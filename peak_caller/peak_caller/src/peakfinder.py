#!/nas3/yeolab/Software/Python_dependencies/bin/python
#We will follow the UCSC genome browser assumption of using a zero based half open cord system
import pysam
import optparse
from optparse import OptionParser, SUPPRESS_HELP
import os
import sys
from subprocess import Popen, PIPE, call
import tempfile
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from deap import dtm
import scipy
from scipy import interpolate, optimize, linalg, misc, stats
import pickle
import math
import time
import pybedtools
from random import sample as rs
from seqTools import *
import gzip
import peaks
import pkg_resources
host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
os.system("echo $PATH")


"""

Checks to make sure a BAM file has an index, if the index does not exist it is created

Usage undefined if file does not exist (check is made earlier in program)
bamfile - a path to a bam file

Returns 1 
TODO make it so a failaure returns 0

"""
def check_for_index(bamfile):
    
    if not os.path.exists(bamfile):
        raise NameError("file %s does not exist" % (bamfile))
    
    if os.path.exists(bamfile + ".bai"):
        return 1
    else:
        print "Index for %s does not exist, indexing bamfile" %(bamfile)
        process = call(["samtools", "index", str(bamfile)])
        
        #error if file is not correct format
        print process
        print type(process)
        if process == -11: 
            raise NameError("file %s not of correct type" % (bamfile))
        
        return 1

def get_FDR_cutoff_mode(readlengths, genelength, iterations=1000, mincut = 2, alpha=.05):
    """
    Find randomized method, as in FOX2ES NSMB paper.
    """
    if readlengths.__len__() < 20: # if you have very few reads on a gene, don't waste time trying to find a cutoff
        return mincut
    cmd = "./peaks"
    bad = 1
    tries=0
    while bad == 1 and tries < 5:
        try:
            process = Popen([cmd, "-f", "stdin", "-L",str(genelength),"-r", str(iterations), "-a", str(alpha)], stdin=PIPE, stdout=PIPE)
            results, err = process.communicate("\n".join(map(str, readlengths)))
            return_val = process.wait()
            bad=0
        except OSError:
            print "Couldn't open a process for thresholding, trying again"
            tries +=1
        
    if bad == 1:
        return "error"
    obs = 0
    cutoff = mincut
    for x in results.split("\n"):
        if x == "":
            continue
        try:
            cut, n_observed = map(int, x.strip().split("\t"))
        except:
            pass
        if n_observed > obs and cut > cutoff:
            obs = n_observed
            cutoff = cut
    if cutoff < mincut:
        cutoff = mincut            
    return int(cutoff)

def get_FDR_cutoff_mean(readlengths, genelength, iterations=1000, mincut = 2, alpha = 0.05):
    """
    Find randomized method, as in FOX2ES NSMB paper.
    MEAN, not MODE
    scatter reads, calcaluate number of reads to pass fdr threshold, takes average observed cutoff
    readlengths -- list of lengths of aligned portions of reads
    genelength -- effective gene length (unalignable regions aren't counted)
    interations -- number of times to repeat FDR thresholding calculation 
    mincut -- min threshold possible to return
    alpha -- FDR alpha 
    
   
    
    Returns an int, the number of reads needed to meet the FDR cutoff
    TODO: Allow the minimum cutoff to be paramaritizied
    TODO: double check math on this
    """
    
    if readlengths.__len__() < 20: # if you have very few reads on a gene, don't waste time trying to find a cutoff
        return mincut
    results = peaks.shuffle(genelength, iterations, 0, .05, readlengths) 
    
    total = 0
    n =0
    
    #parses results from peaks script, calculates mean from peaks results 
    #should document peaks function call return value somewhere around here
    
    for cut, n_observed in enumerate(results):
        total += (cut * n_observed)
        n+=n_observed
        
    #logic for min cutoffs 
    cutoff = total/iterations
    if cutoff < mincut:
        cutoff = mincut
    return int(round(cutoff, 0))

"""

Loads bed file into a dictionary with the key being the name and a string being the value

Input:
BED -- a bed file to load

Return:
A dictionary with the key being the name position of the bed file and the values being the
ordered bed file

TODO: Refactor to used bedtools instead

"""
def build_geneinfo(BED):
    
    #opens bed file, either zipped or unzipped
    try:
        bedfile = gzip.open(BED, "rb")
    except:
        bedfile = open(BED, "r")
        
    GI = dict()
    
    for line in bedfile.readlines():
        chr, start, stop, name, score, signstrand = line.strip().split("\t")
        chr.replace("chr", "")
        GI[name] = "|".join([chr, name, start, stop, str(signstrand)])
    bedfile.close()
    return GI

"""

Builds a dictionary of gene names and lengths of mappable regions in that gene

Input:
A two column file with the first column being the gene name and the second column being the
mappable length of the gene

Return:
A dictionary with the key being the name of the gene and the value being the length

"""
def build_lengths(file):
    FI=open(file,"r")
    LEN = dict()

    for line in FI.readlines():
        name, len = line.strip().split("\t")
        LEN[name]=int(len)

    FI.close()
    return LEN

"""


Finds contiguous (within margin) regions that have reads, the area between covered locations
is defined as regions without any coverage

Input:
data - wiggle track in list form each value is the coverage at that location
margin - distance between sections 

Output:
A list of strings in the form "start_location|stop_location"
ie [50|60, 70|80] 

TODO: Modify to allow for thresholded margins 
"""     
def find_sections(data, margin):

    sections = list()
    section_num =0
    start = 0
    stop = 0
    highlight=False
    gap=0
    margin = int(margin)
    
    #walk along wiggle track until a gap is length of margin, when that happens reset, call that a region
    #and reset
    for i, val in enumerate(data):
        stop =i
        if val > 0:
            gap=0
            if highlight is False:
                start = i - margin
                if start < 0:
                    start =0
            highlight=True
        else:
            gap += 1
        if highlight is True and gap > margin:
            #reset
            highlight=False
            gap=0
            sect = str(start) + "|" + str(stop)
            sections.append(sect)
    
    #catches last section
    #really need to get rid of all this string parsing
    if highlight is True:
        sect = str(start) + "|" + str(stop)
        sections.append(sect)
    return(sections)
        
def plotSpline(spline, data, xvals, threshold=None):
    """
    Plot a smoothing spline and real data
    """
    f = plt.figure()
    ax1 = f.add_subplot(111)
    ax1.plot(xvals, data)
    ax1.plot(xvals, spline.__call__(xvals))
    if threshold is not None:
        ax1.axhline(y=threshold)
    f.show()
    
    pass

"""

Calculate a smoothing spline for ydata in xdata then return complexity-penalized residuals
or return the smoothing spline if resid ==False

Parameters:
x -- Int, range of spline
xdata -- list, wiggle track positions   
ydata -- list, wiggle track coverage 
k -- int, degree of spline

Output: spline object and error if asked for by setting resid, probably want to factor out error calcluation into 
second function

refer to scipy documentation for futher questions.
This functions results are undefined if all requiered argunments are not supplied
"""
def find_univariateSpline(x, xdata, ydata, k, weight=None, resid=True):

    try:
        spline = scipy.interpolate.UnivariateSpline(xdata, ydata, s=x,k=k, w=weight)
        #plotSpline(spline, ydata, xdata, 33)        
        #computationally intensive 
        
        #this section calclauates error of splines based on residuals * number of turns
        #need to explain why this happens Mike
        if resid is True:
            #knots = spline.get_knots()
            func = spline.__call__(xdata)
            turns = sum(abs(diff(sign(diff(func)))))/2 # number of turns in the function
            err = sqrt((spline.get_residual())*(turns**4))
            #may also work.
            ###norm = linalg.norm(func)
            ###err = sqrt((spline.get_residual()) + norm**3) #penalize complexity
            #print str(v) + "\t" + str(turns)
            return(err)
        else:
            return(spline)
    except:
        return(Inf)

def plotSections(wiggle, sections, threshold):
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(wiggle)
    ax.axhline(y=threshold)
    for sect in sections:
        #mark active sections
        positions = list() 
        codes = list()
        start, stop = map(int, sect.split("|"))
        positions.append([start, 0.6])
        codes.append(Path.MOVETO)
        positions.append([stop, 0.6])
        codes.append(Path.LINETO)
        positions.append([stop, 0.05])
        codes.append(Path.LINETO)
        positions.append([start, 0.05])
        codes.append(Path.LINETO)
        positions.append([start, 0.6])
        codes.append(Path.LINETO)
        positions.append([start, 0.6])
        codes.append(Path.CLOSEPOLY)
        path = Path(positions, codes)
        patch = patches.PathPatch(path, lw=1)
        ax.add_patch(patch)
    f.show()

    """
    
    scipy.stats.poisson.cdf
    compute the p-value for a peak of peak_length length with reads_in_peak reads,
    given that the read is gene_length long and has reads_in_gene reads

    If there are fewer than 3 reads expected to fall in the region, assume there's 3 reads
    expected...
    
    Paramaters
    ----------
    reads_in_gene: Integer representing number of reads in gene
    reads_in_peak: Integer reperesnting the number of reads in a specific peak
    gene_length: Integer representing length of gene
    peak_length: Integer representing length of peak
    
    Returns double, the p-value that the peak is significant
    If calcluation fails returns 1
    
    """
def poissonP(reads_in_gene, reads_in_peak, gene_length, peak_length):

    try:
        #lam is estimate of the lambda value
        #poission takes a value and the lambda 
        #this is average number of reads per single site in the gene, but
        #a peak is not a single site, so it the average number gets multipled by the peak 
        #length as an estimator of the mean
        
        #TODO: check with boyko or someone else about this math.  I think its a strong over 
        #estimate of the mean
        lam = (float(reads_in_gene)/(gene_length))*(peak_length)

        if lam < 3:
            lam =3;
        cumP = 1-scipy.stats.poisson.cdf(reads_in_peak, int(lam))
        return cumP
    except:
        return 1

"""

calls peaks for an individual gene 

loc - string of all gene locations
gene_length - effective length of gene
takes bam file or bam file object.  Serial uses object parallel uses location (name)
trim collapses redundant reads (same start and stop position) --might not belong here
margin - space between sections for calling new peaks
FDR_alpha - false discovery rate, p-value bonferoni correct from peaks script (called in setup)
user_threshold - user defined FDR thershold (probably should be factored into FDR_alpha
minreads - min reads in section to try and call peaks
poisson_cutoff - p-value for signifance cut off for number of reads in peak that gets called - might want to use ashifted distribution
plotit - makes figures 
quiet - supresses output
outfile - ???
w_cutoff - width cutoff, peaks narrower than this are discarted 
windowssize - for super local calculation distance left and right to look 
SloP - super local p-value instead of gene-wide p-value
correct_P - boolean bonferoni correction of p-values from poisson

"""
def call_peaks(loc, gene_length, bam_fileobj=None, bam_file=None, trim=False, margin=25, FDR_alpha=0.05,user_threshold=None,
               minreads=20, poisson_cutoff=0.05, plotit=False, quiet=False, outfile=None, w_cutoff=10, windowsize=1000, SloP = False, correct_P = False):
    #setup
    chrom, gene_name, tx_start, tx_end, signstrand = loc.split("|")
    
    #logic reading bam files
    if bam_file is None and bam_fileobj is None:
        #using a file object is faster for serial processing bot doesn't work in parallel
        print "you have to pick either bam file or bam file object, not both"
        exit()
    elif bam_fileobj is None:
        bam_fileobj = pysam.Samfile(bam_file, 'rb')
        
    tx_start, tx_end = map(int, [tx_start, tx_end])
    subset_reads = bam_fileobj.fetch(reference=chrom, start=tx_start,end=tx_end)

    #need to document reads to wiggle
    wiggle, jxns, pos_counts, lengths, allreads =readsToWiggle_pysam(subset_reads,tx_start, tx_end, keepstrand=signstrand, trim=trim)
    
    
    r = peaks_from_info(wiggle,pos_counts, lengths, loc, gene_length, trim, margin, FDR_alpha,user_threshold,minreads, poisson_cutoff, plotit, quiet, outfile, w_cutoff, windowsize, SloP, correct_P)

    return r


"""

same args as before 
wiggle is converted from bam file
pos_counts - one point per read instead of coverage of entire read
lengths - lengths aligned portions of reads 
rest are the same fix later

"""
def peaks_from_info(wiggle, pos_counts, lengths, loc, gene_length, trim=False, margin=25, FDR_alpha=0.05,user_threshold=None,
                                   minreads=20, poisson_cutoff=0.05, plotit=False, quiet=False, outfile=None, w_cutoff=10, windowsize=1000, SloP = False, correct_P = False):


    peakDict = {}
    #these are what is built in this dict, complicated enough that it might be worth turning into an object
    #peakDict['clusters'] = {}
    #peakDict['sections'] = {}
    #peakDict['nreads'] = int()
    #peakDict['threshold'] = int()
    #peakDict['loc'] = loc
    
    #data munging
    chrom, gene_name, tx_start, tx_end, signstrand = loc.split("|")
    tx_start, tx_end = map(int, [tx_start, tx_end])    
    
    #used for poisson calclulation? 
    nreads_in_gene = sum(pos_counts)
    
    #decides FDR calcalation, maybe move getFRDcutoff mean into c code
    if user_threshold is None:
        gene_threshold = get_FDR_cutoff_mean(lengths, gene_length, alpha=FDR_alpha)
    else:
        gene_threshold = user_threshold
    
    if gene_threshold == "error":
        print "I had a hard time with this one: %s.  I think I'll use a threshold of 50" %(loc)
        threshold=50
    peakDict['clusters'] = {}
    peakDict['sections'] = {}
    peakDict['nreads'] = nreads_in_gene
    peakDict['threshold'] = gene_threshold
    peakDict['loc'] = loc
    peakn=1
    tmpsect = {}
    if quiet is not True:
        print "Testing %s" %(loc)
        print "Gene threshold is: %d" %(gene_threshold)
    sections = find_sections(wiggle, margin)


    if plotit is True:
        plotSections(wiggle, sections, gene_threshold)
    bed = list()

    for sect in sections:
        sectstart, sectstop = map(int, sect.split("|"))
        sect_length = sectstop-sectstart+1
        data = wiggle[sectstart:(sectstop+1)]
        cts = pos_counts[sectstart:(sectstop+1)]
        xvals = arange(0, sect_length)
        Nreads = sum(cts)

        #gets random subset of lengths of reads for calculations on a section
        sect_read_lengths = rs(lengths, Nreads) #not exactly the right way to do this but it should be very close.
        peakDict['sections'][sect] ={}
        threshold=int()
        
        #makes sure there are enough reads
        if Nreads < minreads:
            if quiet is not True:
                print "%d is not enough reads, skipping section: %s" %(Nreads, sect)
            continue
        else:
            if quiet is not True:
                print "Analyzing section %s with %d reads" %(sect, Nreads)
        
            
        #sets super-local if requested, might be able to factor this
        if user_threshold is None:
            if SloP is True:
                threshold = get_FDR_cutoff_mean(sect_read_lengths, sect_length, alpha=FDR_alpha)
                if quiet is not True:
                    print "Using super-local threshold %d" %(threshold)
            else:
                threshold= gene_threshold
        else:
            threshold= user_threshold

        #saves threshold for each individual section
        peakDict['sections'][sect]['threshold'] = threshold
        peakDict['sections'][sect]['nreads'] = Nreads

        #if wiggle track never excides threshold
        if max(data) < threshold:
            if quiet is not True:
                print "data not high enough, stopping"
            continue
        
        #fitting splines logic, black magic 
        try:
            degree=3 #cubic spline
            weights = None 
            fo = True #output information about the fitting optimization
            if quiet is True:
                fo=False
            #for very large windows with many reads a large smoothing parameter is required.  test several different options to determine a reasonable inital estimate
            
             #Goal is to find optimnal smooting paramater in multiple steps
           #x1 initial estimate of smoothing paramater 
            #step 1, identify good initial value
            x1=(sectstop-sectstart+1)
            x0=x1
            useme = 1
            
           
            
            #step 2, refine so as not to runinto local minima later, try to come up with a good way of getting optimal paramater
            error  = find_univariateSpline(x1, xvals, data, degree, weights, resid=True)
            for i in range(2, 11):
                x2 = x1*i
                
                #tries find optimal initial smooting paraater in this loop
                err = find_univariateSpline(x2, xvals, data, degree, weights, resid=True)
                if err < error:
                    x0 = x2
                    useme = i

            if quiet is not True:
                print "I'm using (region length) * %d as the initial estimate for the smoothing parameter" %(useme)            
            try:
                #fine optimization of smooting paramater
                cutoff =float(0)
                tries =0
                while cutoff <5:# shouldn't get smoothing coef's this small.. increase the initial estimate and try again. WARNING: BLACK MAGIC
                    tries += 1
                    if tries == 3: # increasing this may improve accuracy, but at the cost of running time.
                        break
                    sp = scipy.optimize.minimize(find_univariateSpline, x0, args=(xvals, data, degree, weights),
                                                 options={'disp':fo}, method="Powell")
                    #fit a smoothing spline using an optimal parameter for smoothing and with weights proportional to the number of reads aligned at each position if weights is set
                    if sp.success is True:
                        cutoff = sp.x
                    else:
                        pass
                    x0 += sect_length
            except:
                print "%s failed spline fitting at section %s with sp: %s" %(loc, sect, str(sp))
                continue

            if quiet is not True:
                print "optimized smoothing parameter"
        #if we are going to save and output as a pickle fi is %s" %(str(cutoff))
        #final fit spline
            spline = find_univariateSpline(cutoff, xvals, data, degree, weights, resid=False)
            if plotit is True:
                plotSpline(spline, data, xvals, threshold)
            
            #function inside function, remove later
            def thresh(threhold):
                return threshold
            #starts = xvals[diff(sign(spline(xvals) - spline(xvals+1))) < 0]
            
            #finds all turns 
            starts = xvals[diff(sign(spline(xvals) - thresh(xvals))) > 0]
            stops = xvals[diff(sign(spline(xvals) - thresh(xvals))) < 0]
            ### important note: for getting values x->y [inclusive] you must index an array as ar[x:(y+1)]|                     or else you end up with one-too-few values, the second index is non-inclusive
            #append local minima:
            local_minima = xvals[diff(sign(spline(xvals) - spline(xvals+1))) < 0]
            
            #append to list any local minima above threshold
            if any(local_minima >= threshold):
                startlist = list(starts)
                stoplist = list(stops)                
                for val in local_minima:
                    if spline(val) >= threshold:
                        startlist.append(val)
                        stoplist.append(val)
                starts = array(sorted(startlist))
                stops = array(sorted(stoplist))

            #make sure that the start is not a minima 
            if spline(xvals)[0] > threshold:# add a 0 to the beginning of "starts"
                l = list(starts) ## this is HACKED... i couldn't figure out a clean way to do it.
                l.append(0)
                for i in starts:
                    l.append(i)
                starts = array(l)
            
            #removes duplicates 
            starts = array(sorted(set(starts)))
            stops = array(sorted(set(stops)))

            #plt.plot(starts)
            #plt.plot(stops)
            #plt.draw()
            
            #walks along spline, and calls peaks along spline
            #for each start, take the next stop and find the peak between the start and the stop
            for p_start in starts: #subsections that are above threshold
                
                try:
                    p_stop = stops[stops > p_start][0]
                except:
                    p_stop = sect_length -1
                try:
                    peaks = map(lambda x: x+p_start, xvals[diff(sign(diff(spline(xvals[p_start:(p_stop+1)]))))<0])  #peaks with-in this subsection, indexed from section (not subsection) start
                    if quiet is not True:
                        print "I found %d peaks" %(len(peaks))
                except:
                    continue
                
                #handles logic if there are multiple peaks between start and stop
                if peaks.__len__() <=0:
                    continue
                if peaks.__len__() is 1:
                    #gets reads in peak
                    Nreads_in_peak = sum(cts[p_start:(p_stop+1)])
                    if quiet is not True:
                        print "Peak %d - %d has %d reads" %(p_start, (p_stop+1), Nreads_in_peak)
                    
                    #makes sure there enough reads
                    if Nreads_in_peak < minreads or max(data[p_start:(p_stop+1)]) < threshold:
                        if quiet is not True:
                            print "skipping peak, %d is not enough reads" %(Nreads_in_peak)
                        continue

                    #formatting of bed track
                    #start and stop for bed track to be created
                    g_start = tx_start + sectstart + p_start
                    g_stop = tx_start + sectstart + p_stop
                    
                    #highest point in start stop
                    peak = tx_start + sectstart + peaks[0]
                    
                    #makes it thicker so we can see on the browser 
                    thick_start = peak-2
                    thick_stop = peak+2
                    
                    #error checking logic to keep bed files from breaking
                    if thick_start < g_start:
                        thick_start = g_start
                    if thick_stop > g_stop:
                        thick_stop = g_stop
                    peak_length = g_stop-g_start+1

                    #
                    if peak_length < w_cutoff:#skip really small peaks
                        continue
                    peak_name = gene_name + "_" + str(peakn) + "_" + str(int(Nreads_in_peak))
                    
                    #super local logic 
                    #error check to make sure area is in area of gene
                    if peak - tx_start - windowsize < 0: #distance from gene start
                        area_start = 0
                    else:  #for super local gets area around peak for calculation
                        area_start = peak - tx_start - windowsize
                        #area_start = sectstart
                        
                    #same thing except for end of gene instead of start
                    if peak + windowsize > tx_end: #distance to gene stop
                        area_stop = tx_start-tx_end + 1
                    else:
                        area_stop = peak - tx_start + windowsize
                        #area_stop = sectstop

                    #use area reads + 1/2 all other reads in gene: area_reads = sum(pos_counts[area_start:area_stop]) + 0.5*(sum(pos_counts) - sum(pos_counts[area_start:area_stop]))
                    #use area reads:
                    area_reads = sum(pos_counts[area_start:area_stop])
                    area_size = area_stop - area_start + 1

                    #area_reads = sum(pos_counts[sectstart:sectstop])
                    #area_size = sect_length

                    #calcluates poisson based of whole gene vs peak
                    gene_poisP = poissonP(nreads_in_gene, Nreads_in_peak, gene_length, peak_length)
                    if SloP is True:
                        #same thing except for based on super local p-value
                        slop_poisP = poissonP(area_reads, Nreads_in_peak, area_size, peak_length)
                        
                    #makes sure spop_poisP is defined, even if its just normal, something to be removed later,
                    #slop should only be used when defined as true
                    else:
                        slop_poisP=gene_poisP
                    
                    
                    if math.isnan(slop_poisP):
                        slop_poisP = 1
                                            
                    #remove later    
                    if slop_poisP > poisson_cutoff:
                        #continue
                        pass
                    
                    #poisP = 1
                    
                    #defines the bedline of a peak for returning                    
                    bedline = "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d" %(chrom, g_start, g_stop, peak_name, slop_poisP, signstrand, thick_start, thick_stop)

                    #metadata for the specific bedline
                    peakDict['clusters'][bedline] = {}
                    peakDict['clusters'][bedline]['GeneP'] = gene_poisP
                    peakDict['clusters'][bedline]['SloP'] = slop_poisP                    
                    peakDict['clusters'][bedline]['Nreads'] = Nreads_in_peak
                    peakDict['clusters'][bedline]['size'] = peak_length
                    

                    peakn += 1
                else:  #there are more than one peaks in this window
                    #this handles peaks within peaks logic
                    valleys = array(map(lambda x:x+p_start, xvals[diff(sign(diff(spline(xvals[p_start:p_stop+1]))))>0]))#local minima in subsection, relative to section start

                    for subpeak in peaks:
                        subpeak_start = int()
                        subpeak_stop = int()
                        if any(valleys < subpeak):
                            subpeak_start = valleys[valleys < subpeak][-1]
                        else:
                            subpeak_start = starts[starts < subpeak][-1]
                        if any(valleys > subpeak):
                            subpeak_stop = valleys[valleys > subpeak][0]
                        else:
                            subpeak_stop = stops[stops > subpeak][0]
                        peak_length = subpeak_stop - subpeak_start + 1
                        if peak_length < w_cutoff:#skip really small peaks
                            continue
                        Nreads_in_peak = sum(cts[subpeak_start:(subpeak_stop+1)])
                        if Nreads_in_peak < minreads or max(data[subpeak_start:(subpeak_stop+1)]) < threshold:
                            continue
                        g_start = tx_start + subpeak_start + sectstart
                        g_stop = tx_start + subpeak_stop + sectstart
                        peak = tx_start + subpeak + sectstart
                        thick_start = peak-2
                        if thick_start < g_start:
                            thick_start = g_start                        
                        thick_stop = peak+2
                        if thick_stop > g_stop:
                            thick_stop = g_stop                        
                        peak_name = gene_name + "_" + str(peakn) + "_" + str(int(Nreads_in_peak))
                        if peak - tx_start - windowsize < 0: #distance from gene start
                            area_start = 0 
                        else:
                            area_start = peak - tx_start - windowsize
                            #area_start = sectstart
                        if peak + windowsize > tx_end: #distance to gene stop
                            area_stop = tx_start-tx_end + 1
                        else:
                            #area_stop = sectstop
                            area_stop = peak - tx_start + windowsize
                        
                        #area_reads = sum(pos_counts[area_start:area_stop]) + 0.5*(sum(pos_counts) - sum(pos_counts[area_start:area_stop])) #all the reads in the area + half the reads in the rest of the gene.
                        area_reads = sum(pos_counts[area_start:area_stop])                            
                        area_size = area_stop - area_start + 1
                        #area_reads = sum(pos_counts[sectstart:sectstop])
                        #area_size = sect_length
                        
                        gene_poisP = poissonP(nreads_in_gene, Nreads_in_peak, gene_length, peak_length)
                        if SloP is True:
                            slop_poisP = poissonP(area_reads, Nreads_in_peak, area_size, peak_length)
                        else:
                            slop_poisP = gene_poisP
                        
                        if math.isnan(slop_poisP):
                            slop_poisP = 1
                        if slop_poisP > poisson_cutoff: #we'll leave these in to allow for BH p-value correction
                            pass
                        
                        #output results again
                        bedline = "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d" %(chrom, g_start, g_stop, peak_name,
                                                                     slop_poisP, signstrand, thick_start, thick_stop)
                        peakDict['clusters'][bedline] = {}                        
                        peakDict['clusters'][bedline]['SloP'] = slop_poisP
                        peakDict['clusters'][bedline]['GeneP'] = gene_poisP
                        peakDict['clusters'][bedline]['Nreads'] = Nreads_in_peak
                        peakDict['clusters'][bedline]['size'] = peak_length
                        peakn+=1
        except:
            print "spline fitting failed for %s" %(loc)
            
    #inflate p-values based on # of comparisons #bonferroni corrected
    if correct_P is True:            
        for peak in peakDict['clusters']:
            peakDict['clusters'][peak]['p'] = peakDict['clusters'][peak]['p'] * peakn  #bonferroni correct p-value for MHT
        
        

    peakDict['Nclusters'] = peakn
    
    return peakDict



"""

Creates a dictionary containing all information needed to perform peak calling calcluations 
for a single species

Paramaters
-----------
species: string currently not used
chrs: list specifying all the chromosomes in a given species
bed: path to a bed file that contains information on genes (custom file *STRUCTURE_genes.BED.gz)
mrna: path to a file that contains mRNA lengths (custom CSV file contains gene names follwed by gene lengths)
premrna: path to a file that contains pre-mRNA lengths (custom CSV file contains gene names follwed by gene lengths_

Returns dict of all items passed to it

TODO:  Add checking to verify that file are actually passed
"""
def add_species(species, chrs, bed, mrna, premrna):
        par = dict()
        
        #this is non-pythonic, should just combine all lists
        par["chrs"] = [item for sublist in chrs for item in sublist] #expand sublists
        par["gene_bed"] = bed
        par["mRNA"] = mrna
        par["premRNA"] = premrna
        return par
        
def main(options):
    bamfile = options.bam
    quiet = options.quiet
    start = True
    if os.path.exists(bamfile):
        bamfile = os.path.abspath(bamfile) #re-set to include the full path to bamfile
        if quiet is not True:
            #sys.stderr.write("bam file is set to %s\n" %bamfile)
            print "bam file is set to %s\n" %(bamfile)
    else:
        print "Bam file not defined"

    species = options.species
    geneBed = options.geneBEDfile
    genemRNA = options.geneMRNAfile
    genePREmRNA = options.genePREMRNAfile

    species_parameters = dict()

    if species is None and geneBed is None:
        print "You must set either \"species\" or \"geneBed\"+\"geneMRNA\"+\"genePREMRNA\""
        exit()

    lenfile = ""
    
    species_parameters["hg19"] = add_species("hg19", [range(1,22), "X", "Y"],
                                             pkg_resources.resource_filename(__name__, "../data/hg19.AS.STRUCTURE_genes.BED.gz"),
                                             pkg_resources.resource_filename(__name__, "../data/hg19.AS.STRUCTURE_mRNA.lengths"),
                                             pkg_resources.resource_filename(__name__, "../data/hg19.AS.STRUCTURE_premRNA.lengths"))
    species_parameters["hg18"] = add_species("hg18", [range(1,22), "X", "Y"],
                                             pkg_resources.resource_filename(__name__, "../data/hg18.AS.STRUCTURE_genes.BED.gz"),
                                             pkg_resources.resource_filename(__name__, "../data/hg18.AS.STRUCTURE_mRNA.lengths"),
                                             pkg_resources.resource_filename(__name__, "../data/hg18.AS.STRUCTURE_premRNA.lengths"))
    species_parameters["mm9"] = add_species("mm9", [range(1,19), "X", "Y"],
                                            pkg_resources.resource_filename(__name__,"../data/mm9.AS.STRUCTURE_genes.BED.gz"),
                                            pkg_resources.resource_filename(__name__,"../data/mm9.AS.STRUCTURE_mRNA.lengths"),
                                            pkg_resources.resource_filename(__name__,"../data/mm9.AS.STRUCTURE_premRNA.lengths"))
    acceptable_species = ",".join(species_parameters.keys())
    
    #error checking
    if species is not None and geneBed is not None:
        print "You shouldn't set both geneBed and species, defaults exist for %s" %(acceptable_species)
        exit()
    if species is not None and species not in species_parameters:
        print "Defaults don't exist for your species: %s. Please choose from: %s or supply \"geneBed\"+\"geneMRNA\"+\"genePREMRNA\"" %(species, acceptable_species)
        exit()
    if species is None:
        species = "custom"
        species_parameters["custom"] = add_species("custom", [range(1,22), "X", "Y"], geneBed, genemRNA, genePREmRNA)

    #error checking done, this does... something.  This is more setup phase  Uses pre-mrnas?
    if options.premRNA is True:
        lenfile=species_parameters[species]["premRNA"]
    else:
        lenfile=species_parameters[species]["mRNA"]
    
    #builds dict to do processing on, 
    lengths = build_lengths(lenfile)
    genes = build_geneinfo(species_parameters[species]["gene_bed"])
    margin = int(options.margin)
    
    #this should be fixed, args should initally be ints if passed
    try:
        
        maxgenes = int(options.maxgenes)
    except:
        pass
    
    #unwrapping the options is really unnessessary, this could be refactored to remove the unwrapping
    job_name = options.job_name
    prefix = options.prefix
    minreads = int(options.minreads)
    plotit = options.plotit
    poisson_cutoff = options.poisson_cutoff

    g = options.gene ##
    
    #Chooses between parallel and serial runs.  Potentally good to factor these two things out into different function calls
    if options.serial is True:
        #this warning looks incorrect, serial is for processing, not lists of genes
        print "WARNING: not using transcriptome-based cutoff because you decided to supply a list of genes"
        bamfileobj = pysam.Samfile(bamfile, 'rb')
        n = 0 #number of genes analyzed
        genelist = list()
        #allpeaks=pybedtools.BedTool(open(options.outfile, 'w'))
        allpeaks = dict()
        
        #this is all super hacky and needs to be cleaned up, this tests for genes or keys, needs to be more explicit
        try:
            if g.__len__() > 0:
                genelist = g
                print "Trying these genes:"
                for gene in g:
                    print gene

        except:
            genelist = genes.keys()
            
        #unused?
        alleaks = dict()

        n=0
        
        #unused?
        n_clusters= 0
        
        #I think this calls peaks for each gene in the gene list, which could be every gene in the genome
        for gene in genelist:
            n= n+1
            
            #again, hacky should be factored to a single if statement, need to be more explicit about code paths
            if options.maxgenes == None:
                pass
            else:
                if n >= maxgenes:
                    break
            geneinfo = genes[gene]
            genelen = lengths[gene]
            
            #There is a better way of doing timing.  
            t=time.strftime('%X %x %Z')
            print geneinfo+ " started:"+str(t)
            trim=options.trim
            
            #maybe just pass options object?  Should this be a new object in itself?
            gene_peaks = call_peaks( geneinfo, genelen, bam_fileobj=bamfileobj, bam_file=None, trim=trim, user_threshold = options.threshold,
                                     margin=margin, minreads=minreads, poisson_cutoff=poisson_cutoff,plotit=plotit, quiet=quiet,outfile=None, SloP=options.SloP, FDR_alpha = options.FDR_alpha)
            
            #unused?
            n_clusters = gene_peaks['Nclusters']
            
            #prints out all peaks for an individual gene, this should probably be its own function
            for i in gene_peaks['clusters'].keys():

                #note: you must set "correctP" in the call_peaks function to "True" if you'd like the p-values to actually be corrected
                corrected_SloP_pval = gene_peaks['clusters'][i]['SloP']
                corrected_Gene_pval = gene_peaks['clusters'][i]['GeneP']                
                #this i should be an object if you are processing it like this
                (chrom, g_start, g_stop, peak_name, geneP, signstrand, thick_start, thick_stop) = i.split("\t")                            
                
                #only prints peak if it passes quality cutoff requiernments.  
                if corrected_SloP_pval < poisson_cutoff or corrected_Gene_pval < poisson_cutoff:
                    min_pval = min([corrected_SloP_pval, corrected_Gene_pval])

                    bedline = "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d" %(chrom, int(g_start), int(g_stop), peak_name, min_pval, signstrand, int(thick_start), int(thick_stop))
                    #this is the biggest hack I've ever seen.  Convert it out then convert it out of string then convert it back for hashing?  This needs to be a hashable object
                    allpeaks[bedline] = 1
                else:
                    if quiet is not True:
                        print "Failed Gene Pvalue: %s Failed SloP Pvalue: %s for cluster %s" %(corrected_Gene_pval, corrected_SloP_pval, i)
                    pass
            t=time.strftime('%X %x %Z')
            print geneinfo+ " finished:"+str(t)
        outbed = options.outfile + ".BED"
        color=options.color

        
        #import code
        #code.interact(local=locals())

        
        #creates a bedtool out of keys, this should really be whats getting appended to in the inital construction step
        tool = pybedtools.BedTool("\n".join(allpeaks.keys()), from_string=True)
        
        #gracefully fails on zero case 
        if len(allpeaks.keys()) > 0:
            tool = tool.sort()
        
        
        tool.saveas(outbed, trackline="track name=\"%s\" visibility=2 colorByStrand=\"%s %s\"" %(outbed, color, color))   
        
        #outputs results as well as printing them, bad form.  We should probably just print to stdout and let the user decide how to save the results.  
        print tool
    else:
        
        #factor out into error checking 
        if quiet is False:
            print "quiet required for parallel runs"
        if plotit is True:
            print "plots not available for parallel runs"
        
        #badly named variables 
        gl = list()
        ll = list()
        g = options.gene
        genelist = list() # list of genes to call peaks on

        #this is redundant code should be factored, also len() __len__ is a protected function!
        try:
            if g.__len__() > 0:
                genelist = g

        except:
            genelist = genes.keys()
        n=0
        transcriptome_size = 0
        
        #redundant code, need to factor
        for i in genelist: #make sure that the lists are in the right order, and if there's maxgenes set then limit the number of genes to analyze
            n+=1
            if options.maxgenes is None:
                pass
            else:
                if n > maxgenes:
                    break
                
            #gene lengths and length length lengths 
            gl.append(genes[i])
            ll.append(lengths[i])
            transcriptome_size += lengths[i]
        print "trying %d genes" %(gl.__len__())
        trim=options.trim
        
        results = dtm.map(call_peaks, gl, ll, bam_file=bamfile, bam_fileobj=None, trim=trim, margin=margin, user_threshold = options.threshold,
                          minreads=minreads, poisson_cutoff=poisson_cutoff, plotit=False, quiet=True, outfile=None, SloP=options.SloP, FDR_alpha=options.FDR_alpha)
        print "finished with calling peaks"
        
        #if we are going to save and output as a pickle file we should output as a pickle file
        #we should factor instead create a method or object to handle all file output
        if options.save_pickle is True:
            pickle_file = open(options.outfile + ".pickle", 'w')
            pickle.dump(results, file = pickle_file)                
        
        #combine results
        allpeaks = {}

        #count total number of reads in transcriptiome
        transcriptome_reads = 0
        
        for gener in results:
            transcriptome_reads += gener['nreads']
        print "Transcriptome size is %d, transcriptome reads are %d" %(transcriptome_size, transcriptome_reads)
        
        #is this a missed indent?
        for gener in results:
            try:
                #how come this logic for printing clusters is different from the serial logic?
                for cluster in gener['clusters'].keys():
                    try:
                        transcriptomeP = poissonP(transcriptome_reads, gener['clusters'][cluster]['Nreads'], transcriptome_size, gener['clusters'][cluster]['size'])
                        if math.isnan(transcriptomeP):
                            print "Transcriptome P is NaN, transcriptome_reads = %d, cluster reads = %d, transcriptome_size = %d, cluster_size = %d" %(transcriptome_reads, gener['clusters'][cluster]['Nreads'], transcriptome_size, gener['clusters'][cluster]['size'])

                        if transcriptomeP > poisson_cutoff:
                            #print "%s\n Failed Transcriptome cutoff with %s reads, pval: %s" %(cluster, gener['clusters'][cluster]['Nreads'], transcriptomeP)
                            
                            continue
                        min_pval = 1

                        corrected_SloP_pval = gener['clusters'][cluster]['SloP']
                        corrected_Gene_pval = gener['clusters'][cluster]['GeneP']
                        

                        if corrected_SloP_pval < poisson_cutoff or corrected_Gene_pval < poisson_cutoff:
                            min_pval = min([corrected_SloP_pval, corrected_Gene_pval])
                        else:
                            #print "Failed Gene Pvalue: %s and failed SloP Pvalue: %s for cluster %s" %(corrected_Gene_pval, corrected_SloP_pval, i)
                            pass
                            continue


                        (chrom, g_start, g_stop, peak_name, geneP, signstrand, thick_start, thick_stop) = cluster.split("\t")                            
                        bedline = "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d" %(chrom, int(g_start), int(g_stop), peak_name, min_pval, signstrand, int(thick_start), int(thick_stop))
                        allpeaks[bedline] = 1

                    except:
                        print "parsing failed"
                        pass
            except:
                print "error handling genes"
                pass
            
        #again redundant code 
        outbed = options.outfile + ".BED"
        color=options.color
        tool = pybedtools.BedTool("\n".join(allpeaks.keys()), from_string=True).sort().saveas(outbed, trackline="track name=\"%s\" visibility=2 colorByStrand=\"%s %s\"" %(outbed, color, color))

        print "wrote peaks to %s" %(options.outfile)

        return 1

if __name__ == "__main__":
    
    usage="\npython peakfinder.py -b <bamfile> -s <hg18/hg19/mm9>\n OR \npython peakfinder.py -b <bamfile> --customBED <BEDfile> --customMRNA <mRNA lengths> --customPREMRNA <premRNA lengths>"
    description="Fitted Accurate Peaks. Michael Lovci 2012. CLIP peakfinder that uses fitted smoothing splines to define clusters of binding.  Computation is performed in parallel using MPI.  You may decide to use the pre-built gene lists by setting the --species parameter or alternatively you can define your own list of genes to test in BED6/12 format and provide files containing the length of each gene in PREMRNA form and MRNA form (both are required). Questions should be directed to michaeltlovci@gmail.com."
    parser = OptionParser(usage=usage, description=description)

    parser.add_option("--bam", "-b", dest="bam", help="A bam file to call peaks on", type="string", metavar="FILE.bam")

    parser.add_option("--species", "-s", dest="species", help="A species for your peak-finding")

    parser.add_option("--customBED", dest="geneBEDfile", help="bed file to call peaks on, must come withOUT species and with customMRNA and customPREMRNA", metavar="BEDFILE")
    parser.add_option("--customMRNA", dest="geneMRNAfile", help="file with mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")
    parser.add_option("--customPREMRNA", dest="genePREMRNAfile", help="file with pre-mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")

    parser.add_option("--outdir", dest="prefix", default=os.getcwd(), help="output directory, default=cwd")    
    parser.add_option("--outfile", dest="outfile", default="fitted_clusters", help="a bed file output, default:%default")

    parser.add_option("--gene", "-g", dest="gene", action="append", help="A specific gene you'd like try", metavar="GENENAME")
    parser.add_option("--plot", "-p", dest="plotit", action="store_true", help="make figures of the fits", default=False)
    parser.add_option("--quiet", "-q", dest="quiet", action="store_true", help="suppress notifications")

    parser.add_option("--minreads", dest="minreads", help="minimum reads required for a section to start the fitting process.  Default:%default", default=3, type="int", metavar="NREADS")
    parser.add_option("--margin", dest="margin", type="int", help="find sections of genes within M bases that have genes and perform fitting. Default:%default", default=15, metavar="NBASES")
    parser.add_option("--trim", dest="trim", action="store_true", default=False, help="Trim reads with the same start/stop to count as 1")
    parser.add_option("--premRNA", dest="premRNA", action="store_true", help="use premRNA length cutoff, default:%default", default=False)
    parser.add_option("--poisson-cutoff", dest="poisson_cutoff", type="float", help="p-value cutoff for poisson test, Default:%default", default=0.05, metavar="P")
    parser.add_option("--FDR", dest="FDR_alpha", type="float", default=0.05, help="FDR cutoff for significant height estimation, default=%default")
    parser.add_option("--threshold", dest="threshold", type="int", default=None, help="Skip FDR calculation and set a threshold yourself")

    parser.add_option("--serial", dest="serial", action="store_true", help="run genes in sequence (not parallel)")
    parser.add_option("--maxgenes", dest="maxgenes", default=None, help="stop computation after this many genes, for testing", metavar="NGENES")
    parser.add_option("--job_name", dest="job_name", default="FAP", help="name for submitted job. Not used with --serial.  default:%default", metavar="NAME")
    parser.add_option("--processors", dest="np", default=32, help="number of processors to use. Not used with --serial.  default:%default", type="int", metavar="NP")
    parser.add_option("--notify", dest="notify", default=None, help="email address to notify of start, errors and completion", metavar="EMAIL")
    parser.add_option("--superlocal", action = "store_true", dest="SloP", default=False, help="Use super-local p-values, counting reads in a 1KB window around peaks")
    parser.add_option("--color", dest="color", default="0,0,0", help="R,G,B Color for BED track output, default:black (0,0,0)")
    parser.add_option("--start", dest="start", default=False, action="store_true", help=SUPPRESS_HELP) #private, don't use
    parser.add_option("--save-pickle", dest="save_pickle", default=False, action = "store_true", help="Save a pickle file containing the analysis")
    (options,args) = parser.parse_args()
    
    #enforces required usage
    if not (options.bam and ((options.species) or (options.geneBEDfile and options.geneMRNAfile and options.genePREMRNAfile))):
        print "to helpful"
        parser.print_help()
        exit()
        
    check_for_index(options.bam)
    
    if options.serial is True:
        print "Starting serial computation"        
        main(options)
    else:

        if options.start is True:
            dtm.start(main, options)
        else:
            scriptName = os.path.join(options.prefix, options.job_name+".runme.sh")
            runerr = os.path.join(options.prefix, options.job_name+ ".err")
            runout = os.path.join(options.prefix, options.job_name+ ".out")
            shScript = open(scriptName, 'w')
            host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()

            if "optiputer" in host or "compute" in host:
                shScript.write("#!/bin/bash\n#$ -N %s\n#$ -S /bin/bash\n#$ -V\n#$ -pe mpi %d\n#$ -cwd\n#$ -o %s\n#$ -e %s\n#$ -l bigmem\n" %(options.job_name, options.np, runout, runerr))
                if options.notify is not None:
                    shScript.write("#$ -notify\n#$ -m abe\n#$ -M %s\n" %(options.notify))
            
                shScript.write("/opt/openmpi/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines python %s --start\n" %(" ".join(sys.argv)))

            elif "tcc" in host or "triton" in host:
                nnodes = 1
                if int(options.np)%8 == 0:
                    nnodes = int(options.np)/8
                else:
                    nnodes = (int(options.np)/8)+1
                    np = nnodes*8
                    print "You should have used a number of processors that is divisible by 8.  You tried %d and I'll actually use %d." %(options.np, np)
                
                shScript.write("#!/bin/sh\n#PBS -N %s\n#PBS -o %s\n#PBS -e %s\n#PBS -V\n#PBS -S /bin/sh\n#PBS -l nodes=%d:ppn=8\n#PBS -q batch\n#PBS -l walltime=00:50:00\n" %(options.job_name, runout, runerr, nnodes))
                
                if options.notify is not None:
                    shScript.write("#PBS -m abe\n#PBS -M %s\n" %(options.notify))


                shScript.write("let np=$PBS_NUM_NODES*$PBS_NUM_PPN\ncd $PBS_O_WORKDIR\n"); 
                shScript.write("/opt/openmpi/bin/mpirun --mca btl_tcp_if_include myri0 -v -machinefile $PBS_NODEFILE -np $np python %s --start\n" %(" ".join(sys.argv)))
            shScript.close()
            call(["qsub", scriptName]) #subprocess.call
            
                

