#!/nas3/yeolab/Software/Python_dependencies/bin/python
#We will follow the UCSC genome browser assumption of using a zero based half open cord system
import pysam
import optparse
from optparse import OptionParser, SUPPRESS_HELP
import os
import sys
from subprocess import Popen, PIPE, call
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
import gzip
import peaks
import pkg_resources
import pp
host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
job_server = pp.Server()
import logging
logging.disable(logging.INFO)

def verboseprint(*args):
        # Print each argument separately so caller doesn't need to
        # stuff everything to be printed into a single string
            for arg in args:
                print arg,
            print
#verboseprint = lambda *a: None
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
        verboseprint("Index for %s does not exist, indexing bamfile" %(bamfile))
        process = call(["samtools", "index", str(bamfile)])
        
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
            print >>sys.stderr,  "Couldn't open a process for thresholding, trying again"
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

def get_FDR_cutoff_mean(readlengths, genelength, iterations=100, mincut = 2, alpha = 0.05):
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

    if len(readlengths) < 20: # if you have very few reads on a gene, don't waste time trying to find a cutoff        
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
def build_lengths(f):
    FI=open(f,"r")
    gene_lengths = {}

    for line in FI.readlines():
        name, gene_length = line.strip().split("\t")
        gene_lengths[name]=int(gene_length)

    FI.close()
    return gene_lengths

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

def find_sections(data, margin):

    sections = list()
    start = 0
    stop = 0
    highlight = False
    gap = 0
    margin = int(margin)
    
    #walk along wiggle track until a gap is length of margin, when that happens reset, call that a region
    #and reset
    for i, val in enumerate(data):
        stop = i
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
    """
def plotSpline(spline, data, xvals, threshold=None):
    """
    Plot a smoothing spline and real data
    """
    
    verboseprint("plotting")
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
            return(err)
        else:
            return(spline)
    except:
        return(Inf)


"""

Plots each section individually, I think
Wiggle is a list representing a wiggle track
sections is a list of strings of format "start|stop" where start and stop are both integers
threshold is an integer 

"""
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
    plt.show()

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
outfile - ???
w_cutoff - width cutoff, peaks narrower than this are discarted 
windowssize - for super local calculation distance left and right to look 
SloP - super local p-value instead of gene-wide p-value
correct_P - boolean bonferoni correction of p-values from poisson

"""
def call_peaks(loc, gene_length, bam_fileobj=None, bam_file=None, margin=25, FDR_alpha=0.05,user_threshold=None,
               minreads=20, poisson_cutoff=0.05, plotit=False,  outfile=None, w_cutoff=10, windowsize=1000, SloP = False, correct_P = False):
    #setup
    chrom, gene_name, tx_start, tx_end, signstrand = loc.split("|")
    print >> sys.stderr, loc
    #logic reading bam files
    if bam_file is None and bam_fileobj is None:
        #using a file opbject is faster for serial processing bot doesn't work in parallel
        verboseprint("you have to pick either bam file or bam file object, not both")
        exit()
    elif bam_fileobj is None:
        bam_fileobj = pysam.Samfile(bam_file, 'rb')
        
    tx_start, tx_end = map(int, [tx_start, tx_end])
    subset_reads = bam_fileobj.fetch(reference=chrom, start=tx_start,end=tx_end)

    #need to document reads to wiggle
    #wiggle, pos_counts, lengths = readsToWiggle_pysam_foo(subset_reads, tx_start, tx_end, signstrand, "center")
    wiggle, pos_counts, lengths = peaks.readsToWiggle_pysam_foo(subset_reads, tx_start, tx_end, signstrand, "center")
    r = peaks_from_info(list(wiggle), pos_counts, lengths, loc, gene_length, margin, FDR_alpha,user_threshold,minreads, poisson_cutoff, plotit, outfile, w_cutoff, windowsize, SloP, correct_P)

    return r


"""

same args as before 
wiggle is converted from bam file
pos_counts - one point per read instead of coverage of entire read
lengths - lengths aligned portions of reads 
rest are the same fix later

"""
def peaks_from_info(wiggle, pos_counts, lengths, loc, gene_length, margin=25, FDR_alpha=0.05,user_threshold=None,
                                   minreads=20, poisson_cutoff=0.05, plotit=False, outfile=None, w_cutoff=10, windowsize=1000, SloP = False, correct_P = False):

    #screwing around with paraplization 
    import peaks
    from numpy import arange, diff, sign, array
    from random import sample as rs
    import math
    from peakfinder import find_univariateSpline, poissonP
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
        verboseprint("I had a hard time with this one: %s.  I think I'll use a threshold of 50" %(loc))
        threshold=50
    peakDict['clusters'] = {}
    peakDict['sections'] = {}
    peakDict['nreads'] = int(nreads_in_gene)
    peakDict['threshold'] = gene_threshold
    peakDict['loc'] = loc
    peakn=1
 
    verboseprint("Testing %s" %(loc))
    verboseprint("Gene threshold is: %d" %(gene_threshold))
    
    sections = peaks.find_sections(wiggle, margin)


    if plotit is True:      
        plotSections(wiggle, sections, gene_threshold)

    for sect in sections:
        sectstart, sectstop = sect
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
            verboseprint("%d is not enough reads, skipping section: %s" %(Nreads, sect))
            continue
        else:
            verboseprint("Analyzing section %s with %d reads" %(sect, Nreads))
        
            
        #sets super-local if requested, might be able to factor this
        if user_threshold is None:
            if SloP is True:
                threshold = get_FDR_cutoff_mean(sect_read_lengths, sect_length, alpha=FDR_alpha)
                verboseprint("Using super-local threshold %d" %(threshold))
                
            else:
                threshold= gene_threshold
        else:
            threshold= user_threshold

        #saves threshold for each individual section
        peakDict['sections'][sect]['threshold'] = threshold
        peakDict['sections'][sect]['nreads'] = int(Nreads)

        #if wiggle track never excides threshold
        if max(data) < threshold:
            verboseprint("data not high enough, stopping")
            continue
        
        #fitting splines logic, black magic 
        try:
            degree=3 #cubic spline
            weights = None 
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

         
            verboseprint("I'm using (region length) * %d as the initial estimate for the smoothing parameter" %(useme))          
            import scipy
            
            try:
                #fine optimization of smooting paramater
                cutoff =float(0)
                tries =0
                while cutoff <5:# shouldn't get smoothing coef's this small.. increase the initial estimate and try again. WARNING: BLACK MAGIC
                    tries += 1
                    if tries == 3: # increasing this may improve accuracy, but at the cost of running time.
                        break
                    
                    #TODO make verbose a global variable so I can display stuff again
                    sp = scipy.optimize.minimize(find_univariateSpline, x0, args=(xvals, data, degree, weights),
                                                 options={'disp':False}, method="Powell")
                    #fit a smoothing spline using an optimal parameter for smoothing and with weights proportional to the number of reads aligned at each position if weights is set
                    if sp.success is True:
                        cutoff = sp.x
                    else:
                        pass
                    x0 += sect_length
            except e:
                print >>sys.stderr, e
                print >>sys.stderr,  "%s failed spline fitting at section %s with sp: %s" %(loc, sect, str(sp))
                continue

       
            verboseprint ("optimized smoothing parameter")
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
                    verboseprint( "I found %d peaks" %(len(peaks)))
                except:
                    continue
                
                #handles logic if there are multiple peaks between start and stop
                if peaks.__len__() <=0:
                    continue
                if peaks.__len__() is 1:
                    #gets reads in peak
                    Nreads_in_peak = sum(cts[p_start:(p_stop+1)])
                    verboseprint("Peak %d - %d has %d reads" %(p_start, (p_stop+1), Nreads_in_peak))
                    
                    #makes sure there enough reads
                    if Nreads_in_peak < minreads or max(data[p_start:(p_stop+1)]) < threshold:
                        verboseprint("skipping peak, %d is not enough reads" %(Nreads_in_peak))
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
        except :
            print >>sys.stderr,  "spline fitting failed for %s" %(loc)
            raise
            
            
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
    if os.path.exists(bamfile):
        bamfile = os.path.abspath(bamfile) #re-set to include the full path to bamfile
        verboseprint("bam file is set to %s\n" %(bamfile))
    else:
        sys.stderr.write("Bam file not defined")
        raise FileNotFoundException

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
    print options.maxgenes
    if options.maxgenes is not None:
        maxgenes = int(options.maxgenes)

    
    #unwrapping the options is really unnessessary, this could be refactored to remove the unwrapping
    minreads = int(options.minreads)
    plotit = options.plotit
    poisson_cutoff = options.poisson_cutoff

    g = options.gene 
    
    gene_list = list()
    
    allpeaks = set([])
    #gets all the genes to call peaks on
    try:
        if len(g) > 0:
            gene_list = g
    except:
        gene_list = genes.keys()
        
    #Chooses between parallel and serial runs.  Potentally good to factor these two things out into different function calls
   
    bamfileobj = pysam.Samfile(bamfile, 'rb')
        
    results = []
        
    transcriptome_size = 0
    #I think this calls peaks for each gene in the gene list, which could be every gene in the genome
    running_list = []
    length_list  = []
    
    for n, gene in enumerate(gene_list):
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
        verboseprint (geneinfo+ " started:"+str(t))
        transcriptome_size += lengths[gene]
        #TODO make it so transcript size isn't always used
        #this is a filter operation, should make it as such
        running_list.append(genes[gene])
        verboseprint(lengths[gene])
        length_list.append(lengths[gene])

    #print running_list[23918]
    combined_list = zip(running_list, length_list)
  
    result = [call_peaks(gene, length, None, bamfile,  margin, options.FDR_alpha, options.threshold, 
                               minreads,  poisson_cutoff,  False,  None, 10, 1000, options.SloP, False) for gene, length in combined_list]
    print result
    jobs = [job_server.submit(call_peaks, 
                              args = (gene, length, None, bamfile,  margin, options.FDR_alpha, options.threshold, 
                               minreads,  poisson_cutoff,  False,  None, 10, 1000, options.SloP, False,), 
                              depfuncs = (peaks_from_info, get_FDR_cutoff_mean, 
                                          verboseprint,),
                              modules = ("pysam", "os", "sys", "scipy", "math", "time", "pybedtools", 
                               "random", "peaks"),) for gene, length in combined_list]

    print "looking at jobs"
    for job in jobs:
        print job.finished
        print job.tid
        job.wait()
        print job.finished
        results.append(job())   
    verboseprint("finished with calling peaks")

    #if we are going to save and output as a pickle file we should output as a pickle file
    #we should factor instead create a method or object to handle all file output
    if options.save_pickle is True:
        pickle_file = open(options.outfile + ".pickle", 'w')
        pickle.dump(results, file = pickle_file)                
    
    #combine results
    allpeaks = set([])

    #count total number of reads in transcriptiome
    transcriptome_reads = 0
    
    for gene_result in results:
        if gene_result is not None:
            verboseprint("nreads", gene_result['nreads'])
            transcriptome_reads += gene_result['nreads']
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

                    if not options.serial and transcriptomeP > poisson_cutoff:
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
                    allpeaks.add(bedline)

                except:
                    print >>sys.stderr,  "parsing failed"
                    pass
        except:
            print >>sys.stderr,  "error handling genes"
            pass
        
    #again redundant code 
    outbed = options.outfile + ".BED"
    color=options.color
    pybedtools.BedTool("\n".join(allpeaks), from_string=True).sort(stream=True).saveas(outbed, trackline="track name=\"%s\" visibility=2 colorByStrand=\"%s %s\"" %(outbed, color, color))
    print "wrote peaks to %s" %(options.outfile)
    "\n".join(allpeaks)
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
    parser.add_option("--session", dest="session", help="Type of cluster to submit multi-threaded job to, enter SGE or PBS", type="string")

    parser.add_option("--outdir", dest="prefix", default=os.getcwd(), help="output directory, default=cwd")    
    parser.add_option("--outfile", dest="outfile", default="fitted_clusters", help="a bed file output, default:%default")

    parser.add_option("--gene", "-g", dest="gene", action="append", help="A specific gene you'd like try", metavar="GENENAME")
    parser.add_option("--plot", "-p", dest="plotit", action="store_true", help="make figures of the fits", default=False)
    parser.add_option("--verbose", "-q", dest="verbose", action="store_true", help="suppress notifications")

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
    
    global varboseprint
    #creates verbose or scilent output mode
    if options.verbose:
        def verboseprint(*args):
        # Print each argument separately so caller doesn't need to
        # stuff everything to be printed into a single string
            for arg in args:
                print arg,
            print
    else:   
        verboseprint = lambda *a: None      # do-nothing function

    
    #enforces required usage    
    if not (options.bam and ((options.species) or (options.geneBEDfile and options.geneMRNAfile and options.genePREMRNAfile))):
        parser.print_help()
        exit()
        
    check_for_index(options.bam)
    
    verboseprint("Starting peak calling")        
    main(options)
    """else:

        if options.start is True:
            dtm.start(main, options)
        else:
                
            #importing drmaa here so as not to thrash serilazation 
            
            
            scriptName = os.path.join(options.prefix, options.job_name+".runme.sh")
            runerr = os.path.join(options.prefix, options.job_name+ ".err")
            runout = os.path.join(options.prefix, options.job_name+ ".out")
            shScript = open(scriptName, 'w')
            
            #trying to make job submission as general as possible, figure out the engine to submit to
            if options.session is not None:
                sessionType = options.session
            else: #use DRMAA to figure out engine type
                try:
                    from drmaa import Session
                except e: 
                    print >>sys.stderr, e 
                    print >>sys.stderr, "you must either specify the type of cluster you are running or have DRMAA installed"
                sessionType = Session.drmsInfo
                
            
            if sessionType.startswith("GE") or sessionType.startswith("SGE"):
                shScript.write("#!/bin/bash\n#$ -N %s\n#$ -S /bin/bash\n#$ -V\n#$ -pe mpi %d\n#$ -cwd\n#$ -o %s\n#$ -e %s\n#$" %(options.job_name, options.np, runout, runerr))
                if options.notify is not None:
                    shScript.write("#$ -notify\n#$ -m abe\n#$ -M %s\n" %(options.notify))
            
                shScript.write("/opt/openmpi/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines python %s --start\n" %(" ".join(sys.argv)))

            elif sessionType.startswith("PBS"):
                nnodes = 1
                if int(options.np)%8 == 0:
                    nnodes = int(options.np)/8
                else:
                    nnodes = (int(options.np)/8)+1
                    np = nnodes*8
                    verboseprint("You should have used a number of processors that is divisible by 8.  You tried %d and I'll actually use %d." %(options.np, np))
                
                shScript.write("#!/bin/sh\n#PBS -N %s\n#PBS -o %s\n#PBS -e %s\n#PBS -V\n#PBS -S /bin/sh\n#PBS -l nodes=%d:ppn=8\n#PBS -q batch\n#PBS -l walltime=00:50:00\n" %(options.job_name, runout, runerr, nnodes))
                
                if options.notify is not None:
                    shScript.write("#PBS -m abe\n#PBS -M %s\n" %(options.notify))


                shScript.write("let np=$PBS_NUM_NODES*$PBS_NUM_PPN\ncd $PBS_O_WORKDIR\n"); 
                shScript.write("/opt/openmpi/bin/mpirun --mca btl_tcp_if_include myri0 -v -machinefile $PBS_NODEFILE -np $np python %s --start\n" %(" ".join(sys.argv)))
            shScript.close()
            call(["qsub", scriptName]) #subprocess.call
        """
                

