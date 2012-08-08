'''
Created on Jul 25, 2012

@author: gabrielp
'''

import peaks
from numpy import Inf
import sys
import pysam

def verboseprint(*args):
        # Print each argument separately so caller doesn't need to
        # stuff everything to be printed into a single string
            for arg in args:
                print arg,
            print
            
def get_FDR_cutoff_mode(readlengths, genelength, iterations=1000, mincut=2, alpha=.05):
    """
    
    Find randomized method, as in FOX2ES NSMB paper.
    
    """
    if readlengths.__len__() < 20: # if you have very few reads on a gene, don't waste time trying to find a cutoff
        return mincut
    cmd = "./peaks"
    bad = 1
    tries = 0
    while bad == 1 and tries < 5:
        try:
            process = Popen([cmd, "-f", "stdin", "-L", str(genelength), "-r", str(iterations), "-a", str(alpha)], stdin=PIPE, stdout=PIPE)
            results, err = process.communicate("\n".join(map(str, readlengths)))
            return_val = process.wait()
            bad = 0
        except OSError:
            print >> sys.stderr, "Couldn't open a process for thresholding, trying again"
            tries += 1
        
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

def get_FDR_cutoff_mean(readlengths, genelength, iterations=100, mincut=2, alpha=0.05):
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
    n = 0
    
    #parses results from peaks script, calculates mean from peaks results 
    #should document peaks function call return value somewhere around here
        
    for cut, n_observed in enumerate(results):
        total += (cut * n_observed)
        n += n_observed
        
    #logic for min cutoffs 
    cutoff = total / iterations
    if cutoff < mincut:
        cutoff = mincut
    return int(round(cutoff, 0))



def plotSpline(spline, data, xvals, section, threshold=None):
    """
    
    plots spline information
    
    spline - spline from scipy
    data - wiggle track to plot
    xvals - where to plot
    threshold - line to draw so peaks don't go below threshold
    
    """
    import matplotlib.pyplot as plt 
    f = plt.figure()
    plt.title(section)
    ax1 = f.add_subplot(111)
    ax1.plot(xvals, data)
    ax1.plot(xvals, spline(xvals))
    if threshold is not None:
        ax1.axhline(y=threshold)
    plt.show()
    



def find_splineResiduals(x, xdata, ydata, k, weight=None):
    
    """
    
    Returns complexity penalized residuals from a smoothing spline
    
    x -- Int, range of spline
    xdata -- list, wiggle track positions   
    ydata -- list, wiggle track coverage 
    k -- int, degree of spline
    
    """
    #imports nessessary things
    from numpy import diff, sign
    from math import sqrt
    
    spline = find_univariateSpline(x, xdata, ydata, k, weight=None)
    
    #catches spline fitting error
    if spline is None:
        return Inf
    
    #calculates residuals
    func = spline(xdata)
    turns = sum(abs(diff(sign(diff(func))))) / 2 # number of turns in the function
    err = sqrt((spline.get_residual()) * (turns ** 4))
    return(err)
    
def find_univariateSpline(x, xdata, ydata, k, weight=None):
    
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
    from scipy import interpolate
    
    try:
        print xdata, ydata, x, k, weight
        spline = interpolate.UnivariateSpline(xdata, ydata, s=x, k=k, w=weight)
        return(spline)
    
    except Exception as e:
        verboseprint("failed to build spline", e)
        return None


def plotSections(wiggle, sections, threshold):
    
    """
    
    Plots each section individually, I think
    Wiggle is a list representing a wiggle track
    sections is a list of strings of format "start|stop" where start and stop are both integers
    threshold is an integer 
    
    """
    
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches
    
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(wiggle)
    ax.axhline(y=threshold)
    for sect in sections:
        #mark active sections
        positions = list() 
        codes = list()
        start, stop = sect
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


def poissonP(reads_in_gene, reads_in_peak, gene_length, peak_length):
    
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
    
    from scipy import stats
    try:
        #lam is estimate of the lambda value
        #poission takes a value and the lambda 
        #this is average number of reads per single site in the gene, but
        #a peak is not a single site, so it the average number gets multipled by the peak 
        #length as an estimator of the mean
        
        #TODO: check with boyko or someone else about this math.  I think its a strong over 
        #estimate of the mean
        lam = (float(reads_in_gene) / (gene_length)) * (peak_length)

        if lam < 3:
            lam = 3;
        cumP = 1 - stats.poisson.cdf(reads_in_peak, int(lam))
        return cumP
    except Exception as e:
        print e
        return 1


def call_peaks(loc, gene_length, bam_fileobj=None, bam_file=None, margin=25, FDR_alpha=0.05, user_threshold=None,
               minreads=20, poisson_cutoff=0.05, plotit=False, w_cutoff=10, windowsize=1000, SloP=False, correct_P=False):
    
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
    
    w_cutoff - width cutoff, peaks narrower than this are discarted 
    windowssize - for super local calculation distance left and right to look 
    SloP - super local p-value instead of gene-wide p-value
    correct_P - boolean bonferoni correction of p-values from poisson
    
    """
    
    #setup
    chrom, gene_name, tx_start, tx_end, signstrand = loc.split("|")

    
    #logic reading bam files
    if bam_file is None and bam_fileobj is None:
        #using a file opbject is faster for serial processing bot doesn't work in parallel
        verboseprint("you have to pick either bam file or bam file object, not both")
        exit()
    elif bam_fileobj is None:
        bam_fileobj = pysam.Samfile(bam_file, 'rb')
        
    tx_start, tx_end = map(int, [tx_start, tx_end])
    subset_reads = bam_fileobj.fetch(reference=chrom, start=tx_start, end=tx_end)

    #need to document reads to wiggle
    #wiggle, pos_counts, lengths = readsToWiggle_pysam_foo(subset_reads, tx_start, tx_end, signstrand, "center")
    wiggle, jxns, pos_counts, lengths, allreads = peaks.readsToWiggle_pysam(subset_reads, tx_start, tx_end, signstrand, "center", False)
    r = peaks_from_info(list(wiggle), pos_counts, lengths, loc, gene_length, margin, FDR_alpha, user_threshold, minreads, poisson_cutoff, plotit, w_cutoff, windowsize, SloP, correct_P)

    return r


def get_start_stop_pairs_above_threshold(threshold, values):
    
    """

    Idea here is to call all regions above a given threshold and return start stop pairs for those regions
    added twist is that when everthere is a local minima above the threshold we will treat that as a breakpoint
    
    generates start and stop positions for calling peaks on.  Helper function that was abstracted 
    from peaks_from_info
    
    threshold -- threshold for what is siginifant peak
    values -- the values (as a numpy array) arranged from 0-length of the section
    sect_length -- the length of the section that we will attempt to call peaks on
    xvals -- the location of all the xvalues we are looking at
    spline -- spline object from find_univarateSpline
    
    returns list of tuples(start, stop) used for calling peaks
    
    """
    
    from numpy import diff, sign, append, insert, array, arange, r_
    
    xlocs = arange(0, len(values))
    #finds all turns, and generate areas to call peaks in, also makes sure starting and stopping above 
    #maxima is caught
    starts = xlocs[r_[True, diff(values >= threshold)] & (values >= threshold)]
    stops = xlocs[r_[diff(values >= threshold), True] & (values >= threshold)]
    
    #error correction incase my logic is wrong here, assuming that starts and stops
    #are always paired, and the only two cases of not being pared are if the spline starts above
    #the cutoff or the spline starts below the cutoff
    assert len(starts) == len(stops)
    
    ### important note: for getting values x->y [inclusive] you must index an array as ar[x:(y+1)]|                     or else you end up with one-too-few values, the second index is non-inclusive
    
    #gets all local minima, function taken from:
    #http://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array
    #Can't have local minima at start or end, that would get caught by previous check, really need to think
    #about that more
    local_minima = r_[False, values[1:] < values[:-1]] & r_[values[:-1] < values[1:], False]
    
    #append to list any local minima above threshold
    for i, minima in enumerate(local_minima):
        if minima and values[i] >= threshold:
            starts = append(starts, i)
            stops = append(stops, i)
    
    starts = array(sorted(set(starts)))
    stops = array(sorted(set(stops)))
    starts_and_stops = []

    #get all contigous start and stops pairs
    while len(starts) > 0:
        stop_list = stops[stops > starts[0]]
        
        #if there are no more stops left exit the loop and return the currently found starts and stops
        if len(stop_list) == 0:
            break 
        stop = stop_list[0]
        starts_and_stops.append((starts[0], stop))
        starts = starts[starts >= stop]
    
    starts = array(map(lambda x: x[0], starts_and_stops))
    stops = array(map(lambda x: x[1], starts_and_stops))
    return starts_and_stops, starts, stops



def peaks_from_info(wiggle, pos_counts, lengths, loc, gene_length, margin=25, FDR_alpha=0.05, user_threshold=None,
                                   minreads=20, poisson_cutoff=0.05, plotit=False, w_cutoff=10, windowsize=1000, SloP=False, correct_P=False):

    """
    
    same args as before 
    wiggle is converted from bam file
    pos_counts - one point per read instead of coverage of entire read
    lengths - lengths aligned portions of reads 
    rest are the same fix later
    
    """

    #this is how things need to be for parallization to work
    import peaks
    from numpy import arange, diff, sign, array
    from random import sample as rs
    import math
    from src.call_peak import find_splineResiduals, find_univariateSpline, poissonP, plotSections, plotSpline, get_start_stop_pairs_above_threshold
    peakDict = {}
    import scipy  
    from scipy import optimize 
    
    
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
        verboseprint("I had a hard time with this one: %s.  I think I'll use a threshold of 50" % (loc))
        threshold = 50
    peakDict['clusters'] = {}
    peakDict['sections'] = {}
    peakDict['nreads'] = int(nreads_in_gene)
    peakDict['threshold'] = gene_threshold
    peakDict['loc'] = loc
    peakn = 1
 
    verboseprint("Testing %s" % (loc))
    verboseprint("Gene threshold is: %d" % (gene_threshold))
    
    #print wiggle
    #print margin
    sections = peaks.find_sections(wiggle, margin)
    if plotit is True:      
        plotSections(wiggle, sections, gene_threshold)

    for sect in sections:
        sectstart, sectstop = sect
        sect_length = sectstop - sectstart + 1
        data = wiggle[sectstart:(sectstop + 1)]
        cts = pos_counts[sectstart:(sectstop + 1)]
        xvals = arange(0, sect_length)
        Nreads = sum(cts)

        #gets random subset of lengths of reads for calculations on a section
        sect_read_lengths = rs(lengths, Nreads) #not exactly the right way to do this but it should be very close.
        peakDict['sections'][sect] = {}
        threshold = int()
        
        #makes sure there are enough reads
        if Nreads < minreads:
            verboseprint("%d is not enough reads, skipping section: %s" % (Nreads, sect))
            continue
        else:
            verboseprint("Analyzing section %s with %d reads" % (sect, Nreads))
        
            
        #sets super-local if requested, might be able to factor this
        if user_threshold is None:
            if SloP is True:
                threshold = get_FDR_cutoff_mean(sect_read_lengths, sect_length, alpha=FDR_alpha)
                verboseprint("Using super-local threshold %d" % (threshold))
                
            else:
                threshold = gene_threshold
        else:
            threshold = user_threshold

        #saves threshold for each individual section
        peakDict['sections'][sect]['threshold'] = threshold
        peakDict['sections'][sect]['nreads'] = int(Nreads)

        #if wiggle track never excides threshold
        if max(data) < threshold:
            verboseprint("data does not excede threshold, stopping")
            continue
        
        #fitting splines logic, black magic 
        try:
            degree = 3 #cubic spline
            weights = None 
            #for very large windows with many reads a large smoothing parameter is required.  test several different options to determine a reasonable inital estimate
            
            #Goal is to find optimnal smooting paramater in multiple steps
            #x1 initial estimate of smoothing paramater 
            #step 1, identify good initial value
            x1 = (sectstop - sectstart + 1)
            x0 = x1
            useme = 1
            
           
            
            #step 2, refine so as not to runinto local minima later, try to come up with a good way of getting optimal paramater
            error = find_splineResiduals(x1, xvals, data, degree, weights)
            for i in range(2, 11):
                x2 = x1 * i
                
                #tries find optimal initial smooting paraater in this loop
                err = find_splineResiduals(x2, xvals, data, degree, weights)
                if err < error:
                    x0 = x2
                    useme = i

         
            verboseprint("I'm using (region length) * %d as the initial estimate for the smoothing parameter" % (useme))          
                        
            try:
                #fine optimization of smooting paramater
                cutoff = float(0)
                tries = 0
                while cutoff < 5:# shouldn't get smoothing coef's this small.. increase the initial estimate and try again. WARNING: BLACK MAGIC
                    tries += 1
                    if tries == 3: # increasing this may improve accuracy, but at the cost of running time.
                        break
                    
                    #TODO make verbose a global variable so I can display stuff again
                    sp = optimize.minimize(find_splineResiduals, x0, args=(xvals, data, degree, weights),
                                                 options={'disp':False}, method="Powell")
                    #fit a smoothing spline using an optimal parameter for smoothing and with weights proportional to the number of reads aligned at each position if weights is set
                    if sp.success is True:
                        cutoff = sp.x
                    else:
                        pass
                    x0 += sect_length
            except Exception as e:
                print >> sys.stderr, e
                raise
                #print >>sys.stderr,  "%s failed spline fitting at section %s with sp: %s" %(loc, sect)
                #continue

       
            verboseprint ("optimized smoothing parameter")
        #if we are going to save and output as a pickle fi is %s" %(str(cutoff))
        #final fit spline

            spline = find_univariateSpline(cutoff, xvals, data, degree, weights)
          
            if plotit is True:
                plotSpline(spline, data, xvals, peakn, threshold)
            
            starts_and_stops, starts, stops = get_start_stop_pairs_above_threshold(threshold, spline(xvals))

            
            #walks along spline, and calls peaks along spline
            #for each start, take the next stop and find the peak between the start and the stop
            #this is where I need to fix, some peaks starts start right after another start, but not on top of it
            #make sure the next start is after the previous stop
            for p_start, p_stop in starts_and_stops: #subsections that are above threshold
                try:
                    peaks = map(lambda x: x + p_start, xvals[diff(sign(diff(spline(xvals[p_start:(p_stop + 1)])))) < 0])  #peaks with-in this subsection, indexed from section (not subsection) start
                    verboseprint("I found %d peaks" % (len(peaks)))
                except:
                    continue
                
                #handles logic if there are multiple peaks between start and stop
                if len(peaks) <= 0:
                    continue
                if len(peaks) is 1:
                    #gets reads in peak
                    Nreads_in_peak = sum(cts[p_start:(p_stop + 1)])
                    verboseprint("Peak %d (%d - %d) has %d reads" % (peakn, p_start, (p_stop + 1), Nreads_in_peak))
                    
                    #makes sure there enough reads
                    if Nreads_in_peak < minreads or max(data[p_start:(p_stop + 1)]) < threshold:
                        verboseprint("skipping peak, %d is not enough reads" % (Nreads_in_peak))
                        continue

                    #formatting of bed track
                    #start and stop for bed track to be created
                    g_start = tx_start + sectstart + p_start
                    g_stop = tx_start + sectstart + p_stop
                    
                    #highest point in start stop
                    peak = tx_start + sectstart + peaks[0]
                    
                    #makes it thicker so we can see on the browser 
                    thick_start = peak - 2
                    thick_stop = peak + 2
                    
                    #error checking logic to keep bed files from breaking
                    if thick_start < g_start:
                        thick_start = g_start
                    if thick_stop > g_stop:
                        thick_stop = g_stop
                    peak_length = g_stop - g_start + 1

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
                        area_stop = tx_start - tx_end + 1
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
                        slop_poisP = gene_poisP
                    
                    
                    if math.isnan(slop_poisP):
                        slop_poisP = 1
                                            
                    #remove later    
                    if slop_poisP > poisson_cutoff:
                        #continue
                        pass
                    
                    #poisP = 1
                    
                    #defines the bedline of a peak for returning                    
                    bedline = "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d" % (chrom, g_start, g_stop, peak_name, slop_poisP, signstrand, thick_start, thick_stop)

                    #metadata for the specific bedline
                    peakDict['clusters'][bedline] = {}
                    peakDict['clusters'][bedline]['GeneP'] = gene_poisP
                    peakDict['clusters'][bedline]['SloP'] = slop_poisP                    
                    peakDict['clusters'][bedline]['Nreads'] = Nreads_in_peak
                    peakDict['clusters'][bedline]['size'] = peak_length
                    

                    peakn += 1
                else:  #there are more than one peaks in this window
                    #this handles peaks within peaks logic
                    valleys = array(map(lambda x:x + p_start, xvals[diff(sign(diff(spline(xvals[p_start:p_stop + 1])))) > 0]))#local minima in subsection, relative to section start

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
                        Nreads_in_peak = sum(cts[subpeak_start:(subpeak_stop + 1)])
                        if Nreads_in_peak < minreads or max(data[subpeak_start:(subpeak_stop + 1)]) < threshold:
                            continue
                        g_start = tx_start + subpeak_start + sectstart
                        g_stop = tx_start + subpeak_stop + sectstart
                        peak = tx_start + subpeak + sectstart
                        thick_start = peak - 2
                        if thick_start < g_start:
                            thick_start = g_start                        
                        thick_stop = peak + 2
                        if thick_stop > g_stop:
                            thick_stop = g_stop                        
                        peak_name = gene_name + "_" + str(peakn) + "_" + str(int(Nreads_in_peak))
                        if peak - tx_start - windowsize < 0: #distance from gene start
                            area_start = 0 
                        else:
                            area_start = peak - tx_start - windowsize
                            #area_start = sectstart
                        if peak + windowsize > tx_end: #distance to gene stop
                            area_stop = tx_start - tx_end + 1
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
                        bedline = "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d" % (chrom, g_start, g_stop, peak_name,
                                                                     slop_poisP, signstrand, thick_start, thick_stop)
                        peakDict['clusters'][bedline] = {}                        
                        peakDict['clusters'][bedline]['SloP'] = slop_poisP
                        peakDict['clusters'][bedline]['GeneP'] = gene_poisP
                        peakDict['clusters'][bedline]['Nreads'] = Nreads_in_peak
                        peakDict['clusters'][bedline]['size'] = peak_length
                        peakn += 1
        except NameError as e:
            print >> sys.stderr, e
            print >> sys.stderr, "spline fitting failed for %s" % (loc)
            raise
            
            
    #inflate p-values based on # of comparisons #bonferroni corrected
    if correct_P is True:            
        for peak in peakDict['clusters']:
            peakDict['clusters'][peak]['p'] = peakDict['clusters'][peak]['p'] * peakn  #bonferroni correct p-value for MHT
        
        

    peakDict['Nclusters'] = peakn
    
    return peakDict
