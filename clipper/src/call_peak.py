'''
Created on Jul 25, 2012
@author: mlovci
@author: gabrielp
'''

from numpy import Inf
import sys
import pysam
from clipper.src.peaks import readsToWiggle_pysam, shuffle, find_sections
from scipy import optimize
from collections import namedtuple
#import pylab

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from numpy import diff, sign, append, array, arange, r_, empty, argmin

from math import sqrt
from scipy import interpolate
from matplotlib.path import Path

from scipy import stats
from random import sample as rs
import math
import logging
#pylab.rcParams['interactive']=True




class Peak(namedtuple('Peak', ['chrom', 
                               'genomic_start', 
                               'genomic_stop', 
                               'gene_name', 
                               'super_local_poisson_p', 
                               'strand',
                               'thick_start',
                               'thick_stop',
                               'peak_number',
                               'number_reads_in_peak',
                               'gene_poisson_p',
                               'size',
                               'p'])):
    pass
    
def get_FDR_cutoff_mode(readlengths, 
                        genelength, 
                        iterations=1000, 
                        mincut=2, 
                        alpha=.05):
   
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
            logging.info("Couldn't open a process for thresholding, trying again")
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

def get_FDR_cutoff_mean(readlengths, 
                        genelength, 
                        iterations=100, 
                        mincut=2, 
                        alpha=0.05):
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
    
    #if you have very few reads on a gene, don't waste time 
    #trying to find a cutoff        
    if len(readlengths) < 20:
        return mincut
    results = shuffle(genelength, iterations, 0, .05, readlengths) 
    total = 0

    
    #parses results from peaks script, calculates mean from peaks results 
    #should document peaks function call return value somewhere around here
        
    for cut, n_observed in enumerate(results):
        total += (cut * n_observed)
        
    #logic for min cutoffs 
    cutoff = total / iterations
    if cutoff < mincut:
        cutoff = mincut
    return int(round(cutoff, 0))

def count_turns(spline):
    
    """
    
    NOT USED (useful function though so I'll keep it around)
    
    """
    func = spline(spline._data[0])
    turns = sum(abs(diff(sign(diff(func))))) / 2
    return turns
    
class PeakGenerator(object):
    
    """
    
    An abstract class used to encapsulate all potental peak calling methods that
    we have.  New peak calling algorithms should inheret from this class
    
    """

    def __init__(self, xRange, yData):
        
        """
        
        All basic peak calling algorithms need a wiggle track and in the form of the range
        of the data, and the value at each location
        
        """
        
        self.xRange = xRange
        self.yData  = yData
    
    def peaks(self, threshold, plotit):
        
        """
        
        Idenitifes peaks given the constructed object
        
        threshold is the minimum threshold to report a peak at
        plotit plots results
        
        function returns 
        fit_values: ??
        starts_and_stops: a list of tuples detailing the start and stop of each peak
        starts: a list of all the starts
        stops: a list of all the stops
        
        
        """
        
        raise("Error abstract class, peaks not implemented")
    
class SmoothingSpline(PeakGenerator):
    """Class to fit data to a smooth curve"""
    
    def __init__(self, xRange, yData, smoothingFactor=None,
                 lossFunction = "get_turn_penalized_residuals"):
        
        """
        
        xRange -- the range to interpolate the spline over, must be monotonically increasing
        yData  -- the yAxis of the spline that corosponds to the xRange
        smoothingFactor -- tradeoff between smoothness of the spline and how well it fits
        lossFunction -- loss function to use to optomize the spline
        
        """
        
        super(SmoothingSpline,self).__init__(xRange, yData)
        
        if smoothingFactor is None:
            smoothingFactor = len(xRange)
        
        self.k = 3 #degree of spline (cubic)
        self.smoothingFactor = smoothingFactor
        self.spline = None

        #Sets loss function
        if lossFunction == "get_turn_penalized_residuals":
            self.lossFunction = self.get_turn_penalized_residuals
        elif lossFunction == "get_norm_penalized_residuals":
            self.lossFunction = self.get_norm_penalized_residuals
        else:
            raise TypeError("loss function not implemented")

    def get_norm_penalized_residuals(self, spline, norm_weight = 10, residual_weight = 1):
        
        """
        
        Returns an error value for the spline.  IN this case the error is calculated by
        a weighted combination of the norm and the residuals
        
        spline -- the smoothing spline to get the weight of
        norm_weight --the weight to apply to the norm
        residual_weight -- the weight to apply to the residuals
        
        """
        
        from scipy.linalg import norm

        #the exponent is a magic number and subject to change
        err = (norm_weight*norm(spline(self.xRange))**5) + (residual_weight*sqrt(spline.get_residual()))
        return err

    def get_turn_penalized_residuals(self, spline):
    
        """
        
        Returns an error value for the spline.  IN this case the error is calculated by
        by the numbers of turns 
        
        spline -- the smoothing spline to get the weight of
        
        """
        
        func = spline(self.xRange)

        turns = sum(abs(diff(sign(diff(func))))) / 2
 
        err = sqrt((spline.get_residual()) * (turns ** 4))

        return err
    
    def fit_univariate_spline(self, smoothingFactor = None, weight = None):
        
        """
        
        fit a spline, return the spline.
             
        (wrapper for UnivariateSpline with error handling and logging)
        
        Parameters:
        smoothingFactor -- parameter for UnivariateSpline
        xRange -- Int, range of spline
        yData -- list, wiggle track positions   
        k -- int, degree of spline
        weight -- spline weight
        
        Output: spline object
        
        """
        
        if smoothingFactor is None:
            smoothingFactor = self.smoothingFactor
                  
        try:
            spline = interpolate.UnivariateSpline(self.xRange, 
                                                  self.yData, 
                                                  s=smoothingFactor, 
                                                  k=self.k, 
                                                  w=weight)
            
        except Exception as error: #This error shouldn't happen anymore
            logging.error("failed to build spline %s, %s, %s, %s, %s, %s" % (error, 
                                                                             self.xRange, 
                                                                             self.yData, 
                                                                             smoothingFactor, 
                                                                             self.k, 
                                                                             weight) )
            raise

        return spline

    def fit_loss(self, smoothingFactor=None, weight = None):
        """fit a curve with a given smoothing parameter, return the result of the loss fxn"""

        if smoothingFactor == None:
            smoothingFactor = self.smoothingFactor

        spline = self.fit_univariate_spline(smoothingFactor=smoothingFactor, weight=weight)
        err = self.lossFunction(spline)
        
        return err
    
    def plot_spline(self, spline, title=None, threshold=None, fig=None, label = "_nolegend_"):
        
        """
        
        plots spline information
        
        spline - spline from scipy
        data - wiggle track to plot
        xvals - where to plot
        threshold - line to draw so peaks don't go below threshold
        
        """

        ax = plt.gca()
        ax.plot(spline._data[0], spline._data[1], "blue", label="_nolegend")
        ax.plot(spline._data[0], spline(spline._data[0]), label=label)
    
        if threshold is not None:
            ax.axhline(y=threshold)
        
    def plot(self, threshold=None, title=None, label="_nolegend_"):
        """plot data and spline"""
        
        spline = self.spline
        
        self.plot_spline(spline, threshold=threshold, title=title, label=label)

    def optimize_fit(self, s_estimate=None, method = 'L-BFGS-B', bounds=((1,None),),
                     weight=None):
        """

        optimize the smoothingFactor for fitting.
        
        """
        import scipy
        from scipy import optimize
        if s_estimate == None:
            s_estimate = self.smoothingFactor

        minOpts = {'disp':False,
                   'maxiter':1000}

        minimizeResult = scipy.optimize.minimize(self.fit_loss, s_estimate,
                                          #args = (weight),
                                          options = minOpts,
                                          method = method,
                                          bounds = bounds,
                                          )

        if minimizeResult.success:
            optimizedSmoothingFactor = minimizeResult.x
        
        else:
            
            #if optimization fails then we revert back to the estimate, probably should log this
            optimizedSmoothingFactor = s_estimate
          
            #logging.error("Problem spline fitting. Here is the message:\n%s" % (minimizeResult.message))
            #raise Exception
        
        optimizedSpline = self.fit_univariate_spline(optimizedSmoothingFactor, weight)
        smoothingFactor = optimizedSmoothingFactor
        self.spline = optimizedSpline
        
        return optimizedSpline 

    def get_regions_above_threshold(self, threshold, values):
        
        """
    
        Idea here is to call all regions above a given threshold and return start 
        stop pairs for those regions added twist is that when everthere is a local
        minima above the threshold we will treat that as a breakpoint
        
        generates start and stop positions for calling peaks on.
        
        threshold -- threshold for what is siginifant peak
        values -- the values (as a numpy array) arranged from 0-length of the section
        
        returns list of tuples(start, stop) used for calling peaks
        
        """
        
        xlocs = arange(0, len(values))
        
        #finds all turns, between above and below threshold
        #and generate areas to call peaks in, also 
        #makes sure starting and stopping above maxima is caught
        #threshold is at or equal to values, need to correct this
        starts = xlocs[r_[True, diff(values >= threshold)] & (values >= threshold)]
        stops = xlocs[r_[diff(values >= threshold), True] & (values >= threshold)]
        stops = stops + 1 #add to fix off by one bug
        

        #error correction incase my logic is wrong here, assuming that starts
        #and stops are always paired, and the only two cases of not being 
        #pared are if the spline starts above the cutoff or the spline starts
        #below the cutoff
        assert len(starts) == len(stops)
        
        ### important note: for getting values x->y [inclusive] 
        #you must index an array as ar[x:(y+1)]|                    
        # or else you end up with one-too-few values, the second 
        #index is non-inclusive
        
        #gets all local minima, function taken from:
        #http://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array
        #Can't have local minima at start or end, that would get caught by 
        #previous check, really need to think about that more

        local_minima = self.find_local_minima(values)

        #append to list any local minima above threshold
        for i, minima in enumerate(local_minima):
            if minima and values[i] >= threshold:
                starts = append(starts, i)
                stops = append(stops, i)
        
        starts = array(sorted(set(starts)))
        stops = array(sorted(set(stops)))
        starts_and_stops = []
        
        #making sure we aren't in some strange state
        assert len(starts) == len(stops)
        
        #get all contigous start and stops pairs         
        while len(starts) > 0:
            stop_list = stops[stops > starts[0]]
            
            #if there are no more stops left exit the loop and return the 
            #currently found starts and stops
            if len(stop_list) == 0:
                break 
            stop = stop_list[0]
            starts_and_stops.append((starts[0], stop))
            starts = starts[starts >= stop]
        
        starts = array([x[0] for x in starts_and_stops])
        stops  = array([x[1] for x in starts_and_stops])
        return starts_and_stops, starts, stops
    
    def find_local_minima(self, arr):
        
        """
        
        Returns a list of boolean values for an array that mark if a value is a local 
        minima or not True for yes false for no
        
        Importantly for ranges of local minima the value in the middle of the range
        is chosen as the minimum value
        
        """
        
        #walks through array, finding local minima ranges
        
        #hacky way to initalize a new array to all false
        minima = (arr == -1)
        min_range_start = 0
        decreasing = False
        for i in range(len(arr[:-1])):
            
            #array needs to be smooth for this to work, otherwise we'll
            #run into odd edge cases
            #update location of minima start until 
            if arr[i] > arr[i + 1]:
                min_range_start = i + 1
                decreasing = True
            
            if (arr[i] < arr[i+1]) and decreasing is True:
                decreasing = False
                #gets the local minima midpoint
                minima[(min_range_start + i) / 2] = True
        
        return minima
        
    def peaks(self, threshold=0, plotit = False):
        
        """
        
        run optimization on spline fitting.
        return peak start/stops
        
        """

        #step 1, identify good initial value
        initial_smoothing_value = self.smoothingFactor
        bestSmoothingEstimate = initial_smoothing_value


        #step 1 naive spline

        spline = self.fit_univariate_spline()
        self.spline = spline
        
        if plotit == True:
            self.plot()
            

            
        #step 2, refine to avoid local minima later
        #high-temp optimize

        best_error = self.lossFunction(spline)


            
        for i in range(2, 50):   #tries find optimal initial smooting parameter in this loop

            
            cur_smoothing_value = initial_smoothing_value * i

            cur_error = self.fit_loss(cur_smoothing_value)
            self.spline = self.fit_univariate_spline(cur_smoothing_value)
            
            if plotit == True:
                self.plot(label=str(cur_smoothing_value))

   
            if cur_error < best_error:
                bestSmoothingEstimate = cur_smoothing_value
                best_error = cur_error

        try:
            #fine optimization of smooting paramater
            #low-temp optimize
            optimizedSpline = self.optimize_fit(s_estimate=bestSmoothingEstimate)
            self.spline = optimizedSpline
            if plotit:
                self.plot(title = "optimized spline", threshold=threshold)

        except Exception as error:
            logging.error("failed spline fitting optimization at section (major crash)")

            raise

        #descretizes the data so it is easy to get regions above a given threshold
        spline_values = array([int(x) for x in optimizedSpline(self.xRange)])
    
        if plotit is True:

            self.plot(title="A fit", threshold=threshold)       
        
        starts_and_stops, starts, stops = self.get_regions_above_threshold(threshold, 
                                                                      spline_values)
    
        return (spline_values, starts_and_stops, starts, stops)

class GaussMix(object):
    from sklearn import mixture as mix
    
    def __init__(self, xvals, data):
        #data should count one "base" per read (not one point for each position)
        #initialize the fitter
        #transform data
        Tdata = list()
        for x, n in zip(xvals, data):
            for i in range(n):
                Tdata.append(x)
        self.data = Tdata #transformed to x-value frequency

    def fit(self):
        #try multiple numbers of components
        tryThisMany = 20
        models = [None for i in range(1,tryThisMany)]

        for nComponents in range(1, tryThisMany):
            #test fits for many possible numbers of components
            models[nComponents] = mix.GMM(nComponents, covariance_type='full').fit(self.data)
        #BIC = [m.bic(d) for m in models]
        AIC = [m.aic(d) for m in models]
        best = np.argmin(AIC)
        self.nComponents = best + 1
        self.GMM = models[best]
        self.AIC = AIC[best]
        
    def predict(self):
        pass

def find_local_maxima(arr):
    
    """
        
        Returns a list of boolean values for an array that mark if a value is a local 
    maxima or not True for yes false for no
    
    Importantly for ranges of local maxima the value in the middle of the range
    is chosen as the minimum value
    
    """
    
    #walks through array, finding local maxima ranges
    
    #to initalize a new array to all false
    maxima = empty(len(arr), dtype='bool')
    maxima.fill(False)
    
    max_range_start = 0
    increasing = True
    for i in range(len(arr[:-1])):
        
        #update location of maxima start until 
        if arr[i] < arr[i + 1]:
    
            max_range_start = i + 1
            increasing = True
        
        if (arr[i] > arr[i+1]) and increasing is True:
            increasing = False
            #gets the local maxima midpoint
            maxima[(max_range_start + i) / 2] = True
    
    #catches last case
    if increasing: 
        maxima[(max_range_start + len(arr) - 1) / 2] = True
        
    return maxima

def plot_sections(wiggle, sections, threshold):
    
    """
    
    Plots each section individually, I think
    Wiggle is a list representing a wiggle track
    sections is a list of strings of format "start|stop" where start and stop are both integers
    threshold is an integer 
    
    """
    
    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.plot(wiggle)
    axis.axhline(y=threshold)
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
        axis.add_patch(patch)
    plt.show()

def negative_binomial(reads_in_gene, reads_in_peak, gene_length, peak_length):
    
    """
    Paramaters
    ----------
    reads_in_gene: Integer representing number of reads in gene
    reads_in_peak: Integer reperesnting the number of reads in a specific peak
    gene_length: Integer representing length of gene
    peak_length: Integer representing length of peak
    
    Returns double, the p-value that the peak is significant
    If calcluation fails returns 1
    
    """
    
    #lambda
    lam = (float(reads_in_gene) / (gene_length)) * (peak_length)

    if lam < 3:
        lam = 3
        
    p = (lam) / (reads_in_peak + lam)
    stats.nbinom.cdf(reads_in_peak, p)
    
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
    
    try:
        #lam is estimate of the lambda value
        #poission takes a value and the lambda 
        #this is average number of reads per single 
        #site in the gene, but
        #a peak is not a single site, so it the average number 
        #gets multipled by the peak 
        #length as an estimator of the mean
        
        lam = 1 + ((float(reads_in_gene) / (gene_length)) * (peak_length)) #expect at least one read.

        cum_p = 1 - stats.poisson.cdf(reads_in_peak, int(lam))

        return cum_p
    
    except Exception as error:
        print error
        return 1

def call_peaks(loc, gene_length, bam_fileobj=None, bam_file=None, 
               margin=25, fdr_alpha=0.05, user_threshold=None,
               minreads=20, poisson_cutoff=0.05, 
               plotit=False, w_cutoff=10, windowsize=1000, 
               SloP=False, correct_p=False):
    
    """

    calls peaks for an individual gene 
    
    loc - string of all gene location
    gene_length - effective length of gene
    takes bam file or bam file object.  Serial uses object parallel uses location (name)
    margin - space between sections for calling new peaks
    fdr_alpha - false discovery rate, p-value bonferoni correct from peaks script (called in setup)
    user_threshold - user defined FDR thershold (probably should be factored into fdr_alpha

    minreads - min reads in section to try and call peaks
    poisson_cutoff - p-value for signifance cut off for number of reads in peak that gets called - might want to use ashifted distribution
    plotit - makes figures 
    
    w_cutoff - width cutoff, peaks narrower than this are discarted 
    windowssize - for super local calculation distance left and right to look 
    SloP - super local p-value instead of gene-wide p-value
    correct_p - boolean bonferoni correction of p-values from poisson
    
    """
    
    #sys.stderr.write("plotit foo" + str(plotit))
    if plotit:
        plt.rcParams['interactive']=True
        pass
    #setup
    chrom, gene_name, tx_start, tx_end, strand = loc
    logging.error("running on gene %s" % (loc))
    #logic reading bam files
    if bam_file is None and bam_fileobj is None:
        #using a file opbject is faster for serial processing 
        #but doesn't work in parallel
        
        logging.error("""you have to pick either bam file or bam file 
                        object, not both""")
        exit()
    elif bam_fileobj is None:
        bam_fileobj = pysam.Samfile(bam_file, 'rb')
        
    tx_start, tx_end = [int(x) for x in [tx_start, tx_end]]
    subset_reads = bam_fileobj.fetch(reference=chrom, start=tx_start, end=tx_end)

    #need to document reads to wiggle
    wiggle, jxns, pos_counts, lengths, allreads = readsToWiggle_pysam(subset_reads, tx_start, tx_end, strand, "center", False)

    #wiggle, pos_counts, lengths = readsToWiggle_pysam(subset_reads, tx_start, tx_end, strand, "center", False)

    #TODO have a check to kill this if there aren't any reads in a region
        
    result = peaks_from_info(list(wiggle), pos_counts, lengths, loc, gene_length, margin, fdr_alpha, user_threshold, minreads, poisson_cutoff, plotit, w_cutoff, windowsize, SloP, correct_p)
    
    return result

def peaks_from_info(wiggle, pos_counts, lengths, loc, gene_length, 
                    margin=25, fdr_alpha=0.05, user_threshold=None,
                    minreads=20, poisson_cutoff=0.05, plotit=False, 
                    width_cutoff=10, windowsize=1000, SloP=False, 
                    correct_p=False):

    """
    
    same args as before 
    wiggle is converted from bam file
    pos_counts - one point per read instead of coverage of entire read
    lengths - lengths aligned portions of reads 
    rest are the same fix later


    calls peaks for an individual gene 
    

    gene_length - effective length of gene
    margin - space between sections for calling new peaks
    fdr_alpha - false discovery rate, p-value bonferoni correct from peaks script (called in setup)
    user_threshold - user defined FDR thershold (probably should be factored into fdr_alpha
    minreads - min reads in section to try and call peaks
    poisson_cutoff - p-value for signifance cut off for number of reads in peak that gets called - might want to use ashifted distribution
    plotit - makes figures 
    
    w_cutoff - width cutoff, peaks narrower than this are discarted 
    windowssize - for super local calculation distance left and right to look 
    SloP - super local p-value instead of gene-wide p-value
    correct_p - boolean bonferoni correction of p-values from poisson
        
    """

    peak_dict = {}
    
    #all the information nessessary to record a peak, used later, but declared outside of loops

    
    #these are what is built in this dict, complicated enough that it might 
    #be worth turning into an object
    #peak_dict['clusters'] = {}
    #peak_dict['sections'] = {}
    #peak_dict['nreads'] = int()
    #peak_dict['threshold'] = int()
    #peak_dict['loc'] = loc
    
    #data munging
    chrom, gene_name, tx_start, tx_end, strand = loc
    tx_start, tx_end = [int(x) for x in [tx_start, tx_end]]    
    
    #used for poisson calclulation? 
    nreads_in_gene = sum(pos_counts)
    
    #decides FDR calcalation, maybe move getFRDcutoff mean into c code
    gene_threshold = 0
    if user_threshold is None:
        gene_threshold = get_FDR_cutoff_mean(lengths, 
                                             gene_length, 
                                             alpha=fdr_alpha)
    else:
        logging.info("using user threshold")
        gene_threshold = user_threshold
    
    if not isinstance(gene_threshold, int):
        raise TypeError
        
    peak_dict['clusters'] = []
    peak_dict['sections'] = {}
    peak_dict['nreads'] = int(nreads_in_gene)
    peak_dict['threshold'] = gene_threshold
    peak_dict['loc'] = loc
    peak_number=1

 
    sections = find_sections(wiggle, margin)
    if plotit is True:      
        plot_sections(wiggle, sections, gene_threshold)

    for sect in sections:
        sectstart, sectstop = sect
        sect_length = sectstop - sectstart + 1
        data = wiggle[sectstart:(sectstop + 1)]
        cts = pos_counts[sectstart:(sectstop + 1)]
        xvals = arange(0, sect_length)
        Nreads = sum(cts)

        peak_dict['sections'][sect] = {}
        threshold = int()
        peak_dict['sections'][sect]['nreads'] = int(Nreads)

        #makes sure there are enough reads
        if Nreads < minreads:
            logging.info("""%d is not enough reads, skipping section: %s""" %(Nreads, sect))
            peak_dict['sections'][sect]['tried'] = False            
            continue
        else:
            logging.info("""Analyzing section %s with %d reads""" %(sect, Nreads))
            pass
        
            
        if user_threshold == None:
            if SloP:
                
                #gets random subset of lengths of reads for calculations on a section
                #not exactly the right way to do this but it should be very close.
                sect_read_lengths = rs(lengths, Nreads) 
                
                #use the minimum FDR cutoff between superlocal and gene-wide calculations
                threshold = min(gene_threshold, get_FDR_cutoff_mean(sect_read_lengths, 
                                                sect_length, 
                                                alpha=fdr_alpha))
                logging.info("Using super-local threshold %d" %(threshold))
                
            else:
                threshold = gene_threshold
        else:
            threshold = user_threshold

        #saves threshold for each individual section
        peak_dict['sections'][sect]['threshold'] = threshold
        peak_dict['sections'][sect]['nreads'] = int(Nreads)
        peak_dict['sections'][sect]['tried'] = True
        peak_dict['sections'][sect]['nPeaks'] = 0
        #if wiggle track never excides threshold
        if max(data) < threshold:
            logging.info("data does not excede threshold, stopping")
            continue

        fitType = "Spline"
        ###gauss mixture model
        if fitType == "Spline":
            
            initial_smoothing_value = (sectstop - sectstart + 1)
            print xvals, data, initial_smoothing_value
            fitter = SmoothingSpline(xvals, data, initial_smoothing_value,
                            lossFunction="get_norm_penalized_residuals")
        
        elif fitType == "Gaussian":
            fitter = GaussMix(xvals, data)
            

        try:
            (fit_values, starts_and_stops, starts, stops) = fitter.peaks(threshold, plotit)
        except Exception as error:
            print gene_name
            raise error
            
         
        #walks along spline, and calls peaks along spline
        #for each start, take the next stop and find the peak 
        #between the start and the stop this is where I need to 
        #fix, some peaks starts start right after another start, 
        #but not on top of it make sure the next start is after the 
        #previous stop

        #subsections that are above threshold
        for p_start, p_stop in starts_and_stops: 

            #peaks with-in this subsection, indexed from section 
            #(not subsection) start
            #find all local maxima
            peaks = [x + p_start for x in xvals[find_local_maxima(fit_values[p_start:(p_stop + 1)])]]

            assert len(peaks) in (0,1) 

            #handles logic if there are multiple peaks between 
            #start and stop
            if len(peaks) <= 0:
                continue
            if len(peaks) is 1:
                #TODO All this formatting logic doesn't belong here 
                #should be simplifed
                #gets reads in peak
                number_reads_in_peak = sum(cts[p_start:(p_stop + 1)])
                logging.info("""Peak %d (%d - %d) has %d 
                                 reads""" %(peak_number,                                             
                                             p_start,
                                             (p_stop + 1),
                                             number_reads_in_peak))

                #makes sure there enough reads
                if (number_reads_in_peak < minreads or 
                    max(data[p_start:(p_stop + 1)]) < threshold):
                    logging.info("""skipping peak, %d is not enough reads"""
                                 %(number_reads_in_peak))
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

                #best_error checking logic to keep bed files from breaking
                if thick_start < g_start:
                    thick_start = g_start
                if thick_stop > g_stop:
                    thick_stop = g_stop

                peak_length = g_stop - g_start + 1

                #skip really small peaks
                if peak_length < width_cutoff:
                    continue
              

                #super local logic 
                #best_error check to make sure area is in area of gene

                #distance from gene start
                if peak - tx_start - windowsize < 0: 
                    area_start = 0

                #for super local gets area around peak for calculation
                else:  
                    area_start = peak - tx_start - windowsize
                    #area_start = sectstart

                #same thing except for end of gene instead of start
                if peak + windowsize > tx_end: #distance to gene stop
                    area_stop = tx_start - tx_end + 1
                else:
                    area_stop = peak - tx_start + windowsize
                    #area_stop = sectstop

                #use area reads + 1/2 all other reads in gene: 
                #area_reads = sum(pos_counts[area_start:area_stop]) + 
                #0.5*(sum(pos_counts) - 
                #sum(pos_counts[area_start:area_stop]))

                #use area reads:
                area_reads = sum(pos_counts[area_start:area_stop])
                area_size = area_stop - area_start + 1

                #area_reads = sum(pos_counts[sectstart:sectstop])
                #area_size = sect_length

                #calcluates poisson based of whole gene vs peak
                gene_pois_p = poissonP(nreads_in_gene, 
                                       number_reads_in_peak, 
                                       gene_length, 
                                       peak_length)
                if SloP is True:
                    #same thing except for based on super local p-value
                    slop_pois_p = poissonP(area_reads, 
                                          number_reads_in_peak, 
                                          area_size, 
                                          peak_length)

                #makes sure spop_poisP is defined, even if its 
                #just normal, something to be removed later,
                #slop should only be used when defined as true
                else:
                    slop_pois_p = gene_pois_p


                if math.isnan(slop_pois_p):
                    slop_pois_p = 1

                #defines the bedline of a peak for returning
                #TODO This should be abstracted out for now... seperate model from view
                
                peak_dict['clusters'].append(Peak(chrom, 
                                                  g_start, 
                                                  g_stop, 
                                                  gene_name, 
                                                  slop_pois_p, 
                                                  strand,
                                                  thick_start,
                                                  thick_stop,
                                                  peak_number,
                                                  number_reads_in_peak,
                                                  gene_pois_p,
                                                  peak_length,
                                                  0
                                                  )
                                             )
   
                peak_number += 1
                peak_dict['sections'][sect]['nPeaks'] +=1
           
    #inflate p-values based on # of comparisons #bonferroni corrected
    if correct_p is True:
        #best I can tell this never executes...            
        for peak in peak_dict['clusters']:
            peak.p = peak.p * peak_number  #bonferroni correct p-value for MHT
        
        

    peak_dict['Nclusters'] = peak_number
    if plotit:
        import sys
        plt.show()
        v = sys.stdin.read(1)
    return peak_dict
