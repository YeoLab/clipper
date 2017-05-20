'''
Created on Jul 25, 2012
@author: mlovci
@author: gabrielp
'''

from collections import namedtuple
import logging
import math


import HTSeq
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from math import sqrt
import numpy
from numpy import diff, sign, append, array, arange, r_, empty
import numpy as np
import pysam
from scipy import stats
from scipy import interpolate
from scipy.stats import binom
import pybedtools

from clipper.src.peaks import shuffle, find_sections
from clipper.src.readsToWiggle import readsToWiggle_pysam
from clipper.src.bam_helpers import Robust_BAM_Reader

class Peak(namedtuple('Peak', ['chrom',
                               'genomic_start',
                               'genomic_stop',
                               'gene_name',
                               'strand',
                               'thick_start',
                               'thick_stop',
                               'peak_number',
                               'number_reads_in_peak',
                               'size',
                               'p',
                               'effective_length',
                               'peak_length',
                               'area_reads',
                               'area_size',
                               'nreads_in_gene'
                               #'nreads_in_input'
                               ])):
    def __repr__(self):
        """bed8 format"""
        return "\t".join(map(str, [self.chrom, self.genomic_start, self.genomic_stop,
                         "_".join(map(str, [self.gene_name, self.peak_number, self.number_reads_in_peak])),
                        self.strand,
                         self.thick_start, self.thick_stop]))

    def __len__(self):
        return self.genomic_stop - self.genomic_start
    pass

def get_FDR_cutoff_binom(readlengths, genelength, alpha, mincut = 2):
    number_reads = len(readlengths)
    
    if number_reads == 0:
        return mincut
    else:
        read_length = numpy.array(readlengths)
        mean_read_length = numpy.mean(read_length)
        prob = float(mean_read_length) / float(genelength)
        if prob > 1:
            raise ValueError("probability of >= 1 read per-base > 1")
        try:
            k = int(binom.ppf(1 - (alpha), number_reads, prob))
            if k < mincut:
                return mincut
            else:
                return k
        except:
            print read_length, mean_read_length, genelength, prob, alpha, number_reads
            raise
        
def get_FDR_cutoff_mode(readlengths, 
                        genelength, 
                        iterations=1000, 
                        mincut=2, 
                        alpha=.05):
   
    """
    
    Find randomized method, as in FOX2 ES NSMB paper.
    
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
    
    results = shuffle(int(genelength), int(iterations), 0, .05, readlengths) 
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
        
        self.xRange = np.array(xRange)
        self.yData = np.array(yData)
    
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
    
    def __init__(self, xRange, yData, smoothing_factor=None,
                 lossFunction="get_turn_penalized_residuals",
                 threshold=0,
                 num_reads=0):
        
        """
        
        xRange -- the range to interpolate the spline over, must be monotonically increasing
        yData  -- the yAxis of the spline that corosponds to the xRange
        smoothingFactor -- tradeoff between smoothness of the spline and how well it fits
        lossFunction -- loss function to use to optomize the spline
        
        """
        
        super(SmoothingSpline, self).__init__(xRange, yData)
        
        if smoothing_factor is None:
            #smoothingFactor = 0.25 * numpy.sum(yData) #
            smoothing_factor = len(xRange)

        #degree of spline (cubic)
        self.k = 3
        self.num_reads = num_reads
        self.smoothing_factor = smoothing_factor
        self.spline = None
        self.threshold = threshold
        
        #Sets loss function
        if lossFunction == "get_turn_penalized_residuals":
            self.lossFunction = self.get_turn_penalized_residuals
        elif lossFunction == "get_norm_penalized_residuals":
            self.lossFunction = self.get_norm_penalized_residuals
        else:
            raise TypeError("loss function not implemented")

    def get_norm_penalized_residuals(self, spline, norm_weight=1, residual_weight=10):

        """

        Returns an error value for the spline.  IN this case the error is calculated by
        a weighted combination of the norm and the residuals

        spline -- the smoothing spline to get the weight of
        norm_weight --the weight to apply to the norm
        residual_weight -- the weight to apply to the residuals

        """

        from scipy.linalg import norm

        #the exponent is a magic number and subject to change
        err = (norm_weight*norm(spline(self.xRange))**2) + (residual_weight*sqrt(spline.get_residual()))
        return err

    def get_turn_penalized_residuals(self, spline, residual_weight=1, turn_weight=1, turn_exp=10000):

        """

        Returns an error value for the spline.  IN this case the error is calculated by
        by the numbers of turns

        spline -- the smoothing spline to get the weight of

        """

        func = spline(self.xRange)
        turns = sum(abs(diff(sign(diff(func))))) / 2
        turn_exp = 4
        residual_weight = 1000
        err = (residual_weight * sqrt((spline.get_residual()))) * (turn_weight * (turns ** turn_exp))

        return err

    def fit_univariate_spline(self, smoothingFactor=None, weight=None):

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
            smoothingFactor = self.smoothing_factor

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

    def fit_loss(self, smoothing_factor=None, weight=None):
        """fit a curve with a given smoothing parameter, return the result of the loss fxn"""

        if smoothing_factor is None:
            smoothing_factor = self.smoothing_factor

        #print "smoothing_factor", smoothing_factor
        spline = self.fit_univariate_spline(smoothingFactor=smoothing_factor, weight=weight)
        err = self.lossFunction(spline)

        return err


    def optimize_fit(self, s_estimate=None, method='bounded', bounds=((1, None),),
                     weight=None):
        """

        optimize the smoothingFactor for fitting.

        """
        import scipy
        from scipy import optimize
        if s_estimate is None:
            s_estimate = self.smoothing_factor

        min_opts = {'disp': False,
                    'maxiter': 5000}

        def con(t):
            return t[0] > 1

        cons = {'type': "ineq",
                'fun': con}

        #I think that the smoothing bounds are linear
        bound_scale = self.num_reads * 20
        bound_scale = 100 if bound_scale <= 1 else bound_scale
        #print bound_scale
        minimize_result = scipy.optimize.minimize_scalar(self.fit_loss,
                                                          s_estimate,
                                                          #args = (weight),
                                                          options=min_opts,
                                                          method=method,
                                                          #constraints=cons,
                                                          bounds=(1, bound_scale),
                                                          )

        if minimize_result.success:
            optimized_smoothing_factor = minimize_result.x

        else:
            #if optimization fails then we revert back to the estimate, probably should log this
            optimized_smoothing_factor = s_estimate
            #print "smoothing factor", self.smoothing_factor, s_estimate
            logging.error("Problem spline fitting. Here is the message:\n%s" % (minimize_result.message))
            #raise Exception

        optimized_spline = self.fit_univariate_spline(optimized_smoothing_factor, weight)
        self.smoothing_factor = optimized_smoothing_factor
        #print "final smoothing factor", str(self.smoothing_factor), bound_scale
        self.spline = optimized_spline
        #print "optimized: %f" % optimizedSmoothingFactor
        self.result = minimize_result
        return optimized_spline

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
        #add to fix off by one bug
        stops += + 1

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
        stops = array([x[1] for x in starts_and_stops])
        return starts_and_stops, starts, stops

    def find_local_maxima(self, arr):

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

    def peaks(self, threshold=0, plotit=False):

        """

        run optimization on spline fitting.
        return peak start/stops

        """

        #step 1, identify good initial value
        initial_smoothing_value = self.smoothing_factor
        best_smoothing_estimate = initial_smoothing_value

        #step 1 naive spline
        spline = self.fit_univariate_spline()
        self.spline = spline

        if plotit:
            self.plot()

        #step 2, refine to avoid local minima later
        #high-temp optimize
        best_error = self.lossFunction(spline)

        #tries to find optimal initial smoothing parameter in this loop
        for i in range(1, 50):

            cur_smoothing_value = initial_smoothing_value * i

            cur_error = self.fit_loss(cur_smoothing_value)
            self.spline = self.fit_univariate_spline(cur_smoothing_value)

            if plotit:
                self.plot(label=str(cur_smoothing_value))

            if cur_error < best_error:
                best_smoothing_estimate = cur_smoothing_value
                best_error = cur_error


        try:
            #fine optimization of smoothing parameter
            #low-temp optimize
            optimizedSpline = self.optimize_fit(s_estimate=best_smoothing_estimate)
            self.spline = optimizedSpline
            if plotit:
                self.plot(title="optimized spline", threshold=self.threshold)

        except Exception as error:
            logging.error("failed spline fitting optimization at section (major crash)")

            raise

        #descretizes the data so it is easy to get regions above a given threshold
        spline_values = array([int(x) for x in optimizedSpline(self.xRange)])
        if plotit is True:
            self.plot()

        starts_and_stops, starts, stops = self.get_regions_above_threshold(self.threshold,
                                                                      spline_values)
        peak_definitions = []
        for peak_start, peak_stop in starts_and_stops:
            peak_center_offset = np.where(self.find_local_maxima(spline_values[peak_start:(peak_stop + 1)]))[0]
            peak_center = [x + peak_start for x in peak_center_offset]
            if len(peak_center) == 1:
                peak_definitions.append((peak_start, peak_stop, peak_center[0]))
        self.peakCalls = peak_definitions
        return peak_definitions

    def plot(self, ax=None):
        self.peakCalls = self.peaks()
        if ax is None:
            ax = plt.gca()
        for peak in self.peakCalls:
            ax.axvline(x=peak[2], color='red', alpha=1, linewidth=4) #peak middle
            ax.axvspan(peak[0]+1, peak[1]-1, facecolor='blue', linewidth=2, alpha=.2)#peak span
        ax.axhline(y=self.threshold)
        ax.plot(self.yData, c='b', alpha=0.7)
        ax.plot(self.spline(self.xRange))


class Classic(PeakGenerator):
    """ Class to reimplement kaseys original peak calling method """   
    def __init__(self, xRange, yData, max_width, min_width, max_gap):
        
        """
        
        xRange -- the range to interpolate the spline over, must be monotonically increasing
        yData  -- the yAxis of the spline that corresponds to the xRange
                
        """
        
        super(Classic, self).__init__(xRange, yData)
        self.max_width = max_width
        self.min_width = min_width
        self.max_gap = max_gap
    
    #TODO:factor to init
    def peaks(self, plotit):
        peak_definitions = []
        
        in_peak = False
        peak_start = 0 
        for x, y in zip(self.xRange, self.yData):
            
            if in_peak is False and y != 0:
                in_peak = True
                peak_start = x
            
            #set peak stop
            if y != 0:
                peak_stop = x + 1
            
            #if the gap has been reached
            if in_peak and x - peak_stop >= self.max_gap:
                if peak_stop - peak_start < self.min_width:
                    #peak_stop = peak_start + self.min_width
                    pass
                #Change peak calculation and p-value min width calculation
                #also visualization min width should be ~10
                peak_center = peak_start + self.yData[peak_start:peak_stop].index(max(self.yData[peak_start:peak_stop]))
                peak_definitions.append((peak_start, peak_stop, peak_center))
                in_peak = False
            
            #if the max width has been reached
            if in_peak and peak_stop - peak_start >= self.max_width:
                peak_center = peak_start + self.yData[peak_start:peak_stop].index(max(self.yData[peak_start:peak_stop]))

                peak_definitions.append((peak_start, peak_stop, peak_center))
                in_peak = False
        
        #catch last case
        if in_peak: 
            if peak_stop - peak_start < self.min_width:
                peak_stop = peak_start + self.min_width
            
            peak_center = peak_start + self.yData[peak_start:peak_stop].index(max(self.yData[peak_start:peak_stop]))
            peak_definitions.append((peak_start, peak_stop, peak_center))
        
        return peak_definitions

from sklearn import mixture

class MyGMM(mixture.GMM):

    def bic(self, X, scoreWeight = 2, complexityWeight=4):
        return (-scoreWeight * self.score(X).sum() +
                complexityWeight * (self._n_parameters() * np.log(X.shape[0])))


class GaussMix(PeakGenerator):
    #Need to identify local minima for BIC, however fitting a gmm is expensive, we limit the total number of components
    #to search through, and if the bic stays at a minimium for max_checks we break out and call that out min
    def fit(self, max_components=50, max_checks=5):
        
        gmm_result = []
        for x, y in zip(self.xRange, self.yData):
            gmm_result += ([x] * y)

        #try multiple numbers of components, stop checking when you've been gradient-ascending too long
        best_model = None
        num_checks = 0
        for n_components in xrange(1, max_components):
            try:
                model = MyGMM(n_components, covariance_type='full').fit(gmm_result)
            except Exception as e:
                print e

            if best_model is None or model.bic(self.yData) < best_model.bic(self.yData):
                best_model = model
                self.n_components = n_components
                num_checks = 0
            else:
                num_checks += 1

            if num_checks >= max_checks:
                break


        self.bounds = 2 * np.sqrt(best_model.covars_).ravel()   # I think this is 2 sigma (stdev) (99% confidence?)
        self.model = best_model
        self.centers = best_model.means_.ravel()
        self.upper_bounds = self.centers + self.bounds
        self.lower_bounds = self.centers - self.bounds

    def peaks(self, plotit=False):
        
        if plotit:
            self.plot()

        peaks = []
        for prob, x, y, m in zip((self.model.score_samples(self.xRange)[1] > 0.8).T, self.lower_bounds, self.upper_bounds,
                                 self.centers):

            #clean up peaks, remove overlapping portions, keep peaks within xrange
            try:
                lowProbBnd = np.min(self.xRange[prob])
            except:
                lowProbBnd = 0
            try:
                highProbBnd = np.max(self.xRange[prob])
            except:
                highProbBnd = np.max(self.xRange)
            newX = np.max([0, lowProbBnd, x])
            
            newY = np.min([np.max(self.xRange), highProbBnd, y])
            #assert newX < m < newY
            newX = int(math.floor(newX))
            newY = int(math.floor(newY))
            m = int(np.round(m, 0))
            peaks.append((newX, newY, m))
            
        return peaks

    def plot(self, ax=None):

        if ax is None:
            ax = plt.gca()

        ax.plot(self.yData / float(sum(self.yData)))
        for mean, upper, lower, bound in zip(self.centers, self.upper_bounds, self.lower_bounds, self.bounds):
            #ax.axvline(mean, color='r')
            #ax.axvline(upper, color='g')
            #ax.axvline(lower, color='g')
            ax.plot(stats.norm.pdf(self.xRange, loc=mean, scale=bound))
            ax.set_xlim(0,)

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
        
        lam = 1 + ((float(reads_in_gene) / float(gene_length)) * float(peak_length)) #expect at least one read.

        cum_p = stats.poisson.sf(reads_in_peak, lam)

        return cum_p
    
    except Exception as error:
        logging.error("Poisson cutoff failled %s " % (error))
        return 1

#################################################################################################
### REMOVED IN V4
#################
#def read_array(reads, start, stop):
#    reads = (read for read in reads if (read.iv.start_d > start) & (read.iv.end_d < stop))
#    set_of_reads = HTSeq.GenomicArrayOfSets("auto", stranded=True)
#    for read in reads:
#        if read.aligned:
#                for cigop in read.cigar:
#                    if cigop.type != "M":
#                        continue
#                    set_of_reads[cigop.ref_iv] += read
#    return set_of_reads
#def get_reads_in_interval(interval, array_of_sets):
#    return set().union(*[baz for bar, baz in array_of_sets[interval].steps()])
#def count_reads_in_interval(interval, array_of_sets):
#    return len(get_reads_in_interval(interval, array_of_sets))
#def get_aligned_read_length(read):
#    return sum(cigar.size for cigar in read.cigar if cigar.type == "M")
#def read_lengths_from_htseq(reads):
#    return [get_aligned_read_length(read) for read in reads]
#################################################################################################
#
#
##################################################################################################
### ADDED IN V4
def get_reads_in_interval_pysam(interval, tx_start, full_read_array):
    start = max(0, interval.start - tx_start)
    stop = max(0, interval.stop - tx_start)
    return set().union(*[read for read in full_read_array[start:stop]])

def count_reads_in_interval_pysam(interval, tx_start, full_read_array):
    return len(get_reads_in_interval_pysam(interval, tx_start, full_read_array))


def get_aligned_read_length_pysam(read):
    return sum([code_len for code, code_len in read.cigartuples if code == 0])


def read_lengths_from_pysam(reads):
    return [get_aligned_read_length_pysam(read) for read in reads]
##################################################################################################


def bed_to_genomic_interval(bed):
    for interval in bed:
        yield HTSeq.GenomicInterval(str(interval.chrom), interval.start, interval.stop + 1, str(interval.strand))


def call_peaks(interval, gene_length, bam_file=None, max_gap=25,
               fdr_alpha=0.05, user_threshold=None, binom_alpha=0.05, method="binomial",
               min_reads=3, poisson_cutoff=0.05,
               plotit=False, w_cutoff=10, windowsize=1000, 
               SloP=False, max_width=None, min_width=None,
               algorithm="spline", reverse_strand=False, exons=None
               ###################################
               ### ADDED IN V4
               ,gene_no=0 ):                 # ADDED logging per gene feature (github call_peaks which does not log each gene)
               ###################################
    
    """

    calls peaks for an individual gene 
    
    interval - gtf interval describing the gene to query 
    takes bam file or bam file object.  Serial uses object parallel uses location (name)
    max_gap - space between sections for calling new peaks
    fdr_alpha - false discovery rate, p-value bonferoni correct from peaks script (called in setup)
    user_threshold - user defined FDR thershold (probably should be factored into fdr_alpha

    minreads - min reads in section to try and call peaks
    poisson_cutoff - p-value for signifance cut off for number of reads in peak that gets called - might want to use ashifted distribution
    plotit - makes figures 

    w_cutoff - width cutoff, peaks narrower than this are discarted 
    windowssize - for super local calculation distance left and right to look 
    SloP - super local p-value instead of gene-wide p-value
    max_width - int maximum with of classic peak calling algorithm peak
    min_width - int min width of classic peak calling algorithm peak
    max_gap   - int max gap of classic peak calling algorithm peak

    """
    ###########################################################################
    # print("starting call_peaks on gene_no:", gene_no, "interval:", interval)
    genecallpeaksloggingperiode = 100
    should_log_gene_call_peaks_this_time = (gene_no % genecallpeaksloggingperiode == 0)
    ###########################################################################
    if should_log_gene_call_peaks_this_time:
        logging.info(" starting call_peaks on gene_no {}".format(gene_no))
    ###########################################################################

    if plotit:
        plt.rcParams['interactive'] = True
        pass

    bam_fileobj = pysam.Samfile(bam_file, 'rb')
    #fixes non-standard chrom file names (without the chr)
    if not interval.chrom.startswith("chr"):
        interval.chrom = "chr" + interval.chrom

    subset_reads = list(bam_fileobj.fetch(reference=str(interval.chrom), start=interval.start, end=interval.stop))
    strand = str(interval.strand)
    if reverse_strand:
        if strand == "+":
            strand = "-"
        elif strand == "-":
            strand = "+"
#################################################################################################
### ADDED IN V4
#################
    (wiggle, jxns, pos_counts,
     lengths, allreads, read_locations) = readsToWiggle_pysam(subset_reads, interval.start,
                                              interval.stop, strand, "start", False)
#################################################################################################

#################################################################################################
### REMOVED IN V4
#################
#    (wiggle, jxns, pos_counts, lengths, allreads) = readsToWiggle_pysam(
#         subset_reads, interval.start, interval.stop, strand, "start", False)
#################################################################################################

#################################################################################################
### REMOVED IN V4
#################
#    #This is the worst of hacks, need to factor out pysam eventually
#    bam_fileobj = Robust_BAM_Reader(bam_file)
#    subset_reads = list(bam_fileobj.fetch(reference=str(interval.chrom), start=interval.start, end=interval.stop))
#    array_of_reads = read_array(subset_reads, interval.start, interval.stop)
#################################################################################################

    nreads_in_gene = sum(pos_counts)
    gene_length = int(gene_length)
    lengths = [gene_length - 1 if read >= gene_length else read for read in lengths]

    #pre-mRNA Threshold
    if user_threshold is None:
        if method == "binomial":  #Uses Binomial Distribution to get cutoff if specified by user
            #print len(lengths), gene_length, binom_alpha
            premRNA_threshold = get_FDR_cutoff_binom(lengths, gene_length, binom_alpha)
            #print premRNA_threshold
        elif method == "random":
            premRNA_threshold = get_FDR_cutoff_mean(readlengths=lengths,
                                                 genelength=gene_length,
                                                 alpha=fdr_alpha)
        else:
            raise ValueError("Method %s does not exist" % (method))
    else:
        logging.info("using user threshold")
        premRNA_threshold = user_threshold

    #mRNA Threshold
    exons = pybedtools.BedTool(exons)
    exons = exons.filter(lambda x: x.name == interval.attrs['gene_id']).saveas()

    total_exonic_reads = []
    total_exonic_length = 0
    htseq_exons = HTSeq.GenomicArrayOfSets(chroms="auto", stranded=False)

    for exon, exon_interval in zip(exons, bed_to_genomic_interval(exons)):
        exon.stop += 1
        exonic_reads = get_reads_in_interval_pysam(exon, interval.start, read_locations)

        exon_read_lengths = read_lengths_from_pysam(exonic_reads)
        exon_read_lengths = [exon_interval.length - 1 if read > exon_interval.length else read for read in exon_read_lengths]
        total_exonic_reads += exon_read_lengths
        total_exonic_length += exon_interval.length
        htseq_exons[exon_interval] += 'exon'

    mRNA_threshold = get_FDR_cutoff_binom(total_exonic_reads, total_exonic_length, binom_alpha)
    if not isinstance(premRNA_threshold, int):
        raise TypeError

    #these are what is built in this dict, complicated enough that it might
    #be worth turning into an object
    peak_dict = {}
    peak_dict['clusters'] = []
    peak_dict['sections'] = {}
    peak_dict['nreads'] = int(nreads_in_gene)
    peak_dict['threshold'] = premRNA_threshold
    peak_dict['loc'] = interval

    peak_number = 0

    sections = find_sections(wiggle, max_gap)
    if plotit:
        plot_sections(wiggle, sections, premRNA_threshold)

    for sect in sections:

        sectstart, sectstop = sect
        sect_length = sectstop - sectstart + 1
        data = wiggle[sectstart:(sectstop + 1)]

        cur_interval = HTSeq.GenomicInterval(str(interval.chrom),
                                             sectstart + interval.start,
                                             sectstop + interval.start + 1,
                                             strand)

        #Logic to use variable thresholds for exons or introns, still superseded by superLocal logic
        overlaps_exon = len(reduce( set.union, ( val for iv, val in htseq_exons[cur_interval].steps()))) > 0
        gene_threshold = mRNA_threshold if overlaps_exon else premRNA_threshold

        #maybe make a function that takes a GI and converts it into a pybedtools interval
        cur_pybedtools_interval = pybedtools.create_interval_from_list([interval.chrom,
                                              sectstart + interval.start,
                                              sectstop + interval.start + 1,
                                              interval.name,
                                              interval.score,
                                              strand])

        Nreads = count_reads_in_interval_pysam(cur_pybedtools_interval, interval.start, read_locations)

        cts = pos_counts[sectstart:(sectstop + 1)]
        xvals = arange(len(data))
        peak_dict['sections'][sect] = {}
        peak_dict['sections'][sect]['nreads'] = int(Nreads)

        #makes sure there are enough reads
        if Nreads < min_reads:
            logging.info("""%d is not enough reads, skipping section: %s""" % (Nreads, sect))
            peak_dict['sections'][sect]['tried'] = False
            continue
        else:
            logging.info("""Analyzing section %s with %d reads""" % (sect, Nreads))
            pass

        if user_threshold is None:
            if SloP:
                half_width = 500
                section_start = max(0, sectstart + interval.start - half_width)
                section_stop = sectstop + interval.start + 1 + half_width
                expanded_sect_length = section_stop - section_start

                cur_pybedtools_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                                section_start,
                                                                                section_stop,
                                                                                interval.name,
                                                                                interval.score,
                                                                                strand])

                expanded_Nreads = get_reads_in_interval_pysam(cur_pybedtools_interval, interval.start, read_locations)
                sect_read_lengths = read_lengths_from_pysam(expanded_Nreads)
                sect_read_lengths = [sect_length - 1 if read > sect_length else read for read in sect_read_lengths]
                peak_dict['sections'][sect]['expanded_Nreads'] = len(expanded_Nreads)

                if method == "binomial":  #Uses Binomial Distribution to get cutoff if specified by user
                    slop_threshold = get_FDR_cutoff_binom(readlengths=sect_read_lengths,
                                                          genelength=expanded_sect_length,
                                                          alpha=binom_alpha)
                elif method == "random":
                    #use the minimum FDR cutoff between superlocal and gene-wide calculations
                    slop_threshold = get_FDR_cutoff_mean(readlengths=sect_read_lengths,
                                                         genelength=expanded_sect_length,
                                                         alpha=fdr_alpha)
                else:
                    raise ValueError("Method %s does not exist" % (method))
                threshold = max(gene_threshold, slop_threshold)

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

        if max(data) < threshold:
            logging.info("data does not excede threshold, stopping")
            continue

        if algorithm == "spline":
            data = map(float, data)
            #Magic number for initial smoothing, but it works
            initial_smoothing_value = ((sectstop - sectstart + 1)**(1/3)) + 10

            peak_dict['sections'][sect]['smoothing_factor'] = initial_smoothing_value

            logging.info("initial smoothing value: %.2f" % initial_smoothing_value)
            fitter = SmoothingSpline(xvals, data, smoothing_factor=initial_smoothing_value,
                            lossFunction="get_turn_penalized_residuals",
                            threshold=threshold,
                            num_reads=Nreads)

        elif algorithm == "gaussian":
            cts = map(float, cts)
            fitter = GaussMix(xvals, cts)

        elif algorithm == "classic":
            data = map(float, data)
            fitter = Classic(xvals, data, max_width, min_width, max_gap)

        try:
            peak_definitions = fitter.peaks()
            logging.info("optimized smoothing value: %.2f" % fitter.smoothing_factor)
            peak_dict['sections'][sect]['final_smoothing_factor'] = fitter.smoothing_factor
            if peak_definitions is None:
                numpeaks = 0
            else:
                numpeaks = len(peak_definitions)
            logging.info("I identified %d potential peaks" % (numpeaks))

        except Exception as error:
            logging.error("peak finding failed:, %s, %s" % (interval.name, error))
            raise error

        #subsections that are above threshold
        #peak center is actually the location where we think binding should
        #occur, not the average of start and stop

        #Need to get all ranges, count number of reads in each range and compute from there
        for peak_start, peak_stop, peak_center in peak_definitions:

            genomic_start = interval.start + sectstart + peak_start
            genomic_stop = interval.start + sectstart + peak_stop

            cur_pybedtools_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                                genomic_start,
                                                                                genomic_stop,
                                                                                interval.name,
                                                                                interval.score,
                                                                                strand])

            number_reads_in_peak = count_reads_in_interval_pysam(cur_pybedtools_interval, interval.start, read_locations)

            peak_length = genomic_stop - genomic_start + 1

            logging.info("""Peak %d (%d - %d) has %d
                          reads""" % (peak_number, peak_start,
                                     (peak_stop + 1), number_reads_in_peak))

            #highest point in start stop
            genomic_center = interval.start + sectstart + peak_center

            #makes it thicker so we can see on the browser
            #error checking logic to keep bed files from breaking
            thick_start = max(genomic_center - 2, genomic_start)
            thick_stop = min(genomic_center + 2, genomic_stop)


            #super local logic
            area_start = max(0, (peak_center + sectstart) - windowsize)
            area_stop = min((peak_center + sectstart) + windowsize, len(wiggle))

            cur_pybedtools_interval = pybedtools.create_interval_from_list([interval.chrom,
                                                                                interval.start + area_start,
                                                                                interval.start + area_stop,
                                                                                interval.name,
                                                                                interval.score,
                                                                                strand])

            number_reads_in_area = count_reads_in_interval_pysam(cur_pybedtools_interval, interval.start, read_locations)
            area_length = area_stop - area_start + 1

            peak_dict['clusters'].append(Peak(chrom=interval.chrom,
                                              genomic_start=genomic_start,
                                              genomic_stop=genomic_stop,
                                              gene_name=interval.attrs['gene_id'],
                                              strand=interval.strand,
                                              thick_start=thick_start,
                                              thick_stop=thick_stop,
                                              peak_number=peak_number,
                                              number_reads_in_peak=number_reads_in_peak,
                                              size=peak_length,
                                              p=0,
                                              effective_length=int(interval.attrs['effective_length']),
                                              peak_length=peak_length,
                                              area_reads=number_reads_in_area,
                                              area_size=area_length,
                                              nreads_in_gene=nreads_in_gene,
                                              #nreads_in_input=input_number_reads_in_peak,
                                              ))

            peak_number += 1
            peak_dict['sections'][sect]['nPeaks'] += 1

    peak_dict['Nclusters'] = peak_number
    if plotit:
        import sys
        plt.show()
        v = sys.stdin.read(1)
    ###################################################
    # print("returning gene_no:", gene_no, "peak_dict:", peak_dict)
    ####################################################

    return peak_dict
