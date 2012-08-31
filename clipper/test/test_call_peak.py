'''
Created on Jul 25, 2012

@author: gabrielp
'''
import unittest
from clipper.src.call_peak import *
from clipper.src.peaks import *
from numpy import *
from numpy.testing import *
from scipy import interpolate

            
class Test(unittest.TestCase):



    def test_get_FDR_cutoff_mode(self):
        
        """
    
        Performs unit tests on get_FDR_cutoff_mode function
        Function is currently deperacated, tests not done
        
        """
        
        pass
    

    def test_get_FDR_cutoff_mean(self):
        
        """
    
        Performs unit tests on get_FDR_cutoff_mean function
        
        Difficult to test because of random sampling
        
        TODO: Figure out how to manually calculate FDR cutoff, right now 
        I just took the old result and am using that.  Math appears to be correct,
        but this is a really bad practice
        
        """
        
        #Case: Not enough reads, expected result: returns the passed min cutoff 
        result = get_FDR_cutoff_mean([], 100)
        assert result == 2
        
        #Case: Not enough reads with different min cutoff value
        
        result = get_FDR_cutoff_mean([], 100, mincut = 10)
        assert result == 10
        
        #Case, enough reads, mean is higher than minimum cutoff, expected result: mean is returned
        #setup, create read lengths 
        
        #mean can also be estimated by emperically (N * L) / G
        #N is number of reads (20)
        #L is number of nucleotides per read (30)
        #G is genome size (in this case an RNA (100)
        #Read depth should be 6
        
        read_lengths = [30] * 100
        result = get_FDR_cutoff_mean(read_lengths, 1000)
        self.assertEqual(result, 7)
        
        #Second similar case 
        read_lengths = [30] * 20
        result = get_FDR_cutoff_mean(read_lengths, 100)
        assert result == 11 
       
        
        #Case: enough reads, mean is lower than minimum cutoff, expected result: minimum cutoff is returned
        result = get_FDR_cutoff_mean(read_lengths, 100, mincut = 20)
        assert result == 20
        

    def test_plotSpline(self):
        
        """
    
        Quality control function, not important to main running of program
        Not tested
        
        """
        
        pass
        #plotSections([5] * 10, ["1|5", "7|9"], 3)
        #assert 1 == 0
    

    def test_plotSections(self):
        
        """
    
        tests plotSections function, appears to have computer specific issues
        

        Plots each section individually, I think
        Wiggle is a list representing a wiggle track
        sections is a list of strings of format "start|stop" where start and stop are both integers
        threshold is an integer 

        """     
        
        pass
            

    def test_find_univariateSpline(self):
        
        """
    
        Performs unit testing on find_univariate_spline
        As this is mostly a wrapper for a scipy function I will not test spline calling
        beyond basics.
        
        These tests will verify that the resid logic works, and eventually these 
        will be factored into two functions.  Its bad form to have to toggle your return
        value in the function
        
        """
        
        #Case null inputs, expected:  Everything goes to hell, not testing
        
        #Case resid is false, expected: returns univariateSpline 
        #spline, just verifies that it is the same as 
        #the scipy result
        
        #setup 
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2
        smoothing = 5 
        #expected
        expected = interpolate.UnivariateSpline(xvals, data, k=3, s=smoothing)
        
        #test
        result = find_univariate_spline(smoothing, xvals, data, 3)
        
        #hacky test, but they should be about the same
        self.assertAlmostEqual(expected.get_residual(), result.get_residual()) 
        
        #tests error mode
        self.assertRaises(TypeError, find_univariate_spline, None, None, None, None)
  
    def test_find_spline_residuals(self):
        
        """
        
        Test for find_spline_residuals function
        
        """
        
        #setup 
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2
        smoothing = 5 
        #expected
        expected = interpolate.UnivariateSpline(xvals, data, k=3, s=smoothing)
        residual = find_spline_residuals(smoothing, xvals, data, 3)
        self.assertAlmostEqual(residual, sqrt(expected.get_residual()))
        
        #tests error mode
        self.assertRaises(TypeError, find_spline_residuals, None, None, None, None)
        
        

    def test_poissonP(self):
        
        """
    
        Performs unit testing on poissonP
        
        Will not test math behind calling poisson fuction more than nessessary as it has already been tested by scipy
        
        """ 
        
        #Case: fewer than 3 reads fall within a peak region: Expected result should return as though there are 3 expected reads
        result = poissonP(50, 3, 50, 2)
        self.assertAlmostEqual(result, (1 - 0.64723188878223115)) #manually calculated on personal computer stats.poisson.cdf(3, 3)
        
        #Case: More than 3 reads fall within a peak region. Expected: Result should return as normal
        result = poissonP(50, 3, 50, 4)
        self.assertAlmostEqual(result, (1 - 0.4334701203667089)) #manually calculated stats.poisson.cdf(3, 4)
    
        #Case: poissonp throws an exception. Expected: Result should return 1
        result = poissonP(None, None, None, None)    
        assert 1 == result
    
    
    def test_call_peaks(self):
        
        """
    
        Performs unit testing on call_peaks
        
        will not test peaks_from_info here, just the error handling of call peaks
        Need to create a dummy call peaks function or stub or something to keep everything from
        getting called
        
        The way this is currently written it is difficult to test logic other than error checking
        
        """ 

        #peaks_from_info = poissonP
        """call_peaks("chr1|foo|5|10|-", 
                   50, 
                   bam_fileobj=None, 
                   bam_file=None)"""
        assert 1 == 1
        #Case: Bam file is passed
        
        #Case: Bam object is passed
        
        #Case: Both bam file and object are passed
    

    def test_peaks_from_info(self):
        
        """
    
        Performs unit testing on peaks_from_info
        
        Testing of this is currently handled in the allup test
        This method is currently to complicated to test as is
        
        """
        
        pass


    def test_get_start_stop_pairs_above_threshold(self):
        
        """    
    
        tests generating start and stops function,
        need better tests and edge cases 
        
        """
        
        #need to add null and empty inputs 
        
        #Test general flow
        values = array([1,2,3,4,5,5,5,3,2,1])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,8)])
        assert_array_equal(starts, [2])
        assert_array_equal(stops, [8])
        
        #Test starting above threshold
        values = array([4,4,5,5,5,3,2,1])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(0,6)])
        assert_array_equal(starts, [0])
        assert_array_equal(stops, [6])

        #Test ending above threshold
        values = array([1,2,3,4,5,5,5,5,5,5])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,10)])
        assert_array_equal(starts, [2])
        assert_array_equal(stops, [10])
        
        #Test local minima 
        values = array([1,2,3,4,5,5,4,5,5,3])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (6,10)])
        assert_array_equal(starts, [2, 6])
        assert_array_equal(stops, [6,10])
        
        #Test Two peaks
        values = array([1,2,3,4,5,5,2,5,5,3])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (7,10)])
        assert_array_equal(starts, [2,7])
        assert_array_equal(stops, [6,10])
        
        values = array([1,2,3,4,5,5,2,2,5,4,4,3])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (8,12)])
        assert_array_equal(starts, [2,8])
        assert_array_equal(stops, [6,12])
        
        #test two peaks starting above 
        values = array([5,5,5,5,5,5,2,2,5,4,4,3])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(0,6), (8,12)])
        assert_array_equal(starts, [0,8])
        assert_array_equal(stops, [6,12])
        
        #test two peaks ending above
        values = array([1,2,3,4,5,5,2,2,5,4,4,4])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (8,12)])
        assert_array_equal(starts, [2,8])
        assert_array_equal(stops, [6,12])
        
        #test two peaks with one local minima
        values = array([5,5,5,4,5,5,2,2,5,4,4,3])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(0,3), (3,6), (8,12)])
        assert_array_equal(starts, [0, 3, 8])
        assert_array_equal(stops, [3, 6,12])
        
        #test two peaks with two local minima
        values = array([1,2,3,4,5,5,2,2,5,4,5,5])
        starts_and_stops, starts, stops = get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (8,9), (9, 12)])
        assert_array_equal(starts, [2,8,9])
        assert_array_equal(stops, [6,9,12])
        
        #more complicated version
        values = array([3,2,1,2,3,4,3,2,1,2,3,4,5,4,3,2,3,4,5,3,2,0,3])
        threshold = 2
        starts_and_stops, starts, stops = get_regions_above_threshold(threshold, values)
        
        #Not sure if I want to have the last value be a possible peak, 
        #I guess filtering happens later so it doesn't matter
        assert_array_equal(starts_and_stops, [(0,2), (3,8), (9, 15), (15,21), (22, 23)])
        assert_array_equal(starts, [0, 3, 9, 15, 22])
        assert_array_equal(stops, [2, 8, 15, 21, 23])
        
        #more complicated version
        values = array([3,2,1,2,3,4,3,2,1,2,3,4,5,4,3,2,3,4,5,3,2,0,2,3])
        threshold = 2
        starts_and_stops, starts, stops = get_regions_above_threshold(threshold, values)
        assert_array_equal(starts_and_stops, [(0,2), (3,8), (9, 15), (15,21), (22, 24)])

        
        #more complicated version
        values = array([0,0,0, 3, 0, 0 ,0])
        threshold = 2
        starts_and_stops, starts, stops = get_regions_above_threshold(threshold, values)
        assert_array_equal(starts_and_stops, [(3,4)])

        
        #Local minima that has a range
        values = array([10,10,10,10,10,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10])
        threshold = 2
        starts_and_stops, starts, stops = get_regions_above_threshold(threshold, values)
        assert_array_equal(starts_and_stops, [(0,8), (8,20)])

        
        #failing on real data
        values = array([   6.00000000e+00,   6.00000000e+00,   6.00000000e+00,   6.00000000e+00,
                           6.00000000e+00,   6.00000000e+00,   8.00000000e+00,   8.00000000e+00,
                           8.00000000e+00,   1.50000000e+01,   1.50000000e+01,   3.30000000e+01,
                           3.30000000e+01,   3.30000000e+01,   3.30000000e+01,   3.30000000e+01,
                           3.70000000e+01,   3.70000000e+01,   3.70000000e+01,   3.70000000e+01,
                           3.70000000e+01,   3.70000000e+01,   3.70000000e+01,   3.70000000e+01,
                           3.70000000e+01,   3.70000000e+01,   3.70000000e+01,   3.70000000e+01,
                           3.70000000e+01,   3.70000000e+01,   4.20000000e+01,   4.20000000e+01,
                           4.20000000e+01,   3.60000000e+01,   3.60000000e+01,   3.60000000e+01,
                           3.60000000e+01,   3.60000000e+01,   3.60000000e+01,   3.40000000e+01,
                           4.60000000e+01,   5.40000000e+01,   5.00000000e+01,   5.00000000e+01,
                           3.20000000e+01,   3.20000000e+01,   3.20000000e+01,   3.20000000e+01,
                           3.20000000e+01,   2.80000000e+01,   2.80000000e+01,   2.80000000e+01,
                           2.80000000e+01,   2.80000000e+01,   4.50000000e+01,   4.50000000e+01,
                           5.70000000e+01,   8.00000000e+01,   8.00000000e+01,   1.79000000e+02,])
        
        starts_and_stops, starts, stops = get_regions_above_threshold(32, values)
        
        assert_array_equal(starts_and_stops, [(11, 39), (39, 49), (54, 60)])
        
        #more real tests

    
    def test_find_local_minima(self):
        
        """
        
        tests find local minima range function
        
        """
        
        #inital tests used in base version of find starts and stops
        
        #Test local minima 
        values = array([1,2,3,4,5,5,4,5,5,3])
        result = find_local_minima(values)
        true = array([False,False,False,False,False,False,True,False,False,False])
        assert_array_equal(true, result)
        
        #test two peaks with one local minima
        values = array([5,5,5,4,5,5,2,2,5,4,4,3])
        result = find_local_minima(values)
        true = array([False,False,False,True,False,False,True,False,False,False, False, False])
        assert_array_equal(true, result)

        
        #test two peaks with two local minima
        values = array([1,2,3,4,5,5,2,2,5,4,5,5])
        result = find_local_minima(values)
        true = array([False,False,False,False,False,False,True, False, False,True, False, False])
        assert_array_equal(true, result)
        
        #Test long array
        values = array([10,10,9,9,9,9,9,9,9,9,10,10])
        result = find_local_minima(values)
        true = array([False,False,False,False,False,True,False, False, False,False, False, False])
        assert_array_equal(true, result)
        
        #more failing stuff
        values = array([32,  32,  32,  28,  28,  28,  28,  28,  45,  45,  57,  80,  80])
        result = find_local_minima(values)
        true = array([False,False,False,False,False,True,False, False, False,False, False, False, False])
        assert_array_equal(true, result)
        
    def test_find_local_maxima(self):
        
        """
        
        Tests find local maxima function
        
        Need to think of better tests for this, but I'm lazy right now...
        
        """
        
        #Test local minima 
        values = array([1,2,3,4,5,5,4,5,5,3])
        result = find_local_maxima(values)
        true = array([False,False,False,False,True,False,False,True,False,False])
        assert_array_equal(true, result)
        
        #test two peaks with one local minima
        values = array([5,5,5,4,5,5,2,2,5,4,4,3])
        result = find_local_maxima(values)
        true = array([False,True,False,False,True,False, False ,False,True, False, False, False])
        assert_array_equal(true, result)

        #test two peaks with two local minima
        values = array([1,2,3,4,5,5,2,2,5,4,5,5])
        result = find_local_maxima(values)
        true = array([False,False,False,False,True,False,False, False, True, False, True, False])
        assert_array_equal(true, result)
        
        #Test long array
        values = array([10,10,9,9,9,9,9,9,9,9,10,10])
        result = find_local_maxima(values)
        true = array([True,False,False,False,False,False,False, False, False,False, True, False])
        assert_array_equal(true, result)
        
        #test something that is breaking in reality 
        values = array([ 32,  84,  85,  85,])
        result = find_local_maxima(values)
        true = array([False, False, True, False])
        assert_array_equal(true, result)
        
        #more failing stuff
        values = array([32,  32,  32,  28,  28,  28,  28,  28,  45,  45,  57,  80,  80])
        result = find_local_maxima(values)
        true = array([False,True,False,False,False,False,False, False, False, False, False, True, False])
        assert_array_equal(true, result)
        
    def test_get_maxima_real(self):
        #more real tests
        values = array([ 201,  201,  201,  201,  201,  201,  201,  201,  201,  201,  201,  201,
          201,  201,  201,  201,  205,  196,   82,   -1])
        result = find_local_maxima(values)
        true = array([False,False,False,False,False,False,False, False, False, False, False, False, False, False, False, False, True, False, False, False])
        assert_array_equal(true, result)
    def test_get_start_stop_pairs_above_threshold_regression(self):
        
        """
        
        regresstion test to make sure I didn't break anything
        tests generating start and stops function, this is kind of a badly written function,
        need better tests and edge cases 
        
        """
        
        #test to make sure stuff doesn't break, should have made regression test...
        cutoff = 71.5879289615 
        xvals = array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                  23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
                  46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68])
        data = array([2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]) 
        degree = 3 
        weights = None
        threshold = 3

        spline = find_univariate_spline(cutoff, xvals, data, degree, weights)
        
        starts_and_stops, starts, stops = get_regions_above_threshold(threshold, spline(xvals))
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
