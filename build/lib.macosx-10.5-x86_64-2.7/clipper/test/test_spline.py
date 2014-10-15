'''
Created on Jan 8, 2013

@author: gabrielp

Spline is a confusing name
TODO: RENAME class 

tests Class to fit data to a smooth curve

will use inheretince scheme, either a gausian pattern or the old style

will need to generalize tests to both subclasses, main class will probably be a 
interface

'''
import unittest
from numpy import *
from numpy.testing import *
from scipy import interpolate
from clipper.src.call_peak import * 
class Test(unittest.TestCase):


    def setUp(self):
        pass


    def tearDown(self):
        pass


    def testName(self):
        pass

    def test_init(self):
        
        """
        
        Tests init
        
        """
        
        #basic case, should work
        smoothing_spline = SmoothingSpline([], [])
        self.assertEqual(smoothing_spline.lossFunction , smoothing_spline.get_turn_penalized_residuals)
        
        smoothing_spline = SmoothingSpline([], [], lossFunction = "get_turn_penalized_residuals")
        self.assertEqual(smoothing_spline.lossFunction , smoothing_spline.get_turn_penalized_residuals)
        
        smoothing_spline = SmoothingSpline([], [], lossFunction = "get_norm_penalized_residuals")
        self.assertEqual(smoothing_spline.lossFunction , smoothing_spline.get_norm_penalized_residuals)
        
        self.assertRaises(TypeError, SmoothingSpline.__init__, [], [], lossFunction = "foo")
        
        smoothing_spline = SmoothingSpline([5,5,5,5], [1,2,3,4])
        self.assertEqual(smoothing_spline.smoothingFactor , 4) 
        
        smoothing_spline = SmoothingSpline([5,5,5,5], [1,2,3,4], smoothingFactor = 6 )
        self.assertEqual(smoothing_spline.smoothingFactor , 6)        
        
    def test_fit_univariate_spline(self):
        
        """
    
        Performs unit testing on find_univariate_spline
        As this is mostly a wrapper for a scipy function I will not test spline calling
        beyond basics.
        
        The return of the result as well as storing it as a class variable is bad form
        I will try to figure out how to get rid of it
        
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
        smoothing_spline = SmoothingSpline(xvals, data)
        
        result = smoothing_spline.fit_loss(smoothing, replace, weight)
        
        #hacky test, but they should be about the same
        self.assertAlmostEqual(expected.get_residual(), result.get_residual()) 
        
        #tests error mode
        #self.assertRaises(TypeError, find_univariate_spline, None, None, None, None)
     
    def test_fit_loss(self):
        
        """
        
        Tests fit loss
        
        """
        
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2
        smoothing = 5 
        #expected
        expected = interpolate.UnivariateSpline(xvals, data, k=3)
        
        smoothing_spline = SmoothingSpline(xvals, data)
        self.assertEqual(smoothing_spline.get_turn_penalized_residuals(expected), smoothing_spline.fit_loss())
        
        smoothing_spline = SmoothingSpline(xvals, data, lossFunction="get_norm_penalized_residuals")
        self.assertEqual(smoothing_spline.get_norm_penalized_residuals(expected), smoothing_spline.fit_loss())

    
    def test_optimize_fit(self):
        
        """
        
        Tests optimize fit
        
        """
        
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2
        smoothing = 5 
        #expected
        expected = interpolate.UnivariateSpline(xvals, data, k=3)
        
        smoothing_spline = SmoothingSpline(xvals, data)
        
        result = smoothing_spline.optimize_fit()
        result.get_residual()
        
        xvals = [ 0,  1,  2, 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
        ydata = [9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 9930, 0]
        
        smoothing_spline = SmoothingSpline(xvals, ydata, lossFunction="get_norm_penalized_residuals")
        smoothing_spline.optimize_fit(225)
        
    def test_peaks(self):
        
        """
        
        Tests peaks
        
        regression test
        
        need to test: 
        
        really the only thing I need to test is the spline values 
        starts and stops - 
        (spline_values, starts_and_stops, starts, stops)
        
        Not testing main function...
        """

        xvals = arange(0, 162)
        data = [16, 16, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 15, 17, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 17, 32, 32, 84, 85, 85, 85, 85, 85, 85, 85, 85, 85, 95, 109, 111, 112, 112, 112, 112, 112, 112, 117, 117, 117, 117, 117, 117, 117, 116, 116, 116, 116, 100, 85, 85, 33, 32, 32, 32, 32, 32, 32, 32, 32, 34, 24, 26, 36, 35, 35, 35, 35, 35, 35, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 28, 28, 12, 0] 
        initial_smoothing_value = 162

        fitter = SmoothingSpline(xvals, data, initial_smoothing_value,
                            lossFunction="get_norm_penalized_residuals",
                            threshold=32)
        
        starts_and_stops = fitter.peaks( False)
        print starts_and_stops
        self.assertListEqual([(0, 8, 0), (67, 148, 114)], starts_and_stops) 
         
    def test_get_norm_penalized_residuals(self):
        
        """
        
        Tests get norm penalized residuals
        
        Should shore this up, and do error checking...
        """
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2
        smoothing = 5 
        #expected
        expected = interpolate.UnivariateSpline(xvals, data, k=3, s=smoothing)
        #(10*norm(expected(xvals))**5) + 1 * sqrt(expected.get_residual())
        #test
        smoothing_spline = SmoothingSpline(xvals, data)
        
        result = smoothing_spline.get_norm_penalized_residuals(expected)
        
        
        #here the final math result
        self.assertAlmostEqual(74590259.948699772, result, 4)
    def test_get_turn_penalized_residuals(self):
        
        """
        
        Tests get turn penalized residuals
        
        This isn't a great test, because I don't know how it should work
        its an arbitray error function though, so it should be alright
        """
        
        
        
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2
        smoothing = 5 
        #expected
        expected = interpolate.UnivariateSpline(xvals, data, k=3, s=smoothing)
        #(10*norm(expected(xvals))**5) + 1 * sqrt(expected.get_residual())
        #test
        smoothing_spline = SmoothingSpline(xvals, data)
        
        turns = sum(abs(diff(sign(diff(expected(xvals)))))) / 2
        self.assertEqual(1, turns, "turn calculation is wrong")
        
        self.assertAlmostEqual(2.235811645548004, smoothing_spline.get_turn_penalized_residuals(expected), 4) 
   
        
    def test_fit_univariate_spline(self):
        
        """
        
        Tests fit univariate spline
        
        """
        
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2
        smoothing = 5 
        #expected
        expected = interpolate.UnivariateSpline(xvals, data, k=3, s=smoothing)
        #(10*norm(expected(xvals))**5) + 1 * sqrt(expected.get_residual())
        #test
        smoothing_spline = SmoothingSpline(xvals, data)
        
        #tests automatic setting of spline smoothing factor
        expected = interpolate.UnivariateSpline(xvals, data, k=3, s=len(xvals))
        result = smoothing_spline.fit_univariate_spline()
        self.assertEqual(expected.get_residual(), result.get_residual())
        
        #tests manual setting of factor
        expected = interpolate.UnivariateSpline(xvals, data, k=3, s=smoothing)
        result = smoothing_spline.fit_univariate_spline(smoothingFactor=smoothing)
        self.assertEqual(expected.get_residual(), result.get_residual())
        
    def test_find_local_minima(self):
        
        """
        
        tests find local minima range function
        
        """
        
        #inital tests used in base version of find starts and stops
        smoothing_spline = SmoothingSpline([], [])
        #Test local minima 
        values = array([1,2,3,4,5,5,4,5,5,3])
        result = smoothing_spline.find_local_minima(values)
        true = array([False,False,False,False,False,False,True,False,False,False])
        assert_array_equal(true, result)
        
        #test two peaks with one local minima
        values = array([5,5,5,4,5,5,2,2,5,4,4,3])
        result = smoothing_spline.find_local_minima(values)
        true = array([False,False,False,True,False,False,True,False,False,False, False, False])
        assert_array_equal(true, result)

        
        #test two peaks with two local minima
        values = array([1,2,3,4,5,5,2,2,5,4,5,5])
        result = smoothing_spline.find_local_minima(values)
        true = array([False,False,False,False,False,False,True, False, False,True, False, False])
        assert_array_equal(true, result)
        
        #Test long array
        values = array([10,10,9,9,9,9,9,9,9,9,10,10])
        result = smoothing_spline.find_local_minima(values)
        true = array([False,False,False,False,False,True,False, False, False,False, False, False])
        assert_array_equal(true, result)
        
        #more failing stuff
        values = array([32,  32,  32,  28,  28,  28,  28,  28,  45,  45,  57,  80,  80])
        result = smoothing_spline.find_local_minima(values)
        true = array([False,False,False,False,False,True,False, False, False,False, False, False, False])
        assert_array_equal(true, result)
        
    def test_get_start_stop_pairs_above_threshold(self):
        
        """    
    
        tests generating start and stops function,
        need better tests and edge cases 
        
        """
        
        #need to add null and empty inputs 
        smoothing_spline = SmoothingSpline([], [])
        #Test general flow
        values = array([1,2,3,4,5,5,5,3,2,1])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,8)])
        assert_array_equal(starts, [2])
        assert_array_equal(stops, [8])
        
        #Test starting above threshold
        values = array([4,4,5,5,5,3,2,1])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(0,6)])
        assert_array_equal(starts, [0])
        assert_array_equal(stops, [6])

        #Test ending above threshold
        values = array([1,2,3,4,5,5,5,5,5,5])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,10)])
        assert_array_equal(starts, [2])
        assert_array_equal(stops, [10])
        
        #Test local minima 
        values = array([1,2,3,4,5,5,4,5,5,3])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (6,10)])
        assert_array_equal(starts, [2, 6])
        assert_array_equal(stops, [6,10])
        
        #Test Two peaks
        values = array([1,2,3,4,5,5,2,5,5,3])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (7,10)])
        assert_array_equal(starts, [2,7])
        assert_array_equal(stops, [6,10])
        
        values = array([1,2,3,4,5,5,2,2,5,4,4,3])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (8,12)])
        assert_array_equal(starts, [2,8])
        assert_array_equal(stops, [6,12])
        
        #test two peaks starting above 
        values = array([5,5,5,5,5,5,2,2,5,4,4,3])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(0,6), (8,12)])
        assert_array_equal(starts, [0,8])
        assert_array_equal(stops, [6,12])
        
        #test two peaks ending above
        values = array([1,2,3,4,5,5,2,2,5,4,4,4])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (8,12)])
        assert_array_equal(starts, [2,8])
        assert_array_equal(stops, [6,12])
        
        #test two peaks with one local minima
        values = array([5,5,5,4,5,5,2,2,5,4,4,3])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(0,3), (3,6), (8,12)])
        assert_array_equal(starts, [0, 3, 8])
        assert_array_equal(stops, [3, 6,12])
        
        #test two peaks with two local minima
        values = array([1,2,3,4,5,5,2,2,5,4,5,5])
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(3, values)
        assert_array_equal(starts_and_stops, [(2,6), (8,9), (9, 12)])
        assert_array_equal(starts, [2,8,9])
        assert_array_equal(stops, [6,9,12])
        
        #more complicated version
        values = array([3,2,1,2,3,4,3,2,1,2,3,4,5,4,3,2,3,4,5,3,2,0,3])
        threshold = 2
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(threshold, values)
        
        #Not sure if I want to have the last value be a possible peak, 
        #I guess filtering happens later so it doesn't matter
        assert_array_equal(starts_and_stops, [(0,2), (3,8), (9, 15), (15,21), (22, 23)])
        assert_array_equal(starts, [0, 3, 9, 15, 22])
        assert_array_equal(stops, [2, 8, 15, 21, 23])
        
        #more complicated version
        values = array([3,2,1,2,3,4,3,2,1,2,3,4,5,4,3,2,3,4,5,3,2,0,2,3])
        threshold = 2
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(threshold, values)
        assert_array_equal(starts_and_stops, [(0,2), (3,8), (9, 15), (15,21), (22, 24)])

        
        #more complicated version
        values = array([0,0,0, 3, 0, 0 ,0])
        threshold = 2
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(threshold, values)
        assert_array_equal(starts_and_stops, [(3,4)])

        
        #Local minima that has a range
        values = array([10,10,10,10,10,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10])
        threshold = 2
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(threshold, values)
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
        
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(32, values)
        
        assert_array_equal(starts_and_stops, [(11, 39), (39, 49), (54, 60)])
        
        #more real tests
    
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

        threshold = 3
        
        smoothing_spline = SmoothingSpline(xvals, data)
        spline = smoothing_spline.fit_univariate_spline(smoothingFactor = cutoff)
        
        starts_and_stops, starts, stops = smoothing_spline.get_regions_above_threshold(threshold, spline(xvals))

    

        
    def test_find_local_maxima(self):
        
        """
        
        Tests find local maxima function
        
        Need to think of better tests for this, but I'm lazy right now...
        
        """
        
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2

        #test
        smoothing_spline = SmoothingSpline(xvals, data)

        #Test local minima 
        values = array([1,2,3,4,5,5,4,5,5,3])
        result = smoothing_spline.find_local_maxima(values)
        true = array([False,False,False,False,True,False,False,True,False,False])
        assert_array_equal(true, result)
        
        #test two peaks with one local minima
        values = array([5,5,5,4,5,5,2,2,5,4,4,3])
        result = smoothing_spline.find_local_maxima(values)
        true = array([False,True,False,False,True,False, False ,False,True, False, False, False])
        assert_array_equal(true, result)

        #test two peaks with two local minima
        values = array([1,2,3,4,5,5,2,2,5,4,5,5])
        result = smoothing_spline.find_local_maxima(values)
        true = array([False,False,False,False,True,False,False, False, True, False, True, False])
        assert_array_equal(true, result)
        
        #Test long array
        values = array([10,10,9,9,9,9,9,9,9,9,10,10])
        result = smoothing_spline.find_local_maxima(values)
        true = array([True,False,False,False,False,False,False, False, False,False, True, False])
        assert_array_equal(true, result)
        
        #test something that is breaking in reality 
        values = array([ 32,  84,  85,  85,])
        result = smoothing_spline.find_local_maxima(values)
        true = array([False, False, True, False])
        assert_array_equal(true, result)
        
        #more failing stuff
        values = array([32,  32,  32,  28,  28,  28,  28,  28,  45,  45,  57,  80,  80])
        result = smoothing_spline.find_local_maxima(values)
        true = array([False,True,False,False,False,False,False, False, False, False, False, True, False])
        assert_array_equal(true, result)
        
    def test_get_maxima_real(self):
        
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        xvals = range(20)
        data = x1 + x2

        #test
        smoothing_spline = SmoothingSpline(xvals, data)
        #more real tests
        values = array([ 201,  201,  201,  201,  201,  201,  201,  201,  201,  201,  201,  201,
          201,  201,  201,  201,  205,  196,   82,   -1])
        result = smoothing_spline.find_local_maxima(values)
        true = array([False,False,False,False,False,False,False, False, False, False, False, False, False, False, False, False, True, False, False, False])
        assert_array_equal(true, result)
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()