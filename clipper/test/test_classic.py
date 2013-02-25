'''
Created on Feb 25, 2013

@author: gabrielp
'''
import unittest

from clipper.src.call_peak import * 
from numpy import arange
class Test(unittest.TestCase):


    def test_complete_peaks(self):
        
        """
        
        Tests to see if entire section is called as a peak
        
        """
        
        data = [1,1,1,1,1,1,1,1,1,1] #len of 10
        xvals = arange(len(data))
        
        fitter = Classic(xvals, data, max_width=20, min_width=3, max_gap=10)
                         
        
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(0,10, 5)])
        
    def test_large_peaks(self):
        
        """
        
        Tests splitting up a section into two large peaks
        
        """

        
        
        data = [1,1,1,1,1,1,1,1,1,1] #len of 10
        xvals = arange(len(data))
        fitter = Classic(xvals, data, max_width=5, min_width=3, max_gap=10)
                         
        
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(0,5, 2), (5,10, 7)])
    
    def test_large_and_small_peaks(self):
        
        """
        
        Tests calling 3 peaks, two too large, the last one at the end being to small
        
        """
        
        data = [1,1,1,1,1,1,1,1,1,1,1] #len of 11
        xvals = arange(len(data))
        fitter = Classic(xvals, data, max_width=5, min_width=3, max_gap=10)
                         
        
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(0,5, 2), (5,10, 7), (10, 13, 11 )])
    def test_small_peak(self):
        
        """
        
        Tests setting one small peak
        
        """
        
        data = [1,1,1,1,1,1,1,1,1,1] #len of 11
        xvals = arange(len(data))
        
        fitter = Classic(xvals, data, max_width=20, min_width=15, max_gap=10)
                         
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(0,15, 7)])
        
    def test_gap_peaks(self):
        
        """
        
        Tests finding gaps between peaks
        
        """
        
        data = [1,1,1,1,0,0,1,1,1,1] #len of 11
        xvals = arange(len(data))
        
        fitter = Classic(xvals, data, max_width=20, min_width=1, max_gap=1)
                         
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(0,4, 2), (6,10, 8)])
    
    def test_gap_and_small_peak(self):
        
        """
        
        Tests finding gaps between peaks with a small peak before the gap
        
        """
        
        data = [1,1,1,0,0,0,1,1,1,1] #len of 11
        xvals = arange(len(data))
        
        fitter = Classic(xvals, data, max_width=20, min_width=4, max_gap=1)
                         
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(0,4,2), (6,10,8)])
    
    def test_gap_and_small_peak(self):
        
        """
        
        Tests 3 gaps between peaks with a small peak before the gap
        
        """
        
        data = [1,1,0,1,1,0,1,1,1,1] #len of 11
        xvals = arange(len(data))
        
        fitter = Classic(xvals, data, max_width=20, min_width=1, max_gap=0)
                         
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(0,2,1), (3,5,4), (6,10,8)])
    
    def test_leading_zero(self):
        
        """
        
        tests leading zeros
        
        """
        
        data = [0,0,1,1,1,1,1,1,1,1] #len of 11
        xvals = arange(len(data))
        
        fitter = Classic(xvals, data, max_width=20, min_width=1, max_gap=1)
                         
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(2,10, 6)])
    
    def test_traling_zero(self):
        
        """
        
        tests traling zeros
        
        """
        
        data = [1,1,1,1,1,1,1,1,0,0] #len of 11
        xvals = arange(len(data))
        
        fitter = Classic(xvals, data, max_width=20, min_width=1, max_gap=1)
                         
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(0,8, 4)])  
        
    def test_gap_and_small_peaks(self):
        
        """
        
        Tests finding gaps and fixing small peaks
        
        """
        
        data = [1,1,1,0,0,0,1,1,1,1] #len of 11
        xvals = arange(len(data))
        
        fitter = Classic(xvals, data, max_width=5, min_width=5, max_gap=1)
                         
        starts_and_stops = fitter.peaks()
        self.assertListEqual(starts_and_stops, [(0,5,2), (6,11,8)])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_peaks']
    unittest.main()