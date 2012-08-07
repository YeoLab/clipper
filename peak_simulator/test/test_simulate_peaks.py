'''
Created on Aug 6, 2012

@author: gabrielp
'''
import unittest
from src.simulate_peaks import distribute_reads, output_bam, output_bed
from numpy import array, ones
from numpy.testing import *

class Test(unittest.TestCase):


    def test_distribute_reads(self):
        """
        
        Tests distribute reads function
        
        """
        
        weights = array([1,1,1,1])
        num_reads = 4
        result = distribute_reads(weights, num_reads)
        assert_array_equal(result, array([1,1,1,1]))
        
        weights = array([2,2,2,2])
        num_reads = 4
        result = distribute_reads(weights, num_reads)
        assert_array_equal(result, array([1,1,1,1]))
        
        weights = array([2,2,0,0])
        num_reads = 4
        result = distribute_reads(weights, num_reads)
        assert_array_equal(result, array([2,2,0,0]))
        
        weights = array([1,2,0,0])
        num_reads = 4
        result = distribute_reads(weights, num_reads)
        assert_array_equal(result, array([1,3,0,0]))
        
        #This test sort of breaks, need to assign exactly the number of reads
        #Need to figure out how to assign exactly correct numner of reads
        weights = array([1,2,1,1])
        num_reads = 4
        result = distribute_reads(weights, num_reads)
        assert_array_equal(result, array([1,2,1,1]))
        
    def test_output_bam(self):
        
        """
        
        tests output function, a bit difficult to test at the moment because I'm not sure what the 
        output should look like
        
        """
        
        reads = ones(60000)
        read_length = 50
        output_name = "foo.bam"
        output_bam(reads, read_length, output_name)
    
    def test_output_bed(self):
        peaks = ((1,2), (5,10))
        output_name = "foo.bed"
        
        output_bed(peaks, output_name)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
    