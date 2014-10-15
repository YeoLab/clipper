'''
Created on Jul 25, 2012

@author: gabrielp
'''
import unittest
from clipper.src.call_peak import *
from clipper.src.peaks import shuffle, find_sections
from clipper.src.readsToWiggle import readsToWiggle_pysam
from numpy import *
from numpy.testing import *
from scipy import interpolate
import os
import clipper
            
class Test(unittest.TestCase):

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
        
        
    def test_poissonP(self):
        
        """
    
        Performs unit testing on poissonP
        
        Will not test math behind calling poisson fuction more than nessessary as it has already been tested by scipy
        
        """ 
        
        #Case: fewer than 3 reads fall within a peak region: Expected result should return as though there are 3 expected reads
        result = poissonP(50, 3, 50, 2)
        self.assertAlmostEqual(result, (1 - 0.64723188878223115)) #manually calculated on personal computer stats.poisson.cdf(3, 3)
        
        #Case: More than 3 reads fall within a peak region. Expected: Result should return as normal (+1 of what ever it should be)
        result = poissonP(50, 3, 50, 4)
        self.assertAlmostEqual(result, (1 - 0.26502591529736164)) #manually calculated stats.poisson.cdf(3, 5)
    
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
        
    def test_peaks_from_info(self):
        """
        
        Tests peak_from_info function, really badly, I'm basically making this so I can
        see if there is a memory leak
        
        """
        reads = pysam.Samfile(os.path.join(clipper.test_dir(), "allup_test.bam")) 
        reads = reads.fetch(region="chr15:91536649-91537641")
        loc = ['chr15', 'bar', 91536649, 91537641, "+"]
        #wiggle, jxns, pos_counts, lengths, allreads = readsToWiggle_pysam(reads, 91537632, 91537675, '-', 'center', False)
        #result = peaks_from_info(wiggle, pos_counts,lengths,loc, 992, 25,.05, None, 3, .05, False, 10, 1000, False, .05)
        #print result

    def test_call_peaks(self):
        pass
        #reads = pysam.Samfile(os.path.join(clipper.test_dir(), "allup_test.bam")) 
        #reads = reads.fetch(region="chr15:91536649-91537641")
        #loc = ['chr15', 'bar', 91536649, 91537641, "+"]
        #result = peaks_from_info(wiggle, pos_counts,lengths,loc, 992, 25,.05, None, 3, .05, False, 10, 1000, False, .05)
                
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
