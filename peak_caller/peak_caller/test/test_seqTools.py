'''
Created on Jul 20, 2012

@author: gabrielp
'''
import unittest
from peak_caller.src.seqTools import *
import pysam
class Test(unittest.TestCase):


    def test_readsToWiggle_pysam(self):

        reads = pysam.Samfile(pkg_resources.resource_filename(__name__, "../test/allup_test.bam"))
        reads = reads.fetch(region="chr15:91536649-91537641")
        wiggle, jxn, pos_counts, lengths, allreads = readsToWiggle_pysam(reads, 91537632, 91537675, "-", 'center')
        
        wiggle_true = [  2. ,  2.,   2. ,  2. ,  2. ,  2.  , 2. ,  2. , 11. , 11.,  11. , 11.  ,11. , 11. , 11.,
  11. , 11.,  11.,  11. , 11.  ,11. , 11. , 11. , 11.,  11. , 11. , 11.  ,11. , 11.  ,11.,
  11. , 11.,  11.,   9. ,  9. ,  9. ,  9. ,  9.,   9. ,  9.,   9. ,  0. ,  0.,   0.]
        
        for true, test in zip(wiggle_true, wiggle):
            assert true == test
        
        pos_counts_true = [ 0. , 0.,  0. , 0.  ,0. , 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0. , 0. , 
                           0. , 0. , 2.,  0., 0. , 0.,  0.,  0.,  0. , 0.,  9.,  0. , 0.,  0. , 0. ,  
                           0. , 0. , 0. , 0. , 0.,  0.,  0., 0. , 0.,  0. , 0. , 0.,  0.,  0. ,  0.]
        
        for true, test in zip(pos_counts_true, pos_counts):
            assert test == true
            
        assert lengths == [33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33]
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()