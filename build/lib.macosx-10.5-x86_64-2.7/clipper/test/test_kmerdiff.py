'''
Created on Feb 13, 2013

Unit tests for kmerdiff
@author: gabrielp
'''
import unittest
import clipper
from kmerdiff import *
from collections import namedtuple
class Test(unittest.TestCase):

    def test_kmer_diff(self):
        
        """
        
        Tests all of kmer diff 
        (note checked for identity to original perl module seperatly)
        
        """
        
        Motif = namedtuple('Motif', ['freq1', 'freq2', 'delta'])
        
        
        file1 = clipper.test_file("compseq_test.fasta")
        file2 = clipper.test_file("compseq_test_2.fasta")
        result, n1, n2 = kmer_diff(file1, file2, 2)
        
        self.assertEqual(n1, 11.0)
        self.assertEqual(n2, 11.0)
        
        true =  {                       
                  'AA' : Motif(2.0, 1.0, 0.621260744197396 ),
                  'AC' : Motif(0.0, 0.0, 0.0),
                  'AG' : Motif(0.0, 0.0, 0.0),
                  'AT' : Motif(1.0, 2.0, -0.621260744197396),
                  'CA' : Motif(0.0, 0.0, 0.0),
                  'CC' : Motif(2.0, 2.0, 0.0),
                  'CG' : Motif(0.0, 0.0, 0.0),
                  'CT' : Motif(0.0, 0.0, 0.0),
                  'GA' : Motif(0.0, 0.0, 0.0),
                  'GC' : Motif(1.0, 1.0, 0.0),
                  'GG' : Motif(2.0, 2.0, 0.0),
                  'GT' : Motif(0.0, 0.0, 0.0),
                  'TA' : Motif(0.0, 1.0, -1.02353263143832), 
                  'TC' : Motif(0.0, 0.0, 0.0),
                  'TG' : Motif(1.0, 1.0, 0.0),
                  'TT' : Motif(2.0, 1.0, 0.621260744197396),
                 }
        #almost equal dict assertion
        self.assertSetEqual(set(result.keys()), set(true.keys()))
        
        for key in result.keys():
            self.assertEqual(true[key].freq1, result[key].freq1, "failed at %s" % (key))
            self.assertEqual(true[key].freq2, result[key].freq2, "failed at %s" % (key))
            self.assertAlmostEqual(true[key].delta, result[key].delta, delta=3)
            
    def test_file_not_found(self):
        
        """
        
        Makes sure kmer diff fails gracefully when a file isn't found
        
        """
        
        file1 = "foo"
        file2 = "bizbuz"
        
        self.assertRaises(IOError, kmer_diff, file1, file2, 10)
        
    def test_parse_compseq(self):
        
        """
        
        Tests the resut of parse compseq 
        
        """
        
        file = clipper.test_file("compseq_test.out")
        
        diff, dict = parse_compseq(file)
        
        self.assertEqual(diff, 11.0)
        self.assertDictEqual(dict, {
                                    'AA' : 2.0,
                                    'AC' : 0.0,
                                    'AG' : 0.0,
                                    'AT' : 1.0,
                                    'CA' : 0.0,
                                    'CC' : 2.0,
                                    'CG' : 0.0,
                                    'CT' : 0.0,
                                    'GA' : 0.0,
                                    'GC' : 1.0,
                                    'GG' : 2.0,
                                    'GT' : 0.0,
                                    'TA' : 0.0,
                                    'TC' : 0.0,
                                    'TG' : 1.0,
                                    'TT' : 2.0,
                                    })

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()