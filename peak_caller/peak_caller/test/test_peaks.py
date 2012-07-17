'''
Created on Jul 17, 2012

@author: gabrielp
'''
import unittest
import peaks

class Test(unittest.TestCase):

    """
    
    Tests shuffle extention function, mostly checks for error handling due to random nature of the algoritm
    I don't test output throughly, just verify that the correct number of results appears
    TODO: fails on uniform case need to fix
    TODO: fails on small inputs.  Need to either thorw errors or handle this better
    """
    def test_shuffle(self):
        
        #Case: fail on null inputs
        self.assertRaises(TypeError, peaks.shuffle, (None, 1, 0, .05, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, None, 0, .05, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, 1, None, .05, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, 1, 0, None, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, 1, 0, .05, None))
            
        #Case: fail on zero input for [] for the reads
        #result = peaks.shuffle(1,1,0,.05, [])
        #self.assertEqual(result, [0] * 1000)
        
        #case fail on zero input for either length or #iterations
        self.assertRaises(TypeError, peaks.shuffle, (0, 1, 0, .05, [2,3,4]))
        self.assertRaises(TypeError, peaks.shuffle, (1, 0, 0, .05, [2,3,4]))
        
        #case succede and check results (need to figure how to lock down random for testing
        result = peaks.shuffle(100, 3, 0,.05, [5] * 50 )
        self.assertEqual(sum(result), 3)
        
        #makes sure it works on edge cases
        #result = peaks.shuffle(100, 3, 0, .05, [2,3,4])
        #self.assertEqual(sum(result), 3)
        
        #reads longer than gene
        self.assertRaises(TypeError, peaks.shuffle, (1, 1, 0, .05, [2,3,4]))
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()