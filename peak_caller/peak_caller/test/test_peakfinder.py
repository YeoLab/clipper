
from peak_caller.src.peakfinder import *
from optparse import OptionParser, SUPPRESS_HELP
import os
import unittest 
import scipy

class test_peakfinder (unittest.TestCase):
    
    def setUp(self):
        pass
    
    def test_allup(self):
        
        usage="\npython peakfinder.py -b <bamfile> -s <hg18/hg19/mm9>\n OR \npython peakfinder.py -b <bamfile> --customBED <BEDfile> --customMRNA <mRNA lengths> --customPREMRNA <premRNA lengths>"
        description="Fitted Accurate Peaks. Michael Lovci 2012. CLIP peakfinder that uses fitted smoothing splines to define clusters of binding.  Computation is performed in parallel using MPI.  You may decide to use the pre-built gene lists by setting the --species parameter or alternatively you can define your own list of genes to test in BED6/12 format and provide files containing the length of each gene in PREMRNA form and MRNA form (both are required). Questions should be directed to michaeltlovci@gmail.com."
        parser = OptionParser(usage=usage, description=description)
    
        parser.add_option("--bam", "-b", dest="bam", help="A bam file to call peaks on", type="string", metavar="FILE.bam")
    
        parser.add_option("--species", "-s", dest="species", help="A species for your peak-finding")
    
        parser.add_option("--customBED", dest="geneBEDfile", help="bed file to call peaks on, must come withOUT species and with customMRNA and customPREMRNA", metavar="BEDFILE")
        parser.add_option("--customMRNA", dest="geneMRNAfile", help="file with mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")
        parser.add_option("--customPREMRNA", dest="genePREMRNAfile", help="file with pre-mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")
    
        parser.add_option("--outdir", dest="prefix", default=os.getcwd(), help="output directory, default=cwd")    
        parser.add_option("--outfile", dest="outfile", default="fitted_clusters", help="a bed file output, default:%default")
    
        parser.add_option("--gene", "-g", dest="gene", action="append", help="A specific gene you'd like try", metavar="GENENAME")
        parser.add_option("--plot", "-p", dest="plotit", action="store_true", help="make figures of the fits", default=False)
        parser.add_option("--quiet", "-q", dest="quiet", action="store_true", help="suppress notifications")
    
        parser.add_option("--minreads", dest="minreads", help="minimum reads required for a section to start the fitting process.  Default:%default", default=3, type="int", metavar="NREADS")
        parser.add_option("--margin", dest="margin", type="int", help="find sections of genes within M bases that have genes and perform fitting. Default:%default", default=15, metavar="NBASES")
        parser.add_option("--trim", dest="trim", action="store_true", default=False, help="Trim reads with the same start/stop to count as 1")
        parser.add_option("--premRNA", dest="premRNA", action="store_true", help="use premRNA length cutoff, default:%default", default=False)
        parser.add_option("--poisson-cutoff", dest="poisson_cutoff", type="float", help="p-value cutoff for poisson test, Default:%default", default=0.05, metavar="P")
        parser.add_option("--FDR", dest="FDR_alpha", type="float", default=0.05, help="FDR cutoff for significant height estimation, default=%default")
        parser.add_option("--threshold", dest="threshold", type="int", default=None, help="Skip FDR calculation and set a threshold yourself")
    
        parser.add_option("--serial", dest="serial", action="store_true", help="run genes in sequence (not parallel)")
        parser.add_option("--maxgenes", dest="maxgenes", default=None, help="stop computation after this many genes, for testing", metavar="NGENES")
        parser.add_option("--job_name", dest="job_name", default="FAP", help="name for submitted job. Not used with --serial.  default:%default", metavar="NAME")
        parser.add_option("--processors", dest="np", default=32, help="number of processors to use. Not used with --serial.  default:%default", type="int", metavar="NP")
        parser.add_option("--notify", dest="notify", default=None, help="email address to notify of start, errors and completion", metavar="EMAIL")
        parser.add_option("--superlocal", action = "store_true", dest="SloP", default=False, help="Use super-local p-values, counting reads in a 1KB window around peaks")
        parser.add_option("--color", dest="color", default="0,0,0", help="R,G,B Color for BED track output, default:black (0,0,0)")
        parser.add_option("--start", dest="start", default=False, action="store_true", help=SUPPRESS_HELP) #private, don't use
        parser.add_option("--save-pickle", dest="save_pickle", default=False, action = "store_true", help="Save a pickle file containing the analysis")
        
        args = ["-b", "/nas3/gpratt/HEK293/GSM782786_NAMF-mapped.hg19.bam.sorted.bam",
                 "-s", "hg19",
                  "-g", "ENSG00000198901", 
                  "--serial", 
                  "--job_name=peak_test",
                   "--outfile=peak_results",
                   "-q"
                ]    
        (options,args) = parser.parse_args(args)
        main(options)
        tested = open("../src/peak_results.BED")
        correct = open("peak_results.BED")
        for test, correct in zip(tested, correct):
            assert test == correct
    
    """
    
    Performs unit tests on check_for_index function
    
    """
    def test_check_for_index(self):
        #Test if string is null, expected result is operation throws file not exist exception 
        f = None
        self.assertRaises(TypeError, check_for_index, f)
        
        #Test if bam file doesn't exist, expected result is operation throws 
        #file does not exist exception
        f = "/foo/bar"
        self.assertRaises(NameError, check_for_index, f)
        
        #Test if file is not bam, but exists expected result is to throw improper file error
        f = os.path.curdir + "/test_peakfinder.py"
        self.assertRaises(NameError, check_for_index, f)
        
        #Test if file is bam and indexed expected result is returns 1 and succedes
        #should also check if file exists, but I'm lazy
        f = os.path.curdir + "/indexed_test.bam"
        result = check_for_index(f)
        assert result == 1
        
        #Test if file is bam and not indexed, expected result is returns one and succedes
        #should also check if file exists, but I'm lazy
        f = os.path.curdir + "/not_indexed_test.bam"
        result = check_for_index(f)
        assert result == 1
        
        #cleanup (should be in taredown)
        os.remove(os.path.curdir + "/not_indexed_test.bam.bai")
    
    """
    
    Performs unit tests on get_FDR_cutoff_mode function
    Function is currently deperacated, tests not done
    
    """
    def test_get_FDR_cutoff_mode(self):
        pass
    
    """
    
    Performs unit tests on get_FDR_cutoff_mean function
    
    Difficult to test because of random sampling
    
    TODO: Figure out how to manually calculate FDR cutoff, right now 
    I just took the old result and am using that.  Math appears to be correct,
    but this is a really bad practice
    
    """
    def test_get_FDR_cutoff_mean(self):
        
        
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
        assert result == 7
        
        #Second similar case 
        read_lengths = [30] * 20
        result = get_FDR_cutoff_mean(read_lengths, 100)
        assert result == 11 
       
        
        #Case: enough reads, mean is lower than minimum cutoff, expected result: minimum cutoff is returned
        result = get_FDR_cutoff_mean(read_lengths, 100, mincut = 20)
        assert result == 20
    
    """
    
    Performs unit testing on build_geneinfo
    
    I'm hopefully going to remove this method soon so no testing for now
    
    """ 
    def test_build_geneinfo(self):
        pass
    
    """
    
    Performs unit testing on build_lengths
    
    I'm hopefully going to remove this method soon so no unit testing for now
    
    """
    def test_build_lengths(self):
        pass
    
    """
    
    Performs unit testing on find_sections

    """
    def test_find_sections(self):
        #setup 
        
        #Null Case
        self.assertRaises(TypeError, find_sections, (None, 0))
        
        #Case with all zero coverage
        wiggle = [0] * 20
        result = find_sections(wiggle, 0)
        assert result == []
        
        #Case with all non-zero coverage
        wiggle = [5] * 20
        result = find_sections(wiggle, 0)
        assert result == ["0|19"]
      
        #Case with one region on margin of one and two regions on margin of two
        
        #returns two segnments
        wiggle = ([5] * 20) + [0] + ([5] * 20)
        result = find_sections(wiggle, 0)
        print result
        assert result == ["0|19", "21|40"]
        
        #returns one segnment
        result = find_sections(wiggle, 1)
        print result
        assert result == ["0|40"]
        
        #second case returns two segnments
        wiggle = ([5] * 9) + [0] + ([5] * 10)
        result = find_sections(wiggle, 0)
        print result
        assert result == ["0|8", "10|19"]
        
        #returns one segnment
        result = find_sections(wiggle, 1)
        assert result == ["0|19"]
        
        #Edge case where margins stop before the end of genes
        wiggle = [0] + ([5] * 10)
        result = find_sections(wiggle, 0)
        assert result == ["1|10"]
        
        #Edge case where margins start after the start of genes
        wiggle = ([5] * 10) + [0] 
        result = find_sections(wiggle, 0)
        assert result == ["0|9"]
    
    """
    
    Quality control function, not important to main running of program
    Not tested
    
    """
    def test_plotSpline(self):
        pass
    
    """
    
    Performs unit testing on find_univariateSpline
    As this is mostly a wrapper for a scipy function I will not test spline calling
    beyond basics.
    
    These tests will verify that the resid logic works, and eventually these 
    will be factored into two functions.  Its bad form to have to toggle your return
    value in the function
    
    """
    def test_find_univariateSpline(self):
        
        #Case null inputs, expected:  Everything goes to hell, not testing
        
        #Case resid is false, expected: returns univariateSpline spline, just verifies that it is the same as 
        #the scipy result
        
        #setup 
        x1 = range(10)
        x2 = range(10)
        x2.reverse()
        x3 = x1 + x2
        y = x3
        smoothing = 5 
        #expected
        expected = scipy.interpolate.UnivariateSpline(x3, y, k=3, s=smoothing)
        
        #test
        result = find_univariateSpline(smoothing, x3, y, 3, resid=False)
        
        #hacky test, but they should be about the same
        self.assertAlmostEqual(expected.get_residual(), result.get_residual()) 
        
        #Case resid is true. Expected: returns residual sum of squared error between the 
        #spline and the actual curve 
        
    """
    
    Performs unit testing on poissonP
    
    """ 
    def test_poissonP(self):
        pass
    
    """
    
    Performs unit testing on call_peaks
    
    will not test peaks_from_info here, just the error handling of call peaks
    Need to create a dummy call peaks function or stub or something to keep everything from
    getting called
    
    """ 
    def test_call_peaks(self):
        pass
    
    """
    
    Performs unit testing on peaks_from_info
    
    Testing of this is currently handled in the allup test
    This method is currently to complicated to test as is
    
    """
    def test_peaks_from_info(self):
        pass
    
    """
    
    Performs unit testing on add_species
        
    """
    def test_add_species(self):
        pass
    
    """
    
    Performs unit testing on main
    
    Mostly testing validation and input here
        
    """
    def test_main(self):
        pass

if __name__ =='__main__':
    unittest.main()
