
from peak_caller.src.peakfinder import *
from optparse import OptionParser, SUPPRESS_HELP
import os
import unittest 
import scipy
import pkg_resources
import pybedtools

class test_peakfinder(unittest.TestCase):
    
    parser = None
    
    """
    
    General setup, currently creates parser for various allup tests
    
    """
    def setUp(self):
        usage="\npython peakfinder.py -b <bamfile> -s <hg18/hg19/mm9>\n OR \npython peakfinder.py -b <bamfile> --customBED <BEDfile> --customMRNA <mRNA lengths> --customPREMRNA <premRNA lengths>"
        description="Fitted Accurate Peaks. Michael Lovci 2012. CLIP peakfinder that uses fitted smoothing splines to define clusters of binding.  Computation is performed in parallel using MPI.  You may decide to use the pre-built gene lists by setting the --species parameter or alternatively you can define your own list of genes to test in BED6/12 format and provide files containing the length of each gene in PREMRNA form and MRNA form (both are required). Questions should be directed to michaeltlovci@gmail.com."
        self.parser = OptionParser(usage=usage, description=description)
    
        self.parser.add_option("--bam", "-b", dest="bam", help="A bam file to call peaks on", type="string", metavar="FILE.bam")
    
        self.parser.add_option("--species", "-s", dest="species", help="A species for your peak-finding")
    
        self.parser.add_option("--customBED", dest="geneBEDfile", help="bed file to call peaks on, must come withOUT species and with customMRNA and customPREMRNA", metavar="BEDFILE")
        self.parser.add_option("--customMRNA", dest="geneMRNAfile", help="file with mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")
        self.parser.add_option("--customPREMRNA", dest="genePREMRNAfile", help="file with pre-mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")
    
        self.parser.add_option("--outdir", dest="prefix", default=os.getcwd(), help="output directory, default=cwd")    
        self.parser.add_option("--outfile", dest="outfile", default="fitted_clusters", help="a bed file output, default:%default")
    
        self.parser.add_option("--gene", "-g", dest="gene", action="append", help="A specific gene you'd like try", metavar="GENENAME")
        self.parser.add_option("--plot", "-p", dest="plotit", action="store_true", help="make figures of the fits", default=False)
        self.parser.add_option("--quiet", "-q", dest="quiet", action="store_true", help="suppress notifications")
    
        self.parser.add_option("--minreads", dest="minreads", help="minimum reads required for a section to start the fitting process.  Default:%default", default=3, type="int", metavar="NREADS")
        self.parser.add_option("--margin", dest="margin", type="int", help="find sections of genes within M bases that have genes and perform fitting. Default:%default", default=15, metavar="NBASES")
        self.parser.add_option("--trim", dest="trim", action="store_true", default=False, help="Trim reads with the same start/stop to count as 1")
        self.parser.add_option("--premRNA", dest="premRNA", action="store_true", help="use premRNA length cutoff, default:%default", default=False)
        self.parser.add_option("--poisson-cutoff", dest="poisson_cutoff", type="float", help="p-value cutoff for poisson test, Default:%default", default=0.05, metavar="P")
        self.parser.add_option("--FDR", dest="FDR_alpha", type="float", default=0.05, help="FDR cutoff for significant height estimation, default=%default")
        self.parser.add_option("--threshold", dest="threshold", type="int", default=None, help="Skip FDR calculation and set a threshold yourself")
    
        self.parser.add_option("--serial", dest="serial", action="store_true", help="run genes in sequence (not parallel)")
        self.parser.add_option("--maxgenes", dest="maxgenes", default=None, help="stop computation after this many genes, for testing", metavar="NGENES")
        self.parser.add_option("--job_name", dest="job_name", default="FAP", help="name for submitted job. Not used with --serial.  default:%default", metavar="NAME")
        self.parser.add_option("--processors", dest="np", default=32, help="number of processors to use. Not used with --serial.  default:%default", type="int", metavar="NP")
        self.parser.add_option("--notify", dest="notify", default=None, help="email address to notify of start, errors and completion", metavar="EMAIL")
        self.parser.add_option("--superlocal", action = "store_true", dest="SloP", default=False, help="Use super-local p-values, counting reads in a 1KB window around peaks")
        self.parser.add_option("--color", dest="color", default="0,0,0", help="R,G,B Color for BED track output, default:black (0,0,0)")
        self.parser.add_option("--start", dest="start", default=False, action="store_true", help=SUPPRESS_HELP) #private, don't use
        self.parser.add_option("--save-pickle", dest="save_pickle", default=False, action = "store_true", help="Save a pickle file containing the analysis")
        

    
    
    """
    
    Performs basic all up test on entire program (except for main)
    
    """
    def test_allup(self):
        
        
        args = ["-b", pkg_resources.resource_filename(__name__, "../test/allup_test.bam"),
                 "-s", "hg19",
                  "-g", "ENSG00000198901", 
                  "--serial", 
                  "--job_name=peak_test",
                   "--outfile=" + pkg_resources.resource_filename(__name__, "../src/peak_results"),
                   "-q"
                ]    
        (options,args) = self.parser.parse_args(args)
        main(options)
        tested = open(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
        correct = open(pkg_resources.resource_filename(__name__, "../test/peak_results.BED"))
        
        #problem with tracks being different
        tested_tool = pybedtools.BedTool(tested)
        correct_tool = pybedtools.BedTool(correct)
        
        #checks to make sure files are equal and there are not exact dups
        for test, correct in zip(tested_tool, correct_tool):
            self.assertEqual(test, correct)
        
        #cleanup
        #os.remove(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
    
    """
    def test_plotting(self):
        args = ["-b", pkg_resources.resource_filename(__name__, "../test/allup_test.bam"),
                 "-s", "hg19",
                  "-g", "ENSG00000198901", 
                  "--serial", 
                  "--job_name=peak_test",
                   "--outfile=" + pkg_resources.resource_filename(__name__, "../src/peak_results"),
                   "-q",
                   "-p",
                ]
               
        (options,args) = self.parser.parse_args(args)
        main(options)

        tested = open(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
        correct = open(pkg_resources.resource_filename(__name__, "../test/peak_results.BED"))
        
        #problem with tracks being different
        tested.next()
        correct.next()
        for test, correct in zip(tested, correct):
            self.assertEqual(test, correct)
        
        #cleanup
        os.remove(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
        """
        
    """
    
    Checks for overlapping results, we don't want this
    
    """
    def test_checkOverlaps(self):
        
        args = ["-b", pkg_resources.resource_filename(__name__, "../test/allup_test.bam"),
                 "-s", "hg19",
                  "-g", "ENSG00000198901", 
                  "--serial", 
                  "--job_name=peak_test",
                   "--outfile=" + pkg_resources.resource_filename(__name__, "../src/peak_results"),
                   "-q"
                ]    
        (options,args) = self.parser.parse_args(args)
        main(options)
        
        #tests to make sure there are no overlaps
        tested = open(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
        tested_tool2 = pybedtools.BedTool(tested)
        result = tested_tool2.merge(n=True)
        self.assertEqual(result, 0, "there are overlaps in the output file") 
        
        #cleanup
        os.remove(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
        
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
        f = pkg_resources.resource_filename(__name__, "test/test_peakfinder.py")
        self.assertRaises(NameError, check_for_index, f)
        
        #Test if file is bam and indexed expected result is returns 1 and succedes
        #should also check if file exists, but I'm lazy
        f = pkg_resources.resource_filename(__name__, "../test/indexed_test.bam")
        result = check_for_index(f)
        assert result == 1
        
        #Test if file is bam and not indexed, expected result is returns one and succedes
        #should also check if file exists, but I'm lazy
        f = pkg_resources.resource_filename(__name__, "../test/not_indexed_test.bam")
        result = check_for_index(f)
        assert result == 1
        
        #cleanup (should be in taredown)
        os.remove(pkg_resources.resource_filename(__name__, "../test/not_indexed_test.bam.bai"))
    
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
    
    Quality control function, not important to main running of program
    Not tested
    
    """
    def test_plotSpline(self):
        pass
        #plotSections([5] * 10, ["1|5", "7|9"], 3)
        #assert 1 == 0
    
    """
    
    tests plotSections function, appears to have computer specific issues
    
    """
    def test_plotSections(self):
        """

        Plots each section individually, I think
        Wiggle is a list representing a wiggle track
        sections is a list of strings of format "start|stop" where start and stop are both integers
        threshold is an integer 

        """     
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
        #TODO: write better tests, this only tests very basic case
        residual = find_univariateSpline(smoothing, x3, y, 3, resid=True)
        self.assertAlmostEqual(expected.get_residual(), result.get_residual()) 
        
        #this works because the number of turns is 1 so the residual * turns is the residual
        self.assertAlmostEqual(residual, sqrt(expected.get_residual()))
        
        #Case: calculation of univarate spline breaks
        #There should be more error handling this is a general catch, we should be more sensetavie
        #to different types of errors
        assert Inf == find_univariateSpline(None, None, None, None, resid=True)
        
    """
    
    Performs unit testing on poissonP
    
    Will not test math behind calling poisson fuction more than nessessary as it has already been tested by scipy
    
    """ 
    def test_poissonP(self):
        #Case: fewer than 3 reads fall within a peak region: Expected result should return as though there are 3 expected reads
        result = poissonP(50, 3, 50, 2)
        self.assertAlmostEqual(result, (1 - 0.64723188878223115)) #manually calculated on personal computer stats.poisson.cdf(3, 3)
        
        #Case: More than 3 reads fall within a peak region. Expected: Result should return as normal
        result = poissonP(50, 3, 50, 4)
        self.assertAlmostEqual(result, (1 - 0.4334701203667089)) #manually calculated stats.poisson.cdf(3, 4)
    
        #Case: poissonp throws an exception. Expected: Result should return 1
        result = poissonP(None, None, None, None)    
        assert 1 == result
    
    """
    
    Performs unit testing on call_peaks
    
    will not test peaks_from_info here, just the error handling of call peaks
    Need to create a dummy call peaks function or stub or something to keep everything from
    getting called
    
    The way this is currently written it is difficult to test logic other than error checking
    
    """ 
    def test_call_peaks(self):
        #peaks_from_info = poissonP
        """call_peaks("chr1|foo|5|10|-", 
                   50, 
                   bam_fileobj=None, 
                   bam_file=None)"""
        assert 1 == 1
        #Case: Bam file is passed
        
        #Case: Bam object is passed
        
        #Case: Both bam file and object are passed
    
    """
    
    Performs unit testing on peaks_from_info
    
    Testing of this is currently handled in the allup test
    This method is currently to complicated to test as is
    
    """
    def test_peaks_from_info(self):
        pass
    
    """
    
    Performs unit testing on add_species
    
    I'll probably refactor this a bit so I won't work to hard on testing this
    """
    def test_add_species(self):
        #Case: object is returned as expected
        result = add_species("hg19", [range(1,22), "X", "Y"],
                                        "foo",
                                        "bar",
                                        "baz")
        
        assert result == {"chrs" : range(1,22) + ["X"] + ["Y"],
                          "gene_bed" : "foo",
                          "mRNA" : "bar",
                          "premRNA" : "baz"}
                        
    """
    
    Performs unit testing on main
    
    Mostly testing validation and input here
        
    """
    def test_main(self):
        pass
    
    def tearDown(self):
        pass
        
if __name__ =='__main__':
    unittest.main()
    os.remove(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
