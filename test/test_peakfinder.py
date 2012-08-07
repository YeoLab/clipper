
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
        self.parser.add_option("--processors", dest="np", default="autodetect", help="number of processors to use. Not used with --serial.  default:%default", type="str", metavar="NP")
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
                   "--outfile=" + os.getcwd() + "/peak_results",
                   "-q"
                ]    
        (options,args) = self.parser.parse_args(args)
        main(options)
        print os.getcwd()
        tested = open(os.getcwd() + "/peak_results.BED")
        correct = open(pkg_resources.resource_filename(__name__, "../test/peak_results_no_overlap.BED"))
        
        #problem with tracks being different
        tested_tool = pybedtools.BedTool(tested)
        correct_tool = pybedtools.BedTool(correct)
        
        #checks to make sure files are equal and there are not exact dups
        self.assertEqual(len(tested_tool), len(correct_tool))
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
                   "--outfile=" + os.getcwd() + "/peak_results",
                   "-q"
                ]    
        (options,args) = self.parser.parse_args(args)
        main(options)
        
        #tests to make sure there are no overlaps
        tested = open(os.getcwd() + "/peak_results.BED")
        tested_tool2 = pybedtools.BedTool(tested).saveas(os.getcwd() + "/foo.bed")
        result = tested_tool2.intersect(tested_tool2)
        self.assertEqual(len(result), len(tested_tool2), "there are overlaps in the output file") 
        
        #cleanup
        #os.remove(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
       
       
    """
    
    Performs unit tests on trim_reads
    
    """ 
    
    def test_trim_reads(self):
        pass
        #does standard test assuming no melformed input
        #test_file = pkg_resources.resource_filename(__name__, "../test/allup_test.bam")
        #print type(test_file)
        #outfile = trim_reads(test_file)
        #correct = pysam.Samfile(pkg_resources.resource_filename(__name__, "../test/rmdup_test.bam"))
        #test = pysam.Samfile(outfile)
        
        #for t, c in zip(correct, test):
        #    assert t == c
            
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
