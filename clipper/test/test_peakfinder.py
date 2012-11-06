import unittest 
from clipper.src.peakfinder import *
import pkg_resources           
import pysam
import filecmp
class test_peakfinder(unittest.TestCase):
    
    parser = None
    def setUp(self):
              
        """
        
        General setup, currently creates parser for various allup tests
        
        """
        
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
        
        self.parser.add_option("--global-cutoff", dest="global_cutoff", action="store_false", help="apply global transcriptome level cutoff to CLIP-seq peaks, Default:%default", default=True, metavar="P")

        self.parser.add_option("--serial", dest="serial", action="store_true", help="run genes in sequence (not parallel)")
        self.parser.add_option("--maxgenes", dest="maxgenes", default=None, help="stop computation after this many genes, for testing", metavar="NGENES")
        self.parser.add_option("--job_name", dest="job_name", default="FAP", help="name for submitted job. Not used with --serial.  default:%default", metavar="NAME")
        self.parser.add_option("--processors", dest="np", default="autodetect", help="number of processors to use. Not used with --serial.  default:%default", type="str", metavar="NP")
        self.parser.add_option("--notify", dest="notify", default=None, help="email address to notify of start, errors and completion", metavar="EMAIL")
        self.parser.add_option("--superlocal", action = "store_true", dest="SloP", default=False, help="Use super-local p-values, counting reads in a 1KB window around peaks")
        self.parser.add_option("--color", dest="color", default="0,0,0", help="R,G,B Color for BED track output, default:black (0,0,0)")
        self.parser.add_option("--start", dest="start", default=False, action="store_true", help=SUPPRESS_HELP) #private, don't use
        self.parser.add_option("--save-pickle", dest="save_pickle", default=False, action = "store_true", help="Save a pickle file containing the analysis")
        self.parser.add_option("--debug", dest="debug", default=False, action="store_true", help="disables multipcoressing in order to get proper error tracebacks")


    
    

    def test_allup(self):
        
        """
    
        Performs basic all up test on entire program (except for main)
        
        """
        
        args = ["-b", pkg_resources.resource_filename(__name__, "../test/allup_test.bam"),
                 "-s", "hg19",
                  "-g", "ENSG00000198901", 
                  "--serial", 
                  "--job_name=peak_test",
                   "--outfile=" + os.getcwd() + "/peak_results.bed",
                   "-q"
                ]    
        (options, args) = self.parser.parse_args(args)
        main(options)
        print os.getcwd()
        tested = open(os.getcwd() + "/peak_results.bed")
        correct = open(pkg_resources.resource_filename(__name__, "../test/peak_results_no_overlap.BED"))
        
        #problem with tracks being different
        tested_tool = pybedtools.BedTool(tested)
        correct_tool = pybedtools.BedTool(correct)
        
        #checks to make sure files are equal and there are not exact dups
        self.assertAlmostEqual(len(tested_tool), len(correct_tool), delta=3)
        for test, correct in zip(tested_tool, correct_tool):
            self.assertEqual(test, correct)
        
        #cleanup
        os.remove(os.getcwd() + "/peak_results.bed")
 
    def test_check_overlaps(self):
        
        """
    
        Checks for overlapping results, we don't want this
        
        """
        
        args = ["-b", pkg_resources.resource_filename(__name__, "../test/allup_test.bam"),
                 "-s", "hg19",
                  "-g", "ENSG00000198901", 
                  "--serial", 
                  "--job_name=peak_test",
                   "--outfile=" + os.getcwd() + "/peak_results.bed",
                   "-q"
                ]    
        (options, args) = self.parser.parse_args(args)
        main(options)
        
        #tests to make sure there are no overlaps
        tested = open(os.getcwd() + "/peak_results.bed")
        tested_tool2 = pybedtools.BedTool(tested).saveas(os.getcwd() + "/foo.bed")
        result = tested_tool2.intersect(tested_tool2)
        self.assertEqual(len(result), len(tested_tool2), 
                         "there are overlaps in the output file") 
        
        #cleanup
        os.remove(os.getcwd() + "/peak_results.bed")
        os.remove(os.getcwd() + "/foo.bed")
       

    
    def test_trim_reads(self):
        
        """
    
        Performs unit tests on trim_reads
        
        """ 

        #does standard test assuming no melformed input
        test_file = pkg_resources.resource_filename(__name__, "../test/allup_test.bam")
        #print type(test_file)
        outfile = trim_reads(test_file)
        correct = pysam.Samfile(pkg_resources.resource_filename(__name__, "../test/rmdup_test.bam"))
        test = pysam.Samfile(outfile)
        
        assert filecmp.cmp(outfile, pkg_resources.resource_filename(__name__, "../test/rmdup_test.bam") )
        #for t, c in zip(correct, test):
        #    assert t == c
            

    def test_check_for_index(self):
        
        """
    
        Performs unit tests on check_for_index function
        
        """
        
        #Test if string is null, expected result is operation 
        #throws file not exist exception 
        handle = None
        self.assertRaises(TypeError, check_for_index, handle)
        
        #Test if bam file doesn't exist, expected result is operation throws 
        #file does not exist exception
        handle = "/foo/bar"
        self.assertRaises(NameError, check_for_index, handle)
        
        #Test if file is not bam, but exists expected 
        #result is to throw improper file error
        handle = pkg_resources.resource_filename(__name__, "test/test_peakfinder.py")
        self.assertRaises(NameError, check_for_index, handle)
        
        #Test if file is bam and indexed expected 
        #result is returns 1 and succedes
        #should also check if file exists, but I'm lazy
        handle = pkg_resources.resource_filename(__name__, "../test/indexed_test.bam")
        result = check_for_index(handle)
        assert result == 1
        
        #Test if file is bam and not indexed, expected 
        #result is returns one and succedes
        #should also check if file exists, but I'm lazy
        handle = pkg_resources.resource_filename(__name__, "../test/not_indexed_test.bam")
        result = check_for_index(handle)
        assert result == 1
        
        #cleanup (should be in taredown)
        os.remove(pkg_resources.resource_filename(__name__, "../test/not_indexed_test.bam.bai"))
    

    def test_build_geneinfo(self):
        self.maxDiff = 10000000
        """
    
        Performs unit testing on build_geneinfo
        
        I'm hopefully going to remove this method soon so no testing for now
        
        """ 
        
        #checks error mode
        self.assertRaises(TypeError, build_geneinfo, None)
        
        self.assertRaises(IOError, build_geneinfo, "foo")

        #checks working mode
        geneinfo = build_geneinfo(
                    clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"))

  
        true_values = {
        "ENSG00000232113" : ["chr1",  "ENSG00000232113",  173604911,      173606273, "+"],
        "ENSG00000228150" : ["chr1",  "ENSG00000228150",  10002980,       10010032,  "+"],
        "ENSG00000223883" : ["chr1",  "ENSG00000223883",  69521580,       69650686,  "+"],
        "ENSG00000135750" : ["chr1",  "ENSG00000135750",  233749749,      233808258, "+"],
        "ENSG00000227280" : ["chr1",  "ENSG00000227280",  145373053,      145375554 ,"-"],
        }


        self.assertDictEqual(geneinfo, true_values)
        
    def test_build_lengths(self):
        
        """
    
        Performs unit testing on build_lengths
        
        I'm hopefully going to remove this method soon so no unit testing for now
        
        """
        
        #Checks error mode
        self.assertRaises(ValueError, build_lengths, None)
        
        self.assertRaises(ValueError, build_lengths, clipper.test_file("foo"))
        
        #checks working mode
        lengths = build_lengths(
                    clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"))

        true = {"ENSG00000232113" : 384,
                "ENSG00000228150" : 323,
                "ENSG00000223883" : 437,
                "ENSG00000135750" : 3141,
                "ENSG00000227280" : 212,
                }       
         
        self.assertDictEqual(lengths, true)
        
    def test_add_species(self):
        
        """
    
        Performs unit testing on add_species
        
        I'll probably refactor this a bit so I won't work to hard on testing this
       
        """
        
        #Case: object is returned as expected
        result = add_species("hg19", [range(1, 22), "X", "Y"], 
                                        "foo", 
                                        "bar", 
                                        "baz")
        
        self.assertEqual(result, {"chrs" : range(1, 22) + ["X"] + ["Y"], 
                          "gene_bed" : "foo", 
                          "mRNA" : "bar", 
                          "premRNA" : "baz"})
    
    def test_get_acceptable_species(self):
        
        """
        
        Test get acceptable species 
        
        """
        
        result = get_acceptable_species()
        
        #just a quick test to make sure it works, probably need to fix this
        #later
        self.assertSetEqual(result, set(["test", "hg19", "mm9", "hg18"]))
    
      
    def test_build_transcript_data(self):
        self.maxDiff = 10000000
        """
        
        Tests building transcript data and returning the proper values
        
        Doesn't assume malformed data
        
        """
    
        #tests error modes    
        self.assertRaises(ValueError, build_transcript_data, None, None, None, None, True)
        
        self.assertRaises(ValueError, build_transcript_data, "foo", "bar", "bar", "bar", True)
        
        self.assertRaises(ValueError, build_transcript_data, "bar", None, None, None, True)
        
        #tests hg19 to make sure its equal to logic
        genes, lengths = build_transcript_data("test", None, None, None, True)
        true_genes = build_geneinfo(clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"))
        true_lengths   = build_lengths(clipper.data_file("test.AS.STRUCTURE_premRNA.lengths"))
    
        self.assertDictEqual(genes, true_genes)
        self.assertDictEqual(lengths, true_lengths)
        
        #tests hg19 on premrna lengths
        genes, lengths = build_transcript_data("test", None, None, None, False)
        true_genes = build_geneinfo(clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"))
        true_lengths  = build_lengths(clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"))
    
        self.assertDictEqual(genes, true_genes)
        self.assertDictEqual(lengths, true_lengths)
        
        #Test custom files 
        #this should all work, 
        self.assertRaises(IOError, build_transcript_data, None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), None, clipper.data_file("test.AS.STRUCTURE_premRNA.lengths"), False)
        build_transcript_data(None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"), clipper.data_file("test.AS.STRUCTURE_premRNA.lengths"), True)
        build_transcript_data(None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"), clipper.data_file("test.AS.STRUCTURE_premRNA.lengths"), False)
        build_transcript_data(None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"), None, False)
        build_transcript_data(None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), None, clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"), True)
    
    def test_transcriptome_filter(self):
        """
        
        Tests transcriptome filter
        not great tests, but good enough to make sure we don't have regressions
        
        """
        cluster = {'Nreads' : 5, "size" : 10}
        transcriptome_size = 1000
        transcriptome_reads = 10000
        poisson_cutoff = .05
        
        result = transcriptome_filter(poisson_cutoff, transcriptome_size, transcriptome_reads, cluster)
        
        self.assertFalse(result) 
        
        cluster = {'Nreads' : 10000, "size" : 100}
        transcriptome_size = 1000
        transcriptome_reads = 10000
        poisson_cutoff = .05
        
        result = transcriptome_filter(poisson_cutoff, transcriptome_size, transcriptome_reads, cluster)
        self.assertTrue(result)
        
        cluster = {'Nreads' : 0, "size" : 0}
        transcriptome_size = 0
        transcriptome_reads = 10000
        poisson_cutoff = .05
        
        result = transcriptome_filter(poisson_cutoff, transcriptome_size, transcriptome_reads, cluster)
        self.assertFalse(result)
        
    def test_filter_results(self):
        
        """
        
        tests filter results
        
        good for regression tests, need to do better verification...
        
        """
        results = [{'loc': ['chr15', 'ENSG00000198901', 91509274, 91537804, '-'], 
          'Nclusters': 24, 
          'nreads': 2086, 
          'threshold': 32, 
          'clusters': {'chr15\t1\t10\tENSG1\t3.44651351902e-09\t-\t50\t60': {'GeneP': .04, 'Nreads': 52, 'SloP': .04, 'size': 32}, 
                       'chr15\t200\t300\tENSG2\t0.0\t-\t140\t160': {'GeneP': .06, 'Nreads': 239, 'SloP': .06, 'size': 45}}}]
        

        transcriptome_size = 10000
        transcriptome_reads = 100000
        
        result = filter_results(results, .07, transcriptome_size, transcriptome_reads, False)
        self.assertSetEqual(set(['chr15\t1\t10\tENSG1\t0.04\t-\t50\t60', 'chr15\t200\t300\tENSG2\t0.06\t-\t140\t160']), result)
        #assert False
        
        result = filter_results(results, .05, transcriptome_size, transcriptome_reads, False)
        self.assertSetEqual(set(['chr15\t1\t10\tENSG1\t0.04\t-\t50\t60']), result)

        result = filter_results(results, .07, transcriptome_size, transcriptome_reads, True)
        self.assertSetEqual(set([]), result)

    def test_count_transcriptome_reads(self):
        """
        
        Tests count_transcriptome_reads
        
        """
        
        results = [
          {'loc': ['chr15', 'ENSG00000198901', 91509274, 91537804, '-'], 
          'Nclusters': 24, 
          'nreads': 200, 
          'threshold': 32, 
          'clusters': {'chr15\t1\t10\tENSG1\t3.44651351902e-09\t-\t50\t60': {'GeneP': 3.4465135190231422e-09, 'Nreads': 100, 'SloP': 3.4465135190231422e-09, 'size': 32}, 
                       'chr15\t200\t300\tENSG2\t0.0\t-\t140\t160': {'GeneP': 0.0, 'Nreads': 100, 'SloP': 0.0, 'size': 45}}}, 
          {'loc': ['chr15', 'ENSG00000198901', 91509274, 91537804, '-'], 
          'Nclusters': 24, 
          'nreads': 200, 
          'threshold': 32, 
          'clusters': {'chr15\t1\t10\tENSG1\t3.44651351902e-09\t-\t50\t60': {'GeneP': 3.4465135190231422e-09, 'Nreads': 100, 'SloP': 3.4465135190231422e-09, 'size': 32}, 
                       'chr15\t200\t300\tENSG2\t0.0\t-\t140\t160': {'GeneP': 0.0, 'Nreads': 100, 'SloP': 0.0, 'size': 45}}},  
                   ]
        result = count_transcriptome_reads(results)
        
        self.assertEqual(400, result)
        
    def test_main(self):
        
        """
    
        Performs unit testing on main
        
        Mostly testing validation and input here
        TODO: fill in test...    
        
        """
        pass
    
    def tearDown(self):
        pass
        
if __name__ == '__main__':
    unittest.main()
    os.remove(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
