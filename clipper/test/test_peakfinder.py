import unittest
import pkg_resources
import pybedtools
# class Test(unittest.TestCase):
    

    
    # def test_classic_allup(self):
    #
    #     """
    #     DON't NEED TO TEST SINCE OVERRIDE method option with spline
    #     Performs basic all up test on entire program (using classic algorithm) (except for main)
    #
    #     """
    #
    #     #self.assertTrue(False, "test is currently disabled output from logging causes it to crash")
    #     args = ["-b", clipper.test_file("allup_test.bam"),
    #              "-s", "hg19",
    #              "-g", "ENSG00000198901",
    #              "--outfile=" + os.getcwd() + "/allup_peak_results_classic.bed",
    #              "--debug",
    #              "--algorithm=classic"
    #             ]
    #
    #     (options, args) = self.parser.parse_args(args)
    #     options = override_options(options)
    #
    #
    #     main(options)
    #
    #     tested = open(os.getcwd() + "/allup_peak_results_classic.bed")
    #     correct = open(clipper.test_file("peak_results_no_overlap.BED"))
    #
    #
    #     #problem with tracks being different
    #     tested_tool = pybedtools.BedTool(tested)
    #     correct_tool = pybedtools.BedTool(correct)
    #
    #     #self.assertAlmostEqual(len(tested_tool), len(correct_tool), delta=3)
    #
    #     """
    #     for test, correct in zip(tested_tool, correct_tool):
    #         self.assertEqual(test, correct)
    #
    #
    #     """
    #
    #     #cleanup
    #     os.remove(os.getcwd() + "/allup_peak_results_classic.bed")

    
    # def test_gtf_allup(self):
    #
    #     """
    #       OUTPUT NO PEAKS,
    #     Performs basic all up test on entire program (using classic algorithm) using gtf file
    #
    #     """
    #
    #     #self.assertTrue(False, "test is currently disabled output from logging causes it to crash")
    #     args = ["-b", clipper.test_file("allup_test.bam"),
    #              "--gtfFile", clipper.test_file("ensembl_test.gtf"),
    #              "-g", "ENSG00000198901",
    #              "--outfile=" + os.getcwd() + "/allup_peak_results_ensembl_test.bed",
    #              "--debug",
    #              "--algorithm=classic"
    #             ]
    #
    #     (options, args) = self.parser.parse_args(args)
    #     options = override_options(options)
    #
    #
    #     main(options)
        
    # def test_filter(self):
    #
    #     """
    #
    #     allup test for transcriptome filter
    #     makes sure special test file
    #     detects only one peak when filter is enabled and detects two peaks when filter is disabled
    #
    #     """
    #     args = ["--bam", clipper.test_file("transcriptome_filter.sort.bam"),
    #             "--species", "hg19",
    #             "--gene", "ENSG00000198901",
    #             "--outfile", os.getcwd() + "/cut_off_included.bed",
    #             "-q",
    #             "--debug",
    #             ]
    #
    #     (options, args) = self.parser.parse_args(args)
    #     options = override_options(options)
    #     print(options.reverse_strand)
    #     main(options)
    #
    #     tested = open(os.getcwd() + "/cut_off_included.bed")
    #
    #     #problem with tracks being different
    #     tested_tool = pybedtools.BedTool(tested)
    #
    #
    #     #checks to make sure files are equal and there are not exact dups
    #     #cutoff stuff is borken and possibly buggy, need to fix later
    #     #self.assertEqual(len(tested_tool), 1)
    #
    #     #cleanup
    #     os.remove(os.getcwd() + "/cut_off_included.bed")
    #
    # def test_cutoff(self):
    #
    #     """
    #
    #     test_cutoff Tests that the cutoff code works if its enabled
    #
    #     """
    #
    #     args = ["-b", clipper.test_file("transcriptome_filter.sort.bam"),
    #              "-s", "hg19",
    #               "-g", "ENSG00000198901",
    #               '-g', "ENSG00000226167",
    #                "--outfile=" + os.getcwd() + "/no_cut_off.bed",
    #                "-q",
    #                "--disable_global_cutoff",
    #                '--debug'
    #             ]
    #     (options, args) = self.parser.parse_args(args)
    #     options = override_options(options)
    #     main(options)
    #
    #     tested = open(os.getcwd() + "/no_cut_off.bed")
    #
    #     #problem with tracks being different
    #     tested_tool = pybedtools.BedTool(tested)
    #
    #
    #     #checks to make sure files are equal and there are not exact dups
    #     #cutoff of stuff is broken and possibly buggy, need to fix later
    #     #self.assertEqual(len(tested_tool), 2)
    #
    #     #cleanup
    #     os.remove(os.getcwd() + "/no_cut_off.bed")
        

       

    
    
#    def test_build_bed(self):
#        
#        """
#        
#        Tests building transcript data from bed with both pre-mrna and mrna data
#        
#        """
#        
#        test = pybedtools.BedTool(clipper.test_file("test_bed_creation.bed"))
#
#        true_genes = {["chr1", "AS_STRUCTURE", "mRNA", 66999065, 67210057, ".", "+", ".", "transcript_id=ENST00000237247;effective_length=3997" ],
#                      ["chr1", "AS_STRUCTURE", "mRNA", 66999274, 66999274, ".", "+", ".", "transcript_id=ENST00000371039;effective_length=4080" ], 
#                      ["chr1", "AS_STRUCTURE", "mRNA", 66999297, 67145425, ".", "+", ".", "transcript_id=ENST00000424320;effective_length=951" ], 
# 
#        
#        genes = build_transcript_data_bed(test, False)
#        
#        self.assertDictEqual(true_genes, genes, "mrna genes not equal")
#        
#        #tests pre-mrna 
#        
#        true_lengths = {"ENST00000237247" : 210992, "ENST00000371039" : 211494, "ENST00000424320" : 146128}
#        genes = build_transcript_data_bed(test, True)
#        
#        self.assertDictEqual(true_genes, genes, "pre-mrna genes not equal")
#        self.assertDictEqual(true_lengths, lengths, "pre-mrna lengths not equal")
#        
#     def test_build_geneinfo(self):
#
#         """
#
#         Performs unit testing on build_geneinfo
#
#         I'm hopefully going to remove this method soon so no testing for now
#
#         """
#
#         #checks error mode
#         self.assertRaises(TypeError, build_geneinfo, None)
#
#         self.assertRaises(IOError, build_geneinfo, "foo")
#
#         #checks working mode
#         geneinfo = build_geneinfo(
#                     clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"))
#
#
#         true_values = {
#         "ENSG00000232113" : ["chr1",  "ENSG00000232113",  173604911,      173606273, "+"],
#         "ENSG00000228150" : ["chr1",  "ENSG00000228150",  10002980,       10010032,  "+"],
#         "ENSG00000223883" : ["chr1",  "ENSG00000223883",  69521580,       69650686,  "+"],
#         "ENSG00000135750" : ["chr1",  "ENSG00000135750",  233749749,      233808258, "+"],
#         "ENSG00000227280" : ["chr1",  "ENSG00000227280",  145373053,      145375554 ,"-"],
#         }
#
#
#         self.assertDictEqual(geneinfo, true_values)
#
#     def test_build_lengths(self):
#
#         """
#
#         Performs unit testing on build_lengths
#
#         I'm hopefully going to remove this method soon so no unit testing for now
#
#         """
#
#         #Checks error mode
#         self.assertRaises(ValueError, build_lengths, None)
#
#         self.assertRaises(ValueError, build_lengths, clipper.test_file("foo"))
#
#         #checks working mode
#         lengths = build_lengths(
#                     clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"))
#
#         true = {"ENSG00000232113" : 384,
#                 "ENSG00000228150" : 323,
#                 "ENSG00000223883" : 437,
#                 "ENSG00000135750" : 3141,
#                 "ENSG00000227280" : 212,
#                 }
#
#         self.assertDictEqual(lengths, true)
#
#     def test_add_species(self):
#
#         """
#
#         Performs unit testing on add_species
#
#         I'll probably refactor this a bit so I won't work to hard on testing this
#
#         """
#
#         #Case: object is returned as expected
#         result = add_species("hg19", [range(1, 22), "X", "Y"],
#                                         "foo",
#                                         "bar",
#                                         "baz")
#
#         self.assertEqual(result, {"chrs" : range(1, 22) + ["X"] + ["Y"],
#                           "gene_bed" : "foo",
#                           "mRNA" : "bar",
#                           "premRNA" : "baz"})
#
#     def test_get_acceptable_species(self):
#
#         """
#
#         Test get acceptable species
#
#         """
#
#         result = get_acceptable_species()
#
#         # make sure some main genomes are in here.
#         self.assertIn("hg19", result)
#         self.assertIn("mm9", result)
#         self.assertIn("mm10", result)
#         self.assertIn("GRCh38", result)
#
#
#     def test_build_transcript_data(self):
#         self.maxDiff = 10000000
#         """
#
#         Tests building transcript data and returning the proper values
#
#         Doesn't assume malformed data
#
#         """
#
#         #tests error modes
#         self.assertRaises(ValueError, build_transcript_data, None, None, None, None, True)
#
#         self.assertRaises(ValueError, build_transcript_data, "foo", "bar", "bar", "bar", True)
#
#         self.assertRaises(ValueError, build_transcript_data, "bar", None, None, None, True)
#
#         #tests hg19 to make sure its equal to logic
#         genes = build_transcript_data("test", None, None, None, True).sort()
#         true_genes = pybedtools.BedTool(
#                 [["chr1", "AS_STRUCTURE", "mRNA", 173604911, 173606273, ".", "+", ".", "gene_id=ENSG00000232113; effective_length=1147" ],
#                 ["chr1", "AS_STRUCTURE", "mRNA", 10002980, 10010032, ".", "+", ".", "gene_id=ENSG00000228150; effective_length=3088" ],
#                 ["chr1", "AS_STRUCTURE", "mRNA", 69521580, 69650686, ".", "+", ".", "gene_id=ENSG00000223883; effective_length=46051" ],
#                 ["chr1", "AS_STRUCTURE", "mRNA", 233749749, 233808258, ".", "+", ".", "gene_id=ENSG00000135750; effective_length=35997" ],
#                 ["chr1", "AS_STRUCTURE", "mRNA", 145373053, 145375554, ".", "-", ".", "gene_id=ENSG00000227280; effective_length=609" ]],
#                                         ).sort()
#
#         self.assertEqual(str(genes), str(true_genes))
#
#         #tests hg19 on premrna lengths
#         genes = build_transcript_data("test", None, None, None, False).sort()
#
#         true_genes = pybedtools.BedTool(
#                 [["chr1", "AS_STRUCTURE", "mRNA", 173604911, 173606273, ".", "+", ".", "gene_id=ENSG00000232113; effective_length=384" ],
#                 ["chr1", "AS_STRUCTURE", "mRNA", 10002980, 10010032, ".", "+", ".", "gene_id=ENSG00000228150; effective_length=323" ],
#                 ["chr1", "AS_STRUCTURE", "mRNA", 69521580, 69650686, ".", "+", ".", "gene_id=ENSG00000223883; effective_length=437" ],
#                 ["chr1", "AS_STRUCTURE", "mRNA", 233749749, 233808258, ".", "+", ".", "gene_id=ENSG00000135750; effective_length=3141" ],
#                 ["chr1", "AS_STRUCTURE", "mRNA", 145373053, 145375554, ".", "-", ".", "gene_id=ENSG00000227280; effective_length=212" ]],
#                                         ).sort()
#
#         self.assertEqual(str(genes), str(true_genes))
#
#         #Test custom files
#         #this should all work,
#         self.assertRaises(IOError, build_transcript_data, None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), None, clipper.data_file("test.AS.STRUCTURE_premRNA.lengths"), False)
#         build_transcript_data(None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"), clipper.data_file("test.AS.STRUCTURE_premRNA.lengths"), True)
#         build_transcript_data(None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"), clipper.data_file("test.AS.STRUCTURE_premRNA.lengths"), False)
#         build_transcript_data(None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"), None, False)
#         build_transcript_data(None, clipper.data_file("test.AS.STRUCTURE_genes.BED.gz"), None, clipper.data_file("test.AS.STRUCTURE_mRNA.lengths"), True)
    
    # def test_transcriptome_filter(self):
    #
    #     """
    #     DEPRECATED, no function named `transcriptome_filter`
    #     Tests transcriptome filter
    #     not great tests, but good enough to make sure we don't have regressions
    #
    #     """
    #
    #
    #     cluster = Peak(0,0,0,0,0,0,0,0,0,5,0,10,0,0,0,0)
    #     #cluster = {'Nreads' : 5, "size" : 10}
    #     transcriptome_size = 1000
    #     transcriptome_reads = 10000
    #     poisson_cutoff = .05
    #
    #     result = transcriptome_filter(poisson_cutoff, transcriptome_size, transcriptome_reads, cluster)
    #
    #     self.assertEqual(result, 1)
    #
    #     #cluster = {'Nreads' : 10000, "size" : 100}
    #     cluster = Peak(0,0,0,0,0,0,0,0,0,10000,0,100,0,0,0,0)
    #
    #     transcriptome_size = 1000
    #     transcriptome_reads = 10000
    #     poisson_cutoff = .05
    #
    #     result = transcriptome_filter(poisson_cutoff, transcriptome_size, transcriptome_reads, cluster)
    #     self.assertEqual(result,0.0)
    #
    #     #cluster = {'Nreads' : 0, "size" : 0}
    #     cluster = Peak(0,0,0,0,0,0,0,0,0,0,0,0,0)
    #     transcriptome_size = 0
    #     transcriptome_reads = 10000
    #     poisson_cutoff = .05
    #
    #     result = transcriptome_filter(poisson_cutoff, transcriptome_size, transcriptome_reads, cluster)
    #     self.assertEqual(result, 1)
        

        

    #
    # def tearDown(self):
    #     pass

#tests for hadoop mapping, currently not used / not worth the time to fix
#    def test_mapper_premrna(self):
#        
#        """
#        
#        tests the mapper to make sure that its not breaking / outputs call_peaks results
#        
#        """
#        
#        args = ["-b", pkg_resources.resource_filename(__name__, "../test/allup_test.bam"),
#                 "-s", "hg19",
#                 "-g", "ENSG00000198901", 
#                 "--outfile=" + os.getcwd() + "/allup_peak_results.bed",
#                ]
#
#        (options, args) = self.parser.parse_args(args)
#        
#        mapper(options, "chr1    66999065    67210057    ENST00000237247    0    +    67000041    67208778    0    27    25,123,64,25,84,57,55,176,12,12,25,52,86,93,75,501,81,128,127,60,112,156,133,203,65,165,1302,    0,863,92464,99687,100697,106394,109427,110161,127130,134147,137612,138561,139898,143621,146295,148486,150724,155765,156807,162051,185911,195881,200365,205952,207275,207889,209690,")
#   
#    def test_mapper_mrna(self):
#        
#        """
#        
#        Tests the pre mrna mapper
#        
#        """
#        
#        args = ["-b", pkg_resources.resource_filename(__name__, "../test/allup_test.bam"),
#                 "-s", "hg19",
#                 "-g", "ENSG00000198901", 
#                 "--outfile=" + os.getcwd() + "/allup_peak_results.bed",
#                ]
#
#        (options, args) = self.parser.parse_args(args)
#        
#        mapper(options, "chr1    66999065    67210057    ENST00000237247    0    +    67000041    67208778    0    27    25,123,64,25,84,57,55,176,12,12,25,52,86,93,75,501,81,128,127,60,112,156,133,203,65,165,1302,    0,863,92464,99687,100697,106394,109427,110161,127130,134147,137612,138561,139898,143621,146295,148486,150724,155765,156807,162051,185911,195881,200365,205952,207275,207889,209690,")

# if __name__ == '__main__':
#     unittest.main()
#     os.remove(pkg_resources.resource_filename(__name__, "../src/peak_results.BED"))
