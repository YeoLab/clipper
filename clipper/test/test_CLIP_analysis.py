'''
Created on Sep 10, 2012

@author: gabrielp
'''
import unittest
import CLIP_analysis 
from CLIP_analysis import *
from optparse import OptionParser, SUPPRESS_HELP 
import os
import pkg_resources       
import clipper
import pybedtools
from numpy.testing import *
from numpy import array
class Test(unittest.TestCase):
    
    """  Runs tests on CLIP_analyis module """
    
    parser = None
    
    def setUp(self):
        self.parser = OptionParser()
    
        self.parser.add_option("--clusters", dest="clusters", help="BED file of clusters", metavar="BED")
        self.parser.add_option("--bam", dest="bam")
        self.parser.add_option("--species", "-s", dest="species", help = "genome version")
        ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
        self.parser.add_option("--reAssign", dest="assign", action="store_true", default=False, help="re-assign clusters, if not set it will re-use existing assigned clusters") 
        ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
        self.parser.add_option("--rePhast", dest="rePhast", action="store_true", default=False, help="re-calculate conservation, must have been done before") 
        self.parser.add_option("--old_motifs", dest="reMotif", action="store_false", default=True, help="use old motif files")
        self.parser.add_option("--motif", dest="motif", action="append", help="Files of motif locations", default=None)
        self.parser.add_option("--homer", dest="homer", action="store_true", help="What does this do?", default=False)
        self.parser.add_option("--conservation", dest="cons", help="what does this do?", action="store_true")
        self.parser.add_option("--structure", dest="structure", help="what does this do?", action="store_true")
        self.parser.add_option("--nrand", dest="nrand", default=3, help="selects number of times to randomly sample genome", type="int")
        self.parser.add_option("--k", dest="k", action="append", help="what does this do?", default=[6])
        self.parser.add_option("--outdir", "-o", dest="outdir", default=os.getcwd(), help="directory for output, default:cwd")
        self.parser.add_option("--run_phast", dest="run_phast", action="store_true", help="what does this do?", default=False)
        self.parser.add_option("--AS_Structure", dest="as_structure", help="Location of AS_Structure directory (chromosme files should be inside)", default=None)
        self.parser.add_option("--genome_location", dest="genome_location", help="location of all.fa file for genome of interest", default=None)
        self.parser.add_option("--homer_path", dest="homer_path", help="path to homer, if not in default path", default=None)
        self.parser.add_option("--phastcons_location", dest="phastcons_location", help="location of phastcons file", default=None)
        self.parser.add_option("--regions_location", dest="regions_location" , help="directory of genomic regions for a species", default=None)
        self.parser.add_option("--motif_location", dest="motif_location", help="directory of pre-computed motifs for analysis", default=os.getcwd())
        self.parser.add_option("--runPhast", dest="runPhast", action="store_true", default=False, help="Run Phastcons ") 


     
    def testName(self):
        pass

    def test_allup(self):
        
        """ runs entire program on small test dataset """

        args = ["--clusters", clipper.test_file("clip_analysis_test_peak_results.bed"),
                "-s", "mm9",
                "--bam", clipper.test_file("allup_test.bam"),
                "--AS_Structure", "/home/gabrielp/bioinformatics/Yeo_Lab/clip_analysis_metadata/mm9data4",
                '--genome_location', '/home/gabrielp/bioinformatics/Yeo_Lab/clip_analysis_metadata/mm9/mm9.fa', 
                #'--regions_location', clipper.test_file("knownGene_sample.gtf"),
                "--regions_location", '/home/gabrielp/bioinformatics/Yeo_Lab/clip_analysis_metadata/regions',
                '--phastcons_location', '/home/gabrielp/bioinformatics/Yeo_Lab/clip_analysis_metadata/phastcons/mm9_phastcons.bw',
                '--motif', 'AAAAAA',
                '--nrand', '1',
                '--rePhast',
                '--runPhast'
                ]    
        (options, args) = self.parser.parse_args(args)
        CLIP_analysis.main(options)
    
   
    def test_CLIP_figure(self):
        pass
    
    def test_intersection(self):
        pass
    
    def test_adjust_offsets(self):
        
        """
        
        Tests adjust offsets, only on bed12 files
        
        """
        
        offsets = {"ENSMUSG00000051951_1_147" : 10, 
                   "ENSG00000198901_2_52" : 10 ,
                   "ENSG00000198901_3_239" : 10, 
                   "ENSG00000198901_4_85" : 10 ,
                   "ENSG00000198901_5_47" : 10 ,
                   "ENSG00000198901_6_119" : 10 ,
                   "ENSG00000198901_7_58" : 10 ,
                   "ENSG00000198901_8_588" : 10 ,
                   "ENSG00000198901_10_92" : 10 ,
                   "ENSG00000198901_11_59" : 10 ,
                   "ENSG00000198901_12_196" : 10 ,
                   "ENSG00000198901_13_36" : 10 ,

                   }
        bedtool = pybedtools.BedTool(clipper.test_file("clip_analysis_test_peak_results.bed"))
        
        
        results = adjust_offsets(bedtool, offsets)
        
        true_results = ((3206126, 3206130),
                   (91513660, 91513664),
                   (91517394, 91517398),
                   (91517935, 91517939),
                   (91522404, 91522408),
                   (91523607, 91523611),
                   (91524250, 91524254),
                   (91525137, 91525141),
                   (91527347, 91527351),
                   (91527937, 91527941),
                   (91528034, 91528038),
                   (91537658, 91537662),
                   )
        for result, true_result in zip(results, true_results):
            self.assertEqual(int(result[6]), true_result[0])
            self.assertEqual(int(result[7]), true_result[1])
        
    def test_adjust_offsets_short(self):
        
        """
        
        Tests adjust offsets on bed files shorter than bed8 bed files 
        
        """
        tool = pybedtools.BedTool("chr15    91512755    91512836    ENSG00000198901_1_147    0    -", from_string=True)
        offsets = {"ENSG00000198901_1_147" : 10}
        results = adjust_offsets(tool, offsets)
        
    def test_count_genomic_types(self):
        
        """
        
        Tests count genomic types
        
        """
        
        result, bed_result = parse_AS_STRUCTURE_dict("test", clipper.test_dir())
        result = count_genomic_types(result)
        
        self.assertDictEqual(result, {"CE:" : 14})
        
    def test_parse_AS_STRUCTURE_dict(self):
        
        """
        
        Tests build AS structure dict
        
        """
        
        true_tool = pybedtools.BedTool("chr10\t127496045\t127555714\tENSMUSG00000040054\t0\t+\n", from_string=True)
        
        result, result_bed = parse_AS_STRUCTURE_dict("test", clipper.test_dir())
        print str(result_bed)
        self.assertEqual(len(true_tool.intersect(result_bed)), 1)
        test_result = result["ENSMUSG00000040054"]
        
        true_exons = {0:'127496045-127496082', 
                    1:'127528690-127528832',
                    2:'127533494-127533579', 
                    3:'127545949-127546087', 
                    4:'127547810-127548404', 
                    5:'127549637-127549823', 
                    6:'127550518-127550737', 
                    7:'127551389-127551839', 
                    8:'127552080-127552141', 
                    9:'127553116-127553225', 
                    10:'127553361-127553463', 
                    11:'127553602-127553813',
                    12:'127555610-127555714'}
        
  
        self.assertDictEqual(test_result["exons"], true_exons)
        
        self.assertDictEqual(test_result["introns"], {0 :'127496083-127528689', 
                                                     1 :'127528833-127533493', 
                                                     2 :'127533580-127545948', 
                                                     3 :'127546088-127547809', 
                                                     4 : '127548405-127549636',
                                                     5 :'127549824-127550517', 
                                                     6 :'127550738-127551388', 
                                                     7 :'127551840-127552079', 
                                                     8 : '127552142-127553115', 
                                                     9 : '127553226-127553360', 
                                                     10 : '127553464-127553601', 
                                                     11 :'127553814-127555609'}, 
                             "introns not equal")

        self.assertDictEqual(test_result["types"], {0 : "CE:", 
                                                    1 : "CE:", 
                                                    2 : "CE:", 
                                                    3 : "CE:",
                                                    4 : "CE:", 
                                                    5 : "CE:", 
                                                    6 : "CE:", 
                                                    7 : "CE:", 
                                                    8 : "CE:", 
                                                    9 : "CE:", 
                                                    10 : "CE:", 
                                                    11 : "CE:", 
                                                    12 : "CE:" }, 
                             "types not equal")

        self.assertEqual(test_result["tx_stop"], 127555714)
        self.assertEqual(test_result["tx_start"], 127496045)
        self.assertEqual(test_result["premRNA_length"], 59670)
        self.assertEqual(test_result["mRNA_length"], 2451)

    
    def test_count_genomic_region_sizes(self):
        
        """
        
        Regression test to make sure we can count all genomic regions
        Need to write better test, for now we'll just make sure it doesn't crash
        
        """
        
        results = count_genomic_region_sizes('/home/gabrielp/bioinformatics/Yeo_Lab/clip_analysis_metadata/regions', "mm9")
        
    def test_assign_to_regions(self):
        
        """
        
        Basic test of assign to regions 
        
        """
        
        #bedtool = pybedtools.BedTool(clipper.test_file("clip_analysis_test_peak_results.bed"))
        
        #assign_to_regions(tool, clusters, speciesFA, regions_dir, regions, species="hg19", nrand = 3, getseq=False):

    
    def test_build_assigned_from_existing(self):
        pass
    
    def test_eliminate_invalid_bed12(self):
        pass
    
    def test_get_offsets_bed12(self):
        
        """
        
        Very lazy test of bed12 offsets
        
        """
        
        bedtool = pybedtools.BedTool(clipper.test_file("clip_analysis_test_peak_results.bed"))
        result = get_offsets_bed12(bedtool)
    
        self.assertDictContainsSubset({"ENSMUSG00000051951_1_147" : 5}, result)
            
    def plot_motif_dist(self):
        pass
    
    def test_RNA_position(self):
        pass
    
    def test_run_kmerdiff(self):
        pass
    
    def test_run_homer(self):
        
        """
        
        tests run_homer, makes sure homer is installed and outputs the correct files
        
        """
        
        #foreground = clipper.test_file("clip_analysis_test_peak_results.bed.all.real.fa")
        #background = clipper.test_file("clip_analysis_test_peak_results.bed.all.random.fa")
        #run_homer(foreground, background)
        
    def test_write_seqs(self):
        pass
    
    def motif_boxplots(self):
        pass
    
    def bedlengths(self):
        pass
    
    def chop(self):
        pass
    
    def test_fa_file(self):
        pass
    
    def make_dir(self):
        pass
    
    def test_count_total_reads(self):
        
        """
        
        Makes sure we are counting the correct number of reads in genes
        
        """
        
        bam = pybedtools.BedTool(clipper.test_file("allup_test.bam"))
        gene_dfn = pybedtools.BedTool(clipper.test_file("hg19_genes.bed"))
        
        result = count_total_reads(bam, gene_dfn)
        
        self.assertEqual(result, 2086)
    
    def test_count_reads_per_cluster(self):
        
        """
        
        Tests count reads per cluster
        
        """
        
        bedtool = pybedtools.BedTool(clipper.test_file("clip_analysis_test_peak_results.bed"))
        
        total_reads, reads_per_cluster = count_reads_per_cluster(bedtool)
        
        self.assertListEqual([147,52, 239, 85, 47, 119, 58, 588, 92, 59, 196, 36], reads_per_cluster)
        self.assertEqual(sum([147,52, 239, 85, 47, 119, 58, 588, 92, 59, 196, 36]), total_reads)
    
    def test_count_reads_per_cluster_merged(self):
        
        """
        
        Tries to count a read in a cluster that has been merged
        
        """
        
        tool = pybedtools.BedTool("chr15    91512755    91512836    ENSMUSG00000025736_1_83;ENSMUSG00000091321_6_83    0    -", from_string=True)
        total_reads, reads_per_cluster = count_reads_per_cluster(tool)
        
        self.assertListEqual([83], reads_per_cluster)
    def test_calculate_kmer_diff(self):
        
        """
        
        Tests calculate kmer diff
        
        Not testing computations (if kmer diff works) here, just that the dir structure is as it should be.  
        
        """
        
        result = calculate_kmer_diff([3,4], ['all'], "clip_analysis_test_peak_results.bed", "/home/gabrielp/bioinformatics/Yeo_Lab/clipper/clipper/test")
        self.assertListEqual(result.keys(), ["all"])
        self.assertListEqual(result['all'].keys(), [3,4])
    
    def test_get_mean_phastcons(self):
        
        """
        
        Tests get_mean phastcons
        
        """
        
        phastcons = '/home/gabrielp/bioinformatics/Yeo_Lab/clip_analysis_metadata/phastcons/mm9_phastcons.bw'
        tool = pybedtools.BedTool("chr15    91512755    91512836    ENSG00000198901_1_147    0    -", from_string=True)
        result = get_mean_phastcons(tool, phastcons)
        
        #pre-computed correct answer
        self.assertAlmostEqual(result[0], 0.01738272, delta=5)
        
    def test_calculate_phastcons(self):
        
        """
        
        Tests calculate phastcons
        
        Only tests that stucture of file is correct
        
        """
        cluster_regions = {"exon" : 
                           {'real' : pybedtools.BedTool(clipper.test_file("clip_analysis_test_peak_results.bed.exon.real.BED")),
                            'rand' : {0: pybedtools.BedTool(clipper.test_file("clip_analysis_test_peak_results.bed.exon.rand.0.BED"))}}
                           }
        regions = (["all", "exon", "UTR3", "UTR5", "proxintron500", "distintron500"])    
        phastcons = '/home/gabrielp/bioinformatics/Yeo_Lab/clip_analysis_metadata/phastcons/mm9_phastcons.bw'

        result = calculate_phastcons(regions, cluster_regions, phastcons)
        self.assertListEqual(result.keys(), ['real', 'rand'])
        self.assertListEqual(result['real'].keys(), ['all', 'exon'])
        
    def test_get_motif_distance(self):
        
        """
        
        Tests get motif distance
        
        """
        
        result = get_motif_distance(pybedtools.BedTool(clipper.test_file("motif_distance_test_cluster.bed")),
                           pybedtools.BedTool(clipper.test_file("motif_distance_test_motif.bed")))
        
        #hand calculated results from motif distance
        self.assertListEqual(result, [-12,112, 12])
        
    def test_calculate_motif_distance(self):
        cluster_regions = {"exon" : 
                           {'real' : pybedtools.BedTool(clipper.test_file("clip_analysis_test_peak_results.bed.exon.real.BED")),
                            'rand' : {0: pybedtools.BedTool(clipper.test_file("clip_analysis_test_peak_results.bed.exon.rand.0.BED"))}}
                           }
        
        region_sizes = {'all' : 30, "exon" : 30}
        #doesn't matter if I output actual results
        result = calculate_motif_distance(cluster_regions, region_sizes,
                                 pybedtools.BedTool(clipper.test_file("motif_distance_test_motif.bed")))
        
        #make sure dict structure is sound
        self.assertListEqual(result.keys(), ['all', 'exon'])
        self.assertListEqual(result['all'].keys(), ['real', 'rand'])
        self.assertListEqual(result['all']['real'].keys(), ['dist', 'size'])
        
    def test_build_genomic_regions(self):
        """
        
        Tests build genomic regions
        
        """

        CDS = pybedtools.BedTool("""chr1\t7700\t7900\tfoo\t0\t+\n
                                   chr1\t7999\t8500\tfoo\t0\t+\n""", from_string = True)
        UTR5 = pybedtools.BedTool("""chr1\t7499\t7700\tfoo\t0\t+\n""", from_string = True)
        UTR3 = pybedtools.BedTool("""chr1\t8500\t9000\tfoo\t0\t+\n""", from_string = True)
        proxintron = pybedtools.BedTool("""chr1\t100\t300\tfoo\t0\t+\n
                                          chr1\t798\t998\tfoo\t0\t+\n
                                          chr1\t2000\t2200\tfoo\t0\t+\n
                                          chr1\t2798\t2998\tfoo\t0\t+\n
                                          chr1\t6000\t6200\tfoo\t0\t+\n
                                          chr1\t6798\t6998\tfoo\t0\t+\n
                                          chr1\t7900\t7998\tfoo\t0\t+\n""", from_string = True
                                          )
        distintron = pybedtools.BedTool("""chr1\t301\t797\tfoo\t0\t+\n
                                           chr1\t2201\t2797\tfoo\t0\t+\n
                                           chr1\t6201\t6797\tfoo\t0\t+\n""", from_string = True)
        
        regions = build_genomic_regions(pybedtools.BedTool(clipper.test_file("test.gtf")), prox_distance=200)        
        
        #print UTR3

        #print regions['UTR3']
        print proxintron
        print regions['proxintron']
        #print regions['distintron']
        
        self.assertEqual(len(CDS.intersect(regions['CDS'], f= 1.0, r = True)), 2)
        self.assertEqual(len(UTR5.intersect(regions['UTR5'], f= 1.0, r = True)), 1)
        self.assertEqual(len(UTR3.intersect(regions['UTR3'], f= 1.0, r = True)), 1)
        self.assertEqual(len(proxintron.intersect(regions['proxintron'], f= 1.0, r = True)), 7)
        self.assertEqual(len(distintron.intersect(regions['distintron'], f= 1.0, r = True)), 3)
        

        
    def test_main(self):
        pass
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()