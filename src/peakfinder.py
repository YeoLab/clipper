#!/nas3/yeolab/Software/Python_dependencies/bin/python
#We will follow the UCSC genome browser assumption of using a zero based half open cord system
from pysam import rmdup
from optparse import OptionParser, SUPPRESS_HELP
import os
import sys
from subprocess import call
import pickle
import time
import pybedtools
import gzip
import pkg_resources
import pp
import math
from call_peak import call_peaks, peaks_from_info, get_FDR_cutoff_mean, poissonP
import logging
#logging.basicConfig(level=logging.INFO)
logging.disable(logging.INFO)

#define verbose printing here for test cases
global varboseprint
def verboseprint(*args):
        # Print each argument separately so caller doesn't need to
        # stuff everything to be printed into a single string
            for arg in args:
                print arg,
            print



def trim_reads(bamfile):
    
    """

    Wrapper to remove PCR duplicate reads from bed file
    
    Input
    bamfile -- location of bamfile on disk
    assumes .bam ending of bam file
    returns bamfile_trimed.bam file
    
    """
    
    if not os.path.exists(bamfile):
        raise NameError("file %s does not exist" % (bamfile))
    
    outfile = ".".join(bamfile.split(".")[:-1])
    outfile += "_trimmed.bam"
    rmdup("-S", bamfile, outfile)
    return outfile

def check_for_index(bamfile, make=True):
    
    """

    Checks to make sure a BAM file has an index, if the index does not exist it is created
    
    Usage undefined if file does not exist (check is made earlier in program)
    bamfile - a path to a bam file
    
    Returns 1 
    TODO make it so a failaure returns 0

    """

    if not os.path.exists(bamfile):
        raise NameError("file %s does not exist" % (bamfile))
    
    if os.path.exists(bamfile + ".bai"):
        return 1
    else:
        verboseprint("Index for %s does not exist, indexing bamfile" % (bamfile))
        process = call(["samtools", "index", str(bamfile)])
        
        if process == -11: 
            raise NameError("file %s not of correct type" % (bamfile))
        
        return 1

def build_geneinfo(BED):
    
    """

    Loads bed file into a dictionary with the key being the name and a string being the value
    
    Input:
    BED -- a bed file to load
    
    Return:
    A dictionary with the key being the name position of the bed file and the values being the
    ordered bed file
    
    TODO: Refactor to used bedtools instead
    
    """
    
    #opens bed file, either zipped or unzipped
    try:
        bedfile = gzip.open(BED, "rb")
    except:
        bedfile = open(BED, "r")
        
    GI = dict()
    
    for line in bedfile.readlines():
        chromosome, start, stop, name, score, signstrand = line.strip().split("\t")
        chromosome.replace("chr", "")
        GI[name] = "|".join([chromosome, name, start, stop, str(signstrand)])
    bedfile.close()
    return GI

def build_lengths(f):
    
    """
    
    Builds a dictionary of gene names and lengths of mappable regions in that gene
    
    Input:
    A two column file with the first column being the gene name and the second column being the
    mappable length of the gene
    
    Return:
    A dictionary with the key being the name of the gene and the value being the length
    
    """
    
    FI = open(f, "r")
    gene_lengths = {}

    for line in FI.readlines():
        name, gene_length = line.strip().split("\t")
        gene_lengths[name] = int(gene_length)

    FI.close()

    return gene_lengths



def add_species(species, chrs, bed, mrna, premrna):
    
    """

    Creates a dictionary containing all information needed to perform peak calling calcluations 
    for a single species
    
    Paramaters
    -----------
    species: string currently not used
    chrs: list specifying all the chromosomes in a given species
    bed: path to a bed file that contains information on genes (custom file *STRUCTURE_genes.BED.gz)
    mrna: path to a file that contains mRNA lengths (custom CSV file contains gene names follwed by gene lengths)
    premrna: path to a file that contains pre-mRNA lengths (custom CSV file contains gene names follwed by gene lengths_
    
    Returns dict of all items passed to it
    
    TODO:  Add checking to verify that file are actually passed
    """
    par = dict()
    
    #this is non-pythonic, should just combine all lists
    par["chrs"] = [item for sublist in chrs for item in sublist] #expand sublists
    par["gene_bed"] = bed
    par["mRNA"] = mrna
    par["premRNA"] = premrna
    return par
        
def main(options):
    
    job_server = pp.Server(ncpus=options.np)
    bamfile = options.bam
    if os.path.exists(bamfile):
        bamfile = os.path.abspath(bamfile) #re-set to include the full path to bamfile
        verboseprint("bam file is set to %s\n" % (bamfile))
    else:
        sys.stderr.write("Bam file not defined")
        raise IOError

    species = options.species
    geneBed = options.geneBEDfile
    genemRNA = options.geneMRNAfile
    genePREmRNA = options.genePREMRNAfile

    species_parameters = dict()

    if species is None and geneBed is None:
        print "You must set either \"species\" or \"geneBed\"+\"geneMRNA\"+\"genePREMRNA\""
        exit()

    lenfile = ""
    
    species_parameters["hg19"] = add_species("hg19", [range(1, 22), "X", "Y"],
                                             pkg_resources.resource_filename(__name__, "../data/hg19.AS.STRUCTURE_genes.BED.gz"),
                                             pkg_resources.resource_filename(__name__, "../data/hg19.AS.STRUCTURE_mRNA.lengths"),
                                             pkg_resources.resource_filename(__name__, "../data/hg19.AS.STRUCTURE_premRNA.lengths"))
    species_parameters["hg18"] = add_species("hg18", [range(1, 22), "X", "Y"],
                                             pkg_resources.resource_filename(__name__, "../data/hg18.AS.STRUCTURE_genes.BED.gz"),
                                             pkg_resources.resource_filename(__name__, "../data/hg18.AS.STRUCTURE_mRNA.lengths"),
                                             pkg_resources.resource_filename(__name__, "../data/hg18.AS.STRUCTURE_premRNA.lengths"))
    species_parameters["mm9"] = add_species("mm9", [range(1, 19), "X", "Y"],
                                             pkg_resources.resource_filename(__name__, "../data/mm9.AS.STRUCTURE_genes.BED.gz"),
                                             pkg_resources.resource_filename(__name__, "../data/mm9.AS.STRUCTURE_mRNA.lengths"),
                                             pkg_resources.resource_filename(__name__, "../data/mm9.AS.STRUCTURE_premRNA.lengths"))
    acceptable_species = ",".join(species_parameters.keys())
    
    #error checking
    if species is not None and geneBed is not None:
        print "You shouldn't set both geneBed and species, defaults exist for %s" % (acceptable_species)
        exit()
    if species is not None and species not in species_parameters:
        print "Defaults don't exist for your species: %s. Please choose from: %s or supply \"geneBed\"+\"geneMRNA\"+\"genePREMRNA\"" % (species, acceptable_species)
        exit()
    if species is None:
        species = "custom"
        species_parameters["custom"] = add_species("custom", [range(1, 22), "X", "Y"], geneBed, genemRNA, genePREmRNA)

    #error checking done, this does... something.  This is more setup phase  Uses pre-mrnas?
    if options.premRNA is True:
        lenfile = species_parameters[species]["premRNA"]
    else:
        lenfile = species_parameters[species]["mRNA"]
    
    #builds dict to do processing on, 
    lengths = build_lengths(lenfile)
    genes = build_geneinfo(species_parameters[species]["gene_bed"])
    margin = int(options.margin)
    
    #this should be fixed, args should initally be ints if passed
    if options.maxgenes is not None:
        maxgenes = int(options.maxgenes)

    minreads = int(options.minreads)
    poisson_cutoff = options.poisson_cutoff

    g = options.gene 
    
    gene_list = list()
    
    #gets all the genes to call peaks on
    try:
        if len(g) > 0:
            gene_list = g
    except:
        gene_list = genes.keys()
                
    
    results = []
    transcriptome_size = 0
    
    #I think this calls peaks for each gene in the gene list, which could be every gene in the genome
    running_list = []
    length_list = []
    
    for n, gene in enumerate(gene_list):
        #again, hacky should be factored to a single if statement, need to be more explicit about code paths
        if options.maxgenes == None:
            pass
        else:
            if n >= maxgenes:
                break
            
        geneinfo = genes[gene]
        
        #There is a better way of doing timing.  
        t = time.strftime('%X %x %Z')
        verboseprint(geneinfo + " started:" + str(t))
        transcriptome_size += lengths[gene]
        #TODO make it so transcript size isn't always used
        #this is a filter operation, should make it as such
        running_list.append(genes[gene])
        length_list.append(lengths[gene])
        
        verboseprint(lengths[gene])
        
        #for debugging purposes, sometimes 
        #call_peaks(genes[gene], lengths[gene], None, bamfile,  margin, options.FDR_alpha, options.threshold, 
        #                       minreads,  poisson_cutoff,  options.plotit, 10, 1000, options.SloP, False,) 

 
    combined_list = zip(running_list, length_list)
  
    jobs = [job_server.submit(call_peaks,
                              args=(gene, length, None, bamfile, margin, options.FDR_alpha, options.threshold,
                               minreads, poisson_cutoff, options.plotit, 10, 1000, options.SloP, False,),
                              depfuncs=(peaks_from_info, get_FDR_cutoff_mean,
                                          verboseprint,),
                              modules=("pysam", "os", "sys", "scipy", "math", "time",
                               "random", "peaks"),) for gene, length in combined_list]

    for job in jobs:
        results.append(job())   
    verboseprint("finished with calling peaks")

    #if we are going to save and output as a pickle file we should output as a pickle file
    #we should factor instead create a method or object to handle all file output
    if options.save_pickle is True:
        pickle_file = open(options.outfile + ".pickle", 'w')
        pickle.dump(results, file=pickle_file)                
    
    #combine results
    allpeaks = set([])

    #count total number of reads in transcriptiome
    transcriptome_reads = 0
    
    for gene_result in results:
        if gene_result is not None:
            verboseprint("nreads", gene_result['nreads'])
            transcriptome_reads += gene_result['nreads']
    print "Transcriptome size is %d, transcriptome reads are %d" % (transcriptome_size, transcriptome_reads)
    
    #is this a missed indent?
    for gener in results:
        if gener['clusters'] is None:
            print >> sys.stderr, gener, "no clusters"
            continue
        
        for cluster in gener['clusters'].keys():
            try:
                transcriptomeP = poissonP(transcriptome_reads, gener['clusters'][cluster]['Nreads'], transcriptome_size, gener['clusters'][cluster]['size'])
                if math.isnan(transcriptomeP):
                    print "Transcriptome P is NaN, transcriptome_reads = %d, cluster reads = %d, transcriptome_size = %d, cluster_size = %d" % (transcriptome_reads, gener['clusters'][cluster]['Nreads'], transcriptome_size, gener['clusters'][cluster]['size'])
                    continue
                
                if transcriptomeP > poisson_cutoff:
                    print "%s\n Failed Transcriptome cutoff with %s reads, pval: %s" % (cluster, gener['clusters'][cluster]['Nreads'], transcriptomeP)
                    continue
                
                min_pval = 1

                corrected_SloP_pval = gener['clusters'][cluster]['SloP']
                corrected_Gene_pval = gener['clusters'][cluster]['GeneP']

                if corrected_SloP_pval < poisson_cutoff or corrected_Gene_pval < poisson_cutoff:
                    min_pval = min([corrected_SloP_pval, corrected_Gene_pval])
                else:
                    verboseprint("Failed Gene Pvalue: %s and failed SloP Pvalue: %s for cluster %s" % (corrected_Gene_pval, corrected_SloP_pval, cluster))
                    continue


                (chrom, g_start, g_stop, peak_name, geneP, signstrand, thick_start, thick_stop) = cluster.split("\t")
                #print >> sys.stderr, cluster                           
                bedline = "%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d" % (chrom, int(g_start), int(g_stop), peak_name, min_pval, signstrand, int(thick_start), int(thick_stop))
                allpeaks.add(bedline)

            except NameError as e:
                print >> sys.stderr, e
                print >> sys.stderr, "parsing failed"
                raise e
        
    #again redundant code 
    outbed = options.outfile + ".BED"
    color = options.color
    pybedtools.BedTool("\n".join(allpeaks), from_string=True).sort(stream=True).saveas(outbed, trackline="track name=\"%s\" visibility=2 colorByStrand=\"%s %s\"" % (outbed, color, color))
    print "wrote peaks to %s" % (options.outfile)
    "\n".join(allpeaks)
    return 1
 

def call_main():
    
    usage = "\npython peakfinder.py -b <bamfile> -s <hg18/hg19/mm9>\n OR \npython peakfinder.py -b <bamfile> --customBED <BEDfile> --customMRNA <mRNA lengths> --customPREMRNA <premRNA lengths>"
    description = "Fitted Accurate Peaks. Michael Lovci 2012. CLIP peakfinder that uses fitted smoothing splines to define clusters of binding.  Computation is performed in parallel using MPI.  You may decide to use the pre-built gene lists by setting the --species parameter or alternatively you can define your own list of genes to test in BED6/12 format and provide files containing the length of each gene in PREMRNA form and MRNA form (both are required). Questions should be directed to michaeltlovci@gmail.com."
    parser = OptionParser(usage=usage, description=description)

    parser.add_option("--bam", "-b", dest="bam", help="A bam file to call peaks on", type="string", metavar="FILE.bam")

    parser.add_option("--species", "-s", dest="species", help="A species for your peak-finding, either hg19 or mm9")
    
    #we don't have custom scripts or documentation to support this right now, removing until those get added in
    #parser.add_option("--customBED", dest="geneBEDfile", help="bed file to call peaks on, must come withOUT species and with customMRNA and customPREMRNA", metavar="BEDFILE")
    #parser.add_option("--customMRNA", dest="geneMRNAfile", help="file with mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")
    #parser.add_option("--customPREMRNA", dest="genePREMRNAfile", help="file with pre-mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")
    parser.add_option("--outfile", "-o", dest="outfile", default="fitted_clusters", help="a bed file output, default:%default")
    parser.add_option("--gene", "-g", dest="gene", action="append", help="A specific gene you'd like try", metavar="GENENAME")
    parser.add_option("--minreads", dest="minreads", help="minimum reads required for a section to start the fitting process.  Default:%default", default=3, type="int", metavar="NREADS")
    parser.add_option("--margin", dest="margin", type="int", help="find sections of genes within M bases that have genes and perform fitting. Default:%default", default=15, metavar="NBASES")
    parser.add_option("--trim", dest="trim", action="store_true", default=False, help="Trim reads with the same start/stop to count as 1")
    parser.add_option("--premRNA", dest="premRNA", action="store_true", help="use premRNA length cutoff, default:%default", default=False)
    parser.add_option("--poisson-cutoff", dest="poisson_cutoff", type="float", help="p-value cutoff for poisson test, Default:%default", default=0.05, metavar="P")
    parser.add_option("--FDR", dest="FDR_alpha", type="float", default=0.05, help="FDR cutoff for significant height estimation, default=%default")
    parser.add_option("--threshold", dest="threshold", type="int", default=None, help="Skip FDR calculation and set a threshold yourself")
    parser.add_option("--maxgenes", dest="maxgenes", default=None, help="stop computation after this many genes, for testing", metavar="NGENES")
    parser.add_option("--processors", dest="np", default="autodetect", help="Number of processors to use. Default: All processors on machine", type="str", metavar="NP")
    parser.add_option("--superlocal", action="store_true", dest="SloP", default=False, help="Use super-local p-values, counting reads in a 1KB window around peaks")
    parser.add_option("--color", dest="color", default="0,0,0", help="R,G,B Color for BED track output, default:black (0,0,0)")
    parser.add_option("--plot", "-p", dest="plotit", action="store_true", help="make figures of the fits", default=False)
    parser.add_option("--verbose", "-q", dest="verbose", action="store_true", help="suppress notifications")
    parser.add_option("--save-pickle", dest="save_pickle", default=False, action="store_true", help="Save a pickle file containing the analysis")
    (options, args) = parser.parse_args()
    
    
    #creates verbose or scilent output mode
    global verboseprint
    if options.verbose:
        def verboseprint(*args):
        # Print each argument separately so caller doesn't need to
        # stuff everything to be printed into a single string
            for arg in args:
                print arg,
            print
    else:   
        verboseprint = lambda *a: None      # do-nothing function

    
    #enforces required usage    
    if not (options.bam and ((options.species) or (options.geneBEDfile and options.geneMRNAfile and options.genePREMRNAfile))):
        parser.print_help()
        exit()
    
    #If triming option is set use pysam to remove duplicate reads for us, trims strictly ignoring paired end and strandness
    if options.trim:
        options.bam = trim_reads(options.bam)
    
    check_for_index(options.bam)
    
    verboseprint("Starting peak calling")        
    main(options)

#so hacky... need to factor some of this out
if __name__ == "__main__":
    call_main()
