from pysam import rmdup
from optparse import OptionParser, SUPPRESS_HELP
import os
import sys
from subprocess import call
import pickle
import time
import pybedtools
import gzip
import math
import multiprocessing 
import clipper
from clipper import data_dir
from clipper.src.call_peak import call_peaks, poissonP
import logging
logging.captureWarnings(True)
    
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
    outfile += ".rmdup.bam"
    rmdup("-S", bamfile, outfile)
    return outfile

def check_for_index(bamfile):
    
    """

    Checks to make sure a BAM file has an index, if the index does not exist it is created
    
    Usage undefined if file does not exist (check is made earlier in program)
    bamfile - a path to a bam file
    
    Returns 1 

    """

    if not os.path.exists(bamfile):
        raise NameError("file %s does not exist" % (bamfile))
    
    if os.path.exists(bamfile + ".bai"):
        return 1
    else:
        logging.error("Index for %s does not exist, indexing bamfile" % (bamfile))
        #result = pysam.index(str(bamfile))
#TODO fix this, indexing is very fragile
        process = call(["samtools", "index", str(bamfile)])
        
        if process == -11: 
            raise NameError("file %s not of correct type" % (bamfile))
        
        return 1

def build_geneinfo(bed):
    
    """

    Loads bed file into a dictionary with the key being the name and a string being the value
    
    Input:
    BED -- a bed file to load
    
    Return:
    A dictionary with the key being the name position of the bed file and the values being the
    ordered bed file
    
    """
    
    #opens bed file, either zipped or unzipped
    try:
        bedfile = gzip.open(bed, "rb")
    except IOError:
        bedfile = open(bed, "r")
        
    gene_info = dict()
    
    for line in bedfile.readlines():
        chromosome, start, stop, name, score, signstrand = line.strip().split()
        gene_info[name] = [chromosome, name, int(start), 
                           int(stop), str(signstrand)]
    
    bedfile.close()
    return gene_info

def build_lengths(length_file):
    
    """
    
    Builds a dictionary of gene names and lengths of mappable regions in that gene
    
    Input:
    A two column file with the first column being the gene name and the second column being the
    mappable length of the gene
    
    Return:
    A dictionary with the key being the name of the gene and the value being the length
    
    """
    
    try:
        handle = open(length_file, "r")
        gene_lengths = {}
    
        for line in handle.readlines():
            name, gene_length = line.strip().split("\t")
            gene_lengths[name] = int(gene_length)
    
        handle.close()
        
    except TypeError:
        raise ValueError("file %s not found" % length_file)
    except ValueError:
        raise ValueError("file not formatted correctly, expects two columns gene<tab>length")
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
    #expand sublists
    par["chrs"] = [item for sublist in chrs for item in sublist] 
    par["gene_bed"] = bed
    par["mRNA"] = mrna
    par["premRNA"] = premrna
    return par
 
def func_star(varables):
    """ covert f([1,2]) to f(1,2) """
    return call_peaks(*varables)


def get_acceptable_species():
    
    """
    
    Finds all species in data directory 
    
    """
    
    acceptable_species = set([])
    for fn in os.listdir(clipper.data_dir()):
        fn = fn.split(".")[0]
        
        if fn == "__init__":
            continue
        
        acceptable_species.add(fn)
    
    return acceptable_species

 
def build_transcript_data_bed(bed_file, pre_mrna):
    
    """
    
    Generates gene lengths and gene names from BED12 file
    
    bed_file - pybedtools object defining genes
    pre_mrna - flag True indicates use pre-mRNA lengths instead of mRNA lengths 
    
    TODO just turn this totally into a pybedtools object instead of converting 
    it to a dictionary
    """
    
    gene_info = {}
    gene_lengths = {}
    
    for line in bed_file:
        
        #builds gene info
        gene_info[line.name] = [line.chrom, line.name, line.start, line.stop, line.strand]
        #builds gene lengths
        
        if pre_mrna:
            gene_lengths[line.name] = line.stop - line.start
        else:
            #Just gets the lengths of the exons (although no mention of cds or not... not important)
            gene_lengths[line.name] = sum([int(x) for x in line[10][:-1].strip().split(",")])
            
    return gene_info, gene_lengths

def build_transcript_data(species, gene_bed, gene_mrna, gene_pre_mrna, pre_mrna):
    
    """
    
    Generates transcript data structures to call peaks on
    
    Allows for either predefined files (from the data directory) 
    or custom files
    
    Accepts species, and genebed, genemrnaand genepremrna options
    
    species - the species to run on
    gene_bed - an abribtary bed file of locations to search for peaks (should be gene locations)
    gene_mrna - the effective length of the mrna of a gene (unmappable regions removed)
    gene_premrna - the effective length of the pre-mrna (unmappable regions removed)
    pre_mrna - flag True indicates use pre-mRNA lengths instead of mRNA lengths
     
    returns genes and lengths dict
    
    """
    
    #error checking 
    
    acceptable_species = get_acceptable_species()
    if (species is None and 
        gene_bed is None and 
        (gene_mrna is None or gene_pre_mrna is None)):
        
        raise ValueError("You must set either \"species\" or \"geneBed\"+\"geneMRNA\"+\"genePREMRNA\"")

    if species is not None and gene_bed is not None:
        raise ValueError("You shouldn't set both geneBed and species, defaults exist for %s" % (acceptable_species))
    
    #Now actually assign values
    if species is not None:
        try:
            gene_bed      = clipper.data_file(species + ".AS.STRUCTURE_genes.BED.gz")
            gene_mrna     = clipper.data_file(species + ".AS.STRUCTURE_mRNA.lengths")
            gene_pre_mrna = clipper.data_file(species + ".AS.STRUCTURE_premRNA.lengths")
            
        except ValueError:
            raise ValueError("Defaults don't exist for your species: %s. Please choose from: %s or supply \"geneBed\"+\"geneMRNA\"+\"genePREMRNA\"" % (species, acceptable_species))

    #Selects mRNA or preMRNA lengths
    if pre_mrna is True:
        lenfile = gene_pre_mrna
    else:
        lenfile = gene_mrna

    if lenfile is None:
        raise IOError("""didn't pass correct mRNA length file option 
                    with given length file""")
        
    #builds dict to do processing on,
    genes = build_geneinfo(gene_bed)
    lengths = build_lengths(lenfile)
    
    
    return genes, lengths



def transcriptome_filter(poisson_cutoff, transcriptome_size, transcriptome_reads, cluster):
    
    """
    filters each cluster by if it passes a transciptome wide cutoff or not, returns true if it passes
    transcriptome cutoff, false if not
    
    poisson_cutoff - float,user set cutoff 
    transcriptome_size - int number of genes in transcriptome
    transcritpmoe_reads - int total number of reads analized
    cluster - named tuple , namedtuple('Peak', ['chrom', 
                                      'genomic_start', 
                                      'genomic_stop', 
                                      'gene_name', 
                                      'super_local_poisson_p', 
                                      'strand',
                                      'thick_start',
                                      'thick_stop',
                                      'peak_number',
                                      'number_reads_in_peak',
                                      'gene_poisson_p',
                                      'size'
                                      'p'
                                      ])
    """
    
    transcriptome_p = poissonP(transcriptome_reads, 
                               cluster.number_reads_in_peak, 
                               transcriptome_size, 
                               cluster.size)
    
    if math.isnan(transcriptome_p):
        logging.info("""Transcriptome P is NaN, transcriptome_reads = %d, cluster reads = %d, transcriptome_size = %d, cluster_size = %d""" % (transcriptome_reads, cluster.number_reads_in_peak, transcriptome_size, cluster.size))
        return False
    
    if transcriptome_p > poisson_cutoff:
        logging.info("%s\n Failed Transcriptome cutoff with %s reads, pval: %s" % (cluster, cluster.number_reads_in_peak, transcriptome_p))

        return False
    
    return True


def count_transcriptome_reads(results):
    
    """ 
    
    Counts number of reads in the entire transcriptome
    
    results -- the result returned back by call_peaks
    
    returns int, the number of reads in the transcriptome
    
    """
    #count total number of reads in transcriptiome
    transcriptome_reads = 0

    for gene_result in results:
        if gene_result is not None:
            logging.info("nreads: %d" % (gene_result['nreads']))
            transcriptome_reads += gene_result['nreads']
    
    
    return transcriptome_reads

def filter_results(results, poisson_cutoff, transcriptome_size, transcriptome_reads, use_global_cutoff):
    
    """
    
    Takes a list of results, filters them based off of various argunments and returns only the filtered
    reads
    
    options - the options object from the initial parsing
    poisson_cutoff - user defined possion cutoff (also from options) that filters reads
    results - list of results generated by call_peaks
    transcriptome_size - number of genes there are in the transcriptome
    
    """
    #combine results
    allpeaks = set([])
        
    for gene_result in results:

        #alert user that there aren't any clusters for specific gene
        if gene_result['clusters'] is None:
            logging.info("%s no clusters" % (gene_result))

        
        for cluster in gene_result['clusters']:
            meets_cutoff = True
            try:
                
                if use_global_cutoff:
                    meets_cutoff = meets_cutoff and transcriptome_filter(poisson_cutoff, 
                                                                         transcriptome_size, 
                                                                         transcriptome_reads,  
                                                                         cluster)
                
                #should factor out this as well, but I'll leave it be until nessessary            
                #does SlOP always get used?  it looks like it does
                corrected_SloP_pval = cluster.super_local_poisson_p
                corrected_gene_pval = cluster.gene_poisson_p
                min_pval = min([corrected_SloP_pval, corrected_gene_pval])
                
                if not (min_pval < poisson_cutoff):
                    logging.info("Failed Gene Pvalue: %s and failed SloP Pvalue: %s for Gene %s %s" % (corrected_gene_pval, corrected_SloP_pval, cluster.gene_name, cluster.peak_number))
                    meets_cutoff = False
                
                if meets_cutoff:

                    #adds beadline to total peaks that worked
                    #the name column is hacky because the ucsc genome browser doesn't like
                    #random columns, will need to do some parsing later to fix (awk it)
                    allpeaks.add("\t".join([str(x) for x in [
                                           cluster.chrom, 
                                           cluster.genomic_start, 
                                           cluster.genomic_stop, 
                                           cluster.gene_name  + "_" + str(cluster.peak_number) + "_" + str(cluster.number_reads_in_peak), 
                                           min_pval, 
                                           cluster.strand, 
                                           cluster.thick_start, 
                                           cluster.thick_stop,
                                           ]]))
        
            except NameError as error:
                logging.error("parsing failed: %s" % (error))
                raise error
            
    return allpeaks

def mapper(options, line):
    bedtool = pybedtools.BedTool(line, from_string=True)

    #loads in the last bedline, because bedtools doesn't have a .next()
    for bedline in bedtool:
        pass
    
    if options.premRNA:
        length = bedline.stop - bedline.start
    else:
        length = sum([int(x) for x in bedline[10][:-1].strip().split(",")]) #Just gets the lengths of the exons (although no mention of cds or not... not important)
    
    print call_peaks([bedline.chrom, bedline.name, bedline.start, bedline.stop, bedline.strand], length, None, options.bam, int(options.margin), options.FDR_alpha, 
        options.threshold, int(options.minreads), options.poisson_cutoff, 
        options.plotit, 10, 1000, options.SloP, False)
    
def hadoop_mapper(options):
    
    """
    
    Expermental mapper to give call peaks (running call_peak function) from hadoop
    
    """
    #being lazy for now will make this integrated eventually
   
    for line in sys.stdin:
        mapper(options, line)
        
def main(options):
    
    check_for_index(options.bam)
    
    if options.np == 'autodetect':
        options.np = multiprocessing.cpu_count()

    pool = multiprocessing.Pool(int(options.np))
    
    #job_server = pp.Server(ncpus=options.np) #old pp stuff
    
    bamfile = options.bam
    
    if os.path.exists(bamfile):
        #re-set to include the full path to bamfile
        bamfile = os.path.abspath(bamfile) 
        logging.info("bam file is set to %s\n" % (bamfile))
    else:
        logging.error("Bam file: %s is not defined" % (bamfile))
        raise IOError
    
    if options.bedFile is not None:
        genes, lengths = build_transcript_data_bed(options.bedFile, options.premRNA)
    else:
        genes, lengths = build_transcript_data(options.species, 
                                               options.geneBEDfile, 
                                               options.geneMRNAfile, 
                                               options.genePREMRNAfile,
                                               options.premRNA)
    
    margin = int(options.margin)
    
    #this should be fixed, args should initally be ints if passed
    if options.maxgenes is not None:
        maxgenes = int(options.maxgenes)

    minreads = int(options.minreads)
    poisson_cutoff = options.poisson_cutoff

    #gets all the genes to call peaks on

    if options.gene is not None and len(options.gene) > 0:
        gene_list = options.gene
    else: #selects all genes
        gene_list = genes.keys()
    

    results = []
    
    #Set up peak calling by gene
    running_list = [genes[gene] for gene in gene_list]
    length_list  = [lengths[gene] for gene in gene_list]
    
    #truncates for max genes
    if options.maxgenes is not None:
        import random

        ind = random.sample(xrange(len(genes.keys())), k=options.maxgenes)
        running_list = running_list[ind]
        length_list  = length_list[ind]
    
    transcriptome_size = sum(length_list)
    #do the parralization
    tasks =  [(gene, length, None, bamfile, margin, options.FDR_alpha, 
               options.threshold, minreads, poisson_cutoff, 
               options.plotit, 10, 1000, options.SloP, False,
               options.max_width, options.min_width, options.max_gap,
               options.algorithm)
              for gene, length in zip(running_list, length_list)]
    
    jobs = []
    if options.debug:
        
        for job in tasks:
            jobs.append(func_star(job))
        
        for job in jobs:
            results.append(job)   
    
    else:
        #sets chunk size to be a fair bit smaller, than total input, but not
        #to small
        chunk_size = len(tasks) // int(options.np) * 10
        if chunk_size < 1:
            chunk_size = 1
        
      
        for job in tasks:
            jobs.append(pool.apply_async(call_peaks, job))
        
        for job in jobs:
            try:
                results.append(job.get(timeout=360))
            except Exception as error:
                logging.error(error)
        #jobs = pool.map(func_star, tasks, chunksize=chunk_size)
        
    pool.close()
    

    logging.info("finished with calling peaks")
    #if we are going to save and output as a pickle file we should 
    #output as a pickle file we should factor instead create a method 
    #or object to handle all file output
    if options.save_pickle is True:
        pickle_file = open(options.outfile + ".pickle", 'w')
        pickle.dump(results, file=pickle_file)                
    
    transcriptome_reads = count_transcriptome_reads(results)
    
    logging.info("""Transcriptome size is %d, transcriptome reads are %d""" % (transcriptome_size, transcriptome_reads))

    filtered_peaks = filter_results(results, 
                              poisson_cutoff, 
                              transcriptome_size,  
                              transcriptome_reads, 
                              options.use_global_cutoff)
    
    
    #need to include scaliable way to see ALL peaks, not just filtered one
        
    outbed = options.outfile
    color = options.color
    pybedtools.BedTool("\n".join(filtered_peaks), from_string=True).sort(stream=True).saveas(outbed, trackline="track name=\"%s\" visibility=2 colorByStrand=\"%s %s\"" % (outbed, color, color))

    logging.info("wrote peaks to %s" % (options.outfile))
    
def call_main():
    
    usage = """\npython peakfinder.py -b <bamfile> -s <hg18/hg19/mm9>\n OR 
    \npython peakfinder.py -b <bamfile> --customBED <BEDfile> --customMRNA 
    <mRNA lengths> --customPREMRNA <premRNA lengths>"""
    description = """CLIPper. Michael Lovci, Gabriel Pratt 2012. 
                     CLIP peakfinder that uses fitted smoothing splines to 
                     define clusters of binding.  Computation is performed in
                     parallel using parallelPython. 
                     Refer to: https://github.com/YeoLab/clipper/wiki for instructions. 
                     Questions should be directed to michaeltlovci@gmail.com."""

    parser = OptionParser(usage=usage, description=description)

    parser.add_option("--bam", "-b", dest="bam", help="A bam file to call peaks on", type="string", metavar="FILE.bam")

    parser.add_option("--species", "-s", dest="species", help="A species for your peak-finding, either hg19 or mm9")
    
    #we don't have custom scripts or documentation to support this right now, removing until those get added in
    parser.add_option("--customBED", dest="geneBEDfile", help="bed file to call peaks on, must come withOUT species and with customMRNA and customPREMRNA", metavar="BEDFILE")
    parser.add_option("--customMRNA", dest="geneMRNAfile", help="file with mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")
    parser.add_option("--customPREMRNA", dest="genePREMRNAfile", help="file with pre-mRNA lengths for your bed file in format: GENENAME<tab>LEN", metavar="FILE")
    parser.add_option("--outfile", "-o", dest="outfile", default="fitted_clusters", help="a bed file output, default:%default")
    parser.add_option("--gene", "-g", dest="gene", action="append", help="A specific gene you'd like try", metavar="GENENAME")
    parser.add_option("--minreads", dest="minreads", help="minimum reads required for a section to start the fitting process.  Default:%default", default=3, type="int", metavar="NREADS")
    parser.add_option("--margin", dest="margin", type="int", help="find sections of genes within M bases that have genes and perform fitting. Default:%default", default=15, metavar="NBASES")
    parser.add_option("--trim", dest="trim", action="store_true", default=False, help="Trim reads with the same start/stop to count as 1")
    parser.add_option("--premRNA", dest="premRNA", action="store_true", help="use premRNA length cutoff, default:%default", default=False)
    parser.add_option("--poisson-cutoff", dest="poisson_cutoff", type="float", help="p-value cutoff for poisson test, Default:%default", default=0.05, metavar="P")
    parser.add_option("--disable_global_cutoff", dest="use_global_cutoff", action="store_false", help="disables global transcriptome level cutoff to CLIP-seq peaks, Default:On", default=True, metavar="P")
    parser.add_option("--FDR", dest="FDR_alpha", type="float", default=0.05, help="FDR cutoff for significant height estimation, default=%default")
    parser.add_option("--threshold", dest="threshold", type="int", default=None, help="Skip FDR calculation and set a threshold yourself")
    parser.add_option("--maxgenes", dest="maxgenes", default=None, help="stop computation after this many genes, for testing", metavar="NGENES")
    parser.add_option("--processors", dest="np", default="autodetect", help="Number of processors to use. Default: All processors on machine", type="str", metavar="NP")
    parser.add_option("--superlocal", action="store_true", dest="SloP", default=False, help="Use super-local p-values, counting reads in a 1KB window around peaks")
    parser.add_option("--color", dest="color", default="0,0,0", help="R,G,B Color for BED track output, default:black (0,0,0)")
    parser.add_option("--plot", "-p", dest="plotit", action="store_true", help="make figures of the fits", default=False)
    parser.add_option("--verbose", "-q", dest="verbose", action="store_true", help="suppress notifications")
    parser.add_option("--save-pickle", dest="save_pickle", default=False, action="store_true", help="Save a pickle file containing the analysis")
    parser.add_option("--debug", dest="debug", default=False, action="store_true", help="disables multipcoressing in order to get proper error tracebacks")
    parser.add_option("--bedFile", dest="bedFile", help="use a bed file instead of the AS structure data")
    parser.add_option("--max_width", dest="max_width", type="int", default=75, help="Defines max width for classic algorithm")
    parser.add_option("--min_width", dest="min_width", type="int", default=50, help="Defines min width for classic algorithm")
    parser.add_option("--max_gap", dest="max_gap",type="int", default=10, help="defines min gap for classic algorithm")
    parser.add_option("--algorithm", dest="algorithm",default="spline", help="Defines algorithm to run, currently Spline, or Classic")
    parser.add_option("--hadoop", dest="hadoop",default=False, action="store_true", help="Run in hadoop mode")
    
    (options, args) = parser.parse_args()
 
    if options.plotit:
        options.debug=True
    
    #enforces required usage    
    if not (options.bam and ((options.species) or (options.geneBEDfile and options.geneMRNAfile and options.genePREMRNAfile))):
        parser.print_help()
        exit()
    
    #If triming option is set use samtools to remove duplicate 
    #reads for us, trims strictly ignoring paired end and strandness
    if options.trim:
        options.bam = trim_reads(options.bam)
    
    logging.info("Starting peak calling")
    
    if options.hadoop:
        hadoop_mapper(options)
    else:
        main(options)


if __name__ == "__main__":
    call_main()
