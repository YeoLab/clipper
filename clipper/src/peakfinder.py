#!/bin/env python

from __future__ import print_function
from collections import defaultdict
import gzip
import logging
import multiprocessing
import numpy as np
from optparse import OptionParser, SUPPRESS_HELP
import os
import pickle
from subprocess import call
import pandas as pd
import pybedtools

import clipper
from clipper.src.call_peak import call_peaks, poissonP

logging.captureWarnings(True)


def check_for_index(bamfile):
    
    """

    Checks to make sure a BAM file has an index, if the index does not exist it is created
    
    Usage undefined if file does not exist (check is made earlier in program)
    bamfile - a path to a bam file
    
    """

    if not os.path.exists(bamfile):
        raise NameError("file %s does not exist" % (bamfile))
    
    if os.path.exists(bamfile + ".bai"):
        return 
    if not bamfile.endswith(".bam"):
        raise NameError("file %s not of correct type" % (bamfile))
    else:
        logging.info("Index for %s does not exist, indexing bamfile" % (bamfile))
        
        process = call(["samtools", "index", str(bamfile)])
        
        if process == -11: 
            raise NameError("file %s not of correct type" % (bamfile))


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
 
##############################
# only used along with --debug, was used by github code, but no longer necessary
##############################
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


def build_transcript_data_gtf_as_structure(species, pre_mrna):
    
    """
    
    gtf_file - gtf file generated from AS_STRUCTURE_gtf ipython notebook 
    pre_mrna - if true uses pre mRNA length instead of mRNA length
    
    """
    bedtoolintervals = []
    x = clipper.data_file(species + ".AS.STRUCTURE.COMPILED.gff")
    gtf_file = pybedtools.BedTool(x)
    for gene in gtf_file:
        effective_length = gene.attrs['premrna_length'] if pre_mrna else gene.attrs['mrna_length']
        attrs = "gene_id=%s;" % (gene.attrs['gene_id'])
        if "transcript_ids" in gene.attrs:
            attrs += "transcript_ids=%s;" % (gene.attrs['transcript_ids']) 
        attrs += "effective_length=%s" % (str(effective_length)) 
        
        bedtoolintervals.append(pybedtools.create_interval_from_list(map(str, [gene['chrom'],
                                                                               "AS_STRUCTURE",
                                                                               "mRNA",
                                                                               str(gene.start + 1),
                                                                               str(gene.stop + 1),
                                                                               "0",
                                                                               gene['strand'],
                                                                               ".",
                                                                               attrs
                                                                               ])))

    return pybedtools.BedTool(bedtoolintervals)


def build_transcript_data_gtf(gtf_file, pre_mrna):
    
    """
    
    Generates GTF file to use when calling genes
    Returns the longest gene from a group of transcripts to call peaks on (this isn't optimal 
    behavior, but until we get a general as structure working its alright)
    
    gtf_file - bedtool from a standard gtf file
    pre_mrna - boolean flag to use pre_mrna instead of mrna 
    
    """
    
    #objects for default dict, no need to test or factor out
    def default_transcript():
        return {'chrom' : None, 'start': np.inf, "stop" : np.NINF, "strand" : None, "gene_id" : None, "mRNA_length" : 0}
    def default_gene():
        return {'start' : 0, 'stop' : 0}
    
    #get all transcripts, their starts, stops and mrna lengths
    transcripts = defaultdict(default_transcript)
    gtf_file = gtf_file.filter(lambda x: x[2] == 'exon')
    for interval in gtf_file:
        cur_transcript = transcripts[interval.attrs['transcript_id']]
        cur_transcript['start'] = min(cur_transcript['start'], interval.start)
        cur_transcript['stop'] = max(cur_transcript['stop'], interval.stop)
        cur_transcript['chrom'] = interval.chrom
        cur_transcript['strand'] = interval.strand
        cur_transcript['gene_id'] = interval.attrs['gene_id']
        cur_transcript['mRNA_length'] += interval.length
        cur_transcript['transcript_id'] = interval.attrs['transcript_id']
        
    #get the longest transcript from each gene group
    longest_genes = defaultdict(default_gene)
    for transcript_name, transcript in transcripts.items():
        cur_gene = transcript['gene_id']
        foo = longest_genes[cur_gene]
        best_length = longest_genes[cur_gene]['stop'] - longest_genes[cur_gene]['start']
        cur_length = transcript['stop'] - transcript['start']
        if  best_length < cur_length:
            longest_genes[cur_gene] = transcript
    
    #convert back into a gtf file 
    bedtoolintervals = []
    for gene in longest_genes.values():
        effective_length = gene['stop'] - gene['start'] if pre_mrna else gene['mRNA_length']
        bedtoolintervals.append(pybedtools.create_interval_from_list([gene['chrom'], 
                                        "AS_STRUCTURE", 
                                        "mRNA", 
                                        str(gene['start']), 
                                        str(gene['stop']),
                                        "0", 
                                        gene['strand'], 
                                        ".",
                                        "gene_id=" + gene['gene_id'] + "; transcript_id=" + gene['transcript_id'] + "; effective_length=" + str(effective_length)]))
    return pybedtools.BedTool(bedtoolintervals)


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
    
    #this is a stopgap until it can be fully factored out, returing a gtf file of 
    #genes and effective lengths, eventually this is the file we want to pass in
    gtf_list = []
    
    for gene in genes.keys():
        gtf_list.append(pybedtools.create_interval_from_list([genes[gene][0], 
                        "AS_STRUCTURE", 
                        "mRNA",
                        str(genes[gene][2]), 
                        str(genes[gene][3]),
                        ".",
                        str(genes[gene][4]),
                        ".",
                        "gene_id=" + gene + "; effective_length=" + str(lengths[gene])]))

    return pybedtools.BedTool(gtf_list)


def count_transcriptome_reads(peaks_dicts):
    
    """ 
    
    Counts number of reads in the entire transcriptome
    
    results -- the result returned back by call_peaks
    
    returns int, the number of reads in the transcriptome
    
    """
    ################################################
    # logging.info(" number of reads per gene result")
    ################################################
    #count total number of reads in transcriptiome
    transcriptome_reads = 0

    for peaks_dict_no, peaks_dict in enumerate(peaks_dicts):
        if peaks_dict is not None:
            if int(peaks_dict['nreads']) > 10 :
                logging.info("     gene_result_no: %s , number_of_reads: %d" % (peaks_dict_no, peaks_dict['nreads']))
            transcriptome_reads += peaks_dict['nreads']

    return transcriptome_reads


def count_transcriptome_length(results):
    transcriptome_length = 0

    for gene_result in results:
        if gene_result is not None:
            transcriptome_length += int(gene_result['loc'].attrs['effective_length'])

    return transcriptome_length


def transcriptome_poissonP(cluster):
    return poissonP(cluster.transcriptome_reads,
                    cluster.number_reads_in_peak,
                    cluster.transcriptome_size,
                    cluster['size'])


def transcript_poissonP(cluster):
    return poissonP(cluster.nreads_in_gene,
                    cluster.number_reads_in_peak,
                    cluster.effective_length,
                    cluster['size'])


def superlocal_poissonP(cluster):
    return poissonP(cluster.area_reads,
                    cluster.number_reads_in_peak,
                    cluster.area_size,
                    cluster['size'])


def write_peak_bedtool_string(cluster):
    cluster_info_list = [
        cluster.chrom,
        cluster.genomic_start,
        cluster.genomic_stop,
        cluster.gene_name  + "_" + str(cluster.peak_number) + "_" + str(cluster.number_reads_in_peak),
        cluster.final_p_value,
        cluster.strand,
        cluster.thick_start,
        cluster.thick_stop,
        ]
    cluster_bedtool_string = "\t".join([str(info) for info in cluster_info_list])
    return cluster_bedtool_string


def dictify(some_named_tuple):
    return dict((s, getattr(some_named_tuple, s)) for s in some_named_tuple._fields)


def make_peak_dataframe(peaks_dicts):
    peaks_list = [dictify(cluster) for peaks_dict in peaks_dicts for cluster in peaks_dict['clusters' ]]
    peaks_dataframe = pd.DataFrame(peaks_list)
    return peaks_dataframe


def bh_correct(df):
    """
    :param df:
    :return: returns dataframe wtih adjusted p-value
    """
    df = df.sort_values("final_p_value")
    df['sort_rank'] = np.arange(1, len(df) + 1)
    df['bh_corrected'] = df.apply(lambda x: min(((len(df) / x.sort_rank) * x.final_p_value), 1), axis=1)
    df['padj'] = df.sort_values("final_p_value", ascending=False).bh_corrected.cummin()
    return df.sort_index()

def filter_peaks_dicts(peaks_dicts, poisson_cutoff, transcriptome_size,
                   transcriptome_reads, use_global_cutoff, 
                   bonferroni_correct, algorithm="spline", superlocal=False, min_width=50,
                   bypassfiltering=False):
    # TODO bonferroni_correct is now always True!
    
    """
    
    Takes a list of peaks_dict's, filters them based off of various argunments and returns only the filtered
    reads
    
    options - the options object from the initial parsing
    poisson_cutoff - user defined possion cutoff (also from options) that filters reads
    results - list of results generated by call_peaks
    transcriptome_size - number of genes there are in the transcriptome
    
    """

    peaks_dataframe = make_peak_dataframe(peaks_dicts)
    total_clusters = len(peaks_dataframe)
    #########################################################################
    # logging.info(" total clusters BEFORE filtering : %s" % (total_clusters))
    #########################################################################
    if total_clusters == 0:
        logging.info(" no peaks detected in dataset")
        return []

    # TODO always False !
    if algorithm == "classic":
        peaks_dataframe['peak_length'] = peaks_dataframe['peak_length'].apply(lambda x: max(x, min_width))
    peaks_dataframe['transcriptome_size'] = transcriptome_size
    peaks_dataframe['transcriptome_reads'] = transcriptome_reads
    peaks_dataframe['transcriptome_poisson_p'] = peaks_dataframe.apply(transcriptome_poissonP, axis=1) if use_global_cutoff else np.nan
    peaks_dataframe['transcript_poisson_p'] = peaks_dataframe.apply(transcript_poissonP, axis=1)
    peaks_dataframe['superlocal_poisson_p'] = peaks_dataframe.apply(superlocal_poissonP, axis=1) if superlocal else np.nan

    if algorithm == "classic":
    # TODO this never happens!
        peaks_dataframe['final_p_value'] = peaks_dataframe[['transcript_poisson_p', 'superlocal_poisson_p']].max(axis=1)
    else:
    # TODO this is always the case!
        peaks_dataframe['final_p_value'] = peaks_dataframe[['transcript_poisson_p', 'superlocal_poisson_p']].min(axis=1)

    if bonferroni_correct:
    # TODO always True!
        peaks_dataframe = bh_correct(peaks_dataframe)
        #peaks_dataframe['final_p_value'] = (peaks_dataframe['final_p_value'] * total_clusters)    # TODO still needed?

    if bypassfiltering :
        filtered_peaks_dataframe = peaks_dataframe
    else:
    #This is a bug I should fix, padj isn't getting printed, the uncorreded p-value is
        filtered_peaks_dataframe = peaks_dataframe[peaks_dataframe['padj'] < poisson_cutoff]

    #########################################################################
    # logging.info(" total clusters AFTER filtering : %s" % len(filtered_peaks_dataframe))
    #########################################################################
    filtered_peak_bedtool_strings = filtered_peaks_dataframe.apply(write_peak_bedtool_string, axis=1).values
    filtered_peak_bedtool_tsv = "\n".join(filtered_peak_bedtool_strings) + "\n"
    return filtered_peak_bedtool_tsv


def get_exon_bed(species):

    short_species = species.split("_")[0]
    return os.path.join(clipper.data_dir(), "regions", "%s_%s.bed" % (short_species, "exons"))


def main(options):
    ##############################################
    # logging.info("options : {}".format(options))
    ##############################################
    check_for_index(options.bam)
    
    if options.np == 'autodetect':
        options.np = multiprocessing.cpu_count()

    pool = multiprocessing.Pool(int(options.np))
        
    bamfile = options.bam
    
    if os.path.exists(bamfile):
        #re-set to include the full path to bamfile
        bamfile = os.path.abspath(bamfile) 
        logging.info("bam file is set to %s\n" % (bamfile))
    else:
        logging.error("Bam file: %s is not defined" % (bamfile))
        raise IOError
    
    if options.gtfFile:
        # TODO always False - no longer an option
        bedtool = build_transcript_data_gtf(pybedtools.BedTool(options.gtfFile), options.premRNA)
    else:
        bedtool = build_transcript_data_gtf_as_structure(options.species, options.premRNA)
    bedtool.saveas()

    #gets a bedtool of all genes to call peaks on
    if options.gene:
        bedtool = bedtool.filter(lambda x: x.attrs['gene_id'] in options.gene)

    # options.maxgenes   # truncates for max bedtool
    if options.maxgenes:
        ########################################################################
        logging.info(" number of genes before maxing : {}".format(len(bedtool)))
        logging.info(" max genes from user input: {}".format(options.maxgenes))
        ########################################################################
        if options.maxgenes < len(bedtool):
            bedtool = bedtool.random_subset(int(options.maxgenes))
        else:
            logging.info(" number of genes <= max genes from user , not truncating genes")
            pass

    exons = get_exon_bed(options.species)

    bedtool = bedtool.saveas()

    tasks = [(bedtool_interval, bedtool_interval.attrs['effective_length'], bamfile, options.max_gap, options.FDR_alpha,
              options.threshold, options.binom, options.method, options.minreads, options.poisson_cutoff,
              options.plotit, 10, 1000, options.SloP, options.max_width,
              options.min_width, options.algorithm,
              # TODO options.algorithm now always "spline" !
              options.reverse_strand, exons
              )
             for gene_no, bedtool_interval in enumerate(bedtool)]
    ##################################
    # print("len(tasks):", len(tasks))
    ##################################

    #jobs = []
    peaks_dicts = []
    # generate list of all peaks_dict's, (one peaks_dict per gene)
    ##############################################################
    if options.debug:
        peaks_dicts = [call_peaks(*task) for task in tasks]
    
    else:        
        jobs = [pool.apply_async(call_peaks, task) for task in tasks]
        
        for job, task in zip(jobs, tasks):
            try:
                peaks_dicts.append(job.get(timeout=options.timeout))
            except multiprocessing.TimeoutError as error:
                print()
                ####################################################################################################################################
                logging.error("gene %s timed out after %s minutes on bedinterval: %s" % (task[0].attrs['gene_id'], options.timeout / 60,  task[0] ))
                ####################################################################################################################################
        
    pool.close()
    #################################################################
    logging.info("finished call_peaks on all genes")
    #################################################################

    ############################################################
    logging.info(" starting adding up transcriptome-wise reads")
    ############################################################
    transcriptome_reads = count_transcriptome_reads(peaks_dicts)
    transcriptome_size = count_transcriptome_length(peaks_dicts)
    ####################################################################################
    logging.info(" transcriptome size in bases: {}".format(transcriptome_size))
    logging.info(" transcriptome total number of reads: {}".format(transcriptome_reads))
    ####################################################################################
    filtered_peak_bedtool_tsv = filter_peaks_dicts(peaks_dicts,
                              options.poisson_cutoff, 
                              transcriptome_size,  
                              transcriptome_reads, 
                              options.use_global_cutoff,
                              options.bonferroni_correct,
                              options.algorithm,
                              options.SloP,
                              options.min_width,
                              bypassfiltering=False)

    ##########################################################
    # logging.info(" 1: {}".format(filtered_peak_bedtool_tsv))
    ##########################################################


    ###############
    # writing files
    ###############

    # options.outfileF, options.save_pickle
    #======================================
    outbedF = options.outfileF
    wether_to_save_pickle = options.save_pickle
    #
    # writing tsv files
    #==================
    with open(outbedF + ".tsv", 'w') as tsvfile:
        tsvfile.write(filtered_peak_bedtool_tsv)
    #
    # writing bed files
    #==================
    pybedtools.BedTool(filtered_peak_bedtool_tsv, from_string=True).sort(stream=True).saveas(outbedF)
    ########################################################
    #logging.info(" wrote filtered peaks to %s" % (outbedF))
    ########################################################
    #
    # writing pickle files
    #=====================
    if wether_to_save_pickle is True:
        with open(outbedF + ".pickle", 'w') as picklefile:                     # TODO Can't pickle save after filtering ? as we have a tsv now, not a peaks_dicts list !?
            pickle.dump(peaks_dicts, file=picklefile)



def call_main():
    
    usage = """
    THIS IS CLIPPER FOR ECLIP VERSION 0.1.4
    python peakfinder.py -b <bamfile> -s <hg18/hg19/mm9> OR 
    python peakfinder.py -b <bamfile> --customBED <BEDfile> --customMRNA 
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
    #parser.add_option("--gtfFile", dest="gtfFile", help="use a gtf file instead of the AS structure data")
    parser.add_option("--outfile", "-o", dest="outfileF", default="fitted_clusters", help="a bed file output, default:%default")
    parser.add_option("--gene", "-g", dest="gene", action="append", help="A specific gene you'd like try", metavar="GENENAME")
    parser.add_option("--minreads", dest="minreads", help="minimum reads required for a section to start the fitting process.  Default:%default", default=3, type="int", metavar="NREADS")
    #parser.add_option("--premRNA", dest="premRNA", action="store_true", help="use premRNA length cutoff, default:%default", default=False)
    parser.add_option("--poisson-cutoff", dest="poisson_cutoff", type="float", help="p-value cutoff for poisson test, Default:%default", default=0.05, metavar="P")
    parser.add_option("--disable_global_cutoff", dest="use_global_cutoff", action="store_false", help="disables global transcriptome level cutoff to CLIP-seq peaks, Default:On", default=True, metavar="P")
    parser.add_option("--FDR", dest="FDR_alpha", type="float", default=0.05, help="FDR cutoff for significant height estimation, default=%default")
    #parser.add_option("--threshold-method", dest="method", default="binomial", help="Method used for determining height threshold, Can use default=random or binomial")
    parser.add_option("--binomial", dest="binom", type="float", default=0.05, help ="Alpha significance threshold for using Binomial distribution for determining height threshold, default=%default")
    parser.add_option("--threshold", dest="threshold", type="int", default=None, help="Skip FDR calculation and set a threshold yourself")
    parser.add_option("--maxgenes", dest="maxgenes", default=None, type="int", help="stop computation after this many genes, for testing", metavar="NGENES")
    parser.add_option("--processors", dest="np", default="autodetect", help="Number of processors to use. Default: All processors on machine", type="str", metavar="NP")
    #parser.add_option("--superlocal", action="store_true", dest="SloP", default=False, help="Use super-local p-values, counting reads in a 1KB window around peaks")
    parser.add_option("--plot", "-p", dest="plotit", action="store_true", help="make figures of the fits", default=False)
    parser.add_option("--verbose", "-v", dest="verbose", action="store_true", default=False)                                            # TODO unused ?
    parser.add_option("--quiet", "-q", dest="quiet", action="store_true", default=False, help="suppress notifications")                 # TODO unused ?
    parser.add_option("--save-pickle", dest="save_pickle", default=False, action="store_true", help="Save a pickle file containing the analysis")
    parser.add_option("--debug", dest="debug", default=False, action="store_true", help="disables multipcoressing in order to get proper error tracebacks")
    #parser.add_option("--max_width", dest="max_width", type="int", default=75, help="Defines max width for classic algorithm, default: %default")
    #parser.add_option("--min_width", dest="min_width", type="int", default=50, help="Defines min width for classic algorithm, default: %default")
    parser.add_option("--max_gap", dest="max_gap",type="int", default=15, help="defines maximum gap between reads before calling a region a new section, default: %default")
    #parser.add_option("--bonferroni", dest="bonferroni_correct",action="store_true", default=False, help="Perform Bonferroni on data before filtering")
    #parser.add_option("--algorithm", dest="algorithm",default="spline", help="Defines algorithm to run, currently spline, classic, gaussian")
    #parser.add_option("--reverse_strand", dest="reverse_strand",default=False, action="store_true", help="adds option to reverse strand")
    parser.add_option("--timeout", dest="timeout", default=None, type=int, help="adds timeout (in seconds) to genes that take too long (useful for debugging only, or if you don't care about higly expressed genes)")


    (options, args) = parser.parse_args()
 
    ##########################################################################
    #logfile = options.outfileF + '.log'
    #logging.basicConfig(filename=logfile ,level=logging.DEBUG)
    #logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
    ##########################################################################
    if options.plotit:
        options.debug=True

    options.premRNA = True
    options.gtfFile = None
    # overrides parser option --threshold-method
    options.method = "binomial"
    # overrides parser option --superlocal
    options.SloP = True
    # overrides parser option --bonferroni
    options.bonferroni_correct = True
    options.algorithm = "spline"
    options.reverse_strand = False
    options.max_width = 75
    options.min_width = 50

    #enforces required usage
    if not (options.bam and ((options.species) or (options.gtfFile))): 
        parser.print_help()
        exit()

    ################################################
    logging.info("starting peakfinder on all genes")
    ################################################
    main(options)
    ################################################
    logging.info(" finished peakfinder on all genes")
    ################################################


if __name__ == "__main__":
    call_main()
