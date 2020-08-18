'''
Created on Jul 25, 2012
@author: mlovci
@author: gabrielp

Refactored on June 29, 2020
@author: algaebrown

'''
from __future__ import print_function

############################################################
####### This files calls the whole pipeline ################
############################################################
from clipper.src.utils import check_for_index, build_transcript_data_gtf_as_structure, get_exon_bed
from clipper.src.call_peak import call_peaks
from clipper.src.filter_peak import count_transcriptome_reads, count_transcriptome_length, filter_peaks_dicts
import multiprocessing
import pybedtools
import os
import pickle
from optparse import OptionParser
import logging

logging.captureWarnings(True)

def main(options):
    """
    Run the whole pipeline
    :rtype: None
    """
    ##############################################
    # logging.info("options : {}".format(options))
    ##############################################

    ############ CHECKING FILE STATUS ############
    check_for_index(options.bam)
    bamfile = options.bam

    if os.path.exists(bamfile):
        # re-set to include the full path to bamfile
        bamfile = os.path.abspath(bamfile)
        logging.info("bam file is set to %s\n" % (bamfile))
    else:
        logging.error("Bam file: %s is not defined" % (bamfile))
        raise IOError

    ########### PREPARE GENE LENGTH ################
    # if options.gtfFile:
    #    # TODO always False - no longer an option
    #   bedtool = build_transcript_data_gtf(pybedtools.BedTool(options.gtfFile), options.premRNA)
    # else:
    bedtool = build_transcript_data_gtf_as_structure(options.species, options.premRNA).saveas()

    # gets a bedtool of all genes to call peaks on
    if options.gene:
        bedtool = bedtool.filter(lambda x: x.attrs['gene_id'].split('.')[0] in options.gene).saveas()  ### bug

    # options.maxgenes   # truncates for max bedtool
    if options.maxgenes:
        logging.info(" number of genes before maxing : {}".format(len(bedtool)))
        logging.info(" max genes from user input: {}".format(options.maxgenes))
        ########################################################################
        if options.maxgenes < len(bedtool):
            bedtool = bedtool.random_subset(int(options.maxgenes)).saveas()
        else:
            logging.info(" number of genes <= max genes from user , not truncating genes")
            pass

    if len(bedtool) == 0:
        raise Warning('Bedtool length is 0; check gene id')

    exons = get_exon_bed(options.species)

    ############### PREPARE MULTIPROCESSING ##############################

    if options.np == 'autodetect':
        options.np = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(int(options.np))

    tasks = [(bedtool_interval, bedtool_interval.attrs['effective_length'], bamfile, options.max_gap, options.FDR_alpha,
              options.threshold, options.binom, options.method, options.minreads, options.poisson_cutoff,
              options.plotit, 10, 1000, options.SloP, options.max_width,
              options.min_width, options.algorithm,
              # TODO options.algorithm now always "spline" !
              options.reverse_strand, exons
              )
             for gene_no, bedtool_interval in enumerate(bedtool)]

    logging.info('Total tasks: {}'.format(len(tasks)))

    ############## CALL PEAKS BY HEIGHT AND CURVE##########################

    # generate list of all peaks_dict's, (one peaks_dict per gene)
    peaks_dicts = []

    if options.debug:
        peaks_dicts = [call_peaks(*task) for task in tasks]

    else:
        jobs = [pool.apply_async(call_peaks, task) for task in tasks]

        for job, task in zip(jobs, tasks):
            try:
                peaks_dicts.append(job.get(timeout=options.timeout))
            except multiprocessing.TimeoutError as error:
                print()

                logging.error("gene %s timed out after %s minutes on bedinterval: %s" % (
                task[0].attrs['gene_id'], options.timeout / 60, task[0]))


    pool.close()
    logging.info("finished call_peaks on all genes")


    ################### FILTER PEAK BY READ #################################
    logging.info(" starting adding up transcriptome-wise reads")

    transcriptome_reads = count_transcriptome_reads(peaks_dicts)
    transcriptome_size = count_transcriptome_length(peaks_dicts)

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

    ############### WRITE TO FILE #####################################

    if type(filtered_peak_bedtool_tsv) == str:
        # with open(outbedF + ".tsv", 'w') as tsvfile:
        #    tsvfile.write(filtered_peak_bedtool_tsv)
        # filtered_peak_bedtool_dataframe.to_csv(tsvfile, sep = '\t')
        pybedtools.BedTool(filtered_peak_bedtool_tsv, from_string=True).sort(stream=True).saveas(options.outfileF)

    if options.save_pickle is True:
        with open(options.outfileF + ".pickle", 'w') as f:
            # TODO Can't pickle save after filtering ? as we have a tsv now, not a peaks_dicts list !?
            pickle.dump(peaks_dicts, file=f)

def option_parser():
    ''' return parser
    :return: OptionParser object
    '''
    usage = """
        THIS IS CLIPPER FOR ECLIP VERSION 2.0.0
        clipper -b YOUR_BAM_FILE.bam -o YOUR_OUT_FILE.bed -s hg19 """
    description = """CLIPper. Michael Lovci, Gabriel Pratt 2012, Hsuan-lin Her 2020.
                         CLIP peakfinder that uses fitted smoothing splines to 
                         define clusters of binding.  Computation is performed in
                         parallel using parallelPython. 
                         Refer to: https://github.com/YeoLab/clipper/ for instructions. 
                         Questions should be directed to hsher@ucsd.edu"""
    parser = OptionParser(usage=usage, description=description)

    parser.add_option("--bam", "-b", dest="bam", help="A bam file to call peaks on", type="string", metavar="FILE.bam")

    parser.add_option("--species", "-s", dest="species", help="A species for your peak-finding, either hg19 or mm9")
    # parser.add_option("--gtfFile", dest="gtfFile", help="use a gtf file instead of the AS structure data")
    parser.add_option("--outfile", "-o", dest="outfileF", default="fitted_clusters",
                      help="a bed file output, default:%default")
    parser.add_option("--gene", "-g", dest="gene", action="append", help="A specific gene you'd like try",
                      metavar="GENENAME")
    parser.add_option("--minreads", dest="minreads",
                      help="minimum reads required for a section to start the fitting process.  Default:%default",
                      default=3, type="int", metavar="NREADS")
    # parser.add_option("--premRNA", dest="premRNA", action="store_true", help="use premRNA length cutoff, default:%default", default=False)
    parser.add_option("--poisson-cutoff", dest="poisson_cutoff", type="float",
                      help="p-value cutoff for poisson test, Default:%default", default=0.05, metavar="P")
    parser.add_option("--disable_global_cutoff", dest="use_global_cutoff", action="store_false",
                      help="disables global transcriptome level cutoff to CLIP-seq peaks, Default:On", default=True,
                      metavar="P")
    parser.add_option("--FDR", dest="FDR_alpha", type="float", default=0.05,
                      help="FDR cutoff for significant height estimation, default=%default")
    # parser.add_option("--threshold-method", dest="method", default="binomial", help="Method used for determining height threshold, Can use default=random or binomial")
    parser.add_option("--binomial", dest="binom", type="float", default=0.05,
                      help="Alpha significance threshold for using Binomial distribution for determining height threshold, default=%default")
    parser.add_option("--threshold", dest="threshold", type="int", default=None,
                      help="Skip FDR calculation and set a threshold yourself")
    parser.add_option("--maxgenes", dest="maxgenes", default=None, type="int",
                      help="stop computation after this many genes, for testing", metavar="NGENES")
    parser.add_option("--processors", dest="np", default="autodetect",
                      help="Number of processors to use. Default: All processors on machine", type="str", metavar="NP")
    # parser.add_option("--superlocal", action="store_true", dest="SloP", default=False, help="Use super-local p-values, counting reads in a 1KB window around peaks")
    parser.add_option("--plot", "-p", dest="plotit", action="store_true", help="make figures of the fits",
                      default=False)
    parser.add_option("--verbose", "-v", dest="verbose", action="store_true", default=False)  # TODO unused ?
    parser.add_option("--quiet", "-q", dest="quiet", action="store_true", default=False,
                      help="suppress notifications")  # TODO unused ?
    parser.add_option("--save-pickle", dest="save_pickle", default=False, action="store_true",
                      help="Save a pickle file containing the analysis")
    parser.add_option("--debug", dest="debug", default=False, action="store_true",
                      help="disables multipcoressing in order to get proper error tracebacks")
    # parser.add_option("--max_width", dest="max_width", type="int", default=75, help="Defines max width for classic algorithm, default: %default")
    # parser.add_option("--min_width", dest="min_width", type="int", default=50, help="Defines min width for classic algorithm, default: %default")
    parser.add_option("--max_gap", dest="max_gap", type="int", default=15,
                      help="defines maximum gap between reads before calling a region a new section, default: %default")
    # parser.add_option("--bonferroni", dest="bonferroni_correct",action="store_true", default=False, help="Perform Bonferroni on data before filtering")
    # parser.add_option("--algorithm", dest="algorithm",default="spline", help="Defines algorithm to run, currently spline, classic, gaussian")
    # parser.add_option("--reverse_strand", dest="reverse_strand",default=False, action="store_true", help="adds option to reverse strand")
    parser.add_option("--timeout", dest="timeout", default=None, type=int,
                      help="adds timeout (in seconds) to genes that take too long (useful for debugging only, or if you don't care about higly expressed genes)")
    return parser

def override_options(options):
    ''' override some options (gtfFile, premRNA, method, SloP, bonferroni, algorithm) are removed and are added here
    :param options: OptionParser object
    :return: OptionParser object
    :rtype: OptionParser object
    '''
    if options.plotit:
        options.debug = True

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

    return options

def call_main():
    """
    Command line interface
    :rtype: None
    """
    parser = option_parser()
    (options, args) = parser.parse_args()

    ##########################################################################
    # logfile = options.outfileF + '.log'
    # logging.basicConfig(filename=logfile ,level=logging.DEBUG)
    # logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
    ##########################################################################

    # override some old options
    options = override_options(options)

    print(options)

    # enforces required usage
    # if not (options.bam and ((options.species) or (options.gtfFile))):
    if not (options.bam and options.species):
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
