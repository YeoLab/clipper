'''
Created on Jul 25, 2012
@author: mlovci
@author: gabrielp

Refactored on June 29, 2020
@author: algaebrown

'''
import logging
import numpy as np
import pandas as pd
from scipy import stats
from clipper.src.utils import write_peak_bedtool_string


############################################################
### This files is for filtering peak by No. reads ##########
############################################################

def count_transcriptome_reads(peaks_dicts):
    """
    Counts number of reads in the entire transcriptome
    :param peaks_dicts: list of peak_dict (dict) containing peaks called in each transcript
    :return: int, the number of reads in the transcriptome
    :rtype: int
    """
    ################################################
    # logging.info(" number of reads per gene result")
    ################################################
    # count total number of reads in transcriptiome
    transcriptome_reads = 0

    for peaks_dict_no, peaks_dict in enumerate(peaks_dicts):
        if peaks_dict is not None:
            if int(peaks_dict['nreads']) > 10:
                logging.info("     gene_result_no: %s , number_of_reads: %d" % (peaks_dict_no, peaks_dict['nreads']))
            transcriptome_reads += peaks_dict['nreads']

    return transcriptome_reads


def count_transcriptome_length(results):
    """
    :param results: list of peak_dicts
    :return: total transcriptome_length
    :rtype: int
    """
    transcriptome_length = 0

    for gene_result in results:
        if gene_result is not None:
            transcriptome_length += int(gene_result['loc'].attrs['effective_length'])

    return transcriptome_length

def dictify(some_named_tuple):
    """
    Convert Named Tuple to dictionary
    :param some_named_tuple: collections.namedtuple()
    :return:dict,
    """
    return dict((s, getattr(some_named_tuple, s)) for s in some_named_tuple._fields)


def make_peak_dataframe(peaks_dicts):
    '''
    make peak_dicts['cluster'](peaks) into dataframe
    :param peaks_dicts: list of peak_dict
    :return: pd.DataFrame containing all peaks
    '''
    peaks_list = [dictify(cluster) for peaks_dict in peaks_dicts for cluster in peaks_dict['clusters']]
    peaks_dataframe = pd.DataFrame(peaks_list)
    return peaks_dataframe


def poissonP(reads_in_gene, reads_in_peak, gene_length, peak_length):
    """

    scipy.stats.poisson.cdf
    compute the p-value for a peak of peak_length length with reads_in_peak reads,
    given that the read is gene_length long and has reads_in_gene reads

    If there are fewer than 3 reads expected to fall in the region, assume there's 3 reads
    expected...

    Paramaters
    ----------
    reads_in_gene: Integer representing number of reads in gene
    reads_in_peak: Integer representing the number of reads in a specific peak
    gene_length: Integer representing length of gene
    peak_length: Integer representing length of peak

    Returns double, the p-value that the peak is significant
    If calculation fails returns 1

    """

    try:
        # lam is estimate of the lambda value
        # poission takes a value and the lambda
        # this is average number of reads per single
        # site in the gene, but
        # a peak is not a single site, so it the average number
        # gets multipled by the peak
        # length as an estimator of the mean

        lam = 1 + ((float(reads_in_gene) / float(gene_length)) * float(peak_length))  # expect at least one read.

        cum_p = stats.poisson.sf(reads_in_peak, lam)

        return cum_p

    except Exception as error:
        logging.error("Poisson cutoff failled %s " % (error))
        return 1

def transcriptome_poissonP(cluster):
    """
    Applied to peaks_dataframe with params in there
    used to determine P value if the number of reads in a peak is significant by comparing *whole transcriptome*
    :param cluster: pd.Series.
    :return: function, PoissionP(transcriptome_reads, no.read in peak, transcriptome_size, cluster size)
    :rtype: function"""
    return poissonP(cluster.transcriptome_reads,
                    cluster.number_reads_in_peak,
                    cluster.transcriptome_size,
                    cluster['size'])


def transcript_poissonP(cluster):
    """
    Applied to peaks_dataframe with params in there
    return P value indicating if the number of reads in peak is significant by comparing whole transcript
    :rtype: function
    :param cluster: pd.Series in peaks_dataframe
    :return: function PoissionP(transcript_read, no.read in peak, transcript effective length, cluster size)
    """
    return poissonP(cluster.nreads_in_gene,
                    cluster.number_reads_in_peak,
                    cluster.effective_length,
                    cluster['size'])


def superlocal_poissonP(cluster):
    """
    Applied to peaks_dataframe with params in there
    return P value indicating if the number of reads in peak is significant by comparing local region (+/- 500 in section)
    :rtype: function
    :param cluster: pd.Series in peaks_dataframe
    :return: function PoissionP(local_read, no.read in peak, local length, cluster size)
    """
    return poissonP(cluster.area_reads,
                    cluster.number_reads_in_peak,
                    cluster.area_size,
                    cluster['size'])


def bh_correct(df):
    """
    Bonferroni correction on peak_dataframe
    :return: pd.DataFrame,  dataframe with adjusted p-value
    :param df: pd.DataFrame with adjusted p-value ['padj']
    """
    df = df.sort_values("final_p_value")
    df['sort_rank'] = np.arange(1, len(df) + 1)
    df['bh_corrected'] = df.apply(lambda x: min(((len(df) / x.sort_rank) * x.final_p_value), 1), axis=1)
    df['padj'] = df.sort_values("final_p_value", ascending=False).bh_corrected.cummin()
    return df.sort_index()

def filter_peaks_dicts(peaks_dicts, poisson_cutoff, transcriptome_size,
                       transcriptome_reads, use_global_cutoff,
                       bonferroni_correct=True, algorithm="spline", superlocal=False, min_width=50,
                       bypassfiltering=False):
    """

    Takes a list of peaks_dicts, filters them based on number of reads in peak, Bonferroni correction

    :param peaks_dicts: list, list of peaks_dicts,
    :param poisson_cutoff: float, p-value cutoff for Poission (number of reads)
    :param transcriptome_size: int, No.genes in transcriptome
    :param transcriptome_reads: int, total no. reads in transcriptome
    :param use_global_cutoff: bool
    :param bonferroni_correct: bool, default True (OVERRIDED)
    :param algorithm: str, default "spline" (OVERRIDED)
    :param superlocal: bool, default TRUE, +/-500 region in each section
    :param min_width: int, min peak width
    :param bypassfiltering: bool, default FALSE, if TRUE, return unfiltered result
    :return: str, Peak in the format of NarrowPeak
    :rtype: str

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
    peaks_dataframe['transcriptome_poisson_p'] = peaks_dataframe.apply(transcriptome_poissonP,
                                                                       axis=1) if use_global_cutoff else np.nan
    peaks_dataframe['transcript_poisson_p'] = peaks_dataframe.apply(transcript_poissonP, axis=1)
    peaks_dataframe['superlocal_poisson_p'] = peaks_dataframe.apply(superlocal_poissonP,
                                                                    axis=1) if superlocal else np.nan

    if algorithm == "classic":
        # TODO this never happens!
        peaks_dataframe['final_p_value'] = peaks_dataframe[['transcript_poisson_p', 'superlocal_poisson_p']].max(axis=1)
    else:
        # TODO this is always the case!
        peaks_dataframe['final_p_value'] = peaks_dataframe[['transcript_poisson_p', 'superlocal_poisson_p']].min(axis=1)

    if bonferroni_correct:
        # TODO always True!
        peaks_dataframe = bh_correct(peaks_dataframe)
        # peaks_dataframe['final_p_value'] = (peaks_dataframe['final_p_value'] * total_clusters)    # TODO still needed?

    if bypassfiltering:
        filtered_peaks_dataframe = peaks_dataframe
    else:
        # This is a bug I should fix, padj isn't getting printed, the uncorreded p-value is
        filtered_peaks_dataframe = peaks_dataframe[peaks_dataframe['padj'] < poisson_cutoff]

    #########################################################################
    # logging.info(" total clusters AFTER filtering : %s" % len(filtered_peaks_dataframe))
    #########################################################################

    filtered_peak_bedtool_strings = filtered_peaks_dataframe.apply(write_peak_bedtool_string, axis=1).values
    filtered_peak_bedtool_tsv = "\n".join(filtered_peak_bedtool_strings) + "\n"
    return filtered_peak_bedtool_tsv
