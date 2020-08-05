'''
Created on Jul 25, 2012
@author: mlovci
@author: gabrielp

Refactored on June 29, 2020
@author: algaebrown

'''

############################################################
####### This files contains many helper function ###########
############################################################
import os
import logging
from subprocess import call
import clipper
import pybedtools

def check_for_index(bamfile):
    """
    Checks to make sure a BAM file has an index, if the index does not exist it is created
    Usage undefined if file does not exist (check is made earlier in program)
    bamfile - a path to a bam file
    :param bamfile: path to bamfile
    :return: None
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

####################### LENGTHS #####################################

def build_transcript_data_gtf_as_structure(species, pre_mrna):
    """
    calculate effective length for each transcript from pre-created gtf file in clipper/data
    Returns Bedtool containing effective length
    :param species: (str) genome name
    :param pre_mrna: (bool) if true uses pre-mRNA length instead of mRNA length
    :return: (pybedtools.Bedtool)
    :rtype: pybedtools.BedTool

    """
    bedtool_intervals = []
    x = clipper.data_file(species + ".AS.STRUCTURE.COMPILED.gff")
    gtf_file = pybedtools.BedTool(x)
    for gene in gtf_file:
        effective_length = gene.attrs['premrna_length'] if pre_mrna else gene.attrs['mrna_length']
        attrs = "gene_id=%s;" % (gene.attrs['gene_id'])
        if "transcript_ids" in gene.attrs:
            attrs += "transcript_ids=%s;" % (gene.attrs['transcript_ids'])
        attrs += "effective_length=%s" % (str(effective_length))

        # add to bedtool_intervals
        to_string = map(str, [gene['chrom'],"AS_STRUCTURE","mRNA",str(gene.start + 1),str(gene.stop + 1),"0",gene['strand'],".",attrs]) # map object
        bedtool_intervals.append(pybedtools.create_interval_from_list(list(to_string)))


    return pybedtools.BedTool(bedtool_intervals)

def get_exon_bed(species):
    short_species = species.split("_")[0]
    return os.path.join(clipper.data_dir(), "regions", "%s_%s.bed" % (short_species, "exons"))

def write_peak_bedtool_string(cluster):
    """
    Format Peak into NarrowBed format
    :param cluster: Peak object
    :return: str, format to NarrowBed format with tab-delimited, [chrom, start, stop, name, pval, strand, thick_start, thick_stop]
    """
    cluster_info_list = [
        cluster.chrom,
        cluster.genomic_start,
        cluster.genomic_stop,
        cluster.gene_name + "_" + str(cluster.peak_number) + "_" + str(cluster.number_reads_in_peak),
        cluster.final_p_value,
        cluster.strand,
        cluster.thick_start,
        cluster.thick_stop,
    ]
    cluster_bedtool_string = "\t".join([str(info) for info in cluster_info_list])
    return cluster_bedtool_string

################################ DEPRECATED UNUSED ###############################################
def build_geneinfo(bed):
    """
    Loads bed file into a dictionary with the key being the name and a string being the value

    :param bed: (str) file path to bed file
    :return: gene_info (dict) key: name position of the bed file; values: ordered bed file
    """

    # opens bed file, either zipped or unzipped
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
    :rtype: dict

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

    Creates a dictionary containing all information needed to perform peak calling calculations
    for a single species

    :parameters
    -----------
    species: string currently not used
    chrs: list specifying all the chromosomes in a given species
    bed: path to a bed file that contains information on genes (custom file *STRUCTURE_genes.BED.gz)
    mrna: path to a file that contains mRNA lengths (custom CSV file contains gene names followed by gene lengths)
    premrna: path to a file that contains pre-mRNA lengths (custom CSV file contains gene names followed by gene lengths_

    Returns dict of all items passed to it

    TODO:  Add checking to verify that file are actually passed
    :rtype: dict
    """
    par = dict()

    # this is non-pythonic, should just combine all lists
    # expand sublists
    par["chrs"] = [item for sublist in chrs for item in sublist]
    par["gene_bed"] = bed
    par["mRNA"] = mrna
    par["premRNA"] = premrna
    return par


def build_transcript_data(species, gene_bed, gene_mrna, gene_pre_mrna, pre_mrna):
    """

    Generates transcript data structures to call peaks on
    Allows for either predefined files (from the data directory)
    or custom files

    :param species: (str) the species genome to run on
    :param gene_bed: an abribtary bed file of locations to search for peaks (should be gene locations)
    :param gene_mrna: the effective length of the mrna of a gene (unmappable regions removed)
    :param gene_pre_mrna: the effective length of the pre-mrna (unmappable regions removed)
    :param pre_mrna: (bool) flag True indicates use pre-mRNA lengths instead of mRNA lengths
    :return:
    :rtype: pyBedtool.BedTool

    """

    # error checking

    acceptable_species = get_acceptable_species()
    if (species is None and
            gene_bed is None and
            (gene_mrna is None or gene_pre_mrna is None)):
        raise ValueError("You must set either \"species\" or \"geneBed\"+\"geneMRNA\"+\"genePREMRNA\"")

    if species is not None and gene_bed is not None:
        raise ValueError("You shouldn't set both geneBed and species, defaults exist for %s" % (acceptable_species))

    # Now actually assign values
    if species is not None:
        try:
            gene_bed = clipper.data_file(species + ".AS.STRUCTURE_genes.BED.gz")
            gene_mrna = clipper.data_file(species + ".AS.STRUCTURE_mRNA.lengths")
            gene_pre_mrna = clipper.data_file(species + ".AS.STRUCTURE_premRNA.lengths")

        except ValueError:
            raise ValueError(
                "Defaults don't exist for your species: %s. Please choose from: %s or supply \"geneBed\"+\"geneMRNA\"+\"genePREMRNA\"" % (
                    species, acceptable_species))

    # Selects mRNA or preMRNA lengths
    if pre_mrna is True:
        lenfile = gene_pre_mrna
    else:
        lenfile = gene_mrna

    if lenfile is None:
        raise IOError("""didn't pass correct mRNA length file option 
                    with given length file""")

    # builds dict to do processing on,
    genes = build_geneinfo(gene_bed)
    lengths = build_lengths(lenfile)

    # this is a stopgap until it can be fully factored out, returing a gtf file of
    # genes and effective lengths, eventually this is the file we want to pass in
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
                                                              "gene_id=" + gene + "; effective_length=" + str(
                                                                  lengths[gene])]))

    return pybedtools.BedTool(gtf_list)


def build_transcript_data_gtf(gtf_file, pre_mrna):
    """
    Generates GTF file to use when calling genes
    Returns the longest gene from a group of transcripts to call peaks on (this isn't optimal
    behavior, but until we get a general as structure working its alright)


    :param gtf_file: (pybedtool.BedTool) bedtool from a standard gtf file
    :param pre_mrna: (bool) use pre_mRNA length instead of mrna
    :return: pyBedtools.BedTool
    """

    # objects for default dict, no need to test or factor out
    def default_transcript():
        return {'chrom': None, 'start': np.inf, "stop": np.NINF, "strand": None, "gene_id": None, "mRNA_length": 0}

    def default_gene():
        return {'start': 0, 'stop': 0}

    # get all transcripts, their starts, stops and mrna lengths
    transcripts = defaultdict(default_transcript)
    gtf_file = gtf_file.filter(lambda x: x[2] == 'exon').saveas()

    for interval in gtf_file:
        cur_transcript = transcripts[interval.attrs['transcript_id']]
        cur_transcript['start'] = min(cur_transcript['start'], interval.start)
        cur_transcript['stop'] = max(cur_transcript['stop'], interval.stop)
        cur_transcript['chrom'] = interval.chrom
        cur_transcript['strand'] = interval.strand
        cur_transcript['gene_id'] = interval.attrs['gene_id']
        cur_transcript['mRNA_length'] += interval.length
        cur_transcript['transcript_id'] = interval.attrs['transcript_id']

    # get the longest transcript from each gene group
    longest_genes = defaultdict(default_gene)
    for transcript_name, transcript in transcripts.items():
        cur_gene = transcript['gene_id']
        foo = longest_genes[cur_gene]
        best_length = longest_genes[cur_gene]['stop'] - longest_genes[cur_gene]['start']
        cur_length = transcript['stop'] - transcript['start']
        if best_length < cur_length:
            longest_genes[cur_gene] = transcript

    # convert back into a gtf file
    bedtool_intervals = []
    for gene in longest_genes.values():
        effective_length = gene['stop'] - gene['start'] if pre_mrna else gene['mRNA_length']
        bedtool_intervals.append(pybedtools.create_interval_from_list([gene['chrom'],
                                                                       "AS_STRUCTURE",
                                                                       "mRNA",
                                                                       str(gene['start']),
                                                                       str(gene['stop']),
                                                                       "0",
                                                                       gene['strand'],
                                                                       ".",
                                                                       "gene_id=" + gene[
                                                                           'gene_id'] + "; transcript_id=" + gene[
                                                                           'transcript_id'] + "; effective_length=" + str(
                                                                           effective_length)]))
    return pybedtools.BedTool(bedtool_intervals)


##############################
# only used along with --debug, was used by github code, but no longer necessary
##############################
def func_star(variables):
    """ covert f([1,2]) to f(1,2) """
    return call_peaks(*variables)


def get_acceptable_species():
    """

    Finds all species in data directory
    :return: acceptable_species (set): string of available genome

    """
    acceptable_species = set([])
    for fn in os.listdir(clipper.data_dir()):
        fn = fn.split(".")[0]
        if fn == "__init__":
            continue
        acceptable_species.add(fn)
    return acceptable_species
