

"""

Analizes CLIP data given a bed file and a bam file

Michael Lovci and Gabriel Pratt

"""
from collections import Counter, OrderedDict, defaultdict
from itertools import izip_longest
from optparse import OptionParser
import os
import pickle
import subprocess

from bx.bbi.bigwig_file import BigWigFile
import HTSeq
import gffutils
import numpy as np
import pandas as pd
import pybedtools
from sklearn.cluster import KMeans
from AS_Structure_tools import parse_AS_STRUCTURE_dict
import clipper
from clipper.src import CLIP_analysis_display
from clipper.src.kmerdiff import kmer_diff
from clipper.src.get_genomic_regions import GenomicFeatures
from gscripts.general.pybedtools_helpers import small_peaks, convert_to_mRNA_position

import matplotlib as mpl
from matplotlib import rc

mpl.rcParams['svg.fonttype'] = 'none'
rc('text', usetex=False)
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})


def name_to_chrom(interval):
    interval.chrom = interval.name
    return interval


def intersection(a, b=None):
    
    """
    
    A : bedtool
    B : bedtool
    Returns A - B and A intersect B 
    A with b and returns everything in a but not b and everything in a but... ???
    
    """
    
    overlapping = a.intersect(b, wa=True, u=True, stream=True).saveas()
    remaining = a.intersect(b, wa=True, v=True, stream=True).saveas()

    return remaining, overlapping


def count_genomic_types(as_structure_dict):
    
    """
    
    Counts number of exons for each splicing type
    CE SE ect...
    
    as_structure_dict - as structure dict format as defined in parse_AS_STRUCTURE_dict
    
    returns dict of counts
    
    """
    
    genomic_types = Counter()
    for gene in as_structure_dict.values():
      
        for exon_type in gene['types'].values():
            genomic_types[exon_type] += 1
    
    return genomic_types


def count_genomic_region_sizes(regions, species="hg19"):
    
    """
    
    Counts the genomic region sizes for the specified regions dir
    (should be changed to GTF)
    
    """
    
    genomic_region_sizes = {}
    #TODO update this to work of GFF file, because something isn't matching up...
    for region in regions:
        region_tool = pybedtools.BedTool(os.path.join(clipper.data_dir(), "regions",  species + "_" + region + ".bed"))
        genomic_region_sizes[region] = region_tool.total_coverage()
    return genomic_region_sizes


def get_sizes(cluster_regions):
    """

    :param cluster_regions: Dict of clusters generatd by assign_to_regions
    :return: dict of {region: size of region}

    """

    return {region: bedtools['real'].total_coverage() for region, bedtools in cluster_regions.items()}


def make_fasta_files_from_regions(cluster_regions, clusters, fasta_dir, speciesFA):
    """
    Makes and saves fasta files given a cluster regions dict
    :param cluster_regions:
    :return:  nothing, but saves lots of fasta files to a dir
    """

    for region in cluster_regions:
        region_fa = fa_file(clusters, region=region, directory=fasta_dir, type="real")
        cluster_regions[region]['real'].sequence(fi=speciesFA, s=True).save_seqs(region_fa)

        if "rand" not in cluster_regions[region]:
            continue

        rand_list = []
        for n_rand in cluster_regions[region]['rand']:
            rand_list.append(cluster_regions[region]['rand'][n_rand].sequence(fi=speciesFA, s=True))

        write_seqs(fa_file(clusters, region=region, directory=fasta_dir, type="random"), rand_list)


def save_bedtools(cluster_regions, clusters, assigned_dir):
    """
    Given cluster regions file saves all bedtools sanely and returns result
    :param cluster_regions:
    :return:
    """
    for region in cluster_regions:
        output_file = "%s.%s.real.BED" % (clusters, region)
        cluster_regions[region]['real'] = cluster_regions[region]['real'].sort().saveas(os.path.join(assigned_dir, output_file))

        if "rand" not in cluster_regions[region]:
            continue

        for n_rand in cluster_regions[region]['rand']:
            output_file = "%s.%s.rand.%s.BED" % (clusters, region, n_rand)
            cluster_regions[region]['rand'][n_rand] = cluster_regions[region]['rand'][n_rand].sort().saveas(os.path.join(assigned_dir, output_file))

    return cluster_regions


def fix_strand(interval):
    strands = list(set(interval.strand.split(",")))
    if len(strands) > 1:
        raise NameError("Both strands are present, something went wrong during the merge")
    interval.strand = strands[0]
    return interval


def fix_shuffled_strand(shuffled_tool, regions_file):
    """
    Fixes strand of shuffled tool if there was an include file used

    Chooses the first overlapping thing as the correct strand, if there is anti-sense transcription
    this will fail
    shuffled_tool: a sorted shuffled bedtool
    shuffled_file: incl file that did the shuffling
    """

    regions_tool = pybedtools.BedTool(regions_file)
    intersected = shuffled_tool.intersect(regions_file, wao=True, stream=True, sorted=True).saveas()

    #Don't think about refactoring this because of leaky files in pybedtools
    shuffled_tool_field_count = shuffled_tool.field_count()
    regions_tool_field_count = regions_tool.field_count()

    total_header_size = shuffled_tool_field_count + regions_tool_field_count + 1

    intersected = intersected.to_dataframe(names=range(total_header_size))
    intersected = intersected.groupby([0,1,2,3]).first()
    if regions_tool.file_type == "bed":
        intersected[5] = intersected[shuffled_tool_field_count + 5]
    if regions_tool.file_type == "gff":
        intersected[5] = intersected[shuffled_tool_field_count + 6]

    values = intersected.apply(lambda x: "\t".join([str(item) for item in list(x.name) + list(x[:shuffled_tool_field_count - 4])]), axis=1).values
    new_bedtool = pybedtools.BedTool("\n".join(values), from_string=True)
    return new_bedtool

def assign_to_regions(tool, clusters, regions, assigned_dir=".", species="hg19", nrand=3):
    
    """
    
    Assigns each cluster to a genic region
    finally saves all generated bed and fasta files for future analysis...
    
    Needs to be refactored desperatly
    
    tool - a bed tool (each line represnting a cluster)
    clusters - name of cluster file
    regions_dir - the directory that has genomic regions already processed
    regions - dict [str] regions to process, not used now but should be after refactoring 
    assigned_dir - location to save files in
    species - str species to segment
    nrand - int number offsets times to shuffle for null hypothesis
    getseq - boolean gets the full sequence to store
    
    shuffling for background, should still be factored out
    
    """

    bedtracks = {}

    #TODO refactor this to use get_genomic_regions
    for region in regions:
        bedtracks[region] = pybedtools.BedTool(os.path.join(clipper.data_dir(), "regions", "%s_%s.bed" % (species,
                                                                                                          region)))
              
    #creates the basics of bed dict
    bed_dict = {'all': {'rand': {}}}

    #can't append None to string so hack a null bed tool
    for n in range(nrand):
        bed_dict['all']['rand'][n] = pybedtools.BedTool("", from_string=True)

    bed_dict['all']['real'] = pybedtools.BedTool("", from_string=True)
    
    offsets = get_offsets_bed12(tool)
    if len(list(tool[0])) < 8:
        tool = tool.sort().merge(s=True, c="4,5,6", o="collapse,collapse,collapse", stream=True).each(fix_strand).saveas()
    else:
        tool = tool.sort().merge(s=True, c="4,5,6,7,8", o="collapse,collapse,collapse,min,min", stream=True).each(fix_strand).saveas()
    remaining_clusters = adjust_offsets(tool, offsets)
    
    print "There are a total %d clusters I'll examine" % (len(tool))
    for region in regions:
        remaining_clusters, overlapping = intersection(remaining_clusters, b=bedtracks[region])

        #if for some reason there isn't a peak in the region skip it
        if len(overlapping) == 0:
            print "ignoring %s " % region
            continue
        
        #sets up bed dict for this region
        bed_dict[region] = {'real': overlapping.sort(stream=True).saveas(),
                            'rand': {}}

        no_overlapping_count = len(remaining_clusters)
        overlapping_count = len(bed_dict[region]['real'])
        print "For region: %s found %d that overlap and %d that don't" % (region,
                                                                          overlapping_count,
                                                                          no_overlapping_count)

        bed_dict['all']['real'] = bed_dict['all']['real'].cat(bed_dict[region]['real'], stream=True, postmerge=False)

        #saves offsets so after shuffling the offsets can be readjusted
        offset_dict = get_offsets_bed12(bed_dict[region]['real'])
        for i in range(nrand):
            random_intervals = bed_dict[region]['real'].shuffle(genome=species, incl=bedtracks[region].fn).sort()
            random_intervals = fix_shuffled_strand(random_intervals, bedtracks[region].fn)
            random_intervals = adjust_offsets(random_intervals, offset_dict)
            bed_dict[region]['rand'][i] = random_intervals.saveas()

            bed_dict['all']['rand'][i] = bed_dict['all']['rand'][i].cat(bed_dict[region]['rand'][i], stream=True, postmerge=False)

        #if there are no more clusters to assign stop trying
        if no_overlapping_count == 0:
            break

    print "After assigning %d un-categorized regions" % len(remaining_clusters)

    if len(remaining_clusters) > 0:
        bed_dict['uncatagorized'] = {'real': remaining_clusters.sort(stream=True).saveas()}

    bed_dict = save_bedtools(bed_dict, clusters, assigned_dir)
    return bed_dict


def get_offsets_bed12(bedtool):
    
    """
    
    Gets offsets for each cluster in a bed12 file (CLIPper formatted bed file)
    Offsets are the difference between the wide peak and the narrow (predicted binding site)
    This will break if we analize different types of data
    tool - bedtool
    
    returns offset for each cluster
    
    """

    try:
        offset_dict = {}
        for interval in bedtool:

            if interval.strand == "+":
                offset = int(interval[6]) - int(interval[1])
            else:
                offset = int(interval[2]) - int(interval[7])

            offset_dict[interval.name] = offset

        return offset_dict
    except:
        print "Not Valid bed12 file, continuing processing, some things may be strange"
        return None


def adjust_offsets(bedtool, offsets=None):

    """

    For finding motiff position relative to center of peaks
    Handles merging overlapping peaks, merge > bed12 -> bed6
    picks first offset from merged peak and assigns tqhat to
    new peaks

    Input:
    tool - a bed tool (bed12) object
    offset - dict of key:peak, value:int
    Adjusts the offsets for each transcript in a bedtool


    """

    if offsets is None:
        return bedtool

    clusters = []
    for interval in bedtool:
        #the ; represents two merged locations
        if "," in interval.name and interval.name not in offsets:
            offset = offsets[interval.name.split(",")[0]]
        else:
            offset = offsets[interval.name]

        if interval.strand == "+":
            thick_start = interval.start + offset
            thick_stop = thick_start + 4
        else:
            thick_stop = interval.stop - offset
            thick_start = thick_stop - 4

        clusters.append("\t".join([str(x) for x in (interval.chrom,
                                                    interval.start,
                                                    interval.stop,
                                                    interval.name,
                                                    interval.score,
                                                    interval.strand,
                                                    thick_start,
                                                    thick_stop)]))

    return pybedtools.BedTool("\n".join(clusters), from_string=True)


def generate_region_dict(bedtool):
    region_dict = defaultdict(list)
    
    for interval in bedtool:
        region_dict[interval.name].append(interval)
    
    for gene in region_dict.keys():
        if region_dict[gene][0].strand == "-":
            region_dict[gene].reverse()
            
    return region_dict


def convert_to_transcript(feature_dict):
    
    """
    
    Converts bed file to be based only off transcripts (transcripts are defined by chromosome)
    
    returns modified dict
        
    """
    return {name: bedtool.each(name_to_chrom).saveas() for name, bedtool in feature_dict.items()}


def convert_to_mrna(feature_dict, exon_dict):
    
    """
    
    converts transcript dict into mRNA locations given a dict of exons
    
    feature_dict - generated from convert to transcript, dict of bedfiles 
    exon_dict - dict of { genes : [exons] }
    
    """

    return {name: bedtool.each(convert_to_mRNA_position, exon_dict).filter(lambda x: x.chrom != "none").saveas() for name, bedtool in feature_dict.items()}


def invert_neg(interval):
    interval[-1] = str(int(interval[-1]) * -1)
    return interval


def get_feature_distances(bedtool, regions, features):
    """

    :param bedtool: bedtools of peaks, both  (only works with clipper defined peaks)
    :param regions: dict of regions, used to convert genomic coordinates to mrna coordinates
    :param features: dict of features to find distance from
    :return:
    """

    exon_dict = generate_region_dict(regions['exons'])
    
    bed_center = bedtool.each(small_peaks).saveas()
    beds_center_transcripts = bed_center.each(name_to_chrom).saveas()
    beds_center_transcripts_mrna = beds_center_transcripts.each(convert_to_mRNA_position, exon_dict).filter(lambda x: x.chrom != "none").saveas()

    features_transcript = convert_to_transcript(features)
    features_mrna = convert_to_mrna(features_transcript, exon_dict)
    
    #for pre-mrna
    features_transcript_closest = {}
    for name, feature in features_transcript.items():
        try:
            features_transcript_closest[name] = {"dist": beds_center_transcripts.closest(feature, s=True, D="b", t="first", id=True).filter(lambda x: x[-2] != ".").saveas()}
        except:
            features_transcript_closest[name] = {"dist": None}
    #for mrna
    features_mrna_closest = {} 
    for name, feature in features_mrna.items():
        try:
            features_mrna_closest[name] = {"dist": beds_center_transcripts_mrna.closest(feature, s=True, D="ref", t="first").filter(lambda x: x[-2] != ".").each(invert_neg).saveas()}
        except:
            features_mrna_closest[name] = {"dist": None}

    return features_transcript_closest, features_mrna_closest


def get_region_distributions(bedtool, regions):
    """
    Gets location of each peak across different regions
    :param bedtool: pybedtool of peaks
    :param regions: dict of regions to
    :return:
    """
    distributions = {}
    for name, region in regions.items():
        region_dict = generate_region_dict(region)
        distributions[name] = get_distributions(bedtool, region_dict)
    return distributions


def get_distributions(bedtool, region_dict):

    """

    Gets distributions from RNA_position function

    bedtool - clipper bed file to generate distributions from
    region_dict - generate_region_dict dict defining regions

    """

    exon_distributions = []
    total_distributions = []
    num_errors = []
    num_missed = []
    bound_in_regions = []
    genes = []
    for interval in bedtool:
        try:
            #will need to redefine this to use intervals
            exon, total, bound_in_region, gene = RNA_position_interval(interval, region_dict)

            if total is not None:
                total_distributions.append(total)
                exon_distributions.append(exon)
                bound_in_regions.append(bound_in_region)
                genes.append(gene)
            else:
                num_missed.append(interval)
        except Exception as e:
            print e
            num_errors.append(interval)

    return {'individual': exon_distributions, 'total': total_distributions, "gene_ids": genes,
            "region_numbers": bound_in_regions, 'errors': num_errors, 'missed': num_missed}


def RNA_position_interval(interval, location_dict):

    """

    makes mrna and pre-mrna peak_center figure
    interval - single interval

    location_dict = dict{gene_id : {strand : "+/-", regions : list((start,stop)
    as_structure_dict - from build AS structure dict

    will return distribution across entire region + just across specific region
    Might be able to use my ribo-seq stuff for genic -> transcriptomic location conversion

    this is based off as structure, which provides sequences ordered with first exon being the first exon on the gene,
    not first in the chromosome (like gff does) THIS WILL NOT WORK WITH RAW GFF DATA

    """

    #think about turning the location_dict into a gff file
    #gets thickstart and stop
    peak_center = (int(interval[6]) + int(interval[7])) / 2

    try:
        gene = interval.name.split("_")[0]
    except:
        #takes first gene if there are multiple overlapping
        gene = interval.name.split(";")[0].split("_")[0]

    if gene not in location_dict:
        raise KeyError(gene + """ not in current region dict, ignoring cluster (not to worry, this error gets thrown if the peak isn't in the region being looked at, or you have your annotations wrong, double check you're using a supported genome or make a regions file yourself""")

    if not interval.strand == location_dict[gene][0].strand:
        raise ValueError("strands not the same, there is some issue with gene annotations")

    total_length = float(sum(region.length for region in location_dict[gene]))

    running_length = 0
    for region_number, region in enumerate(location_dict[gene]):
        length = float(region.length)

        if int(region.start) <= peak_center <= int(region.stop):
            if interval.strand == "+":
                total_location = running_length + (peak_center - region.start)
                total_fraction = np.round((total_location / total_length), 3)

                individual_fraction = (peak_center - region.start) / length

            elif interval.strand == "-":
                total_location = running_length + (region.stop - peak_center)
                total_fraction = total_location / total_length
                individual_fraction = (region.stop - peak_center) / length
            else:
                raise ValueError("Strand not correct strand is %s" % interval.strand)

            #probably not nessessary
            if total_fraction < 0 or total_fraction > 1:
                raise ValueError("total_fraction is bad: %f, gene %s, total_length: %s, total_location: %s" % (total_fraction,
                                                                                                               gene,
                                                                                                              total_length,
                                                                                                               total_location))
            return individual_fraction, total_fraction, region_number + 1, gene

        running_length += length

    #clusters fall outside of regions integrated
    return None, None


def get_closest_exon_types(bedtool, as_structure_dict):
    """

    For each peak gets type of exon nearest to peak (AS, CE ect...)
    :param bedtool: Bedtool of peaks
    :param as_structure_dict: AS Structure dict
    :return:
    """
    #get closest features to peaks
    bedtool_list = []
    for name, gene in as_structure_dict.items():
        for exon, type in zip(gene['exons'].values(), gene['types'].values()):
            start, stop = map(int, exon.split("-"))
            bedtool_list.append([gene['chrom'], start, stop, name, 0, gene['strand'], type])

    feature_tool = pybedtools.BedTool(bedtool_list)
    return Counter([interval[-1] for interval in bedtool.closest(feature_tool, s=True)])

#TODO Start small module for getting read densiites
from HTSeq import SAM_Alignment
import pysam


class Robust_BAM_Reader(HTSeq.BAM_Reader):

    def __iter__( self ):
        sf = pysam.Samfile(self.filename, "rb")
        self.record_no = 0
        for pa in sf:
            try:
                yield SAM_Alignment.from_pysam_AlignedRead( pa, sf )
            except OverflowError:
                pass
            self.record_no += 1


def get_bam_coverage(bamfile):
    """

    Given a bam file returns a properly covered htseq coverage file (this is slow)

    """
    bam = Robust_BAM_Reader(bamfile)
    coverage = HTSeq.GenomicArray("auto", typecode="i", stranded=True)
    for read in bam:
        if read.aligned:
            for cigop in read.cigar:
                if cigop.type != "M":
                    continue
                coverage[cigop.ref_iv] += 1
    return coverage


def get_bam_counts(bamfile):
    """
    Given a bam file returns a coverage file with just count of the number of reads that start at a specific location
    :param bamfile:
    :return HTSeq coverage :
    """
    bam = Robust_BAM_Reader(bamfile)
    coverage = HTSeq.GenomicArray("auto", typecode="i")
    for read in bam:
        read.iv.length = 1
        if read.aligned:
            coverage[read.iv] += 1
    return coverage


def bed_to_genomic_interval(bed):
    for interval in bed:
        yield HTSeq.GenomicInterval(interval.chrom, interval.start, interval.stop, interval.strand)


def get_densities(intervals, coverage):
    for interval in bed_to_genomic_interval(intervals):
        density = np.fromiter(coverage[interval], dtype="i")
        if interval.strand == "-":
            density = density[::-1]
        yield density


def adjust_width(interval, width=250):
    interval.start -= width
    interval.stop += width
    return interval


def cluster_peaks(bedtool, coverage, k=16):
    
    """
    
    Given a bedtool of peaks, positive and negative bigwig files
    gets read densities around peaks, normalizes them and outputs clusters and dataframe 
    
    """

    #TODO change small peaks so they can automatically fall back on center of peak if bed format isn't present
    bed_center = bedtool.each(small_peaks).saveas()
    bedtool_df = pd.DataFrame(list(get_densities(bed_center.each(adjust_width), coverage)))

    bedtool_df_mag_normalized = bedtool_df.div(bedtool_df.sum(axis=1), axis=0)
    bedtool_df_mag_normalized[bedtool_df_mag_normalized > .1] = .1
    #TODO Write Test, tested slowly in ipython notebook, but need to make some sample data to test with
    data = np.array(bedtool_df_mag_normalized)
    
    #fixes bad test case
    if len(data) < k:
        k = len(data) / 2
        
    kmeans_classifier = KMeans(k)
    classes = kmeans_classifier.fit_predict(data)
    return bedtool_df_mag_normalized, classes


def run_homer(foreground, background, k=list([5,6,7,8,9]), outloc=os.getcwd()):
    
    """
    
    runs homer with standard args
    output location is saved
    
    foreground - str, location of fasta file for the foreground distribution
    background - str, location of fasta file for the background distribution
    k - different k-mers to examine
    outloc - directory to output homer results 

    --make optional make work off locations and not fasta files 
    
    """
    #findMotifs.pl clusters.fa fasta outloc -nofacts p 4 -rna -S 10 -len 5,6,7,8,9 -noconvert -nogo -fasta background.fa
    #converts k to a string for use in subprocess
    k = ",".join([str(x) for x in k])
    print "starting Homer"
    
    try:
        with open(os.devnull, 'w') as fnull:
            result = subprocess.call(["findMotifs.pl", foreground, "fasta", 
                                  outloc, "-nofacts", "-p", "4", "-rna", "-S", "20",
                                   "-len", k, "-noconvert", "-nogo", "-fasta", background], shell=False, stdout=fnull)
            
        print "Homer Finished, output here: %s" % outloc
    except OSError:
        print "Homer not installed, ignoring motif generation, install homer for this to work"  
        raise   


def write_seqs(outfile, bedtool_list):
    
    """
    
    outputs bedtools file to another file 
    
    Combines all sequences because its the punitive background that gets used by other analyses
    """
    
    #TODO refactor to just the bedtools save function
    with open(outfile, "w") as fn:
        for bedtool in bedtool_list:
            with open(bedtool.seqfn) as fasta:
                fn.write(fasta.read())


def bedlengths(tool):
    
    """
    
    returns lengths of all lines in a bedtool
    
    """
    
    x = list()
    for line in tool:
        x.append(line.length)
    return x


def chop(chr, start,end, wsize=5):
    
    """
    
    writes some sort of phastcons thing..., not quite sure
    
    For circos plotting, not used add back if we want more output later, ignore for now
    
    """
    
    file = open("phastcons.txt", 'w')
    i = start
    while i < end:
        x = pybedtools.Interval(chr, i, (i+wsize))
        p = get_mean_phastcons(x, species="hg19")
        file.write("\t".join(map(str, [chr, i, (i+wsize-1), p])) + "\n")
        i += wsize
    file.close()


def get_mean_phastcons(bedtool, phastcons_location, sample_size = 1000):
    
    """
    
    Get means phastcons scores for all intervals in a bed tool
    bedtool - bedtool to extract data from
    phastcons_location - location of phastcons file
    
    """
    
    with open(phastcons_location) as bw_file:
        bw = BigWigFile(bw_file)
    
        data = []
        
        for bedline in bedtool.random_subset(min(len(bedtool), sample_size)):
            conservation_values = bw.get_as_array(bedline.chrom, bedline.start, bedline.stop)
            try:
                if len(conservation_values) > 0:
                    mean_phastcons = np.mean(conservation_values)
                else:
                    mean_phastcons = 0
                data.append(mean_phastcons)
            except TypeError:
                pass
    return data


def fa_file(filename, region=None, directory=None, type="real"):
    
    """
    
    Formats fasta file, and returns fully qualified name
    Checks if a fasta file exists returns the file attaced to a region
    or something 
    
    """
    
    if not os.path.exists(directory):
        raise IOError("directory: %s, doesn't exist" % (directory))
    
    if region is not None:
        full_name = filename + "." +  region + "." + type + ".fa"
        return os.path.join(directory, full_name)
    else:
        full_name = filename + "." + type + ".fa"
        return os.path.join(directory, full_name)


def make_dir(dir_name):
    
    """ makes a dir, dir_name if it doesn't exist otherwise does nothing """
    
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
            

def count_total_reads(genomic_regions, counts):
    gas = HTSeq.GenomicArrayOfSets( "auto", stranded=True)

    for region_name, cur_region in genomic_regions.items():
        for interval in bed_to_genomic_interval(cur_region):
            gas[interval] += region_name

    total_counts = defaultdict(int)
    for overlapping_interval, count in counts.steps():
        for step in gas[overlapping_interval].steps():
            for item in step[1]:
                total_counts[item] += count
    return total_counts

def count_reads_per_cluster(bedtool, counts):
    
    """
    
    Counts the number of reads in each cluster
    
    bedtool: bedtool containing the specially formatted CLIPper bedlines
    bamtool: a pybedtools tool for all the reads in the bam file
    returns list(int) each index being a read in the cluster
    
    """
    #generates data for figure 1 and 2
    #gets reads in clusters (figure 1)
    #gets reads per cluster (figure 2)

    reads_per_cluster = []
    try:
        for cluster in bedtool:
            #handles merged clusters (taking the reads from only the first one)
            first_cluster = cluster.name.split(";")[0]
            number_reads_in_peak = first_cluster.split("_")[2]
            reads_per_cluster.append(int(number_reads_in_peak))

    #can't count based on reported number, reverting to HTSeq
    except:
        result = get_densities(bedtool, counts)
        reads_per_cluster = list([sum(item) for item in result])

    reads_in_clusters = sum(reads_per_cluster)
    return reads_in_clusters, reads_per_cluster


def calculate_kmer_diff(kmer_list, regions, clusters, fasta_dir):
    
    """
    
    calculates results of kmer diff and returns the results in a dict
    
    kmer_list - list[int] kmer sizes compute
    regions   - list[str] regions to compute kmers on
    clusters  - str cluster name
    fasta_dir - str directoy of fasta files
    
    """
    
    if kmer_list is None:
        return None
    
    kmer_results = {}    

    for region in regions:

        real_fa = fa_file(clusters, region=region, directory=fasta_dir, type="real")
        rand_fa = fa_file(clusters, region=region, directory=fasta_dir, type="random")

        kmer_results[region] = {}
        
        for k in kmer_list:
            try:
                kmers, n1, n2 = kmer_diff(real_fa, rand_fa, k)
                kmer_results[region][k] = kmers
            except IOError as e:
                print e
                print "Ignoring: %s" % region
 
    return kmer_results


def calculate_homer_motifs(kmer_list, regions, clusters, fasta_dir, homerout):
    
    """
    
    Calculates motiffs for homer
    
    kmer_list - list[int] different kmers to examine
    regions   - list[str] regions to compute kmers on
    clusters - str name out output file
    fasta_dir - str directoy of fasta files
    homerout - str output dir of homer
    
    """
    if kmer_list is None:
        return None
    
    for region in regions:
        #reads nicely named files
        real_fa = fa_file(clusters, region=region, directory=fasta_dir, type="real")
        rand_fa = fa_file(clusters, region=region, directory=fasta_dir, type="random")

        region_homer_out = os.path.join(homerout, region)
        run_homer(real_fa, rand_fa, kmer_list, outloc=region_homer_out)


def calculate_phastcons(regions, cluster_regions, phastcons_location, sample_size=1000):
    
    """
    
    Generates phastcons scores for all clusters
    
    regions - list[str] regions to analize
    cluster_regions - special dict that contains all random and normal cluster info
    phastcons_location - str location of phastcons file
    sample_size - int the number of regions to randomly sample to calculate phastcons (needed because sampling from bw file is slow for phastcons)
    
    """
    
    
    ###conservation --should use multiprocessing to speed this part up!
    #start output_file conservation logic, very slow...
    phast_values = {"real": {}, "rand": {}}
    print "Fetching Phastcons Scores..."
    #phastcons values for all regions except "all"

    #skip "all" combine them later
    for region in regions:
        print "%s..." % region
        #gracefully fail if the region isn't represented
        try:
            phast_values["real"][region] = get_mean_phastcons(cluster_regions[region]['real'], 
                                                              phastcons_location, 
                                                              sample_size=sample_size)
            
            #can't concatanate zero length arrays, so prime it
            randPhast = np.array([])
            for rand_bed in cluster_regions[region]['rand'].values():                
                randPhast = np.concatenate((randPhast, get_mean_phastcons(rand_bed, 
                                                                          phastcons_location,  
                                                                          sample_size=sample_size)), axis=1)
            phast_values["rand"][region] = randPhast
        except KeyError as e:
            print "ignoring: ", e
            continue
        except Exception as e:
            print region, e
            raise


    #get the "all" values
    phast_values["real"]["all"] = np.concatenate(phast_values["real"].values())
    phast_values["rand"]["all"] = np.concatenate(phast_values["rand"].values())

    return phast_values


def calculate_motif_distance(cluster_regions, region_sizes, motif_tool, slopsize = 200):
    
    """
    
    Calculates distance from each cluster to the nearest motif specified in motif 
    tool.  Returns dict used for printing outputs.  Dict format is 
    
    dict{region : {'real': {'size' : int, 'dist' : list[int (distance to motif)},
                   'rand': {'size' : int, 'dist' : list[int (distance to motif)}}
    
    Input: 
    cluster_regions - special dict of clusters
    region_sizes - dict[region] : int count of number of the size of the region
    motif_tool - bedtool of that contains the location of that motif in the genome
        
    """
    
    dist_dict = {}
    #for each region that exists calculate stats
    for region in cluster_regions.keys():
        dist_dict[region] = {}
        dist_dict[region]['real'] = {}
        dist_dict[region]['rand'] = {}
        
        dist_dict[region]['real']['dist'] = get_motif_distance(cluster_regions[region]['real'], motif_tool, slop=slopsize)
        dist_dict[region]['real']['size'] = region_sizes[region]

        nrand = len(cluster_regions[region]['rand'].keys())
        dist_dict[region]['rand']['dist'] = []
        for i in range(nrand):
            dist_dict[region]['rand']['dist'] += get_motif_distance(cluster_regions[region]['rand'][i], motif_tool, slop=slopsize)
        dist_dict[region]['rand']['size'] = region_sizes[region] * nrand
        
        print "%s done" % region
    
    #assigns 'all' information
    dist_dict['all'] = {}
    dist_dict['all']['real'] = {'dist': [], 'size': 0}
    dist_dict['all']['rand'] = {'dist': [], 'size': 0}
    for region in dist_dict:
        dist_dict['all']['real']['size'] += dist_dict[region]['real']['size']
        dist_dict['all']['rand']['size'] += dist_dict[region]['rand']['size']
        dist_dict['all']['real']['dist'] += dist_dict[region]['real']['dist']
        dist_dict['all']['rand']['dist'] += dist_dict[region]['rand']['dist']

    return dist_dict


def get_motif_distance(clusters, motif, slop=500):
    
    """
    
    Compares two bed files and computes distance from center of first (indicated by bed12)
    to center of second (by bed12)
    
    Input:
      
    clusters - bedtool (bed12)
    motif - bedtool (bed12)

    returns distance from clusters to nearest motif 
    
    """
    
    windowed_clusters = clusters.window(motif, w=slop, sm=True)
    distances = []
    for line in windowed_clusters:
        #gets the midpoint of the cluster and the motif and finds the distance between them
        cluster_center = (int(line[7]) + int(line[6])) / 2
        motif_center = (int(line[15]) + int(line[14])) / 2
        distance = motif_center - cluster_center

        if line[5] == "-":
            distance *= -1
        distances.append(distance)

    return distances


def generate_motif_distances(cluster_regions, region_sizes, motifs, motif_location, species):
    
    """
    
    Generates all motif distances for a lsit of motifs 
    returns list[motif_distances]
    
    motif_location - str location that motifs are stored
    species - str species (for finding stored motifs)
    motifs - list of motifs to analize
    cluster_regions - dict from parse clusters 
    
    """
    
    motif_distance_list = []
#given a specific motif in a motif file generate distances from that motif...?
    for motif in motifs:
        mf = "motif_" + motif + ".BED"
        mfgz = "motif_" + motif + ".BED.gz"

        motif_tool = None
        if os.path.exists(os.path.join(motif_location, species, mf)):
            motif_tool = pybedtools.BedTool(os.path.join(motifBASE, species, mf))
        elif os.path.exists(os.path.join(motif_location, species, mfgz)):
            motif_tool = pybedtools.BedTool(os.path.join(motif_location, species, mfgz))
        else:
            print "MOTIF BED FILE for motif: %s is not available, please build it" % (mf)
        
        if motif_tool is not None:
            motif_distance_list.append(calculate_motif_distance(cluster_regions, region_sizes, motif_tool))
    
    return motif_distance_list


def make_unique(bedtool):
        for x, interval in enumerate(bedtool):
            interval.name = interval.name + "_" + str(x)
            yield interval


def visualize(clusters, extension, out_dir):
    with open(os.path.join("%s.clip_analysis.pickle" % clusters)) as out_file:
        clip_viz = CLIP_analysis_display.ClipVisualization(out_file)
    qc_fig = clip_viz.CLIP_QC_figure()
    distribution_fig = clip_viz.plot_distributions()
    qc_fig.savefig(os.path.join(out_dir, clusters + ".qc_fig." + extension))
    distribution_fig.savefig(os.path.join(out_dir, clusters + ".DistFig." + extension))


def adjust_name(interval):
    interval.name = interval[-3]
    return interval


def add_thick_marks(interval):
    interval[6] = str(interval.start)
    interval[7] = str(interval.stop)
    return interval


def smooth_interval(interval):
    """Takes only the first item merged, everything else is tossed"""
    for item in range(3,8):
        interval[item] = interval[item].split(",")[0]
    return interval


def infer_info(bedtool, genes):
    bedtool = bedtool.intersect(genes, s=True, loj=True).\
                      sort(). \
                      each(adjust_name). \
                      each(add_thick_marks). \
                      merge(s=True, c="4,5,6,7,8", o="distinct,distinct,distinct,distinct,distinct"). \
                      each(smooth_interval). \
                      saveas()
    return bedtool


def main(bedtool, bam, species, runPhast=False, motifs=[], k=[6], nrand=3,
         outdir=os.getcwd(), db=None, as_structure=None, genome_location=None, phastcons_location=None,
         regions_location=None, motif_location=os.getcwd(), metrics="CLIP_Analysis.metrics", extension="svg",
         infer=False):
    
    """
    
    Runs all analysies 
    
    one thing to do is make graphs fail gracefully 
    
    """
    print "starting"
    #gets clusters in a bed tools + names species 
    clusters = os.path.basename(bedtool)
    species = species
    out_dict = {}
    #In case names aren't unique make them all unique
    clusters_bed = pybedtools.BedTool(make_unique(pybedtools.BedTool(bedtool))).saveas()
    if len(clusters_bed) <= 1:
        raise IllegalArgunmentException("not enough reads to properly analyze bed file")
    coverage = get_bam_coverage(bam)
    counts = get_bam_counts(bam)

    #makes output file names 
    clusters = str.replace(clusters, ".BED", "")
    k = [int(x) for x in k]
    
    #sets up output dirs
    make_dir(outdir)        

    assigned_dir = os.path.join(outdir, "assigned")
    misc_dir = os.path.join(outdir, "misc")
    fasta_dir = os.path.join(outdir, "fasta")    
    homerout_base = os.path.join(outdir, "homer")
    make_dir(homerout_base)
    homerout = os.path.join(homerout_base, clusters)    

    make_dir(assigned_dir)
    make_dir(misc_dir)
    make_dir(fasta_dir)
    make_dir(homerout)
    
    #lazy refactor, this used to be a list, but the dict acts the same until I don't want it to...
    regions = OrderedDict()
    regions["all"] = "All"
    regions["cds"] = "CDS"
    regions["three_prime_utrs"] = "3' UTR"
    regions["five_prime_utrs"] = "5' UTR"
    regions["proxintron500"] = "Proximal\nIntron"
    regions["distintron500"] = "Distal\nIntron"
    if db is not None:
        db = gffutils.FeatureDB(db)
    else:
        print "gff utils db not defined, this is fine, but falling back onto pre-set region defentions"
        db = None
    features = GenomicFeatures(species, db,  regions_dir=regions_location)
    genomic_regions = features.get_genomic_regions()
    features = features.get_feature_locations()

    clusters_bed = pybedtools.BedTool(bedtool)
    if infer:
        clusters_bed = infer_info(clusters_bed, genomic_regions['genes'])

    clusters_bed = pybedtools.BedTool(make_unique(clusters_bed)).saveas()

     #all catagory would break some analysies, create copy and remove it
    assigned_regions = regions.copy()
    del assigned_regions['all']

    if as_structure is not None:
        genes_dict, genes_bed = parse_AS_STRUCTURE_dict(species, as_structure)
    else:
        print "AS STRUCTURE file not listed, alt-splicing figure will not be generated"

    cluster_regions = assign_to_regions(tool=clusters_bed, clusters=clusters, regions=assigned_regions,
                                        assigned_dir=assigned_dir, species=species, nrand=nrand)

    print "getting cluster sizes"
    region_sizes = get_sizes(cluster_regions)

    print "getting genomic regions sizes"
    genic_region_sizes = count_genomic_region_sizes(assigned_regions, species)
    print "counting reads per cluster"
    reads_in_clusters, reads_per_cluster = count_reads_per_cluster(cluster_regions['all']['real'], counts)

    print "counting total reads in regions"
    region_read_counts = count_total_reads(genomic_regions, counts)
    total_reads = region_read_counts['genes']

    print "clustering peaks"
    read_densities, classes = cluster_peaks(cluster_regions['all']['real'], coverage)

    #generates cluster lengths (figure 3)
    print "getting cluster lengths"
    cluster_lengths = bedlengths(cluster_regions['all']['real'])

    print "getting peak locations"
    features_transcript_closest, features_mrna_closest = get_feature_distances(cluster_regions['all']['real'], genomic_regions, features)
    features_transcript_closest = {name:  pybedtools.BedTool(bedtool['dist'].saveas("%s_%s_transcript.bed" % (clusters, name)).fn) for name, bedtool in features_transcript_closest.items() if bedtool['dist'] is not None}
    features_mrna_closest = {name:  pybedtools.BedTool(bedtool['dist'].saveas("%s_%s_mrna.bed" % (clusters, name)).fn) for name, bedtool in features_mrna_closest.items() if bedtool['dist'] is not None}

    #Distribution counting only works if genes have been pre-assignned
    distributions = get_region_distributions(cluster_regions['all']['real'], genomic_regions)

    if as_structure is not None:
        #also builds figure 10 (exon distances)
        genomic_types = count_genomic_types(genes_dict)
        exon_types = get_closest_exon_types(cluster_regions['all']['real'], genes_dict)

        #generates figure 10 (exon distances)
        type_count = [exon_types["CE:"], exon_types["SE:"], exon_types["MXE:"], exon_types["A5E:"], exon_types["A3E:"]]

        genomic_type_count = [genomic_types["CE:"], genomic_types["SE:"], genomic_types["MXE:"], genomic_types["A5E:"],
                              genomic_types["A3E:"]]

        out_dict["genomic_type_count"] = genomic_type_count
        out_dict["type_count"] = type_count

    if genome_location is not None:
        make_fasta_files_from_regions(cluster_regions, clusters, fasta_dir, genome_location)
        calculate_homer_motifs(k, regions, clusters, fasta_dir, homerout)
        kmer_results = calculate_kmer_diff(k, regions, clusters, fasta_dir)
    else:
        print "No genome fasta file provide, motif identification will not be performed"

    phast_values = None
    print "starting phast"
    if runPhast:
        phast_values = calculate_phastcons(assigned_regions, cluster_regions, phastcons_location)
    print "ending phast"

    motif_distances = []
    try:
        if motifs:
            motif_distances = generate_motif_distances(cluster_regions, region_sizes, motifs, motif_location,
                                                       species)
    except:
        pass

    out_dict['features_transcript_closest'] = features_transcript_closest
    out_dict['features_mrna_closest'] = features_mrna_closest
    out_dict['distributions'] = distributions
    out_dict["region_sizes"] = region_sizes
    out_dict["reads_in_clusters"] = reads_in_clusters
    out_dict["reads_out_clusters"] = (total_reads - reads_in_clusters)
    out_dict["cluster_lengths"] = cluster_lengths
    out_dict["reads_per_cluster"] = reads_per_cluster
    out_dict["genic_region_sizes"] = genic_region_sizes
    out_dict["kmer_results"] = kmer_results
    out_dict["motifs"] = motifs
    out_dict["phast_values"] = phast_values
    out_dict["motif_distances"] = motif_distances
    out_dict['data'] = np.array(read_densities)
    out_dict['classes'] = classes
    out_dict['region_read_counts'] = region_read_counts
    out_dict['homerout'] = homerout
    out_dict['regions'] = regions

    with open(os.path.join("%s.clip_analysis.pickle" % clusters), 'w') as out_file:
        pickle.dump(out_dict, file=out_file)
    print "file saved"

    visualize(clusters, extension, outdir)
      
    #prints distance of clusters from various motifs in a different figure
    #TODO use homer / MEME results to print this out by default?  Make different sub package?
    try:
        if motifs:
            motif_fig = CLIP_analysis_display.plot_motifs(motif_distances)
            motif_fig.savefig(clusters + ".motif_distribution." + extension)
    except:
        pass
    
    with open(metrics, 'w') as outfile:
        outfile.write("FRiP\n")
        outfile.write("\t".join([str(float(reads_in_clusters) / float(total_reads))]) + "\n")
    

def call_main():  
    parser = OptionParser()
    
    parser.add_option("--clusters", dest="clusters", help="BED file of clusters", metavar="BED")
    parser.add_option("--bam", dest="bam", help="The bam file from the CLIP analysis")
    parser.add_option("--species", "-s", dest="species", help = "genome version")
    parser.add_option("--runPhast", dest="runPhast", action="store_true", default=False, help="Run Phastcons ")
    parser.add_option("--motifs", dest="motifs", action="append", help="Motifs to use (files of motifs give must exist in motif_directory directory)", default=[])
    parser.add_option("--k", dest="k", action="append", help="k-mer and homer motif ananlysis", default=[6])
    parser.add_option("--nrand", dest="nrand", default=3, help="selects number of times to randomly sample genome", type="int")
    parser.add_option("--outdir", "-o", dest="outdir", default=os.getcwd(), help="directory for output, default:cwd")
    ##Below here are critical files that always need to be referenced
    parser.add_option("--gff_db", dest="db", help="gff database from gffutils to generate annotations with")
    parser.add_option("--AS_Structure", dest="as_structure",  help="Location of AS_Structure directory (chromosme files should be inside)", default=None)
    parser.add_option("--genome_location", dest="genome_location", help="location of all.fa file for genome of interest", default=None)
    parser.add_option("--phastcons_location", dest="phastcons_location",  help="location of phastcons file", default=None)
    parser.add_option("--regions_location", dest="regions_location",  help="directory of genomic regions for a species (default: clipper defined regions)", default=None)
    parser.add_option("--motif_directory", dest="motif_location",  help="directory of pre-computed motifs for analysis", default=os.getcwd())
    parser.add_option("--metrics", dest="metrics", default="CLIP_Analysis.metrics", help="file name to output metrics to")
    parser.add_option("--extension", dest="extension", default="svg", help="file extension to use (svg, png, pdf...)")
    parser.add_option("--infer", default=False, action="store_true", help="Infer peak center and gene if if peak finding algorithm doesn't report it")

    (options, args) = parser.parse_args()
    
    #error checking
    if options.clusters is None or options.bam is None or options.species is None:
        parser.print_help()
        exit()
    main(bedtool=options.clusters, bam=options.bam, species=options.species,
         runPhast=options.runPhast,
         motifs=options.motifs, k=options.k, nrand=options.nrand, outdir=options.outdir, db=options.db,
         as_structure=options.as_structure,
         genome_location=options.genome_location,
         phastcons_location=options.phastcons_location, regions_location=options.regions_location,
         motif_location=options.motif_location, metrics=options.metrics, extension=options.extension,
         infer=options.infer,
         )

if __name__ == "__main__":
    call_main()
