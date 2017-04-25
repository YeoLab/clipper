

"""

Analizes CLIP data given a bed file and a bam file

Michael Lovci and Gabriel Pratt

"""
from collections import Counter, OrderedDict, defaultdict
import os
import subprocess

from bx.bbi.bigwig_file import BigWigFile
import HTSeq
import numpy as np
import pandas as pd
import pybedtools
from sklearn.cluster import KMeans

import clipper
from clipper.src.kmerdiff import kmer_diff
from clipper.src.bam_helpers import Robust_BAM_Reader


def convert_to_mRNA_position(interval, gene_model):
    """

    Returns distance from nearest feature assuming distance is centered on the feature

    bedtool - a bedtool to find closest feature of
    gene_model - dict of lists of pybedtools

    Returns bed objects mapped to mRNA position instead of genomic position

    Assumes both the bedtools object and the feature are 1bp long so we get the distance from both from their start sites

    Negative strand gets modified to be positive strand like, this will fuck with closest bed
    need to do everything on the positive strand from here on out
    """

    #feature_dict = {feature.name : feature for feature in closest_feature}


    if interval.chrom not in gene_model:
        interval.chrom = "none"
        return interval
        #raise KeyError(interval.chrom + " not in current as stucture dict ignoring cluster ")

    #gene model - dict of list of intervals, always at least length 1
    if not interval.strand == gene_model[interval.chrom][0].strand:
        interval.chrom = "none"
        return interval
        #raise ValueError("strands not the same, there is some issue with gene annotations")

    running_length = 0

    for region in gene_model[interval.chrom]:

        if interval.start >= int(region.start) and interval.start <= int(region.stop):
            if interval.strand == "+":
                tmp_start = running_length + (interval.start - region.start)
                tmp_end = running_length + (interval.end - region.start)

            elif interval.strand == "-": #need the temps for swaping start and end
                tmp_start = running_length + (region.stop - interval.end)
                tmp_end = running_length + (region.stop - interval.start)

            if int(tmp_start) <= 0 or int(tmp_end) <= 0:
                tmp_start = 1
                tmp_end = 1

            interval.start = tmp_start
            interval.end = tmp_end

            return interval
        running_length += region.length
    interval.chrom = "none"
    return interval


def small_peaks(feature):
    """

    feature - pybedtools feature

    returns center of clipper called peak (the middle of the narrow start / stop)

    """

    try:
        feature.start = (int(feature[6]) + int(feature[7])) / 2
        feature.stop = ((int(feature[6]) + int(feature[7])) / 2) + 1
    except IndexError:
        #if I don't have that # fall back on just taking the peak center
        midpoint = (feature.start + feature.stop) / 2
        feature.start = midpoint
        feature.stop = midpoint + 1

    feature.name = feature.name.split("_")[0]
    return feature


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



def fix_strand(interval, warn=True):
    #this only is comptabale with bedtools 2.25.0
    #lst = interval.fields
    #del lst[3]
    #interval = pybedtools.interval_constructor(lst)
    strands = list(set(interval.strand.split(",")))
    if len(strands) > 1 and warn:
        #There is something odd going on between my local box, pybedtools and bedtools.  Merge isn't acting exacly the same, this works on TSCC
        raise NameError("Both strands are present, something went wrong during the merge")
    interval.strand = strands[0]
    return interval

def fix_strand_v26(interval, warn=True):
    #this only is comptabale with bedtools 2.26.0
    interval = list(interval)
    interval.pop(3)
    interval = pybedtools.create_interval_from_list(interval)

    strands = list(set(interval.strand.split(",")))
    
    if len(strands) > 1 and warn:
        #There is something odd going on between my local box, pybedtools and bedtools.  Merge isn't acting exacly the same, this works on TSCC
        print interval
        print interval.strand
        print strands
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


def move_name(interval, original_length):
    interval.name = interval[original_length + 3]
    return interval


def fix_name(interval):
    """
    currently pybedtools crashes if the 1st and 3rd columns is an int and there are more than 11 lines
    This fixes that issue but turing the name column to foo
    """

    interval.name = "foo"
    return interval

def assign_to_regions(tool, clusters=None, assigned_dir=".", species="hg19", nrand=3):
    
    """
    
    Assigns each cluster to a genic region
    finally saves all generated bed and fasta files for future analysis...

    tool - a bed tool (each line represnting a cluster)
    clusters - name of cluster file (optional)
    assigned_dir - location to save files in
    species - str species to segment
    nrand - int number offsets times to shuffle for null hypothesis


    """
    if clusters is None:
        clusters, ext = os.path.splitext(os.path.basename(tool.fn))
    bedtracks = {}

    regions, assigned_regions = regions_generator()
    short_species = species.split("_")[0]
    if short_species == "GRCh38":
        short_species = "hg38"

    for region in regions:
        bedtracks[region] = pybedtools.BedTool(os.path.join(clipper.data_dir(), "regions", "%s_%s.bed" % (species,
                                                                                                          region)))
    #creates the basics of bed dict
    bed_dict = {'all': {'rand': {}}}

    genes = pybedtools.BedTool(os.path.join(clipper.data_dir(), "regions", "%s_genes.bed" % (species)))


    offsets = get_offsets_bed12(tool)
    if tool.field_count() <= 5:
        tool.sort().merge().saveas()
    elif 6 <= tool.field_count() < 8:
        #Hack to get around not having gene name assigned by peak caller, due to overlapping genes this won't be perfect
        #move_name_real = functools.partial(move_name, original_length=len(tool[0].fields))
        #tool = tool.intersect(genes, wo=True, s=True).each(move_name_real).saveas()
        #fix_strand_ok = functools.partial(fix_strand, warn=False)
        tool = tool.sort().merge(s=True, c="4,5,6", o="collapse,collapse,collapse").each(fix_strand_v26).saveas()
    #elif not tool[0][7].isdigit():
    #    tool = tool.sort().merge(s=True, c="4,5,6", o="collapse,collapse,collapse").each(fix_strand).each(fix_name).saveas()
    else: #Clipper, this is ideal we like this technique
        tool = tool.sort().merge(s=True, c="4,5,6,7,8", o="collapse,collapse,collapse,min,min").each(fix_strand_v26).saveas()

    remaining_clusters = adjust_offsets(tool, offsets)

    # print "There are a total %d clusters I'll examine" % (len(tool))
    for region in regions:
        remaining_clusters, overlapping = intersection(remaining_clusters, b=bedtracks[region])

        #if for some reason there isn't a peak in the region skip it
        if len(overlapping) == 0:
            # print "ignoring %s " % region
            continue

        #sets up bed dict for this region
        bed_dict[region] = {'real': overlapping.sort(stream=True).saveas(),
                            'rand': {}}

        no_overlapping_count = len(remaining_clusters)
        overlapping_count = len(bed_dict[region]['real'])
        # print "For region: %s found %d that overlap and %d that don't" % (region,
        #                                                                   overlapping_count,
        #                                                                   no_overlapping_count)

        if 'real' not in bed_dict['all']:
            bed_dict['all']['real'] = bed_dict[region]['real']
        else:
            bed_dict['all']['real'] = bed_dict['all']['real'].cat(bed_dict[region]['real'], stream=True, postmerge=False).saveas()

        #saves offsets so after shuffling the offsets can be readjusted
        offset_dict = get_offsets_bed12(bed_dict[region]['real'])
        for i in range(nrand):
            random_intervals = bed_dict[region]['real'].shuffle(genome=short_species, incl=bedtracks[region].fn).sort()
            random_intervals = fix_shuffled_strand(random_intervals, bedtracks[region].fn)
            random_intervals = adjust_offsets(random_intervals, offset_dict)
            bed_dict[region]['rand'][i] = random_intervals.saveas()

            if i not in bed_dict['all']['rand']:
                bed_dict['all']['rand'][i] = bed_dict[region]['rand'][i]
            else:
                bed_dict['all']['rand'][i] = bed_dict['all']['rand'][i].cat(bed_dict[region]['rand'][i], stream=True, postmerge=False)


        #if there are no more clusters to assign stop trying
        if no_overlapping_count == 0:
            break

    # print "After assigning %d un-categorized regions" % len(remaining_clusters)

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

    if bedtool.field_count() < 8:
        print "Not Valid bed12 file, continuing processing, some things may be strange"
        return None

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
        print "Not Valid bed12 file, continuing processing, some things may be strange, also this will cause a file leak in pybedtools, watch out"
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


def convert_to_transcript(bedtool):
    
    """
    
    Converts bed file to be based only off transcripts (transcripts are defined by chromosome)
    
    returns modified dict
        
    """

    return bedtool.each(name_to_chrom).sort().saveas()


def convert_to_mrna(bedtool, exon_dict):
    
    """
    
    converts transcript dict into mRNA locations given a dict of exons
    
    feature_dict - generated from convert to transcript, dict of bedfiles 
    exon_dict - dict of { genes : [exons] }
    
    """

    return bedtool.each(convert_to_mRNA_position, exon_dict).filter(lambda x: x.chrom != "none").sort().saveas()


def invert_neg(interval):
    interval[-1] = str(int(interval[-1]) * -1)
    return interval


def get_feature_distances(bedtool, features, exons):
    """

    :param bedtool: bedtools of peaks, both  (only works with clipper defined peaks)
    :param features: dict of features to find distance from
    :param exons: bedtool of exons, used to convert genomic coordinates to mrna coordinates

    :return:
    """

    exon_dict = generate_region_dict(exons)
    
    bed_center = bedtool.each(small_peaks).saveas()
    beds_center_transcripts = convert_to_transcript(bed_center)
    beds_center_transcripts_mrna = convert_to_mrna(beds_center_transcripts, exon_dict)

    features_transcript = {name: convert_to_transcript(bedtool) for name, bedtool in features.items()}
    features_mrna = {name: convert_to_mrna(bedtool, exon_dict) for name, bedtool in features_transcript.items()}

    #for pre-mrna
    features_transcript_closest = defaultdict(lambda: None)
    for name, feature in features_transcript.items():
        features_transcript_closest[name] = beds_center_transcripts.closest(feature, s=True, D="b").filter(lambda x: x[-2] != ".").saveas()

    #for mrna
    features_mrna_closest = defaultdict(lambda: None)
    for name, feature in features_mrna.items():
        features_mrna_closest[name] = beds_center_transcripts_mrna.closest(feature, s=True, D="ref").filter(lambda x: x[-2] != ".").each(invert_neg).saveas()

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
    return None, None, None, None


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

    feature_tool = pybedtools.BedTool(bedtool_list).sort()
    return Counter([interval[-1] for interval in bedtool.closest(feature_tool, s=True)])

#TODO Start small module for getting read densiites


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
        yield HTSeq.GenomicInterval(str(interval.chrom), interval.start, interval.stop + 1, str(interval.strand))


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
    #I'm calling peaks in the wrong places?  This is wrong, I should look into it.  I shouldn't need to dropnas
    bedtool_df_mag_normalized = bedtool_df_mag_normalized.dropna()

    #TODO Write Test, tested slowly in ipython notebook, but need to make some sample data to test with
    data = np.array(bedtool_df_mag_normalized)
    #fixes bad test case
    if len(data) < k:
        #In case there are very few peaks in a region the max will keep k at least 1
        k = max(len(data) / 2, 1)
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

def trim_bedtool(bedtool):
    result = []
    for x, interval in enumerate(bedtool):
        interval.name = "{}_{}:{}-{}".format(interval.name, interval.chrom, interval.start, interval.stop)
        result.append(interval[:4])
    return pybedtools.BedTool(result).saveas()

def bigWigAverageOverBed(bw_file, bedtool):

    out_name = bedtool.fn.split(".")[:-1]
    if len(out_name) is not 1:
        out_name = ".".join(out_name)

    outfile = os.path.join(os.getcwd(), out_name + ".tab")
    with open(os.devnull, 'w') as fnull:
        print 'bigWigAverageOverBed', bw_file, trim_bedtool(bedtool).fn, outfile
        subprocess.check_call(" ".join(['bigWigAverageOverBed',
                         bw_file,
                         trim_bedtool(bedtool).fn,
                         outfile]), shell=True)
    return pd.read_table(outfile, index_col=0, header=None, names=['name', 'size', 'covered', 'sum', 'mean0', 'mean'])

def calculate_phastcons(cluster_regions, phastcons_location):
    result = {}
    for region in cluster_regions.keys():
        result[region] = {}
        region_dict = {}
        region_dict['real'] = pd.concat({1: bigWigAverageOverBed(phastcons_location, cluster_regions[region]['real'])})
        if 'rand' in cluster_regions[region]:
            rand = {}
            for rand_num in cluster_regions[region]['rand']:
                rand[rand_num] = bigWigAverageOverBed(phastcons_location, cluster_regions[region]['rand'][rand_num])
            region_dict['rand'] = pd.concat(rand)
        result[region] = pd.concat(region_dict)
    result = pd.concat(result)
    return result


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


def regions_generator():

    """
    returns ordered dict of regions without the "all" region and with the "All" region, in that order
    """

    #lazy refactor, should turn entire thing into object, make this a field
    regions = OrderedDict()
    regions["all"] = "All"
    regions["cds"] = "CDS"
    regions["three_prime_utrs"] = "3' UTR"
    regions["five_prime_utrs"] = "5' UTR"
    regions["proxintron500"] = "Proximal\nIntron"
    regions["distintron500"] = "Distal\nIntron"

    #all catagory would break some analysies, create copy and remove it
    assigned_regions = regions.copy()
    del assigned_regions['all']
    return assigned_regions, regions
