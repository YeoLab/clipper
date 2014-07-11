

"""

Analizes CLIP data given a bed file and a bam file

Michael Lovci and Gabriel Pratt

"""
from collections import Counter, OrderedDict, defaultdict
from itertools import izip_longest
from optparse import OptionParser
import os
import pickle
from random import sample
import subprocess

from bx.bbi.bigwig_file import BigWigFile
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
from gscripts.general.pybedtools_helpers import small_peaks, shuffle_and_adjust, closest_by_feature, get_three_prime_end, get_five_prime_end, convert_to_mRNA_position, adjust_after_shuffle

import matplotlib as mpl
from matplotlib import rc

mpl.rcParams['svg.fonttype'] = 'none'
rc('text', usetex=False)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def name_to_chrom(interval):
    interval.chrom = interval.name
    return interval

def intersection(A, B=None):
    
    """
    
    A : bedtool
    B : bedtool
    Returns A - B and A intersect B 
    A with b and returns everything in a but not b and everything in a but... ???
    
    """
    
    overlapping = A.intersect(B, wa=True, u=True).saveas()
    remaining   = A.intersect(B, wa=True, v=True).saveas()
    #overlapping = A.subtract(B, s=True, A=True).saveas() #without regions that overlap
    #print len(less)
    #remaining = A.subtract(less, s=True, A=True) #only regions that overlap
    return remaining, overlapping 


def adjust_offsets(tool, offsets=None):
    
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
    
    if offsets == None:
        raise Exception, "no offsets"
    
    clusters = []
    for bed_line in tool:

        #handles multiple different types of bed files
        if len(bed_line.fields) < 8:
            thick_start = 0
            thick_stop = 0
        else:
            thick_start = int(bed_line[6])
            thick_stop = int(bed_line[7])
        
        #the ; represents two merged locations 
        if "," in bed_line.name and bed_line.name not in offsets:
            offset = offsets[bed_line.name.split(",")[0]]
        else:
            offset = offsets[bed_line.name]
        
        if bed_line.strand == "+":
            thick_start = bed_line.start + offset
            thick_stop = thick_start + 4
        else:
            thick_stop = bed_line.stop - offset
            thick_start = thick_stop - 4
            
        clusters.append("\t".join([str(x) for x in (bed_line.chrom, 
                                                    bed_line.start, 
                                                    bed_line.stop, 
                                                    bed_line.name, 
                                                    bed_line.score, 
                                                    bed_line.strand, 
                                                    thick_start, 
                                                    thick_stop)]))
    
    return pybedtools.BedTool("\n".join(clusters), from_string=True)
    



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



def count_genomic_region_sizes(regions_dir, regions, species="hg19"):
    
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

def build_genomic_regions(tool, prox_distance=500):
    
    """
    
    Makes dict of seperated genomic regions from gtf file
    
    tool: gtf formatted bedfile with gene definitions.  Must contain defenitions 3UTR, 5UTR, CDS and exon
    prox_distance, the distance to define an intron as proximal
    returns dict of regions 5UTR, 3UTR, CDS, proxintron, distintron
    
    """
    
    prev_line = None
    proximal_introns = ""
    distal_introns   = ""
    
    for line in tool.filter(lambda read: read[2] == "exon"):
        if prev_line is None or not prev_line.name == line.name:
            prev_line = line
            continue
        
        #adjust by one to prevent overlapping regions
        
        if line.start - prev_line.stop < prox_distance * 2:
             proximal_introns += "\t".join([str(x) for x in [prev_line.chrom, "stdin", "proxintron", prev_line.stop + 1, line.start - 1, prev_line.name, prev_line.strand, "gene_id \"%s\";" % (line.name), "\n"]])
        else:
            #May be an off by one bug here, but it doesn't matter to much... (hopefully)
            proximal_introns += "\t".join([str(x) for x in [prev_line.chrom, "stdin", "proxintron", prev_line.stop + 1, prev_line.stop + prox_distance, prev_line.name, prev_line.strand, "gene_id \"%s\";" % (line.name),  "\n"]])
            proximal_introns += "\t".join([str(x) for x in [line.chrom, "stdin", "proxintron", line.start - prox_distance, line.start - 1, line.name, line.strand, "gene_id \"%s\";" % (line.name), "\n"]])
            
            #adjust by 2 bp so we don't overlap the stop / start site of either the exon or the proxintron
            distal_introns += "\t".join([str(x) for x in [line.chrom, "stdin", "proxintron",  prev_line.stop + prox_distance + 2, line.start - prox_distance - 2, line.name, line.strand, "gene_id \"%s\";" % (line.name), "\n"]])

        prev_line = line
        
    proximal_introns = pybedtools.BedTool(proximal_introns, from_string=True)
    distal_introns   = pybedtools.BedTool(distal_introns, from_string=True)
    UTR_5 = tool.filter(lambda read: read[2] == "5UTR").saveas()
    UTR_3 = tool.filter(lambda read: read[2] == "3UTR").saveas()
    CDS  = tool.filter(lambda read: read[2] == "CDS").saveas()
    
    return {"CDS" : CDS, "UTR5" : UTR_5, "UTR3" : UTR_3, "proxintron" : proximal_introns, "distintron" : distal_introns}

def assign_to_regions(tool, clusters, speciesFA, regions_dir, regions, 
                      assigned_dir, fasta_dir, species="hg19", nrand = 3, 
                      getseq=False):
    
    """
    
    Assigns each cluster to a genic region
    finally saves all generated bed and fasta files for future analysis...
    
    Needs to be refactored desperatly
    
    tool - a bed tool (each line represnting a cluster)
    clusters - name of cluster file
    speciesFA - the species fasta file
    regions_dir - the directory that has genomic regions already processed
    regions - dict [str] regions to process, not used now but should be after refactoring 
    assigned_dir - location to save files in
    fasta_dir -location to save fasta files in
    species - str species to segment
    nrand - int number offsets times to shuffle for null hypothesis
    getseq - boolean gets the full sequence to store
    
    shuffling for background, should still be factored out
    
    """
    
    
    
    bedtracks = {}

    for region in regions:
        bedtracks[region] = pybedtools.BedTool(os.path.join(clipper.data_dir(), "regions", species + "_" + region + ".bed"))
              
    #creates the basics of bed dict
    bed_dict = {}
    bed_dict['all'] = {}
    bed_dict['all']['rand'] = {}

    #can't append None to string so hack a null bed tool
    for n in range(nrand):
        bed_dict['all']['rand'][n] = pybedtools.BedTool("", from_string=True)

    bed_dict['all']['real'] = pybedtools.BedTool("", from_string = True)
    
    offsets = get_offsets_bed12(tool)
    tool = tool.merge(s=True, nms=True, scores="max")
    remaining_clusters = adjust_offsets(tool, offsets)
    
    print "There are a total %d clusters I'll examine" % (len(tool))
    sizes = {}
    
    for region in regions:
        output_file = clusters + "." + region + ".real.BED"
        
        #handles case of no more things to assign... I think
 
        remaining_clusters, overlapping = intersection(remaining_clusters, 
                                                  B = bedtracks[region])  

        
        #if for some reason there isn't a peak in the region skip it
        if len(overlapping) == 0:
            print "ignoring %s " % (region)
            continue
        
        #sets up bed dict for this region
        bed_dict[region] = {}
        bed_dict[region]['rand'] = {}
        
        bed_dict[region]['real'] = overlapping.filter(is_valid_bed12).sort().saveas(os.path.join(assigned_dir, output_file))
        remaining_clusters = pybedtools.BedTool(str(remaining_clusters.filter(is_valid_bed12)), from_string=True)
        
        no_overlapping_count = len(remaining_clusters)
        overlapping_count = len(bed_dict[region]['real'])
        print "For region: %s I found %d that overlap and %d that don't" % (region, overlapping_count, no_overlapping_count)
              
        bed_dict['all']['real'] = pybedtools.BedTool(str(bed_dict['all']['real']) + str(bed_dict[region]['real']), from_string=True)    
        
        #this should be factored out
        sizes[region] = bed_dict[region]['real'].total_coverage()
        
        region_fa = fa_file(clusters, region = region, directory=fasta_dir, type="real")
        bed_dict[region]['real'].sequence(fi=speciesFA, s=True).save_seqs(region_fa)

        #saves offsets so after shuffling the offsets can be readjusted    
        offset_dict = get_offsets_bed12(overlapping)
    
        #TODO refactor to different function
        rand_list = []
        for i in range(nrand):
            output_file = clusters + "." + region + ".rand." + str(i) + ".BED"
            
            #for each region shuffles all peaks in that region around the region 
            #then pulls out sequences if requested 
            random_intervals = bed_dict[region]['real'].shuffle(genome=species, incl=bedtracks[region].fn).sort()

            #shuffling doesn't change offsets so we adjust bed 11 and 12 lines here to correct 
            random_intervals = adjust_offsets(random_intervals, offset_dict)
            
            #saves intervals, makes the all interval by appending all regions together
            bed_dict[region]['rand'][i] = random_intervals.saveas(os.path.join(assigned_dir, output_file))
            bed_dict['all']['rand'][i] = pybedtools.BedTool(str(bed_dict['all']['rand'][i]) + str(bed_dict[region]['rand'][i]), from_string=True)
            
            rand_list.append(bed_dict[region]['rand'][i].sequence(fi=speciesFA, s=True))
            
        write_seqs(fa_file(clusters, region = region, directory=fasta_dir, type="random"), rand_list)
        
        #if there are no more clusters to assign stop trying
        if no_overlapping_count == 0:
            break
    print "After assigning, I\'m left with %d un-categorized regions" %(len(remaining_clusters))

    if len(remaining_clusters) > 0:
        bed_dict['uncatagorized'] = remaining_clusters.sort()

    #Save results for the all section
 
    real_fa = fa_file(clusters, region = "all", directory=fasta_dir, type="real")

    bed_dict['all']['real'] = bed_dict['all']['real'].sort()
    bed_dict['all']['real'].sequence(fi=speciesFA,  s=True).save_seqs(real_fa)
    
    rand_list = []
    for i in range(nrand):
        output_file = clusters + ".all.rand." + str(i) + ".BED" 
        bed_dict['all']['rand'][i] = bed_dict['all']['rand'][i].sort()
        rand_list.append(bed_dict['all']['rand'][i].sequence(fi=speciesFA, s=True))
  
    write_seqs(fa_file(clusters, region = "all", directory=fasta_dir, type="random"), rand_list)
    
    return bed_dict, sizes

def build_assigned_from_existing(assigned_dir, clusters, regions, nrand):
    
    """
    
    Loads results produced from above analysis from saved file - can either be switched to
    piclking or removed 
    
    assigned_dir - location of files
    clusters - name of experiment
    regions - list of genic regions 
    nrand - number of shuffled datsets to look for
    
    """
    
    CLUS_regions = {}
    for region in regions:
        CLUS_regions[region]={}
        for n in range(nrand):
            CLUS_regions[region]['rand']={}
    for region, region_descrip in regions:
        bedfile = os.path.join(assigned_dir, "%s.%s.real.BED" %(clusters, region))
        CLUS_regions[region]['real'] = pybedtools.BedTool(bedfile)
        for n in range(nrand):
            randbedfile = os.path.join(assigned_dir, "%s.%s.rand.%s.BED" %(clusters, region, str(n)))                
            CLUS_regions[region]['rand'][n] = pybedtools.BedTool(randbedfile)
    try:
        sizes = pickle.load(open(os.path.join(assigned_dir, "%s.sizes.pickle" %(clusters)), 'rb'))
    except:
        sizes=[1,1,1,1,1]
    try:
        Gsizes = pickle.load(open(os.path.join(assigned_dir, "Gsizes.pickle"), 'rb'))
    except:
        Gsizes=[1,1,1,1,1]

    return CLUS_regions, sizes, Gsizes

def is_valid_bed12(x):
    
    """
    
    Removes clusters that start in invalid locations for bed12 formatted files
    
    x - bedline
    
    returns either true if line is valid 
    
    """
    
    chr, start, stop, name, score, strand, tstart, tstop = str(x).split("\t")
    start, stop, tstart, tstop = map(int, (start, stop, tstart, tstop))
    if start < tstart and stop > tstop:
        return True
    else:
        return False

def get_offsets_bed12(tool):
    
    """
    
    Gets offsets for each cluster in a bed12 file (CLIPper formatted bed file)
    Offsets are the difference between the wide peak and the narrow (predicted binding site)
    This will break if we analize different types of data
    tool - bedtool
    
    returns offset for each cluster
    
    """
    
    offset_dict = {}
    for line in tool:
        if line.strand == "+":
            #get difference between thick start and start
            offset = int(line[6]) - int(line[1]) 
        else:
            #get difference between thick start and start
            offset = int(line[2]) - int(line[7])
        
        #need a primary key...
        offset_dict[line.name] = offset
        
    return offset_dict


def RNA_position_interval(interval, location_dict):
    
    """
    
    makes mrna and pre-mrna peak_center figure 
    interval - single interval
    
    location_dict = dict{gene_id : {strand : "+/-", regions : list((start,stop)
    as_structure_dict - from build AS structure dict 
    
    will return distribution across entire region + just across specific region
    Might be able to use my ribo-seq stuff for genic -> transcriptomic location conversion
    
    this is based off as structure, which provides sequences ordered with first exon being the first exon on the gene, not 
    first in the chromosome (like gff does) THIS WILL NOT WORK WITH RAW GFF DATA
    
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
        raise KeyError(gene + " not in current as stucture dict ignoring cluster ")
    
    if not interval.strand == location_dict[gene][0].strand:
        raise ValueError("strands not the same, there is some issue with gene annotations")
    
    total_length = float(sum(region.length for region in location_dict[gene]))
    
    running_length = 0
    
    #reverses list if negative strand so running count can work
    #if location_dict[gene]['strand'] == "-":
    #    location_dict[gene]['regions'].reverse()
        
    for region in location_dict[gene]:
        length = float(region.length) 

        if peak_center >= int(region.start) and peak_center <= int(region.stop):
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
            return individual_fraction, total_fraction
        
        running_length += length
        
    return None, None #clusters fall outside of regions integrated

def RNA_position(interval, location_dict):
    
    """
    
    makes mrna and pre-mrna peak_center figure 
    interval - single interval
    
    location_dict = dict{gene_id : {strand : "+/-", regions : list((start,stop)
    as_structure_dict - from build AS structure dict 
    
    will return distribution across entire region + just across specific region
    Might be able to use my ribo-seq stuff for genic -> transcriptomic location conversion
    
    this is based off as structure, which provides sequences ordered with first exon being the first exon on the gene, not 
    first in the chromosome (like gff does) THIS WILL NOT WORK WITH RAW GFF DATA
    
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
        raise KeyError(gene + " not in current as stucture dict ignoring cluster ")
    
    if not interval.strand == location_dict[gene]['strand']:
        raise ValueError("strands not the same, there is some issue with gene annotations")
    
    total_length = float(sum(stop - start for start, stop in location_dict[gene]['regions']))
    
    running_length = 0
    
    #reverses list if negative strand so running count can work
    #if location_dict[gene]['strand'] == "-":
    #    location_dict[gene]['regions'].reverse()
        
    for start, stop in location_dict[gene]['regions']:
        length = float(stop - start) 

        if peak_center >= int(start) and peak_center <= int(stop):
            if interval.strand == "+":
                total_location = running_length + (peak_center - start)
                total_fraction = np.round((total_location / total_length), 3)

                individual_fraction = (peak_center - start) / length

            elif interval.strand == "-":
                total_location = running_length + (stop - peak_center)
                total_fraction = total_location / total_length

                individual_fraction = (stop - peak_center) / length
                
            else:
                raise ValueError("Strand not correct strand is %s" % interval.strand)

            #probably not nessessary
            if total_fraction < 0 or total_fraction > 1:
                raise ValueError("total_fraction is bad: %f, gene %s, total_length: %s, total_location: %s" % (total_fraction, 
                                                                                                               gene, 
                                                                                                              total_length,
                                                                                                               total_location))
            return individual_fraction, total_fraction
        
        running_length += length
        
    return None, None #clusters fall outside of regions integrated 

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
    return { name : bedtool.each(name_to_chrom).saveas() for name, bedtool in feature_dict.items()}

def convert_to_mrna(feature_dict, exon_dict):
    
    """
    
    converts transcript dict into mRNA locations given a dict of exons
    
    feature_dict - generated from convert to transcript, dict of bedfiles 
    exon_dict - dict of { genes : [exons] }
    
    """

    return { name : bedtool.each(convert_to_mRNA_position, exon_dict).filter(lambda x: x.chrom != "none").saveas() for name, bedtool in feature_dict.items()}

def invert_neg(interval):
    interval[-1] = str(int(interval[-1]) * -1)
    return interval

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
    for interval in bedtool:
        try: 
            #will need to redefine this to use intervals
            exon, total = RNA_position_interval(interval, region_dict)
  
            if total is not None:
                total_distributions.append(total)
                exon_distributions.append(exon)
            else:
                num_missed.append(interval)
        except Exception as e:
            print e
            num_errors.append(interval)

    return {'individual' : exon_distributions, 'total' : total_distributions, 
            'errors' : num_errors, 'missed' : num_missed}

def calculate_peak_locations(bedtool, genes, regions, features):
    
    """
    
    bedtools - bedtools of peaks, both  (only works with clipper defined peaks)
    genes - as_structure dictionary 
    db - database of gff file from gffutils
    
    Counts peak locations returns mRNA, pre mRNA, exon and intron positions
    and the nearest type of exon (CE SE ect...)
    
    cluster_regions - bedtool describing locations of peaks / clusters
    
    """
    

    mRNA_positions = []
    premRNA_positions = []
    intron_positions = []
    exon_positions = []


    exon_dict = generate_region_dict(regions['exons'])
    
    bed_center = bedtool.each(small_peaks).saveas()
    beds_center_transcripts = bed_center.each(name_to_chrom).saveas()
    beds_center_transcripts_mrna = beds_center_transcripts.each(convert_to_mRNA_position, exon_dict).filter(lambda x: x.chrom != "none").saveas()
    
    #####################
    #Here be feature code
    #####################


    features_transcript = convert_to_transcript(features)
    features_mrna = convert_to_mrna(features_transcript, exon_dict)
    
    #for pre-mrna
    #can't use list comprehensions because if file is empty than error is thrown...
    features_transcript_closest = {}
    for name, feature in features_transcript.items():
        try:
            features_transcript_closest[name] = {"dist" : beds_center_transcripts.closest(feature, s=True, D="b", t="first", id=True).filter(lambda x: x[-2] != ".").saveas()}
        except:
            features_transcript_closest[name] = {"dist" : None}
    #for mrna
    features_mrna_closest = {} 
    for name, feature in features_mrna.items():
        try:
            features_mrna_closest[name] = {"dist" : beds_center_transcripts_mrna.closest(feature, s=True, D="ref", t="first").filter(lambda x: x[-2] != ".").each(invert_neg).saveas()}
        except:
            features_mrna_closest[name] = {"dist" : None}
            
    #I should factor out this building logic...
    exons_dict = {}
    for gene in genes:
        exons_dict[gene] = {}
        exons_dict[gene]["strand"] = genes[gene]['strand']
        exons_dict[gene]['regions'] = [map(int, region.split("-")) for region in genes[gene]['exons'].values()]
        
    introns_dict = {}
    for gene in genes:
        introns_dict[gene] = {}
        introns_dict[gene]["strand"] = genes[gene]['strand']
        introns_dict[gene]['regions'] = [map(int, region.split("-")) for region in genes[gene]['introns'].values()]
    
    #might as well stitch them both together to get pre-mRNA distribution
    premRNA_dict = {}
    for gene in genes:
        premRNA_dict[gene] = {}
        premRNA_dict[gene]['regions'] = []
        premRNA_dict[gene]['strand'] = genes[gene]['strand']
        for exon, intron in izip_longest(genes[gene]['exons'].values(), genes[gene]['introns'].values()):
            premRNA_dict[gene]['regions'].append(map(int, exon.split("-")))
            if intron is not None:
                premRNA_dict[gene]['regions'].append(map(int, intron.split("-")))

#counts nearest exon to peak and gets RNA
#gets rna positon for every line as well
    #make nearest types bed file (custom format for closest bed parsing

    #need to generate CDS lists
    #need to make nearest type function
    
    distributions = {}
    for name, region in regions.items(): 
        region_dict = generate_region_dict(region)
        distributions[name] = get_distributions(bedtool, region_dict)
    
    
    for interval in bedtool:
        try:
            exon_frac, mRNA_frac = RNA_position(interval, exons_dict)
        except:
            exon_frac, mRNA_frac = None, None
        try:
            intron_frac, no_meaning = RNA_position(interval, introns_dict)
        except:
            intron_frac, no_meaning = None, None
        try:
            total_frac, premRNA_frac = RNA_position(interval, premRNA_dict)
        except:
            total_frac, premRNA_frac = None, None

        if mRNA_frac is not None:
            mRNA_positions.append(mRNA_frac)
        if exon_frac is not None:
            exon_positions.append(exon_frac)
        if premRNA_frac is not None:
            premRNA_positions.append(premRNA_frac)
        if intron_frac is not None:
            intron_positions.append(intron_frac)


    #get closest features to peaks
    bedtool_list = []
    for name, gene in genes.items():
        for exon, type in zip(gene['exons'].values(), gene['types'].values()):
            start, stop = map(int, exon.split("-"))
            bedtool_list.append([gene['chrom'], start, stop, name, 0, gene['strand'], type])
    
    feature_tool = pybedtools.BedTool(bedtool_list)
    
    #gets and counts closest features to each peak
    #this doesn't take into account badly placed peaks, we'll know if this is a problem if there a lot of ignoring error messages
    types = Counter([interval[-1] for interval in bedtool.closest(feature_tool, s=True)])
    
    return types, premRNA_positions, mRNA_positions, exon_positions, intron_positions, features_transcript_closest, features_mrna_closest, distributions


#here we start a small module for getting read densities:

def get_peak_wiggle(tool, data_pos, data_neg, size = 1000, num_peaks=100):
    
    """
    
    returns matrix of (n,25) of wiggle and peak information
    tool - bedtool contains location of peaks
    data - bigwig file handle - contains bw data for the peaks

    """
    
    results = []
    index = []
    for line in tool[:num_peaks]:
        peak_center = line.start
        pad_amt = size / 2
        try:
            if line.strand == "+":
                result = data_pos.get_as_array(line.chrom, peak_center - pad_amt , peak_center + pad_amt)
            if line.strand == "-":
                result = data_neg.get_as_array(line.chrom, peak_center - pad_amt , peak_center + pad_amt)[::-1]
                
            result = [0 if np.isnan(x) else x for x in result]
        
            results.append(result)
            index.append(line.name)
            
        except:
            pass
    return pd.DataFrame(results, index=index)

def cluser_peaks(bedtool, bw_pos, bw_neg, k=16):
    
    """
    
    Given a bedtool of peaks, positive and negative bigwig files
    gets read densities around peaks, normalizes them and outputs clusters and dataframe 
    
    """
    bed_center = bedtool.each(small_peaks).saveas()
    bedtool_df = get_peak_wiggle(bed_center, bw_pos, bw_neg, num_peaks=1000000)

    bedtool_df_mag_normalized = bedtool_df.div(bedtool_df.sum(axis=1), axis=0).fillna(0)
    bedtool_df_mag_normalized[bedtool_df_mag_normalized > .1] = .1
    
    data = np.array(bedtool_df_mag_normalized)
    
    #fixes bad test case
    if len(data) < k:
        k = len(data) / 2
        
    kmeansClassifier = KMeans(k)
    classes = kmeansClassifier.fit_predict(data)
    return bedtool_df_mag_normalized, classes

def run_homer(foreground, background, k = list([5,6,7,8,9]), outloc = os.getcwd()):
    
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
            
        print "Homer Finished, output here: %s" %(outloc)
    except OSError:
        print "Homer not installed, ignoring motif generation, install homer for this to work"  
        raise   
    
def write_seqs(outfile, bedtool_list):
    
    """
    
    outputs bedtools file to another file 
    
    Combines all sequences because its the punative background that gets used by other analysies
    """
    
    #TODO refactor to just the bedtools save function
    f = open(outfile, "w")
    for bedtool in bedtool_list:
        f.write(open(bedtool.seqfn).read())
    f.close()
    return

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


def fa_file(filename, region = None, directory=None, type = "real"):
    
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
            

def count_total_reads(bam, gene_definition):
    
    """
    
    Counts total reads in genes 
    
    bam: bedtool (reading in a bam file) all the reads in the original clip-seq experiment
    gene_definition: bedtool, defines locations of all genes to analize
    
    """
    
    return len(bam.intersect(gene_definition, u=True))

def count_reads_per_cluster(bedtool, bamtool):
    
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
    except: #can't count based on reported number, reverting to bedtools 
        cov = bamtool.coverage(bedtool, s=True, counts=True)
        for cluster in cov:
            reads_per_cluster.append(int(cluster[-1]))
            
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
                print "Ignoring: %s" % (region)
 
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
    phast_values = {"real" : {}, "rand" : {}}
    print "Fetching Phastcons Scores..."
    #phastcons values for all regions except "all"
    
    for region in regions: #skip "all" combine them later
        print "%s..." % (region)

        try: #gracefully fail if the region isn't represented
            phast_values["real"][region] = get_mean_phastcons(cluster_regions[region]['real'], 
                                                              phastcons_location, 
                                                              sample_size = sample_size)
            
            #can't concatanate zero length arrays, so prime it
            randPhast = np.array([])
            for rand_bed in cluster_regions[region]['rand'].values():                
                randPhast = np.concatenate((randPhast, get_mean_phastcons(rand_bed, 
                                                                          phastcons_location,  
                                                                          sample_size = sample_size)), axis=1)
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
        
        print "%s done" % (region)
    
    #assigns 'all' information
    dist_dict['all'] = {}
    dist_dict['all']['real'] = {'dist' : [], 'size' : 0}
    dist_dict['all']['rand'] = {'dist' : [], 'size' : 0}
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
            distance = distance * -1
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

def make_unique(feature):
        global uniq_count
        uniq_count +=1
        feature.name = feature.name + "_" + str(uniq_count)
        return feature
    
def main(options):
    
    """
    
    Runs all analysies 
    
    one thing to do is make graphs fail gracefully 
    
    """
    print "starting"
    #gets clusters in a bed tools + names species 
    clusters = os.path.basename(options.clusters)
    species = options.species
    
    #In case names aren't unique make them all unique
    global uniq_count
    uniq_count = 0
    clusters_bed = pybedtools.BedTool(options.clusters).each(make_unique).saveas()
    bamtool = pybedtools.BedTool(options.bam)
    bw_pos = BigWigFile(open(options.bw_pos))
    bw_neg = BigWigFile(open(options.bw_neg))
    
    #makes output file names 
    clusters = str.replace(clusters, ".BED", "")
    options.k = [int(x) for x in options.k]
    
    #all the different motifs that the user specifices 
    motifs = options.motifs
    outdir = options.outdir
    
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
    regions["all" ] = "All"
    regions["cds"] = "CDS"
    regions["three_prime_utrs"] = "3' UTR"
    regions["five_prime_utrs"] = "5' UTR"
    regions["proxintron500"] = "Proximal\nIntron"
    regions["distintron500"] = "Distal\nIntron"
    
    db = gffutils.FeatureDB(options.db)
    
    features = GenomicFeatures(species, db,  regions_dir=options.regions_location)
    genomic_regions = features.get_genomic_regions()
    features = features.get_feature_locations()
    
     #all catagory would break some analysies, create copy and remove it
    assigned_regions = regions.copy()
    del assigned_regions['all']
    
    genes_dict, genes_bed = parse_AS_STRUCTURE_dict(species, options.as_structure)
    
    if not options.assign:
        try:
            cluster_regions, region_sizes, genic_region_sizes = build_assigned_from_existing(assigned_dir, clusters, regions, options.nrand)
            
            print "I used a pre-assigned set output_file BED files... score!"
        except:
            
            print "I had problems retreiving region-assigned BED files from %s, i'll rebuild" % (assigned_dir)
            options.assign = True
    
    if options.assign:
        print "Assigning Clusters to Genic Regions"
        
       
        
        
        cluster_regions, region_sizes = assign_to_regions(clusters_bed,
                                                   clusters,
                                                   options.genome_location, 
                                                   options.regions_location, 
                                                   assigned_regions,
                                                   assigned_dir=assigned_dir,
                                                   fasta_dir=fasta_dir,
                                                   species=species, 
                                                   getseq=True, 
                                                   nrand=options.nrand,
                                                   )
        
        genic_region_sizes = count_genomic_region_sizes(options.regions_location, assigned_regions,
                                                            species)


            
    print "Counting reads in clusters...",
    reads_in_clusters, reads_per_cluster = count_reads_per_cluster(cluster_regions['all']['real'], bamtool)
    
    #might want to actually count genes_dict, not clusters...
    total_reads = count_total_reads(bamtool, genes_bed)
    region_read_counts = {}
    for region_name, cur_region in genomic_regions.items():
        region_read_counts[region_name] = count_total_reads(bamtool, cur_region)  
    
    #one stat is just generated here
    #generates cluster lengths (figure 3)
    cluster_lengths = bedlengths(cluster_regions['all']['real'])
    print "done"
    
    #also builds figure 10 (exon distances)
    genomic_types = count_genomic_types(genes_dict)
    
    #figures 5 and 6, builds pre-mrna, mrna exon and intron distributions
    (types, premRNA_positions, mRNA_positions, exon_positions, 
     intron_positions, features_transcript_closest, 
     features_mrna_closest, distributions) = calculate_peak_locations(cluster_regions['all']['real'], 
                                                                      genes_dict, 
                                                                      genomic_regions, 
                                                                      features)
    
    read_densities, classes = cluser_peaks(cluster_regions['all']['real'], bw_pos, bw_neg)
    
    #gtypes is total genomic content 
    #types is what clusters are
    #generates figure 10 (exon distances)
    
    type_count = [types["CE:"], types["SE:"], 
                  types["MXE:"], types["A5E:"], types["A3E:"]]
    
    genomic_type_count = [genomic_types["CE:"], genomic_types["SE:"], 
                          genomic_types["MXE:"], genomic_types["A5E:"], 
                          genomic_types["A3E:"]]    
    kmer_results = None
    if options.reMotif:
        kmer_results = calculate_kmer_diff(options.k, regions, clusters, fasta_dir)        
    
    if options.homer:
        calculate_homer_motifs(options.k, regions, clusters, fasta_dir, homerout)
    
    phast_values = None
    #loads phastcons values output_file generates them again
    if not options.rePhast:
        try:
            phast_values = pickle.load(open(os.path.join(misc_dir, "%s.phast.pickle" % (clusters))))
        except:
            options.rePhast = True

    print "starting phast"
    if options.runPhast:
        phast_values = calculate_phastcons(assigned_regions, cluster_regions, options.phastcons_location)
    print "ending phast"
    
    motif_distances = []
    try:
        if motifs:
            motif_distances = generate_motif_distances(cluster_regions, region_sizes, motifs, options.motif_location, options.species)
                    
    except:
        pass

    #save all analysies in a pickle dict
    out_dict = {}
    out_dict["region_sizes"] = region_sizes
    out_dict["reads_in_clusters"] = reads_in_clusters
    out_dict["reads_out_clusters"] = (total_reads - reads_in_clusters)
    out_dict["cluster_lengths"] = cluster_lengths
    out_dict["reads_per_cluster"] = reads_per_cluster
    out_dict["premRNA_positions"] = premRNA_positions
    out_dict["mRNA_positions"] = mRNA_positions
    out_dict["exon_positions"] = exon_positions
    out_dict["intron_positions"] = intron_positions
    out_dict["genic_region_sizes"] = genic_region_sizes
    out_dict["genomic_type_count"] = genomic_type_count
    out_dict["type_count"] = type_count
    out_dict["kmer_results"] = kmer_results
    out_dict["motifs"] = motifs
    out_dict["phast_values"] = phast_values
    out_dict["motif_distances"] = motif_distances
    for name, feature in features_transcript_closest.items():
        if feature['dist'] is not None:
            feature['dist'].saveas(clusters + "_" + name + "_transcript.bed")

    for name, feature in features_mrna_closest.items():
        if feature['dist'] is not None:
            feature['dist'].saveas(clusters + "_" + name + "_mrna.bed")
    out_dict['distributions'] = distributions
    out_dict['data'] = np.array(read_densities)
    out_dict['classes'] = classes
    out_dict['region_read_counts'] = region_read_counts
    out_file = open(os.path.join("%s.clip_analysis.pickle" %(clusters)), 'w')
    pickle.dump(out_dict, file=out_file)
    
    print "file saved"
    #build qc figure
    QCfig_params = [reads_in_clusters, (total_reads - reads_in_clusters), 
                    cluster_lengths, reads_per_cluster, distributions['genes']['total'], 
                    distributions['exons']['total'], distributions['exons']['individual'], 
                    distributions['introns']['individual'], genic_region_sizes, region_sizes, 
                    genomic_type_count, type_count, homerout, kmer_results, motifs, phast_values, 
                    regions, np.array(read_densities), classes]
    
    QCfig = CLIP_analysis_display.CLIP_QC_figure(*QCfig_params)
    distribution_fig = CLIP_analysis_display.plot_distributions(features_transcript_closest, features_mrna_closest, distributions)
    QCfig.savefig(os.path.join(outdir, clusters + ".QCfig." + options.extension))
    distribution_fig.savefig(os.path.join(outdir, clusters + ".DistFig." + options.extension))  
      
    #prints distance of clusters from various motifs in a different figure
    
    try:
        if motifs:
            motif_fig = CLIP_analysis_display.plot_motifs(motif_distances)
            motif_fig.savefig(clusters + ".motif_distribution." + options.extension)
    except:
        pass
    
    with open(options.metrics, 'w') as outfile:
        outfile.write("FRiP\n")
        outfile.write("\t".join([str(float(reads_in_clusters) / float(total_reads))]) + "\n")
    
    
    #fin
def call_main():  
    parser = OptionParser()
    
    parser.add_option("--clusters", dest="clusters", help="BED file of clusters", metavar="BED")
    parser.add_option("--bam", dest="bam", help="The bam file from the CLIP analysis")
    parser.add_option("--species", "-s", dest="species", help = "genome version")
    ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
    #parser.add_option("--reAssign", dest="assign", action="store_true", default=False, help="re-assign clusters, if not set it will re-use existing assigned clusters") 
    ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
    parser.add_option("--rePhast", dest="rePhast", action="store_true", default=False, help="re-calculate conservation, must have been done before") 
    parser.add_option("--runPhast", dest="runPhast", action="store_true", default=False, help="Run Phastcons ") 
    parser.add_option("--runMotif", dest="reMotif", action="store_true", default=False, help="Calculate Motif scores")
    parser.add_option("--runHomer", dest="homer", action="store_true", help="Runs homer", default=False)
    parser.add_option("--motifs", dest="motifs", action="append", help="Motifs to use (files of motifs give must exist in motif_directory directory)", default=[])
    parser.add_option("--k", dest="k", action="append", help="k-mer and homer motif ananlysis", default=[6])
    #parser.add_option("--conservation", dest="cons", help="Runs conservation (might not do anything)", action="store_true")
    #parser.add_option("--structure", dest="structure", help="also doesn't do anything gets structure maps", action="store_true")
    parser.add_option("--nrand", dest="nrand", default=3, help="selects number of times to randomly sample genome", type="int")
    parser.add_option("--outdir", "-o", dest="outdir", default=os.getcwd(), help="directory for output, default:cwd")
    ##Below here are critical files that always need to be referenced
    parser.add_option("--gff_db", dest="db", help="gff database from gffutils to generate annotations with")
    parser.add_option("--AS_Structure", dest="as_structure",  help="Location of AS_Structure directory (chromosme files should be inside)", default=None)
    parser.add_option("--genome_location", dest="genome_location", help="location of all.fa file for genome of interest", default=None)
    parser.add_option("--homer_path", dest="homer_path", action="append", help="path to homer, if not in default path", default=None)
    parser.add_option("--phastcons_location", dest="phastcons_location",  help="location of phastcons file", default=None)
    parser.add_option("--regions_location", dest="regions_location",  help="directory of genomic regions for a species", default=None)
    parser.add_option("--motif_directory", dest="motif_location",  help="directory of pre-computed motifs for analysis", default=os.getcwd())
    parser.add_option("--reAssign", dest="assign", action="store_true", default=False, help="re-assign clusters, if not set it will re-use existing assigned clusters")
    parser.add_option("--metrics", dest="metrics", default="CLIP_Analysis.metrics", help="file name to output metrics to")
    parser.add_option("--extension", dest="extension", default="svg", help="file extension to use (svg, png, pdf...)")
    parser.add_option("--bw_pos",  help="bigwig file, on the positive strand")
    parser.add_option("--bw_neg", help="bigwige file on the negative strand")

    (options, args) = parser.parse_args()
    
    #error checking
    if options.clusters is None or options.bam is None or options.species is None:
        parser.print_help()
        exit()
    main(options)

if __name__== "__main__":
    call_main()
