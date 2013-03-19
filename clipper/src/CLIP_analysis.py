
"""

Analizes CLIP data given a bed file and a bam file

Michael Lovci and Gabriel Pratt

"""
import pybedtools
import numpy as np
from optparse import OptionParser
import os
import pickle
import subprocess
from bx.bbi.bigwig_file import BigWigFile
from clipper.src import CLIP_Analysis_Display
from clipper.src.kmerdiff import kmer_diff
from collections import Counter
from random import sample

def intersection(A, B=None):
    
    """
    
    A : bedtool
    B : bedtool
    Returns A - B and A intersect B 
    A with b and returns everything in a but not b and everything in a but... ???
    
    """
    
    less = A.subtract(B, s=True) #without regions that overlap
    more = A.subtract(less, s=True) #only regions that overlap
    return less, more


def adjust_offsets(tool, offsets=None):
    
    """
    
    For finding motiff position relative to center of peaks
    Handles merging overlapping peaks, merge > bed12 -> bed6
    picks first offset from merged peak and assigns that to 
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

        #handles multiple different types of bed files (don't know why)
        if len(bed_line.fields) < 8:
            thick_start = 0
            thick_stop = 0
        else:
            thick_start = int(bed_line[6])
            thick_stop = int(bed_line[7])
        
        #the ; represents two merged locations 
        if ";" in bed_line.name and bed_line.name not in offsets:
            offset = offsets[bed_line.name.split(";")[0]]
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

def parse_AS_STRUCTURE_dict(species, working_dir):
    
    """
    
    Important return values:
    
    info a dict of gene info
    genes - bedtool of gene locations
    
    Parses out all important AS structure - see constructed dict in function
    for information on what is needed...
    
    also returns bed file of genes 
    Should refactor to be a real object, but I'm lazy right now...
    
    """
    

    
    print species
    if species == "hg19":
        chroms = [str(x) for x in range(1, 23)] #1-22
        chroms.append("X")
        chroms.append("Y")
    elif species == "mm9":
        chroms = [str(x) for x in range(1, 20)] #1-19
        chroms.append("X")
        chroms.append("Y")        
    elif species == "test":
        chroms = ["1"]

    info = {}
    bed_string = ""
    
    for chrom in chroms:

        as_file = os.path.join(working_dir, species + ".tx." + chrom + ".AS.STRUCTURE")
        
        f = open(as_file, "r")
        for line in f.readlines():
            if not line.startswith(">"):
                continue
            
            blank, gene, chrom, transcripts, d2, d3, d4, strand, number_of_exons, exonloc, intronloc, exonlen, intronlen, asType, locType = line.strip().split("\t")
            signstrand = "-"
            if int(strand) == 1:
                signstrand = "+"
            number_of_exons = int(number_of_exons)
            info[gene] = {}
            info[gene]['chrom'] = "chr" + str(chrom)
            info[gene]['strand'] = signstrand
            info[gene]['exons'] = {}
            info[gene]['introns'] = {}
            info[gene]['types'] = {}                        
            exons = exonloc.split("|")
            introns = intronloc.split("|")
            types = asType.split("|")
            info[gene]['numEx'] = number_of_exons
            info[gene]['mRNA_length'] = 0
            info[gene]['premRNA_length'] = 0
            tx_start = np.Inf
            tx_stop  = np.NINF
            for i, exon in enumerate(exons):
                if i == number_of_exons: #last exon is empty
                    continue
            
                info[gene]['exons'][i] = exon
                info[gene]['types'][i] = types[i]
                exstart, exstop = [int(x) for x in exon.split("-")]
                tx_start = min(exstart, tx_start)
                tx_stop = max(exstop, tx_stop)
                
                #there is an off by one bug in here somewhere, this if off 
                #from exon and intron lengths column
                info[gene]['mRNA_length'] += exstop-exstart+1
                info[gene]['premRNA_length'] += exstop-exstart+1                
            for i, intron in enumerate(introns):
                if i == number_of_exons-1: #only number_of_exons-1 introns
                    continue
                info[gene]['introns'][i] = intron
                intstart, intstop = [int(x) for x in intron.split("-")]
                info[gene]['premRNA_length'] += intstop-intstart+1
            info[gene]['tx_start'] = tx_start
            info[gene]['tx_stop'] = tx_stop
            bed_string += "\t".join([info[gene]['chrom'], str(info[gene]['tx_start']), str(info[gene]['tx_stop']), gene, "0", info[gene]['strand']]) + "\n"
            
            
    return info, pybedtools.BedTool(bed_string, from_string=True)

def count_genomic_region_sizes(regions_dir, species="hg19"):
    
    """
    
    Counts the genomic region sizes for the specified regions dir
    (should be changed to GTF)
    
    """
    
    genomic_region_sizes = {}
    regions = ["exon", "UTR3", "UTR5", "proxintron500", "distintron500"]
    
    for region in regions:
        region_tool = pybedtools.BedTool(os.path.join(regions_dir, region + "_" + species + "_frea_sorted.withscore"))
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
    speciesFA - the species fasta file
    regions_dir - the directory that has genomic regions already processed
    regions - list [str] regions to process, not used now but should be after refactoring
    species - str species to segment
    nrand - int number offsets times to shuffle for null hypothesis
    getseq - boolean gets the full sequence to store
    
    shuffling for background, should still be factored out
    
    """
    
    region_files = {}
    
    #constructs bed tools for each region
    #TODO fix names
    bedtracks = {}
    #bedtracks = build_genomic_regions(pybedtools.BedTool(regions_dir))
    #regions = bedtracks.keys()
    regions = ["exon", "UTR3", "UTR5", "proxintron500", "distintron500"]
    for region in regions:
        region_files[region] = os.path.join(regions_dir, region + "_" + species + "_frea_sorted.withscore")
        bedtracks[region] = pybedtools.BedTool(region_files[region])
              
    #creates the basics of bed dict
    bed_dict = {}
    bed_dict['all'] = {}
    bed_dict['all']['rand'] = {}

    #can't append None to string so hack a null bed tool
    bed_dict['all']['real'] = pybedtools.BedTool("", from_string = True)
    
    offsets = get_offsets_bed12(tool)
    tool = tool.merge(s=True, nms=True, scores="max")
    remaining_clusters = adjust_offsets(tool, offsets)
    
    print "There are a total %d clusters I'll examine" % (len(tool))
    sizes = {}
    
    for region in regions:
        output_file = clusters + "." + region + ".real.BED"
        
        #check with mike about this
        remaining_clusters, overlapping = intersection(remaining_clusters, 
                                                       B = bedtracks[region])  
        
        #if for some reason there isn't a peak in the region skip it
        if len(overlapping) == 0:
            print "ignoring %s " % (region)
            continue
        
        #sets up bed dict for this region
        bed_dict[region] = {}
        bed_dict[region]['rand'] = {}
        
        #hack to allow saving of filter results
        bed_dict[region]['real'] = overlapping.filter(is_valid_bed12).sort().saveas(os.path.join(assigned_dir, output_file))
        remaining_clusters = pybedtools.BedTool(str(remaining_clusters.filter(is_valid_bed12)), from_string=True)
        
        no_overlapping_count = len(remaining_clusters)
        overlapping_count = len(bed_dict[region]['real'])
        print "For region: %s I found %d that overlap and %d that don't" % (region, overlapping_count, no_overlapping_count)
              
        bed_dict['all']['real'] = pybedtools.BedTool(str(bed_dict['all']['real']) + str(bed_dict[region]['real']), from_string=True)    
        
        #this should be factored out
        sizes[region] = bed_dict[region]['real'].total_coverage()

        bed_dict[region]['real'].sequence(fi=speciesFA, s=True)

        #saves offsets so after shuffling the offsets can be readjusted    
        offset_dict = get_offsets_bed12(overlapping)
    
        #TODO refactor to different function
        for i in range(nrand):
            output_file = clusters + "." + region + ".rand." + str(i) + ".BED"
            bed_dict['all']['rand'][i] = pybedtools.BedTool("", from_string=True)
            
            #for each region shuffles all peaks in that region around the region 
            #then pulls out sequences if requested 
            random_intervals = bed_dict[region]['real'].shuffle(genome=species, incl=bedtracks[region].fn).sort()

            #shuffling doesn't change offsets so we adjust bed 11 and 12 lines here to correct 
            random_intervals = adjust_offsets(random_intervals, offset_dict)
            
            #saves intervals, makes the all interval by appending all regions together
            bed_dict[region]['rand'][i] = random_intervals.saveas(os.path.join(assigned_dir, output_file))
            bed_dict['all']['rand'][i] = pybedtools.BedTool(str(bed_dict['all']['rand'][i]) + str(bed_dict[region]['rand'][i]), from_string=True)


    print "After assigning, I\'m left with %d un-categorized regions" %(len(remaining_clusters))
    try:
        #TODO this is wrong that only gets the non-overlapping somethings...
        bed_dict['uncatagorized'] = remaining_clusters.sort()
    except:
        pass

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
    for region in regions:
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


    
def RNA_position(bedline, as_structure_dict):
    
    """
    
    makes mrna and pre-mrna position figure 
    bedline - single bedline
    as_structure_dict - from build AS structure dict 
    
    Might be able to use my ribo-seq stuff for genic -> transcriptomic location conversion
    
    """
    
    mRNA_pos = 0
    pre_pos = 0
    exon_frac = None
    intron_frac = None
    nearest_type= None
    chrom, start, stop, name, score, strand, thickstart, thickstop = str(bedline).strip().split("\t")
    thickstart, thickstop = [int(x) for x in (thickstart, thickstop)]
    position = int((thickstart + thickstop)/2)
    

    
    try:
        gene, n, reads = name.split("_")
    except:
        #takes first gene if there are multiple overlapping 
        gene, n, reads = name.split(";")[0].split("_")
    
    if gene not in as_structure_dict:
        raise KeyError(gene + " not in current as stucture dict ignoring cluster ")
    
    for exN in range(as_structure_dict[gene]['numEx']):
        exstart, exstop = map(int, as_structure_dict[gene]['exons'][exN].split("-"))
        exlen = exstop-exstart+1
        if position >= exstart and position <= exstop:
            if strand == "+":
                mRNA_pos += position - exstart
                exon_frac = np.round(((position-exstart)/float(exlen)), 3)
                pre_pos += position - exstart                
            else:
                mRNA_pos += exstop - position
                exon_frac = np.round(((exstop - position)/float(exlen)), 3)
                pre_pos += exstop - position

            mRNA_frac = np.round((mRNA_pos/float(as_structure_dict[gene]['mRNA_length'])), 3)
            premRNA_frac = np.round((pre_pos/float(as_structure_dict[gene]['premRNA_length'])), 3)

            if mRNA_frac < 0 or mRNA_frac > 1:
                print "mRNA_frac is bad: %f, gene %s" %(mRNA_frac, gene)

            if premRNA_frac < 0 or premRNA_frac > 1:
                print "premRNA_frac is bad: %f, gene %s" %(premRNA_frac, gene)

            nearest_type = as_structure_dict[gene]['types'][exN]
            return mRNA_frac, premRNA_frac, exon_frac, None, nearest_type
        else:
            mRNA_pos += exlen
            pre_pos += exlen
            if exN < as_structure_dict[gene]['numEx']: #there are only exN - 1 introns
                intstart, intstop = map(int, as_structure_dict[gene]['introns'][exN].split("-"))
                intlen = intstop-intstart+1
                if position >= intstart and position <= intstop:
                    mRNA_pos = None                
                    if strand == "+":
                        pre_pos += position - intstart
                        intron_frac = np.round(((position-intstart)/float(intlen)), 3)
                    else:
                        pre_pos += intstop - position
                        intron_frac = np.round(((intstop - position)/float(intlen)), 3)
                    premRNA_frac = np.round((pre_pos/float(as_structure_dict[gene]['premRNA_length'])), 3)
                    if premRNA_frac > 0.5:
                        nearest_type=as_structure_dict[gene]['types'][exN+1]
                    else:
                        nearest_type=as_structure_dict[gene]['types'][exN]


                    if premRNA_frac < 0 or premRNA_frac > 1:
                        print "premRNA_frac is bad: %f, gene %s" %(premRNA_frac, gene)


            
                    return None, premRNA_frac, None, intron_frac, nearest_type
                else:
                    pre_pos += intlen
                    
    raise ValueError("Clusters fall outside gene")

def calculate_peak_locations(clusters, genes):
    
    """
    
    Counts peak locations returns mRNA, pre mRNA, exon and intron positions
    and the nearest type of exon (CE SE ect...)
    
    cluster_regions - bedtool describing locations of peaks / clusters
    
    """
    types = Counter()
    print "locating clusters within genes"
##This should be abstracted to some sort output_file list or something...
#figures 5 and 6, builds pre-mrna, mrna exon and intron distributions
    mRNA_positions = []
    premRNA_positions = []
    intron_positions = []
    exon_positions = []
#counts nearest exon to peak and gets RNA
#gets rna positon for every line as well
    for line in clusters:
        try:
            mRNA_frac, premRNA_frac, exon_frac, intron_frac, nearest_type = RNA_position(line, genes)
        
            if mRNA_frac is not None:
                mRNA_positions.append(mRNA_frac)
            if exon_frac is not None:
                exon_positions.append(exon_frac)
            if premRNA_frac is not None:
                premRNA_positions.append(premRNA_frac)
            if intron_frac is not None:
                intron_positions.append(intron_frac)
            if nearest_type is not None:
                types[nearest_type] += 1
        except KeyError:
            print "ignoring: " + str(line)
    
    return types, premRNA_positions, mRNA_positions, exon_positions, intron_positions

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
        result = subprocess.call(["findMotifs.pl", foreground, "fasta", outloc, "-p", "4", "-rna", "-S", "20", "-len", k, "-fasta", background])
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
 
def get_mean_phastcons(bedtool, phastcons_location, sample_size = 10000):
    
    """
    
    Get means phastcons scores for all intervals in a bed tool
    bedtool - bedtool to extract data from
    phastcons_location - location of phastcons file
    
    """
    
    f = open(phastcons_location, 'r')
    bw = BigWigFile(file=f)

    #if bedtool
    data = []
    for bedline in sample(bedtool, min(len(bedtool), sample_size)):
              
        conservation_values = bw.get_as_array(bedline.chrom, bedline.start, bedline.stop)
        
        if len(conservation_values) > 0:
            mean_phastcons = np.mean(conservation_values)
        else:
            mean_phastcons = 0
        data.append(mean_phastcons)
        
    return data


def fa_file(filename, region = None, directory=None, type = "real"):
    
    """
    
    Formats fasta file, and returns fully qualified name
    Checks if a fasta file exists returns the file attaced to a region
    or something 
    
    """
    
    if not os.path.exists(directory):
        raise Exception
    
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

def count_reads_per_cluster(bedtool):
    
    """
    
    Counts the number of reads in each cluster
    
    bedtool: bedtool containing the specially formatted CLIPper bedlines
    
    returns list(int) each index being a read in the cluster
    
    """
    #generates data for figure 1 and 2
    #gets reads in clusters (figure 1)
    #gets reads per cluster (figure 2)

    reads_per_cluster = []
    for cluster in bedtool:
        
        #handles merged clusters (taking the reads from only the first one)
        first_cluster = cluster.name.split(";")[0]
        number_reads_in_peak = first_cluster.split("_")[2]
        reads_per_cluster.append(int(number_reads_in_peak))
 
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
                kmer_results[region][k] = kmer_diff(real_fa, rand_fa, k)
            except IOError as e:
                print e
                print "Ignoring: %s" % (region)
                
    return kmer_results
 
def calculate_homer_motifs(kmer_list, regions, run_homer, clusters, fasta_dir, homerout):
    
    """
    
    Calculates motiffs for homer
    
    kmer_list - list[int] different kmers to examine
    regions   - list[str] regions to compute kmers on
    run_homer - boolean, run homer or not
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
        if run_homer is True:
            region_homer_out = os.path.join(homerout, region)
            run_homer(real_fa, rand_fa, kmer_list, outloc=region_homer_out)

def calculate_phastcons(regions, cluster_regions, phastcons_location):
    
    """
    
    Generates phastcons scores for all clusters
    
    regions - list[str] regions to analize
    cluster_regions - special dict that contains all random and normal cluster info
    phastcons_location - str location of phastcons file
    
    
    """
    
    
    ###conservation --should use multiprocessing to speed this part up!
#start output_file conservation logic, very slow...
    phast_values = {"real" : {}, "rand" : {}}
    print "Fetching Phastcons Scores..."
    #phastcons values for all regions except "all"
    for region in regions[1:]: #skip "all" combine them later
        print "%s..." % (region)
        try: #gracefully fail if the region isn't represented
            #think about collecting random sample to study...
            #samplesize = 1000
            print "getting real..." #gets phastcons values real regions
            phast_values["real"][region] = get_mean_phastcons(cluster_regions[region]['real'], phastcons_location)
            
            #can't concatanate zero length arrays, so prime it
            randPhast = np.array([])
            for i in range(len(cluster_regions[region]['rand'])):
                randPhast = np.concatenate((randPhast, get_mean_phastcons(cluster_regions[region]['rand'][i], phastcons_location)), axis=1)
            phast_values["rand"][region] = randPhast
        except KeyError as e:
            print "ignoring: ", e
            continue
    
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

def main(options):
    
    """
    
    Runs all analysies 
    
    one thing to do is make graphs fail gracefully 
    
    """
    print "starting"
    #gets clusters in a bed tools + names species 
    clusters = options.clusters
    species = options.species
    clusters_bed = pybedtools.BedTool(clusters)

    #makes output file names 
    clusters = str.replace(clusters, ".BED", "")
    options.k = [int(x) for x in options.k]
    
    #all the different motifs that the user specifices 
    motifs = list(options.motif)
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

    regions = (["all", "exon", "UTR3", "UTR5", "proxintron500", "distintron500"])    
    
    #Not quite sure whats going on here, but its one logical block
    #either reassigns clusters to genic regions or reads from already
    #made assigned lists
    genes_dict, genes_bed = parse_AS_STRUCTURE_dict(species, options.as_structure)
    
    if options.assign is False:
        try:
            cluster_regions, region_sizes, genic_region_sizes = build_assigned_from_existing(assigned_dir, clusters, regions, options.nrand)
            
            print "I used a pre-assigned set output_file BED files... score!"
        except:
            print "I had problems retreiving region-assigned BED files from %s, i'll rebuild" % (assigned_dir)
            options.assign = True
    
    #This is what I'm working on tomorrow
    #assign to reions / this big chunk that should get factoed out of main...
    if options.assign is True:
        print "Assigning Clusters to Genic Regions"
        #TODO what is region sizes?
        cluster_regions, region_sizes = assign_to_regions(clusters_bed,
                                                   clusters,
                                                   options.genome_location, 
                                                   options.regions_location, 
                                                   regions,
                                                   assigned_dir,
                                                   fasta_dir,
                                                   species=species, 
                                                   getseq=True, 
                                                   nrand=options.nrand,
                                                   )
        
        genic_region_sizes = count_genomic_region_sizes(options.regions_location, 
                                                            species)

    print "Counting reads in clusters...",
    reads_in_clusters, reads_per_cluster = count_reads_per_cluster(cluster_regions['all']['real'])
    
    #might want to actually count genes_dict, not clusters...
    total_reads = count_total_reads(pybedtools.BedTool(options.bam), genes_bed)
    
    #one stat is just generated here
    #generates cluster lengths (figure 3)
    cluster_lengths = bedlengths(cluster_regions['all']['real'])
    print "done"
    

    
    #also builds figure 10 (exon distances)
    
    genomic_types = count_genomic_types(genes_dict)
    types, premRNA_positions, mRNA_positions, exon_positions, intron_positions = calculate_peak_locations(cluster_regions['all']['real'], genes_dict)
        
    #gtypes is total genomic content 
    #types is what clusters are
    #generates figure 10 (exon distances)
    type_count = [types["CE:"], types["SE:"], 
                  types["MXE:"], types["A5E:"], types["A3E:"]]
    
    genomic_type_count = [genomic_types["CE:"], genomic_types["SE:"], 
                          genomic_types["MXE:"], genomic_types["A5E:"], 
                          genomic_types["A3E:"]]    
    kmer_results = []
    if options.reMotif is True:
        kmer_results = calculate_kmer_diff(options.k, regions, clusters, fasta_dir)        
        calculate_homer_motifs(options.k, regions, options.homer, clusters, fasta_dir, homerout)
    
    phast_values = []
    #loads phastcons values output_file generates them again
    if not options.rePhast:
        try:
            phast_values = pickle.load(open(os.path.join(misc_dir, "%s.phast.pickle" % (clusters))))
        except:
            options.rePhast = True

    if options.rePhast and options.runPhast:
        phast_values = calculate_phastcons(regions, cluster_regions, options.phastcons_location)
    
    #build qc figure
    QCfig_params = [reads_in_clusters, (total_reads - reads_in_clusters), 
                    cluster_lengths, reads_per_cluster, premRNA_positions, 
                    mRNA_positions, exon_positions, intron_positions, 
                    genic_region_sizes, region_sizes, genomic_type_count,
                     type_count, homerout, kmer_results, motifs, phast_values]
    
    QCfig = CLIP_Analysis_Display.CLIP_QC_figure(*QCfig_params)
    fn = clusters + ".QCfig.pdf"
    outFig = os.path.join(outdir, fn)
    QCfig.savefig(outFig)
        
    #prints distance of clusters from various motifs in a different figure
    try:
        if motifs is not None:
            motif_distances = generate_motif_distances(cluster_regions, region_sizes, motifs, options.motif_location, options.species)
            motif_fig = CLIP_Analysis_Display.plot_motifs(motif_distances)
            motif_fig.savefig(clusters + ".motif_distribution.pdf")
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

    out_file = open(os.path.join(assigned_dir, "%s.pickle" %(clusters)), 'w')
    pickle.dump(out_dict, file=out_file)

    #fin
if __name__== "__main__":
    parser = OptionParser()
    
    parser.add_option("--clusters", dest="clusters", help="BED file of clusters", metavar="BED")
    parser.add_option("--bam", dest="bam", help="The bam file from the CLIP analysis")
    parser.add_option("--species", "-s", dest="species", help = "genome version")
    ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
    parser.add_option("--reAssign", dest="assign", action="store_true", default=False, help="re-assign clusters, if not set it will re-use existing assigned clusters") 
    ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
    parser.add_option("--rePhast", dest="rePhast", action="store_true", default=False, help="re-calculate conservation, must have been done before") 
    parser.add_option("--runPhast", dest="runPhast", action="store_true", default=False, help="Run Phastcons ") 

    parser.add_option("--runMotif", dest="reMotif", action="store_true", default=False, help="Calculate Motif scores")
    parser.add_option("--motif", dest="motif", action="append", help="Files of motif locations", default=None)
    parser.add_option("--homer", dest="homer", action="store_true", help="Runs homer", default=False)
    parser.add_option("--k", dest="k", action="append", help="k-mer and homer motif ananlysis", default=[6])
    parser.add_option("--conservation", dest="cons", help="Runs conservation (might not do anything)", action="store_true")
    parser.add_option("--structure", dest="structure", help="also doesn't do anything gets structure maps", action="store_true")
    parser.add_option("--nrand", dest="nrand", default=3, help="selects number of times to randomly sample genome", type="int")
    parser.add_option("--outdir", "-o", dest="outdir", default=os.getcwd(), help="directory for output, default:cwd")
    ##Below here are critical files that always need to be referenced
    parser.add_option("--AS_Structure", dest="as_structure",  help="Location of AS_Structure directory (chromosme files should be inside)", default=None)
    parser.add_option("--genome_location", dest="genome_location", help="location of all.fa file for genome of interest", default=None)
    parser.add_option("--homer_path", dest="homer_path", action="append", help="path to homer, if not in default path", default=None)
    parser.add_option("--phastcons_location", dest="phastcons_location",  help="location of phastcons file", default=None)
    parser.add_option("--regions_location", dest="regions_location",  help="directory of genomic regions for a species", default=None)
    parser.add_option("--motif_directory", dest="motif_location",  help="directory of pre-computed motifs for analysis", default=os.getcwd())
   
    (options, args) = parser.parse_args()
    
    #error checking
    if options.clusters is None or options.bam is None or options.species is None:
        parser.print_help()
        exit()
        
    main(options)

