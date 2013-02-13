
"""

Analizes CLIP data given a bed file and a bam file

Michael Lovci and Gabriel Pratt

"""
import pybedtools
import numpy as np
from optparse import OptionParser
import os
import sys
import pickle
import random
#from subprocess import Popen, call, PIPE
import subprocess
from bx.bbi.bigwig_file import BigWigFile
import pysam
import clipper.src.CLIP_Analysis_Display
import pylab
from kmerdiff import kmer_diff

def intersection(A, B=None):
    
    """
    
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
        #TODO need to refactor to avoid this try junk
        #TODO don't change bedtool, use its map function instead
        try:
            chrom, start, stop, name, score, strand, thick_start, thick_stop = str(bed_line).strip().split("\t")
        
        except:
            chrom, start, stop, name, score, strand = str(bed_line).strip().split("\t")
            thick_start = 0
            thick_stop = 0
            
        start, stop, thick_start, thick_stop = [int(x) for x in (start, stop, thick_start, thick_stop)]
        
        #the ; represents two merged locations 
        if ";" in name and name not in offsets:
            o = name.split(";")
            offset = offsets[o[0]]
        else:
            offset = offsets[name]
        
        if strand == "+":
            thick_start = start + offset
            thick_stop = thick_start + 4
        else:
            thick_stop = stop - offset
            thick_start = thick_stop -4
        clusters.append("\t".join([str(x) for x in (chrom, start, stop, name, score, strand, thick_start, thick_stop)]))
    
    return pybedtools.BedTool("\n".join(clusters), from_string=True)
    




def build_AS_STRUCTURE_dict(species, working_dir):
    
    #TODO make more geneticreplace with clipper AS structure dict method
    #totally revamp convert to gtf files
    #these have to go
    """
    
    Important return values:
    
    Parses out all important AS structure - see constructed dict in function
    for information on what is needed...
    
    GTypes - number of each exon type (used for getting background number
    of each type of exon) (can also probably be removed)
    
    """
    
    if species == "hg19":
        chrs = [str(x) for x in range(1,23)] #1-22
        chrs.append("X")
        chrs.append("Y")
    elif species == "mm9":
        chrs = [str(x) for x in range(1,20)] #1-19
        chrs.append("X")
        chrs.append("Y")        

    info = dict()
    Gtypes = dict()
    for chr in chrs:

        ASfile = os.path.join(working_dir, species + ".tx." + chr + ".AS.STRUCTURE")
        
        f = open(ASfile, "r")
        for line in f.readlines():
            if not line.startswith(">"):
                continue
            
            blank, gene, chr, transcripts, d2, d3, d4, strand, numex, exonloc, intronloc, exonlen, intronlen, asType, locType = line.strip().split("\t")
            signstrand = "-"
            if int(strand) == 1:
                signstrand = "+"
            numex = int(numex)
            info[gene] = {}
            info[gene]['strand'] = signstrand
            info[gene]['exons'] = {}
            info[gene]['introns'] = {}
            info[gene]['types'] = {}                        
            exons = exonloc.split("|")
            introns = intronloc.split("|")
            types = asType.split("|")
            info[gene]['numEx'] = numex
            info[gene]['mRNA_length'] = 0
            info[gene]['premRNA_length'] = 0
            tx_start = 1000000000000
            tx_stop = -1000000000000
            for i, exon in enumerate(exons):
                if i == numex: #last exon is empty
                    continue
                
                #TODO refactor to default dict
                try:
                    Gtypes[types[i]] +=1
                except:
                    Gtypes[types[i]] = 1
                info[gene]['exons'][i] = exon
                info[gene]['types'][i] = types[i]
                exstart, exstop = map(int, exon.split("-"))
                tx_start = min(exstart, tx_start)
                tx_stop = max(exstop, tx_stop)
                info[gene]['mRNA_length'] += exstop-exstart+1
                info[gene]['premRNA_length'] += exstop-exstart+1                
            for i, intron in enumerate(introns):
                if i == numex-1: #only numex-1 introns
                    continue
                info[gene]['introns'][i] = intron
                intstart, intstop = map(int, intron.split("-"))
                info[gene]['premRNA_length'] += intstop-intstart+1
            info[gene]['tx_start'] = tx_start
            info[gene]['tx_stop'] = tx_stop
    return info, Gtypes

def assign_to_regions(tool, speciesFA, regions_dir, species="hg19", nrand = 3, getseq=False):
    
    """
    
    Assigns each cluster to a genic region
    
    tool - a bed tool (each line represnting a cluster)
    species - the species to assign to
    nrand - int number offsets times to shuffle for null hypothesis
    getseq - boolean gets the full sequence to store
    
    shuffling for background, should still be factored out
    
    """
      
    UTR3File = os.path.join(regions_dir, "UTR3_" + species + "_frea_sorted.withscore")
    UTR3 = pybedtools.BedTool(UTR3File)
    G_UTR3_size = UTR3.total_coverage()
    
    UTR5File = os.path.join(regions_dir, "UTR5_" + species + "_frea_sorted.withscore")
    UTR5 = pybedtools.BedTool(UTR5File)
    G_UTR5_size = UTR5.total_coverage()
    
    exonFile = os.path.join(regions_dir, "exon_" + species + "_frea_sorted.withscore")
    exon = pybedtools.BedTool(exonFile)
    G_exon_size = exon.total_coverage()    
    
    proxintronFile = os.path.join(regions_dir, "proxintron500_" + species + "_frea_sorted.withscore")
    proxintron = pybedtools.BedTool(proxintronFile)
    G_proxintron_size = proxintron.total_coverage()        
    
    distintronFile = os.path.join(regions_dir, "distintron500_" + species + "_frea_sorted.withscore")
    distintron = pybedtools.BedTool(distintronFile)
    G_distintron_size = distintron.total_coverage()    
    
    all_regions = (["exon", "UTR3", "UTR5", "proxintron", "distintron"])
    all_bedtracks = ([exon, UTR3, UTR5, proxintron, distintron])
    all_regionfiles = ([exonFile, UTR3File, UTR5File, proxintronFile, distintronFile])
    
    bed_dict = {}
    last = 0
    offsets = get_offsets_bed12(tool)
    tool =  tool.merge(s=True, nms=True, scores="max")
    tool = adjust_offsets(tool, offsets)
    Gsizes = [G_exon_size, G_UTR3_size, G_UTR5_size, G_proxintron_size, G_distintron_size]
    sizes = [0, 0, 0, 0, 0]

    print "There are a total offsets %d clusters I'll examine" %(len(tool))

    for j, region in enumerate(all_regions):
        try:
            #portions offsets the regions that overlap with a genic region
            no_overlapping, only_overlapping = intersection(tool, B = all_bedtracks[j])  
        except Exception as e:
            print e
            "prnt intersection failed"
            continue
        only_overlapping = pybedtools.BedTool(str(only_overlapping.filter(eliminate_invalid_bed12)), from_string=True)
        no_overlapping = pybedtools.BedTool(str(no_overlapping.filter(eliminate_invalid_bed12)), from_string=True)
        noCT = len(no_overlapping)
        onlyCT = len(only_overlapping)
        last += onlyCT
        both = last + noCT
        print "For region: %s I found %d that overlap and %d that don't, making a total offsets %d" %(region, onlyCT, noCT, both)
        bed_dict[region]= {}
        bed_dict[region]['real'] = None
        try:
            bed_dict[region]['real'] = only_overlapping.sort()
        except:
            continue
        try:
            bed_dict['all']['real'] = pybedtools.BedTool(str(bed_dict['all']['real']) + str(bed_dict[region]['real']), from_string=True)
        except:
            bed_dict['all']= {}
            bed_dict['all']['rand']= {}
            bed_dict['all']['real'] = bed_dict[region]['real']
        try:
            sizes[j] = bed_dict[region]['real'].total_coverage()
        except:
            continue
        if getseq is True:
            try:
                bed_dict[region]['real'].sequence(fi=speciesFA, s=True)
            except:
                pass
        ofdic = get_offsets_bed12(only_overlapping)
        bed_dict[region]['rand'] = {}
        
        #TODO refactor to different function
        for i in range(nrand):
            try:
                #for each region shuffles all peaks in that region around the region 
                #then pulls out sequences if requested 
                
                random_intervals = bed_dict[region]['real'].shuffle(genome=species, incl=all_regionfiles[j]).sort()
            except:
                continue
            
            #shuffling doesn't change offsets so we adjust bed 11 and 12 lines here to correct 
            random_intervals = adjust_offsets(random_intervals, ofdic)
            if getseq is True:
                pass
                #print speciesFA
                #random_intervals.sequence(fi=speciesFA, s=True)
                
            bed_dict[region]['rand'][i] = random_intervals
            try:
                bed_dict['all']['rand'][i] = pybedtools.BedTool(str(bed_dict['all']['rand'][i]) + str(bed_dict[region]['rand'][i]), from_string=True)
            except:
                bed_dict['all']['rand'][i] = bed_dict[region]['rand'][i]

        tool = no_overlapping

    print "After assigning, I\'m left with %d un-categorized regions" %(len(tool))
    try:
        bed_dict['uncatagorized'] = tool.sort()
    except:
        pass
    bed_dict['all']['real'] = bed_dict['all']['real'].sort()
    bed_dict['all']['real'].sequence(fi=speciesFA, s=True)
    for i in range(nrand):
        bed_dict['all']['rand'][i] = bed_dict['all']['rand'][i].sort()
        bed_dict['all']['rand'][i].sequence(fi=speciesFA, s=True)
    return bed_dict, sizes, Gsizes

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

def eliminate_invalid_bed12(x):
    
    """
    
    Removes clusters that start in invalid locations for bed12 formatted files
    
    x - bedline
    
    returns either true or false if its a valid line or not
    
    """
    
    chr, start, stop, name, score, strand, tstart, tstop = str(x).split("\t")
    start, stop, tstart, tstop = map(int, (start, stop, tstart, tstop))
    if start < tstart and stop > tstop:
        return True
    else:
        return False

def get_offsets_bed12(tool):
    
    """
    
    gets offsets for each location in a bed12 file
    
    tool - bedtool
    
    returns offset for each location
    
    """
    
    of = {}
    for line in tool:
        chr, start, stop, name, score, strand, tstart, tstop = str(line).split("\t")
        start, stop, tstart, tstop = map(int, (start, stop, tstart, tstop))
        if strand == "+":
            offset = tstart -start 
        else:
            offset = stop - tstop
        of[name] = offset
    return of


    
def RNA_position(bedline, GI):
    
    """
    
    makes mrna and pre-mrna position figure 
    bedline - single bedline
    GI - from build AS structure dict 
    
    Might be able to use my ribo-seq stuff for genic -> transcriptomic location conversion
    
    """
    
    mRNA_pos = 0
    pre_pos = 0
    exon_frac = None
    intron_frac = None
    nearest_type=None
    chr, start, stop, name, score, strand, thickstart, thickstop = str(bedline).strip().split("\t")
    thickstart, thickstop = map(int, (thickstart, thickstop))
    position = int((thickstart + thickstop)/2)
    try:
        gene, n, reads = name.split("_")
    except:
        
        #takes first gene if there are multiple overlapping 
        gene, n, reads = name.split(";")[0].split("_")

    for exN in range(GI[gene]['numEx']):
        exstart, exstop = map(int, GI[gene]['exons'][exN].split("-"))
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

            mRNA_frac = np.round((mRNA_pos/float(GI[gene]['mRNA_length'])), 3)
            premRNA_frac = np.round((pre_pos/float(GI[gene]['premRNA_length'])), 3)

            if mRNA_frac < 0 or mRNA_frac > 1:
                print "mRNA_frac is bad: %f, gene %s" %(mRNA_frac, gene)

            if premRNA_frac < 0 or premRNA_frac > 1:
                print "premRNA_frac is bad: %f, gene %s" %(premRNA_frac, gene)

            nearest_type = GI[gene]['types'][exN]
            return mRNA_frac, premRNA_frac, exon_frac, None, nearest_type
        else:
            mRNA_pos += exlen
            pre_pos += exlen
            if exN < GI[gene]['numEx']: #there are only exN - 1 introns
                intstart, intstop = map(int, GI[gene]['introns'][exN].split("-"))
                intlen = intstop-intstart+1
                if position >= intstart and position <= intstop:
                    mRNA_pos = None                
                    if strand == "+":
                        pre_pos += position - intstart
                        intron_frac = np.round(((position-intstart)/float(intlen)), 3)
                    else:
                        pre_pos += intstop - position
                        intron_frac = np.round(((intstop - position)/float(intlen)), 3)
                    premRNA_frac = np.round((pre_pos/float(GI[gene]['premRNA_length'])), 3)
                    if premRNA_frac > 0.5:
                        nearest_type=GI[gene]['types'][exN+1]
                    else:
                        nearest_type=GI[gene]['types'][exN]


                    if premRNA_frac < 0 or premRNA_frac > 1:
                        print "premRNA_frac is bad: %f, gene %s" %(premRNA_frac, gene)


            
                    return None, premRNA_frac, None, intron_frac, nearest_type
                else:
                    pre_pos += intlen
    return None, None, None, None, None #shouldn't ever get here. this would mean that the clusters doesn't fall within the gene
                


def run_kmerdiff(clustersFA, backgroundFA, k=6, outfile=None):
    
    """
    
    Runs kmerdiff
    kmerdiff -computes zscores for each k-mer (consier re-writing in python / biopython)
    
    clustersFA - clusters in fasta format
    background - background locations
    k - size of kmer
    outfile - output
    
    returns sorted output
    
    """
    
    if outfile is None:
        outfile = clustersFA + ".k" + str(k) + ".kmerdiff"
    print "Running Kmer analysis"
    #TODO make more general 
    subprocess.call(["perl", (data_file + "/yeolab/Software/generalscripts/kmerdiff.pl"), "-file1", clustersFA, "-file2", backgroundFA, "-k", str(k), "-o", outfile])
    srtout = outfile + ".sort"
    srt = open(srtout, 'w')
    print "Sorting Results..."
    subprocess.call(["perl", (data_file + "/yeolab/Software/generalscripts/sortkmerdiff.pl"), outfile], stdout=srt)
    print "Kmer analysis done, output here: %s" %(srtout)
    srt.close()
    return srtout
    

def run_homer(foreground, background, k = list([5,6,7,8,9]), outloc = os.getcwd()):
    
    """
    
    runs homer with standard args
    
    output location is saved
    
    --make optional make work off locations and not fasta files 
    
    """
    #findMotifs.pl clusters.fa fasta outloc -nofacts p 4 -rna -S 10 -len 5,6,7,8,9 -noconvert -nogo -fasta background.fa
    ks = ",".join(map(str, k))
    print "starting Homer"
    subprocess.call(["findMotifs.pl", foreground, "fasta", outloc, "-p", "4", "-rna", "-S", "20", "-len", ks, "-fasta", background])
    print "Homer Finished, output here: %s" %(outloc)
    return


def write_seqs(outfile, bedtool_list):
    
    """
    
    outputs bedtools file to another file 
    
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
        p = get_phastcons(x, species="hg19")
        file.write("\t".join(map(str, [chr, i, (i+wsize-1), p])) + "\n")
        i += wsize
    file.close()
 
def get_phastcons(bedtool, phastcons_location, species=None, index=None, ):
    
    """
    
    Get phastcons scores for intervals in a bed tool
    
    """
    
    if species is None and index is None:
        print "Error, must select species or index"
    
    f = open(phastcons_location, 'r')
    bw = BigWigFile(file=f)

    try:
        
        #if its a line
        #for each line fetch bigwig values 
        type(bedtool)
        v = bedtool.chrom #is a single interval
        vals = bw.get(bedtool.chrom, bedtool.start, bedtool.stop)
        consvals = list(v[-1] for v in vals)
        if len(consvals) > 0:
            mean_phastcons = np.mean(consvals)
        else:
            mean_phastcons=0
        data = mean_phastcons
    except:
        
        #if bedtool
        for i, bedline in enumerate(bedtool):
            data = np.ndarray(len(bedtool))        
            vals = bw.get(bedline.chrom, bedline.start, bedline.stop)
            consvals = list(v[-1] for v in vals)
            if len(consvals) > 0:
                mean_phastcons = np.mean(consvals)
            else:
                mean_phastcons=0
            data[i] = mean_phastcons
            
    #returns mean phastcons score for each line 
    #returns inconistant data types, need to convert so it just returns an array 
    return data


def fa_file(filename, region = None, fd=None, type= "real"):
    
    """
    
    Way of organizing fasta file names  
    Checks if a fasta file exists returns the file attaced to a region
    or something 
    
    """
    
    if not os.path.exists(fd):
        raise Exception
    
    if region is not None:
        x =filename+"."+  region+ "."+ type+ ".fa"
        return os.path.join(fd, x)
    else:
        x = filename+ "."+ type + ".fa"
        return os.path.join(fd, x)

def make_dir(dir_name):
    
    """ makes a dir, dir_name if it doesn't exist otherwise does nothing """
    
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)
            
def main(options):
    
    """
    
    Runs all analysies 
    
    one thing to do is make graphs fail gracefully 
    
    """
    
    #from subprocess import Popen, PIPE
    #host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
    
    #gets clusters in a bed tools + names species 
    clusters = options.clusters
    species = options.species
    clusters_bed = pybedtools.BedTool(clusters)

    #makes output file names 
    clusters = str.replace(clusters, ".BED", "")
    options.k = map(int, options.k)
    outdir = options.outdir
    
    #sets up output dirs
    make_dir(outdir)        

    assigned_dir = os.path.join(outdir, "assigned")
    misc_dir = os.path.join(outdir, "misc")
    fastadir = os.path.join(outdir, "fasta")    
    kmerout = os.path.join(outdir, "kmer")
    homerout_base = os.path.join(outdir, "homer")
    make_dir(homerout_base)
    homerout = os.path.join(homerout_base, clusters)    

    make_dir(assigned_dir)
    make_dir(misc_dir)
    make_dir(fastadir)
    make_dir(homerout)
    make_dir(kmerout)

    all_regions = (["all", "exon", "UTR3", "UTR5", "proxintron", "distintron"])    

    #Not quite sure whats going on here, but its one logical block
    #either reassigns clusters to genic regions or reads from already
    #made assigned lists

    if options.assign is False:
        try:
            cluster_regions, sizes, Gsizes = build_assigned_from_existing(assigned_dir, clusters, all_regions, options.nrand)
            print "I used a pre-assigned set output_file BED files... score!"
        except:
            print "I had problems retreiving region-assigned BED files from %s, i'll rebuild" % (assigned_dir)
            options.assign = True
            
    if options.assign is True:
        print "Assigning Clusters to Genic Regions"
        cluster_regions, sizes, Gsizes = assign_to_regions(clusters_bed,options.genome_location, options.regions_location, species=species, getseq=True, nrand=options.nrand)
        print "Done Assigning"
        
        print "Saving BED and Fasta Files...",

        #outputs little files (maybe move inside output_file assign to regions)
        sizes_out = open(os.path.join(assigned_dir, "%s.sizes.pickle" %(clusters)), 'w')
        pickle.dump(sizes, file=sizes_out)
        sizes_out.close()    
        Gsizes_out = open(os.path.join(assigned_dir, "Gsizes.pickle"), 'w')
        pickle.dump(Gsizes, file=Gsizes_out)
        Gsizes_out.close()
        
        #this is where all saving happens for assign to regions
        for region in all_regions:
            output_file = clusters + "." + region+ ".real.BED"
            try:
                cluster_regions[region]['real'].saveas(os.path.join(assigned_dir, output_file))
            except:
                continue
            for n in range(options.nrand):
                output_file = clusters + "." + region+ ".rand." + str(n) + ".BED"
                try:
                    cluster_regions[region]['rand'][n].saveas(os.path.join(assigned_dir, output_file))
                except:
                    continue
                
        print "done"

        #creates pretty file names for all regions
        for region in all_regions:
            try:
                real_fa = fa_file(clusters, region=region, fd = fastadir, type="real")
                rand_fa = fa_file(clusters, region=region, fd = fastadir, type="random")
                cluster_regions[region]['real'].save_seqs(real_fa)

                l = list()#list output_file randoms
                for n in cluster_regions[region]['rand'].keys():
                    l.append(cluster_regions[region]['rand'][n])
                write_seqs(rand_fa, l)        
            except:
                continue            
                                   
    print "Counting reads in clusters...",
    
    #generates data for figure 1 and 2
    #gets reads in clusters (figure 1)
    #gets reads per cluster (figure 2)
    reads_in_clusters = 0
    reads_per_cluster = list()
    for cluster in cluster_regions['all']['real']:
        chr, start, stop, name, score, strand, tstart, tstop = str(cluster).strip().split("\t")
        try:
            gene, n, reads = name.split("_")
        except:
            try:
                gene, n, reads = name.split(";")[0].split("_")
            except:
                pass
        if int(reads)> 1:
            reads_per_cluster.append(int(reads))
        reads_in_clusters += int(reads)
    print "done"
    
    #need to get rid output_file this pickleing busniess, its a waste output_file space and doesn't work with other methods
    #gets total number output_file reads (figure 1)
    #gets total number output_file reads from clipper analysis (Need to make clipper automatically output
    #pickle file
    print "Getting total number output_file reads...",
    total_reads = 0;
    try:
        pickle_file = clusters + ".pickle"
        if os.path.exists(pickle_file):
            pf = pickle.load(open(pickle_file, 'rb'))
        else:
            print "Couldn't find %s" %(pickle_file)
        print "Found %s" %(pickle_file)
        for gene in pf:
            total_reads += gene['nreads']
    
    #if clipper didn't output gets it from flagstat
    except:
        print "Couldn't find a pickled file, resorting to flagstat for total reads. (this includes intergenic reads)"
        flagstats = pysam.flagstat(options.bam)
        total_reads =int(flagstats[2].split(" ")[0])
        
    print "done, there were %d" %(total_reads)
    print "Gathering bed lengths...",
    
    #one stat is just generated here
    #generates cluster lengths (figure 3)
    cluster_lengths = bedlengths(cluster_regions['all']['real'])
    print "done"
    
    ##This should be abstracted to some sort output_file list or something...
    #figures 5 and 6, builds pre-mrna, mrna exon and intron distributions 
    mRNA_positions = list()
    premRNA_positions = list()
    intron_positions = list()
    exon_positions = list()
    
    #also builds figure 10 (exon distances)
    GENES, Gtypes = build_AS_STRUCTURE_dict(species, options.as_structure)
    types = {}
    for type in ["CE:", "SE:", "MXE:", "A5E:", "A3E:"]:
        types[type]=0
    print "locating clusters within genes",
    
    
    try:
        #counts nearest exon to peak and gets RNA 
        #gets rna positon for every line as well
        for line in (cluster_regions['all']['real']):
            mRNA_frac, premRNA_frac, exon_frac, intron_frac, nearest_type = RNA_position(line, GENES)
            if mRNA_frac is not None:
                mRNA_positions.append(mRNA_frac)
            if exon_frac is not None:
                exon_positions.append(exon_frac)
            if premRNA_frac is not None:
                premRNA_positions.append(premRNA_frac)
            if intron_frac is not None:
                intron_positions.append(intron_frac)
            if nearest_type is not None:
                try:
                    types[nearest_type] += 1
                except:
                    types[nearest_type] =1
    except:
        print "there were errors, skipping"
    print "done"
    
    #gtypes is total genomic content 
    #types is what clusters are
    #generates figure 10 (exon distances)
    type_count = [types["CE:"], types["SE:"], types["MXE:"], types["A5E:"], types["A3E:"]]
    Gtype_count = [Gtypes["CE:"], Gtypes["SE:"], Gtypes["MXE:"], Gtypes["A5E:"], Gtypes["A3E:"]]    

    ### write fasta files and run homer and/or kmer analysis if at least one analysis is requested
    #runs kmer and homer analysis
    
    kmer_results = {} 
    if options.reMotif is True:
        for region in all_regions:

            #reads nicely named files 
            real_fa = fa_file(clusters, region=region, fd =  fastadir, type="real")
            rand_fa = fa_file(clusters, region=region, fd =  fastadir, type="random")
            if options.k is not None:
                if options.homer is True:
                    region_homer_out = os.path.join(homerout, region)
                    run_homer(real_fa, rand_fa, options.k,  outloc=region_homer_out)
                for k in options.k:
                    kmer_results[k] = {}
                    kmer_results[k][region] = kmerdiff(real_fa, rand_fa, k)
                    kmerfile = clusters + ".k" + str(k) + "." + region + ".kmerdiff"
                    kmerfile = os.path.join(kmerout, kmerfile)
                    kmer_sorted_output = run_kmerdiff(real_fa, rand_fa, outfile=kmerfile, k=k)

            
    #all the different motifs that the user specifices 
    motifs = list(options.motif)
    kmer_box_params = [kmerout, clusters, options.k, motifs]

    ###conservation --should use multiprocessing to speed this part up!
    #start output_file conservation logic, very slow...
    phast_values = list()
    
    #loads phastcons values output_file generates them again
    if options.rePhast is False:
        try:
            phast_values = pickle.load(open(os.path.join(misc_dir, "%s.phast.pickle" %(clusters))))
        except:
            options.rePhast = True

    #generates again
    if options.rePhast is True:
        print "Fetching Phastcons Scores...",
        
        #phastcons values for all regions except "all"
        for region in all_regions[1:]: #skip "all" combine them later
            print ("%s..." %(region)),
            try:
                samplesize=1000
                
                #because it takes so long to fetch only select 1000 output_file them, not actually
                #implemented
                if len(cluster_regions[region]['real']) > samplesize:
                    R1 = cluster_regions[region]['real']                
                    # R1 = random.sample(cluster_regions[region]['real'], samplesize)
                else:
                    R1 = cluster_regions[region]['real']

                #realPhast = get_phastcons(cluster_regions[region]['real'], species=options.species)
                print "getting real...",
                
                #gets phastcons values real regions 
                realPhast = get_phastcons(R1, options.phastcons_location, species=options.species)
                randPhast = list()
                
                #logic for random stuff (could be precomputed)
                for i in range(options.nrand):
                    if len(cluster_regions[region]['rand'][i]) > samplesize:
                        R2 = cluster_regions[region]['rand'][i]                    
                        #R2 = random.sample(cluster_regions[region]['rand'][i], samplesize)
                    else:
                        R2 = cluster_regions[region]['rand'][i]
                    print ("getting rand %d" %(i)),
                    randPhast.extend(get_phastcons(R2, options.phastcons_location, species=options.species).tolist())
                
                #list output_file lists for real and random for every genic region
                phast_values.append(realPhast)
                phast_values.append(randPhast)
            
            except:
                continue
            
        #hacky selection of real values from phast_values
        all_real = np.concatenate(phast_values[::2])
        
        #hacky selection output_file random values from phast_values
        all_rand = np.concatenate(phast_values[1::2])
        
        #adds back in all and rand to phast_values list
        phast_values.insert(0,all_rand)
        phast_values.insert(0,all_real)
        pickout = open(os.path.join(misc_dir, "%s.phast.pickle" %(clusters)), 'w')
        pickle.dump(phast_values, file = pickout)
    
    
    Zscores = None  #old. remove
    
    #build qc figure
    QCfig_params = [reads_in_clusters, (total_reads - reads_in_clusters), cluster_lengths, reads_per_cluster, premRNA_positions, mRNA_positions, exon_positions, intron_positions, Gsizes, sizes, Gtype_count, type_count, Zscores, homerout, kmer_box_params, phast_values]

    #save results 
    pickout = open(os.path.join(outdir, "misc", "%s.qcfig_params.pickle" %(clusters)), 'w')
    pickle.dump(QCfig_params, file = pickout)
    QCfig = CLIP_Analysis_Display.CLIP_QC_figure(*QCfig_params)
    fn = clusters + ".QCfig.pdf"
    outFig = os.path.join(outdir, fn)
    
    #TODO Fix output output_file file (Don't know why its crashing right now
    print >> sys.stderr, outFig
    QCfig.savefig(outFig)
                    
    ### does something with motifs doesn't appear to work right now
    
    #reads in existing precompiled motif file
    motifs = list(options.motif)
    
    if motifs is not None and False: #TODO hack to get stuff compiling fix soon
        motifBASE  = options.motif_location
        fig = pylab.figure(figsize=(8.5, 11))
        colors = ["red", "orange", "green", "blue", "purple", "brown", "black", "pink", "gray", "cyan", "magenta"]
        for i, motif in enumerate(motifs):
            mf = "motif_" + motif + ".BED"
            mfgz = "motif_" + motif + ".BED.gz"
            print os.path.join(motifBASE,species,mf)
            motifFILE = None

            if os.path.exists(os.path.join(motifBASE,species, mf)):
                motifFILE = os.path.join(motifBASE,species, mf)
            elif os.path.exists(os.path.join(motifBASE,species, mfgz)):
                motifFILE= os.path.join(motifBASE,species, mfgz)
            else:
                print "MOTIF BED FILE for motif: %s is not available, please build it" %(mf)
                continue
            
            #plots motif distance from the precompiled file to the clusters 
            plot_motif_dist(cluster_regions, motifFILE, fig, color = colors[i], species=species, slopsize=200)
        pylab.savefig(clusters + ".motif_distribution.pdf")
        
        #fin
if __name__== "__main__":
    parser = OptionParser()
    
    parser.add_option("--clusters", dest="clusters", help="BED file of clusters", metavar="BED")
    parser.add_option("--bam", dest="bam")
    parser.add_option("--species", "-s", dest="species", help = "genome version")
    ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
    parser.add_option("--reAssign", dest="assign", action="store_true", default=False, help="re-assign clusters, if not set it will re-use existing assigned clusters") 
    ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
    parser.add_option("--rePhast", dest="rePhast", action="store_true", default=False, help="re-calculate conservation, must have been done before") 
    parser.add_option("--old_motifs", dest="reMotif", action="store_false", default=True, help="use old motif files")
    parser.add_option("--motif", dest="motif", action="append", help="Files of motif locations", default=None)
    parser.add_option("--homer", dest="homer", action="store_true", help="Runs homer", default=False)
    parser.add_option("--k", dest="k", action="append", help="k-mer and homer motif ananlysis", default=[6])
    parser.add_option("--conservation", dest="cons", help="Runs conservation (might not do anything)", action="store_true")
    parser.add_option("--structure", dest="structure", help="also doesn't do anything gets structure maps", action="store_true")
    parser.add_option("--nrand", dest="nrand", default=3, help="selects number of times to randomly sample genome", type="int")
    parser.add_option("--outdir", "-o", dest="outdir", default=os.getcwd(), help="directory for output, default:cwd")
    parser.add_option("--run_phast", dest="run_phast", action="store_true", help="re-runs phastcons (not be implemented)", default=False)
    ##Below here are critical files that always need to be referenced
    parser.add_option("--AS_Structure", dest="as_structure",  help="Location of AS_Structure directory (chromosme files should be inside)", default=None)
    parser.add_option("--genome_location", dest="genome_location", help="location of all.fa file for genome of interest", default=None)
    parser.add_option("--homer_path", dest="homer_path", action="append", help="path to homer, if not in default path", default=None)
    parser.add_option("--phastcons_location", dest="phastcons_location",  help="location of phastcons file", default=None)
    parser.add_option("--regions_location", dest="regions_location",  help="directory of genomic regions for a species", default=None)
    parser.add_option("--motif_directory", dest="motif_dir",  help="directory of pre-computed motifs for analysis", default=None)
   
    (options, args) = parser.parse_args()
    
    #error checking
    if options.clusters is None or options.bam is None or options.species is None:
        parser.print_help()
        exit()
        
    main(options)

