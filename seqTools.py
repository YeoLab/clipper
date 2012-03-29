
import pybedtools

from optparse import OptionParser
import os
import sys
import random
import pickle
from numpy import *
import subprocess
from bx.bbi.bigwig_file import BigWigFile
import pysam
from subprocess import Popen, PIPE
host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
if "optiputer" in host or "compute" in host:
    basedir = "/nas/nas0"
elif "tcc" in host or "triton" in host:
    basedir = "/projects"
else:
    print "Where am I?"
    raise Exception


if not os.path.exists("/nas3/scratch/lovci/pybedtools_tmp"):
    pybedtools.set_tempdir("/projects/lovci/projects/tmp/python_tmp")
else:
    pybedtools.set_tempdir("/nas3/scratch/lovci/pybedtools_tmp")
    pass

def build_AS_STRUCTURE_dict(species):
    if species == "hg19":
        chrs = map(str,range(1,23)) #1-22
        chrs.append("X")
        chrs.append("Y")
    elif species == "mm9":
        chrs = map(str,range(1,20)) #1-19
        chrs.append("X")
        chrs.append("Y")        

    info = dict()
    Gtypes = dict()
    for chr in chrs:
        ASfile = basedir + "/yeolab/Genome/ensembl/AS_STRUCTURE/" + species + "data4/" + species + ".tx." + chr + ".AS.STRUCTURE"
        f = open(ASfile, "r")
        for line in f.readlines():
            if not line.startswith(">"):
                continue
            
            blank, gene, chr, transcripts, d2, d3, d4, strand, numex, exonloc, intronloc, exonlen, intronlen, asType, locType = line.strip().split("\t")
            signstrand = "-"
            if int(strand) == 1:
                signstrand = "+"
            numex = int(numex)
            info[gene]={}
            info[gene]['strand'] = signstrand
            info[gene]['exons'] = {}
            info[gene]['introns'] = {}
            info[gene]['types'] = {}                        
            exons = exonloc.split("|")
            introns = intronloc.split("|")
            types = asType.split("|")
            info[gene]['numEx'] = numex
            info[gene]['mRNA_length'] =0
            info[gene]['premRNA_length'] =0
            tx_start = 1000000000000
            tx_stop = -1000000000000
            for i, exon in enumerate(exons):
                if i == numex: #last exon is empty
                    continue
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

def assign_to_regions(tool, species="hg19", nrand = 3, getseq=False):
    speciesFA = ""
    if species =="hg19" or species == "mm9" or species =="hg18":
        speciesFA = basedir + "/yeolab/Genome/ucsc/" +species + "/chromosomes/all.fa"
    else:
        raise Exception, "Unknown species"
    root = basedir + "/lovci/projects/ucscBED/" + species
    UTR3File=os.path.join(root, "UTR3_"+species+"_frea_sorted.withscore")
    UTR3 = pybedtools.BedTool(UTR3File)
    G_UTR3_size = UTR3.total_coverage()
    UTR5File = os.path.join(root,"UTR5_"+species+"_frea_sorted.withscore")
    UTR5 = pybedtools.BedTool(UTR5File)
    G_UTR5_size = UTR5.total_coverage()
    exonFile=os.path.join(root,"exon_"+species+"_frea_sorted.withscore")
    exon = pybedtools.BedTool(exonFile)
    G_exon_size = exon.total_coverage()    
    proxintronFile= os.path.join(root,"proxintron500_"+species+"_frea_sorted.withscore")
    proxintron = pybedtools.BedTool(proxintronFile)
    G_proxintron_size = proxintron.total_coverage()        
    distintronFile = os.path.join(root,"distintron500_"+species+"_frea_sorted.withscore")
    distintron = pybedtools.BedTool(distintronFile)
    G_distintron_size = distintron.total_coverage()    
    all_regions = (["exon", "UTR3", "UTR5", "proxintron", "distintron"])
    all_bedtracks = ([exon, UTR3, UTR5, proxintron, distintron])
    all_regionfiles = ([exonFile, UTR3File, UTR5File, proxintronFile, distintronFile])
    bed_dict = {}
    last = 0
    of = get_offsets_bed12(tool)
    tool =  tool.merge(s=True, nms=True, scores="max")
    tool = adjust_offsets(tool, of)
    Gsizes = [G_exon_size, G_UTR3_size, G_UTR5_size, G_proxintron_size, G_distintron_size]
    sizes = [0, 0, 0, 0, 0]
    print "There are a total of %d clusters I'll examine" %(tool.__len__())
    for j, region in enumerate(all_regions):
        no_overlapping, only_overlapping =intersection(tool, B = all_bedtracks[j])  #portions of the regions that overlap with a genic region
        only_overlapping = pybedtools.BedTool(str(only_overlapping.filter(eliminate_invalid_bed12)), from_string=True)
        no_overlapping = pybedtools.BedTool(str(no_overlapping.filter(eliminate_invalid_bed12)), from_string=True)
        noCT = no_overlapping.__len__()
        onlyCT = only_overlapping.__len__()
        last += onlyCT
        both = last + noCT
        print "For region: %s I found %d that overlap and %d that don't, making a total of %d" %(region, onlyCT, noCT, both)
        bed_dict[region]= {}
        bed_dict[region]['real'] = only_overlapping.sort()
        try:
            bed_dict['all']['real'] = pybedtools.BedTool(str(bed_dict['all']['real']) + str(bed_dict[region]['real']), from_string=True)
        except:
            bed_dict['all']= {}
            bed_dict['all']['rand']= {}
            bed_dict['all']['real'] = bed_dict[region]['real']
        sizes[j] = bed_dict[region]['real'].total_coverage()
        if getseq is True:
            bed_dict[region]['real'].sequence(fi=speciesFA, s=True)
        ofdic = get_offsets_bed12(only_overlapping)
        bed_dict[region]['rand'] = {}
        for i in range(nrand):
            ri = bed_dict[region]['real'].shuffle(genome=species, incl=all_regionfiles[j]).sort()
            ri = adjust_offsets(ri, ofdic)
            if getseq is True:
                ri.sequence(fi=speciesFA, s=True)
                
            bed_dict[region]['rand'][i] = ri
            try:
                bed_dict['all']['rand'][i] = pybedtools.BedTool(str(bed_dict['all']['rand'][i]) + str(bed_dict[region]['rand'][i]), from_string=True)
            except:
                bed_dict['all']['rand'][i] = bed_dict[region]['rand'][i]

        tool = no_overlapping
    #colors = ["#E52C27", "#C3996B", "#3C54A4", "#48843D", "#852882"] #red, tan, blue, green, purple
    #labels = ["Exon", "3'UTR", "5'UTR", "Proximal Intron", "Distal Intron"]
    #pie_data = piefig.add_subplot(121, aspect=1, title="Clusters' Coverage")
    #pie_genomic = piefig.add_subplot(122, aspect=1, title="Genomic Coverage")
    #pie_data.pie(sizes, colors=colors, labels=labels)
    #pie_genomic.pie(Gsizes, colors=colors, labels=labels)
    print "After assigning, I\'m left with %d un-categorized regions" %(tool.__len__())
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
        sizes=0
    try:
        Gsizes = pickle.load(open(os.path.join(assigned_dir, "Gsizes.pickle", 'rb')))
    except:
        Gsizes=0


    
    return CLUS_regions, sizes, Gsizes


def readsToWiggle_pysam(reads, tx_start, tx_end, keepstrand=None, trim=False, usePos='center'):
    wiggle = zeros((tx_end - tx_start + 1), dtype='f') #total read overlapping each position
    pos_counts = zeros((tx_end - tx_start + 1), dtype='f') #each read only gives one point, use usePos to determine whether the point comes from the start, center, or end
    allreads = {}
    jxns = {}
    lengths = list()
    for read in reads:
        if read.is_reverse is True and keepstrand is "+":
            continue
        if read.is_reverse is False and keepstrand is "-":
            continue
        aligned_positions = read.positions
        if aligned_positions[0] < tx_start or aligned_positions[-1] > tx_end:
            continue # skip reads that fall outside the gene bounds
        readpos = str(aligned_positions[0]) +"-" + str(aligned_positions[-1])
        read_start = aligned_positions[0]
        read_stop = aligned_positions[-1]        
        if trim is True and readpos in reads:
            continue
        lengths.append(read.qlen)        
        try:
            if usePos == "center":
                center = (int(aligned_positions[-1]+aligned_positions[0])/2)-tx_start
                pos_counts[center] += 1
            elif usePos == "start":
                if keepstrand == "+":
                    pos_counts[(aligned_positions[0]-tx_start)]+=1
                else:
                    pos_counts[(aligned_positions[-1]-tx_start)]+=1
            elif usePos == "end":
                if keepstrand == "+":
                    pos_counts[(aligned_positions[-1]-tx_start)]+=1
                else:
                    pos_counts[(aligned_positions[0]-tx_start)]+=1
            else:
                raise NameError('PositionType')
        except NameError:
            print "Invalid value for usePos, choose \"center\", \"start\", or \"end\""
        try:
            allreads[read_start][read_stop]+=1
        except:
            try:
                allreads[read_start][read_stop]=1
            except:
                allreads[read_start] = {}
                allreads[read_start][read_stop]=1
                
        for i,pos in enumerate(aligned_positions):

            if pos < tx_start or pos > tx_end:
                continue
            wig_index = pos-tx_start
            #wiggle[wig_index] += 1./read.qlen
            wiggle[wig_index] += 1.
            try:
                #if there is a junction coming up
                if i+1 < len(aligned_positions) and aligned_positions[i+1] > pos + 1:
                    leftss = pos+1
                    rightss= aligned_positions[i+1]+1
                    if leftss > tx_start and leftss < tx_end and rightss > tx_start and rightss < tx_end:
                        jxn = ":".join(map(str, [leftss, rightss]))
                        try:
                            jxns[jxn] += 1 
                        except:
                            jxns[jxn] = 1
            except:
                print "junction parsing broke"
                pass
    return wiggle, jxns, pos_counts, lengths, allreads
