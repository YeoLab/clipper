import matplotlib as mpl
mpl.use('Agg')    
import pybedtools
import numpy as np
from optparse import OptionParser
import os
import sys
import pickle
import random
from subprocess import Popen, call, PIPE
import pylab

host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
if "optiputer" in host or "compute" in host:
    basedir = "/nas/nas0"
elif "tcc" in host or "triton" in host:
    basedir = "/projects"    
print "basedir: %s" %(basedir)


import subprocess
from bx.bbi.bigwig_file import BigWigFile
import pysam
try:
    pybedtools.set_tempdir(basedir + "/scratch/lovci/pybedtools_tmp")
except:
    pybedtools.set_tempdir(basedir + "/lovci/projects/tmp/python_tmp")

def CLIP_QC_figure(reads_in_clusters, reads_out_clusters, cluster_lengths, reads_per_cluster, premRNA, mRNA, exondist, introndist, genomic_locs, clusters_locs, genomic_types, clusters_types, zscores, homer_location, kmer_box_params, phastcons_values):
    import matplotlib as mpl
    import pylab
    from pylab import cm
    import matplotlib.gridspec as gridspec
    import matplotlib.image as mpimg
    import numpy as np
    import random
    mpl.rcParams['interactive']=False

    fig = pylab.figure(figsize=(20,20), facecolor='white')
    gs1 = gridspec.GridSpec(6,4)
    ax_pie_inout = pylab.subplot(gs1[0, 0], title = "Reads in Clusters", aspect=1)
    ax_pie_inout.pie([reads_in_clusters, reads_out_clusters], labels=["In Clusters", "Outside of\nClusters"])
    gs_nreads_lengths_cons = gridspec.GridSpecFromSubplotSpec(1,7,subplot_spec=gs1[0,1:4])
    #ax_nreads = pylab.subplot(gs1[0, 1:3], title="Reads per Cluster")
    ax_nreads = pylab.subplot(gs_nreads_lengths_cons[0:3], title="Reads per Cluster")
    ax_nreads.hist(reads_per_cluster, bins=50, facecolor='#C8D2B0',log=True, range=(10, np.max(reads_per_cluster)))
    ax_nreads.set_xscale('log')
#    ax_nreads.set_yscale('log')        
    ax_nreads.set_xlabel("log10(N reads)")
    ax_nreads.set_ylabel("Frequency")
    #ax_lengths = pylab.subplot(gs1[0, 3], title="Cluster Lengths")
    ax_lengths = pylab.subplot(gs_nreads_lengths_cons[3], title="Cluster Lengths")
    ax_lengths.set_yscale('log')
    ax_lengths.boxplot(random.sample(cluster_lengths, 2000), vert=1)
    #ax_lengths.set_ylabel("log10(length)")
    ax_lengths.set_xticklabels([])
    formal_labels = ["All", "Exon", "3'UTR", "5'UTR", "Proximal\nIntron", "Distal\nIntron"]
    positions=list()
    width=.1
    for i in range(6):
        positions.append(i-.1)
        positions.append(i)
    ax_cons = pylab.subplot(gs_nreads_lengths_cons[4:7], title="PhastCons Values")
    bp = ax_cons.boxplot(phastcons_values, positions=positions, widths=width, sym='k.')
    boxes = bp['boxes']

    for realbox in boxes[::2]:
        realbox.set_color('blue')

    for randbox in boxes[1::2]:
        randbox.set_color('red')

    outliers = bp['fliers']
    for i in outliers:
        i.set_markersize(3.)

    medians = bp['medians']
    for line in medians[::2]:
        line.set_color('blue')
    for line in medians[1::2]:
        line.set_color('red')        
    ax_cons.hold(True)
    start, stop = ax_cons.get_xlim()
    ticks = np.linspace(start, stop, 7)
    starts= ticks[:-1]
    stops = ticks[1:]
    reals = phastcons_values[::2]
    rands = phastcons_values[1::2]

    bins=20
    for n in range(6):
        ax_cons.hold(True)        
        heights1, edges1 = np.histogram(reals[n], normed=True, bins=bins)
        heights2, edges2 = np.histogram(rands[n], normed=True, bins=bins)
        CDF1=np.cumsum(heights1)/np.sum(heights1)
        CDF1 = np.insert(CDF1, 0, 0)
        CDF2=np.cumsum(heights2)/np.sum(heights2)
        CDF2 = np.insert(CDF2, 0, 0)        
        #xvals = np.linspace(starts[n], stops[n], (bins+1))
        yvals = np.linspace(0,1, (bins+1))
        
        if n ==0:
            label1 = "real"
            label2 = "random"
        else:
            label1 = "_nolegend_"
            label2 = "_nolegend_"            
        ax_cons.plot( starts[n]+CDF1, yvals, 'blue', label=label1)
        ax_cons.plot( starts[n]+CDF2, yvals,  'red', label=label2)
    for x in starts[1:]:
        ax_cons.axvline(x, color='k', lw=2)
    

    
    ax_cons.set_xticks([0,1,2,3,4,5])
    ax_cons.set_xticklabels(formal_labels)
    ax_cons.set_ylabel("PhastCons Score")

    gs_line2 = gridspec.GridSpecFromSubplotSpec(1,6, subplot_spec=gs1[1,:])
    #ax_genedist = pylab.subplot(gs1[1, 0:2], title="mRNA and preMRNA Distribution")
    ax_genedist = pylab.subplot(gs_line2[0:2], title="mRNA and pre-mRNA Distribution")
    ax_genedist.set_xlabel("Fraction of region")
    ax_genedist.set_ylabel("Frequency")
    prehist = ax_genedist.hist(premRNA, range=(0,1.0), histtype="step", color="red", label="pre-mRNA", bins=100, normed=True)
    ax_genedist.hold(True)

    mhist = ax_genedist.hist(mRNA, range=(0,1.0), histtype="step", color="blue", label="mRNA", bins=100, normed=True)
    ax_genedist.legend(frameon=False, loc=2, markerscale=0.2, )

    #ax_exondist = pylab.subplot(gs1[1, 2:4], title="Distribution across exons")
    ax_exondist = pylab.subplot(gs_line2[2:4], title="Distribution Across Exons and Introns")
    ax_exondist.hist(exondist, range=(0,1.0), histtype="step", color="blue", label="Intron", bins=100, normed=True)
    ax_exondist.set_xlabel("Fraction of region")
    ax_exondist.set_ylabel("Exon Location Frequency", color='blue')
    #ax_introndist = pylab.subplot(gs1[2, 2:4], title="Distribution across introns")
    #ax_introndist.set_xlabel("Fractionof region")    

    ax_introndist = ax_exondist.twinx()
    ax_introndist.hist(introndist, range=(0,1.0), histtype="step", color="red", label="Intron", bins=100, normed=True)

    ax_introndist.set_ylabel("Intron Location Frequency", color='red')
    for tick in ax_introndist.get_yticklabels():
        tick.set_color('red')
    for tick in ax_exondist.get_yticklabels():
        tick.set_color('blue')
        
        #ax_introndist.set_xticklabels([])
    pylab.tight_layout()

    colors = ["#E52C27", "#C3996B", "#3C54A4", "#48843D", "#852882"] #red, tan, blue, green, purple
    labels = ["Exon", "3'UTR", "5'UTR", "Proximal\nIntron", "Distal\nIntron"]
    gs_pie_nearestType = gridspec.GridSpecFromSubplotSpec(1,3, subplot_spec=gs1[2,0:2])
    ax_pie_genomic = pylab.subplot(gs_pie_nearestType[0], title="Genomic Content", aspect=1)
    ax_pie_clusters = pylab.subplot(gs_pie_nearestType[1], title="Clusters' Content", aspect=1)
    ax_bar_exontypes = pylab.subplot(gs_pie_nearestType[2], title="Nearest Exon Types")    

        #ax_pie_clusters = pylab.subplot(gs1[2, 1], title="Clusters' Content", aspect=1)
        #ax_pie_genomic = pylab.subplot(gs_pies[2, 0], title="Genomic Content", aspect=1)
        #ax_bar_exontypes = pylab.subplot(gs1[3,0], title="Nearest Exon Types")        
    ax_pie_genomic.pie(genomic_locs, colors=colors, labels=labels)
    ax_pie_clusters.pie(clusters_locs, colors=colors, labels=labels)

    ind = np.arange(5)
    clusters_types = 100*np.array(clusters_types, dtype="float")/np.sum(clusters_types)
    genomic_types = 100*np.array(genomic_types, dtype="float")/np.sum(genomic_types)
    difference = clusters_types-genomic_types
    #width=0.35
    #bar1=ax_bar_exontypes.bar(ind, clusters_types, color='r', width=width)
    #bar2=ax_bar_exontypes.bar(ind+width, genomic_types, color='y', width=width)
    ind = ind-.5
    bar2=ax_bar_exontypes.bar(ind, difference, color='y')
    xlabels = ["", "CE", "SE", "MXE", "A5E", "A3E"]
    ax_bar_exontypes.set_xticklabels(xlabels)
    ax_bar_exontypes.set_ylabel('% Difference (clusters-genomic)')
    ax_bar_exontypes.set_xlabel('Exon Type')
    ax_bar_exontypes.axhline(y=0, ls="-",color='k')
    #ax_bar_exontypes.legend( (bar1[0], bar2[0]), ('Clusters', 'Genome'), frameon=False)
    #ax_pie_exontypes_expected = pylab.subplot(gs1[3,0], title="Genome's Exon Types", aspect=1)
    #ax_pie_exontypes_expected.pie(genomic_types, labels=["CE", "SE", "MXE", "A5E", "A3E"], colors=["darkred", "navy", "purple", "green", "orange"])
    #ax_pie_exontypes = pylab.subplot(gs1[3,1], title="Nearest Exon Type", aspect=1)
    #ax_pie_exontypes.pie(clusters_types, labels=["CE", "SE", "MXE", "A5E", "A3E"], colors=["darkred", "navy", "purple", "green", "orange"])
    ax_hist_zscores = pylab.subplot(gs1[2,2:4], title="Motif Z-scores")
    motif_boxplots(*kmer_box_params, subplot=ax_hist_zscores) # * expands the list into individual args

    #ax_hist_zscores.hist(zscores, bins=100, fc=None)
    #ax_hist_zscores.set_ylabel("Frequency")
    #ax_hist_zscores.set_xlabel("Z-score")
    pylab.tight_layout()
    gs2 = gridspec.GridSpecFromSubplotSpec(1,6, subplot_spec=gs1[4:6,:], hspace=0, wspace=0)
    all_regions = (["all", "exon", "UTR3", "UTR5", "proxintron", "distintron"])
    for i, region in enumerate(all_regions):
        gs_homer_motifs = gridspec.GridSpecFromSubplotSpec(8,1, subplot_spec=(gs2[i]))
        for j, space in enumerate(gs_homer_motifs):
            try:
                motifname = "motif" + str(j+1) + ".logo.png"
                formal_labels = ["All", "Exon", "3'UTR", "5'UTR", "Proximal Intron", "Distal Intron"]
                motifFILE = os.path.join(homer_location, region, "homerResults", motifname)
                if os.path.exists(motifFILE):
                    motif = mpimg.imread(motifFILE)
                    if j == 0:
                        title = formal_labels[i]
                        ax = pylab.subplot(space, frameon=False, xticks=[], yticks=[], title=title)                        
                    else:
                        ax = pylab.subplot(space, frameon=False, xticks=[], yticks=[]) #no title
                    ax.imshow(motif)
                else:
                    print "no motif %s" %(motifFILE)
            except:
                pass
    return fig

def intersection(A, B=None):
#    tmp = outfile+ ".tmp"
#    of = open(outfile, 'w')
#    ot = open(tmp, 'w')
    less = A.subtract(B, s=True) #without regions that overlap
    more = A.subtract(less, s=True) #only regions that overlap
    return less, more


def adjust_offsets(tool, offsets=None):
    if offsets == None:
        raise Exception, "no offsets"
    l = list()
    for x in tool:
        try:
            chr, start, stop, name, score, strand, tstart, tstop = x.__str__().strip().split("\t")
        except:
            chr, start, stop, name, score, strand = x.__str__().strip().split("\t")
            tstart=0
            tstop=0
        start, stop, tstart, tstop = map(int, (start, stop, tstart, tstop))
        if ";" in name and name not in offsets:
            o = name.split(";")
            offset = offsets[o[0]]
        else:
            offset = offsets[name]
        if strand== "+":
            tstart = start + offset
            tstop = tstart + 4
        else:
            tstop = stop - offset
            tstart = tstop -4
        l.append("\t".join(map(str, (chr, start, stop, name, score, strand, tstart, tstop))))
    t = pybedtools.BedTool("\n".join(l), from_string=True)
    return t




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
        ASfile = basedir+ "/yeolab/Genome/ensembl/AS_STRUCTURE/" + species + "data4/" + species + ".tx." + chr + ".AS.STRUCTURE"
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
        try:
            no_overlapping, only_overlapping =intersection(tool, B = all_bedtracks[j])  #portions of the regions that overlap with a genic region
        except:
            continue
        only_overlapping = pybedtools.BedTool(str(only_overlapping.filter(eliminate_invalid_bed12)), from_string=True)
        no_overlapping = pybedtools.BedTool(str(no_overlapping.filter(eliminate_invalid_bed12)), from_string=True)
        noCT = no_overlapping.__len__()
        onlyCT = only_overlapping.__len__()
        last += onlyCT
        both = last + noCT
        print "For region: %s I found %d that overlap and %d that don't, making a total of %d" %(region, onlyCT, noCT, both)
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
        for i in range(nrand):
            try:
                ri = bed_dict[region]['real'].shuffle(genome=species, incl=all_regionfiles[j]).sort()
            except:
                continue
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
        sizes=[1,1,1,1,1]
    try:
        Gsizes = pickle.load(open(os.path.join(assigned_dir, "Gsizes.pickle"), 'rb'))
    except:
        Gsizes=[1,1,1,1,1]


    
    return CLUS_regions, sizes, Gsizes

def readsToWiggle_pysam(reads, tx_start, tx_end, strand=None, trim=False, usePos='center'):
    wiggle = zeros((tx_end - tx_start + 1), dtype='f') #total read overlapping each position
    pos_counts = zeros((tx_end - tx_start + 1), dtype='f') #each read only gives one point, use usePos to determine whether the point comes from the start, center, or end
    jxns = {}
    seen = {}
    lengths = list()
    for read in reads:
        #check strand, skip reads on the opposite strand
        if read.is_reverse is True and strand is "+":
            continue
        elif read.is_reverse is False and strand is "-":
            continue


        aligned_positions = read.positions
        if aligned_positions[0] < tx_start or aligned_positions[-1] > tx_end:
            continue # skip reads that fall outside the gene bounds
        lengths.append(read.qlen)
        readpos = str(aligned_positions[0]) +"|" + str(aligned_positions[-1])
        if trim is True and readpos in seen:
            continue
        try:
            if usePos == "center":
                center = (int(aligned_positions[-1]+aligned_positions[0])/2)-tx_start
                pos_counts[center] += 1
            elif usePos == "start":
                if strand == "+":
                    pos_counts[(aligned_positions[0]-tx_start)]+=1
                else:
                    pos_counts[(aligned_positions[-1]-tx_start)]+=1
            elif usePos == "end":
                if strand == "+":
                    pos_counts[(aligned_positions[-1]-tx_start)]+=1
                else:
                    pos_counts[(aligned_positions[0]-tx_start)]+=1
            else:
                raise NameError('PositionType')
        except NameError:
            print "Invalid value for usePos, choose \"center\", \"start\", or \"end\""

        seen[readpos]=1
        for i,pos in enumerate(aligned_positions):
            if pos < tx_start or pos > tx_end:
                continue
            wig_index = pos-tx_start
            #wiggle[wig_index] += 1./read.qlen
            wiggle[wig_index] += 1.
            try:
                #if there is a junction coming up                
                if aligned_positions[i+1] > pos + 1: 
                    leftss = pos+1
                    rightss= aligned_positions[i+1]+1
                    if leftss > tx_start and leftss < tx_end and rightss > tx_start and rightss < tx_end:                      
                        jxn = ":".join(map(str, [leftss, rightss]))
                        try:
                            jxns[jxn] += 1 
                        except:
                            jxns[jxn] = 1
            except:
                pass
    return wiggle, jxns, pos_counts, lengths




def eliminate_invalid_bed12(x):
    chr, start, stop, name, score, strand, tstart, tstop = x.__str__().split("\t")
    start, stop, tstart, tstop = map(int, (start, stop, tstart, tstop))
    if start < tstart and stop > tstop:
        return True
    else:
        return False

def get_offsets_bed12(tool):
    of = {}
    for line in tool:
        chr, start, stop, name, score, strand, tstart, tstop = line.__str__().split("\t")
        start, stop, tstart, tstop = map(int, (start, stop, tstart, tstop))
        if strand == "+":
            offset = tstart -start 
        else:
            offset = stop - tstop
        of[name]=offset
    return of

def get_offsets(clusters, motif, slop=500):
    ov = clusters.window(motif, w=slop, sm=True)
    distances = list()
    for line in ov:
        positions=line.__str__().split("\t")
        cluster_center = int(positions[7])-int(positions[6])/2
        motif_center = int(positions[15]) - int(positions[14])/2
        distance = motif_center - cluster_center
        if positions[5] == "-":
            distance = distance * -1
        distances.append(distance)
    del ov
    return distances
    
def plot_motif_dist(assigned_clusters, motifFILE, figure, nrand=3, color = "red", label=None, species="mm9", slopsize=0, scale='linear'):

    motifBed = pybedtools.BedTool(motifFILE)
    if label is None:
        label=motifFILE
    UTR5dist = get_offsets(assigned_clusters['UTR5']['real'], motifBed, slop=slopsize)
    rand_5UTR = list()

    for i in range(nrand):
        rand_5UTR.extend(get_offsets(assigned_clusters['UTR5']['rand'][i], motifBed, slop=slopsize))

    UTR5size = assigned_clusters['UTR5']['real'].total_coverage()
    UTR5_rand_size = UTR5size*nrand
    print "UTR5 done"
    UTR3dist = get_offsets(assigned_clusters['UTR3']['real'], motifBed, slop=slopsize)
    rand_3UTR = list()
    for i in range(nrand):
        rand_3UTR.extend(get_offsets(assigned_clusters['UTR3']['rand'][i], motifBed, slop=slopsize))

    UTR3size = assigned_clusters['UTR3']['real'].total_coverage()
    UTR3_rand_size = UTR3size*nrand            
    print "UTR3 done"
    exondist = get_offsets(assigned_clusters['exon']['real'], motifBed, slop=slopsize)
    rand_exon = list()
    for i in range(nrand):
        rand_exon.extend(get_offsets(assigned_clusters['exon']['rand'][i], motifBed, slop=slopsize))

    exonsize = assigned_clusters['exon']['real'].total_coverage()
    exon_rand_size = exonsize*nrand            
    print "exon done"

    distintrondist = get_offsets(assigned_clusters['distintron']['real'], motifBed, slop=slopsize)
    rand_distintron = list()
    for i in range(nrand):
        rand_distintron.extend(get_offsets(assigned_clusters['distintron']['rand'][i], motifBed, slop=slopsize))

    distintronsize = assigned_clusters['distintron']['real'].total_coverage()
    distintron_rand_size = distintronsize*nrand                        
    print "distintron done"
    
    proxintrondist = get_offsets(assigned_clusters['proxintron']['real'], motifBed, slop=slopsize)
    rand_proxintron = list()
    for i in range(nrand):
        rand_proxintron.extend(get_offsets(assigned_clusters['proxintron']['rand'][i], motifBed, slop=slopsize))

    proxintronsize = assigned_clusters['proxintron']['real'].total_coverage()
    proxintron_rand_size = proxintronsize*nrand

    print "proxintron done"

    allsize = UTR5size + UTR3size + exonsize + proxintronsize + distintronsize
    all_rand_size = allsize*nrand

    all = list()
    all.extend(UTR5dist)
    all.extend(UTR3dist)
    all.extend(exondist)
    all.extend(proxintrondist)
    all.extend(distintrondist)                    

    all_rand = list()
    all_rand.extend(rand_5UTR)
    all_rand.extend(rand_3UTR)
    all_rand.extend(rand_exon)
    all_rand.extend(rand_distintron)
    all_rand.extend(rand_proxintron)    

    ax_all = figure.add_subplot(321, title="All Clusters")
    ax_all.set_yscale(scale)
    ax_UTR5 = figure.add_subplot(322, title="5'UTR")
    ax_UTR5.set_yscale(scale)
    ax_exon = figure.add_subplot(323, title="Exon")
    ax_exon.set_yscale(scale)
    ax_UTR3 = figure.add_subplot(324, title="3'UTR")
    ax_UTR3.set_yscale(scale)
    ax_proxintron = figure.add_subplot(325, title="Proximal Intron")
    ax_proxintron.set_yscale(scale)
    ax_distintron = figure.add_subplot(326, title="Distal Intron")
    ax_distintron.set_yscale(scale)

    all_hist, all_edges = np.histogram(all, bins=50, range=(-150, 150))
    all_hist = all_hist/(allsize/1000.)        
    all_rand_hist, all_edges_rand = np.histogram(all_rand, bins=50, range=(-150, 150))
    all_rand_hist = all_rand_hist/(all_rand_size/1000.)

    ax_all.plot(all_edges[:-1], all_hist, c=color, linestyle='solid', label=label)
    ax_all.plot(all_edges_rand[:-1], all_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    UTR5_hist, UTR5_edges = np.histogram(UTR5dist, bins=50, range=(-150, 150))
    UTR5_hist = UTR5_hist/(UTR5size/1000.)        
    UTR5_rand_hist, UTR5_edges_rand = np.histogram(rand_5UTR, bins=50, range=(-150, 150))
    UTR5_rand_hist = UTR5_rand_hist/(UTR5_rand_size/1000.)       
    ax_UTR5.plot(UTR5_edges[:-1], UTR5_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_UTR5.hold(True)
    ax_UTR5.plot(UTR5_edges_rand[:-1], UTR5_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    UTR3_hist, UTR3_edges = np.histogram(UTR3dist, bins=50, range=(-150, 150))
    UTR3_hist = UTR3_hist/(UTR3size/1000.)        
    UTR3_rand_hist, UTR3_edges_rand = np.histogram(rand_3UTR, bins=50, range=(-150, 150))
    UTR3_rand_hist = UTR3_rand_hist/(UTR3_rand_size/1000.)       
    ax_UTR3.plot(UTR3_edges[:-1], UTR3_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_UTR3.hold(True)
    ax_UTR3.plot(UTR3_edges_rand[:-1], UTR3_rand_hist, linestyle='dashed', c=color, label="_nolegend_")    

    exon_hist, exon_edges = np.histogram(exondist, bins=50, range=(-150, 150))
    exon_hist = exon_hist/(exonsize/1000.)        
    exon_rand_hist, exon_edges_rand = np.histogram(rand_exon, bins=50, range=(-150, 150))
    exon_rand_hist = exon_rand_hist/(exon_rand_size/1000.)       
    ax_exon.plot(exon_edges[:-1], exon_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_exon.hold(True)
    ax_exon.plot(exon_edges_rand[:-1], exon_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    distintron_hist, distintron_edges = np.histogram(distintrondist, bins=50, range=(-150, 150))
    distintron_hist = distintron_hist/(distintronsize/1000.)        
    distintron_rand_hist, distintron_edges_rand = np.histogram(rand_distintron, bins=50, range=(-150, 150))
    distintron_rand_hist = distintron_rand_hist/(distintron_rand_size/1000.)       
    ax_distintron.plot(distintron_edges[:-1], distintron_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_distintron.hold(True)
    ax_distintron.plot(distintron_edges_rand[:-1], distintron_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    proxintron_hist, proxintron_edges = np.histogram(proxintrondist, bins=50, range=(-150, 150))
    proxintron_hist = proxintron_hist/(proxintronsize/1000.)        
    proxintron_rand_hist, proxintron_edges_rand = np.histogram(rand_proxintron, bins=50, range=(-150, 150))
    proxintron_rand_hist = proxintron_rand_hist/(proxintron_rand_size/1000.)       
    ax_proxintron.plot(proxintron_edges[:-1], proxintron_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_proxintron.hold(True)
    ax_proxintron.plot(proxintron_edges_rand[:-1], proxintron_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    return



def RNA_position(bedline, GI):
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
    if outfile is None:
        outfile = clustersFA + ".k" + str(k) + ".kmerdiff"
    print "Running Kmer analysis"
    subprocess.call(["perl", (basedir + "/yeolab/Software/generalscripts/kmerdiff.pl"), "-file1", clustersFA, "-file2", backgroundFA, "-k", str(k), "-o", outfile])
    srtout = outfile + ".sort"
    srt = open(srtout, 'w')
    print "Sorting Results..."
    subprocess.call(["perl", (basedir + "/yeolab/Software/generalscripts/sortkmerdiff.pl"), outfile], stdout=srt)
    print "Kmer analysis done, output here: %s" %(srtout)
    srt.close()
    return srtout
    

def run_homer(foreground, background, k = list([5,6,7,8,9]), outloc = os.getcwd()):
    #findMotifs.pl clusters.fa fasta outloc -nofacts p 4 -rna -S 10 -len 5,6,7,8,9 -noconvert -nogo -fasta background.fa
    ks = ",".join(map(str, k))
    print "starting Homer"
    subprocess.call(["findMotifs.pl", foreground, "fasta", outloc, "-p", "4", "-rna", "-S", "20", "-len", ks, "-fasta", background])
    print "Homer Finished, output here: %s" %(outloc)
    return


def write_seqs(outfile, bedtool_list):
    f = open(outfile, "w")
    for bedtool in bedtool_list:
        f.write(open(bedtool.seqfn).read())
    f.close()
    return


def motif_boxplots(kmerloc, filename, klengths, highlight_motifs, subplot=None):

    """
    Make bake boxplots of motif z-scores. you must get kmer z-scores first with run_kmerdiff.  up to 11 motifs can be highlighted
    pass a pylab subplot instance to the kwarg \"subplot\" to attach this to a figure, otherwise it will make its own figure
    """
    colorcycle = ["red", "orange", "green", "blue", "purple", "brown", "black", "pink", "gray", "cyan", "magenta"]
    regions = ["all", "exon", "UTR3", "UTR5", "proxintron", "distintron"]
    formal_labels = ["All Regions", "Exon", "3'UTR", "5'UTR", "Proximal Intron", "Distal Intron"]    
    kmers = {}
    all_kmers = set()
    for region in regions:
        kmers[region] = {}
        for k in klengths:
            f = open(os.path.join(kmerloc, "%s.k%s.%s.kmerdiff.sort" %(filename, str(k), region)))  
            for line in f:
                kmer, val = line.strip().split("\t")
                kmer, val = map(str.strip, [kmer, val])
                all_kmers.add(kmer)
                kmers[region][kmer] = float(val)

    for i, m in enumerate(highlight_motifs):
        highlight_motifs[i] = m.lower().replace("u", "t")#hash the motifs, convert to DNA letters
    ak = np.ndarray(shape=(len(all_kmers), len(regions)))
    all_kmers = list(all_kmers)

    for i, kmer in enumerate(all_kmers):
        for j, region in enumerate(regions):            
            ak[i,j] = kmers[region][kmer]
    
    showme=False
    if subplot is None: #
        showme=True
        x = pylab.figure()
        subplot = x.add_subplot(111)

    subplot.boxplot(ak, vert=False, notch=1, sym='k.',  whis=2)
    subplot.set_yticklabels(formal_labels)
    for i, motif in enumerate(highlight_motifs):
        indices= list()
        for ind, k in enumerate(all_kmers):
            if motif in k:
                indices.append(ind)
        y = map(lambda x: x+1, range(len(regions)))
        for m, ind in enumerate(indices):
            if m ==0:
                label=motif
            else:
                label=None
            subplot.plot(ak[ind,:], y, 'o', color=colorcycle[i], label=label, markersize=10)
    subplot.set_xscale('symlog', linthreshx=10)
    subplot.axvline(x=-4)
    subplot.axvline(x=4)
    

    subplot.legend(frameon=False,loc=0, numpoints=1)
    
    subplot.set_xlabel("Z-score")
    if showme is True:
        pylab.show()    

    return ak, all_kmers
def bedlengths(tool):
    x = list()
    for line in tool:
        x.append(line.length)
    return x



def chop(chr, start,end, wsize=5):
    file = open("phastcons.txt", 'w')
    i = start
    while i < end:
        x = pybedtools.Interval(chr, i, (i+wsize))
        p = get_phastcons(x, species="hg19")
        file.write("\t".join(map(str, [chr, i, (i+wsize-1), p])) + "\n")
        i += wsize
    file.close()
 




def get_phastcons(bedtool, species=None, index=None):
    """
    Get phastcons scores for intervals in a bed tool
    """
    if species is None and index is None:
        print "Error, must select species or index"
    if species is not None and index is None:
        if species == "mm9":
            index= basedir + "/yeolab/Conservation/phastCons/mm9_30way/placental/mm9_phastcons.bw"
        elif species == "hg19":
            index = basedir + "/yeolab/Conservation/phastCons/hg19_46way/placentalMammals/reformat/hg19_phastcons.bw"
    f = open(index, 'r')
    bw = BigWigFile(file=f)

    try:
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
        for i, bedline in enumerate(bedtool):
            data = np.ndarray(len(bedtool))        
            vals = bw.get(bedline.chrom, bedline.start, bedline.stop)
            consvals = list(v[-1] for v in vals)
            if len(consvals) > 0:
                mean_phastcons = np.mean(consvals)
            else:
                mean_phastcons=0
            data[i] = mean_phastcons
    return data




def main(options):
    


    
    #
    from subprocess import Popen, PIPE
    host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
    #print host
    #print mpl.get_backend()

    #print mpl.get_backend()
    
    clusters = options.clusters
    species = options.species
    CLUSTERS = pybedtools.BedTool(clusters)

    clusters = str.replace(clusters, ".BED", "")
    options.k= map(int, options.k)
    outdir = options.outdir

    def make_dir(dir_name):
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)

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


    def fa_file(filename, region = None, fd=fastadir, type= "real"):
        if not os.path.exists(fd):
            raise Exception
        if region is not None:
            x =filename+"."+  region+ "."+ type+ ".fa"
            return os.path.join(fd, x)
        else:
            x = filename+ "."+ type + ".fa"
            return os.path.join(fd, x)

    if options.assign is False:
        try:
            CLUS_regions, sizes, Gsizes = build_assigned_from_existing(assigned_dir, clusters, all_regions, options.nrand)
            print "I used a pre-assigned set of BED files... score!"
        except:
            print "I had problems retreiving region-assigned BED files from %s, i'll rebuild" %(assigned_dir)
            options.assign=True
            
    if options.assign is True:
        print "Assigning Clusters to Genic Regions"
        CLUS_regions, sizes, Gsizes = assign_to_regions(CLUSTERS, species=species, getseq=True, nrand=options.nrand)
        print "Done Assigning"

        print "Saving BED and Fasta Files...",

        sizes_out = open(os.path.join(assigned_dir, "%s.sizes.pickle" %(clusters)), 'w')
        pickle.dump(sizes, file=sizes_out)
        sizes_out.close()    
        Gsizes_out = open(os.path.join(assigned_dir, "Gsizes.pickle"), 'w')
        pickle.dump(Gsizes, file=Gsizes_out)
        Gsizes_out.close()

        for region in all_regions:
            of = clusters + "." + region+ ".real.BED"
            try:
                CLUS_regions[region]['real'].saveas(os.path.join(assigned_dir, of))
            except:
                continue
            for n in range(options.nrand):
                of = clusters + "." + region+ ".rand." + str(n) + ".BED"
                try:
                    CLUS_regions[region]['rand'][n].saveas(os.path.join(assigned_dir, of))
                except:
                    continue
                
        print "done"

        for region in all_regions:
            try:
                real_fa = fa_file(clusters, region=region, type="real")
                rand_fa = fa_file(clusters, region=region, type="random")
                CLUS_regions[region]['real'].save_seqs(real_fa)

                l = list()#list of randoms
                for n in CLUS_regions[region]['rand'].keys():
                    l.append(CLUS_regions[region]['rand'][n])
                write_seqs(rand_fa, l)        
            except:
                continue            
                                   
    print "Counting reads in clusters...",
    reads_in_clusters = 0
    reads_per_cluster = list()
    for cluster in CLUS_regions['all']['real']:
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
    #bamfile = pysam.Samfile(options.bam, 'rb')
    print "Getting total number of reads...",
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
            
    except:
        print "Couldn't find a pickled file, resorting to flagstat for total reads. (this includes intergenic reads)"
        flagstats = pysam.flagstat(options.bam)
        total_reads =int(flagstats[2].split(" ")[0])
        
    print "done, there were %d" %(total_reads)
    print "Gathering bed lengths...",
    cluster_lengths = bedlengths(CLUS_regions['all']['real'])
    print "done"
##     
    mRNA_positions = list()
    premRNA_positions = list()
    intron_positions = list()
    exon_positions = list()
    GENES, Gtypes = build_AS_STRUCTURE_dict(species)
    types = {}
    for type in ["CE:", "SE:", "MXE:", "A5E:", "A3E:"]:
        types[type]=0
    print "locating clusters within genes",
    try:
        for line in (CLUS_regions['all']['real']):
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
                                     
    type_count = [types["CE:"], types["SE:"], types["MXE:"], types["A5E:"], types["A3E:"]]
    Gtype_count = [Gtypes["CE:"], Gtypes["SE:"], Gtypes["MXE:"], Gtypes["A5E:"], Gtypes["A3E:"]]    

    ### write fasta files and run homer and/or kmer analysis if at least one analysis is requested
    if options.reMotif is True:
       
        for region in all_regions:
            try:
                real_fa = fa_file(clusters, region=region, type="real")
                rand_fa = fa_file(clusters, region=region, type="random")
                if options.k is not None:
                    if options.homer is True:
                        region_homer_out = os.path.join(homerout, region)
                        run_homer(real_fa, rand_fa, options.k,  outloc=region_homer_out)
                    for k in options.k:                    
                        kmerfile = clusters + ".k" + str(k) + "." + region + ".kmerdiff"
                        kmerfile = os.path.join(kmerout, kmerfile)
                        kmer_sorted_output = run_kmerdiff(real_fa, rand_fa, outfile=kmerfile, k=k)
            except:
                continue

    motifs = list(options.motif)
    kmer_box_params = [kmerout, clusters, options.k, motifs]

    ###conservation --should use multiprocessing to speed this part up!
    phast_values = list()
    if options.rePhast is False:
        try:
            phast_values = pickle.load(open(os.path.join(misc_dir, "%s.phast.pickle" %(clusters))))
        except:
            options.rePhast =True


    if options.rePhast is True:
        print "Fetching Phastcons Scores...",
        for region in all_regions[1:]:#skip "all" combine them later
            print ("%s..." %(region)),
            try:
                samplesize=1000
                if len(CLUS_regions[region]['real']) > samplesize:
                    R1 = CLUS_regions[region]['real']                
                    # R1 = random.sample(CLUS_regions[region]['real'], samplesize)
                else:
                    R1 = CLUS_regions[region]['real']

                #realPhast = get_phastcons(CLUS_regions[region]['real'], species=options.species)
                print "getting real...",
                realPhast = get_phastcons(R1, species=options.species)
                randPhast=list()
                for i in range(options.nrand):
                    if len(CLUS_regions[region]['rand'][i]) > samplesize:
                        R2 = CLUS_regions[region]['rand'][i]                    
                        #R2 = random.sample(CLUS_regions[region]['rand'][i], samplesize)
                    else:
                        R2 = CLUS_regions[region]['rand'][i]
                    print ("getting rand %d" %(i)),
                    randPhast.extend(get_phastcons(R2, species=options.species).tolist())
                phast_values.append(realPhast)
                phast_values.append(randPhast)
            except:
                continue
        all_real = np.concatenate(phast_values[::2])
        all_rand = np.concatenate(phast_values[1::2])
        phast_values.insert(0,all_rand)
        phast_values.insert(0,all_real)
        pickout = open(os.path.join(misc_dir, "%s.phast.pickle" %(clusters)), 'w')
        pickle.dump(phast_values, file = pickout)
    Zscores = None  #old. remove

    QCfig_params = [reads_in_clusters, (total_reads - reads_in_clusters), cluster_lengths, reads_per_cluster, premRNA_positions, mRNA_positions, exon_positions, intron_positions, Gsizes, sizes, Gtype_count, type_count, Zscores, homerout, kmer_box_params, phast_values]

    pickout = open(os.path.join(outdir, "misc", "%s.qcfig_params.pickle" %(clusters)), 'w')
    pickle.dump(QCfig_params, file = pickout)
    QCfig = CLIP_QC_figure(*QCfig_params)
    fn = clusters + ".QCfig.pdf"
    outFig = os.path.join(outdir, fn)
    QCfig.savefig(outFig)
                    
###
    motifs = list(options.motif)
    motifBASE  = basedir + "/lovci/projects/ucscBED"
    if motifs is not None:
        fig = pylab.figure(figsize=(8.5, 11))
        colors = ["red", "orange", "green", "blue", "purple", "brown", "black", "pink", "gray", "cyan", "magenta"]
        for i, motif in enumerate(motifs):
            mf = "motif_" + motif + ".BED"
            mfgz = "motif_" + motif + ".BED.gz"
            print os.path.join(motifBASE,species,mf)
            motifFILE = None
#            import code
#            code.interact(local=locals())
            if os.path.exists(os.path.join(motifBASE,species, mf)):
                motifFILE = os.path.join(motifBASE,species, mf)
            elif os.path.exists(os.path.join(motifBASE,species, mfgz)):
                motifFILE= os.path.join(motifBASE,species, mfgz)
            else:
                print "MOTIF BED FILE for motif: %s is not available, please build it" %(mf)
                continue
            plot_motif_dist(CLUS_regions, motifFILE, fig, color = colors[i], species=species, slopsize=200)
        pylab.savefig(clusters + ".motif_distribution.pdf")

if __name__== "__main__":
    parser = OptionParser()
    
    parser.add_option("--clusters", dest="clusters", help="BED file of clusters", metavar="BED")
    parser.add_option("--reAssign", dest="assign", action="store_true", default=False, help="re-assign clusters, if not set it will re-use existing assigned clusters") ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
    parser.add_option("--rePhast", dest="rePhast", action="store_true", default=False, help="re-calculate conservation, must have been done before") ##to-do. this should be auto-set if the creation date of "clusters" is after creation date fo assigned files
    parser.add_option("--old_motifs", dest="reMotif", action="store_false", default=True, help="use old motif files")
    parser.add_option("--species", dest="species", help = "genome version")
    parser.add_option("--motif", dest="motif", action="append", help="Files of motif locations", default=None)
    parser.add_option("--homer", dest="homer", action="store_true", default=False)
    parser.add_option("--conservation", dest="cons", action="store_true")
    parser.add_option("--structure", dest="structure", action="store_true")
    parser.add_option("--nrand", dest="nrand", default=3, type="int")
    parser.add_option("--k", dest="k", action="append", default=[6])
    parser.add_option("--outdir", dest="outdir", default=os.getcwd(), help="directory for output, default:cwd")
    parser.add_option("--bam", dest="bam")
    
    (options, args) = parser.parse_args()

    main(options)

