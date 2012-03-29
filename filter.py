import pybedtools
import subprocess
import matplotlib as mpl
mpl.use('Agg')
import pylab
import os
import numpy as np
pybedtools.set_tempdir('/nas3/lovci/projects/tmp/python_tmp')
from pybedtools import featurefuncs
from pybedtools.featurefuncs import less_than
from pybedtools.featurefuncs import greater_than
#mpl.rcParams['interactive']=True
def bedlengths(tool):
    x = list()
    for line in tool:
        x.append(line.length)
    return x
def bedscores(tool):
    x = list()
    for line in tool:
        x.append(float(line.score))
    return x


def rewrite(A, B=None, outfile=None):
    tmp = outfile+ ".tmp"
    of = open(outfile, 'w')
    ot = open(tmp, 'w')
    x = A.intersect(B, s=True, wao=True)
#    print x.__len__()
    seen = dict()
    n=0
    for i in x:
        #        n=n+1
        #        if n%1000 == 0:
        #            print n
        chr, start, stop, name, score, strand, chr2, start2, stop2, name2, score2, strand2, ov = str(i).strip().split("\t")
        name = chr + "|" + start + "|" + stop + "|" + name + "|" + strand
        if name in seen:
            continue
        else:
            length = int(stop)-int(start)
            if length == int(ov):
                seen[name]=score
            else:
                seen[name]="redo"
#    print "loaded in the events"
    for event in seen.keys():
        chr, start, stop, name, strand = event.split("|")
        score = seen[event]
        line = chr + "\t"+ start + "\t" + stop + "\t" + name + "\t" + score + "\t" + strand + "\n"
        if score is "redo":
            ot.write(line)
            #            print line,
        else:
            of.write(line)
    ot.close()
    ot = open(tmp, 'r')
    of.close()
    of = open(outfile, 'a')
    subprocess.Popen(["./recalculate_scores.pl"], stdout=of, shell=True, stdin=ot).wait()
    ot.close()
    os.remove(tmp)
    of.close()
    return pybedtools.BedTool(outfile).sort().saveas(outfile)



def intersection(A, B=None):
#    tmp = outfile+ ".tmp"
#    of = open(outfile, 'w')
#    ot = open(tmp, 'w')
    less = A.subtract(B, s=True) #without regions that overlap
    more = A.subtract(less, s=True) #only regions that overlap
    return less, more


def bedCoverage(bed, vs=None):
    #    compute the fraction of bases in bed that overlap with vs retun x if x = "x% of bases in bed are overlapped by vs"
    total_length = np.sum(bedlengths(bed), dtype=float)
    minus_bed = bed.subtract(vs)
    minus_length = np.sum(bedlengths(minus_bed))
    diff_length = total_length - minus_length
    fraction = 100.*(diff_length/total_length)
    return(fraction)

def dumpInFile(array,file=None):
    f = open(file, 'w')
    for i in array:
        for j in i:
            f.write(str(j) + "\t")
        f.write("\n")
    f.close()

def dumpInFile2(vector,file=None):
    f = open(file, 'w')
    for i in vector:
        f.write(str(i) + "\t")
    f.write("\n")
    f.close()



def main():
    length_cutoff = 2000  #remove regions shorter than this.
    UTR3 = pybedtools.BedTool("UTR3_hg19_frea_sorted.withscore")
    UTR5 = pybedtools.BedTool("UTR5_hg19_frea_sorted.withscore")
    exon = pybedtools.BedTool("exon_hg19_frea_sorted.withscore")
    proxintron = pybedtools.BedTool("proxintron500_hg19_frea_sorted.withscore")
    distintron = pybedtools.BedTool("distintron500_hg19_frea_sorted.withscore")

    all_regions = (["exon", "UTR3", "UTR5", "proxintron", "distintron"])
    all_bedtracks = ([exon, UTR3, UTR5, proxintron, distintron])
    low = pybedtools.BedTool("conserved_regions.hg19.margin_2.0-3.BED")
    moderate2 = pybedtools.BedTool("conserved_regions.hg19.margin_2.3-9.BED")
    moderate = pybedtools.BedTool("conserved_regions.hg19.margin_2.3-9.BED")
    high = pybedtools.BedTool("conserved_regions.hg19.margin_2.6-9.BED")
    ultra = pybedtools.BedTool("conserved_regions.hg19.margin_2.9-10.BED")


    high = high.subtract(ultra, s=True)
    tmp = pybedtools.BedTool((str(high)+str(ultra)), from_string=True)
    moderate = moderate.subtract(tmp, s=True)
    moderate2 = moderate2.subtract(ultra, s=True)
    tmp = pybedtools.BedTool((str(high)+str(ultra)+str(moderate) +str(moderate2)), from_string=True)
    low= low.subtract(tmp,s=True)
    tmp.delete_temporary_history(ask=False)
    del tmp

    ultra_pf = ultra.merge(s=True, nms=True, scores="mean").saveas("ultra.prefilter.BED")
    high_pf = high.merge(s=True, nms=True, scores="mean").saveas("high.prefilter.BED")
    moderate_pf = moderate.merge(s=True, nms=True, scores="mean").saveas("moderate.prefilter.BED")
    moderate2_pf =moderate2.merge(s=True, nms=True, scores="mean").saveas("moderate2.prefilter.BED")
    low_pf = low.merge(s=True, nms=True, scores="mean").saveas("low.prefilter.BED")


## #    return ultra, ultra_pf, high, high_pf, moderate, moderate_pf, moderate2, moderate2_pf, low, low_pf




    #recalculate phastcons scores for adjusted regions
    ultra = rewrite(ultra_pf, B=ultra, outfile="ultra.prefilter.fix.BED") 
    high = rewrite(high_pf, B=high, outfile="high.prefilter.fix.BED")
    moderate = rewrite(moderate_pf, B=moderate, outfile="moderate.prefilter.fix.BED")
    moderate2 = rewrite(moderate2_pf, B=moderate2, outfile="moderate2.prefilter.fix.BED")
    low = rewrite(low_pf, B=low, outfile="low.prefilter.fix.BED")

    print "re-examined phastcons scores"
    
    del ultra_pf
    del high_pf
    del moderate_pf
    del moderate2_pf
    del low_pf
    



    low = pybedtools.BedTool("low.prefilter.fix.BED")
    moderate = pybedtools.BedTool("moderate.prefilter.fix.BED")
    moderate2 = pybedtools.BedTool("moderate2.prefilter.fix.BED")
    high = pybedtools.BedTool("high.prefilter.fix.BED")
    ultra = pybedtools.BedTool("ultra.prefilter.fix.BED")
    first_counts = [low.__len__(), moderate2.__len__(), moderate.__len__(),  high.__len__(), ultra.__len__()]
    dumpInFile2(first_counts, file="first_counts.txt")

##     #plot length distributions
    low.lengths = np.log10(bedlengths(low))
    moderate.lengths = np.log10(bedlengths(moderate))
    moderate2.lengths = np.log10(bedlengths(moderate2))
    high.lengths = np.log10(bedlengths(high))
    ultra.lengths = np.log10(bedlengths(ultra))
    lenfigure = pylab.figure()
    pylab.boxplot([low.lengths, moderate.lengths, high.lengths, ultra.lengths])
    pylab.axhline(y=np.log10(length_cutoff), color='r')
    lenfigure.savefig("region_lengths4.pdf")        
    
    lenfigure2 = pylab.figure()
    pylab.boxplot([low.lengths, moderate2.lengths, ultra.lengths])
    pylab.axhline(y=np.log10(length_cutoff), color='r')
    lenfigure2.savefig("region_lengths3.pdf")        


    print "plotted initial length distributions"
    #filter conserved regions < length_cutoff
    
    low = pybedtools.BedTool(str(low.filter(less_than, length_cutoff)), from_string=True)
    moderate = pybedtools.BedTool(str(moderate.filter(less_than, length_cutoff)), from_string=True)
    moderate2 = pybedtools.BedTool(str(moderate2.filter(less_than, length_cutoff)), from_string=True)
    high = pybedtools.BedTool(str(high.filter(less_than, length_cutoff)), from_string=True)
    ultra = pybedtools.BedTool(str(ultra.filter(less_than, length_cutoff)), from_string=True)

    print "removed long ones"
    rmlong_counts = [low.__len__(), moderate2.__len__(),moderate.__len__(),  high.__len__(), ultra.__len__()]
    dumpInFile2(rmlong_counts, file="rmlong_counts.txt")
    #plot conservation scores
    low.cons = bedscores(low)
    moderate.cons = bedscores(moderate)
    moderate2.cons = bedscores(moderate2)
    high.cons = bedscores(high)
    ultra.cons = bedscores(ultra)
    consfigure = pylab.figure()
    pylab.boxplot([low.cons, moderate.cons, high.cons, ultra.cons])
    consfigure.savefig("region_conservation4.pdf")

    consfigure2 = pylab.figure()
    pylab.boxplot([low.cons, moderate2.cons, ultra.cons])
    consfigure2.savefig("region_conservation3.pdf")
    print "plotted conservation scores"

    #lists/containers
    all_cr = ([low, moderate2, moderate, high, ultra])
#    return all_cr
    all_cr_names = (["low", "moderate2", "moderate", "high", "ultra"])
    bed_dict = dict()
    bed_dict['low'] = dict()
    bed_dict['moderate'] = dict()
    bed_dict['moderate2'] = dict()
    bed_dict['high'] = dict()
    bed_dict['ultra'] = dict()

    #allocate conserved regions to genic regions
    for i in range(5):
        for j, region in enumerate(all_regions):
            name = all_cr_names[i] + "_" + region + ".BED"
            no_overlapping, only_overlapping =intersection(all_cr[i], B = all_bedtracks[j])  #portions of the regions that overlap with a genic region
            print str(np.sum(bedlengths(all_cr[i]))) + "\t" + str(np.sum(bedlengths(no_overlapping))) + "\t" + str(np.sum(bedlengths(only_overlapping)))
            bed_dict[all_cr_names[i]][region] = only_overlapping
            all_cr[i] = no_overlapping
            only_overlapping.saveas(name)

    #count the elements in each category
    first_counts_regions = np.ndarray((5,5))


    #remove regions that overlap with annotated miRs
    mir_counts = np.ndarray((5,5))
    mirs = pybedtools.BedTool("/nas3/yeolab/Genome/miRNA/mirbase16/genomes/hsa.BED")
    for i, type in enumerate(all_cr_names):
        for j, region in enumerate(all_regions):
            first_counts_regions[i,j] = bed_dict[type][region].__len__()
            bed_dict[type][region] = bed_dict[type][region].intersect(mirs, v=True)
            mir_counts[i,j] = bed_dict[type][region].__len__()

    dumpInFile(first_counts_regions, file="first_counts_regions.txt")
    dumpInFile(mir_counts, "mir_counts.txt")


    #remove regions that overlap with rmsk elements
    repeats = pybedtools.BedTool("/nas3/yeolab/Genome/ucsc/hg19/database/rmsk.BED")
    repeat_counts = np.ndarray((5,5))
    for i, type in enumerate(all_cr_names):
        for j, region in enumerate(all_regions):

    #        repeat_counts[i,j] = bed_dict[type][region].intersect(repeats).__len__()
            bed_dict[type][region] = bed_dict[type][region].intersect(repeats, v=True)
            repeat_counts[i,j]=bed_dict[type][region].__len__()

    dumpInFile(repeat_counts, file="repeat_counts.txt")


    #remove regions that overlap with TF binding sites from ENCODE
    ENCODE = pybedtools.BedTool("/nas3/yeolab/Genome/ucsc/hg19/database/ENCODE/wgEncodeRegTfbsClustered.bed")
    encode_counts = np.ndarray((5,5))
    for i, type in enumerate(all_cr_names):
        for j, region in enumerate(all_regions):
    #        encode_counts[i,j] = bed_dict[type][region].intersect(ENCODE).__len__()
            name = type + "_" + region + ".filtered_normsk.BED"
            bed_dict[type][region] = bed_dict[type][region].intersect(ENCODE, v=True).saveas(name)
            encode_counts[i,j] = bed_dict[type][region].__len__()

    dumpInFile(encode_counts, file="encode_counts.txt")

    stats = open("loss_stats.txt", 'w')
    for row in range(5):
        first = np.sum(first_counts_regions[row,range(5)])
        mirs = (np.sum(mir_counts[row,range(5)]))
        reps = (np.sum(repeat_counts[row,range(5)]))
        encode = (np.sum(encode_counts[row,range(5)]))
        mirloss = 100* (first - mirs)/ first
        reploss = 100* (mirs - reps) / mirs
        encodeloss = 100 * (reps - encode) /reps
        line =  all_cr_names[row] + "\t" + str(first) + "\t" + str(mirloss) + "\t" + str(reploss) + "\t" + str(encodeloss) + "\t" + str(encode) + "\n"
        stats.write(line)
    stats.close()

    # plot length distributions
    low.lengths = np.log10(bedlengths(low))
    moderate.lengths = np.log10(bedlengths(moderate))
    moderate2.lengths = np.log10(bedlengths(moderate2))
    high.lengths = np.log10(bedlengths(high))
    ultra.lengths = np.log10(bedlengths(ultra))
    lenfigure = pylab.figure()
    pylab.boxplot([low.lengths, moderate2.lengths, moderate.lengths, high.lengths, ultra.lengths])
    pylab.axhline(y=np.log10(length_cutoff), color='r')
    lenfigure.savefig("region_lengths_after.eps")

    CLIPdict = dict()
    CLIPdict["fox2Brain"] = pybedtools.BedTool("/nas3/lovci/projects/FOX2/FOX2_human_brain/CLIP/analysis_gsnap/mike_clusters/sig_clusters.BED")
    CLIPdict["AGO"] = pybedtools.BedTool("CLIPs/hg19/AGO_ALL_par_clusters.hg19.BED")
    CLIPdict["PUM2"] = pybedtools.BedTool("CLIPs/hg19/PUM2_clusters.hg19.BED")
    CLIPdict["QKI"]=pybedtools.BedTool("CLIPs/hg19/QKI_clusters.hg19.BED")
    CLIPdict["TNRC6"]= pybedtools.BedTool("CLIPs/hg19/TNRC6_ALL_clusters.hg19.BED")
    CLIPdict["IGF2BP1"] = pybedtools.BedTool("CLIPs/hg19/IGF2BP1_all.hg19.BED")
    CLIPdict["IGF2BP2"]= pybedtools.BedTool("CLIPs/hg19/IGF2BP2_all.hg19.BED")
    CLIPdict["IGF2BP3"] = pybedtools.BedTool("CLIPs/hg19/IGF2BP3_all.hg19.BED")
    CLIPdict["TDP43"] = pybedtools.BedTool("CLIPs/hg19/TDP43_MP51.hg19.BED")
    CLIPdict["TLS"] = pybedtools.BedTool("CLIPs/hg19/TLS_all.hg19.BED")
    CLIPdict["HA1"] = pybedtools.BedTool("CLIPs/hg19/HnRNPA1all_comb_trim_ingenes_clusters_hg19150.bed")
    CLIPdict["HA2B1"] = pybedtools.BedTool("CLIPs/hg19/HnRNPA2B1_comb_trim_ingenes_clusters_hg19150.bed")
    CLIPdict["HD"] = pybedtools.BedTool("CLIPs/hg19/HnRNPD_comb_trim_ingenes_clusters_hg19150.bed")
    CLIPdict["HF"] = pybedtools.BedTool("CLIPs/hg19/HnRNPF_comb_trim_ingenes_clusters_hg19150.bed")
    CLIPdict["HM"] = pybedtools.BedTool("CLIPs/hg19/HnRNPM_comb_trim_ingenes_clusters_hg19150.bed")
    CLIPdict["HG"] = pybedtools.BedTool("CLIPs/hg19/HnRNPG_comb_trim_ingenes_clusters_hg19150.bed")
    CLIPdict["HU"] = pybedtools.BedTool("CLIPs/hg19/HnRNPU_comb_trim_ingenes_clusters_hg19150.bed")
    CLIPdict["HH"] = pybedtools.BedTool("CLIPs/hg19/HnRNPHall_comb_trim_ingenes_clusters_hg19150.bed")
    CLIPdict["LIN28ES"] = pybedtools.BedTool("CLIPs/hg19/Lin28ES.hg19.BED")
    CLIPdict["LIN28293T"] = pybedtools.BedTool("CLIPs/hg19/Lin28_293T.hg19.BED")

    CLIPdict["ALL"]= pybedtools.BedTool("CLIPs/hg19/ALL_CLIPs.BED")

    f2B =np.ndarray((5,5))
    ago = np.ndarray((5,5))
    pum2 = np.ndarray((5,5))
    qki= np.ndarray((5,5))
    tnrc6 = np.ndarray((5,5))
    igf2bp1 =  np.ndarray((5,5))
    igf2bp2 =  np.ndarray((5,5))
    igf2bp3 =  np.ndarray((5,5))
    all =  np.ndarray((5,5))
    tdp43 = np.ndarray((5,5))
    tls = np.ndarray((5,5))
    ha1 = np.ndarray((5,5))
    ha2b1 = np.ndarray((5,5))
    hd = np.ndarray((5,5))
    hf = np.ndarray((5,5))
    hm = np.ndarray((5,5))
    hg = np.ndarray((5,5))
    hu = np.ndarray((5,5))
    hh = np.ndarray((5,5))
    le = np.ndarray((5,5))
    l2 = np.ndarray((5,5))



    for i, type in enumerate(all_cr_names):
        for j, region in enumerate(all_regions):
            f2B[i,j]=bedCoverage(bed_dict[type][region], vs=CLIPdict["fox2Brain"])
            ago[i,j]=bedCoverage(bed_dict[type][region], vs=CLIPdict["AGO"])
            pum2[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["PUM2"])
            qki[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["QKI"])
            tnrc6[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["TNRC6"])
            igf2bp1[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["IGF2BP1"])
            igf2bp2[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["IGF2BP2"])
            igf2bp3[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["IGF2BP3"])
            all[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["ALL"])
            tdp43[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["TDP43"])
            tls[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["TLS"])
            ha1[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["HA1"])
            ha2b1[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["HA2B1"])
            hd[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["HD"])
            hf[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["HF"])
            hm[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["HM"])
            hu[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["HU"])
            hh[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["HH"])
            le[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["LIN28ES"])
            l2[i,j] = bedCoverage(bed_dict[type][region], vs=CLIPdict["LIN28293T"])


    dumpInFile(f2B, file="fox2.cts.txt")
    dumpInFile(ago, file="ago.cts.txt")
    dumpInFile(pum2, file="pum2.cts.txt")
    dumpInFile(qki, file="qki.cts.txt")
    dumpInFile(tnrc6, file="tnrc6.cts.txt")
    dumpInFile(igf2bp1, file="igf2bp1.cts.txt")
    dumpInFile(igf2bp2, file="igf2bp2.cts.txt")
    dumpInFile(igf2bp3, file="igf2bp3.cts.txt")
    dumpInFile(all, file="all.cts.txt")
    dumpInFile(tdp43, file="tdp43.cts.txt")
    dumpInFile(tls, file="tls.cts.txt")
    dumpInFile(ha1, file="hnrnpa1.cts.txt")
    dumpInFile(ha2b1, file="hnrnpa2b1.cts.txt")
    dumpInFile(hd, file="hnrnpd.cts.txt")
    dumpInFile(hf, file="hnrnpf.cts.txt")
    dumpInFile(hm, file="hnrnpm.cts.txt")
    dumpInFile(hu, file="hnrnpu.cts.txt")
    dumpInFile(hh, file="hnrnph.cts.txt")
    dumpInFile(le, file="lin28es.cts.txt")
    dumpInFile(l2, file="lin28293t.cts.txt")


if __name__ == '__main__':
    print main()
