from subprocess import Popen, call, PIPE
import re
import pybedtools

from CLIP_analysis import get_phastcons
host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
if "optiputer" in host or "compute" in host:
    base = "/nas/nas0"
elif "tcc" in host or "triton" in host:
    base = "/projects"    
print "base directory: %s" %(base)


try:
    pybedtools.set_tempdir(base + "/scratch/lovci/pybedtools_tmp")
except:
    try:
        pybedtools.set_tempdir(base + "/lovci/projects/tmp/pybedtools_tmp")
    except:
        print "No valid pybedtools tmp directory"
        exit()
        





RNAhybridcmd = base + "/yeolab/Software/RNAhybrid/bin/RNAhybrid"
RNAhybrid_parser = base + "/lovci/gscripts/parseRNAhybrid.pl"

def RNAhybrid_hits(target, query):
    xi, theta = map(str, [-20, 10])
    
    rh = Popen([RNAhybridcmd, "-d", (xi + "," + theta), "-n", "70", "-m", "10000", "-e" "-40", str(target), str(query)], stdout = PIPE)
    parse = Popen(["perl", RNAhybrid_parser], stdin=rh.stdout, stdout=PIPE)
    rh.stdout.close()
    results = parse.communicate()[0]
    return results


def build_bed(filename, outfile = "test.bed"):
    f = open(filename, 'r')
    
    proxLen = 500
    data = list()

    out = open(outfile, 'w')
    n=300
    for line in f.readlines():
        if n==0:
            continue
        #n=n-1
        if "Exonloc" in line:
            print "headerskip"
            continue
        id, type, wholeLoc, exonLoc  = line.strip().split("\t")

        shortLine = wholeLoc
        #print shortLine

        chrom, start, stop, strand = re.split(r'\W+', wholeLoc)
        strand = wholeLoc[-1]
        exstart, exstop = exonLoc.split('-')
        upstream_intron = start + "-" + str(int(exstart)- 1)
        downstream_intron = str(int(exstop)+1) + "-" + stop
        if strand == "-":
            upstream_intron, downstream_intron = downstream_intron, upstream_intron

        uI_start , uI_stop = upstream_intron.split('-')
        uI_length = int(uI_stop) - int(uI_start)
        dI_start , dI_stop = downstream_intron.split('-')
        dI_length = int(dI_stop) - int(dI_start)
        
        if uI_length < proxLen:
            #print "upstream intron too short, skipping upstream distal class for %s" %(line)
            uI = "\t".join([chrom , uI_start, uI_stop, "%".join(["ui", shortLine, id, type]), "1", strand])
            uIProx = uI.replace("ui", "uiProx")
            out.write("\n".join([uI,  uIProx]) + "\n")           
        else:
            uI = "\t".join([chrom , uI_start, uI_stop, "%".join(["ui", shortLine, id, type]), "1", strand])
            uIDist = "\t".join([chrom , uI_start, str((int(uI_stop)-proxLen)), "%".join(["uiDist", shortLine, id, type]), "1", strand])
            uIProx = "\t".join([chrom , str(int(uI_stop)-proxLen), uI_stop, "%".join(["uiProx", shortLine, id, type]), "1", strand])    
            out.write("\n".join([uI, uIDist, uIProx]) + "\n")

        if dI_length < proxLen:
            #print "downstream intron too short, skipping downstream distal class for %s" %(line)
            dI = "\t".join([chrom , dI_start, dI_stop, "%".join(["di", shortLine, id, type]), "1", strand])
            dIProx = dI.replace("di", "diProx")
            out.write("\n".join([dI, dIProx]) + "\n")
        else:
            dI = "\t".join([chrom , dI_start, dI_stop, "%".join(["di", shortLine, id, type]), "1", strand])
            dIDist = "\t".join([chrom , str(int(dI_start)+proxLen), dI_stop, "%".join(["diDist", shortLine, id, type]), "1", strand])
            dIProx = "\t".join([chrom , dI_start, str(int(dI_start) + proxLen), "%".join(["diProx", shortLine, id, type]), "1", strand])
            
            out.write("\n".join([dI, dIDist, dIProx]) + "\n")
            
    out.close()
    tool  =pybedtools.BedTool(outfile)
    return(tool)





species = "hg19"
genomeFasta = base+ "/yeolab/Genome/ucsc/" + species + "/chromosomes/all.fa"



b = build_bed(base + "/yeolab/Genome/ensembl/AS_STRUCTURE/hg19data4/hg19.internal_exons.with_intron_flanks.table", outfile = "event_detail.BED")

