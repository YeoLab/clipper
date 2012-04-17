import re
import pysam
from seqTools import *
import sys
import os
from subprocess import Popen, call, PIPE
from optparse import OptionParser, SUPPRESS_HELP
from numpy import *
from deap import dtm
import random

host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
if "optiputer" in host or "compute" in host:
    basedir = "/nas/nas0"
elif "tcc" in host or "triton" in host:
    basedir = "/projects"
else:
    print "Where am I?"
    raise Exception

def assign_reads(gene, splicedict=None, bam_file=None, alignment_slop=10, flip=True, splicetypes=None):
    if splicedict is None or bam_file is None:
        raise Exception
    bam_fileobj = pysam.Samfile(bam_file, 'rb')
    data = {}
    chrom = splicedict[gene]["chromosome"]
    strand = splicedict[gene]["strand"]
    tx_start = splicedict[gene]["tx_start"]
    tx_end = splicedict[gene]["tx_end"]
    signstrand = None
    if flip is not None:
        if flip is True:
            usestrand = strand * -1
        else:
            usestrand = strand
        if usestrand == 1:
            signstrand = "+"
        elif usestrand == -1:
            signstrand = "-"
    subset_reads = bam_fileobj.fetch(reference=chrom, start=tx_start,end=tx_end)


    wig, jxns, nrCounts, readLengths, reads = readsToWiggle_pysam(subset_reads, (tx_start-1000), (tx_end+1000), keepstrand=signstrand)
    data["descriptor"] = gene
    import code
    code.interact(local=locals())
    if "SE" in splicedict[gene] and "SE" in splicetypes:
        data["SE"] = {}
        for loc in splicedict[gene]["SE"]:
            #rangestart = splicedict[gene]["SE"][loc]["rangestart"]
            #rangeend = splicedict[gene]["SE"][loc]["rangeend"]                
            data["SE"][loc] = {}
            data["SE"][loc]["IN"] = 0
            data["SE"][loc]["EX"] = 0
            #import code
            #code.interact(local=locals())
            for structure in splicedict[gene]["SE"][loc]["IN"]:
                if structure.startswith("j"):
                    structurestrip = structure.lstrip("j")
                    if structurestrip in jxns:
                        data["SE"][loc]["IN"] += jxns[structurestrip]
                elif structure.startswith("b"):
                    exstart, exstop = map(int, structure.lstrip("b").split("-"))
                    for position in range(exstart, (exstop+1)):
                        if position in reads:
                            for read_end in reads[position]:
                                if read_end <= exstop:
                                    data["SE"][loc]["IN"] += reads[position][read_end]
                                    
                    #for read in reads:
                    #    rstart, rstop = map(int, read.split("-"))
                    #    if rstart >= (exstart-alignment_slop) and rstop <= (exstop + alignment_slop):
                    #        data["SE"][loc]["IN"] += reads[read]
                    #    else:
                    #        pass
            for structure in splicedict[gene]["SE"][loc]["EX"]:
                if structure.startswith("j"):
                    structurestrip = structure.lstrip("j")
                    if structurestrip in jxns:
                        data["SE"][loc]["EX"] += jxns[structurestrip]

    if "MXE" in splicedict[gene] and "MXE" in splicetypes:
        data["MXE"] = {}
        for loc in splicedict[gene]["MXE"]:
            #rangestart = splicedict[gene]["SE"][loc]["rangestart"]
            #rangeend = splicedict[gene]["SE"][loc]["rangeend"]                
            data["MXE"][loc] = {}
            data["MXE"][loc]["A"] = 0
            data["MXE"][loc]["B"] = 0
            #import code
            #code.interact(local=locals())
            for structure in splicedict[gene]["MXE"][loc]["A"]:
                if structure.startswith("j"):
                    structurestrip = structure.lstrip("j")
                    if structurestrip in jxns:
                        data["MXE"][loc]["A"] += jxns[structurestrip]
                elif structure.startswith("b"):
                    exstart, exstop = map(int, structure.lstrip("b").split("-"))
                    for position in range(exstart, (exstop+1)):
                        if position in reads:
                            for read_end in reads[position]:
                                if read_end <= exstop:
                                    data["MXE"][loc]["A"] += reads[position][read_end]

            for structure in splicedict[gene]["MXE"][loc]["B"]:
                if structure.startswith("j"):
                    structurestrip = structure.lstrip("j")
                    if structurestrip in jxns:
                        data["MXE"][loc]["B"] += jxns[structurestrip]
                elif structure.startswith("b"):
                    exstart, exstop = map(int, structure.lstrip("b").split("-"))
                    for position in range(exstart, (exstop+1)):
                        if position in reads:
                            for read_end in reads[position]:
                                if read_end <= exstop:
                                    data["MXE"][loc]["B"] += reads[position][read_end]                                    
                                    

    return data

def retrieve_splicing(species):
    if species == "hg19":
        chrs = map(str,range(1,23)) #1-22
        chrs.append("X")
        chrs.append("Y")
    elif species == "mm9":
        chrs = map(str,range(1,20)) #1-19
        chrs.append("X")
        chrs.append("Y")        
    info = dict()
    for chr in chrs:
        ASfile = basedir + "/yeolab/Genome/ensembl/AS_STRUCTURE/" + species + "data4/" + species + ".tx." + chr + ".AS.STRUCTURE.flanks"
        f = open(ASfile, "r")
        annotline = f.next()
        eof = False
        while not eof:
            blank, gene, chromosome, transcripts, d2, d3, d4, strand, numex, exonloc, intronloc, exonlen, intronlen, asSplicingType, locSplicingType = annotline.strip().split("\t")
            info[gene] = {}
            info[gene]["chromosome"] = "chr" + chr
            info[gene]["strand"] = int(strand)
            info[gene]["tx_start"] = min(map(int, exonloc.rstrip("|").replace("|", "-").split("-")))
            info[gene]["tx_end"] = max(map(int, exonloc.rstrip("|").replace("|", "-").split("-")))
            line = f.next().strip()
            if line.startswith("<"):
                while True:
                    try:
                        splicing_labels = f.next().strip()
                    except:
                        eof = True
                        break
                    if splicing_labels.startswith(">"):
                        annotline = splicing_labels
                        break
                    loc, splicingType, inorout, labels = splicing_labels.split("\t")

                    if splicingType is "OV" or splicingType is "RI":
                        continue                    # skip RI and OV... not well defined yet

                    if not splicingType in info[gene]:
                        info[gene][splicingType] = {}
                    
#I know the "pass" statements are verbose.
                    if "SE" in splicingType:
                        import code
                        code.interact(local=locals())
                        if not loc in info[gene][splicingType]:
                            info[gene][splicingType][loc] = {}
                            info[gene][splicingType][loc]["rangestart"] = 10000000000000
                            info[gene][splicingType][loc]["rangeend"] = -10000000000000
                            info[gene][splicingType][loc]["bedTrack"] = str()

                            info[gene][splicingType][loc]["IN"] = {}
                            info[gene][splicingType][loc]["EX"] = {}

                            
                        seen = {}
                        versions = labels.rstrip("|").split("|")
                        for v in versions:
                            transcript, vstrand, locUp, locIn, locDown, evidence = v.split(":")

                            if int(vstrand) == -1:
                                locDown, locUp = locUp, locDown


                            locUpx, locUpy = map(str, locUp.split("-"))
                            locDownx, locDowny = map(str, locDown.split("-"))                                                                

                            info[gene][splicingType][loc]["rangestart"] = min(info[gene][splicingType][loc]["rangestart"], (int(locUpx)+1))
                            info[gene][splicingType][loc]["rangeend"] = max(info[gene][splicingType][loc]["rangeend"], (int(locDowny)+1))
                            info[gene][splicingType][loc]["UpExon"]=locUp
                            info[gene][splicingType][loc]["DownExon"]=locDown
                                                                                   
                            if "IN" in inorout:
                                exstart, exstop = map(str, locIn.split("-"))                                
                                try:
                                    info[gene][splicingType][loc][inorout]["b" + locIn] += 1
                                except:
                                    info[gene][splicingType][loc][inorout]["b" + locIn] = 1


                                try:
                                    info[gene][splicingType][loc][inorout]["jup" + locUpy + ":" + str(int(exstart)+1)] += 1#upstream jxn
                                except:
                                    info[gene][splicingType][loc][inorout]["jup" + locUpy + ":" + str(int(exstart)+1)] = 1#upstream jxn
                                try:
                                    info[gene][splicingType][loc][inorout]["jdn" + exstop + ":" + str(int(locDownx)+1)] += 1#dnstream jxn
                                except:
                                    info[gene][splicingType][loc][inorout]["jdn" + exstop + ":" + str(int(locDownx)+1)] = 1 #dnstream jxn
                            else:
                                try:
                                    info[gene][splicingType][loc][inorout]["jex" + locUpy + ":" + str(int(locDownx)+1)] += 1
                                except:
                                    info[gene][splicingType][loc][inorout]["jex" + locUpy + ":" + str(int(locDownx)+1)] = 1
                        if int(strand) == 1:
                            signstrand = "+"
                        else:
                            signstrand = "-"
                        #refine bed track
                        info[gene][splicingType][loc]["bedTrack"] = "\t".join([("chr" + chr), str(info[gene][splicingType][loc]["rangestart"]), str(info[gene][splicingType][loc]["rangeend"]), (gene), "1", signstrand])
                    
                    elif "MXE" in splicingType:
                        #define A or B with varied 5' and 3' exon termini
                        if not loc in info[gene][splicingType]:
                            info[gene][splicingType][loc] = {}
                            info[gene][splicingType][loc]["A"] = {}
                            info[gene][splicingType][loc]["B"] = {}
                            info[gene][splicingType][loc]["rangestart"] = 100000000000
                            info[gene][splicingType][loc]["rangeend"] = -100000000000                            

                        versions = labels.rstrip("|").split("|")
                        for v in versions:
                            transcript, vstrand, locUp, locIn, locDown, evidence = v.split(":")
                            if int(vstrand) == -1:
                                locDown, locUp = locUp, locDown
                                pass

                            locUpx, locUpy = map(str, locUp.split("-"))
                            locDownx, locDowny = map(str, locDown.split("-"))                                                                
                            info[gene][splicingType][loc]["rangestart"] = min(info[gene][splicingType][loc]["rangestart"], int(locUpx))
                            info[gene][splicingType][loc]["rangeend"] = max(info[gene][splicingType][loc]["rangeend"], int(locDowny))                            

                            if "IN" in inorout:
                                exstart, exstop = map(str, locIn.split("-"))                                
                                try:
                                    info[gene][splicingType][loc]["A"]["b" + locIn] += 1#body
                                except:
                                    info[gene][splicingType][loc]["A"]["b" + locIn] = 1#body
                                exstart, exstop = locIn.split("-")
                                try:
                                    info[gene][splicingType][loc]["A"]["jup" + locUpy + ":" + str(int(exstart)+1)] += 1#upstream jxn
                                except:
                                    info[gene][splicingType][loc]["A"]["jup" + locUpy + ":" + str(int(exstart)+1)] = 1 #upstream jxn
                                try:
                                    info[gene][splicingType][loc]["A"]["jdn" + exstop + ":" + str(int(locDownx)+1)] +=1 #dnstream jxn
                                except:
                                    info[gene][splicingType][loc]["A"]["jdn" + exstop + ":" + str(int(locDownx) +1)] =1 #dnstream jxn
                                pass
                            else:
                                #simply the other exon, but IN/EX designations like SE persist.
                                exstart, exstop = map(str, locIn.split("-"))
                                try:
                                    info[gene][splicingType][loc]["B"]["b" + locIn] += 1#body
                                except:
                                    info[gene][splicingType][loc]["B"]["b" + locIn] = 1#body
                                exstart, exstop = locIn.split("-")
                                try:
                                    info[gene][splicingType][loc]["B"]["jup" + locUpy + ":" + str(int(exstart)+1)] +=1#upstream jxn
                                except:
                                    info[gene][splicingType][loc]["B"]["jup" + locUpy + ":" + str(int(exstart)+1)] =1#upstream jxn
                                try:
                                    info[gene][splicingType][loc]["B"]["jdn" + exstop + ":" + str(int(locDownx)+1)] +=1#dnstream jxn
                                except:
                                    info[gene][splicingType][loc]["B"]["jdn" + exstop + ":" + str(int(locDownx)+1)] =1#dnstream jxn
                        if int(strand) == 1:
                            signstrand = "+"
                        else:
                            signstrand = "-"

                            # refine bed track
                           
                        info[gene][splicingType][loc]["bedTrack"] = "\t".join([("chr" + chr), str(info[gene][splicingType][loc]["rangestart"]), str(info[gene][splicingType][loc]["rangeend"]), (gene), "1", signstrand])                                    

                    elif "A5E" in splicingType or "A3E" in splicingType:
                        #not working yet... DANGER!
                        continue
                        if not loc in info[gene][splicingType]: 
                            info[gene][splicingType][loc] = {}                           
                            info[gene][splicingType][loc]['jxns'] = {}
                            info[gene][splicingType][loc]["rangestart"] = 100000000000
                            info[gene][splicingType][loc]["rangeend"] = -100000000000                            

                        versions = labels.rstrip("|").split("|")
                        for v in versions:
                            transcript, vstrand, locUp, locIn, locDown, evidence = v.split(":")
                            if int(vstrand) == -1:
                                locDown, locUp = locUp, locDown
                            info[gene][splicingType][loc]["rangestart"] = min(info[gene][splicingType][loc]["rangestart"], int(locUp))
                            info[gene][splicingType][loc]["rangeend"] = max(info[gene][splicingType][loc]["rangeend"], int(locDown))                                                        
                            exstart, exstop = locIn.split("-")                            
                            jxns = [("jup" + locUp + ":" + str(int(exstart)+1)), ("jdn" + exstop + ":" + str(int(exstop)+1))]
                            for jxn in jxns:
                                try:
                                    info[gene][splicingType][loc]['jxns'][jxn] +=1
                                except:
                                    info[gene][splicingType][loc]['jxns'][jxn] =1

    return info
            

def main(options):
    species = options.species
    bamfile = options.bam
    splicetypes = options.splicetypes
    if splicetypes is None:
        print "you must specify what type of splicing to examine: SE/MXE are implemented now"

    splicing = retrieve_splicing(species)

    if options.gene is not None:
        genes = options.gene
    else:
        genes = splicing.keys()
    
    if options.maxgenes is not None:
        if not options.maxgenes > len(genes):
            genes = random.sample(genes, options.maxgenes)

    for gene in genes:
        x = assign_reads(gene, splicedict=splicing, bam_file=bamfile)
    #data = dtm.map(assign_reads, genes, splicedict=splicing, bam_file=bamfile, splicetypes = splicetypes)
    st = "_".join(splicetypes)
    if options.outfile is None:
        outfile = os.path.join(options.prefix, (bamfile.replace(".bam", ".splices.pickle") + "." + st))
    else:
        outfile = options.outfile
    pickle.dump(data, file=open(outfile, 'w'))

    #AS, splicingTypes  = bulid_AS_STRUCTURE_dict("mm9")
    #SE_list = list()
    #find SE exons, convert to a usable format
    #for gene in AS:
    #    for exonSplicingType in AS[gene]['splicingTypes']:
    #        isoforms
    #        if "SE" in exonSplicingType:
    #other exon splicingTypes...

    
if __name__ == "__main__":
    usage = "python splicing.py --bam <bamfile> --species <species>"
    description = "Given a bam file, count reads to isoforms. comparisons come later"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--bam", '-b', dest="bam", help = "bam file")
    parser.add_option("--species", dest="species")
    parser.add_option("--outfile", dest="outfile", default=None)
    parser.add_option("--gene", dest="gene", default=None, action="append", type="str")
    parser.add_option("--maxgenes", dest="maxgenes", type="int", default=None)    
    parser.add_option("--start", dest="start", action="store_true", default=False, help=SUPPRESS_HELP)
    parser.add_option("--prefix", dest="prefix", default=os.getcwd(), help="output location")
    parser.add_option("--job_name", dest="job_name", default="splice", help="job name")
    parser.add_option("--processors",  dest="np", type="int", default=32, help="number of processors to use")
    parser.add_option("--notify",  dest="notify", default=None, help="email")
    parser.add_option("--wait_to_exit", dest="wait", default=False, action ="store_true")
    parser.add_option("--splicetypes", dest="splicetypes", default=None, action="append")
    (options,args) = parser.parse_args()
    main(options)
    exit()
    if options.start is True:
        dtm.start(main, options)
    else:

        if not os.path.exists(options.prefix):
            try:
                os.mkdir(options.prefix)
            except:
                raise Exception
        scriptName = os.path.join(options.prefix, options.job_name+".runme.sh")
        runerr = os.path.join(options.prefix, options.job_name+ ".err")
        runout = os.path.join(options.prefix, options.job_name+ ".out")
        shScript = open(scriptName, 'w')
        host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()

        if "optiputer" in host or "compute" in host:
            shScript.write("#!/bin/bash\n#$ -N %s\n#$ -S /bin/bash\n#$ -V\n#$ -pe mpi %d\n#$ -cwd\n#$ -o %s\n#$ -e %s\n#$ -l bigmem\n" %(options.job_name, options.np, runout, runerr))
            if options.notify is not None:
                shScript.write("#$ -notify\n#$ -m abe\n#$ -M %s\n" %(options.notify))
            
            shScript.write("/opt/openmpi/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines python %s --start\n" %(" ".join(sys.argv)))

            
        elif "tcc" in host or "triton" in host:
            nnodes = 1
            if int(options.np)%8 == 0:
                nnodes = int(options.np)/8
            else:

                nnodes = (int(options.np)/8)+1
                np = nnodes*8
                print "You should have used a number of processors that is divisible by 8.  You tried %d and I'll actually use %d." %(options.np, np)
            shScript.write("#!/bin/sh\n#PBS -N %s\n#PBS -o %s\n#PBS -e %s\n#PBS -V\n#PBS -S /bin/sh\n#PBS -l nodes=%d:ppn=8\n#PBS -q batch\n#PBS -l walltime=00:50:00\n" %(options.job_name, runout, runerr, nnodes))
                
            if options.notify is not None:
                shScript.write("#PBS -m abe\n#PBS -M %s\n" %(options.notify))


            shScript.write("let np=$PBS_NUM_NODES*$PBS_NUM_PPN\ncd $PBS_O_WORKDIR\n"); 
            shScript.write("/opt/openmpi/bin/mpirun --mca btl_tcp_if_include myri0 -v -machinefile $PBS_NODEFILE -np $np python %s --start\n" %(" ".join(sys.argv)))
            if options.wait is True:
                shScript.write("x=$PBS_JOBID\nJOB_ID=`echo $x | perl -lane '@a = split(/\./,$_); print $a[0]'`\n")
        if options.wait is True:
            shScript.write("touch $JOB_ID.JOBDONE\n")
        shScript.close()                                             
        jobID = re.findall(r'\d+', (Popen(["qsub", scriptName], stdout=PIPE).communicate()[0].strip()))[0] #call qsub and catch the job identifier
        jobDone = "%s.JOBDONE" %(jobID)
        
        st = "_".join(options.splicetypes)
        if options.outfile is  None:
            outfile = os.path.join(options.prefix, (options.bam.replace(".bam", ".splices.") + st + ".pickle"))
        else:
            outfile = options.outfile
        out = os.path.join(options.prefix, outfile)
        print "Job %s submitted.  Output will be here: %s" %(jobID ,out)
        if options.wait is True:
            import spin
            spin.spin(jobDone)
            
                
    
