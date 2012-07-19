from deap import dtm
from subprocess import Popen, call, PIPE
import re
import pybedtools
import os,sys
import pickle
from seqTools import fetchseq, chop
#from CLIP_analysis import get_phastcons
host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
if "optiputer" in host or "compute" in host:
    base = "/nas/nas0"
elif "tcc" in host or "triton" in host:
    base = "/projects"    

try:
    pybedtools.set_tempdir(base + "/scratch/lovci/pybedtools_tmp")
except:
    try:
        pybedtools.set_tempdir(base + "/lovci/projects/tmp/python_tmp")
    except:
        print "No valid pybedtools tmp directory"
        exit()


RNAhybridcmd = base + "/yeolab/Software/RNAhybrid/bin/RNAhybrid"
RNAhybrid_parser = base + "/lovci/gscripts/parseRNAhybrid.pl"

def RNAhybrid_hits(query, target, mfe_cutoff=-40, chi=-20.0, theta=10.0):
    """
    query and target are sequences, other options are parameters for RNAhybrid, see those docs
    output is a parsed list of matches. Parsing Scheme: upper-case = match; lower-case = non-pair; dash=gap
    """
    
    #print " ".join(map(str, [RNAhybridcmd, "-d", (str(chi) + "," + str(theta)), "-n", "70", "-m", "10000", "-e", str(mfe_cutoff), str(target), str(query)]))
    tries = 0
    while tries < 5: #recover from problems opening a Popen instance
        try:
            rh = Popen([RNAhybridcmd, "-d", (str(chi) + "," + str(theta)),
                        "-n", "70", "-m", "10000", "-e", str(mfe_cutoff),
                        str(target), str(query)], stdout = PIPE)
            parse = Popen(["perl", RNAhybrid_parser], stdin=rh.stdout, stdout=PIPE)
            rh.stdout.close()
            output = parse.communicate()
            results = output[0].split("\n")[:-1]
            if tries > 0:
                print "problem fixed"
            tries=100 #made it through, exit the loop
        except:
            print "problem with %s" %(str(query))
            tries = tries+1

    if tries == 5:
        print "I had a problem opening a Popen instance"
        print "error, %s" %(sys.exc_info())
    for i, r in enumerate(results):
        results[i] = r.split("\t")
    return results

def build_bed(filename, outfile = "test.bed",    proxLen = 500):
    """
    Make a bed-file of up/down-stream distal/proximal intronic regions with labels in the 'name' field
    """
    try:
        f = open(filename, 'r')
    except:
        print "couldn't open %s" %(filename)
        exit()
    data = list()
    out = open(outfile, 'w')
    for line in f.readlines():
        id, type, wholeLoc, exonLoc  = line.strip().split("\t")
        shortLine = wholeLoc
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
        exon = "\t".join([chrom , exstart, exstop, "%".join(["exon", shortLine, id, type]), "1", strand])
        out.write(exon +"\n")
        if uI_length < proxLen:
            #print "upstream intron too short, skipping upstream distal class for %s" %(line)
            uI = "\t".join([chrom , uI_start, uI_stop, "%".join(["ui", shortLine, id, type]), "1", strand])
            uIProx = uI.replace("ui", "uiProx")
            out.write("\n".join([uI,  uIProx]) + "\n")           
        else:
            uI = "\t".join([chrom , uI_start, uI_stop, "%".join(["ui", shortLine, id, type]), "1", strand])
            if strand == "+":
                uIProx = "\t".join([chrom , str(int(uI_stop)-proxLen), uI_stop, "%".join(["uiProx", shortLine, id, type]), "1", strand])
                uIDist = "\t".join([chrom , uI_start, str((int(uI_stop)-proxLen)), "%".join(["uiDist", shortLine, id, type]), "1", strand])                
            else:
                uIProx = "\t".join([chrom , uI_start, str(int(uI_start)+proxLen), "%".join(["uiProx", shortLine, id, type]), "1", strand])
                uIDist = "\t".join([chrom , str(int(uI_start)+proxLen),  uI_stop, "%".join(["uiDist", shortLine, id, type]), "1", strand])
            out.write("\n".join([uI, uIDist, uIProx]) + "\n")
        if dI_length < proxLen:
            #print "downstream intron too short, skipping downstream distal class for %s" %(line)
            dI = "\t".join([chrom , dI_start, dI_stop, "%".join(["di", shortLine, id, type]), "1", strand])
            dIProx = dI.replace("di", "diProx")
            out.write("\n".join([dI, dIProx]) + "\n")
        else:
            dI = "\t".join([chrom , dI_start, dI_stop, "%".join(["di", shortLine, id, type]), "1", strand])
            if strand == "+":
                dIDist = "\t".join([chrom , str(int(dI_start)+proxLen), dI_stop, "%".join(["diDist", shortLine, id, type]), "1", strand])
                dIProx = "\t".join([chrom , dI_start, str(int(dI_start) + proxLen), "%".join(["diProx", shortLine, id, type]), "1", strand])
            else:
                dIProx = "\t".join([chrom , str(int(dI_stop) - proxLen), dI_stop , "%".join(["diProx", shortLine, id, type]), "1", strand])                
                dIDist = "\t".join([chrom , dI_start, str(int(dI_stop) - proxLen), "%".join(["diDist", shortLine, id, type]), "1", strand])
            out.write("\n".join([dI, dIDist, dIProx]) + "\n")            
    out.close()
    tool  =pybedtools.BedTool(outfile)
    tool = pybedtools.BedTool(str(tool.remove_invalid()), from_string=True).saveas(outfile)
    return(tool)

def bed12FromMatchedPair(a,b, x,y, chr="chrZ", strand = "+", color="255,0,0", score="1", name="bedline"):
    """
    Make separate bed6 entries a single 'linked' bed12 entry
    input is the x/y coordinates of the two regions to be linked
    output is a single bed12 line
    """
    a,b,x,y = map(int, [a,b,x,y])
    if a<=y and b >=x:
        #overlapping...skip
        return
    if b < x:
        start = a
        stop = y
        block1len = b-a
        block2len = y-x
        block1Start = 0
        block2Start = x-a
    else:
        start = x
        stop = b
        block1len = y-x
        block2len = b-a
        block1Start = 0
        block2Start = a-x
    bedline = "\t".join(map(str, [chr, start, stop, name, score, strand,  start, stop, color, "2", ",".join(map(str, [block1len, block2len])), ",".join(map(str,[block1Start, block2Start]))]))
    return bedline
    
def matched_positions(query, target, querychunks=70, queryOverlap=.6, targetchunks=10000, targetOverlap=.1, mfe_cutoff=-40, species="hg19", quiet =False):
    """
    Run RNAhybrid using a genome range for query and target, step
    along each in (query/target)chunks-sized overlapping windows, target and query are
    bedtools intervals returned positions are 0-based from the start of
    target
    """
    if quiet is not True:
        print "Pairing %s with %s" %(str(query), str(target))

    qSeq = fetchseq(species, query.chrom.replace("chr", ""),
                    query.start, query.stop, query.strand)
    tSeq = fetchseq(species, target.chrom.replace("chr", ""),
                    target.start, target.stop, target.strand)

    tChop = chop(tSeq, chunkSize=targetchunks, chunkOverlap=targetOverlap)
    matches = list()


    for i, T in tChop:
        qChop = chop(qSeq, chunkSize=querychunks, chunkOverlap=queryOverlap)
        for j, Q in qChop:
            #print str(i) + "\t" + str(j)
            #print T + "\t" + Q
            hits = RNAhybrid_hits(Q, T, mfe_cutoff=mfe_cutoff)
            for hit in hits:
                tname, tlen, tpos, qname, qlen, mfe, pval, tseq, qseq = hit
                tMatchLen = len(tseq.replace("-", ""))                
                if target.strand == "+":
                    Tgenome_start = target.start +  i+int(tpos) + 1
                    Tgenome_stop = Tgenome_start + tMatchLen + 1
                    Qgenome_start = query.start + j
                    Qgenome_stop = Qgenome_start + len(Q)
                else:
                    Tgenome_stop = target.stop - i - int(tpos) + 1
                    Tgenome_start = Tgenome_stop - tMatchLen + 1
                    Qgenome_stop = query.stop - j
                    Qgenome_start = Qgenome_stop - len(Q)
                    
                if Qgenome_start< Tgenome_start:
                    color= "255,0,0"
                else:
                    color="0,0,255"

                bedline = bed12FromMatchedPair(Qgenome_start, Qgenome_stop,
                                               Tgenome_start, Tgenome_stop,
                                               chr=target.chrom, strand=target.strand,
                                               color=color, score=-float(mfe))
                #Db_entry = RNApair(bedline)
                
                if bedline is not None:
                    matches.append(bedline)
    if quiet is not True:
        print "Found %s structural matches" %(str(len(matches)))
    return matches


def run_fold(picklefile, dir=None, eventName=None, outfile=None, mfe_cutoff=-40, species = "hg19"):
    """
    given a picklefile (output from parse_event_detail), pair regions
    """
    if eventName ==None:#if a colloquial name is not given, get it from the filename
        eventNameF = picklefile.replace(".pickle", "").split("/")[-1] 
        eventName=eventNameF
    if not os.path.exists(picklefile):
        raise Exception
    info = pickle.load(open(picklefile, 'r'))
    #print "loaded file %s" %(picklefile)

    #if not "SE" in  info['type']:
        #print "skipping %s because it is of type: %s" %(picklefile, info['type'])
        #return

    matches =list()    
    try:
        di = info['di']
        ui = info['ui']
        diProx = info['diProx']
        uiProx = info['uiProx']
        diDist = info['diDist']
        uiDist = info['uiDist']
    except:
        #print "skipping %s because it is missing parts" %(picklefile)
        return matches
    #print "mfe_cutoff = %f" %(mfe_cutoff)
    try:
        matches.extend([i.replace("bedline", (eventName + "%" + "%s" %(info['type']) + "%" + "diProx_diDist")) for i in matched_positions(diProx, diDist, mfe_cutoff=mfe_cutoff, species=species)])
    except:
        print "error. Exception raised was: %s" %(str(sys.exc_info()))
        
    try:
        matches.extend([i.replace("bedline", (eventName + "%" + "%s" %(info['type']) + "%" + "diProx_uiDist")) for i in matched_positions(diProx, uiDist, mfe_cutoff=mfe_cutoff, species=species)])
    except:
        pass
    try:
        matches.extend([i.replace("bedline", (eventName + "%" + "%s" %(info['type']) + "%" + "diProx_uiProx")) for i in matched_positions(diProx, uiProx, mfe_cutoff=mfe_cutoff, species=species)])
    except:
        pass        
    try:
        matches.extend([i.replace("bedline", (eventName + "%" + "%s" %(info['type']) + "%" + "uiProx_uiDist")) for i in matched_positions(uiProx, uiDist, mfe_cutoff=mfe_cutoff, species=species)])
    except:        
        pass
    try:
        matches.extend([i.replace("bedline", (eventName + "%" + "%s" %(info['type']) + "%" + "uiProx_diDist")) for i in matched_positions(uiProx, diDist, mfe_cutoff=mfe_cutoff, species=species)])
    except:
        pass
    try:
        matches.extend([i.replace("bedline", (eventName + "%" + "%s" %(info['type']) + "%" + "uiProx_diProx")) for i in matched_positions(uiProx, diProx, mfe_cutoff=mfe_cutoff, species=species)])
    except:
        pass    
    all = "\n".join(matches)
    if outfile is not None:
        return pybedtools.BedTool(all, from_string=True).saveas(outfile, trackline="track_name='%s.RNAlinks' itemRgb=On" %(eventName))
    else:
        return all

      
def fold_a_dir(eventName, dir="structure_tmp/", outdir = "beds", rewrite=False, mfe_cutoff=-45, species="hg19"):

    #print "trying %s" %(eventName)
    from pybedtools import BedTool as BT
    import glob
    x = ""
    #print "Using output directory:%s" %(outdir)
    out = list()
    #if os.path.exists(os.path.join(outdir, (eventName + ".RNAlinks.bed"))) and not rewrite:
    #    #file exists and user has not specificed a rewrite
    #    return(BT(os.path.join(outdir, (eventName + ".RNAlinks.bed"))))



    #IMPORTANT:Delete locks later!!!

    #make a lock so other processes can't touch it. 
    if os.path.exists(os.path.join(outdir, (eventName + ".lock"))):
        print "skipping locked file"
        return None
    else:
        f = open(os.path.join(outdir, (eventName + ".lock")), 'w')
        f.write("lock")
        f.close()
                     
    
    if os.path.exists(os.path.join(outdir, (eventName + ".RNAlinks.bed"))): #this may lead to a race and data collision in the unlikely circumstance that two computers are running structure.py on the same set of genes to the same output directory
        #check whether the file contains "empty" indicating that the job was killed before a bed-file can be written or whether it shuold be re-tried
        
        f = open(os.path.join(outdir, (eventName + ".RNAlinks.bed")), 'r');
        if str(f.readline()).startswith("empty") or rewrite is True:
            #the file was opened previously, but never finished or, the user requested a re-write
            #print "removing existing files, probably corrupt intermediates"
            [os.remove(i) for i in glob.glob(os.path.join(outdir, (eventName + "*.RNAlinks.bed")))]#did not finish, partial files are probably corrupt
            pass
        else:
            f.close() #this is a legit bed file, i'll leave it be
            return(BT(os.path.join(outdir, (eventName + ".RNAlinks.bed"))))
        
        f.close() #this is a legit bed file, i'll leave it be        
    f = open(os.path.join(outdir, (eventName + ".RNAlinks.bed")), 'w');
    f.write("empty\n")
    f.close()   #reserve this event and go on.
    files =glob.glob(dir + "*%s*.pickle" %(eventName))
    if len (files) ==0:
        print "Nothing to be done for %s" %(eventName)
    if len(files) > 1:
        print "You have multiple files for %s, I'll merge the output" %(eventName)
        pass
    for file in files:
        print file
        try:
            outfile = os.path.join(outdir, "%s" %(file.split("/")[-1].replace(".pickle", ".RNAlinks.bed")))
            if os.path.exists(outfile) and not rewrite:
                out.append(BT(outfile))
                continue
            tool = run_fold(file, outfile=outfile,mfe_cutoff=mfe_cutoff, species=species)
            
            out.append(tool)
        except:
            print "problem with %s" %(file)
    if len(out) ==1:
        print "using a single event for %s" %(eventName)
        print "Wrote to file %s" %(outfile)
    else:
        import operator
        print "Merging multiple events"
        nCmps = len(out)
        for i in out:
            if i is None:
                nCmps = nCmps-1

        if nCmps > 0:
            try:
                newBed = ""
                for event in out:
                    for line in event:
                        if line is not None:
                            newBed = str(newBed) + "\n" + str(line)
                
                tool = BT(newBed, from_string=True).saveas(os.path.join(outdir, (eventName + ".RNAlinks.bed")), trackline="track_name='%s.RNAlinks' itemRgb=On" %(eventName)) #merge all bed lines into one and save with a track line
                print "Wrote to file %s" %(os.path.join(outdir, (eventName + ".RNAlinks.bed")))

            except:
                print sys.exc_info()
                tool = None
        else:
            tool = None
    return tool

def parse_event_detail(tool, tmpdir= "structure_tmp/", skip=True, findMe=None):
    """ makes .pickle files for each event """
    info = {}
    for line in tool:
        chr, start, stop, name, score, strand = str(line).strip().split("\t")
        if findMe is not None and findMe not in name: #for error checking, set findMe to some value to skip all lines not containing findMe
            continue
        relLoc, wholeLoc, eventName, type = name.split("%")
        eventNameF = eventName.replace("|", "_")
        if os.path.exists(os.path.join(tmpdir, (eventNameF + ".pickle"))):
            continue
        #if "SE" not in type:
        #continue
        try:
            info[eventName][relLoc] = line
        except:
            info[eventName] = {}
            info[eventName][relLoc] = line                
            info[eventName]["type"] = type
    for eventName in info:
        eventNameF = eventName.replace("|", "_")
        f = open(os.path.join(tmpdir, (eventNameF + ".pickle")), 'w')
        pickle.dump(info[eventName], file=f)
        f.close()

        return(info.keys())

def build_db(species, outfile="event_detail.BED", tmpdir="structure_tmp/"):
    try:
        regions = pybedtools.BedTool(outfile)
        regions = regions.remove_invalid()
    except:
        print "make a new db"
        regions = build_bed(base + "/yeolab/Genome/ensembl/AS_STRUCTURE/" + species + "data4/" + species + ".internal_exons.with_intron_flanks.table", outfile = outfile)        
    try:
        os.mkdir(tmpdir)
    except:
        print "WARNING: couldn't make %s... it may already exist" %(tmpdir)

    events = parse_event_detail(regions)
    return(events)
        
        
def get_names(db):
    tool = pybedtools.BedTool(db)
    names = list()
    for line in tool:
        chr, start, stop, name, score, strand = str(line).strip().split("\t")
        relLoc, wholeLoc, eventName, type = name.split("%")
        names.append(eventName)
    return names
        
def main(options):
    if options.gene is None:
        #options.gene =  x = list(map(str.strip, open(base + "/lovci/projects/FOX2/FOX2_human_brain/CLIP/analysis_gsnap/bound_genes.slop.p05.t0.intron_only.txt").readlines()))[:-1]
        options.gene =  x = list(map(str.strip, open(base + "/lovci/projects/conservation/hg19/mammal_cons/ultra_allIntron.genes.txt").readlines()))[:-1]
        genelist=options.gene
        #print "which genes to run?... this will take awhile"
        #genelist = get_names(options.db)
    else:
        genelist = options.gene
    #mongoport = "8585"
    #mongo = Popen(["ssh", "-L", ("%s:localhost:%s" %(mongoport, mongoport)), "oolite", "-N"])  #connect to mongo db


    if options.max_genes is not None:
        if not len(genelist) < options.max_genes: #already < max_genes genes in genelist
            import random
            #sample a random subset
            genelist = random.sample(genelist, options.max_genes)
    
    if options.serial is True:
        for gene in genelist:
            fold_a_dir(gene, rewrite=options.rewrite, mfe_cutoff=options.mfe_cutoff, dir=options.dbdir, species=options.species, outdir=options.outdir)
    else:
        dtm.map(fold_a_dir, genelist, rewrite=options.rewrite, mfe_cutoff=options.mfe_cutoff, dir=options.dbdir, species=options.species, outdir=options.outdir)


if __name__ == "__main__":
    print "base directory: %s" %(base)
    import sys
    import optparse
    from optparse import OptionParser, SUPPRESS_HELP
    usage = "None"
    description = "None"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--db",         dest="db",        default="event_detail.BED")
    parser.add_option("--rewrite",    dest="rewrite",   default=False,  action="store_true")    
    parser.add_option("--mfe_cutoff", dest="mfe_cutoff",default=-45)
    parser.add_option("--dbdir",      dest="dbdir",     default="structure_tmp/")
    parser.add_option("--outdir",      dest="outdir",     default="beds/")
    parser.add_option("--max_genes", dest="max_genes", default=None, type="int")
    parser.add_option("--gene",       dest="gene",      default=None,   action="append")
    
    parser.add_option("--SGEprefix",     dest="prefix",    default=os.getcwd(), help="location for sh script and error logs")
    parser.add_option("--start",      dest="start",     default=False,  action="store_true", help=SUPPRESS_HELP) #private, don't use
    parser.add_option("--species",    dest="species",   default="hg19")
    parser.add_option("--serial",     dest="serial",                    action="store_true", help="run genes in sequence (not parallel)")
    parser.add_option("--job_name",   dest="job_name",  default="STRUC", help="name for submitted job. Not used with --serial.  default:%default", metavar="NAME")
    parser.add_option("--processors", dest="np",        default=32, help="number of processors to use. Not used with --serial.  default:%default", type="int", metavar="NP")
    parser.add_option("--notify",     dest="notify",    default=None, help="email address to notify of start, errors and completion", metavar="EMAIL")
    parser.add_option("--wait_to_exit", dest="wait",    default=False, action ="store_true")    
    (options,args) = parser.parse_args()


    if options.serial is True:
        main(options)
        exit()
    if options.start is True:
        dtm.start(main, options)
    else:       
        if os.path.exists(options.db):
            print "Do you want to rebuild the database:%s of exons and introns?" %(options.db)
            userAnswer = sys.stdin.readline()
            if "y" in userAnswer:
                print "plase wait"
                build_db(options.species, outfile=options.db)
        else:
            print "Building database %s, press enter to continue" %(options.db)
            null = sys.stdin.readline()
            build_db(options.species, outfile=options.db)

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
            shScript.write("#!/bin/sh\n#PBS -N %s\n#PBS -o %s\n#PBS -e %s\n#PBS -V\n#PBS -S /bin/sh\n#PBS -l nodes=%d:ppn=8\n#PBS -q batch\n#PBS -l walltime=01:00:00\n" %(options.job_name, runout, runerr, nnodes))
                
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
        

        print "Job %s submitted." %(jobID)
        if options.wait is True:
            import spin
            spin.spin(jobDone)
  
