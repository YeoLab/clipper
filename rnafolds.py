from optparse import OptionParser, SUPPRESS_HELP
from subprocess import Popen, PIPE, call
import pybedtools
from deap import dtm
import pickle
from Bio import SeqIO
import os
import sys


def get_folds_probabilities_single(FastaIO, foldparams="-u 1 -W 80 -L 40"):


    """
    runs RNAplfold and returns a dictionary with the probability that each base is unpaired
    """

    name, seq = FastaIO.id, FastaIO.seq.tostring()
    
    process = Popen(("/nas3/yeolab/Software/ViennaRNA-1.8.4/Progs/RNAplfold %s" %(foldparams)).split(), stdin=PIPE, stdout=PIPE)
    #print ("/nas3/yeolab/Software/ViennaRNA-1.8.4/Progs/RNAplfold %s" %(foldparams)).split()    
    #print ">%s\n%s\n" %(name, seq)
    run_items = process.communicate(">%s\n%s\n" %(name, seq))
    all_analyzed = run_items[0].replace(">", "").split()
    probabilities = {}
    for file in all_analyzed:
        filename = file + "_lunp"
        print filename
        try:
            f = open(filename, 'r')
            probabilities[file] = list()
            for line in f:
                if "#" in line:
                    continue #skip header
                l = line.strip().split("\t")
                probability=l[1]
                probabilities[file].append(probability)
            f.close()
            os.remove(file + "_lunp")
            os.remove(file + "_dp.ps")
        except:
            print "I had problems with %s, I'll skip it" %(filename)
    return probabilities



def main(options):
    fafile = options.fasta

    outpick = fafile.replace(".fa", ".pickle")
    
    fastafile = SeqIO.parse(open(fafile, 'r'), 'fasta')
    r = dtm.map(get_folds_probabilities_single, fastafile, foldparams=(options.fp).replace("_", " "))
            
    f = open(outpick, 'w')
    pickle.dump(r, file=f)


if __name__== "__main__":

    usage="python %s <options> --fasta <fastafile>\n"
    description="Run RNAplfold parallel-ly.  expects a file with a \".fa\" extension and outputs a <fastafile>.pickle file which is a pickled version of the fold output"
        
    parser = OptionParser(usage=usage, description=description)        
    parser.add_option("--fasta", dest="fasta", default=None, help="fasta file")    
    parser.add_option("--job_name", dest="job_name", default="fold", help="name for submitted job.  default:%default", metavar="NAME")
    parser.add_option("--processors", dest="np", default=4, help="number of processors to use.  default:%default", type="int", metavar="NP")
    parser.add_option("--notify", dest="notify", default=None, help="email address to notify of start, errors and completion", metavar="EMAIL")
    parser.add_option("--foldparams", dest="fp", type="str", action="store", default="-W_80_-L_40_-u_1", help="folding parameters default=%default", metavar="PARAMETER STRING")
    parser.add_option("--start", dest="start", default=False, action="store_true", help=SUPPRESS_HELP) #private, don't use

    (options, args) = parser.parse_args()
    options.prefix = os.getcwd()

    #main(options)
    #exit()


    
    if options.start is True:
        dtm.start(main, options)
        #main(options)
    
    else:
        scriptName = os.path.join(options.prefix, options.job_name+".RUNFOLD.sh")
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
        shScript.close()
        call(["qsub", scriptName]) #subprocess.call
            
                







