'''
Created on Aug 2, 2012

@author: gabrielp
'''


#Things to include
#gene expression levels -> coverage for any given gene

#Might be worthwhile to define a genome data structure that can handle repetative elements 
import random
import numpy
from numpy import array, round, arange, cumsum, where, zeros
import pysam
import pybedtools


def assign_peaks(genome, peak_size, num_peaks):
    
    """

    Returns a list of locations for peaks
    randomly assigned to locations in the genome
    
    Input
    genome: genome (will decide on format later)
    peak_size: size of peak to generate
    
    TODO: This assignment does not properly model a true genome, fix to allow true modeling of human genome
    
    """
    
    peaks = set([])
    for x in range(num_peaks):
        start = random.randrange(0, len(genome))
        stop  = start + peak_size
        peaks.add((start, stop))
    
    return peaks



def distribute_background_weights(genome):
    
    """

    Distributes background weights across the entire genome, returns array of weights
    
    genome - array, or list of arrays to represent genome, still optimizing
    
    """
    
    return numpy.random.gamma(5, 10, len(genome))

def distribute_peak_weights(peaks, background_weights, enrichment_coeff):

    """
    
    Returns a genome with both peak weights and background weights calculated as an array
    
    peaks - list of tuples [(start, stop)] that represent the location of peaks in the genome
    background_weights - array that represents the background weight of any read mapping to the genome
    
    """
    
    #peaks should be an ordered list in this situation
    peaks = list(peaks)
    
    #calculate the average background weight
    average_background_weight = numpy.mean(background_weights)
    average_peak_weight = average_background_weight * enrichment_coeff 
    
    #gets the total size of all peaks and the total background weights 
    total_size = 0
    peak_weights = zeros(len(peaks))
    for i, locs in enumerate(peaks):
        start, stop = locs
        total_size += stop - start
        peak_weights[i] = sum(background_weights[start:stop])
    
    #figures the total weight left to assign to peaks
    total_peak_weight = int((total_size * average_peak_weight) - sum(peak_weights))
     
    
    
    #for our purposes weight will be discrete
    #randomly distributes all remaining weights in a powerlaw form to all peak objects
    for x in range(total_peak_weight):
        peak_probs = cumsum(peak_weights * 1. / sum(peak_weights))
        peak = where((peak_probs > random.random()) == True)[0][0]
        peak_weights[peak] += 1
        start, stop = peaks[peak]
        background_weights[random.randint(start, stop)] += 1

    #TODO distribute weights normally instead of unifromily 
    return background_weights
    
        
    
    
        
    #pick a peak based on its total binding weight
    #update weight at peak
    #iterate until weight is done being added
    
    #redistribute peaks as a bionimal distribution 
def distribute_reads(weights, total_reads):
    """
    
    Distributes reads along the weights array returns an array of number of reads at a specific start location
    
    weight - an array of weights to distribute reads along
    num_reads - total number of reads to distribute
    read_length - length of reads to distribute
    
    TODO: print out 
    """
    
    total_weight = sum(weights)
    start_sites = numpy.zeros(len(weights))
    
    #distribute read starts
    for read_location, w in enumerate(weights):
        fractional_weight = float(w) / float(total_weight)
        num_reads = fractional_weight * total_reads
        start_sites[read_location] = round(num_reads)
        
    return start_sites 

def output_bam(reads, read_length, outfile_name):
    
    """
    
    Prints a list of reads to standard out
    
    Inputs
    lst: a list of reads to output as a bam file
    
    """
    
    #human information, will need to abstract to other species eventually
    header = { 'HD': {'VN': '1.0'},
            'SQ': [ {'SN':'chr1',      'LN':249250621},
                    {'SN':'chr2',     'LN':243199373},
                    {'SN':'chr3',     'LN':198022430},
                    {'SN':'chr4',     'LN':191154276},
                    {'SN':'chr5',     'LN':180915260},
                    {'SN':'chr6',     'LN':171115067},
                    {'SN':'chr7',     'LN':159138663},
                    {'SN':'chr8',     'LN':146364022},
                    {'SN':'chr9',     'LN':141213431},
                    {'SN':'chr10',    'LN':135534747},
                    {'SN':'chr11',    'LN':135006516},
                    {'SN':'chr12',    'LN':133851895},
                    {'SN':'chr13',    'LN':115169878},
                    {'SN':'chr14',    'LN':107349540},
                    {'SN':'chr15',    'LN':102531392},
                    {'SN':'chr16',    'LN':90354753},
                    {'SN':'chr17',    'LN':81195210},
                    {'SN':'chr18',    'LN':78077248},
                    {'SN':'chr19',    'LN':59128983},
                    {'SN':'chr20',    'LN':63025520},
                    {'SN':'chr21',    'LN':48129895},
                    {'SN':'chr22',    'LN':51304566},
                    {'SN':'chrX',     'LN':155270560},
                    {'SN':'chrY',     'LN':59373566},
                    {'SN':'chrM',     'LN':16571},] } 
    

    outfile = pysam.Samfile(outfile_name, "wb", header = header)
    
    #prints out all reads in reads
    read_count = 0
    for read_location, num_read in enumerate(reads):
        for i in arange(num_read):
            a = pysam.AlignedRead()
            a.qname = "read_%i" % (read_count)
            a.seq = "A" * read_length
            a.flag = 16
            a.rname = 0 
            a.pos = read_location
            a.mapq = 255
            a.cigar = ( (0,read_length - 1), (1, 1) )
            a.mrnm = 0 
            a.isize = 0
            a.qual = "<" * read_length
            outfile.write(a)
            
            read_count += 1
    outfile.close()
    
    
def output_bed(peaks, outfile):
    
    #hacky way to make peaks into a bedfile, need to change from chr1 at some point
    bedstring = "\n".join(map(lambda x: "chr1\t" + x , map(lambda x: "\t".join(map(str, x)), peaks)))
    
    tool = pybedtools.BedTool(bedstring, from_string= True)
    tool.saveas(outfile
                )
def run():
    peaks = assign_peaks(range(60000), 50, 50)
    background_weights = distribute_background_weights(range(60000))
    total_weights = distribute_peak_weights(peaks, background_weights ,5)
    reads = distribute_reads(total_weights, 60000)
    output_bam(reads, 50, "foo.bam")
    output_bed(peaks, "foo.bed")
    #
    #pass

if __name__ == "__main__":
    run()