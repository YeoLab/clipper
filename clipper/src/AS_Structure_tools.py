'''
Created on Mar 19, 2013

@author: gabrielp
'''
import os
import numpy as np
import pybedtools
def parse_AS_STRUCTURE_dict(species, working_dir):
    
    """
    
    species - acceptable species to parse
    working_dir - location of as_structure files
    
    Important return values:
    
    info a dict of gene info
    genes - bedtool of gene locations
    
    Parses out all important AS structure - see constructed dict in function
    for information on what is needed...
    
    also returns bed file of genes 
    Should refactor to be a real object, but I'm lazy right now...
    
    """

    if species == "hg19" or species == "hg18":
        chroms = [str(x) for x in range(1, 23)] #1-22
        chroms.append("X")
        chroms.append("Y")
    elif species == "mm9":
        chroms = [str(x) for x in range(1, 20)] #1-19
        chroms.append("X")
        chroms.append("Y")
    elif species == "ce10":
        chroms = ["I", "II", "III", "IV", "V", "X", "M"]
    elif species == "test":
        chroms = ["1"]

    info = {}
    bed_string = ""
    
    for chrom in chroms:

        as_file = os.path.join(working_dir, species + ".tx." + chrom + ".AS.STRUCTURE")
        
        with open(as_file, "r") as AS_STRUCTURE_FILE:
            for line in AS_STRUCTURE_FILE:
                if not line.startswith(">"):
                    continue
                
                blank, gene, chrom, transcripts, d2, d3, d4, strand, number_of_exons, exonloc, intronloc, exonlen, intronlen, asType, locType = line.strip().split("\t")
                signstrand = "-"
                if int(strand) == 1:
                    signstrand = "+"
                info[gene] = {}
                info[gene]['chrom'] = "chr" + str(chrom)
                info[gene]['transcripts'] = transcripts
                info[gene]['strand'] = signstrand
                info[gene]['exons'] = {}
                info[gene]['introns'] = {}
                info[gene]['types'] = {}                        
                exons = exonloc.split("|")
                introns = intronloc.split("|")
                types = asType.split("|")
                info[gene]['numEx'] = int(number_of_exons)
                info[gene]['mRNA_length'] = 0
                info[gene]['premRNA_length'] = 0
                tx_start = np.Inf
                tx_stop  = np.NINF
                for i, exon in enumerate(exons):
                    if i == info[gene]['numEx']: #last exon is empty
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
                    if i == info[gene]['numEx']-1: #only number_of_exons-1 introns
                        continue
                    info[gene]['introns'][i] = intron
                    intstart, intstop = [int(x) for x in intron.split("-")]
                    info[gene]['premRNA_length'] += intstop-intstart+1
                info[gene]['tx_start'] = tx_start
                info[gene]['tx_stop'] = tx_stop
                bed_string += "\t".join([info[gene]['chrom'], str(info[gene]['tx_start']), str(info[gene]['tx_stop']), gene, "0", info[gene]['strand']]) + "\n"
            
            
    return info, pybedtools.BedTool(bed_string, from_string=True)

def parse_AS_STRUCTURE_COMPILED(species, working_dir):
    
    """
    
    species - acceptable species to parse
    working_dir - location of as_structure files
    
    Important return values:
    
    info a dict of gene info
    genes - bedtool of gene locations
    
    Parses out all important AS structure - see constructed dict in function
    for information on what is needed...
    
    also returns bed file of genes 
    Should refactor to be a real object, but I'm lazy right now...
    
    """

    if species == "hg19" or species == "hg18":
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
        
        with open(as_file, "r") as AS_STRUCTURE_FILE:
            for line in AS_STRUCTURE_FILE:
                if not line.startswith(">"):
                    continue
                
                blank, gene, chrom, transcripts, d2, d3, d4, strand, number_of_exons, exonloc, intronloc, exonlen, intronlen, asType, locType = line.strip().split("\t")
