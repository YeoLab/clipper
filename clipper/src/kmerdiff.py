'''
Created on Sep 19, 2012

@author: gabrielp

Clone of perl kmer diff program written by gene

#
# Output from 'compseq'
#
# The Expected frequencies are calculated on the (false) assumption that every
# word has equal frequency.
#
# The input sequences are:
#       NR_004433.3,923-1118_1
#       NR_033590.1,120-219_1

# ... et al.


Word size       4
Total count     13838

#
# Word  Obs Count       Obs Frequency   Exp Frequency   Obs/Exp Frequency
#
AAAA    4               0.0002891       0.0000051       56.2165053
AAAC    0               0.0000000       0.0000051       0.0000000
AAAD    0               0.0000000       0.0000051       0.0000000
AAAE    0               0.0000000       0.0000051       0.0000000
AAAF    0               0.0000000       0.0000051       0.0000000

Other    104        0.0075155    0.0000000    10000000000.0000000

'''
import os
import sys
import subprocess
def kmer_diff(file1, file2, k):
    
    """
    ASSUMES EMBOSS IS INSTALLED AND IN YOUR PATH
    (Need to include some checks for if this is not the case)
    
    file1 is a file to find words from
    file2 is another file to find words from
    k is the size of words to search for
    
    """
    
    Motif = namedtype('Motif', ['freq1', 'freq2', 'delta'])
    
    #call compseq
    subprocess.call(['compseq', 
                     '-word', k, 
                     '-sequence', file1,
                     '-outfile', file1 + ".compseq"]
                    )
    
    subprocess.call(['compseq', 
                     '-word', k, 
                     '-sequence', file2,
                     '-outfile', file2 + ".compseq"]
                    )
    
    n1, freq1 = parse_compseq(file1 + ".compseq")
    n2, freq2 = parse_compseq(file2 + ".compseq")
    
    results = {}
    for key in freq1.keys():
        g = (freq1[key] + freq2[key]) / (n1 + n2)
        
        if g == 0:
            delta = 0
        else:
            delta = (freq1[key] / n1) - (freq2[key] / n2) / sqrt(((1/n1 + 1/n2) * g * (1-g)))
        
        results[key] = Motif(freq1[key], freq2[key], delta)
    
    return results, n1, n2
def parse_compseq(file):
    
    """
    
    parses compseq file to return the total number of something...
    
    also the number of occurances for each sequence
    
    """
    fi = open(file)
    
    result = {}
    for line in fi:
        line_split = line.strip().split()
        if line.startswith("Total count"):
            total = float(line_split[1])
        if line.startswith("Other"):
            other = float(line_split[1])
        if not (line.startswith("#") or line.startswith("Word size") or line.startswith("\n")):
            result[line[0]] = float(line[1])
    return total - other, result

if __name__ == "__main__":
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    k = sys.argv[3]
    results, n1, n2kmer_diff(file1, file2, k)
    for key in results.keys():
        result = results[key]
        print "%s\t%s\t%s\t%s" % (key.lower(), str(result.freq1), str(result.freq2), result.delta)
    print "%s: %s words\n %s: %s words\n" %(file1, n1, file2, n2)