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
from collections import namedtuple
from math import sqrt
import pandas as pd


class Motif(namedtuple('Motif', ['freq1', 'freq2', 'delta'])):
    
    """
    
    Class to encapsulate motifs 
    
    """


def run_jellyfish(fn, out_file, k):
    with open(os.devnull, 'w') as fnull:
        db_file = "{}.{}.jellyfish".format(fn, k)
        jellyfish_db = ['jellyfish', 'count',
                        '-m', str(k),
                        '-s', '50M',
                        fn,
                        '-o', db_file]
        jellyfish_dump = ['jellyfish', 'dump',
                          '-c', '-t',
                          '-o', out_file,
                          db_file
                            ]
        try:
            subprocess.call(jellyfish_db,
                                shell=False,
                                stdout=fnull,
                                stderr=fnull
                            )

            subprocess.call(jellyfish_dump,
                                shell=False,
                                stdout=fnull,
                                stderr=fnull)

        except OSError as e:
            print "jellyfish probably not installed, please install"
            raise e
    df = pd.read_table(out_file,
                       index_col=0,
                       header=None,
                       names=['kmer', 'count'])
    return df


def kmer_diff(file1, file2, k):

    """

    ASSUMES EMBOSS IS INSTALLED AND IN YOUR PATH
    (Need to include some checks for if this is not the case)

    file1 is a file to find words from
    file2 is another file to find words from
    k is the size of words to search for

    """



    if not os.path.exists(file1):
        raise IOError(file1 + " does not exist")
    if not os.path.exists(file2):
        raise IOError(file2 + " does not exist")

    freq1 = run_jellyfish(file1,
             file1 + ".counts.txt",
             k)
    n1 = float(freq1.sum()['count'])

    freq2 = run_jellyfish(file2,
              file2 + ".random.counts.txt",
             k)
    n2 = float(freq2.sum()['count'])

    freqs = pd.merge(freq1, freq2, left_index=True, right_index=True, how="outer", suffixes=["_freq1", "_freq2"])
    freqs = freqs.fillna(0)

    results = {}
    for key, row in freqs.iterrows():

        try:
            g = (row.count_freq1 + row.count_freq2) / (n1 + n2)
        except ZeroDivisionError:
            g = 0

        #delta is some sort of strange change measurement, not quite
        if g == 0:
            delta = 0
        else:
            delta = ((row.count_freq1 / n1) - (row.count_freq2 / n2)) / sqrt((1/n1 + 1/n2) * g * (1-g))

        results[key] = Motif(row.count_freq1, row.count_freq2, delta)

    return results, n1, n2

if __name__ == "__main__":
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    k = sys.argv[3]
    results, n1, n2 = kmer_diff(file1, file2, k)
    for key in results.keys():
        result = results[key]
        print "%s\t%s\t%s\t%s" % (key.lower(), str(result.freq1), str(result.freq2), result.delta)
    print "%s: %s words\n %s: %s words\n" %(file1, n1, file2, n2)
