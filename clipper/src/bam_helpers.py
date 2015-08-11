__author__ = 'gpratt'

import HTSeq
from HTSeq import SAM_Alignment
import pysam

class Robust_BAM_Reader(HTSeq.BAM_Reader):

    def __iter__( self ):
        sf = pysam.Samfile(self.filename, "rb")
        self.record_no = 0
        for pa in sf:
            try:
                yield SAM_Alignment.from_pysam_AlignedRead( pa, sf )
            except OverflowError:
                pass
            self.record_no += 1

    def fetch( self, reference = None, start = None, end = None, region = None ):
        sf = pysam.Samfile(self.filename, "rb")
        self.record_no = 0
        try:
           for pa in sf.fetch( reference, start, end, region ):
            try:
                yield SAM_Alignment.from_pysam_AlignedRead( pa, sf )
            except OverflowError:
                pass
            self.record_no += 1

        except ValueError as e:
           if e.message == "fetch called on bamfile without index":
              print "Error: ", e.message
              print "Your bam index file is missing or wrongly named, convention is that file 'x.bam' has index file 'x.bam.bai'!"
           else:
              raise
        except:
           raise