import unittest
from clipper.src.utils import *
import clipper
import os

class Test(unittest.TestCase):

    def test_check_for_index(self):
        """

        Performs unit tests on check_for_index function

        """

        # Test if string is null, expected result is operation
        # throws file not exist exception
        handle = None
        self.assertRaises(TypeError, check_for_index, handle)

        # Test if bam file doesn't exist, expected result is operation throws
        # file does not exist exception
        handle = "/foo/bar"
        self.assertRaises(NameError, check_for_index, handle)

        # Test if file is not bam, but exists expected
        # result is to throw improper file error
        handle = clipper.test_file("test_peakfinder.py")
        self.assertRaises(NameError, check_for_index, handle)

        # Test if file is bam and indexed expected
        # result is returns 1 and succedes
        # should also check if file exists, but I'm lazy
        handle = clipper.test_file("indexed_test.bam")
        result = check_for_index(handle)
        assert result == None

        # Test if file is bam and not indexed, expected
        # result is returns one and succedes
        # should also check if file exists, but I'm lazy
        handle = clipper.test_file("not_indexed_test.bam")
        result = check_for_index(handle)
        assert result == None

        # cleanup (should be in taredown)
        os.remove(clipper.test_file("not_indexed_test.bam.bai"))
    # def test_build_transcript_data_gtf(self):
    #     """
    #      THIS FUNCTION IS DEPRECATED AND UNUSED
    #     Tests build transcript data gtf, tests two genes, with some noise
    #
    #     """
    #
    #     # tests pre-mrna
    #     genes = build_transcript_data_gtf(pybedtools.BedTool(clipper.test_file("data.gtf")), True).sort()
    #     true_genes = pybedtools.BedTool(
    #         [["chrI", "AS_STRUCTURE", "mRNA", 7741935, 7950951, "0", "+", ".",
    #           "gene_id=NR_070240; transcript_id=NR_070240; effective_length=209016"],
    #          ["chrI", "AS_STRUCTURE", "mRNA", 8378298, 8378421, "0", "-", ".",
    #           "gene_id=NM_001129046; transcript_id=NM_001129046; effective_length=123"], ]
    #     ).sort()
    #
    #     self.assertEqual(str(genes), str(true_genes))
    #
    #     # tests mrna lengths
    #     genes = build_transcript_data_gtf(pybedtools.BedTool(clipper.test_file("data.gtf")), False).sort()
    #     true_genes = pybedtools.BedTool(
    #         [["chrI", "AS_STRUCTURE", "mRNA", 7741935, 7950951, "0", "+", ".",
    #           "gene_id=NR_070240; transcript_id=NR_070240; effective_length=30"],
    #          ["chrI", "AS_STRUCTURE", "mRNA", 8378298, 8378421, "0", "-", ".",
    #           "gene_id=NM_001129046; transcript_id=NM_001129046; effective_length=123"], ]
    #     ).sort()
    #
    #     self.assertEqual(str(genes), str(true_genes))
    #
    # def test_build_transcript_data_gtf_longer(self):
    #     """
    #       THIS FUNCTION IS DEPRECATED AND UNUSED
    #     Tests build transcript data, but makes sure it gets the longer of two transcripts with the same gene name
    #
    #     """
    #
    #     genes = build_transcript_data_gtf(pybedtools.BedTool(clipper.test_file("data_2.gtf")), False).sort()
    #     true_genes = pybedtools.BedTool(
    #         [["chrI", "AS_STRUCTURE", "mRNA", 7741935, 7950970, "0", "+", ".",
    #           "gene_id=NR_070240; transcript_id=NR_070241; effective_length=41"]],
    #     ).sort()
    #
    #     self.assertEqual(str(genes), str(true_genes))

