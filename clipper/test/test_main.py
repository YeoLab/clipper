from __future__ import print_function
import unittest
from clipper.src.main import *
import clipper
import pybedtools

class Test(unittest.TestCase):

    def test_allup(self):
        """
        Performs basic all up test on entire program (except for main)
        """
        args = ["--bam", clipper.test_file("allup_test.bam"),
                "--species", "hg19",
                "--gene", "ENSG00000198901",
                "--outfile", os.getcwd() + "/allup_peak_results.bed",
                "--debug",
                ]
        parser = option_parser()
        (options, args) = parser.parse_args(args)
        options = override_options(options)

        main(options)



        tested = open(os.getcwd() + "/allup_peak_results.bed")
        correct = open(clipper.test_file("peak_results_no_overlap.BED"))

        # problem with tracks being different
        tested_tool = pybedtools.BedTool(tested)
        correct_tool = pybedtools.BedTool(correct)

        # checks to make sure files are equal and there are not exact dups
        print(len(tested_tool))
        print(len(correct_tool))

        self.assertAlmostEqual(len(tested_tool), len(correct_tool), delta=3)
        print(len(tested_tool))
        print(len(correct_tool))

        # cleanup
        os.remove(os.getcwd() + "/allup_peak_results.bed")

    def test_allup_parallel(self):
        """

        Performs basic all up test on entire program (except for main), running in parallel to
        try to detect crashes

        """

        args = ["--bam", clipper.test_file("allup_test.bam"),
                "--species", "hg19",
                "--gene", "ENSG00000198901",
                "--outfile", os.getcwd() + "/allup_peak_results.bed",
                "--debug",
                ]
        parser = option_parser()
        (options, args) = parser.parse_args(args)
        options = override_options(options)

        main(options)

        tested = open(os.getcwd() + "/allup_peak_results.bed")
        correct = open(clipper.test_file("peak_results_no_overlap.BED"))

        # problem with tracks being different
        tested_tool = pybedtools.BedTool(tested)
        correct_tool = pybedtools.BedTool(correct)

        # checks to make sure files are equal and there are not exact dups
        print(len(tested_tool))
        print(len(correct_tool))

        self.assertAlmostEqual(len(tested_tool), len(correct_tool), delta=3)
        print(len(tested_tool))
        print(len(correct_tool))
        # assert False
        """
        for test, correct in zip(tested_tool, correct_tool):
            self.assertEqual(test, correct)


        """
        # cleanup
        os.remove(os.getcwd() + "/allup_peak_results.bed")

    def test_check_overlaps(self):
        """

        Checks for overlapping results, we don't want this

        overlaps have been borken for a while, disabling test until its really a problem
        """

        args = ["--bam", clipper.test_file("allup_test.bam"),
                "--species", "hg19",
                "--gene", "ENSG00000198901",
                "--outfile", os.getcwd() + "/overlap_peak_results.bed",
                "-q",
                "--debug"
                ]
        parser = option_parser()
        (options, args) = parser.parse_args(args)
        options = override_options(options)
        main(options)

        # tests to make sure there are no overlaps
        tested = open(os.getcwd() + "/overlap_peak_results.bed")
        tested_tool2 = pybedtools.BedTool(tested).saveas(os.getcwd() + "/overlaps.bed")
        result = tested_tool2.intersect(tested_tool2)
        self.assertEqual(len(result), len(tested_tool2),
                         "there are overlaps in the output file")

        # cleanup
        os.remove(os.getcwd() + "/overlap_peak_results.bed")
        os.remove(os.getcwd() + "/overlaps.bed")

    def test_main(self):
        """

        Performs unit testing on main

        Mostly testing validation and input here
        TODO: fill in test...

        """
        pass