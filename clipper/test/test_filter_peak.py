import unittest
from clipper.src.filter_peak import *
from clipper.src.call_peak import *

class Test(unittest.TestCase):

    def test_poissonP(self):
        """

        Performs unit testing on poissonP

        Will not test math behind calling poisson fuction more than nessessary as it has already been tested by scipy

        """

        # Case: fewer than 3 reads fall within a peak region: Expected result should return as though there are 3 expected reads
        result = poissonP(50, 3, 50, 2)
        self.assertAlmostEqual(result, (
                    1 - 0.64723188878223115))  # manually calculated on personal computer stats.poisson.cdf(3, 3)

        # Case: More than 3 reads fall within a peak region. Expected: Result should return as normal (+1 of what ever it should be)
        result = poissonP(50, 3, 50, 4)
        self.assertAlmostEqual(result, (1 - 0.26502591529736164))  # manually calculated stats.poisson.cdf(3, 5)

        # Case: poissonp throws an exception. Expected: Result should return 1
        result = poissonP(None, None, None, None)
        assert 1 == result

    def test_filter_peaks_dicts(self):
        """

        Tests filter results

        good for regression tests, need to do better verification...

        """

        peak1 = Peak(chrom="chr15",
                     genomic_start=1, genomic_stop=10,
                     gene_name="ENSG1",
                     strand="-",
                     thick_start=50, thick_stop=60,
                     peak_number=1, number_reads_in_peak=52, size=32, p=0,
                     effective_length=425, peak_length=32, area_reads=89, area_size=95, nreads_in_gene=400)
        peak2 = Peak(chrom="chr15",
                     genomic_start=200, genomic_stop=300,
                     gene_name="ENSG2",
                     strand="-",
                     thick_start=140, thick_stop=160,
                     peak_number=2, number_reads_in_peak=239, size=45, p=0,
                     effective_length=200, peak_length=45, area_reads=400, area_size=100, nreads_in_gene=300)

        peak_dict = {'loc': ['chr15', 'ENSG00000198901', 91509274, 91537804, '-'],
                     'Nclusters': 24,
                     'nreads': 2086,
                     'threshold': 32,
                     'clusters': [peak1, peak2]}

        transcriptome_size = 10000
        transcriptome_reads = 100000

        # try different params
        result = filter_peaks_dicts([peak_dict], .07, transcriptome_size, transcriptome_reads,
                                    use_global_cutoff=False, bonferroni_correct=True,
                                    superlocal=True, min_width=50, bypassfiltering=False)

        ans = 'chr15\t1\t10\tENSG1_1_52\t0.000199322700259\t-\t50\t60\nchr15\t200\t300\tENSG2_2_239\t2.24979875151e-58\t-\t140\t160\n'
        self.assertIn('ENSG1_1_52', result)
        self.assertIn('ENSG2_2_239', result)

        # lower poission cutoff
        result = filter_peaks_dicts([peak_dict], .00001, transcriptome_size, transcriptome_reads,
                                    use_global_cutoff=False, bonferroni_correct=True,
                                    superlocal=True, min_width=50, bypassfiltering=False)
        self.assertIn('ENSG2_2_239', result)
        self.assertNotIn('ENSG1_1_52', result)

        # use global cutoff
        result = filter_peaks_dicts([peak_dict], .007, transcriptome_size, transcriptome_reads,
                                    use_global_cutoff=True, bonferroni_correct=True,
                                    superlocal=False, min_width=50, bypassfiltering=False)
        self.assertIn('ENSG1_1_52', result)
        self.assertIn('ENSG2_2_239', result)

    def test_bonferroni_correct_filter_results(self):
        """

        Tests to make sure bonferroni correction works

        """

        transcriptome_size = 10000
        transcriptome_reads = 100000

        peak1 = Peak(chrom="chr15",
                     genomic_start=1, genomic_stop=10,
                     gene_name="ENSG1",
                     strand="-",
                     thick_start=50, thick_stop=60,
                     peak_number=1, number_reads_in_peak=52, size=32, p=0,
                     effective_length=425, peak_length=32, area_reads=89, area_size=95, nreads_in_gene=400)
        peak2 = Peak(chrom="chr15",
                     genomic_start=200, genomic_stop=300,
                     gene_name="ENSG2",
                     strand="-",
                     thick_start=140, thick_stop=160,
                     peak_number=2, number_reads_in_peak=239, size=45, p=0,
                     effective_length=200, peak_length=45, area_reads=400, area_size=100, nreads_in_gene=300)
        peak_dicts = [{'loc': ['chr15', 'ENSG00000198901', 91509274, 91537804, '-'],
                       'Nclusters': 24,
                       'nreads': 2086,
                       'threshold': 32,
                       'clusters': [peak1, peak2]}]

        result = filter_peaks_dicts(peak_dicts, .007, transcriptome_size, transcriptome_reads,
                                    use_global_cutoff=False, bonferroni_correct=True,
                                    superlocal=False, min_width=50, bypassfiltering=False)
        self.assertIn('ENSG1_1_52', result)
        self.assertIn('ENSG2_2_239', result)

    def test_count_transcriptome_reads(self):
        """

        Tests count_transcriptome_reads

        """

        results = [
            {'loc': ['chr15', 'ENSG00000198901', 91509274, 91537804, '-'],
             'Nclusters': 24,
             'nreads': 200,
             'threshold': 32,
             'clusters': {
                 'chr15\t1\t10\tENSG1\t3.44651351902e-09\t-\t50\t60': {'GeneP': 3.4465135190231422e-09, 'Nreads': 100,
                                                                       'SloP': 3.4465135190231422e-09, 'size': 32},
                 'chr15\t200\t300\tENSG2\t0.0\t-\t140\t160': {'GeneP': 0.0, 'Nreads': 100, 'SloP': 0.0, 'size': 45}}},
            {'loc': ['chr15', 'ENSG00000198901', 91509274, 91537804, '-'],
             'Nclusters': 24,
             'nreads': 200,
             'threshold': 32,
             'clusters': {
                 'chr15\t1\t10\tENSG1\t3.44651351902e-09\t-\t50\t60': {'GeneP': 3.4465135190231422e-09, 'Nreads': 100,
                                                                       'SloP': 3.4465135190231422e-09, 'size': 32},
                 'chr15\t200\t300\tENSG2\t0.0\t-\t140\t160': {'GeneP': 0.0, 'Nreads': 100, 'SloP': 0.0, 'size': 45}}},
        ]
        result = count_transcriptome_reads(results)

        self.assertEqual(400, result)

