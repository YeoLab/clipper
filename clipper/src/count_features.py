__author__ = 'gpratt'
from clipper.src import CLIP_analysis
import pybedtools
import pandas as pd
from collections import OrderedDict
import os
import clipper
import HTSeq
import collections
from itertools import izip

#NEED TO TEST AND CHECK THIS CODE BEFORE SUBMITTING THE PAPER
class RegionCounter():
    def __init__(self, species):

        regions = OrderedDict()
        regions["all"] = "All"
        regions["cds"] = "CDS"
        regions["three_prime_utrs"] = "3' UTR"
        regions["five_prime_utrs"] = "5' UTR"
        regions["proxintron500"] = "Proximal\nIntron"
        regions["distintron500"] = "Distal\nIntron"
        regions['exons'] = "Exons"

        assigned_regions = regions.copy()
        del assigned_regions['all']
        self.species = species
        self.assigned_regions = assigned_regions
        self.features = self.make_features()

    def make_features(self):
        Region = collections.namedtuple("Region", ["region", "gene_id"])

        bedtracks = {}
        for region in self.assigned_regions:
            bedtracks[region] = pybedtools.BedTool(os.path.join(clipper.data_dir(),
                                                                "regions", "%s_%s.bed" % (self.species, region)))



        features = HTSeq.GenomicArrayOfSets("auto", stranded=True)
        for region, bedtrack in bedtracks.items():
            for iv, interval in izip(CLIP_analysis.bed_to_genomic_interval(bedtrack), bedtrack):
                features[iv] = set([Region(region, interval.name)])
        return features


    def count_features(self, bam_file):
        bam_file = HTSeq.BAM_Reader(bam_file)
        counts = collections.defaultdict(collections.Counter)
        for x, read in enumerate(bam_file):
            region_ids = set()
            for iv, val in self.features[read.iv].steps():
                region_ids |= val

            gene_ids = {region_id.gene_id for region_id in region_ids}
            if len(gene_ids) == 1:

                cur_regions = {region_id.region for region_id in region_ids}

                for region in self.assigned_regions:
                    if region in cur_regions:
                        break

                gene_id = list(region_ids)[0]

                counts[gene_id.gene_id][region] += 1

            elif len(gene_ids) == 0:
                counts["_no_feature"]['none'] += 1
            else:
                counts["_ambiguous"]['none'] += 1

        return pd.DataFrame(counts)