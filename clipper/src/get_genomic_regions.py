"""

class to get genomic features from gffutils _db

"""

import copy
from collections import defaultdict
import os


import gffutils
import pybedtools
from clipper import data_dir


class GenomicFeatures():
    """

    class to get genomic features from gffutils _db
    
    """
    def __init__(self, species, db=None, regions_dir=None, gencode=False):
        """
        
        creates genomic features function, chooses 
        how to direct creation of features based on _species
        
        regions_dir : str location to create region
        species: str species (hg19, mm9, ce10
        db: gffutils FeatureDb object
        
        """

        if regions_dir == None:
            regions_dir = os.path.join(data_dir(), "regions")
        self._regions_dir = regions_dir
        self._db = db
        self._species = species

        #I'm going to be lazy and leave this here, its needed to make a new genomic features for human genomes
        #engineering so it doesn't take too much time on load will be slightly annoying so just uncomment when you need it
        result = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        for x, feature in enumerate(db.all_features()):
            gene_ids = feature.attributes['gene_id']
            transcript_ids = feature.attributes['transcript_id']
            feature_type = feature.featuretype


            if feature_type == "gene":
                if len(gene_ids) != 1:
                    print gene_ids[0]
                    break

                result[gene_ids[0]]['gene'] = feature
            else:
                for gene_id in gene_ids:
                    for transcript_id in transcript_ids:
                        result[gene_id][transcript_id][feature_type].append(feature)

        self._feature_hash = result

        if species in ["hg19", "mm9", "hg19_v19", "GRCh38_v24", "hb27"] or gencode:
            self._feature_names = {
                             "five_prime_utr" : "five_prime_utr",
                             "three_prime_utr" : "three_prime_utr",
                             "exon" : "exon",
                             "CDS" : "CDS",
                             "gene_id" : "gene_id",
                             "transcript" : "transcript",
                             }
            self._fix_chrom = self._fix_chrom_null
            
        elif species in ["ce10"]:
            self._feature_names = {
                             "five_prime_utr" : "five_prime_UTR",
                             "three_prime_utr" : "three_prime_UTR",
                             "exon" : "exon",
                             "CDS" : "CDS",
                             "gene_id" : "ID",
                             "transcript" : "mRNA",
                             }
            self._fix_chrom = self._fix_chrom_ce10
        if self._db is not None:
            self.featuretypes = list(self._db.featuretypes())
        else:
            self.featuretypes = None
    def _to_bed(self, interval):
        """
        interval: gffutils interval
        
        converts gffutils to bed format, setting id to be the 
        id of the gene
        
        """
        gene_id = interval.attributes[self._feature_names['gene_id']]
        if type(gene_id) is list:
            gene_id = gene_id[0]        
        return (interval.chrom, interval.start, 
                interval.stop, gene_id, "0", interval.strand)
    
    def _fix_chrom_ce10(self, interval):
        """
        
        Adjusts ce10 interval to be in line with ucsc browser format
        
        """
        interval.chrom = "chr" + interval.chrom.replace("tDNA", "")
        return interval 
    
    def _fix_chrom_null(self, interval):
        
        """
        
        null fixer for when things are in ucsc format
        
        """
        return interval
    
    def _get_proximal_distal_introns(self, introns, prox_size=500):
        """
    
        From a given gene returns all its proximal and distal introns, 
        proximal and distal being defined by prox_size 
        gene = iterator of exons belonging to a single gene
        prox_size = int size of proximal introns
    
        returns bedtool(proximal introns), bedtool(distal introns)
    
        """
    
        prox_introns = []
        dist_introns = []

        for intron in introns:
            #want distal introns to have size at least one, otherwise
            if len(intron) <= (prox_size * 2) + 3: 
                #they are proximal
                prox_introns.append(intron)
            else:
                #create prox and dist intron ranges from intron 
                #(this is dangerous, but copying doesn't work
                start_prox_intron = copy.deepcopy(intron)
                start_prox_intron.stop = intron.start + prox_size
                prox_introns.append(start_prox_intron)
            
                stop_prox_intron = copy.deepcopy(intron)
                stop_prox_intron.start = intron.stop - prox_size
                prox_introns.append(stop_prox_intron)
            
                dist_intron = copy.deepcopy(intron)
                dist_intron.start = intron.start + prox_size + 1
                dist_intron.stop = intron.stop - prox_size - 1
                dist_introns.append(dist_intron)
        return prox_introns, dist_introns

    def _rename_regions(self, intervals, gene_id):
        """
        
        renames list of intervals
        
        intervals: gffutils interval list
        gene_id: value to replace id with
        
        returns updated list of intervals
        """
        updated_regions = []
        for interval in intervals:
            interval.attributes[self._feature_names['gene_id']] = gene_id
            updated_regions.append(interval)
        return updated_regions
    
    def _merge_and_rename_regions(self, regions, gene_id):
        """
        
        region: list gffutils features to merge 
        gene_id: str string to rename gene_id as (changes feature names to be sane)
        
        """
        
        if len(regions) > 0:
            regions = self._db.merge(sorted(list(regions), 
                                            key = lambda x: x.start)) 
                                           
        
            return self._rename_regions(regions, gene_id)
        return []
    
    def _get_utrs(self, gene, mrna, cds, feature_types):
        """
        
        mrna: gffutils feature an mRNA
        cds: list[gffutls features] coding sequences belonging to the mRNA
        features_types: list[str] features avaiable for parsing 
        
        returns back list[gffutils feature] for 5' and 3' UTRS
        """
        mrna_five_prime_utrs = []
        mrna_three_prime_utrs = []
        if (self._feature_names['five_prime_utr'] in feature_types and self._feature_names['three_prime_utr'] in feature_types):
            
            mrna_five_prime_utrs = (list(self._db.children(mrna, featuretype=self._feature_names['five_prime_utr'])))
            mrna_three_prime_utrs = (list(self._db.children(mrna, featuretype=self._feature_names['three_prime_utr'])))
        else:
            #If 5' and 3' utr annotations don't exist 
            #generate them from CDS and UTR information (handles gencode case)
            utrs = list(self._feature_hash[gene][mrna]['UTR'])
            if len(cds) == 0:
                return [], []
                
            first_cds, last_cds = cds[0], cds[-1]
            for utr in utrs:
                if utr.strand == "+":
                    if utr.stop < first_cds.start:
                        utr.featuretype = self._feature_names["five_prime_utr"]
                        mrna_five_prime_utrs.append(utr)
                    elif last_cds.stop < utr.start:
                        utr.featuretype = self._feature_names["three_prime_utr"]
                        mrna_three_prime_utrs.append(utr)
                    else:
                        print "something odd"
    
                elif utr.strand == "-":
                    if last_cds.stop < utr.start:
                        utr.featuretype = self._feature_names["five_prime_utr"]
                        mrna_five_prime_utrs.append(utr)
                    elif utr.start < first_cds.start:
                        utr.featuretype = self._feature_names["three_prime_utr"]
                        mrna_three_prime_utrs.append(utr)
                    else:
                        print "odd in the negative strand"
        return mrna_five_prime_utrs, mrna_three_prime_utrs

    def _interval_key(self, interval):
        return interval.start

    def _gene_regions(self, gene_id, prox_size=500):
        gene_five_prime_utrs = []
        gene_three_prime_utrs = []
        gene_cds = []
        gene_exons = []
        gene_introns = []
        gene_dist_introns = []
        gene_prox_introns = []
        for mrna in self._feature_hash[gene_id].keys():
            if mrna == "gene":
                continue
            mrna_cds = []
            if self._feature_names['CDS'] in self._feature_hash[gene_id][mrna]:
                mrna_cds = list(self._feature_hash[gene_id][mrna][self._feature_names['CDS']])
            mrna_exons = []
            if self._feature_names['exon'] in self._feature_hash[gene_id][mrna]:
                mrna_exons = list(self._feature_hash[gene_id][mrna][self._feature_names['exon']])
            gene_exons += mrna_exons
            gene_cds += mrna_cds
            mrna_five_prime_utrs, mrna_three_prime_utrs = self._get_utrs(gene_id, mrna, mrna_cds, self.featuretypes)
            gene_five_prime_utrs += mrna_five_prime_utrs
            gene_three_prime_utrs += mrna_three_prime_utrs
            mrna_introns = list(self._db.interfeatures(self._db.merge(sorted(mrna_exons, key=self._interval_key))))
            gene_introns += mrna_introns
            mrna_prox_introns, mrna_dist_introns = self._get_proximal_distal_introns(mrna_introns, prox_size)
            gene_dist_introns += mrna_dist_introns
            gene_prox_introns += mrna_prox_introns

        gene = self._feature_hash[gene_id]['gene']
        gene_id = gene.attributes[self._feature_names['gene_id']]



        return self._merge_and_rename_regions(gene_cds, gene.id), \
               self._merge_and_rename_regions(gene_dist_introns, gene_id), \
               self._merge_and_rename_regions(gene_exons, gene_id),  \
               self._merge_and_rename_regions(gene_five_prime_utrs, gene_id), \
               self._merge_and_rename_regions(gene_introns, gene_id), \
               self._merge_and_rename_regions(gene_prox_introns, gene_id), \
               self._merge_and_rename_regions(gene_three_prime_utrs, gene_id)

    def get_genomic_regions(self, prox_size=500, limit_genes=False, flush_cashe=False):
        
        """
        
        returns bedtool of all non-overlapping regions in the genome, exons, cds, 3' utrs and 5' utrs
        _species - string of the _species to analyze
        _db - _db handle generated by gtf utils
        
        Potental off by one bug here, need to examine more closely
        """
        region_and_species = os.path.join(self._regions_dir, self._species)
        regions = ["genes", "five_prime_utrs", "three_prime_utrs", "cds", 
                   "exons", "introns", "proxintron", "distintron",
                   ]
        try:
            if flush_cashe:
                raise ValueError
            results = {}
            for region in regions:
                if region in ["proxintron", "distintron"]:
                    results[region] = pybedtools.BedTool("%s_%s%d.bed" % (region_and_species, 
                                                                           region, prox_size))
                else:
                    results[region] = pybedtools.BedTool("%s_%s.bed" % (region_and_species, 
                                                                         region))
            return results
        except ValueError as e:
            print e
            pass
        
        three_prime_utrs = []
        five_prime_utrs = []
        cds = []
        exons = []
        dist_introns = []
        prox_introns = []
        gene_list = []
        introns = []
        for i, gene in enumerate(self._feature_hash.keys()):
            gene_list.append(self._feature_hash[gene]['gene'])
            if i % 2000 == 0:
                print "processed %d genes" % (i)
                if i == 2000 and limit_genes:
                    break


            gene_cds, gene_dist_introns, gene_exons, gene_five_prime_utrs, gene_introns, gene_prox_introns, gene_three_prime_utrs = self._gene_regions(gene)
            three_prime_utrs += gene_three_prime_utrs
            five_prime_utrs += gene_five_prime_utrs
            cds += gene_cds
            exons += gene_exons
            dist_introns += gene_dist_introns
            prox_introns += gene_prox_introns
            introns += gene_introns

        #make exons and introns
        results = {"genes": gene_list,
                   "five_prime_utrs": five_prime_utrs,
                   "three_prime_utrs": three_prime_utrs,
                   "cds": cds,
                   "proxintron": prox_introns,
                   "distintron": dist_introns,
                   "exons": exons,
                   "introns": introns}
        
        for name, intervals in results.items():
            intervals = pybedtools.BedTool(map(self._to_bed, intervals)).remove_invalid().sort().each(self._fix_chrom)
            
            if name in ["proxintron", "distintron"]:
                results[name] = intervals.saveas(region_and_species + "_%s%d.bed" % (name, 
                                                                                     prox_size))
            else:
                results[name] = intervals.saveas(region_and_species + "_%s.bed" % (name))

        return results
        
    def get_feature_locations(self, limit_genes=False, flush_cashe=False):

        
        """
        
        Gets locations of genic features, five prime sites, 3 prime sites, poly a sites stop codons start codons and tss
        based off annotated gtf _db file
        
        _db - _db handle generated by gtf utils
        
        returns dict of bedfiles     { five_prime_ends : bedtool 
                                       three_prime_ends
                                       poly_a_sites
                                       stop_codons
                                       transcription_start_sites 
                                    } 
        
        """

        transcriptome = { "five_prime_ends" : [],
                    "three_prime_ends" : [],
                    "poly_a_sites" : [],
                    "stop_codons" :  [],
                    "start_codons" :  [],
                    "transcription_start_sites" : []}
        
        region_and_species = os.path.join(self._regions_dir, self._species)
        try:
            if flush_cashe:
                raise ValueError
         
            return {region : pybedtools.BedTool("%s_%s.bed" % (region_and_species, 
                                                               region)) for region in transcriptome}
    
        except ValueError:
            pass

        for i, gene_id in enumerate(self._db.features_of_type('gene')):
            if i % 2000 == 0:
                print "processed %d genes" % (i)
                if i == 2000 and limit_genes:
                    break
                
            gene = { "five_prime_ends": [],
                    "three_prime_ends": [],
                    "poly_a_sites": [],
                    "stop_codons":  [],
                    "start_codons":  [],
                    "transcription_start_sites": []}
            try:
                for exon in self._db.children(gene_id, featuretype='exon'):
                    exon_start = copy.deepcopy(exon)
                    exon_start.start = exon.start + 1
   
                    exon_stop = copy.deepcopy(exon)
                    exon_stop.start = exon_stop.stop
                    exon_stop.stop += 1
                    
                    if exon.strand == "-":
                        exon_start, exon_stop = exon_stop, exon_start 
                        
                    gene['five_prime_ends'].append(exon_start)
                    gene['three_prime_ends'].append(exon_stop)
                
                #transcript vs mRNA need to look at the difference
                for transcript in self._db.children(gene_id, featuretype=self._feature_names['transcript']):
                    transcript_start = copy.deepcopy(transcript)
                    transcript_start.stop = transcript.start + 1
                    
                    transcript_stop = copy.deepcopy(transcript)
                    transcript_stop.start = transcript_stop.stop
                    transcript_stop.stop += 1
                 
                    if transcript.strand == "-":
                        transcript_start, transcript_stop = transcript_stop, transcript_start
                        
                    gene['poly_a_sites'].append(transcript_stop)
                    gene['transcription_start_sites'].append(transcript_start)
                if self._species == "ce10": #need to generalize later
                    for transcript in self._db.children(gene_id, featuretype=self._feature_names['transcript']):
                        try:
                            cds = list(self._db.children(transcript, 
                                                         featuretype='CDS'))
                            
                            first_cds, last_cds = cds[0], cds[-1]

                            if first_cds.strand == '-':
                                first_cds, last_cds = last_cds, first_cds
                                
                            start_codon = first_cds
                            start_codon.stop = first_cds.start + 1
                            gene['start_codons'].append(start_codon)
                                
                            stop_codon = last_cds
                            stop_codon.start = stop_codon.stop
                            stop_codon.stop  = stop_codon.stop + 1
                            gene['stop_codons'].append(stop_codon)

                        except:
                            pass
                else: #for hg19 and mm9 gencode 
                    for start_codon in self._db.children(gene_id, featuretype='start_codon'):
                        start_codon.stop = start_codon.start + 1
                        gene['start_codons'].append(start_codon)
                        
                    for stop_codon in self._db.children(gene_id, featuretype='stop_codon'):
                        stop_codon.start = stop_codon.stop
                        stop_codon.stop  = stop_codon.stop + 1
                        gene['stop_codons'].append(stop_codon)
                    
            except IndexError:
                pass
            gene_id = gene_id.attributes[self._feature_names['gene_id']]
            for region in gene:
                transcriptome[region] += self._merge_and_rename_regions(gene[region], gene_id)

        for name, intervals in transcriptome.items():
            transcriptome[name] = pybedtools.BedTool(map(self._to_bed, intervals)).\
                remove_invalid().sort().each(self._fix_chrom).saveas("%s_%s.bed" % (region_and_species, name))

        return transcriptome
