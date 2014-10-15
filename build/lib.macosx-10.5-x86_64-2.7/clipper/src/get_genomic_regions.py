"""

class to get genomic features from gffutils _db

"""

import copy
import os

import gffutils
import pybedtools
from clipper import data_dir
class GenomicFeatures():
    """

    class to get genomic features from gffutils _db
    
    """
    def __init__(self, species, db, regions_dir=None):
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
        
        if species in ["hg19", "mm9"]:
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
                #they are poximal
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
    
    def _get_utrs(self, mrna, cds, feature_types):
        """
        
        mrna: gffutils feature an mRNA
        cds: list[gffutls features] coding sequences belonging to the mRNA
        features_types: list[str] features avaiable for parsing 
        
        returns back list[gffutils feature] for 5' and 3' UTRS
        """
        mrna_five_prime_utrs = []
        mrna_three_prime_utrs = []
        if (self._feature_names['five_prime_utr']  in feature_types and
            self._feature_names['three_prime_utr'] in feature_types):
            
            mrna_five_prime_utrs = (list(self._db.children(mrna, featuretype=self._feature_names['five_prime_utr'])))
            mrna_three_prime_utrs = (list(self._db.children(mrna, featuretype=self._feature_names['three_prime_utr'])))
        else:
            #If 5' and 3' utr annotations don't exist 
            #generate them from CDS and UTR information (handles gencode case)
            utrs = list(self._db.children(mrna, featuretype='UTR'))
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
        feature_types = list(self._db.featuretypes())
        for i, gene in enumerate(self._db.features_of_type('gene')):
            if i % 2000 == 0:
                print "processed %d genes" % (i)
                if i == 2000 and limit_genes:
                    break
                
            gene_list.append(gene)
            gene_five_prime_utrs = []
            gene_three_prime_utrs = []
            gene_cds = []
            gene_exons = []
            
            for mrna in self._db.children(gene, featuretype=self._feature_names['transcript']):
                mrna_cds = list(self._db.children(mrna, featuretype=self._feature_names['CDS']))
                mrna_exons = list(self._db.children(mrna, featuretype=self._feature_names['exon']))
                gene_exons += mrna_exons
                gene_cds += mrna_cds
                mrna_five_prime_utrs, mrna_three_prime_utrs = self._get_utrs(mrna, mrna_cds, feature_types)
                gene_five_prime_utrs += mrna_five_prime_utrs
                gene_three_prime_utrs += mrna_three_prime_utrs
                   
            gene_id = gene.attributes[self._feature_names['gene_id']]
            merged_exons = self._merge_and_rename_regions(gene_exons, gene_id)
            gene_introns = list(self._db.interfeatures(merged_exons))
            gene_prox_introns, gene_dist_introns = self._get_proximal_distal_introns(gene_introns, 
                                                                                     prox_size)
            
            exons += merged_exons
            introns += gene_introns
            prox_introns += self._rename_regions(gene_prox_introns, gene_id)
            dist_introns += self._rename_regions(gene_dist_introns, gene_id)
            cds += self._merge_and_rename_regions(gene_cds, gene.id)       
            five_prime_utrs += self._merge_and_rename_regions(gene_five_prime_utrs, gene_id)
            three_prime_utrs += self._merge_and_rename_regions(gene_three_prime_utrs, gene_id)

        #make exons and introns
        results = {"genes" : gene_list,
                   "five_prime_utrs" : five_prime_utrs,
                   "three_prime_utrs" : three_prime_utrs,
                   "cds" : cds, 
                   "proxintron" : prox_introns, 
                   "distintron": dist_introns,
                   "exons" : exons,
                   "introns" : introns}
        
        for name, intervals in results.items():
            intervals = pybedtools.BedTool(map(self._to_bed, intervals)).remove_invalid().sort().each(self._fix_chrom)
            
            if name in ["prox_introns", "dist_introns"]:
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
        
        five_prime_ends = []
        three_prime_ends = []
        poly_a_sites = []
        stop_codons = []
        start_codons = []
        transcription_start_sites = []
        for i, gene_id in enumerate(self._db.features_of_type('gene')):
            if i % 2000 == 0:
                print "processed %d genes" % (i)
                if i == 2000 and limit_genes:
                    break
                
            gene = { "five_prime_ends" : [],
                    "three_prime_ends" : [],
                    "poly_a_sites" : [],
                    "stop_codons" :  [],
                    "start_codons" :  [],
                    "transcription_start_sites" : []}
            try:
                for exon in self._db.children(gene_id, featuretype='exon'):
                    exon_start = copy.deepcopy(exon)
                    exon_start.start = exon.start + 1
   
                    exon_stop = copy.deepcopy(exon)
                    exon_stop.start = exon_stop.stop
                    exon_stop.stop = exon_stop.stop + 1  
                    
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
                    transcript_stop.stop = transcript_stop.stop + 1
                 
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
