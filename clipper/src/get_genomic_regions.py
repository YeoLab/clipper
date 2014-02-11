"""

class to get genomic features from gffutils _db

"""

import copy
import os

import gffutils
import pybedtools

class GenomicFeatures():
    """

    class to get genomic features from gffutils _db
    
    """
    def __init__(self, regions_dir, species, db):
        """
        
        creates genomic features function, chooses 
        how to direct creation of features based on _species
        
        regions_dir : str location to create region
        species: str species (hg19, mm9, ce10
        db: gffutils FeatureDb object
        
        """
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
                             }
            self._fix_chrom = self._fix_chrom_null
            
        elif species in ["ce10"]:
            self._feature_names = {
                             "five_prime_utr" : "five_prime_UTR",
                             "three_prime_utr" : "three_prime_UTR",
                             "exon" : "exon",
                             "CDS" : "CDS",
                             "gene_id" : "ID",
                             }
            self._fix_chrom = self._fix_chrom_ce10
    
    def _to_bed(self, interval):
        """
        interval: gffutils interval
        
        converts gffutils to bed format, setting id to be the 
        id of the gene
        
        """
        gene_id = interval.attributes[self._feature_names['gene_id']][0]
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
            regions = sorted(list(self._db.merge(regions)), 
                                           key = lambda x: x.start)
        
            return self._rename_regions(regions, gene_id)
        return []
    
    def get_genomic_regions(self, prox_size=500):
        
        """
        
        returns bedtool of all non-overlapping regions in the genome, exons, cds, 3' utrs and 5' utrs
        _species - string of the _species to analyze
        _db - _db handle generated by gtf utils
        
        Potental off by one bug here, need to examine more closely
        """
        region_and_species = os.path.join(self._regions_dir, self._species)
        regions = ["genes",
                   "five_prime_utrs"
                   "three_prime_utrs"
                   "cds",
                   "exons",
                   "introns",
                   "prox_introns",
                   "dist_introns",
                   ]
        try:
            results = {}
            for k, intervals in results.items():
                if k in ["prox_introns", "dist_introns"]:
                    results[k] = pybedtools.BedTool("%s_%s%d.bed" % (region_and_species, 
                                                                     k, prox_size))
                else:
                    results[k] = pybedtools.BedTool("%s_%s.bed" % (region_and_species, 
                                                                   prox_size))
    
        except Exception as e:
            pass
        
        genes = self._db.features_of_type('gene')
        three_prime_utrs = []
        five_prime_utrs = []
        cds = []
        exons = []
        dist_introns = []
        prox_introns = []
        gene_list = []
        introns = []
        for i, gene in enumerate(genes):
            if i % 500 == 0:
                print "processed %d genes" % (i)
                #break
                
            mrnas = list(self._db.children(gene, featuretype='mRNA'))
            #print "gene"
            gene_three_prime_utrs = []
            gene_five_prime_utrs = []
            gene_cds = []
            gene_exons = []
            
            gene_list.append(gene)
            for mrna in mrnas:
                cur_cds = list(self._db.children(mrna, featuretype=self._feature_names['CDS']))
                cur_exons = list(self._db.children(mrna, featuretype=self._feature_names['exon']))
                gene_exons += cur_exons
                gene_cds += cur_cds
                
                try:
                    gene_five_prime_utrs += (list(self._db.children(mrna, featuretype=self._feature_names['five_prime_utr'])))
                    gene_three_prime_utrs += (list(self._db.children(mrna, featuretype=self._feature_names['three_prime_utr'])))
                except:
                    #If 5' and 3' utr annotations don't exist generate them from CDS and UTR information (handles gencode case)
                    utrs =  list(self._db.children(mrna, featuretype='UTR'))
                    if len(cur_cds) == 0:
                        continue
                        
                    first_cds, last_cds = cur_cds[0], cur_cds[-1]
                                        
                    #this is ugly code, can I fix it?
                    for utr in utrs:
                        if utr.strand == "+":
                            if utr.stop < first_cds.start:
                                utr.featuretype = self._feature_names["five_prime_utr"]
                                gene_five_prime_utrs.append(utr)
                            elif last_cds.stop < utr.start:
                                utr.featuretype = self._feature_names["three_prime_utr"]
                                gene_three_prime_utrs.append(utr)
                            else:
                                print "something odd"
    
                        elif utr.strand == "-":
                            if last_cds.stop < utr.start:
                                utr.featuretype = self._feature_names["five_prime_utr"]
                                gene_five_prime_utrs.append(utr)
                            elif utr.start < first_cds.start:
                                utr.featuretype = self._feature_names["three_prime_utr"]
                                gene_three_prime_utrs.append(utr)
                            else:
                                print "odd in the negative strand"
                                
            gene_id = gene.attributes[self._feature_names['gene_id']]
            merged_exons = self._merge_and_rename_regions(gene_exons, gene.id)
            exons += merged_exons
            
            gene_introns = list(self._db.interfeatures(merged_exons))
            introns += gene_introns
            cur_prox_introns, cur_dist_introns = self._get_proximal_distal_introns(gene_introns, 
                                                                                  prox_size)
            prox_introns += self._rename_regions(cur_prox_introns, gene_id)
            dist_introns += self._rename_regions(cur_dist_introns, gene_id)
            
            cds += self._merge_and_rename_regions(gene_cds, gene.id)       
            five_prime_utrs += self._merge_and_rename_regions(gene_five_prime_utrs, gene_id)
            three_prime_utrs += self._merge_and_rename_regions(gene_three_prime_utrs, gene_id)

    
        #make exons and introns
        exons = pybedtools.BedTool(map(self._to_bed, exons)).sort()
        results = {"genes" : gene_list,
                   "five_prime_utrs" : five_prime_utrs,
                   "three_prime_utrs" : three_prime_utrs,
                   "cds" : cds, 
                   "prox_introns" : prox_introns, 
                   "dist_introns": dist_introns,
                   "exons" : exons,
                   "introns" : introns}
        
        for name, intervals in results.items():
            intervals = pybedtools.BedTool(map(self._to_bed, introns)).sort().each(self._fix_chrom)
            
            if name in ["prox_introns", "dist_introns"]:
                results[name] = intervals.saveas(region_and_species + "_%s%d.bed" % (name, 
                                                                                     prox_size))
            else:
                results[name] = intervals.saveas(region_and_species + "_%s.bed" % (name))

        return results
        
    def get_feature_locations(self):
        
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
        regions = ["five_prime_ends",
                   "three_prime_ends",
                   "poly_a_sites",
                   "stop_codons",
                   "start_codons",
                   "transcription_start_sites"]
        try:
            region_and_species = os.path.join(self._regions_dir, self._species)
            return {region : pybedtools.BedTool("%s_%s.bed" % (region_and_species, 
                                                               region)) for region in regions}
    
        except:
            pass
        
        genes = self._db.features_of_type('gene')
        five_prime_ends = []
        three_prime_ends = []
        poly_a_sites = []
        stop_codons = []
        start_codons = []
        transcription_start_sites = []
        for gene in genes:
            try:
                for exon in self._db.children(gene, featuretype='exon'):
                    exon_start = [exon.chrom, exon.start, exon.start + 1, 
                                  gene.id, "0", gene.strand]
                    exon_stop = [exon.chrom, exon.stop, exon.stop, 
                                 gene.id, "0", gene.strand]
                    
                    if exon.strand == "-":
                        exon_start, exon_stop = exon_stop, exon_start 
                        
                    five_prime_ends.append(exon_start)
                    three_prime_ends.append(exon_stop)
                
                #transcript vs mRNA need to look at the difference
                for transcript in self._db.children(gene, featuretype='mRNA'):
                    transcript_start = [transcript.chrom, transcript.start, 
                                        transcript.start + 1, gene.id, 
                                        "0", gene.strand]
                    transcript_stop = [transcript.chrom, transcript.stop, 
                                       transcript.stop + 1, gene.id, 
                                       "0", gene.strand]
                    
                    if transcript.strand == "-":
                        transcript_start, transcript_stop = transcript_stop, transcript_start
                        
                    poly_a_sites.append(transcript_stop)
                    transcription_start_sites.append(transcript_start)
                if self._species == "ce10": #need to generalize later
                    for transcript in self._db.children(gene, featuretype='mRNA'):
                        try:
                            cds = list(self._db.children(transcript, featuretype='CDS'))
                            first, last = cds[0], cds[-1]
                            
                            if first.strand == '-':
                                first, last = last, first
                            start_codons.append([first.chrom, 
                                                     first.start, 
                                                     first.start + 1, 
                                                     gene.id, 
                                                     "0", 
                                                     first.strand])
                            stop_codons.append([last.chrom, 
                                                    last.stop, 
                                                    last.stop + 1, 
                                                    gene.id, 
                                                    "0", 
                                                    last.strand])

                        except:
                            pass
                else: #for hg19 and mm9 gencode 
                    for start_codon in self._db.children(gene, featuretype='start_codon'):
                        start_codons.append([start_codon.chrom, 
                                             start_codon.stop, 
                                             start_codon.stop + 1, 
                                             gene.id, 
                                             "0", 
                                             gene.strand])
                        
                    for stop_codon in self._db.children(gene, featuretype='stop_codon'):
                        stop_codons.append([stop_codon.chrom, 
                                            stop_codon.stop, 
                                            stop_codon.stop + 1, 
                                            gene.id, 
                                            "0", 
                                            gene.strand])
                    
            except IndexError:
                pass
            
        
        results = { "five_prime_ends" : pybedtools.BedTool(five_prime_ends),
                 "three_prime_ends" :pybedtools.BedTool(three_prime_ends),
                 "poly_a_sites" :  pybedtools.BedTool(poly_a_sites),
                 "stop_codons" :  pybedtools.BedTool(stop_codons),
                 "start_codons" :  pybedtools.BedTool(start_codons),
                 "transcription_start_sites" : pybedtools.BedTool(transcription_start_sites)}

        for name, intervals in results.items():
            results[name] = intervals.each(self._fix_chrom).saveas("%s_%s.bed" % (region_and_species, name))

        return results   