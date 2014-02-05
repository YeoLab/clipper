import gffutils, pybedtools 

db = gffutils.FeatureDB('GENES.db')

id_to_gene = {}
mrna_to_id = {}



for x in db.features_of_type('gene'):
    id_to_gene[x.attributes['ID'][0]] = x.attributes['sequence_name'][0]

for ft in db.featuretypes():
    if ft == "gene":
        continue

    for x in db.features_of_type(ft):
        mrna_to_id[x.attributes['ID'][0]] = x.attributes['Parent'][0]




_ = pybedtools.BedTool([("chr" + i.chrom, i.start, i.stop, i.attributes['sequence_name'][0], '0', i.strand) for i in db.features_of_type('gene')]).sort().saveas('ce10_genes.bed')

_ = pybedtools.BedTool([("chr" + i.chrom, i.start, i.stop, id_to_gene[mrna_to_id[i['Parent'][0]]], '0', i.strand) for i in db.features_of_type('three_prime_UTR')]).sort().saveas('ce10_three_prime_utrs.bed')

_ = pybedtools.BedTool([("chr" + i.chrom, i.start, i.stop, id_to_gene[mrna_to_id[i['Parent'][0]]], '0', i.strand) for i in db.features_of_type('five_prime_UTR')]).sort().saveas('ce10_five_prime_utrs.bed')

_ = pybedtools.BedTool([("chr" + i.chrom, i.start, i.stop, id_to_gene[mrna_to_id[i['Parent'][0]]], '0', i.strand) for i in db.features_of_type('CDS')]).sort().saveas('ce10_cds.bed')

exons = pybedtools.BedTool([("chr" + i.chrom, i.start, i.stop, id_to_gene[mrna_to_id[i['Parent'][0]]], '0', i.strand) for i in db.features_of_type('exon')]).sort().saveas('ce10_exons.bed')

from clipper.src.CLIP_analysis import get_introns
introns  = get_introns(exons).sort().saveas("ce10_introns.bed")


    

