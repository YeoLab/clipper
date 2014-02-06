import gffutils, pybedtools 

#usage: python get_genomic_regions.py > ce10_genes.gff3

db = gffutils.FeatureDB('GENES.db')


id_to_gene = {}
mrna_to_id = {}

for x in db.features_of_type('gene'):
    id_to_gene[x.attributes['ID'][0]] = x.attributes['sequence_name'][0]
    mRNAs = db.children(x, featuretype='mRNA')
    premRNA_len = len(x)
    lns = [0,]
    for mRNA in mRNAs:
        exons = db.children(mRNA, featuretype='exon')
        this_len = sum([len(i) for i in exons])
        lns.append(this_len)
    mRNA_len = max(lns)
    if mRNA_len == 0:
        mRNA_len = premRNA_len
    x.chrom = x.chrom.replace('tDNA', '')
    print str(x).replace('sequence_name', 'gene_id') +  ';mrna_length=%d;premrna_length=%d' % (mRNA_len, premRNA_len)


    

for ft in db.featuretypes():
    if ft == "gene":
        continue

    for x in db.features_of_type(ft):
        mrna_to_id[x.attributes['ID'][0]] = x.attributes['Parent'][0]

_ = pybedtools.BedTool([("chr" + i.chrom.replace("tDNA", ""), i.start, i.stop, i.attributes['sequence_name'][0], '0', i.strand) for i in db.features_of_type('gene')]).sort().saveas('ce10_genes.bed')

_ = pybedtools.BedTool([("chr" + i.chrom.replace("tDNA", ""), i.start, i.stop, id_to_gene[mrna_to_id[i['Parent'][0]]], '0', i.strand) for i in db.features_of_type('three_prime_UTR')]).sort().saveas('ce10_three_prime_utrs.bed')

_ = pybedtools.BedTool([("chr" + i.chrom.replace("tDNA", ""), i.start, i.stop, id_to_gene[mrna_to_id[i['Parent'][0]]], '0', i.strand) for i in db.features_of_type('five_prime_UTR')]).sort().saveas('ce10_five_prime_utrs.bed')

_ = pybedtools.BedTool([("chr" + i.chrom.replace("tDNA", ""), i.start, i.stop, id_to_gene[mrna_to_id[i['Parent'][0]]], '0', i.strand) for i in db.features_of_type('CDS')]).sort().saveas('ce10_cds.bed')

def rename(x):
    x.name = ",".join(set(x.name.split(",")))
    return x

exons = pybedtools.BedTool([("chr" + i.chrom.replace("tDNA", ""), i.start, i.stop, id_to_gene[mrna_to_id[i['Parent'][0]]], '0', i.strand) for i in db.features_of_type('exon')]).sort().merge(nms=True, s=True).each(rename).saveas('ce10_exons.bed')

from clipper.src.CLIP_analysis import get_introns
introns  = get_introns(exons).sort().saveas("ce10_introns.bed")

dist = []
prox = []
dist_len = 500

for intron in introns:
    if len(intron) < dist_len * 2:
        prox.append(intron)

    else:
        up = intron
        down = intron
        middle = intron
        up.stop = up.start + dist_len
        down.start = down.stop - dist_len
        middle.start = up.stop
        middle.stop = down.start
        prox.append(up)
        prox.append(down)
        dist.append(middle)

pybedtools.BedTool(prox).sort().saveas("ce10_proxintron500.bed")
pybedtools.BedTool(dist).sort().saveas("ce10_distintron500.bed")

from clipper.src.CLIP_analysis import get_feature_locations_ce10
get_feature_locations_ce10('.', 'ce10', db)
