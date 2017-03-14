__author__ = 'gpratt'

import cPickle as pickle

#This is bad style, but I need to fix import errors with matplotlib and clip analysis display
from clipper.src.CLIP_analysis import *
from clipper.src import CLIP_analysis_display

from AS_Structure_tools import parse_AS_STRUCTURE_dict
from clipper.src.get_genomic_regions import GenomicFeatures
from optparse import OptionParser
import gffutils

import matplotlib as mpl
from matplotlib import rc
mpl.rcParams['svg.fonttype'] = 'none'
rc('text', usetex=False)
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})


def main(bedtool, bam, species, runPhast=False, motifs=[], k=[6], nrand=3,
         outdir=os.getcwd(), db=None, as_structure=None, genome_location=None, phastcons_location=None,
         regions_location=None, motif_location=os.getcwd(), metrics="CLIP_Analysis.metrics", extension="svg",
         infer=False):

    """

    Runs all analysies

    one thing to do is make graphs fail gracefully

    """
    print "starting"
    #gets clusters in a bed tools + names species
    clusters = os.path.basename(bedtool)
    species = species
    short_species = species.split("_")[0]

    out_dict = {}
    #In case names aren't unique make them all unique
    clusters_bed = pybedtools.BedTool(make_unique(pybedtools.BedTool(bedtool))).saveas()
    if len(clusters_bed) <= 1:
        print "No Clusters, killing now to save time"
        return

    coverage = get_bam_coverage(bam)
    counts = get_bam_counts(bam)

    #makes output file names
    clusters = str.replace(clusters, ".BED", "")
    k = [int(x) for x in k]

    #sets up output dirs
    make_dir(outdir)

    assigned_dir = os.path.join(outdir, "assigned")
    misc_dir = os.path.join(outdir, "misc")
    fasta_dir = os.path.join(outdir, "fasta")
    homerout_base = os.path.join(outdir, "homer")
    make_dir(homerout_base)
    homerout = os.path.join(homerout_base, clusters)

    make_dir(assigned_dir)
    make_dir(misc_dir)
    make_dir(fasta_dir)
    make_dir(homerout)

    assigned_regions, regions = regions_generator()

    if db is not None:
        db = gffutils.FeatureDB(db)
    else:
        print "gff utils db not defined, this is fine, but falling back onto pre-set region defentions"
        db = None
    genomic_features = GenomicFeatures(species, db,  regions_dir=regions_location)
    genomic_regions = genomic_features.get_genomic_regions()
    features = genomic_features.get_feature_locations()

    clusters_bed = pybedtools.BedTool(bedtool)
    if infer:
        clusters_bed = infer_info(clusters_bed, genomic_regions['genes'])

    clusters_bed = pybedtools.BedTool(make_unique(clusters_bed)).saveas()



    if as_structure is not None:
        genes_dict, genes_bed = parse_AS_STRUCTURE_dict(short_species, as_structure)
    else:
        print "AS STRUCTURE file not listed, alt-splicing figure will not be generated"

    cluster_regions = assign_to_regions(tool=clusters_bed, clusters=clusters,
                                        assigned_dir=assigned_dir, species=species, nrand=nrand)

    print "getting cluster sizes"
    region_sizes = get_sizes(cluster_regions)

    print "getting genomic regions sizes"
    genic_region_sizes = count_genomic_region_sizes(assigned_regions, species)
    print "counting reads per cluster"
    reads_in_clusters, reads_per_cluster = count_reads_per_cluster(cluster_regions['all']['real'], counts)

    print "counting total reads in regions"
    region_read_counts = count_total_reads(genomic_regions, counts)
    total_reads = region_read_counts['genes']

    print "clustering peaks"
    read_densities, classes = cluster_peaks(cluster_regions['all']['real'], coverage)

    #generates cluster lengths (figure 3)
    print "getting cluster lengths"
    cluster_lengths = bedlengths(cluster_regions['all']['real'])

    print "getting peak locations"
    features_transcript_closest, features_mrna_closest = get_feature_distances(cluster_regions['all']['real'], features,
                                                                               genomic_regions['exons'], )
    features_transcript_closest = {name:  pybedtools.BedTool(bedtool.saveas("%s_%s_transcript.bed" % (clusters, name)).fn) for name, bedtool in features_transcript_closest.items()}
    features_mrna_closest = {name:  pybedtools.BedTool(bedtool.saveas("%s_%s_mrna.bed" % (clusters, name)).fn) for name, bedtool in features_mrna_closest.items()}

    #Distribution counting only works if genes have been pre-assignned
    distributions = get_region_distributions(cluster_regions['all']['real'], genomic_regions)

    if as_structure is not None:
        #also builds figure 10 (exon distances)
        genomic_types = count_genomic_types(genes_dict)
        exon_types = get_closest_exon_types(cluster_regions['all']['real'], genes_dict)

        #generates figure 10 (exon distances)
        type_count = [exon_types["CE:"], exon_types["SE:"], exon_types["MXE:"], exon_types["A5E:"], exon_types["A3E:"]]

        genomic_type_count = [genomic_types["CE:"], genomic_types["SE:"], genomic_types["MXE:"], genomic_types["A5E:"],
                              genomic_types["A3E:"]]

        out_dict["genomic_type_count"] = genomic_type_count
        out_dict["type_count"] = type_count

    if genome_location is not None:
        make_fasta_files_from_regions(cluster_regions, clusters, fasta_dir, genome_location)
        calculate_homer_motifs(k, regions, clusters, fasta_dir, homerout)
        kmer_results = calculate_kmer_diff(k, regions, clusters, fasta_dir)
    else:
        print "No genome fasta file provide, motif identification will not be performed"

    phast_values = None
    print "starting phast"
    if runPhast:
        phast_values = calculate_phastcons(cluster_regions, phastcons_location)
    print "ending phast"

    motif_distances = []
    try:
        if motifs:
            motif_distances = generate_motif_distances(cluster_regions, region_sizes, motifs, motif_location,
                                                       short_species)
    except:
        pass

    out_dict['features_transcript_closest'] = features_transcript_closest
    out_dict['features_mrna_closest'] = features_mrna_closest
    out_dict['distributions'] = distributions
    out_dict["region_sizes"] = region_sizes
    out_dict["reads_in_clusters"] = reads_in_clusters
    out_dict["reads_out_clusters"] = (total_reads - reads_in_clusters)
    out_dict["cluster_lengths"] = cluster_lengths
    out_dict["reads_per_cluster"] = reads_per_cluster
    out_dict["genic_region_sizes"] = genic_region_sizes
    out_dict["kmer_results"] = kmer_results
    out_dict["motifs"] = motifs
    out_dict["phast_values"] = phast_values
    out_dict["motif_distances"] = motif_distances
    out_dict['data'] = np.array(read_densities)
    out_dict['classes'] = classes
    out_dict['region_read_counts'] = region_read_counts
    out_dict['homerout'] = homerout
    out_dict['regions'] = regions

    with open(os.path.join("%s.clip_analysis.pickle" % clusters), 'w') as out_file:
        pickle.dump(out_dict, file=out_file)
    print "file saved"

    visualize(clusters, extension, outdir)

    try:
        if motifs:
            motif_fig = CLIP_analysis_display.plot_motifs(motif_distances)
            motif_fig.savefig(clusters + ".motif_distribution." + extension)
    except:
        pass

    with open(metrics, 'w') as outfile:
        outfile.write("FRiP\n")
        outfile.write("\t".join([str(float(reads_in_clusters) / float(total_reads))]) + "\n")

def visualize(clusters, extension, out_dir):
    with open(os.path.join("%s.clip_analysis.pickle" % clusters)) as out_file:
        clip_viz = CLIP_analysis_display.ClipVisualization(out_file)
    qc_fig = clip_viz.CLIP_QC_figure()
    distribution_fig = clip_viz.plot_distributions()
    qc_fig.savefig(os.path.join(out_dir, clusters + ".qc_fig." + extension))
    distribution_fig.savefig(os.path.join(out_dir, clusters + ".DistFig." + extension))

def call_main():
    parser = OptionParser()

    parser.add_option("--clusters", dest="clusters", help="BED file of clusters", metavar="BED")
    parser.add_option("--bam", dest="bam", help="The bam file from the CLIP analysis")
    parser.add_option("--species", "-s", dest="species", help = "genome version")
    parser.add_option("--runPhast", dest="runPhast", action="store_true", default=False, help="Run Phastcons ")
    parser.add_option("--motifs", dest="motifs", action="append", help="Motifs to use (files of motifs give must exist in motif_directory directory)", default=[])
    parser.add_option("--k", dest="k", action="append", help="k-mer and homer motif ananlysis", default=[6])
    parser.add_option("--nrand", dest="nrand", default=3, help="selects number of times to randomly sample genome", type="int")
    parser.add_option("--outdir", "-o", dest="outdir", default=os.getcwd(), help="directory for output, default:cwd")
    ##Below here are critical files that always need to be referenced
    parser.add_option("--gff_db", dest="db", help="gff database from gffutils to generate annotations with")
    parser.add_option("--AS_Structure", dest="as_structure",  help="Location of AS_Structure directory (chromosme files should be inside)", default=None)
    parser.add_option("--genome_location", dest="genome_location", help="location of all.fa file for genome of interest", default=None)
    parser.add_option("--phastcons_location", dest="phastcons_location",  help="location of phastcons file", default=None)
    parser.add_option("--regions_location", dest="regions_location",  help="directory of genomic regions for a species (default: clipper defined regions)", default=None)
    parser.add_option("--motif_directory", dest="motif_location",  help="directory of pre-computed motifs for analysis", default=os.getcwd())
    parser.add_option("--metrics", dest="metrics", default="CLIP_Analysis.metrics", help="file name to output metrics to")
    parser.add_option("--extension", dest="extension", default="svg", help="file extension to use (svg, png, pdf...)")
    parser.add_option("--infer", default=False, action="store_true", help="Infer peak center and gene if if peak finding algorithm doesn't report it")

    (options, args) = parser.parse_args()

    #error checking
    if options.clusters is None or options.bam is None or options.species is None:
        parser.print_help()
        exit()
    main(bedtool=options.clusters, bam=options.bam, species=options.species,
         runPhast=options.runPhast,
         motifs=options.motifs, k=options.k, nrand=options.nrand, outdir=options.outdir, db=options.db,
         as_structure=options.as_structure,
         genome_location=options.genome_location,
         phastcons_location=options.phastcons_location, regions_location=options.regions_location,
         motif_location=options.motif_location, metrics=options.metrics, extension=options.extension,
         infer=options.infer,
         )

if __name__ == "__main__":
    call_main()
