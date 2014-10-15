'''
Created on Sep 18, 2012

@author: gabrielp
'''

import cPickle as pickle
from collections import OrderedDict
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['interactive'] = False
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import random
import os
import pybedtools
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

class ClipVisualization():

    def __init__(self, pickle_file):
        """

        :param pickle:
        :return:
        """

        #TODO Remove this is a temporary fix, new versions of pickle file should contain this
        #OBJECT
        self.regions = OrderedDict()
        self.regions["all"] = "All"
        self.regions["cds"] = "CDS"
        self.regions["three_prime_utrs"] = "3' UTR"
        self.regions["five_prime_utrs"] = "5' UTR"
        self.regions["proxintron500"] = "Proximal\nIntron"
        self.regions["distintron500"] = "Distal\nIntron"

        data = pickle.load(pickle_file)

        self.reads_in_clusters = data['reads_in_clusters']
        self.reads_out_clusters = data['reads_out_clusters']
        self.cluster_lengths = data['cluster_lengths']
        self.reads_per_cluster = data['reads_per_cluster']
        self.premRNA = data['distributions']['genes']['total']
        self.mRNA = data['distributions']['exons']['total']
        self.exondist = data['distributions']['exons']['individual']
        self.introndist = data['distributions'] ['introns']['individual']
        self.genomic_locs = data['genic_region_sizes']
        self.clusters_locs = data['region_sizes']
        self.genomic_types = data['genomic_type_count']
        self.clusters_types = data['type_count']
        #self.homer_location = data['homerout']
        self.kmer_results = data['kmer_results']
        self.motifs = data['motifs']
        self.phastcons_values = data['phast_values']
        #self.regions = data['regions']
        self.read_densities = np.array(data['data'])
        self.classes = data['classes']
        self.distributions = data['distributions']

        #TODO Need to load in transcripts and mrnas for plotting

        #out_dict["premRNA_positions"] = premRNA_positions
        #out_dict["mRNA_positions"] = mRNA_positions
        #out_dict["exon_positions"] = exon_positions
        #out_dict["intron_positions"] = intron_positions

    def build_reads_in_clusters(self, ax, reads_in_clusters=None, reads_out_clusters=None):

        """

        Draws the number of reads in a cluster and the number of reads
        outside of clusters on a pie chart

        ax - the axis to draw on
        reads_in_clusters - int, the number of reads falling into a cluster
        reads_out_clusters - int the number of reads falling outside of a cluster

        """

        if reads_in_clusters is None:
            reads_in_clusters = self.reads_in_clusters
        if reads_out_clusters is None:
            reads_out_clusters = self.reads_out_clusters

        ax.pie([reads_in_clusters, reads_out_clusters],
               labels=["In Clusters", "Outside of\nClusters"])
        return ax

    def build_reads_per_cluster(self, ax_nreads, reads_per_cluster=None):

        """

        Draws the number of reads per cluster for each cluster

        ax - the axis to draw on
        reads_per_cluster - list, the number of reads in a cluster


        """
        if reads_per_cluster is None:
            reads_per_cluster = self.reads_per_cluster

        ax_nreads.hist(reads_per_cluster,
                       bins=50,
                       facecolor='#C8D2B0',
                       log=True,
                       range=(10, max(11, np.max(reads_per_cluster))))

        ax_nreads.set_xscale('log')
        ax_nreads.set_xlabel("log10(N reads)")
        ax_nreads.set_ylabel("Frequency")
        return ax_nreads

    def build_cluster_lengths(self, ax_lengths, cluster_lengths=None):

        """

        Selects a random sample of all cluster length
        and draws 2000 of them in a boxplot

        ax - the axis to draw on
        cluster_lengths - list, the length of each cluster


        """

        if cluster_lengths is None:
            cluster_lengths = self.cluster_lengths

        ax_lengths.set_yscale('log')
        ax_lengths.boxplot(random.sample(cluster_lengths,
                                         min(2000, len(cluster_lengths))),
                           vert=1)

        ax_lengths.set_xticklabels([])
        return ax_lengths

    def build_gene_distribution(self, ax, premrna=None, mrna=None):

        """

        Draws the distribution of clusters across mRNA and pre-mRNA meta-genes

        ax - the axis to draw on
        premRNA - a list of locations (in percent) that describes where a cluster binds
        in preMRNA

        mRNA - a list of locations (in percent) that describes where a cluster binds
        in mRNA

        """
        if premrna is None:
            premrna = self.premRNA
        if mrna is None:
            mrna = self.mRNA
        ax, alt_ax = self.build_distribution(ax, premrna, mrna)

        ax.set_xlabel("Fraction of region")
        ax.set_ylabel("pre-mRNA Location Frequency", color='red')
        alt_ax.set_ylabel("mRNA Location Frequency", color='blue')
        return ax, alt_ax

    def build_exon_exon_distribution(self, ax, exon=None, intron=None):

        """

        Draws the distribution of clusters across exon and intron meta-genes

        ax - the axis to draw on
        exon - a list of locations (in percent) that describes where a cluster binds
        in exons

        intron - a list of locations (in percent) that describes where a cluster binds
        in introns

        """
        if exon is None:
            exon = self.exondist
        if intron is None:
            intron = self.introndist
        ax, alt_ax = self.build_distribution(ax, exon, intron)

        ax.set_xlabel("Fraction of region")
        ax.set_ylabel("Exon Location Frequency", color='red')
        alt_ax.set_ylabel("Intron Location Frequency", color='blue')
        return ax, alt_ax

    def build_distribution(self, ax, dist1, dist2):

        """

        Plots two histograms on the same axis returning a second axis for labeling purposes.
        ax - the axis to draw on

        dist1 - the first distribution to plot
        dist2 - the second distribution to plot

        """
        #print dist1
        #print dist2

        alt_ax = ax.twinx()

        #error checking in case there is a null distribution for some reasion...
        if len(dist1) > 0:
            count, bins = np.histogram(dist1,
                                       bins=100,
                                       range=(0, 1.0),
                                       normed=True)

            ax.plot([(bins[x] + bins[x+1]) / 2 for x in range(len(bins) - 1)], count, color="red")
            for tick in ax.get_yticklabels():
                tick.set_color('blue')

        if len(dist2) > 0:
            count, bins = np.histogram(dist2,
                                       range=(0, 1.0),
                                       bins=100,
                                       normed=True)

            alt_ax.plot([(bins[x] + bins[x+1]) / 2 for x in range(len(bins) - 1)], count, color="blue")
            for tick in alt_ax.get_yticklabels():
                tick.set_color('red')

        return ax, alt_ax

    def build_genomic_content(self, ax, genomic_locs=None, regions=None):

        """

        builds pie chart describing the total size of things falling inside the genome

        ax - the axis to draw on
        genomic_locs - list of ints 5 with the total size of each genomic region
        stored inside order is Exon, 3'UTR, 5'UTR, Proximal Intron, Distal Intron

        """
        if genomic_locs is None:
            genomic_locs = self.genomic_locs
        if regions is None:
            regions = self.regions

        ax = self.build_pie_chart_content(ax, genomic_locs, regions)
        return ax

    def build_cluster_content(self, ax, cluster_locs=None, regions=None):

        """

        builds pie chart describing the total size of things falling inside the clusters

        ax - the axis to draw on
        cluster_locs - list of ints 5 with the total size of each region in
        clusters stored inside order is Exon, 3'UTR, 5'UTR, Proximal Intron,
        Distal Intron
        """

        if cluster_locs is None:
            cluster_locs = self.clusters_locs
        if regions is None:
            regions = self.regions

        ax = self.build_pie_chart_content(ax, cluster_locs, regions)
        return ax

    def build_pie_chart_content(self, ax, regions_count, regions):

        """

        builds pie chart describing the total size of things falling inside regions

        ax - the axis to draw on
        regions_count - dict of {region : count} of all the regions counted
        regions - dict raw name : good name
        clusters stored inside order is Exon, 3'UTR, 5'UTR, Proximal Intron,
        Distal Intron

        """

        #red, tan, blue, green, purple

        colors = {"cds": "#E52C27",
                  "three_prime_utrs": "#C3996B",
                  "five_prime_utrs": "#3C54A4",
                  "proxintron500": "#48843D",
                  "distintron500": "#852882"}

        ax.pie(regions_count.values(), colors=[colors[region] for region in regions_count.keys()],
               labels=[regions[region] for region in regions_count.keys()])
        return ax

    def build_nearest_exon(self, ax, genomic_types=None, clusters_types=None):

        """

        Creates a bar char showing the the enrichement of the nearest exons relative
        to clusters that are idenitifed

        ax - the axis to draw on
        genomic_types - background distribution of exons in the human genome
        clusters_types - distribuiton of clusters into exons

        order for both lists is CE, SE, MXE, A5E, A3E

        """

        #TODO get rid of this comuptation it should be factored into the model
        if genomic_types is None:
            genomic_types = self.genomic_types
        if clusters_types is None:
            clusters_types = self.clusters_types

        clusters_types = 100 * np.array(clusters_types, dtype="float") / np.sum(clusters_types)
        genomic_types = 100 * np.array(genomic_types, dtype="float") / np.sum(genomic_types)
        difference = clusters_types - genomic_types
        ind = np.arange(5)
        ind -= .5

        ax.bar(ind, difference, color='y')
        ax.set_xticklabels(["", "CE", "SE", "MXE", "A5E", "A3E"])
        ax.set_ylabel('% Difference (clusters-genomic)')
        ax.set_xlabel('Exon Type')
        ax.axhline(y=0, ls="-", color='k')
        return ax

    def build_common_motifs(self, motif_grid, homer_location=None, regions=None):

        """

        Find a the top 8 common motifs in each region as determined by homer
        and places their images on the motif grid

        motif_grid - a grid section to plot the most common motifs onto
        homer_location - the location on the file system where homer is
        (for getting the data, should factor out)

        """
        if homer_location is None:
            homer_location = self.homer_location
        if regions is None:
            regions = self.regions

        for i, region in enumerate(regions.keys()):

            #make a gridspec for the top 8 motifs
            gs_homer_motifs = gridspec.GridSpecFromSubplotSpec(8, 1, subplot_spec=(motif_grid[i]))

            #for each top 8 motifs
            for j, gs in enumerate(gs_homer_motifs):

                #get it from where homer stored the output
                motif_logo = "motif" + str(j + 1) + ".logo.png"
                motif_pvalue = 'motif' + str(j + 1) + '.motif'
                motif_logo_file = os.path.join(homer_location, region, "homerResults", motif_logo)
                motif_pvalue_file = os.path.join(homer_location, region, "homerResults", motif_pvalue)
                if os.path.exists(motif_logo_file) and os.path.exists(motif_pvalue_file):
                    motif = mpimg.imread(motif_logo_file)
                    pvalue = float(open(motif_pvalue_file).readline().strip().split(",")[-1].split(":")[-1])

                    #print title only once
                    if j == 0:
                        ax = plt.subplot(gs, frameon=False,
                                         xticks=[],
                                         yticks=[],
                                         title=regions[region] + "\n" + '{:.2e}'.format(pvalue))
                    else:
                        ax = plt.subplot(gs, frameon=False, xticks=[], yticks=[], title='{:.2e}'.format(pvalue))
                    ax.imshow(motif)
                else:
                    print "no motif %s" % motif_logo_file

    def build_phastcons_values(self, ax, phastcons_values=None, regions=None):

        """

        Builds phastcons CDF plots

        ax - axis to plot phastcons values on
        phastcons_values - clusterfuck of a list that stores CDF for real
        and random values
        regions dict informal name : formal name

        """
        if phastcons_values is None:
            phastcons_values = self.phastcons_values
        if regions is None:
            regions = self.regions

        intersecting_regions = set(phastcons_values['real'].keys()).intersection(set(phastcons_values['rand'].keys()))

        for region in intersecting_regions:
            phastcons_values['real'][region] = [value for value in phastcons_values['real'][region] if not math.isnan(value)]
            phastcons_values['rand'][region] = [value for value in phastcons_values['rand'][region] if not math.isnan(value)]

        box_values = []
        #need to mix the real and random values, one after another
        for region in intersecting_regions:
            box_values.append(phastcons_values['real'][region])
            box_values.append(phastcons_values['rand'][region])

        #positions = np.arange(len(intersecting_regions) * 2)
        width = .1

        #make boxplots of phastcons values for all regions
        bp = ax.boxplot(box_values,
                        #positions = positions,
                        widths=width,
                        sym='k.')

        boxes = bp['boxes']

        #slicing trick, even boxes are true, odd boxses are randomized
        for realbox in boxes[:2]:
            realbox.set_color('blue')

        for randbox in boxes[1::2]:
            randbox.set_color('red')

        #sets all outlighters sizes to 3
        [i.set_markersize(3.) for i in bp['fliers']]

        #sets true medians to blue
        medians = bp['medians']
        for line in medians[:2]:
            line.set_color('blue')

        #sets randomized medians to red
        for line in medians[1::2]:
            line.set_color('red')

        start, stop = ax.get_xlim()
        #start = 0
        #stop += .5
        #ticks to seperate out the different regions...
        ticks = np.linspace(start, stop, len(intersecting_regions) + 1)
        starts = np.array(ticks[:-1])
        bins = 20
        for n, region in enumerate(intersecting_regions):
            heights1, edges1 = np.histogram(phastcons_values['real'][region], normed=True, bins=bins)
            heights2, edges2 = np.histogram(phastcons_values['rand'][region], normed=True, bins=bins)

            #insert 0 at start to anchor CDF
            cdf1 = np.cumsum(heights1) / np.sum(heights1)
            cdf1 = np.insert(cdf1, 0, 0)
            cdf2 = np.cumsum(heights2) / np.sum(heights2)
            cdf2 = np.insert(cdf2, 0, 0)
            yvals = np.linspace(0, 1, (bins + 1))
            if n == 0:
                label1 = "real"
                label2 = "random"
            else:
                label1 = "_nolegend_"
                label2 = "_nolegend_"
            ax.plot((starts[n] + (cdf1 * 2)), yvals, 'blue', label=label1)
            ax.plot((starts[n] + (cdf2 * 2)), yvals, 'red', label=label2)

        for x in starts[1:]:
            ax.axvline(x, color='k', lw=2)

        #sets ticks to be in the middle of box plots
        ax.set_xticks([(ax.get_xticks()[x] + ax.get_xticks()[x+1]) / 2.0 for x in range(len(ax.get_xticks()))[::2]])
        ax.set_xticklabels([regions[region] for region in list(intersecting_regions)])
        ax.set_ylabel("PhastCons Score")
        ax.legend()
        return ax

    def build_motif_boxplots(self, ax, kmer_results=None, highlight_motifs=None, regions=None):

        """

        Make bake boxplots of motif z-scores. you must get kmer z-scores first with run_kmerdiff.
        up to 11 motifs can be highlighted
        pass a pylab subplot instance to the kwarg \"subplot\" to attach this to a figure,
        otherwise it will make its own figure

        kmer_results - dict[region][k][kmer] = motif (motif object defined in kmerdirr)
        the results of calculate_kmer_diff
        highlight_motifs - list[str] motifs to plot

        I'll come back to this, in a way its alrgiht...
        """

        if kmer_results is None:
            kmer_results = self.kmer_results
        if highlight_motifs is None:
            highlight_motifs = self.motifs
        if regions is None:
            regions = self.regions

        colorcycle = ["red", "orange", "green", "blue",
                      "purple", "brown", "black", "pink",
                      "gray", "cyan", "magenta"]

        kmers = {}
        all_kmers = set()

        #might consider converting into pandas data frame
        for region in kmer_results:
            kmers[region] = {}

            #loads all kmers into the allmers set and into a dict
            #print kmer_results
            for k in kmer_results[region].keys():
                for kmer in kmer_results[region][k]:
                    for kmer, motif_data in kmer_results[region][k].items():
                        all_kmers.add(kmer)
                        kmers[region][kmer] = float(motif_data.delta)

        for i, m in enumerate(highlight_motifs):
            #hash the motifs, convert to DNA letters
            highlight_motifs[i] = m.lower().replace("u", "t")

        #creates an ndarray to load all kmer values into
        ak = np.ndarray(shape=(len(all_kmers), len(kmer_results.keys())))

        for i, kmer in enumerate(all_kmers):
            for j, region in enumerate(kmer_results.keys()):
                try:
                    ak[i, j] = kmers[region][kmer]
                #if kmer doesn't exist in specified region its value is zero
                except:
                    ak[i, j] = 0

        #loads everything up and only showed the motifs that are highlighted...
        ax.boxplot(ak, vert=False, notch=1, sym='k.',  whis=2)
        ax.set_yticklabels([regions[region] for region in kmer_results.keys()])
        for i, motif in enumerate(highlight_motifs):
            indices = list()
            for ind, k in enumerate(all_kmers):
                if motif in k:
                    indices.append(ind)
            y = map(lambda x: x+1, range(len(kmer_results.keys())))
            for m, ind in enumerate(indices):
                if m == 0:
                    label = motif
                else:
                    label = None
                ax.plot(ak[ind, :], y, 'o', color=colorcycle[i], label=label, markersize=10)
        ax.set_xscale('symlog', linthreshx=10)
        ax.axvline(x=-4)
        ax.axvline(x=4)

        ax.legend(frameon=False, loc=0, numpoints=1)

        ax.set_xlabel("Z-score")
        return ax

    def plot_distance_and_sn(self, ax, closest_tool, dist_name, normed=False, **kwargs):

        """

        Take a dict of name : bedfile (that has had closest bed run on it) and prints out the distribution 200bp away
        from the closest results

        """
        #TODO Figure out how to fix this, sample labeling needs to be added in and kwargs need to be as well
        #also probably shouldn't do the whole gridspec switching in here, do it somewhere else
        ax.set_title(dist_name)

        sample_distances = np.array([int(result[-1]) for result in closest_tool])
        sample_y = np.histogram(sample_distances, range=(-200, 200), bins=25, normed=normed)[0]
        bins = np.histogram(sample_distances, range=(-200, 200), bins=25)[1]
        x = [(bins[n] + bins[n+1]) / 2 for n in range(len(bins) - 1)]
        ax.set_xticks(np.arange(-200, 201, 100))
        [tick.set_fontsize(10) for tick in ax.get_xticklabels()]
        [tick.set_fontsize(10) for tick in ax.get_yticklabels()]

        ax.plot(x, sample_y, linewidth=3, alpha=.7, **kwargs)
        ax.legend(loc=0)
        sns.despine(ax=ax)
        return ax

    def generate_distribution(self, dist):

        """

        Helper to return his back as a line for easy plotting

        """

        count, bins = np.histogram(dist, range=(0, 1.0), bins=50, density=True)
        return count, [(bins[x] + bins[x+1]) / 2 for x in range(len(bins) - 1)]

    def generate_peak_distribution(self, ax, distribution, total=False, **kwargs):

        """

        ax - axis to plot on
        distribution - distribution to plot (in the form of a list

        """

        if total:
            total_int = 'total'
        else:
            total_int = 'individual'
        sample_count, bins = self.generate_distribution(distribution[total_int])
        ax.plot(bins, sample_count, linewidth=3, alpha=.7, **kwargs)

        ax.legend(loc=0)
        return ax

    def build_peak_densities(self, ax, read_densities=None, classes=None):

        if read_densities is None:
            read_densities = self.read_densities
        if classes is None:
            classes = self.classes

        reordered = read_densities[classes.argsort()]
        cluster = ax.matshow(reordered, aspect='auto', origin='lower')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cluster, cax=cax)

    def plot_distributions(self, ax=None, **kwargs):

        """

        Plots all motifs given in motif distances returns the figure for saving
        motif_distances - list of results from calculate_motif_distance

        tried to make multiple plotting axees, but it turns out this is difficult as gridspec creates new objects,
        which screw with the ability for subplot / add_subplot to get the proper axes
        """
        feature_proper_name_dict = {"transcription_start_sites": "Transcription\nStart Sites",
                                    "start_codons": "Start Codons",
                                    "stop_codons": "Start Codons",
                                    "five_prime_ends": "Five Prime Ends",
                                    "three_prime_ends": "Three Prime Ends",
                                    "poly_a_sites": "Poly A Sites",
                                    'five_prime_utrs': "5' UTR",
                                    'three_prime_utrs': "3' UTR",
                                    'exons': "Exons",
                                    'introns': "Introns",
                                    'cds': 'CDS',
        }

        gs_x = 4
        gs_y = 6
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(11, 8.5))
            full_grid = gridspec.GridSpec(gs_x, gs_y)
        else:
            full_grid = gridspec.GridSpecFromSubplotSpec(gs_x, gs_y, ax.get_subplotspec())
            fig = plt.gcf()

        #premRNA distributions
        gs_premRNA = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=full_grid[0, :])
        for x, (name, region) in enumerate(self.features_transcript_closest.items()):
            print "foo"
            print type(gs_premRNA)
            ax = plt.add_subplot(gs_premRNA[x])
            ax = self.plot_distance_and_sn(ax,
                                           region,
                                           feature_proper_name_dict[name],
                                           **kwargs)

            if ax.is_first_col():
                ax.set_ylabel("Pre-mRNA")
        #mRNA distributions
        gs_mRNA = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=full_grid[1, :])
        for x, (name, region) in enumerate(self.features_mrna_closest.items()):
            ax = plt.subplot(gs_mRNA[x])
            ax = self.plot_distance_and_sn(ax,
                                           region,
                                           feature_proper_name_dict[name],
                                           **kwargs)
            if ax.is_first_col():
                ax.set_ylabel("mRNA")
        #peak distributions across regions
        regions = ['three_prime_utrs', 'five_prime_utrs', 'exons', 'introns', 'cds']
        gs_individual_distribution = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=full_grid[2, :])
        for x, name in enumerate(regions):
            ax = plt.subplot(gs_individual_distribution[x])
            ax.set_title(feature_proper_name_dict[name])
            ax = self.generate_peak_distribution(ax,
                                                 self.distributions[name],
                                                 total=False,
                                                 **kwargs)

            if ax.is_first_col():
                ax.set_ylabel("Individual")

        gs_full_distribution = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=full_grid[3, :])
        for x, name in enumerate(regions):
            ax = plt.subplot(gs_full_distribution[x])
            #"Peaks\n across %s total" %
            ax.set_title(feature_proper_name_dict[name])
            ax = self.generate_peak_distribution(ax,
                                                 self.distributions[name],
                                                 total=True,
                                                 **kwargs)

            if ax.is_first_col():
                ax.set_ylabel("Total")
        fig.tight_layout(pad=3)
        return fig

    def plot_motifs(self, motif_distances):

        """

        Plots all motifs given in motif distances returns the figure for saving
        motif_distances - list of results from calculate_motif_distance

        """

        fig = plt.figure(figsize=(8.5, 11))
        colors = ['red', 'orange', 'green',
                  'blue', 'purple', 'brown',
                  'black', 'pink', 'gray', 'cyan',
                  'magenta']

        for i, motif in enumerate(motif_distances):
            self.plot_motif_dist(motif, fig, color=colors[i], species=species, slopsize=200)

        return fig

    def plot_motif_dist(self, motif_distances, figure, color="red", label=None, scale='linear'):

        """

        Plots distances of motifs from clusters, data structure generated by calculate_motif_distance

        assigned_clusters - dict clusters assigned to specific genic regions + random assignments
        motifFILE - precompiled bed12 /w transcriptome locations for given motif
        figure - output location



        dict{region : {'real': {'size' : int, 'dist' : list[int (distance to motif)},
                       'rand': {'size' : int, 'dist' : list[int (distance to motif)}}

        """

        subplot_number = 320
        for region in motif_distances:
            subplot_number += 1
            ax_region = figure.add_subplot(subplot_number, title=region)
            ax_region.set_yscale(scale)

            region_hist, region_edges = np.histogram(motif_distances['all']['real']['dist'],
                                                     bins=50, range=(-150, 150))
            region_hist /= motif_distances['all']['real']['size'] / 1000.

            region_rand_hist, region_edges_rand = np.histogram(motif_distances['all']['rand']['dist'],
                                                               bins=50, range=(-150, 150))
            region_rand_hist /= motif_distances['all']['rand']['size'] / 1000.

            #plots all motifs on same canvis
            ax_region.plot(region_edges[:-1], region_hist, c=color, linestyle='solid', label=label)
            ax_region.hold(True)
            ax_region.plot(region_edges_rand[:-1], region_rand_hist,
                           linestyle='dashed', c=color, label="_nolegend_")
        #TODO reformat this to taking in a gridspec
        return figure

    def CLIP_QC_figure(self):

        """

        Main Visualization logic

        Input --
        reads_in_clusters -- int, number of reads in a cluster
        reads_out_clusters -- int, number of reads outside of clusters
        cluster_lengths -- list of ints, length of each cluster
        reads_per_cluster -- list of ints, number of reads inside each cluster
        premRNA -- unknown looks like a list of something
        mRNA -- unknown looks like a list of ints for something used for mRNA and pre-mRNA dirstibutions...
        exondist -- distance from exon?
        introndist -- distance from intron?
        genomic_locs -- ?? int
        clusters_locs -- ?? int
        genomic_types --
        clusters_types --
        zscores -- list of zcores for each peak
        homer_location -- location of homer results file
        kmer_results - dict[region][k][(kmer, value) the results of calculate_kmer_diff
        motifs - list[str] motifs to plot
        phastcons_values -- list[float] - conservation scores for each cluster
        regions - dict of short name : printable name
        read_densities - read denisites around peaks, np matrix
        classes - np array, classes resulting from k-means clustering, goes with read denisites

        """

        #First do layout for main figures
        fig = plt.figure(figsize=(20, 20), facecolor='white')

        #The grid for the entire display
        full_grid = gridspec.GridSpec(6, 4)

        #creates Reads in clusters axis
        ax_pie_inout = plt.subplot(full_grid[0, 0], title="Reads in Clusters", aspect=1)

        #Lay out the rest of the first column
        gs_nreads_lengths_cons = gridspec.GridSpecFromSubplotSpec(1, 7, subplot_spec=full_grid[0, 1:4])

        #Creates reads per cluster axis
        ax_nreads = plt.subplot(gs_nreads_lengths_cons[0:3], title="Reads per Cluster")

        #creates cluster lengths axis
        ax_lengths = plt.subplot(gs_nreads_lengths_cons[3], title="Cluster Lengths")

        #Creates phastcons axis (lets modify this a bit to have it be its own thing
        ax_cons = plt.subplot(gs_nreads_lengths_cons[4:7], title="PhastCons Values")

        #Lays out second column
        gs_line2 = gridspec.GridSpecFromSubplotSpec(1, 6, subplot_spec=full_grid[1, :])

        #Creates gene distance axis
        ax_genedist = plt.subplot(gs_line2[0:2], title="mRNA and pre-mRNA Distribution")

        #Creates exon distance axis
        ax_exondist = plt.subplot(gs_line2[2:4], title="Distribution Across Exons and Introns")

        #Creates read densities axis
        ax_read_densities = plt.subplot(gs_line2[4:], title="Read Densities Around Peaks")

        #Creates out 3rd column minus the motif z-scores
        gs_pie_nearestType = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=full_grid[2, 0:2])

        #Creates genomic content axis
        ax_pie_genomic = plt.subplot(gs_pie_nearestType[0], title="Genomic Content", aspect=1)

        #Creates exonic content axis
        ax_pie_exonic = plt.subplot(gs_pie_nearestType[1], title="Exonic Content", aspect=1)

        #Creates cluster content axis
        ax_pie_clusters = plt.subplot(gs_pie_nearestType[2], title="Clusters' Content", aspect=1)

        #Creates nearest exon type axis
        ax_bar_exontypes = plt.subplot(gs_pie_nearestType[3], title="Nearest Exon Types")

        #Creates z-score axis
        ax_hist_zscores = plt.subplot(full_grid[2,2:4], title="Motif Z-scores")

        #The grid for just the motif display section
        motif_grid = gridspec.GridSpecFromSubplotSpec(1, 6,
                                                      subplot_spec=full_grid[4:6, :],
                                                      hspace=0,
                                                      wspace=0)

        #Pass specific axies + data off for plotting
        self.build_reads_in_clusters(ax_pie_inout)
        self.build_reads_per_cluster(ax_nreads)
        self.build_cluster_lengths(ax_lengths)
        if self.phastcons_values is not None:
            self.build_phastcons_values(ax_cons)
        self.build_gene_distribution(ax_genedist)
        self.build_exon_exon_distribution(ax_exondist)
        self.build_peak_densities(ax_read_densities)
        self.build_genomic_content(ax_pie_genomic)

        #filter out intronic regions
        exonic_locs = {name: value for name, value in self.genomic_locs.items() if name not in ['proxintron500',
                                                                                           'distintron500']}
        self.build_genomic_content(ax_pie_exonic, exonic_locs, self.regions)

        self.build_cluster_content(ax_pie_clusters)
        self.build_nearest_exon(ax_bar_exontypes)
        try:
            self.build_common_motifs(motif_grid)
        except AttributeError:
            print "Can't Build Homer Motifs"
            pass

        if self.kmer_results is not None:
            self.build_motif_boxplots(ax_hist_zscores)
        plt.tight_layout()
        return fig
