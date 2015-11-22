'''
Created on Sep 18, 2012

@author: gabrielp
'''


import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['interactive'] = False

import cPickle as pickle
from collections import OrderedDict
import math
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
import statsmodels as sm



class ClipVisualization():

    def __init__(self, pickle_file=None):
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

        if pickle_file is None:
            return

        data = pickle.load(pickle_file)

        self.reads_in_clusters = data['reads_in_clusters'] if "reads_in_clusters" in data else None
        self.reads_out_clusters = data['reads_out_clusters'] if "reads_out_clusters" in data else None
        self.cluster_lengths = data['cluster_lengths'] if "cluster_lengths" in data else None
        self.reads_per_cluster = data['reads_per_cluster'] if "reads_per_cluster" in data else None
        self.premRNA = data['distributions']['genes']['total'] if "distributions" in data else None
        self.mRNA = data['distributions']['exons']['total'] if "distributions" in data else None
        self.exondist = data['distributions']['exons']['individual'] if "distributions" in data else None
        self.introndist = data['distributions']['introns']['individual'] if "distributions" in data else None
        self.distributions = data['distributions'] if "distributions" in data else None
        self.genomic_locs = data['genic_region_sizes'] if "genic_region_sizes" in data else None
        self.clusters_locs = data['region_sizes'] if "region_sizes" in data else None
        self.genomic_types = data['genomic_type_count'] if "genomic_type_count" in data else None
        self.clusters_types = data['type_count'] if "type_count" in data else None
        self.kmer_results = data['kmer_results'] if "kmer_results" in data else None
        self.motifs = data['motifs'] if "motifs" in data else None
        self.phastcons_values = data['phast_values'] if "phast_values" in data else None
        self.read_densities = np.array(data['data']) if "data" in data else None
        self.classes = data['classes'] if "classes" in data else None
        self.features_transcript_closest = data['features_transcript_closest'] if "features_transcript_closest" in data else None
        self.features_mrna_closest = data['features_mrna_closest'] if "features_mrna_closest" in data else None
        self.homer_location = data['homerout']

    def plot_cdf(self, cdf_list, **kwargs):
        cdf = sm.distributions.ECDF(cdf_list)
        cdf_linspace = np.linspace(min(cdf_list), max(cdf_list))
        if kwargs['ax'] is not None:
            ax = kwargs['ax']
            del kwargs['ax']
            ax.plot(cdf_linspace, cdf(cdf_linspace), **kwargs)
            ax.set_ylim((0, 1))
        else:
            plot(cdf_linspace, cdf(cdf_linspace), **kwargs)

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

        if reads_out_clusters is None or reads_out_clusters is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

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

        if reads_per_cluster is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

        sns.kdeplot(np.array(self.reads_per_cluster), ax=ax_nreads)
        [tick.set_rotation(90) for tick in ax_nreads.get_xticklabels()]
        ax_nreads.set_xlim(0,)

        ax_nreads.set_xlabel("N reads)")
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

        if cluster_lengths is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

        sns.kdeplot(np.array(self.cluster_lengths), ax=ax_lengths)
        [tick.set_rotation(90) for tick in ax_lengths.get_xticklabels()]
        ax_lengths.set_xlim(0,)
        ax_lengths.set_ylabel("Frequency")
        ax_lengths.set_xlabel("Length (bp)")
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

        if premrna is None or mrna is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

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

        if exon is None or intron is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

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

        if genomic_locs is None or regions is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

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

        if cluster_locs is None or regions is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

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

        values = []
        color_list = []
        labels = []
        for region in regions_count.viewkeys() & regions.viewkeys() & colors.viewkeys():
            values.append(regions_count[region])
            color_list.append(colors[region])
            labels.append(regions[region])
        ax.pie(values, colors=color_list, labels=labels)
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

        if genomic_types is None or clusters_types is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

        clusters_types = 100 * np.array(clusters_types, dtype="float") / np.sum(clusters_types)
        genomic_types = 100 * np.array(genomic_types, dtype="float") / np.sum(genomic_types)
        difference = clusters_types - genomic_types
        ind = np.arange(5.0)
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

        if homer_location is None or regions is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

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

    def build_phastcons_values(self, gs, phastcons_values=None, regions=None):

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

        if phastcons_values is None or regions is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

        non_null_regions = self.phastcons_values.index.levels[0] & regions.keys()

        for x, (region) in enumerate(non_null_regions):
            print non_null_regions
            ax = plt.subplot(gs[x])
            self.plot_cdf(self.phastcons_values.ix[region, 'real'].mean0, ax=ax, label="Actual")
            self.plot_cdf(self.phastcons_values.ix[region, 'rand'].mean0, ax=ax, label="Random")
            tmp_lim = ax.get_ylim()
            tmp_yticklabels = ax.get_yticks()
            first = (tmp_lim[1] - tmp_lim[0]) * .25
            second = (tmp_lim[1] - tmp_lim[0]) * .75

            bp = ax.boxplot([self.phastcons_values.ix[region, 'real'].mean0,
                             self.phastcons_values.ix[region, 'rand'].mean0],
                            widths=0.05,
                            sym='k.',
                            patch_artist=True,
                            vert=False,
                            positions =[first, second])
            plt.setp(bp['boxes'][0], color='b', alpha=.7)
            plt.setp(bp['boxes'][1], color='g', alpha=.7)
            ax.set_yticks(tmp_yticklabels)
            ax.set_yticklabels(tmp_yticklabels, fontsize=6)
            ax.set_ylim(tmp_lim)
            sns.despine()
            if ax.is_first_col():
                ax.set_ylabel("Cumulative\nFraction", fontsize=8)
            else:
                ax.set_yticklabels([])

            if ax.is_last_col() and ax.is_first_row():
                ax.legend()
            x_label = regions[region].replace("\n", " ")

            if ax.is_last_row():
                x_label += "\nConservation Score"
            ax.set_xlabel(x_label, fontsize=8)


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

        if kmer_results is None or highlight_motifs is None or regions is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

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

        if read_densities is None or classes is None:
            raise NotImplementedError("Pickle file doesn't have data to generate this figure")

        reordered = read_densities[classes.argsort()]
        cluster = ax.matshow(reordered, aspect='auto', origin='lower')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cluster, cax=cax)

    def plot_distributions(self, ax=None, **kwargs):

        """

        Plots all motifs given in motif distances returns the figure for saving
        motif_distances - list of results from calculate_motif_distance

        """
        feature_proper_name_dict = {"transcription_start_sites": "Transcription\nStart Sites",
                                    "start_codons": "Start Codons",
                                    "stop_codons": "Stop Codons",
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
            ax = plt.subplot(gs_premRNA[x])
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
        fig = plt.figure(figsize=(20, 20))

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
        gs_cons = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs_nreads_lengths_cons[4:7])

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
        try:
            self.build_reads_in_clusters(ax_pie_inout)
        except NotImplementedError:
            pass
        try:
            self.build_reads_per_cluster(ax_nreads)
        except NotImplementedError:
            pass
        try:
            self.build_cluster_lengths(ax_lengths)
        except NotImplementedError:
            pass
        try:
            self.build_phastcons_values(gs_cons)
        except NotImplementedError:
            pass
        try:
            self.build_gene_distribution(ax_genedist)
        except NotImplementedError:
            pass
        try:
            self.build_exon_exon_distribution(ax_exondist)
        except NotImplementedError:
            pass
        try:
            self.build_peak_densities(ax_read_densities)
        except NotImplementedError:
            pass
        try:
            self.build_genomic_content(ax_pie_genomic)
        except NotImplementedError:
            pass

        #filter out intronic regions
        exonic_locs = {name: value for name, value in self.genomic_locs.items() if name not in ['proxintron500',
                                                                                           'distintron500']}
        try:
            self.build_genomic_content(ax_pie_exonic, exonic_locs, self.regions)
        except NotImplementedError:
            pass
        try:
            self.build_cluster_content(ax_pie_clusters)
        except NotImplementedError:
            pass
        try:
            self.build_nearest_exon(ax_bar_exontypes)
        except NotImplementedError:
            pass
        try:
            self.build_common_motifs(motif_grid)
        except AttributeError:
            print "Can't Build Homer Motifs"
            pass

        if self.kmer_results is not None:
            self.build_motif_boxplots(ax_hist_zscores)
        plt.tight_layout()
        return fig
