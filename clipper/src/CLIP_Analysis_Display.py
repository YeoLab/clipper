'''
Created on Sep 18, 2012

@author: gabrielp
'''

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

def build_reads_in_clusters(ax, reads_in_clusters, reads_out_clusters):
    
    """
    
    Draws the number of reads in a cluster and the number of reads 
    outside of clusters on a pie chart
    
    ax - the axis to draw on  
    reads_in_clusters - int, the number of reads falling into a cluster
    reads_out_clusters - int the number of reads falling outside of a cluster
    
    """

    ax.pie([reads_in_clusters, reads_out_clusters], 
           labels=["In Clusters", "Outside of\nClusters"])



def build_reads_per_cluster(ax_nreads, reads_per_cluster):

    """
    
    Draws the number of reads per cluster for each cluster 
    
    ax - the axis to draw on  
    reads_per_cluster - list, the number of reads in a cluster

    
    """
    
    ax_nreads.hist(reads_per_cluster, 
                   bins=50, 
                   facecolor='#C8D2B0', 
                   log=True, 
                   range=(10, np.max(reads_per_cluster)))
    
    ax_nreads.set_xscale('log')
    ax_nreads.set_xlabel("log10(N reads)")
    ax_nreads.set_ylabel("Frequency")
    

def build_cluster_lengths(ax_lengths, cluster_lengths):
    
    """
    
    Selects a random sample of all cluster length 
    and draws 2000 of them in a boxplot 
    
    ax - the axis to draw on  
    cluster_lengths - list, the length of each cluster

    
    """
    
    ax_lengths.set_yscale('log')
    ax_lengths.boxplot(random.sample(cluster_lengths, 
                                     min(2000, len(cluster_lengths))), 
                       vert=1)
    
    ax_lengths.set_xticklabels([])


def build_gene_distribution(ax, premRNA, mRNA):
    
    """
    
    Draws the distribution of clusters across mRNA and pre-mRNA meta-genes
    
    ax - the axis to draw on
    premRNA - a list of locations (in percent) that describes where a cluster binds
    in preMRNA 
    
    mRNA - a list of locations (in percent) that describes where a cluster binds
    in mRNA  
    
    """
    
    alt_ax = build_distribution(ax, premRNA, mRNA)
 
    ax.set_xlabel("Fraction of region")
    ax.set_ylabel("pre-mRNA Location Frequency", color = 'blue')
    alt_ax.set_ylabel("mRNA Location Frequency", color = 'red')

def build_exon_exon_distribution(ax, exon, intron):
    
    """
    
    Draws the distribution of clusters across exon and intron meta-genes
    
    ax - the axis to draw on
    exon - a list of locations (in percent) that describes where a cluster binds
    in exons
    
    intron - a list of locations (in percent) that describes where a cluster binds
    in introns  
    
    """
    
    alt_ax = build_distribution(ax, exon, intron)
    
    ax.set_xlabel("Fraction of region")
    ax.set_ylabel("Exon Location Frequency", color='blue')
    alt_ax.set_ylabel("Intron Location Frequency", color='red')


def build_distribution(ax, dist1, dist2):
    
    """
    
    Plots two histograms on the same axis returning a second axis for labeling purposes.  
    ax - the axis to draw on
    
    dist1 - the first distribution to plot
    dist2 - the second distribution to plot
    
    """
    

    alt_ax = ax.twinx()
    ax.hist(dist1, range=(0, 1.0), 
            histtype="step", 
            color="red", 
            bins=100, 
            normed=True)
    
    alt_ax.hist(dist2, 
                range=(0, 1.0), 
                histtype="step", 
                color="blue", 
                bins=100, 
                normed=True)

    for tick in alt_ax.get_yticklabels():
        tick.set_color('red')
    
    for tick in ax.get_yticklabels():
        tick.set_color('blue')
        
    return alt_ax


def build_genomic_content(ax, genomic_locs):
    
    """
    
    builds pie chart describing the total size of things falling inside the genome
    
    ax - the axis to draw on
    genomic_locs - list of ints 5 with the total size of each genomic region 
    stored inside order is Exon, 3'UTR, 5'UTR, Proximal Intron, Distal Intron
    
    """
    
    build_pie_chart_content(ax, genomic_locs)

def build_cluster_content(ax, cluster_locs):
    
    """
    
    builds pie chart describing the total size of things falling inside the clusters
    
    ax - the axis to draw on
    cluster_locs - list of ints 5 with the total size of each region in 
    clusters stored inside order is Exon, 3'UTR, 5'UTR, Proximal Intron, 
    Distal Intron
    """
    
    build_pie_chart_content(ax, cluster_locs)

def build_pie_chart_content(ax, regions):
    
    """
    
    builds pie chart describing the total size of things falling inside regions
    
    ax - the axis to draw on
    clusters_locs - list of ints 5 with the total size of each regions in 
    clusters stored inside order is Exon, 3'UTR, 5'UTR, Proximal Intron, 
    Distal Intron
    
    """

    #red, tan, blue, green, purple
    colors = ["#E52C27", "#C3996B", "#3C54A4", "#48843D", "#852882"] 
    labels = ["Exon", "3'UTR", "5'UTR", "Proximal\nIntron", "Distal\nIntron"]
    
    ax.pie(regions, colors=colors, labels=labels)
    

def build_nearest_exon(ax, genomic_types, clusters_types):
    
    """
    
    Creates a bar char showing the the enrichement of the nearest exons relative 
    to clusters that are idenitifed 
    
    ax - the axis to draw on
    genomic_types - background distribution of exons in the human genome
    clusters_types - distribuiton of clusters into exons 
    
    order for both lists is CE, SE, MXE, A5E, A3E
    
    """
    
    #TODO get rid of this comuptation it should be factored into the model
    clusters_types = 100 * np.array(clusters_types, dtype="float") / np.sum(clusters_types)
    genomic_types  = 100 * np.array(genomic_types, dtype="float") / np.sum(genomic_types)
    difference = clusters_types - genomic_types
    ind = np.arange(5)
    ind = ind - .5
    
    ax.bar(ind, difference, color='y')
    ax.set_xticklabels(["", "CE", "SE", "MXE", "A5E", "A3E"])
    ax.set_ylabel('% Difference (clusters-genomic)')
    ax.set_xlabel('Exon Type')
    ax.axhline(y=0, ls="-", color='k')


def build_common_motifs(motif_grid, homer_location):
    
    """
    
    Find a the top 8 common motifs in each region as determined by homer 
    and places their images on the motif grid
    
    motif_grid - a grid section to plot the most common motifs onto
    homer_location - the location on the file system where homer is 
    (for getting the data, should factor out)
    
    """
    
    all_regions = ["all", 
                    "exon", 
                    "UTR3", 
                    "UTR5", 
                    "proxintron", 
                    "distintron",
                    ]
    
    formal_labels = ["All", 
                     "Exon", 
                     "3'UTR", 
                     "5'UTR", 
                     "Proximal Intron", 
                     "Distal Intron",
                     ]
    
    #for reach region 
    for i, region in enumerate(all_regions):
        
        #make a gridspec for the top 8 motifs 
        gs_homer_motifs = gridspec.GridSpecFromSubplotSpec(8, 1, subplot_spec=(motif_grid[i]))
        
        #for each top 8 motifs
        for j, space in enumerate(gs_homer_motifs):

            #get it from where homer stored the output
            motif_name = "motif" + str(j + 1) + ".logo.png"
            
            motif_file = os.path.join(homer_location, region, "homerResults", motif_name)
            
            if os.path.exists(motif_file):
                motif = mpimg.imread(motif_file)
            
                #print title only once
                if j == 0:
                    title = formal_labels[i]
                    ax = plt.subplot(space, frameon=False, xticks=[], yticks=[], title=title)
                else:
                    ax = plt.subplot(space, frameon=False, xticks=[], yticks=[])
                ax.imshow(motif)
            else:
                print "no motif %s" % (motif_file)



def build_phastcons_values(ax, phastcons_values):
    
    """
    
    Builds phastcons CDF plots
    
    ax - axis to plot phastcons values on
    phastcons_values - clusterfuck of a list that stores CDF for real
    and random values
    
    TODO: Fix bad data structures and get working
    
    """
    
    formal_labels = ["All", 
                     "Exon", 
                     "3'UTR", 
                     "5'UTR", 
                     "Proximal\nIntron", 
                     "Distal\nIntron",
                     ]
    
    phastcons_values = [[value for value in lst if not math.isnan(value)] for lst in phastcons_values]
    
    print len(phastcons_values)
    for value in phastcons_values:
        print len(value)
            
    positions = np.arange(5) - .1
    width = .1
    
  
    #make boxplots of phastcons values for all regions
    bp = ax.boxplot([range(10)], 
                    positions = positions, 
                    widths = width, 
                    sym = 'k.')
    
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
    
    #sets randomized medians to false
    for line in medians[1::2]:
        line.set_color('red')
    
    start, stop = ax.get_xlim()
    ticks = np.linspace(start, stop, 7)
    starts = ticks[:-1]
    stops = ticks[1:]
    reals = phastcons_values[::2]
    rands = phastcons_values[1::2]
    bins = 20
    for n in range(6):
        heights1, edges1 = np.histogram(reals[n], normed=True, bins=bins)
        heights2, edges2 = np.histogram(rands[n], normed=True, bins=bins)
        CDF1 = np.cumsum(heights1) / np.sum(heights1)
        CDF1 = np.insert(CDF1, 0, 0)
        CDF2 = np.cumsum(heights2) / np.sum(heights2)
        CDF2 = np.insert(CDF2, 0, 0)

        yvals = np.linspace(0, 1, (bins + 1))
        if n == 0:
            label1 = "real"
            label2 = "random"
        else:
            label1 = "_nolegend_"
            label2 = "_nolegend_"
        ax.plot(starts[n] + CDF1, yvals, 'blue', label=label1)
        ax.plot(starts[n] + CDF2, yvals, 'red', label=label2)
    
    for x in starts[1:]:
        ax.axvline(x, color='k', lw=2)
    
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(formal_labels)
    ax.set_ylabel("PhastCons Score")

def CLIP_QC_figure(reads_in_clusters, reads_out_clusters, cluster_lengths, 
                   reads_per_cluster, premRNA, mRNA, exondist, introndist, 
                   genomic_locs, clusters_locs, genomic_types, clusters_types,
                   zscores, homer_location, kmer_box_params, phastcons_values):
    
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
    kmer_box_params -- ??
    phastcons_values -- ??
    
    """
    
    #First do layout for main figures
    fig = plt.figure(figsize=(20, 20), facecolor='white')
    
    #The grid for the entire display
    full_grid = gridspec.GridSpec(6, 4)
    
    #creates Reads in clusters axis
    ax_pie_inout = plt.subplot(full_grid[0, 0], title = "Reads in Clusters", aspect=1)
    
    #Lay out the rest of the first column
    gs_nreads_lengths_cons = gridspec.GridSpecFromSubplotSpec(1, 7, subplot_spec=full_grid[0, 1:4])
    
    #Creates reads per cluster axis
    ax_nreads = plt.subplot(gs_nreads_lengths_cons[0:3], title = "Reads per Cluster")
    
    #creates cluster lengths axis
    ax_lengths = plt.subplot(gs_nreads_lengths_cons[3], title = "Cluster Lengths")
    
    #Creates phastcons axis (lets modify this a bit to have it be its own thing
    ax_cons = plt.subplot(gs_nreads_lengths_cons[4:7], title = "PhastCons Values")
    
    #Lays out second column 
    gs_line2 = gridspec.GridSpecFromSubplotSpec(1,6, subplot_spec=full_grid[1,:])
    
    #Creates gene distance axis
    ax_genedist = plt.subplot(gs_line2[0:2], title = "mRNA and pre-mRNA Distribution")
    
    #Creates exon distance axis
    ax_exondist = plt.subplot(gs_line2[2:4], title = "Distribution Across Exons and Introns")
    
    #Creates out 3rd column minus the motif z-scores
    gs_pie_nearestType = gridspec.GridSpecFromSubplotSpec(1,3, subplot_spec=full_grid[2,0:2])
    
    #Creates genomic content axis 
    ax_pie_genomic = plt.subplot(gs_pie_nearestType[0], title = "Genomic Content", aspect=1)
    
    #Creates cluster content axis 
    ax_pie_clusters = plt.subplot(gs_pie_nearestType[1], title = "Clusters' Content", aspect=1)
    
    #Creates nearest exon type axis
    ax_bar_exontypes = plt.subplot(gs_pie_nearestType[2], title = "Nearest Exon Types")  
    
    #Creates z-score axis 
    ax_hist_zscores = plt.subplot(full_grid[2,2:4], title = "Motif Z-scores")
    
    #The grid for just the motif display section
    motif_grid = gridspec.GridSpecFromSubplotSpec(1, 6, 
                                                  subplot_spec = full_grid[4:6,:], 
                                                  hspace = 0, 
                                                  wspace = 0)

    #Pass specific axies + data off for plotting
    build_reads_in_clusters(ax_pie_inout, reads_in_clusters, reads_out_clusters)
    build_reads_per_cluster(ax_nreads, reads_per_cluster)
    build_cluster_lengths(ax_lengths, cluster_lengths)
    #build_phastcons_values(ax_cons, phastcons_values)
    #TODO add in z-score stuff again
    build_gene_distribution(ax_genedist, premRNA, mRNA)
    build_exon_exon_distribution(ax_exondist, exondist, introndist)
    build_genomic_content(ax_pie_genomic, genomic_locs)
    build_cluster_content(ax_pie_clusters, clusters_locs)
    build_nearest_exon(ax_bar_exontypes, genomic_types, clusters_types)
    build_common_motifs(motif_grid, homer_location)
            
    plt.tight_layout()
    return fig

def plot_motif_dist(assigned_clusters, motifFILE, figure, nrand=3, color = "red", label=None, species="mm9", slopsize=0, scale='linear'):
    
    """
    
    Plotting function, don't want touch this with a 10 foot pole, can be refactored easily
    Takes output of get_motif_distance and plots for each different genic region
    
    assigned_clusters - dict clusters assigned to specific genic regions + random assignments 
    motifFILE - precompiled bed12 /w transcriptome locations for given motif  
    figure - output location
    nrand - number of randomization steps
    slopsize - size around clusters to look for other motifs (windowed bed)
    
    """
    
    motifBed = pybedtools.BedTool(motifFILE)
    if label is None:
        label=motifFILE
    UTR5dist = get_motif_distance(assigned_clusters['UTR5']['real'], motifBed, slop=slopsize)
    rand_5UTR = list()

    for i in range(nrand):
        rand_5UTR.extend(get_motif_distance(assigned_clusters['UTR5']['rand'][i], motifBed, slop=slopsize))

    UTR5size = assigned_clusters['UTR5']['real'].total_coverage()
    UTR5_rand_size = UTR5size*nrand
    print "UTR5 done"
    UTR3dist = get_motif_distance(assigned_clusters['UTR3']['real'], motifBed, slop=slopsize)
    rand_3UTR = list()
    for i in range(nrand):
        rand_3UTR.extend(get_motif_distance(assigned_clusters['UTR3']['rand'][i], motifBed, slop=slopsize))

    UTR3size = assigned_clusters['UTR3']['real'].total_coverage()
    UTR3_rand_size = UTR3size*nrand            
    print "UTR3 done"
    exondist = get_motif_distance(assigned_clusters['exon']['real'], motifBed, slop=slopsize)
    rand_exon = list()
    for i in range(nrand):
        rand_exon.extend(get_motif_distance(assigned_clusters['exon']['rand'][i], motifBed, slop=slopsize))

    exonsize = assigned_clusters['exon']['real'].total_coverage()
    exon_rand_size = exonsize*nrand            
    print "exon done"

    distintrondist = get_motif_distance(assigned_clusters['distintron']['real'], motifBed, slop=slopsize)
    rand_distintron = list()
    for i in range(nrand):
        rand_distintron.extend(get_motif_distance(assigned_clusters['distintron']['rand'][i], motifBed, slop=slopsize))

    distintronsize = assigned_clusters['distintron']['real'].total_coverage()
    distintron_rand_size = distintronsize*nrand                        
    print "distintron done"
    
    proxintrondist = get_motif_distance(assigned_clusters['proxintron']['real'], motifBed, slop=slopsize)
    rand_proxintron = list()
    for i in range(nrand):
        rand_proxintron.extend(get_motif_distance(assigned_clusters['proxintron']['rand'][i], motifBed, slop=slopsize))

    proxintronsize = assigned_clusters['proxintron']['real'].total_coverage()
    proxintron_rand_size = proxintronsize*nrand

    print "proxintron done"

    allsize = UTR5size + UTR3size + exonsize + proxintronsize + distintronsize
    all_rand_size = allsize*nrand

    all = list()
    all.extend(UTR5dist)
    all.extend(UTR3dist)
    all.extend(exondist)
    all.extend(proxintrondist)
    all.extend(distintrondist)                    

    all_rand = list()
    all_rand.extend(rand_5UTR)
    all_rand.extend(rand_3UTR)
    all_rand.extend(rand_exon)
    all_rand.extend(rand_distintron)
    all_rand.extend(rand_proxintron)    

    ax_all = figure.add_subplot(321, title="All Clusters")
    ax_all.set_yscale(scale)
    ax_UTR5 = figure.add_subplot(322, title="5'UTR")
    ax_UTR5.set_yscale(scale)
    ax_exon = figure.add_subplot(323, title="Exon")
    ax_exon.set_yscale(scale)
    ax_UTR3 = figure.add_subplot(324, title="3'UTR")
    ax_UTR3.set_yscale(scale)
    ax_proxintron = figure.add_subplot(325, title="Proximal Intron")
    ax_proxintron.set_yscale(scale)
    ax_distintron = figure.add_subplot(326, title="Distal Intron")
    ax_distintron.set_yscale(scale)

    all_hist, all_edges = np.histogram(all, bins=50, range=(-150, 150))
    all_hist = all_hist/(allsize/1000.)        
    all_rand_hist, all_edges_rand = np.histogram(all_rand, bins=50, range=(-150, 150))
    all_rand_hist = all_rand_hist/(all_rand_size/1000.)

    ax_all.plot(all_edges[:-1], all_hist, c=color, linestyle='solid', label=label)
    ax_all.plot(all_edges_rand[:-1], all_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    UTR5_hist, UTR5_edges = np.histogram(UTR5dist, bins=50, range=(-150, 150))
    UTR5_hist = UTR5_hist/(UTR5size/1000.)        
    UTR5_rand_hist, UTR5_edges_rand = np.histogram(rand_5UTR, bins=50, range=(-150, 150))
    UTR5_rand_hist = UTR5_rand_hist/(UTR5_rand_size/1000.)       
    ax_UTR5.plot(UTR5_edges[:-1], UTR5_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_UTR5.hold(True)
    ax_UTR5.plot(UTR5_edges_rand[:-1], UTR5_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    UTR3_hist, UTR3_edges = np.histogram(UTR3dist, bins=50, range=(-150, 150))
    UTR3_hist = UTR3_hist/(UTR3size/1000.)        
    UTR3_rand_hist, UTR3_edges_rand = np.histogram(rand_3UTR, bins=50, range=(-150, 150))
    UTR3_rand_hist = UTR3_rand_hist/(UTR3_rand_size/1000.)       
    ax_UTR3.plot(UTR3_edges[:-1], UTR3_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_UTR3.hold(True)
    ax_UTR3.plot(UTR3_edges_rand[:-1], UTR3_rand_hist, linestyle='dashed', c=color, label="_nolegend_")    

    exon_hist, exon_edges = np.histogram(exondist, bins=50, range=(-150, 150))
    exon_hist = exon_hist/(exonsize/1000.)        
    exon_rand_hist, exon_edges_rand = np.histogram(rand_exon, bins=50, range=(-150, 150))
    exon_rand_hist = exon_rand_hist/(exon_rand_size/1000.)       
    ax_exon.plot(exon_edges[:-1], exon_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_exon.hold(True)
    ax_exon.plot(exon_edges_rand[:-1], exon_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    distintron_hist, distintron_edges = np.histogram(distintrondist, bins=50, range=(-150, 150))
    distintron_hist = distintron_hist/(distintronsize/1000.)        
    distintron_rand_hist, distintron_edges_rand = np.histogram(rand_distintron, bins=50, range=(-150, 150))
    distintron_rand_hist = distintron_rand_hist/(distintron_rand_size/1000.)       
    ax_distintron.plot(distintron_edges[:-1], distintron_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_distintron.hold(True)
    ax_distintron.plot(distintron_edges_rand[:-1], distintron_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    proxintron_hist, proxintron_edges = np.histogram(proxintrondist, bins=50, range=(-150, 150))
    proxintron_hist = proxintron_hist/(proxintronsize/1000.)        
    proxintron_rand_hist, proxintron_edges_rand = np.histogram(rand_proxintron, bins=50, range=(-150, 150))
    proxintron_rand_hist = proxintron_rand_hist/(proxintron_rand_size/1000.)       
    ax_proxintron.plot(proxintron_edges[:-1], proxintron_hist, linestyle='solid', c=color, label="_nolegend_")
    ax_proxintron.hold(True)
    ax_proxintron.plot(proxintron_edges_rand[:-1], proxintron_rand_hist, linestyle='dashed', c=color, label="_nolegend_")

    return

def motif_boxplots(kmerloc, filename, klengths, highlight_motifs, subplot=None):

    """
    
    Make bake boxplots of motif z-scores. you must get kmer z-scores first with run_kmerdiff.  up to 11 motifs can be highlighted
    pass a pylab subplot instance to the kwarg \"subplot\" to attach this to a figure, otherwise it will make its own figure
    
    Fix this up later...
    
    """
    
    colorcycle = ["red", "orange", "green", "blue", "purple", "brown", "black", "pink", "gray", "cyan", "magenta"]
    regions = ["all", "exon", "UTR3", "UTR5", "proxintron", "distintron"]
    formal_labels = ["All Regions", "Exon", "3'UTR", "5'UTR", "Proximal Intron", "Distal Intron"]    
    kmers = {}
    all_kmers = set()
    for region in regions:
        kmers[region] = {}
        for k in klengths:
            if(os.path.exists(os.path.join(kmerloc, "%s.k%s.%s.kmerdiff.sort" %(filename, str(k), region)))):
                f = open(os.path.join(kmerloc, "%s.k%s.%s.kmerdiff.sort" %(filename, str(k), region)))  
                for line in f:
                    kmer, val = line.strip().split("\t")
                    kmer, val = map(str.strip, [kmer, val])
                    all_kmers.add(kmer)
                    kmers[region][kmer] = float(val)

    for i, m in enumerate(highlight_motifs):
        #hash the motifs, convert to DNA letters
        highlight_motifs[i] = m.lower().replace("u", "t")
    ak = np.ndarray(shape=(len(all_kmers), len(regions)))
    all_kmers = list(all_kmers)

    for i, kmer in enumerate(all_kmers):
        for j, region in enumerate(regions):            
            ak[i,j] = kmers[region][kmer]
    
    showme=False
    if subplot is None: #
        showme=True
        x = plt.figure()
        subplot = x.add_subplot(111)

    subplot.boxplot(ak, vert=False, notch=1, sym='k.',  whis=2)
    subplot.set_yticklabels(formal_labels)
    for i, motif in enumerate(highlight_motifs):
        indices= list()
        for ind, k in enumerate(all_kmers):
            if motif in k:
                indices.append(ind)
        y = map(lambda x: x+1, range(len(regions)))
        for m, ind in enumerate(indices):
            if m ==0:
                label=motif
            else:
                label=None
            subplot.plot(ak[ind,:], y, 'o', color=colorcycle[i], label=label, markersize=10)
    subplot.set_xscale('symlog', linthreshx=10)
    subplot.axvline(x=-4)
    subplot.axvline(x=4)
    

    subplot.legend(frameon=False,loc=0, numpoints=1)
    
    subplot.set_xlabel("Z-score")
    if showme is True:
        plt.show()    

    return ak, all_kmers

def get_motif_distance(clusters, motif, slop=500):
    
    """
    
    TODO: This shouldn't be in the visualization side of things, need to factor out
    
    Compares two bed files and computes distance from center of first (indicated by bed12)
    to center of second (by bed12)
    

    Gets offsets for each cluster

    Input:
      
    clusters - bedtool (bed12)
    motif - bedtool (bed12)
    
    returns distance from clusters to nearest motif 
    
    """
    
    ov = clusters.window(motif, w=slop, sm=True)
    distances = list()
    for line in ov:
        positions=str(line).split("\t")
        cluster_center = int(positions[7])-int(positions[6])/2
        motif_center = int(positions[15]) - int(positions[14])/2
        distance = motif_center - cluster_center
        if positions[5] == "-":
            distance = distance * -1
        distances.append(distance)
    del ov
    return distances
