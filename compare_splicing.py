import scipy.stats
from optparse import OptionParser
import pickle
import numpy as np
from splicing import retrieve_splicing


def mergeSamples(samples, splicetypes=["SE"]):
    
    data = {}
    for splicetype in splicetypes:
        data[splicetype] = {}
    for sample in samples:
        sampleFilename, sampleLabel = sample
        sampleData = pickle.load(open(sampleFilename))
        for i, geneItem in enumerate(sampleData):
            gene = geneItem["descriptor"]
            if not gene in data[splicetype]:
                data[splicetype][gene] = {}
            if not splicetype in geneItem:
                continue
            for loc in geneItem[splicetype]:
                if not loc in data[splicetype][gene]:
                    data[splicetype][gene][loc] = {}
                data[splicetype][gene][loc][sampleLabel] = geneItem[splicetype][loc]

    return data


def main(options):
    samples = options.samples
    spliceData = mergeSamples(samples, splicetypes=options.splicetype)
    #import code
    #code.interact(local=locals())
    pval_cutoff = options.pval

    if options.species is not None:
        annotation = retrieve_splicing(options.species)
    
    if "SE" in spliceData:
        if len(options.samples) == 2:
            s1_label = samples[0][1]
            s2_label = samples[1][1]            
            SEoutfile = s1_label + ".vs." + s2_label + ".SEs_comparison"
            SEout = open(SEoutfile, 'w')
            if options.species is not None:
                SEBEDfile = s1_label + ".vs." + s2_label + ".SEs_comparison.BED"
                SEbed = open(SEBEDfile, 'w')

            s1_inlabel = "_".join([s1_label, "IN"])
            s1_exlabel = "_".join([s1_label, "EX"])

            s2_inlabel = "_".join([s2_label, "IN"])
            s2_exlabel = "_".join([s2_label, "EX"])            
            
            header= "\t".join(["Gene", "Exonloc", "p-value", "Test", "Testdetails", "significant?", "direction", s1_inlabel, s1_exlabel, s2_inlabel, s2_exlabel]) + "\n"

                
            SEout.write(header)
            for gene in spliceData["SE"]:
                for loc in spliceData["SE"][gene]:
                    sample1_IN = spliceData["SE"][gene][loc][samples[0][1]]["IN"]
                    sample1_EX = spliceData["SE"][gene][loc][samples[0][1]]["EX"]
                    sample2_IN = spliceData["SE"][gene][loc][samples[1][1]]["IN"]
                    sample2_EX = spliceData["SE"][gene][loc][samples[1][1]]["EX"]
                    sampledata = np.array([[sample1_IN, sample1_EX], [sample2_IN, sample2_EX]])

                    s1_ratio= (sample1_IN+1)/float(sample1_EX+1)
                    s2_ratio= (sample2_IN+1)/float(sample2_EX+1)
                    direction = np.sign(s2_ratio - s1_ratio)

                    if np.any(sampledata < 5):
                        test="fisher_exact"
                        odds, p = scipy.stats.fisher_exact(sampledata)
                        testdetails= "%e" %(odds)
                        
                    else:
                        test = "chi"
                        chi2, p, dof, exp = scipy.stats.chi2_contingency(sampledata, correction=True)

                        testdetails= "%e" %(chi2)
                    issig= "no"

                    if p < pval_cutoff:
                        issig = "yes"
                        if options.species is not None:
                            chr, start, stop, name, score, strand = annotation[gene]["SE"][loc]["bedTrack"].split("\t")
                            if direction <0:
                                color = "0,255,0"
                            else:
                                color = "255,0,0"
                            sci_pval = "%E" %(p) #scientific notation
                            bedline = "\t".join([chr, start, stop, gene, sci_pval, strand, start, stop, color]) + "\n"
                            SEbed.write(bedline)
                    line = "\t".join(map(str, [gene, loc, p, test, testdetails, issig, direction, sample1_IN, sample1_EX, sample2_IN, sample2_EX]))
                    SEout.write(line + "\n")

    if "MXE" in spliceData:
        if len(options.samples) == 2:
            s1_label = samples[0][1]
            s2_label = samples[1][1]            
            MXEoutfile = s1_label + ".vs." + s2_label + ".MXEs_comparison"
            MXEout = open(MXEoutfile, 'w')
            if options.species is not None:
                MXEBEDfile = s1_label + ".vs." + s2_label + ".MXEs_comparison.BED"
                MXEbed = open(MXEBEDfile, 'w')

            s1_inlabel = "_".join([s1_label, "A"])
            s1_exlabel = "_".join([s1_label, "B"])

            s2_inlabel = "_".join([s2_label, "A"])
            s2_exlabel = "_".join([s2_label, "B"])            
            
            header= "\t".join(["Gene", "Exonloc", "p-value", "Test", "Testdetails", "significant?", "direction", s1_inlabel, s1_exlabel, s2_inlabel, s2_exlabel]) + "\n"

                
            MXEout.write(header)
            for gene in spliceData["MXE"]:
                for loc in spliceData["MXE"][gene]:
                    sample1_IN = spliceData["MXE"][gene][loc][samples[0][1]]["A"]
                    sample1_EX = spliceData["MXE"][gene][loc][samples[0][1]]["B"]
                    sample2_IN = spliceData["MXE"][gene][loc][samples[1][1]]["A"]
                    sample2_EX = spliceData["MXE"][gene][loc][samples[1][1]]["B"]
                    sampledata = np.array([[sample1_IN, sample1_EX], [sample2_IN, sample2_EX]])

                    s1_ratio= (sample1_IN+1)/float(sample1_EX+1)
                    s2_ratio= (sample2_IN+1)/float(sample2_EX+1)
                    direction = np.sign(s2_ratio - s1_ratio)

                    if np.any(sampledata < 5):
                        test="fisher_exact"
                        odds, p = scipy.stats.fisher_exact(sampledata)
                        testdetails= "%e" %(odds)
                        
                    else:
                        test = "chi"
                        chi2, p, dof, exp = scipy.stats.chi2_contingency(sampledata, correction=True)

                        testdetails= "%e" %(chi2)
                    issig= "no"

                    if p < pval_cutoff:
                        issig = "yes"
                        if options.species is not None:
                            try:
                                chr, start, stop, name, score, strand = annotation[gene]["MXE"][loc]["bedTrack"].split("\t")
                                if direction <0:
                                    color = "0,255,0"
                                else:
                                    color = "255,0,0"
                                sci_pval = "%E" %(p) #scientific notation
                                bedline = "\t".join([chr, start, stop, gene, sci_pval, strand, start, stop, color]) + "\n"
                                MXEbed.write(bedline)
                            except:
                                print loc + " failed, no bed line created"
                                pass
                    line = "\t".join(map(str, [gene, loc, p, test, testdetails, issig, direction, sample1_IN, sample1_EX, sample2_IN, sample2_EX]))
                    MXEout.write(line + "\n")                    
                    
                    
                
                
    






            

        


if __name__ == "__main__":
    parser = OptionParser()

    
    parser.add_option("--sample", nargs=2, action="append", dest="samples", help="Two values: --sample filename label")

    parser.add_option("--pvalue", dest="pval", default=0.05, help="p-value cutoff for chi2 or fisher exact")

    parser.add_option("--splice_type", dest="splicetype", default=["SE"], action="append")
    parser.add_option("--species", dest="species", default=None)


    options, args = parser.parse_args()
    main(options)





    

    
