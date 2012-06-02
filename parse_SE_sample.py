import sys
f = open(sys.argv[1])
f.readline() #header
for line in f.readlines():
    event, misoMean, ciLow, ciHigh, isoforms, counts, acounts = line.strip().split("\t")
    
    inclusion = -1
    exclusion = -1
    

    if not "," in misoMean:
        inclusion = float(misoMean)
        exclusion = 1-float(misoMean)
    print "\t".join(map(str, [event, inclusion, exclusion]))

