import re
import time
import pysam
from seqTools import *
import sys
import os
from deap import dtm
from optparse import OptionParser

host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
if "optiputer" in host or "compute" in host:
	basedir = "/nas/nas0"
elif "tcc" in host or "triton" in host:
	basedir = "/projects"
else:
	print "Where am I?"
	raise Exception

parser = OptionParser()
parser.add_option("-s", "--species", dest="species")
parser.add_option("-b", "--bam_file", dest="bam_path")
parser.add_option("-f", "--flip", dest="flip")
parser.add_option("-o", "--out_file", dest="out_file")
parser.add_option("-r", "--raw_count_out", dest="raw_count_file")




#bam_path = sys.argv[1]
#flip = sys.argv[2]
#bam_file = pysam.Samfile(bam_path, 'rb')
gene_info = dict()
regions_info = dict()

#out_file = open(sys.argv[3], 'a');
#raw_count_file = open(sys.argv[4], 'a');

def get_wig(keys):

	if (keys not in gene_info):
		print keys
		print "not here"

	try:
		subset_reads = bam_file.fetch(reference = gene_info[keys]['chr'], start = int(gene_info[keys]["start"]), end = int(gene_info[keys]["stop"]))
		
		keep_strand = gene_info[keys]["strand"]
		if (str(flip) == "flip"):
			if (str(keep_strand) == '-'):
				keep_strand = '+'
			elif (str(keep_strand) == '+'):
				keep_strand = '-'

		wig, jxns, nr_counts, read_lengths, reads = readsToWiggle_pysam(subset_reads, int(gene_info[keys]["start"]), int(gene_info[keys]["stop"]), keepstrand = keep_strand)

	except:
		print "part1"
	lines=[]
	region_sum=0
	while (len(regions_info[keys]) >=2 ):

		try:
			start = int(regions_info[keys].pop(0)) - int(gene_info[keys]["start"] )
			stop = int(regions_info[keys].pop(0)) - int(gene_info[keys]["start"] )
		except:
			print "part2"

		sum = 0
		#sum = sum(wig[start:stop])
		for x in range(start, stop):
			sum+=wig[x]	
		
		region_sum+=sum

		read_count=0
		for y in nr_counts:
			read_count+=int(y)

		gene_info[keys]["raw_count"] = read_count

		#lines.append( keys+"\t"+str(int(start)+int(gene_info[keys]["start"] ))+"\t"+str(int(stop)+int(gene_info[keys]["start"]) )+"\t"+str(sum))

		lines.append( keys+"\t"+str(int(start)+int(gene_info[keys]["start"] ))+"\t"+str(int(stop)+int(gene_info[keys]["start"]) )+"\t"+str(read_count)+"\t"+str(sum))

	for i, line in enumerate(lines):
		lines[i] = line+"\t"+str(region_sum)
	
	return lines

def count_to_regions(species):

	if species == "hg19":
		chrs = map(str,range(1,23)) #1-22
		chrs.append("X")
		chrs.append("Y")
	elif species == "mm9":
		chrs = map(str, range(1,19))
		chrs.append("X")
		chrs.append("Y")
	elif species == "ce6":
		chrs = ("I", "II", "III", "IV", "V", "X")

	for chr in chrs:
		regions_file = "/nas3/ppliu/genic_regions/"+species+"/genic_regions_"+species+".chr"+chr
		f = open(regions_file, 'r')
		for line in f:

			line = line.strip()
			chromosome, start, stop, strand, ensembl_id, frea_annot = line.strip().split("\t");
	

			if ( int(strand) == 1):
				strand = '-'
			elif( int(strand) == 0):
				strand = '+'
		
			if (ensembl_id not in gene_info):
				regions_info[ensembl_id] = []

			regions_info[ensembl_id].append(start)
			regions_info[ensembl_id].append(stop)



			if (ensembl_id in gene_info):
				if (int(start) < int(gene_info[ensembl_id]["start"])):
					gene_info[ensembl_id]["start"] = start
				if (int(stop) > int(gene_info[ensembl_id]["stop"])):
					gene_info[ensembl_id]["stop"] = stop

			else:
				gene_info[ensembl_id]={}
				gene_info[ensembl_id]["chr"] = chromosome
				gene_info[ensembl_id]["start"] = start
				gene_info[ensembl_id]["stop"] = stop
				gene_info[ensembl_id]["strand"] = strand
				gene_info[ensembl_id]["frea"] = frea_annot

		
		f.close()


def main():
	


	print "Day, Date :", time.strftime("%A, %B %d, %Y", time.localtime())
	print "Time (12hr) :", time.strftime("%I:%M:%S %p", time.localtime())

	key_list=[]
	for keys in gene_info:
		key_list.append(keys)

	region_lines = dtm.map(get_wig, iter(key_list))
	for info in region_lines:
		out_file.write("\n".join(info))
		out_file.write("\n")
	
	for keys in gene_info:
		raw_count_file.write(gene_info[keys]["chr"]+"\t"+gene_info[keys]["start"]+"\t"+gene_info[keys]["stop"]+"\t"+gene_info[keys]["raw_count"]+"\n")

	print "Day, Date :", time.strftime("%A, %B %d, %Y", time.localtime())
	print "Time (12hr) :", time.strftime("%I:%M:%S %p", time.localtime())

if __name__ == "__main__":


	count_to_regions("hg19")
	dtm.start(main)
	
