__author__ = 'gpratt'

import pysam

from optparse import OptionParser
import os
import pybedtools
import scipy.stats
import numpy as np

def format_line(interval):
    return "{}:{}-{}:{}:{}".format(interval.chrom, interval.start, interval.stop, interval.strand, interval.score)


def chisq_value(ip_count, ip_total, input_count, input_total):
    try:
        g, p, dof, expctd = scipy.stats.chi2_contingency([[ip_count, ip_total], [input_count, input_total]])
    except Exception:
        print ip_count, ip_total, input_count, input_total
        raise

    return g, p


def fisher_exact(ip_count, ip_total, input_count, input_total):
    try:
        oddsratio, p = scipy.stats.fisher_exact([[ip_count, ip_total], [input_count, input_total]])
    except Exception:
        print ip_count, ip_total, input_count, input_total
        raise

    return "F", p


def fisher_or_chisq(ip_count, ip_total, input_count, input_total):

    total = ip_count + ip_total + input_count + input_total

    expected_ip_count = (ip_count + input_count) * (ip_count + ip_total) / total
    expected_ip_total = (ip_total + input_total) * (ip_count + ip_total) / total
    expected_input_count = (ip_count+input_count) * (input_count + input_total) / total
    expected_input_total = (ip_total + input_total) * (input_count + input_total) / total

    direction = "enriched"
    if ip_count < expected_ip_count:
        direction = "depleted"
        return 1, "DEPL", "N", direction
    if any(i < 5 for i in [ip_count, ip_total, input_count, input_total,
                           expected_ip_count, expected_ip_total, expected_input_count, expected_input_total]):
        val, p = fisher_exact(ip_count, ip_total, input_count, input_total)
        analysis_type = "F"
    else:
        val, p = chisq_value(ip_count, ip_total, input_count, input_total)
        analysis_type = "C"
    return p, val, analysis_type, direction


def adjusted_positions(read):
    positions = enumerate(read.positions)
    for t in read.cigartuples:
        value, times = t

        #value 3 is splice junction value 2 is deletion in read
        if value == 3 or value == 1 or value == 4:
            continue

        if value == 2:  #In the case of a deletion yeild from the previous current position, this could be buggy in the case of deletions occuring after anything but a match
            for x in xrange(times):
                cur_position += 1
                yield cur_position

        for x in xrange(times):
            cur_position = positions.next()[1]
            yield cur_position


def stranded_reads(bam, interval, reads_found=None):
    reads = list(bam.fetch(interval.chrom, interval.start, interval.stop))
    #this is a possibly too smart way of getting only reads on the correct strand, using xor, I tested it an it works
    reads = [read for read in reads if (interval.strand == "+") ^ read.is_reverse]
    final_reads = []
    for read in reads:
        for position in adjusted_positions(read):
            if interval.start <= position < interval.stop:
                #I was trying to make this more elegant, handing both a read hash and not one, but I don't know if it worked
                if reads_found is not None:
                    if read.qname in reads_found:
                        break
                    reads_found.add(read.qname)
                final_reads.append(read)
                break

    return final_reads


def calculate_entropy(ip_count, input_count, ip_total, input_total):
    p_ip = float(ip_count) / ip_total
    p_input = float(input_count) / input_total

    try:
        entropy = p_ip * np.log2(p_ip / p_input)
    except FloatingPointError:
        print ip_count, input_count, ip_total, input_total, p_ip, p_input

    return entropy


def calculate_log2fc(ip_count, input_count, ip_mapped, input_mapped):

    #rpr is reads per read, ignoring the million mapped number
    ip_rpm = (ip_count / ip_mapped) * 1000000
    input_rpm = (input_count / input_mapped) * 1000000
    log2fc = np.log2(ip_rpm / input_rpm)

    return log2fc


def write_lines(interval, ip_count, input_count, ip_mapped, input_mapped, out_file, out_full_file):

    #Sometimes "peaks" dont have reads in them, off by one bug, this is the easiest way to fix
    ip_count = 1.0 if ip_count == 0 else ip_count
    ip_count = float(ip_count)

    #add pesudocount of 1 to the input
    input_count = float(input_count) + 1

    line = format_line(interval)
    ip_total = float(ip_mapped - ip_count)
    input_total = float(input_mapped - input_count)

    p, val, analysis_type, direction = fisher_or_chisq(ip_count, ip_total, input_count, input_total)

    log10_p = abs(400.0 if p == 0 else np.log10(p) * -1)
    log2fc = calculate_log2fc(ip_count, input_count, ip_total, input_mapped)
    entropy = calculate_entropy(ip_count, input_count, ip_mapped, input_mapped)

    out_file.write("\t".join(map(str, [interval.chrom, interval.start, interval.stop,
                                           log10_p, log2fc, interval.strand, entropy])) + "\n")

    out_full_file.write("\t".join(map(str, [interval.chrom, interval.start, interval.stop,
                                       line, int(ip_count), int(input_count), p, val, analysis_type, direction,
                                       log10_p, log2fc, entropy])) + "\n")


def input_norm_peaks(peaks, ip_bam, input_bam, out, out_full):
    with open(out_full, 'w') as out_full_file:
        with open(out, 'w') as out_file:
            for interval in peaks:
                ip_count = len(stranded_reads(ip_bam, interval))
                input_count = len(stranded_reads(input_bam, interval))
                write_lines(interval, ip_count, input_count, ip_bam.mapped, input_bam.mapped, out_file, out_full_file)


def em_peaks(peaks, ip_bam, input_bam, out, out_full):
    bedtool = pybedtools.BedTool(peaks)
    bedtool = bedtool.to_dataframe(names=["chrom", "start", "stop", "name", "ip_count", "input_count", "p_val",
                                          "chisq_stat", "test_type", "enriched", "log10_pval", "log2_fc", "entropy"])
    bedtool = bedtool.sort_values("entropy", ascending=False)

    ip_reads_found = set([])
    input_reads_found = set([])
    with open(out, 'w') as out_file:
        with open(out_full, 'w') as out_full_file:
            for name, row in bedtool.iterrows():
                chrom, loc, strand, score = row.values[3].split(":")
                interval = pybedtools.create_interval_from_list([row[0], row[1], row[2], 'foo', score, strand])

                ip_count = len(stranded_reads(ip_bam, interval, ip_reads_found))
                input_count = len(stranded_reads(input_bam, interval, input_reads_found))
                write_lines(interval, ip_count, input_count, ip_bam.mapped, input_bam.mapped, out_file, out_full_file)


def run(peaks, ip_bam, input_bam, out, out_full, just_norm=False):
    out_prefix, out_ext = os.path.splitext(out)
    out_full_prefix, out_full_ext = os.path.splitext(out_full)

    out_version = "{}.v{}{}".format(out_prefix, 0, out_ext)
    out_full_version = "{}.v{}{}".format(out_full_prefix, 0, out_full_ext)
    input_norm_peaks(peaks, ip_bam, input_bam, out_version, out_full_version)

    if just_norm:
        return

    #explicitly define the number of round of input norm to do
    for x in range(1, 4):
        prev_out_full_version = out_full_version
        out_version = "{}.v{}{}".format(out_prefix, x, out_ext)
        out_full_version = "{}.v{}{}".format(out_full_prefix, x, out_full_ext)

        em_peaks(prev_out_full_version, ip_bam, input_bam, out_version, out_full_version)


def call_main():

    description = """CLIPper Input Norm. Michael Lovci, Gabriel Pratt, Eric Van Nostrand 2017.
                     Input Norm for CLIPper.  For results similar to ENCODE peaks on the DCC website
                     (but not identitacal due to numerical percision issues between python and the initial perl implementation
                     use v0 outputs).  The iterative assignment becomes not useful after ~4 rounds.
                     Refer to: https://github.com/YeoLab/clipper/wiki for instructions.
                     Questions should be directed to gpratt@ucsd.edu."""

    parser = OptionParser(description=description)

    parser.add_option("--peaks", "-p", dest="peaks", help="A peak file to input normalize", type="string")
    parser.add_option("--ip_bam", dest="ip_bam", help="a bam file of the IP reads to input normalze against")
    parser.add_option("--input_bam", dest="input_bam", help="a bam file of the input reads to input normalze against")
    parser.add_option("--out", "-o", dest="out", default="norm_clusters.bed",
                      help="a bed file output, default:%default")
    parser.add_option("--out_full", "-f", dest="out_full", default="norm_clusters.full",
                      help="A full file output (used for specific analysies), default:%default")

    (options, args) = parser.parse_args()
    run(pybedtools.BedTool(options.peaks), pysam.Samfile(options.ip_bam), pysam.Samfile(options.input_bam), options.out,
        options.out_full)

if __name__ == "__main__":
    call_main()
