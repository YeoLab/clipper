#import numpy

def readsToWiggle_pysam(reads, int tx_start, int tx_end, keepstrand, usePos, bint fracional_input):
    """
    
    converts pysam to a wiggle vector and some other stuff.
    
    input: (bamfile.fetch obj), tx_start, tx_end, strand, readPos, trim
    bamfile.fetch obj is from pysam.Samfile().fetch    
    tx_start is the genome coordinate of the start position of the window you care about
    tx_stop is genome stop
    strand is + or - to indicate which strand to create a wiggle from
    usePos is the position you'll be calling the read cover from one of: ["center", "start", "end"]
    fractional: boolean output fractional results instead of integers 
    
    output: wiggle array, jxns, positional counts (cover), read lengths, read locations
    
    """
    cdef int cur_pos
    cdef int read_start
    cdef int read_stop
    cdef int next_pos
    cdef float increment_value
    #cdef vector[int] vect

    gene_size = tx_end - tx_start + 1
    lengths = []

    #Something about numpy makes this break, the error is SystemError: Objects/listobject.c:169: bad argument to internal function
    #Probably better to use numpy arrays, but I don't know how much better
    
    wiggle = [0] * gene_size #numpy.zeros(shape = gene_size, dtype = 'i')
    pos_counts = [0] * gene_size #numpy.zeros(shape = gene_size, dtype = 'i')
    explicit_locations = [0] * gene_size #numpy.zeros(shape = gene_size, dtype = 'O')

    for x in xrange(len(explicit_locations)):
        explicit_locations[x] = set([])

    all_reads = set([])
    junctions = {}
    
    for read in reads:
        if read.is_reverse and keepstrand == "+":
            continue
        elif not read.is_reverse and keepstrand == "-":
            continue
            
        read_start = read.positions[0]
        read_stop = read.positions[-1]
        
        if read_start < tx_start or read_stop > tx_end:
            continue
            
        read_len = len(read.positions)
        lengths.append(read_len)
        
        if usePos == "center":
            pos_counts[(((read_stop + read_start) / 2) - tx_start)] += 1
        elif usePos == "start":
            if keepstrand == "+":
                pos_counts[read_start - tx_start] += 1
            else: 
                pos_counts[read_stop - tx_start] += 1
        elif usePos == "end":
             if keepstrand == "-":
                pos_counts[read_start - tx_start] += 1
             else: 
                pos_counts[read_stop - tx_start] += 1
        
        all_reads.add((read_start, read_stop))
        
        increment_value = (1.0 / read_len) if fracional_input else 1.0

        cigops = list(get_full_length_cigar(read))

        for cur_pos, next_pos, cigop in zip(read.positions, read.positions[1:], cigops):
            #if cur is not next to the next position than its a junction
            if cur_pos + 1 != next_pos:
                junction = (cur_pos, next_pos)
                if junction not in junctions:
                    junctions[junction] = 0
                junctions[junction] += 1
                
            wiggle[cur_pos - tx_start] += increment_value
            if cigop == 0: #Exact matches only, doing this because it duplicates HTSeq behavior
                explicit_locations[cur_pos - tx_start].add(read)

        #needed to get last read counted
        wiggle[read.positions[-1] - tx_start] += increment_value
        if cigops[-1] == 0:
            explicit_locations[read.positions[-1] - tx_start].add(read)

    return wiggle, junctions, pos_counts, lengths, all_reads, explicit_locations

def get_full_length_cigar(read):
    for t in read.cigartuples:
        value, times = t
        for x in xrange(times):
            yield value