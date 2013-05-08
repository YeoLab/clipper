def readsToWiggle_pysam(reads, int tx_start, int tx_end, keepstrand, usePos, bint fracional_input):
    """
    
    converts pysam to a wiggle vector and some other stuff.
    
    input: ((bamfile.fetch obj), tx_start, tx_end, strand, readPos, trim
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
    
    lengths = []
    wiggle = [0] * (tx_end - tx_start + 1)
    pos_counts = [0] * (tx_end - tx_start + 1)
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
        
        for cur_pos, next_pos in zip(read.positions, read.positions[1:]):
            #if cur is not next to the next position than its a junction
            if cur_pos + 1 != next_pos:
                junction = (cur_pos, next_pos)
                if junction not in junctions:
                    junctions[junction] = 0
                junctions[junction] += 1
                
            wiggle[cur_pos - tx_start] += increment_value
        
        #needed to get last read counted
        wiggle[read.positions[-1] - tx_start] += increment_value
                
    return wiggle, junctions, pos_counts, lengths, all_reads