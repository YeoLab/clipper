#!/bin/bash
#converts all my gtf files to bed files, forgot why I wanted this, we'll keep it here though.
cat hg19.AS.STRUCTURE.COMPILED.gff | awk 'BEGIN {OFS="\t"} {print $1, $4, $5, substr($9, 9, 17), 0, $7}' > hg19.AS.STRUCTURE.COMPILED.bed
cat mm9.AS.STRUCTURE.COMPILED.gff | awk 'BEGIN {OFS="\t"} {print $1, $4, $5, $2, 0, $7}' > mm9.AS.STRUCTURE.COMPILED.bed
cat mm10.AS.STRUCTURE.COMPILED.gff | awk 'BEGIN {OFS="\t"} {print $1, $4, $5, $2, 0, $7}' > mm10.AS.STRUCTURE.COMPILED.bed
cat ce10.AS.STRUCTURE.COMPILED.gff | awk 'BEGIN {OFS="\t"} {print $1, $4, $5, $2, 0, $7}' > ce10.AS.STRUCTURE.COMPILED.bed
