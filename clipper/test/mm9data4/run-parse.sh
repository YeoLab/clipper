#!/bin/sh

for chr in `seq 1 19` "X" "Y"
do
    perl /nas3/ppliu/projects/ASparser/parser_to_gff.pl mm9.tx.$chr.AS.STRUCTURE > mm9.tx.$chr.AS.STRUCTURE.SE.gff
done


rm -rf mm9.AS.STRUCTURE.gff
for chr in `seq 1 19` "X" "Y"
do
    cat mm9.tx.$chr.AS.STRUCTURE.SE.gff >> mm9.AS.STRUCTURE.SE.gff
    rm -rf mm9.tx.$chr.AS.STRUCTURE.SE.gff
done

