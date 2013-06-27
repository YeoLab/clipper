#!/bin/sh
for i in `ls *sorted`
do
 perl /nas3/lovci/projects/scripts/add_score.pl $i
done
