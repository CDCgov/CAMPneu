#!/bin/bash

stats1x=$1
stats10x=$2
depth=$3
prefix=$4

printf "%s\t%s\t%s\t%s\t%s\n" "sample" "breadth_coverage>10X(%)" "breadth_coverage>0X(%)" "avg_depth" "avg_depth/breadth_coverage>0X" > ${prefix}_cov_depth.tsv

breadth1x=$(grep "percentage of target genome with coverage > 0 (%):" $stats1x | awk 'BEGIN{ OFS="\t" }{ print $11 }') 
breadth10x=$(grep "percentage of target genome with coverage > 10 (%):" $stats10x | awk 'BEGIN{ OFS="\t" }{ print $11 }') 

avg_depth=$(awk '{ cov_sum+=$3 }END{ if( cov_sum > 0 ){ print cov_sum/NR }else{ print 0 } }' $depth)

normalized=$(awk -vB=$breadth1x -vD=$avg_depth 'BEGIN{ b=B/100; print D/b }')

printf "%s\t%s\t%s\t%s\t%s\n" "$prefix" "$breadth10x" "$breadth1x" "$avg_depth" "$normalized" >> ${prefix}_cov_depth.tsv