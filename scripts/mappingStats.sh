#!/bin/bash
# Args 
# $1 name of sam-file
# $2 P-val hreshold 
# $3 human readable type 

UNIQUE=$(cat $1 | /home/frellsen/share/ceprano-tools/samfilter -s $2 -A | wc -l)
MAPPED=$(cat $1 | /home/frellsen/share/ceprano-tools/samfilter -A | wc -l)
ALL=$(cat $1 | grep -vc "^@")

PCTUNMAPPED=$(echo "scale=4; 1 - $MAPPED/$ALL;" | bc -l)
PCTMAPPED=$(echo "scale=4; $MAPPED/$ALL;" | bc -l)
PCTUNIQUE=$(echo "scale=4; $UNIQUE/$ALL;" | bc -l)
PCTMULTIPLE=$(echo "scale=4; $PCTMAPPED - $PCTUNIQUE;" | bc -l)

if [ $3 = 'y' ]; then
    echo "all" $ALL
    echo "probability_threshold",$2
    echo "Mapped          ", 0$PCTMAPPED, $MAPPED
    echo "Unmapped        ", 0$PCTUNMAPPED, $(( $ALL - $MAPPED )) 
    echo "Confidently_mapped ", 0$PCTUNIQUE, $UNIQUE
    echo "Multiple_mapped ", 0$PCTMULTIPLE, $(( $MAPPED - $UNIQUE ))
else
MULTIPLE=$(( $MAPPED - $UNIQUE ))
UNMAPPED=$(( $ALL - $MAPPED ))
echo $UNIQUE $MULTIPLE $ALL $MAPPED $UNMAPPED
fi