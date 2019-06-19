#!/bin/bash

if [ $# != 4 ]
then
    echo ""
    echo "             bash $0 <R1_fastq.gz> <Prefix> <Cell_number> <Out_dir>"
    echo ""
    echo ""
    exit 0
fi

R1=$1
prefix=$2
cellNumber=$3
outDir=$4
if [ ! -d $outDir ]
then
    mkdir -p $outDir
fi
cd $outDir
log=$outDir/$prefix.time.log
echo "*****************************" >> $log
echo "Command: sh $0 $*" >> $log
echo "*****************************" >> $log
echo "Job begin: `date` " >> $log

umi_tools whitelist \
--stdin $R1 \
--extract-method=regex \
--bc-pattern="(?P<cell_1>.{8,11}){s<=2}(?P<discard_1>GTGATTGCTTGTGAC){s<=2}(?P<cell_2>.{8}){s<=2}(?P<umi_1>.{6})T{5}.*" \
--plot-prefix $prefix \
--set-cell-number=$cellNumber \
--log2stderr > $prefix.whitelist.tsv
wait
echo "Job end: `date` " >> $log
