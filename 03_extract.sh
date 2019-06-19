#!/bin/bash

if [ $# != 5 ]
then
    echo ""
    echo "        bash $0 <R1_fastq.gz> <R2_fastq.gz> <whitelist> <prefix> <Outdir>"
    echo ""
    echo ""
    exit 0
fi

R1=$1
R2=$2
whitelist=$3
prefix=$4
outDir=$5
log=$outDir/$prefix.time.log
if [ ! -d $outDir ]
then 
    mkdir -p $outDir
fi
cd $outDir
echo "************************" >> $log
echo "Command: bash $0 $*" >> $log
echo "************************" >> $log
echo "Job begin: `date`" >> $log

umi_tools extract \
--extract-method=regex \
--bc-pattern="(?P<cell_1>.{8,11}){s<=2}(?P<discard_1>GTGATTGCTTGTGAC){s<=2}(?P<cell_2>.{8}){s<=2}(?P<umi_1>.{6})T{5}.*" \
--filter-cell-barcode \
--error-correct-cell \
--whitelist=$whitelist \
--stdin $R1 --stdout $prefix\_extract_R1.fastq.gz \
--read2-in $R2 --read2-out $prefix\_extract_R2.fastq.gz \
-L $prefix.extract.log
wait

python /home/longzhao/scSeq/src/scSeq.py extract_head_seq $prefix 100 $prefix\_extract_R1.fastq.gz $prefix\_extract_R2.fastq.gz

fastp -i $prefix\_head_R2.fastq.gz -o $prefix\_head_R2.clean.fastq.gz \
      -j $prefix\_R2_head.json -h $prefix\_R2_head.html -R $prefix -A --thread 4

python /home/longzhao/scSeq/src/scSeq.py get_single_cell_read_count $prefix $prefix\_head_R2.clean.fastq.gz

python /home/longzhao/scSeq/src/scSeq.py summary4fastp $prefix\_R2_head.json

echo "Job end: `date`" >> $log


