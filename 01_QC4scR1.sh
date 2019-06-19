#!/bin/bash

if [ $# != 3 ]
then 
    echo ""
    echo "        Usage:    <R1_fastq.gz> <sampleName> <outDir>"
    echo ""
    echo "        RefInfo:  https://github.com/OpenGene/fastp"
    echo ""
    exit 1
fi

R1=$1
sample=$2
outDir=$3
Log=$outDir/$sample.time.log

if [ ! -d $outDir ]
then
    mkdir -p $outDir
fi
cd $outDir
echo "*****************************" >> $Log
echo "Command: sh $0 $*" >> $Log
echo "*****************************" >> $Log
echo "Job begin: `date`" >> $Log

fastp -i $R1 -L --length_limit 40 -j $sample.qc4scR1.json -h $sample.scSeqR1.html -R $sample -A --thread 4 
wait

python /home/longzhao/scSeq/src/scSeq.py summary4fastp $sample.qc4scR1.json

echo "Job end: `date`" >> $Log
