#!/bin/bash

if [ $# != 5 ]
then
    echo ""
    echo "        bash $0 <R1_fastq.gz> <R2_fastq.gz> <prefix> <cell_Number> <Outdir>"
    echo ""
    echo ""
    exit 0
fi

R1=$1
R2=$2
prefix=$3
cellNum=$4
outDir=$5

log=$outDir/$prefix.time.log

if [ ! -d $outDir ]
then 
    mkdir -p $outDir
fi
cd $outDir

echo "=============================================" >> $log
echo "Command: $0 $*" >> $log
echo "=============================================" >> $log

echo "[`date`] 01_QC4scR1 begin......" >> $log
bash /home/longzhao/scSeq/01_QC4scR1.sh $R1 $prefix $outDir/01_QC

echo "[`date`] 02_whiteList begin......" >> $log
bash /home/longzhao/scSeq/02_whitelist.sh $R1 $prefix $cellNum $outDir/02_whitelist

echo "[`date`] 03_extract begin......" >> $log
bash /home/longzhao/scSeq/03_extract.sh $R1 $R2 $outDir/02_whitelist/$prefix.whitelist.tsv $prefix $outDir/03_extract

echo "[`date`] 04_mapAndAssign begin......" >> $log
bash /home/longzhao/scSeq/04_mapAndAssign2Gene.sh $outDir/03_extract/$prefix\_head_R2.clean.fastq.gz $prefix $outDir/04_mapAndAssign

echo "[`date`] 05_scSeqAnalysis begin......" >> $log
bash /home/longzhao/scSeq/05_seurat_oneSample.sh $outDir/04_mapAndAssign/$prefix\_all_counts.tsv.gz $outDir/03_extract/$prefix\_cellReadsCounts.tsv $prefix $outDir/05_scSeqAnalysis

echo "[`date`] 06_GenerateReport begin......" >> $log
bash /home/longzhao/scSeq/06_report_oneSample.sh $outDir $outDir/$prefix\_Result $prefix


echo "[`date`] Task finished." >> $log




