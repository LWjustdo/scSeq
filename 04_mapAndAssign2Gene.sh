#!/bin/bash

if [ $# != 3 ]
then 
    echo ""
    echo "           bash $0 <R2_fastq.gz> <prefix> <outDir>"
    echo ""
    echo ""
    exit 0
fi

R2=$1
prefix=$2
outDir=$3
log=$outDir/$prefix.time.log

if [ ! -d $outDir ]
then 
    mkdir -p $outDir
fi
cd $outDir
echo "****************************" >> $log
echo "Command: $0 $* " >> $log
echo "****************************" >> $log
echo "Job begin: `date` " >> $log

##########  mapping #########################
STAR --runMode alignReads --runThreadN 12 \
--genomeDir /home/longzhao/scSeq/star_index_hg38 \
--readFilesIn $R2 \
--readFilesCommand zcat \
--outFileNamePrefix $prefix. \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate 

wait
########## featureCounts ######################
/data/biosoft/subread-1.6.3-Linux-x86_64/bin/featureCounts -a /data/DataBase/GENECODE/GRCh38/gencode.v29.annotation.gtf -R BAM -g gene_name -o $prefix.featureCounts  ./$prefix.Aligned.sortedByCoord.out.bam -T 8 &>$prefix.featureCounts.log

samtools sort $prefix.Aligned.sortedByCoord.out.bam.featureCounts.bam $prefix.assigned_sorted

samtools index $prefix.assigned_sorted.bam

rm $prefix.Aligned.sortedByCoord.out.bam.featureCounts.bam

umi_tools count --per-gene --gene-tag=XT  --per-cell --wide-format-cell-counts -I $prefix.assigned_sorted.bam -S $prefix.counts.tsv.gz

wait

python /home/longzhao/scSeq/src/scSeq.py grep_STAR_map_info $prefix $prefix.Log.final.out

python /home/longzhao/scSeq/src/scSeq.py map_tsv2_all  $prefix $prefix.counts.tsv.gz

echo "Job end: `date` " >> $log




