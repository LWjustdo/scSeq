#!/bin/bash

if [ $# != 4 ]
then
    echo ""
    echo "         bash $0 <Count_tsv> <cell_Reads_count.tsv> <prefix> <outDir>"
    echo ""
    echo ""
    echo "         Note:"
    echo "              1).The <Count_tsv> is output from 04_mapAndAssign2Gene.sh;"
    exit 0
fi

tsv=$1
readstsv=$2
prefix=$3
outDir=$4
log=$outDir/$prefix.time.log
if [ ! -d $outDir ]
then
    mkdir -p $outDir
fi
cd $outDir

echo "**********************************" >> $log
echo "Command: sh $0 $* " >> $log
echo "**********************************" >> $log
echo "Job begin: `date` " >> $log

Rscript /home/longzhao/scSeq/src/seurat_multiSample.R $outDir $prefix 100 200 0.1 TRUE $tsv $prefix
wait 
if [ ! -f $outDir/$prefix\_ClusterHeatmap.png ]
then 
    rm $outDir/*png $outDir/*txt $outDir/*rds $outDir/*csv $outDir/*tsv
    echo "" >> $log
    echo "Redo 2" >>  $log
    echo "Command: sh $0 $outDir $prefix 100 150 0.3 TRUE $tsv $prefix" >> $log
    Rscript /home/longzhao/scSeq/src/seurat_multiSample.R $outDir $prefix 100 150 0.3 TRUE $tsv $prefix
    wait
fi

if [ ! -f $outDir/$prefix\_ClusterHeatmap.png ]
then 
    rm $outDir/*png $outDir/*txt $outDir/*rds $outDir/*csv $outDir/*tsv
    echo "" >> $log
    echo "Redo 3" >>  $log
    echo "Command: sh $0 $outDir $prefix 30 80 0.3 FALSE $tsv $prefix" >> $log
    Rscript /home/longzhao/scSeq/src/seurat_multiSample.R $outDir $prefix 30 80 0.3 FALSE $tsv $prefix
    wait
fi

python /home/longzhao/scSeq/src/scSeq.py sequence_saturation $prefix $prefix\_metaData.csv $readstsv

python /home/longzhao/scSeq/src/scSeq.py grep_seurat_cell_filter $prefix $prefix.seuratAnalysis4Muti.txt

echo "Job end: `date` " >> $log
