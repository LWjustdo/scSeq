#!/bin/bash

if [ $# != 3 ]
then
    echo ""
    echo "        bash $0 <In_Dir> <Out_Dir> <prefix>"
    echo ""
    echo ""
    exit 0
fi

inDir=$1
outDir=$2
prefix=$3
log=$outDir/$prefix.report.log

if [ ! -d $outDir ]
then
    mkdir -p $outDir
fi
cd $outDir
mkdir -p $outDir/01_QC $outDir/02_whitelist $outDir/03_extract $outDir/04_mapAndAssign $outDir/05_scSeqAnalysis $outDir/img
echo "=============================================" >> $log
echo "Command: $0 $*" >> $log
echo "=============================================" >> $log
echo "[`date`] Job begin...." >> $log

#========= 01_QC ===================
#mv $inDir/01_QC/$prefix.R1.html $outDir/01_QC/
cp $inDir/01_QC/$prefix\_QCsummary.xls $outDir/01_QC/
#======== 02_whitelist =============
cp $inDir/02_whitelist/$prefix\_*.png $outDir/02_whitelist/
cp $inDir/02_whitelist/$prefix.whitelist.tsv  $outDir/02_whitelist/
#======== 03_extract ===============
#mv $inDir/03_extract/$prefix\_extract_R*.fastq.gz $outDir/03_extract/
cp $inDir/03_extract/$prefix\_cellReadsCounts.tsv $outDir/03_extract/
cp $inDir/03_extract/$prefix\_cellReadsInfo.tsv $outDir/03_extract/
cp $inDir/03_extract/$prefix\_R2_head_QCsummary.xls $outDir/03_extract/
#======== 04_mapAndAssign ==========
#mv $inDir/04_mapAndAssign/$prefix.assigned_sorted.bam $outDir/04_mapAndAssign/
#mv $inDir/04_mapAndAssign/$prefix.assigned_sorted.bam.bai $outDir/04_mapAndAssign/
cp $inDir/04_mapAndAssign/$prefix\_mapstat.tsv $outDir/04_mapAndAssign/
#======== 05_scSeqAnalysis =========
cp $inDir/04_mapAndAssign/$prefix\_all_counts.tsv.gz $outDir/05_scSeqAnalysis/
cp $inDir/05_scSeqAnalysis/$prefix\_saturation.tsv $outDir/05_scSeqAnalysis/
cp $inDir/05_scSeqAnalysis/$prefix\_*.csv $outDir/05_scSeqAnalysis/
cp $inDir/05_scSeqAnalysis/$prefix\_cellFiltered.tsv $outDir/05_scSeqAnalysis/
cp $inDir/05_scSeqAnalysis/$prefix\_ClusterHeatmap.png $outDir/05_scSeqAnalysis/
cp $inDir/05_scSeqAnalysis/$prefix\_GenePlot.png $outDir/05_scSeqAnalysis/
cp $inDir/05_scSeqAnalysis/$prefix\_QCPlot.png $outDir/05_scSeqAnalysis/
cp $inDir/05_scSeqAnalysis/$prefix\_MitoPlot.png $outDir/05_scSeqAnalysis/
cp $inDir/05_scSeqAnalysis/$prefix\_TSNEPlot.png $outDir/05_scSeqAnalysis/
cp $inDir/05_scSeqAnalysis/$prefix\_VariableGenes.png $outDir/05_scSeqAnalysis/
#
wait
#generate result.json
python /home/longzhao/scSeq/src/scSeq.py generate_json $prefix /home/longzhao/scSeq/src/report/temp/result.json
#prepration log img
cp -r /home/longzhao/scSeq/src/report/img/* $outDir/img/
#generate report
python /home/longzhao/scSeq/src/scSeq.py generate_report $prefix /home/longzhao/scSeq/src/report/temp/scSeq_template.html $prefix\_result.json $outDir
#html2pdf
wkhtmltopdf --header-html /home/longzhao/scSeq/src/report/temp/header.html --footer-html /home/longzhao/scSeq/src/report/temp/footer.html --margin-top 2cm --margin-right 1cm --margin-bottom 1.8cm --margin-left 1cm --encoding utf-8 --header-spacing 12 --footer-spacing 4 $prefix\_scSeqReport.html $prefix\_scSeqReport.pdf
#
mv $prefix\_result.json ../

echo "[`date`] Job finished..." >> $log
mv $log ../








