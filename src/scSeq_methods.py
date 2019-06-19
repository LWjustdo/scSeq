# python3
# coding = utf-8
# Author: liang.wan, 2019-05-20
"""
tools for DragonDropinn single cell sequence pipeline
"""
import os
import json
import gzip
import csv
import re
import time
#from collections import

def fastqIterate(infile):
    '''iterate over contents of fastq file.'''

    def convert2string(b):
        if type(b) == str:
            return b
        else:
            return b.decode("utf-8")

    while 1:
        line1 = convert2string(infile.readline())
        if not line1:
            break
        if not line1.startswith('@'):
            raise Exception("parsing error: expected '@' in line %s" % line1)
        line2 = convert2string(infile.readline())
        line3 = convert2string(infile.readline())
        if not line3.startswith('+'):
            raise Exception("parsing error: expected '+' in line %s" % line3)
        line4 = convert2string(infile.readline())
        # incomplete entry
        if not line4:
            raise Exception("incomplete entry for %s" % line1)

        yield "@%s\n%s\n+\n%s" %(line1[1:-1],line2[:-1], line4[:-1])


def joinedFastqIterate(fastq_iterator1, fastq_iterator2, strict=True):
    '''This will return an iterator that returns tuples of fastq records.
    At each step it will confirm that the first field of the read name
    (before the first whitespace character) is identical between the
    two reads. The response if it is not depends on the value of
    :param:`strict`. If strict is true an error is returned. If strict
    is `False` the second file is advanced until a read that matches
    is found.

    This allows for protocols where read one contains cell barcodes, and these
    reads have been filtered and corrected before processing without regard
    to read2

    '''

    for read1 in fastq_iterator1:
        read2 = next(fastq_iterator2)
        pair_id = read1.identifier.split()[0]
        if not strict:
            while read2.identifier.split()[0] != pair_id:
                read2 = next(fastq_iterator2)
        if not read2.identifier.split()[0] == pair_id:
            raise ValueError("\nRead pairs do not match\n%s != %s" %
                             (pair_id, read2.identifier.split()[0]))
        yield (read1, read2)


def summary4fastp(fastp_json):
    """
    grep information from fastp output json.

    Usage: scSeq.py summary4fastp <fastp_json>

    """
    name = os.path.basename(fastp_json).split(".")[0]
    outfile = open(name + "_QCsummary.xls", "w")
    with open(fastp_json) as f:
        infoDic = json.load(f)
        summaryDic = infoDic["summary"]
        filterDic = infoDic["filtering_result"]
        before_filtering = summaryDic["before_filtering"]
        after_filtering = summaryDic["after_filtering"]

        raw_reads = before_filtering["total_reads"] if before_filtering[
                "total_reads"] else "NA"
        raw_bases = before_filtering["total_bases"] if before_filtering[
                "total_bases"] else "NA"
        rawQ30 = "%.2f" % (100 * float(before_filtering["q30_rate"]))
        rawGC = "%.2f" % (100 * float(before_filtering["gc_content"]))
        clean_reads = after_filtering["total_reads"] if after_filtering[
                "total_reads"] else "NA"
        clean_bases = after_filtering["total_bases"] if after_filtering[
                "total_bases"] else "NA"
        cleanQ30 = "%.2f" % (100 * float(after_filtering["q30_rate"]))
        cleanGC = "%.2f" % (100 * float(after_filtering["gc_content"]))
        passed_filter_reads = "%.2f" % ((
                100 * float(filterDic["passed_filter_reads"]) / float(
                raw_reads)))
        print("sample\traw_reads\traw_bases\traw_Q30\traw_GC\tclean_reads\t\
        clean_bases\tclean_Q30\tclean_GC\tpassed_filter_reads", file=outfile)
        print("%s\t%s\t%s\t%s%%\t%s%%\t%s\t%s\t%s%%\t%s%%\t%s%%" % (
            name, raw_reads, raw_bases, rawQ30, rawGC, clean_reads, clean_bases,
            cleanQ30, cleanGC, passed_filter_reads), file=outfile)
    outfile.close()
    return outfile


def extract_head_seq(prefix, seqnum, read1, read2=None):
    """
    Extract the pre-sequence from the fastq file.

    Usage: scSeq.py extract_head_seq <prefix> <pre_num> <R1_fastq.gz> <R2_fastq.gz=None>

    """
    outfile1 = gzip.open(prefix + "_head_R1.fastq.gz", "wt")
    with gzip.open(read1, "rb") as f1:
        read1_iterate = fastqIterate(f1)
        for read in read1_iterate:
            read = read.split("\n")
            outfile1.write("%s\n%s\n+\n%s\n" % (read[0], read[1][:int(seqnum)],
                                              read[3][:int(seqnum)]))
    outfile1.close()
    if read2:
        outfile2 = gzip.open(prefix + "_head_R2.fastq.gz", "wt")
        with gzip.open(read2, "rb") as f2:
            read2_iterate = fastqIterate(f2)
            for read in read2_iterate:
                read = read.split("\n")
                outfile2.write("%s\n%s\n+\n%s\n" % (read[0], read[1][:int(seqnum)],
                                                  read[3][:int(seqnum)]))

        outfile2.close()
        return outfile1, outfile2
    else:
        return outfile1


def get_single_cell_read_count(prefix, infile):
    """
    Calculate each single cell reads number.

    Usage: scSeq.py get_single_cell_read_count <prefix> <R1/2_fastq.gz>

    """
    scdic = {}
    outfile1 = open(prefix + "_cellReadsCounts.tsv", "w")
    outfile2 = open(prefix + "_cellReadsInfo.tsv", "w")
    outfile1.write("%s\t%s\n" % ("cell_ID", "reads_count"))
    outfile2.write("%s\t%s\t%s\n" % ("all_cells_count", "all_cells_reads",
                                    "average_reads_per_cell"))
    with gzip.open(infile, "rb") as f:
        read_iterate = fastqIterate(f)
        for read in read_iterate:
            barcode = read.split()[0].split("_")[1].strip()
            scdic[barcode] = scdic.get(barcode, 0) + 1
    for i in sorted(scdic.items(), key=lambda item: item[1], reverse=True):
        print("%s\t%s" % (i[0], i[1]), file=outfile1)
    outfile1.close()
    outfile2.write("%i\t%i\t%.2f\n" % (len(scdic.keys()), sum(scdic.values()),
                                       sum(scdic.values())/len(scdic.keys())))
    outfile2.close()
    return outfile1, outfile2


def map_tsv2_all(prefix, intable):
    """
    Map gene in counts.tsv to all hg.38 genome genes, include blank gene.

    Usage: scSeq.py map_tsv2_all <prefix> <intable.tsv.gz>

    """
    _tableDic = {}
    _allgeneList = "/home/longzhao/scSeq/src/geneList.uniq.sorted"
    if os.path.exists(_allgeneList):
        pass
    else:
        raise Exception("allGeneListFoundError: please check and run again!")
    outfile = gzip.open(prefix + "_all_counts.tsv.gz", "wt")
    with gzip.open(intable, "rt") as f:
        for line in f:
            l = line.strip().split("\t")
            _tableDic[l[0]] = line.strip()
    cellNum = len(_tableDic["gene"].strip().split("\t")) - 1
    blankRow = ("\t" + "0\t"*cellNum)[:-1]
    print(_tableDic["gene"], file=outfile)
    with open(_allgeneList, "r") as f:
        for gene in f:
            gene = gene.strip()
            if gene in _tableDic:
                print(_tableDic[gene], file=outfile)
            else:
                print(gene + blankRow, file=outfile)
    return outfile


def sequence_saturation(prefix, metaData, cellReadsCount):
    """
    Calculate sequence saturation from final high quality cell.

    Usage: scSeq.py sequence_saturation <prefix> <metaData.csv> <cellReadsCount.tsv>

    """
    outfile = open(prefix + "_saturation.tsv", "w")
    metaDic = {}
    nGene, nUMI = [], []
    def median(aList):
        aList = sorted(aList)
        Len = len(aList)
        if Len == 0:
            return 0
        elif Len % 2 == 1:
            index = int(Len / 2) + 1
            return aList[index]
        elif Len % 2 == 0:
            index = int(Len / 2)
            return (aList[index] + aList[index + 1]) / 2
    with open(metaData, "r") as f:
        metareader = csv.reader(f)
        for line in metareader:
            if metareader.line_num == 1:
                continue
            else:
                metaDic[line[0].strip().split("_")[1]] = 0
                nGene.append(int(line[1]))
                nUMI.append(int(line[2]))
    with open(cellReadsCount, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            if line[0].strip() in metaDic:
                metaDic[line[0].strip()] = int(line[1])
            else:
                pass
    outfile.write("%s\t%s\n" % ("Item", "Score"))
    outfile.write("%s\t%i\n" % ("final_number_of_cells", len(metaDic.values())))
    outfile.write("%s\t%.2f\n" % (
    "mean_read_per_cell", float(sum(metaDic.values()) / len(metaDic.keys()))))
    outfile.write("%s\t%.2f\n" % ("median_gene_per_cell", median(nGene)))
    outfile.write("%s\t%.2f\n" % ("median_umi_per_cell", median(nUMI)))
    outfile.write("%s\t%.2f%%\n" % (
    "sequencing_saturation", (1 - sum(nUMI) / sum(metaDic.values())) * 100))
    return outfile

def grep_STAR_map_info(prefix, starLog):
    """
    grep Map stat information from STAR.Log.final.out.

    Usage: scSeq.py grep_STAR_map_info <prefix> <STAR.Log.final.out>

    """
    infoList = [
        "Number of input reads",
        "Average input read length",
        "Uniquely mapped reads number",
        "Uniquely mapped reads %",
        "Average mapped length",
        "Number of reads mapped to too many loci",
        "% of reads mapped to too many loci"
    ]
    fileDic = {}
    with open(prefix + "_mapstat.tsv", "w") as f:
        with open(starLog, "r") as f1:
            for line in f1:
                line = line.strip().split("|")
                if len(line) == 2:
                    fileDic[line[0].strip()] = line[1].strip()
        print("Item\tScore", file=f)
        for _ in infoList:
            print("%s\t%s" % (_, fileDic[_]), file=f)
    return 

def grep_seurat_cell_filter(prefix, seuratstdout):
    """
    grep filter information from seurat QC step.

    Usage: scSeq.py grep_seurat_cell_filter <prefix> <seuratstdout.txt>

    """
    pattern = "^(\d+) genes across (\d+) samples.$"
    out = "Item\tall_gene_counts\tall_cell_counts\n"
    tmp = []
    with open(prefix + "_cellFiltered.tsv", "w") as f:
        with open(seuratstdout, "r") as f1:
            for line in f1:
                line = line.strip()
                match = re.compile(pattern).search(line)
                if match:
                    tmp.append([match.group(1), match.group(2)])
                else:
                    pass
        if len(tmp) >= 2:
            out += ("Before filter\t%s\t%s\nAfter filter\t%s\t%s" %
                    (tmp[0][0], tmp[0][1], tmp[1][0], tmp[1][1]))
        else:
            out += "Before filter\tNA\tNA\nAfter filter\tNA\tNA"
        print(out, file=f)
    return


def generate_json(prefix, templateJson):
    """
    generate result.json from work directory

    Usage: scSeq.py generate_json <prefix> <templateJson>

    """
    with open(templateJson, "r") as f:
        tmpjson = f.read()
    now = time.strftime("%Y-%m-%d", time.localtime())
    tmpjson = tmpjson.replace("time", now)
    tmpjson = tmpjson.replace("prefix", prefix)
    with open(prefix + "_result.json", "w") as f:
        f.write(tmpjson)
    return 














