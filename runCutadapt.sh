#! /bin/bash
file1=$1
file2=$2

newFile1=`echo $file1 | sed 's/rawData/cutadapt/'`
newFile2=`echo $file2 | sed 's/rawData/cutadapt/'`

CUTADAPT=/Users/semple/anaconda/bin/cutadapt
$CUTADAPT  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $newFile1 -p $newFile2 $file1 $file2 >> ./cutadapt/cutadapt_log.txt


