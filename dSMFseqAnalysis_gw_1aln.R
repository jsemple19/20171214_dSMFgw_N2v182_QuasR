## mapping with QuasR
####### Align bisulfite sequencing of dSMF amplicons to C elegans genome
# date: 2017-06-28
# author: Jennnifer Semple
# Scripts adapted from Drosophila scripts kindly provided by Arnaud Krebs
#
# input required:
# fasta files from sequencing in rawData directory
# sampleList.txt with FileName1 FileName2 SampleName headings in rawData directory
#
# output files include:
# alignment .bam files (contain a unique automatically generated ref num) in aln/ directory
# QC plot: ./plots/QC_QualityTrimmed_"date".pdf in plots/ directory
# list of aligned files to recall the project in other scripts: QuasR_Aligned.txt in current dir


library(QuasR)
library("BSgenome.Celegans.UCSC.ce11")
library("BSgenome.Ecoli.NCBI.20080805")
# collect citations for packages used
packageBib<-toBibtex(c(citation("QuasR"),
                       citation("BSgenome.Celegans.UCSC.ce11")))

setwd("/Users/semple/Documents/MeisterLab/sequencingData/20171214_dSMFgw_N2v182")

source('./R/callAllCs.r') #to call all Cs
source('./R/useful_functionsV1.r') #load the ranges

# make clusters, needed to run in parallel
cluObj=makeCluster(3)

#setup directory structure this to desired location of your alignments
if (!dir.exists(paste0(path,"/aln"))){
  dir.create(paste0(path,"/aln"))
}
if (!dir.exists(paste0(path,"/tmp"))) {  #for trimmomatic quality trimmed reads
  dir.create(paste0(path,"/tmp"))
}
if (!dir.exists(paste0(path,"/cutadapt"))) {  #for trimmomatic quality trimmed reads
  dir.create(paste0(path,"/cutadapt"))
}
if (!dir.exists(paste0(path,"/rds"))) {
  dir.create(paste0(path,"/rds"))
}
if (!dir.exists(paste0(path,"/plots"))) {
  dir.create(paste0(path,"/plots"))
}
if (!dir.exists(paste0(path,"/bed"))) {
  dir.create(paste0(path,"/bed"))
}
path='.'
my.alignmentsDir=paste(path,'/aln/',sep='')
tmp=paste(path,'/tmp/',sep='')

#create auxiliary file
export(Ecoli,"/tmp/Ecoli.fasta")

AuxInput=as.data.frame(cbind(
  FileName=paste(tmp,"Ecoli.fasta",sep=''),
  AuxName="Ecoli"))

write.table(AuxInput,paste0(path,'QuasR_auxInput.txt'),quote=F,row.names=F,sep='\t')

# create sampleList file with ls *.fastq.gz > sampleList.txt  edit file in VI to add columns
#load the experiments
seqExp=read.table(paste0(path,'rawData/sampleList.txt'),sep='\t',header=T,stringsAsFactors=F)

#create the QuasR Aln file
samples=as.data.frame(cbind(FileName1=paste(path,"rawData/",seqExp$FileName1,sep=''),
                            FileName2=paste(path,"rawData/",seqExp$FileName2,sep=''),
                            SampleName=as.character(seqExp$SampleName)),stringsAsFactors=F)


###############################
# trim adaptors with cutadapt #
###############################
file.remove(paste0(path,"/cutadapt/cutadapt_log.txt"))
for(i in sl(samples[,1])){
  #i=1
  #spID=as.character(samples$SampleName[i])
  #clip the low quality bases #remove adapters
  system(paste(
    path, '/runCutadapt.sh ',
    samples$FileName1[i], ' ', samples$FileName2[i],
    sep='')
  )
}

#create the QuasR Aln table
samples=as.data.frame(cbind(FileName1=paste(path,"/cutadapt/",seqExp$FileName1,sep=''),
                            FileName2=paste(path,"/cutadapt/",seqExp$FileName2,sep=''),
                            SampleName=as.character(seqExp$SampleName)),stringsAsFactors=F)


###########################
#trim the low quality bases
###########################

for(i in sl(samples[,1])){
  #i=1
  spID=as.character(samples$SampleName[i])
  #clip the low quality bases #remove adapters
  system(paste(
    'java -jar $HOME/Trimmomatic-0.36/trimmomatic-0.36.jar PE ',
    samples$FileName1[i],' ', samples$FileName2[i], ' ',
    './tmp/',samples$SampleName[i],'_forward_paired.fq.gz ',
    './tmp/',samples$SampleName[i],'_forward_unpaired.fq.gz ',
    './tmp/',samples$SampleName[i],'_reverse_paired.fq.gz ',
    './tmp/',samples$SampleName[i],'_reverse_unpaired.fq.gz ',
    'ILLUMINACLIP:$HOME/Trimmomatic-0.36/adapters/TruSeq_2-3_PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36',
    sep='')
  )
}


AlnInput=as.data.frame(cbind(
  FileName1=paste(tmp,samples$SampleName,'_forward_paired.fq.gz',sep=''),
  FileName2=paste(tmp,samples$SampleName,'_reverse_paired.fq.gz',sep=''),
  SampleName=as.character(samples$SampleName)))

write.table(AlnInput,paste0(path,'QuasR_input.txt'),quote=F,row.names=F,sep='\t')



###########################
#Align the full length fragments
###########################
QuasRdef='-k 2 --best --strata'

NOMEproj<-qAlign(sampleFile="QuasR_input.txt",
                 genome="BSgenome.Celegans.UCSC.ce11",
                 auxiliaryFile="QuasR_auxInput.txt",
                 aligner="Rbowtie",
                 paired="fr",
                 bisulfite="dir", #for gw
                 projectName="dSMF_gw_N2vF2",
                 alignmentsDir=my.alignmentsDir,
                 clObj=cluObj ,
                 alignmentParameter=paste('-e 150 -X 600 ',QuasRdef,sep=''),
                 # bowtie_usage() for parameters:
                 #"  -k <int>           report up to <int> good alignments per read (default: 1)"
                 #"  --best             hits guaranteed best stratum; ties broken by quality"
                 #"  --strata           hits in sub-optimal strata aren't reported (requires --best)"
                 #"  -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)"
                 #"  -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)"
                 cacheDir = tmp)



#QC of the alignment
#todayDate<-format(Sys.time(), "%Y%m%d")
qQCReport(NOMEproj,paste0(path,'/plots/QC_QualityTrimmed.pdf'),clObj=cluObj)


NOMEproj
alnStats<-as.data.frame(alignmentStats(NOMEproj))
alnStats$perCentMapped<-round(100*alnStats$mapped/(alnStats$mapped+alnStats$unmapped),2)
alnStats

##### with cutadapt and then trimmomatic
# seqlength mapped unmapped perCentMapped
# N2_DE_gwV006:genome 100286401 172494   290962         37.22
# N2_DE_gwV007:genome 100286401 175360   292334         37.49
# F2_DE_gwV008:genome 100286401 109988   359530         23.43
# F2_DE_gwV009:genome 100286401 180370   287544         38.55
# N2_DE_gwV006:Ecoli   64754917   1046   289916          0.36
# N2_DE_gwV007:Ecoli   64754917    974   291360          0.33
# F2_DE_gwV008:Ecoli   64754917   9342   350188          2.60
# F2_DE_gwV009:Ecoli   64754917   2400   285144          0.83


#### with trimmomatic only

#undir
# seqlength mapped unmapped perCentMapped
# N2_DE_gwV006:genome 100286401 168284   109322         60.62
# N2_DE_gwV007:genome 100286401 171744   103214         62.46
# F2_DE_gwV008:genome 100286401 107482   183262         36.97
# F2_DE_gwV009:genome 100286401 177240   120514         59.53
# N2_DE_gwV006:Ecoli   64754917   1036   108286          0.95
# N2_DE_gwV007:Ecoli   64754917    964   102250          0.93
# F2_DE_gwV008:Ecoli   64754917   8980   174282          4.90
# F2_DE_gwV009:Ecoli   64754917   2416   118098          2.00

#dir - doesn't seem to make a difference to mapping...?!
# seqlength mapped unmapped perCentMapped
# N2_DE_gwV006:genome 100286401 168282   109324         60.62
# N2_DE_gwV007:genome 100286401 171740   103218         62.46
# F2_DE_gwV008:genome 100286401 107460   183284         36.96
# F2_DE_gwV009:genome 100286401 177250   120504         59.53
# N2_DE_gwV006:Ecoli   64754917    996   108328          0.91
# N2_DE_gwV007:Ecoli   64754917    922   102296          0.89
# F2_DE_gwV008:Ecoli   64754917   9094   174190          4.96
# F2_DE_gwV009:Ecoli   64754917   2386   118118          1.98

#with new adapter file
# seqlength mapped unmapped perCentMapped
# N2_DE_gwV006:genome 100286401 168292   106908         61.15
# N2_DE_gwV007:genome 100286401 171786   101788         62.79
# F2_DE_gwV008:genome 100286401 107510   179624         37.44
# F2_DE_gwV009:genome 100286401 177256   118928         59.85
# N2_DE_gwV006:Ecoli   64754917   1036   105872          0.97
# N2_DE_gwV007:Ecoli   64754917    998   100790          0.98
# F2_DE_gwV008:Ecoli   64754917   9080   170544          5.06
# F2_DE_gwV009:Ecoli   64754917   2386   116542          2.01


alignments1=as.data.frame(alignments(NOMEproj)$genome) #pulls out name of .bam files created

unlink(c('QuasR_Aligned.txt'))
write.table(alignments1,paste0(path,'QuasR_Aligned.txt'),quote=F,col.names=T,row.names=F,sep='\t',append=T)


