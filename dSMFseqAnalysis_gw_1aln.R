## mapping with QuasR
####### Align bisulfite sequencing of dSMF amplicons to C elegans genome
# date: 2018-02-20
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

args = commandArgs(trailingOnly=TRUE)

setwd(args[1])


## get genome file locations
genomeEcoli<-paste0(args[1],"/../publicData/GenomeVer/Ecoli/Ecoli.fasta")
genomeCelegans<-paste0(args[1],"/../publicData/GenomeVer/WS250/c_elegans.PRJNA13758.WS250.genomic.fa")


source('./R/callAllCs.r') #to call all Cs
source('./R/useful_functionsV1.r') #load the ranges

# make clusters, needed to run in parallel
threadNum=4
cluObj=makeCluster(threadNum)

#setup directory structure this to desired location of your alignments
if (!dir.exists("./aln")){
  dir.create("./aln")
}
if (!dir.exists("./tmp")) {  #for trimmomatic quality trimmed reads
  dir.create("./tmp")
}
if (!dir.exists("./rds")) {
  dir.create("./rds")
}
if (!dir.exists("./plots")) {
  dir.create("./plots")
}
if (!dir.exists("./bed")) {
  dir.create("./bed")
}
path='./'
my.alignmentsDir=paste(path,'aln/',sep='')
tmp=paste(path,'tmp/',sep='')

#create auxiliary file

AuxInput=as.data.frame(cbind(
  FileName=genomeEcoli,
  AuxName="Ecoli"))

write.table(AuxInput,paste0(path,'QuasR_auxInput.txt'),quote=F,row.names=F,sep='\t')

# create sampleList file with ls *.fastq.gz > sampleList.txt  edit file in VI to add columns
#load the experiments
seqExp=read.table(paste0(path,'rawData/sampleList.txt'),sep='\t',header=T,stringsAsFactors=F)

#create the QuasR Aln file
samples=as.data.frame(cbind(FileName1=seqExp$FileName1,
                            FileName2=seqExp$FileName2,
                            SampleName=as.character(seqExp$SampleName)),stringsAsFactors=F)

###########################
#trim the low quality bases
###########################

for(i in sl(samples[,1])){
  #i=1
  spID=as.character(samples$SampleName[i])
  #clip the low quality bases #remove adapters
  system(paste(
    'trimmomatic PE -threads ', threadNum, ' ',
    samples$FileName1[i],' ', samples$FileName2[i], ' ',
    './tmp/',samples$SampleName[i],'_forward_paired.fq.gz ',
    './tmp/',samples$SampleName[i],'_forward_unpaired.fq.gz ',
    './tmp/',samples$SampleName[i],'_reverse_paired.fq.gz ',
    './tmp/',samples$SampleName[i],'_reverse_unpaired.fq.gz ',
    'ILLUMINACLIP:${HOME}/TruSeq_2-3_PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36',
    sep='')
  )
}


#    'java -jar $HOME/Trimmomatic-0.36/trimmomatic-0.36.jar PE ',

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
                 genome=genomeCelegans,
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


alignments1=as.data.frame(alignments(NOMEproj)$genome) #pulls out name of .bam files created

unlink(c('QuasR_Aligned.txt'))
write.table(alignments1,'QuasR_Aligned.txt',quote=F,col.names=T,row.names=F,sep='\t',append=T)


#QC of the alignment
#todayDate<-format(Sys.time(), "%Y%m%d")
qQCReport(NOMEproj,paste0('./plots/QC_QualityTrimmed.pdf'),clObj=cluObj)


NOMEproj
alnStats<-as.data.frame(alignmentStats(NOMEproj))
alnStats$perCentMapped<-round(100*alnStats$mapped/(alnStats$mapped+alnStats$unmapped),2)
alnStats


