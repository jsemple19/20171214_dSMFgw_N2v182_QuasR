#!/usr/bin/make -f
## making small subset of data to test mappers locally



##############################
########### VARIABLES #########
###############################
runDir := /scratch/cluster/monthly/jsemple/20171214_dSMFgw_N2v182_QuasR
scriptDir := /home/jsemple/20171214_dSMFgw_N2v182_QuasR
dataDir := /home/jsemple/archive/20171214_dSMFgw_N2v182
sampleNames := N2_DE_gwV006 N2_DE_gwV007 F2_DE_gwV008 F2_DE_gwV009

# construct lists of sequence file namess
bname :=  $(addprefix 180126_SNK268_A_L001_JIB-, 1 2 3 4)
longbname := $(addsuffix _R1, $(bname)) $(addsuffix _R2, $(bname))


#list of the final output files
objects := $(addsuffix .fastq.gz, $(addprefix ${runDir}/rawData/, $(longbname))) \
	${runDir}/rawData/sampleList.txt ${runDir}/R/callAllCs.r ${runDir}/R/useful_functionsV1.r \
	${runDir}/../publicData/GenomeVer/Ecoli/Ecoli.fasta \
	${runDir}/../publicData/GenomeVer/WS250/c_elegans.PRJNA13758.WS250.genomic.fa

###############################
########### RULES  ############
###############################

all: $(objects)

#use cleanall when you want to force rerunning all the analysis
cleanall:
	rm -f $(objects)

.PHONY: all cleanall 


#######################################################
## get initial read stats                            ##
#######################################################

#run fastqc on downloaded sequences
${runDir}/rawData/%.fastq.gz: ${dataDir}/%.fastq.gz
	mkdir -p ${runDir}/rawData
	cp $^ $@

${runDir}/rawData/sampleList.txt: $(addsuffix .fastq.gz, $(addprefix ${runDir}/rawData/, $(longbname)))
	ls ${runDir}/rawData/*_R1.fastq.gz > fileList1.txt
	ls ${runDir}/rawData/*_R2.fastq.gz > fileList2.txt
	printf '%s\n'  ${sampleNames} > sampleNames.txt
	echo "FileName1	FileName2	SampleName" > ${runDir}/rawData/sampleList.txt
	paste fileList1.txt fileList2.txt sampleNames.txt >> ${runDir}/rawData/sampleList.txt 
	rm fileList?.txt sampleNames.txt 

${runDir}/R/%.r: ${scriptDir}/R/%.r
	mkdir -p ${runDir}/R
	cp $^ $@

${runDir}/../publicData/GenomeVer/Ecoli/Ecoli.fasta: ${dataDir}/../publicData/GenomeVer/Ecoli/Ecoli.fasta
	mkdir -p ${runDir}/../publicData/GenomeVer/Ecoli
	cp $^ $@

${runDir}/../publicData/GenomeVer/WS250/c_elegans.PRJNA13758.WS250.genomic.fa: ${dataDir}/../publicData/GenomeVer/WS250/c_elegans.PRJNA13758.WS250.genomic.fa
	mkdir -p ${runDir}/../publicData/GenomeVer/WS250
	cp $^ $@

