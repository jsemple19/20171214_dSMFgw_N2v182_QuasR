#!/usr/bin/make -f
## making small subset of data to test mappers locally



###############################
########### VARIABLES #########
###############################

#genomefile:= ${HOME}/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa
#trimmomaticDIR := ${HOME}/Trimmomatic-0.36
#methIndGenomeFiles := $(addsuffix ${genomefile}.bwameth.ct2, .sa .amb .ann .pac .bwt)
bname :=  $(addprefix 180126_SNK268_A_L001_JIB-, 1 2 3 4)
longbname := $(addsuffix _R1, $(bname)) $(addsuffix _R2, $(bname))


#list of the final output files
objects := 	$(addsuffix .fastq.gz, $(addprefix rawData/, $(longbname)))


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
rawData/%.fastq.gz: fastq/%.fastq.gz
	mkdir -p rawData
	gunzip -c $^ | head -n 1000000 | gzip > $@ 
