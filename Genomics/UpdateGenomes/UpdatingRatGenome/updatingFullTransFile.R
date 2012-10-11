rm(list=ls())

setwd("/Users/laurasaba/Documents/Affymetrix/SNP Masks/Exon Array/Rat")


#############################
###  Retained Probe Sets  ###
#############################

new.pgf <- read.table(file="masked.pgf.simple.txt",sep="\t",header=TRUE)
retainProbeSet <- unique(new.pgf$probeSetID)

###  Full Transcript Clusters  ###

full.mps <- read.table(file="Source/RaEx-1_0-st-v1.r2.dt1.rn4.full.mps",sep="\t",header=TRUE)

###  Comprehensive Transcript Clusters  ###

comp.mps <- read.table(file="Source/RaEx-1_0-st-v1.r2.dt1.rn4.comprehensive.mps",sep="\t",header=TRUE)


#combine full and comprehensive

all.mps <- rbind(full.mps,comp.mps)
all.mps$numProbeSets <- ceiling(nchar(as.character(all.mps$probeset_list))/8)

## entries that do not have a transcript cluster ID  ##
noTranscript <- all.mps[is.na(all.mps[,"transcript_cluster_id"]),]
noTranscript2 <- noTranscript[noTranscript[,"probeset_id"] %in% retainProbeSet,]

## entries that have a transcript cluster ID  ##
all.mps <- all.mps[!is.na(all.mps[,"transcript_cluster_id"]),]

## remove duplicate transcript clusters

redo <- c()	
for(i in 1:nrow(all.mps)){
	redo <- rbind(redo,cbind(all.mps[i,"transcript_cluster_id"],unlist(strsplit(as.character(all.mps[i,"probeset_list"])," "))))
	}	
	
colnames(redo)<-c("transcript_cluster_id","probeset_id")

#get rid of duplicates
redo <- redo[!duplicated(redo),]

#remove masked probe sets
reduced <- redo[as.character(redo[,2]) %in% retainProbeSet,]

#calculate number of probes per probeset
numPerSet <- as.matrix(table(new.pgf$probeSetID),nc=1)

reduced <- merge(reduced,numPerSet,by.x="probeset_id",by.y="row.names")

#combine for output

colnames(reduced)[colnames(reduced)=="V1"] = "probe_count"
forPython <- rbind(as.matrix(reduced),as.matrix(noTranscript2[,c("probeset_id","transcript_cluster_id","probe_count")]))

forPython <- forPython[order(forPython[,"transcript_cluster_id"]),]
write.table(forPython,file="all.mps",quote=FALSE,sep="\t",na="0",row.names=FALSE,col.names=TRUE)
