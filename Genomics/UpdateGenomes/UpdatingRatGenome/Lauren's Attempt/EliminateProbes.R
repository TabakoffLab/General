rm(list=ls())
setwd("/home/kiemele/HXB/UpdateRatGenome/try20121012")

pgfFile <- read.table(file="newpgffile.pgf",header=TRUE,sep="\t",fill=TRUE)
## 4,104,557 probes

uniqueHits <- read.table(file="snpsInProbes.txt",sep="\t",header=FALSE)

uniqueHits <- uniqueHits[!grepl("random",uniqueHits$V1),]
uniqueHits <- uniqueHits[!(uniqueHits$V1 %in% c("chrM","chrUn")),]

probeStatus <- merge(pgfFile,uniqueHits,by.x="probeID",by.y="V4",all.x=TRUE)
probeStatus$status <- rep("unique Hit",nrow(probeStatus))
probeStatus$status[is.na(probeStatus$V1)] <- "nonUnique"
table(probeStatus$status)
## 3,664,621 probes uniquely hit the genome
## 439,936 probes did NOT uniquely hit the genome


probeStatus$status[probeStatus$V5>0 & probeStatus$status=="unique Hit"] <- "with snp"
## 108,563 uniquelly mapping probes with known SNP

table(probeStatus$status)

#4,104,557 original probes

# 439,936 did not hit the genome or hit multiple places in the genome
# 108,563 probes contain a SNP between SHR and BN
#__________________________________
# 548,499 probes to be removed
 
to.be.removed = probeStatus[probeStatus$status!="unique Hit",]

#########################################
###  Remove Probes from the pgf File  ###
#########################################


removedProbes <- to.be.removed[to.be.removed$type=="main",]
# 497,275 of probes to be removed are "main" probes

##  removing bad probes
new.pgf <- pgfFile[!(pgfFile$probeID %in% removedProbes$probeID),]

##  identifying probe sets that have less than 3 probes
numPerSet <- table(new.pgf$probeSetID[new.pgf$type=="main"])
probeSets.to.remove <- rownames(numPerSet)[numPerSet<3]

##  removing all probes from bad probe sets
new.pgf <- new.pgf[!(new.pgf$probeSetID %in% probeSets.to.remove),]

length(unique(pgfFile$probeSetID)) - length(unique(new.pgf$probeSetID))
##  174,294 probe sets removed  ##
write.table(new.pgf,file="masked.pgf.simple.txt",sep="\t",row.names=FALSE,quote=FALSE)

#############################
###  Retained Probe Sets  ###
#############################
rm(list=ls())
new.pgf <- read.table(file="masked.pgf.simple.txt",sep="\t",header=TRUE)
retainProbeSet <- unique(new.pgf$probeSetID)

##  Probe Set Files
# Using rn4 files because Affymetrix hasn't updated to rn5
core.ps <- read.table(file="RaEx-1_0-st-v1.r2.dt1.rn4.core.ps",sep="\t",comment.char="#", header=TRUE)
reduced.core.ps <- c(as.character(core.ps[1:12,]),as.character(core.ps[core.ps[,1] %in% retainProbeSet,]))

nrow(core.ps) - length(reduced.core.ps)
#14,007 core probe sets removed

write.table(reduced.core.ps,file="/home/kiemele/HXB/UpdateRatGenome/try20121012/results/RaEx-1_0-st-v1.r2.dt1.rn5.core.MASKED.ps",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

extended.ps <- read.table(file="RaEx-1_0-st-v1.r2.dt1.rn4.extended.ps",sep="\t",comment.char="#", header=TRUE)
reduced.extended.ps <- c(as.character(extended.ps[1:12,]),as.character(extended.ps[extended.ps[,1] %in% retainProbeSet,]))

nrow(extended.ps) - length(reduced.extended.ps)
#29,345 core/extended probe sets removed

write.table(reduced.extended.ps,file="/home/kiemele/HXB/UpdateRatGenome/try20121012/results/RaEx-1_0-st-v1.r2.dt1.rn4.extended.MASKED.ps",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

full.ps <- read.table(file="RaEx-1_0-st-v1.r2.dt1.rn4.full.ps",sep="\t",comment.char="#", header=TRUE)
reduced.full.ps <- c(as.character(full.ps[1:12,]),as.character(full.ps[full.ps[,1] %in% retainProbeSet,]))
nrow(full.ps) - length(reduced.full.ps)
#149,405 full probe sets removed

write.table(reduced.full.ps,file="/home/kiemele/HXB/UpdateRatGenome/try20121012/results/RaEx-1_0-st-v1.r2.dt1.rn4.full.MASKED.ps",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

reduced.all.ps <- c(as.character(full.ps[1:12,]),as.character(retainProbeSet))
write.table(reduced.all.ps,file="/home/kiemele/HXB/UpdateRatGenome/try20121012/results/RaEx-1_0-st-v1.r2.dt1.rn4.all.MASKED.ps",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

##  Transcript Clusters

core.mps <- read.table(file="RaEx-1_0-st-v1.r2.dt1.rn4.core.mps",sep="\t",header=TRUE)

## entries that do not have a transcript cluster ID  ##
noTranscript <- core.mps[is.na(core.mps[,"transcript_cluster_id"]),]
noTranscript2 <- noTranscript[noTranscript[,"probeset_id"] %in% retainProbeSet,]

## entries that have a transcript cluster ID  ##
core.mps <- core.mps[!is.na(core.mps[,"transcript_cluster_id"]),]
redo <- c()	
for(i in 1:nrow(core.mps)){
	redo <- rbind(redo,cbind(core.mps[i,"transcript_cluster_id"],unlist(strsplit(as.character(core.mps[i,"probeset_list"])," "))))
	}	
	
colnames(redo)<-c("transcript_cluster_id","probeset_id")
numPerSet <- as.matrix(table(new.pgf$probeSetID),nc=1)

reduced <- redo[as.character(redo[,2]) %in% retainProbeSet,]
reduced <- merge(reduced,numPerSet,by.x="probeset_id",by.y="row.names")

#combine for output

colnames(reduced)[colnames(reduced)=="V1"] = "probe_count"
forPython <- rbind(as.matrix(reduced),as.matrix(noTranscript2[,c("probeset_id","transcript_cluster_id","probe_count")]))

forPython <- forPython[order(forPython[,"transcript_cluster_id"]),]
write.table(forPython,file="/home/kiemele/HXB/UpdateRatGenome/try20121012/results/core.mps",quote=FALSE,sep="\t",na="0",row.names=FALSE,col.names=TRUE)

###  Extended Transcript Clusters  ###

extended.mps <- read.table(file="RaEx-1_0-st-v1.r2.dt1.rn4.extended.mps",sep="\t",header=TRUE)

## entries that do not have a transcript cluster ID  ##
noTranscript <- extended.mps[is.na(extended.mps[,"transcript_cluster_id"]),]
noTranscript2 <- noTranscript[noTranscript[,"probeset_id"] %in% retainProbeSet,]

## entries that have a transcript cluster ID  ##
extended.mps <- extended.mps[!is.na(extended.mps[,"transcript_cluster_id"]),]
redo <- c()	
for(i in 1:nrow(extended.mps)){
	redo <- rbind(redo,cbind(extended.mps[i,"transcript_cluster_id"],unlist(strsplit(as.character(extended.mps[i,"probeset_list"])," "))))
	}	
	
colnames(redo)<-c("transcript_cluster_id","probeset_id")
numPerSet <- as.matrix(table(new.pgf$probeSetID),nc=1)

reduced <- redo[as.character(redo[,2]) %in% retainProbeSet,]
reduced <- merge(reduced,numPerSet,by.x="probeset_id",by.y="row.names")

#combine for output

colnames(reduced)[colnames(reduced)=="V1"] = "probe_count"
forPython <- rbind(as.matrix(reduced),as.matrix(noTranscript2[,c("probeset_id","transcript_cluster_id","probe_count")]))

forPython <- forPython[order(forPython[,"transcript_cluster_id"]),]
write.table(forPython,file="/home/kiemele/HXB/UpdateRatGenome/try20121012/results/extended.mps",quote=FALSE,sep="\t",na="0",row.names=FALSE,col.names=TRUE)

###  Full Transcript Clusters  ###
rm(list=ls())
new.pgf <- read.table(file="masked.pgf.simple.txt",sep="\t",header=TRUE)
retainProbeSet <- unique(new.pgf$probeSetID)

full.mps <- read.table(file="RaEx-1_0-st-v1.r2.dt1.rn4.full.mps",sep="\t",header=TRUE)

## entries that do not have a transcript cluster ID  ##
noTranscript <- full.mps[is.na(full.mps[,"transcript_cluster_id"]),]
noTranscript2 <- noTranscript[noTranscript[,"probeset_id"] %in% retainProbeSet,]

## entries that have a transcript cluster ID  ##
full.mps <- full.mps[!is.na(full.mps[,"transcript_cluster_id"]),]
redo <- c()	
for(i in 1:nrow(full.mps)){
	redo <- rbind(redo,cbind(full.mps[i,"transcript_cluster_id"],unlist(strsplit(as.character(full.mps[i,"probeset_list"])," "))))
	}
	
colnames(redo)<-c("transcript_cluster_id","probeset_id")
numPerSet <- as.matrix(table(new.pgf$probeSetID),nc=1)

reduced <- redo[as.character(redo[,2]) %in% retainProbeSet,]
reduced <- merge(reduced,numPerSet,by.x="probeset_id",by.y="row.names")

#combine for output

colnames(reduced)[colnames(reduced)=="V1"] = "probe_count"
forPython <- rbind(as.matrix(reduced),as.matrix(noTranscript2[,c("probeset_id","transcript_cluster_id","probe_count")]))

forPython <- forPython[order(forPython[,"transcript_cluster_id"]),]
write.table(forPython,file="/home/kiemele/HXB/UpdateRatGenome/try20121012/results/full.mps",quote=FALSE,sep="\t",na="0",row.names=FALSE,col.names=TRUE)

