rm(list=ls())
setwd(ls)

pgfFile <- read.table(file="Masks/MoEx_newpgffile.pgf",header=TRUE,sep="\t",fill=TRUE)
## 4,779,547 probes

ilsHits <- read.table(file="SNPs_Probes/ILS_snpsInProbes.txt",sep="\t",header=FALSE)
ilsHits <- ilsHits[!grepl("random",ilsHits$V1),]
ilsHits <- ilsHits[!(ilsHits$V1 %in% c("chrM","chrUn")),]
#175,132 of 4,149,331 total probes with locations

issHits <- read.table(file="SNPs_Probes/ISS_snpsInProbes.txt",sep="\t",header=FALSE)
issHits <- issHits[!grepl("random",issHits$V1),]
issHits <- issHits[!(issHits$V1 %in% c("chrM","chrUn")),]
#152,905 of 4,171,558 total probes with locations

ProbesWithSNP <- rbind(ilsHits[ilsHits$V5>0,1:5],issHits[issHits$V5>0,1:5])
ProbesWithSNP <- ProbesWithSNP[!(duplicated(ProbesWithSNP[,1:5])),]
dim(ProbesWithSNP)


uniqueHits <- rbind(ilsHits[,1:5],issHits[,1:5])
uniqueHits <- uniqueHits[!(duplicated(uniqueHits[,1:5])),]

probeStatus <- merge(pgfFile,uniqueHits,by.x="probeID",by.y="V4",all.x=TRUE)
probeStatus$status <- rep("unique Hit",nrow(probeStatus))
probeStatus$status[is.na(probeStatus$V1)] <- "nonUnique"
table(probeStatus$status)
## 4,485,273 probes uniquely hit the genome
## 455,085 probes did NOT uniquely hit the genome


probeStatus$status[probeStatus$V5>0 & probeStatus$status=="unique Hit"] <- "with snp"
## 245,574 uniquelly mapping probes with known SNP

table(probeStatus$status)

#4,779,547 original probes

# 455,085 did not hit the genome or hit multiple places in the genome
# 245,574 probes contain a SNP between SHR and BN
#__________________________________
# 700,659 probes to be removed
 
to.be.removed = probeStatus[probeStatus$status!="unique Hit",]

#########################################
###  Remove Probes from the pgf File  ###
#########################################


removedProbes <- to.be.removed[to.be.removed$type=="main",]
# 656,536 of probes to be removed are "main" probes

##  removing bad probes
new.pgf <- pgfFile[!(pgfFile$probeID %in% removedProbes$probeID),]

##  identifying probe sets that have less than 3 probes
numPerSet <- table(new.pgf$probeSetID[new.pgf$type=="main"])
probeSets.to.remove <- rownames(numPerSet)[numPerSet<3]

##  removing all probes from bad probe sets
new.pgf <- new.pgf[!(new.pgf$probeSetID %in% probeSets.to.remove),]

length(unique(pgfFile$probeSetID)) - length(unique(new.pgf$probeSetID))
##  204,720 probe sets removed  ##
write.table(new.pgf,file="Masks/LXS.masked.pgf.simple.txt",sep="\t",row.names=FALSE,quote=FALSE)

#############################
###  Retained Probe Sets  ###
#############################
rm(list=ls())
new.pgf <- read.table(file="Masks/LXS.masked.pgf.simple.txt",sep="\t",header=TRUE)
retainProbeSet <- unique(new.pgf$probeSetID)

##  Probe Set Files

core.ps <- read.table(file="Source/MoEx-1_0-st-v1.r2.dt1.mm9.core.ps",sep="\t",comment.char="")
reduced.core.ps <- c(as.character(core.ps[1:12,]),as.character(core.ps[core.ps[,1] %in% retainProbeSet,]))

nrow(core.ps) - length(reduced.core.ps)
#31,707 core probe sets removed

write.table(reduced.core.ps,file="Masks/LXS.MoEx-1_0-st-v1.r2.dt1.mm10.core.MASKED.ps",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

extended.ps <- read.table(file="Source/MoEx-1_0-st-v1.r2.dt1.mm9.extended.ps",sep="\t",comment.char="")
reduced.extended.ps <- c(as.character(extended.ps[1:12,]),as.character(extended.ps[extended.ps[,1] %in% retainProbeSet,]))

nrow(extended.ps) - length(reduced.extended.ps)
#99,879 core/extended probe sets removed

write.table(reduced.extended.ps,file="Masks/LXS.MoEx-1_0-st-v1.r2.dt1.mm10.extended.MASKED.ps",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

full.ps <- read.table(file="Source/MoEx-1_0-st-v1.r2.dt1.mm9.full.ps",sep="\t",comment.char="")
reduced.full.ps <- c(as.character(full.ps[1:12,]),as.character(full.ps[full.ps[,1] %in% retainProbeSet,]))
nrow(full.ps) - length(reduced.full.ps)

#204,479 full probe sets removed

write.table(reduced.full.ps,file="Masks/LXS.MoEx-1_0-st-v1.r2.dt1.mm10.full.MASKED.ps",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

reduced.all.ps <- c(as.character(full.ps[1:12,]),as.character(retainProbeSet))
write.table(reduced.all.ps,file="Masks/LXS.MoEx-1_0-st-v1.r2.dt1.mm10.all.MASKED.ps",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

##  Transcript Clusters

core.mps <- read.table(file="Source/MoEx-1_0-st-v1.r2.dt1.mm9.core.mps",sep="\t",header=TRUE)

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
write.table(forPython,file="Masks/LXS.core.mps",quote=FALSE,sep="\t",na="0",row.names=FALSE,col.names=TRUE)

###  Extended Transcript Clusters  ###

extended.mps <- read.table(file="Source/MoEx-1_0-st-v1.r2.dt1.mm9.extended.mps",sep="\t",header=TRUE)

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
write.table(forPython,file="Masks/LXS.extended.mps",quote=FALSE,sep="\t",na="0",row.names=FALSE,col.names=TRUE)

###  Full Transcript Clusters  ###
rm(list=ls())
new.pgf <- read.table(file="Masks/LXS.masked.pgf.simple.txt",sep="\t",header=TRUE)
retainProbeSet <- unique(new.pgf$probeSetID)

full.mps <- read.table(file="Source/MoEx-1_0-st-v1.r2.dt1.mm9.full.mps",sep="\t",header=TRUE)

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
write.table(forPython,file="Masks/full.mps",quote=FALSE,sep="\t",na="0",row.names=FALSE,col.names=TRUE)

