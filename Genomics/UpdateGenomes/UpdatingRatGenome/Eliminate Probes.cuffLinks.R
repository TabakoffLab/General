setwd("/Users/laurasaba/Documents/Affymetrix/SNP Masks/Exon Array/Rat/")

#############################
###  Retained Probe Sets  ###
#############################

new.pgf <- read.table(file="Results/masked.pgf.simple.txt",sep="\t",header=TRUE)
retainProbeSet <- unique(new.pgf$probeSetID)
numPerSet <- as.matrix(table(new.pgf$probeSetID),nc=1)

##  Transcript Clusters

## entries that do not have a transcript cluster ID  ##
core.mps <- read.table(file="Source/RaEx-1_0-st-v1.r2.dt1.rn4.core.mps",sep="\t",header=TRUE)
noTranscript <- core.mps[is.na(core.mps[,"transcript_cluster_id"]),]
noTranscript2 <- noTranscript[noTranscript[,"probeset_id"] %in% retainProbeSet,]

## cuffLink transcript clusters  ##
cuffLinks <- read.table(file="Source/cuffLinksTransClust.uniqueToTranscript.txt",sep="\t",header=TRUE)
cuffLinks$transcript_cluster_id <- 1000000 + 10*cuffLinks$geneCuffID + cuffLinks$isoformCuffID
redo <- cuffLinks[,c("transcript_cluster_id","probeset_id")]

reduced <- redo[as.character(redo[,2]) %in% retainProbeSet,]
reduced <- merge(reduced,numPerSet,by.x="probeset_id",by.y="row.names")

#combine for output

colnames(reduced)[colnames(reduced)=="V1"] = "probe_count"
forPython <- rbind(as.matrix(reduced),as.matrix(noTranscript2[,c("probeset_id","transcript_cluster_id","probe_count")]))

forPython <- forPython[order(forPython[,"transcript_cluster_id"]),]
write.table(forPython,file="Results/cuffLinksUnique.mps",quote=FALSE,sep="\t",na="0",row.names=FALSE,col.names=TRUE)

## cuffLink transcript clusters - All OS ##
cuffLinks <- read.table(file="Source/cuffLinksTransClust.allPossible.txt",sep="\t",header=TRUE)
cuffLinks$transcript_cluster_id <- 1000000 + 10*cuffLinks$geneCuffID + cuffLinks$isoformCuffID
redo <- cuffLinks[,c("transcript_cluster_id","probeset_id")]

reduced <- redo[as.character(redo[,2]) %in% retainProbeSet,]
reduced <- merge(reduced,numPerSet,by.x="probeset_id",by.y="row.names")

#combine for output

colnames(reduced)[colnames(reduced)=="V1"] = "probe_count"
forPython <- rbind(as.matrix(reduced),as.matrix(noTranscript2[,c("probeset_id","transcript_cluster_id","probe_count")]))

forPython <- forPython[order(forPython[,"transcript_cluster_id"]),]
write.table(forPython,file="Results/cuffLinksAllPS.mps",quote=FALSE,sep="\t",na="0",row.names=FALSE,col.names=TRUE)

