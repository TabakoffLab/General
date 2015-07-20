rm(list = ls())

setwd('/Volumes/LauraS')
options(stringsAsFactors=FALSE)

##  Identifying probes that have more than one hit on the genome  ##
hits <- read.table(file="Affymetrix/SNP.Masks/Mouse 430v2/Data/Mouse430_NCBI37hits-multiple.tab",header=FALSE)
multiMatch <- hits[-grep("AFFX",hits[,1]),]  # removing control probes
multiMatch = cbind(multiMatch,probeID = unlist(lapply(strsplit(multiMatch,split="-",fixed=TRUE),function(a) paste(a[1:3],collapse="-"))))

##  Identifying probes that do NOT hit the genome  ##
noHits <- read.table(file="Affymetrix/SNP.Masks/Mouse 430v2/Data/Mouse430_NCBI37hits-none.tab",header=FALSE)
noMatch <- noHits[-grep("AFFX", noHits[,1]),]  # removing control probes
noMatch = cbind(noMatch,probeID = unlist(lapply(strsplit(noMatch,split="-",fixed=TRUE),function(a) paste(a[1:3],collapse="-"))))

##  Identifying probes with a known SNP  ##
SNPs <- read.table(file="Affymetrix/SNP.Masks/Mouse 430v2/Data/snpsInMouse430v2Probes.16strains.mm9.txt",sep="\t",header=FALSE)
SNPs$probeID = unlist(lapply(strsplit(SNPs$V4,split=":",fixed=TRUE),function(a) gsub(";","",paste(a[3:5],collapse="-"))))

##  Combining the three sets of probes to be deleted  ##
combined <- c(as.character(multiMatch[,"probeID"]),as.character(noMatch[,"probeID"]),as.character(SNPs$probeID))
combined <- combined[!duplicated(combined)]  # eliminating duplicates

length(combined)
#73,894 probes

##  Separate into three variables (probeID, X, Y)  ##
reduced <- data.frame(ProbeID=unlist(lapply(strsplit(combined,split="-",fixed=TRUE),function(a) a[1])),X=unlist(lapply(strsplit(combined,split="-",fixed=TRUE),function(a) a[2])),Y=unlist(lapply(strsplit(combined,split="-",fixed=TRUE),function(a) a[3])))

##  Identify Probe Sets With Less than 4 Probes Remaining  ##
all.probes <- read.table(file="Affymetrix/SNP.Masks/Mouse 430v2/Data/Mouse430_2_probe.tab",sep='\t',header=TRUE)
remaining.probes <- all.probes[!(paste(all.probes[,1],all.probes[,2],all.probes[,3],sep="-") %in% combined),]

remove.probeset <- table(remaining.probes[,1])[table(remaining.probes[,1])<4]
length(remove.probeset)
#1,446 probe sets removed
sum(table(remove.probeset)*as.numeric(names(table(remove.probeset))))
#2,896 additional probes removed

remove.probeset.xy <- merge(remove.probeset,all.probes,by.x=0,by.y="ProbeSetName")
reduced2 <- as.matrix(remove.probeset.xy[,c(1,3,4)])
reduced2[,2] <- gsub(" ","",reduced2[,2])
reduced2[,3] <- gsub(" ","",reduced2[,3])
colnames(reduced2) <- colnames(reduced)

##  Mismatch probes need to be deleted also  ##
pm.to.be.deleted <- rbind(reduced,reduced2)
pm.to.be.deleted <- pm.to.be.deleted[!duplicated(pm.to.be.deleted),]
#76,790 probes

mm.to.be.deleted <- cbind(pm.to.be.deleted[,1:2],as.numeric(pm.to.be.deleted[,3])+1)
colnames(mm.to.be.deleted) <- colnames(pm.to.be.deleted)
final.to.be.deleted <- rbind(pm.to.be.deleted,mm.to.be.deleted)

##  Percent of perfect match probes that were deleted  ##
(dim(final.to.be.deleted)[1]/2)/dim(all.probes)[1]
# 15.5%

##  Write to txt file  ##
write.table(final.to.be.deleted,file="Affymetrix/SNP.Masks/Mouse 430v2/Masks/mask.v37.16strains.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)



