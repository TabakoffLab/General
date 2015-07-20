
## ILS SNPs ##
ILS <- read.table(file="~/mm10/array/430v2/SNPs_Probes/ILS_snpsInProbes.txt",sep="\t",header=FALSE)
ILS <- ILS[!grepl("random",ILS$V1),]
ILS <- ILS[!(ILS$V1 %in% c("chrM","chrUn")),]
dim(ILS)
sum(ILS$V5==0)
#406976 wo SNPs out of 417834


ISS <- read.table(file="~/mm10/array/430v2/SNPs_Probes/ISS_snpsInProbes.txt",sep="\t",header=FALSE)
ISS <- ISS[!grepl("random", ISS$V1),]
ISS <- ISS[!(ISS$V1 %in% c("chrM","chrUn")),]
dim(ISS)
sum(ISS$V5==0)
#408463 wo SNPs out of 417834

ProbesWithSNP <- rbind(ILS[ILS$V5>0,1:5],ISS[ISS$V5>0,1:5])
dim(ProbesWithSNP)
#20229
ProbesWithSNP <- ProbesWithSNP[!(duplicated(ProbesWithSNP[,1:5])),]
dim(ProbesWithSNP)
#15429

DNAseq <- rbind(ILS[ILS$V5==0,1:5],ISS[ISS$V5==0,1:5])
DNAseq <- DNAseq[!(duplicated(DNAseq[,1:5])),]
dim(DNAseq)
#412883

DNAseq <- DNAseq[!(DNAseq$V4 %in% ProbesWithSNP$V4),]
dim(DNAseq)
#402556

options(scipen=10)
write.table(DNAseq[,c("V1","V2","V3","V4","V5")],file="~/mm10/array/430v2/Aligned/LXS.nonMaskedLoc.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(ProbesWithSNP[,c("V1","V2","V3","V4","V5")],file="~/mm10/array/430v2/SNPs_Probes/LXS.wSNPs.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(ProbesWithSNP[,c("V4")],file="~/mm10/array/430v2/Mask/LXS.430v2.forMask.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)


multi <- read.table(file="~/mm10/Aligned/Mouse430_2.default.multipleMatches.txt",sep="\t")
multi<-multi[!duplicated(multi),1]
snp <- read.table(file="~/mm10/array/430v2/Masks/LXS.430v2.forMask.txt",sep="\t")
snp<-snp[!duplicated(snp),1]
masklist <- c(multi[,1],snp[,1])
masklist<-masklist[!duplicated(masklist)]
mask$pid1<-unlist(lapply(strsplit(as.character(masklist),split=":",fixed=TRUE),function(a) a[1]))
mask$pid2<-unlist(lapply(strsplit(as.character(masklist),split=":",fixed=TRUE),function(a) a[2]))
mask$pid3<-unlist(lapply(strsplit(as.character(masklist),split=":",fixed=TRUE),function(a) a[3]))