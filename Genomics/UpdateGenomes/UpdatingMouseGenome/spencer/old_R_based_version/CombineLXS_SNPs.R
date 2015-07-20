
## ILS SNPs ##
ILS <- read.table(file="~/mm10/array/MoEx/SNPs_Probes/ILS_snpsInProbes.txt",sep="\t",header=FALSE)
ILS <- ILS[!grepl("random",ILS$V1),]
ILS <- ILS[!(ILS$V1 %in% c("chrM","chrUn")),]
dim(ILS)
sum(ILS$V5==0)
#4149331 wo SNPs out of 4324463


ISS <- read.table(file="~/mm10/array/MoEx/SNPs_Probes/ISS_snpsInProbes.txt",sep="\t",header=FALSE)
ISS <- ISS[!grepl("random", ISS$V1),]
ISS <- ISS[!(ISS$V1 %in% c("chrM","chrUn")),]
dim(ISS)
sum(ISS$V5==0)
#4171558 wo SNPs out of 4324463

ProbesWithSNP <- rbind(ILS[ILS$V5>0,1:5],ISS[ISS$V5>0,1:5])
dim(ProbesWithSNP)
#328037
ProbesWithSNP <- ProbesWithSNP[!(duplicated(ProbesWithSNP[,1:5])),]
dim(ProbesWithSNP)
#245547

DNAseq <- rbind(ILS[ILS$V5==0,1:5],ISS[ISS$V5==0,1:5])
DNAseq <- DNAseq[!(duplicated(DNAseq[,1:5])),]
dim(DNAseq)
#4239727
DNAseq <- DNAseq[!(DNAseq$V4 %in% ProbesWithSNP$V4),]
dim(DNAseq)
#4081162

options(scipen=10)
write.table(DNAseq[,c("V1","V2","V3","V4","V5")],file="~/mm10/array/MoEx/SNPs_Probes/LXS_combinedSNPs.BED",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)


