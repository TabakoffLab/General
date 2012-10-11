rm(list=ls())

for(i in 1:9){
	results <- read.table(file=paste("~/BLAT/rat/ExonResults/Set",i,"/combined.psl",sep=""),sep="\t",skip=5)
	results<-results[!duplicated(results),]
	colnames(results)<-c("match","misMatch","repMatch","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd","tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts")
	results$qDiff <- results$qEnd - results$qStart
	perfectBLAT <- results[results$match==25 & results$misMatch==0 & results$qDiff==25,]
	perfectBLAT$probeID <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(perfectBLAT$tName),split=":",fixed=TRUE),function(a) a[3])),split=";",fixed=TRUE),function(a) a[1]))
	perfectBLAT <- perfectBLAT[!(perfectBLAT$probeID %in% as.character(unique(perfectBLAT$probeID[duplicated(perfectBLAT$probeID)]))),] 
	perfectBLAT <- perfectBLAT[,c("qName","qStart","qEnd","probeID")]

	if(i==1) write.table(perfectBLAT,file="~/BLAT/rat/ExonResults/perfectMatches.txt",sep="\t",row.names=FALSE,col.names=FALSE,append=FALSE,quote=FALSE)
	if(i!=1) write.table(perfectBLAT,file="~/BLAT/rat/ExonResults/perfectMatches.txt",sep="\t",row.names=FALSE,col.names=FALSE,append=TRUE,quote=FALSE)
}

rm(list=ls())
snps <- read.table(file="~/Documents/BxH_HxB Rats/Sequencing/SHR_SNPs_with_ssID.txt",sep="\t",header=TRUE)
snps$end <- snps$position+1
snps <- snps[,c("chr","position","end","ref_allele","shr_allele")]
snps$quality <- rep(999,nrow(snps))

## BN-Lx SNPs ##
BNLx <- read.table(file="~/Documents/Affymetrix/SNP Masks/Exon Array/Rat/Source/BNLX.snps.bed",sep="\t",header=FALSE)
BNLx <- BNLx[!grepl("random",BNLx$V1),]
BNLx <- BNLx[!(BNLx$V1 %in% c("chrM","chrUn")),]
dim(BNLx)
#351,543 putative SNPs

sum(BNLx$V6>150)
#91,804 high quality SNPs

SHRH <- read.table(file="~/Documents/Affymetrix/SNP Masks/Exon Array/Rat/Source/SHRH.snps.bed",sep="\t",header=FALSE)
SHRH <- SHRH[!grepl("random", SHRH$V1),]
SHRH <- SHRH[!(SHRH$V1 %in% c("chrM","chrUn")),]
#4,303,450 putative SNPs

sum(SHRH$V6>150)
#2,439,182 high quality SNPs

DNAseq <- rbind(BNLx[BNLx$V6>150,1:6],SHRH[SHRH$V6>150,1:6])
#2,530,986
DNAseq <- DNAseq[!(duplicated(DNAseq[,1:5])),]
#2,472,027


allSNPs <- merge(DNAseq,snps,by.x=paste("V",1:5,sep=""),by.y=c("chr","position","end","ref_allele","shr_allele"),all=TRUE)
allSNPs$start <- allSNPs$V2-1

options(scipen=10)

write.table(allSNPs[,c("V1","start","V2","V4","V5")],file="~/BLAT/rat/ExonResults/combinedSNPs.BED",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)