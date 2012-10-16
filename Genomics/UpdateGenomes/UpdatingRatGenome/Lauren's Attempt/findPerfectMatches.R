#######################################
# Objective #1: 					  #
# Find the Perfect Matched Alignments # 
#######################################

rm(list=ls())
setwd("/home/kiemele/HXB/UpdateRatGenome")

results <- read.table(file="output4.try2.bed",sep="\t")
dim(results)
# 10,190,795 rows 9 columns 
 
results<-results[!duplicated(results),]
# 10190795x9


#colnames(results)<-c("match","misMatch","repMatch","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd","tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts")
colnames(results)<-c("match", "misMatch", "chr", "qStart", "qEnd", "ProbeID", "ProbeX", "ProbeY", "strand")
results$qDiff <- results$qEnd - results$qStart
perfectBLAT <- results[results$match==25 & results$misMatch==0 & results$qDiff==25,]
#7,456,508 perfect BLAT matches

perfectBLAT <- perfectBLAT[!(perfectBLAT$ProbeID %in% as.character(unique(perfectBLAT$ProbeID[duplicated(perfectBLAT$ProbeID)]))),]
#3,668,092 perfect BLAT

perfectBLAT <- perfectBLAT[,c("chr", "qStart","qEnd", "ProbeID")]
write.table(perfectBLAT,file="perfectMatches.txt",sep="\t",row.names=FALSE,col.names=FALSE,append=FALSE,quote=FALSE)

#####################################################################
# Object #2:														#
# Merge SNP Data From DNA Seq 										#
# Note: Can always use the STAR SNPs instead  						#
# Just not as accruate because strains slightly different than ours #
#####################################################################

##SHRH ##
# Snps
setwd("/home/kiemele/HXB/UpdateRatGenome/try20121012")
rm(list=ls())
snps <- read.table(file="/data/helicos/RatDNA/trimmedReads/SHRH.Filtered.Snps.bed",sep="\t",header=FALSE)
colnames(snps) = c('chr', 'position', 'end', 'name', 'alleles', 'score')
# 3,638,975 SHRH snps total 

alleles = unlist(strsplit(as.vector(snps[,"alleles"]), ":", fixed=TRUE))
snps$ref_allele = alleles[seq(1, length(alleles), 2)]
snps$shr_allele = alleles[seq(2, length(alleles), 2)]

sum(snps$score>150)
SHRH.snps = snps[snps$score>150, ]
#2,232,156 high quality SHRH SNPs

SHRH.snps <- SHRH.snps[,c("chr","position","end","ref_allele","shr_allele")]

#Indels
indels <- read.table(file="/data/helicos/RatDNA/trimmedReads/SHRH.Filtered.Indels.bed",sep="\t",header=FALSE)
colnames(indels) = c('chr', 'position', 'end', 'name', 'insertion', 'score')
#1,154,000 SHRH indels 

#Asking Steve about locations, he reccomends the following:
indels$position = indels$position +1 
indels$end = indels$end +1 

insertion =  unlist(strsplit(as.vector(indels[,"insertion"]), ":", fixed=TRUE))
indels$ref_insertion = insertion[seq(1, length(insertion), 2)]
indels$shr_insertion = insertion[seq(2, length(insertion), 2)]

sum(indels$score>150)
SHRH.indels = indels[indels$score>150, ]
#490,098 high quality SHRH indels


##BNLx ##
# Snps
snps <- read.table(file="/data/helicos/RatDNA/trimmedReads/BNLX.Filtered.Snps.bed",sep="\t",header=FALSE)
colnames(snps) = c('chr', 'position', 'end', 'name', 'alleles', 'score')
# 112,159 BNLx snps total 

alleles = unlist(strsplit(as.vector(snps[,"alleles"]), ":", fixed=TRUE))
snps$ref_allele = alleles[seq(1, length(alleles), 2)]
snps$shr_allele = alleles[seq(2, length(alleles), 2)]

sum(snps$score>150)
BNLx.snps = snps[snps$score>150, ]
#20,371 high quality BNLx SNPs

BNLx.snps <- BNLx.snps[,c("chr","position","end","ref_allele","shr_allele")]

#Indels
indels <- read.table(file="/data/helicos/RatDNA/trimmedReads/BNLX.Filtered.Indels.bed",sep="\t",header=FALSE)
colnames(indels) = c('chr', 'position', 'end', 'name', 'insertion', 'score')
#131,420 BNLx indels 

#Asking Steve about locations, he reccomends the following:
indels$position = indels$position +1 
indels$end = indels$end +1 

insertion =  unlist(strsplit(as.vector(indels[,"insertion"]), ":", fixed=TRUE))
indels$ref_insertion = insertion[seq(1, length(insertion), 2)]
indels$shr_insertion = insertion[seq(2, length(insertion), 2)]

sum(indels$score>150)
BNLx.indels = indels[indels$score>150, ]
#35,329 high quality BNLx indels

###########################################
# Combine the SNPs and INDELs into 1 file #
###########################################

Indels.total = rbind(SHRH.indels, BNLx.indels)
snps.total = rbind(SHRH.snps, BNLx.snps)

#Get common colums
Indels.total = Indels.total[,c('chr', 'position', 'end', 'ref_insertion', 'shr_insertion')]
colnames(Indels.total) = c('chr', 'position', 'end', 'ref_allele', 'shr_allele')

#Combine 
all = rbind(snps.total, Indels.total)
options(scipen=10)

write.table(all, file="combinedSNPs.bed", sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)





############LS Original Code 
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


