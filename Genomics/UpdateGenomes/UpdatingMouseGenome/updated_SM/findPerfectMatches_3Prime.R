	rm(list=ls())
	results <- read.table(file="~/mm10/Mouse430_2.default.psl",sep="\t",skip=5)
	results<-results[!duplicated(results),]
	colnames(results)<-c("match","misMatch","repMatch","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd","tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts")
	results$qDiff <- results$qEnd - results$qStart
	#results 927376
	perfectBLAT <- results[results$match==25 & results$misMatch==0 & results$qDiff==25,]
	#perfectMatch 685425
	perfectBLAT$probeID <-unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(perfectBLAT$qName),split=":",fixed=TRUE),function(a) paste(a[3],a[4],a[5],sep=":"))),split=";",fixed=TRUE),function(a) a[1]))
	multiBLAT <- perfectBLAT[(perfectBLAT$probeID %in% as.character(unique(perfectBLAT$probeID[duplicated(perfectBLAT$probeID)]))),]
	dim(multiBLAT)
	perfectBLAT <- perfectBLAT[!(perfectBLAT$probeID %in% as.character(unique(perfectBLAT$probeID[duplicated(perfectBLAT$probeID)]))),]
	dim(perfectBLAT) 
	#unique alignments 417838
	perfectBLAT <- perfectBLAT[,c("tName","tStart","tEnd","probeID")]
	write.table(perfectBLAT,file="~/mm10/array/430v2/Aligned/Mouse430_2.default.perfectMatches.txt",sep="\t",row.names=FALSE,col.names=FALSE,append=FALSE,quote=FALSE)
	multiBLAT <- multiBLAT[,c("probeID")]
	write.table(multiBLAT,file="~/mm10/array/430v2/Aligned/Mouse430_2.default.multipleMatches.txt",sep="\t",row.names=FALSE,col.names=FALSE,append=FALSE,quote=FALSE)