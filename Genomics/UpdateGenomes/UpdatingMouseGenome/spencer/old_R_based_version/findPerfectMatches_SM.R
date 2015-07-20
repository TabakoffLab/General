rm(list=ls())
	
	
	results <- read.table(file="/Users/smahaffey/Desktop/MoEx1_0st_probe.default.psl",sep="\t",skip=5)
	results<-results[!duplicated(results),]
	colnames(results)<-c("match","misMatch","repMatch","nCount","qNumInsert","qBaseInsert","tNumInsert","tBaseInsert","strand","qName","qSize","qStart","qEnd","tName","tSize","tStart","tEnd","blockCount","blockSizes","qStarts","tStarts")
	results$qDiff <- results$qEnd - results$qStart
	perfectBLAT <- results[results$match==25 & results$misMatch==0 & results$qDiff==25,]
	#Exon1.0st
	perfectBLAT$probeID <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(perfectBLAT$qName),split=":",fixed=TRUE),function(a) a[3])),split=";",fixed=TRUE),function(a) a[1]))
	perfectBLAT <- perfectBLAT[!(perfectBLAT$probeID %in% as.character(unique(perfectBLAT$probeID[duplicated(perfectBLAT$probeID)]))),] 
	perfectBLAT <- perfectBLAT[,c("tName","tStart","tEnd","probeID","misMatch","strand")]
	write.table(perfectBLAT,file="/Users/smahaffey/Desktop/MoEx1_0st_probe.default.perfectMatches.wStrand.txt",sep="\t",row.names=FALSE,col.names=FALSE,append=FALSE,quote=FALSE)