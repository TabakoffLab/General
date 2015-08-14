#Initiate Python
#python 

import string

###  open output file  ###
outputFile = "/Volumes/Data/mm10/array/MoEx/Masks/MoEx-1_0-st-v1.r2.dt1.mm10.core.LXS.MASKED.mps"
fout = open(outputFile, "wb")

####################################################
###  Get the header rows from original mps file  ###
####################################################

###  open original pgf file  ###
fname = "/Volumes/Data/mm10/array/MoEx/Source/MoEx-1_0-st-v1.r2.dt1.mm9.core.mps"

###  write header rows to new pgf file  ###
fp = open(fname)
header=0
i=1
while header==0:
	tmp = fp.next()
	header = (tmp.find("probeset_id") > -1)
	fout.write("%s" % (tmp))
	i=i+1


##################################################
###  Reformat to fit Affymetrix PGF standards  ###
##################################################

###  new masked pgf file in simple format  ###
fname = "/Volumes/Data/mm10/array/MoEx/Masks/forPython/LXS.core.mps"

###  count then number of records  ###
fp = open(fname)
num_lines = sum(1 for line in open(fname))

fp = open(fname)
colNames = fp.next()

space = " "
j=0
transcript_cluster_id="a"
probeset_list = "tmp"
probe_count = 0
while j < (num_lines-1):
	j=j+1
	read = fp.next()
	read = string.split(read)
	if read[1]=="0":
		fout.write("%s\t\t%s\t\t\r\n" % (read[0],read[0]))
	if read[1]==transcript_cluster_id: 
		probeset_list = space.join([probeset_list,read[0]])
		probe_count = probe_count + map(int,read[2])
	if read[1]!=transcript_cluster_id:
		if j!=1: fout.write("%s\t%s\t%s\t\t%s\r\n" % (transcript_cluster_id,transcript_cluster_id,probeset_list,sum(probe_count))) 
		transcript_cluster_id=read[1]
		probeset_list=read[0]
		probe_count = map(int,read[2])



fout.close()



