#Initiate Python
#python 

import string

###  open output file  ###
outputFile = "/Volumes/Data/mm10/array/MoEx/Masks/MoEx-1_0-st-v1.r2.LXS.MASKED.pgf"
fout = open(outputFile, "wb")

####################################################
###  Get the header rows from original pgf file  ###
####################################################

###  open original pgf file  ###
fname = "/Volumes/Data/mm10/array/MoEx/Source/MoEx-1_0-st-v1.r2.pgf"

###  write header rows to new pgf file  ###
fp = open(fname)
header=0
i=1
while header==0:
	tmp = fp.next()
	header = (tmp.find("header2") > -1)
	fout.write("%s" % (tmp))
	i=i+1


##################################################
###  Reformat to fit Affymetrix PGF standards  ###
##################################################

###  new masked pgf file in simple format  ###
fname = "/Volumes/Data/mm10/array/MoEx/Masks/forPython/LXS.masked.pgf.simple.txt"

###  count then number of records  ###
fp = open(fname)
num_lines = sum(1 for line in open(fname))

fp = open(fname)
colNames = fp.next()

j=0
probeSetID="a"
IDnum=-1
while j < (num_lines-1):
	j=j+1
	read = fp.next()
	read = string.split(read)
	if read[0]!=probeSetID: 
		fout.write("%s\t%s\r\n" % (read[0],read[1]))
		probeSetID=read[0]
	if read[0]==probeSetID: 
		if IDnum!=read[2]:
			fout.write("\t%s\r\n" % (read[2]))
			IDnum=read[2]
		fout.write("\t\t%s\t%s\t%s\t%s\t%s\t%s\r\n" % (read[3],read[4],read[5],read[6],read[7],read[8]))



fout.close()

