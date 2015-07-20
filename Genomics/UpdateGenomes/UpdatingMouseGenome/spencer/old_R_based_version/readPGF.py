#Initiate Python
#python 

import string

###  open output file  ###
outputFile = "/Volumes/Data/mm10/masks/MoEx_newpgffile.pgf"
fout = open(outputFile, "wb")
fout.write("probeSetID\ttype\tIDnum\tprobeID\ttypeProbe\tgc_count\tprobe_length\tinterrogation_position\tprobe_sequence\n")

###  open input file  ###
fname = "/Volumes/Data/mm10/masks/MoEx-1_0-st-v1.r2.pgf"

###  count then number of records  ###
fp = open(fname)
num_lines = sum(1 for line in open(fname))

fp = open(fname)
header=0
i=1
while header==0:
	tmp = fp.next()
	header = (tmp.find("header2") > -1)
	i=i+1


j=0
while j < (num_lines - i):
	j=j+1
	read = fp.next()
	read = string.split(read)
	if len(read)==2: probeSetID = read[0] 
	if len(read)==2: type = read[1]
	if len(read)==1: IDnum = read[0]
	if len(read)==6: probeID = read[0]
	if len(read)==6: typeProbe = read[1]
	if len(read)==6: gc_count = read[2]
	if len(read)==6: probe_length = read[3]
	if len(read)==6: interrogation_position = read[4]
	if len(read)==6: probe_sequence = read[5]
	if len(read)==6: fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (probeSetID,type,IDnum,probeID,typeProbe,gc_count,probe_length,interrogation_position,probe_sequence))

fout.close()

