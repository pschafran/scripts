#!/usr/bin/python
# 2020 May 13 (C) Peter W Schafran ps997@cornell.edu
#
# Extracts a range of bases out of a FASTA file. If output sequence name not provided in 4th column of input file, output sequence name is 'contigX_positionX-X'
# Requires a tab-delimited file with columns: contig name \t start position \t end position \t (optional) output sequence name
# Usage: trimFastaBases.py sequence.fasta targetFile.txt

import sys
from Bio.Seq import Seq

openSeqFile = open(sys.argv[1], "r")
openOutFile = open("%s_extractedRegions.fasta" % sys.argv[2], "w")
openTargetFile = open(sys.argv[2], "r")

# check targetFile
for line in openTargetFile:
	if len(line.strip("\n").split("\t")) < 3:
		print("ERROR: missing elements in line'n")
		print(line)
		exit(1)
	else:
		if int(line.strip("\n").split("\t")[1]) > int(line.strip("\n").split("\t")[2]):
			pass
			#print("WARNING: Start position comes after end position. Output will be in forward direction\n")
			#print(line)
			#exit(1)
		elif int(line.strip("\n").split("\t")[1]) <= 0 or int(line.strip("\n").split("\t")[2]) <= 0:
			print("ERROR: Start and end position can't be negative or zero\n")
			print(line)
			exit(1)

seqDict = {}
seqList = []
for seqLine in openSeqFile:
	if seqLine.startswith(">"):
		try:
			if len(seqList) > 1:
				seq = "".join(seqList)
			else:
				seq = seqList[0]
			seqDict[seqName] = seq
		except:
			seqList = []
		seqName = seqLine.strip(">\n").split()[0]
		seqList = []
	else:
		seqList.append(seqLine.strip("\n"))
if len(seqList) > 1:
	seq = "".join(seqList)
else:
	seq = seqList[0]
#seqDict[seqName] = seq
seqDict[seqName] = list(seq)
#for i in seqDict.keys():
#	print("%s %s" %(i, len(seqDict[i])))
openTargetFile.seek(0)
for targetLine in openTargetFile:
	reverse = False
	targetContig = targetLine.strip(">\n").split("\t")[0]
	startPos = int(targetLine.strip(">\n").split("\t")[1])
	endPos = int(targetLine.strip(">\n").split("\t")[2])
	if endPos < startPos:
		tmpPos = endPos
		endPos = startPos
		startPos = tmpPos
		reverse = True
	# check for bad positions
	if endPos > len(seqDict[targetContig]):
		print("ERROR: End position is beyond end of sequence!")
		print(targetLine)
		print(len(seqDict[targetContig]))
		exit(1)
	if startPos == 1 and endPos == len(seqDict[targetContig]):
		print("ERROR: Whole sequence selected to be trimmed!")
		print(targetLine)
		exit(1)
	try:
		openOutFile.write(">%s\n" % targetLine.strip("\n").split("\t")[3])
	except:
		openOutFile.write(">%s_pos%s-%s\n" %(targetContig, startPos, endPos))
	baseCount = 0
	writeCount = 0
	seqSubset = seqDict[targetContig][startPos-1:endPos-1]
	if reverse == True:
		dna = Seq("".join(seqSubset))
		seqSubseq = str(dna.reverse_complement())
	openOutFile.write(seqSubset)
	#for base in seqDict[targetContig]:
	#	baseCount += 1
	#	if baseCount >= startPos and baseCount <= endPos:
	#		openOutFile.write(base)
	#		writeCount += 1
	#		if writeCount % 80 == 0:
	#			openOutFile.write("\n")
	openOutFile.write("\n")

openSeqFile.close()
openOutFile.close()
openTargetFile.close()
