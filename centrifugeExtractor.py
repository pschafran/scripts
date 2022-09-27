#! /home/ps997/miniconda3/bin/python

# Usage: centrifugeExtractor.py centrifuge.results taxIDsToReject.txt reads.fastq

import sys
from Bio import SeqIO

centResultFile = sys.argv[1]
taxIDsFile = sys.argv[2]
readsFile = sys.argv[3]

taxIDsList = []
keepReadsList = []
rejectReadsList = []

openTaxIDsFile = open(taxIDsFile, "r")
for line in openTaxIDsFile:
	taxIDsList.append(line.strip("\n"))
openTaxIDsFile.close()


openCentResultFile = open(centResultFile , "r")
for line in openCentResultFile:
	splitline = line.strip("\n").split("\t")
	if splitline[2] in taxIDsList:
		rejectReadsList.append(splitline[0])
	elif splitline[2] not in taxIDsList:
		keepReadsList.append(splitline[0])
openCentResultFile.close()

outputKeepFile = "%s.keep.fastq" % readsFile
outputRejectFile = "%s.reject.fastq" % readsFile

input_seq_iterator = SeqIO.parse(readsFile, "fastq")
keep_seq_iterator = (record for record in input_seq_iterator if record.id in keepReadsList)
SeqIO.write(keep_seq_iterator, outputKeepFile, "fastq")

input_seq_iterator = SeqIO.parse(readsFile, "fastq")
reject_seq_iterator = (record for record in input_seq_iterator if record.id in rejectReadsList)
SeqIO.write(reject_seq_iterator, outputRejectFile, "fastq")
