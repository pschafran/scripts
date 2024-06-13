#! /usr/bin/env python
#
# 2023 August 22. Peter Schafran ps997@cornell.edu

# Usage: convertNexusToFasta.py file.nex

from Bio import SeqIO
import sys

filename = sys.argv[1]
try:
	if (filename.split(".")[-1]) != "nexus" and (filename.split(".")[-1]) != "nex":
		print("ERROR: Expected .nexus or .nex file extension")
		exit(1)
	else:
		fileprefix = ".".join(filename.split(".")[0:-1])
		SeqIO.convert(filename,"nexus","%s.fasta" %(fileprefix),"fasta")
except:
	print("ERROR: Something went wrong...this shouldn't happen")
	exit(1)

