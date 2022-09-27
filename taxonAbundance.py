#! /usr/bin/env python

# script.py assembly_info_abundance.txt Species1.fasta Species2.fasta Species3.fasta ...

import sys

abundanceFile = sys.argv[1]
taxonFiles = sys.argv[2:]

dict = {}
with open(abundanceFile,"r") as infile:
	for line in infile:
		if line.startswith("#"):
			pass
		else:
			splitline = line.strip("\n").split("\t")
			dict[splitline[0]] = {"length" : splitline[1]}
			dict[splitline[0]].update({"coverage" : splitline[2]})
			dict[splitline[0]].update({"abundance" : splitline[3]})

with open("abundance_by_taxon.txt", "w") as outfile:
	for file in taxonFiles:
		with open(file, "r") as infile:
			taxon = file.split("_")[-1].split(".")[0]
			abundance = 0
			for line in infile:
				if line.startswith(">"):
					seqname = line.strip(">\n").split(" ")[0]
					abundance = int(abundance) + int(dict[seqname]["abundance"])
			outfile.write("%s\t%s\n" %(taxon, abundance))

