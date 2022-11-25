#! /usr/bin/env python

# Usage: renameGFF.py -i file.gff --contig-table name_table.tsv
# name_table.tsv is a 2-column tab-separated file with the first item in each line as the original contig name, second item is the new contig name
import sys

assemblyID = ""

for arg in sys.argv:
	if arg in ["-i","--in"]:
		gff = sys.argv[sys.argv.index(arg)+1]
	elif arg == "--contig-table":
		renameTable = open(sys.argv[sys.argv.index(arg)+1], "r")
	elif arg == "--assembly-id":
		assemblyID = sys.argv[sys.argv.index(arg)+1]

gffFilename = ".".join(gff.split(".")[0:-1])
outfile = open("%s_renamed.gff" % gtfFilename, "w")
if len(assemblyID) >= 1:
	outfile.write("#assembly ID = %s\n" % assemblyID)

renameDict = {}
for line in renameTable:
	oldContig = line.strip("\n").split("\t")[0]
	oldContig = oldContig.split(" ")[0]
	newContig = line.strip("\n").split("\t")[1]
	newContig = newContig.split(" ")[0]
	renameDict[oldContig] = newContig

with open(gff,"r") as ingff:
	for line in ingff:
		if line.startswith("#"):
			outfile.write(line)
		else:
			line = line.strip("\n")
			contig = line.split("\t")[0]
			restofline = "\t".join(line.split("\t")[1:])
			newContig = renameDict[contig]
			outfile.write("%s\t%s\n" %(newContig,restofline))
outfile.close()
renameTable.close()
