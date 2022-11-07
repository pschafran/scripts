#! /usr/bin/env python

# Usage: renameGTF.py -i file.gtf --contig-table name_table.tsv --gene-prefix Genus_species --assembly-id 6cde96438c2e713efa5c285e4fe3a62d
# name_table.tsv is a 2-column tab-separated file with the first item in each line as the original contig name, second item is the new contig name
import sys

genePrefix = ""
assemblyID = ""

for arg in sys.argv:
	if arg in ["-i","--in"]:
		ingtf = open(sys.argv[sys.argv.index(arg)+1], "r")
		gtf = sys.argv[sys.argv.index(arg)+1]
	elif arg == "--contig-table":
		renameTable = open(sys.argv[sys.argv.index(arg)+1], "r")
	elif arg == "--gene-prefix":
		genePrefix = sys.argv[sys.argv.index(arg)+1]
	elif arg == "--assembly-id":
		assemblyID = sys.argv[sys.argv.index(arg)+1]

#renameTable = open(sys.argv[1], "r")
#ingtf = open(sys.argv[0], "r")
gtfFilename = ".".join(gtf.split("."))[0:-1]
outfile = open("%s_renamed.gtf" % gtfFilename, "w")
if len(assemblyID) >= 1:
	outfile.write("#assembly ID = %s\n" % assemblyID)

renameDict = {}
for line in renameTable:
	oldContig = line.strip("\n").split("\t")[0]
	newContig = line.strip("\n").split("\t")[1]
	renameDict[oldContig] = newContig

for line in ingtf:
	line = line.strip("\n")
	if not line.startswith("#"):
		contig = line.split("\t")[0]
		source = line.split("\t")[1]
		feature = line.split("\t")[2]
		startPos = line.split("\t")[3]
		endPos = line.split("\t")[4]
		score = line.split("\t")[5]
		strand = line.split("\t")[6]
		frame = line.split("\t")[7]
		description = line.strip(";").split("\t")[8]
		newContig = renameDict[contig]
		if feature != "gene":
			splitDescription = description.split("; ")
			for item in splitDescription:
				splitSplitDescription = splitDescription.split(" ")
				if splitSplitDescription[0] == "transcript_id":
					transcript_id = splitSplitDescription[1].strip('''"''')
				elif splitSplitDescription[1] == "gene_id":
					gene_id = splitSplitDescription[1].strip('''"''')
			if len(genePrefix) >= 1:
				newGene = "genePrefix_%s" % gene_id
				newTranscript = "genePrefix_%s" % transcript_id
			else:
				newGene = gene_id
				newTranscript = transcript_id
			outfile.write('''%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\ttranscript_id "%s"; gene_id "%s"\n''' %(newContig, source, feature, startPos, endPos, score, strand, frame, newTranscript, newGene ))
		else:
			gene_id = description[0]
			if len(genePrefix) >= 1:
				newGene = "genePrefix_%s" % gene_id
				newTranscript = "genePrefix_%s" % transcript_id
			else:
				newGene = gene_id
				newTranscript = transcript_id
			outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(newContig, source, feature, startPos, endPos, score, strand, frame, newGene ))

		#transcript = transcriptSplit2.strip("\"")
		#gtfDict[transcript] = {"contig":contig, "startPos":startPos, "endPos":endPos, "direction":direction}

	else:
		outfile.write("%s\n" % line)
