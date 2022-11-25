#! /usr/bin/env python

# Usage: renameGTF.py -i file.gtf --contig-table name_table.tsv --interval 10 --digits 5 --gene-prefix Genus_species --assembly-id 6cde96438c2e713efa5c285e4fe3a62d
# name_table.tsv is a 2-column tab-separated file with the first item in each line as the original contig name, second item is the new contig name
# Interval is the gap in numbering between adjacent genes e.g. g00010, g00020, g00030 for interval = 10
import sys

genePrefix = ""
assemblyID = ""
interval = 1
digits = 1

for arg in sys.argv:
	if arg in ["-i","--in"]:
		gtf = sys.argv[sys.argv.index(arg)+1]
	elif arg == "--contig-table":
		renameTable = open(sys.argv[sys.argv.index(arg)+1], "r")
	elif arg == "--gene-prefix":
		genePrefix = sys.argv[sys.argv.index(arg)+1]
	elif arg == "--assembly-id":
		assemblyID = sys.argv[sys.argv.index(arg)+1]
	elif arg =="--interval":
		interval = int(sys.argv[sys.argv.index(arg)+1])
	elif arg =="--digits":
		digits = int(sys.argv[sys.argv.index(arg)+1])

#renameTable = open(sys.argv[1], "r")
#ingtf = open(sys.argv[0], "r")
gtfFilename = ".".join(gtf.split(".")[0:-1])
outfile = open("%s_renamed.gtf" % gtfFilename, "w")
if len(assemblyID) >= 1:
	outfile.write("# assembly ID = %s\n" % assemblyID)

renameDict = {}
contigOrderList = []
for line in renameTable:
	oldContig = line.strip("\n").split("\t")[0]
	oldContig = oldContig.split(" ")[0]
	newContig = line.strip("\n").split("\t")[1]
	newContig = newContig.split(" ")[0]
	renameDict[oldContig] = newContig
	contigOrderList.append(oldContig)

gtfDict = {}
#startPosList = []
#sortedGeneList = []
with open(gtf, "r") as ingtf:
	for line in ingtf:
		line = line.strip("\n")
		if not line.startswith("#"):
			contig = line.split("\t")[0]
			if contig not in gtfDict.keys():
				gtfDict.update({contig : {"startPosDict" : {}}})
				gtfDict[contig].update({"startPosList" : []})
				gtfDict[contig].update({"sortedGeneList" : []})
			#source = line.split("\t")[1]
			feature = line.split("\t")[2]
			startPos = int(line.split("\t")[3])
			#endPos = line.split("\t")[4]
			#score = line.split("\t")[5]
			#strand = line.split("\t")[6]
			#frame = line.split("\t")[7]
			description = line.strip(";").split("\t")[8]
			#newContig = renameDict[contig]
			if feature != "gene":
				pass
				#splitDescription = description.split("; ")
				#for item in splitDescription:
				#	splitSplitDescription = item.split(" ")
				#	if splitSplitDescription[0] == "transcript_id":
				#		transcript_id = splitSplitDescription[1].strip('''"''')
				#	elif splitSplitDescription[1] == "gene_id":
				#		gene_id = splitSplitDescription[1].strip('''"''')
				#if len(genePrefix) >= 1:
				#	newGene = "%s_%s" % (genePrefix, gene_id)
				#	newTranscript = "%s_%s" % (genePrefix, transcript_id)
				#else:
				#	newGene = gene_id
				#	newTranscript = transcript_id
				#outfile.write('''%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\ttranscript_id "%s"; gene_id "%s"\n''' %(newContig, source, feature, startPos, endPos, score, strand, frame, newTranscript, newGene ))
			else:
				gene_id = description
				gtfDict[contig]["startPosList"].append(startPos)
				try:
					gtfDict[contig]["startPosDict"][startPos].append(gene_id)
				except:
					gtfDict[contig]["startPosDict"].update({startPos : gene_id})

				#if len(genePrefix) >= 1:
				#	newGene = "%s_%s" % (genePrefix, gene_id)
				#else:
				#	newGene = gene_id

				#outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(newContig, source, feature, startPos, endPos, score, strand, frame, newGene ))

			#transcript = transcriptSplit2.strip("\"")
			#gtfDict[transcript] = {"contig":contig, "startPos":startPos, "endPos":endPos, "direction":direction}

		else:
			pass
		#outfile.write("%s\n" % line)
for contig in gtfDict:
	sortedStartPos = sorted(gtfDict[contig]["startPosList"])
	for item in sortedStartPos:
		gtfDict[contig]["sortedGeneList"].append(gtfDict[contig]["startPosDict"][item])
counter = 0
conversionDict = {}
for contig in gtfDict:
	for sorted_gene in gtfDict[contig]["sortedGeneList"]:
		counter = interval + counter
		formattedCounter = '{num:0{width}}'.format(num=counter, width=digits)
		if len(genePrefix) >= 1:
			newGene = "%s%s" % (genePrefix, formattedCounter)
		else:
			newGene = formattedCounter
		conversionDict.update({sorted_gene : newGene})
with open(gtf, "r") as ingtf:
	for line in ingtf:
		line = line.strip("\n")
		if not line.startswith("#"):
			contig = line.split("\t")[0]
			#if contig not in gtfDict.keys():
			#	gtfDict.update({contig : {"geneList" : []}})
			#	gtfDict.update({contig : {"prevStartPos" : "FIRST"}})
				#gtfDict.update({contig : {"maxStartPos" : 0}})
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
					splitSplitDescription = item.split(" ")
					if splitSplitDescription[0] == "transcript_id":
						transcript_id = splitSplitDescription[1].strip('''"''')
						transcript_number = transcript_id.split(".")[1]
					elif splitSplitDescription[1] == "gene_id":
						gene_id = splitSplitDescription[1].strip('''"''')

				#newGene = gtfDict["conversionDict"][gene_id]
				newGene = conversionDict[gene_id]
				#newTranscript = "%s.%s" % (gtfDict["conversionDict"][gene_id], transcript_number)
				newTranscript = "%s.%s" % (conversionDict[gene_id], transcript_number)

				outfile.write('''%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\ttranscript_id "%s"; gene_id "%s"\n''' %(newContig, source, feature, startPos, endPos, score, strand, frame, newTranscript, newGene ))
			else:
				gene_id = description
				#newGene = gtfDict["conversionDict"][gene_id]
				newGene = conversionDict[gene_id]
				outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(newContig, source, feature, startPos, endPos, score, strand, frame, newGene ))
		outfile.write("%s\n" % line)
renameTable.close()
outfile.close()
