#! /usr/bin/env python

# Usage: filterGenes.py augustus.hints.gtf score
# Score is any number between 0 and 1.
# If all transcripts for a gene fall below the score,the gene is removed from the output.

import sys, re, math

def geomean(xs):
    return math.exp(math.fsum(math.log(x) for x in xs) / len(xs))

gtfFile = sys.argv[1]
gtfFilename = ".".join(gtfFile.split(".")[0:-1])
score = float(sys.argv[2])
if score < 0 or score > 1:
	print("ERROR: Score must be between 0 and 1")

mainDict = {}
with open(gtfFile, "r") as infile:
	for line in infile:
		if not line.startswith("#"):
			splitline = line.strip("\n").split("\t")
			if splitline[2] == "gene":
				mainDict[splitline[8]] = {"gene_gtf" : line, "transcripts" : {}}
			else:
				gene_id_results = re.search('''gene_id "(.*?)"''',line)
				geneID = gene_id_results.group(1)
				transcript_id_results = re.search('''transcript_id "(.*?)"''',line)
				transcriptID = transcript_id_results.group(1)
				if transcriptID not in mainDict[geneID]["transcripts"].keys():
					mainDict[geneID]["transcripts"].update({transcriptID : {"gtf" : [], "scores" : []}})
				mainDict[geneID]["transcripts"][transcriptID]["gtf"].append(line)
				if splitline[2] == "CDS":
					mainDict[geneID]["transcripts"][transcriptID]["scores"].append(float(splitline[5]))

transcriptsAboveScore = []
transcriptsBelowScore = []
genesBelowScore = []
for gene in mainDict:
	allFail = True
	for transcript in mainDict[gene]["transcripts"]:
		if geomean(mainDict[gene]["transcripts"][transcript]["scores"]) < score:
			transcriptsBelowScore.append(transcript)
		else:
			transcriptsAboveScore.append(transcript)
			allFail = False
	if allFail == True:
		genesBelowScore.append(gene)

if len(transcriptsBelowScore) == 0:
	print("All transcripts passed. No output generated.\n")
	exit(0)

with open("%s_score%s.gtf" %(gtfFilename,score), "w") as outfile:
	#outfile.write("# %s\n" % " ".join(sys.argv))
	for gene in mainDict:
		if gene not in genesBelowScore:
			outfile.write(mainDict[gene]["gene_gtf"])
			#print("%s" % mainDict[gene]["gene_gtf"])
		for transcript in mainDict[gene]["transcripts"]:
			if transcript not in transcriptsBelowScore:
				for gtfline in mainDict[gene]["transcripts"][transcript]["gtf"]:
					outfile.write(gtfline)
					#print("%s" % gtfline)
with open("failed_transcripts_%s.txt" % score, "w") as ftOutfile:
	for transcript in transcriptsBelowScore:
		ftOutfile.write("%s\n" % transcript)
with open("failed_genes_%s.txt" % score, "w") as gnOutfile:
	for gene in genesBelowScore:
		gnOutfile.write("%s\n" % gene)

print("Transcripts passed:\t%s" % len(transcriptsAboveScore))
print("Transcripts failed:\t%s" % len(transcriptsBelowScore))
print("Genes passed:\t%s" % (len(mainDict.keys())-len(genesBelowScore)))
print("Genes failed:\t%s" % len(genesBelowScore))
print("Failed transcripts written to failed_transcripts_%s.txt" % score)
print("Failed genes written to failed_genes_%s.txt" % score)
print("GTF without failed transcripts and genes written to %s_%s.gtf" %(gtfFilename,score))
