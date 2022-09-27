#! /usr/bin/env python

dict = {}
with open("assembly_info.txt","r") as infile:
	for line in infile:
		if line.startswith("#"):
			pass
		else:
			splitline = line.strip("\n").split("\t")
			dict[splitline[0]] = {"length" : splitline[1]}
			dict[splitline[0]].update({"coverage" : splitline[2]})
for contig in dict:
	abundance = int(dict[contig]["length"]) * int(dict[contig]["coverage"])
	dict[contig].update({"abundance" : abundance})
with open("assembly_info_abundance.txt","w") as outfile:
	outfile.write("#seq_name\tlength\tcoverage\tabundance\n")
	for contig in dict:
		outfile.write("%s\t%s\t%s\t%s\n" %(contig, dict[contig]["length"], dict[contig]["coverage"], dict[contig]["abundance"]))
