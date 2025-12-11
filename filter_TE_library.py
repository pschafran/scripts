#! /usr/bin/env python


# PWS 2025 12 11

from Bio import SeqIO
import sys
import re
import subprocess

usage='''
filter_TE_library.py mod.EDTA.TElib.fa
'''

def parse_fasta(infile):
	AllSeq = SeqIO.parse(infile, 'fasta')
	return [i for i in AllSeq]

# Default variables
tpases_blastdb = '/media/data/resources/Tpases020812'
plantprot_blastdb = '/media/data/resources/blast_dbs/uniprot_swissprot_plants_2025-12-11.fasta'
#plantprot_blastdb = '/media/data/resources/blast_dbs/uniprot_viridiplantae_reference_proteomes/uniprot_viridiplantae_reference_proteomes'
blastx_evalue = '1e-10'
blastx_numdescriptions = '10'
blastx_numthreads = '12'
blastx_outfmt = '6'
blastx_maxtargetseqs = '1'


# Stage 1 - Separate known and unknown LTRs from TE library
infile = sys.argv[1]

input = parse_fasta(infile)
TElibDict = SeqIO.to_dict(SeqIO.parse(infile, "fasta"))

print("Input TE sequences: %s" % len(TElibDict.keys()))

unknown_TE_list= []
known_TE_list= []

output_unknown = open('%s.LTRlib.unknown.fa' % infile, 'w')
output_known = open('%s.LTRlib.known.fa' % infile, 'w')

for seq in input:
	if str(seq.id).find('nknown') != -1:
		output_unknown.write('>' + str(seq.id) + '\n' + str(seq.seq) + '\n')
		unknown_TE_list.append(str(seq.id))
	else:
		output_known.write('>' + str(seq.id) + '\n' + str(seq.seq) + '\n')
		known_TE_list.append(str(seq.id)) 

print("\tKnown TEs: %s" % len(known_TE_list))
print("\tUnknown TEs: %s" % len(unknown_TE_list))

# Stage 2 - BLAST unknown LTRs against a database of transposases
blastx_query = '%s.LTRlib.unknown.fa' % infile
blastx_outfile = '%s.LTRlib.unknown.blastx.out' % infile
blastx_cmd = 'blastx -query %s -db %s -evalue %s -max_hsps 1 -max_target_seqs %s -num_threads %s -outfmt %s' %(blastx_query, tpases_blastdb, blastx_evalue, blastx_maxtargetseqs, blastx_numthreads, blastx_outfmt)
blastx_results = subprocess.run(blastx_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True, shell=True)

# Stage 3 - Extract the positive hits from BLAST results
identified_hits = '%s.LTRlib.unknown.identified.fa' % infile
unknown_hits = '%s.LTRlib.unknown.unknown.fa' % infile
unknown_fasta = SeqIO.to_dict(SeqIO.parse(blastx_query, "fasta"))
source_hits = []
with open(identified_hits, "w") as identified_outfile, open(unknown_hits, "w") as unknown_outfile:
	for line in blastx_results.stdout.splitlines():
		query_id = line.split("\t")[0]
		hit_id = line.split("\t")[1]
		if query_id not in source_hits:
			#identified_outfile.write(">%s\n%s\n" %(query_id, unknown_fasta[query_id].seq))
			try:
				unknown_TE_list.remove(str(query_id))
			except:
				print("Can't remove %s from unknown list" % query_id)
			try:
				known_TE_list.append(str(query_id))
			except:
				print("Can't add %s to known list" % query_id)
		source_hits.append(query_id)
	#for seq in unknown_fasta:
	#	if seq not in source_hits:
	#		unknown_outfile.write(">%s\n%s\n" %(unknown_fasta[seq].id, unknown_fasta[seq].seq))
known_TE_list = set(known_TE_list)
seq_iterator = SeqIO.parse(infile,"fasta")
known_TE_records = [record for record in seq_iterator if record.id in known_TE_list ]
try:
	SeqIO.write(known_TE_records,"%s.LTRlib.known.2.fa" % infile,"fasta")
except:
	with open("%s.LTRlib.known.2.fa" % infile,"w") as outfile:
		for record in known_TE_records:
			print(TElibDict[record].id)
			print(TElibDict[record].seq)
			outfile.write(">%s\n" % TElibDict[record].id)
			outfile.write(TElibDict[record].seq)
	exit(1)

# Stats
print("\nAfter BLASTing transposases:")
print("\tKnown TEs:%s" %(len(set(known_TE_list))))
print("\tUnknown TEs:%s" %(len(set(unknown_TE_list))))

# Stage 4 - BLAST known LTRs against a database of plant proteins
#stage4_outfile = '%s.LTRlib.known.blastx.out' % infile
blastx_query = '%s.LTRlib.known.2.fa' % infile
blastx_cmd ='blastx -query %s -db %s -evalue %s -max_hsps 1 -max_target_seqs %s -num_threads %s -outfmt %s' %(blastx_query, plantprot_blastdb, blastx_evalue, blastx_maxtargetseqs, blastx_numthreads, blastx_outfmt)
blastx_results = subprocess.run(blastx_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True, shell=True)

# Stage 5 - How to replace ProtExcluder package?
# Need to parse blast output to get the alignment lines between each query-subject pair.
# Calculate the amino acid complexity of each alignment, remove low complexity matches
# For now do a filter of any BLAST hits to UniProtRef
known_fasta = SeqIO.to_dict(SeqIO.parse(blastx_query,"fasta"))
for line in blastx_results.stdout.splitlines():
	query_id = line.split("\t")[0]
	hit_id = line.split("\t")[1]
	try:
		known_TE_list.remove(str(query_id))
	except:
		print("Can't remove %s from known list" % query_id)
	try:
		unknown_TE_list.append(str(query_id))
	except:
		print("Can't add %s to unknown list" % query_id)

known_TE_list = set(known_TE_list)
seq_iterator = SeqIO.parse(infile,"fasta")
known_TE_records = [record for record in seq_iterator if record.id in known_TE_list ]
try:
        SeqIO.write(known_TE_records,"%s.LTRlib.known.final.fa" % infile,"fasta")
except:
	print("Error writing LTRlib.known.final.fa")

print("\nAfter BLASTing UniProtRef plant proteomes:")
print("\tKnown TEs:%s" %(len(set(known_TE_list))))
print("\tUnknown TEs:%s" %(len(set(unknown_TE_list))))
