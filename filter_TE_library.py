#! /usr/bin/env python


# PWS 2025 12 11

from Bio import SeqIO
import sys
import re


usage='''
filter_TE_library.py mod.EDTA.TElib.fa
'''

def parse_fasta(infile):
	AllSeq = SeqIO.parse(infile, 'fasta')
	return [i for i in AllSeq]


# Default variables
tpases_blastdb = '/media/data/resources/Tpases020812'
blastx_evalue = '1e-10'
blastx_numdescriptions = '10'
blastx_numthreads = '12'
blastx_outfmt = '6'
blastx_maxtargetseqs = '1'


# Stage 1 - Separate known and unknown LTRs from TE library
infile = sys.argv[1]

input = parse_fasta(infile)
output_unknown = open('%s.LTRlib.unknown.fa' % infile, 'w')
output_known = open('%s.LTRlib.known.fa' % infile, 'w')

for seq in input:
	if str(seq.id).find('nknown') != -1:
		output_unknown.write('>' + str(seq.id) + '\n' + str(seq.seq) + '\n')
	else:
		output_known.write('>' + str(seq.id) + '\n' + str(seq.seq) + '\n')

# Stage 2 - BLAST unknown LTRs against a database of transposases

blastx_query = '%s.LTRlib.unknown.fa' % infile
blastx_outfile = '%s.LTRlib.unknown.blastx.out' % infile
blastx_cmd = 'blastx -query %s -db %s -evalue %s -max_target_seqs %s -num_threads %s -outfmt %s' %(blastx_query, tpases_blastdb, blastx_evalue, blastx_maxtargetseqs, blastx_numthreads, blastx_outfmt)
blastx_results = subprocess.run(blatsx_cmd, stdout=subprocess.PIPE, universal_newlines=True)

# Stage 3 - Extract the positive hits from BLAST results

identified_hits = '%s.LTRlib.unknown.identified.fa' % infile
unknown_hits = '%s.LTRlib.unknown.unknown.fa' % infile

unknown_fasta = parse_fasta(blastx_query)
source_hits = []
with open(identified_hits, "w") as identified_outfile, open(unknown_hits, "w") as unknown_outfile:
	for line in blastx_results:
		query_id = line.split("\t")[0]
		hit_id = line.split("\t")[1]
		source_hits.append(query_id)
		identified_outfile.write(">%s\n%s\n" %(hit_id, unknown_fasta[query_id].seq)
	for unknown_fasta.id not in source_hits:
		unknown_outfile.write(">%s\n%s\n") %(unknown_fasta.id, unknown_fasta.seq)
