#! /bin/bash

# Define functions
command_exists ()
{
	command -v $1
}

sort_and_rename_assembly()
{
	python3 /home/ps997/scripts/getFastaSeqLengths.py $1
	if [[ ! -s "$1"_sequence_lengths.tmp ]]
		then warning_report "Assembly could not be sorted. Check that assembly.fasta exists."
	fi
	sort -k2,2nr "$1"_sequence_lengths.tmp > "$1"_sequence_lengths.sorted.tmp
	awk -F"\t" -v prefix=$CONTIG_PREFIX -v hash=$uniqueMark -v i="1" '{ print $1"\t"prefix""i++" "hash }' "$1"_sequence_lengths.sorted.tmp > "$1"_sequence_lengths.new_contig_names.tsv
	python3 /home/ps997/scripts/renameFastaAndReorder.py "$1" "$1"_sequence_lengths.new_contig_names.tsv
}

error_exit()
{
	printf "\nERROR: %s\n\n" "$1"
	exit 1
}

warning_report()
{
	printf "\nWARNING: %s\n\n" "$1"
}

error_break()
{
	printf "\nERROR: %s\n\n" "$1"
	break
}

# Read in config file
source $(readlink -e $1)

# Need to source conda in subshell to be able to activate environments
if [ -e $CONDA_PATH ]
	then source $CONDA_PATH
elif [ -e ~/miniconda3/etc/profile.d/conda.sh ]
	then source ~/miniconda3/etc/profile.d/conda.sh
elif [ -e ~/anaconda3/etc/profile.d/conda.sh ]
	then source ~/anaconda3/etc/profile.d/conda.sh
else
	echo "ERROR: Can't initialize conda environment. Check that ~/miniconda3/etc/profile.d/conda.sh or ~/anaconda3/etc/profile.d/conda.sh or the user-supplied path exists."
	exit 1
fi
conda activate base

# Check input files exist and convert relative paths into absolute
if [[ -e $NANOPORE ]]
	then NANOPORE=$(readlink -f $NANOPORE)
else
	error_exit "$NANOPORE file not found"
fi

if [[ -e $ILLUMINA1 ]]
	then ILLUMINA1=$(readlink -f $ILLUMINA1)
else
	warning_report "$ILLUMINA1 file not found"
fi

if [[ -e $ILLUMINA2 ]]
	then ILLUMINA2=$(readlink -f $ILLUMINA2)
else
	warning_report "$ILLUMINA2 file not found"
fi

if [[ -e $RNA1 ]]
	then RNA1=$(readlink -f $RNA1)
else
	warning_report "$RNA1 file not found"
fi

if [[ -e $RNA2 ]]
	then RNA2=$(readlink -f $RNA2)
else
	warning_report "$RNA2 file not found"
fi

if [[ -e $REFPROT ]]
	then REFPROT=$(readlink -f $REFPROT)
else
	warning_report "$REFPROT" file not found
fi

if [[ -e $PLASTOME ]]
	then PLASTOME=$(readlink -f $PLASTOME)
else
	warning_report "$PLASTOME file not found"
fi

if [[ -e $CHONDROME ]]
	then CHONDROME=$(readlink -f $CHONDROME)
else
	warning_report "$CHONDROME file not found"
fi
# Check dependencies
for i in $ASSEMBLER
	do echo "$i"
 	if [ "$i" == "flye" ]
		then if [ ! -x $(command_exists "flye") ]
			then error_exit "Flye couldn't be executed"
		fi
	elif [ "$i" == "miniasm" ]
		then if [ ! -x $(command_exists "minimap2") ]
			then error_exit "Minimap2 couldn't be executed"
		elif [ ! -x $(command_exists "miniasm") ]
			then error_exit "Miniasm couldn't be executed"
		elif [ ! -x $(command_exists "minipolish") ]
			then error_exit "Minipolish couldn't be executed"
		fi
	elif [ "$i" == "raven" ]
		then if [ ! -x $(command_exists "raven") ]
			then error_exit "Raven couldn't be executed"
		fi
	elif [ "$i" == "haslr" ]
		then if [ ! -x $(command_exists haslr.py) ]
			then error_exit "HASLR couldn't be executed"
		fi
	elif [ "$i" == "canu" ]
		then if [ ! -x $(command_exists $CANU_PATH) ]
			then error_exit "CANU counldn't be executed. Did you remember to supply a path in the config file?"
		fi
	elif [ "$i" == "wtdbg2" ]
		then if [ ! -x $(command_exists wtdbg2) ]
			then error_exit "Wtdbg2 counldn't be executed"
		elif [ ! -x $(command_exists minimap2) ]
			then error_exit "Minimap2 couldn't be executed"
		elif [ ! -x $(command_exists samtools) ]
			then error_exit "Samtools couldn't be executed"
		fi
	elif [ "$i" == "masurca-flye" ] || [ "$i" == "masurca-cabog" ]
		then if [ ! -x $(command_exists $MASURCA_PATH) ]
			then error_exit "Masurca couldn't be executed. Did you remember to supply a path in the config file?"
		fi
	fi
done

if [ ! -x $(command_exists $REPEATMASKER_PATH) ]
	then error_exit "RepeatMasker couldn't be executed. Did you remember to supply a path in the config file?"
fi




echo ""
echo "***** RUN CONDITIONS *****"
echo "Sample Name: $SAMPLE_NAME"
echo "Contig prefix: $CONTIG_PREFIX"
printf "Genome Size: %s\n" $GENOME_SIZE
echo "Assemblers: $ASSEMBLER"
echo "Assembly selection based on: $SELECT_CRITERION"
echo "BUSCO lineage: $BUSCO_DB"
echo ""
echo "***** INPUT FILES *****"
printf "NANOPORE file:\t%s\n" $NANOPORE
printf "ILLUMINA files:\t%s\n" $ILLUMINA1
printf "\t\t%s\n" $ILLUMINA2
printf "RNA files:\t%s\n" $RNA1
printf "\t\t%s\n" $RNA2
echo "Reference protein file: $REFPROT"
echo "Chloroplast genome file: $PLASTOME"
echo "Mitocondria genome file: $CHONDROME"
echo ""
echo "***** DEPENDENCY PATHS *****"
echo "MASURCA_PATH: $MASURCA_PATH"
echo "AUGUSTUS_CONFIG_PATH: $AUGUSTUS_CONFIG_PATH"
echo "AUGUSTUS_SCRIPTS_PATH: $AUGUSTUS_SCRIPTS_PATH"
echo "AUGUSTUS_BIN_PATH: $AUGUSTUS_BIN_PATH"
echo "GENEMARK_PATH: $GENEMARK_PATH"
printf "CONDA_PATH:\t%s\n" $CONDA_PATH


# Start script in ~/$SAMPLE_NAME/assemblies/ directory
# Expects a directory structure like this:
#
# Creates this directory structure:
# ./assemblies/
#	|./flye/
#	|./miniasm/
#	|./masurca-flye/
#	|./masurca-cabog/
#	|./raven/
#	|
#	|./repeat-masking/
#	|	|./flye/
#	|	|	|./0_RemoveOrganelles/
#	|	|	|./1a_RepeatModeler/
#	|	|	|./1b_EDTA/
#	|	|	|./2_ProteinExclude/
#	|	|	|./3_RepeatMasker
#	|	|
#	|	|./masurca-cabog/
#	|	|	|./0_RemoveOrganelles/
#	|	|	|./1a_RepeatModeler/
#	|	|	|./1b_EDTA/
#	|	|	|./2_ProteinExclude/
#	|	|	|./3_RepeatMasker
#	|	|
#	|	|./masurca-flye/
#	|	|	|./0_RemoveOrganelles/
#	|	|	|./1a_RepeatModeler/
#	|	|	|./1b_EDTA/
#	|	|	|./2_ProteinExclude/
#	|	|	|./3_RepeatMasker
#	|	|
#	|	|./miniasm/
#	|	|	|./0_RemoveOrganelles/
#	|	|	|./1a_RepeatModeler/
#	|	|	|./1b_EDTA/
#	|	|	|./2_ProteinExclude/
#	|	|	|./3_RepeatMasker
#	|
#	|./gene-prediction/
#	|	|./flye/
#	|	|	|./1a_RepeatModeler/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1b_EDTA/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1c_EDTA_noFilter/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|
#	|	|./masurca-cabog/
#	|	|	|./1a_RepeatModeler/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1b_EDTA/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1c_EDTA_noFilter/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|
#	|	|./masurca-flye/
#	|	|	|./1a_RepeatModeler/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1b_EDTA/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1c_EDTA_noFilter/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|
#	|	|./miniasm/
#	|	|	|./1a_RepeatModeler/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1b_EDTA/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|./1c_EDTA_noFilter/
#	|	|	|	|./braker-RNA/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#	|	|	|	|./braker-RNA+PROT/augustus.hints.[aa|codingseq]  <-- These are probably the final files you want
#		|./$SAMPLE_NAME_combined_proteins.faa
#		|./$SAMPLE_NAME_clustered_proteins.faa <-- All proteins from all assemblies and repeat-maskings clustered with CD-HIT at 0.9 similarity
#


#Cores used by repeat modeler gets inflated 3x so create a new varible 1/3 of user specified number
RM_CORES=$((CORES/3))

# START_DIR holds location where script was executed
START_DIR=$(pwd)

# Create output directory where everything takes place
if [[ ! -d $SAMPLE_NAME ]]
	then mkdir "$SAMPLE_NAME"
fi

# HOMEDIR holds location to base output firectory
HOMEDIR="$START_DIR"/"$SAMPLE_NAME"

# Store the date when script started. Might implement later to indicate when different runs were done
DATE=$(date +%F)

# Remove reads shorter then $MIN_READ_LENGTH
cd "$SAMPLE_NAME"
HOMEDIR=$(pwd)
printf "\n***** REMOVING SHORT ONT READS *****\n"
if [ -s "$HOMEDIR"/nanopore.minlen"$MIN_READ_LENGTH".fastq.gz ]
	then printf "Reads already filtered to correct length. Skipping..."
	NANOPORE="$HOMEDIR"/nanopore.minlen"$MIN_READ_LENGTH".fastq.gz
else
	if [[ $NANOPORE = *.gz ]]
		then awk -v MIN_READ_LENGTH=$MIN_READ_LENGTH 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= MIN_READ_LENGTH ) {print header, seq, qheader, qseq}}' <(zcat $NANOPORE) | gzip > nanopore.minlen"$MIN_READ_LENGTH".fastq.gz
	else
		awk -v MIN_READ_LENGTH=$MIN_READ_LENGTH 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= MIN_READ_LENGTH ) {print header, seq, qheader, qseq}}' < $NANOPORE| gzip > nanopore.minlen"$MIN_READ_LENGTH".fastq.gz
	fi
	NANOPORE="$HOMEDIR"/nanopore.minlen"$MIN_READ_LENGTH".fastq.gz
	printf "Longer reads written to %s\n" $NANOPORE
fi

# Make a hash unique to the input files that will be added to all files produced
uniqueMark=$(tar -cf - $NANOPORE $ILLUMINA1 $ILLUMINA2 $RNA1 $RNA2 $REFPROT $PLASTOME $CHONDROME | md5sum)
echo "Input file hash: $uniqueMark"
# Start assembly stage
printf "\n***** ASSEMBLY *****\n"
if [ ! -d assemblies ]
	then mkdir assemblies
fi
cd assemblies

for i in $ASSEMBLER
	do
		if [ ! -d $i ]
			then mkdir $i
		fi
		cd $i

		if [ $i == "flye" ]
			then if [ -s "$SAMPLE_NAME"_flye_assembly.fasta ]
				then printf "Flye already completed. Skipping...\n"
			else
				if [[ ! -s assembly.fasta ]]
					then printf "Running flye...\n"
					conda activate flye
					flye --nano-raw $NANOPORE -t $CORES -o . 1> flye.out 2> flye.err
					if [ ! -s assembly.fasta ]
						then printf "ERROR: Flye failed to finish. Check flye.log and flye.err\n"
					fi
				else
					printf "Flye raw assembly appears to exist. Will try renaming and sorting...\n"
				fi
				conda activate base
				# rename genome contigs by length
				sort_and_rename_assembly assembly.fasta
				#python3 /home/ps997/scripts/getFastaSeqLengths.py assembly.fasta
				#if [[ ! -s assembly.fasta_sequence_lengths.tmp ]]
				#	then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				#fi
				#sort -k2,2nr assembly.fasta_sequence_lengths.tmp > assembly.fasta_sequence_lengths.sorted.tmp
				#awk -F"\t" -v prefix=$CONTIG_PREFIX -v i="1" '{ print $1"\t"prefix""i++ }' assembly.fasta_sequence_lengths.sorted.tmp > assembly.fasta_sequence_lengths.new_contig_names.tsv
				#python3 /home/ps997/scripts/renameFastaAndReorder.py assembly.fasta assembly.fasta_sequence_lengths.new_contig_names.tsv
				if [[ ! -s assembly_renamed.fasta ]]
					then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				fi
				ln -sf assembly_renamed.fasta "$SAMPLE_NAME"_flye_assembly.fasta
				cd "$HOMEDIR"/assemblies/
				ln -sf "$i"/"$SAMPLE_NAME"_flye_assembly.fasta
				printf "Flye complete.\n"
			fi
		fi

		if [ $i == "miniasm" ]
			then if [ -e "$SAMPLE_NAME"_miniasm_assembly.fasta ]
				then echo "Miniasm already completed. Skipping..."
			else
				printf "Running miniasm..."
				conda activate miniasm
				minimap2 -x ava-ont -t $CORES $NANOPORE $NANOPORE 2> minimap.err | gzip -1 > reads.paf.gz
				miniasm -f $NANOPORE reads.paf.gz > "$SAMPLE_NAME".gfa 2> miniasm.err
				minipolish -t $CORES $NANOPORE "$SAMPLE_NAME".gfa > "$SAMPLE_NAME"_polished.gfa 2> minipolish.err
				awk '/^S/{print ">"$2"\n"$3}' "$SAMPLE_NAME"_polished.gfa | fold > assembly.fasta
				conda activate base
				sort_and_rename_assembly assembly.fasta
				#python3 /home/ps997/scripts/getFastaSeqLengths.py assembly.fasta
				#if [[ ! -s assembly.fasta_sequence_lengths.tmp ]]
				#	then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				#fi
				#sort -k2,2nr assembly.fasta_sequence_lengths.tmp > assembly.fasta_sequence_lengths.sorted.tmp
				#awk -F"\t" -v prefix=$CONTIG_PREFIX -v i="1" '{ print $1"\t"prefix""i++ }' assembly.fasta_sequence_lengths.sorted.tmp > assembly.fasta_sequence_lengths.new_contig_names.tsv
				#python3 /home/ps997/scripts/renameFastaAndReorder.py assembly.fasta assembly.fasta_sequence_lengths.new_contig_names.tsv
				if [[ ! -s assembly_renamed.fasta ]]
					then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				fi
				ln -sf assembly_renamed.fasta "$SAMPLE_NAME"_miniasm_assembly.fasta
				cd "$HOMEDIR"/assemblies/
				ln -sf "$i"/"$SAMPLE_NAME"_miniasm_assembly.fasta
				printf "Miniasm complete.\n"
			fi
		fi

		if [ $i == "raven" ]
			then if [ -e "$SAMPLE_NAME"_raven_assembly.fasta ]
				then echo "Raven already completed. Skipping..."
			else
				printf "Running raven..."
				conda activate raven
				raven -t $CORES $NANOPORE > assembly.fasta 2> raven.err
				conda activate base
				sort_and_rename_assembly assembly.fasta
				#python3 /home/ps997/scripts/getFastaSeqLengths.py assembly.fasta
				#if [[ ! -s assembly.fasta_sequence_lengths.tmp ]]
				#	then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				#fi
				#sort -k2,2nr assembly.fasta_sequence_lengths.tmp > assembly.fasta_sequence_lengths.sorted.tmp
				#awk -F"\t" -v prefix=$CONTIG_PREFIX -v i="1" '{ print $1"\t"prefix""i++ }' assembly.fasta_sequence_lengths.sorted.tmp > assembly.fasta_sequence_lengths.new_contig_names.tsv
				#python3 /home/ps997/scripts/renameFastaAndReorder.py assembly.fasta assembly.fasta_sequence_lengths.new_contig_names.tsv
				if [[ ! -s assembly_renamed.fasta ]]
					then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				fi
				ln -sf assembly_renamed.fasta "$SAMPLE_NAME"_raven_assembly.fasta
				cd "$HOMEDIR"/assemblies/
				ln -sf "$i"/"$SAMPLE_NAME"_raven_assembly.fasta
				printf "Raven complete.\n"
			fi
		fi

		if [ $i == "masurca-flye" ]
			then echo "Running Masurca with flye..."
			conda activate masurca
			printf "DATA\nPE= pe 500 50  %s  %s\n" $ILLUMINA1 $ILLUMINA2 >> masurca_config.txt
			printf "NANOPORE=%s\nEND\n\n" $NANOPORE >> masurca_config.txt
			printf "PARAMETERS\nEXTEND_JUMP_READS=0\nGRAPH_KMER_SIZE = auto\nUSE_LINKING_MATES = 0\nUSE_GRID=0\nGRID_ENGINE=SGE\nGRID_QUEUE=all.q\nGRID_BATCH_SIZE=500000000\nLHE_COVERAGE=25\nLIMIT_JUMP_COVERAGE = 300\nCA_PARAMETERS =  cgwErrorRate=0.15\nCLOSE_GAPS=1\nNUM_THREADS = %s\nJF_SIZE = 200000000\nSOAP_ASSEMBLY=0\nFLYE_ASSEMBLY=1\nEND\n" $CORES >> masurca_config.txt
			"$MASURCA_PATH" masurca_config.txt
			./assemble.sh > masurca.out 2> masurca.err
			echo "Masurca-Flye complete."
			conda activate base
		fi

		if [ $i == "masurca-cabog" ]
			then echo "Running Masurca with CABOG..."
			conda activate masurca
			printf "DATA\nPE= pe 500 50  %s  %s\n" $ILLUMINA1 $ILLUMINA2 >> masurca_config.txt
			printf "NANOPORE=%s\nEND\n\n" $NANOPORE >> masurca_config.txt
			printf "PARAMETERS\nEXTEND_JUMP_READS=0\nGRAPH_KMER_SIZE = auto\nUSE_LINKING_MATES = 0\nUSE_GRID=0\nGRID_ENGINE=SGE\nGRID_QUEUE=all.q\nGRID_BATCH_SIZE=500000000\nLHE_COVERAGE=25\nLIMIT_JUMP_COVERAGE = 300\nCA_PARAMETERS =  cgwErrorRate=0.15\nCLOSE_GAPS=1\nNUM_THREADS = %s\nJF_SIZE = 200000000\nSOAP_ASSEMBLY=0\nFLYE_ASSEMBLY=0\nEND\n" $CORES >> masurca_config.txt
			"$MASURCA_PATH" masurca_config.txt
			./assemble.sh > masurca.out 2> masurca.err
			echo "Masurca-CABOG complete."
			conda activate base
		fi

		if [ $i == "wtdbg2" ]
			then if [ -e "$SAMPLE_NAME"_wtdbg2_assembly.fasta ]
				then printf "Wtdbg2 already completed. Skipping...\n"
			else
				printf "Running wtdbg2..."
				conda activate wtdbg2
				wtdbg2 -x ont -i $NANOPORE -t $CORES -fo dbg > wtdbg2.out 2> wtdbg2.err
				wtpoa-cns -t $CORES -i dbg.ctg.lay.gz -fo dbg.raw.fa > wtpoa-cns.out 2> wtpoa-cns.err
				minimap2 -t $CORES -ax map-ont -r2k dbg.raw.fa reads.fa.gz 2> minimap2.err | samtools sort -@4 >dbg.bam
				samtools view -F0x900 dbg.bam | wtpoa-cns -t $CORES -d dbg.raw.fa -i - -fo assembly.fasta > wtpoa-cns.out 2> wtpoa-cns.err
				conda activate base
				sort_and_rename_assembly assembly.fasta
				#python3 /home/ps997/scripts/getFastaSeqLengths.py assembly.fasta
				#if [[ ! -s assembly.fasta_sequence_lengths.tmp ]]
				#	then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				#fi
				#sort -k2,2nr assembly.fasta_sequence_lengths.tmp > assembly.fasta_sequence_lengths.sorted.tmp
				#awk -F"\t" -v prefix=$CONTIG_PREFIX -v i="1" '{ print $1"\t"prefix""i++ }' assembly.fasta_sequence_lengths.sorted.tmp > assembly.fasta_sequence_lengths.new_contig_names.tsv
				#python3 /home/ps997/scripts/renameFastaAndReorder.py assembly.fasta assembly.fasta_sequence_lengths.new_contig_names.tsv
				if [[ ! -s assembly_renamed.fasta ]]
					then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				fi
				ln -sf assembly_renamed.fasta "$SAMPLE_NAME"_wtdbg2_assembly.fasta
				cd "$HOMEDIR"/assemblies/
				ln -sf "$i"/"$SAMPLE_NAME"_wtdbg2_assembly.fasta
				# do something
				printf "Wtdbg2 complete.\n"
			fi
		fi

		if [ $i == "haslr" ]
			then if [ -e "$SAMPLE_NAME"_haslr_assembly.fasta ]
				then printf "HASLR already complete. Skipping...\n"
			else
				printf "Running HASLR..."
				conda activate haslr
				haslr.py -t $CORES -g $GENOME_SIZE -o . -l $NANOPORE -x nanopore -s $ILLUMINA1 $ILLUMINA2 1> haslr.out 2> haslr.err
				ln -sf asm_contigs_k49_a3_lr25x_b500_s3_sim0.85/asm.final.fa assembly.fasta
				conda activate base
				sort_and_rename_assembly assembly.fasta
				#python3 /home/ps997/scripts/getFastaSeqLengths.py assembly.fasta
				#if [[ ! -s assembly.fasta_sequence_lengths.tmp ]]
				#	then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				#fi
				#sort -k2,2nr assembly.fasta_sequence_lengths.tmp > assembly.fasta_sequence_lengths.sorted.tmp
				#awk -F"\t" -v prefix=$CONTIG_PREFIX -v i="1" '{ print $1"\t"prefix""i++ }' assembly.fasta_sequence_lengths.sorted.tmp > assembly.fasta_sequence_lengths.new_contig_names.tsv
				#python3 /home/ps997/scripts/renameFastaAndReorder.py assembly.fasta assembly.fasta_sequence_lengths.new_contig_names.tsv
				if [[ ! -s assembly_renamed.fasta ]]
					then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				fi
				ln -sf assembly_renamed.fasta "$SAMPLE_NAME"_haslr_assembly.fasta
				cd "$HOMEDIR"/assemblies/
				ln -sf "$i"/"$SAMPLE_NAME"_haslr_assembly.fasta
				printf "HASLR complete.\n"
			fi
		fi

		if [ $i == "canu" ]
			then if [ -e "$SAMPLE_NAME"_canu_assembly.fasta ]
				then printf "CANU already complete. Skipping...\n"
			else
				printf "Running CANU..."
				$CANU_PATH -p $SAMPLE_NAME genomeSize="$GENOME_SIZE" maxThreads="$CORES" -nanopore $NANOPORE 1> canu.out 2> canu.err
				ln -sf "$SAMPLE_NAME".contigs.fasta assembly.fasta
				sort_and_rename_assembly assembly.fasta
				#python3 /home/ps997/scripts/getFastaSeqLengths.py assembly.fasta
				#sort -k2,2nr assembly.fasta_sequence_lengths.tmp > assembly.fasta_sequence_lengths.sorted.tmp
				#awk -F"\t" -v prefix=$CONTIG_PREFIX -v i="1" '{ print $1"\t"prefix""i++ }' assembly.fasta_sequence_lengths.sorted.tmp > assembly.fasta_sequence_lengths.new_contig_names.tsv
				#python3 /home/ps997/scripts/renameFastaAndReorder.py assembly.fasta assembly.fasta_sequence_lengths.new_contig_names.tsv
				if [[ ! -s assembly_renamed.fasta ]]
					then error_break "Assembly could not be sorted. Check that assembly.fasta exists."
				fi
				ln -sf assembly_renamed.fasta "$SAMPLE_NAME"_canu_assembly.fasta
				cd "$HOMEDIR"/assemblies/
				ln -sf "$i"/"$SAMPLE_NAME"_canu_assembly.fasta
				printf "CANU complete.\n"
				conda activate base
			fi
		fi
	cd "$HOMEDIR"/assemblies/
done
cd "$HOMEDIR"/assemblies/

printf "\n***** ASSEMBLY DONE *****\n"
#for i in "$HOMEDIR"/assemblies/*/"$SAMPLE_NAME"*assembly.fasta
#	do if [ -e $i ]
#		then ln -sf $i "$HOMEDIR"/assemblies/
#	fi
#done

printf "\n***** RUNNING BUSCO *****\n"
if [ -e busco_complete.txt ]
	then rm busco_complete.txt
fi
touch busco_complete.txt
conda activate busco
for i in "$SAMPLE_NAME"*assembly.fasta
	do if [ -e $i ]
		then j=$(basename $i)
		if [ -e "$j"_BUSCO/short_summary*txt ]
			then printf "%s BUSCO results already exist. Skipping...\n" $i
			BUSCO_COMPLETE=$(grep "C:" "$j"_BUSCO/short_summary*txt | awk -F":|%" '{print $2}')
			printf "%s\t%s\n" $i $BUSCO_COMPLETE >> busco_complete.txt
		else
			if [ $BUSCO_DB == "AUTO" ]
				then busco -c $CORES -i $i -o "$j"_BUSCO -m genome --auto-lineage &> "$i"_BUSCO.log
				BUSCO_COMPLETE=$(grep "C:" "$j"_BUSCO/short_summary*txt | awk -F":|%" '{print $2}')
				printf "%s\t%s\n" $i $BUSCO_COMPLETE >> busco_complete.txt
			else
				busco -c $CORES -i $i -o "$j"_BUSCO -m genome -l $BUSCO_DB &> "$i"_BUSCO.log
				BUSCO_COMPLETE=$(grep "C:" "$j"_BUSCO/short_summary*txt | awk -F":|%" '{print $2}')
				printf "%s\t%s\n" $i $BUSCO_COMPLETE >> busco_complete.txt
			fi
		fi
	fi
done
conda activate base

# RUN AUN script to calculate area under N curve to determine
python3 /home/ps997/scripts/auN.py *assembly.fasta > area_under_Ncurve.txt
sort -k2,2nr area_under_Ncurve.txt | head -n 1 | cut -f 1 > best_assembly.txt
BEST_ASSEMBLY=$(sort -k2,2nr area_under_Ncurve.txt | head -n 1 | cut -f 1)
BEST_BUSCO=$(sort -k2,2nr busco_complete.txt | head -n 1 | cut -f 1)


ASSEMBLIES_TO_POLISH=()
if [ "$ASSEMBLY_SELECTION_STAGE" = "ASSEMBLY" ]
	then if [ "$SELECT_CRITERION" = "CONTIGUITY" ]
		then j=$(basename $BEST_ASSEMBLY | rev | cut -d"_" -f 2 | rev)
		ASSEMBLIES_TO_POLISH+=("$j")
	elif [ "$SELECT_CRITERION" = "BUSCO" ]
		then j=$(basename $BEST_BUSCO | rev | cut -d"_" -f 2 | rev)
		ASSEMBLIES_TO_POLISH+=("$j")
	else
		if [ -e busco_x_contiguity.txt ]
			then rm busco_x_contiguity.txt
		fi
		touch busco_x_contiguity.txt
		for i in "$SAMPLE_NAME"*assembly.fasta
			do if [ -e $i ]
					then TMP_CONTIGUITY=$(grep $i area_under_Ncurve.txt | cut -f 2)
					TMP_BUSCO=$(grep $i busco_complete.txt | cut -f 2)
					BUSCO_CONTIGUITY=$(printf "%f\t%f" $TMP_CONTIGUITY $TMP_BUSCO | awk -F"\t" '{print 0.0001*$1*$2}')
					printf "%s\t%s\n" $i $BUSCO_CONTIGUITY >> busco_x_contiguity.txt
				fi
			done
			BEST_BUSCO_CONTIGUITY=$(sort -k2,2nr busco_x_contiguity.txt | head -n 1 | cut -f 1)
			j=$(basename $BEST_BUSCO_CONTIGUITY | rev | cut -d"_" -f 2 | rev)
			ASSEMBLIES_TO_POLISH+=("$j")
	printf "Assembly to polish: %s\n" $j
	fi
else
	for i in "$SAMPLE_NAME"*assembly.fasta
		do if [ -e $i ]
			then j=$(basename $i | rev | cut -d"_" -f 2 | rev)
			ASSEMBLIES_TO_POLISH+=("$j")
		fi
	done
	printf "Assemblies to polish: %s\n" "${ASSEMBLIES_TO_POLISH[@]}"
fi

printf "\n***** POLISHING *****\n"
cd "$HOMEDIR"
if [ ! -d polishing ]
	then mkdir polishing
fi
cd polishing
for i in "${ASSEMBLIES_TO_POLISH[@]}"
	do if [ ! -d $i ]
		then mkdir $i
	fi
	if [ -e "$SAMPLE_NAME"_"$i"_polished_assembly.fasta ]
		then printf "%s already polished. Skipping...\n" $i
	else
		cd $i
		printf "%s assembly\n" $i
		ln -sf "$HOMEDIR"/assemblies/"$SAMPLE_NAME"_"$i"_assembly.fasta pilon-iter0.fasta
		PILONCOUNT=1
		while [ $PILONCOUNT -le $POLISHING_ROUNDS ]
			do PRIORCOUNT=$((PILONCOUNT-1))
			printf "\tpolishing round %s\n" $PILONCOUNT
			bwa index pilon-iter"$PRIORCOUNT".fasta &> /dev/null
			bwa mem -t $CORES pilon-iter"$PRIORCOUNT".fasta $ILLUMINA1 $ILLUMINA2 2> /dev/null | samtools sort -o pilon-iter"$PILONCOUNT".illumina.bam &> samtools.out
			minimap2 -ax map-ont -t $CORES pilon-iter"$PRIORCOUNT".fasta $NANOPORE 2> minimap2.err | samtools sort -o pilon-iter"$PILONCOUNT".ont.bam &> samtools.out
			samtools index pilon-iter"$PILONCOUNT".illumina.bam
			samtools index pilon-iter"$PILONCOUNT".ont.bam
			java -Xmx"$JAVAMEM" -jar $PILON_PATH --genome pilon-iter"$PRIORCOUNT".fasta --frags pilon-iter"$PILONCOUNT".illumina.bam --nanopore pilon-iter"$PILONCOUNT".ont.bam --output pilon-iter"$PILONCOUNT" --changes > pilon.out 2> pilon.err
			sed -i 's/_pilon//' pilon-iter"$PILONCOUNT".fasta
			PILONCOUNT=$((PILONCOUNT+1))
		done
		PRIORCOUNT=$((PILONCOUNT-1))
		ln -sf pilon-iter"$PRIORCOUNT".fasta "$SAMPLE_NAME"_"$i"_polished_assembly.fasta
		cd "$HOMEDIR"/polishing/
		ln -sf "$i"/"$SAMPLE_NAME"_"$i"_polished_assembly.fasta
		rm "$i"/*.bam
	fi
	cd "$HOMEDIR"/polishing/
done

printf "\n***** POLISHING DONE *****\n"

for i in "$HOME"/polishing/*/*polished_assembly.fasta
	do if [ -e $i ]
		then ln -sf $i "$HOMEDIR"/polishing/
	fi
done
cd "$HOMEDIR"/polishing/

printf "\n***** RUNNING BUSCO *****\n"
conda activate busco
for i in *polished_assembly.fasta
	do if [ -e $i ]
		then j=$(basename $i)
		if [ -e "$j"_BUSCO/short_summary*txt ]
			then printf "%s BUSCO results already exist. Skipping...\n" $i
		else
			if [ $BUSCO_DB == "AUTO" ]
				then busco -c $CORES -i $i -o "$j"_BUSCO -m genome --auto-lineage
				BUSCO_COMPLETE=$(grep "C:" "$j"_BUSCO/short_summary*txt | awk -F":|%" '{print $2}')
				printf "%s\t%s\n" $i $BUSCO_COMPLETE >> busco_complete.txt
			else
				busco -c $CORES -i $i -o "$j"_BUSCO -m genome -l $BUSCO_DB
				BUSCO_COMPLETE=$(grep "C:" "$j"_BUSCO/short_summary*txt | awk -F":|%" '{print $2}')
				printf "%s\t%s\n" $i $BUSCO_COMPLETE >> busco_complete.txt
			fi
		fi
	fi
done
conda activate base

# RUN AUN script to calculate area under N curve to determine
python3 /home/ps997/scripts/auN.py *polished_assembly.fasta > area_under_Ncurve.txt
sort -k2,2nr area_under_Ncurve.txt | head -n 1 | cut -f 1 > best_assembly.txt
BEST_ASSEMBLY=$(sort -k2,2nr area_under_Ncurve.txt | head -n 1 | cut -f 1)
BEST_BUSCO=$(sort -k2,2nr busco_complete.txt | head -n 1 | cut -f 1)

ASSEMBLIES_TO_MASK=()
if [ "$ASSEMBLY_SELECTION_STAGE" = "POLISHING" ]
	then if [ "$SELECT_CRITERION" = "CONTIGUITY" ]
		then j=$(basename $BEST_ASSEMBLY | rev | cut -d"_" -f 3 | rev)
		ASSEMBLIES_TO_MASK+=("$j")
		printf "Assembly to mask: %s\n" $j
	elif [ "$SELECT_CRITERION" = "BUSCO" ]
		then j=$(basename $BEST_BUSCO | rev | cut -d"_" -f 3 | rev)
		ASSEMBLIES_TO_MASK+=("$j")
		printf "Assembly to mask: %s\n" $j
	else [ "$SELECT_CRITERION" = "BUSCO*CONTIGUITY" ]
		for i in *polished_assembly.fasta
			do if [ -e $i ]
				then TMP_CONTIGUITY=$(grep $i area_under_Ncurve.txt | cut -f 2)
				TMP_BUSCO=$(grep $i busco_complete.txt | cut -f 2)
				BUSCO_CONTIGUITY=$(awk '{print $1*$2}' <<< "${TMP_CONTIGUITY} ${TMP_BUSCO}")
				printf "%s\t%s\n" $i $BUSCO_CONTIGUITY >> busco_x_contiguity.txt
			fi
		done
		BEST_BUSCO_CONTIGUITY=$(sort -k2,2nr busco_x_contiguity.txt | head -n 1 | cut -f 1)
		j=$(basename $BEST_BUSCO_CONTIGUITY | rev | cut -d"_" -f 3 | rev)
		ASSEMBLIES_TO_MASK+=("$j")
		printf "Assembly to mask: %s\n" $j
	fi
else
	for i in *polished_assembly.fasta
		do if [[ -e $i ]]
			then j=$(basename $i | rev | cut -d"_" -f 3 | rev)
			ASSEMBLIES_TO_MASK+=("$j")
		fi
	done
	printf "Assemblies tow mask: %s\n" "${ASSEMBLIES_TO_POLISH[@]}"
fi

printf "\n***** REPEAT MASKING *****\n"
cd $HOMEDIR
if [[ ! -d repeat-masking ]]
	then mkdir repeat-masking
fi
cd repeat-masking
if [[ -s $PLASTOME ]]
 	then if [[ -s $CHONDROME ]]
		then cat $PLASTOME $CHONDROME > organelle_genomes.fa
	else
		cat $PLASTOME > organelle_genomes.fa
	fi
elif [[ -s $CHONDROME ]]
	then cat $CHONDROME > organelle_genomes.fa
else
	printf "No organelle sequences supplied, will not be removed from assembly."
fi

if [[ -s organelle_genomes.fa ]]
	then makeblastdb -in organelle_genomes.fa -dbtype nucl &> makeblastdb.out
fi

for i in "${ASSEMBLIES_TO_MASK[@]}"
	do if [[ ! -d $i ]]
		then mkdir $i
	fi
	cd $i
	if [ ! -d 0_RemoveOrganelles ]
		then mkdir 0_RemoveOrganelles
	fi
	if [ ! -d 1a_RepeatModeler ]
		then mkdir 1a_RepeatModeler
	fi
	if [ ! -d 1b_EDTA ]
		then mkdir 1b_EDTA
	fi
	if [ ! -d 2_ProteinExclude ]
		then mkdir 2_ProteinExclude
	fi
	if [ ! -d 3_RepeatMasker ]
		then mkdir 3_RepeatMasker
	fi

	if [ -s organelle_genomes.fa ]
		then cd "$HOMEDIR"/repeat-masking/"$i"/0_RemoveOrganelles
		if [ ! -e 0_assembly.fa ]
			then ln -sf "$HOMEDIR"/polishing/"$i"/"$SAMPLE_NAME"_"$i"_polished_assembly.fasta 0_assembly.fa
		fi
		#sed -i --follow-symlinks 's/_pilon_pilon_pilon_pilon_pilon//' 0_assembly.fa # Should be fixed with change in pilon loop
		if [ ! -s organelles.blast ]
			then blastn -query 0_assembly.fa -db "$HOMEDIR"/repeat-masking/organelle_genomes.fa -evalue 0.0001 -outfmt "6 qseqid sseqid pident length qcovs bitscore qstart qend sstart send" -out organelles.blast -num_threads "$CORES"
		fi
		awk -F"\t" '{if ($5 > 90) {print $1} }' organelles.blast | sort | uniq > organelle_contigs.txt
		awk -F"\t" '{if ($5 < 90 && $4 > 1500) {print $0} }' organelles.blast | sort -k1,1 -k7,7 | uniq > potential_misassemblies.txt
		/home/ps997/scripts/removeScaffoldsFromFasta.py 0_assembly.fa organelle_contigs.txt
		ln -sf 0_assembly.fa_organelle_contigs.txt-removed.fasta ../1a_RepeatModeler/1_assembly.fa
		ln -sf 0_assembly.fa_organelle_contigs.txt-removed.fasta ../1b_EDTA/1_assembly.fa
	fi

	cd "$HOMEDIR"/repeat-masking/"$i"/1a_RepeatModeler/
	if [[ -s 1_assembly.fa_db-families.fa.LTRlib.known.final.fa ]]
		then printf "RepeatModeler repeat library already exists for %s, skipping..." $i
	else
		if [ ! -e 1_assembly.fa ]
			then ln -s "$HOMEDIR"/polishing/"$i"/"$SAMPLE_NAME"_"$i"_polished_assembly.fasta 1_assembly.fa
		fi
		if [ ! -f 1_assembly.fa_db.nsq ]
			then /home/ps997/bin/RepeatModeler-2.0.1/BuildDatabase -name 1_assembly.fa_db 1_assembly.fa
		fi
		if [ ! -f 1_assembly.fa_db-families.fa ]
			then /home/ps997/bin/RepeatModeler-2.0.1/RepeatModeler -database 1_assembly.fa_db -pa $RM_CORES -LTRStruct > repeatmodeler.out 2> repeatmodeler.err
		fi
		/home/ps997/scripts/extract_unknownLTR.py 1_assembly.fa_db-families.fa
		if [ ! -s 1_assembly.fa_db-families.fa.LTRlib.unknown.fa.blastx ]
			then blastx -query 1_assembly.fa_db-families.fa.LTRlib.unknown.fa -db /home/fay-wei/bin/Tpases020812 -evalue 1e-10 -num_descriptions 10 -out 1_assembly.fa_db-families.fa.LTRlib.unknown.fa.blastx -num_threads "$CORES"
		fi
		perl /home/fay-wei/bin/Custom-Repeat-Library/transposon_blast_parse.pl --blastx *blastx --modelerunknown 1_assembly.fa_db-families.fa.LTRlib.unknown.fa
		cat 1_assembly.fa_db-families.fa.LTRlib.known.fa identified_elements.txt > 1_assembly.fa_db-families.fa.LTRlib.known.final.fa
		mv unknown_elements.txt 1_assembly.fa_db-families.fa.LTRlib.unknown.final.fa
		/home/ps997/scripts/rename_seq.py 1_assembly.fa_db-families.fa.LTRlib.known.final.fa
	fi
	cd "$HOMEDIR"/repeat-masking
done

for i in "${ASSEMBLIES_TO_MASK[@]}"
	do cd "$HOMEDIR"/repeat-masking/"$i"/1b_EDTA/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [[ -s 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.known.final.fa ]]
		then printf "EDTA repeat library exists for %s, skipping..." $i
	else
		conda activate EDTA_v2
		#export PATH="/home/ps997/bin/NINJA-0.95-cluster_only/NINJA:/home/ps997/bin/cd-hit-v4.8.1-2019-0228:/home/ps997/bin/bbmap:/home/ps997/bin/BBMap:/home/ps997/scripts:/home/fay-wei/bin/racon/build/bin:/home/fay-wei/bin/jellyfish-2.2.5/bin:/home/ps997/miniconda3/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/usr/lib/jvm/java-8-oracle/bin:/usr/lib/jvm/java-8-oracle/db/bin:/usr/lib/jvm/java-8-oracle/jre/bin:/opt/dell/srvadmin/bin:/home/ps997/edirect:/home/ps997/edirect"
		if [ ! -e 1_assembly.fa ]
				then ln -sf "$HOMEDIR"/polishing/"$i"/"$SAMPLE_NAME"_"$i"_polished_assembly.fasta 1_assembly.fa
		fi
		if [ ! -s 1_assembly.fa.mod.EDTA.TElib.fa ]
			then ~/bin/EDTA/EDTA.pl --sensitive 1 --anno 1 --evaluate 1 -t $CORES --genome 1_assembly.fa --repeatmasker /home/ps997/bin/RepeatMasker/RepeatMasker 2> edta.out
		fi
		/home/ps997/scripts/extract_unknownLTR.py 1_assembly.fa.mod.EDTA.TElib.fa
		if [ ! -s 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.fa.blastx ]
			then blastx -query 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.fa -db /home/fay-wei/bin/Tpases020812 -evalue 1e-10 -num_descriptions 10 -out 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.fa.blastx -num_threads "$CORES"
		fi
		perl /home/fay-wei/bin/Custom-Repeat-Library/transposon_blast_parse.pl --blastx *blastx --modelerunknown 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.fa
		cat 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.known.fa identified_elements.txt > 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.known.final.fa
		mv unknown_elements.txt 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.unknown.final.fa
		/home/ps997/scripts/rename_seq.py 1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.known.final.fa
	fi
	cd $HOMEDIR/repeat-masking
done

conda activate base

for i in "${ASSEMBLIES_TO_MASK[@]}"
	do cd "$HOMEDIR"/repeat-masking/"$i"/2_ProteinExclude/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -e 2_RepeatModeler.LTRlib.known.fa ]
		then ln -sf ../1a_RepeatModeler/1_assembly.fa_db-families.fa.LTRlib.known.final.fa.renamed.fa 2_RepeatModeler.LTRlib.known.fa
	fi
	if [ ! -e 2_EDTA.LTRLib.known.fa ]
		then ln -sf ../1b_EDTA/1_assembly.fa.mod.EDTA.TElib.fa.LTRlib.known.final.fa.renamed.fa 2_EDTA.LTRLib.known.fa
	fi
	if [[ ! -s 2_RepeatModeler.LTRlib.known.fanoProtFinal ]]
		then blastx -db /home/fay-wei/bin/uniprot/uniprot_sprot_plants.fasta -query 2_RepeatModeler.LTRlib.known.fa -out RepeatModeler2uniprot_plant_blast.out -num_threads $CORES
		perl /home/ps997/scripts/ProtExcluder1.1/ProtExcluder.pl RepeatModeler2uniprot_plant_blast.out 2_RepeatModeler.LTRlib.known.fa
	fi
	if [[ ! -s 2_EDTA.LTRLib.known.fanoProtFinal ]]
		then blastx -db /home/fay-wei/bin/uniprot/uniprot_sprot_plants.fasta -query 2_EDTA.LTRLib.known.fa -out EDTA2uniprot_plant_blast.out -num_threads $CORES
		perl /home/ps997/scripts/ProtExcluder1.1/ProtExcluder.pl EDTA2uniprot_plant_blast.out 2_EDTA.LTRLib.known.fa
	fi

	cd "$HOMEDIR"/repeat-masking
done

for i in "${ASSEMBLIES_TO_MASK[@]}"
	do cd "$HOMEDIR"/repeat-masking/"$i"/3_RepeatMasker
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -d EDTA_out ]
		then mkdir EDTA_out
	fi
	if [ ! -d RepeatModeler_out ]
		then mkdir RepeatModeler_out
	fi
	if [ ! -e 3_RepeatModeler_LTRLib ]
		then ln -sf ../2_ProteinExclude/2_RepeatModeler.LTRlib.known.fanoProtFinal 3_RepeatModeler_LTRLib
	fi
	if [ ! -e 3_EDTA_LTRLib ]
		then ln -sf ../2_ProteinExclude/2_EDTA.LTRLib.known.fanoProtFinal 3_EDTA_LTRLib
	fi
	if [ ! -e 3_assembly.fa ]
		then ln -sf ../1a_RepeatModeler/1_assembly.fa 3_assembly.fa
	fi

	if [ ! -s RepeatModeler_out/3_assembly.fa.masked ]
		then /home/ps997/bin/RepeatMasker/RepeatMasker -noisy -a -gff -u -pa $CORES -lib 3_RepeatModeler_LTRLib 3_assembly.fa
		rsync --remove-source-files 3_assembly.fa.* RepeatModeler_out
	fi
	if [ ! -s EDTA_out/3_assembly.fa.masked ]
		then /home/ps997/bin/RepeatMasker/RepeatMasker -noisy -a -gff -u -pa $CORES -lib 3_EDTA_LTRLib 3_assembly.fa
		rsync --remove-source-files 3_assembly.fa.* EDTA_out
	fi

	if [ ! -s EDTA_detailed_summary ]
		then /home/ps997/bin/RepeatMasker/util/buildSummary.pl -species "$SAMPLE_NAME" -useAbsoluteGenomeSize EDTA_out/3_assembly.fa.out > EDTA_detailed_summary
	fi
	if [ ! -s RepeatModeler_detailed_summary ]
		then /home/ps997/bin/RepeatMasker/util/buildSummary.pl -species "$SAMPLE_NAME" -useAbsoluteGenomeSize RepeatModeler_out/3_assembly.fa.out > RepeatModeler_detailed_summary
	fi
	cd "$HOMEDIR"/repeat-masking
done
cd $HOMEDIR

exit 0

if [ ! -d gene-prediction ]
	then mkdir gene-prediction
fi
cd gene-prediction

for i in $ASSEMBLER
	do if [ ! -d $i ]
		then mkdir $i
	fi
done

for i in $ASSEMBLER
	do cd "$i"
	if [ ! -d 1a_RepeatModeler ]
		then mkdir 1a_RepeatModeler
	fi
	if [ ! -d 1b_EDTA ]
		then mkdir 1b_EDTA
	fi
	if [ ! -d 1c_EDTA_nofilter ]
		then mkdir 1c_EDTA_nofilter
	fi
	for j in 1a_RepeatModeler 1b_EDTA 1c_EDTA_nofilter
		do cd "$j"
		if [ ! -e RNA_R1.fq.gz ]
			then ln -sf $RNA1 RNA_R1.fq.gz
		fi
		if [ ! -e RNA_R2.fq.gz ]
			then ln -sf $RNA2 RNA_R2.fq.gz
		fi
		if [ ! -e Hornwort_orthogroups.faa ]
			then ln -sf $REFPROT Hornwort_orthogroups.faa
		fi
		cd ..
	done

	cd 1a_RepeatModeler/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -e 1_RepeatModeler_masked_assembly.fa ]
		then ln -sf ../../../repeat-masking/"$i"/3_RepeatMasker/RepeatModeler_out/3_assembly.fa.masked 1_RepeatModeler_masked_assembly.fa
	fi
	if [ ! -e RNA_mapped.bam ] || [ $(stat --printf="%s" RNA_mapped.bam) -lt 100 ]
		then hisat2-build -p "$CORES" 1_RepeatModeler_masked_assembly.fa 1_RepeatModeler_masked_assembly >& hisat-build.out
		hisat2 -p "$CORES" -x 1_RepeatModeler_masked_assembly -1 RNA_R1.fq.gz -2 RNA_R2.fq.gz 2> hisat-align.out | samtools view -b | samtools sort > RNA_mapped.bam
		samtools flagstat RNA_mapped.bam > RNA_mapped.bam.flagstats
		samtools stats RNA_mapped.bam > RNA_mapped.bam.samstats
	fi
	if [ ! -e RNA_mapped.bam.alignmentsummary ]
		then picard CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=1_RepeatModeler_masked_assembly.fa INPUT=RNA_mapped.bam OUTPUT=RNA_mapped.bam.alignmentsummary
	fi

	cd ../1b_EDTA/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -e 1_EDTA_masked_assembly.fa ]
		then ln -sf ../../../repeat-masking/"$i"/3_RepeatMasker/EDTA_out/3_assembly.fa.masked 1_EDTA_masked_assembly.fa
	fi
	if [ ! -e RNA_mapped.bam ] || [ $(stat --printf="%s" RNA_mapped.bam) -lt 100 ]
		then hisat2-build -p $CORES 1_EDTA_masked_assembly.fa 1_EDTA_masked_assembly >& hisat-build.out
		hisat2 -p "$CORES" -x 1_EDTA_masked_assembly -1 RNA_R1.fq.gz -2 RNA_R2.fq.gz 2> hisat-align.out | samtools view -b | samtools sort > RNA_mapped.bam
		samtools flagstat RNA_mapped.bam > RNA_mapped.bam.flagstats
		samtools stats RNA_mapped.bam > RNA_mapped.bam.samstats
	fi
	if [ ! -e RNA_mapped.bam.alignmentsummary ]
		then picard CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=1_EDTA_masked_assembly.fa INPUT=RNA_mapped.bam OUTPUT=RNA_mapped.bam.alignmentsummary
	fi

	cd ../1c_EDTA_nofilter
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH

	if [ ! -e 1_assembly.fa.mod.MAKER.masked.fa ]
		then ln -sf ../../../repeat-masking/"$i"/1b_EDTA/1_assembly.fa.mod.MAKER.masked 1_assembly.fa.mod.MAKER.masked.fa
	fi
	if [ ! -e RNA_mapped.bam ] || [ $(stat --printf="%s" RNA_mapped.bam) -lt 100 ]
	then hisat2-build -p $CORES 1_assembly.fa.mod.MAKER.masked.fa 1_assembly.fa.mod.MAKER.masked >& hisat-build.out
		hisat2 -p "$CORES" -x 1_assembly.fa.mod.MAKER.masked -1 RNA_R1.fq.gz -2 RNA_R2.fq.gz 2> hisat-align.out | samtools view -b | samtools sort > RNA_mapped.bam
		samtools flagstat RNA_mapped.bam > RNA_mapped.bam.flagstats
		samtools stats RNA_mapped.bam > RNA_mapped.bam.samstats
	fi
	if [ ! -e RNA_mapped.bam.alignmentsummary ]
		then picard CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=1_assembly.fa.mod.MAKER.masked.fa INPUT=RNA_mapped.bam OUTPUT=RNA_mapped.bam.alignmentsummary
	fi

	cd ../1a_RepeatModeler/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	conda activate braker
	if [ ! -f braker-RNA/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_RepeatModeler_masked_assembly.fa --species Sp_"$k" --bam RNA_mapped.bam --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA/
	fi
	if [ ! -f braker-RNA+PROT/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_RepeatModeler_masked_assembly.fa --species Sp_"$k" --bam RNA_mapped.bam --prot_seq Hornwort_orthogroups.faa --prg=gth --gth2traingenes --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.log 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA+PROT/
	fi

	cd ../1b_EDTA/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -f braker-RNA/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_EDTA_masked_assembly.fa --species Sp_"$k" --bam RNA_mapped.bam --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA/
	fi
	if [ ! -f braker-RNA+PROT/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_EDTA_masked_assembly.fa --species Sp_"$k" --bam RNA_mapped.bam --prot_seq Hornwort_orthogroups.faa --prg=gth --gth2traingenes --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA+PROT/
	fi

	cd ../1c_EDTA_nofilter/
	CURRENT_PATH=$(pwd)
	echo $CURRENT_PATH
	if [ ! -f braker-RNA/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_assembly.fa.mod.MAKER.masked.fa --species Sp_"$k" --bam RNA_mapped.bam --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA/
	fi
	if [ ! -f braker-RNA+PROT/augustus.hints.aa ]
		then k=1
		for k in {1..10000} ; do if [ -d /home/ps997/bin/Augustus/config/species/Sp_"$k" ] ; then true ; else break ; fi ; done
		echo "Running Braker as Sp_$k"
		/home/ps997/miniconda3/envs/braker/bin/braker.pl --genome 1_assembly.fa.mod.MAKER.masked.fa --species Sp_"$k" --bam RNA_mapped.bam --prot_seq Hornwort_orthogroups.faa --prg=gth --gth2traingenes --verbosity 3 --cores $CORES --GENEMARK_PATH=/home/ps997/bin/gmes_linux_64 --AUGUSTUS_CONFIG_PATH=/home/ps997/bin/Augustus/config --AUGUSTUS_BIN_PATH=/home/ps997/miniconda3/envs/braker/bin/ --AUGUSTUS_SCRIPTS_PATH=/home/ps997/bin/Augustus/scripts 1>braker.out 2>braker.err
		rsync -ar --remove-source-files braker/* braker-RNA+PROT/
	fi
	conda activate base
	cd $HOMEDIR/gene-prediction
done
conda activate base
