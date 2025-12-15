#!/bin/bash

while getopts a:r:t:i:x:y: option
	do
		case "${option}"
		in
		a) assembly=${OPTARG};;
		r) reads=${OPTARG};;
		t) threads=${OPTARG};;
		i) iteration_num=${OPTARG};;
		x) illumina1=${OPTARG};;
		y) illumina2=${OPTARG};;
		esac
	done

echo ${assembly} ${reads} ${threads} ${iteration_num} ${illumina1} ${illumina2}

# assembly=
# reads=
# threads=4
# iteration_num=4

iteration=1
while [ $iteration -le $iteration_num ]

	do
	echo Iteration: ${iteration}
	minimap2 -ax lr:hq -t $threads $assembly $reads | samtools sort -o ont.bam
	minimap2 -ax sr -t $threads $assembly $illumina1 $illumina2 | samtools sort -o illumina.bam
	#bwa index ${assembly}
	#bwa mem -t $threads $assembly $illumina1 $illumina2 | samtools sort -o bwa.bam
	samtools index illumina.bam
	samtools index ont.bam
	java -Xmx100G -jar /home/ps997/bin/pilon-1.24.jar --genome $assembly --frags illumina.bam --nanopore ont.bam --output pilon.${iteration}

	assembly=pilon.${iteration}.fasta
	iteration=$((iteration+1))

	done
