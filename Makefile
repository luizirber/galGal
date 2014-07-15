SHELL=/bin/bash

all: moleculo_galGal4

moleculo_galGal4: outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.sorted.bam

inputs/galGal4/galGal4.fa.masked.gz:
	wget -SNc ftp://hgdownload.cse.ucsc.edu/goldenPath/galGal4/bigZips/galGal4.fa.masked.gz -P inputs/galGal4/

outputs/galGal4/galGal4.fa.masked: inputs/galGal4/galGal4.fa.masked.gz
	mkdir -p outputs/galGal4
	cp -a {inputs,outputs}/galGal4/galGal4.fa.masked.gz
	gunzip -f outputs/galGal4/galGal4.fa.masked.gz






outputs/galGal4/galGal4.fa.masked.sa: outputs/galGal4/galGal4.fa.masked
	module try-load bwa
	bwa index outputs/galGal4/galGal4.fa.masked

outputs/galGal4/galGal4.fa.masked.fai: outputs/galGal4/galGal4.fa.masked
	module try-load samtools
	samtools faidx outputs/galGal4/galGal4.fa.masked



outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq: inputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.gz
	mkdir -p outputs/moleculo
	cp -a {inputs,outputs}/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.gz
	gunzip -f outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.gz



outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.sorted.bam: outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq outputs/galGal4/galGal4.fa.masked.sa
	module try-load bwa
	module try-load samtools
	bwa mem outputs/galGal4/galGal4.fa.masked outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq > outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.sam
	samtools import outputs/galGal4/galGal4.fa.masked.fai outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.sam outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.bam
	samtools sort outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.bam outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.sorted
	samtools index outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.sorted.bam
