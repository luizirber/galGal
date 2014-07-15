SHELL=/bin/bash

all: moleculo_galGal4

moleculo_galGal4: outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-2_LongRead_500_1499nt.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-2_LongRead.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-3_LongRead_500_1499nt.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-3_LongRead.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-4_LongRead_500_1499nt.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-4_LongRead.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-5_LongRead_500_1499nt.fastq.sorted.bam \
	outputs/moleculo/LR6000017-DNA_A01-LRAAA-5_LongRead.fastq.sorted.bam

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



outputs/moleculo/%.fastq: inputs/moleculo/%.fastq.gz
	mkdir -p outputs/moleculo
	cp -a $< $@.gz
	gunzip -f $@.gz

%.fastq.sorted.bam: %.fastq outputs/galGal4/galGal4.fa.masked.sa outputs/galGal4/galGal4.fa.masked.fai
	module try-load bwa
	module try-load samtools
	bwa mem outputs/galGal4/galGal4.fa.masked $< > $<.sam
	samtools import outputs/galGal4/galGal4.fa.masked.fai $<.sam $<.bam
	samtools sort $<.bam $<.sorted
	samtools index $@
