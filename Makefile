SHELL=/bin/bash

.PHONY: clean all

all: moleculo_galGal4 moleculo_galGal4_masked

moleculo_galGal4_masked: outputs/moleculo/galGal4.masked.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam

moleculo_galGal4: outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam

inputs/galGal4/galGal4.%.gz:
	wget -SNc ftp://hgdownload.cse.ucsc.edu/goldenPath/galGal4/bigZips/$(@F) -P inputs/galGal4/

inputs/uniprot/uniprot_sprot.fasta.gz:
	wget -SNc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/$(@F) -P inputs/uniprot/

#######################################################################

outputs/uniprot/uniprot_sprot.fasta: inputs/uniprot/uniprot_sprot.fasta.gz
	mkdir -p outputs/uniprot
	cp -a {inputs,outputs}/uniprot/$(<F)
	gunzip -f outputs/uniprot/$(<F)

outputs/moleculo/%-AllReads.sorted.bam: \
  outputs/moleculo/%-1_LongRead_500_1499nt.fastq.sorted.bam \
  outputs/moleculo/%-1_LongRead.fastq.sorted.bam \
  outputs/moleculo/%-2_LongRead_500_1499nt.fastq.sorted.bam \
  outputs/moleculo/%-2_LongRead.fastq.sorted.bam \
  outputs/moleculo/%-3_LongRead_500_1499nt.fastq.sorted.bam \
  outputs/moleculo/%-3_LongRead.fastq.sorted.bam \
  outputs/moleculo/%-4_LongRead_500_1499nt.fastq.sorted.bam \
  outputs/moleculo/%-4_LongRead.fastq.sorted.bam \
  outputs/moleculo/%-5_LongRead_500_1499nt.fastq.sorted.bam \
  outputs/moleculo/%-5_LongRead.fastq.sorted.bam
	module try-load samtools
	samtools merge $@ $^

outputs/galGal4/galGal4.%.sa: outputs/galGal4/galGal4.%
	module try-load bwa
	bwa index $<

outputs/galGal4/galGal4.%.fai: outputs/galGal4/galGal4.%
	module try-load samtools
	samtools faidx $<

outputs/galGal4/galGal4.%: inputs/galGal4/galGal4.%.gz
	mkdir -p outputs/galGal4
	cp {inputs,outputs}/galGal4/$(<F)
	gunzip -fN outputs/galGal4/$(<F)

outputs/moleculo/%.fastq: inputs/moleculo/%.fastq.gz
	mkdir -p outputs/moleculo
	cp -a $< $@.gz
	gunzip -f $@.gz

outputs/moleculo/%.fastq.sorted.bam: outputs/moleculo/%.fastq.bam
	module try-load samtools
	samtools sort $< $(basename $@ .bam)
	samtools index $@

outputs/moleculo/galGal4.masked.%.fastq.bam: outputs/moleculo/%.fastq outputs/galGal4/galGal4.fa.masked.sa outputs/galGal4/galGal4.fa.masked.fai
	module try-load bwa
	module try-load samtools
	bwa mem outputs/galGal4/galGal4.fa.masked $< > $<.sam.galGal4.masked
	samtools import outputs/galGal4/galGal4.fa.masked.fai $<.sam.galGal4.masked $@

outputs/moleculo/galGal4.%.fastq.bam: outputs/moleculo/%.fastq outputs/galGal4/galGal4.fa.sa outputs/galGal4/galGal4.fa.fai
	module try-load bwa
	module try-load samtools
	bwa mem outputs/galGal4/galGal4.fa $< > $<.sam.galGal4
	samtools import outputs/galGal4/galGal4.fa.fai $<.sam.galGal4 $@

clean:
	- rm outputs/moleculo/*.sam.galGal4
	- rm outputs/moleculo/*.sam.galGal4.masked
	- rm outputs/moleculo/*.fastq.bam
	# -rm outputs/galGal4/*.{amb,pac,ann}

.PRECIOUS: %.sorted.bam
