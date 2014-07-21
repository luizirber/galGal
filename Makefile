SHELL=/bin/bash

.PHONY: clean all

all: moleculo_galGal4 moleculo_galGal3 moleculo_galGal5

#######################################################################
# High-level targets
#######################################################################

moleculo_galGal3: outputs/moleculo/galGal3.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam

moleculo_galGal4: outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam

moleculo_galGal5: outputs/moleculo/galGal5.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam

#######################################################################
# Inputs
#######################################################################

inputs/galGal4/galGal4.%.gz:
	wget -SNc ftp://hgdownload.cse.ucsc.edu/goldenPath/galGal4/bigZips/$(@F) -P inputs/galGal4/

inputs/galGal3/galGal3.fa:
	wget -SNc http://hgdownload.soe.ucsc.edu/goldenPath/galGal3/bigZips/galGal3.2bit -P inputs/galGal3/

inputs/uniprot/uniprot_sprot.fasta.gz:
	wget -SNc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/$(@F) -P inputs/uniprot/

#######################################################################
# Outputs
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
	samtools merge $@ $^
	samtools index $@

outputs/galGal%.sa: outputs/galGal%
	bwa index $<

outputs/galGal%.fai: outputs/galGal%
	samtools faidx $<

outputs/galGal4/galGal4.%: inputs/galGal4/galGal4.%.gz
	mkdir -p outputs/galGal4
	cp {inputs,outputs}/galGal4/$(<F)
	gunzip -fN outputs/galGal4/$(<F)

outputs/galGal3/galGal3.fa: inputs/galGal3/galGal3.2bit
	mkdir -p outputs/galGal3
	wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa -P outputs/galGal3/
	#wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa -P outputs/galGal3/
	chmod +x outputs/galGal3/twoBitToFa
	outputs/galGal3/twoBitToFa $< $@
	-rm outputs/galGal3/twoBitToFa

outputs/galGal5/galGal5.fa: inputs/galGal5/Chicken_Chr1_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr2_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr3_PBJelly_14.1.15.QUIVER.fasta \
inputs/galGal5/Chicken_Chr4_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr5_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr6_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr7_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr8_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr9_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr10_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr11_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr12_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr13_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr14_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr15_PBjelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr16_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr17_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr18_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr19_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr20_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr21_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr22_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr23_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr24_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr25_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr26_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr27_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr28_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_Chr32_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_ChrW_PBJelly_14.1.15_QUIVER.fasta \
inputs/galGal5/Chicken_ChrZ_PBJelly_14.1.15_QUIVER.fasta
	mkdir -p outputs/galGal5
	scripts/galGal5_fix.sh "$^" > $@

outputs/moleculo/%.fastq: inputs/moleculo/%.fastq.gz
	mkdir -p outputs/moleculo
	cp -a $< $@.gz
	gunzip -f $@.gz

outputs/moleculo/%.fastq.sorted.bam: outputs/moleculo/%.fastq.bam
	samtools sort $< $(basename $@ .bam)
	samtools index $@

outputs/moleculo/galGal4.%.fastq.bam: outputs/moleculo/%.fastq outputs/galGal4/galGal4.fa.sa outputs/galGal4/galGal4.fa.fai
	bwa mem outputs/galGal4/galGal4.fa $< > $<.sam.galGal4
	samtools import outputs/galGal4/galGal4.fa.fai $<.sam.galGal4 $@

outputs/moleculo/galGal3.%.fastq.bam: outputs/moleculo/%.fastq outputs/galGal3/galGal3.fa.sa outputs/galGal3/galGal3.fa.fai
	bwa mem outputs/galGal3/galGal3.fa $< > $<.sam.galGal3
	samtools import outputs/galGal3/galGal3.fa.fai $<.sam.galGal3 $@

outputs/moleculo/galGal5.%.fastq.bam: outputs/moleculo/%.fastq outputs/galGal5/galGal5.fa.sa outputs/galGal5/galGal5.fa.fai
	bwa mem outputs/galGal5/galGal5.fa $< > $<.sam.galGal5
	samtools import outputs/galGal5/galGal5.fa.fai $<.sam.galGal5 $@

#######################################################################

outputs/uniprot/uniprot.namedb: outputs/uniprot/uniprot_sprot.fasta
	cd outputs/uniprot && formatdb -i uniprot_sprot.fasta -o T -p T
	cd outputs/uniprot && python ../../scripts/make-namedb.py uniprot_sprot.fasta uniprot.namedb
	cd outputs/uniprot && python -m screed.fadbm uniprot_sprot.fasta

clean:
	- rm outputs/moleculo/*.sam.galGal*
	- rm outputs/moleculo/*.fastq.bam
	# -rm outputs/galGal4/*.{amb,pac,ann}

#.PRECIOUS: %.sorted.bam
.SECONDARY:
