include Makefile.inc

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

inputs/msu/Chicker.tar.gz:
	wget -SNc https://dl.dropboxusercontent.com/u/1455804/Chicken.tar.gz -P inputs/msu/

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

#outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq: outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq_screed
#	python -m screed.dump_to_fastq $< $@

outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq: \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead.fastq \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-2_LongRead_500_1499nt.fastq \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-2_LongRead.fastq \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-3_LongRead_500_1499nt.fastq \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-3_LongRead.fastq \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-4_LongRead_500_1499nt.fastq \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-4_LongRead.fastq \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-5_LongRead_500_1499nt.fastq \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-5_LongRead.fastq
	cat $^ > $@

outputs/moleculo/%.fastq: inputs/moleculo/%.fastq.gz
	mkdir -p outputs/moleculo
	cp -a $< $@.gz
	gunzip -f $@.gz

outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq_screed: \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead_500_1499nt.fastq_screed \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-1_LongRead.fastq_screed \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-2_LongRead_500_1499nt.fastq_screed \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-2_LongRead.fastq_screed \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-3_LongRead_500_1499nt.fastq_screed \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-3_LongRead.fastq_screed \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-4_LongRead_500_1499nt.fastq_screed \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-4_LongRead.fastq_screed \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-5_LongRead_500_1499nt.fastq_screed \
  outputs/moleculo/LR6000017-DNA_A01-LRAAA-5_LongRead.fastq_screed
	python scripts/merge_repeated.py $@ $^

outputs/moleculo/%.fastq_screed: outputs/moleculo/%.fastq
	python -m screed.fqdbm $<

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

outputs/coverage/%.pd_df.csv: outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq \
  outputs/galGal4/galGal4.fa outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam \
  workdirs/%/output/output.00500 workdirs/%/output/output.01000 \
  workdirs/%/output/output.01500 workdirs/%/output/output.02000 \
  workdirs/%/output/output.02500 workdirs/%/output/output.03000 \
  workdirs/%/output/output.03500 workdirs/%/output/output.03000 \
  workdirs/%/output/output.03500 workdirs/%/output/output.04000 \
  workdirs/%/output/output.04500 workdirs/%/output/output.05000 \
  workdirs/%/output/output.05500 workdirs/%/output/output.06000 \
  workdirs/%/output/output.06500 workdirs/%/output/output.07000 \
  workdirs/%/output/output.07500 workdirs/%/output/output.08000 \
  workdirs/%/output/output.08500 workdirs/%/output/output.09000 \
  workdirs/%/output/output.09500 workdirs/%/output/output.10000 \
  workdirs/%/output/output.10500 workdirs/%/output/output.11000 \
  workdirs/%/output/output.11500 workdirs/%/output/output.12000 \
  workdirs/%/output/output.12500 workdirs/%/output/output.13000 \
  workdirs/%/output/output.13500 workdirs/%/output/output.14000 \
  workdirs/%/output/output.14500 workdirs/%/output/output.15000 \
  workdirs/%_90/output/output.00500 workdirs/%_90/output/output.01000 \
  workdirs/%_90/output/output.01500 workdirs/%_90/output/output.02000 \
  workdirs/%_90/output/output.02500 workdirs/%_90/output/output.03000 \
  workdirs/%_90/output/output.03500 workdirs/%_90/output/output.03000 \
  workdirs/%_90/output/output.03500 workdirs/%_90/output/output.04000 \
  workdirs/%_90/output/output.04500 workdirs/%_90/output/output.05000 \
  workdirs/%_90/output/output.05500 workdirs/%_90/output/output.06000 \
  workdirs/%_90/output/output.06500 workdirs/%_90/output/output.07000 \
  workdirs/%_90/output/output.07500 workdirs/%_90/output/output.08000 \
  workdirs/%_90/output/output.08500 workdirs/%_90/output/output.09000 \
  workdirs/%_90/output/output.09500 workdirs/%_90/output/output.10000 \
  workdirs/%_90/output/output.10500 workdirs/%_90/output/output.11000 \
  workdirs/%_90/output/output.11500 workdirs/%_90/output/output.12000 \
  workdirs/%_90/output/output.12500 workdirs/%_90/output/output.13000 \
  workdirs/%_90/output/output.13500 workdirs/%_90/output/output.14000 \
  workdirs/%_90/output/output.14500 workdirs/%_90/output/output.15000
	echo python scripts/count_reads_pd.py $< workdirs/$(*F)/output workdirs/$(*F)_90/output $@

workdirs/%.pbs: $(@:.pbs=)
	echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | qsub -l ${COVERAGE_RES} -N cov.${subst output.,,$(@F)} -o $@.out -e $@.err

#workdirs/%_90/output/output.00500: $(subst REF,galGal4,outputs/REF/REF.fa) outputs/moleculo/%.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
#workdirs/galGal4_90/output/output.00500: outputs/galGal4/galGal4.fa outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
workdirs/galGal4_90/output/output.% : outputs/galGal4/galGal4.fa outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 --minlen=0.9 1>$@

#workdirs/galGal4/output/output.00500: inputs/galGal4/galGal4.fa.gz outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq
workdirs/galGal4/output/output.%: outputs/galGal4/galGal4.fa outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 1>$@

outputs/uniprot/uniprot.namedb: outputs/uniprot/uniprot_sprot.fasta
	cd outputs/uniprot && formatdb -i uniprot_sprot.fasta -o T -p T
	cd outputs/uniprot && python ../../scripts/make-namedb.py uniprot_sprot.fasta uniprot.namedb
	cd outputs/uniprot && python -m screed.fadbm uniprot_sprot.fasta

#######################################################################

clean:
	- rm outputs/moleculo/*.sam.galGal*
	- rm outputs/moleculo/*.fastq.bam
	# -rm outputs/galGal4/*.{amb,pac,ann}

#.PRECIOUS: %.sorted.bam
.SECONDARY:
