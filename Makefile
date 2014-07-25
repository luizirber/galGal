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
# PBS stuff
#######################################################################

workdirs/%.pbs:
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${COVERAGE_RES} -N cov.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 60; done
	@grep "galGal PBS job finished: SUCCESS" $@

#######################################################################
# Inputs
#######################################################################

inputs/galGal4/galGal4.%.gz:
	wget -SNc ftp://hgdownload.cse.ucsc.edu/goldenPath/galGal4/bigZips/$(@F) -P inputs/galGal4/

inputs/galGal3/galGal3.fa:
	wget -SNc http://hgdownload.soe.ucsc.edu/goldenPath/galGal3/bigZips/galGal3.2bit -P inputs/galGal3/

inputs/uniprot/uniprot_sprot.fasta.gz:
	wget -SNc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/$(@F) -P inputs/uniprot/

inputs/msu/Chicken.tar.gz:
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

outputs/galGal5/galGal5.fa: $(GALGAL5_SOURCES)
	mkdir -p outputs/galGal5
	scripts/galGal5_fix.sh "$^" > $@

#outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq: outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq_screed
#	python -m screed.dump_to_fastq $< $@

outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq: $(MOLECULO_READS)
	cat $^ > $@

outputs/moleculo/%.fastq: inputs/moleculo/%.fastq.gz
	mkdir -p outputs/moleculo
	cp -a $< $@.gz
	gunzip -f $@.gz

outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq_screed: $(addsuffix _screed, $(MOLECULO_READS))
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

#######################################################################

outputs/coverage/%.pd_df.csv: outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq \
  outputs/galGal4/galGal4.fa outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam \
  outputs/galGal5/galGal5.fa outputs/moleculo/galGal5.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam \
  workdirs/%/output/output.00500.pbs workdirs/%/output/output.01000.pbs \
  workdirs/%/output/output.01500.pbs workdirs/%/output/output.02000.pbs \
  workdirs/%/output/output.02500.pbs workdirs/%/output/output.03000.pbs \
  workdirs/%/output/output.03500.pbs workdirs/%/output/output.03000.pbs \
  workdirs/%/output/output.03500.pbs workdirs/%/output/output.04000.pbs \
  workdirs/%/output/output.04500.pbs workdirs/%/output/output.05000.pbs \
  workdirs/%/output/output.05500.pbs workdirs/%/output/output.06000.pbs \
  workdirs/%/output/output.06500.pbs workdirs/%/output/output.07000.pbs \
  workdirs/%/output/output.07500.pbs workdirs/%/output/output.08000.pbs \
  workdirs/%/output/output.08500.pbs workdirs/%/output/output.09000.pbs \
  workdirs/%/output/output.09500.pbs workdirs/%/output/output.10000.pbs \
  workdirs/%/output/output.10500.pbs workdirs/%/output/output.11000.pbs \
  workdirs/%/output/output.11500.pbs workdirs/%/output/output.12000.pbs \
  workdirs/%/output/output.12500.pbs workdirs/%/output/output.13000.pbs \
  workdirs/%/output/output.13500.pbs workdirs/%/output/output.14000.pbs \
  workdirs/%/output/output.14500.pbs workdirs/%/output/output.15000.pbs \
  workdirs/%_90/output/output.00500.pbs workdirs/%_90/output/output.01000.pbs \
  workdirs/%_90/output/output.01500.pbs workdirs/%_90/output/output.02000.pbs \
  workdirs/%_90/output/output.02500.pbs workdirs/%_90/output/output.03000.pbs \
  workdirs/%_90/output/output.03500.pbs workdirs/%_90/output/output.03000.pbs \
  workdirs/%_90/output/output.03500.pbs workdirs/%_90/output/output.04000.pbs \
  workdirs/%_90/output/output.04500.pbs workdirs/%_90/output/output.05000.pbs \
  workdirs/%_90/output/output.05500.pbs workdirs/%_90/output/output.06000.pbs \
  workdirs/%_90/output/output.06500.pbs workdirs/%_90/output/output.07000.pbs \
  workdirs/%_90/output/output.07500.pbs workdirs/%_90/output/output.08000.pbs \
  workdirs/%_90/output/output.08500.pbs workdirs/%_90/output/output.09000.pbs \
  workdirs/%_90/output/output.09500.pbs workdirs/%_90/output/output.10000.pbs \
  workdirs/%_90/output/output.10500.pbs workdirs/%_90/output/output.11000.pbs \
  workdirs/%_90/output/output.11500.pbs workdirs/%_90/output/output.12000.pbs \
  workdirs/%_90/output/output.12500.pbs workdirs/%_90/output/output.13000.pbs \
  workdirs/%_90/output/output.13500.pbs workdirs/%_90/output/output.14000.pbs \
  workdirs/%_90/output/output.14500.pbs workdirs/%_90/output/output.15000.pbs
	mkdir -p $(@D)
	python scripts/count_reads_pd.py $< workdirs/$(*F)/output workdirs/$(*F)_90/output $@

workdirs/galGal4_90/output/output.%: outputs/galGal4/galGal4.fa outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 --minlen=0.9 1>$@

workdirs/galGal5_90/output/output.%: outputs/galGal5/galGal5.fa outputs/moleculo/galGal5.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 --minlen=0.9 1>$@

workdirs/galGal4/output/output.%: outputs/galGal4/galGal4.fa outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 1>$@

workdirs/galGal5/output/output.%: outputs/galGal5/galGal5.fa outputs/moleculo/galGal5.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 1>$@

#######################################################################

outputs/uniprot/uniprot.namedb: outputs/uniprot/uniprot_sprot.fasta
	cd outputs/uniprot && formatdb -i uniprot_sprot.fasta -o T -p T
	cd outputs/uniprot && python ../../scripts/make-namedb.py uniprot_sprot.fasta uniprot.namedb
	cd outputs/uniprot && python -m screed.fadbm uniprot_sprot.fasta

#######################################################################

clean:
	- rm outputs/moleculo/*.sam.galGal*
	- rm outputs/moleculo/*.fastq.bam
	# -rm outputs/galGal4/*.{amb,pac,ann}

dependencies.png:
	make -Bnd | make2graph > output.dot
	dot -Tpng output.dot > dependencies.png
	-rm output.dot

#.PRECIOUS: %.sorted.bam
.SECONDARY:
