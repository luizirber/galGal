include Makefile.inc

.PHONY: clean all

all: moleculo_galGal4 moleculo_galGal3 moleculo_galGal5

#######################################################################
# High-level targets
#######################################################################

moleculo_galGal3: outputs/moleculo/galGal3.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam

moleculo_galGal4: outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam

moleculo_galGal5: outputs/moleculo/galGal5.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam

msu_galGal3: outputs/msu/galGal3.msu.sorted.bam

msu_galGal4: outputs/msu/galGal4.msu.sorted.bam

msu_galGal5: outputs/msu/galGal5.msu.sorted.bam

msu_latest_galGal4: outputs/msu/galGal4.latest.sorted.bam

pacbio_galGal4: outputs/pacbio/galGal4.pacbio.sorted.bam

#######################################################################
# PBS stuff
#######################################################################

workdirs/galGal%.pbs:
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${COVERAGE_RES} -N cov.${subst output.,,$(@F)} -o $@ -e $@.err -A ged-intel11 | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 60; done
	@grep "galGal PBS job finished: SUCCESS" $@

workdirs/msu%.pbs:
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${COVERAGE_MSU_RES} -N cov.${subst output.,,$(@F)} -o $@ -e $@.err -A ged-intel11 | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 60; done
	@grep "galGal PBS job finished: SUCCESS" $@

workdirs/blat/transc_pacbio.pbs:
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${BLAT_PACBIO_RES} -N blat.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 600; done
	@grep "galGal PBS job finished: SUCCESS" $@

workdirs/blat/transc_%.pbs:
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${BLAT_RES} -N blat.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 600; done
	@grep "galGal PBS job finished: SUCCESS" $@

outputs/pacbio_assembly/chicken_2.fasta.pbs:
	mkdir -p $(@D)
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${PBCR_RES} -N pbcr.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 600; done
	@grep "galGal PBS job finished: SUCCESS" $@

outputs/pacbio_unmapped/chicken_unmapped.fasta.pbs:
	mkdir -p $(@D)
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${PBCR_RES} -N pbcr.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 600; done
	@grep "galGal PBS job finished: SUCCESS" $@

outputs/pacbio/galGal4.sam.pbs:
	mkdir -p $(@D)
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${BLASR_RES} -N blasr.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 600; done
	@grep "galGal PBS job finished: SUCCESS" $@

workdir/results/Moleculo_RNA_BLAST.txt.pbs:
	mkdir -p $(@D)
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${BLAST_MOL_RES} -N blast.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 600; done
	@grep "galGal PBS job finished: SUCCESS" $@

workdir/results/Chick_RNA_BLAST.txt.pbs:
	mkdir -p $(@D)
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${BLAST_REF_RES} -N blast.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 600; done
	@grep "galGal PBS job finished: SUCCESS" $@

workdir/results/%.pbs:
	mkdir -p $(@D)
	JOBID=`echo make $(subst .pbs,,$@) | cat pbs/header.sub - pbs/footer.sub | \
	  qsub -l ${BLAST_RES} -N blast.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 600; done
	@grep "galGal PBS job finished: SUCCESS" $@

#######################################################################
# Inputs
#######################################################################

inputs/reference/galGal4.%.gz:
	mkdir -p $(@D)
	wget -SNc ftp://hgdownload.cse.ucsc.edu/goldenPath/galGal4/bigZips/$(@F) -P inputs/reference/

inputs/reference/galGal3.2bit:
	mkdir -p $(@D)
	wget -SNc http://hgdownload.soe.ucsc.edu/goldenPath/galGal3/bigZips/galGal3.2bit -P inputs/reference/

inputs/uniprot/uniprot_sprot.fasta.gz:
	wget -SNc ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/$(@F) -P inputs/uniprot/

#inputs/msu/Chicken.tar.gz:
#	wget -SNc https://dl.dropboxusercontent.com/u/1455804/Chicken.tar.gz -P inputs/msu/

#inputs/msu/chicken_contig_list.all.txt.gz:
#	wget -SNc https://dl.dropboxusercontent.com/u/1455804/chicken_contig_list.all.txt.gz -P inputs/msu/

#######################################################################
# Outputs
#######################################################################

outputs/uniprot/uniprot_sprot.fasta: inputs/uniprot/uniprot_sprot.fasta.gz
	mkdir -p outputs/uniprot
	cp -a {inputs,outputs}/uniprot/$(<F)
	gunzip -f outputs/uniprot/$(<F)

outputs/uniprot/%.nsq: outputs/uniprot/%
	cd $(@D) && \
	formatdb -i $(<F) -o T -p T

outputs/uniprot/%.namedb: outputs/uniprot/%.fasta
	cd $(@D) && \
	python $(PRJ_ROOT)/scripts/make-namedb.py $(<F) $(@F)

outputs/uniprot/%_screed: outputs/uniprot/%
	python -m screed.fadbm $<

#######################################################################

outputs/chicken_transcripts/global_merged.fa.clean.nr: inputs/chicken_transcripts/global_merged.fa.clean.nr
	mkdir -p $(@D)
	cp -a $< $@

outputs/chicken_transcripts/global_merged.fa.clean.nr_screed: outputs/chicken_transcripts/global_merged.fa.clean.nr
	python -m screed.fadbm $<

outputs/chicken_transcripts/%.part: outputs/chicken_transcripts/%
	cd $(@D) && \
	do-partition.py -x 1e10 -N 4 --threads 4 chick_rna $(<F)

outputs/chicken_transcripts/%.renamed.fasta.gz: outputs/chicken_transcripts/%.part
	cd $(@D) && \
	python $(PRJ_ROOT)/scripts/rename-with-partitions.py chick_rna $*

outputs/chicken_transcripts/%.renamed.fasta: outputs/chicken_transcripts/%.renamed.fasta.gz
	cd $(@D) && \
    gunzip -f $(<F)

outputs/chicken_transcripts/%.nsq: outputs/chicken_transcripts/%
	cd $(@D) && \
	formatdb -i $(<F) -o T -p F

#######################################################################

outputs/rna/pacbio/galGal4/summary.csv: outputs/reference/galGal4.fa_screed outputs/pacbio/galGal4.unmapped_reads workdirs/blat/pacbio_minlen200.h5 outputs/chicken_transcripts/global_merged.fa.clean.nr_screed
	mkdir -p $(@D)
	python scripts/venn.py galGal4 outputs/reference/galGal4.fa_screed /reference_galGal4 outputs/chicken_transcripts/global_merged.fa.clean.nr_screed pacbio workdirs/blat/pacbio_minlen200.h5 outputs/pacbio/galGal4.unmapped_reads

outputs/rna/moleculo/%/summary.csv: outputs/reference/%.fa_screed outputs/moleculo/%.unmapped_reads workdirs/blat/minlen200.h5 outputs/chicken_transcripts/global_merged.fa.clean.nr_screed
	mkdir -p $(@D)
	python scripts/venn.py $* outputs/reference/$*.fa_screed /reference_$* outputs/chicken_transcripts/global_merged.fa.clean.nr_screed moleculo workdirs/blat/minlen200.h5 outputs/moleculo/$*.unmapped_reads

outputs/rna/%/only_rna.fa: outputs/rna/%/summary.csv
outputs/rna/%/only_rna_ref.fa: outputs/rna/%/summary.csv
outputs/rna/%/only_rna_seq.fa: outputs/rna/%/summary.csv
outputs/rna/%/intersection.fa: outputs/rna/%/summary.csv

#######################################################################

outputs/rna/%.part: outputs/rna/%
	cd $(@D) && \
	do-partition.py -x 1e10 -N 4 --threads 4 chick_rna $(<F)

outputs/rna/%.renamed.fasta.gz: outputs/rna/%.part
	cd $(@D) && \
	python $(PRJ_ROOT)/scripts/rename-with-partitions.py chick_rna $*

outputs/rna/%.renamed.fasta: outputs/rna/%.renamed.fasta.gz
	cd $(@D) && \
    gunzip -f $(<F)

outputs/rna/%.nsq: outputs/rna/%
	cd $(@D) && \
	formatdb -i $(<F) -o T -p F

#######################################################################

workdir/results/chick.x.uniprot: outputs/rna/only_rna.fa.renamed.fasta outputs/uniprot/uniprot_sprot.fasta
	mkdir -p $(@D)
	blastall -i $(word 1,$^) -d $(word 2,$^) -e 1e-3 -p blastx -o $@ -a 8 -v 4 -b 4

workdir/results/uniprot.x.chick: outputs/uniprot/uniprot_sprot.fasta outputs/rna/only_rna.fa.renamed.fasta
	mkdir -p $(@D)
	blastall -i $(word 1,$^) -d $(word 2,$^) -e 1e-3 -p tblastn -o $@ -a 8 -v 4 -b 4

outputs/rna/%.homol: outputs/rna/%
	cd $(@D) && \
	python $(PRJ_ROOT)/scripts/make-uni-best-hits.py $(<F) $(@F)

outputs/rna/chick.x.uniprot.ortho: outputs/rna/chick.x.uniprot outputs/rna/uniprot.x.chick
	cd $(@D) && \
	python $(PRJ_ROOT)/scripts/make-reciprocal-best-hits.py chick.x.uniprot uniprot.x.chick $(@F)

outputs/rna/%.annot: outputs/chicken_transcripts/% outputs/rna/chick.x.uniprot.ortho outputs/rna/chick.x.uniprot.homol
	cd $(@D) && \
	python $(PRJ_ROOT)/scripts/annotate-seqs.py $(PRJ_ROOT)/$< chick.x.uniprot.ortho chick.x.uniprot.homol

outputs/rna/ortho.fa: outputs/rna/global_merged.fa.clean.nr.renamed.fasta.annot
	python $(PRJ_ROOT)/scripts/extract.py $< $@ ortho

#######################################################################

outputs/reference/%.nsq: outputs/reference/%
	makeblastdb -in $< -out $< -dbtype nucl

outputs/moleculo/%.00.nsq: outputs/moleculo/%
	makeblastdb -in $< -out $< -dbtype nucl

workdir/results/Chick_RNA_BLAST.txt: outputs/chicken_transcripts/global_merged.fa.clean.nr.renamed.fasta outputs/reference/galGal4.fa.nsq
	blastn -query $< -out $@ -db outputs/reference/galGal4.fa -outfmt 6 -evalue 1e-3 -num_threads 8

workdir/results/Moleculo_RNA_BLAST.txt: outputs/chicken_transcripts/global_merged.fa.clean.nr.renamed.fasta outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fasta.00.nsq
	blastn -query $< -out $@ -db outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fasta -outfmt 6 -evalue 1e-3 -num_threads 16

outputs/rna/match_ortho.rna.ref.txt: outputs/rna/ortho.fa outputs/rna/Chick_RNA_BLAST.txt
	cd $(@D) && \
	python $(PRJ_ROOT)/scripts/find_match_2.py ortho.fa Chick_RNA_BLAST.txt ortho.rna.ref

outputs/rna/match_ortho.rna.moleculo.txt: outputs/rna/ortho.fa outputs/rna/Moleculo_RNA_BLAST.txt
	cd $(@D) && \
	python $(PRJ_ROOT)/scripts/find_match_2.py ortho.fa Moleculo_RNA_BLAST.txt ortho.rna.moleculo

#######################################################################

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

#######################################################################

outputs/reference/%.sa: outputs/reference/%
	bwa index $<

outputs/reference/%.fai: outputs/reference/%
	samtools faidx $<

outputs/reference/galGal4.%: inputs/reference/galGal4.%.gz
	mkdir -p $(@D)
	cp {inputs,outputs}/reference/$(<F)
	gunzip -fN outputs/reference/$(<F)

outputs/reference/galGal3.fa: inputs/reference/galGal3.2bit
	mkdir -p $(@D)
	wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa -P outputs/reference/
	#wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa -P outputs/reference/
	chmod +x outputs/reference/twoBitToFa
	outputs/reference/twoBitToFa $< $@
	-rm outputs/reference/twoBitToFa

outputs/reference/galGal5.fa: $(GALGAL5_SOURCES)
	mkdir -p $(@D)
	scripts/galGal5_fix.sh "$^" > $@

#######################################################################

outputs/msu/msu.fasta: inputs/msu/Chicken.tar.gz
	mkdir -p $(@D)
	tar xf $< -C $(@D)
	gunzip -c $(@D)/Chicken/contig_list.txt.gz > $@.raw
	python scripts/msuseq_fix.py $@.raw $@
	-rm -rf outputs/msu/Chicken/ $@.raw

outputs/msu/%.msu.bam: outputs/msu/msu.fasta outputs/reference/%.fa.sa outputs/reference/%.fa.fai
	mkdir -p $(@D)
	bwa mem outputs/reference/$(*F).fa $< > $<.sam.$(*F)
	samtools import outputs/reference/$(*F).fa.fai $<.sam.$(*F) $@

outputs/msu/%.msu.sorted.bam: outputs/msu/%.msu.bam
	samtools sort $< $(basename $@ .bam)
	samtools index $@

outputs/msu/%.unmapped_reads: outputs/msu/%.msu.sorted.bam
	scripts/extract_reads.sh $< > $@

#######################################################################

outputs/msu/latest.fasta: inputs/msu/chicken_contig_list.all.txt.gz
	mkdir -p $(@D)
	gunzip -c $< > $@.raw
	python scripts/msuseq_fix.py $@.raw $@

outputs/msu/%.latest.bam: outputs/msu/latest.fasta outputs/reference/%.fa.sa outputs/reference/%.fa.fai
	mkdir -p $(@D)
	bwa mem outputs/reference/$(*F).fa $< > $<.sam.$(*F)
	samtools import outputs/reference/$(*F).fa.fai $<.sam.$(*F) $@

outputs/msu/%.latest.sorted.bam: outputs/msu/%.latest.bam
	samtools sort $< $(basename $@ .bam)
	samtools index $@

outputs/msu/%.latest.unmapped_reads: outputs/msu/%.latest.sorted.bam
	scripts/extract_reads.sh $< > $@

#######################################################################

outputs/msu/coverage/%.latest.pd_df.csv: outputs/msu/latest.fasta \
  outputs/reference/galGal4.fa outputs/msu/%.latest.sorted.bam \
  workdirs/msulatest/%/output/output.00500 workdirs/msulatest/%/output/output.01000 \
  workdirs/msulatest/%/output/output.01500 workdirs/msulatest/%/output/output.02000 \
  workdirs/msulatest/%/output/output.02500 workdirs/msulatest/%/output/output.03000 \
  workdirs/msulatest/%/output/output.03500 workdirs/msulatest/%/output/output.03000 \
  workdirs/msulatest/%/output/output.03500 workdirs/msulatest/%/output/output.04000 \
  workdirs/msulatest/%/output/output.04500 workdirs/msulatest/%/output/output.05000 \
  workdirs/msulatest/%/output/output.05500 workdirs/msulatest/%/output/output.06000 \
  workdirs/msulatest/%_90/output/output.00500 workdirs/msulatest/%_90/output/output.01000 \
  workdirs/msulatest/%_90/output/output.01500 workdirs/msulatest/%_90/output/output.02000 \
  workdirs/msulatest/%_90/output/output.02500 workdirs/msulatest/%_90/output/output.03000 \
  workdirs/msulatest/%_90/output/output.03500 workdirs/msulatest/%_90/output/output.03000 \
  workdirs/msulatest/%_90/output/output.03500 workdirs/msulatest/%_90/output/output.04000 \
  workdirs/msulatest/%_90/output/output.04500 workdirs/msulatest/%_90/output/output.05000 \
  workdirs/msulatest/%_90/output/output.05500 workdirs/msulatest/%_90/output/output.06000
	mkdir -p $(@D)
	python scripts/count_reads_pd.py $< workdirs/msulatest/$(*F)/output workdirs/msulatest/$(*F)_90/output $@

workdirs/msulatest/galGal4_90/output/output.%: outputs/reference/galGal4.fa outputs/msu/galGal4.latest.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 --minlen=0.9 1>$@

workdirs/msulatest/galGal4/output/output.%: outputs/reference/galGal4.fa outputs/msu/galGal4.latest.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 1>$@

#######################################################################

outputs/msu/coverage/%.pd_df.csv: outputs/msu/msu.fasta \
  outputs/reference/galGal4.fa outputs/reference/galGal5.fa \
  outputs/msu/%.msu.sorted.bam \
  workdirs/msucov/%/output/output.00500.pbs workdirs/msucov/%/output/output.01000.pbs \
  workdirs/msucov/%/output/output.01500.pbs workdirs/msucov/%/output/output.02000.pbs \
  workdirs/msucov/%/output/output.02500.pbs workdirs/msucov/%/output/output.03000.pbs \
  workdirs/msucov/%/output/output.03500.pbs workdirs/msucov/%/output/output.03000.pbs \
  workdirs/msucov/%/output/output.03500.pbs workdirs/msucov/%/output/output.04000.pbs \
  workdirs/msucov/%/output/output.04500.pbs workdirs/msucov/%/output/output.05000.pbs \
  workdirs/msucov/%/output/output.05500.pbs workdirs/msucov/%/output/output.06000.pbs \
  workdirs/msucov/%_90/output/output.00500.pbs workdirs/msucov/%_90/output/output.01000.pbs \
  workdirs/msucov/%_90/output/output.01500.pbs workdirs/msucov/%_90/output/output.02000.pbs \
  workdirs/msucov/%_90/output/output.02500.pbs workdirs/msucov/%_90/output/output.03000.pbs \
  workdirs/msucov/%_90/output/output.03500.pbs workdirs/msucov/%_90/output/output.03000.pbs \
  workdirs/msucov/%_90/output/output.03500.pbs workdirs/msucov/%_90/output/output.04000.pbs \
  workdirs/msucov/%_90/output/output.04500.pbs workdirs/msucov/%_90/output/output.05000.pbs \
  workdirs/msucov/%_90/output/output.05500.pbs workdirs/msucov/%_90/output/output.06000.pbs
	mkdir -p $(@D)
	python scripts/count_reads_pd.py $< workdirs/msucov/$(*F)/output workdirs/msucov/$(*F)_90/output $@

workdirs/msucov/galGal4_90/output/output.%: outputs/reference/galGal4.fa outputs/msu/galGal4.msu.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 --minlen=0.9 1>$@

workdirs/msucov/galGal5_90/output/output.%: outputs/reference/galGal5.fa outputs/msu/galGal5.msu.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 --minlen=0.9 1>$@

workdirs/msucov/galGal4/output/output.%: outputs/reference/galGal4.fa outputs/msu/galGal4.msu.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 1>$@

workdirs/msucov/galGal5/output/output.%: outputs/reference/galGal5.fa outputs/msu/galGal5.msu.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 1>$@

#######################################################################

#outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq: outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq_screed
#	python -m screed.dump_to_fastq $< $@

outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq: $(MOLECULO_READS)
	cat $^ > $@

outputs/moleculo/%.fasta: outputs/moleculo/%.fastq_screed
	mkdir -p outputs/moleculo
	python -m screed.dump_to_fasta $< $@

outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq_screed: $(addsuffix _screed, $(MOLECULO_READS))
	python scripts/merge_repeated.py $@ $^

outputs/reference/%_screed: outputs/reference/%
	python -m screed.fadbm $<

outputs/moleculo/%_screed: outputs/moleculo/%
	python -m screed.fqdbm $<

outputs/moleculo/%.fastq.sorted.bam: outputs/moleculo/%.fastq.bam
	samtools sort $< $(basename $@ .bam)
	samtools index $@

outputs/moleculo/galGal4.%.fastq.bam: outputs/moleculo/%.fastq outputs/reference/galGal4.fa.sa outputs/reference/galGal4.fa.fai
	bwa mem outputs/reference/galGal4.fa $< > $<.sam.galGal4
	samtools import outputs/reference/galGal4.fa.fai $<.sam.galGal4 $@

outputs/moleculo/galGal3.%.fastq.bam: outputs/moleculo/%.fastq outputs/reference/galGal3.fa.sa outputs/reference/galGal3.fa.fai
	bwa mem outputs/reference/galGal3.fa $< > $<.sam.galGal3
	samtools import outputs/reference/galGal3.fa.fai $<.sam.galGal3 $@

outputs/moleculo/galGal5.%.fastq.bam: outputs/moleculo/%.fastq outputs/reference/galGal5.fa.sa outputs/reference/galGal5.fa.fai
	bwa mem outputs/reference/galGal5.fa $< > $<.sam.galGal5
	samtools import outputs/reference/galGal5.fa.fai $<.sam.galGal5 $@

outputs/moleculo/galGal4.%.unmapped_reads: outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-%.fastq.sorted.bam
	scripts/extract_reads.sh $< > $@

outputs/moleculo/galGal5.%.unmapped_reads: outputs/moleculo/galGal5.LR6000017-DNA_A01-LRAAA-%.fastq.sorted.bam
	scripts/extract_reads.sh $< > $@

outputs/moleculo/%.unmapped_reads: \
  outputs/moleculo/%.1_LongRead_500_1499nt.unmapped_reads \
  outputs/moleculo/%.1_LongRead.unmapped_reads \
  outputs/moleculo/%.2_LongRead_500_1499nt.unmapped_reads \
  outputs/moleculo/%.2_LongRead.unmapped_reads \
  outputs/moleculo/%.3_LongRead_500_1499nt.unmapped_reads \
  outputs/moleculo/%.3_LongRead.unmapped_reads \
  outputs/moleculo/%.4_LongRead_500_1499nt.unmapped_reads \
  outputs/moleculo/%.4_LongRead.unmapped_reads \
  outputs/moleculo/%.5_LongRead_500_1499nt.unmapped_reads \
  outputs/moleculo/%.5_LongRead.unmapped_reads
	cat $^ > $@

#######################################################################

outputs/coverage/%.pd_df.csv: outputs/moleculo/LR6000017-DNA_A01-LRAAA-AllReads.fastq \
  outputs/reference/galGal4.fa outputs/reference/galGal5.fa \
  outputs/moleculo/%.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam \
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

workdirs/galGal4_90/output/output.%: outputs/reference/galGal4.fa outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 --minlen=0.9 1>$@

workdirs/galGal5_90/output/output.%: outputs/reference/galGal5.fa outputs/moleculo/galGal5.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 --minlen=0.9 1>$@

workdirs/galGal4/output/output.%: outputs/reference/galGal4.fa outputs/moleculo/galGal4.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 1>$@

workdirs/galGal5/output/output.%: outputs/reference/galGal5.fa outputs/moleculo/galGal5.LR6000017-DNA_A01-LRAAA-AllReads.sorted.bam
	mkdir -p ${@D}
	bioinfo bam_coverage $^ ${subst output.,,$(@F)} --mapq=30 1>$@

#######################################################################

workdirs/blat/pacbio_minlen200.h5: workdirs/blat/transc_reference_galGal4.pbs \
	workdirs/blat/transc_pacbio.pbs
	mkdir -p $(@D)
	python scripts/blat_merge_outputs.py $@ $(foreach inp,$^,$(subst .pbs,,${inp}))

workdirs/blat/msu_latest_minlen200.h5: workdirs/blat/transc_reference_galGal4.pbs \
  workdirs/blat/transc_msu_latest.pbs
	mkdir -p $(@D)
	python scripts/blat_merge_outputs.py $@ $(foreach inp,$^,$(subst .pbs,,${inp}))

workdirs/blat/msu_minlen200.h5: workdirs/blat/transc_reference_galGal4.pbs \
  workdirs/blat/transc_reference_galGal5.pbs \
  workdirs/blat/transc_msu.pbs
	mkdir -p $(@D)
	python scripts/blat_merge_outputs.py $@ $(foreach inp,$^,$(subst .pbs,,${inp}))

workdirs/blat/minlen200.h5: workdirs/blat/transc_reference_galGal4.pbs \
  workdirs/blat/transc_reference_galGal5.pbs \
  workdirs/blat/transc_moleculo_1_LongRead_500_1499nt.pbs \
  workdirs/blat/transc_moleculo_1_LongRead.pbs \
  workdirs/blat/transc_moleculo_2_LongRead_500_1499nt.pbs \
  workdirs/blat/transc_moleculo_2_LongRead.pbs \
  workdirs/blat/transc_moleculo_3_LongRead_500_1499nt.pbs \
  workdirs/blat/transc_moleculo_3_LongRead.pbs \
  workdirs/blat/transc_moleculo_4_LongRead_500_1499nt.pbs \
  workdirs/blat/transc_moleculo_4_LongRead.pbs \
  workdirs/blat/transc_moleculo_5_LongRead_500_1499nt.pbs \
  workdirs/blat/transc_moleculo_5_LongRead.pbs
	mkdir -p $(@D)
	python scripts/blat_merge_outputs.py $@ $(foreach inp,$^,$(subst .pbs,,${inp}))

workdirs/blat/transc_reference_galGal4: outputs/reference/galGal4.fa outputs/chicken_transcripts/global_merged.fa.clean.nr
	mkdir -p $(@D)
	blat -out=blast8 $^ $@

workdirs/blat/transc_reference_galGal5: outputs/reference/galGal5.fa outputs/chicken_transcripts/global_merged.fa.clean.nr
	mkdir -p $(@D)
	blat -out=blast8 $^ $@

workdirs/blat/transc_moleculo_%: outputs/moleculo/LR6000017-DNA_A01-LRAAA-%.fasta outputs/chicken_transcripts/global_merged.fa.clean.nr
	mkdir -p $(@D)
	blat -out=blast8 $^ $@

workdirs/blat/transc_msu_latest: outputs/msu/latest.fasta outputs/chicken_transcripts/global_merged.fa.clean.nr
	mkdir -p $(@D)
	blat -out=blast8 $^ $@

workdirs/blat/transc_msu: outputs/msu/msu.fasta outputs/chicken_transcripts/global_merged.fa.clean.nr
	mkdir -p $(@D)
	blat -out=blast8 $^ $@

workdirs/blat/transc_pacbio: outputs/pacbio_unmapped/chicken_unmapped.unique.fasta outputs/chicken_transcripts/global_merged.fa.clean.nr
	mkdir -p $(@D)
	blat -out=blast8 $^ $@
	#cat $< | parallel -j 8 --round-robin --pipe --recstart ">" "blat -out=blast8 -noHead $(word 2,$^) stdin >(cat) >&2" > $@

#######################################################################

outputs/rna/msu/galGal4/%.fa.gz: outputs/rna/msu/galGal4/%.fa
	gzip -c $< > $@

outputs/rna/moleculo/galGal4/%.fa.gz: outputs/rna/moleculo/galGal4/%.fa
	gzip -c $< > $@

publish_mrnaseq_msu: $(addprefix outputs/rna/msu/galGal4/,only_rna.fa.gz only_rna_msu.fa.gz only_rna_ref.fa.gz intersection.fa.gz)
	ssh athyra "mkdir -p public_html/mrnaseq/msu/galGal4"
	scp $^ athyra:public_html/mrnaseq/msu/galGal4

publish_mrnaseq_filtered: $(addprefix outputs/rna/cobb/,A-ref_rna.fa B-ref_rna_mol.fa C-mol_rna.fa D-only_rna.fa ortho.fa)
	ssh athyra "mkdir -p public_html/ortho"
	scp $^ athyra:public_html/ortho/

publish_mrnaseq_moleculo: $(addprefix outputs/rna/moleculo/galGal4/,only_rna.fa.gz only_rna_mol.fa.gz only_rna_ref.fa.gz intersection.fa.gz)
	ssh athyra "mkdir -p public_html/mrnaseq/moleculo/galGal4"
	scp $^ athyra:public_html/mrnaseq/moleculo/galGal4

#######################################################################

outputs/uniprot/uniprot.namedb: outputs/uniprot/uniprot_sprot.fasta
	cd outputs/uniprot && formatdb -i uniprot_sprot.fasta -o T -p T
	cd outputs/uniprot && python ../../scripts/make-namedb.py uniprot_sprot.fasta uniprot.namedb
	cd outputs/uniprot && python -m screed.fadbm uniprot_sprot.fasta

#######################################################################

#outputs/pacbio/%.pacbio.bam: outputs/pacbio_assembly/Chicken_10Kb20Kb_40X_Filtered_Subreads.fastq outputs/reference/galGal4.fa.sab
#	mkdir -p $(@D)
#	blasr $< outputs/reference/$(*F).fa -sa outputs/reference/$(*F).fa.sab -sam -nproc 16 -out $<.sam.$(*F)
#	samtools import outputs/reference/$(*F).fa.fai $<.sam.$(*F) $@

#outputs/reference/galGal4.fa.sab: outputs/reference/galGal4.fa
#	sawriter $@ $<

outputs/pacbio/%.pacbio.bam: outputs/pacbio_assembly/Chicken_10Kb20Kb_40X_Filtered_Subreads.fastq outputs/reference/%.fa.sa outputs/reference/%.fa.fai
	mkdir -p $(@D)
	bwa mem -x pacbio outputs/reference/$(*F).fa $< > $<.sam.$(*F)
	samtools import outputs/reference/$(*F).fa.fai $<.sam.$(*F) $@

outputs/pacbio/%.pacbio.sorted.bam: outputs/pacbio/%.pacbio.bam
	samtools sort $< $(basename $@ .bam)
	samtools index $@

outputs/pacbio/%.unmapped_reads: outputs/pacbio/%.pacbio.sorted.bam
	scripts/extract_reads.sh $< > $@

outputs/pacbio/%.unmapped_reads: outputs/pacbio/%.pacbio.sorted.bam
	scripts/extract_reads.py -o $@ $<

#######################################################################

outputs/pacbio_assembly/%_screed: outputs/pacbio_assembly/%
	python -m screed.fqdbm $<

outputs/pacbio/%.fasta: outputs/pacbio_assembly/%.fastq_screed
	mkdir -p outputs/pacbio
	python -m screed.dump_to_fasta $< $@

#######################################################################

outputs/pacbio_assembly/PBcR_Specfile_mer_14.txt: inputs/pacbio_assembly/PBcR_Specfile_mer_14.txt
	mkdir -p $(@D)
	cp -a $< $@

outputs/pacbio_assembly/Chicken_10Kb20Kb_40X_Filtered_Subreads.fastq: inputs/pacbio_assembly/Chicken_10Kb20Kb_40X_Filtered_Subreads.fastq
	mkdir -p $(@D)
	cp -a $< $@

outputs/pacbio_assembly/chicken_2.fasta: outputs/pacbio_assembly/PBcR_Specfile_mer_14.txt \
                                         outputs/pacbio_assembly/Chicken_10Kb20Kb_40X_Filtered_Subreads.fastq
	source hpcc.modules ; \
    cd $(@D) ; \
	PBcR -l chicken_2 -s PBcR_Specfile_mer_14.txt -fastq Chicken_10Kb20Kb_40X_Filtered_Subreads.fastq -maxCoverage 0 genomeSize=1046932099 > >(tee $(@F).stdout.log) 2> >(tee $(@F).stderr.log >&2)

#######################################################################

outputs/pacbio_unmapped/PBcR_Specfile_mer_14.txt: inputs/pacbio_assembly/PBcR_Specfile_mer_14.txt
	mkdir -p $(@D)
	cp -a $< $@

outputs/pacbio_unmapped/galGal4.unmapped_reads: outputs/pacbio/galGal4.unmapped_reads
	mkdir -p $(@D)
	cp -a $< $@

outputs/pacbio_unmapped/chicken_unmapped.fasta: outputs/pacbio_unmapped/PBcR_Specfile_mer_14.txt \
                                                outputs/pacbio_unmapped/galGal4.unmapped_reads
	source hpcc.modules ; \
    cd $(@D) ; \
	PBcR -l chicken_unmapped -s PBcR_Specfile_mer_14.txt -fastq galGal4.unmapped_reads -maxCoverage 0 genomeSize=52346605 -noclean > >(tee $(@F).stdout.log) 2> >(tee $(@F).stderr.log >&2)
#	PBcR -l chicken_unmapped -s PBcR_Specfile_mer_14.txt -fastq galGal4.unmapped_reads -maxCoverage 0 genomeSize=52346605 -noclean -pbcns=0 falconcns=1 > >(tee $(@F).stdout.log) 2> >(tee $(@F).stderr.log >&2)

outputs/pacbio_unmapped/chicken_unmapped.unique.fasta: outputs/pacbio_unmapped/chicken_unmapped.fasta
	awk '{if ($$0 ~ /^>/) {print $$0 "/" NR} else {print $$0 }}' $< > $@


#######################################################################

clean:
	- rm outputs/moleculo/*.sam.galGal*
	- rm outputs/moleculo/*.fastq.bam
	# -rm outputs/galGal4/*.{amb,pac,ann}
	find . -iname "*.pyc" -delete

dependencies.png:
	#make -Bnd | make2graph --gexf > output.gexf
	make -Bnd | make2graph | circo -Tpng -o dependencies.png

#.PRECIOUS: %.sorted.bam
.SECONDARY:

#.PRECIOUS: %.sorted.bam
.SECONDARY:
