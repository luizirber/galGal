PBCR_RES=walltime=1:00:00,nodes=1:ppn=32,mem=20gb -A ged

K12.fasta.pbs:
	JOBID=`echo make $(subst .pbs,,$@) | cat header.sub - footer.sub | \
	  qsub -l ${PBCR_RES} -N pbcr.${subst output.,,$(@F)} -o $@ -e $@.err | cut -d"." -f1` ; \
	while [ -n "$$(qstat -a |grep $${JOBID})" ]; do sleep 600; done
	@grep "galGal PBS job finished: SUCCESS" $@

K12.fasta:
	PBcR -l K12 -s pacbio.spec -fastq ecoli_filtered.fastq genomeSize=4650000
