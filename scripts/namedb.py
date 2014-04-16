import cPickle
import os
import screed

fp = open('../uniprot/uniprot.namedb')
is_ncbi = cPickle.load(fp)
uniprot_names = cPickle.load(fp)
fp.close()

uniprot_fullname = cPickle.load(open('../uniprot/uniprot.namedb.fullname'))
uniprot_seqs = screed.ScreedDB('../uniprot/uniprot_sprot.fasta')
