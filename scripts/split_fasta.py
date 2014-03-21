#! /usr/bin/env python

import screed

db = screed.read_fasta_sequences('galGal4.fa.masked')

keys = [k for k in db.keys() if "_" not in k and "Un" not in k]

filtered = screed.fasta.FASTA_Writer('galGal4.fa.masked.filtered')
for k in keys:
    record = db[k]
    filtered.write(record)

for k in keys:
    record = db[k]
    filtered = screed.fasta.FASTA_Writer(k)
    filtered.write(record)
