#! /usr/bin/env python

from __future__ import print_function

import os
import sys

import numpy as np
import pandas as pd
import screed


#query = "LR6000017-DNA_A01-LRAAA-AllReads.fastq"
query = sys.argv[1]


query_lengths = np.fromiter((len(r.sequence) for r in screed.open(query)), 'i8')

data = {}
total_ref_bases_covered = {}
total_ref_bases_covered_90 = {}

for i in range(500, 15001, 500):
    data[i] = {}
    reads_over_threshold = query_lengths[query_lengths >= i]

    data[i]['# reads >= minimum'] = len(reads_over_threshold)
    data[i]['# bases in reads >= minimum'] = sum(reads_over_threshold)

    with open(os.path.join(sys.argv[2], 'output.%05d' % i), 'r') as f:
        total_ref_bases_covered[i] = None
        for line in f:
            if 'covered' in line:
                total_ref_bases_covered[i] = int(line.strip().split(':')[-1])

    with open(os.path.join(sys.argv[3], 'output.%05d' % i), 'r') as f:
        total_ref_bases_covered_90[i] = None
        for line in f:
            if 'covered' in line:
                total_ref_bases_covered_90[i] = int(line.strip().split(':')[-1])


df = pd.DataFrame.from_dict(data=data, orient='index')
hist = np.histogram(query_lengths, range(500, 15501, 500))
df['histogram'] = pd.Series(hist[0], hist[1][:-1])
df['total ref bases covered'] = pd.Series(total_ref_bases_covered, index=total_ref_bases_covered.keys())
df['total ref bases covered, alignment length > 90% read'] = pd.Series(total_ref_bases_covered_90, index=total_ref_bases_covered_90.keys())

df.to_csv(sys.argv[4])
