#! /usr/bin/env python

from __future__ import print_function

import sys

import numpy as np
import screed
from tqdm import tqdm


query = "LR6000017-DNA_A01-LRAAA-AllReads.fastq"


query_lengths = np.fromiter((len(r.sequence)
                             for r in tqdm(screed.open(query), leave=True)),
                            'i8')
print()

for i in range(500, 15001, 500):
    reads_over_threshold = query_lengths[query_lengths > i]
    print(len(reads_over_threshold), sum(reads_over_threshold))

hist = np.histogram(query_lengths, range(500, 15501, 500))
for i in hist[0]:
    print(i)
