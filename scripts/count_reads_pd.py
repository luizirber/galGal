#! /usr/bin/env python

from __future__ import print_function

import sys

import numpy as np
import pandas as pd
import screed
from tqdm import tqdm


#query = "LR6000017-DNA_A01-LRAAA-AllReads.fastq"
query = sys.argv[1]


query_lengths = np.fromiter((len(r.sequence)
                             for r in tqdm(screed.open(query), leave=True)),
                            'i8')
print()

data = {}
for i in range(500, 15001, 500):
    data[i] = {}
    reads_over_threshold = query_lengths[query_lengths >= i]

    data[i]['# reads >= minimum'] = len(reads_over_threshold)
    data[i]['# bases in reads >= minimum'] = sum(reads_over_threshold)
df = pd.DataFrame.from_dict(data=data, orient='index')

hist = np.histogram(query_lengths, range(500, 15501, 500))
df['histogram'] = pd.Series(hist[0], hist[1][:-1])

# grep "ref bases covered" bam_cov.*
total_ref_bases_covered = [
  861340070, 856041875, 851206517, 846848062, 842270459, 837297964, 831939901,
  826275521, 820383036, 813324316, 805880858, 797349471, 786997143, 773834101,
  752631701, 700732257, 622487206, 519975737, 397122647, 271267193, 161094016,
  76707337, 25066690, 4830387, 856471, 461040, 421499, 328581, 260590, 191072
]

df['total ref bases covered'] = total_ref_bases_covered

# grep "ref bases covered" output_0.9/output.*|cut -d ":" -f3
total_ref_bases_covered_90 = [
 822849049, 813355067, 804452160, 795848482, 786925399, 777548083, 767605790,
 757371322, 746935144, 735442127, 723596727, 710879582, 696277887, 679037408,
 651879672, 587672522, 501621674, 401902417, 295127743, 194892528, 112797988,
 52316978, 16741243, 3133169, 510501, 270244, 230703, 161312, 125694, 85222
]

df['total ref bases covered, alignment length > 90% read'] = total_ref_bases_covered_90

df.to_csv('pd_df.csv')
