#! /usr/bin/env python

from glob import glob
import os

import pandas as pd


COLUMN_NAMES = [s.strip()
                for s in "Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end , e-value, bit score".split(',')]

output = pd.DataFrame()

for output_file in glob("../workdir/results/transc_*.filtered"):
    table = pd.read_table(output_file, names=COLUMN_NAMES)
    table['source'] = "_".join(os.path.splitext(os.path.basename(output_file))[0].split('_')[1:])
    output = output.append(table)

#output.to_csv('merged.csv')
output.to_hdf('merged.h5', keys=output.keys(), mode='w', complib='blosc', complevel=9)
