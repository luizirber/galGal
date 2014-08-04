#! /usr/bin/env python

from glob import glob
import sys
import os

import pandas as pd


COLUMN_NAMES = [s.strip()
                for s in "Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end , e-value, bit score".split(',')]

def blat_merge_outputs(output_file, results):
    for result in results:
        print 'reading file', result
        table = pd.read_table(result, names=COLUMN_NAMES)
        source = "_".join(os.path.splitext(os.path.basename(result))[0].split('_')[1:])
        table['source'] = source
        table[table['alignment length'] > 200].to_hdf(output_file, source)


if __name__ == "__main__":
    output = sys.argv[1]
    inputs = sys.argv[2:]
    blat_merge_outputs(output, inputs)
