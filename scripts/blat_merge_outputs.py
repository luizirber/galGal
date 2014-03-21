#! /usr/bin/env python

from glob import glob
import os

import pandas as pd


COLUMN_NAMES = [s.strip()
                for s in "Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end , e-value, bit score".split(',')]

output = pd.DataFrame()

for output_file in glob("../workdir/results/transc_*.filtered"):
    print 'reading file', output_file
    table = pd.read_table(output_file, names=COLUMN_NAMES)
    source = "_".join(os.path.splitext(os.path.basename(output_file))[0].split('_')[1:])
    table['source'] = source
#    print 'appending source', source
#    output = pd.concat([output, table])
#    print 'finished appending source', source
    table[table['alignment length'] > 200].to_hdf("minlen200.h5", source)
    #table.to_hdf("merged.h5", source)#, mode='w', complib='blosc', complevel=9)


#output.to_csv('merged.csv')
#print 'writing HDF5'
#output.to_hdf('merged.h5', 'merged.h5', mode='w', complib='blosc', complevel=9)
