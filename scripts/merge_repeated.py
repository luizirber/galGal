#! /usr/bin/env python

from itertools import chain
import os
import sys

import screed
from screed import fastq

def iter_filename(original, filename=None):
    for name in original:
        record = original[name]
        if filename:
            f = "=".join(['File', filename])
            record['name'] = "-".join([record['name'], f])
        yield record


def read_fastq_repeated_names(output, inputs):
    iters = []
    for file in inputs:
        db = screed.ScreedDB(file)

        filename = os.path.basename(file).replace('LR6000017-DNA_A01-LRAAA-', "")
        if "500" in filename:
            filename = filename[0] + "_500"
        else:
            filename = filename[0]

        iters.append(iter_filename(db, filename))

    screed.create_db(output, fastq.FieldTypes, chain(*iters))


if __name__ == '__main__':
    output = sys.argv[1]
    inputs = sys.argv[2:]
    read_fastq_repeated_names(output, inputs)
