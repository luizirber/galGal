from glob import glob
from itertools import chain

import screed
from screed import fastq
from tqdm import tqdm

def iter_filename(original, filename=None):
    for name in original:
        record = original[name]
        if filename:
            f = "=".join(['File', filename])
            record['name'] = "-".join([record['name'], f])
        yield record


def read_fastq_repeated_names():
    iters = []
    for file in glob('LR6000017-DNA_A01-LRAAA-*LongRead*.fastq_screed'):
        db = screed.ScreedDB(file)

        filename = file.replace('LR6000017-DNA_A01-LRAAA-', "")
        if "500" in filename:
            filename = filename[0] + "_500"
        else:
            filename = filename[0]

        iters.append(tqdm(iter_filename(db, filename)))
        #iterfunc = openscreed.open(filename)
        #iters.append(iter_filename(iterfunc, file))

    screed.create_db('LR6000017-DNA_A01-LRAAA-AllReads.fastq_screed', fastq.FieldTypes, chain(*iters))

read_fastq_repeated_names()
