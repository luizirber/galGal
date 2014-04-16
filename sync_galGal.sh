#!/bin/bash -exu

# Reference:
# http://www.gnu.org/software/parallel/man.html#example__parallelizing_rsync

export SRCDIR="/mnt/home/irberlui/biodata/galGal"
export DESTDIR="/mnt/scratch/tg/irberlui/biodata/galGal"

cd $SRCDIR; find . -type f -path ./misc -prune -o -size +100000 | parallel -j 10 -v mkdir -p $DESTDIR/{//}\;rsync -avP {} $DESTDIR/{}

rsync -av $SRCDIR/ $DESTDIR/ -f"- .git/"

rsync -av $DESTDIR/ $SRCDIR/ -f"- .git/"
