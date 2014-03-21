#!/bin/bash -exu

rsync -azPv ~/biodata/galGal/ /mnt/scratch/irberlui/biodata/galGal/ --exclude={*/.git}
rsync -azPv /mnt/scratch/irberlui/biodata/galGal/ ~/biodata/galGal/ --exclude={*/.git}
