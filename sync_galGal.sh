#!/bin/bash -exu

rsync -azPv ~/biodata/galGal/ /mnt/scratch/irberlui/biodata/galGal/ -f"- .git/"
rsync -azPv /mnt/scratch/irberlui/biodata/galGal/ ~/biodata/galGal/ -f"- .git/"
