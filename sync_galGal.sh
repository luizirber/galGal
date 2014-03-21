#!/bin/bash -exu

rsync -azPv ~/biodata/galGal/ /mnt/scratch/irberlui/biodata/galGal/
rsync -azPv /mnt/scratch/irberlui/biodata/galGal/ ~/biodata/galGal/
