merSize=16
mhap=-k 16 --num-hashes 512 --num-min-matches 3 --threshold 0.04

ovlMemory=32
ovlStoreMemory=32000
threads=32
ovlConcurrency=1
ovlRefBlockSize=20000
ovlCorrBatchSize = 100000
cnsConcurrency=8
merylThreads=32
merylMemory=32000
frgCorrThreads = 16
frgCorrBatchSize = 100000

useGrid=0
scriptOnGrid=0
ovlCorrOnGrid=0
frgCorrOnGrid=0

gridEngine=PBS
gridOptionsScript = -l nodes=1:ppn=1,mem=2GB,walltime=03:00:00 -A ged
gridOptionsOverlapCorrection = -l nodes=1:ppn=1,mem=16GB,walltime=03:00:00 -A ged
gridOptionsConsensus = -l nodes=1:ppn=8,mem=16GB,walltime=03:00:00 -A ged
gridOptionsOverlap = -l nodes=1:ppn=15,mem=16GB,walltime=03:00:00 -A ged
gridOptionsCorrection = -l nodes=1:ppn=15,mem=16GB,walltime=03:00:00 -A ged
gridOptionsFragmentCorrection = -l nodes=1:ppn=16,mem=16GB,walltime=03:00:00 -A ged

#localStaging="/mnt/scratch/tg/irberlui/biodata/galGal/inputs/ecoli/tempK12"
