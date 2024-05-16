#!/bin/bash

#SBATCH -t 1:00:00
#SBATCH -A bhatele-lab-cmsc
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
# #SBATCH -o human-1.out
#SBATCH -o human-2.out

export PARLAY_NUM_THREADS=128

#../bin/caps_sa Homo_sapiens.GRCh38.dna.toplevel.fa human_genome.out
../bin/caps_sa GRCh38_latest_genomic.fna.gz human_genome2.out
