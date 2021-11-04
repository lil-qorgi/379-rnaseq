#!/bin/bash
#SBATCH --job-name=tophat_WTC2 --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=netID@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=40G

# Tophat

module load tophat

tophat -o WTC2_tophat \
-G GCF_000182965.3_ASM18296v3_genomic.gff \
Calbicans \
WTC2_1.trPE.fq WTC2_2.trPE.fq

module unload tophat