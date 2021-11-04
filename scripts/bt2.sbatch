#!/bin/bash
#SBATCH --job-name=bt2_WTC2 --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=netID@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

# Bowtie2

module load bowtie2/2.4.4

bowtie2 -x Calbicans -1 WTC2_1.trPE.fq -2 WTC2_2.trPE.fq -S WTC2.sam

module unload bowtie2/2.4.4