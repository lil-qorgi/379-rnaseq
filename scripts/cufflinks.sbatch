#!/bin/bash
#SBATCH --job-name=cufflinks_WTC2 --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=qz108@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=40G

# Cufflinks

module load cufflinks

cufflinks -o cufflinks_output --GTF GCF_000182965.3_ASM18296v3_genomic.gff accepted_hits.bam

module unload cufflinks