#!/bin/bash
#SBATCH --job-name=cuffmerge --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=qz108@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=40G

# Variables
output_dir=cuffmerge_output
reference_anno=GCF_000182965.3_ASM18296v3_genomic.gff
gtf_list=transcripts_gtf.txt

# Cuffmerge
module load anaconda2
module load cufflinks

cuffmerge \
-o $output_dir \
-g $reference_anno \
$gtf_list

module unload cufflinks
module unload anaconda2 