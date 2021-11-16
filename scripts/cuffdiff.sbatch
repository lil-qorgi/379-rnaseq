#!/bin/bash
#SBATCH --job-name=cuffdiff --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=qz108@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=40G

# Variables
input_dir=cuffdiff_input
output_dir=cuffdiff_output
bam_suffix=accepted_hits.bam

# Cuffmerge
module load cufflinks

cuffdiff \
-o $output_dir \
$input_dir/merged.gtf \
$input_dir/WTA1_$bam_suffix,$input_dir/WTB1_$bam_suffix,$input_dir/WTC1_$bam_suffix \
$input_dir/WTA2_$bam_suffix,$input_dir/WTB2_$bam_suffix,$input_dir/WTC2_$bam_suffix

module unload cufflinks
