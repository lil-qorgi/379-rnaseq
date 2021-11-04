#!/bin/bash
#SBATCH --job-name=trim_WTC2 --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=<netID>@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

# Author: QZ            #

# Set up variables	#
trim_prog=/home/qz108/Trimmomatic/trimmomatic-0.39.jar

in_R1=WTC2_1.fq.gz
in_R2=WTC2_2.fq.gz

out_R1_PE=WTC2_1.trPE.fq.gz
out_R1_SE=WTC2_1.trSE.fq.gz

out_R2_PE=WTC2_2.trPE.fq.gz
out_R2_SE=WTC2_2.trSE.fq.gz

adapters=/home/qz108/Trimmomatic/adapters/TruSeq3-PE.fa

# Run Trimmomatic	#

java -jar $trim_prog PE \
$in_R1 \
$in_R2 \
$out_R1_PE \
$out_R1_SE \
$out_R2_PE \
$out_R2_SE \
ILLUMINACLIP:$adapters:2:30:10 \
HEADCROP:10 \
TRAILING:10 \
SLIDINGWINDOW:4:15 \
MINLEN:75