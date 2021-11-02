# *C. albicans* RNAseq Workflow Notes 

Name&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Q Zhang
## Version
v0.1.0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2021.10.12-10.31
- Entries logged using Google Docs.

v1.0.0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2021.11.02-??
- Entries transferred to and logged using GitHub.
- Repo: lil-qorgi:379-rnaseq

## *Requirements for this document
- Document everything that someone would need to follow to recapitulate your analysis:
    - downloading files, software, commands, scripts, analyses runs, etc.
    - explain all file names
    - [see example notebook](https://docs.google.com/document/d/1KJAZZpvzMqlzfekco8q_vjajGqBRnzcg/edit)

# Menu
- [2021.10.13 - Read counting](#20211013---read-counting)
- [2021.10.19 - Trimmomatic](#20211019---trimmomatic)
- [2021.10.21 - Obtain reference genome](#20211021---obtain-reference-genome)
- [2021.10.26 - bowtie2 sequence alignment](#20211026---bowtie2-sequence-alignment)
- [2021.10.28 - Tophat alignment](#20211028---tophat-alignment)
- [2021.11.02 - Transfer to GitHub](#20211102---transfer-to-github)
- [&lt;Template&gt; yyyy.mm.dd - Title](#template-yyyymmdd---title)

# 2021.10.13 - Read counting
[Back to menu](#menu)

## Summary description of objectives / goals of the analysis.
Goal is to count total number of reads in fastq files before cleaning for WTC2 R1 and R2 fastq files. 

## Files involved.
WTC2_1.fq.gz, WTC2_2.fq.gz
- WT = Wildtype Candida albicans
- C = one of three biological replicates
- C2 = Thi– = sampled from yeast grown from environment without thiamine
- _1, _2 = forward & reverse reads (?)

## Specific commands used in analyses. (Copy & paste.)
```
# unzip the .gz files
$ gunzip WTC*

# count the number of lines in _1.fq & divide by 4
$ wc -l WTC2_1.fq
81630776
$ bc -l <<< '81630776/4'
20407694

# count the number of lines in _2.fq & divide by 4
$ wc -l WTC2_2.fq
81630776
$ bc -l <<< '81630776/4'
20407694
```

## Summary description of results & interpretation.
**There are 20,407,694 reads in each of the two WTC2 fastq files.**

# 2021.10.19 - Trimmomatic
[Back to menu](#menu)

## Summary description of objectives / goals of the analysis.
The goal is to use Trimmomatic to trim the raw reads so as to remove the problematic first ten bases of each read, reduce adapter content, improve reverse read quality, and do an overall trimming via sliding window. 

Review trimmed PE files in FastQC to ensure Trim achieved aims as intended.

## Files involved.
Input files into Trimmomatic:
- WTC2_1.fq.gz, WTC2_2.fq.gz
- For descriptions, see previous entry.

Output files from Trimmomatic:
- Paired-end and Single-end trimmed sequence files for forward reads (R1)
    - WTC2_1.trPE.fq.gz
    - WTC2_1.trSE.fq.gz
- Paired-end and Single-end trimmed sequence files for reverse reads (R2)
    - WTC2_2.trPE.fq.gz
    - WTC2_2.trSE.fq.gz
- *Note that the trSE versions will not be used in subsequent steps.


## Specific commands used in analyses. (Copy & paste.)
### Create trim.sbatch
```
$ nano trim.sbatch
```
Inside **trim.sbatch**
```
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
```

### Run sbatch Trimmomatic command
```
$ sbatch trim.sbatch

# to check on sbatch status
$ squeue -u qz108 
# job id was 24813

# wait for sbatch to finish…
# see low long it took
$ sacct -j 24813 --format=Elapsed
39:05
```

### Observe output on GCP
```
# observe output
$ less z01.trim_WTC2
Input Read Pairs: 20407694 Both Surviving: 19459631 (95.35%) Forward Only Surviving: 599344 (2.94%) Reverse Only Surviving: 218664 (1.07%) Dropped: 130055 (0.64%)
# This already gives the answer of 19459631 reads. Let's check by unzipping.

# copy & store the .gz files in a local folder for later use.
# unzip output files and check number of reads
$ gunzip WTC2*trPE*.gz
$ wc -l WTC2_1.trPE.fq
77838524

$ bc -l <<< '/4'
19459631
```

### Observe output on local
```
# download (with custom alias)
get_hpc /home/qz108/RNA_seq_workflow/WTC2*trPE*
# open with FastQC & observe.
```

## Summary description of results & interpretation.
There are 19,459,631 reads in the trimmed paired-end fastq files.
**95.35% surviving.** 

Per Base Sequence Content greatly improved from before, but still has some visible wobbliness in the first 3 bases. 

High-Duplication sequences do not appear to have been trimmed, which is fine considering this is RNA-seq data.

# 2021.10.21 - Obtain reference genome
[Back to menu](#menu)
## Summary description of objectives / goals of the analysis.
Goal is to download refseq C. albicans genome assembly from Entrez genome into a separate folder.

Species: [Candida albicans SC5314 (budding yeasts)](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=237561&lvl=3&lin=f&keep=1&srchmode=1&unlock)
Link at Entrez Genome: [GCA_000182965.3](https://www.ncbi.nlm.nih.gov/assembly/GCA_000182965.3)

## Files involved.
## Specific commands used in analyses. (Copy & paste.)
```
$ mkdir refseq_GCF_000182965.3
$ cd refseq_GCF_000182965.3
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz
```

## Summary description of results & interpretation.
File obtained:
GCF_000182965.3_ASM18296v3_genomic.fna.gz

# 2021.10.26 - bowtie2 sequence alignment
[Back to menu](#menu)

## Summary description of objectives / goals of the analysis.
Run a bowtie2 (non-spliced) alignment of trimmed reads to the reference genome for C. albicans. Mostly, this bowtie2 run results will be replaced by the tophat results later on. This is mostly to familiarize with bowtie2 usage.

## Files involved.
### —bowtie2-build—
### *input files*
- GCF_000182965.3_ASM18296v3_genomic.fna.gz
    - It's the C. albicans reference genome that will be used.
### *output files*
Index files
- Calbicans.1.bt2
- Calbicans.2.bt2
- Calbicans.3.bt2
- Calbicans.4.bt2
- Calbicans.rev.1.bt2
- Calbicans.rev.2.bt2

### —bowtie2—
### *input files*
Index files
- Calbicans.1.bt2
- Calbicans.2.bt2
- Calbicans.3.bt2
- Calbicans.4.bt2
- Calbicans.rev.1.bt2
- Calbicans.rev.2.bt2
Trimmed FastQ paired-end files (from Trimmomatic) to align
- WTC2_1.trPE.fq
    - forward PE cleaned reads
- WTC2_2.trPE.fq
    - reverse PE cleaned reads
### *output files*
Sequence alignment map
- WTC2.sam

## Specific commands used in analyses. (Copy & paste.)
```
# Create & navigate to directory
$ mkdir 1021-26_bowtie2_alignment
$ cd 1021-26_bowtie2_alignment
```
### —bowtie2-build—
```
# This is to build an index of the reference genome, which is then used by bowtie2 for alignment.

# Enter a compute node
$ srun --pty bash

# Do index run on reference genome
$ bowtie2-build GCF_000182965.3_ASM18296v3_genomic.fna.gz Calbicans
```
### —bowtie2—
```
# Copy over the WTC2 trPE fastq files from the outputs of the Trimmomatic run.
$ cp ../1019_Trimmomatic/WTC2*trPE.fq .

# Create sbatch file for bowtie2

$ nano bt2.sbatch
```
—bt2.sbatch—
```
#!/bin/bash
#SBATCH --job-name=bt2_WTC2 --output=z01.%x
#SBATCH --mail-type=END,FAIL --mail-user=netID@georgetown.edu
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --time=72:00:00
#SBATCH --mem=4G

# Bowtie2

module load bowtie2/2.4.4

bowtie2 -x Calbicans -1 WTC2_1.trPE.fq -2 WTC2_2.trPE.fq -S WTC2.sam

module unload bowtie2/2.4.4
```
——
```
# Run bowtie2 for alignment via slurm script
$ sbatch bt2.sbatch

# Wait for completion...

# Check run stats.
$ 31143 --format=jobid,jobname,account,state,elapsed

# Check z01 output:
$ less z01.bt2_WTC2
```
—z01.bt2_WTC2—
```
19459631 reads; of these:
  19459631 (100.00%) were paired; of these:
    942014 (4.84%) aligned concordantly 0 times
    17615558 (90.52%) aligned concordantly exactly 1 time
    902059 (4.64%) aligned concordantly >1 times
    ----
    942014 pairs aligned concordantly 0 times; of these:
      300558 (31.91%) aligned discordantly 1 time
    ----
    641456 pairs aligned 0 times concordantly or discordantly; of these:
      1282912 mates make up the pairs; of these:
        791009 (61.66%) aligned 0 times
        448352 (34.95%) aligned exactly 1 time
        43551 (3.39%) aligned >1 times
97.97% overall alignment rate
```
——


## Summary description of results & interpretation.
90.52% of the input (paired) reads aligned concordantly exactly 1 time.

The bowtie2 job took 01:37:18 (1.5+ hours) to complete.

[See class's alignment results.](https://docs.google.com/spreadsheets/d/1GZOHIUhfYA523kkOdoF9APQEkRsberddvEsEZOkIIvA/edit)


# 2021.10.28 - Tophat alignment
[Back to menu](#menu)

## Summary description of objectives / goals of the analysis.
Perform spliced alignment using tophat. 

## Files involved.
### —Tophat—
### *input files*
Annotation file (unzipped) for identifying splice sites
- GCF_000182965.3_ASM18296v3_genomic.gff

Index files built by bt2
- Calbicans.&lt;x&gt;.bt2 
    - (see previous entry)

Paired-end inputs
- WTC2_1.trPE.fq
- WTC2_2.trPE.fq

### *output files*
Folder containing all outputs
- WTC2_tophat/
    - **accepted_hits.bam**
    - **align_summary.txt**
    - deletions.bed  
    - insertions.bed   
    - junctions.bed  
    - logs/
    - prep_reads.info
    - unmapped.bam

**Bolded files** are important.

## Specific commands used in analyses. (Copy & paste.)
```
$ gunzip GCF_000182965.3_ASM18296v3_genomic.gff.gz

$ nano tophat_align.sbatch
```
—tophat_align.sbatch—
```
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
```
——
```
$ sbatch tophat_align.sbatch
jobid: 37967

# wait for some hours…

$ less z01.tophat_WTC2
# last two lines:
[2021-10-28 13:03:03] A summary of the alignment counts can be found in WTC2_tophat/align_summary.txt
[2021-10-28 13:03:03] Run complete: 02:52:28 elapsed

$ cd WTC2_tophat

$ less align_summary.txt
```
—align_summary.txt—
```
Left reads:
          Input     :  19459631
           Mapped   :  18359897 (94.3% of input)
            of these:    418573 ( 2.3%) have multiple alignments (216 have >20)
Right reads:
          Input     :  19459631
           Mapped   :  18213032 (93.6% of input)
            of these:    414252 ( 2.3%) have multiple alignments (215 have >20)
94.0% overall read mapping rate.

Aligned pairs:  17584504
     of these:    401012 ( 2.3%) have multiple alignments
                   15884 ( 0.1%) are discordant alignments
90.3% concordant pair alignment rate.
```
——


## Summary description of results & interpretation.
Tophat spliced alignment resulted in 90.3% concordant pair alignment rate. 

[See class's alignment results.](https://docs.google.com/spreadsheets/d/1GZOHIUhfYA523kkOdoF9APQEkRsberddvEsEZOkIIvA/edit)

# 2021.11.02 - Transfer to GitHub
[Back to menu](#menu)

## Summary description of objectives / goals of the analysis.
Transfer notes to GitHub.

## Files involved.
Repository
- lil-qorgi:379-rnaseq

READ<span>ME.md</span>
- Notes document for RNAseq workflow.

## Specific commands used in analyses. (Copy & paste.)
```
# Create local repo.
# (in parent folder:)
$ git clone https://github.com/lil-qorgi/379-rnaseq.git RNAseq_notes
$ cd RNAseq_notes
# make changes to README.md
# add & commit changes
$ git add .
$ git commit -m "Transferred journal to readme."
$ git push
```

## Summary description of results & interpretation.

<br></br>
# &lt;Template&gt; yyyy.mm.dd - Title

## Summary description of objectives / goals of the analysis.

## Files involved.

## Specific commands used in analyses. (Copy & paste.)

## Summary description of results & interpretation.