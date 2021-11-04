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
- [2021.11.04 - Infer transcripts using cufflinks](#20211104---infer-transcripts-using-cufflinks)
- [2021.11.09 - Merge transcript annotation files using cuffmerge](#20211109---merge-transcript-annotation-files-using-cuffmerge)
- [&lt;Template&gt; yyyy.mm.dd - Title](#template-yyyymmdd---title)

# 2021.10.13 - Read counting
[Back to menu](#menu)

## Objectives of the analysis.
Goal is to count total number of reads in fastq files before cleaning for WTC2 R1 and R2 fastq files. 

## Files involved.
WTC2_1.fq.gz, WTC2_2.fq.gz
- WT = Wildtype Candida albicans
- C = one of three biological replicates
- C2 = Thi– = sampled from yeast grown from environment without thiamine
- _1, _2 = forward & reverse reads (?)

## Specific commands used in the analysis.
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

## Results & interpretation.
**There are 20,407,694 reads in each of the two WTC2 fastq files.**

# 2021.10.19 - Trimmomatic
[Back to menu](#menu)

## Objectives of the analysis.
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


## Specific commands used in the analysis.
### Create trim.sbatch
```
$ nano trim.sbatch
```
- [**trim.sbatch**](scripts/trim.sbatch)

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
```
- [**z01.trim_WTC2**](slurm_outputs/z01.trim_WTC2)

```
# Highlighted output:
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

## Results & interpretation.
There are 19,459,631 reads in the trimmed paired-end fastq files.
**95.35% surviving.** 

Per Base Sequence Content greatly improved from before, but still has some visible wobbliness in the first 3 bases. 

High-Duplication sequences do not appear to have been trimmed, which is fine considering this is RNA-seq data.

# 2021.10.21 - Obtain reference genome
[Back to menu](#menu)
## Objectives of the analysis.
Goal is to download refseq C. albicans genome assembly from Entrez genome into a separate folder.

Species: [Candida albicans SC5314 (budding yeasts)](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=237561&lvl=3&lin=f&keep=1&srchmode=1&unlock)
Link at Entrez Genome: [GCA_000182965.3](https://www.ncbi.nlm.nih.gov/assembly/GCA_000182965.3)

## Files involved.
## Specific commands used in the analysis.
```
$ mkdir refseq_GCF_000182965.3
$ cd refseq_GCF_000182965.3
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz
```

## Results & interpretation.
File obtained:
GCF_000182965.3_ASM18296v3_genomic.fna.gz

# 2021.10.26 - bowtie2 sequence alignment
[Back to menu](#menu)

## Objectives of the analysis.
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

## Specific commands used in the analysis.
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
- [**bt2.sbatch**](scripts/bt2.sbatch)

```
# Run bowtie2 for alignment via slurm script
$ sbatch bt2.sbatch

# Wait for completion...

# Check run stats.
$ 31143 --format=jobid,jobname,account,state,elapsed

# Check z01 output:
$ less z01.bt2_WTC2
```
- [**z01.bt2_WTC2**](slurm_outputs/z01.bt2_WTC2)

## Results & interpretation.
90.52% of the input (paired) reads aligned concordantly exactly 1 time.

97.97% overall alignment rate

The bowtie2 job took 01:37:18 (1.5+ hours) to complete.

[See class's alignment results.](https://docs.google.com/spreadsheets/d/1GZOHIUhfYA523kkOdoF9APQEkRsberddvEsEZOkIIvA/edit)


# 2021.10.28 - Tophat alignment
[Back to menu](#menu)

## Objectives of the analysis.
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
        - The sequences that were successfully aligned by tophat to the reference genome.
    - **align_summary.txt**
        - A summary of the percentage of input reads that were aligned, either concordantly or discordantly, to the reference genome.
    - deletions.bed  
    - insertions.bed   
    - junctions.bed  
    - logs/
    - prep_reads.info
    - unmapped.bam

**Bolded files** are important.

## Specific commands used in the analysis.
```
$ gunzip GCF_000182965.3_ASM18296v3_genomic.gff.gz

$ nano tophat_align.sbatch
```
- [**tophat_align.sbatch**](scripts/tophat_align.sbatch)

```
$ sbatch tophat_align.sbatch
jobid: 37967

# wait for some hours…

$ less z01.tophat_WTC2
# last two lines of z01.tophat_WTC2 are:

[2021-10-28 13:03:03] A summary of the alignment counts can be found in WTC2_tophat/align_summary.txt
[2021-10-28 13:03:03] Run complete: 02:52:28 elapsed

$ cd WTC2_tophat

$ less align_summary.txt
```
- [**align_summary.txt**](summary_outputs/align_summary.txt)

## Results & interpretation.
Tophat spliced alignment resulted in 90.3% concordant pair alignment rate. 
Reads aligned to the reference genome were collected in the accepted_hits.bam file.

[See class's alignment results.](https://docs.google.com/spreadsheets/d/1GZOHIUhfYA523kkOdoF9APQEkRsberddvEsEZOkIIvA/edit)

# 2021.11.02 - Transfer to GitHub
[Back to menu](#menu)

## Objectives of the analysis.
Transfer notes to GitHub.

## Files involved.
Repository
- lil-qorgi:379-rnaseq

READ<span>ME.md</span>
- Notes document for RNAseq workflow.

## Specific commands used in the analysis.
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

## Results & interpretation.
The journal has been transferred onto Github.


# 2021.11.04 - Infer transcripts using cufflinks
[Back to menu](#menu)

## Objectives of the analysis.
The goal is to use cufflinks to take in the tophat read alignment results to infer the transcripts that are found in the sequenced samples. 

See [cufflinks manual](http://cole-trapnell-lab.github.io/cufflinks/manual/)

## Files involved.
### —cufflinks—
### *input files*
- GCF_000182965.3_ASM18296v3_genomic.gff
    - The reference genome annotation file
        - Note that, even though cufflinks manual says it requires .gtf., .gff works for cufflinks as well because the two file formats are very similar.
- accepted_hits.bam
    - The sequences that were successfully aligned by tophat to the reference genome from the previous step.
### *output files*
- cufflinks_output/
    - genes.fpkm_tracking
    - isoforms.fpkm_tracking
    - skipped.gtf
    - **transcripts.gtf**
        - this is the annotations file containing the transcripts inferred from the accepted_hits.bam alignment file.

## Specific commands used in the analysis.
```
# First, pull all required files together into the same folder for easier usage.

# Then write the slurm script for cufflinks
$ nano cufflinks.sbatch
```
- [**cufflinks.sbatch**](scripts/cufflinks.sbatch)

```
# Run the sbatch script
$ sbatch cufflinks.sbatch

$ squeue -u qz108

# Wait for job to complete.

$ sacct -j 44135 --format=jobid,jobname,account,state,elapsed
# The job took 05:48 minutes to complete.

# See the sbatch output
$ less z01.cufflinks_WTC2
```
- [**z01.cufflinks_WTC2**](slurm_outputs/z01.cufflinks_WTC2)

```
# One of the last lines says:
Processed 6131 loci.

$ cd cufflinks_output/
$ ls -l

-rw-r--r-- 1 qz108 users  739513 Nov  2 19:59 genes.fpkm_tracking
-rw-r--r-- 1 qz108 users  756069 Nov  2 19:59 isoforms.fpkm_tracking
-rw-r--r-- 1 qz108 users       0 Nov  2 19:55 skipped.gtf
-rw-r--r-- 1 qz108 users 2958558 Nov  2 19:59 transcripts.gtf

$ less transcripts.gft
# The file looks good.

$ wc -l transcripts.gtf
# There are 12994 lines in the file. There is no header in the file.
```

## Results & interpretation.
The cufflinks program successfully inferred transcripts based on the tophat alignment results. 

# 2021.11.09 - Merge transcript annotation files using cuffmerge
[Back to menu](#menu)

## Objectives of the analysis.
The goal is to merge the transcript annotations from all biological replicates to obtain the merged annotation file merged.gtf.

## Files involved.
### —cuffmerge—
### *input files*
- GCF_000182965.3_ASM18296v3_genomic.gff
    - Note that, similar to cufflinks, even though cuffmerge  requires .gtf., .gff works as well because the two file formats are very similar.
- transcript gtf's, resulting from running previous steps from Trimmomatic to tophat (not including bowtie2) over each original raw reads file from the sequencing run.
    - WTA1_transcripts.gtf
    - WTA2_transcripts.gtf
    - WTB1_transcripts.gtf
    - WTB2_transcripts.gtf
    - WTC1_transcripts.gtf
    - WTC2_transcripts.gtf
- file specifying the transcript gtfs
    - transcripts_gtf.txt
        ```
        WTA1_transcripts.gtf
        WTA2_transcripts.gtf
        WTB1_transcripts.gtf
        WTB2_transcripts.gtf
        WTC1_transcripts.gtf
        WTC2_transcripts.gtf
        ```

### *output files*
- cuffmerge_output/
    - logs/
        - run.log
    - **merged.gtf**
        - this is the annotation file that merges the transcripts.gtf file from the previous (cufflinks) step obtained for all biological replicates: WTA1, A2, B1, B2, C1, C2. 

## Specific commands used in the analysis.
```
# First, move all required input files into the same folder.
# Obtain all transcripts.gtf files from Google Bucket:
$ gsutil cp gs://gu-biology-dept-class/*.gtf 

# Create a text file specifying the transcript.gtf files
$ nano transcripts_gtf.txt
```
```
WTA1_transcripts.gtf
WTA2_transcripts.gtf
WTB1_transcripts.gtf
WTB2_transcripts.gtf
WTC1_transcripts.gtf
WTC2_transcripts.gtf
```
```
# Create an sbatch file for cuffmerge.
$ nano cuffmerge.sbatch
```
- [**cuffmerge.sbatch**](scripts/cuffmerge.sbatch)

```
# check folder contents
$ ls -1
cuffmerge.sbatch
GCF_000182965.3_ASM18296v3_genomic.gff
transcripts_gtf.txt
WTA1_transcripts.gtf
WTA2_transcripts.gtf
WTB1_transcripts.gtf
WTB2_transcripts.gtf
WTC1_transcripts.gtf
WTC2_transcripts.gtf

# run sbatch script
sbatch cuffmerge.sbatch

# wait for sbatch script to finish

# check the sbatch output file
$ less z01.cuffmerge
```
- [**z01.cuffmerge**](slurm_outputs/z01.cuffmerge)

```
# go to cuffmerge_output
$ cd cuffmerge_output

# check contents
$ ll
drwxr-xr-x 2 qz108 users       0 Nov  4 10:34 logs
-rw-r--r-- 1 qz108 users 1577746 Nov  4 10:35 merged.gtf

# check job elapsed time
$ sacct -j 47686 --format=jobname,jobid,user,elapsed
# the process took 35 seconds.
```

## Results & interpretation.
The process took 35 seconds.

Cuffmerge successfully merged the transcript annotation (.gtf) files and produced the merged.gtf file.

<br></br>
# —Template—
<br></br>

# yyyy.mm.dd - Title
[Back to menu](#menu)

## Objectives of the analysis.

## Files involved.

## Specific commands used in the analysis.

## Results & interpretation.
