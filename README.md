# *C. albicans* RNAseq Workflow Notes 

Name&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Q Zhang

Course&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;BIOL 379

## Version
v0.1.0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2021.10.12-10.31
- Entries logged using Google Docs.

v1.0.0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2021.11.02-??
- Entries transferred to and logged using GitHub.
- Repo: lil-qorgi:379-rnaseq

## (Requirements for this document)
- Document everything that someone would need to follow to recapitulate your analysis:
    - downloading files, software, commands, scripts, analyses runs, etc.
    - explain all file names
    - [see example notebook (private link)](https://docs.google.com/document/d/1KJAZZpvzMqlzfekco8q_vjajGqBRnzcg/edit)

# Introduction & Background

The Rolfes lab is conducting a study to identify possible genes that regulate thiamine biosynthesis in the yeast *Candida albicans*. *Candida albicans* is a normally commensal species commonly found on humans, but can become pathogenic under certain environmental conditions, leading to the candidiasis disease. By examining the genes that are significantly differentially expressed after removing thiamine from the growing medium, the team hopes to better understand the role of thiamine in the normal growth as well as the pathogenesis of *C. albicans*. 

The Rolfes lab has sequenced wildtype *Candida albicans* yeast cells grown in two media differing only in thiamine presence — with the control having thiamine (Thi+) and the treatment absent thiamine (Thi-). For each environmental condition, three biological replicates, A, B, and C, were grown and their cellular RNA sequenced via next-gen sequencing techniques. 

The aim of this bioinformatics workflow is to identify *C. albicans* genes that were differentially expressed between Thi+ and Thi- conditions. To do so, we will first clean up the raw sequences, then assemble identified transcripts, and, finally, compare gene expressions between the two treatment conditions across all three replicates. 

## Starting files overview
The two treatment conditions and three biological replicates per condition resulted in six distinct RNA-sequencing samples, each being a paired-end sequencing sample contained in a pair of forward- and reverse-reads FASTQ files. The naming convention to be maintained throughout the workflow is as follows:

The raw RNA-seq files of paired-end sequencing data all derive from the same wildtype (WT) cell line.

| Thiamine+ environment             | Thiamine- environment                 |
| -----------                       | -----------                           |
| **WTA1**: wildtype, replicate A, Thi+ | **WTA2**: wildtype, replicate A, Thi-     |
| **WTB1**: wildtype, replicate B, Thi+ | **WTB2**: wildtype, replicate B, Thi-     |
| **WTC1**: wildtype, replicate C, Thi+ | **WTC2**: wildtype, replicate C, Thi-     |

---

<br></br>

# Menu
- [2021.10.13 - Obtain reads and read counting](#20211013---obtain-reads-and-read-counting)
- [2021.10.19 - Clean reads via Trimmomatic](#20211019---clean-reads-via-trimmomatic)
- [2021.10.21 - Obtain reference genome](#20211021---obtain-reference-genome)
- [2021.10.26 - bowtie2 sequence alignment (mostly for practice)](#20211026---bowtie2-sequence-alignment-mostly-for-practice)
- [2021.10.28 - Tophat alignment](#20211028---tophat-alignment)
- [2021.11.02 - Transfer to GitHub](#20211102---transfer-to-github)
- [2021.11.04 - Infer transcripts using cufflinks](#20211104---infer-transcripts-using-cufflinks)
- [2021.11.09 - Merge transcript annotation files using cuffmerge](#20211109---merge-transcript-annotation-files-using-cuffmerge)
- [2021.11.11 - Identify differentially expressed genes using cuffdiff](#20211111---identify-differentially-expressed-genes-using-cuffdiff)
- [2021.11.16 - Building summary table](##20211116---building-summary-table)
- [2021.11.18 - cummeRbund]()
- [Template](#template)

---

<br></br>

# 2021.10.13 - Obtain reads and read counting
[Back to menu](#menu)

## Objective(s) of this step of analysis.
The goal in this step is to obtain the raw RNA sequencing reads files (focusing on the WTC2 biological replicate) and count the total number of reads in the R1 and R2 fastq files before cleaning. 

## Files involved.
As the class used a divide & conquer strategy on these files, I focused on WTC2's paired-end raw data files.
- WTC2_1.fq.gz
    - forward reads of wildtype, Thi-, replicate C.
- WTC2_2.fq.gz
    - reverse reads of wildtype, Thi-, replicate C.


## Specific commands used in the analysis.
```bash
# obtain the .gz files from Google Bucket
$ gsutil cp gs://gu-biology-dept-class/WTC2*.gz

# unzip the .gz files
$ gunzip *.gz

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
There are **20,407,694** reads in both the forward and the reverse read files for WTC2.

The agreement in read numbers between the forward and reverse read files makes sense because they are obtained from the same paired-end sequencing run.

---

<br></br>

# 2021.10.19 - Clean reads via Trimmomatic
[Back to menu](#menu)

## Objective(s) of this step of analysis.
The goal is to use Trimmomatic to trim the raw reads so as to remove the problematic first ten bases of each read, reduce adapter content, improve reverse read quality, and do an overall trimming via sliding window. 

Review trimmed PE files in FastQC to ensure Trim achieved aims as intended.

## Files involved.
### —Trimmomatic—
### *input files*
- WTC2_1.fq.gz, WTC2_2.fq.gz
    - Raw data files for wildtype C2 treatment. (Yeast grown without thiamine.) (see also previous entry)

### *output files*
- Paired-end and Single-end trimmed sequence files for forward reads (R1)
    - WTC2_1.trPE.fq.gz
    - WTC2_1.trSE.fq.gz
- Paired-end and Single-end trimmed sequence files for reverse reads (R2)
    - WTC2_2.trPE.fq.gz
    - WTC2_2.trSE.fq.gz
- *Note that the trSE versions will not be used in subsequent steps.


## Specific commands used in the analysis.
### Create trim.sbatch
```bash
$ nano trim.sbatch
```
- [**trim.sbatch**](scripts/trim.sbatch)

### Run sbatch Trimmomatic command
```bash
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
```bash
# observe output
$ less z01.trim_WTC2
```
- [**z01.trim_WTC2**](slurm_outputs/z01.trim_WTC2)

```bash
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
```bash
# download (with custom alias)
get_hpc /home/qz108/RNA_seq_workflow/WTC2*trPE*

# open with FastQC & observe.
```

## Results & interpretation.
There are 19,459,631 reads in the trimmed paired-end fastq files.
**95.35% surviving.** 

Per Base Sequence Content greatly improved from before, but still has some visible wobbliness in the first 3 bases. 

High-Duplication sequences do not appear to have been trimmed, which is fine considering this is RNA-seq data.

---

<br></br>

# 2021.10.21 - Obtain reference genome
[Back to menu](#menu)
## Objective(s) of this step of analysis.
Goal is to download refseq C. albicans genome assembly from Entrez genome into a separate folder.

Species: [Candida albicans SC5314 (budding yeasts)](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=237561&lvl=3&lin=f&keep=1&srchmode=1&unlock)
Link at Entrez Genome: [GCA_000182965.3](https://www.ncbi.nlm.nih.gov/assembly/GCA_000182965.3)

## Files involved.
## Specific commands used in the analysis.
```bash
$ mkdir refseq_GCF_000182965.3
$ cd refseq_GCF_000182965.3
$ wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz
```

## Results & interpretation.
File obtained:
GCF_000182965.3_ASM18296v3_genomic.fna.gz

---

<br></br>

# 2021.10.26 - bowtie2 sequence alignment (mostly for practice)
[Back to menu](#menu)

## Objective(s) of this step of analysis.
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
```bash
# Create & navigate to directory
$ mkdir 1021-26_bowtie2_alignment
$ cd 1021-26_bowtie2_alignment
```
### —bowtie2-build—
```bash
# This is to build an index of the reference genome, which is then used by bowtie2 for alignment.

# Enter a compute node
$ srun --pty bash

# Do index run on reference genome
$ bowtie2-build GCF_000182965.3_ASM18296v3_genomic.fna.gz Calbicans
```
### —bowtie2—
```bash
# Copy over the WTC2 trPE fastq files from the outputs of the Trimmomatic run.
$ cp ../1019_Trimmomatic/WTC2*trPE.fq .

# Create sbatch file for bowtie2

$ nano bt2.sbatch
```
- [**bt2.sbatch**](scripts/bt2.sbatch)

```bash
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

---

<br></br>

# 2021.10.28 - Tophat alignment
[Back to menu](#menu)

## Objective(s) of this step of analysis.
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
```bash
$ gunzip GCF_000182965.3_ASM18296v3_genomic.gff.gz

$ nano tophat_align.sbatch
```
- [**tophat_align.sbatch**](scripts/tophat_align.sbatch)

```bash
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
- [**align_summary.txt**](summary_outputs/tophat_align_summary.txt)

## Results & interpretation.
Tophat spliced alignment resulted in 90.3% concordant pair alignment rate. 
Reads aligned to the reference genome were collected in the accepted_hits.bam file.

[See class's alignment results.](https://docs.google.com/spreadsheets/d/1GZOHIUhfYA523kkOdoF9APQEkRsberddvEsEZOkIIvA/edit)

---

<br></br>

# 2021.11.02 - Transfer to GitHub
[Back to menu](#menu)

## Objective(s) of this step of analysis.
Transfer notes to GitHub.

## Files involved.
Repository
- lil-qorgi:379-rnaseq

READ<span>ME.md</span>
- Notes document for RNAseq workflow.

## Specific commands used in the analysis.
```bash
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

---

<br></br>

# 2021.11.04 - Infer transcripts using cufflinks
[Back to menu](#menu)

## Objective(s) of this step of analysis.
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
```bash
# First, pull all required files together into the same folder for easier usage.

# Then write the slurm script for cufflinks
$ nano cufflinks.sbatch
```
- [**cufflinks.sbatch**](scripts/cufflinks.sbatch)

```bash
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

```bash
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

---

<br></br>

# 2021.11.09 - Merge transcript annotation files using cuffmerge
[Back to menu](#menu)

## Objective(s) of this step of analysis.
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
        - this is the annotation file that merges the transcripts.gtf file from the previous (cufflinks) step obtained for all biological replicates and experimental conditions: WTA1, A2, B1, B2, C1, C2. 

## Specific commands used in the analysis.
```bash
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
```bash
# Create an sbatch file for cuffmerge.
$ nano cuffmerge.sbatch
```
- [**cuffmerge.sbatch**](scripts/cuffmerge.sbatch)

```bash
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

```bash
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

---


<br></br>

# 2021.11.11 - Identify differentially expressed genes using cuffdiff
[Back to menu](#menu)

## Objective(s) of this step of analysis.
The goal is to identify genes that are significantly differentially expressed between the control (Thi+) and treatment (Thi-) groups ...

i.e., between WTX1 (Thi+) vs. WTX2 (Thi-), where X is one of {A, B, C} 

... using cuffdiff.

## Files involved.
### —cuffdiff—
### *input files*
Place all input files for cuffdiff into a new folder named "cuffdiff_input"
- cuffdiff_input/
    - **merged.gtf**
        - the feature annotation file from the previous step that merged all biological replicates (A1 through C2) via cuffmerge.
    - &ast;*The following files are the  accepted hits for all replicates under Thi+ and Thi- treatments*:
        - **WTA1_accepted_hits.bam**
        - **WTB1_accepted_hits.bam**
        - **WTC1_accepted_hits.bam**
        - **WTA2_accepted_hits.bam**
        - **WTB2_accepted_hits.bam**
        - **WTC2_accepted_hits.bam**
        - these accepted hits .bam files were obtained from tophat alignment over Trimmomatic-cleaned RNAseq data for each of the samples. In other words, we repeated the previous analysis steps from the first step (2021.10.13) through to the tophat alignment step (2021.10.28) for each raw sequence reads file (excluding the bowtie2 step, which was for practice). The links below reference those steps.
            1. [Obtain the reference genome](#20211021---obtain-reference-genome) for *Candida albicans*. It only needs to be obtained once.
            1. [Obtain the raw RNA-seq paired-end fastq files](#20211013---obtain-reads-and-read-counting) for the biological replicate (e.g. WTB1)
            1. Clean the reads via [Trimmomatic]((#20211019---clean-reads-via-trimmomatic))
            1. Run the paired-end reads against the reference genome using the [Tophat spliced alignment program](#20211028---tophat-alignment) to obtain accepted_hits.bam for the replicate. 
                - After obtaining the file, prepend accepted_hits.bam with the replicate name (e.g. WTB1_accepted_hits.bam), then copy the file into cuffdiff_input/

### *output files*
- cuffdiff_output/
    - **gene_exp.diff**
    - ... (many other cuffdiff results)
- [**z01.cuffdiff**](slurm_outputs/z01.cuffdiff)

## Specific commands used in the analysis.
```bash
$ mkdir cuffdiff_input
# copy merged.gtf and all of the <sample>.bam files into cuffdiff_input/

# edit cuffdiff.sbatch
$ nano cuffdiff.sbatch
```
[cuffdiff.sbatch](scripts/cuffdiff.sbatch)

```bash
# run sbatch script. Note that cuffdiff_inputs/ is parallel to cuffdiff.sbatch
$ sbatch cuffdiff.sbatch

# wait for job to finish...

# check job status
$ sacct -j <jobid> --format=jobid,jobname,account,state,elapsed

# it took 01:20:09 (~80 minutes) for the job to finish.

# observe the slurm output 
$ less z01.cuffdiff
```
[**z01.cuffdiff**](slurm_outputs/z01.cuffdiff)
```bash
# first change into the output folder, then check the output files
$ cd cuffdiff_output

$ ls -1
bias_params.info
cds.count_tracking
cds.diff
cds_exp.diff
cds.fpkm_tracking
cds.read_group_tracking
gene_exp.diff
genes.count_tracking
genes.fpkm_tracking
genes.read_group_tracking
isoform_exp.diff
isoforms.count_tracking
isoforms.fpkm_tracking
isoforms.read_group_tracking
promoters.diff
read_groups.info
run.info
splicing.diff
tss_group_exp.diff
tss_groups.count_tracking
tss_groups.fpkm_tracking
tss_groups.read_group_tracking
var_model.info

# check whether the gene_exp.diff file was properly generated 
$ less gene_exp.diff

# the file appears to have 14 tab-delimited columns and contents. This file appears to be properly generated.
```

## Results & interpretation.
Cuffdiff differential expression analysis output appeared successful. We can now move on to building the RNAseq summary table of differentially genes. 


---

<br></br>

# 2021.11.16 - Building summary table
[Back to menu](#menu)

## Objective(s) of this step of analysis.
Using the results stored in gene_exp.diff from the previous cuffdiff run as well as Uniprot and Entrez Nucleotide websites (for protein functional annotation), we will build a summary table of the genes that were (statistically) significantly differentially expressed in the absence of thiamine when compared with in the presence of thiamine.

The method is to link the genes along with their FPKMs, log2 FC, and q-value from cuffdiff with NCBI protein IDs and known descriptions of their functions from Uniprot and Entrez Nucleotide.

## Files involved.
- **gene_exp.diff**
    - value_1
        - FPKM for the gene in Thiamine+ samples
    - value_2
        - FPKM for the gene in Thiamine- samples 



## Specific commands used in the analysis.
Columns we'll need in our summary table:
- TUXEDO SUITE ID
- NCBI ID
- GENE LABEL (acronym)
- FPKM_THI+
- FPKM_THI-
- LOG2_FC (THI- relative THI+)
- Q-VALUE (FDR CORRECTED)
- NOTES INTERPRETATION (from uniport and entrez nucleotide)

To begin, we'll copy merged.gtf and gene_exp.diff into a separate folder. 

In this case, I named that folder 1116_cuffdiff_results_interpretation/

We'll first obtain the columns we want from gene_exp.diff.
```bash
# The columns that we need from gene_exp.diff are: 1, 3, 8, 9, 10, 13
$ grep "yes" gene_exp.diff | cut -f1,3,8-10,13 > partial_summary.txt
```

We'll next obtain the NCBI IDs that correspond to the TUXEDO SUITE IDs (XLOC IDs), which are in column 1 of partial_summary.txt, by grepping the XLOC IDs against merged.gtf, which contains the NCBI IDs.
```bash
# We first extract only the XLOC ID column from the "yes," i.e., differentially expressed, genes. This is similar to the previous command, but includes only the XLOC ID column.
$ grep "yes" gene_exp.diff | cut -f1 > signif_xlocIDs

# We grep these significantly regulated xlocIDs (signif_xlocIDs) against merged.gtf.
$ grep -wFf signif_xlocIDs merged.gtf > results_summary

# From results_summary, we extract the information representing xlocIDs, gene names, and NCBI IDs. The gene names column is used for double-checking whether it lines up with the gene names from partial_summary.txt
$ cut -f9 results_summary|cut -d " " -f2,8,10 > results_summary2.txt
```

Now that we have both the file containing gene_exp.diff's relevant columns and file containing their corresponding NCBI IDs, we'll first download them to the local computer.
```bash
# The following commands are run on local, inside my RNA-seq local folder.

# Use gcloud compute scp to copy partial_summary.txt and results_summary2.txt to a new subfolder.

# Rename the files to more clearly represent what they are used for. 
$ mv partial_summary.txt partial_from_gene_exp.txt
$ mv results_summary2.txt partial_from_merged.txt
```

Create 


export them both into Excel, sort both, then check if they align. If so, we will then join the tables by introducing .


## Results & interpretation.

---



<br></br>

# —Template—

<copy-from-here>

<br></br>

# yyyy.mm.dd - Title
[Back to menu](#menu)

## Objective(s) of this step of analysis.

## Files involved.

## Specific commands used in the analysis.
```bash
```

## Results & interpretation.

---

<copy-to-here>
