# ShapeJumper_V1.0
Identifies and processes deletions from illumina sequencing reads.  
Tom Christy 2021
twchrist@email.unc.edu

### Requirements:  
Python 2.7  including modules numpy, matplotlib, argparse, re, 
ShapeMapper-V2, available at https://github.com/Weeks-UNC/shapemapper2
Flash, available at https://sourceforge.net/projects/flashpage/  
BWA, available at http://bio-bwa.sourceforge.net/
All programs should be installed and added to the path

### Installation:
clone this package to your home directory with this command:
git clone https://github.com/Weeks-UNC/ShapeJumper_V1.0.git

Within the now downloaded ShapeJumper_V1.0 folder is the script ShapeJumper.sh
Copy this script to your working directory.

### Required inputs:
**Fasta file** - a text file containing one or more reference DNA sequences of the RNAs probed in SHAPE-JuMP experiment.
Example Format:
>gene1  
AGCTAGA  
>gene2  
AAAAGGGTTGACGAGA  

**Fastq file** - text files obtained from sequencer. Fastq files can be paired or unpaired. Fastq file for both crosslink and control sample required.
Fastq file format:
Read Name
Sequence
+
Quality (PHRED scores)

### Execution Instructions:
Copy script ShapeJumper.sh to your current working directory. The directory you run this script in is where output will be generated.
**Paired End Reads**
bash ShapeJumper.sh referenceSequence.fasta crosslinkRead1.fastq crosslinkread2.fastq controlRead1.fastq controlRead2.fastq"
**Unpaired Reads**
bash ShapeJumper.sh referenceSequence.fasta crosslinkUnPairedReads.fastq controlUnPairedReads.fastq

### Example files
**Reference Fasta File**

**Paired End Reads**

**Unpaired Reads**
## Output Description
