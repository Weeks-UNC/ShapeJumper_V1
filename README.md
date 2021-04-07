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
*bash ShapeJumper.sh referenceSequence.fasta crosslinkRead1.fastq crosslinkread2.fastq controlRead1.fastq controlRead2.fastq*
**Unpaired Reads**  
*bash ShapeJumper.sh referenceSequence.fasta crosslinkUnPairedReads.fastq controlUnPairedReads.fastq*

### Example files  
Provided example files from a SHAPE-JuMP experiment on RNase P Catalytic domain, a small RNA with available sructure on pdb: 3DHS.  
RNA was crosslinked with TBIA and IA was used as a mono-adduct control.  

**Reference Fasta File**
*RNaseP_WithStructureCassette.fa*  
Contains the reference 268 nucleotides of DNA sequence for the RNase P Catalytic domain.  
Also included are the 5' and 3' structure cassettes. Most _in vitro_ studies of small RNAs use transcripts with structure cassettes to aid in library prep.  
The 5' structure casette is 14 nucleotides long, the 3' cassette is 43 nucleotides.  

See https://doi.org/10.1021/ja043822v for more in depth explanation.  

**Paired End Reads**
*TBIA-RNaseP_S1_L001_R1_001.fastq  
TBIA-RNaseP_S1_L001_R2_001.fastq  
IA-RNaseP_S2_L001_R1_001.fastq  
IA-RNaseP_S2_L001_R2_001.fastq*  

**Unpaired Reads**
*TBIA-RNaseP.extendedFrags.fastq  
IA-RNaseP.extendedFrags.fastq*  

## Output Description
