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

Input file names: The fastq files are assumed to have Illumina generated file names for paired end reads and FLASH generated filenames for Unpaired.  
Paired end read file names should contain a sample name followed by *\_S##\_L001\_R#\_001.fastq* where ## is an ilumina generated sample number and # is the read number, 1 or 2.  
Unpaired reads should contain a sample name followed by *.extendedFrags.fastq*.  
Failure to follow naming conventions will result in truncated and misnamed output files.

### Example files  
Provided example files from a SHAPE-JuMP experiment on RNase P Catalytic domain, a small RNA with available sructure on pdb: 3DHS.  
RNA was crosslinked with TBIA and IA was used as a mono-adduct control.  

**Reference Fasta File**  
*RNaseP_WithStructureCassette.fa*  
Contains the reference 268 nucleotides of DNA sequence for the RNase P catalytic domain.  
Also included are the 5' and 3' structure cassettes. Most _in vitro_ studies of small RNAs use transcripts with structure cassettes to aid in library prep.  
The 5' structure casette is 14 nucleotides long, the 3' cassette is 43 nucleotides.  

See https://doi.org/10.1021/ja043822v for more in depth explanation.  

**Paired End Reads**  
*TBIA-RNaseP_S1_L001_R1_001.fastq  
TBIA-RNaseP_S1_L001_R2_001.fastq  
IA-RNaseP_S2_L001_R1_001.fastq  
IA-RNaseP_S2_L001_R2_001.fastq*  

100,000 paired end reads from a SHAPE-JuMP experiment on RNase P catalytic domain. TBIA samples were crosslinked, IA is the mono-adduct control.  

**Unpaired Reads**  
*TBIA-RNaseP.extendedFrags.fastq  
IA-RNaseP.extendedFrags.fastq*  

Again 100,000 reads from a SHAPE-JuMP experiment on RNase P, but these have already been merged with FLASH to create paired end reads.

## Output Description  

**Deletion Text File**  
The final output from succesful execution of ShapeJumper will be stored in a text file ending in *_Merged_ProcessedDeletions.txt* and starting with the name of the crosslinked sample.  
EXAMPLE OUTPUT:  
Total Reads Aligned:437994      Total Deletions:27806.0  
rnasep  123     184     0.0015252557  
rnasep  113     184     0.0012315027  

The first line is a header, denoting total reads aligned in the crosslinked sample. Total Deletions are the raw count of deletions longer than 10 nucleotides observed in the crosslinked sample.  
Subsequent lines follow the same 4 column format:  
- Column 1 = Reference name from fasta file matching alignment. Samples with multiple reference sequences may contain multiple names.
- Column 2 = 
- Column 3 =
- Column 4 = 
