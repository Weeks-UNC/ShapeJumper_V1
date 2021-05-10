# ShapeJumper_V1.0
Identifies and processes deletions from Illumina sequencing reads.  
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
\>ReferenceSequenceName1  
AGCTAGA  
\>ReferenceSequenceName2  
AAAAGGGTTGACGAGA  

**Fastq file** - text files obtained from sequencer. Fastq files can be paired or unpaired. Fastq file for both crosslink and control sample required.
Fastq file format:
@Read Name
Sequence
+
Quality (PHRED scores)
*See https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf for an explanation of Quality scores generated by Illumina sequencing  
See https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm to decode ASCII symbols of PHRED scores*  


### Execution Instructions:
Copy script ShapeJumper.sh to your current working directory. The directory you run this script in is where output will be generated.  
**Paired End Reads**  
`bash ShapeJumper.sh referenceSequence.fasta crosslinkRead1.fastq crosslinkread2.fastq controlRead1.fastq controlRead2.fastq`  
**Unpaired Reads**  
`bash ShapeJumper.sh referenceSequence.fasta crosslinkUnPairedReads.fastq controlUnPairedReads.fastq`  

Input file names: The fastq files are assumed to have Illumina generated file names for paired end reads and FLASH generated filenames for Unpaired.  
Paired end read file names should contain a sample name followed by `\_S##\_L001\_R#\_001.fastq` where ## is an ilumina generated sample number and # is the read number, 1 or 2.  
Unpaired reads should contain the sample name followed by `.extendedFrags.fastq`.  
Failure to follow naming conventions will result in truncated and misnamed output files.

### Example files  
Provided example files from a SHAPE-JuMP experiment on RNase P Catalytic domain, a small RNA with an available structure in the pdb: 3DHS.  
RNA was crosslinked with TBIA and IA was used as a mono-adduct control.  
*Note:* All read files have been compressed. To extract files use command `tar -xzf fileName.tar.gz`

**Reference Fasta File**  
`RNaseP_WithStructureCassette.fa`  
Contains the reference 268 nucleotides of DNA sequence for the RNase P catalytic domain.  
Also included are the 5' and 3' structure cassettes. Most _in vitro_ studies of small RNAs use transcripts with structure cassettes to aid in library prep.  
The 5' structure casette is 14 nucleotides long, the 3' cassette is 43 nucleotides.  

See https://doi.org/10.1021/ja043822v for a more in depth explanation of structure cassettes.  

**Paired End Reads**  
`TBIA-RNaseP_S1_L001_R1_001.fastq  
TBIA-RNaseP_S1_L001_R2_001.fastq`  
`IA-RNaseP_S2_L001_R1_001.fastq  
IA-RNaseP_S2_L001_R2_001.fastq`  

25,000 paired end reads from a SHAPE-JuMP experiment on RNase P catalytic domain. TBIA samples were crosslinked, IA is the mono-adduct control.  

**Unpaired Reads**  
`TBIA-RNaseP.extendedFrags.fastq`  
`IA-RNaseP.extendedFrags.fastq`  

25,000 reads from the same SHAPE-JuMP experiment on RNase P, but these have already been merged with FLASH to create paired end reads.

## Output Description  

**Deletion Text File**  
The final output from succesful execution of ShapeJumper will be stored in a text file ending in `_Merged_ProcessedDeletions.txt` and starting with the name of the crosslinked sample.  
Deletions in this file have been fully processed: Rates have been normalized and subtracted by control deletion rates. Ambiguous deletions and and those with inserts in the deletion longer than 10 nucleotides have been removed. Exact edge matching at deletion sites has been enforced. The 5' deletion start sites have been shifted 2 nucleotides downstream. If selected, numbering has been adjusted for to account for structure cassettes.

EXAMPLE OUTPUT:  
Total Reads Aligned:11450       Total Deletions:728.0
rnasep  123     184     0.0017623086
rnasep  81      94      0.0013324450

The first line is a header, denoting total reads aligned in the crosslinked sample. Total Deletions are the raw count of deletions longer than 10 nucleotides observed in the crosslinked sample, regardless of downstream filtering by ambiguous deletions or so on.  
Subsequent lines follow the same 4 column format:  
- Column 1 = Reference name from fasta file matching alignment. Samples with multiple reference sequences may contain multiple names.
- Column 2 = Deletion start site. Numbering is relative to reference fasta sequence.
- Column 3 = Deletion stop site.
- Column 4 = Deletion rate frequency. This value is normalized by read depth.

## Intermediate Output Files  
A folder will be generated during execution, ShapeJumperIntermediateFiles. In it are contained all the files generated during ShapeJumper execution. Files are generated for both crosslink and control inputs. All files are text files.  

File Name Guide:
- `.extendedFrags.fastq` = Pair mate merged reads post FLASH
- `.notCombined_1.fastq` = Reads that were unable to be merged by FLASH. A \_2 file is generated for read 2.
- `.sam` = Alignments generated by BWA-MEM from extendedFrags fastq files.
- `_NoFlash.sam` = Alignments of reads that did not merge by FLASH.   
- `_Merged.sam` = The alignment files from both FLASH merged and unmerged reads concatenated together into one file.  
- `_deletions.txt` = Set of deletions identified in Merged.sam file. Each deletion found is recorded with nucleotide coordinates of the deletion and the frequency it is found.
- `_normalizedDels.txt` = Deletions with counts normalized by read depth. If the option to renumber deletion sites to account for structure cassettes was selected, that is implementd here.  
- `_Subtracted.txt` = Set of crosslinked sample deletions with normalized rates subtracted by normalized rates of control sample.  
- `_DelReadNames.txt` = Stores every read name of sequences containing one or more identified deletions. Can be useful for follow up/in-depth analysis.
- `.fa.amb, .fa.ann, .fa.bwt, .fa.pac, .fa.sa` = index files generated by BWA-MEM


*Note:* Depending on options selected, some files will not be present.
