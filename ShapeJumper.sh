#!/bin/bash
#ShapeJumper.sh

#June 29, 2020
#takes in a reference fasta file and reads. Reads may be paired end or unpaired.
#Process all reads through ShapeJumper pipeline and outputs a list of all deletions
#identified and their rates.

set -e # exit on first error (if any)

#handle arguments, including optional arguments
script="ShapeJumper"
#Declare the number of mandatory args
margs=3

#examples of usage
function example {
    echo -e "example usage - unpaired reads: \nbash ShapeJumper.sh -f referenceSequence.fasta -t1 crosslinkUnPairedReads.fastq -c1 controlUnPairedReads.fastq"
	echo "OR"
    echo -e "example usage - paired end reads: \nbash ShapeJumper.sh -f referenceSequence.fasta -t1 crosslinkRead1.fastq -t2 crosslinkread2.fastq -c1 controlRead1.fastq -c2 controlRead2.fastq"	
}

function usage {
    echo -e "usage: bash $script -f FASTA -t1 TEST -c1 CONTROL [OPTIONALARGUMENTS]\n"
}

function help {
  usage
    echo -e "Mandatory Arguments:"
    echo -e "  -f, --fasta  REFERENCEFASTA.FASTA  \n\"Fasta file containing sequences of all target RNAs.\"\n"
    echo -e "  -t1, --test1  TESTREADS.FASTQ  \n\"Fastq file containing reads from test/crosslinked sample.\"\n"
    echo -e "  -c1, --control1  CONTROLREADS.FASTQ  \n\"Fastq file containing reads from control/monoAdduct sample.\"\n"
    echo -e "Optional Arguments:"
    echo -e "  -t2, --test2   TESTREADS2.FASTQ   \n\"Fastq file containing paired reads 2 from test/crosslinked sample.\"\n"
    echo -e "  -c2, --control2   CONTROLREADS2.FASTQ  \n\"Fastq file containing paired reads 2 from control/monoAdduct sample.\"\n"
    echo -e "  -sc, --structureCassette \n\"This option renumbers the output deletions to remove structure cassettes.\"\n"
    echo -e "  -s, --shift5Prime SHIFTINTEGER \n\"Enter a custom shift to the 5\' end rather than the default 2.\"\n"
    echo -e "  -a, --keepAmbiguousDeletions \n\"Keeps ambiguous deletions in the final output rather than removing them by default.\"\n"
    echo -e "  -n, --normalizeByBothEnds \n\"Normalize deletions by read depths at the 5\' and 3\' ends rather than just 3\'.\"\n"
    echo -e "  -h,  --help             \n\"Prints this help\"\n"
  example
}

# Ensures that the number of passed args are at least equals
# to the declared number of mandatory args.
# It also handles the special case of the -h or --help arg.
function margs_precheck {
	if [ $2 ] && [ $1 -lt $margs ]; then
		if [ $2 == "--help" ] || [ $2 == "-h" ]; then
			help
			exit
		else
	    	usage
			example
	    	exit 1 # error
		fi
	fi
}

# Ensures that all the mandatory args are not empty
function margs_check {
	if [ $# -lt $margs ]; then
	    usage
	  	example
	    exit 1 # error
	fi
}

# Main method
#make sure mandatory arguments were input
margs_precheck $# $1

#set empty variables for mandatory arguments
fasta=
testFastq1=
controlFastq1=
#set default values for optional arguments
testFastq2=
controlFastq2=
sc=
nt5Shift=2
ambig="--NoAmbiguousDeletions"
norm=

# Args while-loop, parses command line arguments and stores them in above variables
while [ "$1" != "" ];
do
   case $1 in
   -f  | --fasta )  shift
                          fasta=$1
                		  ;;
   -t1  | --test1 )  shift
   						  testFastq1=$1
			              ;;
   -c1  | --control1  )  shift
   							controlFastq1=$1
                          ;;
   -t2  | --test2  )  shift
                          testFastq2=$1
                          ;;
   -c2  | --control2  )  shift
                          controlFastq2=$1
                          ;;
   -s  | --shift5Prime  )  shift
                          nt5Shift=$1
                          ;;
   -sc  | --structureCassette  )  sc="--sc"
                          ;;
   -a  | --keepAmbiguousDeletions  )  ambig=""
                          ;;  
   -n  | --normalizeByBothEnds  )  norm="--bothEnds"
                          ;;                                            
   -h   | --help )        help
                          exit
                          ;;
   *)                     
                          echo "$script: illegal option $1"
                          usage
						  example
						  exit 1 # error
                          ;;
    esac
    shift
done


#process reads if they are unparied
if [ ! -e "$testFastq2" ] && [ ! -e "$controlFastq2" ]; then
	#trim reads
	mkdir -p ShapeJumperIntermediateFiles
	echo "Trimming Reads:"
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $testFastq1 --out ShapeJumperIntermediateFiles/$testFastq1
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $controlFastq1 --out ShapeJumperIntermediateFiles/$controlFastq1
	#copy the fasta file into the intermediate directory
	cp $fasta ShapeJumperIntermediateFiles
	#all analysis will take place in this intermediate directory
	cd ShapeJumperIntermediateFiles
	
	#Run BWA-MEM to align reads
	echo "Aligning Reads:"
	bwa index $fasta
	testSam=${testFastq1:0:${#testFastq1} - 6}".sam"
	bwa mem $fasta $testFastq1 -O 2 -B 2 -k 10 -T 15 >${testSam}
	controlSam=${controlFastq1:0:${#controlFastq1} - 6}".sam"
	bwa mem $fasta $controlFastq1 -O 2 -B 2 -k 10 -T 15 >${controlSam}

else
	#process reads if they are paired
	#trim reads
	mkdir -p ShapeJumperIntermediateFiles
	echo "Trimming Reads:"
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $testFastq1 --out ShapeJumperIntermediateFiles/$testFastq1
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $testFastq2 --out ShapeJumperIntermediateFiles/$testFastq2
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $controlFastq1 --out ShapeJumperIntermediateFiles/$controlFastq1
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $controlFastq2 --out ShapeJumperIntermediateFiles/$controlFastq2

	#copy the fasta file into the intermediate directory
	cp $fasta ShapeJumperIntermediateFiles
	#all analysis will take place in this intermediate directory
	cd ShapeJumperIntermediateFiles
	
	#Flash the trimmed reads together
	echo "Mate Pair Merging Reads:"
	flash $testFastq1 $testFastq2 -o ${testFastq1:0:${#testFastq1} - 21}
	flash $controlFastq1 $controlFastq2 -o ${controlFastq1:0:${#controlFastq1} - 21}
	#flash generates unneccessary histograms, remove them
	rm *.hist
	rm *.histogram
	testFlashed=${testFastq1:0:${#testFastq1} - 21}".extendedFrags.fastq"
	controlFlashed=${controlFastq1:0:${#controlFastq1} - 21}".extendedFrags.fastq"
	
	#Run BWA MEM on Flashed reads
	echo "Aligning Reads:"
	testSam=${testFlashed:0:${#testFlashed} - 20}".sam"
	controlSam=${controlFlashed:0:${#controlFlashed} - 20}".sam"
	bwa index $fasta
	bwa mem $fasta $testFlashed -O 2 -B 2 -k 10 -T 15 >$testSam
	bwa mem $fasta $controlFlashed -O 2 -B 2 -k 10 -T 15 >$controlSam
	
	#Run BWA MEM on the reads that could not be Flashed together
	testNotCombined1=${testFastq1:0:${#testFastq1} - 21}".notCombined_1.fastq"
	testNotCombined2=${testFastq2:0:${#testFastq2} - 21}".notCombined_2.fastq"
	testNoFlashSam=${testFastq1:0:${#testFastq1} - 21}"_NoFlash.sam"
	controlNotCombined1=${controlFastq1:0:${#controlFastq1} - 21}".notCombined_1.fastq"
	controlNotCombined2=${controlFastq2:0:${#controlFastq2} - 21}".notCombined_2.fastq"
	controlNoFlashSam=${controlFastq1:0:${#controlFastq1} - 21}"_NoFlash.sam"
	bwa mem $fasta $testNotCombined1 $testNotCombined2 -O 2 -B 2 -k 10 -T 15 >$testNoFlashSam
	bwa mem $fasta $controlNotCombined1 $controlNotCombined2 -O 2 -B 2 -k 10 -T 15 >$controlNoFlashSam
	
	#Merge the flashed and not flashed sam files together
	name=${testSam:0:${#testSam} - 4}
	header=$testSam
	files="$testSam $testNoFlashSam"
	testMergedSam=${name}_Merged.sam
	echo $testSam
	echo $testNoFlashSam
	echo $testMergedSam
	(grep ^@ $header; for f in $files; do grep -v ^@ $f; done) > $testMergedSam
	name=${controlSam:0:${#controlSam} - 4}
	header=$controlSam
	files="$controlSam $controlNoFlashSam"
	controlMergedSam=${name}_Merged.sam
	(grep ^@ $header; for f in $files; do grep -v ^@ $f; done) > $controlMergedSam
	
	#rename merged sam variable for continuity
	testSam=$testMergedSam
	controlSam=$controlMergedSam
	
fi
#find deletions in both samples now that reads have been aligned
echo "Finding Deletions:"
python ~/ShapeJumper_V1.0/ShapeJumperScripts/identifyDeletions.py $testSam $fasta -rb $ambig -e -mi 10 
python ~/ShapeJumper_V1.0/ShapeJumperScripts/identifyDeletions.py $controlSam $fasta -rb $ambig -e -mi 10 

#normalize deletions
echo "Normalizing Deletions:"
testDeletions=${testSam:0:${#testSam} - 4}"*_deletions.txt"
controlDeletions=${controlSam:0:${#controlSam} - 4}"*_deletions.txt"
python ~/ShapeJumper_V1.0/ShapeJumperScripts/normalizeDeletionRates.py $testDeletions $testSam $fasta $sc $norm
python ~/ShapeJumper_V1.0/ShapeJumperScripts/normalizeDeletionRates.py $controlDeletions $controlSam $fasta $sc $norm

echo "Subtracting normalized control deletions from normalized crosslinked deletions"
testNormDeletions=${testSam:0:${#testSam} - 4}"*_*ormalizedDels.txt"
controlNormDeletions=${controlSam:0:${#controlSam} - 4}"*_*ormalizedDels.txt"
python ~/ShapeJumper_V1.0/ShapeJumperScripts/subtractDeletionRates.py $testNormDeletions $controlNormDeletions --s

echo "Shift 5\' nucleotide position by $nt5Shift"
subtractedDeletions=${testNormDeletions:0:${#testNormDeletions} - 4}"_Subtracted.txt"
python ~/ShapeJumper_V1.0/ShapeJumperScripts/shiftDeletionSites.py $subtractedDeletions --shift5nt $nt5Shift

#rename shifted deletion file
if [ $nt5Shift -eq 0 ]; then
	shiftedDeletions=${subtractedDeletions}
else
	shiftedDeletions=${subtractedDeletions:0:${#subtractedDeletions} - 4}"_shift*.txt"
fi
if [ ! -e "$testFastq2" ] && [ ! -e "$controlFastq2" ]; then
	mv $shiftedDeletions ../${testSam:0:${#testSam} - 4}"_ProcessedDeletions.txt"
else
	mv $shiftedDeletions ../${testSam:0:${#testSam} - 11}"_ProcessedDeletions.txt"
fi
cd ..
echo "Done!"
