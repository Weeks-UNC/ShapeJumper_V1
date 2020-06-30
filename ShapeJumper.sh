#!bin/bash
#ShapeJumper.sh

#June 29, 2020
#takes in a reference fasta file and reads. Reads may be paired end or unpaired.
#Process all reads through ShapeJumper pipeline and outputs a list of all deletions
#identified and their rates.

#check for argument count
if [ $# -eq 3 ]; then
	fasta=$1
	testFastq1=$2
	controlFastq1=$3
elif [ $# -eq 5 ]; then
	fasta=$1
	testFastq1=$2
	testFastq2=$3
	controlFastq1=$4
	controlFastq2=$5

else
	echo "Wrong number of arguments"
	echo "USAGE: bash ShapeJumper.sh referenceSequence.fasta crosslinkRead1.fastq crosslinkread2.fastq controlRead1.fastq controlRead2.fastq"
	echo "OR"
	echo "bash ShapeJumper.sh referenceSequence.fasta crosslinkUnPairedReads.fastq controlUnPairedReads.fastq"
	exit 1
fi
	
#check the files to make sure they exist
if [ ! -e "$fasta" ]
then
	echo "ERROR: In $0, the input fasta file $fasta does not exist"
	exit 1
fi

if [ ! -e "$testFastq1" ]
then
	echo "ERROR: In $0, the input fastq file $testFastq1 does not exist"
	exit 1
fi

if [ ! -e "$controlFastq1" ]
then
	echo "ERROR: In $0, the input fastq file $controlFastq1 does not exist"
	exit 1
fi


#process reads if they are unparied
if [ $# -eq 3 ]; then
	#trim reads
	mkdir ShapeJumperIntermediateFiles
	echo "Trimming Reads:"
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $testFastq1 --out ShapeJumperIntermediateFiles/$testFastq1
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $controlFastq1 --out ShapeJumperIntermediateFiles/$controlFastq1
	#copy the fasta file into the intermediate directory
	cp $fasta ShapeJumperIntermediateFiles
	#all analysis will take place in this intermeditate directory
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
	mkdir ShapeJumperIntermediateFiles
	echo "Trimming Reads:"
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $testFastq1 --out ShapeJumperIntermediateFiles/$testFastq1
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $testFastq2 --out ShapeJumperIntermediateFiles/$testFastq2
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $controlFastq1 --out ShapeJumperIntermediateFiles/$controlFastq1
	shapemapper_read_trimmer --min_phred 20 --min_length 25 --window_size 5 --in $controlFastq2 --out ShapeJumperIntermediateFiles/$controlFastq2

	#move original reads into their own directory, move trimmed reads into current directory
	#copy the fasta file into the intermediate directory
	cp $fasta ShapeJumperIntermediateFiles
	#all analysis will take place in this intermeditate directory
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
python ~/ShapeJumper_V1.0/ShapeJumperScripts/identifyDeletions.py $testSam $fasta -rb -a -e -mi 10 -i 1 -m
python ~/ShapeJumper_V1.0/ShapeJumperScripts/identifyDeletions.py $controlSam $fasta -rb -a -e -mi 10 -i 1 -m

#normalize deletions
echo "Normalizing Deletions:"
testDeletions=${testSam:0:${#testSam} - 4}"*_deletions.txt"
controlDeletions=${controlSam:0:${#controlSam} - 4}"*_deletions.txt"
python ~/ShapeJumper_V1.0/ShapeJumperScripts/normalizeDeletionRates.py $testDeletions $testSam $fasta --sc
python ~/ShapeJumper_V1.0/ShapeJumperScripts/normalizeDeletionRates.py $controlDeletions $controlSam $fasta --sc

echo "Subtracting normalized control deletions from normalized crosslinked deletions"
testNormDeletions=${testSam:0:${#testSam} - 4}"*_normalizedDels.txt"
controlNormDeletions=${controlSam:0:${#controlSam} - 4}"*_normalizedDels.txt"
python ~/ShapeJumper_V1.0/ShapeJumperScripts/subtractDeletionRates.py $testNormDeletions $controlNormDeletions --s

echo "Shift 5\' nucleotide position by 2"
subtractedDeletions=${testNormDeletions:0:${#testNormDeletions} - 4}"_Subtracted.txt"
python ~/ShapeJumper_V1.0/ShapeJumperScripts/shiftDeletionSites.py $subtractedDeletions --shift5nt 2

#rename shifted deletion file
shiftedDeletions=${subtractedDeletions:0:${#subtractedDeletions} - 4}"_shift*.txt"
mv $shiftedDeletions ../${testSam:0:${#testSam} - 4}"_ProcessedDeletions.txt"
cd ..
echo "Done!"
