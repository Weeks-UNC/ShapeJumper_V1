#normalizeDels.oy

#October 24, 2018

#This script is an improvement over the previous normalize deletions. It is able to handle paired end reads without falsely double counting nucleotides that have overlap
import sys
import numpy as np
import re
import argparse
import matplotlib
matplotlib.use("PDF")
from matplotlib import pyplot as plt

def parseArgs():
	prs = argparse.ArgumentParser()
	prs.add_argument("deletions", type=str, help="Unique Deletion Locations")
	prs.add_argument("sam", type=str, help="SAM file from BWAMEM Alignment")
	prs.add_argument("fasta", type=str, help="fasta file of target sequence")
	prs.add_argument('--sc',action='store_true',default=False,help='Corrects deletion locations for structure cassettes and removes all deletions that occur within structure cassettes')
	prs.add_argument('--bothEnds',action='store_true',default=False,help='Determine the median depth of the 5nts at BOTH ends of the deletion and divide the deletion Count by square rooted product of these medians.')
	prs.add_argument('--HiC',action='store_true',default=False,help='HiC Normalization: delcount is divided by the square root of the products of the read depths of the start and stop NTs of the deletion')
	prs.add_argument('--plotReadDepth',action="store_true",help="Plot the read depth per nucleotide and save the figure. Limits are set to make visualizing fluctuations easier.")
	o = prs.parse_args()
	return o

# returns array with start and stop sites of the alignment relative to the reference sequence
def pullSeqEnds(array):
	start = int(array[3].strip())
	end = start - 1
	cigar = array[5].strip()
	s = re.findall('(\d+[IDMN])', cigar)
	for pair in s:
		# pull out number and type
		p = re.compile('(\d+)([MIDN])')
		m = p.search(pair)
		num = int(m.group(1))
		type = m.group(2)
		# run through cigar options
		if (type == 'D' or type == 'M' or type == 'N'):
			end = end + num
	return [start, end]
# returns array with start and stop sites of the alignment relative to the read sequence
def pullReadEnds(array):
	cigar = array[5].strip()

	start = 1
	#find any clips at the beginning of the cigar, use to move the start
	if re.match('\d+[SH]', cigar):
		m = re.match('(\d+)[SH]', cigar)
		num = int(m.group(1))
		start += num
	end = start - 1
	#run through the cigar string and add together all pertinent operation types to discern where the read stops aligning
	s = re.findall('(\d+[IM])', cigar)
	for pair in s:
		# pull out number and type
		p = re.compile('(\d+)([MI])')
		m = p.search(pair)
		num = int(m.group(1))
		type = m.group(2)
		# run through cigar options
		if (type == 'I' or type == 'M'):
			end = end + num
	return [start,end]


# takes in a target tuple and a query tuple and determines if any of the query is within the target.
def inBetween(target, query):
	inB = False
	if (target[0] <= query[0] <= target[1]):
		inB = True
	if (target[0] <= query[1] <= target[1]):
		inB = True
	return inB


#This function takes in a sam flag an reports if the flag denotes that the alignment if
#of a reverse read strand
def isReverseStrand(samFlag):
	if int(samFlag[-5]) == 1:
		return True
	else:
		return False

#Takes in a sam flag and reports if the alignment is a main alignment
def isPrimaryAlignment(samFlag):
	if int(samFlag[-12]) == 0 and int(samFlag[-9]) == 0:
		return True
	else:
		return False

def isSecondRead(samFlag):
	if int(samFlag[-8]) == 1:
		return True
	else:
		return False

#make sure that the alignments aren't significantly longer than the length of the the read
def checkReadLen(seqCoords, readLength):
	alignedLen = 0
	for i in range(len(seqCoords)):
		a = seqCoords[i]
		alignedLen += a[1][1] - a[1][0]
	if alignedLen > readLength + 10:
		print "HEY youre alignments are using more than the whole read"
		print readName
		return False
	else:
		return True

#takes in a start site, cigar string, and ntList and adds the depth from that cigar string to the list
def addDepth(start, cigar, siteList):
	s = re.findall('(\d+[IDM])',cigar)
	#end = start -1
	for pair in s:
		#pull out number and type
		p = re.compile('(\d+)([MID])')
		m = p.search(pair)
		num = int(m.group(1))
		type = m.group(2)
		#run through cigar options
		if(type == 'M' or type =='D'):
			end = start + num -1
			#now add to depth
			siteList[start:end+1] = np.logical_or(siteList[start:end+1],1)
			start = end + 1
			#end = start - 1
	return siteList


if __name__ == '__main__':

	args = parseArgs()


	#load in fasta sequence and get sequence length
	inF2 = open(args.fasta,'r')
	refGeneName = inF2.readline().strip()[1:]
	refSeqLen = len(inF2.readline().strip())
	print refGeneName
	print "Reference sequence length is "+str(refSeqLen)
	inF2.close()

	# import sam file
	inputFileName = args.sam
	inF = open(inputFileName, 'r')
	lines = inF.readlines()
	inF.close()

	#create an empty array that will hold how often each nucleotide in the reference sequence has been aligned to
	ntList = np.zeros(refSeqLen + 1)
	# note I am treating this array as if it was not 0 indexed

	#store the first index each new read is encountered in a list
	readIndices = list()
	prevRead = ""
	for i in range(len(lines)):
		line = lines[i]
		#make sure we're not looking at the header
		if line[0] != "@":
			cols = line.strip().split()
			if cols[0] != prevRead:
				readIndices.append(i)
				prevRead = cols[0]
	#add in the last index so we can grab the last set of alignments too
	readIndices.append(len(lines))

	# use this list of indices to process all the alignments for each read in batches.
	for i in range(len(readIndices) - 1):
		# instantiate variable for each read
		# grab the range of indices of all alignment for the current read
		startIndex = readIndices[i]
		stopIndex = readIndices[i + 1]
		anyAlignmentsFound = False
		currentReadNtList = np.zeros(refSeqLen + 1)
		for index in range(startIndex, stopIndex):
			alignment = lines[index].strip().split()
			# convert flag into 16 bit binary
			flag = format(int(alignment[1]), '016b')
			cigarString = alignment[5]
			# if we're on the first alignment, determine the gene we're aligning to
			if index == startIndex:
				if alignment[2] == '*':
					continue
				else:
					geneName = alignment[2]
					if geneName != refGeneName:
						break
					readName = alignment[0]

			# skip this alignment if it does not align to a gene
			if alignment[2] == '*' or alignment[5] == '*':
				continue
			else:
				anyAlignmentsFound = True

			# if alignment doesn't match gene name of first alignment, skip the alignment
			if geneName != alignment[2]:
				continue
			# if the alignment is a primary alignment, note the direction. I'm assuming that the first alignment will be a primary alignment
			if isPrimaryAlignment(flag):
				reverseStrand = isReverseStrand(flag)
				readLen = len(alignment[9])
			# if the alignment is secondary and disagrees with the primary alignment direction, skip the alignment
			if reverseStrand != isReverseStrand(flag):
				continue

			#now that we've made absolutely sure the current alignment aligns to the right gene and is not a secondary alignment reusing part of the read from the first alignment
			#process the cigar string and add 1 depth to any nucleotide covered by this cigar string that wasn't covered before by the same read
			seqStart = int(alignment[3])
			currentReadNtList = addDepth(seqStart, cigarString, currentReadNtList)

		#now that every alignment has been checked, if any of them aligned to the target sequence add the depth from
		#that alignment to the total depth of the whole sequence
		if anyAlignmentsFound:
			ntList = ntList + currentReadNtList

	#now that all alignments depths have been totalled, plot the read depth, if requested
	if args.plotReadDepth:
		name = args.sam[:-4] + refGeneName
		plt.figure(name+" Read Depth", figsize=(15, 10))
		plt.plot(ntList[1:])
		plt.title(args.sam[:-4]+" "+refGeneName+" Read Depth")
		plt.ylabel('Read Depth')
		plt.xlabel('Nucleotide')
		plt.tight_layout()
		plt.xlim(1,refSeqLen)
		plt.ylim(min(ntList[1:])- (max(ntList[1:])/10) ,max(ntList[1:]) + (max(ntList[1:])/10) )
		#print plt.ylim()
		#plt.show()
		plt.savefig(name+"_Read_Depth.pdf")
	#print ntList[10:20]
	#load in the deletions, and print them out after dividing their rate by the depth.
	#Also take into account the type of normalization selected and structure cassette filtering

	delF = open(args.deletions, 'r')
	delLines = delF.readlines()
	delF.close()

	newDelList = list()
	# loop through file
	for i in range(1, len(delLines)):
		# load in line and pull out end and del count
		line = delLines[i].strip().split()
		name = line[0]
		if (name == refGeneName):
			start = int(line[1])
			end = int(line[2])
			count = float(line[3])
			# handle HiC Case
			if (args.HiC):
				# grap the 5' end depth
				startDepth = ntList[start]
				# grap the 3' end depth
				endDepth = ntList[end]
				# modify the count by the square of the product
				count = count / np.sqrt(startDepth * endDepth)
			elif args.bothEnds:
				# pull the depths of the 5nts at the start and end of the deletion
				startSubList = ntList[start - 4:start + 1]
				endSubList = ntList[end:end + 5]
				# get the median depth of these sets of 5 nt depths
				startMed = float(np.median(startSubList))

				endMed = float(np.median(endSubList))
				# modify the count by the square root of the product of the median depths
				count = count / np.sqrt(startMed * endMed)
			else:
				# handle the standard read depth normalization case.
				# grab 5 nts after end of deletion and get median read depth
				endSubList = ntList[end:end + 5]
				endMed = float(np.median(endSubList))
				# normalize read depth and put that back into deletion line
				count = count / endMed
			line[3] = count
			# correct location for structure cassettes if aplicable
			if (args.sc):
				line[1] = int(line[1]) - 14
				line[2] = int(line[2]) - 14
				if (line[1] < 1 or line[2] > refSeqLen-14-43):
					# do not include deletions that occur in the 5' SC
					continue
			newDelList.append(line)
	# sort the new Del List
	#print "\nDeletions before sorting"
	#print newDelList
	newDelList = sorted(newDelList, key=lambda x: x[3], reverse=True)
	#print "\nDeletions after sorting"
	#print newDelList
	# add in header line
	newDelList.insert(0, delLines[0])
	# print new normalized deletion list
	if (args.HiC):
		outName = args.deletions[:-14] + "_HiCNormalizedDels.txt"
	elif args.bothEnds:
		outName = args.deletions[:-14] + "_bothEndsNormalizedDels.txt"
	else:
		outName = args.deletions[:-14] + "_normalizedDels.txt"
	outF = open(outName, 'w')
	outF.write(str(newDelList[0]))
	for i in range(1, len(newDelList)):
		line = newDelList[i]
		outF.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str(line[2]) + '\t' + str(line[3]) + "\n")
	outF.close()
