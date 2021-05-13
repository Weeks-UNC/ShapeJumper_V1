#FindDeletions.py

#Tom Christy
#October 11 2018

#This script is designed to take in a sam file and report the deletions within. It also takes in a fasta file that can have
#multiple sequences, to report which deletion belongs to which sequence.

#The reason this script is different from the rest of my find deletions script is two-fold
# 1) It will work on both paired end and flashed reads
# 2) It will now find all the deletions in a read > 10 nt rather than just the longest.
# Bonus: There are now options for checking against fake deletion reads, outputting intramolecular merged reads and isolating deletions with insertions.

import sys
import argparse
import re
import numpy as np


def parseArgs():
	prs = argparse.ArgumentParser()
	prs.add_argument("SamFile", type=str, help="Input SAM file from a BWAMEM alignment")
	prs.add_argument("fasta",type=str,help="Input Fasta file with reference sequence.")
	prs.add_argument("-f","--FakeDeletions", action='store_true',help="Output read name with every deletion to check found deletion vs actual sites stated in the read name.")
	prs.add_argument("-i","--IndelOutput",type=int, default=0,help="Input the minimum insert size to for an indel. Print out a specific deletions file for all the deletions containing indels")
	prs.add_argument("-mut","--mutationOutput", nargs=2, default=[-1,-1], help="Input primer lengths. When printing out deletions by each read, add a column for the mutation rate of the read. Primer length is used to determine regions to ignore.")
	prs.add_argument("-m","--MergedOutput",action='store_true',help="Print out any reads that appear contain intraMolecular Merges.")
	prs.add_argument("-w","--writeNotFound",action='store_true',help="Write to ouput files all reads that did not align and all reads without deletions.")
	prs.add_argument("-mi","--MaxInsertLimit",type=int, default = -1, help="Input a number, only deletions with an insert less than or equal to this number will be output.")
	prs.add_argument("-li","--MinInsertLimit",type=int, default = -1, help="Input a number, only deletions with an insert greater than or equal to this number will be output.")
	prs.add_argument("-a","--NoAmbiguousDeletions",action='store_true',help="Remove any deletions where one of the sites has ambiguity as to where it should be.")
	prs.add_argument("-o","--OnlyAmbiguousDeletions",action='store_true',help="Output only deletions where one or both sites have ambiguity.")
	prs.add_argument("-e","--ExactMatchAtDelSite",action='store_true',help="Requires that the 3 nucleotides at the start/stop of the deletion exactly match the reference. This is to mitigate error from inserts.")
	prs.add_argument("-n","--InsertIdentity",default="",type=str,help="Input a nucleotide. Will return deletions that contain inserts contating only that nucleotide")
	prs.add_argument("-rb","--removeBackwardDels",action="store_true",help="If selected, any deletions found where the start site comes after the stop site will be removed. This will remove all deletions that detect 5\' to 3'\ circularization as well as some odd non template activity.")
	o = prs.parse_args()
	return o

#takes in the start and stop coordinates of two reads to make sure the alignments aren't re-using nucleotides
def inside(a1nt1,a1nt2,a2nt1,a2nt2):
	problem = False
	if a1nt1 >= a2nt1 and a1nt1 <= a2nt2 and a1nt2 >= a2nt1 and a1nt2 <= a2nt2:
		problem = True
	if a2nt1 >= a1nt1 and a2nt1 <= a1nt2 and a2nt2 >= a1nt1 and a2nt2 <= a1nt2:
		problem = True
	return problem

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

#takes in a reference sequence index and determines if the nt does not occur within a found deletion
def validSite(index, delSiteList):
	valid = True
	for d in delSiteList:
		#index is not valid if it falls within the deletion
		if index > d[0] and index < d[1]:
			valid = False
	return valid
#takes in a read and outputs the number of mutations and length of reference sequence analyzed
def parseCigarString(cigarString, read, refSeq, startIndex, ntSiteList, delSiteList):
	# define regular expression patterns for parsing CIGAR string
	mismatchCount = 0.0
	smallDelCount = 0.0
	length = 0.0
	item = re.compile(r"[0-9]+[M|I|D|S|H]")
	number = re.compile(r"([0-9]+)")
	splitCigarString = item.findall(cigarString.strip())
	cigarList = [number.split(x)[1:] for x in splitCigarString]
	# print "CIGAR: %s\nsplit CIGAR: %s"%(cigarString,str(cigarList))
	readCopy = read
	refIndex = startIndex - 1  # begin at leftmost matching position
	for i in range(len(cigarList)):
		region = cigarList[i]
		regionLength = int(region[0])
		regionType = region[1]
		# alignment match (could be match or mismatch, use raw read to resolve)
		if regionType == "M":

			matchRegion = read[:regionLength]
			read = read[regionLength:]
			for regionIndex in range(len(matchRegion)):
				nuc = matchRegion[regionIndex]
				#make sure that the currently analyzed nt has not already been analyzed
				if ntSiteList[refIndex] == 0 and validSite(refIndex, delSiteList):
					#update ntSiteList now that this index has been considered
					ntSiteList[refIndex] = 1
					length += 1
					if nuc != refSeq[refIndex]:
						# nuc mismatch found
						mismatchCount += 1

				refIndex += 1
		# insertion
		elif regionType == "I":
			read = read[regionLength:]
			#okay this might be a bit hacky but I'm going to use the reference sequence site to determine if this nucleotide has already been counted
			if ntSiteList[refIndex+1] == 0 and validSite(refIndex+1, delSiteList) and regionLength <= 20:
				mismatchCount += regionLength
		# deletion
		elif regionType == "D":
			for deletionIndex in range(regionLength):
				if (regionLength < 6):
					if ntSiteList[refIndex] == 0 and validSite(refIndex, delSiteList):
						ntSiteList[refIndex] = 1
						mismatchCount += 1
						smallDelCount += 1
						length += 1
				refIndex += 1

		# padding
		elif (regionType == "S"):
			read = read[regionLength:]
	# nothing is required to handle H events

	# return the del rates for the 5' and then 3' end of the read.
	# if the deletion is split between two alignments, the first entry will be 0
	#return [smallDelCount, length]
	return [mismatchCount, length]

#takes in two sequences of the same length, determines if they are the same. Used to make
#sure the ends of the deletions have reads that perfectly match the reference
def sequencesMatch(seq1, seq2):
	same = True
	try:
		for n in range(len(seq1)):
			if seq1[n] != seq2[n]:
				same = False
	except:
		print "Out of bounds at "+readName
		print "seq1 is "+seq1
		print "seq2 is "+seq2
		sys.exit()
	return same

#takes in a deletion and some other pertinent info and returns the deletion such that the
#deletion start and stop sites exactly match
def makeDelEndsExactMatch(deletion, readNT, readSeq, refSeq):
	global badMatchCount
	ref1NT = deletion[0]
	ref2NT = deletion[1]
	insert = deletion[2]
	insertNts = deletion[3]
	shiftCounter = 0
	# handle the 5' deletion start site
	exactMatchAtStartSiteFound = False
	#print "Examing 5\' side"
	while exactMatchAtStartSiteFound == False:
		# grab the 5' sequence from the reference and the read
		refRegion = refSeq[ref1NT - 3:ref1NT]
		readRegion = readSeq[readNT - 3:readNT]
		#print "refNt is "+str(ref1NT)
		#print refRegion
		#print "readNT is " + str(readNT)
		#print readRegion
		#make sure that the read region is not empty and longer than 3 nts
		if not readRegion or len(readRegion) < 3:
			shiftCounter = 12
			allMatch = False
			break
		# determine if the two segments match
		allMatch = sequencesMatch(refRegion, readRegion)
		# if they don't match move the read and referenct nts back one as well as the length of the deletion
		if allMatch == False:
			ref1NT = ref1NT - 1
			readNT = readNT - 1
			insert = insert + 1
			insertNts = readSeq[readNT] + insertNts
			shiftCounter+=1
			if shiftCounter > 10:
				break
		# if the sequences do match, break out of this loop, but after updating the deletion coordinates
		exactMatchAtStartSiteFound = allMatch

	# handle 3' side of the deletion now
	exactMatchAtStopSiteFound = False
	readNT = readNT + insert
	#print "Examing 3\' side"
	while exactMatchAtStopSiteFound == False:
		# grab the 3' sequence from the reference and the read
		refRegion = refSeq[ref2NT-1:ref2NT + 2]
		readRegion = readSeq[readNT:readNT + 3]
		#print "refNt is " + str(ref2NT)
		#print refRegion
		#print "readNT is " + str(readNT)
		#print readRegion
		# make sure that the read region is not empty and longer than 3 nts
		if not readRegion or len(readRegion) < 3:
			shiftCounter = 12
			allMatch = False
			break
		# determine if the two segments match
		allMatch = sequencesMatch(refRegion, readRegion)
		# if they don't match move the read and reference nt up one as well as the length of the deletion
		if allMatch == False:
			insertNts = insertNts + readSeq[readNT]
			ref2NT = ref2NT + 1
			readNT = readNT + 1
			insert = insert + 1
			shiftCounter+=1
			if shiftCounter > 10:
				break
		# if the sequences do match, break out of this loop, but after updating the deletion coordinates
		exactMatchAtStopSiteFound = allMatch

		deletion = [ref1NT, ref2NT, insert, insertNts]
	if shiftCounter > 10:
		print readName+ " had to shift "+str(shiftCounter)+" nucleotides and still did not find exact match, skipping read "+readName
		badMatchCount += 1
		deletion = [-1,-1,-1, '']
	return deletion
#make sure that the alignments aren't significantly longer than the length of the the read
def checkReadLen(seqCoords, readLength):
	global overlappingAlignmentCount
	alignedLen = 0
	for i in range(len(seqCoords)):
		a = seqCoords[i]
		alignedLen += a[1][1] - a[1][0]
	if alignedLen > readLength + 10:
		print "Alignments are using more than the whole read at read "+readName+" it is being skipped."
		overlappingAlignmentCount += 1
		return False
	else:
		return True
#this function determines if an array containing a deletion matches to deletion that already exists in the current list
#does not takes insert size into acount, also if both ends of the deletion are within 2 nts of an existing deletion, that is considered a match
def freshDeletion(newDeletion, currentDelList):
	fresh = True
	newDelStart = newDeletion[0]
	newDelStop = newDeletion[1]
	for currentDel in currentDelList:
		currentDelStart = currentDel[0]
		currentDelStop = currentDel[1]
		if abs(newDelStart - currentDelStart) <= 5 or abs(newDelStop - currentDelStop) <= 5:
			fresh = False
	return fresh


#run through the nucleotides in the insert and determine if they all match the given nt
def checkInsertIdentity(insertNucleotides, nt):
	allMatchNt = True
	if insertNucleotides:
		for i in range(len(insertNucleotides)):
			if insertNucleotides[i] != nt:
				allMatchNt = False
	else:
		allMatchNt = False
	return allMatchNt

#take in an existing deletion list and an alignment and add all within alignment deletions longer than 10 nt to the deletion list
def findInAlignmentDels(alignmentCols, delList, referenceSequence, readSequence):
	#print referenceSequence
	#print readSequence
	insertSize = 0

	start = int(alignmentCols[3].strip())
	currentNt = start - 1
	currentReadNt = 0
	cigar = alignmentCols[5].strip()
	cigarIndex = 0
	#first search for two long deletions with an insertion in the middle cause these can look like two separate deletions.
	s = re.findall('(\d+[SHIDMN])', cigar)
	skipCount = 0
	for pair in s:
		#so this is if I found a deletion that was spread across multiple pairs in the cigar string
		#the other pairs that were used to make the deletion are skipped
		#print "Processing "+pair
		insertNts = ''
		if skipCount> 0:
			skipCount = skipCount -1
			continue
		# pull out number and type
		p = re.compile('(\d+)([SHMIDN])')
		m = p.search(pair)
		num = int(m.group(1))
		type = m.group(2)
		# run through cigar options
		if type == 'S' or type == 'H' or type == 'I':
			cigarIndex = cigarIndex + len(m.group(1)) + 1
			if type != 'H':
				currentReadNt += num
		if (type == 'M' or type == 'N'):
			currentNt = currentNt + num
			currentReadNt += num
			cigarIndex = cigarIndex + len(m.group(1)) + 1
		if (type == 'D'):
			#look for a case where an insert has been split up within a bigger deletion
			m3 = re.match('(\d+)D(\d+)M(\d+)D(\d+)M(\d+)D', cigar[cigarIndex:])
			validTrip = False
			if m3:
				del1 = int(m3.group(1))
				match1 = int(m3.group(2))
				del2 = int(m3.group(3))
				match2 = int(m3.group(4))
				del3 = int(m3.group(5))
				insertNts = readSequence[currentReadNt: currentReadNt + match1 + match2]
				newDel = [currentNt, currentNt + del1 + match1 + del2 + match2+ del3 + 1, match1 + match2, insertNts]
				if del1 >= 3 and del2 >= 3 and del3 >= 3 and del1 + del2 + del3 > 10 and match1 + match2 <= 10 and freshDeletion(newDel, delList):
					validTrip = True
					skipCount = 4

					# if requested, make sure nucleotides at deletion start and stop sites are exact matches
					if args.ExactMatchAtDelSite:
						newDel = makeDelEndsExactMatch(newDel, currentReadNt, readSequence, referenceSequence)
					if args.NoAmbiguousDeletions or args.OnlyAmbiguousDeletions:
						fiveRefNuc = referenceSequence[newDel[0] - 1]
						threeRefNuc = referenceSequence[newDel[1] - 1]
						fiveSeqNuc = referenceSequence[newDel[0]]
						threeSeqNuc = referenceSequence[newDel[1] - 2]
						if (fiveSeqNuc != threeRefNuc and fiveRefNuc != threeSeqNuc) and args.NoAmbiguousDeletions:
							delList.append(newDel)
						if (fiveSeqNuc == threeRefNuc or fiveRefNuc == threeSeqNuc) and args.OnlyAmbiguousDeletions:
							delList.append(newDel)
					else:
						delList.append(newDel)
					currentNt = newDel[1] -1
					currentReadNt = currentReadNt + match1 + match2
					cigarIndex = cigarIndex + len(m3.group(1)) + len(m3.group(2)) + len(m3.group(3)) + len(m3.group(4)) + len(m3.group(5)) + 5
					continue

			#see if this deletion fits pattern of two deletions interrupted by an insert
			m2 = re.match('(\d+)D(\d+)M(\d+)D',cigar[cigarIndex:])
			validDub = False
			if m2 and validTrip == False:
				del1 = int(m2.group(1))
				match = int(m2.group(2))
				del2 = int(m2.group(3))
				insertNts = readSequence[currentReadNt: currentReadNt + match]
				newDel = [currentNt, currentNt + del1 + match + del2 + 1, match, insertNts]
				if del1 >= 3 and del2 >= 3 and del1 + del2 >10 and match <= 10 and freshDeletion(newDel, delList):
					validDub = True
					skipCount = 2
					# if requested, make sure nucleotides at deletion start and stop sites are exact matches
					if args.ExactMatchAtDelSite:
						newDel = makeDelEndsExactMatch(newDel, currentReadNt, readSequence, referenceSequence)

					if args.NoAmbiguousDeletions or args.OnlyAmbiguousDeletions:
						fiveRefNuc = referenceSequence[newDel[0] - 1]
						threeRefNuc = referenceSequence[newDel[1]-1]
						fiveSeqNuc = referenceSequence[newDel[0]]
						threeSeqNuc = referenceSequence[newDel[1]-2]
						if (fiveSeqNuc != threeRefNuc and fiveRefNuc != threeSeqNuc) and args.NoAmbiguousDeletions:
							delList.append(newDel)
						if (fiveSeqNuc == threeRefNuc or fiveRefNuc == threeSeqNuc) and args.OnlyAmbiguousDeletions:
							delList.append(newDel)
					else:
						delList.append(newDel)
					currentNt = newDel[1] -1
					currentReadNt = currentReadNt + match
					cigarIndex = cigarIndex + len(m2.group(1)) + len(m2.group(2)) + len(m2.group(3)) + 3
					continue

			#store deletions larger than 10 that haven't already been found
			#deletions are stored as the last nucleotide to align and the first new nucleotide to align, hence the +1
			newDel = [currentNt,currentNt+num+1,insertSize, '']
			if num > 10 and freshDeletion(newDel,delList) and validTrip == False and validDub == False:
				#if requested, make sure nucleotides at deletion start and stop sites are exact matches
				if args.ExactMatchAtDelSite:
					newDel = makeDelEndsExactMatch(newDel, currentReadNt, readSequence, referenceSequence)

				if args.NoAmbiguousDeletions or args.OnlyAmbiguousDeletions:
					fiveRefNuc = referenceSequence[newDel[0] - 1]
					threeRefNuc = referenceSequence[newDel[1] - 1]
					fiveSeqNuc = referenceSequence[newDel[0]]
					threeSeqNuc = referenceSequence[newDel[1] - 2]
					if (fiveSeqNuc != threeRefNuc and fiveRefNuc != threeSeqNuc) and args.NoAmbiguousDeletions:
						delList.append(newDel)
					if (fiveSeqNuc == threeRefNuc or fiveRefNuc == threeSeqNuc) and args.OnlyAmbiguousDeletions:
						delList.append(newDel)
				else:
					delList.append(newDel)
			currentNt = currentNt + num
			cigarIndex = cigarIndex + len(m.group(1)) + 1
		#print "Current NT is "+str(currentNt)
		#print "Current Read NT is "+str(currentReadNt)
	return delList

#takes in a list of the sequence and read coordinates of all the alignments. Use these to find any deletions from
#the multi alignments.
def findMultiAlignDels(seqCoords, delList, readLength, referenceSequence, readSequence):
	global weirdDoubleSequences
	#set the amount of read overlap allowed
	if args.NoAmbiguousDeletions or args.OnlyAmbiguousDeletions:
		overlap = 0
	else:
		overlap = -2

	#if there is only one alignment for the read, exit
	if len(seqCoords) == 1:
		return delList
	else:
		#sort the seq coordinates so that start points of the reads are ordered from lowest to highest
		#this should mean that sequences coordinates will also be ordered this way, except in the case of
		#circularization detection where the 3' end of the sequence is at the 5' end of the read
		seqCoords = sorted(seqCoords, key=lambda x: x[:][1])
		insertNts = ''
		#compare each set of coordinates to its next neighbor
		for i in range(len(seqCoords)-1):
			a1 = seqCoords[i]
			a2 = seqCoords[i+1]
			#print "a1 is "
			#print a1
			#print "a2 is "
			#print a2
			#print seqCoords
			#grab the start/stop sites of both alignments
			a1SeqStart = a1[0][0]
			a1SeqStop = a1[0][1]
			a1ReadStart = a1[1][0]
			a1ReadStop = a1[1][1]
			a2SeqStart = a2[0][0]
			a2SeqStop = a2[0][1]
			a2ReadStart = a2[1][0]
			a2ReadStop = a2[1][1]
			#yell real loud if the alignemnts are using the same portion of the read
			if inside(a1ReadStart,a1ReadStop,a2ReadStart,a2ReadStop):
				print readName + " ALIGNMENTS ARE REUSING THE SAME NTS FROM THE READ, worth a more in depth look"
			#run through all previous sets of seq coordinates and determine if they match the current reference seq coordinates. this is to prevent any issues from the same sequence being duplicated
			for a0 in seqCoords[:i+1]:
				a0SeqStart = a0[0][0]
				a0SeqStop = a0[0][1]
				if inside(a0SeqStart, a0SeqStop, a2SeqStart, a2SeqStop):
					print readName + " ALIGNMENTS ARE REUSING THE SAME NTS FROM THE READ FROM THE SEQUENCE at: "
					weirdDoubleSequences += 1
					print alignment
					print "skipping this read"
					#if different parts of the read are aligning to the same sequence, something weird is happening, so set the deletion to -1. This will cause the deletion to be skipped.
					return [[-1, -1, -1,""]]
			#make sure that the reads aren't overlapping
			if(a2ReadStart - a1ReadStop > overlap):
				insertSize = a2ReadStart - a1ReadStop - 1
				if insertSize < 0:
					insertSize = 0
				else:
					insertNts = readSequence[a1ReadStop:a2ReadStart - 1]
				# only add the deletions if the sequence distance is > 10 or handle circularization case where read starts 3' of the end of the read (in terms of sequence space) because of jumping from one end to the other
				newDel = [a1SeqStop, a2SeqStart, insertSize, insertNts]
				if(a2SeqStart - a1SeqStop > 10) and freshDeletion(newDel,delList)  or (a2SeqStart - a1SeqStop < -10) and (a2SeqStop - a1SeqStart < -10) and freshDeletion(newDel, delList):
					newDel = [a1SeqStop, a2SeqStart, insertSize, insertNts]
					if args.ExactMatchAtDelSite:
						# handle the 5' deletion start site first, then do the 3' stop site
						exactMatchAtStartSiteFound = False
						while exactMatchAtStartSiteFound == False:
							# grab the 5' sequence from the reference and the read
							refRegion = referenceSequence[a1SeqStop - 3:a1SeqStop]
							readRegion = readSequence[a1ReadStop - 3:a1ReadStop]
							# determine if the two segments match
							allMatch = sequencesMatch(refRegion, readRegion)
							# if they don't match move the read and referenct nts back one as well as the length of the deletion
							if allMatch == False:
								a1SeqStop = a1SeqStop - 1
								a1ReadStop = a1ReadStop - 1
								insertSize = insertSize + 1
							# if the sequences do match, break out of this loop, but after updating the deletion coordinates
							exactMatchAtStartSiteFound = allMatch
							newDel = [a1SeqStop, a2SeqStart, insertSize]

						exactMatchAtStopSiteFound = False
						while exactMatchAtStopSiteFound == False:
							# grab the 5' sequence from the reference and the read
							refRegion = referenceSequence[a2SeqStart:a2SeqStart+3]
							readRegion = readSequence[a2ReadStart:a2ReadStart+3]
							# determine if the two segments match
							allMatch = sequencesMatch(refRegion, readRegion)
							# if they don't match move the read and referenct nts back one as well as the length of the deletion
							if allMatch == False:
								a2SeqStart = a2SeqStart + 1
								a2ReadStart = a2ReadStart + 1
								insertSize = insertSize + 1
							# if the sequences do match, break out of this loop, but after updating the deletion coordinates
							exactMatchAtStopSiteFound = allMatch
							newDel = [a1SeqStop, a2SeqStart, insertSize]
					insertNts = readSequence[a1ReadStop:a2ReadStart-1]
					newDel = [a1SeqStop, a2SeqStart, insertSize, insertNts]

					if args.NoAmbiguousDeletions or args.OnlyAmbiguousDeletions:
						fiveRefNuc = referenceSequence[a1SeqStop-1]
						threeRefNuc = referenceSequence[a2SeqStart-1]
						try:
							fiveSeqNuc = referenceSequence[a1SeqStop]
						except:
							#print "Hey I'm trying to check an ambiguous deletion start site at the end of the sequence"
							#when the start site of the deletion is at the end of the reference sequence, the deletion cannot be ambiguous in this direction.
							#So i'm assuming the deletion is real (caused by circularization) and using the try statement to bypass the normal error that would be caused by this check.
							fiveSeqNuc = ""
							#sys.exit()
						threeSeqNuc = referenceSequence[a2SeqStart-2]
						#print newDel
						#print "fiveSeqNuc " + fiveSeqNuc + " at "+str(a1SeqStop+1)
						#print "threeRefNuc " + threeRefNuc + " at "+str(a2SeqStart)
						#print "fiveRefNuc " + fiveRefNuc + " at "+str(a1SeqStop)
						#print "threeSeqNuc " + threeSeqNuc + " at "+str(a2SeqStart-1)
						
						if (fiveSeqNuc != threeRefNuc and fiveRefNuc != threeSeqNuc) and args.NoAmbiguousDeletions:
							delList.append(newDel)
						if (fiveSeqNuc == threeRefNuc or fiveRefNuc == threeSeqNuc) and args.OnlyAmbiguousDeletions:
							delList.append(newDel)
					else:
						delList.append(newDel)
				#handle the merged read case of the same RNA species interacting with itself
				if a2SeqStart - a1SeqStop < -10 and (a2SeqStop - a1SeqStart > 0) and checkReadLen(seqCoords, readLength):
					print readName + " might be a merged read."
					if args.MergedOutput:
						outF6.write(readName+"\t"+str(a1SeqStart)+"\t"+str(a1SeqStop)+"\t"+str(a2SeqStart)+"\t"+str(a2SeqStop)+"\n")
					#I don't trust any deletions found in a merged read
					return [[-1,-1]]

		return delList



if __name__ == '__main__':

	args = parseArgs()
	
	#redirect print statements to log file
	old_stdout = sys.stdout
	log_file = open(args.SamFile[:-4]+"_DeletionReport.log","w")
	sys.stdout = log_file
	
	# import sam file
	inputFileName = args.SamFile
	inF = open(args.SamFile, 'r')
	lines = inF.readlines()
	inF.close()

	 # load in reference fasta file
	inF2 = open(args.fasta, 'r')
	fastaLines = inF2.readlines()
	inF2.close()

	# load fasta file sequences into a dict indexed by the gene name
	refFasta = {}
	for i in range(len(fastaLines)):
		line = fastaLines[i]
		if line[0] == '>':
			#stores the gene name linked to the sequence. The gene name doesn't have the > symbol or any text after the first space stored
			#because this is how the aligner parses out the gene name from a fasta file.
			refFasta[line[1:].rstrip().split()[0]] = fastaLines[i + 1]

	# keep track of deletions  and alignments not found
	if args.writeNotFound:
		outputFileName = inputFileName[:-4] + "_NoAlign.txt"
		outF2 = open(outputFileName, 'w')
		outF2.write("\#Contains all alignments that did not align to a gene\n")
		outputFileName = inputFileName[:-4] + "_NotFound.txt"
		outF3 = open(outputFileName, 'w')
		outF3.write("\#Contains all alignments of reads where no deletions was found\n")
	outputFileName = inputFileName[:-4] + "_DelReadNames.txt"
	outF4 = open(outputFileName, 'w')
	if args.FakeDeletions:
		outputFileName = inputFileName[:-4] + "_DelPerRead.txt"
		outF5 = open(outputFileName, 'w')
		outF5.write("ReadName\tdelStart\tdelStop\tMutationRate\n")
		delByReadList = []
	if args.MergedOutput:
		outputFileName = inputFileName[:-4] + "_MergedReads.txt"
		outF6 = open(outputFileName, 'w')
		outF6.write("GeneName\tStart1\tStop1\tStart2\tStop2\t\#These are cases where the same rna has overlapping alignments, looking like an interRNA jump on the same sequence\n")
	if args.IndelOutput > 0:
		outputFileName = inputFileName[:-4] + "_Indels.txt"
		outF7 = open(outputFileName, 'w')
		outputFileName = inputFileName[:-4] + "_IndelPerRead.txt"
		outF8 = open(outputFileName, 'w')
		outF8.write("GeneName\tReadName\tDelStart\tDelStop\tInsertSize\tInsertedNucleotides\n")
	# set up counters and data storage
	totalReadsAligned = 0
	totalDeletions = 0.0
	totalIndels = 0.0
	multiDelReads = 0.0
	backwardsDels = 0.0
	weirdDoubleSequences = 0
	badMatchCount = 0
	overlappingAlignmentCount = 0
	deletions = list()
	twoStepDeletions = list()
	delDict = {}
	indelDict = {}
	#keeps track on the actual nucleotides making up the inserts
	insertDict = {}
	insertDict['None'] = 0
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
	refSeq = ""
	# use this list of indices to process all the alignments for each read in batches.
	for i in range(len(readIndices)-1):
		#instantiate variable for each read
		deletions = []
		seqCoordinates = []
		#grab the range of indices of all alignment for the current read
		startIndex = readIndices[i]
		stopIndex = readIndices[i+1]
		anyAlignmentsFound = False
		#store all reads that get analyzed, used later when calculating mutation rates
		pertinentReads = list()
		for index in range(startIndex,stopIndex):
			alignment = lines[index].strip().split()
			#convert flag into 16 bit binary
			flag = format(int(alignment[1]),'016b')

			# if we're on the first alignment, determine the gene we're aligning to
			if index == startIndex:
				if alignment[2] == '*':
					if args.writeNotFound:
						outF2.write(lines[index])
					continue
				else:
					geneName = alignment[2]
					#if requested, make sure that gene matches one of the genes in the fasta, otherwise, don't process it
					if refFasta.has_key(geneName) == False:
						break
					else:
						refSeq = refFasta[geneName]
						refSeqLen = len(refSeq)
					readName = alignment[0]
					totalReadsAligned += 1

			#skip this alignment if it does not align to a gene
			if alignment[2] == '*' or alignment[5] == '*':
				if args.writeNotFound:
					outF2.write(lines[index])
				continue
			else:
				anyAlignmentsFound = True


			#if alignment doesn't match gene name of first alignment, skip the alignment
			if geneName != alignment[2]:
				continue
			pertinentReads.append(alignment)
			#if this alignment is the first alignment of the second read, process all the seq coordinates of the first read
			#to find deletions and then reset the set of seq coordinates to nothing
			if isSecondRead(flag) and isPrimaryAlignment(flag):

				checkReadLen(seqCoordinates, readLen)
				deletions = findMultiAlignDels(seqCoordinates, deletions, readLen, refSeq, readSeq)
				seqCoordinates = []
			#if the alignment is a primary alignment, note the direction. I'm assuming that the first alignment will be a primary alignment
			if isPrimaryAlignment(flag):
				reverseStrand = isReverseStrand(flag)
				readLen = len(alignment[9])
				readSeq = alignment[9]
			#if the alignment is secondary and disagrees with the primary alignment direction, skip the alignment
			if reverseStrand != isReverseStrand(flag):
				continue

			#find any within alignment deletions
			deletions = findInAlignmentDels(alignment, deletions, refSeq, readSeq)
			#record the start and stop sites of the alignment and of the read
			seqCoordinates.append([pullSeqEnds(alignment),pullReadEnds(alignment)])
		#process all possible found deletions if any of the alignments aligned to the reference
		if anyAlignmentsFound and checkReadLen(seqCoordinates, readLen):
			#Use the set of readAlignments find all deletions invovling multiple alignments
			deletions = findMultiAlignDels(seqCoordinates, deletions, readLen, refSeq, readSeq)
			#see if -1 is in the set of deletions, if it is that means that one of the alignments returned a merged read or had some other issue
			#and thus I don't trust any deletions that come out of that alignment so we'll skip the analysis
			#print seqCoordinates
			#print deletions
			if -1 in (x[0] for x in deletions):
				continue

			#write out readname if deletions were found

			if len(deletions) >= 2:
				#make sure that the deletions don't overlap with each other and far enough away from each other
				overlap = False
				for m in range(len(deletions)-1):
					del1 = deletions[m]
					for n in range(m+1,len(deletions)):
						del2 = deletions[n]
						if inBetween(del1, del2):
							overlap = True
						if del2[0] - del1[1] <= 15:
							overlap = True
				if overlap:
					continue
				multiDelReads += 1
				print readName+" has "+str(len(deletions))+" deletions."
			if len(deletions) > 0:
				outF4.write(readName+"\n")
				#detrmine mutation rate, if requested, of reads that have a deletion
				if int(args.mutationOutput[0]) > -1:
					mutsFound = 0
					sequencedNTs = 0
					# list the length of the reference used to keep track and make sure I'm not double counting
					# mutations at regions where the reads overlap
					siteList = np.zeros(refSeqLen + 1)
					#set primer regions to 1
					for nt in range(1,int(args.mutationOutput[0])+1):
						siteList[nt] = 1
					for nt in range(refSeqLen+1	 - int(args.mutationOutput[1]),refSeqLen+1):
						siteList[nt] = 1
					for alignment in pertinentReads:
						[m, l] = parseCigarString(alignment[5], alignment[9], refSeq, int(alignment[3]), siteList, deletions)
						mutsFound += m
						sequencedNTs += l
					#print "Mutations Found = "+str(mutsFound)
					#print "Read Length = "+str(sequencedNTs)
					mutRate = round(float(mutsFound)/float(sequencedNTs),4)
					#add mutation rate to each deletion
					for d in range(len(deletions)):
						currentDel = deletions[d]
						deletions[d] = [currentDel[0],currentDel[1],currentDel[2],mutRate]



			#if no deletions found write out the alignments
			if len(deletions) == 0 and args.writeNotFound:
				for index in range(startIndex, stopIndex):
					outF3.write(lines[index])


			#now that all deletions in the alignment have been found, add them to the deletion dictionary that totals up all deletions at every position
			for deletion in deletions:
				insertLen = deletion[2]
				#if selected, don't add a deletion where the start site occurs after the stop site
				if args.removeBackwardDels:
					if deletion[0] > deletion[1]:
						#print "remove backwards deletions removed deletion"
						#print deletion
						backwardsDels += 1
						continue

				#if selected, determine if inset matches the desired nucleotide identity
				if args.InsertIdentity:
					if checkInsertIdentity(deletion[3],args.InsertIdentity) == False:
						continue
				if (args.MaxInsertLimit == -1 or insertLen <= args.MaxInsertLimit) and (args.MinInsertLimit == -1 or insertLen >= args.MinInsertLimit):
					if args.FakeDeletions:
						if int(args.mutationOutput[0]) > -1:
							outF5.write(alignment[0] + "\t" + str(deletion[0]) + "\t" + str(deletion[1]) + "\t" + str(deletion[3]) + "\n")
						else:
							outF5.write(alignment[0]+"\t"+str(deletion[0])+"\t"+str(deletion[1])+"\n")
						if args.IndelOutput == 0:
							delByReadList.append([readName,geneName,deletion[0],deletion[1]])
						else:
							delByReadList.append([readName, geneName, deletion[0], deletion[1], deletion[2], deletion[3]])
					totalDeletions += 1
					if (delDict.has_key((geneName, deletion[0], deletion[1]))):
						delDict[geneName, deletion[0], deletion[1]] += 1
					else:
						delDict[geneName, deletion[0], deletion[1]] = 1
					#if requested, output deletions that have an insertion above the specified cutoff
					#and store these indel locations for future output
					if args.IndelOutput > 0:
						#only work with inserts larger than user specified cutoff
						#keeps a running tally of how many indels occur at each location and the total insert length found
						if insertLen >= args.IndelOutput:
							totalIndels += 1.0
							# write every indel read to file
							outF8.write(geneName + "\t" + readName + "\t" + str(deletion[0]) + "\t" + str(deletion[1]) + "\t" + str(insertLen) + "\t" + str(deletion[3]) + "\n")
							if (indelDict.has_key((geneName, deletion[0], deletion[1]))):
								indelDict[geneName, deletion[0], deletion[1]][0] += 1
								indelDict[geneName, deletion[0], deletion[1]][0] += insertLen
							else:
								indelDict[geneName, deletion[0], deletion[1]] = [1, insertLen]

							#add the actual nucleotides of the insert to the dictionary that keeps track of the inserts found
							if insertDict.has_key(deletion[3]):
								insertDict[deletion[3]] = insertDict[deletion[3]] +1
							else:
								insertDict[deletion[3]] = 1
							# if there is no insert mark that down in the dictionary
						if deletion[3] == '':
							insertDict['None'] = insertDict['None'] + 1



								#now that all alignments have been analyzed for deletions, load the deletion dictionary into a list
	uniqueDeletions = list()
	for key in delDict.keys():
		uniqueDeletions.append([key[0], key[1], key[2], delDict[key]])
	#sort the new deletion list by the frequency of deletions
	uniqueDeletions = sorted(uniqueDeletions, key=lambda x: x[3], reverse=True)

	# print deletions to output file
	outputFileName = inputFileName[:-4]
	if args.removeBackwardDels:
		outputFileName += "_NoCircDels"
	if args.MaxInsertLimit != -1:
		outputFileName += "_MaxInsertLimit"+str(args.MaxInsertLimit)
	if args.MinInsertLimit != -1:
		outputFileName += "_MinInsertLimit"+str(args.MinInsertLimit)
	if args.InsertIdentity:
		outputFileName += "_Only"+args.InsertIdentity+"Insert"
	if args.NoAmbiguousDeletions:
		outputFileName += "_NoAmbig"
	if args.OnlyAmbiguousDeletions:
		outputFileName += "_OnlyAmbig"
	if args.ExactMatchAtDelSite:
		outputFileName += "_DelSiteExactMatch"
	outputFileName = outputFileName+ "_deletions.txt"
	outF = open(outputFileName, 'w')
	outF.write("Total Reads Aligned:" + str(totalReadsAligned) + "	Total Deletions:" + str(totalDeletions) + "\n")

	if args.FakeDeletions:
		if args.IndelOutput == 0:
			for line in delByReadList:
				outF.write(line[0] + "\t" + line[1] + "\t" + str(line[2]) + "\t" + str(line[3]) + "\n")
		else:
			for line in delByReadList:
				outF.write(line[0] + "\t" + line[1] + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(line[4]) + "\t" + line[5] + "\n")
	else:
		for line in uniqueDeletions:
			outF.write(line[0] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\n")  # for reporting deletions by sequence location

	#redirect print statements to standard output
	sys.stdout = old_stdout
	log_file.close()
	
	#report deletion statistics
	print args.SamFile[:-4]
	print str(multiDelReads) +" Reads have more than one deletion."
	if args.removeBackwardDels:
		print str(backwardsDels)+" backwards deletions found and removed"
	if args.ExactMatchAtDelSite:
		print str(badMatchCount) + " deletions removed due to a shift of more than 11 nts before finding exact match at deletion site"
	print str(overlappingAlignmentCount) + "or fewer reads had alignments that used overlapping segments from the read and were discounted"
	print str(weirdDoubleSequences) + " reads contain the same sequence more than once and were discounted from consideration"
	#if we're doing indels, write out the final file here
	if args.IndelOutput > 0:
		print "Insert identities and frequencies"
		#convert counts of types to percentages, and determine how many of each length I see
		count1 = 0.0
		count2 = 0.0
		count3 = 0.0
		count4 = 0.0
		count5More = 0.0
		for key in insertDict.keys():
			count = insertDict[key]
			if len(key) == 1:
				count1 += count
			elif len(key) == 2:
				count2 += count
			elif len(key) == 3:
				count3 += count
			elif len(key) == 4:
				if key != "None":
					count4 += count
			else:
				count5More += count
			#prevent a division by zero error by making sure that inserts exist
			if totalIndels > 0:
				insertDict[key] = round(count/totalIndels * 100.0,3)
			else:
				insertDict[key] = 0.00
		print "Frequencies by insert length:"
		if totalIndels > 0:
			print "Length 1 = "+str(round(count1/totalIndels*100.0,2))+"%, Length 2 = "+str(round(count2/totalIndels*100.0,2))+"%, Length 3 = "+str(round(count3/totalIndels*100.0,2))+"%, Length 4 = "+str(round(count4/totalIndels*100.0,2))+"%, Length >=5 = "+str(round(count5More/totalIndels*100.0,2))+"%"
		else:
			print "No Inserts Found"
		print "Frequencies by Identity, in Percent of Total Indels"
		print sorted(insertDict.items(), key=lambda x: (-x[1], x[0]))
		uniqueIndels = list()
		for key in indelDict.keys():
			#divide the total insert length of all indels by the frequency of the deletion to get an average indel rate
			uniqueIndels.append([key[0], key[1], key[2], indelDict[key][0],	 round(float(indelDict[key][1])/indelDict[key][0],1)])
		# sort the new deletion list by the frequency of deletions
		uniqueIndels = sorted(uniqueIndels, key=lambda x: x[3], reverse=True)

		# print deletions to output file

		outF7.write("Total Reads Aligned:" + str(totalReadsAligned) + " Total Indels:" + str(totalIndels) + "\n")

		for line in uniqueIndels:
			# outF.write(line[0]+"\t"+line[1]+"\t"+str(line[2])+"\t"+str(line[3])+"\n") #for when report each reads dels
			outF7.write(line[0] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(line[4]) + "\n")	 # for reporting deletions by sequence location