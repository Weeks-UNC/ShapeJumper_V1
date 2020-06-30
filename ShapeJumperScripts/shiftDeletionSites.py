#shiftDels.py

#Tom Christy

#June 15, 2017

#this script takes in a deletion file, either top dels or all dels, and shifts the 5' or 3' end by a user specified
#amount. While both shifts are optional, you must select at least one to actually have the deletions shift

#if, after the shift a deletion is less than 10 nt long, it is removed from the data set.

import argparse
import re

def parseArgs():
	prs = argparse.ArgumentParser()
	prs.add_argument("dels",type=str, help = "Name of deletion file that you want to shift")
	prs.add_argument("--shift5nt",type = int, default = 0, help = "Number of nucleotides by which you want to shift the 5' end of the deletion, can be positive or negative")
	prs.add_argument("--shift3nt",type = int, default = 0, help = "Number of nucleotides by which you want to shift the 3' end of the deletion, can be positive or negative")
	o = prs.parse_args()
	return o

args = parseArgs()

#check and see if the first line in test
#load the deletions into a list and iterate through.
inF = open(args.dels,'r')
delLines = inF.readlines()
inF.close()

#determine title and open up output file
title = args.dels[:-4]
if(args.shift5nt != 0):
	title += "_shift5nt"+str(args.shift5nt)

if(args.shift3nt != 0):
	title += "_shift3nt"+str(args.shift3nt)
title += ".txt"
outF = open(title, 'w')
outF.write(delLines.pop(0))

for line in delLines:
	cols = line.strip().split()
	if(len(cols) == 4):
		nt1 = int(cols[1])
		nt2 = int(cols[2])
		#determine which nt is 5'
		nt1 = nt1 + args.shift5nt
		nt2 = nt2 + args.shift3nt
		if(abs(nt2 - nt1) > 10 and nt1 > 0):
			cols[1] = str(nt1)
			cols[2] = str(nt2)
			s = '\t'.join(cols)
			outF.write(s+"\n")
	elif(len(cols) == 3):
		nt1 = int(cols[0])
		nt2 = int(cols[1])
		#determine which nt is 5'
		nt1 = nt1 + args.shift5nt
		nt2 = nt2 + args.shift3nt
		if(abs(nt2 - nt1) > 10 and nt1 > 0):
			cols[0] = str(nt1)
			cols[1] = str(nt2)
			s = '\t'.join(cols)
			outF.write(s+"\n")
	else:
		outF.write(line)
	
	

