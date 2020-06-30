#filterDels.py

#March 9 2015

#This script takes in a crosslinked deletions list and its control. In general, crosslinked deletions
#are filtered by the control. If there are no options, any deletion that occurs in the control is removed from
#the crosslinked set. If the --s option is set, overlapping deletions are subtracted by their normalized del count.

import numpy
import argparse

def parseArgs():
    prs = argparse.ArgumentParser()
    prs.add_argument("test", type=str, help="Crosslinked deletions")
    prs.add_argument("control", type=str, help="deletions from a control like DMSO or NMIA")
    prs.add_argument('--s',action='store_true',default=False,help='Control deletions are subtracted rather than filtered')
    prs.add_argument('--s2',action='store_true',default=False,help='Control deletions are subtracted at twice their actual rate')
    prs.add_argument('--sd',action='store_true',default=False,help='Control deletions are subtracted and then divided by the experimental rate')
    prs.add_argument('--sdc',action='store_true',default=False,help='Control deletions are subtracted and then divided by the control rate')
    prs.add_argument('--rankOrder',action='store_true',help="Subtract the ranks of the deletions rather than their rates. Unique Experimental rates get top ranking.")
    prs.add_argument('--unionRankOrder',action='store_true',help="Similar to rank order but the ranks of the deletions are determined after a union of the experimental and control deletions.")
    o = prs.parse_args()
    return o

def rankOrderDeletions(deletions):
	#sort deletions by rate, highest at top
	deletions = sorted(deletions,key=lambda x: float(x[3]), reverse=True)
	order = range(len(deletions),0,-1)
	for i in range(len(deletions)):
		deletions[i][3] = order[i]
	return deletions

args = parseArgs()

#load in crosslinked dels
inCr = open(args.test,'r')
crossLines = inCr.readlines()
inCr.close()

#load in control dels
inCo = open(args.control,'r')
controlLines = inCo.readlines()
inCo.close()

#set up output file
if(args.s):
	outF = open(args.test[:-4]+"_Subtracted.txt",'w')
elif(args.sd):
	outF = open(args.test[:-4]+"_SubDiv.txt",'w')
elif(args.sdc):
	outF = open(args.test[:-4]+"_SubDivCon.txt",'w')
elif(args.s2):
	outF = open(args.test[:-4]+"_SubtractedDouble.txt",'w')
elif(args.rankOrder):
	outF = open(args.test[:-4]+"_RankOrderSub.txt",'w')
elif(args.unionRankOrder):
	outF = open(args.test[:-4]+"_UnionRankOrderSub.txt",'w')
else:
	outF = open(args.test[:-4]+"_Filtered.txt",'w')

newDelList = list()

if args.rankOrder:
	#load cross link deletion into array
	crossDels = list()
	for i in range(1,len(crossLines)):
		crossLine = crossLines[i].strip().split()
		crossLine.append(0)
		crossDels.append(crossLine)
		
	#do the same for control
	#load cross link deletion into array
	controlDels = list()
	for i in range(1,len(controlLines)):
		controlLine = controlLines[i].strip().split()
		controlLine.append(0)
		controlDels.append(controlLine)
	
	#compare the two deletion sets. If sort crosslinked dels into unique and intersecting
	#remove all control deletions that don't exist in the crosslinked
	crossDelsUnique = list()
	crossDelsInter = list()
	controlDelsInter = list()
	for crossLine in crossDels:
		for controlLine in controlDels:
			if (crossLine[1] == controlLine[1] and crossLine[2] == controlLine[2]):
				crossLine[4] = 1
				controlLine[4] = 1
	

	for crossLine in crossDels:
		if(crossLine[4] == 1):
			crossDelsInter.append(crossLine)
		else:
			crossDelsUnique.append(crossLine)
	for controlLine in controlDels:
		if controlLine[4] == 1:
			controlDelsInter.append(controlLine)
	
	#rank order the deletions
	crossDelsInter = rankOrderDeletions(crossDelsInter)
	controlDelsInter = rankOrderDeletions(controlDelsInter)
		
	#run through the intersection deletions, and subtract their ranks
	for i in range(len(crossDelsInter)):
		crossStart = crossDelsInter[i][1]
		crossStop = crossDelsInter[i][2]		
		for j in range(len(controlDelsInter)):
			controlStart = controlDelsInter[j][1]
			controlStop = controlDelsInter[j][2]
			if(crossStart == controlStart and crossStop == controlStop):
				crossDelsInter[i][3] = crossDelsInter[i][3] - controlDelsInter[j][3]

	#sort the results of the subtraction by their new rate
	crossDelsInter = sorted(crossDelsInter,key=lambda x: float(x[3]), reverse=True)
	#sort the unique dels, merge them with the ordered subtracted dels and apply a new order
	crossDelsUnique = sorted(crossDelsUnique,key=lambda x: float(x[3]), reverse=True)
	

	newDelList = crossDelsUnique + crossDelsInter
	order = range(len(newDelList),0,-1)
	for i in range(len(newDelList)):
		newDelList[i][3] = order[i]

elif(args.unionRankOrder):
	#load cross link deletion into array joint array
	#identify Crosslink deletions with a 1 in column 4
	unionDels = list()
	for i in range(1,len(crossLines)):
		crossLine = crossLines[i].strip().split()
		crossLine.append(1)
		unionDels.append(crossLine)
		
	#do the same for control
	#load cross link deletion into array joint array
	#identify control deletions with a 0 in column 4
	controlDels = list()
	for i in range(1,len(controlLines)):
		controlLine = controlLines[i].strip().split()
		controlLine.append(0)
		unionDels.append(controlLine)
	
	#rank order the union set of deletions
	unionDels = rankOrderDeletions(unionDels)

	#compare all union deletions to each other, if a TBIA and NMIA share coordinates, subtract the crosslink order number from the control
	for i in range(len(unionDels)):
		if(unionDels[i][4] == 1):
			for j in range(len(unionDels)):
				if unionDels[i][1] == unionDels[j][1] and unionDels[i][2] == unionDels[j][2] and unionDels[j][4] == 0:
					unionDels[i][3] = unionDels[i][3] - unionDels[j][3]
					break

	
	#run through the union deletion set and pull all Crosslink deletions into a new deletions list to be printed
	
	for line in unionDels:
		if(line[4] == 1):
			newDelList.append(line)
	#rerank this list of crosslink deletions
	newDelList = rankOrderDeletions(newDelList)


else:
	for i in range(1,len(crossLines)):
	#rund through all crosslink deletions and load data
		crossLine = crossLines[i].strip().split()
		geneName = crossLine[0]
		crossStart = crossLine[1]
		crossStop = crossLine[2]
		crossFreq = float(crossLine[3])
		found = False
		for j in range(1,len(controlLines)):
			#load in control deletion data and compare to crosslinked del
			controlLine = controlLines[j].strip().split()
			controlStart = controlLine[1]
			controlStop = controlLine[2]
			controlFreq = float(controlLine[3])
			if(crossStart == controlStart and crossStop == controlStop):
				found = True
				#adjust frequency accordingly
				if(args.s):
					crossLine[3] = crossFreq - controlFreq
					newDelList.append(crossLine)
				if(args.sd):
					crossLine[3] = (crossFreq - controlFreq)/crossFreq
					newDelList.append(crossLine)
				if(args.sdc):
					crossLine[3] = (crossFreq - controlFreq)/controlFreq
					newDelList.append(crossLine)
				if(args.s2):
					crossLine[3] = crossFreq - (2 * controlFreq)
					newDelList.append(crossLine)
				break
		if(found == False):
			newDelList.append(crossLine)
#re sort deletions by frequency
newDelList = sorted(newDelList,key=lambda x: float(x[3]), reverse=True)
#add in header
newDelList.insert(0,crossLines[0])
#print out altered deletions
outF.write(str(newDelList[0]))
for i in range(1,len(newDelList)):
	line = newDelList[i]
	outF.write(str(line[0])+'\t'+str(line[1])+'\t'+str(line[2])+'\t'+"{:1.10f}".format(float(line[3]))+"\n")
outF.close()