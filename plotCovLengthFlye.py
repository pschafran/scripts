#! /home/ps997/miniconda3/bin/python
# Needs to run with python 3 on our server!!!


# Usage: plotCovLengthFlye.py assembly_info.txt FigureTitle(also output file name -- no spaces!)

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from matplotlib import colors
from numpy import log


assembly_info = open(sys.argv[1], "r")
title = sys.argv[2]
contigDict = {}
covList = []
lengthList = []
logCovList = []
logLengthList = []
circularList = []
repeatList = []
multiplicityList = []

linecount = 0
for line in assembly_info:
	if linecount > 0:
		splitline = line.strip("\n").split("\t")
		seqName = splitline[0]
		length = float(splitline[1])
		cov = float(splitline[2])
		circular = splitline[3]
		repeat = splitline[4]
		multiplicity = float(splitline[5])
		graphPath = splitline[6]
		
		contigDict[seqName] = {"length":length}
		contigDict[seqName].update({"cov":cov})
		contigDict[seqName].update({"circular":circular})
		contigDict[seqName].update({"repeat":repeat})
		contigDict[seqName].update({"multiplicity":multiplicity})
		contigDict[seqName].update({"graphPath":graphPath})
		
		lengthList.append(length/1000000)
		covList.append(cov)
		logLengthList.append(log(length))
		logCovList.append(log(cov))
		if circular == "+":
			circularList.append(1)
		elif circular == "-":
			circularList.append(0)
		if repeat == "+":
			repeatList.append(1)
		elif repeat == "-":
			repeatList.append(0)
		multiplicityList.append(multiplicity)
		
		linecount += 1
	else:
		linecount += 1

assembly_info.close()


assemblySize = np.sum(lengthList)
sizeCounter = 0
n50 = assemblySize*0.75
for i in sorted(lengthList, reverse=True):
	if sizeCounter < n50:
		n50size = i
		sizeCounter += i
	else:
		pass

n50covList = []
for key in contigDict.keys():
	if contigDict[key]["length"]/1000000 >= n50size:
		n50covList.append(contigDict[key]["cov"])

lowerQuartile = np.quantile(n50covList, 0.25)
median = np.median(n50covList)
upperQuartile = np.quantile(n50covList, 0.75)
IQR = upperQuartile - lowerQuartile
lowerOutlier = IQR - (1.5*IQR)
upperOutlier = IQR + (1.5*IQR)

standardDev = np.std(n50covList)
upperOutlier = median + (3*standardDev)
lowerOutlier = median - (3*standardDev)


ax = plt.gca()
plt.grid(which='major',axis='y',linestyle='--')
ax.scatter(lengthList, covList, c = 'blue', alpha = 0.5, s=5, marker = 'o')
ax.set_yscale('log')
ax.set_xlabel('Length (Mb)')
ax.set_ylabel('Coverage (X)')
ax.set_title(title)
plt.axvline(x=0.25, ymin = 0, ymax = 1, linewidth=0.333, color='black')
plt.axhline(y=lowerOutlier, xmin = 0, xmax = 1, linewidth=0.333, color='red', alpha = 0.5)
#plt.axhline(y=median, xmin = 0, xmax = 1, linewidth=0.333, color='black', alpha = 0.5)
plt.axhline(y=upperOutlier, xmin = 0, xmax = 1, linewidth=0.333, color='red', alpha = 0.5)
for index in range(1,len(lengthList)):
	if lengthList[index-1] < 0.25:
		pass
		
	elif lengthList[index-1] > 0.25 and covList[index-1] < upperOutlier and covList[index-1] > lowerOutlier:
		pass
	else:
		plt.annotate(s = " ",xy = (lengthList[index-1],covList[index-1]), xytext = (15,0), textcoords ='offset points',ha='center', va='center', arrowprops=dict(facecolor='black', headwidth = 4, headlength = 4, width = 1, shrink = 0.15, alpha = 0.5))


ax.set_axisbelow(True)
plt.savefig("%s.pdf" %(title),format = "pdf")
plt.close()