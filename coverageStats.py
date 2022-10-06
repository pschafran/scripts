#! /home/ps997/miniconda3/bin/python
# Needs to run with python 3 on our server!!!

# Usage: coverageStats.py depth_output_from_samtools.txt

import numpy as np
import sys
from scipy import stats
from scipy import signal
from sklearn.neighbors import KernelDensity

def running_mean(x, N):
	cumsum = np.cumsum(np.insert(x, 0, 0))
	return (cumsum[N:] - cumsum[:-N]) / float(N)

file = open(sys.argv[1],"r")

contigList = []
posList = []
depthList = []
contigDict = {}
assemblySize = 0
contigLengthList = []
indexer = 0
print('Parsing depth file...')
for line in file:
	splitline = line.strip("\n").split("\t")
	contig = splitline[0]
	pos = int(splitline[1])
	depth = float(splitline[2])
	depthList.append(depth)
	try:
		contigDict[contig][0].append(int(pos))
		contigDict[contig][1].append(float(depth))
	except:
		contigDict[contig] = [[],[]]
		contigDict[contig][0].append(int(pos))
		contigDict[contig][1].append(float(depth))

for key in contigDict.keys():
	hiCovList = []
	assemblySize += len(contigDict[key][1])
	contigLengthList.append(len(contigDict[key][1]))
	contigDict[key].append(np.median(contigDict[key][1])) # contigDict[key][2]
	contigDict[key].append(np.std(contigDict[key][1])) # contigDict[key][3]
	contigDict[key].append(np.min(contigDict[key][1])) # contigDict[key][4]
	contigDict[key].append(np.max(contigDict[key][1])) # contigDict[key][5]


depthMedian = np.median(depthList)
depthStd = np.std(depthList)
depthMin = np.min(depthList)
depthMax = np.max(depthList)
depthMode = stats.mode(depthList)[0][0]
depthRange = depthMax - depthMin
depthLowerBound = depthMedian - 100
if depthLowerBound < 0:
	depthLowerBound = 0
depthUpperBound = depthMedian + 100
depthBins = int(depthRange/5)
if depthBins < 100:
	depthBins = int(depthRange)
	if depthBins == 0:
		depthBins = 1

n50size = 0
n50list = []
for length in sorted(contigLengthList, reverse=True):
	if n50size < assemblySize*0.5:
		n50list.append(length)
		n50size += length
	else:
		n50list.append(length)
		n50size += length
n50 = np.min(n50list)

hist = np.histogram(depthList)
kde = KernelDensity(kernel='gaussian').fit(hist)
peaks, peakProperties = signal.find_peaks(kde)
prominences, peakBounds = signal.peak_prominences(kde, peaks)


"Assembly size: %s\n" (assemblySize)
"N50: %s\n" (n50)
"Median depth (total sites): %s\n" %(depthMedian)
"Minimum depth (total sites): %s\n" %(depthMin)
"Maximum depth (total sites): %s\n" %(depthMax)
