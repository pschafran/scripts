#! /usr/bin/env python3
# Needs to run with python 3!

# Usage: plotCovAlongContigs.py -i depth_output_from_samtools.txt [options]

help ='''Usage:
plotCovAlongContigs.py -i depth_output_from_samtools.txt [options]

Options:
-i, --input	text file containing depth data in 3-column format (i.e. output from samtools depth)
-w, --window-size	size of sliding window to average depth over (default: 0.05% of genome size)
-y, --y-min	minimum y-axis value for coverage plots (default: 0)
-Y, --y-max	maximum y-axis value for coverage plots (default: maximum value of the data)
-x, --x-min	minimum x-axis value for density plots (default: 0)
-X, --x-max	maximum x-axis value for density plots (default: median coverage + 100)
-b, --bin-width	histogram bin size (default: 5)
'''

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from scipy import stats
from matplotlib import colors
from numpy import log
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
import random

def running_mean(x, N):
	cumsum = np.cumsum(np.insert(x, 0, 0))
	return (cumsum[N:] - cumsum[:-N]) / float(N)

if "-h" in sys.argv:
	print(help)
	exit(0)
if len(sys.argv) < 2:
	print("ERROR: too few arguments")
	print(help)
	exit(1)

input = False
avgingWindowSize = False
yMin = False
yMax = False
xMin = False
xMax = False
binWidth = 5

index = 0
for item in sys.argv:
	if item in ["-i","--input"]:
		input = sys.argv[index+1]
	elif item in ["-w", "--window-size"]:
		avgingWindowSize = int(sys.argv[index+1])
	elif item in ["-y","--y-min"]:
		yMin = int(sys.argv[index+1])
	elif item in ["-Y","--y-max"]:
		yMax = int(sys.argv[index+1])
	elif item in ["-x","--x-min"]:
		xMin = int(sys.argv[index+1])
	elif item in ["-X","--x-max"]:
		xMax = int(sys.argv[index+1])
	elif item in ["-b","--bin-width"]:
		binWidth = float(sys.argv[index+1])
	index+=1

if input == False:
	print("ERROR: No input file provided")
	print(help)
	exit(1)
filename = input.split(".txt")[0]

file = open(input,"r")
#if smoothingFactor < 1:
#	print("WARNING: Smoothing was set to less than 1. Now set to 1.")
#	smoothingFactor = 1
#if smoothingFactor > 100:
#	print("WARNING: Smoothing was set to greater than 100. Now set to 100.")
#	smoothingFactor = 100

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
depthListRand = random.sample(depthList, k=10000)
#print('Finding high coverage outlier regions...')
#for key in contigDict.keys():
#	peaks, _ = find_peaks(contigDict[key][1], prominence = 1)
#	baseline = np.array(contigDict[key][1])
#	contigDict[key].append(peaks)
for key in contigDict.keys():
#	hiCovList = []
	assemblySize += len(contigDict[key][1])
	contigLengthList.append(len(contigDict[key][1]))
	contigDict[key].append(np.median(contigDict[key][1])) # contigDict[key][2]
	contigDict[key].append(np.std(contigDict[key][1])) # contigDict[key][3]
	contigDict[key].append(np.min(contigDict[key][1])) # contigDict[key][4]
	contigDict[key].append(np.max(contigDict[key][1])) # contigDict[key][5]
#	if any(x > contigDict[key][2]+(3*contigDict[key][3]) for x in contigDict[key][1]):
		#print('	%s' %(key))
		#outfile = open("%s_hiCovRegions.txt" %(key), "w")
		#outfile.write("position\tcoverage\n")
#		for i in range(0, len(contigDict[key][0])):
#			if contigDict[key][1][i] > contigDict[key][2]+(3*contigDict[key][3]):
#				hiCovList.append(int(contigDict[key][0][i]))
		#		outfile.write("%d\t%f\n" %(contigDict[key][0][i-1], contigDict[key][1][i-1]))
		#outfile.close()
#		outliers = np.array(hiCovList)
#		contigDict[key].append(contigDict[key][2]+(3*contigDict[key][3])) # contigDict[key][6]

depthMedian = np.median(depthList)
depthStd = np.std(depthList)
depthMin = np.min(depthList)
depthMax = np.max(depthList)
depthMode = stats.mode(depthList)[0][0]
#depthRange = depthMax - depthMin
if xMin == False:
	depthLowerBound = depthMedian - 100
	if depthLowerBound < 0:
		depthLowerBound = 0
else:
	depthLowerBound = xMin
if xMax == False:
	depthUpperBound = depthMedian + 100
else:
	depthUpperBound = xMax
depthRange = depthMax- depthMin
#depthBins = int(depthRange/5)
#if depthBins < 100:
#	depthBins = int(depthRange)
#	if depthBins == 0:
#		depthBins = 1

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
#n50 = 0
print("Assembly Size: %s" %(assemblySize))
if avgingWindowSize == False:
	avgingWindowSize = int(assemblySize * 0.005)
movAvgDict = {}
print("Window Size: %s" %(avgingWindowSize))
print('Smoothing data...')
runningyMax = 0
runningyMin = 0
for key in contigDict.keys():
	#smoothingBases = int(len(contigDict[key][1]) * (smoothingFactor/100))
	#movAvgDict[key] = [np.array(running_mean(contigDict[key][0], avgingWindowSize)).tolist(), np.array(running_mean(contigDict[key][1], avgingWindowSize)).tolist()]
	movAvgList = np.convolve(contigDict[key][1], avgingWindowSize, mode = "valid").tolist()
	movAvgDict[key] = [np.array(contigDict[key][0], movAvgList)]
	if yMin == False and runningyMin > min(movAvgDict[key][1]):
		runningyMin = int(min(movAvgDict[key][1]))
	if yMax == False and runningyMax < max(movAvgDict[key][1]):
		runningyMax = int(max(movAvgDict[key][1]))
if yMin == False:
	yMin = runningyMin
if yMax == False:
	yMax = runningyMax

print('Plotting charts...')
print("xMin: %s" % depthLowerBound)
print("xMax: %s" % depthUpperBound)
print("yMin: %s" % yMin)
print("yMax: %s" % yMax)
fig, ax = plt.subplots(1,1)
ax.hist(depthListRand, bins = np.arange(int(min(depthListRand)), int(max(depthListRand)) + binWidth , binWidth), density = True)
ax.set_xlabel("Read Depth")
ax.set_ylabel("Density")
ax.set_title("Read Depth Across All Sites")
plt.grid(which='major',axis='both',linestyle='dashed')
ax.set_axisbelow(True)
ax.set_xlim(depthLowerBound, depthUpperBound)
ax.text(0.95,0.95, "Median: %d" %(depthMedian), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
ax.text(0.95,0.9, "Mode: %d" %(depthMode), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
plt.savefig("%s_depth_hist.pdf" %(filename),format = "pdf")
plt.close()

for key in contigDict.keys():
	baseline = np.array(contigDict[key][1])
	fig, ax = plt.subplots(1,1)
	ax.plot(movAvgDict[key][0], movAvgDict[key][1], color = "blue", linewidth = 0.5)
	#ax.plot(movAvgDict[key][0], movAvgDict[key][1], color = "blue")
	ax.hlines(contigDict[key][2]+contigDict[key][3], xmin = np.min(contigDict[key][0]), xmax = np.max(contigDict[key][0]), color = "yellow", label = "1 SD")
	ax.hlines(contigDict[key][2]+(2*contigDict[key][3]), xmin = np.min(contigDict[key][0]), xmax = np.max(contigDict[key][0]), color = "orange", label = "2 SD")
	ax.hlines(contigDict[key][2]+(3*contigDict[key][3]), xmin = np.min(contigDict[key][0]), xmax = np.max(contigDict[key][0]), color = "red", label = "3 SD")
	ax.set_xlabel("Position")
	ax.set_ylabel("Read Depth")
	ax.set_title("%s" %(key))
	plt.grid(which='major',axis='both',linestyle='dashed')
	#ax.set_ylim(ymin=(contigDict[key][4] - 100), ymax = (contigDict[key][5] + 100))
	#ax.set_yscale('log')
	ax.set_axisbelow(True)
	ax.set_ylim(ymin = yMin,ymax = yMax)
	#ax.text(0.7, 1.1, "Median Depth: %d" %(contigDict[key][2]), horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
	#ax.text(0.7, 1.05, "Standard Deviation: %d" %(contigDict[key][3]), horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
	plt.savefig("%s_%skbp.pdf" %(key,(avgingWindowSize/1000)),format = "pdf")
	plt.close()

for key in contigDict.keys():
	contigMin = min(contigDict[key][1])
	contigMax = max(contigDict[key][1])
	fig, ax = plt.subplots(1,1)
	ax.hist(random.sample(contigDict[key][1],k=10000), bins = np.arange(int(contigMin), int(contigMax + binWidth), binWidth), density = True)
	ax.set_xlabel("Read Depth")
	ax.set_ylabel("Density")
	ax.set_title("%s" %(key))
	plt.grid(which='major',axis='both',linestyle='dashed')
	ax.set_axisbelow(True)
	ax.set_xlim(xmin = depthLowerBound, xmax = depthUpperBound)
	ax.text(0.95,0.95, "Median: %d" %(np.median(contigDict[key][1])), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
	ax.text(0.95,0.9, "Mode: %d" %(stats.mode(contigDict[key][1])[0][0]), horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
	plt.savefig("%s_depth_hist.pdf" %(key),format = "pdf")
	plt.close()

if len(contigDict.keys()) > 1:
	fig, ax = plt.subplots(len(contigDict.keys()), 1, sharex=True, sharey=True, tight_layout=True, figsize = (8,len(contigDict.keys())))
	index = 0
	for key in movAvgDict.keys():
		halfway = int(len(contigDict.keys())/2)
		#if len(contigDict[key][0]) >= n50:
		ax[index].plot(movAvgDict[key][0], movAvgDict[key][1])
		ax[index].text(1.01, 0.5, key, horizontalalignment='left', verticalalignment='center',transform=ax[index].transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
		ax[index].spines['right'].set_visible(False)
		ax[index].spines['top'].set_visible(False)
		ax[index].set_ylim(ymin = yMin, ymax = yMax)
		if index == halfway:
			ax[index].set_ylabel("Read Depth")
		index += 1
	ax[index-1].set_xlabel("Position")
	plt.savefig("%s_all_contigs_%skbp.pdf" %(filename, (avgingWindowSize/1000)),format = "pdf")
	plt.close()

	#fig, ax = plt.subplots(len(n50list), 1, sharex=True, sharey=True, tight_layout=True, figsize = (8,len(n50list)))
	#index = 0
	#for key in movAvgDict.keys():
	#	halfway = int(len(n50list)/2)
	#	if len(contigDict[key][0]) >= n50:
	#		baseline = np.array(movAvgDict[key][1])
	#		ax[index].plot(movAvgDict[key][0], movAvgDict[key][1])
	#		#ax[index].plot(contigDict[key][6], baseline[contigDict[key][6]], color = "red")
	#		ax[index].text(1.01, 0.5, key, horizontalalignment='left', verticalalignment='center',transform=ax[index].transAxes, bbox=dict(facecolor='white', edgecolor='white', alpha=0.5))
	#		ax[index].spines['right'].set_visible(False)
	#		ax[index].spines['top'].set_visible(False)
	#		ax[index].set_yscale('log')
	#		if index == halfway:
	#			ax[index].set_ylabel("Read Depth")
	#		index += 1
	#ax[index-1].set_xlabel("Position")
	#plt.savefig("%s_all_contigs_logCov_%skb.pdf" %((sys.argv[1].split(".txt")[0]), (avgingWindowSize/1000)),format = "pdf")
	#plt.close()
        #
	#fig = plt.figure()
	#ax = fig.gca(projection='3d')
	#index = 0
	#for key in movAvgDict.keys():
	#	if len(contigDict[key][0]) >= n50:
	#		ax.plot(movAvgDict[key][0], movAvgDict[key][1], zs=index, zdir='y', label= "%s" %(key))
	#		index -= 1
	#plt.savefig("%s_3Dn50contigs_%skb.pdf" %(sys.argv[1].split(".txt")[0], (avgingWindowSize/1000)),format = "pdf")
	#plt.close()
