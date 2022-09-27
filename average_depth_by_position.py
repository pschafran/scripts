#! /usr/bin/env python

import sys
import statistics

infile = sys.argv[1]

positionDict = {}

with open(infile,"r") as openinfile:
  for line in openinfile:
    splitline = line.strip("\n").split("\t")
    position = splitline[2]
    depth = splitline[3]
    try:
      positionDict[position].append(depth)
    except:
      positionDict = {position : [depth]}

with open("%s_average_by_position.txt" % infile , "r") as outfile:
  for key in positionDict:
    average = numpy.mean(positionDict[key])
    outfile.write("%s\t%s\n" %(key,average))
