import numpy as np
from Bio import SeqIO
import csv
import os
import re
import math
import matplotlib.pyplot as plt
from collections import OrderedDict

# Setwd
os.chdir("/mywd")

def readFasta(pathToFile, seqCounter, seq_per_protein):
    tmpCounter = 0
    for seq_record in SeqIO.parse(pathToFile, "fasta"):
        prot_id = re.split('\:', (re.split('\|', seq_record.id)[0]))[0]
        if (prot_id in inputDict):
            inputDict[prot_id].append(seq_record)
        else:
            inputDict[prot_id] = [seq_record]
        seqCounter = seqCounter + 1
        tmpCounter += 1
    if (tmpCounter in seq_per_protein):
        seq_per_protein[tmpCounter] += 1
    else:
        seq_per_protein[tmpCounter] = 1
    return seqCounter
    # print(seq_record.id)
    # print(repr(seq_record.seq))
    # print(len(seq_record))


def readOutput(pathToFile, lengthSeq, id_prot):
    outputDict[id_prot] = np.memmap(pathToFile, dtype=np.float32, mode='r', shape=lengthSeq)


pathToInput = "/mywd/bdb_bval/input/"
pathToOutput = "/mywd/bdb_bval/output/"

inputDict = {}
outputDict = {}

seq_per_protein = {}
seqCounter = 0
protCounter = 0

# read input
print("Reading input")
for file in os.listdir(path=pathToInput):
    protCounter = protCounter + 1
    seqCounter = readFasta(pathToInput + file, seqCounter, seq_per_protein)

print("Number of proteins: ", protCounter)
print("Number of chains: ", seqCounter)

# read output
print("Reading output")
for prot_id in enumerate(inputDict.keys()):
    prot_id_lower = str.lower(prot_id[1])
    # print (prot_id[1])
    lenSeq = len(inputDict[prot_id[1]][0])
    readOutput(pathToOutput + prot_id_lower + ".bdb.memmap", lengthSeq=lenSeq, id_prot=prot_id[1])

# plotting
print("plotting")
ordered = OrderedDict(sorted(seq_per_protein.items(), key=lambda t: t[0]))
plt.bar(range(len(ordered)), list(ordered.values()), align='center')
plt.xticks(range(len(ordered)), list(ordered.keys()))
plt.xlabel("Number of chains per protein")
plt.ylabel("Frequency of number of chains per protein")
plt.title("Distribution of chains per protein")
plt.savefig('./chains_per_prot.png')
plt.show()

# histogram for feature distribution
print("plotting histogram for feature distribution")
from itertools import chain

allOutput = list(chain.from_iterable(outputDict.values()))
allOutputSorted = np.sort(allOutput)
x = allOutputSorted[np.logical_not(np.isnan(allOutputSorted))]

plt.hist(x, bins=20)
plt.title("Histogram of B-values")
plt.xlabel("B-value")
plt.ylabel("Frequency")
plt.savefig('./bvals.png')
plt.show()

# histogram for feature distribution of one protein
print("Plotting histogram for feature distribution of one protein")
removedNAN = outputDict['1A1X'][np.logical_not(np.isnan(outputDict['1A1X']))]
plt.hist(removedNAN, bins=10)
plt.title("Histogram for B-values of the protein 1A1X")
plt.xlabel("B-values")
plt.ylabel("Frequency")
plt.savefig('./bval_1A1X.png')
plt.show()


#Pie chart plot of features
weak = 0
strong = 0
normal = 0

for bval in x:
    if bval < -1:
        weak += 1
    elif bval > 1:
        strong += 1
    else:
        normal += 1


labels = 'Normal', 'Weak', 'Strong'
colors = ['darkolivegreen', 'indianred', 'goldenrod']
sizes = normal, weak, strong
plt.pie(sizes, labels=labels, colors=colors, startangle=50, autopct='%1.1f%%')
plt.axis('equal')
#   plt.title("Distribution of B-values")
plt.savefig('./bval_pie.png')
plt.show()

############### Import Prot tVec ######################

protVec = {}

with open('protVec_100d_3grams.csv', newline='') as protVec_file:
    prot_reader = csv.reader(protVec_file, delimiter='\t')

    # Skips first line
    next(prot_reader)

    for row in prot_reader:
        three_mer = row[0]
        variables = row[1:]
        if three_mer in protVec:
            print("ERROR: 3mer occurs multiple times!!!")
            print(three_mer)
        else:
            protVec[three_mer] = variables


