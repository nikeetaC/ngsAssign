import random
import re
import random
import os
import gzip
from Bio import  SeqIO
import sys

if(len(sys.argv) < 1):
    exit("Usage: python3.6 generateFastq.py genomeFile")
genomeFile = sys.argv[1]

#Reading reference genome
print("Reading genome sequence")
genomeDict = dict()
for record in SeqIO.parse(genomeFile, "fasta"):
    genomeDict[record.id] = str(record.seq)

#Creating .fastq file from reference genome
print("Making fastq file of 50bp reads")
dictPhredRef = {20:'5'}
if(os.path.isfile('output.fastq') == True):
    os.remove('output.fastq')
f=open('output.fastq', 'a')

#number of contigs present in reference
noOfContigs = len(genomeDict.keys())
#loop over
for i in range(100000):
    randomContig = list(genomeDict.keys())[random.randint(0,noOfContigs-1)] #select random contig
    contigSeq = str(genomeDict[randomContig]) #contigseq
    startPos = random.randrange(0,len(contigSeq)-50) #random start
    endPos = startPos + 50 # endpos
    strId = '@seqID:' + randomContig + ":" + str(startPos) + ':' + str(endPos) + '\n'
    f.write(strId)
    f.write(contigSeq[startPos:endPos] + '\n')
    f.write('+\n')
    for j in range(0, 50):
        f.write(dictPhredRef[20])
    f.write('\n')

print("Mapping mock data on reference genome")
cmdBWA2 = 'bwa mem ' + genomeFile + ' output.fastq -o alignOut.sam'
os.system(cmdBWA2)

print("Checking Error rate:")
errorCount = 0
with open("alignOut.sam","r") as f:
    for i in f.readlines():
        if(i.startswith("seqID")):
            queryStart = int(i.split("\t")[0].split(":")[2]) + 1
            alignStart = i.split('\t')[3]
            if(str(queryStart) != alignStart):
                errorCount+=1

print("Error Rate:" + str(errorCount/100000) + '\n' + 'Error Count: ' + str(errorCount))
