#!/usr/bin/env python

import os,re
from collections import OrderedDict

numForAnalysis = 5 #125
#sequenceFileDir = "querySeqFiles52"

file_sbatch = '''#!/bin/bash
#SBATCH --job-name=OTS
#SBATCH --partition=compute
#SBATCH --time=23:00:00
#SBATCH --mail-user="yourEmail@XXX.jp"
#SBATCH --mail-type=FAIL,END
#SBATCH --input=none
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --array=0-NUM
#SBATCH --chdir=WORKDIR

query="$(printf 910_list${SLURM_ARRAY_TASK_ID}).out"
#  echo $query > "910_list${SLURM_ARRAY_TASK_ID}.out"

for fq in `cat ${query}`
do
  ./orthoscope_star.py $fq
done
'''


def readFasta_dict(dirAddressFN, InfileNameFN):
    #print("dirAddressFN + InfileNameFN", dirAddressFN + InfileNameFN)
    #exit()
    Infile = open(dirAddressFN + InfileNameFN, "r")
    seqDictFN  = OrderedDict()
    for Line in Infile:
        Line = Line.rstrip("\n")
        if not Line:
            continue
        elif Line[0] == ">":
            Name            = Line
            seqDictFN[Name] = ""
        else:
            Line = Line.replace("\n", "")
            Line = Line.replace("\r", "")
            #if InfileNameFN == "000_cds_assigned_by_ID.txt":
            #    Line = re.sub("-", "", Line)
            seqDictFN[Name] += Line.upper()
    Infile.close()
    return seqDictFN
#############

f = open("list_geneIDs.txt")
geneIDs = list(f)
f.close()

countInfile = 0
countOutfile = 0
stock = []
prefix_outfile = "910_list"
for geneID in geneIDs:
    geneID = geneID.rstrip("\n")
    if countInfile % numForAnalysis == 0 and countInfile > 0:
        out = open(prefix_outfile + str(countOutfile) + ".out", "w")
        for file1 in stock:
            out.write(file1 + "\n")
        out.close()
        countOutfile += 1;
        stock = [];
        stock.append(geneID)
    else:
        stock.append(geneID)
    countInfile +=1
out = open(prefix_outfile + str(countOutfile) + ".out", "w")
for geneID in stock:
    geneID = geneID.rstrip("\n")
    out.write(geneID + "\n")
out.close()

#outlist = open("030_list.txt", "w")
#for file in geneIDs:
#    outlist.write(file + "\n")
#outlist.close()

out2 = open("920_arrayJobAA.slurm", "w")
file_sbatch1 = re.sub("NUM", str(countOutfile), file_sbatch)
file_sbatch1 = re.sub("WORKDIR", os.getcwd(), file_sbatch1)
out2.write(file_sbatch1)
out2.close()

