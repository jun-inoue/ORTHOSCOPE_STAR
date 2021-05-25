#!/usr/local/bin/python3

import os,re
from collections import OrderedDict

numAnalysis_eachNode = 5 # 125

################
file_sbatch = '''#!/bin/bash
#SBATCH --job-name=yuraiP
#SBATCH --partition=compute
#SBATCH --time=23:00:00
#SBATCH --mail-user="yourEmail@xxxx.ac.jp"
#SBATCH --mail-type=FAIL,END
#SBATCH --input=none
#SBATCH --mem=8G          # memory for tree search
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --array=0-NUM
#SBATCH --chdir=WORKDIR

query="$(printf 010_list${SLURM_ARRAY_TASK_ID}).out"
#  echo $query > "010_list${SLURM_ARRAY_TASK_ID}.out"

for fq in `cat ${query}`
do
  ./orthoscope_star.py $fq
done
'''
#######
def make_file(count_file_fn, stock_geneID_fn):
    out = open("010_list" + str(count_file_fn) + ".txt", "w")
    for geneID in stock_geneID_fn:
        out.write(geneID + "\n")
    out.close()

#######

f = open ("list_geneIDs.txt")
geneIDs = list(f)
f.close()

count_geneID = 1
count_file = 0
stock_geneID = []
for geneID in geneIDs:
    geneID = geneID.rstrip("\n")
    #print(count_geneID)
    stock_geneID.append(geneID)
    if count_geneID == 1:
        count_geneID += 1
    elif count_geneID % numAnalysis_eachNode == 0:
        make_file(count_file, stock_geneID)
        stock_geneID = []
        count_file += 1
        count_geneID += 1
    else:
        count_geneID += 1
#print(stock_geneID)
make_file(count_file, stock_geneID)

out2 = open("920_arrayJobAA.slurm", "w")
file_sbatch1 = re.sub("NUM", str(count_file), file_sbatch)
out2.write(file_sbatch1)
out2.close()
