#$ -S /bin/bash
#$ -cwd
#$ -t 1-7
#$ -tc 7
#$ -l s_vmem=32G
#$ -l mem_req=32G
#$ -N an_array_job

echo running on `hostname`
echo starting at
date
echo -e ""
#echo SGE_TASK_ID $SGE_TASK_ID
#SGE_TASK_ID=`expr $SGE_TASK_ID - 1`
echo SGE_TASK_ID $SGE_TASK_ID

PADDING_N=$(printf "%03d" ${SGE_TASK_ID})
echo PADDING_N $PADDING_N

for gene in `cat ./list_geneIDs_split/list_geneIDs-${PADDING_N}.txt`
do
        echo ########### $gene ###########
        python3 /home/jun-inoue/ORTHOSCOPE_STAR-1.1.5/orthoscope_star_v1.1.6.py ${gene}
done

echo -e ""
echo ending at
date

