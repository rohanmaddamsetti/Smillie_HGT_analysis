#!/bin/bash
#PBS -j oe
#PBS -l walltime=04:00:00,nodes=1:ppn=1,feature=gbe
#PBS -t 1-2
#PBS -N rohan_blastrun

##NOTE: I never actually used this script for anything! It was easier just to run
## my command-line blast on kingshorses rather than on the HPC.
## I've kept this script just for reference in case I want to modify this script
## to run future jobs on the HPC using a script that looks like this.

module load powertools
module load BLAST+/2.2.25
cd ${PBS_O_WORKDIR}
##sleep ${PBS_ARRAYID} ##This is for preventing collisions.

INPUT=./input/*.fasta
for file in $INPUT ##loop through input files.
do
    name=`basename $file .fasta`
    if [ ! -f ./output/$name.out ]
    then
        touch ./output/$name.out
	blastn -task megablast -db REL606db -query $file -out ./output/$name.out	  
    fi
     ##do some kind of check on the remaining walltime
     ###safetyfactor = 10*60 sec
    walltime=`qstat -f ${PBS_JOBID} | grep used.walltime | cut -d "=" -f 2`
    echo $walltime

    #currenttime=`timeconvert $walltime`
    #max=$((4*60*60-10*60))  
    #echo $walltime
#    echo $currenttime
    #if [ $currentime -ge $max ] 
    #then
    #    qsub -t ${PBS_ARRAYID} -N $PBS_JOBNAME $0
#	exit
#    fi
done

qstat -f ${PBS_JOBID} ##prints out statistics on the job run.
