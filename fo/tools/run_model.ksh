#!/bin/ksh

 diri="/home/ps98/scratch/OUT_GM_CTRL/"
 rundir=memls_out
 rm -rf $rundir
 mkdir $rundir
 cp out2mat.m $rundir/out2mat_copy.m
 cp matlab.sbatch $rundir
 cd $rundir
 chmod +rwx out2mat_copy.m
 sed "s,__diri__,${diri}," -i out2mat_copy.m 
 sed "s,__mfile__,out2mat_copy," -i matlab.sbatch
 sbatch matlab.sbatch
