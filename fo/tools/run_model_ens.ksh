#!/bin/ksh

diri="/home/ps98/scratch/mpdaf_CCI_run/out/ens_"
diri="/home/ps98/scratch/OUT_GM_OL/ens_"
for iens in {1..48}
do
 rundir=ens_${iens}
# rm -rf $rundir
 mkdir $rundir
 cp out2mat.m $rundir/out2mat_copy.m
 cp matlab.sbatch $rundir
 cd $rundir
 chmod +rwx out2mat_copy.m
 sed "s,__diri__,${diri}${iens}," -i out2mat_copy.m 
 sed "s,__mfile__,out2mat_copy," -i matlab.sbatch
 sbatch matlab.sbatch
 cd ..
done
