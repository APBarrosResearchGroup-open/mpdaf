#!/bin/ksh
#Usage: ./run_obs.ksh
diri="/projects/illinois/eng/cee/barros/snowpits/2017pit/"
rundir=memls_out
#
fils=(`ls ${diri}*.txt`)
#
#rm -rf $rundir
mkdir $rundir
cd $rundir
ctr=0
for fil in ${fils[*]}
do
  echo $fil
  let ctr=${ctr}+1
  cp ../txt2mat.m txt2mat_copy_${ctr}.m
  cp ../matlab.sbatch ./matlab.sbatch
  chmod +w matlab.sbatch
  chmod +rwx txt2mat_copy_${ctr}.m
  sed "s,__fil__,${fil}," -i txt2mat_copy_${ctr}.m 
  sed "s,__mfile__,txt2mat_copy_${ctr}," -i matlab.sbatch
  sbatch matlab.sbatch
done
