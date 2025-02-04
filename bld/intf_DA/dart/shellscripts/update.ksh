#!/bin/ksh

#Usage: sbatch update.ksh NENS
#NENS - total ensemble number

#Job Submission to CCI
#
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=job_update
##SBATCH --partition=secondary 
#SBATCH --partition=IllinoisComputes
#SBATCH --account=ps98-ic
#SBATCH --output=job_update.o%j
#SBATCH --error=job_update.e%j

#User Settings

export NENS=$1
export WORK="$HOME/scratch/mpdaf_CCI_run"
#
rm job_mpdaf-* 
rm myjob.*

cd $WORK

#Step 1: update 
#----------------------------------------

for instance in {1..$(($NENS))}
do
  echo $instance
  cd mpdaf_instance_$instance


 #rm *_*.out
  
  # Update namelist for restart run
  istart_old=`grep istart namelist.nml | cut -d'=' -f2 | cut -d',' -f1`
  restart_freq=`grep itime namelist.nml | cut -d'=' -f2 | cut -d',' -f1`
  istart_new=$(( $istart_old + $restart_freq ))
  sed "s,istart=${istart_old},istart=${istart_new}," -i namelist.nml
  
  cd ..
done

exit 0
#wait
