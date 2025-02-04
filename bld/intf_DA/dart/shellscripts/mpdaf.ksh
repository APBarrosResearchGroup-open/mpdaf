#!/bin/ksh

#Usage: ./mpdaf.ksh NENS
#Update any namelist in mpdaf ensembles

#User Settings
export WORK="$HOME/scratch/mpdaf_CCI_run"
#
export NENS=$1

#
cd $WORK
# Clean run directory

# Step 1: Update namelist
#----------------------------------------
for instance in {1..$(($NENS))}
do
  echo $instance
  cd mpdaf_instance_$instance
  sed "s,istart=0,istart=336," -i namelist.nml
  sed "s,itime=336,itime=48," -i namelist.nml   #daily
#  sed "s,itime=336,itime=144," -i namelist.nml   #3 day
#  sed "s,itime=336,itime=240," -i namelist.nml   #5 day
#  sed "s,itime=336,itime=480," -i namelist.nml   #10 day
  #istart_old=`grep istart namelist.nml | cut -d'=' -f2 | cut -d',' -f1`
  #echo $istart_old
  cd ..
done

wait

exit 0
