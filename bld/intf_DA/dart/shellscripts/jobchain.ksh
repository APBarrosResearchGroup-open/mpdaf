#!/bin/ksh
#Usage: sbatch jobchain.ksh NENS STEP

#SBATCH --job-name="mpdafDART"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=job_mpdaf-out.%j
#SBATCH --error=job_mpdaf-err.%j
#SBATCH --time=00:05:00
##SBATCH --partition=secondary
#SBATCH --partition=IllinoisComputes
#SBATCH --account=ps98-ic

#User Settings
#-------------------------------------------------
SPATH="$HOME/mpdaf/bld/intf_DA/dart/shellscripts"
WPATH="$HOME/scratch/mpdaf_CCI_run"
JOBSCRIPT0="jobchain.ksh"
JOBSCRIPT1="mpdaf_run.ksh"
JOBSCRIPT2="dart.ksh"
#-------------------------------------------------


NENS=$1
STEP=$2

if [[ $STEP == "run" ]] then
  echo " MPDAF run ...."
  echo " "
  JOBID=$(sbatch $WPATH/$JOBSCRIPT1 2>&1 | awk '{print $(NF)}')
  echo $JOBID
  echo "sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $NENS filter"
  JOBID=$(sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $NENS "filter" 2>&1 | awk '{print $(NF)}')
fi

if [[ $STEP == "filter" ]] then
  echo " Running DART filter ..."
  echo " "
  JOBID=$(sbatch $SPATH/$JOBSCRIPT2 $NENS 2>&1 | awk '{print $(NF)}')
  echo $JOBID
  echo "sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $NENS run"
  JOBID=$(sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $NENS "run" 2>&1 | awk '{print $(NF)}')
fi

exit 0
