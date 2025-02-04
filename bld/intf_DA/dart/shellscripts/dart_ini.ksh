#!/bin/ksh

#Usage: sbatch dart_ini.ksh NENS
#NENS - total ensemble number

#Job Submission to CCI
#
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=job_dart
#SBATCH --partition=IllinoisComputes
##SBATCH --partition=secondary
#SBATCH --account=ps98-ic
#SBATCH --output=job_dart.o%j
#SBATCH --error=job_dart.e%j

#User Settings

export NENS=$1
export WORK="$HOME/scratch/mpdaf_CCI_run"
export DART_DIR="$HOME/DART/models/mshm/work"
export DART_OBS_DIR="$HOME/DART/observations/obs_converters/text_snex/work"
#

cd $WORK

# Clean run directory
rm dart_*.nc 
rm *.out
rm *.nml
rm filter*
rm dart_to_mshm
rm mshm_to_dart

# Copy DART executables and namelist
cp $DART_DIR/input.nml .
cp $DART_DIR/dart_to_mshm .
cp $DART_DIR/mshm_to_dart .
cp $DART_DIR/filter .
cp $DART_DIR/filter_input_list.txt .
cp $DART_DIR/filter_output_list.txt .
cp $DART_DIR/advance_model.csh .

# Copy observations in DART space
cp $DART_OBS_DIR/*.out .

# For inflation
sed "s,inf_initial_from_restart    = .true.,inf_initial_from_restart    = .false.," -i input.nml
sed "s,inf_sd_initial_from_restart = .true.,inf_sd_initial_from_restart = .false.," -i input.nml

#Step 1: MSHM to DART
#----------------------------------------
for instance in {1..$(($NENS))}
do
  echo $instance
  cd mpdaf_instance_$instance

  rm input.nml
  rm restartMSHM
  #rm dart_prior*.nc

  ln -s ../input.nml .
  mshm_restart=`ls -1 restartMSHM* | tail -n -1`
  echo "Ens no " $instance " Using " $mshm_restart
  ln -s $mshm_restart restartMSHM
  ../mshm_to_dart || exit 1
  dfilename="dart_prior_"$instance".nc"
  #This is needed for the 1st initial run, will be overwritten by filter next time
  mv dart_assim_time.txt ../
  #
  mv dart_prior.nc ../${dfilename}
  cd ..
done

wait

#Step 2: filter
#----------------------------------------

date
rm dart_prior.nc
rm obs_seq.out

#ofilename="data_snex_"`cat dart_assim_time.txt`".out"
ofilename="data_snex_"`cat dart_assim_time.txt`".out"
#Need this file to read geo data
ln -s ${dfilename} dart_prior.nc
ln -s ${ofilename} obs_seq.out
#
mpirun ./filter || exit 3  >> log_file 2>> err_file
date

wait

#Step 3: DART to MSHM
#----------------------------------------

for instance in {1..$(($NENS))}
do
  echo $instance
  cd mpdaf_instance_$instance

#  mshm_new_restart="mshmRestart_"$instance
# If we want to preserve old restart

  rm dart_prior.nc
  rm dart_posterior.nc
  rm ${mshm_new_restart}
  
  dfilename="dart_prior_"$instance".nc"
  ln -s ../${dfilename} dart_prior.nc
  dfilename="dart_posterior_"$instance".nc"
  ln -s ../${dfilename} dart_posterior.nc

  ../dart_to_mshm || exit 1
  mv mshm_update_restart ${mshm_restart} 

  cd ..
done


exit 0
#wait
