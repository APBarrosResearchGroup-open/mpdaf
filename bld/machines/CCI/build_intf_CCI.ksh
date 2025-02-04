#!/bin/ksh

if [[ ${numEns} = 1 ]] ; then

cat << EndOfText >mpdaf_run.ksh
#!/bin/ksh
#Job Submission to CCI
#
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=myjob
#SBATCH --partition=secondary
##SBATCH --partition=IllinoisComputes
##SBATCH --account=ps98-ic
#SBATCH --output=myjob.o%j
#SBATCH --error=myjob.e%j
#
# Run
cd $rundir
source $rundir/loadenvs
date

./a.out
exit 0

EndOfText

else

npn=16
let nodes=${numEns}/${npn}

cat << EndOfText >mpdaf_run.ksh
#!/bin/ksh
#Job Submission to CCI
#
#SBATCH --time=01:00:00
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${npn}
#SBATCH --job-name=myjob
#SBATCH --partition=secondary 
##SBATCH --partition=IllinoisComputes
##SBATCH --account=ps98-ic
###SBATCH --array=1-${numEns}
#SBATCH --output=myjob.o%j
#SBATCH --error=myjob.e%j

#
# Run
##cd $origrundir/mpdaf_instance_\$SLURM_ARRAY_TASK_ID
##source $origrundir/loadenvs
date

##./a.out
##exit 0

for (( i=1; i<=${numEns}; i++))
do
#  cd /scratch/users/ps98//mpdaf_CCI_run/mpdaf_instance_\$i
#  source /scratch/users/ps98//mpdaf_CCI_run/loadenvs
  cd $origrundir/mpdaf_instance_\$i
  source $origrundir/loadenvs
  ./a.out &
done
wait

EndOfText

fi

cp $rootdir/bld/machines/${platform}/loadenvs  ${origrundir}
mv mpdaf_run.ksh ${origrundir}
