#!/bin/bash
###############################################################################
##                                                                           ##
##                           Campus Cluster                                  ##
##                    Sample MATLAB Job Batch Script                         ##
##                                                                           ##
## SLURM Options (To view, run the following command below)                  ##
##                                                                           ##
##     man sbatch                                                            ##
##                                                                           ##
###############################################################################
#                                                                            
#SBATCH --time=00:30:00                  # Job run time (hh:mm:ss)
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks-per-node=16             # Number of task (cores/ppn) per node
#SBATCH --job-name=matlab_job            # Name of batch job
##SBATCH --partition=secondary            # Partition (queue)           
#SBATCH --partition=IllinoisComputes
#SBATCH --account=ps98-ic
#SBATCH --output=matlab.o%j              # Name of batch job output file
##SBATCH --error=matlab.e%j              # Name of batch job error file
##SBATCH --mail-user=NetID@illinois.edu  # Send email notifications
##SBATCH --mail-type=BEGIN,END           # Type of email notifications to send
#                                                                            
###############################################################################

# Change to the directory from which the batch job was submitted
# Note: SLURM defaults to running jobs in the directory where
# they are submitted, no need for cd'ing to $SLURM_SUBMIT_DIR

#cd ${SLURM_SUBMIT_DIR}
rm *.o*
rm *.mat
#rm *.dat

# Clear the value set in the DISPLAY environment variable
# to run the CLI version of MATLAB
unset DISPLAY

# Load MATLAB module (Enable MATLAB in user environment)
module load matlab

# Run matlab code and redirect output into a file whose
# name includes the JobID.

matlab -nodisplay -r __mfile__ >& ${SLURM_JOB_NAME}.${SLURM_JOB_ID}.oe

