Log: P. Shrestha (Oct 6 2022)
   
**MEMLS** should be _downloaded_ within this folder. The script in tools will access the matlab codes to process
the MSHM ascii files.

#### To run MEMLS offline:  
**/tools** folder contains
 
* _out2mat.m_   - for MSHM outputs
* _txt2mat.m_   - for snow pit profiles

Input parameters are:

* number of grid points, number of layers, number of time steps
* MEMLS parameters (m, q, sh, sv, ssh, ssv, iangle, frequency)

Scripts:

* matlab.sbatch - job submission script for CCI cluster
* run_obs.ksh - run script for snow pit measurements
    needs vertical profile of observations
* run_model.ksh - run script for single model runs
    needs model outputs
* run_model_ens.ksh - run script for ensemble runs
    needs model outputs

#### For running
 * Check the paths, time loop and parameters in .m file
 * ./run_XXXX.ksh
   It will generate memls_out directory and execute MEMLS

 * For data with np>1
  * Band aids
    nt=nt1, np=np1  
    use np = 1  
    and nt = np1 X nt1  
    and for big time loop, use actual it0:nt1:it1  

