#!/bin/ksh

#Usage: ./archive.ksh
export RUN_DIR="$HOME/scratch/mpdaf_CCI_run/"
NENS=$1

#DART run
cd ${RUN_DIR}
mkdir nc
mkdir dp
mkdir pm
mkdir am 
mkdir om

mv dart_*.nc dp/
mv analysis_*.nc am/
mv preassim_* pm
mv output_* om

mv dp nc/
mv am nc/
mv pm nc/
mv om nc/

#Model run
mkdir out
for instance in {1..$(($NENS))}
do
  echo $instance
  cd mpdaf_instance_$instance
  mkdir ens_$instance
  mv zongSWE_* ens_$instance
  mv zongSD_* ens_$instance
  mv SWE_*  ens_$instance
  mv Depth_* ens_$instance
  mv Density_* ens_$instance
  mv Dsnow_* ens_$instance
  mv LW_* ens_$instance
  mv Tsnow_* ens_$instance
  mv Tsoil_* ens_$instance
  mv ens_$instance ../out
  cd ..
done

#Save dart diagnostics
mkdir dartdiag
mv obs_seq.final_* dartdiag/
mv dartdiag out/

exit 0
