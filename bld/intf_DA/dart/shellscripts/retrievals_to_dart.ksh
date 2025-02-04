#!/bin/ksh

# Convert text data to DART space, one snapshot data at multiple locations 
# Users need to choose the observation errors, 
# Also generates retrievals_out.txt file which contains the list of files generated 

export RTV_DIR=$HOME/DART/observations/obs_converters/text_snex/
export RTV_DATA=${RTV_DIR}/data
export RTV_WORK=${RTV_DIR}/work

#Cleanup
rm ${RTV_WORK}/*.out
rm ${RTV_WORK}/retrievals_out.txt

cd ${RTV_DATA}
# Choose which data to be used
cp retrievals_10/data_retrievals* .
#
for file in data_retrievals*
do
  echo $file
  rm text_input_file
  ln -s $file text_input_file
  cd ${RTV_WORK}
  ./text_to_obs 
  mv obs_seq.out $file".out" 
  echo $file".out" >> retrievals_out.txt
  cd ${RTV_DATA}
done
rm data_retrievals*
