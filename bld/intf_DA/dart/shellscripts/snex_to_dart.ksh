#!/bin/ksh

# Convert text data to DART space, data are at daily frequency
# Users need to choose the observation errors, eg., snex 10, snex 20 etc
# Also generates snex_out.txt file which contains the list of files generated 

export SNEX_DIR=$HOME/DART/observations/obs_converters/text_snex/
export SNEX_DATA=${SNEX_DIR}/data
export SNEX_WORK=${SNEX_DIR}/work

#Cleanup
rm ${SNEX_WORK}/*.out
rm ${SNEX_WORK}/snex_out.txt

cd ${SNEX_DATA}
# Choose which data to be used
#cp snex_10/data_snex* .
cp snex_test/data_snex* .
#
for file in data_snex*
do
  echo $file
  rm text_input_file
  ln -s $file text_input_file
  cd ${SNEX_WORK}
  ./text_to_obs 
  mv obs_seq.out $file".out" 
  echo $file".out" >> snex_out.txt
  cd ${SNEX_DATA}
done
rm data_snex*
