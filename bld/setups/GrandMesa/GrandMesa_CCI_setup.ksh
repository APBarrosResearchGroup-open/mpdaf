#!/bin/ksh

 defForcingDir="/scsf"
 defStartDate="2016-10-01 00"
 defRunHours=5832

 cd ${scratchdir}
 if [ -d "${scratchdir}/inputdata" ]; then
   echo "Removing current inputdata link ..."
   rm inputdata
 fi
 echo "Linking ",${scratchdir}/data/${testcase} 
 ln -s ${scratchdir}/data/${testcase}  ${scratchdir}/inputdata

 #cd $rundir
