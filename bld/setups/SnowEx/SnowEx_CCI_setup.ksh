#!/bin/ksh

 defForcingDir="/scsf"
 defStartDate="2016-10-11 00"
 defRunHours=8750

 cd ${scratchdir}
 if [ -d "${scratchdir}/inputdata" ]; then
   echo "Removing Inputdata link ..."
   rm inputdata
 fi
 echo "Linking inputdata ", ${scratchdir}/data/${testcase}
 ln -s ${scratchdir}/data/${testcase}  ${scratchdir}/inputdata

 #cd $rundir
