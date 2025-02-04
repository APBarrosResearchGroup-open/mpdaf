#!/bin/ksh

# User Settings
#------------------------------------------------------------------------
# Model Namelist Parameters
#------------------------------------------------------------------------
inputdatadir="/home/$USER/scratch/inputdata/"
ibegin=0                              #Start time step 
ntsp=17035                            #Total time step
nlevs=30                              #Number of vertical levels
ngridpts=1                            #Number of Grid Pts
dt=1800                               #Time step in seconds 
ptb1=0.4                              #Multiplicative random uniform Precipitation perturbation
ptb2=2.0                              #Additive random uniform Temperature perturbation
ptb3=0.05                             #Multiplicative random uniform SW perturbation
ptb4=0.1                              #Multiplicative random uniform LW perturbation
ptb5=0.25                             #Additive random uniform Kscale perturbation
kfrac=1.0                             #Melt fraction of top layer for rain on snow
fdate="2016 10 11 00 00 00"           #Start day of run (yyyy mm dd hh mm ss)
platform="CCI"                        #Select Machine
testcase="SnowEx"                     #Select Testcase
numEns=1                              #Select ensemble size
scratchdir="/scratch/${USER}/"  #Select rundir
#------------------------------------------------------------------------
#Do not change below
#------------------------------------------------------------------------

getRoot(){
  #automatically determine root dir
  cpwd=`pwd`
  rootdir=`echo $cpwd | sed 's@^/@@'`
  rootdir=`echo $cpwd | sed 's@/bld@@'` #remove bld to get rootpath
}
getRoot
print "Root Directory is ......................." $rootdir

#Check date
date=`date +%d%m%y-%H%M%S`

#Select rundir
rundir="${scratchdir}/mpdaf_${platform}_run"

mkdir -p $rundir
if [[ -n $(ls -A $rundir) ]] ; then
     mv ${rundir} ${rundir}_backup_${date}
     mkdir $rundir
fi

#Check for ensemble run
origrundir=${rundir}
for inst in {1..${numEns}}
do
  if [[ ${numEns} > 1 ]] ; then
     rundir=${origrundir}/mpdaf_instance_${inst}  
     mkdir -p ${rundir}
  fi
  #Copy executables
  cp $rootdir/bin/a.out $rundir
  #Copy namelist
  cp $rootdir/bld/setups/${testcase}/namelist.nml ${rundir}
  sed "s,__inputdatadir__,${inputdatadir}," -i ${rundir}/namelist.nml 
  sed "s,__ibegin__,${ibegin}," -i ${rundir}/namelist.nml
  sed "s,__ntsp__,${ntsp}," -i ${rundir}/namelist.nml
  sed "s,__inl__,${nlevs}," -i ${rundir}/namelist.nml
  sed "s,__inp__,${ngridpts}," -i ${rundir}/namelist.nml
  sed "s,__rdt__,${dt}," -i ${rundir}/namelist.nml
  sed "s,__fdate__,${fdate}," -i ${rundir}/namelist.nml
  sed "s,__ptb1__,${ptb1}," -i ${rundir}/namelist.nml
  sed "s,__ptb2__,${ptb2}," -i ${rundir}/namelist.nml
  sed "s,__ptb3__,${ptb3}," -i ${rundir}/namelist.nml
  sed "s,__ptb4__,${ptb4}," -i ${rundir}/namelist.nml
  sed "s,__ptb5__,${ptb5}," -i ${rundir}/namelist.nml
  sed "s,__kfrac__,${kfrac}," -i ${rundir}/namelist.nml
  sed "s,__eid__,${inst}," -i ${rundir}/namelist.nml
  sed "s,__nens__,${numEns}," -i ${rundir}/namelist.nml

  #Setup test case
  . ${rootdir}/bld/setups/${testcase}/${testcase}_${platform}_setup.ksh
done

#Setup machine
  . $rootdir/bld/machines/${platform}/build_intf_${platform}.ksh


print "setup script finished sucessfully"
print "Run Directory is ......................." $origrundir
print "Number of ensembles is ........" ${numEns}
ls -l $rundir
