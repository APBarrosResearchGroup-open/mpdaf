#!/bin/ksh

#------------------------------------------------------------------------

if [[ $# -eq 0 ]];then
   print "Usage: ./build.ksh <0,1>"
   exit
fi

getRoot(){
  #automatically determine root dir
  cpwd=`pwd`
  rootdir=`echo $cpwd | sed 's@^/@@'`
  rootdir=`echo $cpwd | sed 's@/bld@@'` #remove bld to get rootpath
}
getRoot
print "Root Directory is ......................." $rootdir

tmp_dir=$rootdir"/tmp_mshm"

if [[ $1 == 0 ]] then
  print "Clean build ......."
  mkdir -p $tmp_dir
  cp $rootdir/mshm/mshm.f90 $tmp_dir
  cd $tmp_dir
else
  cd $tmp_dir
fi

gfortran mshm.f90
mv $tmp_dir/a.out $rootdir/bin

print "build script finished sucessfully"

