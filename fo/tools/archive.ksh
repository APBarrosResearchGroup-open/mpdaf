#!/bin/ksh

odir="DAT_GM_OL"
nens=48

mkdir ${odir}
for iens in {1..${nens}}
do
 echo ${ens}
 rundir=ens_${iens}
 edir=ens_${iens}
 mkdir ${odir}/${edir}
 cp ${rundir}/memls*.dat ${odir}/${edir}
done
