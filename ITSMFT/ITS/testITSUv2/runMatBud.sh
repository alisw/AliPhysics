#!/bin/bash

# This simple procedure produces all geometry_*root files
# needed for an analysis of the material budget using the
# macro GetMaterialBudget.C and MakeMatBudPlots.C
# It takes one optional parameter, the maximum  number of build levels.
# It only requires CreateITSUv2.C_template as input
# M.Sitta 11 Nov 2014, 05 May 2015

if [ "x$1" != "x" ]; then
  MAXBUILDLEVEL=$1
else
  MAXBUILDLEVEL=6
fi

rm -Rf geometry_[0-5].root

for i in `seq 0 $MAXBUILDLEVEL`; do
  sed "s/BUILDLEVEL/$i/" CreateITSUv2.C_template >CreateITSUv2.C
  aliroot -b -q sim.C\(1\)
  mv geometry.root geometry_$i.root
done

ls -l geometry_*root
