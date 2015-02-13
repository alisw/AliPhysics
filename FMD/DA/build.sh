#!/bin/sh

tgt=$1 ; shift
if test "x$tgt" = "x" ;then 
    tgt=FMDPedestalda
fi


dateCFlags=`/opt/date/.commonScripts/date-config --cflags | sed "s/long' 'long/long\\ long/"`
dateLibs=`/opt/date/.commonScripts/date-config --monitorlibs`
amoreCFlags=`/opt/amore/bin/amore-config --includes`
amoreLibs=`/opt/amore/bin/amore-config --auxlibs`
rootCFlags=`root-config --cflags`
rootLibs=`root-config --libs`
dateCFlags=`echo $dateCFlags | sed 's/-Dlong64="long long"//'`
aliSrc=$HOME/source
aliBuild=$HOME/build
aliLibs="
    -I${aliSrc}/FMD/FMDrec -L${aliBuild}/FMD/FMDrec -lFMDrec \
    -I${aliSrc}/FMD/FMDsim -L${aliBuild}/FMD/FMDsim -lFMDsim \
    -I${aliSrc}/FMD/FMDutil -L${aliBuild}/FMD/FMDutil -lFMDutil \
    -I${aliSrc}/FMD/FMDbase -L${aliBuild}/FMD/FMDbase -lFMDbase \
    -I${aliSrc}/ANALYSIS/ANALYSIS -L${aliBuild}/ANALYSIS/ANALYSIS -lANALYSIS \
    -I${aliSrc}/RAW/RAWDatarecOnline  -L${aliBuild}/RAW/RAWDatarecOnline  -lRAWDatarecOnline  \
    -I${aliSrc}/RAW/RAWDatarec  -L${aliBuild}/RAW/RAWDatarec  -lRAWDatarec  \
    -I${aliSrc}/RAW/RAWDatasim  -L${aliBuild}/RAW/RAWDatasim  -lRAWDatasim  \
    -I${aliSrc}/RAW/RAWDatabase -L${aliBuild}/RAW/RAWDatabase -lRAWDatabase \
    -I${aliSrc}/HLT/BASE/HOMER -L${aliBuild}/HLT/BASE/HOMER -lAliHLTHOMER \
    -I${aliSrc}/HLTb/BASE -L${aliBuild}/HLT/BASE -lHLTbase \
    -I${aliSrc}/STEER/STEER -L${aliBuild}/STEER/STEER -lSTEER \
    -I${aliSrc}/STEER/CDB  -L${aliBuild}/STEER/CDB  -lCDB  \
    -I${aliSrc}/STEER/ESD -L${aliBuild}/STEER/ESD -lESD \
    -I${aliSrc}/STEER/STEERBase -L${aliBuild}/STEER/STEERBase -lSTEERBase \
"
echo "dateCFlags=$dateCFlags"
set -x
g++ -DALI_AMORE \
    ${dateCFlags} \
    -Dlong64="long long" \
    ${amoreCFlags} \
    ${rootCFlags} \
    -I/opt/daqDA-lib/ -L /opt/daqDA-lib \
    -I$ALICE_ROOT/include -L$ALICE_ROOT/lib  \
    ${aliLibs} \
    -lRMySQL \
    ${rootLibs} \
    $@ \
    ${tgt}.cxx \
    ${amoreLibs} \
    ${dateLibs} \
    -ldaqDA \
    -o ${tgt}_dyn

