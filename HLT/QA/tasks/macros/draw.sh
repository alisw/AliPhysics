#!/bin/bash

if [ -d perfImg ] ; then
    rm -rf perfImg
fi

if [ -d perfRoot ] ; then
    rm -rf perfRoot
fi

mkdir perfImg
mkdir -p perfImg/qa
#mkdir -p perfImg/eff
#mkdir -p perfImg/res
mkdir perfRoot

nr1=0
nr2=35

for ii in {0..5} ; do
    if [ ! -d perfImg/$ii ] ; then
	mkdir -p perfImg/qa/$ii
    fi
done
 
if [ ! -d perfImg/event ] ; then
    mkdir -p perfImg/qa/event
fi

if [ $1 ] ; then
    FOLDER=$1
else
    FOLDER="./"
fi

aliroot -b -l -q drawPerformanceTPCQAofflineHLT.C'('\"${FOLDER}\"')'
#aliroot -b -l -q drawPerformanceTPCEff.C'('\"${FOLDER}\"')'
#aliroot -b -l -q drawPerformanceTPCRes.C'('\"${FOLDER}\"')'

for ii in {0..5} ; do
    name=`cat drawPerformanceTPCQAofflineHLT.C | grep "if (cuts == $ii"`
    title=`echo $name | awk -F'// ' '{ printf $2 }' | awk -F' --' '{ printf $1 }'`

    for jj in `ls perfImg/qa/$ii` ; do
	mv perfImg/qa/$ii/$jj perfImg/qa/$ii/${title}_$jj
    done

    mv perfImg/qa/$ii perfImg/qa/$title
done

if [ -f perfImg.tar.bz2 ] ; then
    rm perfImg.tar.bz2
fi

tar cjf perfQA.tar.bz2 perfImg/ perfRoot/

if [ "$FOLDER" != "./" ] ; then
    mv perfImg ${FOLDER}/perfImg
    mv perfRoot ${FOLDER}/perfRoot
fi