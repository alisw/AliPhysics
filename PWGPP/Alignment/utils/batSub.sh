#!/bin/bash

jdl="algJDL50.jdl"
split='50'
ttl='28000'
pref='out'
ddir='/alice/cern.ch/user/s/shahoian/algTest'

Usage() {
    echo "Usage: mrgfl.sh  inpFile1 inpFile2 ...."
    echo "options:"
    echo "       [-j <JDL=${jdl}>] "
    echo "       [-s <split=${split}>] "
    echo "       [-t <TTL=${ttl}>] "
    echo "       [-p <outPrefix=${pref}>] "
    echo "       [-d <xmldir=${ddir}>]"
    exit 1
}

if [ $# -lt 1 ] ; then Usage ;fi

nxml=0

shift 1 ;

#process argumens
while [[ $# > 0 ]]
do
key="$1"

case $key in
    -j|--jdl)
	shift 1 ;
	if [ $# -lt 1 ] ; then Usage ;fi	
	jdl="$1"
	;;
#
    -s|--split)
	shift 1 ;
	if [ $# -lt 1 ] ; then Usage ;fi	
	split="$1"
	;;
#
    -t|--ttl)
	shift 1 ;
	if [ $# -lt 1 ] ; then Usage ;fi	
	ttl="$1"
	;;
#
    -p|--pref)
	shift 1 ;
	if [ $# -lt 1 ] ; then Usage ;fi	
	pref="$1"
	;;
#
    -d|--dir)
	shift 1 ;
	if [ $# -lt 1 ] ; then Usage ;fi	
	ddir="$1"
	;;	
#
    *)
	flx=$1
	extension="${flx##*.}"
	if [ "$extension" != "xml" ] ; then
	    echo "argument $flx has not .xml extension"
	    Usage
	fi	
	nxml=$((nxml+1))
	xmlList[$nxml]=$1
	;;
esac
shift # past argument or value
done

echo "jdl   = $jdl"
echo "Split = $split"
echo "TTL   = $ttl"
echo "Dir   = $ddir"
echo "$nxml input xml files"

if  [ $nxml -lt 1 ] ; then 
    echo "No xml collections supplied"
    Usage ;
fi

curd='/alice/cern.ch/user/s/shahoian/algTest/v0'

for coll in "${xmlList[@]}"
do
   echo "$coll"
   outDir=$(basename "$coll")
   outDir="${pref}${outDir%.*}"
   echo 'submitting ' ${ddir}/${coll} $curd $outDir $jdl $ttl $split 
   submit $jdl ${ddir}/${coll} $curd $outDir $jdl $ttl $split 
   ps

done
