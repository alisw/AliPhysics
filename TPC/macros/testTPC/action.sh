#!/bin/bash
#aguments
#1 TString jobID, 
#2 TString inputData
#3 TString outputDir
#4 TString   action

echo INPUT ARGUMENT FOR ACTION
echo $1
echo $2
echo $3
echo $4
#
# 
#
olddir=`pwd`
mkdir -p $jobhome/$1
mkdir -p $jobhome/$1/$4
cd       $jobhome/$1/$4
ln -sf $olddir home
echo HOME DIR `pwd`
#
#
#
echo GETTING DATA - START
date

CISIZE=0
if [  $AGENTINPUTTYPE -eq 2 ]; then 
    CINAME=`echo $2| sed s_root://voalice04.cern.ch:1094/__`
    echo CASTOR COPY $CINAME ; 
    CISIZE=`nsls -l $CINAME | gawk '{print $5}'`
    if [  $CISIZE -gt 1000000 ]; then 
        rfcp  $CINAME data.root 
    else
       echo FILE TOO SMALL
       exit
    fi;
fi;
if [  $AGENTINPUTTYPE -eq 0 ]; then 
    echo XRD COPY $2 ;
    xrdcp -np $2 data.root
fi;
echo GETTING DATA - STOP
date
echo LS DATA
ls -al
 
CISIZE=`ls -l data.root | gawk '{print $5}'`
if [  $CISIZE -lt 100000 ]; then 
    exit
fi;
#
#
#
echo BEGIN ACTION
date
echo aliroot -b -q $olddir/macros/$4.C
aliroot -b -q $olddir/macros/$4.C
echo END ACTION $1
date
#
#
#
rm data.root
rm TPCsignal.root
ls -al `pwd`
echo CREATING ZIP FILE
zip -n root $4 *.root *.log


isxrd=`echo $3 | grep root://`
mkdir $olddir/filelists
flist=$olddir/filelists/$1.list
touch $flist

#
# ALIEN - PROBLEM WITH HOST CERTIFICATE
#
if [  ${#alien_HOME} -gt 1 ]; then
    dirname=`echo $alien_HOME$3/$4 | sed s/.root//g  `
    echo alien_mkdir $dirname
    alien_mkdir -p $dirname
    echo ALIEN COPY DATA START 
    echo  alien_cp -n $4.zip alien:${dirname}/$4_se.zip@ALICE::GSI::SE
    #echo  alien_cp -n  $4.zip alien:${dirname}$4.zip@ALICE::ALICE::CASTOR2 
     for name in `ls *.root`; do
	echo alien_cp -n $name alien:${dirname}/$name@ALICE::GSI::SE 
	alien_cp -n $name alien:${dirname}/$name@ALICE::GSI::SE
	
    done
    alien_cp  -n  $4.zip  alien:${dirname}/$4_se.zip@ALICE::GSI::SE  
    #alien_cp  $4.zip -n alien:${dirname}$4.zip@ALICE::ALICE::CASTOR2  
    echo END OF ALIEN COPY
fi;



if [  ${#isxrd} -lt 1 ]; then
   dirname=`echo $CASTOR_HOME$3/$4 | sed s/.root//g  `
   echo CASTOR COPPING DATA - START
   echo DIRNAME  - $dirname
   echo XRD PATH = root://voalice04.cern.ch:1094/$dirname
   echo     
   rfmkdir -p  $dirname
   #
   for name in `ls *.root`; do
	echo  rfcp  $name $dirname/$name
	rfcp $name $dirname/$name
	echo root://voalice04.cern.ch:1094/$dirname/$name >>$flist
    done
    echo rfcp $4.zip $dirname/$4.zip
    rfcp $4.zip $dirname/$4.zip
    echo root://voalice04.cern.ch:1094/$dirname/$4.zip >>$flist
    echo END OF COPY    
else
    echo XRD COPY DATA START
    dirname=`echo $3/$4 | sed s/.root//g` 
    echo DIRNAME  - $dirname
    for name in `ls *.root`; do
	echo  xrdcp $name $dirname/$name
	xrdcp -np $name $dirname/$name
	echo $dirname/$name >>$flist
    done
    echo xrdcp -np  $4.zip  $dirname.zip
    xrdcp -np $4.zip  $dirname/$4.zip 
    echo $dirname/$4.zip >>$flist
    echo END OF COPY
fi;



cd $olddir

