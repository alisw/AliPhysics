#! /bin/bash	
# 1. Merge all the QAData files from a given Run
# 2. Run rge QA Checker on the merged file
PRODYEAR=$2
PRODCYCLE=$3
RUNNUMBER=$4
BASEDIR=""
MODE=""
XMLCOLL="QA"_$RUNNUMBER".xml"
PATTERN="Merged.QA.Data.root"
MERGEDOUT="Merged.QA.Data.$RUNNUMBER.root"
dialog=/sw/bin/dialog
Help(){
echo "Usage: MergeAndCheck.sh -h|s|r|c [xxxx, yyyy, zzzz]"
echo "                        -h : get this list"
echo "                        -s <production year> <production cycle> <run number> : processes MC data"
echo "                        -r <year> <production cycle> <run number> : processes RAW data"
echo "                        -c <run number> : processes RAW data only"
}

SetupAliEnEnv(){
echo "...............Setting up AliEn environment"
if [ ! -e "$HOME/.globus/usercert.pem" ]; then
  echo "No grid certificate found in $HOME/.globus"
  exit 1
fi  
if [ ! -e "/tmp/gclient_env_$UID" ]; then
  alien-token-init
  if [ ! "$?" -eq "0" ]; then
    echo "alien-token-init failed"
    exit 1
  fi
else
  source /tmp/gclient_env_$UID
fi
}

MakeXML(){
echo "...............Make XML for" $BASEDIR
gbbox find -x MergedQACollection $BASEDIR $PATTERN > $XMLCOLL
if [ ! "$?" -eq "0" ]; then
  echo "MergeXML failed"
  exit 1
fi
}

Merge(){
echo "...............Make Merge of" $BASEDIR
aliroot -b << EOF
  AliQAManager * qam = AliQAManager::QAManager($MODE) ; 
  qam->MergeXML("$XMLCOLL", NULL, "$MERGEDOUT")
 .q
EOF
if [ ! "$?" -eq "0" ]; then
  echo "Merge failed"
  exit 1
fi
}

Check(){
echo "...............Make Check of $MERGEDOUT"
if [ -e "QA.$RUNNUMBER.root" ] ; then 
  rm QA.$RUNNUMBER.root
fi  
if [ -e "QA.root" ] ; then 
  rm QA.root
fi
aliroot -b << EOF
  AliQAv1::SetQARefStorage("local://$ALICE_ROOT/QAref") ;
  AliQAChecker::Instance()->Run("$MERGEDOUT")
 .q
EOF
mv QA.root QA.$RUNNUMBER.root
if [ ! "$?" -eq "0" ]; then
  echo "Check failed"
  exit 1
fi
}

Save(){
echo "...............Save $MERGEDOUT and QA.$RUNNUMBER.root in alien://$BASEDIR"
gbbox cp file:$MERGEDOUT $BASEDIR
gbbox cp file:QA$.RUNNUMBER.root $BASEDIR 
if [ ! "$?" -eq "0" ]; then
  echo "Save failed"
  exit 1
fi
}

# Start
if [ $# -eq 0 ] ; then
 Help 
 exit 1
fi
# help
if [ "$1" == "-h" ] ; then
 Help
 exit 0
fi
# run
if [ "$1" == "-s" -o "$1" == "-r" -o "$1" == -a ] ; then
  if [ $# -lt 4 ] ; then
    echo "Missing data !"
    Help
    exit 1
  fi
  if [ "$1" == "-s" ] ; then 
    BASEDIR="/alice/sim/"$PRODYEAR"/"$PRODCYCLE"/"$RUNNUMBER 
    MODE="AliQAv1::kSIMMODE"
   fi
  if [ "$1" == "-r" ] ; then 
    BASEDIR="/alice/data/"$PRODYEAR"/"$PRODCYCLE"/"$RUNNUMBER 
    MODE="AliQAv1::kRECMODE"
  fi
  SetupAliEnEnv 
  MakeXML
  Merge
  Check
  Save
fi
if [ "$1" == "-c" ] ; then 
  if [ $# -lt 2 ] ; then
    echo "Missing data !"
    Help
    exit 1
  fi
  RUNNUMBER=$2
  MERGEDOUT="Merged.QA.Data.$RUNNUMBER.root"
  Check
fi

# Done 
if [ -e $XMLCOLL ] ; then 
  rm $XMLCOLL
fi
exit 0