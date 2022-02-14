#!/bin/bash

compileOnly=0
verbose=0

optList="cv"
while getopts $optList option
do
case $option in
  c ) compileOnly=1;;
  v ) verbose=1;;
  * ) echo "Unimplemented option chosen."
  EXIT=1
  ;;
  esac
done

shift $(($OPTIND - 1))

# needs 3 arguments
if [[ ( $compileOnly -eq 0 && $# -ne 3 ) || "$EXIT" -eq 1 ]]; then
  echo "Usage: $(basename $0) (-v) <processDir> <numevts> <seed>"
  echo "       $(basename $0) (-v) -c [ <processDir> ]"
  echo "       -c compile only"
  echo "       -v verbose"
  exit 4
fi


processDir=${1-"all"}
numevts=${2-0}
seed=${3-0}

if [ $verbose -eq 1 ]; then
  echo "Software info"
  uname -a
  if [ `which lsb_release` ]; then
    lsb_release -a
  fi
  gcc --version
  gfortran --version
fi

lhapdfDir="LHAPDF"
lhapdfConfig="lhapdf-config"
lhapdfTar="${lhapdfDir}.tgz"
pdfSetDir="PDFsets"
pdfSetTar="${pdfSetDir}.tgz"
powhegDir="POWHEG-BOX"
powhegTar="${powhegDir}.tgz"

executable="pwhg_main"
inputFile="powheg.input"
baseInputFile="base_$inputFile"


if [ $verbose -eq 1 ]; then
  echo "Path variables:"
  echo $PATH
  echo $LD_LIBRARY_PATH
fi

############################################
# Strip cvmfs from path (last resource...) #
############################################
function stripPath()
{
  echo $1 | awk -F ":" '{
    stripped="";
    for(ipath=0; ipath<NF;ipath++) {
      if ( index($ipath,"cvmfs") > 0 ) continue;
      if ( length(stripped)>0 ) stripped=stripped ":";
      stripped=stripped $ipath;
    } print stripped;
  }'
}

orig_PATH="$PATH"
orig_LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
#####################################################

currDir=`pwd`

if [ $compileOnly -eq 0 ]; then
  if [ $numevts -le 0 ]; then
    echo "Number of events mst be >0"
    exit
  fi
  if [ ! -e "$baseInputFile" ]; then
    echo "Cannot find ${baseInputFile}. Nothing done";
    exit
  fi
fi

if [ ! -e "$powhegTar" ]; then
  echo "Cannot find ${powhegTar}. Nothing done";
  exit
fi
tar -xzvf "$powhegTar"


if [ -e "$pdfSetTar" ]; then
  tar -xzvf $pdfSetTar
  export LHAPATH="$pdfSetDir"
elif [ -e "$lhapdfTar" ]; then
  if [ ! -z "$ALICE_ROOT" ]; then
    echo "Use PDF in aliroot"
    export LHAPATH="$ALICE_ROOT/LHAPDF/PDFsets"
  else
    echo "Please provide $pdfSetTar with PDFsets"
    echo "This is needed for LHAPDF!"
    exit
  fi
elif [ -e "$lhapdfConfig" ]; then
  export LHAPATH="$ALICE_ROOT/LHAPDF/PDFsets"
else
  echo "Assume you're using the native pdf"
  pdfFilename="cteq6m"
  cp -pi $powhegDir/pdfdata/${pdfFilename}.tbl $pdfFilename
  useLhapdf=0
fi

if [ -e "$lhapdfTar" ]; then

  export PATH=$(stripPath $PATH)
  export LD_LIBRARY_PATH=$(stripPath $LD_LIBRARY_PATH)

  if [ $verbose -eq 1 ]; then
    echo "Changing path (will come back to original one at end of script)"
    echo $PATH
    echo $LD_LIBRARY_PATH
  fi

  lhapdfInstallDir="$currDir/install_LHAPDF"
  tar -xzvf $lhapdfTar
  if [ "$(find $lhapdfDir -name libLHAPDF.*)" ]; then
    mv $lhapdfDir $lhapdfInstallDir
  else
    cd $lhapdfDir
    make clean
    ./configure --prefix="$lhapdfInstallDir" --disable-octave --disable-old-ccwrap --disable-doxygen --disable-octave --enable-static=no
    make install
    cd $currDir
#    cp $lhapdfDir/config.log . #REMEMBER TO CUT
  fi
  export PATH="$PATH:$lhapdfInstallDir/bin"
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$lhapdfInstallDir/lib"

elif [ -e "$lhapdfConfig" ]; then
  # Used to link with LHAPDF in aliroot
  echo "Linking with LHAPDF in aliroot"
  chmod u+x $lhapdfConfig
  export PATH="$currDir:$PATH"
fi


processPath="$powhegDir/$processDir"
if [ "$processDir" = "all" ]; then
  processPath="$powhegDir"
fi

makefileList=$(find "$processPath" -name "Makefile")
for ifile in $makefileList; do
  execFile=${ifile//"Makefile"/"$executable"}
  if [ ! -e "$execFile" ]; then
    cd $(dirname $ifile)
    make $executable
  fi
  cd $currDir
done

if [ $compileOnly -eq 0 ]; then
  awk -v lseed="$seed" -v lnumevts="$numevts" '{
    if ( index($0,"iseed") ) gsub("12345",lseed);
    if ( index($0,"numevts") ) gsub("12345",lnumevts);
    print $0; }' $baseInputFile > $inputFile

  $processPath/pwhg_main
fi

export PATH="$orig_PATH"
export LD_LIBRARY_PATH="$orig_LD_LIBRARY_PATH"
