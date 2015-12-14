#!/bin/bash

#REMEMBER TO CHANGE INPUT HANDLERS IN RUN MACROS FOR AODs!
# ./runTrig.sh -d '/alice/data/2015/LHC15f/<RUN09>/pass1/%.%/;runList=225106,226225,226220,226170'
# ./runTrig.sh -d '/alice/data/2015/LHC15f/225106,226225/pass1/%.%/'
# defaults
PROOFDATASET="Find;BasePath=/alice/data/2015/LHC15f/000226062/pass2/%.%/;FileName=root_archive.zip;Anchor=AliESDs.root;Tree=/esdTree;Mode=remote;"
isMC="kFALSE"
NEV=2000
FIRSTEV=0
#ADDTASKMACRO="AddAnalysisTaskdNdEtaPP13(\\\"outfile.root\\\",\\\"%s\\\",%f,%f,%f,%f,\\\"%s\\\",%f,%f,%d,%d,%d,%d,%f,%f,%d,%f,%f,%f,%f,%d,%f,%f,%d)"
ADDTASKMACRO="AddAnalysisTaskdNdEtaPP13(<outname>,\\\"%s\\\",%f,%f,%f,%f,\\\"%s\\\",%f,%f,%d,%d,%d,%d,%f,%f,%d,%f,%f,%f,%f,%d,%f,%f,%d,%u)"
#outname is replaced dynamically below, chance the code if you want to change the filename policy
# is mc is replaced below
USEPHYSICSSELECTION=kTRUE
OUTFNAME="AnalysisResults.root"
DOMULTSELTREE=kFALSE

# GRID STUFF
DATADIR="/alice/sim/2015/LHC15g3c2/"
RUNLIST="226062"
DATAPATTERN="/*/AliESDs.root"
GRIDWORKINGDIR="dNdeta_LHC15g3c2_MC"

RUNPROOF="NO"
RUNLOCAL="NO"
RUNTESTDEST="NO"
RUNGRID="NO"
RUNCORRECTION="NO"

OADBMCMULT="OADB-LHC15g3c2_plus.root"
GRIDMODE="full"

CENTRVAR="V0M"

UNIQUENAMECORR="dndeta"

give_help() {

    echo "Available Options: "
    echo " -h                    Give Help"
    echo " -d <dataset>          Set the dataset to be processed on proof"
    echo "                       The format on VAF id the following: "
    echo "                        \"/alice/data/2015/LHC15f/000225768/pass1/%.%/;\""
    echo "                       if you want to specify and AOD processing you need to edit the runProof macro"
    echo "                       (Default: $PROOFDATASET)"
    echo " -l <files>            Set files to be processed locally"
    echo " -r proof|local|dset|grid "
    echo "                       Run on proof, locally or on grid. The dset mode simply checks the dataset creation on VAF"
    echo "                       for grid, you also have to specify the additional options below"
    echo " -n <nev>              Set number of events (default: $NEV)"
    echo " -u                    Show FIXME and TODO from classes and macros"
    echo " -s kTRUE|kFALSE       Use physics selection (default: $USEPHYSICSSELECTION)"
    echo " -m                    Use this flag if processing MC"
    echo " -t                    Enable multiplicity selection trees (not yet supported on proof, flag needs to be propagated)"
    echo " -o fname              Set out filename (default $OUTFNAME)"
    echo " -a oadbMCMult         Set the OADB to be used for the mult selection in the MC. Can either be a period name "
    echo "                       (e.g. LHC15f), to use the same absolute boundaries used in data, or a root file "
    echo "                       (e.g. OADB-LHC15g3c2_plus.root) which contains the calibration histos to be shipped with "
    echo "                       the analysis task. The root file can be obtained with the macro"
    echo "                       OADB/COMMON/MULTIPLICITY/macros/CalibratePeriodMC.C"
    echo "                       (Default: $OADBMCMULT)"
    echo " -v centravar          Centrality variable (Default $CENTRVAR)"
    echo ""
    echo "Other GRID Options"
    echo " -b <basedir>          Set the base data dir for GRID processing (Default: $DATADIR)"
    echo " -p <pattern>          Specify the data pattern, e.g. \"/muon_calo_pass1/*/AliESDs.root\""
    echo "                       Default: $DATAPATTERN"
    echo " -w <workdir>          Set The workdir (Default $GRIDWORKINGDIR)"
    echo " -x <RUNLIST>          Runlist, comma separated (Default $RUNLIST)"
    echo " -g <gridmode>         Typically, full or test (Default: $GRIDMODE)"
    echo ""
    echo "Correction:"
    echo " -c filedata.root filemc.root    Runs the correction macro"
    echo " -o <tag>                        Unique name prepended to corrected histos output (Default $UNIQUENAMECORR)"
    echo ""
    echo "VAF Cheatsheet"
    echo "vaf-enter"
    echo "vafctl start"
    echo "vafreq <NUM_OF_WORKERS>"
    echo "vafcount"
    echo "vafctl stop #release resources"

}

while getopts "hd:l:mr:g:n:s:o:p:w:x:b:ca:tv:" opt; do
  case $opt in
      h) 
	  give_help
	  exit 0
	  ;;
      d) 
	  PROOFDATASET="Find;BasePath=$OPTARG;FileName=root_archive.zip;Anchor=AliESDs.root;Tree=/esdTree;Mode=remote;"
	  ;;
      l)
	  rm  files.txt
	  for i in $OPTARG
	  do 
	      echo $i >> files.txt
	  done
	  ;;
      m)
	  isMC=kTRUE;
	  ;;
      v)
	  CENTRVAR=$OPTARG;
	  ;;
      n) 
	  NEV=$OPTARG
	  ;;
      a) 
	  OADBMCMULT=$OPTARG
	  ;;
      o) 
	  OUTFNAME=$OPTARG
          UNIQUENAMECORR=$OPTARG
	  ;;
      p) 
	  PATTERN=$OPTARG
	  ;;
      w) 
	  GRIDWORKINGDIR=$OPTARG
	  ;;
      b) 
	  DATADIR=$OPTARG
	  ;;
      x) 
	  RUNLIST=$OPTARG
	  ;;
      s) 
	  USEPHYSICSSELECTION=$OPTARG
	  ;;
      t)
          DOMULTSELTREE=kTRUE
          ;;
      g)
          GRIDMODE=$OPTARG
          ;;

      r)
	  if [ "$OPTARG" = "proof" ]
	      then
	      RUNPROOF=YES
	  fi
	  if [ "$OPTARG" = "local" ]
	      then
	      RUNLOCAL=YES
	  fi
	  if [ "$OPTARG" = "dset" ]
	      then
	      RUNTESTDSET=YES
          fi
	  if [ "$OPTARG" = "grid" ]
	      then
	      RUNGRID=YES
	  fi
	  ;;
      u) 
	  for i in `ls *.{cxx,h,C}`
	  do
	      # the first grep is just used for the return code (stored in $?) so that we only pring relevant files
	      grep 'FIXME\|TODO' $i > /dev/null 
	      if [ "$?" -eq "0" ]
		  then
		  echo -e "\033[01;31m===== $i =====\033[0m"
		  grep --color=auto -i -n -B1 -A5 'FIXME\|TODO' $i
		  echo -e "\033[01;31m========================================\033[0m"
	      fi
	  done
	  exit 1
	  ;;
      c)
          RUNCORRECTION=YES
          ;;

      \?)
	  echo "Invalid option: -$OPTARG" >&2
	  exit 1
	  ;;
      :)
	  echo "Option -$OPTARG requires an argument." >&2
	  exit 1
	  ;;

  esac
done

if [ $OPTIND -eq 1 ]; then give_help; fi

if [ "$RUNPROOF" = "YES" ]
    then   
        #    ADDTASKMACRO=${ADDTASKMACRO/<outname>/\\\"out_`basename $PROOFDATASET`.root\\\"}
    ADDTASKMACRO=${ADDTASKMACRO/<outname>/\\\"$OUTFNAME\\\"}
    root -b -q runProofdNdeta.C\(\"$PROOFDATASET\",$USEPHYSICSSELECTION,$isMC,$NEV,$FIRSTEV,\"$OADBMCMULT\",\"$ADDTASKMACRO\"\)
fi

if [ "$RUNLOCAL" = "YES" ]
    then   
        echo "RUNLOCAL Not supported"
fi

if [ "$RUNGRID" = "YES" ]
    then   
    ADDTASKMACRO=${ADDTASKMACRO/<outname>/\\\"$OUTFNAME\\\"}
    root -b -q runGridEsd.C\(\"$DATADIR\",\"$RUNLIST\",\"$DATAPATTERN\",\"$GRIDWORKINGDIR\",$USEPHYSICSSELECTION,$isMC,$DOMULTSELTREE,\"$OADBMCMULT\",\"$GRIDMODE\",\"$CENTRVAR\",\"$ADDTASKMACRO\"\)
fi

if [ "$RUNTESTDSET" = "YES" ]
    then   
    root -b -q CreateVAFDataset.C\(\"$PROOFDATASET\"\)
fi


if [ "$RUNCORRECTION" = "YES" ]
then
    shift $((OPTIND -1)) # allows remaining parameters to be parsed as $1 ...
    root ../macros/ SaveCanvas.C+ CorrectSpectraMultiMCBG.C+\(\"$1\",\"$2\",\"clist\",\"$UNIQUENAMECORR\"\)
fi
  
