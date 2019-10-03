# QA processing for the FMD
#
# it defines two functions: runLevelQA and periodLevelQA
# 
#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber    e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
fwd=${ALICE_PHYSICS}/PWGLF/FORWARD/analysis2
# fwd=$ANA_SRC

runLevelQA()
{
  #full path of QAresults.root is provided
  ${fwd}/qa/runQA.sh "$1" "$dataType" "$year" "$period" "$pass" "$runNumber"
}

periodLevelQA()
{
  #path of the merged period trending.root is provided
  ${fwd}/qa/periodQA.sh "$1" 

  if test x$outputDir = x ; then 
      outputDir=/afs/cern.ch/work/a/aliqafmd/www
  fi
  if test ! -d $outputDir ; then 
      echo "No output directory specified, index not recreated"
      exit 0
  fi
  (cd ${outputDir} && \
      pwd && \
      # if test -d FMD ; then cd FMD ; fi ; \
      ${fwd}/makeIndex.sh \
	  -t "QA for the FMD" \
	  -d "For more see <a href='://cern.ch/go/6Bwz'>TWiki pages</a>." \
	  -m 5 \
	  -o index.html \
	  -l \
	  -f ; \
      cp ${fwd}/qa/style.css . ; \
      cp ${fwd}/qa/script.js . ; \
      cp ${fwd}/qa/fmd_favicon.png . ; \
      cp ${fwd}/qa/fmd_logo.png . ; \
      )
}
#
# EOF
#

