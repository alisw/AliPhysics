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

runLevelQA()
{
  aliroot -l -b -q ${fwd}/qa/RunQA.C\(\"${qaFile}\",\"${dataType}\",${year},\"${period}\",\"${pass}\",${runNumber}\)
}

periodLevelQA()

{
 
 #path of the merged period trending.root is provided
 
 # 1: file
 # 2: data type 
 # 3: year
 # 4: period
 # 5: pass 
 aliroot -l -b -q ${fwd}/qa/PeriodQA.C\(\"$1\",\"\",0,\"\",\"\"\)

  (cd ${outputDir} && \
      pwd && \
      if test -d FMD ; then cd FMD ; fi ; \
      ${fwd}/qa/makeIndex.sh -t "QA for the FMD" -d "QA" -m 5 \
         -o index.html -l ; \
      cp ${fwd}/qa/style.css . ; \
      cp ${fwd}/qa/script.js . ; \
      cp ${fwd}/qa/fmd_favicon.png . ; \
      cp ${fwd}/qa/fmd_logo.png . ; \
      )
}
#
# EOF
#

