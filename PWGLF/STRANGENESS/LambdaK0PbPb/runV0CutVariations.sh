#
prod=LHC11a10b_plus

for dir in 0005 0010 1020 2040 4060 6080 8090
do
  s=${dir:0:2}; e=${dir:2:2}
  echo $dir $s $e
  echo "root.exe -q run.C\\($s,$e,kTRUE,kTRUE\\)"  # MC, select non-injected
     root.exe -q runV0CutVariations.C\($s,$e,kTRUE,kTRUE\) >& MC_nonInj.log
  echo "root.exe -q run.C\\($s,$e,kTRUE,kFALSE\\)"  # MC
     root.exe -q runV0CutVariations.C\($s,$e,kTRUE,kFALSE\) >& MC.log
  mkdir $dir/$prod
  mv AliV0CutVariationsMC_nonInj.root AliV0CutVariationsMC.root MC_nonInj.log MC.log $dir/$prod
done
