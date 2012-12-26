#
centrality="0005 0010 1020 2040 4060 6080 8090"

production="LHC11a10b_plus LHC11a10b_bis"


#### Real data
fn=LHC10h_pass2/Merged.root
for dir in $centrality
do
  s=${dir:0:2}; e=${dir:2:2}
  echo $dir $s $e
  root.exe -q runV0CutVariations.C\($s,$e,kFALSE,kFALSE,\"$fn\"\) >& real.log
  mkdir $dir
  mv AliV0CutVariations.root real.log $dir
done

echo
echo

#### Injected MC productions
for prod in $production
do
    echo "Injected MC production $prod"
    fn=$prod/Merged.root
    for dir in $centrality
    do
        s=${dir:0:2}; e=${dir:2:2}
        echo $dir $s $e
        echo "Non-injected"
      root.exe -q runV0CutVariations.C\($s,$e,kTRUE,kTRUE,\"$fn\"\) >& MC_nonInj.log
        echo "All"
      root.exe -q runV0CutVariations.C\($s,$e,kTRUE,kFALSE,\"$fn\"\) >& MC.log
        mkdir $dir/$prod
        mv AliV0CutVariationsMC_nonInj.root AliV0CutVariationsMC.root MC_nonInj.log MC.log $dir/$prod
    done
done

echo
echo 

#### Merge efficiency input
for dir in $centrality
do
  s=${dir:0:2}; e=${dir:2:2}
  echo $dir $s $e
  cd $dir
  hadd AliV0CutVariationsMC.root `find . -name "AliV0CutVariationsMC.root"`
  hadd AliV0CutVariationsMC_nonInj.root `find . -name "AliV0CutVariationsMC_nonInj.root"`
  cd ..
done

