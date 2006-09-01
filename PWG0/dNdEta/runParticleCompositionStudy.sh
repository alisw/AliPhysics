/* $Id$ */

# This script runs the dN/dEta analysis with different correction maps to gather systematics
# It runs with the "normal map", and with 4 other different cases where particle species are enhanced
# or reduced.
# The normal map is expected in correction_map.root, created by AlidNdEtaCorrectionSelector
# The others in new_compositions.root in the folders (K|p)(Boosted|Reduced), created
#   by AlidNdEtaSystematicsSelector and Composition() out of drawSystematics.C


function run
{
  root -l -q testAnalysis2.C\(\"analysisInput.txt\",10000,0,kFALSE,kFALSE,kTRUE,\"$1\",\"$2\"\)
  mv analysis_esd.root $3

  if [ "$?" -ne "0" ]
  then
    echo "$3 failed"
    exit
  fi

  sleep 5
}

rm analysis_esd.root

#run correction_map.root dndeta_correction systematics_dndeta_reference.root
#run new_compositions.root KBoosted systematics_dndeta_KBoosted.root
#run new_compositions.root KReduced systematics_dndeta_KReduced.root
run new_compositions.root pBoosted systematics_dndeta_pBoosted.root
run new_compositions.root pReduced systematics_dndeta_pReduced.root

