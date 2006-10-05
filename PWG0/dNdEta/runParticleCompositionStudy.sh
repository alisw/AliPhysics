# /* $Id$ */

# This script runs the dN/dEta analysis with different correction maps to gather systematics
# It runs with the "normal map", and with 4 other different cases where particle species are enhanced
# or reduced.
# The normal map is expected in correction_map.root, created by AlidNdEtaCorrectionSelector
# The others in new_compositions.root in the folders (K|p)(Boosted|Reduced), created
#   by AlidNdEtaSystematicsSelector and Composition() out of drawSystematics.C


function run
{
  root -l -q testAnalysis2.C\(\"analysisInputMerged.txt\",10000,0,kFALSE,kFALSE,kTRUE,\"$1\",\"$2\"\)
  mv analysis_esd.root $3

  if [ "$?" -ne "0" ]
  then
    echo "$3 failed"
    exit
  fi

  sleep 5
}

rm analysis_esd.root

## this runs particle composition study
#run correction_map.root dndeta_correction systematics_dndeta_reference.root
#run new_compositions.root KBoosted systematics_dndeta_KBoosted.root
#run new_compositions.root KReduced systematics_dndeta_KReduced.root
#run new_compositions.root pBoosted systematics_dndeta_pBoosted.root
#run new_compositions.root pReduced systematics_dndeta_pReduced.root

## this runs vertex reco study

#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vertexreco_pythia systematics_vtxreco_pythia.root
#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vertexreco_ddmore systematics_vtxreco_ddmore.root
#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vertexreco_ddless systematics_vtxreco_ddless.root
#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vertexreco_sdmore systematics_vtxreco_sdmore.root
#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vertexreco_sdless systematics_vtxreco_sdless.root
#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vertexreco_dmore systematics_vtxreco_dmore.root
#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vertexreco_dless systematics_vtxreco_dless.root


## this runs trigger bias study

#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_trigger_pythia systematics_trigger_pythia.root
#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_trigger_ddmore systematics_trigger_ddmore.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_trigger_ddless systematics_trigger_ddless.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_trigger_sdmore systematics_trigger_sdmore.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_trigger_sdless systematics_trigger_sdless.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_trigger_dmore systematics_trigger_dmore.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_trigger_dless systematics_trigger_dless.root


## this runs trigger bias and vertex reco study

#run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vtxtrigger_pythia systematics_vtxtrigger_pythia.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vtxtrigger_ddmore systematics_vtxtrigger_ddmore.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vtxtrigger_ddless systematics_vtxtrigger_ddless.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vtxtrigger_sdmore systematics_vtxtrigger_sdmore.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vtxtrigger_sdless systematics_vtxtrigger_sdless.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vtxtrigger_dmore systematics_vtxtrigger_dmore.root
run systematics_vtxtrigger_compositions.root dndeta_correction_syst_vtxtrigger_dless systematics_vtxtrigger_dless.root
