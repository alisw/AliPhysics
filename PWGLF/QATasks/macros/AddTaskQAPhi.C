AliAnalysisTask *AddTaskQAPhi(const char *dataType = "Rsn_pp_ESD_MINI",const char *extra="10.0,-1,0,-1,1,0,10,1,-1,-1,0,1")
{

   /////////////////////////////////////////////////////
   // How to set the parameter for the different systems
   // pp ESD -> "Rsn_pp_ESD_MINI"
   // pp AOD -> "Rsn_pp_AOD_MINI"
   // PbPb ESD -> "Rsn_PbPb_ESD_MINI"
   // PbPb AOD -> "Rsn_PbPb_AOD_MINI"
   // pPb ESD -> "Rsn_pPb_ESD_MINI"
   // pPb AOD -> "Rsn_pPb_AOD_MINI"
   /////////////////////////////////////////////////////

   /////////////////////////////////////////////////////
   //
   // extra - string contains comma separated arguments for macro
   //         $ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train/RsnTrainSettingsExtra.C
   //
   /////////////////////////////////////////////////////

   gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/RESONANCES/macros/lego_train/AddRsnTaskTrain.C"));
   return AddRsnTaskTrain(dataType,"Phi","PhiNsigma:KTPCnsig30","","RsnTrainSettingsExtra.C",extra);
}
