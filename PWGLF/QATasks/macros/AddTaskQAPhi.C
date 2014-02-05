AliAnalysisTask *AddTaskQAPhi(const char *dataType = "Rsn_pp_ESD_MINI")
{
   gROOT->LoadMacro(gSystem->ExpandPathName("$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train/AddRsnTaskTrain.C"));
   return AddRsnTaskTrain(dataType,"Phi","PhiNsigma:KTPCnsig30","","RsnTrainSettingsExtra.C","10.0,-1,0,-1,1,0,10,1,-1,-1,0,1");
}
