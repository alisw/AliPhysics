AliAnalysisTask *AddRsnTaskTrain(const char *commonStr,const char *rsnStr,const char *rsnCutStr) {
   // rsnStr -> <Name>
   // rsnCutStr -> <CutName>
   // This will use AddRsnPairs<Name>.C
   // and for cuts AddRsnDaughterCuts<CutName>.C
   // and <opt> string is passed to AddRsnDaughterCuts<CutName>.C
   // so you can control different cut settings
   // string "<Name>:mon" means that it will add monitoring histograms from cuts
   // Note : for now you have to set gRsnUseMiniPackage = 0 to have mon histograms
   //    return AddRsnTask("<Name>:mon","<CutName>:<opt>","");
   // or like we are using it now
   //    return AddRsnTask(rsnStr,rsnCutStr,"");

   Bool_t valid;
   AliRsnTrainManager::GetGlobalStr("LegoTrainPath",valid);
   if (!valid) {
      TString legoTrainPath = "$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train";
      AliRsnTrainManager::SetGlobalStr("LegoTrainPath",legoTrainPath.Data());
   }

   // Creating Rsn Train Manager
   AliRsnTrainManager *rsnMgr = new AliRsnTrainManager();

   if (!RsnLoadMacroTrain("RsnTrainCommonSettings.C")) return kFALSE;
   RsnTrainCommonSettings(commonStr);

   rsnMgr->Print();

   if (!RsnLoadMacroTrain("AddRsnTask.C")) return kFALSE;
   return AddRsnTask(rsnStr,rsnCutStr,"");
}

Bool_t RsnLoadMacroTrain(TString macro,TString path="") {

   Bool_t valid;
   TString lego_path = AliAnalysisManager::GetGlobalStr("RsnLegoTrainPath",valid);
   if (!valid) lego_path = "$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train";

   if (!gSystem->AccessPathName(macro.Data())) {
      gROOT->LoadMacro(macro.Data());
      Printf("Macro loaded from %s/%s ...",gSystem->pwd(),macro.Data());
      return kTRUE;
   }

   if (!gSystem->AccessPathName(gSystem->ExpandPathName(Form("%s/%s",lego_path.Data(),macro.Data())))) {
      gROOT->LoadMacro(gSystem->ExpandPathName(Form("%s/%s",lego_path.Data(),macro.Data())));
      Printf("Macro loaded from %s ...",gSystem->ExpandPathName(Form("%s/%s",lego_path.Data(),macro.Data())));
      return kTRUE;
   }

   Printf("Error loading %s",macro.Data());

   return kFALSE;
}