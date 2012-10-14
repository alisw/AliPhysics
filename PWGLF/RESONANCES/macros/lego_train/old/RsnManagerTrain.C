TList *RsnManagerTrain(TString addRsnManager="AddRsnToManagerTrain.C",
                       Int_t isESD=0,
                       Int_t isMC=0,
                       Int_t isPP=1,
                       Int_t useRsnMini = 1,
                       Int_t useMixing = 0,
                       Int_t numMix = 10,
                       Int_t fullOutput = 1)
{

   // sets Rsn Lego train path
   TString legoTrainPath = "$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train";

   // creates list
   TList *listRsn = new TList();

   TString commonCutOption="";
   commonCutOption = "mon";

   Printf("Adding RsnManger : %s",addRsnManager.Data());
   AddResonanceToRsnManager(listRsn,addRsnManager.Data(),legoTrainPath.Data());

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Printf("Error[RsnManager] mgr is null !!!");
      return 0;
   }

   Bool_t valid;
   AliAnalysisManager::GetGlobalStr("rsnLegoTrainPath",valid);
   if (valid) {
      return list;
   }

   Printf("Setting up RSN variables ...");
   AliAnalysisManager::SetGlobalStr("rsnLegoTrainPath",legoTrainPath.Data());
   AliAnalysisManager::SetGlobalInt("rsnIsPP",isPP);
   AliAnalysisManager::SetGlobalInt("rsnUseMC",isMC);
   AliAnalysisManager::SetGlobalInt("rsnUseMiniPackage",useRsnMini);

   // mixing setting
   AliAnalysisManager::SetGlobalInt("rsnUseMixing",useMixing);
   AliAnalysisManager::SetGlobalInt("rsnNumMix",numMix);

   // oputput settings
   AliAnalysisManager::SetGlobalInt("rsnOutputFull",fullOutput);

   // expert options (don't change)
   AliAnalysisManager::SetGlobalInt("rsnMixPrintRefresh",-1);

   AliAnalysisManager::SetGlobalStr("rsnLegoTrainCommonCutOption",commonCutOption.Data());

   return listRsn;
}

Bool_t AddResonanceToRsnManager(TList *listRsn,TString rsnAddMacro="AddRsnToManagerPhi.C",TString path="") {
   if (!listRsn) return kFALSE;

   RsnManagerLoadMacro(rsnAddMacro,path);
   rsnAddMacro.ReplaceAll(".C","");
   gROOT->ProcessLine(TString::Format("%s((TList*)%p)",rsnAddMacro.Data(),listRsn).Data());

   return kTRUE;
}

Bool_t RsnManagerLoadMacro(TString macro,TString path="") {

   TString lego_path=path;

   if (lego_path.IsNull()) {
      Bool_t valid;
      lego_path = AliAnalysisManager::GetGlobalStr("rsnLegoTrainPath",valid);
      if (!valid) lego_path = "$ALICE_ROOT/PWG2/RESONANCES/macros/lego_train";
   }
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

