TList *RsnManager() {

   Int_t isESD             = 0;
   
   Int_t isPP              = 1;
   Int_t useRsnMini        = 1;

   Int_t useMixing         = 0;
   Int_t numMix            = 10;

   Int_t fullOutput        = 1;
   Int_t mcMomentum        = 0;
   Int_t mcMon             = 0;
   
   Int_t usePhysSel        = 0;
   Int_t useCentralityTask = 0;

   Int_t useEventMixPar    = 0;
   Int_t useRsnPar         = 0;
   Int_t useRsnParDev      = -1;

   TString rootver = "v5-32-01";
   TString alirootver = "";
//      alirootver = "v5-03-07-AN";

   TString legoTrainPath = "$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train";
//    legoTrainPath = "$HOME/git/AliRsn/PWGLF/RESONANCES/macros/lego_train";


   TList *listRsn = new TList();

   // This will use AddRsnPairs<Name>.C
   // and for cuts AddRsnDaughterCuts<CutName>.C
   // and <opt> string is passed to AddRsnDaughterCuts<CutName>.C
   // so you can control different cut settings
   // string "<Name>:mon" means that it will add monitoring histograms from cuts
   // Note : for now you have to set gRsnUseMiniPackage = 0 to have mon histograms
   //    listRsn->Add(new TNamed("<Name>:mon","<CutName>:<opt>"));

   TString commonCutOption="";
   commonCutOption = "mon";
//    commonCutOption += "_eta";

   AddResonanceToRsnManager(listRsn,"AddRsnToManagerPhi.C",legoTrainPath.Data());
//    AddResonanceToRsnManager(listRsn,"AddRsnToManagerKStar.C",legoTrainPath.Data());
//    AddResonanceToRsnManager(listRsn,"AddRsnToManagerRho.C",legoTrainPath.Data());
//    AddResonanceToRsnManager(listRsn,"AddRsnToManagerLambda.C",legoTrainPath.Data());


   
   //============= ONLY for GRID ====================
   TString dsConfig;

   //   isPP = 0;
   //   dsConfig = "datasets-grid/LHC10h_p2_ESD.txt";
   //   dsConfig = "datasets-grid/LHC10h_p2_AOD049.txt";
   //   dsConfig = "datasets-grid/LHC10h_p2_AOD073.txt";

   //   dsConfig = "datasets-grid/LHC11a10b_AOD080.txt";

   //      isPP = 1;
   //      dsConfig = "datasets-grid/LHC10b_p2_ESD.txt";
   //      dsConfig = "datasets-grid/LHC10b_p2_AOD038.txt";

   // pp 2.76 TeV data
   isPP = 1;
   // data
   dsConfig = "datasets-grid/LHC11a_AOD072.txt";
   // mc
   dsConfig = "datasets-grid/LHC11h5b_AOD079.txt";
   dsConfig = "datasets-grid/LHC11e3a_AOD074.txt";



   //================================================

   ///////////////////////////////////////////
   // don't edit next lines (EXPERTS ONLY)
   ///////////////////////////////////////////
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Printf("Error[RsnManager] mgr is null !!!");
      return 0;
   }

   AliAnalysisManager::SetGlobalStr("rsnLegoTrainPath",legoTrainPath.Data());

   AliAnalysisManager::SetGlobalInt("rsnIsPP",isPP);

   // use parfiles instead of libs
   AliAnalysisManager::SetGlobalInt("rsnUseEventMixingPar",useEventMixPar);
   AliAnalysisManager::SetGlobalInt("rsnUseRSNPar",useRsnPar);
   AliAnalysisManager::SetGlobalInt("rsnUseRSNParDev",useRsnParDev);
   // common options
   AliAnalysisManager::SetGlobalInt("rsnUsePhysSel",usePhysSel);
   AliAnalysisManager::SetGlobalInt("rsnUseCentralityTask",useCentralityTask);
   AliAnalysisManager::SetGlobalInt("rsnUsePIDResponse",1);
   // rsn options

   AliAnalysisManager::SetGlobalInt("rsnUseMixing",useMixing);
   if (useRsnMini) {
      AliAnalysisManager::SetGlobalInt("rsnUseMiniPackage",1);
      AliAnalysisManager::SetGlobalInt("rsnUseRsnInputHandler",0);
      AliAnalysisManager::SetGlobalInt("rsnSplitMgrByTasks",1);
   } else  {
      AliAnalysisManager::SetGlobalInt("rsnUseMiniPackage",0);
      AliAnalysisManager::SetGlobalInt("rsnUseRsnInputHandler",1);
      AliAnalysisManager::SetGlobalInt("rsnSplitMgrByTasks",0);
//       AliAnalysisManager::SetGlobalInt("rsnUseMixing",0);
   }

   // mixing setting
   AliAnalysisManager::SetGlobalInt("rsnNumMix",numMix);

   // oputput settings
   AliAnalysisManager::SetGlobalInt("rsnOutputFull",fullOutput);
   AliAnalysisManager::SetGlobalInt("rsnUseMCMomentum",mcMomentum);
   AliAnalysisManager::SetGlobalInt("rsnUseMCMonitoring",mcMon);


   // expert options (don't change)
   AliAnalysisManager::SetGlobalInt("rsnMixPrintRefresh",-1);

   // RSN train settings for GRID
   AliAnalysisManager::SetGlobalStr("rsnTrainDSConfig",dsConfig.Data());

   // root and aliroot version
   AliAnalysisManager::SetGlobalStr("rsnLegoTrainROOTversion",rootver.Data());
   AliAnalysisManager::SetGlobalStr("rsnLegoTrainAliROOTversion",alirootver.Data());

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

