TList *RsnManager() {

   Int_t isESD             = 0;
   
   Int_t isPP              = 0;
   Int_t useRsnMini        = 1;

   Int_t useMixing         = 1;
   Int_t numMix            = 5;

   Int_t fullOutput        = 1;
   Int_t mcMomentum        = 0;
   Int_t mcMon             = 0;
   
   Int_t usePhysSel        = -1;
//   usePhysSel              = AliVEvent::kMB;
   Int_t useCentralityTask = 0;
   Int_t useEventPlaneTask = 0;
   
   Double_t eventCutVertex = 10.0;
   Int_t useYNotEta        = 0;
   
   // useCommonQualityCut=-1  -> Defaultcuts for 2010
   Int_t useCommonQualityCut = -1;
   useCommonQualityCut = 5;

   Int_t useEventMixPar    = 0;
   Int_t useRsnPar         = 0;
   Int_t useRsnParDev      = -1;

   TString rootver = "v5-33-02b";
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
   TString dsConfig = "datasets-grid/LHC11e3a_AOD074.txt";
   Int_t globalTrainID=0;
   Int_t numRuns = 1000;
   Int_t numRunsSkip = 0;

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
   AliAnalysisManager::SetGlobalInt("rsnUseEventPlaneTask",useEventPlaneTask);
   
   AliAnalysisManager::SetGlobalInt("rsnUsePIDResponse",1);
   
   AliAnalysisManager::SetGlobalInt("rsnCommonQualityCut",useCommonQualityCut);
   AliAnalysisManager::SetGlobalInt("rsnUseRapidity",useYNotEta);
   
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
   AliAnalysisManager::SetGlobalDbl("rsnEventCutVertex",eventCutVertex);

   // expert options (don't change)
   AliAnalysisManager::SetGlobalInt("rsnMixPrintRefresh",-1);

   // RSN train settings for GRID
   AliAnalysisManager::SetGlobalStr("rsnTrainDSConfig",dsConfig.Data());
   AliAnalysisManager::SetGlobalInt("rsnGlobalTrainID",globalTrainID);
   AliAnalysisManager::SetGlobalInt("rsnGridNumRuns",numRuns);
   AliAnalysisManager::SetGlobalInt("rsnGridNumRunsSkip",numRunsSkip);

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
      if (!valid) lego_path = "$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train";
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

