TList *RsnManager() {

   Int_t isPP           = 1; // in GRID case it will be overwriten
   Int_t useRsnMini     = 1;

   Int_t useMixing      = 1;
   Int_t numMix         = 10;

   Int_t fullOutput     = 1;
   Int_t mcMomentum     = 0;

   Int_t useEventMixPar = 0;
   Int_t useRsnPar      = 0;
   
   TString rootver = "v5-30-06";
   TString alirootver = "";
//   alirootver = "v5-02-17-AN";
   
   TString legoTrainPath = "$ALICE_ROOT/PWGLF/RESONANCES/macros/lego_train";
//   legoTrainPath = "$HOME/git/PWG2resonances/PWG2/RESONANCES/macros/lego_train";


   TList *listRsn = new TList();

   // This will use AddRsnPairs<Name>.C
   // and for cuts AddRsnDaughterCuts<CutName>.C
   // and <opt> string is passed to AddRsnDaughterCuts<CutName>.C
   // so you can control different cut settings
   // string "<Name>:mon" means that it will add monitoring histograms from cuts
   // Note : for now you have to set gRsnUseMiniPackage = 0 to have mon histograms
   //    listRsn->Add(new TNamed("<Name>:mon","<CutName>:<opt>"));


   //    listRsn->Add(new TNamed("Phi","Phi2010"));
   listRsn->Add(new TNamed("Phi","Phi2010:mon"));
//   listRsn->Add(new TNamed("Phi","Phi2010:tpconly_sigma1_mon"));
//   listRsn->Add(new TNamed("Phi","Phi2010:tpconly_sigma2_mon"));
//   listRsn->Add(new TNamed("Phi","Phi2010:tpconly_sigma3_mon"));
//   listRsn->Add(new TNamed("Phi","Phi2010:tofonly_sigma1_mon"));
//   listRsn->Add(new TNamed("Phi","Phi2010:tofonly_sigma2_mon"));
//   listRsn->Add(new TNamed("Phi","Phi2010:tofonly_sigma3_mon"));
//   listRsn->Add(new TNamed("Phi","BPID:mon"));
//   listRsn->Add(new TNamed("Phi","Phi2010:qualityonly_mon"));
   //
   //    // in case you have MC
   //    listRsn->Add(new TNamed("Phi","PDG:mon"));
   //
   //    listRsn->Add(new TNamed("KStar","KStar2010:mon"));
   //    listRsn->Add(new TNamed("KStar","BPID:mon"));

   //    listRsn->Add(new TNamed("KStar","KStar:mon"));
   //    listRsn->Add(new TNamed("KStar","KStar:TPCTOFpidDefaultKstarPP2010_mon"));
   //    listRsn->Add(new TNamed("KStar","KStar:FastTPCpid1point5sigma_mon"));
   //    listRsn->Add(new TNamed("KStar","KStar:FastTPCpid2sigma_mon"));


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

   isPP = 1;
   dsConfig = "datasets-grid/LHC11a_AOD072.txt";
   dsConfig = "datasets-grid/LHC11h5b_AOD079.txt";


   //================================================

   ///////////////////////////////////////////
   // don't edit next lines (EXPERTS ONLY)
   ///////////////////////////////////////////

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) { Printf("Error[RsnManager] mgr is null !!!"); return 0; }

   AliAnalysisManager::SetGlobalStr("rsnLegoTrainPath",legoTrainPath.Data());

   AliAnalysisManager::SetGlobalInt("rsnIsPP",isPP);

   // use parfiles instead of libs
   AliAnalysisManager::SetGlobalInt("rsnUseEventMixingPar",useEventMixPar);
   AliAnalysisManager::SetGlobalInt("rsnUseRSNPar",useRsnPar);

   // common options
   AliAnalysisManager::SetGlobalInt("rsnUsePhysSel",0);
   AliAnalysisManager::SetGlobalInt("rsnUsePIDResponse",1);
   // rsn options

   if (useRsnMini) {
      AliAnalysisManager::SetGlobalInt("rsnUseMiniPackage",1);
      AliAnalysisManager::SetGlobalInt("rsnUseRsnInputHandler",0);
      AliAnalysisManager::SetGlobalInt("rsnSplitMgrByTasks",1);
      AliAnalysisManager::SetGlobalInt("rsnUseMixing",useMixing);

   } else  {
      AliAnalysisManager::SetGlobalInt("rsnUseMiniPackage",0);
      AliAnalysisManager::SetGlobalInt("rsnUseRsnInputHandler",1);
      AliAnalysisManager::SetGlobalInt("rsnSplitMgrByTasks",0);
      AliAnalysisManager::SetGlobalInt("rsnUseMixing",0);
   }

   // mixing setting
   AliAnalysisManager::SetGlobalInt("rsnNumMix",numMix);

   // oputput settings
   AliAnalysisManager::SetGlobalInt("rsnOutputFull",fullOutput);
   AliAnalysisManager::SetGlobalInt("rsnUseMCMomentum",mcMomentum);

   // expert options (don't change)
   AliAnalysisManager::SetGlobalInt("rsnMixPrintRefresh",-1);

   // RSN train settings for GRID
   AliAnalysisManager::SetGlobalStr("rsnTrainDSConfig",dsConfig.Data());

   // root and aliroot version
   AliAnalysisManager::SetGlobalStr("rsnLegoTrainROOTversion",rootver.Data());
   AliAnalysisManager::SetGlobalStr("rsnLegoTrainAliROOTversion",alirootver.Data());

   return listRsn;
}
