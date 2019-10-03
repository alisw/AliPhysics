TList *RsnPrivateTrainBaseSettings() {

   Bool_t valid;
   TString legoTrainPath = AliAnalysisManager::GetGlobalStr("RsnLegoTrainPath",valid);
   Int_t usePhysSel        = -1;
//   usePhysSel              = AliVEvent::kMB;

   Int_t usePIDResponseTask      = 1;
   Int_t useCentralityTask       = 0;
   Int_t useEventPlaneTask       = 1;
   Int_t useVZEROEPSelectionTask = 1;
   Int_t usePIDqa                = 0;

   Int_t useEventMixPar    = 0;
   Int_t useRsnPar         = 0;
   Int_t useRsnParDev      = -1;
//   useRsnParDev            = 1;

   TString rootver = "v5-34-02-1";
   TString alirootver = "";
//      alirootver = "v5-03-07-AN";

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
      return;
   }

   // use parfiles instead of libs
   AliAnalysisManager::SetGlobalInt("rsnUseEventMixingPar",useEventMixPar);
   AliAnalysisManager::SetGlobalInt("rsnUseRSNPar",useRsnPar);
   AliAnalysisManager::SetGlobalInt("rsnUseRSNParDev",useRsnParDev);

   // common options
   AliAnalysisManager::SetGlobalInt("rsnUsePIDResponse",usePIDResponseTask);
   AliAnalysisManager::SetGlobalInt("rsnUsePhysSel",usePhysSel);
   AliAnalysisManager::SetGlobalInt("rsnUseCentralityTask",useCentralityTask);
   AliAnalysisManager::SetGlobalInt("rsnUseEventPlaneTask",useEventPlaneTask);
   AliAnalysisManager::SetGlobalInt("rsnUseVZEROEPSelection",useVZEROEPSelectionTask);

   AliAnalysisManager::SetGlobalInt("rsnUsePIDqa",usePIDqa);

   // RSN train settings for GRID
   AliAnalysisManager::SetGlobalStr("rsnTrainDSConfig",dsConfig.Data());
   AliAnalysisManager::SetGlobalInt("rsnGlobalTrainID",globalTrainID);
   AliAnalysisManager::SetGlobalInt("rsnGridNumRuns",numRuns);
   AliAnalysisManager::SetGlobalInt("rsnGridNumRunsSkip",numRunsSkip);

   // root and aliroot version
   AliAnalysisManager::SetGlobalStr("rsnLegoTrainROOTversion",rootver.Data());
   AliAnalysisManager::SetGlobalStr("rsnLegoTrainAliROOTversion",alirootver.Data());

}
