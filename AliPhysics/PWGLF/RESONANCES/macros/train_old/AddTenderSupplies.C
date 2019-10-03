AliAnalysisTask *AddTenderSupplies
(
   Float_t tofres       = 80,
   Bool_t  corrExpTimes = kFALSE,
   Bool_t  applyT0      = kFALSE
)
{
   // get the current analysis manager
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      Error("AddTask_tender_Tender", "No analysis manager found.");
      exit(0);
      return 0;
   }

   //
   // === Add tender to the ANALYSIS manager and set default storage =====
   //
   AliTender *tender = new AliTender("AnalysisTender");
   tender->SetCheckEventSelection(kFALSE);
   //tender->SetDefaultCDBStorage("raw://");
   tender->SetDefaultCDBStorage("alien://folder=/alice/data/2010/OCDB");
   mgr->AddTask(tender);

   //
   // === Attach VZERO supply ============================================
   //
   AliVZEROTenderSupply *VZEROtender = new AliVZEROTenderSupply("VZEROtender");
   tender->AddSupply(VZEROtender);

   //
   // === Attach TPC supply ==============================================
   //
   AliTPCTenderSupply *TPCtender = new AliTPCTenderSupply("TPCtender");
   tender->AddSupply(TPCtender);

   //
   // === Attach TOF supply ==============================================
   //
   AliTOFTenderSupply *TOFtender = new AliTOFTenderSupply("TOFtender");
   TOFtender->SetTOFres(tofres);
   TOFtender->SetApplyT0(applyT0);
   TOFtender->SetCorrectExpTimes(corrExpTimes);
   tender->AddSupply(TOFtender);

   //
   // === Define output containers, please use 'username'_'somename' =====
   //
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tender_event", AliESDEvent::Class(), AliAnalysisManager::kExchangeContainer, "default_tender");
   mgr->ConnectInput(tender,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(tender,  1, coutput1);

   return tender;
}
