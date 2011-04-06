typedef enum AliESDpid::EStartTimeType_t TOFtimeType;

TString     rsnTenderStorage         ("alien://folder=/alice/data/2010/OCDB");
TString     rsnTenderOptions         ("VZERO+TPC+TOF");
TOFtimeType rsnTenderTOFtime         = AliESDpid::kTOF_T0;
Double_t    rsnTenderTOFres          = 80.0; 
Bool_t      rsnTenderTOFcorrExpTimes = kFALSE;
Bool_t      rsnTenderTOFapplyT0      = kFALSE;

//__________________________________________________________________________________________________
//
// Special function to add tender when no multi handler is used
//
void AddTender()
{
   Info("Setup", "Adding tender directly to manager");
   
   AliTender *tender = new AliTender("AnalysisTender");
   tender->SetCheckEventSelection(rsnTenderOptions.Contains("SEL"));
   tender->SetDefaultCDBStorage(rsnTenderStorage.Data());
   AddTenderSupplies(tender);
   
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   mgr->AddTask(tender);
}

//__________________________________________________________________________________________________
//
// Special function to add tender handler when multi handler is used
//
void AddTenderHandler(AliMultiInputEventHandler *multiInputHandler)
{
   Info("Setup", "Adding tender handler");
   if (!multiInputHandler) return;
   
   // works only with ESDs
   AliESDInputHandler*esdIH = dynamic_cast<AliESDInputHandler*>(multiInputHandler->GetFirstInputEventHandler());
   if (!esdIH) {
      Error("Setup", "No ESD handler found");
      return;
   }   
   
   // add tender handler and configure tender inside
   AliTenderInputEventHandler *tenderIH = new AliTenderInputEventHandler();
   AliTender *tender = tenderIH->GetTender();
   AddTenderSupplies(tender);
   
   // add handler to event handler
   multiInputHandler->AddInputEventHandler(tenderIH);
}

//__________________________________________________________________________________________________
//
// Special function to add tender supplies
//
void AddTenderSupplies(AliTender *tender)
{
   if (!tender) return;
   
   Bool_t useV0  = rsnTenderOptions.Contains("V0");
   Bool_t useTPC = rsnTenderOptions.Contains("TPC");
   Bool_t useTOF = rsnTenderOptions.Contains("TOF");
   Bool_t useTRD = rsnTenderOptions.Contains("TRD");
   Bool_t usePID = rsnTenderOptions.Contains("PID");
   Bool_t useVTX = rsnTenderOptions.Contains("PrimVtx");
   Bool_t evSel  = rsnTenderOptions.Contains("SEL");
   
   tender->SetCheckEventSelection(evSel);
   tender->SetDefaultCDBStorage("raw://");

   // VZERO
   if (useV0) {
      Info("Setup", "Adding tender supply for VZERO");
      AliVZEROTenderSupply *vzeroSupply = new AliVZEROTenderSupply("VZEROtender");
      vzeroSupply->SetDebug(kFALSE);
      tender->AddSupply(vzeroSupply);
   }
   
   // TPC
   if (useTPC) {
      Info("Setup", "Adding tender supply for TPC");
      AliTPCTenderSupply *tpcSupply = new AliTPCTenderSupply("TPCtender");
      tpcSupply->SetDebugLevel(2);
      //tpcSupply->SetMip(50.);
      tender->AddSupply(tpcSupply);
   }
   
   // TOF
   if (useTOF) {
      Info("Setup", "Adding tender supply for TOF");
      AliTOFTenderSupply *tofTender = new AliTOFTenderSupply("TOFtender");
      tofTender->SetTimeZeroType(rsnTenderTOFtime);
      tofTender->SetTOFres(rsnTenderTOFres);
      tofTender->SetApplyT0(rsnTenderTOFapplyT0);
      tofTender->SetCorrectExpTimes(rsnTenderTOFcorrExpTimes);
      tender->AddSupply(tofTender);
   }
   
   // TRD
   if (useTRD) {
      Info("Setup", "Adding tender supply for TRD");
      AliTRDTenderSupply *trdSupply = new AliTRDTenderSupply("TRDtender");
      tender->AddSupply(trdSupply);
   }
   
   // PID
   if (usePID) {
      Info("Setup", "Adding tender supply for PID");
      tender->AddSupply(new AliPIDTenderSupply("PIDtender"));
   }
   
   // Primary Vertex
   if (useVTX) {
      Info("Setup", "Adding tender supply for Primary Vertex");
      tender->AddSupply(new AliVtxTenderSupply("PriVtxtender"));
   }
}
