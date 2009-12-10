void Plot(Int_t run);
void AliTRDrunCalib(const char*filename,Int_t runnr)
{
 
 
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTRDcalib");
  
 
  AliLog::SetGlobalLogLevel(AliLog::kError);
  gROOT->LoadMacro(Form("%s/PWG0/CreateESDChain.C", gSystem->ExpandPathName("$ALICE_ROOT")));
  TChain *chain = CreateESDChain((const char*)filename, -1);
  chain->SetBranchStatus("*FMD*",0);
  chain->SetBranchStatus("*Calo*",0);
  chain->SetBranchStatus("Tracks", 1);
  chain->SetBranchStatus("ESDfriend*",1);
  chain->Lookup();
  chain->GetListOfFiles()->Print();
  printf("\n ----> CHAIN HAS %d ENTRIES <----\n\n", (Int_t)chain->GetEntries());
    

  TGrid::Connect("alien://");
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("alien://folder=/alice/data/2009/OCDB?cacheFold=/lustre/alice/local/alice/data/2009/OCDB");
  man->SetCacheFlag(kTRUE);
  man->SetRun(runnr);

  AliTRDcalibDB *calib = AliTRDcalibDB::Instance();
  const AliTRDCalDet* caldet = calib->GetGainFactorDet();

  AliAnalysisManager *TRDqa = new AliAnalysisManager("TRD QA/Calibration");
  TRDqa->SetInputEventHandler(new AliESDInputHandler);
  
  AliTRDCalibTask *calibTask = new AliTRDCalibTask();
  calibTask->SetHisto2d(kTRUE);
  calibTask->SetVector2d(kTRUE);
  calibTask->SetVdriftLinear(kTRUE);
  calibTask->SetNz(0,0);
  calibTask->SetNrphi(0,0);
  calibTask->SetNz(0,1);
  calibTask->SetNrphi(0,1);
  calibTask->SetNz(0,2);
  calibTask->SetNrphi(0,2);
  calibTask->SetLow(0);
  calibTask->SetHigh(30);
  calibTask->SetFillZero(kFALSE);
  calibTask->SetDebug(3);
  calibTask->SetNbTimeBins(30);
  //calibTask->SetMaxEvent(20);
  calibTask->SetRequirePrimaryVertex(kTRUE);
  calibTask->SetMinNbOfContributors(1);
  calibTask->SetMaxCluster(100.0);
  calibTask->SetNbMaxCluster(2);
  calibTask->SetCalDetGain(caldet);

  /////////////////////////////
  // Track cuts
  /////////////////////////////
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("trackcuts","trackcuts");
  trackCuts->SetMinNClustersTPC(70);
  trackCuts->SetMaxChi2PerClusterTPC(3.5);
  //trackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  trackCuts->SetRequireTPCRefit(kTRUE);
  //trackCuts->SetRequireITSRefit(kTRUE);
  //trackCuts->SetMinNsigmaToVertex(10);
  //trackCuts->SetRequireSigmaToVertex(kFALSE);
  trackCuts->SetAcceptKinkDaughters(kFALSE);
  trackCuts->SetMaxDCAToVertexZ(10.0);
  //trackCuts->SetMaxDCAToVertexXY(4.0);


  calibTask->SetESDtrackCuts(trackCuts);

  TRDqa->AddTask(calibTask);
  
  TRDqa->ConnectInput( calibTask, 0, TRDqa->GetCommonInputContainer());
  TRDqa->ConnectOutput(calibTask, 0, TRDqa->CreateContainer(calibTask->GetName(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("TRDcalibration%d.root", runnr)));

  if(TRDqa->InitAnalysis()){
    TRDqa->Print();
    TRDqa->StartAnalysis("local", chain);
  }
  
  delete calibTask;

  Plot(runnr);
}

void Plot(Int_t run){
  TFile *calibration=new TFile(Form("TRDcalibration%d.root",run));
  TList *lister = (TList *) calibration->Get("AliTRDCalibTask");
  TH2I *hCH2d = (TH2I*)lister->FindObject("CH2d");
  TProfile2D *hPH2d = (TProfile2D*)lister->FindObject("PH2d"); // tracks/event
  TProfile2D *hPRF2d = (TProfile2D*)lister->FindObject("PRF2d");
  
  calibration->Close();

  TCanvas *tcalib = new TCanvas ("CalibrationMonitor", "CalibrationMonitor",50,50,600,900);
  tcalib->Divide(3,1);
  tcalib->cd(1);
  if(hCH2d) hCH2d->Draw("lego2"); 
  tcalib->cd(2);
  if(hPH2d) hPH2d->Draw("lego2"); 
  tcalib->cd(3);
  if(hPRF2d) hPRF2d->Draw("lego2"); 

  tcalib->SaveAs(Form("CalibrationMon_run%d.gif",run));

}
