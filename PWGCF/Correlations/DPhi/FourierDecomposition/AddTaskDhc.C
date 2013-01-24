AliDhcTask *AddTaskDhc(Int_t iAna = 1)
{
  const char *nTracks     = "PicoTracks";
  const char *inputTracks = "HybridTracks";

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskDhc", "No analysis manager found.");
    return;
  }

  // Track Cuts
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTpcTrack.C");
  AliEmcalEsdTpcTrackTask *hybTask = AddTaskEmcalEsdTpcTrack(inputTracks,"Hybrid_LHC11h");

  // Pico Tracks
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(nTracks, inputTracks, "LHC11h");

  // Binning
  Double_t arPt[5] = {0.5, 1.0, 2.0, 4.0};
  TAxis *axPt = new TAxis(3,arPt);
  Double_t arCent[5] = {0.0, 20.0, 40.0, 60.0, 100.0};
  TAxis *axCent = new TAxis(4,arCent);
  TAxis *axZvtx = new TAxis(1,-10.0,10.0);
  Double_t arCentMix[9] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0};
  TAxis *axCentMix = new TAxis(8,arCentMix);
  TAxis *axZvtxMix = new TAxis(8,-10.0,10.0);

  // Efficiency correction files
  TFile *fiHEff  = 0;
  TFile *fiMuEff = 0;
  fiHEff = TFile::Open("alien:///alice/cern.ch/user/t/tschuste/correction_hybrid_nulled.root","OLD");
  if (!fiHEff){
    cout << "Requested file:" << fiHEff << " was not opened. ABORT." << endl;
    return;
  }
  THnF* hHEff = (THnF*) fiHEff->Get("correction");
  
  fiMuEff = TFile::Open("alien:///alice/cern.ch/user/t/tschuste/correction_muon.root","OLD");
  if (!fiMuEff){
    cout << "Requested file:" << fiMuEff << " was not opened. ABORT." << endl;
    return;
  }
  THnF* hMuEff = (THnF*) fiMuEff->Get("correction");
  

  AliDhcTask *dhcTask_NW = new AliDhcTask("Task_tschuste_Dhc_NW");
  if (iAna==1) { // h-h
    Int_t nDetaBins = 40;
    Int_t nDPhiBins = 72;
    dhcTask_NW->SetAnaMode(AliDhcTask::kHH);
    dhcTask_NW->SetHEffT(hHEff);
    dhcTask_NW->SetHEffA(hHEff);
    dhcTask_NW->SetEtaMax(1.2);
  } else if (iAna==2) { // mu-h
    Int_t nDetaBins = 100;
    Int_t nDPhiBins = 36;
    dhcTask_NW->SetAnaMode(AliDhcTask::kMuH);
    dhcTask_NW->SetHEffT(hMuEff);
    dhcTask_NW->SetHEffA(hHEff);
    dhcTask_NW->SetEtaMax(5.0);
  }
  dhcTask_NW->SetTracksName(nTracks);
  dhcTask_NW->SetDoWeights(kFALSE);
  dhcTask_NW->SetCentMethod("V0M");
  dhcTask_NW->SetDEtaDPhiBins(nDetaBins,nDPhiBins);
  dhcTask_NW->SetPtTBins(axPt);
  dhcTask_NW->SetPtABins(axPt);
  dhcTask_NW->SetCentBins(axCent);
  dhcTask_NW->SetZVtxBins(axZvtx);
  dhcTask_NW->SetCentMixBins(axCentMix);
  dhcTask_NW->SetZVtxMixBins(axZvtxMix);
  dhcTask_NW->SelectCollisionCandidates(AliVEvent::kINT7);
  dhcTask_NW->SetVerbosity(10);
  mgr->AddTask(dhcTask_NW);
  
  AliAnalysisDataContainer *co_Dhc_NW = mgr->CreateContainer("Cont_tschuste_DhcAna_NW", 
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             Form("%s:PWGCF.outDhc_NW_%d.root", AliAnalysisManager::GetCommonFileName(), iAna));
  mgr->ConnectInput(dhcTask_NW,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(dhcTask_NW,1,co_Dhc_NW);
  
  return dhcTask_NW;

}
