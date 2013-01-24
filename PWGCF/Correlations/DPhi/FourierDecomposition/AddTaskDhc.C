AliDhcTask *AddTaskDhc(Int_t iAna = 1, TString chTaskFile = "alien:///alice/cern.ch/user/t/tschuste/LEGO_DhcTask.root")
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
  
  AliDhcTask *dhcTask = 0x0;
  
  if (iAna==99) { // load task from file
    TFile *fiDhcTask = 0x0;
    fiDhcTask = TFile::Open(chTaskFile,"OLD");
    if (!fiDhcTask){
      cout << "Requested file:" << fiDhcTask << " was not opened. ABORT." << endl;
      return;
    }
    dhcTask = (AliDhcTask*) fiDhcTask->Get("Task_tschuste_Dhc");
  } else { // create a new task
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
    TFile *fiHEff  = 0x0;
    TFile *fiMuEff = 0x0;
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
    
    
    dhcTask = new AliDhcTask("Task_tschuste_Dhc");
    if (iAna==1) { // h-h
      Int_t nDetaBins = 40;
      Int_t nDPhiBins = 72;
      dhcTask->SetAnaMode(AliDhcTask::kHH);
      dhcTask->SetHEffT(hHEff);
      dhcTask->SetHEffA(hHEff);
      dhcTask->SetEtaMax(1.2);
    } else if (iAna==2) { // mu-h
      Int_t nDetaBins = 100;
      Int_t nDPhiBins = 36;
      dhcTask->SetAnaMode(AliDhcTask::kMuH);
      dhcTask->SetHEffT(hMuEff);
      dhcTask->SetHEffA(hHEff);
      dhcTask->SetEtaMax(5.0);
    }
    dhcTask->SetTracksName(nTracks);
    dhcTask->SetDoWeights(kFALSE);
    dhcTask->SetCentMethod("V0M");
    dhcTask->SetDEtaDPhiBins(nDetaBins,nDPhiBins);
    dhcTask->SetPtTBins(axPt);
    dhcTask->SetPtABins(axPt);
    dhcTask->SetCentBins(axCent);
    dhcTask->SetZVtxBins(axZvtx);
    dhcTask->SetCentMixBins(axCentMix);
    dhcTask->SetZVtxMixBins(axZvtxMix);
    dhcTask->SelectCollisionCandidates(AliVEvent::kINT7);
    dhcTask->SetVerbosity(10);
    mgr->AddTask(dhcTask);
    
    AliAnalysisDataContainer *co_Dhc = mgr->CreateContainer("Cont_tschuste_DhcAna",
                                                               TList::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               Form("%s:PWGCF.outDhc_%d.root", AliAnalysisManager::GetCommonFileName(), iAna));
    mgr->ConnectInput(dhcTask,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(dhcTask,1,co_Dhc);
  }
  
  return dhcTask;
  
}
