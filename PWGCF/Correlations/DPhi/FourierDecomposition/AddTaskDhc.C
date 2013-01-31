AliDhcTask *AddTaskDhc(Int_t iAna = 2, TString chHEffFile = "alien:///alice/cern.ch/user/t/tschuste/correction_hybrid_nulled.root", TString chMuEffFile = "", TString chTaskFile = "", TString chTaskName = "")
{
  const char *nTracks     = "PicoTracks";
  const char *inputTracks = "HybridTracks";
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskDhc", "No analysis manager found.");
    return;
  }

  TString chIsESD("ESD");
  if (chIsESD.EqualTo(mgr->GetInputEventHandler()->GetDataType())) {
    Info("AddTaskDhc","adding ESD track selection tasks ...");
    // ESD Track Cuts
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTpcTrack.C");
    AliEmcalEsdTpcTrackTask *hybTask = AddTaskEmcalEsdTpcTrack(inputTracks,"Hybrid_LHC11h",kFALSE);
    // Pico Tracks
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(nTracks, inputTracks, "LHC11h");
  } else {
    Info("AddTaskDhc","AOD analysis, no extra track selection tasks required.");
    sprintf(nTracks,"tracks");
  }  

  AliDhcTask *dhcTask = 0x0;
  
  if (!chTaskFile.EqualTo("")) { // if string is given, load the task from file
    iAna=999;
    TFile *fiDhcTask = 0x0;
    fiDhcTask = TFile::Open(chTaskFile,"OLD");
    if (!fiDhcTask){
      cout << "Requested file:" << fiDhcTask << " was not opened. ABORT." << endl;
      return;
    }
    dhcTask = (AliDhcTask*) fiDhcTask->Get(chTaskName);
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
    THnF* hHEff    = 0x0;
    THnF* hMuEff   = 0x0;
    if (!chHEffFile.EqualTo("")) {
      fiHEff = TFile::Open(chHEffFile,"OLD");
      if (!fiHEff){
        cout << "Requested file:" << fiHEff << " was not opened. ABORT." << endl;
        return;
      }
      hHEff = (THnF*) fiHEff->Get("correction");
    }
    
    if (!chMuEffFile.EqualTo("")) {
      fiMuEff = TFile::Open(chMuEffFile,"OLD");
      if (!fiMuEff){
        cout << "Requested file:" << fiMuEff << " was not opened. ABORT." << endl;
        return;
      }
      hMuEff = (THnF*) fiMuEff->Get("correction");
    }
    
    dhcTask = new AliDhcTask(Form("Task_tschuste_Dhc_%d",iAna));
    if (iAna==1) { // h-h
      Int_t nDetaBins = 40;
      Int_t nDPhiBins = 72;
      dhcTask->SetAnaMode(AliDhcTask::kHH);
      dhcTask->SetHEffT(hHEff);
      dhcTask->SetHEffA(hHEff);
      dhcTask->SetEtaMax(1.2);
      dhcTask->SetPtTACrit(kTRUE);
    } else if (iAna==2) { // mu-h
      Int_t nDetaBins = 100;
      Int_t nDPhiBins = 36;
      dhcTask->SetAnaMode(AliDhcTask::kMuH);
      dhcTask->SetHEffT(hMuEff);
      dhcTask->SetHEffA(hHEff);
      dhcTask->SetEtaMax(5.0);
      dhcTask->SetPtTACrit(kFALSE);
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
    dhcTask->SetVerbosity(0);
  }
  if (!dhcTask) {
    cout << "no DhcTask---return" << endl;
    return 0x0;
  }

  mgr->AddTask(dhcTask);    
  AliAnalysisDataContainer *co_Dhc = mgr->CreateContainer(Form("Cont_tschuste_DhcAna_%d",iAna),
                                                          TList::Class(),
                                                          AliAnalysisManager::kOutputContainer,
                                                          Form("%s:PWGCF_outDhc_%d", AliAnalysisManager::GetCommonFileName(), iAna));
  mgr->ConnectInput(dhcTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(dhcTask,1,co_Dhc);
  
  return dhcTask;
  
}
