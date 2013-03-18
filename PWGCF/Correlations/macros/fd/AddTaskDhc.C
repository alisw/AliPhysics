// $Id$

// File config location
// alien:///alice/cern.ch/user/t/tschuste/correction_hybrid_nulled.root

AliDhcTask *AddTaskDhc(
  Int_t iAna          = 2, /*1=h-h, 2=mu-h, 4=mu-mu*/ 
  TString chUName     = "", 
  TString chHEffFile  = "", 
  TString chMuEffFile = "", 
  TString chTaskFile  = "", 
  TString chTaskName  = "", 
  TString chNTracks   = "PicoTracks",
  TString centSel     = "V0M",
  TString className   = "",
  Bool_t  doMassCut   = kFALSE,
  Bool_t  doFillSame  = kFALSE,
  UInt_t  trigsel     = AliVEvent::kINT7
)
{
  Char_t chExtraName[256];
  
  // Get the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskDhc", "No analysis manager found.");
    return;
  }

  // ESD or AOD? write into the task name etc.
  sprintf(chExtraName,"_%s",mgr->GetInputEventHandler()->GetDataType());

  AliDhcTask *dhcTask = 0x0;

  // if string chTaskFile is given, load a pre-configured task from file
  if (!chTaskFile.EqualTo("")) {
    iAna=999;
    TFile *fiDhcTask = 0x0;
    fiDhcTask = TFile::Open(chTaskFile,"OLD");
    if (!fiDhcTask){
      Error("AddTaskDhc",Form("Requested file: %s was not opened. ABORT.",chTaskFile));
      return;
    }
    dhcTask = (AliDhcTask*) fiDhcTask->Get(chTaskName);
  }
  else { // create a new task
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
        Error("AddTaskDhc",Form("Requested file: %s was not opened. ABORT.",chHEffFile));
        return;
      }
      hHEff = (THnF*) fiHEff->Get("correction");
    }
    if (!chMuEffFile.EqualTo("")) {
      fiMuEff = TFile::Open(chMuEffFile,"OLD");
      if (!fiMuEff){
        Error("AddTaskDhc",Form("Requested file: %s was not opened. ABORT.",chMuEffFile));
        return;
      }
      hMuEff = (THnF*) fiMuEff->Get("correction");
    }
    
    dhcTask = new AliDhcTask("Task_Dhc_Temp_Name");
    if (iAna==1) { // h-h
      Int_t nDetaBins = 40;
      Int_t nDPhiBins = 72;
      dhcTask->SetAnaMode(AliDhcTask::kHH);
      dhcTask->SetHEffT(hHEff);
      dhcTask->SetHEffA(hHEff);
      dhcTask->SetEtaMax(1.2);
      dhcTask->SetPtTACrit(kTRUE);
      sprintf(chExtraName,"%s_HH",chExtraName);
      if (hHEff) {
        sprintf(chExtraName,"%s_corrH",chExtraName);
      }
    } else if (iAna==2) { // mu-h
      Int_t nDetaBins = 100;
      Int_t nDPhiBins = 36;
      dhcTask->SetAnaMode(AliDhcTask::kMuH);
      dhcTask->SetHEffT(hMuEff);
      dhcTask->SetHEffA(hHEff);
      dhcTask->SetEtaMax(5.0);
      dhcTask->SetPtTACrit(kFALSE);
      sprintf(chExtraName,"%s_MuH",chExtraName);
      if (hMuEff) {
        sprintf(chExtraName,"%s_corrMu",chExtraName);
      }
      if (hHEff) {
        sprintf(chExtraName,"%s_corrH",chExtraName);
      }
    } else if (iAna==4) { // mu-mu
      Int_t nDetaBins = 60;
      Int_t nDPhiBins = 36;
      dhcTask->SetAnaMode(AliDhcTask::kMuMu);
      dhcTask->SetHEffT(hMuEff);
      dhcTask->SetHEffA(hMuEff);
      dhcTask->SetEtaMax(5.0);
      dhcTask->SetPtTACrit(kFALSE);
      sprintf(chExtraName,"%s_MuMu",chExtraName);
      if (hMuEff) {
        sprintf(chExtraName,"%s_corrMu",chExtraName);
      }
    } else {
      Error("AddTaskDhc", Form("iAna %d not known", iAna));
    }
    dhcTask->SetTracksName(chNTracks);
    dhcTask->SetDoWeights(kFALSE);
    dhcTask->SetDoFillSame(doFillSame);
    dhcTask->SetDoMassCut(doMassCut);
    dhcTask->SetClassName(className);
    dhcTask->SetCentMethod(centSel);
    dhcTask->SetDEtaDPhiBins(nDetaBins,nDPhiBins);
    dhcTask->SetPtTBins(axPt);
    dhcTask->SetPtABins(axPt);
    dhcTask->SetCentBins(axCent);
    dhcTask->SetZVtxBins(axZvtx);
    dhcTask->SetCentMixBins(axCentMix);
    dhcTask->SetZVtxMixBins(axZvtxMix);
    dhcTask->SelectCollisionCandidates(trigsel);
    dhcTask->SetVerbosity(0);
  }
  if (!dhcTask) {
    Error("AddTaskDhc","no dhcTask");
    return 0x0;
  }

  // make a unique task name
  Char_t chNewTaskName[256];
  if (chTaskName.EqualTo("")) {
    sprintf(chNewTaskName,"Task_Dhc%s_%s_%s",chExtraName,centSel.Data(),chUName.Data());
  } else {
    sprintf(chNewTaskName,"%s",chTaskName.Data());
  }

  AliDhcTask *mgrTask = mgr->GetTask(chNewTaskName);
  if (mgrTask)
    sprintf(chNewTaskName,"%s_bis%04.0f",chNewTaskName,10000*gRandom->Rndm());
  
  dhcTask->SetName(chNewTaskName);
  Info("AddTaskDhc",Form("DHC Analysis, adding task %s",dhcTask->GetName()));
  mgr->AddTask(dhcTask);
  AliAnalysisDataContainer *co_Dhc = mgr->CreateContainer(Form("Cont_%s",chNewTaskName),
                                                          TList::Class(),
                                                          AliAnalysisManager::kOutputContainer,
                                                          Form("%s:PWGCF_out_%s", AliAnalysisManager::GetCommonFileName(), chNewTaskName));
  mgr->ConnectInput(dhcTask,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(dhcTask,1,co_Dhc);
  
  return dhcTask;
}
