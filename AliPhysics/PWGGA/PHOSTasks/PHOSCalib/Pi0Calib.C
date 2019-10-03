void Pi0Calib(// const char* dataset="/default/polishch/LHC12c_*_ESDs_pass1_phos_PHOS_calib",
              const char* dataset="/default/polishch/LHC12b_000*_ESDs_pass1_phos",
	      const Int_t nEvents=-1, 
	      const Int_t nEventsSkip=0)
{
  // Enabling authentication with grid certificate
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");

  const TString proofCluster="alice-caf.cern.ch";
  const TString alirootVersion = "VO_ALICE@AliRoot::v5-04-25-AN";

  // Connecting to proof cluster
  TProof::Open(proofCluster.Data());
  //   TProof::Open(proofCluster.Data(), "workers=1x");
  //   TProof::Open(proofCluster.Data(), "masteronly");
  if (!gProof) {
    Error("runAAF.C", "Connection to proof cluster failed!");
    return;
  }

  TList* list = new TList();
  list->Add(new TNamed("ALIROOT_MODE", "ALIROOT"));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS","ANALYSIS:OADB:ANALYSISalice"));
  list->Add(new TNamed("ALIROOT_EXTRA_LIBS","ANALYSIS:OADB:ANALYSISalice:PWGGAPHOSTasks"));

  gProof->EnablePackage(alirootVersion.Data(),list);

  cout << "PHOSCalib: processing dataset " << dataset << endl;
//   AliLog::SetGlobalLogLevel(AliLog::kError);
  // A task can be compiled dynamically with AClic
  gProof->Load("AliCaloPhotonC.cxx++g");
  gProof->Load("AliAnalysisTaskPi0Calib.cxx++g");
  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0Spectrum");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);

  
  // Debug level
  mgr->SetDebugLevel(0);

  // Add physics selection
//  gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
//  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();

  // Add my task
  AliAnalysisTaskPi0Calib *task1 = new AliAnalysisTaskPi0Calib("Pi0Spectrum");
//   task1->SelectCollisionCandidates(AliVEvent::kAnyINT);


/*  
  TFile *fBadMap = TFile::Open("alien:///alice/cern.ch/user/p/prsnko/BadMaps/BadMap_LHC10h_period1.root");
  //  TFile *fBadMap = TFile::Open("BadMap_LHC10h_period1.root");
  if(fBadMap->IsOpen()){
    printf("\n\n...Adding PHOS bad channel map \n") ;
    gROOT->cd();
    char key[55] ;
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key,"PHOS_BadMap_mod%d",mod) ;
      TH2I * h = (TH2I*)fBadMap->Get(key) ;
      if(h)
        task1->SetPHOSBadMap(mod,h) ;
    }
    fBadMap->Close() ;
  }
*/
/*
  TFile *fTOF = TFile::Open("calibTOF.root");
  if(fTOF->IsOpen()){
    printf("\n\n...Adding PHOS TOF recalibration map \n") ;
    gROOT->cd();
    char key[55] ;
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key,"T0_mod%d",mod) ;
      TH2D * h = (TH2D*)fTOF->Get(key) ;
      if(h)
        task1->SetPHOSTOF(mod,h) ;
      sprintf(key,"BadMapTOF_mod%d",mod) ;
      TH2D * h = (TH2D*)fTOF->Get(key) ;
      if(h)
        task1->SetPHOSTOFBadMap(mod,h) ;

    }
    fTOF->Close() ;
  }
*/
  TFile *fAmp = TFile::Open("Calibration_pass1.root");
  if(fAmp->IsOpen()){
    printf("\n\n...Adding PHOS Amp recalibration map \n") ;
    gROOT->cd();
    char key[55] ;
    for(Int_t mod=1;mod<4; mod++){
      sprintf(key,"Mass_mod%d",mod) ;
      TH2D * h = (TH2D*)fAmp->Get(key) ;
      if(h)
        task1->SetPHOSCalib(mod,h) ;
    }
    fAmp->Close() ;
  }

  
  mgr->AddTask(task1);

  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histESD",TList::Class(),AliAnalysisManager::kOutputContainer,"histos.root");
  
  // Connect input/output
  mgr->ConnectInput(task1 , 0, cinput);
  mgr->ConnectOutput(task1, 1, coutput1);
  
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("proof", dataset, nEvents, nEventsSkip);
  }
  
}
