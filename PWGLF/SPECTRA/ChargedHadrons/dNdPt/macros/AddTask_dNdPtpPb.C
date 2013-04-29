void AddTask_dNdPtpPb()
{
/*
CheckLoadLibrary("libPWG0base");
CheckLoadLibrary("libPWG0dep");
CheckLoadLibrary("libPWG0selectors");
*/

  //Get current Analysis Manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_dNdPtAnalysis_TPCITS", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
//  AliLog::SetGlobalLogLevel(AliLog::kError);
//  mgr->SetDebugLevel(0);

  //
  // Create event cuts
  //
  Float_t zvWindow = 30. ;

  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kTRUE);
  //evtCuts->SetTriggerRequired(kFALSE);

  //
  // Create geom. acceptance cuts
  //
  //Float_t etaWindow = 1.0;
  eta1 = -0.765409; eta2 = -0.165409; 
  Float_t ptMin = 0.1 ;

  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(eta1,eta2);
  accCuts->SetPtRange(ptMin,1.e10);

  //
  // Create standard esd track cuts
  //
  Int_t cutMode = 2014;

  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");
  //gROOT->LoadMacro("./CreatedNdPtTrackCuts.C");
  AliESDtrackCuts* esdTrackCuts = CreatedNdPtTrackCuts(cutMode);
  if (!esdTrackCuts) {
    printf("ERROR: esdTrackCuts could not be created\n")
    return;
  } else {
    //esdTrackCuts->SetHistogramsOn(kTRUE);
    esdTrackCuts->SetHistogramsOn(kTRUE);
  }

  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //
  // Create task
  //
  AlidNdPtTask *task = new AlidNdPtTask("AlidNdPtTask_TPCITS");
  task->SetUseMCInfo(hasMC);

  // trigger  
  task->SelectCollisionCandidates(AliVEvent::kINT7); 

  //
  // set analysis options from the Helper here !!!
  //

  //AliTriggerAnalysis::Trigger trigger = AliTriggerAnalysis::kMB1;
  AlidNdPtHelper::AnalysisMode analysisMode = AlidNdPtHelper::kTPCITS;
  AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kAllPart ;

  //
  // Create analysis object
  //

  AlidNdPtAnalysispPb *fdNdPtAnalysis = new AlidNdPtAnalysispPb("dNdPtAnalysis_TPCITS","dN/dPt Analysis with TPC-ITS tracking");//Now AnalysispPb before: Analzsiis
  fdNdPtAnalysis->SetEventCuts(evtCuts);
  fdNdPtAnalysis->SetAcceptanceCuts(accCuts);
  fdNdPtAnalysis->SetTrackCuts(esdTrackCuts);
  //fdNdPtAnalysis->SetBackgroundCuts(backCuts);
  fdNdPtAnalysis->SetAnalysisMode(analysisMode); 
  fdNdPtAnalysis->SetParticleMode(particleMode); 
  
  //j kCINT5 should work for 2012 data, kINT7 for 2013
  
  //fdNdPtAnalysis->SetTrigger(trigger);
 // fdNdPtAnalysis->SetTriggerMask(AliVEvent::kCINT5);
  fdNdPtAnalysis->SetTriggerMask(AliVEvent::kINT7);
  //fdNdPtAnalysis->SetTriggerMask(AliVEvent::kEMC1);  
  if(hasMC) 
  {
    //physTrigSel->SetAnalyzeMC();
    //fdNdPtAnalysis->SetPhysicsTriggerSelection(physTrigSel); 

    fdNdPtAnalysis->SetUseMCInfo(kTRUE);
    fdNdPtAnalysis->SetHistogramsOn(kTRUE);
    //fdNdPtAnalysis->SetHistogramsOn(kFALSE);
  }
  else { // online trigger
//     physTrigSel->SetUseBXNumbers();
//     physTrigSel->SetComputeBG();
//     fdNdPtAnalysis->SetPhysicsTriggerSelection(physTrigSel); 
  }
  
    // change binning
    Int_t multNbins = 152;  
    Double_t binsMult[153];
    for (int i=0; i<=multNbins; i++) { binsMult[i] = -0.5 + i; }
    binsMult[152] = 1000.;
    // change binning
    const Int_t ptNbins = 81;
    Double_t bins[82] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0};
    Double_t* binsPt = new Double_t[82];
    for (int i=0; i<82; i++) {binsPt[i] = bins[i];}
    fdNdPtAnalysis->SetBinsPt(ptNbins, binsPt);
    fdNdPtAnalysis->SetBinsPtCorr(ptNbins, binsPt);  
    fdNdPtAnalysis->SetBinsMult(multNbins, binsMult);
    
    
       // y shift -0.465409
    const Int_t EtaNbins = 30;
    Double_t binse[31] = {-1.465409,-1.365409,-1.265409,-1.165409,-1.065409,-0.965409,-0.865409,-0.765409,-0.665409,-0.565409,-0.465409,-0.365409,-0.265409,-0.165409,-0.065409,0.034591,0.134591,0.234591,0.334591,0.434591,0.534591,0.634591,0.734591,0.834591,0.934591,1.034591,1.134591,1.234591,1.334591,1.434591,1.534591};
        Double_t* binsEta = new Double_t[31];
    for (int i=0; i<31; i++) {binsEta[i] = binse[i];}
    fdNdPtAnalysis->SetBinsEta(EtaNbins,binsEta);     
    
    
    // set centrality estimator
   // fdNdPtAnalysis->SetCentralityEstimators("NPA");
  

  // Add analysis object
  task->AddAnalysisObject( fdNdPtAnalysis );
	
  // Add task
  mgr->AddTask(task);

  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
 
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

//  AliAnalysisDataContainer *coutput = mgr->CreateContainer("jgronef_dNdPtpPb_TPCITS", TList::Class(), AliAnalysisManager::kOutputContainer, "jgronef_dNdPtpPb_TPCITS.root");   //    <-- Old Way works on Batch
 
   AliAnalysisDataContainer *coutput = mgr->CreateContainer("dNdPtpPb", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:dNdPtHistos", mgr->GetCommonFileName())); //    <-- New Way changed to work on Grid
  
  
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
    
  
  /*
   //                                                         Aus Phillips Macro kopiert:
  
   // Add task
   mgr->AddTask(task);
   
   // Create containers for input
   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   AliAnalysisDataContainer *coutput  = mgr->CreateContainer("dNdPtPbPb", 
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer, 	
							     Form("%s:dNdPtHistos", mgr->GetCommonFileName()));
							     
							     mgr->ConnectInput(task, 0, cinput);
							     mgr->ConnectOutput(task, 1, coutput);
							     
							    */ 
  
  
  

}

