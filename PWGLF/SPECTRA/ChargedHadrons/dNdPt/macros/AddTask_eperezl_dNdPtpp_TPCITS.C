/*
 * 
 * Last modified:  12/02/2016
 * By: Edgar Perez Lezama <eperezle@cern.ch>, GSI-Darmstadt
 *  
 * 
 */




void AddTask_eperezl_dNdPtpp_TPCITS(const int cutMode=223, const char *taskMode = "AllPart kINT7", char *suffix="")
{

 /*
CheckLoadLibrary("libPWG0base");
CheckLoadLibrary("libPWG0dep");
CheckLoadLibrary("libPWG0selectors");
*/

  TString TaskMode(taskMode);
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_dNdPtAnalysis_TPCITS", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);

  //
  // Create physics trigger selection class
  //
  //AliPhysicsSelection *physTrigSel =  new AliPhysicsSelection();
  //physTrigSel->AddBackgroundIdentification(new AliBackgroundSelection());

  //
  // Create event cuts
  //
  Float_t zvWindow =30. ;

  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kTRUE);

  //
  // Create geom. acceptance cuts
  //
  Float_t etaWindow = 1.0 ;
  Float_t ptMin = 0.10;

  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);

  //
  // Create standard esd track cuts
  //
  //   Int_t cutMode = 222;
 gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");
  AliESDtrackCuts* esdTrackCuts = CreatedNdPtTrackCuts(cutMode);
  if (!esdTrackCuts) {
    printf("ERROR: esdTrackCuts could not be created\n");
    return;
  } else {
    esdTrackCuts->SetHistogramsOn(kTRUE);
  }

  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //
  // Create task
  //
  AlidNdPtTask *task = new AlidNdPtTask("AlidNdPtTask_TPCITS");  
  task->SetUseMCInfo(hasMC);

  // trigger  
  if(TaskMode.Contains("MB"))    task->SelectCollisionCandidates(AliVEvent::kMB); 
  if(TaskMode.Contains("INT7"))  task->SelectCollisionCandidates(AliVEvent::kINT7); 
//    task->SelectCollisionCandidates(AliVEvent::kINT7);

  //
  // set analysis options from the Helper here !!!
  //

  //AliTriggerAnalysis::Trigger trigger = AliTriggerAnalysis::kMB1;
  AlidNdPtHelper::AnalysisMode analysisMode = AlidNdPtHelper::kTPCITS;
  
  if(TaskMode.Contains("AllPart"))    AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kAllPart;
  if(TaskMode.Contains("Pion"))       AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kMCPion;
  if(TaskMode.Contains("Proton"))     AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kMCProton;
  if(TaskMode.Contains("Kaon"))       AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kMCKaon;
  if(TaskMode.Contains("Rest"))       AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kMCRest;
  if(TaskMode.Contains("Plus"))       AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kPlus;
  if(TaskMode.Contains("Minus"))      AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kMinus;

  
  //
  // Create analysis object
  //

  AlidNdPtAnalysis *fdNdPtAnalysis = new AlidNdPtAnalysis("dNdPtAnalysis","dN/dPt Analysis with TPC-ITS tracking");
  fdNdPtAnalysis->SetEventCuts(evtCuts);
  fdNdPtAnalysis->SetAcceptanceCuts(accCuts);
  fdNdPtAnalysis->SetTrackCuts(esdTrackCuts);
  fdNdPtAnalysis->SetAnalysisMode(analysisMode); 
  fdNdPtAnalysis->SetParticleMode(particleMode); 
  if(TaskMode.Contains("MB"))     fdNdPtAnalysis->SetTriggerMask(AliVEvent::kMB);
  if(TaskMode.Contains("INT7"))   fdNdPtAnalysis->SetTriggerMask(AliVEvent::kINT7);
//   fdNdPtAnalysis->SetTriggerMask(AliVEvent::kINT7);
  //fdNdPtAnalysis->SetTriggerMask(AliVEvent::kEMC1);  

  fdNdPtAnalysis->SetUseMCInfo(hasMC);
  fdNdPtAnalysis->SetHistogramsOn(hasMC);

    // change binning
    Int_t multNbins = 252;  
    Double_t binsMult[253];
    for (int i=0; i<=multNbins; i++) { binsMult[i] = -0.5 + i; }
    binsMult[252] = 1000.;
    // change binning
    const Int_t ptNbins = 81;
    Double_t bins[82] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0};
    Double_t* binsPt = new Double_t[82];
    for (int i=0; i<82; i++) {binsPt[i] = bins[i];}
    fdNdPtAnalysis->SetBinsPt(ptNbins, binsPt);
    fdNdPtAnalysis->SetBinsPtCorr(ptNbins, binsPt);  
    fdNdPtAnalysis->SetBinsMult(multNbins, binsMult);
    


  // Add analysis object
  task->AddAnalysisObject( fdNdPtAnalysis );
	
  // Add task
  mgr->AddTask(task);

  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

//   AliAnalysisDataContainer *coutput = mgr->CreateContainer("eperezl_dNdPtpp_TPCITS",
//                                                            TList::Class(),
//                                                            AliAnalysisManager::kOutputContainer,
//                                                            Form("%s:dNdPtHistos", mgr->GetCommonFileName())); 

   TString outputFile("AnalysisResults.root");

  AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("dNdPtpp_%d",cutMode),
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   Form("%s:dNdPt_eperezl", outputFile.Data())); 
  
  mgr->ConnectOutput(task, 1, coutput);

}

