
void AddTask_mkoehler_dNdPtPbPb()
{

	CheckLoadLibrary("libPWG0base");
	CheckLoadLibrary("libPWG0dep");
	CheckLoadLibrary("libPWG0selectors");

 AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (!mgr) {
    Error("AddTask_mkoehler_dNdPtPbPb", "No analysis manager found.");
    return 0;
  }

  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);

  //
  // Create physics trigger selection class
  //
  AliPhysicsSelection *physTrigSel =  new AliPhysicsSelection();

  //
  // Create event cuts
  //
  Float_t zvWindow = 20. ;

  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kTRUE);

  //
  // Create geom. acceptance cuts
  //
  Float_t etaWindow = 1. ;
  Float_t ptMin = 0.15 ;

  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);
  accCuts->SetMaxDCAr(3.0);
  accCuts->SetMaxDCAz(30.0);

  //
  // Create standard esd track cuts
  //
  Int_t cutMode = 23;

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/dNdPt/macros/CreatedNdPtTrackCuts.C");
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
  AlidNdPtTask *task = new AlidNdPtTask("AlidNdPtTask");
  task->SetUseMCInfo(hasMC);

 //
 // set analysis options from the Helper here !!!
 //
 AlidNdPtHelper::OutputObject outputObject = AlidNdPtHelper::kAnalysisPbPb;
 AlidNdPtHelper::AnalysisMode analysisMode = AlidNdPtHelper::kTPC ;
 AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kAllPart ;


  //
  // Create cut analysis object
  //
  if(outputObject==AlidNdPtHelper::kAnalysisPbPb){

    AlidNdPtAnalysisPbPb *fdNdPtAnalysisPbPb = new AlidNdPtAnalysisPbPb("dNdPtAnalysisPbPb","dN/dPt Analysis");
    fdNdPtAnalysisPbPb->SetEventCuts(evtCuts);
    fdNdPtAnalysisPbPb->SetAcceptanceCuts(accCuts);
    fdNdPtAnalysisPbPb->SetTrackCuts(esdTrackCuts);
    fdNdPtAnalysisPbPb->SetAnalysisMode(analysisMode); 
    fdNdPtAnalysisPbPb->SetParticleMode(particleMode); 


    if(hasMC) 
    {
       physTrigSel->SetAnalyzeMC();
       fdNdPtAnalysisPbPb->SetPhysicsTriggerSelection(physTrigSel);

       fdNdPtAnalysisPbPb->SetUseMCInfo(kTRUE);
       fdNdPtAnalysisPbPb->SetHistogramsOn(kTRUE);
       //fdNdPtAnalysisPbPb->SetHistogramsOn(kFALSE);
    }else { // online trigger
    fdNdPtAnalysisPbPb->SetPhysicsTriggerSelection(physTrigSel); 
    }

    task->AddAnalysisObject( fdNdPtAnalysisPbPb );
  }

  // create analysis object
  if(outputObject==AlidNdPtHelper::kAnalysis) 
  {
    AlidNdPtAnalysis *fdNdPtAnalysis = new AlidNdPtAnalysis("dNdPtAnalysis","dN/dPt Analysis");
    fdNdPtAnalysis->SetEventCuts(evtCuts);
    fdNdPtAnalysis->SetAcceptanceCuts(accCuts);
    fdNdPtAnalysis->SetTrackCuts(esdTrackCuts);
    fdNdPtAnalysis->SetAnalysisMode(analysisMode); 
    fdNdPtAnalysis->SetParticleMode(particleMode); 
    if(hasMC) 
    {
       physTrigSel->SetAnalyzeMC();
       fdNdPtAnalysis->SetPhysicsTriggerSelection(physTrigSel);

       fdNdPtAnalysis->SetUseMCInfo(kTRUE);
       fdNdPtAnalysis->SetHistogramsOn(kTRUE);
       //fdNdPtAnalysis->SetHistogramsOn(kFALSE);
    }
    else { // online trigger
    fdNdPtAnalysis->SetPhysicsTriggerSelection(physTrigSel); 
    }

    task->AddAnalysisObject( fdNdPtAnalysis );
  }

 // Add task
  mgr->AddTask(task);

  // Create containers for input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(task, 0, cinput);

  AliAnalysisDataContainer *coutput = mgr->CreateContainer("mkoehler_dNdPtPbPb", TList::Class(), AliAnalysisManager::kOutputContainer, "mkoehler_dNdPtPbPb.root");
  mgr->ConnectOutput(task, 1, coutput);



}

