AliAnalysisTaskSE *AddTask_mkrueger_CutTest(Int_t cutMode = 223){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_mkrueger_CutTest", "No analysis manager found.");
    return 0;
  }

    //
  // Create event cuts
  //
  Float_t zvWindow = 10.;

  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kTRUE);

  //
  // Create geom. acceptance cuts
  //
  Float_t etaWindow = 0.8;
  Float_t ptMin = 0.0;
  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);

  string cutModeString = std::to_string(cutMode);
  TMacro createTrackCutsMacro(gSystem->ExpandPathName("$ALICE_PHYSICS/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C"));
  AliESDtrackCuts* esdTrackCuts = reinterpret_cast<AliESDtrackCuts *>(createTrackCutsMacro.Exec(cutModeString.c_str()));
  if(!esdTrackCuts) cout << "ERROR: Could not find esdTrackCuts." << endl;

   esdTrackCuts->SetHistogramsOn(kFALSE);

    // check for mc
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

  //
  // Create task
  //
//  gROOT->LoadMacro("cuttest/AliAnalysisTaskCutTest.cxx+g");
  AliAnalysisTaskCutTest *task = new AliAnalysisTaskCutTest("cuttestANDresolution");

  task->SetUseMCInfo(hasMC);
  //task->SelectCollisionCandidates(AliVEvent::kMB & AliVEvent::kCINT5 & AliVEvent::kINT7);
//  task->SelectCollisionCandidates(AliVEvent::kMB);
  task->SelectCollisionCandidates(AliVEvent::kINT7);
  task->SetTrackCuts(esdTrackCuts);
  task->SetAcceptanceCuts(accCuts);
  task->SetEventCuts(evtCuts);

  mgr->AddTask(task);


  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  AliAnalysisDataContainer *coutput1 =
  mgr->CreateContainer("mkrueger_CutTest", TObjArray::Class(),
			     AliAnalysisManager::kOutputContainer,"AnalysisResults.root");

  //           connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput1);

  return task;
}
