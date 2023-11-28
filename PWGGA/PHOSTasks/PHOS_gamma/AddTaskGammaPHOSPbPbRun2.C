AliAnalysisTaskGammaPHOSPbPbRun2* AddTaskGammaPHOSPbPbRun2 (TString name = "PHOSGammaPbPbRun2",
					    TString mode = "2468",
					    UInt_t offlineTriggerMask = AliVEvent::kINT7,
					    Int_t harmonics = 2
					    )
{
  //Add a task AliAnalysisTaskGammaPHOSPbPbRun2 to the analysis train
  //Author: Dmitri Peressounko
  /* $Id$ */

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSGammaFlow", "No analysis manager to connect to");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSGammaFlow", "This task requires an input event handler");
    return NULL;
  }
/*
  if(harmonics!=2 && harmonics!=3){
    ::Error("AddTaskPHOSGammaFlow", Form("Only harmonics 2 and 3 allowed, you use %d",harmonics));
    return NULL;    
  }
  */
  AliAnalysisTaskGammaPHOSPbPbRun2* task = new AliAnalysisTaskGammaPHOSPbPbRun2(Form("%s_%s", name.Data(), mode.Data()));

  mgr->AddTask(task);

  task->SetDistCut(kTRUE) ;
  task->SetHarmonics(harmonics) ;  
  task->SelectCollisionCandidates(offlineTriggerMask);

  //std::vector<Double_t> 
  //mode  = "2468";               // 0) 1% 1) 5% 2) 10% 3) 15% 4) 20% 5) 30% 6) 40% 7) 50% 8) 60% 9) 70% a) 80% b) 90%
				 // examples: 
				 // centrMode = "1", fCenBinEdges = {0, 5, 100}
				 // centrMode = "2468, fCenBinEdges = {0, 10, 20, 40, 60, 100}
				 // centrMode = "468a", fCenBinEdges = {0, 20, 40, 60, 80, 100}
				 // centrMode = "", fCenBinEdges = {0, 100}
  task->SetCentralityIntervals(mode);

//  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer() );
  
  AliAnalysisDataContainer *cinput   = mgr->GetCommonInputContainer(); 
  mgr->ConnectInput(task , 0, cinput);

  TString contName = Form("%s:GammaPHOSPbPb_%s", AliAnalysisManager::GetCommonFileName(), mode.Data());

  AliAnalysisDataContainer *coutput1 =
    mgr->CreateContainer("QA", THashList::Class(), AliAnalysisManager::kOutputContainer, contName);
  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer("Photons", THashList::Class(), AliAnalysisManager::kOutputContainer, contName);
  AliAnalysisDataContainer *coutput3 =
    mgr->CreateContainer("Tracks", THashList::Class(), AliAnalysisManager::kOutputContainer, contName );
  AliAnalysisDataContainer *coutput4 =
    mgr->CreateContainer("EventPlane", THashList::Class(), AliAnalysisManager::kOutputContainer, contName);
  AliAnalysisDataContainer *coutput5 =
    mgr->CreateContainer("MC", THashList::Class(), AliAnalysisManager::kOutputContainer, contName);
  AliAnalysisDataContainer *coutput6 =
    mgr->CreateContainer("Pi0", THashList::Class(), AliAnalysisManager::kOutputContainer, contName);


  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);
  mgr->ConnectOutput(task, 3, coutput3);
  mgr->ConnectOutput(task, 4, coutput4);
  mgr->ConnectOutput(task, 5, coutput5);
  mgr->ConnectOutput(task, 6, coutput6);
  
  return task;
}
