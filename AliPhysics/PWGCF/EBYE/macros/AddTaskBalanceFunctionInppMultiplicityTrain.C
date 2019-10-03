//_________________________________________________________//
AliAnalysisTaskBF *AddTaskBalanceFunctionInppMultiplicityTrain(Int_t nMultMin = 0,
							       Int_t nMultMax = 100,
							       Double_t vertexZ=10.,
							       Double_t DCAxy=-1,
							       Double_t DCAz=-1,
							       Double_t ptMin=0.3,
							       Double_t ptMax=1.5,
							       Double_t etaMin=-0.8,
							       Double_t etaMax=0.8,
							       TString fileNameBase="AnalysisResults") {
  // Creates a balance function analysis task and adds it to the analysis manager.
  // Get the pointer to the existing analysis manager via the static access method.
  TString outputFileName(fileNameBase);
  outputFileName.Append(".root");

  
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //===========================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskBF", "This task requires an input event handler");
    return NULL;
  }
  TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/EBYE/macros/configBalanceFunctionAnalysis.C");
  AliBalance *bf  = 0;  // Balance Function object
  AliBalance *bfs = 0;  // shuffled Balance function object

  if (analysisType=="ESD"){
    bf  = GetBalanceFunctionObject("ESD",0,nMultMin,nMultMax);
    bfs = GetBalanceFunctionObject("ESD",0,nMultMin,nMultMax,kTRUE);
  }
  else if (analysisType=="AOD"){
    bf  = GetBalanceFunctionObject("AOD",0,nMultMin,nMultMax);
    bfs = GetBalanceFunctionObject("AOD",0,nMultMin,nMultMax,kTRUE);
  }
  else{
    bf  = GetBalanceFunctionObject("MC",0,nMultMin,nMultMax);
    bfs = GetBalanceFunctionObject("MC",0,nMultMin,nMultMax,kTRUE);
  }

  // Create the task, add it to manager and configure it.
  //===========================================================================
  AliAnalysisTaskBF *taskBF = new AliAnalysisTaskBF("TaskBF");
  taskBF->SetAnalysisObject(bf);
  taskBF->SetShufflingObject(bfs);
  if(analysisType == "ESD") {
    AliESDtrackCuts *trackCuts = GetTrackCutsObject();
    taskBF->SetAnalysisCutObject(trackCuts);
    
    // offline trigger selection (AliVEvent.h)
    // taskBF->UseOfflineTrigger(); // NOT used (selection is done with the AliAnalysisTaskSE::SelectCollisionCandidates()) 
    // with this only selected events are analyzed (first 2 bins in event QA histogram are the same))
    taskBF->SelectCollisionCandidates(AliVEvent::kMB);
  }
  else if(analysisType == "AOD") {
    // pt and eta cut (pt_min, pt_max, eta_min, eta_max)
    taskBF->SetAODtrackCutBit(128);
    taskBF->SetKinematicsCutsAOD(ptMin,ptMax,etaMin,etaMax);
    //taskBF->SetExtraDCACutsAOD(DCAxy,DCAz);
    taskBF->SetMultiplicityRange(nMultMin,nMultMax);
  }
  
  // vertex cut (x,y,z)
  taskBF->SetVertexDiamond(.3,.3,vertexZ);
  
  //bf->PrintAnalysisSettings();
  mgr->AddTask(taskBF);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //======================================================================
  TString listQAName = "listQA"; listQAName += (Int_t) (nMultMin);
  listQAName += "-"; listQAName += (Int_t) (nMultMax);
  AliAnalysisDataContainer *coutQA = mgr->CreateContainer(listQAName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
 
  TString listBFName = "listBF"; listBFName += (Int_t) (nMultMin);
  listBFName += "-"; listBFName += (Int_t) (nMultMax);
  AliAnalysisDataContainer *coutBF = mgr->CreateContainer(listBFName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  
  TString listBFshuffledName = "listBFshuffled"; listBFshuffledName += (Int_t) (nMultMin);
  listBFshuffledName += "-"; listBFshuffledName += (Int_t) (nMultMax);
  AliAnalysisDataContainer *coutBFS= mgr->CreateContainer(listBFshuffledName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(taskBF, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskBF, 1, coutQA);
  mgr->ConnectOutput(taskBF, 2, coutBF);
  mgr->ConnectOutput(taskBF, 3, coutBFS);

  return taskBF;
}
