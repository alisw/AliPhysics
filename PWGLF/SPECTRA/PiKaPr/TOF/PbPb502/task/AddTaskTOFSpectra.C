AliAnalysisTaskTOFSpectra *AddTaskTOFSpectra(Int_t optTree = kTRUE, Bool_t readMC = kFALSE, Bool_t HeavyIon = kTRUE, Bool_t ChannelMismatch = kTRUE, Bool_t CutVariation = kFALSE, Int_t SimpleCutMode = -1, TString prefix = "", TString tname = "TOFSpectra"){
  // Creates, configures and attaches to the train the task for pi, K , p spectra
  // Get the pointer to the existing analysis manager via the static access method.
  //============================================================================== 

  
  Info("AddTaskTOFSpectra","Adding a new task %s with this settings optTree = %i, readMC = %i, HeavyIon = %i, ChannelMismatch = %i, CutVariation = %i, SimpleCutMode = %i", tname.Data(), optTree,readMC, HeavyIon, ChannelMismatch, CutVariation, SimpleCutMode);  
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTOFSpectra", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskTOFSpectra", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    Error("AddTaskTOFSpectra", "This task requires to run on ESD");
    return NULL;
  }
  
  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  
  // Create and configure the task
  AliAnalysisTaskTOFSpectra *taskTOF = new AliAnalysisTaskTOFSpectra(tname, HeavyIon, readMC, optTree, ChannelMismatch, CutVariation, SimpleCutMode);
  taskTOF->SetTreeFlag(optTree);
  taskTOF->SetHeavyIonFlag(HeavyIon);
  taskTOF->SetMCFlag(readMC);
  taskTOF->SetChannelFlag(ChannelMismatch);

  mgr->AddTask(taskTOF);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
//   TString outputFileName = AliAnalysisManager::GetCommonFileName();
//   outputFileName += ":PWGLFSpectraTOF";
//   Info("AddTaskTOFSpectra","The results of this task will be found in: %s",outputFileName.Data());  
  
  //Create and attach input
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); 
  mgr->ConnectInput(taskTOF, 0,cinput);

  //Create and attach output
  TString ListFileName = "TListTOF";
  if(readMC) ListFileName += "_MC";
  if(ChannelMismatch) ListFileName += "_Mismatch";
  AliAnalysisDataContainer *cOutput1 = mgr->CreateContainer(Form("cOutputList%s",prefix.Data()),TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s.root",ListFileName.Data()));
  mgr->ConnectOutput(taskTOF, 1, cOutput1);
  
  AliAnalysisDataContainer *cOutput2;
  if(optTree){
    TString TreeFileName = "TreeTOF";
    if(readMC) TreeFileName += "_MC";
    cOutput2 = mgr->CreateContainer("cOutputTree",TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s.root",TreeFileName.Data()));
    cOutput2->SetSpecialOutput();
    mgr->ConnectOutput(taskTOF, 2, cOutput2); 

  }
  
  //   if(optTree){
  //     AliAnalysisDataContainer *cNsigma = mgr->CreateContainer("Special", TTree::Class(),AliAnalysisManager::kOutputContainer,"Special.root");
  //     cNsigma->SetSpecialOutput();
  //     mgr->ConnectOutput(taskTOF,3,cNsigma);
  //     if(readMC){ 
  //       AliAnalysisDataContainer *cMC = mgr->CreateContainer("SpecialMC", TTree::Class(),AliAnalysisManager::kOutputContainer,"SpecialMC.root");
  //       cMC->SetSpecialOutput();
  //       mgr->ConnectOutput(taskTOF,4,cMC);
  //     }
  //   }
  
  
  return taskTOF;
}   
