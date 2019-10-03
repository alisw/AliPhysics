AliAnalysisTaskSEITSsaSpectraMultiplicity *AddTaskITSsaSpectraMultiplicity(Int_t optTree=0,Bool_t readMC=0,Bool_t useV0=0,Double_t lowmult=-1, Double_t upmult=-1,Int_t hi=0, TString prefix = ""){
  // Creates, configures and attaches to the train the task for pi, K , p spectra
  // with ITS standalone tracks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  
  TString multest[2] = {"RefMult","V0Mult"};
  TString multinter = "";
  if(useV0==kFALSE) multinter = Form("%.0fto%.0f",lowmult,upmult);    
  else{
    if(lowmult<0) multinter += Form("%.0f",lowmult);
    else if(lowmult <0.1) multinter += Form("%.2f",lowmult);
    else if(lowmult <1) multinter += Form("%.1f",lowmult);
    else multinter += Form("%.0f",lowmult);
    multinter+="to";
    if(upmult<0) multinter += Form("%.0f",upmult);
    else if(upmult <0.1) multinter += Form("%.2f",upmult);
    else if(upmult <1) multinter += Form("%.1f",upmult);
    else multinter += Form("%.0f",upmult);
  }
  Int_t estindex = 0;
  if(useV0) estindex = 1;
  
  ::Info("AddTaskITSsaSpectraMultiplicity","Adding a new task with this settings optTree = %i, readMC = %i, useV0 = %i Multiplicity Interval = %s hi = %i",optTree,readMC,useV0,multinter.Data(),hi);  
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskITSsaSpectraMultiplicity", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskITSsaSpectraMultiplicity", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskITSsaSpectraMultiplicity", "This task requires to run on ESD");
    return NULL;
  }
  
  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  // Create and configure the task
  
  AliAnalysisTaskSEITSsaSpectraMultiplicity *taskits = new AliAnalysisTaskSEITSsaSpectraMultiplicity();
  taskits->SelectCollisionCandidates();
  taskits->SetReadMC(readMC);
  taskits->SetFillTree(optTree);
  taskits->SetV0Estimator(useV0);
  
  
  if(hi)
  {
    Float_t lowcencut=(Float_t)lowmult;
    Float_t upcencut=(Float_t)upmult;
    taskits->SetCentralityCut(lowcencut,upcencut);
    taskits->SetHImode();
  }
  else{
    taskits->SetMultBin(lowmult,upmult);
  }
  
  taskits->PrintStatus();
  
  mgr->AddTask(taskits);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWG2SpectraITSsa";
  ::Info("AddTaskITSsaSpectraMultiplicity","The results of this task will be found in: %s",outputFileName.Data());  
  
  AliAnalysisDataContainer *coutput = 0x0;
  AliAnalysisDataContainer *coutputCuts = 0x0;
  
  if(hi){
    coutput = mgr->CreateContainer(Form("clistITSsaCent%.0fto%.0f",lowmult,upmult),
				   TList::Class(),
				   AliAnalysisManager::kOutputContainer,
				   outputFileName );
  }
  else{
        
    coutput = mgr->CreateContainer(Form("clistITSsa%s%s%s",multest[estindex].Data(),multinter.Data(),prefix.Data()),
				   TList::Class(),
				   AliAnalysisManager::kOutputContainer,
				   outputFileName );    
    
    coutputCuts = mgr->CreateContainer(Form("DCACut%s%s%s",multest[estindex].Data(),multinter.Data(),prefix.Data()),
				       TList::Class(),
				       AliAnalysisManager::kParamContainer,
				       outputFileName );
  }
  

  mgr->ConnectInput(taskits, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskits, 1, coutput);
  mgr->ConnectOutput(taskits, 2, coutputCuts);
  
  if(optTree){
    AliAnalysisDataContainer *cNsigma = mgr->CreateContainer("SpecialNSigmaOutput", TTree::Class(),AliAnalysisManager::kOutputContainer,"SpecialNSigmaOutput.root");
    cNsigma->SetSpecialOutput();
    mgr->ConnectOutput(taskits,3,cNsigma);
    if(readMC){ 
      AliAnalysisDataContainer *cMC = mgr->CreateContainer("SpecialMCOutput", TTree::Class(),AliAnalysisManager::kOutputContainer,"SpecialMCOutput.root");
      cMC->SetSpecialOutput();
      mgr->ConnectOutput(taskits,4,cMC);
    }
  }
  
  
  return taskits;
}   
