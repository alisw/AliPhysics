AliAnalysisTaskSEITSsaSpectra *AddTaskITSsaSpectra(Int_t optNtuple=0,Int_t readMC=0,Int_t lowcut=-1,Int_t upcut=-1,Int_t hi=0){
  // Creates, configures and attaches to the train the task for pi, K , p spectra
  // with ITS standalone tracks
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskITSsaSpectra", "No analysis manager to connect to.");
    return NULL;
  }   
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskITSsaSpectra", "This task requires an input event handler");
    return NULL;
  }   
  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskITSsaSpectra", "This task requires to run on ESD");
    return NULL;
  }
  
  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = new AliMCEventHandler;
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  // Create and configure the task
  
  AliAnalysisTaskSEITSsaSpectra *taskits = new AliAnalysisTaskSEITSsaSpectra();
  taskits->SelectCollisionCandidates();
  taskits->SetReadMC(readMC);
  taskits->SetFillNtuple(optNtuple);
  
  if(hi)
    {
      Float_t lowcencut=(Float_t)lowcut;
      Float_t upcencut=(Float_t)upcut;
      taskits->SetCentralityCut(lowcencut,upcencut);
      taskits->SetHImode();
    }
  else{
    taskits->SetMultBin(lowcut,upcut);
  }
  
  mgr->AddTask(taskits);
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWG2SpectraITSsa";
  
  AliAnalysisDataContainer *coutput =0x0;
  
  if(hi)
    {
      coutput = mgr->CreateContainer(Form("clistITSsaCent%ito%i",lowcut,upcut),
				      TList::Class(),
				      AliAnalysisManager::kOutputContainer,
				      outputFileName );
    }
  else
    {
      coutput = mgr->CreateContainer(Form("clistITSsaMult%ito%i",lowcut,upcut),
				     TList::Class(),
				      AliAnalysisManager::kOutputContainer,
				      outputFileName );    
    }
  
  mgr->ConnectInput(taskits, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskits, 1, coutput);
  return taskits;
}   
