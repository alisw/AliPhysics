// 
// // AddTaskSimSpectraLFpPb
// 
AliAnalysisTask *AddTaskSimSpectraLFpPb(TString suffixName ="",Float_t ly, Float_t hy){
  
  // Create the task, add it to manager and configure it
  //===========================================================================
  
  AliAnalysisTaskSimSpectraLFpPb* taskSpectraLFMC = new  AliAnalysisTaskSimSpectraLFpPb("AliAnalysisTaskSpectraLFMC");
  taskSpectraLFMC -> SetYRange(ly, hy);
  
  // Get the pointer to the existing analysis manager via the static access method
  //===========================================================================
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){ Printf("AliAnalysisTaskSimSpectraLF: No analysis manager to connect to."); return NULL; }

  // Check the analysis type using the event handlers connected to the analysis manager
  //===========================================================================

  if(!mgr->GetMCtruthEventHandler()){ Printf("AliAnalysisTaskSimSpectraLFpPb: This task requires an input MC event handler."); return NULL; }

  // ADD the task
  //===========================================================================
  mgr -> AddTask(taskSpectraLFMC);  

  // Create containers for input/output

  TString finDirname	= "";
  TString inname	= "cinput";
  TString outBasic	= "cList";

  finDirname	+= suffixName.Data();
  inname	+= finDirname.Data();
  outBasic	+= finDirname.Data();

  
  // Input and Output Slots
  //===========================================================================

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWGLF_SimSpectraLF";
  
  AliAnalysisDataContainer *coutSim = mgr->CreateContainer(outBasic,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  //  AliAnalysisDataContainer *coutSim1 = mgr->CreateContainer(outBasic,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  mgr->ConnectInput (taskSpectraLFMC, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskSpectraLFMC, 1, coutSim);
  // mgr->ConnectOutput(taskSpectraLFMC, 2, coutSim1);


  // create containers for input/output                                                                                                                              
  
  
  return taskSpectraLFMC;
  
}

