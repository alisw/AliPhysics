AliAnalysisTaskSEImpParRes *AddTaskImpParRes(Bool_t readMC=kFALSE,Int_t selPdg=-1,Bool_t diamond=kTRUE)
{
  //
  // Configuration for the study of the impact parameter resolution
  //
  // xianbao.yuan@pd.infn.it, andrea.dainese@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskImpParRes", "No analysis manager to connect to.");
    return NULL;
  }   

  
  // Aanalysis task    
  AliAnalysisTaskSEImpParRes *d0ResTask = new AliAnalysisTaskSEImpParRes("d0ResAnalysis");
  d0ResTask->SetDebugLevel(2);
  d0ResTask->SetReadMC(readMC);
  d0ResTask->SetSelectedPdg(selPdg);
  d0ResTask->SetUseDiamond(diamond);
  mgr->AddTask(d0ResTask);
 
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputd0distr = mgr->CreateContainer("cinputd0distr",TChain::Class(), 
								 AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputd0ITSpureSARec = mgr->CreateContainer("coutputd0ITSpureSARec",TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   "ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputd0ITSpureSASkip = mgr->CreateContainer("coutputd0ITSpureSASkip",TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   "ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputd0allPointRec = mgr->CreateContainer("coutputd0allPointRec",TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   "ImpParRes.Performance.root");
  
  AliAnalysisDataContainer *coutputd0allPointSkip = mgr->CreateContainer("coutputd0allPointSkip",TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   "ImpParRes.Performance.root");
 
  AliAnalysisDataContainer *coutputd0partPointRec = mgr->CreateContainer("coutputd0partPointRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputd0partPointSkip = mgr->CreateContainer("coutputd0partPointSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputd0onepointSPDRec = mgr->CreateContainer("coutputd0onepointSPDRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputd0onepointSPDSkip = mgr->CreateContainer("coutputd0onepointSPDSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputd0postvTracRec = mgr->CreateContainer("coutputd0postvTracRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");
 
 
  AliAnalysisDataContainer *coutputd0postvTracSkip = mgr->CreateContainer("coutputd0postvTracSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");
 
  AliAnalysisDataContainer *coutputd0negtvTracRec = mgr->CreateContainer("coutputd0negtvTracRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");
 
  AliAnalysisDataContainer *coutputd0negtvTracSkip = mgr->CreateContainer("coutputd0negtvTracSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");


 

  AliAnalysisDataContainer *coutputd0pullAllpointRec = mgr->CreateContainer("coutputd0pullAllpointRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");


  AliAnalysisDataContainer *coutputd0pullAllpointSkip = mgr->CreateContainer("coutputd0pullAllpointSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputd0onlyRefitRec = mgr->CreateContainer("coutputd0onlyRefitRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputd0onlyRefitSkip = mgr->CreateContainer("coutputd0onlyRefitSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputd0sinThetaRec = mgr->CreateContainer("coutputd0sinThetaRec",TList::Class(),
 									AliAnalysisManager::kOutputContainer,
 									"ImpParRes.Performance.root");
  
  AliAnalysisDataContainer *coutputd0sinThetaSkip = mgr->CreateContainer("coutputd0sinThetaSkip",TList::Class(),
									 AliAnalysisManager::kOutputContainer,
									 "ImpParRes.Performance.root");
  
  AliAnalysisDataContainer *coutputd0Pt = mgr->CreateContainer("coutputd0Pt",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									"ImpParRes.Performance.root");

  
  AliAnalysisDataContainer *coutputNentries = mgr->CreateContainer("coutputNentries",TH1F::Class(),
								     AliAnalysisManager::kOutputContainer, 
								   "ImpParRes.Performance.root");

  AliAnalysisDataContainer *coutputEstimVtx = mgr->CreateContainer("coutputEstimVtx",TH1F::Class(),
								     AliAnalysisManager::kOutputContainer, 
								   "ImpParRes.Performance.root");
 

 
  mgr->ConnectInput(d0ResTask,0,mgr->GetCommonInputContainer()); 
  //mgr->ConnectOutput(d0ResTask,0,mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(d0ResTask,1,coutputd0ITSpureSARec);
  mgr->ConnectOutput(d0ResTask,2,coutputd0ITSpureSASkip);
  mgr->ConnectOutput(d0ResTask,3,coutputd0allPointRec);
  mgr->ConnectOutput(d0ResTask,4,coutputd0allPointSkip);
  mgr->ConnectOutput(d0ResTask,5,coutputd0partPointRec);
  mgr->ConnectOutput(d0ResTask,6,coutputd0partPointSkip);
  mgr->ConnectOutput(d0ResTask,7,coutputd0onepointSPDRec);
  mgr->ConnectOutput(d0ResTask,8,coutputd0onepointSPDSkip);
  mgr->ConnectOutput(d0ResTask,9,coutputd0postvTracRec); 
  mgr->ConnectOutput(d0ResTask,10,coutputd0postvTracSkip);
  mgr->ConnectOutput(d0ResTask,11,coutputd0negtvTracRec);
  mgr->ConnectOutput(d0ResTask,12,coutputd0negtvTracSkip);
  mgr->ConnectOutput(d0ResTask,13,coutputd0pullAllpointRec);
  mgr->ConnectOutput(d0ResTask,14,coutputd0pullAllpointSkip);
  mgr->ConnectOutput(d0ResTask,15,coutputd0onlyRefitRec);
  mgr->ConnectOutput(d0ResTask,16,coutputd0onlyRefitSkip);
  mgr->ConnectOutput(d0ResTask,17,coutputd0sinThetaRec);
  mgr->ConnectOutput(d0ResTask,18,coutputd0sinThetaSkip);
  mgr->ConnectOutput(d0ResTask,19,coutputd0Pt);
  mgr->ConnectOutput(d0ResTask,20,coutputNentries);
  mgr->ConnectOutput(d0ResTask,21,coutputEstimVtx);

  return d0ResTask;
}
