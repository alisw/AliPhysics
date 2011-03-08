AliAnalysisTaskSEImpParRes *AddTaskImpParRes(Bool_t readMC=kFALSE,
					     Int_t selPdg=-1,
					     Bool_t diamond=kTRUE,
					     Bool_t skipTrack=kTRUE,
					     Int_t minmult=0,
					     Int_t maxmult=1000000)
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
  d0ResTask->SetSkipTrack(skipTrack);
  d0ResTask->SetMultiplicityRange(minmult,maxmult);
  mgr->AddTask(d0ResTask);

  TString fname=Form("%s:ImpParRes_Performance",mgr->GetCommonFileName());
  if(selPdg>0) {fname+=selPdg;}

 
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputd0distr = mgr->CreateContainer("cinputd0distr",TChain::Class(), 
								 AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputd0ITSpureSARec = mgr->CreateContainer("coutputd0ITSpureSARec",TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());

  AliAnalysisDataContainer *coutputd0ITSpureSASkip = mgr->CreateContainer("coutputd0ITSpureSASkip",TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());

  AliAnalysisDataContainer *coutputd0allPointRec = mgr->CreateContainer("coutputd0allPointRec",TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());
  
  AliAnalysisDataContainer *coutputd0allPointSkip = mgr->CreateContainer("coutputd0allPointSkip",TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());
 
  AliAnalysisDataContainer *coutputd0partPointRec = mgr->CreateContainer("coutputd0partPointRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0partPointSkip = mgr->CreateContainer("coutputd0partPointSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onepointSPDRec = mgr->CreateContainer("coutputd0onepointSPDRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onepointSPDSkip = mgr->CreateContainer("coutputd0onepointSPDSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0postvTracRec = mgr->CreateContainer("coutputd0postvTracRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
 
  AliAnalysisDataContainer *coutputd0postvTracSkip = mgr->CreateContainer("coutputd0postvTracSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0negtvTracRec = mgr->CreateContainer("coutputd0negtvTracRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0negtvTracSkip = mgr->CreateContainer("coutputd0negtvTracSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0pullAllpointRec = mgr->CreateContainer("coutputd0pullAllpointRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


  AliAnalysisDataContainer *coutputd0pullAllpointSkip = mgr->CreateContainer("coutputd0pullAllpointSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onlyRefitRec = mgr->CreateContainer("coutputd0onlyRefitRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onlyRefitSkip = mgr->CreateContainer("coutputd0onlyRefitSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


  AliAnalysisDataContainer *coutputd0sinThetaRec = mgr->CreateContainer("coutputd0sinThetaRec",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


  AliAnalysisDataContainer *coutputd0sinThetaSkip = mgr->CreateContainer("coutputd0sinThetaSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


 AliAnalysisDataContainer *coutputd0allPointTrue = mgr->CreateContainer("coutputd0allPointTrue",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0postvTracTrue = mgr->CreateContainer("coutputd0postvTracTrue",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0negtvTracTrue = mgr->CreateContainer("coutputd0negtvTracTrue",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0pullAllpointTrue = mgr->CreateContainer("coutputd0pullAllpointTrue",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


 AliAnalysisDataContainer *coutputd0phiAllpointSkip = mgr->CreateContainer("coutputd0phiAllpointSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0phiPostvtracSkip = mgr->CreateContainer("coutputd0phiPostvtracSkip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
 AliAnalysisDataContainer *coutputd0phiNegtvtracSkip = mgr->CreateContainer("coutputd0phiNegtvtracSkip",TList::Class(),
									    AliAnalysisManager::kOutputContainer,
									  fname.Data());

 
 AliAnalysisDataContainer *coutputd0clusterTypeSPD01Skip = mgr->CreateContainer("coutputd0clusterTypeSPD01Skip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD02Skip = mgr->CreateContainer("coutputd0clusterTypeSPD02Skip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD03Skip = mgr->CreateContainer("coutputd0clusterTypeSPD03Skip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD11Skip = mgr->CreateContainer("coutputd0clusterTypeSPD11Skip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD12Skip = mgr->CreateContainer("coutputd0clusterTypeSPD12Skip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD13Skip = mgr->CreateContainer("coutputd0clusterTypeSPD13Skip",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0PID = mgr->CreateContainer("coutputd0PID",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0Pt = mgr->CreateContainer("coutputd0Pt",TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputNentries = mgr->CreateContainer("coutputNentries",TH1F::Class(),
								     AliAnalysisManager::kOutputContainer, 
								   fname.Data());

  AliAnalysisDataContainer *coutputEstimVtx = mgr->CreateContainer("coutputEstimVtx",TH1F::Class(),
								     AliAnalysisManager::kOutputContainer, 
								   fname.Data());

  // Attach input  
  mgr->ConnectInput(d0ResTask,0,mgr->GetCommonInputContainer()); 
  // Attack output
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
  mgr->ConnectOutput(d0ResTask,19,coutputd0allPointTrue);
  mgr->ConnectOutput(d0ResTask,20,coutputd0postvTracTrue);
  mgr->ConnectOutput(d0ResTask,21,coutputd0negtvTracTrue);
  mgr->ConnectOutput(d0ResTask,22,coutputd0pullAllpointTrue);
  mgr->ConnectOutput(d0ResTask,23,coutputd0phiAllpointSkip);
  mgr->ConnectOutput(d0ResTask,24,coutputd0phiPostvtracSkip);
  mgr->ConnectOutput(d0ResTask,25,coutputd0phiNegtvtracSkip);
  mgr->ConnectOutput(d0ResTask,26,coutputd0clusterTypeSPD01Skip);
  mgr->ConnectOutput(d0ResTask,27,coutputd0clusterTypeSPD02Skip);
  mgr->ConnectOutput(d0ResTask,28,coutputd0clusterTypeSPD03Skip);
  mgr->ConnectOutput(d0ResTask,29,coutputd0clusterTypeSPD11Skip);
  mgr->ConnectOutput(d0ResTask,30,coutputd0clusterTypeSPD12Skip);
  mgr->ConnectOutput(d0ResTask,31,coutputd0clusterTypeSPD13Skip);
  mgr->ConnectOutput(d0ResTask,32,coutputd0PID);
  mgr->ConnectOutput(d0ResTask,33,coutputd0Pt);
  mgr->ConnectOutput(d0ResTask,34,coutputNentries);
  mgr->ConnectOutput(d0ResTask,35,coutputEstimVtx);

  return d0ResTask;
}
