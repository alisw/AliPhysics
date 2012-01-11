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

  AliAnalysisDataContainer *coutputd0ITSpureSARec = mgr->CreateContainer(Form("coutputd0ITSpureSARec_%d_%d",minmult,maxmult),TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());

  AliAnalysisDataContainer *coutputd0ITSpureSASkip = mgr->CreateContainer(Form("coutputd0ITSpureSASkip_%d_%d",minmult,maxmult),TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());

  AliAnalysisDataContainer *coutputd0allPointRec = mgr->CreateContainer(Form("coutputd0allPointRec_%d_%d",minmult,maxmult),TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());
  
  AliAnalysisDataContainer *coutputd0allPointSkip = mgr->CreateContainer(Form("coutputd0allPointSkip_%d_%d",minmult,maxmult),TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());
 
  AliAnalysisDataContainer *coutputd0partPointRec = mgr->CreateContainer(Form("coutputd0partPointRec_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0partPointSkip = mgr->CreateContainer(Form("coutputd0partPointSkip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onepointSPDRec = mgr->CreateContainer(Form("coutputd0onepointSPDRec_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onepointSPDSkip = mgr->CreateContainer(Form("coutputd0onepointSPDSkip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0postvTracRec = mgr->CreateContainer(Form("coutputd0postvTracRec_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
 
  AliAnalysisDataContainer *coutputd0postvTracSkip = mgr->CreateContainer(Form("coutputd0postvTracSkip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0negtvTracRec = mgr->CreateContainer(Form("coutputd0negtvTracRec_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0negtvTracSkip = mgr->CreateContainer(Form("coutputd0negtvTracSkip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0pullAllpointRec = mgr->CreateContainer(Form("coutputd0pullAllpointRec_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


  AliAnalysisDataContainer *coutputd0pullAllpointSkip = mgr->CreateContainer(Form("coutputd0pullAllpointSkip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onlyRefitRec = mgr->CreateContainer(Form("coutputd0onlyRefitRec_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onlyRefitSkip = mgr->CreateContainer(Form("coutputd0onlyRefitSkip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


  AliAnalysisDataContainer *coutputd0sinThetaRec = mgr->CreateContainer(Form("coutputd0sinThetaRec_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


  AliAnalysisDataContainer *coutputd0sinThetaSkip = mgr->CreateContainer(Form("coutputd0sinThetaSkip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


 AliAnalysisDataContainer *coutputd0allPointTrue = mgr->CreateContainer(Form("coutputd0allPointTrue_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0postvTracTrue = mgr->CreateContainer(Form("coutputd0postvTracTrue_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0negtvTracTrue = mgr->CreateContainer(Form("coutputd0negtvTracTrue_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0pullAllpointTrue = mgr->CreateContainer(Form("coutputd0pullAllpointTrue_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


 AliAnalysisDataContainer *coutputd0phiAllpointSkip = mgr->CreateContainer(Form("coutputd0phiAllpointSkip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0phiPostvtracSkip = mgr->CreateContainer(Form("coutputd0phiPostvtracSkip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
 AliAnalysisDataContainer *coutputd0phiNegtvtracSkip = mgr->CreateContainer(Form("coutputd0phiNegtvtracSkip_%d_%d",minmult,maxmult),TList::Class(),
									    AliAnalysisManager::kOutputContainer,
									  fname.Data());

 
 AliAnalysisDataContainer *coutputd0clusterTypeSPD01Skip = mgr->CreateContainer(Form("coutputd0clusterTypeSPD01Skip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD02Skip = mgr->CreateContainer(Form("coutputd0clusterTypeSPD02Skip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD03Skip = mgr->CreateContainer(Form("coutputd0clusterTypeSPD03Skip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD11Skip = mgr->CreateContainer(Form("coutputd0clusterTypeSPD11Skip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD12Skip = mgr->CreateContainer(Form("coutputd0clusterTypeSPD12Skip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0clusterTypeSPD13Skip = mgr->CreateContainer(Form("coutputd0clusterTypeSPD13Skip_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0PID = mgr->CreateContainer(Form("coutputd0PID_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0Pt = mgr->CreateContainer(Form("coutputd0Pt_%d_%d",minmult,maxmult),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputNentries = mgr->CreateContainer(Form("coutputNentries_%d_%d",minmult,maxmult),TH1F::Class(),
								     AliAnalysisManager::kOutputContainer, 
								   fname.Data());

  AliAnalysisDataContainer *coutputEstimVtx = mgr->CreateContainer(Form("coutputEstimVtx_%d_%d",minmult,maxmult),TH1F::Class(),
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
