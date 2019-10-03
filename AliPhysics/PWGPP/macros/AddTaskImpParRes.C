/*AliAnalysisTaskSEImpParRes *AddTaskImpParRes(TString dirName="",
					     Bool_t readMC=kFALSE,
					     Bool_t isAOD=kTRUE,
					     Int_t SPDreq=1,
					     Int_t selPdg=-1,
					     Bool_t diamond=kTRUE,
					     Bool_t skipTrack=kTRUE,
					     Int_t minmult=0,
					     Int_t maxmult=1000000)*/
AliAnalysisTaskSEImpParRes *AddTaskImpParRes(Bool_t readMC=kFALSE,
					     Int_t selPdg=-1,
					     Bool_t diamond=kTRUE,
					     Bool_t skipTrack=kTRUE,
					     Int_t minmult=0,
					     Int_t maxmult=1000000,
					     Int_t checkSDDIsIn=1,
					     TString dirName="",
					     Bool_t isAOD=kFALSE,
					     Int_t SPDreq=1)
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

  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("d0ResAnalysisESDTrackCuts");
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  switch(SPDreq){
  case(1):
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    break;
  case(2):
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    break;
  case(3):
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    break;
  default:
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    break;
    }
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  //esdTrackCuts->SetMinNCrossedRowsTPC(70);
  esdTrackCuts->SetMinNClustersTPC(70);

  
  // Aanalysis task    
  AliAnalysisTaskSEImpParRes *d0ResTask = new AliAnalysisTaskSEImpParRes("d0ResAnalysis");
  d0ResTask->SetDebugLevel(2);
  d0ResTask->SetReadMC(readMC);
  d0ResTask->SetIsAOD(isAOD);
  d0ResTask->SetSelectedPdg(selPdg);
  d0ResTask->SetUseDiamond(diamond);
  d0ResTask->SetSkipTrack(skipTrack);
  d0ResTask->SetMultiplicityRange(minmult,maxmult);
  d0ResTask->SetESDtrackCuts(esdTrackCuts);
  mgr->AddTask(d0ResTask);

  TString fname=Form("%s:ImpParRes_Performance",mgr->GetCommonFileName());
  if(selPdg>0) {fname+=selPdg;}
  fname += dirName.Data();

  TString name=dirName;
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputd0distr = mgr->CreateContainer(Form("cinputd0distr%s",name.Data()),TChain::Class(), 
								 AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputd0ITSpureSARec = mgr->CreateContainer(Form("coutputd0ITSpureSARec_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());

  AliAnalysisDataContainer *coutputd0ITSpureSASkip = mgr->CreateContainer(Form("coutputd0ITSpureSASkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());

  AliAnalysisDataContainer *coutputd0allPointRec = mgr->CreateContainer(Form("coutputd0allPointRec_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());
  
  AliAnalysisDataContainer *coutputd0allPointSkip = mgr->CreateContainer(Form("coutputd0allPointSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(), 
								   AliAnalysisManager::kOutputContainer,
								   fname.Data());
 
  AliAnalysisDataContainer *coutputd0partPointRec = mgr->CreateContainer(Form("coutputd0partPointRec_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0partPointSkip = mgr->CreateContainer(Form("coutputd0partPointSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onepointSPDRec = mgr->CreateContainer(Form("coutputd0onepointSPDRec_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onepointSPDSkip = mgr->CreateContainer(Form("coutputd0onepointSPDSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0postvTracRec = mgr->CreateContainer(Form("coutputd0postvTracRec_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
 
  AliAnalysisDataContainer *coutputd0postvTracSkip = mgr->CreateContainer(Form("coutputd0postvTracSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0negtvTracRec = mgr->CreateContainer(Form("coutputd0negtvTracRec_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0negtvTracSkip = mgr->CreateContainer(Form("coutputd0negtvTracSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0pullAllpointRec = mgr->CreateContainer(Form("coutputd0pullAllpointRec_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


  AliAnalysisDataContainer *coutputd0pullAllpointSkip = mgr->CreateContainer(Form("coutputd0pullAllpointSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onlyRefitRec = mgr->CreateContainer(Form("coutputd0onlyRefitRec_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

  AliAnalysisDataContainer *coutputd0onlyRefitSkip = mgr->CreateContainer(Form("coutputd0onlyRefitSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


  AliAnalysisDataContainer *coutputd0sinThetaRec = mgr->CreateContainer(Form("coutputd0sinThetaRec_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


  AliAnalysisDataContainer *coutputd0sinThetaSkip = mgr->CreateContainer(Form("coutputd0sinThetaSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


 AliAnalysisDataContainer *coutputd0allPointTrue = mgr->CreateContainer(Form("coutputd0allPointTrue_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0postvTracTrue = mgr->CreateContainer(Form("coutputd0postvTracTrue_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0negtvTracTrue = mgr->CreateContainer(Form("coutputd0negtvTracTrue_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0pullAllpointTrue = mgr->CreateContainer(Form("coutputd0pullAllpointTrue_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());


 AliAnalysisDataContainer *coutputd0phiAllpointSkip = mgr->CreateContainer(Form("coutputd0phiAllpointSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());

 AliAnalysisDataContainer *coutputd0phiPostvtracSkip = mgr->CreateContainer(Form("coutputd0phiPostvtracSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
 AliAnalysisDataContainer *coutputd0phiNegtvtracSkip = mgr->CreateContainer(Form("coutputd0phiNegtvtracSkip_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									    AliAnalysisManager::kOutputContainer,
									  fname.Data());

 
  AliAnalysisDataContainer *coutputd0PID = mgr->CreateContainer(Form("coutputd0PID_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputd0Pt = mgr->CreateContainer(Form("coutputd0Pt_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
									AliAnalysisManager::kOutputContainer,
									fname.Data());
 
  AliAnalysisDataContainer *coutputNentries = mgr->CreateContainer(Form("coutputNentries_%d_%d%s",minmult,maxmult,name.Data()),TH1F::Class(),
								     AliAnalysisManager::kOutputContainer, 
								   fname.Data());

  AliAnalysisDataContainer *coutputEstimVtx = mgr->CreateContainer(Form("coutputEstimVtx_%d_%d%s",minmult,maxmult,name.Data()),TH1F::Class(),
								     AliAnalysisManager::kOutputContainer, 
								   fname.Data());

  AliAnalysisDataContainer *coutputd0withESDTC = mgr->CreateContainer(Form("coutputd0withESDTC_%d_%d%s",minmult,maxmult,name.Data()),TList::Class(),
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
  mgr->ConnectOutput(d0ResTask,26,coutputd0PID);
  mgr->ConnectOutput(d0ResTask,27,coutputd0Pt);
  mgr->ConnectOutput(d0ResTask,28,coutputNentries);
  mgr->ConnectOutput(d0ResTask,29,coutputEstimVtx);
  mgr->ConnectOutput(d0ResTask,30,coutputd0withESDTC);

  return d0ResTask;
}
