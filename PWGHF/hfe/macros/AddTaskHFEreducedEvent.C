AliAnalysisTask *AddTaskHFEreducedEvent(UInt_t trigger=131073,Int_t minnbTPC=30,Int_t minnbTPCPID=80,Int_t minnbITS=2){

  //
  // Produce reduced events
  //
  

  // Name
  TString appendixx("HFEreducedEvent");
 
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  // task
  AliHFEreducedEventCreatorAOD *reducedEventCreator = new AliHFEreducedEventCreatorAOD("HFEreducedEventCreator");
  reducedEventCreator->SetMinNclustersTPC(minnbTPC);
  reducedEventCreator->SetMinNclustersTPCPID(minnbTPCPID);
  reducedEventCreator->SetMinNclustersITS(minnbITS);
  reducedEventCreator->SelectCollisionCandidates(trigger); 

  //AliHFEpidTPC *tpcpid = reducedEventCreator->GetTPCResponse();

  /*
  TF1 *etaCorrection = GetEtaCorrection();
  if(etaCorrection){
    tpcpid->SetEtaCorrection(etaCorrection);
  }
  */
  /*
  TF1 *centralityCorrection = new TF1("centralityCorrection", "pol1", 0., 10000.);
  centralityCorrection->SetParameter(0, 1.0);
  centralityCorrection->SetParameter(1, -0.00002);
  tpcpid->SetCentralityCorrection(centralityCorrection);
  */

  mgr->AddTask(reducedEventCreator);
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendixx.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(reducedEventCreator,1, mgr->CreateContainer(Form("list_%s",appendixx.Data()), TTree::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(reducedEventCreator,0, cinput );    

  return NULL;

  
}
