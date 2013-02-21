AliAnalysisTask *AddTaskHFEreducedEvent(UInt_t trigger=131073,Int_t minnbTPC=70,Int_t minnbTPCPID=80,Int_t minnbITS=3,Float_t nbOfSigmaTOF=3.){

  //
  // Produce reduced events
  //
  

  // Name
  TString appendixx(TString::Format("HFEreducedEventt%dTPCcl%dpidcl%dITScl%dTOFsigma%d",(Int_t)trigger,(Int_t) minnbTPC,(Int_t) minnbTPCPID,(Int_t) minnbITS,(Int_t) nbOfSigmaTOF));
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  // task
  AliHFEreducedEventCreatorAOD *reducedEventCreator = new AliHFEreducedEventCreatorAOD("HFEreducedEventCreator");
  reducedEventCreator->SetMinNclustersTPC(minnbTPC);
  reducedEventCreator->SetMinNclustersTPCPID(minnbTPCPID);
  reducedEventCreator->SetMinNclustersITS(minnbITS);
  reducedEventCreator->SetNbOfTOFSigma(nbOfSigmaTOF);
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
