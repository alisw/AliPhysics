AliAnalysisTask *AddTaskHFEreducedEvent(Bool_t MCthere=kFALSE, Int_t TRDtrigger=0,Int_t minnbTPC=70,Int_t minnbTPCPID=80,Int_t minnbITS=3,
					Bool_t isRemoveFirstEvent=kTRUE){

  //
  // Produce reduced events
  //
  

  // Name
    TString appendixx(TString::Format("HFEreducedEventt%dTPCcl%dpidcl%dITScl%d",(Int_t)TRDtrigger,(Int_t) minnbTPC,(Int_t) minnbTPCPID,(Int_t) minnbITS,
				      (Int_t) isRemoveFirstEvent));
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

  // task
  AliHFEreducedEventCreatorESD *reducedEventCreator = new AliHFEreducedEventCreatorESD("HFEreducedEventCreator");
  reducedEventCreator->SetMinNclustersTPC(minnbTPC);
  reducedEventCreator->SetMinNclustersTPCPID(minnbTPCPID);
  reducedEventCreator->SetMinNclustersITS(minnbITS);
  if(isRemoveFirstEvent) reducedEventCreator->SetRemoveFirstEventFromChunk();

  if(TRDtrigger==0) reducedEventCreator->SelectCollisionCandidates(AliVEvent::kINT7);
  else reducedEventCreator->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kTRD);

  mgr->AddTask(reducedEventCreator);
  
  TString containerName = mgr->GetCommonFileName();
  containerName += ":";
  containerName += appendixx.Data();

  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  mgr->ConnectOutput(reducedEventCreator,1, mgr->CreateContainer(Form("list_%s",appendixx.Data()), TTree::Class(),AliAnalysisManager::kOutputContainer,containerName.Data()));
  mgr->ConnectInput(reducedEventCreator,0, cinput );    

  return NULL;

  
}
