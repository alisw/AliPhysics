AliAnalysisTaskSE *AddTaskNuclMult(Bool_t kAOD=kTRUE){

  //for ESDs
  if(!kAOD) {
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
    
    //To set the track cuts (useful only for ESDs)
    AliESDtrackCuts* esdTrackCutsStd2010 = new AliESDtrackCuts("AliESDtrackCuts", "Standard2010");
    esdTrackCutsStd2010->GetStandardITSTPCTrackCuts2010(kFALSE,0);
    esdTrackCutsStd2010->SetMaxDCAToVertexXY(2.4);
    esdTrackCutsStd2010->SetMaxDCAToVertexZ(3.2);
    esdTrackCutsStd2010->SetDCAToVertex2D(kTRUE);
    
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsStd2010);
    //task->SetTrackFilter(trackFilter);
  }
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  const Int_t Nmult=15;//number of multiplicity classes
  const Int_t Ntask=Nmult+1;//+1 = +analysis on over all Minimum Bias collisions
  
  Char_t mytaskName[100];
  snprintf(mytaskName,100,"AliAnalysisNuclMult");
  
  AliAnalysisNuclMult *task[Ntask];
  for(Int_t i=0;i<Ntask;i++) {
    task[i] = new AliAnalysisNuclMult(mytaskName);
    mgr->AddTask(task[i]);
  }
  
  Int_t multiplicityRanges[Nmult+1]={1,4,7,10,15,20,25,30,40,50,60,70,80,90,100,999};
  
  for(Int_t i=0;i<Ntask;i++) {
    if(!kAOD) task[i]->SetTrackFilter(trackFilter);
    
    if(i<Nmult) task[i]->SetMultiplicityRange(multiplicityRanges[i],multiplicityRanges[i+1]);
    else task[i]->SetMultiplicityRange(-999,999);//over all Minimum Bias collisions
  }
  
  AliAnalysisDataContainer *cinput[Ntask];
  AliAnalysisDataContainer *cOutputL[Ntask];
  
  for(Int_t i=0;i<Ntask;i++) {
    //task input
    Char_t name[1000];
    if(i<Nmult) snprintf(name,1000,"cchain1_kAOD=%i_multMin=%03i_multMax=%03i",kAOD,multiplicityRanges[i],multiplicityRanges[i+1]);
    else snprintf(name,1000,"cchain1_kAOD=%i_MinimumBias",kAOD);//over all Minimum Bias collisions
    cinput[i] = mgr->CreateContainer(name,TChain::Class(),AliAnalysisManager::kInputContainer);
    mgr->ConnectInput(task[i],0,mgr->GetCommonInputContainer());
    
    // task output  
    if(i<Nmult) snprintf(name,1000,"Results_kAOD=%i_multMin=%03i_multMax=%03i",kAOD,multiplicityRanges[i],multiplicityRanges[i+1]);
    else snprintf(name,1000,"Results_kAOD=%i_MinimumBias",kAOD);
    cOutputL[i]= mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task[i],1,cOutputL[i]);
  }
  
  return task[0];
}
