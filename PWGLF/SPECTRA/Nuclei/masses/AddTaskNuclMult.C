AliAnalysisTaskSE *AddTaskNuclMult(Bool_t kAOD=kTRUE){

  //for ESDs
  if(!kAOD) {
    //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
    
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
  
  const Int_t Nmult_0=12;//number of multiplicity classes (VZERO Amplitude Estimator)
  const Int_t Nmult_1=15;//number of multiplicity classes (Mid-pseudorapidity Estimator)
  
  const Int_t Ntask_0=Nmult_0+1;//+1 = +analysis on over all Minimum Bias collisions
  const Int_t Ntask_1=Nmult_1+1;//+1 = +analysis on over all Minimum Bias collisions
  
  Char_t mytaskName[100];
  snprintf(mytaskName,100,"AliAnalysisNuclMult");
  
  AliAnalysisNuclMult *task_0[Ntask_0];
  AliAnalysisNuclMult *task_1[Ntask_1];
  for(Int_t i=0;i<Ntask_0;i++) {
    task_0[i] = new AliAnalysisNuclMult(mytaskName);
    mgr->AddTask(task_0[i]);
  }
  for(Int_t i=0;i<Ntask_1;i++) {
    task_1[i] = new AliAnalysisNuclMult(mytaskName);
    mgr->AddTask(task_1[i]);
  }
  
  Float_t multiplicityRanges_0[Nmult_0+1]={0.00,0.01,0.10,1,5,10,15,20,30,40,50,70,100};
  Float_t multiplicityRanges_1[Nmult_1+1]={1,4,7,10,15,20,25,30,40,50,60,70,80,90,100,999};

  Bool_t iMultEstimator=0;
  
  for(Int_t i=0;i<Ntask_0;i++) {
    if(!kAOD) task_0[i]->SetTrackFilter(trackFilter);
    task_0[i]->SetMultEstimator(iMultEstimator);
    
    if(i<Nmult_0) task_0[i]->SetMultiplicityRange(multiplicityRanges_0[i],multiplicityRanges_0[i+1]);
    else task_0[i]->SetMultiplicityRange(-999,999);//over all Minimum Bias collisions
  }
  
  iMultEstimator=1;//pay attention to that;
  for(Int_t i=0;i<Ntask_1;i++) {
    if(!kAOD) task_1[i]->SetTrackFilter(trackFilter);
    task_1[i]->SetMultEstimator(iMultEstimator);
    
    if(i<Nmult_1) task_1[i]->SetMultiplicityRange(multiplicityRanges_1[i],multiplicityRanges_1[i+1]);
    else task_1[i]->SetMultiplicityRange(-999,999);//over all Minimum Bias collisions
  }
  
  iMultEstimator=0;//pay attention to that;
  AliAnalysisDataContainer *cinput_0[Ntask_0];
  AliAnalysisDataContainer *cOutputL_0[Ntask_0];
  for(Int_t i=0;i<Ntask_0;i++) {
    //task input
    Char_t name[1000];
    if(i<Nmult_0) snprintf(name,1000,"cchain1_kAOD=%i_isMultEs=%i_multMin=%.02f_multMax=%.02f",kAOD,iMultEstimator,multiplicityRanges_0[i],multiplicityRanges_0[i+1]);
    else snprintf(name,1000,"cchain1_kAOD=%i_MinimumBias_isMultEs=%i",kAOD,iMultEstimator);//over all Minimum Bias collisions
    cinput_0[i] = mgr->CreateContainer(name,TChain::Class(),AliAnalysisManager::kInputContainer);
    mgr->ConnectInput(task_0[i],0,mgr->GetCommonInputContainer());
    
    // task output  
    if(i<Nmult_0) snprintf(name,1000,"Results_kAOD=%i_isMultEs=%i_multMin=%.02f_multMax=%.02f",kAOD,iMultEstimator,multiplicityRanges_0[i],multiplicityRanges_0[i+1]);
    else snprintf(name,1000,"Results_kAOD=%i_MinimumBias_isMultEs=%i",kAOD,iMultEstimator);
    cOutputL_0[i] = mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task_0[i],1,cOutputL_0[i]);
  }
  
  iMultEstimator=1;//pay attention to that;
  AliAnalysisDataContainer *cinput_1[Ntask_1];
  AliAnalysisDataContainer *cOutputL_1[Ntask_1];
  for(Int_t i=0;i<Ntask_1;i++) {
    //task input
    Char_t name[1000];
    if(i<Nmult_1) snprintf(name,1000,"cchain1_kAOD=%i_isMultEs=%i_multMin=%03.0f_multMax=%03.0f",kAOD,iMultEstimator,multiplicityRanges_1[i],multiplicityRanges_1[i+1]);
    else snprintf(name,1000,"cchain1_kAOD=%i_MinimumBias_isMultEs=%i",kAOD,iMultEstimator);//over all Minimum Bias collisions
    cinput_1[i] = mgr->CreateContainer(name,TChain::Class(),AliAnalysisManager::kInputContainer);
    mgr->ConnectInput(task_1[i],0,mgr->GetCommonInputContainer());
    
    // task output  
    if(i<Nmult_1) snprintf(name,1000,"Results_kAOD=%i_isMultEs=%i_multMin=%03.0f_multMax=%03.0f",kAOD,iMultEstimator,multiplicityRanges_1[i],multiplicityRanges_1[i+1]);
    else snprintf(name,1000,"Results_kAOD=%i_MinimumBias_isMultEs=%i",kAOD,iMultEstimator);
    cOutputL_1[i] = mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task_1[i],1,cOutputL_1[i]);
  }
  
  return task_0[0];
}
