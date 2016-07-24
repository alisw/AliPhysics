AliAnalysisTaskSE *AddTaskMCNuclMult(Bool_t isMC=kTRUE){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  //for ESDs
  //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

  //To set the track cuts (useful only for ESDs)
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("esdTrackCuts", "Standard2010");
  //esdTrackCuts=esdTrackCuts->GetStandardITSTPCTrackCuts2010(kFALSE,0);
  //corresponding to:
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  esdTrackCuts->SetMaxDCAToVertexZ(2);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetMaxChi2PerClusterITS(36);
  //DCA cuts set in the task
    
  AliPPVsMultUtils *fAliPPVsMultUtils = new AliPPVsMultUtils();

  const Int_t Nmult_0 = 12;//0;//number of multiplicity bins (VZERO Amplitude Estimator)
  const Int_t Ntask_0 = Nmult_0+1;//Nmult_0;//+1 for integrated multiplicity
  Float_t multiplicityRanges_0[Nmult_0+1] = {0,5,10,15,20,25,30,35,40,45,50,70,100};//{0,100};

  AliAnalysisMCNuclMult *task_0[Ntask_0];
  for(Int_t i=0;i<Ntask_0;i++) {
    task_0[i] = new AliAnalysisMCNuclMult("AliAnalysisMCNuclMult");
    mgr->AddTask(task_0[i]);
  }
  
  for(Int_t i=0;i<Ntask_0;i++) {
    task_0[i]->SetESDtrackCutsObj(esdTrackCuts);
    task_0[i]->SetPPVsMultUtilsObj(fAliPPVsMultUtils);
    //task_0[i]->SetMultiplicityRange(multiplicityRanges_0[i],multiplicityRanges_0[i+1]);
    if(i<Nmult_0) task_0[i]->SetMultiplicityRange(multiplicityRanges_0[i],multiplicityRanges_0[i+1]);
    else task_0[i]->SetMultiplicityRange(0,100);
  }
  
  AliAnalysisDataContainer *cinput_0[Ntask_0];
  AliAnalysisDataContainer *cOutputL_0[Ntask_0];
  for(Int_t i=0;i<Ntask_0;i++) {
    //input
    Char_t name[1000];
    //snprintf(name,1000,"cchain1_multMin=%.02f_multMax=%.02f",multiplicityRanges_0[i],multiplicityRanges_0[i+1]);
    if(i<Nmult_0) snprintf(name,1000,"cchain1_multMin=%.02f_multMax=%.02f",multiplicityRanges_0[i],multiplicityRanges_0[i+1]);
    else snprintf(name,1000,"cchain1_multMin=%.02f_multMax=%.02f",0.,100.);
    cinput_0[i] = mgr->CreateContainer(name,TChain::Class(),AliAnalysisManager::kInputContainer);
    mgr->ConnectInput(task_0[i],0,mgr->GetCommonInputContainer());
    
    //output  
    //snprintf(name,1000,"multMin=%.02f_multMax=%.02f",multiplicityRanges_0[i],multiplicityRanges_0[i+1]);
    if(i<Nmult_0) snprintf(name,1000,"multMin=%.02f_multMax=%.02f",multiplicityRanges_0[i],multiplicityRanges_0[i+1]);
    else snprintf(name,1000,"multMin=%.02f_multMax=%.02f",0.,100.);
    cOutputL_0[i] = mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task_0[i],1,cOutputL_0[i]);
  }
  
  return task_0[0];
}
