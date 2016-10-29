AliAnalysisTaskSE *AddTaskNuclMult(){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  //for ESDs
  //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

  //To set the track cuts (useful only for ESDs)
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("esdTrackCuts", "Standard2010");
  //esdTrackCuts=esdTrackCuts->GetStandardITSTPCTrackCuts2010(kFALSE,0);
  //corresponding to (with SetMaxDCAToVertexZ(2) also):
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  //esdTrackCuts->SetMaxDCAToVertexZ(2);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetMaxChi2PerClusterITS(36);
  //DCA cuts set in the task
    
  AliPPVsMultUtils *fAliPPVsMultUtils = new AliPPVsMultUtils();

  const Int_t Ntask = 2;
  
  Float_t DCAzMax[Ntask] = {1.,10.};

  AliAnalysisNuclMult *task[Ntask];
  for(Int_t i=0;i<Ntask;i++) {
    task[i] = new AliAnalysisNuclMult("AliAnalysisNuclMult");
    mgr->AddTask(task[i]);
  }
  
  for(Int_t i=0;i<Ntask;i++) {
    task[i]->SetESDtrackCutsObj(esdTrackCuts);
    task[i]->SetPPVsMultUtilsObj(fAliPPVsMultUtils);
    task[i]->SetMultiplicityRange(0,100);
    task[i]->SetDCAzMax(DCAzMax[i]);
  }
  
  AliAnalysisDataContainer *cinput[Ntask];
  AliAnalysisDataContainer *cOutputL[Ntask];
  for(Int_t i=0;i<Ntask;i++) {
    //input
    Char_t name[1000];
    snprintf(name,1000,"cchain1_out_DCAzMax=%.1f",DCAzMax[i]);
    cinput[i] = mgr->CreateContainer(name,TChain::Class(),AliAnalysisManager::kInputContainer);
    mgr->ConnectInput(task[i],0,mgr->GetCommonInputContainer());
    
    //output  
    snprintf(name,1000,"out_DCAzMax=%.1f",DCAzMax[i]);
    cOutputL[i] = mgr->CreateContainer(name,TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
    mgr->ConnectOutput(task[i],1,cOutputL[i]);
  }
  
  return task[0];
}
