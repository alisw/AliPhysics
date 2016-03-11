AliAnalysisTaskJetsEvshape* AddTaskJetsEvshape(const char *ntracks            = "Tracks",
					       const char *nclusters          = "CaloClusters",
					       const char *njets              = "Jets",
					       const char *name               = "jets_evshape")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if(!mgr){
    ::Error("AddTaskJetsEvshape", "No analysis manager to connect to.");
    return 0x0;
  }
  if(!mgr->GetInputEventHandler()){
    ::Error("AddTaskJetsEvshape", "This task requires an input event handler.");
    return 0x0;
  }

  AliAnalysisTaskJetsEvshape *task = new AliAnalysisTaskJetsEvshape(name);

  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer(Form("emcal_%s", name), TList::Class(), AliAnalysisManager::kOutputContainer,
                         Form("%s:PWGJE_jets_evshape", AliAnalysisManager::GetCommonFileName()));

  AliAnalysisDataContainer *coutput2 =
    mgr->CreateContainer(Form("spec_%s", name), TList::Class(), AliAnalysisManager::kOutputContainer,
                         Form("%s:PWGJE_jets_evshape", AliAnalysisManager::GetCommonFileName()));

  if (!coutput || !coutput2) {
    ::Error("AddTaskJetsEvshape", "output containers could not be created");
    return 0x0;
  }

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  if(trackCont) trackCont->SetClassName("AliVTrack");
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);

  AliJetContainer *jetCont = task->AddJetContainer(njets, "TPC", 0.4);
  printf("just added jet container %p\n", jetCont);
  if(jetCont) {
    jetCont->SetRhoName("");
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
    jetCont->SetZLeadingCut(0.98,0.98);
    // jetCont->SetPercAreaCut(0.6);
    jetCont->SetJetPtCut(1.);
    jetCont->SetLeadingHadronType(0);
  }

  mgr->AddTask(task);

  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
  if (mgr->GetCommonOutputContainer())
    mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(task, 1, coutput);
  mgr->ConnectOutput(task, 2, coutput2);

  return task;
}
