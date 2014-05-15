AliAnalysisTaskPPJetSpectra *jetAnal = 0;

AliAnalysisTaskPPJetSpectra *AddTaskPPJetSpectra(TString name = "", TString recBranch = "", TString genBranch = "", TString recBckg = "", TString genBckg= "", UInt_t fEvtSelMask=AliVEvent::kMB, Float_t fEvtClass = -1,  Int_t nVtxContCut=0, Double_t fVtxZcut = 10, Double_t fVtxRcut = 1, UInt_t nTrackFilter = 272, Double_t trackPtMin = 0.15, Double_t trackPtMax = 100., Double_t trackEtaAbs = 0.9, Double_t fJetPtMin= 5, Double_t fJetEta = 0.9, Double_t fJetZ = 0.99)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) 
  {
    ::Error("AddTaskJetAnalysis", "No Analysis Manager found");
    return NULL;
  }
  AliAnalysisTaskPPJetSpectra *task = 0;

  if(name.Length() == 0) name = Form("JetAnalysis-%s-%s-%s-%s",recBranch.Data(),genBranch.Data(), recBckg.Data(), genBckg.Data() );
  task = new AliAnalysisTaskPPJetSpectra(name.Data());
  task->SetRecJetBranch(recBranch);
  task->SetGenJetBranch(genBranch);
  task->SetRecBckgBranch(recBckg);
  task->SetGenBckgBranch(genBckg);
  task->SetEventSelectionMask(fEvtSelMask);
  task->SetEventClass(fEvtClass);
  task->SetTrackFilter(nTrackFilter);
  task->SetVertexCuts(nVtxContCut, fVtxZcut, fVtxRcut);
  task->SetTrackCuts(trackPtMin, trackPtMax, trackEtaAbs);
  task->SetJetCuts(fJetPtMin, fJetEta, fJetZ);

  if(recBranch.Contains("MC") || genBranch.Contains("MC")) task->UseMC(kTRUE);
  else task->UseMC(kFALSE);

  if(recBranch.Contains("MC2")) task->SetTrackType(3);
  else if(recBranch.Contains("MC")) task->SetTrackType(2);
  else task->SetTrackType(1);
  
  if(genBranch.Contains("MC2")) task->SetParticleType(3);
  else if(genBranch.Contains("MC")) task->SetParticleType(2);
  else task->SetParticleType(1);

  mgr->AddTask(task);

  AliAnalysisDataContainer *container = mgr->CreateContainer(name.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:JETANALYSIS_%s", AliAnalysisManager::GetCommonFileName(), name.Data()));

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(task, 1, container);

  mgr->Print();
  cout<<__LINE__<<endl;
  
  return task;
}
