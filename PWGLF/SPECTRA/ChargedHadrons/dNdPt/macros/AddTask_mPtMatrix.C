AliAnalysisTask *AddTask_mmarquar_MatrixPbPb(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_mmarquar_MatrixPbPb", "No analysis manager found.");
    return 0;
  }

  // --- Check for MC ---
  AliMCEventHandler  *mcH = static_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());


//   //============= Set Task Name ===================
//   TString taskName=("AliAnalysisTaskMatrix.cxx+");
//   //===============================================
//   //            Load the task
//   gROOT->LoadMacro(taskName.Data());
//   if (gProof){
//     TString taskSO=gSystem->pwd();
//     taskSO+="/";
//     taskSO+=taskName(0,taskName.First('.'))+"_cxx.so";
//     gProof->Exec(Form("gSystem->Load(\"%s\")",taskSO.Data()),kTRUE);
//   }

  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskMatrix *task = new AliAnalysisTaskMatrix("mmarquar_MatrixPbPb");

  // --- Set ESD track Cuts ---

  //gSystem->Load("libANALYSISalice.so");
  //gSystem->Load("libANALYSIS.so");
  //gSystem->Load("libPWGUDbase.so");
  //gSystem->Load("libPWGLFspectra.so");

  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");

  // use dNdPt cuts
  
  Int_t cutMode = 200;
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");
  AliESDtrackCuts* esdTrackCuts = CreatedNdPtTrackCuts(cutMode);

  esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);
  esdTrackCuts->SetMaxChi2PerClusterITS(36.);

  Double_t ptmin = 0.15;
  Double_t ptmax = 1.e10;
  Double_t etamax = 0.3;
  Double_t etamin = -0.3;
  esdTrackCuts->SetPtRange(ptmin, ptmax);
  esdTrackCuts->SetEtaRange(etamin, etamax);

  task->SetAliESDtrackCuts(esdTrackCuts);
  task->SetUseCentrality(1);
  task->SetCutsGenParticle(ptmin, ptmax, etamin, etamax);
  Float_t zvWindow = 30.;
  task->SetMaxVertexZ(zvWindow);
  task->SetCentLimit(10);
  // if mc - use mc
  if (mcH) task->SetUseMC(kTRUE);
  else task->SetUseMC(kFALSE);
  task->SetSubsample(0);

  task->SetTrigger(AliTriggerAnalysis::kAcceptAll);

  mgr->AddTask(task);

  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  //this is the old way!!!
  //AliAnalysisDataContainer *cinput  = (AliAnalysisDataContainer*)mgr->GetContainers()->FindObject("cAUTO_INPUT");

  //            define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("mmarquar_MatrixPbPb", TList::Class(),
	AliAnalysisManager::kOutputContainer,"mmarquar_MatrixPbPb.root");

  //           connect containers
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  1, coutput1);

  return task;
}

