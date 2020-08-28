AliAnalysisTaskFilterHe3 *AddTask_FilterHe3(AliPID::EParticleType ParticleType = AliPID::kTriton)
{
  
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_FilterHe3", "No analysis manager found.");
    return 0;
  }
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask_FilterHe3", "This task requires an input event handler");
    return NULL;
  }  
  
  TString name = "standard";
  //========= Add task to the ANALYSIS manager =====
  AliAnalysisTaskFilterHe3 *task = new AliAnalysisTaskFilterHe3(Form("akalweitTaskFilterHe3_%s", name.Data()));
  task->SetParticleType(ParticleType);
  //task->SelectCollisionCandidates(AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral);

  switch((Int_t) ParticleType) {

    case (Int_t)AliPID::kHe3 : {
      //keep default setters
      break;
    }

    case (Int_t)AliPID::kTriton : {
       task->SetMinNsig(-4.);
       task->SetMaxNsig(4.);
       task->SetPtotPosMin(0.5);
       task->SetPtotNegMin(0.5);
       task->SetPtotPosMax(1.5);
       task->SetPtotNegMax(1.5);
       task->SetMinMass(2.);
       task->SetMaxMass(4.);
      break;
    }
  }

  TString ParticleName;
  switch((Int_t) ParticleType) {
  case (Int_t)AliPID::kHe3 : ParticleName = "He3"; break; 
  case (Int_t)AliPID::kTriton : ParticleName = "Triton"; break; 
  }

  //================================================
  //              data containers
  //================================================

  //dumm output container
  AliAnalysisDataContainer *coutput0 =
    mgr->CreateContainer(Form("akalweit_treeFilter%s",ParticleName.Data()),
			 TTree::Class(),
			 AliAnalysisManager::kExchangeContainer,
			 "akalweit_default");
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("akalweit_filter%s_hist", ParticleName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:Filter%s_%s", ParticleName.Data(), name.Data()));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("akalweit_filter%s_names", ParticleName.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, Form("AnalysisResults.root:Filter%s_List_%s", ParticleName.Data(), name.Data()));

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
  mgr->ConnectOutput(task, 2, coutput2);


 
  return task;
}
