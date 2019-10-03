AliAnalysisTaskRhoMassScale* AddTaskRhoMassScale(const char *rhomNeutral,
						 const char *rhomChEmcal,
						 const char *rhomCh2xEmcal,
						 const char *ntracks,
						 const Double_t kMinTrackPt = 0.15,
						 Int_t       pSel        = AliVEvent::kAny,
						 const char * njetsNeutral  = "",
						 const char * njetsCharged   = "",
						 const Double_t R        = 0.4,
						 const char *type        = "EMCAL",					     
						 TString     tag         = "") {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskRhoMassScale","No analysis manager found.");
      return 0;
    }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskRhoMassScale", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName = Form("RhoMassScale%s",tag.Data());
  TString strType(type);

  //Configure jet mass detector response task
  AliAnalysisTaskRhoMassScale *task = new AliAnalysisTaskRhoMassScale(wagonName.Data());

  task->SetNCentBins(1);
  //task->SetVzRange(-10.,10.);

  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);

  task->SetRhoMNeutralName(rhomNeutral);
  task->SetRhoMChargedEmcalName(rhomChEmcal);
  task->SetRhoMCharged2xEmcalName(rhomCh2xEmcal);

  AliParticleContainer *pcont = task->AddParticleContainer(ntracks);
  pcont->SetParticlePtCut(kMinTrackPt);
  pcont->SetParticleEtaLimits(-0.9,0.9);
  pcont->SetCharge(AliParticleContainer::kPositiveCharge);

  task->SetJetContainerNeutral(0);
  // task->SetJetContainerCharged(1);

  AliJetContainer *jetContNeutral = task->AddJetContainer(njetsNeutral,strType,R);
  if(jetContNeutral) {
    jetContNeutral->SetPercAreaCut(0.6);
  }

  // AliJetContainer *jetContCharged = task->AddJetContainer(njetsCharged,strType,R);
  // if(jetContCharged) {
  //   jetContCharged->SetPercAreaCut(0.6);
  // }

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName(wagonName);
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);

  return task;  
}

