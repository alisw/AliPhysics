EmcalHJetMassAnalysis::AliAnalysisTaskEmcalHJetMass* AddTaskEmcalHJetMass(const char * njetsBase,
						   const Double_t R,
						   const char * nrhoBase,
						   const char * ntracks,
						   const char * nclusters,
						   const char * type,
						   Double_t     ptMinH         = 5.,
						   Double_t     maxDPhi        = 0.6,
						   Int_t        pSel,
						   const char * CentEst        = "V0M",
						   TString      trigClass      = "",
						   TString      kEmcalTriggers = "",
						   TString      nJetsUnsub     = "",
						   TString      tag            = "") {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalHJetMass","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalHJetMass", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName = Form("HJetMass_%s_TC%s%s",njetsBase,trigClass.Data(),tag.Data());

  //Configure jet tagger task
  EmcalHJetMassAnalysis::AliAnalysisTaskEmcalHJetMass *task = new EmcalHJetMassAnalysis::AliAnalysisTaskEmcalHJetMass(wagonName.Data());
  task->SetNCentBins(4);
  task->SetMaxDeltaPhi(maxDPhi);
  
  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  trackCont->SetParticlePtCut(ptMinH);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);

  task->SetJetContainerBase(0);

  TString strType(type);
  AliJetContainer *jetContBase = task->AddJetContainer(njetsBase,strType,R);
  if(jetContBase) {
    jetContBase->SetRhoName(nrhoBase);
    jetContBase->ConnectParticleContainer(trackCont);
    jetContBase->ConnectClusterContainer(clusterCont);
    jetContBase->SetPercAreaCut(0.6);
  }
  if(!nJetsUnsub.IsNull()) {
    AliJetContainer *jetContUS = task->AddJetContainer(nJetsUnsub.Data(),strType,R);
    if(jetContUS) {
      jetContUS->SetRhoName(nrhoBase);
      jetContUS->ConnectParticleContainer(trackCont);
      jetContUS->ConnectClusterContainer(clusterCont);
      jetContUS->SetPercAreaCut(0.6);
    }
  }

  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);

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

