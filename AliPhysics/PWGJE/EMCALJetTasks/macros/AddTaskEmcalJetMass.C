AliAnalysisTaskEmcalJetMass* AddTaskEmcalJetMass(const char * njetsBase,
						 const Double_t R,
						 const char * nrhoBase,
						 const char * ntracks, 
						 const char * nclusters,
						 const char *type,					     
						 const char *CentEst,
						 Int_t       pSel,
						 TString     trigClass      = "",
						 TString     kEmcalTriggers = "",
						 TString     nJetsUnsub     = "",
						 const char * nrhoMassBase  = "",
						 TString     optionName     = "") {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalJetMass","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalJetMass", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName = Form("JetMass_%s_TC%s%s",njetsBase,trigClass.Data(),optionName.Data());

  //Configure jet tagger task
  AliAnalysisTaskEmcalJetMass *task = new AliAnalysisTaskEmcalJetMass(wagonName.Data());

  task->SetNCentBins(4);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);

  task->SetJetContainerBase(0);

  TString strType(type);
  AliJetContainer *jetContBase = task->AddJetContainer(njetsBase,strType,R);
  if(jetContBase) {
    jetContBase->SetRhoName(nrhoBase);
    jetContBase->SetRhoMassName(nrhoMassBase);
    jetContBase->ConnectParticleContainer(trackCont);
    jetContBase->ConnectClusterContainer(clusterCont);
    jetContBase->SetPercAreaCut(0.6);
  }
  if(!nJetsUnsub.IsNull()) {
    AliJetContainer *jetContUS = task->AddJetContainer(nJetsUnsub.Data(),strType,R);
    if(jetContUS) {
      jetContUS->SetRhoName(nrhoBase);
      jetContUS->SetRhoMassName(nrhoMassBase);
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

