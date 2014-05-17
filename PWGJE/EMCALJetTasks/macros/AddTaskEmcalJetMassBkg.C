AliAnalysisTaskEmcalJetMassBkg* AddTaskEmcalJetMassBkg(const char * njetsBase, 
						       const Double_t R,
						       const char * nrhoBase, 
						       const char * nrhoMass,
						       const char * ntracks, 
						       const char * nclusters,	
						       const char *type,					     
						       const char *CentEst,
						       Int_t       pSel,
						       TString     trigClass,
						       TString     kEmcalTriggers) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalJetMassBkg","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalJetMassBkg", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName = Form("JetMassBkg_%s_TC%s",njetsBase,trigClass.Data());

  //Configure jet tagger task
  AliAnalysisTaskEmcalJetMassBkg *task = new AliAnalysisTaskEmcalJetMassBkg(wagonName.Data());

  task->SetNCentBins(4);
  task->SetConeRadius(R);

  if (strcmp(type,"TPC")==0)
    task->SetConeEtaPhiTPC();
  else if (strcmp(type,"EMCAL")==0)
    task->SetConeEtaPhiEMCAL();

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);

  TString strJetsBase(njetsBase);
  if(strJetsBase.Contains("JetPythia")) {
    if(trackCont)   trackCont->SetTrackBitMap(TObject::kBitMask);
    if(clusterCont) clusterCont->SetClusterBitMap(TObject::kBitMask);
  }

  task->SetJetContainerBase(0);

  TString strType(type);
  AliJetContainer *jetContBase = task->AddJetContainer(njetsBase,strType,R);
  if(jetContBase) {
    jetContBase->SetRhoName(nrhoBase);
    jetContBase->SetRhoMassName(nrhoMass);
    jetContBase->ConnectParticleContainer(trackCont);
    jetContBase->ConnectClusterContainer(clusterCont);
    //  jetContBase->SetZLeadingCut(0.98,0.98);
    jetContBase->SetPercAreaCut(0.6);
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
