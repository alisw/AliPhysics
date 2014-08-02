AliAnalysisTaskJetShapeGR *AddTaskJetShapeGR(const char * njetsBase,
					     const char * njetsSub,
					     const char * njetsTrue,
					     const Double_t R,
					     const char * nrhoBase,
					     const char * nrhoMass,
					     const char * ntracks,
					     const char * nclusters,
					     const char * ntracksTrue,
					     const char * type           = "TPC",
					     const char * CentEst        = "V0M",
					     Int_t        pSel           = AliVEvent::kAny,
					     TString      trigClass      = "",
					     TString      kEmcalTriggers = "",
					     TString      tag            = "MCMatch")
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskJetShapeGR","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskJetShapeGR", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName = Form("JetShapeGR_%s_TC%s%s",njetsBase,trigClass.Data(),tag.Data());

  //Configure jet tagger task
  AliAnalysisTaskJetShapeGR *task = new AliAnalysisTaskJetShapeGR(wagonName.Data());

  task->SetNCentBins(4);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
  AliParticleContainer *trackContTrue  = task->AddParticleContainer(ntracksTrue);

  Int_t ic = 0;

  TString strType(type);
  AliJetContainer *jetContBase = task->AddJetContainer(njetsBase,strType,R);
  if(jetContBase) {
    jetContBase->SetRhoName(nrhoBase);
    jetContBase->SetRhoMassName(nrhoMass);
    jetContBase->ConnectParticleContainer(trackCont);
    jetContBase->ConnectClusterContainer(clusterCont);
    jetContBase->SetPercAreaCut(0.6);
    task->SetJetContainerBase(ic);
    ic++;
  }
  AliJetContainer *jetContSub = task->AddJetContainer(njetsSub,strType,R);
  if(jetContSub) {
    jetContSub->SetRhoName(nrhoBase);
    jetContSub->SetRhoMassName(nrhoMass);
    jetContSub->ConnectParticleContainer(trackCont);
    jetContSub->ConnectClusterContainer(clusterCont);
    jetContSub->SetPercAreaCut(0.6);
    jetContSub->SetJetPtCut(-1e6);
    task->SetJetContainerSub(ic);
    ic++;
  }
  AliJetContainer *jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
  if(jetContTrue) {
    jetContTrue->ConnectParticleContainer(trackContTrue);
    jetContTrue->SetPercAreaCut(0.6);
    jetContTrue->SetJetPtCut(0.1);
    task->SetJetContainerTrue(ic);
    ic++;
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
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("%sTree",contName.Data()), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);

  return task;



}
