AliAnalysisTaskJetShapeDeriv *AddTaskJetShapeDeriv(const char * njetsBase,
						   const char * njetsNoEmb,
						   const Double_t R,
						   const char * nrhoBase,
						   const char * nrhoMass,
						   const char * ntracks,
						   const char * nclusters,
						   const char * type           = "TPC",
						   const char * CentEst        = "V0M",
						   Int_t        pSel           = AliVEvent::kAny,
						   TString      trigClass      = "",
						   TString      kEmcalTriggers = "",
						   Bool_t       bCreateTree    = kFALSE,
						   Bool_t       removeoverlap  = kFALSE,
						   const char * njetsOverl     = "",
						   const char * ntmptracksOvlJ = "",
						   Double_t     sigJetpTCut    = 5.,
						   TString      tag            = "MCMatch")
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskJetShapeDeriv","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskJetShapeDeriv", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName = Form("JetShapeDeriv_%s_TC%s%s%s",njetsBase,trigClass.Data(),tag.Data(), removeoverlap ? "noOvl" : "");

  //Configure jet tagger task
  AliAnalysisTaskJetShapeDeriv *task = new AliAnalysisTaskJetShapeDeriv(wagonName.Data());

  task->SetNCentBins(4);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliParticleContainer *trackContO = task->AddParticleContainer(ntmptracksOvlJ);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);

  Int_t contindx = 0;
  TString strType(type);
  AliJetContainer *jetContBase = task->AddJetContainer(njetsBase,strType,R);
  if(jetContBase) {
    jetContBase->SetRhoName(nrhoBase);
    jetContBase->SetRhoMassName(nrhoMass);
    jetContBase->ConnectParticleContainer(trackCont);
    jetContBase->ConnectClusterContainer(clusterCont);
    jetContBase->SetPercAreaCut(0.6);
    task->SetJetContainerBase(contindx);
    contindx+=1;
  }

  AliJetContainer *jetContNoEmb = task->AddJetContainer(njetsNoEmb,strType,R);
  if(jetContNoEmb) {
    jetContNoEmb->SetRhoName(nrhoBase);
    jetContNoEmb->SetRhoMassName(nrhoMass);
    jetContNoEmb->ConnectParticleContainer(trackCont);
    jetContNoEmb->ConnectClusterContainer(clusterCont);
    jetContNoEmb->SetPercAreaCut(0.6);
    jetContNoEmb->SetJetPtCut(-1e6);
    task->SetJetContainerNoEmb(contindx);
    contindx+=1;
  }

  AliJetContainer *jetContOverl = 0x0;
  if(removeoverlap) {
     jetContOverl = task->AddJetContainer(njetsOverl,strType,R);
     task->SetRemoveOverlapTrackJet(kTRUE, R);
     if(jetContOverl) {
     	jetContOverl->SetRhoName(nrhoBase);
     	jetContOverl->SetRhoMassName(nrhoMass);
     	jetContOverl->ConnectParticleContainer(trackContO);
     	jetContOverl->ConnectClusterContainer(clusterCont);
     	jetContOverl->SetPercAreaCut(0.6);
     	jetContOverl->SetJetPtCut(sigJetpTCut);
     	task->SetJetContainerOverlap(contindx);
     	contindx+=1;
     }
  }
  
  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
  task->SetCentralityEstimator(CentEst);
  task->SelectCollisionCandidates(pSel);
  task->SetUseAliAnaUtils(kFALSE);
  task->SetCreateTree(bCreateTree);

  mgr->AddTask(task);

  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName(wagonName);
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);
  if(bCreateTree) {
    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(Form("%sTree",contName.Data()), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    mgr->ConnectOutput(task,2,coutput2);
  }

  return task;
}
