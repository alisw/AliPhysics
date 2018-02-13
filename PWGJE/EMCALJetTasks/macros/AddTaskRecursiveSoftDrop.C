AliAnalysisTaskRecursiveSoftDrop* AddTaskRecursiveSoftDrop(         const char * njetsData, //data jets
								    const char * njetsTrue, //Pyhthia Particle Level
								    const char * njetsDet,
								    const char * njetsHybridUs,
								    const char * njetsHybridS,
								    const Double_t R,
								    const char * nrhoBase, 
								    const char * ntracksData,
                                                                    const char * ntracksTrue,
                                                                    const char * ntracksDet, 
								    const char * ntracksHybridUs,
								    const char * ntracksHybridS,
								    const char *type,				      
								    const char *CentEst,
								    Int_t       pSel,
								    TString     trigClass      = "",
								    TString     kEmcalTriggers = "",
								    TString     tag            = "",
								    AliAnalysisTaskRecursiveSoftDrop::JetShapeSub jetShapeSub,
								    AliAnalysisTaskRecursiveSoftDrop::JetShapeSub fjetType
								    ) {
  
  
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskRecursiveSoftDrop","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AliAnalysisTaskRecursiveSoftDrop", "This task requires an input event handler");
      return NULL;
    }

  if(fJetType==AliAnalysisTaskRecursiveSoftDrop::kData){
    TString wagonName1 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%s",njetsData,trigClass.Data(),tag.Data());
    TString wagonName2 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_Det",njetsData,trigClass.Data(),tag.Data());
    TString wagonName3 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_True",njetsData,trigClass.Data(),tag.Data());
  }
  if(fJetType==AliAnalysisTaskRecursiveSoftDrop::kEmb){
    TString wagonName1 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%s",njetsHybridS,trigClass.Data(),tag.Data());
    TString wagonName2 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_Det",njetsHybridS,trigClass.Data(),tag.Data());
    TString wagonName3 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_True",njetsHybridS,trigClass.Data(),tag.Data());
  }
  //Configure jet tagger task
  AliAnalysisTaskRecursiveSoftDrop *task = new AliAnalysisTaskRecursiveSoftDrop(wagonName1.Data());

  task->SetJetShapeSub(jetShapeSub);
  task->SetJetType(fjetType);
  task->SetJetRadius(R);
  
  AliParticleContainer *trackContData=0x0;
  AliParticleContainer *trackContDet=0x0;
  AliParticleContainer *trackContTrue=0x0;

  trackContData = task->AddParticleContainer(ntracksData);
  if(fJetType==kEmb){
    trackContDet = task->AddParticleContainer(ntracksDet);
    trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  }
  AliJetContainer *JetContData=0x0;
  AliJetContainer *JetContTrue=0x0;
  AliJetContainer *JetContDet=0x0;
  AliJetContainer *JetContHybridUs=0x0;
  AliJetContainer *JetContHybridS=0x0;
  TString strType(type);

  if(fJetType==kData){
    JetContData = task->AddJetContainer(njetsData,strType,R); //Data
    if(JetContData) {
      JetContData->SetRhoName(nrhoBase);
      JetContData->ConnectParticleContainer(trackContData);
      JetContData->SetPercAreaCut(0.6);
      JetContData->SetJetRadius(R);
      JetContData->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContData->SetAreaEmcCut(-2);
    }  
  }
  if(fJetType==kEmb){
    JetContHybridS = task->AddJetContainer(njetsHybridS,strType,R); //HybridS
    if(JetContHybridS) {
      JetContHybridS->SetRhoName(nrhoBase);
      JetContHybridS->ConnectParticleContainer(trackContHybridS);
      JetContHybridS->SetPercAreaCut(0.6);
      JetContHybridS->SetJetRadius(R);
      JetContHybridS->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContHybridS->SetAreaEmcCut(-2);
    }
    JetContHybridUs = task->AddJetContainer(njetsHybridUs,strType,R); //HybridUs
    if(JetContHybridUs) {
      JetContHybridUs->SetRhoName(nrhoBase);
      JetContHybridUs->ConnectParticleContainer(trackContHybridUs);
      JetContHybridUs->SetPercAreaCut(0.6);
      JetContHybridUs->SetJetRadius(R);
      JetContHybridUs->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContHybridUs->SetAreaEmcCut(-2);
    }
    JetContDet = task->AddJetContainer(njetsDet,strType,R); //Det
    if(JetContDet) {
      JetContDet->SetRhoName(nrhoBase);
      JetContDet->ConnectParticleContainer(trackContDet);
      JetContDet->SetPercAreaCut(0.6);
      JetContDet->SetJetRadius(R);
      JetContDet->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContDet->SetAreaEmcCut(-2);
    }
    JetContTrue = task->AddJetContainer(njetsTrue,strType,R); //True
    if(JetContTrue) {
      JetContTrue->SetRhoName(nrhoBase);
      JetContTrue->ConnectParticleContainer(trackContTrue);
      JetContTrue->SetPercAreaCut(0.6);
      JetContTrue->SetJetRadius(R);
      JetContTrue->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContTrue->SetAreaEmcCut(-2);
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
  TString contName1(wagonName1);
  TString contName2(wagonName2);
  TString contName3(wagonName3);

  if (jetShapeType == AliAnalysisTaskRecoilJetYield::kEmb){
    contName1 += "_Embedded";
    contName2 += "_Embedded";
    contName2 += "_Embedded";

  }

 
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(contName3.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,3,coutput3);

  return task;  

}

