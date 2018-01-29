AliAnalysisTaskRecursiveSoftDrop* AddTaskRecursiveSoftDrop(                                 const char * njetsData, //data jets
								    const Double_t R,
								    const char * nrhoBase, 
								    const char * ntracksData,
								    const char *type,				      
								    const char *CentEst,
								    Int_t       pSel,
								    TString     trigClass      = "",
								    TString     kEmcalTriggers = "",
								    TString     tag            = "",
								    AliAnalysisTaskRecursiveSoftDrop::JetShapeSub jetShapeSub
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

    TString wagonName1 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%s",njetsData,trigClass.Data(),tag.Data());
    TString wagonName2 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_Det",njetsData,trigClass.Data(),tag.Data());
    TString wagonName3 = Form("AliAnalysisTaskRecursiveSoftDrop_%s_TC%s%sTree_True",njetsData,trigClass.Data(),tag.Data());

  //Configure jet tagger task
  AliAnalysisTaskRecursiveSoftDrop *task = new AliAnalysisTaskRecursiveSoftDrop(wagonName1.Data());

  task->SetJetShapeSub(jetShapeSub);
  task->SetJetRadius(R);
  
  AliParticleContainer *trackContData=0x0; 
  trackContData = task->AddParticleContainer(ntracksData);
 
  AliJetContainer *JetContData=0x0;

  TString strType(type);



    JetContData = task->AddJetContainer(njetsData,strType,R); //Data
    if(JetContData) {
      JetContData->SetRhoName(nrhoBase);
      JetContData->ConnectParticleContainer(trackContData);
      JetContData->SetPercAreaCut(0.6);
      JetContData->SetJetRadius(R);
      JetContData->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecursiveSoftDrop::kConstSub) JetContData->SetAreaEmcCut(-2);
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



 
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(contName3.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,3,coutput3);

  return task;  

}

