AliAnalysisTaskJetMassResponseDet* AddTaskJetMassResponseDet(const char * njetsPart,
							     const char * njetsDet,
							     const Double_t R,
							     const char *type,					     
							     Int_t       pSel,
							     TString     kEmcalTriggers = "",
							     TString     tag            = "") {

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

  TString wagonName = Form("JetMassResponseDet_%s%s",njetsDet,tag.Data());
  TString strType(type);

  //Configure jet mass detector response task
  AliAnalysisTaskJetMassResponseDet *task = new AliAnalysisTaskJetMassResponseDet(wagonName.Data());

  task->SetNCentBins(1);
  //task->SetVzRange(-10.,10.);

  task->SetJetContainerPart(0);
  task->SetJetContainerDet(1);

  AliJetContainer *jetContPart = task->AddJetContainer(njetsPart,strType,R);
  if(jetContPart) {
    // jetContPart->SetPercAreaCut(0.6);
  }

  AliJetContainer *jetContDet = task->AddJetContainer(njetsDet,strType,R);
  if(jetContDet) {
    jetContDet->SetPercAreaCut(0.6);
  }

  task->SetCaloTriggerPatchInfoName(kEmcalTriggers.Data());
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

