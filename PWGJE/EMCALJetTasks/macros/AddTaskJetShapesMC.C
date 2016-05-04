AliAnalysisTaskEmcalJetShapesMC* AddTaskJetShapesMC(const char * njetsBase,
                                                    const Double_t jetradius,
                                                    const Double_t subjetradius,
                                                    const char *ntracksPartLevel,
                                                    const char *type,
                                                    const char *CentEst,
                                                    Int_t       pSel,
                                                    TString     trigClass      = "",
                                                    TString     kEmcalTriggers = "",
                                                    TString     tag            = "",
                                                    const char *rhoName = "",
                                                    AliAnalysisTaskEmcalJetShapesMC::JetShapeType jetShapeType,
                                                    AliAnalysisTaskEmcalJetShapesMC::JetShapeSub jetShapeSub,
                                                    AliAnalysisTaskEmcalJetShapesMC::JetSelectionType jetSelection,
                                                    Float_t minpTHTrigger =0.,  Float_t maxpTHTrigger =0.,
                                                    AliAnalysisTaskEmcalJetShapesMC::DerivSubtrOrder derivSubtrOrder = 0 ) {
  

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AliAnalysisTaskEmcalJetShapesMC","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalQGTagging", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName1 = Form("JetShapesMC_%s_Histos%s%s",njetsBase,trigClass.Data(),tag.Data());
  TString wagonName2 = Form("JetShapesMC_%s_Tree%s%s",njetsBase,trigClass.Data(),tag.Data());

  //Configure jet tagger task
  AliAnalysisTaskEmcalJetShapesMC *task = new AliAnalysisTaskEmcalJetShapesMC(Form("JetShapesMC_%s", njetsBase));

  //task->SetNCentBins(4);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);
  task->SetJetRadius(jetradius);
  task->SetSubjetRadius(subjetradius);
  
  if (jetSelection == AliAnalysisTaskEmcalJetShapesMC::kRecoil) task->SetPtTriggerSelections(minpTHTrigger, maxpTHTrigger);

  TString thename(njetsBase);

  AliParticleContainer *trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
  
  AliJetContainer *jetContBase=0x0;
  TString strType(type);

  jetContBase = task->AddJetContainer(njetsBase,strType,jetradius);
  if(jetContBase) {
    jetContBase->SetRhoName(rhoName);
    jetContBase->ConnectParticleContainer(trackContPartLevel);
    //jetContBase->ConnectClusterContainer(clusterCont);
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
  TString contName1(wagonName1);
  TString contName2(wagonName2);
  
  if (jetShapeType == AliAnalysisTaskEmcalJetShapesMC::kGenShapes) {
    contName1 += "_GenShapes";
    contName2 += "_GenShapes";
  }
  
  switch (jetShapeSub) {
      
    case AliAnalysisTaskEmcalJetShapesMC::kNoSub:
      contName1 += "_NoSub";
      contName2 += "_NoSub";
      break;
      
    case AliAnalysisTaskEmcalJetShapesMC::kConstSub:
      contName1 += "_ConstSub";
      contName2 += "_ConstSub";
      break;
      
    case AliAnalysisTaskEmcalJetShapesMC::kDerivSub:
      contName1 += "_DerivSub";
      contName2 += "_DerivSub";
      break;


  }
  
  switch (jetSelection) {
    case AliAnalysisTaskEmcalJetShapesMC::kInclusive:
      contName1 += "_Incl";
      contName2 += "_Incl";
      break;


    case AliAnalysisTaskEmcalJetShapesMC::kRecoil:
      TString recoilTriggerString = Form("_Recoil_%.0f_%0.f", minpTHTrigger, maxpTHTrigger);
      contName1 += recoilTriggerString;
      contName2 += recoilTriggerString;

      break;
    
  }

  
  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
 
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(),TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);

  return task;  

}

