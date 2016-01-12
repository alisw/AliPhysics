AliAnalysisTaskFakeJets* AddTaskEmcalFakeJets(const char * njetsBase,
                                                     const char * njetsUS,
						     const char * njetsTrue,
                                                     const char * njetsPartLevel,
						     const Double_t R,
						     const char * nrhoBase, 
						     const char * ntracks, 
                                                     const char * ntracksUS,
                                                     const char *ntracksPartLevel,
						     const char * nclusters,
						     const char * ntracksTrue,
						     const char *type,				      
						     const char *CentEst,
						     Int_t       pSel,
						     TString     trigClass      = "",
						     TString     kEmcalTriggers = "",
						     TString     tag            = "",
						     AliAnalysisTaskFakeJets::JetShapeType jetShapeType,
						     AliAnalysisTaskFakeJets::JetShapeSub jetShapeSub,
						     AliAnalysisTaskFakeJets::JetSelectionType jetSelection,
                 Float_t minpTHTrigger =0.,  Float_t maxpTHTrigger =0., AliAnalysisTaskFakeJets::DerivSubtrOrder derivSubtrOrder = 0 ) {
 

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalFakeJets","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskEmcalFakeJets", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName = Form("FakeJets_%s_TC%s%s",njetsBase,trigClass.Data(),tag.Data());

  //Configure jet tagger task
  AliAnalysisTaskFakeJets *task = new AliAnalysisTaskFakeJets(wagonName.Data());

  //task->SetNCentBins(4);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);
  
  if (jetSelection == AliAnalysisTaskFakeJets::kRecoil) task->SetPtTriggerSelections(minpTHTrigger, maxpTHTrigger);

  TString thename(njetsBase);
  //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  //Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer *trackContUS  = task->AddParticleContainer(ntracksUS);
  //Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer *trackContTrue  = task->AddParticleContainer(ntracksTrue);
  //Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
  AliParticleContainer *trackContPartLevel  = task->AddParticleContainer(ntracksPartLevel);

  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
  
  AliJetContainer *jetContBase=0x0;
  AliJetContainer *jetContUS=0x0;
  AliJetContainer *jetContTrue=0x0;
  AliJetContainer *jetContPart=0x0;
  TString strType(type);

  if (jetShapeType==AliAnalysisTaskFakeJets::kTrue) {
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      jetContBase->SetPythiaInfoName("PythiaInfo");
    }
  }
  
  if (jetShapeType==AliAnalysisTaskFakeJets::kTrueDet){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      jetContBase->SetPythiaInfoName("PythiaInfo");
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(0.6); 
      jetContTrue->SetPythiaInfoName("PythiaInfo");
    }
  }
  
  if (jetShapeType==AliAnalysisTaskFakeJets::kData){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      if(jetShapeSub==AliAnalysisTaskFakeJets::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }    
  }
  
  if (jetShapeType==AliAnalysisTaskFakeJets::kDetEmb){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      jetContBase->SetPythiaInfoName("PythiaInfo");
      if(jetShapeSub==AliAnalysisTaskFakeJets::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(0.6); 
      jetContTrue->SetPythiaInfoName("PythiaInfo");
          }
    
    if(jetShapeSub==AliAnalysisTaskFakeJets::kConstSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(0.6);
        jetContUS->SetPythiaInfoName("PythiaInfo");
      }
    }
  }

  if ((jetShapeType==AliAnalysisTaskFakeJets::kDetEmbPart)||(jetShapeType==AliAnalysisTaskFakeJets::kDetEmbPartPythia)){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      jetContBase->SetPythiaInfoName("PythiaInfo");
      if(jetShapeSub==AliAnalysisTaskFakeJets::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(0.6); 
      jetContTrue->SetPythiaInfoName("PythiaInfo");
    }
    
    if(jetShapeSub==AliAnalysisTaskFakeJets::kConstSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(0.6);
        jetContUS->SetPythiaInfoName("PythiaInfo");
      }
    }
 
     jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
      if(jetContPart) {
      jetContPart->SetRhoName(nrhoBase);
      jetContPart->ConnectParticleContainer(trackContPartLevel);
      jetContPart->SetPercAreaCut(0.6); 
      jetContPart->SetPythiaInfoName("PythiaInfo");
    }




    }













  if (jetShapeType==AliAnalysisTaskFakeJets::kPythiaDef){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      //jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
    }
    
    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      //  jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(0.6);
     
    }
    
    if(jetShapeSub==AliAnalysisTaskFakeJets::kConstSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
       // jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(0.6);
       
      }
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
  TString contName1(wagonName);

  if (jetShapeType == AliAnalysisTaskFakeJets::kTrue) contName1 += "_True"; 
  if (jetShapeType == AliAnalysisTaskFakeJets::kTrueDet) contName1 += "_TrueDet"; 
  if (jetShapeType == AliAnalysisTaskFakeJets::kData) contName1 += "_Data"; 
  if (jetShapeType == AliAnalysisTaskFakeJets::kDetEmb) contName1 += "_DetEmb"; 
   if (jetShapeType == AliAnalysisTaskFakeJets::kDetEmbPart) contName1 += "_DetEmbPart"; 
  if (jetShapeType == AliAnalysisTaskFakeJets::kPythiaDef) contName1 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskFakeJets::kNoSub) contName1 += "_NoSub"; 
  if (jetShapeSub == AliAnalysisTaskFakeJets::kConstSub) contName1 += "_ConstSub"; 
  if (jetShapeSub == AliAnalysisTaskFakeJets::kDerivSub) contName1 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskFakeJets::kInclusive) contName1 += "_Incl";
  if (jetSelection == AliAnalysisTaskFakeJets::kRecoil) {

  TString recoilTriggerString = Form("_Recoil_%.0f_%0.f", minpTHTrigger, maxpTHTrigger);
  contName1 += recoilTriggerString;
  }



  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    
  mgr->ConnectOutput(task,1,coutput1);

  return task;  

}

