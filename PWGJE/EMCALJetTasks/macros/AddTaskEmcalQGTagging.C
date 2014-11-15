AliAnalysisTaskEmcalQGTagging* AddTaskEmcalQGTagging(const char * njetsBase,
						     const char * njetsTrue,
						     const Double_t R,
						     const char * nrhoBase, 
						     const char * ntracks, 
						     const char * nclusters,
						     const char * ntracksTrue,
						     const char *type,				      
						     const char *CentEst,
						     Int_t       pSel,
						     TString     trigClass      = "",
						     TString     kEmcalTriggers = "",
						     TString     tag            = "",
						     AliAnalysisTaskEmcalQGTagging::JetShapeType jetShapeType, AliAnalysisTaskEmcalQGTagging::JetShapeSub jetShapeSub ) {
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskEmcalQGTagging","No analysis manager found.");
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

  TString wagonName = Form("JetQGTaggings_%s_TC%s%s",njetsBase,trigClass.Data(),tag.Data());

  //Configure jet tagger task
  AliAnalysisTaskEmcalQGTagging *task = new AliAnalysisTaskEmcalQGTagging(wagonName.Data());

  //task->SetNCentBins(4);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  
  TString thename(njetsBase);
  //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont  = task->AddParticleContainer(ntracks);
  AliParticleContainer *trackContTrue  = task->AddParticleContainer(ntracksTrue);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);

  AliJetContainer *jetContBase=0x0;
  AliJetContainer *jetContTrue=0x0;

  TString strType(type);

  if (jetShapeType==AliAnalysisTaskEmcalQGTagging::kTrue) {
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      jetContBase->SetPartonInfoName("PartonsInfo");
    }
  }
  
  if (jetShapeType==AliAnalysisTaskEmcalQGTagging::kTrueDet){
    
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      jetContBase->SetPartonInfoName("PartonsInfo");
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(0.6); 
      jetContTrue->SetPartonInfoName("PartonsInfo");
    }
  }  

  if (jetShapeType==AliAnalysisTaskEmcalQGTagging::kData){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
    }    
  }
  
  if (jetShapeType==AliAnalysisTaskEmcalQGTagging::kDetEmb){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      jetContBase->SetPartonInfoName("PartonsInfo");
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(0.6); 
      jetContTrue->SetPartonInfoName("PartonsInfo");
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

  if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kTrue) contName1 += "_True"; 
  if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kTrueDet) contName1 += "_TrueDet"; 
  if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kData) contName1 += "_Data"; 
  if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kDetEmb) contName1 += "_DetEmb"; 
 
  if (jetShapeSub == AliAnalysisTaskEmcalQGTagging::kNoSub) contName1 += "_NoSub"; 
  if (jetShapeSub == AliAnalysisTaskEmcalQGTagging::kConstSub) contName1 += "_ConstSub"; 
  if (jetShapeSub == AliAnalysisTaskEmcalQGTagging::kDerivSub) contName1 += "_DerivSub"; 


  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
    
  mgr->ConnectOutput(task,1,coutput1);

  return task;  

}

