AliAnalysisTaskEmcalQGTagging* AddTaskEmcalQGTagging(const char * njetsBase,
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
						     AliAnalysisTaskEmcalQGTagging::JetShapeType jetShapeType,
						     AliAnalysisTaskEmcalQGTagging::JetShapeSub jetShapeSub,
						     AliAnalysisTaskEmcalQGTagging::JetSelectionType jetSelection,
                 Float_t minpTHTrigger =0.,  Float_t maxpTHTrigger =0., AliAnalysisTaskEmcalQGTagging::DerivSubtrOrder derivSubtrOrder = 0 ) {
 

  
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

  TString wagonName1 = Form("JetQGTaggings_%s_TC%s%s",njetsBase,trigClass.Data(),tag.Data());
  TString wagonName2 = Form("JetQGTaggings_%s_TC%s%sTree",njetsBase,trigClass.Data(),tag.Data());
  //Configure jet tagger task
  AliAnalysisTaskEmcalQGTagging *task = new AliAnalysisTaskEmcalQGTagging(wagonName1.Data());

  //task->SetNCentBins(4);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);
  
  if (jetSelection == AliAnalysisTaskEmcalQGTagging::kRecoil) task->SetPtTriggerSelections(minpTHTrigger, maxpTHTrigger);

  TString thename(njetsBase);
  //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont;// = task->AddTrackContainer(ntracks);
 
  if ((jetShapeSub==AliAnalysisTaskEmcalQGTagging::kConstSub) && ((jetShapeType==AliAnalysisTaskEmcalQGTagging::kData) || (jetShapeType==AliAnalysisTaskEmcalQGTagging::kDetEmbPartPythia) || (jetShapeType==AliAnalysisTaskEmcalQGTagging::kPythiaDef))){
    trackCont = task->AddParticleContainer(ntracks);}
  else trackCont = task->AddTrackContainer(ntracks);

  
  //Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer *trackContUS  = task->AddTrackContainer(ntracksUS);
  //Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer *trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  //Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
  
  AliParticleContainer *trackContPartLevel=0;
  
  if ((jetShapeSub==AliAnalysisTaskEmcalQGTagging::kConstSub) && ((jetShapeType==AliAnalysisTaskEmcalQGTagging::kMCTrue) || (jetShapeType==AliAnalysisTaskEmcalQGTagging::kPythiaDef))){
    trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
  }
  else trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
  
   //Printf("ntracksPartLevel() = %s, trackContPartLevel=%p ", ntracksPartLevel, trackContPartLevel);
  

  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
  
  AliJetContainer *jetContBase=0x0;
  AliJetContainer *jetContUS=0x0;
  AliJetContainer *jetContTrue=0x0;
  AliJetContainer *jetContPart=0x0;
  TString strType(type);

  if ((jetShapeType==AliAnalysisTaskEmcalQGTagging::kMCTrue || (jetShapeType==AliAnalysisTaskEmcalQGTagging::kGenOnTheFly))) {
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackContPartLevel);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
    }
  }
  
  if (jetShapeType==AliAnalysisTaskEmcalQGTagging::kData){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
      if(jetShapeSub==AliAnalysisTaskEmcalQGTagging::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }    
  }
  

  if (jetShapeType==AliAnalysisTaskEmcalQGTagging::kDetEmbPartPythia){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
     
      if(jetShapeSub==AliAnalysisTaskEmcalQGTagging::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(0.6); 
    
    }
    
    if(jetShapeSub==AliAnalysisTaskEmcalQGTagging::kConstSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(0.6);
       
      }
    }
 
     jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
      if(jetContPart) {
        jetContPart->SetRhoName(nrhoBase);
        jetContPart->ConnectParticleContainer(trackContPartLevel);
        jetContPart->SetPercAreaCut(0.6);
        
      }
  }
  
  if (jetShapeType==AliAnalysisTaskEmcalQGTagging::kPythiaDef){
    
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(0.6);
    }
    
    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(0.6);
      
    }
    
    if(jetShapeSub==AliAnalysisTaskEmcalQGTagging::kConstSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(0.6);
        
      }
    }
    
    jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
    if(jetContPart) {
      jetContPart->SetRhoName(nrhoBase);
      jetContPart->ConnectParticleContainer(trackContPartLevel);
      jetContPart->SetPercAreaCut(0.6);
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

  if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kMCTrue) contName1 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kData) contName1 += "_Data"; 
 
  if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kPythiaDef) contName1 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskEmcalQGTagging::kNoSub) contName1 += "_NoSub"; 
  if (jetShapeSub == AliAnalysisTaskEmcalQGTagging::kConstSub) contName1 += "_ConstSub"; 
  if (jetShapeSub == AliAnalysisTaskEmcalQGTagging::kDerivSub) contName1 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskEmcalQGTagging::kInclusive) contName1 += "_Incl";
  if (jetSelection == AliAnalysisTaskEmcalQGTagging::kRecoil) {

  TString recoilTriggerString = Form("_Recoil_%.0f_%0.f", minpTHTrigger, maxpTHTrigger);
  contName1 += recoilTriggerString;
  }


    if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kMCTrue) contName2 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kData) contName2 += "_Data"; 
 
  if (jetShapeType == AliAnalysisTaskEmcalQGTagging::kPythiaDef) contName2 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskEmcalQGTagging::kNoSub) contName2 += "_NoSub"; 
  if (jetShapeSub == AliAnalysisTaskEmcalQGTagging::kConstSub) contName2 += "_ConstSub"; 
  if (jetShapeSub == AliAnalysisTaskEmcalQGTagging::kDerivSub) contName2 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskEmcalQGTagging::kInclusive) contName2 += "_Incl";
  if (jetSelection == AliAnalysisTaskEmcalQGTagging::kRecoil) {

  TString recoilTriggerString = Form("_Recoil_%.0f_%0.f", minpTHTrigger, maxpTHTrigger);
  contName2 += recoilTriggerString;
  }




  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   mgr->ConnectOutput(task,1,coutput1);


   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);
   
  return task;  

}

