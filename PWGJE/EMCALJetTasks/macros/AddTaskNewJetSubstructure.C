AliAnalysisTaskNewJetSubstructure* AddTaskNewJetSubstructure(const char * njetsBase,
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
						     AliAnalysisTaskNewJetSubstructure::JetShapeType jetShapeType = AliAnalysisTaskNewJetSubstructure::kMCTrue,
						     AliAnalysisTaskNewJetSubstructure::JetShapeSub jetShapeSub = AliAnalysisTaskNewJetSubstructure::kNoSub,
						     AliAnalysisTaskNewJetSubstructure::JetSelectionType jetSelection = AliAnalysisTaskNewJetSubstructure::kInclusive,
						     Float_t minpTHTrigger =0.,  Float_t maxpTHTrigger =0., Float_t acut =0.6, AliAnalysisTaskNewJetSubstructure::DerivSubtrOrder derivSubtrOrder = AliAnalysisTaskNewJetSubstructure::kSecondOrder ) {
 

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskNewJetSubstructure","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskNewJetSubstructure", "This task requires an input event handler");
      return NULL;
    }

  TString wagonName1 = Form("JetSubstructure_%s_TC%s%s",njetsBase,trigClass.Data(),tag.Data());
  TString wagonName2 = Form("JetSubstructure_%s_TC%s%sTree",njetsBase,trigClass.Data(),tag.Data());
  //Configure jet tagger task
  AliAnalysisTaskNewJetSubstructure *task = new AliAnalysisTaskNewJetSubstructure(wagonName1.Data());

  //task->SetNCentBins(4);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);
  
 

  TString thename(njetsBase);
  //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont;// = task->AddTrackContainer(ntracks);
 
  if ((jetShapeSub==AliAnalysisTaskNewJetSubstructure::kConstSub || jetShapeSub==AliAnalysisTaskNewJetSubstructure::kEventSub ) && ((jetShapeType==AliAnalysisTaskNewJetSubstructure::kData) || (jetShapeType==AliAnalysisTaskNewJetSubstructure::kDetEmbPartPythia) || (jetShapeType==AliAnalysisTaskNewJetSubstructure::kPythiaDef))){
    trackCont = task->AddParticleContainer(ntracks);}
  else trackCont = task->AddTrackContainer(ntracks);

  
  //Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer *trackContUS  = task->AddTrackContainer(ntracksUS);
  //Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer *trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  //Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
   if (jetShapeType==AliAnalysisTaskNewJetSubstructure::kDetEmbPartPythia) trackContTrue->SetIsEmbedding(true);
  AliParticleContainer *trackContPartLevel=0;
  
  if ((jetShapeSub==AliAnalysisTaskNewJetSubstructure::kConstSub) && ((jetShapeType==AliAnalysisTaskNewJetSubstructure::kMCTrue) || (jetShapeType==AliAnalysisTaskNewJetSubstructure::kPythiaDef))){
    trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
  }
  else trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
    if (jetShapeType==AliAnalysisTaskNewJetSubstructure::kDetEmbPartPythia) trackContPartLevel->SetIsEmbedding(true);
   //Printf("ntracksPartLevel() = %s, trackContPartLevel=%p ", ntracksPartLevel, trackContPartLevel);
  

  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
  
  AliJetContainer *jetContBase=0x0;
  AliJetContainer *jetContUS=0x0;
  AliJetContainer *jetContTrue=0x0;
  AliJetContainer *jetContPart=0x0;
  TString strType(type);

  if ((jetShapeType==AliAnalysisTaskNewJetSubstructure::kMCTrue || (jetShapeType==AliAnalysisTaskNewJetSubstructure::kGenOnTheFly))) {
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackContPartLevel);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }
  }
  
  if (jetShapeType==AliAnalysisTaskNewJetSubstructure::kData){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
      if(jetShapeSub==AliAnalysisTaskNewJetSubstructure::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }    
  }
  

  if (jetShapeType==AliAnalysisTaskNewJetSubstructure::kDetEmbPartPythia){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
     
      if(jetShapeSub==AliAnalysisTaskNewJetSubstructure::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut); 
    
    }
    
    if(jetShapeSub==AliAnalysisTaskNewJetSubstructure::kConstSub || jetShapeSub==AliAnalysisTaskNewJetSubstructure::kEventSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(acut);
       
      }
    }
 
     jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
      if(jetContPart) {

        jetContPart->SetRhoName(nrhoBase);
        jetContPart->ConnectParticleContainer(trackContPartLevel);
        jetContPart->SetPercAreaCut(acut);
        
      }
  }
  
  if (jetShapeType==AliAnalysisTaskNewJetSubstructure::kPythiaDef){
    
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }
    
    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);
      
    }
    
    if(jetShapeSub==AliAnalysisTaskNewJetSubstructure::kConstSub){
      jetContUS=task->AddJetContainer(njetsUS,strType,R);
      if(jetContUS) {
        jetContUS->SetRhoName(nrhoBase);
        jetContUS->ConnectParticleContainer(trackContUS);
        jetContUS->SetPercAreaCut(acut);
        
      }
    }
    
    jetContPart = task->AddJetContainer(njetsPartLevel,strType,R);
    if(jetContPart) {
      jetContPart->SetRhoName(nrhoBase);
      jetContPart->ConnectParticleContainer(trackContPartLevel);
      jetContPart->SetPercAreaCut(acut);
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

  if (jetShapeType == AliAnalysisTaskNewJetSubstructure::kMCTrue) contName1 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskNewJetSubstructure::kData) contName1 += "_Data"; 
 
  if (jetShapeType == AliAnalysisTaskNewJetSubstructure::kPythiaDef) contName1 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskNewJetSubstructure::kNoSub) contName1 += "_NoSub"; 
  if (jetShapeSub == AliAnalysisTaskNewJetSubstructure::kConstSub) contName1 += "_ConstSub";
   if (jetShapeSub == AliAnalysisTaskNewJetSubstructure::kEventSub) contName1 += "_EventSub"; 
  if (jetShapeSub == AliAnalysisTaskNewJetSubstructure::kDerivSub) contName1 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskNewJetSubstructure::kInclusive) contName1 += "_Incl";
 


    if (jetShapeType == AliAnalysisTaskNewJetSubstructure::kMCTrue) contName2 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskNewJetSubstructure::kData) contName2 += "_Data"; 
 
  if (jetShapeType == AliAnalysisTaskNewJetSubstructure::kPythiaDef) contName2 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskNewJetSubstructure::kNoSub) contName2 += "_NoSub"; 
  if (jetShapeSub == AliAnalysisTaskNewJetSubstructure::kConstSub) contName2 += "_ConstSub";
  if (jetShapeSub == AliAnalysisTaskNewJetSubstructure::kEventSub) contName2 += "_EventSub"; 
  if (jetShapeSub == AliAnalysisTaskNewJetSubstructure::kDerivSub) contName2 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskNewJetSubstructure::kInclusive) contName2 += "_Incl";




  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   mgr->ConnectOutput(task,1,coutput1);


   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);
   
  return task;  

}

