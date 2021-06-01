AliAnalysisTaskJetPull* AddTaskJtPull(const char * njetsBase,
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
						     AliAnalysisTaskJetPull::JetShapeType jetShapeType = AliAnalysisTaskJetPull::kMCTrue,
						     AliAnalysisTaskJetPull::JetShapeSub jetShapeSub = AliAnalysisTaskJetPull::kNoSub,
						     AliAnalysisTaskJetPull::JetSelectionType jetSelection = AliAnalysisTaskJetPull::kInclusive,
						     Float_t minpTHTrigger =0.,  Float_t maxpTHTrigger =0., Float_t acut =0.6, AliAnalysisTaskJetPull::DerivSubtrOrder derivSubtrOrder = AliAnalysisTaskJetPull::kSecondOrder ) {
 

  
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
  AliAnalysisTaskJetPull *task = new AliAnalysisTaskJetPull(wagonName1.Data());

  //task->SetNCentBins(4);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);
  
 

  TString thename(njetsBase);
  //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont;// = task->AddTrackContainer(ntracks);
 
  if ((jetShapeSub==AliAnalysisTaskJetPull::kConstSub || jetShapeSub==AliAnalysisTaskJetPull::kEventSub ) && ((jetShapeType==AliAnalysisTaskJetPull::kData) || (jetShapeType==AliAnalysisTaskJetPull::kDetEmbPartPythia) || (jetShapeType==AliAnalysisTaskJetPull::kPythiaDef))){
    trackCont = task->AddParticleContainer(ntracks);}
  else trackCont = task->AddTrackContainer(ntracks);

  
  //Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer *trackContUS  = task->AddTrackContainer(ntracksUS);
  //Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer *trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  //Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
   if (jetShapeType==AliAnalysisTaskJetPull::kDetEmbPartPythia) trackContTrue->SetIsEmbedding(true);
  AliParticleContainer *trackContPartLevel=0;
  
  if ((jetShapeSub==AliAnalysisTaskJetPull::kConstSub) && ((jetShapeType==AliAnalysisTaskJetPull::kMCTrue) || (jetShapeType==AliAnalysisTaskJetPull::kPythiaDef))){
    trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
  }
  else trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
    if (jetShapeType==AliAnalysisTaskJetPull::kDetEmbPartPythia) trackContPartLevel->SetIsEmbedding(true);
   //Printf("ntracksPartLevel() = %s, trackContPartLevel=%p ", ntracksPartLevel, trackContPartLevel);
  

  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
  
  AliJetContainer *jetContBase=0x0;
  AliJetContainer *jetContUS=0x0;
  AliJetContainer *jetContTrue=0x0;
  AliJetContainer *jetContPart=0x0;
  TString strType(type);

  if ((jetShapeType==AliAnalysisTaskJetPull::kMCTrue || (jetShapeType==AliAnalysisTaskJetPull::kGenOnTheFly))) {
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackContPartLevel);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }
  }
  
  if (jetShapeType==AliAnalysisTaskJetPull::kData){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
      if(jetShapeSub==AliAnalysisTaskJetPull::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }    
  }
  

  if (jetShapeType==AliAnalysisTaskJetPull::kDetEmbPartPythia){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
     
      if(jetShapeSub==AliAnalysisTaskJetPull::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut); 
    
    }
    
    if(jetShapeSub==AliAnalysisTaskJetPull::kConstSub || jetShapeSub==AliAnalysisTaskJetPull::kEventSub){
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
  
  if (jetShapeType==AliAnalysisTaskJetPull::kPythiaDef){
    
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
    
    if(jetShapeSub==AliAnalysisTaskJetPull::kConstSub){
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

  if (jetShapeType == AliAnalysisTaskJetPull::kMCTrue) contName1 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskJetPull::kData) contName1 += "_Data"; 
 
  if (jetShapeType == AliAnalysisTaskJetPull::kPythiaDef) contName1 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskJetPull::kNoSub) contName1 += "_NoSub"; 
  if (jetShapeSub == AliAnalysisTaskJetPull::kConstSub) contName1 += "_ConstSub";
   if (jetShapeSub == AliAnalysisTaskJetPull::kEventSub) contName1 += "_EventSub"; 
  if (jetShapeSub == AliAnalysisTaskJetPull::kDerivSub) contName1 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskJetPull::kInclusive) contName1 += "_Incl";
 


    if (jetShapeType == AliAnalysisTaskJetPull::kMCTrue) contName2 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskJetPull::kData) contName2 += "_Data"; 
 
  if (jetShapeType == AliAnalysisTaskJetPull::kPythiaDef) contName2 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskJetPull::kNoSub) contName2 += "_NoSub"; 
  if (jetShapeSub == AliAnalysisTaskJetPull::kConstSub) contName2 += "_ConstSub";
  if (jetShapeSub == AliAnalysisTaskJetPull::kEventSub) contName2 += "_EventSub"; 
  if (jetShapeSub == AliAnalysisTaskJetPull::kDerivSub) contName2 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskJetPull::kInclusive) contName2 += "_Incl";




  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   mgr->ConnectOutput(task,1,coutput1);


   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);
   
  return task;  

}

