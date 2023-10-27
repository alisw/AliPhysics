AliAnalysisTaskJetsEEC* AddTaskJetsEEC(const char * njetsBase,
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
                             AliAnalysisTaskJetsEEC::JetShapeType jetShapeType = AliAnalysisTaskJetsEEC::kMCTrue,
                             AliAnalysisTaskJetsEEC::JetShapeSub jetShapeSub = AliAnalysisTaskJetsEEC::kNoSub,
                             AliAnalysisTaskJetsEEC::JetSelectionType jetSelection = AliAnalysisTaskJetsEEC::kInclusive,
                             Float_t minpTHTrigger =0.,  Float_t maxpTHTrigger =0., Float_t acut =0.6, AliAnalysisTaskJetsEEC::DerivSubtrOrder derivSubtrOrder = AliAnalysisTaskJetsEEC::kSecondOrder ) {
 

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskJetsEEC","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskJetsEEC", "This task requires an input event handler");
      return NULL;
    }


  TString wagonName1 = Form("JetsEEC_%s_TC%s%s",njetsBase,trigClass.Data(),tag.Data());
  TString wagonName2 = Form("JetsEEC_%s_TC%s%sTree",njetsBase,trigClass.Data(),tag.Data());
  //Configure jet tagger task
  AliAnalysisTaskJetsEEC *task = new AliAnalysisTaskJetsEEC(wagonName1.Data());

  //task->SetNCentBins(4);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);
  
 

  TString thename(njetsBase);
  //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackCont;// = task->AddTrackContainer(ntracks);
 
  if ((jetShapeSub==AliAnalysisTaskJetsEEC::kConstSub || jetShapeSub==AliAnalysisTaskJetsEEC::kEventSub ) && ((jetShapeType==AliAnalysisTaskJetsEEC::kData) || (jetShapeType==AliAnalysisTaskJetsEEC::kDetEmbPartPythia) || (jetShapeType==AliAnalysisTaskJetsEEC::kPythiaDef))){
    trackCont = task->AddParticleContainer(ntracks);}
  else trackCont = task->AddTrackContainer(ntracks);

  
  //Printf("tracks() = %s, trackCont =%p", ntracks, trackCont);
  AliParticleContainer *trackContUS  = task->AddTrackContainer(ntracksUS);
  //Printf("tracksUS() = %s", ntracksUS);
  AliParticleContainer *trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  //Printf("ntracksTrue() = %s, trackContTrue=%p ", ntracksTrue, trackContTrue);
   if (jetShapeType==AliAnalysisTaskJetsEEC::kDetEmbPartPythia) trackContTrue->SetIsEmbedding(true);
  AliParticleContainer *trackContPartLevel=0;
  
  if ((jetShapeSub==AliAnalysisTaskJetsEEC::kConstSub) && ((jetShapeType==AliAnalysisTaskJetsEEC::kMCTrue) || (jetShapeType==AliAnalysisTaskJetsEEC::kPythiaDef))){
    trackContPartLevel = task->AddParticleContainer(ntracksPartLevel);
  }
  else trackContPartLevel = task->AddMCParticleContainer(ntracksPartLevel);
    if (jetShapeType==AliAnalysisTaskJetsEEC::kDetEmbPartPythia) trackContPartLevel->SetIsEmbedding(true);
   //Printf("ntracksPartLevel() = %s, trackContPartLevel=%p ", ntracksPartLevel, trackContPartLevel);
  

  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
  
  AliJetContainer *jetContBase=0x0;
  AliJetContainer *jetContUS=0x0;
  AliJetContainer *jetContTrue=0x0;
  AliJetContainer *jetContPart=0x0;
  TString strType(type);

  if ((jetShapeType==AliAnalysisTaskJetsEEC::kMCTrue || (jetShapeType==AliAnalysisTaskJetsEEC::kGenOnTheFly))) {
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackContPartLevel);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
    }
  }
  
  if (jetShapeType==AliAnalysisTaskJetsEEC::kData){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
      if(jetShapeSub==AliAnalysisTaskJetsEEC::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }
  }
  

  if (jetShapeType==AliAnalysisTaskJetsEEC::kDetEmbPartPythia){
    jetContBase = task->AddJetContainer(njetsBase,strType,R);
    if(jetContBase) {
      jetContBase->SetRhoName(nrhoBase);
      jetContBase->ConnectParticleContainer(trackCont);
      jetContBase->ConnectClusterContainer(clusterCont);
      jetContBase->SetPercAreaCut(acut);
     
      if(jetShapeSub==AliAnalysisTaskJetsEEC::kConstSub) jetContBase->SetAreaEmcCut(-2);
    }

    jetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(jetContTrue) {
      
      jetContTrue->SetRhoName(nrhoBase);
      jetContTrue->ConnectParticleContainer(trackContTrue);
      jetContTrue->SetPercAreaCut(acut);
    
    }
    
    if(jetShapeSub==AliAnalysisTaskJetsEEC::kConstSub || jetShapeSub==AliAnalysisTaskJetsEEC::kEventSub){
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
  
  if (jetShapeType==AliAnalysisTaskJetsEEC::kPythiaDef){
    
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
    
    if(jetShapeSub==AliAnalysisTaskJetsEEC::kConstSub){
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

  if (jetShapeType == AliAnalysisTaskJetsEEC::kMCTrue) contName1 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskJetsEEC::kData) contName1 += "_Data";
  if (jetShapeType == AliAnalysisTaskJetsEEC::kGenOnTheFly) contName1 += "_GenOnTheFly";
  if (jetShapeType == AliAnalysisTaskJetsEEC::kPythiaDef) contName1 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskJetsEEC::kNoSub) contName1 += "_NoSub";
  if (jetShapeSub == AliAnalysisTaskJetsEEC::kConstSub) contName1 += "_ConstSub";
   if (jetShapeSub == AliAnalysisTaskJetsEEC::kEventSub) contName1 += "_EventSub";
  if (jetShapeSub == AliAnalysisTaskJetsEEC::kDerivSub) contName1 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskJetsEEC::kInclusive) contName1 += "_Incl";
 


    if (jetShapeType == AliAnalysisTaskJetsEEC::kMCTrue) contName2 += "_MCTrue";
  if (jetShapeType == AliAnalysisTaskJetsEEC::kData) contName2 += "_Data";
  if (jetShapeType == AliAnalysisTaskJetsEEC::kGenOnTheFly) contName2 += "_GenOnTheFly";
  if (jetShapeType == AliAnalysisTaskJetsEEC::kPythiaDef) contName2 +="_PythiaDef";
  if (jetShapeSub == AliAnalysisTaskJetsEEC::kNoSub) contName2 += "_NoSub";
  if (jetShapeSub == AliAnalysisTaskJetsEEC::kConstSub) contName2 += "_ConstSub";
  if (jetShapeSub == AliAnalysisTaskJetsEEC::kEventSub) contName2 += "_EventSub";
  if (jetShapeSub == AliAnalysisTaskJetsEEC::kDerivSub) contName2 += "_DerivSub";
  
  if (jetSelection == AliAnalysisTaskJetsEEC::kInclusive) contName2 += "_Incl";




  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
   mgr->ConnectOutput(task,1,coutput1);


   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);

  return task;

}

