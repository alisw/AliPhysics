AliAnalysisTaskRecoilJetYield* AddTaskRecoilJetYield(const char * njetsData, //data jets
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
								    const char * nclusters,
								    const char *type,				      
								    const char *CentEst,
								    Double_t fSharedFractionPtMin,
								    Int_t SubJetAlgorithm,
								    Float_t SubJetRadius,
								    Float_t SubJetMinPt,
								    Int_t       pSel,
								    TString     trigClass      = "",
								    TString     kEmcalTriggers = "",
								    TString     tag            = "",
								    AliAnalysisTaskRecoilJetYield::JetShapeType jetShapeType,
								    AliAnalysisTaskRecoilJetYield::JetShapeSub jetShapeSub,
								    AliAnalysisTaskRecoilJetYield::JetSelectionType jetSelection,
								    Float_t minpTHTrigger =0.,  Float_t maxpTHTrigger =0., AliAnalysisTaskRecoilJetYield::DerivSubtrOrder derivSubtrOrder = 0  ) {
  
  
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskAliAnalysisTaskRecoilJetYield","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AliAnalysisTaskRecoilJetYield", "This task requires an input event handler");
      return NULL;
    }

  if (jetShapeType==AliAnalysisTaskRecoilJetYield::kData || jetShapeType==AliAnalysisTaskRecoilJetYield::kSim){
    TString wagonName1 = Form("AliAnalysisTaskRecoilJetYield_%s_TC%s%s",njetsData,trigClass.Data(),tag.Data());
    TString wagonName2 = Form("AliAnalysisTaskRecoilJetYield_%s_TC%s%sTree",njetsData,trigClass.Data(),tag.Data());
  }
  if (jetShapeType==AliAnalysisTaskRecoilJetYield::kTrue || jetShapeType==AliAnalysisTaskRecoilJetYield::kTrueDet || jetShapeType==AliAnalysisTaskRecoilJetYield::kGenOnTheFly){
    TString wagonName1 = Form("AliAnalysisTaskRecoilJetYield_%s_TC%s%s",njetsTrue,trigClass.Data(),tag.Data());
    TString wagonName2 = Form("AliAnalysisTaskRecoilJetYield_%s_TC%s%sTree",njetsTrue,trigClass.Data(),tag.Data());
  }
  if (jetShapeType==AliAnalysisTaskRecoilJetYield::kDetEmbPart){
    TString wagonName1 = Form("AliAnalysisTaskRecoilJetYield_%s_TC%s%s",njetsHybridS,trigClass.Data(),tag.Data());
    TString wagonName2 = Form("AliAnalysisTaskRecoilJetYield_%s_TC%s%sTree",njetsHybridS,trigClass.Data(),tag.Data());
  }
  //Configure jet tagger task
  AliAnalysisTaskRecoilJetYield *task = new AliAnalysisTaskRecoilJetYield(wagonName1.Data());

  //  task->SetNCentBins(4);
  task->SetJetShapeType(jetShapeType);
  task->SetJetShapeSub(jetShapeSub);
  task->SetJetSelection(jetSelection);
  task->SetSubJetAlgorithm(SubJetAlgorithm);
  task->SetSubJetRadius(SubJetRadius);
  task->SetSubJetMinPt(SubJetMinPt);
  task->SetJetRadius(R);
  task->SetSharedFractionPtMin(fSharedFractionPtMin);
  task->SetDerivativeSubtractionOrder(derivSubtrOrder);
  if (jetSelection == AliAnalysisTaskRecoilJetYield::kRecoil) task->SetPtTriggerSelections(minpTHTrigger, maxpTHTrigger);

  // TString thename(njetsBase);
  //if(thename.Contains("Sub")) task->SetIsConstSub(kTRUE);
  //task->SetVzRange(-10.,10.);

  AliParticleContainer *trackContData=0x0; 
  AliParticleContainer *trackContDet=0x0;
  AliParticleContainer *trackContTrue=0x0;

  if (jetShapeSub == AliAnalysisTaskRecoilJetYield::kConstSub){
    trackContData = task->AddParticleContainer(ntracksData);
    trackContDet = task->AddParticleContainer(ntracksDet);
  }
  else{
    trackContData = task->AddTrackContainer(ntracksData);
    trackContDet = task->AddMCParticleContainer(ntracksDet);
  }
  if ((jetShapeType==AliAnalysisTaskRecoilJetYield::kData || jetShapeType==AliAnalysisTaskRecoilJetYield::kSim) && (jetShapeSub == AliAnalysisTaskRecoilJetYield::kConstSub)){
    trackContTrue = task->AddTrackContainer(ntracksTrue); //Unsubtracted tracks branch
  }
  else trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  AliParticleContainer *trackContHybridUs = task->AddParticleContainer(ntracksHybridUs);
  AliParticleContainer *trackContHybridS = task->AddParticleContainer(ntracksHybridS);
  AliClusterContainer *clusterCont = task->AddClusterContainer(nclusters);
    
  AliJetContainer *JetContData=0x0;
  AliJetContainer *JetContTrue=0x0;
  AliJetContainer *JetContDet=0x0;
  AliJetContainer *JetContHybridUs=0x0;
  AliJetContainer *JetContHybridS=0x0;

  TString strType(type);
  
  ///////kTrue///////
  if (jetShapeType==AliAnalysisTaskRecoilJetYield::kTrue) {
    JetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(JetContTrue) {
      JetContTrue->SetRhoName(nrhoBase);
      JetContTrue->ConnectParticleContainer(trackContTrue);
      JetContTrue->ConnectClusterContainer(clusterCont);
      JetContTrue->SetPercAreaCut(0.6);
      JetContTrue->SetJetRadius(R);
      JetContTrue->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
    }
  }
  ///////kGenOnTheFly///////
  if (jetShapeType==AliAnalysisTaskRecoilJetYield::kGenOnTheFly) {
    JetContTrue = task->AddJetContainer(njetsTrue,strType,R);
    if(JetContTrue) {
      JetContTrue->SetRhoName(nrhoBase);
      JetContTrue->ConnectParticleContainer(trackContTrue);
      JetContTrue->ConnectClusterContainer(clusterCont);
      JetContTrue->SetJetRadius(R);
      JetContTrue->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      JetContTrue->SetPercAreaCut(0.6);
    }
    if(jetShapeSub==AliAnalysisTaskRecoilJetYield::kConstSub){                                                                             
      JetContDet=task->AddJetContainer(njetsDet,strType,R);     //So we can access the unsubtracted particle container                                                                            
      if(JetContDet) {
        JetContDet->SetRhoName(nrhoBase);
        JetContDet->ConnectParticleContainer(trackContDet);
        JetContDet->SetPercAreaCut(0.6);
	JetContDet->SetJetRadius(R);
	JetContDet->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      }
    }  
  }

  ///////kTrueDet///////
  if (jetShapeType==AliAnalysisTaskRecoilJetYield::kTrueDet){
    
    JetContDet = task->AddJetContainer(njetsDet,strType,R);  //detector level MC
    if(JetContDet) {
      JetContDet->SetRhoName(nrhoBase);
      JetContDet->ConnectParticleContainer(trackContDet);
      JetContDet->ConnectClusterContainer(clusterCont);
      JetContDet->SetPercAreaCut(0.6);
      JetContDet->SetJetRadius(R);
      JetContDet->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
    }

    JetContTrue = task->AddJetContainer(njetsTrue,strType,R); //Particle Level MC
    if(JetContTrue) {
      JetContTrue->SetRhoName(nrhoBase);
      JetContTrue->ConnectParticleContainer(trackContTrue);
      JetContTrue->SetPercAreaCut(0.6);
      JetContTrue->SetJetRadius(R);
      JetContTrue->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
    }
  }  

  ///////kData||kSim///////
  if (jetShapeType==AliAnalysisTaskRecoilJetYield::kData || jetShapeType==AliAnalysisTaskRecoilJetYield::kSim){
    JetContData = task->AddJetContainer(njetsData,strType,R); //Data
    if(JetContData) {
      JetContData->SetRhoName(nrhoBase);
      JetContData->ConnectParticleContainer(trackContData);
      JetContData->ConnectClusterContainer(clusterCont);
      JetContData->SetPercAreaCut(0.6);
      JetContData->SetJetRadius(R);
      JetContData->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecoilJetYield::kConstSub) JetContData->SetAreaEmcCut(-2);
    }
    if(jetShapeSub==AliAnalysisTaskRecoilJetYield::kConstSub){                                                                             
      JetContTrue=task->AddJetContainer(njetsTrue,strType,R);     //So we can access the unsubtracted particle container                                                                            
      if(JetContTrue) {
        JetContTrue->SetRhoName(nrhoBase);
        JetContTrue->ConnectParticleContainer(trackContTrue);
        JetContTrue->SetPercAreaCut(0.6);
	JetContTrue->SetJetRadius(R);
	JetContTrue->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      }
    }   
  }

  ///////kDetEmbPart///////
  if (jetShapeType==AliAnalysisTaskRecoilJetYield::kDetEmbPart){
    JetContHybridS = task->AddJetContainer(njetsHybridS,strType,R);  //Subtracted Hybrid (Pb+Pyhthia Det Level)                                                               
    if(JetContHybridS) {
      JetContHybridS->SetRhoName(nrhoBase);
      JetContHybridS->ConnectParticleContainer(trackContHybridS);
      JetContHybridS->ConnectClusterContainer(clusterCont);
      JetContHybridS->SetPercAreaCut(0.6);
      JetContHybridS->SetJetRadius(R);
      JetContHybridS->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
      if(jetShapeSub==AliAnalysisTaskRecoilJetYield::kConstSub) JetContHybridS->SetAreaEmcCut(-2); //??????????                                                             
    }
    //  if(jetShapeSub==AliAnalysisTaskRecoilJetYield::kConstSub){ // In the case of Non Const Subtracted fill the Unsubtracted and Subtracted jet branches with the same hybrid branch. But you need to fill Both! For simplicity of coding                                                                              
    JetContHybridUs=task->AddJetContainer(njetsHybridUs,strType,R);     //Unsubtracted Hybrid                                                                              
    if(JetContHybridUs) {
      JetContHybridUs->SetRhoName(nrhoBase);
      JetContHybridUs->ConnectParticleContainer(trackContHybridUs);
      JetContHybridUs->SetPercAreaCut(0.6);
      JetContHybridUs->SetJetRadius(R);
      JetContHybridUs->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
    }
    //    }
    JetContDet = task->AddJetContainer(njetsDet,strType,R); //Pythia Detector Level                                                                                        
    if(JetContDet) {
      JetContDet->SetRhoName(nrhoBase);
      JetContDet->ConnectParticleContainer(trackContDet);
      JetContDet->SetPercAreaCut(0.6);
      JetContDet->SetJetRadius(R);
      JetContDet->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
    }
    JetContTrue = task->AddJetContainer(njetsTrue,strType,R); //Pyhthia Particle Level                                                                                
    if(JetContTrue) {
      JetContTrue->SetRhoName(nrhoBase);
      JetContTrue->ConnectParticleContainer(trackContTrue);
      JetContTrue->SetPercAreaCut(0.6);
      JetContTrue->SetJetRadius(R);
      JetContTrue->SetJetAcceptanceType(AliEmcalJet::kTPCfid);
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
  
  if (jetShapeType == AliAnalysisTaskRecoilJetYield::kTrue){
    contName1 += "_True";
    contName2 += "_True";
  }
  if (jetShapeType == AliAnalysisTaskRecoilJetYield::kTrueDet){
    contName1 += "_TrueDet";
    contName2 += "_TrueDet";
  }
  if (jetShapeType == AliAnalysisTaskRecoilJetYield::kData){
    contName1 += "_Data";
    contName2 += "_Data";
  }
  if (jetShapeType == AliAnalysisTaskRecoilJetYield::kSim){
    contName1 += "_Sim";
    contName2 += "_Sim";
  }
  if (jetShapeType == AliAnalysisTaskRecoilJetYield::kDetEmbPart){
    contName1 += "_DetEmbPart";
    contName2 += "_DetEmbPart";
  }
  if (jetShapeType == AliAnalysisTaskRecoilJetYield::kGenOnTheFly){
    contName1 += "_GenOnTheFly";
    contName2 += "_GenOnTheFly";
  }

  if (jetShapeSub == AliAnalysisTaskRecoilJetYield::kNoSub){
    contName1 += "_NoSub";
    contName2 += "_NoSub";
  }
  if (jetShapeSub == AliAnalysisTaskRecoilJetYield::kConstSub){
    contName1 += "_ConstSub";
    contName2 += "_ConstSub";
  }
  if (jetShapeSub == AliAnalysisTaskRecoilJetYield::kDerivSub && derivSubtrOrder == 0){
    contName1 += "_DerivSubSecondOrder";
    contName2 += "_DerivSubSecondOrder";
  }
  if (jetShapeSub == AliAnalysisTaskRecoilJetYield::kDerivSub && derivSubtrOrder == 1){
    contName1 += "_DerivSubFirstOrder";
    contName2 += "_DerivSubFirstOrder";
  }
  if (jetSelection == AliAnalysisTaskRecoilJetYield::kInclusive){
    contName1 += "_Incl";
    contName2 += "_Incl";
  }
  if (jetSelection == AliAnalysisTaskRecoilJetYield::kRecoil) {
  TString recoilTriggerString = Form("_Recoil_%.0f_%0.f", minpTHTrigger, maxpTHTrigger);
  contName1 += recoilTriggerString;
  contName2 += recoilTriggerString;
  }
  TString SubJetRadiusString = Form("_SubJetRadius_%f", SubJetRadius);
  contName1 += SubJetRadiusString;
  contName2 += SubJetRadiusString;
  TString SubJetMinPtString = Form("_SubJetMinPt_%f", SubJetMinPt);
  contName1 += SubJetMinPtString;
  contName2 += SubJetMinPtString;
 



  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);

  return task;  

}

