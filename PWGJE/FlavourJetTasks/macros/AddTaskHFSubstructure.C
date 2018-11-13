AliAnalysisTaskHFSubstructure* AddTaskAliAnalysisTaskHFSubstructure(const char * ntracksData,
                                                                    const char * ntracksDet,
                                                                    const char * ntracksTrue,
								    const Double_t R,
								    AliAnalysisTaskHFSubstructure::ECandidateType_t ECandidateType = AliAnalysisTaskHFSubstructure::kD0toKpi,
								    AliAnalysisTaskHFSubstructure::JetShapeType jetShapeType = AliAnalysisTaskHFSubstructure::kData) {
  
  
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
    {
      Error("AddTaskAliAnalysisTaskHFSubstructure","No analysis manager found.");
      return 0;
    }
  Bool_t ismc=kFALSE;
  ismc = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
    {
      ::Error("AliAnalysisTaskHFSubstructure", "This task requires an input event handler");
      return NULL;
    }
  TString wagonName1, wagonName2,wagonName3;
  TString tag="";
  if (ECandidateType == AliAnalysisTaskHFSubstructure::kD0toKpi) tag="kD0toKpi";
  if (jetShapeType==AliAnalysisTaskHFSubstructure::kData){
    wagonName1 = Form("AliAnalysisTaskHFSubstructure_%s_TC%s",ntracksData,tag.Data());
    wagonName2 = Form("AliAnalysisTaskHFSubstructure_%s_TC%sTree",ntracksData,tag.Data());
    wagonName3 = Form("AliAnalysisTaskHFSubstructure_%s_TC%sTreeSplittings",ntracksData,tag.Data());
  }
  if (jetShapeType == AliAnalysisTaskHFSubstructure::kDetSignal || jetShapeType == AliAnalysisTaskHFSubstructure::kDetBackground || jetShapeType == AliAnalysisTaskHFSubstructure::kDetReflection || jetShapeType == AliAnalysisTaskHFSubstructure::kTrueDet){
    wagonName1 = Form("AliAnalysisTaskHFSubstructure_%s_TC%s",ntracksDet,tag.Data());
    wagonName2 = Form("AliAnalysisTaskHFSubstructure_%s_TC%sTree",ntracksDet,tag.Data());
    wagonName3 = Form("AliAnalysisTaskHFSubstructure_%s_TC%sTreeSplittings",ntracksDet,tag.Data());
  }
  if (jetShapeType==AliAnalysisTaskHFSubstructure::kTrue){
    wagonName1 = Form("AliAnalysisTaskHFSubstructure_%s_TC%s",ntracksTrue,tag.Data());
    wagonName2 = Form("AliAnalysisTaskHFSubstructure_%s_TC%sTree",ntracksTrue,tag.Data());
    wagonName3 = Form("AliAnalysisTaskHFSubstructure_%s_TC%sTreeSplittings",ntracksTrue,tag.Data());
  }
  //Configure jet tagger task
  AliAnalysisTaskHFSubstructure *task = new AliAnalysisTaskHFSubstructure(wagonName1);


  task->SetECandidateType_t(ECandidateType);
  task->SetJetShapeType(jetShapeType);
  task->SetJetRadius(R);


  AliParticleContainer *trackContData=0x0;  //why not track containers?
  AliParticleContainer *trackContDet=0x0;
  AliParticleContainer *trackContTrue=0x0;


  if (jetShapeType == AliAnalysisTaskHFSubstructure::kData) trackContData = task->AddTrackContainer(ntracksData);
  else if (jetShapeType == AliAnalysisTaskHFSubstructure::kDetSignal || jetShapeType == AliAnalysisTaskHFSubstructure::kDetBackground || jetShapeType == AliAnalysisTaskHFSubstructure::kDetReflection || jetShapeType == AliAnalysisTaskHFSubstructure::kTrueDet){
    trackContDet = task->AddTrackContainer(ntracksDet);
    trackContTrue = task->AddMCParticleContainer(ntracksTrue);
  }
  else if (jetShapeType == AliAnalysisTaskHFSubstructure::kTrue) trackContTrue = task->AddMCParticleContainer(ntracksTrue);

  


  task->SetUseAliAnaUtils(kFALSE);

  mgr->AddTask(task);
  
  //Connnect input
  mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer() );

  //Connect output
  TString contName1(wagonName1);
  TString contName2(wagonName2);
  TString contName3(wagonName3);
  
  contName2 += "_Splittings";
  

  if (jetShapeType == AliAnalysisTaskHFSubstructure::kData){
    contName1 += "_Data";
    contName2 += "_Data";
  }

  if (jetShapeType == AliAnalysisTaskHFSubstructure::kDetSignal){
    contName1 += "_DetSignal";
    contName2 += "_DetSignal";
  }
  
  if (jetShapeType == AliAnalysisTaskHFSubstructure::kDetBackground){
    contName1 += "_DetBackground";
    contName2 += "_DetBackgroundl";
  }
  
  if (jetShapeType == AliAnalysisTaskHFSubstructure::kDetReflection){
    contName1 += "_DetReflection";
    contName2 += "_DetReflection";
  }
  
  if (jetShapeType == AliAnalysisTaskHFSubstructure::kTrueDet){
    contName1 += "_TrueDet";
    contName2 += "_TrueDet";
  }
  
  if (jetShapeType == AliAnalysisTaskHFSubstructure::kTrue){
    contName1 += "_True";
    contName2 += "_True";
  }


  TString outputfile = Form("%s",AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,1,coutput1);
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,2,coutput2);
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(contName3.Data(), TTree::Class(),AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectOutput(task,3,coutput3);

  return task;  

}

