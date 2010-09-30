AliAnalysisTask *AddTaskHFE(){
  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFE", "No analysis manager found.");
    return NULL;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskHFE", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if (type=="AOD"){
    ::Error("AddTaskHFE", "The tasks exits because AODs are in input");
    return NULL;
  }
  Bool_t MCthere=kTRUE;
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
  if(!mcH){
    MCthere=kFALSE;
  }

  //============= Set Task Name ===================
  //TString taskName=("AliAnalysisTaskHFE.cxx+");
  //===============================================

  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCuts","HFE Standard Cuts");
  hfecuts->CreateStandardCuts();
  hfecuts->SetMinNClustersTPC(110);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);
  hfecuts->SetSigmaToVertex(10);
  hfecuts->SetQAOn();
  //hfecuts->SetMinNTrackletsTRD(5);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE("HFEanalysis");
  task->SetESDAnalysis();
  if (MCthere)
    task->SetHasMCData(kTRUE);
  else{
    task->SetHasMCData(kFALSE);
    task->SelectCollisionCandidates();
  }

  task->SetPIDStrategy(6);
  task->SetHFECuts(hfecuts);
  if(!MCthere){
    TF1 *hBackground = new TF1("hadronicBackgroundFunction", "[0]+[1]*TMath::Erf([2]*x+[3])", 0, 20);
    hBackground->SetParameter(0, 0.2194);
    hBackground->SetParameter(1, 0.217);
    hBackground->SetParameter(2, 0.6829);
    hBackground->SetParameter(3, -2.697);
    task->SetBackGroundFactorsFunction(hBackground);
  }

  // kPIDqa needs to be off for flat pT spectra !!!
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);

  mgr->AddTask(task);

  //----------------------
  //create data containers
  //----------------------
 
  //find input container
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  TString containerName = mgr->GetCommonFileName();
  containerName += ":PWG3_hfe";
  
  task->ConnectOutput(1, mgr->CreateContainer("HFE_nEvents", TH1I::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  task->ConnectOutput(2, mgr->CreateContainer("HFE_Results", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  task->ConnectOutput(3, mgr->CreateContainer("HFE_QA", TList::Class(),
					      AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectInput  (task,  0, cinput );
  
  return task;
}
