AliAnalysisTask *AddTaskTrackingUncertAOT(Bool_t readMC = kFALSE,
                                          TString trigClass = "CINT1B",
                                          AliVEvent::EOfflineTriggerTypes trigMask = AliVEvent::kMB,
                                          AliAnalysisTrackingUncertaintiesAOT::ESpecies_t specie=(AliAnalysisTrackingUncertaintiesAOT::kSpecPion|AliAnalysisTrackingUncertaintiesAOT::kSpecKaon),
                                          TString outputSuffix = "",
                                          Bool_t doCutV0multTPCout = kFALSE,
                                          AliAnalysisTrackingUncertaintiesAOT::ECentrality centrSel = AliAnalysisTrackingUncertaintiesAOT::kCentOff,
                                          Double_t minCentrality = 0.,
                                          Double_t maxCentrality = 100.,
                                          Double_t MaxDCAxy = 2.4,
                                          Double_t MaxDCAz  = 3.2,
                                          Double_t MaxEta   = 0.8,
                                          Double_t CrossRowsOverFndCltTPC = 0.8,
                                          AliESDtrackCuts::ITSClusterRequirement spdReq=AliESDtrackCuts::kAny,
                                          Bool_t useGenPt = kFALSE) {
    
  //
  // add task of tracking uncertainty
  //
  //
  //get the current analysis manager
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTrackingUncertAOT", "No analysis manager found.");
    return 0;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskImpParDistrib", "This task requires an input event handler");
    return NULL;
  }
    
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AddTaskImpParDistrib", "This task requires to run on ESD");
    return NULL;
  }
  //
  //========= Add task for standard analysis to the ANALYSIS manager ====
  //
  AliAnalysisTrackingUncertaintiesAOT *task    = new AliAnalysisTrackingUncertaintiesAOT("trackingUncertainty");
  //
  task->SetReadMC(readMC);
  task->SetTriggerClass(trigClass.Data());
  task->SetTriggerMask(trigMask);
  task->SetSpecie(specie);
  task->SetMaxDCAxy(MaxDCAxy);
  task->SetMaxDCAz(MaxDCAz);
  task->SetEtaRange(MaxEta);
  task->SetCrossRowsOverFndCltTPC(CrossRowsOverFndCltTPC);
  task->SetUseCutV0multVsTPCout(doCutV0multTPCout);
  task->SetSPDRequirement(spdReq);
  task->SetUseCentrality(centrSel);
  task->SetMinCentrality(minCentrality);
  task->SetMaxCentrality(maxCentrality);
  task->SetUseGeneratedPt(useGenPt);
    
  mgr->AddTask(task);
  ULong64_t SPeciee = task->GetSpecie();
  TString suffix = "_";
  if(SPeciee&AliAnalysisTrackingUncertaintiesAOT::kSpecElectron) suffix += "e";
  if(SPeciee&AliAnalysisTrackingUncertaintiesAOT::kSpecPion)     suffix += "Pi";
  if(SPeciee&AliAnalysisTrackingUncertaintiesAOT::kSpecKaon)     suffix += "K";
  if(SPeciee&AliAnalysisTrackingUncertaintiesAOT::kSpecProton)   suffix += "p";
  if(SPeciee&AliAnalysisTrackingUncertaintiesAOT::kAll)          suffix  = "_All";
  suffix += outputSuffix;
  //
  //
  //======================================================================
  //              data containers
  //======================================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
    
  //define output containers
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(Form("trackingUncert%s",suffix.Data()), TList::Class(),AliAnalysisManager::kOutputContainer,"AnalysisResults.root");
    
  //connect containers
  mgr->ConnectInput  (task, 0, cinput );
  mgr->ConnectOutput (task, 1, coutput1);
  //
  //
  //
  return task;
    
}

