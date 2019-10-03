AliEPSelectionTask *AddTaskEventplane(Bool_t useEtaGap=kTRUE,Float_t etaGap=0.4,Bool_t posTPCAOD=kFALSE,TString containername = "EPStat")
{
  // Macro to connect an event plane selection task to an existing analysis manager.

  if(useEtaGap && posTPCAOD){
    ::Error("AddTaskEventplane", "eta-splitting of events and one side of TPC not possible at same time!");
    return NULL;
  }
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEventplane", "No analysis manager to connect to.");
    return NULL;
  }      
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEventplane", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 
  AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection");
  eventplaneTask->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral | AliVEvent::kCentral);
  if (inputDataType == "AOD"){
    eventplaneTask->SetInput("AOD");
  }
  eventplaneTask->SetTrackType("TPC");
  eventplaneTask->SetUsePtWeight();
  eventplaneTask->SetUsePhiWeight();
  eventplaneTask->SetSaveTrackContribution();
  if(useEtaGap){
    eventplaneTask->SetSubeventsSplitMethod(AliEPSelectionTask::kEta); 
    eventplaneTask->SetEtaGap(etaGap); 
  }
  if(posTPCAOD){
    eventplaneTask->SetPersonalAODtrackCuts(128,0.,0.8,0.15,20.);
    eventplaneTask->SetSubeventsSplitMethod(AliEPSelectionTask::kRandom);
  }
   // Cuts on primary tracks
  AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Standard");
  esdTrackCutsL->SetMinNClustersTPC(50);
  esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
//  esdTrackCutsL->SetMaxCovDiagonalElements(2, 2, 0.5, 0.5, 2);
//  esdTrackCutsL->SetRequireTPCRefit(kTRUE);
  esdTrackCutsL->SetMaxDCAToVertexXY(3.0);
  esdTrackCutsL->SetMaxDCAToVertexZ(3.0);
  esdTrackCutsL->SetDCAToVertex2D(kTRUE);
//  esdTrackCutsL->SetRequireSigmaToVertex(kFALSE);
//  esdTrackCutsL->SetAcceptKinkDaughters(kFALSE);
  eventplaneTask->SetPersonalESDtrackCuts(esdTrackCutsL);


  
  mgr->AddTask(eventplaneTask);

  AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containername,
                TList::Class(), AliAnalysisManager::kOutputContainer,
                "EventStat_temp.root");
  
  mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(eventplaneTask,1,coutput1);

  return eventplaneTask;
}
