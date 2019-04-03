/*  created by fbellini@cern.ch on 29/04/2013 */
/*  last modified by fbellini   on 29/04/2013 */

// UInt_t kTriggerMuonAll = AliVEvent::kMUL7 | AliVEvent::kMUSH7 | AliVEvent::kMUU7 | AliVEvent::kMUS7;
// UInt_t kTriggerMuonBarell = AliVEvent::kMUU7;
// UInt_t kTriggerEMC   = AliVEvent::kEMC7;
// UInt_t kTriggerHM   = AliVEvent::kHighMult;
// UInt_t kTriggerInt = AliVEvent::kAnyINT;
// UInt_t kTriggerMask = kTriggerInt;

AliAnalysisTaskTOFqaID * AddTaskTOFqaID(Bool_t  flagEnableAdvancedCheck = kFALSE, 
				   UInt_t  triggerMask = AliVEvent::kAnyINT, 
				   Int_t   trackCutSetTOFqa = AliAnalysisTaskTOFqaID::kStd2011crossedRows, 
				   Bool_t  flagEnableChargeSplit = kFALSE,
				   TString cutName = "",
				   Bool_t  isMC = kFALSE, 
				   Short_t absPdgCode = 0,
				   UInt_t  runN = 0,
				   const char *cdb     = "raw://")
{
  // Task for checking TOF QA
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
   ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }   
  runN = mgr->GetRunFromPath();
  Printf(":::::::  HERE I TRY TO GET RUN  NUMBER = %i", runN);

  // Check the analysis type using the event handlers connected to the analysis manager.
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 
  if (isMC && (absPdgCode<0)) {
    ::Error("AddTask","Invalid selected PDG code for MC. Set absPdgCode=0 for no specie selection.");
  }
 
  // Create the task
  AliAnalysisTaskTOFqaID *task = new AliAnalysisTaskTOFqaID(Form("taskTOFqaID_%i",absPdgCode));
  task->SetOCDBInfo(cdb, runN);
  task->EnableAdvancedCheck(flagEnableAdvancedCheck);
  task->EnableChargeSplit(flagEnableChargeSplit);
  task->SetSelectMCspecies(isMC, absPdgCode);
  task->SelectCollisionCandidates(triggerMask);
  //AliLog::SetClassDebugLevel("AliAnalysisTaskTOFqaID",4);
  mgr->AddTask(task);

  //Define track cut set
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "TrackCutsTOFqa"); 
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  
  if ( (trackCutSetTOFqa<0) || (trackCutSetTOFqa>=AliAnalysisTaskTOFqaID::kNCutSetTOFqa) ) trackCutSetTOFqa = 0;

  Bool_t selPrimaries=kTRUE;
  Int_t clusterCut = 1; //set to 1 to use crossed rows and to 0 to use clusters
  Bool_t cutAcceptanceEdges = kTRUE;
  Bool_t removeDistortedRegions = kFALSE;




  Printf(":::: Setting cut scheme #%i", trackCutSetTOFqa);
  if (trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kRun1Cuts) {
    //use track cuts used for QA during run1 (before July 2014)
    esdTrackCuts->SetMinNClustersTPC(70); 
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE); 
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					   AliESDtrackCuts::kAny);
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//selects primaries
    esdTrackCuts->SetMaxDCAToVertexZ(2);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  }
  if (trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kStd2010)
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selPrimaries, 0);
  else if (trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kStd2010crossedRows)
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selPrimaries, 1);
  else if (trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kStd2011 )
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(selPrimaries, 0);
  else if (trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kStd2011crossedRows )
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(selPrimaries, 1);
  else if (trackCutSetTOFqa == AliAnalysisTaskTOFqaID::kStd2015crossedRows )
    esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb(selPrimaries, clusterCut, cutAcceptanceEdges, removeDistortedRegions);
  else {
    //use track cuts used for QA during run1 (before July 2014)
    esdTrackCuts->SetMinNClustersTPC(70); 
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE); 
    esdTrackCuts->SetRequireTPCRefit(kTRUE);
    esdTrackCuts->SetRequireITSRefit(kTRUE);
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					   AliESDtrackCuts::kAny);
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//selects primaries
    esdTrackCuts->SetMaxDCAToVertexZ(2);
    esdTrackCuts->SetDCAToVertex2D(kFALSE);
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  }
  esdTrackCuts->Dump();
  trackFilter->AddCuts(esdTrackCuts);
  task->SetTrackFilter(trackFilter); 
  
  TString partName(task->GetSpeciesName(absPdgCode));

  // Create containers for input/output
  AliAnalysisDataContainer* cInputTOFqa = mgr->CreateContainer("cInputTOFqa", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer* cGeneralTOFqa = mgr->CreateContainer(Form("base_%s%s", partName.Data(), cutName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TOF", mgr->GetCommonFileName()));
  AliAnalysisDataContainer* cTimeZeroTOFqa = mgr->CreateContainer(Form("timeZero_%s%s", partName.Data(), cutName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TOF", mgr->GetCommonFileName()));
  AliAnalysisDataContainer* cPIDTOFqa = mgr->CreateContainer(Form("pid_%s%s", partName.Data(), cutName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TOF", mgr->GetCommonFileName()));
  AliAnalysisDataContainer* cTRDcheckTOFqa = mgr->CreateContainer(Form("trd_%s%s", partName.Data(), cutName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TOF", mgr->GetCommonFileName()));
  AliAnalysisDataContainer* cTriggerTOFqa = mgr->CreateContainer(Form("trigger_%s%s", partName.Data(), cutName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:TOF", mgr->GetCommonFileName()));

  // Attach i/o
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cGeneralTOFqa);
  mgr->ConnectOutput(task, 2, cTimeZeroTOFqa);
  mgr->ConnectOutput(task, 3, cPIDTOFqa);
  mgr->ConnectOutput(task, 4, cTRDcheckTOFqa);
  mgr->ConnectOutput(task, 5, cTriggerTOFqa);
  return task;
}
