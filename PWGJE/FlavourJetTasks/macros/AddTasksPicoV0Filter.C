const Bool_t bAnaInfoMC     = kTRUE;
const Bool_t bIsInfoAOD     = kTRUE;
const Bool_t bIsPhysSel     = kFALSE;
const Bool_t bIsCentSel     = kFALSE;
const Bool_t bIsEvPnSel     = kFALSE;
const Bool_t bIsRespPID     = kTRUE;
const UInt_t wTriggerMask   = AliVEvent::kINT7;
const UInt_t wCollisionType = AliPicoHeaderCJ::kPA;

const TString sCentEsti     = "V0A";  // "V0M"; "V0A"; "V0C"
//=============================================================================

Bool_t AddTasksPicoV0Filter()
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  if (bIsInfoAOD) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
    AliAODInputHandler *aodIH = AddAODHandler();
  } else {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    AliESDInputHandler *esdIH = AddESDHandler();
//  esdIH->SetReadFriends(kFALSE);
  }

  if (bAnaInfoMC && !bIsInfoAOD) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C");
    AliMCEventHandler *mctEH = AddMCHandler(kFALSE);
  }

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODOutputHandler.C");
  AliAODHandler *aodH = AddAODOutputHandler();
  aodH->SetOutputFileName("AliAOD.PicoV0s.root");
  aodH->SetFillAOD(kTRUE);
  aodH->SetCreateNonStandardAOD();
  mgr->SetOutputEventHandler(aodH);
//=============================================================================

  if (bIsPhysSel) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *taskPhysSel = AddTaskPhysicsSelection(bAnaInfoMC);
  }

  if (bIsCentSel) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentSel = AddTaskCentrality(kTRUE, bIsAOD);
    if (wTriggerMask) taskCentSel->SelectCollisionCandidates(wTriggerMask);
    if (bAnaInfoMC)   taskCentSel->SetMCInput();
  }

  if (bIsEvPnSel) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
    AliEPSelectionTask *taskEventPlane = AddTaskEventplane();
    if (wTriggerMask) taskEventPlane->SelectCollisionCandidates(wTriggerMask);
    if (bAnaInfoMC)   taskEventPlane->SetUseMCRP();
  }

  if (bIsRespPID) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskSE *taskRespPID = AddTaskPIDResponse(bAnaInfoMC, kTRUE,
                                                        kFALSE, 2,
                                                        kFALSE, "",
                                                        kTRUE, kTRUE,
                                                        2);
    if (wTriggerMask) taskRespPID->SelectCollisionCandidates(wTriggerMask);
  }
//=============================================================================

  AliAnalysisTaskSEPicoV0Maker *taskPicoV0Maker = new AliAnalysisTaskSEPicoV0Maker("AliAnalysisTaskSEPicoV0Maker", bAnaInfoMC);
  if (wTriggerMask) {
    taskPicoV0Maker->SetTriggerMask(wTriggerMask);
    taskPicoV0Maker->SelectCollisionCandidates(wTriggerMask);
  }

  taskPicoV0Maker->SetCollitionType(wCollisionType);
  taskPicoV0Maker->SetCentralityEstimator(sCentEsti.Data());
  taskPicoV0Maker->SetVertexContributorN(1);
//taskPicoV0Maker->SetRefitV0ESD();
//taskPicoV0Maker->SetSkipFastOnly();
//taskPicoV0Maker->SetDMPjetMC();

  mgr->AddTask(taskPicoV0Maker);
  mgr->ConnectInput(taskPicoV0Maker,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPicoV0Maker, 1, mgr->CreateContainer("listPicoV0MakerEH",
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              AliAnalysisManager::GetCommonFileName()));

  if (bAnaInfoMC) mgr->ConnectOutput(taskPicoV0Maker, 2, mgr->CreateContainer("listPicoV0MakerMC",
                                                                              TList::Class(),
                                                                              AliAnalysisManager::kOutputContainer,
                                                                              AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  AliAnalysisTaskSEPicoV0Filter *taskPicoV0Filter = new AliAnalysisTaskSEPicoV0Filter("AliAnalysisTaskSEPicoV0Filter");
  if (wTriggerMask) taskPicoV0Filter->SelectCollisionCandidates(wTriggerMask);
  taskPicoV0Filter->SetAnaInfoMC(bAnaInfoMC);

  mgr->AddTask(taskPicoV0Filter);
  mgr->ConnectInput(taskPicoV0Filter, 0, mgr->GetCommonInputContainer());
//=============================================================================

/*if (bAnaInfoMC) {
    if (!bIsInfoAOD) {
      AliAnalysisTaskExtractPerformanceV0 *taskExtractV0 = new AliAnalysisTaskExtractPerformanceV0("AliAnalysisTaskExtractPerformanceV0");
      if (wTriggerMask) taskExtractV0->SelectCollisionCandidates(wTriggerMask);

      taskExtractV0->SetIsNuclear(kTRUE);
      taskExtractV0->SetINT7Trigger(kTRUE);
      taskExtractV0->SetUseOnTheFly(kFALSE);
      taskExtractV0->SetTakeAllTracks(kFALSE);
      taskExtractV0->SetpARapidityShift();
      taskExtractV0->SetCentralityEstimator(sCentEsti);
      taskExtractV0->SetLightWeightAnalysis(kFALSE);
      taskExtractV0->SetpAVertexSelection(kTRUE);
      taskExtractV0->SetSpecialExecution(kFALSE);
      taskExtractV0->SetSaveAssociatedOnly(kTRUE);
      taskExtractV0->SetSkipTrigger(kTRUE);
      taskExtractV0->SetDoNotCallTPCdEdx(kFALSE);
//    taskExtractV0->SetDiffractiveOnly(kTRUE);

      mgr->AddTask(taskExtractV0);
      mgr->ConnectInput(taskExtractV0,  0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(taskExtractV0, 1, mgr->CreateContainer("clistV0MC",
                                                                TList::Class(),
                                                                AliAnalysisManager::kOutputContainer,
                                                                AliAnalysisManager::GetCommonFileName()));

      AliAnalysisDataContainer *cV0Tree = mgr->CreateContainer("cTreeMC",
                                                               TTree::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               AliAnalysisManager::GetCommonFileName());
      cV0Tree->SetSpecialOutput();
      mgr->ConnectOutput(taskExtractV0, 2, cV0Tree);
    }
  } else {
    if (bIsInfoAOD) {
      AliAnalysisTaskExtractV0AOD *taskExtractV0 = new AliAnalysisTaskExtractV0AOD("AliAnalysisTaskExtractV0AOD");
      if (wTriggerMask) taskExtractV0->SelectCollisionCandidates(wTriggerMask);

      taskExtractV0->SetIsNuclear(kTRUE);
      taskExtractV0->SetIsLowEnergyPP(kFALSE);
      taskExtractV0->SetUseOnTheFly(kFALSE);
      taskExtractV0->SetTriggerMask("kINT7");

      mgr->AddTask(taskExtractV0);
      mgr->ConnectInput(taskExtractV0,  0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(taskExtractV0, 1, mgr->CreateContainer("cListV0",
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              AliAnalysisManager::GetCommonFileName()));

      AliAnalysisDataContainer *cV0Tree = mgr->CreateContainer("cTree",
                                                               TTree::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               AliAnalysisManager::GetCommonFileName());
      cV0Tree->SetSpecialOutput();
      mgr->ConnectOutput(taskExtractV0, 2, cV0Tree);
    } else {
      AliAnalysisTaskExtractV0 *taskExtractV0 = new AliAnalysisTaskExtractV0("AliAnalysisTaskExtractV0");
      if (wTriggerMask) taskExtractV0->SelectCollisionCandidates(wTriggerMask);

      taskExtractV0->SetIsNuclear(kTRUE);
      taskExtractV0->SetINT7Trigger(kTRUE);
      taskExtractV0->SetUseOnTheFly(kFALSE);
      taskExtractV0->SetTakeAllTracks(kFALSE);
      taskExtractV0->SetCentralityEstimator(sCentEsti);
      taskExtractV0->SetLightWeightAnalysis(kFALSE);
      taskExtractV0->SetpAVertexSelection(kTRUE);
      taskExtractV0->SetSkipTrigger(kFALSE);
      taskExtractV0->SetTPCdEdxSelection(kFALSE);

      mgr->AddTask(taskExtractV0);
      mgr->ConnectInput(taskExtractV0,  0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput(taskExtractV0, 1, mgr->CreateContainer("clistV0",
                                                                TList::Class(),
                                                                AliAnalysisManager::kOutputContainer,
                                                                AliAnalysisManager::GetCommonFileName()));

      AliAnalysisDataContainer *cV0Tree = mgr->CreateContainer("cTree",
                                                               TTree::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               AliAnalysisManager::GetCommonFileName());
      cV0Tree->SetSpecialOutput();
      mgr->ConnectOutput(taskExtractV0, 2, cV0Tree);
    }
  }*/

  return kFALSE;
}
