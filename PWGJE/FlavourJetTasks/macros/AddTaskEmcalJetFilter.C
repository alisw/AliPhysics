const Bool_t bIsInfoAOD = kTRUE;
const Bool_t bIsPhysSel = kFALSE;
const Bool_t bIsCentSel = kFALSE;
const Bool_t bIsEvPnSel = kFALSE;
const Bool_t bIsRespPID = kFALSE;

const Bool_t bAnaJetR02 = kTRUE;
const Bool_t bAnaJetR03 = kTRUE;
const Bool_t bAnaJetR04 = kTRUE;
const Bool_t bAnaInfoMC = kFALSE;

const UInt_t wTriggerMask   = AliVEvent::kINT7;
//=============================================================================

const TString sPeriodIn = "LHC13b";
const TString sCentEsti = "V0A";  // "V0M"; "V0A"; "V0C"
const Double_t dCentMin =   0.;
const Double_t dCentMax = 100.;

const TString sInputTrkRD = (bIsInfoAOD ? "tracks"          : "Tracks");
const TString sInputClsRD = (bIsInfoAOD ? "caloClusters"    : "CaloClusters");
const TString sFilterTrks = (bIsInfoAOD ? "AODFilterTracks" : "ESDFilterTracks");

const TString sUsedTrksRD = "PicoTracks";
const TString sUsedClusRD = "CaloClustersCorr";

const TString sUsedTrksMC = "MCParticles";
const TString sUsedClusMC = "";

const TString sUsedRhoRD02  = "RhoRD02";  // "RhoRD"
const TString sUsedRhoRD03  = "RhoRD03";  // "RhoRD"
const TString sUsedRhoRD04  = "RhoRD04";  // "RhoRD"

const TString sUsedRhoMC02  = "";         // "RhoMC"
const TString sUsedRhoMC03  = "";         // "RhoMC"
const TString sUsedRhoMC04  = "";         // "RhoMC"

const TString sAnaType = "TPC";  // "TPC"; "EMCAL"; "USER"
const Int_t   nLeading = 0;      // 0: charged; 1: neutral; 2: both
//=============================================================================

const Int_t nJetAlgo = 1;  // 0: KT; 1: AKT
const Int_t nJetType = 1;  // 0: FullJet; 1: ChargedJet; 2: NeutralJet

const Double_t dJetPtCut   = 8.0;
const Double_t dJetAreaCut = 0.6;

const Double_t dTrkPtCut  = 0.15;
const Double_t dCluEnCut  = 0.30;
//=============================================================================

Bool_t AddTasksEmcalJetFilter()
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

  if (bAnaInfoMC && (!bIsInfoAOD)) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C");
    AliMCEventHandler *mctEH = AddMCHandler(kTRUE);
  }

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODOutputHandler.C");
  AliAODHandler *aodH = AddAODOutputHandler();
  aodH->SetOutputFileName("AliAOD.PicoJets.root");
  aodH->SetFillAOD(kTRUE);
  aodH->SetCreateNonStandardAOD();
  mgr->SetOutputEventHandler(aodH);
//=============================================================================

  if (bIsPhysSel) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
    AliEmcalPhysicsSelectionTask *taksPhysSel = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE, wTriggerMask, 5., 5., 10., kTRUE, -1, -1, -1, -1);
    if (bAnaInfoMC) {
      AliEmcalPhysicsSelection *pPhysSel = static_cast<AliEmcalPhysicsSelection*>(taksPhysSel->GetPhysicsSelection());
      if (!pPhysSel) return kTRUE; pPhysSel->SetAnalyzeMC();
    }
  }

  if (bIsCentSel) {
    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
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
    AliAnalysisTaskSE *taskRespPID = AddTaskPIDResponse(bAnaInfoMC);
    if (wTriggerMask)  taskRespPID->SelectCollisionCandidates(wTriggerMask);
  }
//=============================================================================

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
  taskCDB->SetFallBackToRaw(kTRUE);

  if (bIsInfoAOD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
    AliEmcalAodTrackFilterTask *taskTrkFilterAOD = AddTaskEmcalAodTrackFilter(sFilterTrks.Data(),
                                                                              sInputTrkRD.Data(),
                                                                              sPeriodIn.Data());
    if (wTriggerMask) taskTrkFilterAOD->SelectCollisionCandidates(wTriggerMask);
  } else {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
    AliEmcalEsdTrackFilterTask *taskTrkFilterESD = AddTaskEmcalEsdTrackFilter(sFilterTrks.Data(),
                                                                              Form("Hybrid_%s",sPeriodIn.Data()));
    if (wTriggerMask) taskTrkFilterESD->SelectCollisionCandidates(wTriggerMask);
    esdfilter->SetDoPropagation(kTRUE);
    esdfilter->SetDist(440.);
  }

  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  AliEmcalPicoTrackMaker *taskPicoTrackRD = AddTaskEmcalPicoTrackMaker(sUsedTrksRD.Data(),
                                                                       sFilterTrks.Data());
  if (wTriggerMask) taskPicoTrackRD->SelectCollisionCandidates(wTriggerMask);
  if (bAnaInfoMC)   taskPicoTrackRD->SetMC(kTRUE);
//=============================================================================

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");

  if (bAnaJetR02) {
    AliEmcalJetTask *taskAktRD02 = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                   ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                   nJetAlgo,
                                                   0.2,
                                                   nJetType,
                                                   dTrkPtCut,
                                                   dCluEnCut);
    if (wTriggerMask) taskAktRD02->SelectCollisionCandidates(wTriggerMask);

    if (!sUsedRhoRD02.IsNull()) {
      AliEmcalJetTask *taskKtRD02 = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                    ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                    0,  // KT
                                                    0.2,
                                                    nJetType);
      if (wTriggerMask) taskKtRD02->SelectCollisionCandidates(wTriggerMask);
      taskKtRD02->SetMinJetPt(0.);

      AliAnalysisTaskRhoSparse *taskRhoRD02 = AddTaskRhoSparse(taskKtRD02->GetName(),
                                                               taskAktRD02->GetName(),
                                                               ((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                               ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                               sUsedRhoRD02.Data(),
                                                               0.2,
                                                               sAnaType.Data(),
                                                               0.01,                 // jet area cut
                                                               0.15,                 // jet pT cut
                                                               0.,                   // EMC area
                                                               0x0,                  // scale fxn
                                                               0,                    // excluded leadings
                                                               kTRUE,                // output histogram
                                                               sUsedRhoRD02.Data(),  // task name
                                                               kTRUE);               // use CMS rho
      if (wTriggerMask) taskRhoRD02->SelectCollisionCandidates(wTriggerMask);
      taskRhoRD02->SetCentRange(dCentMin, dCentMax);
      taskRhoRD02->SetCentralityEstimator(sCentEsti.Data());
    }
  }

  if (bAnaJetR03) {
    AliEmcalJetTask *taskAktRD03 = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                   ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                   nJetAlgo,
                                                   0.3,
                                                   nJetType,
                                                   dTrkPtCut,
                                                   dCluEnCut);
    if (wTriggerMask) taskAktRD03->SelectCollisionCandidates(wTriggerMask);

    if (!sUsedRhoRD03.IsNull()) {
      AliEmcalJetTask *taskKtRD03 = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                    ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                    0,  // KT
                                                    0.3,
                                                    nJetType);
      if (wTriggerMask) taskKtRD03->SelectCollisionCandidates(wTriggerMask);
      taskKtRD03->SetMinJetPt(0.);

      AliAnalysisTaskRhoSparse *taskRhoRD03 = AddTaskRhoSparse(taskKtRD03->GetName(),
                                                               taskAktRD03->GetName(),
                                                               ((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                               ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                               sUsedRhoRD03.Data(),
                                                               0.3,
                                                               sAnaType.Data(),
                                                               0.01,                 // jet area cut
                                                               0.15,                 // jet pT cut
                                                               0.,                   // EMC area
                                                               0x0,                  // scale fxn
                                                               0,                    // excluded leadings
                                                               kTRUE,                // output histogram
                                                               sUsedRhoRD03.Data(),  // task name
                                                               kTRUE);               // use CMS rho
      if (wTriggerMask) taskRhoRD03->SelectCollisionCandidates(wTriggerMask);
      taskRhoRD03->SetCentRange(dCentMin, dCentMax);
      taskRhoRD03->SetCentralityEstimator(sCentEsti.Data());
    }
  }

  if (bAnaJetR04) {
    AliEmcalJetTask *taskAktRD04 = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                   ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                   nJetAlgo,
                                                   0.4,
                                                   nJetType,
                                                   dTrkPtCut,
                                                   dCluEnCut);
    if (wTriggerMask) taskAktRD04->SelectCollisionCandidates(wTriggerMask);

    if (!sUsedRhoRD04.IsNull()) {
      AliEmcalJetTask *taskKtRD04 = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                    ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                    0,  // KT
                                                    0.4,
                                                    nJetType);
      if (wTriggerMask) taskKtRD04->SelectCollisionCandidates(wTriggerMask);
      taskKtRD04->SetMinJetPt(0.);

      AliAnalysisTaskRhoSparse *taskRhoRD04 = AddTaskRhoSparse(taskKtRD04->GetName(),
                                                               taskAktRD04->GetName(),
                                                               ((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                               ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                               sUsedRhoRD04.Data(),
                                                               0.4,
                                                               sAnaType.Data(),
                                                               0.01,                 // jet area cut
                                                               0.15,                 // jet pT cut
                                                               0.,                   // EMC area
                                                               0x0,                  // scale fxn
                                                               0,                    // excluded leadings
                                                               kTRUE,                // output histogram
                                                               sUsedRhoRD04.Data(),  // task name
                                                               kTRUE);               // use CMS rho
      if (wTriggerMask) taskRhoRD04->SelectCollisionCandidates(wTriggerMask);
      taskRhoRD04->SetCentRange(dCentMin, dCentMax);
      taskRhoRD04->SetCentralityEstimator(sCentEsti.Data());
    }
  }
//=============================================================================

  if (bAnaInfoMC) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *taskPicoTrackMC = AddTaskMCTrackSelector(sUsedTrksMC.Data(),
                                                                      kFALSE,  // NK
                                                                      kTRUE);  // CH
    if (wTriggerMask) taskPicoTrackMC->SelectCollisionCandidates(wTriggerMask);
  }
//=============================================================================

  AliAnalysisTaskEmcalJetFilter *taskEmcalJetFilter = new AliAnalysisTaskEmcalJetFilter("AliAnalysisTaskEmcalJetFilter");
  if (wTriggerMask) taskEmcalJetFilter->SelectCollisionCandidates(wTriggerMask);

  taskEmcalJetFilter->SetCentRange(dCentMin, dCentMax);
  taskEmcalJetFilter->SetCentralityEstimator(sCentEsti.Data());
//taskEmcalJetFilter->SetForceBeamType(0);
//taskEmcalJetFilter->SetIsPythia(kTRUE);

  if (bAnaJetR02) {
    AliJetContainer *pContJetsRD02 = taskEmcalJetFilter->AddJetContainer(taskAktRD02->GetName(), sAnaType.Data(), 0.2);
    pContJetsRD02->SetPercAreaCut(dJetAreaCut);
    pContJetsRD02->SetJetPtCut(dJetPtCut);
    pContJetsRD02->SetRhoName(sUsedRhoRD02.Data());
//  pContJetsRD02->SetLocalRhoName();
    pContJetsRD02->SetLeadingHadronType(nLeading);
    pContJetsRD02->ConnectParticleContainer(taskEmcalJetFilter->AddParticleContainer((nJetType!=2) ? sUsedTrksRD.Data() : ""));
    pContJetsRD02->ConnectClusterContainer( taskEmcalJetFilter->AddClusterContainer( (nJetType!=1) ? sUsedClusRD.Data() : ""));
    pContJetsRD02->SetNameTitle("AliEMcalJetContainerRD02", "AliEMcalJetContainerRD02");
    taskEmcalJetFilter->SetNameJetRD02("AliEMcalJetContainerRD02");
  }

  if (bAnaJetR03) {
    AliJetContainer *pContJetsRD03 = taskEmcalJetFilter->AddJetContainer(taskAktRD03->GetName(), sAnaType.Data(), 0.3);
    pContJetsRD03->SetPercAreaCut(dJetAreaCut);
    pContJetsRD03->SetJetPtCut(dJetPtCut);
    pContJetsRD03->SetRhoName(sUsedRhoRD03.Data());
//  pContJetsRD03->SetLocalRhoName();
    pContJetsRD03->SetLeadingHadronType(nLeading);
    pContJetsRD03->ConnectParticleContainer(taskEmcalJetFilter->AddParticleContainer((nJetType!=2) ? sUsedTrksRD.Data() : ""));
    pContJetsRD03->ConnectClusterContainer( taskEmcalJetFilter->AddClusterContainer( (nJetType!=1) ? sUsedClusRD.Data() : ""));
    pContJetsRD03->SetNameTitle("AliEMcalJetContainerRD03", "AliEMcalJetContainerRD03");
    taskEmcalJetFilter->SetNameJetRD03("AliEMcalJetContainerRD03");
  }

  if (bAnaJetR04) {
    AliJetContainer *pContJetsRD04 = taskEmcalJetFilter->AddJetContainer(taskAktRD04->GetName(), sAnaType.Data(), 0.4);
    pContJetsRD04->SetPercAreaCut(dJetAreaCut);
    pContJetsRD04->SetJetPtCut(dJetPtCut);
    pContJetsRD04->SetRhoName(sUsedRhoRD04.Data());
//  pContJetsRD04->SetLocalRhoName();
    pContJetsRD04->SetLeadingHadronType(nLeading);
    pContJetsRD04->ConnectParticleContainer(taskEmcalJetFilter->AddParticleContainer((nJetType!=2) ? sUsedTrksRD.Data() : ""));
    pContJetsRD04->ConnectClusterContainer( taskEmcalJetFilter->AddClusterContainer( (nJetType!=1) ? sUsedClusRD.Data() : ""));
    pContJetsRD04->SetNameTitle("AliEMcalJetContainerRD04", "AliEMcalJetContainerRD04");
    taskEmcalJetFilter->SetNameJetRD04("AliEMcalJetContainerRD04");
  }

  mgr->AddTask(taskEmcalJetFilter);
  mgr->ConnectInput(taskEmcalJetFilter,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEmcalJetFilter, 1, mgr->CreateContainer("listEmcalJet",
                                                                 TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 AliAnalysisManager::GetCommonFileName()));

  return kFALSE;
}
