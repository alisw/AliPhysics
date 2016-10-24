const Bool_t bIsInfoAOD = kTRUE;
const Bool_t bIsPhysSel = kFALSE;
const Bool_t bIsCentSel = kFALSE;
const Bool_t bIsEvPnSel = kFALSE;
const Bool_t bIsRespPID = kTRUE;

const Bool_t bAnaInfoMC = kTRUE;
const Bool_t bAnaJetsMC = kFALSE;

const UInt_t wTriggerMask   = AliVEvent::kINT7;
const UInt_t wCollisionType = AliPicoHeaderCJ::kPA;
//=============================================================================

const TString sPeriodIn = "LHC13b4";
const TString sCentEsti = "V0A";  // "V0M"; "V0A"; "V0C"
const Double_t dCentMin =   0.;
const Double_t dCentMax = 100.;

const TString sInputTrkRD = (bIsInfoAOD ? "tracks"       : "Tracks");
const TString sInputClsRD = (bIsInfoAOD ? "caloClusters" : "CaloClusters");

const TString sUsedTrksRD = "PicoTracks";
const TString sUsedClusRD = "CaloClustersCorr";

const TString sUsedTrksMC = "MCParticles";
const TString sUsedClusMC = "";

const TString sUsedRhoRD  = "RhoRD";  // "RhoRD"
const TString sUsedRhoMC  = "";       // "RhoMC"

const TString sAnaType = "TPC";  // "TPC"; "EMCAL"; "USER"
const Int_t   nLeading = 0;      // 0: charged; 1: neutral; 2: both
//=============================================================================

const Double_t dV0Cuts[] = {
  0.5,    // min V0 radius
  200.,   // max V0 radius
  1.,     // max V0 daus DCA
  0.06,   // min tr DCA to PV
  70.,    // min tr X rows TPC
  0.8,    // min tr X/F rows TPC
  5.,     // max Ka PID sigma TPC
  0.97,   // min Ka cosPA
  20.,    // max Ka ctau
  0.005,  // min Ka deltaM
  5.,     // max La PID sigma TPC
  0.995,  // min La cosPA
  30.,    // max La ctau
  0.01    // min La deltaM
};

const Double_t dV0EtaMin = -0.75;
const Double_t dV0EtaMax =  0.75;

const TString sFileV0InvM = "/global/scratch2/sd/xmzhang/local/AnaCorrJets20140301/outputs/JE_V0A_000_010/AnalysisOutputs_FitInvMrd.root";
//=============================================================================

const Int_t nJetAlgo = 1;  // 0: KT; 1: AKT
const Int_t nJetType = 1;  // 0: FullJet; 1: ChargedJet; 2: NeutralJet

const Double_t dTrkPtCut  = 0.15;
const Double_t dCluEnCut  = 0.30;

const Double_t dJetRadius  = 0.4;
const Double_t dJetPtCut   = 10.;
const Double_t dJetAreaCut = 0.6;

const Double_t dJetEtaMin = dV0EtaMin + dJetRadius;
const Double_t dJetEtaMax = dV0EtaMax - dJetRadius;
//=============================================================================

Bool_t AddTasksEmcalJetV0CF()
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
//=============================================================================

  if (bIsPhysSel) {
/*  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
    AliEmcalPhysicsSelectionTask *taksPhysSel = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE, wTriggerMask, 5., 5., 10., kTRUE, -1, -1, -1, -1);
    if (bAnaInfoMC) {
      AliEmcalPhysicsSelection *pPhysSel = static_cast<AliEmcalPhysicsSelection*>(taksPhysSel->GetPhysicsSelection());
      if (!pPhysSel) return kTRUE; pPhysSel->SetAnalyzeMC();
    }*/

    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *taskPhysSel = AddTaskPhysicsSelection(bAnaInfoMC);
  }

  if (bIsCentSel) {
    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentSel = AddTaskCentrality(kTRUE, bIsInfoAOD);
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
  if (wTriggerMask)  taskCDB->SelectCollisionCandidates(wTriggerMask);
  taskCDB->SetFallBackToRaw(kTRUE);

  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  AliEmcalPicoTrackMaker *taskPicoTrackRD = AddTaskEmcalPicoTrackMaker(sUsedTrksRD.Data(),
                                                                       sInputTrkRD.Data(),
                                                                       sPeriodIn.Data());
  if (wTriggerMask) taskPicoTrackRD->SelectCollisionCandidates(wTriggerMask);
  if (bAnaInfoMC)   taskPicoTrackRD->SetMC(kTRUE);

  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  AliEmcalJetTask *taskAktRD = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                               ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                               nJetAlgo,
                                               dJetRadius,
                                               nJetType,
                                               dTrkPtCut,
                                               dCluEnCut);
  if (wTriggerMask) taskAktRD->SelectCollisionCandidates(wTriggerMask);

  if (!sUsedRhoRD.IsNull()) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
    AliEmcalJetTask *taskKtRD = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                0,  // KT
                                                dJetRadius,
                                                nJetType);
    if (wTriggerMask) taskKtRD->SelectCollisionCandidates(wTriggerMask);
    taskKtRD->SetMinJetPt(0.);

    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
    AliAnalysisTaskRhoSparse *taskRhoRD = AddTaskRhoSparse(taskKtRD->GetName(),
                                                           taskAktRD->GetName(),
                                                           ((nJetType!=2) ? sUsedTrksRD.Data() : ""),
                                                           ((nJetType!=1) ? sUsedClusRD.Data() : ""),
                                                           sUsedRhoRD.Data(),
                                                           dJetRadius,
                                                           sAnaType.Data(),
                                                           0.01,               // jet area cut
                                                           0.15,               // jet pT cut
                                                           0.,                 // EMC area
                                                           0x0,                // scale fxn
                                                           0,                  // excluded leadings
                                                           kTRUE,              // output histogram
                                                           sUsedRhoRD.Data(),  // task name
                                                           kTRUE);             // use CMS rho
    if (wTriggerMask) taskRhoRD->SelectCollisionCandidates(wTriggerMask);
    taskRhoRD->SetCentRange(dCentMin, dCentMax);
    taskRhoRD->SetCentralityEstimator(sCentEsti.Data());
  }
//=============================================================================

  if (bAnaJetsMC) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *taskPicoTrackMC = AddTaskMCTrackSelector(sUsedTrksMC.Data(),
                                                                      kFALSE,  // NK
                                                                      kTRUE);  // CH
    if (wTriggerMask) taskPicoTrackMC->SelectCollisionCandidates(wTriggerMask);

    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
    AliEmcalJetTask *taskAktMC = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksMC.Data() : ""),
                                                 ((nJetType!=1) ? sUsedClusMC.Data() : ""),
                                                 nJetAlgo,
                                                 dJetRadius,
                                                 nJetType,
                                                 dTrkPtCut,
                                                 dCluEnCut);
    if (wTriggerMask) taskAktMC->SelectCollisionCandidates(wTriggerMask);

    if (!sUsedRhoMC.IsNull()) {
      gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
      AliEmcalJetTask *taskKtMC = AddTaskEmcalJet(((nJetType!=2) ? sUsedTrksMC.Data() : ""),
                                                  ((nJetType!=1) ? sUsedClusMC.Data() : ""),
                                                  0,  // KT
                                                  dJetRadius,
                                                  nJetType);
      if (wTriggerMask) taskKtMC->SelectCollisionCandidates(wTriggerMask);
      taskKtMC->SetMinJetPt(0.);

      gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");
      AliAnalysisTaskRhoSparse *taskRhoMC = AddTaskRhoSparse(taskKtMC->GetName(),
                                                             taskAktMC->GetName(),
                                                             ((nJetType!=2) ? sUsedTrksMC.Data() : ""),
                                                             ((nJetType!=1) ? sUsedClusMC.Data() : ""),
                                                             sUsedRhoMC.Data(),
                                                             dJetRadius,
                                                             sAnaType.Data(),
                                                             0.01,               // jet area cut
                                                             0.15,               // jet pT cut
                                                             0.,                 // EMC area
                                                             0x0,                // scale fxn
                                                             0,                  // excluded leadings
                                                             kTRUE,              // output histogram
                                                             sUsedRhoMC.Data(),  // task name
                                                             kTRUE);             // use CMS rho
      if (wTriggerMask) taskRhoMC->SelectCollisionCandidates(wTriggerMask);
      taskRhoMC->SetCentRange(dCentMin, dCentMax);
      taskRhoMC->SetCentralityEstimator(sCentEsti.Data());
    }
  }
//=============================================================================

  AliAnalysisTaskSEPicoV0MakerMC *taskPicoV0Maker = new AliAnalysisTaskSEPicoV0MakerMC("AliAnalysisTaskSEPicoV0MakerMC");
  if (wTriggerMask) {
    taskPicoV0Maker->SetTriggerMask(wTriggerMask);
    taskPicoV0Maker->SelectCollisionCandidates(wTriggerMask);
  }

  taskPicoV0Maker->SetCollitionType(wCollisionType);
  taskPicoV0Maker->SetCentralityEstimator(sCentEsti.Data());

//taskPicoV0Maker->SetRefitV0ESD();
//taskPicoV0Maker->SetSkipFastOnly();
//taskPicoV0Maker->SetDMPjetMC();

  taskPicoV0Maker->SetV0Cuts(dV0Cuts);
  taskPicoV0Maker->SetDauEtaRange(-0.8, 0.8);

  mgr->AddTask(taskPicoV0Maker);
  mgr->ConnectInput(taskPicoV0Maker,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskPicoV0Maker, 1, mgr->CreateContainer("listPicoV0MakerMC",
                                                              TList::Class(),
                                                              AliAnalysisManager::kOutputContainer,
                                                              AliAnalysisManager::GetCommonFileName()));
//=============================================================================

  TFile *file = TFile::Open(sFileV0InvM.Data(), "READ");
  TH1D *hKshortInvM = (TH1D*)((TList*)file->Get("listFitInvMrd_Kshort_Default"))->FindObject("hFitPtInvM"); hKshortInvM->SetName("hKshortInvM");
  TH1D *hLambdaInvM = (TH1D*)((TList*)file->Get("listFitInvMrd_Lambda_Default"))->FindObject("hFitPtInvM"); hLambdaInvM->SetName("hLambdaInvM");
  TH1D *hAntiLaInvM = (TH1D*)((TList*)file->Get("listFitInvMrd_AntiLa_Default"))->FindObject("hFitPtInvM"); hAntiLaInvM->SetName("hAntiLaInvM");
  file->Close();

  AliAnalysisTaskEmcalJetV0CF *taskEmcalJetV0Filter = new AliAnalysisTaskEmcalJetV0CF("AliAnalysisTaskEmcalJetV0CF");
  if (wTriggerMask) taskEmcalJetV0Filter->SelectCollisionCandidates(wTriggerMask);

  taskEmcalJetV0Filter->SetCentRange(dCentMin, dCentMax);
  taskEmcalJetV0Filter->SetCentralityEstimator(sCentEsti.Data());

  taskEmcalJetV0Filter->SetHistoKshortInvM(hKshortInvM);
  taskEmcalJetV0Filter->SetHistoLambdaInvM(hLambdaInvM);
  taskEmcalJetV0Filter->SetHistoAntiLaInvM(hAntiLaInvM);
  taskEmcalJetV0Filter->SetV0EtaRange(dV0EtaMin, dV0EtaMax);

//taskEmcalJetV0Filter->SetForceBeamType(0);
//taskEmcalJetV0Filter->SetIsPythia(kTRUE);

  AliJetContainer *pContJetsRD = taskEmcalJetV0Filter->AddJetContainer(taskAktRD->GetName(), "USER", dJetRadius);
  pContJetsRD->SetPercAreaCut(dJetAreaCut);
  pContJetsRD->SetJetPtCut(dJetPtCut);
  pContJetsRD->SetJetEtaLimits(dJetEtaMin, dJetEtaMax);
  pContJetsRD->SetRhoName(sUsedRhoRD.Data());
//pContJetsRD->SetLocalRhoName();
  pContJetsRD->SetLeadingHadronType(nLeading);
  pContJetsRD->ConnectParticleContainer(taskEmcalJetV0Filter->AddParticleContainer((nJetType!=2) ? sUsedTrksRD.Data() : ""));
  pContJetsRD->ConnectClusterContainer( taskEmcalJetV0Filter->AddClusterContainer( (nJetType!=1) ? sUsedClusRD.Data() : ""));

  if (bAnaJetsMC) {
    AliJetContainer *pContJetsMC = taskEmcalJetV0Filter->AddJetContainer(taskAktMC->GetName(), sAnaType.Data(), dJetRadius);
    pContJetsMC->SetPercAreaCut(dJetAreaCut);
    pContJetsMC->SetJetPtCut(dJetPtCut);
    pContJetsMC->SetJetEtaLimits(dJetEtaMin, dJetEtaMax);
    pContJetsMC->SetRhoName(sUsedRhoMC.Data());
//  pContJetsMC->SetLocalRhoName();
    pContJetsMC->SetLeadingHadronType(nLeading);
    pContJetsMC->ConnectParticleContainer(taskEmcalJetV0Filter->AddParticleContainer((nJetType!=2) ? sUsedTrksMC.Data() : ""));
    pContJetsMC->ConnectClusterContainer( taskEmcalJetV0Filter->AddClusterContainer( (nJetType!=1) ? sUsedClusMC.Data() : ""));
  }

  mgr->AddTask(taskEmcalJetV0Filter);
  mgr->ConnectInput(taskEmcalJetV0Filter,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEmcalJetV0Filter, 1, mgr->CreateContainer("listEmcalJetTask",
                                                                   TList::Class(),
                                                                   AliAnalysisManager::kOutputContainer,
                                                                   AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectOutput(taskEmcalJetV0Filter, 2, mgr->CreateContainer("listEmcalJetV0CF",
                                                                   TList::Class(),
                                                                   AliAnalysisManager::kOutputContainer,
                                                                   AliAnalysisManager::GetCommonFileName()));

  return kFALSE;
}
