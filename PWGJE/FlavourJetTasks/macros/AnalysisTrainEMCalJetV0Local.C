const TString sPeriod  = "LHC11a";
const TString sDataset = Form("%s.txt",sPeriod.Data());

const Bool_t bAOD    = kTRUE;
const Bool_t bMC     = kFALSE;
const Bool_t bOwnAOD = kTRUE;

const UInt_t wEvMask = AliVEvent::kMB; // kINT7

const Bool_t bMultSel = kFALSE;
const Double_t dMinMult = 0.;
const Double_t dMaxMult = 10.;

const Double_t dJetR = 0.4;
const Double_t dCutTrkPt = 0.15;
const Double_t dCutClusE = 0.30;
const Double_t dCutJetPt = 1.;

const Bool_t bJetBkg  = kFALSE;
const Bool_t bFullJet = kFALSE;
//=============================================================================

void AnalysisTrainEMCalJetV0Local()
{
  const TString sType(bFullJet ? "EMCAL" : "TPC");

  const AliAnalysisTaskEmcal::BeamType wBeamType(CheckBeamType());
  const Double_t dGhostArea(wBeamType==(AliAnalysisTaskEmcal::kpp) ? 0.01 : 0.005);

  const Bool_t bRUNII((sPeriod.Length()==6) && sPeriod.BeginsWith("LHC15"));

  AliTrackContainer::SetDefTrackCutsPeriod(sPeriod);
//=============================================================================

  TChain *chain = CreateChain();
  if (!chain) {
    ::Error("AnalysisTrainEMCalJetV0Local.C::AnalysisTrainEMCalJetV0Local", "Creating input chain failed!");
    return;
  }
//=============================================================================

  AliAnalysisManager *mgr = new AliAnalysisManager("AliAnaTrainEMCalJetV0Local", "Analysis Train EMCal Jet V0 Local");
//mgr->SetDebugLevel(3);
//=============================================================================

  if (bAOD) {
    AliAODInputHandler* pAODHandler = AliAnalysisTaskEmcal::AddAODHandler();
  } else {
    AliESDInputHandler* pESDHandler = AliAnalysisTaskEmcal::AddESDHandler();
  }

  if (bMC && (!bAOD)) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C");
    AliMCEventHandler *mctEH = AddMCHandler(kTRUE);
  }

  if (bOwnAOD) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODOutputHandler.C");
    AliAODHandler *aodH = AddAODOutputHandler();
    aodH->SetOutputFileName("AliAOD.Jets.root");
    aodH->SetCreateNonStandardAOD();
    aodH->SetFillAOD(kTRUE);
  }
//=============================================================================

  if (!bAOD) {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *taskPhysSel = AddTaskPhysicsSelection();
  }
//=============================================================================

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTask *taskPID = AddTaskPIDResponse(bMC, kTRUE, kTRUE, "2", kFALSE, "", kTRUE, kFALSE);
//=============================================================================

  if (bMultSel) {
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *taskMult = AddTaskMultSelection(kFALSE);
  }
//=============================================================================

  if (bFullJet && (!bAOD)) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
    taskCDB->SetFallBackToRaw(kTRUE);

//TODO hadronic correction
  }
//=============================================================================

  const TString sRho(bFullJet ? "Rho_Scaled" : "Rho");

  if (bJetBkg) {
    AliEmcalJetTask *taskKtJet = AliEmcalJetTask::AddTaskEmcalJet("usedefault", "",
      AliJetContainer::kt_algorithm, 0.4, AliJetContainer::kChargedJet,
      dCutTrkPt, dCutClusE, dGhostArea, AliJetContainer::pt_scheme,
      "Jet", 0., kFALSE, kFALSE);

    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C");
    AliAnalysisTaskRho* taskRho = AddTaskRhoNew("usedefault", "usedefault", sRho.Data(), 0.4);
    taskRho->SetExcludeLeadJets(2);

    if (bFullJet) {
      TString sFuncPath = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root";
      TString sFuncName = "LHC11h_HadCorr20_ClustersV2";
      taskRho->LoadRhoFunction(sFuncPath, sFuncName);
    }
  }
//=============================================================================

  AliEmcalJetTask *taskAktJet = AliEmcalJetTask::AddTaskEmcalJet("usedefault",
    bFullJet ? "usedefault" : "", AliJetContainer::antikt_algorithm, dJetR,
    bFullJet ? (AliJetContainer::kFullJet) : (AliJetContainer::kChargedJet),
    dCutTrkPt, dCutClusE, dGhostArea, AliJetContainer::pt_scheme,
    "Jet", dCutJetPt, kFALSE, kFALSE);

  if (bFullJet) taskAktJet->GetClusterContainer()->SetDefaultClusterEnergy(AliVCluster::kHadCorr);

  const TString sJets(taskAktJet->GetName());
//const TString sTrks(taskAktJet->GetParticleContainer()->GetName());
//const TString sClus(bFullJet ? taskAktJet->GetClusterContainer()->GetName() : "");
//=============================================================================

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskPicoV0Maker.C");
  AliAnalysisTaskSEPicoV0Maker *taskPicoV0 = AddTaskPicoV0Maker(bMC);

  taskPicoV0->SetUseAnaUtils(kTRUE);
  taskPicoV0->SetTriggerMask(wEvMask);
  if (bMultSel) taskPicoV0->SetMultRange(dMinMult, dMaxMult);
  if (wBeamType==(AliAnalysisTaskEmcal::kpp)) taskPicoV0->SetCollitionType(AliPicoBase::kPP);
  if (wBeamType==(AliAnalysisTaskEmcal::kpA)) taskPicoV0->SetCollitionType(AliPicoBase::kPA);
  if (wBeamType==(AliAnalysisTaskEmcal::kAA)) taskPicoV0->SetCollitionType(AliPicoBase::kAA);
//=============================================================================

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskEmcalJetV0Filter.C");
  AliAnalysisTaskEmcalJetV0Filter *taskJetV0Filter = AddTaskEmcalJetV0Filter("usedefault",
                                                                  bFullJet ? "usedefault" : "",
                                                                  bFullJet ? "usedefault" : "");

  taskJetV0Filter->SetHistoBins(600, 0., 300.);
  AliParticleContainer *cPar(taskJetV0Filter->GetParticleContainer());
  AliClusterContainer *clus(taskJetV0Filter->GetClusterContainer());
  AliJetContainer *cJet = taskJetV0Filter->AddJetContainer(sJets.Data(),sType.Data(),dJetR);

  if (cPar) {
    cPar->SetParticlePtCut(dCutTrkPt);
  }

  if (clus) {
    clus->SetClusECut(0.);
    clus->SetClusPtCut(0.);
    clus->SetClusNonLinCorrEnergyCut(0.);
    clus->SetClusHadCorrEnergyCut(0.3);
    clus->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  if (cJet) {
    cJet->SetPercAreaCut(0.6);
    cJet->ConnectParticleContainer(cPar);
    cJet->ConnectClusterContainer(clus);

    if (bJetBkg) {
      cJet->SetRhoName(sRho.Data());
    }
  }
//=============================================================================


  TIter next(mgr->GetTasks());
  AliAnalysisTaskSE *task(0x0);
  while ((taskSE = dynamic_cast<AliAnalysisTaskSE*>(next()))) {
    if (wEvMask) if (!(taskSE->InheritsFrom("AliPhysicsSelectionTask"))) {
      taskSE->SelectCollisionCandidates(wEvMask);
    }

    if (taskSE->InheritsFrom("AliAnalysisTaskEmcal")) {
      AliAnalysisTaskEmcal *taskEmcal = static_cast<AliAnalysisTaskEmcal*>(taskSE);

      taskEmcal->SetForceBeamType(wBeamType);
      taskEmcal->SetUseNewCentralityEstimation(bRUNII);

      taskEmcal->SetUseAliAnaUtils(kTRUE);
      taskEmcal->SetUseSPDTrackletVsClusterBG(kTRUE);
      if (bMultSel) taskEmcal->SetCentRange(dMinMult, dMaxMult);
    }
  }
//=============================================================================

  if (mgr->InitAnalysis()) { mgr->PrintStatus(); mgr->StartAnalysis("local",chain); }
  return;
}

//_____________________________________________________________________________
AliAnalysisTaskEmcal::BeamType CheckBeamType()
{
  if ((sPeriod=="LHC10h") || (sPeriod=="LHC11h") || (sPeriod=="LHC15o")) {
    return (AliAnalysisTaskEmcal::kAA);
  }

  if ((sPeriod=="LHC12g") || (sPeriod=="LHC13b") || (sPeriod=="LHC13c") ||
      (sPeriod=="LHC13d") || (sPeriod=="LHC13e") || (sPeriod=="LHC13f") ||
      (sPeriod=="LHC16q") || (sPeriod=="LHC16r") || (sPeriod=="LHC16s") ||
      (sPeriod=="LHC16t")) {
    return (AliAnalysisTaskEmcal::kpA);
  }

  return (AliAnalysisTaskEmcal::kpp);
}

//_____________________________________________________________________________
TChain *CreateChain()
{
  TChain *chain = 0;
  if (bAOD)
    chain = new TChain("aodTree");
  else
    chain = new TChain("esdTree");

  if (gSystem->AccessPathName(sDataset.Data())) {
    ::Error("AnalysisTrainEMCalJetV0Local.C::CreateChain","Dataset %s does not exist!",sDataset.Data());
    return NULL;
  }

  TString dataFile;
  ifstream dataList(sDataset.Data(), ios::in);
  while (!dataList.eof()) {
    dataFile.ReadLine(dataList,kFALSE);
    if (!dataFile.EndsWith(".root")) continue;
    if (!gSystem->AccessPathName(dataFile.Data())) chain->Add(dataFile.Data());
  } dataList.close();

  return chain;
}

//_____________________________________________________________________________
TChain *CreateAODFriendChain()
{
  if (!sDataset.EndsWith(".txt")) return 0;

  TChain *chain = new TChain("aodTree");
  TChain *cFrid = new TChain("aodTree");

  TString dataFile;
  ifstream dataList(sDataset.Data(), ios::in);
  while (!dataList.eof()) {
    dataFile.ReadLine(dataList,kFALSE);
    if (!dataFile.EndsWith("AliAOD.root")) continue;
    if (!gSystem->AccessPathName(dataFile.Data())) chain->Add(dataFile.Data());

    dataFile.ReplaceAll("AliAOD.root","AliAOD.VertexingHF.root");
    if (!gSystem->AccessPathName(dataFile.Data())) cFrid->Add(dataFile.Data());
  } dataList.close();

  chain->AddFriend(cFrid);
  return chain;
}
