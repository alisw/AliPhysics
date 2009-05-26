//
// Macro to create the full analysis manager for Resonances
//

static Double_t  cov11 = 2;
static Double_t  cov22 = 2;
static Double_t  cov33 = 0.5;
static Double_t  cov44 = 0.5;
static Double_t  cov55 = 2;
static Double_t  nSigmaToVertex = 4;
static Double_t  dcaToVertex = 3.0;
static Double_t  maxChi2PerClusterTPC = 3.5;
static Bool_t    requireTPCRefit = kTRUE;
static Bool_t    requireSigmaToVertex = kTRUE;
static Bool_t    acceptKinkDaughters = kFALSE;
static Int_t     minNClustersTPC = 50;

Bool_t AddAnalysisTaskRsn(const char *outFile = "rsn.root", Bool_t useKine = kTRUE)
{
  // Retrieve analysis manager.
  // Since it is usually created in the steering macro,
  // then we don't ever need to initialize a new one
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  // no need to  create handlers for input and MonteCarlo
  // Done in the general steering macro
  //  AliMCEventHandler  *mcHandler  = new AliMCEventHandler();
  //  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  //AliAODHandler      *aodHandler = new AliAODHandler();
  //aodHandler->SetOutputFileName("AliAOD.root");
  //  if (useKine) mgr->SetMCtruthEventHandler(mcHandler);
  //  mgr->SetInputEventHandler(esdHandler);
  //mgr->SetOutputEventHandler(aodHandler);

  // initialize task
  AliRsnAnalysisSE *task = new AliRsnAnalysisSE("AliRsnAnalysisSE");

  // set prior probabilities for PID
  task->SetPriorProbability(AliPID::kElectron, 0.02);
  task->SetPriorProbability(AliPID::kMuon,     0.02);
  task->SetPriorProbability(AliPID::kPion,     0.83);
  task->SetPriorProbability(AliPID::kKaon,     0.07);
  task->SetPriorProbability(AliPID::kProton,   0.06);
  task->DumpPriors();

  // initialize analysis manager with pairs from config
  AliRsnAnalysisManager *anaMgr = task->GetAnalysisManager("MyAnalysisSE");
  anaMgr->Add(CreatePairsNoPID("PHI_NoPID_0"   , 333, AliPID::kKaon, AliPID::kKaon, 10000.0));
  anaMgr->Add(CreatePairsNoPID("PHI_NoPID_10"  , 333, AliPID::kKaon, AliPID::kKaon,     0.1));
  anaMgr->Add(CreatePairsNoPID("PHI_NoPID_20"  , 333, AliPID::kKaon, AliPID::kKaon,     0.2));
  anaMgr->Add(CreatePairsPID  ("PHI_PID"       , 333, AliPID::kKaon, AliPID::kKaon));
  anaMgr->Add(CreatePairsNoPID("KSTAR_NoPID_0" , 313, AliPID::kKaon, AliPID::kPion, 10000.0));
  anaMgr->Add(CreatePairsNoPID("KSTAR_NoPID_10", 313, AliPID::kKaon, AliPID::kPion,     0.1));
  anaMgr->Add(CreatePairsNoPID("KSTAR_NoPID_20", 313, AliPID::kKaon, AliPID::kPion,     0.2));
  anaMgr->Add(CreatePairsPID  ("KSTAR_PID"     , 313, AliPID::kKaon, AliPID::kPion));

  // setup cuts
  AliESDtrackCuts *esdCuts = new AliESDtrackCuts;
  esdCuts->SetMaxCovDiagonalElements(cov11, cov22, cov33, cov44, cov55);
  esdCuts->SetRequireSigmaToVertex(requireSigmaToVertex);
  if (requireSigmaToVertex) esdCuts->SetMaxNsigmaToVertex(nSigmaToVertex);
  else
  {
    esdCuts->SetDCAToVertexZ(dcaToVertex);
    esdCuts->SetDCAToVertexXY(dcaToVertex);
  }
  esdCuts->SetRequireTPCRefit(requireTPCRefit);
  esdCuts->SetAcceptKinkDaughters(acceptKinkDaughters);
  esdCuts->SetMinNClustersTPC(minNClustersTPC);
  esdCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
  task->SetESDtrackCuts(esdCuts);

  // add the task to manager
  mgr->AddTask(task);

  // initialize container for the output
  // (input container is common)
  AliAnalysisDataContainer *out = mgr->CreateContainer("RSN", TList::Class(), AliAnalysisManager::kOutputContainer, outFile);

  // connect input and output slots for this task
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, out);
}

void AddRsnFunctionsToPair(AliRsnPair *pair)
{
//
// In general, for all processed pairs in one analysis the same functions are computed.
// Then, they are defined separately here and added in the same way to all pairs.
//

  // define histogram templates
  AliRsnHistoDef *hdIM   = new AliRsnHistoDef(1000,  0.0,   2.0);
  AliRsnHistoDef *hdPt   = new AliRsnHistoDef(  10,  0.0,  10.0);
  AliRsnHistoDef *hdMult = new AliRsnHistoDef(   8,  0.0, 200.0);

  // define functions:
  // -- one computed w.r. to Pt and Mult
  // -- one computed w.r. to Eta and Mult
  AliRsnFunction *fcnIM1 = new AliRsnFunction(AliRsnFunction::kInvMass, AliRsnFunction::kBinPairPt, AliRsnFunction::kBinEventMult, hdIM, hdPt, hdMult);

  // add to pair
  pair->AddFunction(fcnIM1);
}

AliRsnPairManager* CreatePairsNoPID
(
  const char            *pairMgrName,    // name for the pair manager
  Int_t                  resonancePDG,   // PDG code of resonance (for true pairs)
  AliPID::EParticleType  type1,          // PID of one member of decay (+)
  AliPID::EParticleType  type2,          // PID of other member of decay (-)
  Double_t               bbCut = 10000.0 // cut on Bethe-Bloch
)
{
//
// Creates an AliRsnPairMgr for a specified resonance, which contains:
// - signal (inv. mass)
// - event mixing (inv. mass)
// - like-sign (inv. mass)
// - true pairs (inv. mass, resolution)
//
// For all pairs, a binning in Pt and Eta is provided, and a cut in multiplicity
// which defines a multiplicity bin where the analysis is computed.
//
// Arguments define how the pair manager must be built, and are explained above
//

  AliRsnPairManager  *pairMgr  = new AliRsnPairManager(pairMgrName);

  // === PAIR DEFINITIONS =========================================================================

  // if particle #1 and #2 are different, two histograms must be built
  // for each scheme (signal, true, mixed, like-sign) exchanging both particles and signs
  Int_t i, j, nArray = 1;
  if (type1 != type2) nArray = 2;

  AliRsnPairDef *defUnlike[2] = {0, 0};
  AliRsnPairDef *defLikePP[2] = {0, 0};
  AliRsnPairDef *defLikeMM[2] = {0, 0};

  defUnlike[0] = new AliRsnPairDef(type1, '+', type2, '-', resonancePDG);
  defLikePP[0] = new AliRsnPairDef(type1, '+', type2, '+', resonancePDG);
  defLikeMM[0] = new AliRsnPairDef(type1, '-', type2, '-', resonancePDG);

  defUnlike[1] = new AliRsnPairDef(type2, '+', type1, '-', resonancePDG);
  defLikePP[1] = new AliRsnPairDef(type2, '+', type1, '+', resonancePDG);
  defLikeMM[1] = new AliRsnPairDef(type2, '-', type1, '-', resonancePDG);

  // === PAIR ANALYSIS ENGINES ====================================================================

  // define null (dummy) objects and initialize only the ones which are needed,
  // depending again on particle types;
  // array is organized as follows:
  // [0] - true pairs
  // [1] - signal
  // [2] - like PP
  // [3] - like MM
  AliRsnPair *noPID[2][4] = {0,0,0,0,0,0,0,0};

  for (i = 0; i < nArray; i++) {
    noPID[i][0] = new AliRsnPair(AliRsnPair::kNoPID, defUnlike[i]);
    noPID[i][1] = new AliRsnPair(AliRsnPair::kNoPID, defUnlike[i]);
    noPID[i][2] = new AliRsnPair(AliRsnPair::kNoPID, defLikePP[i]);
    noPID[i][3] = new AliRsnPair(AliRsnPair::kNoPID, defLikeMM[i]);
  }

  // === CUTS =====================================================================================

  // cuts for tracks:
  // -- Bethe-Bloch & required kaon PID
  AliRsnCutStd *cutReqKaon = new AliRsnCutStd("cutReqKaon", AliRsnCutStd::kRequiredPID, (Int_t)AliPID::kKaon);
  AliRsnCutBetheBloch *cutKaonBB = new AliRsnCutBetheBloch("cutKaon", bbCut, AliPID::kKaon);
  cutKaonBB->SetCalibConstant(0, 0.76176e-1);
  cutKaonBB->SetCalibConstant(1, 10.632);
  cutKaonBB->SetCalibConstant(2, 0.13279e-4);
  cutKaonBB->SetCalibConstant(3, 1.8631);
  cutKaonBB->SetCalibConstant(4, 1.9479);

  // cuts on pairs:
  // -- true daughters of a phi resonance (only for true pairs histogram)cutSetPairTrue->AddCut(cutTrue);
  AliRsnCutStd *cutTruePair = new AliRsnCutStd("cutTrue", AliRsnCutStd::kTruePair, resonancePDG);

  // cut set definition for all pairs
  AliRsnCutSet *cutSetParticle = new AliRsnCutSet("trackCuts");
  cutSetParticle->AddCut(cutKaonBB);
  cutSetParticle->AddCut(cutReqKaon);
  cutSetParticle->SetCutScheme("(cutReqKaon&cutKaonBB)|(!cutReqKaon)");

  // cut set definition for true pairs
  AliRsnCutSet *cutSetPairTrue = new AliRsnCutSet("truePairs");
  cutSetPairTrue->AddCut(cutTruePair);
  cutSetPairTrue->SetCutScheme("cutTrue");

  // cut manager for all pairs
  AliRsnCutMgr *cutMgrAll = new AliRsnCutMgr("std", "All");
  cutMgrAll->SetCutSet(AliRsnCut::kParticle, cutSetParticle);

  // cut manager for all pairs
  AliRsnCutMgr *cutMgrTrue = new AliRsnCutMgr("true", "True");
  cutMgrTrue->SetCutSet(AliRsnCut::kParticle, cutSetParticle);
  cutMgrTrue->SetCutSet(AliRsnCut::kPair, cutSetPairTrue);

  for (i = 0; i < nArray; i++) {
    noPID[i][0]->SetCutMgr(cutMgrTrue);
    for (j = 1; j < 4; j++) {
      noPID[i][j]->SetCutMgr(cutMgrAll);
    }
  }

  // === FUNCTIONS ================================================================================

  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 4; j++) {
      AddRsnFunctionsToPair(noPID[i][j]);
    }
  }

  // === ADD TO PAIR MANAGER ======================================================================

  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 4; j++) {
      pairMgr->AddPair(noPID[i][j]);
    }
  }

  return pairMgr;
}

AliRsnPairManager* CreatePairsPID
(
  const char            *pairMgrName,    // name for the pair manager
  Int_t                  resonancePDG,   // PDG code of resonance (for true pairs)
  AliPID::EParticleType  type1,          // PID of one member of decay (+)
  AliPID::EParticleType  type2           // PID of other member of decay (-)
)
{
//
// Creates an AliRsnPairMgr for a specified resonance, which contains:
// - signal (inv. mass)
// - event mixing (inv. mass)
// - like-sign (inv. mass)
// - true pairs (inv. mass, resolution)
//
// For all pairs, a binning in Pt and Eta is provided, and a cut in multiplicity
// which defines a multiplicity bin where the analysis is computed.
//
// Arguments define how the pair manager must be built, and are explained above
//

  AliRsnPairManager  *pairMgr  = new AliRsnPairManager(pairMgrName);

  // === PAIR DEFINITIONS =========================================================================

  // if particle #1 and #2 are different, two histograms must be built
  // for each scheme (signal, true, mixed, like-sign) exchanging both particles and signs
  Int_t i, j, nArray = 1;
  if (type1 != type2) nArray = 2;

  AliRsnPairDef *defUnlike[2] = {0, 0};
  AliRsnPairDef *defLikePP[2] = {0, 0};
  AliRsnPairDef *defLikeMM[2] = {0, 0};

  defUnlike[0] = new AliRsnPairDef(type1, '+', type2, '-', resonancePDG);
  defLikePP[0] = new AliRsnPairDef(type1, '+', type2, '+', resonancePDG);
  defLikeMM[0] = new AliRsnPairDef(type1, '-', type2, '-', resonancePDG);

  defUnlike[1] = new AliRsnPairDef(type2, '+', type1, '-', resonancePDG);
  defLikePP[1] = new AliRsnPairDef(type2, '+', type1, '+', resonancePDG);
  defLikeMM[1] = new AliRsnPairDef(type2, '-', type1, '-', resonancePDG);

  // === PAIR ANALYSIS ENGINES ====================================================================

  // define null (dummy) objects and initialize only the ones which are needed,
  // depending again on particle types;
  // array is organized as follows:
  // [0] - true pairs
  // [1] - signal
  // [2] - like PP
  // [3] - like MM
  AliRsnPair *perfectPID[2][4]   = {0,0,0,0,0,0,0,0};
  AliRsnPair *realisticPID[2][4] = {0,0,0,0,0,0,0,0};

  for (i = 0; i < nArray; i++) {
    perfectPID[i][0] = new AliRsnPair(AliRsnPair::kPerfectPID, defUnlike[i]);
    perfectPID[i][1] = new AliRsnPair(AliRsnPair::kPerfectPID, defUnlike[i]);
    perfectPID[i][2] = new AliRsnPair(AliRsnPair::kPerfectPID, defLikePP[i]);
    perfectPID[i][3] = new AliRsnPair(AliRsnPair::kPerfectPID, defLikeMM[i]);

    realisticPID[i][0] = new AliRsnPair(AliRsnPair::kRealisticPID, defUnlike[i]);
    realisticPID[i][1] = new AliRsnPair(AliRsnPair::kRealisticPID, defUnlike[i]);
    realisticPID[i][2] = new AliRsnPair(AliRsnPair::kRealisticPID, defLikePP[i]);
    realisticPID[i][3] = new AliRsnPair(AliRsnPair::kRealisticPID, defLikeMM[i]);
  }

  // === CUTS =====================================================================================

  // cuts for tracks:
  // -- nothing

  // cuts on pairs:
  // -- true daughters of a phi resonance (only for true pairs histogram)cutSetPairTrue->AddCut(cutTrue);
  AliRsnCutStd *cutTruePair = new AliRsnCutStd("cutTrue", AliRsnCutStd::kTruePair, resonancePDG);

  // cut set definition for true pairs
  AliRsnCutSet *cutSetPairTrue = new AliRsnCutSet("truePairs");
  cutSetPairTrue->AddCut(cutTruePair);
  cutSetPairTrue->SetCutScheme("cutTrue");

  // cut manager for all pairs
  AliRsnCutMgr *cutMgrAll = new AliRsnCutMgr("std", "All");

  // cut manager for all pairs
  AliRsnCutMgr *cutMgrTrue = new AliRsnCutMgr("true", "True");
  cutMgrTrue->SetCutSet(AliRsnCut::kPair, cutSetPairTrue);

  for (i = 0; i < nArray; i++) {
    perfectPID[i][0]->SetCutMgr(cutMgrTrue);
    realisticPID[i][0]->SetCutMgr(cutMgrTrue);
    for (j = 1; j < 4; j++) {
      perfectPID[i][j]->SetCutMgr(cutMgrAll);
      realisticPID[i][j]->SetCutMgr(cutMgrAll);
    }
  }

  // === FUNCTIONS ================================================================================

  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 4; j++) {
      AddRsnFunctionsToPair(perfectPID[i][j]);
      AddRsnFunctionsToPair(realisticPID[i][j]);
    }
  }

  // === ADD TO PAIR MANAGER ======================================================================

  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 4; j++) {
      pairMgr->AddPair(perfectPID[i][j]);
      pairMgr->AddPair(realisticPID[i][j]);
    }
  }

  return pairMgr;
}

