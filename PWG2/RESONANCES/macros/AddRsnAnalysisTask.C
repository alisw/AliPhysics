//=============================================================================
//
// *** AliRsnAnalysisTask.C ***
//
// This macro initialize a complete AnalysisTask object for resonances.
// All the possibilities are made available through enumerations,
// with the same style of the Config.C file used for simulation.
//
//=============================================================================



static AliRsnAnalysisTaskSEBase::EInputType inputType      = AliRsnAnalysisTaskSEBase::kESDMC;
static const char*                          outputFileName = "rsn.root";

static Bool_t    useAutoHandler = kFALSE;
static Int_t     bufferSize     = 3000;
static Int_t     pidArraySize   = 2000;
static Int_t     nMixedEvents   = 10;

static Bool_t    useESDTrackCuts = kTRUE;
static Double_t  cov11 = 2000000.0;              // disabled (was 2.0)
static Double_t  cov22 = 2000000.0;              // disabled (was 2.0)
static Double_t  cov33 = 0.5;
static Double_t  cov44 = 0.5;
static Double_t  cov55 = 2;
static Double_t  nSigmaToVertex = 4;
static Double_t  dcaToVertex = 3.0;
static Double_t  maxChi2PerClusterTPC = 3.5;
static Bool_t    requireTPCRefit = kFALSE;      // disabled (was kTRUE)
static Bool_t    requireSigmaToVertex = kTRUE;
static Bool_t    acceptKinkDaughters = kFALSE;
static Int_t     minNClustersTPC = 50;

Int_t AddRsnAnalysisTask()
{
//
// Core method of the macro.
// Creates and initialized the analysis task object and inserts it
// into the AnalysisManager object passed as argument.
//

  // initialize task objects and retrieve managers for PID and reading
  AliRsnAnalysisSE *task   = new AliRsnAnalysisSE("RsnAnalysisTask", bufferSize);
  AliRsnReader     *reader = task->GetReader();
  AliRsnPID        *pid    = task->GetPID();

  // set the input type
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddRsnAnalysisTask", "NO AnalysisManager found!!!");
    return 0x0;
  }
  task->SetInputType(inputType, mgr, useAutoHandler);

  // define event-mixing cuts
  AliRsnCut    *cutVz     = new AliRsnCut("cutVz" , "", AliRsnCut::kVzDifference, 0.0, 1.0);
  AliRsnCut    *cutMult   = new AliRsnCut("cutMult", "", AliRsnCut::kMultiplicityDifference, 0, 10);
  AliRsnCutSet *cutMixing = new AliRsnCutSet("cutMixing");
  cutMixing->AddCut(cutVz);
  cutMixing->AddCut(cutMult);
  cutMixing->SetCutScheme("cutVz&cutMult");
  task->SetMixingCut(cutMixing);
  task->SetMixingNum(nMixedEvents);

  // ESD settings
  reader->SetCheckSplit(kFALSE);
  reader->SetPIDArraysSize(pidArraySize);
  pid->SetPriorProbability(AliRsnPID::kElectron, 0.02);
  pid->SetPriorProbability(AliRsnPID::kMuon,     0.02);
  pid->SetPriorProbability(AliRsnPID::kPion,     0.83);
  pid->SetPriorProbability(AliRsnPID::kKaon,     0.07);
  pid->SetPriorProbability(AliRsnPID::kProton,   0.06);
  pid->SetMaxPt(1000.0);
  pid->SetMinProb(0.0);

  // ESD primary track cuts
  AliESDtrackCuts *esdCuts = reader->GetESDTrackCuts();
  esdCuts->SetMaxCovDiagonalElements(cov11, cov22, cov33, cov44, cov55);
  esdCuts->SetRequireSigmaToVertex(requireSigmaToVertex);
  if (requireSigmaToVertex) {
    esdCuts->SetMaxNsigmaToVertex(nSigmaToVertex);
  }
  else {
    esdCuts->SetDCAToVertexZ(dcaToVertex);
    esdCuts->SetDCAToVertexXY(dcaToVertex);
  }
  esdCuts->SetRequireTPCRefit(requireTPCRefit);
  esdCuts->SetAcceptKingDaughters(acceptKinkDaughters);
  esdCuts->SetMinNClustersTPC(minNClustersTPC);
  esdCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
  // set usage flag
  reader->SetUseESDTrackCuts(useESDTrackCuts);

  // add configs for PHI and KSTAR in multiplicity bins
  Int_t mult[6] = {0, 25, 50, 75, 100, 200};
  for (Int_t im = 0; im < 5; im++) {
    task->AddPairMgr(CreatePairs(Form("PHI_%d-%d"  , mult[im]+1, mult[im+1]), 333, AliRsnPID::kKaon, AliRsnPID::kKaon, mult[im]+1, mult[im+1]));
    task->AddPairMgr(CreatePairs(Form("KSTAR_%d-%d", mult[im]+1, mult[im+1]), 313, AliRsnPID::kPion, AliRsnPID::kKaon, mult[im]+1, mult[im+1]));
  }
  // integrated
  task->AddPairMgr(CreatePairs("PHI"  , 333, AliRsnPID::kKaon, AliRsnPID::kKaon, 0, 0));
  task->AddPairMgr(CreatePairs("KSTAR", 313, AliRsnPID::kPion, AliRsnPID::kKaon, 0, 0));

  // initialize containers
  AliAnalysisDataContainer *input = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *dummy = mgr->CreateContainer("dummy", TTree::Class(), AliAnalysisManager::kOutputContainer, "default");
  AliAnalysisDataContainer *out   = mgr->CreateContainer("RSN", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

  // connect containers to AnalysisManager
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, input);
  mgr->ConnectOutput(task, 0, dummy);
  mgr->ConnectOutput(task, 1, out);
}

AliRsnPairMgr* CreatePairs
(
  const char       *pairMgrName,    // name for the pair manager
  Int_t             resonancePDG,   // PDG code of resonance (for true pairs)
  AliRsnPID::EType  type1,          // PID of one member of decay (+)
  AliRsnPID::EType  type2,          // PID of other member of decay (-)
  Int_t             multMin,        // lower edge of multiplicity cut
  Int_t             multMax         // upper edge of multiplicity cut
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

  AliRsnPairMgr  *pairMgr  = new AliRsnPairMgr(pairMgrName);

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
  // [2] - mixing
  // [3] - like PP
  // [4] - like MM
  AliRsnPair *noPID[2][5]        = {0,0,0,0,0,0,0,0,0,0};
  AliRsnPair *perfectPID[2][5]   = {0,0,0,0,0,0,0,0,0,0};

  for (i = 0; i < nArray; i++) {
    noPID[i][0] = new AliRsnPair(AliRsnPair::kNoPID, defUnlike[i]);
    noPID[i][1] = new AliRsnPair(AliRsnPair::kNoPID, defUnlike[i]);
    noPID[i][2] = new AliRsnPair(AliRsnPair::kNoPIDMix, defUnlike[i]);
    noPID[i][3] = new AliRsnPair(AliRsnPair::kNoPID, defLikePP[i]);
    noPID[i][4] = new AliRsnPair(AliRsnPair::kNoPID, defLikeMM[i]);

    perfectPID[i][0] = new AliRsnPair(AliRsnPair::kPerfectPID, defUnlike[i]);
    perfectPID[i][1] = new AliRsnPair(AliRsnPair::kPerfectPID, defUnlike[i]);
    perfectPID[i][2] = new AliRsnPair(AliRsnPair::kPerfectPIDMix, defUnlike[i]);
    perfectPID[i][3] = new AliRsnPair(AliRsnPair::kPerfectPID, defLikePP[i]);
    perfectPID[i][4] = new AliRsnPair(AliRsnPair::kPerfectPID, defLikeMM[i]);
  }

  // === CUTS =====================================================================================

  // cuts for tracks:
  // - nothing

  // cuts on pairs:
  // - minimum transverse momentum of the pair (for all)
  // - true daughters of the defined resonance (only for true pairs histogram)
  AliRsnCut    *cutTrue        = new AliRsnCut("cutTrue", "", AliRsnCut::kIsTruePair, resonancePDG);
  AliRsnCutSet *cutSetPairTrue = new AliRsnCutSet("truePairs");
  cutSetPairTrue->AddCut(cutTrue);
  cutSetPairTrue->SetCutScheme("cutTrue");

  // cuts on events:
  // - multiplicity bin
  AliRsnCut    *cutEventMult = new AliRsnCut("cutMult", "", AliRsnCut::kTrueMultiplicity, (Int_t)multMin, (Int_t)multMax);
  AliRsnCutSet *cutSetEvent  = new AliRsnCutSet("multiplicity");
  cutSetEvent->AddCut(cutEventMult);
  cutSetEvent->SetCutScheme("cutMult");

  // define cut manager for NOT true pairs
  AliRsnCutMgr *cutMgr = new AliRsnCutMgr("default", "All pairs");
  if (multMin != multMax) cutMgr->SetCutSet(AliRsnCut::kEvent, cutSetEvent);

  // define cut manager for true pairs
  AliRsnCutMgr *cutMgrTrue = new AliRsnCutMgr("true", "True pairs");
  cutMgrTrue->SetCutSet(AliRsnCut::kPair, cutSetPairTrue);
  if (multMin != multMax) cutMgrTrue->SetCutSet(AliRsnCut::kEvent, cutSetEvent);

  // add cuts to pair analysis
  for (i = 0; i < nArray; i++) {
    noPID[i][0]->SetCutMgr(cutMgrTrue);
    perfectPID[i][0]->SetCutMgr(cutMgrTrue);
    for (j = 1; j < 5; j++) {
      noPID[i][j]->SetCutMgr(cutMgr);
      perfectPID[i][j]->SetCutMgr(cutMgr);
    }
  }

  // === FUNCTIONS ================================================================================

  // define histogram templates
  AliRsnHistoDef *hdIM  = new AliRsnHistoDef(1000, 0.0, 2.0);
  AliRsnFunction *fcnIM = new AliRsnFunction(AliRsnFunction::kInvMass, hdIM);

  // binning in pt
  //Double_t mom[36];
  //mom[0] = 0.2;
  //for (Int_t i = 1; i < 25; i++) mom[i] = mom[i-1] + 0.2;
  //for (Int_t i = 25; i <= 35; i++) mom[i] = mom[i-1] + 0.5;
  //fcnIM->SetBinningCut(AliRsnCut::kTransMomentum, 36, mom, 0);
  //fcnIM->SetBinningCut(AliRsnCut::kTransMomentum, 0.0, 10.0, 0.5, 0);

  // binning in eta
  //fcnIM->SetBinningCut(AliRsnCut::kEta, -1.0, 1.0, 0.2, 1);

  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 5; j++) {
      noPID[i][j]->AddFunction(fcnIM);
      perfectPID[i][j]->AddFunction(fcnIM);
    }
  }

  // === ADD TO PAIR MANAGER ======================================================================

  for (i = 0; i < nArray; i++) {
    for (j = 0; j < 5; j++) {
      pairMgr->AddPair(noPID[i][j]);
      pairMgr->AddPair(perfectPID[i][j]);
    }
  }

  return pairMgr;
}
