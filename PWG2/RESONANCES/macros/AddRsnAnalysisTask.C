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
static Int_t     nMixedEvents   = 5;

static Bool_t    useRsnTrackCuts = kFALSE;
static Double_t  etaRange = 10.0;
static Double_t  rsnImpactParam = 3.0;
static Double_t  ptMin = 0.0;

static Bool_t    useESDTrackCuts = kFALSE;
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

Int_t AddRsnAnalysisTask()
{
//
// Core method of the macro.
// Creates and initialized the analysis task object and inserts it
// into the AnalysisManager object passed as argument.
//

  // initialize task objects and retrieve managers for PID and reading
  AliRsnAnalysisSE *task    = new AliRsnAnalysisSE("RsnAnalysisTask", bufferSize);
  AliRsnReader     *reader  = task->GetReader();
  AliRsnPID        *pid     = task->GetPID();

  // set the input type
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  task->SetInputType(inputType, mgr, useAutoHandler);

  // define event-mixing cuts
  AliRsnCut    *cutVz     = new AliRsnCut("cutVz" , "", AliRsnCut::kVzDifference, 0.0, 0.5);
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
  pid->SetMaxPt(100.0);
  pid->SetMinProb(0.0);

  // default cuts for ITS+TPC
  AliESDtrackCuts *esdCuts = reader->GetESDTrackCuts();
  esdCuts->SetMaxCovDiagonalElements(cov11, cov22, cov33, cov44, cov55);
  esdCuts->SetRequireSigmaToVertex(requireSigmaToVertex);
  if (requireSigmaToVertex) esdCuts->SetMaxNsigmaToVertex(nSigmaToVertex);
  else
  {
    esdCuts->SetDCAToVertexZ(dcaToVertex);
    esdCuts->SetDCAToVertexXY(dcaToVertex);
  }
  esdCuts->SetRequireTPCRefit(requireTPCRefit);
  esdCuts->SetAcceptKingDaughters(acceptKinkDaughters);
  esdCuts->SetMinNClustersTPC(minNClustersTPC);
  esdCuts->SetMaxChi2PerClusterTPC(maxChi2PerClusterTPC);
  reader->SetUseESDTrackCuts(useESDTrackCuts);

  // default cuts from package for tracks
  AliRsnCut *cutPt = new AliRsnCut("cutPt", "", AliRsnCut::kTransMomentumMC, ptMin, 1000000.0);
  AliRsnCut *cutEta = new AliRsnCut("cutEtaMC", "", AliRsnCut::kEtaMC, -etaRange, etaRange);
  AliRsnCut *cutDr = new AliRsnCut("cutDr", "", AliRsnCut::kRadialImpactParam, 0.0, rsnImpactParam);
  AliRsnCutSet *rsnCuts = reader->GetRsnTrackCuts();
  rsnCuts->AddCut(cutPt);
  rsnCuts->AddCut(cutEta);
  rsnCuts->AddCut(cutDr);
  rsnCuts->SetCutScheme("cutPt&cutEta&cutDr");
  reader->SetUseRsnTrackCuts(useRsnTrackCuts);

  // add configs for all impact parameters
  task->AddPairMgrFromConfig("CreatePairsPhi.C");

  // initialize containers
  AliAnalysisDataContainer *input = mgr->CreateContainer("in", TChain::Class(), AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *dummy = mgr->CreateContainer("dummy1", TTree::Class(), AliAnalysisManager::kOutputContainer, "default");
  AliAnalysisDataContainer *out   = mgr->CreateContainer("RSN", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);

  // connect containers to AnalysisManager
  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, input);
  mgr->ConnectOutput(task, 0, dummy);
  mgr->ConnectOutput(task, 1, out);
}
