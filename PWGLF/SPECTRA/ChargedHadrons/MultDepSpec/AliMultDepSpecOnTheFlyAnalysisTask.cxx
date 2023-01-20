#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"

#include "TSystem.h"
#include "TList.h"

#include "AliVEventHandler.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliMultSelectionTask.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliEventCuts.h"

#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliVHeader.h"

#include "AliMultDepSpecOnTheFlyAnalysisTask.h"

using std::array;
using std::string;
using std::vector;

//**************************************************************************************************
/**
 * Define axis that can be used in the histograms.
 */
//**************************************************************************************************
void AliMultDepSpecOnTheFlyAnalysisTask::SetAxis(unsigned int dim, const std::string name, const std::string title,
                                                 const std::vector<double>& binEdges, int nBins)
{
  if (binEdges.size() != 2 && nBins != 0) {
    AliError("Specifying the number of bins is only required for fixed bin widths.");
    nBins = 0;
  }
  fAxes[dim] = {name, title, binEdges, nBins};
}

//**************************************************************************************************
/**
 * Define default axis properties.
 */
//**************************************************************************************************
void AliMultDepSpecOnTheFlyAnalysisTask::DefineDefaultAxes()
{
  std::vector<double> ptBins = {0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                                0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5,
                                6.0, 6.5, 7.0, 8.0, 9.0, 10.0};
  if (fHighPtMode == 1) {
    // simple extension of pt range to 50 GeV/c
    std::vector<double> highPtBins = {20., 30., 40., 50.};
    ptBins.insert(ptBins.end(), highPtBins.begin(), highPtBins.end());
  } else if (fHighPtMode == 2) {
    // simple extension of pt range to 50 GeV/c with less bins
    std::vector<double> highPtBins = {20., 50.};
    ptBins.insert(ptBins.end(), highPtBins.begin(), highPtBins.end());
  } else if (fHighPtMode == 3) {
    // binning for improved RAA reference up to 100 GeV/c
    std::vector<double> highPtBins = {11., 12., 13., 14., 15., 16., 18., 20., 22., 24., 26., 30., 34., 40., 50., 60., 80., 100.};
    ptBins.insert(ptBins.end(), highPtBins.begin(), highPtBins.end());
  }

  int nBinsMultTrue = fMaxMultTrue + 1;
  SetAxis(mult_true, "mult_true", "#it{N}_{ch}", {-0.5, nBinsMultTrue - 0.5}, nBinsMultTrue);
  SetAxis(pt_true, "pt_true", "#it{p}_{T} (GeV/#it{c})", ptBins);
}

//**************************************************************************************************
/**
 * Book a histogram with specified dimensions.
 */
//**************************************************************************************************
template <typename T>
void AliMultDepSpecOnTheFlyAnalysisTask::BookHistogram(Hist<T>& histContainer, const std::string& histName, const std::vector<unsigned int>& dimensions)
{
  for (auto& dim : dimensions) {
    if (fAxes.find(dim) == fAxes.end()) {
      AliFatal(Form("Not all axes for histogram %s were specified properly!", histName.data()));
      return;
    }
    histContainer.AddAxis(fAxes[dim]);
  }
  fOutputList->Add(histContainer.GenerateHist(histName));
}

//**************************************************************************************************
/**
 * Function executed once before the event loop. Create histograms here.
 */
//**************************************************************************************************
void AliMultDepSpecOnTheFlyAnalysisTask::BookHistograms()
{
  BookHistogram(fHist_multDist_evt_gen, "multDist_evt_gen", {mult_true});
  BookHistogram(fHist_multPtSpec_prim_gen, "multPtSpec_prim_gen", {mult_true, pt_true});
}

//**************************************************************************************************
/**
 * Function executed once before the event loop.
 */
//**************************************************************************************************
void AliMultDepSpecOnTheFlyAnalysisTask::UserCreateOutputObjects()
{
  OpenFile(1, "recreate");
  fOutputList.reset(new TList());
  fOutputList->SetOwner();
  DefineDefaultAxes();
  BookHistograms();
  PostData(1, fOutputList.get());
}

//**************************************************************************************************
/**
 * Function executed for each event.
 */
//**************************************************************************************************
void AliMultDepSpecOnTheFlyAnalysisTask::UserExec(Option_t*)
{
  if (!InitEvent()) return;
  if (fMCAcceptEvent) {
    fHist_multDist_evt_gen.Fill(fMultTrue);
  }
  LoopTrue();
  PostData(1, fOutputList.get());
}

//**************************************************************************************************
/**
 * Initialize event quantities. Sets fEvent, fMCEvent, fMeasMult, fTrueMult
 */
//**************************************************************************************************
bool AliMultDepSpecOnTheFlyAnalysisTask::InitEvent()
{
  fMCEvent = MCEvent();
  if (!fMCEvent) {
    AliError("fMCEvent not available\n");
    return false;
  }
  fMCVtxZ = fMCEvent->GetPrimaryVertex()->GetZ();
  LoopTrue(true); // set true multiplicity fTrueMult
  fMCIsGoodZPos = std::abs(fMCVtxZ) <= 10;
  fMCIsGoodEventClass = (fMultTrue > 0);
  fMCAcceptEvent = fMCIsGoodEventClass && fMCIsGoodZPos;
  return true;
}

//**************************************************************************************************
/**
 * Loop over generated mc particles. Can either count multiplicity or fill histograms.
 */
//**************************************************************************************************
void AliMultDepSpecOnTheFlyAnalysisTask::LoopTrue(bool count)
{
  if (count) {
    fMultTrue = 0;
  }
  for (int particleID = 0; particleID < fMCEvent->GetNumberOfTracks(); ++particleID) {
    if (!InitParticle<AliMCParticle>(particleID)) continue;
    if (!fMCIsChargedPrimary) continue;
    
    if (count) {
      ++fMultTrue;
    } else {
      if (fMCAcceptEvent) {
        fHist_multPtSpec_prim_gen.Fill(fMultTrue, fMCPt);
      }
    }
  }
}

//**************************************************************************************************
/**
 * Initializes particle properties and returns false if not charged primary or secondary or if
 * particle is of kinematic range range. Works for AliMCParticle (ESD) and AliAODMCParticle (AOD).
 */
//**************************************************************************************************
template <typename Particle_t>
bool AliMultDepSpecOnTheFlyAnalysisTask::InitParticle(int particleID)
{
  Particle_t* particle = static_cast<Particle_t*>(fMCEvent->GetTrack(particleID));

  if (!particle) {
    AliFatal("Particle not found\n");
    return false;
  }

  if (!(TMath::Abs(particle->Charge()) > 0.01)) return false; // reject all neutral particles

  fMCIsChargedPrimary = fMCEvent->IsPhysicalPrimary(particleID);
  bool isChargedSecondaryFromWeakDecay = (fMCIsChargedPrimary) ? false : fMCEvent->IsSecondaryFromWeakDecay(particleID);
  fMCIsChargedSecondary = (fMCIsChargedPrimary) ? false : (isChargedSecondaryFromWeakDecay || fMCEvent->IsSecondaryFromMaterial(particleID));

  // not interested in anything non-final
  if (!(fMCIsChargedPrimary || fMCIsChargedSecondary)) return false;

  fMCPt = particle->Pt();
  fMCEta = particle->Eta();

  if ((fMCPt <= fMinPt + PRECISION) || (fMCPt >= fMaxPt - PRECISION) || (fMCEta <= fMinEta + PRECISION) || (fMCEta >= fMaxEta - PRECISION)) {
    return false;
  }
  return true;
}

//**************************************************************************************************
/**
 * Function to add this task to a train.
 */
//**************************************************************************************************
AliMultDepSpecOnTheFlyAnalysisTask* AliMultDepSpecOnTheFlyAnalysisTask::AddTaskMultDepSpecOnTheFlyMC(const string& dataSet, TString options)
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMultDepSpec", "No analysis manager found.");
    return nullptr;
  }
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMultDepSpec", "No input event handler found.");
    return nullptr;
  }

  char taskName[100] = "";
  sprintf(taskName, "%s_onTheFlyMC", dataSet.data());
  AliMultDepSpecOnTheFlyAnalysisTask* task = new AliMultDepSpecOnTheFlyAnalysisTask(taskName);
  task->InitTask(dataSet, options);

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(taskName, TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));
  return task;
}

//**************************************************************************************************
/**
 Initialise task.
 */
//**************************************************************************************************
bool AliMultDepSpecOnTheFlyAnalysisTask::InitTask(string dataSet, TString options)
{
  if (!SetupTask(dataSet, options)) return false;
  string colSys = dataSet.substr(0, dataSet.find("_"));
  return true;
}

//**************************************************************************************************
/**
 * Apply default and data set specific settings like triggers and mulitplicity binning. Override defaults with user options.
 */
//**************************************************************************************************
bool AliMultDepSpecOnTheFlyAnalysisTask::SetupTask(string dataSet, TString options)
{
  vector<string> dataSets = {
    "pp_2TeV",
    "pp_5TeV",
    "pp_7TeV",
    "pp_8TeV",
    "pp_13TeV",
    "pPb_5TeV",
    "pPb_8TeV",
    "XeXe_5TeV",
    "PbPb_2TeV",
    "PbPb_5TeV",
  };

  if (std::find(dataSets.begin(), dataSets.end(), dataSet) == dataSets.end()) {
    AliError("Settings for specified dataset are not defined!\n");
    return false;
  }

  if (dataSet.find("pPb") != string::npos) {
    fMaxMultTrue = 200;
  } else if (dataSet.find("PbPb") != string::npos || dataSet.find("XeXe") != string::npos) {
    fMaxMultTrue = 3800;
  }
  if (!options.Contains("consistentBinningAA")) {
    if (dataSet.find("PbPb_2TeV") != string::npos) {
      fMaxMultTrue = 3200;
    }
    if (dataSet.find("XeXe") != string::npos) {
      fMaxMultTrue = 3200;
    }
  }

  // kinematic cuts:
  SetEtaRange(-0.8, 0.8);
  SetPtRange(0.15, 10.0);

  if (options.Contains("highPtMode::50")) {
    fHighPtMode = 1;
    SetPtRange(0.15, 50.0);
  }
  if (options.Contains("highPtMode::50coarse")) {
    fHighPtMode = 2;
    SetPtRange(0.15, 50.0);
  }
  if (options.Contains("highPtMode::100")) {
    fHighPtMode = 3;
    SetPtRange(0.15, 100.0);
  }

  return true;
}
