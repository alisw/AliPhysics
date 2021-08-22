#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "TGeoGlobalMagField.h"

#include "TSystem.h"
#include "TList.h"
#include "TRandom3.h"

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

#include "AliMultDepSpecAnalysisTask.h"
#include "AliMCSpectraWeights.h"

using std::array;
using std::string;
using std::vector;

//**************************************************************************************************
/**
 * Define axis that can be used in the histograms.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTask::SetAxis(unsigned int dim, const std::string name, const std::string title,
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
void AliMultDepSpecAnalysisTask::DefineDefaultAxes(int maxMultMeas, int maxMultTrue)
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

  SetAxis(pt_meas, "pt_meas", "#it{p}^{ meas}_{T} (GeV/#it{c})", ptBins);

  int nBinsMultMeas = maxMultMeas + 1;
  SetAxis(mult_meas, "mult_meas", "#it{N}^{ meas}_{ch}", {-0.5, nBinsMultMeas - 0.5}, nBinsMultMeas);

  if (fIsMC) {
    int nBinsMultTrue = maxMultTrue + 1;
    SetAxis(mult_true, "mult_true", "#it{N}_{ch}", {-0.5, nBinsMultTrue - 0.5}, nBinsMultTrue);
    SetAxis(pt_true, "pt_true", "#it{p}_{T} (GeV/#it{c})", ptBins);
  }
}

//**************************************************************************************************
/**
 * Book a histogram with specified dimensions.
 */
//**************************************************************************************************
template <typename T>
void AliMultDepSpecAnalysisTask::BookHistogram(Hist::Hist<T>& histContainer, const std::string& histName, const std::vector<unsigned int>& dimensions)
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
 * Fill QA histograms with some relevant track properties.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTask::FillTrackQA(AliESDtrack* track)
{
  fHist_signed1Pt.Fill(track->GetSigned1Pt());
  fHist_eta.Fill(track->Eta());
  fHist_phi.Fill(track->Phi());

  fHist_itsFoundClusters.Fill(track->GetITSclusters(0));
  for (unsigned int i = 0; i < 6; i++) {
    if (track->GetITSClusterMap() & (1 << i)) fHist_itsHits.Fill(i + 1);
  }
  fHist_itsChi2PerCluster.Fill((track->GetITSclusters(0) > 0.) ? (track->GetITSchi2() / track->GetITSclusters(0)) : -1.);

  fHist_tpcFindableClusters.Fill(track->GetTPCNclsF());
  fHist_tpcFoundClusters.Fill(track->GetTPCclusters(0));
  fHist_tpcCrossedRows.Fill(track->GetTPCCrossedRows());
  fHist_tpcCrossedRowsOverFindableClusters.Fill((track->GetTPCNclsF() > 0) ? track->GetTPCCrossedRows() / track->GetTPCNclsF() : -1.);
  fHist_tpcChi2PerCluster.Fill((track->GetTPCclusters(0) > 0.) ? track->GetTPCchi2() / track->GetTPCclusters(0) : -1.);
  fHist_tpcGoldenChi2.Fill(track->GetChi2TPCConstrainedVsGlobal(static_cast<AliESDVertex*>(const_cast<AliVVertex*>(fEvent->GetPrimaryVertex()))));
  if (track->GetInnerParam()) fHist_tpcGeomLength.Fill(track->GetLengthInActiveZone(1, 3., 220, fEvent->GetMagneticField()));

  float dca[2];
  float dcaCov[3];
  track->GetImpactParameters(dca, dcaCov);
  fHist_dcaXY.Fill(fPt, dca[0]);
  fHist_dcaZ.Fill(dca[1]);
}

//**************************************************************************************************
/**
 * Function executed once before the event loop. Create histograms here.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTask::BookHistograms()
{
  if (fIsNominalSetting) {
    fQAList = new TList(); // owned by output list
    fQAList->SetOwner();
    fQAList->SetName("QA");
    fOutputList->Add(fQAList);

    fQAList->Add(fHist_trainInfo.GenerateHist("trainInfo"));
    fQAList->Add(fHist_runStatistics.GenerateHist("runStatistics"));

    fHist_eventSelection.AddAxis("selectionStages", "event selection stages", 8, -0.5, 7.5);
    auto evSelHist = fHist_eventSelection.GenerateHist("eventSelection");
    evSelHist->GetXaxis()->SetBinLabel(1, "all");
    evSelHist->GetXaxis()->SetBinLabel(2, "physics selection");
    evSelHist->GetXaxis()->SetBinLabel(3, "no pileup");
    evSelHist->GetXaxis()->SetBinLabel(4, "triggered");
    evSelHist->GetXaxis()->SetBinLabel(5, "has vertex");
    evSelHist->GetXaxis()->SetBinLabel(6, "vertex quality");
    evSelHist->GetXaxis()->SetBinLabel(7, "event quality");
    evSelHist->GetXaxis()->SetBinLabel(8, "vertex position");
    fQAList->Add(evSelHist);

    if (fIsMC) {
      fHist_eventSelectionMC.AddAxis("selectionStages", "event selection stages", 10, -0.5, 9.5);
      auto eventSelectionHistMC = fHist_eventSelectionMC.GenerateHist("eventSelectionMC");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(1, "INEL");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(2, "INEL > 0");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(3, "fiducial");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(4, "physics selection");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(5, "no pileup");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(6, "triggered");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(7, "has vertex");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(8, "vertex quality");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(9, "event quality");
      eventSelectionHistMC->GetXaxis()->SetBinLabel(10, "vertex position");
      fQAList->Add(eventSelectionHistMC);

      fHist_zVtxGen.AddAxis("zv", "V_{z} (cm)", 61, -30.5, 30.5);
      fQAList->Add(fHist_zVtxGen.GenerateHist("zVtxGen"));

      fHist_deltaPt.AddAxis(fAxes[pt_meas]);
      fHist_deltaPt.AddAxis("delta_pt", "#Delta(#it{p}_{T}) / #it{p}^{ meas}_{T}", 100, 0., 0.3);
      fQAList->Add(fHist_deltaPt.GenerateHist("deltaPt"));
    }

    fHist_zVtxMeas.AddAxis("zv", "V^{ meas}_{z} (cm)", 61, -30.5, 30.5);
    fQAList->Add(fHist_zVtxMeas.GenerateHist("zVtxMeas"));

    fHist_sigmaPt.AddAxis(fAxes[pt_meas]);
    fHist_sigmaPt.AddAxis("sigmaPt", "#sigma(#it{p}^{ meas}_{T}) / #it{p}^{ meas}_{T}", 100, 0., 0.3);
    fQAList->Add(fHist_sigmaPt.GenerateHist("sigmaPt"));

    fHist_signed1Pt.AddAxis("signed1Pt", "q/p_{T}", 200, -8, 8);
    fQAList->Add(fHist_signed1Pt.GenerateHist("signed1Pt"));

    fHist_eta.AddAxis("eta", "#eta", 101, -1.0, 1.0);
    fQAList->Add(fHist_eta.GenerateHist("eta"));

    fHist_phi.AddAxis("phi", "#phi [rad]", 90, 0., 2 * M_PI); // 5 bins for each of the 18 sectors
    fQAList->Add(fHist_phi.GenerateHist("phi"));

    fHist_itsFoundClusters.AddAxis("itsFoundClusters", "# clusters ITS", 7, -0.5, 6.5);
    fQAList->Add(fHist_itsFoundClusters.GenerateHist("itsFoundClusters"));

    fHist_itsHits.AddAxis("itsHits", "layer ITS", 6, 0.5, 6.5);
    fQAList->Add(fHist_itsHits.GenerateHist("itsHits"));

    fHist_itsChi2PerCluster.AddAxis("itsChi2PerCluster", "chi2 per cluster ITS", 100, 0, 50);
    fQAList->Add(fHist_itsChi2PerCluster.GenerateHist("itsChi2PerCluster"));

    fHist_tpcFindableClusters.AddAxis("tpcFindableClusters", "# findable clusters TPC", 81, 79.5, 160.5);
    fQAList->Add(fHist_tpcFindableClusters.GenerateHist("tpcFindableClusters"));

    fHist_tpcFoundClusters.AddAxis("foundClustersTPC", "# found clusters TPC", 81, 79.5, 160.5);
    fQAList->Add(fHist_tpcFoundClusters.GenerateHist("tpcFoundClusters"));

    fHist_tpcCrossedRows.AddAxis("tpcCrossedRows", "# crossed rows TPC", 81, 79.5, 160.5);
    fQAList->Add(fHist_tpcCrossedRows.GenerateHist("tpcCrossedRows"));

    fHist_tpcCrossedRowsOverFindableClusters.AddAxis("tpcCrossedRowsOverFindableClusters", "crossed rows / findable clusters TPC", 60, 0.7, 1.3);
    fQAList->Add(fHist_tpcCrossedRowsOverFindableClusters.GenerateHist("tpcCrossedRowsOverFindableClusters"));

    fHist_tpcChi2PerCluster.AddAxis("chi2PerClusterTPC", "chi2 / cluster TPC", 140, 0., 7.);
    fQAList->Add(fHist_tpcChi2PerCluster.GenerateHist("tpcChi2PerCluster"));

    fHist_tpcGoldenChi2.AddAxis("goldenChi2TPC", "chi2 global vs. TPC constrained", 41, -0.5, 40.5);
    fQAList->Add(fHist_tpcGoldenChi2.GenerateHist("tpcGoldenChi2"));

    fHist_tpcGeomLength.AddAxis("geomLengthTPC", "geometric length in TPC", 51, 111.5, 162.5);
    fQAList->Add(fHist_tpcGeomLength.GenerateHist("tpcGeomLength"));

    fHist_dcaXY.AddAxis(fAxes[pt_meas]);
    fHist_dcaXY.AddAxis("dcaXY", "dca XY (cm)", 100, -0.5, 0.5);
    fQAList->Add(fHist_dcaXY.GenerateHist("dcaXY"));

    fHist_dcaZ.AddAxis("dcaZ", "dca Z (cm)", 100, -0.5, 0.5);
    fQAList->Add(fHist_dcaZ.GenerateHist("dcaZ"));
  }

  BookHistogram(fHist_multDist_evt_meas, "multDist_evt_meas", {mult_meas});
  BookHistogram(fHist_multPtSpec_trk_meas, "multPtSpec_trk_meas", {mult_meas, pt_meas});

  if (fIsMC) {
    BookHistogram(fHist_multDist_evt_gen, "multDist_evt_gen", {mult_true});
    BookHistogram(fHist_multDist_evt_gen_trig, "multDist_evt_gen_trig", {mult_true});
    BookHistogram(fHist_multCorrel_evt, "multCorrel_evt", {mult_meas, mult_true});
    BookHistogram(fHist_multCorrel_prim, "multCorrel_prim", {mult_meas, mult_true});
    BookHistogram(fHist_ptCorrel_prim, "ptCorrel_prim", {pt_meas, pt_true});
    BookHistogram(fHist_multPtSpec_prim_gen, "multPtSpec_prim_gen", {mult_true, pt_true});
    BookHistogram(fHist_multPtSpec_prim_gen_sigloss, "multPtSpec_prim_gen_sigloss", {mult_true, pt_true});
    BookHistogram(fHist_multPtSpec_prim_gen_notrig, "multPtSpec_prim_gen_notrig", {mult_true, pt_true});
    BookHistogram(fHist_multPtSpec_prim_meas, "multPtSpec_prim_meas", {mult_true, pt_true});
    BookHistogram(fHist_multPtSpec_trk_prim_meas, "multPtSpec_trk_prim_meas", {mult_meas, pt_meas});
    BookHistogram(fHist_multPtSpec_trk_sec_meas, "multPtSpec_trk_sec_meas", {mult_meas, pt_meas});
  }

  // check required memory
  double requiredMemory =
    fHist_eventSelection.GetSize() +
    fHist_eventSelectionMC.GetSize() +
    fHist_zVtxMeas.GetSize() +
    fHist_multDist_evt_meas.GetSize() +
    fHist_multPtSpec_trk_meas.GetSize() +
    fHist_sigmaPt.GetSize() +
    fHist_zVtxGen.GetSize() +
    fHist_multDist_evt_gen.GetSize() +
    fHist_multDist_evt_gen_trig.GetSize() +
    fHist_deltaPt.GetSize() +
    fHist_multCorrel_evt.GetSize(0.045) +
    fHist_multCorrel_prim.GetSize(0.045) +
    fHist_ptCorrel_prim.GetSize() +
    fHist_multPtSpec_prim_gen.GetSize() +
    fHist_multPtSpec_prim_gen_sigloss.GetSize() +
    fHist_multPtSpec_prim_gen_notrig.GetSize() +
    fHist_multPtSpec_prim_meas.GetSize() +
    fHist_multPtSpec_trk_prim_meas.GetSize() +
    fHist_multPtSpec_trk_sec_meas.GetSize() +
    fHist_signed1Pt.GetSize() +
    fHist_eta.GetSize() +
    fHist_phi.GetSize() +
    fHist_itsFoundClusters.GetSize() +
    fHist_itsHits.GetSize() +
    fHist_itsChi2PerCluster.GetSize() +
    fHist_tpcFindableClusters.GetSize() +
    fHist_tpcFoundClusters.GetSize() +
    fHist_tpcCrossedRows.GetSize() +
    fHist_tpcCrossedRowsOverFindableClusters.GetSize() +
    fHist_tpcChi2PerCluster.GetSize() +
    fHist_tpcGoldenChi2.GetSize() +
    fHist_tpcGeomLength.GetSize();

  // max allowed memory per train job: 8 GiB
  AliError(Form("\n\nEstimated memory usage of histograms: %.2f MiB. For all 23 systematic variations: %.2f MiB\n", requiredMemory / 1048576, 23 * requiredMemory / 1048576));
}

//**************************************************************************************************
/**
 * Function executed once before the event loop.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTask::UserCreateOutputObjects()
{
  OpenFile(1, "recreate");
  fOutputList.reset(new TList());
  fOutputList->SetOwner();
  fRand.reset(new TRandom3());

  // book user histograms
  BookHistograms();

  // configure event selection
  fEventCuts.reset(new AliEventCuts());
  fEventCuts->OverrideCentralityFramework();                   // disables min/max cent cuts in AA as well as cent estimator correlation cuts (which are disabled by default anyhow)
  fEventCuts->OverrideAutomaticTriggerSelection(fTriggerMask); // by default AliEventCuts selects MB and centrality triggered events

  PostData(1, fOutputList.get());
}

//**************************************************************************************************
/**
 * Function executed for each event.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTask::UserExec(Option_t*)
{
  if (!InitEvent()) return;

  if (fIsNominalSetting) {
    // all events -> after physics selection (with pileup rejection) -> no pileup -> triggered -> has vertex -> vertex quality -> event quality -> vertex position
    fHist_eventSelection.Fill(0);
    if (fEventPassesPhysSel) {
      fHist_eventSelection.Fill(1);
      if (fEventCuts->PassedCut(AliEventCuts::kPileUp) && fEventCuts->PassedCut(AliEventCuts::kTPCPileUp)) {
        fHist_eventSelection.Fill(2);
        if (fEventCuts->PassedCut(AliEventCuts::kTrigger)) {
          fHist_eventSelection.Fill(3);
          if (fEventCuts->PassedCut(AliEventCuts::kVertex)) {
            fHist_eventSelection.Fill(4);
            if (fEventCuts->CheckNormalisationMask(AliEventCuts::kHasReconstructedVertex)) {
              fHist_eventSelection.Fill(5);
              if (fEventCuts->CheckNormalisationMask(AliEventCuts::kPassesNonVertexRelatedSelections)) {
                fHist_eventSelection.Fill(6);
                if (fEventCuts->CheckNormalisationMask(AliEventCuts::kPassesAllCuts)) {
                  fHist_eventSelection.Fill(7);
                }
              }
            }
          }
        }
      }
    }

    if (fIsMC) {
      // inel -> inelgt0 -> fiducial -> after physics selection (with pileup rejection) -> no pileup -> triggered -> has vertex -> vertex quality -> event quality -> vertex position
      fHist_eventSelectionMC.Fill(0);
      if (fMCIsINELGT0) {
        fHist_eventSelectionMC.Fill(1);
        if (fMultTrue > 0) {
          fHist_eventSelectionMC.Fill(2);
          if (fEventPassesPhysSel) {
            fHist_eventSelectionMC.Fill(3);
            if (fEventCuts->PassedCut(AliEventCuts::kPileUp) && fEventCuts->PassedCut(AliEventCuts::kTPCPileUp)) {
              fHist_eventSelectionMC.Fill(4);
              if (fEventCuts->PassedCut(AliEventCuts::kTrigger)) {
                fHist_eventSelectionMC.Fill(5);
                if (fEventCuts->PassedCut(AliEventCuts::kVertex)) {
                  fHist_eventSelectionMC.Fill(6);
                  if (fEventCuts->CheckNormalisationMask(AliEventCuts::kHasReconstructedVertex)) {
                    fHist_eventSelectionMC.Fill(7);
                    if (fEventCuts->CheckNormalisationMask(AliEventCuts::kPassesNonVertexRelatedSelections)) {
                      fHist_eventSelectionMC.Fill(8);
                      if (fEventCuts->CheckNormalisationMask(AliEventCuts::kPassesAllCuts)) {
                        fHist_eventSelectionMC.Fill(9);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (fEventCuts->CheckNormalisationMask(AliEventCuts::kHasReconstructedVertex)) {
      fHist_zVtxMeas.Fill(fVtxZ);
    }
  }

  if (fIsMC) {
    if (fIsNominalSetting && fMCIsGoodEventClass) {
      fHist_zVtxGen.Fill(fMCVtxZ);
    }
    if (fMCAcceptEvent) {
      fHist_multDist_evt_gen.Fill(fMultTrue);
      if (fIsTriggered) fHist_multDist_evt_gen_trig.Fill(fMultTrue);
    }
    LoopTrue();
  }

  if (!fAcceptEvent) return;

  if (fIsNominalSetting) {
    fHist_runStatistics.Fill(std::to_string(fRunNumber).data());
    // fill meta train info once per computing job
    if (fIsFirstEventInJob) {
      string trainInfo;
      string period = AliMultSelectionTask::GetPeriodNameByRunNumber(fRunNumber).Data();
      string path{CurrentFileName()};
      if (fIsMC) {
        string prodName = AliMultSelectionTask::GetPeriodNameByGenericPath(path.data()).Data();
        if (!prodName.empty()) trainInfo += prodName + " (" + period + ")";
      } else {
        string recoPass;
        recoPass = path.substr(path.find("pass"));
        recoPass = recoPass.substr(0, recoPass.find("/"));
        recoPass = recoPass.substr(0, recoPass.find("__"));
        recoPass = recoPass.substr(0, recoPass.find("_root_archive_")); // for local testing
        if (!recoPass.empty()) trainInfo += period + "_" + recoPass;
      }
      trainInfo += ", " + fTrainMetadata;
      if (fMCEventClass == EventClass::triggered) {
        trainInfo += ", triggered";
      } else if (fMCEventClass == EventClass::fiducial) {
        trainInfo += ", fiducial";
      } else if (fMCEventClass == EventClass::inelgt0) {
        trainInfo += ", INEL>0";
      } else {
        trainInfo += ", all";
      }
      trainInfo += " events";

      if (!fMCEnableDDC || !fMCSpectraWeights) {
        trainInfo += ", no DDCs";
      }

      fHist_trainInfo.Fill(trainInfo.data());
      fIsFirstEventInJob = false;
    }
  }

  fHist_multDist_evt_meas.Fill(fMultMeas);

  if (fIsMC && fMCAcceptEvent) {
    fHist_multCorrel_evt.Fill(fMultMeas, fMultTrue);
  }
  LoopMeas();

  PostData(1, fOutputList.get());
}

//**************************************************************************************************
/**
 * Function to get data-driven secondary scaling weights.
 */
//**************************************************************************************************
double AliMultDepSpecAnalysisTask::GetSecScalingFactor(AliVParticle* particle)
{
  // FIXME: add PCC handler
  return 1.0;
}

//**************************************************************************************************
/**
 * Function to get data-driven particle composition weights.
 */
//**************************************************************************************************
double AliMultDepSpecAnalysisTask::GetParticleWeight(AliVParticle* particle)
{
  if (fMCSpectraWeights) {
    return fMCSpectraWeights->GetMCSpectraWeight(particle->Particle(), fMCPCCMode);
  }
  return 1.0;
}

//**************************************************************************************************
/**
 * Initialize event quantities. Sets fEvent, fMCEvent, fMeasMult, fTrueMult
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTask::InitEvent()
{
  fEvent = InputEvent();
  if (!fEvent) {
    AliError("fEvent not available\n");
    return false;
  }
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  fEventPassesPhysSel = bool(handl->IsEventSelected());
  fIsTriggered = handl->IsEventSelected() & fTriggerMask;

  if (fIsMC) {
    fMCEvent = MCEvent();
    if (!fMCEvent) {
      AliError("fMCEvent not available\n");
      return false;
    }
    if (fIsNewReco && AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(fMCEvent)) {
      // reject all the rare events with simulated in-bunch pileup
      // those events would have a biased true multiplicity as we cannot identify the particles actually coming from the trigger event
      return false;
    }

    fMCVtxZ = fMCEvent->GetPrimaryVertex()->GetZ();

    if (fMCEnableDDC) {
      // get mc spectra weights object for data-driven corrections
      AliMCSpectraWeightsHandler* mcWeightsHandler = static_cast<AliMCSpectraWeightsHandler*>(fEvent->FindListObject("fMCSpectraWeights"));
      fMCSpectraWeights = (mcWeightsHandler) ? mcWeightsHandler->fMCSpectraWeight : nullptr;
    }
  }

  // event info for random seeds
  fRunNumber = fEvent->GetRunNumber();
  fTimeStamp = fEvent->GetTimeStamp();
  fEventNumber = fEvent->GetEventNumberInFile();

  // check if event is accepted (dataset dependent)
  fAcceptEvent = fEventCuts->AcceptEvent(fEvent);

  fVtxZ = fEvent->GetPrimaryVertex()->GetZ();
  // for the golden chi2 cut, the magnetic field needs to be initialized
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    fEvent->InitMagneticField();
  }

  // first determine mc truth properties
  if (fIsMC) {
    fMCIsINELGT0 = false; // will be set to true in case a charged particle within eta +-1 is found
    LoopTrue(true);       // set true multiplicity fTrueMult
    fMCIsGoodZPos = std::abs(fMCVtxZ) <= 10;

    if (fMCEventClass == EventClass::triggered) {
      fMCIsGoodEventClass = fIsTriggered;
    } else if (fMCEventClass == EventClass::fiducial) {
      fMCIsGoodEventClass = (fMultTrue > 0);
    } else if (fMCEventClass == EventClass::inelgt0) {
      fMCIsGoodEventClass = fMCIsINELGT0;
    } else {
      fMCIsGoodEventClass = true; // fallback to everything the generator produces
    }
    fMCAcceptEvent = fMCIsGoodEventClass && fMCIsGoodZPos;
  }
  LoopMeas(true); // set measured multiplicity fMeasMult

  return true;
}

//**************************************************************************************************
/**
 * Loop over measured tracks. Can either count multiplicity or fill histograms.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTask::LoopMeas(bool count)
{
  if (count) {
    fMultMeas = 0;
  }
  AliVTrack* track = nullptr;
  for (int i = 0; i < fEvent->GetNumberOfTracks(); ++i) {
    track = dynamic_cast<AliVTrack*>(fEvent->GetTrack(i));
    // set track properties and check if track is in kin range and has good quality
    if (!InitTrack(track)) continue;

    // initialize particle properties
    bool isValidParticle = false;
    if (fIsMC) {
      // mc label corresponding to measured track (negative label indicates bad quality track)
      int mcLabel = std::abs(track->GetLabel());

      // set mc particle properties and check if it is charged prim/sec and in kin range
      if (fIsESD) {
        isValidParticle = InitParticle((AliMCParticle*)fMCEvent->GetTrack(mcLabel));
      } else {
        isValidParticle = InitParticle((AliAODMCParticle*)fMCEvent->GetTrack(mcLabel));
      }
      if (fMCIsPileupParticle) continue; // skip tracks from pileup in mc
    }

    if (count) {
      fMultMeas += fNRepetitions;
    } else {

      if (fIsESD && fIsNominalSetting) FillTrackQA(dynamic_cast<AliESDtrack*>(track));

      for (int i = 0; i < fNRepetitions; ++i) {

        fHist_multPtSpec_trk_meas.Fill(fMultMeas, fPt);
        if (fIsNominalSetting) fHist_sigmaPt.Fill(fPt, fSigmaPt);

        if (fIsMC && fMCAcceptEvent && isValidParticle) {

          if (fIsNominalSetting) fHist_deltaPt.Fill(fPt, TMath::Abs(fPt - fMCPt) / fPt);

          if (fMCIsChargedPrimary) {
            fHist_multPtSpec_trk_prim_meas.Fill(fMultMeas, fPt);

            fHist_multPtSpec_prim_meas.Fill(fMultTrue, fMCPt);

            fHist_ptCorrel_prim.Fill(fPt, fMCPt);
            fHist_multCorrel_prim.Fill(fMultMeas, fMultTrue);

          } else {
            fHist_multPtSpec_trk_sec_meas.Fill(fMultMeas, fPt);
          }
        }
      }
    }
  }
}

//**************************************************************************************************
/**
 * Loop over generated mc particles. Can either count multiplicity or fill histograms.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTask::LoopTrue(bool count)
{
  if (count) {
    fMultTrue = 0;
  }

  for (int i = 0; i < fMCEvent->GetNumberOfTracks(); ++i) {
    // sets fMCPt, fMCEta, ... and checks if particle in kin range
    if (fIsESD) {
      if (!InitParticle((AliMCParticle*)fMCEvent->GetTrack(i))) continue;
    } else {
      if (!InitParticle((AliAODMCParticle*)fMCEvent->GetTrack(i))) continue;
    }
    if (!fMCIsChargedPrimary) continue;

    if (count) {
      fMultTrue += fNRepetitions;
    } else {
      for (int i = 0; i < fNRepetitions; ++i) {
        if (fMCAcceptEvent) {
          fHist_multPtSpec_prim_gen.Fill(fMultTrue, fMCPt);
          if (!fAcceptEvent) {
            fHist_multPtSpec_prim_gen_sigloss.Fill(fMultTrue, fMCPt);
          }
        }
        if (fMCIsGoodZPos && !fIsTriggered) {
          fHist_multPtSpec_prim_gen_notrig.Fill(fMultTrue, fMCPt);
        }
      }
    }
  }
}

//**************************************************************************************************
/**
 * Initializes track properties and returns false if track is not available, has bad quality or is not in kinematic range.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTask::InitTrack(AliVTrack* track)
{
  if (!track) {
    AliFatal("Track not found\n");
    return false;
  }
  // reset the two mc variables that have side effects on measured loop in mc
  fNRepetitions = 1;
  fMCIsPileupParticle = false;

  // temporary solution to remove tracks from background events (corresponds to fMCIsPileupParticle = true)
  // FIXME: we may actually not want to reject those tracks but count them as contamination instead
  if (fIsMC && fIsNewReco && std::abs(track->GetLabel()) >= AliMCEvent::BgLabelOffset()) return false;

  fPt = track->Pt();
  fEta = track->Eta();
  if ((fPt <= fMinPt + PRECISION) || (fPt >= fMaxPt - PRECISION) || (fEta <= fMinEta + PRECISION) || (fEta >= fMaxEta - PRECISION)) {
    return false;
  }
  if (!fTrackCuts->IsSelected(track)) return false;

  if (fIsESD) {
    fSigmaPt = fPt * TMath::Sqrt(dynamic_cast<AliESDtrack*>(track)->GetSigma1Pt2());
  } else {
    // for AODs this is only possible with massive overhead
    // (cov matrix entries defined in pxpypz space need to be converted back to sigma 1/pt)
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(track);
    double cov[21] = {0};
    double pxpypz[3] = {0};
    double xyz[3] = {0};
    AliExternalTrackParam exParam;
    aodTrack->GetCovMatrix(cov);
    aodTrack->PxPyPz(pxpypz);
    aodTrack->GetXYZ(xyz);
    exParam.Set(xyz, pxpypz, cov, aodTrack->Charge());
    fSigmaPt = fPt * TMath::Sqrt(exParam.GetSigma1Pt2());
  }
  return true;
}

//**************************************************************************************************
/**
 * Initializes particle properties and returns false if not charged primary or secondary or if
 * particle is of kinematic range range. Works for AliMCParticles (ESD) and AliAODMCParticles (AOD).
 */
//**************************************************************************************************
template <typename Particle_t>
bool AliMultDepSpecAnalysisTask::InitParticle(Particle_t* particle)
{
  if (!particle) {
    AliFatal("Particle not found\n");
    return false;
  }
  fMCLabel = particle->GetLabel();
  fMCIsPileupParticle = false;
  // reject all particles and tracks that come from simulated out-of-bunch pileup
  if (fIsNewReco && AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(fMCLabel, fMCEvent)) {
    fMCIsPileupParticle = true; // store this info as it is relevant for track loop as well
    return false;
  }

  if (!(TMath::Abs(particle->Charge()) > 0.01)) return false; // reject all neutral particles

  fMCIsChargedPrimary = fMCEvent->IsPhysicalPrimary(fMCLabel);
  fMCIsChargedSecondary = (fMCIsChargedPrimary) ? false : (fMCEvent->IsSecondaryFromWeakDecay(fMCLabel) || fMCEvent->IsSecondaryFromMaterial(fMCLabel));

  // not interested in anything non-final
  if (!(fMCIsChargedPrimary || fMCIsChargedSecondary)) return false;

  fMCPt = particle->Pt();
  fMCEta = particle->Eta();

  // event class is INEL>0 in case it has a charged particle in abs(eta) < 1
  fMCIsINELGT0 = fMCIsINELGT0 || (fMCIsChargedPrimary && (std::abs(fMCEta) < 1.));

  fMCParticleWeight = 1.0;
  fMCSecScaleWeight = 1.0;
  fNRepetitions = 1;

  if (fMCEnableDDC) {
    if (fMCIsChargedPrimary) {
      fMCParticleWeight = GetParticleWeight(static_cast<AliVParticle*>(particle));
      fNRepetitions = GetNRepetitons(fMCParticleWeight);
    } else {
      fMCSecScaleWeight = GetSecScalingFactor(static_cast<AliVParticle*>(particle));
      fNRepetitions = GetNRepetitons(fMCSecScaleWeight);
    }
  }

  if ((fMCPt <= fMinPt + PRECISION) || (fMCPt >= fMaxPt - PRECISION) || (fMCEta <= fMinEta + PRECISION) || (fMCEta >= fMaxEta - PRECISION)) {
    return false;
  }
  return true;
}

//**************************************************************************************************
/**
 * Decide how often to repeat particle in MC to match data.
 */
//**************************************************************************************************
int AliMultDepSpecAnalysisTask::GetNRepetitons(double scalingFactor)
{
  int nRepetitions = (int)scalingFactor;
  double rest = scalingFactor - nRepetitions;

  fRand->SetSeed(GetSeed());
  nRepetitions += (fRand->Rndm() <= rest) ? 1 : 0;
  return nRepetitions;
}

//**************************************************************************************************
/**
 * Define random (but reproducable) seed.
 */
//**************************************************************************************************
unsigned long AliMultDepSpecAnalysisTask::GetSeed()
{
  if (fUseRandomSeed) {
    return 0;
  }

  unsigned long seed = fEventNumber;
  seed <<= 5;
  seed += fRunNumber;
  seed <<= 2;
  seed += fMCLabel;
  seed <<= 3;
  seed += fTimeStamp;
  return seed;
}

//**************************************************************************************************
/**
 * Function to initialise centrality.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTask::InitCentrality()
{
  AliMultSelection* multSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if (multSelection) {
    fCent = multSelection->GetMultiplicityPercentile("V0M");
  } else {
    // check for legacy cent selection used in old PbPb data
    AliCentrality* centSelection = fEvent->GetCentrality();
    if (centSelection) {
      fCent = centSelection->GetCentralityPercentile("V0M");
    } else {
      AliError("Centrality wagon not found!");
      fCent = 999.;
      return false;
    }
  }
  return true;
}

//**************************************************************************************************
/**
 * Function to add this task to a train.
 */
//**************************************************************************************************
AliMultDepSpecAnalysisTask* AliMultDepSpecAnalysisTask::AddTaskMultDepSpec(const string& dataSet, int cutModeLow, int cutModeHigh, TString options, bool isMC)
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
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  bool isAOD = false;
  if (type.Contains("AOD")) {
    isAOD = true;
  } else {
    // for ESDs isMC can be determined automatically
    isMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler() != nullptr);
  }
  string mode = (isMC) ? "MC" : "Data";

  AliMultDepSpecAnalysisTask* returnTask = nullptr;
  char taskName[100] = "";
  for (int cutMode = cutModeLow; cutMode <= cutModeHigh; cutMode++) {
    sprintf(taskName, "%s_%s_cutMode_%d", dataSet.data(), mode.data(), cutMode);
    AliMultDepSpecAnalysisTask* task = new AliMultDepSpecAnalysisTask(taskName);
    if (!task->InitTask(isMC, isAOD, dataSet, options, cutMode)) {
      delete task;
      task = nullptr;
      continue;
    }
    if (!returnTask) returnTask = task; // return one of the tasks
    task->SaveTrainMetadata();

    // hang task in train
    mgr->AddTask(task);
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer(taskName, TList::Class(), AliAnalysisManager::kOutputContainer, "AnalysisResults.root"));
  }
  return returnTask;
}

//**************************************************************************************************
/**
 * Function to extract metadata string for current train run.
 */
//**************************************************************************************************
void AliMultDepSpecAnalysisTask::SaveTrainMetadata()
{
  // Save train metadata
  string trainID = (gSystem->Getenv("TRAIN_ID")) ? gSystem->Getenv("TRAIN_ID") : "";
  string trainRun = (gSystem->Getenv("TRAIN_RUN_ID")) ? gSystem->Getenv("TRAIN_RUN_ID") : "";
  string aliPhysTag = (gSystem->Getenv("ALIROOT_VERSION")) ? gSystem->Getenv("ALIROOT_VERSION") : "";
  if (aliPhysTag.find("::") != string::npos) aliPhysTag = aliPhysTag.substr(aliPhysTag.find("::") + 2);
  std::map<string, string> trainIdNames = {
    // ESD trains
    {"36", "LF_pp"},
    {"37", "LF_pp_MC"},
    {"51", "LF_pPb"},
    {"53", "LF_pPb_MC"},
    {"26", "LF_PbPb"},
    {"27", "LF_PbPb_MC"},
    // AOD trains
    {"39", "LF_pp_AOD"},
    {"38", "LF_pp_MC_AOD"},
    {"72", "LF_pPb_AOD"},
    {"73", "LF_pPb_MC_AOD"},
    {"20", "LF_PbPb_AOD"},
    {"21", "LF_PbPb_MC_AOD"},
  };
  string trainName = (trainIdNames.find(trainID) == trainIdNames.end()) ? "-" : trainIdNames[trainID];
  fTrainMetadata = trainName + "#" + trainRun + " @ " + aliPhysTag;
}

//**************************************************************************************************
/**
 * Initialize task for specific cut mode. Creates ESD track cuts and particle compositon objects with correct settings.
 * Call this funktion in the AddTask function. The resulting settings will be streamed.
 * In AODs the track selections are based on the best possible conversion of the specified esd cuts (see AliESDTrackCuts::AcceptVTrack()).
 * Due to the nature of the AOD format it is not possible to have an exact 1-to-1conversion for all of the cuts (e.g. geom length cut).
 * In addition, some cut variables are not available in early AOD productions: e.g. the golden chi2 cut can only be applied since August 2016.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTask::InitTask(bool isMC, bool isAOD, string dataSet, TString options, int cutMode)
{
  if ((!isMC && cutMode > 119) || cutMode < 100 || cutMode > 122) return false;
  if (cutMode == 100) fIsNominalSetting = true;
  fIsMC = isMC;
  fIsESD = !isAOD;
  if (!SetupTask(dataSet, options)) return false;
  string colSys = dataSet.substr(0, dataSet.find("_"));

  fTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

  // these are the default track selections (cutMode == 100)
  fTrackCuts->SetRequireTPCRefit(true);
  fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  if (fIsNewReco) {
    // all reconstructions since mid 2020 require a tighter cut
    fTrackCuts->SetMaxChi2PerClusterTPC(2.5);
  } else {
    fTrackCuts->SetMaxChi2PerClusterTPC(4.0);
  }
  fTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
  fTrackCuts->SetRequireITSRefit(true);
  fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fTrackCuts->SetMaxChi2PerClusterITS(36.);
  fTrackCuts->SetDCAToVertex2D(false);
  fTrackCuts->SetRequireSigmaToVertex(false);
  fTrackCuts->SetMaxDCAToVertexZ(2.0);
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); // 7 sigma cut
  fTrackCuts->SetAcceptKinkDaughters(false);
  fTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);     // golden chi2 cut
  fTrackCuts->SetCutGeoNcrNcl(3, 130, 1.5, 0.85, 0.7); // geometrical length cut

  // cut-variations:
  if (cutMode == 101) {
    fTrackCuts->SetMaxChi2PerClusterITS(25.);
  } else if (cutMode == 102) {
    fTrackCuts->SetMaxChi2PerClusterITS(49.);
  } else if (cutMode == 103) {
    if (fIsNewReco) {
      fTrackCuts->SetMaxChi2PerClusterTPC(1.7);
    } else {
      fTrackCuts->SetMaxChi2PerClusterTPC(3.0);
    }
  } else if (cutMode == 104) {
    if (fIsNewReco) {
      fTrackCuts->SetMaxChi2PerClusterTPC(3.4);
    } else {
      fTrackCuts->SetMaxChi2PerClusterTPC(5.0);
    }
  } else if (cutMode == 105) {
    fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
  } else if (cutMode == 106) {
    fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
  } else if (cutMode == 107) {
    fTrackCuts->SetMaxFractionSharedTPCClusters(0.2);
  } else if (cutMode == 108) {
    fTrackCuts->SetMaxFractionSharedTPCClusters(1.0);
  } else if (cutMode == 109) {
    fTrackCuts->SetMaxChi2TPCConstrainedGlobal(25.);
  } else if (cutMode == 110) {
    fTrackCuts->SetMaxChi2TPCConstrainedGlobal(49.);
  } else if (cutMode == 111) {
    fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0104+0.0200/pt^1.01");
  } else if (cutMode == 112) {
    fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0260+0.0500/pt^1.01");
  } else if (cutMode == 113) {
    fTrackCuts->SetMaxDCAToVertexZ(1.0);
  } else if (cutMode == 114) {
    fTrackCuts->SetMaxDCAToVertexZ(5.0);
  } else if (cutMode == 115) {
    fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
  } else if (cutMode == 116) {
    fTrackCuts->SetCutGeoNcrNcl(3, 120, 1.5, 0.85, 0.7);
  } else if (cutMode == 117) {
    fTrackCuts->SetCutGeoNcrNcl(3, 140, 1.5, 0.85, 0.7);
  } else if (cutMode == 118) {
    fTrackCuts->SetCutGeoNcrNcl(4, 130, 1.5, 0.85, 0.7);
  } else if (cutMode == 119) {
    fTrackCuts->SetCutGeoNcrNcl(2, 130, 1.5, 0.85, 0.7);
  }

  fMCEnableDDC = true;
  if (options.Contains("noDDC")) {
    fMCEnableDDC = false;
  }

  // in MC we always apply data driven corrections to account for wrong particle composition in the generator
  if (isMC && fMCEnableDDC) {
    fMCPCCMode = 0;
    fMCSecScalingMode = 0;

    if (cutMode == 120) {
      fMCPCCMode = 1; // shift correction factor up within its systematics
    } else if (cutMode == 121) {
      fMCPCCMode = -1; // shift correction factor down within its systematics
    } else if (cutMode == 122) {
      fMCSecScalingMode = 1; // use 3 template dca fits as variation
    }
  }
  return true;
}

//**************************************************************************************************
/**
 * Apply default and data set specific settings like triggers and mulitplicity binning. Override defaults with user options.
 */
//**************************************************************************************************
bool AliMultDepSpecAnalysisTask::SetupTask(string dataSet, TString options)
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

  // dataset specific trigger settings
  unsigned int triggerMask = AliVEvent::kINT7; // V0and used in Run2 (non single diffractive!)
  if (dataSet.find("pp_2TeV") != string::npos || dataSet.find("pp_7TeV") != string::npos || dataSet.find("PbPb_2TeV") != string::npos) {
    triggerMask = AliVEvent::kMB; // V0or used in Run1
  }
  SetTriggerMask(triggerMask);

  unsigned int eventCass = EventClass::fiducial; // by default normalize to fiducial event class
  if (options.Contains("EventClass::triggered")) {
    eventCass = EventClass::triggered;
  } else if (options.Contains("EventClass::inelgt0")) {
    eventCass = EventClass::inelgt0;
  }
  SetEventClass(eventCass);

  // for reconstructions after mid 2020 chi2/clusterTPC cut moves from 4 -> 2.5 and MCs contain out-of-bunch pileup
  if (!options.Contains("oldReco") &&
      (dataSet.find("pp_13TeV") != string::npos || dataSet.find("pPb") != string::npos || dataSet.find("PbPb_5TeV") != string::npos)) {
    SetIsNewReco();
  }

  int maxMultMeas = 100;
  int maxMultTrue = 100;

  if (dataSet.find("pPb") != string::npos) {
    maxMultMeas = 200;
    maxMultTrue = 200;
  } else if (dataSet.find("PbPb") != string::npos || dataSet.find("XeXe") != string::npos) {
    maxMultMeas = 2500;
    maxMultTrue = 3800;
  }

  if (!options.Contains("consistentBinningAA")) {
    if (dataSet.find("PbPb_2TeV") != string::npos) {
      maxMultMeas = 2000;
      maxMultTrue = 3200;
    }
    if (dataSet.find("XeXe") != string::npos) {
      maxMultMeas = 1800;
      maxMultTrue = 2800;
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

  DefineDefaultAxes(maxMultMeas, maxMultTrue);
  return true;
}
