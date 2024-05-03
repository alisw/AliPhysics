/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.
 *                                                                        *
 * Authors: Omar Vazquez (omar.vazquez.rueda@cern.ch)                     *
 **************************************************************************/

/* This source code yields the histograms for the measurement of the
 * speed of sound using very central Pb-Pb collisions.
 */

class TTree;

class AliPPVsMultUtils;
class AliESDtrackCuts;

#include <AliAnalysisFilter.h>
#include <AliESDVertex.h>
#include <AliHeader.h>
#include <AliMultiplicity.h>
#include <TBits.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>

#include <iostream>
#include <vector>

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDUtils.h"
#include "AliESDVZERO.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
// #include "AliEventCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
// #include "AliMCParticle.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TParticle.h"
#include "TProfile.h"
#include "TVector3.h"

using std::cout;
using std::endl;
using std::vector;

static constexpr int v0m_Nbins{1};
// static constexpr double v0m_bins[v0m_Nbins + 1] = {0.0, 5.0};
static constexpr double uc_v0m_bins_high[v0m_Nbins] = {5.0};
static constexpr double uc_v0m_bins_low[v0m_Nbins] = {0.0};
static const char* uc_v0m_bins_name[v0m_Nbins] = {"0_5"};

#include "AliAnalysisTaskDataSpeedOfSound.h"

class AliAnalysisTaskDataSpeedOfSound;  // your analysis class

ClassImp(AliAnalysisTaskDataSpeedOfSound)  // classimp: necessary for root

    AliAnalysisTaskDataSpeedOfSound::AliAnalysisTaskDataSpeedOfSound()
    : AliAnalysisTaskSE(),
      fESD(0),
      fEventCuts(0x0),
      fMCStack(0),
      fMC(0),
      fUseMC(kFALSE),
      fIsSystematics(true),
      fVaryVtxZPos(false),
      fMinVtxZPos(-5.0),
      fMaxVtxZPos(5.0),
      fSystematic(1),
      fTrigger(AliVEvent::kCentral),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fOutputList(0),
      fEtaCut(0.8),
      fEtaCutForpTwSPDGap(0.4),
      fEtaCutSPDGapMin(0.7),
      fEtaCutSPDGapMax(1.4),
      fEtaCutForpTwTPCGap(0.3),
      fEtaCutTPCGapMin(0.5),
      fEtaCutTPCGapMax(0.8),
      fPtMin(0.15),
      fPtMinCent(0.15),
      fPtMaxCent(10.0),
      fV0Mmin(0.0),
      fV0Mmax(80.0),
      fHMCut(10.0),
      ftrackmult08(0),
      fv0mpercentile(0),
      fv0mamplitude(0),
      fTrackletsEtaGap(0),
      fTracksEtaGapTPC(0),
      fMultSelection(0x0),
      hNch(0),
      hNchvsV0MAmp(0),
      hV0MvsV0MAmp(0),
      pV0MAmpChannel(0),
      hV0MAmplitude(0),
      hV0Percentile(0),
      hPtvsNch(0),
      pPtvsNch(0),
      pPtvsV0MAmp(0),
      hPtvsV0MAmp(0),
      hPhiEtaGapTPC(0),
      hPtWithCutForCent(0),
      hPhiEtaGapSPD(0),
      hBestVtxZ(0),
      hTrackletsEtaGap(0),
      hTracksEtaGapTPC(0),
      hPtvsTrackletsEtaGap(0),
      hPtvsTracksEtaGapTPC(0),
      pPtvsTrackletsEtaGap(0),
      pPtvsTracksEtaGapTPC(0),
      hSPDFull(0),
      hSPDEtaAdj(0),
      hSPDEtaGapW(0),
      hSPDEtaGapWW(0),
      hEtFull(0),
      hEtEtaGap(0),
      hPtvsSPDFull(0),
      hPtvsSPDEtaAdj(0),
      hPtvsSPDEtaGapW(0),
      hPtvsSPDEtaGapWW(0),
      hPtvsEtFull(0),
      hPtvsEtEtaGap(0),
      hPtvsTPCEtaGapWidepT(0),
      hPtvsEtEtaGapWidepT(0),
      pZVtxvsSPDClus(0),
      pSPDClusvsEta(0),
      fSPDVtxCut(3.0),
      fSPDFull(0),
      fSPDEtaAdj(0),
      fSPDEtaGapW(0),
      fSPDEtaGapWW(0),
      fZDC(-999.0),
      pZDCvsV0MAmp(0),
      pZDCvsTPCFull(0),
      pZDCvsTPCEtaGap(0),
      pZDCvsSPDFull(0),
      pZDCvsSPDEtaGap(0),
      pZDCvsSPDEtaAdj(0),
      pZDCvsSPDEtaGapW(0),
      pZDCvsEtFull(0),
      pZDCvsEtEtaGap(0) {
  for (int i = 0; i < v0m_Nbins; ++i) {
    hDCAxyData[i] = 0;
  }
}
//_____________________________________________________________________________
AliAnalysisTaskDataSpeedOfSound::AliAnalysisTaskDataSpeedOfSound(
    const char* name)
    : AliAnalysisTaskSE(name),
      fESD(0),
      fEventCuts(0x0),
      fMCStack(0),
      fMC(0),
      fUseMC(kFALSE),
      fIsSystematics(true),
      fVaryVtxZPos(false),
      fMinVtxZPos(-5.0),
      fMaxVtxZPos(5.0),
      fSystematic(1),
      fTrigger(AliVEvent::kCentral),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fOutputList(0),
      fEtaCut(0.8),
      fEtaCutForpTwSPDGap(0.4),
      fEtaCutSPDGapMin(0.7),
      fEtaCutSPDGapMax(1.4),
      fEtaCutForpTwTPCGap(0.3),
      fEtaCutTPCGapMin(0.5),
      fEtaCutTPCGapMax(0.8),
      fPtMin(0.15),
      fPtMinCent(0.15),
      fPtMaxCent(10.0),
      fV0Mmin(0.0),
      fV0Mmax(80.0),
      fHMCut(10.0),
      ftrackmult08(0),
      fv0mpercentile(0),
      fv0mamplitude(0),
      fTrackletsEtaGap(0),
      fTracksEtaGapTPC(0),
      fMultSelection(0x0),
      hNch(0),
      hNchvsV0MAmp(0),
      hV0MvsV0MAmp(0),
      pV0MAmpChannel(0),
      hV0MAmplitude(0),
      hV0Percentile(0),
      hPtvsNch(0),
      pPtvsNch(0),
      pPtvsV0MAmp(0),
      hPtvsV0MAmp(0),
      hPhiEtaGapTPC(0),
      hPtWithCutForCent(0),
      hPhiEtaGapSPD(0),
      hBestVtxZ(0),
      hTrackletsEtaGap(0),
      hTracksEtaGapTPC(0),
      hPtvsTrackletsEtaGap(0),
      hPtvsTracksEtaGapTPC(0),
      pPtvsTrackletsEtaGap(0),
      pPtvsTracksEtaGapTPC(0),
      hSPDFull(0),
      hSPDEtaAdj(0),
      hSPDEtaGapW(0),
      hSPDEtaGapWW(0),
      hEtFull(0),
      hEtEtaGap(0),
      hPtvsSPDFull(0),
      hPtvsSPDEtaAdj(0),
      hPtvsSPDEtaGapW(0),
      hPtvsSPDEtaGapWW(0),
      hPtvsEtFull(0),
      hPtvsEtEtaGap(0),
      hPtvsTPCEtaGapWidepT(0),
      hPtvsEtEtaGapWidepT(0),
      pZVtxvsSPDClus(0),
      pSPDClusvsEta(0),
      fSPDVtxCut(3.0),
      fSPDFull(0),
      fSPDEtaAdj(0),
      fSPDEtaGapW(0),
      fSPDEtaGapWW(0),
      fZDC(-999.0),
      pZDCvsV0MAmp(0),
      pZDCvsTPCFull(0),
      pZDCvsTPCEtaGap(0),
      pZDCvsSPDFull(0),
      pZDCvsSPDEtaGap(0),
      pZDCvsSPDEtaAdj(0),
      pZDCvsSPDEtaGapW(0),
      pZDCvsEtFull(0),
      pZDCvsEtEtaGap(0) {
  for (int i = 0; i < v0m_Nbins; ++i) {
    hDCAxyData[i] = 0;
  }
  DefineInput(0, TChain::Class());  // define the input of the analysis: in this
                                    // case you take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about
  // it, does its work automatically
  DefineOutput(1, TList::Class());  // define the ouptut of the analysis: in
                                    // this case it's a list of histograms
}
//_____________________________________________________________________________
AliAnalysisTaskDataSpeedOfSound::~AliAnalysisTaskDataSpeedOfSound() {
  // destructor
  if (fOutputList) {
    delete fOutputList;  // at the end of your task, it is deleted from memory
                         // by calling this function
    fOutputList = 0x0;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound::UserCreateOutputObjects() {
  if (!fTrackFilter) {
    fTrackFilter = new AliAnalysisFilter("trackFilter2015");
    AliESDtrackCuts* fCuts = new AliESDtrackCuts();
    fCuts->SetMaxFractionSharedTPCClusters(0.4);
    fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
    fCuts->SetMaxChi2PerClusterTPC(4);
    fCuts->SetAcceptKinkDaughters(kFALSE);
    fCuts->SetRequireTPCRefit(kTRUE);
    fCuts->SetRequireITSRefit(kTRUE);
    fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                    AliESDtrackCuts::kAny);
    fCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    fCuts->SetMaxChi2TPCConstrainedGlobal(36);
    fCuts->SetMaxDCAToVertexZ(2);
    fCuts->SetDCAToVertex2D(kFALSE);
    fCuts->SetRequireSigmaToVertex(kFALSE);
    fCuts->SetMaxChi2PerClusterITS(36);
    fCuts->SetEtaRange(-0.8, 0.8);

    if (fIsSystematics) {
      ChangeCut(fCuts);
    }
    fTrackFilter->AddCuts(fCuts);
  }

  // track cuts to find contamination via DCA distribution
  if (!fTrackFilterwoDCA) {
    fTrackFilterwoDCA = new AliAnalysisFilter("trackFilter2015");
    AliESDtrackCuts* fCuts3 = new AliESDtrackCuts();
    fCuts3->SetMaxFractionSharedTPCClusters(0.4);                //
    fCuts3->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);  //
    fCuts3->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);           //
    fCuts3->SetMaxChi2PerClusterTPC(4);                          //
    fCuts3->SetAcceptKinkDaughters(kFALSE);                      //
    fCuts3->SetRequireTPCRefit(kTRUE);                           //
    fCuts3->SetRequireITSRefit(kTRUE);                           //
    fCuts3->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                     AliESDtrackCuts::kAny);  //
    // fCuts3->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//
    // fCuts3->SetMaxChi2TPCConstrainedGlobal(36);//
    fCuts3->SetMaxDCAToVertexZ(2);            //
    fCuts3->SetDCAToVertex2D(kFALSE);         //
    fCuts3->SetRequireSigmaToVertex(kFALSE);  //
    fCuts3->SetMaxChi2PerClusterITS(36);      //
    fCuts3->SetEtaRange(-0.8, 0.8);

    if (fIsSystematics) {
      ChangeCut(fCuts3);
    }
    fTrackFilterwoDCA->AddCuts(fCuts3);
  }

  if (fVaryVtxZPos) {
    fEventCuts.SetManualMode();  //! Enable manual mode
    fEventCuts.fMinVtz = fMinVtxZPos;
    fEventCuts.fMaxVtz = fMaxVtxZPos;
  }

  // create output objects
  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  const int nPtbins = 12;
  double Ptbins[nPtbins + 1] = {0.15, 0.2, 0.4, 0.6,  0.8,  1.0, 1.5,
                                2.0,  3.0, 5.0, 10.0, 20.0, 50.0};

  constexpr int pt_Nbins{210};
  double pt_bins[pt_Nbins + 1] = {0};
  for (int i = 0; i <= pt_Nbins; ++i) {
    pt_bins[i] = 0.15 + (i * 0.05);
  }

  // Nch (|eta|<0.8)
  constexpr int nch_Nbins{1250};
  double nch_bins[nch_Nbins + 1] = {0};
  for (int i = 0; i <= nch_Nbins; ++i) {
    nch_bins[i] = 0.5 + (2.0 * i);
  }

  // Nch (0.5<|eta|<0.8)
  constexpr int nchEtaGapTPC_Nbins{500};
  double nchEtaGapTPC_bins[nchEtaGapTPC_Nbins + 1] = {0};
  for (int i = 0; i <= nchEtaGapTPC_Nbins; ++i) {
    nchEtaGapTPC_bins[i] = 0.5 + (2.0 * i);
  }

  // SPD Tracklets (|eta|<0.8)
  constexpr int SPD0p8_Nbins{1550};
  double SPD0p8_bins[SPD0p8_Nbins + 1] = {0};
  for (int i = 0; i <= SPD0p8_Nbins; ++i) {
    SPD0p8_bins[i] = 0.5 + (2.0 * i);
  }

  // SPD Tracklets (|eta|<0.4)
  constexpr int SPD0p4_Nbins{850};
  double SPD0p4_bins[SPD0p4_Nbins + 1] = {0};
  for (int i = 0; i <= SPD0p4_Nbins; ++i) {
    SPD0p4_bins[i] = 0.5 + (2.0 * i);
  }

  // nTracklets (0.5<|eta|<0.8)
  constexpr int SPDEtaGap_Nbins{700};
  double SPDEtaGap_bins[SPDEtaGap_Nbins + 1] = {0};
  for (int i = 0; i <= SPDEtaGap_Nbins; ++i) {
    SPDEtaGap_bins[i] = 0.5 + (2.0 * i);
  }

  // Et binning
  constexpr int Et_Nbins{1000};
  double Et_bins[Et_Nbins + 1] = {0};
  for (int i = 0; i <= Et_Nbins; ++i) {
    Et_bins[i] = 0.0 + (2.0 * i);
  }

  constexpr int EtEtaGap_Nbins{580};
  double EtEtaGap_bins[EtEtaGap_Nbins + 1] = {0};
  for (int i = 0; i <= EtEtaGap_Nbins; ++i) {
    EtEtaGap_bins[i] = 0.0 + (2.0 * i);
  }

  constexpr double v0mAmp_width{25.0};
  constexpr int v0mAmp_Nbins{1720};
  double v0mAmp_bins[v0mAmp_Nbins + 1] = {0};
  for (int i = 0; i <= v0mAmp_Nbins; ++i) {
    v0mAmp_bins[i] = 0.0 + i * v0mAmp_width;
  }

  constexpr int dcaxy_Nbins{100};
  double dcaxy_bins[dcaxy_Nbins + 1] = {0};
  for (int i = 0; i <= dcaxy_Nbins; ++i) {
    dcaxy_bins[i] = -3.0 + (0.06 * i);
  }

  constexpr int v0m_Nbins080{6};
  constexpr double v0m_bins080[v0m_Nbins080 + 1] = {0.0,  1.0,  5.0, 10.0,
                                                    20.0, 50.0, 80.0};

  hV0Percentile =
      new TH1F("hV0Percentile", ";V0M (%);Entries", v0m_Nbins080, v0m_bins080);

  hNch = new TH1F("hTPCFull", ";#it{N}_{ch} (|#eta|#leq0.8); Counts;",
                  nch_Nbins, nch_bins);

  hNchvsV0MAmp =
      new TH2F("hNchvsV0MAmp", ";#it{N}_{ch} (|#eta|#leq0.8); V0M amplitude",
               nch_Nbins, nch_bins, v0mAmp_Nbins, v0mAmp_bins);

  hV0MvsV0MAmp = new TH2F("hV0MvsV0MAmp", ";V0M amplitude; V0M percentile",
                          v0mAmp_Nbins, v0mAmp_bins, v0m_Nbins080, v0m_bins080);

  pV0MAmpChannel =
      new TProfile("pV0MAmpChannel", ";V0 Channel; Amplitude;", 64, -0.5, 63.5);

  hV0MAmplitude =
      new TH1F("hV0MAmp", ";V0M amplitude; Counts", v0mAmp_Nbins, v0mAmp_bins);

  hPtvsNch = new TH2D("hPtvsTPCFull",
                      "; #it{N}_{ch} (|#eta|#leq0.8); #it{p}_{T} "
                      "(|#eta|#leq0.8, GeV/#it{c});",
                      nch_Nbins, nch_bins, pt_Nbins, pt_bins);

  pPtvsNch = new TProfile("pPtvsTPCFull",
                          "; #it{N}_{ch} (|#eta|#leq0.8); #LT#it{p}_{T}#GT "
                          "(|#eta|#leq0.8, GeV/#it{c})",
                          nch_Nbins, nch_bins);

  pPtvsV0MAmp = new TProfile(
      "pPtvsV0MAmp",
      "; V0M Amplitude; #LT#it{p}_{T}#GT (|#eta|#leq0.8, GeV/#it{c})",
      v0mAmp_Nbins, v0mAmp_bins);

  hPtvsV0MAmp = new TH2D(
      "hPtvsV0MAmp", ";V0M amplitude; #it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      v0mAmp_Nbins, v0mAmp_bins, pt_Nbins, pt_bins);

  hPhiEtaGapTPC =
      new TH2F("hPhiEtaGapTPC", ";#varphi; #eta (0.5#leq|#eta|#leq0.8)", 30, 0,
               2 * TMath::Pi(), 20, -1.0, 1.0);

  hPtWithCutForCent = new TH1F(
      "hPtWithCutForCent",
      Form("pT cut: %f - %f; pT (GeV/c); Counts", fPtMinCent, fPtMaxCent),
      pt_Nbins, pt_bins);

  hPhiEtaGapSPD =
      new TH2F("hPhiEtaGapSPD", ";#varphi; #eta (0.5#leq|#eta|#leq0.8) ;", 30,
               0, 2 * TMath::Pi(), 20, -1.0, 1.0);

  hBestVtxZ =
      new TH1F("hBestVtxZ", ";Vertex_{#it{z}} (cm); Counts;", 400, -11, 11);

  hTrackletsEtaGap = new TH1F(
      "hSPDEtaGap", "; #it{N}_{tracklet} (0.5#leq|#eta|#leq0.8); Entries",
      SPDEtaGap_Nbins, SPDEtaGap_bins);

  hTracksEtaGapTPC =
      new TH1F("hTPCEtaGap", "; #it{N}_{ch} (0.5#leq|#eta|#leq0.8); Entries",
               nchEtaGapTPC_Nbins, nchEtaGapTPC_bins);

  hPtvsTrackletsEtaGap =
      new TH2D("hPtvsSPDEtaGap",
               "; #it{N}_{tracklet} (0.5#leq|#eta|#leq0.8); #it{p}_{T} "
               "(|#eta|#leq0.4, GeV/#it{c})",
               SPDEtaGap_Nbins, SPDEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsTracksEtaGapTPC = new TH2D(
      "hPtvsTPCEtaGap",
      "; #it{N}_{ch} (0.5#leq#eta#leq0.8); #it{p}_{T} (|#eta|<0.3, GeV/#it{c})",
      nchEtaGapTPC_Nbins, nchEtaGapTPC_bins, pt_Nbins, pt_bins);

  pPtvsTrackletsEtaGap =
      new TProfile("pPtvsSPDEtaGap",
                   "; #it{N}_{tracklet} (0.5leq|#eta|#leq0.8); "
                   "#LT#it{p}_{T}#GT (|#eta|#leq0.4, GeV/#it{c})",
                   SPDEtaGap_Nbins, SPDEtaGap_bins);

  pPtvsTracksEtaGapTPC =
      new TProfile("pPtvsTPCEtaGap",
                   "; #it{N}_{nch} (0.5#leq|#eta|#leq0.8); "
                   "#LT#it{p}_{T}#GT (|#eta|#leq0.3, GeV/#it{c})",
                   nchEtaGapTPC_Nbins, nchEtaGapTPC_bins);

  hSPDFull = new TH1F("hSPDFull", "; #it{N}_{ch} (|#eta|#leq0.8); Entries",
                      SPD0p8_Nbins, SPD0p8_bins);

  hSPDEtaAdj =
      new TH1F("hSPDEtaAdj", "; #it{N}_{ch} (0.4<|#eta|#leq0.8); Entries",
               SPD0p4_Nbins, SPD0p4_bins);

  hSPDEtaGapW =
      new TH1F("hSPDEtaGapW", "; #it{N}_{ch} (0.7#leq|#eta|#leq1); Entries",
               SPDEtaGap_Nbins, SPDEtaGap_bins);

  hSPDEtaGapWW =
      new TH1F("hSPDEtaGapWW", "; #it{N}_{ch} (1#leq|#eta|#leq1.3); Entries",
               SPDEtaGap_Nbins, SPDEtaGap_bins);

  hEtFull = new TH1F("hEtFull", ";#it{E}_{T} (|#eta|#leq0.8); Entries",
                     Et_Nbins, Et_bins);

  hEtEtaGap =
      new TH1F("hEtEtaGap", ";#it{E}_{T} (0.5#leq|#eta|#leq0.8);Entries",
               EtEtaGap_Nbins, EtEtaGap_bins);

  hPtvsSPDFull = new TH2D("hPtvsSPDFull",
                          "; #it{N}_{ch} (|#eta|#leq0.8); #it{p}_{T} "
                          "(|#eta|#leq0.8, GeV/#it{c})",
                          SPD0p8_Nbins, SPD0p8_bins, pt_Nbins, pt_bins);

  hPtvsSPDEtaAdj = new TH2D("hPtvsSPDEtaAdj",
                            "; #it{N}_{ch} (0.4<|#eta|#leq0.8); #it{p}_{T} "
                            "(|#eta|#leq0.4, GeV/#it{c})",
                            SPD0p4_Nbins, SPD0p4_bins, pt_Nbins, pt_bins);

  hPtvsSPDEtaGapW =
      new TH2D("hPtvsSPDEtaGapW",
               "; #it{N}_{ch} (0.7#leq|#eta|#leq1); #it{p}_{T} "
               "(|#eta|#leq0.3, GeV/#it{c})",
               SPDEtaGap_Nbins, SPDEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsSPDEtaGapWW =
      new TH2D("hPtvsSPDEtaGapWW",
               "; #it{N}_{ch} (1#leq|#eta|#leq1.3); #it{p}_{T} "
               "(|#eta|#leq0.3, GeV/#it{c})",
               SPDEtaGap_Nbins, SPDEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsEtFull = new TH2D(
      "hPtvsEtFull",
      ";#it{E}_{T} (|#eta|#leq0.8); #it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      Et_Nbins, Et_bins, pt_Nbins, pt_bins);

  hPtvsEtEtaGap = new TH2D("hPtvsEtEtaGap",
                           ";#it{E}_{T} (0.5#leq|#eta|#leq0.8); #it{p}_{T} "
                           "(|#eta|#leq0.3, GeV/#it{c})",
                           EtEtaGap_Nbins, EtEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsTPCEtaGapWidepT = new TH2F(
      "hPtvsTPCEtaGapWidepT",
      "; #it{N}_{ch} (0.5#leq#eta#leq0.8); #it{p}_{T} (|#eta|<0.3, GeV/#it{c})",
      nchEtaGapTPC_Nbins, nchEtaGapTPC_bins, nPtbins, Ptbins);

  hPtvsEtEtaGapWidepT =
      new TH2F("hPtvsEtEtaGapWidepT",
               "; #it{E}_{T} (0.5#leq|#eta|#leq0.8); #it{p}_{T} "
               "(|#eta|#leq0.3, GeV/#it{c});",
               EtEtaGap_Nbins, EtEtaGap_bins, nPtbins, Ptbins);

  pZVtxvsSPDClus =
      new TProfile("pZVtxvsSPDClus", ";Z_{SPD} Vertex (cm); <SPD clusters>;",
                   40, -10.0, 10.0);

  pSPDClusvsEta =
      new TProfile("pSPDClusvsEta", ";#eta; <SPD clusters>;", 30, -1.5, 1.5);

  pZDCvsV0MAmp = new TProfile("pZDCvsV0MAmp", ";V0M Amp;<ZDC>;", v0mAmp_Nbins,
                              v0mAmp_bins);
  pZDCvsTPCFull =
      new TProfile("pZDCvsTPCFull", ";Nch;<ZDC>;", nch_Nbins, nch_bins);
  pZDCvsTPCEtaGap =
      new TProfile("pZDCvsTPCEtaGap", ";Nch;<ZDC>;", nch_Nbins, nch_bins);
  pZDCvsSPDFull =
      new TProfile("pZDCvsSPDFull", ";Nch;<ZDC>;", SPD0p8_Nbins, SPD0p8_bins);
  pZDCvsSPDEtaGap = new TProfile("pZDCvsSPDEtaGap", ";Nch;<ZDC>;",
                                 SPDEtaGap_Nbins, SPDEtaGap_bins);
  pZDCvsSPDEtaAdj =
      new TProfile("pZDCvsSPDEtaAdj", ";Nch;<ZDC>;", SPD0p4_Nbins, SPD0p4_bins);
  pZDCvsSPDEtaGapW = new TProfile("pZDCvsSPDEtaGapW", ";Nch;<ZDC>;",
                                  SPDEtaGap_Nbins, SPDEtaGap_bins);
  pZDCvsEtFull =
      new TProfile("pZDCvsEtFull", ";E_{T};<ZDC>;", Et_Nbins, Et_bins);
  pZDCvsEtEtaGap = new TProfile("pZDCvsEtEtaGap", ";E_{T};<ZDC>;",
                                EtEtaGap_Nbins, EtEtaGap_bins);

  fOutputList->Add(hBestVtxZ);
  fOutputList->Add(pZVtxvsSPDClus);
  fOutputList->Add(pSPDClusvsEta);
  fOutputList->Add(hV0Percentile);
  fOutputList->Add(hNchvsV0MAmp);
  fOutputList->Add(hV0MvsV0MAmp);

  fOutputList->Add(hV0MAmplitude);
  fOutputList->Add(hPtvsV0MAmp);
  fOutputList->Add(pPtvsV0MAmp);

  fOutputList->Add(hNch);
  fOutputList->Add(hPtvsNch);
  fOutputList->Add(pPtvsNch);

  fOutputList->Add(hTrackletsEtaGap);
  fOutputList->Add(hPtvsTrackletsEtaGap);
  fOutputList->Add(pPtvsTrackletsEtaGap);

  fOutputList->Add(hTracksEtaGapTPC);
  fOutputList->Add(hPtvsTracksEtaGapTPC);
  fOutputList->Add(pPtvsTracksEtaGapTPC);

  fOutputList->Add(hPhiEtaGapSPD);
  fOutputList->Add(hPhiEtaGapTPC);
  fOutputList->Add(pV0MAmpChannel);
  fOutputList->Add(hPtWithCutForCent);

  fOutputList->Add(hSPDFull);
  fOutputList->Add(hPtvsSPDFull);

  fOutputList->Add(hSPDEtaAdj);
  fOutputList->Add(hPtvsSPDEtaAdj);

  fOutputList->Add(hSPDEtaGapW);
  fOutputList->Add(hPtvsSPDEtaGapW);

  fOutputList->Add(hSPDEtaGapWW);
  fOutputList->Add(hPtvsSPDEtaGapWW);

  fOutputList->Add(hEtFull);
  fOutputList->Add(hPtvsEtFull);

  fOutputList->Add(hEtEtaGap);
  fOutputList->Add(hPtvsEtEtaGap);

  fOutputList->Add(hPtvsTPCEtaGapWidepT);
  fOutputList->Add(hPtvsEtEtaGapWidepT);

  fOutputList->Add(pZDCvsV0MAmp);
  fOutputList->Add(pZDCvsTPCFull);
  fOutputList->Add(pZDCvsTPCEtaGap);
  fOutputList->Add(pZDCvsSPDFull);
  fOutputList->Add(pZDCvsSPDEtaGap);
  fOutputList->Add(pZDCvsSPDEtaAdj);
  fOutputList->Add(pZDCvsSPDEtaGapW);
  fOutputList->Add(pZDCvsEtFull);
  fOutputList->Add(pZDCvsEtEtaGap);

  for (int i = 0; i < v0m_Nbins; ++i) {
    hDCAxyData[i] = new TH2F(Form("hDCAxyData_%s", uc_v0m_bins_name[i]),
                             ";DCA_{xy} (cm);#it{p}_{T} (GeV/#it{c});",
                             dcaxy_Nbins, dcaxy_bins, pt_Nbins, pt_bins);
    fOutputList->Add(hDCAxyData[i]);
  }

  fEventCuts.AddQAplotsToList(fOutputList);
  PostData(1, fOutputList);  // postdata will notify the analysis manager of
                             // changes / updates to the
}
//_____________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound::UserExec(Option_t*) {
  AliVEvent* event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  fESD = dynamic_cast<AliESDEvent*>(event);

  if (!fESD) {
    Printf("%s:%d ESDEvent not found in Input Manager", (char*)__FILE__,
           __LINE__);
    this->Dump();
    return;
  }

  if (fUseMC) {
    //      E S D
    fMC = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!fMC) {
      Printf("%s:%d MCEvent not found in Input Manager", (char*)__FILE__,
             __LINE__);
      this->Dump();
      return;
    }
    fMCStack = fMC->Stack();
  }

  ftrackmult08 = -999.0;
  fv0mpercentile = -999.0;
  fv0mamplitude = -999.0;

  fMultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
  fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  ftrackmult08 = AliESDtrackCuts::GetReferenceMultiplicity(
      fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);

  //! Analyze only the 0--80 % V0M range
  if (!(fv0mpercentile >= fV0Mmin && fv0mpercentile < fV0Mmax)) {
    return;
  }

  //! Trigger selection
  bool isEventTriggered{false};
  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  isEventTriggered = fSelectMask & fTrigger;
  if (!isEventTriggered) {
    return;
  }

  // Good events
  if (!fEventCuts.AcceptEvent(event)) {
    PostData(1, fOutputList);
    return;
  }

  // Good vertex
  bool hasRecVertex{false};
  hasRecVertex = HasRecVertex();
  if (!hasRecVertex) {
    return;
  }

  VertexPosition();

  //! Get calibrated V0 amplitude
  GetCalibratedV0Amplitude();

  //! Get SPD tracklets multiplicity
  GetSPDMultiplicity();

  //! Get ZDC Centrality
  GetZDC();

  //! DCAxy templates MC and Data
  DCAxyDistributions();

  //! Data Multiplicity distributions
  MultiplicityDistributions();

  PostData(1, fOutputList);
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSound::Terminate(Option_t*) {}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSound::VertexPosition() {
  //! best primary vertex available
  const AliVVertex* vtx = fEventCuts.GetPrimaryVertex();
  hBestVtxZ->Fill(vtx->GetZ());
}

//______________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound::GetZDC() {
  AliESDZDC* esdZDC = fESD->GetESDZDC();
  if (!esdZDC) {
    return;
  }

  fZDC = -999.0;
  double zp = esdZDC->GetZDCP1Energy() + esdZDC->GetZDCP2Energy();
  double zn = esdZDC->GetZDCN1Energy() + esdZDC->GetZDCN2Energy();
  zp *= 0.001;
  zn *= 0.001;
  fZDC = zp + zn;
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSound::GetCalibratedV0Amplitude() {
  float mV0M{0.0};
  for (int i = 0; i < 64; i++) {
    mV0M += fESD->GetVZEROEqMultiplicity(i);
    pV0MAmpChannel->Fill(i, fESD->GetVZEROEqMultiplicity(i));
  }
  fv0mamplitude = mV0M;
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSound::GetSPDMultiplicity() {
  fTrackletsEtaGap = 0;

  fSPDFull = 0;
  fSPDEtaAdj = 0;
  fSPDEtaGapW = 0;
  fSPDEtaGapWW = 0;

  int nTracklets = 0;
  float spdVtxZ = -999.0;
  AliMultiplicity* SPDptr = fESD->GetMultiplicity();
  if (!SPDptr) {
    return;
  }
  const AliESDVertex* spdVtx = fESD->GetPrimaryVertexSPD();
  if (!spdVtx) {
    return;
  }
  if (spdVtx->GetNContributors() <= 2) {
    return;
  }
  spdVtxZ = spdVtx->GetZ();
  if (TMath::Abs(spdVtxZ) > fSPDVtxCut) {
    return;
  }

  pZVtxvsSPDClus->Fill(spdVtxZ, SPDptr->GetNumberOfSPDClusters());

  nTracklets = SPDptr->GetNumberOfTracklets();
  for (auto it = 0; it < nTracklets; it++) {
    double eta = SPDptr->GetEta(it);
    double phi = SPDptr->GetPhi(it);

    pSPDClusvsEta->Fill(eta, SPDptr->GetNumberOfSPDClusters());

    if (TMath::Abs(eta) <= fEtaCut) {
      fSPDFull++;
    }

    if (TMath::Abs(eta) > 0.4 && TMath::Abs(eta) <= fEtaCut) {
      fSPDEtaAdj++;
    }

    if (TMath::Abs(eta) > 0.7 && TMath::Abs(eta) <= 1.0) {
      fSPDEtaGapW++;
    }

    if (TMath::Abs(eta) > 1.0 && TMath::Abs(eta) <= 1.3) {
      fSPDEtaGapWW++;
    }

    if (TMath::Abs(eta) >= fEtaCutSPDGapMin &&
        TMath::Abs(eta) <= fEtaCutSPDGapMax) {
      fTrackletsEtaGap++;
      hPhiEtaGapSPD->Fill(phi, eta);
    }
  }
}
//______________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound::MultiplicityDistributions() {
  int rec_nch{0};
  double etfull{0.0};
  double etetagap{0.0};
  fTracksEtaGapTPC = 0;
  const double masspi{0.13957};
  const int n_tracks{fESD->GetNumberOfTracks()};

  for (int i = 0; i < n_tracks; ++i) {
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
    if (!track) {
      continue;
    }
    if (!fTrackFilter->IsSelected(track)) {
      continue;
    }
    if (track->Pt() < fPtMinCent || track->Pt() > fPtMaxCent) {
      continue;
    }
    if (track->Charge() == 0) {
      continue;
    }
    if (TMath::Abs(track->Eta()) > fEtaCut) {
      continue;
    }

    if (TMath::Abs(track->Eta()) >= fEtaCutTPCGapMin &&
        TMath::Abs(track->Eta()) <= fEtaCutTPCGapMax) {
      fTracksEtaGapTPC++;
      hPhiEtaGapTPC->Fill(track->Phi(), track->Eta());
      etetagap += TMath::Sqrt(TMath::Power(track->Pt(), 2.0) +
                              TMath::Power(masspi, 2.0));
    }
    hPtWithCutForCent->Fill(track->Pt());
    etfull +=
        TMath::Sqrt(TMath::Power(track->Pt(), 2.0) + TMath::Power(masspi, 2.0));
    rec_nch++;
  }

  for (int i = 0; i < n_tracks; ++i) {
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
    if (!track) {
      continue;
    }
    if (!fTrackFilter->IsSelected(track)) {
      continue;
    }
    if (track->Pt() < fPtMin) {
      continue;
    }
    if (track->Charge() == 0) {
      continue;
    }
    if (TMath::Abs(track->Eta()) > fEtaCut) {
      continue;
    }
    double pt = track->Pt();

    //! Nch |eta|<=0.8 and Spectra |eta|<=0.8
    hPtvsNch->Fill(rec_nch, pt);
    pPtvsNch->Fill(rec_nch, pt);

    //! Nch in V0 and Spectra |eta|<=0.8
    hPtvsV0MAmp->Fill(fv0mamplitude, pt);
    pPtvsV0MAmp->Fill(fv0mamplitude, pt);

    hPtvsSPDFull->Fill(fSPDFull, pt);
    hPtvsEtFull->Fill(etfull, pt);

    //! SPD Eta Adjacent
    if (TMath::Abs(track->Eta()) <= 0.4) {
      hPtvsSPDEtaAdj->Fill(fSPDEtaAdj, pt);
    }

    //! SPD Eta Gap
    if (TMath::Abs(track->Eta()) <= fEtaCutForpTwSPDGap) {
      hPtvsTrackletsEtaGap->Fill(fTrackletsEtaGap, pt);
      pPtvsTrackletsEtaGap->Fill(fTrackletsEtaGap, pt);
      hPtvsSPDEtaGapW->Fill(fSPDEtaGapW, pt);
      hPtvsSPDEtaGapWW->Fill(fSPDEtaGapWW, pt);
    }

    //! Nch 0.5<=|eta|<=0.8 and Spectra |eta|<=0.3
    if (TMath::Abs(track->Eta()) <= fEtaCutForpTwTPCGap) {
      hPtvsTracksEtaGapTPC->Fill(fTracksEtaGapTPC, pt);
      pPtvsTracksEtaGapTPC->Fill(fTracksEtaGapTPC, pt);
      hPtvsTPCEtaGapWidepT->Fill(fTracksEtaGapTPC, pt);
      hPtvsEtEtaGap->Fill(etetagap, pt);
      hPtvsEtEtaGapWidepT->Fill(etetagap, pt);
    }
  }

  hV0Percentile->Fill(fv0mpercentile);
  hV0MAmplitude->Fill(fv0mamplitude);
  hNchvsV0MAmp->Fill(rec_nch, fv0mamplitude);
  hV0MvsV0MAmp->Fill(fv0mamplitude, fv0mpercentile);
  hNch->Fill(rec_nch);

  hTrackletsEtaGap->Fill(fTrackletsEtaGap);
  hTracksEtaGapTPC->Fill(fTracksEtaGapTPC);

  hSPDFull->Fill(fSPDFull);
  hSPDEtaAdj->Fill(fSPDEtaAdj);
  hSPDEtaGapW->Fill(fSPDEtaGapW);
  hSPDEtaGapWW->Fill(fSPDEtaGapWW);

  hEtFull->Fill(etfull);
  hEtEtaGap->Fill(etetagap);

  pZDCvsV0MAmp->Fill(fv0mamplitude, fZDC);
  pZDCvsTPCFull->Fill(rec_nch, fZDC);
  pZDCvsTPCEtaGap->Fill(fTracksEtaGapTPC, fZDC);
  pZDCvsSPDFull->Fill(fSPDFull, fZDC);
  pZDCvsSPDEtaGap->Fill(fTrackletsEtaGap, fZDC);
  pZDCvsSPDEtaAdj->Fill(fSPDEtaAdj, fZDC);
  pZDCvsSPDEtaGapW->Fill(fSPDEtaGapW, fZDC);
  pZDCvsEtFull->Fill(etfull, fZDC);
  pZDCvsEtEtaGap->Fill(etetagap, fZDC);
}

//____________________________________________________________

void AliAnalysisTaskDataSpeedOfSound::DCAxyDistributions() const {
  int index{-1};
  for (int i = 0; i < v0m_Nbins; ++i) {
    if (fv0mpercentile >= uc_v0m_bins_low[i] &&
        fv0mpercentile < uc_v0m_bins_high[i]) {
      index = i;
      break;
    }
  }

  if (index < 0) {
    return;
  }

  const int n_tracks{fESD->GetNumberOfTracks()};
  for (int i = 0; i < n_tracks; ++i) {
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
    if (!track) {
      continue;
    }
    if (!fTrackFilterwoDCA->IsSelected(track)) {
      continue;
    }
    if (track->Pt() < fPtMin) {
      continue;
    }
    if (TMath::Abs(track->Eta()) > fEtaCut) {
      continue;
    }
    if (track->Charge() == 0) {
      continue;
    }

    float dcaxy = -999;
    float dcaz = -999;
    track->GetImpactParameters(dcaxy, dcaz);
    hDCAxyData[index]->Fill(dcaxy, track->Pt());
  }
}

//____________________________________________________________

Bool_t AliAnalysisTaskDataSpeedOfSound::HasRecVertex() {
  float fMaxDeltaSpdTrackAbsolute = 0.5f;
  float fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  float fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  float fMaxResolutionSPDvertex = 0.25f;
  float fMaxDispersionSPDvertex = 1.e14f;

  Bool_t fRequireTrackVertex = true;
  unsigned long fFlag;
  fFlag = BIT(AliEventCuts::kNoCuts);

  const AliVVertex* vtTrc = fESD->GetPrimaryVertex();
  bool isTrackV = true;
  if (vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ()) isTrackV = false;
  const AliVVertex* vtSPD = fESD->GetPrimaryVertexSPD();

  if (vtSPD->GetNContributors() > 0) fFlag |= BIT(AliEventCuts::kVertexSPD);

  if (vtTrc->GetNContributors() > 1 && isTrackV)
    fFlag |= BIT(AliEventCuts::kVertexTracks);

  if (((fFlag & BIT(AliEventCuts::kVertexTracks)) || !fRequireTrackVertex) &&
      (fFlag & BIT(AliEventCuts::kVertexSPD)))
    fFlag |= BIT(AliEventCuts::kVertex);

  const AliVVertex*& vtx =
      bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
  AliVVertex* fPrimaryVertex = const_cast<AliVVertex*>(vtx);
  if (!fPrimaryVertex) return kFALSE;

  /// Vertex quality cuts
  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = bool(fFlag & AliEventCuts::kVertexSPD) &&
                      bool(fFlag & AliEventCuts::kVertexTracks)
                  ? vtTrc->GetZ() - vtSPD->GetZ()
                  : 0.;  /// If one of the two vertices is not available
                         /// this cut is always passed.
  double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
  double errTrc =
      bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  /// vertex dispersion for run1, only for ESD, AOD code to be added here
  const AliESDVertex* vtSPDESD = dynamic_cast<const AliESDVertex*>(vtSPD);
  double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
  if ((TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute &&
       nsigTot <= fMaxDeltaSpdTrackNsigmaSPD &&
       nsigTrc <= fMaxDeltaSpdTrackNsigmaTrack) &&  // discrepancy
                                                    // track-SPD vertex
      (!vtSPD->IsFromVertexerZ() ||
       TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
      (!vtSPD->IsFromVertexerZ() ||
       vtSPDdispersion <= fMaxDispersionSPDvertex)  /// vertex dispersion cut
                                                    /// for run1, only for ESD
      )  // quality cut on vertexer SPD z
    fFlag |= BIT(AliEventCuts::kVertexQuality);

  Bool_t hasVtx = (TESTBIT(fFlag, AliEventCuts::kVertex)) &&
                  (TESTBIT(fFlag, AliEventCuts::kVertexQuality));

  return hasVtx;
}

//____________________________________________________________

void AliAnalysisTaskDataSpeedOfSound::ChangeCut(AliESDtrackCuts* fCuts) {
  cout << "Changing track cut (systematic variation): " << fSystematic << '\n';
  switch (fSystematic) {
    case 0:
      fCuts->SetMaxDCAToVertexZ(1);
      break;
    case 1:
      fCuts->SetMaxDCAToVertexZ(5);
      break;
    case 2:
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
      break;
    case 3:
      fCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
      break;
    case 4:
      fCuts->SetMaxFractionSharedTPCClusters(0.2);
      break;
    case 5:
      fCuts->SetMaxFractionSharedTPCClusters(1);
      break;
    case 6:
      fCuts->SetMaxChi2PerClusterTPC(3);
      break;
    case 7:
      fCuts->SetMaxChi2PerClusterTPC(5);
      break;
    case 8:
      fCuts->SetMaxChi2PerClusterITS(25);
      break;
    case 9:
      fCuts->SetMaxChi2PerClusterITS(49);
      break;
    case 10:
      fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                      AliESDtrackCuts::kNone);
      break;
    case 11:
      fCuts->SetCutGeoNcrNcl(2., 130., 1.5, 0.85, 0.7);
      break;
    case 12:
      fCuts->SetCutGeoNcrNcl(4., 130., 1.5, 0.85, 0.7);
      break;
    case 13:
      fCuts->SetCutGeoNcrNcl(3., 120., 1.5, 0.85, 0.7);
      break;
    case 14:
      fCuts->SetCutGeoNcrNcl(3., 140., 1.5, 0.85, 0.7);
      break;
    case 15:
      fCuts->SetMaxChi2TPCConstrainedGlobal(25);
      break;
    case 16:
      fCuts->SetMaxChi2TPCConstrainedGlobal(49);
      break;
    default:
      cout << "fSystematic not defined!" << '\n';
  }
}
