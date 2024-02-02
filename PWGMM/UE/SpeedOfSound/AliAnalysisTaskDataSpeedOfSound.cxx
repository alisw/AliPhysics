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
      fSystematic(1),
      fTrigger(AliVEvent::kCentral),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fOutputList(0),
      fEtaCut(0.8),
      fEtaMin(-0.8),
      fEtaMax(0.8),
      fEtaGappT(0.4),
      fEtaGapNchMin(0.7),
      fEtaGapNchMax(1.4),
      fEtaGapTPCpT(0.3),
      fEtaGapTPCNchMin(0.5),
      fEtaGapTPCNchMax(0.8),
      fPtMin(0.15),
      fV0Mmin(0.0),
      fV0Mmax(80.0),
      fHMCut(10.0),
      ftrackmult08(0),
      fv0mpercentile(0),
      fv0mamplitude(0),
      fTracklets14(0),
      fTracklets10(0),
      fTrackletsEtaGap(0),
      fTracksEtaGapTPC(0),
      fza(0),
      fzc(0),
      fzn(0),
      fdcaxy(-999),
      fdcaz(-999),
      fMultSelection(0x0),
      hNchvsV0M(0),
      hNchvsV0MAmp(0),
      hV0MvsV0MAmp(0),
      pV0MAmpChannel(0),
      hV0MAmplitude(0),
      hV0Mmult(0),
      pPtvsNch(0),
      pPtEtaNegvsNchEtaPos(0),
      pPtEtaPosvsNchEtaNeg(0),
      pPtvsV0MAmp(0),
      hPtvsV0MAmp(0),
      hNchEtaPosvsNchEtaNeg(0),
      hPtEtaNegvsNchEtaPos(0),
      hPtEtaPosvsNchEtaNeg(0),
      hPhiEtaSPD(0),
      hPhiEtaGapTPC(0),
      hEtaGapSPD(0),
      hVtxZvsTracklets(0),
      hTrackletsvsV0MAmp14(0),
      hTrackletsvsV0MAmp10(0),
      hTrackletsEtaGap(0),
      hTracksEtaGapTPC(0),
      hPtvsTracklets14(0),
      hPtvsTracklets10(0),
      hPtvsTrackletsEtaGap(0),
      hPtvsTracksEtaGapTPC(0),
      pPtvsTracklets14(0),
      pPtvsTracklets10(0),
      pPtvsTrackletsEtaGap(0),
      pPtvsTracksEtaGapTPC(0) {
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
      fSystematic(1),
      fTrigger(AliVEvent::kCentral),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fOutputList(0),
      fEtaCut(0.8),
      fEtaMin(-0.8),
      fEtaMax(0.8),
      fEtaGappT(0.4),
      fEtaGapNchMin(0.7),
      fEtaGapNchMax(1.4),
      fEtaGapTPCpT(0.3),
      fEtaGapTPCNchMin(0.5),
      fEtaGapTPCNchMax(0.8),
      fPtMin(0.15),
      fV0Mmin(0.0),
      fV0Mmax(80.0),
      fHMCut(10.0),
      ftrackmult08(0),
      fv0mpercentile(0),
      fv0mamplitude(0),
      fTracklets14(0),
      fTracklets10(0),
      fTrackletsEtaGap(0),
      fTracksEtaGapTPC(0),
      fza(0),
      fzc(0),
      fzn(0),
      fdcaxy(-999),
      fdcaz(-999),
      fMultSelection(0x0),
      hNchvsV0M(0),
      hNchvsV0MAmp(0),
      hV0MvsV0MAmp(0),
      pV0MAmpChannel(0),
      hV0MAmplitude(0),
      hV0Mmult(0),
      pPtvsNch(0),
      pPtEtaNegvsNchEtaPos(0),
      pPtEtaPosvsNchEtaNeg(0),
      pPtvsV0MAmp(0),
      hPtvsV0MAmp(0),
      hNchEtaPosvsNchEtaNeg(0),
      hPtEtaNegvsNchEtaPos(0),
      hPtEtaPosvsNchEtaNeg(0),
      hPhiEtaSPD(0),
      hPhiEtaGapTPC(0),
      hEtaGapSPD(0),
      hVtxZvsTracklets(0),
      hTrackletsvsV0MAmp14(0),
      hTrackletsvsV0MAmp10(0),
      hTrackletsEtaGap(0),
      hTracksEtaGapTPC(0),
      hPtvsTracklets14(0),
      hPtvsTracklets10(0),
      hPtvsTrackletsEtaGap(0),
      hPtvsTracksEtaGapTPC(0),
      pPtvsTracklets14(0),
      pPtvsTracklets10(0),
      pPtvsTrackletsEtaGap(0),
      pPtvsTracksEtaGapTPC(0) {
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
    fTrackFilterwoDCA->AddCuts(fCuts3);
  }

  // create output objects
  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  constexpr int pt_Nbins{210};
  double pt_bins[pt_Nbins + 1] = {0};
  for (int i = 0; i <= pt_Nbins; ++i) {
    pt_bins[i] = 0.15 + (i * 0.05);
  }

  // Nch (|eta|<0.8) up to 2500
  constexpr double nch_width{2.0};
  constexpr int nch_Nbins{1250};
  double nch_bins[nch_Nbins + 1] = {0};
  for (int i = 0; i <= nch_Nbins; ++i) {
    nch_bins[i] = nch_width * i;
  }

  // nTracklets (|eta|<1.4) up to 10000
  constexpr double tracklets_width{5.0};
  constexpr int tracklets_Nbins{1200};
  double tracklets_bins[tracklets_Nbins + 1] = {0};
  for (int i = 0; i <= tracklets_Nbins; ++i) {
    tracklets_bins[i] = tracklets_width * i;
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

  hV0Mmult =
      new TH1F("hV0Mmult", ";V0M (%);Entries", v0m_Nbins080, v0m_bins080);
  fOutputList->Add(hV0Mmult);

  hNchvsV0M = new TH2F("hNchvsV0M", ";#it{N}_{ch}; V0M (%)", nch_Nbins,
                       nch_bins, v0m_Nbins080, v0m_bins080);

  hNchvsV0MAmp = new TH2F("hNchvsV0MAmp", ";#it{N}_{ch}; V0M Amp", nch_Nbins,
                          nch_bins, v0mAmp_Nbins, v0mAmp_bins);

  hV0MvsV0MAmp = new TH2F("hV0MvsV0MAmp", ";V0M Ampl; V0M (%)", v0mAmp_Nbins,
                          v0mAmp_bins, v0m_Nbins080, v0m_bins080);

  pV0MAmpChannel =
      new TProfile("pV0MAmpChannel", ";V0 Channel; Amplitude;", 64, -0.5, 63.5);

  hV0MAmplitude =
      new TH1F("hV0MAmp", ";V0M Amplitude; Entries", v0mAmp_Nbins, v0mAmp_bins);

  pPtvsNch = new TProfile("pPtvsNch",
                          "; #it{N}_{ch}^{rec}; #LT#it{p}_{T}#GT GeV/#it{c}",
                          nch_Nbins, nch_bins);
  pPtEtaNegvsNchEtaPos =
      new TProfile("pPtEtaNegvsNchEtaPos",
                   "; #it{N}_{ch} (0#leq#eta#leq0.8); #it{p}_{T} (GeV/#it{c})"
                   "(-0.8#leq#eta<0)",
                   nch_Nbins, nch_bins);
  pPtEtaPosvsNchEtaNeg =
      new TProfile("pPtEtaPosvsNchEtaNeg",
                   "; #it{N}_{ch} (-0.8#leq#eta<0); #it{p}_{T} (GeV/#it{c})"
                   "(0#geq#eta#leq0.8)",
                   nch_Nbins, nch_bins);
  pPtvsV0MAmp = new TProfile("pPtvsV0MAmp",
                             "; V0M Amplitude; #LT#it{p}_{T}#GT GeV/#it{c}",
                             v0mAmp_Nbins, v0mAmp_bins);

  hPtvsV0MAmp = new TH2D("hPtvsV0MAmp", ";V0M Amp; #it{p}_{T} (GeV/#it{c})",
                         v0mAmp_Nbins, v0mAmp_bins, pt_Nbins, pt_bins);

  hNchEtaPosvsNchEtaNeg =
      new TH2D("hNchEtaPosvsNchEtaNeg",
               ";#it{N}_{ch} (0#leq#eta#leq0.8); #it{N}_{ch} (-0.8#leq#eta<0)",
               nch_Nbins, nch_bins, nch_Nbins, nch_bins);

  hPtEtaNegvsNchEtaPos =
      new TH2D("hPtEtaNegvsNchEtaPos",
               "; #it{N}_{ch} (0#leq#eta#leq0.8); #it{p}_{T} (GeV/#it{c})"
               "(-0.8#leq#eta<0)",
               nch_Nbins, nch_bins, pt_Nbins, pt_bins);

  hPtEtaPosvsNchEtaNeg =
      new TH2D("hPtEtaPosvsNchEtaNeg",
               "; #it{N}_{ch} (-0.8#leq#eta<0); #it{p}_{T} (GeV/#it{c})"
               "(0#geq#eta#leq0.8)",
               nch_Nbins, nch_bins, pt_Nbins, pt_bins);

  hPhiEtaSPD = new TH2F("hPhiEtaSPD", ";#varphi; #eta (|#eta|#leq1.4)", 80, 0,
                        2 * TMath::Pi(), 80, -1.4, 1.4);

  hPhiEtaGapTPC =
      new TH2F("hPhiEtaGapTPC", ";#varphi; #eta (0.5#leq|#eta|#leq0.8)", 30, 0,
               2 * TMath::Pi(), 80, -1.4, 1.4);

  hEtaGapSPD = new TH1F("hEtaGapSPD", ";#eta (0.7#leq|#eta|#leq1.4) ;Counts",
                        300, -1.5, 1.5);

  hVtxZvsTracklets =
      new TH2F("hVtxZvsTracklets", ";#it{Z}_{vtz} (cm); #it{N}_{tracklets}", 50,
               -10, 10, tracklets_Nbins, tracklets_bins);

  hTrackletsvsV0MAmp14 = new TH2F(
      "hTrackletsvsV0MAmp14", "; #it{N}_{tracklet} (|#eta|#leq1.4); V0M Amp",
      tracklets_Nbins, tracklets_bins, v0mAmp_Nbins, v0mAmp_bins);

  hTrackletsvsV0MAmp10 = new TH2F(
      "hTrackletsvsV0MAmp10", "; #it{N}_{tracklet} (|#eta|#leq1); V0M Amp",
      tracklets_Nbins, tracklets_bins, v0mAmp_Nbins, v0mAmp_bins);

  hTrackletsEtaGap = new TH1F(
      "hTrackletsEtaGap", "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); Entries",
      tracklets_Nbins, tracklets_bins);

  hTracksEtaGapTPC = new TH1F("hTracksEtaGapTPC",
                              "; #it{N}_{ch} (0.5#leq|#eta|#leq0.8); Entries",
                              nch_Nbins, nch_bins);

  hPtvsTracklets14 =
      new TH2D("hPtvsTracklets14",
               "; #it{N}_{tracklet} (|#eta|#leq1.4); #it{p}_{T} (GeV/#it{c})",
               tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);

  hPtvsTracklets10 =
      new TH2D("hPtvsTracklets10",
               "; #it{N}_{tracklet} (|#eta|#leq1); #it{p}_{T} (GeV/#it{c})",
               tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);

  hPtvsTrackletsEtaGap = new TH2D(
      "hPtvsTrackletsEtaGap",
      "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); #it{p}_{T} (GeV/#it{c})",
      tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);

  hPtvsTracksEtaGapTPC = new TH2D(
      "hPtvsTracksEtaGapTPC",
      "; #it{N}_{ch} (0.5#leq#eta#leq0.8); #it{p}_{T} (|#eta|<0.3, GeV/#it{c})",
      nch_Nbins, nch_bins, pt_Nbins, pt_bins);

  pPtvsTracklets14 = new TProfile(
      "pPtvsTracklets14",
      "; #it{N}_{tracklet} (|#eta|#leq1.4); #LT#it{p}_{T}#GT (GeV/#it{c})",
      tracklets_Nbins, tracklets_bins);

  pPtvsTracklets10 = new TProfile(
      "pPtvsTracklets10",
      "; #it{N}_{tracklet} (|#eta|#leq1); #LT#it{p}_{T}#GT (GeV/#it{c})",
      tracklets_Nbins, tracklets_bins);

  pPtvsTrackletsEtaGap =
      new TProfile("pPtvsTrackletsEtaGap",
                   "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); "
                   "#LT#it{p}_{T}#GT (GeV/#it{c})",
                   tracklets_Nbins, tracklets_bins);

  pPtvsTracksEtaGapTPC =
      new TProfile("pPtvsTracksEtaGapTPC",
                   "; #it{N}_{tracklet} (0.5#leq|#eta|#leq0.8); "
                   "#LT#it{p}_{T}#GT (|#eta|#leq0.3, GeV/#it{c})",
                   tracklets_Nbins, tracklets_bins);

  fOutputList->Add(hNchvsV0M);
  fOutputList->Add(hNchvsV0MAmp);
  fOutputList->Add(hV0MvsV0MAmp);
  fOutputList->Add(pV0MAmpChannel);
  fOutputList->Add(hV0MAmplitude);
  fOutputList->Add(pPtvsNch);
  fOutputList->Add(pPtEtaNegvsNchEtaPos);
  fOutputList->Add(pPtEtaPosvsNchEtaNeg);
  fOutputList->Add(pPtvsV0MAmp);
  fOutputList->Add(hPtvsV0MAmp);
  fOutputList->Add(hNchEtaPosvsNchEtaNeg);
  fOutputList->Add(hPtEtaNegvsNchEtaPos);
  fOutputList->Add(hPtEtaPosvsNchEtaNeg);
  fOutputList->Add(hPhiEtaSPD);
  fOutputList->Add(hPhiEtaGapTPC);
  fOutputList->Add(hEtaGapSPD);
  fOutputList->Add(hVtxZvsTracklets);
  fOutputList->Add(hTrackletsvsV0MAmp14);
  fOutputList->Add(hTrackletsvsV0MAmp10);
  fOutputList->Add(hTrackletsEtaGap);
  fOutputList->Add(hTracksEtaGapTPC);
  fOutputList->Add(hPtvsTracklets14);
  fOutputList->Add(pPtvsTracklets14);
  fOutputList->Add(hPtvsTracklets10);
  fOutputList->Add(pPtvsTracklets10);
  fOutputList->Add(hPtvsTrackletsEtaGap);
  fOutputList->Add(pPtvsTrackletsEtaGap);
  fOutputList->Add(hPtvsTracksEtaGapTPC);
  fOutputList->Add(pPtvsTracksEtaGapTPC);

  for (int i = 0; i < v0m_Nbins; ++i) {
    hDCAxyData[i] = new TH2F(Form("hDCAxyData_%s", uc_v0m_bins_name[i]),
                             ";DCA_{xy} (cm);#it{p}_{T} GeV/#it{c}",
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
  bool hasRecVertex = false;
  hasRecVertex = HasRecVertex();
  if (!hasRecVertex) {
    return;
  }

  //! Get calibrated V0 amplitude
  GetCalibratedV0Amplitude();

  //! Get SPD tracklets multiplicity
  GetSPDMultiplicity();

  //! Get ZDC Centrality
  // GetZDCCentrality();

  //! DCAxy templates MC and Data
  DCAxyDistributions();

  //! Data Multiplicity distributions
  MultiplicityDistributions();

  PostData(1, fOutputList);
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSound::Terminate(Option_t*) {}

//______________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound::GetZDCCentrality() {
  AliESDZDC* esdZDC = fESD->GetESDZDC();
  if (!esdZDC) {
    return;
  }

  double zc = -1.0;
  double za = -1.0;
  double zn = -1.0;
  fza = -1;
  fzc = -1;
  fzn = -1;
  zc = esdZDC->GetZDCN1Energy();
  za = esdZDC->GetZDCN2Energy();
  zc = zc * 0.001;
  za = za * 0.001;
  zn = zc + za;
  fza = za;
  fzc = zc;
  fzn = zn;
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
  fTracklets14 = 0;
  fTracklets10 = 0;
  fTrackletsEtaGap = 0;
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
  if (TMath::Abs(spdVtxZ) > 10.0) {
    return;
  }

  nTracklets = SPDptr->GetNumberOfTracklets();

  for (auto it = 0; it < nTracklets; it++) {
    double eta = SPDptr->GetEta(it);
    double phi = SPDptr->GetPhi(it);

    if (TMath::Abs(eta) <= 1.0) {
      fTracklets10++;
    }

    if (fEtaGapNchMin <= TMath::Abs(eta) && TMath::Abs(eta) <= fEtaGapNchMax) {
      fTrackletsEtaGap++;
      hEtaGapSPD->Fill(eta);
    }

    if (TMath::Abs(eta) <= 1.4) {
      hPhiEtaSPD->Fill(phi, eta);
      fTracklets14++;
    }
  }

  hVtxZvsTracklets->Fill(spdVtxZ, fTracklets14);
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSound::MultiplicityDistributions() {
  int rec_nch{0};
  int rec_nch_neg_eta{0};
  int rec_nch_pos_eta{0};
  fTracksEtaGapTPC = 0;
  const int n_tracks{fESD->GetNumberOfTracks()};

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
    if (fEtaMin <= track->Eta() && track->Eta() < 0.0) {
      rec_nch_neg_eta++;
    }
    if (track->Eta() >= 0.0 && track->Eta() <= fEtaMax) {
      rec_nch_pos_eta++;
    }
    if (fEtaGapTPCNchMin <= TMath::Abs(track->Eta()) &&
        TMath::Abs(track->Eta()) <= fEtaGapTPCNchMax) {
      fTracksEtaGapTPC++;
      hPhiEtaGapTPC->Fill(track->Phi(), track->Eta());
    }
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
    //! pT Spectra with NEGATIVE eta
    if (fEtaMin <= track->Eta() && track->Eta() < 0.0) {
      hPtEtaNegvsNchEtaPos->Fill(rec_nch_pos_eta, pt);
      pPtEtaNegvsNchEtaPos->Fill(rec_nch_pos_eta, pt);
    }
    //! pT Spectra with POSITIVE eta
    if (track->Eta() >= 0.0 && track->Eta() <= fEtaMax) {
      hPtEtaPosvsNchEtaNeg->Fill(rec_nch_neg_eta, pt);
      pPtEtaPosvsNchEtaNeg->Fill(rec_nch_neg_eta, pt);
    }
    pPtvsNch->Fill(rec_nch, pt);
    pPtvsV0MAmp->Fill(fv0mamplitude, pt);
    hPtvsV0MAmp->Fill(fv0mamplitude, pt);
    if (fTracklets10 > 0) {
      hPtvsTracklets10->Fill(fTracklets10, pt);
      pPtvsTracklets10->Fill(fTracklets10, pt);
    }
    if (fTracklets14 > 0) {
      hPtvsTracklets14->Fill(fTracklets14, pt);
      pPtvsTracklets14->Fill(fTracklets14, pt);
    }
    if (fTrackletsEtaGap > 0) {
      if (TMath::Abs(track->Eta()) <= fEtaGappT) {
        hPtvsTrackletsEtaGap->Fill(fTrackletsEtaGap, pt);
        pPtvsTrackletsEtaGap->Fill(fTrackletsEtaGap, pt);
      }
    }
    if (fTracksEtaGapTPC > 0) {
      if (TMath::Abs(track->Eta()) <= fEtaGapTPCpT) {
        hPtvsTracksEtaGapTPC->Fill(fTracksEtaGapTPC, pt);
        pPtvsTracksEtaGapTPC->Fill(fTracksEtaGapTPC, pt);
      }
    }
  }

  hV0Mmult->Fill(fv0mpercentile);
  hV0MAmplitude->Fill(fv0mamplitude);
  hNchvsV0M->Fill(rec_nch, fv0mpercentile);
  hNchvsV0MAmp->Fill(rec_nch, fv0mamplitude);
  hV0MvsV0MAmp->Fill(fv0mamplitude, fv0mpercentile);
  hNchEtaPosvsNchEtaNeg->Fill(rec_nch_pos_eta, rec_nch_neg_eta);
  if (fTracklets10 > 0) {
    hTrackletsvsV0MAmp10->Fill(fTracklets10, fv0mamplitude);
  }
  if (fTracklets14 > 0) {
    hTrackletsvsV0MAmp14->Fill(fTracklets14, fv0mamplitude);
  }
  if (fTrackletsEtaGap > 0) {
    hTrackletsEtaGap->Fill(fTrackletsEtaGap);
  }
  if (fTracksEtaGapTPC > 0) {
    hTracksEtaGapTPC->Fill(fTracksEtaGapTPC);
  }
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
