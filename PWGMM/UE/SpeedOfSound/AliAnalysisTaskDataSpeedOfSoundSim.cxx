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

static const int nSpecies{8};
static const int nParticleSource{3};
static const char* Species[nSpecies] = {"Ch",     "Pion",   "Kaon",   "Proton",
                                        "Sigmap", "Sigmam", "Lambda", "Rest"};
static const char* ParticleSource[nParticleSource] = {"Pri", "WeDe", "MaIn"};

#include "AliAnalysisTaskDataSpeedOfSoundSim.h"

class AliAnalysisTaskDataSpeedOfSoundSim;  // your analysis class

ClassImp(AliAnalysisTaskDataSpeedOfSoundSim)  // classimp: necessary for root

    AliAnalysisTaskDataSpeedOfSoundSim::AliAnalysisTaskDataSpeedOfSoundSim()
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
      fEtaMin(-0.8),
      fEtaMax(0.8),
      fEtaGapSPDpT(0.4),
      fEtaGapSPDNchMin(0.7),
      fEtaGapSPDNchMax(1.4),
      fEtaGapTPCpT(0.3),
      fPtMinCent(0.15),
      fPtMaxCent(10.0),
      fEtaGapTPCNchMin(0.5),
      fEtaGapTPCNchMax(0.8),
      fPtMin(0.15),
      fV0Mmin(0.0),
      fV0Mmax(80.0),
      fHMCut(10.0),
      fRandomNumberCut(0.5),
      fv0mpercentile(0),
      fv0mamplitude(0),
      fRecNch(0),
      fTrueNch(0),
      fTrueNch14(0),
      fTrueNch10(0),
      fTrueNchEtaGap(0),
      fTrueNchEtaGapTPC(0),
      fTrueNchEtaPos(0),
      fTrueNchEtaNeg(0),
      fTrueV0(0),
      fTracklets14(0),
      fTracklets10(0),
      fTrackletsEtaGap(0),
      fTracksEtaGapTPC(0),
      fTracksEtaNeg(0),
      fTracksEtaPos(0),
      fTracksEta08(0),
      fdcaxy(-999),
      fdcaz(-999),
      fMultSelection(0x0),
      hBestVtxZ(0),
      hNchvsV0M(0),
      hNchvsV0MAmp(0),
      hV0MvsV0MAmp(0),
      pV0MAmpChannel(0),
      hV0MAmplitude(0),
      hCounter(0),
      hV0Mmult(0),
      hPtvsV0MAmp(0),
      hNchEtaPosvsNchEtaNeg(0),
      hPtEtaNegvsNchEtaPos(0),
      hPtEtaPosvsNchEtaNeg(0),
      hTrueVtxZ(0),
      hTrueNch(0),
      hTrueV0MAmp(0),
      hTruePtvsTrueNch(0),
      hTrueNchEtaPosvsTrueNchEtaNeg(0),
      hTruePtEtaNegvsTrueNchEtaPos(0),
      hTruePtEtaPosvsTrueNchEtaNeg(0),
      hTrueNchEtaGap(0),
      hTrueNchEtaGapTPC(0),
      hTruePtvsTrueNchEtaGap(0),
      hTruePtvsTrueNchEtaGapTPC(0),
      hTrueNchvsMeasNchEtaGapSPD(0),
      hTrueNchvsMeasNchEtaGapTPC(0),
      hTrueNchvsMeasNchHalfTPC(0),
      hPtOutAll_ch(0),
      hPhiEta(0),
      hPhiEtaGap_SPD(0),
      hPhiEtaGap_TPC(0),
      hTrackletsEtaGap(0),
      hTracksEtaGapTPC(0),
      hPtvsTrackletsEtaGap(0),
      hPtvsTracksEtaGapTPC(0) {
  for (int i = 0; i < nParticleSource; ++i) {
    hDCAxy[i] = 0;
  }
  for (int i = 0; i < nSpecies; ++i) {
    hPtInPrim[i] = 0;
    hPtOutPrim[i] = 0;
    hPtInPrim_EtaGap_SPD[i] = 0;
    hPtOutPrim_EtaGap_SPD[i] = 0;
    hPtInPrim_EtaGap_TPC[i] = 0;
    hPtOutPrim_EtaGap_TPC[i] = 0;
    hPtInPrim_HalfEta[i] = 0;
    hPtOutPrim_HalfEta[i] = 0;
  }
}
//_____________________________________________________________________________
AliAnalysisTaskDataSpeedOfSoundSim::AliAnalysisTaskDataSpeedOfSoundSim(
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
      fEtaMin(-0.8),
      fEtaMax(0.8),
      fEtaGapSPDpT(0.4),
      fEtaGapSPDNchMin(0.7),
      fEtaGapSPDNchMax(1.4),
      fEtaGapTPCpT(0.3),
      fPtMinCent(0.15),
      fPtMaxCent(10.0),
      fEtaGapTPCNchMin(0.5),
      fEtaGapTPCNchMax(0.8),
      fPtMin(0.15),
      fV0Mmin(0.0),
      fV0Mmax(80.0),
      fHMCut(10.0),
      fRandomNumberCut(0.5),
      fv0mpercentile(0),
      fv0mamplitude(0),
      fRecNch(0),
      fTrueNch(0),
      fTrueNch14(0),
      fTrueNch10(0),
      fTrueNchEtaGap(0),
      fTrueNchEtaGapTPC(0),
      fTrueNchEtaPos(0),
      fTrueNchEtaNeg(0),
      fTrueV0(0),
      fTracklets14(0),
      fTracklets10(0),
      fTrackletsEtaGap(0),
      fTracksEtaGapTPC(0),
      fTracksEtaNeg(0),
      fTracksEtaPos(0),
      fTracksEta08(0),
      fdcaxy(-999),
      fdcaz(-999),
      fMultSelection(0x0),
      hBestVtxZ(0),
      hNchvsV0M(0),
      hNchvsV0MAmp(0),
      hV0MvsV0MAmp(0),
      pV0MAmpChannel(0),
      hV0MAmplitude(0),
      hCounter(0),
      hV0Mmult(0),
      hPtvsV0MAmp(0),
      hNchEtaPosvsNchEtaNeg(0),
      hPtEtaNegvsNchEtaPos(0),
      hPtEtaPosvsNchEtaNeg(0),
      hTrueVtxZ(0),
      hTrueNch(0),
      hTrueV0MAmp(0),
      hTruePtvsTrueNch(0),
      hTrueNchEtaPosvsTrueNchEtaNeg(0),
      hTruePtEtaNegvsTrueNchEtaPos(0),
      hTruePtEtaPosvsTrueNchEtaNeg(0),
      hTrueNchEtaGap(0),
      hTrueNchEtaGapTPC(0),
      hTruePtvsTrueNchEtaGap(0),
      hTruePtvsTrueNchEtaGapTPC(0),
      hTrueNchvsMeasNchEtaGapSPD(0),
      hTrueNchvsMeasNchEtaGapTPC(0),
      hTrueNchvsMeasNchHalfTPC(0),
      hPtOutAll_ch(0),
      hPhiEta(0),
      hPhiEtaGap_SPD(0),
      hPhiEtaGap_TPC(0),
      hTrackletsEtaGap(0),
      hTracksEtaGapTPC(0),
      hPtvsTrackletsEtaGap(0),
      hPtvsTracksEtaGapTPC(0) {
  for (int i = 0; i < nParticleSource; ++i) {
    hDCAxy[i] = 0;
  }
  for (int i = 0; i < nSpecies; ++i) {
    hPtInPrim[i] = 0;
    hPtOutPrim[i] = 0;
    hPtInPrim_EtaGap_SPD[i] = 0;
    hPtOutPrim_EtaGap_SPD[i] = 0;
    hPtInPrim_EtaGap_TPC[i] = 0;
    hPtOutPrim_EtaGap_TPC[i] = 0;
    hPtInPrim_HalfEta[i] = 0;
    hPtOutPrim_HalfEta[i] = 0;
  }
  DefineInput(0, TChain::Class());  // define the input of the analysis: in this
                                    // case you take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about
  // it, does its work automatically
  DefineOutput(1, TList::Class());  // define the ouptut of the analysis: in
                                    // this case it's a list of histograms
}
//_____________________________________________________________________________
AliAnalysisTaskDataSpeedOfSoundSim::~AliAnalysisTaskDataSpeedOfSoundSim() {
  // destructor
  if (fOutputList) {
    delete fOutputList;  // at the end of your task, it is deleted from memory
                         // by calling this function
    fOutputList = 0x0;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSoundSim::UserCreateOutputObjects() {
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

  constexpr int pt_Nbins{210};
  double pt_bins[pt_Nbins + 1] = {0};
  for (int i = 0; i <= pt_Nbins; ++i) {
    pt_bins[i] = 0.15 + (i * 0.05);
  }

  // Nch (|eta|<0.8) up to 2500
  constexpr double nch_width{2.0};
  constexpr int nch_Nbins{2500};
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

  constexpr int v0mAmp_Nbins_true{800};
  double v0mAmp_bins_true[v0mAmp_Nbins_true + 1] = {0};
  for (int i = 0; i <= v0mAmp_Nbins_true; ++i) {
    v0mAmp_bins_true[i] = 0.0 + i * v0mAmp_width;
  }

  constexpr int dcaxy_Nbins{100};
  double dcaxy_bins[dcaxy_Nbins + 1] = {0};
  for (int i = 0; i <= dcaxy_Nbins; ++i) {
    dcaxy_bins[i] = -3.0 + (0.06 * i);
  }

  constexpr int v0m_Nbins080{6};
  constexpr double v0m_bins080[v0m_Nbins080 + 1] = {0.0,  1.0,  5.0, 10.0,
                                                    20.0, 50.0, 80.0};

  hCounter = new TH1F("hCounter", ";Counts;", 2, 0, 1);
  hCounter->GetXaxis()->SetBinLabel(1, "Spectra");
  hCounter->GetXaxis()->SetBinLabel(2, "Corrections");
  fOutputList->Add(hCounter);

  hBestVtxZ =
      new TH1F("hBestVtxZ", ";Vertex_{#it{z}} (cm); Counts;", 400, -11, 11);
  fOutputList->Add(hBestVtxZ);

  hV0Mmult =
      new TH1F("hV0Mmult", ";V0M (%);Entries", v0m_Nbins080, v0m_bins080);
  fOutputList->Add(hV0Mmult);

  hNchvsV0M = new TH2D("hNchvsV0M", ";#it{N}_{ch}; V0M (%)", nch_Nbins,
                       nch_bins, v0m_Nbins080, v0m_bins080);

  hNchvsV0MAmp = new TH2D("hNchvsV0MAmp", ";#it{N}_{ch}; V0M Amp", nch_Nbins,
                          nch_bins, v0mAmp_Nbins, v0mAmp_bins);

  hV0MvsV0MAmp = new TH2D("hV0MvsV0MAmp", ";V0M Ampl; V0M (%)", v0mAmp_Nbins,
                          v0mAmp_bins, v0m_Nbins080, v0m_bins080);

  pV0MAmpChannel =
      new TProfile("pV0MAmpChannel", ";V0 Channel; Amplitude;", 64, -0.5, 63.5);

  hV0MAmplitude =
      new TH1D("hV0MAmp", ";V0M Amplitude; Entries", v0mAmp_Nbins, v0mAmp_bins);

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
  hPhiEta = new TH2F("hPhiEta", ";#varphi; #eta (|#eta|#leq0.8)", 80, 0,
                     2 * TMath::Pi(), 80, -1.4, 1.4);
  hPhiEtaGap_SPD =
      new TH2F("hPhiEtaGap_SPD", ";#varphi; #eta (0.7#leq|#eta|#leq1.4)", 80, 0,
               2 * TMath::Pi(), 80, -1.4, 1.4);
  hPhiEtaGap_TPC =
      new TH2F("hPhiEtaGap_TPC", ";#varphi; #eta (0.5#leq|#eta|#leq0.8)", 80, 0,
               2 * TMath::Pi(), 80, -1.4, 1.4);

  hTrackletsEtaGap = new TH1F(
      "hTrackletsEtaGap", ";#it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); Counts",
      tracklets_Nbins, tracklets_bins);

  hTracksEtaGapTPC = new TH1F("hTracksEtaGapTPC",
                              ";#it{N}_{ch} (0.5#leq|#eta|#leq0.8); Counts",
                              nch_Nbins, nch_bins);

  hPtvsTrackletsEtaGap =
      new TH2D("hPtvsTrackletsEtaGap",
               "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); #it{p}_{T} "
               "(|#eta|<0.4, GeV/#it{c})",
               tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);

  hPtvsTracksEtaGapTPC = new TH2D("hPtvsTracksEtaGapTPC",
                                  "; #it{N}_{ch} (0.5#leq|#eta|#leq0.8); "
                                  "#it{p}_{T} (|#eta|<0.3, GeV/#it{c})",
                                  nch_Nbins, nch_bins, pt_Nbins, pt_bins);

  fOutputList->Add(hNchvsV0M);
  fOutputList->Add(hNchvsV0MAmp);
  fOutputList->Add(hV0MvsV0MAmp);
  fOutputList->Add(pV0MAmpChannel);
  fOutputList->Add(hV0MAmplitude);
  fOutputList->Add(hPtvsV0MAmp);
  fOutputList->Add(hNchEtaPosvsNchEtaNeg);
  fOutputList->Add(hPtEtaNegvsNchEtaPos);
  fOutputList->Add(hPtEtaPosvsNchEtaNeg);
  fOutputList->Add(hPhiEta);
  fOutputList->Add(hPhiEtaGap_SPD);
  fOutputList->Add(hPhiEtaGap_TPC);
  fOutputList->Add(hTrackletsEtaGap);
  fOutputList->Add(hTracksEtaGapTPC);
  fOutputList->Add(hPtvsTrackletsEtaGap);
  fOutputList->Add(hPtvsTracksEtaGapTPC);

  hTrueVtxZ =
      new TH1F("hTrueVtxZ", ";z-vertex position;Entries", 200, -10.0, 10.0);

  hTrueNch = new TH1F("hTrueNch", "; #it{N}_{ch}^{true} (|#eta|<0.8); Entries",
                      nch_Nbins, nch_bins);

  hTrueV0MAmp = new TH1F("hTrueV0MAmp", "; V0M Amplitude; Entries",
                         v0mAmp_Nbins_true, v0mAmp_bins_true);

  hTrueNchEtaGap = new TH1F("hTrueNchEtaGap",
                            "; #it{N}_{ch} (0.7#leq|#eta|#leq1.4); Entries ",
                            tracklets_Nbins, tracklets_bins);

  hTrueNchEtaGapTPC = new TH1F("hTrueNchEtaGapTPC",
                               "; #it{N}_{ch} (0.5#leq|#eta|#leq0.8); Entries ",
                               nch_Nbins, nch_bins);

  hTruePtvsTrueNch =
      new TH2D("hTruePtvsTrueNch", ";#it{N}_{ch}^{true}; #it{p}_{T} GeV/#it{c}",
               nch_Nbins, nch_bins, pt_Nbins, pt_bins);

  hTrueNchEtaPosvsTrueNchEtaNeg =
      new TH2D("hTrueNchEtaPosvsTrueNchEtaNeg",
               ";#it{N}_{ch} (0#leq#eta#leq0.8); #it{N}_{ch} (-0.8#leq#eta<0)",
               nch_Nbins, nch_bins, nch_Nbins, nch_bins);

  hTruePtEtaNegvsTrueNchEtaPos =
      new TH2D("hTruePtEtaNegvsTrueNchEtaPos",
               "; #it{N}_{ch} (0#leq#eta#leq0.8); #it{p}_{T} "
               "(GeV/#it{c}) (-0.8#leq#eta<0)",
               nch_Nbins, nch_bins, pt_Nbins, pt_bins);

  hTruePtEtaPosvsTrueNchEtaNeg =
      new TH2D("hTruePtEtaPosvsTrueNchEtaNeg",
               "; #it{N}_{ch} (-0.8#leq#eta<0); #it{p}_{T} "
               "(GeV/#it{c}) (0#geq#eta#leq0.8)",
               nch_Nbins, nch_bins, pt_Nbins, pt_bins);

  hPtOutAll_ch = new TH1F("hPtOutAll_Ch", ";#it{p}_{T} (GeV/#it{c}); Entries",
                          pt_Nbins, pt_bins);

  hTruePtvsTrueNchEtaGap =
      new TH2D("hTruePtvsTrueNchEtaGap",
               "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); #it{p}_{T} "
               "(|#eta|<0.4, GeV/#it{c})",
               tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);

  hTruePtvsTrueNchEtaGapTPC =
      new TH2D("hTruePtvsTrueNchEtaGapTPC",
               "; #it{N}_{ch} (0.5#leq|#eta|#leq0.8); #it{p}_{T} (|#eta|<0.3, "
               "GeV/#it{c})",
               nch_Nbins, nch_bins, pt_Nbins, pt_bins);

  hTrueNchvsMeasNchEtaGapSPD = new TH2D(
      "hTrueNchvsMeasNchEtaGapSPD",
      "; #it{N}_{ch}^{true} (0.7#leq|#eta|#leq1.4); #it{N}_{ch}^{rec} "
      "(0.7#leq|#eta|#leq1.4)",
      tracklets_Nbins, tracklets_bins, tracklets_Nbins, tracklets_bins);

  hTrueNchvsMeasNchEtaGapTPC =
      new TH2D("hTrueNchvsMeasNchEtaGapTPC",
               "; #it{N}_{ch}^{true} (0.5#leq|#eta|#leq0.8); #it{N}_{ch}^{rec} "
               "(0.5#leq|#eta|#leq0.8)",
               nch_Nbins, nch_bins, nch_Nbins, nch_bins);

  hTrueNchvsMeasNchHalfTPC =
      new TH2D("hTrueNchvsMeasNchHalfTPC",
               "; #it{N}_{ch}^{true} (0#leq#eta#leq0.8); #it{N}_{ch}^{rec} "
               "(0#leq#eta#leq0.8)",
               nch_Nbins, nch_bins, nch_Nbins, nch_bins);

  fOutputList->Add(hTrueVtxZ);
  fOutputList->Add(hTrueNch);
  fOutputList->Add(hTrueV0MAmp);
  fOutputList->Add(hTruePtvsTrueNch);
  fOutputList->Add(hTrueNchEtaPosvsTrueNchEtaNeg);
  fOutputList->Add(hTruePtEtaNegvsTrueNchEtaPos);
  fOutputList->Add(hTruePtEtaPosvsTrueNchEtaNeg);
  fOutputList->Add(hTruePtvsTrueNchEtaGap);
  fOutputList->Add(hTruePtvsTrueNchEtaGapTPC);
  fOutputList->Add(hTrueNchEtaGap);
  fOutputList->Add(hTrueNchEtaGapTPC);
  fOutputList->Add(hTrueNchvsMeasNchEtaGapSPD);
  fOutputList->Add(hTrueNchvsMeasNchEtaGapTPC);
  fOutputList->Add(hTrueNchvsMeasNchHalfTPC);
  fOutputList->Add(hPtOutAll_ch);

  for (int i = 0; i < nSpecies; ++i) {
    hPtInPrim[i] = new TH1F(Form("hPtInPrim_%s", Species[i]),
                            "With 0#leqV0#leq5 cut (official "
                            "framework);#it{p}_{T} (GeV/#it{c}); Entries",
                            pt_Nbins, pt_bins);
    hPtInPrim_EtaGap_SPD[i] =
        new TH2D(Form("hPtInPrim_EtaGap_SPD_%s", Species[i]),
                 ";#it{N}_{ch}^{true} (0.7#leq|#eta|#leq1.4);#it{p}_{T} "
                 "(|#eta|#leq0.4, GeV/#it{c});",
                 tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);
    hPtInPrim_EtaGap_TPC[i] =
        new TH2D(Form("hPtInPrim_EtaGap_TPC_%s", Species[i]),
                 ";#it{N}_{ch}^{true} (0.5#leq|#eta|#leq0.8);#it{p}_{T} "
                 "(|#eta|#leq0.3, GeV/#it{c});",
                 nch_Nbins, nch_bins, pt_Nbins, pt_bins);
    hPtInPrim_HalfEta[i] =
        new TH2D(Form("hPtInPrim_HalfEta_%s", Species[i]),
                 ";#it{N}_{ch}^{true} (0#leq#eta#leq0.8);#it{p}_{T} "
                 "(-0.8#leq#eta<0, GeV/#it{c});",
                 nch_Nbins, nch_bins, pt_Nbins, pt_bins);
    hPtOutPrim[i] = new TH1F(Form("hPtOutPrim_%s", Species[i]),
                             "With 0#leqV0#leq5 cut (official framework);"
                             "#it{p}_{T} (|#eta|<0.8, GeV/#it{c}); Entries",
                             pt_Nbins, pt_bins);
    hPtOutPrim_EtaGap_SPD[i] =
        new TH2D(Form("hPtOutPrim_EtaGap_SPD_%s", Species[i]),
                 ";#it{N}_{ch}^{rec} (0.7#leq|#eta|#leq1.4);#it{p}_{T} "
                 "(|#eta|#leq0.4, GeV/#it{c});",
                 tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);
    hPtOutPrim_EtaGap_TPC[i] =
        new TH2D(Form("hPtOutPrim_EtaGap_TPC_%s", Species[i]),
                 ";#it{N}_{ch}^{rec} (0.5#leq|#eta|#leq0.8);#it{p}_{T} "
                 "(|#eta|#leq0.3, GeV/#it{c});",
                 nch_Nbins, nch_bins, pt_Nbins, pt_bins);
    hPtOutPrim_HalfEta[i] =
        new TH2D(Form("hPtOutPrim_HalfEta_%s", Species[i]),
                 ";#it{N}_{ch}^{rec} (0#leq#eta#leq0.8);#it{p}_{T} "
                 "(-0.8#leq#eta<0, GeV/#it{c});",
                 nch_Nbins, nch_bins, pt_Nbins, pt_bins);

    fOutputList->Add(hPtInPrim[i]);
    fOutputList->Add(hPtOutPrim[i]);
    fOutputList->Add(hPtInPrim_EtaGap_SPD[i]);
    fOutputList->Add(hPtOutPrim_EtaGap_SPD[i]);
    fOutputList->Add(hPtInPrim_EtaGap_TPC[i]);
    fOutputList->Add(hPtOutPrim_EtaGap_TPC[i]);
    fOutputList->Add(hPtInPrim_HalfEta[i]);
    fOutputList->Add(hPtOutPrim_HalfEta[i]);
  }

  for (int i = 0; i < nParticleSource; ++i) {
    hDCAxy[i] = new TH2F(Form("hDCAxy_%s", ParticleSource[i]),
                         "With 0#leqV0#leq5 cut (official framework);DCA_{xy} "
                         "(cm);#it{p}_{T} (|#eta|<0.8, GeV/#it{c})",
                         dcaxy_Nbins, dcaxy_bins, pt_Nbins, pt_bins);
    fOutputList->Add(hDCAxy[i]);
  }

  fEventCuts.AddQAplotsToList(fOutputList);
  PostData(1, fOutputList);  // postdata will notify the analysis manager of
                             // changes / updates to the
}
//_____________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSoundSim::UserExec(Option_t*) {
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

  double random_number = -1.0;
  bool fill_corrections{false};
  gRandom->SetSeed(0);
  random_number = gRandom->Uniform(0.0, 1.0);
  // if random_number < 0.5 --> Multiplicity Distributions
  // if random_number >= 0.5 --> Detector Response & Corrections
  if (random_number >= fRandomNumberCut) {
    fill_corrections = true;
  }

  fv0mpercentile = -999.0;
  fv0mamplitude = -999.0;

  fMultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
  fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  // ftrackmult08 = AliESDtrackCuts::GetReferenceMultiplicity(
  //     fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);

  //! Analyze only the 0--80 % V0M range
  if (!(fv0mpercentile >= fV0Mmin && fv0mpercentile < fV0Mmax)) {
    return;
  }

  bool isGoodVtxPosMC{false};
  AliHeader* headerMC = fMC->Header();
  AliGenEventHeader* genHeader = headerMC->GenEventHeader();
  TArrayF vtxMC(3);  // primary vertex  MC
  vtxMC[0] = 9999;
  vtxMC[1] = 9999;
  vtxMC[2] = 9999;  // initialize with dummy
  if (genHeader) {
    genHeader->PrimaryVertex(vtxMC);
  }
  if (TMath::Abs(vtxMC[2]) <= 10.0) {
    isGoodVtxPosMC = true;
  }

  float vtx_z{999};
  vtx_z = vtxMC[2];

  if (!isGoodVtxPosMC) {
    return;
  }

  ReadMCEvent();

  if (!(fTrueNch10 > 0 && fTrueV0 > 0)) {
    return;
  }

  hCounter->Fill(random_number);
  hTrueVtxZ->Fill(vtx_z);

  if (!fill_corrections) {
    TrueMultiplicityDistributions();
  }

  bool isEventTriggered{false};
  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  isEventTriggered = fSelectMask & fTrigger;
  if (!isEventTriggered) {
    return;
  }

  if (!fEventCuts.AcceptEvent(event)) {
    PostData(1, fOutputList);
    return;
  }

  bool hasRecVertex = false;
  hasRecVertex = HasRecVertex();
  if (!hasRecVertex) {
    return;
  }

  VertexPosition();
  GetCalibratedV0Amplitude();
  GetSPDMultiplicity();
  GetTPCMultiplicity();

  if (isGoodVtxPosMC) {
    if (!fill_corrections) {
      MultiplicityDistributions();
    } else {
      TrackingEfficiency();
      DCA();
    }
  }

  PostData(1, fOutputList);
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::Terminate(Option_t*) {}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::VertexPosition() {
  //! best primary vertex available
  const AliVVertex* vtx = fEventCuts.GetPrimaryVertex();
  hBestVtxZ->Fill(vtx->GetZ());
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::GetCalibratedV0Amplitude() {
  float mV0M{0.0};
  for (int i = 0; i < 64; i++) {
    mV0M += fESD->GetVZEROEqMultiplicity(i);
    pV0MAmpChannel->Fill(i, fESD->GetVZEROEqMultiplicity(i));
  }
  fv0mamplitude = mV0M;
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::GetSPDMultiplicity() {
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
    if (TMath::Abs(eta) <= 1.4) {
      fTracklets14++;
    }
    if (TMath::Abs(eta) >= fEtaGapSPDNchMin &&
        TMath::Abs(eta) <= fEtaGapSPDNchMax) {
      hPhiEtaGap_SPD->Fill(phi, eta);
      fTrackletsEtaGap++;
    }
  }
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::GetTPCMultiplicity() {
  fTracksEtaGapTPC = 0;
  fTracksEtaNeg = 0;
  fTracksEtaPos = 0;
  fTracksEta08 = 0;
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
    fTracksEta08++;
    hPhiEta->Fill(track->Phi(), track->Eta());
    if (track->Eta() >= fEtaMin && track->Eta() < 0.0) {
      fTracksEtaNeg++;
    }
    if (track->Eta() >= 0.0 && track->Eta() <= fEtaMax) {
      fTracksEtaPos++;
    }
    if (TMath::Abs(track->Eta()) >= fEtaGapTPCNchMin &&
        TMath::Abs(track->Eta()) <= fEtaGapTPCNchMax) {
      hPhiEtaGap_TPC->Fill(track->Phi(), track->Eta());
      fTracksEtaGapTPC++;
    }
  }
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::MultiplicityDistributions() {
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
    double pt = track->Pt();
    //! pT Spectra with NEGATIVE eta
    if (track->Eta() >= fEtaMin && track->Eta() < 0.0) {
      hPtEtaNegvsNchEtaPos->Fill(fTracksEtaPos, pt);
    }
    //! pT Spectra with POSITIVE eta
    if (track->Eta() >= 0.0 && track->Eta() <= fEtaMax) {
      hPtEtaPosvsNchEtaNeg->Fill(fTracksEtaNeg, pt);
    }
    //! pT Spectra with SPD eta gap
    if (TMath::Abs(track->Eta()) <= fEtaGapSPDpT) {
      hPtvsTrackletsEtaGap->Fill(fTrackletsEtaGap, pt);
    }
    //! pT Spectra with TPC eta gap
    if (TMath::Abs(track->Eta()) <= fEtaGapTPCpT) {
      hPtvsTracksEtaGapTPC->Fill(fTracksEtaGapTPC, pt);
    }
    hPtvsV0MAmp->Fill(fv0mamplitude, pt);
  }

  hV0Mmult->Fill(fv0mpercentile);
  hV0MAmplitude->Fill(fv0mamplitude);
  hNchvsV0M->Fill(fTracksEta08, fv0mpercentile);
  hNchvsV0MAmp->Fill(fTracksEta08, fv0mamplitude);
  hV0MvsV0MAmp->Fill(fv0mamplitude, fv0mpercentile);
  hNchEtaPosvsNchEtaNeg->Fill(fTracksEtaPos, fTracksEtaNeg);
  hTrackletsEtaGap->Fill(fTrackletsEtaGap);
  hTracksEtaGapTPC->Fill(fTracksEtaGapTPC);
}

//____________________________________________________________
void AliAnalysisTaskDataSpeedOfSoundSim::DCA() {
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

    int label = -1;
    label = TMath::Abs(track->GetLabel());

    TParticle* particle = fMC->GetTrack(label)->Particle();
    if (!particle) {
      continue;
    }
    if (track->Charge() == 0) {
      continue;
    }

    float dcaxy = -999.0;
    float dcaz = -999.0;
    track->GetImpactParameters(dcaxy, dcaz);

    if (fv0mpercentile >= 0.0 && fv0mpercentile <= 5.0) {
      if (TMath::Abs(track->Eta()) <= fEtaCut) {
        if (fMC->IsPhysicalPrimary(label)) {
          hDCAxy[0]->Fill(dcaxy, track->Pt());
        } else if (fMC->IsSecondaryFromWeakDecay(label)) {
          hDCAxy[1]->Fill(dcaxy, track->Pt());
        } else if (fMC->IsSecondaryFromMaterial(label)) {
          hDCAxy[2]->Fill(dcaxy, track->Pt());
        } else {
          continue;
        }
      }
    }
  }
}
//____________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::TrackingEfficiency() {
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
    if (TMath::Abs(track->Eta()) > fEtaCut) {
      continue;
    }

    int label = -1;
    label = TMath::Abs(track->GetLabel());

    TParticle* particle = fMC->GetTrack(label)->Particle();
    if (!particle) {
      continue;
    }
    if (track->Charge() == 0) {
      continue;
    }

    if (fv0mpercentile >= 0.0 && fv0mpercentile <= 5.0) {
      hPtOutAll_ch->Fill(track->Pt());
    }

    if (!fMC->IsPhysicalPrimary(label)) {
      continue;
    }

    int partPDG = particle->GetPdgCode();
    int pidCode = GetPidCode(partPDG);

    if (fv0mpercentile >= 0.0 && fv0mpercentile <= 5.0) {
      if (TMath::Abs(track->Eta()) <= fEtaCut) {
        hPtOutPrim[0]->Fill(track->Pt());
        hPtOutPrim[pidCode]->Fill(track->Pt());
      }
    }

    if (TMath::Abs(track->Eta()) <= fEtaGapSPDpT) {
      hPtOutPrim_EtaGap_SPD[0]->Fill(fTrackletsEtaGap, track->Pt());
      hPtOutPrim_EtaGap_SPD[pidCode]->Fill(fTrackletsEtaGap, track->Pt());
    }

    if (TMath::Abs(track->Eta()) <= fEtaGapTPCpT) {
      hPtOutPrim_EtaGap_TPC[0]->Fill(fTracksEtaGapTPC, track->Pt());
      hPtOutPrim_EtaGap_TPC[pidCode]->Fill(fTracksEtaGapTPC, track->Pt());
    }

    if (track->Eta() >= fEtaMin && track->Eta() < 0.0) {
      hPtOutPrim_HalfEta[0]->Fill(fTracksEtaPos, track->Pt());
      hPtOutPrim_HalfEta[pidCode]->Fill(fTracksEtaPos, track->Pt());
    }
  }

  const int n_particles{fMC->GetNumberOfTracks()};
  for (int i = 0; i < n_particles; ++i) {
    AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
    if (!particle) {
      continue;
    }
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC)) {
      continue;
    }
    if (TMath::Abs(particle->Eta()) > fEtaCut) {
      continue;
    }
    if (particle->Pt() < fPtMin) {
      continue;
    }
    if (!fMC->IsPhysicalPrimary(i)) {
      continue;
    }

    int partPDG = particle->PdgCode();
    int pidCode = GetPidCode(partPDG);

    if (particle->Charge() == 0) {
      if (pidCode == 6) {
        if (fv0mpercentile >= 0.0 && fv0mpercentile <= 5.0) {
          if (TMath::Abs(particle->Eta()) <= fEtaCut) {
            hPtInPrim[pidCode]->Fill(particle->Pt());
          }
        }

        if (TMath::Abs(particle->Eta()) <= fEtaGapSPDpT) {
          hPtInPrim_EtaGap_SPD[pidCode]->Fill(fTrueNchEtaGap, particle->Pt());
        }

        if (TMath::Abs(particle->Eta()) <= fEtaGapTPCpT) {
          hPtInPrim_EtaGap_TPC[pidCode]->Fill(fTrueNchEtaGapTPC,
                                              particle->Pt());
        }

        if (particle->Eta() >= fEtaMin && particle->Eta() < 0.0) {
          hPtInPrim_HalfEta[pidCode]->Fill(fTrueNchEtaPos, particle->Pt());
        }
      }
      continue;
    }

    if (fv0mpercentile >= 0.0 && fv0mpercentile <= 5.0) {
      if (TMath::Abs(particle->Eta()) <= fEtaCut) {
        hPtInPrim[0]->Fill(particle->Pt());
        hPtInPrim[pidCode]->Fill(particle->Pt());
      }
    }

    if (TMath::Abs(particle->Eta()) <= fEtaGapSPDpT) {
      hPtInPrim_EtaGap_SPD[0]->Fill(fTrueNchEtaGap, particle->Pt());
      hPtInPrim_EtaGap_SPD[pidCode]->Fill(fTrueNchEtaGap, particle->Pt());
    }

    if (TMath::Abs(particle->Eta()) <= fEtaGapTPCpT) {
      hPtInPrim_EtaGap_TPC[0]->Fill(fTrueNchEtaGapTPC, particle->Pt());
      hPtInPrim_EtaGap_TPC[pidCode]->Fill(fTrueNchEtaGapTPC, particle->Pt());
    }

    if (particle->Eta() >= fEtaMin && particle->Eta() < 0.0) {
      hPtInPrim_HalfEta[0]->Fill(fTrueNchEtaPos, particle->Pt());
      hPtInPrim_HalfEta[pidCode]->Fill(fTrueNchEtaPos, particle->Pt());
    }
  }
  hTrueNchvsMeasNchEtaGapSPD->Fill(fTrueNchEtaGap, fTrackletsEtaGap);
  hTrueNchvsMeasNchEtaGapTPC->Fill(fTrueNchEtaGapTPC, fTracksEtaGapTPC);
  hTrueNchvsMeasNchHalfTPC->Fill(fTrueNchEtaPos, fTracksEtaPos);
}

//____________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::ReadMCEvent() {
  fTrueNch = 0;
  fTrueV0 = 0;
  fTrueNch10 = 0;
  fTrueNch14 = 0;
  fTrueNchEtaPos = 0;
  fTrueNchEtaNeg = 0;
  fTrueNchEtaGap = 0;
  fTrueNchEtaGapTPC = 0;

  const int n_particles{fMC->GetNumberOfTracks()};
  for (int i = 0; i < n_particles; ++i) {
    AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
    if (!particle) {
      continue;
    }
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC)) {
      continue;
    }
    if (particle->Charge() == 0.0) {
      continue;
    }
    if (!fMC->IsPhysicalPrimary(i)) {
      continue;
    }
    if ((2.8 < particle->Eta() && particle->Eta() < 5.1) ||
        (-3.7 < particle->Eta() && particle->Eta() < -1.7)) {
      fTrueV0++;
    }
    if (TMath::Abs(particle->Eta()) < 1.4) {
      fTrueNch14++;
    }
    if (TMath::Abs(particle->Eta()) < 1.0) {
      fTrueNch10++;
    }
    if (TMath::Abs(particle->Eta()) >= fEtaGapSPDNchMin &&
        TMath::Abs(particle->Eta()) <= fEtaGapSPDNchMax) {
      fTrueNchEtaGap++;
    }
    if (particle->Pt() < fPtMinCent || particle->Pt() > fPtMaxCent) {
      continue;
    }
    if (TMath::Abs(particle->Eta()) > fEtaCut) {
      continue;
    }
    if (particle->Eta() >= fEtaMin && particle->Eta() < 0.0) {
      fTrueNchEtaNeg++;
    }
    if (particle->Eta() >= 0.0 && particle->Eta() <= fEtaMax) {
      fTrueNchEtaPos++;
    }
    if (TMath::Abs(particle->Eta()) >= fEtaGapTPCNchMin &&
        TMath::Abs(particle->Eta()) <= fEtaGapTPCNchMax) {
      fTrueNchEtaGapTPC++;
    }
    fTrueNch++;
  }
}

//____________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::TrueMultiplicityDistributions() {
  const int n_particles{fMC->GetNumberOfTracks()};
  for (int i = 0; i < n_particles; ++i) {
    AliMCParticle* particle = (AliMCParticle*)fMC->GetTrack(i);
    if (!particle) {
      continue;
    }
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC)) {
      continue;
    }
    if (particle->Charge() == 0.0) {
      continue;
    }
    if (particle->Pt() < fPtMin) {
      continue;
    }
    if (TMath::Abs(particle->Eta()) > fEtaCut) {
      continue;
    }
    if (!fMC->IsPhysicalPrimary(i)) {
      continue;
    }

    hTruePtvsTrueNch->Fill(fTrueNch, particle->Pt());
    if (TMath::Abs(particle->Eta()) <= fEtaGapSPDpT) {
      hTruePtvsTrueNchEtaGap->Fill(fTrueNchEtaGap, particle->Pt());
    }
    if (TMath::Abs(particle->Eta()) <= fEtaGapTPCpT) {
      hTruePtvsTrueNchEtaGapTPC->Fill(fTrueNchEtaGapTPC, particle->Pt());
    }
    if (particle->Eta() >= fEtaMin && particle->Eta() < 0.0) {
      hTruePtEtaNegvsTrueNchEtaPos->Fill(fTrueNchEtaPos, particle->Pt());
    }
    if (particle->Eta() >= 0.0 && particle->Eta() <= fEtaMax) {
      hTruePtEtaPosvsTrueNchEtaNeg->Fill(fTrueNchEtaNeg, particle->Pt());
    }
  }

  hTrueNch->Fill(fTrueNch);
  hTrueV0MAmp->Fill(fTrueV0);
  hTrueNchEtaGap->Fill(fTrueNchEtaGap);
  hTrueNchEtaGapTPC->Fill(fTrueNchEtaGapTPC);
  hTrueNchEtaPosvsTrueNchEtaNeg->Fill(fTrueNchEtaPos, fTrueNchEtaNeg);
}

//____________________________________________________________

Bool_t AliAnalysisTaskDataSpeedOfSoundSim::HasRecVertex() {
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
                  : 0.;  /// If one of the two vertices is not available this
                         /// cut is always passed.
  double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
  double errTrc =
      bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  /// vertex dispersion for run1, only for ESD, AOD code to be added here
  const AliESDVertex* vtSPDESD = dynamic_cast<const AliESDVertex*>(vtSPD);
  double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
  if ((TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute &&
       nsigTot <= fMaxDeltaSpdTrackNsigmaSPD &&
       nsigTrc <=
           fMaxDeltaSpdTrackNsigmaTrack) &&  // discrepancy track-SPD vertex
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

void AliAnalysisTaskDataSpeedOfSoundSim::ChangeCut(AliESDtrackCuts* fCuts) {
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
//_____________________________________________________________________________
int AliAnalysisTaskDataSpeedOfSoundSim::GetPidCode(int pdgCode) const {
  // return our internal code for pions, kaons, and protons

  int pidCode{1};
  switch (TMath::Abs(pdgCode)) {
    case 211:
      pidCode = 1;  // pion
      break;
    case 321:
      pidCode = 2;  // kaon
      break;
    case 2212:
      pidCode = 3;  // proton
      break;
    case 3222:
      pidCode = 4;  // sigma plus
      break;
    case 3112:
      pidCode = 5;  // sigma minus
      break;
    case 3122:
      pidCode = 6;  // lambda
      break;
    default:
      pidCode = 7;  // something else?
  };

  return pidCode;
}
