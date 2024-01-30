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
static constexpr double uc_v0m_bins_high[v0m_Nbins] = {5.0};
static constexpr double uc_v0m_bins_low[v0m_Nbins] = {0.0};
static const char* uc_v0m_bins_name[v0m_Nbins] = {"0_5"};

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
      fSystematic(1),
      fTrigger(AliVEvent::kCentral),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fOutputList(0),
      fEtaCut(0.8),
      fEtaMin(-0.8),
      fEtaMax(0.8),
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
      fTrueNchEtaPos(0),
      fTrueNchEtaNeg(0),
      fTrueV0(0),
      fTracklets14(0),
      fTracklets10(0),
      fTrackletsEtaGap(0),
      fdcaxy(-999),
      fdcaz(-999),
      fMultSelection(0x0),
      hNchvsV0M(0),
      hNchvsV0MAmp(0),
      hV0MvsV0MAmp(0),
      pV0MAmpChannel(0),
      hV0MAmplitude(0),
      hCounter(0),
      hV0Mmult(0),
      pPtvsNch(0),
      pPtEtaNegvsNchEtaPos(0),
      pPtEtaPosvsNchEtaNeg(0),
      pPtvsV0MAmp(0),
      hPtvsV0MAmp(0),
      hNchEtaPosvsNchEtaNeg(0),
      hPtEtaNegvsNchEtaPos(0),
      hPtEtaPosvsNchEtaNeg(0),
      hTrueVtxZ(0),
      hTrueNch(0),
      hTrueV0MAmp(0),
      hNchResponse(0),
      hTruePtvsTrueNch(0),
      hTrueNchEtaPosvsTrueNchEtaNeg(0),
      hTruePtEtaNegvsTrueNchEtaPos(0),
      hTruePtEtaPosvsTrueNchEtaNeg(0),
      hTrueNchEtaGap(0),
      hTruePtvsTrueNchEtaGap(0),
      hPtInPrim_ch(0),
      hPtInPrim_pion(0),
      hPtInPrim_kaon(0),
      hPtInPrim_proton(0),
      hPtInPrim_sigmap(0),
      hPtInPrim_sigmam(0),
      hPtInPrim_lambda(0),
      hPtInPrim_rest(0),
      hPtOutAll_ch(0),
      hPtOutPrim_ch(0),
      hPtOutPrim_pion(0),
      hPtOutPrim_kaon(0),
      hPtOutPrim_proton(0),
      hPtOutPrim_sigmap(0),
      hPtOutPrim_sigmam(0),
      hPtOutPrim_rest(0),
      hPhiEtaSPD(0),
      hVtxZvsTracklets(0),
      hTrackletsEtaGapvsV0MAmp(0),
      hPtvsTrackletsEtaGap(0) {
  for (int i = 0; i < v0m_Nbins; ++i) {
    hDCAxyPri[i] = 0;
    hDCAxyWeDe[i] = 0;
    hDCAxyMaIn[i] = 0;
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
      fSystematic(1),
      fTrigger(AliVEvent::kCentral),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fOutputList(0),
      fEtaCut(0.8),
      fEtaMin(-0.8),
      fEtaMax(0.8),
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
      fTrueNchEtaPos(0),
      fTrueNchEtaNeg(0),
      fTrueV0(0),
      fTracklets14(0),
      fTracklets10(0),
      fTrackletsEtaGap(0),
      fdcaxy(-999),
      fdcaz(-999),
      fMultSelection(0x0),
      hNchvsV0M(0),
      hNchvsV0MAmp(0),
      hV0MvsV0MAmp(0),
      pV0MAmpChannel(0),
      hV0MAmplitude(0),
      hCounter(0),
      hV0Mmult(0),
      pPtvsNch(0),
      pPtEtaNegvsNchEtaPos(0),
      pPtEtaPosvsNchEtaNeg(0),
      pPtvsV0MAmp(0),
      hPtvsV0MAmp(0),
      hNchEtaPosvsNchEtaNeg(0),
      hPtEtaNegvsNchEtaPos(0),
      hPtEtaPosvsNchEtaNeg(0),
      hTrueVtxZ(0),
      hTrueNch(0),
      hTrueV0MAmp(0),
      hNchResponse(0),
      hTruePtvsTrueNch(0),
      hTrueNchEtaPosvsTrueNchEtaNeg(0),
      hTruePtEtaNegvsTrueNchEtaPos(0),
      hTruePtEtaPosvsTrueNchEtaNeg(0),
      hTrueNchEtaGap(0),
      hTruePtvsTrueNchEtaGap(0),
      hPtInPrim_ch(0),
      hPtInPrim_pion(0),
      hPtInPrim_kaon(0),
      hPtInPrim_proton(0),
      hPtInPrim_sigmap(0),
      hPtInPrim_sigmam(0),
      hPtInPrim_lambda(0),
      hPtInPrim_rest(0),
      hPtOutAll_ch(0),
      hPtOutPrim_ch(0),
      hPtOutPrim_pion(0),
      hPtOutPrim_kaon(0),
      hPtOutPrim_proton(0),
      hPtOutPrim_sigmap(0),
      hPtOutPrim_sigmam(0),
      hPtOutPrim_rest(0),
      hPhiEtaSPD(0),
      hVtxZvsTracklets(0),
      hTrackletsEtaGapvsV0MAmp(0),
      hPtvsTrackletsEtaGap(0) {
  for (int i = 0; i < v0m_Nbins; ++i) {
    hDCAxyPri[i] = 0;
    hDCAxyWeDe[i] = 0;
    hDCAxyMaIn[i] = 0;
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

  hPhiEtaSPD = new TH2D("hPhiEtaSPD", ";#varphi; #eta", 80, 0, 2 * TMath::Pi(),
                        80, -1.4, 1.4);

  hVtxZvsTracklets =
      new TH2D("hVtxZvsTracklets", ";#it{Z}_{vtz} (cm); #it{N}_{tracklets}", 50,
               -10, 10, tracklets_Nbins, tracklets_bins);

  hTrackletsEtaGapvsV0MAmp =
      new TH2D("hTrackletsEtaGapvsV0MAmp",
               "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); V0M Amp",
               tracklets_Nbins, tracklets_bins, v0mAmp_Nbins, v0mAmp_bins);

  hPtvsTrackletsEtaGap = new TH2D(
      "hPtvsTrackletsEtaGap",
      "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); #it{p}_{T} (GeV/#it{c})",
      tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);

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
  fOutputList->Add(hVtxZvsTracklets);
  fOutputList->Add(hTrackletsEtaGapvsV0MAmp);
  fOutputList->Add(hPtvsTrackletsEtaGap);

  hTrueVtxZ =
      new TH1F("hTrueVtxZ", ";z-vertex position;Entries", 200, -10.0, 10.0);

  hNchResponse =
      new TH2D("hNchResponse", ";#it{N}_{ch}^{rec}; #it{N}_{ch}^{true};",
               nch_Nbins, nch_bins, nch_Nbins, nch_bins);

  hTrueNch = new TH1F("hTrueNch", "; #it{N}_{ch}^{true} (|#eta|<0.8); Entries",
                      nch_Nbins, nch_bins);

  hTrueV0MAmp = new TH1F("hTrueV0MAmp", "; V0M Amplitude; Entries",
                         v0mAmp_Nbins_true, v0mAmp_bins_true);

  hTrueNchEtaGap = new TH1F("hTrueNchEtaGap",
                            "; #it{N}_{ch} (0.7#leq|#eta|#leq1.4); Entries ",
                            tracklets_Nbins, tracklets_bins);

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

  hPtOutPrim_ch = new TH1F("hPtOutPrim_Ch", ";#it{p}_{T} (GeV/#it{c}); Entries",
                           pt_Nbins, pt_bins);

  hPtOutPrim_pion =
      new TH1F("hPtOutPrim_Pion", ";#it{p}_{T} (GeV/#it{c}); Entries", pt_Nbins,
               pt_bins);

  hPtOutPrim_kaon =
      new TH1F("hPtOutPrim_Kaon", ";#it{p}_{T} (GeV/#it{c}); Entries", pt_Nbins,
               pt_bins);

  hPtOutPrim_proton =
      new TH1F("hPtOutPrim_Proton", ";#it{p}_{T} (GeV/#it{c}); Entries",
               pt_Nbins, pt_bins);

  hPtOutPrim_sigmap =
      new TH1F("hPtOutPrim_Sigmap", ";#it{p}_{T} (GeV/#it{c}); Entries",
               pt_Nbins, pt_bins);

  hPtOutPrim_sigmam =
      new TH1F("hPtOutPrim_Sigmam", ";#it{p}_{T} (GeV/#it{c}); Entries",
               pt_Nbins, pt_bins);

  hPtOutPrim_rest =
      new TH1F("hPtOutPrim_Rest", ";#it{p}_{T} (GeV/#it{c}); Entries", pt_Nbins,
               pt_bins);

  hPtInPrim_ch = new TH1F("hPtInPrim_Ch", ";#it{p}_{T} (GeV/#it{c}); Entries",
                          pt_Nbins, pt_bins);

  hPtInPrim_pion = new TH1F(
      "hPtInPrim_Pion", ";#it{p}_{T} (GeV/#it{c}); Entries", pt_Nbins, pt_bins);

  hPtInPrim_kaon = new TH1F(
      "hPtInPrim_Kaon", ";#it{p}_{T} (GeV/#it{c}); Entries", pt_Nbins, pt_bins);

  hPtInPrim_proton =
      new TH1F("hPtInPrim_Proton", ";#it{p}_{T} (GeV/#it{c}); Entries",
               pt_Nbins, pt_bins);

  hPtInPrim_sigmap =
      new TH1F("hPtInPrim_Sigmap", ";#it{p}_{T} (GeV/#it{c}); Entries",
               pt_Nbins, pt_bins);

  hPtInPrim_sigmam =
      new TH1F("hPtInPrim_Sigmam", ";#it{p}_{T} (GeV/#it{c}); Entries",
               pt_Nbins, pt_bins);

  hPtInPrim_lambda =
      new TH1F("hPtInPrim_Lambda", ";#it{p}_{T} (GeV/#it{c}); Entries",
               pt_Nbins, pt_bins);

  hPtInPrim_rest = new TH1F(
      "hPtInPrim_Rest", ";#it{p}_{T} (GeV/#it{c}); Entries", pt_Nbins, pt_bins);

  hTruePtvsTrueNchEtaGap = new TH2D(
      "hTruePtvsTrueNchEtaGap",
      "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); #it{p}_{T} (GeV/#it{c})",
      tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);

  for (int i = 0; i < v0m_Nbins; ++i) {
    hDCAxyPri[i] = new TH2F(Form("hDCAxyPri_%s", uc_v0m_bins_name[i]),
                            ";DCA_{xy} (cm);#it{p}_{T} GeV/#it{c}", dcaxy_Nbins,
                            dcaxy_bins, pt_Nbins, pt_bins);
    hDCAxyWeDe[i] = new TH2F(Form("hDCAxyWeDe_%s", uc_v0m_bins_name[i]),
                             ";DCA_{xy} (cm);#it{p}_{T} GeV/#it{c}",
                             dcaxy_Nbins, dcaxy_bins, pt_Nbins, pt_bins);
    hDCAxyMaIn[i] = new TH2F(Form("hDCAxyMaIn_%s", uc_v0m_bins_name[i]),
                             ";DCA_{xy} (cm);#it{p}_{T} GeV/#it{c}",
                             dcaxy_Nbins, dcaxy_bins, pt_Nbins, pt_bins);
  }

  fOutputList->Add(hTrueVtxZ);
  fOutputList->Add(hTrueNch);
  fOutputList->Add(hTrueV0MAmp);
  fOutputList->Add(hNchResponse);
  fOutputList->Add(hPtInPrim_ch);
  fOutputList->Add(hPtInPrim_pion);
  fOutputList->Add(hPtInPrim_kaon);
  fOutputList->Add(hPtInPrim_proton);
  fOutputList->Add(hPtInPrim_sigmap);
  fOutputList->Add(hPtInPrim_sigmam);
  fOutputList->Add(hPtInPrim_lambda);
  fOutputList->Add(hPtInPrim_rest);
  fOutputList->Add(hPtOutAll_ch);
  fOutputList->Add(hPtOutPrim_ch);
  fOutputList->Add(hPtOutPrim_pion);
  fOutputList->Add(hPtOutPrim_kaon);
  fOutputList->Add(hPtOutPrim_proton);
  fOutputList->Add(hPtOutPrim_sigmap);
  fOutputList->Add(hPtOutPrim_sigmam);
  fOutputList->Add(hPtOutPrim_rest);
  fOutputList->Add(hTruePtvsTrueNch);
  fOutputList->Add(hTrueNchEtaPosvsTrueNchEtaNeg);
  fOutputList->Add(hTruePtEtaNegvsTrueNchEtaPos);
  fOutputList->Add(hTruePtEtaPosvsTrueNchEtaNeg);
  fOutputList->Add(hTruePtvsTrueNchEtaGap);
  fOutputList->Add(hTrueNchEtaGap);

  for (int i = 0; i < v0m_Nbins; ++i) {
    fOutputList->Add(hDCAxyPri[i]);
    fOutputList->Add(hDCAxyWeDe[i]);
    fOutputList->Add(hDCAxyMaIn[i]);
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

  GetCalibratedV0Amplitude();
  GetSPDMultiplicity();

  if (isGoodVtxPosMC) {
    if (!fill_corrections) {
      MultiplicityDistributions();
    } else {
      DetectorResponse();
      TrackingEfficiency();
      DCAxyDistributions();
    }
  }

  PostData(1, fOutputList);
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::Terminate(Option_t*) {}

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

    if (TMath::Abs(eta) < 1.0) {
      fTracklets10++;
    }
    if (TMath::Abs(eta) < 1.4) {
      hPhiEtaSPD->Fill(phi, eta);
      fTracklets14++;
    }
    if (0.7 <= TMath::Abs(eta) && TMath::Abs(eta) <= 1.4) {
      fTrackletsEtaGap++;
    }
  }

  hVtxZvsTracklets->Fill(spdVtxZ, fTracklets14);
}

//______________________________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::MultiplicityDistributions() {
  int rec_nch{0};
  int rec_nch_neg_eta{0};
  int rec_nch_pos_eta{0};
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
    rec_nch++;
  }

  hV0Mmult->Fill(fv0mpercentile);
  hV0MAmplitude->Fill(fv0mamplitude);
  hNchvsV0M->Fill(rec_nch, fv0mpercentile);
  hNchvsV0MAmp->Fill(rec_nch, fv0mamplitude);
  hV0MvsV0MAmp->Fill(fv0mamplitude, fv0mpercentile);
  hNchEtaPosvsNchEtaNeg->Fill(rec_nch_pos_eta, rec_nch_neg_eta);
  if (fTrackletsEtaGap > 0) {
    hTrackletsEtaGapvsV0MAmp->Fill(fTrackletsEtaGap, fv0mamplitude);
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
    if (fTrackletsEtaGap > 0) {
      if (TMath::Abs(track->Eta()) <= 0.4) {
        hPtvsTrackletsEtaGap->Fill(fTrackletsEtaGap, pt);
      }
    }
  }
}

//____________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::DCAxyDistributions() {
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

    int label = -1;
    label = TMath::Abs(track->GetLabel());
    if (fMC->IsPhysicalPrimary(label)) {
      hDCAxyPri[index]->Fill(dcaxy, track->Pt());
    } else if (fMC->IsSecondaryFromWeakDecay(label)) {
      hDCAxyWeDe[index]->Fill(dcaxy, track->Pt());
    } else if (fMC->IsSecondaryFromMaterial(label)) {
      hDCAxyMaIn[index]->Fill(dcaxy, track->Pt());
    } else {
      continue;
    }
  }
}

//____________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::TrackingEfficiency() {
  if (fv0mpercentile < 0.0 || fv0mpercentile > 5.0) {
    return;
  }

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
    hPtOutAll_ch->Fill(track->Pt());
    if (!fMC->IsPhysicalPrimary(label)) {
      continue;
    }

    hPtOutPrim_ch->Fill(track->Pt());
    int partPDG = TMath::Abs(particle->GetPdgCode());
    if (partPDG == 211) {
      hPtOutPrim_pion->Fill(track->Pt());  // pions
    } else if (partPDG == 321) {
      hPtOutPrim_kaon->Fill(track->Pt());  // kaons
    } else if (partPDG == 2212) {
      hPtOutPrim_proton->Fill(track->Pt());  // protons
    } else if (partPDG == 3222) {
      hPtOutPrim_sigmap->Fill(track->Pt());  // sigma plus
    } else if (partPDG == 3112) {
      hPtOutPrim_sigmam->Fill(track->Pt());  // sigma minus
    } else {
      hPtOutPrim_rest->Fill(track->Pt());  // rest of the charged particles
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

    int partPDG = TMath::Abs(particle->PdgCode());

    if (partPDG == 3122) {
      hPtInPrim_lambda->Fill(particle->Pt());  // Lambda
    }
    if (particle->Charge() == 0) {
      continue;
    }

    hPtInPrim_ch->Fill(particle->Pt());
    if (partPDG == 211) {
      hPtInPrim_pion->Fill(particle->Pt());  // pions
    } else if (partPDG == 321) {
      hPtInPrim_kaon->Fill(particle->Pt());  // kaons
    } else if (partPDG == 2212) {
      hPtInPrim_proton->Fill(particle->Pt());  // protons
    } else if (partPDG == 3222) {
      hPtInPrim_sigmap->Fill(particle->Pt());  // sigma plus
    } else if (partPDG == 3112) {
      hPtInPrim_sigmam->Fill(particle->Pt());  // sigma minus
    } else {
      hPtInPrim_rest->Fill(particle->Pt());  // rest of the charged particles
    }
  }
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
    if ((2.8 < particle->Eta() && particle->Eta() < 5.1) ||
        (-3.7 < particle->Eta() && particle->Eta() < -1.7)) {
      if (fMC->IsPhysicalPrimary(i)) {
        fTrueV0++;
      }
    }
    if (TMath::Abs(particle->Eta()) < 1.4) {
      if (fMC->IsPhysicalPrimary(i)) {
        fTrueNch14++;
      }
    }
    if (TMath::Abs(particle->Eta()) < 1.0) {
      if (fMC->IsPhysicalPrimary(i)) {
        fTrueNch10++;
      }
    }
    if (0.7 <= TMath::Abs(particle->Eta()) &&
        TMath::Abs(particle->Eta()) <= 1.4) {
      if (fMC->IsPhysicalPrimary(i)) {
        fTrueNchEtaGap++;
      }
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

    fTrueNch++;
    if (fEtaMin <= particle->Eta() && particle->Eta() < 0.0) {
      fTrueNchEtaNeg++;
    }
    if (particle->Eta() >= 0.0 && particle->Eta() <= fEtaMax) {
      fTrueNchEtaPos++;
    }
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
    if (TMath::Abs(particle->Eta()) <= 0.4) {
      hTruePtvsTrueNchEtaGap->Fill(fTrueNchEtaGap, particle->Pt());
    }
    if (fEtaMin <= particle->Eta() && particle->Eta() < 0.0) {
      hTruePtEtaNegvsTrueNchEtaPos->Fill(fTrueNchEtaPos, particle->Pt());
    }
    if (particle->Eta() >= 0.0 && particle->Eta() <= fEtaMax) {
      hTruePtEtaPosvsTrueNchEtaNeg->Fill(fTrueNchEtaNeg, particle->Pt());
    }
  }

  hTrueNch->Fill(fTrueNch);
  hTrueV0MAmp->Fill(fTrueV0);
  hTrueNchEtaGap->Fill(fTrueNchEtaGap);
  hTrueNchEtaPosvsTrueNchEtaNeg->Fill(fTrueNchEtaPos, fTrueNchEtaNeg);
}

//____________________________________________________________

void AliAnalysisTaskDataSpeedOfSoundSim::DetectorResponse() {
  fRecNch = 0;
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
    fRecNch++;
  }
  hNchResponse->Fill(fRecNch, fTrueNch);
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
