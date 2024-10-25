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
      fTowerEnergy(true),
      fUseZEMcut(true),
      fIsZNhit(false),
      fIsZEMhit(false),
      fZEMcut(800.0),
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
      fZNC(0.0),
      fZNA(0.0),
      fZPC(0.0),
      fZPA(0.0),
      fZDC(0.0),
      fZEM(0.0),
      fMultSelection(0x0),
      hNch(0),
      pNchvsV0MAmp(0),
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
      hEtFull(0),
      hEtEtaGap(0),
      hPtvsSPDFull(0),
      hPtvsSPDEtaAdj(0),
      hPtvsSPDEtaGapW(0),
      hPtvsEtFull(0),
      hPtvsEtEtaGap(0),
      hZNC(0),
      hZNA(0),
      hZPC(0),
      hZPA(0),
      hZDC(0),
      hZNCvspT(0),
      hZNAvspT(0),
      hZPCvspT(0),
      hZPAvspT(0),
      hZDCvspT(0),
      hZDCvsZEM(0),
      hV0MvsZDC(0),
      hV0MvsZNC(0),
      hV0MvsZNA(0),
      hV0MvsZPC(0),
      hV0MvsZPA(0),
      hZNhits(0),
      hZEMhits(0),
      hZNZEMhits(0),
      pZDCvspT(0),
      pZDCvsZEM(0),
      pZDCvsNch(0),
      pZVtxvsSPDClus(0),
      pSPDClusvsEta(0),
      pZNCvsV0Amp(0),
      pZNAvsV0Amp(0),
      pZPCvsV0Amp(0),
      pZPAvsV0Amp(0),
      pZNCvsNch(0),
      pZNAvsNch(0),
      pZPCvsNch(0),
      pZPAvsNch(0),
      pZNCvsEt(0),
      pZNAvsEt(0),
      pZPCvsEt(0),
      pZPAvsEt(0),
      pZNCvsSPD(0),
      pZNAvsSPD(0),
      pZPCvsSPD(0),
      pZPAvsSPD(0),
      fSPDVtxCut(3.0),
      fSPDFull(0),
      fSPDEtaAdj(0),
      fSPDEtaGapW(0),
      fSPDEtaGapWW(0) {
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
      fTowerEnergy(true),
      fUseZEMcut(true),
      fIsZNhit(false),
      fIsZEMhit(false),
      fZEMcut(800.0),
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
      fZNC(0.0),
      fZNA(0.0),
      fZPC(0.0),
      fZPA(0.0),
      fZDC(0.0),
      fZEM(0.0),
      fMultSelection(0x0),
      hNch(0),
      pNchvsV0MAmp(0),
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
      hEtFull(0),
      hEtEtaGap(0),
      hPtvsSPDFull(0),
      hPtvsSPDEtaAdj(0),
      hPtvsSPDEtaGapW(0),
      hPtvsEtFull(0),
      hPtvsEtEtaGap(0),
      hZNC(0),
      hZNA(0),
      hZPC(0),
      hZPA(0),
      hZDC(0),
      hZNCvspT(0),
      hZNAvspT(0),
      hZPCvspT(0),
      hZPAvspT(0),
      hZDCvspT(0),
      hZDCvsZEM(0),
      hV0MvsZDC(0),
      hV0MvsZNC(0),
      hV0MvsZNA(0),
      hV0MvsZPC(0),
      hV0MvsZPA(0),
      hZNhits(0),
      hZEMhits(0),
      hZNZEMhits(0),
      pZDCvspT(0),
      pZDCvsZEM(0),
      pZDCvsNch(0),
      pZVtxvsSPDClus(0),
      pSPDClusvsEta(0),
      pZNCvsV0Amp(0),
      pZNAvsV0Amp(0),
      pZPCvsV0Amp(0),
      pZPAvsV0Amp(0),
      pZNCvsNch(0),
      pZNAvsNch(0),
      pZPCvsNch(0),
      pZPAvsNch(0),
      pZNCvsEt(0),
      pZNAvsEt(0),
      pZPCvsEt(0),
      pZPAvsEt(0),
      pZNCvsSPD(0),
      pZNAvsSPD(0),
      pZPCvsSPD(0),
      pZPAvsSPD(0),
      fSPDVtxCut(3.0),
      fSPDFull(0),
      fSPDEtaAdj(0),
      fSPDEtaGapW(0),
      fSPDEtaGapWW(0) {
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

  /*const int nPtbins = 11;*/
  /*double Ptbins[nPtbins + 1] = {0.15, 0.2, 0.4, 0.6, 0.8,  1.0,*/
  /*                              1.5,  2.0, 3.0, 5.0, 10.0, 50.0};*/

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

  // Nch (0.5<eta<0.8)
  constexpr int NchOneSide_Nbins{250};
  double NchOneSide_bins[NchOneSide_Nbins + 1] = {0};
  for (int i = 0; i <= NchOneSide_Nbins; ++i) {
    NchOneSide_bins[i] = 0.5 + (2.0 * i);
  }

  // SPD Tracklets (|eta|<0.8)
  constexpr int SPD0p8_Nbins{1650};
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

  constexpr int EtEtaGap_Nbins{400};
  double EtEtaGap_bins[EtEtaGap_Nbins + 1] = {0};
  for (int i = 0; i <= EtEtaGap_Nbins; ++i) {
    EtEtaGap_bins[i] = 0.0 + (2.0 * i);
  }

  // Et (0.5<eta<0.8)
  constexpr int EtOneSide_Nbins{200};
  double EtOneSide_bins[EtOneSide_Nbins + 1] = {0};
  for (int i = 0; i <= EtOneSide_Nbins; ++i) {
    EtOneSide_bins[i] = 0.0 + (2.0 * i);
  }

  constexpr double v0mAmp_width{25.0};
  constexpr int v0mAmp_Nbins{1720};
  double v0mAmp_bins[v0mAmp_Nbins + 1] = {0};
  for (int i = 0; i <= v0mAmp_Nbins; ++i) {
    v0mAmp_bins[i] = 0.0 + i * v0mAmp_width;
  }

  constexpr int ZN_Nbins{1500};
  double ZN_bins[ZN_Nbins + 1] = {0};
  for (int i = 0; i <= ZN_Nbins; ++i) {
    ZN_bins[i] = 0.0 + (0.1 * i);
  }

  constexpr int ZP_Nbins{510};
  double ZP_bins[ZP_Nbins + 1] = {0};
  for (int i = 0; i <= ZP_Nbins; ++i) {
    ZP_bins[i] = 0.0 + (0.1 * i);
  }

  constexpr int ZDC_Nbins{2550};
  double ZDC_bins[ZDC_Nbins + 1] = {0};
  for (int i = 0; i <= ZDC_Nbins; ++i) {
    ZDC_bins[i] = 0.0 + (0.1 * i);
  }

  constexpr int dcaxy_Nbins{100};
  double dcaxy_bins[dcaxy_Nbins + 1] = {0};
  for (int i = 0; i <= dcaxy_Nbins; ++i) {
    dcaxy_bins[i] = -3.0 + (0.06 * i);
  }

  constexpr int v0m_Nbins0100{100};
  double v0m_bins0100[v0m_Nbins0100 + 1] = {0.};
  for (int i = 0; i <= v0m_Nbins0100; ++i) {
    v0m_bins0100[i] = (double)i;
  }
  std::string title{"Common PMT:TowerEnergy[0]"};
  if (!fTowerEnergy) title = "All five PMTs: fZDCEnergy";

  hV0Percentile = new TH1F("hV0Percentile", ";V0M (%);Entries", v0m_Nbins0100,
                           v0m_bins0100);

  hNch = new TH1F("hTPCFull", ";#it{N}_{ch} (|#eta|#leq0.8); Counts;",
                  nch_Nbins, nch_bins);

  pNchvsV0MAmp = new TProfile(
      "pNchvsV0MAmp", ";V0M amplitude; #LT#it{N}_{ch}#GT (|#eta|#leq0.8);",
      v0mAmp_Nbins, v0mAmp_bins);

  pV0MAmpChannel =
      new TProfile("pV0MAmpChannel", ";V0 Channel; Amplitude;", 64, -0.5, 63.5);

  hV0MAmplitude =
      new TH1F("hV0MAmp", ";V0M amplitude;Counts", v0mAmp_Nbins, v0mAmp_bins);

  hPtvsNch = new TH2F("hPtvsTPCFull",
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

  hPtvsV0MAmp = new TH2F(
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
      new TH2F("hPtvsSPDEtaGap",
               "; #it{N}_{tracklet} (0.5#leq|#eta|#leq0.8); #it{p}_{T} "
               "(|#eta|#leq0.4, GeV/#it{c})",
               SPDEtaGap_Nbins, SPDEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsTracksEtaGapTPC = new TH2F(
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
      new TH1F("hSPDEtaAdj", "; #it{N}_{ch} (0.3<|#eta|#leq0.6); Entries",
               SPD0p4_Nbins, SPD0p4_bins);

  hSPDEtaGapW =
      new TH1F("hSPDEtaGapW", "; #it{N}_{ch} (0.7#leq|#eta|#leq1); Entries",
               SPDEtaGap_Nbins, SPDEtaGap_bins);

  hEtFull = new TH1F("hEtFull", ";#it{E}_{T} (|#eta|#leq0.8); Entries",
                     Et_Nbins, Et_bins);

  hEtEtaGap =
      new TH1F("hEtEtaGap", ";#it{E}_{T} (0.5#leq|#eta|#leq0.8);Entries",
               EtEtaGap_Nbins, EtEtaGap_bins);

  hPtvsSPDFull = new TH2F("hPtvsSPDFull",
                          "; #it{N}_{ch} (|#eta|#leq0.8); #it{p}_{T} "
                          "(|#eta|#leq0.8, GeV/#it{c})",
                          SPD0p8_Nbins, SPD0p8_bins, pt_Nbins, pt_bins);

  hPtvsSPDEtaAdj = new TH2F("hPtvsSPDEtaAdj",
                            "; #it{N}_{ch} (0.3<|#eta|#leq0.6); #it{p}_{T} "
                            "(|#eta|#leq0.3, GeV/#it{c})",
                            SPD0p4_Nbins, SPD0p4_bins, pt_Nbins, pt_bins);

  hPtvsSPDEtaGapW =
      new TH2F("hPtvsSPDEtaGapW",
               "; #it{N}_{ch} (0.7#leq|#eta|#leq1); #it{p}_{T} "
               "(|#eta|#leq0.3, GeV/#it{c})",
               SPDEtaGap_Nbins, SPDEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsEtFull = new TH2F(
      "hPtvsEtFull",
      ";#it{E}_{T} (|#eta|#leq0.8); #it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      Et_Nbins, Et_bins, pt_Nbins, pt_bins);

  hZNC = new TH1F("hZNC", ";ZNC/2.511 [TeV]; Entries;", ZN_Nbins, ZN_bins);
  hZNA = new TH1F("hZNA", ";ZNA/2.511 [TeV]; Entries;", ZN_Nbins, ZN_bins);
  hZPC = new TH1F("hZPC", ";ZPC/2.511 [TeV]; Entries;", ZP_Nbins, ZP_bins);
  hZPA = new TH1F("hZPA", ";ZPA/2.511 [TeV]; Entries;", ZP_Nbins, ZP_bins);
  hZDC = new TH1F("hZDC", ";ZDC/2.511 [TeV]; Entries;", ZDC_Nbins, ZDC_bins);

  hZNZEMhits =
      new TH1F("hZNZEMhits", "Events with ZEMs and ZNs hits;V0M percentile; ",
               v0m_Nbins0100, v0m_bins0100);
  hZEMhits =
      new TH2F("hZEMhits",
               "Bin 1 = YES, Bin 2 = NO;V0M percentile; Hits in ZEM1 OR ZEM2?",
               v0m_Nbins0100, v0m_bins0100, 2, 0., 2.);
  hZNhits =
      new TH2F("hZNhits",
               "Bin 1 = YES, Bin 2 = NO;V0M percentile; Hits in ZNC AND ZNA?;",
               v0m_Nbins0100, v0m_bins0100, 2, 0., 2.);

  hZNCvspT = new TH2F("hZNCvspT", ";ZNC/2.511 [TeV]; p_{T} [GeV/c];", ZN_Nbins,
                      ZN_bins, pt_Nbins, pt_bins);
  hZNAvspT = new TH2F("hZNAvspT", ";ZNA/2.511 [TeV]; p_{T} [GeV/c];", ZN_Nbins,
                      ZN_bins, pt_Nbins, pt_bins);
  hZPCvspT = new TH2F("hZPCvspT", ";ZPC/2.511 [TeV]; p_{T} [GeV/c];", ZP_Nbins,
                      ZP_bins, pt_Nbins, pt_bins);
  hZPAvspT = new TH2F("hZPAvspT", ";ZPA/2.511 [TeV]; p_{T} [GeV/c];", ZP_Nbins,
                      ZP_bins, pt_Nbins, pt_bins);
  hZDCvspT = new TH2F("hZDCvspT", ";ZDC/2.511 [TeV]; p_{T} [GeV/c];", ZDC_Nbins,
                      ZDC_bins, pt_Nbins, pt_bins);
  hZDCvsZEM = new TH2F("hZDCvsZEM", ";ZEM (a.u.); ZDC/2.511 [TeV];", 162, 0.,
                       1620., 255, 0., 255.);
  hV0MvsZDC =
      new TH2F("hV0MvsZDC", "With ZEM cut;V0M centrality;ZDC/2.511 [TeV];",
               v0m_Nbins0100, v0m_bins0100, 255, 0., 255.);

  hV0MvsZNC =
      new TH2F("hV0MvsZNC", "With ZEM cut;V0M centrality;ZNC/2.511 [TeV];",
               v0m_Nbins0100, v0m_bins0100, 150, 0., 150.);
  hV0MvsZNA =
      new TH2F("hV0MvsZNA", "With ZEM cut;V0M centrality;ZNA/2.511 [TeV];",
               v0m_Nbins0100, v0m_bins0100, 150, 0., 150.);
  hV0MvsZPC =
      new TH2F("hV0MvsZPC", "With ZEM cut;V0M centrality;ZPC/2.511 [TeV];",
               v0m_Nbins0100, v0m_bins0100, 150, 0., 51.);
  hV0MvsZPA =
      new TH2F("hV0MvsZPA", "With ZEM cut;V0M centrality;ZPA/2.511 [TeV];",
               v0m_Nbins0100, v0m_bins0100, 150, 0., 51.);

  pZDCvspT = new TProfile("pZDCvspT", ";ZDC/2.511 [TeV]; <pT> [GeV/c];", 255,
                          0., 255.);
  pZDCvsNch = new TProfile("pZDCvsNch",
                           "With ZEM cut;Nch (|#eta|<0.8);<ZDC>/2.511 [TeV];",
                           nch_Nbins, nch_bins);
  pZDCvsZEM =
      new TProfile("pZDCvsZEM", "Without ZEM cut;ZEM (a.u.);<ZDC>/2.511 (TeV);",
                   162, 0., 1620.);

  hPtvsEtEtaGap = new TH2F("hPtvsEtEtaGap",
                           ";#it{E}_{T} (0.5#leq|#eta|#leq0.8); #it{p}_{T} "
                           "(|#eta|#leq0.3, GeV/#it{c})",
                           EtEtaGap_Nbins, EtEtaGap_bins, pt_Nbins, pt_bins);
  pZVtxvsSPDClus =
      new TProfile("pZVtxvsSPDClus", ";Z_{SPD} Vertex (cm); <SPD clusters>;",
                   40, -10.0, 10.0);
  pSPDClusvsEta =
      new TProfile("pSPDClusvsEta", ";#eta; <SPD clusters>;", 30, -1.5, 1.5);

  pZNCvsV0Amp = new TProfile(
      "pZNCvsV0Amp", Form("%s;V0 Amp;#it{E}_{ZNC} [TeV];", title.c_str()),
      v0mAmp_Nbins, v0mAmp_bins, 0., 251.);
  pZNAvsV0Amp = new TProfile(
      "pZNAvsV0Amp", Form("%s;V0 Amp;#it{E}_{ZNA} [TeV];", title.c_str()),
      v0mAmp_Nbins, v0mAmp_bins, 0., 252.);
  pZPCvsV0Amp = new TProfile(
      "pZPCvsV0Amp", Form("%s;V0M Per;#it{E}_{ZPC} [TeV];", title.c_str()),
      v0mAmp_Nbins, v0mAmp_bins, 0., 76.);
  pZPAvsV0Amp = new TProfile(
      "pZPAvsV0Amp", Form("%s;V0M Per;#it{E}_{ZPA} [TeV];", title.c_str()),
      v0mAmp_Nbins, v0mAmp_bins, 0., 76.);
  pZNCvsNch = new TProfile("pZNCvsNch",
                           Form("%s;Nch;#it{E}_{ZNC} [TeV];", title.c_str()),
                           nch_Nbins, nch_bins, 0., 251.);
  pZNAvsNch = new TProfile("pZNAvsNch",
                           Form("%s;Nch;#it{E}_{ZNA} [TeV];", title.c_str()),
                           nch_Nbins, nch_bins, 0., 252.);
  pZPCvsNch = new TProfile("pZPCvsNch",
                           Form("%s;Nch;#it{E}_{ZPC} [TeV];", title.c_str()),
                           nch_Nbins, nch_bins, 0., 76.);
  pZPAvsNch = new TProfile("pZPAvsNch",
                           Form("%s;Nch;#it{E}_{ZPA} [TeV];", title.c_str()),
                           nch_Nbins, nch_bins, 0., 76.);
  pZNCvsEt =
      new TProfile("pZNCvsEt", Form("%s;Et;#it{E}_{ZNC} [TeV];", title.c_str()),
                   Et_Nbins, Et_bins, 0., 251.);
  pZNAvsEt =
      new TProfile("pZNAvsEt", Form("%s;Et;#it{E}_{ZNA} [TeV];", title.c_str()),
                   Et_Nbins, Et_bins, 0., 252.);
  pZPCvsEt =
      new TProfile("pZPCvsEt", Form("%s;Et;#it{E}_{ZPC} [TeV];", title.c_str()),
                   Et_Nbins, Et_bins, 0., 76.);
  pZPAvsEt =
      new TProfile("pZPAvsEt", Form("%s;Et;#it{E}_{ZPA} [TeV];", title.c_str()),
                   Et_Nbins, Et_bins, 0., 76.);

  pZNCvsSPD = new TProfile("pZNCvsSPD",
                           Form("%s;Nch;#it{E}_{ZNC} [TeV];", title.c_str()),
                           SPD0p8_Nbins, SPD0p8_bins, 0., 251.);
  pZNAvsSPD = new TProfile("pZNAvsSPD",
                           Form("%s;Nch;#it{E}_{ZNA} [TeV];", title.c_str()),
                           SPD0p8_Nbins, SPD0p8_bins, 0., 252.);
  pZPCvsSPD = new TProfile("pZPCvsSPD",
                           Form("%s;Nch;#it{E}_{ZPC} [TeV];", title.c_str()),
                           SPD0p8_Nbins, SPD0p8_bins, 0., 76.);
  pZPAvsSPD = new TProfile("pZPAvsSPD",
                           Form("%s;Nch;#it{E}_{ZPA} [TeV];", title.c_str()),
                           SPD0p8_Nbins, SPD0p8_bins, 0., 76.);

  fOutputList->Add(hBestVtxZ);
  fOutputList->Add(pZVtxvsSPDClus);
  fOutputList->Add(pSPDClusvsEta);
  fOutputList->Add(pNchvsV0MAmp);
  fOutputList->Add(hV0Percentile);

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

  fOutputList->Add(hEtFull);
  fOutputList->Add(hPtvsEtFull);

  fOutputList->Add(hEtEtaGap);
  fOutputList->Add(hPtvsEtEtaGap);

  fOutputList->Add(hZNC);
  fOutputList->Add(hZNCvspT);
  fOutputList->Add(hZNA);
  fOutputList->Add(hZNAvspT);
  fOutputList->Add(hZPC);
  fOutputList->Add(hZPCvspT);
  fOutputList->Add(hZPA);
  fOutputList->Add(hZPAvspT);
  fOutputList->Add(hZDC);
  fOutputList->Add(hZDCvspT);
  fOutputList->Add(pZDCvspT);

  fOutputList->Add(hV0MvsZNC);
  fOutputList->Add(hV0MvsZNA);
  fOutputList->Add(hV0MvsZPC);
  fOutputList->Add(hV0MvsZPA);
  fOutputList->Add(hV0MvsZDC);
  fOutputList->Add(hZDCvsZEM);
  fOutputList->Add(pZDCvsZEM);
  fOutputList->Add(pZDCvsNch);
  fOutputList->Add(hZNZEMhits);
  fOutputList->Add(hZNhits);
  fOutputList->Add(hZEMhits);

  fOutputList->Add(pZNCvsV0Amp);
  fOutputList->Add(pZNAvsV0Amp);
  fOutputList->Add(pZNCvsNch);
  fOutputList->Add(pZNAvsNch);
  fOutputList->Add(pZNCvsEt);
  fOutputList->Add(pZNAvsEt);
  fOutputList->Add(pZNCvsSPD);
  fOutputList->Add(pZNAvsSPD);
  fOutputList->Add(pZPCvsV0Amp);
  fOutputList->Add(pZPAvsV0Amp);
  fOutputList->Add(pZPCvsNch);
  fOutputList->Add(pZPAvsNch);
  fOutputList->Add(pZPCvsEt);
  fOutputList->Add(pZPAvsEt);
  fOutputList->Add(pZPCvsSPD);
  fOutputList->Add(pZPAvsSPD);

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

  //! DCAxy templates MC and Data
  // DCAxyDistributions();

  //! ZDC
  GetZDC();

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
void AliAnalysisTaskDataSpeedOfSound::GetCalibratedV0Amplitude() {
  float mV0M{0.0};
  for (int i = 0; i < 64; i++) {
    mV0M += fESD->GetVZEROEqMultiplicity(i);
    pV0MAmpChannel->Fill(i, fESD->GetVZEROEqMultiplicity(i));
  }
  fv0mamplitude = mV0M;
}
//______________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound::GetZDC() {
  AliESDZDC* esdZDC = fESD->GetESDZDC();
  if (!esdZDC) {
    return;
  }

  fZNC = -999.;
  fZNA = -999.;
  fZPA = -999.;
  fZPC = -999.;
  fZDC = -999.;
  fZEM = -999.;
  fIsZNhit = false;
  fIsZEMhit = false;

  const double eA{2.511};
  const double gev2tev{1. / 1000.};
  const double zNC{esdZDC->GetZDCN1Energy()};
  const double zNA{esdZDC->GetZDCN2Energy()};
  const double zPC{esdZDC->GetZDCP1Energy()};
  const double zPA{esdZDC->GetZDCP2Energy()};

  const double* towZNC{esdZDC->GetZN1TowerEnergy()};
  const double* towZPC{esdZDC->GetZP1TowerEnergy()};
  const double* towZNA{esdZDC->GetZN2TowerEnergy()};
  const double* towZPA{esdZDC->GetZP2TowerEnergy()};
  double zem{esdZDC->GetZEM1Energy() + esdZDC->GetZEM2Energy()};

  if (esdZDC->IsZNChit() && esdZDC->IsZNAhit()) {
    hZNhits->Fill(fv0mpercentile, 0.5);
    fIsZNhit = true;
  } else {
    hZNhits->Fill(fv0mpercentile, 1.5);
    fIsZNhit = false;
  }

  if (esdZDC->IsZEM1hit() || esdZDC->IsZEM2hit()) {
    hZEMhits->Fill(fv0mpercentile, 0.5);
    fIsZEMhit = true;
  } else {
    hZEMhits->Fill(fv0mpercentile, 1.5);
    fIsZEMhit = false;
  }

  double znc{999.};
  double zna{999.};
  double zpc{999.};
  double zpa{999.};
  if (fTowerEnergy) {
    znc = towZNC[0] * gev2tev;
    zna = towZNA[0] * gev2tev;
    zpc = towZPC[0] * gev2tev;
    zpa = towZPA[0] * gev2tev;
  } else {
    znc = zNC * gev2tev;
    zna = zNA * gev2tev;
    zpc = zPC * gev2tev;
    zpa = zPA * gev2tev;
  }

  znc /= eA;
  zna /= eA;
  zpc /= eA;
  zpa /= eA;

  //! Efficiency corrected
  //! https://alice-notes.web.cern.ch/system/files/notes/analysis/1062/2021-06-03-ALICE_AN_ZDC_efficiency.pdf
  fZNA = zna;
  fZPA = zpa;
  fZNC = znc;
  fZPC = zpc;
  fZDC = (fZNC + fZNA + fZPC + fZPA);
  fZEM = zem;

  if (fUseZEMcut && fIsZEMhit && fIsZNhit) {
    //! Use ZEM cut to select low ZDC events
    if (fZEM > fZEMcut) {
      hZNC->Fill(fZNC);
      hZNA->Fill(fZNA);
      hZPC->Fill(fZPC);
      hZPA->Fill(fZPA);
      hZDC->Fill(fZDC);
      hZDCvsZEM->Fill(fZEM, fZDC);
    } else {
      hZNC->Fill(149.);
      hZNA->Fill(149.);
      hZPC->Fill(50.);
      hZPA->Fill(50.);
      hZDC->Fill(254.);
      hZDCvsZEM->Fill(fZEM, 254.);
    }
  }

  if (fIsZEMhit && fIsZNhit) {
    hZNZEMhits->Fill(fv0mpercentile);
    pZDCvsZEM->Fill(fZEM, fZDC);
    pZNCvsV0Amp->Fill(fv0mamplitude, fZNC * eA);
    pZNAvsV0Amp->Fill(fv0mamplitude, fZNA * eA);
    pZPCvsV0Amp->Fill(fv0mamplitude, fZPC * eA);
    pZPAvsV0Amp->Fill(fv0mamplitude, fZPA * eA);
    pZNCvsSPD->Fill(fSPDFull, fZNC * eA);
    pZNAvsSPD->Fill(fSPDFull, fZNA * eA);
    pZPCvsSPD->Fill(fSPDFull, fZPC * eA);
    pZPAvsSPD->Fill(fSPDFull, fZPA * eA);
  }
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

    if (TMath::Abs(eta) > 0.3 && TMath::Abs(eta) <= 0.6) {
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
  const double eA{2.511};
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

    if (fUseZEMcut && fIsZEMhit && fIsZNhit) {
      //! Use ZEM cut to select low ZDC events
      if (fZEM > fZEMcut) {
        hZNCvspT->Fill(fZNC, pt);
        hZNAvspT->Fill(fZNA, pt);
        hZPCvspT->Fill(fZPC, pt);
        hZPAvspT->Fill(fZPA, pt);
        hZDCvspT->Fill(fZDC, pt);
        pZDCvspT->Fill(fZDC, pt);
      } else {
        hZNCvspT->Fill(149., pt);
        hZNAvspT->Fill(149., pt);
        hZPCvspT->Fill(50., pt);
        hZPAvspT->Fill(50., pt);
        hZDCvspT->Fill(254., pt);
      }
    }

    //! SPD Eta Adjacent
    if (TMath::Abs(track->Eta()) <= fEtaCutForpTwTPCGap) {
      hPtvsSPDEtaAdj->Fill(fSPDEtaAdj, pt);
    }

    //! SPD Eta Gap
    if (TMath::Abs(track->Eta()) <= fEtaCutForpTwSPDGap) {
      hPtvsTrackletsEtaGap->Fill(fTrackletsEtaGap, pt);
      pPtvsTrackletsEtaGap->Fill(fTrackletsEtaGap, pt);
      hPtvsSPDEtaGapW->Fill(fSPDEtaGapW, pt);
    }

    //! Nch 0.5<=|eta|<=0.8 and Spectra |eta|<=0.3
    if (TMath::Abs(track->Eta()) <= fEtaCutForpTwTPCGap) {
      hPtvsTracksEtaGapTPC->Fill(fTracksEtaGapTPC, pt);
      pPtvsTracksEtaGapTPC->Fill(fTracksEtaGapTPC, pt);
      hPtvsEtEtaGap->Fill(etetagap, pt);
    }
  }

  hV0Percentile->Fill(fv0mpercentile);
  hV0MAmplitude->Fill(fv0mamplitude);
  pNchvsV0MAmp->Fill(fv0mamplitude, rec_nch);
  hNch->Fill(rec_nch);

  if (fUseZEMcut && fIsZEMhit && fIsZNhit) {
    //! Use ZEM cut to select low ZDC events
    if (fZEM > fZEMcut) {
      pZDCvsNch->Fill(rec_nch, fZDC);
      hV0MvsZDC->Fill(fv0mpercentile, fZDC);
      hV0MvsZNC->Fill(fv0mpercentile, fZNC);
      hV0MvsZNA->Fill(fv0mpercentile, fZNA);
      hV0MvsZPC->Fill(fv0mpercentile, fZPC);
      hV0MvsZPA->Fill(fv0mpercentile, fZPA);
    }
  }

  hTrackletsEtaGap->Fill(fTrackletsEtaGap);
  hTracksEtaGapTPC->Fill(fTracksEtaGapTPC);

  hSPDFull->Fill(fSPDFull);
  hSPDEtaAdj->Fill(fSPDEtaAdj);
  hSPDEtaGapW->Fill(fSPDEtaGapW);

  hEtFull->Fill(etfull);
  hEtEtaGap->Fill(etetagap);

  if (fIsZEMhit && fIsZNhit) {
    pZNCvsNch->Fill(rec_nch, fZNC * eA);
    pZNAvsNch->Fill(rec_nch, fZNA * eA);
    pZPCvsNch->Fill(rec_nch, fZPC * eA);
    pZPAvsNch->Fill(rec_nch, fZPA * eA);
    pZNCvsEt->Fill(etfull, fZNC * eA);
    pZNAvsEt->Fill(etfull, fZNA * eA);
    pZPCvsEt->Fill(etfull, fZPC * eA);
    pZPAvsEt->Fill(etfull, fZPA * eA);
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
