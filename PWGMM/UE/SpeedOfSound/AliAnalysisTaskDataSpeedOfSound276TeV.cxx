/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
 * provided "as is" without express or implied warranty.                  *
 * **************************************************************************/

#include "AliAnalysisTaskDataSpeedOfSound276TeV.h"

// ROOT includes
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TList.h>
#include <TMath.h>
#include <TParticle.h>
#include <TProfile.h>
#include <TTree.h>

// AliRoot includes
#include <AliAODMCHeader.h>
#include <AliAODPid.h>
#include <AliAODTrack.h>
#include <AliAODVZERO.h>
#include <AliAODVertex.h>
#include <AliAnalysisFilter.h>
#include <AliAnalysisManager.h>
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDVZERO.h>
#include <AliESDVertex.h>
#include <AliESDtrackCuts.h>
#include <AliESDv0.h>
#include <AliExternalTrackParam.h>
#include <AliGenDPMjetEventHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliHeader.h>
#include <AliKFVertex.h>
#include <AliLog.h>
#include <AliMCEvent.h>
#include <AliMCEventHandler.h>
#include <AliStack.h>
#include <TTreeStream.h>

#include "AliCentrality.h"

// STL includes
#include <iostream>
using namespace std;

//
// Responsible:
// Omar Vazquez (U. Houston)

ClassImp(AliAnalysisTaskDataSpeedOfSound276TeV)
    //_____________________________________________________________________________
    AliAnalysisTaskDataSpeedOfSound276TeV::
        AliAnalysisTaskDataSpeedOfSound276TeV()
    : AliAnalysisTaskSE(),
      fESD(0x0),
      fAOD(0x0),
      fMC(0x0),
      fMCStack(0x0),
      fMCArray(0x0),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fCentEst("V0M"),
      fAnalysisType("ESD"),
      fAnalysisMC(kFALSE),
      fAnalysisPbPb(kFALSE),
      ftrigBit(0x0),
      fRandom(0x0),
      fPileUpRej(kFALSE),
      fVtxCut(10.0),
      fEtaCut(0.9),
      fEtaMin(-0.8),
      fEtaMax(0.8),
      fPtMin(0.15),
      fEtaGappT(0.4),
      fMinCent(0.0),
      fMaxCent(100.0),
      fStoreMcIn(kFALSE),
      fEtaGapNchMin(0.7),
      fEtaGapNchMax(1.4),
      fMcProcessType(-999),
      fTriggeredEventMB(-999),
      fVtxStatus(-999),
      fZvtx(-999),
      fZvtxMC(-999),
      fRun(-999),
      fEventId(-999),
      fv0mamplitude(0.0),
      fv0mpercentile(0.0),
      fTracklets14(-999),
      fTracklets10(-999),
      fTrackletsEtaGap(-999),
      fListOfObjects(0),
      fEvents(0x0),
      fVtx(0x0),
      fVtxMC(0x0),
      fVtxBeforeCuts(0x0),
      fVtxAfterCuts(0x0),
      fn1(0x0),
      fcent(0x0),
      hPhi(0x0),
      hPhiEtaSPD(0x0),
      hPhiEtaGapSPD(0x0),
      hVtxZvsTracklets(0x0),
      pV0MAmpChannel(0x0),
      hPtEtaNegvsNchEtaPos(0x0),
      pPtEtaNegvsNchEtaPos(0x0),
      hPtEtaPosvsNchEtaNeg(0x0),
      pPtEtaPosvsNchEtaNeg(0x0),
      pPtvsNch(0x0),
      pPtvsV0MAmp(0x0),
      hPtvsV0MAmp(0x0),
      hPtvsTracklets10(0x0),
      pPtvsTracklets10(0x0),
      hPtvsTracklets14(0x0),
      pPtvsTracklets14(0x0),
      hPtvsTrackletsEtaGap(0x0),
      pPtvsTrackletsEtaGap(0x0),
      hV0Mmult(0x0),
      hV0MAmplitude(0x0),
      hNchvsV0M(0x0),
      hNchvsV0MAmp(0x0),
      hV0MvsV0MAmp(0x0),
      hNchEtaPosvsNchEtaNeg(0x0),
      hTrackletsvsV0MAmp10(0x0),
      hTrackletsvsV0MAmp14(0x0),
      hTrackletsEtaGap(0x0),
      hTrackletsvsV0MAmpEtaGap(0x0),
      hDaDCAxy(0x0)

{
  // default constructor
  for (Int_t pid = 0; pid < 8; ++pid) {
    hMcIn[pid] = 0;
    hMcOut[pid] = 0;
  }

  for (Int_t i = 0; i < 3; ++i) {
    hDCAxy[i] = 0;
  }
}

AliAnalysisTaskDataSpeedOfSound276TeV::AliAnalysisTaskDataSpeedOfSound276TeV(
    const char* name)
    : AliAnalysisTaskSE(name),
      fESD(0x0),
      fAOD(0x0),
      fMC(0x0),
      fMCStack(0x0),
      fMCArray(0x0),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fCentEst("V0M"),
      fAnalysisType("ESD"),
      fAnalysisMC(kFALSE),
      fAnalysisPbPb(kFALSE),
      ftrigBit(0x0),
      fRandom(0x0),
      fPileUpRej(kFALSE),
      fVtxCut(10.0),
      fEtaCut(0.9),
      fEtaMin(-0.8),
      fEtaMax(0.8),
      fPtMin(0.15),
      fEtaGappT(0.4),
      fMinCent(0.0),
      fMaxCent(100.0),
      fStoreMcIn(kFALSE),
      fEtaGapNchMin(0.7),
      fEtaGapNchMax(1.4),
      fMcProcessType(-999),
      fTriggeredEventMB(-999),
      fVtxStatus(-999),
      fZvtx(-999),
      fZvtxMC(-999),
      fRun(-999),
      fEventId(-999),
      fv0mamplitude(0.0),
      fv0mpercentile(0.0),
      fListOfObjects(0),
      fEvents(0x0),
      fVtx(0x0),
      fVtxMC(0x0),
      fVtxBeforeCuts(0x0),
      fVtxAfterCuts(0x0),
      fn1(0x0),
      fcent(0x0),
      hPhi(0x0),
      hPhiEtaSPD(0x0),
      hPhiEtaGapSPD(0x0),
      hVtxZvsTracklets(0x0),
      pV0MAmpChannel(0x0),
      hPtEtaNegvsNchEtaPos(0x0),
      pPtEtaNegvsNchEtaPos(0x0),
      hPtEtaPosvsNchEtaNeg(0x0),
      pPtEtaPosvsNchEtaNeg(0x0),
      pPtvsNch(0x0),
      pPtvsV0MAmp(0x0),
      hPtvsV0MAmp(0x0),
      hPtvsTracklets10(0x0),
      pPtvsTracklets10(0x0),
      hPtvsTracklets14(0x0),
      pPtvsTracklets14(0x0),
      hPtvsTrackletsEtaGap(0x0),
      pPtvsTrackletsEtaGap(0x0),
      hV0Mmult(0x0),
      hV0MAmplitude(0x0),
      hNchvsV0M(0x0),
      hNchvsV0MAmp(0x0),
      hV0MvsV0MAmp(0x0),
      hNchEtaPosvsNchEtaNeg(0x0),
      hTrackletsvsV0MAmp10(0x0),
      hTrackletsvsV0MAmp14(0x0),
      hTrackletsEtaGap(0x0),
      hTrackletsvsV0MAmpEtaGap(0x0),
      hDaDCAxy(0x0)

{
  // Default constructor (should not be used)
  for (Int_t pid = 0; pid < 8; ++pid) {
    hMcIn[pid] = 0;
    hMcOut[pid] = 0;
  }

  for (Int_t i = 0; i < 3; ++i) {
    hDCAxy[i] = 0;
  }
  DefineOutput(1, TList::Class());  // esto es nuevo
}

AliAnalysisTaskDataSpeedOfSound276TeV::
    ~AliAnalysisTaskDataSpeedOfSound276TeV() {
  //
  // Destructor
  //
}

//______________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::UserCreateOutputObjects() {
  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested
  // We also create the random generator here so it might get different seeds...
  fRandom = new TRandom(0);  // 0 means random seed

  // OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  if (!fTrackFilter) {
    // Same as in the Nch vs. mult in pp, p-Pb and Pb-Pb
    fTrackFilter = new AliAnalysisFilter("trackFilter2015");
    AliESDtrackCuts* fCuts2_1 = new AliESDtrackCuts();
    fCuts2_1->SetMaxFractionSharedTPCClusters(0.4);                //
    fCuts2_1->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);  //
    fCuts2_1->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);           //
    fCuts2_1->SetMaxChi2PerClusterTPC(4);                          //
    fCuts2_1->SetAcceptKinkDaughters(kFALSE);                      //
    fCuts2_1->SetRequireTPCRefit(kTRUE);                           //
    fCuts2_1->SetRequireITSRefit(kTRUE);                           //
    fCuts2_1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                       AliESDtrackCuts::kAny);    //
    fCuts2_1->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");  //
    fCuts2_1->SetMaxChi2TPCConstrainedGlobal(36);                 //
    fCuts2_1->SetMaxDCAToVertexZ(2);                              //
    fCuts2_1->SetDCAToVertex2D(kFALSE);                           //
    fCuts2_1->SetRequireSigmaToVertex(kFALSE);                    //
    fCuts2_1->SetMaxChi2PerClusterITS(36);                        //
    fCuts2_1->SetEtaRange(-0.8, 0.8);
    fTrackFilter->AddCuts(fCuts2_1);
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

  //
  // Histograms
  //
  fEvents =
      new TH1I("fEvents", "Number of analyzed events; Events; Counts", 3, 0, 3);
  fListOfObjects->Add(fEvents);

  fn1 = new TH1F("fn1", "fn1", 11, -0.5, 10.5);
  fListOfObjects->Add(fn1);

  fcent = new TH1F("fcent", "fcent", 104, -2, 102);
  fListOfObjects->Add(fcent);

  fVtx = new TH1I("fVtx", "Vtx info (0=no, 1=yes); Vtx; Counts", 2, -0.5, 1.5);
  // fListOfObjects->Add(fVtx);

  fVtxBeforeCuts = new TH1F(
      "fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts",
      120, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);

  fVtxAfterCuts = new TH1F("fVtxAfterCuts",
                           "Vtx distribution (before cuts); Vtx z [cm]; Counts",
                           120, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);

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

  constexpr int v0m_Nbins080{6};
  constexpr double v0m_bins080[v0m_Nbins080 + 1] = {0.0,  1.0,  5.0, 10.0,
                                                    20.0, 50.0, 80.0};

  const Int_t nPtBinsV0s = 25;
  Double_t ptBinsV0s[nPtBinsV0s + 1] = {
      0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,  1.2,  1.4,
      1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 9.0, 12.0, 15.0, 20.0};

  constexpr int dcaxy_Nbins{100};
  double dcaxy_bins[dcaxy_Nbins + 1] = {0};
  for (int i = 0; i <= dcaxy_Nbins; ++i) {
    dcaxy_bins[i] = -3.0 + (0.06 * i);
  }

  const char* Pid[8] = {"Ch",     "Pion",   "Kaon",   "Proton",
                        "Sigmap", "Sigmam", "Lambda", "Rest"};

  const char* ParticleType[3] = {"Pri", "WeDe", "MaIn"};

  hVtxZvsTracklets =
      new TH2F("hVtxZvsTracklets", ";#it{Z}_{vtz} (cm); #it{N}_{tracklets}", 50,
               -10, 10, tracklets_Nbins, tracklets_bins);
  fListOfObjects->Add(hVtxZvsTracklets);

  hPhi =
      new TH2D("histPhi", "pt; #phi'", nPtBinsV0s, ptBinsV0s, 90, -0.05, 0.4);
  // dE/dx vs phi, pions at the MIP
  // fListOfObjects->Add(hPhi);

  hPhiEtaSPD = new TH2F("hPhiEtaSPD", ";#varphi; #eta (|#eta|#leq1.4)", 80, 0,
                        2 * TMath::Pi(), 80, -1.4, 1.4);
  fListOfObjects->Add(hPhiEtaSPD);

  hPhiEtaGapSPD =
      new TH2F("hPhiEtaGapSPD", ";#varphi; #eta (0.7#leq|#eta|#leq1.4)", 80, 0,
               2 * TMath::Pi(), 80, -1.4, 1.4);
  fListOfObjects->Add(hPhiEtaGapSPD);

  pV0MAmpChannel =
      new TProfile("pV0MAmpChannel", ";V0 Channel; Amplitude;", 64, -0.5, 63.5);
  fListOfObjects->Add(pV0MAmpChannel);

  hV0Mmult =
      new TH1F("hV0Mmult", ";V0M (%);Entries", v0m_Nbins080, v0m_bins080);
  fListOfObjects->Add(hV0Mmult);

  hV0MAmplitude =
      new TH1F("hV0MAmp", ";V0M Amplitude; Entries", v0mAmp_Nbins, v0mAmp_bins);
  fListOfObjects->Add(hV0MAmplitude);

  hNchvsV0M = new TH2F("hNchvsV0M", ";#it{N}_{ch} (|#eta|<0.8); V0M (%)",
                       nch_Nbins, nch_bins, v0m_Nbins080, v0m_bins080);
  fListOfObjects->Add(hNchvsV0M);

  hNchvsV0MAmp = new TH2F("hNchvsV0MAmp", ";#it{N}_{ch} (|#eta|<0.8); V0M Amp",
                          nch_Nbins, nch_bins, v0mAmp_Nbins, v0mAmp_bins);
  fListOfObjects->Add(hNchvsV0MAmp);

  hV0MvsV0MAmp = new TH2F("hV0MvsV0MAmp", ";V0M Ampl; V0M (%)", v0mAmp_Nbins,
                          v0mAmp_bins, v0m_Nbins080, v0m_bins080);
  fListOfObjects->Add(hV0MvsV0MAmp);

  hNchEtaPosvsNchEtaNeg =
      new TH2D("hNchEtaPosvsNchEtaNeg",
               ";#it{N}_{ch} (0#leq#eta#leq0.8); #it{N}_{ch} (-0.8#leq#eta<0)",
               nch_Nbins, nch_bins, nch_Nbins, nch_bins);
  fListOfObjects->Add(hNchEtaPosvsNchEtaNeg);

  hPtEtaNegvsNchEtaPos = new TH2D("hPtEtaNegvsNchEtaPos",
                                  "; #it{N}_{ch} (0#leq#eta#leq0.8); "
                                  "#it{p}_{T} (GeV/#it{c}, -0.8#leq#eta<0)",
                                  nch_Nbins, nch_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtEtaNegvsNchEtaPos);

  pPtEtaNegvsNchEtaPos = new TProfile("pPtEtaNegvsNchEtaPos",
                                      "; #it{N}_{ch} (0#leq#eta#leq0.8); "
                                      "#it{p}_{T} (GeV/#it{c}, -0.8#leq#eta<0)",
                                      nch_Nbins, nch_bins);
  fListOfObjects->Add(pPtEtaNegvsNchEtaPos);

  hPtEtaPosvsNchEtaNeg = new TH2D("hPtEtaPosvsNchEtaNeg",
                                  "; #it{N}_{ch} (-0.8#leq#eta<0); #it{p}_{T} "
                                  "(GeV/#it{c}, 0#geq#eta#leq0.8)",
                                  nch_Nbins, nch_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtEtaPosvsNchEtaNeg);

  pPtEtaPosvsNchEtaNeg =
      new TProfile("pPtEtaPosvsNchEtaNeg",
                   "; #it{N}_{ch} (-0.8#leq#eta<0); #it{p}_{T} (GeV/#it{c}, "
                   "0#geq#eta#leq0.8)",
                   nch_Nbins, nch_bins);
  fListOfObjects->Add(pPtEtaPosvsNchEtaNeg);

  pPtvsNch = new TProfile("pPtvsNch",
                          "; #it{N}_{ch}^{rec} (|#eta|<0.8); #LT#it{p}_{T}#GT "
                          "(GeV/#it{c}, |#eta|<0.8)",
                          nch_Nbins, nch_bins);
  fListOfObjects->Add(pPtvsNch);

  hPtvsV0MAmp =
      new TH2D("hPtvsV0MAmp", ";V0M Amp; #it{p}_{T} (GeV/#it{c}, |#eta|<0.8)",
               v0mAmp_Nbins, v0mAmp_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsV0MAmp);

  pPtvsV0MAmp = new TProfile(
      "pPtvsV0MAmp",
      "; V0M Amplitude; #LT#it{p}_{T}#GT (GeV/#it{c}, |#eta|<0.8) ",
      v0mAmp_Nbins, v0mAmp_bins);
  fListOfObjects->Add(pPtvsV0MAmp);

  hTrackletsvsV0MAmp14 = new TH2F(
      "hTrackletsvsV0MAmp14", "; #it{N}_{tracklet} (|#eta|#leq1.4); V0M Amp",
      tracklets_Nbins, tracklets_bins, v0mAmp_Nbins, v0mAmp_bins);
  fListOfObjects->Add(hTrackletsvsV0MAmp14);

  hTrackletsvsV0MAmp10 = new TH2F(
      "hTrackletsvsV0MAmp10", "; #it{N}_{tracklet} (|#eta|#leq1); V0M Amp",
      tracklets_Nbins, tracklets_bins, v0mAmp_Nbins, v0mAmp_bins);
  fListOfObjects->Add(hTrackletsvsV0MAmp10);

  hTrackletsvsV0MAmpEtaGap =
      new TH2F("hTrackletsvsV0MAmpEtaGap",
               "; #it{N}_{tracklet} (0.7<|#eta|#leq1.4); V0M Amp",
               tracklets_Nbins, tracklets_bins, v0mAmp_Nbins, v0mAmp_bins);
  fListOfObjects->Add(hTrackletsvsV0MAmpEtaGap);

  hTrackletsEtaGap = new TH1F(
      "hTrackletsEtaGap", "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); Entries",
      tracklets_Nbins, tracklets_bins);
  fListOfObjects->Add(hTrackletsEtaGap);

  hPtvsTracklets14 =
      new TH2D("hPtvsTracklets14",
               "; #it{N}_{tracklet} (|#eta|#leq1.4); #it{p}_{T} (GeV/#it{c}, "
               "|#eta|<0.8)",
               tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsTracklets14);

  pPtvsTracklets14 = new TProfile("pPtvsTracklets14",
                                  "; #it{N}_{tracklet} (|#eta|#leq1.4); "
                                  "#LT#it{p}_{T}#GT (GeV/#it{c}, |#eta|<0.8)",
                                  tracklets_Nbins, tracklets_bins);
  fListOfObjects->Add(pPtvsTracklets14);

  hPtvsTracklets10 = new TH2D(
      "hPtvsTracklets10",
      "; #it{N}_{tracklet} (|#eta|#leq1); #it{p}_{T} (GeV/#it{c}, |#eta|<0.8)",
      tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsTracklets10);

  pPtvsTracklets10 = new TProfile("pPtvsTracklets10",
                                  "; #it{N}_{tracklet} (|#eta|#leq1); "
                                  "#LT#it{p}_{T}#GT (GeV/#it{c}, |#eta|<0.8)",
                                  tracklets_Nbins, tracklets_bins);
  fListOfObjects->Add(pPtvsTracklets10);

  hPtvsTrackletsEtaGap =
      new TH2D("hPtvsTrackletsEtaGap",
               "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); #it{p}_{T} "
               "(GeV/#it{c}, |#eta|<0.4)",
               tracklets_Nbins, tracklets_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsTrackletsEtaGap);

  pPtvsTrackletsEtaGap =
      new TProfile("pPtvsTrackletsEtaGap",
                   "; #it{N}_{tracklet} (0.7#leq|#eta|#leq1.4); "
                   "#LT#it{p}_{T}#GT (GeV/#it{c}, |#eta|<0.4)",
                   tracklets_Nbins, tracklets_bins);
  fListOfObjects->Add(pPtvsTrackletsEtaGap);

  hDaDCAxy = new TH2D("hDCAxyData_0_5", ";DCAxy (cm); pT", dcaxy_Nbins,
                      dcaxy_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hDaDCAxy);

  if (fAnalysisMC) {
    for (Int_t pid = 0; pid < 8; pid++) {
      hMcIn[pid] =
          new TH1D(Form("hPtInPrim_%s", Pid[pid]),
                   Form("MC in (pid %s)", Pid[pid]), pt_Nbins, pt_bins);
      fListOfObjects->Add(hMcIn[pid]);

      hMcOut[pid] =
          new TH1D(Form("hPtOutPrim_%s", Pid[pid]),
                   Form("MC out (pid %s)", Pid[pid]), pt_Nbins, pt_bins);
      fListOfObjects->Add(hMcOut[pid]);
    }

    for (Int_t i = 0; i < 3; ++i) {
      hDCAxy[i] =
          new TH2D(Form("hDCAxy%s_0_5", ParticleType[i]), ";DCAxy (cm); pT",
                   dcaxy_Nbins, dcaxy_bins, pt_Nbins, pt_bins);
      fListOfObjects->Add(hDCAxy[i]);
    }

    fVtxMC = new TH1F("fVtxMC", "mc vtx", 120, -30, 30);
    fListOfObjects->Add(fVtxMC);
  }

  // Post output data.
  PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::UserExec(Option_t*) {
  // Main loop

  //
  // First we make sure that we have valid input(s)!
  //

  AliVEvent* event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  if (fAnalysisType == "ESD") {
    fESD = dynamic_cast<AliESDEvent*>(event);
    if (!fESD) {
      Printf("%s:%d ESDEvent not found in Input Manager", (char*)__FILE__,
             __LINE__);
      this->Dump();
      return;
    }
  } else {
    fAOD = dynamic_cast<AliAODEvent*>(event);
    if (!fAOD) {
      Printf("%s:%d AODEvent not found in Input Manager", (char*)__FILE__,
             __LINE__);
      this->Dump();
      return;
    }
  }

  if (fAnalysisMC) {
    if (fAnalysisType == "ESD") {
      fMC = dynamic_cast<AliMCEvent*>(MCEvent());
      if (!fMC) {
        Printf("%s:%d MCEvent not found in Input Manager", (char*)__FILE__,
               __LINE__);
        this->Dump();
        return;
      }

      fMCStack = fMC->Stack();

      if (!fMCStack) {
        Printf("%s:%d MCStack not found in Input Manager", (char*)__FILE__,
               __LINE__);
        this->Dump();
        return;
      }
    } else {  // AOD

      fMC = dynamic_cast<AliMCEvent*>(MCEvent());
      if (fMC) fMC->Dump();

      fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
      if (!fMCArray) {
        Printf("%s:%d AOD MC array not found in Input Manager", (char*)__FILE__,
               __LINE__);
        this->Dump();
        return;
      }
    }
  }

  // Get trigger decision
  fTriggeredEventMB = 0;  // init

  fn1->Fill(0);

  if (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()
                                   ->GetInputEventHandler()))
          ->IsEventSelected() &
      ftrigBit) {
    fTriggeredEventMB = 1;  // event triggered as minimum bias
  }

  // Get process type for MC
  fMcProcessType = 0;  // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD

  // real data that are not triggered we skip
  if (!fAnalysisMC && !fTriggeredEventMB) {
    return;
  }

  fn1->Fill(1);

  if (fAnalysisMC) {
    if (fAnalysisType == "ESD") {
      AliHeader* headerMC = fMC->Header();
      if (headerMC) {
        AliGenEventHeader* genHeader = headerMC->GenEventHeader();
        TArrayF vtxMC(3);  // primary vertex  MC
        vtxMC[0] = 9999;
        vtxMC[1] = 9999;
        vtxMC[2] = 9999;  // initialize with dummy
        if (genHeader) {
          genHeader->PrimaryVertex(vtxMC);
        }
        fZvtxMC = vtxMC[2];

        // PYTHIA:
        AliGenPythiaEventHeader* pythiaGenHeader =
            dynamic_cast<AliGenPythiaEventHeader*>(headerMC->GenEventHeader());
        if (pythiaGenHeader) {  // works only for pythia
          fMcProcessType =
              GetPythiaEventProcessType(pythiaGenHeader->ProcessType());
        }
        // PHOJET:
        AliGenDPMjetEventHeader* dpmJetGenHeader =
            dynamic_cast<AliGenDPMjetEventHeader*>(headerMC->GenEventHeader());
        if (dpmJetGenHeader) {
          fMcProcessType =
              GetDPMjetEventProcessType(dpmJetGenHeader->ProcessType());
        }
      }
    } else {  // AOD

      AliAODMCHeader* mcHeader =
          dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject("mcHeader"));

      if (mcHeader) {
        fZvtxMC = mcHeader->GetVtxZ();

        if (strstr(mcHeader->GetGeneratorName(), "Pythia")) {
          fMcProcessType = GetPythiaEventProcessType(mcHeader->GetEventType());
        } else {
          fMcProcessType = GetDPMjetEventProcessType(mcHeader->GetEventType());
        }
      }
    }
  }

  if (fAnalysisType == "ESD") {
    const AliESDVertex* vtxESD = fESD->GetPrimaryVertexTracks();
    // # of tracklets/tracks used for the estimate
    if (vtxESD->GetNContributors() < 1) {
      // SPD vertex
      vtxESD = fESD->GetPrimaryVertexSPD();
      /* quality checks on SPD-vertex */
      TString vertexType = vtxESD->GetTitle();
      if (vertexType.Contains("vertexer: Z") &&
          (vtxESD->GetDispersion() > 0.04 || vtxESD->GetZRes() > 0.25))
        fZvtx = -1599;  // vertex = 0x0; //
      else if (vtxESD->GetNContributors() < 1)
        fZvtx = -999;  // vertex = 0x0; //
      else
        fZvtx = vtxESD->GetZ();
    } else {
      fZvtx = vtxESD->GetZ();
    }
  } else  // AOD
  {
    fZvtx = GetVertex(fAOD);
  }

  fVtxBeforeCuts->Fill(fZvtx);

  // cut on the z position of vertex
  if (TMath::Abs(fZvtx) > fVtxCut) {
    return;
  }
  fn1->Fill(2);

  Float_t centrality = -10.0;

  // only analyze triggered events
  if (fTriggeredEventMB) {
    if (fAnalysisType == "ESD") {
      if (fAnalysisPbPb) {
        AliCentrality* centObject = fESD->GetCentrality();
        centrality = centObject->GetCentralityPercentile(fCentEst);
        if ((centrality > fMaxCent) || (centrality < fMinCent)) {
          return;
        }
        fv0mpercentile = centrality;
      }
      fcent->Fill(centrality);
      fn1->Fill(3);
      if (fAnalysisMC) {
        if (TMath::Abs(fZvtxMC) < fVtxCut) {
          ProcessMCTruthESD();
          fVtxMC->Fill(fZvtxMC);
        }
      }
      AnalyzeESD(fESD);
    } else {  // AOD
      if (fAnalysisPbPb) {
        AliCentrality* centObject = fAOD->GetCentrality();
        if (centObject) {
          centrality = centObject->GetCentralityPercentile(fCentEst);
        }
        // cout<<"centrality="<<centrality<<endl;
        if ((centrality > fMaxCent) || (centrality < fMinCent)) return;
      }
      fcent->Fill(centrality);
      fn1->Fill(3);
      if (fAnalysisMC) {
        if (TMath::Abs(fZvtxMC) < fVtxCut) {
          // ProcessMCTruthAOD();
          fVtxMC->Fill(fZvtxMC);
        }
      }
      AnalyzeAOD(fAOD);
    }
  }

  fVtxAfterCuts->Fill(fZvtx);

  // Post output data.
  PostData(1, fListOfObjects);
}

//________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::AnalyzeESD(AliESDEvent* esdEvent) {
  fRun = esdEvent->GetRunNumber();
  fEventId = 0;
  if (esdEvent->GetHeader()) {
    fEventId = GetEventIdAsLong(esdEvent->GetHeader());
  }

  Bool_t isPileup = esdEvent->IsPileupFromSPD();
  if (fPileUpRej) {
    if (isPileup) {
      return;
    }
    fn1->Fill(4);
  }

  //  Int_t     event     = esdEvent->GetEventNumberInFile();
  // UInt_t    time      = esdEvent->GetTimeStamp();
  //  ULong64_t trigger   = esdEvent->GetTriggerMask();
  // magf = esdEvent->GetMagneticField();

  if (fTriggeredEventMB) {  // Only MC case can we have not triggered events

    // accepted event
    fEvents->Fill(0);

    // Change, 10/04/13. Now accept all events to do a correct normalization
    // if(fVtxStatus!=1) return; // accepted vertex
    // Int_t nESDTracks = esdEvent->GetNumberOfTracks();

    GetCalibratedV0Amplitude();

    //! Get SPD tracklets multiplicity
    GetSPDMultiplicity();

    //! pT spectra
    MultiplicityDistributions();

    DCAxyDistributions();

    if (fAnalysisMC) {
      TrackingEfficiency();
    }

    // ProduceArrayTrksESD(esdEvent);  // produce array with global track
    // parameters

    fEvents->Fill(1);

  }  // end if triggered
}
//________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::MultiplicityDistributions() {
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
    hTrackletsvsV0MAmpEtaGap->Fill(fTrackletsEtaGap, fv0mamplitude);
  }
}
//________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::GetSPDMultiplicity() {
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
      hPhiEtaGapSPD->Fill(phi, eta);
    }

    if (TMath::Abs(eta) <= 1.4) {
      fTracklets14++;
      hPhiEtaSPD->Fill(phi, eta);
    }
  }

  hVtxZvsTracklets->Fill(spdVtxZ, fTracklets14);
}
//________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::GetCalibratedV0Amplitude() {
  float mV0M{0.0};
  for (int i = 0; i < 64; i++) {
    mV0M += fESD->GetVZEROEqMultiplicity(i);
    pV0MAmpChannel->Fill(i, fESD->GetVZEROEqMultiplicity(i));
  }
  fv0mamplitude = mV0M;
}
//________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::AnalyzeAOD(AliAODEvent* aodEvent) {
  fRun = aodEvent->GetRunNumber();
  fEventId = 0;
  if (aodEvent->GetHeader()) fEventId = GetEventIdAsLong(aodEvent->GetHeader());

  // UInt_t    time      = 0; // Missing AOD info? aodEvent->GetTimeStamp();
  // magf = aodEvent->GetMagneticField();

  // Int_t     trackmult = 0; // no pt cuts
  // Int_t     nadded    = 0;

  Bool_t isPileup = aodEvent->IsPileupFromSPD();
  if (fPileUpRej)
    if (isPileup) return;
  fn1->Fill(4);

  if (fTriggeredEventMB) {  // Only MC case can we have not triggered events

    // accepted event
    fEvents->Fill(0);

    // if(fVtxStatus!=1) return; // accepted vertex
    // Int_t nAODTracks = aodEvent->GetNumberOfTracks();

    // ProduceArrayTrksAOD(aodEvent);

    fEvents->Fill(1);

  }  // end if triggered
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskDataSpeedOfSound276TeV::GetVertex(
    const AliVEvent* event) const {
  Float_t zvtx = -999;

  const AliVVertex* primaryVertex = event->GetPrimaryVertex();

  if (primaryVertex->GetNContributors() > 0) zvtx = primaryVertex->GetZ();

  return zvtx;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskDataSpeedOfSound276TeV::GetPidCode(Int_t pdgCode) const {
  // return our internal code for pions, kaons, and protons

  Short_t pidCode = 6;

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
      pidCode = 7;  // rest
  };

  return pidCode;
}
//_____________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::ProcessMCTruthESD() {
  // Fill the special MC histogram with the MC truth info

  if (fv0mpercentile < 0.0 || fv0mpercentile > 5.0) {
    return;
  }

  const Int_t nTracksMC = fMCStack->GetNtrack();
  for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
    // Cuts
    if (!(fMCStack->IsPhysicalPrimary(iTracks))) continue;

    TParticle* trackMC = fMCStack->Particle(iTracks);

    TParticlePDG* pdgPart = trackMC->GetPDG();
    Double_t chargeMC = pdgPart->Charge();

    if (TMath::Abs(trackMC->Eta()) > fEtaCut) continue;

    Short_t pidCodeMC = 0;
    pidCodeMC = GetPidCode(trackMC->GetPdgCode());

    // lambda
    if (TMath::Abs(trackMC->GetPdgCode()) == 3122) {
      hMcIn[6]->Fill(trackMC->Pt());
    }

    if (chargeMC == 0) continue;

    hMcIn[0]->Fill(trackMC->Pt());
    hMcIn[pidCodeMC]->Fill(trackMC->Pt());

  }  // MC track loop
}
//_____________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::TrackingEfficiency() {
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

    TParticle* particle = fMCStack->Particle(label);
    if (!particle) {
      continue;
    }
    if (track->Charge() == 0) {
      continue;
    }
    if (!fMCStack->IsPhysicalPrimary(label)) {
      continue;
    }

    hMcOut[0]->Fill(track->Pt());
    Int_t pidCode = GetPidCode(particle->GetPdgCode());

    hMcOut[pidCode]->Fill(track->Pt());
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskDataSpeedOfSound276TeV::DCAxyDistributions() {
  if (fv0mpercentile < 0.0 || fv0mpercentile > 5.0) {
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

    hDaDCAxy->Fill(dcaxy, track->Pt());
  }

  if (fAnalysisMC) {
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
//_____________________________________________________________________________
Short_t AliAnalysisTaskDataSpeedOfSound276TeV::GetPythiaEventProcessType(
    Int_t pythiaType) {
  //
  // Get the process type of the event.  PYTHIA
  //
  // source PWG0   dNdpt

  Short_t globalType = -1;  // init

  if (pythiaType == 92 || pythiaType == 93) {
    globalType = 2;  // single diffractive
  } else if (pythiaType == 94) {
    globalType = 3;  // double diffractive
  }
  // else if(pythiaType != 91){ // also exclude elastic to be sure... CKB??
  else {
    globalType = 1;  // non diffractive
  }
  return globalType;
}

//_____________________________________________________________________________
Short_t AliAnalysisTaskDataSpeedOfSound276TeV::GetDPMjetEventProcessType(
    Int_t dpmJetType) {
  //
  // get the process type of the event.  PHOJET
  //
  // source PWG0   dNdpt
  // can only read pythia headers, either directly or from cocktalil header
  Short_t globalType = -1;

  if (dpmJetType == 1 ||
      dpmJetType == 4) {  // explicitly inelastic plus central diffraction
    globalType = 1;
  } else if (dpmJetType == 5 || dpmJetType == 6) {
    globalType = 2;
  } else if (dpmJetType == 7) {
    globalType = 3;
  }
  return globalType;
}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskDataSpeedOfSound276TeV::GetEventIdAsLong(
    AliVHeader* header) const {
  // To have a unique id for each event in a run!
  // Modified from AliRawReader.h
  return ((ULong64_t)header->GetBunchCrossNumber() +
          (ULong64_t)header->GetOrbitNumber() * 3564 +
          (ULong64_t)header->GetPeriodNumber() * 16777215 * 3564);
}

//____________________________________________________________________
TParticle* AliAnalysisTaskDataSpeedOfSound276TeV::FindPrimaryMother(
    AliStack* stack, Int_t label) {
  //
  // Finds the first mother among the primary particles of the particle
  // identified by <label>, i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  Int_t motherLabel = FindPrimaryMotherLabel(stack, label);
  if (motherLabel < 0) return 0;

  return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskDataSpeedOfSound276TeV::FindPrimaryMotherLabel(
    AliStack* stack, Int_t label) {
  //
  // Finds the first mother among the primary particles of the particle
  // identified by <label>, i.e. the primary that "caused" this particle
  //
  // returns its label
  //
  // Taken from AliPWG0Helper class
  //
  const Int_t nPrim = stack->GetNprimary();

  while (label >= nPrim) {
    // printf("Particle %d (pdg %d) is not a primary. Let's check its mother
    // %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

    TParticle* particle = stack->Particle(label);
    if (!particle) {
      AliDebugGeneral(
          "FindPrimaryMotherLabel", AliLog::kError,
          Form("UNEXPECTED: particle with label %d not found in stack.",
               label));
      return -1;
    }

    // find mother
    if (particle->GetMother(0) < 0) {
      AliDebugGeneral(
          "FindPrimaryMotherLabel", AliLog::kError,
          Form("UNEXPECTED: Could not find mother of secondary particle %d.",
               label));
      return -1;
    }

    label = particle->GetMother(0);
  }

  return label;
}

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskDataSpeedOfSound276TeV::FindPrimaryMotherAOD(
    AliAODMCParticle* startParticle) {
  //
  // Finds the first mother among the primary particles of the particle
  // identified by <label>, i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  AliAODMCParticle* mcPart = startParticle;

  while (mcPart) {
    if (mcPart->IsPrimary()) return mcPart;

    Int_t mother = mcPart->GetMother();

    mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
  }

  return 0;
}
//____________________________________________________________________
TParticle* AliAnalysisTaskDataSpeedOfSound276TeV::FindPrimaryMotherV0(
    AliStack* stack, Int_t label) {
  //
  // Finds the first mother among the primary particles of the particle
  // identified by <label>, i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  Int_t nSteps = 0;

  Int_t motherLabel = FindPrimaryMotherLabelV0(stack, label, nSteps);
  if (motherLabel < 0) return 0;

  return stack->Particle(motherLabel);
}

//____________________________________________________________________
Int_t AliAnalysisTaskDataSpeedOfSound276TeV::FindPrimaryMotherLabelV0(
    AliStack* stack, Int_t label, Int_t& nSteps) {
  //
  // Finds the first mother among the primary particles of the particle
  // identified by <label>, i.e. the primary that "caused" this particle
  //
  // returns its label
  //
  // Taken from AliPWG0Helper class
  //
  nSteps = 0;
  const Int_t nPrim = stack->GetNprimary();

  while (label >= nPrim) {
    // printf("Particle %d (pdg %d) is not a primary. Let's check its mother
    // %d\n", label, mother->GetPdgCode(), mother->GetMother(0));

    nSteps++;  // 1 level down

    TParticle* particle = stack->Particle(label);
    if (!particle) {
      AliDebugGeneral(
          "FindPrimaryMotherLabelV0", AliLog::kError,
          Form("UNEXPECTED: particle with label %d not found in stack.",
               label));
      return -1;
    }

    // find mother
    if (particle->GetMother(0) < 0) {
      AliDebugGeneral(
          "FindPrimaryMotherLabelV0", AliLog::kError,
          Form("UNEXPECTED: Could not find mother of secondary particle %d.",
               label));
      return -1;
    }

    label = particle->GetMother(0);
  }

  return label;
}

//____________________________________________________________________
AliAODMCParticle* AliAnalysisTaskDataSpeedOfSound276TeV::FindPrimaryMotherAODV0(
    AliAODMCParticle* startParticle, Int_t& nSteps) {
  //
  // Finds the first mother among the primary particles of the particle
  // identified by <label>, i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  nSteps = 0;

  AliAODMCParticle* mcPart = startParticle;

  while (mcPart) {
    if (mcPart->IsPrimary()) return mcPart;

    Int_t mother = mcPart->GetMother();

    mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
    nSteps++;  // 1 level down
  }

  return 0;
}
//__________________________________________________________________________
Bool_t AliAnalysisTaskDataSpeedOfSound276TeV::PhiCut(Double_t pt, Double_t phi,
                                                     Double_t q, Float_t mag,
                                                     TF1* phiCutLow,
                                                     TF1* phiCutHigh) {
  if (pt < 2.0) return kTRUE;

  // TF1* cutLow = new TF1("StandardPhiCutLow", "0.1/x/x+pi/18.0-0.025", 0, 50);
  // TF1* cutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0,
  // 50);

  // Double_t phi = track->Phi();
  if (mag < 0)  // for negatve polarity field
    phi = TMath::TwoPi() - phi;
  if (q < 0)  // for negatve charge
    phi = TMath::TwoPi() - phi;

  phi += TMath::Pi() / 18.0;  // to center gap in the middle
  phi = fmod(phi, TMath::Pi() / 9.0);

  if (phi < phiCutHigh->Eval(pt) && phi > phiCutLow->Eval(pt))
    return kFALSE;  // reject track

  hPhi->Fill(pt, phi);

  return kTRUE;
}
