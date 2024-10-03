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

#include "AliAnalysisTaskZNZP.h"

class AliAnalysisTaskZNZP;  // your analysis class

ClassImp(AliAnalysisTaskZNZP)  // classimp: necessary for root

    AliAnalysisTaskZNZP::AliAnalysisTaskZNZP()
    : AliAnalysisTaskSE(),
      fESD(0),
      fEventCuts(0x0),
      fMCStack(0),
      fMC(0),
      fUseMC(kFALSE),
      fTowerEnergy(true),
      fTrigger(AliVEvent::kCentral),
      fMultSelection(0x0),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fOutputList(0),
      fEtaCut(0.8),
      fPtMin(0.15),
      fV0Mmin(0.0),
      fV0Mmax(80.0),
      ftrackmult08(0),
      fv0mpercentile(0),
      fv0mamplitude(0),
      fSPD(0),
      fNchTPC(0),
      fET(0.),
      fZNCvsV0M(0),
      fZNAvsV0M(0),
      pV0MAmpChannel(0),
      hV0Percentile(0),
      hBestVtxZ(0),
      hZNAvsV0M(0),
      hZNCvsV0M(0),
      hAsyN(0),
      hZPAvsV0M(0),
      hZPCvsV0M(0),
      hZNCNorm(0),
      hZNANorm(0),
      pZNChannel(0),
      pZPChannel(0),
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
      hZNCTowvsNCEn(0),
      hZPATowvsPAEn(0) {}
//_____________________________________________________________________________
AliAnalysisTaskZNZP::AliAnalysisTaskZNZP(const char* name)
    : AliAnalysisTaskSE(name),
      fESD(0),
      fEventCuts(0x0),
      fMCStack(0),
      fMC(0),
      fUseMC(kFALSE),
      fTowerEnergy(true),
      fTrigger(AliVEvent::kCentral),
      fMultSelection(0x0),
      fTrackFilter(0x0),
      fTrackFilterwoDCA(0x0),
      fOutputList(0),
      fEtaCut(0.8),
      fPtMin(0.15),
      fV0Mmin(0.0),
      fV0Mmax(80.0),
      ftrackmult08(0),
      fv0mpercentile(0),
      fv0mamplitude(0),
      fSPD(0),
      fNchTPC(0),
      fET(0.),
      fZNCvsV0M(0),
      fZNAvsV0M(0),
      pV0MAmpChannel(0),
      hV0Percentile(0),
      hBestVtxZ(0),
      hZNAvsV0M(0),
      hZNCvsV0M(0),
      hAsyN(0),
      hZPAvsV0M(0),
      hZPCvsV0M(0),
      hZNCNorm(0),
      hZNANorm(0),
      pZNChannel(0),
      pZPChannel(0),
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
      hZNCTowvsNCEn(0),
      hZPATowvsPAEn(0) {
  DefineInput(0, TChain::Class());  // define the input of the analysis: in this
                                    // case you take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about
  // it, does its work automatically
  DefineOutput(1, TList::Class());  // define the ouptut of the analysis: in
                                    // this case it's a list of histograms
}
//_____________________________________________________________________________
AliAnalysisTaskZNZP::~AliAnalysisTaskZNZP() {
  // destructor
  if (fOutputList) {
    delete fOutputList;  // at the end of your task, it is deleted from memory
                         // by calling this function
    fOutputList = 0x0;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskZNZP::UserCreateOutputObjects() {
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

  //! <ZNC> and <ZNA> obtained from the run 296621
  fZNCvsV0M = new TF1("fZNCvsV0M", "pol4", 0., 100.);
  fZNCvsV0M->SetParameter(0, 8.09862);
  fZNCvsV0M->SetParameter(1, 2.36909);
  fZNCvsV0M->SetParameter(2, -0.0537093);
  fZNCvsV0M->SetParameter(3, 0.000549338);
  fZNCvsV0M->SetParameter(4, -2.60388e-06);

  fZNAvsV0M = new TF1("fZNAvsV0M", "pol4", 0., 100.);
  fZNAvsV0M->SetParameter(0, 7.10818);
  fZNAvsV0M->SetParameter(1, 2.20761);
  fZNAvsV0M->SetParameter(2, -0.0492212);
  fZNAvsV0M->SetParameter(3, 0.000514146);
  fZNAvsV0M->SetParameter(4, -2.56141e-06);

  // create output objects
  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  constexpr double v0mAmp_width{25.0};
  constexpr int v0mAmp_Nbins{1720};
  double v0mAmp_bins[v0mAmp_Nbins + 1] = {0};
  for (int i = 0; i <= v0mAmp_Nbins; ++i) {
    v0mAmp_bins[i] = 0.0 + i * v0mAmp_width;
  }

  // SPD Tracklets (|eta|<0.8)
  constexpr int SPD0p8_Nbins{1650};
  double SPD0p8_bins[SPD0p8_Nbins + 1] = {0};
  for (int i = 0; i <= SPD0p8_Nbins; ++i) {
    SPD0p8_bins[i] = 0.5 + (2.0 * i);
  }

  // Nch (|eta|<0.8)
  constexpr int nch_Nbins{1250};
  double nch_bins[nch_Nbins + 1] = {0};
  for (int i = 0; i <= nch_Nbins; ++i) {
    nch_bins[i] = 0.5 + (2.0 * i);
  }

  constexpr int Et_Nbins{1000};
  double Et_bins[Et_Nbins + 1] = {0};
  for (int i = 0; i <= Et_Nbins; ++i) {
    Et_bins[i] = 0.0 + (2.0 * i);
  }

  const int nBinsV0M090{80};
  double BinsV0M090[nBinsV0M090 + 1] = {0.0};
  for (int i = 0; i <= nBinsV0M090; ++i) {
    BinsV0M090[i] = 0.0 + (double)i;
  }

  hV0Percentile = new TH1F("hV0M", ";V0M (%);Entries", nBinsV0M090, BinsV0M090);
  pV0MAmpChannel =
      new TProfile("pV0Channel", ";Channel; Amplitude;", 64, -0.5, 63.5);
  hBestVtxZ =
      new TH1F("hBestVtxZ", ";Vertex_{#it{z}} (cm); Counts;", 400, -11, 11);

  std::string title{"Common PMT:TowerEnergy[0]"};
  if (!fTowerEnergy) title = "All five PMTs: fZDCEnergy";

  hZNAvsV0M = new TH2F(
      "hZNAvsV0M", Form("%s;V0M Per;#it{E}_{ZNA} [TeV]/2.511;", title.c_str()),
      nBinsV0M090, BinsV0M090, 100, 0.0, 100.0);
  hZNCvsV0M = new TH2F(
      "hZNCvsV0M", Form("%s;V0M Per;#it{E}_{ZNC} [TeV]/2.511;", title.c_str()),
      nBinsV0M090, BinsV0M090, 100, 0.0, 100.0);
  hZPAvsV0M = new TH2F(
      "hZPAvsV0M", Form("%s;V0M Per;#it{E}_{ZPA} [TeV]/2.511;", title.c_str()),
      nBinsV0M090, BinsV0M090, 30, 0.0, 30.0);
  hZPCvsV0M = new TH2F(
      "hZPCvsV0M", Form("%s;V0M Per;#it{E}_{ZPC} [TeV]/2.511;", title.c_str()),
      nBinsV0M090, BinsV0M090, 30, 0.0, 30.0);
  hAsyN =
      new TH2F("hAsyN", "Neutron asymmetry;V0M Per; N_{C}-N_{A}/N_{C}+N_{A};",
               nBinsV0M090, BinsV0M090, 50, -1.0, 1.0);

  hZNCNorm = new TH2F("hZNCNorm", ";V0M;<ZNC>/<ZNC>;", nBinsV0M090, BinsV0M090,
                      40, 0., 2.);

  hZNANorm = new TH2F("hZNANorm", ";V0M;<ZNA>/<ZNA>;", nBinsV0M090, BinsV0M090,
                      40, 0., 2.);
  hZNCTowvsNCEn =
      new TH2F("hZNCTowvsNCEn",
               ";fZDCN1Energy [TeV/2.511];fZN1TowerEnergy[0] [TeV/2.511];", 100,
               0., 100., 100, 0., 100.);
  hZPATowvsPAEn =
      new TH2F("hZPATowvsPAEn",
               ";fZDCN2Energy [TeV/2.511];fZN2TowerEnergy[0] [TeV/2.511];", 30,
               0., 30., 30, 0., 30.);
  pZNChannel =
      new TProfile("pZNChannel", ";Channel; <ZN> [TeV/2.511];", 10, 0., 10.);
  pZPChannel =
      new TProfile("pZPChannel", ";Channel; <ZP> [TeV/2.511];", 10, 0., 10.);

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
  fOutputList->Add(hV0Percentile);
  fOutputList->Add(pZNChannel);
  fOutputList->Add(hZNAvsV0M);
  fOutputList->Add(hZNCvsV0M);
  fOutputList->Add(hZNCTowvsNCEn);
  fOutputList->Add(hZNCNorm);
  fOutputList->Add(hZNANorm);
  fOutputList->Add(pZNCvsV0Amp);
  fOutputList->Add(pZNAvsV0Amp);
  fOutputList->Add(pZNCvsNch);
  fOutputList->Add(pZNAvsNch);
  fOutputList->Add(pZNCvsEt);
  fOutputList->Add(pZNAvsEt);
  fOutputList->Add(pZNCvsSPD);
  fOutputList->Add(pZNAvsSPD);
  fOutputList->Add(hAsyN);
  fOutputList->Add(pZPChannel);
  fOutputList->Add(hZPAvsV0M);
  fOutputList->Add(hZPCvsV0M);
  fOutputList->Add(hZPATowvsPAEn);
  fOutputList->Add(pZPCvsV0Amp);
  fOutputList->Add(pZPAvsV0Amp);
  fOutputList->Add(pZPCvsNch);
  fOutputList->Add(pZPAvsNch);
  fOutputList->Add(pZPCvsEt);
  fOutputList->Add(pZPAvsEt);
  fOutputList->Add(pZPCvsSPD);
  fOutputList->Add(pZPAvsSPD);

  fEventCuts.AddQAplotsToList(fOutputList);
  PostData(1, fOutputList);  // postdata will notify the analysis manager of
                             // changes / updates to the
}
//_____________________________________________________________________________
void AliAnalysisTaskZNZP::UserExec(Option_t*) {
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

  SPDActivity();

  //! Nch & Et with TPC
  TPCActivity();

  //! Get ZDC Centrality
  GetZDC();

  PostData(1, fOutputList);
}
//____________________________________________________________
void AliAnalysisTaskZNZP::Terminate(Option_t*) {}
//____________________________________________________________
void AliAnalysisTaskZNZP::VertexPosition() {
  //! best primary vertex available
  const AliVVertex* vtx = fEventCuts.GetPrimaryVertex();
  hBestVtxZ->Fill(vtx->GetZ());
}
//____________________________________________________________
void AliAnalysisTaskZNZP::GetZDC() {
  AliESDZDC* esdZDC = fESD->GetESDZDC();
  if (!esdZDC) {
    return;
  }

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

  hZNCvsV0M->Fill(fv0mpercentile, znc / eA);
  hZNAvsV0M->Fill(fv0mpercentile, zna / eA);
  hZPCvsV0M->Fill(fv0mpercentile, zpc / eA);
  hZPAvsV0M->Fill(fv0mpercentile, zpa / eA);

  pZNCvsV0Amp->Fill(fv0mamplitude, znc);
  pZNAvsV0Amp->Fill(fv0mamplitude, zna);
  pZPCvsV0Amp->Fill(fv0mamplitude, zpc);
  pZPAvsV0Amp->Fill(fv0mamplitude, zpa);

  pZNCvsSPD->Fill(fSPD, znc);
  pZNAvsSPD->Fill(fSPD, zna);
  pZPCvsSPD->Fill(fSPD, zpc);
  pZPAvsSPD->Fill(fSPD, zpa);

  pZNCvsNch->Fill(fNchTPC, znc);
  pZNAvsNch->Fill(fNchTPC, zna);
  pZPCvsNch->Fill(fNchTPC, zpc);
  pZPAvsNch->Fill(fNchTPC, zpa);

  pZNCvsEt->Fill(fET, znc);
  pZNAvsEt->Fill(fET, zna);
  pZPCvsEt->Fill(fET, zpc);
  pZPAvsEt->Fill(fET, zpa);

  hZNCTowvsNCEn->Fill(esdZDC->GetZDCN1Energy() * (gev2tev / eA),
                      towZNC[0] * (gev2tev / eA));
  hZPATowvsPAEn->Fill(esdZDC->GetZDCP2Energy() * (gev2tev / eA),
                      towZPA[0] * (gev2tev / eA));

  for (int i = 0; i < 5; ++i) {
    pZNChannel->Fill(i + 0.5, towZNC[i] * gev2tev / eA);
    pZNChannel->Fill(i + 5.5, towZNA[i] * gev2tev / eA);
    pZPChannel->Fill(i + 0.5, towZPC[i] * gev2tev / eA);
    pZPChannel->Fill(i + 5.5, towZPA[i] * gev2tev / eA);
  }

  //! Asymmetry
  if (!(fv0mpercentile > 0. && fv0mpercentile < 80.)) {
    return;
  }

  //! meanZNC = <EZC/eA> = <EZC>/2.511
  double meanZNC{fZNCvsV0M->Eval(fv0mpercentile)};
  double meanZNA{fZNAvsV0M->Eval(fv0mpercentile)};
  hZNCNorm->Fill(fv0mpercentile, (znc / eA) / meanZNC);
  hZNANorm->Fill(fv0mpercentile, (zna / eA) / meanZNA);
  const double znc_sca{(znc / eA) / meanZNC};
  const double zna_sca{(zna / eA) / meanZNA};
  hAsyN->Fill(fv0mpercentile, (znc_sca - zna_sca) / (znc_sca + zna_sca));
}
//____________________________________________________________
void AliAnalysisTaskZNZP::SPDActivity() {
  fSPD = 0;

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
  if (TMath::Abs(spdVtx->GetZ()) > 10.) {
    return;
  }
  for (auto it = 0; it < SPDptr->GetNumberOfTracklets(); it++) {
    if (TMath::Abs(SPDptr->GetEta(it)) > fEtaCut) {
      continue;
    }

    fSPD++;
  }
}
//____________________________________________________________
void AliAnalysisTaskZNZP::GetCalibratedV0Amplitude() {
  float mV0M{0.0};
  for (int i = 0; i < 64; i++) {
    mV0M += fESD->GetVZEROEqMultiplicity(i);
    pV0MAmpChannel->Fill(i, fESD->GetVZEROEqMultiplicity(i));
  }
  fv0mamplitude = mV0M;
  hV0Percentile->Fill(fv0mpercentile);
  /*hV0MAmplitude->Fill(fv0mamplitude);*/
}
//____________________________________________________________
void AliAnalysisTaskZNZP::TPCActivity() {
  const double masspi{0.13957};
  const int n_tracks{fESD->GetNumberOfTracks()};

  fNchTPC = 0;
  fET = 0.;

  for (int i = 0; i < n_tracks; ++i) {
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(i));
    if (!track) {
      continue;
    }
    if (!fTrackFilter->IsSelected(track)) {
      continue;
    }
    if (track->Pt() < 0.15 || track->Pt() > 50.) {
      continue;
    }
    if (track->Charge() == 0) {
      continue;
    }
    if (TMath::Abs(track->Eta()) > fEtaCut) {
      continue;
    }

    double pt{track->Pt()};
    fNchTPC++;
    fET += TMath::Sqrt(TMath::Power(pt, 2.0) + TMath::Power(masspi, 2.0));
  }
}
//____________________________________________________________
void AliAnalysisTaskZNZP::DCAxyDistributions() const {
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
  }
}

//____________________________________________________________

Bool_t AliAnalysisTaskZNZP::HasRecVertex() {
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
/*bool AliAnalysisTaskZNZP::MeanSigmaZN(double& mean, double& sigma,*/
/*                                      const std::string& ZN) {*/
/*  const std::vector<double> meanZNC{*/
/*      {9.41689, 11.4554, 13.6935, 15.899,  17.9693, 19.8913, 21.8775, 23.6385,*/
/*       25.2471, 26.8672, 28.3384, 29.8036, 31.0197, 32.3755, 33.5516, 34.5927,*/
/*       35.6032, 36.6078, 37.5322, 38.3336, 39.1421, 39.8723, 40.5853, 41.1441,*/
/*       41.7837, 42.2702, 42.8783, 43.321,  43.7049, 44.1573, 44.5504, 44.8895,*/
/*       45.1275, 45.3678, 45.7246, 45.7786, 45.9834, 46.1794, 46.2849, 46.2833,*/
/*       46.3337, 46.4907, 46.4214, 46.344,  46.3888, 46.2107, 46.1971, 46.0917,*/
/*       45.893,  45.6904, 45.4781, 45.4185, 45.0878, 44.8286, 44.5048, 44.2761,*/
/*       43.9231, 43.421,  43.1319, 42.6159, 42.1982, 41.7332, 41.1261, 40.8089,*/
/*       40.1382, 39.5972, 39.127,  38.4715, 37.7199, 36.9609}};*/
/**/
/*  const std::vector<double> sigmaZNC{*/
/*      {3.50677, 3.82505, 4.11002, 4.24869, 4.46803, 4.59024, 4.90368, 5.07824,*/
/*       5.17144, 5.19741, 5.30018, 5.51713, 5.46297, 5.66198, 5.76174, 5.80704,*/
/*       5.83463, 5.97203, 6.03367, 6.11339, 6.12927, 6.17397, 6.24563, 6.32831,*/
/*       6.45077, 6.45699, 6.46789, 6.52385, 6.57878, 6.67186, 6.77811, 6.72446,*/
/*       6.85183, 6.85365, 6.96893, 6.95236, 6.9604,  7.09092, 7.10086, 7.13674,*/
/*       7.2995,  7.24527, 7.29297, 7.38196, 7.45292, 7.48775, 7.59025, 7.67292,*/
/*       7.73854, 7.77752, 7.84931, 8.00492, 7.91375, 8.13183, 8.19829, 8.33288,*/
/*       8.42185, 8.57176, 8.7218,  8.87617, 9.03015, 9.14355, 9.28268, 9.51969,*/
/*       9.68484, 9.92729, 10.0658, 10.3063, 10.4252, 10.6264}};*/
/**/
/*  const std::vector<double> meanZNA{*/
/*      {7.64853, 9.31038, 11.1556, 13.0141, 14.6987, 16.4405, 17.8888, 19.4125,*/
/*       20.8443, 22.1771, 23.4108, 24.5767, 25.6076, 26.716,  27.7061, 28.6364,*/
/*       29.557,  30.2871, 31.0859, 31.7577, 32.4302, 33.0959, 33.7153, 34.2552,*/
/*       34.8104, 35.1898, 35.7727, 36.1279, 36.4649, 36.8997, 37.2122, 37.3977,*/
/*       37.6971, 37.9211, 38.1646, 38.2855, 38.5034, 38.6686, 38.7706, 38.8934,*/
/*       38.9517, 38.975,  39.0276, 39.0661, 38.9757, 39.0779, 38.8923, 38.8756,*/
/*       38.7662, 38.6167, 38.527,  38.2397, 38.1209, 37.9127, 37.763,  37.4913,*/
/*       37.2047, 36.9349, 36.5257, 36.207,  35.8965, 35.5434, 35.0563, 34.6988,*/
/*       34.1538, 33.8234, 33.3074, 32.82,   32.1863, 31.6347}};*/
/**/
/*  const std::vector<double> sigmaZNA{*/
/*      {3.10314, 3.13242, 3.48659, 3.77875, 3.7878,  3.915,   4.03227, 4.20045,*/
/*       4.31136, 4.40799, 4.56356, 4.61643, 4.65824, 4.74357, 4.78761, 4.92246,*/
/*       4.99679, 5.05502, 5.0036,  5.12261, 5.17488, 5.19715, 5.26379, 5.20759,*/
/*       5.38565, 5.41547, 5.42096, 5.47782, 5.47648, 5.65329, 5.61113, 5.61817,*/
/*       5.70361, 5.62584, 5.73757, 5.73626, 5.75184, 5.85641, 5.92166, 5.9685,*/
/*       5.99215, 6.08534, 6.02079, 6.1867,  6.17472, 6.26609, 6.2535,  6.34106,*/
/*       6.32414, 6.48443, 6.57763, 6.57414, 6.70035, 6.78231, 6.80588, 6.87396,*/
/*       7.06222, 7.13683, 7.23114, 7.44464, 7.49467, 7.61375, 7.7578,  7.96836,*/
/*       8.12227, 8.33368, 8.3856,  8.5707,  8.92543, 8.95762}};*/
/**/
/*  bool goodCent{false};*/
/*  if (fv0mpercentile >= 0. && fv0mpercentile < 70.) goodCent = true;*/
/**/
/*  int bin{CentBin()};*/
/*  if (!goodCent || bin < 0) return false;*/
/**/
/*  if (ZN == "ZNC") {*/
/*    mean = meanZNC.at(bin);*/
/*    sigma = sigmaZNC.at(bin);*/
/*  }*/
/*  if (ZN == "ZNA") {*/
/*    mean = meanZNA.at(bin);*/
/*    sigma = sigmaZNA.at(bin);*/
/*  }*/
/**/
/*  return goodCent;*/
/*}*/
//____________________________________________________________
/*int AliAnalysisTaskZNZP::CentBin() {*/
/*  const std::vector<double> low{*/
/*      {0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12., 13.,*/
/*       14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27.,*/
/*       28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41.,*/
/*       42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55.,*/
/*       56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69.}};*/
/*  const std::vector<double> high{*/
/*      {1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12., 13., 14.,*/
/*       15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28.,*/
/*       29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42.,*/
/*       43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55., 56.,*/
/*       57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69., 70.}};*/
/**/
/*  int bin{-999};*/
/*  for (int i = 0; i < 70; ++i) {*/
/*    if (fv0mpercentile >= low.at(i) && fv0mpercentile < high.at(i)) {*/
/*      bin = i;*/
/*      break;*/
/*    }*/
/*  }*/
/**/
/*  return bin;*/
/*}*/
