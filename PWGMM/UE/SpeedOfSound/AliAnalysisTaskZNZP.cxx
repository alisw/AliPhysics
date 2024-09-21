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
      fIsSystematics(true),
      fVaryVtxZPos(false),
      fSaveAsy(false),
      fMinVtxZPos(-5.0),
      fMaxVtxZPos(5.0),
      fSystematic(1),
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
      pV0MAmpChannel(0),
      hV0Percentile(0),
      hBestVtxZ(0),
      hZNvsV0MPer(0),
      hZNvsV0M(0),
      hZNAvsV0M(0),
      hZNCvsV0M(0),
      hAsyN(0),
      hZPvsV0MPer(0),
      hZPvsV0M(0),
      hZPAvsV0M(0),
      hZPCvsV0M(0),
      /*hAsyP(0),*/
      hZNCpmc(0),
      hZNApmc(0),
      hZPCpmc(0),
      hZPApmc(0),
      hZNCNorm(0),
      hZNANorm(0) {}
//_____________________________________________________________________________
AliAnalysisTaskZNZP::AliAnalysisTaskZNZP(const char* name)
    : AliAnalysisTaskSE(name),
      fESD(0),
      fEventCuts(0x0),
      fMCStack(0),
      fMC(0),
      fUseMC(kFALSE),
      fIsSystematics(true),
      fVaryVtxZPos(false),
      fSaveAsy(false),
      fMinVtxZPos(-5.0),
      fMaxVtxZPos(5.0),
      fSystematic(1),
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
      pV0MAmpChannel(0),
      hV0Percentile(0),
      hBestVtxZ(0),
      hZNvsV0MPer(0),
      hZNvsV0M(0),
      hZNAvsV0M(0),
      hZNCvsV0M(0),
      hAsyN(0),
      hZPvsV0MPer(0),
      hZPvsV0M(0),
      hZPAvsV0M(0),
      hZPCvsV0M(0),
      /*hAsyP(0),*/
      hZNCpmc(0),
      hZNApmc(0),
      hZPCpmc(0),
      hZPApmc(0),
      hZNCNorm(0),
      hZNANorm(0) {
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

    /*if (fIsSystematics) {*/
    /*  ChangeCut(fCuts);*/
    /*}*/
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

    /*if (fIsSystematics) {*/
    /*  ChangeCut(fCuts3);*/
    /*}*/
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

  /*constexpr double v0mAmp_width{25.0};*/
  /*constexpr int v0mAmp_Nbins{1720};*/
  /*double v0mAmp_bins[v0mAmp_Nbins + 1] = {0};*/
  /*for (int i = 0; i <= v0mAmp_Nbins; ++i) {*/
  /*  v0mAmp_bins[i] = 0.0 + i * v0mAmp_width;*/
  /*}*/

  const int nBinsV0M090{70};
  double BinsV0M090[nBinsV0M090 + 1] = {0.0};
  for (int i = 0; i <= nBinsV0M090; ++i) {
    BinsV0M090[i] = 0.0 + (double)i;
  }

  /*constexpr int v0m_Nbins080{6};*/
  /*constexpr double v0m_bins080[v0m_Nbins080 + 1] = {0.0,  1.0,  5.0, 10.0,*/
  /*                                                  20.0, 50.0, 80.0};*/

  hV0Percentile = new TH1F("hV0M", ";V0M (%);Entries", nBinsV0M090, BinsV0M090);
  pV0MAmpChannel =
      new TProfile("pV0Channel", ";Channel; Amplitude;", 64, -0.5, 63.5);
  /*hV0MAmplitude =*/
  /*    new TH1F("hV0MAmp", ";V0M amplitude;Counts", v0mAmp_Nbins,
   * v0mAmp_bins);*/
  hBestVtxZ =
      new TH1F("hBestVtxZ", ";Vertex_{#it{z}} (cm); Counts;", 400, -11, 11);

  hZNvsV0M = new TH2F("hZNvsV0M", "ZNA+ZNC;V0M Per; #it{E}_{ZN} [TeV]/2.511;",
                      nBinsV0M090, BinsV0M090, 200, 0.0, 200.0);
  hZPvsV0M = new TH2F("hZPvsV0M", "ZPA+ZPC;V0M Per; #it{E}_{ZP} [TeV]/2.511;",
                      nBinsV0M090, BinsV0M090, 60, 0.0, 60.0);

  hZNAvsV0M = new TH2F("hZNAvsV0M", "ZNA;V0M Per; #it{E}_{ZN} [TeV]/2.511;",
                       nBinsV0M090, BinsV0M090, 100, 0.0, 100.0);
  hZNCvsV0M = new TH2F("hZNCvsV0M", "ZNC;V0M Per; #it{E}_{ZN} [TeV]/2.511;",
                       nBinsV0M090, BinsV0M090, 100, 0.0, 100.0);

  hZPAvsV0M = new TH2F("hZPAvsV0M", "ZPA;V0M Per; #it{E}_{ZP} [TeV]/2.511;",
                       nBinsV0M090, BinsV0M090, 30, 0.0, 30.0);
  hZPCvsV0M = new TH2F("hZPCvsV0M", "ZPC;V0M Per; #it{E}_{ZP} [TeV]/2.511;",
                       nBinsV0M090, BinsV0M090, 30, 0.0, 30.0);
  hAsyN =
      new TH2F("hAsyN", "Neutron asymmetry;V0M Per; N_{C}-N_{A}/N_{C}+N_{A};",
               nBinsV0M090, BinsV0M090, 50, -1.0, 1.0);
  /*hAsyP =*/
  /*    new TH2F("hAsyP", "Proton asymmetry;V0M Per;
   * P_{C}-P_{A}/P_{C}+P_{A};",*/
  /*             nBinsV0M090, BinsV0M090, 50, -1.0, 1.0);*/
  hZNvsV0MPer =
      new TH2F("hZNvsV0MPer", "(ZNA+ZNC)/2;V0M Per; #it{E}_{ZN} [TeV]/2.511;",
               nBinsV0M090, BinsV0M090, 100, 0.0, 100.0);
  hZPvsV0MPer =
      new TH2F("hZPvsV0MPer", "(ZPA+ZPC)/2;V0M Per; #it{E}_{ZP} [TeV]/2.511;",
               nBinsV0M090, BinsV0M090, 30, 0.0, 30.0);

  hZNCNorm = new TH2F("hZNCNorm", ";V0M;<ZNC>/<ZNC>;", nBinsV0M090, BinsV0M090,
                      40, 0., 2.);

  hZNANorm = new TH2F("hZNANorm", ";V0M;<ZNA>/<ZNA>;", nBinsV0M090, BinsV0M090,
                      40, 0., 2.);

  /*hZNCNormSca = new TH2F("hZNCNormSca", ";V0M;<ZNC>/<ZNC> with #sigma
   * scaled;",*/
  /*                       nBinsV0M090, BinsV0M090, 40, 0., 2.);*/

  /*hZNANormSca = new TH2F("hZNANormSca", ";V0M;<ZNA>/<ZNA> with #sigma scaled;",*/
  /*                       nBinsV0M090, BinsV0M090, 40, 0., 2.);*/

  hZNCpmc = new TH1F("hZNCpmc", "ZNC PMC;ZNC energy;Entries", 520, 0., 130.);
  hZNApmc = new TH1F("hZNApmc", "ZNA PMC;ZNA energy;Entries", 520, 0., 130.);
  hZPCpmc = new TH1F("hZPCpmc", "ZPC PMC;ZPC energy;Entries", 120, 0., 30.);
  hZPApmc = new TH1F("hZPApmc", "ZPA PMC;ZPA energy;Entries", 120, 0., 30.);

  fOutputList->Add(hBestVtxZ);
  fOutputList->Add(hV0Percentile);
  fOutputList->Add(hZNvsV0MPer);
  fOutputList->Add(hZNvsV0M);
  fOutputList->Add(hZNAvsV0M);
  fOutputList->Add(hZNCvsV0M);
  fOutputList->Add(hZNCpmc);
  fOutputList->Add(hZNApmc);
  fOutputList->Add(hZNCNorm);
  fOutputList->Add(hZNANorm);
  if (fSaveAsy) {
    fOutputList->Add(hAsyN);
    /*fOutputList->Add(hZNCNormSca);*/
    /*fOutputList->Add(hZNANormSca);*/
  }

  fOutputList->Add(hZPvsV0MPer);
  fOutputList->Add(hZPvsV0M);
  fOutputList->Add(hZPAvsV0M);
  fOutputList->Add(hZPCvsV0M);
  /*fOutputList->Add(hAsyP);*/
  fOutputList->Add(hZPCpmc);
  fOutputList->Add(hZPApmc);

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
  /*GetCalibratedV0Amplitude();*/

  //! Get ZDC Centrality
  hV0Percentile->Fill(fv0mpercentile);
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

  const double znc{esdZDC->GetZDCN1Energy() / 1000.0};
  const double zna{esdZDC->GetZDCN2Energy() / 1000.0};
  const double zpc{esdZDC->GetZDCP1Energy() / 1000.0};
  const double zpa{esdZDC->GetZDCP2Energy() / 1000.0};

  double meanZNC{1.0};
  double sigmaZNC{1.0};
  double meanZNA{1.0};
  double sigmaZNA{1.0};
  bool gCent1{MeanSigmaZN(meanZNC, sigmaZNC, "ZNC")};
  bool gCent2{MeanSigmaZN(meanZNA, sigmaZNA, "ZNA")};
  if (!gCent1 && !gCent2) {
    return;
  }

  hZNCNorm->Fill(fv0mpercentile, znc / (meanZNC * 2.511));
  hZNANorm->Fill(fv0mpercentile, zna / (meanZNA * 2.511));
  const double znc_sca{znc / (meanZNC * 2.511)};
  const double zna_sca{zna / (meanZNA * 2.511)};

  // Non-average ZN & ZP
  hZNvsV0M->Fill(fv0mpercentile, (znc + zna) / 2.511);
  hZPvsV0M->Fill(fv0mpercentile, (zpc + zpa) / 2.511);

  hZNCvsV0M->Fill(fv0mpercentile, znc / 2.511);
  hZNAvsV0M->Fill(fv0mpercentile, zna / 2.511);

  hZPCvsV0M->Fill(fv0mpercentile, zpc / 2.511);
  hZPAvsV0M->Fill(fv0mpercentile, zpa / 2.511);

  hAsyN->Fill(fv0mpercentile, (znc_sca - zna_sca) / (znc_sca + zna_sca));
  /*hAsyP->Fill(fv0mpercentile, (zpc - zpa) / (zpc + zpa));*/

  const double* towZNC{esdZDC->GetZN1TowerEnergy()};
  const double* towZPC{esdZDC->GetZP1TowerEnergy()};
  const double* towZNA{esdZDC->GetZN2TowerEnergy()};
  const double* towZPA{esdZDC->GetZP2TowerEnergy()};
  hZNApmc->Fill(towZNA[0] / 1000.0);
  hZNCpmc->Fill(towZNC[0] / 1000.0);
  hZPApmc->Fill(towZPA[0] / 1000.0);
  hZPCpmc->Fill(towZPC[0] / 1000.0);

  /*const double* towZNCLG{esdZDC->GetZN1TowerEnergyLR()};*/
  /*const double* towZNALG{esdZDC->GetZN2TowerEnergyLR()};*/
  /*cout << "ZNC LG= " << towZNCLG[0] << " | ZNA LG= " << towZNALG[0] <<
   * '\n';*/

  // Average the energy detected in each calorimeter
  hZNvsV0MPer->Fill(fv0mpercentile, (znc + zna) / (2.0 * 2.511));
  /*hZNvsV0MAmp->Fill(fv0mamplitude, (znc + zna) / (2.0 * 2.511));*/
  /*pZNvsV0MAmp->Fill(fv0mamplitude, (znc + zna) / 2.0);*/

  hZPvsV0MPer->Fill(fv0mpercentile, (zpc + zpa) / (2.0 * 2.511));
  /*hZPvsV0MAmp->Fill(fv0mamplitude, (zpc + zpa) / (2.0 * 2.511));*/
  /*pZPvsV0MAmp->Fill(fv0mamplitude, (zpc + zpa) / 2.0);*/
}
//____________________________________________________________
void AliAnalysisTaskZNZP::GetCalibratedV0Amplitude() {
  float mV0M{0.0};
  for (int i = 0; i < 64; i++) {
    mV0M += fESD->GetVZEROEqMultiplicity(i);
    pV0MAmpChannel->Fill(i, fESD->GetVZEROEqMultiplicity(i));
  }
  fv0mamplitude = mV0M;

  /*hV0MAmplitude->Fill(fv0mamplitude);*/
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
bool AliAnalysisTaskZNZP::MeanSigmaZN(double& mean, double& sigma,
                                      const std::string& ZN) {
  const std::vector<double> meanZNC{
      {9.41689, 11.4554, 13.6935, 15.899,  17.9693, 19.8913, 21.8775, 23.6385,
       25.2471, 26.8672, 28.3384, 29.8036, 31.0197, 32.3755, 33.5516, 34.5927,
       35.6032, 36.6078, 37.5322, 38.3336, 39.1421, 39.8723, 40.5853, 41.1441,
       41.7837, 42.2702, 42.8783, 43.321,  43.7049, 44.1573, 44.5504, 44.8895,
       45.1275, 45.3678, 45.7246, 45.7786, 45.9834, 46.1794, 46.2849, 46.2833,
       46.3337, 46.4907, 46.4214, 46.344,  46.3888, 46.2107, 46.1971, 46.0917,
       45.893,  45.6904, 45.4781, 45.4185, 45.0878, 44.8286, 44.5048, 44.2761,
       43.9231, 43.421,  43.1319, 42.6159, 42.1982, 41.7332, 41.1261, 40.8089,
       40.1382, 39.5972, 39.127,  38.4715, 37.7199, 36.9609}};

  const std::vector<double> sigmaZNC{
      {3.50677, 3.82505, 4.11002, 4.24869, 4.46803, 4.59024, 4.90368, 5.07824,
       5.17144, 5.19741, 5.30018, 5.51713, 5.46297, 5.66198, 5.76174, 5.80704,
       5.83463, 5.97203, 6.03367, 6.11339, 6.12927, 6.17397, 6.24563, 6.32831,
       6.45077, 6.45699, 6.46789, 6.52385, 6.57878, 6.67186, 6.77811, 6.72446,
       6.85183, 6.85365, 6.96893, 6.95236, 6.9604,  7.09092, 7.10086, 7.13674,
       7.2995,  7.24527, 7.29297, 7.38196, 7.45292, 7.48775, 7.59025, 7.67292,
       7.73854, 7.77752, 7.84931, 8.00492, 7.91375, 8.13183, 8.19829, 8.33288,
       8.42185, 8.57176, 8.7218,  8.87617, 9.03015, 9.14355, 9.28268, 9.51969,
       9.68484, 9.92729, 10.0658, 10.3063, 10.4252, 10.6264}};

  const std::vector<double> meanZNA{
      {7.64853, 9.31038, 11.1556, 13.0141, 14.6987, 16.4405, 17.8888, 19.4125,
       20.8443, 22.1771, 23.4108, 24.5767, 25.6076, 26.716,  27.7061, 28.6364,
       29.557,  30.2871, 31.0859, 31.7577, 32.4302, 33.0959, 33.7153, 34.2552,
       34.8104, 35.1898, 35.7727, 36.1279, 36.4649, 36.8997, 37.2122, 37.3977,
       37.6971, 37.9211, 38.1646, 38.2855, 38.5034, 38.6686, 38.7706, 38.8934,
       38.9517, 38.975,  39.0276, 39.0661, 38.9757, 39.0779, 38.8923, 38.8756,
       38.7662, 38.6167, 38.527,  38.2397, 38.1209, 37.9127, 37.763,  37.4913,
       37.2047, 36.9349, 36.5257, 36.207,  35.8965, 35.5434, 35.0563, 34.6988,
       34.1538, 33.8234, 33.3074, 32.82,   32.1863, 31.6347}};

  const std::vector<double> sigmaZNA{
      {3.10314, 3.13242, 3.48659, 3.77875, 3.7878,  3.915,   4.03227, 4.20045,
       4.31136, 4.40799, 4.56356, 4.61643, 4.65824, 4.74357, 4.78761, 4.92246,
       4.99679, 5.05502, 5.0036,  5.12261, 5.17488, 5.19715, 5.26379, 5.20759,
       5.38565, 5.41547, 5.42096, 5.47782, 5.47648, 5.65329, 5.61113, 5.61817,
       5.70361, 5.62584, 5.73757, 5.73626, 5.75184, 5.85641, 5.92166, 5.9685,
       5.99215, 6.08534, 6.02079, 6.1867,  6.17472, 6.26609, 6.2535,  6.34106,
       6.32414, 6.48443, 6.57763, 6.57414, 6.70035, 6.78231, 6.80588, 6.87396,
       7.06222, 7.13683, 7.23114, 7.44464, 7.49467, 7.61375, 7.7578,  7.96836,
       8.12227, 8.33368, 8.3856,  8.5707,  8.92543, 8.95762}};

  bool goodCent{false};
  if (fv0mpercentile >= 0. && fv0mpercentile < 70.) goodCent = true;

  int bin{CentBin()};
  if (!goodCent || bin < 0) return false;

  if (ZN == "ZNC") {
    mean = meanZNC.at(bin);
    sigma = sigmaZNC.at(bin);
  }
  if (ZN == "ZNA") {
    mean = meanZNA.at(bin);
    sigma = sigmaZNA.at(bin);
  }

  return goodCent;
}
//____________________________________________________________
int AliAnalysisTaskZNZP::CentBin() {
  const std::vector<double> low{
      {0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12., 13.,
       14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27.,
       28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41.,
       42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55.,
       56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69.}};
  const std::vector<double> high{
      {1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 11., 12., 13., 14.,
       15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28.,
       29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42.,
       43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55., 56.,
       57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69., 70.}};

  int bin{-999};
  for (int i = 0; i < 70; ++i) {
    if (fv0mpercentile >= low.at(i) && fv0mpercentile < high.at(i)) {
      bin = i;
      break;
    }
  }

  return bin;
}
