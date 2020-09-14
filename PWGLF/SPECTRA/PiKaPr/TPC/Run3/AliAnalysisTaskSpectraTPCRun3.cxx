/**************************************************************************
*  Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                         *
*  Author: The ALICE Off-line Project.                                    *
*  Contributors are mentioned in the code where appropriate.              *
*                                                                         *
*  Permission to use, copy, modify and distribute this software and its   *
*  documentation strictly for non-commercial purposes is hereby granted   *
*  without fee, provided that the above copyright notice appears in all   *
*  copies and that both the copyright notice and this permission notice   *
*  appear in the supporting documentation. The authors make no claims     *
*  about the suitability of this software for any purpose. It is          *
*  provided "as is" without express or implied warranty.                  *
**************************************************************************/

///
/// \author Nicolo' Jacazio
/// \file   AliAnalysisTaskSpectraTPCRun3.cxx
/// \brief  Task to produce a TPC spectra with O2 standards
/// \since  14/09/2020
///

#include "AliAnalysisTaskSpectraTPCRun3.h"

ClassImp(AliAnalysisTaskSpectraTPCRun3);

AliAnalysisTaskSpectraTPCRun3::AliAnalysisTaskSpectraTPCRun3()
    : AliAnalysisTaskSE("")
    , fEventCut()
{
  // default constructor, nothing to do
}

AliAnalysisTaskSpectraTPCRun3::AliAnalysisTaskSpectraTPCRun3(const char* name)
    : AliAnalysisTaskSE(name)
    , fEventCut()
{
  bbparam[0] = 0.0330656;
  bbparam[1] = 19.9759;
  bbparam[2] = -7.31532e-09;
  bbparam[3] = 2.72023;
  bbparam[4] = 6.0812;
  bbparam[5] = 51.5071;
  bbparam[6] = 2.3;
  bbresoparam[0] = 0.07;
  bbresoparam[1] = 0.;
  //
  // main constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskSpectraTPCRun3::~AliAnalysisTaskSpectraTPCRun3()
{
  //
  // destructor
  //
  if (fOutputList) {
    delete fOutputList;
  }
}

void AliAnalysisTaskSpectraTPCRun3::UserCreateOutputObjects()
{
  //
  // create the output objects
  //
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
// Histograms
#define DOTH1F(OBJ, ...)               \
  {                                    \
    OBJ = new TH1F(#OBJ, __VA_ARGS__); \
    fOutputList->Add(OBJ);             \
  }
#define DOTH2F(OBJ, ...)               \
  {                                    \
    OBJ = new TH2F(#OBJ, __VA_ARGS__); \
    fOutputList->Add(OBJ);             \
  }

#define makelogaxis(h)                                            \
  {                                                               \
    const Int_t nbins = h->GetNbinsX();                           \
    double binp[nbins + 1];                                       \
    double max = h->GetXaxis()->GetBinUpEdge(nbins);              \
    double min = h->GetXaxis()->GetBinLowEdge(1);                 \
    double lmin = TMath::Log10(min);                              \
    double ldelta = (TMath::Log10(max) - lmin) / ((double)nbins); \
    for (int i = 0; i < nbins; i++) {                             \
      binp[i] = TMath::Exp(TMath::Log(10) * (lmin + i * ldelta)); \
    }                                                             \
    binp[nbins] = max + 1;                                        \
    h->GetXaxis()->Set(nbins, binp);                              \
  }

#define BIN_AXIS 1000, 0.001, 20, 1000, 0, 1000

  DOTH2F(htpcsignal, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);
  DOTH2F(hexpEl, ";#it{p} (GeV/#it{c});TPC expected signal e;Tracks", BIN_AXIS);
  DOTH2F(hexpMu, ";#it{p} (GeV/#it{c});TPC expected signal #mu;Tracks", BIN_AXIS);
  DOTH2F(hexpPi, ";#it{p} (GeV/#it{c});TPC expected signal #pi;Tracks", BIN_AXIS);
  DOTH2F(hexpKa, ";#it{p} (GeV/#it{c});TPC expected signal K;Tracks", BIN_AXIS);
  DOTH2F(hexpPr, ";#it{p} (GeV/#it{c});TPC expected signal p;Tracks", BIN_AXIS);
  DOTH2F(hexpDe, ";#it{p} (GeV/#it{c});TPC expected signal d;Tracks", BIN_AXIS);
  DOTH2F(hexpTr, ";#it{p} (GeV/#it{c});TPC expected signal t;Tracks", BIN_AXIS);
  DOTH2F(hexpHe, ";#it{p} (GeV/#it{c});TPC expected signal ^{3}He;Tracks", BIN_AXIS);
  DOTH2F(hexpAl, ";#it{p} (GeV/#it{c});TPC expected signal #alpha;Tracks", BIN_AXIS);

  makelogaxis(htpcsignal);
  makelogaxis(hexpEl);
  makelogaxis(hexpMu);
  makelogaxis(hexpPi);
  makelogaxis(hexpKa);
  makelogaxis(hexpPr);
  makelogaxis(hexpDe);
  makelogaxis(hexpTr);
  makelogaxis(hexpHe);
  makelogaxis(hexpAl);
#undef BIN_AXIS

#define BIN_AXIS 1000, 0.001, 20, 1000, -10, 10
  DOTH2F(hnsigmaEl, ";#it{p} (GeV/#it{c});TPC N_{sigma e};Tracks", BIN_AXIS);
  DOTH2F(hnsigmaMu, ";#it{p} (GeV/#it{c});TPC N_{sigma #mu};Tracks", BIN_AXIS);
  DOTH2F(hnsigmaPi, ";#it{p} (GeV/#it{c});TPC N_{sigma #pi};Tracks", BIN_AXIS);
  DOTH2F(hnsigmaKa, ";#it{p} (GeV/#it{c});TPC N_{sigma K};Tracks", BIN_AXIS);
  DOTH2F(hnsigmaPr, ";#it{p} (GeV/#it{c});TPC N_{sigma p};Tracks", BIN_AXIS);
  DOTH2F(hnsigmaDe, ";#it{p} (GeV/#it{c});TPC N_{sigma d};Tracks", BIN_AXIS);
  DOTH2F(hnsigmaTr, ";#it{p} (GeV/#it{c});TPC N_{sigma t};Tracks", BIN_AXIS);
  DOTH2F(hnsigmaHe, ";#it{p} (GeV/#it{c});TPC N_{sigma ^{3}He};Tracks", BIN_AXIS);
  DOTH2F(hnsigmaAl, ";#it{p} (GeV/#it{c});TPC N_{sigma #alpha};Tracks", BIN_AXIS);

  makelogaxis(hnsigmaEl);
  makelogaxis(hnsigmaMu);
  makelogaxis(hnsigmaPi);
  makelogaxis(hnsigmaKa);
  makelogaxis(hnsigmaPr);
  makelogaxis(hnsigmaDe);
  makelogaxis(hnsigmaTr);
  makelogaxis(hnsigmaHe);
  makelogaxis(hnsigmaAl);
#undef BIN_AXIS

  if (fMakeTOFPlots) {
#define BIN_AXIS 1000, 0.001, 20, 1000, 0, 1000

    DOTH2F(htpcsignalEl, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);
    DOTH2F(htpcsignalMu, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);
    DOTH2F(htpcsignalPi, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);
    DOTH2F(htpcsignalKa, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);
    DOTH2F(htpcsignalPr, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);
    DOTH2F(htpcsignalDe, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);
    DOTH2F(htpcsignalTr, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);
    DOTH2F(htpcsignalHe, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);
    DOTH2F(htpcsignalAl, ";#it{p} (GeV/#it{c});TPC Signal;Tracks", BIN_AXIS);

    makelogaxis(htpcsignalEl);
    makelogaxis(htpcsignalMu);
    makelogaxis(htpcsignalPi);
    makelogaxis(htpcsignalKa);
    makelogaxis(htpcsignalPr);
    makelogaxis(htpcsignalDe);
    makelogaxis(htpcsignalTr);
    makelogaxis(htpcsignalHe);
    makelogaxis(htpcsignalAl);
#undef BIN_AXIS
  }

  // Pt
#define TIT ";#it{p}_{T} (GeV/#it{c});Tracks"
  DOTH1F(hpt_El, TIT, 100, 0, 20);
  DOTH1F(hpt_Pi, TIT, 100, 0, 20);
  DOTH1F(hpt_Ka, TIT, 100, 0, 20);
  DOTH1F(hpt_Pr, TIT, 100, 0, 20);
#undef TIT
  // P
#define TIT ";#it{p} (GeV/#it{c});Tracks"
  DOTH1F(hp_El, TIT, 100, 0, 20);
  DOTH1F(hp_Pi, TIT, 100, 0, 20);
  DOTH1F(hp_Ka, TIT, 100, 0, 20);
  DOTH1F(hp_Pr, TIT, 100, 0, 20);
#undef TIT

#undef makelogaxis
  TList* lev = new TList();
  lev->SetOwner(kTRUE);
  fEventCut.AddQAplotsToList(lev);
  fOutputList->Add(lev);
  //
  PostData(1, fOutputList);
}

void AliAnalysisTaskSpectraTPCRun3::UserExec(Option_t*)
{
  //
  // loop over events
  //
  // AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  // if (man) {
  //   AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
  //   if (inputHandler)
  //     fPIDResponse = inputHandler->GetPIDResponse();
  // }

  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent)
    return;
  if (!fEventCut.AcceptEvent(fEvent)) {
    PostData(1, fOutputList);
    return;
  }

  for (Int_t itrk = 0; itrk < fEvent->GetNumberOfTracks(); itrk++) {

    /* get track */
    fTrack = static_cast<AliAODTrack*>(fEvent->GetTrack(itrk));
    if (!fTrack)
      continue;
    /* check accept track */
    if (!fTrack->TestFilterBit(32))
      continue;

    const float mom = fTrack->GetTPCmomentum();
    htpcsignal->Fill(mom, fTrack->GetTPCsignal());
    tpcExpSignalEl = ExpectedSignal(AliPID::kElectron);
    tpcExpSignalMu = ExpectedSignal(AliPID::kMuon);
    tpcExpSignalPi = ExpectedSignal(AliPID::kPion);
    tpcExpSignalKa = ExpectedSignal(AliPID::kKaon);
    tpcExpSignalPr = ExpectedSignal(AliPID::kProton);
    tpcExpSignalDe = ExpectedSignal(AliPID::kDeuteron);
    tpcExpSignalTr = ExpectedSignal(AliPID::kTriton);
    tpcExpSignalHe = ExpectedSignal(AliPID::kHe3);
    tpcExpSignalAl = ExpectedSignal(AliPID::kAlpha);

    hexpEl->Fill(mom, tpcExpSignalEl);
    hexpMu->Fill(mom, tpcExpSignalMu);
    hexpPi->Fill(mom, tpcExpSignalPi);
    hexpKa->Fill(mom, tpcExpSignalKa);
    hexpPr->Fill(mom, tpcExpSignalPr);
    hexpDe->Fill(mom, tpcExpSignalDe);
    hexpTr->Fill(mom, tpcExpSignalTr);
    hexpHe->Fill(mom, tpcExpSignalHe);
    hexpAl->Fill(mom, tpcExpSignalAl);

    if (ExpectedReso() > 0) {
      tpcNSigmaEl = (fTrack->GetTPCsignal() - tpcExpSignalEl) / ExpectedReso();
      tpcNSigmaMu = (fTrack->GetTPCsignal() - tpcExpSignalMu) / ExpectedReso();
      tpcNSigmaPi = (fTrack->GetTPCsignal() - tpcExpSignalPi) / ExpectedReso();
      tpcNSigmaKa = (fTrack->GetTPCsignal() - tpcExpSignalKa) / ExpectedReso();
      tpcNSigmaPr = (fTrack->GetTPCsignal() - tpcExpSignalPr) / ExpectedReso();
      tpcNSigmaDe = (fTrack->GetTPCsignal() - tpcExpSignalDe) / ExpectedReso();
      tpcNSigmaTr = (fTrack->GetTPCsignal() - tpcExpSignalTr) / ExpectedReso();
      tpcNSigmaHe = (fTrack->GetTPCsignal() - tpcExpSignalHe) / ExpectedReso();
      tpcNSigmaAl = (fTrack->GetTPCsignal() - tpcExpSignalAl) / ExpectedReso();
    } else {
      tpcNSigmaEl = 0;
      tpcNSigmaMu = 0;
      tpcNSigmaPi = 0;
      tpcNSigmaKa = 0;
      tpcNSigmaPr = 0;
      tpcNSigmaDe = 0;
      tpcNSigmaTr = 0;
      tpcNSigmaHe = 0;
      tpcNSigmaAl = 0;
    }

    hnsigmaEl->Fill(fTrack->P(), tpcNSigmaEl);
    hnsigmaMu->Fill(fTrack->P(), tpcNSigmaMu);
    hnsigmaPi->Fill(fTrack->P(), tpcNSigmaPi);
    hnsigmaKa->Fill(fTrack->P(), tpcNSigmaKa);
    hnsigmaPr->Fill(fTrack->P(), tpcNSigmaPr);
    hnsigmaDe->Fill(fTrack->P(), tpcNSigmaDe);
    hnsigmaTr->Fill(fTrack->P(), tpcNSigmaTr);
    hnsigmaHe->Fill(fTrack->P(), tpcNSigmaHe);
    hnsigmaAl->Fill(fTrack->P(), tpcNSigmaAl);

    if (fMakeTOFPlots) { // Check of PID with TOF information
      htpcsignalEl->Fill(mom, fTrack->GetTPCsignal());
      htpcsignalMu->Fill(mom, fTrack->GetTPCsignal());
      htpcsignalPi->Fill(mom, fTrack->GetTPCsignal());
      htpcsignalKa->Fill(mom, fTrack->GetTPCsignal());
      htpcsignalPr->Fill(mom, fTrack->GetTPCsignal());
      htpcsignalDe->Fill(mom, fTrack->GetTPCsignal());
      htpcsignalTr->Fill(mom, fTrack->GetTPCsignal());
      htpcsignalHe->Fill(mom, fTrack->GetTPCsignal());
      htpcsignalAl->Fill(mom, fTrack->GetTPCsignal());
    }

    if (TMath::Abs(tpcNSigmaEl) < 3) {
      hp_El->Fill(fTrack->P());
      hpt_El->Fill(fTrack->Pt());
    }
    if (TMath::Abs(tpcNSigmaPi) < 3) {
      hp_Pi->Fill(fTrack->P());
      hpt_Pi->Fill(fTrack->Pt());
    }
    if (TMath::Abs(tpcNSigmaKa) < 3) {
      hp_Ka->Fill(fTrack->P());
      hpt_Ka->Fill(fTrack->Pt());
    }
    if (TMath::Abs(tpcNSigmaPr) < 3) {
      hp_Pr->Fill(fTrack->P());
      hpt_Pr->Fill(fTrack->Pt());
    }
  }
  //
  // post the data and end the event loop
  //
  PostData(1, fOutputList);
}

void AliAnalysisTaskSpectraTPCRun3::Terminate(Option_t*) {}
