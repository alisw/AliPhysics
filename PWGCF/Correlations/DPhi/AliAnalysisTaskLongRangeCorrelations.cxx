// -*- C++ -*-
/**************************************************************************
 * Copyright(c) 2012,      ALICE Experiment at CERN, All rights reserved. *
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
 **************************************************************************/
// $Id: AliAnalysisTaskLongRangeCorrelations.cxx 217 2012-11-06 10:19:42Z cmayer $

#include <TChain.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TH2.h>

#include "AliEventPoolManager.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

#include "AliAnalysisTaskLongRangeCorrelations.h"


ClassImp(AliAnalysisTaskLongRangeCorrelations);
ClassImp(LRCParticle);

AliAnalysisTaskLongRangeCorrelations::AliAnalysisTaskLongRangeCorrelations(const char *name)
  : AliAnalysisTaskSE(name)
  , fOutputList(NULL)
  , fRunMixing(kFALSE)
  , fPoolMgr(NULL)
  , fTrackFilter(128)
  , fCentMin(0), fCentMax(20)
  , fPtMin(0.2), fPtMax(1e10)
  , fPhiMin(0.), fPhiMax(TMath::TwoPi())
  , fnBinsCent( 220), fnBinsPt(400), fnBinsPhi(4),              fnBinsEta(120)
  , fxMinCent( -5.0), fxMinPt( 0.0), fxMinPhi( 0.0),            fxMinEta(-1.5)
  , fxMaxCent(105.0), fxMaxPt( 4.0), fxMaxPhi( TMath::TwoPi()), fxMaxEta( 1.5) {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskLongRangeCorrelations::~AliAnalysisTaskLongRangeCorrelations() {
  if (NULL != fOutputList) {
    delete fOutputList;
    fOutputList = NULL;
  }
}

void AliAnalysisTaskLongRangeCorrelations::UserCreateOutputObjects() {
  fOutputList = new THashList;
  fOutputList->SetOwner(kTRUE);
  fOutputList->SetName(GetOutputListName());

  // Event statistics
  const char* eventStatLabels[] = { 
    "All Events",
    "Physics Selection",
    "Centrality Selection",
    "Analyzed Events"
  };
  const size_t nEventStat(sizeof(eventStatLabels)/sizeof(const char*));
  TH1* hStats(new TH1D("histEventStats", "histEventStats", nEventStat, -0.5, nEventStat-0.5));
  for (size_t i=0; i<nEventStat; ++i)
    hStats->GetXaxis()->SetBinLabel(1+i, eventStatLabels[i]);
  fOutputList->Add(hStats);

  // QA histograms
  fOutputList->Add(new TH2D("histQACentPt", "histQACentPt;centrality V0M(%);p_{T} (GeV/c);",
			    fnBinsCent, fxMinCent, fxMaxCent, fnBinsPt, fxMinPt, fxMaxPt));
  fOutputList->Add(new TH2D("histQAPhiEta", "histQAPhiEta;#phi (rad);#eta;",
			    200, 0.0, TMath::TwoPi(), 300, -1.5, 1.5));

  // Moments
  fOutputList->Add(MakeHistSparsePhiEta("histMoment1PhiEta_1"));
  fOutputList->Add(MakeHistSparsePhiEta("histMoment1PhiEta_2"));
  fOutputList->Add(MakeHistSparsePhiEtaPhiEta("histMoment2PhiEtaPhiEta_11"));
  fOutputList->Add(MakeHistSparsePhiEtaPhiEta("histMoment2PhiEtaPhiEta_12"));
  fOutputList->Add(MakeHistSparsePhiEtaPhiEta("histMoment2PhiEtaPhiEta_22"));

  if (fRunMixing)
    SetupForMixing();

  PostData(1, fOutputList);
}

void AliAnalysisTaskLongRangeCorrelations::UserExec(Option_t* ) {
  AliAnalysisManager* pManager(AliAnalysisManager::GetAnalysisManager());
  if (NULL == pManager) return;

  AliInputEventHandler* pInputHandler(dynamic_cast<AliInputEventHandler*>(pManager->GetInputEventHandler()));
  if (NULL == pInputHandler) return;
  
  AliAODEvent* pAOD(dynamic_cast<AliAODEvent*>(InputEvent()));
  if (NULL == pAOD) return;

  AliAODHeader *pAODHeader = pAOD->GetHeader();
  if (NULL == pAODHeader) return;

  Fill("histEventStats", 0.); // all events

  if (!pInputHandler->IsEventSelected()) return;
  Fill("histEventStats", 1.); // events passing physics selection

  const Double_t centrality(pAODHeader->GetCentralityP()->GetCentralityPercentile("V0M"));
  if (centrality < fCentMin || centrality >= fCentMax) return;
  Fill("histEventStats", 2.); // events passing centrality selection

  // event is accepted
  TObjArray* tracksMain(GetAcceptedTracks(pAOD, pAODHeader, centrality));

  if (fRunMixing) {  
    AliEventPool* pEventPool(fPoolMgr->GetEventPool(centrality, pAOD->GetPrimaryVertex()->GetZ()));
    if (NULL == pEventPool)
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality, pAOD->GetPrimaryVertex()->GetZ()));
    
//     pEventPool->PrintInfo();
    if (pEventPool->IsReady()
	|| pEventPool->NTracksInPool() > fMixingTracks/10
	|| pEventPool->GetCurrentNEvents() >= 5) {
      const Int_t nMix(pEventPool->GetCurrentNEvents());
      Fill("histEventStats", 3.); // analyzed events
      for (Int_t i(0); i<nMix; ++i) {
	TObjArray* tracksMixed(pEventPool->GetEvent(i));
	CalculateMoments(tracksMain, tracksMixed, 1./nMix);
      }
    }
    // Update the Event pool
    pEventPool->UpdatePool(tracksMain);
  } else { // no mixing
    Fill("histEventStats", 3.); // analyzed events
    CalculateMoments(tracksMain, tracksMain, 1.);
    delete tracksMain;
  }
}

void AliAnalysisTaskLongRangeCorrelations::Terminate(Option_t* ) {
  //
}

TString AliAnalysisTaskLongRangeCorrelations::GetOutputListName() const {
  TString listName("listLRC");
  listName += TString::Format("_%smix",         fRunMixing ? "" : "no");
  listName += TString::Format("_trackFilt%d",   fTrackFilter);
  listName += TString::Format("_cent%.0fT%.0f", fCentMin, fCentMax);
  listName += TString::Format("_ptMin%.0fMeV",  1e3*fPtMin);
  listName += TString::Format("_phi%.0fT%.0f",  TMath::RadToDeg()*fPhiMin, TMath::RadToDeg()*fPhiMax);
  return listName;
}

void AliAnalysisTaskLongRangeCorrelations::SetupForMixing() {
  const Int_t trackDepth(fMixingTracks);
  const Int_t poolsize(1000); // Maximum number of events

  Double_t centralityBins[] = { // centrality bins
    0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,90.,100.
//     0.,20.,100.
  };
  const Int_t nCentralityBins(sizeof(centralityBins)/sizeof(Double_t) - 1);

  Double_t vertexBins[] = {  // Zvtx bins
    -10., -7., -5., -3., -1., 1., 3., 5., 7., 10.
  };
  const Int_t nVertexBins(sizeof(vertexBins)/sizeof(Double_t) - 1);

  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth,
				     nCentralityBins, centralityBins,
				     nVertexBins,     vertexBins);
}

THnSparse* AliAnalysisTaskLongRangeCorrelations::MakeHistSparsePhiEta(const char* name) const {
  const Int_t   nBinsM1[] = { fnBinsPhi, fnBinsEta };
  const Double_t xMinM1[] = {  fxMinPhi,  fxMinEta };
  const Double_t xMaxM1[] = {  fxMaxPhi,  fxMaxEta };
  const TString title(TString(name)
		      +";#phi;#eta;");
  return new THnSparseD(name, title.Data(), 2, nBinsM1, xMinM1, xMaxM1);
}
THnSparse* AliAnalysisTaskLongRangeCorrelations::MakeHistSparsePhiEtaPhiEta(const char* name) const {
  const Int_t   nBinsM2[] = {  fnBinsPhi, fnBinsEta,  fnBinsPhi, fnBinsEta };
  const Double_t xMinM2[] = {  fxMinPhi,  fxMinEta,   fxMinPhi,  fxMinEta  };
  const Double_t xMaxM2[] = {  fxMaxPhi,  fxMaxEta,   fxMaxPhi,  fxMaxEta  };
  const TString title(TString(name)
		      +";#phi_{1};#eta_{1}"
		      +";#phi_{2};#eta_{2};");
  return new THnSparseD(name, title.Data(), 4, nBinsM2, xMinM2, xMaxM2);
}

TObjArray* AliAnalysisTaskLongRangeCorrelations::GetAcceptedTracks(AliAODEvent* pAOD,
								   AliAODHeader* /* pAODHeader */,
								   Double_t centrality) {
  TObjArray* tracks= new TObjArray;
  tracks->SetOwner(kTRUE);

  AliDebug(5, Form("#tracks= %6d %f", pAOD->GetNumberOfTracks(), centrality));
  for (Long64_t i(0); i<pAOD->GetNumberOfTracks(); ++i) {
    AliAODTrack* pAODTrack(dynamic_cast<AliAODTrack *>(pAOD->GetTrack(i)));
    if (NULL == pAODTrack) continue;    
    if (!pAODTrack->TestFilterBit(fTrackFilter)) continue;

    Fill("histQACentPt", centrality, pAODTrack->Pt());
    Fill("histQAPhiEta", pAODTrack->Phi(), pAODTrack->Eta());
    if (pAODTrack->Phi() < fPhiMin || pAODTrack->Phi() > fPhiMax) continue;
    if (pAODTrack->Pt()  < fPtMin  || pAODTrack->Pt()  > fPtMax)  continue;

    tracks->Add(new LRCParticle(pAODTrack->Eta(), pAODTrack->Phi()));
  } // next track
  return tracks;
}

void AliAnalysisTaskLongRangeCorrelations::CalculateMoments(TObjArray* tracks1,
							    TObjArray* tracks2,
							    Double_t weight) {
  // <n_1>
  THnSparse* hN1ForThisEvent(ComputeNForThisEvent(tracks1, "hN1"));
  THnSparse* hN1(dynamic_cast<THnSparse*>(fOutputList->FindObject("histMoment1PhiEta_1")));
  hN1->Add(hN1ForThisEvent, weight);

  // <n_2>
  THnSparse* hN2ForThisEvent(ComputeNForThisEvent(tracks2, "hN2"));
  THnSparse* hN2(dynamic_cast<THnSparse*>(fOutputList->FindObject("histMoment1PhiEta_2")));
  hN2->Add(hN2ForThisEvent, weight);

  // <n_1 n_1>
  ComputeN2ForThisEvent(hN1ForThisEvent, hN1ForThisEvent, "histMoment2PhiEtaPhiEta_11", weight);
  // <n_1 n_2>
  ComputeN2ForThisEvent(hN1ForThisEvent, hN2ForThisEvent, "histMoment2PhiEtaPhiEta_12", weight);
  // <n_2 n_2>
  ComputeN2ForThisEvent(hN2ForThisEvent, hN2ForThisEvent, "histMoment2PhiEtaPhiEta_22", weight);

  // clean up
  delete hN1ForThisEvent;
  delete hN2ForThisEvent;
}

THnSparse* AliAnalysisTaskLongRangeCorrelations::ComputeNForThisEvent(TObjArray* tracks, const char* histName) const {
  THnSparse* hN(MakeHistSparsePhiEta(histName));
  const Long64_t nTracks(tracks->GetEntries());
  for (Long64_t i(0); i<nTracks; ++i) {
    const LRCParticle* p(dynamic_cast<LRCParticle*>(tracks->At(i)));
    const Double_t x[] = { p->Phi(), p->Eta() };
    hN->Fill(x);
  }
  return hN;
}

void AliAnalysisTaskLongRangeCorrelations::ComputeN2ForThisEvent(THnSparse* hN1, THnSparse* hN2,
								 const char* histName, Double_t weight) {
  if (NULL == fOutputList) return;
  THnSparse* hs(dynamic_cast<THnSparse*>(fOutputList->FindObject(histName)));
  if (hs == NULL) return;

  for (Long64_t i(0); i<hN1->GetNbins(); ++i) {
    Double_t x[] = {0,0, 0,0};
    Int_t coord1[] = {0,0};
    const Double_t n1(hN1->GetBinContent(i, coord1));
    for (Int_t k(0); k<2; ++k)
      x[k] = hN1->GetAxis(k)->GetBinCenter(coord1[k]);

    for (Long64_t j(0); j<hN2->GetNbins(); ++j) {
      Int_t coord2[] = {0,0};
      const Double_t n2(hN2->GetBinContent(j, coord2));
      for (Int_t k(0); k<2; ++k)
	x[2+k] = hN2->GetAxis(k)->GetBinCenter(coord2[k]);

      hs->Fill(x, weight*n1*n2);
    }
  }
}

void AliAnalysisTaskLongRangeCorrelations::Fill(const char* histName, Double_t x) {
  if (NULL == fOutputList) return;
  TH1* h = dynamic_cast<TH1*>(fOutputList->FindObject(histName));
  if (h == NULL) return;
  h->Fill(x);
}
void AliAnalysisTaskLongRangeCorrelations::Fill(const char* histName, Double_t x, Double_t y) {
  if (NULL == fOutputList) return;
  TH2* h = dynamic_cast<TH2*>(fOutputList->FindObject(histName));
  if (h == NULL) return;
  h->Fill(x, y);
}
void AliAnalysisTaskLongRangeCorrelations::Fill(const char* histName, const Double_t* x, Double_t w) {
  if (NULL == fOutputList) return;
  THnSparse* h = dynamic_cast<THnSparse*>(fOutputList->FindObject(histName));
  if (h == NULL) return;
  h->Fill(x, w);
}
