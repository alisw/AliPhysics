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
// $Id: AliAnalysisTaskLongRangeCorrelations.cxx 241 2012-12-19 19:45:59Z cmayer $

#include <numeric>
#include <functional>

#include <TChain.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include "AliEventPoolManager.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

#include "AliAnalysisTaskLongRangeCorrelations.h"


ClassImp(AliAnalysisTaskLongRangeCorrelations);
ClassImp(LRCParticle);

AliAnalysisTaskLongRangeCorrelations::AliAnalysisTaskLongRangeCorrelations(const char *name)
  : AliAnalysisTaskSE(name)
  , fOutputList(NULL)
  , fVertexZ(NULL)
  , fRunMixing(kFALSE)
  , fPoolMgr(NULL)
  , fMixingTracks(50000)
  , fTrackFilter(128)
  , fCentMin(0), fCentMax(20)
  , fPtMin(0.2), fPtMax(1e10)
  , fPhiMin(0.), fPhiMax(TMath::TwoPi())
  , fMaxAbsVertexZ(10.)
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

  const Double_t vertexBins[] = {  // Zvtx bins
    -10., -7., -5., -3., -1., 1., 3., 5., 7., 10.
  };
  fVertexZ = new TAxis(sizeof(vertexBins)/sizeof(Double_t)-1, vertexBins);
  fVertexZ->SetName("VertexZAxis");
  fOutputList->Add(fVertexZ);

  // Event statistics
  const char* eventStatLabels[] = { 
    "All Events",
    "Physics Selection",
    "Centrality Selection",
    "Vertex Selection",
    "Analyzed Events"
  };
  const size_t nEventStat(sizeof(eventStatLabels)/sizeof(const char*));
  TH1* hStats(new TH1D("histEventStats", "histEventStats", nEventStat, -0.5, nEventStat-0.5));
  for (size_t i(0); i<nEventStat; ++i)
    hStats->GetXaxis()->SetBinLabel(1+i, eventStatLabels[i]);
  fOutputList->Add(hStats);

  // QA histograms
  fOutputList->Add(new TH1D("histQAVertexZ", ";vertex-z (cm)",
			    800, -40, 40));
  fOutputList->Add(new TH3D("histQACentPt", ";charge;centrality V0M(%);p_{T} (GeV/c)",
			    2, -0.5, 1.5, fnBinsCent, fxMinCent, fxMaxCent, fnBinsPt, fxMinPt, fxMaxPt));
  fOutputList->Add(new TH3D("histQAPhiEta", ";charge;#phi (rad);#eta",
			    2, -0.5, 1.5, 200, 0.0, TMath::TwoPi(), 300, -1.5, 1.5));

  // N(eta) distributions with different binnings
  fOutputList->Add(new TH2D("histNEta_300", ";#eta;N", 300, -1.5, 1.5, 1000, 0., 1000.));  // 0.01
  fOutputList->Add(new TH2D("histNEta_120", ";#eta;N", 120, -1.5, 1.5, 1000, 0., 1000.));  // 0.025
  fOutputList->Add(new TH2D("histNEta__30", ";#eta;N",  30, -1.5, 1.5, 1000, 0., 1000.));  // 0.1
  fOutputList->Add(new TH2D("histNEta__15", ";#eta;N",  15, -1.5, 1.5, 1000, 0., 1000.));  // 0.2

  // Moments
  fOutputList->Add(MakeHistSparsePhiEta("histMoment1PhiEta_1"));
  fOutputList->Add(MakeHistSparsePhiEta("histMoment1PhiEta_2"));
  fOutputList->Add(MakeHistSparsePhiEtaPhiEta("histMoment2PhiEtaPhiEta_11"));
  fOutputList->Add(MakeHistSparsePhiEtaPhiEta("histMoment2PhiEtaPhiEta_12"));
  fOutputList->Add(MakeHistSparsePhiEtaPhiEta("histMoment2PhiEtaPhiEta_22"));
  
  // add MC Histograms
  const Int_t N(fOutputList->GetEntries());
  for (Int_t i(0); i<N; ++i)
    fOutputList->Add(fOutputList->At(i)->Clone(TString("MC_")
					       + fOutputList->At(i)->GetName()));

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

  AliAODMCHeader *pAODMCHeader(dynamic_cast<AliAODMCHeader*>(pAOD->FindListObject(AliAODMCHeader::StdBranchName())));
  const Bool_t isMC(NULL != pAODMCHeader);

  TClonesArray *arrayMC(NULL);
  if (isMC) {
    arrayMC = dynamic_cast<TClonesArray*>(pAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
    if (NULL == arrayMC) return;
  }

  // --- event cuts MC ---
  if (isMC) {
    Fill("MC_histEventStats", 0.); // all events
    Fill("MC_histEventStats", 1.); // events passing physics selection
    Fill("MC_histEventStats", 2.); // events passing centrality selection
    Fill("MC_histQAVertexZ", pAODMCHeader->GetVtxZ());

    // vertex cut in MC
    if (TMath::Abs(pAODMCHeader->GetVtxZ()) > fMaxAbsVertexZ) return;    
    Fill("MC_histEventStats", 3.); // events passing vertex selection
  }

  // --- event cuts data ---
  Fill("histEventStats", 0.); // all events

  if (!pInputHandler->IsEventSelected()) return;
  Fill("histEventStats", 1.); // events passing physics selection

  const Double_t centrality(pAODHeader->GetCentralityP()->GetCentralityPercentile("V0M"));
  AliDebug(3, Form("centrality=%f", centrality));
  if (centrality < fCentMin || centrality >= fCentMax) return;
  Fill("histEventStats", 2.); // events passing centrality selection

  // vertex selection
  const Int_t nVertex(pAOD->GetNumberOfVertices());
  if (0 == nVertex) return;
  const AliAODVertex* pVertex(pAOD->GetPrimaryVertex());
  if (NULL == pVertex) return;
  const Int_t nTracksPrimary(pVertex->GetNContributors());
  if (nTracksPrimary < 1) return;

  const Double_t zVertex(pVertex->GetZ());
  Fill("histQAVertexZ", zVertex);
  if (TMath::Abs(zVertex) > fMaxAbsVertexZ) return;

  Fill("histEventStats", 3.); // events passing vertex selection

  // event is accepted
  TObjArray* tracksMain(GetAcceptedTracks(pAOD, centrality));

  if (fRunMixing) {  
    AliEventPool* pEventPool(fPoolMgr->GetEventPool(centrality, pAOD->GetPrimaryVertex()->GetZ()));
    if (NULL == pEventPool)
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality, pAOD->GetPrimaryVertex()->GetZ()));
    
//     pEventPool->PrintInfo();
    if (pEventPool->IsReady()
	|| pEventPool->NTracksInPool()     > fMixingTracks/10
	|| pEventPool->GetCurrentNEvents() >= 5) {
      Fill("histEventStats", 4.); // analyzed events
      const Int_t nMix(pEventPool->GetCurrentNEvents());
      for (Int_t i(0); i<nMix; ++i) {
 	TObjArray* tracksMixed(pEventPool->GetEvent(i));
 	CalculateMoments("", tracksMain, tracksMixed, zVertex, 1./nMix);
      }
    }
    // Update the Event pool
    pEventPool->UpdatePool(tracksMain);
  } else { // no mixing
    Fill("histEventStats", 4.); // analyzed events
    CalculateMoments("", tracksMain, tracksMain, zVertex, 1.);
    delete tracksMain;
  }
  
  if (isMC) {
    TObjArray* tracksMC(GetAcceptedTracks(arrayMC, centrality));
    Fill("MC_histEventStats", 4.); // analyzed MC events
    CalculateMoments("MC_", tracksMC, tracksMC, zVertex, 1.);
    delete tracksMC;
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

  Double_t centralityBins[] = { fCentMin, fCentMax };
  const Int_t nCentralityBins(sizeof(centralityBins)/sizeof(Double_t) - 1);

  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth,
				     nCentralityBins, centralityBins,
				     fVertexZ->GetNbins()+1, 
				     const_cast<Double_t*>(fVertexZ->GetXbins()->GetArray()));
}

THnSparse* AliAnalysisTaskLongRangeCorrelations::MakeHistSparsePhiEta(const char* name) const {
  const Int_t nVertexZBins=fVertexZ->GetNbins();
  const Int_t   nBinsM1[] = { fnBinsPhi, fnBinsEta, nVertexZBins     };
  const Double_t xMinM1[] = {  fxMinPhi,  fxMinEta, 0.5              };
  const Double_t xMaxM1[] = {  fxMaxPhi,  fxMaxEta, nVertexZBins+0.5 };
  const TString title(TString(name)
		      +";#phi;#eta;vertex Z (bin);");
  return new THnSparseD(name, title.Data(), 3, nBinsM1, xMinM1, xMaxM1);
}
THnSparse* AliAnalysisTaskLongRangeCorrelations::MakeHistSparsePhiEtaPhiEta(const char* name) const {
  const Int_t nVertexZBins=fVertexZ->GetNbins();
  const Int_t   nBinsM2[] = {  fnBinsPhi, fnBinsEta,  fnBinsPhi, fnBinsEta, nVertexZBins     };
  const Double_t xMinM2[] = {  fxMinPhi,  fxMinEta,   fxMinPhi,  fxMinEta,  0.5              };
  const Double_t xMaxM2[] = {  fxMaxPhi,  fxMaxEta,   fxMaxPhi,  fxMaxEta,  nVertexZBins+0.5 };
  const TString title(TString(name)
		      +";#phi_{1};#eta_{1}"
		      +";#phi_{2};#eta_{2};vertex Z (bin);");
  return new THnSparseD(name, title.Data(), 5, nBinsM2, xMinM2, xMaxM2);
}

TObjArray* AliAnalysisTaskLongRangeCorrelations::GetAcceptedTracks(AliAODEvent* pAOD,
								   Double_t centrality) {
  TObjArray* tracks= new TObjArray;
  tracks->SetOwner(kTRUE);

  const Long64_t N(pAOD->GetNumberOfTracks());
  AliDebug(5, Form("#tracks= %6lld %f", N, centrality));
  for (Long64_t i(0); i<N; ++i) {
    AliAODTrack* pAODTrack(dynamic_cast<AliAODTrack*>(pAOD->GetTrack(i)));
    if (NULL == pAODTrack) continue;

    // track filter selection
    if (!pAODTrack->TestFilterBit(fTrackFilter)) continue;

    // select only primary tracks
    if (pAODTrack->GetType() != AliAODTrack::kPrimary) continue;

    // select only charged tracks
    if (pAODTrack->Charge() == 0) continue;

    Fill("histQACentPt", pAODTrack->Charge()>0, centrality,       pAODTrack->Pt());
    Fill("histQAPhiEta", pAODTrack->Charge()>0, pAODTrack->Phi(), pAODTrack->Eta());
    if (pAODTrack->Phi() < fPhiMin || pAODTrack->Phi() > fPhiMax) continue;
    if (pAODTrack->Pt()  < fPtMin  || pAODTrack->Pt()  > fPtMax)  continue;

    tracks->Add(new LRCParticle(pAODTrack->Eta(), pAODTrack->Phi()));
  } // next track
  return tracks;
}

TObjArray* AliAnalysisTaskLongRangeCorrelations::GetAcceptedTracks(TClonesArray* tracksMC,
								   Double_t centrality) {
  TObjArray* tracks= new TObjArray;
  tracks->SetOwner(kTRUE);

  const Long64_t N(tracksMC->GetEntriesFast());
  AliDebug(5, Form("#tracks= %6lld %f", N, centrality));
  for (Long64_t i(0); i<N; ++i) {
    AliAODMCParticle* pMCTrack(dynamic_cast<AliAODMCParticle*>(tracksMC->At(i)));
    if (NULL == pMCTrack) continue;    

    // no track filter selection for MC tracks    

    // select only primary tracks
    if (kFALSE == pMCTrack->IsPhysicalPrimary()) continue;

    // select only charged tracks
    if (pMCTrack->Charge() == 0) continue;

    Fill("MC_histQACentPt", pMCTrack->Charge()>0, centrality,      pMCTrack->Pt());
    Fill("MC_histQAPhiEta", pMCTrack->Charge()>0, pMCTrack->Phi(), pMCTrack->Eta());

    if (pMCTrack->Phi() < fPhiMin || pMCTrack->Phi() > fPhiMax) continue;
    if (pMCTrack->Pt()  < fPtMin  || pMCTrack->Pt()  > fPtMax)  continue;

    tracks->Add(new LRCParticle(pMCTrack->Eta(), pMCTrack->Phi()));
  } // next track
  return tracks;
}

void AliAnalysisTaskLongRangeCorrelations::CalculateMoments(TString prefix,
							    TObjArray* tracks1,
							    TObjArray* tracks2,
							    Double_t vertexZ,
							    Double_t weight) {
  THnSparse* hN1(dynamic_cast<THnSparse*>(fOutputList->FindObject(prefix+"histMoment1PhiEta_1")));
  THnSparse* hN2(dynamic_cast<THnSparse*>(fOutputList->FindObject(prefix+"histMoment1PhiEta_2")));
  if (NULL == hN1) return;
  if (NULL == hN2) return;

  // <n_1>
  THnSparse* hN1ForThisEvent(ComputeNForThisEvent(tracks1, "hN1", vertexZ));
  hN1->Add(hN1ForThisEvent, weight);

  // <n_2>
  THnSparse* hN2ForThisEvent(ComputeNForThisEvent(tracks2, "hN2", vertexZ));
  hN2->Add(hN2ForThisEvent, weight);

  // n(eta) distributions
  FillNEtaHist(prefix+"histNEta_300", hN1ForThisEvent, weight);
  FillNEtaHist(prefix+"histNEta_120", hN1ForThisEvent, weight);
  FillNEtaHist(prefix+"histNEta__30", hN1ForThisEvent, weight);
  FillNEtaHist(prefix+"histNEta__15", hN1ForThisEvent, weight);

  TObjArray* hNs(new TObjArray);

  // <n_1 n_1>
  hNs->AddAt(hN1ForThisEvent, 0);
  hNs->AddAt(hN1ForThisEvent, 1);
  ComputeNXForThisEvent(hNs, prefix+"histMoment2PhiEtaPhiEta_11", vertexZ, weight);

  // <n_1 n_2>
  hNs->AddAt(hN2ForThisEvent, 1);
  ComputeNXForThisEvent(hNs, prefix+"histMoment2PhiEtaPhiEta_12", vertexZ, weight);

  // <n_2 n_2>
  hNs->AddAt(hN2ForThisEvent, 0);
  ComputeNXForThisEvent(hNs, prefix+"histMoment2PhiEtaPhiEta_22", vertexZ, weight);

  // clean up
  delete hNs;
  delete hN1ForThisEvent;
  delete hN2ForThisEvent;
}

void AliAnalysisTaskLongRangeCorrelations::FillNEtaHist(TString name,
							THnSparse* hs,
							Double_t weight) {

  TH2* hSum(dynamic_cast<TH2*>(fOutputList->FindObject(name)));
  if (NULL == hSum) return;

  TH2* hPerEvent(dynamic_cast<TH2*>(hSum->Clone("hPerEvent")));
  if (NULL == hPerEvent) return;
  hPerEvent->Reset();  

  // fill hPerEvent
  const Long64_t N(hs->GetNbins());
  for (Long64_t i(0); i<N; ++i) {
    Int_t coord[2] = { 0, 0 };
    const Double_t n(hs->GetBinContent(i, coord));
    const Double_t eta(hs->GetAxis(1)->GetBinCenter(coord[1]));
    hPerEvent->Fill(eta, n);
  }

  // add zero counts for those eta bins with zero tracks
  TH1* h(hPerEvent->ProjectionX());
  for (Int_t i(1); i<=h->GetNbinsX(); ++i)
    hPerEvent->SetBinContent(i,1, Double_t(h->GetBinContent(i) == 0));

  hSum->Add(hPerEvent, weight);

  delete h;
  delete hPerEvent;
}

THnSparse* AliAnalysisTaskLongRangeCorrelations::ComputeNForThisEvent(TObjArray* tracks,
								      const char* histName,
								      Double_t vertexZ) const {
  const Double_t vertexZBin(fVertexZ->FindBin(vertexZ));
  THnSparse* hN(MakeHistSparsePhiEta(histName));
  const Long64_t nTracks(tracks->GetEntriesFast());
  for (Long64_t i(0); i<nTracks; ++i) {
    const LRCParticle* p(dynamic_cast<LRCParticle*>(tracks->At(i)));
    if (NULL == p) continue;
    const Double_t x[] = { p->Phi(), p->Eta(), vertexZBin };
    hN->Fill(x);
  }
  return hN;
}

class MultiDimIterator {
public:
  MultiDimIterator(TObjArray* _fHs, Double_t vertexZBin)
    : fHs(_fHs)
    , fN(fHs->GetEntriesFast())
    , fDims(fN,  0)
    , fIdxs(fN,  0)
    , fNs  (fN,  0)
    , fX (2*fN+1, 0)
    , fj(0) {
    for (Long64_t i=0; i<fN; ++i) {
      THnSparse* hNi(reinterpret_cast<THnSparse*>(fHs->At(i)));
      if (NULL == hNi)
	AliFatal("NULL == hNi");
      fDims[i] = hNi->GetNbins();
      AliDebug(3, Form("%lld %s %lld", i, hNi->GetName(), fDims[i]));
      update(i);
    }
    fX[2*fN] = vertexZBin;
  }

  const char* ClassName() const { return "MultiDimIterator"; }
  const Double_t* GetX() const { return &fX.front(); }
  Double_t        GetN() const { // returns the product of all multiplicities
    return std::accumulate(fNs.begin(), fNs.end(), Double_t(1), std::multiplies<Double_t>());
  }
  Bool_t end() const { return fj == fN; }

  MultiDimIterator& operator++() {
    Long64_t j(0);
    for (; j<fN; ++j) {
      ++fIdxs[j];
      update(j);
      if (fIdxs[j] < fDims[j]) break;
      fIdxs[j] = 0;
      update(j);
    }
    fj = j;
    return *this;
  }

protected:
  void update(Long64_t j) {
    THnSparse* hs(reinterpret_cast<THnSparse*>(fHs->At(j)));
    Int_t coord[] = { 0, 0 };
    fNs[j] = hs->GetBinContent(fIdxs[j], coord);
    for (Int_t k(0); k<2; ++k)
      fX[2*j+k] = hs->GetAxis(k)->GetBinCenter(coord[k]);
  }
private:
  MultiDimIterator(const MultiDimIterator&);
  MultiDimIterator& operator=(const MultiDimIterator&);

  TObjArray*            fHs;   // array of THnSparse histograms
  const Long64_t        fN;    // number of histograms
  std::vector<Long64_t> fDims; // number of filled bins for each THnSparse
  std::vector<Long64_t> fIdxs; // indices
  std::vector<Double_t> fNs;   // number of tracks
  std::vector<Double_t> fX;    // coordinate
  Long64_t              fj;    // state
} ;

void AliAnalysisTaskLongRangeCorrelations::ComputeNXForThisEvent(TObjArray* hNs,
								 const char* histName,
								 Double_t vertexZ,
								 Double_t weight) {
  if (NULL == fOutputList) return;

  THnSparse* hs(dynamic_cast<THnSparse*>(fOutputList->FindObject(histName)));
  if (hs == NULL) return;

  for (MultiDimIterator mdi(hNs, fVertexZ->FindBin(vertexZ)); !mdi.end(); ++mdi)
    hs->Fill(mdi.GetX(), mdi.GetN()*weight);
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
void AliAnalysisTaskLongRangeCorrelations::Fill(const char* histName, Double_t x, Double_t y, Double_t z) {
  if (NULL == fOutputList) return;
  TH3* h = dynamic_cast<TH3*>(fOutputList->FindObject(histName));
  if (h == NULL) return;
  h->Fill(x, y, z);
}
void AliAnalysisTaskLongRangeCorrelations::Fill(const char* histName, const Double_t* x, Double_t w) {
  if (NULL == fOutputList) return;
  THnSparse* h = dynamic_cast<THnSparse*>(fOutputList->FindObject(histName));
  if (h == NULL) return;
  h->Fill(x, w);
}
