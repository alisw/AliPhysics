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
// $Id: AliAnalysisTaskLongRangeCorrelations.cxx 407 2014-03-21 11:55:57Z cmayer $

#include <numeric>
#include <functional>
#include <set>

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
#include "AliTHn.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
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
  , fSelectPrimaryMCParticles(0)
  , fSelectPrimaryMCDataParticles(0)
  , fNMin(-1)
  , fNMax(-1)
  , fDeltaEta(-1)
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
  const Int_t nVertexBins(sizeof(vertexBins)/sizeof(Double_t)-1);
  fVertexZ = new TAxis(nVertexBins, vertexBins);
  fVertexZ->SetName("VertexZAxis");
  fOutputList->Add(fVertexZ);

  fOutputList->Add(new TH1D("histVertexZ", ";vertex Z (cm)", nVertexBins, vertexBins));

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
  fOutputList->Add(new TH2D("histQACentrality", ";centrality V0M(%);selected;",
			    fnBinsCent,fxMinCent,fxMaxCent, 2, -0.5, 1.5));
  fOutputList->Add(new TH2D("histQAVertexZ", ";vertex-z (cm);selected;",
			    800, -40., 40., 2, -0.5, 1.5));
  fOutputList->Add(new TH1D("histQAMultiplicityBeforeCuts",
			    ";tracks", 10000, 0, 10000));
  fOutputList->Add(new TH1D("histQAMultiplicityAfterCuts",
			    ";selected tracks", 10000, 0, 10000));
  fOutputList->Add(new TH3D("histQACentPt", ";charge;centrality V0M(%);p_{T} (GeV/c)",
			    2, -0.5, 1.5, fnBinsCent, fxMinCent, fxMaxCent, fnBinsPt, fxMinPt, fxMaxPt));
  fOutputList->Add(new TH3D("histQAPhiEta", ";charge;#phi (rad);#eta",
			    2, -0.5, 1.5, 200, 0.0, TMath::TwoPi(), 300, -1.5, 1.5));

  // N(eta) distributions with different binnings
  fOutputList->Add(new TH2D("histNEta_120", ";#eta;N", 120, -1.5, 1.5, 1000, -.5, 1000-.5));  // 0.025
  fOutputList->Add(new TH2D("histNEta__60", ";#eta;N",  60, -1.5, 1.5, 1000, -.5, 1000-.5));  // 0.05
  fOutputList->Add(new TH2D("histNEta__30", ";#eta;N",  30, -1.5, 1.5, 1000, -.5, 1000-.5));  // 0.1
  fOutputList->Add(new TH2D("histNEta__15", ";#eta;N",  15, -1.5, 1.5, 1000, -.5, 1000-.5));  // 0.2

  // Moments
  fOutputList->Add(MakeHistPhiEta("histMoment1PhiEta_1"));
  if (fRunMixing)
    fOutputList->Add(MakeHistPhiEta("histMoment1PhiEta_2"));
  fOutputList->Add(MakeHistPhiEtaPhiEta("histMoment2PhiEtaPhiEta_11"));
  if (fRunMixing) {
    fOutputList->Add(MakeHistPhiEtaPhiEta("histMoment2PhiEtaPhiEta_12"));
    fOutputList->Add(MakeHistPhiEtaPhiEta("histMoment2PhiEtaPhiEta_22"));
  }
  
  // add MC Histograms
  const Int_t N(fOutputList->GetEntries());
  for (Int_t i(0); i<N; ++i)
    fOutputList->Add(fOutputList->At(i)->Clone(TString("MC_") + fOutputList->At(i)->GetName()));

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

  AliAODHeader *pAODHeader = dynamic_cast<AliAODHeader*>(pAOD->GetHeader());
  if(!pAODHeader) AliFatal("Not a standard AOD");
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
  }

  // --- event cuts data ---
  Fill("histEventStats", 0.); // all events

  // physics selection
  if (!pInputHandler->IsEventSelected()) return;
  Fill("histEventStats", 1.); // events passing physics selection

  // centrality selection
  const Double_t centrality(pAODHeader->GetCentralityP()->GetCentralityPercentile("V0M"));
  AliDebug(3, Form("centrality=%f", centrality));
  const Bool_t centralitySelected(centrality > fCentMin && centrality < fCentMax);
  Fill("histQACentrality", centrality, centralitySelected);
  if (!centralitySelected) return;
  Fill("histEventStats", 2.); // events passing centrality selection
  if (isMC)
    Fill("MC_histEventStats", 2.); // events passing centrality selection

  // vertex selection -- data
  const Int_t nVertex(pAOD->GetNumberOfVertices());
  if (0 == nVertex) return;
  const AliAODVertex* pVertex(pAOD->GetPrimaryVertex());
  if (NULL == pVertex) return;
  const Int_t nTracksPrimary(pVertex->GetNContributors());
  if (nTracksPrimary < 1) return;

  const Double_t zVertex(pVertex->GetZ());
  const Bool_t vertexSelectedData(TMath::Abs(zVertex) < fMaxAbsVertexZ);

  // vertex selection -- MC
  Bool_t vertexSelectedMC(kTRUE);
  if (isMC)
    vertexSelectedMC = (TMath::Abs(pAODMCHeader->GetVtxZ()) < fMaxAbsVertexZ);

  // combined vertex selection (data and MC)
  const Bool_t vertexSelected(vertexSelectedData && vertexSelectedMC);
  Fill("histQAVertexZ", zVertex, vertexSelected);
  if (isMC)
    Fill("MC_histQAVertexZ", pAODMCHeader->GetVtxZ(), vertexSelected);
  if (!vertexSelected) return;
  Fill("histEventStats",    3.); // events passing vertex selection
  if (isMC)
    Fill("MC_histEventStats", 3.); // events passing vertex selection

  // ------------------
  // event is accepted
  // ------------------

  TObjArray* tracksMain(GetAcceptedTracks(pAOD, arrayMC, centrality));
  Fill("histQAMultiplicityBeforeCuts", pAOD->GetNumberOfTracks());
  Fill("histQAMultiplicityAfterCuts", tracksMain->GetEntriesFast());
  
  if (fRunMixing) {  
    AliEventPool* pEventPool(fPoolMgr->GetEventPool(centrality, pAOD->GetPrimaryVertex()->GetZ()));
    if (NULL == pEventPool)
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality, pAOD->GetPrimaryVertex()->GetZ()));
    
//     pEventPool->PrintInfo();
    if (pEventPool->IsReady()
	|| pEventPool->NTracksInPool()     > fMixingTracks/10
	|| pEventPool->GetCurrentNEvents() >= 5) {
      Fill("histEventStats", 4.); // analyzed events
      Fill("histVertexZ", zVertex);
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
    Fill("histVertexZ", zVertex);
    CalculateMoments("", tracksMain, tracksMain, zVertex, 1.);
    delete tracksMain;
  }
  
  if (isMC) {
    TObjArray* tracksMC(GetAcceptedTracks(arrayMC, centrality));
    const Double_t x[2] = {
      static_cast<Double_t>(arrayMC->GetEntriesFast()),
      static_cast<Double_t>(tracksMC->GetEntriesFast())
    };
    Fill("MC_histQAMultiplicity", x);
    Fill("MC_histEventStats", 4.); // analyzed MC events
    Fill("MC_histVertexZ", pAOD->GetPrimaryVertex()->GetZ());
    CalculateMoments("MC_", tracksMC, tracksMC, zVertex, 1.);
    delete tracksMC;
  }
}

void AliAnalysisTaskLongRangeCorrelations::Terminate(Option_t* ) {
  //
  fOutputList = dynamic_cast<TList*>(GetOutputData(1));
  if (NULL == fOutputList) {
    AliFatal("NULL == fOutputList");
    return; // needed to avoid coverity warning
  }
}

TString AliAnalysisTaskLongRangeCorrelations::GetOutputListName() const {
  TString listName("listLRC");
  listName += TString::Format("_%smix",         fRunMixing ? "" : "no");
  listName += TString::Format("_trackFilt%d",   fTrackFilter);
  listName += TString::Format("_cent%.0fT%.0f", fCentMin, fCentMax);
  listName += TString::Format("_ptMin%.0fMeV",  1e3*fPtMin);
  listName += TString::Format("_phi%.0fT%.0f",  TMath::RadToDeg()*fPhiMin, TMath::RadToDeg()*fPhiMax);
  if ( 1 == fSelectPrimaryMCParticles)
    listName += "_selPrimMC";
  if (-1 == fSelectPrimaryMCParticles)
    listName += "_selNonPrimMC";
  if ( 1 == fSelectPrimaryMCDataParticles)
    listName += "_selPrimMCData";
  if (-1 == fSelectPrimaryMCDataParticles)
    listName += "_selNonPrimMCData";
  if (-1 != fNMin)
    listName += TString::Format("_nMin%d", fNMin);
  if (-1 != fNMax)
    listName += TString::Format("_nMax%d", fNMax);
  if (fDeltaEta >= 0)
    listName += TString::Format("_deltaEta%02.0f", 10*fDeltaEta);
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

AliTHn* AliAnalysisTaskLongRangeCorrelations::MakeHistPhiEta(const char* name) const {
  const Int_t nVertexZBins=fVertexZ->GetNbins();
  const Int_t   nBinsM1[] = { fnBinsPhi, fnBinsEta, nVertexZBins     };
  const Double_t xMinM1[] = {  fxMinPhi,  fxMinEta, 0.5              };
  const Double_t xMaxM1[] = {  fxMaxPhi,  fxMaxEta, nVertexZBins+0.5 };
  const TString title(TString(name)
		      +";#phi;#eta;vertex Z (bin);");
  AliTHn* h(new AliTHn(name, title.Data(), 1, 3, nBinsM1));
  for (Int_t i=0; i<3; ++i)
    h->SetBinLimits(i, xMinM1[i], xMaxM1[i]);
  return h;
}

AliTHn* AliAnalysisTaskLongRangeCorrelations::MakeHistPhiEtaPhiEta(const char* name) const {
  const Int_t nVertexZBins=fVertexZ->GetNbins();
  const Int_t   nBinsM2[] = {  fnBinsPhi, fnBinsEta,  fnBinsPhi, fnBinsEta, nVertexZBins     };
  const Double_t xMinM2[] = {  fxMinPhi,  fxMinEta,   fxMinPhi,  fxMinEta,  0.5              };
  const Double_t xMaxM2[] = {  fxMaxPhi,  fxMaxEta,   fxMaxPhi,  fxMaxEta,  nVertexZBins+0.5 };
  const TString title(TString(name)
		      +";#phi_{1};#eta_{1}"
		      +";#phi_{2};#eta_{2};vertex Z (bin);");
  AliTHn* h(new AliTHn(name, title.Data(), 1, 5, nBinsM2));
  for (Int_t i=0; i<5; ++i)
    h->SetBinLimits(i, xMinM2[i], xMaxM2[i]);
  return h;
}

TObjArray* AliAnalysisTaskLongRangeCorrelations::GetAcceptedTracks(AliAODEvent*  pAOD,
								   TClonesArray* arrayMC,
								   Double_t      centrality) {
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

    if (NULL != arrayMC) {
      const Int_t label(pAODTrack->GetLabel());
      AliAODMCParticle* mcParticle((label >= 0) 
				   ? static_cast<AliAODMCParticle*>(arrayMC->At(label))
				   : NULL);
      if (label >=0 && NULL == mcParticle)
	AliFatal("MC particle not found");

//       const Bool_t isPhysicalPrimary(mcParticle->IsPhysicalPrimary());
//       Printf("mcParticle->IsPhysicalPrimary() %d", isPhysicalPrimary);

      switch (fSelectPrimaryMCDataParticles) {
      case -1:
	if (label < 0) continue;
	if (kTRUE  == mcParticle->IsPhysicalPrimary()) continue;
	break;
      case  0:
	// NOP, take all tracks
	break;
      case  1:
	if (label < 0) continue;
	if (kFALSE == mcParticle->IsPhysicalPrimary()) continue;
	break;
      default:
	AliFatal("fSelectPrimaryMCDataParticles != {-1,0,1}");
      }            
    }

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
  // for keeping track of MC labels
  std::set<Int_t> labelSet;

  TObjArray* tracks= new TObjArray;
  tracks->SetOwner(kTRUE);

  const Long64_t N(tracksMC->GetEntriesFast());
  AliDebug(5, Form("#tracks= %6lld %f", N, centrality));
  for (Long64_t i(0); i<N; ++i) {
    AliAODMCParticle* pMCTrack(dynamic_cast<AliAODMCParticle*>(tracksMC->At(i)));
    if (NULL == pMCTrack) continue;    

    // no track filter selection for MC tracks    

    if (labelSet.find(pMCTrack->Label()) != labelSet.end()) {
      Printf("Duplicate Label= %3d", pMCTrack->Label());
      continue;
    }
    labelSet.insert(pMCTrack->Label());    
    
//     Printf("isPrim = %d", pMCTrack->IsPhysicalPrimary());

    switch (fSelectPrimaryMCParticles) {
    case -1:
      if (kTRUE == pMCTrack->IsPhysicalPrimary()) continue;
      break;
    case  0:
      // NOP, take all MC tracks
      break;
    case  1:
      if (kFALSE == pMCTrack->IsPhysicalPrimary()) continue;
      break;
    default:
      AliFatal("fSelectPrimaryMCParticles != {-1,0,1}");
    }

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

Bool_t AddTHnSparseToAliTHn(AliTHn* h, THnSparse* hs, Double_t weight) {
  if (h->GetNVar() != hs->GetNdimensions())
    return kFALSE;

  const size_t nDim(hs->GetNdimensions());
  
  Int_t    *coord = new Int_t[nDim];
  Double_t *x     = new Double_t[nDim];

  const Long64_t N(hs->GetNbins());
  for (Long64_t i(0); i<N; ++i) {
    const Double_t n(hs->GetBinContent(i, coord));
    for (size_t j=0; j<nDim; ++j) 
      x[j] = hs->GetAxis(j)->GetBinCenter(coord[j]);
    h->Fill(x, 0, n*weight);
  }

  delete[] coord;
  delete[] x;

  return kTRUE;
}

void AliAnalysisTaskLongRangeCorrelations::CalculateMoments(TString prefix,
							    TObjArray* tracks1,
							    TObjArray* tracks2,
							    Double_t vertexZ,
							    Double_t weight) {
  THnSparse* hN1ForThisEvent(ComputeNForThisEvent(tracks1, "hN1", vertexZ));

  if (fDeltaEta >= 0) {
    hN1ForThisEvent->GetAxis(1)->SetRangeUser( fDeltaEta/2+0.00001,  fDeltaEta/2+0.19999);
    TH1 *hTemp = hN1ForThisEvent->Projection(0);
    const Long64_t ncPlus = Long64_t(hTemp->GetEntries());
    delete hTemp;
    
    hN1ForThisEvent->GetAxis(1)->SetRangeUser(-fDeltaEta/2-0.19999, -fDeltaEta/2-0.00001);
    hTemp = hN1ForThisEvent->Projection(0);
    const Long64_t ncMinus = Long64_t(hTemp->GetEntries());
    delete hTemp;
    
    // restore full axis range
    hN1ForThisEvent->GetAxis(1)->SetRange(0, -1);

    if (fNMin != -1 && ncPlus < fNMin) return;
    if (fNMax != -1 && ncPlus > fNMax) return;
    
    if (fNMin != -1 && ncMinus < fNMin) return;
    if (fNMax != -1 && ncMinus > fNMax) return;
  }

  AliTHn* hN1(dynamic_cast<AliTHn*>(fOutputList->FindObject(prefix+"histMoment1PhiEta_1")));
  if (NULL == hN1) return;

  // <n_1>
  AddTHnSparseToAliTHn(hN1,hN1ForThisEvent, weight); 
//   hN1->GetGrid(0)->GetGrid()->Add(hN1ForThisEvent, weight);

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

  if (fRunMixing) {
    AliTHn* hN2(dynamic_cast<AliTHn*>(fOutputList->FindObject(prefix+"histMoment1PhiEta_2")));
    if (NULL == hN2) return;
    // <n_2>
    THnSparse* hN2ForThisEvent(ComputeNForThisEvent(tracks2, "hN2", vertexZ));
    AddTHnSparseToAliTHn(hN2,hN2ForThisEvent, weight); 
//     hN2->GetGrid(0)->GetGrid()->Add(hN2ForThisEvent, weight);

    // <n_1 n_2>
    hNs->AddAt(hN2ForThisEvent, 1);
    ComputeNXForThisEvent(hNs, prefix+"histMoment2PhiEtaPhiEta_12", vertexZ, weight);
    
    // <n_2 n_2>
    hNs->AddAt(hN2ForThisEvent, 0);
    ComputeNXForThisEvent(hNs, prefix+"histMoment2PhiEtaPhiEta_22", vertexZ, weight);

    // clean up
    delete hN2ForThisEvent;
  }

  // clean up
  delete hNs;
  delete hN1ForThisEvent;
}

void AliAnalysisTaskLongRangeCorrelations::FillNEtaHist(TString name,
							THnSparse* hs,
							Double_t weight) {

  TH2* hSum(dynamic_cast<TH2*>(fOutputList->FindObject(name)));
  if (NULL == hSum) return;

  TH2* hPerEvent(dynamic_cast<TH2*>(hSum->Clone("hPerEvent")));
  if (NULL == hPerEvent) return;
  hPerEvent->Reset();  

  TH1* h1PerEvent(hPerEvent->ProjectionX());

  // fill h1PerEvent
  const Long64_t N(hs->GetNbins());
  for (Long64_t i(0); i<N; ++i) {
    Int_t coord[2] = { 0, 0 };
    const Double_t n(hs->GetBinContent(i, coord));
    const Double_t eta(hs->GetAxis(1)->GetBinCenter(coord[1]));
    h1PerEvent->Fill(eta, n);
  }
 
  for (Int_t i(1); i<=h1PerEvent->GetNbinsX(); ++i)
    hPerEvent->Fill(h1PerEvent->GetXaxis()->GetBinCenter(i),
		    h1PerEvent->GetBinContent(i));

  hSum->Add(hPerEvent, weight);

  delete hPerEvent;
  delete h1PerEvent;
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

  AliTHn* hs(dynamic_cast<AliTHn*>(fOutputList->FindObject(histName)));
  if (hs == NULL) return;

  for (MultiDimIterator mdi(hNs, fVertexZ->FindBin(vertexZ)); !mdi.end(); ++mdi)
    hs->Fill(mdi.GetX(), 0, mdi.GetN()*weight);

//   hs->FillParent();
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
