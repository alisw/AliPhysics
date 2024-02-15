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
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Analysis task to produce trees of lightweight events
// evgeny.kryshen@cern.ch

// aliroot
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliExternalTrackParam.h"
#include "AliCentrality.h"
#include "AliAnalysisFilter.h"
#include "AliVMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliEventPoolManager.h"
#include "AliVParticle.h"
#include "AliCFParticle.h"
// root
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TChain.h"
#include "TTree.h"

#include "AliAnalysisTaskThreePartCorr.h"

ClassImp(AliAnalysisTaskThreePartCorr)

AliAnalysisTaskThreePartCorr::AliAnalysisTaskThreePartCorr() :
AliAnalysisTaskSE(),
  fESD(0x0),
  fAOD(0x0),
  fPIDResponse(0x0),
  fUtils(0x0),
  fListOfHistos(0x0),
  fQAList(0x0),
  fCorrList(0x0),
  hEventStatistics(0x0),
  Centaxis(0x0),
  Zvtxaxis(0x0),
  Triggptaxis(0x0),
  Ptaxis(0x0),
  LambdaInvMassaxis(0x0),
  XiInvMassaxis(0x0),
  fPoolMgr(0x0),
  hTrackPt(0x0),
  hTrackEta(0x0),
  hTrackPhi(0x0),
  hEventCent(0x0),
  hEventZvtx(0x0),
  hTrackdEdx(0x0),
  hNSigmaPion(0x0),
  hNSigmaKaon(0x0),
  hNSigmaProton(0x0),
  hInvMassLambda(0x0),
  hInvMassXi(0x0),
  hSameLambda_SGNL(0x0),
  hSameLambda_SB(0x0),
  hSameXi(0x0), 
  hMixLambda_SGNL(0x0),
  hMixLambda_SB(0x0),
  hMixXi(0x0),
  fCentV0M(-1),
  fIsAOD(kTRUE),
  Nmaxmixevents(-1),
  Nmaxmixtracks(50000)
{
  // empty constructor
}

//==================================================================================================================================================================================================================

AliAnalysisTaskThreePartCorr::AliAnalysisTaskThreePartCorr(const char *name) :
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fAOD(0x0),
  fPIDResponse(0x0),
  fUtils(0x0),
  fListOfHistos(0x0),
  fQAList(0x0),
  fCorrList(0x0),
  hEventStatistics(0x0),
  Centaxis(0x0),
  Zvtxaxis(0x0),
  Triggptaxis(0x0),
  Ptaxis(0x0),
  LambdaInvMassaxis(0x0),
  XiInvMassaxis(0x0),
  fPoolMgr(0x0),
  hTrackPt(0x0),
  hTrackEta(0x0),
  hTrackPhi(0x0),
  hEventCent(0x0),
  hEventZvtx(0x0),
  hTrackdEdx(0x0),
  hNSigmaPion(0x0),
  hNSigmaKaon(0x0),
  hNSigmaProton(0x0),
  hInvMassLambda(0x0),
  hInvMassXi(0x0),
  hSameLambda_SGNL(0x0),
  hSameLambda_SB(0x0),
  hSameXi(0x0), 
  hMixLambda_SGNL(0x0),
  hMixLambda_SB(0x0),
  hMixXi(0x0),
  fCentV0M(-1),
  fIsAOD(kTRUE),
  Nmaxmixevents(-1),
  Nmaxmixtracks(5000)
{
  Info("AliAnalysisTaskThreePartCorr", "Calling Constructor");
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//==================================================================================================================================================================================================================

AliAnalysisTaskThreePartCorr::~AliAnalysisTaskThreePartCorr()
{
  if(fListOfHistos) { delete fListOfHistos; }
}

//==================================================================================================================================================================================================================

void AliAnalysisTaskThreePartCorr::UserCreateOutputObjects() {  
  fListOfHistos = new TList();
  fListOfHistos->SetOwner();

  fQAList = new TList();
  fCorrList = new TList();
  fQAList->SetName("QAList");
  fCorrList->SetName("CorrList");
  
  TList *fSEHistos = new TList();
  TList *fMEHistos = new TList();

  hEventStatistics = new TH1I("hEventStatistics", "", 10, 0, 10);
  //hEventStatistics->SetBit(TH1::kCanRebin);
  fListOfHistos->Add(hEventStatistics);

  const Int_t Ncentbins = 9, Nzvtxbins = 8, Ntriggptbins = 3;
  Double_t Centbins[Ncentbins+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
  Double_t Zvtxbins[Nzvtxbins+1] = {-7, -5, -3, -1, 0, 1, 3, 5, 7};
  Double_t Triggptbins[Ntriggptbins+1] = {0, 4, 8, 12};
  Centaxis = new TAxis(Ncentbins, Centbins);
  Zvtxaxis = new TAxis(Nzvtxbins, Zvtxbins);
  Triggptaxis = new TAxis(Ntriggptbins, Triggptbins);
  fPoolMgr = new AliEventPoolManager(Nmaxmixevents, Nmaxmixtracks, Ncentbins, Centbins, Nzvtxbins, Zvtxbins);
  fPoolMgr->SetTargetValues(Nmaxmixtracks, 0.1, 5);

  const Int_t Nptbins = 120;
  Double_t Ptbins[Nptbins+1];
  Ptaxis = new TAxis(Nptbins, 0, 12);
  for(Int_t i = 0; i < Nptbins+1; i++) { Ptbins[i] = Ptaxis->GetBinLowEdge(i+1); }
  hTrackPt = new TH1D("hTrackPt", "hTrackPt", 100, 0, 4);
  hTrackEta = new TH1D("hTrackEta", "hTrackEta", 100, -1, 1);
  hTrackPhi = new TH1D("hTrackPhi", "hTrackPhi", 100, (-1./2)*TMath::Pi(), (5./2)*TMath::Pi());
  hEventCent = new TH1D("hEventCentrality", "hEventCentrality", Ncentbins, Centbins);
  hEventZvtx = new TH1D("hEventZvtx", "hEventZvtx", Nzvtxbins, Zvtxbins);
  hTrackdEdx = new TH2D("hTrackdEdx", "hTrackdEdx", 200, 0.1, 20, 200, 20, 600);
  hNSigmaPion = new TH2D("hNSigmaPion", "hNSigmaPion", 5, 0.5, 3, 40, -4, 4);
  hNSigmaKaon = new TH2D("hNSigmaKaon", "hNSigmaKaon", 5, 0.5, 3, 40, -4, 4);
  hNSigmaProton = new TH2D("hNSigmaProton", "hNSigmaProton", 5, 0.5, 3, 40, -4, 4);
  fQAList->Add(hTrackPt);
  fQAList->Add(hTrackEta);
  fQAList->Add(hTrackPhi);
  fQAList->Add(hEventCent);
  fQAList->Add(hEventZvtx);
  fQAList->Add(hTrackdEdx);
  fQAList->Add(hNSigmaPion);
  fQAList->Add(hNSigmaKaon);
  fQAList->Add(hNSigmaProton);

  const Int_t Ninvmassbins = 100;
  Double_t LambdaInvMassbins[Ninvmassbins+1];
  Double_t XiInvMassbins[(2*Ninvmassbins)+1];
  LambdaInvMassaxis = new TAxis(Ninvmassbins, 1.08, 1.16);
  XiInvMassaxis = new TAxis(2*Ninvmassbins, 1.28, 1.36);
  for(Int_t i = 0; i < Ninvmassbins+1; i++) { LambdaInvMassbins[i] = LambdaInvMassaxis->GetBinLowEdge(i+1); }
  for(Int_t i = 0; i < (2*Ninvmassbins)+1; i++) { XiInvMassbins[i] = XiInvMassaxis->GetBinLowEdge(i+1); }  
  hInvMassLambda = new TH3D("hInvMassLambda", "hInvMassLambda", Ninvmassbins, LambdaInvMassbins, Nptbins, Ptbins, Ncentbins, Centbins);
  hInvMassXi = new TH3D("hInvMassXi", "hInvMassXi", 2*Ninvmassbins, XiInvMassbins, Nptbins, Ptbins, Ncentbins, Centbins);
  fQAList->Add(hInvMassLambda);
  fQAList->Add(hInvMassXi);

  const char* Track;
  hSameLambda_SGNL = new TH3D****[3];
  hSameLambda_SB = new TH3D****[3];
  hSameXi = new TH3D****[3];
  hMixLambda_SGNL = new TH3D****[3];
  hMixLambda_SB = new TH3D****[3];
  hMixXi = new TH3D****[3];  
  for(Int_t A = 0; A < 3; A++) { // Associates
    if(A == 0) { Track = "#pi"; }
    else if(A == 1) { Track = "K"; }
    else if(A == 2) { Track = "p"; }    
    hSameLambda_SGNL[A] = new TH3D***[Ncentbins];
    hSameLambda_SB[A] = new TH3D***[Ncentbins];
    hSameXi[A] = new TH3D***[Ncentbins];
    hMixLambda_SGNL[A] = new TH3D***[Ncentbins];
    hMixLambda_SB[A] = new TH3D***[Ncentbins];
    hMixXi[A] = new TH3D***[Ncentbins];
    
    for(Int_t C = 0; C < Ncentbins; C++) { // Centrality bins
      hSameLambda_SGNL[A][C] = new TH3D**[Nzvtxbins];
      hSameLambda_SB[A][C] = new TH3D**[Nzvtxbins];
      hSameXi[A][C] = new TH3D**[Nzvtxbins];
      hMixLambda_SGNL[A][C] = new TH3D**[Nzvtxbins];
      hMixLambda_SB[A][C] = new TH3D**[Nzvtxbins];
      hMixXi[A][C] = new TH3D**[Nzvtxbins];
      
      for(Int_t Z = 0; Z < Nzvtxbins; Z++) { // Z-vertex bins
	hSameLambda_SGNL[A][C][Z] = new TH3D*[Ntriggptbins];
	hSameLambda_SB[A][C][Z] = new TH3D*[Ntriggptbins];
	hSameXi[A][C][Z] = new TH3D*[Ntriggptbins];
	hMixLambda_SGNL[A][C][Z] = new TH3D*[Ntriggptbins];
	hMixLambda_SB[A][C][Z] = new TH3D*[Ntriggptbins];
	hMixXi[A][C][Z] = new TH3D*[Ntriggptbins];
	
	for(Int_t P = 0; P < Ntriggptbins; P++) { // Trigger transverse momentum bins	
	  hSameLambda_SGNL[A][C][Z][P] = new TH3D(Form("hSameLambda_SGNL[%i][%i][%i][%i]", A, C, Z, P), Form("Same-event #Lambda - %s correlator (SGNL region)", Track), 36, (-1./2)*TMath::Pi(), (3./2)*TMath::Pi(), 32, -1.52, 1.52, 5, 0.5, 3);
	  hSameLambda_SB[A][C][Z][P] = new TH3D(Form("hSameLambda_SB[%i][%i][%i][%i]", A, C, Z, P), Form("Same-event #Lambda - %s correlator (SB region)", Track), 36, (-1./2)*TMath::Pi(), (3./2)*TMath::Pi(), 32, -1.52, 1.52, 5, 0.5, 3);
	  hSameXi[A][C][Z][P] = new TH3D(Form("hSameXi[%i][%i][%i][%i]", A, C, Z, P), Form("Same-event #Xi - %s correlator", Track), 36, (-1./2)*TMath::Pi(), (3./2)*TMath::Pi(), 30, -1.5, 1.5, 5, 0.5, 3);	  
	  hMixLambda_SGNL[A][C][Z][P] = new TH3D(Form("hMixLambda_SGNL[%i][%i][%i][%i]", A, C, Z, P), Form("Mixed-event #Lambda - %s correlator (SGNL region)", Track), 36, (-1./2)*TMath::Pi(), (3./2)*TMath::Pi(), 32, -1.52, 1.52, 5, 0.5, 3);
	  hMixLambda_SB[A][C][Z][P] = new TH3D(Form("hMixLambda_SB[%i][%i][%i][%i]", A, C, Z, P), Form("Mixed-event #Lambda - %s correlator (SB region)", Track), 36, (-1./2)*TMath::Pi(), (3./2)*TMath::Pi(), 32, -1.52, 1.52, 5, 0.5, 3);
	  hMixXi[A][C][Z][P] = new TH3D(Form("hMixXi[%i][%i][%i][%i]", A, C, Z, P), Form("Mixed-event #Xi - %s correlator", Track), 36, (-1./2)*TMath::Pi(), (3./2)*TMath::Pi(), 30, -1.5, 1.5, 5, 0.5, 3);

	  fSEHistos->Add(hSameLambda_SGNL[A][C][Z][P]);
	  fSEHistos->Add(hSameLambda_SB[A][C][Z][P]);
	  fSEHistos->Add(hSameXi[A][C][Z][P]);
	  fMEHistos->Add(hMixLambda_SGNL[A][C][Z][P]);
	  fMEHistos->Add(hMixLambda_SB[A][C][Z][P]);
	  fMEHistos->Add(hMixXi[A][C][Z][P]);
	}
      }
    }
  }
  
  fCorrList->AddAll(fSEHistos);
  fCorrList->AddAll(fMEHistos);
  fSEHistos->Clear();
  fMEHistos->Clear();
  fListOfHistos->Add(fQAList);
  fListOfHistos->Add(fCorrList);
  
  //fUtils = new AliAnalysisUtils();
  //fUtils->SetMinPlpContribSPD(3);
  //fUtils->SetMinPlpContribMV(3);

  PostData(1, fListOfHistos);
}

//==================================================================================================================================================================================================================

void AliAnalysisTaskThreePartCorr::UserExec(Option_t *) {
  hEventStatistics->Fill("before cuts", 1);
  
  if(!fInputEvent) { return; }
  hEventStatistics->Fill("after event check", 1);

  if(fIsAOD) { // AOD
    fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
    if(!fAOD) { return; }
    hEventStatistics->Fill("after aod check", 1); }
  else { // ESD
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if(!fESD) { return; }
    hEventStatistics->Fill("after esd check", 1); }
  
  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse) { return; }
  hEventStatistics->Fill("after pid check", 1);
  
  if(!(fInputHandler->IsEventSelected())) { return; }
  hEventStatistics->Fill("physics selection", 1);

  AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent->FindListObject("MultSelection");
  if(!MultSelection) { return; }
  hEventStatistics->Fill("found MultSelection object", 1);

  Double_t fVtx[3];
  AliVVertex *Vertex = (AliVVertex*)fInputEvent->GetPrimaryVertex();
  if(!Vertex) { return; }
  fVtx[0] = Vertex->GetX();
  fVtx[1] = Vertex->GetY();
  fVtx[2] = Vertex->GetZ();
  hEventStatistics->Fill("found primary vertex", 1);
  if(TMath::Abs(fVtx[2]) > 7) { return; }
  hEventStatistics->Fill("z-vertex selection", 1);

  fCentV0M = MultSelection->GetMultiplicityPercentile("V0M");
  if(fCentV0M < 0 || fCentV0M > 90) { return; }
  hEventStatistics->Fill("centrality selection", 1);
  
  // Loop over reconstructed tracks
  Int_t nTracks = 0, nV0s = 0, nXis = 0; 
  if(fIsAOD) {
    nTracks = fAOD->GetNumberOfTracks();
    nV0s = fAOD->GetNumberOfV0s();
    nXis = fAOD->GetNumberOfCascades(); } // AOD
  else {
    nTracks = fESD->GetNumberOfTracks();
    nV0s = fESD->GetNumberOfV0s();
    nXis = fESD->GetNumberOfCascades(); } // ESD

  Int_t Centindex = (Centaxis->FindBin(fCentV0M)) - 1;
  Int_t Zvtxindex = (Zvtxaxis->FindBin(fVtx[2])) - 1;
  Int_t Ptindex;
  AliEventPool *Pool = fPoolMgr->GetEventPool(fCentV0M,fVtx[2]);

  TObjArray *TrackList = new TObjArray();
  TObjArray *ReducedList = new TObjArray();
  TObjArray *MixTracks = 0x0;
  AliAODTrack *aodtrack = 0x0;
  AliAODTrack *aodassociate = 0x0;
  AliAODTrack *aodmix = 0x0;
  AliAODv0 *aodV0trigger = 0x0;
  AliAODcascade *aodXitrigger = 0x0;

  Double_t XiPhi = 0, XiEta = 0;
  Double_t DeltaPhi = 0, DeltaRap = 0;

  hEventCent->Fill(fCentV0M);
  hEventZvtx->Fill(fVtx[2]);
  
  //================================================================================================================================================================================================================

  if(fIsAOD) { // AOD
    for(Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
      aodtrack = (AliAODTrack*)fAOD->GetTrack(iTrack);
      if(AssociateCuts(aodtrack)) {
	TrackList->Add(aodtrack);
	hTrackPt->Fill(aodtrack->Pt());
	hTrackEta->Fill(aodtrack->Eta());
	hTrackPhi->Fill(aodtrack->Phi());
	if(TrackPID(aodtrack) != AliPID::kUnknown) { hTrackdEdx->Fill(aodtrack->P(), fPIDResponse->GetTPCsignalTunedOnData(aodtrack)); }
	if(TrackPID(aodtrack) == AliPID::kPion) { hNSigmaPion->Fill(aodtrack->Pt(), fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kPion)); }
	else if(TrackPID(aodtrack) == AliPID::kKaon) { hNSigmaKaon->Fill(aodtrack->Pt(), fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kKaon)); }
	else if(TrackPID(aodtrack) == AliPID::kProton) { hNSigmaProton->Fill(aodtrack->Pt(), fPIDResponse->NumberOfSigmasTPC(aodtrack, AliPID::kProton)); }
      }
    }

    ReducedList = CloneAndReduceTrackList(TrackList);
    Pool->UpdatePool(ReducedList);	

    // Lambda correlations 
    for(Int_t iV0 = 0; iV0 < nV0s; iV0++) {         
      aodV0trigger = fAOD->GetV0(iV0);
      if(V0TriggerCuts(aodV0trigger, fVtx)) {

	Ptindex = Triggptaxis->FindBin(aodV0trigger->Pt()) - 1;
	hInvMassLambda->Fill(aodV0trigger->MassLambda(), aodV0trigger->Pt(), fCentV0M);
	for(Int_t iAssociate = 0; iAssociate < ReducedList->GetEntriesFast(); iAssociate++) {
	  aodassociate = (AliAODTrack*)ReducedList->UncheckedAt(iAssociate);

	  DeltaPhi = PhiCorrection(aodV0trigger->Phi(), aodassociate->Phi());
	  DeltaRap = aodV0trigger->Eta() - aodassociate->Eta(); // Pseudorapidity instead of rapidity
	  if(aodV0trigger->MassLambda() >= (massLambda - 3*DGaussSigma) && aodV0trigger->MassLambda() <= (massLambda + 3*DGaussSigma)) {
	    if(TrackPID(aodassociate, kTRUE) == AliPID::kPion) { hSameLambda_SGNL[0][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodassociate->Pt()); }	      
	    else if(TrackPID(aodassociate, kTRUE) == AliPID::kKaon) { hSameLambda_SGNL[1][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodassociate->Pt()); }	    
	    else if(TrackPID(aodassociate, kTRUE) == AliPID::kProton) { hSameLambda_SGNL[2][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodassociate->Pt()); }	    
	  }	  
	  else {
	    if(TrackPID(aodassociate, kTRUE) == AliPID::kPion) { hSameLambda_SB[0][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodassociate->Pt()); }	      
	    else if(TrackPID(aodassociate, kTRUE) == AliPID::kKaon) { hSameLambda_SB[1][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodassociate->Pt()); }	    
	    else if(TrackPID(aodassociate, kTRUE) == AliPID::kProton) { hSameLambda_SB[2][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodassociate->Pt()); }
	  }
	}
	
	if(Pool->IsReady()) {
	  for(Int_t iMixEvent = 0; iMixEvent < Pool->GetCurrentNEvents(); iMixEvent++) {
	    MixTracks = Pool->GetEvent(iMixEvent);
	    for(Int_t iMixTrack = 0; iMixTrack < MixTracks->GetEntriesFast(); iMixTrack++) {
	      aodmix = (AliAODTrack*)MixTracks->UncheckedAt(iMixTrack);	      
		
	      DeltaPhi = PhiCorrection(aodV0trigger->Phi(), aodmix->Phi());
	      DeltaRap = aodV0trigger->Eta() - aodmix->Eta(); // Pseudorapidity instead of rapidity	      
	      if(aodV0trigger->MassLambda() >= (massLambda - 3*DGaussSigma) && aodV0trigger->MassLambda() <= (massLambda + 3*DGaussSigma)) {		
		if(TrackPID(aodmix, kTRUE) == AliPID::kPion) { hMixLambda_SGNL[0][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodmix->Pt()); }		  
		else if(TrackPID(aodmix, kTRUE) == AliPID::kKaon) { hMixLambda_SGNL[1][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodmix->Pt()); }
		else if(TrackPID(aodmix, kTRUE) == AliPID::kProton) { hMixLambda_SGNL[2][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodmix->Pt()); }
	      }	      
	      else {
		if(TrackPID(aodmix, kTRUE) == AliPID::kPion) { hMixLambda_SB[0][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodmix->Pt()); }		  
		else if(TrackPID(aodmix, kTRUE) == AliPID::kKaon) { hMixLambda_SB[1][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodmix->Pt()); }
		else if(TrackPID(aodmix, kTRUE) == AliPID::kProton) { hMixLambda_SB[2][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodmix->Pt()); }
	      }
	    }
	  }
	}
      }
    }

    // Xi correlations
    for(Int_t iXi = 0; iXi < nXis; iXi++) {
      aodXitrigger = fAOD->GetCascade(iXi);
      if(XiTriggerCuts(aodXitrigger)) {

	Ptindex = Triggptaxis->FindBin(TMath::Sqrt(aodXitrigger->Pt2Xi())) - 1;
	hInvMassXi->Fill(aodXitrigger->MassXi(), TMath::Sqrt(aodXitrigger->Pt2Xi()), fCentV0M);
	XiPhi = TMath::ATan((aodXitrigger->MomXiY())/(aodXitrigger->MomXiX()));
	XiEta = TMath::ATanH((aodXitrigger->MomXiZ())/(TMath::Sqrt(aodXitrigger->Ptot2Xi())));
	for(Int_t iAssociate = 0; iAssociate < ReducedList->GetEntriesFast(); iAssociate++) {
	  aodassociate = (AliAODTrack*)ReducedList->UncheckedAt(iAssociate);

	  DeltaPhi = PhiCorrection(XiPhi, aodassociate->Phi());
	  DeltaRap = XiEta - aodassociate->Eta(); // Pseudorapidity instead of rapidity
	  if(TrackPID(aodassociate, kTRUE) == AliPID::kPion) { hSameXi[0][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodassociate->Pt()); }
	  else if(TrackPID(aodassociate, kTRUE) == AliPID::kKaon) { hSameXi[1][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodassociate->Pt()); }	  
	  else if(TrackPID(aodassociate, kTRUE) == AliPID::kProton) { hSameXi[2][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodassociate->Pt()); }
	}
	
	if(Pool->IsReady()) {
	  for(Int_t iMixEvent = 0; iMixEvent < Pool->GetCurrentNEvents(); iMixEvent++) {
	    MixTracks = Pool->GetEvent(iMixEvent);
	    for(Int_t iMixTrack = 0; iMixTrack < MixTracks->GetEntriesFast(); iMixTrack++) {
	      aodmix = (AliAODTrack*)MixTracks->UncheckedAt(iMixTrack);

	      DeltaPhi = PhiCorrection(XiPhi, aodmix->Phi());
	      DeltaRap = XiEta - aodmix->Eta(); // Pseudorapidity instead of rapidity   
	      if(TrackPID(aodmix, kTRUE) == AliPID::kPion) { hMixXi[0][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodmix->Pt()); }
	      else if(TrackPID(aodmix, kTRUE) == AliPID::kKaon) { hMixXi[1][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodmix->Pt()); }	  
	      else if(TrackPID(aodmix, kTRUE) == AliPID::kProton) { hMixXi[2][Centindex][Zvtxindex][Ptindex]->Fill(DeltaPhi, DeltaRap, aodmix->Pt()); }
	    }
	  }
	}
      }
    }
  }

  PostData(1, fListOfHistos);
}

//==================================================================================================================================================================================================================

Bool_t AliAnalysisTaskThreePartCorr::V0TriggerCuts(AliAODv0 *V0, Double_t PrimVertex[3]) {
  //AliAODVertex *TrigVertex = (AliAODVertex*)V0->GetSecondaryVtx();
  //AliAODTrack *DTrackPos = (AliAODTrack*)TrigVertex->GetDaughter(0);
  //AliAODTrack *DTrackNeg = (AliAODTrack*)TrigVertex->GetDaughter(1);
  //if(DTrackPos->Pt() < 0.15 || DTrackNeg->Pt() < 0.15 || DTrackPos->Pt() > 20 || DTrackNeg->Pt() > 20) { return kFALSE; }
  //if(TMath::Abs(DTrackPos->Eta()) > 0.8 || TMath::Abs(DTrackNeg->Eta()) > 0.8) { return kFALSE; }
  
  if(V0->Pt() < 1.0 || V0->Pt() > 12.0) { return kFALSE; }
  if(TMath::Abs(V0->Eta()) > 0.72) { return kFALSE; }
  //if(V0->DcaV0Daughters() > 1.0) { return kFALSE; }
  //if(V0->CosPointingAngle(PrimVertex) < 0.995) { return kFALSE; }

  return kTRUE;
}

Bool_t AliAnalysisTaskThreePartCorr::XiTriggerCuts(AliAODcascade *Xi) {
  Double_t Pt = TMath::Sqrt(Xi->Pt2Xi());
  Double_t Eta = TMath::ATanH((Xi->MomXiZ())/(TMath::Sqrt(Xi->Ptot2Xi())));
  if(Pt < 1.2 || Pt > 12.0) { return kFALSE; }
  if(TMath::Abs(Eta) > 0.7) { return kFALSE; }

  return kTRUE;
}

Bool_t AliAnalysisTaskThreePartCorr::AssociateCuts(AliAODTrack *Track) {
  if(Track->Pt() < 0.2 || Track->Pt() > 3.0) { return kFALSE; }
  if(TMath::Abs(Track->Eta()) > 0.8) { return kFALSE; }
  if(!(Track->TestFilterBit(768))) { return kFALSE; }

  Double_t minNSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Track, TrackPID(Track)));
  if(minNSigma > 4) { return kFALSE; }

  return kTRUE;
} 

Double_t AliAnalysisTaskThreePartCorr::PhiCorrection(Double_t TriggerPhi, Double_t AssociatePhi) {
  Double_t dPhi = TriggerPhi - AssociatePhi;
  
  if(dPhi < (-1./2)*TMath::Pi()) { dPhi = dPhi + 2*TMath::Pi(); }
  else if(dPhi > (3./2)*TMath::Pi()) { dPhi = dPhi - 2*TMath::Pi(); }
  return dPhi;
}

AliPID::EParticleType AliAnalysisTaskThreePartCorr::TrackPID(AliAODTrack *Track, Bool_t FillHist) {
  Double_t NSigma[3];
  if(FillHist) {
    AliCFParticle *CFParticle = (AliCFParticle*)Track;
    NSigma[0] = CFParticle->GetAt(0);
    NSigma[1] = CFParticle->GetAt(1);
    NSigma[2] = CFParticle->GetAt(2);
  }
  else {     
    NSigma[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Track, AliPID::kPion));
    NSigma[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Track, AliPID::kKaon));
    NSigma[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Track, AliPID::kProton));
  }
    
  if(NSigma[0] < std::min(NSigma[1], NSigma[2])) { return AliPID::kPion; }
  else if(NSigma[1] < std::min(NSigma[0], NSigma[2])) { return AliPID::kKaon; }
  else if(NSigma[2] < std::min(NSigma[0], NSigma[1])) { return AliPID::kProton; }

  return AliPID::kUnknown;
}

TObjArray *AliAnalysisTaskThreePartCorr::CloneAndReduceTrackList(TObjArray *Tracks) {
  TObjArray *CFArray = new TObjArray;
  CFArray->SetOwner(kTRUE);

  for(Int_t i = 0; i < Tracks->GetEntriesFast(); i++) {
    AliVParticle *Particle = (AliVParticle*)Tracks->UncheckedAt(i);
    AliCFParticle *CF = new AliCFParticle(Particle->Pt(), Particle->Eta(), Particle->Phi(), Particle->Charge(), 0, 3);
    CF->SetAt(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Particle, AliPID::kPion)), 0);
    CF->SetAt(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Particle, AliPID::kKaon)), 1);
    CF->SetAt(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Particle, AliPID::kProton)), 2);
    CF->SetUniqueID(Particle->GetUniqueID());
    CFArray->Add(CF);
  }

  return CFArray;
}   

//==================================================================================================================================================================================================================

void AliAnalysisTaskThreePartCorr::Terminate(Option_t *) {
  
}

//==================================================================================================================================================================================================================
