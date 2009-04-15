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

/* $Id: $ */
#include "AliAnaESDSpectraQA.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TList.h"
#include "TChain.h"
#include "TDirectory.h"

#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"


const Int_t AliAnaESDSpectraQA::fgkNPtBins=38;
const Float_t AliAnaESDSpectraQA::fgkPtMin=2;
const Float_t AliAnaESDSpectraQA::fgkPtMax=40;
const Int_t AliAnaESDSpectraQA::fgkNPhiBins=18;

ClassImp(AliAnaESDSpectraQA)

AliAnaESDSpectraQA::AliAnaESDSpectraQA(): AliAnalysisTask("AliAnaESDSpectraQA", ""), 
  fESD(0), 
  fTrackCuts(0), 
  fNEvent(0), // just to avoid warnings, inititialized in InitPointers too
  fPtAll(0),  //
  fPtSel(0),  //
  fHistList(0)
{
  InitHistPointers();
}

AliAnaESDSpectraQA::AliAnaESDSpectraQA(const char *name): 
  AliAnalysisTask(name, ""), 
  fESD(0),
  fTrackCuts(new AliESDtrackCuts),
  fNEvent(0), // just to avoid warnings, inititialized in InitPointers too
  fPtAll(0),  //
  fPtSel(0),  // 
  fHistList(0) {
  // Input slot #0 works with a TChain ESD
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList
  DefineOutput(0, TList::Class());
  InitHistPointers();
  //fTrackCuts = new AliESDtrackCuts;
  fTrackCuts->SetAcceptKingDaughters(kFALSE);
  fTrackCuts->SetRequireTPCRefit(kTRUE);
  fTrackCuts->SetEtaRange(-1,1);
  // Add chisq criterium to reject 'bad' tracks that might not make it as prim
}

//________________________________________________________________________
void AliAnaESDSpectraQA::ConnectInputData(Option_t *) 
{
  // Connect ESD here
  // Called once
  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  
  if( handler && handler->InheritsFrom("AliESDInputHandler") ) {
     fESD  =  ((AliESDInputHandler*)handler)->GetEvent();
  } 
  else {
    AliFatal("I can't get any ESD Event Handler");
  }
}

void AliAnaESDSpectraQA::InitHistPointers() {
  fNEvent = 0;
  fPtAll = 0;
  fPtSel = 0;
  for (Int_t i=0; i< 4; i++) {
    fHists[i].PhiPtNPointTPC = 0;
    fHists[i].PhiPtNPointITS = 0;
    fHists[i].PhiPtChisqC = 0;
    fHists[i].PhiPtChisqTPC = 0;
    fHists[i].PhiPtDCAR = 0;
    fHists[i].PhiPtDCAZ = 0;
  }
}

void AliAnaESDSpectraQA::CreateOutputObjects() {
  fHistList = new TList();
  //fDirectory = new TDirectory("trig_hists","Trigger histos");
  //fDirectory->cd();
  TString labels[4];
  labels[kNegA]="NegA";
  labels[kPosA]="PosA";
  labels[kNegC]="NegC";
  labels[kPosC]="PosC";

  static const Float_t kMinPhi = 0;
  static const Float_t kMaxPhi = 2*TMath::Pi();
  fNEvent = new TH1F("NEvent","NEvent",1,-0.5,0.5);
  fHistList->Add(fNEvent);
  fPtAll = new TH1F("PtAll","PtAll",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtAll);
  fPtSel = new TH1F("PtSel","PtSel",fgkNPtBins, fgkPtMin, fgkPtMax);
  fHistList->Add(fPtSel);
  for (Int_t iSide = 0; iSide < 4; iSide++) {
    TString hname="PhiPtNpointTPC";
    hname += labels[iSide];
    fHists[iSide].PhiPtNPointTPC = new TH3F(hname,hname+";#phi;p_{T} (GeV);N_{point,TPC}",fgkNPhiBins,kMinPhi,kMaxPhi,fgkNPtBins,fgkPtMin,fgkPtMax,160,0.5,160.5);
    fHistList->Add(fHists[iSide].PhiPtNPointTPC);

    hname="PhiPtNpointITS";
    hname += labels[iSide];
    fHists[iSide].PhiPtNPointITS = new TH3F(hname,hname+";#phi;p_{T} (GeV);N_{point,ITS}",fgkNPhiBins,kMinPhi,kMaxPhi,fgkNPtBins,fgkPtMin,fgkPtMax,9,-0.5,8.5);
    fHistList->Add(fHists[iSide].PhiPtNPointITS);

    hname="PhiPtChisqC";
    hname += labels[iSide];
    fHists[iSide].PhiPtChisqC = new TH3F(hname,hname+";#phi;p_{T} (GeV);#Chi^{2}/NDF",fgkNPhiBins,kMinPhi,kMaxPhi,fgkNPtBins,fgkPtMin,fgkPtMax,160,0,80);
    fHistList->Add(fHists[iSide].PhiPtChisqC);

    hname="PhiPtChisqTPC";
    hname += labels[iSide];
    fHists[iSide].PhiPtChisqTPC = new TH3F(hname,hname+";#phi;p_{T} (GeV);#Chi^{2}/NDF",fgkNPhiBins,kMinPhi,kMaxPhi,fgkNPtBins,fgkPtMin,fgkPtMax,50,0,5);
    fHistList->Add(fHists[iSide].PhiPtChisqTPC);

    hname="PhiPtDCAR";
    hname += labels[iSide];
    fHists[iSide].PhiPtDCAR = new TH3F(hname,hname+";#phi;p_{T} (GeV);DCAR",fgkNPhiBins,kMinPhi,kMaxPhi,fgkNPtBins,fgkPtMin,fgkPtMax,200,-1,1);
    fHistList->Add(fHists[iSide].PhiPtDCAR);

    hname="PhiPtDCAZ";
    hname += labels[iSide];
    fHists[iSide].PhiPtDCAZ = new TH3F(hname,hname+";#phi;p_{T} (GeV);DCAZ",fgkNPhiBins,kMinPhi,kMaxPhi,fgkNPtBins,fgkPtMin,fgkPtMax,200,-2,2);
    fHistList->Add(fHists[iSide].PhiPtDCAZ);

    hname="PhiPtSigmaToVertex";
    hname += labels[iSide];
    fHists[iSide].PhiPtSigmaToVertex = new TH3F(hname,hname+";#phi;p_{T} (GeV);n#sigma to vtx",fgkNPhiBins,kMinPhi,kMaxPhi,fgkNPtBins,fgkPtMin,fgkPtMax,50,0,5);
    fHistList->Add(fHists[iSide].PhiPtSigmaToVertex);
  }
}

void AliAnaESDSpectraQA::Exec(Option_t *) {  
  // Main loop
  // Called for each event

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  const AliESDVertex *vtx = fESD->GetPrimaryVertex();

  // Need vertex cut
  if (vtx->GetNContributors() < 2)
    return;

  printf("Vertex title %s, status %d, nCont %d\n",vtx->GetTitle(), vtx->GetStatus(), vtx->GetNContributors());
  // Need to keep track of evts without vertex

  Int_t nTracks = fESD->GetNumberOfTracks();
  AliDebug(2,Form("nTracks %d", nTracks));
  printf("nTracks %d\n", nTracks);
  static Int_t fMult = 0;
  fMult = 0;   // Need extra init bc of static
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {
    AliESDtrack *track = fESD->GetTrack(iTrack);
    hists *curTypeHists = 0;

    if (fTrackCuts->AcceptTrack(track)) {
      fMult++;

      Float_t dca2D, dcaZ;
      track->GetImpactParameters(dca2D,dcaZ);
      fPtAll->Fill(track->Pt());

      if (track->Eta() > 0) { // A side (crude for now)
	if (track->Charge() > 0) 
	  curTypeHists = &fHists[kPosA];
	else
	  curTypeHists = &fHists[kNegA];	
      }
      else { // C side
	if (track->Charge() > 0) 
	  curTypeHists = &fHists[kPosC];
	else
	  curTypeHists = &fHists[kNegC];	
      }

      Float_t phi = track->Phi();
      Float_t pt = track->Pt();

      UChar_t itsMap = track->GetITSClusterMap();
      Int_t nPointITS = 0;
      for (Int_t i=0; i < 6; i++) {
	if (itsMap & (1 << i))
	  nPointITS ++;
      }

      Float_t sigToVertex = fTrackCuts->GetSigmaToVertex(track);
      Float_t chisqC = 1000;
      if (track->GetConstrainedParam())
	chisqC = track->GetConstrainedChi2();

      curTypeHists->PhiPtNPointTPC->Fill(phi,pt,track->GetTPCNcls());
      curTypeHists->PhiPtNPointITS->Fill(phi,pt,nPointITS);
      curTypeHists->PhiPtChisqC->Fill(phi,pt,chisqC);
      if(track->GetTPCNclsF()>5){
	curTypeHists->PhiPtChisqTPC->Fill(phi,pt,track->GetTPCchi2()/(track->GetTPCNclsF()-5));
      }      
      curTypeHists->PhiPtDCAR->Fill(phi,pt,dca2D);
      curTypeHists->PhiPtDCAZ->Fill(phi,pt,dcaZ);
      curTypeHists->PhiPtSigmaToVertex->Fill(phi,pt,sigToVertex);
    }
  }
  
  // Post output data
  PostData(0, fHistList);
}
