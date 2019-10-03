/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////
// Basic QA task to monitor the observables related to
// flow and balance function analysis. It is intended
// for both pp and PbPb data.
//
// Mainteiners:
//   Carlos Perez (cperez@cern.ch)
//   Alis Rodriguez (alisrm@nikhef.nl)
//////////////////////////////////////////////////////

#include "TChain.h"
#include "TList.h"
#include "TH2D.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVZERO.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliGlobalFBFqa.h"

ClassImp(AliGlobalFBFqa)

// C O N S T R U C T O R =======================================================
AliGlobalFBFqa::AliGlobalFBFqa() :
  AliAnalysisTaskSE(), fDebugger(kFALSE), fOutputList(NULL), fEvents(NULL),
  fPhiRinVZERO(NULL), fPhiEtaTPC50(NULL), fPhiEtaITSSA(NULL), 
  fPsiCenVZERO(NULL), fPsiCenTPC50(NULL), fPsiCenITSSA(NULL) {
  // default constructor
  if(fDebugger) printf("AliGlobalFBFqa:: Default Constructor\n");
}
// C O N S T R U C T O R =======================================================
AliGlobalFBFqa::AliGlobalFBFqa(const char *name) :
  AliAnalysisTaskSE(name), fDebugger(kFALSE), fOutputList(NULL), fEvents(NULL),
  fPhiRinVZERO(NULL), fPhiEtaTPC50(NULL), fPhiEtaITSSA(NULL), 
  fPsiCenVZERO(NULL), fPsiCenTPC50(NULL), fPsiCenITSSA(NULL) {
  // named constructor
  if(fDebugger) printf("AliGlobalFBFqa:: Named Constructor\n");
  DefineInput( 0,TChain::Class());
  DefineOutput(1,TList::Class());
}
// D E S T R U C T O R =========================================================
AliGlobalFBFqa::~AliGlobalFBFqa() {
  // destructor
  if(fDebugger) printf("AliGlobalFBFqa:: Destructor\n");
  if(fOutputList) delete fOutputList;
}
// U S E R   C R E A T E   O U T P U T   O B J E C T S =========================
void AliGlobalFBFqa::UserCreateOutputObjects() {
  // user create output object
  if(fDebugger) printf("AliGlobalFBFqa:: UserCreateOutputObjects\n");
  fOutputList = new TList();
  fOutputList->SetOwner();
  // number of events
  TList *tQAEvents = new TList();
  tQAEvents->SetName("Events");
  tQAEvents->SetOwner();
    fEvents = new TH2D("Events","Events;TRK;V0M", 20,0,100, 20,0,100);
    tQAEvents->Add(fEvents);
  fOutputList->Add(tQAEvents);
  TList *tQAPhaseSpace = new TList();
  tQAPhaseSpace->SetOwner();
  tQAPhaseSpace->SetName("PhaseSpace");
    fPhiRinVZERO = new TH2D("PhiRinVZERO","PhiRinVZERO;#phi (rad);Ring Number",60,0,TMath::Pi(),8,0,8);
    tQAPhaseSpace->Add(fPhiRinVZERO);
    fPhiEtaTPC50 = new TH2D("PhiEtaTPC50","PhiEtaTPC50;#phi (rad);#eta",60,0,TMath::Pi(),16,-1.6,1.6);
    tQAPhaseSpace->Add(fPhiEtaTPC50);
    fPhiEtaITSSA = new TH2D("PhiEtaITSSA","PhiEtaITSSA;#phi (rad);#eta",60,0,TMath::Pi(),20,-2,2);
    tQAPhaseSpace->Add(fPhiEtaITSSA);
  fOutputList->Add(tQAPhaseSpace);
  TList *tQAEventPlane = new TList();
  tQAEventPlane->SetOwner();
  tQAEventPlane->SetName("EventPlane");
    fPsiCenVZERO = new TH2D("PsiCenVZERO","PsiCenVZERO;#phi (rad);Centrality TRK",60,0,TMath::Pi(),20,0,100);
    tQAEventPlane->Add(fPsiCenVZERO);
    fPsiCenTPC50 = new TH2D("PsiTPCCen50","PsiCenTPC50;#phi (rad);Centrality V0M",60,0,TMath::Pi(),20,0,100);
    tQAEventPlane->Add(fPsiCenTPC50);
    fPsiCenITSSA = new TH2D("PsiCenITSSA","PsiCenITSSA;#phi (rad);Centrality V0M",60,0,TMath::Pi(),20,0,100);
    tQAEventPlane->Add(fPsiCenITSSA);
  fOutputList->Add(tQAEventPlane);
  //
  PostData(1,fOutputList);
}
// U S E R   E X E C ===========================================================
void AliGlobalFBFqa::UserExec(Option_t *) {
  // user exec
  if(fDebugger) printf("AliGlobalFBFqa:: UserExec\n");
  AliESDEvent *myESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!myESD) {
    if(fDebugger) printf("AliGlobalFBFqa:: UserExec ESDEvent not found\n");
    return;
  }
  AliESDVZERO *myVZero = myESD->GetVZEROData();
  if(!myVZero) {
    if(fDebugger) printf("AliGlobalFBFqa:: UserExec VZERO not found\n");
    return;
  }
  AliESDVertex *myVertex = (AliESDVertex*) myESD->GetPrimaryVertex();
  if(!myVertex) {
    if(fDebugger) printf("AliGlobalFBFqa:: UserExec Vertex not found\n");
    return;
  }
  if(myVertex->GetNContributors() < 2) {
    if(fDebugger) printf("AliGlobalFBFqa:: UserExec poor vertex\n");
    return;
  }
  Double_t ccTRK=-1, ccV0M=-1;
  AliCentrality *myCentrality = myESD->GetCentrality();
  if(myCentrality) {
    ccTRK = myCentrality->GetCentralityPercentileUnchecked("TRK");
    ccV0M = myCentrality->GetCentralityPercentileUnchecked("V0M");
  } else
    if(fDebugger) printf("AliGlobalFBFqa:: UserExec no centrality object\n");
  //============================================================================
  if(fDebugger) 
    printf("AliGlobalFBFqa:: UserExec cc says %.1f and %.1f\n", ccTRK,ccV0M);
  fEvents->Fill( ccTRK, ccV0M );
  Double_t dPhi, dEta;
  Double_t cosTPC=0, sinTPC=0;
  Double_t cosITS=0, sinITS=0;
  Double_t cosVZE=0, sinVZE=0;
  // Global tracks QA
  if( (ccV0M>5)&&(ccV0M<90) ) { // cutting out edges of multiplicity
    Int_t nTracks = myESD->GetNumberOfTracks();
    AliESDtrack *track;
    AliESDtrackCuts *myTPC=CutsTPC50Generic();
    AliESDtrackCuts *myITS=CutsITSSAGeneric();
    for(Int_t i=0; i!=nTracks; ++i) {
      track = (AliESDtrack*) myESD->GetTrack( i );
      dPhi = track->Phi();
      dEta = track->Eta();
      if( myTPC->IsSelected(track) ) { // TPC50
        cosTPC += TMath::Cos(2*dPhi);
        sinTPC += TMath::Sin(2*dPhi);
        fPhiEtaTPC50->Fill( dPhi, dEta );
      }
      if( myITS->IsSelected(track) ) { // ITSSA
        cosITS += TMath::Cos(2*dPhi);
        sinITS += TMath::Sin(2*dPhi);
        fPhiEtaITSSA->Fill( dPhi, dEta );
      }
    }
    double dPsiTPC = 0.5*TMath::ATan2( sinTPC, cosTPC ) + TMath::PiOver2();
    double dPsiITS = 0.5*TMath::ATan2( sinITS, cosITS ) + TMath::PiOver2();
    fPsiCenTPC50->Fill( dPsiTPC, ccV0M );
    fPsiCenITSSA->Fill( dPsiITS, ccV0M );
  }
  // VZERO QA
  if( (ccTRK>5)&&(ccTRK<90) ) { // cutting out edges of multiplicity
    for(int i=0;i!=64;++i) {
      dPhi = TMath::Pi()*(i%8)/4;
      dEta = i/8;
      cosVZE += myVZero->GetMultiplicity(i)*TMath::Cos(2*dPhi);
      sinVZE += myVZero->GetMultiplicity(i)*TMath::Sin(2*dPhi);
      fPhiRinVZERO->Fill( dPhi, dEta );
    }
    double dPsiVZE = 0.5*TMath::ATan2( sinVZE, cosVZE ) + TMath::PiOver2();
    fPsiCenVZERO->Fill( dPsiVZE, ccTRK );
  }
  PostData(1,fOutputList);
  return;
}

AliESDtrackCuts* AliGlobalFBFqa::CutsITSSAGeneric() {
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  esdTrackCuts->SetRequireITSStandAlone(kTRUE);
  esdTrackCuts->SetRequireITSPureStandAlone(kFALSE);
  esdTrackCuts->SetRequireITSRefit(kTRUE); 
  esdTrackCuts->SetMinNClustersITS(4);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  esdTrackCuts->SetMaxChi2PerClusterITS(2.5);
  // 7*(0.0033+0.0045/pt^1.3)
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0231+0.0315/pt^1.3");
  esdTrackCuts->SetRequireITSPid(kTRUE);
  esdTrackCuts->SetMaxNOfMissingITSPoints(1);
  return esdTrackCuts;
}

AliESDtrackCuts* AliGlobalFBFqa::CutsTPC50Generic() {
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  esdTrackCuts->SetMinNClustersTPC(50);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetMaxDCAToVertexZ(3.2);
  esdTrackCuts->SetMaxDCAToVertexXY(2.4);
  esdTrackCuts->SetDCAToVertex2D(kTRUE);
  return esdTrackCuts;
}

