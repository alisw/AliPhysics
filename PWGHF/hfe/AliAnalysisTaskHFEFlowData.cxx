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
//
// Flow task
// 
// Authors:
//   Raphaelle Bailhache <R.Bailhache@gsi.de>
//
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TVector2.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "AliVEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliPIDResponse.h"
#include "AliESDVZERO.h"
#include "AliESDUtils.h"
#include "AliMCParticle.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"

#include "AliFlowCandidateTrack.h"
#include "AliFlowEvent.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowVector.h"
#include "AliFlowCommonConstants.h"

#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliHFEVZEROEventPlane.h"

#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskHFEFlowData.h"


//____________________________________________________________________
AliAnalysisTaskHFEFlowData::AliAnalysisTaskHFEFlowData() :
  AliAnalysisTaskSE(),
  fListHist(0x0), 
  fAODAnalysis(kFALSE),
  fUseFlagAOD(kFALSE),
  fApplyCut(kTRUE),
  fFlags(1<<4),
  fVZEROEventPlane(kFALSE),
  fVZEROEventPlaneA(kFALSE),
  fVZEROEventPlaneC(kFALSE),
  fSubEtaGapTPC(kFALSE),
  fEtaGap(0.0),
  fDebugLevel(0),
  fHFECuts(0),
  fPID(0),
  fPIDqa(0),
  fHistEV(0),
  fEventPlane(0x0),
  fCosResabc(0x0),
  fCosRes(0x0),
  fDeltaPhiMaps(0x0),
  fCosPhiMaps(0x0)
{
  // Constructor

  
}
//______________________________________________________________________________
AliAnalysisTaskHFEFlowData:: AliAnalysisTaskHFEFlowData(const char *name) :
  AliAnalysisTaskSE(name),
  fListHist(0x0),
  fAODAnalysis(kFALSE),
  fUseFlagAOD(kFALSE),
  fApplyCut(kTRUE),
  fFlags(1<<4), 
  fVZEROEventPlane(kFALSE),
  fVZEROEventPlaneA(kFALSE),
  fVZEROEventPlaneC(kFALSE),
  fSubEtaGapTPC(kFALSE),
  fEtaGap(0.0),
  fDebugLevel(0),
  fHFECuts(0),
  fPID(0),
  fPIDqa(0),
  fHistEV(0),
  fEventPlane(0x0),
  fCosResabc(0x0),
  fCosRes(0x0),
  fDeltaPhiMaps(0x0),
  fCosPhiMaps(0x0)
{
  //
  // named ctor
  //
  
  fPID = new AliHFEpid("hfePid");
  fPIDqa = new AliHFEpidQAmanager;

  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
   
}
//____________________________________________________________
AliAnalysisTaskHFEFlowData::AliAnalysisTaskHFEFlowData(const AliAnalysisTaskHFEFlowData &ref):
  AliAnalysisTaskSE(ref),
  fListHist(0x0),
  fAODAnalysis(ref.fAODAnalysis), 
  fUseFlagAOD(ref.fUseFlagAOD),
  fApplyCut(ref.fApplyCut),
  fFlags(ref.fFlags),
  fVZEROEventPlane(ref.fVZEROEventPlane),
  fVZEROEventPlaneA(ref.fVZEROEventPlaneA),
  fVZEROEventPlaneC(ref.fVZEROEventPlaneC),
  fSubEtaGapTPC(ref.fSubEtaGapTPC),
  fEtaGap(ref.fEtaGap),
  fDebugLevel(ref.fDebugLevel),
  fHFECuts(0),
  fPID(0),
  fPIDqa(0),
  fHistEV(0),
  fEventPlane(0x0),
  fCosResabc(0x0),
  fCosRes(0x0),
  fDeltaPhiMaps(0x0),
  fCosPhiMaps(0x0)
{
  //
  // Copy Constructor
  //
  ref.Copy(*this);
}

//____________________________________________________________
AliAnalysisTaskHFEFlowData &AliAnalysisTaskHFEFlowData::operator=(const AliAnalysisTaskHFEFlowData &ref){
  //
  // Assignment operator
  //
  if(this == &ref) 
    ref.Copy(*this);
  return *this;
}

//____________________________________________________________
void AliAnalysisTaskHFEFlowData::Copy(TObject &o) const {
  // 
  // Copy into object o
  //
  AliAnalysisTaskHFEFlowData &target = dynamic_cast<AliAnalysisTaskHFEFlowData &>(o);
  target.fAODAnalysis = fAODAnalysis;
  target.fUseFlagAOD = fUseFlagAOD;
  target.fApplyCut = fApplyCut;
  target.fFlags = fFlags;
  target.fVZEROEventPlane = fVZEROEventPlane;
  target.fVZEROEventPlaneA = fVZEROEventPlaneA;
  target.fVZEROEventPlaneC = fVZEROEventPlaneC;
  target.fSubEtaGapTPC = fSubEtaGapTPC;
  target.fEtaGap = fEtaGap;
  target.fDebugLevel = fDebugLevel;
  target.fHFECuts = fHFECuts;
  target.fPID = fPID;
  target.fPIDqa = fPIDqa;
  
}
//____________________________________________________________
AliAnalysisTaskHFEFlowData::~AliAnalysisTaskHFEFlowData(){
  //
  // Destructor
  //
  if(fListHist) delete fListHist;
  if(fHFECuts) delete fHFECuts;
  if(fPID) delete fPID;
  if(fPIDqa) delete fPIDqa;
 

}
//________________________________________________________________________
void AliAnalysisTaskHFEFlowData::UserCreateOutputObjects()
{

  //********************
  // Create histograms
  //********************

  //**************
  // Cuts
  //**************

  AliInfo("AliAnalysisTaskHFEFlowData: create output objects");

  // HFE cuts

  if(!fHFECuts){
    fHFECuts = new AliHFEcuts;
    fHFECuts->CreateStandardCuts();
  }
  fHFECuts->Initialize();
  if(fAODAnalysis) fHFECuts->SetAOD();  

  AliInfo("AliAnalysisTaskHFEFlowData: HFE cuts initialize");

  // PID HFE
  //fPID->SetHasMCData(HasMCData());
  if(!fPID->GetNumberOfPIDdetectors()) fPID->AddDetector("TPC", 0);
  fPID->InitializePID();
  fPIDqa->Initialize(fPID);
  fPID->SortDetectors();

  AliInfo("AliAnalysisTaskHFEFlowData: pid and pidqa");
  
  //**************************
  // Bins for the THnSparse
  //**************************

  Int_t nBinsPt = 24;
  Double_t binLimPt[25] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2,
			   1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 3., 3.5, 4., 5.,
			   6.};
  
  Int_t nBinsEtaLess = 2;
  Double_t minEta = -0.8;
  Double_t maxEta = 0.8;
  Double_t binLimEtaLess[nBinsEtaLess+1];
  for(Int_t i=0; i<=nBinsEtaLess; i++) binLimEtaLess[i]=(Double_t)minEta + (maxEta-minEta)/nBinsEtaLess*(Double_t)i ;
 
  Int_t nBinsCos = 50;
  Double_t minCos = -1.0;
  Double_t maxCos = 1.0;
  Double_t binLimCos[nBinsCos+1];
  for(Int_t i=0; i<=nBinsCos; i++) binLimCos[i]=(Double_t)minCos + (maxCos-minCos)/nBinsCos*(Double_t)i ;
 
  Int_t nBinsC = 11;
  Double_t minC = 0.0;
  Double_t maxC = 11.0;
  Double_t binLimC[nBinsC+1];
  for(Int_t i=0; i<=nBinsC; i++) binLimC[i]=(Double_t)minC + (maxC-minC)/nBinsC*(Double_t)i ;

  Int_t nBinsCMore = 20;
  Double_t minCMore = 0.0;
  Double_t maxCMore = 20.0;
  Double_t binLimCMore[nBinsCMore+1];
  for(Int_t i=0; i<=nBinsCMore; i++) binLimCMore[i]=(Double_t)minCMore + (maxCMore-minCMore)/nBinsCMore*(Double_t)i ;

  Int_t nBinsPhi = 25;
  Double_t minPhi = 0.0;
  Double_t maxPhi = TMath::Pi();
  Double_t binLimPhi[nBinsPhi+1];
  for(Int_t i=0; i<=nBinsPhi; i++) {
    binLimPhi[i]=(Double_t)minPhi + (maxPhi-minPhi)/nBinsPhi*(Double_t)i ;
  }

  Int_t nBinsCharge = 2;
  Double_t minCharge = -1.0;
  Double_t maxCharge = 1.0;
  Double_t binLimCharge[nBinsCharge+1];
  for(Int_t i=0; i<=nBinsCharge; i++) binLimCharge[i]=(Double_t)minCharge + (maxCharge-minCharge)/nBinsCharge*(Double_t)i ;
  
  AliInfo("AliAnalysisTaskHFEFlowData: array of bins");

  //******************
  // Histograms
  //******************
    
  fListHist = new TList();
  fListHist->SetOwner();

  AliInfo("AliAnalysisTaskHFEFlowData: list created");
  
  // Histos
  fHistEV = new TH2D("fHistEV", "events", 3, 0, 3, 3, 0,3);
  
  AliInfo("AliAnalysisTaskHFEFlowData: fHistEv");

  // Event plane as function of phiep, centrality
  const Int_t nDima=5;
  Int_t nBina[nDima] = {nBinsPhi,nBinsPhi,nBinsPhi,nBinsPhi,nBinsC};
  fEventPlane = new THnSparseF("EventPlane","EventPlane",nDima,nBina);
  fEventPlane->SetBinEdges(0,binLimPhi);
  fEventPlane->SetBinEdges(1,binLimPhi);
  fEventPlane->SetBinEdges(2,binLimPhi);
  fEventPlane->SetBinEdges(3,binLimPhi);
  fEventPlane->SetBinEdges(4,binLimC);
  fEventPlane->Sumw2();

  AliInfo("AliAnalysisTaskHFEFlowData: fEventPlane");
  
  // Resolution cosres_abc centrality
  const Int_t nDimfbis=4;
  Int_t nBinfbis[nDimfbis] = {nBinsCos,nBinsCos,nBinsCos,nBinsCMore};
  fCosResabc = new THnSparseF("CosRes_abc","CosRes_abc",nDimfbis,nBinfbis);
  fCosResabc->SetBinEdges(0,binLimCos);
  fCosResabc->SetBinEdges(1,binLimCos);
  fCosResabc->SetBinEdges(2,binLimCos);
  fCosResabc->SetBinEdges(3,binLimCMore);
  fCosResabc->Sumw2();

  AliInfo("AliAnalysisTaskHFEFlowData: fCosResabc");

  // Resolution cosres centrality
  const Int_t nDimf=2;
  Int_t nBinf[nDimf] = {nBinsCos, nBinsCMore};
  fCosRes = new THnSparseF("CosRes","CosRes",nDimf,nBinf);
  fCosRes->SetBinEdges(0,binLimCos);
  fCosRes->SetBinEdges(1,binLimCMore);
  fCosRes->Sumw2();

  AliInfo("AliAnalysisTaskHFEFlowData: fCosRes");
  
  // Maps delta phi
  const Int_t nDimg=5;
  Int_t nBing[nDimg] = {nBinsPhi,nBinsC,nBinsPt, nBinsCharge,nBinsEtaLess};
  fDeltaPhiMaps = new THnSparseF("DeltaPhiMaps","DeltaPhiMaps",nDimg,nBing);
  fDeltaPhiMaps->SetBinEdges(0,binLimPhi);
  fDeltaPhiMaps->SetBinEdges(1,binLimC);
  fDeltaPhiMaps->SetBinEdges(2,binLimPt);
  fDeltaPhiMaps->SetBinEdges(3,binLimCharge);
  fDeltaPhiMaps->SetBinEdges(4,binLimEtaLess);
  fDeltaPhiMaps->Sumw2();  

  AliInfo("AliAnalysisTaskHFEFlowData: fDeltaPhiMaps");

  // Maps cos phi
  const Int_t nDimh=5;
  Int_t nBinh[nDimh] = {nBinsCos,nBinsC,nBinsPt,nBinsCharge,nBinsEtaLess};
  fCosPhiMaps = new THnSparseF("CosPhiMaps","CosPhiMaps",nDimh,nBinh);
  fCosPhiMaps->SetBinEdges(0,binLimCos);
  fCosPhiMaps->SetBinEdges(1,binLimC);
  fCosPhiMaps->SetBinEdges(2,binLimPt);
  fCosPhiMaps->SetBinEdges(3,binLimCharge);
  fCosPhiMaps->SetBinEdges(4,binLimEtaLess);
  fCosPhiMaps->Sumw2();

  AliInfo("AliAnalysisTaskHFEFlowData: fCosPhiMaps");

  //**************************
  // Add to the list
  //******************************

  //fListHist->Add(qaCutsRP);
  fListHist->Add(fPIDqa->MakeList("HFEpidQA"));
  fListHist->Add(fHistEV);
  fListHist->Add(fEventPlane);
  fListHist->Add(fCosRes);
  fListHist->Add(fCosResabc);
  fListHist->Add(fDeltaPhiMaps);
  fListHist->Add(fCosPhiMaps);

  AliInfo("AliAnalysisTaskHFEFlowData: added to the list");
  
  
  PostData(1, fListHist);

  AliInfo("AliAnalysisTaskHFEFlowData: Post Data");

}
   
//________________________________________________________________________
void AliAnalysisTaskHFEFlowData::UserExec(Option_t */*option*/)
{
  //
  // Loop over event
  //

  AliInfo("AliAnalysisTaskHFEFlowData: UserExec");
   
  Float_t cntr = 0.0;
  Double_t binct = 11.5;
  Double_t binctMore = 20.5;
  Float_t binctt = -1.0;
  
  Double_t valuensparsea[5];
  Double_t valuensparsefbis[4];
  Double_t valuensparsef[2];
  Double_t valuensparseg[5];
  Double_t valuensparseh[5];

  AliDebug(1, "Variable initialized");

  
  /////////////////
  // centrality
  /////////////////
  
  //AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  //if(!esd) return;
  AliCentrality *centrality = fInputEvent->GetCentrality();
  //printf("Got the centrality\n");
  if(!centrality) {
    PostData(1, fListHist);
    return;
  }
  cntr = centrality->GetCentralityPercentile("V0M");
  if((0.0< cntr) && (cntr<5.0)) binct = 0.5;
  if((5.0< cntr) && (cntr<10.0)) binct = 1.5;
  if((10.0< cntr) && (cntr<20.0)) binct = 2.5;
  if((20.0< cntr) && (cntr<30.0)) binct = 3.5;
  if((30.0< cntr) && (cntr<40.0)) binct = 4.5;
  if((40.0< cntr) && (cntr<50.0)) binct = 5.5;
  if((50.0< cntr) && (cntr<60.0)) binct = 6.5;
  if((60.0< cntr) && (cntr<70.0)) binct = 7.5;
  if((70.0< cntr) && (cntr<80.0)) binct = 8.5;
  if((80.0< cntr) && (cntr<90.0)) binct = 9.5;
  if((90.0< cntr) && (cntr<100.0)) binct = 10.5;
  
  if((0.< cntr) && (cntr < 20.)) binctt = 0.5;
  if((20.< cntr) && (cntr < 40.)) binctt = 1.5;
  if((40.< cntr) && (cntr < 80.)) binctt = 2.5;

  if((0.0< cntr) && (cntr<5.0)) binctMore = 0.5;
  if((5.0< cntr) && (cntr<10.0)) binctMore = 1.5;
  if((10.0< cntr) && (cntr<15.0)) binctMore = 2.5;
  if((15.0< cntr) && (cntr<20.0)) binctMore = 3.5;
  if((20.0< cntr) && (cntr<25.0)) binctMore = 4.5;
  if((25.0< cntr) && (cntr<30.0)) binctMore = 5.5;
  if((30.0< cntr) && (cntr<35.0)) binctMore = 6.5;
  if((35.0< cntr) && (cntr<40.0)) binctMore = 7.5;
  if((40.0< cntr) && (cntr<45.0)) binctMore = 8.5;
  if((45.0< cntr) && (cntr<50.0)) binctMore = 9.5;
  if((50.0< cntr) && (cntr<55.0)) binctMore = 10.5;
  if((55.0< cntr) && (cntr<60.0)) binctMore = 11.5;
  if((60.0< cntr) && (cntr<65.0)) binctMore = 12.5;
  if((65.0< cntr) && (cntr<70.0)) binctMore = 13.5;
  if((70.0< cntr) && (cntr<75.0)) binctMore = 14.5;
  if((75.0< cntr) && (cntr<80.0)) binctMore = 15.5;
  if((80.0< cntr) && (cntr<85.0)) binctMore = 16.5;
  if((85.0< cntr) && (cntr<90.0)) binctMore = 17.5;
  if((90.0< cntr) && (cntr<95.0)) binctMore = 18.5;
  if((95.0< cntr) && (cntr<100.0)) binctMore = 19.5;

  
  if(binct > 11.0) {
    PostData(1, fListHist);
    return;
  }
 
  AliDebug(1, "Centrality");

  // centrality
  valuensparsea[4] = binct;  
  valuensparsef[1] = binctMore;  
  valuensparsefbis[3] = binctMore;  
  valuensparseg[1] = binct;
  valuensparseh[1] = binct; 
  
  //////////////////////
  // run number
  //////////////////////

  Int_t runnumber = fInputEvent->GetRunNumber();
  
  if(!fPID->IsInitialized()){
    fPID->InitializePID(runnumber);
  }

  AliDebug(1, "Run number");

  //////////
  // PID
  //////////
 
  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    PostData(1, fListHist);
    return;
  }
  fPID->SetPIDResponse(pidResponse);

  AliDebug(1, "PID");

  fHistEV->Fill(binctt,0.0);
 
  AliDebug(1, "fHistEv");

  //////////////////
  // Event cut
  //////////////////
  if(!fHFECuts->CheckEventCuts("fEvRecCuts", fInputEvent)) {
    PostData(1, fListHist);
    return;
  }

  AliDebug(1, "Event cut");

  fHistEV->Fill(binctt,1.0);

  ////////////////////////////////////  
  // First method event plane
  ////////////////////////////////////

  AliEventplane* vEPa = fInputEvent->GetEventplane();
  Float_t eventPlanea = 0.0;
  Float_t eventPlaneTPC = 0.0;
  Float_t eventPlaneV0A = 0.0;
  Float_t eventPlaneV0C = 0.0;
  Float_t eventPlaneV0 = 0.0;
  TVector2 *standardQ = 0x0;
  TVector2 *qsub1a = 0x0;
  TVector2 *qsub2a = 0x0;

  // V0

  eventPlaneV0 = TVector2::Phi_0_2pi(vEPa->GetEventplane("V0", fInputEvent,2));
  if(eventPlaneV0 > TMath::Pi()) eventPlaneV0 = eventPlaneV0 - TMath::Pi();
  eventPlaneV0A = TVector2::Phi_0_2pi(vEPa->GetEventplane("V0A", fInputEvent,2));
  if(eventPlaneV0A > TMath::Pi()) eventPlaneV0A = eventPlaneV0A - TMath::Pi();
  eventPlaneV0C = TVector2::Phi_0_2pi(vEPa->GetEventplane("V0C", fInputEvent,2));
  if(eventPlaneV0C > TMath::Pi()) eventPlaneV0C = eventPlaneV0C - TMath::Pi();
  
  AliDebug(1, "V0 event plane");
  
  // TPC

  standardQ = vEPa->GetQVector(); 
  Double_t qx = -1.0;
  Double_t qy = -1.0;
  if(standardQ) {
    qx = standardQ->X();
    qy = standardQ->Y();
  }  
  TVector2 qVectorfortrack;
  qVectorfortrack.Set(qx,qy);
  eventPlaneTPC = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.; 

  AliDebug(1, "TPC event plane");

  // Choose the one used for v2

  if(fVZEROEventPlane) eventPlanea = eventPlaneV0;
  if(fVZEROEventPlaneA) eventPlanea = eventPlaneV0A;
  if(fVZEROEventPlaneC) eventPlanea = eventPlaneV0C;
  if(!fVZEROEventPlane) eventPlanea = eventPlaneTPC;

  Float_t eventPlanesub1a = -100.0;
  Float_t eventPlanesub2a = -100.0;
  Double_t diffsub1sub2a = -100.0;
  Double_t diffsubasubb = -100.0;
  Double_t diffsubasubc = -100.0;
  Double_t diffsubbsubc = -100.0;
  
  diffsubasubb = TMath::Cos(2.*(eventPlaneV0A - eventPlaneV0C));
  diffsubasubc = TMath::Cos(2.*(eventPlaneV0A - eventPlaneTPC));
  diffsubbsubc = TMath::Cos(2.*(eventPlaneV0C - eventPlaneTPC));
  
  qsub1a = vEPa->GetQsub1();
  qsub2a = vEPa->GetQsub2();
  if(qsub1a) eventPlanesub1a = TVector2::Phi_0_2pi(qsub1a->Phi())/2.;
  if(qsub2a) eventPlanesub2a = TVector2::Phi_0_2pi(qsub2a->Phi())/2.;
  if(qsub1a && qsub2a) {
    diffsub1sub2a = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
  }

  AliDebug(1, "Diff");
  
  /////////////////////////////////////////////////////////
  // Cut for event with event plane reconstructed by all
  ////////////////////////////////////////////////////////
  
  if((!standardQ) || (!qsub1a) || (!qsub2a)) {
    PostData(1, fListHist);
    return;
  }

  AliDebug(1, "Number of tracks");
  
  Int_t nbtracks = fInputEvent->GetNumberOfTracks();
  
  //////////////////////
  // Fill Histos
  //////////////////////

  fHistEV->Fill(binctt,2.0);
  
  // Fill
  valuensparsea[0] = eventPlaneV0A;
  valuensparsea[1] = eventPlaneV0C;
  valuensparsea[2] = eventPlaneTPC;
  valuensparsea[3] = eventPlaneV0;  
  fEventPlane->Fill(&valuensparsea[0]);
  
  if(!fVZEROEventPlane) {
    valuensparsef[0] = diffsub1sub2a;
    fCosRes->Fill(&valuensparsef[0]);
  }
  else {
    valuensparsefbis[0] = diffsubasubb;
    valuensparsefbis[1] = diffsubbsubc;
    valuensparsefbis[2] = diffsubasubc;
    fCosResabc->Fill(&valuensparsefbis[0]);
  }
    
  
  //////////////////////////
  // Loop over ESD track
  //////////////////////////
 
  AliDebug(1, "Loop tracks");

  for(Int_t k = 0; k < nbtracks; k++){
    
    AliVTrack *track = (AliVTrack *) fInputEvent->GetTrack(k);
    if(!track) continue;

    if(fAODAnalysis) {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
      if(!aodtrack){
	AliError("AOD track is not there");
	PostData(1, fListHist);
	return;
      }  
      //printf("Find AOD track on\n");
      if(fUseFlagAOD){
	if(aodtrack->GetFlags() != fFlags) continue;  // Only process AOD tracks where the HFE is set
      }
    }
    
    if(fApplyCut) {
      Bool_t survived = kTRUE;
      for(Int_t icut = AliHFEcuts::kStepRecKineITSTPC; icut <= AliHFEcuts::kStepHFEcutsTRD; icut++){
	if(!fHFECuts->CheckParticleCuts(icut + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)){
	  survived = kFALSE;
	  break;
	}
      }
      if(!survived) continue;
    }
    
    // Apply PID
    AliHFEpidObject hfetrack;
    if(!fAODAnalysis) hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    else hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
    hfetrack.SetRecTrack(track);
    hfetrack.SetCentrality((Int_t)binct);
    hfetrack.SetPbPb();
    if(!fPID->IsSelected(&hfetrack,0x0,"recTrackCont",fPIDqa)) {
      continue;
    }
    
    /////////////////////////////////////////////////////////
    // Subtract electron candidate from TPC event plane
    ////////////////////////////////////////////////////////
    Float_t eventplanesubtracted = 0.0;    

    if(!fVZEROEventPlane) {
      // Subtract the tracks from the event plane
      Double_t qX = standardQ->X() - vEPa->GetQContributionX(track);  //Modify the components: subtract the track you want to look at with your analysis
      Double_t qY = standardQ->Y() - vEPa->GetQContributionY(track);  //Modify the components: subtract the track you want to look at with your analysis
      TVector2 newQVectorfortrack;
      newQVectorfortrack.Set(qX,qY);
      eventplanesubtracted = TVector2::Phi_0_2pi(newQVectorfortrack.Phi())/2; 
    }
    else eventplanesubtracted = eventPlanea;

    ////////////////////////////////////////
    // Fill pt and eta for the THnSparseF
    ///////////////////////////////////////

    valuensparseg[2] = track->Pt();
    valuensparseh[2] = track->Pt();
    if(track->Charge() > 0.0) {
      valuensparseg[3] = 0.2;
      valuensparseh[3] = 0.2;
    }
    else {
      valuensparseg[3] = -0.2;
      valuensparseh[3] = -0.2;
    }
    valuensparseh[4] = track->Eta();
    valuensparseg[4] = track->Eta();

    ///////////////////////////////
    // Event plane without track
    /////////////////////////////
    Bool_t fillEventPlane = kTRUE;
    if(!fVZEROEventPlane){
      if((!qsub1a) || (!qsub2a)) fillEventPlane = kFALSE;
      if(fSubEtaGapTPC) {
	if(track->Eta() < (- fEtaGap/2.)) eventplanesubtracted = eventPlanesub1a;
	else if(track->Eta() > (fEtaGap/2.)) eventplanesubtracted = eventPlanesub2a;
	else fillEventPlane = kFALSE;
      }
    }
    
    ///////////////////////
    // Calculate deltaphi
    ///////////////////////
    Double_t phitrack = track->Phi();  
    Double_t deltaphi = TVector2::Phi_0_2pi(phitrack - eventplanesubtracted);
    if(deltaphi > TMath::Pi()) deltaphi = deltaphi - TMath::Pi();
   
    /////////////////////
    // Fill THnSparseF
    /////////////////////

    valuensparseg[0] = deltaphi;
    if(fillEventPlane) fDeltaPhiMaps->Fill(&valuensparseg[0]);
    
    //
    valuensparseh[0] = TMath::Cos(2*TVector2::Phi_mpi_pi(phitrack-eventplanesubtracted));
    if(fillEventPlane) {
      fCosPhiMaps->Fill(&valuensparseh[0]);
    }
    
  }


  AliDebug(1, "Post data");
  
  PostData(1, fListHist);


 
}
