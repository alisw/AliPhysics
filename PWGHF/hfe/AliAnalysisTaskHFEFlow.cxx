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

#include "AliAnalysisTaskSE.h"

#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliPIDResponse.h"
#include "AliESDVZERO.h"
#include "AliESDUtils.h"
#include "AliMCParticle.h"

#include "AliFlowCandidateTrack.h"
#include "AliFlowEvent.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowVector.h"
#include "AliFlowCommonConstants.h"


#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"



#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskHFEFlow.h"


//____________________________________________________________________
AliAnalysisTaskHFEFlow::AliAnalysisTaskHFEFlow() :
  AliAnalysisTaskSE(),
  fListHist(), 
  fVZEROEventPlane(kFALSE),
  fVZEROEventPlaneA(kFALSE),
  fVZEROEventPlaneC(kFALSE),
  fSubEtaGapTPC(kFALSE),
  fEtaGap(0.0),
  fNbBinsCentralityQCumulant(4),
  fNbBinsPtQCumulant(15),
  fMinPtQCumulant(0.0),
  fMaxPtQCumulant(6.0),
  fAfterBurnerOn(kFALSE),
  fNonFlowNumberOfTrackClones(0),
  fV1(0.),
  fV2(0.),
  fV3(0.),
  fV4(0.),
  fV5(0.),
  fMaxNumberOfIterations(100),
  fPrecisionPhi(0.001),
  fUseMCReactionPlane(kFALSE),
  fMCPID(kFALSE),
  fDebugLevel(0),
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fPID(0),
  fPIDqa(0),
  fflowEvent(0x0),
  fHistEV(0),
  fEventPlane(0x0),
  fEventPlaneaftersubtraction(0x0),
  fCosSin2phiep(0x0),
  fCos2phie(0x0),
  fSin2phie(0x0),
  fCos2phiep(0x0),
  fSin2phiep(0x0),
  fSin2phiephiep(0x0),
  fCosResabc(0x0),
  fProfileCosResab(0x0),
  fProfileCosResac(0x0),
  fProfileCosResbc(0x0),
  fCosRes(0x0),
  fProfileCosRes(0x0),
  fDeltaPhiMaps(0x0),
  fCosPhiMaps(0x0),
  fProfileCosPhiMaps(0x0)
{
  // Constructor

  for(Int_t k = 0; k < 10; k++) {
    fBinCentralityLess[k] = 0.0;
  }
  
}
//______________________________________________________________________________
AliAnalysisTaskHFEFlow:: AliAnalysisTaskHFEFlow(const char *name) :
  AliAnalysisTaskSE(name),
  fListHist(), 
  fVZEROEventPlane(kFALSE),
  fVZEROEventPlaneA(kFALSE),
  fVZEROEventPlaneC(kFALSE),
  fSubEtaGapTPC(kFALSE),
  fEtaGap(0.0),
  fNbBinsCentralityQCumulant(4),
  fNbBinsPtQCumulant(15),
  fMinPtQCumulant(0.0),
  fMaxPtQCumulant(6.0),
  fAfterBurnerOn(kFALSE),
  fNonFlowNumberOfTrackClones(0),
  fV1(0.),
  fV2(0.),
  fV3(0.),
  fV4(0.),
  fV5(0.),
  fMaxNumberOfIterations(100),
  fPrecisionPhi(0.001),
  fUseMCReactionPlane(kFALSE),
  fMCPID(kFALSE),
  fDebugLevel(0),
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fPID(0),
  fPIDqa(0),
  fflowEvent(0x0),
  fHistEV(0),
  fEventPlane(0x0),
  fEventPlaneaftersubtraction(0x0),
  fCosSin2phiep(0x0),
  fCos2phie(0x0),
  fSin2phie(0x0),
  fCos2phiep(0x0),
  fSin2phiep(0x0),
  fSin2phiephiep(0x0),
  fCosResabc(0x0),
  fProfileCosResab(0x0),
  fProfileCosResac(0x0),
  fProfileCosResbc(0x0),
  fCosRes(0x0),
  fProfileCosRes(0x0),
  fDeltaPhiMaps(0x0),
  fCosPhiMaps(0x0),
  fProfileCosPhiMaps(0x0)
{
  //
  // named ctor
  //
  
  for(Int_t k = 0; k < 10; k++) {
    fBinCentralityLess[k] = 0.0;
  }
  fBinCentralityLess[0] = 0.0;
  fBinCentralityLess[1] = 20.0;
  fBinCentralityLess[2] = 40.0;
  fBinCentralityLess[3] = 60.0;
  fBinCentralityLess[4] = 80.0;
  
  fPID = new AliHFEpid("hfePid");
  fPIDqa = new AliHFEpidQAmanager;

  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
  for(Int_t bincless = 0; bincless < fNbBinsCentralityQCumulant; bincless++) {
    DefineOutput(bincless+2,AliFlowEventSimple::Class()); 
  }
  
}
//________________________________________________________________________
void AliAnalysisTaskHFEFlow::UserCreateOutputObjects()
{

  //********************
  // Create histograms
  //********************

  //**************
  // Cuts
  //**************

  //---------Data selection----------
  //kMC, kGlobal, kTPCstandalone, kSPDtracklet, kPMD
  //AliFlowTrackCuts::trackParameterType rptype = AliFlowTrackCuts::kGlobal;
  //AliFlowTrackCuts::trackParameterType poitype = AliFlowTrackCuts::kGlobal;

  //---------Parameter mixing--------
  //kPure - no mixing, kTrackWithMCkine, kTrackWithMCPID, kTrackWithMCpt
  //AliFlowTrackCuts::trackParameterMix rpmix = AliFlowTrackCuts::kPure;
  //AliFlowTrackCuts::trackParameterMix poimix = AliFlowTrackCuts::kPure;

  // RP TRACK CUTS:
  fcutsRP = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();
  fcutsRP->SetName("StandartTPC");
  fcutsRP->SetEtaRange(-0.9,0.9);
  fcutsRP->SetQA(kTRUE);
  TList *qaCutsRP = fcutsRP->GetQA();
  qaCutsRP->SetName("QA_StandartTPC_RP");

  //POI TRACK CUTS:
  fcutsPOI = new AliFlowTrackCuts("dummy");
  fcutsPOI->SetParamType(AliFlowTrackCuts::kGlobal);
  fcutsPOI->SetPtRange(+1,-1); // select nothing QUICK
  fcutsPOI->SetEtaRange(+1,-1); // select nothing VZERO

  // Flow
  AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(10000);
  cc->SetMultMin(0);
  cc->SetMultMax(10000.);
  cc->SetNbinsPt(fNbBinsPtQCumulant);
  cc->SetPtMin(fMinPtQCumulant);
  cc->SetPtMax(fMaxPtQCumulant);
  cc->SetNbinsPhi(180);
  cc->SetPhiMin(0.0);
  cc->SetPhiMax(TMath::TwoPi());
  cc->SetNbinsEta(200);
  cc->SetEtaMin(-0.9);
  cc->SetEtaMax(+0.9);
  cc->SetNbinsQ(500);
  cc->SetQMin(0.0);
  cc->SetQMax(3.0);

  
  // HFE cuts

  if(!fHFECuts){
    fHFECuts = new AliHFEcuts;
    fHFECuts->CreateStandardCuts();
  }
  fHFECuts->Initialize();
  
  // PID HFE
  //fPID->SetHasMCData(HasMCData());
  if(!fPID->GetNumberOfPIDdetectors()) fPID->AddDetector("TPC", 0);
  fPID->InitializePID();
  fPIDqa->Initialize(fPID);
  fPID->SortDetectors();
  
  //**************************
  // Bins for the THnSparse
  //**************************

  Int_t nBinsPt = 25;
  Double_t minPt = 0.001;
  Double_t maxPt = 10.0;
  Double_t binLimLogPt[nBinsPt+1];
  Double_t binLimPt[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++) binLimLogPt[i]=(Double_t)TMath::Log10(minPt) + (TMath::Log10(maxPt)-TMath::Log10(minPt))/nBinsPt*(Double_t)i ;
  for(Int_t i=0; i<=nBinsPt; i++) binLimPt[i]=(Double_t)TMath::Power(10,binLimLogPt[i]);

  Int_t nBinsPtPlus = 15;
  Double_t minPtPlus = 0.0;
  Double_t maxPtPlus = 6.0;
  Double_t binLimPtPlus[nBinsPtPlus+1];
  for(Int_t i=0; i<=nBinsPtPlus; i++) binLimPtPlus[i]=(Double_t)minPtPlus + (maxPtPlus-minPtPlus)/nBinsPtPlus*(Double_t)i ;

  Int_t nBinsEta = 8;
  Double_t minEta = -0.8;
  Double_t maxEta = 0.8;
  Double_t binLimEta[nBinsEta+1];
  for(Int_t i=0; i<=nBinsEta; i++) binLimEta[i]=(Double_t)minEta + (maxEta-minEta)/nBinsEta*(Double_t)i ;
 
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
    //printf("bin phi is %f for %d\n",binLimPhi[i],i);
  }
  
  //******************
  // Histograms
  //******************
    
  fListHist = new TList();

  // Histos
  fHistEV = new TH2D("fHistEV", "events", 3, 0, 3, 3, 0,3);
  
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
  
  // Event Plane after subtraction as function of phiep, centrality, pt, eta
  const Int_t nDimb=2;
  Int_t nBinb[nDimb] = {nBinsPhi, nBinsC};
  fEventPlaneaftersubtraction = new THnSparseF("EventPlane_aftersubtraction","EventPlane_aftersubtraction",nDimb,nBinb);
  fEventPlaneaftersubtraction->SetBinEdges(0,binLimPhi);
  fEventPlaneaftersubtraction->SetBinEdges(1,binLimC);
  fEventPlaneaftersubtraction->Sumw2();

  // Monitoring of the event Plane cos(2phi) sin(2phi) centrality
  const Int_t nDimi=3;
  Int_t nBini[nDimi] = {nBinsCos, nBinsCos, nBinsCMore};
  fCosSin2phiep = new THnSparseF("CosSin2phiep","CosSin2phiep",nDimi,nBini);
  fCosSin2phiep->SetBinEdges(0,binLimCos);
  fCosSin2phiep->SetBinEdges(1,binLimCos);
  fCosSin2phiep->SetBinEdges(2,binLimCMore);
  fCosSin2phiep->Sumw2();

  // Monitoring Event plane after subtraction of the track
  const Int_t nDime=4;
  Int_t nBine[nDime] = {nBinsCos, nBinsC, nBinsPt, nBinsEta};
  fCos2phie = new THnSparseF("cos2phie","cos2phie",nDime,nBine);
  fCos2phie->SetBinEdges(2,binLimPt);
  fCos2phie->SetBinEdges(3,binLimEta);
  fCos2phie->SetBinEdges(0,binLimCos);
  fCos2phie->SetBinEdges(1,binLimC);
  fCos2phie->Sumw2();
  fSin2phie = new THnSparseF("sin2phie","sin2phie",nDime,nBine);
  fSin2phie->SetBinEdges(2,binLimPt);
  fSin2phie->SetBinEdges(3,binLimEta);
  fSin2phie->SetBinEdges(0,binLimCos);
  fSin2phie->SetBinEdges(1,binLimC);
  fSin2phie->Sumw2();
  fCos2phiep = new THnSparseF("cos2phiep","cos2phiep",nDime,nBine);
  fCos2phiep->SetBinEdges(2,binLimPt);
  fCos2phiep->SetBinEdges(3,binLimEta);
  fCos2phiep->SetBinEdges(0,binLimCos);
  fCos2phiep->SetBinEdges(1,binLimC);
  fCos2phiep->Sumw2();
  fSin2phiep = new THnSparseF("sin2phiep","sin2phiep",nDime,nBine);
  fSin2phiep->SetBinEdges(2,binLimPt);
  fSin2phiep->SetBinEdges(3,binLimEta);
  fSin2phiep->SetBinEdges(0,binLimCos);
  fSin2phiep->SetBinEdges(1,binLimC);
  fSin2phiep->Sumw2();
  fSin2phiephiep = new THnSparseF("sin2phie_phiep","sin2phie_phiep",nDime,nBine);
  fSin2phiephiep->SetBinEdges(2,binLimPt);
  fSin2phiephiep->SetBinEdges(3,binLimEta);
  fSin2phiephiep->SetBinEdges(0,binLimCos);
  fSin2phiephiep->SetBinEdges(1,binLimC);
  fSin2phiephiep->Sumw2();  

  // Resolution cosres_abc centrality
  const Int_t nDimfbis=4;
  Int_t nBinfbis[nDimfbis] = {nBinsCos,nBinsCos,nBinsCos,nBinsCMore};
  fCosResabc = new THnSparseF("CosRes_abc","CosRes_abc",nDimfbis,nBinfbis);
  fCosResabc->SetBinEdges(0,binLimCos);
  fCosResabc->SetBinEdges(1,binLimCos);
  fCosResabc->SetBinEdges(2,binLimCos);
  fCosResabc->SetBinEdges(3,binLimCMore);
  fCosResabc->Sumw2();

  // Profile cosres centrality with 3 subevents
  fProfileCosResab = new TProfile("ProfileCosRes_a_b","ProfileCosRes_a_b",nBinsCMore,binLimCMore);
  fProfileCosResab->Sumw2();
  fProfileCosResac = new TProfile("ProfileCosRes_a_c","ProfileCosRes_a_c",nBinsCMore,binLimCMore);
  fProfileCosResac->Sumw2();
  fProfileCosResbc = new TProfile("ProfileCosRes_b_c","ProfileCosRes_b_c",nBinsCMore,binLimCMore);
  fProfileCosResbc->Sumw2();
  
  // Resolution cosres centrality
  const Int_t nDimf=2;
  Int_t nBinf[nDimf] = {nBinsCos, nBinsCMore};
  fCosRes = new THnSparseF("CosRes","CosRes",nDimf,nBinf);
  fCosRes->SetBinEdges(0,binLimCos);
  fCosRes->SetBinEdges(1,binLimCMore);
  fCosRes->Sumw2();

  // Profile cosres centrality
  fProfileCosRes = new TProfile("ProfileCosRes","ProfileCosRes",nBinsCMore,binLimCMore);
  fProfileCosRes->Sumw2();
  
  // Maps delta phi
  const Int_t nDimg=3;
  Int_t nBing[nDimg] = {nBinsPhi,nBinsC,nBinsPt};
  fDeltaPhiMaps = new THnSparseF("DeltaPhiMaps","DeltaPhiMaps",nDimg,nBing);
  fDeltaPhiMaps->SetBinEdges(0,binLimPhi);
  fDeltaPhiMaps->SetBinEdges(1,binLimC);
  fDeltaPhiMaps->SetBinEdges(2,binLimPt);
  fDeltaPhiMaps->Sumw2();  

  // Maps cos phi
  const Int_t nDimh=3;
  Int_t nBinh[nDimh] = {nBinsCos,nBinsC,nBinsPt};
  fCosPhiMaps = new THnSparseF("CosPhiMaps","CosPhiMaps",nDimh,nBinh);
  fCosPhiMaps->SetBinEdges(0,binLimCos);
  fCosPhiMaps->SetBinEdges(1,binLimC);
  fCosPhiMaps->SetBinEdges(2,binLimPt);
  fCosPhiMaps->Sumw2();

  // Profile Maps cos phi
  fProfileCosPhiMaps = new TProfile2D("ProfileCosPhiMaps","ProfileCosPhiMaps",fNbBinsCentralityQCumulant,&fBinCentralityLess[0],nBinsPtPlus,binLimPtPlus);
  fProfileCosPhiMaps->Sumw2();


  //**************************
  // Add to the list
  //******************************

  fListHist->Add(qaCutsRP);
  fListHist->Add(fPIDqa->MakeList("HFEpidQA"));
  fListHist->Add(fHistEV);
  fListHist->Add(fProfileCosRes);
  fListHist->Add(fProfileCosResab);
  fListHist->Add(fProfileCosResac);
  fListHist->Add(fProfileCosResbc);
  fListHist->Add(fCosSin2phiep);
  fListHist->Add(fEventPlane);
  fListHist->Add(fEventPlaneaftersubtraction);
  fListHist->Add(fCos2phie);
  fListHist->Add(fSin2phie);
  fListHist->Add(fCos2phiep);
  fListHist->Add(fSin2phiep);
  fListHist->Add(fSin2phiephiep);
  fListHist->Add(fCosRes);
  fListHist->Add(fCosResabc);
  fListHist->Add(fDeltaPhiMaps);
  fListHist->Add(fCosPhiMaps);
  fListHist->Add(fProfileCosPhiMaps);
  

  PostData(1, fListHist);


}
   
//________________________________________________________________________
void AliAnalysisTaskHFEFlow::UserExec(Option_t */*option*/)
{
  //
  // Loop over event
  //
   
  Double_t massElectron = 0.000511;
  Double_t mcReactionPlane = 0.0;
  Bool_t   eventplanedefined = kTRUE;

  Float_t cntr = 0.0;
  Double_t binct = 11.5;
  Double_t binctMore = 20.5;
  Double_t binctLess = -0.5;
  Float_t binctt = -1.0;
  
  Double_t valuecossinephiep[3];
  Double_t valuensparsea[5];
  Double_t valuensparseabis[5];
  Double_t valuensparsee[4];
  Double_t valuensparsef[2];
  Double_t valuensparsefbis[4];
  Double_t valuensparseg[3];
  Double_t valuensparseh[3];
  Double_t valuensparsehprofile[3];

  AliMCEvent *mcEvent = MCEvent();
  AliMCParticle *mctrack = NULL;
  
  /////////////////
  // centrality
  /////////////////
  
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!esd) return;
  AliCentrality *centrality = esd->GetCentrality();
  if(!centrality) return;
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

  binctLess = cntr;
  
  //for(Int_t bincless = 0; bincless < fNbBinsCentralityQCumulant; bincless++) {
  //if((fBinCentralityLess[bincless]< cntr) && (cntr < fBinCentralityLess[bincless+1])) {
  // binctLess = bincless+0.5;
  //printf("fBinCentralityLess[bincless] %f and binctLess %f\n",fBinCentralityLess[bincless],binctLess);
  //}
  //}
    
  
  if(binct > 11.0) return;
 
  // centrality
  valuensparsea[4] = binct;  
  valuensparseabis[1] = binct;  
  valuensparsee[1] = binct;    
  valuensparsef[1] = binctMore;  
  valuensparsefbis[3] = binctMore;  
  valuensparseg[1] = binct;
  valuensparseh[1] = binct; 
  valuensparsehprofile[1] = binctLess; 
  valuecossinephiep[2] = binctMore;

  //////////////////////
  // run number
  //////////////////////

  Int_t runnumber = esd->GetRunNumber();
   
  if(!fPID->IsInitialized()){
    // Initialize PID with the given run number
    fPID->InitializePID(runnumber);
  }


  //////////
  // PID
  //////////
 
  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    //printf("No PID response set\n");
    return;
  }
  fPID->SetPIDResponse(pidResponse);

  fHistEV->Fill(binctt,0.0);
 

  //////////////////
  // Event cut
  //////////////////
  if(!fHFECuts->CheckEventCuts("fEvRecCuts", esd)) {
    PostData(1, fListHist);
    return;
  }

  fHistEV->Fill(binctt,1.0);

  ////////////////////////////////////  
  // First method event plane
  ////////////////////////////////////

  AliEventplane* esdEPa = esd->GetEventplane();
  Float_t eventPlanea = 0.0;
  Float_t eventPlaneTPC = 0.0;
  Float_t eventPlaneV0A = 0.0;
  Float_t eventPlaneV0C = 0.0;
  Float_t eventPlaneV0 = 0.0;
  TVector2 *standardQ = 0x0;
  TVector2 *qsub1a = 0x0;
  TVector2 *qsub2a = 0x0;

  // V0
 
  eventPlaneV0 = TVector2::Phi_0_2pi(esdEPa->GetEventplane("V0", esd,2));
  if(eventPlaneV0 > TMath::Pi()) eventPlaneV0 = eventPlaneV0 - TMath::Pi();
  eventPlaneV0A = TVector2::Phi_0_2pi(esdEPa->GetEventplane("V0A", esd,2));
  if(eventPlaneV0A > TMath::Pi()) eventPlaneV0A = eventPlaneV0A - TMath::Pi();
  eventPlaneV0C = TVector2::Phi_0_2pi(esdEPa->GetEventplane("V0C", esd,2));
  if(eventPlaneV0C > TMath::Pi()) eventPlaneV0C = eventPlaneV0C - TMath::Pi();
  
  // TPC

  standardQ = esdEPa->GetQVector(); 
  Double_t qx = -1.0;
  Double_t qy = -1.0;
  if(standardQ) {
    qx = standardQ->X();
    qy = standardQ->Y();
  }  
  TVector2 qVectorfortrack;
  qVectorfortrack.Set(qx,qy);
  eventPlaneTPC = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.; 

  // Choose the one used for v2

  if(fVZEROEventPlane) eventPlanea = eventPlaneV0;
  if(fVZEROEventPlaneA) eventPlanea = eventPlaneV0A;
  if(fVZEROEventPlaneC) eventPlanea = eventPlaneV0C;
  if(!fVZEROEventPlane) eventPlanea = eventPlaneTPC;

  valuecossinephiep[0] = TMath::Cos(2*eventPlanea);
  valuecossinephiep[1] = TMath::Sin(2*eventPlanea);

  Float_t eventPlanesub1a = -100.0;
  Float_t eventPlanesub2a = -100.0;
  Double_t diffsub1sub2a = -100.0;
  Double_t diffsubasubb = -100.0;
  Double_t diffsubasubc = -100.0;
  Double_t diffsubbsubc = -100.0;

  //if(fVZEROEventPlane) {
  diffsubasubb = TMath::Cos(2.*(eventPlaneV0A - eventPlaneV0C));
  diffsubasubc = TMath::Cos(2.*(eventPlaneV0A - eventPlaneTPC));
  diffsubbsubc = TMath::Cos(2.*(eventPlaneV0C - eventPlaneTPC));
  //}
  //else {
  qsub1a = esdEPa->GetQsub1();
  qsub2a = esdEPa->GetQsub2();
  if(qsub1a) eventPlanesub1a = TVector2::Phi_0_2pi(qsub1a->Phi())/2.;
  if(qsub2a) eventPlanesub2a = TVector2::Phi_0_2pi(qsub2a->Phi())/2.;
  if(qsub1a && qsub2a) diffsub1sub2a = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
  //}
  

  /////////////////////////////////////////////////////////
  // Cut for event with event plane reconstructed by all
  ////////////////////////////////////////////////////////
  
  //if(!fVZEROEventPlane) {
  // if(!standardQ) {
  //  eventplanedefined = kFALSE;
      //PostData(1, fListHist);
      //return;
  //}
  //}

  if((!standardQ) || (!qsub1a) || (!qsub2a)) return;

  ///////////////////////
  // AliFlowEvent  
  //////////////////////

  Int_t nbtracks = esd->GetNumberOfTracks();
  //printf("Number of tracks %d\n",nbtracks);

  fcutsRP->SetEvent( InputEvent(), MCEvent());
  fcutsPOI->SetEvent( InputEvent(), MCEvent());
  if( fflowEvent ){ 
    fflowEvent->~AliFlowEvent();
    new(fflowEvent) AliFlowEvent(fcutsRP,fcutsPOI);
  }else fflowEvent = new AliFlowEvent(fcutsRP,fcutsPOI);
  if(mcEvent && mcEvent->GenEventHeader()) {
    fflowEvent->SetMCReactionPlaneAngle(mcEvent);
    //if reaction plane not set from elsewhere randomize it before adding flow
    //if (!fflowEvent->IsSetMCReactionPlaneAngle()) fflowEvent->SetMCReactionPlaneAngle(gRandom->Uniform(0.0,TMath::TwoPi()));
    mcReactionPlane = TVector2::Phi_0_2pi(fflowEvent->GetMCReactionPlaneAngle());
    if(mcReactionPlane > TMath::Pi()) mcReactionPlane = mcReactionPlane - TMath::Pi();
    //printf("MC reaction plane %f\n",mcReactionPlane);
  }
  fflowEvent->SetReferenceMultiplicity( nbtracks );
  fflowEvent->DefineDeadZone(0,0,0,0);
  //fflowEvent.TagSubeventsInEta(-0.8,-0.1,0.1,0.8);

  ////////////////
  // MC
  ///////////////
  if(fUseMCReactionPlane) {
    eventPlanea = mcReactionPlane;
    diffsub1sub2a = 0.0;
  }

  
  //////////////////////
  // Fill Histos
  //////////////////////

  if(eventplanedefined) {
    
    fHistEV->Fill(binctt,2.0);
    
    // Fill
    valuensparsea[0] = eventPlaneV0A;
    valuensparsea[1] = eventPlaneV0C;
    valuensparsea[2] = eventPlaneTPC;
    valuensparsea[3] = eventPlaneV0;  
    fEventPlane->Fill(&valuensparsea[0]);

    // Fill
    fCosSin2phiep->Fill(&valuecossinephiep[0]);
    
    if(!fVZEROEventPlane) {
      valuensparsef[0] = diffsub1sub2a;
      fCosRes->Fill(&valuensparsef[0]);
      if(fDebugLevel > 0) {
	fProfileCosRes->Fill(valuensparsef[1],valuensparsef[0]);
      }
    }
    else {
      valuensparsefbis[0] = diffsubasubb;
      valuensparsefbis[1] = diffsubbsubc;
      valuensparsefbis[2] = diffsubasubc;
      fCosResabc->Fill(&valuensparsefbis[0]);
      if(fDebugLevel > 0) {
	fProfileCosResab->Fill(valuensparsefbis[3],valuensparsefbis[0]);
	fProfileCosResac->Fill(valuensparsefbis[3],valuensparsefbis[1]);
	fProfileCosResbc->Fill(valuensparsefbis[3],valuensparsefbis[2]);
      }
    }
    
  }
  
  //////////////////////////
  // Loop over ESD track
  //////////////////////////
 

  for(Int_t k = 0; k < nbtracks; k++){
    
    AliESDtrack *track = esd->GetTrack(k);
    if(!track) continue;

    Bool_t survived = kTRUE;
    for(Int_t icut = AliHFEcuts::kStepRecKineITSTPC; icut <= AliHFEcuts::kStepHFEcutsTRD; icut++){
      if(!fHFECuts->CheckParticleCuts(icut + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)){
	survived = kFALSE;
	break;
      }
      //printf("Pass the cut %d\n",icut);
    }
    
    if(!survived) continue;

    // Apply PID for Data
    if(!fMCPID) {
      AliHFEpidObject hfetrack;
      hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
      hfetrack.SetRecTrack(track);
      hfetrack.SetCentrality((Int_t)binct);
      //printf("centrality %f and %d\n",binct,hfetrack.GetCentrality());
      hfetrack.SetPbPb();
      if(!fPID->IsSelected(&hfetrack,0x0,"",fPIDqa)) {
	continue;
      }
    }
    else {
      if(!mcEvent) continue;
      if(!(mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel()))))) continue;
      //printf("PdgCode %d\n",TMath::Abs(mctrack->Particle()->GetPdgCode()));
      if(TMath::Abs(mctrack->Particle()->GetPdgCode())!=11) continue;
    }


    /////////////////////////////////////////////////////////////////////////////
    // Add candidate to AliFlowEvent for POI and subtract from RP if needed
    ////////////////////////////////////////////////////////////////////////////
    Int_t idtrack = static_cast<AliVTrack*>(track)->GetID();
    Bool_t found = kFALSE;
    Int_t numberoffound = 0;
    //printf("A: Number of tracks %d\n",fflowEvent->NumberOfTracks());
    for(Int_t iRPs=0; iRPs< fflowEvent->NumberOfTracks(); iRPs++) {
      AliFlowTrack *iRP = (AliFlowTrack*) (fflowEvent->GetTrack(iRPs));
      //if(!iRP->InRPSelection()) continue;
      if( TMath::Abs(idtrack) == TMath::Abs(iRP->GetID()) ) {
	iRP->SetForPOISelection(kTRUE);
	found = kTRUE;
	numberoffound ++;
      }
    }
    //printf("Found %d mal\n",numberoffound);
    if(!found) {
      AliFlowCandidateTrack *sTrack = (AliFlowCandidateTrack*) MakeTrack(massElectron,track->Pt(),track->Phi(), track->Eta());
      sTrack->SetID(idtrack);
      fflowEvent->AddTrack(sTrack);
      //printf("Add the track\n");
    }
    //printf("B: Number of tracks %d\n",fflowEvent->NumberOfTracks());
    
    
    /////////////////////////////////////////////////////////
    // Subtract electron candidate from TPC event plane
    ////////////////////////////////////////////////////////
    Float_t eventplanesubtracted = 0.0;    

    //if(eventplanedefined && (!fVZEROEventPlane)) {
    if(!fVZEROEventPlane) {
      // Subtract the tracks from the event plane
      Double_t qX = standardQ->X() - esdEPa->GetQContributionX(track);  //Modify the components: subtract the track you want to look at with your analysis
      Double_t qY = standardQ->Y() - esdEPa->GetQContributionY(track);  //Modify the components: subtract the track you want to look at with your analysis
      TVector2 newQVectorfortrack;
      newQVectorfortrack.Set(qX,qY);
      eventplanesubtracted = TVector2::Phi_0_2pi(newQVectorfortrack.Phi())/2; 
    }
    else eventplanesubtracted = eventPlanea;

    ////////////////////////////////////////
    // Fill pt and eta for the THnSparseF
    ///////////////////////////////////////

    valuensparsee[2] = track->Pt();
    valuensparsee[3] = track->Eta();    
    valuensparseg[2] = track->Pt();
    valuensparseh[2] = track->Pt();
    valuensparsehprofile[2] = track->Pt();

    Bool_t fillEventPlane = kTRUE;
    if(!fVZEROEventPlane){
      //if((!qsub1a) || (!qsub2a) || (!eventplanedefined)) fillEventPlane = kFALSE;
      if((!qsub1a) || (!qsub2a)) fillEventPlane = kFALSE;
      if(fSubEtaGapTPC) {
	if(track->Eta() < (- fEtaGap/2.)) eventplanesubtracted = eventPlanesub1a;
	else if(track->Eta() > (fEtaGap/2.)) eventplanesubtracted = eventPlanesub2a;
	else fillEventPlane = kFALSE;
      }
    }
    
    
    ///////////////
    // MC
    //////////////
    if(fUseMCReactionPlane) {
      eventplanesubtracted = mcReactionPlane;
      fillEventPlane = kTRUE;
    }
    
    //////////////////////////////////////////////////////////////////////////////
    ///////////////////////////AFTERBURNER
    Double_t phitrack = track->Phi();    
    if (fAfterBurnerOn)
      {
	phitrack = GetPhiAfterAddV2(track->Phi(),mcReactionPlane);
      }
    //////////////////////////////////////////////////////////////////////////////


    ///////////////////////
    // Calculate deltaphi
    ///////////////////////
    
    // Suppose phi track is between 0.0 and phi
    Double_t deltaphi = TVector2::Phi_0_2pi(phitrack - eventplanesubtracted);
    if(deltaphi > TMath::Pi()) deltaphi = deltaphi - TMath::Pi();
   
    /////////////////////
    // Fill THnSparseF
    /////////////////////

    //
    valuensparseabis[0] = eventplanesubtracted;
    if(fillEventPlane) fEventPlaneaftersubtraction->Fill(&valuensparseabis[0]);
    

    if(fDebugLevel > 1) 
      {
	//
	valuensparsee[0] = TMath::Cos(2*phitrack);
	fCos2phie->Fill(&valuensparsee[0]);
	valuensparsee[0] = TMath::Sin(2*phitrack);
	fSin2phie->Fill(&valuensparsee[0]);
	//
	valuensparsee[0] = TMath::Cos(2*eventplanesubtracted);
	if(fillEventPlane) fCos2phiep->Fill(&valuensparsee[0]);
	valuensparsee[0] = TMath::Sin(2*eventplanesubtracted);
	if(fillEventPlane) fSin2phiep->Fill(&valuensparsee[0]);
	valuensparsee[0] = TMath::Sin(2*TVector2::Phi_mpi_pi(phitrack-eventplanesubtracted));
	if(fillEventPlane) fSin2phiephiep->Fill(&valuensparsee[0]);
        //
      }

    // 
    valuensparseg[0] = deltaphi;
    if(fillEventPlane) fDeltaPhiMaps->Fill(&valuensparseg[0]);
    
    //
    valuensparseh[0] = TMath::Cos(2*TVector2::Phi_mpi_pi(phitrack-eventplanesubtracted));
    if(fillEventPlane) {
      fCosPhiMaps->Fill(&valuensparseh[0]);
      if(fDebugLevel > 0) {
	fProfileCosPhiMaps->Fill(valuensparsehprofile[1],valuensparsehprofile[2],valuensparseh[0]);
      }
    }
    
    
    
  }

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////AFTERBURNER
  if (fAfterBurnerOn)
    {
      fflowEvent->AddFlow(fV1,fV2,fV3,fV4,fV5);     //add flow
      fflowEvent->CloneTracks(fNonFlowNumberOfTrackClones); //add nonflow by cloning tracks
    }
  //////////////////////////////////////////////////////////////////////////////



  for(Int_t bincless = 0; bincless < fNbBinsCentralityQCumulant; bincless++) {
    if((fBinCentralityLess[bincless]< cntr) && (cntr < fBinCentralityLess[bincless+1])) PostData(bincless+2,fflowEvent);
  }
  
  PostData(1, fListHist);


 
}
//______________________________________________________________________________
AliFlowCandidateTrack *AliAnalysisTaskHFEFlow::MakeTrack( Double_t mass, 
                          Double_t pt, Double_t phi, Double_t eta) {
  //
  //  Make Track (Not needed actually)
  //

  AliFlowCandidateTrack *sTrack = new AliFlowCandidateTrack();
  sTrack->SetMass(mass);
  sTrack->SetPt(pt);
  sTrack->SetPhi(phi);
  sTrack->SetEta(eta);
  sTrack->SetForPOISelection(kTRUE);
  sTrack->SetForRPSelection(kFALSE);
  return sTrack;
}
//_________________________________________________________________________________ 
Double_t AliAnalysisTaskHFEFlow::GetPhiAfterAddV2(Double_t phi,Double_t reactionPlaneAngle) const
{
  //
  // Adds v2, uses Newton-Raphson iteration
  //
  Double_t phiend=phi;
  Double_t phi0=phi;
  Double_t f=0.;
  Double_t fp=0.;
  Double_t phiprev=0.;

  for (Int_t i=0; i<fMaxNumberOfIterations; i++)
  {
    phiprev=phiend; //store last value for comparison
    f =  phiend-phi0+fV2*TMath::Sin(2.*(phiend-reactionPlaneAngle));
    fp = 1.0+2.0*fV2*TMath::Cos(2.*(phiend-reactionPlaneAngle)); //first derivative
    phiend -= f/fp;
    if (TMath::AreEqualAbs(phiprev,phiend,fPrecisionPhi)) break;
  }
  return phiend;
}
