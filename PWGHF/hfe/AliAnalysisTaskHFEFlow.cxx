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
  fSubEtaGapTPC(kFALSE),
  fEtaGap(0.0),
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
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fPID(0),
  fPIDqa(0),
  fflowEvent(0x0),
  fHistEV(0),
  fEventPlane(0x0),
  fEventPlaneaftersubtraction(0x0),
  fCos2phie(0x0),
  fSin2phie(0x0),
  fCos2phiep(0x0),
  fSin2phiep(0x0),
  fSin2phiephiep(0x0),
  fCosRes(0x0),
  fDeltaPhiMaps(0x0),
  fCosPhiMaps(0x0)
{
  // Constructor

  
  
}
//______________________________________________________________________________
AliAnalysisTaskHFEFlow:: AliAnalysisTaskHFEFlow(const char *name) :
  AliAnalysisTaskSE(name),
  fListHist(), 
  fSubEtaGapTPC(kFALSE),
  fEtaGap(0.0),
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
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fPID(0),
  fPIDqa(0),
  fflowEvent(0x0),
  fHistEV(0),
  fEventPlane(0x0),
  fEventPlaneaftersubtraction(0x0),
  fCos2phie(0x0),
  fSin2phie(0x0),
  fCos2phiep(0x0),
  fSin2phiep(0x0),
  fSin2phiephiep(0x0),
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
  DefineOutput(2,AliFlowEventSimple::Class()); // first band 0-20%
  DefineOutput(3,AliFlowEventSimple::Class()); // first band 20-40%
  DefineOutput(4,AliFlowEventSimple::Class()); // first band 40-60%
  DefineOutput(5,AliFlowEventSimple::Class()); // first band 60-80%
  
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
 
  Int_t nBinsPhi = 50;
  Double_t minPhi = 0.0;
  Double_t maxPhi = TMath::Pi();
  Double_t binLimPhi[nBinsPhi+1];
  for(Int_t i=0; i<=nBinsPhi; i++) {
    binLimPhi[i]=(Double_t)minPhi + (maxPhi-minPhi)/nBinsPhi*(Double_t)i ;
    //printf("bin phi is %f for %d\n",binLimPhi[i],i);
  }
 
  Int_t nBinsQ = 200;
  Double_t minQ = -300.0;
  Double_t maxQ = 300.0;
  Double_t binLimQ[nBinsQ+1];
  for(Int_t i=0; i<=nBinsQ; i++) binLimQ[i]=(Double_t)minQ + (maxQ-minQ)/nBinsQ*(Double_t)i ;

  Int_t nBinsM = 250;
  Double_t minM = 0.0;
  Double_t maxM = 800.0;
  Double_t binLimM[nBinsM+1];
  for(Int_t i=0; i<=nBinsM; i++) binLimM[i]=(Double_t)minM + (maxM-minM)/nBinsM*(Double_t)i ;

  Int_t nBinsCELL = 64;
  Double_t minCELL = 0.0;
  Double_t maxCELL = 64.0;
  Double_t binLimCELL[nBinsCELL+1];
  for(Int_t i=0; i<=nBinsCELL; i++) binLimCELL[i]=(Double_t)minCELL + (maxCELL-minCELL)/nBinsCELL*(Double_t)i ;
  

  //******************
  // Histograms
  //******************
    
  fListHist = new TList();

  // Histos
  fHistEV = new TH2D("fHistEV", "events", 3, 0, 3, 3, 0,3);
  
  // Event plane as function of phiep, centrality
  const Int_t nDima=2;
  Int_t nBina[nDima] = {nBinsPhi, nBinsC};
  fEventPlane = new THnSparseF("EventPlane","EventPlane",nDima,nBina);
  fEventPlane->SetBinEdges(0,binLimPhi);
  fEventPlane->SetBinEdges(1,binLimC);
  
  // Event Plane after subtraction as function of phiep, centrality, pt, eta
  const Int_t nDimb=2;
  Int_t nBinb[nDimb] = {nBinsPhi, nBinsC};
  fEventPlaneaftersubtraction = new THnSparseF("EventPlane_aftersubtraction","EventPlane_aftersubtraction",nDimb,nBinb);
  fEventPlaneaftersubtraction->SetBinEdges(0,binLimPhi);
  fEventPlaneaftersubtraction->SetBinEdges(1,binLimC);

  // Monitoring Event plane after subtraction of the track
  const Int_t nDime=4;
  Int_t nBine[nDime] = {nBinsCos, nBinsC, nBinsPt, nBinsEta};
  fCos2phie = new THnSparseF("cos2phie","cos2phie",nDime,nBine);
  fCos2phie->SetBinEdges(2,binLimPt);
  fCos2phie->SetBinEdges(3,binLimEta);
  fCos2phie->SetBinEdges(0,binLimCos);
  fCos2phie->SetBinEdges(1,binLimC);
  fSin2phie = new THnSparseF("sin2phie","sin2phie",nDime,nBine);
  fSin2phie->SetBinEdges(2,binLimPt);
  fSin2phie->SetBinEdges(3,binLimEta);
  fSin2phie->SetBinEdges(0,binLimCos);
  fSin2phie->SetBinEdges(1,binLimC);
  fCos2phiep = new THnSparseF("cos2phiep","cos2phiep",nDime,nBine);
  fCos2phiep->SetBinEdges(2,binLimPt);
  fCos2phiep->SetBinEdges(3,binLimEta);
  fCos2phiep->SetBinEdges(0,binLimCos);
  fCos2phiep->SetBinEdges(1,binLimC);
  fSin2phiep = new THnSparseF("sin2phiep","sin2phiep",nDime,nBine);
  fSin2phiep->SetBinEdges(2,binLimPt);
  fSin2phiep->SetBinEdges(3,binLimEta);
  fSin2phiep->SetBinEdges(0,binLimCos);
  fSin2phiep->SetBinEdges(1,binLimC);
  fSin2phiephiep = new THnSparseF("sin2phie_phiep","sin2phie_phiep",nDime,nBine);
  fSin2phiephiep->SetBinEdges(2,binLimPt);
  fSin2phiephiep->SetBinEdges(3,binLimEta);
  fSin2phiephiep->SetBinEdges(0,binLimCos);
  fSin2phiephiep->SetBinEdges(1,binLimC);
  
  // Resolution cosres centrality
  const Int_t nDimf=2;
  Int_t nBinf[nDimf] = {nBinsCos, nBinsC};
  fCosRes = new THnSparseF("CosRes","CosRes",nDimf,nBinf);
  fCosRes->SetBinEdges(0,binLimCos);
  fCosRes->SetBinEdges(1,binLimC);
  
  // Maps delta phi
  const Int_t nDimg=3;
  Int_t nBing[nDimg] = {nBinsPhi,nBinsC,nBinsPt};
  fDeltaPhiMaps = new THnSparseF("DeltaPhiMaps","DeltaPhiMaps",nDimg,nBing);
  fDeltaPhiMaps->SetBinEdges(0,binLimPhi);
  fDeltaPhiMaps->SetBinEdges(1,binLimC);
  fDeltaPhiMaps->SetBinEdges(2,binLimPt);
  
  // Maps cos phi
  const Int_t nDimh=3;
  Int_t nBinh[nDimh] = {nBinsCos,nBinsC,nBinsPt};
  fCosPhiMaps = new THnSparseF("CosPhiMaps","CosPhiMaps",nDimh,nBinh);
  fCosPhiMaps->SetBinEdges(0,binLimCos);
  fCosPhiMaps->SetBinEdges(1,binLimC);
  fCosPhiMaps->SetBinEdges(2,binLimPt);
  
  //**************************
  // Add to the list
  //******************************

  fListHist->Add(fPIDqa->MakeList("HFEpidQA"));
	
  fListHist->Add(fHistEV);

  fListHist->Add(fEventPlane);
  fListHist->Add(fEventPlaneaftersubtraction);
  fListHist->Add(fCos2phie);
  fListHist->Add(fSin2phie);
  fListHist->Add(fCos2phiep);
  fListHist->Add(fSin2phiep);
  fListHist->Add(fSin2phiephiep);
  fListHist->Add(fCosRes);
  fListHist->Add(fDeltaPhiMaps);
  fListHist->Add(fCosPhiMaps);
  

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

  Float_t cntr = 0.0;
  Double_t binct = 11.5;
  Float_t binctt = -1.0;
  
  Double_t valuensparsea[2];
  Double_t valuensparsee[4];
  Double_t valuensparsef[2];
  Double_t valuensparseg[3];
  Double_t valuensparseh[3];

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

  if(binct > 11.0) return;
 
  // centrality
  valuensparsea[1] = binct;  
  valuensparsee[1] = binct;    
  valuensparsef[1] = binct;  
  valuensparseg[1] = binct;
  valuensparseh[1] = binct; 
  

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
  // First method TPC event plane
  ////////////////////////////////////

  AliEventplane* esdEPa = esd->GetEventplane();
  TVector2 *standardQ = esdEPa->GetQVector();  // This is the "standard" Q-Vector
  Double_t qx = -1.0;
  Double_t qy = -1.0;
  if(standardQ) {
    qx = standardQ->X();
    qy = standardQ->Y();
  }  
  TVector2 qVectorfortrack;
  qVectorfortrack.Set(qx,qy);
  Float_t eventPlanea = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.;   
  
  TVector2 *qsub1a = esdEPa->GetQsub1();
  TVector2 *qsub2a = esdEPa->GetQsub2();
  Float_t eventPlanesub1a = -100.0;
  Float_t eventPlanesub2a = -100.0;
  if(qsub1a) eventPlanesub1a = TVector2::Phi_0_2pi(qsub1a->Phi())/2.;
  if(qsub2a) eventPlanesub2a = TVector2::Phi_0_2pi(qsub2a->Phi())/2.;
  Double_t diffsub1sub2a = -100.0;
  if(qsub1a && qsub2a) {
    diffsub1sub2a = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
  }

  //////////////////////////////////////////
  // Cut for event with TPC event defined
  /////////////////////////////////////////
  
  if(!standardQ) {
    PostData(1, fListHist);
    return;
  }

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

  fHistEV->Fill(binctt,2.0);

  // Fill
  valuensparsea[0] = eventPlanea;  
  fEventPlane->Fill(&valuensparsea[0]);

  valuensparsef[0] = diffsub1sub2a;
  fCosRes->Fill(&valuensparsef[0]);
  
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

    // Subtract the tracks from the event plane
    Double_t qX = standardQ->X() - esdEPa->GetQContributionX(track);  //Modify the components: subtract the track you want to look at with your analysis
    Double_t qY = standardQ->Y() - esdEPa->GetQContributionY(track);  //Modify the components: subtract the track you want to look at with your analysis
    TVector2 newQVectorfortrack;
    newQVectorfortrack.Set(qX,qY);
    Float_t eventplanesubtracted = TVector2::Phi_0_2pi(newQVectorfortrack.Phi())/2; 

    ////////////////////////////////////////
    // Fill pt and eta for the THnSparseF
    ///////////////////////////////////////

    valuensparsee[2] = track->Pt();
    valuensparsee[3] = track->Eta();    
    valuensparseg[2] = track->Pt();
    valuensparseh[2] = track->Pt();

    Bool_t fillTPC = kTRUE;
    if((!qsub1a) || (!qsub2a)) fillTPC = kFALSE;

    if(fSubEtaGapTPC) {
      if(track->Eta() < (- fEtaGap/2.)) eventplanesubtracted = eventPlanesub1a;
      else if(track->Eta() > (fEtaGap/2.)) eventplanesubtracted = eventPlanesub2a;
      else fillTPC = kFALSE;
    }

    ///////////////
    // MC
    //////////////
    if(fUseMCReactionPlane) {
      eventplanesubtracted = mcReactionPlane;
      fillTPC = kTRUE;
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
    valuensparsea[0] = eventplanesubtracted;
    if(fillTPC) fEventPlaneaftersubtraction->Fill(&valuensparsea[0]);
    
    //
    valuensparsee[0] = TMath::Cos(2*phitrack);
    fCos2phie->Fill(&valuensparsee[0]);
    valuensparsee[0] = TMath::Sin(2*phitrack);
    fSin2phie->Fill(&valuensparsee[0]);
    
    //
    valuensparsee[0] = TMath::Cos(2*eventplanesubtracted);
    if(fillTPC) fCos2phiep->Fill(&valuensparsee[0]);
    valuensparsee[0] = TMath::Sin(2*eventplanesubtracted);
    if(fillTPC) fSin2phiep->Fill(&valuensparsee[0]);
    valuensparsee[0] = TMath::Sin(2*TVector2::Phi_mpi_pi(phitrack-eventplanesubtracted));
    if(fillTPC) fSin2phiephiep->Fill(&valuensparsee[0]);
    
    // 
    valuensparseg[0] = deltaphi;
    if(fillTPC) fDeltaPhiMaps->Fill(&valuensparseg[0]);
    
    //
    valuensparseh[0] = TMath::Cos(2*TVector2::Phi_mpi_pi(phitrack-eventplanesubtracted));
    if(fillTPC) fCosPhiMaps->Fill(&valuensparseh[0]);
    
  }

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////AFTERBURNER
  if (fAfterBurnerOn)
    {
      fflowEvent->AddFlow(fV1,fV2,fV3,fV4,fV5);     //add flow
      fflowEvent->CloneTracks(fNonFlowNumberOfTrackClones); //add nonflow by cloning tracks
    }
  //////////////////////////////////////////////////////////////////////////////



  if((0.0< cntr) && (cntr<20.0))  PostData(2,fflowEvent);
  if((20.0< cntr) && (cntr<40.0))  PostData(3,fflowEvent);
  if((40.0< cntr) && (cntr<60.0))  PostData(4,fflowEvent);
  if((60.0< cntr) && (cntr<80.0))  PostData(5,fflowEvent);
  
 
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
