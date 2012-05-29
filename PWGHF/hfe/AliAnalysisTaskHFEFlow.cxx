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
#include "TLorentzVector.h"
#include "TParticle.h"

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
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODTrack.h"
#include "AliStack.h"
#include "AliMCEvent.h"

#include "AliFlowCandidateTrack.h"
#include "AliFlowEvent.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowVector.h"
#include "AliFlowCommonConstants.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliHFEVZEROEventPlane.h"

#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskHFEFlow.h"


//____________________________________________________________________
AliAnalysisTaskHFEFlow::AliAnalysisTaskHFEFlow() :
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
  fNbBinsCentralityQCumulant(4),
  fNbBinsPtQCumulant(12),
  fMinPtQCumulant(0.2),
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
  fNoPID(kFALSE),
  fChi2OverNDFCut(3.0),
  fMaxdca(3.0),
  fMaxopeningtheta(0.02),
  fMaxopeningphi(0.1),
  fMaxopening3D(0.1),
  fMaxInvmass(0.1),
  fSetMassConstraint(kFALSE),
  fDebugLevel(0),
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fPID(0),
  fPIDqa(0),
  fflowEvent(NULL),
  fHFEBackgroundCuts(0),
  fPIDBackground(0),
  fPIDBackgroundqa(0),
  fAlgorithmMA(kTRUE),
  fArraytrack(NULL),
  fCounterPoolBackground(0),
  fHFEVZEROEventPlane(0x0),
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
  fSinResabc(0x0),
  fProfileCosResab(0x0),
  fProfileCosResac(0x0),
  fProfileCosResbc(0x0),
  fCosRes(0x0),
  fSinRes(0x0),
  fProfileCosRes(0x0),
  fTrackingCuts(0x0),
  fDeltaPhiMapsBeforePID(0x0),
  fCosPhiMapsBeforePID(0x0),
  fDeltaPhiMaps(0x0),
  fCosPhiMaps(0x0),
  fProfileCosPhiMaps(0x0),
  fDeltaPhiMapsTaggedPhotonic(0x0),
  fCosPhiMapsTaggedPhotonic(0x0),
  fDeltaPhiMapsTaggedNonPhotonic(0x0),
  fCosPhiMapsTaggedNonPhotonic(0x0),
  fDeltaPhiMapsTaggedPhotonicLS(0x0),
  fCosPhiMapsTaggedPhotonicLS(0x0),
  fMCSourceDeltaPhiMaps(0x0),
  fOppSignDeltaPhiMaps(0x0),
  fSameSignDeltaPhiMaps(0x0),
  fOppSignAngle(0x0),
  fSameSignAngle(0x0)
{
  // Constructor

  for(Int_t k = 0; k < 10; k++) {
    fBinCentralityLess[k] = 0.0;
  }
  
}
//______________________________________________________________________________
AliAnalysisTaskHFEFlow:: AliAnalysisTaskHFEFlow(const char *name) :
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
  fNoPID(kFALSE),
  fChi2OverNDFCut(3.0),
  fMaxdca(3.0),
  fMaxopeningtheta(0.02),
  fMaxopeningphi(0.1),
  fMaxopening3D(0.1),
  fMaxInvmass(0.1),
  fSetMassConstraint(kFALSE),
  fDebugLevel(0),
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fPID(0),
  fPIDqa(0),
  fflowEvent(NULL),
  fHFEBackgroundCuts(0),
  fPIDBackground(0),
  fPIDBackgroundqa(0),
  fAlgorithmMA(kTRUE),  
  fArraytrack(NULL),
  fCounterPoolBackground(0),
  fHFEVZEROEventPlane(0x0),
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
  fSinResabc(0x0),
  fProfileCosResab(0x0),
  fProfileCosResac(0x0),
  fProfileCosResbc(0x0),
  fCosRes(0x0),
  fSinRes(0x0),
  fProfileCosRes(0x0),
  fTrackingCuts(0x0),
  fDeltaPhiMapsBeforePID(0x0),
  fCosPhiMapsBeforePID(0x0),
  fDeltaPhiMaps(0x0),
  fCosPhiMaps(0x0),
  fProfileCosPhiMaps(0x0),
  fDeltaPhiMapsTaggedPhotonic(0x0),
  fCosPhiMapsTaggedPhotonic(0x0),
  fDeltaPhiMapsTaggedNonPhotonic(0x0),
  fCosPhiMapsTaggedNonPhotonic(0x0),
  fDeltaPhiMapsTaggedPhotonicLS(0x0),
  fCosPhiMapsTaggedPhotonicLS(0x0),
  fMCSourceDeltaPhiMaps(0x0),
  fOppSignDeltaPhiMaps(0x0),
  fSameSignDeltaPhiMaps(0x0),
  fOppSignAngle(0x0),
  fSameSignAngle(0x0)
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

  fPIDBackground = new AliHFEpid("hfePidBackground");
  fPIDBackgroundqa = new AliHFEpidQAmanager;

  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
  for(Int_t bincless = 0; bincless < fNbBinsCentralityQCumulant; bincless++) {
    DefineOutput(bincless+2,AliFlowEventSimple::Class()); 
  }
  
}
//____________________________________________________________
AliAnalysisTaskHFEFlow::AliAnalysisTaskHFEFlow(const AliAnalysisTaskHFEFlow &ref):
  AliAnalysisTaskSE(ref),
  fListHist(NULL),
  fAODAnalysis(ref.fAODAnalysis), 
  fUseFlagAOD(ref.fUseFlagAOD),
  fApplyCut(ref.fApplyCut),
  fFlags(ref.fFlags),
  fVZEROEventPlane(ref.fVZEROEventPlane),
  fVZEROEventPlaneA(ref.fVZEROEventPlaneA),
  fVZEROEventPlaneC(ref.fVZEROEventPlaneC),
  fSubEtaGapTPC(ref.fSubEtaGapTPC),
  fEtaGap(ref.fEtaGap),
  fNbBinsCentralityQCumulant(ref.fNbBinsCentralityQCumulant),
  fNbBinsPtQCumulant(ref.fNbBinsPtQCumulant),
  fMinPtQCumulant(ref.fMinPtQCumulant),
  fMaxPtQCumulant(ref.fMaxPtQCumulant),
  fAfterBurnerOn(ref.fAfterBurnerOn),
  fNonFlowNumberOfTrackClones(ref.fNonFlowNumberOfTrackClones),
  fV1(ref.fV1),
  fV2(ref.fV2),
  fV3(ref.fV3),
  fV4(ref.fV4),
  fV5(ref.fV5),
  fMaxNumberOfIterations(ref.fMaxNumberOfIterations),
  fPrecisionPhi(ref.fPrecisionPhi),
  fUseMCReactionPlane(ref.fUseMCReactionPlane),
  fMCPID(ref.fMCPID),
  fNoPID(ref.fNoPID),
  fChi2OverNDFCut(ref.fChi2OverNDFCut),
  fMaxdca(ref.fMaxdca),
  fMaxopeningtheta(ref.fMaxopeningtheta),
  fMaxopeningphi(ref.fMaxopeningphi),
  fMaxopening3D(ref.fMaxopening3D),
  fMaxInvmass(ref.fMaxInvmass),
  fSetMassConstraint(ref.fSetMassConstraint),
  fDebugLevel(ref.fDebugLevel),
  fcutsRP(NULL),
  fcutsPOI(NULL),
  fHFECuts(NULL),
  fPID(NULL),
  fPIDqa(NULL),
  fflowEvent(NULL),
  fHFEBackgroundCuts(NULL),
  fPIDBackground(NULL),
  fPIDBackgroundqa(NULL),
  fAlgorithmMA(ref.fAlgorithmMA),
  fArraytrack(NULL),
  fCounterPoolBackground(ref.fCounterPoolBackground),
  fHFEVZEROEventPlane(NULL),
  fHistEV(NULL),
  fEventPlane(NULL),
  fEventPlaneaftersubtraction(NULL),
  fCosSin2phiep(NULL),
  fCos2phie(NULL),
  fSin2phie(NULL),
  fCos2phiep(NULL),
  fSin2phiep(NULL),
  fSin2phiephiep(NULL),
  fCosResabc(NULL),
  fSinResabc(NULL),
  fProfileCosResab(NULL),
  fProfileCosResac(NULL),
  fProfileCosResbc(NULL),
  fCosRes(NULL),
  fSinRes(NULL),
  fProfileCosRes(NULL),
  fTrackingCuts(NULL),
  fDeltaPhiMapsBeforePID(NULL),
  fCosPhiMapsBeforePID(NULL),
  fDeltaPhiMaps(NULL),
  fCosPhiMaps(NULL),
  fProfileCosPhiMaps(NULL),
  fDeltaPhiMapsTaggedPhotonic(NULL),
  fCosPhiMapsTaggedPhotonic(NULL),
  fDeltaPhiMapsTaggedNonPhotonic(NULL),
  fCosPhiMapsTaggedNonPhotonic(NULL),
  fDeltaPhiMapsTaggedPhotonicLS(NULL),
  fCosPhiMapsTaggedPhotonicLS(NULL),
  fMCSourceDeltaPhiMaps(NULL),
  fOppSignDeltaPhiMaps(NULL),
  fSameSignDeltaPhiMaps(NULL),
  fOppSignAngle(NULL),
  fSameSignAngle(NULL)
{
  //
  // Copy Constructor
  //
  ref.Copy(*this);
}

//____________________________________________________________
AliAnalysisTaskHFEFlow &AliAnalysisTaskHFEFlow::operator=(const AliAnalysisTaskHFEFlow &ref){
  //
  // Assignment operator
  //
  if(this == &ref) 
    ref.Copy(*this);
  return *this;
}

//____________________________________________________________
void AliAnalysisTaskHFEFlow::Copy(TObject &o) const {
  // 
  // Copy into object o
  //
  AliAnalysisTaskHFEFlow &target = dynamic_cast<AliAnalysisTaskHFEFlow &>(o);
  target.fAODAnalysis = fAODAnalysis;
  target.fUseFlagAOD = fUseFlagAOD;
  target.fApplyCut = fApplyCut;
  target.fFlags = fFlags;
  target.fVZEROEventPlane = fVZEROEventPlane;
  target.fVZEROEventPlaneA = fVZEROEventPlaneA;
  target.fVZEROEventPlaneC = fVZEROEventPlaneC;
  target.fSubEtaGapTPC = fSubEtaGapTPC;
  target.fEtaGap = fEtaGap;
  target.fNbBinsCentralityQCumulant = fNbBinsCentralityQCumulant;
  target.fNbBinsPtQCumulant = fNbBinsPtQCumulant;
  target.fMinPtQCumulant = fMinPtQCumulant;
  target.fMaxPtQCumulant = fMaxPtQCumulant;
  target.fAfterBurnerOn = fAfterBurnerOn;
  target.fNonFlowNumberOfTrackClones = fNonFlowNumberOfTrackClones;
  target.fV1 = fV1;
  target.fV2 = fV2;
  target.fV3 = fV3;
  target.fV4 = fV4;
  target.fV5 = fV5;
  target.fMaxNumberOfIterations = fMaxNumberOfIterations;
  target.fPrecisionPhi = fPrecisionPhi;
  target.fUseMCReactionPlane = fUseMCReactionPlane;
  target.fMCPID = fMCPID;
  target.fNoPID = fNoPID;
  target.fChi2OverNDFCut = fChi2OverNDFCut;
  target.fMaxdca = fMaxdca;
  target.fMaxopeningtheta = fMaxopeningtheta;
  target.fMaxopeningphi = fMaxopeningphi;
  target.fMaxopening3D = fMaxopening3D;
  target.fMaxInvmass = fMaxInvmass;
  target.fSetMassConstraint =  fSetMassConstraint;
  target.fAlgorithmMA = fAlgorithmMA;
  target.fCounterPoolBackground = fCounterPoolBackground;
  target.fDebugLevel = fDebugLevel;
  target.fcutsRP = fcutsRP;
  target.fcutsPOI = fcutsPOI;
  target.fHFECuts = fHFECuts;
  target.fPID = fPID;
  target.fPIDqa = fPIDqa;
  target.fHFEVZEROEventPlane = fHFEVZEROEventPlane;
 
}
//____________________________________________________________
AliAnalysisTaskHFEFlow::~AliAnalysisTaskHFEFlow(){
  //
  // Destructor
  //
  if(fArraytrack) delete fArraytrack;
  if(fListHist) delete fListHist;
  if(fcutsRP) delete fcutsRP;
  if(fcutsPOI) delete fcutsPOI;
  if(fHFECuts) delete fHFECuts;
  if(fPID) delete fPID;
  if(fPIDqa) delete fPIDqa;
  if(fflowEvent) delete fflowEvent;
  if(fHFEBackgroundCuts) delete fHFEBackgroundCuts;
  if(fPIDBackground) delete fPIDBackground;
  if(fPIDBackgroundqa) delete fPIDBackgroundqa;
  //if(fHFEVZEROEventPlane) delete fHFEVZEROEventPlane;
 

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
 
  // AOD or ESD
  AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
    SetAODAnalysis(kTRUE);
    //printf("Put AOD analysis on\n");
  } else {
    SetAODAnalysis(kFALSE);
  }


  // RP TRACK CUTS:
  fcutsRP = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();
  fcutsRP->SetName("StandartTPC");
  fcutsRP->SetEtaRange(-0.9,0.9);
  fcutsRP->SetQA(kTRUE);
  //TList *qaCutsRP = fcutsRP->GetQA();
  //qaCutsRP->SetName("QA_StandartTPC_RP");

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
  if(fAODAnalysis) fHFECuts->SetAOD();  

  // PID HFE
  //fPID->SetHasMCData(HasMCData());
  if(!fPID->GetNumberOfPIDdetectors()) fPID->AddDetector("TPC", 0);
  fPID->InitializePID();
  fPIDqa->Initialize(fPID);
  fPID->SortDetectors();

  // HFE Background cuts

  if(!fHFEBackgroundCuts){
     fHFEBackgroundCuts = new AliESDtrackCuts();
     fHFEBackgroundCuts->SetName("nackgroundcuts");
     //Configure Default Track Cuts
     fHFEBackgroundCuts->SetAcceptKinkDaughters(kFALSE);
     fHFEBackgroundCuts->SetRequireTPCRefit(kTRUE);
     fHFEBackgroundCuts->SetEtaRange(-0.9,0.9);
     fHFEBackgroundCuts->SetRequireSigmaToVertex(kTRUE);
     fHFEBackgroundCuts->SetMaxChi2PerClusterTPC(4.0);
     fHFEBackgroundCuts->SetMinNClustersTPC(50);
     fHFEBackgroundCuts->SetPtRange(0.3,1e10);
  }
  
  // PID background HFE
  if(!fPIDBackground->GetNumberOfPIDdetectors()) fPIDBackground->AddDetector("TPC", 0);
  fPIDBackground->InitializePID();
  fPIDBackgroundqa->Initialize(fPIDBackground);
  fPIDBackground->SortDetectors();
  


  //**************************
  // Bins for the THnSparse
  //**************************

  Int_t nBinsPt = 44;
  Double_t minPt = 0.1;
  Double_t maxPt = 20.0;
  Double_t binLimLogPt[nBinsPt+1];
  Double_t binLimPt[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++) binLimLogPt[i]=(Double_t)TMath::Log10(minPt) + (TMath::Log10(maxPt)-TMath::Log10(minPt))/nBinsPt*(Double_t)i ;
  for(Int_t i=0; i<=nBinsPt; i++) binLimPt[i]=(Double_t)TMath::Power(10,binLimLogPt[i]);

  Int_t nBinsPtPlus = fNbBinsPtQCumulant;
  Double_t minPtPlus = fMinPtQCumulant;
  Double_t maxPtPlus = fMaxPtQCumulant;
  Double_t binLimPtPlus[nBinsPtPlus+1];
  for(Int_t i=0; i<=nBinsPtPlus; i++) binLimPtPlus[i]=(Double_t)minPtPlus + (maxPtPlus-minPtPlus)/nBinsPtPlus*(Double_t)i ;

  Int_t nBinsEta = 8;
  Double_t minEta = -0.8;
  Double_t maxEta = 0.8;
  Double_t binLimEta[nBinsEta+1];
  for(Int_t i=0; i<=nBinsEta; i++) binLimEta[i]=(Double_t)minEta + (maxEta-minEta)/nBinsEta*(Double_t)i ;

  Int_t nBinsStep = 7;
  Double_t minStep = 0.;
  Double_t maxStep = 7.;
  Double_t binLimStep[nBinsStep+1];
  for(Int_t i=0; i<=nBinsStep; i++) binLimStep[i]=(Double_t)minStep + (maxStep-minStep)/nBinsStep*(Double_t)i ;

  Int_t nBinsEtaLess = 2;
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

  Int_t nBinsPhi = 12;
  Double_t minPhi = 0.0;
  Double_t maxPhi = TMath::Pi();
  Double_t binLimPhi[nBinsPhi+1];
  for(Int_t i=0; i<=nBinsPhi; i++) {
    binLimPhi[i]=(Double_t)minPhi + (maxPhi-minPhi)/nBinsPhi*(Double_t)i ;
    //printf("bin phi is %f for %d\n",binLimPhi[i],i);
  }

  Int_t nBinsAngle = 40;
  Double_t minAngle = 0.0;
  Double_t maxAngle = 1.0;
  Double_t binLimAngle[nBinsAngle+1];
  for(Int_t i=0; i<=nBinsAngle; i++) {
    binLimAngle[i]=(Double_t)minAngle + (maxAngle-minAngle)/nBinsAngle*(Double_t)i ;
    //printf("bin phi is %f for %d\n",binLimPhi[i],i);
  }

  Int_t nBinsCharge = 2;
  Double_t minCharge = -1.0;
  Double_t maxCharge = 1.0;
  Double_t binLimCharge[nBinsCharge+1];
  for(Int_t i=0; i<=nBinsCharge; i++) binLimCharge[i]=(Double_t)minCharge + (maxCharge-minCharge)/nBinsCharge*(Double_t)i ;

  Int_t nBinsSource = 10;
  Double_t minSource = 0.;
  Double_t maxSource = 10.;
  Double_t binLimSource[nBinsSource+1];
  for(Int_t i=0; i<=nBinsSource; i++) binLimSource[i]=(Double_t)minSource + (maxSource-minSource)/nBinsSource*(Double_t)i ;

  Int_t nBinsInvMass = 50;
  Double_t minInvMass = 0.;
  Double_t maxInvMass = 0.3;
  Double_t binLimInvMass[nBinsInvMass+1];
  for(Int_t i=0; i<=nBinsInvMass; i++) binLimInvMass[i]=(Double_t)minInvMass + (maxInvMass-minInvMass)/nBinsInvMass*(Double_t)i ;
  
  //******************
  // Histograms
  //******************
    
  fListHist = new TList();
  fListHist->SetOwner();

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

  const Int_t nDimfbiss=4;
  Int_t nBinfbiss[nDimfbiss] = {nBinsCos,nBinsCos,nBinsCos,nBinsC};
  fSinResabc = new THnSparseF("SinRes_abc","SinRes_abc",nDimfbiss,nBinfbiss);
  fSinResabc->SetBinEdges(0,binLimCos);
  fSinResabc->SetBinEdges(1,binLimCos);
  fSinResabc->SetBinEdges(2,binLimCos);
  fSinResabc->SetBinEdges(3,binLimC);
  fSinResabc->Sumw2();

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

  const Int_t nDimff=2;
  Int_t nBinff[nDimff] = {nBinsCos, nBinsC};
  fSinRes = new THnSparseF("SinRes","SinRes",nDimff,nBinff);
  fSinRes->SetBinEdges(0,binLimCos);
  fSinRes->SetBinEdges(1,binLimC);
  fSinRes->Sumw2();

  // Profile cosres centrality
  fProfileCosRes = new TProfile("ProfileCosRes","ProfileCosRes",nBinsCMore,binLimCMore);
  fProfileCosRes->Sumw2();

  // Debugging tracking steps
  const Int_t nDimTrStep=2;
  Int_t nBinTrStep[nDimTrStep] = {nBinsPt,nBinsStep};
  fTrackingCuts = new THnSparseF("TrackingCuts","TrackingCuts",nDimTrStep,nBinTrStep);
  fTrackingCuts->SetBinEdges(0,binLimPt);
  fTrackingCuts->SetBinEdges(1,binLimStep);
  fTrackingCuts->Sumw2();
  
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

  const Int_t nDimgb=3;
  Int_t nBingb[nDimgb] = {nBinsPhi,nBinsC,nBinsPt};

  fDeltaPhiMapsBeforePID = new THnSparseF("DeltaPhiMapsBeforePID","DeltaPhiMapsBeforePID",nDimgb,nBingb);
  fDeltaPhiMapsBeforePID->SetBinEdges(0,binLimPhi);
  fDeltaPhiMapsBeforePID->SetBinEdges(1,binLimC);
  fDeltaPhiMapsBeforePID->SetBinEdges(2,binLimPt);
  fDeltaPhiMapsBeforePID->Sumw2();  

 
  fDeltaPhiMapsTaggedPhotonic = new THnSparseF("DeltaPhiMapsTaggedPhotonic","DeltaPhiMapsTaggedPhotonic",nDimgb,nBingb);
  fDeltaPhiMapsTaggedPhotonic->SetBinEdges(0,binLimPhi);
  fDeltaPhiMapsTaggedPhotonic->SetBinEdges(1,binLimC);
  fDeltaPhiMapsTaggedPhotonic->SetBinEdges(2,binLimPt);
  fDeltaPhiMapsTaggedPhotonic->Sumw2();  

  fDeltaPhiMapsTaggedNonPhotonic = new THnSparseF("DeltaPhiMapsTaggedNonPhotonic","DeltaPhiMapsTaggedNonPhotonic",nDimgb,nBingb);
  fDeltaPhiMapsTaggedNonPhotonic->SetBinEdges(0,binLimPhi);
  fDeltaPhiMapsTaggedNonPhotonic->SetBinEdges(1,binLimC);
  fDeltaPhiMapsTaggedNonPhotonic->SetBinEdges(2,binLimPt);
  fDeltaPhiMapsTaggedNonPhotonic->Sumw2();  

  fDeltaPhiMapsTaggedPhotonicLS = new THnSparseF("DeltaPhiMapsTaggedPhotonicLS","DeltaPhiMapsTaggedPhotonicLS",nDimgb,nBingb);
  fDeltaPhiMapsTaggedPhotonicLS->SetBinEdges(0,binLimPhi);
  fDeltaPhiMapsTaggedPhotonicLS->SetBinEdges(1,binLimC);
  fDeltaPhiMapsTaggedPhotonicLS->SetBinEdges(2,binLimPt);
  fDeltaPhiMapsTaggedPhotonicLS->Sumw2();  

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

  const Int_t nDimhb=3;
  Int_t nBinhb[nDimhb] = {nBinsCos,nBinsC,nBinsPt};

  fCosPhiMapsBeforePID = new THnSparseF("CosPhiMapsBeforePID","CosPhiMapsBeforePID",nDimhb,nBinhb);
  fCosPhiMapsBeforePID->SetBinEdges(0,binLimCos);
  fCosPhiMapsBeforePID->SetBinEdges(1,binLimC);
  fCosPhiMapsBeforePID->SetBinEdges(2,binLimPt);
  fCosPhiMapsBeforePID->Sumw2();

  fCosPhiMapsTaggedPhotonic = new THnSparseF("CosPhiMapsTaggedPhotonic","CosPhiMapsTaggedPhotonic",nDimhb,nBinhb);
  fCosPhiMapsTaggedPhotonic->SetBinEdges(0,binLimCos);
  fCosPhiMapsTaggedPhotonic->SetBinEdges(1,binLimC);
  fCosPhiMapsTaggedPhotonic->SetBinEdges(2,binLimPt);
  fCosPhiMapsTaggedPhotonic->Sumw2();

  fCosPhiMapsTaggedNonPhotonic = new THnSparseF("CosPhiMapsTaggedNonPhotonic","CosPhiMapsTaggedNonPhotonic",nDimhb,nBinhb);
  fCosPhiMapsTaggedNonPhotonic->SetBinEdges(0,binLimCos);
  fCosPhiMapsTaggedNonPhotonic->SetBinEdges(1,binLimC);
  fCosPhiMapsTaggedNonPhotonic->SetBinEdges(2,binLimPt);
  fCosPhiMapsTaggedNonPhotonic->Sumw2();
  
  fCosPhiMapsTaggedPhotonicLS = new THnSparseF("CosPhiMapsTaggedPhotonicLS","CosPhiMapsTaggedPhotonicLS",nDimhb,nBinhb);
  fCosPhiMapsTaggedPhotonicLS->SetBinEdges(0,binLimCos);
  fCosPhiMapsTaggedPhotonicLS->SetBinEdges(1,binLimC);
  fCosPhiMapsTaggedPhotonicLS->SetBinEdges(2,binLimPt);
  fCosPhiMapsTaggedPhotonicLS->Sumw2();

  // Profile Maps cos phi
  fProfileCosPhiMaps = new TProfile2D("ProfileCosPhiMaps","ProfileCosPhiMaps",nBinsC,binLimC,nBinsPt,binLimPt);
  fProfileCosPhiMaps->Sumw2();

  // Background study
  const Int_t nDimMCSource=3;
  Int_t nBinMCSource[nDimMCSource] = {nBinsC,nBinsPt,nBinsSource};
  fMCSourceDeltaPhiMaps = new THnSparseF("MCSourceDeltaPhiMaps","MCSourceDeltaPhiMaps",nDimMCSource,nBinMCSource);
  fMCSourceDeltaPhiMaps->SetBinEdges(0,binLimC);
  fMCSourceDeltaPhiMaps->SetBinEdges(1,binLimPt);
  fMCSourceDeltaPhiMaps->SetBinEdges(2,binLimSource);
  fMCSourceDeltaPhiMaps->Sumw2();

  // Maps invmass opposite
  const Int_t nDimOppSign=5;
  Int_t nBinOppSign[nDimOppSign] = {nBinsPhi,nBinsC,nBinsPt,nBinsInvMass,nBinsSource};
  fOppSignDeltaPhiMaps = new THnSparseF("OppSignDeltaPhiMaps","OppSignDeltaPhiMaps",nDimOppSign,nBinOppSign);
  fOppSignDeltaPhiMaps->SetBinEdges(0,binLimPhi);
  fOppSignDeltaPhiMaps->SetBinEdges(1,binLimC);
  fOppSignDeltaPhiMaps->SetBinEdges(2,binLimPt);
  fOppSignDeltaPhiMaps->SetBinEdges(3,binLimInvMass);
  fOppSignDeltaPhiMaps->SetBinEdges(4,binLimSource);
  fOppSignDeltaPhiMaps->Sumw2();

  // Maps invmass same sign
  const Int_t nDimSameSign=5;
  Int_t nBinSameSign[nDimSameSign] = {nBinsPhi,nBinsC,nBinsPt,nBinsInvMass,nBinsSource};
  fSameSignDeltaPhiMaps = new THnSparseF("SameSignDeltaPhiMaps","SameSignDeltaPhiMaps",nDimSameSign,nBinSameSign);
  fSameSignDeltaPhiMaps->SetBinEdges(0,binLimPhi);
  fSameSignDeltaPhiMaps->SetBinEdges(1,binLimC);
  fSameSignDeltaPhiMaps->SetBinEdges(2,binLimPt);
  fSameSignDeltaPhiMaps->SetBinEdges(3,binLimInvMass);
  fSameSignDeltaPhiMaps->SetBinEdges(4,binLimSource);
  fSameSignDeltaPhiMaps->Sumw2();

  // Maps angle same sign
  const Int_t nDimAngleSameSign=3;
  Int_t nBinAngleSameSign[nDimAngleSameSign] = {nBinsAngle,nBinsC,nBinsSource};
  fSameSignAngle = new THnSparseF("SameSignAngleMaps","SameSignAngleMaps",nDimAngleSameSign,nBinAngleSameSign);
  fSameSignAngle->SetBinEdges(0,binLimAngle);
  fSameSignAngle->SetBinEdges(1,binLimC);
  fSameSignAngle->SetBinEdges(2,binLimSource);
  fSameSignAngle->Sumw2();

  // Maps angle opp sign
  const Int_t nDimAngleOppSign=3;
  Int_t nBinAngleOppSign[nDimAngleOppSign] = {nBinsAngle,nBinsC,nBinsSource};
  fOppSignAngle = new THnSparseF("OppSignAngleMaps","OppSignAngleMaps",nDimAngleOppSign,nBinAngleOppSign);
  fOppSignAngle->SetBinEdges(0,binLimAngle);
  fOppSignAngle->SetBinEdges(1,binLimC);
  fOppSignAngle->SetBinEdges(2,binLimSource);
  fOppSignAngle->Sumw2();


  //**************************
  // Add to the list
  //******************************

  //fListHist->Add(qaCutsRP);
  fListHist->Add(fPIDqa->MakeList("HFEpidQA"));
  fListHist->Add(fPIDBackgroundqa->MakeList("HFEpidBackgroundQA"));
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
  fListHist->Add(fSinRes);
  fListHist->Add(fSinResabc);
  fListHist->Add(fTrackingCuts);
  fListHist->Add(fDeltaPhiMapsBeforePID);
  fListHist->Add(fCosPhiMapsBeforePID);
  fListHist->Add(fDeltaPhiMaps);
  fListHist->Add(fCosPhiMaps);
  fListHist->Add(fProfileCosPhiMaps);
  fListHist->Add(fDeltaPhiMapsTaggedPhotonic);
  fListHist->Add(fCosPhiMapsTaggedPhotonic);
  fListHist->Add(fDeltaPhiMapsTaggedNonPhotonic);
  fListHist->Add(fCosPhiMapsTaggedNonPhotonic);
  fListHist->Add(fDeltaPhiMapsTaggedPhotonicLS);
  fListHist->Add(fCosPhiMapsTaggedPhotonicLS);
  fListHist->Add(fMCSourceDeltaPhiMaps);
  fListHist->Add(fOppSignDeltaPhiMaps);
  fListHist->Add(fSameSignDeltaPhiMaps);
  fListHist->Add(fSameSignAngle);
  fListHist->Add(fOppSignAngle);


  if(fHFEVZEROEventPlane && (fDebugLevel > 2)) fListHist->Add(fHFEVZEROEventPlane->GetOutputList());
  

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
  Double_t valuensparsefsin[2];
  Double_t valuensparsefbis[4];
  Double_t valuensparsefbissin[4];
  Double_t valuensparseg[5];
  Double_t valuensparseh[5];
  Double_t valuensparsehprofile[3];
  Double_t valuensparseMCSourceDeltaPhiMaps[3];
  Double_t valuetrackingcuts[2];
   
  AliMCEvent *mcEvent = MCEvent();
  AliMCParticle *mctrack = NULL;
    
  /////////////////
  // centrality
  /////////////////
  
  //AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  //if(!esd) return;
  AliCentrality *centrality = fInputEvent->GetCentrality();
  //printf("Got the centrality\n");
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
   
  
  if(binct > 11.0) return;
 
  // centrality
  valuensparsea[4] = binct;  
  valuensparseabis[1] = binct;  
  valuensparsee[1] = binct;    
  valuensparsef[1] = binctMore;  
  valuensparsefsin[1] = binct;  
  valuensparsefbis[3] = binctMore;  
  valuensparsefbissin[3] = binct;  
  valuensparseg[1] = binct;
  valuensparseh[1] = binct; 
  valuensparsehprofile[1] = binct; 
  valuecossinephiep[2] = binctMore;
  valuensparseMCSourceDeltaPhiMaps[0] = binct;
 
  //////////////////////
  // run number
  //////////////////////

  Int_t runnumber = fInputEvent->GetRunNumber();
  //printf("Run number %d\n",runnumber);
   
  if(!fPID->IsInitialized()){
    // Initialize PID with the given run number
    fPID->InitializePID(runnumber);
  }

  if(!fPIDBackground->IsInitialized()){
    // Initialize PID with the given run number
    fPIDBackground->InitializePID(runnumber);
  }

  fHFECuts->SetRecEvent(fInputEvent);
  if(mcEvent) fHFECuts->SetMCEvent(mcEvent);


  //////////
  // PID
  //////////
 
  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    //printf("No PID response set\n");
    return;
  }
  fPID->SetPIDResponse(pidResponse);
  fPIDBackground->SetPIDResponse(pidResponse);

  fHistEV->Fill(binctt,0.0);
 

  //////////////////
  // Event cut
  //////////////////
  if(!fHFECuts->CheckEventCuts("fEvRecCuts", fInputEvent)) {
    //printf("Do not pass the event cut\n");
    PostData(1, fListHist);
    return;
  }

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

  if(fHFEVZEROEventPlane && (!fAODAnalysis)){

    //AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
    //if(!esd) return;
    
    fHFEVZEROEventPlane->ProcessEvent(fInputEvent);
    
    if(TMath::Abs(fHFEVZEROEventPlane->GetEventPlaneV0A()+100) < 0.0000001) eventPlaneV0A = -100.0;
    else {
      eventPlaneV0A = TVector2::Phi_0_2pi(fHFEVZEROEventPlane->GetEventPlaneV0A());
      if(eventPlaneV0A > TMath::Pi()) eventPlaneV0A = eventPlaneV0A - TMath::Pi();
    }
    
    if(TMath::Abs(fHFEVZEROEventPlane->GetEventPlaneV0C()+100) < 0.0000001) eventPlaneV0C = -100.0;
    else {
      eventPlaneV0C = TVector2::Phi_0_2pi(fHFEVZEROEventPlane->GetEventPlaneV0C());
      if(eventPlaneV0C > TMath::Pi()) eventPlaneV0C = eventPlaneV0C - TMath::Pi();
    }

    if(TMath::Abs(fHFEVZEROEventPlane->GetEventPlaneV0()+100) < 0.0000001) eventPlaneV0 = -100.0;
    else {
      eventPlaneV0 = TVector2::Phi_0_2pi(fHFEVZEROEventPlane->GetEventPlaneV0());
      if(eventPlaneV0 > TMath::Pi()) eventPlaneV0 = eventPlaneV0 - TMath::Pi();
    }
    
  }
  else {
    
    eventPlaneV0 = TVector2::Phi_0_2pi(vEPa->GetEventplane("V0", fInputEvent,2));
    //printf("eventPlaneV0 %f\n",eventPlaneV0);
    if(eventPlaneV0 > TMath::Pi()) eventPlaneV0 = eventPlaneV0 - TMath::Pi();
    //printf("eventPlaneV0 %f\n",eventPlaneV0);
    eventPlaneV0A = TVector2::Phi_0_2pi(vEPa->GetEventplane("V0A", fInputEvent,2));
    if(eventPlaneV0A > TMath::Pi()) eventPlaneV0A = eventPlaneV0A - TMath::Pi();
    eventPlaneV0C = TVector2::Phi_0_2pi(vEPa->GetEventplane("V0C", fInputEvent,2));
    if(eventPlaneV0C > TMath::Pi()) eventPlaneV0C = eventPlaneV0C - TMath::Pi();
  
  }

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
  Double_t diffsub1sub2asin = -100.0;
  Double_t diffsubasubb = -100.0;
  Double_t diffsubasubc = -100.0;
  Double_t diffsubbsubc = -100.0;
  Double_t diffsubasubbsin = -100.0;
  Double_t diffsubasubcsin = -100.0;
  Double_t diffsubbsubcsin = -100.0;

  //if(fVZEROEventPlane) {
  diffsubasubb = TMath::Cos(2.*(eventPlaneV0A - eventPlaneV0C));
  diffsubasubc = TMath::Cos(2.*(eventPlaneV0A - eventPlaneTPC));
  diffsubbsubc = TMath::Cos(2.*(eventPlaneV0C - eventPlaneTPC));

  diffsubasubbsin = TMath::Sin(2.*(eventPlaneV0A - eventPlaneV0C));
  diffsubasubcsin = TMath::Sin(2.*(eventPlaneV0A - eventPlaneTPC));
  diffsubbsubcsin = TMath::Sin(2.*(eventPlaneV0C - eventPlaneTPC));
  //}
  //else {
  qsub1a = vEPa->GetQsub1();
  qsub2a = vEPa->GetQsub2();
  if(qsub1a) eventPlanesub1a = TVector2::Phi_0_2pi(qsub1a->Phi())/2.;
  if(qsub2a) eventPlanesub2a = TVector2::Phi_0_2pi(qsub2a->Phi())/2.;
  if(qsub1a && qsub2a) {
    diffsub1sub2a = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
    diffsub1sub2asin = TMath::Sin(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
  }
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

  if((!standardQ) || (!qsub1a) || (!qsub2a)) {
    //printf("No event plane\n");
    return;
  }

  ///////////////////////
  // AliFlowEvent  
  //////////////////////

  Int_t nbtracks = fInputEvent->GetNumberOfTracks();
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
    if(fDebugLevel > 5) fCosSin2phiep->Fill(&valuecossinephiep[0]);
    
    if(!fVZEROEventPlane) {
      valuensparsef[0] = diffsub1sub2a;
      fCosRes->Fill(&valuensparsef[0]);
      valuensparsefsin[0] = diffsub1sub2asin;
      if(fDebugLevel > 5) fSinRes->Fill(&valuensparsefsin[0]);
      if(fDebugLevel > 5) {
	fProfileCosRes->Fill(valuensparsef[1],valuensparsef[0]);
      }
    }
    else {
      valuensparsefbis[0] = diffsubasubb;
      valuensparsefbis[1] = diffsubbsubc;
      valuensparsefbis[2] = diffsubasubc;
      fCosResabc->Fill(&valuensparsefbis[0]);
      valuensparsefbissin[0] = diffsubasubbsin;
      valuensparsefbissin[1] = diffsubbsubcsin;
      valuensparsefbissin[2] = diffsubasubcsin;
      if(fDebugLevel > 5) fSinResabc->Fill(&valuensparsefbissin[0]);
      if(fDebugLevel > 5) {
	fProfileCosResab->Fill(valuensparsefbis[3],valuensparsefbis[0]);
	fProfileCosResac->Fill(valuensparsefbis[3],valuensparsefbis[1]);
	fProfileCosResbc->Fill(valuensparsefbis[3],valuensparsefbis[2]);
      }
    }
    
  }
  
  ////////////////////////////////////////
  // Loop to determine pool background
  /////////////////////////////////////////
  if(  fArraytrack ){ 
     fArraytrack->~TArrayI();
     new(fArraytrack) TArrayI(nbtracks);
  }
  else {  
    fArraytrack = new TArrayI(nbtracks);
  }
  fCounterPoolBackground = 0;
  for(Int_t k = 0; k < nbtracks; k++){
    
    AliVTrack *track = (AliVTrack *) fInputEvent->GetTrack(k);
    if(!track) continue;
    
    // Track cuts
    Bool_t survivedbackground = kTRUE;
    if(fAODAnalysis) {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
      if(aodtrack) {
	AliESDtrack esdTrack(aodtrack);
	// set the TPC cluster info
	esdTrack.SetTPCClusterMap(aodtrack->GetTPCClusterMap());
	esdTrack.SetTPCSharedMap(aodtrack->GetTPCSharedMap());
	esdTrack.SetTPCPointsF(aodtrack->GetTPCNclsF());
	// needed to calculate the impact parameters
	AliAODEvent *aodeventu = dynamic_cast<AliAODEvent *>(fInputEvent);
	if(aodeventu) {
	  AliAODVertex *vAOD = aodeventu->GetPrimaryVertex();
	  Double_t bfield = aodeventu->GetMagneticField();
	  Double_t pos[3],cov[6];
	  vAOD->GetXYZ(pos);
	  vAOD->GetCovarianceMatrix(cov);
	  const AliESDVertex vESD(pos,cov,100.,100);
	  esdTrack.RelateToVertex(&vESD,bfield,3.);
	} 
	if(!fHFEBackgroundCuts->IsSelected(&esdTrack)) {
	  survivedbackground = kFALSE;
	}
      }
    }
    else {
      AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
      if(esdtrack) {
	if(!fHFEBackgroundCuts->IsSelected(esdtrack)) survivedbackground = kFALSE;
      }
    }
    // PID
    if(survivedbackground) {
      // PID track cuts
      AliHFEpidObject hfetrack2;
      if(!fAODAnalysis) hfetrack2.SetAnalysisType(AliHFEpidObject::kESDanalysis);
      else hfetrack2.SetAnalysisType(AliHFEpidObject::kAODanalysis);
      hfetrack2.SetRecTrack(track);
      hfetrack2.SetCentrality((Int_t)binct);
      //printf("centrality %f and %d\n",binct,hfetrack.GetCentrality());
      hfetrack2.SetPbPb();
      if(fPIDBackground->IsSelected(&hfetrack2,0x0,"recTrackCont",fPIDBackgroundqa)) {
	fArraytrack->AddAt(k,fCounterPoolBackground);
	fCounterPoolBackground++;
	//printf("fCounterPoolBackground %d, track %d\n",fCounterPoolBackground,k);
      }
    }
  }

  // Look at kink mother in case of AOD
  Int_t numberofvertices = 1;
  AliAODEvent *aodevent = NULL;
  Int_t numberofmotherkink = 0;  
  if(fAODAnalysis) {
    aodevent = dynamic_cast<AliAODEvent *>(fInputEvent);
    if(aodevent) {
      numberofvertices = aodevent->GetNumberOfVertices();
    }
  }
  Double_t listofmotherkink[numberofvertices];
  if(fAODAnalysis && aodevent) {
    for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
      AliAODVertex *aodvertex = aodevent->GetVertex(ivertex);
      if(!aodvertex) continue;
      if(aodvertex->GetType()==AliAODVertex::kKink) {
	AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
	if(!mother) continue;
	Int_t idmother = mother->GetID();
	listofmotherkink[numberofmotherkink] = idmother;
	//printf("ID %d\n",idmother);
	numberofmotherkink++;
      }
    }
  }
  
  //////////////////////////
  // Loop over track
  //////////////////////////
  
  for(Int_t k = 0; k < nbtracks; k++){
      
    AliVTrack *track = (AliVTrack *) fInputEvent->GetTrack(k);
    if(!track) continue;
    
    if(fAODAnalysis) {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
      if(!aodtrack){
	AliError("AOD track is not there");
	continue;
      }  
      //printf("Find AOD track on\n");
      if(fUseFlagAOD){
	if(aodtrack->GetFlags() != fFlags) continue;  // Only process AOD tracks where the HFE is set
      }
    }
    
    if(fApplyCut) {

      valuetrackingcuts[0] = track->Pt(); 
      valuetrackingcuts[1] = 0;

      // RecKine: ITSTPC cuts  
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      if(fDebugLevel > 3) fTrackingCuts->Fill(&valuetrackingcuts[0]);
      
      // Reject kink mother
      if(fAODAnalysis) {
	Bool_t kinkmotherpass = kTRUE;
	for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
	  if(track->GetID() == listofmotherkink[kinkmother]) {
	    kinkmotherpass = kFALSE;
	    continue;
	  }
	}
	if(!kinkmotherpass) continue;
      }
      else {
	AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(track);
	if(esdtrack){  
	  if(esdtrack->GetKinkIndex(0) != 0) continue; 
	} // Quick and dirty fix to reject both kink mothers and daughters
      }
            
      valuetrackingcuts[1] = 1; 
       if(fDebugLevel > 3) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      // RecPrim
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      valuetrackingcuts[1] = 2; 
       if(fDebugLevel > 3) fTrackingCuts->Fill(&valuetrackingcuts[0]);    

      // HFEcuts: ITS layers cuts
       if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
       valuetrackingcuts[1] = 3; 
       if(fDebugLevel > 3) fTrackingCuts->Fill(&valuetrackingcuts[0]);     

      // HFE cuts: TOF PID and mismatch flag
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTOF + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      valuetrackingcuts[1] = 4; 
       if(fDebugLevel > 3) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      
      // HFE cuts: TPC PID cleanup
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTPC + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      valuetrackingcuts[1] = 5; 
       if(fDebugLevel > 3) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      
      // HFEcuts: Nb of tracklets TRD0
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      valuetrackingcuts[1] = 6; 
       if(fDebugLevel > 3) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      
    }

    //printf("Survived\n");

    /////////////////////////////////////////////////////////
    // Subtract candidate from TPC event plane
    ////////////////////////////////////////////////////////
    Float_t eventplanesubtracted = 0.0;    

    //if(eventplanedefined && (!fVZEROEventPlane)) {
    if(!fVZEROEventPlane) {
      // Subtract the tracks from the event plane
      Double_t qX = standardQ->X() - vEPa->GetQContributionX(track);  //Modify the components: subtract the track you want to look at with your analysis
      Double_t qY = standardQ->Y() - vEPa->GetQContributionY(track);  //Modify the components: subtract the track you want to look at with your analysis
      TVector2 newQVectorfortrack;
      newQVectorfortrack.Set(qX,qY);
      eventplanesubtracted = TVector2::Phi_0_2pi(newQVectorfortrack.Phi())/2; 
    }
    else eventplanesubtracted = eventPlanea;

    ///////////////////////////////////////////
    // Event plane
    //////////////////////////////////////////
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

    ////////////////////////////////////////
    // Define variables
    ///////////////////////////////////////

    valuensparsee[2] = track->Pt();
    valuensparsee[3] = track->Eta();    
    valuensparseg[2] = track->Pt();
    valuensparseh[2] = track->Pt();
    valuensparsehprofile[2] = track->Pt();
    valuensparseMCSourceDeltaPhiMaps[1] = track->Pt();
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

    //printf("charge %d\n",(Int_t)track->Charge());

    ////////////////////////
    // Fill before PID
    ///////////////////////
    
    if(fDebugLevel > 4) { 
      
      valuensparseg[0] = deltaphi;
      if(fillEventPlane) fDeltaPhiMapsBeforePID->Fill(&valuensparseg[0]);
      
      //
      valuensparseh[0] = TMath::Cos(2*TVector2::Phi_mpi_pi(phitrack-eventplanesubtracted));
      if(fillEventPlane) {
	fCosPhiMapsBeforePID->Fill(&valuensparseh[0]);
      }
    }
    
    ////////////////////////
    // Apply PID
    ////////////////////////
    if(!fNoPID) {
      // Apply PID for Data
      if(!fMCPID) {
	AliHFEpidObject hfetrack;
	if(!fAODAnalysis) hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
	else hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
	hfetrack.SetRecTrack(track);
	hfetrack.SetCentrality((Int_t)binct);
	//printf("centrality %f and %d\n",binct,hfetrack.GetCentrality());
	hfetrack.SetPbPb();
	if(!fPID->IsSelected(&hfetrack,0x0,"recTrackCont",fPIDqa)) {
	  continue;
	}
      }
      else {
	if(!mcEvent) continue;
	if(!(mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel()))))) continue;
	//printf("PdgCode %d\n",TMath::Abs(mctrack->Particle()->GetPdgCode()));
	if(TMath::Abs(mctrack->Particle()->GetPdgCode())!=11) continue;
      }
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
    
    
  
    /////////////////////
    // Fill THnSparseF
    /////////////////////

    //
    valuensparseabis[0] = eventplanesubtracted;
    if((fillEventPlane) && (fDebugLevel > 5)) fEventPlaneaftersubtraction->Fill(&valuensparseabis[0]);
    

    if(fDebugLevel > 5) 
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
    
    if(fDebugLevel > 1) {
      // background
      Int_t source = 0;
      Int_t indexmother = -1;
      source = FindMother(TMath::Abs(track->GetLabel()),mcEvent, indexmother);
      valuensparseMCSourceDeltaPhiMaps[2] = source;
      if(mcEvent) fMCSourceDeltaPhiMaps->Fill(&valuensparseMCSourceDeltaPhiMaps[0]);
      Int_t taggedvalue = LookAtNonHFE(k,track,fInputEvent,mcEvent,binct,deltaphi,source,indexmother);
      if(fillEventPlane) {
	// No opposite charge partner found in the invariant mass choosen
	if((taggedvalue!=2) && (taggedvalue!=6)) {
	  fDeltaPhiMapsTaggedNonPhotonic->Fill(&valuensparseg[0]);
	  fCosPhiMapsTaggedNonPhotonic->Fill(&valuensparseh[0]);
	}
	// One opposite charge partner found in the invariant mass choosen
	if((taggedvalue==2) || (taggedvalue==6)) {
	  fDeltaPhiMapsTaggedPhotonic->Fill(&valuensparseg[0]);
	  fCosPhiMapsTaggedPhotonic->Fill(&valuensparseh[0]);
	}
	// One same charge partner found in the invariant mass choosen
	if((taggedvalue==4) || (taggedvalue==6)) {
	  fDeltaPhiMapsTaggedPhotonicLS->Fill(&valuensparseg[0]);
	  fCosPhiMapsTaggedPhotonicLS->Fill(&valuensparseh[0]);
	}
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

  if(fArraytrack) {
    delete fArraytrack;
    fArraytrack = NULL;
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
//_____________________________________________________________________________________________
Int_t AliAnalysisTaskHFEFlow::LookAtNonHFE(Int_t iTrack1, AliVTrack *track1, AliVEvent *vEvent, AliMCEvent *mcEvent,Int_t binct,Double_t deltaphi,Int_t source,Int_t indexmother)
{	
  //
  // Look At Non HFE
  //

  // return -1 if nothing
  // return 2 if opposite charge within the mass range found
  // return 4 if like charge within the mass range found
  // return 6 if opposite charge and like charge within the mass range found
  //

  Int_t taggedphotonic = -1;

  Bool_t oppositetaggedphotonic = kFALSE;
  Bool_t sametaggedphotonic = kFALSE;

  //printf("fCounterPoolBackground %d in LookAtNonHFE!!!\n",fCounterPoolBackground);
  if(!fArraytrack) return taggedphotonic;
  //printf("process track %d\n",iTrack1);
  
  TVector3 v3D1;
  TVector3 v3D2;

  Double_t valuensparseDeltaPhiMaps[5];
  Double_t valueangle[3];

  valuensparseDeltaPhiMaps[1] = binct;
  valuensparseDeltaPhiMaps[2] = track1->Pt();
  valuensparseDeltaPhiMaps[0] = deltaphi;
  valuensparseDeltaPhiMaps[4] = source;
  
  valueangle[2] = source;
  valueangle[1] = binct;

  //Magnetic Field
  Double_t bfield = vEvent->GetMagneticField();

  // Get Primary vertex
  const AliVVertex *pVtx = vEvent->GetPrimaryVertex();
  
  for(Int_t idex = 0; idex < fCounterPoolBackground; idex++) 
    {

      Int_t iTrack2 = fArraytrack->At(idex);
      //printf("track %d\n",iTrack2);
      AliVTrack* track2 = (AliVTrack *) vEvent->GetTrack(iTrack2);
      if (!track2) 
	{
	  printf("ERROR: Could not receive track %d\n", iTrack2);
	  continue;
	}
      if(iTrack2==iTrack1) continue;
      //printf("Different\n");

      // track cuts and PID already done

      // if MC look
      if(mcEvent) {
	Int_t source2 = 0;
	Int_t indexmother2 = -1;
	source2 = FindMother(TMath::Abs(track2->GetLabel()),mcEvent, indexmother2);
	if(source2 >=0 ) {
	  if((indexmother2 == indexmother) && (source == source2)) {
	    if(source == kElectronfromconversion) {
	      valueangle[2] = kElectronfromconversionboth;
	      valuensparseDeltaPhiMaps[4] = kElectronfromconversionboth;
	    }
	    if(source == kElectronfrompi0) {
	      valueangle[2] = kElectronfrompi0both;
	      valuensparseDeltaPhiMaps[4] = kElectronfrompi0both;
	    }
	    if(source == kElectronfrometa) {
	      valueangle[2] = kElectronfrometaboth;
	      valuensparseDeltaPhiMaps[4] = kElectronfrometaboth;
	    }
	  }
	}
      }
      
      if(fAlgorithmMA && (!fAODAnalysis))
	{
	  // tracks
	  AliESDtrack *esdtrack2 = dynamic_cast<AliESDtrack *>(track2);   
	  AliESDtrack *esdtrack1 = dynamic_cast<AliESDtrack *>(track1);      
	  if((!esdtrack2) || (!esdtrack1)) continue;

	  //Variables
	  Double_t p1[3];
	  Double_t p2[3];
	  Double_t xt1; //radial position track 1 at the DCA point
	  Double_t xt2; //radial position track 2 at the DCA point
	  //DCA track1-track2
	  Double_t dca12 = esdtrack2->GetDCA(esdtrack1,bfield,xt2,xt1);

	  // Cut dca
	  if(dca12 > fMaxdca) continue;
	  
	  //Momento of the track extrapolated to DCA track-track	
	  //Track1
	  Bool_t hasdcaT1 = esdtrack1->GetPxPyPzAt(xt1,bfield,p1);
	  //Track2
	  Bool_t hasdcaT2 = esdtrack2->GetPxPyPzAt(xt2,bfield,p2);
	  
	  if(!hasdcaT1 || !hasdcaT2) AliWarning("It could be a problem in the extrapolation");
	  
	  //track1-track2 Invariant Mass
	  Double_t eMass = 0.000510998910; //Electron mass in GeV
	  Double_t pP1 = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]); //Track 1 momentum
	  Double_t pP2 = sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]); //Track 2 momentum
	  Double_t eE1 = TMath::Sqrt(pP1*pP1+eMass*eMass);
	  Double_t eE2 = TMath::Sqrt(pP2*pP2+eMass*eMass);
	  
	  //TLorentzVector v1(p1[0],p1[1],p1[2],sqrt(eMass*eMass+pP1*pP1));
	  //TLorentzVector v2(p2[0],p2[1],p2[2],sqrt(eMass*eMass+pP2*pP2));
	  //Double_t imass = (v1+v2).M(); //Invariant Mass
	  //Double_t angle3D = v1.Angle(v2.Vect()); //Opening Angle (Total Angle)
	  
	  // daughter
	  v3D1.SetXYZ(p1[0],p1[1],p1[2]);
	  v3D2.SetXYZ(p2[0],p2[1],p2[2]);
	  Double_t openingangle = TVector2::Phi_0_2pi(v3D2.Angle(v3D1));
	  
	  // mother
	  TVector3 motherrec = v3D1 + v3D2;
	  Double_t invmass = TMath::Sqrt((eE1+eE2)*(eE1+eE2)-(motherrec.Px()*motherrec.Px()+motherrec.Py()*motherrec.Py()+motherrec.Pz()*motherrec.Pz()));
	  
	  // xy
	  //TVector3 vectordiff = v3D1 - v3D2;
	  //Double_t diffphi = TVector2::Phi_0_2pi(vectordiff.Phi());
	  //Double_t massxy = TMath::Sqrt((eE1+eE2)*(eE1+eE2)-(pP1*pP1+pP2*pP2+2*pP1*pP2*TMath::Cos(diffphi)));

	  // rz
	  //Double_t difftheta = TVector2::Phi_0_2pi(vectordiff.Eta());
	  //Double_t massrz = TMath::Sqrt((eE1+eE2)*(eE1+eE2)-(pP1*pP1+pP2*pP2+2*pP1*pP2*TMath::Cos(difftheta)));
  

	  Float_t fCharge1 = track1->Charge();
	  Float_t fCharge2 = track2->Charge();

	  // Fill Histo
	  //valueangle[0] = diffphi;
	  //valueangle[1] = difftheta;
	  valueangle[0] = openingangle;
	  if((fCharge1*fCharge2)>0.0) fSameSignAngle->Fill(&valueangle[0]);
	  else fOppSignAngle->Fill(&valueangle[0]);

	  // Cut
	  if(openingangle > fMaxopening3D) continue;
	  //if(difftheta > fMaxopeningtheta) continue;
	  //if(diffphi > fMaxopeningphi) continue;

	  // Invmass
	  valuensparseDeltaPhiMaps[3] = invmass;
	  if((fCharge1*fCharge2)>0.0) fSameSignDeltaPhiMaps->Fill(&valuensparseDeltaPhiMaps[0]);
	  else fOppSignDeltaPhiMaps->Fill(&valuensparseDeltaPhiMaps[0]);
	  
	  // Cut
	  if(invmass < fMaxInvmass) {
	    if((fCharge1*fCharge2)<0.0) oppositetaggedphotonic=kTRUE;
	    if((fCharge1*fCharge2)>0.0) sametaggedphotonic=kTRUE;
	  }


	}
      else 
	{
	  Int_t fPDGtrack1 = 11; 
	  Int_t fPDGtrack2 = 11;
	  
	  Float_t fCharge1 = track1->Charge();
	  Float_t fCharge2 = track2->Charge();
	  
	  if(fCharge1>0) fPDGtrack1 = -11;
	  if(fCharge2>0) fPDGtrack2 = -11;
	  
	  AliKFParticle ktrack1(*track1, fPDGtrack1);
	  AliKFParticle ktrack2(*track2, fPDGtrack2);
	  AliKFParticle recoGamma(ktrack1, ktrack2);
	  
	  //Reconstruction Cuts
	  if(recoGamma.GetNDF()<1) continue;
	  Double_t chi2OverNDF = recoGamma.GetChi2()/recoGamma.GetNDF();
	  if(TMath::Sqrt(TMath::Abs(chi2OverNDF))>fChi2OverNDFCut) continue;
	  
	  // if set mass constraint
	  if(fSetMassConstraint && pVtx) {
	    AliKFVertex primV(*pVtx);
	    primV += recoGamma;
	    recoGamma.SetProductionVertex(primV);
	    recoGamma.SetMassConstraint(0,0.0001);
	  }    

	  //Invariant Mass
	  Double_t imass; 
	  Double_t width;
	  recoGamma.GetMass(imass,width);
	  
	  //Opening Angle (Total Angle)
	  Double_t angle = ktrack1.GetAngle(ktrack2);
	  valueangle[0] = angle;
	  if((fCharge1*fCharge2)>0.0) fSameSignAngle->Fill(&valueangle[0]);
	  else fOppSignAngle->Fill(&valueangle[0]);

	  // Cut
	  if(angle > fMaxopening3D) continue;	  

	  // Invmass
	  valuensparseDeltaPhiMaps[3] = imass;
	  if((fCharge1*fCharge2)>0.0) fSameSignDeltaPhiMaps->Fill(&valuensparseDeltaPhiMaps[0]);
	  else fOppSignDeltaPhiMaps->Fill(&valuensparseDeltaPhiMaps[0]);
	  
	  
	  // Cut
	  if(imass < fMaxInvmass) {
	    if((fCharge1*fCharge2)<0.0) oppositetaggedphotonic=kTRUE;
	    if((fCharge1*fCharge2)>0.0) sametaggedphotonic=kTRUE;
	  }
	
	}
    }
  
  if(oppositetaggedphotonic && sametaggedphotonic){
    taggedphotonic = 6;
  }

  if(!oppositetaggedphotonic && sametaggedphotonic){
    taggedphotonic = 4;
  }

  if(oppositetaggedphotonic && !sametaggedphotonic){
    taggedphotonic = 2;
  }

  
  return taggedphotonic;
}

//_________________________________________________________________________
Int_t AliAnalysisTaskHFEFlow::FindMother(Int_t tr, AliMCEvent *mcEvent, Int_t &indexmother){
  //
  // Find the mother if MC
  //

  if(!mcEvent) return 0;

  Int_t pdg = CheckPdg(tr,mcEvent);
  if(TMath::Abs(pdg)!= 11) {
    indexmother = -1;
    return kNoElectron;
  }
  
  indexmother = IsMotherGamma(tr,mcEvent);
  if(indexmother > 0) return kElectronfromconversion;
  indexmother = IsMotherPi0(tr,mcEvent);
  if(indexmother > 0) return kElectronfrompi0;
  indexmother = IsMotherC(tr,mcEvent);
  if(indexmother > 0) return kElectronfromC;
  indexmother = IsMotherB(tr,mcEvent);
  if(indexmother > 0) return kElectronfromB;
  indexmother = IsMotherEta(tr,mcEvent);
  if(indexmother > 0) return kElectronfrometa;
  
  return kElectronfromother;


}
//____________________________________________________________________________________________________________
Int_t AliAnalysisTaskHFEFlow::CheckPdg(Int_t tr, AliMCEvent* mcEvent) {

  //
  // Return the pdg of the particle
  //

  Int_t pdgcode = -1;

  if(tr < 0) return pdgcode;
  AliVParticle *mctrack = mcEvent->GetTrack(tr);
 
  
  if(mctrack->IsA() == AliMCParticle::Class()) {
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return pdgcode;
    pdgcode = mctrackesd->PdgCode();
  }

  if(mctrack->IsA() == AliAODMCParticle::Class()) {
    AliAODMCParticle *mctrackaod = NULL;
    if(!(mctrackaod = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return pdgcode;
    pdgcode = mctrackaod->GetPdgCode();
  }
  
  return pdgcode;

 
}
//____________________________________________________________________________________________________________
Int_t AliAnalysisTaskHFEFlow::IsMotherGamma(Int_t tr, AliMCEvent* mcEvent) {

  //
  // Return the lab of gamma mother or -1 if not gamma
  //

  if(tr < 0) return -1;
  AliVParticle *mctrack = mcEvent->GetTrack(tr);
  
  if(mctrack->IsA() == AliMCParticle::Class()) {
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();
    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother(); 
    if(imother < 0) return -1;  
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;
    // Check gamma    
    Int_t pdg = mother->GetPdgCode();
    if(TMath::Abs(pdg) == 22) return imother;
    if(TMath::Abs(pdg) == 11) {
      return IsMotherGamma(imother,mcEvent);
    }
    return -1;
  }

  if(mctrack->IsA() == AliAODMCParticle::Class()) {
    AliAODMCParticle *mctrackaod = NULL;
    if(!(mctrackaod = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0) return -1;  
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    // Check gamma    
    Int_t pdg = mothertrack->GetPdgCode();
    if(TMath::Abs(pdg) == 22) return imother;
    if(TMath::Abs(pdg) == 11) {
      return IsMotherGamma(imother,mcEvent);
    }
    return -1;

  }
  
  return -1;

 
}
//
//____________________________________________________________________________________________________________
Int_t AliAnalysisTaskHFEFlow::IsMotherPi0(Int_t tr, AliMCEvent* mcEvent) {

  //
  // Return the lab of pi0 mother or -1 if not pi0
  //

  if(tr < 0) return -1;
  AliVParticle *mctrack = mcEvent->GetTrack(tr);
  
  if(mctrack->IsA() == AliMCParticle::Class()) {
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();
    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother(); 
    if(imother < 0) return -1;  
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;
    // Check gamma    
    Int_t pdg = mother->GetPdgCode();
    if(TMath::Abs(pdg) == 111) return imother;
    if(TMath::Abs(pdg) == 11) {
      return IsMotherPi0(imother,mcEvent);
    }
    return -1;
  }

  if(mctrack->IsA() == AliAODMCParticle::Class()) {
    AliAODMCParticle *mctrackaod = NULL;
    if(!(mctrackaod = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0) return -1;  
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    // Check gamma    
    Int_t pdg = mothertrack->GetPdgCode();
    if(TMath::Abs(pdg) == 111) return imother;
    if(TMath::Abs(pdg) == 11) {
      return IsMotherPi0(imother,mcEvent);
    }
    return -1;
  }

  return -1;
 
}
//____________________________________________________________________________________________________________
Int_t AliAnalysisTaskHFEFlow::IsMotherC(Int_t tr, AliMCEvent* mcEvent) {

  //
  // Return the lab of signal mother or -1 if not signal
  //

  if(tr < 0) return -1;
  AliVParticle *mctrack = mcEvent->GetTrack(tr);
  
  if(mctrack->IsA() == AliMCParticle::Class()) {
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();
    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother(); 
    if(imother < 0) return -1;  
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;
    // Check gamma    
    Int_t pdg = mother->GetPdgCode();
    if((TMath::Abs(pdg)==411) || (TMath::Abs(pdg)==421) || (TMath::Abs(pdg)==431) || (TMath::Abs(pdg)==4122) || (TMath::Abs(pdg)==4132) || (TMath::Abs(pdg)==4232) || (TMath::Abs(pdg)==43320)) return imother;
    if(TMath::Abs(pdg) == 11) {
      return IsMotherC(imother,mcEvent);
    }
    return -1;
  }

  if(mctrack->IsA() == AliAODMCParticle::Class()) {
    AliAODMCParticle *mctrackaod = NULL;
    if(!(mctrackaod = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0) return -1;  
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    // Check gamma    
    Int_t pdg = mothertrack->GetPdgCode();
    if((TMath::Abs(pdg)==411) || (TMath::Abs(pdg)==421) || (TMath::Abs(pdg)==431) || (TMath::Abs(pdg)==4122) || (TMath::Abs(pdg)==4132) || (TMath::Abs(pdg)==4232) || (TMath::Abs(pdg)==43320)) return imother;
    if(TMath::Abs(pdg) == 11) {
      return IsMotherC(imother,mcEvent);
    }
    return -1;
  }

  return -1;

}
//____________________________________________________________________________________________________________
Int_t AliAnalysisTaskHFEFlow::IsMotherB(Int_t tr, AliMCEvent* mcEvent) {

  //
  // Return the lab of signal mother or -1 if not signal
  //

  if(tr < 0) return -1;
  AliVParticle *mctrack = mcEvent->GetTrack(tr);
  
  if(mctrack->IsA() == AliMCParticle::Class()) {
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();
    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother(); 
    if(imother < 0) return -1;  
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;
    // Check gamma    
    Int_t pdg = mother->GetPdgCode();
    if((TMath::Abs(pdg)==511) || (TMath::Abs(pdg)==521) || (TMath::Abs(pdg)==531) || (TMath::Abs(pdg)==5122) || (TMath::Abs(pdg)==5132) || (TMath::Abs(pdg)==5232) || (TMath::Abs(pdg)==53320)) return imother; 
    if(TMath::Abs(pdg) == 11) {
      return IsMotherB(imother,mcEvent);
    }
    return -1;
  }

  if(mctrack->IsA() == AliAODMCParticle::Class()) {
    AliAODMCParticle *mctrackaod = NULL;
    if(!(mctrackaod = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0) return -1;  
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    // Check gamma    
    Int_t pdg = mothertrack->GetPdgCode();
    if((TMath::Abs(pdg)==511) || (TMath::Abs(pdg)==521) || (TMath::Abs(pdg)==531) || (TMath::Abs(pdg)==5122) || (TMath::Abs(pdg)==5132) || (TMath::Abs(pdg)==5232) || (TMath::Abs(pdg)==53320)) return imother;
    if(TMath::Abs(pdg) == 11) {
      return IsMotherB(imother,mcEvent);
    }
    return -1;
  }

  return -1;

}
//____________________________________________________________________________________________________________
Int_t AliAnalysisTaskHFEFlow::IsMotherEta(Int_t tr, AliMCEvent* mcEvent) {

  //
  // Return the lab of pi0 mother or -1 if not pi0
  //

 if(tr < 0) return -1;
  AliVParticle *mctrack = mcEvent->GetTrack(tr);
  
  if(mctrack->IsA() == AliMCParticle::Class()) {
    AliMCParticle *mctrackesd = NULL;
    if(!(mctrackesd = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    TParticle *particle = 0x0;
    particle = mctrackesd->Particle();
    // Take mother
    if(!particle) return -1;
    Int_t imother   = particle->GetFirstMother(); 
    if(imother < 0) return -1;  
    AliMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    TParticle * mother = mothertrack->Particle();
    if(!mother) return -1;
    // Check gamma    
    Int_t pdg = mother->GetPdgCode();
    if(TMath::Abs(pdg) == 221) return imother;
    if(TMath::Abs(pdg) == 11) {
      return IsMotherEta(imother,mcEvent);
    }
    return -1;
  }

  if(mctrack->IsA() == AliAODMCParticle::Class()) {
    AliAODMCParticle *mctrackaod = NULL;
    if(!(mctrackaod = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(tr))))) return -1;
    // Take mother
    Int_t imother = mctrackaod->GetMother();
    if(imother < 0) return -1;  
    AliAODMCParticle *mothertrack = NULL;
    if(!(mothertrack = dynamic_cast<AliAODMCParticle *>(mcEvent->GetTrack(TMath::Abs(imother))))) return -1;
    // Check gamma    
    Int_t pdg = mothertrack->GetPdgCode();
    if(TMath::Abs(pdg) == 221) return imother;
    if(TMath::Abs(pdg) == 11) {
      return IsMotherEta(imother,mcEvent);
    }
    return -1;
  }

  return -1;
  
}
