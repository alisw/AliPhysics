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
//   Theodor Rascanu <trascanu@stud.uni-frankfurt.de>
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
#include "TF1.h"
#include <TArrayD.h>

#include <TDirectory.h>
#include <TTreeStream.h>

#include "AliVEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliPID.h"
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
#include "AliGenCocktailEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliHFEcollection.h"

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
#include "AliAnalysisTaskFlowTPCTOFEPSP.h"
#include "AliAODMCHeader.h"
#include "TClonesArray.h"
#include "AliHFENonPhotonicElectron.h"
#include "AliHFEmcQA.h"

ClassImp(AliAnalysisTaskFlowTPCTOFEPSP)

//____________________________________________________________________
AliAnalysisTaskFlowTPCTOFEPSP::AliAnalysisTaskFlowTPCTOFEPSP() :
AliAnalysisTaskSE(),
  fListHist(0x0),
  fHistMCQA(0x0),
  fAODAnalysis(kFALSE),
  fFilter(1<<4),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fBackgroundSubtraction(NULL),
  fMCQA(NULL),
  fSelectGenerator(-1),
  fVZEROEventPlane(kFALSE),
  fVZEROEventPlaneA(kFALSE),
  fVZEROEventPlaneC(kFALSE),
  fSubEtaGapTPC(kFALSE),
  fEtaGap(0.0),
  fPtBinning(),
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
  fSP(kFALSE),
  fVariableMultiplicity(0),
  fTriggerUsed(0),
  fMCPID(kFALSE),
  fNoPID(kFALSE),
  fMonitorEventPlane(kFALSE),
  fMonitorContamination(kFALSE),
  fMonitorPhotonic(kFALSE),
  fMonitorWithoutPID(kFALSE),
  fMonitorTrackCuts(kFALSE),
  fMonitorQCumulant(kFALSE),
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fPID(0),
  fPIDTOFOnly(0),
  fPIDqa(0),
  fflowEvent(NULL),
  fAsFunctionOfP(kTRUE),
  fHFEVZEROEventPlane(0x0),
  fHistEV(0),
  fHistPileUp(0),
  fPileUpCut(kFALSE),
  fEventPlane(0x0),
  fEventPlaneaftersubtraction(0x0),
  fFractionContamination(0x0),
  fContaminationv2(0x0),
  fContaminationmeanpt(0x0),
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
  fDeltaPhiMapsContamination(0x0),
  fCosPhiMaps(0x0),
  fProfileCosPhiMaps(0x0),
  fDebugStreamer(0)
{
  // Constructor

  for(Int_t k = 0; k < 10; k++) {
    fBinCentralityLess[k] = 0.0;
  }
  for(Int_t k = 0; k < 11; k++) {
    fContamination[k] = NULL;
    fv2contamination[k] = NULL;
  }
  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
   
}
//______________________________________________________________________________
AliAnalysisTaskFlowTPCTOFEPSP:: AliAnalysisTaskFlowTPCTOFEPSP(const char *name) :
  AliAnalysisTaskSE(name),
  fListHist(0x0),
  fHistMCQA(0x0),
  fAODAnalysis(kFALSE),
  fFilter(1<<4), 
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fBackgroundSubtraction(NULL),
  fMCQA(NULL),
  fSelectGenerator(-1),
  fVZEROEventPlane(kFALSE),
  fVZEROEventPlaneA(kFALSE),
  fVZEROEventPlaneC(kFALSE),
  fSubEtaGapTPC(kFALSE),
  fEtaGap(0.0),
  fPtBinning(),
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
  fSP(kFALSE),
  fVariableMultiplicity(0),
  fTriggerUsed(0),  
  fMCPID(kFALSE),
  fNoPID(kFALSE),
  fMonitorEventPlane(kFALSE),
  fMonitorContamination(kFALSE),
  fMonitorPhotonic(kFALSE),
  fMonitorWithoutPID(kFALSE),
  fMonitorTrackCuts(kFALSE),
  fMonitorQCumulant(kFALSE),
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fPID(0),
  fPIDTOFOnly(0),
  fPIDqa(0),
  fflowEvent(NULL),
  fAsFunctionOfP(kTRUE),
  fHFEVZEROEventPlane(0x0),
  fHistEV(0),
  fHistPileUp(0),
  fPileUpCut(kFALSE),
  fEventPlane(0x0),
  fEventPlaneaftersubtraction(0x0),
  fFractionContamination(0x0),
  fContaminationv2(0x0),
  fContaminationmeanpt(0x0),
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
  fDeltaPhiMapsContamination(0x0),
  fCosPhiMaps(0x0),
  fProfileCosPhiMaps(0x0),
  fDebugStreamer(0)
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

  for(Int_t k = 0; k < 11; k++) {
    fContamination[k] = NULL;
    fv2contamination[k] = NULL;
  }
  
  fPID = new AliHFEpid("hfePid\n");
  fPIDqa = new AliHFEpidQAmanager;

  fPIDTOFOnly = new AliHFEpid("hfePidTOFOnly\n");

  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));

  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
  //for(Int_t bincless = 0; bincless < fNbBinsCentralityQCumulant; bincless++) {
  //  DefineOutput(bincless+2,AliFlowEventSimple::Class()); 
  //}
  
}
//____________________________________________________________
AliAnalysisTaskFlowTPCTOFEPSP::AliAnalysisTaskFlowTPCTOFEPSP(const AliAnalysisTaskFlowTPCTOFEPSP &ref):
  AliAnalysisTaskSE(ref),
  fListHist(NULL),
  fHistMCQA(0x0),
  fAODAnalysis(ref.fAODAnalysis), 
  fFilter(ref.fFilter),
  fAODMCHeader(ref.fAODMCHeader),
  fAODArrayMCInfo(ref.fAODArrayMCInfo),
  fBackgroundSubtraction(ref.fBackgroundSubtraction),
  fMCQA(NULL),
  fSelectGenerator(ref.fSelectGenerator),
  fVZEROEventPlane(ref.fVZEROEventPlane),
  fVZEROEventPlaneA(ref.fVZEROEventPlaneA),
  fVZEROEventPlaneC(ref.fVZEROEventPlaneC),
  fSubEtaGapTPC(ref.fSubEtaGapTPC),
  fEtaGap(ref.fEtaGap),
  fPtBinning(ref.fPtBinning),
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
  fSP(ref.fSP),
  fVariableMultiplicity(ref.fVariableMultiplicity),
  fTriggerUsed(ref.fTriggerUsed),
  fMCPID(ref.fMCPID),
  fNoPID(ref.fNoPID),
  fMonitorEventPlane(ref.fMonitorEventPlane),
  fMonitorContamination(ref.fMonitorContamination),
  fMonitorPhotonic(ref.fMonitorPhotonic),
  fMonitorWithoutPID(ref.fMonitorWithoutPID),
  fMonitorTrackCuts(ref.fMonitorTrackCuts),
  fMonitorQCumulant(ref.fMonitorQCumulant),
  fcutsRP(NULL),
  fcutsPOI(NULL),
  fHFECuts(NULL),
  fPID(NULL),
  fPIDTOFOnly(NULL),
  fPIDqa(NULL),
  fflowEvent(NULL),
  fAsFunctionOfP(ref.fAsFunctionOfP),
  fHFEVZEROEventPlane(NULL),
  fHistEV(NULL),
  fHistPileUp(NULL),
  fPileUpCut(kFALSE),
  fEventPlane(NULL),
  fEventPlaneaftersubtraction(NULL),
  fFractionContamination(NULL),
  fContaminationv2(NULL),
  fContaminationmeanpt(0x0),
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
  fDeltaPhiMapsContamination(NULL),
  fCosPhiMaps(NULL),
  fProfileCosPhiMaps(NULL),
  fDebugStreamer(0)
{
  //
  // Copy Constructor
  //

  for(Int_t k = 0; k < 10; k++) {
    fBinCentralityLess[k] = 0.0;
  }
  for(Int_t k = 0; k < 11; k++) {
    fContamination[k] = NULL;
    fv2contamination[k] = NULL;
  }
   
  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
  
  ref.Copy(*this);
}

//____________________________________________________________
AliAnalysisTaskFlowTPCTOFEPSP &AliAnalysisTaskFlowTPCTOFEPSP::operator=(const AliAnalysisTaskFlowTPCTOFEPSP &ref){
  //
  // Assignment operator
  //
  if(this == &ref) 
    ref.Copy(*this);
  return *this;
}

//____________________________________________________________
void AliAnalysisTaskFlowTPCTOFEPSP::Copy(TObject &o) const {
  // 
  // Copy into object o
  //
  AliAnalysisTaskFlowTPCTOFEPSP &target = dynamic_cast<AliAnalysisTaskFlowTPCTOFEPSP &>(o);
  target.fListHist = fListHist;
  target.fHistMCQA = fHistMCQA;
  target.fAODAnalysis = fAODAnalysis;
  target.fFilter = fFilter;
  target.fAODMCHeader = fAODMCHeader;
  target.fAODArrayMCInfo = fAODArrayMCInfo;
  target.fBackgroundSubtraction = fBackgroundSubtraction;
  target.fMCQA = fMCQA;
  target.fSelectGenerator = fSelectGenerator;
  target.fVZEROEventPlane = fVZEROEventPlane;
  target.fVZEROEventPlaneA = fVZEROEventPlaneA;
  target.fVZEROEventPlaneC = fVZEROEventPlaneC;
  target.fSubEtaGapTPC = fSubEtaGapTPC;
  target.fEtaGap = fEtaGap;
  target.fPtBinning = fPtBinning;
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
  target.fSP = fSP;
  target.fVariableMultiplicity = fVariableMultiplicity;
  target.fTriggerUsed = fTriggerUsed;
  target.fMCPID = fMCPID;
  target.fNoPID = fNoPID;
  target.fMonitorEventPlane = fMonitorEventPlane;
  target.fMonitorContamination = fMonitorContamination;
  target.fMonitorPhotonic = fMonitorPhotonic;
  target.fMonitorWithoutPID = fMonitorWithoutPID;
  target.fMonitorTrackCuts = fMonitorTrackCuts;
  target.fMonitorQCumulant = fMonitorQCumulant;
  target.fcutsRP = fcutsRP;
  target.fcutsPOI = fcutsPOI;
  target.fHFECuts = fHFECuts;
  target.fPID = fPID;
  target.fPIDTOFOnly = fPIDTOFOnly;
  target.fPIDqa = fPIDqa;
  target.fflowEvent = fflowEvent;
  target.fAsFunctionOfP = fAsFunctionOfP;	 	 	 
  target.fHFEVZEROEventPlane = fHFEVZEROEventPlane; 	
  target.fHistEV=fHistEV;       		 	 
  target.fHistPileUp=fHistPileUp;   		 	 
  target.fPileUpCut=fPileUpCut;         		 	 
  target.fEventPlane=fEventPlane;     		 	 
  target.fEventPlaneaftersubtraction=fEventPlaneaftersubtraction; 		 	 
  target.fFractionContamination=fFractionContamination;     		 	 
  target.fContaminationv2=fContaminationv2;           
  target.fContaminationmeanpt=fContaminationmeanpt;           		 	 
  target.fCosSin2phiep=fCosSin2phiep;         		 	 
  target.fCos2phie=fCos2phie;   		 	 
  target.fSin2phie=fSin2phie;   		 	 
  target.fCos2phiep=fCos2phiep; 		 	 
  target.fSin2phiep=fSin2phiep; 		 	 
  target.fSin2phiephiep=fSin2phiephiep;   		 	 
  target.fCosResabc=fCosResabc; 		 	 
  target.fSinResabc=fSinResabc; 		 	 
  target.fProfileCosResab=fProfileCosResab; 		 	 
  target.fProfileCosResac=fProfileCosResac; 		 	 
  target.fProfileCosResbc=fProfileCosResbc; 		 	 
  target.fCosRes=fCosRes; 		 	 
  target.fSinRes=fSinRes; 		 	 
  target.fProfileCosRes=fProfileCosRes; 		 	 
  target.fTrackingCuts=fTrackingCuts; 		 	 
  target.fDeltaPhiMapsBeforePID=fDeltaPhiMapsBeforePID; 		 	 
  target.fCosPhiMapsBeforePID=fCosPhiMapsBeforePID; 		 	 
  target.fDeltaPhiMaps=fDeltaPhiMaps; 		 	 
  target.fDeltaPhiMapsContamination=fDeltaPhiMapsContamination; 		 	 
  target.fCosPhiMaps=fCosPhiMaps;         		 	 
  target.fProfileCosPhiMaps=fProfileCosPhiMaps;   		 	 
  
  for(Int_t k = 0; k < 10; k++) {
    target.fBinCentralityLess[k] = fBinCentralityLess[k];
  }
  for(Int_t k = 0; k < 11; k++) {
    target.fContamination[k] = fContamination[k];
    target.fv2contamination[k] = fv2contamination[k];
  }
  target.fDebugStreamer=fDebugStreamer;
}
//____________________________________________________________
AliAnalysisTaskFlowTPCTOFEPSP::~AliAnalysisTaskFlowTPCTOFEPSP(){
  //
  // Destructor
  //
  

  if(fListHist) delete fListHist;
  if(fcutsRP) delete fcutsRP;
  if(fcutsPOI) delete fcutsPOI;
  if(fHFECuts) delete fHFECuts;
  if(fPID) delete fPID;
  if(fPIDTOFOnly) delete fPIDTOFOnly;
  //if(fPIDqa) delete fPIDqa;
  if(fflowEvent) delete fflowEvent;
  if ( fDebugStreamer ) delete fDebugStreamer;
  if(fMCQA) delete fMCQA;
  

}
//________________________________________________________________________
void AliAnalysisTaskFlowTPCTOFEPSP::UserCreateOutputObjects()
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

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: User create output objects\n");


  // AOD or ESD
  AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
    SetAODAnalysis(kTRUE);
    //printf("Put AOD analysis on\n");
  } else {
    SetAODAnalysis(kFALSE);
  }

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: AOD ESD\n");

  // RP TRACK CUTS:
  fcutsRP = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();
  fcutsRP->SetName("StandartTPC\n");
  fcutsRP->SetEtaRange(-0.9,0.9);
  fcutsRP->SetQA(kTRUE);
  //TList *qaCutsRP = fcutsRP->GetQA();
  //qaCutsRP->SetName("QA_StandartTPC_RP\n");

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: cutsRP\n");

  //POI TRACK CUTS:
  fcutsPOI = new AliFlowTrackCuts("dummy\n");
  fcutsPOI->SetParamType(AliFlowTrackCuts::kGlobal);
  fcutsPOI->SetPtRange(+1,-1); // select nothing QUICK
  fcutsPOI->SetEtaRange(+1,-1); // select nothing VZERO

  if( fflowEvent ){ 
    fflowEvent->~AliFlowEvent();
    new(fflowEvent) AliFlowEvent(fcutsRP,fcutsPOI);
  }
  else fflowEvent = new AliFlowEvent(fcutsRP,fcutsPOI);
    
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: cutsPOI\n");

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

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: common constants\n");

  // HFE cuts

  if(!fHFECuts){
    fHFECuts = new AliHFEcuts;
    fHFECuts->CreateStandardCuts();
  }
  fHFECuts->Initialize();
  if(fAODAnalysis) fHFECuts->SetAOD();  

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: HFE cuts\n");


  // PID HFE
  //fPID->SetHasMCData(HasMCData());
  if(!fPID) {
    fPID =new AliHFEpid("hfePid\n");
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: pid init 0\n");
  }
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: pid init 1\n");
  if(!fPID->GetNumberOfPIDdetectors()) fPID->AddDetector("TPC", 0);
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: GetNumber of PID detectors %d\n",fPID->GetNumberOfPIDdetectors());
  fPID->InitializePID();
  fPIDqa->Initialize(fPID);
  fPID->SortDetectors();
  
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: pid and pidqa\n");

  if(!fPIDTOFOnly->GetNumberOfPIDdetectors()) {
    fPIDTOFOnly->AddDetector("TOF", 0);
    fPIDTOFOnly->ConfigureTOF(3.);
  }
  fPIDTOFOnly->InitializePID();
  fPIDTOFOnly->SortDetectors();

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: pidtof\n");

 
  //**************************
  // Bins for the THnSparse
  //**************************

  /*
  Int_t nBinsPt = 44;
  Double_t minPt = 0.1;
  Double_t maxPt = 20.0;
  Double_t binLimLogPt[nBinsPt+1];
  Double_t binLimPt[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++) binLimLogPt[i]=(Double_t)TMath::Log10(minPt) + (TMath::Log10(maxPt)-TMath::Log10(minPt))/nBinsPt*(Double_t)i ;
  for(Int_t i=0; i<=nBinsPt; i++) binLimPt[i]=(Double_t)TMath::Power(10,binLimLogPt[i]);
  */

  Int_t nBinsPt = 20;
  Double_t binLimPt[21] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2,
			   1.3, 1.4, 1.5, 2., 2.5, 3., 4., 6.};

  if(!fPtBinning.GetSize()) fPtBinning.Set(nBinsPt+1, binLimPt);

  //Int_t nBinsPtPlus = fNbBinsPtQCumulant;
  //Double_t minPtPlus = fMinPtQCumulant;
  //Double_t maxPtPlus = fMaxPtQCumulant;
  //Double_t binLimPtPlus[nBinsPtPlus+1];
  //for(Int_t i=0; i<=nBinsPtPlus; i++) binLimPtPlus[i]=(Double_t)minPtPlus + (maxPtPlus-minPtPlus)/nBinsPtPlus*(Double_t)i ;

  Int_t nBinsEta = 8;
  Double_t minEta = -0.8;
  Double_t maxEta = 0.8;
  Double_t binLimEta[nBinsEta+1];
  for(Int_t i=0; i<=nBinsEta; i++) binLimEta[i]=(Double_t)minEta + (maxEta-minEta)/nBinsEta*(Double_t)i ;

  Int_t nBinsStep = 6;
  Double_t minStep = 0.;
  Double_t maxStep = 6.;
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

  // Int_t nBinsCosSP = 50;
  // Double_t minCosSP = -100.0;
  // Double_t maxCosSP = 100.0;
  // Double_t binLimCosSP[nBinsCosSP+1];
  // for(Int_t i=0; i<=nBinsCosSP; i++) binLimCosSP[i]=(Double_t)minCosSP + (maxCosSP-minCosSP)/nBinsCosSP*(Double_t)i ;
 
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

  Int_t nBinsCEvenMore = 100;
  Double_t minCEvenMore = 0.0;
  Double_t maxCEvenMore = 100.0;
  Double_t binLimCEvenMore[nBinsCEvenMore+1];
  for(Int_t i=0; i<=nBinsCEvenMore; i++) binLimCEvenMore[i]=(Double_t)minCEvenMore + (maxCEvenMore-minCEvenMore)/nBinsCEvenMore*(Double_t)i ;

  Int_t nBinsPhi = 8;
  Double_t minPhi = 0.0;
  Double_t maxPhi = TMath::Pi();
  Double_t binLimPhi[nBinsPhi+1];
  for(Int_t i=0; i<=nBinsPhi; i++) {
    binLimPhi[i]=(Double_t)minPhi + (maxPhi-minPhi)/nBinsPhi*(Double_t)i ;
    //printf("bin phi is %f for %d\n",binLimPhi[i],i);
  }

  Int_t nBinsPhiLess = 2.0;
  Double_t minPhiLess = 0.0;
  Double_t maxPhiLess = 2.0;
  Double_t binLimPhiLess[nBinsPhiLess+1];
  for(Int_t i=0; i<=nBinsPhiLess; i++) {
    binLimPhiLess[i]=(Double_t)minPhiLess + (maxPhiLess-minPhiLess)/nBinsPhiLess*(Double_t)i ;
  }

  Int_t nBinsTPCdEdx = 140;
  Double_t minTPCdEdx = -12.0;
  Double_t maxTPCdEdx = 12.0;
  Double_t binLimTPCdEdx[nBinsTPCdEdx+1];
  for(Int_t i=0; i<=nBinsTPCdEdx; i++) {
    binLimTPCdEdx[i]=(Double_t)minTPCdEdx + (maxTPCdEdx-minTPCdEdx)/nBinsTPCdEdx*(Double_t)i ;
  }

  Int_t nBinsAngle = 40;
  Double_t minAngle = 0.0;
  Double_t maxAngle = 1.0;
  Double_t binLimAngle[nBinsAngle+1];
  for(Int_t i=0; i<=nBinsAngle; i++) {
    binLimAngle[i]=(Double_t)minAngle + (maxAngle-minAngle)/nBinsAngle*(Double_t)i ;
    //printf("bin phi is %f for %d\n",binLimAngle[i],i);
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

  Int_t nBinsMult = 100;
  Double_t minMult = 0.;
  Double_t maxMult = 3000;
  Double_t binLimMult[nBinsMult+1];
  //for(Int_t i=0; i<=nBinsMult; i++) binLimMult[i]=TMath::Power((Double_t)minMult + (TMath::Sqrt(maxMult)-TMath::Sqrt(minMult))/nBinsMult*(Double_t)i,2);
  for(Int_t i=0; i<=nBinsMult; i++) binLimMult[i]=(Double_t)minMult + (maxMult-minMult)/nBinsMult*(Double_t)i;

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: variables\n");
  
  //******************
  // Histograms
  //******************
    
  fListHist = new TList();
  fListHist->SetOwner();

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: list\n");

  // Minimum histos

  // Histos
  fHistEV = new TH2D("fHistEV", "events", 3, 0, 3, 3, 0,3);
  
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: histev\n");

  // V0 multiplicity vs # of tracks vs centraliy
  const Int_t nDimPU=4;
  Int_t nBinPU[nDimPU] = {nBinsCEvenMore,nBinsCEvenMore,nBinsMult,nBinsMult};
  fHistPileUp = new THnSparseF("PileUp","PileUp",nDimPU,nBinPU);
  fHistPileUp->SetBinEdges(0,binLimCEvenMore);
  fHistPileUp->SetBinEdges(1,binLimCEvenMore);
  fHistPileUp->SetBinEdges(2,binLimMult);
  fHistPileUp->SetBinEdges(3,binLimMult);
  fHistPileUp->Sumw2();
  
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: eventplane\n");
  
  // Event plane as function of phiep, centrality
  const Int_t nDima=4;
  Int_t nBina[nDima] = {nBinsPhi,nBinsPhi,nBinsPhi,nBinsC};
  fEventPlane = new THnSparseF("EventPlane","EventPlane",nDima,nBina);
  fEventPlane->SetBinEdges(0,binLimPhi);
  fEventPlane->SetBinEdges(1,binLimPhi);
  fEventPlane->SetBinEdges(2,binLimPhi);
  fEventPlane->SetBinEdges(3,binLimC);
  fEventPlane->Sumw2();

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: eventplane\n");

  // Fraction of contamination, centrality
  const Int_t nDimcont=2;
  Int_t nBincont[nDimcont] = {fPtBinning.GetSize()-1,nBinsC};
  fFractionContamination = new THnSparseF("Contamination","Contamination",nDimcont,nBincont);
  fFractionContamination->SetBinEdges(0,fPtBinning.GetArray());
  fFractionContamination->SetBinEdges(1,binLimC);
  fFractionContamination->Sumw2();
  //  
  fContaminationv2 = new TProfile2D("Contaminationv2","",nBinsC,binLimC,fPtBinning.GetSize()-1,fPtBinning.GetArray());
  fContaminationv2->Sumw2();
  //  
  fContaminationmeanpt = new TProfile2D("Contaminationmeanpt","",nBinsC,binLimC,fPtBinning.GetSize()-1,fPtBinning.GetArray());
  fContaminationmeanpt->Sumw2();

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: fraction of contamination\n");

  // Resolution cosres_abc centrality
  const Int_t nDimfbis=4;
  Int_t nBinfbis[nDimfbis] = {nBinsCos,nBinsCos,nBinsCos,nBinsCMore};
  fCosResabc = new THnSparseF("CosRes_abc","CosRes_abc",nDimfbis,nBinfbis);
  fCosResabc->SetBinEdges(0,binLimCos);
  fCosResabc->SetBinEdges(1,binLimCos);
  fCosResabc->SetBinEdges(2,binLimCos);
  fCosResabc->SetBinEdges(3,binLimCMore);
  fCosResabc->Sumw2();

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: cosresabc\n");

  // Resolution cosres centrality
  const Int_t nDimf=2;
  Int_t nBinf[nDimf] = {nBinsCos, nBinsCMore};
  fCosRes = new THnSparseF("CosRes","CosRes",nDimf,nBinf);
  fCosRes->SetBinEdges(0,binLimCos);
  fCosRes->SetBinEdges(1,binLimCMore);
  fCosRes->Sumw2();

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: cosres\n");

  // Maps delta phi
  const Int_t nDimg=5;
  Int_t nBing[nDimg] = {nBinsPhi,nBinsC,fPtBinning.GetSize()-1, nBinsCharge,nBinsEtaLess};
  fDeltaPhiMaps = new THnSparseF("DeltaPhiMaps","DeltaPhiMaps",nDimg,nBing);
  fDeltaPhiMaps->SetBinEdges(0,binLimPhi);
  fDeltaPhiMaps->SetBinEdges(1,binLimC);
  fDeltaPhiMaps->SetBinEdges(2,fPtBinning.GetArray());
  fDeltaPhiMaps->SetBinEdges(3,binLimCharge);
  fDeltaPhiMaps->SetBinEdges(4,binLimEtaLess);
  fDeltaPhiMaps->Sumw2();  

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: deltaphimaps\n");

  // Maps cos phi
  const Int_t nDimh=5;
  Int_t nBinh[nDimh] = {nBinsCos,nBinsC,fPtBinning.GetSize()-1,nBinsCharge,nBinsEtaLess};
  fCosPhiMaps = new THnSparseF("CosPhiMaps","CosPhiMaps",nDimh,nBinh);
  fCosPhiMaps->SetBinEdges(0,binLimCos);
  fCosPhiMaps->SetBinEdges(1,binLimC);
  fCosPhiMaps->SetBinEdges(2,fPtBinning.GetArray());
  fCosPhiMaps->SetBinEdges(3,binLimCharge);
  fCosPhiMaps->SetBinEdges(4,binLimEtaLess);
  fCosPhiMaps->Sumw2();

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: cosphimaps\n");

  //
  // fMonitorEventPlane
  //
  //

  if(fMonitorEventPlane) {  
    // Event Plane after subtraction as function of phiep, centrality, pt, eta
    const Int_t nDimb=2;
    Int_t nBinb[nDimb] = {nBinsPhi, nBinsC};
    fEventPlaneaftersubtraction = new THnSparseF("EventPlane_aftersubtraction","EventPlane_aftersubtraction",nDimb,nBinb);
    fEventPlaneaftersubtraction->SetBinEdges(0,binLimPhi);
    fEventPlaneaftersubtraction->SetBinEdges(1,binLimC);
    fEventPlaneaftersubtraction->Sumw2();

    //printf("AliAnalysisTaskFlowTPCTOFEPSP: eventplane after sub\n");
    
    // Monitoring of the event Plane cos(2phi) sin(2phi) centrality
    const Int_t nDimi=3;
    Int_t nBini[nDimi] = {nBinsCos, nBinsCos, nBinsCMore};
    fCosSin2phiep = new THnSparseF("CosSin2phiep","CosSin2phiep",nDimi,nBini);
    fCosSin2phiep->SetBinEdges(0,binLimCos);
    fCosSin2phiep->SetBinEdges(1,binLimCos);
    fCosSin2phiep->SetBinEdges(2,binLimCMore);
    fCosSin2phiep->Sumw2();

    //printf("AliAnalysisTaskFlowTPCTOFEPSP: cossin2phiep\n");
    
    // Monitoring Event plane after subtraction of the track
    const Int_t nDime=4;
    Int_t nBine[nDime] = {nBinsCos, nBinsC, fPtBinning.GetSize()-1, nBinsEta};
    fCos2phie = new THnSparseF("cos2phie","cos2phie",nDime,nBine);
    fCos2phie->SetBinEdges(2,fPtBinning.GetArray());
    fCos2phie->SetBinEdges(3,binLimEta);
    fCos2phie->SetBinEdges(0,binLimCos);
    fCos2phie->SetBinEdges(1,binLimC);
    fCos2phie->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: cos2phie\n");
    fSin2phie = new THnSparseF("sin2phie","sin2phie",nDime,nBine);
    fSin2phie->SetBinEdges(2,fPtBinning.GetArray());
    fSin2phie->SetBinEdges(3,binLimEta);
    fSin2phie->SetBinEdges(0,binLimCos);
    fSin2phie->SetBinEdges(1,binLimC);
    fSin2phie->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: sin2phie\n");
    fCos2phiep = new THnSparseF("cos2phiep","cos2phiep",nDime,nBine);
    fCos2phiep->SetBinEdges(2,fPtBinning.GetArray());
    fCos2phiep->SetBinEdges(3,binLimEta);
    fCos2phiep->SetBinEdges(0,binLimCos);
    fCos2phiep->SetBinEdges(1,binLimC);
    fCos2phiep->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: cos2phiep\n");
    fSin2phiep = new THnSparseF("sin2phiep","sin2phiep",nDime,nBine);
    fSin2phiep->SetBinEdges(2,fPtBinning.GetArray());
    fSin2phiep->SetBinEdges(3,binLimEta);
    fSin2phiep->SetBinEdges(0,binLimCos);
    fSin2phiep->SetBinEdges(1,binLimC);
    fSin2phiep->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: sin2phiep\n");
    fSin2phiephiep = new THnSparseF("sin2phie_phiep","sin2phie_phiep",nDime,nBine);
    fSin2phiephiep->SetBinEdges(2,fPtBinning.GetArray());
    fSin2phiephiep->SetBinEdges(3,binLimEta);
    fSin2phiephiep->SetBinEdges(0,binLimCos);
    fSin2phiephiep->SetBinEdges(1,binLimC);
    fSin2phiephiep->Sumw2();  
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: sin2phiephiep\n");
    
    const Int_t nDimfbiss=4;
    Int_t nBinfbiss[nDimfbiss] = {nBinsCos,nBinsCos,nBinsCos,nBinsC};
    fSinResabc = new THnSparseF("SinRes_abc","SinRes_abc",nDimfbiss,nBinfbiss);
    fSinResabc->SetBinEdges(0,binLimCos);
    fSinResabc->SetBinEdges(1,binLimCos);
    fSinResabc->SetBinEdges(2,binLimCos);
    fSinResabc->SetBinEdges(3,binLimC);
    fSinResabc->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: sinresabc\n");
    
    // Profile cosres centrality with 3 subevents
    fProfileCosResab = new TProfile("ProfileCosRes_a_b","ProfileCosRes_a_b",nBinsCMore,binLimCMore);
    fProfileCosResab->Sumw2();
    fProfileCosResac = new TProfile("ProfileCosRes_a_c","ProfileCosRes_a_c",nBinsCMore,binLimCMore);
    fProfileCosResac->Sumw2();
    fProfileCosResbc = new TProfile("ProfileCosRes_b_c","ProfileCosRes_b_c",nBinsCMore,binLimCMore);
    fProfileCosResbc->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: profilecosresbc\n");
    
    //
    const Int_t nDimff=2;
    Int_t nBinff[nDimff] = {nBinsCos, nBinsC};
    fSinRes = new THnSparseF("SinRes","SinRes",nDimff,nBinff);
    fSinRes->SetBinEdges(0,binLimCos);
    fSinRes->SetBinEdges(1,binLimC);
    fSinRes->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: sinres\n");
    
    // Profile cosres centrality
    fProfileCosRes = new TProfile("ProfileCosRes","ProfileCosRes",nBinsCMore,binLimCMore);
    fProfileCosRes->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: profilecosres\n");
    
    // Profile Maps cos phi
    fProfileCosPhiMaps = new TProfile2D("ProfileCosPhiMaps","ProfileCosPhiMaps",nBinsC,binLimC,fPtBinning.GetSize()-1,fPtBinning.GetArray());
    fProfileCosPhiMaps->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: profilecosphimaps\n");

  }
  //
  // fMonitorTrackCuts
  //

  if(fMonitorTrackCuts) {
    // Debugging tracking steps
    const Int_t nDimTrStep=2;
    Int_t nBinTrStep[nDimTrStep] = {fPtBinning.GetSize()-1,nBinsStep};
    fTrackingCuts = new THnSparseF("TrackingCuts","TrackingCuts",nDimTrStep,nBinTrStep);
    fTrackingCuts->SetBinEdges(0,fPtBinning.GetArray());
    fTrackingCuts->SetBinEdges(1,binLimStep);
    fTrackingCuts->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: trackingcuts\n");
  }

  //
  // fMonitorContamination
  //

  if(fMonitorContamination) { 
    // Maps delta phi contamination
    const Int_t nDimgcont=4;
    Int_t nBingcont[nDimgcont] = {nBinsPhiLess,nBinsC,fPtBinning.GetSize()-1, nBinsTPCdEdx};
    fDeltaPhiMapsContamination = new THnSparseF("DeltaPhiMapsContamination","DeltaPhiMapsContamination",nDimgcont,nBingcont);
    fDeltaPhiMapsContamination->SetBinEdges(0,binLimPhiLess);
    fDeltaPhiMapsContamination->SetBinEdges(1,binLimC);
    fDeltaPhiMapsContamination->SetBinEdges(2,fPtBinning.GetArray());
    fDeltaPhiMapsContamination->SetBinEdges(3,binLimTPCdEdx);
    fDeltaPhiMapsContamination->Sumw2();  
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: deltaphimapscontamination\n");

  }
  //
  // fMonitorWithoutPID
  //

  if(fMonitorWithoutPID) {
    //
    const Int_t nDimgb=3;
    Int_t nBingb[nDimgb] = {nBinsPhi,nBinsC,fPtBinning.GetSize()-1};
    
    fDeltaPhiMapsBeforePID = new THnSparseF("DeltaPhiMapsBeforePID","DeltaPhiMapsBeforePID",nDimgb,nBingb);
    fDeltaPhiMapsBeforePID->SetBinEdges(0,binLimPhi);
    fDeltaPhiMapsBeforePID->SetBinEdges(1,binLimC);
    fDeltaPhiMapsBeforePID->SetBinEdges(2,fPtBinning.GetArray());
    fDeltaPhiMapsBeforePID->Sumw2();  
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: deltaphimapsbeforepid\n");
    
    const Int_t nDimhb=3;
    Int_t nBinhb[nDimhb] = {nBinsCos,nBinsC,fPtBinning.GetSize()-1};
    
    fCosPhiMapsBeforePID = new THnSparseF("CosPhiMapsBeforePID","CosPhiMapsBeforePID",nDimhb,nBinhb);
    fCosPhiMapsBeforePID->SetBinEdges(0,binLimCos);
    fCosPhiMapsBeforePID->SetBinEdges(1,binLimC);
    fCosPhiMapsBeforePID->SetBinEdges(2,fPtBinning.GetArray());
    fCosPhiMapsBeforePID->Sumw2();
    //printf("AliAnalysisTaskFlowTPCTOFEPSP: cosphimapsbeforepid\n");
  }
  //
  // fMonitorPhotonic
  //
  if (fMonitorPhotonic) {
    if(!fBackgroundSubtraction) fBackgroundSubtraction = new AliHFENonPhotonicElectron();
    if(fAODAnalysis) fBackgroundSubtraction->SetAOD(kTRUE);  
    fBackgroundSubtraction->Init();
    // mcQA----------------------------------
    AliInfo("MC QA on");
    if(!fMCQA) fMCQA = new AliHFEmcQA;
    if(!fHistMCQA) fHistMCQA = new TList();
    fHistMCQA->SetOwner();
    fMCQA->SetPbPb();
    //if(TestBit(kWeightHist)){
    //  fMCQA->EnableGetWeightHist();
    //}
    fMCQA->CreatDefaultHistograms(fHistMCQA);
    fMCQA->SetBackgroundWeightFactor(fElecBackgroundFactor[0][0][0],fBinLimit);
    fListHist->Add(fHistMCQA);
  }
  


  //**************************
  // Add to the list
  //******************************

  fListHist->Add(fHistEV);
  fListHist->Add(fHistPileUp);
  fListHist->Add(fEventPlane);
  fListHist->Add(fFractionContamination);
  fListHist->Add(fCosRes);
  fListHist->Add(fCosResabc);
  fListHist->Add(fCosPhiMaps);
  fListHist->Add(fDeltaPhiMaps);
  fListHist->Add(fPIDqa->MakeList("HFEpidQA"));
  fListHist->Add(fContaminationv2);
  fListHist->Add(fContaminationmeanpt);
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: add default\n");

  if(fMonitorEventPlane) {
    fListHist->Add(fProfileCosRes);
    fListHist->Add(fProfileCosResab);
    fListHist->Add(fProfileCosResac);
    fListHist->Add(fProfileCosResbc);
    fListHist->Add(fCosSin2phiep);
    fListHist->Add(fCos2phie);
    fListHist->Add(fSin2phie);
    fListHist->Add(fCos2phiep);
    fListHist->Add(fSin2phiep);
    fListHist->Add(fSin2phiephiep);
    fListHist->Add(fSinRes);
    fListHist->Add(fSinResabc);
    fListHist->Add(fProfileCosPhiMaps);
  }
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: add monitor\n");

  if(fMonitorTrackCuts) fListHist->Add(fTrackingCuts);

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: add monitortrackcuts\n");

  if(fMonitorContamination) {
    fListHist->Add(fDeltaPhiMapsContamination);
  }
  
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: add deltaphimapscontamination\n");

  if(fMonitorWithoutPID) {
    fListHist->Add(fDeltaPhiMapsBeforePID);
    fListHist->Add(fCosPhiMapsBeforePID);
  }

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: add without pid\n");

  if(fMonitorPhotonic) {
    fListHist->Add(fBackgroundSubtraction->GetListOutput());
  }

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: add photonic\n");

  if(fHFEVZEROEventPlane && fMonitorEventPlane) fListHist->Add(fHFEVZEROEventPlane->GetOutputList());
  
  //printf("AliAnalysisTaskFlowTPCTOFEPSP: add event plane\n");

  fListHist->Print();

  PostData(1, fListHist);
 

  //printf("AliAnalysisTaskFlowTPCTOFEPSP: post\n");


}
   
//________________________________________________________________________
void AliAnalysisTaskFlowTPCTOFEPSP::UserExec(Option_t */*option*/)
{
  //
  // Loop over event
  //
   
 
  Double_t massElectron = 0.000511;
  Double_t mcReactionPlane = 0.0;

  Float_t cntr = 0.0;
  Double_t binct = 11.5;
  Double_t binctMore = 20.5;
  Double_t binctLess = -0.5;
  Float_t binctt = -1.0;
  
  Double_t valuecossinephiep[3];
  Double_t valuensparsea[4];
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
  Double_t valuedeltaphicontamination[4];
  Double_t valuefractioncont[2];
   
  AliMCEvent *mcEvent = MCEvent();
  AliMCParticle *mctrack = NULL;
  AliAODMCParticle *mctrackaod = NULL;

  // MC info
  Bool_t mcthere = kTRUE;
  if(fAODAnalysis) {
    AliAODEvent *aodE = dynamic_cast<AliAODEvent *>(fInputEvent);
    if(!aodE){
      //        printf("testd\n\n");
      AliError("No AOD Event\n");
      return;
    }
    fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fInputEvent->FindListObject(AliAODMCHeader::StdBranchName()));
    if(!fAODMCHeader){ 
      mcthere = kFALSE;
    }
    fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fAODArrayMCInfo){ 
      mcthere = kFALSE;
    }
    else {
      fHFECuts->SetMCEvent(aodE);
      if(fMonitorPhotonic) {
	fBackgroundSubtraction->SetAODArrayMCInfo(fAODArrayMCInfo);
	//printf("Set the AOD Array MC Info for BackgroundSubtraction\n");
	fMCQA->SetMCArray(fAODArrayMCInfo);
      }
    }
  }
  else {
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH) mcthere = kFALSE;
    else {
      if(fMonitorPhotonic) {
	fBackgroundSubtraction->SetMCEvent(fMCEvent);
	fMCQA->SetMCEvent(fMCEvent);
	fMCQA->SetGenEventHeader(fMCEvent->GenEventHeader());
      }
    }
  }

  /////////////////////
  // Trigger selection
  ////////////////////

  UInt_t isEventSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(fTriggerUsed==0){
    
    // MB, semi-central and central
    
    if ( !((isEventSelected & AliVEvent::kCentral) |
	   (isEventSelected & AliVEvent::kSemiCentral) |
	   (isEventSelected & AliVEvent::kMB)) ) return;
    
  }
  else if(fTriggerUsed==1){
    
    // semi-central Ionut

    if ( !((isEventSelected & AliVEvent::kCentral) |
	   (isEventSelected & AliVEvent::kSemiCentral) |
	   (isEventSelected & AliVEvent::kMB)) ) return;
    
    Bool_t isMB = (InputEvent()->GetTriggerMask() & (ULong64_t(1)<<1));
    //Bool_t isCentral = (InputEvent()->GetTriggerMask() & (ULong64_t(1)<<4));
    Bool_t isSemiCentral = (InputEvent()->GetTriggerMask() & (ULong64_t(1)<<7));
    
    if(!(isSemiCentral | isMB)) return;
    
  }
  else if(fTriggerUsed==2){

    // semi-central Andrea and Muons
    
    if ( !(isEventSelected & AliVEvent::kAny) ) return;
    
    //TString firedTriggerClasses = static_cast<const AliAODEvent*>(InputEvent())->GetFiredTriggerClasses();
    TString firedTriggerClasses = InputEvent()->GetFiredTriggerClasses();
        
    if ( ! ( firedTriggerClasses.Contains("CVLN_B2-B-NOPF-ALLNOTRD") || firedTriggerClasses.Contains("CVLN_R1-B-NOPF-ALLNOTRD") || firedTriggerClasses.Contains("CSEMI_R1-B-NOPF-ALLNOTRD") ) ) return;
  }


  /////////////////
  // centrality
  /////////////////
  
  //AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  //if(!esd) return;
  AliCentrality *centrality = fInputEvent->GetCentrality();
  if(!centrality) return;
  cntr = centrality->GetCentralityPercentile("V0M");
  //printf("Got the centrality %f\n",cntr);
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
  valuensparsea[3] = binct;  
  valuensparseabis[1] = binct;  
  valuensparsee[1] = binct;    
  valuensparsef[1] = binctMore;  
  valuensparsefsin[1] = binct;  
  valuensparsefbis[3] = binctMore;  
  valuensparsefbissin[3] = binct;  
  valuensparseg[1] = binct;
  valuensparseh[1] = binct; 
  valuefractioncont[1] = binct;
  valuensparsehprofile[1] = binct; 
  valuecossinephiep[2] = binctMore;
  valuensparseMCSourceDeltaPhiMaps[0] = binct;
  valuedeltaphicontamination[1] = binct;
 
  //////////////////////
  // run number
  //////////////////////

  Int_t runnumber = fInputEvent->GetRunNumber();
  //printf("Run number %d\n",runnumber);
   
  if(!fPID->IsInitialized()){
    // Initialize PID with the given run number
    fPID->InitializePID(runnumber);
  }
  if(!fPIDTOFOnly->IsInitialized()){
    // Initialize PID with the given run number
    fPIDTOFOnly->InitializePID(runnumber);
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
  fPIDTOFOnly->SetPIDResponse(pidResponse);
  if(fMonitorPhotonic) fBackgroundSubtraction->InitRun(fInputEvent,pidResponse);

  fHistEV->Fill(binctt,0.0);

  //////////////////
  // Event cut
  //////////////////
  if(!fHFECuts->CheckEventCuts("fEvRecCuts", fInputEvent)) {
    //printf("Does not pass the event cut\n");
    PostData(1, fListHist);
    return;
  }

  fHistEV->Fill(binctt,1.0);


  ///////////////////////////////////////////////////////////
  // PileUpCut
  ///////////////////////////////////////////////////////////

  Float_t multTPC(0.); // tpc mult estimate
  Float_t multGlob(0.); // global multiplicity
  const Int_t nGoodTracks = fInputEvent->GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill tpc mult
    AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(iTracks));
    if (!trackAOD) continue;
    if (!(trackAOD->TestFilterBit(1))) continue;
    if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70)  || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2)) continue;
    multTPC++;
  }
  for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) { // fill global mult
    AliAODTrack* trackAOD = dynamic_cast<AliAODTrack*>(fInputEvent->GetTrack(iTracks));
    if (!trackAOD) continue;
    if (!(trackAOD->TestFilterBit(16))) continue;
    if ((trackAOD->Pt() < .2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > .8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1)) continue;
    Double_t b[2] = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
    if (!(trackAOD->PropagateToDCA(fInputEvent->GetPrimaryVertex(), fInputEvent->GetMagneticField(), 100., b, bCov))) continue;
    if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3)) continue;
    multGlob++;
  } //track loop

  Double_t pileup[4];
  pileup[0]=fInputEvent->GetCentrality()->GetCentralityPercentile("V0M\n");
  pileup[1]=fInputEvent->GetCentrality()->GetCentralityPercentile("TRK\n");
  pileup[2]=multTPC;
  pileup[3]=multGlob;
  fHistPileUp->Fill(pileup);

  if(fPileUpCut){
    if (TMath::Abs(pileup[0]-pileup[1]) > 5) {
      //printf("Does not pass the centrality correlation cut\n");
      return;
    }
    if(multTPC < (-36.81+1.48*multGlob) && multTPC > (63.03+1.78*multGlob)){
      //printf("Does not pass the multiplicity correlation cut\n");
      return;
    }
  }
 
  // AliVVZERO* vzeroData=fInputEvent->GetVZEROData();
  // Double_t mult[3],multV0A(0),multV0C(0);
  // for(Int_t i=0; i<32; ++i) {
  //   multV0A += vzeroData->GetMultiplicityV0A(i);
  //   multV0C += vzeroData->GetMultiplicityV0C(i);
  // }

  // int ntrk=0;
  // for(Int_t k = 0; k < fInputEvent->GetNumberOfTracks(); k++){
  //   AliVTrack *track = (AliVTrack *) fInputEvent->GetTrack(k);
  //   if(!track) continue;
  //   if(!(track->GetStatus()&AliVTrack::kITSrefit)) continue;
  //   if(!(track->GetStatus()&AliVTrack::kTPCrefit)) continue;
  //   ntrk++;
  // }
    
  // mult[0]=fInputEvent->GetNumberOfTracks();
  // mult[1]=multV0A+multV0C;
  // mult[2]=binctMore;
  // fHistPileUp->Fill(mult);

  // if(fUpperPileUpCut&&fLowerPileUpCut){
  //   if((mult[0]<fLowerPileUpCut->Eval(mult[1])) || 
  //      (mult[0]>fUpperPileUpCut->Eval(mult[1]))){
  //     printf("Does not pass the pileup cut\n");
  //     PostData(1, fListHist);
  //     return;
  //   }
  // }

  ////////////////////////////////////  
  // First method event plane
  ////////////////////////////////////

  AliEventplane* vEPa = fInputEvent->GetEventplane();
  Float_t eventPlanea = 0.0;
  Float_t eventPlaneTPC = 0.0;
  Float_t eventPlaneV0A = 0.0;
  Float_t eventPlaneV0C = 0.0;
  Float_t eventPlaneV0 = 0.0;
  TVector2 *qTPC = 0x0;
  TVector2 *qsub1a = 0x0;
  TVector2 *qsub2a = 0x0;
  TVector2 qV0A,qV0C,qV0,*qAna;
  
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

    Double_t qVx, qVy;  //TR: info
    eventPlaneV0 = TVector2::Phi_0_2pi(vEPa->CalculateVZEROEventPlane(fInputEvent,10,2,qVx,qVy));
    if(eventPlaneV0 > TMath::Pi()) eventPlaneV0 = eventPlaneV0 - TMath::Pi();
    qV0.Set(qVx,qVy);
    eventPlaneV0A = TVector2::Phi_0_2pi(vEPa->CalculateVZEROEventPlane(fInputEvent,8,2,qVx,qVy));
    if(eventPlaneV0A > TMath::Pi()) eventPlaneV0A = eventPlaneV0A - TMath::Pi();
    qV0A.Set(qVx,qVy);
    eventPlaneV0C = TVector2::Phi_0_2pi(vEPa->CalculateVZEROEventPlane(fInputEvent,9,2,qVx,qVy));
    if(eventPlaneV0C > TMath::Pi()) eventPlaneV0C = eventPlaneV0C - TMath::Pi();
    qV0C.Set(qVx,qVy);

    if(eventPlaneV0<-900) return;
    if(eventPlaneV0A<-900) return;
    if(eventPlaneV0C<-900) return;


    eventPlaneV0=TVector2::Phi_0_2pi(eventPlaneV0);
    eventPlaneV0A=TVector2::Phi_0_2pi(eventPlaneV0A);
    eventPlaneV0C=TVector2::Phi_0_2pi(eventPlaneV0C);
  }


  // TPC

  qTPC = vEPa->GetQVector(); 
  Double_t qx = -1.0;
  Double_t qy = -1.0;
  if(qTPC) {
    qx = qTPC->X();
    qy = qTPC->Y();
  }  
  TVector2 qVectorfortrack;
  qVectorfortrack.Set(qx,qy);
  eventPlaneTPC = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.; 

  // Choose the one used for v2

  if(fVZEROEventPlane){ //TR: info
    eventPlanea = eventPlaneV0;
    qAna = &qV0;
  }
  if(fVZEROEventPlaneA){
    eventPlanea = eventPlaneV0A;
    qAna = &qV0A;
  }
  if(fVZEROEventPlaneC){
    eventPlanea = eventPlaneV0C;
    qAna = &qV0C;
  }
  if(!fVZEROEventPlane){
    eventPlanea = eventPlaneTPC;
    qAna = &qV0C;
  }

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

  // two sub event TPC
  qsub1a = vEPa->GetQsub1();
  qsub2a = vEPa->GetQsub2();

  /////////////////////////////////////////////////////////
  // Cut for event with event plane reconstructed by all
  ////////////////////////////////////////////////////////
  
  if((!qTPC) || (!qsub1a) || (!qsub2a)) {
    //printf("No event plane\n");
    return;
  }

  eventPlanesub1a = TVector2::Phi_0_2pi(qsub1a->Phi())/2.;
  eventPlanesub2a = TVector2::Phi_0_2pi(qsub2a->Phi())/2.;
  diffsub1sub2a = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
  diffsub1sub2asin = TMath::Sin(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));


  // if ( !fDebugStreamer ) {
  //   //debug stream
  //   TDirectory *backup = gDirectory;
  //   fDebugStreamer = new TTreeSRedirector("TaskFlowTPCTOFEPSPdebug.root\n");
  //   if ( backup ) backup->cd();  //we don't want to be cd'd to the debug streamer
  // }     

  // {

  //   double v0nrom = TMath::Sqrt(qV0.X()*qV0.X()+qV0.Y()*qV0.Y());
  //   double v0Anrom = TMath::Sqrt(qV0A.X()*qV0A.X()+qV0A.Y()*qV0A.Y());
  //   double v0Cnrom = TMath::Sqrt(qV0C.X()*qV0C.X()+qV0C.Y()*qV0C.Y());
  //   double sub1nrom = TMath::Sqrt(qsub1a->X()*qsub1a->X()+qsub1a->Y()*qsub1a->Y());
  //   double sub2nrom = TMath::Sqrt(qsub2a->X()*qsub2a->X()+qsub2a->Y()*qsub2a->Y());

  //   (* fDebugStreamer) << "UserExec" <<
  //     "binct="<<binct<<
  //     "qV0="<<v0nrom<<
  //     "qV0A="<<v0Anrom<<
  //     "qV0C="<<v0Cnrom<<
  //     "qsub1a="<<sub1nrom<<
  //     "qsub2a="<<sub2nrom<<
  //     "\n";
  // }

  // three sub events in case of VZEROA and VZEROC
  if(!fSP){
    diffsubasubb = TMath::Cos(2.*(eventPlaneV0A - eventPlaneV0C));  //TR: 
    diffsubasubc = TMath::Cos(2.*(eventPlaneV0A - eventPlaneTPC));  //TR: 
    diffsubbsubc = TMath::Cos(2.*(eventPlaneV0C - eventPlaneTPC));  //TR: 
  }
  else{
    if(fVZEROEventPlaneA){
      diffsubasubb = qV0A.X()*qV0C.X()+qV0A.Y()*qV0C.Y();
      diffsubasubc = qV0A.X()*qTPC->X()+qV0A.Y()*qTPC->Y();
      diffsubbsubc = qV0C.X()*qTPC->X()+qV0C.Y()*qTPC->Y();
    }
    else if(fVZEROEventPlaneC){
      diffsubasubb = qV0C.X()*qV0A.X()+qV0C.Y()*qV0A.Y();
      diffsubasubc = qV0C.X()*qTPC->X()+qV0C.Y()*qTPC->Y();
      diffsubbsubc = qV0A.X()*qTPC->X()+qV0A.Y()*qTPC->Y();
    }
  }

  diffsubasubbsin = TMath::Sin(2.*(eventPlaneV0A - eventPlaneV0C));
  diffsubasubcsin = TMath::Sin(2.*(eventPlaneV0A - eventPlaneTPC));
  diffsubbsubcsin = TMath::Sin(2.*(eventPlaneV0C - eventPlaneTPC));
  // three sub events in case of VZERO all
  if(fVZEROEventPlane && (!fVZEROEventPlaneA) && (!fVZEROEventPlaneC)) {
    if(!fSP){
      diffsubasubb = TMath::Cos(2.*(eventPlaneV0 - eventPlanesub1a));     //TR: 
      diffsubasubc = TMath::Cos(2.*(eventPlaneV0 - eventPlanesub2a));     //TR: 
      diffsubbsubc = TMath::Cos(2.*(eventPlanesub1a - eventPlanesub2a));  //TR: 
    }
    else{
      diffsubasubb = qV0.X()*qsub1a->X()+qV0.Y()*qsub1a->Y();	       
      diffsubasubc = qV0.X()*qsub2a->X()+qV0.Y()*qsub2a->Y();	 
      diffsubbsubc = qsub1a->X()*qsub2a->X()+qsub1a->Y()*qsub2a->Y();
    }
    
    diffsubasubbsin = TMath::Sin(2.*(eventPlaneV0 - eventPlanesub1a));
    diffsubasubcsin = TMath::Sin(2.*(eventPlaneV0 - eventPlanesub2a));
    diffsubbsubcsin = TMath::Sin(2.*(eventPlanesub1a - eventPlanesub2a));
  }
  
  //////////////////////////////////////
  // AliFlowEvent  and MC event plane
  /////////////////////////////////////

  Int_t nbtracks = fInputEvent->GetNumberOfTracks();
  //printf("Number of tracks %d\n",nbtracks);

  if(fMonitorQCumulant) {

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
  }

  
  //////////////////////
  // Fill Histos
  //////////////////////

  fHistEV->Fill(binctt,2.0);
    
  // Fill
  valuensparsea[0] = eventPlaneV0A;
  valuensparsea[1] = eventPlaneV0C;
  valuensparsea[2] = eventPlaneTPC;
  if(fVZEROEventPlane && (!fVZEROEventPlaneA) && (!fVZEROEventPlaneC)) {
    // case VZERO all
    valuensparsea[0] = eventPlaneV0;
    valuensparsea[1] = eventPlanesub1a;
    valuensparsea[2] = eventPlanesub2a;
  } 
  fEventPlane->Fill(&valuensparsea[0]);

  // Fill
  if(fMonitorEventPlane) fCosSin2phiep->Fill(&valuecossinephiep[0]);
    
  if(!fVZEROEventPlane) {
    valuensparsef[0] = diffsub1sub2a;
    fCosRes->Fill(&valuensparsef[0]);
    valuensparsefsin[0] = diffsub1sub2asin;
    if(fMonitorEventPlane) fSinRes->Fill(&valuensparsefsin[0]);
    if(fMonitorEventPlane) {
      fProfileCosRes->Fill(valuensparsef[1],valuensparsef[0]);
    }
  }
  else {
    valuensparsefbis[0] = diffsubasubb;
    valuensparsefbis[1] = diffsubasubc;
    valuensparsefbis[2] = diffsubbsubc;
    fCosResabc->Fill(&valuensparsefbis[0]); //TR: info
    valuensparsefbissin[0] = diffsubasubbsin;
    valuensparsefbissin[1] = diffsubbsubcsin;
    valuensparsefbissin[2] = diffsubasubcsin;
    if(fMonitorEventPlane) fSinResabc->Fill(&valuensparsefbissin[0]);
    if(fMonitorEventPlane) {
      fProfileCosResab->Fill(valuensparsefbis[3],valuensparsefbis[0]);
      fProfileCosResac->Fill(valuensparsefbis[3],valuensparsefbis[1]);
      fProfileCosResbc->Fill(valuensparsefbis[3],valuensparsefbis[2]);
    }
  }
  
  ////////////////////////////////////////
  // Loop to determine pool background
  /////////////////////////////////////////
  if(fMonitorPhotonic) {
    fBackgroundSubtraction->FillPoolAssociatedTracks(fInputEvent,binct);
    if(mcthere) {
      fMCQA->SetCentrality(binct);
      fMCQA->SetPercentrality(static_cast<Int_t>(cntr));
      fMCQA->SetPbPb();
      //printf("Init fMCQA\n");
      fMCQA->Init();
      //printf("GetMesonKine\n");
      //if(!fAODAnalysis) fMCQA->GetMesonKine();
      // Light variant of fMCQA->GetMesonKine();
      GetMesonKine(binct);
    }
  }
    
  
  
  //////////////////////////
  // Loop over track
  //////////////////////////
  //printf("Loop over tracks\n");
  for(Int_t k = 0; k < nbtracks; k++){
      
    AliVTrack *track = (AliVTrack *) fInputEvent->GetTrack(k);
    if(!track) continue;

    // mc track as well
    if(mcthere) {
      mctrackaod = NULL;
      mctrack = NULL;
      if(fAODAnalysis){
	if(fAODArrayMCInfo){
	  Int_t label = TMath::Abs(track->GetLabel());
	  if(label && label < fAODArrayMCInfo->GetEntriesFast())
	    mctrackaod = dynamic_cast<AliAODMCParticle *>(fAODArrayMCInfo->At(label));
	}
      }
      else {
	if(mcEvent) mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));
      }
    }

    // aod filtering
    if(fAODAnalysis) {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
      if(!aodtrack){
	AliError("AOD track is not there\n");
	continue;
      }  
      //printf("Find AOD track on\n");
      if(!(aodtrack->TestFilterBit(fFilter))) continue;  // Only process AOD tracks where the HFE is set
    }
    
    valuetrackingcuts[0] = track->Pt(); 
    valuetrackingcuts[1] = 0;

    // RecKine: ITSTPC cuts  
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);
    
    // RecPrim
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    valuetrackingcuts[1] = 1; 
    if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
    
    // HFEcuts: ITS layers cuts
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    valuetrackingcuts[1] = 2; 
    if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);     
    
    // HFE cuts: TOF and mismatch flag
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTOF + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    valuetrackingcuts[1] = 3; 
    if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
    
    // HFE cuts: TPC PID cleanup
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTPC + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    valuetrackingcuts[1] = 4; 
    if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
    
    // HFEcuts: Nb of tracklets TRD0
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    valuetrackingcuts[1] = 5; 
    if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      
    //printf("Survived\n");

    /////////////////////////////////////////////////////////
    // Subtract candidate from TPC event plane
    ////////////////////////////////////////////////////////
    Float_t eventplanesubtracted = 0.0;    

    if(!fVZEROEventPlane) {
      // Subtract the tracks from the event plane
      Double_t qX = qTPC->X() - vEPa->GetQContributionX(track);  //Modify the components: subtract the track you want to look at with your analysis
      Double_t qY = qTPC->Y() - vEPa->GetQContributionY(track);  //Modify the components: subtract the track you want to look at with your analysis
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

    ////////////////////////////////
    // Determine the deltaphi bin
    ///////////////////////////////

    // in-plane
    if(((deltaphi<(TMath::Pi()/4.)) && (deltaphi>0.0)) || ((deltaphi>(3*TMath::Pi()/4.)) && (deltaphi<TMath::Pi()))) valuedeltaphicontamination[0] = 0.5;
    // out-of-plane
    if((deltaphi>(TMath::Pi()/4.)) && (deltaphi<(3*TMath::Pi()/4.))) valuedeltaphicontamination[0] = 1.5;

    ////////////////////////////////////////
    // Define variables
    ///////////////////////////////////////


    valuedeltaphicontamination[2] = track->Pt();
    valuensparsee[2] = track->Pt();
    valuensparsee[3] = track->Eta();    
    valuensparseg[2] = track->Pt();
    valuensparseh[2] = track->Pt();
    valuefractioncont[0] = track->Pt();
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
    
    if(fMonitorWithoutPID) { 
      
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
	// pid object
	AliHFEpidObject hfetrack;
	if(!fAODAnalysis){
	  hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
	  if(fVariableMultiplicity==0) 
	    hfetrack.SetMulitplicity(cntr);
	  if(fVariableMultiplicity==1)
	    hfetrack.SetMulitplicity(((AliESDEvent*)fInputEvent)->GetNumberOfESDTracks()/8.);
	  if(fVariableMultiplicity==2)
	    hfetrack.SetMulitplicity(((AliESDEvent*)fInputEvent)->GetPrimaryVertexSPD()->GetNContributors());
	}else{
	  hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
	  if(fVariableMultiplicity==0) 
	    hfetrack.SetMulitplicity(cntr);
	  if(fVariableMultiplicity==1)
	    hfetrack.SetMulitplicity(((AliAODEvent*)fInputEvent)->GetNumberOfESDTracks()/8.);
	  if(fVariableMultiplicity==2)
	    hfetrack.SetMulitplicity(((AliAODEvent*)fInputEvent)->GetPrimaryVertexSPD()->GetNContributors());
	}
	hfetrack.SetRecTrack(track);
	hfetrack.SetCentrality((Int_t)binct);
	hfetrack.SetPbPb();

	// Only TOF PID
	if(fMonitorContamination) {
	  if(fPIDTOFOnly->IsSelected(&hfetrack,0x0,"recTrackCont",0x0)) {
	    Float_t nsigma = pidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
	    valuedeltaphicontamination[3] = nsigma;
	    fDeltaPhiMapsContamination->Fill(&valuedeltaphicontamination[0]);
	  }
	}

	// Complete PID TOF+TPC
	if(!fPID->IsSelected(&hfetrack,0x0,"recTrackCont",fPIDqa)) {
	  continue;
	}
	//printf("Pass the TOF-TPC PID\n");
	
      }
      else {
	if(fAODAnalysis) {
	  if(!mctrackaod) continue;
	  if(TMath::Abs(mctrackaod->GetPdgCode())!=11) continue;
	}
	else {
	  if(!mctrack) continue;
	  //printf("PdgCode %d\n",TMath::Abs(mctrack->Particle()->GetPdgCode()));
	  if(TMath::Abs(mctrack->Particle()->GetPdgCode())!=11) continue;
	}
      }
    }


    /////////////////////////////////////////////////////////////////////////////
    // Add candidate to AliFlowEvent for POI and subtract from RP if needed
    ////////////////////////////////////////////////////////////////////////////
    if(fMonitorQCumulant) {
      Int_t idtrack = static_cast<AliVTrack*>(track)->GetID();
      Bool_t found = kFALSE;
      Int_t numberoffound = 0;
      //printf("A: Number of tracks %d\n",fflowEvent->NumberOfTracks());
      for(Int_t iRPs=0; iRPs< fflowEvent->NumberOfTracks(); iRPs++) {
	AliFlowTrack *iRP = (AliFlowTrack*) (fflowEvent->GetTrack(iRPs));
	if(!iRP) continue;
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
    }
    
  
    /////////////////////
    // Fill THnSparseF
    /////////////////////

    //
    valuensparseabis[0] = eventplanesubtracted;
    if((fillEventPlane) && (fMonitorEventPlane)) fEventPlaneaftersubtraction->Fill(&valuensparseabis[0]);
    

    if(fMonitorEventPlane) 
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
      fCosPhiMaps->Fill(&valuensparseh[0]); //TR: fCosPhiQSum+=valuensparseh[0]*TMath:Sqrt(qAna->X()*qAna->X()+qAna->Y()*qAna->Y()); fCosPhiQN++;
      //printf("Filled fCosPhiMaps\n");
      if((valuefractioncont[1] >=0) && (valuefractioncont[1] < 11)){
	if(fContamination[((Int_t)valuefractioncont[1])]){
	  //printf("Contamination\n");
	  Double_t weight = 1.;
	  if(fAsFunctionOfP) weight = fContamination[((Int_t)valuefractioncont[1])]->Eval(track->P());
	  else weight = fContamination[((Int_t)valuefractioncont[1])]->Eval(track->Pt());
	  if(weight<0.0) weight=0.0;
	  if(weight>1.0) weight=1.0;
	  fFractionContamination->Fill(&valuefractioncont[0],weight);
	  if(fv2contamination[((Int_t)valuefractioncont[1])]){
	    Double_t v2 =  fv2contamination[((Int_t)valuefractioncont[1])]->Eval(track->Pt());
	    //printf("value v2 %f, contamination %f and pt %f centrality %d\n",v2,weight,track->Pt(),(Int_t)valuefractioncont[1]);
	    //printf("Check for centrality 3: value v2 %f, contamination %f\n",fv2contamination[3]->Eval(track->Pt()),fContamination[3]->Eval(track->P()));
	    //printf("Check for centrality 4: value v2 %f, contamination %f\n",fv2contamination[4]->Eval(track->Pt()),fContamination[4]->Eval(track->P()));
	    //printf("Check for centrality 5: value v2 %f, contamination %f\n",fv2contamination[5]->Eval(track->Pt()),fContamination[5]->Eval(track->P()));
	    fContaminationv2->Fill(valuefractioncont[1],valuefractioncont[0],v2,weight);
	  }
	  fContaminationmeanpt->Fill(valuefractioncont[1],valuefractioncont[0],TMath::Abs(track->Pt()));
	}     
      }
      if(fMonitorEventPlane) {
	if(fSP)
	  valuensparseh[0] *= TMath::Sqrt(qAna->X()*qAna->X()+qAna->Y()*qAna->Y());
	fProfileCosPhiMaps->Fill(valuensparsehprofile[1],valuensparsehprofile[2],valuensparseh[0]);  //TR: info
      }
    }

    
    if(fMonitorPhotonic) {
      //printf("Background subtraction\n");
      Int_t indexmother = -1;
      Int_t source = -1;
      Int_t mcQAsource = -1;
      Double_t weightNonPhotonicFactor = 1.;
      Bool_t fillphotonic = kTRUE;
      if(mcthere) {
	if(fAODAnalysis) {
	  //printf("AOD\n");
	  source = fBackgroundSubtraction->FindMother(TMath::Abs(track->GetLabel()),indexmother);
	  //printf("source %d\n",source);
	  if(fBackgroundSubtraction->GetLevelBack()>=0) {
	    //printf("Background level %d\n",fBackgroundSubtraction->GetLevelBack());
	    // weights
	    if(fMCQA) {
	      mcQAsource = fMCQA->GetElecSource(mctrackaod, kTRUE);
	      weightNonPhotonicFactor = TMath::Abs(fMCQA->GetWeightFactor(mctrackaod, fBackgroundSubtraction->GetLevelBack())); // positive:conversion e, negative: nonHFE
	      //printf("weight %f\n",weightNonPhotonicFactor);
	    }
	  }
	   if(fSelectGenerator>=0) {
	     //printf("SelectGenerator %d\n",fSelectGenerator);
	     // select generator
	     Int_t Prim = GetPrimary(TMath::Abs(track->GetLabel()));
	     //printf("Prim %d\n",Prim);
	     if(Prim>=0) {
	       AliAODMCParticle *AODMCtrackPrim = (AliAODMCParticle*)fAODArrayMCInfo->At(TMath::Abs(Prim));
	       Int_t trkIndexPrim = AODMCtrackPrim->GetLabel();//gives index of the particle in original MCparticle array
	       if(IsFromHijing(trkIndexPrim)) {
		 if(fSelectGenerator==1) fillphotonic = kFALSE;
		 //printf("Is from Hijing\n");
	       }
	       else {
		 if(fSelectGenerator==0) fillphotonic = kFALSE;
	       }
	     }
	     else fillphotonic = kFALSE;
	   }
	}
	else {
	  //printf("ESD\n");
	  if(mctrack) source = fBackgroundSubtraction->FindMother(mctrack->GetLabel(),indexmother);
	  if(fBackgroundSubtraction->GetLevelBack()>=0) {
	    // weights
	    if(fMCQA) {
	      mcQAsource = fMCQA->GetElecSource(mctrack, kTRUE);
	      weightNonPhotonicFactor = TMath::Abs(fMCQA->GetWeightFactor(mctrack, fBackgroundSubtraction->GetLevelBack())); // positive:conversion e, negative: nonHFE 
	    }
	  }
	  if(fSelectGenerator>=0) {
	    // select generator
	    Int_t Prim = GetPrimary(mctrack->GetLabel());
	    if(Prim>=0) {
	      AliMCParticle *MCtrackPrim = (AliMCParticle *)(fMCEvent->GetTrack(TMath::Abs(Prim)));
	      Int_t trkIndexPrim = MCtrackPrim->GetLabel();//gives index of the particle in original MCparticle array
	      if(IsFromHijing(trkIndexPrim)) {
		if(fSelectGenerator==1) fillphotonic = kFALSE;
	      }
	      else {
		if(fSelectGenerator==0) fillphotonic = kFALSE;
	      }
	    }
	    else fillphotonic = kFALSE;
	  }
	}
      }
      //printf("source %d and weight %f\n",source,weightNonPhotonicFactor);
      if(fillphotonic) {
	//printf("Fill histos\n");
	fBackgroundSubtraction->LookAtNonHFE(k, track, fInputEvent, weightNonPhotonicFactor, binct, deltaphi, source, indexmother, mcQAsource);
      }
      //printf("Look At Non HFE\n");
    }
    
    
  }
  
  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////AFTERBURNER
  if (fAfterBurnerOn &  fMonitorQCumulant)
    {
      fflowEvent->AddFlow(fV1,fV2,fV3,fV4,fV5);     //add flow
      fflowEvent->CloneTracks(fNonFlowNumberOfTrackClones); //add nonflow by cloning tracks
    }
  //////////////////////////////////////////////////////////////////////////////
  
  
  
  //for(Int_t bincless = 0; bincless < fNbBinsCentralityQCumulant; bincless++) {
  //  if((fBinCentralityLess[bincless]< cntr) && (cntr < fBinCentralityLess[bincless+1])) PostData(bincless+2,fflowEvent);
  //}
  
  
  
  if(fMonitorPhotonic) fBackgroundSubtraction->CountPoolAssociated(fInputEvent,binct);
  
  //printf("Finish\n");
  
  PostData(1, fListHist);
  
  
  
}
//______________________________________________________________________________
AliFlowCandidateTrack *AliAnalysisTaskFlowTPCTOFEPSP::MakeTrack( Double_t mass, 
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
Double_t AliAnalysisTaskFlowTPCTOFEPSP::GetPhiAfterAddV2(Double_t phi,Double_t reactionPlaneAngle) const
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
//___________________________________________________________________________
Bool_t AliAnalysisTaskFlowTPCTOFEPSP::IsFromHijing(Int_t Index)
{
  //
  //Check if the particle is from Hijing or Enhanced event
  //

  Int_t nBG =-1;
  
  if(fAODAnalysis) {
    AliAODMCHeader *mcHeader;
   
    
    mcHeader = dynamic_cast<AliAODMCHeader*>(fInputEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
      //printf("Could not find MC Header in AOD\n");
        return kFALSE;
    }
    
    TList *List = mcHeader->GetCocktailHeaders();
    AliGenHijingEventHeader* hijingH = dynamic_cast<AliGenHijingEventHeader*>(List->FindObject("Hijing"));
    if (!hijingH){
      //printf("no GenHijing header");
        return kFALSE;
    }
    nBG = hijingH->NProduced();
  }
  else {

    if(!fMCEvent) return kFALSE;
    
    TString genname=fMCEvent->GenEventHeader()->ClassName();
    //Int_t typeHF=-1;
    TList* List=0x0;
    if(genname.Contains("CocktailEventHeader")){
      AliGenCocktailEventHeader *cockhead=(AliGenCocktailEventHeader*)fMCEvent->GenEventHeader();
      List=cockhead->GetHeaders();
      AliGenHijingEventHeader* hijingH = dynamic_cast<AliGenHijingEventHeader*>(List->FindObject("Hijing"));
      if (!hijingH){
        //printf("no GenHijing header");
        return kFALSE;
      }
      nBG = hijingH->NProduced();
    }
    else {
      //printf("No cocktail generators\n");
      return kFALSE;
    }

  }
  
    return (Index < nBG);
}
//_________________________________________

Int_t AliAnalysisTaskFlowTPCTOFEPSP::GetPrimary(Int_t id)
{
    
  //
  // Return number of primary that has generated track
  //

  if(fAODAnalysis && (!fAODArrayMCInfo)) return -1;
  if((!fAODAnalysis) && (!fMCEvent)) return -1;
  
  int current, parent;
    parent=id;
    while (1) {
        current=parent;
	if(fAODAnalysis) {
	  AliAODMCParticle *Part = (AliAODMCParticle*)fAODArrayMCInfo->At(current);
	  parent=Part->GetMother();
	}
	else {
	  AliMCParticle *Part = (AliMCParticle *) fMCEvent->GetTrack(current);
	  parent = ((TParticle *)Part->Particle())->GetFirstMother();
	}
        //  cout << "GetPartArr momid :"  << parent << endl;
        if(parent<0) return current;
    }
}
//__________________________________________
void AliAnalysisTaskFlowTPCTOFEPSP::GetMesonKine(Int_t centrality) 
{
  //
  // Get meson pt spectra
  //

  if(!fMCQA) return;
  if(fAODAnalysis && (!fAODArrayMCInfo)) return;
  if((!fAODAnalysis) && (!fMCEvent)) return;

  TList *qalist = fMCQA->GetList();
  //qalist->Print();
  TList *mcQACollectionlist = (TList *) qalist->FindObject("list_TaskMCQA");
  if(!mcQACollectionlist) {
    //printf("No MC QA collection list\n");
    return;
  }
  //mcQACollectionlist->Print();

  if(centrality>=11) {
    //printf("Centrality out of histogram array limits: %d", centrality);
    return;
  }


  if(!fAODAnalysis) {

     AliVParticle *mctrack2 = NULL;
     AliMCParticle *mctrack0 = NULL;
     TH2F *histo = 0x0;
     TH1F *histo1D = 0x0;
    
     for(Int_t imc = 0; imc <fMCEvent->GetNumberOfPrimaries(); imc++){
       if(!(mctrack2 = fMCEvent->GetTrack(imc))) continue;
       TParticle* mcpart0 = fMCEvent->Stack()->Particle(imc);
       if(!mcpart0) continue;
       mctrack0 = dynamic_cast<AliMCParticle *>(mctrack2);
       if(!mctrack0) continue;
       if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))>0.8) continue;
       // mc source
       Float_t mcsource = fMCQA->GetElecSource(mctrack0, kFALSE);
       // generator
       Bool_t fillphotonic = kTRUE;
       if(fSelectGenerator>=0) {
	 // select generator
	 Int_t Prim = GetPrimary(mctrack2->GetLabel());
	 if(Prim>=0) {
	   AliMCParticle *MCtrackPrim = (AliMCParticle *)(fMCEvent->GetTrack(TMath::Abs(Prim)));
	   Int_t trkIndexPrim = MCtrackPrim->GetLabel();//gives index of the particle in original MCparticle array
	   if(IsFromHijing(trkIndexPrim)) {
	     if(fSelectGenerator==1) fillphotonic = kFALSE;
	   }
	   else {
	     if(fSelectGenerator==0) fillphotonic = kFALSE;
	   }
	 }
	 else fillphotonic = kFALSE;
       }
       if(!fillphotonic) continue;

       
       if(TMath::Abs(mctrack0->PdgCode()) == 111) // pi0 
       {
	 histo = (TH2F *) mcQACollectionlist->FindObject(Form("pionspectraLog2D_centrbin%i",centrality));
	 if(histo) histo->Fill(mcsource,mctrack0->Pt()); 
       }
       else if(TMath::Abs(mctrack0->PdgCode()) == 221) // eta 
	 {
	   histo = (TH2F *) mcQACollectionlist->FindObject(Form("etaspectraLog2D_centrbin%i",centrality));
	   if(histo) histo->Fill(mcsource,mctrack0->Pt());
	 }
       else if(TMath::Abs(mctrack0->PdgCode()) == 223) // omega
	 {
	   histo = (TH2F *) mcQACollectionlist->FindObject(Form("omegaspectraLog2D_centrbin%i",centrality));
	   if(histo) histo->Fill(mcsource,mctrack0->Pt());
	 }
       else if(TMath::Abs(mctrack0->PdgCode()) == 333) // phi 
	 {
	   histo = (TH2F *) mcQACollectionlist->FindObject(Form("phispectraLog2D_centrbin%i",centrality));
	   if(histo) histo->Fill(mcsource,mctrack0->Pt());
         }
       else if(TMath::Abs(mctrack0->PdgCode()) == 331) // eta prime
	 {
	   histo = (TH2F *) mcQACollectionlist->FindObject(Form("etapspectraLog2D_centrbin%i",centrality));
	   if(histo) histo->Fill(mcsource,mctrack0->Pt());
	 }
       else if(TMath::Abs(mctrack0->PdgCode()) == 113) // rho
	 {
	   histo = (TH2F *) mcQACollectionlist->FindObject(Form("rhospectraLog2D_centrbin%i",centrality));
	   if(histo) histo->Fill(mcsource,mctrack0->Pt());
	 }
       else if(TMath::Abs(mctrack0->PdgCode()) == 321) // kaon+-
       {
	 histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("kaonspectraLog_centrbin%i",centrality));
	 if(histo1D) {
	   histo1D->Fill(mctrack0->Pt());
	   //printf("Fill histo %s\n",Form("kaonspectraLog_centrbin%i",centrality));
	 }
       }
       else if(TMath::Abs(mctrack0->PdgCode()) == 130) // k0L
	 {
	   histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("k0LspectraLog_centrbin%i",centrality));
	   if(histo1D) {
	     histo1D->Fill(mctrack0->Pt());
	     //printf("Fill histo %s\n",Form("k0LpspectraLog_centrbin%i",centrality));
	   }
	 }
       else if(TMath::Abs(mctrack0->PdgCode()) == 310) // k0S
	 {
	   histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("k0SspectraLog_centrbin%i",centrality));
	   if(histo1D) {
	     histo1D->Fill(mctrack0->Pt());
	     //printf("Fill histo %s\n",Form("k0SpspectraLog_centrbin%i",centrality));
	   }
	 }
       else if(TMath::Abs(mctrack0->PdgCode()) == 3122) // lamda
	 {
	   histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("lamdaspectraLog_centrbin%i",centrality));
	   if(histo1D) {
	     histo1D->Fill(mctrack0->Pt());
	     //printf("Fill histo %s\n",Form("lamdaspectraLog_centrbin%i",centrality));
	   }
	 }
       else if(TMath::Abs(mctrack0->PdgCode()) == 3222) // sigma
	 {
	   histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("sigmaspectraLog_centrbin%i",centrality));
	   if(histo1D) {
	     histo1D->Fill(mctrack0->Pt());
	     //printf("Fill histo %s\n",Form("sigmaspectraLog_centrbin%i",centrality));
	   }	 
	 }
     }  
  } else {

 
    AliAODMCParticle *mctrack0 = NULL;
    TH2F *histo = 0x0;
    TH1F *histo1D = 0x0;
    
    for(Int_t imc=0; imc< fAODArrayMCInfo->GetEntries(); imc++){
      mctrack0 = (AliAODMCParticle*)fAODArrayMCInfo->At(imc);
      // IsPrimary
      Bool_t primMC = mctrack0->IsPrimary();
      if(!primMC) continue;
      // Eta cut
      if(TMath::Abs(mctrack0->Eta()) > 0.8) continue;
      // mc source
      Float_t mcsource = fMCQA->GetElecSource(mctrack0, kFALSE);
      //printf("PdgCode %d\n",mctrack0->PdgCode());
      // Generator
      Bool_t fillphotonic = kTRUE;
      if(fSelectGenerator>=0) {
	//printf("SelectGenerator %d\n",fSelectGenerator);
	// select generator
	Int_t Prim = GetPrimary(imc);
	//printf("Prim %d\n",Prim);
	if(Prim>=0) {
	  AliAODMCParticle *AODMCtrackPrim = (AliAODMCParticle*)fAODArrayMCInfo->At(TMath::Abs(Prim));
	  Int_t trkIndexPrim = AODMCtrackPrim->GetLabel();//gives index of the particle in original MCparticle array
	  if(IsFromHijing(trkIndexPrim)) {
	    if(fSelectGenerator==1) fillphotonic = kFALSE;
	    //printf("Is from Hijing\n");
	  }
	  else {
	    if(fSelectGenerator==0) fillphotonic = kFALSE;
	  }
	}
	else fillphotonic = kFALSE;
      }
      if(!fillphotonic) continue;
      //printf("Fill\n");
      // Fill
      if(TMath::Abs(mctrack0->PdgCode()) == 111) // pi0 
       { 
	 histo = (TH2F *) mcQACollectionlist->FindObject(Form("pionspectraLog2D_centrbin%i",centrality));
	 if(histo) {
	   histo->Fill(mcsource,mctrack0->Pt());
	   //printf("Fill histo %s\n",Form("pionspectraLog2D_centrbin%i",centrality));
	 }
       }
      else if(TMath::Abs(mctrack0->PdgCode()) == 221) // eta 
	{
	  histo = (TH2F *) mcQACollectionlist->FindObject(Form("etaspectraLog2D_centrbin%i",centrality));
	  if(histo) {
	    histo->Fill(mcsource,mctrack0->Pt());
	    //printf("Fill histo %s\n",Form("etaspectraLog2D_centrbin%i",centrality));
	  }
	}
      else if(TMath::Abs(mctrack0->PdgCode()) == 223) // omega
	 {
	   histo = (TH2F *) mcQACollectionlist->FindObject(Form("omegaspectraLog2D_centrbin%i",centrality));
	   if(histo) {
	     histo->Fill(mcsource,mctrack0->Pt());
	     //printf("Fill histo %s\n",Form("omegaspectraLog2D_centrbin%i",centrality));
	   }
	 }
      else if(TMath::Abs(mctrack0->PdgCode()) == 333) // phi 
	{
	  histo = (TH2F *) mcQACollectionlist->FindObject(Form("phispectraLog2D_centrbin%i",centrality));
	  if(histo) {
	    histo->Fill(mcsource,mctrack0->Pt());
	    //printf("Fill histo %s\n",Form("phispectraLog2D_centrbin%i",centrality));
	  }
	}
      else if(TMath::Abs(mctrack0->PdgCode()) == 331) // eta prime
	{
	  histo = (TH2F *) mcQACollectionlist->FindObject(Form("etapspectraLog2D_centrbin%i",centrality));
	  if(histo) {
	    histo->Fill(mcsource,mctrack0->Pt());
	    //printf("Fill histo %s\n",Form("etapspectraLog2D_centrbin%i",centrality));
	  }
	}
      else if(TMath::Abs(mctrack0->PdgCode()) == 113) // rho
	{
	  histo = (TH2F *) mcQACollectionlist->FindObject(Form("rhospectraLog2D_centrbin%i",centrality));
	  if(histo) {
	    histo->Fill(mcsource,mctrack0->Pt());
	    //printf("Fill histo %s\n",Form("rhospectraLog2D_centrbin%i",centrality));
	  }
	}
      else if(TMath::Abs(mctrack0->PdgCode()) == 321) // kaon+-
       {
	 histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("kaonspectraLog_centrbin%i",centrality));
	 if(histo1D) {
	   histo1D->Fill(mctrack0->Pt());
	   //printf("Fill histo %s\n",Form("kaonspectraLog_centrbin%i",centrality));
	 }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 130) // k0L
       {
	 histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("k0LspectraLog_centrbin%i",centrality));
	 if(histo1D) {
	   histo1D->Fill(mctrack0->Pt());
	   //printf("Fill histo %s\n",Form("k0LpspectraLog_centrbin%i",centrality));
	 }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 310) // k0S
       {
	 histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("k0SspectraLog_centrbin%i",centrality));
	 if(histo1D) {
	   histo1D->Fill(mctrack0->Pt());
	   //printf("Fill histo %s\n",Form("k0SpspectraLog_centrbin%i",centrality));
	 }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 3122) // lamda
       {
	 histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("lamdaspectraLog_centrbin%i",centrality));
	 if(histo1D) {
	   histo1D->Fill(mctrack0->Pt());
	   //printf("Fill histo %s\n",Form("lamdaspectraLog_centrbin%i",centrality));
	 }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 3222) // sigma
       {
	 histo1D = (TH1F *) mcQACollectionlist->FindObject(Form("sigmaspectraLog_centrbin%i",centrality));
	 if(histo1D) {
	   histo1D->Fill(mctrack0->Pt());
	   //printf("Fill histo %s\n",Form("sigmaspectraLog_centrbin%i",centrality));
	 }	 
       }
    }  
  }
}
