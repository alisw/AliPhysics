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


//____________________________________________________________________
AliAnalysisTaskFlowTPCTOFEPSP::AliAnalysisTaskFlowTPCTOFEPSP() :
  AliAnalysisTaskSE(),
  fListHist(0x0), 
  fAODAnalysis(kFALSE),
  fUseFilterAOD(kFALSE),
  fApplyCut(kTRUE),
  fFilter(1<<4),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fBackgroundSubtraction(NULL),
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
  fSP(kFALSE),
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
  fMonitorEventPlane(kFALSE),
  fMonitorContamination(kFALSE),
  fMonitorPhotonic(kFALSE),
  fMonitorWithoutPID(kFALSE),
  fMonitorTrackCuts(kFALSE),
  fMonitorQCumulant(kFALSE),
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fRejectKinkMother(kFALSE),
  fPID(0),
  fPIDTOFOnly(0),
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
  fHistPileUp(0),
  fPileUpCut(kFALSE),
  fEventPlane(0x0),
  fEventPlaneaftersubtraction(0x0),
  fFractionContamination(0x0),
  fContaminationv2(0x0),
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
  fDeltaPhiMapsTaggedPhotonic(0x0),
  //fCosPhiMapsTaggedPhotonic(0x0),
  fDeltaPhiMapsTaggedNonPhotonic(0x0),
  //fCosPhiMapsTaggedNonPhotonic(0x0),
  fDeltaPhiMapsTaggedPhotonicLS(0x0),
  //fCosPhiMapsTaggedPhotonicLS(0x0),
  fMCSourceDeltaPhiMaps(0x0),
  fOppSignDeltaPhiMaps(0x0),
  fSameSignDeltaPhiMaps(0x0),
  fOppSignAngle(0x0),
  fSameSignAngle(0x0),
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
   
}
//______________________________________________________________________________
AliAnalysisTaskFlowTPCTOFEPSP:: AliAnalysisTaskFlowTPCTOFEPSP(const char *name) :
  AliAnalysisTaskSE(name),
  fListHist(0x0),
  fAODAnalysis(kFALSE),
  fUseFilterAOD(kFALSE),
  fApplyCut(kTRUE),
  fFilter(1<<4), 
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fBackgroundSubtraction(NULL),
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
  fSP(kFALSE),
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
  fMonitorEventPlane(kFALSE),
  fMonitorContamination(kFALSE),
  fMonitorPhotonic(kFALSE),
  fMonitorWithoutPID(kFALSE),
  fMonitorTrackCuts(kFALSE),
  fMonitorQCumulant(kFALSE),
  fcutsRP(0),
  fcutsPOI(0),
  fHFECuts(0),
  fRejectKinkMother(kFALSE),
  fPID(0),
  fPIDTOFOnly(0),
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
  fHistPileUp(0),
  fPileUpCut(kFALSE),
  fEventPlane(0x0),
  fEventPlaneaftersubtraction(0x0),
  fFractionContamination(0x0),
  fContaminationv2(0x0),
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
  fDeltaPhiMapsTaggedPhotonic(0x0),
  //fCosPhiMapsTaggedPhotonic(0x0),
  fDeltaPhiMapsTaggedNonPhotonic(0x0),
  //fCosPhiMapsTaggedNonPhotonic(0x0),
  fDeltaPhiMapsTaggedPhotonicLS(0x0),
  //fCosPhiMapsTaggedPhotonicLS(0x0),
  fMCSourceDeltaPhiMaps(0x0),
  fOppSignDeltaPhiMaps(0x0),
  fSameSignDeltaPhiMaps(0x0),
  fOppSignAngle(0x0),
  fSameSignAngle(0x0),
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
  
  fPID = new AliHFEpid("hfePid");
  fPIDqa = new AliHFEpidQAmanager;

  fPIDBackground = new AliHFEpid("hfePidBackground");
  fPIDBackgroundqa = new AliHFEpidQAmanager;

  fPIDTOFOnly = new AliHFEpid("hfePidTOFOnly");

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
  fAODAnalysis(ref.fAODAnalysis), 
  fUseFilterAOD(ref.fUseFilterAOD),
  fApplyCut(ref.fApplyCut),
  fFilter(ref.fFilter),
  fAODMCHeader(ref.fAODMCHeader),
  fAODArrayMCInfo(ref.fAODArrayMCInfo),
  fBackgroundSubtraction(ref.fBackgroundSubtraction),
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
  fSP(ref.fSP),
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
  fMonitorEventPlane(ref.fMonitorEventPlane),
  fMonitorContamination(ref.fMonitorContamination),
  fMonitorPhotonic(ref.fMonitorPhotonic),
  fMonitorWithoutPID(ref.fMonitorWithoutPID),
  fMonitorTrackCuts(ref.fMonitorTrackCuts),
  fMonitorQCumulant(ref.fMonitorQCumulant),
  fcutsRP(NULL),
  fcutsPOI(NULL),
  fHFECuts(NULL),
  fRejectKinkMother(ref.fRejectKinkMother),
  fPID(NULL),
  fPIDTOFOnly(NULL),
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
  fHistPileUp(NULL),
  fPileUpCut(kFALSE),
  fEventPlane(NULL),
  fEventPlaneaftersubtraction(NULL),
  fFractionContamination(NULL),
  fContaminationv2(NULL),
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
  fDeltaPhiMapsTaggedPhotonic(NULL),
  //fCosPhiMapsTaggedPhotonic(NULL),
  fDeltaPhiMapsTaggedNonPhotonic(NULL),
  //fCosPhiMapsTaggedNonPhotonic(NULL),
  fDeltaPhiMapsTaggedPhotonicLS(NULL),
  //fCosPhiMapsTaggedPhotonicLS(NULL),
  fMCSourceDeltaPhiMaps(NULL),
  fOppSignDeltaPhiMaps(NULL),
  fSameSignDeltaPhiMaps(NULL),
  fOppSignAngle(NULL),
  fSameSignAngle(NULL),
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
  target.fAODAnalysis = fAODAnalysis;
  target.fUseFilterAOD = fUseFilterAOD;
  target.fApplyCut = fApplyCut;
  target.fFilter = fFilter;
  target.fAODMCHeader = fAODMCHeader;
  target.fAODArrayMCInfo = fAODArrayMCInfo;
  target.fBackgroundSubtraction = fBackgroundSubtraction;
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
  target.fSP = fSP;
  target.fMCPID = fMCPID;
  target.fNoPID = fNoPID;
  target.fChi2OverNDFCut = fChi2OverNDFCut;
  target.fMaxdca = fMaxdca;
  target.fMaxopeningtheta = fMaxopeningtheta;
  target.fMaxopeningphi = fMaxopeningphi;
  target.fMaxopening3D = fMaxopening3D;
  target.fMaxInvmass = fMaxInvmass;
  target.fSetMassConstraint =  fSetMassConstraint;
  target.fDebugLevel = fDebugLevel;
  target.fMonitorEventPlane = fMonitorEventPlane;
  target.fMonitorContamination = fMonitorContamination;
  target.fMonitorPhotonic = fMonitorPhotonic;
  target.fMonitorWithoutPID = fMonitorWithoutPID;
  target.fMonitorTrackCuts = fMonitorTrackCuts;
  target.fMonitorQCumulant = fMonitorQCumulant;
  target.fcutsRP = fcutsRP;
  target.fcutsPOI = fcutsPOI;
  target.fHFECuts = fHFECuts;
  target.fRejectKinkMother = fRejectKinkMother;
  target.fPID = fPID;
  target.fPIDTOFOnly = fPIDTOFOnly;
  target.fPIDqa = fPIDqa;
  target.fflowEvent = fflowEvent;
  target.fHFEBackgroundCuts = fHFEBackgroundCuts;  	 
  target.fPIDBackground = fPIDBackground; 		
  target.fPIDBackgroundqa = fPIDBackgroundqa; 		 	 
  target.fAlgorithmMA = fAlgorithmMA; 		 	 
  target.fArraytrack = fArraytrack; 		 	 
  target.fCounterPoolBackground = fCounterPoolBackground; 		 	 
  target.fHFEVZEROEventPlane = fHFEVZEROEventPlane; 	
  target.fHistEV=fHistEV;       		 	 
  target.fHistPileUp=fHistPileUp;   		 	 
  target.fPileUpCut=fPileUpCut;         		 	 
  target.fEventPlane=fEventPlane;     		 	 
  target.fEventPlaneaftersubtraction=fEventPlaneaftersubtraction; 		 	 
  target.fFractionContamination=fFractionContamination;     		 	 
  target.fContaminationv2=fContaminationv2;           		 	 
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
  target.fDeltaPhiMapsTaggedPhotonic=fDeltaPhiMapsTaggedPhotonic; 		 	 
  target.fDeltaPhiMapsTaggedNonPhotonic=fDeltaPhiMapsTaggedNonPhotonic; 		     
  target.fDeltaPhiMapsTaggedPhotonicLS=fDeltaPhiMapsTaggedPhotonicLS; 		 	 
  target.fMCSourceDeltaPhiMaps=fMCSourceDeltaPhiMaps; 		 	 
  target.fOppSignDeltaPhiMaps=fOppSignDeltaPhiMaps;   		 	 
  target.fSameSignDeltaPhiMaps=fSameSignDeltaPhiMaps; 		 	 
  target.fOppSignAngle=fOppSignAngle;         		 	 
  target.fSameSignAngle=fSameSignAngle;   		 	 
  
  
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
  if(fArraytrack) delete fArraytrack;
  if(fListHist) delete fListHist;
  if(fcutsRP) delete fcutsRP;
  if(fcutsPOI) delete fcutsPOI;
  if(fHFECuts) delete fHFECuts;
  if(fPID) delete fPID;
  if(fPIDTOFOnly) delete fPIDTOFOnly;
  //if(fPIDqa) delete fPIDqa;
  if(fflowEvent) delete fflowEvent;
  if(fHFEBackgroundCuts) delete fHFEBackgroundCuts;
  if(fPIDBackground) delete fPIDBackground;
  if(fBackgroundSubtraction) delete fBackgroundSubtraction;
  //if(fPIDBackgroundqa) delete fPIDBackgroundqa;
  //if(fHFEVZEROEventPlane) delete fHFEVZEROEventPlane;
  if ( fDebugStreamer ) delete fDebugStreamer;

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

  AliDebug(2,"test");

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: User create output objects");
 
  // AOD or ESD
  AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
    SetAODAnalysis(kTRUE);
    AliDebug(2,"Put AOD analysis on");
  } else {
    SetAODAnalysis(kFALSE);
  }

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: AOD ESD");

  // RP TRACK CUTS:
  fcutsRP = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();
  fcutsRP->SetName("StandartTPC");
  fcutsRP->SetEtaRange(-0.9,0.9);
  fcutsRP->SetQA(kTRUE);
  //TList *qaCutsRP = fcutsRP->GetQA();
  //qaCutsRP->SetName("QA_StandartTPC_RP");

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cutsRP");

  //POI TRACK CUTS:
  fcutsPOI = new AliFlowTrackCuts("dummy");
  fcutsPOI->SetParamType(AliFlowTrackCuts::kGlobal);
  fcutsPOI->SetPtRange(+1,-1); // select nothing QUICK
  fcutsPOI->SetEtaRange(+1,-1); // select nothing VZERO

  if( fflowEvent ){ 
    fflowEvent->~AliFlowEvent();
    new(fflowEvent) AliFlowEvent(fcutsRP,fcutsPOI);
  }
  else fflowEvent = new AliFlowEvent(fcutsRP,fcutsPOI);
    
  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cutsPOI");

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

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: common constants");

  
  // HFE cuts

  if(!fHFECuts){
    fHFECuts = new AliHFEcuts;
    fHFECuts->CreateStandardCuts();
  }
  fHFECuts->Initialize();
  if(fAODAnalysis) fHFECuts->SetAOD();  

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: HFE cuts");


  // PID HFE
  //fPID->SetHasMCData(HasMCData());
  if(!fPID) {
    fPID =new AliHFEpid("hfePid");
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: pid init 0");
  }
  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: pid init 1");
  if(!fPID->GetNumberOfPIDdetectors()) fPID->AddDetector("TPC", 0);
  AliDebug(2,Form("AliAnalysisTaskFlowTPCTOFEPSP: GetNumber of PID detectors %d",fPID->GetNumberOfPIDdetectors()));
  fPID->InitializePID();
  AliDebug(2,"Init ");
  fPIDqa->Initialize(fPID);
  AliDebug(2,"Init qa");
  fPID->SortDetectors();
  AliDebug(2,"Sort detectors");

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: pid and pidqa");

  if(!fPIDTOFOnly->GetNumberOfPIDdetectors()) fPIDTOFOnly->AddDetector("TPC", 0);
  fPIDTOFOnly->InitializePID();
  fPIDTOFOnly->SortDetectors();

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: pidtof");

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
  
  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: hfe background");

  // PID background HFE
  if(!fPIDBackground->GetNumberOfPIDdetectors()) fPIDBackground->AddDetector("TPC", 0);
  fPIDBackground->InitializePID();
  fPIDBackgroundqa->Initialize(fPIDBackground);
  fPIDBackground->SortDetectors();

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: pid background");

  if (fMonitorPhotonic) {
    if(!fBackgroundSubtraction) fBackgroundSubtraction = new AliHFENonPhotonicElectron();
    fBackgroundSubtraction->Init();
  }
  


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

  Int_t nBinsPt = 24;
  Double_t binLimPt[25] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2,
			   1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 3., 3.5, 4., 5.,
			   6.};


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
    AliDebug(2,Form("bin phi is %f for %d",binLimPhi[i],i));
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
    AliDebug(2,Form("bin phi is %f for %d",binLimPhi[i],i));
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

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: variables");
  
  //******************
  // Histograms
  //******************
    
  fListHist = new TList();
  fListHist->SetOwner();

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: list");

  // Minimum histos

  // Histos
  fHistEV = new TH2D("fHistEV", "events", 3, 0, 3, 3, 0,3);
  
  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: histev");

  // V0 multiplicity vs # of tracks vs centraliy
  const Int_t nDimPU=4;
  Int_t nBinPU[nDimPU] = {nBinsCEvenMore,nBinsCEvenMore,nBinsMult,nBinsMult};
  fHistPileUp = new THnSparseF("PileUp","PileUp",nDimPU,nBinPU);
  fHistPileUp->SetBinEdges(0,binLimCEvenMore);
  fHistPileUp->SetBinEdges(1,binLimCEvenMore);
  fHistPileUp->SetBinEdges(2,binLimMult);
  fHistPileUp->SetBinEdges(3,binLimMult);
  fHistPileUp->Sumw2();
  
  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: eventplane");
  
  // Event plane as function of phiep, centrality
  const Int_t nDima=4;
  Int_t nBina[nDima] = {nBinsPhi,nBinsPhi,nBinsPhi,nBinsC};
  fEventPlane = new THnSparseF("EventPlane","EventPlane",nDima,nBina);
  fEventPlane->SetBinEdges(0,binLimPhi);
  fEventPlane->SetBinEdges(1,binLimPhi);
  fEventPlane->SetBinEdges(2,binLimPhi);
  fEventPlane->SetBinEdges(3,binLimC);
  fEventPlane->Sumw2();

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: eventplane");

  // Fraction of contamination, centrality
  const Int_t nDimcont=2;
  Int_t nBincont[nDimcont] = {nBinsPt,nBinsC};
  fFractionContamination = new THnSparseF("Contamination","Contamination",nDimcont,nBincont);
  fFractionContamination->SetBinEdges(0,binLimPt);
  fFractionContamination->SetBinEdges(1,binLimC);
  fFractionContamination->Sumw2();
  //  
  fContaminationv2 = new TProfile2D("Contaminationv2","",nBinsC,binLimC,nBinsPt,binLimPt);
  fContaminationv2->Sumw2();

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: fraction of contamination");

  // Resolution cosres_abc centrality
  const Int_t nDimfbis=4;
  Int_t nBinfbis[nDimfbis] = {nBinsCos,nBinsCos,nBinsCos,nBinsCMore};
  fCosResabc = new THnSparseF("CosRes_abc","CosRes_abc",nDimfbis,nBinfbis);
  fCosResabc->SetBinEdges(0,binLimCos);
  fCosResabc->SetBinEdges(1,binLimCos);
  fCosResabc->SetBinEdges(2,binLimCos);
  fCosResabc->SetBinEdges(3,binLimCMore);
  fCosResabc->Sumw2();

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cosresabc");

  // Resolution cosres centrality
  const Int_t nDimf=2;
  Int_t nBinf[nDimf] = {nBinsCos, nBinsCMore};
  fCosRes = new THnSparseF("CosRes","CosRes",nDimf,nBinf);
  fCosRes->SetBinEdges(0,binLimCos);
  fCosRes->SetBinEdges(1,binLimCMore);
  fCosRes->Sumw2();

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cosres");

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

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: deltaphimaps");

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

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cosphimaps");

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

    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: eventplane after sub");
    
    // Monitoring of the event Plane cos(2phi) sin(2phi) centrality
    const Int_t nDimi=3;
    Int_t nBini[nDimi] = {nBinsCos, nBinsCos, nBinsCMore};
    fCosSin2phiep = new THnSparseF("CosSin2phiep","CosSin2phiep",nDimi,nBini);
    fCosSin2phiep->SetBinEdges(0,binLimCos);
    fCosSin2phiep->SetBinEdges(1,binLimCos);
    fCosSin2phiep->SetBinEdges(2,binLimCMore);
    fCosSin2phiep->Sumw2();

    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cossin2phiep");
    
    // Monitoring Event plane after subtraction of the track
    const Int_t nDime=4;
    Int_t nBine[nDime] = {nBinsCos, nBinsC, nBinsPt, nBinsEta};
    fCos2phie = new THnSparseF("cos2phie","cos2phie",nDime,nBine);
    fCos2phie->SetBinEdges(2,binLimPt);
    fCos2phie->SetBinEdges(3,binLimEta);
    fCos2phie->SetBinEdges(0,binLimCos);
    fCos2phie->SetBinEdges(1,binLimC);
    fCos2phie->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cos2phie");
    fSin2phie = new THnSparseF("sin2phie","sin2phie",nDime,nBine);
    fSin2phie->SetBinEdges(2,binLimPt);
    fSin2phie->SetBinEdges(3,binLimEta);
    fSin2phie->SetBinEdges(0,binLimCos);
    fSin2phie->SetBinEdges(1,binLimC);
    fSin2phie->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: sin2phie");
    fCos2phiep = new THnSparseF("cos2phiep","cos2phiep",nDime,nBine);
    fCos2phiep->SetBinEdges(2,binLimPt);
    fCos2phiep->SetBinEdges(3,binLimEta);
    fCos2phiep->SetBinEdges(0,binLimCos);
    fCos2phiep->SetBinEdges(1,binLimC);
    fCos2phiep->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cos2phiep");
    fSin2phiep = new THnSparseF("sin2phiep","sin2phiep",nDime,nBine);
    fSin2phiep->SetBinEdges(2,binLimPt);
    fSin2phiep->SetBinEdges(3,binLimEta);
    fSin2phiep->SetBinEdges(0,binLimCos);
    fSin2phiep->SetBinEdges(1,binLimC);
    fSin2phiep->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: sin2phiep");
    fSin2phiephiep = new THnSparseF("sin2phie_phiep","sin2phie_phiep",nDime,nBine);
    fSin2phiephiep->SetBinEdges(2,binLimPt);
    fSin2phiephiep->SetBinEdges(3,binLimEta);
    fSin2phiephiep->SetBinEdges(0,binLimCos);
    fSin2phiephiep->SetBinEdges(1,binLimC);
    fSin2phiephiep->Sumw2();  
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: sin2phiephiep");
    
    const Int_t nDimfbiss=4;
    Int_t nBinfbiss[nDimfbiss] = {nBinsCos,nBinsCos,nBinsCos,nBinsC};
    fSinResabc = new THnSparseF("SinRes_abc","SinRes_abc",nDimfbiss,nBinfbiss);
    fSinResabc->SetBinEdges(0,binLimCos);
    fSinResabc->SetBinEdges(1,binLimCos);
    fSinResabc->SetBinEdges(2,binLimCos);
    fSinResabc->SetBinEdges(3,binLimC);
    fSinResabc->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: sinresabc");
    
    // Profile cosres centrality with 3 subevents
    fProfileCosResab = new TProfile("ProfileCosRes_a_b","ProfileCosRes_a_b",nBinsCMore,binLimCMore);
    fProfileCosResab->Sumw2();
    fProfileCosResac = new TProfile("ProfileCosRes_a_c","ProfileCosRes_a_c",nBinsCMore,binLimCMore);
    fProfileCosResac->Sumw2();
    fProfileCosResbc = new TProfile("ProfileCosRes_b_c","ProfileCosRes_b_c",nBinsCMore,binLimCMore);
    fProfileCosResbc->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: profilecosresbc");
    
    //
    const Int_t nDimff=2;
    Int_t nBinff[nDimff] = {nBinsCos, nBinsC};
    fSinRes = new THnSparseF("SinRes","SinRes",nDimff,nBinff);
    fSinRes->SetBinEdges(0,binLimCos);
    fSinRes->SetBinEdges(1,binLimC);
    fSinRes->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: sinres");
    
    // Profile cosres centrality
    fProfileCosRes = new TProfile("ProfileCosRes","ProfileCosRes",nBinsCMore,binLimCMore);
    fProfileCosRes->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: profilecosres");
    
    // Profile Maps cos phi
    fProfileCosPhiMaps = new TProfile2D("ProfileCosPhiMaps","ProfileCosPhiMaps",nBinsC,binLimC,nBinsPt,binLimPt);
    fProfileCosPhiMaps->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: profilecosphimaps");

  }
  //
  // fMonitorTrackCuts
  //

  if(fMonitorTrackCuts) {
    // Debugging tracking steps
    const Int_t nDimTrStep=2;
    Int_t nBinTrStep[nDimTrStep] = {nBinsPt,nBinsStep};
    fTrackingCuts = new THnSparseF("TrackingCuts","TrackingCuts",nDimTrStep,nBinTrStep);
    fTrackingCuts->SetBinEdges(0,binLimPt);
    fTrackingCuts->SetBinEdges(1,binLimStep);
    fTrackingCuts->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: trackingcuts");
  }

  //
  // fMonitorContamination
  //

  if(fMonitorContamination) { 
    // Maps delta phi contamination
    const Int_t nDimgcont=4;
    Int_t nBingcont[nDimgcont] = {nBinsPhiLess,nBinsC,nBinsPt, nBinsTPCdEdx};
    fDeltaPhiMapsContamination = new THnSparseF("DeltaPhiMapsContamination","DeltaPhiMapsContamination",nDimgcont,nBingcont);
    fDeltaPhiMapsContamination->SetBinEdges(0,binLimPhiLess);
    fDeltaPhiMapsContamination->SetBinEdges(1,binLimC);
    fDeltaPhiMapsContamination->SetBinEdges(2,binLimPt);
    fDeltaPhiMapsContamination->SetBinEdges(3,binLimTPCdEdx);
    fDeltaPhiMapsContamination->Sumw2();  
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: deltaphimapscontamination");

  }
  //
  // fMonitorWithoutPID
  //

  if(fMonitorWithoutPID) {
    //
    const Int_t nDimgb=3;
    Int_t nBingb[nDimgb] = {nBinsPhi,nBinsC,nBinsPt};
    
    fDeltaPhiMapsBeforePID = new THnSparseF("DeltaPhiMapsBeforePID","DeltaPhiMapsBeforePID",nDimgb,nBingb);
    fDeltaPhiMapsBeforePID->SetBinEdges(0,binLimPhi);
    fDeltaPhiMapsBeforePID->SetBinEdges(1,binLimC);
    fDeltaPhiMapsBeforePID->SetBinEdges(2,binLimPt);
    fDeltaPhiMapsBeforePID->Sumw2();  
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: deltaphimapsbeforepid");
    
    const Int_t nDimhb=3;
    Int_t nBinhb[nDimhb] = {nBinsCos,nBinsC,nBinsPt};
    
    fCosPhiMapsBeforePID = new THnSparseF("CosPhiMapsBeforePID","CosPhiMapsBeforePID",nDimhb,nBinhb);
    fCosPhiMapsBeforePID->SetBinEdges(0,binLimCos);
    fCosPhiMapsBeforePID->SetBinEdges(1,binLimC);
    fCosPhiMapsBeforePID->SetBinEdges(2,binLimPt);
    fCosPhiMapsBeforePID->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cosphimapsbeforepid");
  }
  //
  // fMonitorPhotonic
  //

  if(fMonitorPhotonic) {
    
    const Int_t nDimgbp=3;
    Int_t nBingbp[nDimgbp] = {nBinsPhi,nBinsC,nBinsPt};
    
    fDeltaPhiMapsTaggedPhotonic = new THnSparseF("DeltaPhiMapsTaggedPhotonic","DeltaPhiMapsTaggedPhotonic",nDimgbp,nBingbp);
    fDeltaPhiMapsTaggedPhotonic->SetBinEdges(0,binLimPhi);
    fDeltaPhiMapsTaggedPhotonic->SetBinEdges(1,binLimC);
    fDeltaPhiMapsTaggedPhotonic->SetBinEdges(2,binLimPt);
    fDeltaPhiMapsTaggedPhotonic->Sumw2();  
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: deltaphimapstaggedphotonic");
    
    fDeltaPhiMapsTaggedNonPhotonic = new THnSparseF("DeltaPhiMapsTaggedNonPhotonic","DeltaPhiMapsTaggedNonPhotonic",nDimgbp,nBingbp);
    fDeltaPhiMapsTaggedNonPhotonic->SetBinEdges(0,binLimPhi);
    fDeltaPhiMapsTaggedNonPhotonic->SetBinEdges(1,binLimC);
    fDeltaPhiMapsTaggedNonPhotonic->SetBinEdges(2,binLimPt);
    fDeltaPhiMapsTaggedNonPhotonic->Sumw2();  
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: deltaphimapstaggednonphotonic");
    
    fDeltaPhiMapsTaggedPhotonicLS = new THnSparseF("DeltaPhiMapsTaggedPhotonicLS","DeltaPhiMapsTaggedPhotonicLS",nDimgbp,nBingbp);
    fDeltaPhiMapsTaggedPhotonicLS->SetBinEdges(0,binLimPhi);
    fDeltaPhiMapsTaggedPhotonicLS->SetBinEdges(1,binLimC);
    fDeltaPhiMapsTaggedPhotonicLS->SetBinEdges(2,binLimPt);
    fDeltaPhiMapsTaggedPhotonicLS->Sumw2();  
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: deltaphimapstaggedphotonicls");    

    /*
    const Int_t nDimhbp=3;
    Int_t nBinhbp[nDimhbp] = {nBinsCos,nBinsC,nBinsPt};
    
    fCosPhiMapsTaggedPhotonic = new THnSparseF("CosPhiMapsTaggedPhotonic","CosPhiMapsTaggedPhotonic",nDimhbp,nBinhbp);
    fCosPhiMapsTaggedPhotonic->SetBinEdges(0,binLimCos);
    fCosPhiMapsTaggedPhotonic->SetBinEdges(1,binLimC);
    fCosPhiMapsTaggedPhotonic->SetBinEdges(2,binLimPt);
    fCosPhiMapsTaggedPhotonic->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cosphimapstaggedphotonic");
    
    fCosPhiMapsTaggedNonPhotonic = new THnSparseF("CosPhiMapsTaggedNonPhotonic","CosPhiMapsTaggedNonPhotonic",nDimhbp,nBinhbp);
    fCosPhiMapsTaggedNonPhotonic->SetBinEdges(0,binLimCos);
    fCosPhiMapsTaggedNonPhotonic->SetBinEdges(1,binLimC);
    fCosPhiMapsTaggedNonPhotonic->SetBinEdges(2,binLimPt);
    fCosPhiMapsTaggedNonPhotonic->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cosphimapstaggednonphotonic");
    
    fCosPhiMapsTaggedPhotonicLS = new THnSparseF("CosPhiMapsTaggedPhotonicLS","CosPhiMapsTaggedPhotonicLS",nDimhbp,nBinhbp);
    fCosPhiMapsTaggedPhotonicLS->SetBinEdges(0,binLimCos);
    fCosPhiMapsTaggedPhotonicLS->SetBinEdges(1,binLimC);
    fCosPhiMapsTaggedPhotonicLS->SetBinEdges(2,binLimPt);
    fCosPhiMapsTaggedPhotonicLS->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: cosphimapstaggedphotonicls");
    */
    const Int_t nDimMCSource=3;
    Int_t nBinMCSource[nDimMCSource] = {nBinsC,nBinsPt,nBinsSource};
    fMCSourceDeltaPhiMaps = new THnSparseF("MCSourceDeltaPhiMaps","MCSourceDeltaPhiMaps",nDimMCSource,nBinMCSource);
    fMCSourceDeltaPhiMaps->SetBinEdges(0,binLimC);
    fMCSourceDeltaPhiMaps->SetBinEdges(1,binLimPt);
    fMCSourceDeltaPhiMaps->SetBinEdges(2,binLimSource);
    fMCSourceDeltaPhiMaps->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: mcsourcedeltaphimaps");
    
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
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: oppsigndeltaphimaps");
    
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
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: samesigndeltaphimaps");
    
    // Maps angle same sign
    const Int_t nDimAngleSameSign=3;
    Int_t nBinAngleSameSign[nDimAngleSameSign] = {nBinsAngle,nBinsC,nBinsSource};
    fSameSignAngle = new THnSparseF("SameSignAngleMaps","SameSignAngleMaps",nDimAngleSameSign,nBinAngleSameSign);
    fSameSignAngle->SetBinEdges(0,binLimAngle);
    fSameSignAngle->SetBinEdges(1,binLimC);
    fSameSignAngle->SetBinEdges(2,binLimSource);
    fSameSignAngle->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: samesignangle");
    
    // Maps angle opp sign
    const Int_t nDimAngleOppSign=3;
    Int_t nBinAngleOppSign[nDimAngleOppSign] = {nBinsAngle,nBinsC,nBinsSource};
    fOppSignAngle = new THnSparseF("OppSignAngleMaps","OppSignAngleMaps",nDimAngleOppSign,nBinAngleOppSign);
    fOppSignAngle->SetBinEdges(0,binLimAngle);
    fOppSignAngle->SetBinEdges(1,binLimC);
    fOppSignAngle->SetBinEdges(2,binLimSource);
    fOppSignAngle->Sumw2();
    AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: oppsignangle");

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
  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: add default");

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
  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: add monitor");

  if(fMonitorTrackCuts) fListHist->Add(fTrackingCuts);

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: add monitortrackcuts");

  if(fMonitorContamination) {
    fListHist->Add(fDeltaPhiMapsContamination);
  }
  
  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: add deltaphimapscontamination");

  if(fMonitorWithoutPID) {
    fListHist->Add(fDeltaPhiMapsBeforePID);
    fListHist->Add(fCosPhiMapsBeforePID);
  }

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: add without pid");

  if(fMonitorPhotonic) {
  fListHist->Add(fPIDBackgroundqa->MakeList("HFEpidBackgroundQA"));
  fListHist->Add(fDeltaPhiMapsTaggedPhotonic);
  //fListHist->Add(fCosPhiMapsTaggedPhotonic);
  fListHist->Add(fDeltaPhiMapsTaggedNonPhotonic);
  //fListHist->Add(fCosPhiMapsTaggedNonPhotonic);
  fListHist->Add(fDeltaPhiMapsTaggedPhotonicLS);
  //fListHist->Add(fCosPhiMapsTaggedPhotonicLS);
  fListHist->Add(fMCSourceDeltaPhiMaps);
  fListHist->Add(fOppSignDeltaPhiMaps);
  fListHist->Add(fSameSignDeltaPhiMaps);
  fListHist->Add(fSameSignAngle);
  fListHist->Add(fOppSignAngle);
  fListHist->Add(fBackgroundSubtraction->GetListOutput());
  }

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: add photonic");

  if(fHFEVZEROEventPlane && fMonitorEventPlane) fListHist->Add(fHFEVZEROEventPlane->GetOutputList());
  
  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: add event plane");

  PostData(1, fListHist);
  //for(Int_t bincless = 0; bincless < fNbBinsCentralityQCumulant; bincless++) {
  // PostData(bincless+2,fflowEvent); 
  //}

  AliDebug(2,"AliAnalysisTaskFlowTPCTOFEPSP: post");


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

  // MC info
  Bool_t mcthere = kTRUE;
  if(fAODAnalysis) {
    AliAODEvent *aodE = dynamic_cast<AliAODEvent *>(fInputEvent);
    if(!aodE){
      //        printf("testd\n");
      AliError("No AOD Event");
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
      if(fMonitorPhotonic) fBackgroundSubtraction->SetAODArrayMCInfo(fAODArrayMCInfo);
    }
  }
  else {
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH) mcthere = kFALSE;
    else {
      if(fMonitorPhotonic) fBackgroundSubtraction->SetMCEvent(fMCEvent);
    }
  }

    
  /////////////////
  // centrality
  /////////////////
  
  //AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  //if(!esd) return;
  AliCentrality *centrality = fInputEvent->GetCentrality();
  //AliDebug(2,"Got the centrality");
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
  AliDebug(2,Form("Run number %d",runnumber));
   
  if(!fPID->IsInitialized()){
    // Initialize PID with the given run number
    fPID->InitializePID(runnumber);
  }
  if(!fPIDTOFOnly->IsInitialized()){
    // Initialize PID with the given run number
    fPIDTOFOnly->InitializePID(runnumber);
  }

  //
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
    AliDebug(2,"No PID response set");
    return;
  }
  fPID->SetPIDResponse(pidResponse);
  fPIDTOFOnly->SetPIDResponse(pidResponse);
  fPIDBackground->SetPIDResponse(pidResponse);
  if(fMonitorPhotonic) fBackgroundSubtraction->InitRun(fInputEvent,pidResponse);

  fHistEV->Fill(binctt,0.0);

  //////////////////
  // Event cut
  //////////////////
  if(!fHFECuts->CheckEventCuts("fEvRecCuts", fInputEvent)) {
    AliDebug(2,"Does not pass the event cut");
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
  pileup[0]=fInputEvent->GetCentrality()->GetCentralityPercentile("V0M");
  pileup[1]=fInputEvent->GetCentrality()->GetCentralityPercentile("TRK");
  pileup[2]=multTPC;
  pileup[3]=multGlob;
  fHistPileUp->Fill(pileup);

  if(fPileUpCut){
    if (TMath::Abs(pileup[0]-pileup[1]) > 5) {
      AliDebug(2,"Does not pass the centrality correlation cut");
      return;
    }
    if(multTPC < (-36.81+1.48*multGlob) && multTPC > (63.03+1.78*multGlob)){
      AliDebug(2,"Does not pass the multiplicity correlation cut");
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
  //     AliDebug(2,"Does not pass the pileup cut");
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
    AliDebug(2,"No event plane");
    return;
  }

  eventPlanesub1a = TVector2::Phi_0_2pi(qsub1a->Phi())/2.;
  eventPlanesub2a = TVector2::Phi_0_2pi(qsub2a->Phi())/2.;
  diffsub1sub2a = TMath::Cos(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));
  diffsub1sub2asin = TMath::Sin(2.*TVector2::Phi_0_2pi(qsub1a->Phi()/2.- qsub2a->Phi()/2.));


  // if ( !fDebugStreamer ) {
  //   //debug stream
  //   TDirectory *backup = gDirectory;
  //   fDebugStreamer = new TTreeSRedirector("TaskFlowTPCTOFEPSPdebug.root");
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
  AliDebug(2,Form("Number of tracks %d",nbtracks));

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
      AliDebug(2,Form("MC reaction plane %f",mcReactionPlane));
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
  //printf("%f %f %f\n",valuensparsea[0],valuensparsea[1],valuensparsea[2]);
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
	AliDebug(2,Form("centrality %f and %d",binct,hfetrack2.GetCentrality()));
	hfetrack2.SetPbPb();
	if(fPIDBackground->IsSelected(&hfetrack2,0x0,"recTrackCont",fPIDBackgroundqa)) {
	  fArraytrack->AddAt(k,fCounterPoolBackground);
	  fCounterPoolBackground++;
	  AliDebug(2,Form("fCounterPoolBackground %d, track %d",fCounterPoolBackground,k));
	}
      }
    }
  }


  ///////////////////////////////////
  // Loop for kink mother AOD
  //////////////////////////////////
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
	AliDebug(2,Form("ID %d",idmother));
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
      AliDebug(2,"Find AOD track on");
      if(fUseFilterAOD){
	if(!(aodtrack->TestFilterBit(fFilter))) continue;  // Only process AOD tracks where the HFE is set
      }
    }

    if(fApplyCut) {

      valuetrackingcuts[0] = track->Pt(); 
      valuetrackingcuts[1] = 0;
 
      // RecKine: ITSTPC cuts  
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);

      // Reject kink mother
      if(fRejectKinkMother) {
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
      }      

      valuetrackingcuts[1] = 1; 
      if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      // RecPrim
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      valuetrackingcuts[1] = 2; 
      if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      
      // HFEcuts: ITS layers cuts
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
       valuetrackingcuts[1] = 3; 
       if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);     

      // HFE cuts: TOF PID and mismatch flag
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTOF + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      valuetrackingcuts[1] = 4; 
      if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      
      // HFE cuts: TPC PID cleanup
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTPC + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      valuetrackingcuts[1] = 5; 
      if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      
      // HFEcuts: Nb of tracklets TRD0
      if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
      valuetrackingcuts[1] = 6; 
      if(fMonitorTrackCuts) fTrackingCuts->Fill(&valuetrackingcuts[0]);    
      
    }
    AliDebug(2,"Survived");

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

    AliDebug(2,Form("charge %d",(Int_t)track->Charge()));

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
	if(!fAODAnalysis) hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
	else hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
	hfetrack.SetRecTrack(track);
	hfetrack.SetCentrality((Int_t)binct);
	AliDebug(2,Form("centrality %f and %d",binct,hfetrack.GetCentrality()));
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
      }
      else {
	if(!mcEvent) continue;
	if(!(mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel()))))) continue;
	AliDebug(2,Form("PdgCode %d",TMath::Abs(mctrack->Particle()->GetPdgCode())));
	if(TMath::Abs(mctrack->Particle()->GetPdgCode())!=11) continue;
      }
    }


    /////////////////////////////////////////////////////////////////////////////
    // Add candidate to AliFlowEvent for POI and subtract from RP if needed
    ////////////////////////////////////////////////////////////////////////////
    if(fMonitorQCumulant) {
      Int_t idtrack = static_cast<AliVTrack*>(track)->GetID();
      Bool_t found = kFALSE;
      Int_t numberoffound = 0;
      AliDebug(2,Form("A: Number of tracks %d",fflowEvent->NumberOfTracks()));
      for(Int_t iRPs=0; iRPs< fflowEvent->NumberOfTracks(); iRPs++) {
	AliFlowTrack *iRP = (AliFlowTrack*) (fflowEvent->GetTrack(iRPs));
	//if(!iRP->InRPSelection()) continue;
	if( TMath::Abs(idtrack) == TMath::Abs(iRP->GetID()) ) {
	  iRP->SetForPOISelection(kTRUE);
	  found = kTRUE;
	  numberoffound ++;
	}
      }
      AliDebug(2,Form("Found %d mal",numberoffound));
      if(!found) {
	AliFlowCandidateTrack *sTrack = (AliFlowCandidateTrack*) MakeTrack(massElectron,track->Pt(),track->Phi(), track->Eta());
	sTrack->SetID(idtrack);
	fflowEvent->AddTrack(sTrack);
	AliDebug(2,"Add the track");
      }
      AliDebug(2,Form("B: Number of tracks %d",fflowEvent->NumberOfTracks()));
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
      if((valuefractioncont[1] >=0) && (valuefractioncont[1] < 11)){
	if(fContamination[((Int_t)valuefractioncont[1])]){
	  Double_t weight = fContamination[((Int_t)valuefractioncont[1])]->Eval(track->P());
	  if(weight<0.0) weight=0.0;
	  if(weight>1.0) weight=1.0;
	  fFractionContamination->Fill(&valuefractioncont[0],weight);
	  if(fv2contamination[((Int_t)valuefractioncont[1])]){
	    Double_t v2 =  fv2contamination[((Int_t)valuefractioncont[1])]->Eval(track->Pt());
	    AliDebug(2,Form("value v2 %f, contamination %f and pt %f centrality %d\n",v2,weight,track->Pt(),(Int_t)valuefractioncont[1]));
	    AliDebug(2,Form("Check for centrality 3: value v2 %f, contamination %f\n",fv2contamination[3]->Eval(track->Pt()),fContamination[3]->Eval(track->P())));
	    AliDebug(2,Form("Check for centrality 4: value v2 %f, contamination %f\n",fv2contamination[4]->Eval(track->Pt()),fContamination[4]->Eval(track->P())));
	    AliDebug(2,Form("Check for centrality 5: value v2 %f, contamination %f\n",fv2contamination[5]->Eval(track->Pt()),fContamination[5]->Eval(track->P())));
	    fContaminationv2->Fill(valuefractioncont[1],valuefractioncont[0],v2,weight);
	  }
	}     
      }
      if(fMonitorEventPlane) {
	if(fSP)
	  valuensparseh[0] *= TMath::Sqrt(qAna->X()*qAna->X()+qAna->Y()*qAna->Y());
	fProfileCosPhiMaps->Fill(valuensparsehprofile[1],valuensparsehprofile[2],valuensparseh[0]);  //TR: info
      }
    }
    
    if(fMonitorPhotonic) {
      Int_t indexmother = -1;
      Int_t source = 1;
      if(mcthere) source = fBackgroundSubtraction->FindMother(mctrack->GetLabel(),indexmother);
      fBackgroundSubtraction->LookAtNonHFE(k, track, fInputEvent, 1, binct, deltaphi, source, indexmother);
      
      if((!fAODAnalysis && mcthere) || !mcthere) {
	// background
	source = 0;
	indexmother = -1;
	source = FindMother(TMath::Abs(track->GetLabel()),mcEvent, indexmother);
	valuensparseMCSourceDeltaPhiMaps[2] = source;
	if(mcEvent) fMCSourceDeltaPhiMaps->Fill(&valuensparseMCSourceDeltaPhiMaps[0]);
	//LookAtNonHFE(k,track,fInputEvent,mcEvent,binct,deltaphi,source,indexmother);
	Int_t taggedvalue = LookAtNonHFE(k,track,fInputEvent,mcEvent,binct,deltaphi,source,indexmother);
	if(fMonitorPhotonic) {
	  // No opposite charge partner found in the invariant mass choosen
	  if((taggedvalue!=2) && (taggedvalue!=6)) {
	    //fDeltaPhiMapsTaggedNonPhotonic->Fill(&valuensparseg[0]);
	    //fCosPhiMapsTaggedNonPhotonic->Fill(&valuensparseh[0]);
	  }
	  // One opposite charge partner found in the invariant mass choosen
	  if((taggedvalue==2) || (taggedvalue==6)) {
	    fDeltaPhiMapsTaggedPhotonic->Fill(&valuensparseg[0]);
	    //fCosPhiMapsTaggedPhotonic->Fill(&valuensparseh[0]);
	  }
	  // One same charge partner found in the invariant mass choosen
	  if((taggedvalue==4) || (taggedvalue==6)) {
	    fDeltaPhiMapsTaggedPhotonicLS->Fill(&valuensparseg[0]);
	    //fCosPhiMapsTaggedPhotonicLS->Fill(&valuensparseh[0]);
	  }
	}
      }
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

  if(fMonitorPhotonic) {
    if(fArraytrack) {
      delete fArraytrack;
      fArraytrack = NULL;
    }
  }
  
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
//_____________________________________________________________________________________________
Int_t AliAnalysisTaskFlowTPCTOFEPSP::LookAtNonHFE(Int_t iTrack1, AliVTrack *track1, AliVEvent *vEvent, AliMCEvent *mcEvent,Int_t binct,Double_t deltaphi,Int_t source,Int_t indexmother)
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

  AliDebug(2,Form("fCounterPoolBackground %d in LookAtNonHFE!!!",fCounterPoolBackground));
  if(!fArraytrack) return taggedphotonic;
  AliDebug(2,Form("process track %d",iTrack1));
  
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

  // Pdg code
  Int_t pdg1 = CheckPdg(TMath::Abs(track1->GetLabel()),mcEvent);
  Int_t numberfound = 0;

  //Magnetic Field
  Double_t bfield = vEvent->GetMagneticField();

  // Get Primary vertex
  const AliVVertex *pVtx = vEvent->GetPrimaryVertex();
  
  for(Int_t idex = 0; idex < fCounterPoolBackground; idex++) 
    {

      Int_t iTrack2 = fArraytrack->At(idex);
      AliDebug(2,Form("track %d",iTrack2));
      AliVTrack* track2 = (AliVTrack *) vEvent->GetTrack(iTrack2);
      if (!track2) 
	{
	  printf("ERROR: Could not receive track %d", iTrack2);
	  continue;
	}
      if(iTrack2==iTrack1) continue;
      AliDebug(2,"Different");

      // Reset the MC info
      valueangle[2] = source;
      valuensparseDeltaPhiMaps[4] = source;

      // track cuts and PID already done

      // if MC look
      Int_t pdg2 = -100;
      if(mcEvent) {
	Int_t source2 = 0;
	Int_t indexmother2 = -1;
	source2 = FindMother(TMath::Abs(track2->GetLabel()),mcEvent, indexmother2);
	pdg2 = CheckPdg(TMath::Abs(track2->GetLabel()),mcEvent);
	if(source2 >=0 ) {
	  if((indexmother2 == indexmother) && (source == source2) && ((pdg1*pdg2)<0.0)) {
	    if(source == kElectronfromconversion) {
	      valueangle[2] = kElectronfromconversionboth;
	      valuensparseDeltaPhiMaps[4] = kElectronfromconversionboth;
	      numberfound++;
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

	  // DCA
	  //Double_t dca12 = ktrack1.GetDistanceFromParticle(ktrack2);
	  //if(dca12 > fMaxdca) continue;	  

	  // if set mass constraint
	  if(fSetMassConstraint && pVtx) {
	    AliKFVertex primV(*pVtx);
	    primV += recoGamma;
	    primV -= ktrack1;
	    primV -= ktrack2;
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
	  else {
	    fOppSignDeltaPhiMaps->Fill(&valuensparseDeltaPhiMaps[0]);
	    /*
	    if(valueangle[2] == kElectronfromconversionboth) {
	      printf("Reconstructed charge1 %f, charge 2 %f and invmass %f",fCharge1,fCharge2,imass);
	      printf("MC charge1 %d, charge 2 %d",pdg1,pdg2);
	      printf("DCA %f",dca12);
	      printf("Number of found %d",numberfound);
	    }
	    */
	  }
	  
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
Int_t AliAnalysisTaskFlowTPCTOFEPSP::FindMother(Int_t tr, AliMCEvent *mcEvent, Int_t &indexmother){
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
Int_t AliAnalysisTaskFlowTPCTOFEPSP::CheckPdg(Int_t tr, AliMCEvent* mcEvent) {

  //
  // Return the pdg of the particle
  //


  Int_t pdgcode = -1;
  if(tr < 0) return pdgcode;

  if(!mcEvent) return pdgcode;

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
Int_t AliAnalysisTaskFlowTPCTOFEPSP::IsMotherGamma(Int_t tr, AliMCEvent* mcEvent) {

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
Int_t AliAnalysisTaskFlowTPCTOFEPSP::IsMotherPi0(Int_t tr, AliMCEvent* mcEvent) {

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
Int_t AliAnalysisTaskFlowTPCTOFEPSP::IsMotherC(Int_t tr, AliMCEvent* mcEvent) {

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
Int_t AliAnalysisTaskFlowTPCTOFEPSP::IsMotherB(Int_t tr, AliMCEvent* mcEvent) {

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
Int_t AliAnalysisTaskFlowTPCTOFEPSP::IsMotherEta(Int_t tr, AliMCEvent* mcEvent) {

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
