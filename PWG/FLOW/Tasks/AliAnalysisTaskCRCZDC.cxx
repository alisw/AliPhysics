/**************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in thce supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/**********************************
* analysis task for CRC with ZDC *
*                                *
* author: Jacopo Margutti        *
*         (margutti@nikhef.nl)   *
**********************************/

#include "Riostream.h" //needed as include
#include "TChain.h"
#include "TList.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TTimeStamp.h"
#include "TStopwatch.h"
#include "TProfile.h"
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3D.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TParticle.h>
#include "TProfile2D.h"
#include "TProfile3D.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskCRCZDC.h"
#include "AliMultSelection.h"
#include "AliLumiTools.h"

// ALICE Correction Framework
#include "AliCFManager.h"

// Interface to Event generators to get Reaction Plane Angle
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenGeVSimEventHeader.h"
#include "AliGenEposEventHeader.h"

// Interface to Load short life particles
#include "TObjArray.h"
#include "AliFlowCandidateTrack.h"

// Interface to make the Flow Event Simple used in the flow analysis methods
#include "AliFlowEvent.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliFlowCommonConstants.h"

//NanoAODs
#include "AliNanoAODTrack.h"
#include "AliNanoAODHeader.h"

ClassImp(AliAnalysisTaskCRCZDC)

//________________________________________________________________________
AliAnalysisTaskCRCZDC::AliAnalysisTaskCRCZDC():
AliAnalysisTaskSE(""),
fAnalysisType(kAUTOMATIC),
fRPType(""),
fCFManager1(NULL),
fCFManager2(NULL),
fCutsEvent(NULL),
fCutsRP(NULL),
fCutsPOI(NULL),
fCutContainer(new TList()),
fQAList(NULL),
fMinMult(0),
fMaxMult(10000000),
fAnalysisUtil(NULL),
fMinA(-1.0),
fMaxA(-0.01),
fMinB(0.01),
fMaxB(1.0),
fGenHeader(NULL),
fPythiaGenHeader(NULL),
fHijingGenHeader(NULL),
fFlowTrack(NULL),
fQAon(kFALSE),
fLoadCandidates(kFALSE),
fNbinsMult(10000),
fNbinsPt(100),
fNbinsPhi(100),
fNbinsEta(200),
fNbinsQ(500),
fNbinsMass(1),
fMultMin(0.),
fMultMax(10000.),
fPtMin(0.),
fPtMax(10.),
fPhiMin(0.),
fPhiMax(TMath::TwoPi()),
fEtaMin(-5.),
fEtaMax(5.),
fQMin(0.),
fQMax(3.),
fMassMin(-1.),
fMassMax(0.),
fHistWeightvsPhiMin(0.),
fHistWeightvsPhiMax(3.),
fExcludedEtaMin(0.),
fExcludedEtaMax(0.),
fExcludedPhiMin(0.),
fExcludedPhiMax(0.),
fAfterburnerOn(kFALSE),
fNonFlowNumberOfTrackClones(0),
fV1(0.),
fV2(0.),
fV3(0.),
fV4(0.),
fV5(0.),
fDifferentialV2(0),
fFlowEvent(NULL),
fShuffleTracks(kFALSE),
fMyTRandom3(NULL),
fAnalysisInput(kAOD),
fIsMCInput(kFALSE),
fUseMCCen(kTRUE),
fRejectPileUp(kTRUE),
fRejectPileUpTight(kFALSE),
fResetNegativeZDC(kFALSE),
fCentrLowLim(0.),
fCentrUpLim(100.),
fCentrEstimator(kV0M),
fOutput(0x0),
fOutputRecenter1(0x0),
fOutputRecenter2(0x0),
fOutputRecenter3(0x0),
fhZNCvsZNA(0x0),
fhZDCCvsZDCCA(0x0),
fhZNCvsZPC(0x0),
fhZNAvsZPA(0x0),
fhZNvsZP(0x0),
fhZNvsVZERO(0x0),
fhZDCvsVZERO(0x0),
fhAsymm(0x0),
fhZNAvsAsymm(0x0),
fhZNCvsAsymm(0x0),
fhZNCvscentrality(0x0),
fhZNAvscentrality(0x0),
fhZPCvscentrality(0x0),
fhZPAvscentrality(0x0),
fZPAvsZNASignal(0x0),
fZPCvsZNCSignal(0x0),
fZNenergyBeforeCalibration(0x0),
fZNenergyAfterCalibration(0x0),
fCRCnRun(0),
fZDCGainAlpha(0.395),
fDataSet(kAny),
fStack(0x0),
fCutTPC(kFALSE),
fCenDis(0x0),
fVZEROMult(0x0),
fMultSelection(0x0),
fPileUpCount(0x0),
fPileUpMultSelCount(0x0),
fMultTOFLowCut(0x0),
fMultTOFHighCut(0x0),
fUseTowerEq(kFALSE),
fFillZNCenDisRbR(kFALSE),
fTowerEqList(NULL),
fZDCCalibList(NULL),
fZDCCalibListStep3CommonPart(NULL),
fZDCCalibListStep3RunByRun(NULL),
fUseBadTowerCalib(kFALSE),
fBadTowerCalibList(NULL),
fVZEROGainEqList(NULL),
fVZEROQVecRecList(NULL),
fUseZDCSpectraCorr(kFALSE),
fCorrectPhiTracklets(kFALSE),
fZDCSpectraCorrList(NULL),
fSpectraMCList(NULL),
fTrackQAList(NULL),
fBadTowerStuffList(NULL),
fVZEROStuffList(NULL),
fVZEROGainEqHist(NULL),
fMinRingVZC(1),
fMaxRingVZC(4),
fMinRingVZA(5),
fMaxRingVZA(8),
fCachedRunNum(0),
fhZNSpectra(0x0),
fhZNSpectraCor(0x0),
fhZNSpectraPow(0x0),
fhZPSpectra(0x0),
fhZNBCCorr(0x0),
fQATrackTPCNcls(NULL),
fQATrackITSNcls(NULL),
fQATrackTPCchi2(NULL),
fQATrackITSchi2(NULL),
fQATrackTPCScls(NULL),
fQATrackITSScls(NULL)
{ 
  for(int i=0; i<5; i++){
    fhZNCPM[i] = 0x0;
    fhZNAPM[i] = 0x0;
    //@Shi add fhZPCPM and fhZPAPM
    fhZPCPM[i] = 0x0;
    fhZPAPM[i] = 0x0;
  }
  for(int i=0; i<4; i++){
    fhZNCPMQiPMC[i] = 0x0;
    fhZNAPMQiPMC[i] = 0x0;
    //@Shi add fhZPCPMQiPMC fhZPAPMQiPMC
    fhZPCPMQiPMC[i] = 0x0;
    fhZPAPMQiPMC[i] = 0x0;
  }
  for(Int_t r=0; r<fCRCMaxnRun; r++) {
    fRunList[r] = 0;
  }
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] =  NULL;
    }
  }
  for(Int_t c=0; c<100; c++) {
    fBadTowerCalibHist[c] = NULL;
  }
  for (Int_t k=0; k<fkVZEROnHar; k++) {
    //    fVZEROQVectorRecQx[k] = NULL;
    //    fVZEROQVectorRecQy[k] = NULL;
    fVZEROQVectorRecQxStored[k] = NULL;
    fVZEROQVectorRecQyStored[k] = NULL;
    for (Int_t t=0; t<fkVZEROnQAplots; t++) {
      fVZEROQVectorRecFinal[k][t] = NULL;
    }
  }
  for(Int_t i=0; i<8; i++) {
    SpecCorMu1[i] = NULL;
    SpecCorMu2[i] = NULL;
    SpecCorSi[i] = NULL;
    SpecCorAv[i] = NULL;
  }
  //@shi initialize histograms for recentering ZDC (begin)
  for(Int_t i=0; i<4; i++) {
    fAvr_Run_CentQ[i] = NULL;
    fAvr_Run_VtxXYZQ[i] = NULL;
    fAvr_Cent_VtxXYZQ[i] = NULL;
  }  
  //@shi initialize histograms for recentering ZDC (end)
  this->InitializeRunArrays();
  fMyTRandom3 = new TRandom3(1);
  gRandom->SetSeed(fMyTRandom3->Integer(65539));
  for(Int_t j=0; j<2; j++) {
    for(Int_t c=0; c<10; c++) {
      fPtSpecGen[j][c] = NULL;
      fPtSpecFB32[j][c] = NULL;
      fPtSpecFB96[j][c] = NULL;
      fPtSpecFB128[j][c] = NULL;
      fPtSpecFB768[j][c] = NULL;
    }
  }
  for (Int_t c=0; c<2; c++) {
    fhZNCenDis[c] = NULL;
  }
  for(Int_t fb=0; fb<fKNFBs; fb++){
    for (Int_t i=0; i<4; i++) {
      fTrackQADCAxy[fb][i] = NULL;
      fTrackQADCAz[fb][i] = NULL;
      fTrackQApT[fb][i] = NULL;
      fTrackQADphi[fb][i] = NULL;
      fEbEQRe[fb][i] = NULL;
      fEbEQIm[fb][i] = NULL;
      fEbEQMu[fb][i] = NULL;
    }
  }
}

//________________________________________________________________________
AliAnalysisTaskCRCZDC::AliAnalysisTaskCRCZDC(const char *name, TString RPtype, Bool_t on, UInt_t iseed, Bool_t bCandidates, Int_t StepZDCRecenter):
AliAnalysisTaskSE(name),
fAnalysisType(kAUTOMATIC),
fRPType(RPtype),
fCFManager1(NULL),
fCFManager2(NULL),
fCutsEvent(NULL),
fCutsRP(NULL),
fCutsPOI(NULL),
fCutContainer(new TList()),
fQAList(NULL),
fMinMult(0),
fMaxMult(10000000),
fAnalysisUtil(NULL),
fMinA(-1.0),
fMaxA(-0.01),
fMinB(0.01),
fMaxB(1.0),
fQAon(on),
fLoadCandidates(bCandidates),
fNbinsMult(10000),
fNbinsPt(100),
fNbinsPhi(100),
fNbinsEta(200),
fNbinsQ(500),
fNbinsMass(1),
fMultMin(0.),
fMultMax(10000.),
fPtMin(0.),
fPtMax(10.),
fPhiMin(0.),
fPhiMax(TMath::TwoPi()),
fEtaMin(-5.),
fEtaMax(5.),
fQMin(0.),
fQMax(3.),
fMassMin(-1.),
fMassMax(0.),
fHistWeightvsPhiMin(0.),
fHistWeightvsPhiMax(3.),
fExcludedEtaMin(0.),
fExcludedEtaMax(0.),
fExcludedPhiMin(0.),
fExcludedPhiMax(0.),
fAfterburnerOn(kFALSE),
fNonFlowNumberOfTrackClones(0),
fV1(0.),
fV2(0.),
fV3(0.),
fV4(0.),
fV5(0.),
fDifferentialV2(0),
fFlowEvent(NULL),
fShuffleTracks(kFALSE),
fMyTRandom3(NULL),
fAnalysisInput(kAOD),
fIsMCInput(kFALSE),
fUseMCCen(kTRUE),
fRejectPileUp(kTRUE),
fRejectPileUpTight(kFALSE),
fResetNegativeZDC(kFALSE),
fCentrLowLim(0.),
fCentrUpLim(100.),
fCentrEstimator(kV0M),
fOutput(0x0),
fOutputRecenter1(0x0),
fOutputRecenter2(0x0),
fOutputRecenter3(0x0),
fhZNCvsZNA(0x0),
fhZDCCvsZDCCA(0x0),
fhZNCvsZPC(0x0),
fhZNAvsZPA(0x0),
fhZNvsZP(0x0),
fhZNvsVZERO(0x0),
fhZDCvsVZERO(0x0),
fhAsymm(0x0),
fhZNAvsAsymm(0x0),
fhZNCvsAsymm(0x0),
fhZNCvscentrality(0x0),
fhZNAvscentrality(0x0),
fhZPCvscentrality(0x0),
fhZPAvscentrality(0x0),
fZPAvsZNASignal(0x0),
fZPCvsZNCSignal(0x0),
fZNenergyBeforeCalibration(0x0),
fZNenergyAfterCalibration(0x0),
fDataSet(kAny),
fCRCnRun(0),
fZDCGainAlpha(0.395),
fGenHeader(NULL),
fPythiaGenHeader(NULL),
fHijingGenHeader(NULL),
fFlowTrack(NULL),
fStack(0x0),
fCutTPC(kFALSE),
fCenDis(0x0),
fVZEROMult(0x0),
fMultSelection(0x0),
fPileUpCount(0x0),
fPileUpMultSelCount(0x0),
fMultTOFLowCut(0x0),
fMultTOFHighCut(0x0),
fUseTowerEq(kFALSE),
fFillZNCenDisRbR(kFALSE),
fTowerEqList(NULL),
fZDCCalibList(NULL),
fZDCCalibListStep3CommonPart(NULL),
fZDCCalibListStep3RunByRun(NULL),
fUseBadTowerCalib(kFALSE),
fBadTowerCalibList(NULL),
fVZEROGainEqList(NULL),
fVZEROQVecRecList(NULL),
fUseZDCSpectraCorr(kFALSE),
fCorrectPhiTracklets(kFALSE),
fZDCSpectraCorrList(NULL),
fSpectraMCList(NULL),
fTrackQAList(NULL),
fBadTowerStuffList(NULL),
fVZEROStuffList(NULL),
fVZEROGainEqHist(NULL),
fMinRingVZC(1),
fMaxRingVZC(4),
fMinRingVZA(5),
fMaxRingVZA(8),
fCachedRunNum(0),
fhZNSpectra(0x0),
fhZNSpectraCor(0x0),
fhZNSpectraPow(0x0),
fhZPSpectra(0x0),
fhZNBCCorr(0x0),
fQATrackTPCNcls(NULL),
fQATrackITSNcls(NULL),
fQATrackTPCchi2(NULL),
fQATrackITSchi2(NULL),
fQATrackTPCScls(NULL),
fQATrackITSScls(NULL)
{
  for(int i=0; i<5; i++){
    fhZNCPM[i] = 0x0;
    fhZNAPM[i] = 0x0;
    //@Shi add fhZPCPM and fhZPAPM
    fhZPCPM[i] = 0x0;
    fhZPAPM[i] = 0x0;
  }
  for(int i=0; i<4; i++){
    fhZNCPMQiPMC[i] = 0x0;
    fhZNAPMQiPMC[i] = 0x0;
    //@Shi add fhZPCPMQiPMC fhZPAPMQiPMC
    fhZPCPMQiPMC[i] = 0x0;
    fhZPAPMQiPMC[i] = 0x0;
  }
  for(Int_t r=0; r<fCRCMaxnRun; r++) {
    fRunList[r] = 0;
  }
  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] =  NULL;
    }
  }
  for(Int_t c=0; c<100; c++) {
    fBadTowerCalibHist[c] = NULL;
  }
  for (Int_t k=0; k<fkVZEROnHar; k++) {
    //    fVZEROQVectorRecQx[k] = NULL;
    //    fVZEROQVectorRecQy[k] = NULL;
    fVZEROQVectorRecQxStored[k] = NULL;
    fVZEROQVectorRecQyStored[k] = NULL;
    for (Int_t t=0; t<fkVZEROnQAplots; t++) {
      fVZEROQVectorRecFinal[k][t] = NULL;
    }
  }
  for(Int_t i=0; i<8; i++) {
    SpecCorMu1[i] = NULL;
    SpecCorMu2[i] = NULL;
    SpecCorSi[i] = NULL;
    SpecCorAv[i] = NULL;
  }
  //@shi initialize histograms for recentering ZDC (begin)
  for(Int_t i=0; i<4; i++) {
    fAvr_Run_CentQ[i] = NULL;
    fAvr_Run_VtxXYZQ[i] = NULL;
    fAvr_Cent_VtxXYZQ[i] = NULL;
  }  
  //@shi initialize histograms for recentering ZDC (end)
  this->InitializeRunArrays();
  fMyTRandom3 = new TRandom3(iseed);
  gRandom->SetSeed(fMyTRandom3->Integer(65539));

  DefineInput(0, TChain::Class());
  // Define output slots here
  // Define here the flow event output
  DefineOutput(1, AliFlowEventSimple::Class());
  DefineOutput(2, TList::Class());
  
  if (StepZDCRecenter >= 0) {
	DefineOutput(3, TList::Class());
	DefineOutput(4, TList::Class());
	DefineOutput(5, TList::Class());
  }

  for(Int_t j=0; j<2; j++) {
    for(Int_t c=0; c<10; c++) {
      fPtSpecGen[j][c] = NULL;
      fPtSpecFB32[j][c] = NULL;
      fPtSpecFB96[j][c] = NULL;
      fPtSpecFB128[j][c] = NULL;
      fPtSpecFB768[j][c] = NULL;
    }
  }
  for (Int_t c=0; c<2; c++) {
    fhZNCenDis[c] = NULL;
  }
  for(Int_t fb=0; fb<fKNFBs; fb++){
    for (Int_t i=0; i<4; i++) {
      fTrackQADCAxy[fb][i] = NULL;
      fTrackQADCAz[fb][i] = NULL;
      fTrackQApT[fb][i] = NULL;
      fTrackQADphi[fb][i] = NULL;
      fEbEQRe[fb][i] = NULL;
      fEbEQIm[fb][i] = NULL;
      fEbEQMu[fb][i] = NULL;
    }
  }
}

//________________________________________________________________________
AliAnalysisTaskCRCZDC::~AliAnalysisTaskCRCZDC()
{
  // Destructor
  if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fOutput; fOutput=0;
  }
  //@Shi add destructor for fOutputRecenter1 and fOutputRecenter2 and fOutputRecenter3
  if(fStepZDCRecenter >= 0) {
    if(fOutputRecenter1 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
      delete fOutputRecenter1; fOutputRecenter1=0;
    }
    if(fOutputRecenter2 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
      delete fOutputRecenter2; fOutputRecenter2=0;
    }
    if(fOutputRecenter3 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
      delete fOutputRecenter3; fOutputRecenter3=0;
    }
  }
  delete fMyTRandom3;
  delete fFlowEvent;
  delete fFlowTrack;
  delete fCutsEvent;
  if (fTowerEqList)  delete fTowerEqList;
  if (fBadTowerCalibList) delete fBadTowerCalibList;
  if (fVZEROGainEqList) delete fVZEROGainEqList;
  if (fVZEROQVecRecList) delete fVZEROQVecRecList;
  if (fZDCSpectraCorrList) delete fZDCSpectraCorrList;
  if (fAnalysisUtil) delete fAnalysisUtil;
  if (fQAList) delete fQAList;
  if (fCutContainer) fCutContainer->Delete(); delete fCutContainer;
  if (fZDCCalibList) delete fZDCCalibList; //@shi calibration file for ZDC recentering
  if (fZDCCalibListStep3CommonPart) delete fZDCCalibListStep3CommonPart; 
  if (fZDCCalibListStep3RunByRun) delete fZDCCalibListStep3RunByRun; 
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::InitializeRunArrays()
{
  for(Int_t r=0;r<fCRCMaxnRun;r++) {
    fCRCQVecListRun[r] = NULL;
    //@Shi Add initializing run by run ZN centroid vs centrality
    fRecenter1ListRunbyRun[r] = NULL; 
    fRecenter2ListRunbyRun[r] = NULL; 
    fRecenter3ListRunbyRun[r] = NULL; 
    for (Int_t c=0; c<2; c++) {
      fhZNCenDisRbR[r][c] = NULL;
    }
    for(Int_t k=0;k<fCRCnTow;k++) {
      fZNCTower[r][k] = NULL;
      fZNATower[r][k] = NULL;
      fZPCTower[r][k] = NULL;
      fZPATower[r][k] = NULL;
    }
    //    fhZNSpectraRbR[r] = NULL;
    
    for (Int_t c=0; c<4; c++) {
		fRun_VtxXQPreCalib[r][c] = NULL;
		fRun_VtxYQPreCalib[r][c] = NULL;
		fRun_VtxZQPreCalib[r][c] = NULL;
		fRun_VtxXQCalibStep1[r][c] = NULL;
		fRun_VtxYQCalibStep1[r][c] = NULL;
		fRun_VtxZQCalibStep1[r][c] = NULL;
		fRun_VtxXQCalibStep2[r][c] = NULL;
		fRun_VtxYQCalibStep2[r][c] = NULL;
		fRun_VtxZQCalibStep2[r][c] = NULL;
		fRun_CentQCalib[r][c] = NULL;
		fRun_VtxXQCalib[r][c] = NULL;
		fRun_VtxYQCalib[r][c] = NULL;
		fRun_VtxZQCalib[r][c] = NULL;
		fRun_CentQCalib2[r][c] = NULL;
		fRun_CentQ[r][c] = NULL;
		fRun_VtxXYZQ[r][c] = NULL;
	}  
  }
  fCorrQAReCRe = NULL;
  fCorrQAReCIm = NULL;
  fCorrQAImCRe = NULL;
  fCorrQAImCIm = NULL;
  for(Int_t r=0;r<fnCentBinForRecentering;r++) {
	for (Int_t c=0; c<4; c++) {
      fCent_VtxXYZQ[r][c] = NULL;
	}
  }
  
  fAve_VtxX = NULL;
  fAve_VtxY = NULL;
  fAve_VtxZ = NULL;
  
  //   for(Int_t i=0;i<fnCen;i++) {
  //     fPtPhiEtaRbRFB128[r][i] = NULL;
  //     fPtPhiEtaRbRFB768[r][i] = NULL;
  //   }
  // }
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::UserCreateOutputObjects()
{
  // Create the output containers
  //set the common constants
  AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(fNbinsMult);
  cc->SetNbinsPt(fNbinsPt);
  cc->SetNbinsPhi(fNbinsPhi);
  cc->SetNbinsEta(fNbinsEta);
  cc->SetNbinsQ(fNbinsQ);
  cc->SetNbinsMass(fNbinsMass);
  cc->SetMultMin(fMultMin);
  cc->SetMultMax(fMultMax);
  cc->SetPtMin(fPtMin);
  cc->SetPtMax(fPtMax);
  cc->SetPhiMin(fPhiMin);
  cc->SetPhiMax(fPhiMax);
  cc->SetEtaMin(fEtaMin);
  cc->SetEtaMax(fEtaMax);
  cc->SetQMin(fQMin);
  cc->SetQMax(fQMax);
  cc->SetMassMin(fMassMin);
  cc->SetMassMax(fMassMax);
  cc->SetHistWeightvsPhiMax(fHistWeightvsPhiMax);
  cc->SetHistWeightvsPhiMin(fHistWeightvsPhiMin);

  fFlowEvent = new AliFlowEvent(20000);
  fFlowTrack = new AliFlowTrack();

  //printf("  AliAnalysisTaskCRCZDC::UserCreateOutputObjects()\n\n");
  fOutput = new TList();
  fOutput->SetOwner(kTRUE);
  //fOutput->SetName("output");
  
  //@Shi add fOutputRecenter1 & fOutputRecenter2 & fOutputRecenter3
  if (fStepZDCRecenter >= 0) {
    fOutputRecenter1 = new TList();
    fOutputRecenter1->SetOwner(kTRUE);
    fOutputRecenter2 = new TList();
    fOutputRecenter2->SetOwner(kTRUE);
    fOutputRecenter3 = new TList();
    fOutputRecenter3->SetOwner(kTRUE);
  }
  if (fQAon) {
    fQAList = new TList();
    fQAList->SetOwner(kTRUE);
    fQAList->SetName("AliFlowEventCuts QA");
    if (fCutsEvent->GetQA()) fQAList->Add(fCutsEvent->GetQA()); //0
    if (fCutsRP->GetQA()) fQAList->Add(fCutsRP->GetQA());       //1
    if (fCutsPOI->GetQA())fQAList->Add(fCutsPOI->GetQA());      //2
    fOutput->Add(fQAList);
  }

  fVZEROStuffList = new TList();
  fVZEROStuffList->SetOwner(kTRUE);
  fVZEROStuffList->SetName("VZERO stuff");
  fOutput->Add(fVZEROStuffList);

  fBadTowerStuffList = new TList();
  fBadTowerStuffList->SetOwner(kTRUE);
  fBadTowerStuffList->SetName("BadTowerCalib");
  fOutput->Add(fBadTowerStuffList);

  fCenDis = new TH1F("fCenDis", "fCenDis", 100, 0., 100.);
  fOutput->Add(fCenDis);
  fPileUpCount = new TH1F("fPileUpCount", "fPileUpCount", 9, 0., 9.);
  fPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultiplicityComb08");
  fPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>7.5");
  fPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(8,"multESDTPCDif");
  fPileUpCount->GetXaxis()->SetBinLabel(9,"extraPileUpMultSel");
  fOutput->Add(fPileUpCount);
  fPileUpMultSelCount = new TH1F("fPileUpMultSelCount", "fPileUpMultSelCount", 8, 0., 8.);
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(1,"IsNotPileup");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(2,"IsNotPileupMV");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(3,"IsNotPileupInMultBins");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(4,"InconsistentVertices");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(5,"TrackletVsCluster");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(6,"AsymmetricInVZERO");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(7,"IncompleteDAQ");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(8,"GoodVertex2016");
  fOutput->Add(fPileUpMultSelCount);
  
  // Add two histograms to record the number of occurance of Negative EZNA and EZNC value (Shi)
  fRecordNegativeEZNA = new TH1F("fRecordNegativeEZNA", "fRecordNegativeEZNA", 2, 0., 2.);
  fRecordNegativeEZNA->GetXaxis()->SetBinLabel(1,"Positive EZNA (okay)");
  fRecordNegativeEZNA->GetXaxis()->SetBinLabel(2,"Negative EZNA (problematic)");
  fOutput->Add(fRecordNegativeEZNA);
  fRecordNegativeEZNC = new TH1F("fRecordNegativeEZNC", "fRecordNegativeEZNC", 2, 0., 2.);
  fRecordNegativeEZNC->GetXaxis()->SetBinLabel(1,"Positive EZNC (okay)");
  fRecordNegativeEZNC->GetXaxis()->SetBinLabel(2,"Negative EZNC (problematic)");
  fOutput->Add(fRecordNegativeEZNC);

  fMultTOFLowCut = new TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
  fMultTOFLowCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
  fOutput->Add(fMultTOFLowCut);
  fMultTOFHighCut = new TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
  fMultTOFHighCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
  fOutput->Add(fMultTOFHighCut);

  for(Int_t c=0; c<2; c++) {
    for(Int_t i=0; i<5; i++) {
      fTowerGainEq[c][i] = new TH1D();
      fOutput->Add(fTowerGainEq[c][i]);
    }
  }

  for (Int_t c=0; c<2; c++) {
    fhZNCenDis[c] = new TH3D(Form("fhZNCenDis[%d]",c), Form("fhZNCenDis[%d]",c), 100, 0., 100., 100, -2., 2. , 100., -2., 2.);
    fOutput->Add(fhZNCenDis[c]);
  }

  if(fBadTowerCalibList) {
    for(Int_t c=0; c<100; c++) {
      fBadTowerCalibHist[c] = (TH2D*)fBadTowerCalibList->FindObject(Form("TH2Resp[%d]",c));
      fBadTowerStuffList->Add(fBadTowerCalibHist[c]);
    }
  }
  if(fZDCSpectraCorrList) {
    for(Int_t i=0; i<8; i++) {
      SpecCorMu1[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorMu1[%d]",i));
      fOutput->Add(SpecCorMu1[i]);
      SpecCorMu2[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorMu2[%d]",i));
      fOutput->Add(SpecCorMu2[i]);
      SpecCorAv[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorAv[%d]",i));
      fOutput->Add(SpecCorAv[i]);
      SpecCorSi[i] = (TH1D*)fZDCSpectraCorrList->FindObject(Form("SpecCorSi[%d]",i));
      fOutput->Add(SpecCorSi[i]);
    }
  }

  fhZNSpectra = new TH3D("fhZNSpectra","fhZNSpectra",100,0.,100.,8,0.,8.,1000,0.,1.E5);
  fOutput->Add(fhZNSpectra);
  fhZNSpectraCor = new TH3D("fhZNSpectraCor","fhZNSpectraCor",100,0.,100.,8,0.,8.,1000,0.,1.E5);
  fOutput->Add(fhZNSpectraCor);
  fhZNSpectraPow = new TH3D("fhZNSpectraPow","fhZNSpectraPow",100,0.,100.,8,0.,8.,1000,0.,TMath::Power(1.E5,fZDCGainAlpha));
  fOutput->Add(fhZNSpectraPow);
  fhZNBCCorr = new TH3D("fhZNBCCorr","fhZNBCCorr",100,0.,100.,500,0.,1.E5,500,0.,1.E5);
  fOutput->Add(fhZNBCCorr);
  
  //@Shi add fhZPSpectra
  fhZPSpectra = new TH3D("fhZPSpectra","fhZPSpectra",100,0.,100.,8,0.,8.,1000,0.,1.E5);
  fOutput->Add(fhZPSpectra);

  fQATrackTPCNcls = new TH3D("fQATrackTPCNcls","fQATrackTPCNcls",50,0.,TMath::TwoPi(),16,-0.8,0.8,50,50.,150.);
  fOutput->Add(fQATrackTPCNcls);
  fQATrackITSNcls = new TH3D("fQATrackITSNcls","fQATrackITSNcls",50,0.,TMath::TwoPi(),16,-0.8,0.8,6,0.,6.);
  fOutput->Add(fQATrackITSNcls);
  fQATrackTPCchi2 = new TH3D("fQATrackTPCchi2","fQATrackTPCchi2",50,0.,TMath::TwoPi(),16,-0.8,0.8,50,0.,5.);
  fOutput->Add(fQATrackTPCchi2);
  fQATrackITSchi2 = new TH3D("fQATrackITSchi2","fQATrackITSchi2",50,0.,TMath::TwoPi(),16,-0.8,0.8,50,0.,50.);
  fOutput->Add(fQATrackITSchi2);
  fQATrackTPCScls = new TH3D("fQATrackTPCScls","fQATrackTPCScls",50,0.,TMath::TwoPi(),16,-0.8,0.8,50,0.,1.);
  fOutput->Add(fQATrackTPCScls);
  fQATrackITSScls = new TH3D("fQATrackITSScls","fQATrackITSScls",50,0.,TMath::TwoPi(),16,-0.8,0.8,50,0.,1.);
  fOutput->Add(fQATrackITSScls);

  if(fAnalysisType == kMCAOD) {

    fSpectraMCList = new TList();
    fSpectraMCList->SetOwner(kTRUE);
    fSpectraMCList->SetName("Spectra");
    fOutput->Add(fSpectraMCList);

    Double_t xmin[] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.33,2.66,3.,3.5,4.,5.,6.,9.,20.};
    for(Int_t j=0; j<2; j++) {
      for(Int_t c=0; c<10; c++) {
        fPtSpecGen[j][c] = new TH1F(Form("fPtSpecGen[%d][%d]",j,c), Form("fPtSpecGen[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecGen[j][c]);
        fPtSpecFB32[j][c] = new TH1F(Form("fPtSpecFB32[%d][%d]",j,c), Form("fPtSpecFB32[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecFB32[j][c]);
        fPtSpecFB96[j][c] = new TH1F(Form("fPtSpecFB96[%d][%d]",j,c), Form("fPtSpecFB96[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecFB96[j][c]);
        fPtSpecFB128[j][c] = new TH1F(Form("fPtSpecFB128[%d][%d]",j,c), Form("fPtSpecFB128[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecFB128[j][c]);
        fPtSpecFB768[j][c] = new TH1F(Form("fPtSpecFB768[%d][%d]",j,c), Form("fPtSpecFB768[%d][%d]",j,c), 23, xmin);
        fSpectraMCList->Add(fPtSpecFB768[j][c]);
      }
    }
  }

  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);

  for(int i=0; i<5; i++){
    char hname[20];
    sprintf(hname,"hZNCPM%d",i);
    fhZNCPM[i] = new TH1F(hname, hname, 200, -50., 140000);
    fOutput->Add(fhZNCPM[i]);
    //
    sprintf(hname,"hZNAPM%d",i);
    fhZNAPM[i] = new TH1F(hname, hname, 200, -50., 140000);
    fOutput->Add(fhZNAPM[i]);
    //
    //@Shi add fhZPCPM and fhZPAPM
    sprintf(hname,"hZPCPM%d",i);
    fhZPCPM[i] = new TH1F(hname, hname, 200, -50., 140000);
    fOutput->Add(fhZPCPM[i]);
    //
    sprintf(hname,"hZPAPM%d",i);
    fhZPAPM[i] = new TH1F(hname, hname, 200, -50., 140000);
    fOutput->Add(fhZPAPM[i]);
    //
    if(i<4){
      //
      char hnamenc[20];
      sprintf(hnamenc, "hZNCPMQ%dPMC",i+1);
      fhZNCPMQiPMC[i] = new TH1F(hnamenc, hnamenc, 100, 0., 1.);
      fOutput->Add(fhZNCPMQiPMC[i]);
      //
      char hnamena[20];
      sprintf(hnamena, "hZNAPMQ%dPMC",i+1);
      fhZNAPMQiPMC[i] = new TH1F(hnamena, hnamena, 100, 0., 1.);
      fOutput->Add(fhZNAPMQiPMC[i]);
      //
      //@Shi add fhZPCPMQiPMC fhZPAPMQiPMC
      char hnamepc[20];
      sprintf(hnamepc, "hZPCPMQ%dPMC",i+1);
      fhZPCPMQiPMC[i] = new TH1F(hnamepc, hnamepc, 100, 0., 1.);
      fOutput->Add(fhZPCPMQiPMC[i]);
      //
      char hnamepa[20];
      sprintf(hnamepa, "hZPAPMQ%dPMC",i+1);
      fhZPAPMQiPMC[i] = new TH1F(hnamepa, hnamepa, 100, 0., 1.);
      fOutput->Add(fhZPAPMQiPMC[i]);
    }
  }

  fhZNCvsZNA = new TH2F("hZNCvsZNA","hZNCvsZNA",200,-50.,2.E5,200,-50.,2.E5);
  fOutput->Add(fhZNCvsZNA);
  fhZDCCvsZDCCA = new TH2F("hZDCCvsZDCCA","hZDCCvsZDCCA",200,0.,2.8E5,200,0.,2.8E5);
  fOutput->Add(fhZDCCvsZDCCA);
  fhZNCvsZPC = new TH2F("hZNCvsZPC","hZNCvsZPC",200,-50.,0.8E5,200,-50.,2.E5);
  fOutput->Add(fhZNCvsZPC);
  fhZNAvsZPA = new TH2F("hZNAvsZPA","hZNAvsZPA",200,-50.,0.8E5,200,-50.,2.E5);
  fOutput->Add(fhZNAvsZPA);
  fhZNvsZP = new TH2F("hZNvsZP","hZNvsZP",200,-50.,1.6E5,200,-50.,4.E5);
  fOutput->Add(fhZNvsZP);
  fhZNvsVZERO = new TH2F("hZNvsVZERO","hZNvsVZERO",250,0.,35000.,200,0.,4.E5);
  fOutput->Add(fhZNvsVZERO);
  fhZDCvsVZERO = new TH2F("hZDCvsVZERO","hZDCvsVZERO",250,0.,35000.,250,0.,5.6E5);
  fOutput->Add(fhZDCvsVZERO);

  fhAsymm = new TH1F("hAsymm" , "Asimmetry ",200,-1.,1.);
  fOutput->Add(fhAsymm);
  fhZNAvsAsymm = new TH2F("hZNAvsAsymm","ZNA vs. asymm.",200,-1.,1.,200,-50.,2.E5);
  fOutput->Add(fhZNAvsAsymm);
  fhZNCvsAsymm = new TH2F("hZNCvsAsymm","ZNC vs. asymm.",200,-1.,1.,200,-50.,2.E5);
  fOutput->Add(fhZNCvsAsymm);

  fhZNCvscentrality = new TH2F("hZNCvscentrality","hZNCvscentrality",100,0.,100.,200,-50.,2.E5);
  fOutput->Add(fhZNCvscentrality);
  fhZNAvscentrality = new TH2F("hZNAvscentrality","hZNAvscentrality",100,0.,100.,200,-50.,2.E5);
  fOutput->Add(fhZNAvscentrality);
  fhZPCvscentrality = new TH2F("hZPCvscentrality","hZPCvscentrality",100,0.,100.,200,-50.,0.8E5);
  fOutput->Add(fhZPCvscentrality);
  fhZPAvscentrality = new TH2F("hZPAvscentrality","hZPAvscentrality",100,0.,100.,200,-50.,0.8E5);
  fOutput->Add(fhZPAvscentrality);

  //@Shi add ZN and ZP corelation hists (begin)
  fZPAvsZNASignal = new TH3D("fZPAvsZNASignal","fZPAvsZNASignal",251,-1.,250.,251,-1.,250.,100,0,100);
  fZPAvsZNASignal->GetXaxis()->SetTitle("ZNA signal (a.u.) #times 10^{3}");
  fZPAvsZNASignal->GetYaxis()->SetTitle("ZPA signal (a.u.) #times 10^{3}");
  fZPAvsZNASignal->GetZaxis()->SetTitle("Centrality (%)");
  fZPAvsZNASignal->Sumw2();
  fOutput->Add(fZPAvsZNASignal);
  
  fZPCvsZNCSignal = new TH3D("fZPCvsZNCSignal","fZPCvsZNCSignal",251,-1.,250.,251,-1.,250.,100,0,100);
  fZPCvsZNCSignal->GetXaxis()->SetTitle("ZNC signal (a.u.) #times 10^{3}");
  fZPCvsZNCSignal->GetYaxis()->SetTitle("ZPC signal (a.u.) #times 10^{3}");
  fZPCvsZNCSignal->GetZaxis()->SetTitle("Centrality (%)");
  fZPCvsZNCSignal->Sumw2();
  fOutput->Add(fZPCvsZNCSignal);
  //@Shi add ZN and ZP corelation hists (end)
  
  //@Shi Add test histogram for gain equalization (begin)
  fZNenergyBeforeCalibration = new TProfile("fZNenergyBeforeCalibration","fZNenergyBeforeCalibration",10,0,10,"s");
  fZNenergyBeforeCalibration->SetStats(kFALSE);
  fOutput->Add(fZNenergyBeforeCalibration);
  
  fZNenergyAfterCalibration = new TProfile("fZNenergyAfterCalibration","fZNenergyAfterCalibration",10,0,10,"s");
  fZNenergyAfterCalibration->SetStats(kFALSE);
  fOutput->Add(fZNenergyAfterCalibration);
  //@Shi Add test histogram for gain equalization (end)
  //********************************************************************

  Int_t dRun10h[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};

  Int_t dRun11h[] = {167902, 167903, 167915, 167920, 167985, 167987, 167988, 168066, 168068, 168069, 168076, 168104, 168105, 168107, 168108, 168115, 168212, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168461, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168984, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169143, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169956, 169965, 170027, 170036,170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593};

  Int_t dRun15o[] = {244917, 244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683, 246148}; // @Shi add 246148

  Int_t dRun15ov6[] = {244918, 244975, 244980, 244982, 244983, 245064, 245066, 245068, 246390, 246391, 246392, 246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};

  Int_t dRun15opidfix[] = {245145, 245146, 245151, 245152, 245231, 245232, 245259, 245343, 245345, 245346, 245347, 245349, 245353, 245396, 245397, 245401, 245407, 245409, 245441, 245446, 245450, 245454, 245496, 245497, 245501, 245504, 245505, 245507, 245535, 245540, 245542, 245543, 245544, 245545, 245554};

  //@Shi start
  Int_t dRun18r[] = {297317, 297311, 297310, 297278, 297222, 297221, 297218, 297196, 297195, 297193, 297133, 297132, 297129, 297128, 297124, 297123, 297119, 297118, 297117, 297085, 297035, 297031, 296966, 296941, 296938, 296935, 296934};
  
  Double_t dVtxPosX15o[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0760396, 0.0761476, 0.0754612, 0.0760119, 0.0771416, 0.0758959, 0.0758307, 0.0764312, 0.0712992, 0.0740006, 0.076035, 0.0795941, 0.0758193, 0.0753836, 0.0759469, 0.0753271, 0.0748559, 0.0755779, 0.0747179, 0.0743789, 0.074226, 0.0738555, 0.0741127, 0.0755049, 0.079727, 0.0754529, 0.0747599, 0.0744282, 0.0742795, 0.0750923, 0.0765961, 0.0762358, 0.0765928, 0.0752035, 0.0767834, 0.0759724, 0.0758235, 0.0690952, 0.0693622, 0.0695388, 0.0704506, 0.070026, 0.0703322, 0.0702859, 0.0695319, 0.0684041, 0.0683909, 0.0696078, 0.0699702, 0.0689661, 0.0677066, 0.0689856, 0.0714685, 0.0690362, 0.0703379, 0.0692874, 0.0702451, 0.0693919, 0.0693631, 0.0702106, 0.0703336, 0.0696804, 0.0668393, 0.0696303, 0.0684486, 0.0693902, 0.0682269, 0.0686902, 0.0688619, 0.069442, 0.0705462, 0.0695982, 0.069336, 0.0685833, 0.0677059, 0.0690834, 0.0691257, 0.0690399, 0.0695431};
  
  Double_t dVtxPosY15o[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.336358, 0.336082, 0.33601, 0.336299, 0.338571, 0.336182, 0.336955, 0.33598, 0.337921, 0.334241, 0.335462, 0.337443, 0.334584, 0.336416, 0.335418, 0.336588, 0.338577, 0.335504, 0.336177, 0.336241, 0.338079, 0.338119, 0.337332, 0.336716, 0.340298, 0.337025, 0.337512, 0.337696, 0.336138, 0.338704, 0.336543, 0.337053, 0.335586, 0.335519, 0.335771, 0.334203, 0.335871, 0.32961, 0.329341, 0.328825, 0.330096, 0.328709, 0.329233, 0.329063, 0.329943, 0.330227, 0.329343, 0.330058, 0.32979, 0.330226, 0.330673, 0.330379, 0.325801, 0.329745, 0.327493, 0.329334, 0.329097, 0.331733, 0.330179, 0.329786, 0.330113, 0.327863, 0.331576, 0.329589, 0.329758, 0.32966, 0.329914, 0.329771, 0.330217, 0.327307, 0.32939, 0.329085, 0.329112, 0.331714, 0.327878, 0.331697, 0.330765, 0.331914, 0.33046};
  
  Double_t dVtxPosZ15o[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.575355, 0.371561, 0.513012, 0.471856, 0.406659, 0.284534, 0.257454, 0.199246, 0.325792, 0.212815, 0.26548, 0.348875, 0.368373, 0.431977, 0.528148, 0.508048, 0.478282, 0.436153, 0.330369, 0.406381, 0.40795, 0.411249, 0.445349, 0.412348, 0.391552, 0.353029, 0.338251, 0.251904, 0.293615, 0.544099, 0.352431, 0.221797, 0.232368, 0.35809, 0.234556, 0.300599, 0.375358, 0.418464, 0.476625, 0.385246, 0.333402, 0.314478, 0.326505, 0.375008, 0.289914, 0.410377, 0.33794, 0.331634, 0.347134, 0.343325, 0.367387, 0.400036, 0.307101, 0.300977, 0.357842, 0.377861, 0.401782, 0.432738, 0.446801, 0.43286, 0.416691, 0.423076, 0.398294, 0.479613, 0.422342, 0.443408, 0.455862, 0.656827, 0.704932, 0.289011, 0.392294, 0.419466, 0.396562, 0.377537, 0.347602, 0.296413, 0.405798, 0.462996, 0.440022};
  //@Shi end
  
  if(fDataSet==k2010) {fCRCnRun=92;}
  if(fDataSet==k2011) {fCRCnRun=119;}
  if(fDataSet==k2015) {
	  fCRCnRun=91;
	  fAvVtxPosX=TArrayD(fCRCnRun,dVtxPosX15o);
      fAvVtxPosY=TArrayD(fCRCnRun,dVtxPosY15o);
      fAvVtxPosZ=TArrayD(fCRCnRun,dVtxPosZ15o);
  } // @Shi add 246148 which makes it become 91
  if(fDataSet==k2015v6) {fCRCnRun=91;}
  if(fDataSet==k2015pidfix) {fCRCnRun=35;}
  if(fDataSet==kAny) {fCRCnRun=1;}
  //@Shi start
  if(fDataSet==k2018r) {fCRCnRun=27;}
  //@Shi end
  Int_t d=0;
  for(Int_t r=0; r<fCRCnRun; r++) {
    if(fDataSet==k2010)   fRunList[d] = dRun10h[r];
    if(fDataSet==k2011)   fRunList[d] = dRun11h[r];
    if(fDataSet==k2015)   fRunList[d] = dRun15o[r];
    if(fDataSet==k2015v6) fRunList[d] = dRun15ov6[r];
    if(fDataSet==k2015pidfix) fRunList[d] = dRun15opidfix[r];
    if(fDataSet==kAny) fRunList[d] = 1;
    //@Shi start
    if(fDataSet==k2018r) fRunList[d] = dRun18r[r];
    //@Shi end
    d++;
  }

  fVZEROMult = new TProfile2D("fVZEROMult","fVZEROMult",fCRCnRun,0.,1.*fCRCnRun,64,0.,64.);
  for (Int_t i=0; i<fCRCnRun; i++) {
    fVZEROMult->GetXaxis()->SetBinLabel(i+1,Form("%d",fRunList[i]));
  }
  fVZEROStuffList->Add(fVZEROMult);

  if(fVZEROGainEqList) {
    fVZEROGainEqHist = (TH2D*)fVZEROGainEqList->FindObject("VZEROEqGain");
    fVZEROStuffList->Add(fVZEROGainEqHist);
  }
  if(fVZEROQVecRecList) {
    for (Int_t k=0; k<fkVZEROnHar; k++) {
      fVZEROQVectorRecQxStored[k] = (TProfile3D*)fVZEROQVecRecList->FindObject(Form("fVZEROQVectorRecQx[%d]",k));
      fVZEROQVectorRecQxStored[k]->SetTitle(Form("fVZEROQVectorRecQxStored[%d]",k));
      fVZEROQVectorRecQxStored[k]->SetName(Form("fVZEROQVectorRecQxStored[%d]",k));
      fVZEROStuffList->Add(fVZEROQVectorRecQxStored[k]);
      fVZEROQVectorRecQyStored[k] = (TProfile3D*)fVZEROQVecRecList->FindObject(Form("fVZEROQVectorRecQy[%d]",k));
      fVZEROQVectorRecQyStored[k]->SetTitle(Form("fVZEROQVectorRecQyStored[%d]",k));
      fVZEROQVectorRecQyStored[k]->SetName(Form("fVZEROQVectorRecQyStored[%d]",k));
      fVZEROStuffList->Add(fVZEROQVectorRecQyStored[k]);
      for (Int_t t=0; t<fkVZEROnQAplots; t++) {
        fVZEROQVectorRecFinal[k][t] = new TProfile2D(Form("fVZEROQVectorRecFinal[%d][%d]",k,t),Form("fVZEROQVectorRecFinal[%d][%d]",k,t),fCRCnRun,0.,1.*fCRCnRun,100,0.,100.,"s");
        fVZEROQVectorRecFinal[k][t]->Sumw2();
        fVZEROStuffList->Add(fVZEROQVectorRecFinal[k][t]);
      }
    }
  }

  //  for (Int_t k=0; k<fkVZEROnHar; k++) {
  //    fVZEROQVectorRecQx[k] = new TProfile3D(Form("fVZEROQVectorRecQx[%d]",k),Form("fVZEROQVectorRecQx[%d]",k),fCRCnRun,0.,1.*fCRCnRun,100,0.,100.,8,0.,8.,"s");
  //    fVZEROQVectorRecQx[k]->Sumw2();
  //    fVZEROStuffList->Add(fVZEROQVectorRecQx[k]);
  //    fVZEROQVectorRecQy[k] = new TProfile3D(Form("fVZEROQVectorRecQy[%d]",k),Form("fVZEROQVectorRecQy[%d]",k),fCRCnRun,0.,1.*fCRCnRun,100,0.,100.,8,0.,8.,"s");
  //    fVZEROQVectorRecQy[k]->Sumw2();
  //    fVZEROStuffList->Add(fVZEROQVectorRecQy[k]);
  //  }

  // track QA

  if(fAnalysisType == kTrackQA) {

    fTrackQAList = new TList();
    fTrackQAList->SetOwner(kTRUE);
    fTrackQAList->SetName("TrackQA");
    fOutput->Add(fTrackQAList);

    Int_t nBinsEta=10;
    Double_t dBinsEta[]={-0.8,-0.64,-0.48,-0.32,-0.16,0.,0.16,0.32,0.48,0.64,0.8};
    Int_t nBinsPt = 26;
    Double_t dBinsPt[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.33,2.66,3.,3.5,4.,5.,6.,8.,10.,14.,20.,30.,50.};

    for(Int_t fb=0; fb<fKNFBs; fb++){
      for (Int_t i=0; i<4; i++) {
        // DCA
        fTrackQADCAxy[fb][i] = new TH3D(Form("fTrackQADCAxy[%d][%d]",fb,i),";#eta;p_{T} [GeV/c];DCAxy [cm]",10,-0.8,0.8,10,0.,5.,480,-2.4,2.4);
        fTrackQAList->Add(fTrackQADCAxy[fb][i]);
        fTrackQADCAz[fb][i] = new TH3D(Form("fTrackQADCAz[%d][%d]",fb,i),";#eta;p_{T} [GeV/c];DCAz [cm]",10,-0.8,0.8,10,0.,5.,640,-3.2,3.2);
        fTrackQAList->Add(fTrackQADCAz[fb][i]);
        // pT
        fTrackQApT[fb][i] = new TH2D(Form("fTrackQApT[%d][%d]",fb,i),";#eta;p_{T} [GeV/c]",nBinsEta,dBinsEta,nBinsPt,dBinsPt);
        fTrackQAList->Add(fTrackQApT[fb][i]);
        // Dphi
        fTrackQADphi[fb][i] = new TProfile2D(Form("fTrackQADphi[%d][%d]",fb,i),";#eta;p_{T} [GeV/c];cos(#Delta#phi)",nBinsEta,dBinsEta,nBinsPt,dBinsPt);
        fTrackQAList->Add(fTrackQADphi[fb][i]);
        // EbE
        fEbEQRe[fb][i] = new TH2D(Form("fEbEQRe[%d][%d]",fb,i),";#eta;p_{T} [GeV/c];QRe(EbE)",nBinsEta,dBinsEta,nBinsPt,dBinsPt);
        fTrackQAList->Add(fEbEQRe[fb][i]);
        fEbEQIm[fb][i] = new TH2D(Form("fEbEQIm[%d][%d]",fb,i),";#eta;p_{T} [GeV/c];QRe(EbE)",nBinsEta,dBinsEta,nBinsPt,dBinsPt);
        fTrackQAList->Add(fEbEQIm[fb][i]);
        fEbEQMu[fb][i] = new TH2D(Form("fEbEQMu[%d][%d]",fb,i),";#eta;p_{T} [GeV/c];QRe(EbE)",nBinsEta,dBinsEta,nBinsPt,dBinsPt);
        fTrackQAList->Add(fEbEQMu[fb][i]);
      }
    }

  }

  // run-by-run stuff

  if(!fUseTowerEq) {
    Double_t ptmin[] = {0.2,0.4,0.6,0.8,1.,1.2,1.4,1.8,2.2,3.,4.,6.,8.,12.,20.};
    Double_t phimin[] = {0.,TMath::Pi()/8.,2*TMath::Pi()/8.,3*TMath::Pi()/8.,4*TMath::Pi()/8.,5*TMath::Pi()/8.,6*TMath::Pi()/8.,7*TMath::Pi()/8.,8*TMath::Pi()/8.,9*TMath::Pi()/8.,10*TMath::Pi()/8.,11*TMath::Pi()/8.,12*TMath::Pi()/8.,13*TMath::Pi()/8.,14*TMath::Pi()/8.,15*TMath::Pi()/8.,16*TMath::Pi()/8.};
    Double_t etamin[] = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};

    for(Int_t r=0;r<fCRCnRun;r++) {
      fCRCQVecListRun[r] = new TList();
      fCRCQVecListRun[r]->SetName(Form("Run %d",fRunList[r]));
      fCRCQVecListRun[r]->SetOwner(kTRUE);
      fOutput->Add(fCRCQVecListRun[r]);
      
      //@Shi add List Run 
      if (fStepZDCRecenter>=0) {
        if (r<30) {
          fRecenter1ListRunbyRun[r] = new TList();
          fRecenter1ListRunbyRun[r]->SetName(Form("Run %d",fRunList[r]));
          fRecenter1ListRunbyRun[r]->SetOwner(kTRUE);
          fOutputRecenter1->Add(fRecenter1ListRunbyRun[r]);
	    } else if (r<60) {
		  fRecenter2ListRunbyRun[r-30] = new TList(); // r-30 so that the index starts at 0
          fRecenter2ListRunbyRun[r-30]->SetName(Form("Run %d",fRunList[r]));
          fRecenter2ListRunbyRun[r-30]->SetOwner(kTRUE);
          fOutputRecenter2->Add(fRecenter2ListRunbyRun[r-30]); 
	    } else {
		  fRecenter3ListRunbyRun[r-60] = new TList(); // r-60 so that the index starts at 0
          fRecenter3ListRunbyRun[r-60]->SetName(Form("Run %d",fRunList[r]));
          fRecenter3ListRunbyRun[r-60]->SetOwner(kTRUE);
          fOutputRecenter3->Add(fRecenter3ListRunbyRun[r-60]); 
		}
	  }

	  //@Shi Add run by run ZN centroid vs centrality (begin)
      for (Int_t c=0; c<2; c++) {
        fhZNCenDisRbR[r][c] = new TH3D(Form("fhZNCenDisRbR[%d][%d]",fRunList[r],c), Form("fhZNCenDisRbR[%d][%d]",fRunList[r],c), 100, 0., 100., 100, -2., 2. , 100., -2., 2.);
        fCRCQVecListRun[r]->Add(fhZNCenDisRbR[r][c]);
	  }
	  //@Shi Add run by run ZN centroid vs centrality (end)
	  
      for(Int_t k=0;k<fCRCnTow;k++) {
        fZNCTower[r][k] = new TProfile(Form("fZNCTower[%d][%d]",fRunList[r],k),Form("fZNCTower[%d][%d]",fRunList[r],k),100,0.,100.,"s");
        fZNCTower[r][k]->Sumw2();
        fCRCQVecListRun[r]->Add(fZNCTower[r][k]);
        fZNATower[r][k] = new TProfile(Form("fZNATower[%d][%d]",fRunList[r],k),Form("fZNATower[%d][%d]",fRunList[r],k),100,0.,100.,"s");
        fZNATower[r][k]->Sumw2();
        fCRCQVecListRun[r]->Add(fZNATower[r][k]);
        //@Shi add ZPCTower and ZPATower (begin)
        fZPCTower[r][k] = new TProfile(Form("fZPCTower[%d][%d]",fRunList[r],k),Form("fZPCTower[%d][%d]",fRunList[r],k),100,0.,100.,"s");
        fZPCTower[r][k]->Sumw2();
        fCRCQVecListRun[r]->Add(fZPCTower[r][k]);
        fZPATower[r][k] = new TProfile(Form("fZPATower[%d][%d]",fRunList[r],k),Form("fZPATower[%d][%d]",fRunList[r],k),100,0.,100.,"s");
        fZPATower[r][k]->Sumw2();
        fCRCQVecListRun[r]->Add(fZPATower[r][k]);
        //@Shi add ZPCTower and ZPATower (end)
      }

      //    fhZNSpectraRbR[r] = new TH3D(Form("fhZNSpectraRbR[%d]",fRunList[r]),Form("fhZNSpectraRbR[%d]",fRunList[r]),50,0.,100.,8,0.,8.,100,0.,1.E5);
      //    fCRCQVecListRun[r]->Add(fhZNSpectraRbR[r]);

      //   for(Int_t i=0;i<fnCen;i++) {
      //     fPtPhiEtaRbRFB128[r][i] = new TH3F(Form("fPtPhiEtaRbRFB128[%d][%d]",r,i),Form("fPtPhiEtaRbRFB128[%d][%d]",r,i),14, ptmin, 16, phimin, 16, etamin);
      //     fCRCQVecListRun[r]->Add(fPtPhiEtaRbRFB128[r][i]);
      //     fPtPhiEtaRbRFB768[r][i] = new TH3F(Form("fPtPhiEtaRbRFB768[%d][%d]",r,i),Form("fPtPhiEtaRbRFB768[%d][%d]",r,i),14, ptmin, 16, phimin, 16, etamin);
      //     fCRCQVecListRun[r]->Add(fPtPhiEtaRbRFB768[r][i]);
      //   }
    }
    
    //@Shi Add run by run recentering histograms for ZDC (begin) (if !fUseTowerEq = kTRUE)
    
    if (fStepZDCRecenter>=0) {
		fAve_VtxX = new TProfile("fAve_VtxX", "fAve_VtxX", fCRCnRun, 0, fCRCnRun, "s");
		fAve_VtxY = new TProfile("fAve_VtxY", "fAve_VtxY", fCRCnRun, 0, fCRCnRun, "s");
		fAve_VtxZ = new TProfile("fAve_VtxZ", "fAve_VtxZ", fCRCnRun, 0, fCRCnRun, "s");
		fOutputRecenter1->Add(fAve_VtxX);
		fOutputRecenter1->Add(fAve_VtxY);
		fOutputRecenter1->Add(fAve_VtxZ);
		if (fStoreCalibZDCRecenter){
		  for (Int_t r=0;r<fCRCnRun;r++) {
			for (Int_t c=0; c<4; c++) {
			  if (fStepZDCRecenter>=0) {
				// vertex_x: range: [0.08, 0.1], bins: 40
				if (fDataSet==k2015) {
				  fRun_VtxXQPreCalib[r][c] = new TProfile(Form("fRun_VtxXQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxXQPreCalib[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxXQPreCalib[r][c] = new TProfile(Form("fRun_VtxXQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxXQPreCalib[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, "s");
				}
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxXQPreCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXQPreCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXQPreCalib[r][c]);
				}
				// vertex y: range: [0.36, 0.38], bins: 40
				if (fDataSet==k2015) {
				  fRun_VtxYQPreCalib[r][c] = new TProfile(Form("fRun_VtxYQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxYQPreCalib[%d][%d]",fRunList[r],c), 50, -0.010, 0.010, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxYQPreCalib[r][c] = new TProfile(Form("fRun_VtxYQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxYQPreCalib[%d][%d]",fRunList[r],c), 40, 0.36, 0.38, "s");
				}
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxYQPreCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxYQPreCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxYQPreCalib[r][c]);
				}
				// vertex z: range: [-10, 10], bins: 40
				fRun_VtxZQPreCalib[r][c] = new TProfile(Form("fRun_VtxZQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxZQPreCalib[%d][%d]",fRunList[r],c), 40, -10, 10, "s");
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxZQPreCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxZQPreCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxZQPreCalib[r][c]);
				}
			  } 
			  if (fStepZDCRecenter>=1) {
				// centrality: 1% range: [0., 100.], bins: 100
				fRun_CentQCalib[r][c] = new TProfile(Form("fRun_CentQCalib[%d][%d]",fRunList[r],c), Form("fRun_CentQCalib[%d][%d]",fRunList[r],c), 100, 0., 100.,"s");
				if (fDataSet==k2015) {
				  fRun_VtxXQCalibStep1[r][c] = new TProfile(Form("fRun_VtxXQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalibStep1[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, "s");
				  fRun_VtxYQCalibStep1[r][c] = new TProfile(Form("fRun_VtxYQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalibStep1[%d][%d]",fRunList[r],c), 50, -0.010, 0.010, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxXQCalibStep1[r][c] = new TProfile(Form("fRun_VtxXQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalibStep1[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, "s");
				  fRun_VtxYQCalibStep1[r][c] = new TProfile(Form("fRun_VtxYQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalibStep1[%d][%d]",fRunList[r],c), 40, 0.36, 0.38, "s");
				}
				fRun_VtxZQCalibStep1[r][c] = new TProfile(Form("fRun_VtxZQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxZQCalibStep1[%d][%d]",fRunList[r],c), 40, -10, 10, "s");

				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_CentQCalib[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxXQCalibStep1[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxYQCalibStep1[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxZQCalibStep1[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_CentQCalib[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXQCalibStep1[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxYQCalibStep1[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxZQCalibStep1[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_CentQCalib[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXQCalibStep1[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxYQCalibStep1[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxZQCalibStep1[r][c]);
				}
			  } 
				
			  if (fStepZDCRecenter>=2) {
				if (fDataSet==k2015) {
				  fRun_VtxXQCalibStep2[r][c] = new TProfile(Form("fRun_VtxXQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalibStep2[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, "s");
				  fRun_VtxYQCalibStep2[r][c] = new TProfile(Form("fRun_VtxYQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalibStep2[%d][%d]",fRunList[r],c), 50, -0.010, 0.010, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxXQCalibStep2[r][c] = new TProfile(Form("fRun_VtxXQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalibStep2[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, "s");
				  fRun_VtxYQCalibStep2[r][c] = new TProfile(Form("fRun_VtxYQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalibStep2[%d][%d]",fRunList[r],c), 40, 0.36, 0.38, "s");
				}
				fRun_VtxZQCalibStep2[r][c] = new TProfile(Form("fRun_VtxZQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxZQCalibStep2[%d][%d]",fRunList[r],c), 40, -10, 10, "s");

				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxXQCalibStep2[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxYQCalibStep2[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxZQCalibStep2[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXQCalibStep2[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxYQCalibStep2[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxZQCalibStep2[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXQCalibStep2[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxYQCalibStep2[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxZQCalibStep2[r][c]);
				}
			  } 
				
			  if (fStepZDCRecenter>=3) {
				// vertex_x: range: [0.08, 0.1], bins: 40
				if (fDataSet==k2015) {
				  fRun_VtxXQCalib[r][c] = new TProfile(Form("fRun_VtxXQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalib[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxXQCalib[r][c] = new TProfile(Form("fRun_VtxXQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalib[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, "s");
				}
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxXQCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXQCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXQCalib[r][c]);
				}
				// vertex y: range: [0.36, 0.38], bins: 40
				if (fDataSet==k2015) {
				  fRun_VtxYQCalib[r][c] = new TProfile(Form("fRun_VtxYQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalib[%d][%d]",fRunList[r],c), 50, -0.010, 0.010, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxYQCalib[r][c] = new TProfile(Form("fRun_VtxYQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalib[%d][%d]",fRunList[r],c), 40, 0.36, 0.38, "s");
				}
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxYQCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxYQCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxYQCalib[r][c]);
				}
				// vertex z: range: [-10, 10], bins: 40
				fRun_VtxZQCalib[r][c] = new TProfile(Form("fRun_VtxZQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxZQCalib[%d][%d]",fRunList[r],c), 40, -10, 10, "s");
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxZQCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxZQCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxZQCalib[r][c]);
				}
				// centrality: 1% range: [0., 100.], bins: 100
				fRun_CentQCalib2[r][c] = new TProfile(Form("fRun_CentQCalib2[%d][%d]",fRunList[r],c), Form("fRun_CentQCalib2[%d][%d]",fRunList[r],c), 100, 0., 100.,"s");
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_CentQCalib2[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_CentQCalib2[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_CentQCalib2[r][c]);
				}
			  }
			}
		  }
		  fCorrQAReCRe = new TProfile("CorrQAReCRe", "Correlation QARe x QCRe", 100, 0., 100.);
		  fCorrQAReCIm = new TProfile("CorrQAReCIm", "Correlation QARe x QCIm", 100, 0., 100.);
		  fCorrQAImCRe = new TProfile("CorrQAImCRe", "Correlation QAIm x QCRe", 100, 0., 100.);
		  fCorrQAImCIm = new TProfile("CorrQAImCIm", "Correlation QAIm x QCIm", 100, 0., 100.);
		  fOutputRecenter1->Add(fCorrQAReCRe);
		  fOutputRecenter1->Add(fCorrQAReCIm);
		  fOutputRecenter1->Add(fCorrQAImCRe);
		  fOutputRecenter1->Add(fCorrQAImCIm);
		}

		for (Int_t r=0;r<fCRCnRun;r++) {
		  for (Int_t c=0; c<4; c++) {
			if (fStepZDCRecenter >= 0) {
			  // centrality: 1% range: [0., 100.], bins: 100
			  fRun_CentQ[r][c] = new TProfile(Form("fRun_CentQ[%d][%d]",fRunList[r],c), Form("fRun_CentQ[%d][%d]",fRunList[r],c), 100,0.,100.,"s");
			  if (r<30) {
				fRecenter1ListRunbyRun[r]->Add(fRun_CentQ[r][c]);
			  } else if (r<60) {
				fRecenter2ListRunbyRun[r-30]->Add(fRun_CentQ[r][c]);
			  } else {
				fRecenter3ListRunbyRun[r-60]->Add(fRun_CentQ[r][c]);
			  }
			}
			  
			if (fStepZDCRecenter >= 2) {
			  if (fDataSet==k2015) {
				fRun_VtxXYZQ[r][c] = new TProfile3D(Form("fRun_VtxXYZQ[%d][%d]",fRunList[r],c), Form("fRun_VtxXYZQ[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, 50, -0.010, 0.010, 40, -10, 10, "s");
			  } else if (fDataSet==k2018r) {
				fRun_VtxXYZQ[r][c] = new TProfile3D(Form("fRun_VtxXYZQ[%d][%d]",fRunList[r],c), Form("fRun_VtxXYZQ[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, 40, 0.36, 0.38, 40, -10, 10, "s");
			  }
			  if (r<30) {
				fRecenter1ListRunbyRun[r]->Add(fRun_VtxXYZQ[r][c]);
			  } else if (r<60) {
				fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXYZQ[r][c]);
			  } else {
				fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXYZQ[r][c]);
			  }
			}
			  
			/*if (fStepZDCRecenter >= 3) {
			  // not step 3 hist to be saved
			}*/
		  }
		}
		
		//const Int_t fnCentBinForRecentering = 20; // this means that a wider centrality bin is used {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100}
		for (Int_t r=0;r<fnCentBinForRecentering;r++) {
		  for (Int_t c=0; c<4; c++) {
			if (fStepZDCRecenter >= 1) {
			  // vertex_x, y, z
			  if (fDataSet==k2015) {
				fCent_VtxXYZQ[r][c] = new TProfile3D(Form("fCent_VtxXYZQ[%d][%d]",r,c), Form("fCent_VtxXYZQ[%d][%d]",r,c), 50, -0.01, 0.01, 50, -0.010, 0.010, 40, -10, 10, "s");
			  } else if (fDataSet==k2018r) {
				fCent_VtxXYZQ[r][c] = new TProfile3D(Form("fCent_VtxXYZQ[%d][%d]",r,c), Form("fCent_VtxXYZQ[%d][%d]",r,c), 40, 0.08, 0.1, 40, 0.36, 0.38, 40, -10, 10, "s");
			  }
			  fOutputRecenter1->Add(fCent_VtxXYZQ[r][c]);
			}
		  }
		}
	}
    //@Shi Add run by run recentering histograms for ZDC (end)
    
  } else { //@Shi add for !fUseTowerEq = False (begin)
	for(Int_t r=0;r<fCRCnRun;r++) {
	    fCRCQVecListRun[r] = new TList();
        fCRCQVecListRun[r]->SetName(Form("Run %d",fRunList[r]));
        fCRCQVecListRun[r]->SetOwner(kTRUE);
        fOutput->Add(fCRCQVecListRun[r]);
	
		//@Shi add List Run 
		if (fStepZDCRecenter>=0) {
          if (r<30) {
            fRecenter1ListRunbyRun[r] = new TList();
            fRecenter1ListRunbyRun[r]->SetName(Form("Run %d",fRunList[r]));
            fRecenter1ListRunbyRun[r]->SetOwner(kTRUE);
            fOutputRecenter1->Add(fRecenter1ListRunbyRun[r]);
	      } else if (r<60) {
		    fRecenter2ListRunbyRun[r-30] = new TList(); // r-30 so that the index starts at 0
            fRecenter2ListRunbyRun[r-30]->SetName(Form("Run %d",fRunList[r]));
            fRecenter2ListRunbyRun[r-30]->SetOwner(kTRUE);
            fOutputRecenter2->Add(fRecenter2ListRunbyRun[r-30]); 
	      } else {
		    fRecenter3ListRunbyRun[r-60] = new TList(); // r-60 so that the index starts at 0
            fRecenter3ListRunbyRun[r-60]->SetName(Form("Run %d",fRunList[r]));
            fRecenter3ListRunbyRun[r-60]->SetOwner(kTRUE);
            fOutputRecenter3->Add(fRecenter3ListRunbyRun[r-60]); 
		  }
		}
		  
		if (fFillZNCenDisRbR) {
		  //@Shi Add run by run ZN centroid vs centrality
		  for(Int_t r=0;r<fCRCnRun;r++) {
			for (Int_t c=0; c<2; c++) {
			  fhZNCenDisRbR[r][c] = new TH3D(Form("fhZNCenDisRbR[%d][%d]",fRunList[r],c), Form("fhZNCenDisRbR[%d][%d]",fRunList[r],c), 100, 0., 100., 100, -2., 2. , 100., -2., 2.);
			  fCRCQVecListRun[r]->Add(fhZNCenDisRbR[r][c]);
			}
		  }
		}
    }
    
    //@Shi Add run by run recentering histograms for ZDC (begin) (if !fUseTowerEq = kTRUE)
    if (fStepZDCRecenter>=0) {
		fAve_VtxX = new TProfile("fAve_VtxX", "fAve_VtxX", fCRCnRun, 0, fCRCnRun, "s");
		fAve_VtxY = new TProfile("fAve_VtxY", "fAve_VtxY", fCRCnRun, 0, fCRCnRun, "s");
		fAve_VtxZ = new TProfile("fAve_VtxZ", "fAve_VtxZ", fCRCnRun, 0, fCRCnRun, "s");
		fOutputRecenter1->Add(fAve_VtxX);
		fOutputRecenter1->Add(fAve_VtxY);
		fOutputRecenter1->Add(fAve_VtxZ);
		
		if (fStoreCalibZDCRecenter){
		  for (Int_t r=0;r<fCRCnRun;r++) {
			for (Int_t c=0; c<4; c++) {
			  if (fStepZDCRecenter>=0) {
				// vertex_x: range: [0.08, 0.1], bins: 40
				if (fDataSet==k2015) {
				  fRun_VtxXQPreCalib[r][c] = new TProfile(Form("fRun_VtxXQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxXQPreCalib[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxXQPreCalib[r][c] = new TProfile(Form("fRun_VtxXQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxXQPreCalib[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, "s");
				}
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxXQPreCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXQPreCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXQPreCalib[r][c]);
				}
				// vertex y: range: [0.36, 0.38], bins: 40
				if (fDataSet==k2015) {
				  fRun_VtxYQPreCalib[r][c] = new TProfile(Form("fRun_VtxYQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxYQPreCalib[%d][%d]",fRunList[r],c), 50, -0.010, 0.010, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxYQPreCalib[r][c] = new TProfile(Form("fRun_VtxYQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxYQPreCalib[%d][%d]",fRunList[r],c), 40, 0.36, 0.38, "s");
				}
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxYQPreCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxYQPreCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxYQPreCalib[r][c]);
				}
				// vertex z: range: [-10, 10], bins: 40
				fRun_VtxZQPreCalib[r][c] = new TProfile(Form("fRun_VtxZQPreCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxZQPreCalib[%d][%d]",fRunList[r],c), 40, -10, 10, "s");
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxZQPreCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxZQPreCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxZQPreCalib[r][c]);
				}
			  } 
			  if (fStepZDCRecenter>=1) {
				// centrality: 1% range: [0., 100.], bins: 100
				fRun_CentQCalib[r][c] = new TProfile(Form("fRun_CentQCalib[%d][%d]",fRunList[r],c), Form("fRun_CentQCalib[%d][%d]",fRunList[r],c), 100, 0., 100.,"s");
				if (fDataSet==k2015) {
				  fRun_VtxXQCalibStep1[r][c] = new TProfile(Form("fRun_VtxXQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalibStep1[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, "s");
				  fRun_VtxYQCalibStep1[r][c] = new TProfile(Form("fRun_VtxYQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalibStep1[%d][%d]",fRunList[r],c), 50, -0.010, 0.010, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxXQCalibStep1[r][c] = new TProfile(Form("fRun_VtxXQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalibStep1[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, "s");
				  fRun_VtxYQCalibStep1[r][c] = new TProfile(Form("fRun_VtxYQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalibStep1[%d][%d]",fRunList[r],c), 40, 0.36, 0.38, "s");
				}
				fRun_VtxZQCalibStep1[r][c] = new TProfile(Form("fRun_VtxZQCalibStep1[%d][%d]",fRunList[r],c), Form("fRun_VtxZQCalibStep1[%d][%d]",fRunList[r],c), 40, -10, 10, "s");

				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_CentQCalib[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxXQCalibStep1[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxYQCalibStep1[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxZQCalibStep1[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_CentQCalib[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXQCalibStep1[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxYQCalibStep1[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxZQCalibStep1[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_CentQCalib[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXQCalibStep1[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxYQCalibStep1[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxZQCalibStep1[r][c]);
				}
			  } 
				
			  if (fStepZDCRecenter>=2) {
				if (fDataSet==k2015) {
				  fRun_VtxXQCalibStep2[r][c] = new TProfile(Form("fRun_VtxXQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalibStep2[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, "s");
				  fRun_VtxYQCalibStep2[r][c] = new TProfile(Form("fRun_VtxYQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalibStep2[%d][%d]",fRunList[r],c), 50, -0.010, 0.010, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxXQCalibStep2[r][c] = new TProfile(Form("fRun_VtxXQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalibStep2[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, "s");
				  fRun_VtxYQCalibStep2[r][c] = new TProfile(Form("fRun_VtxYQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalibStep2[%d][%d]",fRunList[r],c), 40, 0.36, 0.38, "s");
				}
				fRun_VtxZQCalibStep2[r][c] = new TProfile(Form("fRun_VtxZQCalibStep2[%d][%d]",fRunList[r],c), Form("fRun_VtxZQCalibStep2[%d][%d]",fRunList[r],c), 40, -10, 10, "s");

				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxXQCalibStep2[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxYQCalibStep2[r][c]);
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxZQCalibStep2[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXQCalibStep2[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxYQCalibStep2[r][c]);
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxZQCalibStep2[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXQCalibStep2[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxYQCalibStep2[r][c]);
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxZQCalibStep2[r][c]);
				}
			  } 
				
			  if (fStepZDCRecenter>=3) {
				// vertex_x: range: [0.08, 0.1], bins: 40
				if (fDataSet==k2015) {
				  fRun_VtxXQCalib[r][c] = new TProfile(Form("fRun_VtxXQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalib[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxXQCalib[r][c] = new TProfile(Form("fRun_VtxXQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxXQCalib[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, "s");
				}
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxXQCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXQCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXQCalib[r][c]);
				}
				// vertex y: range: [0.36, 0.38], bins: 40
				if (fDataSet==k2015) {
				  fRun_VtxYQCalib[r][c] = new TProfile(Form("fRun_VtxYQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalib[%d][%d]",fRunList[r],c), 50, -0.010, 0.010, "s");
				} else if (fDataSet==k2018r) {
				  fRun_VtxYQCalib[r][c] = new TProfile(Form("fRun_VtxYQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxYQCalib[%d][%d]",fRunList[r],c), 40, 0.36, 0.38, "s");
				}
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxYQCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxYQCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxYQCalib[r][c]);
				}
				// vertex z: range: [-10, 10], bins: 40
				fRun_VtxZQCalib[r][c] = new TProfile(Form("fRun_VtxZQCalib[%d][%d]",fRunList[r],c), Form("fRun_VtxZQCalib[%d][%d]",fRunList[r],c), 40, -10, 10, "s");
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_VtxZQCalib[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxZQCalib[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxZQCalib[r][c]);
				}
				// centrality: 1% range: [0., 100.], bins: 100
				fRun_CentQCalib2[r][c] = new TProfile(Form("fRun_CentQCalib2[%d][%d]",fRunList[r],c), Form("fRun_CentQCalib2[%d][%d]",fRunList[r],c), 100, 0., 100.,"s");
				if (r<30) {
				  fRecenter1ListRunbyRun[r]->Add(fRun_CentQCalib2[r][c]);
				} else if (r<60) {
				  fRecenter2ListRunbyRun[r-30]->Add(fRun_CentQCalib2[r][c]);
				} else {
				  fRecenter3ListRunbyRun[r-60]->Add(fRun_CentQCalib2[r][c]);
				}
			  }
			}
		  }
		  fCorrQAReCRe = new TProfile("CorrQAReCRe", "Correlation QARe x QCRe", 100, 0., 100.);
		  fCorrQAReCIm = new TProfile("CorrQAReCIm", "Correlation QARe x QCIm", 100, 0., 100.);
		  fCorrQAImCRe = new TProfile("CorrQAImCRe", "Correlation QAIm x QCRe", 100, 0., 100.);
		  fCorrQAImCIm = new TProfile("CorrQAImCIm", "Correlation QAIm x QCIm", 100, 0., 100.);
		  fOutputRecenter1->Add(fCorrQAReCRe);
		  fOutputRecenter1->Add(fCorrQAReCIm);
		  fOutputRecenter1->Add(fCorrQAImCRe);
		  fOutputRecenter1->Add(fCorrQAImCIm);
		}

		for (Int_t r=0;r<fCRCnRun;r++) {
		  for (Int_t c=0; c<4; c++) {
			if (fStepZDCRecenter >= 0) {
			  // centrality: 1% range: [0., 100.], bins: 100
			  fRun_CentQ[r][c] = new TProfile(Form("fRun_CentQ[%d][%d]",fRunList[r],c), Form("fRun_CentQ[%d][%d]",fRunList[r],c), 100,0.,100.,"s");
			  if (r<30) {
				fRecenter1ListRunbyRun[r]->Add(fRun_CentQ[r][c]);
			  } else if (r<60) {
				fRecenter2ListRunbyRun[r-30]->Add(fRun_CentQ[r][c]);
			  } else {
				fRecenter3ListRunbyRun[r-60]->Add(fRun_CentQ[r][c]);
			  }
			}
			  
			if (fStepZDCRecenter >= 2) {
			  if (fDataSet==k2015) {
				fRun_VtxXYZQ[r][c] = new TProfile3D(Form("fRun_VtxXYZQ[%d][%d]",fRunList[r],c), Form("fRun_VtxXYZQ[%d][%d]",fRunList[r],c), 50, -0.01, 0.01, 50, -0.010, 0.010, 40, -10, 10, "s");
			  } else if (fDataSet==k2018r) {
				fRun_VtxXYZQ[r][c] = new TProfile3D(Form("fRun_VtxXYZQ[%d][%d]",fRunList[r],c), Form("fRun_VtxXYZQ[%d][%d]",fRunList[r],c), 40, 0.08, 0.1, 40, 0.36, 0.38, 40, -10, 10, "s");  
			  }
			  if (r<30) {
				fRecenter1ListRunbyRun[r]->Add(fRun_VtxXYZQ[r][c]);
			  } else if (r<60) {
				fRecenter2ListRunbyRun[r-30]->Add(fRun_VtxXYZQ[r][c]);
			  } else {
				fRecenter3ListRunbyRun[r-60]->Add(fRun_VtxXYZQ[r][c]);
			  }
			}
			  
			/*if (fStepZDCRecenter >= 3) {
			  // not step 3 hist to be saved
			}*/
		  }
		}
		
		//const Int_t fnCentBinForRecentering = 20; // this means that a wider centrality bin is used {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100}
		for (Int_t r=0;r<fnCentBinForRecentering;r++) {
		  for (Int_t c=0; c<4; c++) {
			if (fStepZDCRecenter >= 1) {
			  // vertex_x, y, z
			  if (fDataSet==k2015) {
				fCent_VtxXYZQ[r][c] = new TProfile3D(Form("fCent_VtxXYZQ[%d][%d]",r,c), Form("fCent_VtxXYZQ[%d][%d]",r,c), 50, -0.01, 0.01, 50, -0.010, 0.010, 40, -10, 10, "s");
			  } else if (fDataSet==k2018r) {
				fCent_VtxXYZQ[r][c] = new TProfile3D(Form("fCent_VtxXYZQ[%d][%d]",r,c), Form("fCent_VtxXYZQ[%d][%d]",r,c), 40, 0.08, 0.1, 40, 0.36, 0.38, 40, -10, 10, "s");
			  }
			  fOutputRecenter1->Add(fCent_VtxXYZQ[r][c]);
			}
		  }
		}
	}
    //@Shi Add run by run recentering histograms for ZDC (end)
    
  } //@Shi add for !fUseTowerEq = False (end)
  
  PostData(2, fOutput);
  if (fStepZDCRecenter>=0) {
    PostData(3, fOutputRecenter1);
    PostData(4, fOutputRecenter2);
    PostData(5, fOutputRecenter3);
  }
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  AliMCEvent* McEvent = MCEvent();
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  // AliMultiplicity* myTracklets = NULL;
  // AliESDPmdTrack* pmdtracks = NULL;
  // int availableINslot=1;

  if (!(fCutsRP&&fCutsPOI&&fCutsEvent)) {
    AliError("cuts not set");
    return;
  }

  Int_t numTracks =  aod->GetNumberOfTracks();
  if (numTracks==0) return;

  TObject *head = aod->GetHeader();
  Int_t RunBin=-1, bin=0, RunNum=-1;

  if(!head->InheritsFrom("AliNanoAODStorage")){ //no nanoAOD
    RunNum = aod->GetRunNumber();
    for(Int_t c=0;c<fCRCnRun;c++) {
      if(fRunList[c]==RunNum) RunBin=bin;
      else bin++;
    }
    if(RunBin==-1) return;
    if(fDataSet==kAny) RunBin=0;
  }
  
  //DEFAULT - automatically takes care of everything
  if (fAnalysisType == kAUTOMATIC || fAnalysisType == kTracklets) {

    // get centrality
    Double_t centrV0M=300, centrCL1=300, centrCL0=300, centrTRK=300;
    if(!head->InheritsFrom("AliNanoAODStorage")){
      if(fDataSet!=k2015 && fDataSet!=k2015v6 &&  fDataSet!=k2015pidfix && fDataSet!=k2018r) {  //@Shi test code fDataSet!=2018r
        centrV0M = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
        centrCL1 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");
        centrCL0 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL0");
        centrTRK = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("TRK");
      } else {
        fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
        if(!fMultSelection) {
          //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
          AliWarning("AliMultSelection object not found!");
        } else {
          centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
          centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
          centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
          centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");
        }
      }
    }else{

      AliNanoAODHeader *nanoAodHeader = (AliNanoAODHeader*) head;
      RunNum = nanoAodHeader->GetRunNumber();
      for(Int_t c=0;c<fCRCnRun;c++) {
        if(fRunList[c]==RunNum) RunBin=bin;
        else bin++;
      }

      if(RunBin==-1) return;
      if(fDataSet==kAny) RunBin=0;

      centrV0M = nanoAodHeader->GetCentr("V0M");
      centrTRK = nanoAodHeader->GetCentr("TRK");
      centrCL1 = nanoAodHeader->GetCentr("CL1");
      centrCL0 = nanoAodHeader->GetCentr("CL0");

    }
    //check event cuts
    if (InputEvent()) {
      if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
      if(fRejectPileUp) {
        Bool_t IsPileUp = SelectPileup(aod); //@Shi not really working for 2018
        if(IsPileUp) return;
      }
    }

    //first attach all possible information to the cuts
    fCutsRP->SetEvent( InputEvent(), MCEvent() );  //attach event
    fCutsPOI->SetEvent( InputEvent(), MCEvent() );

    //then make the event
    fFlowEvent->Fill( fCutsRP, fCutsPOI );

    if(fAnalysisType == kTracklets) {
      // fill with tracklets
      AliAODTracklets* anInputTracklets = (AliAODTracklets*)aod->GetTracklets();
      Int_t multSPD = anInputTracklets->GetNumberOfTracklets();
      Double_t BField = aod->GetMagneticField();
      //loop over tracklets
      for (Int_t itracklet=0; itracklet<multSPD; ++itracklet) {
        Float_t thetaTr= anInputTracklets->GetTheta(itracklet);
        Float_t phiTr= anInputTracklets->GetPhi(itracklet);
        // calculate eta
        Float_t etaTr = -TMath::Log(TMath::Tan(thetaTr/2.));
        //make new AliFLowTrackSimple
        AliFlowTrack* pTrack = new AliFlowTrack();
        pTrack->SetPt(0.5);
        pTrack->SetEta(etaTr);
        // set charge: "according to Ruben, with positive magnetic field, the positive tracks rotate clockwise. Since the angle stored in AOD is the difference between the hit on the first and second layers of SPD, in positive mag field, positive delta_phi -> positive track charge."
        Double_t DeltaPhi = anInputTracklets->GetDeltaPhi(itracklet);
        if(BField>0. && DeltaPhi>0.) pTrack->SetCharge(1);
        if(BField>0. && DeltaPhi<0.) pTrack->SetCharge(-1);
        if(BField<0. && DeltaPhi>0.) pTrack->SetCharge(-1);
        if(BField<0. && DeltaPhi<0.) pTrack->SetCharge(1);
        // correction of phi
        if(fCorrectPhiTracklets) {
          phiTr += 39./34.*DeltaPhi;
          if (phiTr < 0.)  phiTr += 2.*TMath::Pi();
          if (phiTr > 2.*TMath::Pi()) phiTr -= 2.*TMath::Pi();
        }
        pTrack->SetPhi(phiTr);
        //marking the particles as POI type 2:
        fFlowEvent->IncrementNumberOfPOIs(2);
        pTrack->SetPOItype(2,kTRUE);
        pTrack->SetSource(AliFlowTrack::kFromTracklet);
        //Add the track to the flowevent
        fFlowEvent->AddTrack(pTrack);
      }
    }

    fFlowEvent->SetReferenceMultiplicity(fCutsEvent->GetReferenceMultiplicity(InputEvent(),McEvent));

    if(fCentrEstimator==kV0M) fFlowEvent->SetCentrality(centrV0M);
    if(fCentrEstimator==kCL0) fFlowEvent->SetCentrality(centrCL0);
    if(fCentrEstimator==kCL1) fFlowEvent->SetCentrality(centrCL1);
    if(fCentrEstimator==kTRK) fFlowEvent->SetCentrality(centrTRK);
    fFlowEvent->SetCentralityCL1(centrCL1);
    fFlowEvent->SetCentralityTRK(centrTRK);
    //   fFlowEvent->SetNITSCL1(((AliVAODHeader*)aod->GetHeader())->GetNumberOfITSClusters(1));

    Double_t SumV0=0.;
    UInt_t period, orbit24;

    if(!head->InheritsFrom("AliNanoAODStorage")){
      for(Int_t i=0; i<64; i++) {
        if(std::isfinite(aod->GetVZEROEqMultiplicity(i))) SumV0 += aod->GetVZEROEqMultiplicity(i);
      }
      period = aod->GetPeriodNumber();
      orbit24 = aod->GetOrbitNumber(); // wrapped down to 24 bits
    }else{
      AliNanoAODHeader *nanoAodHeader = (AliNanoAODHeader*) head;
      SumV0 = nanoAodHeader->GetVar(nanoAodHeader->GetVarIndex("cstV0"));
      period = nanoAodHeader->GetVar(nanoAodHeader->GetVarIndex("cstPeriod"));
      orbit24 = nanoAodHeader->GetVar(nanoAodHeader->GetVarIndex("cstOrbit")); // wrapped down to 24 bits
    }
    fFlowEvent->SetNITSCL1(SumV0);

    // set absolute orbit number

    if (period > 255) { // 8 bits
      period = 255;
      orbit24 = (1<<24) - 1;
    }
    if (orbit24 >= (1<<24)) { // 24 bits
      period = 255;
      orbit24 = (1<<24) - 1;
    }
    UInt_t orbit  = period * (1<<24) + orbit24;
    fFlowEvent->SetAbsOrbit(orbit);

    Double_t vtxpos[3]={0.,0.,0.};
    vtxpos[0] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetX();
    vtxpos[1] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetY();
    vtxpos[2] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetZ();
    fFlowEvent->SetVertexPosition(vtxpos);

    if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);

    // run-by-run QA
    //   for(Int_t jTracks = 0; jTracks<aod->GetNumberOfTracks(); jTracks++){
    //     AliAODTrack* track = (AliAODTrack*)aod->GetTrack(jTracks);
    //     if(!track) continue;
    //     // general kinematic & quality cuts
    //     if (track->Pt() < .2 || track->Pt() > 20. || TMath::Abs(track->Eta()) > .8 || track->GetTPCNcls() < 70)  continue;
    //     if (track->TestFilterBit(128)) fPtPhiEtaRbRFB128[RunBin][CenBin]->Fill(track->Pt(),track->Phi(),track->Eta());
    //     if (track->TestFilterBit(768)) fPtPhiEtaRbRFB768[RunBin][CenBin]->Fill(track->Pt(),track->Phi(),track->Eta());
    //   }
    fCenDis->Fill(centrV0M);

  }

  if (fAnalysisType == kTrackQA) {

    // empty flow event
    fFlowEvent->ClearFast();

    // get centrality
    Double_t centrV0M=300;
    if(fDataSet!=k2015 && fDataSet!=k2015v6 && fDataSet!=k2015pidfix) {
      centrV0M = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
    } else {
      fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      if(!fMultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
      } else {
        centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
      }
    }
    if(centrV0M<10. || centrV0M>40.) return;

    //check event cuts
    if (InputEvent()) {
      if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
      if(fRejectPileUp) {
        Bool_t IsPileUp = SelectPileup(aod);
        if(IsPileUp) return;
      }
    }

    Double_t BField = aod->GetMagneticField();

    for(Int_t jTracks = 0; jTracks<aod->GetNumberOfTracks(); jTracks++){

      AliAODTrack* track = (AliAODTrack*)aod->GetTrack(jTracks);
      if(!track) continue;

      Double_t dPt = track->Pt();
      Double_t dEta = track->Eta();
      Double_t dPhi = track->Phi();
      Double_t dCharge = track->Charge();
      Int_t cw = 0;
      if(BField>0. && dCharge>0.) cw=0;
      if(BField>0. && dCharge<0.) cw=1;
      if(BField<0. && dCharge>0.) cw=2;
      if(BField<0. && dCharge<0.) cw=3;

      // general kinematic cuts
      if (dPt < .2 || dPt > 50. || TMath::Abs(dEta) > 0.8) continue;

      // cut on DCA
      Double_t DCAxy = track->DCA();
      Double_t DCAz = track->ZAtDCA();
      if (std::abs((Int_t)DCAxy)==999 || std::abs((Int_t)DCAz)==999) {
        // re-evaluate the dca as it seems to not be natively present
        // allowed only for tracks inside the beam pipe
        Double_t pos[3] = {-99., -99., -99.};
        track->GetPosition(pos);
        if(pos[0]*pos[0]+pos[1]*pos[1] <= 3.*3.) {
          AliAODTrack copy(*track);       // stack copy
          Double_t b[2] = {-99., -99.};
          Double_t bCov[3] = {-99., -99., -99.};
          if(copy.PropagateToDCA(aod->GetPrimaryVertex(), aod->GetMagneticField(), 100., b, bCov)) {
            DCAxy = b[0];
            DCAz = b[1];
          }
        }
      }
      if(fabs(DCAxy)>2.4 || fabs(DCAz)>3.2) continue;

      // various cuts on TPC clusters
      if (track->GetTPCNcls() < 70) continue;
      Double_t chi2_per_tpc = track->Chi2perNDF();
      if (chi2_per_tpc < 0.1 || chi2_per_tpc > 4.) continue;
      Double_t fraction_shared_tpccls = 1.*track->GetTPCnclsS()/track->GetTPCncls();
      if (fraction_shared_tpccls > 0.4) continue;

      Int_t FBarray[4] = {32,96,128,768};

      // test filter bits
      for(Int_t fb=0; fb<fKNFBs; fb++){
        if (track->TestFilterBit(FBarray[fb])) {
          fTrackQADCAxy[fb][cw]->Fill(dEta,dPt,DCAxy);
          fTrackQADCAz[fb][cw]->Fill(dEta,dPt,DCAz);
          fTrackQApT[fb][cw]->Fill(dEta,dPt);
          fEbEQRe[fb][cw]->Fill(dEta,dPt,TMath::Cos(dPhi));
          fEbEQIm[fb][cw]->Fill(dEta,dPt,TMath::Sin(dPhi));
          fEbEQMu[fb][cw]->Fill(dEta,dPt);
        }
      }
    }

    // compute cos(#Delta#phi)
    for(Int_t bx=1; bx<=fEbEQRe[0][0]->GetXaxis()->GetNbins(); bx++) {
      for(Int_t by=1; by<=fEbEQRe[0][0]->GetYaxis()->GetNbins(); by++) {

        Double_t dEta = fEbEQRe[0][0]->GetXaxis()->GetBinCenter(bx);
        Double_t dPt  = fEbEQRe[0][0]->GetYaxis()->GetBinCenter(by);

        for(Int_t fb=0; fb<fKNFBs; fb++){
          for(Int_t c=0; c<4; c++){

            Double_t QRe = fEbEQRe[fb][c]->GetBinContent(bx,by);
            Double_t QIm = fEbEQIm[fb][c]->GetBinContent(bx,by);
            Double_t M   = 1.*fEbEQMu[fb][c]->GetBinContent(bx,by);

            if(M>1.) {
              Double_t c1 = (QRe*QRe+QIm*QIm-M)/(M*(M-1.));
              fTrackQADphi[fb][c]->Fill(dEta,dPt,c1);
            }
          }
        }

      }
    }

    // reset EbE histograms
    for(Int_t fb=0; fb<fKNFBs; fb++){
      for(Int_t c=0; c<4; c++){
        fEbEQRe[fb][c]->Reset();
        fEbEQIm[fb][c]->Reset();
        fEbEQMu[fb][c]->Reset();
      }
    }

  }

  if (fAnalysisType == kMCAOD) {

    //check event cuts
    if (InputEvent()) {
      if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
      if(fRejectPileUp && fAnalysisUtil->IsPileUpEvent(InputEvent())) return;
    }

    fFlowEvent->ClearFast();

    if(!McEvent) {
      AliError("ERROR: Could not retrieve MCEvent");
      return;
    }
    fStack = (TClonesArray*)aod->FindListObject(AliAODMCParticle::StdBranchName());
    if(!fStack){
      AliError("ERROR: Could not retrieve MCStack");
      return;
    }

    // get centrality (from AliMultSelection or AliCentrality)
    Double_t centr = 300;
    if(fDataSet==k2015 || fDataSet==k2015v6 || fDataSet==k2015pidfix) {
      fMultSelection = (AliMultSelection*)aod->FindListObject("MultSelection");
      if(!fMultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
      } else {
        centr = fMultSelection->GetMultiplicityPercentile("V0M");
      }
    } else {
      centr = (((AliVAODHeader*)aod->GetHeader())->GetCentralityP())->GetCentralityPercentile("V0M");
    }
    // centrality bin
    if (centr<fCentrLowLim || centr>=fCentrUpLim ) return;
    Int_t CenBin = -1;
    CenBin = GetCenBin(centr);
    if(CenBin==-1) return;
    fCenDis->Fill(centr);

    AliMCEvent* McEventFake = NULL;
    //first attach all possible information to the cuts
    fCutsRP->SetEvent( InputEvent(), McEventFake );  //attach event
    fCutsPOI->SetEvent( InputEvent(), McEventFake );

    //then make the event
    fFlowEvent->Fill( fCutsRP, fCutsPOI );

    fFlowEvent->SetCentrality(centr);
    fFlowEvent->SetCentralityCL1(centr);
    fFlowEvent->SetCentralityTRK(centr);
    fFlowEvent->SetReferenceMultiplicity(fFlowEvent->GetNumberOfRPs());

    Double_t SumV0=0.;
    for(Int_t i=0; i<64; i++) {
      if(std::isfinite(aod->GetVZEROEqMultiplicity(i))) SumV0 += aod->GetVZEROEqMultiplicity(i);
    }
    fFlowEvent->SetNITSCL1(SumV0);

    Double_t vtxpos[3]={0.,0.,0.};
    vtxpos[0] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetX();
    vtxpos[1] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetY();
    vtxpos[2] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetZ();
    fFlowEvent->SetVertexPosition(vtxpos);

    // reconstructed
    for(Int_t jTracks = 0; jTracks<aod->GetNumberOfTracks(); jTracks++){

      AliAODTrack* track = (AliAODTrack*)aod->GetTrack(jTracks);
      if(!track) continue;

      // to select primaries
      Int_t lp = TMath::Abs(track->GetLabel());

      // general kinematic cuts
      if (track->Pt() < .2 || track->Pt() > 20. || TMath::Abs(track->Eta()) > 0.8) continue;

      // cut on DCA
      Double_t DCAxy = track->DCA();
      Double_t DCAz = track->ZAtDCA();
      if(fabs(DCAxy)>2.4 || fabs(DCAz)>3.2) continue;

      // various cuts on TPC clusters
      if (track->GetTPCNcls() < 70) continue;
      Double_t chi2_per_tpc = track->Chi2perNDF();
      if (chi2_per_tpc < 0.1 || chi2_per_tpc > 4.) continue;
      Double_t fraction_shared_tpccls = 1.*track->GetTPCnclsS()/track->GetTPCncls();
      if (fraction_shared_tpccls > 0.4) continue;

      // test filter bits
      if (((AliAODMCParticle*)fStack->At(lp))->IsPhysicalPrimary()) {
        if (track->TestFilterBit(32))  fPtSpecFB32[0][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(96))  fPtSpecFB96[0][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(128)) fPtSpecFB128[0][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(768)) fPtSpecFB768[0][CenBin]->Fill(track->Pt());
      } else {
        if (track->TestFilterBit(32))  fPtSpecFB32[1][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(96))  fPtSpecFB96[1][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(128)) fPtSpecFB128[1][CenBin]->Fill(track->Pt());
        if (track->TestFilterBit(768)) fPtSpecFB768[1][CenBin]->Fill(track->Pt());
      }

    }

    // generated (physical primaries)
    for(Int_t jTracks = 0; jTracks<fStack->GetEntriesFast(); jTracks++) {
      AliAODMCParticle *MCpart = (AliAODMCParticle*)fStack->At(jTracks);
      if (!MCpart) {
        printf("ERROR: Could not receive MC track %d\n", jTracks);
        continue;
      }

      // kinematic cuts
      if ( MCpart->Pt() < .2 || MCpart->Pt() > 20. || TMath::Abs(MCpart->Eta()) > .8 ) continue;
      // select charged primaries
      if ( MCpart->Charge() == 0. || !MCpart->IsPhysicalPrimary()) continue;

      fPtSpecGen[0][CenBin]->Fill(MCpart->Pt());
    }

    //    fGenHeader = McEvent->GenEventHeader();
    //    if(fGenHeader) fPythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(fGenHeader);
    //  printf("#reconstructed : %d (rejected from cuts %d), #MC primaries : %d (rejected from cuts %d) \n",AODPOIs,AODbads,MCPrims,MCSecos);
    fFlowEvent->SetReferenceMultiplicity(aod->GetNumberOfTracks());
    fFlowEvent->SetCentrality(centr);
    //    if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);
    fFlowEvent->SetRun(RunNum);
    //  printf("Run : %d, RefMult : %d, Cent : %f \n",fFlowEvent->GetRun(),fFlowEvent->GetReferenceMultiplicity(),fFlowEvent->GetCentrality());
  }

  if(fAnalysisType == kMCESD) {

    fFlowEvent->ClearFast();

    if(!esd) {
      AliError("ERROR: Could not retrieve ESDEvent");
      return;
    }
    if(!McEvent) {
      AliError("ERROR: Could not retrieve MCEvent");
      return;
    }
    AliStack* fStack = fMCEvent->Stack();
    if(!fStack) {
      AliError("ERROR: Could not retrieve MCStack");
      return;
    }

    AliESDVertex *vertex = (AliESDVertex*) esd->GetPrimaryVertex();
    if (!vertex) return;
    if (TMath::Abs(vertex->GetZ()) > 10. ) return;
    if (vertex->GetNContributors() < 1 ) return;
    AliCentrality *centrality = esd->GetCentrality();
    if (!centrality) return;
    Double_t centr = centrality->GetCentralityPercentile("V0M");
    if (centr<fCentrLowLim || centr>=fCentrUpLim ) return;
    Int_t CenBin = -1;
    if (centr>0. && centr<5.) CenBin=0;
    if (centr>5. && centr<10.) CenBin=1;
    if (centr>10. && centr<20.) CenBin=2;
    if (centr>20. && centr<30.) CenBin=3;
    if (centr>30. && centr<40.) CenBin=4;
    if (centr>40. && centr<50.) CenBin=5;
    if (centr>50. && centr<60.) CenBin=6;
    if (centr>60. && centr<70.) CenBin=7;
    if (centr>70. && centr<80.) CenBin=8;
    if (centr>80. && centr<90.) CenBin=9;
    if(CenBin==-1) return;

    //Generated
    Int_t MCPrims = 0;
    for ( Int_t i=0 ; i<fStack->GetNtrack() ; i++ )  {

      //Primaries Selection
      TParticle *particle = (TParticle*)fStack->Particle(i);
      if (!particle) continue;
      if (!fStack->IsPhysicalPrimary(i)) continue;
      if ( particle->GetPDG()->Charge() == 0.) continue;

      //Kinematic Cuts
      if ( particle->Pt()<0.2 || particle->Pt()>10. ) continue;
      if ( TMath::Abs(particle->Eta())>0.8 ) continue;

      fFlowTrack->SetPhi(particle->Phi());
      fFlowTrack->SetEta(particle->Eta());
      fFlowTrack->SetPt(particle->Pt());
      fFlowTrack->SetSource(AliFlowTrack::kFromMC);
      fFlowTrack->SetForRPSelection(kTRUE);
      fFlowEvent->IncrementNumberOfPOIs(0);
      fFlowTrack->SetForPOISelection(kFALSE);
      fFlowEvent->InsertTrack(fFlowTrack);
      MCPrims++;

      fPtSpecGen[0][CenBin]->Fill(particle->Pt());

    }

    //Reconstructed
    Int_t ESDPrims = 0;
    for (Int_t i=0 ; i<esd->GetNumberOfTracks() ; i++)  {

      //Get reconstructed track
      AliVTrack *vtrack = static_cast<AliVTrack*>(esd->GetTrack(i));
      AliESDtrack *track = dynamic_cast<AliESDtrack*>(vtrack);
      if (!track) continue;

      //Primaries selection
      Int_t lp = TMath::Abs(track->GetLabel());
      if (!fStack->IsPhysicalPrimary(lp)) continue;
      TParticle *particle = (TParticle*)fStack->Particle(lp);
      if (!particle) continue;
      if (particle->GetPDG()->Charge() == 0.) continue;

      //   if(!fCutsPOI->PassesESDcuts(track)) continue;

      Bool_t pass = kTRUE;

      if(fCutTPC) {
        //    printf("******* cutting TPC ******** \n");
        UShort_t ntpccls = track->GetTPCNcls();
        Double_t tpcchi2 = track->GetTPCchi2();
        if (tpcchi2<0.2 || tpcchi2 >=4.) {
          //     printf("TPCchi2 : %e %e ",tpcchi2,track->GetTPCchi2Iter1());
          pass=kFALSE;
        }
        if (ntpccls < 70) {
          //     printf("#TPCcluster : %u %u %u %u ",ntpccls,track->GetTPCNclsF(),track->GetTPCNclsFIter1(),track->GetTPCNclsIter1());
          pass=kFALSE;
        }
      }

      Float_t dcaxy=0.0;
      Float_t dcaz=0.0;
      track->GetImpactParameters(dcaxy,dcaz);
      if (dcaxy > 0.3 || dcaz > 0.3) {
        //    printf("DCA : %e %e ",dcaxy,dcaz);
        pass=kFALSE;
      }
      if(!pass) continue;

      //Kinematic Cuts
      if ( track->Pt()<0.2 || track->Pt()>10. ) continue;
      if ( TMath::Abs(track->Eta())>0.8 ) continue;

      fFlowTrack->SetPhi(track->Phi());
      fFlowTrack->SetEta(track->Eta());
      fFlowTrack->SetPt(track->Pt());
      fFlowTrack->SetSource(AliFlowTrack::kFromESD);
      fFlowTrack->SetForRPSelection(kFALSE);
      fFlowTrack->SetForPOISelection(kTRUE);
      fFlowEvent->IncrementNumberOfPOIs(1);
      fFlowEvent->InsertTrack(fFlowTrack);
      ESDPrims++;

    }

    //  printf("#reconstructed : %d , #MC primaries : %d \n",ESDPrims,MCPrims);
    fFlowEvent->SetReferenceMultiplicity(esd->GetNumberOfTracks());
    fFlowEvent->SetCentrality(centr);
    if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);
    fFlowEvent->SetRun(esd->GetRunNumber());
    //  printf("Run : %d, RefMult : %d, Cent : %f \n",fFlowEvent->GetRun(),fFlowEvent->GetReferenceMultiplicity(),fFlowEvent->GetCentrality());

  } // end of if(fAnalysisType ==  kMCESD)

  if(fAnalysisType == kMCkine) {

    fFlowEvent->ClearFast();

    AliInputEventHandler* McHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!McHandler) {
      AliError("ERROR: Could not retrieve MCtruthEventHandler");
      return;
    }
    McEvent = McHandler->MCEvent();
    if(!McEvent) {
      AliError("ERROR: Could not retrieve MC event");
      return;
    }

    Int_t nTracks = McEvent->GetNumberOfTracks();
    //  Int_t nPrimTr = McEvent->GetNumberOfPrimaries();

    //loop over tracks
    for (Int_t itrkN=0; itrkN<nTracks; itrkN++) {
      //get input particle
      AliMCParticle* pParticle = dynamic_cast<AliMCParticle*>(McEvent->GetTrack(itrkN));
      if (!pParticle) continue;

      //check if track passes the cuts
      if (McEvent->IsPhysicalPrimary(itrkN) && pParticle->Charge()!=0) {
        fFlowTrack->Set(pParticle);
        fFlowTrack->SetSource(AliFlowTrack::kFromMC);
        fFlowTrack->SetForRPSelection(kTRUE);
        fFlowEvent->IncrementNumberOfPOIs(0);
        fFlowTrack->SetForPOISelection(kTRUE);
        fFlowEvent->IncrementNumberOfPOIs(1);
        fFlowEvent->InsertTrack(fFlowTrack);
      }
    }// for all tracks

    // if monte carlo event get reaction plane from monte carlo (depends on generator)
    if (McEvent && McEvent->GenEventHeader()) fFlowEvent->SetMCReactionPlaneAngle(McEvent);
    // set reference multiplicity
    fFlowEvent->SetReferenceMultiplicity(McEvent->GetNumberOfTracks());
    // tag subevents
    fFlowEvent->TagSubeventsInEta(fMinA,fMaxA,fMinB,fMaxB);
    // set centrality from impact parameter
    Double_t ImpPar=0., CenPer=0.;
    fGenHeader = McEvent->GenEventHeader();
    if(fGenHeader){
      fPythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader*>(fGenHeader);
      if(fPythiaGenHeader) ImpPar = fPythiaGenHeader->GetImpactParameter();
      fHijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(fGenHeader);
      if(fHijingGenHeader) ImpPar = fHijingGenHeader->ImpactParameter();
      if(ImpPar) CenPer = 0.3859796743103508*pow(ImpPar,2.);
      if(CenPer>0. && CenPer<100.) fFlowEvent->SetCentrality(CenPer);
      else return;
      fFlowEvent->SetRun(1);
    }

  } // end of if(fAnalysisType == kMCkine)

  if (!fFlowEvent) return; //shuts up coverity

  //check final event cuts
  Int_t mult = fFlowEvent->NumberOfTracks();
  //  AliInfo(Form("FlowEvent has %i tracks",mult));
  if (mult<fMinMult || mult>fMaxMult) {
    AliWarning("FlowEvent cut on multiplicity"); return;
  }

  //define dead zone
  fFlowEvent->DefineDeadZone(fExcludedEtaMin, fExcludedEtaMax, fExcludedPhiMin, fExcludedPhiMax );

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////AFTERBURNER
  if (fAfterburnerOn)
  {
    //if reaction plane not set from elsewhere randomize it before adding flow
    if (!fFlowEvent->IsSetMCReactionPlaneAngle())
    fFlowEvent->SetMCReactionPlaneAngle(gRandom->Uniform(0.0,TMath::TwoPi()));

    if(fDifferentialV2)
    fFlowEvent->AddV2(fDifferentialV2);
    else
    fFlowEvent->AddFlow(fV1,fV2,fV3,fV4,fV5);     //add flow
    fFlowEvent->CloneTracks(fNonFlowNumberOfTrackClones); //add nonflow by cloning tracks
  }
  //////////////////////////////////////////////////////////////////////////////

  //tag subEvents
  fFlowEvent->TagSubeventsInEta(fMinA,fMaxA,fMinB,fMaxB);

  //do we want to serve shullfed tracks to everybody?
  fFlowEvent->SetShuffleTracks(fShuffleTracks);

  // associate the mother particles to their daughters in the flow event (if any)
  fFlowEvent->FindDaughters();

  //fListHistos->Print();
  //fOutputFile->WriteObject(fFlowEvent,"myFlowEventSimple");

  //********************************************************************************************************************************

  if(fAnalysisType == kAUTOMATIC || fAnalysisType == kTracklets) {

    // PHYSICS SELECTION
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *hdr = (AliInputEventHandler*)am->GetInputEventHandler();

    if(hdr->IsEventSelected()==0 && !head->InheritsFrom("AliNanoAODStorage")) return;
    //if(hdr->IsEventSelected() && AliVEvent::kAny) {

    Double_t centrperc = fFlowEvent->GetCentrality();
    Int_t cenb = (Int_t)centrperc;

    Int_t nTracklets ;
    if(!head->InheritsFrom("AliNanoAODStorage")){
      AliAODTracklets *trackl = aod->GetTracklets();
      nTracklets = trackl->GetNumberOfTracklets();
    }else{
      AliNanoAODHeader *nanoAodHeader = (AliNanoAODHeader*) head;
      nTracklets = nanoAodHeader->GetVar(nanoAodHeader->GetVarIndex("cstNTrackelets"));
    }


    // get VZERO data
    AliAODVZERO *vzeroAOD = aod->GetVZEROData();
    Double_t multV0A = vzeroAOD->GetMTotV0A();
    Double_t multV0C = vzeroAOD->GetMTotV0C();

    // set VZERO Q-vectors
    if(fDataSet==k2015 || fDataSet==k2015v6) {   //@Shi no fVZEROQVectorRecQxStored and fVZEROQVectorRecQyStored available for 2018
      Int_t CachednRing = 1;
      Double_t QxTot[fkVZEROnHar] = {0.}, QyTot[fkVZEROnHar] = {0.};
      Double_t denom = 0.;
      Double_t V0TotQC[fkVZEROnHar][2] = {{0.}}, V0TotQA[fkVZEROnHar][2] = {{0.}};
      Double_t MultC[fkVZEROnHar] = {0.}, MultA[fkVZEROnHar] = {0.};

      for(Int_t i=0; i<64; i++) {

        // correct multiplicity per channel
        Double_t mult = vzeroAOD->GetMultiplicity(i);
        if(fVZEROGainEqHist) {
          Double_t EqFactor = fVZEROGainEqHist->GetBinContent(RunBin+1,i+1);
          if(EqFactor>0.) mult *= EqFactor;
        }
        fVZEROMult->Fill(RunBin+0.5,i+0.5,mult);

        // build Q-vector per ring
        Int_t nRing = (Int_t)i/8 + 1;
        Double_t ChPhi = TMath::PiOver4()*(0.5+i%8);

        if(i == 63) {
          for (Int_t k=0; k<fkVZEROnHar; k++) {
            QxTot[k] += mult*TMath::Cos((k+1.)*ChPhi);
            QyTot[k] += mult*TMath::Sin((k+1.)*ChPhi);
          }
          denom += mult;
          nRing++;
        }

        if(nRing!=CachednRing && denom!=0) {      //@Shi denom has to be nonzero
          for (Int_t k=0; k<fkVZEROnHar; k++) {
            Double_t QxRec = QxTot[k]/denom;
            Double_t QyRec = QyTot[k]/denom;
            // store values for re-centering
            //            fVZEROQVectorRecQx[k]->Fill(RunBin+0.5,centrperc,CachednRing-0.5,QxRec);
            //            fVZEROQVectorRecQy[k]->Fill(RunBin+0.5,centrperc,CachednRing-0.5,QyRec);
            // do re-centering
            if(fVZEROQVectorRecQxStored[k]) {
              if(!std::isnan(fVZEROQVectorRecQxStored[k]->GetBinContent(fVZEROQVectorRecQxStored[k]->FindBin(RunBin+0.5,centrperc,CachednRing-0.5)))) QxRec -= fVZEROQVectorRecQxStored[k]->GetBinContent(fVZEROQVectorRecQxStored[k]->FindBin(RunBin+0.5,centrperc,CachednRing-0.5));
              if(!std::isnan(fVZEROQVectorRecQyStored[k]->GetBinContent(fVZEROQVectorRecQyStored[k]->FindBin(RunBin+0.5,centrperc,CachednRing-0.5)))) QyRec -= fVZEROQVectorRecQyStored[k]->GetBinContent(fVZEROQVectorRecQyStored[k]->FindBin(RunBin+0.5,centrperc,CachednRing-0.5));
            }
            // sum of Q-vectors over all rings (total V0 Q-vector)
            if (CachednRing >= fMinRingVZC && CachednRing <= fMaxRingVZC) {
              V0TotQC[k][0] += QxRec*denom;
              V0TotQC[k][1] += QyRec*denom;
              MultC[k] += denom;
            }
            if (CachednRing >= fMinRingVZA && CachednRing <= fMaxRingVZA) {
              V0TotQA[k][0] += QxRec*denom;
              V0TotQA[k][1] += QyRec*denom;
              MultA[k] += denom;
            }
            QxTot[k] = 0.;
            QyTot[k] = 0.;
          }
          denom = 0.;
          CachednRing = nRing;
        }
        for (Int_t k=0; k<fkVZEROnHar; k++) {
          QxTot[k] += mult*TMath::Cos((k+1.)*ChPhi);
          QyTot[k] += mult*TMath::Sin((k+1.)*ChPhi);
        }
        denom += mult;
      }

      for (Int_t k=0; k<fkVZEROnHar; k++) {
        if(MultC[k]>0. && MultA[k]>0.) {
          Double_t QCx = V0TotQC[k][0]/MultC[k], QCy = V0TotQC[k][1]/MultC[k], QAx = V0TotQA[k][0]/MultA[k], QAy = V0TotQA[k][1]/MultA[k];
          if(!std::isnan(QCx) && !std::isnan(QCy) && !std::isnan(QAx) && !std::isnan(QAy)) {
            fFlowEvent->SetV02Qsub(QCx,QCy,MultC[k],QAx,QAy,MultA[k],k+1);
            fVZEROQVectorRecFinal[k][0]->Fill(RunBin+0.5,centrperc,QCx);
            fVZEROQVectorRecFinal[k][1]->Fill(RunBin+0.5,centrperc,QCy);
            fVZEROQVectorRecFinal[k][2]->Fill(RunBin+0.5,centrperc,QAx);
            fVZEROQVectorRecFinal[k][3]->Fill(RunBin+0.5,centrperc,QAy);
            fVZEROQVectorRecFinal[k][4]->Fill(RunBin+0.5,centrperc,QCx*QAx);
            fVZEROQVectorRecFinal[k][5]->Fill(RunBin+0.5,centrperc,QCy*QAy);
            fVZEROQVectorRecFinal[k][6]->Fill(RunBin+0.5,centrperc,QCx*QAy);
            fVZEROQVectorRecFinal[k][7]->Fill(RunBin+0.5,centrperc,QCy*QAx);
          } else {
            fFlowEvent->SetV02Qsub(0.,0.,0.,0.,0.,0.,k+1);
          }
        } else {
          fFlowEvent->SetV02Qsub(0.,0.,0.,0.,0.,0.,k+1);
        }
      }
    }

    //      AliAODForwardMult* aodForward = static_cast<AliAODForwardMult*>(aodEvent->FindListObject("Forward"));
    //      const TH2D& d2Ndetadphi = aodForward->GetHistogram();
    //      Int_t nEta = d2Ndetadphi.GetXaxis()->GetNbins();
    //      Int_t nPhi = d2Ndetadphi.GetYaxis()->GetNbins();
    //      Double_t ret = 0.;
    //      // Loop over eta
    //      for (Int_t iEta = 1; iEta <= nEta; iEta++) {
    //        Int_t valid = d2Ndetadphi.GetBinContent(iEta, 0);
    //        if (!valid) continue; // No data expected for this eta
    //        // Loop over phi
    //        for (Int_t iPhi = 1; i <= nPhi; i++) {
    //          ret = d2Ndetadphi.GetBinContent(iEta, iPhi);
    //          printf("eta %e phi %e : %e \n",d2Ndetadphi.GetXaxis()->GetBinCenter(iEta),d2Ndetadphi.GetYaxis()->GetBinCenter(iPhi),ret);
    //        }
    //      }

    AliAODZDC *aodZDC = aod->GetZDCData();

    const Double_t * towZNCraw = aodZDC->GetZNCTowerEnergy();
    const Double_t * towZNAraw = aodZDC->GetZNATowerEnergy();
    
    const Double_t * towZPCraw = aodZDC->GetZPCTowerEnergy(); //@Shi add ZPC
    const Double_t * towZPAraw = aodZDC->GetZPATowerEnergy(); //@Shi add ZPA
	
	Double_t towZNCrawEnergy[5]={0.}, towZNArawEnergy[5]={0.};
	for(Int_t i=0; i<5; i++) {
		towZNCrawEnergy[i] = towZNCraw[i];
		towZNArawEnergy[i] = towZNAraw[i];
	}
	
	fFlowEvent->SetTowZNCraw(towZNCrawEnergy);
	fFlowEvent->SetTowZNAraw(towZNArawEnergy);
    // Get centroid from ZDCs *******************************************************

    Double_t Enucl = (RunNum < 209122 ? 1380. : 2511.);
    Double_t xyZNC[2]={0.,0.}, xyZNA[2]={0.,0.};
    Double_t towZNC[5]={0.}, towZNA[5]={0.};

    Double_t ZNCcalib=1., ZNAcalib=1.;
    
    Double_t xyZPC[2]={0.,0.}, xyZPA[2]={0.,0.}; //@Shi add ZPC ZPA xy
    Double_t towZPC[5]={0.}, towZPA[5]={0.}; //@Shi add towZPC towZPA
    
    if(fUseTowerEq) {
      if(RunNum!=fCachedRunNum) {
        for(Int_t i=0; i<5; i++) {
          fTowerGainEq[0][i] = (TH1D*)(fTowerEqList->FindObject(Form("fZNCTower[%d][%d]",RunNum,i)));
          fTowerGainEq[1][i] = (TH1D*)(fTowerEqList->FindObject(Form("fZNATower[%d][%d]",RunNum,i)));
        }
      }
      
      for(Int_t i=0; i<5; i++) {
        if(fTowerGainEq[0][i]) towZNC[i] = towZNCraw[i]*fTowerGainEq[0][i]->GetBinContent(fTowerGainEq[0][i]->FindBin(centrperc));
        if(fTowerGainEq[1][i]) towZNA[i] = towZNAraw[i]*fTowerGainEq[1][i]->GetBinContent(fTowerGainEq[1][i]->FindBin(centrperc));
        if(fResetNegativeZDC) {
          if(towZNC[i]<0.) towZNC[i] = 0.;
          if(towZNA[i]<0.) towZNA[i] = 0.;
        }
        //@Shi Add test histogram for gain equalization
		fZNenergyBeforeCalibration->Fill(i+0.5,towZNCraw[i]);
		fZNenergyBeforeCalibration->Fill(i+5.5,towZNAraw[i]);
        fZNenergyAfterCalibration->Fill(i+0.5,towZNC[i]);
		fZNenergyAfterCalibration->Fill(i+5.5,towZNA[i]);
      }
    } else { //@Shi uncalibrated
      for(Int_t i=0; i<5; i++) {
        towZNC[i] = towZNCraw[i];
        towZNA[i] = towZNAraw[i];
        towZPC[i] = towZPCraw[i]; //@Shi add towZPC
        towZPA[i] = towZPAraw[i]; //@Shi add towZPA
        if(fResetNegativeZDC) {
          if(towZNC[i]<0.) towZNC[i] = 0.;
          if(towZNA[i]<0.) towZNA[i] = 0.;
          if(towZPC[i]<0.) towZPC[i] = 0.; //@Shi add towZPC
          if(towZPA[i]<0.) towZPA[i] = 0.; //@Shi add towZPA
        }
        fZNCTower[RunBin][i]->Fill(centrperc,towZNC[i]);
        fZNATower[RunBin][i]->Fill(centrperc,towZNA[i]);
        fZPCTower[RunBin][i]->Fill(centrperc,towZPC[i]); //@Shi add filling fZPCTower
        fZPATower[RunBin][i]->Fill(centrperc,towZPA[i]); //@Shi add filling fZPATower
      }
    }

	//@Shi Fill ZP vs ZN correlation using the common channel (begin)
	fZPAvsZNASignal->Fill(towZPA[0]/1000,towZNA[0]/1000,centrperc);
    fZPCvsZNCSignal->Fill(towZPC[0]/1000,towZNC[0]/1000,centrperc);
    //@Shi Fill ZP vs ZN correlation using the common channel (end)

    if(RunNum>=245829 && RunNum<=246994) towZNA[2] = 0.;  //@Shi Not sure why. Added RunNum<=246994 so that this line only applies to LHC15o
    Double_t zncEnergy=0., znaEnergy=0.;
    Double_t zpcEnergy=0., zpaEnergy=0.; //@Shi add zpcEnergy and zpaEnergy
    for(Int_t i=0; i<5; i++){
      zncEnergy += towZNC[i];
      znaEnergy += towZNA[i];
      zpcEnergy += towZPC[i]; //@Shi add zpcEnergy
      zpaEnergy += towZPA[i]; //@Shi add znaEnergy
    }
    if(RunNum>=245829 && RunNum<=246994) znaEnergy *= 8./7.; //@Shi Not sure why. Added RunNum<=246994 so that this line only applies to LHC15o
    fFlowEvent->SetZNCQ0(towZNC[0]);
    fFlowEvent->SetZNAQ0(towZNA[0]);

    Double_t energyZNC =0;
    Double_t energyZNA =0;
    Double_t energyZPC =0;
    Double_t energyZPA =0;
    if(!head->InheritsFrom("AliNanoAODStorage")){
      energyZNC = ((AliVAODHeader*)aod->GetHeader())->GetZDCN1Energy();
      energyZNA = ((AliVAODHeader*)aod->GetHeader())->GetZDCN2Energy();
      energyZPC = ((AliVAODHeader*)aod->GetHeader())->GetZDCP1Energy();
      energyZPA = ((AliVAODHeader*)aod->GetHeader())->GetZDCP2Energy();
    }else{
      AliNanoAODHeader *nanoAodHeader = (AliNanoAODHeader*) head;
      energyZNC = nanoAodHeader->GetVar(nanoAodHeader->GetVarIndex("cstEnergyZNC")); //@Shi always slightly larger than towZNCraw[0]. Why diff?
      energyZNA = nanoAodHeader->GetVar(nanoAodHeader->GetVarIndex("cstEnergyZNA"));
      energyZPC = nanoAodHeader->GetVar(nanoAodHeader->GetVarIndex("cstEnergyZPC"));
      energyZPA = nanoAodHeader->GetVar(nanoAodHeader->GetVarIndex("cstEnergyZPA"));
    }

    fFlowEvent->SetZNCEnergy(energyZNC);
    fFlowEvent->SetZNAEnergy(energyZNA);

    fFlowEvent->SetZPCEnergy(energyZPC);
    fFlowEvent->SetZPAEnergy(energyZPA);

    const Double_t x[4] = {-1.75, 1.75, -1.75, 1.75};
    const Double_t y[4] = {-1.75, -1.75, 1.75, 1.75};
    const Double_t ZDCphi[4] = {TMath::PiOver4(), TMath::PiOver4()*3., TMath::PiOver4()*5., TMath::PiOver4()*7.}; //@Shi 
    Double_t numXZNC=0., numYZNC=0., denZNC=0., cZNC, wZNC, EZNC, SumEZNC=0.;
    Double_t numXZNA=0., numYZNA=0., denZNA=0., cZNA, wZNA, EZNA, SumEZNA=0., BadChOr;
    Bool_t fAllChONZNC=kTRUE, fAllChONZNA=kTRUE;

    Double_t EZPC; //@Shi add some vaiables for ZPC
    Double_t EZPA; //@Shi add some vaiables for ZPA
    
    if (fUseMCCen) {
      for(Int_t i=0; i<4; i++){

        // get energy
        EZNC = towZNC[i+1];
        fhZNSpectra->Fill(centrperc,i+0.5,EZNC);
        //          fhZNSpectraRbR[RunBin]->Fill(centrperc,i+0.5,EZNC);
        if(fUseZDCSpectraCorr && EZNC>0.) { //@Shi fUseZDCSpectraCorr is false by default
          Double_t mu1 = SpecCorMu1[i]->Interpolate(centrperc);
          Double_t mu2 = SpecCorMu2[i]->Interpolate(centrperc);
          Double_t av = SpecCorAv[i]->Interpolate(centrperc);
          Double_t cor1 = SpecCorSi[i]->Interpolate(centrperc);
          EZNC = exp( (log(EZNC) - mu1 + mu2*cor1)/cor1 ) + av;
          fhZNSpectraCor->Fill(centrperc,i+0.5,EZNC);
        }
        if(fUseZDCSpectraCorr && EZNC<=0.) fAllChONZNC=kFALSE;

        SumEZNC += EZNC;

        // build centroid
        if (EZNC < 0) {
			fRecordNegativeEZNC->Fill(1.5);
			EZNC = 0; // @Shi protect negative EZNC value to screw up Power(EZNC, fZDCGainAlpha)
        } else {
			fRecordNegativeEZNC->Fill(0.5);
		}
        
        wZNC = TMath::Power(EZNC, fZDCGainAlpha);
        numXZNC += x[i]*wZNC;
        numYZNC += y[i]*wZNC;
        denZNC += wZNC;
        fhZNSpectraPow->Fill(centrperc,i+0.5,wZNC);

        // get energy
        if(fDataSet==k2015 || fDataSet==k2015v6) { //@Shi bad tower is fixed for 2018
          if(i==1) {
            EZNA = towZNA[0]-towZNA[1]-towZNA[3]-towZNA[4];
            if(fUseBadTowerCalib && fBadTowerCalibHist[cenb]) {
              EZNA = GetBadTowerResp(EZNA, fBadTowerCalibHist[cenb]);
            }
          } else {
            EZNA = towZNA[i+1];
          }
        } else {
          EZNA = towZNA[i+1];
        }
        fhZNSpectra->Fill(centrperc,i+4.5,EZNA);
        //          fhZNSpectraRbR[RunBin]->Fill(centrperc,i+4.5,EZNA);
        if(fUseZDCSpectraCorr && EZNA>0.) {
          Double_t mu1 = SpecCorMu1[i+4]->Interpolate(centrperc);
          Double_t mu2 = SpecCorMu2[i+4]->Interpolate(centrperc);
          Double_t av = SpecCorAv[i+4]->Interpolate(centrperc);
          Double_t cor1 = SpecCorSi[i+4]->Interpolate(centrperc);
          EZNA = exp( (log(EZNA) - mu1 + mu2*cor1)/cor1 ) + av;
          fhZNSpectraCor->Fill(centrperc,i+4.5,EZNA);
        }
        if(fUseZDCSpectraCorr && EZNA<=0.) fAllChONZNA=kFALSE;
        SumEZNA += EZNA;

        // build centroid
        if (EZNA < 0) {
			fRecordNegativeEZNA->Fill(1.5);
			EZNA = 0; // Shi protect negative EZNA value to screw up Power(EZNC, fZDCGainAlpha)
        } else {
			fRecordNegativeEZNA->Fill(0.5);
		}
        
        wZNA = TMath::Power(EZNA, fZDCGainAlpha);
        numXZNA += x[i]*wZNA;
        numYZNA += y[i]*wZNA;
        denZNA += wZNA;
        fhZNSpectraPow->Fill(centrperc,i+4.5,wZNA);
        
        //@Shi add ZP part (begin)////////////
        // get energy for ZPC
        EZPC = towZPC[i+1];
        fhZPSpectra->Fill(centrperc,i+0.5,EZPC);
        // get energy for ZPA
        EZPA = towZPA[i+1];
		fhZPSpectra->Fill(centrperc,i+0.5,EZPA);
        //@Shi add ZP part (end)
      }
      // store distribution for unfolding
      if(RunNum<245829) {
        Double_t recoE = towZNA[0]-towZNA[1]-towZNA[3]-towZNA[4];
        Double_t trueE = towZNA[2];
        fhZNBCCorr->Fill(centrperc,trueE,recoE);
      }
      if(denZNC>0.){
        Double_t nSpecnC = SumEZNC/Enucl;
        cZNC = 1.89358-0.71262/(nSpecnC+0.71789);
        xyZNC[0] = cZNC*numXZNC/denZNC;
        xyZNC[1] = cZNC*numYZNC/denZNC;
        denZNC *= cZNC;
      }
      else{
        xyZNC[0] = xyZNC[1] = 0.;
      }
      if(denZNA>0.){
        Double_t nSpecnA = SumEZNA/Enucl;
        cZNA = 1.89358-0.71262/(nSpecnA+0.71789);
        xyZNA[0] = cZNA*numXZNA/denZNA;
        xyZNA[1] = cZNA*numYZNA/denZNA;
        denZNA *= cZNA;
      }
      else{
        xyZNA[0] = xyZNA[1] = 0.;
      }
    } else {
      for(Int_t i=0; i<4; i++) {
        if(towZNC[i+1]>0.) {
          wZNC = TMath::Power(towZNC[i+1], fZDCGainAlpha);
          numXZNC += x[i]*wZNC;
          numYZNC += y[i]*wZNC;
          denZNC += wZNC;
        }
        if(towZNA[i+1]>0.) {
          wZNA = TMath::Power(towZNA[i+1], fZDCGainAlpha);
          numXZNA += x[i]*wZNA;
          numYZNA += y[i]*wZNA;
          denZNA += wZNA;
        }
      }
      if(denZNC!=0) {
        xyZNC[0] = numXZNC/denZNC;
        xyZNC[1] = numYZNC/denZNC;
      }
      else{
        xyZNC[0] = xyZNC[1] = 999.;
        zncEnergy = 0.;
      }
      if(denZNA!=0) {
        xyZNA[0] = numXZNA/denZNA;
        xyZNA[1] = numYZNA/denZNA;
      }
      else{
        xyZNA[0] = xyZNA[1] = 999.;
        znaEnergy = 0.;
      }
    }

    if(!fAllChONZNC) denZNC=-1.;
    if(!fAllChONZNA) denZNA=-1.;

    if(denZNC>0. && pow(xyZNC[0]*xyZNC[0]+xyZNC[1]*xyZNC[1],0.5)>1.E-6) fhZNCenDis[0]->Fill(centrperc,xyZNC[0],xyZNC[1]);
    if(denZNA>0. && pow(xyZNA[0]*xyZNA[0]+xyZNA[1]*xyZNA[1],0.5)>1.E-6) fhZNCenDis[1]->Fill(centrperc,-xyZNA[0], xyZNA[1]);
    
    //@Shi fill run by run ZN centroid vs. centrality
    if(fFillZNCenDisRbR) {
		if(denZNC>0. && pow(xyZNC[0]*xyZNC[0]+xyZNC[1]*xyZNC[1],0.5)>1.E-6) fhZNCenDisRbR[RunBin][0]->Fill(centrperc,xyZNC[0],xyZNC[1]);
		if(denZNA>0. && pow(xyZNA[0]*xyZNA[0]+xyZNA[1]*xyZNA[1],0.5)>1.E-6) fhZNCenDisRbR[RunBin][1]->Fill(centrperc,-xyZNA[0], xyZNA[1]);
	}

    fFlowEvent->SetZDC2Qsub(xyZNC,denZNC,xyZNA,denZNA);
    // ******************************************************************************
    // @shi add the code to fill the histograms necessary for recentering
    //const Int_t fnCentBinForRecentering = 20; // this means that a wider centrality bin is used {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100}
    Int_t CentBin = Int_t(centrperc/(100/fnCentBinForRecentering)); // 0 to 19
    if (centrperc == 100) CentBin = 19; // centrality cannot be larger than 100 and when it is exactly 100, use bin 95-100
    
    // Step 0; The profiles for recentring step 1 are filled
    Double_t vtxpos[3]={0.,0.,0.};
    vtxpos[0] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetX();
    vtxpos[1] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetY();
    vtxpos[2] = ((AliAODVertex*)aod->GetPrimaryVertex())->GetZ();
    
    if (fStepZDCRecenter >= 0){
		fAve_VtxX->Fill(RunBin, vtxpos[0]);
		fAve_VtxY->Fill(RunBin, vtxpos[1]);
		fAve_VtxZ->Fill(RunBin, vtxpos[2]);
	}
    // re-centered around zer (implemented only for run2)
	if(fDataSet==k2015 || fDataSet==k2015v6 || fDataSet==k2015pidfix) {
		if(fAvVtxPosX[RunBin]) fVtxPosCor[0] = vtxpos[0]-fAvVtxPosX[RunBin];
		if(fAvVtxPosY[RunBin]) fVtxPosCor[1] = vtxpos[1]-fAvVtxPosY[RunBin];
		if(fAvVtxPosZ[RunBin]) fVtxPosCor[2] = vtxpos[2]-fAvVtxPosZ[RunBin];
	} else {
		fVtxPosCor[0] = vtxpos[0];
		fVtxPosCor[1] = vtxpos[1];
		fVtxPosCor[2] = vtxpos[2];
	}
  
    if (fStepZDCRecenter >= 0){
      if(fStoreCalibZDCRecenter){
		fRun_VtxXQPreCalib[RunBin][0]->Fill(fVtxPosCor[0],xyZNC[0]);
		fRun_VtxXQPreCalib[RunBin][1]->Fill(fVtxPosCor[0],xyZNC[1]);
		fRun_VtxXQPreCalib[RunBin][2]->Fill(fVtxPosCor[0],-xyZNA[0]);
		fRun_VtxXQPreCalib[RunBin][3]->Fill(fVtxPosCor[0],xyZNA[1]);
		
		fRun_VtxYQPreCalib[RunBin][0]->Fill(fVtxPosCor[1],xyZNC[0]);
		fRun_VtxYQPreCalib[RunBin][1]->Fill(fVtxPosCor[1],xyZNC[1]);
		fRun_VtxYQPreCalib[RunBin][2]->Fill(fVtxPosCor[1],-xyZNA[0]);
		fRun_VtxYQPreCalib[RunBin][3]->Fill(fVtxPosCor[1],xyZNA[1]);
		
		fRun_VtxZQPreCalib[RunBin][0]->Fill(fVtxPosCor[2],xyZNC[0]);
		fRun_VtxZQPreCalib[RunBin][1]->Fill(fVtxPosCor[2],xyZNC[1]);
		fRun_VtxZQPreCalib[RunBin][2]->Fill(fVtxPosCor[2],-xyZNA[0]);
		fRun_VtxZQPreCalib[RunBin][3]->Fill(fVtxPosCor[2],xyZNA[1]);
		
        //FillProfiles(&fRun_VtxXQPreCalib,fEvInfo.fRunNum, fEvInfo.fVtxX);
        //FillProfiles(&fRun_VtxYQPreCalib,fEvInfo.fRunNum, fEvInfo.fVtxY);
        //FillProfiles(&fRun_VtxZQPreCalib,fEvInfo.fRunNum, fEvInfo.fVtxZ);
      }
      fRun_CentQ[RunBin][0]->Fill(centrperc, xyZNC[0]); // 1% interval for centrality
      fRun_CentQ[RunBin][1]->Fill(centrperc, xyZNC[1]); // 1% interval for centrality
      fRun_CentQ[RunBin][2]->Fill(centrperc, xyZNA[0]); // 1% interval for centrality
      fRun_CentQ[RunBin][3]->Fill(centrperc, xyZNA[1]); // 1% interval for centrality 0-100 100 bins
      //FillProfiles(&fRun_CentQ, fEvInfo.fRunNum, fEvInfo.fCent);
    }
    
    // Step 1; The first recentring step is performed. Afterwards the
    // profiles for recentring step 2 are filled.

    if (fStepZDCRecenter >= 1) {
	  //auto means = GetMeansProfiles(fAvr_Run_CentQ, fEvInfo.fCent);
      //SubtractMeanFromQVectors(means, fEvInfo.fRunNum);
      
      // Load Calib hist
      if (fStepZDCRecenter < 3) {
        fAvr_Run_CentQ[0] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fRun_CentQ[%d][%d]",RunNum,0)));
        fAvr_Run_CentQ[1] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fRun_CentQ[%d][%d]",RunNum,1)));
        fAvr_Run_CentQ[2] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fRun_CentQ[%d][%d]",RunNum,2)));
        fAvr_Run_CentQ[3] = (TProfile*)(fZDCCalibList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fRun_CentQ[%d][%d]",RunNum,3)));
      } else if (fStepZDCRecenter >= 3) {
		fAvr_Run_CentQ[0] = (TProfile*)(fZDCCalibListStep3RunByRun->FindObject(Form("Run %d",RunNum))->FindObject(Form("fRun_CentQ[%d][%d]",RunNum,0)));
        fAvr_Run_CentQ[1] = (TProfile*)(fZDCCalibListStep3RunByRun->FindObject(Form("Run %d",RunNum))->FindObject(Form("fRun_CentQ[%d][%d]",RunNum,1)));
        fAvr_Run_CentQ[2] = (TProfile*)(fZDCCalibListStep3RunByRun->FindObject(Form("Run %d",RunNum))->FindObject(Form("fRun_CentQ[%d][%d]",RunNum,2)));
        fAvr_Run_CentQ[3] = (TProfile*)(fZDCCalibListStep3RunByRun->FindObject(Form("Run %d",RunNum))->FindObject(Form("fRun_CentQ[%d][%d]",RunNum,3)));
	  }
      if (fAvr_Run_CentQ[0]) {
        Double_t AvQCRe = fAvr_Run_CentQ[0]->GetBinContent(fAvr_Run_CentQ[0]->FindBin(centrperc));
        Double_t AvQCIm = fAvr_Run_CentQ[1]->GetBinContent(fAvr_Run_CentQ[1]->FindBin(centrperc));

        Double_t AvQARe = fAvr_Run_CentQ[2]->GetBinContent(fAvr_Run_CentQ[2]->FindBin(centrperc));
        Double_t AvQAIm = fAvr_Run_CentQ[3]->GetBinContent(fAvr_Run_CentQ[3]->FindBin(centrperc));

        if (AvQCRe && AvQCIm && sqrt(xyZNC[0]*xyZNC[0]+xyZNC[1]*xyZNC[1])>1.E-6) {
          xyZNC[0] = xyZNC[0]-AvQCRe;
          xyZNC[1] = xyZNC[1]-AvQCIm;
        }
        
        if(AvQARe && AvQAIm && sqrt(xyZNA[0]*xyZNA[0]+xyZNA[1]*xyZNA[1])>1.E-6) {
          xyZNA[0] = xyZNA[0]-AvQARe;
          xyZNA[1] = xyZNA[1]-AvQAIm;
        }
      }
      
      if (fStoreCalibZDCRecenter) {
		fRun_CentQCalib[RunBin][0]->Fill(centrperc, xyZNC[0]); 
        fRun_CentQCalib[RunBin][1]->Fill(centrperc, xyZNC[1]); 
        fRun_CentQCalib[RunBin][2]->Fill(centrperc, -xyZNA[0]); 
        fRun_CentQCalib[RunBin][3]->Fill(centrperc, xyZNA[1]);
        
        fRun_VtxXQCalibStep1[RunBin][0]->Fill(fVtxPosCor[0],xyZNC[0]);
		fRun_VtxXQCalibStep1[RunBin][1]->Fill(fVtxPosCor[0],xyZNC[1]);
		fRun_VtxXQCalibStep1[RunBin][2]->Fill(fVtxPosCor[0],-xyZNA[0]);
		fRun_VtxXQCalibStep1[RunBin][3]->Fill(fVtxPosCor[0],xyZNA[1]);
		
		fRun_VtxYQCalibStep1[RunBin][0]->Fill(fVtxPosCor[1],xyZNC[0]);
		fRun_VtxYQCalibStep1[RunBin][1]->Fill(fVtxPosCor[1],xyZNC[1]);
		fRun_VtxYQCalibStep1[RunBin][2]->Fill(fVtxPosCor[1],-xyZNA[0]);
		fRun_VtxYQCalibStep1[RunBin][3]->Fill(fVtxPosCor[1],xyZNA[1]);
		
		fRun_VtxZQCalibStep1[RunBin][0]->Fill(fVtxPosCor[2],xyZNC[0]);
		fRun_VtxZQCalibStep1[RunBin][1]->Fill(fVtxPosCor[2],xyZNC[1]);
		fRun_VtxZQCalibStep1[RunBin][2]->Fill(fVtxPosCor[2],-xyZNA[0]);
		fRun_VtxZQCalibStep1[RunBin][3]->Fill(fVtxPosCor[2],xyZNA[1]);
		//FillProfiles(&fRun_CentQCalib, fEvInfo.fRunNum, fEvInfo.fCent);
	  }
	  
      fCent_VtxXYZQ[CentBin][0]->Fill(fVtxPosCor[0], fVtxPosCor[1], fVtxPosCor[2], xyZNC[0]); // a TProfile3D Tprofile
      fCent_VtxXYZQ[CentBin][1]->Fill(fVtxPosCor[0], fVtxPosCor[1], fVtxPosCor[2], xyZNC[1]); 
      fCent_VtxXYZQ[CentBin][2]->Fill(fVtxPosCor[0], fVtxPosCor[1], fVtxPosCor[2], xyZNA[0]);
      fCent_VtxXYZQ[CentBin][3]->Fill(fVtxPosCor[0], fVtxPosCor[1], fVtxPosCor[2], xyZNA[1]); // cent 0-100 20 bins, vertex X 0.08-0.1 40 bins, vertex Y 0.36-0.38 40 bins, vertex Z -10-10 40 bins
      //FillProfiles(&fCent_VtxXYZQ, fEvInfo.fCentBin, fEvInfo.fVtxX, fEvInfo.fVtxY, fEvInfo.fVtxZ);
    }

    // Step 2; The second recentring step is performed. Afterwards the
    // profiles for recentring step 3 are filled.

    if (fStepZDCRecenter >= 2) {
      //auto means = GetMeansProfiles(fAvr_Cent_VtxXYZQ, fEvInfo.fVtxX, fEvInfo.fVtxY, fEvInfo.fVtxZ);
      //SubtractMeanFromQVectors(means, fEvInfo.fCentBin);
      Bool_t withinvtx = kTRUE;
      
      if (fStepZDCRecenter < 3) { // if the step is 3, the calib file has to be splitted to run-by-run
        for(Int_t k=0; k<4; k++) {
          fAvr_Cent_VtxXYZQ[k] = (TProfile3D*)(fZDCCalibList->FindObject(Form("fCent_VtxXYZQ[%d][%d]",CentBin,k)));
        }
      } else if (fStepZDCRecenter >= 3) { // at step 3, the calib file is loaded separately, pass the run independent calib using fZDCCalibListStep3CommonPart
		for(Int_t k=0; k<4; k++) {
          fAvr_Cent_VtxXYZQ[k] = (TProfile3D*)(fZDCCalibListStep3CommonPart->FindObject(Form("fCent_VtxXYZQ[%d][%d]",CentBin,k)));
        }
	  }
      
      if(fVtxPosCor[0] < fAvr_Cent_VtxXYZQ[0]->GetXaxis()->GetXmin() || fVtxPosCor[0] > fAvr_Cent_VtxXYZQ[0]->GetXaxis()->GetXmax()) withinvtx = kFALSE;
      if(fVtxPosCor[1] < fAvr_Cent_VtxXYZQ[0]->GetYaxis()->GetXmin() || fVtxPosCor[1] > fAvr_Cent_VtxXYZQ[0]->GetYaxis()->GetXmax()) withinvtx = kFALSE;
      if(fVtxPosCor[2] < fAvr_Cent_VtxXYZQ[0]->GetZaxis()->GetXmin() || fVtxPosCor[2] > fAvr_Cent_VtxXYZQ[0]->GetZaxis()->GetXmax()) withinvtx = kFALSE;

      if(withinvtx) {
        xyZNC[0] -= fAvr_Cent_VtxXYZQ[0]->GetBinContent(fAvr_Cent_VtxXYZQ[0]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        xyZNC[1] -= fAvr_Cent_VtxXYZQ[1]->GetBinContent(fAvr_Cent_VtxXYZQ[1]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));

        xyZNA[0] -= fAvr_Cent_VtxXYZQ[2]->GetBinContent(fAvr_Cent_VtxXYZQ[2]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));
        xyZNA[1] -= fAvr_Cent_VtxXYZQ[3]->GetBinContent(fAvr_Cent_VtxXYZQ[3]->FindBin(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]));

      } else {
        Double_t vx = fVtxPosCor[0];
        Double_t vy = fVtxPosCor[1];
        Double_t vz = fVtxPosCor[2];
        if(fVtxPosCor[0] < fAvr_Cent_VtxXYZQ[0]->GetXaxis()->GetXmin()) vx = fAvr_Cent_VtxXYZQ[0]->GetXaxis()->GetBinCenter(1);
        if(fVtxPosCor[0] > fAvr_Cent_VtxXYZQ[0]->GetXaxis()->GetXmax()) vx = fAvr_Cent_VtxXYZQ[0]->GetXaxis()->GetBinCenter(fAvr_Cent_VtxXYZQ[0]->GetNbinsX());
        if(fVtxPosCor[1] < fAvr_Cent_VtxXYZQ[0]->GetYaxis()->GetXmin()) vy = fAvr_Cent_VtxXYZQ[0]->GetYaxis()->GetBinCenter(1);
        if(fVtxPosCor[1] > fAvr_Cent_VtxXYZQ[0]->GetYaxis()->GetXmax()) vy = fAvr_Cent_VtxXYZQ[0]->GetYaxis()->GetBinCenter(fAvr_Cent_VtxXYZQ[0]->GetNbinsY());
        if(fVtxPosCor[2] < fAvr_Cent_VtxXYZQ[0]->GetYaxis()->GetXmin()) vz = fAvr_Cent_VtxXYZQ[0]->GetZaxis()->GetBinCenter(1);
        if(fVtxPosCor[2] > fAvr_Cent_VtxXYZQ[0]->GetYaxis()->GetXmax()) vz = fAvr_Cent_VtxXYZQ[0]->GetZaxis()->GetBinCenter(fAvr_Cent_VtxXYZQ[0]->GetNbinsZ());
		  
        xyZNC[0] -= fAvr_Cent_VtxXYZQ[0]->GetBinContent(fAvr_Cent_VtxXYZQ[0]->FindBin(vx,vy,vz));
        xyZNC[1] -= fAvr_Cent_VtxXYZQ[1]->GetBinContent(fAvr_Cent_VtxXYZQ[1]->FindBin(vx,vy,vz));

        xyZNA[0] -= fAvr_Cent_VtxXYZQ[2]->GetBinContent(fAvr_Cent_VtxXYZQ[2]->FindBin(vx,vy,vz));
        xyZNA[1] -= fAvr_Cent_VtxXYZQ[3]->GetBinContent(fAvr_Cent_VtxXYZQ[3]->FindBin(vx,vy,vz));
      }
      
      if(fStoreCalibZDCRecenter){
		fRun_VtxXQCalibStep2[RunBin][0]->Fill(fVtxPosCor[0],xyZNC[0]);
		fRun_VtxXQCalibStep2[RunBin][1]->Fill(fVtxPosCor[0],xyZNC[1]);
		fRun_VtxXQCalibStep2[RunBin][2]->Fill(fVtxPosCor[0],-xyZNA[0]);
		fRun_VtxXQCalibStep2[RunBin][3]->Fill(fVtxPosCor[0],xyZNA[1]);
		
		fRun_VtxYQCalibStep2[RunBin][0]->Fill(fVtxPosCor[1],xyZNC[0]);
		fRun_VtxYQCalibStep2[RunBin][1]->Fill(fVtxPosCor[1],xyZNC[1]);
		fRun_VtxYQCalibStep2[RunBin][2]->Fill(fVtxPosCor[1],-xyZNA[0]);
		fRun_VtxYQCalibStep2[RunBin][3]->Fill(fVtxPosCor[1],xyZNA[1]);
		
		fRun_VtxZQCalibStep2[RunBin][0]->Fill(fVtxPosCor[2],xyZNC[0]);
		fRun_VtxZQCalibStep2[RunBin][1]->Fill(fVtxPosCor[2],xyZNC[1]);
		fRun_VtxZQCalibStep2[RunBin][2]->Fill(fVtxPosCor[2],-xyZNA[0]);
		fRun_VtxZQCalibStep2[RunBin][3]->Fill(fVtxPosCor[2],xyZNA[1]);
      }

      fRun_VtxXYZQ[RunBin][0]->Fill(fVtxPosCor[0], fVtxPosCor[1], fVtxPosCor[2], xyZNC[0]); // a TProfile3D Tprofile
      fRun_VtxXYZQ[RunBin][1]->Fill(fVtxPosCor[0], fVtxPosCor[1], fVtxPosCor[2], xyZNC[1]);
      fRun_VtxXYZQ[RunBin][2]->Fill(fVtxPosCor[0], fVtxPosCor[1], fVtxPosCor[2], xyZNA[0]);
      fRun_VtxXYZQ[RunBin][3]->Fill(fVtxPosCor[0], fVtxPosCor[1], fVtxPosCor[2], xyZNA[1]);
      
      //FillProfiles(&fRun_VtxXYZQ, fEvInfo.fRunNum, fEvInfo.fVtxX, fEvInfo.fVtxY, fEvInfo.fVtxZ);
    }

    // Step 3; The third recentring step is performed.

    if (fStepZDCRecenter >= 3) {
      //auto mean_step2 = GetMeansProfiles(fAvr_Run_VtxXYZQ, fEvInfo.fVtxX, fEvInfo.fVtxY, fEvInfo.fVtxZ);
      //SubtractMeanFromQVectors(mean_step2, fEvInfo.fRunNum);
      for(Int_t k=0; k<4; k++) {
		fAvr_Run_VtxXYZQ[k] = (TProfile3D*)(fZDCCalibListStep3RunByRun->FindObject(Form("Run %d",RunNum))->FindObject(Form("fRun_VtxXYZQ[%d][%d]",RunNum,k)));
	  }
      
      // check if possible to interpolate
      Bool_t bInterp = kTRUE;
      Int_t bx = fAvr_Run_VtxXYZQ[0]->GetXaxis()->FindBin(fVtxPosCor[0]);
      Int_t by = fAvr_Run_VtxXYZQ[0]->GetYaxis()->FindBin(fVtxPosCor[1]);
      Int_t bz = fAvr_Run_VtxXYZQ[0]->GetZaxis()->FindBin(fVtxPosCor[2]);
      if(bx==1 || bx==fAvr_Run_VtxXYZQ[0]->GetXaxis()->GetNbins()) bInterp = kFALSE;
      if(by==1 || by==fAvr_Run_VtxXYZQ[0]->GetYaxis()->GetNbins()) bInterp = kFALSE;
      if(bz==1 || bz==fAvr_Run_VtxXYZQ[0]->GetZaxis()->GetNbins()) bInterp = kFALSE;
      if(bInterp) {
        xyZNC[0] -= fAvr_Run_VtxXYZQ[0]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
        xyZNC[1] -= fAvr_Run_VtxXYZQ[1]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
        xyZNA[0] -= fAvr_Run_VtxXYZQ[2]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
        xyZNA[1] -= fAvr_Run_VtxXYZQ[3]->Interpolate(fVtxPosCor[0],fVtxPosCor[1],fVtxPosCor[2]);
      } else {
        xyZNC[0] -= fAvr_Run_VtxXYZQ[0]->GetBinContent(bx,by,bz);
        xyZNC[1] -= fAvr_Run_VtxXYZQ[1]->GetBinContent(bx,by,bz);
        xyZNA[0] -= fAvr_Run_VtxXYZQ[2]->GetBinContent(bx,by,bz);
        xyZNA[1] -= fAvr_Run_VtxXYZQ[3]->GetBinContent(bx,by,bz);
      }
    
      
      if (fStoreCalibZDCRecenter){
		fRun_VtxXQCalib[RunBin][0]->Fill(fVtxPosCor[0], xyZNC[0]);
		fRun_VtxXQCalib[RunBin][1]->Fill(fVtxPosCor[0], xyZNC[1]);
		fRun_VtxXQCalib[RunBin][2]->Fill(fVtxPosCor[0], xyZNA[0]);
		fRun_VtxXQCalib[RunBin][3]->Fill(fVtxPosCor[0], xyZNA[1]);
		
		fRun_VtxYQCalib[RunBin][0]->Fill(fVtxPosCor[1], xyZNC[0]);
		fRun_VtxYQCalib[RunBin][1]->Fill(fVtxPosCor[1], xyZNC[1]);
		fRun_VtxYQCalib[RunBin][2]->Fill(fVtxPosCor[1], xyZNA[0]);
		fRun_VtxYQCalib[RunBin][3]->Fill(fVtxPosCor[1], xyZNA[1]);
		
		fRun_VtxZQCalib[RunBin][0]->Fill(fVtxPosCor[2], xyZNC[0]);
		fRun_VtxZQCalib[RunBin][1]->Fill(fVtxPosCor[2], xyZNC[1]);
		fRun_VtxZQCalib[RunBin][2]->Fill(fVtxPosCor[2], xyZNA[0]);
		fRun_VtxZQCalib[RunBin][3]->Fill(fVtxPosCor[2], xyZNA[1]);
		
		fRun_CentQCalib2[RunBin][0]->Fill(centrperc, xyZNC[0]); 
        fRun_CentQCalib2[RunBin][1]->Fill(centrperc, xyZNC[1]); 
        fRun_CentQCalib2[RunBin][2]->Fill(centrperc, -xyZNA[0]); 
        fRun_CentQCalib2[RunBin][3]->Fill(centrperc, xyZNA[1]);
        //FillProfiles(&fRun_VtxXQCalib,fEvInfo.fRunNum, fEvInfo.fVtxX);
        //FillProfiles(&fRun_VtxYQCalib,fEvInfo.fRunNum, fEvInfo.fVtxY);
        //FillProfiles(&fRun_VtxZQCalib,fEvInfo.fRunNum, fEvInfo.fVtxZ);
        //FillProfiles(&fRun_CentQCalib2, fEvInfo.fRunNum, fEvInfo.fCent);
      }
    }
    
    if (fStepZDCRecenter >=0 && fStoreCalibZDCRecenter){
      fCorrQAReCRe->Fill(centrperc,-xyZNA[0]*xyZNC[0]); // -QAReR*QCReR
      fCorrQAReCIm->Fill(centrperc,-xyZNA[0]*xyZNC[1]); // -QAReR*QCImR
      fCorrQAImCRe->Fill(centrperc,xyZNA[1]*xyZNC[0]); // QAImR*QCReR
      fCorrQAImCIm->Fill(centrperc,xyZNA[1]*xyZNC[1]); // QAImR*QCImR
    } 
    
    // ******************************************************************************

    for(int i=0; i<5; i++){
      fhZNCPM[i]->Fill(towZNC[i]);
      if((i<4) && (towZNC[0]>0.)) fhZNCPMQiPMC[i]->Fill(towZNC[i+1]/towZNC[0]);
    }
    for(int i=0; i<5; i++){
      fhZNAPM[i]->Fill(towZNA[i]);
      if(((i<4) && towZNA[0]>0.)) fhZNAPMQiPMC[i]->Fill(towZNA[i+1]/towZNA[0]);
    }
    
    //@Shi add fhZPCPM fhZPAPM (begin)
    for(int i=0; i<5; i++){
      fhZPCPM[i]->Fill(towZPC[i]);
      if((i<4) && (towZPC[0]>0.)) fhZPCPMQiPMC[i]->Fill(towZPC[i+1]/towZPC[0]);
    }
    for(int i=0; i<5; i++){
      fhZPAPM[i]->Fill(towZPA[i]);
      if(((i<4) && towZPA[0]>0.)) fhZPAPMQiPMC[i]->Fill(towZPA[i+1]/towZPA[0]);
    }
	//@Shi add fhZPCPM fhZPAPM (end)
	
    fhZNCvsZNA->Fill(energyZNA, energyZNC);
    fhZDCCvsZDCCA->Fill(energyZNA+energyZPA, energyZNC+energyZPC);
    fhZNCvsZPC->Fill(energyZPC, energyZNC);
    fhZNAvsZPA->Fill(energyZPA, energyZNA);
    fhZNvsZP->Fill(energyZPA+energyZPC, energyZNA+energyZNC);
    fhZNvsVZERO->Fill(multV0A+multV0C, energyZNC+energyZNA);
    fhZDCvsVZERO->Fill(multV0A+multV0C, energyZNA+energyZPA+energyZNC+energyZPC);

    Double_t asymmetry = -999.;
    if((energyZNC+energyZNA)>0.) asymmetry = (energyZNC-energyZNA)/(energyZNC+energyZNA);
    fhAsymm->Fill(asymmetry);
    fhZNAvsAsymm->Fill(asymmetry, energyZNA);
    fhZNCvsAsymm->Fill(asymmetry, energyZNC);

    fhZNCvscentrality->Fill(centrperc, energyZNC);
    fhZNAvscentrality->Fill(centrperc, energyZNA);
    fhZPCvscentrality->Fill(centrperc, energyZPC);
    fhZPAvscentrality->Fill(centrperc, energyZPA);

    // } // PHYSICS SELECTION

  }

  // p) cache run number
  fCachedRunNum = fFlowEvent->GetRun();

  //  printf("debug: NoRPs %e, NoPOIs %e, RunNum %d, Cen %e \n",fFlowEvent->GetNumberOfRPs(),fFlowEvent->GetNumberOfPOIs(),fCachedRunNum,fFlowEvent->GetCentrality());

  PostData(1, fFlowEvent);

  PostData(2, fOutput);
  if (fStepZDCRecenter >= 0) {
    PostData(3, fOutputRecenter1);
    PostData(4, fOutputRecenter2);
    PostData(5, fOutputRecenter3);
  }
}
//________________________________________________________________________

Bool_t AliAnalysisTaskCRCZDC::SelectPileup(AliAODEvent *aod)
{
  Bool_t BisPileup=kFALSE;

  TObject *head =aod->GetHeader();
  if (head->InheritsFrom("AliNanoAODStorage")){

    AliNanoAODHeader * nanohead = (AliNanoAODHeader*)head;
    Int_t pileupIndex = nanohead->GetVarIndex("cstPileUp");
    if (nanohead->GetVar(pileupIndex)==0) BisPileup=kFALSE;
    if (nanohead->GetVar(pileupIndex)==1) BisPileup=kTRUE;

  } else {

    Double_t centrV0M=300., centrCL1=300.;

    if(fDataSet!=k2015 && fDataSet!=k2015v6 && fDataSet!=k2015pidfix && fDataSet!=k2018r) { //@shi add 2018 

      // pileup for LHC10h and LHC11h

      centrV0M = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
      centrCL1 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");

      // check anyway pileup
      if (plpMV(aod)) {
        fPileUpCount->Fill(0.5);
        BisPileup=kTRUE;
      }

      Short_t isPileup = aod->IsPileupFromSPD(3);
      if (isPileup != 0) {
        fPileUpCount->Fill(1.5);
        BisPileup=kTRUE;
      }

      if (((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
        fPileUpCount->Fill(2.5);
        BisPileup=kTRUE;
      }

      if (aod->IsIncompleteDAQ())  {
        fPileUpCount->Fill(3.5);
        BisPileup=kTRUE;
      }

      // check vertex consistency
      const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
      const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();

      if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1)  {
        fPileUpCount->Fill(5.5);
        BisPileup=kTRUE;
      }

      double covTrc[6], covSPD[6];
      vtTrc->GetCovarianceMatrix(covTrc);
      vtSPD->GetCovarianceMatrix(covSPD);

      double dz = vtTrc->GetZ() - vtSPD->GetZ();

      double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
      double errTrc = TMath::Sqrt(covTrc[5]);
      double nsigTot = dz/errTot;
      double nsigTrc = dz/errTrc;

      if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
        fPileUpCount->Fill(6.5);
        BisPileup=kTRUE;
      }

      if (fAnalysisUtil->IsPileUpEvent(InputEvent())) {
        fPileUpCount->Fill(7.5);
        BisPileup=kTRUE;
      }

    }
    else {

      // pileup for LHC15o, using AliMultSelection

      if(fMultSelection) {
        centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
        centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
      } else {
        BisPileup=kTRUE;
      }

      // pileup from AliMultSelection
      if(!fMultSelection->GetThisEventIsNotPileup()) fPileUpMultSelCount->Fill(0.5);
      if(!fMultSelection->GetThisEventIsNotPileupMV()) fPileUpMultSelCount->Fill(1.5);
      if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) fPileUpMultSelCount->Fill(2.5);
      if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) fPileUpMultSelCount->Fill(3.5);
      if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) fPileUpMultSelCount->Fill(4.5);
      if(!fMultSelection->GetThisEventIsNotAsymmetricInVZERO()) fPileUpMultSelCount->Fill(5.5);
      if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) fPileUpMultSelCount->Fill(6.5);
      if(!fMultSelection->GetThisEventHasGoodVertex2016()) fPileUpMultSelCount->Fill(7.5);

      // pile-up a la Dobrin for LHC15o
      if (plpMV(aod)) {
        fPileUpCount->Fill(0.5);
        BisPileup=kTRUE;
      }

      Short_t isPileup = aod->IsPileupFromSPD(3);
      if (isPileup != 0) {
        fPileUpCount->Fill(1.5);
        BisPileup=kTRUE;
      }

      if (((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
        fPileUpCount->Fill(2.5);
        BisPileup=kTRUE;
      }

      if (aod->IsIncompleteDAQ())  {
        fPileUpCount->Fill(3.5);
        BisPileup=kTRUE;
      }

      if(fabs(centrV0M-centrCL1)>7.5)  {
        fPileUpCount->Fill(4.5);
        BisPileup=kTRUE;
      }

      // check vertex consistency
      const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
      const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();

      if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1)  {
        fPileUpCount->Fill(5.5);
        BisPileup=kTRUE;
      }

      double covTrc[6], covSPD[6];
      vtTrc->GetCovarianceMatrix(covTrc);
      vtSPD->GetCovarianceMatrix(covSPD);

      double dz = vtTrc->GetZ() - vtSPD->GetZ();

      double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
      double errTrc = TMath::Sqrt(covTrc[5]);
      double nsigTot = dz/errTot;
      double nsigTrc = dz/errTrc;

      if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
        fPileUpCount->Fill(6.5);
        BisPileup=kTRUE;
      }

      // cuts on tracks
      const Int_t nTracks = aod->GetNumberOfTracks();
      Int_t multEsd = ((AliAODHeader*)aod->GetHeader())->GetNumberOfESDTracks();

      Int_t multTrk = 0;
      Int_t multTrkBefC = 0;
      Int_t multTrkTOFBefC = 0;
      Int_t multTPC = 0;

      for (Int_t it = 0; it < nTracks; it++) {

        AliAODTrack* aodTrk = (AliAODTrack*)aod->GetTrack(it);
        if (!aodTrk){
          delete aodTrk;
          continue;
        }

        //      if (aodTrk->TestFilterBit(32)){
        //        multTrkBefC++;
        //
        //        if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
        //          multTrkTOFBefC++;
        //
        //        if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2) && (aodTrk->Pt() < 20.))
        //          multTrk++;
        //      }

        if (aodTrk->TestFilterBit(128))
        multTPC++;

        if (centrV0M<10. && aodTrk->TestFilterBit(768) && aodTrk->Pt()>0.2) {

          // cut on # TPC clusters
          Int_t ntpccls = aodTrk->GetTPCNcls();
          fQATrackTPCNcls->Fill(aodTrk->Phi(),aodTrk->Eta(),ntpccls);

          // cut on # ITS clusters
          Int_t nitscls = aodTrk->GetITSNcls();
          fQATrackITSNcls->Fill(aodTrk->Phi(),aodTrk->Eta(),nitscls);

          // cut on chi2 / # TPC clusters
          Double_t chi2tpc = 0.;
          if(ntpccls>0) {
            chi2tpc = aodTrk->Chi2perNDF();
          }
          fQATrackTPCchi2->Fill(aodTrk->Phi(),aodTrk->Eta(),chi2tpc);

          // cut on chi2 / # ITS clusters
          Double_t chi2its = 0.;
          if(nitscls>0) {
            chi2its = aodTrk->GetITSchi2()/aodTrk->GetITSNcls();
          }
          fQATrackITSchi2->Fill(aodTrk->Phi(),aodTrk->Eta(),chi2its);

          // cut on fraction shared TPC clusters
          Double_t fshtpccls = 0.;
          if(ntpccls>0) {
            Int_t ntpcclsS = aodTrk->GetTPCnclsS();
            fshtpccls = 1.*ntpcclsS/ntpccls;
          }
          fQATrackTPCScls->Fill(aodTrk->Phi(),aodTrk->Eta(),fshtpccls);

          // cut on fraction shared ITS clusters
          Double_t fshitscls = 0.;
          Int_t nshcl = 0;
          if(nitscls>0) {
            for (Int_t i=0; i<6; i++) {
              if(aodTrk->HasSharedPointOnITSLayer(i)) nshcl++;
            }
            fshitscls = 1.*nshcl/nitscls;
          }
          fQATrackITSScls->Fill(aodTrk->Phi(),aodTrk->Eta(),fshitscls);

        }

      } // end of for (Int_t it = 0; it < nTracks; it++)

      Double_t multTPCn = multTPC;
      Double_t multEsdn = multEsd;
      Double_t multESDTPCDif = multEsdn - multTPCn*3.38;

      if (multESDTPCDif > (fRejectPileUpTight?700.:15000.)) {
        fPileUpCount->Fill(7.5);
        BisPileup=kTRUE;
      }

      if(fRejectPileUpTight) {
        if(BisPileup==kFALSE) {
          if(!fMultSelection->GetThisEventIsNotPileup()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
          if(BisPileup) fPileUpCount->Fill(8.5);
        }
      }
    }
  }


  return BisPileup;
}

//________________________________________________________________________

Double_t AliAnalysisTaskCRCZDC::GetBadTowerResp(Double_t Et, TH2D* BadTowerCalibHist)
{
  Double_t EtC = BadTowerCalibHist->ProjectionY("",BadTowerCalibHist->GetXaxis()->FindBin(Et),BadTowerCalibHist->GetXaxis()->FindBin(Et))->GetRandom();
  return EtC;
}

//________________________________________________________________________

Int_t AliAnalysisTaskCRCZDC::GetCenBin(Double_t Centrality)
{
  Int_t CenBin=-1;
  if (Centrality>0. && Centrality<5.) CenBin=0;
  if (Centrality>5. && Centrality<10.) CenBin=1;
  if (Centrality>10. && Centrality<20.) CenBin=2;
  if (Centrality>20. && Centrality<30.) CenBin=3;
  if (Centrality>30. && Centrality<40.) CenBin=4;
  if (Centrality>40. && Centrality<50.) CenBin=5;
  if (Centrality>50. && Centrality<60.) CenBin=6;
  if (Centrality>60. && Centrality<70.) CenBin=7;
  if (Centrality>70. && Centrality<80.) CenBin=8;
  if (Centrality>80. && Centrality<90.) CenBin=9;
  if (CenBin>=fnCen) CenBin=-1;
  if (fnCen==1) CenBin=0;
  return CenBin;
} // end of AliFlowAnalysisCRC::GetCRCCenBin(Double_t Centrality)
//_____________________________________________________________________________

Double_t AliAnalysisTaskCRCZDC::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    printf("One of vertices is not valid\n");
    return 0;
  }
  static TMatrixDSym vVb(3);
  double dist = -1;
  double dx = v0->GetX()-v1->GetX();
  double dy = v0->GetY()-v1->GetY();
  double dz = v0->GetZ()-v1->GetZ();
  double cov0[6],cov1[6];
  v0->GetCovarianceMatrix(cov0);
  v1->GetCovarianceMatrix(cov1);
  vVb(0,0) = cov0[0]+cov1[0];
  vVb(1,1) = cov0[2]+cov1[2];
  vVb(2,2) = cov0[5]+cov1[5];
  vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
  vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
  vVb.InvertFast();
  if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
  +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;
}
//________________________________________________________________________

Bool_t AliAnalysisTaskCRCZDC::plpMV(const AliAODEvent* aod)
{
  // check for multi-vertexer pile-up

  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2 = 5.0;
  const double kMinWDist = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;
  int nPlp = 0;

  if ( !(nPlp=aod->GetNumberOfPileupVerticesTracks()) ) return kFALSE;
  vtPrm = aod->GetPrimaryVertex();
  if (vtPrm == aod->GetPrimaryVertexSPD()) return kTRUE; // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for (int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)aod->GetPileupVertexTracks(ipl);
    //
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF() > kMaxPlpChi2) continue;
    //  int bcPlp = vtPlp->GetBC();
    //  if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other BC
    //
    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist) continue;
    //
    return kTRUE; // pile-up: well separated vertices
  }

  return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::SetCutsRP(AliFlowTrackCuts* cutsRP) {
  fCutContainer->AddAt(cutsRP,0); fCutsRP=cutsRP; cutsRP->SetPOItype(0);
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::SetCutsPOI(AliFlowTrackCuts* cutsPOI) {
  fCutContainer->AddAt(cutsPOI,1); fCutsPOI=cutsPOI; cutsPOI->SetPOItype(1);
}

//________________________________________________________________________
void AliAnalysisTaskCRCZDC::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  /*  if(fDebug > 1) printf(" **** AliAnalysisTaskCRCZDC::Terminate() \n");

  //fOutput = dynamic_cast<TList*> (GetOutputData(1));
  //if(!fOutput) printf("ERROR: fOutput not available\n");
  */
}

void AliAnalysisTaskCRCZDC::NotifyRun()
{
  //open file
  TGrid::Connect("alien://");
  if (fStepZDCRecenter >= 3) {
    TString ZDCRecenterFileName = Form("alien:///alice/cern.ch/user/s/sqiu/15o_ZDCRunByRunCalib/15o_ZDCcalibVar_%d.root",fCurrentRunNumber);
    TFile* ZDCRecenterFileRunByRun = TFile::Open(ZDCRecenterFileName, "READ");
    if(fStepZDCRecenter > 0) {
      if(ZDCRecenterFileRunByRun) {
        TList* ZDCRecenterListRunByRun = (TList*)(ZDCRecenterFileRunByRun->FindObjectAny("Q Vectors")); // hardcoded TList Q Vectors
        if(ZDCRecenterListRunByRun) {
	      SetZDCCalibListStep3RunByRun(ZDCRecenterListRunByRun);
	    } else {
          std::cout << "ERROR: ZDCRecenterList do not exist!" << std::endl;
          exit(1);
        }
      } else {
	    std::cout << "ERROR: if fStepZDCRecenter larger than 0, ZDCRecenterFile should exist!" << std::endl;
        exit(1);
      }
    }
    delete ZDCRecenterFileRunByRun;
  }

}
