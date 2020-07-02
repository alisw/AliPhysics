/****************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.   *
 *                                                                          *
 * Authors: Friederike Bock                                                 *
 * Version 1.0                                                              *
 *                                                                          *
 * Permission to use, copy, modify and distribute this software and its     *
 * documentation strictly for non-commercial purposes is hereby granted     *
 * without fee, provided that the above copyright notice appears in all     *
 * copies and that both the copyright notice and this permission notice     *
 * appear in the supporting documentation. The authors make no claims       *
 * about the suitability of this software for any purpose. It is            *
 * provided "as is" without express or implied warranty.                    *
 ***************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class handling photon selection cuts for
// Gamma Conversion analysis
//---------------------------------------------
////////////////////////////////////////////////

#include "AliConversionPhotonCuts.h"

#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliDataFile.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "AliMCEvent.h"
#include "AliAODConversionMother.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliOADBContainer.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliTriggerAnalysis.h"
#include "AliV0ReaderV1.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliTRDTriggerAnalysis.h"
#include "AliDalitzAODESDMC.h"
#include "AliDalitzEventMC.h"

class iostream;

using namespace std;

/// \cond CLASSIMP
ClassImp(AliConversionPhotonCuts)
/// \endcond

const char* AliConversionPhotonCuts::fgkCutNames[AliConversionPhotonCuts::kNCuts] = {
  "V0FinderType",           // 0
  "EtaCut",                 // 1
  "MinRCut",                // 2
  "EtaForPhiCut",           // 3
  "MinPhiCut",              // 4
  "MaxPhiCut",              // 5
  "SinglePtCut",            // 6
  "ClsTPCCut",              // 7
  "ededxSigmaCut",          // 8
  "pidedxSigmaCut",         // 9
  "piMomdedxSigmaCut",      // 10
  "piMaxMomdedxSigmaCut",   // 11
  "LowPRejectionSigmaCut",  // 12
  "TOFelectronPID",         // 13
  "ITSelectronPID",         // 14 -- new ITS PID
  "TRDelectronPID",         // 15 -- new TRD PID
  "QtMaxCut",               // 16
  "Chi2GammaCut",           // 17
  "PsiPair",                // 18
  "DoPhotonAsymmetryCut",   // 19
  "CosinePointingAngle",    // 20
  "SharedElectronCuts",     // 21
  "RejectToCloseV0s",       // 22
  "DcaRPrimVtx",            // 23
  "DcaZPrimVtx",            // 24
  "EventPlane"              // 25
};

//________________________________________________________________________
AliConversionPhotonCuts::AliConversionPhotonCuts(const char *name,const char *title) :
  AliAnalysisCuts(name,title),
  fHistograms(NULL),
  fPIDResponse(NULL),
  fDoLightOutput(0),
  fDoPlotTrackPID(kFALSE),
  fV0ReaderName("V0ReaderV1"),
  fMaxR(200),
  fMinR(0),
  fEtaCut(0.9),
  fEtaCutMin(-0.1),
  fEtaForPhiCutMin(-10.),
  fEtaForPhiCutMax(10.),
  fMinPhiCut(0.),
  fMaxPhiCut(100.),
  fDoShrinkTPCAcceptance(0),
  fPtCut(0.02),
  fSinglePtCut(0),
  fSinglePtCut2(0),
  fDoAsymPtCut(kFALSE),
  fMaxZ(1000),
  fMinClsTPC(0.),
  fMinClsTPCToF(0.),
  fMaxTPCChi2NDF(0.),
  fLineCutZRSlope(0.),
  fLineCutZValue(0),
  fLineCutZRSlopeMin(0.),
  fLineCutZValueMin(0),
  fChi2CutConversion(1000),
  fChi2CutConversionExpFunc(-1),
  fPIDProbabilityCutNegativeParticle(0),
  fPIDProbabilityCutPositiveParticle(0),
  fDodEdxSigmaCut(kTRUE),
  fDoTOFsigmaCut(kFALSE),
  fPIDTRDEfficiency(1),
  fDoTRDPID(kFALSE),
  fPIDnSigmaAboveElectronLine(100),
  fPIDnSigmaBelowElectronLine(-100),
  fTofPIDnSigmaAboveElectronLine(100),
  fTofPIDnSigmaBelowElectronLine(-100),
  fPIDnSigmaAbovePionLine(0),
  fPIDnSigmaAbovePionLineHighPt(-100),
  fPIDMinPnSigmaAbovePionLine(0),
  fPIDMaxPnSigmaAbovePionLine(0),
  fDoKaonRejectionLowP(kFALSE),
  fDoProtonRejectionLowP(kFALSE),
  fDoPionRejectionLowP(kFALSE),
  fPIDnSigmaAtLowPAroundKaonLine(0),
  fPIDnSigmaAtLowPAroundProtonLine(0),
  fPIDnSigmaAtLowPAroundPionLine(0),
  fPIDMinPKaonRejectionLowP(1.5),
  fPIDMinPProtonRejectionLowP(2),
  fPIDMinPPionRejectionLowP(0),
  fDoQtGammaSelection(1),
  fDo2DQt(kFALSE),
  fQtMax(100),
  fQtPtMax(100),
  fNSigmaMass(0.),
  fUseEtaMinCut(kFALSE),
  fUseOnFlyV0Finder(kTRUE),
  fUseOnFlyV0FinderSameSign(0),
  fUseBDTPhotonCuts(0),
  fDoPhotonAsymmetryCut(kTRUE),
  fDoPhotonPDependentAsymCut(kFALSE),
  fFAsymmetryCut(0),
  fMinPPhotonAsymmetryCut(100.),
  fMinPhotonAsymmetry(0.),
  fMaxPhotonAsymmetry(0.95),
  fUseCorrectedTPCClsInfo(kFALSE),
  fUseTOFpid(kFALSE),
  fUseTOFtiming(kFALSE),
  fTOFtimeMin(-1000),
  fTOFtimeMax(1000),
  fTOFtimingBothLegs(kFALSE),
  fUseTOFpidMinMom(kFALSE),
  fTofPIDMinMom(0.4),
  fOpeningAngle(0.005),
  fPsiPairCut(10000),
  fDo2DPsiPairChi2(0),
  fIncludeRejectedPsiPair(kFALSE),
  fCosPAngleCut(10000),
  fDoToCloseV0sCut(kFALSE),
  fminV0Dist(200.),
  fDoSharedElecCut(kFALSE),
  fDoPhotonQualitySelectionCut(kFALSE),
  fDoPhotonQualityRejectionCut(kFALSE),
  fPhotonQualityCut(0),
  fPhotonQualityCutTRD(0),
  fPhotonQualityCutTOF(0),
  fRandom(0),
  fElectronArraySize(500),
  fElectronLabelArray(NULL),
  fDCAZPrimVtxCut(1000),
  fDCARPrimVtxCut(1000),
  fInPlaneOutOfPlane(0),
  fConversionPointXArray(0.0),
  fConversionPointYArray(0.0),
  fConversionPointZArray(0.0),
  fCutString(NULL),
  fCutStringRead(""),
  fIsHeavyIon(0),
  fUseITSpid(kFALSE),
  fITSPIDnSigmaAboveElectronLine(100),
  fITSPIDnSigmaBelowElectronLine(-100),
  fMaxPtPIDITS(1.5),
  fTRDPIDAboveCut(100),
  fTRDPIDBelowCut(-100),
  fDoDoubleCountingCut(kFALSE),
  fMinRDC(0.),
  fDeltaR(0.),
  fOpenAngle(0.),
  fSwitchToKappa(kFALSE),
  fKappaMinCut(-1),
  fKappaMaxCut(1000),
  fDoElecDeDxPostCalibration(kFALSE),
  fIsRecalibDepTPCCl(kTRUE),
  fHistoEtaDistV0s(NULL),
  fHistoEtaDistV0sAfterdEdxCuts(NULL),
  fHistodEdxCuts(NULL),
  fHistoTPCdEdxbefore(NULL),
  fHistoTPCdEdxafter(NULL),
  fHistoTPCdEdxSigbefore(NULL),
  fHistoTPCdEdxSigafter(NULL),
  fHistoTPCChi2NDFBefore(NULL),
  fHistoTPCChi2NDFAfter(NULL),
  fHistoTPCChi2NDF2D(NULL),
  fHistoKappaafter(NULL),
  fHistoTOFbefore(NULL),
  fHistoTOFSigbefore(NULL),
  fHistoTOFSigafter(NULL),
  fHistoITSSigbefore(NULL),
  fHistoITSSigafter(NULL),
  fHistoPsiPairDeltaPhiafter(NULL),
  fHistoTrackCuts(NULL),
  fHistoTrackPID(NULL),
  fHistoPhotonCuts(NULL),
  fHistoInvMassbefore(NULL),
  fHistoArmenterosbefore(NULL),
  fHistoInvMassafter(NULL),
  fHistoArmenterosafter(NULL),
  fHistoAsymmetrybefore(NULL),
  fHistoAsymmetryafter(NULL),
  fHistoAcceptanceCuts(NULL),
  fHistoCutIndex(NULL),
  fHistoTOFtimeVSMomentum(NULL),
  fHistoEventPlanePhi(NULL),
  fPreSelCut(kFALSE),
  fProcessAODCheck(kFALSE),
  fMaterialBudgetWeightsInitialized(kFALSE),
  fProfileContainingMaterialBudgetWeights(NULL),
  fFileNameElecDeDxPostCalibration(""),
  fElecDeDxPostCalibrationInitialized(kFALSE),
  fRecalibCurrentRun(-1),
  fnRBins(4),
  fHistoEleMapRecalib(NULL),
  fHistoPosMapRecalib(NULL),
  fGoodRegionCMin(0),
  fGoodRegionAMin(0),
  fBadRegionCMin(0),
  fBadRegionAMin(0),
  fGoodRegionCMax(0),
  fGoodRegionAMax(0),
  fBadRegionCMax(0),
  fBadRegionAMax(0),
  fExcludeMinR(180.),
  fExcludeMaxR(250.)
{
  InitPIDResponse();
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=0;}
  fCutString=new TObjString((GetCutNumber()).Data());

  fElectronLabelArray = new Int_t[fElectronArraySize];

  fHistoEleMapRecalib  = new TH2S*[fnRBins];
  fHistoPosMapRecalib  = new TH2S*[fnRBins];

  for (Int_t i = 0; i < fnRBins; i++) {
    fHistoEleMapRecalib[i] = NULL;
    fHistoPosMapRecalib[i] = NULL;
  }
}

//________________________________________________________________________
AliConversionPhotonCuts::AliConversionPhotonCuts(const AliConversionPhotonCuts &ref) :
  AliAnalysisCuts(ref),
  fHistograms(NULL),
  fPIDResponse(NULL),
  fDoLightOutput(ref.fDoLightOutput),
  fDoPlotTrackPID(ref.fDoPlotTrackPID),
  fV0ReaderName("V0ReaderV1"),
  fMaxR(ref.fMaxR),
  fMinR(ref.fMinR),
  fEtaCut(ref.fEtaCut),
  fEtaCutMin(ref.fEtaCutMin),
  fEtaForPhiCutMin(ref.fEtaForPhiCutMin),
  fEtaForPhiCutMax(ref.fEtaForPhiCutMax),
  fMinPhiCut(ref.fMinPhiCut),
  fMaxPhiCut(ref.fMaxPhiCut),
  fDoShrinkTPCAcceptance(ref.fDoShrinkTPCAcceptance),
  fPtCut(ref.fPtCut),
  fSinglePtCut(ref.fSinglePtCut),
  fSinglePtCut2(ref.fSinglePtCut2),
  fDoAsymPtCut(ref.fDoAsymPtCut),
  fMaxZ(ref.fMaxZ),
  fMinClsTPC(ref.fMinClsTPC),
  fMinClsTPCToF(ref.fMinClsTPCToF),
  fMaxTPCChi2NDF(ref.fMaxTPCChi2NDF),
  fLineCutZRSlope(ref.fLineCutZRSlope),
  fLineCutZValue(ref.fLineCutZValue),
  fLineCutZRSlopeMin(ref.fLineCutZRSlopeMin),
  fLineCutZValueMin(ref.fLineCutZValueMin),
  fChi2CutConversion(ref.fChi2CutConversion),
  fChi2CutConversionExpFunc(ref.fChi2CutConversionExpFunc),
  fPIDProbabilityCutNegativeParticle(ref.fPIDProbabilityCutNegativeParticle),
  fPIDProbabilityCutPositiveParticle(ref.fPIDProbabilityCutPositiveParticle),
  fDodEdxSigmaCut(ref. fDodEdxSigmaCut),
  fDoTOFsigmaCut(ref.fDoTOFsigmaCut),
  fPIDTRDEfficiency(ref.fPIDTRDEfficiency),
  fDoTRDPID(ref.fDoTRDPID),
  fPIDnSigmaAboveElectronLine(ref.fPIDnSigmaAboveElectronLine),
  fPIDnSigmaBelowElectronLine(ref.fPIDnSigmaBelowElectronLine),
  fTofPIDnSigmaAboveElectronLine(ref.fTofPIDnSigmaAboveElectronLine),
  fTofPIDnSigmaBelowElectronLine(ref.fTofPIDnSigmaBelowElectronLine),
  fPIDnSigmaAbovePionLine(ref.fPIDnSigmaAbovePionLine),
  fPIDnSigmaAbovePionLineHighPt(ref.fPIDnSigmaAbovePionLineHighPt),
  fPIDMinPnSigmaAbovePionLine(ref.fPIDMinPnSigmaAbovePionLine),
  fPIDMaxPnSigmaAbovePionLine(ref.fPIDMaxPnSigmaAbovePionLine),
  fDoKaonRejectionLowP(ref.fDoKaonRejectionLowP),
  fDoProtonRejectionLowP(ref.fDoProtonRejectionLowP),
  fDoPionRejectionLowP(ref.fDoPionRejectionLowP),
  fPIDnSigmaAtLowPAroundKaonLine(ref.fPIDnSigmaAtLowPAroundKaonLine),
  fPIDnSigmaAtLowPAroundProtonLine(ref.fPIDnSigmaAtLowPAroundProtonLine),
  fPIDnSigmaAtLowPAroundPionLine(ref.fPIDnSigmaAtLowPAroundPionLine),
  fPIDMinPKaonRejectionLowP(ref.fPIDMinPKaonRejectionLowP),
  fPIDMinPProtonRejectionLowP(ref.fPIDMinPProtonRejectionLowP),
  fPIDMinPPionRejectionLowP(ref.fPIDMinPPionRejectionLowP),
  fDoQtGammaSelection(ref.fDoQtGammaSelection),
  fDo2DQt(ref.fDo2DQt),
  fQtMax(ref.fQtMax),
  fQtPtMax(ref.fQtPtMax),
  fNSigmaMass(ref.fNSigmaMass),
  fUseEtaMinCut(ref.fUseEtaMinCut),
  fUseOnFlyV0Finder(ref.fUseOnFlyV0Finder),
  fUseOnFlyV0FinderSameSign(ref.fUseOnFlyV0FinderSameSign),
  fUseBDTPhotonCuts(ref.fUseBDTPhotonCuts),
  fDoPhotonAsymmetryCut(ref.fDoPhotonAsymmetryCut),
  fDoPhotonPDependentAsymCut(ref.fDoPhotonPDependentAsymCut),
  fFAsymmetryCut(ref.fFAsymmetryCut),
  fMinPPhotonAsymmetryCut(ref.fMinPPhotonAsymmetryCut),
  fMinPhotonAsymmetry(ref.fMinPhotonAsymmetry),
  fMaxPhotonAsymmetry(ref.fMaxPhotonAsymmetry),
  fUseCorrectedTPCClsInfo(ref.fUseCorrectedTPCClsInfo),
  fUseTOFpid(ref.fUseTOFpid),
  fUseTOFtiming(ref.fUseTOFtiming),
  fTOFtimeMin(ref.fTOFtimeMin),
  fTOFtimeMax(ref.fTOFtimeMax),
  fTOFtimingBothLegs(ref.fTOFtimingBothLegs),
  fUseTOFpidMinMom(ref.fUseTOFpidMinMom),
  fTofPIDMinMom(ref.fTofPIDMinMom),
  fOpeningAngle(ref.fOpeningAngle),
  fPsiPairCut(ref.fPsiPairCut),
  fDo2DPsiPairChi2(ref.fDo2DPsiPairChi2),
  fIncludeRejectedPsiPair(ref.fIncludeRejectedPsiPair),
  fCosPAngleCut(ref.fCosPAngleCut),
  fDoToCloseV0sCut(ref.fDoToCloseV0sCut),
  fminV0Dist(ref.fminV0Dist),
  fDoSharedElecCut(ref.fDoSharedElecCut),
  fDoPhotonQualitySelectionCut(ref.fDoPhotonQualitySelectionCut),
  fDoPhotonQualityRejectionCut(ref.fDoPhotonQualityRejectionCut),
  fPhotonQualityCut(ref.fPhotonQualityCut),
  fPhotonQualityCutTRD(ref.fPhotonQualityCutTRD),
  fPhotonQualityCutTOF(ref.fPhotonQualityCutTOF),
  fRandom(ref.fRandom),
  fElectronArraySize(ref.fElectronArraySize),
  fElectronLabelArray(NULL),
  fDCAZPrimVtxCut(ref.fDCAZPrimVtxCut),
  fDCARPrimVtxCut(ref.fDCAZPrimVtxCut),
  fInPlaneOutOfPlane(ref.fInPlaneOutOfPlane),
  fConversionPointXArray(ref.fConversionPointXArray),
  fConversionPointYArray(ref.fConversionPointYArray),
  fConversionPointZArray(ref.fConversionPointZArray),
  fCutString(NULL),
  fCutStringRead(""),
  fIsHeavyIon(ref.fIsHeavyIon),
  fUseITSpid(ref.fUseITSpid),
  fITSPIDnSigmaAboveElectronLine(ref.fITSPIDnSigmaAboveElectronLine),
  fITSPIDnSigmaBelowElectronLine(ref.fITSPIDnSigmaBelowElectronLine),
  fMaxPtPIDITS(ref.fMaxPtPIDITS),
  fTRDPIDAboveCut(ref.fTRDPIDAboveCut),
  fTRDPIDBelowCut(ref.fTRDPIDBelowCut),
  fDoDoubleCountingCut(ref.fDoDoubleCountingCut),
  fMinRDC(ref.fMinRDC),
  fDeltaR(ref.fDeltaR),
  fOpenAngle(ref.fOpenAngle),
  fSwitchToKappa(ref.fSwitchToKappa),
  fKappaMinCut(ref.fKappaMinCut),
  fKappaMaxCut(ref.fKappaMaxCut),
  fDoElecDeDxPostCalibration(ref.fDoElecDeDxPostCalibration),
  fIsRecalibDepTPCCl(ref.fIsRecalibDepTPCCl),
  fHistoEtaDistV0s(NULL),
  fHistoEtaDistV0sAfterdEdxCuts(NULL),
  fHistodEdxCuts(NULL),
  fHistoTPCdEdxbefore(NULL),
  fHistoTPCdEdxafter(NULL),
  fHistoTPCdEdxSigbefore(NULL),
  fHistoTPCdEdxSigafter(NULL),
  fHistoTPCChi2NDFBefore(NULL),
  fHistoTPCChi2NDFAfter(NULL),
  fHistoTPCChi2NDF2D(NULL),
  fHistoKappaafter(NULL),
  fHistoTOFbefore(NULL),
  fHistoTOFSigbefore(NULL),
  fHistoTOFSigafter(NULL),
  fHistoITSSigbefore(NULL),
  fHistoITSSigafter(NULL),
  fHistoPsiPairDeltaPhiafter(NULL),
  fHistoTrackCuts(NULL),
  fHistoTrackPID(NULL),
  fHistoPhotonCuts(NULL),
  fHistoInvMassbefore(NULL),
  fHistoArmenterosbefore(NULL),
  fHistoInvMassafter(NULL),
  fHistoArmenterosafter(NULL),
  fHistoAsymmetrybefore(NULL),
  fHistoAsymmetryafter(NULL),
  fHistoAcceptanceCuts(NULL),
  fHistoCutIndex(NULL),
  fHistoTOFtimeVSMomentum(NULL),
  fHistoEventPlanePhi(NULL),
  fPreSelCut(ref.fPreSelCut),
  fProcessAODCheck(ref.fProcessAODCheck),
  fMaterialBudgetWeightsInitialized(ref.fMaterialBudgetWeightsInitialized),
  fProfileContainingMaterialBudgetWeights(ref.fProfileContainingMaterialBudgetWeights),
  fFileNameElecDeDxPostCalibration(ref.fFileNameElecDeDxPostCalibration),
  fElecDeDxPostCalibrationInitialized(ref.fElecDeDxPostCalibrationInitialized),
  fRecalibCurrentRun(ref.fRecalibCurrentRun),
  fnRBins(ref.fnRBins),
  fHistoEleMapRecalib(ref.fHistoEleMapRecalib),
  fHistoPosMapRecalib(ref.fHistoPosMapRecalib),
  fGoodRegionCMin(ref.fGoodRegionCMin),
  fGoodRegionAMin(ref.fGoodRegionAMin),
  fBadRegionCMin(ref.fBadRegionCMin),
  fBadRegionAMin(ref.fBadRegionAMin),
  fGoodRegionCMax(ref.fGoodRegionCMax),
  fGoodRegionAMax(ref.fGoodRegionAMax),
  fBadRegionCMax(ref.fBadRegionCMax),
  fBadRegionAMax(ref.fBadRegionAMax),
  fExcludeMinR(ref.fExcludeMinR),
  fExcludeMaxR(ref.fExcludeMaxR)
{
  // Copy Constructor
  for(Int_t jj=0;jj<kNCuts;jj++){fCuts[jj]=ref.fCuts[jj];}
  fCutString=new TObjString((GetCutNumber()).Data());
  fElectronLabelArray = new Int_t[fElectronArraySize];
  // dont copy histograms (if you like histograms, call InitCutHistograms())
}


//________________________________________________________________________
AliConversionPhotonCuts::~AliConversionPhotonCuts() {
  // Destructor
  //Deleting fHistograms leads to seg fault it it's added to output collection of a task
  // if(fHistograms)
  //    delete fHistograms;
  // fHistograms = NULL;
  if(fCutString != NULL){
    delete fCutString;
    fCutString = NULL;
  }
  if(fElectronLabelArray){
    delete fElectronLabelArray;
    fElectronLabelArray = NULL;
  }

  if(fFAsymmetryCut != NULL){
    delete fFAsymmetryCut;
    fFAsymmetryCut = NULL;
  }
  if(fProfileContainingMaterialBudgetWeights){
      delete fProfileContainingMaterialBudgetWeights;
      fProfileContainingMaterialBudgetWeights = 0x0;
  }

  // if( fHistoEleMapRecalib != NULL){
  //   delete fHistoEleMapRecalib;
  //   fHistoEleMapRecalib =NULL;
  // }

  for (Int_t i = 0; i < fnRBins; i++) {
    if( fHistoEleMapRecalib[i]  != NULL){
      delete fHistoEleMapRecalib[i] ;
      fHistoEleMapRecalib[i]  =NULL;
    }

    if( fHistoPosMapRecalib[i]  != NULL){
      delete fHistoPosMapRecalib[i] ;
      fHistoPosMapRecalib[i]  =NULL;
    }

  }


}

//________________________________________________________________________
void AliConversionPhotonCuts::InitCutHistograms(TString name, Bool_t preCut){

  // Initialize Cut Histograms for QA (only initialized and filled if function is called)
  TH1::AddDirectory(kFALSE);


  if(fHistograms != NULL){
    delete fHistograms;
    fHistograms=NULL;
  }

  if(fDoLightOutput==2) {
      AliInfo("Minimal output chosen");
      return;
  }

  if(fHistograms==NULL){
    fHistograms=new TList();
    fHistograms->SetOwner(kTRUE);
    if(name=="")fHistograms->SetName(Form("ConvCuts_%s",GetCutNumber().Data()));
    else fHistograms->SetName(Form("%s_%s",name.Data(),GetCutNumber().Data()));
  }

  // IsPhotonSelected
  fHistoCutIndex=new TH1F(Form("IsPhotonSelected %s",GetCutNumber().Data()),"IsPhotonSelected",12,-0.5,11.5);
  fHistoCutIndex->GetXaxis()->SetBinLabel(kPhotonIn+1,"in");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kOnFly+1,"onfly");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kNoTracks+1,"no tracks");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kNoV0+1,"miss. V0 in AOD");
  if (!fSwitchToKappa)fHistoCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"PID");
  else fHistoCutIndex->GetXaxis()->SetBinLabel(kdEdxCuts+1,"Kappa+[TOF,ITS,TRD] PID");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kTrackCuts+1,"Track cuts");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kConvPointFail+1,"ConvPoint fail");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kPhotonCuts+1,"PhotonCuts");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kEventPlane+1,"EventPlane");
  fHistoCutIndex->GetXaxis()->SetBinLabel(kPhotonOut+1,"out");
  fHistograms->Add(fHistoCutIndex);

  // Track Cuts
  fHistoTrackCuts=new TH1F(Form("TrackCuts %s",GetCutNumber().Data()),"TrackCuts",10,-0.5,9.5);
  fHistoTrackCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(2,"likesign");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(3,"ntpccl");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(4,"acceptance");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(5,"singlept");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(6,"TPCrefit");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(7,"kink");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(8,"TPCChi2");
  fHistoTrackCuts->GetXaxis()->SetBinLabel(9,"out");
  fHistograms->Add(fHistoTrackCuts);

  // Photon Cuts
  fHistoPhotonCuts=new TH2F(Form("PhotonCuts %s",GetCutNumber().Data()),"PhotonCuts vs p_{T,#gamma}",15,-0.5,14.5,250,0,50);
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(2,"qtcut");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(3,"chi2");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(4,"acceptance");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(5,"asymmetry");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(6,"pidprob");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(7,"cortpcclinfo");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(8,"PsiPair");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(9,"CosPAngle");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(10,"DCA R");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(11,"DCA Z");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(12,"Photon Quality");
  fHistoPhotonCuts->GetXaxis()->SetBinLabel(13,"out");
  fHistograms->Add(fHistoPhotonCuts);

  if (fProfileContainingMaterialBudgetWeights){
      fProfileContainingMaterialBudgetWeights->SetName("InputMaterialBudgetWeightsPerGamma");
      fHistograms->Add(fProfileContainingMaterialBudgetWeights);
  }

  if (fDoElecDeDxPostCalibration || fElecDeDxPostCalibrationInitialized){
    for (Int_t i = 0; i < fnRBins; i++) {
      if( fHistoEleMapRecalib[i] ){
        fHistograms->Add(fHistoEleMapRecalib[i]);
      }
      if(fHistoPosMapRecalib[i]){
        fHistograms->Add(fHistoPosMapRecalib[i]);
      }
    }
  }

  if(fUseTOFtiming){
    fHistoTOFtimeVSMomentum=new TH2F(Form("TOFtime_Momentum %s",GetCutNumber().Data()),"TOFtime_Momentum",400,-150,250,400,0,20.);
    fHistograms->Add(fHistoTOFtimeVSMomentum);
  }

  if(!fDoLightOutput){

    if(preCut){
      fHistoInvMassbefore=new TH1F(Form("InvMass_before %s",GetCutNumber().Data()),"InvMass_before",1000,0,0.3);
      fHistograms->Add(fHistoInvMassbefore);
      fHistoArmenterosbefore=new TH2F(Form("Armenteros_before %s",GetCutNumber().Data()),"Armenteros_before",200,-1,1,1000,0,1.);
      fHistograms->Add(fHistoArmenterosbefore);
      fHistoEtaDistV0s = new TH1F(Form("Eta_before %s",GetCutNumber().Data()),"Eta_before",2000,-2,2);
      fHistograms->Add(fHistoEtaDistV0s);
      fHistoAsymmetrybefore=new TH2F(Form("Asymmetry_before %s",GetCutNumber().Data()),"Asymmetry_before",150,0.03,20.,200,0,1.);
      fHistograms->Add(fHistoAsymmetrybefore);
    }
    fHistoInvMassafter=new TH1F(Form("InvMass_after %s",GetCutNumber().Data()),"InvMass_after",1000,0,0.3);
    fHistograms->Add(fHistoInvMassafter);
    fHistoArmenterosafter=new TH2F(Form("Armenteros_after %s",GetCutNumber().Data()),"Armenteros_after",200,-1,1,250,0,0.25);
    fHistograms->Add(fHistoArmenterosafter);
    // AM - save always to see distribution after selections cut
    //    if(fDoPhotonAsymmetryCut){
    fHistoAsymmetryafter=new TH2F(Form("Asymmetry_after %s",GetCutNumber().Data()),"Asymmetry_after",150,0.03,20.,200,0,1.);
    fHistograms->Add(fHistoAsymmetryafter);
      //    }
  }

  if (fDoPlotTrackPID) {
      fHistoTrackPID=new TH2F(Form("TrackPID %s",GetCutNumber().Data()),"TrackPID",15,-0.5,14.5,250,0,25);
      fHistoTrackPID->GetXaxis()->SetBinLabel(1,"Electron");
      fHistoTrackPID->GetXaxis()->SetBinLabel(2,"Muon");
      fHistoTrackPID->GetXaxis()->SetBinLabel(3,"Pion");
      fHistoTrackPID->GetXaxis()->SetBinLabel(4,"Kaon");
      fHistoTrackPID->GetXaxis()->SetBinLabel(5,"Proton");
      fHistoTrackPID->GetXaxis()->SetBinLabel(6,"Deuteron");
      fHistoTrackPID->GetXaxis()->SetBinLabel(7,"Triton");
      fHistoTrackPID->GetXaxis()->SetBinLabel(8,"He3");
      fHistoTrackPID->GetXaxis()->SetBinLabel(9,"Alpha");
      fHistoTrackPID->GetXaxis()->SetBinLabel(10,"Photon");
      fHistoTrackPID->GetXaxis()->SetBinLabel(11,"Pi0");
      fHistoTrackPID->GetXaxis()->SetBinLabel(12,"Neutron");
      fHistoTrackPID->GetXaxis()->SetBinLabel(13,"Kaon0");
      fHistoTrackPID->GetXaxis()->SetBinLabel(14,"EleCon");
      fHistoTrackPID->GetXaxis()->SetBinLabel(15,"Unknown");
      fHistograms->Add(fHistoTrackPID);
  }

  fHistoAcceptanceCuts=new TH2F(Form("PhotonAcceptanceCuts %s",GetCutNumber().Data()),"PhotonAcceptanceCuts vs p_{T,#gamma}",12,-0.5,11.5,250,0,50);
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(2,"maxR");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(3,"minR");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(4,"ExcludeR");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(5,"line");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(6,"maxZ");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(7,"eta");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(8,"phisector");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(9,"minpt");
  fHistoAcceptanceCuts->GetXaxis()->SetBinLabel(10,"out");
  fHistograms->Add(fHistoAcceptanceCuts);

  // dEdx Cuts
  fHistodEdxCuts=new TH2F(Form("dEdxCuts %s",GetCutNumber().Data()),"dEdxCuts vs p_{T,e}",11,-0.5,10.5,250,0,50);
  fHistodEdxCuts->GetXaxis()->SetBinLabel(1,"in");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(2,"TPCelectron");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(3,"TPCpion");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(4,"TPCpionhighp");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(5,"TPCkaonlowprej");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(6,"TPCprotonlowprej");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(7,"TPCpionlowprej");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(8,"TOFelectron");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(9,"ITSelectron");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(10,"TRDelectron");
  fHistodEdxCuts->GetXaxis()->SetBinLabel(11,"out");
  fHistograms->Add(fHistodEdxCuts);

  if(!fDoLightOutput){
    TAxis *AxisBeforedEdx = NULL;
    TAxis *AxisBeforedEdxSig = NULL;
    TAxis *AxisBeforeTOF = NULL;
    TAxis *AxisBeforeTOFSig = NULL;
    TAxis *AxisBeforeITSSig = NULL;
    TAxis *AxisBeforeAsymmetry = NULL;

    if(preCut){
      fHistoTPCdEdxbefore=new TH2F(Form("Gamma_dEdx_before %s",GetCutNumber().Data()),"dEdx Gamma before" ,150,0.03,20,800,0,200);
      fHistograms->Add(fHistoTPCdEdxbefore);
      AxisBeforedEdx = fHistoTPCdEdxbefore->GetXaxis();
      fHistoTPCdEdxSigbefore=new TH2F(Form("Gamma_dEdxSig_before %s",GetCutNumber().Data()),"dEdx Sigma Gamma before" ,150,0.03,20,400,-10,10);
      fHistograms->Add(fHistoTPCdEdxSigbefore);
      AxisBeforedEdxSig = fHistoTPCdEdxSigbefore->GetXaxis();

      fHistoTOFbefore=new TH2F(Form("Gamma_TOF_before %s",GetCutNumber().Data()),"TOF Gamma before" ,150,0.03,20,11000,-1000,10000);
      fHistograms->Add(fHistoTOFbefore);
      AxisBeforeTOF = fHistoTOFbefore->GetXaxis();
      fHistoTOFSigbefore=new TH2F(Form("Gamma_TOFSig_before %s",GetCutNumber().Data()),"TOF Sigma Gamma before" ,150,0.03,20,400,-6,10);
      fHistograms->Add(fHistoTOFSigbefore);
      AxisBeforeTOFSig = fHistoTOFSigbefore->GetXaxis();

      fHistoITSSigbefore=new TH2F(Form("Gamma_ITSSig_before %s",GetCutNumber().Data()),"ITS Sigma Gamma before" ,150,0.03,20,400,-10,10);
      fHistograms->Add(fHistoITSSigbefore);
      AxisBeforeITSSig = fHistoITSSigbefore->GetXaxis();

      AxisBeforeAsymmetry = fHistoAsymmetrybefore->GetXaxis();
    }

    fHistoTPCdEdxSigafter=new TH2F(Form("Gamma_dEdxSig_after %s",GetCutNumber().Data()),"dEdx Sigma Gamma after" ,150,0.03,20,400, -10,10);
    fHistograms->Add(fHistoTPCdEdxSigafter);

    fHistoTPCdEdxafter=new TH2F(Form("Gamma_dEdx_after %s",GetCutNumber().Data()),"dEdx Gamma after" ,150,0.03,20,800,0,200);
    fHistograms->Add(fHistoTPCdEdxafter);

    fHistoKappaafter=new TH2F(Form("Gamma_Kappa_after %s",GetCutNumber().Data()),"Kappa Gamma after" ,150,0.03,20,200,-20,20);
    fHistograms->Add(fHistoKappaafter);

    fHistoTOFSigafter=new TH2F(Form("Gamma_TOFSig_after %s",GetCutNumber().Data()),"TOF Sigma Gamma after" ,150,0.03,20,400,-6,10);
    fHistograms->Add(fHistoTOFSigafter);

    fHistoITSSigafter=new TH2F(Form("Gamma_ITSSig_after %s",GetCutNumber().Data()),"ITS Sigma Gamma after" ,150,0.03,20,400,-10,10);
    fHistograms->Add(fHistoITSSigafter);

    fHistoEtaDistV0sAfterdEdxCuts = new TH1F(Form("Eta_afterdEdx %s",GetCutNumber().Data()),"Eta_afterdEdx",2000,-2,2);
    fHistograms->Add(fHistoEtaDistV0sAfterdEdxCuts);

    fHistoPsiPairDeltaPhiafter=new TH2F(Form("Gamma_PsiPairDeltaPhi_after %s",GetCutNumber().Data()),"Psi Pair vs Delta Phi Gamma after" ,200,-2,2,200,-2,2);
    fHistograms->Add(fHistoPsiPairDeltaPhiafter);

    fHistoTPCChi2NDFBefore = new TH1F(Form("TPCChi2NDF_before %s",GetCutNumber().Data()),"TPCChi2NDF before cut",120,-2,10);
    fHistograms->Add(fHistoTPCChi2NDFBefore);
    fHistoTPCChi2NDFAfter = new TH1F(Form("TPCChi2NDF_after %s",GetCutNumber().Data()),"TPCChi2NDF after cut",120,-2,10);
    fHistograms->Add(fHistoTPCChi2NDFAfter);
    fHistoTPCChi2NDF2D = new TH2F(Form("TPCChi2NDF2D_before %s",GetCutNumber().Data()),"TPCChi2NDF neg vs pos track before cut",120,-2,10,120,-2,10);
    fHistograms->Add(fHistoTPCChi2NDF2D);

    TAxis *AxisAfter = fHistoTPCdEdxSigafter->GetXaxis();
    Int_t bins = AxisAfter->GetNbins();
    Double_t from = AxisAfter->GetXmin();
    Double_t to = AxisAfter->GetXmax();
    Double_t *newBins = new Double_t[bins+1];
    newBins[0] = from;
    Double_t factor = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];
    AxisAfter->Set(bins, newBins);
    AxisAfter = fHistoTOFSigafter->GetXaxis();
    AxisAfter->Set(bins, newBins);
    AxisAfter = fHistoTPCdEdxafter->GetXaxis();
    AxisAfter->Set(bins, newBins);
    AxisAfter = fHistoKappaafter->GetXaxis();
    AxisAfter->Set(bins, newBins);
    AxisAfter = fHistoITSSigafter->GetXaxis();
    AxisAfter->Set(bins, newBins);
    //    if(fDoPhotonAsymmetryCut){
    AxisAfter = fHistoAsymmetryafter->GetXaxis();
    AxisAfter->Set(bins, newBins);
      //    }
    if(preCut){
      AxisBeforedEdx->Set(bins, newBins);
      AxisBeforeTOF->Set(bins, newBins);
      AxisBeforedEdxSig->Set(bins, newBins);
      AxisBeforeTOFSig->Set(bins, newBins);
      AxisBeforeITSSig->Set(bins, newBins);
      AxisBeforeAsymmetry->Set(bins, newBins);
    }
    delete [] newBins;

    // Event Cuts and Info
    if(!preCut){
      fHistoEventPlanePhi=new TH1F(Form("EventPlaneMinusPhotonAngle %s",GetCutNumber().Data()),"EventPlaneMinusPhotonAngle",360,-TMath::Pi(),TMath::Pi());
      fHistograms->Add(fHistoEventPlanePhi);
    }
  }

  TH1::AddDirectory(kTRUE);
}

//________________________________________________________________________
Bool_t AliConversionPhotonCuts::InitPIDResponse(){
  // Set Pointer to AliPIDResponse

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    if(fPIDResponse)return kTRUE;

  }


  return kFALSE;
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::InitializeElecDeDxPostCalibration(TString filename) {
  AliInfo("Entering loading of correction map for post calibration");

  TFile *file = TFile::Open(filename.Data());
  if(!file){
    AliError(Form("file for electron dEdx post calibration %s not found",filename.Data()));
    return kFALSE;
  }else{
    AliInfo(Form("found %s ",filename.Data()));
  }

  for(Int_t i=0;i<fnRBins;i++){
   if (fIsRecalibDepTPCCl){
    fHistoEleMapRecalib[i]  = (TH2S*)file->Get(Form("Ele_Cl%d_recalib",i));
    fHistoPosMapRecalib[i]  = (TH2S*)file->Get(Form("Pos_Cl%d_recalib",i));
   }else{
    fHistoEleMapRecalib[i]  = (TH2S*)file->Get(Form("Ele_R%d_recalib",i));
    fHistoPosMapRecalib[i]  = (TH2S*)file->Get(Form("Pos_R%d_recalib",i));
   }
  }

  if (fHistoEleMapRecalib[0] == NULL || fHistoEleMapRecalib[1] == NULL ||
      fHistoEleMapRecalib[2] == NULL || fHistoEleMapRecalib[3] == NULL  ){
    AliFatal("Histograms for dedx post calibration not found in %s despite being requested!");
    return kFALSE;// code must break if histograms are not found!
  }
  for(Int_t i=0;i<fnRBins;i++){
    fHistoEleMapRecalib[i]  ->SetDirectory(0);
    fHistoPosMapRecalib[i]  ->SetDirectory(0);
  }

  file->Close();
  delete file;
  fElecDeDxPostCalibrationInitialized=kTRUE;
  return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::LoadElecDeDxPostCalibration(Int_t runNumber) {

  if(runNumber==fRecalibCurrentRun || fElecDeDxPostCalibrationInitialized)
    return kTRUE;
  else
    fRecalibCurrentRun=runNumber;

  // Else choose the one in the $ALICE_PHYSICS directory (or EOS)
  AliInfo("LoadingdEdx recalibration maps from OADB/PWGGA/TPCdEdxRecalibOADB.root");

  TFile *fileRecalib = TFile::Open(AliDataFile::GetFileNameOADB("PWGGA/TPCdEdxRecalibOADB.root").data(),"read");
  if (!fileRecalib || fileRecalib->IsZombie())
    {
      AliWarning("OADB/PWGGA/TPCdEdxRecalibOADB.root was not found");
      return kFALSE;
    }
  if (fileRecalib) delete fileRecalib;
  AliOADBContainer *contRecalibTPC = new AliOADBContainer("");
  contRecalibTPC->InitFromFile(AliDataFile::GetFileNameOADB("PWGGA/TPCdEdxRecalibOADB.root").data(),"AliTPCdEdxRecalib");

  TObjArray *arrayTPCRecalib=(TObjArray*)contRecalibTPC->GetObject(runNumber);
  if (!arrayTPCRecalib)
    {
      AliWarning(Form("No TPC dEdx recalibration found for run number: %d", runNumber));
      delete contRecalibTPC;
      return kFALSE;
    }

  for(Int_t i=0;i<fnRBins;i++){
    if (fIsRecalibDepTPCCl){
      fHistoEleMapRecalib[i]  = (TH2S*)arrayTPCRecalib->FindObject(Form("Ele_Cl%d_recalib",i));
      fHistoPosMapRecalib[i]  = (TH2S*)arrayTPCRecalib->FindObject(Form("Pos_Cl%d_recalib",i));
    } else {
      fHistoEleMapRecalib[i]  = (TH2S*)arrayTPCRecalib->FindObject(Form("Ele_R%d_recalib",i));
      fHistoPosMapRecalib[i]  = (TH2S*)arrayTPCRecalib->FindObject(Form("Pos_R%d_recalib",i));
    }
  }
  if (fHistoEleMapRecalib[0] == NULL || fHistoEleMapRecalib[1] == NULL ||
      fHistoEleMapRecalib[2] == NULL || fHistoEleMapRecalib[3] == NULL  ){
    AliWarning("Histograms for dedx post calibration not found in %s despite being requested!");
    return kFALSE;// code must break if histograms are not found!
  }
  for(Int_t i=0;i<fnRBins;i++){
    fHistoEleMapRecalib[i]  ->SetDirectory(0);
    fHistoPosMapRecalib[i]  ->SetDirectory(0);
  }
  delete contRecalibTPC;
  AliInfo(Form("dEdx recalibration maps successfully loaded from OADB/PWGGA/TPCdEdxRecalibOADB.root for run %d",runNumber));
  return kTRUE;

}

//_________________________________________________________________________
Double_t AliConversionPhotonCuts::GetCorrectedElectronTPCResponse(Short_t charge, Double_t nsig, Double_t P, Double_t Eta, Double_t TPCCl, Double_t R){

  Double_t Charge  = charge;
  Double_t CornSig = nsig;
  Double_t mean          = 1.;
  Double_t width         = 1.;
  //X axis 12 Y axis 18  ... common for all R slice
  Int_t BinP    = 4;   //  default value
  Int_t BinEta  = 9;  //  default value
  Int_t BinCorr   = -1;
  Double_t arrayTPCCl[5]    = {0, 60, 100, 150, 180};
  Double_t arrayR[5]        = {0, 33.5, 72, 145, 180};

  for (Int_t i = 0; i<4; i++){
    if (fIsRecalibDepTPCCl) {
      if (TPCCl > arrayTPCCl[i] && TPCCl <= arrayTPCCl[i+1]){
        BinCorr = i;
      }
    } else {
      if (R > arrayR[i] && R <= arrayR[i+1]){
        BinCorr = i;
      }
    }
  }
  if (BinCorr == -1 || BinCorr >  3){
    if (fIsRecalibDepTPCCl) cout<< " no valid TPC cluster number ..., not recalibrating"<< endl;
    else cout<< " no valid R bin number ..., not recalibrating"<< endl;
    return CornSig;// do nothing if correction map is not avaible
  }

  if(Charge<0){
    if (fHistoEleMapRecalib[BinCorr] == NULL ){
      cout<< " histograms are null..., going out"<< endl;
      return CornSig;// do nothing if correction map is not avaible
    }

    BinP    = fHistoEleMapRecalib[BinCorr]->GetXaxis()->FindBin(P);
    BinEta  = fHistoEleMapRecalib[BinCorr]->GetYaxis()->FindBin(Eta);

    if(P>0. && P<10.){
      mean  = (Double_t)fHistoEleMapRecalib[BinCorr]->GetBinContent(BinP,BinEta)/1000;
      width = (Double_t)fHistoEleMapRecalib[BinCorr]->GetBinError(BinP,BinEta)/1000;
    }else if(P>=10.){// use bin edge value
      mean  = (Double_t)fHistoEleMapRecalib[BinCorr]->GetBinContent(fHistoEleMapRecalib[BinCorr]->GetNbinsX(),BinEta)/1000;
      width = (Double_t)fHistoEleMapRecalib[BinCorr]->GetBinError(fHistoEleMapRecalib[BinCorr]->GetNbinsX(),BinEta)/1000;
    }

  }else{
    if (fHistoPosMapRecalib[BinCorr] == NULL ){
      cout<< " histograms are null..., going out"<< endl;
      return CornSig;// do nothing if correction map is not avaible
    }

    BinP    = fHistoPosMapRecalib[BinCorr]->GetXaxis()->FindBin(P);
    BinEta  = fHistoPosMapRecalib[BinCorr]->GetYaxis()->FindBin(Eta);

    if(P>0. && P<10.){
      mean  = (Double_t)fHistoPosMapRecalib[BinCorr]->GetBinContent(BinP,BinEta)/1000;
      width = (Double_t)fHistoPosMapRecalib[BinCorr]->GetBinError(BinP,BinEta)/1000;
    }else if(P>=10.){// use bin edge value
      mean  = (Double_t)fHistoPosMapRecalib[BinCorr]->GetBinContent(fHistoPosMapRecalib[BinCorr]->GetNbinsX(),BinEta)/1000;
      width = (Double_t)fHistoPosMapRecalib[BinCorr]->GetBinError(fHistoPosMapRecalib[BinCorr]->GetNbinsX(),BinEta)/1000;
    }
  }
  if (width!=0.){
    CornSig = (nsig - mean) / width;
  }
  return CornSig;
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PhotonIsSelectedMC(TParticle *particle,AliMCEvent *mcEvent,Bool_t checkForConvertedGamma){
  // MonteCarlo Photon Selection

  if(!mcEvent)return kFALSE;

  if (particle->GetPdgCode() == 22){


    if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) )
      return kFALSE;
        if(fEtaCutMin>-0.1){
            if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) )
                return kFALSE;
    }

    if(particle->GetMother(0) >-1 && mcEvent->Particle(particle->GetMother(0))->GetPdgCode() == 22){
      return kFALSE; // no photon as mothers!
    }

    // removed, decision on primary and secondary taken in main task
// 		if(particle->GetMother(0) >= mcEvent->GetNumberOfPrimaries()){
// 			return kFALSE; // the gamma has a mother, and it is not a primary particle
// 		}

    if(!checkForConvertedGamma) return kTRUE; // return in case of accepted gamma

    // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
    TParticle* ePos = NULL;
    TParticle* eNeg = NULL;

    if(particle->GetNDaughters() >= 2){
      for(Int_t daughterIndex=particle->GetFirstDaughter();daughterIndex<=particle->GetLastDaughter();daughterIndex++){
        if(daughterIndex<0) continue;
        TParticle *tmpDaughter = mcEvent->Particle(daughterIndex);
        if(tmpDaughter->GetUniqueID() == 5){
        if(tmpDaughter->GetPdgCode() == 11){
          eNeg = tmpDaughter;
        } else if(tmpDaughter->GetPdgCode() == -11){
          ePos = tmpDaughter;
        }
        }
      }
    }

    if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
      return kFALSE;
    }

    if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ||
      eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) )
      return kFALSE;

    if(fEtaCutMin > -0.1){
      if( (ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin)) ||
        (eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin)) )
        return kFALSE;
    }

    if(ePos->R()>fMaxR){
      return kFALSE; // cuts on distance from collision point
    }

    if(TMath::Abs(ePos->Vz()) > fMaxZ){
      return kFALSE;  // outside material
    }
    if(TMath::Abs(eNeg->Vz()) > fMaxZ){
      return kFALSE;  // outside material
    }

    if( ePos->R() <= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE;  // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   ePos->R() >= ((TMath::Abs(ePos->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    if( eNeg->R() <= ((TMath::Abs(eNeg->Vz()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE; // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   eNeg->R() >= ((TMath::Abs(eNeg->Vz()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    return kTRUE;
    //if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
  }
  return kFALSE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PhotonIsSelectedAODMC(AliAODMCParticle *particle,TClonesArray *aodmcArray,Bool_t checkForConvertedGamma){
  // MonteCarlo Photon Selection

  if(!aodmcArray)return kFALSE;

  if (particle->GetPdgCode() == 22){
    if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) )
      return kFALSE;
    if(fEtaCutMin>-0.1){
      if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) )
        return kFALSE;
    }

    if(particle->GetMother() > -1 && (static_cast<AliAODMCParticle*>(aodmcArray->At(particle->GetMother())))->GetPdgCode() == 22){
        return kFALSE; // no photon as mothers!
    }

    if(!checkForConvertedGamma) return kTRUE; // return in case of accepted gamma

    // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
    AliAODMCParticle* ePos = NULL;
    AliAODMCParticle* eNeg = NULL;

    if(particle->GetNDaughters() >= 2){
      for(Int_t daughterIndex=particle->GetDaughterLabel(0);daughterIndex<=particle->GetDaughterLabel(1);daughterIndex++){
        AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(aodmcArray->At(daughterIndex));
        if(!tmpDaughter) continue;
        if(((tmpDaughter->GetMCProcessCode())) == 5){    // STILL A BUG IN ALIROOT >>8 HAS TPO BE REMOVED AFTER FIX
        if(tmpDaughter->GetPdgCode() == 11){
          eNeg = tmpDaughter;
        } else if(tmpDaughter->GetPdgCode() == -11){
          ePos = tmpDaughter;
        }
        }
      }
    }

    if(ePos == NULL || eNeg == NULL){ // means we do not have two daughters from pair production
      return kFALSE;
    }

    if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ||
      eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) )
      return kFALSE;

    if(fEtaCutMin > -0.1){
      if( (ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin)) ||
        (eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin)) )
        return kFALSE;
    }

    Double_t rPos = sqrt( (ePos->Xv()*ePos->Xv()) + (ePos->Yv()*ePos->Yv()) );
    Double_t rNeg = sqrt( (eNeg->Xv()*eNeg->Xv()) + (eNeg->Yv()*eNeg->Yv()) );

    if(rPos>fMaxR){
      return kFALSE; // cuts on distance from collision point
    }
    if(TMath::Abs(ePos->Zv()) > fMaxZ){
      return kFALSE;  // outside material
    }
    if(TMath::Abs(eNeg->Zv()) > fMaxZ){
      return kFALSE;  // outside material
    }

    if( rPos <= ((TMath::Abs(ePos->Zv()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE;  // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   rPos >= ((TMath::Abs(ePos->Zv()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    if( rNeg <= ((TMath::Abs(eNeg->Zv()) * fLineCutZRSlope) - fLineCutZValue)){
      return kFALSE; // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   rNeg >= ((TMath::Abs(eNeg->Zv()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
      return kFALSE;
    }

    return kTRUE;
    //if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
  }
  return kFALSE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PhotonIsSelectedMCAODESD(AliDalitzAODESDMC* particle,AliDalitzEventMC *mcEvent,Bool_t checkForConvertedGamma) const{
// MonteCarlo Photon Selection
    if(!mcEvent)return kFALSE;
    if (particle->GetPdgCodeG() == 22){
        if( particle->EtaG() > (fEtaCut) || particle->EtaG() < (-fEtaCut) )
            return kFALSE;
        if(fEtaCutMin>-0.1){
        if( particle->EtaG() < (fEtaCutMin) && particle->EtaG() > (-fEtaCutMin) )
            return kFALSE;
        }
        std::unique_ptr<AliDalitzAODESDMC> Templeak;
        Templeak = std::unique_ptr<AliDalitzAODESDMC>(mcEvent->Particle(particle->GetMotherG()));

        if(particle->GetMotherG() >-1 && Templeak->GetPdgCodeG() == 22){
            return kFALSE; // no photon as mothers!
        }
        if(!checkForConvertedGamma) return kTRUE; // return in case of accepted gamma
        // looking for conversion gammas (electron + positron from pairbuilding (= 5) )
        std::unique_ptr<AliDalitzAODESDMC> ePos=0x0;
        std::unique_ptr<AliDalitzAODESDMC> eNeg=0x0;
        if(particle->GetNDaughtersG() >= 2){
        //      cout<<particle->GetNDaughtersG()<<endl;
            for(Int_t daughterIndex=particle->GetFirstDaughterG();daughterIndex<=particle->GetLastDaughterG();daughterIndex++){
                if(daughterIndex<0) continue;
                    std::unique_ptr<AliDalitzAODESDMC> tmpDaughter = std::unique_ptr<AliDalitzAODESDMC>(mcEvent->Particle(daughterIndex));
                    // cout<<tmpDaughter->GetPdgCodeG()<<endl;
                    //NOTE 8 Marzo problem here. never and ID 5
                    if(tmpDaughter->GetUniqueIDG() == 5){
                        if(tmpDaughter->GetPdgCodeG() == 11){
                        eNeg = std::unique_ptr<AliDalitzAODESDMC>(mcEvent->Particle(daughterIndex));
                    } else if(tmpDaughter->GetPdgCodeG() == -11){
                    ePos = std::unique_ptr<AliDalitzAODESDMC>(mcEvent->Particle(daughterIndex));
                        }
                    }
            }
        }
    if(ePos.get() == NULL || eNeg.get() == NULL){ // means we do not have two daughters from pair production
        return kFALSE;
    }
    if( ePos->EtaG() > (fEtaCut) || ePos->EtaG() < (-fEtaCut) ||
        eNeg->EtaG() > (fEtaCut) || eNeg->EtaG() < (-fEtaCut) )
        return kFALSE;
    if(fEtaCutMin > -0.1){
        if( (ePos->EtaG() < (fEtaCutMin) && ePos->EtaG() > (-fEtaCutMin)) ||
            (eNeg->EtaG() < (fEtaCutMin) && eNeg->EtaG() > (-fEtaCutMin)) )
            return kFALSE;
            }
    if(ePos->GetRatioVxyG()>fMaxR){
        return kFALSE; // cuts on distance from collision point
    }
        //cout<<" Paso RadioXY "<<endl;
    if(TMath::Abs(ePos->VertexOnZ()) > fMaxZ){
        return kFALSE;  // outside material
    }
    if(TMath::Abs(eNeg->VertexOnZ()) > fMaxZ){
        return kFALSE;  // outside material
    }
    if( ePos->GetRatioVxyG() <= ((TMath::Abs(ePos->VertexOnZ()) * fLineCutZRSlope) - fLineCutZValue)){
        return kFALSE;  // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   ePos->GetRatioVxyG() >= ((TMath::Abs(ePos->VertexOnZ()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
        return kFALSE;
    }
    if( eNeg->GetRatioVxyG() <= ((TMath::Abs(eNeg->VertexOnZ()) * fLineCutZRSlope) - fLineCutZValue)){
        return kFALSE; // line cut to exclude regions where we do not reconstruct
    } else if ( fEtaCutMin != -0.1 &&   eNeg->GetRatioVxyG() >= ((TMath::Abs(eNeg->VertexOnZ()) * fLineCutZRSlopeMin) - fLineCutZValueMin)){
        return kFALSE;
    }
        return kTRUE;
    //if(AcceptanceCut(particle,ePos,eNeg))return kTRUE;
    }

    return kFALSE;
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PhotonCuts(AliConversionPhotonBase *photon,AliVEvent *event){   // Specific Photon Cuts

  Int_t cutIndex = 0;
  if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt());
  cutIndex++;

  AliVTrack * electronCandidate = GetTrack(event,photon->GetTrackLabelNegative());
  AliVTrack * positronCandidate = GetTrack(event,photon->GetTrackLabelPositive());

  // Fill Histos before Cuts
  if(fHistoInvMassbefore)fHistoInvMassbefore->Fill(photon->GetMass());
  if(fHistoArmenterosbefore)fHistoArmenterosbefore->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());
  if(fHistoAsymmetrybefore){
    if(photon->GetPhotonP()!=0 && electronCandidate->P()!=0)fHistoAsymmetrybefore->Fill(photon->GetPhotonP(),electronCandidate->P()/photon->GetPhotonP());
  }
  // Gamma selection based on QT from Armenteros
  if(fDoQtGammaSelection == 1 || fDoQtGammaSelection == 2){
    if(!ArmenterosQtCut(photon)){
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //1
      return kFALSE;
    }
  }
  cutIndex++; //2

  // Chi Cut
  if(photon->GetChi2perNDF() > fChi2CutConversion || photon->GetChi2perNDF() <=0){
    {
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //2
      return kFALSE;
    }
  }
  cutIndex++;//3

  // Reconstruction Acceptance Cuts
  if(!AcceptanceCuts(photon)){
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //3
    return kFALSE;
  }

  cutIndex++; //4
  // Asymmetry Cut
  if(fDoPhotonAsymmetryCut == kTRUE){
    if(!AsymmetryCut(photon,event)){
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //4
      return kFALSE;
    }
  }

  //Check the pid probability
  cutIndex++; //5
  if(!PIDProbabilityCut(photon, event)) {
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //5
    return kFALSE;
  }

  cutIndex++; //6
  if(!CorrectedTPCClusterCut(photon, event)) {
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //6
    return kFALSE;
  }

  Double_t magField = event->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  } else {
    magField =  -1.0;
  }

  Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->Phi()-positronCandidate->Phi());

  cutIndex++; //7
  if(!PsiPairCut(photon)) {
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //7
    return kFALSE;
  }

  cutIndex++; //8
  if(!CosinePAngleCut(photon, event)) {
    if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //8
    return kFALSE;
  }

  AliAODConversionPhoton* photonAOD = dynamic_cast<AliAODConversionPhoton*>(photon);
  if (photonAOD){
    photonAOD->CalculateDistanceOfClossetApproachToPrimVtx(event->GetPrimaryVertex());

    cutIndex++; //9
    if(photonAOD->GetDCArToPrimVtx() > fDCARPrimVtxCut) { //DCA R cut of photon to primary vertex
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //9
      return kFALSE;
    }

    cutIndex++; //10
    if(TMath::Abs(photonAOD->GetDCAzToPrimVtx()) > fDCAZPrimVtxCut) { //DCA Z cut of photon to primary vertex
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //10
      return kFALSE;
    }
  } else {
    cutIndex++; //9
    cutIndex++; //10
  }
  cutIndex++; //11

  if (photonAOD){
    UChar_t photonQuality = 0;
    UChar_t photonQualityTOF = 0;
    UChar_t photonQualityTRD = 0;
    AliAODEvent * aodEvent = dynamic_cast<AliAODEvent*>(event);
    if(aodEvent) {
      photonQuality = DeterminePhotonQualityAOD(photonAOD, event);
    } else {
      photonQuality = photonAOD->GetPhotonQuality();
    }
    if (fPhotonQualityCutTOF) {
      photonQualityTOF = DeterminePhotonQualityTOF(photonAOD, event);
    }
    if (fPhotonQualityCutTRD) {
      photonQualityTRD = DeterminePhotonQualityTRD(photonAOD, event);
    }


    // If fPhotonQualityCutTRD == 0, the TRD part has no effect. Otherwise, selection takes place according to fPhotonQualityCutTRD
    if (fDoPhotonQualitySelectionCut && !(photonQuality == fPhotonQualityCut &&
                                        (!fPhotonQualityCutTRD || (fPhotonQualityCutTRD && photonQualityTRD == fPhotonQualityCutTRD)) &&
                                        (!fPhotonQualityCutTOF || (fPhotonQualityCutTOF && photonQualityTOF == fPhotonQualityCutTOF)))){
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //11
      return kFALSE;
    }
    if (fDoPhotonQualityRejectionCut && (photonQuality == fPhotonQualityCut &&
                                        (!fPhotonQualityCutTRD || (fPhotonQualityCutTRD && photonQualityTRD == fPhotonQualityCutTRD)) &&
                                        (!fPhotonQualityCutTOF || (fPhotonQualityCutTOF && photonQualityTOF == fPhotonQualityCutTOF)))){
      if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //11
      return kFALSE;
    }
  }
  cutIndex++; //12
  if(fHistoPhotonCuts)fHistoPhotonCuts->Fill(cutIndex, photon->GetPhotonPt()); //11

  // Histos after Cuts
  if(fHistoInvMassafter)fHistoInvMassafter->Fill(photon->GetMass());
  if(fHistoArmenterosafter)fHistoArmenterosafter->Fill(photon->GetArmenterosAlpha(),photon->GetArmenterosQt());
  if(fHistoPsiPairDeltaPhiafter)fHistoPsiPairDeltaPhiafter->Fill(deltaPhi,photon->GetPsiPair());
  if(fHistoKappaafter)fHistoKappaafter->Fill(photon->GetPhotonPt(), GetKappaTPC(photon, event));
  if(fHistoAsymmetryafter){
    if(photon->GetPhotonP()!=0 && electronCandidate->P()!=0)fHistoAsymmetryafter->Fill(photon->GetPhotonP(),electronCandidate->P()/photon->GetPhotonP());
  }
  return kTRUE;

}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::CorrectedTPCClusterCut(AliConversionPhotonBase *photon, AliVEvent * event){   //Cut on corrected TPC Cluster Info

  AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
  AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

  if(!negTrack||!posTrack)return kFALSE;

  Double_t negclsToF=0;

  if (!fUseCorrectedTPCClsInfo ){
    if(negTrack->GetTPCNclsF()!=0){
      negclsToF = (Double_t)negTrack->GetNcls(1)/(Double_t)negTrack->GetTPCNclsF();}// Ncluster/Nfindablecluster
  }
  else {
    negclsToF = negTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
  }

  Double_t posclsToF = 0.;
  if (!fUseCorrectedTPCClsInfo ){
    if(posTrack->GetTPCNclsF()!=0){
      posclsToF = (Double_t)posTrack->GetNcls(1)/(Double_t)posTrack->GetTPCNclsF();
    }
  }else{
    posclsToF = posTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(photon->GetConversionRadius()));
  }

  if( negclsToF < fMinClsTPCToF || posclsToF < fMinClsTPCToF ){
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::TrackIsSelected(AliConversionPhotonBase *photon, AliVEvent * event){
  //Selection of Reconstructed Photons

  if(event->IsA()==AliESDEvent::Class()) {
    if(!SelectV0Finder( ( ((AliESDEvent*)event)->GetV0(photon->GetV0Index()))->GetOnFlyStatus() ) ){
      return kFALSE;
    }
  }

  // Get Tracks
  AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
  AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

  if(!negTrack || !posTrack) {
    return kFALSE;
  }

  // check if V0 from AliAODGammaConversion.root is actually contained in AOD by checking if V0 exists with same tracks
  if(event->IsA()==AliAODEvent::Class() && fPreSelCut && ( fIsHeavyIon != 1 )) {
    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(event);

    Bool_t bFound = kFALSE;
    Int_t v0PosID = posTrack->GetID();
    Int_t v0NegID = negTrack->GetID();
    AliAODv0* v0 = NULL;
    for(Int_t iV=0; iV<aodEvent->GetNumberOfV0s(); iV++){
      v0 = aodEvent->GetV0(iV);
      if(!v0) continue;
      if( (v0PosID == v0->GetPosID() && v0NegID == v0->GetNegID()) || (v0PosID == v0->GetNegID() && v0NegID == v0->GetPosID()) ){
        bFound = kTRUE;
        break;
      }
    }
    if(!bFound){
      return kFALSE;
    }
  }

  photon->DeterminePhotonQuality(negTrack,posTrack);

  // Track Cuts
  if(!TracksAreSelected(negTrack, posTrack)){
    return kFALSE;
  }

  // Photon passed cuts
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PhotonIsSelected(AliConversionPhotonBase *photon, AliVEvent * event){
  //Selection of Reconstructed Photons

  FillPhotonCutIndex(kPhotonIn);

  if(event->IsA()==AliESDEvent::Class()) {
    if(!SelectV0Finder( ( ((AliESDEvent*)event)->GetV0(photon->GetV0Index()))->GetOnFlyStatus() ) ){
      FillPhotonCutIndex(kOnFly);
      return kFALSE;
    }
  }

  // Get Tracks
  AliVTrack * negTrack = GetTrack(event, photon->GetTrackLabelNegative());
  AliVTrack * posTrack = GetTrack(event, photon->GetTrackLabelPositive());

  if(!negTrack || !posTrack) {
    FillPhotonCutIndex(kNoTracks);
    return kFALSE;
  }

  // check if V0 from AliAODGammaConversion.root is actually contained in AOD by checking if V0 exists with same tracks
  if(event->IsA()==AliAODEvent::Class() && fPreSelCut && ( fIsHeavyIon != 1 || (fIsHeavyIon == 1 && fProcessAODCheck) )) {
    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(event);

    Bool_t bFound = kFALSE;
    Int_t v0PosID = posTrack->GetID();
    Int_t v0NegID = negTrack->GetID();
    AliAODv0* v0 = NULL;
    for(Int_t iV=0; iV<aodEvent->GetNumberOfV0s(); iV++){
      v0 = aodEvent->GetV0(iV);
      if(!v0) continue;
      if( (v0PosID == v0->GetPosID() && v0NegID == v0->GetNegID()) || (v0PosID == v0->GetNegID() && v0NegID == v0->GetPosID()) ){
        bFound = kTRUE;
        break;
      }
    }
    if(!bFound){
      FillPhotonCutIndex(kNoV0);
      return kFALSE;
    }
  }

  photon->DeterminePhotonQuality(negTrack,posTrack);

  // Track Cuts
  if(!TracksAreSelected(negTrack, posTrack)){
    FillPhotonCutIndex(kTrackCuts);
    return kFALSE;
  }
  if (fHistoEtaDistV0s)fHistoEtaDistV0s->Fill(photon->GetPhotonEta());

  // dEdx Cuts
  if(!KappaCuts(photon, event) || !dEdxCuts(negTrack,photon) || !dEdxCuts(posTrack,photon)) {
    FillPhotonCutIndex(kdEdxCuts);
    return kFALSE;
  }

  if (fHistoEtaDistV0sAfterdEdxCuts)fHistoEtaDistV0sAfterdEdxCuts->Fill(photon->GetPhotonEta());
  // Photon Cuts
  if(!PhotonCuts(photon,event)){
    FillPhotonCutIndex(kPhotonCuts);
    return kFALSE;
  }

  if(fDoPlotTrackPID && fHistoTrackPID){
    Int_t negpidForTracking = (Int_t)negTrack->GetPIDForTracking();
    Int_t pospidForTracking = (Int_t)posTrack->GetPIDForTracking();
    // cout << "PID:  " <<  negpidForTracking << ",     "<< pospidForTracking << endl;
    fHistoTrackPID->Fill((Float_t)negpidForTracking,negTrack->Pt());
    fHistoTrackPID->Fill((Float_t)pospidForTracking,posTrack->Pt());
  }

  // Photon passed cuts
  FillPhotonCutIndex(kPhotonOut);
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::ArmenterosQtCut(AliConversionPhotonBase *photon){   // Armenteros Qt Cut
  if(fDo2DQt){
    if(fDoQtGammaSelection==1){
      if ( !(TMath::Power(photon->GetArmenterosAlpha()/fMaxPhotonAsymmetry,2)+TMath::Power(photon->GetArmenterosQt()/fQtMax,2) < 1) ){
        return kFALSE;
      }
    } else if(fDoQtGammaSelection==2){
      Float_t qtMaxPtDep = fQtPtMax*photon->GetPhotonPt();
      if (qtMaxPtDep > fQtMax)
        qtMaxPtDep      = fQtMax;
      if ( !(TMath::Power(photon->GetArmenterosAlpha()/fMaxPhotonAsymmetry,2)+TMath::Power(photon->GetArmenterosQt()/qtMaxPtDep,2) < 1) ){
        return kFALSE;
      }
    }
  } else {
    if(fDoQtGammaSelection==1){
      if(photon->GetArmenterosQt()>fQtMax){
        return kFALSE;
      }
    } else if(fDoQtGammaSelection==2){
      Float_t qtMaxPtDep = fQtPtMax*photon->GetPhotonPt();
      if (qtMaxPtDep > fQtMax)
        qtMaxPtDep      = fQtMax;
      if(photon->GetArmenterosQt()>qtMaxPtDep){
        return kFALSE;
      }
    }
  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::AcceptanceCuts(AliConversionPhotonBase *photon) {
  // Exclude certain areas for photon reconstruction

  Int_t cutIndex=0;
  if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
  cutIndex++;

  if(photon->GetConversionRadius()>fMaxR){ // cuts on distance from collision point
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;

  if(photon->GetConversionRadius()<fMinR){ // cuts on distance from collision point
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;

  if(photon->GetConversionRadius()>fExcludeMinR && photon->GetConversionRadius()<fExcludeMaxR ){ // cuts on distance from collision point
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;


  if(photon->GetConversionRadius() <= ((TMath::Abs(photon->GetConversionZ())*fLineCutZRSlope)-fLineCutZValue)){
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  else if (fUseEtaMinCut &&  photon->GetConversionRadius() >= ((TMath::Abs(photon->GetConversionZ())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;

  if(TMath::Abs(photon->GetConversionZ()) > fMaxZ ){ // cuts out regions where we do not reconstruct
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;


  if( photon->GetPhotonEta() > (fEtaCut)    || photon->GetPhotonEta() < (-fEtaCut) ){
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  if(fEtaCutMin>-0.1){
    if( photon->GetPhotonEta() < (fEtaCutMin) && photon->GetPhotonEta() > (-fEtaCutMin) ){
      if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
      return kFALSE;
    }
  }
  cutIndex++;

  if (fDoShrinkTPCAcceptance == 1){
    if(photon->GetPhotonEta() > fEtaForPhiCutMin && photon->GetPhotonEta() < fEtaForPhiCutMax ){
      if (fMinPhiCut < fMaxPhiCut){
        if( photon->GetPhotonPhi() > fMinPhiCut && photon->GetPhotonPhi() < fMaxPhiCut ) {
          if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
          return kFALSE;
        }
      } else {
        Double_t photonPhi = photon->GetPhotonPhi();
        if (photon->GetPhotonPhi() < TMath::Pi()) photonPhi = photon->GetPhotonPhi() + 2*TMath::Pi();
        if( photonPhi > fMinPhiCut && photonPhi < fMaxPhiCut+2*TMath::Pi() ) {
          if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
          return kFALSE;
        }
      }
    }
  } else if (fDoShrinkTPCAcceptance == 2){  // accept only photons in 'good region'
      Double_t photonPhi = photon->GetPhotonPhi();
      GetPhiRegions();
      if( photon->GetPhotonEta()>0 && photon->GetPhotonEta()<fEtaCut ){        // A side
          //cout << "A side, eta=" << photon->GetPhotonEta() <<  endl;
          if(!(photonPhi>fGoodRegionAMin && photonPhi<fGoodRegionAMax)){
              //cout  << "photonPhi=" << photonPhi << " excluded" << endl;
              if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
              return kFALSE;
          } //else  cout  << "photonPhi=" << photonPhi << " accepted" << endl;
      } else if(photon->GetPhotonEta()<0 && photon->GetPhotonEta()>-fEtaCut){  // C side
          //cout << "C side, eta=" << photon->GetPhotonEta() <<  endl;
          if (!(photonPhi>fGoodRegionCMin && photonPhi<fGoodRegionCMax)){
              //cout  << "photonPhi=" << photonPhi << " excluded" << endl;
              if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
              return kFALSE;
          } //else  cout  << "photonPhi=" << photonPhi << " accepted" << endl;
      }
  } else if (fDoShrinkTPCAcceptance == 3){   // accept only photons in 'bad region'
      Double_t photonPhi = photon->GetPhotonPhi();
      GetPhiRegions();
      if( photon->GetPhotonEta()>0 && photon->GetPhotonEta()<fEtaCut ){        // A side
          //cout << "A side, eta=" << photon->GetPhotonEta() <<  endl;
          if(!(photonPhi>fBadRegionAMin && photonPhi<fBadRegionAMax)){
              //cout  << "photonPhi=" << photonPhi << " excluded" << endl;
              if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
              return kFALSE;
          } // else cout  << "photonPhi=" << photonPhi << " accepted" << endl;
      } else if(photon->GetPhotonEta()<0 && photon->GetPhotonEta()>-fEtaCut){  // C side
          //cout << "C side, eta=" << photon->GetPhotonEta() <<  endl;
          if (!(photonPhi>fBadRegionCMin && photonPhi<fBadRegionCMax)){
              //cout  << "photonPhi=" << photonPhi << " excluded" << endl;
              if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
              return kFALSE;
          } // else cout  << "photonPhi=" << photonPhi << " accepted" << endl;
      }
  }
  cutIndex++;



  if(photon->GetPhotonPt()<fPtCut){
    if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());
    return kFALSE;
  }
  cutIndex++;

  if(fHistoAcceptanceCuts)fHistoAcceptanceCuts->Fill(cutIndex, photon->GetPhotonPt());

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SpecificTrackCuts(AliAODTrack * negTrack, AliAODTrack * posTrack,Int_t &cutIndex) {
  // Track Cuts which require AOD/ESD specific implementation

  if( !negTrack->IsOn(AliESDtrack::kTPCrefit)  || !posTrack->IsOn(AliESDtrack::kTPCrefit)   )  {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  AliAODVertex * NegVtxType=negTrack->GetProdVertex();
  AliAODVertex * PosVtxType=posTrack->GetProdVertex();
  if( (NegVtxType->GetType())==AliAODVertex::kKink || (PosVtxType->GetType())==AliAODVertex::kKink) {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  // TPC Chi2 cut
  Double_t tpcNClsNeg = negTrack->GetTPCNcls();     // number of TPC clusters
  Double_t tpcNClsPos = posTrack->GetTPCNcls();
  Double_t tpcChi2NDFNeg = negTrack->Chi2perNDF();  // TPC-Chi2 / (nCls -5)  if nCls > 5 otherwise -1
  Double_t tpcChi2NDFPos = posTrack->Chi2perNDF();
  Double_t tpcChi2NDFNegCorr = (tpcNClsNeg>5)?tpcChi2NDFNeg*(tpcNClsNeg-5)/tpcNClsNeg:-1.;    // TPC-Chi2 / nCls  if nCls > 5 otherwise -1
  Double_t tpcChi2NDFPosCorr = (tpcNClsPos>5)?tpcChi2NDFPos*(tpcNClsPos-5)/tpcNClsPos:-1.;    // condition?IfYesThenDoThis:IfNoThenDoThis

  if(fHistoTPCChi2NDFBefore) fHistoTPCChi2NDFBefore->Fill(tpcChi2NDFNegCorr);
  if(fHistoTPCChi2NDF2D)     fHistoTPCChi2NDF2D->Fill(tpcChi2NDFNegCorr, tpcChi2NDFPosCorr);
  if(fMaxTPCChi2NDF>0){ // apply cut
      if(tpcChi2NDFNegCorr > fMaxTPCChi2NDF || tpcChi2NDFPosCorr > fMaxTPCChi2NDF || tpcChi2NDFNegCorr < 0.2 || tpcChi2NDFPosCorr < 0.2){
          if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
          return kFALSE;
      }
  }
  if(fHistoTPCChi2NDFAfter) fHistoTPCChi2NDFAfter->Fill(tpcChi2NDFNegCorr);
  // cutindex is incremented in TracksAreSelected after SpecificTrackCuts was called

  return kTRUE;

}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SpecificTrackCuts(AliESDtrack * negTrack, AliESDtrack * posTrack,Int_t &cutIndex) {
  // Track Cuts which require AOD/ESD specific implementation

  if( !negTrack->IsOn(AliESDtrack::kTPCrefit)  || !posTrack->IsOn(AliESDtrack::kTPCrefit)   )  {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  if(negTrack->GetKinkIndex(0) > 0  || posTrack->GetKinkIndex(0) > 0 ) {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  // TPC Chi2 cut
  Double_t tpcChi2Neg = (Double_t) negTrack->GetTPCchi2();   // TPC-Chi2
  Double_t tpcChi2Pos = (Double_t) posTrack->GetTPCchi2();
  Double_t tpcNClsNeg = negTrack->GetTPCNcls();              // number of TPC clusters
  Double_t tpcNClsPos = posTrack->GetTPCNcls();
  Double_t tpcChi2NDFNegCorr = (tpcNClsNeg>5)?tpcChi2Neg/tpcNClsNeg:-1.;    // TPC-Chi2 / nCls  if nCls > 5 otherwise -1
  Double_t tpcChi2NDFPosCorr = (tpcNClsPos>5)?tpcChi2Pos/tpcNClsPos:-1.;    // condition?IfYesThenDoThis:IfNoThenDoThis

  if(fHistoTPCChi2NDFBefore) fHistoTPCChi2NDFBefore->Fill(tpcChi2NDFNegCorr);
  if(fHistoTPCChi2NDF2D) fHistoTPCChi2NDF2D->Fill(tpcChi2NDFNegCorr, tpcChi2NDFPosCorr);
  if(fMaxTPCChi2NDF>0){ // do cut
      if(tpcChi2NDFNegCorr > fMaxTPCChi2NDF || tpcChi2NDFPosCorr > fMaxTPCChi2NDF || tpcChi2NDFNegCorr < 0.2 || tpcChi2NDFPosCorr < 0.2){
          if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
          return kFALSE;
      }
  }
  if(fHistoTPCChi2NDFAfter) fHistoTPCChi2NDFAfter->Fill(tpcChi2NDFNegCorr);
  // cutindex is incremented in TracksAreSelected after SpecificTrackCuts was called

  return kTRUE;
}



///________________________________________________________________________
Bool_t AliConversionPhotonCuts::TracksAreSelected(AliVTrack * negTrack, AliVTrack * posTrack) {
  // Track Selection for Photon Reconstruction

  Int_t cutIndex=0;
  if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
  cutIndex++;

  // avoid like sign
  if(fUseOnFlyV0FinderSameSign==0){
    if(negTrack->Charge() == posTrack->Charge()) {
      if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //1
      return kFALSE;
    }
  }else if(fUseOnFlyV0FinderSameSign==1){
    if(negTrack->Charge() != posTrack->Charge()) {
      if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //1
      return kFALSE;
    }
  }
  cutIndex++;

  // Number of TPC Clusters


  if( negTrack->GetNcls(1) < fMinClsTPC || posTrack->GetNcls(1) < fMinClsTPC ) {
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //2
    return kFALSE;
  }
  cutIndex++;

  // Acceptance
  if( posTrack->Eta() > (fEtaCut) || posTrack->Eta() < (-fEtaCut) ||
    negTrack->Eta() > (fEtaCut) || negTrack->Eta() < (-fEtaCut) ){
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //3
    return kFALSE;
  }
  if(fEtaCutMin>-0.1){
    if( (posTrack->Eta() < (fEtaCutMin) && posTrack->Eta() > (-fEtaCutMin)) ||
      (negTrack->Eta() < (fEtaCutMin) && negTrack->Eta() > (-fEtaCutMin)) ){
      if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //3
      return kFALSE;
    }
  }
  cutIndex++;

  // Single Pt Cut
  if(fDoAsymPtCut){
    if((posTrack->Pt()<fSinglePtCut || negTrack->Pt()<fSinglePtCut2) && (posTrack->Pt()<fSinglePtCut2 || negTrack->Pt()<fSinglePtCut) ){
      if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //4
      return kFALSE;
    }
  } else {
    if(posTrack->Pt()<fSinglePtCut || negTrack->Pt()<fSinglePtCut){
      if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //4
      return kFALSE;
    }
  }
  // TOF timing cut (fill at same place as single pT)
  if(fUseTOFtiming){
    if(fTOFtimingBothLegs){
      if( !((posTrack->GetStatus()&AliVTrack::kTOFout) && (negTrack->GetStatus()&AliVTrack::kTOFout)) ){ // no timing on both legs
        if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //4
        return kFALSE;
      } else if(
        (((posTrack->GetStatus()&AliVTrack::kTOFout) && (posTrack->GetTOFsignal()/1000 > fTOFtimeMax)) || posTrack->GetTOFsignal()/1000 < fTOFtimeMin) ||
        (((negTrack->GetStatus()&AliVTrack::kTOFout) && (negTrack->GetTOFsignal()/1000 > fTOFtimeMax)) || negTrack->GetTOFsignal()/1000 < fTOFtimeMin)
        ){ // timing outside of cut windows on either leg that has timing information
        if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //4
        return kFALSE;
      }
    } else {
      if( !((posTrack->GetStatus()&AliVTrack::kTOFout) || (negTrack->GetStatus()&AliVTrack::kTOFout)) ){ // timing on at least one leg
        if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //4
        return kFALSE;
      } else if(
        (((posTrack->GetStatus()&AliVTrack::kTOFout) && (posTrack->GetTOFsignal()/1000 > fTOFtimeMax)) || posTrack->GetTOFsignal()/1000 < fTOFtimeMin) ||
        (((negTrack->GetStatus()&AliVTrack::kTOFout) && (negTrack->GetTOFsignal()/1000 > fTOFtimeMax)) || negTrack->GetTOFsignal()/1000 < fTOFtimeMin)
        ){ // timing outside of cut windows on either leg that has timing information
        if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex); //4
        return kFALSE;
      }
    }
    // fill histogram with timing info (in ns) versus momentum
    if(posTrack->GetStatus()&AliVTrack::kTOFout){
      fHistoTOFtimeVSMomentum->Fill(posTrack->GetTOFsignal()/1000,posTrack->Pt());
    }
    if(negTrack->GetStatus()&AliVTrack::kTOFout){
      fHistoTOFtimeVSMomentum->Fill(negTrack->GetTOFsignal()/1000,negTrack->Pt());
    }
  }
  cutIndex++;

  // AOD ESD specific cuts
  Bool_t passCuts = kTRUE;

  if(negTrack->IsA()==AliAODTrack::Class()) {
    passCuts = SpecificTrackCuts(static_cast<AliAODTrack*>(negTrack), static_cast<AliAODTrack*>(posTrack),cutIndex); //5,6,7
  } else {
    passCuts = SpecificTrackCuts(static_cast<AliESDtrack*>(negTrack), static_cast<AliESDtrack*>(posTrack),cutIndex);
  }

  if(!passCuts){
    if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);
    return kFALSE;
  }
  cutIndex++;

  if(fHistoTrackCuts)fHistoTrackCuts->Fill(cutIndex);  // out

  return kTRUE;

}

///________________________________________________________________________
Float_t AliConversionPhotonCuts::GetKappaTPC(AliConversionPhotonBase *gamma, AliVEvent * event){

  if(!fPIDResponse){InitPIDResponse();}// Try to reinitialize PID Response
  if(!fPIDResponse){AliError("No PID Response"); return kTRUE;}// if still missing fatal error

  AliVTrack * negTrack = GetTrack(event, gamma->GetTrackLabelNegative());
  AliVTrack * posTrack = GetTrack(event, gamma->GetTrackLabelPositive());

  Double_t CentrnSig[2] ={-1.,-1.};//negative, positive
  Double_t P[2]         ={-1.,-1.};
  Double_t Eta[2]       ={-1.,-1.};

  Float_t KappaPlus, KappaMinus, Kappa;
  if(fDoElecDeDxPostCalibration){
    CentrnSig[0]=fPIDResponse->NumberOfSigmasTPC(negTrack,AliPID::kElectron);
    CentrnSig[1]=fPIDResponse->NumberOfSigmasTPC(posTrack,AliPID::kElectron);
    P[0]        =negTrack->P();
    P[1]        =posTrack->P();
    Eta[0]      =negTrack->Eta();
    Eta[1]      =posTrack->Eta();
    KappaMinus = GetCorrectedElectronTPCResponse(negTrack->Charge(),CentrnSig[0],P[0],Eta[0],negTrack->GetTPCNcls(),gamma->GetConversionRadius());
    KappaPlus =  GetCorrectedElectronTPCResponse(posTrack->Charge(),CentrnSig[1],P[1],Eta[1],posTrack->GetTPCNcls(),gamma->GetConversionRadius());
  }else{
    KappaMinus = fPIDResponse->NumberOfSigmasTPC(negTrack, AliPID::kElectron);
    KappaPlus =  fPIDResponse->NumberOfSigmasTPC(posTrack, AliPID::kElectron);
  }
  Kappa = ( TMath::Abs(KappaMinus) + TMath::Abs(KappaPlus) ) / 2.0 + 2.0*(KappaMinus+KappaPlus);

  return Kappa;

}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::GetBDTVariableValues(AliConversionPhotonBase *gamma, AliVEvent * event, Float_t* values){

  if(!fPIDResponse){InitPIDResponse();}// Try to reinitialize PID Response
  if(!fPIDResponse){AliError("No PID Response"); return kFALSE;}// if still missing fatal error

  AliVTrack * negTrack = GetTrack(event, gamma->GetTrackLabelNegative());
  AliVTrack * posTrack = GetTrack(event, gamma->GetTrackLabelPositive());
  Int_t nPosClusterITS = 0;
  Int_t nNegClusterITS = 0;
  for(Int_t itsLayer = 0; itsLayer<6;itsLayer++){
    if(TESTBIT(negTrack->GetITSClusterMap(),itsLayer)){
      nNegClusterITS++;
    }
    if(TESTBIT(posTrack->GetITSClusterMap(),itsLayer)){
      nPosClusterITS++;
    }
  }

  //"dEdxElectronITS + dEdxPositronITS"
  if (nPosClusterITS > 0 ){
    values[0] =  posTrack->GetITSsignal();
  } else {
    values[0] =  1000;
  }
  if (nNegClusterITS > 0 ){
    values[0] +=  negTrack->GetITSsignal();
  } else {
    values[0] +=  1000;
  }

  values[1]= (Float_t)posTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(gamma->GetConversionRadius())); //"fracClsTPCPositron"
  values[2]= (Float_t)negTrack->GetTPCClusterInfo(2,0,GetFirstTPCRow(gamma->GetConversionRadius())); //"fracClsTPCElectron"
  values[3]= nPosClusterITS; //"clsITSPositron"
  values[4]= nNegClusterITS; //"clsITSElectron"
  values[5]=fPIDResponse->NumberOfSigmasTPC(negTrack,AliPID::kElectron); //"nSigmaTPCElectron"
  values[6]=fPIDResponse->NumberOfSigmasTPC(posTrack,AliPID::kElectron); //"nSigmaTPCPositron"

  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::dEdxCuts(AliVTrack *fCurrentTrack,AliConversionPhotonBase* photon){
  // Supposed to use post calibration
  // Electron Identification Cuts for Photon reconstruction

  if(!fPIDResponse){InitPIDResponse();}// Try to reinitialize PID Response
  if(!fPIDResponse){AliError("No PID Response"); return kTRUE;}// if still missing fatal error

  Short_t Charge    = fCurrentTrack->Charge();
  Double_t electronNSigmaTPC = fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kElectron);
  Double_t electronNSigmaTPCCor=0.;
  Double_t P=0.;
  Double_t Eta=0.;

  if(fDoElecDeDxPostCalibration){
    P = fCurrentTrack->P();
    Eta = fCurrentTrack->Eta();
    electronNSigmaTPCCor = GetCorrectedElectronTPCResponse(Charge,electronNSigmaTPC,P,Eta,fCurrentTrack->GetTPCNcls(),photon->GetConversionRadius());
  }

  Int_t cutIndex=0;
  if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
  if(fHistoTPCdEdxSigbefore)fHistoTPCdEdxSigbefore->Fill(fCurrentTrack->P(), electronNSigmaTPC);
  if(fHistoTPCdEdxbefore)fHistoTPCdEdxbefore->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());
  cutIndex++; //1
  if(fDodEdxSigmaCut == kTRUE && !fSwitchToKappa){
    // TPC Electron Line
    if(fDoElecDeDxPostCalibration){
      if( electronNSigmaTPCCor < fPIDnSigmaBelowElectronLine ||  electronNSigmaTPCCor >fPIDnSigmaAboveElectronLine ){
        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    } else{
      if( electronNSigmaTPC < fPIDnSigmaBelowElectronLine || electronNSigmaTPC > fPIDnSigmaAboveElectronLine){
        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
    cutIndex++; //2
    // TPC Pion Line
    if( fCurrentTrack->P()>fPIDMinPnSigmaAbovePionLine && fCurrentTrack->P()<fPIDMaxPnSigmaAbovePionLine ){
      if(fDoElecDeDxPostCalibration){
        if( electronNSigmaTPCCor >fPIDnSigmaBelowElectronLine && electronNSigmaTPCCor < fPIDnSigmaAboveElectronLine && fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){
          if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
          return kFALSE;
        }
      } else{
        if( electronNSigmaTPC > fPIDnSigmaBelowElectronLine && electronNSigmaTPC < fPIDnSigmaAboveElectronLine && fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLine){
          if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
          return kFALSE;
        }
      }
    }
    cutIndex++; //3

    // High Pt Pion rej
    if( fCurrentTrack->P()>fPIDMaxPnSigmaAbovePionLine ){
      if(fDoElecDeDxPostCalibration){
        if( electronNSigmaTPCCor > fPIDnSigmaBelowElectronLine && electronNSigmaTPCCor < fPIDnSigmaAboveElectronLine && fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineHighPt){
          if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
          return kFALSE;
        }
      } else{
        if( electronNSigmaTPC > fPIDnSigmaBelowElectronLine && electronNSigmaTPC < fPIDnSigmaAboveElectronLine && fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion)<fPIDnSigmaAbovePionLineHighPt){
          if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
          return kFALSE;
        }
      }
    }
    cutIndex++; //4
  }
  else{cutIndex+=3;} //4

  if(fDoKaonRejectionLowP == kTRUE && !fSwitchToKappa){
    if(fCurrentTrack->P()<fPIDMinPKaonRejectionLowP ){
      if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kKaon))<fPIDnSigmaAtLowPAroundKaonLine){
        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
  }
  cutIndex++; //5

  if(fDoProtonRejectionLowP == kTRUE && !fSwitchToKappa){
    if( fCurrentTrack->P()<fPIDMinPProtonRejectionLowP ){
      if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kProton))<fPIDnSigmaAtLowPAroundProtonLine){
        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
  }
  cutIndex++; //6

  if(fDoPionRejectionLowP == kTRUE && !fSwitchToKappa){
    if( fCurrentTrack->P()<fPIDMinPPionRejectionLowP ){
      if( TMath::Abs(fPIDResponse->NumberOfSigmasTPC(fCurrentTrack,AliPID::kPion))<fPIDnSigmaAtLowPAroundPionLine){
        if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
        return kFALSE;
      }
    }
  }
  cutIndex++; //7

  //  if((fCurrentTrack->GetStatus() & AliESDtrack::kTOFpid ) && !(fCurrentTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
  if((fCurrentTrack->GetStatus() & AliVTrack::kTOFout ) && (fCurrentTrack->GetStatus() & AliVTrack::kTIME)){ // check for TOF signal
    if(fHistoTOFbefore){
      Double_t t0 = fPIDResponse->GetTOFResponse().GetStartTime(fCurrentTrack->P());
      Double_t  times[AliPID::kSPECIESC];
      fCurrentTrack->GetIntegratedTimes(times,AliPID::kSPECIESC);
      Double_t TOFsignal = fCurrentTrack->GetTOFsignal();
      Double_t dT = TOFsignal - t0 - times[0];
      fHistoTOFbefore->Fill(fCurrentTrack->P(),dT);
    }
    if(fHistoTOFSigbefore) fHistoTOFSigbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
    if(fUseTOFpid){
        if(!fUseTOFpidMinMom || (fUseTOFpidMinMom && fCurrentTrack->Pt() > fTofPIDMinMom)){
            if(fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)>fTofPIDnSigmaAboveElectronLine ||
               fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron)<fTofPIDnSigmaBelowElectronLine ){
                if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
                return kFALSE;
            }
        }
    }
    if(fHistoTOFSigafter)fHistoTOFSigafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasTOF(fCurrentTrack, AliPID::kElectron));
  }
  cutIndex++; //8

  if((fCurrentTrack->GetStatus() & AliESDtrack::kITSpid)){
    if(fHistoITSSigbefore) fHistoITSSigbefore->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron));
    if(fUseITSpid){
      if(fCurrentTrack->Pt()<=fMaxPtPIDITS){
        if(fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron)>fITSPIDnSigmaAboveElectronLine || fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron)<fITSPIDnSigmaBelowElectronLine ){
          if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
          return kFALSE;
        }
      }
    }
    if(fHistoITSSigafter)fHistoITSSigafter->Fill(fCurrentTrack->P(),fPIDResponse->NumberOfSigmasITS(fCurrentTrack, AliPID::kElectron));
  }

  cutIndex++; //9

  // Apply TRD PID
  if(fDoTRDPID){
    if(!fPIDResponse->IdentifiedAsElectronTRD(fCurrentTrack,fPIDTRDEfficiency)){
      if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
      return kFALSE;
    }
  }
  cutIndex++; //10

  if(fHistodEdxCuts)fHistodEdxCuts->Fill(cutIndex,fCurrentTrack->Pt());
  if(fDoElecDeDxPostCalibration){
    if(fHistoTPCdEdxSigafter)fHistoTPCdEdxSigafter->Fill(fCurrentTrack->P(),electronNSigmaTPCCor);
  }else{
    if(fHistoTPCdEdxSigafter)fHistoTPCdEdxSigafter->Fill(fCurrentTrack->P(),electronNSigmaTPC);
  }
  if(fHistoTPCdEdxafter)fHistoTPCdEdxafter->Fill(fCurrentTrack->P(),fCurrentTrack->GetTPCsignal());

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::KappaCuts(AliConversionPhotonBase * photon,AliVEvent *event) {
  // abort if Kappa selection not enabled
  if (!fSwitchToKappa) return kTRUE;

  Float_t kappa = GetKappaTPC(photon, event);
  if (kappa < fKappaMinCut) return kFALSE;
  if (kappa > fKappaMaxCut) return kFALSE;
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::AsymmetryCut(AliConversionPhotonBase * photon,AliVEvent *event) {
  // Cut on Energy Asymmetry

  for(Int_t ii=0;ii<2;ii++){

    AliVTrack *track=GetTrack(event,photon->GetTrackLabel(ii));

    if(fDoPhotonPDependentAsymCut){
      Double_t trackNegAsy=0;
      if (photon->GetPhotonP()!=0.){
          trackNegAsy= track->P()/photon->GetPhotonP();
      }

      if( trackNegAsy > fFAsymmetryCut->Eval(photon->GetPhotonP()) || trackNegAsy < 1.-fFAsymmetryCut->Eval(photon->GetPhotonP()) ){
        return kFALSE;
      }

    } else {
      if( track->P() > fMinPPhotonAsymmetryCut ){
        Double_t trackNegAsy=0;
        if (photon->GetPhotonP()!=0.){
          trackNegAsy= track->P()/photon->GetPhotonP();
        }

        if( trackNegAsy<fMinPhotonAsymmetry ||trackNegAsy>(1.- fMinPhotonAsymmetry)){
          return kFALSE;
        }
      }
    }

  }
  return kTRUE;
}

///________________________________________________________________________
AliVTrack *AliConversionPhotonCuts::GetTrack(AliVEvent * event, Int_t label){
  //Returns pointer to the track with given ESD label
  //(Important for AOD implementation, since Track array in AOD data is different
  //from ESD array, but ESD tracklabels are stored in AOD Tracks)

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);
  if(esdEvent) {
    if(label > event->GetNumberOfTracks() ) return NULL;
    AliESDtrack * track = esdEvent->GetTrack(label);
    return track;

  } else {
    if(label == -999999) return NULL; // if AOD relabelling goes wrong, immediately return NULL
    AliVTrack * track = 0x0;
    if(AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()) && ((AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()))->AreAODsRelabeled()){
      if(event->GetTrack(label)) track = dynamic_cast<AliVTrack*>(event->GetTrack(label));
      return track;
    }
    else{
      for(Int_t ii=0; ii<event->GetNumberOfTracks(); ii++) {
        if(event->GetTrack(ii)) track = dynamic_cast<AliVTrack*>(event->GetTrack(ii));
        if(track){
          if(track->GetID() == label) {
            return track;
          }
        }
      }
    }
  }
  //AliDebug(5,(Form("track not found %d %d",label,event->GetNumberOfTracks()));
  return NULL;
}

///________________________________________________________________________
AliESDtrack *AliConversionPhotonCuts::GetESDTrack(AliESDEvent * event, Int_t label){
  //Returns pointer to the track with given ESD label
  //(Important for AOD implementation, since Track array in AOD data is different
  //from ESD array, but ESD tracklabels are stored in AOD Tracks)

  if(event) {
    if(label > event->GetNumberOfTracks() ) return NULL;
    AliESDtrack * track = event->GetTrack(label);
    return track;
  }
  //AliDebug(5,(Form("track not found %d %d",label,event->GetNumberOfTracks()));
  return NULL;
}



///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PIDProbabilityCut(AliConversionPhotonBase *photon, AliVEvent * event){
  // Cut on Electron Probability for Photon Reconstruction

  AliESDEvent * esdEvent = dynamic_cast<AliESDEvent*>(event);

  if(esdEvent){

    Bool_t iResult=kFALSE;

    Double_t posProbArray[AliPID::kSPECIES];
    Double_t negProbArray[AliPID::kSPECIES];

    AliESDtrack* negTrack   = esdEvent->GetTrack(photon->GetTrackLabelNegative());
    AliESDtrack* posTrack   = esdEvent->GetTrack(photon->GetTrackLabelPositive());

    if(negTrack && posTrack){

      negTrack->GetTPCpid(negProbArray);
      posTrack->GetTPCpid(posProbArray);

      if(negProbArray[AliPID::kElectron]>=fPIDProbabilityCutNegativeParticle && posProbArray[AliPID::kElectron]>=fPIDProbabilityCutPositiveParticle){
        iResult=kTRUE;
      }
    }

    return iResult;

  } else {
    ///Not possible for AODs
    return kTRUE;
  }
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg){
  // MC Acceptance Cuts
  //(Certain areas were excluded for photon reconstruction)

  if(particle->R()>fMaxR){
    return kFALSE;}

  if(ePos->R()>fMaxR){
    return kFALSE;
  }

  if(ePos->R()<fMinR){
    return kFALSE;
  }

  if( ePos->R() <= ((TMath::Abs(ePos->Vz())*fLineCutZRSlope)-fLineCutZValue)){
    return kFALSE;
  }
  else if (fUseEtaMinCut &&  ePos->R() >= ((TMath::Abs(ePos->Vz())*fLineCutZRSlopeMin)-fLineCutZValueMin )){
    return kFALSE;
  }

  if(TMath::Abs(eNeg->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
    return kFALSE;
  }

  if(eNeg->Vz()!=ePos->Vz()||eNeg->R()!=ePos->R()){
    return kFALSE;
  }

  if(TMath::Abs(ePos->Vz()) > fMaxZ ){ // cuts out regions where we do not reconstruct
    return kFALSE;
  }


  if( particle->Eta() > (fEtaCut) || particle->Eta() < (-fEtaCut) ){
    return kFALSE;
  }
  if( ePos->Eta() > (fEtaCut) || ePos->Eta() < (-fEtaCut) ){
    return kFALSE;
  }
  if( eNeg->Eta() > (fEtaCut) || eNeg->Eta() < (-fEtaCut) ){
    return kFALSE;
  }
  if(fEtaCutMin>-0.1){
    if( particle->Eta() < (fEtaCutMin) && particle->Eta() > (-fEtaCutMin) ){
      return kFALSE;
    }
    if( ePos->Eta() < (fEtaCutMin) && ePos->Eta() > (-fEtaCutMin) ){
      return kFALSE;
    }
    if( eNeg->Eta() < (fEtaCutMin) && eNeg->Eta() > (-fEtaCutMin) ){
      return kFALSE;
    }
  }

  if(fDoAsymPtCut){
      if((ePos->Pt()<fSinglePtCut || eNeg->Pt()<fSinglePtCut2) && (ePos->Pt()<fSinglePtCut2 || eNeg->Pt()<fSinglePtCut) ){
        return kFALSE;
      }
    } else {
      if(ePos->Pt()<fSinglePtCut || eNeg->Pt()<fSinglePtCut){
        return kFALSE;
      }
    }

  if(particle->Pt()<fPtCut){
    return kFALSE;
  }

  return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::UpdateCutString() {
  ///Update the cut string (if it has been created yet)

  if(fCutString && fCutString->GetString().Length() == kNCuts) {
    fCutString->SetString(GetCutNumber());
  } else {
    return kFALSE;
  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::InitializeCutsFromCutString(const TString analysisCutSelection ) {
  fCutStringRead = Form("%s",analysisCutSelection.Data());

  // Initialize Cuts from a given Cut string
  AliInfo(Form("Set Photoncut Number: %s",analysisCutSelection.Data()));
  if(analysisCutSelection.Length()!=kNCuts) {
    AliError(Form("Cut selection has the wrong length! size is %d, number of cuts is %d", analysisCutSelection.Length(), kNCuts));
    return kFALSE;
  }
  if(!analysisCutSelection.IsAlnum()){
    AliError("Cut selection is not alphanumeric");
    return kFALSE;
  }

  TString analysisCutSelectionLowerCase = Form("%s",analysisCutSelection.Data());
  analysisCutSelectionLowerCase.ToLower();
  const char *cutSelection = analysisCutSelectionLowerCase.Data();
  #define ASSIGNARRAY(i)  fCuts[i] = ((int)cutSelection[i]>=(int)'a') ? cutSelection[i]-'a'+10 : cutSelection[i]-'0'
  for(Int_t ii=0;ii<kNCuts;ii++){
    ASSIGNARRAY(ii);
  }

  // Set Individual Cuts
  for(Int_t ii=0;ii<kNCuts;ii++){
    if(!SetCut(cutIds(ii),fCuts[ii]))return kFALSE;
  }

  PrintCutsWithValues();

  return kTRUE;
}
///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetCut(cutIds cutID, const Int_t value) {
  ///Set individual cut ID

  switch (cutID) {

    case kv0FinderType:
      if( SetV0Finder(value)) {
        fCuts[kv0FinderType] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case ketaCut:
      if( SetEtaCut(value)) {
        fCuts[ketaCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kRCut:
      if( SetRCut(value)) {
        fCuts[kRCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kEtaForPhiSector:
      if( SetEtaForPhiCut(value)) {
        fCuts[kEtaForPhiSector] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    case kMinPhiSector:
      if( SetMinPhiSectorCut(value)) {
        fCuts[kMinPhiSector] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;
    case kMaxPhiSector:
      if( SetMaxPhiSectorCut(value)) {
        fCuts[kMaxPhiSector] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case ksinglePtCut:
      if( SetSinglePtCut(value)) {
        fCuts[ksinglePtCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kclsTPCCut:
      if( SetTPCClusterCut(value)) {
        fCuts[kclsTPCCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kededxSigmaCut:
      if (!fSwitchToKappa){
        if( SetTPCdEdxCutElectronLine(value)) {
          fCuts[kededxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        if( SetKappaTPCCut(value)) {
          fCuts[kededxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      }
    case kpidedxSigmaCut:
      if (!fSwitchToKappa){
        if( SetTPCdEdxCutPionLine(value)) {
          fCuts[kpidedxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        fCuts[kpidedxSigmaCut] = 0;
        return kTRUE;
      }
    case kpiMomdedxSigmaCut:
      if (!fSwitchToKappa){
        if( SetMinMomPiondEdxCut(value)) {
          fCuts[kpiMomdedxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        fCuts[kpiMomdedxSigmaCut] = 0;
        return kTRUE;
      }
    case kpiMaxMomdedxSigmaCut:
      if (!fSwitchToKappa){
        if( SetMaxMomPiondEdxCut(value)) {
          fCuts[kpiMaxMomdedxSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        fCuts[kpiMaxMomdedxSigmaCut] = 0;
        return kTRUE;
      }
    case kLowPRejectionSigmaCut:
      if (!fSwitchToKappa){
        if( SetLowPRejectionCuts(value)) {
          fCuts[kLowPRejectionSigmaCut] = value;
          UpdateCutString();
          return kTRUE;
        } else return kFALSE;
      } else {
        fCuts[kLowPRejectionSigmaCut] = 0;
        return kTRUE;
      }
    case kTOFelectronPID:
      if( SetTOFElectronPIDCut(value)) {
        fCuts[kTOFelectronPID] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kQtMaxCut:
      if( SetQtMaxCut(value)) {
        fCuts[kQtMaxCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;


    case kchi2GammaCut:
      if( SetChi2GammaCut(value)) {
        fCuts[kchi2GammaCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kPsiPair:
      if( SetPsiPairCut(value)) {
        fCuts[kPsiPair] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kdoPhotonAsymmetryCut:
      if( SetPhotonAsymmetryCut(value)) {
        fCuts[kdoPhotonAsymmetryCut] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kCosPAngle:
      if( SetCosPAngleCut(value)) {
        fCuts[kCosPAngle] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kElecShare:
      if( SetSharedElectronCut(value)) {
        fCuts[kElecShare] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kToCloseV0s:
      if( SetToCloseV0sCut(value)) {
        fCuts[kToCloseV0s] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kDcaRPrimVtx:
      if( SetDCARPhotonPrimVtxCut(value)) {
        fCuts[kDcaRPrimVtx] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kDcaZPrimVtx:
      if( SetDCAZPhotonPrimVtxCut(value)) {
        fCuts[kDcaZPrimVtx] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kInPlaneOutOfPlane:
      if( SetInPlaneOutOfPlane(value)) {
        fCuts[kInPlaneOutOfPlane] = value;
        UpdateCutString();
        return kTRUE;
      } else return kFALSE;

    case kITSelectronPID:
    if( SetITSElectronPIDCut(value)) {
      fCuts[kITSelectronPID] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

    case kTRDelectronPID:
    if( SetTRDElectronPIDCut(value)) {
      fCuts[kTRDelectronPID] = value;
      UpdateCutString();
      return kTRUE;
    } else return kFALSE;

    case kNCuts:
      AliError("Cut id out of range");
      return kFALSE;
  }

  AliError("Cut id %d not recognized");
  return kFALSE;
}
///________________________________________________________________________
void AliConversionPhotonCuts::PrintCuts() {
   // Print out current Cut Selection
   for(Int_t ic = 0; ic < kNCuts; ic++) {
      printf("%-30s : %d \n", fgkCutNames[ic], fCuts[ic]);
   }
}

void AliConversionPhotonCuts::PrintCutsWithValues() {
  // Print out current Cut Selection with value
  printf("\nConversion cutnumber \n");
  for(Int_t ic = 0; ic < kNCuts; ic++) {
      printf("%d",fCuts[ic]);  // careful with output: cannot hande characters in cut string
  }
  printf("\n\n");
  printf("Electron cuts & Secondary Track Cuts - only track from secondaries enter analysis: \n");
  printf("\t no like sign pairs from V0s \n");
  if (!fUseCorrectedTPCClsInfo) printf("\t # TPC clusters > %3.2f \n", fMinClsTPC);
  if (fEtaCutMin > -0.1) printf("\t %3.2f < eta_{e} < %3.2f\n", fEtaCutMin, fEtaCut );
    else printf("\t eta_{e} < %3.2f\n", fEtaCut );
  if(fDoShrinkTPCAcceptance == 1) printf("\t reject: %3.2f < phi < %3.2f with %3.2f < eta < %3.2f  \n", fMinPhiCut, fMaxPhiCut, fEtaForPhiCutMin, fEtaForPhiCutMax);
  else if (fDoShrinkTPCAcceptance == 2 ) printf("\t Only use photons in phi regions without distortions\n");
  else if (fDoShrinkTPCAcceptance == 3 ) printf("\t Only use photons in phi regions with strong distortions\n");
  else printf("\t No phi cut\n");
  if(fDoAsymPtCut)
    printf("\t Asymmetric cut: p_{T,e1} > %3.2f and p_{T,e2} > %3.2f\n", fSinglePtCut, fSinglePtCut2 );
  else
    printf("\t p_{T,e} > %3.2f\n", fSinglePtCut );
  printf("\t TPC refit \n");
  printf("\t no kinks \n");
  if (!fSwitchToKappa){
    printf("\t accept: %3.2f < n sigma_{e,TPC} < %3.2f\n", fPIDnSigmaBelowElectronLine, fPIDnSigmaAboveElectronLine );
    printf("\t reject: %3.2f < p_{e,T} < %3.2f, n sigma_{pi,TPC} < %3.2f\n", fPIDMinPnSigmaAbovePionLine, fPIDMaxPnSigmaAbovePionLine, fPIDnSigmaAbovePionLine );
    printf("\t reject: p_{e,T} > %3.2f, n sigma_{pi,TPC} < %3.2f\n", fPIDMaxPnSigmaAbovePionLine, fPIDnSigmaAbovePionLineHighPt );
    if (fDoPionRejectionLowP) printf("\t reject: p_{e,T} < %3.2f, -%3.2f < n sigma_{pi,TPC} < %3.2f\n", fPIDMinPPionRejectionLowP, fPIDnSigmaAtLowPAroundPionLine, fPIDnSigmaAtLowPAroundPionLine );
    if (fDoKaonRejectionLowP) printf("\t reject: -%3.2f < n sigma_{K,TPC} < %3.2f\n", fPIDnSigmaAtLowPAroundKaonLine, fPIDnSigmaAtLowPAroundKaonLine );
    if (fDoProtonRejectionLowP) printf("\t reject: -%3.2f < n sigma_{p,TPC} < %3.2f\n", fPIDnSigmaAtLowPAroundProtonLine, fPIDnSigmaAtLowPAroundProtonLine );
  } else {
    printf("\t accept: %3.2f <= Kappa_{TPC} < %3.2f\n", fKappaMinCut, fKappaMaxCut );
  }
  if (fUseTOFtiming){
    if(fTOFtimingBothLegs) printf("\t requiring TOF timing information on both electrons\n");
    else printf("\t requiring TOF timing information on single electron\n");
  }
  if (fUseTOFpid) printf("\t accept: %3.2f < n sigma_{e,TOF} < %3.2f\n", fTofPIDnSigmaBelowElectronLine, fTofPIDnSigmaAboveElectronLine);
  if (fUseITSpid) printf("\t accept: %3.2f < n sigma_{e,ITS} < %3.2f\n -- up to pT %3.2f", fITSPIDnSigmaBelowElectronLine, fITSPIDnSigmaAboveElectronLine, fMaxPtPIDITS);

  printf("Photon cuts: \n");
  if (fUseOnFlyV0Finder) printf("\t using Onfly V0 finder \n");
  else printf("\t using Offline V0 finder \n");
  if (fDo2DQt){
    printf("\t 2 dimensional q_{T} cut applied with maximum of %3.2f \n", fQtMax );
  } else {
    printf("\t 1 dimensional q_{T} cut applied with maximum of %3.2f \n", fQtMax );
  }
  if (fDo2DPsiPairChi2==1){
    printf("\t 2 dimensional triangle chi^{2} and psi_{pair} cut applied with maximum of chi^{2} = %3.2f and |psi_{pair}| = %3.2f \n", fChi2CutConversion, fPsiPairCut );
  } else if (fDo2DPsiPairChi2==2){
    printf("\t exponential psi_{pair} cut depending on chi^{2} applied with |psi_{pair}| < %3.2f*exp(%3.2f*chi^{2}) \n", fPsiPairCut, fChi2CutConversionExpFunc );
  } else {
    printf("\t chi^{2} max cut chi^{2} < %3.2f \n", fChi2CutConversion );
    printf("\t psi_{pair} max cut |psi_{pair}| < %3.2f \n", fPsiPairCut );
  }
  printf("\t %3.2f < R_{conv} < %3.2f\n", fMinR, fMaxR );
  printf("\t Z_{conv} < %3.2f\n", fMaxZ );
  if (fEtaCutMin > -0.1) printf("\t %3.2f < eta_{conv} < %3.2f\n", fEtaCutMin, fEtaCut );
    else printf("\t eta_{conv} < %3.2f\n", fEtaCut );
  if (fDoPhotonAsymmetryCut) printf("\t for p_{T,track} > %3.2f,  A_{gamma} < %3.2f \n", fMinPPhotonAsymmetryCut, fMinPhotonAsymmetry  );
  if (fDoPhotonPDependentAsymCut && fDoPhotonAsymmetryCut) printf("\t p-dependent asymmetry cut \n");
  if (fUseCorrectedTPCClsInfo) printf("\t #cluster TPC/ #findable clusters TPC (corrected for radius) > %3.2f\n", fMinClsTPCToF );
  if(fMaxTPCChi2NDF>0) printf("\t TPC Chi2 < %3.2f \n", fMaxTPCChi2NDF);
  printf("\t p_{T,gamma} > %3.2f\n", fPtCut );
  printf("\t cos(Theta_{point}) > %3.2f \n", fCosPAngleCut );
  printf("\t dca_{R} < %3.2f \n", fDCARPrimVtxCut );
  printf("\t dca_{Z} < %3.2f \n", fDCAZPrimVtxCut );
  if (fDoPhotonQualitySelectionCut) printf("\t selection based on photon quality with quality %d \n", fPhotonQualityCut );
  if (fDoPhotonQualityRejectionCut) printf("\t rejection based on photon quality with quality %d \n", fPhotonQualityCut );
  if (fPhotonQualityCutTRD || fPhotonQualityCutTOF) printf("\t TRD quality: %d, TOF quality: %d\n", fPhotonQualityCutTRD, fPhotonQualityCutTOF);
  if (fDoDoubleCountingCut) printf("\t Reject doubly counted photons with R > %3.2f, DeltaR < %3.2f, OpenAngle < %3.2f  \n", fMinRDC, fDeltaR,fOpenAngle );

}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetV0Finder(Int_t v0FinderType){   // Set Cut
  switch (v0FinderType){
  case 0:  // on fly V0 finder
    cout << "have chosen onfly V0" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=0;
    break;
  case 1:  // offline V0 finder
    cout << "have chosen offline V0" << endl;
    fUseOnFlyV0Finder=kFALSE;
    fUseOnFlyV0FinderSameSign=0;
    break;
  case 2:  // on fly V0 finder with same signs
    cout << "have chosen onfly V0 same sign pairing" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=1;
    break;
  case 3:  // on fly V0 finder with unlike signs and same signs
    cout << "have chosen onfly V0 unlike sign and same signs pairing" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=2;
    break;
  case 4:  // on fly V0 finder using BDT cuts on the rest of the variables
    cout << "have chosen onfly V0" << endl;
    fUseOnFlyV0Finder=kTRUE;
    fUseOnFlyV0FinderSameSign=0;
    fUseBDTPhotonCuts=kTRUE;
    break;
  default:
    AliError(Form(" v0FinderType not defined %d",v0FinderType));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetEtaCut(Int_t etaCut){   // Set Cut

  //Set Standard LineCutZValues
  fLineCutZValueMin = -2;
  fLineCutZValue = 7.;

  switch(etaCut){
  case 0: // 0.9
    fEtaCut     = 0.9;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 1:  // 0.6  // changed from 1.2 to 0.6 on 2013.06.10
    fEtaCut     = 0.6;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 2:  // 1.4
    fEtaCut     = 1.4;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 3: // 0.65
    fEtaCut     = 0.65;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 4: // 0.75
    fEtaCut     = 0.75;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 5: // 0.5
    fEtaCut     = 0.5;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 6: // 5.
    fEtaCut     = 5.;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 7:
    if (fIsHeavyIon==1){
      fEtaCut     = 0.7;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
    } else {
      fEtaCut     = 0.3;
      fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
      fEtaCutMin     = -0.1;
      fLineCutZRSlopeMin = 0.;
      break;
    }
  // case 8: // 0.1 - 0.8
  //    fEtaCut     = 0.9;
  //    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
  //    fEtaCutMin     = 0.1;
  //    fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
  //    break;
  case 8: // 0.4
    fEtaCut     = 0.4;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 9: // 10
    fEtaCut     = 10;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 10: // a - 0.2-0.9
    fEtaCut     = 0.9;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = 0.2;
    fLineCutZRSlopeMin = 0.;
    break;
  case 11: // b - 0.2-0.9
    fEtaCut     = 0.9;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = 0.2;
    fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin)));
    break;
  case 12: // c - 0.85
    fEtaCut     = 0.85;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  case 13: // d - 0.8
    fEtaCut     = 0.8;
    fLineCutZRSlope = tan(2*atan(exp(-fEtaCut)));
    fEtaCutMin     = -0.1;
    fLineCutZRSlopeMin = 0.;
    break;
  default:
    AliError(Form(" EtaCut not defined %d",etaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetRCut(Int_t RCut){
  // Set Cut
  switch(RCut){
  case 0:
    fMinR=0;
    fMaxR = 180.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 1:
    fMinR=2.8;
    fMaxR = 180.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 2:
    fMinR=5.;
    fMaxR = 180.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 3:
    fMaxR = 70.;
    fMinR = 10.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 4:
    fMaxR = 70.;
    fMinR = 5.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 5:
    fMaxR = 180.;
    fMinR = 10.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 6:
    fMaxR = 180.;
    fMinR = 20.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 7:
    fMaxR = 180.;
    fMinR = 35.; //old 26.
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 8:
    fMaxR = 180.;
    fMinR = 12.5;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 9:
    fMaxR = 180.;
    fMinR = 7.5;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 10:  //a
    fMaxR = 33.5;
    fMinR = 5.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 11:  //b
    fMaxR = 72.;
    fMinR = 33.5;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 12: //c
    fMaxR = 180.;
    fMinR = 72.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 13: //d
    fMaxR = 55.;
    fMinR = 5.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 14: //e
    fMaxR = 180.;
    fMinR = 55.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 15: //f
    fMaxR = 72.;
    fMinR = 5.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 16: //g
    fMaxR = 180.;
    fMinR = 95.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 17: //h
    fMaxR = 13.;
    fMinR = 5.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 18: //i
    fMaxR = 33.5;
    fMinR = 13.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 19:  // j
    fMaxR = 55.;
    fMinR = 33.5;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 20: //k
    fMaxR = 72.;
    fMinR = 55.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 21: //l
    fMaxR = 95.;
    fMinR = 72.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 22: //m
    fMaxR = 180.;
    fMinR = 5.;
    fExcludeMinR = 55.;
    fExcludeMaxR = 72.;
    break;
  case 23: //n
    fMaxR = 180.;
    fMinR = 10.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 24: //o
    fMaxR = 180.;
    fMinR = 15.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 25: //p
    fMaxR = 180.;
    fMinR = 20.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 26: //q
    fMaxR = 95.;
    fMinR = 5.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 27: //r
    fMaxR = 95.;
    fMinR = 10.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 28: //s
    fMaxR = 95.;
    fMinR = 15.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;
  case 29: //t
    fMaxR = 95.;
    fMinR = 20.;
    fExcludeMinR = 180.;
    fExcludeMaxR = 250.;
    break;

  default:
    AliError("RCut not defined");
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetEtaForPhiCut(Int_t etaPhiCut) {

  switch(etaPhiCut) {
  case 0: //no specific eta range selected, full eta range
    fEtaForPhiCutMin = -fEtaCut;
    fEtaForPhiCutMax = fEtaCut;
    break;
  case 1:  //eta < 0 only
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fEtaForPhiCutMin = -fEtaCut;
    fEtaForPhiCutMax = 0.;
    break;
  case 2://eta > 0 only
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fEtaForPhiCutMin = 0.;
    fEtaForPhiCutMax = fEtaCut;
    break;
  case 3:  // distortions cut on A and C side
      fDoShrinkTPCAcceptance = 2; // Only use photons in phi regions without distortions
      break;
  case 4:  // distortions cut on A and C side
      fDoShrinkTPCAcceptance = 3; // Only use photons in phi regions with strong distortions
      break;
  default:
    AliError(Form("EtaForPhiCut not defined %d",etaPhiCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
// This are exclusion cuts, everything between fMinPhiCut & fMaxPhiCut will be excluded
Bool_t AliConversionPhotonCuts::SetMinPhiSectorCut(Int_t minPhiCut) {

  switch(minPhiCut) {
  case 0:
    fDoShrinkTPCAcceptance = 0;
    fMinPhiCut = 0;
    break;
  case 1:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMinPhiCut = 1.7; //OROC C08 large cut
    break;
  case 2:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMinPhiCut = 4.4; //EMCal
    break;
  case 3:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMinPhiCut = 1.0; //PHOS
    break;
  case 4:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMinPhiCut = 3.4; //EMCal tight
    break;
  case 5:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMinPhiCut = 2.0; //OROC C08 medium cut
    break;
  case 6:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMinPhiCut = 2.2; //OROC C08 small cut
    break;
  case 7:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMinPhiCut = 2.4; //OROC C08 tightest cut
    break;
  case 8:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMinPhiCut = 4.54; //PHOS phi
    break;
  case 9:  // distortions cut on A and C side
      fDoShrinkTPCAcceptance = 2; // Only use photons in phi regions without distortions
      break;
  case 10: // distortions cut on A and C side
      fDoShrinkTPCAcceptance = 3; // Only use photons in phi regions with strong distortions
      break;
  case 11:  // b
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMinPhiCut = 0.; // to calculate MBW for PHOS pi0 region
    break;

  default:
    AliError(Form("MinPhiCut not defined %d",minPhiCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
// This are exclusion cuts, everything between fMinPhiCut & fMaxPhiCut will be excluded
Bool_t AliConversionPhotonCuts::SetMaxPhiSectorCut(Int_t maxPhiCut) {

  switch(maxPhiCut) {
  case 0:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 0;
    fMaxPhiCut = 2*TMath::Pi()+0.00001;
    break;
  case 1:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMaxPhiCut = 4.3; //OROC C08 large cut
    break;
  case 2:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMaxPhiCut = 5.8; //EMCal
    break;
  case 3:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMaxPhiCut = 3.0; //PHOS
    break;
  case 4:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMaxPhiCut = 1.; //EMCal
    break;
  case 5:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMaxPhiCut = 4.0; //OROC C08 medium cut
    break;
  case 6:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMaxPhiCut = 3.8; //OROC C08 small cut
    break;
  case 7:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMaxPhiCut = 3.6; //OROC C08 tighest cut
    break;
  case 8:
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMaxPhiCut = 5.59; //PHOS phi
    break;
  case 9:   // distortions cut on A and C side
      fDoShrinkTPCAcceptance = 2; // Only use photons in phi regions without distortions
      break;
  case 10: // distortions cut on A and C side
      fDoShrinkTPCAcceptance = 3; // Only use photons in phi regions with strong distortions
      break;
  case 11:  // b
    if (!fDoShrinkTPCAcceptance) fDoShrinkTPCAcceptance = 1;
    fMaxPhiCut = 3.16; // to calculate MBW for PHOS pi0 region
    break;

  default:
    AliError(Form("MaxPhiCut not defined %d",maxPhiCut));
    return kFALSE;
  }

  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetSinglePtCut(Int_t singlePtCut){   // Set Cut
  switch(singlePtCut){
  case 0: // 0.050 GeV + min gamma pT cut of 20 MeV
    fSinglePtCut = 0.050;
    fPtCut       = 0.02;
    break;
  case 1:  // 0.100 GeV  + min gamma pT cut of 20 MeV
    fSinglePtCut = 0.100;
    fPtCut       = 0.02;
    break;
  case 2:  // 0.150 GeV  + min gamma pT cut of 20 MeV
    fSinglePtCut = 0.150;
    fPtCut       = 0.02;
    break;
  case 3:  // 0.200 GeV  + min gamma pT cut of 20 MeV
    fSinglePtCut = 0.200;
    fPtCut       = 0.02;
    break;
  case 4:  // 0.075 GeV  + min gamma pT cut of 20 MeV
    fSinglePtCut = 0.075;
    fPtCut       = 0.02;
    break;
  case 5:  // 0.125 GeV  + min gamma pT cut of 20 MeV
    fSinglePtCut = 0.125;
    fPtCut       = 0.02;
    break;
  case 6:  // 0.04 GeV  + min gamma pT cut of 10 MeV
    fSinglePtCut = 0.040;
    fPtCut       = 0.01;
    break;
  case 7:  // 0.0 GeV  + min gamma pT cut of 0 MeV
    fSinglePtCut = 0.0;
    fPtCut       = 0.0;
    break;
  case 8:  // 0.02 GeV + min gamma pT cut of 20 MeV,  equivalent to .05 for the low B field runs
    fSinglePtCut = 0.02;
    fPtCut       = 0.01;
    break;
  case 9:  // 0.050 GeV + min gamma pT cut of 100 MeV
    fSinglePtCut = 0.050;
    fPtCut       = 0.100;
    break;
  case 10:  //a: 0.050 GeV + min gamma pT cut of 150 MeV
    fSinglePtCut = 0.050;
    fPtCut       = 0.150;
    break;
  case 11:  //b: 0.050 GeV + min gamma pT cut of 200 MeV
    fSinglePtCut = 0.050;
    fPtCut       = 0.200;
    break;
  case 12:  //c: 0.060 GeV
    fSinglePtCut = 0.060;
    break;
  case 13:  //d: 0.060 GeV + min gamma pT cut of 100 MeV
    fSinglePtCut = 0.060;
    fPtCut       = 0.100;
    break;
  case 14:  //e: 0.060 GeV + min gamma pT cut of 150 MeV
    fSinglePtCut = 0.060;
    fPtCut       = 0.150;
    break;
  case 15:  //f: 0.060 GeV + min gamma pT cut of 200 MeV
    fSinglePtCut = 0.060;
    fPtCut       = 0.200;
    break;
  case 16:  //g: 0.075 GeV + min gamma pT cut of 150 MeV
    fSinglePtCut = 0.075;
    fPtCut       = 0.150;
    break;
  case 17:  //h: 0.100 GeV + min gamma pT cut of 200 MeV
    fSinglePtCut = 0.100;
    fPtCut       = 0.200;
    break;
  case 18:  //i: 0.150 GeV + min gamma pT cut of 300 MeV
    fSinglePtCut = 0.150;
    fPtCut       = 0.300;
    break;
  case 19:  //j: asym: 0.100 GeV and 0.075 GeV
    fSinglePtCut = 0.100;
    fDoAsymPtCut = kTRUE;
    fSinglePtCut2= 0.075;
    break;
  case 20:  //k: asym: 0.150 GeV and 0.075 GeV
    fSinglePtCut = 0.150;
    fDoAsymPtCut = kTRUE;
    fSinglePtCut2= 0.075;
    break;
  case 21:  //l: asym: 0.200 GeV and 0.075 GeV
    fSinglePtCut = 0.200;
    fDoAsymPtCut = kTRUE;
    fSinglePtCut2= 0.075;
    break;
  case 22: // m: 0.080 GeV + min gamma pT cut of 20 MeV
    fSinglePtCut = 0.080;
    fPtCut       = 0.02;
    break;
  case 23: // n: 0.090 GeV + min gamma pT cut of 20 MeV
    fSinglePtCut = 0.090;
    fPtCut       = 0.02;
    break;
  case 24: // o: 0.024 GeV + min gamma pT cut of 20 MeV ; equiv. 0.06 for lowB
    fSinglePtCut = 0.024;
    fPtCut       = 0.01;
    break;
  case 25: // p: 0.030 GeV + min gamma pT cut of 20 MeV ; equiv. 0.075 for lowB
    fSinglePtCut = 0.030;
    fPtCut       = 0.01;
    break;
  case 26: // q: 0.032 GeV + min gamma pT cut of 20 MeV ; equiv. 0.08 for lowB
    fSinglePtCut = 0.032;
    fPtCut       = 0.01;
    break;
  case 27: // r: 0.036 GeV + min gamma pT cut of 20 MeV ; equiv. 0.09 for lowB
    fSinglePtCut = 0.036;
    fPtCut       = 0.01;
    break;
  case 28: // s: 0.040 GeV + min gamma pT cut of 20 MeV ; equiv. 0.075 for lowB
    fSinglePtCut = 0.040;
    fPtCut       = 0.01;
    break;

  default:
    AliError(Form("singlePtCut not defined %d",singlePtCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTPCClusterCut(Int_t clsTPCCut){   // Set Cut
  switch(clsTPCCut){
  case 0: // 0
    fMinClsTPC= 0.;
    break;
  case 1:  // 60
    fMinClsTPC= 60.;
    break;
  case 2:  // 80
    fMinClsTPC= 80.;
    break;
  case 3:  // 100
    fMinClsTPC= 100.;
    break;
  case 4:  // 95% of findable clusters
    fMinClsTPCToF= 0.95;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 5:  // 0% of findable clusters
    fMinClsTPCToF= 0.0;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 6:  // 70% of findable clusters
    fMinClsTPCToF= 0.7;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 7:  // 0% of findable clusters
    fMinClsTPCToF= 0.35;
    fUseCorrectedTPCClsInfo=0;
    break;
  case 8:  // 35% of findable clusters
      fMinClsTPCToF= 0.35;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 9:  // 60% of findable clusters
    fMinClsTPCToF= 0.6;
    fUseCorrectedTPCClsInfo=1;
    break;
  case 10:  // 60% of findable clusters and TPC track Chi2<4
      fMinClsTPCToF= 0.6;
      fUseCorrectedTPCClsInfo=1;
      fMaxTPCChi2NDF=4.0;
      break;
  case 11:  // 60% of findable clusters and TPC track Chi2<3
      fMinClsTPCToF= 0.6;
      fUseCorrectedTPCClsInfo=1;
      fMaxTPCChi2NDF=3.0;
      break;
  case 12:  // 60% of findable clusters and TPC track Chi2<2.5
      fMinClsTPCToF= 0.6;
      fUseCorrectedTPCClsInfo=1;
      fMaxTPCChi2NDF=2.5;
      break;
  default:
    AliError(Form("Warning: clsTPCCut not defined %d",clsTPCCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut){   // Set Cut
  switch(ededxSigmaCut){
  case 0: // -10,10
    fPIDnSigmaBelowElectronLine=-10;
    fPIDnSigmaAboveElectronLine=10;
    break;
  case 1: // -5,5
    fPIDnSigmaBelowElectronLine=-5;
    fPIDnSigmaAboveElectronLine=5;
    break;
  case 2: // -3,5
    fPIDnSigmaBelowElectronLine=-3;
    fPIDnSigmaAboveElectronLine=5;
    break;
  case 3: // -4,5
    fPIDnSigmaBelowElectronLine=-4;
    fPIDnSigmaAboveElectronLine=5;
    break;
  case 4: // -6,7
    fPIDnSigmaBelowElectronLine=-6;
    fPIDnSigmaAboveElectronLine=7;
    break;
  case 5: // -4,4
    fPIDnSigmaBelowElectronLine=-4;
    fPIDnSigmaAboveElectronLine=4;
    break;
  case 6: // -2.5,4
    fPIDnSigmaBelowElectronLine=-2.5;
    fPIDnSigmaAboveElectronLine=4;
    break;
  case 7: // -2,3.5
    fPIDnSigmaBelowElectronLine=-2;
    fPIDnSigmaAboveElectronLine=3.5;
    break;
  case 8: // -2.5,3.
    fPIDnSigmaBelowElectronLine=-2.5;
    fPIDnSigmaAboveElectronLine=3;
    break;
  case 9: // -2.5,5.
    fPIDnSigmaBelowElectronLine=-2.5;
    fPIDnSigmaAboveElectronLine=5;
    break;
  case 10: //a -3,3.
    fPIDnSigmaBelowElectronLine=-3;
    fPIDnSigmaAboveElectronLine=3;
    break;
  case 11: //b -3.2,3.2.
    fPIDnSigmaBelowElectronLine=-3.2;
    fPIDnSigmaAboveElectronLine=3.2;
    break;
  case 12: //c -2.8,2.8
    fPIDnSigmaBelowElectronLine=-2.8;
    fPIDnSigmaAboveElectronLine=2.8;
    break;
  case 13: //d -1E9,1E9
    fPIDnSigmaBelowElectronLine=-1E9;
    fPIDnSigmaAboveElectronLine=1E9;
    break;
  case 14: //e -7,1
    fPIDnSigmaBelowElectronLine=-7;
    fPIDnSigmaAboveElectronLine=1;
    break;
  case 15: //f -3,4
    fPIDnSigmaBelowElectronLine=-3;
    fPIDnSigmaAboveElectronLine=4;
    break;
  default:
    AliError("TPCdEdxCutElectronLine not defined");
    return kFALSE;

  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut){   // Set Cut

  switch(pidedxSigmaCut){
  case 0:  // -10
    fPIDnSigmaAbovePionLine=-10;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 1:   // 0
    fPIDnSigmaAbovePionLine=0;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 2:  // 1
    fPIDnSigmaAbovePionLine=1;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 3:  // 1
    fPIDnSigmaAbovePionLine=2.5;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 4:  // 3.0sigma, 1.0 sigma at high pt
    fPIDnSigmaAbovePionLine=3.0;
    fPIDnSigmaAbovePionLineHighPt=1.;
    break;
  case 5:  // 1
    fPIDnSigmaAbovePionLine=2.;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 6:  // 1
    fPIDnSigmaAbovePionLine=2.;
    fPIDnSigmaAbovePionLineHighPt=0.5;
    break;
  case 7:  // 1
    fPIDnSigmaAbovePionLine=3.5;
    fPIDnSigmaAbovePionLineHighPt=-10;
    break;
  case 8:  // 1
    fPIDnSigmaAbovePionLine=2.;
    fPIDnSigmaAbovePionLineHighPt=1.;
    break;
  case 9:
    fPIDnSigmaAbovePionLine=1; // We need a bit less tight cut on dE/dx
    fPIDnSigmaAbovePionLineHighPt=0.5;
    break;
  case 10: //a
    fPIDnSigmaAbovePionLine=-3; // We need a bit less tight cut on dE/dx
    fPIDnSigmaAbovePionLineHighPt=-14;
  case 11: //b
    fPIDnSigmaAbovePionLine=3;
    fPIDnSigmaAbovePionLineHighPt=2;
    break;
  default:
    AliError(Form("Warning: pidedxSigmaCut not defined %d",pidedxSigmaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetMinMomPiondEdxCut(Int_t piMomdedxSigmaCut){   // Set Cut
  switch(piMomdedxSigmaCut){
  case 0:  // 0.5 GeV
    fPIDMinPnSigmaAbovePionLine=0.5;
    break;
  case 1:  // 1. GeV
    fPIDMinPnSigmaAbovePionLine=1.;
    break;
  case 2:  // 1.5 GeV
    fPIDMinPnSigmaAbovePionLine=1.5;
    break;
  case 3:  // 20.0 GeV
    fPIDMinPnSigmaAbovePionLine=20.;
    break;
  case 4:  // 50.0 GeV
    fPIDMinPnSigmaAbovePionLine=50.;
    break;
  case 5:  // 0.3 GeV
    fPIDMinPnSigmaAbovePionLine=0.3;
    break;
  case 6:  // 0.25 GeV
    fPIDMinPnSigmaAbovePionLine=0.25;
    break;
  case 7:  // 0.4 GeV
    fPIDMinPnSigmaAbovePionLine=0.4;
    break;
  case 8:  // 0.2 GeV
    fPIDMinPnSigmaAbovePionLine=0.2;
    break;
  default:
    AliError(Form("piMomdedxSigmaCut not defined %d",piMomdedxSigmaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetMaxMomPiondEdxCut(Int_t piMaxMomdedxSigmaCut){   // Set Cut
  switch(piMaxMomdedxSigmaCut){
  case 0:  // 100. GeV
    fPIDMaxPnSigmaAbovePionLine=100.;
    break;
  case 1:  // 5. GeV
    fPIDMaxPnSigmaAbovePionLine=5.;
    break;
  case 2:  // 4. GeV
    fPIDMaxPnSigmaAbovePionLine=4.;
    break;
  case 3:  // 3.5 GeV
    fPIDMaxPnSigmaAbovePionLine=3.5;
    break;
  case 4:  // 3. GeV
    fPIDMaxPnSigmaAbovePionLine=3.;
    break;
  case 5:  // 7. GeV
    fPIDMaxPnSigmaAbovePionLine=7.;
    break;
  case 6:  // 2. GeV
    fPIDMaxPnSigmaAbovePionLine=2.;
    break;
  case 7:  // 8. GeV
    fPIDMaxPnSigmaAbovePionLine=8.;
    break;
  default:
    AliError(Form("piMaxMomdedxSigmaCut not defined %d",piMaxMomdedxSigmaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut){   // Set Cut
  switch(LowPRejectionSigmaCut){
  case 0:  //
    fPIDnSigmaAtLowPAroundKaonLine=0;
    fPIDnSigmaAtLowPAroundProtonLine=0;
    fPIDnSigmaAtLowPAroundPionLine=0;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kFALSE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 1:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.5;
    fPIDnSigmaAtLowPAroundProtonLine=0.5;
    fPIDnSigmaAtLowPAroundPionLine=0.5;
    fDoKaonRejectionLowP = kTRUE;
    fDoProtonRejectionLowP = kTRUE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 2:  //
    fPIDnSigmaAtLowPAroundKaonLine=1;
    fPIDnSigmaAtLowPAroundProtonLine=1;
    fPIDnSigmaAtLowPAroundPionLine=1;
    fDoKaonRejectionLowP = kTRUE;
    fDoProtonRejectionLowP = kTRUE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 3:  //
    fPIDnSigmaAtLowPAroundKaonLine=2.;
    fPIDnSigmaAtLowPAroundProtonLine=2.;
    fPIDnSigmaAtLowPAroundPionLine=2.;
    fDoKaonRejectionLowP = kTRUE;
    fDoProtonRejectionLowP = kTRUE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 4:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.;
    fPIDnSigmaAtLowPAroundPionLine=1;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 5:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.;
    fPIDnSigmaAtLowPAroundPionLine=1.5;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 6:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.;
    fPIDnSigmaAtLowPAroundPionLine=2.;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 7:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.;
    fPIDnSigmaAtLowPAroundPionLine=0.5;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kFALSE;
    fDoPionRejectionLowP = kTRUE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  case 8:  //
    fPIDnSigmaAtLowPAroundKaonLine=0.;
    fPIDnSigmaAtLowPAroundProtonLine=0.5;
    fPIDnSigmaAtLowPAroundPionLine=0.;
    fDoKaonRejectionLowP = kFALSE;
    fDoProtonRejectionLowP = kTRUE;
    fDoPionRejectionLowP = kFALSE;
    fPIDMinPPionRejectionLowP = fPIDMinPnSigmaAbovePionLine;
    break;
  default:
    AliError(Form("LowPRejectionSigmaCut not defined %d",LowPRejectionSigmaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetKappaTPCCut(Int_t kappaCut){   // Set Cut
  switch(kappaCut){
  case 0: // completely open
    fKappaMaxCut=200;
    fKappaMinCut=-200;
    break;
  case 1: // mainly pi pi
    fKappaMaxCut=-13;
    fKappaMinCut=-20;
    break;
  case 2: // mainly pi e
    fKappaMaxCut=-6;
    fKappaMinCut=-11;
    break;
  case 3: // signal
    fKappaMaxCut=5;
    fKappaMinCut=-3;
    break;
  case 4: // remaining
    fKappaMaxCut=20;
    fKappaMinCut=11;
    break;
  case 5: // -5-10 full signal peak(including background)
    fKappaMaxCut=10;
    fKappaMinCut=-5;
    break;
  case 6: //
    fKappaMaxCut=10;
    fKappaMinCut=-3;
    break;
  case 7: //
    fKappaMaxCut=10;
    fKappaMinCut=0;
    break;
  default:
    AliError("KappaTPCCut not defined");
    return kFALSE;

  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTOFElectronPIDCut(Int_t TOFelectronPID){
  // Set Cut
  switch(TOFelectronPID){
  case 0: // no cut
    fUseTOFpid = kFALSE;
    fTofPIDnSigmaBelowElectronLine=-100;
    fTofPIDnSigmaAboveElectronLine=100;
    break;
  case 1: // -7,7
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-7;
    fTofPIDnSigmaAboveElectronLine=7;
    break;
  case 2: // -5,5
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-5;
    fTofPIDnSigmaAboveElectronLine=5;
    break;
  case 3: // -3,5
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-3;
    fTofPIDnSigmaAboveElectronLine=5;
    break;
  case 4: // -2,3
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-2;
    fTofPIDnSigmaAboveElectronLine=3;
    break;
  case 5: // -3,3
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-3;
    fTofPIDnSigmaAboveElectronLine=3;
    break;
  case 6: // TOF timing one leg
    fUseTOFpid = kFALSE;
    fUseTOFtiming = kTRUE;
    fTOFtimingBothLegs = kFALSE;
    break;
  case 7: // TOF timing both legs
    fUseTOFpid = kFALSE;
    fUseTOFtiming = kTRUE;
    fTOFtimingBothLegs = kTRUE;
    break;
  case 8: // TOF timing one leg and within 100ns
    fUseTOFpid = kFALSE;
    fUseTOFtiming = kTRUE;
    fTOFtimeMin = -100;
    fTOFtimeMax = 100;
    fTOFtimingBothLegs = kFALSE;
    break;
  case 9: // TOF timing both legs and within 100ns
    fUseTOFpid = kFALSE;
    fUseTOFtiming = kTRUE;
    fTOFtimeMin = -100;
    fTOFtimeMax = 100;
    fTOFtimingBothLegs = kTRUE;
    break;
  case 10: // a  -10,6
    fUseTOFpid = kTRUE;
    fTofPIDnSigmaBelowElectronLine=-10;
    fTofPIDnSigmaAboveElectronLine=6;
    break;
  case 11: // b -4,4 but only if the track momenta are above 0.4GeV/c to cope with large TOF mismatch in central AA collisions at low pT
      fUseTOFpid = kTRUE;
      fTofPIDnSigmaBelowElectronLine=-4;
      fTofPIDnSigmaAboveElectronLine=4;
      fUseTOFpidMinMom = kTRUE;
      fTofPIDMinMom = 0.4;
  default:
    AliError(Form("TOFElectronCut not defined %d",TOFelectronPID));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetITSElectronPIDCut(Int_t ITSelectronPID){
  // Set Cut
  switch(ITSelectronPID){
  case 0: // no cut
    fUseITSpid = kFALSE;
    fITSPIDnSigmaBelowElectronLine=-100;
    fITSPIDnSigmaAboveElectronLine=100;
    fMaxPtPIDITS = 1.5;
    break;
  case 1: // -3,3
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-3;
    fITSPIDnSigmaAboveElectronLine=3;
    fMaxPtPIDITS = 1.5;
    break;
  case 2: // -2,2
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-2;
    fITSPIDnSigmaAboveElectronLine=2;
    fMaxPtPIDITS = 1.5;
    break;
  case 3: // -1,1
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-1;
    fITSPIDnSigmaAboveElectronLine=1;
    fMaxPtPIDITS = 1.5;
    break;
  case 4: // -3,5
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-3;
    fITSPIDnSigmaAboveElectronLine=5;
    fMaxPtPIDITS = 1.5;
    break;
  case 5: // -5,5
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-5;
    fITSPIDnSigmaAboveElectronLine=5;
    fMaxPtPIDITS = 1.5;
    break;
  case 6: // -3,3
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-3;
    fITSPIDnSigmaAboveElectronLine=3;
    fMaxPtPIDITS = 2;
    break;
  case 7: // -2,2
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-2;
    fITSPIDnSigmaAboveElectronLine=2;
    fMaxPtPIDITS = 2;
    break;
  case 8: // -1,1
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-1;
    fITSPIDnSigmaAboveElectronLine=1;
    fMaxPtPIDITS = 2;
    break;
  case 9: // -3,5
    fUseITSpid = kTRUE;
    fITSPIDnSigmaBelowElectronLine=-3;
    fITSPIDnSigmaAboveElectronLine=5;
    fMaxPtPIDITS = 2;
    break;
  default:
    AliError(Form("ITSelectronPID not defined %d",ITSelectronPID));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTRDElectronPIDCut(Int_t TRDelectronPID){
  // Set Cut
  switch(TRDelectronPID){
  case 0: // no cut
    fDoTRDPID = kFALSE;
    fTRDPIDBelowCut=-100;
    fTRDPIDAboveCut=100;
    break;
  default:
    AliError(Form("TRDelectronPID not defined %d",TRDelectronPID));
    return kFALSE;
  }
  return kTRUE;
}


///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetQtMaxCut(Int_t QtMaxCut){   // Set Cut
  switch(QtMaxCut){
  case 0: //
    fQtMax=1.;
    fDoQtGammaSelection=0;
    fDo2DQt=kFALSE;
    break;
  case 1:
    fQtMax=0.1;
    fDo2DQt=kFALSE;
    break;
  case 2:
    fQtMax=0.06;
    fDo2DQt=kTRUE;
    break;
  case 3:
    fQtMax=0.05;
    fDo2DQt=kFALSE;
    break;
  case 4:
    fQtMax=0.03;
    fDo2DQt=kFALSE;
    break;
  case 5:
    fQtMax=0.02;
    fDo2DQt=kFALSE;
    break;
  case 6:
    fQtMax=0.02;
    fDo2DQt=kTRUE;
    break;
  case 7:
    fQtMax=0.15;
    fDo2DQt=kFALSE;
    break;
  case 8:
    fQtMax=0.05;
    fDo2DQt=kTRUE;
    break;
  case 9:
    fQtMax=0.03;
    fDo2DQt=kTRUE;
    break;
  case 10:  //a
    fQtPtMax=0.11;
    fQtMax=0.040;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  case 11: //b
    fQtPtMax=0.125;
    fQtMax=1;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  case 12:  //c
    fQtPtMax=0.125;
    fQtMax=1;
    fDoQtGammaSelection=2;
    fDo2DQt=kFALSE;
    break;
  case 13:  //d
    fQtPtMax=0.125;
    fQtMax=0.050;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  case 14:  //e
    fQtPtMax=0.14;
    fQtMax=0.060;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  case 15:  //f
    fQtPtMax=0.16;
    fQtMax=0.070;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  case 16:  //g
    fQtPtMax=0.180;
    fQtMax=0.032;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  case 17:  //h
    fQtPtMax=0.20;
    fQtMax=1;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  case 18:  //i
    fQtPtMax=0.20;
    fQtMax=0.035;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  case 19:  //j
    fQtPtMax=0.25;
    fQtMax=0.040;
    fDoQtGammaSelection=2;
    fDo2DQt=kFALSE;
    break;
  case 20:  //k
    fQtPtMax=0.30;
    fQtMax=0.045;
    fDoQtGammaSelection=2;
    fDo2DQt=kFALSE;
    break;
  case 21:  //l
    fQtPtMax=0.11;
    fQtMax=0.030;
    fDoQtGammaSelection=2;
    fDo2DQt=kTRUE;
    break;
  default:
    AliError(Form("Warning: QtMaxCut not defined %d",QtMaxCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetChi2GammaCut(Int_t chi2GammaCut){   // Set Cut

  switch(chi2GammaCut){
  case 0: // 100
    fChi2CutConversion = 100.;
    break;
  case 1:  // 50
    fChi2CutConversion = 50.;
    break;
  case 2:  // 30
    fChi2CutConversion = 30.;
    break;
  case 3:
    fChi2CutConversion = 200.;
    break;
  case 4:
    if (fIsHeavyIon==1){
      fChi2CutConversion = 7.;
    } else {
      fChi2CutConversion = 500.;
    }
    break;
  case 5:
    fChi2CutConversion = 100000.;
    break;
  case 6:
    fChi2CutConversion = 5.;
    break;
  case 7:
    fChi2CutConversion = 10.;
    break;
  case 8:
    fChi2CutConversion = 20.;
    break;
  case 9:
    fChi2CutConversion = 15.;
    break;
  case 10: //a
    fChi2CutConversion = 25.;
    break;
  case 11: //b
    fChi2CutConversion = 35.;
    break;
  case 12: //c
    fChi2CutConversion = 40.;
    break;
  case 13: //d
    fChi2CutConversion = 45.;
    break;
  case 14: //e
    fChi2CutConversion = 55.;
    break;
  case 15: //f for exp cut (fDo2DPsiPairChi2 = 2)
    fChi2CutConversion = 50.;
    fChi2CutConversionExpFunc = -0.065;
    break;
  case 16: //g for exp cut (fDo2DPsiPairChi2 = 2)
    fChi2CutConversion = 50.;
    fChi2CutConversionExpFunc = -0.055;
    break;
  case 17: //h for exp cut (fDo2DPsiPairChi2 = 2)
    fChi2CutConversion = 50.;
    fChi2CutConversionExpFunc = -0.050;
    break;
  case 18: //i for exp cut (fDo2DPsiPairChi2 = 2) low B
    fChi2CutConversion = 50.;
    fChi2CutConversionExpFunc = -0.075;
    break;
  case 19: //j for exp cut (fDo2DPsiPairChi2 = 2) low B
    fChi2CutConversion = 50.;
    fChi2CutConversionExpFunc = -0.085;
    break;
  case 20: //k for exp cut (fDo2DPsiPairChi2 = 2)
    fChi2CutConversion = 20.;
    fChi2CutConversionExpFunc = -0.055;
    break;
  case 21: //l for exp cut (fDo2DPsiPairChi2 = 2)
    fChi2CutConversion = 30.;
    fChi2CutConversionExpFunc = -0.11;
    break;
  default:
    AliError(Form("Warning: Chi2GammaCut not defined %d",chi2GammaCut));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetPsiPairCut(Int_t psiCut) {
  switch(psiCut) {
  case 0:
    fPsiPairCut = 10000; //
    break;
  case 1:
    fPsiPairCut = 0.1; //
    break;
  case 2:
    fPsiPairCut = 0.05; // Standard
    break;
  case 3:
    fPsiPairCut = 0.035; //
    break;
  case 4:
    fPsiPairCut = 0.2; //
    break;
  case 5:
    fPsiPairCut = 0.1; //
    fDo2DPsiPairChi2 = 1;
    break;
  case 6:
    fPsiPairCut = 0.05; //
    fDo2DPsiPairChi2 = 1;
    break;
  case 7:
    if (fIsHeavyIon==1){
      fPsiPairCut = 0.07; //
    } else {
      fPsiPairCut = 0.035; //
    }
    fDo2DPsiPairChi2 = 1;
    break;
  case 8:
    fPsiPairCut = 0.2; //
    fDo2DPsiPairChi2 = 1; //
    break;
  case 9:
    //   if (fIsHeavyIon==1){ //AM 2016-05-13
      fPsiPairCut = 0.1; //
      fDo2DPsiPairChi2 = 1;
      fIncludeRejectedPsiPair = kTRUE;
      break;
  case 10: //a
    fPsiPairCut = 0.25; //
    fDo2DPsiPairChi2 = 1; //
    break;
  case 11: //b
    fPsiPairCut = 0.3; //
    fDo2DPsiPairChi2 = 1; //
    break;
  case 12: //c
    fPsiPairCut = 0.15; //
    fDo2DPsiPairChi2 = 1; //
    break;
  case 13: //d
    fPsiPairCut = 0.15; //
    fDo2DPsiPairChi2 = 2; //
    break;
  case 14: //e
    fPsiPairCut = 0.18; //
    fDo2DPsiPairChi2 = 2; //
    break;
  case 15: //f
    fPsiPairCut = 0.20; //
    fDo2DPsiPairChi2 = 2; //
    break;
  case 16: //g
    fPsiPairCut = 0.3; //
    fDo2DPsiPairChi2 = 2; //
    break;
  case 17: //h
    fPsiPairCut = 0.35; //
    fDo2DPsiPairChi2 = 2; //
    break;
  case 18: //i
    fPsiPairCut = 0.40; //
    fDo2DPsiPairChi2 = 2; //
    break;
  default:
    AliError(Form("PsiPairCut not defined %d",psiCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetPhotonAsymmetryCut(Int_t doPhotonAsymmetryCut){
  // Set Cut
  switch(doPhotonAsymmetryCut){
  case 0:
    fDoPhotonAsymmetryCut=0;
    fMinPPhotonAsymmetryCut=100.;
    fMinPhotonAsymmetry=0.;
    break;
  case 1:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=3.5;
    fMinPhotonAsymmetry=0.04;
    break;
  case 2:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=3.5;
    fMinPhotonAsymmetry=0.06;
    break;
  case 3:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=0.0;
    fMinPhotonAsymmetry=0.05;
    break;
  case 4:
    fDoPhotonAsymmetryCut=1;
    fDoPhotonPDependentAsymCut=1;
    fFAsymmetryCut = new TF1("fFAsymmetryCut","[0] + [1]*tanh(2*TMath::Power(x,[2]))",0.,100.);
    fFAsymmetryCut->SetParameter(0,0.3);
    fFAsymmetryCut->SetParameter(1,0.66);
    fFAsymmetryCut->SetParameter(2,0.7);
    fMinPPhotonAsymmetryCut=0.0;
    fMinPhotonAsymmetry=0.;
    break;
  case 5:
    fDoPhotonAsymmetryCut=1;
    fDoPhotonPDependentAsymCut=1;
    fFAsymmetryCut = new TF1("fFAsymmetryCut","[0] + [1]*tanh(2*TMath::Power(x,[2]))",0.,100.);
    fFAsymmetryCut->SetParameter(0,0.14);
    fFAsymmetryCut->SetParameter(1,0.66);
    fFAsymmetryCut->SetParameter(2,0.5);
    fMinPPhotonAsymmetryCut=0.0;
    fMinPhotonAsymmetry=0.;
    break;
  case 6:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=6.;
    fMinPhotonAsymmetry=0.05;
    break;
  case 7:
    fDoPhotonAsymmetryCut=1;
    fMinPPhotonAsymmetryCut=8.;
    fMinPhotonAsymmetry=0.05;
    break;
  case 8:
    fDoPhotonAsymmetryCut=1;
    fDoPhotonPDependentAsymCut=1;
    fFAsymmetryCut = new TF1("fFAsymmetryCut","[0] + [1]*tanh(2*TMath::Power(x,[2]))",0.,100.);
    fFAsymmetryCut->SetParameter(0,0.5);
    fFAsymmetryCut->SetParameter(1,0.46);
    fFAsymmetryCut->SetParameter(2,0.7);
    fMinPPhotonAsymmetryCut=0.0;
    fMinPhotonAsymmetry=0.;
    break;
  case 9:
    fDoPhotonAsymmetryCut=0;
    fMinPPhotonAsymmetryCut=100.;
    fMinPhotonAsymmetry=0.;
    fMaxPhotonAsymmetry=0.99;
    break;
  case 10:
    fDoPhotonAsymmetryCut=0;
    fMinPPhotonAsymmetryCut=100.;
    fMinPhotonAsymmetry=0.;
    fMaxPhotonAsymmetry=1.0;
    break;
  default:
    AliError(Form("PhotonAsymmetryCut not defined %d",doPhotonAsymmetryCut));
    return kFALSE;
  }
  fCuts[kdoPhotonAsymmetryCut]=doPhotonAsymmetryCut;
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetCosPAngleCut(Int_t cosCut) {

  switch(cosCut){
  case 0:
    fCosPAngleCut = -1;
    break;
  case 1:
    fCosPAngleCut = 0;
    break;
  case 2:
    fCosPAngleCut = 0.5;
    break;
  case 3:
    fCosPAngleCut = 0.75;
    break;
  case 4:
    fCosPAngleCut = 0.85;
    break;
  case 5:
    fCosPAngleCut = 0.88;
    break;
  case 6:
    fCosPAngleCut = 0.9;
    break;
  case 7:
    fCosPAngleCut = 0.95;
    break;
  case 8:
    fCosPAngleCut = 0.98;
    break;
  case 9:
    fCosPAngleCut = 0.99;
    break;
  case 10://a
    fCosPAngleCut = 0.995;
    break;
  case 11://b
    fCosPAngleCut = 0.985;
    break;
  case 12://c
    fCosPAngleCut = 0.996;
    break;
  case 13://d
    fCosPAngleCut = 0.997;
    break;
  case 14://e
    fCosPAngleCut = 0.998;
    break;
  case 15://f
    fCosPAngleCut = 0.999;
    break;
  default:
    AliError(Form("Cosine Pointing Angle cut not defined %d",cosCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetSharedElectronCut(Int_t sharedElec) {

    switch(sharedElec){
    case 0:
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kFALSE;
      fDoPhotonQualityRejectionCut = kFALSE;
      fPhotonQualityCut = 0;
      fPhotonQualityCutTRD = 0;
      fPhotonQualityCutTOF = 0;
      break;
    case 1:
      fDoSharedElecCut = kTRUE;
      fDoPhotonQualitySelectionCut = kFALSE;
      fDoPhotonQualityRejectionCut = kFALSE;
      fPhotonQualityCut = 0;
      fPhotonQualityCutTRD = 0;
      fPhotonQualityCutTOF = 0;
      break;
    case 2:
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kTRUE;
      fDoPhotonQualityRejectionCut = kFALSE;
      fPhotonQualityCut = 1;
      fPhotonQualityCutTRD = 0;
      fPhotonQualityCutTOF = 0;
      break;
    case 3:
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kTRUE;
      fDoPhotonQualityRejectionCut = kFALSE;
      fPhotonQualityCut = 2;
      fPhotonQualityCutTRD = 0;
      fPhotonQualityCutTOF = 0;
      break;
    case 4:
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kTRUE;
      fDoPhotonQualityRejectionCut = kFALSE;
      fPhotonQualityCut = 3;
      fPhotonQualityCutTRD = 0;
      fPhotonQualityCutTOF = 0;
      break;
    case 5://Cat1 rejection
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kFALSE;
      fDoPhotonQualityRejectionCut = kTRUE;
      fPhotonQualityCut = 1;
      fPhotonQualityCutTRD = 0;
      fPhotonQualityCutTOF = 0;
      break;
    case 6: // reject TPC-only photons with TOF only
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kFALSE;
      fDoPhotonQualityRejectionCut = kTRUE;
      fPhotonQualityCut = 1;
      fPhotonQualityCutTRD = 0;
      fPhotonQualityCutTOF = 1;
      break;
    case 7: // reject TPC-only photons with TRD only
      fDoSharedElecCut = kFALSE;
      fDoPhotonQualitySelectionCut = kFALSE;
      fDoPhotonQualityRejectionCut = kTRUE;
      fPhotonQualityCut = 1;
      fPhotonQualityCutTRD = 1;
      fPhotonQualityCutTOF = 0;
      break;
    default:
      AliError(Form("Shared Electron Cut not defined %d",sharedElec));
      return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetToCloseV0sCut(Int_t toClose) {

  switch(toClose){
  case 0:
    fDoToCloseV0sCut = kFALSE;
    fminV0Dist = 250;
    break;
  case 1:
    fDoToCloseV0sCut = kTRUE;
    fminV0Dist = 1;
    break;
  case 2:
    fDoToCloseV0sCut = kTRUE;
    fminV0Dist = 2;
    break;
  case 3:
    fDoToCloseV0sCut = kTRUE;
    fminV0Dist = 3;
    break;
  case 4:
    fDoToCloseV0sCut = kTRUE;
    fDoDoubleCountingCut = kTRUE;
    fMinRDC=0.;
    fDeltaR=6.;
    fOpenAngle=0.02;
    break;
  case 5:
    fDoToCloseV0sCut = kTRUE;
    fDoDoubleCountingCut = kTRUE;
    fMinRDC=0.;
    fDeltaR=6.;
    fOpenAngle=0.03;
    break;
  case 6:
    fDoToCloseV0sCut = kTRUE;
    fDoDoubleCountingCut = kTRUE;
    fMinRDC=0.;
    fDeltaR=6.;
    fOpenAngle=0.04;
    break;

  default:
    AliError(Form("Shared Electron Cut not defined %d",toClose));
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetTRDElectronCut(Int_t TRDElectronCut){   // Set Cut
  switch(TRDElectronCut){
  case 0:
    fDoTRDPID=kFALSE;
    break;
  case 1:
    fDoTRDPID=kTRUE;
    fPIDTRDEfficiency=0.1;
    break;
  case 8:
    fDoTRDPID=kTRUE;
    fPIDTRDEfficiency=0.8;
    break;
  case 9:
    fDoTRDPID=kTRUE;
    fPIDTRDEfficiency=0.9;
    break;
  default:
    AliError(Form("TRDElectronCut not defined %d",TRDElectronCut));
    return kFALSE;
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetDCAZPhotonPrimVtxCut(Int_t DCAZPhotonPrimVtx){
  // Set Cut
  switch(DCAZPhotonPrimVtx){
  case 0:  //
    fDCAZPrimVtxCut   = 1000;
    break;
  case 1:  //
    fDCAZPrimVtxCut   = 10;
    break;
  case 2:  //
    fDCAZPrimVtxCut   = 5;
    break;
  case 3:  //
    fDCAZPrimVtxCut   = 4;
    break;
  case 4:  //
    fDCAZPrimVtxCut   = 3;
    break;
  case 5:  //
    fDCAZPrimVtxCut   = 2.5;
    break;
  case 6:  //
    fDCAZPrimVtxCut   = 2;
    break;
  case 7:  //
    fDCAZPrimVtxCut   = 1.5;
    break;
  case 8:  //
    fDCAZPrimVtxCut   = 1;
    break;
  case 9:  //
    fDCAZPrimVtxCut   = 0.5;
    break;
  default:
    cout<<"Warning: DCAZPhotonPrimVtx not defined "<<DCAZPhotonPrimVtx<<endl;
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetDCARPhotonPrimVtxCut(Int_t DCARPhotonPrimVtx){
  // Set Cut
  switch(DCARPhotonPrimVtx){
  case 0:  //
    fDCARPrimVtxCut   = 1000;
    break;
  case 1:  //
    fDCARPrimVtxCut   = 10;
    break;
  case 2:  //
    fDCARPrimVtxCut   = 5;
    break;
  case 3:  //
    fDCARPrimVtxCut   = 4;
    break;
  case 4:  //
    fDCARPrimVtxCut   = 3;
    break;
  case 5:  //
    fDCARPrimVtxCut   = 2.5;
    break;
  case 6:  //
    fDCARPrimVtxCut   = 2;
    break;
  case 7:  //
    fDCARPrimVtxCut   = 1.5;
    break;
  case 8:  //
    fDCARPrimVtxCut   = 1;
    break;
  case 9:  //
    fDCARPrimVtxCut   = 0.5;
    break;
  default:
    cout<<"Warning: DCARPhotonPrimVtx not defined "<<DCARPhotonPrimVtx<<endl;
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::SetInPlaneOutOfPlane(Int_t inOutPlane){
  // Set Cut
  switch(inOutPlane){
  case 0:  //
    fInPlaneOutOfPlane = 0; // No Event Plane
    break;
  case 1:  //
    fInPlaneOutOfPlane = 1; // In-Plane
    break;
  case 2:  //
    fInPlaneOutOfPlane = 2; // Out-Of-Plane
    break;
  default:
    cout<<"Warning: In-Plane or Out-Of-Plane not defined "<<inOutPlane<<endl;
    return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
void AliConversionPhotonCuts::GetPhiRegions(){
    fGoodRegionCMin = 5.0; fGoodRegionCMax = 6.2;
    fGoodRegionAMin = 3.5; fGoodRegionAMax = 6.2;
    fBadRegionCMin  = 0.0; fBadRegionCMax  = 1.5;
    fBadRegionAMin  = 0.0; fBadRegionAMax  = 2.5;
}

///________________________________________________________________________
Int_t AliConversionPhotonCuts::GetFirstTPCRow(Double_t radius){
  // Get first TPC row
  Int_t firstTPCRow = 0;
  Double_t radiusI = 84.8;
  Double_t radiusO = 134.6;
  Double_t radiusOB = 198.;
  Double_t rSizeI = 0.75;
  Double_t rSizeO = 1.;
  Double_t rSizeOB = 1.5;
  Int_t nClsI = 63;
  Int_t nClsIO = 127;

  if(radius <= radiusI){
    return firstTPCRow;
  }
  if(radius>radiusI && radius<=radiusO){
    firstTPCRow = (Int_t)((radius-radiusI)/rSizeI);
  }
  if(radius>radiusO && radius<=radiusOB){
    firstTPCRow = (Int_t)(nClsI+(radius-radiusO)/rSizeO);
  }

  if(radius>radiusOB){
    firstTPCRow =(Int_t)(nClsIO+(radius-radiusOB)/rSizeOB);
  }

  return firstTPCRow;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::CosinePAngleCut(const AliConversionPhotonBase * photon, AliVEvent * event) const {
  ///Check if passes cosine of pointing angle cut
  if(GetCosineOfPointingAngle(photon, event) < fCosPAngleCut){
    return kFALSE;
  }
  return kTRUE;
}

///________________________________________________________________________
Double_t AliConversionPhotonCuts::GetCosineOfPointingAngle( const AliConversionPhotonBase * photon, AliVEvent * event) const{
  // calculates the pointing angle of the recalculated V0

  Double_t momV0[3] = {0,0,0};
  if(event->IsA()==AliESDEvent::Class()){
    AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(event);
    if(!esdEvent) return -999;
    AliESDv0 *v0 = esdEvent->GetV0(photon->GetV0Index());
    if(!v0) return -999;
    v0->GetPxPyPz(momV0[0],momV0[1],momV0[2]);
  }
  if(event->IsA()==AliAODEvent::Class()){
    momV0[0] = photon->GetPx();
    momV0[1] = photon->GetPy();
    momV0[2] = photon->GetPz();
  }

  //Double_t momV0[3] = { photon->GetPx(), photon->GetPy(), photon->GetPz() }; //momentum of the V0
  Double_t PosV0[3] = { photon->GetConversionX() - event->GetPrimaryVertex()->GetX(),
              photon->GetConversionY() - event->GetPrimaryVertex()->GetY(),
              photon->GetConversionZ() - event->GetPrimaryVertex()->GetZ() }; //Recalculated V0 Position vector

  Double_t momV02 = momV0[0]*momV0[0] + momV0[1]*momV0[1] + momV0[2]*momV0[2];
  Double_t PosV02 = PosV0[0]*PosV0[0] + PosV0[1]*PosV0[1] + PosV0[2]*PosV0[2];


  Double_t cosinePointingAngle = -999;
  if(momV02*PosV02 > 0.0)
    cosinePointingAngle = (PosV0[0]*momV0[0] +  PosV0[1]*momV0[1] + PosV0[2]*momV0[2] ) / TMath::Sqrt(momV02 * PosV02);

  return cosinePointingAngle;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::PsiPairCut(const AliConversionPhotonBase * photon) const {
  if (fDo2DPsiPairChi2==1){
    if(fIncludeRejectedPsiPair){
      if (TMath::Abs(photon->GetPsiPair()) < -fPsiPairCut/fChi2CutConversion*photon->GetChi2perNDF() + fPsiPairCut || (photon->GetPsiPair()) == 4){
        return kTRUE;
      } else {
        return kFALSE;
      }

    } else {
      if (TMath::Abs(photon->GetPsiPair()) < -fPsiPairCut/fChi2CutConversion*photon->GetChi2perNDF() + fPsiPairCut ){
        return kTRUE;
      } else {
        return kFALSE;
      }
    }
  } else if (fDo2DPsiPairChi2==2){
    if(fIncludeRejectedPsiPair){
      if (TMath::Abs(photon->GetPsiPair()) < fPsiPairCut*TMath::Exp(fChi2CutConversionExpFunc*photon->GetChi2perNDF()) || (photon->GetPsiPair()) == 4){
        return kTRUE;
      } else {
        return kFALSE;
      }

    } else {
      if (TMath::Abs(photon->GetPsiPair()) < fPsiPairCut*TMath::Exp(fChi2CutConversionExpFunc*photon->GetChi2perNDF()) ){
        return kTRUE;
      } else {
        return kFALSE;
      }
    }
  } else {
    if(fIncludeRejectedPsiPair){
      if(TMath::Abs(photon->GetPsiPair()) > fPsiPairCut || (photon->GetPsiPair()) != 4){
        return kFALSE;
      } else {
        return kTRUE;
      }
    } else {
      if(TMath::Abs(photon->GetPsiPair()) > fPsiPairCut){
        return kFALSE;
      } else {
        return kTRUE;
      }
    }
  }
}

///________________________________________________________________________
TString AliConversionPhotonCuts::GetCutNumber(){
  // returns TString with current cut number
  return fCutStringRead;
}

///________________________________________________________________________
void AliConversionPhotonCuts::FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0){

  Int_t posLabel = photon->GetTrackLabelPositive();
  Int_t negLabel = photon->GetTrackLabelNegative();

  fElectronLabelArray[nV0*2] = posLabel;
  fElectronLabelArray[(nV0*2)+1] = negLabel;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s){

  Int_t posLabel = photon->GetTrackLabelPositive();
  Int_t negLabel = photon->GetTrackLabelNegative();

  for(Int_t i = 0; i<nV0s*2;i++){
    if(i==nV0*2)     continue;
    if(i==(nV0*2)+1) continue;
    if(fElectronLabelArray[i] == posLabel){
      return kFALSE;}
    if(fElectronLabelArray[i] == negLabel){
      return kFALSE;}
  }

  return kTRUE;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0){

  if (fDoDoubleCountingCut && photon->GetConversionRadius() < fMinRDC) return kTRUE;

  Double_t posX = photon->GetConversionX();
  Double_t posY = photon->GetConversionY();
  Double_t posZ = photon->GetConversionZ();

  for(Int_t i = 0;i<photons->GetEntries();i++){
    if(nV0 == i) continue;
    AliAODConversionPhoton *photonComp = (AliAODConversionPhoton*) photons->At(i);
    Double_t posCompX = photonComp->GetConversionX();
    Double_t posCompY = photonComp->GetConversionY();
    Double_t posCompZ = photonComp->GetConversionZ();

    if (!fDoDoubleCountingCut){
      Double_t dist = pow((posX - posCompX),2)+pow((posY - posCompY),2)+pow((posZ - posCompZ),2);

      if(dist < fminV0Dist*fminV0Dist){
        if(photon->GetChi2perNDF() > photonComp->GetChi2perNDF()) return kFALSE;
      }
    }else{
      TVector3 v1(photon->Px(),photon->Py(),photon->Pz());
      TVector3 v2(photonComp->Px(),photonComp->Py(),photonComp->Pz());
      Double_t OpeningAngle=v1.Angle(v2);
      if( OpeningAngle < fOpenAngle && TMath::Abs(photon->GetConversionRadius()-photonComp->GetConversionRadius()) < fDeltaR){
        if(photon->GetChi2perNDF() > photonComp->GetChi2perNDF()) return kFALSE;
      }
    }

  }
  return kTRUE;
}


///________________________________________________________________________
AliConversionPhotonCuts* AliConversionPhotonCuts::GetStandardCuts2010PbPb(){
  //Create and return standard 2010 PbPb cuts
  AliConversionPhotonCuts *cuts=new AliConversionPhotonCuts("StandardCuts2010PbPb","StandardCuts2010PbPb");
  if(!cuts->InitializeCutsFromCutString("04209297002322000000")){
    cout<<"Warning: Initialization of Standardcuts2010PbPb failed"<<endl;}
  return cuts;
}

///________________________________________________________________________
AliConversionPhotonCuts* AliConversionPhotonCuts::GetStandardCuts2010pp(){
  //Create and return standard 2010 PbPb cuts
  AliConversionPhotonCuts *cuts=new AliConversionPhotonCuts("StandardCuts2010pp","StandardCuts2010pp");
  if(!cuts->InitializeCutsFromCutString("00209366300380000000")){
    cout<<"Warning: Initialization of Standardcuts2010pp failed"<<endl;}
  return cuts;
}

///________________________________________________________________________
Bool_t AliConversionPhotonCuts::InPlaneOutOfPlaneCut(Double_t photonPhi, Double_t eventPlaneAngle, Bool_t fill){

  //GetPhotonPhi() 0-2 Pi  //eventPlaneAngle -1pi-1pi
  eventPlaneAngle=eventPlaneAngle+TMath::Pi();
  Double_t gammaToEPAngle = eventPlaneAngle-photonPhi;
  if(gammaToEPAngle < 0) gammaToEPAngle=gammaToEPAngle+2*TMath::Pi();
  gammaToEPAngle = gammaToEPAngle-TMath::Pi(); // angle from -pi +pi

  if(!fInPlaneOutOfPlane){
    if(fill&&fHistoEventPlanePhi)fHistoEventPlanePhi->Fill(gammaToEPAngle);
    return kTRUE;
  }
  else if(fInPlaneOutOfPlane == 1){
    if(TMath::Abs(gammaToEPAngle)<=0.25*TMath::Pi() || TMath::Abs(gammaToEPAngle)>=0.75*TMath::Pi()){
      if(fill&&fHistoEventPlanePhi)fHistoEventPlanePhi->Fill(gammaToEPAngle);
      return kTRUE;
    }
    else return kFALSE;
  }
  else if(fInPlaneOutOfPlane == 2){
    if(TMath::Abs(gammaToEPAngle)>0.25*TMath::Pi() && TMath::Abs(gammaToEPAngle)<0.75*TMath::Pi()){
      if(fill&&fHistoEventPlanePhi)fHistoEventPlanePhi->Fill(gammaToEPAngle);
      return kTRUE;
    }
    else return kFALSE;
  }
  return kFALSE;
}

///________________________________________________________________________
UChar_t AliConversionPhotonCuts::DeterminePhotonQualityAOD(AliAODConversionPhoton* photon, AliVEvent* eventDummy){

  AliAODTrack * negTrack = static_cast<AliAODTrack*>(GetTrack(eventDummy, photon->GetTrackLabelNegative()));
  AliAODTrack * posTrack = static_cast<AliAODTrack*>(GetTrack(eventDummy, photon->GetTrackLabelPositive()));

  if(!negTrack || !posTrack) {
      return 0;
  }
  if(negTrack->Charge() == posTrack->Charge()){
      return 0;
  }
  Int_t nClusterITSneg = negTrack->GetITSNcls();
  Int_t nClusterITSpos = posTrack->GetITSNcls();
  //    cout << nClusterITSneg << "\t" << nClusterITSpos <<endl;

  if (nClusterITSneg > 1 && nClusterITSpos > 1){
    return 3;
  } else if (nClusterITSneg > 1 || nClusterITSpos > 1){
    return 2;
  } else {
    return 1;
  }
  return 0;
}

///________________________________________________________________________
UChar_t AliConversionPhotonCuts::DeterminePhotonQualityTRD(AliAODConversionPhoton* photon, AliVEvent* eventDummy){

  AliVTrack* negTrack = GetTrack(eventDummy, photon->GetTrackLabelNegative());
  AliVTrack* posTrack = GetTrack(eventDummy, photon->GetTrackLabelPositive());

  if(!negTrack || !posTrack) {
      return 0;
  }
  if(negTrack->Charge() == posTrack->Charge()){
      return 0;
  }

  Int_t nClusterTRDneg = negTrack->GetNcls(2);
  Int_t nClusterTRDpos = posTrack->GetNcls(2);
  
  if (nClusterTRDneg > 1 && nClusterTRDpos > 1){
    return 3;
  } else if (nClusterTRDneg > 1 || nClusterTRDpos > 1){
    return 2;
  } else {
    return 1;
  }
}


///________________________________________________________________________
UChar_t AliConversionPhotonCuts::DeterminePhotonQualityTOF(AliAODConversionPhoton* photon, AliVEvent* eventDummy){

  AliVTrack* negTrack = GetTrack(eventDummy, photon->GetTrackLabelNegative());
  AliVTrack* posTrack = GetTrack(eventDummy, photon->GetTrackLabelPositive());

  if(!negTrack || !posTrack) {
      return 0;
  }
  if(negTrack->Charge() == posTrack->Charge()){
      return 0;
  }

  Bool_t negTOFSignal = (negTrack->GetStatus() & AliVTrack::kTOFout ) && (negTrack->GetStatus() & AliVTrack::kTIME);
  Bool_t posTOFSignal = (posTrack->GetStatus() & AliVTrack::kTOFout ) && (posTrack->GetStatus() & AliVTrack::kTIME);

  if (negTOFSignal && posTOFSignal){
    return 3;
  } else if (negTOFSignal || posTOFSignal){
    return 2;
  } else {
    return 1;
  }
}

///__________________________________________________________________________________________
Bool_t AliConversionPhotonCuts::InitializeMaterialBudgetWeights(Int_t flag, TString filename){

    TString nameProfile;
    if      (flag==1){
                nameProfile = "profile2DContainingMaterialBudgetWeights_fewRadialBins";}
    else if (flag==2){
                nameProfile = "profile2DContainingMaterialBudgetWeights_manyRadialBins";}
    else {
        AliError(Form("%d not a valid flag for InitMaterialBudgetWeightingOfPi0Candidates()",flag));
        return kFALSE;
    }
    TFile* file = TFile::Open(filename.Data());
    if (!file) {
        AliError(Form("File %s for materialbudgetweights not found",filename.Data()));
        return kFALSE;
    }
    fProfileContainingMaterialBudgetWeights = (TProfile2D*)file->Get(nameProfile.Data());
    if (!fProfileContainingMaterialBudgetWeights){
        AliError(Form("Histogram %s not found in file",nameProfile.Data()));
        return kFALSE;
    }
    fProfileContainingMaterialBudgetWeights->SetDirectory(0);
    file->Close();
    delete file;

    fMaterialBudgetWeightsInitialized = kTRUE;
    AliInfo(Form("MaterialBudgetWeightingOfPi0Candidates initialized with flag %d. This means %d radial bins will be used for the weighting. File used: %s.",flag, fProfileContainingMaterialBudgetWeights->GetNbinsX(), filename.Data()));
    return kTRUE;
}

///___________________________________________________________________________________________________
 Float_t AliConversionPhotonCuts::GetMaterialBudgetCorrectingWeightForTrueGamma(AliAODConversionPhoton* gamma, Double_t magField){

    Float_t weight = 1.0;
    Float_t gammaConversionRadius = gamma->GetConversionRadius();
    Float_t scalePt=1.;
    Float_t nomMagField = 5.;
    if(magField!=0) 
      scalePt = nomMagField/(TMath::Abs(magField));
    
    // AM:  Scale the pT for correction in case of lowB field
    //    cout<< "scalePt::"<< scalePt<< "    " <<  magField<< endl;

    //AM.  the Omega correction for pT > 0.4 is flat and at high pT the statistics reduces. 
    // So take the correction  at pT=0.5 if pT is > 0.7 GeV/c
    Float_t maxPtForCor = 0.7;  
    Float_t defaultPtForCor = 0.5;  
    Float_t gammaPt = scalePt * gamma->Pt();


    Int_t binX = fProfileContainingMaterialBudgetWeights->GetXaxis()->FindBin(gammaConversionRadius+0.001);
    Int_t binY;

    if (gammaPt < maxPtForCor){
      binY = fProfileContainingMaterialBudgetWeights->GetYaxis()->FindBin(gammaPt+0.001);
    }  else{
      binY = fProfileContainingMaterialBudgetWeights->GetYaxis()->FindBin(defaultPtForCor+0.001);
    }
    if (  (binX > 0 && binX <= fProfileContainingMaterialBudgetWeights->GetNbinsX()) &&
	  (binY > 0 && binY <= fProfileContainingMaterialBudgetWeights->GetNbinsY())){
      weight = fProfileContainingMaterialBudgetWeights->GetBinContent(binX,binY);
    }
    //    cout << gammaConversionRadius<< " " << gammaPt << " " << binX<< " " << binY << " "<<  weight<< endl;
    return weight;
}
