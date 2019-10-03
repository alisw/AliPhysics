
//ERRORs: no member SetCutIsCowboy(int i)

#include <stdio.h>

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TProfile.h"
//#include "TFriendElement.h"
#include "TLorentzVector.h"

//not in $ALICE_PHYSICS/include
//#include "AliMultSelection.h"
#include "AliMultSelection.h"
#include "AliLightV0vertexer.h"
#include "AliLightCascadeVertexer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "AliV0Result.h"
#include "AliCascadeResult.h"
#include "AliCFContainer.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliOADBContainer.h"
#include "AliTriggerIR.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"

//#include "AliRunningCascadeEvent.h"
//#include "AliRunningCascadeTrack.h"
#include "AliAnalysisTaskStrangeCascadesDiscrete.h"

//ClassImp(AliRunningCascadeTrack)
//ClassImp(AliRunningCascadeEvent)
ClassImp(AliAnalysisTaskStrangeCascadesDiscrete)


AliAnalysisTaskStrangeCascadesDiscrete::AliAnalysisTaskStrangeCascadesDiscrete()
: AliAnalysisTaskSE(), fListHist(0), fListK0Short(0), fListLambda(0), fListAntiLambda(0),
fListXiMinus(0), fListXiPlus(0), fListOmegaMinus(0), fListOmegaPlus(0),fEventList(0),
fTreeEvent(0), // fTreeV0(0), //
fTreeCascade(0),
Cascade_Track(0),
Cascade_Event(0),
fTreeCascadeAsEvent(0),
fPIDResponse(0), fESDtrackCuts(0),
fESDtrackCutsITSsa2010(0), fESDtrackCutsGlobal2015(0),
fUtils(0),
number(0),
fEventCuts(0),
fEventCutsStrictAntipileup(0),

fRand(0),
//---> Flags controlling Event Tree output
fkSaveEventTree    ( kTRUE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
fkSaveV0Tree       ( kTRUE ),
fkDownScaleV0      ( kTRUE  ),
fDownScaleFactorV0 ( 0.001  ),
fkPreselectDedx ( kTRUE ),
fkUseOnTheFlyV0Cascading( kFALSE ),
fkDebugWrongPIDForTracking ( kFALSE ),
fkDebugBump(kFALSE),
fkDebugOOBPileup(kFALSE),
fkDoExtraEvSels(kTRUE),
fkPileupRejectionMode(0),
fkUseOldCentrality ( kFALSE ) ,

//---> Flags controlling Cascade TTree output
fkSaveCascadeTree       ( kTRUE  ),
fkDownScaleCascade      ( kTRUE  ),
fDownScaleFactorCascade ( 0.001  ),

//---> Flags controlling Vertexers
fkRunVertexers    ( kFALSE ),
fkUseLightVertexer ( kTRUE ),
fkDoV0Refit       ( kTRUE ),
fkExtraCleanup    ( kTRUE ),

//---> Flag controlling trigger selection
fTrigType(AliVEvent::kMB),

//---> Min pT to save candidate
fMinPtToSave( 0.0 ),
fMaxPtToSave( 100.00 ),

//---> Sandbox switch
fkSandboxMode ( kFALSE ),

//---> Variables for fTreeEvent
fCentrality(0),
fMVPileupFlag(kFALSE),
fOOBPileupFlag(kFALSE),
fNTOFClusters(-1),
fNTOFMatches(-1),
fNTracksITSsa2010(-1),
fNTracksGlobal2015(-1),
fNTracksGlobal2015TriggerPP(-1),
fAmplitudeV0A(-1.),
fAmplitudeV0C(-1.),
fClosestNonEmptyBC(-1),

//---> Variables for fTreeV0
fTreeVariableChi2V0(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0ToPrimVertex(0),
fTreeVariableDcaPosToPrimVertex(0),
fTreeVariableDcaNegToPrimVertex(0),
fTreeVariableV0CosineOfPointingAngle(0),
fTreeVariableV0Radius(0),
fTreeVariablePt(0),
fTreeVariableRapK0Short(0),
fTreeVariableRapLambda(0),
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariableAlphaV0(0),
fTreeVariablePtArmV0(0),
fTreeVariableNegEta(0),
fTreeVariablePosEta(0),

fTreeVariableNSigmasPosProton(0),
fTreeVariableNSigmasPosPion(0),
fTreeVariableNSigmasNegProton(0),
fTreeVariableNSigmasNegPion(0),

fTreeVariableDistOverTotMom(0),
fTreeVariableLeastNbrCrossedRows(0),
fTreeVariableLeastRatioCrossedRowsOverFindable(0),
fTreeVariableMaxChi2PerCluster(0),
fTreeVariableMinTrackLength(0),

fTreeVariablePosPIDForTracking(-1),
fTreeVariableNegPIDForTracking(-1),
fTreeVariablePosdEdx(-1),
fTreeVariableNegdEdx(-1),
fTreeVariablePosInnerP(-1),
fTreeVariableNegInnerP(-1),
fTreeVariableNegTrackStatus(0),
fTreeVariablePosTrackStatus(0),
fTreeVariableNegDCAz(-1),
fTreeVariablePosDCAz(-1),

fTreeVariablePosITSClusters0(0),
fTreeVariablePosITSClusters1(0),
fTreeVariablePosITSClusters2(0),
fTreeVariablePosITSClusters3(0),
fTreeVariablePosITSClusters4(0),
fTreeVariablePosITSClusters5(0),

fTreeVariableNegITSClusters0(0),
fTreeVariableNegITSClusters1(0),
fTreeVariableNegITSClusters2(0),
fTreeVariableNegITSClusters3(0),
fTreeVariableNegITSClusters4(0),
fTreeVariableNegITSClusters5(0),

fTreeVariablePosITSSharedClusters0(0),
fTreeVariablePosITSSharedClusters1(0),
fTreeVariablePosITSSharedClusters2(0),
fTreeVariablePosITSSharedClusters3(0),
fTreeVariablePosITSSharedClusters4(0),
fTreeVariablePosITSSharedClusters5(0),

fTreeVariableNegITSSharedClusters0(0),
fTreeVariableNegITSSharedClusters1(0),
fTreeVariableNegITSSharedClusters2(0),
fTreeVariableNegITSSharedClusters3(0),
fTreeVariableNegITSSharedClusters4(0),
fTreeVariableNegITSSharedClusters5(0),

fTreeVariableIsCowboy(0),

fTreeVariableNegTOFExpTDiff(99999),
fTreeVariablePosTOFExpTDiff(99999),
fTreeVariableNegTOFSignal(99999),
fTreeVariablePosTOFSignal(99999),
fTreeVariableNegTOFBCid(99999),
fTreeVariablePosTOFBCid(99999),
fTreeVariableAmplitudeV0A(-1.),
fTreeVariableAmplitudeV0C(-1.),
fTreeVariableClosestNonEmptyBC(-1),

fTreeVariableCentrality(0),
fTreeVariableMVPileupFlag(kFALSE),
fTreeVariableOOBPileupFlag(kFALSE),

fTreeVariablePrimVertexX(0),
fTreeVariablePrimVertexY(0),
fTreeVariablePrimVertexZ(0),

fTreeVariablePosTrack(0x0),
fTreeVariableNegTrack(0x0),

fTreeVariableMagneticField(0x0),

//---> Variables for fTreeCascade
fTreeCascVarCharge(0),
fTreeCascVarMassAsXi(0),
fTreeCascVarMassAsOmega(0),
fTreeCascVarPt(0),
fTreeCascVarP(0),
fTreeCascVarPz(0),
fTreeCascVarRapXi(0),
fTreeCascVarRapOmega(0),
fTreeCascVarNegEta(0),
fTreeCascVarPosEta(0),
fTreeCascVarBachEta(0),
fTreeCascVarDCACascDaughters(0),
fTreeCascVarDCABachToPrimVtx(0),
fTreeCascVarDCABachToPrimVtxZ(0),
fTreeCascVarDCAV0Daughters(0),
fTreeCascVarDCAV0ToPrimVtx(0),
fTreeCascVarDCAPosToPrimVtx(0),
fTreeCascVarDCAPosToPrimVtxZ(0),
fTreeCascVarDCANegToPrimVtx(0),
fTreeCascVarDCANegToPrimVtxZ(0),
fTreeCascVarCascCosPointingAngle(0),
fTreeCascVarCascDCAtoPVxy(0),
fTreeCascVarCascDCAtoPVz(0),
fTreeCascVarCascRadius(0),
fTreeCascVarV0Mass(0),
fTreeCascVarV0MassLambda(0),
fTreeCascVarV0MassAntiLambda(0),
fTreeCascVarV0CosPointingAngle(0),
fTreeCascVarV0CosPointingAngleSpecial(0),
fTreeCascVarV0Radius(0),
fTreeCascVarDCABachToBaryon(0),
fTreeCascVarWrongCosPA(0),
fTreeCascVarLeastNbrClusters(0),
fTreeCascVarDistOverTotMom(0),
fTreeCascVarMaxChi2PerCluster(0),
fTreeCascVarMinTrackLength(0),

fTreeCascVarNegNSigmaPion(0),
fTreeCascVarNegNSigmaProton(0),
fTreeCascVarPosNSigmaPion(0),
fTreeCascVarPosNSigmaProton(0),
fTreeCascVarBachNSigmaPion(0),
fTreeCascVarBachNSigmaKaon(0),

fTreeCascVarNegTOFNSigmaPion(0),
fTreeCascVarNegTOFNSigmaProton(0),
fTreeCascVarPosTOFNSigmaPion(0),
fTreeCascVarPosTOFNSigmaProton(0),
fTreeCascVarBachTOFNSigmaPion(0),
fTreeCascVarBachTOFNSigmaKaon(0),

fTreeCascVarNegITSNSigmaPion(0),
fTreeCascVarNegITSNSigmaProton(0),
fTreeCascVarPosITSNSigmaPion(0),
fTreeCascVarPosITSNSigmaProton(0),
fTreeCascVarBachITSNSigmaPion(0),
fTreeCascVarBachITSNSigmaKaon(0),


fTreeCascVarChiSquareV0(1e+3),
fTreeCascVarChiSquareCascade(1e+3),

fTreeCascVarBachDCAPVSigmaX2(0),
fTreeCascVarBachDCAPVSigmaY2(0),
fTreeCascVarBachDCAPVSigmaZ2(0),
fTreeCascVarPosDCAPVSigmaX2(0),
fTreeCascVarPosDCAPVSigmaY2(0),
fTreeCascVarPosDCAPVSigmaZ2(0),
fTreeCascVarNegDCAPVSigmaX2(0),
fTreeCascVarNegDCAPVSigmaY2(0),
fTreeCascVarNegDCAPVSigmaZ2(0),

fTreeCascVarPosPIDForTracking(-1),
fTreeCascVarNegPIDForTracking(-1),
fTreeCascVarBachPIDForTracking(-1),
fTreeCascVarNegInnerP(-1),
fTreeCascVarPosInnerP(-1),
fTreeCascVarBachInnerP(-1),
fTreeCascVarNegdEdx(-1),
fTreeCascVarPosdEdx(-1),
fTreeCascVarBachdEdx(-1),

fTreeCascVarNegTrackStatus(0), //!
fTreeCascVarPosTrackStatus(0), //!
fTreeCascVarBachTrackStatus(0), //!
fTreeCascVarNegDCAz(-1),
fTreeCascVarPosDCAz(-1),
fTreeCascVarBachDCAz(-1),
//fTreeCascVarPosTotMom(-1),
//fTreeCascVarNegTotMom(-1),
//fTreeCascVarBachTotMom(-1),

fTreeCascVarNegPx(0),
fTreeCascVarNegPy(0),
fTreeCascVarNegPz(0),
fTreeCascVarPosPx(0),
fTreeCascVarPosPy(0),
fTreeCascVarPosPz(0),
fTreeCascVarBachPx(0),
fTreeCascVarBachPy(0),
fTreeCascVarBachPz(0),
fTreeCascVarV0DecayX(0),
fTreeCascVarV0DecayY(0),
fTreeCascVarV0DecayZ(0),
fTreeCascVarCascadeDecayX(0),
fTreeCascVarCascadeDecayY(0),
fTreeCascVarCascadeDecayZ(0),
fTreeCascVarV0Lifetime(0),
//Track Labels (check for duplicates, etc)
fTreeCascVarNegIndex(0),
fTreeCascVarPosIndex(0),
fTreeCascVarBachIndex(0),
//Event Number (check same-event index mixups)
fTreeCascVarEventNumber(0),

fTreeCascVarPosITSClusters0(0),
fTreeCascVarPosITSClusters1(0),
fTreeCascVarPosITSClusters2(0),
fTreeCascVarPosITSClusters3(0),
fTreeCascVarPosITSClusters4(0),
fTreeCascVarPosITSClusters5(0),

fTreeCascVarNegITSClusters0(0),
fTreeCascVarNegITSClusters1(0),
fTreeCascVarNegITSClusters2(0),
fTreeCascVarNegITSClusters3(0),
fTreeCascVarNegITSClusters4(0),
fTreeCascVarNegITSClusters5(0),

fTreeCascVarBachITSClusters0(0),
fTreeCascVarBachITSClusters1(0),
fTreeCascVarBachITSClusters2(0),
fTreeCascVarBachITSClusters3(0),
fTreeCascVarBachITSClusters4(0),
fTreeCascVarBachITSClusters5(0),

fTreeCascVarPosITSSharedClusters0(0),
fTreeCascVarPosITSSharedClusters1(0),
fTreeCascVarPosITSSharedClusters2(0),
fTreeCascVarPosITSSharedClusters3(0),
fTreeCascVarPosITSSharedClusters4(0),
fTreeCascVarPosITSSharedClusters5(0),

fTreeCascVarNegITSSharedClusters0(0),
fTreeCascVarNegITSSharedClusters1(0),
fTreeCascVarNegITSSharedClusters2(0),
fTreeCascVarNegITSSharedClusters3(0),
fTreeCascVarNegITSSharedClusters4(0),
fTreeCascVarNegITSSharedClusters5(0),

fTreeCascVarBachITSSharedClusters0(0),
fTreeCascVarBachITSSharedClusters1(0),
fTreeCascVarBachITSSharedClusters2(0),
fTreeCascVarBachITSSharedClusters3(0),
fTreeCascVarBachITSSharedClusters4(0),
fTreeCascVarBachITSSharedClusters5(0),

//Variables for debugging the invariant mass bump
//Full momentum information

fTreeCascVarNegTOFExpTDiff(99999),
fTreeCascVarPosTOFExpTDiff(99999),
fTreeCascVarBachTOFExpTDiff(99999),
fTreeCascVarNegTOFSignal(99999),
fTreeCascVarPosTOFSignal(99999),
fTreeCascVarBachTOFSignal(99999),
fTreeCascVarNegTOFBCid(99999),
fTreeCascVarPosTOFBCid(99999),
fTreeCascVarBachTOFBCid(99999),
fTreeCascVarAmplitudeV0A(-1.),
fTreeCascVarAmplitudeV0C(-1.),
fTreeCascVarClosestNonEmptyBC(-1),

fTreeCascVarCentrality(0),
fTreeCascVarMVPileupFlag(kFALSE),
fTreeCascVarOOBPileupFlag(kFALSE),
//Kink tagging
fTreeCascVarBachIsKink(kFALSE),
fTreeCascVarPosIsKink(kFALSE),
fTreeCascVarNegIsKink(kFALSE),
fTreeCascVarIsCowboy(kFALSE),
fTreeCascVarCowboyness(-2),
fTreeCascVarIsCascadeCowboy(kFALSE),
fTreeCascVarCascadeCowboyness(-2),
fkSelectCharge(0),
fTreeCascVarPrimVertexX(0),
fTreeCascVarPrimVertexY(0),
fTreeCascVarPrimVertexZ(0),
fTreeCascVarBachTrack(0x0),
fTreeCascVarPosTrack(0x0),
fTreeCascVarNegTrack(0x0),
fTreeCascVarCosineProtonRestV0RestOmegaMinus(-100),
fTreeCascVarCosineProtonRestV0RestOmegaPlus(-100),
fTreeCascVarDzetaFromMomentaOmegaMinus(1e8),
fTreeCascVarMagneticField(0.0),
//extra
lEffMassXi(0.),
lNegTrackSign(1),
lPosTrackSign(1),
lBachTrackSign(1),
lEvSelCode(100),
lMinV0PA(0.0),
lMaxV0CosPA(1.0),
lMinV0CosPA(0.9),
lMaxCascCosPA(1.0),
lMinCascPA(0.0),
lMaxBBCosPA(1.0),
lMinBBPA(0.0),
dca(0.0),
fTreeCascVarMassXiScratch(0),
fTreeCascVarMassOmegaScratch(0),
fTreeCascVarBackgroundMassOmega(0),
fTreeCascVarBackgroundMassXi(0),
fTreeVariableAlphaXi(0),
fTreeVariablePtArmXi(0),

fkEsdTrackMultiplicity(0),
fkEsdEventTriggerWord(0),
fTreeCascVarBachP(0),
fTreeCascVarBachPt(0),
fTreeCascVarNegP(0),
fTreeCascVarNegPt(0),
fTreeCascVarPosP(0),
fTreeCascVarPosPt(0),
fkMomentaNeeded(kFALSE),
fkITSPID(kFALSE),
fkMassesAndBackgroundFromScratch(kFALSE),
fkDiscreteSymmetryInfo(kFALSE),
fkMultiplicityInfo(kFALSE),
fkExtraInfoCascadeKinematics(kFALSE),
fRunNum(0),
fNtracks(0),


//Histos
fHistEventCounter(0),
fHistEventCounterDifferential(0),
fHistCentrality(0),
fHistCentralityFast(0),
fHistZVertexEvent(0)
//------------------------------------------------
// Tree Variables
{
    
    
}

AliAnalysisTaskStrangeCascadesDiscrete::AliAnalysisTaskStrangeCascadesDiscrete(Bool_t lSaveEventTree,
                                                                                         Bool_t lSaveV0Tree,
                                                                                         Bool_t lSaveCascadeTree,
                                                                                         Bool_t lRunVertexers, //kFALSE per default
                                                                                         Bool_t lUseLightVertexer, //kFALSE per def
                                                                                         Bool_t lMomentaNeeded, //9 branches of p,pz,pt of daughters
                                                                                         Bool_t lExtraInfoCascadeKinematics,
                                                                                         Bool_t lDiscreteSymmetryInfo, //extra branch in cascade tree
                                                                                         Bool_t lMassesAndBackgroundFromScratch, //extra branch in cascade tree
                                                                                         Double_t lCascaderMaxChi2,
                                                                                         Double_t lCascaderV0MinImpactParam, //0.050 per def
                                                                                         Double_t lCascaderV0MassWindow, //0.010 per def
                                                                                         Double_t lCascaderBachMinImpactParam, //0.03 per def
                                                                                         Double_t lCascaderMaxDCAV0andBach, //2.0 per def
                                                                                         Double_t lCascaderMinCosAngle, //0.95 per def
                                                                                         Double_t lCascaderMinRadius, //0.4 per def
                                                                                         Double_t lCascaderMaxRadius, //100. per def
                                                                                         const char *name, TString lExtraOptions
                                                                                         )
: AliAnalysisTaskSE(name), fListHist(0), fListK0Short(0), fListLambda(0), fListAntiLambda(0),
fListXiMinus(0), fListXiPlus(0), fListOmegaMinus(0), fListOmegaPlus(0),fEventList(0),
fTreeEvent(0), //fTreeV0(0),//
fTreeCascade(0),
Cascade_Track(0),
Cascade_Event(0),
fTreeCascadeAsEvent(0),
fPIDResponse(0), fESDtrackCuts(0),
fESDtrackCutsITSsa2010(0), fESDtrackCutsGlobal2015(0),
fUtils(0),
number(0),
fEventCuts(0),
fEventCutsStrictAntipileup(0),
fRand(0),

//---> Flags controlling Event Tree output
fkSaveEventTree    ( kTRUE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
fkSaveV0Tree       ( kFALSE ),
fkDownScaleV0      ( kTRUE  ),
fDownScaleFactorV0 ( 0.001  ),
fkPreselectDedx ( kTRUE ),
fkUseOnTheFlyV0Cascading( kFALSE ),
fkDebugWrongPIDForTracking ( kFALSE ), //also for cascades...
fkDebugBump( kFALSE ),
fkDebugOOBPileup(kFALSE),
fkDoExtraEvSels(kTRUE),
fkPileupRejectionMode(0),
fkUseOldCentrality ( kFALSE ) ,

//---> Flags controlling Cascade TTree output
fkSaveCascadeTree       ( kTRUE  ),
fkDownScaleCascade      ( kTRUE  ),
fDownScaleFactorCascade ( 0.001  ),

//---> Flags controlling Vertexers
fkRunVertexers    ( kFALSE ),
fkUseLightVertexer ( kTRUE ),
fkDoV0Refit       ( kTRUE ),
fkExtraCleanup    ( kTRUE ),

//---> Flag controlling trigger selection
fTrigType(AliVEvent::kMB),

//---> Min pT to save candidate
fMinPtToSave( 0.0 ),
fMaxPtToSave( 100.00 ),

//---> Sandbox switch
fkSandboxMode ( kFALSE ),

//---> Variables for fTreeEvent
fCentrality(0),
fMVPileupFlag(kFALSE),
fOOBPileupFlag(kFALSE),
fNTOFClusters(-1),
fNTOFMatches(-1),
fNTracksITSsa2010(-1),
fNTracksGlobal2015(-1),
fNTracksGlobal2015TriggerPP(-1),
fAmplitudeV0A(-1.),
fAmplitudeV0C(-1.),
fClosestNonEmptyBC(-1),

//---> Variables for fTreeV0
fTreeVariableChi2V0(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0ToPrimVertex(0),
fTreeVariableDcaPosToPrimVertex(0),
fTreeVariableDcaNegToPrimVertex(0),
fTreeVariableV0CosineOfPointingAngle(0),
fTreeVariableV0Radius(0),
fTreeVariablePt(0),
fTreeVariableRapK0Short(0),
fTreeVariableRapLambda(0),
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariableAlphaV0(0),
fTreeVariablePtArmV0(0),
fTreeVariableNegEta(0),
fTreeVariablePosEta(0),

fTreeVariableNSigmasPosProton(0),
fTreeVariableNSigmasPosPion(0),
fTreeVariableNSigmasNegProton(0),
fTreeVariableNSigmasNegPion(0),

fTreeVariableDistOverTotMom(0),
fTreeVariableLeastNbrCrossedRows(0),
fTreeVariableLeastRatioCrossedRowsOverFindable(0),
fTreeVariableMaxChi2PerCluster(0),
fTreeVariableMinTrackLength(0),

fTreeVariablePosPIDForTracking(-1),
fTreeVariableNegPIDForTracking(-1),
fTreeVariablePosdEdx(-1),
fTreeVariableNegdEdx(-1),
fTreeVariablePosInnerP(-1),
fTreeVariableNegInnerP(-1),
fTreeVariableNegTrackStatus(0),
fTreeVariablePosTrackStatus(0),
fTreeVariableNegDCAz(-1),
fTreeVariablePosDCAz(-1),

fTreeVariablePosITSClusters0(0),
fTreeVariablePosITSClusters1(0),
fTreeVariablePosITSClusters2(0),
fTreeVariablePosITSClusters3(0),
fTreeVariablePosITSClusters4(0),
fTreeVariablePosITSClusters5(0),

fTreeVariableNegITSClusters0(0),
fTreeVariableNegITSClusters1(0),
fTreeVariableNegITSClusters2(0),
fTreeVariableNegITSClusters3(0),
fTreeVariableNegITSClusters4(0),
fTreeVariableNegITSClusters5(0),

fTreeVariablePosITSSharedClusters0(0),
fTreeVariablePosITSSharedClusters1(0),
fTreeVariablePosITSSharedClusters2(0),
fTreeVariablePosITSSharedClusters3(0),
fTreeVariablePosITSSharedClusters4(0),
fTreeVariablePosITSSharedClusters5(0),

fTreeVariableNegITSSharedClusters0(0),
fTreeVariableNegITSSharedClusters1(0),
fTreeVariableNegITSSharedClusters2(0),
fTreeVariableNegITSSharedClusters3(0),
fTreeVariableNegITSSharedClusters4(0),
fTreeVariableNegITSSharedClusters5(0),

fTreeVariableIsCowboy(0),

fTreeVariableNegTOFExpTDiff(99999),
fTreeVariablePosTOFExpTDiff(99999),
fTreeVariableNegTOFSignal(99999),
fTreeVariablePosTOFSignal(99999),
fTreeVariableNegTOFBCid(99999),
fTreeVariablePosTOFBCid(99999),
fTreeVariableAmplitudeV0A(-1.),
fTreeVariableAmplitudeV0C(-1.),
fTreeVariableClosestNonEmptyBC(-1),

fTreeVariableCentrality(0),
fTreeVariableMVPileupFlag(kFALSE),
fTreeVariableOOBPileupFlag(kFALSE),

fTreeVariablePrimVertexX(0),
fTreeVariablePrimVertexY(0),
fTreeVariablePrimVertexZ(0),

fTreeVariablePosTrack(0x0),
fTreeVariableNegTrack(0x0),

fTreeVariableMagneticField(0x0),

//---> Variables for fTreeCascade
fTreeCascVarCharge(0),
fTreeCascVarMassAsXi(0),
fTreeCascVarMassAsOmega(0),
fTreeCascVarPt(0),
fTreeCascVarP(0),
fTreeCascVarPz(0),
fTreeCascVarRapXi(0),
fTreeCascVarRapOmega(0),
fTreeCascVarNegEta(0),
fTreeCascVarPosEta(0),
fTreeCascVarBachEta(0),
fTreeCascVarDCACascDaughters(0),
fTreeCascVarDCABachToPrimVtx(0),
fTreeCascVarDCABachToPrimVtxZ(0),
fTreeCascVarDCAV0Daughters(0),
fTreeCascVarDCAV0ToPrimVtx(0),
fTreeCascVarDCAPosToPrimVtx(0),
fTreeCascVarDCAPosToPrimVtxZ(0),
fTreeCascVarDCANegToPrimVtx(0),
fTreeCascVarDCANegToPrimVtxZ(0),
fTreeCascVarCascCosPointingAngle(0),
fTreeCascVarCascDCAtoPVxy(0),
fTreeCascVarCascDCAtoPVz(0),
fTreeCascVarCascRadius(0),
fTreeCascVarV0Mass(0),
fTreeCascVarV0MassLambda(0),
fTreeCascVarV0MassAntiLambda(0),
fTreeCascVarV0CosPointingAngle(0),
fTreeCascVarV0CosPointingAngleSpecial(0),
fTreeCascVarV0Radius(0),
fTreeCascVarDCABachToBaryon(0),
fTreeCascVarWrongCosPA(0),
fTreeCascVarLeastNbrClusters(0),
fTreeCascVarDistOverTotMom(0),
fTreeCascVarMaxChi2PerCluster(0),
fTreeCascVarMinTrackLength(0),

fTreeCascVarNegNSigmaPion(0),
fTreeCascVarNegNSigmaProton(0),
fTreeCascVarPosNSigmaPion(0),
fTreeCascVarPosNSigmaProton(0),
fTreeCascVarBachNSigmaPion(0),
fTreeCascVarBachNSigmaKaon(0),

fTreeCascVarNegTOFNSigmaPion(0),
fTreeCascVarNegTOFNSigmaProton(0),
fTreeCascVarPosTOFNSigmaPion(0),
fTreeCascVarPosTOFNSigmaProton(0),
fTreeCascVarBachTOFNSigmaPion(0),
fTreeCascVarBachTOFNSigmaKaon(0),

fTreeCascVarNegITSNSigmaPion(0),
fTreeCascVarNegITSNSigmaProton(0),
fTreeCascVarPosITSNSigmaPion(0),
fTreeCascVarPosITSNSigmaProton(0),
fTreeCascVarBachITSNSigmaPion(0),
fTreeCascVarBachITSNSigmaKaon(0),

fTreeCascVarChiSquareV0(1e+3),
fTreeCascVarChiSquareCascade(1e+3),

fTreeCascVarBachDCAPVSigmaX2(0),
fTreeCascVarBachDCAPVSigmaY2(0),
fTreeCascVarBachDCAPVSigmaZ2(0),
fTreeCascVarPosDCAPVSigmaX2(0),
fTreeCascVarPosDCAPVSigmaY2(0),
fTreeCascVarPosDCAPVSigmaZ2(0),
fTreeCascVarNegDCAPVSigmaX2(0),
fTreeCascVarNegDCAPVSigmaY2(0),
fTreeCascVarNegDCAPVSigmaZ2(0),

fTreeCascVarPosPIDForTracking(-1),
fTreeCascVarNegPIDForTracking(-1),
fTreeCascVarBachPIDForTracking(-1),
fTreeCascVarNegInnerP(-1),
fTreeCascVarPosInnerP(-1),
fTreeCascVarBachInnerP(-1),
fTreeCascVarNegdEdx(-1),
fTreeCascVarPosdEdx(-1),
fTreeCascVarBachdEdx(-1),

fTreeCascVarNegTrackStatus(0), //!
fTreeCascVarPosTrackStatus(0), //!
fTreeCascVarBachTrackStatus(0), //!
fTreeCascVarNegDCAz(-1),
fTreeCascVarPosDCAz(-1),
fTreeCascVarBachDCAz(-1),
//fTreeCascVarPosTotMom(-1),
//fTreeCascVarNegTotMom(-1),
//fTreeCascVarBachTotMom(-1),

fTreeCascVarNegPx(0),
fTreeCascVarNegPy(0),
fTreeCascVarNegPz(0),
fTreeCascVarPosPx(0),
fTreeCascVarPosPy(0),
fTreeCascVarPosPz(0),
fTreeCascVarBachPx(0),
fTreeCascVarBachPy(0),
fTreeCascVarBachPz(0),
fTreeCascVarV0DecayX(0),
fTreeCascVarV0DecayY(0),
fTreeCascVarV0DecayZ(0),
fTreeCascVarCascadeDecayX(0),
fTreeCascVarCascadeDecayY(0),
fTreeCascVarCascadeDecayZ(0),
fTreeCascVarV0Lifetime(0),
//Track Labels (check for duplicates, etc)
fTreeCascVarNegIndex(0),
fTreeCascVarPosIndex(0),
fTreeCascVarBachIndex(0),
//Event Number (check same-event index mixups)
fTreeCascVarEventNumber(0),

fTreeCascVarPosITSClusters0(0),
fTreeCascVarPosITSClusters1(0),
fTreeCascVarPosITSClusters2(0),
fTreeCascVarPosITSClusters3(0),
fTreeCascVarPosITSClusters4(0),
fTreeCascVarPosITSClusters5(0),

fTreeCascVarNegITSClusters0(0),
fTreeCascVarNegITSClusters1(0),
fTreeCascVarNegITSClusters2(0),
fTreeCascVarNegITSClusters3(0),
fTreeCascVarNegITSClusters4(0),
fTreeCascVarNegITSClusters5(0),

fTreeCascVarBachITSClusters0(0),
fTreeCascVarBachITSClusters1(0),
fTreeCascVarBachITSClusters2(0),
fTreeCascVarBachITSClusters3(0),
fTreeCascVarBachITSClusters4(0),
fTreeCascVarBachITSClusters5(0),

fTreeCascVarPosITSSharedClusters0(0),
fTreeCascVarPosITSSharedClusters1(0),
fTreeCascVarPosITSSharedClusters2(0),
fTreeCascVarPosITSSharedClusters3(0),
fTreeCascVarPosITSSharedClusters4(0),
fTreeCascVarPosITSSharedClusters5(0),

fTreeCascVarNegITSSharedClusters0(0),
fTreeCascVarNegITSSharedClusters1(0),
fTreeCascVarNegITSSharedClusters2(0),
fTreeCascVarNegITSSharedClusters3(0),
fTreeCascVarNegITSSharedClusters4(0),
fTreeCascVarNegITSSharedClusters5(0),

fTreeCascVarBachITSSharedClusters0(0),
fTreeCascVarBachITSSharedClusters1(0),
fTreeCascVarBachITSSharedClusters2(0),
fTreeCascVarBachITSSharedClusters3(0),
fTreeCascVarBachITSSharedClusters4(0),
fTreeCascVarBachITSSharedClusters5(0),

//Variables for debugging the invariant mass bump
//Full momentum information
fTreeCascVarNegTOFExpTDiff(99999),
fTreeCascVarPosTOFExpTDiff(99999),
fTreeCascVarBachTOFExpTDiff(99999),
fTreeCascVarNegTOFSignal(99999),
fTreeCascVarPosTOFSignal(99999),
fTreeCascVarBachTOFSignal(99999),
fTreeCascVarNegTOFBCid(99999),
fTreeCascVarPosTOFBCid(99999),
fTreeCascVarBachTOFBCid(99999),
fTreeCascVarAmplitudeV0A(-1.),
fTreeCascVarAmplitudeV0C(-1.),
fTreeCascVarClosestNonEmptyBC(-1),

fTreeCascVarCentrality(0),
fTreeCascVarMVPileupFlag(kFALSE),
fTreeCascVarOOBPileupFlag(kFALSE),
//Kink tagging
fTreeCascVarBachIsKink(kFALSE),
fTreeCascVarPosIsKink(kFALSE),
fTreeCascVarNegIsKink(kFALSE),
fTreeCascVarIsCowboy(kFALSE),
fTreeCascVarCowboyness(-2),
fTreeCascVarIsCascadeCowboy(kFALSE),
fTreeCascVarCascadeCowboyness(-2),
fkSelectCharge(0),
fTreeCascVarPrimVertexX(0),
fTreeCascVarPrimVertexY(0),
fTreeCascVarPrimVertexZ(0),
fTreeCascVarBachTrack(0x0),
fTreeCascVarPosTrack(0x0),
fTreeCascVarNegTrack(0x0),
fTreeCascVarCosineProtonRestV0RestOmegaMinus(-100),
fTreeCascVarCosineProtonRestV0RestOmegaPlus(-100),
fTreeCascVarDzetaFromMomentaOmegaMinus(1e8),
fTreeCascVarMagneticField(0.0),
//extra
lEffMassXi(0.),
lNegTrackSign(1),
lPosTrackSign(1),
lBachTrackSign(1),
lEvSelCode(100),
lMinV0PA(0.0),
lMaxV0CosPA(1.0),
lMinV0CosPA(0.9),
lMaxCascCosPA(1.0),
lMinCascPA(0.0),
lMaxBBCosPA(1.0),
lMinBBPA(0.0),
dca(0.0),
fTreeCascVarMassXiScratch(0),
fTreeCascVarMassOmegaScratch(0),
fTreeCascVarBackgroundMassOmega(0),
fTreeCascVarBackgroundMassXi(0),
fTreeVariableAlphaXi(0),
fTreeVariablePtArmXi(0),
fkEsdTrackMultiplicity(0),
fkEsdEventTriggerWord(0),
fTreeCascVarBachP(0),
fTreeCascVarBachPt(0),
fTreeCascVarNegP(0),
fTreeCascVarNegPt(0),
fTreeCascVarPosP(0),
fTreeCascVarPosPt(0),
fkMomentaNeeded(kFALSE),
fkITSPID(kFALSE),
fkMassesAndBackgroundFromScratch(kFALSE),
fkDiscreteSymmetryInfo(kFALSE),
fkMultiplicityInfo(kFALSE),
fkExtraInfoCascadeKinematics(kFALSE),
fRunNum(0),
fNtracks(0),

//Histos
fHistEventCounter(0),
fHistEventCounterDifferential(0),
fHistCentrality(0),
fHistCentralityFast(0),
fHistZVertexEvent(0)

{
    
    //Re-vertex: Will only apply for cascade candidates
    
    fkDiscreteSymmetryInfo = lDiscreteSymmetryInfo;
    fkMassesAndBackgroundFromScratch = lMassesAndBackgroundFromScratch;
    fkMomentaNeeded = lMomentaNeeded;
    fkExtraInfoCascadeKinematics = lExtraInfoCascadeKinematics;
    
    fkRunVertexers = lRunVertexers;
    fkUseLightVertexer = lUseLightVertexer;
    
    fV0VertexerSels[0] =  33.  ;  // max allowed chi2
    fV0VertexerSels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
    fV0VertexerSels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
    fV0VertexerSels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
    fV0VertexerSels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
    fV0VertexerSels[5] =   1.0 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
    fV0VertexerSels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
    
    fCascadeVertexerSels[0] = lCascaderMaxChi2  ;  // max allowed chi2 (same as PDC07) (33)
    fCascadeVertexerSels[1] = lCascaderV0MinImpactParam ;  // min allowed V0 impact parameter (0.05)                   (PDC07 : 0.05   / LHC09a4 : 0.025 )
    fCascadeVertexerSels[2] = lCascaderV0MassWindow;  // "window" around the Lambda mass  (0.010)                  (PDC07 : 0.008  / LHC09a4 : 0.010 )
    fCascadeVertexerSels[3] = lCascaderBachMinImpactParam;  // min allowed bachelor's impact parameter (0.03)            (PDC07 : 0.035  / LHC09a4 : 0.025 )
    fCascadeVertexerSels[4] = lCascaderMaxDCAV0andBach;  // max allowed DCA between the V0 and the bachelor (2.0)    (PDC07 : 0.1    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[5] = lCascaderMinCosAngle;  // min allowed cosine of the cascade pointing angle (0.95)  (PDC07 : 0.9985 / LHC09a4 : 0.998 )
    fCascadeVertexerSels[6] = lCascaderMinRadius;  // min radius of the fiducial volume (0.4)                 (PDC07 : 0.9    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[7] = lCascaderMaxRadius;  // max radius of the fiducial volume (100.)                 (PDC07 : 100    / LHC09a4 : 100   )
    
    
    /*fCascadeVertexerSels[0] = 33.;  // max allowed chi2 (same as PDC07) (33)
     fCascadeVertexerSels[1] = 0.05 ;  // min allowed V0 impact parameter (0.05)                   (PDC07 : 0.05   / LHC09a4 : 0.025 )
     fCascadeVertexerSels[2] = 0.010;  // "window" around the Lambda mass  (0.010)                  (PDC07 : 0.008  / LHC09a4 : 0.010 )
     fCascadeVertexerSels[3] = 0.03;  // min allowed bachelor's impact parameter (0.03)            (PDC07 : 0.035  / LHC09a4 : 0.025 )
     fCascadeVertexerSels[4] = 2.0;  // max allowed DCA between the V0 and the bachelor (2.0)    (PDC07 : 0.1    / LHC09a4 : 0.2   )
     fCascadeVertexerSels[5] = 0.95;  // min allowed cosine of the cascade pointing angle (0.95)  (PDC07 : 0.9985 / LHC09a4 : 0.998 )
     fCascadeVertexerSels[6] = 0.4;  // min radius of the fiducial volume (0.4)                 (PDC07 : 0.9    / LHC09a4 : 0.2   )
     fCascadeVertexerSels[7] = 100.;  // max radius of the fiducial volume (100.)                 (PDC07 : 100    / LHC09a4 : 100   )
     */
    
    //[0]+[1]*TMath::Exp([2]*x)+[3]*TMath::Exp([4]*x)
    fLambdaMassMean[0]=1.116; //standard fixed
    fLambdaMassMean[1]=0.0;
    fLambdaMassMean[2]=0.0;
    fLambdaMassMean[3]=0.0;
    fLambdaMassMean[4]=0.0;
    
    //[0]+[1]*x+[2]*TMath::Exp([3]*x)
    fLambdaMassSigma[0]=0.002; //standard at roughly the integ val
    fLambdaMassSigma[1]=0.0;
    fLambdaMassSigma[2]=0.0;
    fLambdaMassSigma[3]=0.0;
    
    fkSaveEventTree    = lSaveEventTree;
    fkSaveV0Tree       = lSaveV0Tree;
    fkSaveCascadeTree  = lSaveCascadeTree;
    
    
    DefineOutput(1, TTree::Class()); // Cascade Tree output (as event)
    // DefineOutput(2, TTree::Class()); //Cascade Tree separate, heavy.
    
    
    
    //Special Debug Options (more to be added as needed)
    // A - Study Wrong PID for tracking bug
    // B - Study invariant mass *B*ump
    // C - Study OOB pileup in pp 2016 data
    // S - Add sandbox mode information, please
    if ( lExtraOptions.Contains("A") ) fkDebugWrongPIDForTracking = kTRUE;
    if ( lExtraOptions.Contains("B") ) fkDebugBump                = kTRUE;
    if ( lExtraOptions.Contains("C") ) fkDebugOOBPileup           = kTRUE;
    if ( lExtraOptions.Contains("S") ) fkSandboxMode              = kTRUE;
}


AliAnalysisTaskStrangeCascadesDiscrete::~AliAnalysisTaskStrangeCascadesDiscrete()
{
    
    
    delete[] fListHist;
    delete[] fListK0Short;
    delete[] fListLambda;
    delete[] fListAntiLambda;
    delete[] fListXiMinus;
    delete[] fListXiPlus;
    delete[] fListOmegaMinus;
    delete[] fListOmegaPlus;
    delete[] fEventList;
    delete[] fTreeEvent;
    //  delete[] fTreeCascade;
    //  delete[] Cascade_Track;
    //  delete[] Cascade_Event;
    // delete[] fTreeCascadeAsEvent;
    
    //not sure if needed
    delete[] fPIDResponse;
    delete[] fESDtrackCuts;
    delete[] fESDtrackCutsITSsa2010;
    delete[] fESDtrackCutsGlobal2015;
    delete[] fUtils;
    delete[] fRand;
    
    delete[] fTreeVariablePosTrack;
    delete[] fTreeVariableNegTrack;
    
    delete[] fTreeCascVarBachTrack;
    delete[] fTreeCascVarPosTrack;
    delete[] fTreeCascVarNegTrack;
    //delete histograms
    delete[] fHistEventCounter;
    delete[] fHistEventCounterDifferential;
    delete[] fHistCentrality;
    delete[] fHistCentralityFast;
    delete[] fHistZVertexEvent;
    
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    
    //Destroy output objects if present
    /*   if (fListHist) {
     delete fListHist;
     //  fListHist = 0x0;
     }
     if (fListK0Short) {
     delete fListK0Short;
     //  fListK0Short = 0x0;
     }
     if (fListLambda) {
     delete fListLambda;
     //   fListLambda = 0x0;
     }
     if (fListAntiLambda) {
     delete fListAntiLambda;
     //   fListAntiLambda = 0x0;
     }
     if (fListXiMinus) {
     delete fListXiMinus;
     //   fListXiMinus = 0x0;
     }
     if (fListXiPlus) {
     delete fListXiPlus;
     //   fListXiPlus = 0x0;
     }
     if (fListOmegaMinus) {
     delete fListOmegaMinus;
     //   fListOmegaMinus = 0x0;
     }
     if (fListOmegaPlus) {
     delete fListOmegaPlus;
     //  fListOmegaPlus = 0x0;
     }
     //  if (fTreeEvent) {
     delete fTreeEvent;
     //   fTreeEvent = 0x0;
     //  }
     //  if (fTreeV0) {
     delete fTreeV0;
     //  fTreeV0 = 0x0;
     //  }
     //  if (fTreeCascade) {
     delete fTreeCascade;
     //  fTreeCascade = 0x0;
     //   }
     if (fUtils) {
     delete fUtils;
     //  fUtils = 0x0;
     }
     if (fRand) {
     delete fRand;
     //  fRand = 0x0;
     }   */
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::UserCreateOutputObjects()
{
    
    //fNumberEventTree serves as a test, in order to figure out whether the programm does fill the trees or just ignores them
    // Int_t number(0);
    //  fNumberEventTree =new TTree("fNumberEventTree", "Number of events");
    //  fNumberEventTree->Branch("number", &number, "number/I");
    //------------------------------------------------
    // fTreeEvent: EbyE information
    //------------------------------------------------
    // if(fkSaveEventTree){
    
    fTreeEvent = new TTree("fTreeEvent","Event");
    fTreeEvent->Branch("fCentrality",&fCentrality,"fCentrality/F");
    fTreeEvent->Branch("fMVPileupFlag",&fMVPileupFlag,"fMVPileupFlag/O");
    fTreeEvent->Branch("fTreeVariablePrimVertexX", &fTreeVariablePrimVertexX, "fTreeVariablePrimVertexX/F" );
    fTreeEvent->Branch("fTreeVariablePrimVertexY", &fTreeVariablePrimVertexY, "fTreeVariablePrimVertexY/F" );
    fTreeEvent->Branch("fTreeVariablePrimVertexZ", &fTreeVariablePrimVertexZ, "fTreeVariablePrimVertexZ/F" );
    fTreeEvent->Branch("fkEsdTrackMultiplicity", &fkEsdTrackMultiplicity, "fkEsdTrackMultiplicity/I");
    fTreeEvent->Branch("fkEsdEventTriggerWord", &fkEsdEventTriggerWord, "fkEsdEventTriggerWord/l"); // l  for ULong64_t, i  for UInt_t
    //
    if ( fkDebugOOBPileup ){
        fTreeEvent->Branch("fOOBPileupFlag",&fOOBPileupFlag,"fOOBPileupFlag/O");
        fTreeEvent->Branch("fNTOFClusters",&fNTOFClusters,"fNTOFClusters/I");
        fTreeEvent->Branch("fNTOFMatches",&fNTOFMatches,"fNTOFMatches/I");
        fTreeEvent->Branch("fNTracksITSsa2010",&fNTracksITSsa2010,"fNTracksITSsa2010/I");
        fTreeEvent->Branch("fNTracksGlobal2015",&fNTracksGlobal2015,"fNTracksGlobal2015/I");
        fTreeEvent->Branch("fNTracksGlobal2015TriggerPP",&fNTracksGlobal2015TriggerPP,"fNTracksGlobal2015TriggerPP/I");
        fTreeEvent->Branch("fAmplitudeV0A",&fAmplitudeV0A,"fAmplitudeV0A/F");
        fTreeEvent->Branch("fAmplitudeV0C",&fAmplitudeV0C,"fAmplitudeV0C/F");
        fTreeEvent->Branch("fClosestNonEmptyBC",&fClosestNonEmptyBC,"fClosestNonEmptyBC/I");
    }
    
    
    //------------------------------------------------
    // fTreeCascade Branch definitions - Cascade Tree
    //------------------------------------------------
    //  if(fkSaveCascadeTree){
    //Create Cascade output tree
    fTreeCascade = new TTree("fTreeCascade","CascadeCandidates");
    //-----------BASIC-INFO---------------------------
    fTreeCascade->Branch("fTreeCascVarCharge",&fTreeCascVarCharge,"fTreeCascVarCharge/I");
    fTreeCascade->Branch("fTreeCascVarMassAsXi",&fTreeCascVarMassAsXi,"fTreeCascVarMassAsXi/F");
    fTreeCascade->Branch("fTreeCascVarMassAsOmega",&fTreeCascVarMassAsOmega,"fTreeCascVarMassAsOmega/F");
    fTreeCascade->Branch("fTreeCascVarPt",&fTreeCascVarPt,"fTreeCascVarPt/F");
    //------------------------------------------------------------added on 12.06.19
    if ( fkExtraInfoCascadeKinematics)
    {
        fTreeCascade->Branch("fTreeCascVarP",&fTreeCascVarP,"fTreeCascVarP/F");
        //   fTreeCascade->Branch("fTreeCascVarPz",&fTreeCascVarPz,"fTreeCascVarPz/F");
        
    }
    fTreeCascade->Branch("fTreeCascVarRapXi",&fTreeCascVarRapXi,"fTreeCascVarRapXi/F");
    fTreeCascade->Branch("fTreeCascVarRapOmega",&fTreeCascVarRapOmega,"fTreeCascVarRapOmega/F");
    fTreeCascade->Branch("fTreeCascVarNegEta",&fTreeCascVarNegEta,"fTreeCascVarNegEta/F");
    fTreeCascade->Branch("fTreeCascVarPosEta",&fTreeCascVarPosEta,"fTreeCascVarPosEta/F");
    fTreeCascade->Branch("fTreeCascVarBachEta",&fTreeCascVarBachEta,"fTreeCascVarBachEta/F");
    //-----------INFO-FOR-CUTS------------------------
    fTreeCascade->Branch("fTreeCascVarDCACascDaughters",&fTreeCascVarDCACascDaughters,"fTreeCascVarDCACascDaughters/F");
    fTreeCascade->Branch("fTreeCascVarDCABachToPrimVtx",&fTreeCascVarDCABachToPrimVtx,"fTreeCascVarDCABachToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarDCABachToPrimVtxZ",&fTreeCascVarDCABachToPrimVtxZ,"fTreeCascVarDCABachToPrimVtxZ/F");
    fTreeCascade->Branch("fTreeCascVarDCAV0Daughters",&fTreeCascVarDCAV0Daughters,"fTreeCascVarDCAV0Daughters/F");
    fTreeCascade->Branch("fTreeCascVarDCAV0ToPrimVtx",&fTreeCascVarDCAV0ToPrimVtx,"fTreeCascVarDCAV0ToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarDCAPosToPrimVtx",&fTreeCascVarDCAPosToPrimVtx,"fTreeCascVarDCAPosToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarDCAPosToPrimVtxZ",&fTreeCascVarDCAPosToPrimVtxZ,"fTreeCascVarDCAPosToPrimVtxZ/F");
    fTreeCascade->Branch("fTreeCascVarDCANegToPrimVtx",&fTreeCascVarDCANegToPrimVtx,"fTreeCascVarDCANegToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarDCANegToPrimVtxZ",&fTreeCascVarDCANegToPrimVtxZ,"fTreeCascVarDCANegToPrimVtxZ/F");
    fTreeCascade->Branch("fTreeCascVarCascCosPointingAngle",&fTreeCascVarCascCosPointingAngle,"fTreeCascVarCascCosPointingAngle/F");
    fTreeCascade->Branch("fTreeCascVarCascDCAtoPVxy",&fTreeCascVarCascDCAtoPVxy,"fTreeCascVarCascDCAtoPVxy/F");
    fTreeCascade->Branch("fTreeCascVarCascDCAtoPVz",&fTreeCascVarCascDCAtoPVz,"fTreeCascVarCascDCAtoPVz/F");
    
    fTreeCascade->Branch("fTreeCascVarCascRadius",&fTreeCascVarCascRadius,"fTreeCascVarCascRadius/F");
    fTreeCascade->Branch("fTreeCascVarV0Mass",&fTreeCascVarV0Mass,"fTreeCascVarV0Mass/F");
    fTreeCascade->Branch("fTreeCascVarV0CosPointingAngle",&fTreeCascVarV0CosPointingAngle,"fTreeCascVarV0CosPointingAngle/F");
    fTreeCascade->Branch("fTreeCascVarV0CosPointingAngleSpecial",&fTreeCascVarV0CosPointingAngleSpecial,"fTreeCascVarV0CosPointingAngleSpecial/F");
    fTreeCascade->Branch("fTreeCascVarV0Radius",&fTreeCascVarV0Radius,"fTreeCascVarV0Radius/F");
    fTreeCascade->Branch("fTreeCascVarDCABachToBaryon",&fTreeCascVarDCABachToBaryon,"fTreeCascVarDCABachToBaryon/F");
    fTreeCascade->Branch("fTreeCascVarWrongCosPA",&fTreeCascVarWrongCosPA,"fTreeCascVarWrongCosPA/F");
    fTreeCascade->Branch("fTreeCascVarLeastNbrClusters",&fTreeCascVarLeastNbrClusters,"fTreeCascVarLeastNbrClusters/I");
    fTreeCascade->Branch("fTreeCascVarMaxChi2PerCluster",&fTreeCascVarMaxChi2PerCluster,"fTreeCascVarMaxChi2PerCluster/F");
    fTreeCascade->Branch("fTreeCascVarMinTrackLength",&fTreeCascVarMinTrackLength,"fTreeCascVarMinTrackLength/F");
    
    //-----------DECAY-LENGTH-INFO--------------------
    fTreeCascade->Branch("fTreeCascVarDistOverTotMom",&fTreeCascVarDistOverTotMom,"fTreeCascVarDistOverTotMom/F");
    //------------------------------------------------
    fTreeCascade->Branch("fTreeCascVarNegNSigmaPion",&fTreeCascVarNegNSigmaPion,"fTreeCascVarNegNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarNegNSigmaProton",&fTreeCascVarNegNSigmaProton,"fTreeCascVarNegNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarPosNSigmaPion",&fTreeCascVarPosNSigmaPion,"fTreeCascVarPosNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarPosNSigmaProton",&fTreeCascVarPosNSigmaProton,"fTreeCascVarPosNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarBachNSigmaPion",&fTreeCascVarBachNSigmaPion,"fTreeCascVarBachNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarBachNSigmaKaon",&fTreeCascVarBachNSigmaKaon,"fTreeCascVarBachNSigmaKaon/F");
    
    fTreeCascade->Branch("fTreeCascVarNegTOFNSigmaPion",&fTreeCascVarNegTOFNSigmaPion,"fTreeCascVarNegTOFNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarNegTOFNSigmaProton",&fTreeCascVarNegTOFNSigmaProton,"fTreeCascVarTOFNegNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarPosTOFNSigmaPion",&fTreeCascVarPosTOFNSigmaPion,"fTreeCascVarPosTOFNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarPosTOFNSigmaProton",&fTreeCascVarPosTOFNSigmaProton,"fTreeCascVarPosTOFNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarBachTOFNSigmaPion",&fTreeCascVarBachTOFNSigmaPion,"fTreeCascVarBachTOFNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarBachTOFNSigmaKaon",&fTreeCascVarBachTOFNSigmaKaon,"fTreeCascVarBachTOFNSigmaKaon/F");
    
    fTreeCascade->Branch("fTreeCascVarIsCowboy",&fTreeCascVarIsCowboy,"fTreeCascVarIsCowboy/O");
    fTreeCascade->Branch("fTreeCascVarIsCascadeCowboy",&fTreeCascVarIsCascadeCowboy,"fTreeCascVarIsCascadeCowboy/O");
    fTreeCascade->Branch("fTreeVariableAlphaXi", &fTreeVariableAlphaXi, "fTreeVariableAlphaXi/D"); //in float!
    fTreeCascade->Branch("fTreeVariablePtArmXi", &fTreeVariablePtArmXi, "fTreeVariablePtArmXi/D"); //in float!
    
    //------------------------------------------------
    // fTreeCascade->Branch("fTreeCascVarChiSquareV0",&fTreeCascVarChiSquareV0,"fTreeCascVarChiSquareV0/F");
    // fTreeCascade->Branch("fTreeCascVarChiSquareCascade",&fTreeCascVarChiSquareCascade,"fTreeCascVarChiSquareCascade/F");
    
    //------------------------------------------------------------added on 11.06.19
    if ( fkMomentaNeeded) {
        fTreeCascade->Branch("fTreeCascVarPosP", &fTreeCascVarPosP, "fTreeCascVarPosP/F");
        fTreeCascade->Branch("fTreeCascVarPosPt", &fTreeCascVarPosPt, "fTreeCascVarPosPt/F");
        fTreeCascade->Branch("fTreeCascVarPosPz", &fTreeCascVarPosPz, "fTreeCascVarPosPz/F");
        
        fTreeCascade->Branch("fTreeCascVarNegP", &fTreeCascVarNegP, "fTreeCascVarNegP/F");
        fTreeCascade->Branch("fTreeCascVarNegPt", &fTreeCascVarNegPt, "fTreeCascVarNegPt/F");
        fTreeCascade->Branch("fTreeCascVarNegPz", &fTreeCascVarNegPz, "fTreeCascVarNegPz/F");
        
        fTreeCascade->Branch("fTreeCascVarBachP", &fTreeCascVarBachP, "fTreeCascVarBachP/F");
        fTreeCascade->Branch("fTreeCascVarBachPt", &fTreeCascVarBachPt, "fTreeCascVarBachPt/F");
        fTreeCascade->Branch("fTreeCascVarBachPz", &fTreeCascVarBachPz, "fTreeCascVarBachPz/F");
        
    }
    
    
    if ( fkDiscreteSymmetryInfo) {
        //--------SYMMETRY PARAMETERS----------------------------
        fTreeCascade->Branch("fTreeCascVarCosineProtonRestV0RestOmegaMinus", &fTreeCascVarCosineProtonRestV0RestOmegaMinus, "fTreeCascVarCosineProtonRestV0RestOmegaMinus/D");
        fTreeCascade->Branch("fTreeCascVarCosineProtonRestV0RestOmegaPlus", &fTreeCascVarCosineProtonRestV0RestOmegaPlus, "fTreeCascVarCosineProtonRestV0RestOmegaPlus/D");
        fTreeCascade->Branch("fTreeCascVarDzetaFromMomentaOmegaMinus",&fTreeCascVarDzetaFromMomentaOmegaMinus,"fTreeCascVarDzetaFromMomentaOmegaMinus/D");
    }
    
    if ( fkITSPID) {
        fTreeCascade->Branch("fTreeCascVarNegITSNSigmaPion",&fTreeCascVarNegITSNSigmaPion,"fTreeCascVarNegITSNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarNegITSNSigmaProton",&fTreeCascVarNegITSNSigmaProton,"fTreeCascVarITSNegNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarPosITSNSigmaPion",&fTreeCascVarPosITSNSigmaPion,"fTreeCascVarPosITSNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarPosITSNSigmaProton",&fTreeCascVarPosITSNSigmaProton,"fTreeCascVarPosITSNSigmaProton/F");
        fTreeCascade->Branch("fTreeCascVarBachITSNSigmaPion",&fTreeCascVarBachITSNSigmaPion,"fTreeCascVarBachITSNSigmaPion/F");
        fTreeCascade->Branch("fTreeCascVarBachITSNSigmaKaon",&fTreeCascVarBachITSNSigmaKaon,"fTreeCascVarBachITSNSigmaKaon/F");
    }
    
    if ( fkMassesAndBackgroundFromScratch ) {
        fTreeCascade->Branch("fTreeCascVarV0MassLambda",&fTreeCascVarV0MassLambda,"fTreeCascVarV0MassLambda/F"); //from scratch
        fTreeCascade->Branch("fTreeCascVarV0MassAntiLambda",&fTreeCascVarV0MassAntiLambda,"fTreeCascVarV0MassAntiLambda/F"); //from scratch
        fTreeCascade->Branch("fTreeCascVarMassXiScratch", &fTreeCascVarMassXiScratch, "fTreeCascVarMassXiScratch/F");
        fTreeCascade->Branch("fTreeCascVarMassOmegaScratch", &fTreeCascVarMassOmegaScratch, "fTreeCascVarMassOmegaScratch/F");
        fTreeCascade->Branch( "fTreeCascVarBackgroundMassOmega", &fTreeCascVarBackgroundMassOmega, "fTreeCascVarBackgroundMassOmega/F");
        fTreeCascade->Branch("fTreeCascVarBackgroundMassXi", &fTreeCascVarBackgroundMassXi,"fTreeCascVarBackgroundMassXi/F");
    }
    
    if ( fkMultiplicityInfo) {
        fTreeCascade->Branch("fTreeCascVarCentrality",&fTreeCascVarCentrality,"fTreeCascVarCentrality/F");
        fTreeCascade->Branch("fTreeCascVarMVPileupFlag",&fTreeCascVarMVPileupFlag,"fTreeCascVarMVPileupFlag/O");
    }
    
    //------------------------------------------------
    if ( fkDebugWrongPIDForTracking ){
        fTreeCascade->Branch("fTreeCascVarPosPIDForTracking",&fTreeCascVarPosPIDForTracking,"fTreeCascVarPosPIDForTracking/I");
        fTreeCascade->Branch("fTreeCascVarNegPIDForTracking",&fTreeCascVarNegPIDForTracking,"fTreeCascVarNegPIDForTracking/I");
        fTreeCascade->Branch("fTreeCascVarBachPIDForTracking",&fTreeCascVarBachPIDForTracking,"fTreeCascVarBachPIDForTracking/I");
        fTreeCascade->Branch("fTreeCascVarPosdEdx",&fTreeCascVarPosdEdx,"fTreeCascVarPosdEdx/F");
        fTreeCascade->Branch("fTreeCascVarNegdEdx",&fTreeCascVarNegdEdx,"fTreeCascVarNegdEdx/F");
        fTreeCascade->Branch("fTreeCascVarBachdEdx",&fTreeCascVarBachdEdx,"fTreeCascVarBachdEdx/F");
        fTreeCascade->Branch("fTreeCascVarPosInnerP",&fTreeCascVarPosInnerP,"fTreeCascVarPosInnerP/F");
        fTreeCascade->Branch("fTreeCascVarNegInnerP",&fTreeCascVarNegInnerP,"fTreeCascVarNegInnerP/F");
        fTreeCascade->Branch("fTreeCascVarBachInnerP",&fTreeCascVarBachInnerP,"fTreeCascVarBachInnerP/F");
        fTreeCascade->Branch("fTreeCascVarNegTrackStatus",&fTreeCascVarNegTrackStatus,"fTreeCascVarNegTrackStatus/l");
        fTreeCascade->Branch("fTreeCascVarPosTrackStatus",&fTreeCascVarPosTrackStatus,"fTreeCascVarPosTrackStatus/l");
        fTreeCascade->Branch("fTreeCascVarBachTrackStatus",&fTreeCascVarBachTrackStatus,"fTreeCascVarBachTrackStatus/l");
        fTreeCascade->Branch("fTreeCascVarNegDCAz",&fTreeCascVarNegDCAz,"fTreeCascVarNegDCAz/F");
        fTreeCascade->Branch("fTreeCascVarPosDCAz",&fTreeCascVarPosDCAz,"fTreeCascVarPosDCAz/F");
        fTreeCascade->Branch("fTreeCascVarBachDCAz",&fTreeCascVarBachDCAz,"fTreeCascVarBachDCAz/F");
    }
    //------------------------------------------------
    if ( fkDebugBump ){
        //Variables for debugging the invariant mass bump
        //Full momentum information
        fTreeCascade->Branch("fTreeCascVarPosPx",&fTreeCascVarPosPx,"fTreeCascVarPosPx/F");
        fTreeCascade->Branch("fTreeCascVarPosPy",&fTreeCascVarPosPy,"fTreeCascVarPosPy/F");
        fTreeCascade->Branch("fTreeCascVarPosPz",&fTreeCascVarPosPz,"fTreeCascVarPosPz/F");
        fTreeCascade->Branch("fTreeCascVarNegPx",&fTreeCascVarNegPx,"fTreeCascVarNegPx/F");
        fTreeCascade->Branch("fTreeCascVarNegPy",&fTreeCascVarNegPy,"fTreeCascVarNegPy/F");
        fTreeCascade->Branch("fTreeCascVarNegPz",&fTreeCascVarNegPz,"fTreeCascVarNegPz/F");
        fTreeCascade->Branch("fTreeCascVarBachPx",&fTreeCascVarBachPx,"fTreeCascVarBachPx/F");
        fTreeCascade->Branch("fTreeCascVarBachPy",&fTreeCascVarBachPy,"fTreeCascVarBachPy/F");
        fTreeCascade->Branch("fTreeCascVarBachPz",&fTreeCascVarBachPz,"fTreeCascVarBachPz/F");
        // Decay positions
        fTreeCascade->Branch("fTreeCascVarV0DecayX",&fTreeCascVarV0DecayX,"fTreeCascVarV0DecayX/F");
        fTreeCascade->Branch("fTreeCascVarV0DecayY",&fTreeCascVarV0DecayY,"fTreeCascVarV0DecayY/F");
        fTreeCascade->Branch("fTreeCascVarV0DecayZ",&fTreeCascVarV0DecayZ,"fTreeCascVarV0DecayZ/F");
        fTreeCascade->Branch("fTreeCascVarCascadeDecayX",&fTreeCascVarCascadeDecayX,"fTreeCascVarCascadeDecayX/F");
        fTreeCascade->Branch("fTreeCascVarCascadeDecayY",&fTreeCascVarCascadeDecayY,"fTreeCascVarCascadeDecayY/F");
        fTreeCascade->Branch("fTreeCascVarCascadeDecayZ",&fTreeCascVarCascadeDecayZ,"fTreeCascVarCascadeDecayZ/F");
        
        fTreeCascade->Branch("fTreeCascVarV0Lifetime",&fTreeCascVarV0Lifetime,"fTreeCascVarV0Lifetime/F");
        
        //Track Labels (check for duplicates, etc)
        fTreeCascade->Branch("fTreeCascVarNegIndex",&fTreeCascVarNegIndex,"fTreeCascVarNegIndex/I");
        fTreeCascade->Branch("fTreeCascVarPosIndex",&fTreeCascVarPosIndex,"fTreeCascVarPosIndex/I");
        fTreeCascade->Branch("fTreeCascVarBachIndex",&fTreeCascVarBachIndex,"fTreeCascVarBachIndex/I");
        //Event Number (check same-event index mixups)
        fTreeCascade->Branch("fTreeCascVarEventNumber",&fTreeCascVarEventNumber,"fTreeCascVarEventNumber/l");
    }
    if ( fkSandboxMode ){
        //Full track info for DCA minim optimization
        fTreeCascade->Branch("fTreeCascVarBachTrack", &fTreeCascVarBachTrack,16000,99);
        fTreeCascade->Branch("fTreeCascVarPosTrack", &fTreeCascVarPosTrack,16000,99);
        fTreeCascade->Branch("fTreeCascVarNegTrack", &fTreeCascVarNegTrack,16000,99);
        
        //for sandbox mode
        fTreeCascade->Branch("fTreeCascVarMagneticField",&fTreeCascVarMagneticField,"fTreeCascVarMagneticField/F");
        
        //Cascade decay position calculation metrics
        fTreeCascade->Branch("fTreeCascVarPrimVertexX",&fTreeCascVarPrimVertexX,"fTreeCascVarPrimVertexX/F");
        fTreeCascade->Branch("fTreeCascVarPrimVertexY",&fTreeCascVarPrimVertexY,"fTreeCascVarPrimVertexY/F");
        fTreeCascade->Branch("fTreeCascVarPrimVertexZ",&fTreeCascVarPrimVertexZ,"fTreeCascVarPrimVertexZ/F");
    }
    if ( fkDebugOOBPileup ) {
        fTreeCascade->Branch("fTreeCascVarPosITSClusters0",&fTreeCascVarPosITSClusters0,"fTreeCascVarPosITSClusters0/O");
        fTreeCascade->Branch("fTreeCascVarPosITSClusters1",&fTreeCascVarPosITSClusters1,"fTreeCascVarPosITSClusters1/O");
        fTreeCascade->Branch("fTreeCascVarPosITSClusters2",&fTreeCascVarPosITSClusters2,"fTreeCascVarPosITSClusters2/O");
        fTreeCascade->Branch("fTreeCascVarPosITSClusters3",&fTreeCascVarPosITSClusters3,"fTreeCascVarPosITSClusters3/O");
        fTreeCascade->Branch("fTreeCascVarPosITSClusters4",&fTreeCascVarPosITSClusters4,"fTreeCascVarPosITSClusters4/O");
        fTreeCascade->Branch("fTreeCascVarPosITSClusters5",&fTreeCascVarPosITSClusters5,"fTreeCascVarPosITSClusters5/O");
        
        fTreeCascade->Branch("fTreeCascVarNegITSClusters0",&fTreeCascVarNegITSClusters0,"fTreeCascVarNegITSClusters0/O");
        fTreeCascade->Branch("fTreeCascVarNegITSClusters1",&fTreeCascVarNegITSClusters1,"fTreeCascVarNegITSClusters1/O");
        fTreeCascade->Branch("fTreeCascVarNegITSClusters2",&fTreeCascVarNegITSClusters2,"fTreeCascVarNegITSClusters2/O");
        fTreeCascade->Branch("fTreeCascVarNegITSClusters3",&fTreeCascVarNegITSClusters3,"fTreeCascVarNegITSClusters3/O");
        fTreeCascade->Branch("fTreeCascVarNegITSClusters4",&fTreeCascVarNegITSClusters4,"fTreeCascVarNegITSClusters4/O");
        fTreeCascade->Branch("fTreeCascVarNegITSClusters5",&fTreeCascVarNegITSClusters5,"fTreeCascVarNegITSClusters5/O");
        
        fTreeCascade->Branch("fTreeCascVarBachITSClusters0",&fTreeCascVarBachITSClusters0,"fTreeCascVarBachITSClusters0/O");
        fTreeCascade->Branch("fTreeCascVarBachITSClusters1",&fTreeCascVarBachITSClusters1,"fTreeCascVarBachITSClusters1/O");
        fTreeCascade->Branch("fTreeCascVarBachITSClusters2",&fTreeCascVarBachITSClusters2,"fTreeCascVarBachITSClusters2/O");
        fTreeCascade->Branch("fTreeCascVarBachITSClusters3",&fTreeCascVarBachITSClusters3,"fTreeCascVarBachITSClusters3/O");
        fTreeCascade->Branch("fTreeCascVarBachITSClusters4",&fTreeCascVarBachITSClusters4,"fTreeCascVarBachITSClusters4/O");
        fTreeCascade->Branch("fTreeCascVarBachITSClusters5",&fTreeCascVarBachITSClusters5,"fTreeCascVarBachITSClusters5/O");
        
        fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters0",&fTreeCascVarPosITSSharedClusters0,"fTreeCascVarPosITSSharedClusters0/O");
        fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters1",&fTreeCascVarPosITSSharedClusters1,"fTreeCascVarPosITSSharedClusters1/O");
        fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters2",&fTreeCascVarPosITSSharedClusters2,"fTreeCascVarPosITSSharedClusters2/O");
        fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters3",&fTreeCascVarPosITSSharedClusters3,"fTreeCascVarPosITSSharedClusters3/O");
        fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters4",&fTreeCascVarPosITSSharedClusters4,"fTreeCascVarPosITSSharedClusters4/O");
        fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters5",&fTreeCascVarPosITSSharedClusters5,"fTreeCascVarPosITSSharedClusters5/O");
        
        fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters0",&fTreeCascVarNegITSSharedClusters0,"fTreeCascVarNegITSSharedClusters0/O");
        fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters1",&fTreeCascVarNegITSSharedClusters1,"fTreeCascVarNegITSSharedClusters1/O");
        fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters2",&fTreeCascVarNegITSSharedClusters2,"fTreeCascVarNegITSSharedClusters2/O");
        fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters3",&fTreeCascVarNegITSSharedClusters3,"fTreeCascVarNegITSSharedClusters3/O");
        fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters4",&fTreeCascVarNegITSSharedClusters4,"fTreeCascVarNegITSSharedClusters4/O");
        fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters5",&fTreeCascVarNegITSSharedClusters5,"fTreeCascVarNegITSSharedClusters5/O");
        
        fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters0",&fTreeCascVarBachITSSharedClusters0,"fTreeCascVarBachITSSharedClusters0/O");
        fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters1",&fTreeCascVarBachITSSharedClusters1,"fTreeCascVarBachITSSharedClusters1/O");
        fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters2",&fTreeCascVarBachITSSharedClusters2,"fTreeCascVarBachITSSharedClusters2/O");
        fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters3",&fTreeCascVarBachITSSharedClusters3,"fTreeCascVarBachITSSharedClusters3/O");
        fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters4",&fTreeCascVarBachITSSharedClusters4,"fTreeCascVarBachITSSharedClusters4/O");
        fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters5",&fTreeCascVarBachITSSharedClusters5,"fTreeCascVarBachITSSharedClusters5/O");
        
        fTreeCascade->Branch("fTreeCascVarNegTOFExpTDiff",&fTreeCascVarNegTOFExpTDiff,"fTreeCascVarNegTOFExpTDiff/F");
        fTreeCascade->Branch("fTreeCascVarPosTOFExpTDiff",&fTreeCascVarPosTOFExpTDiff,"fTreeCascVarPosTOFExpTDiff/F");
        fTreeCascade->Branch("fTreeCascVarBachTOFExpTDiff",&fTreeCascVarBachTOFExpTDiff,"fTreeCascVarBachTOFExpTDiff/F");
        fTreeCascade->Branch("fTreeCascVarNegTOFSignal",&fTreeCascVarNegTOFSignal,"fTreeCascVarNegTOFSignal/F");
        fTreeCascade->Branch("fTreeCascVarPosTOFSignal",&fTreeCascVarPosTOFSignal,"fTreeCascVarPosTOFSignal/F");
        fTreeCascade->Branch("fTreeCascVarBachTOFSignal",&fTreeCascVarBachTOFSignal,"fTreeCascVarBachTOFSignal/F");
        fTreeCascade->Branch("fTreeCascVarNegTOFBCid",&fTreeCascVarNegTOFBCid,"fTreeCascVarNegTOFBCid/I");
        fTreeCascade->Branch("fTreeCascVarPosTOFBCid",&fTreeCascVarPosTOFBCid,"fTreeCascVarPosTOFBCid/I");
        fTreeCascade->Branch("fTreeCascVarBachTOFBCid",&fTreeCascVarBachTOFBCid,"fTreeCascVarBachTOFBCid/I");
        // Event info
        fTreeCascade->Branch("fTreeCascVarOOBPileupFlag",&fTreeCascVarOOBPileupFlag,"fTreeCascVarOOBPileupFlag/O");
        fTreeCascade->Branch("fTreeCascVarAmplitudeV0A",&fTreeCascVarAmplitudeV0A,"fTreeCascVarAmplitudeV0A/F");
        fTreeCascade->Branch("fTreeCascVarAmplitudeV0C",&fTreeCascVarAmplitudeV0C,"fTreeCascVarAmplitudeV0C/F");
        fTreeCascade->Branch("fTreeCascVarClosestNonEmptyBC",&fTreeCascVarClosestNonEmptyBC,"fTreeCascVarClosestNonEmptyBC/I");
        
    }
    
    //Kink tagging
    // fTreeCascade->Branch("fTreeCascVarBachIsKink",&fTreeCascVarBachIsKink,"fTreeCascVarBachIsKink/O");
    // fTreeCascade->Branch("fTreeCascVarPosIsKink",&fTreeCascVarPosIsKink,"fTreeCascVarPosIsKink/O");
    // fTreeCascade->Branch("fTreeCascVarNegIsKink",&fTreeCascVarNegIsKink,"fTreeCascVarNegIsKink/O");
    
    
    // fTreeCascade->Branch("fTreeCascVarCowboyness",&fTreeCascVarCowboyness,"fTreeCascVarCowboyness/F");
    //  fTreeCascade->Branch("fTreeCascVarCascadeCowboyness",&fTreeCascVarCascadeCowboyness,"fTreeCascVarCascadeCowboyness/F");
    // }
    
    //}
    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    //-----------------------------------------------------
    
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    if (!inputHandler) {
        AliWarning("--------------------------------------- NO INPUT HANDLER FOUND---------------------------------------- \n");
        return;
    }
    else AliWarning("--------------------------------------- INPUT HANDLER FOUND indeed---------------------------------------- \n");
    
    
    //fPIDResponse = dynamic_cast<AliPIDResponse*>(inputHandler->GetPIDResponse());
    fPIDResponse = inputHandler->GetPIDResponse();
    // fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
    // fPIDResponse = dynamic_cast<AliPIDResponse*>(inputHandler->GetPIDResponse());
    if (!fPIDResponse) {
        AliWarning("-------------------------------------------- NO FPIDRESPONSE FOUND!----------------------------------------- \n");
        return;
    }
    else AliWarning("-------------------------------------------- FPIDRESPONSE FOUND indeed!----------------------------------------- \n");
    
    
    
    inputHandler->SetNeedField();
    
    // Multiplicity
    if(! fESDtrackCuts ) {
        fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,kFALSE);
        fESDtrackCuts->SetPtRange(0.15);  // adding pt cut
        fESDtrackCuts->SetEtaRange(-1.0, 1.0);
    }
    //Analysis Utils
    if(! fUtils ) {
        fUtils = new AliAnalysisUtils();
    }
    if(! fRand ){
        fRand = new TRandom3();
        // From TRandom3 reference:
        // if seed is 0 (default value) a TUUID is generated and
        // used to fill the first 8 integers of the seed array
        fRand->SetSeed(0);
    }
    
    // OOB Pileup in pp 2016
    if( !fESDtrackCutsGlobal2015 && fkDebugOOBPileup ) {
        fESDtrackCutsGlobal2015 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb(kTRUE,kFALSE);
        //Initial set of cuts - to be adjusted
        fESDtrackCutsGlobal2015->SetPtRange(0.15);
        fESDtrackCutsGlobal2015->SetEtaRange(-1.0, 1.0);
    }
    if( !fESDtrackCutsITSsa2010 && fkDebugOOBPileup ) {
        fESDtrackCutsITSsa2010 = AliESDtrackCuts::GetStandardITSSATrackCuts2010();
    }
    
    //------------------------------------------------
    // V0 Multiplicity Histograms
    //------------------------------------------------
    
    // Create histograms
    fListHist = new TList();
    fListHist->SetOwner("kTRUE");  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
    
    //Event Cuts with strong anti-pilep cuts
    fEventCutsStrictAntipileup.fUseStrongVarCorrelationCut = true;
    fEventCutsStrictAntipileup.fUseVariablesCorrelationCuts = true;
    
    //Add QA plots
    if( fkPileupRejectionMode == 0 ){
        fEventCuts.AddQAplotsToList(fListHist,kTRUE);
    }else{
        fEventCutsStrictAntipileup.AddQAplotsToList(fListHist,kTRUE);
    }
    
    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
        fListHist->Add(fHistEventCounter);
    }
    
    if(! fHistEventCounterDifferential ) {
        //Histogram Output: Event-by-Event
        fHistEventCounterDifferential = new TH1D( "fHistEventCounterDifferential", ";Evt. Sel. Step;Count",2,0,2);
        fHistEventCounterDifferential->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounterDifferential->GetXaxis()->SetBinLabel(2, "Selected");
        fListHist->Add(fHistEventCounterDifferential);
    }
    
    if(! fHistCentrality ) {
        //Histogram Output: Event-by-Event
        fHistCentrality = new TH1D( "fHistCentrality", "WARNING: no pileup rejection applied!;Centrality;Event Count",100,0,100);
        fListHist->Add(fHistCentrality);
    }
    
    //Superlight mode output
    
    // cout << "(fListK0Short, fListLambda,fListAntiLambda, fListXiMinus,fListXiPlus,fListOmegaMinus,fListOmegaPlus) = " << fListK0Short << fListLambda << fListAntiLambda
    //     <<fListXiMinus << fListXiPlus <<fListOmegaMinus << fListOmegaPlus << endl;
    //list are created since there is no error about their connection to output slots (postdata); check
    if ( !fListK0Short    ){ fListK0Short    = new TList();    fListK0Short->SetOwner();    }
    if ( !fListLambda     ){ fListLambda     = new TList();    fListLambda->SetOwner();     }
    if ( !fListAntiLambda ){ fListAntiLambda = new TList();    fListAntiLambda->SetOwner(); }
    
    if ( !fListXiMinus ){ fListXiMinus = new TList();    fListXiMinus->SetOwner(); }
    if ( !fListXiPlus  ){ fListXiPlus = new TList();     fListXiPlus->SetOwner();  }
    if ( !fListOmegaMinus ){ fListOmegaMinus = new TList();    fListOmegaMinus->SetOwner(); }
    if ( !fListOmegaPlus  ){ fListOmegaPlus = new TList();     fListOmegaPlus->SetOwner();  }
    
    //Initialize user objects of cascade type
    Int_t lNbrConfigs = 0, lTotalCfgs = 0;
    AliCascadeResult *lCscRslt = 0x0;
    AliV0Result *lV0Rslt = 0x0;
    
    //K0Short configurations
    lNbrConfigs = fListK0Short->GetEntries();
    lTotalCfgs += lNbrConfigs;
    for(Int_t lcfg=0; lcfg<lNbrConfigs; lcfg++){
        lV0Rslt = (AliV0Result*) fListK0Short->At(lcfg);
        lV0Rslt->InitializeHisto();
    }
    //Lambda configurations
    lNbrConfigs = fListLambda->GetEntries();
    lTotalCfgs += lNbrConfigs;
    for(Int_t lcfg=0; lcfg<lNbrConfigs; lcfg++){
        lV0Rslt = (AliV0Result*) fListLambda->At(lcfg);
        lV0Rslt->InitializeHisto();
        lV0Rslt->InitializeProtonProfile();
        lV0Rslt->InitializeFeeddownMatrix();
    }
    //Lambda configurations
    lNbrConfigs = fListAntiLambda->GetEntries();
    lTotalCfgs += lNbrConfigs;
    for(Int_t lcfg=0; lcfg<lNbrConfigs; lcfg++){
        lV0Rslt = (AliV0Result*) fListAntiLambda->At(lcfg);
        lV0Rslt->InitializeHisto();
        lV0Rslt->InitializeProtonProfile();
        lV0Rslt->InitializeFeeddownMatrix();
    }
    //XiMinus configurations
    lNbrConfigs = fListXiMinus->GetEntries();
    lTotalCfgs += lNbrConfigs;
    for(Int_t lcfg=0; lcfg<lNbrConfigs; lcfg++){
        lCscRslt = (AliCascadeResult*) fListXiMinus->At(lcfg);
        lCscRslt->InitializeHisto();
        lCscRslt->InitializeProtonProfile();
    }
    //XiPlus configurations
    lNbrConfigs = fListXiPlus->GetEntries();
    lTotalCfgs += lNbrConfigs;
    for(Int_t lcfg=0; lcfg<lNbrConfigs; lcfg++){
        lCscRslt = (AliCascadeResult*) fListXiPlus->At(lcfg);
        lCscRslt->InitializeHisto();
        lCscRslt->InitializeProtonProfile();
    }
    //OmegaMinus configurations
    lNbrConfigs = fListOmegaMinus->GetEntries();
    lTotalCfgs += lNbrConfigs;
    for(Int_t lcfg=0; lcfg<lNbrConfigs; lcfg++){
        lCscRslt = (AliCascadeResult*) fListOmegaMinus->At(lcfg);
        lCscRslt->InitializeHisto();
        lCscRslt->InitializeProtonProfile();
    }
    //OmegaPlus configurations
    lNbrConfigs = fListOmegaPlus->GetEntries();
    lTotalCfgs += lNbrConfigs;
    for(Int_t lcfg=0; lcfg<lNbrConfigs; lcfg++){
        lCscRslt = (AliCascadeResult*) fListOmegaPlus->At(lcfg);
        lCscRslt->InitializeHisto();
        lCscRslt->InitializeProtonProfile();
    }
    
    AliWarning( Form("Initialized %i cascade output objects!", lTotalCfgs));
    
    
    
    
    
    //--------------------extra light alternative event list------------------------------------------------------
    fEventList = new TList();
    fEventList -> SetOwner("kTRUE");
    
    if(! fHistCentralityFast ) {
        //Histogram Output: Event-by-Event
        fHistCentralityFast = new TH1D( "fHistCentralityFast", "WARNING: no pileup rejection applied!;Centrality;Event Count",100,0,100);
        fEventList->Add(fHistCentralityFast);
    }
    
    if(! fHistZVertexEvent ) {
        //Histogram Output: Event-by-Event
        fHistZVertexEvent = new TH1D( "fHistZVertexEvent", "z vertex",100,-11,11);
        fEventList->Add(fHistZVertexEvent);
    }
    //--------------------end of extra light alternative event list------------------------------------------------------
    
    
    
    
    //-----store the cascade tree as an event tree-----------------------------
    // Cascade_Track = NULL;
    Cascade_Track = new AliRunningCascadeTrack();
    Cascade_Event = new AliRunningCascadeEvent();
    
    fTreeCascadeAsEvent = NULL;
    fTreeCascadeAsEvent = new TTree("fTreeCascadeAsEvent", "cascade tree as event");
    fTreeCascadeAsEvent->Branch("fTreeCascadeAsEvent_branch", "as event branch", Cascade_Event);
    
    
    
    PostData(1, fTreeCascadeAsEvent);
    //PostData(2, fTreeCascade);
    
    
    
    
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::UserExec(Option_t *)
{
    
    // Main loop
    // Called for each event;
    Double_t cosAAA(-2);
    Double_t dzeta(1e8);
    
    AliESDEvent *lESDevent = 0x0;
    
    
    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of V0s in it.
    
    // Appropriate for ESD analysis!
    
    lESDevent = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return;
    }
    
    //Get VZERO Information for multiplicity later
    AliVVZERO* esdV0 = lESDevent->GetVZEROData();
    if (!esdV0) {
        AliError("AliVVZERO not available");
        return;
    }
    
    Double_t lMagneticField = -10;
    lMagneticField = lESDevent->GetMagneticField( );
    //run number
    
    
    //------------------------------------------------
    // Retrieving IR info for OOB Pileup rejection
    //------------------------------------------------
    if( fkDebugOOBPileup ) {
        fClosestNonEmptyBC = 10*3564; // start with an isolated event
        AliESDHeader* lESDHeader = (AliESDHeader*)lESDevent->GetHeader();
        Int_t    nIRs       = lESDHeader->GetTriggerIREntries();
        Long64_t lThisOrbit = lESDHeader->GetOrbitNumber();
        Int_t    lThisBC    = lESDHeader->GetBunchCrossNumber();
        
        const AliTriggerIR* lIR;
        for(Int_t i=0; i<nIRs; i++) {
            
            lIR = lESDHeader->GetTriggerIR(i);
            Long64_t lOrbit     = lIR->GetOrbit();
            UInt_t   lNWord     = lIR->GetNWord();
            UShort_t *lBCsForIR = lIR->GetBCs();
            Bool_t   *lInt1     = lIR->GetInt1s();
            Bool_t   *lInt2     = lIR->GetInt2s();
            
            for(UInt_t j=0; j<lNWord; j++) {
                
                if( (lInt1[j]) || (lInt2[j]) ) {
                    
                    Int_t lBC = lBCsForIR[j];
                    
                    if((lOrbit == lThisOrbit) && (lBC != lThisBC)) {
                        Int_t lClosestNonEmptyBC = lBC - lThisBC;
                        if(TMath::Abs(lClosestNonEmptyBC)<TMath::Abs(fClosestNonEmptyBC)) fClosestNonEmptyBC = lClosestNonEmptyBC;
                    }
                    
                    if(lOrbit == (lThisOrbit+1)) {
                        Int_t lClosestNonEmptyBC = (lBC+3564) - lThisBC;
                        if(TMath::Abs(lClosestNonEmptyBC)<TMath::Abs(fClosestNonEmptyBC)) fClosestNonEmptyBC = lClosestNonEmptyBC;
                    }
                    
                    if(lOrbit == (lThisOrbit-1)) {
                        Int_t lClosestNonEmptyBC = (lBC-3564) - lThisBC;
                        if(TMath::Abs(lClosestNonEmptyBC)<TMath::Abs(fClosestNonEmptyBC)) fClosestNonEmptyBC = lClosestNonEmptyBC;
                    }
                }
            }
        }
    }
    // done with IR ----------------------------------
    
    //------------------------------------------------
    // Event Selection ---
    //  --- Performed entirely via AliPPVsMultUtils
    // (except removal of incomplete events and SPDClusterVsTracklets cut)
    //------------------------------------------------
    
    //Copy-paste of steps done in AliAnalysisTaskSkeleton
    
    fHistEventCounter->Fill(0.5);
    
    //------------------------------------------------
    // Primary Vertex Requirements Section:
    //  ---> pp: has vertex, |z|<10cm
    //------------------------------------------------
    
    //classical Proton-proton like selection
    const AliESDVertex *lPrimaryBestESDVtx     = lESDevent->GetPrimaryVertex();
    // const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    //  const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();
    
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    
    //sandbox this info, please
    
    fTreeVariablePrimVertexX = lBestPrimaryVtxPos[0];
    fTreeVariablePrimVertexY = lBestPrimaryVtxPos[1];
    fTreeVariablePrimVertexZ = lBestPrimaryVtxPos[2];
    
    Bool_t fAcceptZvrt = kFALSE;
    if (TMath::Abs(fTreeVariablePrimVertexZ) <= 10.)fAcceptZvrt = kTRUE;
    else fAcceptZvrt = kFALSE;
    if (fAcceptZvrt == kFALSE) return;
    
    
    //ITS+TPC Esd track multiplicity-----------------------------------------------------------------------------
    AliESDtrackCuts* esdtrackcuts = new AliESDtrackCuts("esdtrackcuts", "esdtrackcuts");
    // enum trackType = esdtrackcuts->MultEstTrackType(kTrackletsITSTPC);
    //enum AliESDtrackCuts::MultEstTrackType trackType;
    fkEsdTrackMultiplicity = esdtrackcuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTracklets, 0.8, 0.);
    if (fkEsdTrackMultiplicity < 0.) return;
    
    //  ULong64_t fkEsdEventTriggerWord1 = lESDevent->GetTriggerMask();
    
    fkEsdEventTriggerWord =lESDevent->GetTriggerMask();
    
    // TString tw = lESDevent->GetFiredTriggerClasses();
    // cout << "Trigger classes are " << tw << endl;
    //  cout << "Trigger mask = " << fkEsdEventTriggerWord1 <<endl;
    //if (fkEsdEventTriggerWord < 0) return;
    //-----------------------------------------------------------------------------------------------------------
    
    
    fHistZVertexEvent -> Fill(fTreeVariablePrimVertexZ);
    
    fTreeVariableMagneticField = lMagneticField;
    
    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------
    
    Float_t lPercentile = 500;
    // Int_t lEvSelCode = 100;
    AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    } else {
        //V0M Multiplicity Percentile
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
        //cout << "lPercentile = " << lPercentile << endl;  //199 for pp 2.46
        //Event Selection Code
        lEvSelCode = MultSelection->GetEvSelCode();
    }
    
    //just ask AliMultSelection. It will know.
    fMVPileupFlag = kFALSE;
    fMVPileupFlag = MultSelection->GetThisEventIsNotPileupMV();
    
    fCentrality = lPercentile;
    
    
    //===================================================================
    //Override centrality with equivalent run 1 info if requested, please
    if (fkUseOldCentrality) {
        AliCentrality* centrality;
        centrality = lESDevent->GetCentrality();
        if ( centrality ) {
            fCentrality = centrality->GetCentralityPercentile( "V0M" );
        }
    }
    //===================================================================
    
    AliVEvent *ev = InputEvent();
    
    fHistEventCounter->Fill(1.5);
    
    //Bookkeep event number for debugging
    fTreeCascVarEventNumber =
    ( ( ((ULong64_t)lESDevent->GetPeriodNumber() ) << 36 ) |
     ( ((ULong64_t)lESDevent->GetOrbitNumber () ) << 12 ) |
     ((ULong64_t)lESDevent->GetBunchCrossNumber() )  );
    
    //Save info for pileup study (high multiplicity triggers pp 13 TeV - 2016 data)
    if( fkDebugOOBPileup ) {
        fOOBPileupFlag     = !fUtils->IsOutOfBunchPileUp(ev);
        fNTOFClusters      = lESDevent->GetESDTOFClusters()->GetEntriesFast();
        fNTOFMatches       = lESDevent->GetESDTOFMatches()->GetEntriesFast();
        fNTracksITSsa2010  = 0;
        fNTracksGlobal2015 = 0;
        fNTracksGlobal2015TriggerPP = 0;
        //Count tracks with various selections
        for(Long_t itrack = 0; itrack<lESDevent->GetNumberOfTracks(); itrack++) {
            AliVTrack *track = lESDevent -> GetVTrack( itrack );
            if( !track ) continue;
            //Only ITSsa tracks
            if( fESDtrackCutsITSsa2010->AcceptVTrack(track) ) fNTracksITSsa2010++;
            if( !fESDtrackCutsGlobal2015->AcceptVTrack(track) ) continue;
            //Only for accepted tracks
            fNTracksGlobal2015++;
            //Count accepted + TOF time window (info from Alberica)
            //Warning: 12.5 is appropriate for pp (for Pb-Pb use 30)
            if( TMath::Abs( track->GetTOFExpTDiff() ) < 12.5 ) fNTracksGlobal2015TriggerPP++;
        }
        
        //VZERO info
        AliMultEstimator *fEstV0A = 0x0, *fEstV0C = 0x0;
        fEstV0A = (AliMultEstimator*)MultSelection->GetEstimator("V0A");
        fEstV0C = (AliMultEstimator*)MultSelection->GetEstimator("V0C");
        if ( fEstV0A ) fAmplitudeV0A = fEstV0A->GetValue();
        if ( fEstV0C ) fAmplitudeV0C = fEstV0C->GetValue();
        
    }
    
    //Fill centrality histogram
    fHistCentrality->Fill(fCentrality);
    fHistCentralityFast->Fill(fCentrality);
    // cout << "Centrality value = " << fCentrality << endl;
    
    //Event-level fill
    //if (fkSaveEventTree)
    
    //fTreeEvent->Fill();
    
    
    fNtracks  = lESDevent ->GetNumberOfTracks();
    fRunNum = lESDevent ->GetRunNumber();
    
    
    //fill the information about the events
    Cascade_Event-> ClearTrackList();
    Cascade_Event-> setx( fTreeVariablePrimVertexX );
    Cascade_Event-> sety( fTreeVariablePrimVertexY );
    Cascade_Event-> setz( fTreeVariablePrimVertexZ );
    Cascade_Event-> setid(fRunNum); //to do!
    Cascade_Event-> setN_tracks(fNtracks); //to do!
    Cascade_Event-> setcentrality(fCentrality);
    Cascade_Event-> setMVPPileUpFlag(fMVPileupFlag);
    Cascade_Event-> setmultiplicity(fkEsdTrackMultiplicity);
    Cascade_Event-> settrigger_word(fkEsdEventTriggerWord);
    Cascade_Event-> setmagfield(lMagneticField);
    
    
    
    
    //  else Printf("could not be written");
    
    //STOP HERE if skipping event selections (no point in doing the rest...)
    
    //------------------------------------------------
    
    //------------------------------------------------
    // Fill V0 Tree as needed
    //------------------------------------------------
    
    //Variable definition
    Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;
    Double_t lChi2V0 = 0;
    Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
    Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
    Double_t lV0CosineOfPointingAngle = 0;
    Double_t lV0Radius = 0, lPt = 0;
    Double_t lRapK0Short = 0, lRapLambda = 0;
    Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
    Double_t lAlphaV0 = 0, lPtArmV0 = 0;
    
    Double_t fMinV0Pt = 0;
    Double_t fMaxV0Pt = 100;
    
    //------------------------------------------------
    // Rerun V0 Vertexer!
    // WARNING: this will only work if the
    // special "use on the fly cascading" flag
    // is disabled!
    //------------------------------------------------
    
    if( fkRunVertexers && !fkUseOnTheFlyV0Cascading ) {
        //Only reset if not using on-the-fly (or else nothing passes)
        lESDevent->ResetV0s();
        
        //Decide between regular and light vertexer (default: light)
        if ( ! fkUseLightVertexer ){
            //Instantiate vertexer object
            AliV0vertexer lV0vtxer;
            //Set Cuts
            lV0vtxer.SetDefaultCuts(fV0VertexerSels);
            lV0vtxer.SetCuts(fV0VertexerSels);
            //Redo
            lV0vtxer.Tracks2V0vertices(lESDevent);
        } else {
            //Instantiate vertexer object
            AliLightV0vertexer lV0vtxer;
            //Set do or don't do V0 refit for improved precision
            lV0vtxer.SetDoRefit( kFALSE );
            if (fkDoV0Refit) lV0vtxer.SetDoRefit(kTRUE);
            //Set Cuts
            lV0vtxer.SetDefaultCuts(fV0VertexerSels);
            lV0vtxer.SetCuts(fV0VertexerSels);
            //Redo
            lV0vtxer.Tracks2V0vertices(lESDevent);
        }
    }
    
    Int_t nv0s = 0;
    nv0s = lESDevent->GetNumberOfV0s();
    
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
    {   // This is the begining of the V0 loop
        AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
        if (!v0) continue;
        
        CheckChargeV0( v0 );
        //Remove like-sign (will not affect offline V0 candidates!)
        if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() > 0 ){
            continue;
        }
        if( v0->GetParamN()->Charge() < 0 && v0->GetParamP()->Charge() < 0 ){
            continue;
        }
        
        Double_t tDecayVertexV0[3];
        v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);
        
        Double_t tV0mom[3];
        v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
        Double_t lV0TotalMomentum = TMath::Sqrt(
                                                tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );
        
        lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);
        
        lPt = v0->Pt();
        lRapK0Short = v0->RapK0Short();
        lRapLambda  = v0->RapLambda();
        if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;
        
        UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());
        
        Double_t lMomPos[3];
        v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
        Double_t lMomNeg[3];
        v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
        
        //Provisions for cowboy/sailor check
        Double_t lModp1 = TMath::Sqrt( lMomPos[0]*lMomPos[0] + lMomPos[1]*lMomPos[1] );
        Double_t lModp2 = TMath::Sqrt( lMomNeg[0]*lMomNeg[0] + lMomNeg[1]*lMomNeg[1] );
        
        //Calculate vec prod with momenta projected to xy plane
        Double_t lVecProd = (lMomPos[0]*lMomNeg[1] - lMomPos[1]*lMomNeg[0]) / (lModp1*lModp2);
        
        if ( lMagneticField < 0 ) lVecProd *= -1; //invert sign
        
        fTreeVariableIsCowboy = kFALSE;
        if (lVecProd < 0) fTreeVariableIsCowboy = kTRUE;
        
        AliESDtrack *pTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyPos);
        AliESDtrack *nTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyNeg);
        
        fTreeVariablePosTrack = pTrack;
        fTreeVariableNegTrack = nTrack;
        
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
        fTreeVariablePosPIDForTracking = pTrack->GetPIDForTracking();
        fTreeVariableNegPIDForTracking = nTrack->GetPIDForTracking();
        
        //Check its clusters
        fTreeVariablePosITSClusters0 = pTrack->HasPointOnITSLayer(0);
        fTreeVariablePosITSClusters1 = pTrack->HasPointOnITSLayer(1);
        fTreeVariablePosITSClusters2 = pTrack->HasPointOnITSLayer(2);
        fTreeVariablePosITSClusters3 = pTrack->HasPointOnITSLayer(3);
        fTreeVariablePosITSClusters4 = pTrack->HasPointOnITSLayer(4);
        fTreeVariablePosITSClusters5 = pTrack->HasPointOnITSLayer(5);
        
        fTreeVariableNegITSClusters0 = nTrack->HasPointOnITSLayer(0);
        fTreeVariableNegITSClusters1 = nTrack->HasPointOnITSLayer(1);
        fTreeVariableNegITSClusters2 = nTrack->HasPointOnITSLayer(2);
        fTreeVariableNegITSClusters3 = nTrack->HasPointOnITSLayer(3);
        fTreeVariableNegITSClusters4 = nTrack->HasPointOnITSLayer(4);
        fTreeVariableNegITSClusters5 = nTrack->HasPointOnITSLayer(5);
        
        //Check its clusters, shared
        fTreeVariablePosITSSharedClusters0 = pTrack->HasSharedPointOnITSLayer(0);
        fTreeVariablePosITSSharedClusters1 = pTrack->HasSharedPointOnITSLayer(1);
        fTreeVariablePosITSSharedClusters2 = pTrack->HasSharedPointOnITSLayer(2);
        fTreeVariablePosITSSharedClusters3 = pTrack->HasSharedPointOnITSLayer(3);
        fTreeVariablePosITSSharedClusters4 = pTrack->HasSharedPointOnITSLayer(4);
        fTreeVariablePosITSSharedClusters5 = pTrack->HasSharedPointOnITSLayer(5);
        
        fTreeVariableNegITSSharedClusters0 = nTrack->HasSharedPointOnITSLayer(0);
        fTreeVariableNegITSSharedClusters1 = nTrack->HasSharedPointOnITSLayer(1);
        fTreeVariableNegITSSharedClusters2 = nTrack->HasSharedPointOnITSLayer(2);
        fTreeVariableNegITSSharedClusters3 = nTrack->HasSharedPointOnITSLayer(3);
        fTreeVariableNegITSSharedClusters4 = nTrack->HasSharedPointOnITSLayer(4);
        fTreeVariableNegITSSharedClusters5 = nTrack->HasSharedPointOnITSLayer(5);
        
        const AliExternalTrackParam *innernegv0=nTrack->GetInnerParam();
        const AliExternalTrackParam *innerposv0=pTrack->GetInnerParam();
        Float_t lThisPosInnerP = -1;
        Float_t lThisNegInnerP = -1;
        Float_t lThisPosInnerPt = -1;
        Float_t lThisNegInnerPt = -1;
        if(innerposv0)  { lThisPosInnerP  = innerposv0 ->GetP(); }
        if(innernegv0)  { lThisNegInnerP  = innernegv0 ->GetP(); }
        if(innerposv0)  { lThisPosInnerPt  = innerposv0 ->Pt(); }
        if(innernegv0)  { lThisNegInnerPt  = innernegv0 ->Pt(); }
        Float_t lThisPosdEdx = pTrack -> GetTPCsignal();
        Float_t lThisNegdEdx = nTrack -> GetTPCsignal();
        
        fTreeVariablePosdEdx = lThisPosdEdx;
        fTreeVariableNegdEdx = lThisNegdEdx;
        
        fTreeVariablePosInnerP = lThisPosInnerP;
        fTreeVariableNegInnerP = lThisNegInnerP;
        
        //Daughter Eta for Eta selection, afterwards
        fTreeVariableNegEta = nTrack->Eta();
        fTreeVariablePosEta = pTrack->Eta();
        
        if ( fkExtraCleanup ){
            if( TMath::Abs(fTreeVariableNegEta)>0.8 || TMath::Abs(fTreeVariableNegEta)>0.8 ) continue;
            if( TMath::Abs(lRapK0Short        )>0.5 && TMath::Abs(lRapLambda         )>0.5 ) continue;
        }
        
        // Filter like-sign V0 (next: add counter and distribution)
        if ( pTrack->GetSign() == nTrack->GetSign()) {
            continue;
        }
        
        //________________________________________________________________________
        // Track quality cuts
        Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
        Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
        fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
        if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
            fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
        
        // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
        if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        
        //Get status flags
        fTreeVariablePosTrackStatus = pTrack->GetStatus();
        fTreeVariableNegTrackStatus = nTrack->GetStatus();
        
        fTreeVariablePosDCAz = GetDCAz(pTrack);
        fTreeVariableNegDCAz = GetDCAz(nTrack);
        
        //GetKinkIndex condition
        if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;
        
        //Findable clusters > 0 condition
        if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;
        
        //Compute ratio Crossed Rows / Findable clusters
        //Note: above test avoids division by zero!
        Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF()));
        Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF()));
        
        fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
        if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable )
            fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
        
        //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
        //if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;
        
        //Extra track quality: Chi2/cluster for cross-checks
        Float_t lBiggestChi2PerCluster = -1;
        
        Float_t lPosChi2PerCluster = 1000;
        Float_t lNegChi2PerCluster = 1000;
        
        if( pTrack->GetTPCNcls() > 0 ) lPosChi2PerCluster = pTrack->GetTPCchi2() / ((Float_t)pTrack->GetTPCNcls());
        if( nTrack->GetTPCNcls() > 0 ) lNegChi2PerCluster = nTrack->GetTPCchi2() / ((Float_t)nTrack->GetTPCNcls());
        
        if ( lPosChi2PerCluster  > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lPosChi2PerCluster;
        if ( lNegChi2PerCluster  > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lNegChi2PerCluster;
        
        fTreeVariableMaxChi2PerCluster = lBiggestChi2PerCluster;
        
        //Extra track quality: min track length
        Float_t lSmallestTrackLength = 1000;
        Float_t lPosTrackLength = -1;
        Float_t lNegTrackLength = -1;
        
        if (pTrack->GetInnerParam()) lPosTrackLength = pTrack->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if (nTrack->GetInnerParam()) lNegTrackLength = nTrack->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        
        if ( lPosTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lPosTrackLength;
        if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
        
        fTreeVariableMinTrackLength = lSmallestTrackLength;
        
        if ( ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) && lSmallestTrackLength<80 && fkExtraCleanup ) continue;
        
        //End track Quality Cuts
        //________________________________________________________________________
        
        //returns only the transverse impact parameter: standard
        lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lBestPrimaryVtxPos[0],
                                                      lBestPrimaryVtxPos[1],
                                                      lMagneticField) );
        
        
        
        lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lBestPrimaryVtxPos[0],
                                                      lBestPrimaryVtxPos[1],
                                                      lMagneticField));
        
        
        
        
        lOnFlyStatus = v0->GetOnFlyStatus();
        lChi2V0 = v0->GetChi2V0();
        lDcaV0Daughters = v0->GetDcaV0Daughters();
        lDcaV0ToPrimVertex = v0->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        fTreeVariableV0CosineOfPointingAngle=lV0CosineOfPointingAngle;
        
        // Getting invariant mass infos directly from ESD
        v0->ChangeMassHypothesis(310);
        lInvMassK0s = v0->GetEffMass();
        v0->ChangeMassHypothesis(3122);
        lInvMassLambda = v0->GetEffMass();
        v0->ChangeMassHypothesis(-3122);
        lInvMassAntiLambda = v0->GetEffMass();
        lAlphaV0 = v0->AlphaV0();
        lPtArmV0 = v0->PtArmV0();
        
        fTreeVariableMVPileupFlag = fMVPileupFlag;
        
        fTreeVariablePt = v0->Pt();
        fTreeVariableChi2V0 = lChi2V0;
        fTreeVariableDcaV0ToPrimVertex = lDcaV0ToPrimVertex;
        fTreeVariableDcaV0Daughters = lDcaV0Daughters;
        fTreeVariableV0CosineOfPointingAngle = lV0CosineOfPointingAngle;
        fTreeVariableV0Radius = lV0Radius;
        fTreeVariableDcaPosToPrimVertex = lDcaPosToPrimVertex;
        fTreeVariableDcaNegToPrimVertex = lDcaNegToPrimVertex;
        fTreeVariableInvMassK0s = lInvMassK0s;
        fTreeVariableInvMassLambda = lInvMassLambda;
        fTreeVariableInvMassAntiLambda = lInvMassAntiLambda;
        fTreeVariableRapK0Short = lRapK0Short;
        fTreeVariableRapLambda = lRapLambda;
        fTreeVariableAlphaV0 = lAlphaV0;
        fTreeVariablePtArmV0 = lPtArmV0;
        
        //Official means of acquiring N-sigmas
        fTreeVariableNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton);
        fTreeVariableNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion);
        fTreeVariableNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton);
        fTreeVariableNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion);
        
        //This requires an Invariant Mass Hypothesis afterwards
        fTreeVariableDistOverTotMom = TMath::Sqrt(
                                                  TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
                                                  TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
                                                  TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
                                                  );
        fTreeVariableDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure
        
        //Copy Multiplicity information
        fTreeVariableCentrality = fCentrality;
        
        //Info for pileup studies
        fTreeVariableNegTOFExpTDiff = nTrack->GetTOFExpTDiff( lESDevent->GetMagneticField() );
        fTreeVariablePosTOFExpTDiff = pTrack->GetTOFExpTDiff( lESDevent->GetMagneticField() );
        fTreeVariableNegTOFSignal = nTrack->GetTOFsignal() * 1.e-3; // in ns
        fTreeVariablePosTOFSignal = pTrack->GetTOFsignal() * 1.e-3; // in ns
        fTreeVariableNegTOFBCid = nTrack->GetTOFBunchCrossing( lESDevent->GetMagneticField() );
        fTreeVariablePosTOFBCid = pTrack->GetTOFBunchCrossing( lESDevent->GetMagneticField() );
        //Copy OOB pileup flag for this event
        fTreeVariableOOBPileupFlag = fOOBPileupFlag;
        //Copy VZERO information for this event
        fTreeVariableAmplitudeV0A = fAmplitudeV0A;
        fTreeVariableAmplitudeV0C = fAmplitudeV0C;
        //Copy IR information for this event
        fTreeVariableClosestNonEmptyBC = fClosestNonEmptyBC;
        
        
        //------------------------------------------------
        // Fill Tree!
        //------------------------------------------------
        
        // The conditionals are meant to decrease excessive
        // memory usage!
        
        //First Selection: Reject OnFly
        /*   if( lOnFlyStatus == 0 ) {
         //Second Selection: rough 20-sigma band, parametric.
         //K0Short: Enough to parametrize peak broadening with linear function.
         Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeVariablePt;
         Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeVariablePt;
         //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
         //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
         Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeVariablePt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeVariablePt);
         Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeVariablePt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeVariablePt);
         //Do Selection
         if(
         //Case 1: Lambda Selection
         (fTreeVariableInvMassLambda <lUpperLimitLambda  && fTreeVariableInvMassLambda  > lLowerLimitLambda &&
         (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasPosProton) < 7.0 && TMath::Abs(fTreeVariableNSigmasNegPion) < 7.0) )
         )
         ||
         //Case 2: AntiLambda Selection
         (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda &&
         (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasNegProton) < 7.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 7.0) )
         )
         ||
         //Case 3: K0Short Selection
         (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short &&
         (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasNegPion)   < 7.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 7.0) )
         ) ) {
         //Pre-selection in case this is AA...
         
         //Random denial
         
         Bool_t lKeepV0 = kTRUE;
         // The random denial was left away because it is not known what its purpose is. 22.01.19
         //  if(fkDownScaleV0 && ( fRand->Uniform() > fDownScaleFactorV0 )) lKeepV0 = kFALSE;
         
         //pT window
         if( fTreeVariablePt < fMinPtToSave ) lKeepV0 = kFALSE;
         if( fTreeVariablePt > fMaxPtToSave ) lKeepV0 = kFALSE;
         
         //  cout << "The vector (NegEta, PosEta, SaveVOtree, keepVO) = (" << TMath::Abs(fTreeVariableNegEta)<<" , "<<
         //    TMath::Abs(fTreeVariablePosEta) << " , " << fkSaveV0Tree << " , " << lKeepV0 << ")" << endl;
         // lKeepV0 is always kFALSE :(
         
         if ( TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8 && fkSaveV0Tree && lKeepV0 )
         {
         fTreeV0->Fill();
         
         }
         }
         } */
        
        //------------------------------------------------
        // Fill V0 tree over.
        //------------------------------------------------
        
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        // Superlight adaptive output mode
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
        //AliWarning(Form("[V0 Analyses] Processing different configurations (%i detected)",lNumberOfConfigurations));
        TH3F *histoout         = 0x0;
        AliV0Result *lV0Result = 0x0;
        
        //pointers to valid results
        AliV0Result *lPointers[50000];
        Long_t lValidConfigurations=0;
        
        for( Int_t icfg=0; icfg<fListK0Short->GetEntries(); icfg++ ){
            lPointers[lValidConfigurations] = (AliV0Result*) fListK0Short->At(icfg);
            lValidConfigurations++;
        }
        for( Int_t icfg=0; icfg<fListLambda->GetEntries(); icfg++ ){
            lPointers[lValidConfigurations] = (AliV0Result*) fListLambda->At(icfg);
            lValidConfigurations++;
        }
        for( Int_t icfg=0; icfg<fListAntiLambda->GetEntries(); icfg++ ){
            lPointers[lValidConfigurations] = (AliV0Result*) fListAntiLambda->At(icfg);
            lValidConfigurations++;
        }
        
        for(Int_t lcfg=0; lcfg<lValidConfigurations; lcfg++){
            lV0Result = lPointers[lcfg];
            histoout  = lV0Result->GetHistogram();
            
            Float_t lMass = 0;
            Float_t lRap  = 0;
            Float_t lPDGMass = -1;
            Float_t lNegdEdx = 100;
            Float_t lPosdEdx = 100;
            Float_t lBaryonMomentum = -0.5;
            Float_t lBaryonPt = -0.5;
            Float_t lBaryondEdxFromProton = 0;
            
            //========================================================================
            //Setting up: Variable V0 CosPA
            Float_t lV0CosPACut = lV0Result -> GetCutV0CosPA();
            Float_t lVarV0CosPApar[5];
            lVarV0CosPApar[0] = lV0Result->GetCutVarV0CosPAExp0Const();
            lVarV0CosPApar[1] = lV0Result->GetCutVarV0CosPAExp0Slope();
            lVarV0CosPApar[2] = lV0Result->GetCutVarV0CosPAExp1Const();
            lVarV0CosPApar[3] = lV0Result->GetCutVarV0CosPAExp1Slope();
            lVarV0CosPApar[4] = lV0Result->GetCutVarV0CosPAConst();
            Float_t lVarV0CosPA = TMath::Cos(
                                             lVarV0CosPApar[0]*TMath::Exp(lVarV0CosPApar[1]*fTreeVariablePt) +
                                             lVarV0CosPApar[2]*TMath::Exp(lVarV0CosPApar[3]*fTreeVariablePt) +
                                             lVarV0CosPApar[4]);
            if( lV0Result->GetCutUseVarV0CosPA() ){
                //Only use if tighter than the non-variable cut
                if( lVarV0CosPA > lV0CosPACut ) lV0CosPACut = lVarV0CosPA;
            }
            //========================================================================
            
            if ( lV0Result->GetMassHypothesis() == AliV0Result::kK0Short     ){
                lMass    = fTreeVariableInvMassK0s;
                lRap     = fTreeVariableRapK0Short;
                lPDGMass = 0.497;
                lNegdEdx = fTreeVariableNSigmasNegPion;
                lPosdEdx = fTreeVariableNSigmasPosPion;
            }
            if ( lV0Result->GetMassHypothesis() == AliV0Result::kLambda      ){
                lMass = fTreeVariableInvMassLambda;
                lRap = fTreeVariableRapLambda;
                lPDGMass = 1.115683;
                lNegdEdx = fTreeVariableNSigmasNegPion;
                lPosdEdx = fTreeVariableNSigmasPosProton;
                lBaryonMomentum = fTreeVariablePosInnerP;
                lBaryonPt = lThisPosInnerPt;
                lBaryondEdxFromProton = fTreeVariableNSigmasPosProton;
            }
            if ( lV0Result->GetMassHypothesis() == AliV0Result::kAntiLambda  ){
                lMass = fTreeVariableInvMassAntiLambda;
                lRap = fTreeVariableRapLambda;
                lPDGMass = 1.115683;
                lNegdEdx = fTreeVariableNSigmasNegProton;
                lPosdEdx = fTreeVariableNSigmasPosPion;
                lBaryonMomentum = fTreeVariableNegInnerP;
                lBaryonPt = lThisNegInnerPt;
                lBaryondEdxFromProton = fTreeVariableNSigmasNegProton;
            }
            
            if (
                //Check 1: Offline Vertexer
                lOnFlyStatus == lV0Result->GetUseOnTheFly() &&
                
                //Check 2: Basic Acceptance cuts
                lV0Result->GetCutMinEtaTracks() < fTreeVariableNegEta && fTreeVariableNegEta < lV0Result->GetCutMaxEtaTracks() &&
                lV0Result->GetCutMinEtaTracks() < fTreeVariablePosEta && fTreeVariablePosEta < lV0Result->GetCutMaxEtaTracks() &&
                lRap > lV0Result->GetCutMinRapidity() &&
                lRap < lV0Result->GetCutMaxRapidity() &&
                
                //Check 3: Topological Variables
                fTreeVariableV0Radius > lV0Result->GetCutV0Radius() &&
                fTreeVariableV0Radius < lV0Result->GetCutMaxV0Radius() &&
                fTreeVariableDcaNegToPrimVertex > lV0Result->GetCutDCANegToPV() &&
                fTreeVariableDcaPosToPrimVertex > lV0Result->GetCutDCAPosToPV() &&
                fTreeVariableDcaV0Daughters < lV0Result->GetCutDCAV0Daughters() &&
                fTreeVariableV0CosineOfPointingAngle > lV0CosPACut &&
                fTreeVariableDistOverTotMom*lPDGMass < lV0Result->GetCutProperLifetime() &&
                fTreeVariableLeastNbrCrossedRows > lV0Result->GetCutLeastNumberOfCrossedRows() &&
                fTreeVariableLeastRatioCrossedRowsOverFindable > lV0Result->GetCutLeastNumberOfCrossedRowsOverFindable() &&
                
                //Check 4: Minimum momentum of baryon daughter
                ( lV0Result->GetMassHypothesis() == AliV0Result::kK0Short || lBaryonMomentum > lV0Result->GetCutMinBaryonMomentum() ) &&
                
                //Check 5: TPC dEdx selections
                TMath::Abs(lNegdEdx)<lV0Result->GetCutTPCdEdx() &&
                TMath::Abs(lPosdEdx)<lV0Result->GetCutTPCdEdx() &&
                
                //Check 6: Armenteros-Podolanski space cut (for K0Short analysis)
                ( ( lV0Result->GetCutArmenteros() == kFALSE || lV0Result->GetMassHypothesis() != AliV0Result::kK0Short ) || ( fTreeVariablePtArmV0>lV0Result->GetCutArmenterosParameter()*TMath::Abs(fTreeVariableAlphaV0) ) ) &&
                
                //Check 7: kITSrefit track selection if requested
                (
                 ( (fTreeVariableNegTrackStatus & AliESDtrack::kITSrefit) &&
                  (fTreeVariablePosTrackStatus & AliESDtrack::kITSrefit) )
                 ||
                 !lV0Result->GetCutUseITSRefitTracks()
                 )&&
                
                //Check 8: Max Chi2/Clusters if not absurd
                ( lV0Result->GetCutMaxChi2PerCluster()>1e+3 ||
                 fTreeVariableMaxChi2PerCluster < lV0Result->GetCutMaxChi2PerCluster()
                 ) &&
                //Check 9: Min Track Length if positive
                ( lV0Result->GetCutMinTrackLength()<0 || //this is a bit paranoid...
                 fTreeVariableMinTrackLength > lV0Result->GetCutMinTrackLength()
                 )&&
                
                //Check 10: Special 2.76TeV-like dedx
                // Logic: either not requested, or K0Short, or high-pT baryon daughter, or passes cut!
                ( !lV0Result->GetCut276TeVLikedEdx() ||
                 ( lV0Result->GetMassHypothesis() == AliV0Result::kK0Short ||
                  ( lBaryonPt > 1.0 || TMath::Abs(lBaryondEdxFromProton)<3.0 )
                  )
                 )&&
                
                //Check 14: has at least one track with some TOF info, please (reject pileup)
                //          warning: this is still to be studied in more detail!
                (
                 lV0Result->GetCutAtLeastOneTOF() == kFALSE ||
                 (
                  TMath::Abs(fTreeVariableNegTOFSignal) < 100 ||
                  TMath::Abs(fTreeVariablePosTOFSignal) < 100
                  )
                 )//&&
                
                //Check 15: cowboy/sailor for V0
                /*  (
                 lV0Result->GetCutIsCowboy()==0 ||
                 (lV0Result->GetCutIsCowboy()== 1 && fTreeVariableIsCowboy==kTRUE ) ||
                 (lV0Result->GetCutIsCowboy()==-1 && fTreeVariableIsCowboy==kFALSE)
                 )//end cowboy/sailor */
                )//end major if
            {
                //This satisfies all my conditionals! Fill histogram
                histoout -> Fill ( fCentrality, fTreeVariablePt, lMass );
            }
        }
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        // End Superlight adaptive output mode
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
    }// This is the end of the V0 loop
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //------------------------------------------------
    // Rerun cascade vertexer!
    //------------------------------------------------
    
    // fkRunVertexers = kTRUE;
    
    
    
    if( fkRunVertexers ) {
        //Remove existing cascades
        lESDevent->ResetCascades();
        
        //Decide between regular and light vertexer (default: light)
        if ( ! fkUseLightVertexer ){
            //Instantiate vertexer object
            AliCascadeVertexer lCascVtxer;
            lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);
            lCascVtxer.SetCuts(fCascadeVertexerSels);
            lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
        } else {
            AliLightCascadeVertexer lCascVtxer;
            lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);
            lCascVtxer.SetCuts(fCascadeVertexerSels);
            if( fkUseOnTheFlyV0Cascading ) lCascVtxer.SetUseOnTheFlyV0(kTRUE);
            lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
        }
    }
    
    //------------------------------------------------
    // MAIN CASCADE LOOP STARTS HERE
    //------------------------------------------------
    // Code Credit: Antonin Maire (thanks^100)
    // ---> This is an adaptation
    
    Long_t ncascades = 0;
    ncascades = lESDevent->GetNumberOfCascades();
    for (Int_t iXi = 0; iXi < ncascades; iXi++)
    {
        
        //------------------------------------------------
        // Initializations
        //------------------------------------------------
        //Double_t lTrkgPrimaryVtxRadius3D = -500.0;
        //Double_t lBestPrimaryVtxRadius3D = -500.0;
        
        // - 1st part of initialisation : variables needed to store AliESDCascade data members
        //Double_t lEffMassXi      = 0. ;
        //Double_t lChi2Xi         = -1. ;
        Double_t lDcaXiDaughters = -1. ;
        Double_t lXiCosineOfPointingAngle = -1. ;
        Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
        Double_t lXiRadius = -1000. ;
        
        // - 2nd part of initialisation : Nbr of clusters within TPC for the 3 daughter cascade tracks
        Int_t    lPosTPCClusters    = -1; // For ESD only ...//FIXME : wait for availability in AOD
        Int_t    lNegTPCClusters    = -1; // For ESD only ...
        Int_t    lBachTPCClusters   = -1; // For ESD only ...
        
        // - 3rd part of initialisation : about V0 part in cascades
        Double_t lInvMassLambdaAsCascDghter = 0.;
        //Double_t lV0Chi2Xi         = -1. ;
        Double_t lDcaV0DaughtersXi = -1.;
        
        Double_t lDcaBachToPrimVertexXi = -1.;
        Double_t lDcaBachToPrimVertexXiZ = -1.;
        Double_t lDcaV0ToPrimVertexXi = -1.;
        Double_t lDcaPosToPrimVertexXi  = -1.;
        Double_t lDcaPosToPrimVertexXiZ  = -1.;
        Double_t lDcaNegToPrimVertexXi  = -1.;
        Double_t lDcaNegToPrimVertexXiZ  = -1.;
        Double_t lV0CosineOfPointingAngleXi = -1. ;
        Double_t lV0CosineOfPointingAngleXiSpecial = -1. ;
        Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
        Double_t lV0RadiusXi = -1000.0;
        Double_t lV0quality  = 0.;
        
        // - 4th part of initialisation : Effective masses
        Double_t lInvMassXiMinus    = 0.;
        Double_t lInvMassXiPlus     = 0.;
        Double_t lInvMassOmegaMinus = 0.;
        Double_t lInvMassOmegaPlus  = 0.;
        //for the estimation of the background shape for the inv mass distribution
        Double_t lInvMassOmegaPlusBackground = 0.;
        Double_t lInvMassOmegaMinusBackground = 0.;
        Double_t lInvMassXiPlusBackground = 0.;
        Double_t lInvMassXiMinusBackground = 0.;
        
        fTreeCascVarChiSquareV0      = 1e+3;
        fTreeCascVarChiSquareCascade = 1e+3;
        
        // - 6th part of initialisation : extra info for QA
        Double_t lXiMomX       = 0. , lXiMomY = 0., lXiMomZ = 0.;
        //  Double_t lXiMomVec[3] = {0.,0.,0.};
        Double_t lXiTransvMom  = 0. ;
        //Double_t lXiTransvMomMC= 0. ;
        Double_t lXiTotMom     = 0. ;
        
        Double_t lBachMomX       = 0., lBachMomY  = 0., lBachMomZ   = 0.;
        //Double_t lBachTransvMom  = 0.;
        //Double_t lBachTotMom     = 0.;
        
        fTreeCascVarNegNSigmaPion   = -100;
        fTreeCascVarNegNSigmaProton = -100;
        fTreeCascVarPosNSigmaPion   = -100;
        fTreeCascVarPosNSigmaProton = -100;
        fTreeCascVarBachNSigmaPion  = -100;
        fTreeCascVarBachNSigmaKaon  = -100;
        
        fTreeCascVarNegTOFNSigmaPion   = -100;
        fTreeCascVarNegTOFNSigmaProton = -100;
        fTreeCascVarPosTOFNSigmaPion   = -100;
        fTreeCascVarPosTOFNSigmaProton = -100;
        fTreeCascVarBachTOFNSigmaPion  = -100;
        fTreeCascVarBachTOFNSigmaKaon  = -100;
        
        fTreeCascVarBachIsKink = kFALSE;
        fTreeCascVarPosIsKink = kFALSE;
        fTreeCascVarNegIsKink = kFALSE;
        
        
        Short_t  lChargeXi = -2;
        //Double_t lV0toXiCosineOfPointingAngle = 0. ;
        
        Double_t lRapXi   = -20.0, lRapOmega = -20.0; //  lEta = -20.0, lTheta = 360., lPhi = 720. ;
        Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
        
        // -------------------------------------
        // II.ESD - Calculation Part dedicated to Xi vertices (ESD)
        
        AliESDcascade *xi = lESDevent->GetCascade(iXi);
        if (!xi) continue;
        
        // - II.Step 2 : Assigning the necessary variables for specific AliESDcascade data members (ESD)
        //-------------
        lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay
        
        lEffMassXi              = xi->GetEffMassXi();
        
        //ChiSquare implementation
        fTreeCascVarChiSquareV0      = xi->GetChi2V0();
        fTreeCascVarChiSquareCascade = xi->GetChi2Xi();
        
        lDcaXiDaughters     = xi->GetDcaXiDaughters();
        lXiCosineOfPointingAngle = xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                       lBestPrimaryVtxPos[1],
                                                                       lBestPrimaryVtxPos[2] );
        // Take care : the best available vertex should be used (like in AliCascadeVertexer)
        
        xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] );
        lXiRadius = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
        
        fTreeCascVarCascadeDecayX = lPosXi[0];
        fTreeCascVarCascadeDecayY = lPosXi[1];
        fTreeCascVarCascadeDecayZ = lPosXi[2];
        
        // - II.Step 3 : around the tracks : Bach + V0 (ESD)
        // ~ Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
        //-------------
        
        UInt_t lIdxPosXi     = (UInt_t) TMath::Abs( xi->GetPindex() );
        UInt_t lIdxNegXi     = (UInt_t) TMath::Abs( xi->GetNindex() );
        UInt_t lBachIdx     = (UInt_t) TMath::Abs( xi->GetBindex() );
        
        // Care track label can be negative in MC production (linked with the track quality)
        // However = normally, not the case for track index ...
        
        // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
        if(lBachIdx == lIdxNegXi) {
            AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!");
            continue;
        }
        if(lBachIdx == lIdxPosXi) {
            AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!");
            continue;
        }
        
        AliESDtrack *pTrackXi        = lESDevent->GetTrack( lIdxPosXi );
        AliESDtrack *nTrackXi        = lESDevent->GetTrack( lIdxNegXi );
        AliESDtrack *bachTrackXi    = lESDevent->GetTrack( lBachIdx );
        
        
        
        Bool_t goodpTrackXi = kFALSE;
        Bool_t goodnTrackXi = kFALSE;
        Bool_t goodbachTrackXi = kFALSE;
        goodpTrackXi = GoodTrack(pTrackXi);
        goodnTrackXi = GoodTrack(nTrackXi);
        goodbachTrackXi = GoodTrack(bachTrackXi);
        
        
        //Sandbox information: always, regardless of status
        fTreeCascVarBachTrack = bachTrackXi;
        fTreeCascVarPosTrack = pTrackXi;
        fTreeCascVarNegTrack = nTrackXi;
        
        fTreeCascVarNegIndex  = lIdxNegXi;
        fTreeCascVarPosIndex  = lIdxPosXi;
        fTreeCascVarBachIndex = lBachIdx;
        
        if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
            AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
            continue;
        }
        
        fTreeCascVarPosEta = pTrackXi->Eta();
        fTreeCascVarNegEta = nTrackXi->Eta();
        fTreeCascVarBachEta = bachTrackXi->Eta();
        
        
        //momenta
        fTreeCascVarPosPt = pTrackXi -> Pt();
        fTreeCascVarPosPz = pTrackXi->Pz();
        fTreeCascVarPosP = pTrackXi->P();
        
        fTreeCascVarNegPt = nTrackXi -> Pt();
        fTreeCascVarNegPz = nTrackXi->Pz();
        fTreeCascVarNegP = nTrackXi->P();
        
        fTreeCascVarBachPt = pTrackXi -> Pt();
        fTreeCascVarBachPz = pTrackXi->Pz();
        fTreeCascVarBachP = pTrackXi->P();
        
        //Check its clusters
        fTreeCascVarPosITSClusters0 = pTrackXi->HasPointOnITSLayer(0);
        fTreeCascVarPosITSClusters1 = pTrackXi->HasPointOnITSLayer(1);
        fTreeCascVarPosITSClusters2 = pTrackXi->HasPointOnITSLayer(2);
        fTreeCascVarPosITSClusters3 = pTrackXi->HasPointOnITSLayer(3);
        fTreeCascVarPosITSClusters4 = pTrackXi->HasPointOnITSLayer(4);
        fTreeCascVarPosITSClusters5 = pTrackXi->HasPointOnITSLayer(5);
        
        fTreeCascVarNegITSClusters0 = nTrackXi->HasPointOnITSLayer(0);
        fTreeCascVarNegITSClusters1 = nTrackXi->HasPointOnITSLayer(1);
        fTreeCascVarNegITSClusters2 = nTrackXi->HasPointOnITSLayer(2);
        fTreeCascVarNegITSClusters3 = nTrackXi->HasPointOnITSLayer(3);
        fTreeCascVarNegITSClusters4 = nTrackXi->HasPointOnITSLayer(4);
        fTreeCascVarNegITSClusters5 = nTrackXi->HasPointOnITSLayer(5);
        
        fTreeCascVarBachITSClusters0 = bachTrackXi->HasPointOnITSLayer(0);
        fTreeCascVarBachITSClusters1 = bachTrackXi->HasPointOnITSLayer(1);
        fTreeCascVarBachITSClusters2 = bachTrackXi->HasPointOnITSLayer(2);
        fTreeCascVarBachITSClusters3 = bachTrackXi->HasPointOnITSLayer(3);
        fTreeCascVarBachITSClusters4 = bachTrackXi->HasPointOnITSLayer(4);
        fTreeCascVarBachITSClusters5 = bachTrackXi->HasPointOnITSLayer(5);
        
        //Check its clusters, shared
        fTreeCascVarPosITSSharedClusters0 = pTrackXi->HasSharedPointOnITSLayer(0);
        fTreeCascVarPosITSSharedClusters1 = pTrackXi->HasSharedPointOnITSLayer(1);
        fTreeCascVarPosITSSharedClusters2 = pTrackXi->HasSharedPointOnITSLayer(2);
        fTreeCascVarPosITSSharedClusters3 = pTrackXi->HasSharedPointOnITSLayer(3);
        fTreeCascVarPosITSSharedClusters4 = pTrackXi->HasSharedPointOnITSLayer(4);
        fTreeCascVarPosITSSharedClusters5 = pTrackXi->HasSharedPointOnITSLayer(5);
        
        fTreeCascVarNegITSSharedClusters0 = nTrackXi->HasSharedPointOnITSLayer(0);
        fTreeCascVarNegITSSharedClusters1 = nTrackXi->HasSharedPointOnITSLayer(1);
        fTreeCascVarNegITSSharedClusters2 = nTrackXi->HasSharedPointOnITSLayer(2);
        fTreeCascVarNegITSSharedClusters3 = nTrackXi->HasSharedPointOnITSLayer(3);
        fTreeCascVarNegITSSharedClusters4 = nTrackXi->HasSharedPointOnITSLayer(4);
        fTreeCascVarNegITSSharedClusters5 = nTrackXi->HasSharedPointOnITSLayer(5);
        
        fTreeCascVarBachITSSharedClusters0 = bachTrackXi->HasSharedPointOnITSLayer(0);
        fTreeCascVarBachITSSharedClusters1 = bachTrackXi->HasSharedPointOnITSLayer(1);
        fTreeCascVarBachITSSharedClusters2 = bachTrackXi->HasSharedPointOnITSLayer(2);
        fTreeCascVarBachITSSharedClusters3 = bachTrackXi->HasSharedPointOnITSLayer(3);
        fTreeCascVarBachITSSharedClusters4 = bachTrackXi->HasSharedPointOnITSLayer(4);
        fTreeCascVarBachITSSharedClusters5 = bachTrackXi->HasSharedPointOnITSLayer(5);
        
        //GetKinkIndex condition
        if( bachTrackXi->GetKinkIndex(0)>0 ) fTreeCascVarBachIsKink = kTRUE;
        if( pTrackXi->GetKinkIndex(0)>0 ) fTreeCascVarPosIsKink = kTRUE;
        if( nTrackXi->GetKinkIndex(0)>0 ) fTreeCascVarNegIsKink = kTRUE;
        
        //Get track uncertainties
        //WARNING: THIS REFERS TO THE UNCERTAINTIES CLOSEST TO THE PV
        fTreeCascVarNegDCAPVSigmaX2 = TMath::Power(TMath::Sin(nTrackXi->GetAlpha()),2)*nTrackXi->GetSigmaY2();
        fTreeCascVarNegDCAPVSigmaY2 = TMath::Power(TMath::Cos(nTrackXi->GetAlpha()),2)*nTrackXi->GetSigmaY2();
        fTreeCascVarNegDCAPVSigmaZ2 = nTrackXi->GetSigmaZ2();
        
        fTreeCascVarPosDCAPVSigmaX2 = TMath::Power(TMath::Sin(pTrackXi->GetAlpha()),2)*pTrackXi->GetSigmaY2();
        fTreeCascVarPosDCAPVSigmaY2 = TMath::Power(TMath::Cos(pTrackXi->GetAlpha()),2)*pTrackXi->GetSigmaY2();
        fTreeCascVarPosDCAPVSigmaZ2 = pTrackXi->GetSigmaZ2();
        
        fTreeCascVarBachDCAPVSigmaX2 = TMath::Power(TMath::Sin(bachTrackXi->GetAlpha()),2)*bachTrackXi->GetSigmaY2();
        fTreeCascVarBachDCAPVSigmaY2 = TMath::Power(TMath::Cos(bachTrackXi->GetAlpha()),2)*bachTrackXi->GetSigmaY2();
        fTreeCascVarBachDCAPVSigmaZ2 = bachTrackXi->GetSigmaZ2();
        
        //momenta of positive, negative and bachelor tracks
        Double_t lBMom[3], lNMom[3], lPMom[3];
        xi->GetBPxPyPz( lBMom[0], lBMom[1], lBMom[2] );
        xi->GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
        xi->GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );
        
        //fTreeCascVarBachTotMom = TMath::Sqrt( lBMom[0]*lBMom[0] + lBMom[1]*lBMom[1] + lBMom[2]*lBMom[2] );
        //fTreeCascVarPosTotMom  = TMath::Sqrt( lPMom[0]*lPMom[0] + lPMom[1]*lPMom[1] + lPMom[2]*lPMom[2] );
        //fTreeCascVarNegTotMom  = TMath::Sqrt( lNMom[0]*lNMom[0] + lNMom[1]*lNMom[1] + lNMom[2]*lNMom[2] );
        
        fTreeCascVarNegPx = lNMom[0];
        fTreeCascVarNegPy = lNMom[1];
        fTreeCascVarNegPz = lNMom[2];
        fTreeCascVarPosPx = lPMom[0];
        fTreeCascVarPosPy = lPMom[1];
        fTreeCascVarPosPz = lPMom[2];
        fTreeCascVarBachPx = lBMom[0];
        fTreeCascVarBachPy = lBMom[1];
        fTreeCascVarBachPz = lBMom[2];
        
        //Int_t lNegTrackSign = 1;
        //  Int_t lPosTrackSign = 1;
        // Int_t lBachTrackSign = 1;
        
        if( nTrackXi->GetSign() < 0 ) lNegTrackSign = -1;
        if( nTrackXi->GetSign() > 0 ) lNegTrackSign = +1;
        
        if( pTrackXi->GetSign() < 0 ) lPosTrackSign = -1;
        if( pTrackXi->GetSign() > 0 ) lPosTrackSign = +1;
        
        if( bachTrackXi->GetSign() < 0 ) lBachTrackSign = -1;
        if( bachTrackXi->GetSign() > 0 ) lBachTrackSign = +1;
        
        //------------------------------------------------
        // TPC dEdx information
        //------------------------------------------------
        fTreeCascVarNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kPion   );
        fTreeCascVarNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kProton );
        fTreeCascVarPosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kPion );
        fTreeCascVarPosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kProton );
        fTreeCascVarBachNSigmaPion  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachNSigmaKaon  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kKaon );
        
        //------------------------------------------------
        // TOF info (no correction for weak decay traj.)
        //------------------------------------------------
        fTreeCascVarNegTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( nTrackXi, AliPID::kPion   );
        fTreeCascVarNegTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( nTrackXi, AliPID::kProton );
        fTreeCascVarPosTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( pTrackXi, AliPID::kPion );
        fTreeCascVarPosTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( pTrackXi, AliPID::kProton );
        fTreeCascVarBachTOFNSigmaPion  = fPIDResponse->NumberOfSigmasTOF( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachTOFNSigmaKaon  = fPIDResponse->NumberOfSigmasTOF( bachTrackXi, AliPID::kKaon );
        
        fTreeCascVarNegITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( nTrackXi, AliPID::kPion   );
        fTreeCascVarNegITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( nTrackXi, AliPID::kProton );
        fTreeCascVarPosITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( pTrackXi, AliPID::kPion );
        fTreeCascVarPosITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( pTrackXi, AliPID::kProton );
        fTreeCascVarBachITSNSigmaPion  = fPIDResponse->NumberOfSigmasITS( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachITSNSigmaKaon  = fPIDResponse->NumberOfSigmasITS( bachTrackXi, AliPID::kKaon );
        
        //------------------------------------------------
        // Raw TPC dEdx + PIDForTracking information
        //------------------------------------------------
        
        //Step 1: Acquire TPC inner wall total momentum
        const AliExternalTrackParam *innerneg=nTrackXi->GetInnerParam();
        const AliExternalTrackParam *innerpos=pTrackXi->GetInnerParam();
        const AliExternalTrackParam *innerbach=bachTrackXi->GetInnerParam();
        fTreeCascVarPosInnerP = -1;
        fTreeCascVarNegInnerP = -1;
        fTreeCascVarBachInnerP = -1;
        if(innerpos)  { fTreeCascVarPosInnerP  = innerpos ->GetP(); }
        if(innerneg)  { fTreeCascVarNegInnerP  = innerneg ->GetP(); }
        if(innerbach) { fTreeCascVarBachInnerP = innerbach->GetP(); }
        
        //Step 2: Acquire TPC Signals
        fTreeCascVarPosdEdx = pTrackXi->GetTPCsignal();
        fTreeCascVarNegdEdx = nTrackXi->GetTPCsignal();
        fTreeCascVarBachdEdx = bachTrackXi->GetTPCsignal();
        
        //Step 3: Acquire PID For Tracking
        fTreeCascVarPosPIDForTracking = pTrackXi->GetPIDForTracking();
        fTreeCascVarNegPIDForTracking = nTrackXi->GetPIDForTracking();
        fTreeCascVarBachPIDForTracking = bachTrackXi->GetPIDForTracking();
        
        //------------------------------------------------
        // TPC Number of clusters info
        // --- modified to save the smallest number
        // --- of TPC clusters for the 3 tracks
        //------------------------------------------------
        
        lPosTPCClusters   = pTrackXi->GetTPCNcls();
        lNegTPCClusters   = nTrackXi->GetTPCNcls();
        lBachTPCClusters  = bachTrackXi->GetTPCNcls();
        
        // 1 - Poor quality related to TPCrefit
        ULong_t pStatus    = pTrackXi->GetStatus();
        ULong_t nStatus    = nTrackXi->GetStatus();
        ULong_t bachStatus = bachTrackXi->GetStatus();
        
        //fTreeCascVarkITSRefitBachelor = kTRUE;
        //fTreeCascVarkITSRefitNegative = kTRUE;
        //fTreeCascVarkITSRefitPositive = kTRUE;
        
        if ((pStatus&AliESDtrack::kTPCrefit)    == 0) {
            AliDebug(1, "Pb / V0 Pos. track has no TPCrefit ... continue!");
            continue;
        }
        if ((nStatus&AliESDtrack::kTPCrefit)    == 0) {
            AliDebug(1, "Pb / V0 Neg. track has no TPCrefit ... continue!");
            continue;
        }
        if ((bachStatus&AliESDtrack::kTPCrefit) == 0) {
            AliDebug(1, "Pb / Bach.   track has no TPCrefit ... continue!");
            continue;
        }
        
        //Get status flags
        fTreeCascVarPosTrackStatus = pTrackXi->GetStatus();
        fTreeCascVarNegTrackStatus = nTrackXi->GetStatus();
        fTreeCascVarBachTrackStatus = bachTrackXi->GetStatus();
        
        //Sandbox mode information, please
        fTreeCascVarMagneticField = lESDevent->GetMagneticField();
        fTreeCascVarPosDCAz = GetDCAz(pTrackXi);
        fTreeCascVarNegDCAz = GetDCAz(nTrackXi);
        fTreeCascVarBachDCAz = GetDCAz(bachTrackXi);
        
        Float_t lPosChi2PerCluster = pTrackXi->GetTPCchi2() / ((Float_t) lPosTPCClusters);
        Float_t lNegChi2PerCluster = nTrackXi->GetTPCchi2() / ((Float_t) lNegTPCClusters);
        Float_t lBachChi2PerCluster = bachTrackXi->GetTPCchi2() / ((Float_t) lBachTPCClusters);
        
        Int_t leastnumberofclusters = 1000;
        Float_t lBiggestChi2PerCluster = -1;
        
        //Pick minimum
        if( lPosTPCClusters < leastnumberofclusters ) leastnumberofclusters = lPosTPCClusters;
        if( lNegTPCClusters < leastnumberofclusters ) leastnumberofclusters = lNegTPCClusters;
        if( lBachTPCClusters < leastnumberofclusters ) leastnumberofclusters = lBachTPCClusters;
        
        //Pick maximum
        if( lPosChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lPosChi2PerCluster;
        if( lNegChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lNegChi2PerCluster;
        if( lBachChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lBachChi2PerCluster;
        
        //Extra track quality: min track length
        Float_t lSmallestTrackLength = 1000;
        Float_t lPosTrackLength = -1;
        Float_t lNegTrackLength = -1;
        Float_t lBachTrackLength = -1;
        
        if (pTrackXi->GetInnerParam()) lPosTrackLength = pTrackXi->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if (nTrackXi->GetInnerParam()) lNegTrackLength = nTrackXi->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if (bachTrackXi->GetInnerParam()) lBachTrackLength = bachTrackXi->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        
        if ( lPosTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lPosTrackLength;
        if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
        if ( lBachTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lBachTrackLength;
        
        fTreeCascVarMinTrackLength = lSmallestTrackLength;
        
        // 2 - Poor quality related to TPC clusters: lowest cut of 70 clusters
        if(lPosTPCClusters  < 70 && lSmallestTrackLength < 80 && fkExtraCleanup) {
            AliDebug(1, "Pb / V0 Pos. track has less than 70 TPC clusters ... continue!");
            continue;
        }
        if(lNegTPCClusters  < 70 && lSmallestTrackLength < 80 && fkExtraCleanup) {
            AliDebug(1, "Pb / V0 Neg. track has less than 70 TPC clusters ... continue!");
            continue;
        }
        if(lBachTPCClusters < 70 && lSmallestTrackLength < 80 && fkExtraCleanup) {
            AliDebug(1, "Pb / Bach.   track has less than 70 TPC clusters ... continue!");
            continue;
        }
        
        lInvMassLambdaAsCascDghter    = xi->GetEffMass();
        // This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
        lDcaV0DaughtersXi         = xi->GetDcaV0Daughters();
        //lV0Chi2Xi             = xi->GetChi2V0();
        
        
        lV0CosineOfPointingAngleXi     = xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                        lBestPrimaryVtxPos[1],
                                                                        lBestPrimaryVtxPos[2] );
        //Modification: V0 CosPA wrt to Cascade decay vertex
        lV0CosineOfPointingAngleXiSpecial     = xi->GetV0CosineOfPointingAngle( lPosXi[0],
                                                                               lPosXi[1],
                                                                               lPosXi[2] );
        
        lDcaV0ToPrimVertexXi  = xi->GetD( lBestPrimaryVtxPos[0],
                                         lBestPrimaryVtxPos[1],
                                         lBestPrimaryVtxPos[2]);
        
        lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD(    lBestPrimaryVtxPos[0],
                                                               lBestPrimaryVtxPos[1],
                                                               lMagneticField  ) );
        
        //xy and z impact for the bach.track
        Float_t xy_z_bach_to_prim[2];
        bachTrackXi->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2],lMagneticField, xy_z_bach_to_prim);
        lDcaBachToPrimVertexXi =xy_z_bach_to_prim[0];
        lDcaBachToPrimVertexXiZ = xy_z_bach_to_prim[1];
        
        // Note : AliExternalTrackParam::GetD returns an algebraic value ...
        
        xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] );
        lV0RadiusXi    = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
        
        fTreeCascVarPrimVertexX = lBestPrimaryVtxPos[0];
        fTreeCascVarPrimVertexY = lBestPrimaryVtxPos[1];
        fTreeCascVarPrimVertexZ = lBestPrimaryVtxPos[2];
        
        fTreeCascVarV0DecayX = lPosV0Xi[0];
        fTreeCascVarV0DecayY = lPosV0Xi[1];
        fTreeCascVarV0DecayZ = lPosV0Xi[2];
        
        //========================================================================================
        //Calculate V0 lifetime for adaptive decay radius cut
        //3D Distance travelled by the V0 in the cascade
        Float_t lV0DistanceTrav =  TMath::Sqrt(  TMath::Power( lPosV0Xi[0]-lPosXi[0] , 2)
                                               + TMath::Power( lPosV0Xi[1]-lPosXi[1] , 2)
                                               + TMath::Power( lPosV0Xi[2]-lPosXi[2] , 2) );
        
        //Total V0 momentum
        Float_t lV0TotMomentum = TMath::Sqrt(  TMath::Power( lNMom[0]+lPMom[0] , 2)
                                             + TMath::Power( lNMom[1]+lPMom[1] , 2)
                                             + TMath::Power( lNMom[2]+lPMom[2] , 2) );
        
        //
        //  Double_t lV0TotMomentumVec[3] = {lNMom[0]+lPMom[0], lNMom[1]+lPMom[1], lNMom[2]+lPMom[2]};
        
        //V0 transverse momentum
        Float_t lV0Pt = TMath::Sqrt(  TMath::Power( lNMom[0]+lPMom[0] , 2)
                                    + TMath::Power( lNMom[1]+lPMom[1] , 2) );
        
        //Calculate V0 lifetime: mL/p
        if( TMath::Abs(lV0TotMomentum)>1e-5 ){
            fTreeCascVarV0Lifetime = 1.115683*lV0DistanceTrav / lV0TotMomentum;
        }else{
            fTreeCascVarV0Lifetime = -1;
        }
        //========================================================================================
        //Cowboy/sailor info regarding V0 inside cascade
        //Calculate vec prod with momenta projected to xy plane
        //Provisions for cowboy/sailor check
        Double_t lModp1 = TMath::Sqrt( lPMom[0]*lPMom[0] + lPMom[1]*lPMom[1] );
        Double_t lModp2 = TMath::Sqrt( lNMom[0]*lNMom[0] + lNMom[1]*lNMom[1] );
        
        Double_t lVecProd = (lPMom[0]*lNMom[1] - lPMom[1]*lNMom[0]) / (lModp1*lModp2);
        
        if ( lMagneticField < 0 ) lVecProd *= -1; //invert sign
        
        fTreeCascVarIsCowboy = kFALSE;
        if (lVecProd < 0) fTreeCascVarIsCowboy = kTRUE;
        
        fTreeCascVarCowboyness = lVecProd;
        
        Double_t lBachMod = TMath::Sqrt(lBMom[0]*lBMom[0]+lBMom[1]*lBMom[1]);
        Double_t lV0px = lPMom[0] + lNMom[0];
        Double_t lV0py = lPMom[1] + lNMom[1];
        Double_t lVecProdXi = (lV0px*lBMom[1] - lV0py*lBMom[0]) / (lV0Pt*lBachMod);
        
        if ( lMagneticField < 0 ) lVecProdXi *= -1; //invert sign
        
        fTreeCascVarIsCascadeCowboy = kFALSE;
        if (lVecProdXi < 0) fTreeCascVarIsCascadeCowboy = kTRUE;
        
        fTreeCascVarCascadeCowboyness = lVecProdXi;
        //========================================================================================
        
        lDcaPosToPrimVertexXi     = TMath::Abs( pTrackXi    ->GetD(    lBestPrimaryVtxPos[0],
                                                                   lBestPrimaryVtxPos[1],
                                                                   lMagneticField  )     );
        
        lDcaNegToPrimVertexXi     = TMath::Abs( nTrackXi    ->GetD(    lBestPrimaryVtxPos[0],
                                                                   lBestPrimaryVtxPos[1],
                                                                   lMagneticField  )     );
        
        //----------when transverse[0] and z-impact parameter[1] needed------------------------
        Float_t xy_z_pos_to_prim[2];
        pTrackXi->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2],lMagneticField, xy_z_pos_to_prim);
        lDcaPosToPrimVertexXi = xy_z_pos_to_prim[0];
        lDcaPosToPrimVertexXiZ = xy_z_pos_to_prim[1];
        //-------------------------------------------------------------------------------------
        
        //xy and z impact parameters for the neg.track
        Float_t xy_z_neg_to_prim[2];
        nTrackXi->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2],lMagneticField, xy_z_neg_to_prim);
        lDcaNegToPrimVertexXi = xy_z_neg_to_prim[0];
        lDcaNegToPrimVertexXiZ = xy_z_neg_to_prim[1];
        
        
        // - II.Step 4 : around effective masses (ESD)
        // ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
        
        if( bachTrackXi->Charge() < 0 )    {
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , 3312);
            // Calculate the effective mass of the Xi- candidate.
            // pdg code 3312 = Xi-
            lInvMassXiMinus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , 3334);
            // Calculate the effective mass of the Xi- candidate.
            // pdg code 3334 = Omega-
            lInvMassOmegaMinus = xi->GetEffMassXi();
            
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , 3312);     // Back to default hyp.
        }// end if negative bachelor
        
        
        if( bachTrackXi->Charge() >  0 ) {
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , -3312);
            // Calculate the effective mass of the Xi+ candidate.
            // pdg code -3312 = Xi+
            lInvMassXiPlus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , -3334);
            // Calculate the effective mass of the Xi+ candidate.
            // pdg code -3334  = Omega+
            lInvMassOmegaPlus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , -3312);     // Back to "default" hyp.
        }// end if positive bachelor
        
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        //Recalculate from scratch, please
        //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
        //WARNING: This will not be checked for correctness (charge-wise, etc)
        //         It will be up to the user to use the correct variable whenever needed!
        
        //+-+ Recalculate Xi Masses from scratch: will not change with lambda mass as
        //the perfect lambda mass is always assumed
        //   Double_t m_proton = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();
        //   Double_t m_pion = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
        //   Double_t m_kaon = TDatabasePDG::Instance()->GetParticle(321)->Mass(); //K+
        //  cout << "mass of proton / pion / kaon is = " << m_proton << " , " << m_pion << " , " << m_kaon << " , " << endl;
        // >  mass of proton / pion / kaon is = 0.938272 , 0.13957 , 0.493677
        
        Double_t lInvMassOmegaMinusScratch = 0.;
        Double_t lInvMassOmegaPlusScratch = 0.;
        Double_t lInvMassXiMinusScratch = 0.;
        Double_t lInvMassXiPlusScratch = 0.;
        
        // - II.Step 6 : extra info for QA (ESD)
        // miscellaneous pieces of info that may help regarding data quality assessment.
        //-------------
        xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
        //Needed for the calculation of GetCosOfProtonLambdaRestOmegaRest
        // lXiMomVec[0] = lXiMomX; lXiMomVec[1] = lXiMomY; lXiMomVec[2] = lXiMomZ;
        //---------------------
        lXiTransvMom      = TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
        lXiTotMom      = TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
        
        xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
        //lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
        //lBachTotMom      = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );
        
        lChargeXi = xi->Charge();
        
        //lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
        
        lRapXi    = xi->RapXi();
        lRapOmega = xi->RapOmega();
        //lEta      = xi->Eta();
        //lTheta    = xi->Theta() *180.0/TMath::Pi();
        //lPhi      = xi->Phi()   *180.0/TMath::Pi();
        lAlphaXi  = xi->AlphaXi();
        lPtArmXi  = xi->PtArmXi();
        //  cout << "lAlphaXi = " << lAlphaXi << endl;
        
        //----------------------------------------
        // Calculate Cascade DCA to PV, please
        //----------------------------------------
        
        Int_t lChargeCascade = fTreeCascVarCharge;
        
        //cascade properties to get started
        Double_t xyzCascade[3], pxpypzCascade[3], cvCascade[21];
        for(Int_t ii=0;ii<21;ii++) cvCascade[ii]=0.0; //something small
        
        xi->GetXYZcascade( xyzCascade[0],  xyzCascade[1], xyzCascade[2] );
        xi->GetPxPyPz( pxpypzCascade[0], pxpypzCascade[1], pxpypzCascade[2] );
        
        AliExternalTrackParam lCascTrajObject(xyzCascade,pxpypzCascade,cvCascade,lChargeCascade), *hCascTraj = &lCascTrajObject;
        
        Double_t lCascDCAtoPVxy = TMath::Abs(hCascTraj->GetD(lBestPrimaryVtxPos[0],
                                                             lBestPrimaryVtxPos[1],
                                                             lMagneticField) );
        Float_t dzcascade[2];
        hCascTraj->GetDZ(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2], lMagneticField, dzcascade );
        Double_t lCascDCAtoPVz = dzcascade[1];
        
        //assign TTree values
        fTreeCascVarCascDCAtoPVxy = lCascDCAtoPVxy;
        fTreeCascVarCascDCAtoPVz  = lCascDCAtoPVz;
        
        //----------------------------------------
        // Bump studies: perform propagation
        //----------------------------------------
        
        AliESDtrack *lBaryonTrack = 0x0;
        AliESDtrack *lBachelorTrack = 0x0;
        if ( lChargeXi == -1 ){
            lBaryonTrack = pTrackXi;
            lBachelorTrack = bachTrackXi;
        }
        if ( lChargeXi == +1 ){
            lBaryonTrack = nTrackXi;
            lBachelorTrack = bachTrackXi;
        }
        
        fTreeCascVarDCABachToBaryon = -100;
        
        Double_t bMag = lESDevent->GetMagneticField();
        Double_t xn, xp;
        
        //Care has to be taken here
        if ( lBaryonTrack && lBachelorTrack ){
            //Attempt zero: Calculate DCA between bachelor and baryon daughter
            fTreeCascVarDCABachToBaryon = lBaryonTrack->GetDCA(lBachelorTrack, bMag, xn, xp);
        }
        
        fTreeCascVarWrongCosPA = -1;
        if( bachTrackXi->Charge() < 0 )
            fTreeCascVarWrongCosPA = GetCosPA( bachTrackXi , pTrackXi, lESDevent );
        if( bachTrackXi->Charge() > 0 )
            fTreeCascVarWrongCosPA = GetCosPA( bachTrackXi , nTrackXi, lESDevent );
        
        
        //------------------------------------------------
        // Set Variables for adding to tree
        //------------------------------------------------
        
        
        //The member GetCosOfProtonLambdaRestOmegaRest cannot work with pointers, so we need to define new variables
        //with assigned values of pointers
        //  Double_t output[3] = {0.,0.,0.};
        
        fTreeCascVarMVPileupFlag = fMVPileupFlag;
        
        fTreeCascVarCharge    = lChargeXi;
        if (lChargeXi < 0 ){
            //determined with AliCascader
            fTreeCascVarMassAsXi    = lInvMassXiMinus;
            fTreeCascVarMassAsOmega = lInvMassOmegaMinus;
            //determined from scratch using the info of the charged tracks
            fTreeCascVarMassOmegaScratch = lInvMassOmegaMinusScratch;
            fTreeCascVarMassXiScratch = lInvMassXiMinusScratch;
            //combinatorial backgrounds
            fTreeCascVarBackgroundMassOmega = lInvMassOmegaMinusBackground;
            fTreeCascVarBackgroundMassXi = lInvMassXiMinusBackground;
            
            cosAAA = 10;
            //  cosAAA = GetCosOfProtonLambdaRestOmegaRest(pXi, fTreeCascVarMassAsOmega, pL, fTreeCascVarV0MassLambda, pp, m_proton);
            // cout << "The cosine value is = " << cosAAA << endl;
            fTreeCascVarCosineProtonRestV0RestOmegaMinus = cosAAA;
            // dzeta = DzetaFromMomenta(pXi, pL, pp);
            fTreeCascVarDzetaFromMomentaOmegaMinus = dzeta;
            // cout << dzeta << endl;
            //  output[0] = GetThetaOfDaughterInRestOfMotherWorld(pXi, fTreeCascVarMassAsOmega, pL, fTreeCascVarV0MassLambda)[0];
            //  cout << "[phi, omega] = [ " << output[0]*180/TMath::Pi() << " , " << output[1]*180/TMath::Pi() << " ]" << endl;
        }
        if (lChargeXi > 0 ){
            fTreeCascVarMassAsXi    = lInvMassXiPlus;
            fTreeCascVarMassAsOmega = lInvMassOmegaPlus;
            //determined from scratch using the info of the charged tracks
            fTreeCascVarMassOmegaScratch = lInvMassOmegaPlusScratch;
            fTreeCascVarMassXiScratch = lInvMassXiPlusScratch;
            //include combinatorial background for the study of shape
            fTreeCascVarBackgroundMassOmega = lInvMassOmegaPlusBackground;
            fTreeCascVarBackgroundMassXi = lInvMassXiPlusBackground;
            
            
            //anti-proton is handled as a negative track
            cosAAA = 10;
            //  cosAAA = GetCosOfProtonLambdaRestOmegaRest(pXi, fTreeCascVarMassAsOmega, pL, fTreeCascVarV0MassAntiLambda, pn, m_proton);
            fTreeCascVarCosineProtonRestV0RestOmegaPlus = cosAAA;
        }
        fTreeCascVarPt = lXiTransvMom;
        fTreeCascVarP = lXiTotMom;
        fTreeCascVarPz = lXiMomZ;
        fTreeCascVarRapXi = lRapXi ;
        fTreeCascVarRapOmega = lRapOmega ;
        fTreeCascVarDCACascDaughters = lDcaXiDaughters;
        fTreeCascVarDCABachToPrimVtx = lDcaBachToPrimVertexXi;
        fTreeCascVarDCABachToPrimVtxZ = lDcaBachToPrimVertexXiZ; //additionaly the z component
        fTreeCascVarDCAV0Daughters = lDcaV0DaughtersXi;
        fTreeCascVarDCAV0ToPrimVtx = lDcaV0ToPrimVertexXi;
        fTreeCascVarDCAPosToPrimVtx = lDcaPosToPrimVertexXi;
        fTreeCascVarDCAPosToPrimVtxZ = lDcaPosToPrimVertexXiZ;
        fTreeCascVarDCANegToPrimVtx = lDcaNegToPrimVertexXi;
        fTreeCascVarDCANegToPrimVtxZ = lDcaNegToPrimVertexXiZ;
        fTreeCascVarCascCosPointingAngle = lXiCosineOfPointingAngle* bachTrackXi->Charge(); //so that we do not an extra branch in the tree.
        fTreeCascVarCascRadius = lXiRadius;
        fTreeCascVarV0Mass = lInvMassLambdaAsCascDghter;
        fTreeCascVarV0CosPointingAngle = lV0CosineOfPointingAngleXi ;
        fTreeCascVarV0CosPointingAngleSpecial = lV0CosineOfPointingAngleXiSpecial;
        fTreeCascVarV0Radius = lV0RadiusXi;
        fTreeCascVarLeastNbrClusters = leastnumberofclusters;
        fTreeCascVarMaxChi2PerCluster = lBiggestChi2PerCluster;
        
        //Copy Multiplicity information
        fTreeCascVarCentrality = fCentrality;
        
        fTreeCascVarDistOverTotMom = TMath::Sqrt(
                                                 TMath::Power( lPosXi[0] - lBestPrimaryVtxPos[0] , 2) +
                                                 TMath::Power( lPosXi[1] - lBestPrimaryVtxPos[1] , 2) +
                                                 TMath::Power( lPosXi[2] - lBestPrimaryVtxPos[2] , 2)
                                                 );
        fTreeCascVarDistOverTotMom /= (lXiTotMom+1e-13)*bachTrackXi->Charge(); //so that we dont need an extra branch for the cascade charge
        
        //Info for pileup studies
        fTreeCascVarBachTOFExpTDiff = bachTrackXi->GetTOFExpTDiff( bMag );
        fTreeCascVarNegTOFExpTDiff = nTrackXi->GetTOFExpTDiff( bMag );
        fTreeCascVarPosTOFExpTDiff = pTrackXi->GetTOFExpTDiff( bMag );
        fTreeCascVarBachTOFSignal = bachTrackXi->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarNegTOFSignal = nTrackXi->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarPosTOFSignal = pTrackXi->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarBachTOFBCid = bachTrackXi->GetTOFBunchCrossing( bMag );
        fTreeCascVarNegTOFBCid = nTrackXi->GetTOFBunchCrossing( bMag );
        fTreeCascVarPosTOFBCid = pTrackXi->GetTOFBunchCrossing( bMag );
        //Copy OOB pileup flag for this event
        fTreeCascVarOOBPileupFlag = fOOBPileupFlag;
        //Copy VZERO information for this event
        fTreeCascVarAmplitudeV0A = fAmplitudeV0A;
        fTreeCascVarAmplitudeV0C = fAmplitudeV0C;
        //Copy IR information for this event
        fTreeCascVarClosestNonEmptyBC = fClosestNonEmptyBC;
        
        
        //set to test the outcome 14.02.19
        // fkPreselectDedx = kTRUE;
        
        //Armenteros-Podolanski for the Cascades.
        fTreeVariableAlphaXi =lAlphaXi;
        fTreeVariablePtArmXi = lPtArmXi;
        
        
        Bool_t lPassesPreFilterdEdx = kFALSE;
        
        if ( fkExtraCleanup )
        { //set to kTRUE
            //Valid or not valid
            //Meant to provide extra level of cleanup
            if( TMath::Abs(fTreeCascVarPosEta)>0.8 || TMath::Abs(fTreeCascVarNegEta)>0.8 || TMath::Abs(fTreeCascVarBachEta)>0.8 ) return;
            if( TMath::Abs(fTreeCascVarRapXi)>0.5 && TMath::Abs(fTreeCascVarRapOmega)>0.5 ) return;
            if ( fkPreselectDedx ){ //set to kTRUE
                Double_t lWindow = 0.11;
                //XiMinus Pre-selection
                if( fTreeCascVarMassAsXi<1.322+lWindow&&fTreeCascVarMassAsXi>1.322-lWindow &&
                   TMath::Abs(fTreeCascVarPosNSigmaProton) < 4.0 &&
                   TMath::Abs(fTreeCascVarNegNSigmaPion) < 4.0 &&
                   TMath::Abs(fTreeCascVarBachNSigmaPion) < 4.0 &&
                   fTreeCascVarCharge == -1)
                {
                    if (( (!CheckTOFstatus(pTrackXi)) ||  (CheckTOFstatus(pTrackXi) && TMath::Abs(fTreeCascVarPosTOFNSigmaProton)<4.0)  ) &&
                        ( (!CheckTOFstatus(nTrackXi)) ||  (CheckTOFstatus(nTrackXi) && TMath::Abs(fTreeCascVarNegTOFNSigmaPion)<4.0)  ) &&
                        ( (!CheckTOFstatus(bachTrackXi)) ||  (CheckTOFstatus(bachTrackXi) && TMath::Abs(fTreeCascVarBachTOFNSigmaPion)<4.0) ))
                    {
                        
                        lPassesPreFilterdEdx = kTRUE;
                    }
                    
                }
                else if( fTreeCascVarMassAsXi<1.322+lWindow&&fTreeCascVarMassAsXi>1.322-lWindow &&
                        TMath::Abs(fTreeCascVarPosNSigmaPion) < 4.0 &&
                        TMath::Abs(fTreeCascVarNegNSigmaProton) < 4.0 &&
                        TMath::Abs(fTreeCascVarBachNSigmaPion) < 4.0 &&
                        fTreeCascVarCharge == +1 )
                {
                    if (( (!CheckTOFstatus(pTrackXi)) ||  (CheckTOFstatus(pTrackXi) && TMath::Abs(fTreeCascVarPosTOFNSigmaPion)<4.0)  ) &&
                        ( (!CheckTOFstatus(nTrackXi)) ||  (CheckTOFstatus(nTrackXi) && TMath::Abs(fTreeCascVarNegTOFNSigmaProton)<4.0)  ) &&
                        ( (!CheckTOFstatus(bachTrackXi)) ||  (CheckTOFstatus(bachTrackXi) && TMath::Abs(fTreeCascVarBachTOFNSigmaPion)<4.0) ))
                    {
                        
                        lPassesPreFilterdEdx = kTRUE;
                    }
                    
                }
                else if(fTreeCascVarMassAsOmega<1.672+lWindow&&fTreeCascVarMassAsOmega>1.672-lWindow &&
                        TMath::Abs(fTreeCascVarPosNSigmaProton) < 4.0 &&
                        TMath::Abs(fTreeCascVarNegNSigmaPion) < 4.0 &&
                        TMath::Abs(fTreeCascVarBachNSigmaKaon) < 4.0 &&
                        fTreeCascVarCharge == -1  )
                {
                    
                    if (( (!CheckTOFstatus(pTrackXi)) ||  (CheckTOFstatus(pTrackXi) && TMath::Abs(fTreeCascVarPosTOFNSigmaProton)<4.0)  ) &&
                        ( (!CheckTOFstatus(nTrackXi)) ||  (CheckTOFstatus(nTrackXi) && TMath::Abs(fTreeCascVarNegTOFNSigmaPion)<4.0)  ) &&
                        ( (!CheckTOFstatus(bachTrackXi)) ||  (CheckTOFstatus(bachTrackXi) && TMath::Abs(fTreeCascVarBachTOFNSigmaKaon)<4.0) ))
                    {
                        
                        lPassesPreFilterdEdx = kTRUE;
                    }
                }
                else if(fTreeCascVarMassAsOmega<1.672+lWindow&&fTreeCascVarMassAsOmega>1.672-lWindow &&
                        TMath::Abs(fTreeCascVarPosNSigmaPion) < 4.0 &&
                        TMath::Abs(fTreeCascVarNegNSigmaProton) < 4.0 &&
                        TMath::Abs(fTreeCascVarBachNSigmaKaon) < 4.0 &&
                        fTreeCascVarCharge == +1)
                {
                    
                    if (( (!CheckTOFstatus(pTrackXi)) ||  (CheckTOFstatus(pTrackXi) && TMath::Abs(fTreeCascVarPosTOFNSigmaPion)<4.0)  ) &&
                        ( (!CheckTOFstatus(nTrackXi)) ||  (CheckTOFstatus(nTrackXi) && TMath::Abs(fTreeCascVarNegTOFNSigmaProton)<4.0)  ) &&
                        ( (!CheckTOFstatus(bachTrackXi)) ||  (CheckTOFstatus(bachTrackXi) && TMath::Abs(fTreeCascVarBachTOFNSigmaKaon)<4.0) ))
                    {
                        lPassesPreFilterdEdx = kTRUE;
                    }
                    
                }
                else lPassesPreFilterdEdx = kFALSE;
                
            }
            
        }
        
        
        //All vars not specified here: specified elsewhere!
        
        //------------------------------------------------
        // Fill Tree!
        //------------------------------------------------
        
        // The conditional is meant to decrease excessive
        // memory usage! Be careful when loosening the
        // cut!
        
        //Xi    Mass window: 150MeV wide
        //Omega mass window: 150MeV wide
        
        //Random denial
        Bool_t lKeepCascade = kTRUE;
        
        //selection window, 100 MeV
        Bool_t massXiRange = kFALSE;
        if((fTreeCascVarMassAsXi<1.32+0.050&&fTreeCascVarMassAsXi>1.32-0.050)||(fTreeCascVarMassAsOmega<1.68+0.050&&fTreeCascVarMassAsOmega>1.68-0.050)){
            massXiRange = kTRUE;}
        
        //status of the chi2 of V0 and Xi from aliESDcascade
        Bool_t chi2LF = kTRUE;
        //0  and 1024 - predefined by aliesdcascade
        if((fTreeCascVarChiSquareV0 != 0) || (fTreeCascVarChiSquareCascade !=1024)) chi2LF = kFALSE;
        
        
        TLorentzVector LVpos;
        TLorentzVector LVneg;
        TLorentzVector LVbach;
        LVpos.SetPxPyPzE(pTrackXi->Px(), pTrackXi->Py(), pTrackXi->Pz(), pTrackXi->E());
        LVneg.SetPxPyPzE(nTrackXi->Px(), nTrackXi->Py(), nTrackXi->Pz(), nTrackXi->E());
        LVbach.SetPxPyPzE(bachTrackXi->Px(), bachTrackXi->Py(), bachTrackXi->Pz(), bachTrackXi->E());
        
        
        Cascade_Track = Cascade_Event -> createTrack();
        
        if(lPassesPreFilterdEdx &&
           lKeepCascade &&
           massXiRange &&
           goodnTrackXi &&
           goodpTrackXi &&
           goodbachTrackXi &&
           chi2LF //&&
           //NsigmaALL
           )
        {
            //lorentz vectors
            Cascade_Track -> set_TLV_pos(LVpos);
            Cascade_Track -> set_TLV_neg(LVneg);
            Cascade_Track -> set_TLV_bach(LVbach);
            //dcas
            Cascade_Track -> set_dca_pos_to_prim(fTreeCascVarDCAPosToPrimVtx, fTreeCascVarDCAPosToPrimVtxZ);
            Cascade_Track -> set_dca_neg_to_prim(fTreeCascVarDCANegToPrimVtx, fTreeCascVarDCANegToPrimVtxZ);
            Cascade_Track -> set_dca_bach_to_prim(fTreeCascVarDCABachToPrimVtx, fTreeCascVarDCABachToPrimVtxZ);
            Cascade_Track -> set_dca_V0_to_prim(fTreeCascVarDCAV0ToPrimVtx);
            Cascade_Track -> set_dca_Omega_to_prim(fTreeCascVarCascDCAtoPVxy, fTreeCascVarCascDCAtoPVz);
            Cascade_Track -> set_dca_pos_to_neg(fTreeCascVarDCAV0Daughters);
            Cascade_Track -> set_dca_bach_to_Lambda(fTreeCascVarDCACascDaughters);
            Cascade_Track -> set_dca_bach_to_pos(fTreeCascVarDCABachToBaryon);
            //dEdx TPC and TOF
            Cascade_Track -> set_nSigma_dEdx_pos(fTreeCascVarPosNSigmaProton, fTreeCascVarPosNSigmaPion); //[0] Casc charge = -1, [1] Casc +
            Cascade_Track -> set_nSigma_dEdx_neg(fTreeCascVarNegNSigmaPion,  fTreeCascVarNegNSigmaProton);
            Cascade_Track -> set_nSigma_dEdx_bach(fTreeCascVarBachNSigmaKaon, fTreeCascVarBachNSigmaPion);
            Cascade_Track -> set_nSigma_TOF_pos(fTreeCascVarPosTOFNSigmaProton, fTreeCascVarPosTOFNSigmaPion);
            Cascade_Track -> set_nSigma_TOF_neg(fTreeCascVarNegTOFNSigmaPion, fTreeCascVarNegTOFNSigmaProton);
            Cascade_Track -> set_nSigma_TOF_bach(fTreeCascVarBachTOFNSigmaKaon, fTreeCascVarBachTOFNSigmaPion);
            Cascade_Track -> set_CosPointingAngle(fTreeCascVarCascCosPointingAngle);
            //fill with the positions of the V0 and Xi.
            Cascade_Track -> set_CascadeDecayPos(fTreeCascVarCascadeDecayX, fTreeCascVarCascadeDecayY, fTreeCascVarCascadeDecayZ);
            Cascade_Track -> set_V0fromCascadePos(fTreeCascVarV0DecayX, fTreeCascVarV0DecayY, fTreeCascVarV0DecayZ);
            
            fTreeCascadeAsEvent ->Fill();
        }
        
        //   cout << "filling the tree!" << endl;
        
        // fTreeCascade->Fill();
        
    }// end of the Cascade loop (ESD or AOD)
    
    
    
    PostData(1, fTreeCascadeAsEvent);
    // PostData(2, fTreeCascade);
    
    
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityRun2 : ouput data container list not available\n");
        return;
    }
    
    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityRun2 : fHistEventCounter not available");
        return;
    }
    
    TCanvas *canCheck = new TCanvas("AliAnalysisTaskStrangenessVsMultiplicityRun2","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();
    
    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//________________________________________________________________________
Double_t AliAnalysisTaskStrangeCascadesDiscrete::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::AddConfiguration( AliV0Result *lV0Result )
{
    if ( lV0Result->GetMassHypothesis () == AliV0Result::kK0Short ){
        if (!fListK0Short){
            Printf("fListK0Short does not exist. Creating...");
            fListK0Short = new TList();
            fListK0Short->SetOwner();
        }
        fListK0Short->Add(lV0Result);
    }
    if ( lV0Result->GetMassHypothesis () == AliV0Result::kLambda ){
        if (!fListLambda){
            Printf("fListLambda does not exist. Creating...");
            fListLambda = new TList();
            fListLambda->SetOwner();
        }
        fListLambda->Add(lV0Result);
    }
    if ( lV0Result->GetMassHypothesis () == AliV0Result::kAntiLambda ){
        if (!fListAntiLambda){
            Printf("fListAntiLambda does not exist. Creating...");
            fListAntiLambda = new TList();
            fListAntiLambda->SetOwner();
        }
        fListAntiLambda->Add(lV0Result);
    }
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::AddConfiguration( AliCascadeResult *lCascadeResult )
{
    if ( lCascadeResult->GetMassHypothesis () == AliCascadeResult::kXiMinus ){
        if (!fListXiMinus){
            Printf("fListXiMinus does not exist. Creating...");
            fListXiMinus = new TList();
            fListXiMinus->SetOwner();
        }
        fListXiMinus->Add(lCascadeResult);
    }
    if ( lCascadeResult->GetMassHypothesis () == AliCascadeResult::kXiPlus ){
        if (!fListXiPlus){
            Printf("fListXiPlus does not exist. Creating...");
            fListXiPlus = new TList();
            fListXiPlus->SetOwner();
        }
        fListXiPlus->Add(lCascadeResult);
    }
    if ( lCascadeResult->GetMassHypothesis () == AliCascadeResult::kOmegaMinus ){
        if (!fListOmegaMinus){
            Printf("fListOmegaMinus does not exist. Creating...");
            fListOmegaMinus = new TList();
            fListOmegaMinus->SetOwner();
        }
        fListOmegaMinus->Add(lCascadeResult);
    }
    if ( lCascadeResult->GetMassHypothesis () == AliCascadeResult::kOmegaPlus ){
        if (!fListOmegaPlus){
            Printf("fListOmegaPlus does not exist. Creating...");
            fListOmegaPlus = new TList();
            fListOmegaPlus->SetOwner();
        }
        fListOmegaPlus->Add(lCascadeResult);
    }
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::SetupStandardVertexing()
//Meant to store standard re-vertexing configuration
{
    //Tell the task to re-run vertexers
    SetRunVertexers(kTRUE);
    SetDoV0Refit(kTRUE);
    
    //V0-Related topological selections
    SetV0VertexerDCAFirstToPV(0.05);
    SetV0VertexerDCASecondtoPV(0.05);
    SetV0VertexerDCAV0Daughters(1.20);
    SetV0VertexerCosinePA(0.98);
    SetV0VertexerMinRadius(0.9);
    SetV0VertexerMaxRadius(200);
    
    //Cascade-Related topological selections
    SetCascVertexerMinV0ImpactParameter(0.05);
    SetCascVertexerV0MassWindow(0.006);
    SetCascVertexerDCABachToPV(0.02);
    SetCascVertexerDCACascadeDaughters(1.2);
    SetCascVertexerCascadeMinRadius(.8);
    SetCascVertexerCascadeCosinePA(.98);
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::SetupLooseVertexing()
//Meant to store standard re-vertexing configuration
{
    //Tell the task to re-run vertexers
    SetRunVertexers(kTRUE);
    SetDoV0Refit(kTRUE);
    
    //V0-Related topological selections
    SetV0VertexerDCAFirstToPV(0.1);
    SetV0VertexerDCASecondtoPV(0.1);
    SetV0VertexerDCAV0Daughters(1.40);
    SetV0VertexerCosinePA(0.95);
    SetV0VertexerMinRadius(0.9);
    SetV0VertexerMaxRadius(200);
    
    //Cascade-Related topological selections
    SetCascVertexerMinV0ImpactParameter(0.05);
    SetCascVertexerV0MassWindow(0.006);
    SetCascVertexerDCABachToPV(0.02);
    SetCascVertexerDCACascadeDaughters(1.4);
    SetCascVertexerCascadeMinRadius(.5);
    SetCascVertexerCascadeCosinePA(.95);
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::AddTopologicalQAV0(Int_t lRecNumberOfSteps)
//Add all configurations to do QA of topological variables for the V0 analysis
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimits[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
        2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 15};
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    
    //lPtbinlimitsCascade[40] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 14, 17, 20};
    for (int i = 0; i<30; i++){
        if(i < 20)lPtbinlimitsCascade[i] = i/10.;
        if ((i >= 20)&&(i<30))lPtbinlimitsCascade[i] = 2.0 + (i-20)*0.2;
    }
    lPtbinlimitsCascade[30] = 4.5;
    lPtbinlimitsCascade[31] = 5.0;
    lPtbinlimitsCascade[32] = 5.5;
    lPtbinlimitsCascade[33] = 6.5;
    lPtbinlimitsCascade[34] = 8.0;
    lPtbinlimitsCascade[35] = 10;
    lPtbinlimitsCascade[36] = 12;
    lPtbinlimitsCascade[37] = 14;
    lPtbinlimitsCascade[38] = 17;
    lPtbinlimitsCascade[39] = 20;
    
    // Long_t lPtbinnumbCascade = sizeof(lPtbinlimitsCascade)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimits[] = {0, 10};
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleName[] = {"K0Short", "Lambda",  "AntiLambda"};
    
    //STEP 3: Creation of output objects
    
    //Map to mass hypothesis
    AliV0Result::EMassHypo lMassHypoV0[3];
    lMassHypoV0[0] = AliV0Result::kK0Short;
    lMassHypoV0[1] = AliV0Result::kLambda;
    lMassHypoV0[2] = AliV0Result::kAntiLambda;
    
    // Float_t lLifetimeCut[3];
    lLifetimeCut[0] = 20.0;
    lLifetimeCut[1] = 30.0;
    lLifetimeCut[2] = 30.0;
    
    // Float_t lMass[3];
    lMass[0] = 0.497;
    lMass[1] = 1.116;
    lMass[2] = 1.116;
    
    // Float_t lMWindow[3];
    lMWindow[0] = 0.075;
    lMWindow[1] = 0.050;
    lMWindow[2] = 0.050;
    
    //Array of results
    AliV0Result *lV0Result[5000];
    Long_t lNV0 = 0;
    
    //Central results: Stored in indices 0, 1, 2 (careful!)
    for(Int_t i = 0 ; i < 3 ; i ++){
        //Central result, customized binning: the one to use, usually
        lV0Result[lNV0] = new AliV0Result( Form("%s_Central",lParticleName[i].Data() ),lMassHypoV0[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits, 100,lMass[i]-lMWindow[i],lMass[i]+lMWindow[i]);
        //if ( i>0 ) not neeed for real data
        //    lV0Result[lNV0]->InitializeFeeddownMatrix( lPtbinnumb, lPtbinlimits, lPtbinnumbCascade, lPtbinlimitsCascade, lCentbinnumb, lCentbinlimits );
        
        //Setters for V0 Cuts
        lV0Result[lNV0]->SetCutDCANegToPV            ( 0.05 ) ;
        lV0Result[lNV0]->SetCutDCAPosToPV            ( 0.05 ) ;
        lV0Result[lNV0]->SetCutDCAV0Daughters        ( 1.2 ) ;
        lV0Result[lNV0]->SetCutV0CosPA               ( 0.98 ) ;
        lV0Result[lNV0]->SetCutV0Radius              ( 0.9 ) ;
        
        //Miscellaneous
        lV0Result[lNV0]->SetCutProperLifetime        ( lLifetimeCut[i] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( 70 ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( 0.8 ) ;
        lV0Result[lNV0]->SetCutTPCdEdx               ( 4 ) ;
        
        //Add result to pool
        lNV0++;
    }
    
    //Will now proceed to sweep individual variables
    //Number of steps used for the variable sweeps
    const Int_t lNumberOfSteps = lRecNumberOfSteps;
    
    //________________________________________________________
    // Variable 1: DCA Neg to PV
    Float_t lMaxDCANegToPV = 20.00;
    
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCANegToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCANegToPV / ((Float_t) lNumberOfSteps) ;
            lV0Result[lNV0] -> SetCutDCANegToPV ( lThisCut );
            lNV0++;
        }
    }
    //________________________________________________________
    // Variable 2: DCA Pos to PV
    Float_t lMaxDCAPosToPV = 20.00;
    
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAPosToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAPosToPV / ((Float_t) lNumberOfSteps) ;
            lV0Result[lNV0] -> SetCutDCAPosToPV ( lThisCut );
            lNV0++;
        }
    }
    //________________________________________________________
    // Variable 3: DCA V0 Daughters
    Float_t lMaxDCAV0Daughters = 1.20;
    
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAV0DaughtersSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAV0Daughters / ((Float_t) lNumberOfSteps) ;
            lV0Result[lNV0] -> SetCutDCAV0Daughters ( lThisCut );
            lNV0++;
        }
    }
    //________________________________________________________
    // Variable 4: V0 CosPA
    lMinV0CosPA = 0.98;
    lMaxV0CosPA = 1.00;
    Double_t lV0CosPAVals[lNumberOfSteps];
    // Double_t lMinV0PA = 1.0;
    //lMinV0PA = 0.0;
    Double_t lMaxV0PA = TMath::ACos(lMinV0CosPA);
    Double_t lDeltaV0PA = lMaxV0PA / ((Double_t)(lNumberOfSteps));
    for(Int_t iStep = 0; iStep<lNumberOfSteps; iStep++){
        lV0CosPAVals[iStep] = TMath::Cos( ((Float_t)(iStep+1))*lDeltaV0PA );
    }
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"V0CosPASweep",icut) );
            //Add result to pool
            lV0Result[lNV0] -> SetCutV0CosPA ( lV0CosPAVals[icut] );
            lNV0++;
        }
    }
    //________________________________________________________
    // Variable 5: V0 Radius
    Float_t lMinV0Radius = 2.0;
    Float_t lMaxV0Radius = 20.00;
    for(Int_t i = 0 ; i < 3 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%i",lParticleName[i].Data(),"V0RadiusSweep",icut) );
            //Add result to pool
            Float_t lThisCut = lMinV0Radius + (lMaxV0Radius-lMinV0Radius)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            lV0Result[lNV0] -> SetCutV0Radius ( lThisCut );
            lNV0++;
        }
    }
    for (Int_t iconf = 0; iconf<lNV0; iconf++)
        AddConfiguration(lV0Result[iconf]);
    
    //  cout<<"Added "<<lNV0<<" V0 configurations to output."<<endl;
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::AddTopologicalQACascade(Int_t lRecNumberOfSteps)
//Add all configurations to do QA of topological variables for the V0 analysis
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    //Double_t lPtbinlimits[] = {0.4, 0.5, 0.6,
    //    0.7,0.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
    //    2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,
    //    4.4,4.6,4.8,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.,11.,12.};
    //Double_t lPtbinlimits[] = {0.2,0.3, 0.4, 0.5, 0.6,
    //    0.7,0.8,.9,1.0,1.2, 1.4, 1.6, 1.8 ,2.0,
    //    2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
    //    4.4,4.8,5.0,6.0,7.0,8.0,9.0,10.,11.,12.};
    
    Double_t lPtbinlimits[] = {0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,3.2,4.0,5.0,8.};
    Double_t lPtbinlimitsOmega[] = {1.2,1.6,2.0,2.4,3.2,4.0,5.0,8.};
    
    //Double_t lPtbinlimits[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2,
    //3.6, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.5, 8.5, 10, 12};
    
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    Long_t lPtbinnumbOmega = sizeof(lPtbinlimitsOmega)/sizeof(Double_t) - 1;
    
    // centrality binning
    //Double_t lCentbinlimits[] = {0, 10}; //optimize in 0-10%
    Double_t lCentbinlimits[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    //Just a counter and one array, please. Nothing else needed
    AliCascadeResult *lCascadeResult[10000];
    Long_t lN = 0;
    
    //Map to mass hypothesis
    AliCascadeResult::EMassHypo lMassHypo[4];
    lMassHypo[0] = AliCascadeResult::kXiMinus;
    lMassHypo[1] = AliCascadeResult::kXiPlus;
    lMassHypo[2] = AliCascadeResult::kOmegaMinus;
    lMassHypo[3] = AliCascadeResult::kOmegaPlus;
    
    //lLifetimeCut[4];
    lLifetimeCut[0] = 15.0;
    lLifetimeCut[1] = 15.0;
    lLifetimeCut[2] = 12.0;
    lLifetimeCut[3] = 12.0;
    
    //Float_t lMass[4];
    lMass[0] = 1.322;
    lMass[1] = 1.322;
    lMass[2] = 1.672;
    lMass[3] = 1.672;
    
    TString lParticleName[] = {"XiMinus", "XiPlus",  "OmegaMinus", "OmegaPlus"};
    
    //Number of steps used for the variable sweeps
    const Int_t lNumberOfSteps = lRecNumberOfSteps;
    
    //Central results: Stored in indices 0, 1, 2, 3 (careful!)
    for(Int_t i = 0 ; i < 4 ; i ++){
        //Central result, customized binning: the one to use, usually
        if( i<2 ){
            lCascadeResult[lN] = new AliCascadeResult( Form("%s_AnalysisLevel",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits,100,lMass[i]-0.050,lMass[i]+0.050);
        }else{
            lCascadeResult[lN] = new AliCascadeResult( Form("%s_AnalysisLevel",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumbOmega,lPtbinlimitsOmega,100,lMass[i]-0.050,lMass[i]+0.050);
        }
        
        //Default cuts: use vertexer level ones
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( 0.2 ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( 0.2 ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        (  1. ) ;
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarV0CosPA            (TMath::Exp(10.853),
                                                         -25.0322,
                                                         TMath::Exp(-0.843948),
                                                         -0.890794,
                                                         0.057553);
        
        lCascadeResult[lN]->SetCutV0Radius              (  3 ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( 0.1 ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( 0.005 ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( 0.1 ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( 1.0) ;
        lCascadeResult[lN]->SetCutCascRadius            ( 1.2 ) ;
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.95 ) ; //+variable
        if(i < 2){
            lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(4.86664),
                                                             -10.786,
                                                             TMath::Exp(-1.33411),
                                                             -0.729825,
                                                             0.0695724);
        }
        if(i >= 2){
            lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(   12.8752),
                                                             -21.522,
                                                             TMath::Exp( -1.49906),
                                                             -0.813472,
                                                             0.0480962);
        }
        lCascadeResult[lN]->SetCutDCACascadeToPV        ( 0.8 );
        if(i >= 2){
            lCascadeResult[lN]->SetCutDCACascDaughters  ( 0.6 ) ;
            lCascadeResult[lN]->SetCutCascRadius        ( 1.0 ) ;
            lCascadeResult[lN]->SetCutDCACascadeToPV    ( 0.6 ) ;
        }
        lCascadeResult[lN]->SetCutVarDCACascDau ( TMath::Exp(0.0470076), -0.917006, 0, 1, 0.5 );
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( lLifetimeCut[i] ) ;
        lCascadeResult[lN]->SetCutMaxV0Lifetime         ( 30.0  );
        lCascadeResult[lN]->SetCutMinTrackLength        ( 90.0  );
        lCascadeResult[lN]->SetCutLeastNumberOfClusters( -1 );
        lCascadeResult[lN]->SetCutTPCdEdx               ( 3.0 ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( 0.008 ) ;
        lCascadeResult[lN]->SetCutBachBaryonCosPA       ( TMath::Cos(0.04) ) ; //+variable
        lCascadeResult[lN]->SetCutVarBBCosPA            (TMath::Exp(-2.29048),
                                                         -20.2016,
                                                         TMath::Exp(-2.9581),
                                                         -0.649153,
                                                         0.00526455);
        //Add result to pool
        lN++;
    }
    
    //Will now proceed to sweep individual variables
    
    //________________________________________________________
    // Variable 1: DCA Neg to PV
    Float_t lMaxDCANegToPV = 1.5;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCANegToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCANegToPV / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCANegToPV ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 2: DCA Pos to PV
    Float_t lMaxDCAPosToPV = 1.5;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAPosToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAPosToPV / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCAPosToPV ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 3: DCA V0 Daughters
    Float_t lMaxDCAV0Daughters = 1.40;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAV0DaughtersSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAV0Daughters / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCAV0Daughters ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 4: V0 CosPA
    lMinV0CosPA = 0.95;
    // Float_t lMaxV0CosPA = 1.01;
    lMaxV0CosPA = 1.00;
    Double_t lV0CosPAVals[lNumberOfSteps];
    // Double_t lMinV0PA = 0.1;
    // lMinV0PA = 0.0;
    Double_t lMaxV0PA = TMath::ACos(lMinV0CosPA);
    Double_t lDeltaV0PA = lMaxV0PA / ((Double_t)(lNumberOfSteps));
    for(Int_t iStep = 0; iStep<lNumberOfSteps; iStep++){
        lV0CosPAVals[iStep] = TMath::Cos( ((Float_t)(iStep+1))*lDeltaV0PA );
    }
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"V0CosPASweep",icut) );
            //Add result to pool
            lCascadeResult[lN] -> SetCutUseVarV0CosPA( kFALSE );
            lCascadeResult[lN] -> SetCutV0CosPA ( lV0CosPAVals[icut] );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 5: V0 Radius
    Float_t lMinV0Radius = 0.0;
    Float_t lMaxV0Radius = 20.00;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"V0RadiusSweep",icut) );
            //Add result to pool
            Float_t lThisCut = lMinV0Radius + (lMaxV0Radius-lMinV0Radius)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            lCascadeResult[lN] -> SetCutV0Radius ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 6: DCA V0 To PV
    Float_t lMaxDCAV0ToPV = 0.5;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCAV0ToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCAV0ToPV / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCAV0ToPV ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 7: DCA Bach To PV
    Float_t lMaxDCABachToPV = 0.5;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCABachToPVSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCABachToPV / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutDCABachToPV ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 8: DCA Casc Daughters
    Float_t lMaxDCACascDaughters = 1.40;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCACascDaughtersSweep",icut) );
            //Add result to pool
            Float_t lThisCut = ((Float_t)icut+1)*lMaxDCACascDaughters / ((Float_t) lNumberOfSteps) ;
            lCascadeResult[lN] -> SetCutUseVarDCACascDau ( kFALSE );
            lCascadeResult[lN] -> SetCutDCACascDaughters ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 9: Cascade Radius
    Float_t lMinCascRadius = 0.5;
    Float_t lMaxCascRadius = 7.0;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"CascRadiusSweep",icut) );
            //Add result to pool
            Float_t lThisCut = lMinCascRadius + (lMaxCascRadius-lMinCascRadius)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            lCascadeResult[lN] -> SetCutCascRadius ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 10: Cascade CosPA
    Float_t lMinCascCosPA = 0.95;
    // lMaxCascCosPA = 1.00;
    Double_t lCascCosPAVals[lNumberOfSteps];
    //lMinCascPA = 0.0;
    Double_t lMaxCascPA = TMath::ACos(lMinCascCosPA);
    Double_t lDeltaCascPA = lMaxCascPA / ((Double_t)(lNumberOfSteps));
    for(Int_t iStep = 0; iStep<lNumberOfSteps; iStep++){
        lCascCosPAVals[iStep] = TMath::Cos( ((Float_t)(iStep+1))*lDeltaCascPA );
    }
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"CascCosPASweep",icut) );
            //Add result to pool
            lCascadeResult[lN] -> SetCutUseVarCascCosPA( kFALSE );
            lCascadeResult[lN] -> SetCutCascCosPA ( lCascCosPAVals[icut] );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 11: Bach-Baryon CosPA
    Float_t lMinBBCosPA = TMath::Cos(0.1);
    //lMaxBBCosPA = 1.000;
    Double_t lBBCosPAVals[lNumberOfSteps];
    //Double_t lMinBBPA = 0.0;
    Double_t lMaxBBPA = TMath::ACos(lMinBBCosPA);
    Double_t lDeltaBBPA = lMaxBBPA / ((Double_t)(lNumberOfSteps));
    for(Int_t iStep = 0; iStep<lNumberOfSteps; iStep++){
        lBBCosPAVals[iStep] = TMath::Cos( ((Float_t)(iStep+1))*lDeltaBBPA );
    }
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"BBCosPASweep",icut) );
            //Add result to pool
            lCascadeResult[lN] -> SetCutUseVarBBCosPA( kFALSE );
            lCascadeResult[lN] -> SetCutBachBaryonCosPA ( lBBCosPAVals[icut] );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 12: Cascade Lifetime Sweep
    
    Int_t lLifetimeSteps = 15;
    for(Int_t i = 0 ; i < 4 ; i ++){
        Float_t lMinLifetime = 5.00;
        Float_t lMaxLifetime = 20.00;
        for(Int_t icut = 0; icut<lLifetimeSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"CascLifetimeSweep",icut) );
            Float_t lThisCut = lMinLifetime + (lMaxLifetime-lMinLifetime)*(((Float_t)icut)+1)/((Float_t)lLifetimeSteps);
            //Add result to pool
            lCascadeResult[lN] -> SetCutProperLifetime ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 13: V0 Lifetime Sweep
    Float_t lMinV0Lifetime = 8.00;
    Float_t lMaxV0Lifetime = 40.00;
    Int_t lV0LifetimeSteps = 32;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lV0LifetimeSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"MaxV0LifetimeSweep",icut) );
            Float_t lThisCut = lMinV0Lifetime + (lMaxV0Lifetime-lMinV0Lifetime)*(((Float_t)icut)+1)/((Float_t)lV0LifetimeSteps);
            //Add result to pool
            lCascadeResult[lN] -> SetCutMaxV0Lifetime ( lThisCut );
            lN++;
        }
    }
    //________________________________________________________
    // Variable 14: Xi Rejection Sweep (Omega only)
    Float_t lMaxXiRejection = 0.015;
    for(Int_t i = 2 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"XiRejectionSweep",icut) );
            Float_t lThisCut = ((Float_t)icut+1)*lMaxXiRejection / ((Float_t) lNumberOfSteps) ;
            //Add result to pool
            lCascadeResult[lN] -> SetCutXiRejection ( lThisCut );
            lN++;
        }
    }
    
    //________________________________________________________
    // Variable 15: V0 Mass Window
    Float_t lMinV0Mass = 0.002;
    Float_t lMaxV0Mass = 0.010;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"V0MassWindow",icut) );
            Float_t lThisCut = lMinV0Mass + (lMaxV0Mass-lMinV0Mass)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            //Add result to pool
            lCascadeResult[lN] -> SetCutV0Mass ( lThisCut );
            lN++;
        }
    }
    
    //________________________________________________________
    // Variable 16: DCA Cascade To PV
    Float_t lMinDCACascToPV = 0.2;
    Float_t lMaxDCACascToPV = 1.2;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"DCACascadeToPV",icut) );
            Float_t lThisCut = lMinDCACascToPV + (lMaxDCACascToPV-lMinDCACascToPV)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            //Add result to pool
            lCascadeResult[lN] -> SetCutDCACascadeToPV ( lThisCut );
            lN++;
        }
    }
    
    //________________________________________________________
    // Variable 17: Minimum Track Length
    Float_t lMinTrackLength = 60;
    Float_t lMaxTrackLength = 120;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"MinimumTrackLength",icut) );
            Float_t lThisCut = lMinTrackLength + (lMaxTrackLength-lMinTrackLength)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            //Add result to pool
            lCascadeResult[lN] -> SetCutMinTrackLength ( lThisCut );
            lN++;
        }
    }
    
    //________________________________________________________
    // Variable 18: TPC dE/dx
    Float_t lMinTPCnSigma = 2.5;
    Float_t lMaxTPCnSigma = 5.0;
    for(Int_t i = 0 ; i < 4 ; i ++){
        for(Int_t icut = 0; icut<lNumberOfSteps; icut++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%i",lParticleName[i].Data(),"TPCnSigma",icut) );
            Float_t lThisCut = lMinTPCnSigma + (lMaxTPCnSigma-lMinTPCnSigma)*(((Float_t)icut)+1)/((Float_t)lNumberOfSteps);
            //Add result to pool
            lCascadeResult[lN] -> SetCutTPCdEdx ( lThisCut );
            lN++;
        }
    }
    
    for (Int_t iconf = 0; iconf<lN; iconf++)
        AddConfiguration(lCascadeResult[iconf]);
    
    // cout<<"Added "<<lN<<" Cascade configurations to output."<<endl;
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::AddStandardV0Configuration(Bool_t lUseFull, Bool_t lDoSweepLooseTight, Int_t lSweepFullNumb)
//Meant to add some standard V0 analysis Configuration + its corresponding systematics
{
    //======================================================
    // V0 Configurations To use
    //======================================================
    
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimitsV0[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 14, 15, 17, 20};
    Long_t lPtbinnumbV0 = sizeof(lPtbinlimitsV0)/sizeof(Double_t) - 1;
    Double_t lPtbinlimitsXi[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 14, 16, 19, 22, 25};
    Long_t lPtbinnumbXi = sizeof(lPtbinlimitsXi)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimitsV0[] = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90};
    Long_t lCentbinnumbV0 = sizeof(lCentbinlimitsV0)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleNameV0[] =
    {
        "K0Short",
        "Lambda",
        "AntiLambda"
    };
    const Int_t lNPart = 3;
    TString lConfNameV0[] =
    {
        "Loose",
        "Central",
        "Tight"
    };
    const Int_t lNConf = 3;
    TString lCutNameV0[] =
    {
        "DCANegToPV",
        "DCAPosToPV",
        "DCAV0Daughters",
        "V0CosPA",
        "V0Radius",
        "ProperLifetime",
        "TrackLength",
        "LeastNbrCrsOvFind",
        "TPCdEdx",
        "APParameter",
        "V0RadiusMax",
        "LeastNbrCrsRows"
    };
    const Int_t lNCutsForSyst = 10;
    
    // STEP 2: Decide on a set of selections
    
    //1st index: Particle Species
    //2nd index: Loose / Central / Tight
    //3rd index: Number of selection (as ordered above)
    Double_t lcutsV0[lNPart][lNConf][lNCutsForSyst];
    
    //1st index: Particle Species: K0Short, Lambda, AntiLambda
    //2nd index: Loose / Central / Tight: 2%, 5%, 10% signal loss
    Double_t parExp0Const[lNPart][lNConf];
    Double_t parExp0Slope[lNPart][lNConf];
    Double_t parExp1Const[lNPart][lNConf];
    Double_t parExp1Slope[lNPart][lNConf];
    Double_t parConst[lNPart][lNConf];
    
    //=============================================================================================
    // K0SHORT V0 COS PA PARAMETRIZATION
    //---------------------------------------------------------------------------------------------
    //                       LOOSE                         CENTRAL                           TIGHT
    parExp0Const[0][0] =  0.20428;  parExp0Const[0][1] =  0.22692;  parExp0Const[0][2] =  0.28814;
    parExp0Slope[0][0] = -0.73728;  parExp0Slope[0][1] = -1.59317;  parExp0Slope[0][2] = -2.27069;
    parExp1Const[0][0] =  0.09887;  parExp1Const[0][1] =  0.05994;  parExp1Const[0][2] =  0.04320;
    parExp1Slope[0][0] = -0.02822;  parExp1Slope[0][1] = -0.26997;  parExp1Slope[0][2] = -0.29839;
    parConst[0][0] = -0.05302;      parConst[0][1] =  0.00907;      parConst[0][2] =  0.00704;
    //=============================================================================================
    
    //=============================================================================================
    // LAMBDA V0 COS PA PARAMETRIZATION
    //---------------------------------------------------------------------------------------------
    
    
    //                       LOOSE                         CENTRAL                           TIGHT
    parExp0Const[1][0] =  0.22775;  parExp0Const[1][1] =  0.36284;  parExp0Const[1][2] =  0.54877;
    parExp0Slope[1][0] = -1.11579;  parExp0Slope[1][1] = -1.87960;  parExp0Slope[1][2] = -2.72912;
    parExp1Const[1][0] =  0.06266;  parExp1Const[1][1] =  0.04543;  parExp1Const[1][2] =  0.03411;
    parExp1Slope[1][0] = -0.17086;  parExp1Slope[1][1] = -0.20447;  parExp1Slope[1][2] = -0.26965;
    parConst[1][0] =  0.01489;      parConst[1][1] =  0.01085;      parConst[1][2] =  0.00889;
    //=============================================================================================
    
    //=============================================================================================
    // ANTILAMBDA V0 COS PA PARAMETRIZATION
    //---------------------------------------------------------------------------------------------
    //                       LOOSE                         CENTRAL                           TIGHT
    parExp0Const[2][0] =  0.22667;  parExp0Const[2][1] =  0.35809;  parExp0Const[2][2] =  0.54114;
    parExp0Slope[2][0] = -0.93618;  parExp0Slope[2][1] = -1.93860;  parExp0Slope[2][2] = -2.71000;
    parExp1Const[2][0] =  0.06857;  parExp1Const[2][1] =  0.05306;  parExp1Const[2][2] =  0.03664;
    parExp1Slope[2][0] = -0.07015;  parExp1Slope[2][1] = -0.24518;  parExp1Slope[2][2] = -0.28124;
    parConst[2][0] = -0.00707;      parConst[2][1] =  0.01213;      parConst[2][2] =  0.00905;
    //=============================================================================================
    
    //================================================================================
    // K0SHORT SELECTIONS
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcutsV0[0][0][ 0] = 0.05;    lcutsV0[0][1][ 0] =  0.10; lcutsV0[0][2][0] = 0.17; //DCANegToPV
    lcutsV0[0][0][ 1] = 0.05;    lcutsV0[0][1][ 1] =  0.10; lcutsV0[0][2][1] = 0.17; //DCAPosToPV
    lcutsV0[0][0][ 2] = 0.95;    lcutsV0[0][1][ 2] =   0.8; lcutsV0[0][2][2] =  0.7; //DCAV0Daughters
    lcutsV0[0][0][ 3] = 0.95;    lcutsV0[0][1][ 3] =  0.95; lcutsV0[0][2][3] = 0.95; //V0CosPA
    lcutsV0[0][0][ 4] = 4.50;    lcutsV0[0][1][ 4] =  5.00; lcutsV0[0][2][4] = 5.50; //V0Radius
    lcutsV0[0][0][ 5] =   25;    lcutsV0[0][1][ 5] =    20; lcutsV0[0][2][5] =   15; //Proper Lifetime (in cm)
    lcutsV0[0][0][ 6] =   80;    lcutsV0[0][1][ 6] =    90; lcutsV0[0][2][6] =  100; //Track Length
    lcutsV0[0][0][ 7] =  0.7;    lcutsV0[0][1][ 7] =   0.8; lcutsV0[0][2][7] = 0.85; //Least Ratio CrdRows/Findable
    lcutsV0[0][0][ 8] =  4.0;    lcutsV0[0][1][ 8] =   3.0; lcutsV0[0][2][8] =  2.5; //TPC dE/dx
    lcutsV0[0][0][ 9] = 0.18;    lcutsV0[0][1][ 9] =  0.20; lcutsV0[0][2][9] = 0.22; //AP Parameter
    //================================================================================
    
    //================================================================================
    // LAMBDA SELECTIONS
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcutsV0[1][0][ 0] = 0.10;    lcutsV0[1][1][ 0] =  0.25; lcutsV0[1][2][0] = 0.40; //DCANegToPV
    lcutsV0[1][0][ 1] = 0.08;    lcutsV0[1][1][ 1] =  0.10; lcutsV0[1][2][1] = 0.13; //DCAPosToPV
    lcutsV0[1][0][ 2] =  1.0;    lcutsV0[1][1][ 2] =   0.8; lcutsV0[1][2][2] = 0.65; //DCAV0Daughters
    lcutsV0[1][0][ 3] = 0.97;    lcutsV0[1][1][ 3] =  0.98; lcutsV0[1][2][3] = 0.99; //V0CosPA
    lcutsV0[1][0][ 4] = 4.00;    lcutsV0[1][1][ 4] =  5.00; lcutsV0[1][2][4] = 6.00; //V0Radius
    lcutsV0[1][0][ 5] =   30;    lcutsV0[1][1][ 5] =    25; lcutsV0[1][2][5] =   20; //Proper Lifetime (in cm)
    lcutsV0[1][0][ 6] =   80;    lcutsV0[1][1][ 6] =    90; lcutsV0[1][2][6] =  100; //Track Length
    lcutsV0[1][0][ 7] =  0.7;    lcutsV0[1][1][ 7] =   0.8; lcutsV0[1][2][7] = 0.85; //Least Ratio CrdRows/Findable
    lcutsV0[1][0][ 8] =  4.0;    lcutsV0[1][1][ 8] =   3.0; lcutsV0[1][2][8] =  2.5; //TPC dE/dx
    lcutsV0[1][0][ 9] = 0.18;    lcutsV0[1][1][ 9] =  0.20; lcutsV0[1][2][9] = 0.22; //AP Parameter
    //================================================================================
    
    //================================================================================
    // ANTILAMBDA SELECTIONS
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcutsV0[2][0][ 0] = 0.08;    lcutsV0[2][1][ 0] =  0.10; lcutsV0[2][2][0] = 0.13; //DCANegToPV
    lcutsV0[2][0][ 1] = 0.10;    lcutsV0[2][1][ 1] =  0.25; lcutsV0[2][2][1] = 0.40; //DCAPosToPV
    lcutsV0[2][0][ 2] =  1.0;    lcutsV0[2][1][ 2] =   0.8; lcutsV0[2][2][2] = 0.65; //DCAV0Daughters
    lcutsV0[2][0][ 3] = 0.97;    lcutsV0[2][1][ 3] =  0.98; lcutsV0[2][2][3] = 0.99; //V0CosPA
    lcutsV0[2][0][ 4] = 4.00;    lcutsV0[2][1][ 4] =  5.00; lcutsV0[2][2][4] = 6.00; //V0Radius
    lcutsV0[2][0][ 5] =   30;    lcutsV0[2][1][ 5] =    25; lcutsV0[2][2][5] =   20; //Proper Lifetime (in cm)
    lcutsV0[2][0][ 6] =   80;    lcutsV0[2][1][ 6] =    90; lcutsV0[2][2][6] =  100; //Track Length
    lcutsV0[2][0][ 7] =  0.7;    lcutsV0[2][1][ 7] =   0.8; lcutsV0[2][2][7] = 0.85; //Least Ratio CrdRows/Findable
    lcutsV0[2][0][ 8] =  4.0;    lcutsV0[2][1][ 8] =   3.0; lcutsV0[2][2][8] =  2.5; //TPC dE/dx
    lcutsV0[2][0][ 9] = 0.18;    lcutsV0[2][1][ 9] =  0.20; lcutsV0[2][2][9] = 0.22; //AP Parameter
    //================================================================================
    
    //STEP 3: Creation of output objects
    
    //Map to mass hypothesis
    AliV0Result::EMassHypo lMassHypoV0[lNPart];
    lMassHypoV0[0] = AliV0Result::kK0Short;
    lMassHypoV0[1] = AliV0Result::kLambda;
    lMassHypoV0[2] = AliV0Result::kAntiLambda;
    
    //Array of results
    AliV0Result *lV0Result[1000];
    Long_t lNV0 = 0;
    
    //Central results: Stored in indices 0, 1, 2 (careful!)
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        //Central result, customized binning: the one to use, usually
        lV0Result[lNV0] = new AliV0Result( Form("%s_Central",lParticleNameV0[i].Data() ),lMassHypoV0[i],"",lCentbinnumbV0,lCentbinlimitsV0, lPtbinnumbV0,lPtbinlimitsV0);
        
        //feeddown matrix
        if ( i!=0 ) lV0Result[lNV0] -> SetupFeeddownMatrix(lPtbinnumbXi,lPtbinlimitsXi);
        
        //Setters for V0 Cuts
        lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][1][ 0] ) ;
        lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][1][ 1] ) ;
        lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][1][ 2] ) ;
        lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][1][ 3] ) ;
        //Set Variable cut
        lV0Result[lNV0]->SetCutVarV0CosPA               ( parExp0Const[i][1], parExp0Slope[i][1], parExp1Const[i][1], parExp1Slope[i][1], parConst[i][1] ) ;
        
        lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][1][ 4] ) ;
        
        //Miscellaneous
        lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][1][ 5] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( -1 ) ; //no cut here
        lV0Result[lNV0]->SetCutMinTrackLength ( lcutsV0[i][1][ 6] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][1][ 7] ) ;
        lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][1][ 8] ) ;
        lV0Result[lNV0]->SetCutArmenterosParameter               ( lcutsV0[i][1][ 9] ) ;
        
        //Add result to pool
        lNV0++;
    }
    
    //Central full results
    if (lUseFull){
        for(Int_t i = 0 ; i < lNPart ; i ++)
        {
            //Central Result, Full: No rebinning. Will use a significant amount of memory,
            //not to be replicated several times for systematics!
            lV0Result[lNV0] = new AliV0Result( Form("%s_Central_Full",lParticleNameV0[i].Data() ),lMassHypoV0[i]);
            
            //feeddown matrix
            if ( i!=0 ) lV0Result[lNV0] -> SetupFeeddownMatrix(lPtbinnumbXi,lPtbinlimitsXi);
            
            //Setters for V0 Cuts
            lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][1][ 0] ) ;
            lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][1][ 1] ) ;
            lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][1][ 2] ) ;
            lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][1][ 3] ) ;
            //Set Variable cut
            lV0Result[lNV0]->SetCutVarV0CosPA               ( parExp0Const[i][1], parExp0Slope[i][1], parExp1Const[i][1], parExp1Slope[i][1], parConst[i][1] ) ;
            
            lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][1][ 4] ) ;
            
            //Miscellaneous
            lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][1][ 5] ) ;
            lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( -1 ) ; //no cut here
            lV0Result[lNV0]->SetCutMinTrackLength ( lcutsV0[i][1][ 6] ) ;
            lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][1][ 7] ) ;
            lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][1][ 8] ) ;
            
            //Add result to pool
            lNV0++;
        }
    }
    
    //Rapidity sweep
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        for(Int_t ir=0; ir<12; ir++)
        {
            Float_t lLoRap = -0.6 +   (ir)*0.1;
            Float_t lHiRap = -0.6 + (ir+1)*0.1;
            
            //Create a new object from default
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_RapiditySweep_%.1f_%.1f",lParticleNameV0[i].Data(),lLoRap, lHiRap) );
            lV0Result[lNV0] -> SetCutMinRapidity( lLoRap );
            lV0Result[lNV0] -> SetCutMaxRapidity( lHiRap );
            
            //Add result to pool
            lNV0++;
        }
    }
    
    //Number of crossed rows cut
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s",lParticleNameV0[i].Data(),"NCrossedRowsCut") );
        
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( 70 ) ;
        lV0Result[lNV0]->SetCutMinTrackLength ( -1 ) ;
        
        //Add result to pool
        lNV0++;
    }
    
    //No Armenteros-Podolanski cut
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s",lParticleNameV0[i].Data(),"NoAP") );
        lV0Result[lNV0]->SetCutArmenterosParameter(0.0);
        
        //Add result to pool
        lNV0++;
    }
    
    //Explore Use MC Properties vs Use Reconstructed Properties
    for(Int_t i = 0 ; i < lNPart ; i ++){
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_MCUseRecoProp",lParticleNameV0[i].Data() ) );
        
        lV0Result[lNV0]->SetCutMCUseMCProperties(kFALSE);
        
        //Add result to pool
        lNV0++;
    }
    
    //Explore TOF information use
    for(Int_t i = 0 ; i < lNPart ; i ++){
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_AtLeastOneTOF",lParticleNameV0[i].Data() ) );
        lV0Result[lNV0]->SetCutAtLeastOneTOF(kTRUE);
        
        //Add result to pool
        lNV0++;
    }
    
    //Explore no TPC dedx
    for(Int_t i = 0 ; i < lNPart ; i ++){
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_NoTPCdEdx",lParticleNameV0[i].Data() ) );
        lV0Result[lNV0]->SetCutTPCdEdx(1e+6);
        
        //Add result to pool
        lNV0++;
    }
    
    //Explore ITS refit requirement
    for(Int_t i = 0 ; i < lNPart ; i ++){
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_ITSRefitTracks",lParticleNameV0[i].Data() ) );
        lV0Result[lNV0]->SetCutUseITSRefitTracks(kTRUE);
        
        //Add result to pool
        lNV0++;
    }
    
    //================================================================================
    //Set up cut values, tight and loose versions
    //--------------------------------------------------------------------------------
    const Int_t lNCutsForSweep = 12;
    //1st index: Particle Species
    //2nd index: Number of selection
    Double_t lCutsTight[lNPart][lNCutsForSweep];
    Double_t lCutsLoose[lNPart][lNCutsForSweep];
    
    Double_t lMeanLifetime[] = {2.6844, 7.89, 7.89};
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        lCutsTight[i][ 0] =   0.1; //DCANegToPV
        lCutsTight[i][ 1] =   0.1; //DCAPosToPV
        lCutsTight[i][ 2] =     1; //DCAV0Daughters
        lCutsTight[i][ 3] = 0.998; //V0CosPA
        lCutsTight[i][ 4] =     5; //V0Radius
        lCutsTight[i][ 5] = 3 * lMeanLifetime[i]; //Proper Lifetime (in cm)
        lCutsTight[i][ 6] =    -1; //Track Length
        lCutsTight[i][ 7] = -0.01; //Least Ratio CrdRows/Findable
        lCutsTight[i][ 8] =     8; //TPC dE/dx
        lCutsTight[i][ 9] =   0.2; //AP Parameter
        lCutsTight[i][10] =   100; //V0Radius max
        lCutsTight[i][11] =    70; //Least number of CrdRows
        
        for(Int_t j = 0; j < lNCutsForSyst; j++)
        {
            lCutsLoose[i][j] = lcutsV0[i][1][j];
        }
        lCutsLoose[i][10] =   200; //V0Radius max
        lCutsLoose[i][11] =    -1; //Least number of CrdRows
    }
    //================================================================================
    
    //2.76 TeV analysis cuts
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_276Cuts",lParticleNameV0[i].Data() ) );
        
        //Setters for V0 Cuts
        lV0Result[lNV0]->SetCutDCANegToPV            ( lCutsTight[i][0] ) ;
        lV0Result[lNV0]->SetCutDCAPosToPV            ( lCutsTight[i][1] ) ;
        lV0Result[lNV0]->SetCutDCAV0Daughters        ( lCutsTight[i][2] ) ;
        lV0Result[lNV0]->SetCutV0CosPA               ( lCutsTight[i][3] ) ;
        //Use constant cos(PA) cut
        lV0Result[lNV0]->SetCutUseVarV0CosPA         ( kFALSE ) ;
        
        lV0Result[lNV0]->SetCutV0Radius              ( lCutsTight[i][4] ) ;
        lV0Result[lNV0]->SetCutMaxV0Radius           ( lCutsTight[i][10] ) ;
        
        //Miscellaneous
        lV0Result[lNV0]->SetCutProperLifetime        ( lCutsTight[i][5] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( lCutsTight[i][11] ) ;
        lV0Result[lNV0]->SetCutMinTrackLength ( lCutsTight[i][6] ) ; //no cut here
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lCutsTight[i][7] ) ;
        lV0Result[lNV0]->SetCutTPCdEdx               ( 1e6 ) ; //no cut here
        lV0Result[lNV0]->SetCut276TeVLikedEdx        ( kTRUE ) ;
        lV0Result[lNV0]->SetCutArmenterosParameter               ( lCutsTight[i][9] ) ;
        
        //Add result to pool
        lNV0++;
    }
    
    //================================================================================
    // Topo sweeps for tests
    if ( lDoSweepLooseTight )
    { // begin sweeps between loose and tight versions of cuts
        // centrality binning for sweeps
        Double_t lSweepCentBinLimits[] = {0, 90};
        Long_t lSweepCentBinNumb = sizeof(lSweepCentBinLimits)/sizeof(Double_t) - 1;
        
        //Mass binning for sweeps
        Double_t lNMassBins[] = {400, 400, 400};
        Double_t lMass[] = {0.498, 1.116, 1.116};
        Double_t lMassWindow[] = {0.15,0.1,0.1};
        
        //Loose cuts for sweeps
        Int_t lLooseForSweepIndex = lNV0;
        for(Int_t i = 0 ; i < lNPart ; i++)
        {
            //Central result for sweeps
            lV0Result[lNV0] = new AliV0Result( Form("%s_Central_ForSweep",lParticleNameV0[i].Data() ),lMassHypoV0[i],"",lSweepCentBinNumb,lSweepCentBinLimits, lPtbinnumbV0,lPtbinlimitsV0,lNMassBins[i],lMass[i]-lMassWindow[i],lMass[i]+lMassWindow[i]);
            
            //feeddown matrix
            if ( i!=0 ) lV0Result[lNV0] -> SetupFeeddownMatrix(lPtbinnumbXi,lPtbinlimitsXi);
            
            //Setters for V0 Cuts
            lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][1][ 0] ) ;
            lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][1][ 1] ) ;
            lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][1][ 2] ) ;
            lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][1][ 3] ) ;
            //Set Variable cut
            lV0Result[lNV0]->SetCutVarV0CosPA               ( parExp0Const[i][1], parExp0Slope[i][1], parExp1Const[i][1], parExp1Slope[i][1], parConst[i][1] ) ;
            
            lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][1][ 4] ) ;
            
            //Miscellaneous
            lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][1][ 5] ) ;
            lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( -1 ) ; //no cut here
            lV0Result[lNV0]->SetCutMinTrackLength ( lcutsV0[i][1][ 6] ) ;
            lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][1][ 7] ) ;
            lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][1][ 8] ) ;
            lV0Result[lNV0]->SetCutArmenterosParameter               ( lcutsV0[i][1][ 9] ) ;
            
            //Add result to pool
            lNV0++;
        }
        
        //Tightening cuts one by one
        for(Int_t i = 0 ; i < lNPart ; i++)
        {
            for(Int_t iCut = 0; iCut < lNCutsForSweep; iCut++)
            {
                //only proceed if cuts are actually different
                if( ( TMath::Abs( lCutsTight[i][iCut]-lCutsLoose[i][iCut] ) / lCutsLoose[i][iCut] < 0.01 ) && ( iCut != 3 ) ) continue;
                
                Int_t lNSweep = 12;
                for(Int_t iSweep = 1; iSweep <= lNSweep; iSweep++)
                {
                    Double_t lCutValue = lCutsLoose[i][iCut] + (iSweep/(Double_t)lNSweep)*(lCutsTight[i][iCut]-lCutsLoose[i][iCut]);
                    
                    lV0Result[lNV0] = new AliV0Result( lV0Result[lLooseForSweepIndex+i], Form( "%s_Central_%s_%d", lParticleNameV0[i].Data(), lCutNameV0[iCut].Data(), iSweep ) );
                    
                    if(iCut ==  0 ) lV0Result[lNV0]->SetCutDCANegToPV            ( lCutValue ) ;
                    if(iCut ==  1 ) lV0Result[lNV0]->SetCutDCAPosToPV            ( lCutValue ) ;
                    if(iCut ==  2 ) lV0Result[lNV0]->SetCutDCAV0Daughters        ( lCutValue ) ;
                    if(iCut ==  3 )
                    {
                        lV0Result[lNV0]->SetCutV0CosPA               ( lCutValue ) ;
                        lV0Result[lNV0]->SetCutVarV0CosPA            ( parExp0Const[i][1]*(1-iSweep/(Double_t)lNSweep), parExp0Slope[i][1], parExp1Const[i][1]*(1-iSweep/(Double_t)lNSweep), parExp1Slope[i][1], parConst[i][1] + (iSweep/(Double_t)lNSweep)*(TMath::ACos(lCutsTight[i][iCut])-parConst[i][1]) ) ;
                    }
                    if(iCut ==  4 ) lV0Result[lNV0]->SetCutV0Radius              ( lCutValue ) ;
                    
                    //Miscellaneous
                    if(iCut ==  5 ) lV0Result[lNV0]->SetCutProperLifetime        ( lCutValue ) ;
                    if(iCut ==  6 ) lV0Result[lNV0]->SetCutMinTrackLength ( lCutValue ) ;
                    if(iCut ==  7 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lCutValue ) ;
                    if(iCut ==  8 ) lV0Result[lNV0]->SetCutTPCdEdx               ( lCutValue ) ;
                    if(iCut ==  9 ) lV0Result[lNV0]->SetCutArmenterosParameter   ( lCutValue ) ;
                    if(iCut == 10 ) lV0Result[lNV0]->SetCutMaxV0Radius           ( lCutValue ) ;
                    if(iCut == 11 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( lCutValue ) ;
                    
                    //Print this variation, add to pool
                    lV0Result[lNV0]->Print();
                    lNV0++;
                }
            }
        }
        
        //Tight cuts for sweeps
        Int_t lTightForSweepIndex = lNV0;
        for(Int_t i = 0 ; i < lNPart ; i++)
        {
            lV0Result[lNV0] = new AliV0Result( lV0Result[lLooseForSweepIndex+i], Form("%s_276Cuts_ForSweep",lParticleNameV0[i].Data() ) );
            
            //Setters for V0 Cuts
            lV0Result[lNV0]->SetCutDCANegToPV            ( lCutsTight[i][0] ) ;
            lV0Result[lNV0]->SetCutDCAPosToPV            ( lCutsTight[i][1] ) ;
            lV0Result[lNV0]->SetCutDCAV0Daughters        ( lCutsTight[i][2] ) ;
            lV0Result[lNV0]->SetCutV0CosPA               ( lCutsTight[i][3] ) ;
            //Use constant cos(PA) cut
            lV0Result[lNV0]->SetCutUseVarV0CosPA         ( kFALSE ) ;
            
            lV0Result[lNV0]->SetCutV0Radius              ( lCutsTight[i][4] ) ;
            lV0Result[lNV0]->SetCutMaxV0Radius           ( lCutsTight[i][10] ) ;
            
            //Miscellaneous
            lV0Result[lNV0]->SetCutProperLifetime        ( lCutsTight[i][5] ) ;
            lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( lCutsTight[i][11] ) ;
            lV0Result[lNV0]->SetCutMinTrackLength ( lCutsTight[i][6] ) ; //no cut here
            lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lCutsTight[i][7] ) ;
            lV0Result[lNV0]->SetCutTPCdEdx               ( 1e6 ) ; //no cut here
            lV0Result[lNV0]->SetCut276TeVLikedEdx        ( kTRUE ) ;
            lV0Result[lNV0]->SetCutArmenterosParameter               ( lCutsTight[i][9] ) ;
            
            //Add result to pool
            lNV0++;
        }
        
        //Loosening cuts one by one
        for(Int_t i = 0 ; i < lNPart ; i ++)
        {
            for(Int_t iCut = 0; iCut < lNCutsForSweep; iCut++)
            {
                //only proceed if cuts are actually different
                if( ( TMath::Abs( lCutsTight[i][iCut]-lCutsLoose[i][iCut] ) / lCutsLoose[i][iCut] < 0.01 ) && ( iCut != 3 ) ) continue;
                
                Int_t lNSweep = 12;
                for(Int_t iSweep = 1; iSweep <= lNSweep; iSweep++)
                {
                    Double_t lCutValue = lCutsTight[i][iCut] + (iSweep/(Double_t)lNSweep)*(lCutsLoose[i][iCut]-lCutsTight[i][iCut]);
                    
                    lV0Result[lNV0] = new AliV0Result( lV0Result[lTightForSweepIndex+i], Form( "%s_276Cuts_%s_%d", lParticleNameV0[i].Data(), lCutNameV0[iCut].Data(), iSweep ) );
                    
                    if(iCut ==  0 ) lV0Result[lNV0]->SetCutDCANegToPV            ( lCutValue ) ;
                    if(iCut ==  1 ) lV0Result[lNV0]->SetCutDCAPosToPV            ( lCutValue ) ;
                    if(iCut ==  2 ) lV0Result[lNV0]->SetCutDCAV0Daughters        ( lCutValue ) ;
                    if(iCut ==  3 )
                    {
                        lV0Result[lNV0]->SetCutV0CosPA               ( lCutValue ) ;
                        lV0Result[lNV0]->SetCutVarV0CosPA            ( parExp0Const[i][1]*(iSweep/(Double_t)lNSweep), parExp0Slope[i][1], parExp1Const[i][1]*(iSweep/(Double_t)lNSweep), parExp1Slope[i][1], TMath::ACos(lCutsTight[i][iCut]) + (iSweep/(Double_t)lNSweep)*(parConst[i][1]-TMath::ACos(lCutsTight[i][iCut])) ) ;
                    }
                    if(iCut ==  4 ) lV0Result[lNV0]->SetCutV0Radius              ( lCutValue ) ;
                    
                    //Miscellaneous
                    if(iCut ==  5 ) lV0Result[lNV0]->SetCutProperLifetime        ( lCutValue ) ;
                    if(iCut ==  6 ) lV0Result[lNV0]->SetCutMinTrackLength ( lCutValue ) ;
                    if(iCut ==  7 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lCutValue ) ;
                    if(iCut ==  8 )
                    {
                        lV0Result[lNV0]->SetCut276TeVLikedEdx        ( kTRUE ) ;
                        lV0Result[lNV0]->SetCutTPCdEdx               ( lCutValue ) ;
                    }
                    if(iCut ==  9 ) lV0Result[lNV0]->SetCutArmenterosParameter   ( lCutValue ) ;
                    if(iCut == 10 ) lV0Result[lNV0]->SetCutMaxV0Radius           ( lCutValue ) ;
                    if(iCut == 11 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( lCutValue ) ;
                    
                    //Print this variation, add to pool
                    //lV0Result[lNV0]->Print();
                    lNV0++;
                }
            }
        }
    } // end sweeps between loose and tight versions of cuts
    
    if(lSweepFullNumb > 0)
    {
        // centrality binning for sweeps
        Double_t lSweepCentBinLimits[] = {0, 10, 40, 60, 80, 90};
        Long_t lSweepCentBinNumb = sizeof(lSweepCentBinLimits)/sizeof(Double_t) - 1;
        
        //Mass binning for sweeps
        Double_t lNMassBins[] = {200, 200, 200};
        Double_t lMass[] = {0.498, 1.116, 1.116};
        Double_t lMassWindow[] = {0.15,0.1,0.1};
        
        //Loose cuts for sweeps
        Int_t lCentralForSweepIndex = lNV0;
        for(Int_t i = 0 ; i < lNPart ; i++)
        {
            //Central result for sweeps
            lV0Result[lNV0] = new AliV0Result( Form("%s_Central_ForFullSweep",lParticleNameV0[i].Data() ),lMassHypoV0[i],"",lSweepCentBinNumb,lSweepCentBinLimits, lPtbinnumbV0,lPtbinlimitsV0,lNMassBins[i],lMass[i]-lMassWindow[i],lMass[i]+lMassWindow[i]);
            
            //feeddown matrix
            if ( i!=0 ) lV0Result[lNV0] -> SetupFeeddownMatrix(lPtbinnumbXi,lPtbinlimitsXi);
            
            //Setters for V0 Cuts
            lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][1][ 0] ) ;
            lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][1][ 1] ) ;
            lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][1][ 2] ) ;
            lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][1][ 3] ) ;
            //Set Variable cut
            lV0Result[lNV0]->SetCutVarV0CosPA               ( parExp0Const[i][1], parExp0Slope[i][1], parExp1Const[i][1], parExp1Slope[i][1], parConst[i][1] ) ;
            
            lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][1][ 4] ) ;
            
            //Miscellaneous
            lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][1][ 5] ) ;
            lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( -1 ) ; //no cut here
            lV0Result[lNV0]->SetCutMinTrackLength ( lcutsV0[i][1][ 6] ) ;
            lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][1][ 7] ) ;
            lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][1][ 8] ) ;
            lV0Result[lNV0]->SetCutArmenterosParameter               ( lcutsV0[i][1][ 9] ) ;
            
            //Add result to pool
            lNV0++;
        }
        
        //Variable: V0 CosPA
        lMinV0CosPA = 0.97;
        lMaxV0CosPA = 1.0;
        //  lMinV0PA = 0;
        Double_t lMaxV0PA = TMath::ACos(lMinV0CosPA);
        
        for(Int_t i = 0 ; i < lNPart ; i ++)
        {
            for(Int_t iSweep = 0; iSweep <= lSweepFullNumb; iSweep++)
            {
                Double_t lCutValue = TMath::Cos( lMinV0PA + (iSweep/(Double_t)lSweepFullNumb)*(lMaxV0PA-lMinV0PA) );
                
                lV0Result[lNV0] = new AliV0Result( lV0Result[lCentralForSweepIndex+i], Form( "%s_%s_%d", lParticleNameV0[i].Data(), "V0CosPASweep", iSweep ) );
                
                lV0Result[lNV0]->SetCutUseVarV0CosPA( kFALSE ) ;
                lV0Result[lNV0]->SetCutV0CosPA( lCutValue );
                
                //Print this variation, add to pool
                lV0Result[lNV0]->Print();
                lNV0++;
            }
        }
    }
    
    // STEP 4: Creation of objects to be used in systematics
    // Optimized via use of copy constructors
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        for(Int_t iCut = 0 ; iCut < lNCutsForSyst ; iCut ++)
        {
            //LOOSE VARIATIONS
            //Create a new object from default
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%s",lParticleNameV0[i].Data(),lCutNameV0[iCut].Data(),lConfNameV0[0].Data()) );
            
            if(iCut ==  0 ) lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  1 ) lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  2 ) lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  3 )
            {
                lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][0][iCut] ) ;
                lV0Result[lNV0]->SetCutVarV0CosPA               ( parExp0Const[i][0], parExp0Slope[i][0], parExp1Const[i][0], parExp1Slope[i][0], parConst[i][0] ) ;
            }
            if(iCut ==  4 ) lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][0][iCut] ) ;
            
            //Miscellaneous
            if(iCut ==  5 ) lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  6 ) lV0Result[lNV0]->SetCutMinTrackLength ( lcutsV0[i][0][ iCut] ) ;
            if(iCut ==  7 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  8 ) lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][0][iCut] ) ;
            if(iCut ==  9 ) lV0Result[lNV0]->SetCutArmenterosParameter               ( lcutsV0[i][0][iCut] ) ;
            
            //Print this variation, add to pool
            lV0Result[lNV0]->Print();
            lNV0++;
            
            //TIGHT VARIATIONS
            //Create a new object from default
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_%s_%s",lParticleNameV0[i].Data(),lCutNameV0[iCut].Data(),lConfNameV0[2].Data()) );
            
            if(iCut ==  0 ) lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  1 ) lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  2 ) lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  3 )
            {
                lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][2][iCut] ) ;
                lV0Result[lNV0]->SetCutVarV0CosPA               ( parExp0Const[i][2], parExp0Slope[i][2], parExp1Const[i][2], parExp1Slope[i][2], parConst[i][2] ) ;
            }
            if(iCut ==  4 ) lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][2][iCut] ) ;
            
            //Miscellaneous
            if(iCut ==  5 ) lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  6 ) lV0Result[lNV0]->SetCutMinTrackLength ( lcutsV0[i][2][ iCut] ) ;
            if(iCut ==  7 ) lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  8 ) lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][2][iCut] ) ;
            if(iCut ==  9 ) lV0Result[lNV0]->SetCutArmenterosParameter               ( lcutsV0[i][2][iCut] ) ;
            
            //Print this variation, add to pool
            //lV0Result[lNV0]->Print();
            lNV0++;
        }
    }
    
    for (Int_t iconf = 0; iconf<lNV0; iconf++){
        //   cout<<"Adding config named"<<lV0Result[iconf]->GetName()<<endl;
        AddConfiguration(lV0Result[iconf]);
    }
    
    // cout<<"Added "<<lNV0<<" V0 configurations to output."<<endl;
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::AddStandardV0RadiusSweep()
//Meant to do radius sweep for debugging purposes
{
    //======================================================
    // V0 Configurations To use
    //======================================================
    
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimitsV0[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12};
    Long_t lPtbinnumbV0 = sizeof(lPtbinlimitsV0)/sizeof(Double_t) - 1;
    Double_t lPtbinlimitsXi[] ={0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10, 12, 14, 16, 19, 22, 25};
    Long_t lPtbinnumbXi = sizeof(lPtbinlimitsXi)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimitsV0[] = {0, 1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 85, 90};
    Long_t lCentbinnumbV0 = sizeof(lCentbinlimitsV0)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleNameV0[] =
    {
        "K0Short",
        "Lambda",
        "AntiLambda"
    };
    const Int_t lNPart = 3;
    TString lConfNameV0[] =
    {
        "Loose",
        "Central",
        "Tight"
    };
    const Int_t lNConf = 3;
    TString lCutNameV0[] =
    {
        "DCANegToPV",
        "DCAPosToPV",
        "DCAV0Daughters",
        "V0CosPA",
        "V0Radius",
        "ProperLifetime",
        "TrackLength",
        "LeastNbrCrsOvFind",
        "TPCdEdx",
        "APParameter",
        "V0RadiusMax",
        "LeastNbrCrsRows"
    };
    const Int_t lNCutsForSyst = 10;
    
    // STEP 2: Decide on a set of selections
    
    //1st index: Particle Species
    //2nd index: Loose / Central / Tight
    //3rd index: Number of selection (as ordered above)
    Double_t lcutsV0[lNPart][lNConf][lNCutsForSyst];
    
    //1st index: Particle Species: K0Short, Lambda, AntiLambda
    //2nd index: Loose / Central / Tight: 2%, 5%, 10% signal loss
    Double_t parExp0Const[lNPart][lNConf];
    Double_t parExp0Slope[lNPart][lNConf];
    Double_t parExp1Const[lNPart][lNConf];
    Double_t parExp1Slope[lNPart][lNConf];
    Double_t parConst[lNPart][lNConf];
    
    //=============================================================================================
    // K0SHORT V0 COS PA PARAMETRIZATION
    //---------------------------------------------------------------------------------------------
    //                       LOOSE                         CENTRAL                           TIGHT
    parExp0Const[0][0] =  0.20428;  parExp0Const[0][1] =  0.22692;  parExp0Const[0][2] =  0.28814;
    parExp0Slope[0][0] = -0.73728;  parExp0Slope[0][1] = -1.59317;  parExp0Slope[0][2] = -2.27069;
    parExp1Const[0][0] =  0.09887;  parExp1Const[0][1] =  0.05994;  parExp1Const[0][2] =  0.04320;
    parExp1Slope[0][0] = -0.02822;  parExp1Slope[0][1] = -0.26997;  parExp1Slope[0][2] = -0.29839;
    parConst[0][0] = -0.05302;      parConst[0][1] =  0.00907;      parConst[0][2] =  0.00704;
    //=============================================================================================
    
    //=============================================================================================
    // LAMBDA V0 COS PA PARAMETRIZATION
    //---------------------------------------------------------------------------------------------
    //                       LOOSE                         CENTRAL                           TIGHT
    parExp0Const[1][0] =  0.22775;  parExp0Const[1][1] =  0.36284;  parExp0Const[1][2] =  0.54877;
    parExp0Slope[1][0] = -1.11579;  parExp0Slope[1][1] = -1.87960;  parExp0Slope[1][2] = -2.72912;
    parExp1Const[1][0] =  0.06266;  parExp1Const[1][1] =  0.04543;  parExp1Const[1][2] =  0.03411;
    parExp1Slope[1][0] = -0.17086;  parExp1Slope[1][1] = -0.20447;  parExp1Slope[1][2] = -0.26965;
    parConst[1][0] =  0.01489;      parConst[1][1] =  0.01085;      parConst[1][2] =  0.00889;
    //=============================================================================================
    
    //=============================================================================================
    // ANTILAMBDA V0 COS PA PARAMETRIZATION
    //---------------------------------------------------------------------------------------------
    //                       LOOSE                         CENTRAL                           TIGHT
    parExp0Const[2][0] =  0.22667;  parExp0Const[2][1] =  0.35809;  parExp0Const[2][2] =  0.54114;
    parExp0Slope[2][0] = -0.93618;  parExp0Slope[2][1] = -1.93860;  parExp0Slope[2][2] = -2.71000;
    parExp1Const[2][0] =  0.06857;  parExp1Const[2][1] =  0.05306;  parExp1Const[2][2] =  0.03664;
    parExp1Slope[2][0] = -0.07015;  parExp1Slope[2][1] = -0.24518;  parExp1Slope[2][2] = -0.28124;
    parConst[2][0] = -0.00707;      parConst[2][1] =  0.01213;      parConst[2][2] =  0.00905;
    //=============================================================================================
    
    //================================================================================
    // K0SHORT SELECTIONS
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcutsV0[0][0][ 0] = 0.05;    lcutsV0[0][1][ 0] =  0.10; lcutsV0[0][2][0] = 0.17; //DCANegToPV
    lcutsV0[0][0][ 1] = 0.05;    lcutsV0[0][1][ 1] =  0.10; lcutsV0[0][2][1] = 0.17; //DCAPosToPV
    lcutsV0[0][0][ 2] = 0.95;    lcutsV0[0][1][ 2] =   0.8; lcutsV0[0][2][2] =  0.7; //DCAV0Daughters
    lcutsV0[0][0][ 3] = 0.95;    lcutsV0[0][1][ 3] =  0.95; lcutsV0[0][2][3] = 0.95; //V0CosPA
    lcutsV0[0][0][ 4] = 4.50;    lcutsV0[0][1][ 4] =  5.00; lcutsV0[0][2][4] = 5.50; //V0Radius
    lcutsV0[0][0][ 5] =   25;    lcutsV0[0][1][ 5] =    20; lcutsV0[0][2][5] =   15; //Proper Lifetime (in cm)
    lcutsV0[0][0][ 6] =   80;    lcutsV0[0][1][ 6] =    90; lcutsV0[0][2][6] =  100; //Track Length
    lcutsV0[0][0][ 7] =  0.7;    lcutsV0[0][1][ 7] =   0.8; lcutsV0[0][2][7] = 0.85; //Least Ratio CrdRows/Findable
    lcutsV0[0][0][ 8] =  4.0;    lcutsV0[0][1][ 8] =   3.0; lcutsV0[0][2][8] =  2.5; //TPC dE/dx
    lcutsV0[0][0][ 9] = 0.18;    lcutsV0[0][1][ 9] =  0.20; lcutsV0[0][2][9] = 0.22; //AP Parameter
    //================================================================================
    
    //================================================================================
    // LAMBDA SELECTIONS
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcutsV0[1][0][ 0] = 0.10;    lcutsV0[1][1][ 0] =  0.25; lcutsV0[1][2][0] = 0.40; //DCANegToPV
    lcutsV0[1][0][ 1] = 0.08;    lcutsV0[1][1][ 1] =  0.10; lcutsV0[1][2][1] = 0.13; //DCAPosToPV
    lcutsV0[1][0][ 2] =  1.0;    lcutsV0[1][1][ 2] =   0.8; lcutsV0[1][2][2] = 0.65; //DCAV0Daughters
    lcutsV0[1][0][ 3] = 0.97;    lcutsV0[1][1][ 3] =  0.98; lcutsV0[1][2][3] = 0.99; //V0CosPA
    lcutsV0[1][0][ 4] = 4.00;    lcutsV0[1][1][ 4] =  5.00; lcutsV0[1][2][4] = 6.00; //V0Radius
    lcutsV0[1][0][ 5] =   30;    lcutsV0[1][1][ 5] =    25; lcutsV0[1][2][5] =   20; //Proper Lifetime (in cm)
    lcutsV0[1][0][ 6] =   80;    lcutsV0[1][1][ 6] =    90; lcutsV0[1][2][6] =  100; //Track Length
    lcutsV0[1][0][ 7] =  0.7;    lcutsV0[1][1][ 7] =   0.8; lcutsV0[1][2][7] = 0.85; //Least Ratio CrdRows/Findable
    lcutsV0[1][0][ 8] =  4.0;    lcutsV0[1][1][ 8] =   3.0; lcutsV0[1][2][8] =  2.5; //TPC dE/dx
    lcutsV0[1][0][ 9] = 0.18;    lcutsV0[1][1][ 9] =  0.20; lcutsV0[1][2][9] = 0.22; //AP Parameter
    //================================================================================
    
    //================================================================================
    // ANTILAMBDA SELECTIONS
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcutsV0[2][0][ 0] = 0.08;    lcutsV0[2][1][ 0] =  0.10; lcutsV0[2][2][0] = 0.13; //DCANegToPV
    lcutsV0[2][0][ 1] = 0.10;    lcutsV0[2][1][ 1] =  0.25; lcutsV0[2][2][1] = 0.40; //DCAPosToPV
    lcutsV0[2][0][ 2] =  1.0;    lcutsV0[2][1][ 2] =   0.8; lcutsV0[2][2][2] = 0.65; //DCAV0Daughters
    lcutsV0[2][0][ 3] = 0.97;    lcutsV0[2][1][ 3] =  0.98; lcutsV0[2][2][3] = 0.99; //V0CosPA
    lcutsV0[2][0][ 4] = 4.00;    lcutsV0[2][1][ 4] =  5.00; lcutsV0[2][2][4] = 6.00; //V0Radius
    lcutsV0[2][0][ 5] =   30;    lcutsV0[2][1][ 5] =    25; lcutsV0[2][2][5] =   20; //Proper Lifetime (in cm)
    lcutsV0[2][0][ 6] =   80;    lcutsV0[2][1][ 6] =    90; lcutsV0[2][2][6] =  100; //Track Length
    lcutsV0[2][0][ 7] =  0.7;    lcutsV0[2][1][ 7] =   0.8; lcutsV0[2][2][7] = 0.85; //Least Ratio CrdRows/Findable
    lcutsV0[2][0][ 8] =  4.0;    lcutsV0[2][1][ 8] =   3.0; lcutsV0[2][2][8] =  2.5; //TPC dE/dx
    lcutsV0[2][0][ 9] = 0.18;    lcutsV0[2][1][ 9] =  0.20; lcutsV0[2][2][9] = 0.22; //AP Parameter
    //================================================================================
    
    //STEP 3: Creation of output objects
    
    //Map to mass hypothesis
    AliV0Result::EMassHypo lMassHypoV0[lNPart];
    lMassHypoV0[0] = AliV0Result::kK0Short;
    lMassHypoV0[1] = AliV0Result::kLambda;
    lMassHypoV0[2] = AliV0Result::kAntiLambda;
    
    //Array of results
    AliV0Result *lV0Result[3000];
    Long_t lNV0 = 0;
    
    //Central results: Stored in indices 0, 1, 2 (careful!)
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        //Central result, customized binning: the one to use, usually
        lV0Result[lNV0] = new AliV0Result( Form("%s_Central",lParticleNameV0[i].Data() ),lMassHypoV0[i],"",lCentbinnumbV0,lCentbinlimitsV0, lPtbinnumbV0,lPtbinlimitsV0);
        
        //feeddown matrix
        if ( i!=0 ) lV0Result[lNV0] -> SetupFeeddownMatrix(lPtbinnumbXi,lPtbinlimitsXi);
        
        //Setters for V0 Cuts
        lV0Result[lNV0]->SetCutDCANegToPV            ( lcutsV0[i][1][ 0] ) ;
        lV0Result[lNV0]->SetCutDCAPosToPV            ( lcutsV0[i][1][ 1] ) ;
        lV0Result[lNV0]->SetCutDCAV0Daughters        ( lcutsV0[i][1][ 2] ) ;
        lV0Result[lNV0]->SetCutV0CosPA               ( lcutsV0[i][1][ 3] ) ;
        //Set Variable cut
        lV0Result[lNV0]->SetCutVarV0CosPA               ( parExp0Const[i][1], parExp0Slope[i][1], parExp1Const[i][1], parExp1Slope[i][1], parConst[i][1] ) ;
        
        lV0Result[lNV0]->SetCutV0Radius              ( lcutsV0[i][1][ 4] ) ;
        
        //Miscellaneous
        lV0Result[lNV0]->SetCutProperLifetime        ( lcutsV0[i][1][ 5] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRows ( -1 ) ; //no cut here
        lV0Result[lNV0]->SetCutMinTrackLength ( lcutsV0[i][1][ 6] ) ;
        lV0Result[lNV0]->SetCutLeastNumberOfCrossedRowsOverFindable               ( lcutsV0[i][1][ 7] ) ;
        lV0Result[lNV0]->SetCutTPCdEdx               ( lcutsV0[i][1][ 8] ) ;
        lV0Result[lNV0]->SetCutArmenterosParameter               ( lcutsV0[i][1][ 9] ) ;
        
        //Add result to pool
        lNV0++;
    }
    
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_Cowboy", lParticleNameV0[i].Data() ) );
        
        
        //lV0Result[lNV0] -> SetCutIsCowboy(1); //!!!!!!!!!
        //Add result to pool
        lNV0++;
        
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_Sailor", lParticleNameV0[i].Data() ) );
        //  lV0Result[lNV0] -> SetCutIsCowboy(-1); //!!!!!!!!
        //Add result to pool
        lNV0++;
    }
    
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_LooseRadius", lParticleNameV0[i].Data() ) );
        lV0Result[lNV0] -> SetCutV0Radius   ( 0.9 );
        lV0Result[lNV0] -> SetCutMaxV0Radius( 5.0 );
        //Add result to pool
        lNV0++;
        
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_LooseRadius_Cowboy", lParticleNameV0[i].Data() ) );
        lV0Result[lNV0] -> SetCutV0Radius   ( 0.9 );
        lV0Result[lNV0] -> SetCutMaxV0Radius( 5.0 );
        // lV0Result[lNV0] -> SetCutIsCowboy(1);  //!!!!!!
        //Add result to pool
        lNV0++;
        
        //Create a new object from default
        lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_LooseRadius_Sailor", lParticleNameV0[i].Data() ) );
        lV0Result[lNV0] -> SetCutV0Radius   ( 0.9 );
        lV0Result[lNV0] -> SetCutMaxV0Radius( 5.0 );
        //  lV0Result[lNV0] -> SetCutIsCowboy(-1);
        //Add result to pool
        lNV0++;
    }
    
    //Rapidity sweep
    Double_t lRadii[] = {0.9, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10,15,20,30};
    Int_t lNRadii = sizeof(lRadii)/sizeof(Double_t);
    
    for(Int_t i = 0 ; i < lNPart ; i ++)
    {
        for(Int_t ir=0; ir<lNRadii-1; ir++)
        {
            //Create a new object from default
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_RadiusSweep_%.1f_%.1f",lParticleNameV0[i].Data(),lRadii[ir], lRadii[ir+1]) );
            lV0Result[lNV0] -> SetCutV0Radius   ( lRadii[ir  ] );
            lV0Result[lNV0] -> SetCutMaxV0Radius( lRadii[ir+1] );
            //Add result to pool
            lNV0++;
            
            //Create a new object from default
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_RadiusSweep_Cowboy_%.1f_%.1f",lParticleNameV0[i].Data(),lRadii[ir], lRadii[ir+1]) );
            lV0Result[lNV0] -> SetCutV0Radius   ( lRadii[ir  ] );
            lV0Result[lNV0] -> SetCutMaxV0Radius( lRadii[ir+1] );
            //   lV0Result[lNV0] -> SetCutIsCowboy(1);
            //Add result to pool
            lNV0++;
            
            //Create a new object from default
            lV0Result[lNV0] = new AliV0Result( lV0Result[i], Form("%s_RadiusSweep_Sailor_%.1f_%.1f",lParticleNameV0[i].Data(),lRadii[ir], lRadii[ir+1]) );
            lV0Result[lNV0] -> SetCutV0Radius   ( lRadii[ir  ] );
            lV0Result[lNV0] -> SetCutMaxV0Radius( lRadii[ir+1] );
            //  lV0Result[lNV0] -> SetCutIsCowboy(-1);
            //Add result to pool
            lNV0++;
        }
    }
    
    for (Int_t iconf = 0; iconf<lNV0; iconf++){
        //   cout<<"Radius sweep: adding config named "<<lV0Result[iconf]->GetName()<<endl;
        AddConfiguration(lV0Result[iconf]);
    }
    //cout<<"Added "<<lNV0<<" V0 configurations to output."<<endl;
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::AddStandardCascadeConfiguration(Bool_t lUseFull, Bool_t lDoSystematics)
//Meant to add some standard cascade analysis Configuration + its corresponding systematics
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimits[] = {0.4, 0.5, 0.6,
        0.7,0.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
        2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,
        4.4,4.5,4.6,4.8,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.,11.,12.};
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimits[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleName[] = {"XiMinus", "XiPlus",  "OmegaMinus", "OmegaPlus"};
    TString lConfName[]     = {"Loose",   "Central", "Tight"     };
    TString lCutName[]      = {
        "DCANegToPV", //1
        "DCAPosToPV", //2
        "DCAV0Daughters", //3
        "V0Radius", //4
        "DCAV0ToPV", //5
        "V0Mass", //6
        "DCABachToPV", //7
        "DCACascDaughters", //8
        "CascRadius", //9
        "ProperLifetime", //10
        "ProperLifetimeV0", //11
        "MinLength", //12
        "TPCdEdx", //13
        "Competing", //14
        "DCA3DCascToPV" //15
    };
    
    // STEP 2: Decide on a set of selections
    
    //1st index: Particle Species
    //2nd index: Loose / Central / Tight
    //3rd index: Number of selection (as ordered above)
    Double_t lcuts[4][3][15];
    
    //N.B.: These are mostly symmetric, except for the proper lifetime, which is different
    //      for the two particle species. Another potential improvement could be asymmetric
    //      DCA selections for the Neg / Pos tracks for the (anti)Lambda decay, as decay
    //      kinematics would prefer having the pion daughter with a larger DCA.
    
    Int_t lIdx = 0;
    
    //================================================================================
    // XIMINUS SELECTIONS
    lIdx = 0; //Master XiMinus Index
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcuts[lIdx][0][ 0] = 0.10;    lcuts[lIdx][1][ 0] =  0.20; lcuts[lIdx][2][ 0] =  0.30; //DCANegToPV 1
    lcuts[lIdx][0][ 1] = 0.10;    lcuts[lIdx][1][ 1] =  0.20; lcuts[lIdx][2][ 1] =  0.30; //DCAPostToPV 2
    lcuts[lIdx][0][ 2] =  1.2;    lcuts[lIdx][1][ 2] =   1.0; lcuts[lIdx][2][ 2] =   0.8; //DCAV0Daughters 3
    lcuts[lIdx][0][ 3] = 2.00;    lcuts[lIdx][1][ 3] =  3.00; lcuts[lIdx][2][ 3] =   4.0; //V0Radius 4
    lcuts[lIdx][0][ 4] = 0.05;    lcuts[lIdx][1][ 4] =   0.1; lcuts[lIdx][2][ 4] =  0.15; //DCAV0ToPV 5
    lcuts[lIdx][0][ 5] =0.006;    lcuts[lIdx][1][ 5] = 0.005; lcuts[lIdx][2][ 5] = 0.004; //V0Mass 6
    lcuts[lIdx][0][ 6] = 0.05;    lcuts[lIdx][1][ 6] =  0.10; lcuts[lIdx][2][ 6] =  0.15; //DCABachToPV 7
    lcuts[lIdx][0][ 7] = 1.20;    lcuts[lIdx][1][ 7] =   1.0; lcuts[lIdx][2][ 7] =   0.8; //DCACascDaughters 8
    lcuts[lIdx][0][ 8] =  0.8;    lcuts[lIdx][1][ 8] =   1.2; lcuts[lIdx][2][ 8] =  1.60; //CascRadius 9
    lcuts[lIdx][0][ 9] = 17.5;    lcuts[lIdx][1][ 9] =  15.0; lcuts[lIdx][2][ 9] =  12.5; //ProperLifetime 10
    lcuts[lIdx][0][10] = 40.0;    lcuts[lIdx][1][10] =  30.0; lcuts[lIdx][2][10] =  20.0; //ProperLifetimeV0 11
    lcuts[lIdx][0][11] =   80;    lcuts[lIdx][1][11] =    90; lcuts[lIdx][2][11] =   100; //MinimumTrackLength 12
    lcuts[lIdx][0][12] =    4;    lcuts[lIdx][1][12] =     3; lcuts[lIdx][2][12] =   2.5; //TPCdEdx 13
    lcuts[lIdx][0][13] =  0.0;    lcuts[lIdx][1][13] = 0.008; lcuts[lIdx][2][13] = 0.010; //Competing 14
    lcuts[lIdx][0][14] =  1.2;    lcuts[lIdx][1][14] =   0.8; lcuts[lIdx][2][14] =   0.6; //3D DCA Cascade To PV
    //================================================================================
    
    //================================================================================
    // XIPLUS SELECTIONS
    lIdx = 1; //Master XiPlus Index
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcuts[lIdx][0][ 0] = 0.10;    lcuts[lIdx][1][ 0] =  0.20; lcuts[lIdx][2][ 0] =  0.30; //DCANegToPV 1
    lcuts[lIdx][0][ 1] = 0.10;    lcuts[lIdx][1][ 1] =  0.20; lcuts[lIdx][2][ 1] =  0.30; //DCAPostToPV 2
    lcuts[lIdx][0][ 2] =  1.2;    lcuts[lIdx][1][ 2] =   1.0; lcuts[lIdx][2][ 2] =   0.8; //DCAV0Daughters 3
    lcuts[lIdx][0][ 3] = 2.00;    lcuts[lIdx][1][ 3] =  3.00; lcuts[lIdx][2][ 3] =   4.0; //V0Radius 4
    lcuts[lIdx][0][ 4] = 0.05;    lcuts[lIdx][1][ 4] =   0.1; lcuts[lIdx][2][ 4] =  0.15; //DCAV0ToPV 5
    lcuts[lIdx][0][ 5] =0.006;    lcuts[lIdx][1][ 5] = 0.005; lcuts[lIdx][2][ 5] = 0.004; //V0Mass 6
    lcuts[lIdx][0][ 6] = 0.05;    lcuts[lIdx][1][ 6] =  0.10; lcuts[lIdx][2][ 6] =  0.15; //DCABachToPV 7
    lcuts[lIdx][0][ 7] = 1.20;    lcuts[lIdx][1][ 7] =   1.0; lcuts[lIdx][2][ 7] =   0.8; //DCACascDaughters 8
    lcuts[lIdx][0][ 8] =  0.8;    lcuts[lIdx][1][ 8] =   1.2; lcuts[lIdx][2][ 8] =  1.60; //CascRadius 9
    lcuts[lIdx][0][ 9] = 17.5;    lcuts[lIdx][1][ 9] =  15.0; lcuts[lIdx][2][ 9] =  12.5; //ProperLifetime 10
    lcuts[lIdx][0][10] = 40.0;    lcuts[lIdx][1][10] =  30.0; lcuts[lIdx][2][10] =  20.0; //ProperLifetimeV0 11
    lcuts[lIdx][0][11] =   80;    lcuts[lIdx][1][11] =    90; lcuts[lIdx][2][11] =   100; //MinimumTrackLength 12
    lcuts[lIdx][0][12] =    4;    lcuts[lIdx][1][12] =     3; lcuts[lIdx][2][12] =   2.5; //TPCdEdx 13
    lcuts[lIdx][0][13] =  0.0;    lcuts[lIdx][1][13] = 0.008; lcuts[lIdx][2][13] = 0.010; //Competing 14
    lcuts[lIdx][0][14] =  1.2;    lcuts[lIdx][1][14] =   0.8; lcuts[lIdx][2][14] =   0.6; //3D DCA Cascade To PV
    //================================================================================
    
    //================================================================================
    // OMEGAMINUS SELECTIONS
    lIdx = 2; //Master OmegaMinus Index
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcuts[lIdx][0][ 0] = 0.10;    lcuts[lIdx][1][ 0] =  0.20; lcuts[lIdx][2][ 0] =  0.30; //DCANegToPV 1
    lcuts[lIdx][0][ 1] = 0.10;    lcuts[lIdx][1][ 1] =  0.20; lcuts[lIdx][2][ 1] =  0.30; //DCAPostToPV 2
    lcuts[lIdx][0][ 2] =  1.2;    lcuts[lIdx][1][ 2] =   1.0; lcuts[lIdx][2][ 2] =   0.8; //DCAV0Daughters 3
    lcuts[lIdx][0][ 3] = 2.00;    lcuts[lIdx][1][ 3] =  3.00; lcuts[lIdx][2][ 3] =   4.0; //V0Radius 4
    lcuts[lIdx][0][ 4] = 0.05;    lcuts[lIdx][1][ 4] =   0.1; lcuts[lIdx][2][ 4] =  0.15; //DCAV0ToPV 5
    lcuts[lIdx][0][ 5] =0.006;    lcuts[lIdx][1][ 5] = 0.005; lcuts[lIdx][2][ 5] = 0.004; //V0Mass 6
    lcuts[lIdx][0][ 6] = 0.05;    lcuts[lIdx][1][ 6] =  0.10; lcuts[lIdx][2][ 6] =  0.15; //DCABachToPV 7
    lcuts[lIdx][0][ 7] = 1.00;    lcuts[lIdx][1][ 7] =   0.6; lcuts[lIdx][2][ 7] =   0.5; //DCACascDaughters 8
    lcuts[lIdx][0][ 8] =  0.6;    lcuts[lIdx][1][ 8] =   1.0; lcuts[lIdx][2][ 8] =  1.40; //CascRadius 9
    lcuts[lIdx][0][ 9] = 14.0;    lcuts[lIdx][1][ 9] =  12.0; lcuts[lIdx][2][ 9] =  10.0; //ProperLifetime 10
    lcuts[lIdx][0][10] = 40.0;    lcuts[lIdx][1][10] =  30.0; lcuts[lIdx][2][10] =  20.0; //ProperLifetimeV0 11
    lcuts[lIdx][0][11] =   80;    lcuts[lIdx][1][11] =    90; lcuts[lIdx][2][11] =   100; //MinimumTrackLength 12
    lcuts[lIdx][0][12] =    4;    lcuts[lIdx][1][12] =     3; lcuts[lIdx][2][12] =   2.5; //TPCdEdx 13
    lcuts[lIdx][0][13] =  0.0;    lcuts[lIdx][1][13] = 0.008; lcuts[lIdx][2][13] = 0.010; //Competing 14
    lcuts[lIdx][0][14] =  0.8;    lcuts[lIdx][1][14] =   0.6; lcuts[lIdx][2][14] =   0.5; //3D DCA Cascade To PV
    //================================================================================
    
    //================================================================================
    // OMEGAPLUS SELECTIONS
    lIdx = 3; //Master OmegaPlus Index
    //--------------------------------------------------------------------------------
    //                  LOOSE                        CENTRAL                   TIGHT
    lcuts[lIdx][0][ 0] = 0.10;    lcuts[lIdx][1][ 0] =  0.20; lcuts[lIdx][2][ 0] =  0.30; //DCANegToPV 1
    lcuts[lIdx][0][ 1] = 0.10;    lcuts[lIdx][1][ 1] =  0.20; lcuts[lIdx][2][ 1] =  0.30; //DCAPostToPV 2
    lcuts[lIdx][0][ 2] =  1.2;    lcuts[lIdx][1][ 2] =   1.0; lcuts[lIdx][2][ 2] =   0.8; //DCAV0Daughters 3
    lcuts[lIdx][0][ 3] = 2.00;    lcuts[lIdx][1][ 3] =  3.00; lcuts[lIdx][2][ 3] =   4.0; //V0Radius 4
    lcuts[lIdx][0][ 4] = 0.05;    lcuts[lIdx][1][ 4] =   0.1; lcuts[lIdx][2][ 4] =  0.15; //DCAV0ToPV 5
    lcuts[lIdx][0][ 5] =0.006;    lcuts[lIdx][1][ 5] = 0.005; lcuts[lIdx][2][ 5] = 0.004; //V0Mass 6
    lcuts[lIdx][0][ 6] = 0.05;    lcuts[lIdx][1][ 6] =  0.10; lcuts[lIdx][2][ 6] =  0.15; //DCABachToPV 7
    lcuts[lIdx][0][ 7] = 1.00;    lcuts[lIdx][1][ 7] =   0.6; lcuts[lIdx][2][ 7] =   0.5; //DCACascDaughters 8
    lcuts[lIdx][0][ 8] =  0.6;    lcuts[lIdx][1][ 8] =   1.0; lcuts[lIdx][2][ 8] =  1.40; //CascRadius 9
    lcuts[lIdx][0][ 9] = 14.0;    lcuts[lIdx][1][ 9] =  12.0; lcuts[lIdx][2][ 9] =  10.0; //ProperLifetime 10
    lcuts[lIdx][0][10] = 40.0;    lcuts[lIdx][1][10] =  30.0; lcuts[lIdx][2][10] =  20.0; //ProperLifetimeV0 11
    lcuts[lIdx][0][11] =   80;    lcuts[lIdx][1][11] =    90; lcuts[lIdx][2][11] =   100; //MinimumTrackLength 12
    lcuts[lIdx][0][12] =    4;    lcuts[lIdx][1][12] =     3; lcuts[lIdx][2][12] =   2.5; //TPCdEdx 13
    lcuts[lIdx][0][13] =  0.0;    lcuts[lIdx][1][13] = 0.008; lcuts[lIdx][2][13] = 0.010; //Competing 14
    lcuts[lIdx][0][14] =  0.8;    lcuts[lIdx][1][14] =   0.6; lcuts[lIdx][2][14] =   0.5; //3D DCA Cascade To PV
    //================================================================================
    
    //STEP 3: Creation of output objects
    
    //Just a counter and one array, please. Nothing else needed
    AliCascadeResult *lCascadeResult[600];
    Long_t lN = 0;
    
    //Map to mass hypothesis
    AliCascadeResult::EMassHypo lMassHypo[4];
    lMassHypo[0] = AliCascadeResult::kXiMinus;
    lMassHypo[1] = AliCascadeResult::kXiPlus;
    lMassHypo[2] = AliCascadeResult::kOmegaMinus;
    lMassHypo[3] = AliCascadeResult::kOmegaPlus;
    
    
    //Central results: Stored in indices 0, 1, 2, 3 (careful!)
    for(Int_t i = 0 ; i < 4 ; i ++){
        //Central result, customized binning: the one to use, usually
        lCascadeResult[lN] = new AliCascadeResult( Form("%s_Central",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits);
        
        //This is MC: generate profile for G3/F (if ever needed)
        lCascadeResult[lN] -> InitializeProtonProfile();
        
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( lcuts[i][1][ 0] ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( lcuts[i][1][ 1] ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        ( lcuts[i][1][ 2] ) ;
        lCascadeResult[lN]->SetCutV0Radius              ( lcuts[i][1][ 3] ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( lcuts[i][1][ 4] ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( lcuts[i][1][ 5] ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( lcuts[i][1][ 6] ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( lcuts[i][1][ 7] ) ;
        lCascadeResult[lN]->SetCutVarDCACascDau ( TMath::Exp(0.0470076), -0.917006, 0, 1, 0.5 );
        lCascadeResult[lN]->SetCutCascRadius            ( lcuts[i][1][ 8] ) ;
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( lcuts[i][1][ 9] ) ;
        lCascadeResult[lN]->SetCutMaxV0Lifetime         ( lcuts[i][1][10] ) ;
        lCascadeResult[lN]->SetCutMinTrackLength        ( lcuts[i][1][11] ) ;
        lCascadeResult[lN]->SetCutLeastNumberOfClusters( -1 );
        lCascadeResult[lN]->SetCutTPCdEdx               ( lcuts[i][1][12] ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( lcuts[i][1][13] ) ;
        lCascadeResult[lN]->SetCutDCACascadeToPV        ( lcuts[i][1][14] ) ;
        
        //Parametric angle cut initializations
        //V0 cosine of pointing angle
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarV0CosPA            (TMath::Exp(10.853),
                                                         -25.0322,
                                                         TMath::Exp(-0.843948),
                                                         -0.890794,
                                                         0.057553);
        
        //Cascade cosine of pointing angle
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.95 ) ; //+variable
        if(i < 2){
            lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(4.86664),
                                                             -10.786,
                                                             TMath::Exp(-1.33411),
                                                             -0.729825,
                                                             0.0695724);
        }
        if(i >= 2){
            lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(   12.8752),
                                                             -21.522,
                                                             TMath::Exp( -1.49906),
                                                             -0.813472,
                                                             0.0480962);
        }
        
        //BB cosine of pointing angle
        lCascadeResult[lN]->SetCutBachBaryonCosPA       ( TMath::Cos(0.04) ) ; //+variable
        lCascadeResult[lN]->SetCutVarBBCosPA            (TMath::Exp(-2.29048),
                                                         -20.2016,
                                                         TMath::Exp(-2.9581),
                                                         -0.649153,
                                                         0.00526455);
        
        //Add result to pool
        lN++;
    }
    if ( lUseFull ) {
        //Central Full results: Stored in indices 4, 5, 6, 7 (careful!)
        for(Int_t i = 0 ; i < 4 ; i ++){
            lCascadeResult[lN] = new AliCascadeResult( Form("%s_Central_Full",lParticleName[i].Data() ),lMassHypo[i]);
            
            //This is MC: generate profile for G3/F (if ever needed)
            lCascadeResult[lN] -> InitializeProtonProfile();
            
            //Setters for V0 Cuts
            lCascadeResult[lN]->SetCutDCANegToPV            ( lcuts[i][1][ 0] ) ;
            lCascadeResult[lN]->SetCutDCAPosToPV            ( lcuts[i][1][ 1] ) ;
            lCascadeResult[lN]->SetCutDCAV0Daughters        ( lcuts[i][1][ 2] ) ;
            lCascadeResult[lN]->SetCutV0Radius              ( lcuts[i][1][ 3] ) ;
            //Setters for Cascade Cuts
            lCascadeResult[lN]->SetCutDCAV0ToPV             ( lcuts[i][1][ 4] ) ;
            lCascadeResult[lN]->SetCutV0Mass                ( lcuts[i][1][ 5] ) ;
            lCascadeResult[lN]->SetCutDCABachToPV           ( lcuts[i][1][ 6] ) ;
            lCascadeResult[lN]->SetCutDCACascDaughters      ( lcuts[i][1][ 7] ) ;
            lCascadeResult[lN]->SetCutVarDCACascDau ( TMath::Exp(0.0470076), -0.917006, 0, 1, 0.5 );
            lCascadeResult[lN]->SetCutCascRadius            ( lcuts[i][1][ 8] ) ;
            //Miscellaneous
            lCascadeResult[lN]->SetCutProperLifetime        ( lcuts[i][1][ 9] ) ;
            lCascadeResult[lN]->SetCutMaxV0Lifetime         ( lcuts[i][1][10] ) ;
            lCascadeResult[lN]->SetCutMinTrackLength        ( lcuts[i][1][11] ) ;
            lCascadeResult[lN]->SetCutLeastNumberOfClusters( -1 );
            lCascadeResult[lN]->SetCutTPCdEdx               ( lcuts[i][1][12] ) ;
            lCascadeResult[lN]->SetCutXiRejection           ( lcuts[i][1][13] ) ;
            lCascadeResult[lN]->SetCutDCACascadeToPV        ( lcuts[i][1][14] ) ;
            
            //Parametric angle cut initializations
            //V0 cosine of pointing angle
            lCascadeResult[lN]->SetCutV0CosPA               ( 0.95 ) ; //+variable
            lCascadeResult[lN]->SetCutVarV0CosPA            (TMath::Exp(10.853),
                                                             -25.0322,
                                                             TMath::Exp(-0.843948),
                                                             -0.890794,
                                                             0.057553);
            
            //Cascade cosine of pointing angle
            lCascadeResult[lN]->SetCutCascCosPA             ( 0.95 ) ; //+variable
            if(i < 2){
                lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(4.86664),
                                                                 -10.786,
                                                                 TMath::Exp(-1.33411),
                                                                 -0.729825,
                                                                 0.0695724);
            }
            if(i >= 2){
                lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(   12.8752),
                                                                 -21.522,
                                                                 TMath::Exp( -1.49906),
                                                                 -0.813472,
                                                                 0.0480962);
            }
            
            //BB cosine of pointing angle
            lCascadeResult[lN]->SetCutBachBaryonCosPA       ( TMath::Cos(0.04) ) ; //+variable
            lCascadeResult[lN]->SetCutVarBBCosPA            (TMath::Exp(-2.29048),
                                                             -20.2016,
                                                             TMath::Exp(-2.9581),
                                                             -0.649153,
                                                             0.00526455);
            
            //Add result to pool
            lN++;
        }
    }
    
    //Explore restricted rapidity range check
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_OpenCosPA",lParticleName[i].Data() ) );
        
        lCascadeResult[lN] -> SetCutUseVarCascCosPA(kFALSE);
        lCascadeResult[lN] -> SetCutUseVarV0CosPA(kFALSE);
        
        //Add result to pool
        lN++;
    }
    
    //Explore restricted rapidity range check
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_Central_y03",lParticleName[i].Data() ) );
        
        lCascadeResult[lN] -> SetCutMinRapidity(-0.3);
        lCascadeResult[lN] -> SetCutMaxRapidity(+0.3);
        
        //Add result to pool
        lN++;
    }
    
    //Explore No PDG association
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_NoMCInfoUsed",lParticleName[i].Data() ) );
        
        lCascadeResult[lN] -> SetCutMCUseMCProperties(kFALSE);
        lCascadeResult[lN] -> SetCutMCPhysicalPrimary(kFALSE);
        lCascadeResult[lN] -> SetCutMCPDGCodeAssociation(kFALSE);
        
        //Add result to pool
        lN++;
    }
    
    //Explore TOF info use
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_AtLeastOneTOF",lParticleName[i].Data() ) );
        lCascadeResult[lN] -> SetCutAtLeastOneTOF(kTRUE);
        //Add result to pool
        lN++;
    }
    
    //Explore no TPC dEdx use
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_NodEdx",lParticleName[i].Data() ) );
        lCascadeResult[lN] -> SetCutTPCdEdx(1e+6);
        //Add result to pool
        lN++;
    }
    
    //Require ITS refit (will lose tons of signal)
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitTracks",lParticleName[i].Data() ) );
        lCascadeResult[lN] -> SetCutUseITSRefitTracks(kTRUE);
        //Add result to pool
        lN++;
    }
    
    //Require ITS refit (will lose tons of signal)
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitTracks_NoAssoc",lParticleName[i].Data() ) );
        lCascadeResult[lN] -> SetCutUseITSRefitTracks(kTRUE);
        lCascadeResult[lN] -> SetCutMCUseMCProperties(kFALSE);
        lCascadeResult[lN] -> SetCutMCPhysicalPrimary(kFALSE);
        lCascadeResult[lN] -> SetCutMCPDGCodeAssociation(kFALSE);
        //Add result to pool
        lN++;
    }
    
    if( lDoSystematics ){
        //ITS refit requirement map:
        //  [NPB]
        //1  100
        //2  010
        //3  001
        //4  110
        //5  101
        //6  011
        
        //Require ITS refit (will lose tons of signal)
        for(Int_t i = 0 ; i < 4 ; i ++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB100",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitNegative(kTRUE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB010",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitPositive(kTRUE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB001",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitBachelor(kTRUE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB110",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitNegative(kTRUE);
            lCascadeResult[lN] -> SetCutUseITSRefitPositive(kTRUE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB101",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitNegative(kTRUE);
            lCascadeResult[lN] -> SetCutUseITSRefitBachelor(kTRUE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB011",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitPositive(kTRUE);
            lCascadeResult[lN] -> SetCutUseITSRefitBachelor(kTRUE);
            lN++;
        }
        
        //Require ITS refit (will lose tons of signal)
        //No perfect MC association
        for(Int_t i = 0 ; i < 4 ; i ++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB100_NoAssoc",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitNegative(kTRUE);
            lCascadeResult[lN] -> SetCutMCUseMCProperties(kFALSE);
            lCascadeResult[lN] -> SetCutMCPhysicalPrimary(kFALSE);
            lCascadeResult[lN] -> SetCutMCPDGCodeAssociation(kFALSE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB010_NoAssoc",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitPositive(kTRUE);
            lCascadeResult[lN] -> SetCutMCUseMCProperties(kFALSE);
            lCascadeResult[lN] -> SetCutMCPhysicalPrimary(kFALSE);
            lCascadeResult[lN] -> SetCutMCPDGCodeAssociation(kFALSE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB001_NoAssoc",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitBachelor(kTRUE);
            lCascadeResult[lN] -> SetCutMCUseMCProperties(kFALSE);
            lCascadeResult[lN] -> SetCutMCPhysicalPrimary(kFALSE);
            lCascadeResult[lN] -> SetCutMCPDGCodeAssociation(kFALSE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB110_NoAssoc",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitNegative(kTRUE);
            lCascadeResult[lN] -> SetCutUseITSRefitPositive(kTRUE);
            lCascadeResult[lN] -> SetCutMCUseMCProperties(kFALSE);
            lCascadeResult[lN] -> SetCutMCPhysicalPrimary(kFALSE);
            lCascadeResult[lN] -> SetCutMCPDGCodeAssociation(kFALSE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB101_NoAssoc",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitNegative(kTRUE);
            lCascadeResult[lN] -> SetCutUseITSRefitBachelor(kTRUE);
            lCascadeResult[lN] -> SetCutMCUseMCProperties(kFALSE);
            lCascadeResult[lN] -> SetCutMCPhysicalPrimary(kFALSE);
            lCascadeResult[lN] -> SetCutMCPDGCodeAssociation(kFALSE);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_ITSRefitNPB011_NoAssoc",lParticleName[i].Data() ) );
            lCascadeResult[lN] -> SetCutUseITSRefitPositive(kTRUE);
            lCascadeResult[lN] -> SetCutUseITSRefitBachelor(kTRUE);
            lCascadeResult[lN] -> SetCutMCUseMCProperties(kFALSE);
            lCascadeResult[lN] -> SetCutMCPhysicalPrimary(kFALSE);
            lCascadeResult[lN] -> SetCutMCPDGCodeAssociation(kFALSE);
            lN++;
        }
        
        Float_t lLowRap = -0.6;
        Float_t lHighRap = -0.5;
        for(Int_t i=0;i<4;i++){
            lLowRap = -0.6;
            lHighRap = -0.5;
            for(Int_t irapbin=0;irapbin<12;irapbin++){
                lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%f_%f",lParticleName[i].Data(),"DefaultRapiditySweep",lLowRap,lHighRap ) );
                lCascadeResult[lN]->SetCutMinRapidity(lLowRap);
                lCascadeResult[lN]->SetCutMaxRapidity(lHighRap);
                lN++;
                lLowRap+=0.1;
                lHighRap+=0.1;
            }
        }
        
        // STEP 4: Creation of objects to be used in systematics
        // Optimized via use of copy constructors
        for(Int_t i = 0 ; i < 4 ; i ++){
            for(Int_t iCut = 0 ; iCut < 15 ; iCut ++){
                
                //LOOSE VARIATIONS
                //Create a new object from default
                lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),lCutName[iCut].Data(),lConfName[0].Data()) );
                
                if(iCut ==  0 ) lCascadeResult[lN]->SetCutDCANegToPV            ( lcuts[i][0][iCut] ) ;
                if(iCut ==  1 ) lCascadeResult[lN]->SetCutDCAPosToPV            ( lcuts[i][0][iCut] ) ;
                if(iCut ==  2 ) lCascadeResult[lN]->SetCutDCAV0Daughters        ( lcuts[i][0][iCut] ) ;
                if(iCut ==  3 ) lCascadeResult[lN]->SetCutV0Radius              ( lcuts[i][0][iCut] ) ;
                
                //Setters for Cascade Cuts
                if(iCut ==  4 ) lCascadeResult[lN]->SetCutDCAV0ToPV             ( lcuts[i][0][iCut] ) ;
                if(iCut ==  5 ) lCascadeResult[lN]->SetCutV0Mass                ( lcuts[i][0][iCut] ) ;
                if(iCut ==  6 ) lCascadeResult[lN]->SetCutDCABachToPV           ( lcuts[i][0][iCut] ) ;
                if(iCut ==  7 ){
                    lCascadeResult[lN]->SetCutDCACascDaughters      ( lcuts[i][0][iCut] ) ;
                    lCascadeResult[lN]->SetCutVarDCACascDau ( 1.2 * TMath::Exp(0.0470076), -0.917006, 0, 1, 1.2 * 0.5 );
                }
                if(iCut ==  8 ) lCascadeResult[lN]->SetCutCascRadius            ( lcuts[i][0][iCut] ) ;
                
                //Miscellaneous
                if(iCut ==  9 ) lCascadeResult[lN]->SetCutProperLifetime        ( lcuts[i][0][iCut] ) ;
                if(iCut == 10 ) lCascadeResult[lN]->SetCutMaxV0Lifetime         ( lcuts[i][0][iCut] ) ;
                if(iCut == 11 ) lCascadeResult[lN]->SetCutMinTrackLength        ( lcuts[i][0][iCut] ) ;
                if(iCut == 12 ) lCascadeResult[lN]->SetCutTPCdEdx               ( lcuts[i][0][iCut] ) ;
                if(iCut == 13 ) lCascadeResult[lN]->SetCutXiRejection           ( lcuts[i][0][iCut] ) ;
                if(iCut == 14 ) lCascadeResult[lN]->SetCutDCACascadeToPV        ( lcuts[i][0][iCut] ) ;
                
                //Print this variation, add to pool
                //lCascadeResult[lN]->Print();
                lN++;
                
                //TIGHT VARIATIONS
                //Create a new object from default
                lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),lCutName[iCut].Data(),lConfName[2].Data()) );
                
                if(iCut ==  0 ) lCascadeResult[lN]->SetCutDCANegToPV            ( lcuts[i][2][iCut] ) ;
                if(iCut ==  1 ) lCascadeResult[lN]->SetCutDCAPosToPV            ( lcuts[i][2][iCut] ) ;
                if(iCut ==  2 ) lCascadeResult[lN]->SetCutDCAV0Daughters        ( lcuts[i][2][iCut] ) ;
                if(iCut ==  3 ) lCascadeResult[lN]->SetCutV0Radius              ( lcuts[i][2][iCut] ) ;
                
                //Setters for Cascade Cuts
                if(iCut ==  4 ) lCascadeResult[lN]->SetCutDCAV0ToPV             ( lcuts[i][2][iCut] ) ;
                if(iCut ==  5 ) lCascadeResult[lN]->SetCutV0Mass                ( lcuts[i][2][iCut] ) ;
                if(iCut ==  6 ) lCascadeResult[lN]->SetCutDCABachToPV           ( lcuts[i][2][iCut] ) ;
                if(iCut ==  7 ){
                    lCascadeResult[lN]->SetCutDCACascDaughters      ( lcuts[i][2][iCut] ) ;
                    lCascadeResult[lN]->SetCutVarDCACascDau ( 0.8 * TMath::Exp(0.0470076), -0.917006, 0, 1, 0.8 * 0.5 );
                }
                if(iCut ==  8 ) lCascadeResult[lN]->SetCutCascRadius            ( lcuts[i][2][iCut] ) ;
                
                //Miscellaneous
                if(iCut ==  9 ) lCascadeResult[lN]->SetCutProperLifetime        ( lcuts[i][2][iCut] ) ;
                if(iCut == 10 ) lCascadeResult[lN]->SetCutMaxV0Lifetime         ( lcuts[i][2][iCut] ) ;
                if(iCut == 11 ) lCascadeResult[lN]->SetCutMinTrackLength        ( lcuts[i][2][iCut] ) ;
                if(iCut == 12 ) lCascadeResult[lN]->SetCutTPCdEdx               ( lcuts[i][2][iCut] ) ;
                if(iCut == 13 ) lCascadeResult[lN]->SetCutXiRejection           ( lcuts[i][2][iCut] ) ;
                if(iCut == 14 ) lCascadeResult[lN]->SetCutDCACascadeToPV        ( lcuts[i][2][iCut] ) ;
                
                //Print this variation, add to pool
                //lCascadeResult[lN]->Print();
                lN++;
            }
        }
        
        //STEP 5: re-parametrization of cosines for tight and loose variations (done manually)
        for(Int_t i = 0 ; i < 4 ; i ++){
            //======================================================
            //V0CosPA Variations
            //======================================================
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"V0CosPA","Loose") );
            lCascadeResult[lN]->SetCutVarV0CosPA(TMath::Exp(  -1.77429),
                                                 -0.692453,
                                                 TMath::Exp( -2.01938),
                                                 -0.201574,
                                                 0.0776465);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"V0CosPA","Tight") );
            lCascadeResult[lN]->SetCutVarV0CosPA(TMath::Exp(  -1.21892),
                                                 -41.8521,
                                                 TMath::Exp(   -1.278),
                                                 -0.894064,
                                                 0.0303932);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"V0CosPA","VeryTight") );
            lCascadeResult[lN]->SetCutVarV0CosPA(TMath::Exp(   12.8077),
                                                 -21.2944,
                                                 TMath::Exp( -1.53357),
                                                 -0.920017,
                                                 0.0262315);
            
            lN++;
            //======================================================
            //CascCosPA Variations
            //======================================================
            if( i < 2 ){
                lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"CascCosPA","Loose") );
                lCascadeResult[lN]->SetCutVarCascCosPA(TMath::Exp(  -1.77429),
                                                       -0.692453,
                                                       TMath::Exp( -2.01938),
                                                       -0.201574,
                                                       0.0776465);
                lN++;
                lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"CascCosPA","Tight") );
                lCascadeResult[lN]->SetCutVarCascCosPA(TMath::Exp(   12.8752),
                                                       -21.522,
                                                       TMath::Exp( -1.49906),
                                                       -0.813472,
                                                       0.0480962);
                lN++;
                lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"CascCosPA","VeryTight") );
                lCascadeResult[lN]->SetCutVarCascCosPA(TMath::Exp(    12.801),
                                                       -21.6157,
                                                       TMath::Exp( -1.66297),
                                                       -0.889246,
                                                       0.0346838);
                lN++;
            }
            if( i >= 2 ){
                lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"CascCosPA","Loose") );
                lCascadeResult[lN]->SetCutVarCascCosPA(TMath::Exp(4.86664),
                                                       -10.786,
                                                       TMath::Exp(-1.33411),
                                                       -0.729825,
                                                       0.0695724);
                lN++;
                lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"CascCosPA","Tight") );
                lCascadeResult[lN]->SetCutVarCascCosPA(TMath::Exp(    12.801),
                                                       -21.6157,
                                                       TMath::Exp( -1.66297),
                                                       -0.889246,
                                                       0.0346838);
                lN++;
                lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"CascCosPA","VeryTight") );
                lCascadeResult[lN]->SetCutCascCosPA             ( 0.9992 );
                lN++;
            }
            //======================================================
            //BBCosPA Variations
            //======================================================
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"BBCosPA","VeryLoose") );
            lCascadeResult[lN]->SetCutBachBaryonCosPA        ( 2 );
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"BBCosPA","Loose") );
            lCascadeResult[lN]->SetCutBachBaryonCosPA        ( TMath::Cos(0.03) ) ;
            lCascadeResult[lN]->SetCutVarBBCosPA(TMath::Exp(    -2.8798),
                                                 -20.9876,
                                                 TMath::Exp(  -3.10847),
                                                 -0.73045,
                                                 0.00235147);
            lN++;
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%s",lParticleName[i].Data(),"BBCosPA","Tight") );
            lCascadeResult[lN]->SetCutBachBaryonCosPA        ( TMath::Cos(0.05) ) ;
            lCascadeResult[lN]->SetCutVarBBCosPA(TMath::Exp(   12.4606),
                                                 -20.578,
                                                 TMath::Exp( -2.41442),
                                                 -0.709588,
                                                 0.01079);
            lN++;
        }
        
        //STEP 6: V0 Mass sweep
        //for(Int_t i = 0 ; i < 4 ; i ++){
        //    for(Int_t isweep=0; isweep<20;isweep++){
        //        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_V0MassSweep_%i",lParticleName[i].Data(),isweep) );
        //        lCascadeResult[lN]->SetCutV0MassSigma( ((Double_t)(isweep)/4000.0)); //in GeV/c^2
        //        lN++;
        //    }
        //}
        
        //Float_t lLifetimeCut[4];
        lLifetimeCut[0] = 15.0;
        lLifetimeCut[1] = 15.0;
        lLifetimeCut[2] = 12.0;
        lLifetimeCut[3] = 12.0;
        
        // Float_t lMass[4];
        lMass[0] = 1.322;
        lMass[1] = 1.322;
        lMass[2] = 1.672;
        lMass[3] = 1.672;
        /*
         //Old vertexer-level configuration for cross-checks
         for(Int_t i = 0 ; i < 4 ; i ++){
         //Central result, customized binning: the one to use, usually
         lCascadeResult[lN] = new AliCascadeResult( Form("%s_VertexerLevel",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits,100,lMass[i]-0.050,lMass[i]+0.050);
         
         
         //This is MC: generate profile for G3/F (if ever needed)
         lCascadeResult[lN] -> InitializeProtonProfile();
         
         //Default cuts: use vertexer level ones
         //Setters for V0 Cuts
         lCascadeResult[lN]->SetCutDCANegToPV            ( 0.2 ) ;
         lCascadeResult[lN]->SetCutDCAPosToPV            ( 0.2 ) ;
         lCascadeResult[lN]->SetCutDCAV0Daughters        (  1. ) ;
         lCascadeResult[lN]->SetCutV0CosPA               ( 0.98 ) ;
         lCascadeResult[lN]->SetCutV0Radius              (  3 ) ;
         //Setters for Cascade Cuts
         lCascadeResult[lN]->SetCutDCAV0ToPV             ( 0.1 ) ;
         lCascadeResult[lN]->SetCutV0Mass                ( 0.006 ) ;
         lCascadeResult[lN]->SetCutDCABachToPV           ( 0.03 ) ;
         lCascadeResult[lN]->SetCutDCACascDaughters      ( 1. ) ;
         lCascadeResult[lN]->SetCutCascRadius            ( 1.2 ) ;
         if(i==2||i==3)
         lCascadeResult[lN]->SetCutCascRadius            ( 1.0 ) ; //omega case
         lCascadeResult[lN]->SetCutCascCosPA             ( 0.98 ) ;
         //Miscellaneous
         lCascadeResult[lN]->SetCutProperLifetime        ( lLifetimeCut[i] ) ;
         lCascadeResult[lN]->SetCutMinTrackLength           ( 90.0 ) ;
         lCascadeResult[lN]->SetCutTPCdEdx               ( 4.0 ) ;
         lCascadeResult[lN]->SetCutXiRejection           ( 0.008 ) ;
         lCascadeResult[lN]->SetCutBachBaryonCosPA        ( TMath::Cos(0.006) ) ;
         //Add result to pool
         lN++;
         }
         */
    }
    for (Int_t iconf = 0; iconf<lN; iconf++){
        //      cout<<"["<<iconf<<"/"<<lN<<"] Adding config named "<<lCascadeResult[iconf]->GetName()<<endl;
        AddConfiguration(lCascadeResult[iconf]);
    }
    // cout<<"Added "<<lN<<" Cascade configurations to output."<<endl;
}


//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::AddCascadeConfiguration276TeV()
//Adds 2.76 TeV cascade analysis configuration
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimits[] = {0.4, 0.5, 0.6,
        0.7,0.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
        2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,
        4.4,4.6,4.8,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.,11.,12.};
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimits[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleName[] = {"XiMinus", "XiPlus",  "OmegaMinus", "OmegaPlus"};
    
    //Just a counter and one array, please. Nothing else needed
    AliCascadeResult *lCascadeResult[100];
    Long_t lN = 0;
    
    //Map to mass hypothesis
    AliCascadeResult::EMassHypo lMassHypo[4];
    lMassHypo[0] = AliCascadeResult::kXiMinus;
    lMassHypo[1] = AliCascadeResult::kXiPlus;
    lMassHypo[2] = AliCascadeResult::kOmegaMinus;
    lMassHypo[3] = AliCascadeResult::kOmegaPlus;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        //2.76 Config result, customized binning: the one to use, usually
        lCascadeResult[lN] = new AliCascadeResult( Form("%s_276TeV",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits);
        lCascadeResult[lN] -> InitializeProtonProfile();
        
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( 0.1    ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( 0.1    ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        ( 0.8    ) ;
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.95   ) ; // + variable
        lCascadeResult[lN]->SetCutUse276TeVV0CosPA      ( kTRUE  ) ;
        lCascadeResult[lN]->SetCutV0Radius              ( 3.0    ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( 0.1    ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( 0.005  ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( 0.03   ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( 0.3    ) ;
        lCascadeResult[lN]->SetCutCascRadius            ( 1.5    ) ;
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.9992 ) ;
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( 15.0   ) ;
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( 70     ) ;
        lCascadeResult[lN]->SetCutTPCdEdx               ( 4      ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( 0.008  ) ;
        lCascadeResult[lN]->SetCutDCABachToBaryon       ( 0      ) ;
        
        if(i > 1){
            lCascadeResult[lN]->SetCutCascRadius            ( 1.0 ) ;
            lCascadeResult[lN]->SetCutProperLifetime        ( 8.0 ) ;
        }
        
        //Add result to pool
        lN++;
    }
    
    //Explore restricted rapidity range check
    for(Int_t i = 0 ; i < 4 ; i ++){
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_276TeV_y03",lParticleName[i].Data() ) );
        
        lCascadeResult[lN] -> SetCutMinRapidity(-0.3);
        lCascadeResult[lN] -> SetCutMaxRapidity(+0.3);
        
        //Add result to pool
        lN++;
    }
    
    Float_t lLowRap = -0.6;
    Float_t lHighRap = -0.5;
    for(Int_t i=0;i<4;i++){
        lLowRap = -0.6;
        lHighRap = -0.5;
        for(Int_t irapbin=0;irapbin<12;irapbin++){
            lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i], Form("%s_%s_%f_%f",lParticleName[i].Data(),"276TeVRapiditySweep",lLowRap,lHighRap ) );
            lCascadeResult[lN]->SetCutMinRapidity(lLowRap);
            lCascadeResult[lN]->SetCutMaxRapidity(lHighRap);
            lN++;
            lLowRap+=0.1;
            lHighRap+=0.1;
        }
    }
    
    for (Int_t iconf = 0; iconf<lN; iconf++){
        //      cout<<"["<<iconf<<"/"<<lN<<"] Adding config named "<<lCascadeResult[iconf]->GetName()<<endl;
        AddConfiguration(lCascadeResult[iconf]);
    }
    
    //  cout<<"Added "<<lN<<" cascade configurations to output (corresponding to 2.76 TeV analysis cuts)"<<endl;
}

//________________________________________________________________________
//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::AddCascadeConfigurationPreliminaryCrosscheck()
//Adds cross-checks needed for preliminary validation
{
    // STEP 1: Decide on binning (needed to improve on memory consumption)
    
    // pT binning
    Double_t lPtbinlimits[] = {0.4, 0.5, 0.6,
        0.7,0.8,.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,
        2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,
        4.4,4.6,4.8,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.,11.,12.};
    Long_t lPtbinnumb = sizeof(lPtbinlimits)/sizeof(Double_t) - 1;
    
    // centrality binning
    Double_t lCentbinlimits[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
    Long_t lCentbinnumb = sizeof(lCentbinlimits)/sizeof(Double_t) - 1;
    
    // TStrings for output names
    TString lParticleName[] = {"XiMinus", "XiPlus",  "OmegaMinus", "OmegaPlus"};
    
    //Just a counter and one array, please. Nothing else needed
    AliCascadeResult *lCascadeResult[100];
    Long_t lN = 0;
    
    //Map to mass hypothesis
    AliCascadeResult::EMassHypo lMassHypo[4];
    lMassHypo[0] = AliCascadeResult::kXiMinus;
    lMassHypo[1] = AliCascadeResult::kXiPlus;
    lMassHypo[2] = AliCascadeResult::kOmegaMinus;
    lMassHypo[3] = AliCascadeResult::kOmegaPlus;
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        //Preliminary configuration
        lCascadeResult[lN] = new AliCascadeResult( Form("%s_PXL_Prelim",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits);
        lCascadeResult[lN] -> InitializeProtonProfile();
        
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( 0.2    ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( 0.2    ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        ( 1.0    ) ;
        lCascadeResult[lN]->SetCutV0Radius              ( 3.0    ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( 0.1    ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( 0.005  ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( 0.10   ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( 1.0    ) ;
        lCascadeResult[lN]->SetCutCascRadius            ( 1.2    ) ;
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( 15.0   ) ;
        lCascadeResult[lN]->SetCutMaxV0Lifetime         ( 30.0   ) ;
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( 70     ) ;
        lCascadeResult[lN]->SetCutMinTrackLength        ( -1     ) ;
        lCascadeResult[lN]->SetCutTPCdEdx               ( 4      ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( 0.008  ) ;
        
        if(i > 1){
            lCascadeResult[lN]->SetCutDCACascDaughters      ( 0.7 ) ;
            lCascadeResult[lN]->SetCutCascRadius            ( 1.0 ) ;
            lCascadeResult[lN]->SetCutProperLifetime        (12.0 ) ;
        }
        
        //Parametric angle cut initializations
        //V0 cosine of pointing angle
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarV0CosPA            (TMath::Exp(10.853),
                                                         -25.0322,
                                                         TMath::Exp(-0.843948),
                                                         -0.890794,
                                                         0.057553);
        
        //Cascade cosine of pointing angle
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.95 ) ; //+variable
        if(i < 2){
            lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(4.86664),
                                                             -10.786,
                                                             TMath::Exp(-1.33411),
                                                             -0.729825,
                                                             0.0695724);
        }
        if(i >= 2){
            lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(   12.8752),
                                                             -21.522,
                                                             TMath::Exp( -1.49906),
                                                             -0.813472,
                                                             0.0480962);
        }
        
        //BB cosine of pointing angle
        lCascadeResult[lN]->SetCutBachBaryonCosPA       ( TMath::Cos(0.04) ) ; //+variable
        lCascadeResult[lN]->SetCutVarBBCosPA            (TMath::Exp(-2.29048),
                                                         -20.2016,
                                                         TMath::Exp(-2.9581),
                                                         -0.649153,
                                                         0.00526455);
        
        //Add result to pool
        lN++;
    }
    
    for(Int_t i = 0 ; i < 4 ; i ++){
        //Central result, customized binning: the one to use, usually
        lCascadeResult[lN] = new AliCascadeResult( Form("%s_PXL_Central",lParticleName[i].Data() ),lMassHypo[i],"",lCentbinnumb,lCentbinlimits, lPtbinnumb,lPtbinlimits);
        
        //This is MC: generate profile for G3/F (if ever needed)
        lCascadeResult[lN] -> InitializeProtonProfile();
        
        //Setters for V0 Cuts
        lCascadeResult[lN]->SetCutDCANegToPV            ( 0.20 ) ;
        lCascadeResult[lN]->SetCutDCAPosToPV            ( 0.20 ) ;
        lCascadeResult[lN]->SetCutDCAV0Daughters        ( 1.0  ) ;
        lCascadeResult[lN]->SetCutV0Radius              ( 3.0  ) ;
        //Setters for Cascade Cuts
        lCascadeResult[lN]->SetCutDCAV0ToPV             ( 0.1 ) ;
        lCascadeResult[lN]->SetCutV0Mass                ( 0.005 ) ;
        lCascadeResult[lN]->SetCutDCABachToPV           ( 0.1 ) ;
        lCascadeResult[lN]->SetCutDCACascDaughters      ( 1.0 ) ;
        lCascadeResult[lN]->SetCutVarDCACascDau ( TMath::Exp(0.0470076), -0.917006, 0, 1, 0.5 );
        lCascadeResult[lN]->SetCutCascRadius            ( 1.2 ) ;
        //Miscellaneous
        lCascadeResult[lN]->SetCutProperLifetime        ( 15.0 ) ;
        lCascadeResult[lN]->SetCutMaxV0Lifetime         ( 30 ) ;
        lCascadeResult[lN]->SetCutMinTrackLength        ( 90 ) ;
        lCascadeResult[lN]->SetCutLeastNumberOfClusters( -1 );
        lCascadeResult[lN]->SetCutTPCdEdx               ( 3 ) ;
        lCascadeResult[lN]->SetCutXiRejection           ( 0.008 ) ;
        lCascadeResult[lN]->SetCutDCACascadeToPV        ( 0.8 ) ;
        
        if(i > 1){
            lCascadeResult[lN]->SetCutDCACascDaughters      ( 0.7 ) ;
            lCascadeResult[lN]->SetCutCascRadius            ( 1.0 ) ;
            lCascadeResult[lN]->SetCutProperLifetime        (12.0 ) ;
            lCascadeResult[lN]->SetCutDCACascadeToPV        ( 0.6 ) ;
        }
        //Parametric angle cut initializations
        //V0 cosine of pointing angle
        lCascadeResult[lN]->SetCutV0CosPA               ( 0.95 ) ; //+variable
        lCascadeResult[lN]->SetCutVarV0CosPA            (TMath::Exp(10.853),
                                                         -25.0322,
                                                         TMath::Exp(-0.843948),
                                                         -0.890794,
                                                         0.057553);
        
        //Cascade cosine of pointing angle
        lCascadeResult[lN]->SetCutCascCosPA             ( 0.95 ) ; //+variable
        if(i < 2){
            lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(4.86664),
                                                             -10.786,
                                                             TMath::Exp(-1.33411),
                                                             -0.729825,
                                                             0.0695724);
        }
        if(i >= 2){
            lCascadeResult[lN]->SetCutVarCascCosPA          (TMath::Exp(   12.8752),
                                                             -21.522,
                                                             TMath::Exp( -1.49906),
                                                             -0.813472,
                                                             0.0480962);
        }
        
        //BB cosine of pointing angle
        lCascadeResult[lN]->SetCutBachBaryonCosPA       ( TMath::Cos(0.04) ) ; //+variable
        lCascadeResult[lN]->SetCutVarBBCosPA            (TMath::Exp(-2.29048),
                                                         -20.2016,
                                                         TMath::Exp(-2.9581),
                                                         -0.649153,
                                                         0.00526455);
        
        //Add result to pool
        lN++;
    }
    
    //Add incremental changes. Preliminary -> Latest
    for (Int_t i = 0 ; i<4; i++ ){
        
        // TPC n-sigma BB : 4 -> 3
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i],Form("%s_PXL_Prelim_TPCnsigma",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutTPCdEdx ( 3 );
        lN++; //Add result to pool
        
        // Cluster -> Length
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i],Form("%s_PXL_Prelim_ClusLeng",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( -1 ) ;
        lCascadeResult[lN]->SetCutMinTrackLength        ( 90 ) ;
        lN++; //Add result to pool
        
        // Add 3D Cascade to PV
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i],Form("%s_PXL_Prelim_CascadeDCA3D",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutDCACascadeToPV( 0.8 );
        if( i > 1 ) lCascadeResult[lN]->SetCutDCACascadeToPV( 0.6 );
        lN++; //Add result to pool
        
        // Add Casc. Dau. vs pT
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i],Form("%s_PXL_Prelim_CascadeDCA3D",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutDCACascadeToPV( 1.0 );
        if( i > 1 ) lCascadeResult[lN]->SetCutDCACascadeToPV( 0.6 );
        lCascadeResult[lN]->SetCutVarDCACascDau ( TMath::Exp(0.0470076), -0.917006, 0, 1, 0.5 );
        lN++;
        
        // All changes
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[i],Form("%s_PXL_Prelim_AllChanges",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutTPCdEdx ( 3 );
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( -1 ) ;
        lCascadeResult[lN]->SetCutMinTrackLength        ( 90 ) ;
        lCascadeResult[lN]->SetCutDCACascadeToPV( 1.0 );
        if( i > 1 ) lCascadeResult[lN]->SetCutDCACascadeToPV( 0.6 );
        lCascadeResult[lN]->SetCutVarDCACascDau ( TMath::Exp(0.0470076), -0.917006, 0, 1, 0.5 );
        lN++;
    }
    
    
    //Add incremental changes: Latest -> Preliminary
    for (Int_t i = 0 ; i<4; i++ ){
        
        // TPC n-sigma BB : 3 -> 4
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[4+i],Form("%s_PXL_Central_TPCnsigma",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutTPCdEdx ( 4 );
        lN++; //Add result to pool
        
        // Length -> Cluster
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[4+i],Form("%s_PXL_Central_ClusLeng",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( 70 ) ;
        lCascadeResult[lN]->SetCutMinTrackLength        ( -1 ) ;
        lN++; //Add result to pool
        
        // Remove 3D Cascade to PV
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[4+i],Form("%s_PXL_Central_CascadeDCA3D",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutDCACascadeToPV( 1e3 );
        lN++; //Add result to pool
        
        // Remove Casc. Dau. vs pT
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[4+i],Form("%s_PXL_Central_CascadeDCA3D",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutDCACascadeToPV( 1.0 );
        if( i > 1 ) lCascadeResult[lN]->SetCutDCACascadeToPV( 0.7 );
        lCascadeResult[lN]->SetCutUseVarDCACascDau ( kFALSE );
        lN++;
        
        // All changes
        lCascadeResult[lN] = new AliCascadeResult( lCascadeResult[4+i],Form("%s_PXL_Central_AllChanges",lParticleName[i].Data() ) );
        lCascadeResult[lN]->SetCutTPCdEdx ( 4 );
        lCascadeResult[lN]->SetCutLeastNumberOfClusters ( 70 ) ;
        lCascadeResult[lN]->SetCutMinTrackLength        ( -1 ) ;
        lCascadeResult[lN]->SetCutDCACascadeToPV( 1e3 );
        lCascadeResult[lN]->SetCutDCACascadeToPV( 1.0 );
        if( i > 1 ) lCascadeResult[lN]->SetCutDCACascadeToPV( 0.7 );
        lCascadeResult[lN]->SetCutUseVarDCACascDau ( kFALSE );
        lN++;
    }
    
    
    for (Int_t iconf = 0; iconf<lN; iconf++){
        //  cout<<"["<<iconf<<"/"<<lN<<"] Adding config named "<<lCascadeResult[iconf]->GetName()<<endl;
        AddConfiguration(lCascadeResult[iconf]);
    }
    
    //cout<<"Added "<<lN<<" cascade configurations to output (corresponding to cross-checks to preliminary version)"<<endl;
}


//________________________________________________________________________
Float_t AliAnalysisTaskStrangeCascadesDiscrete::GetDCAz(AliESDtrack *lTrack)
//Encapsulation of DCAz calculation
{
    Float_t b[2];
    Float_t bCov[3];
    lTrack->GetImpactParameters(b,bCov);
    if (bCov[0]<=0 || bCov[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal to zero!");
        bCov[0]=0; bCov[2]=0;
    }
    //Float_t dcaToVertexXY = b[0];
    Float_t dcaToVertexZ = b[1];
    
    return dcaToVertexZ;
}


//________________________________________________________________________
Float_t AliAnalysisTaskStrangeCascadesDiscrete::GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent)
//Encapsulation of CosPA calculation (warning: considers AliESDtrack clones)
{
    Float_t lCosPA = -1;
    
    //Get Magnetic field and primary vertex
    Double_t b=lEvent->GetMagneticField();
    const AliESDVertex *vtxT3D=lEvent->GetPrimaryVertex();
    Double_t xPrimaryVertex=vtxT3D->GetX();
    Double_t yPrimaryVertex=vtxT3D->GetY();
    Double_t zPrimaryVertex=vtxT3D->GetZ();
    
    //Copy AliExternalParam for handling
    AliExternalTrackParam nt(*lNegTrack), pt(*lPosTrack), *lNegClone=&nt, *lPosClone=&pt;
    
    //Find DCA
    Double_t xn, xp;
    dca=lNegClone->GetDCA(lPosClone,b,xn,xp);
    
    //Propagate to it
    nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);
    
    //Create V0 object to do propagation
    AliESDv0 vertex(nt,1,pt,2); //Never mind indices, won't use
    
    //Get CosPA
    lCosPA=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
    
    //Return value
    return lCosPA;
}

//________________________________________________________________________
void AliAnalysisTaskStrangeCascadesDiscrete::CheckChargeV0(AliESDv0 *v0)
{
    // This function checks charge of negative and positive daughter tracks.
    // If incorrectly defined (onfly vertexer), swaps out.
    if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
        //V0 daughter track swapping is required! Note: everything is swapped here... P->N, N->P
        Long_t lCorrectNidx = v0->GetPindex();
        Long_t lCorrectPidx = v0->GetNindex();
        Double32_t    lCorrectNmom[3];
        Double32_t    lCorrectPmom[3];
        v0->GetPPxPyPz( lCorrectNmom[0], lCorrectNmom[1], lCorrectNmom[2] );
        v0->GetNPxPyPz( lCorrectPmom[0], lCorrectPmom[1], lCorrectPmom[2] );
        
        AliExternalTrackParam    lCorrectParamN(
                                                v0->GetParamP()->GetX() ,
                                                v0->GetParamP()->GetAlpha() ,
                                                v0->GetParamP()->GetParameter() ,
                                                v0->GetParamP()->GetCovariance()
                                                );
        AliExternalTrackParam    lCorrectParamP(
                                                v0->GetParamN()->GetX() ,
                                                v0->GetParamN()->GetAlpha() ,
                                                v0->GetParamN()->GetParameter() ,
                                                v0->GetParamN()->GetCovariance()
                                                );
        lCorrectParamN.SetMostProbablePt( v0->GetParamP()->GetMostProbablePt() );
        lCorrectParamP.SetMostProbablePt( v0->GetParamN()->GetMostProbablePt() );
        
        //Get Variables___________________________________________________
        Double_t lDcaV0Daughters = v0 -> GetDcaV0Daughters();
        Double_t lCosPALocal     = v0 -> GetV0CosineOfPointingAngle();
        Bool_t lOnFlyStatusLocal = v0 -> GetOnFlyStatus();
        
        //Create Replacement Object_______________________________________
        AliESDv0 *v0correct = new AliESDv0(lCorrectParamN,lCorrectNidx,lCorrectParamP,lCorrectPidx);
        v0correct->SetDcaV0Daughters          ( lDcaV0Daughters   );
        v0correct->SetV0CosineOfPointingAngle ( lCosPALocal       );
        v0correct->ChangeMassHypothesis       ( kK0Short          );
        v0correct->SetOnFlyStatus             ( lOnFlyStatusLocal );
        
        //Reverse Cluster info..._________________________________________
        v0correct->SetClusters( v0->GetClusters( 1 ), v0->GetClusters ( 0 ) );
        
        *v0 = *v0correct;
        //Proper cleanup..._______________________________________________
        v0correct->Delete();
        v0correct = 0x0;
        
        //Just another cross-check and output_____________________________
        if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
            AliWarning("Found Swapped Charges, tried to correct but something FAILED!");
        } else {
            //AliWarning("Found Swapped Charges and fixed.");
        }
        //________________________________________________________________
    } else {
        //Don't touch it! ---
        //Printf("Ah, nice. Charges are already ordered...");
    }
    return;
}
//---------------------------extra implemented------------------------------------------------------------------
//-------Here: Calculate the cosine of angle from directions of (1) proton in rest frame of Lambda and (2) direction of Lambda
// in the rest frame of Xi (Omega)
Double_t AliAnalysisTaskStrangeCascadesDiscrete::GetCosOfProtonLambdaRestOmegaRest(Double_t pXi[3],Double_t mXi,
                                                                                        Double_t pL[3], Double_t mL,
                                                                                        Double_t pp[3], Double_t mp)

{
    Double_t eXi(0), eL(0), ep(0);
    Double_t betaXi[3], betaLambda[3];
    Double_t L_mom_inXirest[3], p_mom_inLrest[3];
    
    //Determine the Energies
    eXi = TMath::Sqrt(mXi*mXi +pXi[0]*pXi[0] +pXi[1]*pXi[1] +pXi[2]*pXi[2]);
    if(!eXi)AliWarning("---------The energy of XI could not be calculated---------");
    eL =TMath::Sqrt(mL*mL+ pL[0]*pL[0] +pL[1]*pL[1] +pL[2]*pL[2]);
    if(!eL)AliWarning("---------The energy of Lambda could not be calculated---------");
    ep =TMath::Sqrt(mp*mp +pp[0]*pp[0] +pp[1]*pp[1] +pp[2]*pp[2]);
    if(!ep)AliWarning("---------The energy of proton could not be calculated---------");
    
    //Define the lorentz vectors
    TLorentzVector* PXi_lab = new TLorentzVector(pXi[0],pXi[1],pXi[2],eXi);
    TLorentzVector* PLambda_lab = new TLorentzVector(pL[0],pL[1],pL[2],eL);
    TLorentzVector* Pproton_lab = new TLorentzVector(pp[0],pp[1],pp[2],ep);
    
    //   TLorentzVector* PLambda_inXi_restframe = new TLorentzVector();
    //   TLorentzVector* Pproton_inLambda_restframe = new TLorentzVector();
    //defined for the checks of boosts
    //  TLorentzVector* Xi_Check = new TLorentzVector();
    //  TLorentzVector* L_Check = new TLorentzVector();
    
    if (!PLambda_lab || !Pproton_lab || !PXi_lab) AliWarning("LorentzVectors could not be defined");
    
    
    //Define the boost parameters; boost depends on velocities beta_i = p_i/E
    for (int i = 0;  i<3; i++)
    {
        //mind the minus sign! Because the boost has to bring the momentum back (into rest frame)
        betaXi[i] = - pXi[i]/eXi;
        betaLambda[i] = - pL[i]/eL;
    }
    
    //Boost the Lambda lab 4 vector into the rest frame of the Xi (e.g. Omega)
    PLambda_lab -> Boost(betaXi[0], betaXi[1], betaXi[2]);
    if(!PLambda_lab)AliWarning("The 4 momentum of Lambda in Xi rest frame could not be determined");
    //Check that Xi goes into its restframe:
    PXi_lab -> Boost(betaXi[0], betaXi[1], betaXi[2]);
    //get the array of momentum
    L_mom_inXirest[0] = PLambda_lab -> Px();
    L_mom_inXirest[1] = PLambda_lab -> Py();
    L_mom_inXirest[2] = PLambda_lab -> Pz();
    
    //boost the vector of Lambda back again into the Lab frame!
    PLambda_lab -> Boost(-betaXi[0], -betaXi[1], -betaXi[2]);
    
    //  if (!L_mom_inXirest)AliWarning("The Double_t array of Lambda momentum could not be got");
    
    
    //Repeat the same for the proton, i.e. boost the lab proton into the restframe of Lambda
    Pproton_lab->Boost(betaLambda[0],betaLambda[1],betaLambda[2]);
    if(!Pproton_lab)AliWarning("The 4 momentum of Proton in Lambda rest frame could not be determined");
    //Check that Lambda goes into its rest frame
    PLambda_lab->Boost(betaLambda[0],betaLambda[1],betaLambda[2]);
    
    //get the array
    p_mom_inLrest[0] = Pproton_lab->Px();
    p_mom_inLrest[1] = Pproton_lab->Py();
    p_mom_inLrest[2] = Pproton_lab->Pz();
    
    // if(!p_mom_inLrest)AliWarning("The Double_t array of Proton momentum could not be got");
    
    //Calculate the cosine of the angle
    Double_t p_tot(0), L_tot(0), vProd(0);
    Double_t cos(-2);
    
    p_tot = TMath::Sqrt(p_mom_inLrest[0]*p_mom_inLrest[0]+p_mom_inLrest[1]*p_mom_inLrest[1]+p_mom_inLrest[2]*p_mom_inLrest[2]);
    L_tot = TMath::Sqrt(L_mom_inXirest[0]*L_mom_inXirest[0]+L_mom_inXirest[1]*L_mom_inXirest[1]+L_mom_inXirest[2]*L_mom_inXirest[2]);
    if((p_tot == 0)||(L_tot==0))AliWarning("Check the total momentum! Should be bigger than 0");
    
    vProd = p_mom_inLrest[0]*L_mom_inXirest[0] + p_mom_inLrest[1]*L_mom_inXirest[1] + p_mom_inLrest[2]*L_mom_inXirest[2];
    cos = vProd/(p_tot*L_tot);
    
    
    //  TLorentzVector* v = new TLorentzVector(3.,0.,0.,5.);
    //  v -> Boost(-0.6,0.,0.);
    //  cout << "v -> E() = " <<  v -> E() << endl;
    
    return cos;
    
}

Double_t AliAnalysisTaskStrangeCascadesDiscrete::DzetaFromMomenta(Double_t p1[3], Double_t p2[3], Double_t pp[3])
{
    //Calculate the vector product of two particles in question, e.g. Omega and Lambda(daughter)
    //This will create an orthogonal vector to the "reaction plane" and is an axial vector (does not change under boosts within the reaction plane)
    Double_t vector_product[3] = {0.,0.,0.,};
    vector_product[0] = p1[1]*p2[2] - p1[2]*p2[1];
    vector_product[1] = p1[2]*p2[0] - p1[0]*p2[2];
    vector_product[2] = p1[0]*p2[1] - p1[1]*p2[0];
    
    Double_t dzeta = 0.;
    dzeta = vector_product[0]*pp[0] + vector_product[1]*pp[1] + vector_product[2]*pp[2];
    if(!dzeta)AliWarning("Dzeta parameter could not be calculated, something wrong with vectors of the chosen particles");
    return dzeta;
    
}
//-------------------------------------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskStrangeCascadesDiscrete::GoodTrack(AliESDtrack* trackESD)
{
    //check the TPC refit.
    ULong_t status = trackESD->GetStatus();
    if((status & AliVTrack::kTPCrefit) == 0) return kFALSE;
    
    //check the ITS status
    Int_t statusITS1 = GetITSstatus(trackESD, 0);
    Int_t statusITS2 = GetITSstatus(trackESD, 1);
    if((statusITS1!=1)&&(statusITS2!=1)) return kFALSE; //-> the hit is in the first OR second ITS layer
    
    //check the ITS innermost layers
    Bool_t lits1 = trackESD -> HasPointOnITSLayer(0);
    Bool_t lits2 = trackESD -> HasPointOnITSLayer(1);
    if ((lits1==0)&&(lits2==0)) return kFALSE;
    
    //TPC clusters
    Int_t tpcnclus=trackESD->GetTPCNcls();
    if(tpcnclus<70) return kFALSE;
    
    //maximal chi2 per TPC cluster should be smaller than 4.
    Float_t tpcchi2 = 1000;
    if(tpcnclus!=0) tpcchi2=(Float_t)trackESD->GetTPCchi2()/tpcnclus;
    if(tpcchi2 > 4) return kFALSE;
    
    //rapidity
    Double_t y = trackESD->Y();
    if (!y || TMath::Abs(y) > 0.5) return kFALSE;
    
    //pseudorapity
    Double_t eta = trackESD->Eta();
    if (!eta || TMath::Abs(eta) > 0.8) return kFALSE;
    
    //kink index
    if(trackESD->GetKinkIndex(0)>0 || trackESD -> GetKinkIndex(1) > 0) return kFALSE;
    
    
    return kTRUE;
    
    
}
//-------------------------------------------------------------------------------------------------------------------------------

Int_t AliAnalysisTaskStrangeCascadesDiscrete::GetITSstatus(const AliVTrack * const track, Int_t layer) const {
    //
    // Check ITS layer status
    //
    Int_t status = 0;
    if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
        Int_t det;
        Float_t xloc, zloc;
        const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
        if(esdtrack) esdtrack->GetITSModuleIndexInfo(layer, det, status, xloc, zloc);
    }
    return status;
}

//-------------------------------------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskStrangeCascadesDiscrete::CheckTOFstatus(AliESDtrack* trackESD){
    
    UInt_t TOFstatus = AliVTrack::kTOFpid;
    if ((trackESD->GetStatus() && TOFstatus) != TOFstatus) return kFALSE;
    return kTRUE;
}






//--------------------------------------------------------------------------------------------------------------------------------
//get a polar angle of the daughter in respect to the world coordinate system
//in a rest frame of the mother
/*Double_t AliAnalysisTaskStrangeCascadesDiscrete::GetThetaOfDaughterInRestOfMotherWorld(Double_t pMother[3],Double_t mMother,
 Double_t pDaughter[3], Double_t mDaughter)
 {
 Double_t eMother(0), eDaughter(0);
 eDaughter = TMath::Sqrt(mDaughter*mDaughter +pDaughter[0]*pDaughter[0] +pDaughter[1]*pDaughter[1] +pDaughter[2]*pDaughter[2]);
 if(!eDaughter)AliWarning("---------The energy of the daughter could not be calculated---------");
 eMother =TMath::Sqrt(mMother*mMother+ pMother[0]*pMother[0] +pMother[1]*pMother[1] +pMother[2]*pMother[2]);
 if(!eMother)AliWarning("---------The energy of the mother could not be calculated---------");
 
 TLorentzVector* PMother = new TLorentzVector(pMother[0],pMother[1],pMother[2],eMother);
 TLorentzVector* PDaughter = new TLorentzVector(pDaughter[0],pDaughter[1],pDaughter[2],eDaughter);
 
 Double_t betaMother[3];
 for (int i = 0;  i<3; i++)
 {
 //mind the minus sign! Because the boost has to bring the momentum back (into rest frame)
 betaMother[i] = - pMother[i]/eMother;
 }
 PDaughter -> Boost(betaMother[0], betaMother[1], betaMother[2]);
 PMother -> Boost(betaMother[0], betaMother[1], betaMother[2]);
 
 Double_t output[3];
 output[0] = PDaughter -> Phi();
 if(TMath::Abs(PMother -> E() - mMother)>1e-4){
 AliWarning("The masses in the rest frame do not match; check the boost!");
 output[1] = -1000;
 output[2] = -1000;
 }
 else{
 output[1] = PDaughter -> Theta();
 output[2] = PDaughter -> CosTheta();
 }
 return output;
 
 }*/
//--------------------------------------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskStrangeCascadesDiscrete::EventTriggerWord(AliESDEvent* lESDevent){
    ULong64_t triggerword = lESDevent->GetTriggerMask();
    //UChar_t triggerword = lESDevent->GetTriggerCluster();
    // if (!lESDevent) break;
    Int_t r(-10.);
    if (triggerword == AliVEvent::kINT7) r = 1;
    else if (triggerword == AliVEvent::kMB) r = 0;
    else if (triggerword == AliVEvent::kMUON) r = 2;
    else if (triggerword == AliVEvent::kHighMult) r = 3;
    else if (triggerword == AliVEvent::kEMC1) r = 4;
    else if (triggerword == AliVEvent::kCINT5) r = 5;
    else if (triggerword == AliVEvent::kCMUS5) r = 6;
    else if (triggerword == AliVEvent::kMUSPB) r = 7;
    else if (triggerword == AliVEvent::kMUSH7) r = 8;
    else if (triggerword == AliVEvent::kMUSHPB) r = 9;
    else if (triggerword == AliVEvent::kMUL7) r = 10;
    else if (triggerword == AliVEvent::kMuonLikePB) r = 11;
    else if (triggerword == AliVEvent::kEMC7) r = 12;
    else if (triggerword == AliVEvent::kEMC8) r = 13;
    else if (triggerword == AliVEvent::kMUS7) r = 14;
    else if (triggerword == AliVEvent::kPHI1) r = 15;
    else if (triggerword == AliVEvent::kPHI7) r = 16;
    else if (triggerword == AliVEvent::kPHI8) r = 17;
    else if (triggerword == AliVEvent::kPHOSPb) r = 18;
    else if (triggerword == AliVEvent::kEMCEJE) r = 19;
    else if (triggerword == AliVEvent::kEMCEGA) r = 20;
    
    else if (triggerword == AliVEvent::kCentral) r = 21;
    else if (triggerword == AliVEvent::kSemiCentral) r = 22;
    else if (triggerword == AliVEvent::kDG5) r = 23;
    else if (triggerword == AliVEvent::kZED) r = 24;
    else if (triggerword == AliVEvent::kSPI7) r = 25;
    
    else if (triggerword == AliVEvent::kSPI) r = 26;
    else if (triggerword == AliVEvent::kINT8) r = 27;
    else if (triggerword == AliVEvent::kMuonSingleLowPt8) r = 28;
    else if (triggerword == AliVEvent::kMuonSingleHighPt8) r = 29;
    
    else if (triggerword == AliVEvent::kMuonLikeLowPt8) r = 30;
    else if (triggerword == AliVEvent::kMuonUnlikeLowPt8) r = 31;
    else if (triggerword == AliVEvent::kMuonUnlikeLowPt0) r = 32;
    
    else if (triggerword == AliVEvent::kUserDefined) r = 33;
    else if (triggerword == AliVEvent::kTRD) r = 34;
    else if (triggerword == AliVEvent::kFastOnly) r = 35;
    //   else if (triggerword == AliVEvent::kAny) r = 36;
    //   else if (triggerword == AliVEvent::kAnyINT) r = 37;
    
    else r = -1;
    
    return r;
}


