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

/* AliAnalysisTaskStrangeCascadesTriggerAODRun2
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TList.h"
#include "TTree.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultEstimator.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliPIDResponse.h"

#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDUtils.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include "AliAnalysisTaskStrangeCascadesTriggerAODRun2.h"

class AliAnalysisTaskStrangeCascadesTriggerAODRun2;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskStrangeCascadesTriggerAODRun2) // classimp: necessary for root

AliAnalysisTaskStrangeCascadesTriggerAODRun2::AliAnalysisTaskStrangeCascadesTriggerAODRun2() : 
AliAnalysisTaskSE(), fOutputList(0), 
fTreeCascade(0), fTreeV0(0), fTreeRsn(0), fTreeRsnBkg(0), fTreePrimTrack(0), fDummyTree(0),
fPIDResponse(0), fUtils(0), 
fRsnEvent(0), fRsnMiniEvent(0), fResonanceFinders(0), fRsnTrackCuts(0), fRsnPairCuts(0),
fPairName(""), fMotherMass(0),

//---> Flags controlling TTree outputs
fSaveV0s(kTRUE), 
fSaveRsn(kTRUE), 
fSavePrimaries(kTRUE), 
fComputationType(""),

//---> Variables controlling PV selections
fkMaxPVR2D(10.), 
fkMaxPVZ(10.), 

//---> Variables controlling cascade and V0 default selections
fMinNbrCrossedRows(0), 
fMinPtToSave(0), 
fMaxPtToSave(20), 
fMaxAbsEta(0.8), 
fMaxAbsRap(0.5),
fRejectCascKink(kFALSE),
fRejectV0Kink(kFALSE),

//---> Flags controlling cascade and V0 custom selections
fCascSaveAddConfig(kFALSE), 
fV0SaveAddConfig(kFALSE),
fPrimariesSaveAddConfig(kFALSE),
fCascadeResult(0), 
fV0Result(0),
fPrimTrackCuts(0),

//---> Variables controlling event mixing
fNMix(0), 
fMaxDiffVz(0), 
fMaxDiffAngle(0), 
fMaxDiffMult(0),

//===========================================================================================
//   Variables for Cascade Tree
//===========================================================================================
//-----------BASIC-INFO---------------------------
fTreeCascVarCharge(0),
fTreeCascVarMassAsXi(0),
fTreeCascVarMassAsOmega(0),
fTreeCascVarPtot(0),
fTreeCascVarV0Ptot(0),
fTreeCascVarPt(0),
fTreeCascVarV0Pt(0),
fTreeCascVarRapXi(0),
fTreeCascVarRapOmega(0),
fTreeCascVarBachEta(0),
fTreeCascVarPosEta(0),
fTreeCascVarNegEta(0),
fTreeCascVarPhi(0),
fTreeCascVarTheta(0),
//-----------INFO-FOR-CUTS------------------------
fTreeCascVarAlpha(0),
fTreeCascVarPtArm(0),
fTreeCascVarAlphaV0(0),
fTreeCascVarPtArmV0(0),
fTreeCascVarDCACascDau(0),
fTreeCascVarDCABachToPV(0),
fTreeCascVarDCAV0Dau(0),
fTreeCascVarDCAV0ToPV(0),
fTreeCascVarDCAPosToPV(0),
fTreeCascVarDCANegToPV(0),
fTreeCascVarCascCosPA(0),
fTreeCascVarCascCosPASpecial(0),
    
fTreeCascVarCascRadius(0),
fTreeCascVarV0Mass(0),
fTreeCascVarV0MassAsLambda(0),
fTreeCascVarV0MassAsAntiLambda(0),
fTreeCascVarV0Radius(0),
fTreeCascVarV0CosPA(0),
fTreeCascVarWrongCosPA(0),
fTreeCascVarDCABachToBaryon(0),
fTreeCascVarLeastNbrCrossedRows(0),
fTreeCascVarLeastNbrClusters(0),
fTreeCascVarLeastRatioCrossedRowsOverFindable(0),
fTreeCascVarNbrCrossedRowsOverLength(0),
fTreeCascVarMaxChi2PerCluster(0),
fTreeCascVarMinTrackLength(0),
//-----------MULTIPLICITY-INFO--------------------
fTreeCascVarCentrality(0),          //! VZERO Percentile             - AliMultSelection
fTreeCascVarVZEROMultSel(0),        //! VZERO Multiplicity           - AliMultSelection
fTreeCascVarVZEROMultSig(0),        //! VZERO Multiplicity           - VZER0 signal
fTreeCascVarVZEROMultSigCorr(0),    //! Corrected VZERO Multiplicity - VZER0 signal
fTreeCascVarSPDMult(0),             //! SPD Multiplicity             - SPD Tracklets
fTreeCascVar_TriggerMask(0),
fTreeCascVarIsIncompleteDAQ(0),
fTreeCascVarSPDPileupFlag(0),
fTreeCascVarMVPileupFlag(0),
fTreeCascVarOOBPileupFlag(0),
//-----------DECAY-LENGTH-INFO--------------------
fTreeCascVarDistOverTotMom(0),
fTreeCascVarV0DistOverTotMom(0),
//------------------------------------------------
//---------------PID-TPC-INFO---------------------
fTreeCascVarBachNSigmaPion(0),
fTreeCascVarBachNSigmaKaon(0),
fTreeCascVarPosNSigmaProton(0),
fTreeCascVarPosNSigmaPion(0),
fTreeCascVarNegNSigmaProton(0),
fTreeCascVarNegNSigmaPion(0),
//---------------PID-TOF-INFO---------------------
fTreeCascVarBachTOFNSigmaPion(0),
fTreeCascVarBachTOFNSigmaKaon(0),
fTreeCascVarPosTOFNSigmaProton(0),
fTreeCascVarPosTOFNSigmaPion(0),
fTreeCascVarNegTOFNSigmaProton(0),
fTreeCascVarNegTOFNSigmaPion(0),
//---------------PID-ITS-INFO---------------------
fTreeCascVarBachITSNSigmaPion(0),
fTreeCascVarBachITSNSigmaKaon(0),
fTreeCascVarPosITSNSigmaProton(0),
fTreeCascVarPosITSNSigmaPion(0),
fTreeCascVarNegITSNSigmaProton(0),
fTreeCascVarNegITSNSigmaPion(0),
//---Raw TPC dEdx + PIDForTracking information----
fTreeCascVarPosdEdx(0),
fTreeCascVarNegdEdx(0),
fTreeCascVarBachdEdx(0),
fTreeCascVarPosPIDForTracking(0),
fTreeCascVarNegPIDForTracking(0),
fTreeCascVarBachPIDForTracking(0),
//------------------------------------------------
fTreeCascVarChi2Cascade(0),
fTreeCascVarChi2CascadePerNDF(0),
fTreeCascVarChi2V0(0),
//------------------------------------------------
fTreeCascVarBachTrackStatus(0),
fTreeCascVarPosTrackStatus(0),
fTreeCascVarNegTrackStatus(0),
//------------------------------------------------
fTreeCascVarBachDCAz(0),
fTreeCascVarPosDCAz(0),
fTreeCascVarNegDCAz(0),
//------------FULL-MOMENTUM-INFO------------------
fTreeCascVarBachPx(0),
fTreeCascVarBachPy(0),
fTreeCascVarBachPz(0),
fTreeCascVarPosPx(0),
fTreeCascVarPosPy(0),
fTreeCascVarPosPz(0),
fTreeCascVarNegPx(0),
fTreeCascVarNegPy(0),
fTreeCascVarNegPz(0),
//------------------------------------------------
fTreeCascVarPrimVertexX(0),
fTreeCascVarPrimVertexY(0),
fTreeCascVarPrimVertexZ(0),
fTreeCascVarCascadeDecayX(0),
fTreeCascVarCascadeDecayY(0),
fTreeCascVarCascadeDecayZ(0),
fTreeCascVarV0DecayX(0),
fTreeCascVarV0DecayY(0),
fTreeCascVarV0DecayZ(0),
//------------------------------------------------
fTreeCascVarBachIndex(0),
fTreeCascVarPosIndex(0),
fTreeCascVarNegIndex(0),
//---------------CLUSTER-INFO---------------------
fTreeCascVarBachITSClusters0(kFALSE),
fTreeCascVarBachITSClusters1(kFALSE),
fTreeCascVarBachITSClusters2(kFALSE),
fTreeCascVarBachITSClusters3(kFALSE),
fTreeCascVarBachITSClusters4(kFALSE),
fTreeCascVarBachITSClusters5(kFALSE),
    
fTreeCascVarPosITSClusters0(kFALSE),
fTreeCascVarPosITSClusters1(kFALSE),
fTreeCascVarPosITSClusters2(kFALSE),
fTreeCascVarPosITSClusters3(kFALSE),
fTreeCascVarPosITSClusters4(kFALSE),
fTreeCascVarPosITSClusters5(kFALSE),
    
fTreeCascVarNegITSClusters0(kFALSE),
fTreeCascVarNegITSClusters1(kFALSE),
fTreeCascVarNegITSClusters2(kFALSE),
fTreeCascVarNegITSClusters3(kFALSE),
fTreeCascVarNegITSClusters4(kFALSE),
fTreeCascVarNegITSClusters5(kFALSE),
    
//------------------------------------------------
fTreeCascVarPosITSSharedClusters0(kFALSE),
fTreeCascVarPosITSSharedClusters1(kFALSE),
fTreeCascVarPosITSSharedClusters2(kFALSE),
fTreeCascVarPosITSSharedClusters3(kFALSE),
fTreeCascVarPosITSSharedClusters4(kFALSE),
fTreeCascVarPosITSSharedClusters5(kFALSE),

fTreeCascVarNegITSSharedClusters0(kFALSE),
fTreeCascVarNegITSSharedClusters1(kFALSE),
fTreeCascVarNegITSSharedClusters2(kFALSE),
fTreeCascVarNegITSSharedClusters3(kFALSE),
fTreeCascVarNegITSSharedClusters4(kFALSE),
fTreeCascVarNegITSSharedClusters5(kFALSE),

fTreeCascVarBachITSSharedClusters0(kFALSE),
fTreeCascVarBachITSSharedClusters1(kFALSE),
fTreeCascVarBachITSSharedClusters2(kFALSE),
fTreeCascVarBachITSSharedClusters3(kFALSE),
fTreeCascVarBachITSSharedClusters4(kFALSE),
fTreeCascVarBachITSSharedClusters5(kFALSE),

//---------------OOB-PILEUP-INFO---------------------
fTreeCascVarBachTOFExpTDiff(0),
fTreeCascVarPosTOFExpTDiff(0),
fTreeCascVarNegTOFExpTDiff(0),

fTreeCascVarBachTOFSignal(0),
fTreeCascVarPosTOFSignal(0),
fTreeCascVarNegTOFSignal(0),
    
fTreeCascVarBachTOFBCid(0),
fTreeCascVarPosTOFBCid(0),
fTreeCascVarNegTOFBCid(0),

//Event info
fTreeCascVarAmplitudeV0A(0),
fTreeCascVarAmplitudeV0C(0),
fTreeCascVarClosestNonEmptyBC(0),
    
//Kink tagging
fTreeCascVarBachIsKink(kFALSE),
fTreeCascVarPosIsKink(kFALSE),
fTreeCascVarNegIsKink(kFALSE),
    
//Cowboy/sailor studies
fTreeCascVarIsCowboy(kFALSE),   //store if V0 is cowboy-like or sailor-like in XY plane
fTreeCascVarCowboyness(0), //negative -> cowboy, positive -> sailor
fTreeCascVarIsCascadeCowboy(kFALSE),   //store if V0 is cowboy-like or sailor-like in XY plane
fTreeCascVarCascadeCowboyness(0), //negative -> cowboy, positive -> sailor
    
fTreeCascVarMagField(0),
fTreeCascVarRunNumber(0),
fTreeCascVarEventNumber(0),
//===========================================================================================
//   Variables for V0 Tree
//===========================================================================================
//-----------BASIC-INFO---------------------------
fTreeV0VarMassAsK0s(0), 
fTreeV0VarMassAsLambda(0), 
fTreeV0VarMassAsAntiLambda(0), 
fTreeV0VarRapK0Short(0), 
fTreeV0VarRapLambda(0), 
fTreeV0VarPosEta(0), 
fTreeV0VarNegEta(0), 
fTreeV0VarPtot(0), 
fTreeV0VarPt(0), 
fTreeV0VarPhi(0),
fTreeV0VarTheta(0),
//-----------INFO-FOR-CUTS------------------------
fTreeV0VarAlpha(0), 
fTreeV0VarPtArm(0),
fTreeV0VarDCAV0Dau(0), 
fTreeV0VarDCAV0ToPV(0), 
fTreeV0VarDCAPosToPV(0), 
fTreeV0VarDCANegToPV(0), 
fTreeV0VarCosPA(0), 
fTreeV0VarRadius(0), 

fTreeV0VarLeastNbrCrossedRows(0),
fTreeV0VarLeastNbrClusters(0),
fTreeV0VarLeastRatioCrossedRowsOverFindable(0),
fTreeV0VarMaxChi2PerCluster(0), 
fTreeV0VarMinTrackLength(0), 
//-----------DECAY-LENGTH-INFO--------------------
fTreeV0VarDistOverTotMom(0),
//---------------PID-TPC-INFO---------------------
fTreeV0VarPosNSigmaProton(0),
fTreeV0VarPosNSigmaPion(0),
fTreeV0VarNegNSigmaProton(0),
fTreeV0VarNegNSigmaPion(0),
//---------------PID-TOF-INFO---------------------
fTreeV0VarPosTOFNSigmaProton(0),
fTreeV0VarPosTOFNSigmaPion(0),
fTreeV0VarNegTOFNSigmaProton(0),
fTreeV0VarNegTOFNSigmaPion(0),
//---------------PID-ITS-INFO---------------------
fTreeV0VarPosITSNSigmaProton(0),
fTreeV0VarPosITSNSigmaPion(0),
fTreeV0VarNegITSNSigmaProton(0),
fTreeV0VarNegITSNSigmaPion(0),
//---Raw TPC dEdx + PIDForTracking information----
fTreeV0VarPosdEdx(0),
fTreeV0VarNegdEdx(0),
fTreeV0VarPosPIDForTracking(0), 
fTreeV0VarNegPIDForTracking(0),
//------------------------------------------------
fTreeV0VarChi2V0(0),         
//------------------------------------------------
fTreeV0VarPosTrackStatus(0), 
fTreeV0VarNegTrackStatus(0), 
//------------------------------------------------
fTreeV0VarPosDCAz(0), 
fTreeV0VarNegDCAz(0), 
//------------FULL-MOMENTUM-INFO------------------
fTreeV0VarPosPx(0), 
fTreeV0VarPosPy(0), 
fTreeV0VarPosPz(0), 
fTreeV0VarNegPx(0), 
fTreeV0VarNegPy(0), 
fTreeV0VarNegPz(0), 
//------------------------------------------------
fTreeV0VarPrimVertexX(0), 
fTreeV0VarPrimVertexY(0), 
fTreeV0VarPrimVertexZ(0), 
fTreeV0VarV0DecayX(0), 
fTreeV0VarV0DecayY(0), 
fTreeV0VarV0DecayZ(0), 
//---------------CLUSTER-INFO---------------------
fTreeV0VarPosITSClusters0(kFALSE),
fTreeV0VarPosITSClusters1(kFALSE),
fTreeV0VarPosITSClusters2(kFALSE),
fTreeV0VarPosITSClusters3(kFALSE),
fTreeV0VarPosITSClusters4(kFALSE),
fTreeV0VarPosITSClusters5(kFALSE),

fTreeV0VarNegITSClusters0(kFALSE),
fTreeV0VarNegITSClusters1(kFALSE),
fTreeV0VarNegITSClusters2(kFALSE),
fTreeV0VarNegITSClusters3(kFALSE),
fTreeV0VarNegITSClusters4(kFALSE),
fTreeV0VarNegITSClusters5(kFALSE),
//------------------------------------------------
fTreeV0VarPosITSSharedClusters0(kFALSE),
fTreeV0VarPosITSSharedClusters1(kFALSE),
fTreeV0VarPosITSSharedClusters2(kFALSE),
fTreeV0VarPosITSSharedClusters3(kFALSE),
fTreeV0VarPosITSSharedClusters4(kFALSE),
fTreeV0VarPosITSSharedClusters5(kFALSE),

fTreeV0VarNegITSSharedClusters0(kFALSE),
fTreeV0VarNegITSSharedClusters1(kFALSE),
fTreeV0VarNegITSSharedClusters2(kFALSE),
fTreeV0VarNegITSSharedClusters3(kFALSE),
fTreeV0VarNegITSSharedClusters4(kFALSE),
fTreeV0VarNegITSSharedClusters5(kFALSE),
//---------------OOB-PILEUP-INFO---------------------
fTreeV0VarNegTOFExpTDiff(0), 
fTreeV0VarPosTOFExpTDiff(0), 
fTreeV0VarNegTOFSignal(0), 
fTreeV0VarPosTOFSignal(0), 
fTreeV0VarNegTOFBCid(0), 
fTreeV0VarPosTOFBCid(0),  

//Event info
fTreeV0VarAmplitudeV0A(0), 
fTreeV0VarAmplitudeV0C(0), 
fTreeV0VarClosestNonEmptyBC(0), 

//Kink tagging
fTreeV0VarPosIsKink(kFALSE),
fTreeV0VarNegIsKink(kFALSE),

//Cowboy/sailor studies
fTreeV0VarIsCowboy(kFALSE), 

fTreeV0VarMagField(0),
fTreeV0VarRunNumber(0),
fTreeV0VarEventNumber(0), 
//===========================================================================================
//   Variables for Resonance Tree
//===========================================================================================
fTreeRsnVarCutIDrsn(0),
fTreeRsnVarPx(0),
fTreeRsnVarPy(0),
fTreeRsnVarPz(0),
fTreeRsnVarInvMass(0),
fTreeRsnVarPassesOOBPileupCut(kFALSE),
fTreeRsnVarEventNumber(0),
//===========================================================================================
//   Variables for Mixed Resonance Tree
//===========================================================================================
fTreeRsnFoundMixEvts(0),
fTreeRsnBkgVarCutIDrsn(0),
fTreeRsnBkgVarPx(0),
fTreeRsnBkgVarPy(0),
fTreeRsnBkgVarPz(0),
fTreeRsnBkgVarInvMass(0),
fTreeRsnBkgVarPassesOOBPileupCut(kFALSE),
fTreeRsnBkgVarEventNumber(0),
fDummyVarEventNumber(0), 
//===========================================================================================
//   Variables for Primary tracks Tree
//===========================================================================================
//-----------BASIC-INFO---------------------------
fTreePrimVarCharge(0),
fTreePrimVarRapPion(0), 
fTreePrimVarRapProton(0), 
fTreePrimVarRapKaon(0),
fTreePrimVarEta(0),
fTreePrimVarTheta(0),
fTreePrimVarPhi(0),
fTreePrimVarPtot(0), 
fTreePrimVarPt(0), 

fTreePrimVarDCAxyToPV(0),
fTreePrimVarDCAzToPV(0),

fTreePrimVarNbrCrossedRows(0),
fTreePrimVarNbrClusters(0),
fTreePrimVarRatioCrossedRowsOverFindable(0),
fTreePrimVarNbrCrossedRowsOverLength(0),
fTreePrimVarFractionSharedTPCClusters(0),
fTreePrimVarITSChi2PerCluster(0),
fTreePrimVarTPCChi2PerCluster(0),
fTreePrimVarTrackLength(0),
//---------------PID-TPC-INFO---------------------
fTreePrimVarNSigmaPion(0),
fTreePrimVarNSigmaKaon(0),
fTreePrimVarNSigmaProton(0),
//---------------PID-TOF-INFO---------------------
fTreePrimVarTOFNSigmaPion(0),
fTreePrimVarTOFNSigmaKaon(0),
fTreePrimVarTOFNSigmaProton(0),
//---------------PID-ITS-INFO---------------------
fTreePrimVarITSNSigmaPion(0),
fTreePrimVarITSNSigmaKaon(0),
fTreePrimVarITSNSigmaProton(0),
//---Raw TPC dEdx + PIDForTracking information----
fTreePrimVardEdx(0),
fTreePrimVarPIDForTracking(0),
//------------------------------------------------
fTreePrimVarTrackStatus(0),
//------------FULL-MOMENTUM-INFO------------------
fTreePrimVarPx(0),
fTreePrimVarPy(0),
fTreePrimVarPz(0),
//---------------CLUSTER-INFO---------------------
fTreePrimVarITSClusters0(0),
fTreePrimVarITSClusters1(0),
fTreePrimVarITSClusters2(0),
fTreePrimVarITSClusters3(0),
fTreePrimVarITSClusters4(0),
fTreePrimVarITSClusters5(0),

//------------------------------------------------
fTreePrimVarITSSharedClusters0(0),
fTreePrimVarITSSharedClusters1(0),
fTreePrimVarITSSharedClusters2(0),
fTreePrimVarITSSharedClusters3(0),
fTreePrimVarITSSharedClusters4(0),
fTreePrimVarITSSharedClusters5(0),

//---------------OOB-PILEUP-INFO---------------------
fTreePrimVarTOFExpTDiff(0),
fTreePrimVarTOFSignal(0),
fTreePrimVarTOFBCid(0),

//Kink tagging
fTreePrimVarIsKink(kFALSE),

fTreePrimVarRunNumber(0),
fTreePrimVarEventNumber(0),
//Histos
fHistEventCounter(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskStrangeCascadesTriggerAODRun2::AliAnalysisTaskStrangeCascadesTriggerAODRun2(const char* name, Bool_t lSaveV0s, Bool_t lSaveRsn, Bool_t lSavePrimaries, TString lcompType) : 
AliAnalysisTaskSE(name), fOutputList(0), 
fTreeCascade(0), fTreeV0(0), fTreeRsn(0), fTreeRsnBkg(0), fTreePrimTrack(0), fDummyTree(0),
fPIDResponse(0), fUtils(0), 
fRsnEvent(0), fRsnMiniEvent(0), fResonanceFinders(0), fRsnTrackCuts(0), fRsnPairCuts(0),
fPairName(""), fMotherMass(0),

//---> Flags controlling TTree outputs
fSaveV0s(kTRUE), 
fSaveRsn(kTRUE), 
fSavePrimaries(kTRUE), 
fComputationType(""),

//---> Variables controlling PV selections
fkMaxPVR2D(10.), 
fkMaxPVZ(10.), 

//---> Variables controlling cascade and V0 default selections
fMinNbrCrossedRows(0), 
fMinPtToSave(0), 
fMaxPtToSave(20), 
fMaxAbsEta(0.8), 
fMaxAbsRap(0.5),
fRejectCascKink(kFALSE),
fRejectV0Kink(kFALSE),

//---> Flags controlling cascade and V0 custom selections
fCascSaveAddConfig(kFALSE), 
fV0SaveAddConfig(kFALSE),
fPrimariesSaveAddConfig(kFALSE),
fCascadeResult(0), 
fV0Result(0),
fPrimTrackCuts(0), 

//---> Variables controlling event mixing
fNMix(0), 
fMaxDiffVz(0), 
fMaxDiffAngle(0), 
fMaxDiffMult(0),

//===========================================================================================
//   Variables for Cascade Tree
//===========================================================================================
//-----------BASIC-INFO---------------------------
fTreeCascVarCharge(0),
fTreeCascVarMassAsXi(0),
fTreeCascVarMassAsOmega(0),
fTreeCascVarPtot(0),
fTreeCascVarV0Ptot(0),
fTreeCascVarPt(0),
fTreeCascVarV0Pt(0),
fTreeCascVarRapXi(0),
fTreeCascVarRapOmega(0),
fTreeCascVarBachEta(0),
fTreeCascVarPosEta(0),
fTreeCascVarNegEta(0),
fTreeCascVarPhi(0),
fTreeCascVarTheta(0),
//-----------INFO-FOR-CUTS------------------------
fTreeCascVarAlpha(0),
fTreeCascVarPtArm(0),
fTreeCascVarAlphaV0(0),
fTreeCascVarPtArmV0(0),
fTreeCascVarDCACascDau(0),
fTreeCascVarDCABachToPV(0),
fTreeCascVarDCAV0Dau(0),
fTreeCascVarDCAV0ToPV(0),
fTreeCascVarDCAPosToPV(0),
fTreeCascVarDCANegToPV(0),
fTreeCascVarCascCosPA(0),
fTreeCascVarCascCosPASpecial(0),
    
fTreeCascVarCascRadius(0),
fTreeCascVarV0Mass(0),
fTreeCascVarV0MassAsLambda(0),
fTreeCascVarV0MassAsAntiLambda(0),
fTreeCascVarV0Radius(0),
fTreeCascVarV0CosPA(0),
fTreeCascVarWrongCosPA(0),
fTreeCascVarDCABachToBaryon(0),
fTreeCascVarLeastNbrCrossedRows(0),
fTreeCascVarLeastNbrClusters(0),
fTreeCascVarLeastRatioCrossedRowsOverFindable(0),
fTreeCascVarNbrCrossedRowsOverLength(0),
fTreeCascVarMaxChi2PerCluster(0),
fTreeCascVarMinTrackLength(0),
//-----------MULTIPLICITY-INFO--------------------
fTreeCascVarCentrality(0),          //! VZERO Percentile             - AliMultSelection
fTreeCascVarVZEROMultSel(0),        //! VZERO Multiplicity           - AliMultSelection
fTreeCascVarVZEROMultSig(0),        //! VZERO Multiplicity           - VZER0 signal
fTreeCascVarVZEROMultSigCorr(0),    //! Corrected VZERO Multiplicity - VZER0 signal
fTreeCascVarSPDMult(0),             //! SPD Multiplicity             - SPD Tracklets
fTreeCascVar_TriggerMask(0),
fTreeCascVarIsIncompleteDAQ(0),
fTreeCascVarSPDPileupFlag(0),
fTreeCascVarMVPileupFlag(0),
fTreeCascVarOOBPileupFlag(0),
//-----------DECAY-LENGTH-INFO--------------------
fTreeCascVarDistOverTotMom(0),
fTreeCascVarV0DistOverTotMom(0),
//------------------------------------------------
//---------------PID-TPC-INFO---------------------
fTreeCascVarBachNSigmaPion(0),
fTreeCascVarBachNSigmaKaon(0),
fTreeCascVarPosNSigmaProton(0),
fTreeCascVarPosNSigmaPion(0),
fTreeCascVarNegNSigmaProton(0),
fTreeCascVarNegNSigmaPion(0),
//---------------PID-TOF-INFO---------------------
fTreeCascVarBachTOFNSigmaPion(0),
fTreeCascVarBachTOFNSigmaKaon(0),
fTreeCascVarPosTOFNSigmaProton(0),
fTreeCascVarPosTOFNSigmaPion(0),
fTreeCascVarNegTOFNSigmaProton(0),
fTreeCascVarNegTOFNSigmaPion(0),
//---------------PID-ITS-INFO---------------------
fTreeCascVarBachITSNSigmaPion(0),
fTreeCascVarBachITSNSigmaKaon(0),
fTreeCascVarPosITSNSigmaProton(0),
fTreeCascVarPosITSNSigmaPion(0),
fTreeCascVarNegITSNSigmaProton(0),
fTreeCascVarNegITSNSigmaPion(0),
//---Raw TPC dEdx + PIDForTracking information----
fTreeCascVarPosdEdx(0),
fTreeCascVarNegdEdx(0),
fTreeCascVarBachdEdx(0),
fTreeCascVarPosPIDForTracking(0),
fTreeCascVarNegPIDForTracking(0),
fTreeCascVarBachPIDForTracking(0),
//------------------------------------------------
fTreeCascVarChi2Cascade(0),
fTreeCascVarChi2CascadePerNDF(0),
fTreeCascVarChi2V0(0),
//------------------------------------------------
fTreeCascVarBachTrackStatus(0),
fTreeCascVarPosTrackStatus(0),
fTreeCascVarNegTrackStatus(0),
//------------------------------------------------
fTreeCascVarBachDCAz(0),
fTreeCascVarPosDCAz(0),
fTreeCascVarNegDCAz(0),
//------------FULL-MOMENTUM-INFO------------------
fTreeCascVarBachPx(0),
fTreeCascVarBachPy(0),
fTreeCascVarBachPz(0),
fTreeCascVarPosPx(0),
fTreeCascVarPosPy(0),
fTreeCascVarPosPz(0),
fTreeCascVarNegPx(0),
fTreeCascVarNegPy(0),
fTreeCascVarNegPz(0),
//------------------------------------------------
fTreeCascVarPrimVertexX(0),
fTreeCascVarPrimVertexY(0),
fTreeCascVarPrimVertexZ(0),
fTreeCascVarCascadeDecayX(0),
fTreeCascVarCascadeDecayY(0),
fTreeCascVarCascadeDecayZ(0),
fTreeCascVarV0DecayX(0),
fTreeCascVarV0DecayY(0),
fTreeCascVarV0DecayZ(0),
//------------------------------------------------
fTreeCascVarBachIndex(0),
fTreeCascVarPosIndex(0),
fTreeCascVarNegIndex(0),
//---------------CLUSTER-INFO---------------------
fTreeCascVarBachITSClusters0(kFALSE),
fTreeCascVarBachITSClusters1(kFALSE),
fTreeCascVarBachITSClusters2(kFALSE),
fTreeCascVarBachITSClusters3(kFALSE),
fTreeCascVarBachITSClusters4(kFALSE),
fTreeCascVarBachITSClusters5(kFALSE),
    
fTreeCascVarPosITSClusters0(kFALSE),
fTreeCascVarPosITSClusters1(kFALSE),
fTreeCascVarPosITSClusters2(kFALSE),
fTreeCascVarPosITSClusters3(kFALSE),
fTreeCascVarPosITSClusters4(kFALSE),
fTreeCascVarPosITSClusters5(kFALSE),
    
fTreeCascVarNegITSClusters0(kFALSE),
fTreeCascVarNegITSClusters1(kFALSE),
fTreeCascVarNegITSClusters2(kFALSE),
fTreeCascVarNegITSClusters3(kFALSE),
fTreeCascVarNegITSClusters4(kFALSE),
fTreeCascVarNegITSClusters5(kFALSE),
    
//------------------------------------------------
fTreeCascVarPosITSSharedClusters0(kFALSE),
fTreeCascVarPosITSSharedClusters1(kFALSE),
fTreeCascVarPosITSSharedClusters2(kFALSE),
fTreeCascVarPosITSSharedClusters3(kFALSE),
fTreeCascVarPosITSSharedClusters4(kFALSE),
fTreeCascVarPosITSSharedClusters5(kFALSE),

fTreeCascVarNegITSSharedClusters0(kFALSE),
fTreeCascVarNegITSSharedClusters1(kFALSE),
fTreeCascVarNegITSSharedClusters2(kFALSE),
fTreeCascVarNegITSSharedClusters3(kFALSE),
fTreeCascVarNegITSSharedClusters4(kFALSE),
fTreeCascVarNegITSSharedClusters5(kFALSE),

fTreeCascVarBachITSSharedClusters0(kFALSE),
fTreeCascVarBachITSSharedClusters1(kFALSE),
fTreeCascVarBachITSSharedClusters2(kFALSE),
fTreeCascVarBachITSSharedClusters3(kFALSE),
fTreeCascVarBachITSSharedClusters4(kFALSE),
fTreeCascVarBachITSSharedClusters5(kFALSE),

//---------------OOB-PILEUP-INFO---------------------
fTreeCascVarBachTOFExpTDiff(0),
fTreeCascVarPosTOFExpTDiff(0),
fTreeCascVarNegTOFExpTDiff(0),

fTreeCascVarBachTOFSignal(0),
fTreeCascVarPosTOFSignal(0),
fTreeCascVarNegTOFSignal(0),
    
fTreeCascVarBachTOFBCid(0),
fTreeCascVarPosTOFBCid(0),
fTreeCascVarNegTOFBCid(0),

//Event info
fTreeCascVarAmplitudeV0A(0),
fTreeCascVarAmplitudeV0C(0),
fTreeCascVarClosestNonEmptyBC(0),
    
//Kink tagging
fTreeCascVarBachIsKink(kFALSE),
fTreeCascVarPosIsKink(kFALSE),
fTreeCascVarNegIsKink(kFALSE),
    
//Cowboy/sailor studies
fTreeCascVarIsCowboy(kFALSE),   //store if V0 is cowboy-like or sailor-like in XY plane
fTreeCascVarCowboyness(0), //negative -> cowboy, positive -> sailor
fTreeCascVarIsCascadeCowboy(kFALSE),   //store if V0 is cowboy-like or sailor-like in XY plane
fTreeCascVarCascadeCowboyness(0), //negative -> cowboy, positive -> sailor
    
fTreeCascVarMagField(0),
fTreeCascVarRunNumber(0),
fTreeCascVarEventNumber(0),
//===========================================================================================
//   Variables for V0 Tree
//===========================================================================================
//-----------BASIC-INFO---------------------------
fTreeV0VarMassAsK0s(0), 
fTreeV0VarMassAsLambda(0), 
fTreeV0VarMassAsAntiLambda(0), 
fTreeV0VarRapK0Short(0), 
fTreeV0VarRapLambda(0), 
fTreeV0VarPosEta(0), 
fTreeV0VarNegEta(0), 
fTreeV0VarPhi(0),
fTreeV0VarTheta(0),
fTreeV0VarPtot(0), 
fTreeV0VarPt(0), 
//-----------INFO-FOR-CUTS------------------------
fTreeV0VarAlpha(0), 
fTreeV0VarPtArm(0),
fTreeV0VarDCAV0Dau(0), 
fTreeV0VarDCAV0ToPV(0), 
fTreeV0VarDCAPosToPV(0), 
fTreeV0VarDCANegToPV(0), 
fTreeV0VarCosPA(0), 
fTreeV0VarRadius(0), 

fTreeV0VarLeastNbrCrossedRows(0),
fTreeV0VarLeastNbrClusters(0),
fTreeV0VarLeastRatioCrossedRowsOverFindable(0),
fTreeV0VarMaxChi2PerCluster(0), 
fTreeV0VarMinTrackLength(0), 
//-----------DECAY-LENGTH-INFO--------------------
fTreeV0VarDistOverTotMom(0),
//---------------PID-TPC-INFO---------------------
fTreeV0VarPosNSigmaProton(0),
fTreeV0VarPosNSigmaPion(0),
fTreeV0VarNegNSigmaProton(0),
fTreeV0VarNegNSigmaPion(0),
//---------------PID-TOF-INFO---------------------
fTreeV0VarPosTOFNSigmaProton(0),
fTreeV0VarPosTOFNSigmaPion(0),
fTreeV0VarNegTOFNSigmaProton(0),
fTreeV0VarNegTOFNSigmaPion(0),
//---------------PID-ITS-INFO---------------------
fTreeV0VarPosITSNSigmaProton(0),
fTreeV0VarPosITSNSigmaPion(0),
fTreeV0VarNegITSNSigmaProton(0),
fTreeV0VarNegITSNSigmaPion(0),
//---Raw TPC dEdx + PIDForTracking information----
fTreeV0VarPosdEdx(0),
fTreeV0VarNegdEdx(0),
fTreeV0VarPosPIDForTracking(0), 
fTreeV0VarNegPIDForTracking(0),
//------------------------------------------------
fTreeV0VarChi2V0(0),         
//------------------------------------------------
fTreeV0VarPosTrackStatus(0), 
fTreeV0VarNegTrackStatus(0), 
//------------------------------------------------
fTreeV0VarPosDCAz(0), 
fTreeV0VarNegDCAz(0), 
//------------FULL-MOMENTUM-INFO------------------
fTreeV0VarPosPx(0), 
fTreeV0VarPosPy(0), 
fTreeV0VarPosPz(0), 
fTreeV0VarNegPx(0), 
fTreeV0VarNegPy(0), 
fTreeV0VarNegPz(0), 
//------------------------------------------------
fTreeV0VarPrimVertexX(0), 
fTreeV0VarPrimVertexY(0), 
fTreeV0VarPrimVertexZ(0), 
fTreeV0VarV0DecayX(0), 
fTreeV0VarV0DecayY(0), 
fTreeV0VarV0DecayZ(0), 
//---------------CLUSTER-INFO---------------------
fTreeV0VarPosITSClusters0(kFALSE),
fTreeV0VarPosITSClusters1(kFALSE),
fTreeV0VarPosITSClusters2(kFALSE),
fTreeV0VarPosITSClusters3(kFALSE),
fTreeV0VarPosITSClusters4(kFALSE),
fTreeV0VarPosITSClusters5(kFALSE),

fTreeV0VarNegITSClusters0(kFALSE),
fTreeV0VarNegITSClusters1(kFALSE),
fTreeV0VarNegITSClusters2(kFALSE),
fTreeV0VarNegITSClusters3(kFALSE),
fTreeV0VarNegITSClusters4(kFALSE),
fTreeV0VarNegITSClusters5(kFALSE),
//------------------------------------------------
fTreeV0VarPosITSSharedClusters0(kFALSE),
fTreeV0VarPosITSSharedClusters1(kFALSE),
fTreeV0VarPosITSSharedClusters2(kFALSE),
fTreeV0VarPosITSSharedClusters3(kFALSE),
fTreeV0VarPosITSSharedClusters4(kFALSE),
fTreeV0VarPosITSSharedClusters5(kFALSE),

fTreeV0VarNegITSSharedClusters0(kFALSE),
fTreeV0VarNegITSSharedClusters1(kFALSE),
fTreeV0VarNegITSSharedClusters2(kFALSE),
fTreeV0VarNegITSSharedClusters3(kFALSE),
fTreeV0VarNegITSSharedClusters4(kFALSE),
fTreeV0VarNegITSSharedClusters5(kFALSE),
//---------------OOB-PILEUP-INFO---------------------
fTreeV0VarNegTOFExpTDiff(0), 
fTreeV0VarPosTOFExpTDiff(0), 
fTreeV0VarNegTOFSignal(0), 
fTreeV0VarPosTOFSignal(0), 
fTreeV0VarNegTOFBCid(0), 
fTreeV0VarPosTOFBCid(0),  

//Event info
fTreeV0VarAmplitudeV0A(0), 
fTreeV0VarAmplitudeV0C(0), 
fTreeV0VarClosestNonEmptyBC(0), 

//Kink tagging
fTreeV0VarPosIsKink(kFALSE),
fTreeV0VarNegIsKink(kFALSE),

//Cowboy/sailor studies
fTreeV0VarIsCowboy(kFALSE), 

fTreeV0VarMagField(0),
fTreeV0VarRunNumber(0),
fTreeV0VarEventNumber(0), 
//===========================================================================================
//   Variables for Resonance Tree
//===========================================================================================
fTreeRsnVarCutIDrsn(0),
fTreeRsnVarPx(0),
fTreeRsnVarPy(0),
fTreeRsnVarPz(0),
fTreeRsnVarInvMass(0),
fTreeRsnVarPassesOOBPileupCut(kFALSE),
fTreeRsnVarEventNumber(0),
//===========================================================================================
//   Variables for Mixed Resonance Tree
//===========================================================================================
fTreeRsnFoundMixEvts(0),
fTreeRsnBkgVarCutIDrsn(0),
fTreeRsnBkgVarPx(0),
fTreeRsnBkgVarPy(0),
fTreeRsnBkgVarPz(0),
fTreeRsnBkgVarInvMass(0),
fTreeRsnBkgVarPassesOOBPileupCut(kFALSE),
fTreeRsnBkgVarEventNumber(0),
fDummyVarEventNumber(0),
//===========================================================================================
//   Variables for Primary tracks Tree
//===========================================================================================
//-----------BASIC-INFO---------------------------
fTreePrimVarCharge(0),
fTreePrimVarRapPion(0), 
fTreePrimVarRapProton(0), 
fTreePrimVarRapKaon(0),
fTreePrimVarEta(0),
fTreePrimVarTheta(0),
fTreePrimVarPhi(0),
fTreePrimVarPtot(0), 
fTreePrimVarPt(0), 

fTreePrimVarDCAxyToPV(0),
fTreePrimVarDCAzToPV(0),

fTreePrimVarNbrCrossedRows(0),
fTreePrimVarNbrClusters(0),
fTreePrimVarRatioCrossedRowsOverFindable(0),
fTreePrimVarNbrCrossedRowsOverLength(0),
fTreePrimVarFractionSharedTPCClusters(0),
fTreePrimVarITSChi2PerCluster(0),
fTreePrimVarTPCChi2PerCluster(0),
fTreePrimVarTrackLength(0),
//---------------PID-TPC-INFO---------------------
fTreePrimVarNSigmaPion(0),
fTreePrimVarNSigmaKaon(0),
fTreePrimVarNSigmaProton(0),
//---------------PID-TOF-INFO---------------------
fTreePrimVarTOFNSigmaPion(0),
fTreePrimVarTOFNSigmaKaon(0),
fTreePrimVarTOFNSigmaProton(0),
//---------------PID-ITS-INFO---------------------
fTreePrimVarITSNSigmaPion(0),
fTreePrimVarITSNSigmaKaon(0),
fTreePrimVarITSNSigmaProton(0),
//---Raw TPC dEdx + PIDForTracking information----
fTreePrimVardEdx(0),
fTreePrimVarPIDForTracking(0),
//------------------------------------------------
fTreePrimVarTrackStatus(0),
//------------FULL-MOMENTUM-INFO------------------
fTreePrimVarPx(0),
fTreePrimVarPy(0),
fTreePrimVarPz(0),
//---------------CLUSTER-INFO---------------------
fTreePrimVarITSClusters0(0),
fTreePrimVarITSClusters1(0),
fTreePrimVarITSClusters2(0),
fTreePrimVarITSClusters3(0),
fTreePrimVarITSClusters4(0),
fTreePrimVarITSClusters5(0),

//------------------------------------------------
fTreePrimVarITSSharedClusters0(0),
fTreePrimVarITSSharedClusters1(0),
fTreePrimVarITSSharedClusters2(0),
fTreePrimVarITSSharedClusters3(0),
fTreePrimVarITSSharedClusters4(0),
fTreePrimVarITSSharedClusters5(0),

//---------------OOB-PILEUP-INFO---------------------
fTreePrimVarTOFExpTDiff(0),
fTreePrimVarTOFSignal(0),
fTreePrimVarTOFBCid(0),

//Kink tagging
fTreePrimVarIsKink(kFALSE),

fTreePrimVarRunNumber(0),
fTreePrimVarEventNumber(0),
//-----------BASIC-INFO---------------------------
//Histos
fHistEventCounter(0)
{
    fSaveV0s            = lSaveV0s;
    fSaveRsn            = lSaveRsn;
    fSavePrimaries      = lSavePrimaries;
    fComputationType    = lcompType;
    
    // constructor
    DefineInput(0, TChain::Class());    
    
    DefineOutput(1, TList::Class());    
    DefineOutput(2, TTree::Class());//cascades
    DefineOutput(3, TTree::Class());//V0s
    DefineOutput(4, TTree::Class());//Resonances
    DefineOutput(5, TTree::Class());//Mixed events resonances
    DefineOutput(6, TTree::Class());//Primaries
}
//_____________________________________________________________________________
AliAnalysisTaskStrangeCascadesTriggerAODRun2::~AliAnalysisTaskStrangeCascadesTriggerAODRun2()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     
    }
    if(fTreeV0) {
        delete fTreeV0;    
    }
    if(fTreeCascade) {
        delete fTreeCascade;    
    }
    if(fTreeRsn) {
        delete fTreeRsn;    
    }
    if(fTreeRsnBkg) {
        delete fTreeRsnBkg;    
    }
    if(fTreePrimTrack) {
        delete fTreePrimTrack;    
    }
    if(fDummyTree) {
        delete fDummyTree;    
    }
    if(fUtils) {
        delete fUtils;    
    }
    if(fRsnEvent) {
        delete fRsnEvent;    
    }
    if(fRsnMiniEvent) {
        delete fRsnMiniEvent;    
    }
    if(fPIDResponse) {
        delete fPIDResponse;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskStrangeCascadesTriggerAODRun2::UserCreateOutputObjects()
{
    // check computation type
    Bool_t okComp = kFALSE;
    if (fComputationType == "")         okComp = kTRUE;
    if (fComputationType == "MIX")      okComp = kTRUE;
    if (fComputationType == "ROTATE1")  okComp = kTRUE;
    if (fComputationType == "ROTATE2")  okComp = kTRUE;
    if (!okComp) {
        AliError(Form("[%s] Unknown computation type", GetName()));
        return;
    }
    
    //===========================================================================================
    //   Create Cascade output tree
    //===========================================================================================
    fTreeCascade = new TTree("fTreeCascade","CascadeCandidates");
    //-----------BASIC-INFO---------------------------
    fTreeCascade->Branch("fTreeCascVarCharge"           , &fTreeCascVarCharge           , "fTreeCascVarCharge/I");
    fTreeCascade->Branch("fTreeCascVarMassAsXi"         , &fTreeCascVarMassAsXi         , "fTreeCascVarMassAsXi/D");
    fTreeCascade->Branch("fTreeCascVarMassAsOmega"      , &fTreeCascVarMassAsOmega      , "fTreeCascVarMassAsOmega/D");
    fTreeCascade->Branch("fTreeCascVarPtot"             , &fTreeCascVarPtot             , "fTreeCascVarPtot/D");
    fTreeCascade->Branch("fTreeCascVarPt"               , &fTreeCascVarPt               , "fTreeCascVarPt/D");
    fTreeCascade->Branch("fTreeCascVarRapXi"            , &fTreeCascVarRapXi            , "fTreeCascVarRapXi/D");
    fTreeCascade->Branch("fTreeCascVarRapOmega"         , &fTreeCascVarRapOmega         , "fTreeCascVarRapOmega/D");
    fTreeCascade->Branch("fTreeCascVarBachEta"          , &fTreeCascVarBachEta          , "fTreeCascVarBachEta/D");
    fTreeCascade->Branch("fTreeCascVarPosEta"           , &fTreeCascVarPosEta           , "fTreeCascVarPosEta/D");
    fTreeCascade->Branch("fTreeCascVarNegEta"           , &fTreeCascVarNegEta           , "fTreeCascVarNegEta/D");
    fTreeCascade->Branch("fTreeCascVarPhi"              , &fTreeCascVarPhi              , "fTreeCascVarPhi/D");
    fTreeCascade->Branch("fTreeCascVarTheta"            , &fTreeCascVarTheta            , "fTreeCascVarTheta/D");
    //-----------INFO-FOR-CUTS------------------------
    fTreeCascade->Branch("fTreeCascVarAlpha"            , &fTreeCascVarAlpha            , "fTreeCascVarAlpha/D");
    fTreeCascade->Branch("fTreeCascVarPtArm"            , &fTreeCascVarPtArm            , "fTreeCascVarPtArm/D");
    fTreeCascade->Branch("fTreeCascVarAlphaV0"          , &fTreeCascVarAlphaV0          , "fTreeCascVarAlphaV0/D");
    fTreeCascade->Branch("fTreeCascVarPtArmV0"          , &fTreeCascVarPtArmV0          , "fTreeCascVarPtArmV0/D");
    fTreeCascade->Branch("fTreeCascVarDCACascDau"       , &fTreeCascVarDCACascDau       , "fTreeCascVarDCACascDau/D");
    fTreeCascade->Branch("fTreeCascVarDCABachToPV"      , &fTreeCascVarDCABachToPV      , "fTreeCascVarDCABachToPV/D");
    fTreeCascade->Branch("fTreeCascVarDCAV0Dau"         , &fTreeCascVarDCAV0Dau         , "fTreeCascVarDCAV0Dau/D");
    fTreeCascade->Branch("fTreeCascVarDCAV0ToPV"        , &fTreeCascVarDCAV0ToPV        , "fTreeCascVarDCAV0ToPV/D");
    fTreeCascade->Branch("fTreeCascVarDCAPosToPV"       , &fTreeCascVarDCAPosToPV       , "fTreeCascVarDCAPosToPV/D");
    fTreeCascade->Branch("fTreeCascVarDCANegToPV"       , &fTreeCascVarDCANegToPV       , "fTreeCascVarDCANegToPV/D");
    fTreeCascade->Branch("fTreeCascVarCascCosPA"        , &fTreeCascVarCascCosPA        , "fTreeCascVarCascCosPA/D");
    fTreeCascade->Branch("fTreeCascVarCascCosPASpecial" , &fTreeCascVarCascCosPASpecial , "fTreeCascVarCascCosPASpecial/D");
    
    fTreeCascade->Branch("fTreeCascVarCascRadius"       , &fTreeCascVarCascRadius       , "fTreeCascVarCascRadius/D");
    fTreeCascade->Branch("fTreeCascVarV0Mass"           , &fTreeCascVarV0Mass           , "fTreeCascVarV0Mass/D");
    fTreeCascade->Branch("fTreeCascVarV0MassAsLambda"   , &fTreeCascVarV0MassAsLambda   , "fTreeCascVarV0MassAsLambda/D");
    fTreeCascade->Branch("fTreeCascVarV0MassAsLambdaAnti", &fTreeCascVarV0MassAsAntiLambda, "fTreeCascVarV0MassAsAntiLambda/D");
    fTreeCascade->Branch("fTreeCascVarV0Radius"          , &fTreeCascVarV0Radius          , "fTreeCascVarV0Radius/D");
    fTreeCascade->Branch("fTreeCascVarV0CosPA"          , &fTreeCascVarV0CosPA          , "fTreeCascVarV0CosPA/D");
    fTreeCascade->Branch("fTreeCascVarWrongCosPA"       , &fTreeCascVarWrongCosPA       , "fTreeCascVarWrongCosPA/D");
    fTreeCascade->Branch("fTreeCascVarDCABachToBaryon"  , &fTreeCascVarDCABachToBaryon  , "fTreeCascVarDCABachToBaryon/D");
    fTreeCascade->Branch("fTreeCascVarLeastNbrCrossedRows" , &fTreeCascVarLeastNbrCrossedRows  , "fTreeCascVarLeastNbrCrossedRows/I");
    fTreeCascade->Branch("fTreeCascVarLeastRatioCrossedRowsOverFindable"   , &fTreeCascVarLeastRatioCrossedRowsOverFindable   , "fTreeCascVarLeastRatioCrossedRowsOverFindable/D");
    fTreeCascade->Branch("fTreeCascVarLeastNbrClusters" , &fTreeCascVarLeastNbrClusters  , "fTreeCascVarLeastNbrClusters/I");
    fTreeCascade->Branch("fTreeCascVarNbrCrossedRowsOverLength" , &fTreeCascVarNbrCrossedRowsOverLength  , "fTreeCascVarNbrCrossedRowsOverLength/D");
    fTreeCascade->Branch("fTreeCascVarMaxChi2PerCluster", &fTreeCascVarMaxChi2PerCluster, "fTreeCascVarMaxChi2PerCluster/D");
    fTreeCascade->Branch("fTreeCascVarMinTrackLength"   , &fTreeCascVarMinTrackLength   , "fTreeCascVarMinTrackLength/D");
    //-----------MULTIPLICITY-INFO--------------------
    fTreeCascade->Branch("fTreeCascVarCentrality"       , &fTreeCascVarCentrality       , "fTreeCascVarCentrality/D");
    fTreeCascade->Branch("fTreeCascVarVZEROMultSel"     , &fTreeCascVarVZEROMultSel     , "fTreeCascVarVZEROMultSel/D");
    fTreeCascade->Branch("fTreeCascVarVZEROMultSig"     , &fTreeCascVarVZEROMultSig     , "fTreeCascVarVZEROMultSig/I");
    fTreeCascade->Branch("fTreeCascVarVZEROMultSigCorr" , &fTreeCascVarVZEROMultSigCorr , "fTreeCascVarVZEROMultSigCorr/I");
    fTreeCascade->Branch("fTreeCascVarSPDMult"          , &fTreeCascVarSPDMult          , "fTreeCascVarSPDMult/I");
    fTreeCascade->Branch("fTreeCascVar_TriggerMask"     , &fTreeCascVar_TriggerMask     , "fTreeCascVar_TriggerMask/i");
    fTreeCascade->Branch("fTreeCascVarIsIncompleteDAQ"  , &fTreeCascVarIsIncompleteDAQ  , "fTreeCascVarIsIncompleteDAQ/O");
    fTreeCascade->Branch("fTreeCascVarSPDPileupFlag"    , &fTreeCascVarSPDPileupFlag    , "fTreeCascVarSPDPileupFlag/O");
    fTreeCascade->Branch("fTreeCascVarMVPileupFlag"     , &fTreeCascVarMVPileupFlag     , "fTreeCascVarMVPileupFlag/O");
    fTreeCascade->Branch("fTreeCascVarOOBPileupFlag"    , &fTreeCascVarOOBPileupFlag    , "fTreeCascVarOOBPileupFlag/O");
    //-----------DECAY-LENGTH-INFO--------------------
    fTreeCascade->Branch("fTreeCascVarDistOverTotMom"   , &fTreeCascVarDistOverTotMom   , "fTreeCascVarDistOverTotMom/D");
    fTreeCascade->Branch("fTreeCascVarV0DistOverTotMom" , &fTreeCascVarV0DistOverTotMom , "fTreeCascVarV0DistOverTotMom/D");
    //---------------PID-TPC-INFO---------------------
    fTreeCascade->Branch("fTreeCascVarBachNSigmaPion"   , &fTreeCascVarBachNSigmaPion   , "fTreeCascVarBachNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarBachNSigmaKaon"   , &fTreeCascVarBachNSigmaKaon   , "fTreeCascVarBachNSigmaKaon/F");
    fTreeCascade->Branch("fTreeCascVarPosNSigmaProton"  , &fTreeCascVarPosNSigmaProton  , "fTreeCascVarPosNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarPosNSigmaPion"    , &fTreeCascVarPosNSigmaPion    , "fTreeCascVarPosNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarNegNSigmaProton"  , &fTreeCascVarNegNSigmaProton  , "fTreeCascVarNegNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarNegNSigmaPion"    , &fTreeCascVarNegNSigmaPion    , "fTreeCascVarNegNSigmaPion/F");
    //---------------PID-TOF-INFO---------------------
    fTreeCascade->Branch("fTreeCascVarBachTOFNSigmaPion", &fTreeCascVarBachTOFNSigmaPion , "fTreeCascVarBachTOFNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarBachTOFNSigmaKaon", &fTreeCascVarBachTOFNSigmaKaon , "fTreeCascVarBachTOFNSigmaKaon/F");
    fTreeCascade->Branch("fTreeCascVarPosTOFNSigmaProton" , &fTreeCascVarPosTOFNSigmaProton  , "fTreeCascVarPosTOFNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarPosTOFNSigmaPion" , &fTreeCascVarPosTOFNSigmaPion  , "fTreeCascVarPosTOFNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarNegTOFNSigmaProton" , &fTreeCascVarNegTOFNSigmaProton  , "fTreeCascVarNegTOFNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarNegTOFNSigmaPion" , &fTreeCascVarNegTOFNSigmaPion  , "fTreeCascVarNegTOFNSigmaPion/F");
    //---------------PID-ITS-INFO---------------------
    fTreeCascade->Branch("fTreeCascVarBachITSNSigmaPion", &fTreeCascVarBachITSNSigmaPion , "fTreeCascVarBachITSNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarBachITSNSigmaKaon", &fTreeCascVarBachITSNSigmaKaon , "fTreeCascVarBachITSNSigmaKaon/F");
    fTreeCascade->Branch("fTreeCascVarPosITSNSigmaProton" , &fTreeCascVarPosITSNSigmaProton  , "fTreeCascVarPosITSNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarPosITSNSigmaPion" , &fTreeCascVarPosITSNSigmaPion  , "fTreeCascVarPosITSNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarNegITSNSigmaProton" , &fTreeCascVarNegITSNSigmaProton  , "fTreeCascVarNegITSNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarNegITSNSigmaPion" , &fTreeCascVarNegITSNSigmaPion  , "fTreeCascVarNegITSNSigmaPion/F");
    //------------------------------------------------
    fTreeCascade->Branch("fTreeCascVarChi2Cascade"      , &fTreeCascVarChi2Cascade      , "fTreeCascVarChi2Cascade/D");
    fTreeCascade->Branch("fTreeCascVarChi2CascadePerNDF", &fTreeCascVarChi2CascadePerNDF, "fTreeCascVarChi2CascadePerNDF/D");
    fTreeCascade->Branch("fTreeCascVarChi2V0"           , &fTreeCascVarChi2V0           , "fTreeCascVarChi2V0/D");
    //------------------------------------------------
    fTreeCascade->Branch("fTreeCascVarBachdEdx"         , &fTreeCascVarBachdEdx         , "fTreeCascVarBachdEdx/D");
    fTreeCascade->Branch("fTreeCascVarPosdEdx"          , &fTreeCascVarPosdEdx          , "fTreeCascVarPosdEdx/D");
    fTreeCascade->Branch("fTreeCascVarNegdEdx"          , &fTreeCascVarNegdEdx          , "fTreeCascVarNegdEdx/D");
    fTreeCascade->Branch("fTreeCascVarBachPIDForTracking"   , &fTreeCascVarBachPIDForTracking   , "fTreeCascVarBachPIDForTracking/D");
    fTreeCascade->Branch("fTreeCascVarPosPIDForTracking"    , &fTreeCascVarPosPIDForTracking    , "fTreeCascVarPosPIDForTracking/D");
    fTreeCascade->Branch("fTreeCascVarNegPIDForTracking"    , &fTreeCascVarNegPIDForTracking    , "fTreeCascVarNegPIDForTracking/D");
    //------------------------------------------------
    fTreeCascade->Branch("fTreeCascVarBachTrackStatus"  , &fTreeCascVarBachTrackStatus  , "fTreeCascVarBachTrackStatus/l");
    fTreeCascade->Branch("fTreeCascVarPosTrackStatus"   , &fTreeCascVarPosTrackStatus   , "fTreeCascVarPosTrackStatus/l");
    fTreeCascade->Branch("fTreeCascVarNegTrackStatus"   , &fTreeCascVarNegTrackStatus   , "fTreeCascVarNegTrackStatus/l");
    //------------------------------------------------
    fTreeCascade->Branch("fTreeCascVarBachDCAz"         , &fTreeCascVarBachDCAz         , "fTreeCascVarBachDCAz/D");
    fTreeCascade->Branch("fTreeCascVarPosDCAz"          , &fTreeCascVarPosDCAz          , "fTreeCascVarPosDCAz/D");
    fTreeCascade->Branch("fTreeCascVarNegDCAz"          , &fTreeCascVarNegDCAz          , "fTreeCascVarNegDCAz/D");
    //------------FULL-MOMENTUM-INFO------------------
    fTreeCascade->Branch("fTreeCascVarBachPx"           , &fTreeCascVarBachPx           , "fTreeCascVarBachPx/D");
    fTreeCascade->Branch("fTreeCascVarBachPy"           , &fTreeCascVarBachPy           , "fTreeCascVarBachPy/D");
    fTreeCascade->Branch("fTreeCascVarBachPz"           , &fTreeCascVarBachPz           , "fTreeCascVarBachPz/D");
    fTreeCascade->Branch("fTreeCascVarPosPx"            , &fTreeCascVarPosPx            , "fTreeCascVarPosPx/D");
    fTreeCascade->Branch("fTreeCascVarPosPy"            , &fTreeCascVarPosPy            , "fTreeCascVarPosPy/D");
    fTreeCascade->Branch("fTreeCascVarPosPz"            , &fTreeCascVarPosPz            , "fTreeCascVarPosPz/D");
    fTreeCascade->Branch("fTreeCascVarNegPx"            , &fTreeCascVarNegPx            , "fTreeCascVarNegPx/D");
    fTreeCascade->Branch("fTreeCascVarNegPy"            , &fTreeCascVarNegPy            , "fTreeCascVarNegPy/D");
    fTreeCascade->Branch("fTreeCascVarNegPz"            , &fTreeCascVarNegPz            , "fTreeCascVarNegPz/D");
    //------------------------------------------------
    fTreeCascade->Branch("fTreeCascVarPrimVertexX"      , &fTreeCascVarPrimVertexX      , "fTreeCascVarPrimVertexX/D");
    fTreeCascade->Branch("fTreeCascVarPrimVertexY"      , &fTreeCascVarPrimVertexY      , "fTreeCascVarPrimVertexY/D");
    fTreeCascade->Branch("fTreeCascVarPrimVertexZ"      , &fTreeCascVarPrimVertexZ      , "fTreeCascVarPrimVertexZ/D");
    fTreeCascade->Branch("fTreeCascVarCascadeDecayX"    , &fTreeCascVarCascadeDecayX    , "fTreeCascVarCascadeDecayX/D");
    fTreeCascade->Branch("fTreeCascVarCascadeDecayY"    , &fTreeCascVarCascadeDecayY    , "fTreeCascVarCascadeDecayY/D");
    fTreeCascade->Branch("fTreeCascVarCascadeDecayZ"    , &fTreeCascVarCascadeDecayZ    , "fTreeCascVarCascadeDecayZ/D");
    fTreeCascade->Branch("fTreeCascVarV0DecayX"         , &fTreeCascVarV0DecayX         , "fTreeCascVarV0DecayX/D");
    fTreeCascade->Branch("fTreeCascVarV0DecayY"         , &fTreeCascVarV0DecayY         , "fTreeCascVarV0DecayY/D");
    fTreeCascade->Branch("fTreeCascVarV0DecayZ"         , &fTreeCascVarV0DecayZ         , "fTreeCascVarV0DecayZ/D");
    //------------------------------------------------
    //fTreeCascade->Branch("fTreeCascVarBachIndex"        , &fTreeCascVarBachIndex        , "fTreeCascVarBachIndex/I");
    //fTreeCascade->Branch("fTreeCascVarPosIndex"         , &fTreeCascVarPosIndex         , "fTreeCascVarPosIndex/I");
    //fTreeCascade->Branch("fTreeCascVarNegIndex"         , &fTreeCascVarNegIndex         , "fTreeCascVarNegIndex/I");
    //Event Number (check same-event index mixups)
    fTreeCascade->Branch("fTreeCascVarEventNumber"      , &fTreeCascVarEventNumber      , "fTreeCascVarEventNumber/l");
    //---------------CLUSTER-INFO---------------------
    fTreeCascade->Branch("fTreeCascVarBachITSClusters0" , &fTreeCascVarBachITSClusters0 , "fTreeCascVarBachITSClusters0/O");
    fTreeCascade->Branch("fTreeCascVarBachITSClusters1" , &fTreeCascVarBachITSClusters1 , "fTreeCascVarBachITSClusters1/O");
    fTreeCascade->Branch("fTreeCascVarBachITSClusters2" , &fTreeCascVarBachITSClusters2 , "fTreeCascVarBachITSClusters2/O");
    fTreeCascade->Branch("fTreeCascVarBachITSClusters3" , &fTreeCascVarBachITSClusters3 , "fTreeCascVarBachITSClusters3/O");
    fTreeCascade->Branch("fTreeCascVarBachITSClusters4" , &fTreeCascVarBachITSClusters4 , "fTreeCascVarBachITSClusters4/O");
    fTreeCascade->Branch("fTreeCascVarBachITSClusters5" , &fTreeCascVarBachITSClusters5 , "fTreeCascVarBachITSClusters5/O");
    
    fTreeCascade->Branch("fTreeCascVarPosITSClusters0"  , &fTreeCascVarPosITSClusters0  , "fTreeCascVarPosITSClusters0/O");
    fTreeCascade->Branch("fTreeCascVarPosITSClusters1"  , &fTreeCascVarPosITSClusters1  , "fTreeCascVarPosITSClusters1/O");
    fTreeCascade->Branch("fTreeCascVarPosITSClusters2"  , &fTreeCascVarPosITSClusters2  , "fTreeCascVarPosITSClusters2/O");
    fTreeCascade->Branch("fTreeCascVarPosITSClusters3"  , &fTreeCascVarPosITSClusters3  , "fTreeCascVarPosITSClusters3/O");
    fTreeCascade->Branch("fTreeCascVarPosITSClusters4"  , &fTreeCascVarPosITSClusters4  , "fTreeCascVarPosITSClusters4/O");
    fTreeCascade->Branch("fTreeCascVarPosITSClusters5"  , &fTreeCascVarPosITSClusters5  , "fTreeCascVarPosITSClusters5/O");
    
    fTreeCascade->Branch("fTreeCascVarNegITSClusters0"  , &fTreeCascVarNegITSClusters0  , "fTreeCascVarNegITSClusters0/O");
    fTreeCascade->Branch("fTreeCascVarNegITSClusters1"  , &fTreeCascVarNegITSClusters1  , "fTreeCascVarNegITSClusters1/O");
    fTreeCascade->Branch("fTreeCascVarNegITSClusters2"  , &fTreeCascVarNegITSClusters2  , "fTreeCascVarNegITSClusters2/O");
    fTreeCascade->Branch("fTreeCascVarNegITSClusters3"  , &fTreeCascVarNegITSClusters3  , "fTreeCascVarNegITSClusters3/O");
    fTreeCascade->Branch("fTreeCascVarNegITSClusters4"  , &fTreeCascVarNegITSClusters4  , "fTreeCascVarNegITSClusters4/O");
    fTreeCascade->Branch("fTreeCascVarNegITSClusters5"  , &fTreeCascVarNegITSClusters5  , "fTreeCascVarNegITSClusters5/O");
    //------------------------------------------------
    fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters0" , &fTreeCascVarBachITSSharedClusters0  , "fTreeCascVarBachITSSharedClusters0/O");
    fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters1" , &fTreeCascVarBachITSSharedClusters1  , "fTreeCascVarBachITSSharedClusters1/O");
    fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters2" , &fTreeCascVarBachITSSharedClusters2  , "fTreeCascVarBachITSSharedClusters2/O");
    fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters3" , &fTreeCascVarBachITSSharedClusters3  , "fTreeCascVarBachITSSharedClusters3/O");
    fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters4" , &fTreeCascVarBachITSSharedClusters4  , "fTreeCascVarBachITSSharedClusters4/O");
    fTreeCascade->Branch("fTreeCascVarBachITSSharedClusters5" , &fTreeCascVarBachITSSharedClusters5  , "fTreeCascVarBachITSSharedClusters5/O");
    
    fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters0" , &fTreeCascVarPosITSSharedClusters0  , "fTreeCascVarPosITSSharedClusters0/O");
    fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters1" , &fTreeCascVarPosITSSharedClusters1  , "fTreeCascVarPosITSSharedClusters1/O");
    fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters2" , &fTreeCascVarPosITSSharedClusters2  , "fTreeCascVarPosITSSharedClusters2/O");
    fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters3" , &fTreeCascVarPosITSSharedClusters3  , "fTreeCascVarPosITSSharedClusters3/O");
    fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters4" , &fTreeCascVarPosITSSharedClusters4  , "fTreeCascVarPosITSSharedClusters4/O");
    fTreeCascade->Branch("fTreeCascVarPosITSSharedClusters5" , &fTreeCascVarPosITSSharedClusters5  , "fTreeCascVarPosITSSharedClusters5/O");
    
    fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters0" , &fTreeCascVarNegITSSharedClusters0  , "fTreeCascVarNegITSSharedClusters0/O");
    fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters1" , &fTreeCascVarNegITSSharedClusters1  , "fTreeCascVarNegITSSharedClusters1/O");
    fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters2" , &fTreeCascVarNegITSSharedClusters2  , "fTreeCascVarNegITSSharedClusters2/O");
    fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters3" , &fTreeCascVarNegITSSharedClusters3  , "fTreeCascVarNegITSSharedClusters3/O");
    fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters4" , &fTreeCascVarNegITSSharedClusters4  , "fTreeCascVarNegITSSharedClusters4/O");
    fTreeCascade->Branch("fTreeCascVarNegITSSharedClusters5" , &fTreeCascVarNegITSSharedClusters5  , "fTreeCascVarNegITSSharedClusters5/O");

    //---------------OOB-PILEUP-INFO---------------------
    fTreeCascade->Branch("fTreeCascVarBachTOFExpTDiff"  , &fTreeCascVarBachTOFExpTDiff  , "fTreeCascVarBachTOFExpTDiff/D");
    fTreeCascade->Branch("fTreeCascVarPosTOFExpTDiff"   , &fTreeCascVarPosTOFExpTDiff   , "fTreeCascVarPosTOFExpTDiff/D");
    fTreeCascade->Branch("fTreeCascVarNegTOFExpTDiff"   , &fTreeCascVarNegTOFExpTDiff   , "fTreeCascVarNegTOFExpTDiff/D");

    fTreeCascade->Branch("fTreeCascVarBachTOFSignal"    , &fTreeCascVarBachTOFSignal    , "fTreeCascVarBachTOFSignal/D");
    fTreeCascade->Branch("fTreeCascVarPosTOFSignal"     , &fTreeCascVarPosTOFSignal     , "fTreeCascVarPosTOFSignal/D");
    fTreeCascade->Branch("fTreeCascVarNegTOFSignal"     , &fTreeCascVarNegTOFSignal     , "fTreeCascVarNegTOFSignal/D");
    
    fTreeCascade->Branch("fTreeCascVarBachTOFBCid"      , &fTreeCascVarBachTOFBCid      , "fTreeCascVarBachTOFBCid/I");
    fTreeCascade->Branch("fTreeCascVarPosTOFBCid"       , &fTreeCascVarPosTOFBCid       , "fTreeCascVarPosTOFBCid/I");
    fTreeCascade->Branch("fTreeCascVarNegTOFBCid"       , &fTreeCascVarNegTOFBCid       , "fTreeCascVarNegTOFBCid/I");
    //Event info
    fTreeCascade->Branch("fTreeCascVarAmplitudeV0A"     , &fTreeCascVarAmplitudeV0A     , "fTreeCascVarAmplitudeV0A/D");
    fTreeCascade->Branch("fTreeCascVarAmplitudeV0C"     , &fTreeCascVarAmplitudeV0C     , "fTreeCascVarAmplitudeV0C/D");
    fTreeCascade->Branch("fTreeCascVarClosestNonEmptyBC", &fTreeCascVarClosestNonEmptyBC, "fTreeCascVarClosestNonEmptyBC/I");
    
    //Kink tagging
    if(!fRejectCascKink)
    {
        fTreeCascade->Branch("fTreeCascVarBachIsKink"       , &fTreeCascVarBachIsKink       , "fTreeCascVarBachIsKink/O");
        fTreeCascade->Branch("fTreeCascVarPosIsKink"        , &fTreeCascVarPosIsKink        , "fTreeCascVarPosIsKink/O");
        fTreeCascade->Branch("fTreeCascVarNegIsKink"        , &fTreeCascVarNegIsKink        , "fTreeCascVarNegIsKink/O");
    }
    
    //Cowboy/sailor studies
    fTreeCascade->Branch("fTreeCascVarIsCowboy"         , &fTreeCascVarIsCowboy         , "fTreeCascVarIsCowboy/O");
    fTreeCascade->Branch("fTreeCascVarCowboyness"       , &fTreeCascVarCowboyness       , "fTreeCascVarCowboyness/D");
    fTreeCascade->Branch("fTreeCascVarIsCascadeCowboy"  , &fTreeCascVarIsCascadeCowboy  , "fTreeCascVarIsCascadeCowboy/O");
    fTreeCascade->Branch("fTreeCascVarCascadeCowboyness", &fTreeCascVarCascadeCowboyness, "fTreeCascVarCascadeCowboyness/D");
    
    fTreeCascade->Branch("fTreeCascVarMagField"         , &fTreeCascVarMagField         , "fTreeCascVarMagField/D");
    fTreeCascade->Branch("fTreeCascVarRunNumber"        , &fTreeCascVarRunNumber        , "fTreeCascVarRunNumber/I");
    
    
    //Optional outputs
    //===========================================================================================
    //   Create V0 output tree
    //===========================================================================================
    fTreeV0 = new TTree("fTreeV0","V0Candidates");
    if( fSaveV0s )
    {
        //-----------BASIC-INFO---------------------------
        fTreeV0->Branch("fTreeV0VarMassAsK0s"           , &fTreeV0VarMassAsK0s              , "fTreeV0VarMassAsK0s/F");
        fTreeV0->Branch("fTreeV0VarMassAsLambda"        , &fTreeV0VarMassAsLambda           , "fTreeV0VarMassAsLambda/F");
        fTreeV0->Branch("fTreeV0VarMassAsAntiLambda"    , &fTreeV0VarMassAsAntiLambda       , "fTreeV0VarMassAsAntiLambda/F");
        fTreeV0->Branch("fTreeV0VarRapK0Short"          , &fTreeV0VarRapK0Short             , "fTreeV0VarRapK0Short/F");
        fTreeV0->Branch("fTreeV0VarRapLambda"           , &fTreeV0VarRapLambda              , "fTreeV0VarRapLambda/F");
        fTreeV0->Branch("fTreeV0VarPosEta"              , &fTreeV0VarPosEta                 , "fTreeV0VarPosEta/F");
        fTreeV0->Branch("fTreeV0VarNegEta"              , &fTreeV0VarNegEta                 , "fTreeV0VarNegEta/F");
        fTreeV0->Branch("fTreeV0VarPhi"                 , &fTreeV0VarPhi                    , "fTreeV0VarPhi/D");
        fTreeV0->Branch("fTreeV0VarTheta"               , &fTreeV0VarTheta                  , "fTreeV0VarTheta/D");
        fTreeV0->Branch("fTreeV0VarPtot"                , &fTreeV0VarPtot                   , "fTreeV0VarPtot/D");
        fTreeV0->Branch("fTreeV0VarPt"                  , &fTreeV0VarPt                     , "fTreeV0VarPt/D");
        //-----------INFO-FOR-CUTS------------------------
        fTreeV0->Branch("fTreeV0VarAlpha"               , &fTreeV0VarAlpha                  , "fTreeV0VarAlpha/F");
        fTreeV0->Branch("fTreeV0VarPtArm"               , &fTreeV0VarPtArm                  , "fTreeV0VarPtArm/F");
        fTreeV0->Branch("fTreeV0VarDCAV0Dau"            , &fTreeV0VarDCAV0Dau               , "fTreeV0VarDCAV0Dau/F");
        fTreeV0->Branch("fTreeV0VarDCAV0ToPV"           , &fTreeV0VarDCAV0ToPV              , "fTreeV0VarDCAV0ToPV/F");
        fTreeV0->Branch("fTreeV0VarDCAPosToPV"          , &fTreeV0VarDCAPosToPV             , "fTreeV0VarDCAPosToPV/F");
        fTreeV0->Branch("fTreeV0VarDCANegToPV"          , &fTreeV0VarDCANegToPV             , "fTreeV0VarDCANegToPV/F");
        fTreeV0->Branch("fTreeV0VarCosPA"               , &fTreeV0VarCosPA                  , "fTreeV0VarCosPA/F");
        fTreeV0->Branch("fTreeV0VarRadius"              , &fTreeV0VarRadius                 , "fTreeV0VarRadius/F");
        
        fTreeV0->Branch("fTreeV0VarLeastNbrCrossedRows" , &fTreeV0VarLeastNbrCrossedRows    , "fTreeV0VarLeastNbrCrossedRows/I");
        fTreeV0->Branch("fTreeV0VarLeastNbrClusters"    , &fTreeV0VarLeastNbrClusters       , "fTreeV0VarLeastNbrClusters/I");
        fTreeV0->Branch("fTreeV0VarLeastRatioCrossedRowsOverFindable"      , &fTreeV0VarLeastRatioCrossedRowsOverFindable           , "fTreeV0VarLeastRatioCrossedRowsOverFindable/D");
        fTreeV0->Branch("fTreeV0VarMaxChi2PerCluster"   , &fTreeV0VarMaxChi2PerCluster      , "fTreeV0VarMaxChi2PerCluster/F");
        fTreeV0->Branch("fTreeV0VarMinTrackLength"      , &fTreeV0VarMinTrackLength         , "fTreeV0VarMinTrackLength/F");
        //-----------DECAY-LENGTH-INFO--------------------
        fTreeV0->Branch("fTreeV0VarDistOverTotMom"      , &fTreeV0VarDistOverTotMom         , "fTreeV0VarDistOverTotMom/F");
        //---------------PID-TPC-INFO---------------------
        fTreeV0->Branch("fTreeV0VarPosNSigmaProton"     , &fTreeV0VarPosNSigmaProton        , "fTreeV0VarPosNSigmaProton/F");
        fTreeV0->Branch("fTreeV0VarPosNSigmaPion"       , &fTreeV0VarPosNSigmaPion          , "fTreeV0VarPosNSigmaPion/F");
        fTreeV0->Branch("fTreeV0VarNegNSigmaProton"     , &fTreeV0VarNegNSigmaProton        , "fTreeV0VarNegNSigmaProton/F");
        fTreeV0->Branch("fTreeV0VarNegNSigmaPion"       , &fTreeV0VarNegNSigmaPion          , "fTreeV0VarNegNSigmaPion/F");
        //---------------PID-TOF-INFO---------------------
        fTreeV0->Branch("fTreeV0VarPosTOFNSigmaProton"  , &fTreeV0VarPosTOFNSigmaProton     , "fTreeV0VarPosTOFNSigmaProton/F");
        fTreeV0->Branch("fTreeV0VarPosTOFNSigmaPion"    , &fTreeV0VarPosTOFNSigmaPion       , "fTreeV0VarPosTOFNSigmaPion/F");
        fTreeV0->Branch("fTreeV0VarNegTOFNSigmaProton"  , &fTreeV0VarNegTOFNSigmaProton     , "fTreeV0VarNegTOFNSigmaProton/F");
        fTreeV0->Branch("fTreeV0VarNegTOFNSigmaPion"    , &fTreeV0VarNegTOFNSigmaPion       , "fTreeV0VarNegTOFNSigmaPion/F");
        //---------------PID-ITS-INFO---------------------
        fTreeV0->Branch("fTreeV0VarPosITSNSigmaProton"  , &fTreeV0VarPosITSNSigmaProton     , "fTreeV0VarPosITSNSigmaProton/F");
        fTreeV0->Branch("fTreeV0VarPosITSNSigmaPion"    , &fTreeV0VarPosITSNSigmaPion       , "fTreeV0VarPosITSNSigmaPion/F");
        fTreeV0->Branch("fTreeV0VarNegITSNSigmaProton"  , &fTreeV0VarNegITSNSigmaProton     , "fTreeV0VarNegITSNSigmaProton/F");
        fTreeV0->Branch("fTreeV0VarNegITSNSigmaPion"    , &fTreeV0VarNegITSNSigmaPion       , "fTreeV0VarNegITSNSigmaPion/F");
        //---Raw TPC dEdx + PIDForTracking information----
        fTreeV0->Branch("fTreeV0VarPosdEdx"             , &fTreeV0VarPosdEdx                , "fTreeV0VarPosdEdx/D");
        fTreeV0->Branch("fTreeV0VarNegdEdx"             , &fTreeV0VarNegdEdx                , "fTreeV0VarNegdEdx/D");
        fTreeV0->Branch("fTreeV0VarPosPIDForTracking"   , &fTreeV0VarPosPIDForTracking      , "fTreeV0VarPosPIDForTracking/D");
        fTreeV0->Branch("fTreeV0VarNegPIDForTracking"   , &fTreeV0VarNegPIDForTracking      , "fTreeV0VarNegPIDForTracking/D");
        //------------------------------------------------
        fTreeV0->Branch("fTreeV0VarChi2V0"              , &fTreeV0VarChi2V0                 , "fTreeV0VarChi2V0/F");
        //------------------------------------------------
        fTreeV0->Branch("fTreeV0VarPosTrackStatus"      , &fTreeV0VarPosTrackStatus         , "fTreeV0VarPosTrackStatus/l");
        fTreeV0->Branch("fTreeV0VarNegTrackStatus"      , &fTreeV0VarNegTrackStatus         , "fTreeV0VarNegTrackStatus/l");
        //------------------------------------------------
        fTreeV0->Branch("fTreeV0VarPosDCAz"             , &fTreeV0VarPosDCAz                , "fTreeV0VarPosDCAz/D");
        fTreeV0->Branch("fTreeV0VarNegDCAz"             , &fTreeV0VarNegDCAz                , "fTreeV0VarNegDCAz/D");
        //------------FULL-MOMENTUM-INFO------------------
        fTreeV0->Branch("fTreeV0VarPosPx"               , &fTreeV0VarPosPx                  , "fTreeV0VarPosPx/D");
        fTreeV0->Branch("fTreeV0VarPosPy"               , &fTreeV0VarPosPy                  , "fTreeV0VarPosPy/D");
        fTreeV0->Branch("fTreeV0VarPosPz"               , &fTreeV0VarPosPz                  , "fTreeV0VarPosPz/D");
        fTreeV0->Branch("fTreeV0VarNegPx"               , &fTreeV0VarNegPx                  , "fTreeV0VarNegPx/D");
        fTreeV0->Branch("fTreeV0VarNegPy"               , &fTreeV0VarNegPy                  , "fTreeV0VarNegPy/D");
        fTreeV0->Branch("fTreeV0VarNegPz"               , &fTreeV0VarNegPz                  , "fTreeV0VarNegPz/D");
        //------------------------------------------------
        fTreeV0->Branch("fTreeV0VarV0DecayX"            , &fTreeV0VarV0DecayX               , "fTreeV0VarV0DecayX/D");
        fTreeV0->Branch("fTreeV0VarV0DecayY"            , &fTreeV0VarV0DecayY               , "fTreeV0VarV0DecayY/D");
        fTreeV0->Branch("fTreeV0VarV0DecayZ"            , &fTreeV0VarV0DecayZ               , "fTreeV0VarV0DecayZ/D");
        //---------------CLUSTER-INFO---------------------
        fTreeV0->Branch("fTreeV0VarPosITSClusters0"     , &fTreeV0VarPosITSClusters0        , "fTreeV0VarPosITSClusters0/O");
        fTreeV0->Branch("fTreeV0VarPosITSClusters1"     , &fTreeV0VarPosITSClusters1        , "fTreeV0VarPosITSClusters1/O");
        fTreeV0->Branch("fTreeV0VarPosITSClusters2"     , &fTreeV0VarPosITSClusters2        , "fTreeV0VarPosITSClusters2/O");
        fTreeV0->Branch("fTreeV0VarPosITSClusters3"     , &fTreeV0VarPosITSClusters3        , "fTreeV0VarPosITSClusters3/O");
        fTreeV0->Branch("fTreeV0VarPosITSClusters4"     , &fTreeV0VarPosITSClusters4        , "fTreeV0VarPosITSClusters4/O");
        fTreeV0->Branch("fTreeV0VarPosITSClusters5"     , &fTreeV0VarPosITSClusters5        , "fTreeV0VarPosITSClusters5/O");
        
        fTreeV0->Branch("fTreeV0VarNegITSClusters0"     , &fTreeV0VarNegITSClusters0        , "fTreeV0VarNegITSClusters0/O");
        fTreeV0->Branch("fTreeV0VarNegITSClusters1"     , &fTreeV0VarNegITSClusters1        , "fTreeV0VarNegITSClusters1/O");
        fTreeV0->Branch("fTreeV0VarNegITSClusters2"     , &fTreeV0VarNegITSClusters2        , "fTreeV0VarNegITSClusters2/O");
        fTreeV0->Branch("fTreeV0VarNegITSClusters3"     , &fTreeV0VarNegITSClusters3        , "fTreeV0VarNegITSClusters3/O");
        fTreeV0->Branch("fTreeV0VarNegITSClusters4"     , &fTreeV0VarNegITSClusters4        , "fTreeV0VarNegITSClusters4/O");
        fTreeV0->Branch("fTreeV0VarNegITSClusters5"     , &fTreeV0VarNegITSClusters5        , "fTreeV0VarNegITSClusters5/O");
        //------------------------------------------------
        fTreeV0->Branch("fTreeV0VarPosITSSharedClusters0"  , &fTreeV0VarPosITSSharedClusters0   , "fTreeV0VarPosITSSharedClusters0/O");
        fTreeV0->Branch("fTreeV0VarPosITSSharedClusters1"  , &fTreeV0VarPosITSSharedClusters1   , "fTreeV0VarPosITSSharedClusters1/O");
        fTreeV0->Branch("fTreeV0VarPosITSSharedClusters2"  , &fTreeV0VarPosITSSharedClusters2   , "fTreeV0VarPosITSSharedClusters2/O");
        fTreeV0->Branch("fTreeV0VarPosITSSharedClusters3"  , &fTreeV0VarPosITSSharedClusters3   , "fTreeV0VarPosITSSharedClusters3/O");
        fTreeV0->Branch("fTreeV0VarPosITSSharedClusters4"  , &fTreeV0VarPosITSSharedClusters4   , "fTreeV0VarPosITSSharedClusters4/O");
        fTreeV0->Branch("fTreeV0VarPosITSSharedClusters5"  , &fTreeV0VarPosITSSharedClusters5   , "fTreeV0VarPosITSSharedClusters5/O");
        
        fTreeV0->Branch("fTreeV0VarNegITSSharedClusters0"  , &fTreeV0VarNegITSSharedClusters0   , "fTreeV0VarNegITSSharedClusters0/O");
        fTreeV0->Branch("fTreeV0VarNegITSSharedClusters1"  , &fTreeV0VarNegITSSharedClusters1   , "fTreeV0VarNegITSSharedClusters1/O");
        fTreeV0->Branch("fTreeV0VarNegITSSharedClusters2"  , &fTreeV0VarNegITSSharedClusters2   , "fTreeV0VarNegITSSharedClusters2/O");
        fTreeV0->Branch("fTreeV0VarNegITSSharedClusters3"  , &fTreeV0VarNegITSSharedClusters3   , "fTreeV0VarNegITSSharedClusters3/O");
        fTreeV0->Branch("fTreeV0VarNegITSSharedClusters4"  , &fTreeV0VarNegITSSharedClusters4   , "fTreeV0VarNegITSSharedClusters4/O");
        fTreeV0->Branch("fTreeV0VarNegITSSharedClusters5"  , &fTreeV0VarNegITSSharedClusters5   , "fTreeV0VarNegITSSharedClusters5/O");
        //---------------OOB-PILEUP-INFO---------------------
        fTreeV0->Branch("fTreeV0VarNegTOFExpTDiff"      , &fTreeV0VarNegTOFExpTDiff         , "fTreeV0VarNegTOFExpTDiff/F");
        fTreeV0->Branch("fTreeV0VarPosTOFExpTDiff"      , &fTreeV0VarPosTOFExpTDiff         , "fTreeV0VarPosTOFExpTDiff/F");
        fTreeV0->Branch("fTreeV0VarNegTOFSignal"        , &fTreeV0VarNegTOFSignal           , "fTreeV0VarNegTOFSignal/F");
        fTreeV0->Branch("fTreeV0VarPosTOFSignal"        , &fTreeV0VarPosTOFSignal           , "fTreeV0VarPosTOFSignal/F");
        fTreeV0->Branch("fTreeV0VarNegTOFBCid"          , &fTreeV0VarNegTOFBCid             , "fTreeV0VarNegTOFBCid/I");
        fTreeV0->Branch("fTreeV0VarPosTOFBCid"          , &fTreeV0VarPosTOFBCid             , "fTreeV0VarPosTOFBCid/I");
        //Kink tagging
        if(!fRejectV0Kink)
        {
            fTreeV0->Branch("fTreeV0VarPosIsKink"        , &fTreeV0VarPosIsKink        , "fTreeV0VarPosIsKink/O");
            fTreeV0->Branch("fTreeV0VarNegIsKink"        , &fTreeV0VarNegIsKink        , "fTreeV0VarNegIsKink/O");
        }
        
        //Cowboy/sailor studies
        fTreeV0->Branch("fTreeV0VarIsCowboy"            , &fTreeV0VarIsCowboy               , "fTreeV0VarIsCowboy/O");
        
        fTreeV0->Branch("fTreeV0VarMagField"            , &fTreeV0VarMagField               , "fTreeV0VarMagField/F");
        fTreeV0->Branch("fTreeV0VarRunNumber"           , &fTreeV0VarRunNumber              , "fTreeV0VarRunNumber/I");
        fTreeV0->Branch("fTreeV0VarEventNumber"         , &fTreeV0VarEventNumber            , "fTreeV0VarEventNumber/l");
    }
    
    //===========================================================================================
    //   Create Resonance output tree
    //===========================================================================================
    fTreeRsn = new TTree("fTreeRsn","ResonanceCandidates");
    if( fSaveRsn )
    {
        fTreeRsn->Branch("fTreeRsnVarCutIDrsn"           , &fTreeRsnVarCutIDrsn           , "fTreeRsnVarCutIDrsn/I");
        fTreeRsn->Branch("fTreeRsnVarPx"                 , &fTreeRsnVarPx                 , "fTreeRsnVarPx/F");
        fTreeRsn->Branch("fTreeRsnVarPy"                 , &fTreeRsnVarPy                 , "fTreeRsnVarPy/F");
        fTreeRsn->Branch("fTreeRsnVarPz"                 , &fTreeRsnVarPz                 , "fTreeRsnVarPz/F");
        fTreeRsn->Branch("fTreeRsnVarInvMass"            , &fTreeRsnVarInvMass            , "fTreeRsnVarInvMass/D");
        fTreeRsn->Branch("fTreeRsnVarPassesOOBPileupCut" , &fTreeRsnVarPassesOOBPileupCut , "fTreeRsnVarPassesOOBPileupCut/O");
        fTreeRsn->Branch("fTreeRsnVarEventNumber"        , &fTreeRsnVarEventNumber        , "fTreeRsnVarEventNumber/l");
    }
    
    //===========================================================================================
    //   Create Mixed-events resonance output tree
    //===========================================================================================
    fTreeRsnBkg = new TTree("fTreeRsnBkg","MixedResonanceCandidates");
    if( fSaveRsn && ( fComputationType == "MIX" || fComputationType == "ROTATE1" || fComputationType == "ROTATE2" ) )
    {
        fTreeRsnBkg->Branch("fTreeRsnBkgVarCutIDrsn"           , &fTreeRsnBkgVarCutIDrsn           , "fTreeRsnBkgVarCutIDrsn/I");
        fTreeRsnBkg->Branch("fTreeRsnBkgVarPx"                 , &fTreeRsnBkgVarPx                 , "fTreeRsnBkgVarPx/F");
        fTreeRsnBkg->Branch("fTreeRsnBkgVarPy"                 , &fTreeRsnBkgVarPy                 , "fTreeRsnBkgVarPy/F");
        fTreeRsnBkg->Branch("fTreeRsnBkgVarPz"                 , &fTreeRsnBkgVarPz                 , "fTreeRsnBkgVarPz/F");
        fTreeRsnBkg->Branch("fTreeRsnBkgVarInvMass"            , &fTreeRsnBkgVarInvMass            , "fTreeRsnBkgVarInvMass/D");
        fTreeRsnBkg->Branch("fTreeRsnBkgVarPassesOOBPileupCut" , &fTreeRsnBkgVarPassesOOBPileupCut , "fTreeRsnBkgVarPassesOOBPileupCut/O");
        fTreeRsnBkg->Branch("fTreeRsnBkgVarEventNumber"   , &fTreeRsnBkgVarEventNumber   , "fTreeRsnBkgVarEventNumber/l");
        
        if( fComputationType == "MIX" )
            fTreeRsnBkg->Branch("fTreeRsnFoundMixEvts" , &fTreeRsnFoundMixEvts, "fTreeRsnFoundMixEvts/I");
        
        fDummyTree = new TTree("fDummyTree","ResonanceCandidates");
        fDummyTree->Branch("fRsnMiniEvent"          , &fRsnMiniEvent            , 16000, 0);
        fDummyTree->Branch("fDummyVarEventNumber"   , &fDummyVarEventNumber     , "fDummyVarEventNumber/l");
    }
    
    //===========================================================================================
    //   Create Primary tracks output tree
    //===========================================================================================
    fTreePrimTrack = new TTree("fTreePrimTrack","PrimaryTrackCandidates");
    if( fSavePrimaries )
    {
        //-----------BASIC-INFO---------------------------
        fTreePrimTrack->Branch("fTreePrimVarCharge"             , &fTreePrimVarCharge               , "fTreePrimVarCharge/I");
        fTreePrimTrack->Branch("fTreePrimVarRapPion"            , &fTreePrimVarRapPion              , "fTreePrimVarRapPion/F");
        fTreePrimTrack->Branch("fTreePrimVarRapProton"          , &fTreePrimVarRapProton            , "fTreePrimVarRapProton/F");
        fTreePrimTrack->Branch("fTreePrimVarRapKaon"            , &fTreePrimVarRapKaon              , "fTreePrimVarRapKaon/F");
        fTreePrimTrack->Branch("fTreePrimVarEta"                , &fTreePrimVarEta                  , "fTreePrimVarEta/F");
        fTreePrimTrack->Branch("fTreePrimVarTheta"              , &fTreePrimVarTheta                , "fTreePrimVarTheta/F");
        fTreePrimTrack->Branch("fTreePrimVarPhi"                , &fTreePrimVarPhi                  , "fTreePrimVarPhi/F");
        fTreePrimTrack->Branch("fTreePrimVarPtot"               , &fTreePrimVarPtot                 , "fTreePrimVarPtot/F");
        fTreePrimTrack->Branch("fTreePrimVarPt"                 , &fTreePrimVarPt                   , "fTreePrimVarPt/F");
        
        fTreePrimTrack->Branch("fTreePrimVarDCAxyToPV"          , &fTreePrimVarDCAxyToPV            , "fTreePrimVarDCAxyToPV/D");
        fTreePrimTrack->Branch("fTreePrimVarDCAzToPV"           , &fTreePrimVarDCAzToPV             , "fTreePrimVarDCAzToPV/D");
        
        fTreePrimTrack->Branch("fTreePrimVarNbrCrossedRows"     , &fTreePrimVarNbrCrossedRows       , "fTreePrimVarNbrCrossedRows/I");
        fTreePrimTrack->Branch("fTreePrimVarNbrClusters"        , &fTreePrimVarNbrClusters          , "fTreePrimVarNbrClusters/I");
        fTreePrimTrack->Branch("fTreePrimVarRatioCrossedRowsOverFindable" , &fTreePrimVarRatioCrossedRowsOverFindable, "fTreePrimVarRatioCrossedRowsOverFindable/D");
        fTreePrimTrack->Branch("fTreePrimVarNbrCrossedRowsOverLength" , &fTreePrimVarNbrCrossedRowsOverLength, "fTreePrimVarNbrCrossedRowsOverLength/D");
        fTreePrimTrack->Branch("fTreePrimVarFractionSharedTPCClusters", &fTreePrimVarFractionSharedTPCClusters, "fTreePrimVarFractionSharedTPCClusters/D");
        fTreePrimTrack->Branch("fTreePrimVarITSChi2PerCluster"  , &fTreePrimVarITSChi2PerCluster    , "fTreePrimVarITSChi2PerCluster/D");
        fTreePrimTrack->Branch("fTreePrimVarTPCChi2PerCluster"  , &fTreePrimVarTPCChi2PerCluster    , "fTreePrimVarTPCChi2PerCluster/D");
        fTreePrimTrack->Branch("fTreePrimVarTrackLength"        , &fTreePrimVarTrackLength          , "fTreePrimVarTrackLength/D");
        //---------------PID-TPC-INFO---------------------
        fTreePrimTrack->Branch("fTreePrimVarNSigmaPion"         , &fTreePrimVarNSigmaPion           , "fTreePrimVarNSigmaPion/F");
        fTreePrimTrack->Branch("fTreePrimVarNSigmaKaon"         , &fTreePrimVarNSigmaKaon           , "fTreePrimVarNSigmaKaon/F");
        fTreePrimTrack->Branch("fTreePrimVarNSigmaProton"       , &fTreePrimVarNSigmaProton         , "fTreePrimVarNSigmaProton/F");
        //---------------PID-TOF-INFO---------------------
        fTreePrimTrack->Branch("fTreePrimVarTOFNSigmaPion"      , &fTreePrimVarTOFNSigmaPion        , "fTreePrimVarTOFNSigmaPion/F");
        fTreePrimTrack->Branch("fTreePrimVarTOFNSigmaKaon"      , &fTreePrimVarTOFNSigmaKaon        , "fTreePrimVarTOFNSigmaKaon/F");
        fTreePrimTrack->Branch("fTreePrimVarTOFNSigmaProton"    , &fTreePrimVarTOFNSigmaProton      , "fTreePrimVarTOFNSigmaProton/F");
        //---------------PID-ITS-INFO---------------------
        fTreePrimTrack->Branch("fTreePrimVarITSNSigmaPion"      , &fTreePrimVarITSNSigmaPion        , "fTreePrimVarITSNSigmaPion/F");
        fTreePrimTrack->Branch("fTreePrimVarITSNSigmaKaon"      , &fTreePrimVarITSNSigmaKaon        , "fTreePrimVarITSNSigmaKaon/F");
        fTreePrimTrack->Branch("fTreePrimVarITSNSigmaProton"    , &fTreePrimVarITSNSigmaProton      , "fTreePrimVarITSNSigmaProton/F");
        //---Raw TPC dEdx + PIDForTracking information----    
        fTreePrimTrack->Branch("fTreePrimVardEdx"               , &fTreePrimVardEdx                 , "fTreePrimVardEdx/D");
        fTreePrimTrack->Branch("fTreePrimVarPIDForTracking"     , &fTreePrimVarPIDForTracking       , "fTreePrimVarPIDForTracking/D");
        //------------------------------------------------
        fTreePrimTrack->Branch("fTreePrimVarTrackStatus"        , &fTreePrimVarTrackStatus          , "fTreePrimVarTrackStatus/l");
        //------------FULL-MOMENTUM-INFO------------------
        fTreePrimTrack->Branch("fTreePrimVarPx"                 , &fTreePrimVarPx                   , "fTreePrimVarPx/D");
        fTreePrimTrack->Branch("fTreePrimVarPy"                 , &fTreePrimVarPy                   , "fTreePrimVarPy/D");
        fTreePrimTrack->Branch("fTreePrimVarPz"                 , &fTreePrimVarPz                   , "fTreePrimVarPz/D");
        //---------------CLUSTER-INFO---------------------
        fTreePrimTrack->Branch("fTreePrimVarITSClusters0"       , &fTreePrimVarITSClusters0         , "fTreePrimVarITSClusters0/O");
        fTreePrimTrack->Branch("fTreePrimVarITSClusters1"       , &fTreePrimVarITSClusters1         , "fTreePrimVarITSClusters1/O");
        fTreePrimTrack->Branch("fTreePrimVarITSClusters2"       , &fTreePrimVarITSClusters2         , "fTreePrimVarITSClusters2/O");
        fTreePrimTrack->Branch("fTreePrimVarITSClusters3"       , &fTreePrimVarITSClusters3         , "fTreePrimVarITSClusters3/O");
        fTreePrimTrack->Branch("fTreePrimVarITSClusters4"       , &fTreePrimVarITSClusters4         , "fTreePrimVarITSClusters4/O");
        fTreePrimTrack->Branch("fTreePrimVarITSClusters5"       , &fTreePrimVarITSClusters5         , "fTreePrimVarITSClusters5/O");
        //------------------------------------------------
        fTreePrimTrack->Branch("fTreePrimVarITSSharedClusters0" , &fTreePrimVarITSSharedClusters0   , "fTreePrimVarITSSharedClusters0/O");
        fTreePrimTrack->Branch("fTreePrimVarITSSharedClusters1" , &fTreePrimVarITSSharedClusters1   , "fTreePrimVarITSSharedClusters1/O");
        fTreePrimTrack->Branch("fTreePrimVarITSSharedClusters2" , &fTreePrimVarITSSharedClusters2   , "fTreePrimVarITSSharedClusters2/O");
        fTreePrimTrack->Branch("fTreePrimVarITSSharedClusters3" , &fTreePrimVarITSSharedClusters3   , "fTreePrimVarITSSharedClusters3/O");
        fTreePrimTrack->Branch("fTreePrimVarITSSharedClusters4" , &fTreePrimVarITSSharedClusters4   , "fTreePrimVarITSSharedClusters4/O");
        fTreePrimTrack->Branch("fTreePrimVarITSSharedClusters5" , &fTreePrimVarITSSharedClusters5   , "fTreePrimVarITSSharedClusters5/O");
        //---------------OOB-PILEUP-INFO---------------------
        fTreePrimTrack->Branch("fTreePrimVarTOFExpTDiff"        , &fTreePrimVarTOFExpTDiff          , "fTreePrimVarTOFExpTDiff/D");
        fTreePrimTrack->Branch("fTreePrimVarTOFSignal"          , &fTreePrimVarTOFSignal            , "fTreePrimVarTOFSignal/D");
        fTreePrimTrack->Branch("fTreePrimVarTOFBCid"            , &fTreePrimVarTOFBCid              , "fTreePrimVarTOFBCid/I");
        
        if(!fPrimariesSaveAddConfig || (fPrimTrackCuts && fPrimTrackCuts->GetAcceptKinkDaughters()) )
            fTreePrimTrack->Branch("fTreePrimVarIsKink"         , &fTreePrimVarIsKink               , "fTreePrimVarIsKink/O");
        
        fTreePrimTrack->Branch("fTreePrimVarRunNumber"          , &fTreePrimVarRunNumber            , "fTreePrimVarRunNumber/I");
        fTreePrimTrack->Branch("fTreePrimVarEventNumber"        , &fTreePrimVarEventNumber          , "fTreePrimVarEventNumber/l");
        
    }
    
    
    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    // Rsn event and mini-event
    if(! fRsnEvent ) fRsnEvent = new AliRsnEvent();
    
    if(! fRsnMiniEvent ) fRsnMiniEvent = new AliRsnMiniEvent();
    
    //Analysis Utils
    if(! fUtils ) fUtils = new AliAnalysisUtils();
    
    fOutputList = new TList();   
    fOutputList->SetOwner(kTRUE);
    
    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
        fOutputList->Add(fHistEventCounter);
    }
    
    PostData(1, fOutputList); 
    PostData(2, fTreeCascade);
    PostData(3, fTreeV0);
    PostData(4, fTreeRsn);
    PostData(5, fTreeRsnBkg);
    PostData(6, fTreePrimTrack);        
}
//_____________________________________________________________________________
void AliAnalysisTaskStrangeCascadesTriggerAODRun2::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    AliAODEvent *lAODevent = 0x0;
    // Connect to the InputEvent
    lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
    if (!lAODevent) {
        AliWarning("ERROR: lAODevent not available \n");
        return;
    }
    fHistEventCounter->Fill(0.5);
    
    Double_t lMagField = lAODevent->GetMagneticField( );
    
    //------------------------------------------------
    // Retrieving IR info for OOB Pileup rejection
    //------------------------------------------------
    Int_t lClosestNonEmptyBC = 10*3564; // start with an isolated event
    
    //------------------------------------------------
    // Primary Vertex Requirements Section:
    //  ---> pp: has vertex, |z|<10cm
    //------------------------------------------------
    
    const AliAODVertex *lPrimaryBestAODVtx     = lAODevent->GetPrimaryVertex();
    //const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    //const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();
    
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );
    
    //Optional cut on the primary vertex
    if ( TMath::Sqrt( TMath::Power(lBestPrimaryVtxPos[0],2)+TMath::Power(lBestPrimaryVtxPos[1],2)) > fkMaxPVR2D ) return;
    if ( TMath::Abs( lBestPrimaryVtxPos[2] ) > fkMaxPVZ ) return;
    
    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------
    
    //Get VZERO Information for multiplicity later
    AliAODVZERO* AODVZERO = (AliAODVZERO*) lAODevent->GetVZEROData();
    if (!AODVZERO) {
        AliError("AliAODVZERO not available");
        return;
    }
    
    //--------------VZERO-MULTIPLICITY----------------
    Double_t vzeroA = AODVZERO->GetMTotV0A();
    Double_t vzeroC = AODVZERO->GetMTotV0C();
    //Double_t vzeroAEq = AliVertexingHFUtils::GetVZEROAEqualizedMultiplicity(aod);
    //Double_t vzeroCEq = AliVertexingHFUtils::GetVZEROCEqualizedMultiplicity(aod);
    Int_t fnV0A = static_cast<Int_t>(vzeroA);
    Int_t fnV0M = static_cast<Int_t>(vzeroA + vzeroC);
    //fnV0MEq = static_cast<Int_t>(vzeroAEq + vzeroCEq);
    Int_t fnV0MCorr = -1;
    //fnV0MEqCorr = -1;
    fnV0MCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroA, lBestPrimaryVtxPos[2]) + AliESDUtils::GetCorrV0C(vzeroC, lBestPrimaryVtxPos[2]));
    //    fnV0MEqCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(vzeroAEq, vtx->GetZ()) + AliESDUtils::GetCorrV0C(vzeroCEq, vtx->GetZ()));
    
    // multiplicity percentiles
    AliMultSelection *MultSelection = (AliMultSelection*) lAODevent -> FindListObject("MultSelection");
    Double_t fPercV0M = MultSelection ? MultSelection->GetMultiplicityPercentile("V0M") : -1.;
    // multiplicity from mult selection task
    AliMultEstimator* MultEst = MultSelection ? MultSelection->GetEstimator("V0M") : nullptr;
    Double_t fMultV0M = MultEst ? MultEst->GetValue() : -1.;
    
    //--------------SPD-MULTIPLICITY-----------------
    //n tracklets
    AliAODTracklets* tracklets= lAODevent->GetTracklets();
    Int_t nTr = tracklets->GetNumberOfTracklets();
    Int_t NbrTrackletEta1=0;
    for(Int_t iTr=0; iTr<nTr; iTr++){
        Double_t theta=tracklets->GetTheta(iTr);
        Double_t Eta=-TMath::Log(TMath::Tan(theta/2.));
        if( TMath::Abs(Eta) < 1.0 ) NbrTrackletEta1++;//count at central rapidity
    }
    Int_t fnTracklets = NbrTrackletEta1;
    Double_t fnTrackletsCorr = -1.;

    /*TProfile *estimatorAvg = fMultEstimatorAvg[GetPeriod(lAODevent)];
    if (estimatorAvg)
        fnTrackletsCorr = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg, countTreta1, vtx->GetZ(), fRefMult));*/
    /*TProfile *estimatorAvgSHM = fMultEstimatorAvgSHM[GetPeriod(aod)];
    if (fCorrNtrVtx && estimatorAvgSHM)
        fnTrackletsCorrSHM = static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvgSHM, countTreta1, vtx->GetZ(), fRefMultSHM));
    */
    
    fTreeCascVarIsIncompleteDAQ     = !MultSelection->GetThisEventIsNotIncompleteDAQ();
    fTreeCascVarSPDPileupFlag       = fUtils->IsPileUpSPD(InputEvent());
    fTreeCascVarMVPileupFlag        = !MultSelection->GetThisEventIsNotPileupMV();
    fTreeCascVarOOBPileupFlag       = fUtils->IsOutOfBunchPileUp(InputEvent());
    fTreeCascVarCentrality          = fPercV0M;
    fTreeCascVarVZEROMultSel        = fMultV0M;
    fTreeCascVarVZEROMultSig        = fnV0M;
    fTreeCascVarVZEROMultSigCorr    = fnV0MCorr;
    fTreeCascVarSPDMult             = fnTracklets;
    
    //Implementation to do trigger selection a posteriori
    UInt_t fEvSel_TriggerMask =
        ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected(); //for full checks later
        
    fTreeCascVar_TriggerMask = fEvSel_TriggerMask;
    
    //VZERO info
    AliMultEstimator *fEstV0A = 0x0, *fEstV0C = 0x0;
    fEstV0A = (AliMultEstimator*)MultSelection->GetEstimator("V0A");
    fEstV0C = (AliMultEstimator*)MultSelection->GetEstimator("V0C");
    Double_t lAmplitudeV0A = -1.;
    Double_t lAmplitudeV0C = -1.;
    if ( fEstV0A ) lAmplitudeV0A = fEstV0A->GetValue();
    if ( fEstV0C ) lAmplitudeV0C = fEstV0C->GetValue();
    
    fHistEventCounter->Fill(1.5);
    
    Bool_t IsEvtWithCascade = kFALSE;
    
    //------------------------------------------------
    // MAIN CASCADE LOOP STARTS HERE
    //------------------------------------------------

    Long_t ncascades = 0;
    ncascades = lAODevent->GetNumberOfCascades();
    for (Int_t iCasc = 0; iCasc < ncascades; iCasc++){
        
        //------------------------------------------------
        // Initializations
        //------------------------------------------------
        Bool_t   lOnFlyStatus           = kFALSE; // if kTRUE, then this V0 is recontructed "on fly" during the tracking
        Double_t lPosCasc[3]            = {0., 0., 0.};// Position of cascade decay vtx
        Double_t lPosV0[3]              = {0.,0.,0.}; // Position of V0 coming from cascade
        
        Double_t lMassAsXiMinus         = 0.;
        Double_t lMassAsXiPlus          = 0.;
        Double_t lMassAsOmegaMinus      = 0.;
        Double_t lMassAsOmegaPlus       = 0.;
        Double_t lV0MassAsLambda        = 0.;
        Double_t lV0MassAsAntiLambda    = 0.;
         
        AliAODcascade *cascade = lAODevent->GetCascade(iCasc);
        if (!cascade) continue;
        
        lPosCasc[0] = cascade->DecayVertexXiX();
        lPosCasc[1] = cascade->DecayVertexXiY();
        lPosCasc[2] = cascade->DecayVertexXiZ();
        
        lPosV0[0] = cascade->DecayVertexV0X();
        lPosV0[1] = cascade->DecayVertexV0Y();
        lPosV0[2] = cascade->DecayVertexV0Z();
        
        //Gather tracks informations
        AliAODTrack *bTrack = dynamic_cast<AliAODTrack*>( cascade->GetDecayVertexXi()->GetDaughter(0) );
        AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>( cascade->GetDaughter(0) );
        AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>( cascade->GetDaughter(1) );
        
        if (!bTrack || !pTrack || !nTrack ) {
            AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
            continue;
        }
        
        //Findable clusters > 0 condition
        if( bTrack->GetTPCNclsF()<=0 || pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;
        
        //Test different mass hypothesis : Xi-/Xi+ and Omega-/Omega+
        if(bTrack->Charge() < 0) lMassAsXiMinus     = cascade->MassXi();
        if(bTrack->Charge() > 0) lMassAsXiPlus      = cascade->MassXi();
        if(bTrack->Charge() < 0) lMassAsOmegaMinus  = cascade->MassOmega();
        if(bTrack->Charge() > 0) lMassAsOmegaPlus   = cascade->MassOmega();
        
        lV0MassAsLambda         = cascade->MassLambda();
        lV0MassAsAntiLambda     = cascade->MassAntiLambda();
        
        Double_t lBachPx = 0. ; Double_t lBachPy = 0.; Double_t lBachPz = 0.;
        Double_t lPosPx = 0.; Double_t lPosPy = 0.; Double_t lPosPz = 0.;
        Double_t lNegPx = 0.; Double_t lNegPy = 0.; Double_t lNegPz = 0.;
        lBachPx = bTrack->Px() ; lBachPy = bTrack->Py() ; lBachPz = bTrack->Pz() ;
        lPosPx  = pTrack->Px() ; lPosPy  = pTrack->Py() ; lPosPz  = pTrack->Pz() ;
        lNegPx  = nTrack->Px() ; lNegPy  = nTrack->Py() ; lNegPz  = nTrack->Pz() ;
        
        Double_t lV0Px = lPosPx + lNegPx;
        Double_t lV0Py = lPosPy + lNegPy;
        Double_t lV0Pz = lPosPz + lNegPz;

        //------------------------------------------------
        // Calculation of the variables related to casc.
        //------------------------------------------------
        lOnFlyStatus                    = cascade->GetOnFlyStatus();
        fTreeCascVarChi2Cascade         = cascade->Chi2Xi();
        fTreeCascVarChi2CascadePerNDF   = ((AliAODVertex*)cascade->GetDecayVertexXi())->GetChi2perNDF();
        fTreeCascVarChi2V0              = cascade->Chi2V0();
        
        fTreeCascVarCharge      = cascade->ChargeXi();
        if(fTreeCascVarCharge > 0){
            fTreeCascVarV0Mass = lV0MassAsAntiLambda;
            fTreeCascVarMassAsXi = lMassAsXiPlus;
            fTreeCascVarMassAsOmega = lMassAsOmegaPlus;
        }
        if(fTreeCascVarCharge < 0){
            fTreeCascVarV0Mass = lV0MassAsLambda;
            fTreeCascVarMassAsXi = lMassAsXiMinus;
            fTreeCascVarMassAsOmega = lMassAsOmegaMinus;
        }
        fTreeCascVarV0MassAsLambda      = lV0MassAsLambda;
        fTreeCascVarV0MassAsAntiLambda  = lV0MassAsAntiLambda;
        
        fTreeCascVarPtot        = cascade->P();
        fTreeCascVarV0Ptot      = TMath::Sqrt( TMath::Power(lPosPx+lNegPx,2) 
                                             + TMath::Power(lPosPy+lNegPy,2) 
                                             + TMath::Power(lPosPz+lNegPz,2) );
        fTreeCascVarPt          = cascade->Pt();
        fTreeCascVarV0Pt        = TMath::Sqrt( TMath::Power(lPosPx+lNegPx,2) 
                                             + TMath::Power(lPosPy+lNegPy,2) );
        fTreeCascVarRapXi       = cascade->RapXi();
        fTreeCascVarRapOmega    = cascade->RapOmega();
        fTreeCascVarBachEta     = bTrack->Eta();
        fTreeCascVarPosEta      = pTrack->Eta();
        fTreeCascVarNegEta      = nTrack->Eta();
        
        fTreeCascVarPhi         = MyPhi(lBachPx+lPosPx+lNegPx, lBachPy+lPosPy+lNegPy);
        fTreeCascVarTheta       = MyTheta(fTreeCascVarPt, lBachPz+lPosPz+lNegPz);
        
        fTreeCascVarAlpha   = cascade->AlphaXi();
        fTreeCascVarPtArm   = cascade->PtArmXi();
        fTreeCascVarAlphaV0 = cascade->AlphaV0();
        fTreeCascVarPtArmV0 = cascade->PtArmV0();
        
        fTreeCascVarDCACascDau  = cascade->DcaXiDaughters();
        fTreeCascVarDCAV0Dau    = cascade->DcaV0Daughters();
        fTreeCascVarDCAV0ToPV   = cascade->DcaV0ToPrimVertex();
        
        fTreeCascVarDCABachToPV = cascade->DcaBachToPrimVertex();
        fTreeCascVarDCAPosToPV  = cascade->DcaPosToPrimVertex();
        fTreeCascVarDCANegToPV  = cascade->DcaNegToPrimVertex();
        
        
        fTreeCascVarCascRadius  = TMath::Sqrt( lPosCasc[0]*lPosCasc[0]  +  lPosCasc[1]*lPosCasc[1] );
        fTreeCascVarV0Radius    = TMath::Sqrt( lPosV0[0]*lPosV0[0]  +  lPosV0[1]*lPosV0[1] );
        fTreeCascVarCascCosPA   = cascade->CosPointingAngleXi(lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
        fTreeCascVarV0CosPA     = cascade->CosPointingAngle(lBestPrimaryVtxPos);
        
        //FIX ME
        //if( bTrack->Charge() < 0 ) fTreeCascVarWrongCosPA = GetCosPA( bTrack , pTrack, lAODevent);
        //if( bTrack->Charge() > 0 ) fTreeCascVarWrongCosPA = GetCosPA( bTrack , nTrack, lAODevent);
        fTreeCascVarWrongCosPA         = cascade->BachBaryonCosPA();
        fTreeCascVarCascCosPASpecial   = cascade->CosPointingAngle(lPosCasc);
        
        //Reject on-the-fly V0s
        if( lOnFlyStatus ) continue;
            
        // Filter like-sign V0
        if ( pTrack->GetSign() == nTrack->GetSign() ) continue;
        
        //________________________________________________________________________
        // Track quality cuts
        //Least Nbr of Crossed Rows
        Int_t lBachNbrCrossedRows = bTrack->GetTPCClusterInfo(2,1);
        Int_t lPosNbrCrossedRows  = pTrack->GetTPCClusterInfo(2,1);
        Int_t lNegNbrCrossedRows  = nTrack->GetTPCClusterInfo(2,1);
        
        Int_t lLeastNbrCrossedRows = (Int_t) lPosNbrCrossedRows;
        if( lNegNbrCrossedRows < lLeastNbrCrossedRows ) lLeastNbrCrossedRows = (Int_t) lNegNbrCrossedRows;
        if( lBachNbrCrossedRows < lLeastNbrCrossedRows ) lLeastNbrCrossedRows = (Int_t) lBachNbrCrossedRows;
        
        fTreeCascVarLeastNbrCrossedRows = lLeastNbrCrossedRows;
        
        //Compute ratio Crossed Rows / Findable clusters
        //Note: above test avoids division by zero!
        Double_t lBachTrackCrossedRowsOverFindable = lBachNbrCrossedRows / ((Double_t)(bTrack->GetTPCNclsF()));
        Double_t lPosTrackCrossedRowsOverFindable = lPosNbrCrossedRows / ((Double_t)(pTrack->GetTPCNclsF()));
        Double_t lNegTrackCrossedRowsOverFindable = lNegNbrCrossedRows / ((Double_t)(nTrack->GetTPCNclsF()));
        
        fTreeCascVarLeastRatioCrossedRowsOverFindable = lBachTrackCrossedRowsOverFindable;
        if( lPosTrackCrossedRowsOverFindable < fTreeCascVarLeastRatioCrossedRowsOverFindable ) 
            fTreeCascVarLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
        if( lNegTrackCrossedRowsOverFindable < fTreeCascVarLeastRatioCrossedRowsOverFindable ) 
            fTreeCascVarLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
        
        // Least Nbr of clusters        
        Int_t lPosTPCClusters        = 0.;
        Int_t lNegTPCClusters        = 0.;
        Int_t lBachTPCClusters       = 0.;

        lPosTPCClusters   = pTrack->GetTPCNcls();
        lNegTPCClusters   = nTrack->GetTPCNcls();
        lBachTPCClusters  = bTrack->GetTPCNcls();
        
        Int_t lLeastNbrOfClusters    = lBachTPCClusters;
        if( lPosTPCClusters < lLeastNbrOfClusters ) lLeastNbrOfClusters = lPosTPCClusters;
        if( lNegTPCClusters < lLeastNbrOfClusters ) lLeastNbrOfClusters = lNegTPCClusters;
        
        fTreeCascVarLeastNbrClusters = lLeastNbrOfClusters;
        
        //Min track length
        Double_t lPosTrackLength        = 0.;
        Double_t lNegTrackLength        = 0.;
        Double_t lBachTrackLength       = 0.;
        
        lBachTrackLength   = GetLengthInActiveZone(bTrack, 2.0, 220.0, lMagField);
        lPosTrackLength    = GetLengthInActiveZone(pTrack, 2.0, 220.0, lMagField);
        lNegTrackLength    = GetLengthInActiveZone(nTrack, 2.0, 220.0, lMagField);
        
        Double_t lSmallestTrackLength = lBachTrackLength;
        if ( lPosTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lPosTrackLength;
        if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
        
        fTreeCascVarMinTrackLength = lSmallestTrackLength;
        
        //Nbr Of Crossed Rows Over Length
        Double_t lBachTrackNbrCrOverLength = bTrack->GetTPCClusterInfo(2,1)/(lBachTrackLength-TMath::Max(fTreeCascVarCascRadius-85.,0.) + 1e-5);
        Double_t lPosTrackNbrCrOverLength = pTrack->GetTPCClusterInfo(2,1)/(lPosTrackLength-TMath::Max(fTreeCascVarV0Radius-85.,0.) + 1e-5);
        Double_t lNegTrackNbrCrOverLength = nTrack->GetTPCClusterInfo(2,1)/(lNegTrackLength-TMath::Max(fTreeCascVarV0Radius-85.,0.) + 1e-5);

        Double_t lLeastNbrCrOverLength = (Double_t) lPosTrackNbrCrOverLength;
        if( lNegTrackNbrCrOverLength < lLeastNbrCrOverLength ) lLeastNbrCrOverLength = (Float_t) lNegTrackNbrCrOverLength;
        if( lBachTrackNbrCrOverLength < lLeastNbrCrOverLength ) lLeastNbrCrOverLength = (Float_t) lBachTrackNbrCrOverLength;
        
        fTreeCascVarNbrCrossedRowsOverLength = lLeastNbrCrOverLength;
        
        //Max Chi2 per cluster
        Double_t lBiggestChi2PerCluster = 0.;
        Double_t lPosChi2PerCluster = pTrack->GetTPCchi2() / ((Float_t) lPosTPCClusters);
        Double_t lNegChi2PerCluster = nTrack->GetTPCchi2() / ((Float_t) lNegTPCClusters);
        Double_t lBachChi2PerCluster = bTrack->GetTPCchi2() / ((Float_t) lBachTPCClusters);
        
        if( lPosChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lPosChi2PerCluster;
        if( lNegChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lNegChi2PerCluster;
        if( lBachChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lBachChi2PerCluster;
        
        fTreeCascVarMaxChi2PerCluster = lBiggestChi2PerCluster;
        
        //________________________________________________________________________
        // Decay length info
        //Cascade
        fTreeCascVarDistOverTotMom  = TMath::Sqrt( TMath::Power( lPosCasc[0] - lBestPrimaryVtxPos[0] , 2) 
                                                 + TMath::Power( lPosCasc[1] - lBestPrimaryVtxPos[1] , 2) 
                                                 + TMath::Power( lPosCasc[2] - lBestPrimaryVtxPos[2] , 2) );
        fTreeCascVarDistOverTotMom /= fTreeCascVarPtot;
        
        //V0
        fTreeCascVarV0DistOverTotMom    = TMath::Sqrt( TMath::Power( lPosV0[0] - lPosCasc[0] , 2) 
                                          + TMath::Power( lPosV0[1] - lPosCasc[1] , 2) 
                                          + TMath::Power( lPosV0[2] - lPosCasc[2] , 2) );
        fTreeCascVarV0DistOverTotMom /= fTreeCascVarV0Ptot;
        
        
        //------------------------------------------------
        // TPC dEdx information
        //------------------------------------------------
        fTreeCascVarBachNSigmaPion  = fPIDResponse->NumberOfSigmasTPC( bTrack, AliPID::kPion );
        fTreeCascVarBachNSigmaKaon  = fPIDResponse->NumberOfSigmasTPC( bTrack, AliPID::kKaon );
        fTreeCascVarPosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
        fTreeCascVarPosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
        fTreeCascVarNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion   );
        fTreeCascVarNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
        
        //------------------------------------------------
        // ITS info 
        //------------------------------------------------
        fTreeCascVarBachITSNSigmaPion  = fPIDResponse->NumberOfSigmasITS( bTrack, AliPID::kPion );
        fTreeCascVarBachITSNSigmaKaon  = fPIDResponse->NumberOfSigmasITS( bTrack, AliPID::kKaon );
        fTreeCascVarPosITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( pTrack, AliPID::kPion );
        fTreeCascVarPosITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( pTrack, AliPID::kProton );
        fTreeCascVarNegITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( nTrack, AliPID::kPion   );
        fTreeCascVarNegITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( nTrack, AliPID::kProton );
        
        //------------------------------------------------
        // TOF info 
        //------------------------------------------------
        fTreeCascVarBachTOFNSigmaPion  = fPIDResponse->NumberOfSigmasTOF( bTrack, AliPID::kPion );
        fTreeCascVarBachTOFNSigmaKaon  = fPIDResponse->NumberOfSigmasTOF( bTrack, AliPID::kKaon );
        fTreeCascVarPosTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( pTrack, AliPID::kPion );
        fTreeCascVarPosTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( pTrack, AliPID::kProton );
        fTreeCascVarNegTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( nTrack, AliPID::kPion   );
        fTreeCascVarNegTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( nTrack, AliPID::kProton );
        
        //------------------------------------------------
        // Raw TPC dEdx + PIDForTracking information
        //------------------------------------------------
        
        //Step 2: Acquire TPC Signals
        fTreeCascVarBachdEdx    = bTrack->GetTPCsignal();
        fTreeCascVarPosdEdx     = pTrack->GetTPCsignal();
        fTreeCascVarNegdEdx     = nTrack->GetTPCsignal();
        
        
        //Step 3: Acquire PID For Tracking
        fTreeCascVarBachPIDForTracking  = bTrack->GetPIDForTracking();
        fTreeCascVarPosPIDForTracking   = pTrack->GetPIDForTracking();
        fTreeCascVarNegPIDForTracking   = nTrack->GetPIDForTracking();
        
        
        //________________________________________________________________________
        // Track status
        fTreeCascVarBachTrackStatus     = bTrack->GetStatus();
        fTreeCascVarPosTrackStatus      = pTrack->GetStatus();
        fTreeCascVarNegTrackStatus      = nTrack->GetStatus();
        
        if ((fTreeCascVarPosTrackStatus & AliAODTrack::kTPCrefit)  == 0) {
            AliDebug(1, "Pb / V0 Pos. track has no TPCrefit ... continue!");
            continue;
        }
        if ((fTreeCascVarNegTrackStatus & AliAODTrack::kTPCrefit)  == 0) {
            AliDebug(1, "Pb / V0 Neg. track has no TPCrefit ... continue!");
            continue;
        }
        if ((fTreeCascVarBachTrackStatus & AliAODTrack::kTPCrefit) == 0) {
            AliDebug(1, "Pb / Bach.   track has no TPCrefit ... continue!");
            continue;
        }

        //________________________________________________________________________
        // Track DCAz
        fTreeCascVarBachDCAz = GetDCAz(bTrack);
        fTreeCascVarPosDCAz = GetDCAz(pTrack);
        fTreeCascVarNegDCAz = GetDCAz(nTrack);
        
        //________________________________________________________________________
        // Momentum info
        fTreeCascVarBachPx = lBachPx ;   fTreeCascVarBachPy = lBachPy ;    fTreeCascVarBachPz = lBachPz ;
        fTreeCascVarPosPx  = lPosPx  ;   fTreeCascVarPosPy  = lPosPy  ;    fTreeCascVarPosPz  = lPosPz ;
        fTreeCascVarNegPx  = lNegPx  ;   fTreeCascVarNegPy  = lNegPy  ;    fTreeCascVarNegPz  = lNegPz ;
        
        //________________________________________________________________________
        // Decay vtx info
        fTreeCascVarCascadeDecayX = lPosCasc[0] ; 
        fTreeCascVarCascadeDecayY = lPosCasc[1] ;  
        fTreeCascVarCascadeDecayZ = lPosCasc[2] ; 
        
        fTreeCascVarV0DecayX = lPosV0[0] ; 
        fTreeCascVarV0DecayY = lPosV0[1] ;  
        fTreeCascVarV0DecayZ = lPosV0[2] ; 
        
        fTreeCascVarPrimVertexX = lBestPrimaryVtxPos[0];
        fTreeCascVarPrimVertexY = lBestPrimaryVtxPos[1];
        fTreeCascVarPrimVertexZ = lBestPrimaryVtxPos[2];
        
        //________________________________________________________________________
        // Particle index
        fTreeCascVarBachIndex   = 0;
        fTreeCascVarPosIndex    = 0;
        fTreeCascVarNegIndex    = 0;
        
        //________________________________________________________________________
        //Check clusters
        fTreeCascVarPosITSClusters0 = pTrack->HasPointOnITSLayer(0);
        fTreeCascVarPosITSClusters1 = pTrack->HasPointOnITSLayer(1);
        fTreeCascVarPosITSClusters2 = pTrack->HasPointOnITSLayer(2);
        fTreeCascVarPosITSClusters3 = pTrack->HasPointOnITSLayer(3);
        fTreeCascVarPosITSClusters4 = pTrack->HasPointOnITSLayer(4);
        fTreeCascVarPosITSClusters5 = pTrack->HasPointOnITSLayer(5);
        
        fTreeCascVarNegITSClusters0 = nTrack->HasPointOnITSLayer(0);
        fTreeCascVarNegITSClusters1 = nTrack->HasPointOnITSLayer(1);
        fTreeCascVarNegITSClusters2 = nTrack->HasPointOnITSLayer(2);
        fTreeCascVarNegITSClusters3 = nTrack->HasPointOnITSLayer(3);
        fTreeCascVarNegITSClusters4 = nTrack->HasPointOnITSLayer(4);
        fTreeCascVarNegITSClusters5 = nTrack->HasPointOnITSLayer(5);
        
        fTreeCascVarBachITSClusters0 = bTrack->HasPointOnITSLayer(0);
        fTreeCascVarBachITSClusters1 = bTrack->HasPointOnITSLayer(1);
        fTreeCascVarBachITSClusters2 = bTrack->HasPointOnITSLayer(2);
        fTreeCascVarBachITSClusters3 = bTrack->HasPointOnITSLayer(3);
        fTreeCascVarBachITSClusters4 = bTrack->HasPointOnITSLayer(4);
        fTreeCascVarBachITSClusters5 = bTrack->HasPointOnITSLayer(5);
        
        //________________________________________________________________________
        //Check its clusters, shared
        fTreeCascVarPosITSSharedClusters0 = pTrack->HasSharedPointOnITSLayer(0);
        fTreeCascVarPosITSSharedClusters1 = pTrack->HasSharedPointOnITSLayer(1);
        fTreeCascVarPosITSSharedClusters2 = pTrack->HasSharedPointOnITSLayer(2);
        fTreeCascVarPosITSSharedClusters3 = pTrack->HasSharedPointOnITSLayer(3);
        fTreeCascVarPosITSSharedClusters4 = pTrack->HasSharedPointOnITSLayer(4);
        fTreeCascVarPosITSSharedClusters5 = pTrack->HasSharedPointOnITSLayer(5);
        
        fTreeCascVarNegITSSharedClusters0 = nTrack->HasSharedPointOnITSLayer(0);
        fTreeCascVarNegITSSharedClusters1 = nTrack->HasSharedPointOnITSLayer(1);
        fTreeCascVarNegITSSharedClusters2 = nTrack->HasSharedPointOnITSLayer(2);
        fTreeCascVarNegITSSharedClusters3 = nTrack->HasSharedPointOnITSLayer(3);
        fTreeCascVarNegITSSharedClusters4 = nTrack->HasSharedPointOnITSLayer(4);
        fTreeCascVarNegITSSharedClusters5 = nTrack->HasSharedPointOnITSLayer(5);
        
        fTreeCascVarBachITSSharedClusters0 = bTrack->HasSharedPointOnITSLayer(0);
        fTreeCascVarBachITSSharedClusters1 = bTrack->HasSharedPointOnITSLayer(1);
        fTreeCascVarBachITSSharedClusters2 = bTrack->HasSharedPointOnITSLayer(2);
        fTreeCascVarBachITSSharedClusters3 = bTrack->HasSharedPointOnITSLayer(3);
        fTreeCascVarBachITSSharedClusters4 = bTrack->HasSharedPointOnITSLayer(4);
        fTreeCascVarBachITSSharedClusters5 = bTrack->HasSharedPointOnITSLayer(5);
        
        //________________________________________________________________________
        //GetKinkIndex condition
        fTreeCascVarBachIsKink  = kFALSE;
        fTreeCascVarPosIsKink   = kFALSE;
        fTreeCascVarNegIsKink   = kFALSE;
        if( bTrack->GetKinkIndex(0)>0 ) fTreeCascVarBachIsKink  = kTRUE;
        if( pTrack->GetKinkIndex(0)>0 ) fTreeCascVarPosIsKink   = kTRUE;
        if( nTrack->GetKinkIndex(0)>0 ) fTreeCascVarNegIsKink   = kTRUE;
        
        //________________________________________________________________________
        //TOF info
        fTreeCascVarBachTOFExpTDiff     = bTrack->GetTOFExpTDiff( lMagField );
        fTreeCascVarPosTOFExpTDiff      = pTrack->GetTOFExpTDiff( lMagField );
        fTreeCascVarNegTOFExpTDiff      = nTrack->GetTOFExpTDiff( lMagField );
        
        fTreeCascVarBachTOFSignal       = bTrack->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarPosTOFSignal        = pTrack->GetTOFsignal() * 1.e-3; // in ns
        fTreeCascVarNegTOFSignal        = nTrack->GetTOFsignal() * 1.e-3; // in ns
        
        fTreeCascVarBachTOFBCid         = bTrack->GetTOFBunchCrossing( lMagField );
        fTreeCascVarPosTOFBCid          = pTrack->GetTOFBunchCrossing( lMagField );
        fTreeCascVarNegTOFBCid          = nTrack->GetTOFBunchCrossing( lMagField );
        
        //Copy VZERO information for this event
        fTreeCascVarAmplitudeV0A = lAmplitudeV0A;
        fTreeCascVarAmplitudeV0C = lAmplitudeV0C;
        //Copy IR information for this event
        fTreeCascVarClosestNonEmptyBC = lClosestNonEmptyBC;
        
        
        //________________________________________________________________________
        //Cowboy/sailor info regarding V0 inside cascade
        //Calculate vec prod with momenta projected to xy plane
        //Provisions for cowboy/sailor check
        Double_t lModp1 = TMath::Sqrt( lPosPx*lPosPx + lPosPy*lPosPy );
        Double_t lModp2 = TMath::Sqrt( lNegPx*lNegPx + lNegPy*lNegPy );
        
        Double_t lVecProd = (lPosPx*lNegPy - lPosPy*lNegPx) / (lModp1*lModp2);
        
        if ( lMagField < 0 ) lVecProd *= -1; //invert sign
        
        fTreeCascVarCowboyness = lVecProd;
        
        fTreeCascVarIsCowboy = kFALSE;
        if (lVecProd < 0) fTreeCascVarIsCowboy = kTRUE;
        
        Double_t lBachMod = TMath::Sqrt(lBachPx*lBachPx + lBachPy*lBachPy);
        Double_t lVecProdCasc = (lV0Px*lBachPy - lV0Py*lBachPx) / (TMath::Sqrt(lV0Px*lV0Px + lV0Py*lV0Py)*lBachMod);
        
        if ( lMagField < 0 ) lVecProdCasc *= -1; //invert sign
        
        fTreeCascVarCascadeCowboyness = lVecProdCasc;
        
        fTreeCascVarIsCascadeCowboy = kFALSE;
        if (lVecProdCasc < 0) fTreeCascVarIsCascadeCowboy = kTRUE;
        
        //________________________________________________________________________
        fTreeCascVarMagField = lMagField;
        fTreeCascVarRunNumber = lAODevent->GetRunNumber();
        fTreeCascVarEventNumber = ( ( ((ULong64_t)lAODevent->GetPeriodNumber() ) << 36 ) |
                                    ( ((ULong64_t)lAODevent->GetOrbitNumber () ) << 12 ) |
                                      ((ULong64_t)lAODevent->GetBunchCrossNumber() )  );
        
        //------------------------------------------------
        // Fill Tree!
        //------------------------------------------------
        //Apply rough selections
        if( fTreeCascVarPt < fMinPtToSave ) continue;
        if( fTreeCascVarPt > fMaxPtToSave ) continue;
        
        if(  fRejectCascKink  && 
            (
                fTreeCascVarBachIsKink || 
                fTreeCascVarPosIsKink  || 
                fTreeCascVarNegIsKink  
            ) 
          ) continue;
        
        
        if( !fCascSaveAddConfig                                                         &&
            
            TMath::Abs(fTreeCascVarBachEta)       < fMaxAbsEta                          && 
            TMath::Abs(fTreeCascVarPosEta )       < fMaxAbsEta                          && 
            TMath::Abs(fTreeCascVarNegEta )       < fMaxAbsEta                          &&
            
            (TMath::Abs(fTreeCascVarRapXi)        < fMaxAbsRap                          || 
             TMath::Abs(fTreeCascVarRapOmega)     < fMaxAbsRap                         )&&
            
            fTreeCascVarLeastNbrCrossedRows       > fMinNbrCrossedRows                  &&
            
            //Do Selection : reject if the both mass selections are simultenously false
            (   //START XI SELECTIONS
                (1.32-0.075 < fTreeCascVarMassAsXi      && fTreeCascVarMassAsXi     < 1.32+0.075)    
                ||
                //START OMEGA SELECTIONS
                (1.68-0.075 < fTreeCascVarMassAsOmega   && fTreeCascVarMassAsOmega  < 1.68+0.075)
            )
        )
        {
            IsEvtWithCascade = kTRUE;
            fTreeCascade->Fill();
        }
            
        
        if( fCascSaveAddConfig )
        {
            Bool_t lCascadeIsSaved = kFALSE;
            
            Double_t lLambdaPDGMass = 1.115683;
            Double_t lXiPDGMass     = 1.32171;
            Double_t lOmegaPDGMass  = 1.67245;
            
            Int_t       lCascResultCharge;
            Bool_t      lCascResultIsOmega;
            
            Double_t    lInvMass;
            Double_t    lPDGMass;
            Double_t    lV0Mass;
            Double_t    lRap;
            Double_t    lBachTPCNSigma;
            Double_t    lPosTPCNSigma;
            Double_t    lNegTPCNSigma;
            
            Int_t nCascRes = fCascadeResult.GetEntries();
            
            for (Int_t iCascRes = 0 ; iCascRes < nCascRes ; iCascRes++) 
            {
                // once the candidate has passed the selection of 
                // one of the cascade type (Xi+, Xi-, Omega+, Omega-), 
                // we can skip the other cut sets
                if( lCascadeIsSaved ) break;
                
                AliCascadeResult *lCascadeResult = (AliCascadeResult *)fCascadeResult[iCascRes];
                
                //========================================================================
                if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiMinus     )
                {
                    lCascResultCharge       = -1;
                    lCascResultIsOmega      = kFALSE;
                    
                    lPDGMass                = lXiPDGMass;
                    lInvMass                = fTreeCascVarMassAsXi;
                    lV0Mass                 = fTreeCascVarV0MassAsLambda;
                    lRap                    = fTreeCascVarRapXi;
                    
                    lBachTPCNSigma          = fTreeCascVarBachNSigmaPion;
                    lPosTPCNSigma           = fTreeCascVarPosNSigmaProton;
                    lNegTPCNSigma           = fTreeCascVarNegNSigmaPion;
                }
                if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kXiPlus      )
                {
                    lCascResultCharge       = +1;
                    lCascResultIsOmega      = kFALSE;
                    
                    lPDGMass                = lXiPDGMass;
                    lInvMass                = fTreeCascVarMassAsXi;
                    lV0Mass                 = fTreeCascVarV0MassAsAntiLambda;
                    lRap                    = fTreeCascVarRapXi;
                    
                    lBachTPCNSigma          = fTreeCascVarBachNSigmaPion;
                    lPosTPCNSigma           = fTreeCascVarPosNSigmaPion;
                    lNegTPCNSigma           = fTreeCascVarNegNSigmaProton;
                }
                if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaMinus     )
                {
                    lCascResultCharge       = -1;
                    lCascResultIsOmega      = kTRUE;
                    
                    lPDGMass                = lOmegaPDGMass;
                    lInvMass                = fTreeCascVarMassAsOmega;
                    lV0Mass                 = fTreeCascVarV0MassAsLambda;
                    lRap                    = fTreeCascVarRapOmega;
                    
                    lBachTPCNSigma          = fTreeCascVarBachNSigmaKaon;
                    lPosTPCNSigma           = fTreeCascVarPosNSigmaProton;
                    lNegTPCNSigma           = fTreeCascVarNegNSigmaPion;
                }
                if ( lCascadeResult->GetMassHypothesis() == AliCascadeResult::kOmegaPlus      )
                {
                    lCascResultCharge       = +1;
                    lCascResultIsOmega      = kTRUE;
                    
                    lPDGMass                = lOmegaPDGMass;
                    lInvMass                = fTreeCascVarMassAsOmega;
                    lV0Mass                 = fTreeCascVarV0MassAsAntiLambda;
                    lRap                    = fTreeCascVarRapOmega;
                    
                    lBachTPCNSigma          = fTreeCascVarBachNSigmaKaon;
                    lPosTPCNSigma           = fTreeCascVarPosNSigmaPion;
                    lNegTPCNSigma           = fTreeCascVarNegNSigmaProton;
                }
                
                if (
                    //Check 1 : Charge
                    fTreeCascVarCharge == lCascResultCharge                                                     &&
                    
                    //Check 2: Basic Acceptance cuts
                    // Bachelor
                    fTreeCascVarBachEta                     > lCascadeResult->GetCutMinEtaTracks()              && 
                    fTreeCascVarBachEta                     < lCascadeResult->GetCutMaxEtaTracks()              &&
                    // Positive
                    fTreeCascVarPosEta                      > lCascadeResult->GetCutMinEtaTracks()              && 
                    fTreeCascVarPosEta                      < lCascadeResult->GetCutMaxEtaTracks()              &&
                    // Negative
                    fTreeCascVarNegEta                      > lCascadeResult->GetCutMinEtaTracks()              &&
                    fTreeCascVarNegEta                      < lCascadeResult->GetCutMaxEtaTracks()              &&
                    
                    // Corresponding rapidity
                    lRap                                    > lCascadeResult->GetCutMinRapidity()               &&
                    lRap                                    < lCascadeResult->GetCutMaxRapidity()               &&
                    
                    //Check 3: Topological Variables
                    // - V0 Selections
                    fTreeCascVarDCAPosToPV                  > lCascadeResult->GetCutDCAPosToPV()                &&
                    fTreeCascVarDCANegToPV                  > lCascadeResult->GetCutDCANegToPV()                &&
                    fTreeCascVarDCAV0Dau                    < lCascadeResult->GetCutDCAV0Daughters()            &&
                    fTreeCascVarDCAV0ToPV                   > lCascadeResult->GetCutDCAV0ToPV()                 &&
                    fTreeCascVarV0Radius                    > lCascadeResult->GetCutV0Radius()                  &&
                    fTreeCascVarV0CosPA                     > lCascadeResult->GetCutV0CosPA()                   &&
                    // - Cascade Selections
                    TMath::Abs(lV0Mass-lLambdaPDGMass)      < lCascadeResult->GetCutV0Mass()                    &&
                    fTreeCascVarDCABachToPV                 > lCascadeResult->GetCutDCABachToPV()               &&
                    fTreeCascVarDCACascDau                  < lCascadeResult->GetCutDCACascDaughters()          &&
                    fTreeCascVarCascRadius                  > lCascadeResult->GetCutCascRadius()                &&
                    fTreeCascVarCascCosPA                   > lCascadeResult->GetCutCascCosPA()                 &&
                    fTreeCascVarLeastNbrCrossedRows         > lCascadeResult->GetCutLeastNumberOfCrossedRows()  &&
                    lPDGMass*fTreeCascVarDistOverTotMom     < lCascadeResult->GetCutProperLifetime()            &&
                    // - Miscellaneous
                    fTreeCascVarLeastNbrClusters            > lCascadeResult->GetCutLeastNumberOfClusters()     &&
                    
                    //Check 4: TPC dEdx selections
                    TMath::Abs(lBachTPCNSigma)              < lCascadeResult->GetCutTPCdEdx()                   &&
                    TMath::Abs(lPosTPCNSigma )              < lCascadeResult->GetCutTPCdEdx()                   &&
                    TMath::Abs(lNegTPCNSigma )              < lCascadeResult->GetCutTPCdEdx()                   &&
                    
                    //Check 5: Min/Max V0 Lifetime cut
                    lLambdaPDGMass*fTreeCascVarV0DistOverTotMom > lCascadeResult->GetCutMinV0Lifetime()         &&
                    lLambdaPDGMass*fTreeCascVarV0DistOverTotMom < lCascadeResult->GetCutMaxV0Lifetime()         &&
                    
                    //Check 6: Min Track Length 
                    fTreeCascVarMinTrackLength              > lCascadeResult->GetCutMinTrackLength()            &&
                    
                    //Check 7: modern track quality selections
                    fTreeCascVarNbrCrossedRowsOverLength    > lCascadeResult->GetCutMinCrossedRowsOverLength()  &&

                    //Check 8: Max Chi2/Clusters 
                    fTreeCascVarMaxChi2PerCluster           < lCascadeResult->GetCutMaxChi2PerCluster()         &&
                    
                    //Check 9: Experimental DCA Bachelor to Baryon cut
                    fTreeCascVarDCABachToBaryon             > lCascadeResult->GetCutDCABachToBaryon()           &&
                    
                    //Check 10: Cut on the Bach Baryon pointing angle (not CosPA)
                    TMath::ACos(fTreeCascVarWrongCosPA)  > TMath::ACos(lCascadeResult->GetCutBachBaryonCosPA()) &&
                    
                    //Check 11a: Xi rejection for Omega analysis
                    (   //Is not Omega, so skip this cut
                        !lCascResultIsOmega   ||
                        //Is Omega --> cut on the Xi mass
                        (lCascResultIsOmega   && 
                        TMath::Abs( fTreeCascVarMassAsXi - lXiPDGMass ) > lCascadeResult->GetCutXiRejection() ) 
                        
                    )                                                                                           &&
                    
                    //Check 11b: Mass cut rejection
                    TMath::Abs( lInvMass - lPDGMass )       < 0.075                                             &&       
                    
                    //Check 12: kITSrefit track selection if requested
                    (   //do not need ITS refit tracks --> skip this cut
                        !( lCascadeResult->GetCutUseITSRefitTracks() )                  ||
                        //otherwise : need to have ITS refit tracks
                        (   ( fTreeCascVarBachTrackStatus   & AliAODTrack::kITSrefit )  &&
                            ( fTreeCascVarPosTrackStatus    & AliAODTrack::kITSrefit )  &&
                            ( fTreeCascVarNegTrackStatus    & AliAODTrack::kITSrefit )  
                        )
                    )                                                                                           &&
                    
                    //Check 13: cowboy/sailor for V0
                    ( lCascadeResult->GetCutIsCowboy()  ==  0                                            ||
                    (lCascadeResult->GetCutIsCowboy()   ==  1 && fTreeCascVarIsCowboy           )        ||
                    (lCascadeResult->GetCutIsCowboy()   == -1 && !fTreeCascVarIsCowboy          )      )        &&//end cowboy/sailor
                    
                    //Check 14: cowboy/sailor for cascade
                    ( lCascadeResult->GetCutIsCascadeCowboy()   ==  0                                      ||
                    (lCascadeResult->GetCutIsCascadeCowboy()    ==  1 && fTreeCascVarIsCascadeCowboy   )   ||
                    (lCascadeResult->GetCutIsCascadeCowboy()    == -1 && !fTreeCascVarIsCascadeCowboy  ) )      &&//end cowboy/sailor
                    
                    //Check 15: ITS or TOF required 
                    (
                    !( lCascadeResult->GetCutITSorTOF() ) || 
                        (   //One of the daughter track has ITSrefit 
                            (   (fTreeCascVarBachTrackStatus  & AliAODTrack::kITSrefit)  ||
                                (fTreeCascVarPosTrackStatus   & AliAODTrack::kITSrefit)  ||
                                (fTreeCascVarNegTrackStatus   & AliAODTrack::kITSrefit) 
                            )||
                            //Or one of the daughter track has a TOF signal
                            (   (TMath::Abs(fTreeCascVarBachTOFExpTDiff+2500. )     > 1e-6 ) ||
                                (TMath::Abs(fTreeCascVarPosTOFExpTDiff+2500.  )     > 1e-6 ) ||
                                (TMath::Abs(fTreeCascVarNegTOFExpTDiff+2500.  )     > 1e-6 ) 
                            )
                        )
                        
                    )
                    
                )//end major if
                {
                    lCascadeIsSaved = kTRUE;
                    IsEvtWithCascade = kTRUE;
                    fTreeCascade->Fill();
                }
                    
            }
                
        }  
    
    }//end loop over the cascades
    
    AliWarning(Form("Number of cascades saved: %lli", fTreeCascade->GetEntries()));
    
    //------------------------------------------------
    // MAIN V0 LOOP STARTS HERE
    //------------------------------------------------
    
    if(IsEvtWithCascade && fSaveV0s){//Adaptation of AliRsnMiniAnalysisTask
        
        Long_t nV0s = 0;
        nV0s = lAODevent->GetNumberOfV0s();
        
        for (Int_t iV0 = 0; iV0 < nV0s; iV0++){
            
            //------------------------------------------------
            // Initializations
            //------------------------------------------------
            Double_t lPosV0[3] = {0.,0.,0.}; // Position of V0 
        
            Bool_t   lOnFlyStatus           = kFALSE; // if kTRUE, then this V0 is recontructed "on fly" during the tracking
            Double_t lV0MassAsK0s           = 0.;
            Double_t lV0MassAsLambda        = 0.;
            Double_t lV0MassAsAntiLambda    = 0.;
            
            AliAODv0 *v0 = lAODevent->GetV0(iV0);
            if (!v0) continue;
            
            //CheckChargeV0( v0 ); //FIXME this won't work for AOD as is
            v0->GetXYZ(lPosV0);
            
            //Gather tracks informations
            AliAODTrack *pTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) ); //0->Positive Daughter
            AliAODTrack *nTrack = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) ); //1->Negative Daughter
            
            if (!pTrack || !nTrack) {
                AliWarning("ERROR: Could not retrieve one of the daughter tracks");
                continue;
            }
            
            //Findable clusters > 0 condition
            if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;
            
            //Test different mass hypothesis : K0s, Lambda, AntiLambda
            lV0MassAsK0s            = v0->MassK0Short();
            lV0MassAsLambda         = v0->MassLambda();
            lV0MassAsAntiLambda     = v0->MassAntiLambda();
            
            Double_t lPosPx = 0.; Double_t lPosPy = 0.; Double_t lPosPz = 0.;
            Double_t lNegPx = 0.; Double_t lNegPy = 0.; Double_t lNegPz = 0.;
            lPosPx  = pTrack->Px() ; lPosPy  = pTrack->Py() ; lPosPz  = pTrack->Pz() ;
            lNegPx  = nTrack->Px() ; lNegPy  = nTrack->Py() ; lNegPz  = nTrack->Pz() ;
            
            //------------------------------------------------
            // Calculation of the variables related to V0
            //------------------------------------------------
            lOnFlyStatus                    = v0->GetOnFlyStatus();
            fTreeV0VarChi2V0                = v0->Chi2V0();
            fTreeV0VarMassAsK0s             = lV0MassAsK0s;
            fTreeV0VarMassAsLambda          = lV0MassAsLambda;
            fTreeV0VarMassAsAntiLambda      = lV0MassAsAntiLambda;
        
            fTreeV0VarPtot                  = TMath::Sqrt( v0->Ptot2V0() );
            fTreeV0VarPt                    = TMath::Sqrt( v0->Pt2V0() );

            fTreeV0VarRapK0Short            = v0->RapK0Short();
            fTreeV0VarRapLambda             = v0->RapLambda();
            fTreeV0VarPosEta                = pTrack->Eta();
            fTreeV0VarNegEta                = nTrack->Eta();
            
            fTreeV0VarPhi                   = MyPhi(lPosPx+lNegPx, lPosPy+lNegPy);
            fTreeV0VarTheta                 = MyTheta(fTreeV0VarPt, lPosPz+lNegPz);
            
            fTreeV0VarAlpha         = v0->AlphaV0();
            fTreeV0VarPtArm         = v0->PtArmV0();
            fTreeV0VarDCAV0Dau      = v0->DcaV0Daughters();
            fTreeV0VarDCAV0ToPV     = v0->DcaV0ToPrimVertex();
            fTreeV0VarDCAPosToPV    = v0->DcaPosToPrimVertex();
            fTreeV0VarDCANegToPV    = v0->DcaNegToPrimVertex();
        
            fTreeV0VarRadius    = TMath::Sqrt( lPosV0[0]*lPosV0[0]  +  lPosV0[1]*lPosV0[1] );
            fTreeV0VarCosPA     = v0->CosPointingAngle(lBestPrimaryVtxPos);
            
            // Filter like-sign V0 
            if ( pTrack->GetSign() == nTrack->GetSign() ) continue;

            
            //________________________________________________________________________
            // Track quality cuts
            //Least Nbr of Crossed Rows
            Int_t lPosNbrCrossedRows  = pTrack->GetTPCClusterInfo(2,1);
            Int_t lNegNbrCrossedRows  = nTrack->GetTPCClusterInfo(2,1);
            
            Int_t lLeastNbrCrossedRows = (Int_t) lPosNbrCrossedRows;
            if( lNegNbrCrossedRows < lLeastNbrCrossedRows ) lLeastNbrCrossedRows = (Int_t) lNegNbrCrossedRows;
            
            fTreeV0VarLeastNbrCrossedRows = lLeastNbrCrossedRows;
            
            //Compute ratio Crossed Rows / Findable clusters
            //Note: above test avoids division by zero!
            Double_t lPosTrackCrossedRowsOverFindable = lPosNbrCrossedRows / ((Double_t)(pTrack->GetTPCNclsF()));
            Double_t lNegTrackCrossedRowsOverFindable = lNegNbrCrossedRows / ((Double_t)(nTrack->GetTPCNclsF()));
            
            fTreeV0VarLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
            if( lNegTrackCrossedRowsOverFindable < fTreeV0VarLeastRatioCrossedRowsOverFindable ) 
                fTreeV0VarLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
            
            // Least Nbr of clusters        
            Int_t lPosTPCClusters        = 0.;
            Int_t lNegTPCClusters        = 0.;

            lPosTPCClusters   = pTrack->GetTPCNcls();
            lNegTPCClusters   = nTrack->GetTPCNcls();
            
            Int_t lLeastNbrOfClusters    = lPosTPCClusters;
            if( lNegTPCClusters < lLeastNbrOfClusters ) lLeastNbrOfClusters = lNegTPCClusters;
            
            fTreeV0VarLeastNbrClusters = lLeastNbrOfClusters;
            
            //Min track length
            Double_t lPosTrackLength        = 0.;
            Double_t lNegTrackLength        = 0.;
            
            lPosTrackLength    = GetLengthInActiveZone(pTrack, 2.0, 220.0, lMagField);
            lNegTrackLength    = GetLengthInActiveZone(nTrack, 2.0, 220.0, lMagField);
            
            Double_t lSmallestTrackLength = lPosTrackLength;
            if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
            
            fTreeV0VarMinTrackLength = lSmallestTrackLength;
            
            //Max Chi2 per cluster
            Double_t lBiggestChi2PerCluster = 0.;
            Double_t lPosChi2PerCluster = pTrack->GetTPCchi2() / ((Float_t) lPosTPCClusters);
            Double_t lNegChi2PerCluster = nTrack->GetTPCchi2() / ((Float_t) lNegTPCClusters);
            
            if( lPosChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lPosChi2PerCluster;
            if( lNegChi2PerCluster > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lNegChi2PerCluster;
        
            fTreeV0VarMaxChi2PerCluster = lBiggestChi2PerCluster;
            
            //________________________________________________________________________
            // Decay length info
            fTreeV0VarDistOverTotMom  = TMath::Sqrt(   TMath::Power( lPosV0[0] - lBestPrimaryVtxPos[0] , 2) 
                                                    + TMath::Power( lPosV0[1] - lBestPrimaryVtxPos[1] , 2) 
                                                    + TMath::Power( lPosV0[2] - lBestPrimaryVtxPos[2] , 2) );
            fTreeV0VarDistOverTotMom /= fTreeV0VarPtot;
            
            
            //------------------------------------------------
            // TPC dEdx information
            //------------------------------------------------
            fTreeV0VarPosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
            fTreeV0VarPosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
            fTreeV0VarNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion   );
            fTreeV0VarNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
            
            //------------------------------------------------
            // ITS info 
            //------------------------------------------------
            fTreeV0VarPosITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( pTrack, AliPID::kPion );
            fTreeV0VarPosITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( pTrack, AliPID::kProton );
            fTreeV0VarNegITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( nTrack, AliPID::kPion   );
            fTreeV0VarNegITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( nTrack, AliPID::kProton );
            
            //------------------------------------------------
            // TOF info 
            //------------------------------------------------
            fTreeV0VarPosTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( pTrack, AliPID::kPion );
            fTreeV0VarPosTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( pTrack, AliPID::kProton );
            fTreeV0VarNegTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( nTrack, AliPID::kPion   );
            fTreeV0VarNegTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( nTrack, AliPID::kProton );
            
            //------------------------------------------------
            // Raw TPC dEdx + PIDForTracking information
            //------------------------------------------------
            
            //Acquire TPC Signals
            fTreeV0VarPosdEdx = pTrack->GetTPCsignal();
            fTreeV0VarNegdEdx = nTrack->GetTPCsignal();
            
            //Acquire PID For Tracking : uses AliPID EParticleType code (0=electron, 1=muon, 2=pion, etc)
            fTreeV0VarPosPIDForTracking = pTrack->GetPIDForTracking();
            fTreeV0VarNegPIDForTracking = nTrack->GetPIDForTracking();
            
            //________________________________________________________________________
            // Track status
            fTreeV0VarPosTrackStatus      = pTrack->GetStatus();
            fTreeV0VarNegTrackStatus      = nTrack->GetStatus();
            
            if ((fTreeV0VarPosTrackStatus & AliAODTrack::kTPCrefit)  == 0) {
                AliDebug(1, "Pb w/ V0 Pos. track has no TPCrefit ... continue!");
                continue;
            }
            if ((fTreeV0VarNegTrackStatus & AliAODTrack::kTPCrefit)  == 0) {
                AliDebug(1, "Pb w/ V0 Neg. track has no TPCrefit ... continue!");
                continue;
            }

            //________________________________________________________________________
            // Track DCAz
            fTreeV0VarPosDCAz = GetDCAz(pTrack);
            fTreeV0VarNegDCAz = GetDCAz(nTrack);
            
            //________________________________________________________________________
            // Momentum info
            fTreeV0VarPosPx  = lPosPx  ;   fTreeV0VarPosPy  = lPosPy  ;    fTreeV0VarPosPz  = lPosPz ;
            fTreeV0VarNegPx  = lNegPx  ;   fTreeV0VarNegPy  = lNegPy  ;    fTreeV0VarNegPz  = lNegPz ;
            
            //________________________________________________________________________
            // Decay vtx info
            fTreeV0VarV0DecayX = lPosV0[0] ; 
            fTreeV0VarV0DecayY = lPosV0[1] ;  
            fTreeV0VarV0DecayZ = lPosV0[2] ; 
            
            fTreeV0VarPrimVertexX = lBestPrimaryVtxPos[0];
            fTreeV0VarPrimVertexY = lBestPrimaryVtxPos[1];
            fTreeV0VarPrimVertexZ = lBestPrimaryVtxPos[2];
            
            //________________________________________________________________________
            //Check clusters
            fTreeV0VarPosITSClusters0 = pTrack->HasPointOnITSLayer(0);
            fTreeV0VarPosITSClusters1 = pTrack->HasPointOnITSLayer(1);
            fTreeV0VarPosITSClusters2 = pTrack->HasPointOnITSLayer(2);
            fTreeV0VarPosITSClusters3 = pTrack->HasPointOnITSLayer(3);
            fTreeV0VarPosITSClusters4 = pTrack->HasPointOnITSLayer(4);
            fTreeV0VarPosITSClusters5 = pTrack->HasPointOnITSLayer(5);
            
            fTreeV0VarNegITSClusters0 = nTrack->HasPointOnITSLayer(0);
            fTreeV0VarNegITSClusters1 = nTrack->HasPointOnITSLayer(1);
            fTreeV0VarNegITSClusters2 = nTrack->HasPointOnITSLayer(2);
            fTreeV0VarNegITSClusters3 = nTrack->HasPointOnITSLayer(3);
            fTreeV0VarNegITSClusters4 = nTrack->HasPointOnITSLayer(4);
            fTreeV0VarNegITSClusters5 = nTrack->HasPointOnITSLayer(5);
            
            //________________________________________________________________________
            //Check its clusters, shared
            fTreeV0VarPosITSSharedClusters0 = pTrack->HasSharedPointOnITSLayer(0);
            fTreeV0VarPosITSSharedClusters1 = pTrack->HasSharedPointOnITSLayer(1);
            fTreeV0VarPosITSSharedClusters2 = pTrack->HasSharedPointOnITSLayer(2);
            fTreeV0VarPosITSSharedClusters3 = pTrack->HasSharedPointOnITSLayer(3);
            fTreeV0VarPosITSSharedClusters4 = pTrack->HasSharedPointOnITSLayer(4);
            fTreeV0VarPosITSSharedClusters5 = pTrack->HasSharedPointOnITSLayer(5);
            
            fTreeV0VarNegITSSharedClusters0 = nTrack->HasSharedPointOnITSLayer(0);
            fTreeV0VarNegITSSharedClusters1 = nTrack->HasSharedPointOnITSLayer(1);
            fTreeV0VarNegITSSharedClusters2 = nTrack->HasSharedPointOnITSLayer(2);
            fTreeV0VarNegITSSharedClusters3 = nTrack->HasSharedPointOnITSLayer(3);
            fTreeV0VarNegITSSharedClusters4 = nTrack->HasSharedPointOnITSLayer(4);
            fTreeV0VarNegITSSharedClusters5 = nTrack->HasSharedPointOnITSLayer(5);
            
            //________________________________________________________________________
            //GetKinkIndex condition
            fTreeV0VarPosIsKink   = kFALSE;
            fTreeV0VarNegIsKink   = kFALSE;
            if( pTrack->GetKinkIndex(0)>0 ) fTreeV0VarPosIsKink   = kTRUE;
            if( nTrack->GetKinkIndex(0)>0 ) fTreeV0VarNegIsKink   = kTRUE;
            
            //________________________________________________________________________
            //TOF info
            fTreeV0VarPosTOFExpTDiff      = pTrack->GetTOFExpTDiff( lMagField );
            fTreeV0VarNegTOFExpTDiff      = nTrack->GetTOFExpTDiff( lMagField );
            
            fTreeV0VarPosTOFSignal        = pTrack->GetTOFsignal() * 1.e-3; // in ns
            fTreeV0VarNegTOFSignal        = nTrack->GetTOFsignal() * 1.e-3; // in ns
            
            fTreeV0VarPosTOFBCid          = pTrack->GetTOFBunchCrossing( lMagField );
            fTreeV0VarNegTOFBCid          = nTrack->GetTOFBunchCrossing( lMagField );
    
            //________________________________________________________________________
            //Cowboy/sailor info regarding V0 
            //Calculate vec prod with momenta projected to xy plane
            //Provisions for cowboy/sailor check
            Double_t lModp1 = TMath::Sqrt( lPosPx*lPosPx + lPosPy*lPosPy );
            Double_t lModp2 = TMath::Sqrt( lNegPx*lNegPx + lNegPy*lNegPy );
            
            //Calculate vec prod with momenta projected to xy plane
            Double_t lVecProd = (lPosPx*lNegPy - lPosPy*lNegPx) / (lModp1*lModp2);
            
            if ( lMagField < 0 ) lVecProd *= -1; //invert sign
            
            fTreeV0VarIsCowboy = kFALSE;
            if (lVecProd < 0) fTreeV0VarIsCowboy = kTRUE;
            
            //________________________________________________________________________
            fTreeV0VarMagField = lMagField;
            fTreeV0VarRunNumber = lAODevent->GetRunNumber();
            fTreeV0VarEventNumber = ( ( ((ULong64_t)lAODevent->GetPeriodNumber() ) << 36 ) |
                                        ( ((ULong64_t)lAODevent->GetOrbitNumber () ) << 12 ) |
                                        ((ULong64_t)lAODevent->GetBunchCrossNumber() )  );
            
            //------------------------------------------------
            // Fill Tree!
            //------------------------------------------------
            //Apply rough selections
            //pT window
            if( fTreeV0VarPt < fMinPtToSave ) continue;
            if( fTreeV0VarPt > fMaxPtToSave ) continue;
            
            // Reject Kink topology
            if(  fRejectV0Kink  && ( fTreeV0VarPosIsKink || fTreeV0VarNegIsKink ) ) continue;
            
            //Second Selection: rough 20-sigma band, parametric.
            //K0Short: Enough to parametrize peak broadening with linear function.
            Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeV0VarPt;
            Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeV0VarPt;
            //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
            //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
            Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeV0VarPt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeV0VarPt);
            Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeV0VarPt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeV0VarPt);
            
            if( !fV0SaveAddConfig )
            {
                
                if( //Reject on-the-fly V0s
                    !lOnFlyStatus                                                       &&
                    
                    TMath::Abs(fTreeV0VarPosEta)            < fMaxAbsEta                && 
                    TMath::Abs(fTreeV0VarNegEta)            < fMaxAbsEta                &&
                    
                    (   
                        TMath::Abs(fTreeV0VarRapK0Short)    < fMaxAbsRap  || 
                        TMath::Abs(fTreeV0VarRapLambda)     < fMaxAbsRap    
                    )                                                                   &&
                    
                    fTreeV0VarLeastNbrCrossedRows           > fMinNbrCrossedRows        &&
                    //Do Selection : reject if the three mass selections are false
                    (   //Case 1: Lambda Selection
                        (lLowerLimitLambda   < fTreeV0VarMassAsLambda      && fTreeV0VarMassAsLambda       < lUpperLimitLambda ) ||
                        //Case 2: AntiLambda Selection
                        (lLowerLimitLambda   < fTreeV0VarMassAsAntiLambda  && fTreeV0VarMassAsAntiLambda   < lUpperLimitLambda ) ||
                        //Case 3: K0Short Selection
                        (lLowerLimitK0Short  < fTreeV0VarMassAsK0s         && fTreeV0VarMassAsK0s          < lUpperLimitK0Short)
                    )
                )
                {
                    fTreeV0->Fill();
                }
                
            }
            
            
            if( fV0SaveAddConfig )
            {
                Bool_t lV0IsSaved = kFALSE; 
                
                Double_t lK0sPDGMass    = 0.493677;
                Double_t lLambdaPDGMass = 1.115683;
                
                Bool_t      lIsK0s;
                Double_t    lInvMass;
                Double_t    lPDGMass;
                Double_t    lV0Mass;
                Double_t    lRap;
                Double_t    lPosTPCNSigma;
                Double_t    lNegTPCNSigma;
                
                Int_t nV0Res = fV0Result.GetEntries();
                
                for (Int_t iV0Res = 0 ; iV0Res < nV0Res ; iV0Res++) 
                {
                    // once the candidate has passed the selection of 
                    // one of the V0 type (K0s, Lambda, AntiLambda), 
                    // we can skip the other cut sets
                    if( lV0IsSaved ) break;
                
                    AliV0Result *lV0Result = (AliV0Result *)fV0Result[iV0Res];
                    
                    //========================================================================
                    if ( lV0Result->GetMassHypothesis() == AliV0Result::kK0Short        )
                    {
                        lIsK0s                  = kTRUE;
                        lInvMass                = fTreeV0VarMassAsK0s;
                        lPDGMass                = lK0sPDGMass;
                        lRap                    = fTreeV0VarRapK0Short;
                        
                        lPosTPCNSigma           = fTreeV0VarPosNSigmaPion;
                        lNegTPCNSigma           = fTreeV0VarNegNSigmaPion;
                    }
                    if ( lV0Result->GetMassHypothesis() == AliV0Result::kLambda         )
                    {
                        lIsK0s                  = kFALSE;
                        lInvMass                = fTreeV0VarMassAsLambda;
                        lPDGMass                = lLambdaPDGMass;
                        lRap                    = fTreeV0VarRapLambda;
                        
                        lPosTPCNSigma           = fTreeV0VarPosNSigmaProton;
                        lNegTPCNSigma           = fTreeV0VarNegNSigmaPion;
                    }
                    if ( lV0Result->GetMassHypothesis() == AliV0Result::kAntiLambda     )
                    {
                        lIsK0s                  = kFALSE;
                        lInvMass                = fTreeV0VarMassAsAntiLambda;
                        lPDGMass                = lLambdaPDGMass;
                        lRap                    = fTreeV0VarRapLambda;
                        
                        lPosTPCNSigma           = fTreeV0VarPosNSigmaPion;
                        lNegTPCNSigma           = fTreeV0VarNegNSigmaProton;
                    }
                    
                    if (
                        //Check 1: Offline Vertexer
                        lOnFlyStatus == lV0Result->GetUseOnTheFly() &&
                        
                        //Check 2: Basic Acceptance cuts
                        // Positive
                        fTreeV0VarPosEta                                > lV0Result->GetCutMinEtaTracks()                           && 
                        fTreeV0VarPosEta                                < lV0Result->GetCutMaxEtaTracks()                           &&
                        // Negative
                        fTreeV0VarNegEta                                > lV0Result->GetCutMinEtaTracks()                           &&
                        fTreeV0VarNegEta                                < lV0Result->GetCutMaxEtaTracks()                           &&
                        
                        // Corresponding rapidity
                        lRap                                            > lV0Result->GetCutMinRapidity()                            &&
                        lRap                                            < lV0Result->GetCutMaxRapidity()                            &&
                        
                        //Check 3: Topological Variables
                        fTreeV0VarDCAPosToPV                            > lV0Result->GetCutDCAPosToPV()                             &&
                        fTreeV0VarDCANegToPV                            > lV0Result->GetCutDCANegToPV()                             &&
                        fTreeV0VarDCAV0Dau                              < lV0Result->GetCutDCAV0Daughters()                         &&
                        fTreeV0VarRadius                                > lV0Result->GetCutV0Radius()                               &&
                        fTreeV0VarCosPA                                 > lV0Result->GetCutV0CosPA()                                &&
                        fTreeV0VarLeastNbrCrossedRows                   > lV0Result->GetCutLeastNumberOfCrossedRows()               &&
                        lPDGMass*fTreeV0VarDistOverTotMom               < lV0Result->GetCutProperLifetime()                         &&
                        // - Miscellaneous
                        fTreeV0VarLeastRatioCrossedRowsOverFindable     > lV0Result->GetCutLeastNumberOfCrossedRowsOverFindable()   &&
                        
                        //Check 4: TPC dEdx selections
                        TMath::Abs(lPosTPCNSigma )                      < lV0Result->GetCutTPCdEdx()                                &&
                        TMath::Abs(lNegTPCNSigma )                      < lV0Result->GetCutTPCdEdx()                                &&
                        
                        //Check 5: Min Track Length 
                        fTreeV0VarMinTrackLength                        > lV0Result->GetCutMinTrackLength()                         &&

                        //Check 6: Max Chi2/Clusters 
                        fTreeV0VarMaxChi2PerCluster                     < lV0Result->GetCutMaxChi2PerCluster()                      &&
                        
                        //Check 7a: K0s/Lambda mass rejection
                        (   //Is K0s --> cut on the Lambda mass
                            (
                                lIsK0s && 
                                TMath::Abs( fTreeV0VarMassAsLambda - lLambdaPDGMass ) > lV0Result->GetCutCompetingV0Rejection() 
                            ) 
                            ||
                            //Is Lambda --> cut on the K0s mass
                            (
                                !lIsK0s && 
                                TMath::Abs( fTreeV0VarMassAsK0s - lK0sPDGMass       ) > lV0Result->GetCutCompetingV0Rejection() 
                            )
                        )                                                                                                          &&
                        
                        //Check 7b: K0s/Lambda mass cut
                        (   //Is K0s
                            (
                                lIsK0s && 
                                lLowerLimitK0Short  < lInvMass  &&  lInvMass    < lUpperLimitK0Short
                            ) 
                            ||
                            //Is Lambda
                            (
                                !lIsK0s && 
                                lLowerLimitLambda   < lInvMass  &&  lInvMass    < lUpperLimitLambda
                            )
                        )                                                                                                          &&
                        
                        //Check 8: kITSrefit track selection if requested
                        (   //do not need ITS refit tracks --> skip this cut
                            !( lV0Result->GetCutUseITSRefitTracks() )                           ||
                            //otherwise : need to have ITS refit tracks
                            (   ( fTreeV0VarPosTrackStatus    & AliAODTrack::kITSrefit )        &&
                                ( fTreeV0VarNegTrackStatus    & AliAODTrack::kITSrefit )  
                            )
                        )                                                                                                           &&
                        
                        //Check 9: cowboy/sailor for V0
                        ( lV0Result->GetCutIsCowboy()  ==  0                                    ||
                        (lV0Result->GetCutIsCowboy()   ==  1 && fTreeV0VarIsCowboy           )  ||
                        (lV0Result->GetCutIsCowboy()   == -1 && !fTreeV0VarIsCowboy          ) )                                    &&
                        
                        //Check 10: ITS or TOF required 
                        (
                          !( lV0Result->GetCutITSorTOF() )                                      || 
                            (   //One of the daughter track has ITSrefit 
                                (   (fTreeV0VarPosTrackStatus   & AliAODTrack::kITSrefit)       ||
                                    (fTreeV0VarNegTrackStatus   & AliAODTrack::kITSrefit) 
                                )||
                                //Or one of the daughter track has a TOF signal
                                (   (TMath::Abs(fTreeV0VarPosTOFExpTDiff+2500.  )     > 1e-6 )  ||
                                    (TMath::Abs(fTreeV0VarNegTOFExpTDiff+2500.  )     > 1e-6 ) 
                                )
                            )
                            
                        )
                        
                    )//end major if
                    {
                        lV0IsSaved = kTRUE;
                        fTreeV0->Fill();
                    }
                        
                }
                    
            }
            
        }
        
    }
    
    if(IsEvtWithCascade && fSaveRsn)//Adaptation of AliRsnMiniAnalysisTask
    {
        // check current event
        fRsnEvent->SetRef(fInputEvent);
        fRsnEvent->SetPIDResponse(fPIDResponse);
        
        //------------------------------------------------
        // Fill the mini-event
        //------------------------------------------------
        fRsnMiniEvent->Clear();
        // assign event-related values
        fRsnMiniEvent->SetRef(fRsnEvent->GetRef());
        fRsnMiniEvent->Vz()             = lBestPrimaryVtxPos[2];
        fRsnMiniEvent->Spherocity()     = ComputeSpherocity();
        fRsnMiniEvent->Angle()          = ComputeAngle();
        fRsnMiniEvent->Mult()           = fPercV0M;
        fRsnMiniEvent->Tracklets()      = fnTracklets;
        fRsnMiniEvent->ID()             = fDummyTree ? fDummyTree->GetEntries() : 0;
        
        // loop on daughters and assign track-related values
        Int_t nPos = 0;
        Int_t nNeg = 0;
        Int_t nNeu = 0;
        
        Int_t nCuts = fRsnTrackCuts.GetEntries();
        Int_t nPart = fRsnEvent->GetAbsoluteSum();//Sum of tracks, V0s, cascades
        
        AliRsnDaughter cursor;
        AliRsnMiniParticle *miniParticle = 0x0;
        
        Int_t Idx = 0;
        for( Int_t iPart = 0 ; iPart < nPart ; iPart++)//loop on the nbr of particles
        {
            fRsnEvent->SetDaughter(cursor, iPart);
            miniParticle = fRsnMiniEvent->AddParticle();
            miniParticle->CopyDaughter(&cursor);
            //miniParticle->Index() = iPart;
            miniParticle->Index() = Idx;

            AliAODTrack* aodtrack = cursor.Ref2AODtrack();
            if(aodtrack)
            {
                //miniParticle->Index()               = aodtrack->GetID();
                // local implementation of AliRsnMiniParticle::TrackPassesOOBPileupCut
                miniParticle->PassesOOBPileupCut()  = ( aodtrack->GetStatus() & AliAODTrack::kITSrefit                      ) ||
                                                      ( aodtrack->GetTOFExpTDiff(fTreeCascVarMagField, kTRUE) + 2500 > 1e-6 );
            }
            
            for (Int_t iCut = 0; iCut < nCuts; iCut++) 
            {
                AliRsnCutSet *cuts = (AliRsnCutSet *)fRsnTrackCuts[iCut];
                if (cuts->IsSelected(&cursor)) miniParticle->SetCutBit(iCut);
            }

            
            // if a track passes at least one track cut, it is added to the pool
            if (miniParticle->CutBits()) 
            {
                if (miniParticle->Charge() == '+') nPos++;
                else if (miniParticle->Charge() == '-') nNeg++;
                else nNeu++;
                Idx++;
            }
            else 
            {
                TClonesArray &arr = fRsnMiniEvent->Particles();
                // Printf("B %d",arr.GetEntries());
                arr.RemoveAt(arr.GetEntries()-1);
                // Printf("A %d",arr.GetEntries());
            }
            
        }
        
        // get number of accepted tracks
        AliDebugClass(1, Form("Event %6llu: total = %5d, accepted = %4d (pos %4d, neg %4d, neu %4d)", fTreeCascVarEventNumber, nPart, (Int_t)fRsnMiniEvent->Particles().GetEntriesFast(), nPos, nNeg, nNeu));
        
        //------------------------------------------------
        // Run Resonance Finder
        //------------------------------------------------
        Int_t nPair = 0;
        Bool_t lHasAtLeastOneResonance = kFALSE;
        
        Int_t lCutIDrsn;       // ID of cut set that identifies the resonance candidates
        
        for (Int_t i = 0 ; i < fResonanceFinders.GetEntries() ; i++)
        {
            AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
            
            // Get the relevant parameters for resonance finding
            lCutIDrsn       = f->GetResonanceCutID();
            fMotherMass     = f->GetResonanceMass();
            fPairName       = f->GetName();
            
            fCharge[0]      = f->GetCharge(0);
            fCharge[1]      = f->GetCharge(1);
            
            fDaughter[0]    = f->GetDaughter(0);
            fDaughter[1]    = f->GetDaughter(1);
            
            fCutID[0]       = f->GetCutID(0);
            fCutID[1]       = f->GetCutID(1);
            
            // Find all the pairs in the event
            TObjArray lPairList;
            lPairList.SetOwner(kTRUE);
            FillPair(fRsnMiniEvent, fRsnMiniEvent, lPairList, "");
            
            // Gather the pair informations and store it in fTreeRsn
            nPair = lPairList.GetEntriesFast();
            AliRsnMiniPair *lPair = 0x0; // Pair of particles 1 and 2
            for( Int_t iPair = 0 ; iPair < nPair ; iPair++)
            {
                //________________________________________________________________________
                //Get the pair info
                lPair = (AliRsnMiniPair*)lPairList[iPair];
                
                fTreeRsnVarCutIDrsn            = lCutIDrsn;
                
                fTreeRsnVarPx                  = lPair->Sum(kFALSE).X();
                fTreeRsnVarPy                  = lPair->Sum(kFALSE).Y();
                fTreeRsnVarPz                  = lPair->Sum(kFALSE).Z();
                
                fTreeRsnVarInvMass             = lPair->InvMass(kFALSE);
                
                fTreeRsnVarPassesOOBPileupCut  = lPair->PassesOOBPileupCut();
                
                fTreeRsnVarEventNumber         = fTreeCascVarEventNumber;
                
                //------------------------------------------------
                // Fill Tree!
                //------------------------------------------------
                
                fTreeRsn->Fill();
                
                lHasAtLeastOneResonance = kTRUE;
                
            }
            
            lPairList.Delete();
            // message
            AliDebugClass(1, Form("Event %lld: def = '%15s' -- fills = %5d", Entry(), fPairName.Data(), nPair));
        }
        
        // if a resonance is found in the event
        // save the AliRsnMiniEvent for event-mixing
        fDummyVarEventNumber = fTreeCascVarEventNumber;
        if( fSaveRsn && lHasAtLeastOneResonance &&
            ( fComputationType == "MIX" || fComputationType == "ROTATE1" || fComputationType == "ROTATE2" )
        ) fDummyTree->Fill();
        
    }
    
    
    if(IsEvtWithCascade && fSavePrimaries){
        
        Int_t ntracksLoop = lAODevent->GetNumberOfTracks();
        
        for(Int_t iTr = 0 ; iTr < ntracksLoop ; iTr++){
            
            AliAODTrack* primTrack = static_cast<AliAODTrack*>(lAODevent->GetTrack(iTr));
            
            if (!primTrack ) {
                AliWarning("ERROR: Could not retrieve AOD primary track ...");
                continue;
            }
            
            //________________________________________________________________________
            // Primary track cuts
            
            Double_t lPx = 0.; Double_t lPy = 0.; Double_t lPz = 0.;
            lPx  = primTrack->Px() ; lPy  = primTrack->Py() ; lPz  = primTrack->Pz() ;
            
            Float_t lDCA[2]    ; Float_t lCovDCA[3]; 
            primTrack->GetImpactParameters(lDCA, lCovDCA);
            Double_t lDCAxy = 0.; Double_t lDCAz = 0.;
            lDCAxy  = lDCA[0]   ; lDCAz   = GetDCAz(primTrack);
            
            //------------------------------------------------
            // Calculation of the variables related to primaries
            //------------------------------------------------
            fTreePrimVarCharge          = primTrack->Charge();
            fTreePrimVarRapPion         = primTrack->Y( AliAODTrack::kPion );
            fTreePrimVarRapProton       = primTrack->Y( AliAODTrack::kProton );
            fTreePrimVarRapKaon         = primTrack->Y( AliAODTrack::kKaon );
            fTreePrimVarEta             = primTrack->Eta();
            fTreePrimVarTheta           = primTrack->Theta();
            fTreePrimVarPhi             = primTrack->Phi();
            fTreePrimVarPtot            = primTrack->P();
            fTreePrimVarPt              = primTrack->Pt();
            
            fTreePrimVarDCAxyToPV       = lDCAxy;
            fTreePrimVarDCAzToPV        = lDCAz;
        
            //________________________________________________________________________
            // Track quality cuts
            //Least Nbr of Crossed Rows
            fTreePrimVarNbrCrossedRows = primTrack->GetTPCClusterInfo(2,1);
            
            // Least Nbr of clusters      
            fTreePrimVarNbrClusters   = primTrack->GetTPCNcls();
            
            //Compute ratio Crossed Rows / Findable clusters
            //Note: this test avoids division by zero!
            if( primTrack->GetTPCNclsF() <= 0. ) continue;
            fTreePrimVarRatioCrossedRowsOverFindable = fTreePrimVarNbrCrossedRows / ((Double_t)primTrack->GetTPCNclsF());
            
            //Min track length
            fTreePrimVarTrackLength   = GetLengthInActiveZone(primTrack, 2.0, 220.0, lMagField);
            
            //Nbr Of Crossed Rows Over Length
            fTreePrimVarNbrCrossedRowsOverLength = primTrack->GetTPCClusterInfo(2,1)/(fTreePrimVarTrackLength + 1e-5);
            
            // Fraction of shared TPC clusters
            fTreePrimVarFractionSharedTPCClusters = Float_t(primTrack->GetTPCnclsS())/Float_t(fTreePrimVarNbrClusters + 1e-5);
            
            //ITS Chi2 per cluster
            fTreePrimVarITSChi2PerCluster = primTrack->GetITSchi2() / ((Float_t) primTrack->GetITSNcls() + 1e-5);
            
            //TPC Chi2 per cluster
            fTreePrimVarTPCChi2PerCluster = primTrack->GetTPCchi2() / ((Float_t) fTreePrimVarNbrClusters + 1e-5);
            
            //------------------------------------------------
            // TPC dEdx information
            //------------------------------------------------
            fTreePrimVarNSigmaPion      = fPIDResponse->NumberOfSigmasTPC( primTrack, AliPID::kPion );
            fTreePrimVarNSigmaKaon      = fPIDResponse->NumberOfSigmasTPC( primTrack, AliPID::kKaon );
            fTreePrimVarNSigmaProton    = fPIDResponse->NumberOfSigmasTPC( primTrack, AliPID::kProton );
            
            //------------------------------------------------
            // ITS info 
            //------------------------------------------------
            fTreePrimVarITSNSigmaPion   = fPIDResponse->NumberOfSigmasITS( primTrack, AliPID::kPion );
            fTreePrimVarITSNSigmaKaon   = fPIDResponse->NumberOfSigmasITS( primTrack, AliPID::kKaon );
            fTreePrimVarITSNSigmaProton = fPIDResponse->NumberOfSigmasITS( primTrack, AliPID::kProton );
            
            //------------------------------------------------
            // TOF info 
            //------------------------------------------------
            fTreePrimVarTOFNSigmaPion   = fPIDResponse->NumberOfSigmasTOF( primTrack, AliPID::kPion );
            fTreePrimVarTOFNSigmaKaon   = fPIDResponse->NumberOfSigmasTOF( primTrack, AliPID::kKaon );
            fTreePrimVarTOFNSigmaProton = fPIDResponse->NumberOfSigmasTOF( primTrack, AliPID::kProton );
            
            //------------------------------------------------
            // Raw TPC dEdx + PIDForTracking information
            //------------------------------------------------
            //Acquire TPC Signals
            fTreePrimVardEdx = primTrack->GetTPCsignal();
            //Acquire PID For Tracking
            fTreePrimVarPIDForTracking = primTrack->GetPIDForTracking();
            
            //________________________________________________________________________
            // Track status
            fTreePrimVarTrackStatus     = primTrack->GetStatus();
            
            if ((fTreePrimVarTrackStatus & AliAODTrack::kTPCrefit) == 0) {
                AliDebug(1, "Pb / Prim. track has no TPCrefit ... continue!");
                continue;
            }
            
            //________________________________________________________________________
            // Momentum info
            fTreePrimVarPx = lPx ;  fTreePrimVarPy = lPy  ;  fTreePrimVarPz = lPz  ; 
            
            //________________________________________________________________________
            //Check clusters
            fTreePrimVarITSClusters0 = primTrack->HasPointOnITSLayer(0);
            fTreePrimVarITSClusters1 = primTrack->HasPointOnITSLayer(1);
            fTreePrimVarITSClusters2 = primTrack->HasPointOnITSLayer(2);
            fTreePrimVarITSClusters3 = primTrack->HasPointOnITSLayer(3);
            fTreePrimVarITSClusters4 = primTrack->HasPointOnITSLayer(4);
            fTreePrimVarITSClusters5 = primTrack->HasPointOnITSLayer(5);
    
            //________________________________________________________________________
            //Check its clusters, shared
            fTreePrimVarITSSharedClusters0 = primTrack->HasSharedPointOnITSLayer(0);
            fTreePrimVarITSSharedClusters1 = primTrack->HasSharedPointOnITSLayer(1);
            fTreePrimVarITSSharedClusters2 = primTrack->HasSharedPointOnITSLayer(2);
            fTreePrimVarITSSharedClusters3 = primTrack->HasSharedPointOnITSLayer(3);
            fTreePrimVarITSSharedClusters4 = primTrack->HasSharedPointOnITSLayer(4);
            fTreePrimVarITSSharedClusters5 = primTrack->HasSharedPointOnITSLayer(5);
            
            //________________________________________________________________________
            //GetKinkIndex condition
            fTreePrimVarIsKink = kFALSE;
            if( primTrack->GetKinkIndex(0)>0 ) fTreePrimVarIsKink = kTRUE;
            
            //________________________________________________________________________
            //TOF info
            fTreePrimVarTOFExpTDiff     = primTrack->GetTOFExpTDiff( lMagField );
            fTreePrimVarTOFSignal       = primTrack->GetTOFsignal() * 1.e-3; // in ns
            fTreePrimVarTOFBCid         = primTrack->GetTOFBunchCrossing( lMagField );
            
            //________________________________________________________________________
            fTreePrimVarRunNumber = fTreeCascVarRunNumber;
            fTreePrimVarEventNumber = fTreeCascVarEventNumber;
            
            //------------------------------------------------
            // Fill Tree!
            //------------------------------------------------
            
            // Filter
            if( !primTrack->TestFilterBit(5)            ) continue;
            
            if( !fPrimariesSaveAddConfig )
            {
                if( fTreePrimVarCharge             == 0     ) continue;// Charge != 0
                if( fTreePrimVarDCAxyToPV           > 0.5   ) continue;// DCAxy < 0.5
                if( fTreePrimVarDCAzToPV            > 2     ) continue;// DCAz < 2
                if( fTreePrimVarTPCChi2PerCluster   > 4.    ) continue;// TPC Chi2/Cluster < 4
                if( fTreePrimVarPt                  < 0.15  ) continue;// Min momentum
                if( TMath::Abs(fTreePrimVarEta)     > 0.8   ) continue;// |eta| < 0.8
                
                if( fTreePrimVarNbrCrossedRows                  < 70    ) continue; // Nbr of Crossed Rows > 70
                if( fTreePrimVarTrackLength                     < 80    ) continue; // Track length > 80cm
                if( fTreePrimVarRatioCrossedRowsOverFindable    < 0.8   ) continue; // Ratio Crossed Rows/Findable > 0.8
                
                if( TMath::Abs(fTreePrimVarNSigmaPion)      > 4. ) continue;
                if( TMath::Abs(fTreePrimVarNSigmaKaon)      > 4. ) continue;
                if( TMath::Abs(fTreePrimVarNSigmaProton)    > 4. ) continue;
                
                fTreePrimTrack->Fill();
            }
            
            if( fPrimariesSaveAddConfig )
            {
                Double_t lCutMinDCAToVertexXY   = 0.;
                Double_t lCutMinDCAToVertexZ    = 0.;   
                Double_t lCutMaxDCAToVertexXY   = 0.;
                Double_t lCutMaxDCAToVertexZ    = 0.;
                
                lCutMinDCAToVertexXY    = fPrimTrackCuts->GetMinDCAToVertexXY();
                lCutMinDCAToVertexZ     = fPrimTrackCuts->GetMinDCAToVertexZ();
                lCutMaxDCAToVertexXY    = fPrimTrackCuts->GetMaxDCAToVertexXY();
                lCutMaxDCAToVertexZ     = fPrimTrackCuts->GetMaxDCAToVertexZ();
                
                TString lMinDCAXYFormula    ( fPrimTrackCuts->GetMinDCAToVertexXYPtDep()    );
                TString lMinDCAZFormula     ( fPrimTrackCuts->GetMinDCAToVertexZPtDep()     );
                TString lMaxDCAXYFormula    ( fPrimTrackCuts->GetMaxDCAToVertexXYPtDep()    );
                TString lMaxDCAZFormula     ( fPrimTrackCuts->GetMaxDCAToVertexZPtDep()     );
                
                if( lMinDCAXYFormula.CompareTo("") )
                {
                    lMinDCAXYFormula.ReplaceAll("pt","x");
                    TFormula* lCutMinDCAToVertexXYPtDep = new TFormula("lCutMinDCAToVertexXYPtDep", lMinDCAXYFormula.Data());
                    
                    lCutMinDCAToVertexXY = lCutMinDCAToVertexXYPtDep->Eval(fTreePrimVarPt);
                }
                if( lMinDCAZFormula.CompareTo("") )
                {
                    lMinDCAZFormula.ReplaceAll("pt","x");
                    TFormula* lCutMinDCAToVertexZPtDep = new TFormula("lCutMinDCAToVertexZPtDep", lMinDCAZFormula.Data());
                    
                    lCutMinDCAToVertexZ = lCutMinDCAToVertexZPtDep->Eval(fTreePrimVarPt);
                }
                if( lMaxDCAXYFormula.CompareTo(""))
                {
                    lMaxDCAXYFormula.ReplaceAll("pt","x");
                    TFormula* lCutMaxDCAToVertexXYPtDep = new TFormula("lCutMaxDCAToVertexXYPtDep", lMaxDCAXYFormula.Data());
                    
                    lCutMaxDCAToVertexXY = lCutMaxDCAToVertexXYPtDep->Eval(fTreePrimVarPt);
                }
                if( lMaxDCAZFormula.CompareTo("") )
                {
                    lMaxDCAZFormula.ReplaceAll("pt","x");
                    TFormula* lCutMaxDCAToVertexZPtDep = new TFormula("lCutMaxDCAToVertexZPtDep", lMaxDCAZFormula.Data());
                    
                    lCutMaxDCAToVertexZ = lCutMaxDCAToVertexZPtDep->Eval(fTreePrimVarPt);
                }
                
                Float_t lCutMinPtTracks     = 0.;
                Float_t lCutMaxPtTracks     = 0.;
                
                Float_t lCutMinEtaTracks    = 0.;
                Float_t lCutMaxEtaTracks    = 0.;
                
                Float_t lCutMinRapTracks    = 0.;
                Float_t lCutMaxRapTracks    = 0.;
                
                fPrimTrackCuts->GetPtRange  (lCutMinPtTracks    , lCutMaxPtTracks   );
                fPrimTrackCuts->GetEtaRange (lCutMinEtaTracks   , lCutMaxEtaTracks  );
                fPrimTrackCuts->GetRapRange (lCutMinRapTracks   , lCutMaxRapTracks  );
                
                if (
                    //Check 1: Charge selection
                    fTreePrimVarCharge                            != 0                                                          &&
                    
                    //Check 2: Transverse momentum selection
                    lCutMinPtTracks     < fTreePrimVarPt    && fTreePrimVarPt   < lCutMaxPtTracks                               && 
                    
                    //Check 3: Basic Acceptance cuts
                    // Pseudo-rapidity
                    lCutMinEtaTracks    < fTreePrimVarEta   && fTreePrimVarEta  < lCutMaxEtaTracks                              && 
                    // Rapidity
                    (
                        ( lCutMinRapTracks  < fTreePrimVarRapPion   && fTreePrimVarRapPion      < lCutMaxRapTracks ) ||
                        ( lCutMinRapTracks  < fTreePrimVarRapKaon   && fTreePrimVarRapKaon      < lCutMaxRapTracks ) ||
                        ( lCutMinRapTracks  < fTreePrimVarRapProton && fTreePrimVarRapProton    < lCutMaxRapTracks ) 
                    )                                                                                                           &&
    
                    //Check 4: Topological Variables
                     (   
                        (// if DCAToVertex2D is off (default)
                            !fPrimTrackCuts->GetDCAToVertex2D() && 
                            (  // Classical cut on the DCAxy 
                                TMath::Abs( fTreePrimVarDCAxyToPV ) > lCutMinDCAToVertexXY   && 
                                TMath::Abs( fTreePrimVarDCAxyToPV ) < lCutMaxDCAToVertexXY   &&
                                // Classical cut on the DCAz 
                                TMath::Abs( fTreePrimVarDCAzToPV )  > lCutMinDCAToVertexZ    && 
                                TMath::Abs( fTreePrimVarDCAzToPV )  < lCutMaxDCAToVertexZ
                            )
                        ) 
                        ||
                        (// if DCAToVertex2D is on 
                            fPrimTrackCuts->GetDCAToVertex2D()  &&
                            (  // sqrt( (DCAxy/MaxCutDCAxy)**2 + (DCAz/MaxCutDCAz)**2 ) > 1 IF MaxCutDCAxy & MaxCutDCAz > 0
                                ( 
                                    ( lCutMaxDCAToVertexXY <= 0 || lCutMaxDCAToVertexZ <= 0 ) ||
                                    TMath::Sqrt( TMath::Power(fTreePrimVarDCAxyToPV/lCutMaxDCAToVertexXY, 2) + 
                                                 TMath::Power(fTreePrimVarDCAzToPV /lCutMaxDCAToVertexZ , 2) ) > 1 
                                ) &&
                                // sqrt( (DCAxy/MinCutDCAxy)**2 + (DCAz/MinCutDCAz)**2 ) < 1 IF MinCutDCAxy & MinCutDCAz > 0
                                ( 
                                    ( lCutMinDCAToVertexXY <= 0 || lCutMinDCAToVertexZ <= 0 ) ||
                                    TMath::Sqrt( TMath::Power(fTreePrimVarDCAxyToPV/lCutMinDCAToVertexXY, 2) + 
                                                 TMath::Power(fTreePrimVarDCAzToPV /lCutMinDCAToVertexZ , 2) ) < 1 
                                ) 
                            )
                        )                                                                                                       
                     )                                                                                                          &&
                    
                    // - Miscellaneous
                    fTreePrimVarNbrCrossedRows                > fPrimTrackCuts->GetMinNCrossedRowsTPC()                         &&
                    fTreePrimVarRatioCrossedRowsOverFindable  > fPrimTrackCuts->GetMinRatioCrossedRowsOverFindableClustersTPC() &&
                    fTreePrimVarTrackLength                   > fPrimTrackCuts->GetMinLengthActiveVolumeTPC()                   &&
                    fTreePrimVarFractionSharedTPCClusters     < fPrimTrackCuts->GetMaxFractionSharedTPCClusters()               &&
                    
                    //Check 5: TPC dEdx selections
                    TMath::Abs(fTreePrimVarNSigmaPion)        < 4.                                                              &&
                    TMath::Abs(fTreePrimVarNSigmaKaon)        < 4.                                                              &&
                    TMath::Abs(fTreePrimVarNSigmaProton)      < 4.                                                              &&

                    //Check 6: Max Chi2/Clusters 
                    fTreePrimVarITSChi2PerCluster             < fPrimTrackCuts->GetMaxChi2PerClusterITS()                       &&
                    fTreePrimVarTPCChi2PerCluster             < fPrimTrackCuts->GetMaxChi2PerClusterTPC()                       &&
                    
                    //Check 7 : Chi2 TPC Constrained Vs Global cut
                    primTrack->GetChi2TPCConstrainedVsGlobal() < fPrimTrackCuts->GetMaxChi2TPCConstrainedGlobal()               &&
                    
                    //Check 8: kITSrefit track selection if requested
                    (   //do not need ITS refit tracks --> skip this cut
                        !( fPrimTrackCuts->GetRequireITSRefit() )                           ||
                        //otherwise : need to have ITS refit tracks
                        ( fTreePrimVarTrackStatus    & AliAODTrack::kITSrefit ) 
                    )                                                                                                           &&

                    //Check 9: Kink rejection
                    (   //do not reject kink topology            OR  reject them if required
                        fPrimTrackCuts->GetAcceptKinkDaughters() || fTreePrimVarIsKink == kFALSE 
                    )                                                                                                           &&
                    
                    //Check 10 : ITS Cluster requirements
                    ( 
                      CheckITSClusterRequirement( fPrimTrackCuts->GetClusterRequirementITS( AliESDtrackCuts::kSPD ),      
                                                  fTreePrimVarITSClusters0, 
                                                  fTreePrimVarITSClusters1)  &&
                        
                      CheckITSClusterRequirement( fPrimTrackCuts->GetClusterRequirementITS( AliESDtrackCuts::kSDD ),
                                                  fTreePrimVarITSClusters2, 
                                                  fTreePrimVarITSClusters3)  &&
                        
                      CheckITSClusterRequirement( fPrimTrackCuts->GetClusterRequirementITS( AliESDtrackCuts::kSSD ),
                                                  fTreePrimVarITSClusters4, 
                                                  fTreePrimVarITSClusters5) 
                    )                                                                                                           
                    
                )//end major if
                {
                    fTreePrimTrack->Fill();
                }
                
            }
        
            
        }
        
        
    }
    
    // Post output data.
    PostData(1, fOutputList); 
    PostData(2, fTreeCascade);
    PostData(3, fTreeV0);
    PostData(4, fTreeRsn);
    PostData(5, fTreeRsnBkg);                           
    PostData(6, fTreePrimTrack);
    
}

//_____________________________________________________________________________
// 
void AliAnalysisTaskStrangeCascadesTriggerAODRun2::FinishTaskOutput()
{
    if( !fSaveRsn || !( fComputationType == "MIX" || fComputationType == "ROTATE1" || fComputationType == "ROTATE2" ) ) return;
    
    fDummyTree->SetBranchAddress("fRsnMiniEvent", &fRsnMiniEvent);
    fDummyTree->SetBranchAddress("fDummyVarEventNumber", &fDummyVarEventNumber);
    
    Int_t nPair = 0;
    Int_t nEntries = fDummyTree->GetEntries();

    printf("Computing background with %d events...\n", nEntries);
    
    // step 0 : check that everything is ready for event mixing
    if( nEntries == 0 || ( fComputationType == "MIX" && fNMix < 1 ) ){
        printf("Stopping here, Resonance tree empty or NMix equal to 0\n");
        return;
    }
    
    
    if( fComputationType == "MIX" )
    {
        // step 1 : construct the pool of candidates --> search for good matchings
        Int_t iEvt2 = 0;
        Int_t NbrOfMatchedEvt[nEntries];
        TString *MatchedEvt = new TString[nEntries];
        //Initialize mixing event counter
        for(Int_t iEvt = 0 ; iEvt < nEntries ; iEvt++){
            NbrOfMatchedEvt[iEvt] = 0;
            MatchedEvt[iEvt] = "|";
        }
        
        //Be careful about the indices !
        for(Int_t iEvt1 = 0 ; iEvt1 < nEntries ; iEvt1++)//loop on event 1
        {
            fDummyTree->GetEntry(iEvt1);
            // make a copy of the mini-event
            AliRsnMiniEvent *miniEventRef = new AliRsnMiniEvent(*fRsnMiniEvent);
            
            if( NbrOfMatchedEvt[iEvt1] >= fNMix )
            {
                AliDebugClass(1, Form("Matches for event %5d = %d [%s] ", miniEventRef->ID(), NbrOfMatchedEvt[iEvt1], MatchedEvt[iEvt1].Data()));
                continue;
            }
            
            for(Int_t iMix = 1 ; iMix < nEntries ; iMix++)//loop on event 2
            {
                iEvt2 = iEvt1 + iMix;
                
                if( iEvt2 == iEvt1) continue; //we look for other events
                if( iEvt2 >= nEntries ) iEvt2 -= nEntries; //look at events before iEvt1
                
                fDummyTree->GetEntry(iEvt2);
                // skip if events do not match
                if (!EventsMatch(miniEventRef, fRsnMiniEvent)) continue;
                // skip if events already in the pool
                if( MatchedEvt[iEvt1].Contains(Form("|%d|", iEvt2)) ) continue;
                if( MatchedEvt[iEvt2].Contains(Form("|%d|", iEvt1)) ) continue;
                
                //cout <<"event 1 : " << iEvt1 << " ; event 2 : " << iEvt2 << endl;
                //add new mixing candidate
                MatchedEvt[iEvt1].Append(Form("%d|", iEvt2));
                NbrOfMatchedEvt[iEvt1]++;
                //if event1 and event2 match, and event2 dont have the required nbr of events yet, add event 1 to his pool
                if( NbrOfMatchedEvt[iEvt2] < fNMix ){
                    MatchedEvt[iEvt2].Append(Form("%d|", iEvt1)); 
                    NbrOfMatchedEvt[iEvt2]++;
                }
                
                if( NbrOfMatchedEvt[iEvt1] >= fNMix ) break;
                
            }//end loop event2
            
            AliDebugClass(1, Form("Matches for event %5d = %d [%s] ", miniEventRef->ID(), NbrOfMatchedEvt[iEvt1], MatchedEvt[iEvt1].Data()));
            
        }//end loop event1
        
        //step 2 : perform event mixing
        TObjArray* EvtMixingPool = 0x0;
        TObjString *os = 0x0;
        
        for(Int_t iEvt1 = 0 ; iEvt1 < nEntries ; iEvt1++)// loop event 1
        {
            fDummyTree->GetEntry(iEvt1);
            
            AliRsnMiniEvent* miniEventRef = new AliRsnMiniEvent(*fRsnMiniEvent);
            ULong64_t CurrentEvtNbr = fDummyVarEventNumber;
            
            EvtMixingPool = MatchedEvt[iEvt1].Tokenize("|");
            TObjArrayIter EvtMixingPoolIterator(EvtMixingPool);
            
            fTreeRsnFoundMixEvts = EvtMixingPool->GetEntries();
            
            while( ( os = (TObjString*)EvtMixingPoolIterator() ) )// loop on evt 2 = matching evt with evt 1
            {
                iEvt2 = (Int_t)( os->GetString().Atoi() );
                
                fDummyTree->GetEntry(iEvt2);
                
                Bool_t sameEvent = (miniEventRef->ID() == fRsnMiniEvent->ID());
                if( sameEvent )
                {
                    AliDebugClass(1, "Skipping same events");
                    continue;
                }
                
                Int_t lCutIDrsn;
                for (Int_t i = 0 ; i < fResonanceFinders.GetEntries() ; i++)
                {
                    AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
                
                    // Get the relevant parameters for resonance finding
                    lCutIDrsn       = f->GetResonanceCutID();
                    fMotherMass     = f->GetResonanceMass();
                    fPairName       = f->GetName();
                    
                    fCharge[0]      = f->GetCharge(0);
                    fCharge[1]      = f->GetCharge(1);
                    
                    fDaughter[0]    = f->GetDaughter(0);
                    fDaughter[1]    = f->GetDaughter(1);
                    
                    fCutID[0]       = f->GetCutID(0);
                    fCutID[1]       = f->GetCutID(1);
                    
                    Bool_t IsSymmetric = (fCharge[0] == fCharge[1]) && (fDaughter[0] == fDaughter[1]);
                    
                    //
                    // Find all the pairs in the event
                    // if the two daughters are not symmetric (different particle or charge)
                    // Also fill the symmetric pair
                    //
                    TObjArray lPairList;
                    lPairList.SetOwner(kTRUE);
                    //Here fComputationType == "MIX"
                    FillPair(miniEventRef, fRsnMiniEvent, lPairList, fComputationType);
                    if( !IsSymmetric ) FillPair(fRsnMiniEvent, miniEventRef, lPairList, fComputationType);
                    
                    // Gather the pair informations and store it in fTreeRsn
                    nPair = lPairList.GetEntriesFast();
                    AliRsnMiniPair *lPair = 0x0; // Pair of particles 1 and 2
                    for( Int_t iPair = 0 ; iPair < nPair ; iPair++)
                    {
                        //________________________________________________________________________
                        //Get the pair info
                        lPair = (AliRsnMiniPair*)lPairList[iPair];
                        
                        fTreeRsnBkgVarCutIDrsn            = lCutIDrsn;
                        
                        fTreeRsnBkgVarPx                  = lPair->Sum(kFALSE).X();
                        fTreeRsnBkgVarPy                  = lPair->Sum(kFALSE).Y();
                        fTreeRsnBkgVarPz                  = lPair->Sum(kFALSE).Z();
                        
                        fTreeRsnBkgVarInvMass             = lPair->InvMass(kFALSE);
                        
                        fTreeRsnBkgVarPassesOOBPileupCut  = lPair->PassesOOBPileupCut();
                        
                        fTreeRsnBkgVarEventNumber         = CurrentEvtNbr;
                        
                        //------------------------------------------------
                        // Fill Tree!
                        //------------------------------------------------
                        
                        fTreeRsnBkg->Fill();
                        
                    }
                    lPairList.Delete();
                    
                }
                
            }
            
            PostData(5, fTreeRsnBkg); 
            
            EvtMixingPool->Clear();
            
        }
        
        printf("Event mixing : done !\n");
        
        delete [] MatchedEvt;
        
    }
    
    if( fComputationType == "ROTATE1" || fComputationType == "ROTATE2" )
    {
        
        for(Int_t iEvt = 0 ; iEvt < nEntries ; iEvt++)//loop on event 
        {
            fDummyTree->GetEntry(iEvt);
            ULong64_t CurrentEvtNbr = fDummyVarEventNumber;
            
            Int_t lCutIDrsn;
            for (Int_t i = 0 ; i < fResonanceFinders.GetEntries() ; i++)
            {
                AliRsnMiniResonanceFinder* f = (AliRsnMiniResonanceFinder*) fResonanceFinders[i];
                
                // Get the relevant parameters for resonance finding
                lCutIDrsn       = f->GetResonanceCutID();
                fMotherMass     = f->GetResonanceMass();
                fPairName       = f->GetName();
                
                fCharge[0]      = f->GetCharge(0);
                fCharge[1]      = f->GetCharge(1);
                
                fDaughter[0]    = f->GetDaughter(0);
                fDaughter[1]    = f->GetDaughter(1);
                
                fCutID[0]       = f->GetCutID(0);
                fCutID[1]       = f->GetCutID(1);
                
                // Find all the pairs in the event
                TObjArray lPairList;
                lPairList.SetOwner(kTRUE);
                //Here fcompType == "ROTATE1" or fcompType == "ROTATE2"
                FillPair(fRsnMiniEvent, fRsnMiniEvent, lPairList, fComputationType);
                
                // Gather the pair informations and store it in fTreeRsn
                nPair = lPairList.GetEntriesFast();
                AliRsnMiniPair *lPair = 0x0; // Pair of particles 1 and 2
                for( Int_t iPair = 0 ; iPair < nPair ; iPair++)
                {
                    //________________________________________________________________________
                    //Get the pair info
                    lPair = (AliRsnMiniPair*)lPairList[iPair];
                    
                    fTreeRsnBkgVarCutIDrsn            = lCutIDrsn;
                    
                    fTreeRsnBkgVarPx                  = lPair->Sum(kFALSE).X();
                    fTreeRsnBkgVarPy                  = lPair->Sum(kFALSE).Y();
                    fTreeRsnBkgVarPz                  = lPair->Sum(kFALSE).Z();
                    
                    fTreeRsnBkgVarInvMass             = lPair->InvMass(kFALSE);
                    
                    fTreeRsnBkgVarPassesOOBPileupCut  = lPair->PassesOOBPileupCut();
                    
                    fTreeRsnBkgVarEventNumber         = CurrentEvtNbr;
                    
                    //------------------------------------------------
                    // Fill Tree!
                    //------------------------------------------------
                    
                    fTreeRsnBkg->Fill();
                    
                }
                
                lPairList.Delete();
                
            }
            
            PostData(5, fTreeRsnBkg); 
            
        }
        
        printf("Rotate1 : done !\n");
        
    }
    
}

//_____________________________________________________________________________
void AliAnalysisTaskStrangeCascadesTriggerAODRun2::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
    
    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        AliError("ERROR - AliAnalysisTaskStrangeCascadesTriggerAODRun2 : ouput data container list not available\n");
        return;
    }
    
    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        AliError("ERROR - AliAnalysisTaskStrangeCascadesTriggerAODRun2 : fHistEventCounter not available");
        return;
    }
    
    TCanvas *canCheck = new TCanvas("AliAnalysisTaskStrangeCascadesTriggerAODRun2","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();
    
    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskStrangeCascadesTriggerAODRun2::GetLengthInActiveZone( AliAODTrack *gt, Float_t deltaY, Float_t deltaZ, Float_t b ){
    AliESDtrack esdTrack( gt );
    esdTrack.SetESDEvent((AliESDEvent*) gt->GetEvent() );
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(gt);
    esdTrack.ResetTrackParamIp(&etp);
    return esdTrack.GetLengthInActiveZone(1, deltaY, deltaZ, b);
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskStrangeCascadesTriggerAODRun2::GetDCAz(AliAODTrack *lTrack)
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


//_____________________________________________________________________________
Double_t AliAnalysisTaskStrangeCascadesTriggerAODRun2::ComputeSpherocity()
//Encapsulation of spherocity calculation
{
    AliVEvent * evTypeS = InputEvent();
    Int_t ntracksLoop = evTypeS->GetNumberOfTracks();
    Int_t GoodTracks = 0;
    Float_t pFull = 0;
    Float_t Spherocity = 2;
    Float_t pt[ntracksLoop],phi[ntracksLoop];
    
    //computing total pt
    Float_t sumapt = 0;
    for(Int_t i1 = 0; i1 < ntracksLoop; i1++){
        AliVTrack   *track = (AliVTrack *)evTypeS->GetTrack(i1);
        AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
        //AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
        if (aodt) if (!aodt->TestFilterBit(5)) continue;
        //if (esdt) if (!fTrackFilter->IsSelected(esdt)) continue;
        if (track->Pt() < 0.15) continue;
        if(TMath::Abs(track->Eta()) > 0.8) continue;
        //pt[i1] = track->Pt();
        pt[i1] = 1.0;
        sumapt += pt[i1];
        GoodTracks++;
    }
    if (GoodTracks <= 9) return -10.0;
    //Getting thrust
    for(Int_t i = 0; i < 360/0.1; i++){
        Float_t numerador = 0;
        Float_t phiparam  = 0;
        Float_t nx = 0;
        Float_t ny = 0;
        phiparam=( (TMath::Pi()) * i * 0.1 ) / 180; // parametrization of the angle
        nx = TMath::Cos(phiparam);            // x component of an unitary vector n
        ny = TMath::Sin(phiparam);            // y component of an unitary vector n
        for(Int_t i1 = 0; i1 < ntracksLoop; i1++){
            AliVTrack   *track = (AliVTrack *)evTypeS->GetTrack(i1);
            AliAODTrack *aodt  = dynamic_cast<AliAODTrack *>(track);
            //AliESDtrack *esdt  = dynamic_cast<AliESDtrack *>(track);
            if (aodt) if (!aodt->TestFilterBit(5)) continue;
            //if (esdt) if (!fTrackFilter->IsSelected(esdt)) continue;
            if (track->Pt() < 0.15) continue;
            if(TMath::Abs(track->Eta()) > 0.8) continue;
            //pt[i1] = track->Pt();
            pt[i1] = 1.0;
            phi[i1] = track->Phi();
            Float_t pxA = pt[i1] * TMath::Cos( phi[i1] );
            Float_t pyA = pt[i1] * TMath::Sin( phi[i1] );
            numerador += TMath::Abs( ny * pxA - nx * pyA );//product between p  proyection in XY plane and the unitary vector
        }
        pFull=TMath::Power( (numerador /(sumapt + 1e-6) ), 2);
        if(pFull < Spherocity)//maximization of pFull
        {
            Spherocity = pFull;
        }
    }
    if (GoodTracks > 9) return ((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;
    else return -10.0;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskStrangeCascadesTriggerAODRun2::ComputeAngle()
//Encapsulation of event plane angle calculation
{
//
// Get the plane angle
//

   AliEventplane *plane = 0x0;

   if (fInputEvent->InheritsFrom(AliESDEvent::Class()))
      plane = fInputEvent->GetEventplane();
   else if (fInputEvent->InheritsFrom(AliAODEvent::Class())) {
      AliAODEvent *aodEvent = (AliAODEvent *)fInputEvent;
      plane = ((AliVAODHeader*)aodEvent->GetHeader())->GetEventplaneP();
   }

   if (plane)
      return plane->GetEventplane("Q");
   else {
      AliWarning("No event plane defined");
      return 1E20;
   }
}

//_____________________________________________________________________________
/// Check if events are compatible 
/// This is true if differences in vz, mult and angle are smaller than
/// the specified values.
Bool_t AliAnalysisTaskStrangeCascadesTriggerAODRun2::EventsMatch(AliRsnMiniEvent* event1, AliRsnMiniEvent* event2)
{
    if (!event1 || !event2) return kFALSE;
    Double_t dVz, dMult, dAngle;
    
    dVz     = TMath::Abs( event1->Vz()      - event2->Vz()      );
    dMult   = TMath::Abs( event1->Mult()    - event2->Mult()    );
    dAngle  = TMath::Abs( event1->Angle()   - event2->Angle()   );
    
    if(dVz      > fMaxDiffVz    ) return kFALSE;
    if(dMult    > fMaxDiffMult  ) return kFALSE;
    if(dAngle   > fMaxDiffAngle ) return kFALSE;
    
    return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskStrangeCascadesTriggerAODRun2::MyPhi(Double_t lPx, Double_t lPy) const
{
    // Local calculation for azimuthal angle
    Double_t phi = 0;
    phi = TMath::ATan2(lPy,lPx);
    if( phi < 0 ) phi += 2*TMath::Pi();
    return phi;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskStrangeCascadesTriggerAODRun2::MyTheta(Double_t lPt, Double_t lPz) const
{
    // Local calculation for azimuthal angle
    Double_t theta = 0.;
    theta = TMath::ATan( lPt/lPz );
    if( theta < 0 ) theta += TMath::Pi();
    return theta;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskStrangeCascadesTriggerAODRun2::CheckITSClusterRequirement(AliESDtrackCuts::ITSClusterRequirement req, 
                                                         Bool_t clusterL1, Bool_t clusterL2)
{
    //
    // Encapsulation of ITS cluster requirement 
    // (from https://github.com/alisw/AliRoot/blob/master/ANALYSIS/ANALYSISalice/AliESDtrackCuts.cxx#L2268)
    //
    // checks if the cluster requirement is fullfilled (in this case: return kTRUE)
    //
    switch (req)
    {
        case AliESDtrackCuts::kOff          :  return kTRUE;
        case AliESDtrackCuts::kNone         :  return !clusterL1 && !clusterL2;
        case AliESDtrackCuts::kAny          :  return clusterL1 || clusterL2;
        case AliESDtrackCuts::kFirst        :  return clusterL1;
        case AliESDtrackCuts::kOnlyFirst    :  return clusterL1 && !clusterL2;
        case AliESDtrackCuts::kSecond       :  return clusterL2;
        case AliESDtrackCuts::kOnlySecond   :  return clusterL2 && !clusterL1;
        case AliESDtrackCuts::kBoth         :  return clusterL1 && clusterL2;
    }

    return kFALSE;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskStrangeCascadesTriggerAODRun2::FillPair(AliRsnMiniEvent *miniEventRef, AliRsnMiniEvent *fRsnMiniEvent, TObjArray &lPairList, TString lcompType)
{
    //
    // Loops on the passed mini-event, and for each pair of particles
    // which satisfy the charge and cut requirements defined here, add an entry.
    // Returns the number of successful fillings.
    // Last argument tells if the type of pairing :
    // -- ""        --> pair particles
    // -- "ROTATE1" --> pair particles but rotate first particle
    // -- "ROTATE2" --> pair particles but rotate second particle
    
    // check computation type
    Bool_t okComp = kFALSE;
    if (lcompType == "")         okComp = kTRUE;
    if (lcompType == "MIX")      okComp = kTRUE;
    if (lcompType == "ROTATE1")  okComp = kTRUE;
    if (lcompType == "ROTATE2")  okComp = kTRUE;
    if (!okComp) {
        AliError(Form("[%s] Unknown computation type", GetName()));
        //return 0;
    }
    
    Int_t nadded = 0;
    
    AliRsnMiniParticle *p1  = 0x0;  // daughter part 1
    AliRsnMiniParticle *p2  = 0x0;  // daughter part 2
    Double_t m1             = AliRsnDaughter::SpeciesMass(fDaughter[0]);   // mass of part 1
    Double_t m2             = AliRsnDaughter::SpeciesMass(fDaughter[1]);   // mass of part 2
    
    AliRsnMiniPair *lPair = new AliRsnMiniPair(); // Pair of particles 1 and 2
    Int_t start;
    
    Bool_t sameCriteria = ((fCharge[0] == fCharge[1]) && (fDaughter[0] == fDaughter[1]));
    Bool_t sameEvent = (miniEventRef->ID() == fRsnMiniEvent->ID());
    
    TString selList1  = "";
    TString selList2  = "";
    TArrayI          fSel1;             // list of selected particles for definition 1
    TArrayI          fSel2;             // list of selected particles for definition 2            
    Int_t   n1 = miniEventRef->CountParticles(fSel1, fCharge[0], fCutID[0]);
    Int_t   n2 = fRsnMiniEvent->CountParticles(fSel2, fCharge[1], fCutID[1]);
    
    for(Int_t i1 = 0; i1 < n1; i1++) selList1.Append(Form("%d ", fSel1[i1]));
    for(Int_t i2 = 0; i2 < n2; i2++) selList2.Append(Form("%d ", fSel2[i2]));
    AliDebugClass(1, Form("[%10s] Part #1: [%s] -- evID %6d -- charge = %c -- cut ID = %d --> %4d tracks (%s)", fPairName.Data(), (miniEventRef == fRsnMiniEvent ? "def" : "mix"), miniEventRef->ID(), fCharge[0], fCutID[0], n1, selList1.Data()));
    AliDebugClass(1, Form("[%10s] Part #2: [%s] -- evID %6d -- charge = %c -- cut ID = %d --> %4d tracks (%s)", fPairName.Data(), (miniEventRef == fRsnMiniEvent ? "def" : "mix"), fRsnMiniEvent->ID(), fCharge[1], fCutID[1], n2, selList2.Data()));
    if (!n1 || !n2) {
        AliDebugClass(1, "No pairs to mix");
        return 0;
    }
    
    // external loop
    for(Int_t i1 = 0 ; i1 < n1 ; i1++) 
    {
        p1 = miniEventRef->GetParticle(fSel1[i1]);
        // define starting point for inner loop
        // if daughter selection criteria (charge, cuts) are the same
        // and the two events coincide, internal loop must start from
        // the first track *after* current one;
        // otherwise it starts from the beginning
        start = ((sameEvent && sameCriteria) ? i1 + 1 : 0);
        AliDebugClass(2, Form("Start point = %d", start));
        // internal loop
        for(Int_t i2 = start ; i2 < n2 ; i2++) 
        {
            p2 = fRsnMiniEvent->GetParticle(fSel2[i2]);
            
            // avoid to mix a particle with itself
            if ( sameEvent && (p1->Index() == p2->Index()) ) 
            {
                AliDebugClass(2, "Skipping same index");
                continue;
            }
            
            // avoid to mix resonance particles
            if( p1->IsResonance() || p2->IsResonance() )
            {
                AliDebugClass(2, "Skipping resonance particles");
                continue;
            }
            
            lPair->Fill(p1, p2, m1, m2, fMotherMass);
            
            // do rotation if needed
            if( lcompType == "ROTATE1" ) lPair->InvertP(kTRUE);
            if( lcompType == "ROTATE2" ) lPair->InvertP(kFALSE);
           
            // check pair against cuts
            if (fRsnPairCuts && !fRsnPairCuts->IsSelected(lPair)) continue;
            
            lPairList.AddLast( new AliRsnMiniPair() );
            ( (AliRsnMiniPair*)lPairList.Last() )->Fill(p1, p2, m1, m2, fMotherMass);
            // do rotation if needed
            if( lcompType == "ROTATE1" ) ( (AliRsnMiniPair*)lPairList.Last() )->InvertP(kTRUE);
            if( lcompType == "ROTATE2" ) ( (AliRsnMiniPair*)lPairList.Last() )->InvertP(kFALSE);
            nadded++;
        }
        
    }
    
    return nadded;
}

