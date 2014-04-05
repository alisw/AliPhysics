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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// This task is meant to explore the possibility of using a VZERO amplitude
// based multiplicity estimator for proton-proton collisions. For this, two 
// main operation methods for this task are foreseen: 
// 
//  1) (under development) it should act as an auxiliary task and provide a 
//     calibrated estimator 
//
//  2) "Debug mode" which will also create a ROOT TTree object with event 
//     by event info potentially used for exploration / calibration. This 
//     includes the following info: 
//    
//      --- All VZERO Amplitudes (saved as Float_t) 
//      --- (optional) time for each channel
//      --- (optional) time width for each channel 
//      --- GetReferenceMultiplicity Estimator, |eta|<0.5 
//      --- GetReferenceMultiplicity Estimator, |eta|<0.8 
//      --- (if MC) True Multiplicity, |eta|<0.5
//      --- (if MC) True Multiplicity,  2.8 < eta < 5.8 (VZEROA region)
//      --- (if MC) True Multiplicity, -3.7 < eta <-1.7 (VZEROC region)
//      --- Run Number
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

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
//#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskStrangenessVsMultiplicity.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStrangenessVsMultiplicity)

AliAnalysisTaskStrangenessVsMultiplicity::AliAnalysisTaskStrangenessVsMultiplicity()
  : AliAnalysisTaskSE(), fListHist(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), 
  fkSaveV0Tree      ( kFALSE ),
  fkSaveCascadeTree ( kTRUE  ),
  fkRunVertexers    ( kTRUE  ), 
  //---> Variables for fTreeEvent
  fAmplitude_V0A   (0),   
  fAmplitude_V0C   (0),   
  fAmplitude_V0AEq (0),   
  fAmplitude_V0CEq (0),  
  fCentrality_V0A(0), 
  fCentrality_V0C(0), 
  fCentrality_V0M(0), 
  fCentrality_V0AEq(0), 
  fCentrality_V0CEq(0), 
  fCentrality_V0MEq(0), 
  fRefMultEta5(0),
  fRefMultEta8(0),
  fRunNumber(0),
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
	
  fTreeVariableCentV0A(0),
  fTreeVariableCentV0C(0),
  fTreeVariableCentV0M(0),
  fTreeVariableCentV0AEq(0),
  fTreeVariableCentV0CEq(0),
  fTreeVariableCentV0MEq(0),
  fTreeVariableAmpV0A(0),
  fTreeVariableAmpV0C(0),
  fTreeVariableAmpV0AEq(0),
  fTreeVariableAmpV0CEq(0),
  fTreeVariableRefMultEta8(0),
  fTreeVariableRefMultEta5(0),
  fTreeVariableRunNumber(0),
  //---> Variables for fTreeCascade
  fTreeCascVarCharge(0), 
	fTreeCascVarMassAsXi(0),
	fTreeCascVarMassAsOmega(0),
	fTreeCascVarPt(0),
	fTreeCascVarRapXi(0),
	fTreeCascVarRapOmega(0),
	fTreeCascVarNegEta(0),
	fTreeCascVarPosEta(0),
	fTreeCascVarBachEta(0),
	fTreeCascVarDCACascDaughters(0),
	fTreeCascVarDCABachToPrimVtx(0),
	fTreeCascVarDCAV0Daughters(0),
	fTreeCascVarDCAV0ToPrimVtx(0),
	fTreeCascVarDCAPosToPrimVtx(0),
	fTreeCascVarDCANegToPrimVtx(0),
	fTreeCascVarCascCosPointingAngle(0),
	fTreeCascVarCascRadius(0),
	fTreeCascVarV0Mass(0),
	fTreeCascVarV0CosPointingAngle(0),
	fTreeCascVarV0CosPointingAngleSpecial(0),
	fTreeCascVarV0Radius(0),
  fTreeCascVarLeastNbrClusters(0),
	fTreeCascVarDistOverTotMom(0),
	fTreeCascVarNegNSigmaPion(0),
	fTreeCascVarNegNSigmaProton(0),
	fTreeCascVarPosNSigmaPion(0),
	fTreeCascVarPosNSigmaProton(0),
	fTreeCascVarBachNSigmaPion(0),
	fTreeCascVarBachNSigmaKaon(0),
	fTreeCascVarCentV0A(0),
	fTreeCascVarCentV0C(0),
	fTreeCascVarCentV0M(0),
	fTreeCascVarCentV0AEq(0),
	fTreeCascVarCentV0CEq(0),
	fTreeCascVarCentV0MEq(0),
	fTreeCascVarAmpV0A(0),
	fTreeCascVarAmpV0C(0),
	fTreeCascVarAmpV0AEq(0),
	fTreeCascVarAmpV0CEq(0),
  fTreeCascVarRefMultEta8(0),
  fTreeCascVarRefMultEta5(0),
  fTreeCascVarRunNumber(0), 
  //Histos 
  fHistEventCounter(0)
//------------------------------------------------
// Tree Variables 
{

}

AliAnalysisTaskStrangenessVsMultiplicity::AliAnalysisTaskStrangenessVsMultiplicity(const char *name) 
  : AliAnalysisTaskSE(name), fListHist(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), 
  fkSaveV0Tree      ( kFALSE ),
  fkSaveCascadeTree ( kTRUE  ), 
  fkRunVertexers    ( kTRUE  ),
  //---> Variables for fTreeEvent
  fAmplitude_V0A (0),   
  fAmplitude_V0C (0), 
  fAmplitude_V0AEq (0),   
  fAmplitude_V0CEq (0), 
  fCentrality_V0A(0), 
  fCentrality_V0C(0), 
  fCentrality_V0M(0), 
  fCentrality_V0AEq(0), 
  fCentrality_V0CEq(0), 
  fCentrality_V0MEq(0), 
  fRefMultEta5(0),
  fRefMultEta8(0),
  fRunNumber(0),
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
	
  fTreeVariableCentV0A(0),
  fTreeVariableCentV0C(0),
  fTreeVariableCentV0M(0),
  fTreeVariableCentV0AEq(0),
  fTreeVariableCentV0CEq(0),
  fTreeVariableCentV0MEq(0),
  fTreeVariableAmpV0A(0),
  fTreeVariableAmpV0C(0),
  fTreeVariableAmpV0AEq(0),
  fTreeVariableAmpV0CEq(0),
  fTreeVariableRefMultEta8(0),
  fTreeVariableRefMultEta5(0),
  fTreeVariableRunNumber(0),
  //---> Variables for fTreeCascade
  fTreeCascVarCharge(0), 
	fTreeCascVarMassAsXi(0),
	fTreeCascVarMassAsOmega(0),
	fTreeCascVarPt(0),
	fTreeCascVarRapXi(0),
	fTreeCascVarRapOmega(0),
	fTreeCascVarNegEta(0),
	fTreeCascVarPosEta(0),
	fTreeCascVarBachEta(0),
	fTreeCascVarDCACascDaughters(0),
	fTreeCascVarDCABachToPrimVtx(0),
	fTreeCascVarDCAV0Daughters(0),
	fTreeCascVarDCAV0ToPrimVtx(0),
	fTreeCascVarDCAPosToPrimVtx(0),
	fTreeCascVarDCANegToPrimVtx(0),
	fTreeCascVarCascCosPointingAngle(0),
	fTreeCascVarCascRadius(0),
	fTreeCascVarV0Mass(0),
	fTreeCascVarV0CosPointingAngle(0),
	fTreeCascVarV0CosPointingAngleSpecial(0),
	fTreeCascVarV0Radius(0),
  fTreeCascVarLeastNbrClusters(0),
	fTreeCascVarDistOverTotMom(0),
	fTreeCascVarNegNSigmaPion(0),
	fTreeCascVarNegNSigmaProton(0),
	fTreeCascVarPosNSigmaPion(0),
	fTreeCascVarPosNSigmaProton(0),
	fTreeCascVarBachNSigmaPion(0),
	fTreeCascVarBachNSigmaKaon(0),
	fTreeCascVarCentV0A(0),
	fTreeCascVarCentV0C(0),
	fTreeCascVarCentV0M(0),
	fTreeCascVarCentV0AEq(0),
	fTreeCascVarCentV0CEq(0),
	fTreeCascVarCentV0MEq(0),
	fTreeCascVarAmpV0A(0),
	fTreeCascVarAmpV0C(0),
	fTreeCascVarAmpV0AEq(0),
	fTreeCascVarAmpV0CEq(0),
  fTreeCascVarRefMultEta8(0),
  fTreeCascVarRefMultEta5(0),
  fTreeCascVarRunNumber(0), 
  //Histos 
  fHistEventCounter(0)
{

  //Re-vertex: Will only apply for cascade candidates

  fV0VertexerSels[0] =  33.  ;  // max allowed chi2
  fV0VertexerSels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
  fV0VertexerSels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
  fV0VertexerSels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
  fV0VertexerSels[5] =   1.0 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
  fV0VertexerSels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
        
  fCascadeVertexerSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
  fCascadeVertexerSels[1] =   0.05 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
  fCascadeVertexerSels[2] =   0.010;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
  fCascadeVertexerSels[3] =   0.03 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
  fCascadeVertexerSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
  fCascadeVertexerSels[5] =   0.95 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
  fCascadeVertexerSels[6] =   0.4  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
  fCascadeVertexerSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
        

  DefineOutput(1, TList::Class()); // Event Counter Histo
  DefineOutput(2, TTree::Class()); // Event Tree
  DefineOutput(3, TTree::Class()); // V0 Tree
  DefineOutput(4, TTree::Class()); // Cascade Tree
}


AliAnalysisTaskStrangenessVsMultiplicity::~AliAnalysisTaskStrangenessVsMultiplicity()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------

   if (fListHist){
      delete fListHist;
      fListHist = 0x0;
   }
   if (fTreeEvent){
      delete fTreeEvent;
      fTreeEvent = 0x0;
   }
   if (fTreeV0){
      delete fTreeV0;
      fTreeV0 = 0x0;
   }
   if (fTreeCascade){
      delete fTreeCascade;
      fTreeCascade = 0x0;
   }
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicity::UserCreateOutputObjects()
{

   OpenFile(2);	
   // Called once

//------------------------------------------------

   fTreeEvent = new TTree("fTreeEvent","Event");

//------------------------------------------------
// fTree Branch definitions - Event by Event info
//------------------------------------------------

//-----------BASIC-INFO---------------------------

  //--- VZERO Data (Integrated)
  fTreeEvent->Branch("fAmplitude_V0A",&fAmplitude_V0A,"fAmplitude_V0A/F");
  fTreeEvent->Branch("fAmplitude_V0C",&fAmplitude_V0C,"fAmplitude_V0C/F");
  fTreeEvent->Branch("fAmplitude_V0AEq",&fAmplitude_V0AEq,"fAmplitude_V0AEq/F");
  fTreeEvent->Branch("fAmplitude_V0CEq",&fAmplitude_V0CEq,"fAmplitude_V0CEq/F");

  //Info from AliCentrality (not necessarily 'centrality' per se) 
  fTreeEvent->Branch("fCentrality_V0A",&fCentrality_V0A,"fCentrality_V0A/F");
  fTreeEvent->Branch("fCentrality_V0C",&fCentrality_V0C,"fCentrality_V0C/F");
  fTreeEvent->Branch("fCentrality_V0M",&fCentrality_V0A,"fCentrality_V0M/F");
  fTreeEvent->Branch("fCentrality_V0AEq",&fCentrality_V0AEq,"fCentrality_V0AEq/F");
  fTreeEvent->Branch("fCentrality_V0CEq",&fCentrality_V0CEq,"fCentrality_V0CEq/F");
  fTreeEvent->Branch("fCentrality_V0MEq",&fCentrality_V0AEq,"fCentrality_V0MEq/F");
  
  //Official GetReferenceMultiplicity
  fTreeEvent->Branch("fRefMultEta5",&fRefMultEta5,"fRefMultEta5/I");
  fTreeEvent->Branch("fRefMultEta8",&fRefMultEta8,"fRefMultEta8/I");

  //Run Number
  fTreeEvent->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");

  //Create Basic V0 Output Tree
  fTreeV0 = new TTree( "fTreeV0", "V0 Candidates");

//------------------------------------------------
// fTreeV0 Branch definitions
//------------------------------------------------

//-----------BASIC-INFO---------------------------
  fTreeV0->Branch("fTreeVariableChi2V0",&fTreeVariableChi2V0,"fTreeVariableChi2V0/F");
  fTreeV0->Branch("fTreeVariableDcaV0Daughters",&fTreeVariableDcaV0Daughters,"fTreeVariableDcaV0Daughters/F");
  fTreeV0->Branch("fTreeVariableDcaV0ToPrimVertex",&fTreeVariableDcaV0ToPrimVertex,"fTreeVariableDcaV0ToPrimVertex/F");
  fTreeV0->Branch("fTreeVariableDcaPosToPrimVertex",&fTreeVariableDcaPosToPrimVertex,"fTreeVariableDcaPosToPrimVertex/F");
  fTreeV0->Branch("fTreeVariableDcaNegToPrimVertex",&fTreeVariableDcaNegToPrimVertex,"fTreeVariableDcaNegToPrimVertex/F");
  fTreeV0->Branch("fTreeVariableV0Radius",&fTreeVariableV0Radius,"fTreeVariableV0Radius/F");
  fTreeV0->Branch("fTreeVariablePt",&fTreeVariablePt,"fTreeVariablePt/F");
  fTreeV0->Branch("fTreeVariableRapK0Short",&fTreeVariableRapK0Short,"fTreeVariableRapK0Short/F");
  fTreeV0->Branch("fTreeVariableRapLambda",&fTreeVariableRapLambda,"fTreeVariableRapLambda/F");
  fTreeV0->Branch("fTreeVariableInvMassK0s",&fTreeVariableInvMassK0s,"fTreeVariableInvMassK0s/F");
  fTreeV0->Branch("fTreeVariableInvMassLambda",&fTreeVariableInvMassLambda,"fTreeVariableInvMassLambda/F");
  fTreeV0->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
  fTreeV0->Branch("fTreeVariableV0CosineOfPointingAngle",&fTreeVariableV0CosineOfPointingAngle,"fTreeVariableV0CosineOfPointingAngle/F");
  fTreeV0->Branch("fTreeVariableAlphaV0",&fTreeVariableAlphaV0,"fTreeVariableAlphaV0/F");
  fTreeV0->Branch("fTreeVariablePtArmV0",&fTreeVariablePtArmV0,"fTreeVariablePtArmV0/F");
  fTreeV0->Branch("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastNbrCrossedRows,"fTreeVariableLeastNbrCrossedRows/I");
  fTreeV0->Branch("fTreeVariableLeastRatioCrossedRowsOverFindable",&fTreeVariableLeastRatioCrossedRowsOverFindable,"fTreeVariableLeastRatioCrossedRowsOverFindable/F");
  fTreeV0->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
  fTreeV0->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
  fTreeV0->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
  fTreeV0->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
  fTreeV0->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
  fTreeV0->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
  fTreeV0->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
//-----------MULTIPLICITY-INFO--------------------
  fTreeV0->Branch("fTreeVariableCentV0A",&fTreeVariableCentV0A,"fTreeVariableCentV0A/F");
  fTreeV0->Branch("fTreeVariableCentV0C",&fTreeVariableCentV0C,"fTreeVariableCentV0C/F");
  fTreeV0->Branch("fTreeVariableCentV0M",&fTreeVariableCentV0M,"fTreeVariableCentV0M/F");
  fTreeV0->Branch("fTreeVariableCentV0AEq",&fTreeVariableCentV0AEq,"fTreeVariableCentV0AEq/F");
  fTreeV0->Branch("fTreeVariableCentV0CEq",&fTreeVariableCentV0CEq,"fTreeVariableCentV0CEq/F");
  fTreeV0->Branch("fTreeVariableCentV0MEq",&fTreeVariableCentV0MEq,"fTreeVariableCentV0MEq/F");
  fTreeV0->Branch("fTreeVariableAmpV0A",&fTreeVariableAmpV0A,"fTreeVariableAmpV0A/F");
  fTreeV0->Branch("fTreeVariableAmpV0C",&fTreeVariableAmpV0C,"fTreeVariableAmpV0C/F");
  fTreeV0->Branch("fTreeVariableAmpV0AEq",&fTreeVariableAmpV0AEq,"fTreeVariableAmpV0AEq/F");
  fTreeV0->Branch("fTreeVariableAmpV0CEq",&fTreeVariableAmpV0CEq,"fTreeVariableAmpV0CEq/F");
  fTreeV0->Branch("fTreeVariableRefMultEta8",&fTreeVariableRefMultEta8,"fTreeVariableRefMultEta8/I");
  fTreeV0->Branch("fTreeVariableRefMultEta5",&fTreeVariableRefMultEta5,"fTreeVariableRefMultEta5/I");
  fTreeV0->Branch("fTreeVariableRunNumber",&fTreeVariableRunNumber,"fTreeVariableRunNumber/I");
//------------------------------------------------

  //Create Cascade output tree
  fTreeCascade = new TTree("fTreeCascade","CascadeCandidates");

//------------------------------------------------
// fTreeCascade Branch definitions - Cascade Tree
//------------------------------------------------

//-----------BASIC-INFO---------------------------
	fTreeCascade->Branch("fTreeCascVarCharge",&fTreeCascVarCharge,"fTreeCascVarCharge/I");	
  fTreeCascade->Branch("fTreeCascVarMassAsXi",&fTreeCascVarMassAsXi,"fTreeCascVarMassAsXi/F");
  fTreeCascade->Branch("fTreeCascVarMassAsOmega",&fTreeCascVarMassAsOmega,"fTreeCascVarMassAsOmega/F");
  fTreeCascade->Branch("fTreeCascVarPt",&fTreeCascVarPt,"fTreeCascVarPt/F");
  fTreeCascade->Branch("fTreeCascVarRapXi",&fTreeCascVarRapXi,"fTreeCascVarRapXi/F");
  fTreeCascade->Branch("fTreeCascVarRapOmega",&fTreeCascVarRapOmega,"fTreeCascVarRapOmega/F");
  fTreeCascade->Branch("fTreeCascVarNegEta",&fTreeCascVarNegEta,"fTreeCascVarNegEta/F");
  fTreeCascade->Branch("fTreeCascVarPosEta",&fTreeCascVarPosEta,"fTreeCascVarPosEta/F");
  fTreeCascade->Branch("fTreeCascVarBachEta",&fTreeCascVarBachEta,"fTreeCascVarBachEta/F");
//-----------INFO-FOR-CUTS------------------------
  fTreeCascade->Branch("fTreeCascVarDCACascDaughters",&fTreeCascVarDCACascDaughters,"fTreeCascVarDCACascDaughters/F");
  fTreeCascade->Branch("fTreeCascVarDCABachToPrimVtx",&fTreeCascVarDCABachToPrimVtx,"fTreeCascVarDCABachToPrimVtx/F");
  fTreeCascade->Branch("fTreeCascVarDCAV0Daughters",&fTreeCascVarDCAV0Daughters,"fTreeCascVarDCAV0Daughters/F");
  fTreeCascade->Branch("fTreeCascVarDCAV0ToPrimVtx",&fTreeCascVarDCAV0ToPrimVtx,"fTreeCascVarDCAV0ToPrimVtx/F");
  fTreeCascade->Branch("fTreeCascVarDCAPosToPrimVtx",&fTreeCascVarDCAPosToPrimVtx,"fTreeCascVarDCAPosToPrimVtx/F");
  fTreeCascade->Branch("fTreeCascVarDCANegToPrimVtx",&fTreeCascVarDCANegToPrimVtx,"fTreeCascVarDCANegToPrimVtx/F");
  fTreeCascade->Branch("fTreeCascVarCascCosPointingAngle",&fTreeCascVarCascCosPointingAngle,"fTreeCascVarCascCosPointingAngle/F");
  fTreeCascade->Branch("fTreeCascVarCascRadius",&fTreeCascVarCascRadius,"fTreeCascVarCascRadius/F");
  fTreeCascade->Branch("fTreeCascVarV0Mass",&fTreeCascVarV0Mass,"fTreeCascVarV0Mass/F");
  fTreeCascade->Branch("fTreeCascVarV0CosPointingAngle",&fTreeCascVarV0CosPointingAngle,"fTreeCascVarV0CosPointingAngle/F");
  fTreeCascade->Branch("fTreeCascVarV0CosPointingAngleSpecial",&fTreeCascVarV0CosPointingAngleSpecial,"fTreeCascVarV0CosPointingAngleSpecial/F");
  fTreeCascade->Branch("fTreeCascVarV0Radius",&fTreeCascVarV0Radius,"fTreeCascVarV0Radius/F");
  fTreeCascade->Branch("fTreeCascVarLeastNbrClusters",&fTreeCascVarLeastNbrClusters,"fTreeCascVarLeastNbrClusters/I");
//-----------MULTIPLICITY-INFO--------------------
  fTreeCascade->Branch("fTreeCascVarCentV0A",&fTreeCascVarCentV0A,"fTreeCascVarCentV0A/F");
  fTreeCascade->Branch("fTreeCascVarCentV0C",&fTreeCascVarCentV0C,"fTreeCascVarCentV0C/F");
  fTreeCascade->Branch("fTreeCascVarCentV0M",&fTreeCascVarCentV0M,"fTreeCascVarCentV0M/F");
  fTreeCascade->Branch("fTreeCascVarCentV0AEq",&fTreeCascVarCentV0AEq,"fTreeCascVarCentV0AEq/F");
  fTreeCascade->Branch("fTreeCascVarCentV0CEq",&fTreeCascVarCentV0CEq,"fTreeCascVarCentV0CEq/F");
  fTreeCascade->Branch("fTreeCascVarCentV0MEq",&fTreeCascVarCentV0MEq,"fTreeCascVarCentV0MEq/F");
  fTreeCascade->Branch("fTreeCascVarAmpV0A",&fTreeCascVarAmpV0A,"fTreeCascVarAmpV0A/F");
  fTreeCascade->Branch("fTreeCascVarAmpV0C",&fTreeCascVarAmpV0C,"fTreeCascVarAmpV0C/F");
  fTreeCascade->Branch("fTreeCascVarAmpV0AEq",&fTreeCascVarAmpV0AEq,"fTreeCascVarAmpV0AEq/F");
  fTreeCascade->Branch("fTreeCascVarAmpV0CEq",&fTreeCascVarAmpV0CEq,"fTreeCascVarAmpV0CEq/F");
  fTreeCascade->Branch("fTreeCascVarRefMultEta8",&fTreeCascVarRefMultEta8,"fTreeCascVarRefMultEta8/I");
  fTreeCascade->Branch("fTreeCascVarRefMultEta5",&fTreeCascVarRefMultEta5,"fTreeCascVarRefMultEta5/I");
  fTreeCascade->Branch("fTreeCascVarRunNumber",&fTreeCascVarRunNumber,"fTreeCascVarRunNumber/I");
//-----------DECAY-LENGTH-INFO--------------------
  fTreeCascade->Branch("fTreeCascVarDistOverTotMom",&fTreeCascVarDistOverTotMom,"fTreeCascVarDistOverTotMom/F");
//------------------------------------------------
  fTreeCascade->Branch("fTreeCascVarNegNSigmaPion",&fTreeCascVarNegNSigmaPion,"fTreeCascVarNegNSigmaPion/F");
  fTreeCascade->Branch("fTreeCascVarNegNSigmaProton",&fTreeCascVarNegNSigmaProton,"fTreeCascVarNegNSigmaProton/F");
  fTreeCascade->Branch("fTreeCascVarPosNSigmaPion",&fTreeCascVarPosNSigmaPion,"fTreeCascVarPosNSigmaPion/F");
  fTreeCascade->Branch("fTreeCascVarPosNSigmaProton",&fTreeCascVarPosNSigmaProton,"fTreeCascVarPosNSigmaProton/F");
  fTreeCascade->Branch("fTreeCascVarBachNSigmaPion",&fTreeCascVarBachNSigmaPion,"fTreeCascVarBachNSigmaPion/F");
  fTreeCascade->Branch("fTreeCascVarBachNSigmaKaon",&fTreeCascVarBachNSigmaKaon,"fTreeCascVarBachNSigmaKaon/F");

//------------------------------------------------
// Particle Identification Setup
//------------------------------------------------
  
   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   fPIDResponse = inputHandler->GetPIDResponse();

  // Multiplicity
  if(! fESDtrackCuts ){
    fESDtrackCuts = new AliESDtrackCuts();
  }

//------------------------------------------------
// V0 Multiplicity Histograms
//------------------------------------------------

   // Create histograms
   OpenFile(1);
   fListHist = new TList();
   fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

   if(! fHistEventCounter ) {
    //Histogram Output: Event-by-Event
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",5,0,5); 
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
    fHistEventCounter->GetXaxis()->SetBinLabel(2, "Phys-Sel");  
    fHistEventCounter->GetXaxis()->SetBinLabel(3, "Has Vtx");  
    fHistEventCounter->GetXaxis()->SetBinLabel(4, "Vtx |z|<10cm");  
    fHistEventCounter->GetXaxis()->SetBinLabel(5, "Isn't Pileup");
    fListHist->Add(fHistEventCounter); 
   }

   //List of Histograms: Normal
   PostData(1, fListHist);

   //TTree Object: Saved to base directory. Should cache to disk while saving. 
   //(Important to avoid excessive memory usage, particularly when merging)
   PostData(2, fTreeEvent);
   PostData(3, fTreeV0);
   PostData(4, fTreeCascade);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicity::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

   AliESDEvent *lESDevent = 0x0;
  
  // Connect to the InputEvent   
  // After these lines, we should have an ESD/AOD event + the number of V0s in it.

   // Appropriate for ESD analysis! 
      
   lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
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
        
   fRunNumber = lESDevent->GetRunNumber();

  Double_t lMagneticField = -10; 
  lMagneticField = lESDevent->GetMagneticField( );

//------------------------------------------------
// Variable Definition
//------------------------------------------------

//------------------------------------------------
// Physics Selection
//------------------------------------------------
  
  fHistEventCounter->Fill(0.5); 

  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
  
  //Standard Min-Bias Selection
  if ( ! isSelected ) {
    PostData(1, fListHist);
    PostData(2, fTreeEvent);
    PostData(3, fTreeV0);
    PostData(4, fTreeCascade);
    return;
  }

  fHistEventCounter->Fill(1.5);
 
  //------------------------------------------------
  // Primary Vertex Requirements Section:
  //  ---> pp: has vertex, |z|<10cm
  //------------------------------------------------
  
  //classical Proton-proton like selection 
  const AliESDVertex *lPrimaryBestESDVtx     = lESDevent->GetPrimaryVertex();	
  const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
  const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();

  Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
  lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

  //Only accept if Tracking or SPD vertex is fine 
  if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtx->GetStatus() ){
    AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
    PostData(1, fListHist); 
    PostData(2, fTreeEvent);
    PostData(3, fTreeV0);
    PostData(4, fTreeCascade);
    return;
  }

  //Has SPD or Tracking Vertex
  fHistEventCounter -> Fill(2.5); 

  //Always do Primary Vertex Selection 
  if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0) {
    AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
    PostData(1, fListHist); 
    PostData(2, fTreeEvent);
    PostData(3, fTreeV0);
    PostData(4, fTreeCascade);
    return;
  }

  //Fill Event selected counter
  fHistEventCounter -> Fill(3.5);

  //------------------------------------------------
  // Check if this isn't pileup
  //------------------------------------------------

  if(lESDevent->IsPileupFromSPDInMultBins() ){
    // minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  
    //-> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
    AliWarning("Pb / Event tagged as pile-up by SPD... return !"); 
    PostData(1, fListHist); 
    PostData(2, fTreeEvent);
    PostData(3, fTreeV0);
    PostData(4, fTreeCascade);
    return; 
  }
  //Fill Event isn't pileup counter
  fHistEventCounter -> Fill(4.5);

//------------------------------------------------
// Multiplicity Information Acquistion
//------------------------------------------------

  //Standard GetReferenceMultiplicity Estimator (0.5 and 0.8)
  fRefMultEta5 = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,0.5);
  fRefMultEta8 = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,0.8);

  // VZERO PART
  Float_t  multV0A  = 0;            //  multiplicity from V0 reco side A
  Float_t  multV0C  = 0;            //  multiplicity from V0 reco side C
  Float_t  multV0AEq  = 0;          //  multiplicity from V0 reco side A
  Float_t  multV0CEq  = 0;          //  multiplicity from V0 reco side C
  Float_t  multV0ACorr  = 0;            //  multiplicity from V0 reco side A
  Float_t  multV0CCorr  = 0;            //  multiplicity from V0 reco side C

  //Non-Equalized Signal: copy of multV0ACorr and multV0CCorr from AliCentralitySelectionTask
  //Getters for uncorrected multiplicity  
  multV0A=esdV0->GetMTotV0A();
  multV0C=esdV0->GetMTotV0C();

  //Get Z vertex position of SPD vertex (why not Tracking if available?) 
  Float_t zvtx = lPrimarySPDVtx->GetZ(); 

  //Acquire Corrected multV0A 
  multV0ACorr = AliESDUtils::GetCorrV0A(multV0A,zvtx);    
  multV0CCorr = AliESDUtils::GetCorrV0C(multV0C,zvtx);   
    
  //Copy to Event Tree for extra information 
  fAmplitude_V0A = multV0ACorr; 
  fAmplitude_V0C = multV0CCorr; 

  // Equalized signals // From AliCentralitySelectionTask
  for(Int_t iCh = 4; iCh < 7; ++iCh) {
    Double_t mult = lESDevent->GetVZEROEqMultiplicity(iCh);
    multV0AEq += mult;
  }
  for(Int_t iCh = 0; iCh < 3; ++iCh) {
    Double_t mult = lESDevent->GetVZEROEqMultiplicity(iCh);
    multV0CEq += mult;
  }
  fAmplitude_V0AEq = multV0AEq; 
  fAmplitude_V0CEq = multV0CEq; 

  fCentrality_V0A   = -100; 
  fCentrality_V0C   = -100; 
  fCentrality_V0M   = -100; 
  fCentrality_V0AEq = -100; 
  fCentrality_V0CEq = -100; 
  fCentrality_V0MEq = -100; 

  //AliCentrality... Check if working? 
  AliCentrality* centrality;
  centrality = lESDevent->GetCentrality();
  if ( !(centrality->GetQuality()>1) ){ 
    fCentrality_V0A   = centrality->GetCentralityPercentile( "V0A"   ); 
    fCentrality_V0C   = centrality->GetCentralityPercentile( "V0C"   ); 
    fCentrality_V0M   = centrality->GetCentralityPercentile( "V0M"   ); 
    fCentrality_V0AEq = centrality->GetCentralityPercentile( "V0AEq" ); 
    fCentrality_V0CEq = centrality->GetCentralityPercentile( "V0CEq" ); 
    fCentrality_V0MEq = centrality->GetCentralityPercentile( "V0MEq" ); 
  }
  
  //Event-level fill 
  fTreeEvent->Fill() ;
  
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

   Int_t nv0s = 0;
   nv0s = lESDevent->GetNumberOfV0s();

   for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
   {// This is the begining of the V0 loop
      AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
      if (!v0) continue;

      Double_t tDecayVertexV0[3]; v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]); 

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

      Double_t lMomPos[3]; v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
      Double_t lMomNeg[3]; v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);

      AliESDtrack *pTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyPos);
      AliESDtrack *nTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyNeg);
      if (!pTrack || !nTrack) {
         Printf("ERROR: Could not retreive one of the daughter track");
         continue;
      }

      //Daughter Eta for Eta selection, afterwards
      fTreeVariableNegEta = nTrack->Eta();
      fTreeVariablePosEta = pTrack->Eta();

      // Filter like-sign V0 (next: add counter and distribution)
      if ( pTrack->GetSign() == nTrack->GetSign()){
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
     
  
      if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;
	
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
      if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;

      //End track Quality Cuts
      //________________________________________________________________________

      lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lBestPrimaryVtxPos[0],
							lBestPrimaryVtxPos[1],
							lMagneticField) );

      lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lBestPrimaryVtxPos[0],
							lBestPrimaryVtxPos[1],
							lMagneticField) );

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
      fTreeVariableNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
      fTreeVariableNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
      fTreeVariableNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
      fTreeVariableNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );

//This requires an Invariant Mass Hypothesis afterwards
      fTreeVariableDistOverTotMom = TMath::Sqrt(
						TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
					);
      fTreeVariableDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure

      //Copy Multiplicity information 
      fTreeVariableCentV0A = fCentrality_V0A; 
      fTreeVariableCentV0C = fCentrality_V0C; 
      fTreeVariableCentV0M = fCentrality_V0M; 
      fTreeVariableCentV0AEq = fCentrality_V0AEq; 
      fTreeVariableCentV0CEq = fCentrality_V0CEq; 
      fTreeVariableCentV0MEq = fCentrality_V0MEq; 
      fTreeVariableAmpV0A = fAmplitude_V0A; 
      fTreeVariableAmpV0C = fAmplitude_V0C; 
      fTreeVariableAmpV0AEq = fAmplitude_V0AEq; 
      fTreeVariableAmpV0CEq = fAmplitude_V0CEq; 
      fTreeVariableRefMultEta8 = fRefMultEta8;
      fTreeVariableRefMultEta5 = fRefMultEta5;
      fTreeVariableRunNumber = fRunNumber; 

//------------------------------------------------
// Fill Tree! 
//------------------------------------------------
     
     // The conditionals are meant to decrease excessive
     // memory usage!
     
     //First Selection: Reject OnFly
     if( lOnFlyStatus == 0 ){
       //Second Selection: rough 20-sigma band, parametric.
       //K0Short: Enough to parametrize peak broadening with linear function.
       Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeVariablePt;
       Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeVariablePt;
       //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
       //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
       Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeVariablePt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeVariablePt);
       Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeVariablePt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeVariablePt);
       //Do Selection
       if( (fTreeVariableInvMassLambda    < lUpperLimitLambda  && fTreeVariableInvMassLambda     > lLowerLimitLambda     ) ||
          (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda     ) ||
          (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short    ) ){
         //Pre-selection in case this is AA...
         if ( TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8 && fkSaveV0Tree ) fTreeV0->Fill();
       }
     }
   }// This is the end of the V0 loop

//------------------------------------------------
// Fill V0 tree over.
//------------------------------------------------



//------------------------------------------------
// Rerun cascade vertexer! 
//------------------------------------------------
    
  if( fkRunVertexers ){ 
    lESDevent->ResetCascades();
    lESDevent->ResetV0s();

    AliV0vertexer lV0vtxer;
    AliCascadeVertexer lCascVtxer;
                  
    lV0vtxer.SetDefaultCuts(fV0VertexerSels);
    lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);

    lV0vtxer.Tracks2V0vertices(lESDevent);
    lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
  }

//------------------------------------------------
// MAIN CASCADE LOOP STARTS HERE
//------------------------------------------------
// Code Credit: Antonin Maire (thanks^100)
// ---> This is an adaptation

  Long_t ncascades = 0;
	ncascades = lESDevent->GetNumberOfCascades();
  
  for (Int_t iXi = 0; iXi < ncascades; iXi++){
    //------------------------------------------------
    // Initializations
    //------------------------------------------------	
	  //Double_t lTrkgPrimaryVtxRadius3D = -500.0;
	  //Double_t lBestPrimaryVtxRadius3D = -500.0;

	  // - 1st part of initialisation : variables needed to store AliESDCascade data members
	  Double_t lEffMassXi      = 0. ;
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
		
	  Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.;
	  Double_t lDcaPosToPrimVertexXi  = -1.;
	  Double_t lDcaNegToPrimVertexXi  = -1.;
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
    
	  // - 6th part of initialisation : extra info for QA
	  Double_t lXiMomX       = 0. , lXiMomY = 0., lXiMomZ = 0.;
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
	
	  Short_t  lChargeXi = -2;
	  //Double_t lV0toXiCosineOfPointingAngle = 0. ;
	
	  Double_t lRapXi   = -20.0, lRapOmega = -20.0; //  lEta = -20.0, lTheta = 360., lPhi = 720. ;
	  //Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
	    
    // -------------------------------------
    // II.ESD - Calculation Part dedicated to Xi vertices (ESD)
    
	  AliESDcascade *xi = lESDevent->GetCascade(iXi);
	  if (!xi) continue;
	
		// - II.Step 2 : Assigning the necessary variables for specific AliESDcascade data members (ESD)	
		//-------------
	  lV0quality = 0.;
	  xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay

	  lEffMassXi  			= xi->GetEffMassXi();
	  //lChi2Xi 			    = xi->GetChi2Xi();
	  lDcaXiDaughters 	= xi->GetDcaXiDaughters();
	  lXiCosineOfPointingAngle 	            = xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                                 lBestPrimaryVtxPos[1],
                                                                                 lBestPrimaryVtxPos[2] );
		  // Take care : the best available vertex should be used (like in AliCascadeVertexer)
	
	  xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] ); 
	  lXiRadius			= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );		

		// - II.Step 3 : around the tracks : Bach + V0 (ESD)
		// ~ Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
		//-------------
		
        UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
        UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
        UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
                // Care track label can be negative in MC production (linked with the track quality)
                // However = normally, not the case for track index ...
          
	  // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
	  if(lBachIdx == lIdxNegXi) {
		  AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
	  }
    if(lBachIdx == lIdxPosXi) {
    	AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
	  }
          
	  AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
	  AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
	  AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx );

	  if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
		  AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
		  continue;
	  }

      fTreeCascVarPosEta = pTrackXi->Eta();
      fTreeCascVarNegEta = nTrackXi->Eta();
      fTreeCascVarBachEta = bachTrackXi->Eta();
      
      Double_t lBMom[3], lNMom[3], lPMom[3];
      xi->GetBPxPyPz( lBMom[0], lBMom[1], lBMom[2] );
      xi->GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
      xi->GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );
      
      //fTreeCascVarBachTransMom = TMath::Sqrt( lBMom[0]*lBMom[0] + lBMom[1]*lBMom[1] );
      //fTreeCascVarPosTransMom  = TMath::Sqrt( lPMom[0]*lPMom[0] + lPMom[1]*lPMom[1] );
      //fTreeCascVarNegTransMom  = TMath::Sqrt( lNMom[0]*lNMom[0] + lNMom[1]*lNMom[1] );
  
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

    if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
    if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
    if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }

	  // 2 - Poor quality related to TPC clusters: lowest cut of 70 clusters
    if(lPosTPCClusters  < 70) { AliWarning("Pb / V0 Pos. track has less than 70 TPC clusters ... continue!"); continue; }
	  if(lNegTPCClusters  < 70) { AliWarning("Pb / V0 Neg. track has less than 70 TPC clusters ... continue!"); continue; }
	  if(lBachTPCClusters < 70) { AliWarning("Pb / Bach.   track has less than 70 TPC clusters ... continue!"); continue; }
	  Int_t leastnumberofclusters = 1000; 
	  if( lPosTPCClusters < leastnumberofclusters ) leastnumberofclusters = lPosTPCClusters;
	  if( lNegTPCClusters < leastnumberofclusters ) leastnumberofclusters = lNegTPCClusters;
	  if( lBachTPCClusters < leastnumberofclusters ) leastnumberofclusters = lBachTPCClusters;

	  lInvMassLambdaAsCascDghter	= xi->GetEffMass();
	  // This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
	  lDcaV0DaughtersXi 		= xi->GetDcaV0Daughters(); 
	  //lV0Chi2Xi 			= xi->GetChi2V0();
	
	  lV0CosineOfPointingAngleXi 	= xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
									    lBestPrimaryVtxPos[1],
									    lBestPrimaryVtxPos[2] );
    //Modification: V0 CosPA wrt to Cascade decay vertex
	  lV0CosineOfPointingAngleXiSpecial 	= xi->GetV0CosineOfPointingAngle( lPosXi[0],
									    lPosXi[1],
									    lPosXi[2] );

	  lDcaV0ToPrimVertexXi 		= xi->GetD( lBestPrimaryVtxPos[0], 
						      lBestPrimaryVtxPos[1], 
						      lBestPrimaryVtxPos[2] );
		
	  lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD(	lBestPrimaryVtxPos[0], 
						       		lBestPrimaryVtxPos[1], 
						       		lMagneticField  ) ); 
					  // Note : AliExternalTrackParam::GetD returns an algebraic value ...
		
	  xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] ); 
	  lV0RadiusXi		= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
	
	  lDcaPosToPrimVertexXi 	= TMath::Abs( pTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
						   		lBestPrimaryVtxPos[1], 
						   		lMagneticField  )     ); 
	
	  lDcaNegToPrimVertexXi 	= TMath::Abs( nTrackXi	->GetD(	lBestPrimaryVtxPos[0], 
					      			lBestPrimaryVtxPos[1], 
					      			lMagneticField  )     ); 
		
	  // - II.Step 4 : around effective masses (ESD)
	  // ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		
	  if( bachTrackXi->Charge() < 0 )	{
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
		  xi->ChangeMassHypothesis(lV0quality , 3312); 	// Back to default hyp.
	  }// end if negative bachelor
	
	
	  if( bachTrackXi->Charge() >  0 ){
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
		  xi->ChangeMassHypothesis(lV0quality , -3312); 	// Back to "default" hyp.
	  }// end if positive bachelor
		  // - II.Step 6 : extra info for QA (ESD)
		  // miscellaneous pieces of info that may help regarding data quality assessment.
		  //-------------

	  xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
		  lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
		  lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
		
	  xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
		  //lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
		  //lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );

	  lChargeXi = xi->Charge();

	  //lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
	
	  lRapXi    = xi->RapXi();
	  lRapOmega = xi->RapOmega();
	  //lEta      = xi->Eta();
	  //lTheta    = xi->Theta() *180.0/TMath::Pi();
	  //lPhi      = xi->Phi()   *180.0/TMath::Pi();
	  //lAlphaXi  = xi->AlphaXi();
	  //lPtArmXi  = xi->PtArmXi();

  //------------------------------------------------
  // Set Variables for adding to tree
  //------------------------------------------------		
	
          fTreeCascVarCharge	= lChargeXi;
          if(lInvMassXiMinus!=0)    fTreeCascVarMassAsXi = lInvMassXiMinus;
          if(lInvMassXiPlus!=0)     fTreeCascVarMassAsXi = lInvMassXiPlus;
          if(lInvMassOmegaMinus!=0) fTreeCascVarMassAsOmega = lInvMassOmegaMinus;
          if(lInvMassOmegaPlus!=0)  fTreeCascVarMassAsOmega = lInvMassOmegaPlus;
          fTreeCascVarPt = lXiTransvMom;
          fTreeCascVarRapXi = lRapXi ;
          fTreeCascVarRapOmega = lRapOmega ;
          fTreeCascVarDCACascDaughters = lDcaXiDaughters;
          fTreeCascVarDCABachToPrimVtx = lDcaBachToPrimVertexXi;
          fTreeCascVarDCAV0Daughters = lDcaV0DaughtersXi;
          fTreeCascVarDCAV0ToPrimVtx = lDcaV0ToPrimVertexXi;
          fTreeCascVarDCAPosToPrimVtx = lDcaPosToPrimVertexXi;
          fTreeCascVarDCANegToPrimVtx = lDcaNegToPrimVertexXi;
          fTreeCascVarCascCosPointingAngle = lXiCosineOfPointingAngle;
          fTreeCascVarCascRadius = lXiRadius;
          fTreeCascVarV0Mass = lInvMassLambdaAsCascDghter;
          fTreeCascVarV0CosPointingAngle = lV0CosineOfPointingAngleXi;
          fTreeCascVarV0CosPointingAngleSpecial = lV0CosineOfPointingAngleXiSpecial;
          fTreeCascVarV0Radius = lV0RadiusXi;
          fTreeCascVarLeastNbrClusters = leastnumberofclusters;

          //Copy Multiplicity information 
          fTreeCascVarCentV0A = fCentrality_V0A; 
          fTreeCascVarCentV0C = fCentrality_V0C; 
          fTreeCascVarCentV0M = fCentrality_V0M; 
          fTreeCascVarCentV0AEq = fCentrality_V0AEq; 
          fTreeCascVarCentV0CEq = fCentrality_V0CEq; 
          fTreeCascVarCentV0MEq = fCentrality_V0MEq; 
          fTreeCascVarAmpV0A = fAmplitude_V0A; 
          fTreeCascVarAmpV0C = fAmplitude_V0C; 
          fTreeCascVarAmpV0AEq = fAmplitude_V0AEq; 
          fTreeCascVarAmpV0CEq = fAmplitude_V0CEq; 
          fTreeCascVarRefMultEta8 = fRefMultEta8;
          fTreeCascVarRefMultEta5 = fRefMultEta5;
          fTreeCascVarRunNumber = fRunNumber; 

          fTreeCascVarDistOverTotMom = TMath::Sqrt(
						TMath::Power( lPosXi[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( lPosXi[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( lPosXi[2] - lBestPrimaryVtxPos[2] , 2)
					);
          fTreeCascVarDistOverTotMom /= (lXiTotMom+1e-13);

//All vars not specified here: specified elsewhere!

//------------------------------------------------
// Fill Tree! 
//------------------------------------------------

// The conditional is meant to decrease excessive
// memory usage! Be careful when loosening the 
// cut!

  //Xi    Mass window: 150MeV wide
  //Omega mass window: 150MeV wide

  if( fkSaveCascadeTree && ( (fTreeCascVarMassAsXi<1.32+0.075&&fTreeCascVarMassAsXi>1.32-0.075) ||
      (fTreeCascVarMassAsOmega<1.68+0.075&&fTreeCascVarMassAsOmega>1.68-0.075) ) ){
      fTreeCascade->Fill();
  }

//------------------------------------------------
// Fill tree over.
//------------------------------------------------

	}// end of the Cascade loop (ESD or AOD)

  // Post output data.
  PostData(1, fListHist); 
  PostData(2, fTreeEvent);
  PostData(3, fTreeV0);
  PostData(4, fTreeCascade);
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicity::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   TList *cRetrievedList = 0x0;
   cRetrievedList = (TList*)GetOutputData(1);
   if(!cRetrievedList){
      Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicity : ouput data container list not available\n");
      return;
   }	
	
   fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
   if (!fHistEventCounter) {
      Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicity : fHistEventCounter not available");
      return;
   }
  
   TCanvas *canCheck = new TCanvas("AliAnalysisTaskStrangenessVsMultiplicity","V0 Multiplicity",10,10,510,510);
   canCheck->cd(1)->SetLogy();

   fHistEventCounter->SetMarkerStyle(22);
   fHistEventCounter->DrawCopy("E");
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskStrangenessVsMultiplicity::MyRapidity(Double_t rE, Double_t rPz) const
{
   // Local calculation for rapidity
   Double_t ReturnValue = -100;
   if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ){ 
      ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
   }
   return ReturnValue;
} 
