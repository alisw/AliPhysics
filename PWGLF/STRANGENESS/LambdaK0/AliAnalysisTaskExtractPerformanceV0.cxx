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
// Modified version of AliAnalysisTaskCheckCascade.cxx.
// This is a 'hybrid' output version, in that it uses a classic TTree
// ROOT object to store the candidates, plus a couple of histograms filled on
// a per-event basis for storing variables too numerous to put in a tree. 
//
// --- Adapted to look for lambdas as well, using code from 
//        AliAnalysisTaskCheckPerformanceStrange.cxx
//
//  --- Algorithm Description 
//   1. Loop over primaries in stack to acquire generated charged Xi
//   2. Loop over stack to find V0s, fill TH3Fs "PrimRawPt"s for Efficiency
//   3. Perform Physics Selection
//   4. Perform Primary Vertex |z|<10cm selection
//   5. Perform Primary Vertex NoTPCOnly vertexing selection (>0 contrib.)
//   6. Perform Pileup Rejection
//   7. Analysis Loops: 
//    7a. Fill TH3Fs "PrimAnalysisPt" for control purposes only
//    7b. Fill TTree object with V0 information, candidates
//
//  Please Report Any Bugs! 
//
//   --- David Dobrigkeit Chinellato
//        (david.chinellato@gmail.com)
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
#include "AliLog.h"

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

#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"

#include "AliAnalysisTaskExtractPerformanceV0.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskExtractPerformanceV0)

AliAnalysisTaskExtractPerformanceV0::AliAnalysisTaskExtractPerformanceV0() 
  : AliAnalysisTaskSE(), fListHistV0(0), fTree(0), fPIDResponse(0), fESDtrackCuts(0),
   fkIsNuclear   ( kFALSE ), 
   fkSwitchINT7  ( kFALSE ),
   fkUseOnTheFly ( kFALSE ),
   fkTakeAllTracks ( kFALSE ),
   fpArapidityShift ( 0.465 ),
  fCentralityEstimator("V0M"),
  fkLightWeight  ( kFALSE ),
  fkFastOnly     ( "" ),
//------------------------------------------------
// Tree Variables 

  fTreeVariablePrimaryStatus(0),
  fTreeVariablePrimaryStatusMother(0),
  fTreeVariableChi2V0(0),
	fTreeVariableDcaV0Daughters(0),
	fTreeVariableDcaV0ToPrimVertex(0),
	fTreeVariableDcaPosToPrimVertex(0),
	fTreeVariableDcaNegToPrimVertex(0),
	fTreeVariableV0CosineOfPointingAngle(0),
	fTreeVariableV0Radius(0),
	fTreeVariablePt(0),
	fTreeVariablePtMC(0),
	fTreeVariableRapK0Short(0),
	fTreeVariableRapLambda(0),
	fTreeVariableRapMC(0),
	fTreeVariableInvMassK0s(0),
	fTreeVariableInvMassLambda(0),
	fTreeVariableInvMassAntiLambda(0),
	fTreeVariableAlphaV0(0),
	fTreeVariablePtArmV0(0),
	fTreeVariableNegTotMomentum(0),
	fTreeVariablePosTotMomentum(0),
	fTreeVariableNegTransvMomentum(0),
	fTreeVariablePosTransvMomentum(0),
	fTreeVariableNegTransvMomentumMC(0),
	fTreeVariablePosTransvMomentumMC(0),
   
	fTreeVariableNSigmasPosProton(0),
	fTreeVariableNSigmasPosPion(0),
	fTreeVariableNSigmasNegProton(0),
	fTreeVariableNSigmasNegPion(0),

	fTreeVariablePtMother(0),
	fTreeVariableV0CreationRadius(0),
  fTreeVariablePID(0),
  fTreeVariablePIDPositive(0),
  fTreeVariablePIDNegative(0),
  fTreeVariablePIDMother(0),
  fTreeVariableIndexStatus(0),
  fTreeVariableIndexStatusMother(0),

  fTreeVariableRunNumber(0),
  fTreeVariableEventNumber(0),

	fTreeVariableDistOverTotMom(0),

	fTreeVariablePosEta(0),
	fTreeVariableNegEta(0),

	fTreeVariableVertexZ(0),

  fTreeVariableLeastNbrCrossedRows(0),
  fTreeVariableLeastRatioCrossedRowsOverFindable(0),

  fTreeVariableMultiplicity(0),
  fTreeVariableMultiplicityV0A(0),
  fTreeVariableMultiplicityZNA(0),
  fTreeVariableMultiplicityTRK(0),
  fTreeVariableMultiplicitySPD(0),
  fTreeVariableMultiplicityMC(0),

  fTreeVariableV0x(0),
  fTreeVariableV0y(0),
  fTreeVariableV0z(0),

  fTreeVariableV0Px(0),
  fTreeVariableV0Py(0),
  fTreeVariableV0Pz(0),

  fTreeVariableMCV0x(0),
  fTreeVariableMCV0y(0),
  fTreeVariableMCV0z(0),

  fTreeVariableMCV0Px(0),
  fTreeVariableMCV0Py(0),
  fTreeVariableMCV0Pz(0),

  fTreeVariablePVx(0),
  fTreeVariablePVy(0),
  fTreeVariablePVz(0),

  fTreeVariableMCPVx(0),
  fTreeVariableMCPVy(0),
  fTreeVariableMCPVz(0),

  fTreeVariableIsNonInjected(0),

  fTreeVariableNegTrackStatus(0),
  fTreeVariablePosTrackStatus(0),

//------------------------------------------------
// HISTOGRAMS
// --- Filled on an Event-by-event basis
//------------------------------------------------
   fHistV0MultiplicityBeforeTrigSel(0),
   fHistV0MultiplicityForTrigEvt(0), 
   fHistV0MultiplicityForSelEvt(0),
   fHistV0MultiplicityForSelEvtNoTPCOnly(0),
   fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup(0),
   fHistMultiplicityBeforeTrigSel(0),
   fHistMultiplicityForTrigEvt(0),
   fHistMultiplicity(0),
   fHistMultiplicityNoTPCOnly(0),
   fHistMultiplicityNoTPCOnlyNoPileup(0),

//V0A Centrality
fHistMultiplicityV0ABeforeTrigSel(0),
fHistMultiplicityV0AForTrigEvt(0),
fHistMultiplicityV0A(0),
fHistMultiplicityV0ANoTPCOnly(0),
fHistMultiplicityV0ANoTPCOnlyNoPileup(0),

//ZNA Centrality
fHistMultiplicityZNABeforeTrigSel(0),
fHistMultiplicityZNAForTrigEvt(0),
fHistMultiplicityZNA(0),
fHistMultiplicityZNANoTPCOnly(0),
fHistMultiplicityZNANoTPCOnlyNoPileup(0),

//TRK Centrality
fHistMultiplicityTRKBeforeTrigSel(0),
fHistMultiplicityTRKForTrigEvt(0),
fHistMultiplicityTRK(0),
fHistMultiplicityTRKNoTPCOnly(0),
fHistMultiplicityTRKNoTPCOnlyNoPileup(0),

//SPD Centrality
fHistMultiplicitySPDBeforeTrigSel(0),
fHistMultiplicitySPDForTrigEvt(0),
fHistMultiplicitySPD(0),
fHistMultiplicitySPDNoTPCOnly(0),
fHistMultiplicitySPDNoTPCOnlyNoPileup(0),

	f2dHistMultiplicityVsTrueBeforeTrigSel(0),
	f2dHistMultiplicityVsTrueForTrigEvt(0),
	f2dHistMultiplicityVsTrue(0),
	f2dHistMultiplicityVsTrueNoTPCOnly(0),
	f2dHistMultiplicityVsTrueNoTPCOnlyNoPileup(0),

  //Raw Data for Vertex Z position estimator change
	f2dHistMultiplicityVsVertexZBeforeTrigSel(0),
	f2dHistMultiplicityVsVertexZForTrigEvt(0),
	f2dHistMultiplicityVsVertexZ(0),
	f2dHistMultiplicityVsVertexZNoTPCOnly(0),
	f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup(0),

  fHistGenVertexZBeforeTrigSel(0),     
  fHistGenVertexZForTrigEvt(0),
  fHistGenVertexZ(0),
  fHistGenVertexZNoTPCOnly(0),
  fHistGenVertexZNoTPCOnlyNoPileup(0),

//------------------------------------------------
// PARTICLE HISTOGRAMS
// --- Filled on a Particle-by-Particle basis
//------------------------------------------------
   f3dHistPrimAnalysisPtVsYVsMultLambda(0),
   f3dHistPrimAnalysisPtVsYVsMultAntiLambda(0),
   f3dHistPrimAnalysisPtVsYVsMultK0Short(0),
   f3dHistPrimAnalysisPtVsYCMSVsMultLambda(0),
   f3dHistPrimAnalysisPtVsYCMSVsMultAntiLambda(0),
   f3dHistPrimAnalysisPtVsYCMSVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsMultLambda(0),
   f3dHistPrimRawPtVsYVsMultAntiLambda(0),
   f3dHistPrimRawPtVsYVsMultK0Short(0),
   f3dHistPrimRawPtVsYCMSVsMultLambda(0),
   f3dHistPrimRawPtVsYCMSVsMultAntiLambda(0),
   f3dHistPrimRawPtVsYCMSVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsMultNonInjLambda(0),
   f3dHistPrimRawPtVsYVsMultNonInjAntiLambda(0),
   f3dHistPrimRawPtVsYVsMultNonInjK0Short(0),
   f3dHistPrimRawPtVsYVsMultMCLambda(0),
   f3dHistPrimRawPtVsYVsMultMCAntiLambda(0),
   f3dHistPrimRawPtVsYVsMultMCK0Short(0),
   f3dHistPrimRawPtVsYVsVertexZLambda(0),
   f3dHistPrimRawPtVsYVsVertexZAntiLambda(0),
   f3dHistPrimRawPtVsYVsVertexZK0Short(0),
   f3dHistPrimCloseToPVPtVsYVsMultLambda(0),
   f3dHistPrimCloseToPVPtVsYVsMultAntiLambda(0),
   f3dHistPrimCloseToPVPtVsYVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsDecayLengthLambda(0),
   f3dHistPrimRawPtVsYVsDecayLengthAntiLambda(0),
   f3dHistPrimRawPtVsYVsDecayLengthK0Short(0),
   f3dHistGenPtVsYVsMultXiMinus(0),
   f3dHistGenPtVsYVsMultXiPlus(0),
   f3dHistGenPtVsYVsMultOmegaMinus(0),
   f3dHistGenPtVsYVsMultOmegaPlus(0),
   f3dHistGenSelectedPtVsYVsMultXiMinus(0),
   f3dHistGenSelectedPtVsYVsMultXiPlus(0),
   f3dHistGenSelectedPtVsYVsMultOmegaMinus(0),
   f3dHistGenSelectedPtVsYVsMultOmegaPlus(0),
   f3dHistGenPtVsYCMSVsMultXiMinus(0),
   f3dHistGenPtVsYCMSVsMultXiPlus(0),
   f3dHistGenPtVsYCMSVsMultOmegaMinus(0),
   f3dHistGenPtVsYCMSVsMultOmegaPlus(0),
   f3dHistGenSelectedPtVsYCMSVsMultXiMinus(0),
   f3dHistGenSelectedPtVsYCMSVsMultXiPlus(0),
   f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus(0),
   f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus(0),
   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0),
   fHistPVxAnalysisHasHighPtLambda(0),
   fHistPVyAnalysisHasHighPtLambda(0),
   fHistPVzAnalysisHasHighPtLambda(0),
   fHistSwappedV0Counter(0)
{
  // Dummy Constructor
}

AliAnalysisTaskExtractPerformanceV0::AliAnalysisTaskExtractPerformanceV0(const char *name) 
  : AliAnalysisTaskSE(name), fListHistV0(0), fTree(0), fPIDResponse(0), fESDtrackCuts(0),
   fkIsNuclear   ( kFALSE ), 
   fkSwitchINT7  ( kFALSE ),
   fkUseOnTheFly ( kFALSE ),
   fkTakeAllTracks ( kFALSE ),
   fpArapidityShift ( 0.465 ),
  fCentralityEstimator("V0M"),
  fkLightWeight  ( kFALSE ),
  fkFastOnly     ( "" ),
//------------------------------------------------
// Tree Variables 

  fTreeVariablePrimaryStatus(0),
  fTreeVariablePrimaryStatusMother(0),
  fTreeVariableChi2V0(0),
	fTreeVariableDcaV0Daughters(0),
	fTreeVariableDcaV0ToPrimVertex(0),
	fTreeVariableDcaPosToPrimVertex(0),
	fTreeVariableDcaNegToPrimVertex(0),
	fTreeVariableV0CosineOfPointingAngle(0),
	fTreeVariableV0Radius(0),
	fTreeVariablePt(0),
	fTreeVariablePtMC(0),
	fTreeVariableRapK0Short(0),
	fTreeVariableRapLambda(0),
	fTreeVariableRapMC(0),
	fTreeVariableInvMassK0s(0),
	fTreeVariableInvMassLambda(0),
	fTreeVariableInvMassAntiLambda(0),
	fTreeVariableAlphaV0(0),
	fTreeVariablePtArmV0(0),
	fTreeVariableNegTotMomentum(0),
	fTreeVariablePosTotMomentum(0),
	fTreeVariableNegTransvMomentum(0),
	fTreeVariablePosTransvMomentum(0),
	fTreeVariableNegTransvMomentumMC(0),
	fTreeVariablePosTransvMomentumMC(0),
   
	fTreeVariableNSigmasPosProton(0),
	fTreeVariableNSigmasPosPion(0),
	fTreeVariableNSigmasNegProton(0),
	fTreeVariableNSigmasNegPion(0),

	fTreeVariablePtMother(0),
	fTreeVariableV0CreationRadius(0),
  fTreeVariablePID(0),
  fTreeVariablePIDPositive(0),
  fTreeVariablePIDNegative(0),
  fTreeVariablePIDMother(0),
  fTreeVariableIndexStatus(0),
  fTreeVariableIndexStatusMother(0),

  fTreeVariableRunNumber(0),
  fTreeVariableEventNumber(0),

	fTreeVariableDistOverTotMom(0),

	fTreeVariablePosEta(0),
	fTreeVariableNegEta(0),

	fTreeVariableVertexZ(0),

  fTreeVariableLeastNbrCrossedRows(0),
  fTreeVariableLeastRatioCrossedRowsOverFindable(0),
  fTreeVariableMultiplicity(0),
  fTreeVariableMultiplicityV0A(0),
  fTreeVariableMultiplicityZNA(0),
  fTreeVariableMultiplicityTRK(0),
  fTreeVariableMultiplicitySPD(0),
  fTreeVariableMultiplicityMC(0),

  fTreeVariableV0x(0),
  fTreeVariableV0y(0),
  fTreeVariableV0z(0),

  fTreeVariableV0Px(0),
  fTreeVariableV0Py(0),
  fTreeVariableV0Pz(0),

  fTreeVariableMCV0x(0),
  fTreeVariableMCV0y(0),
  fTreeVariableMCV0z(0),

  fTreeVariableMCV0Px(0),
  fTreeVariableMCV0Py(0),
  fTreeVariableMCV0Pz(0),

  fTreeVariablePVx(0),
  fTreeVariablePVy(0),
  fTreeVariablePVz(0),

  fTreeVariableMCPVx(0),
  fTreeVariableMCPVy(0),
  fTreeVariableMCPVz(0),

  fTreeVariableIsNonInjected(0),

  fTreeVariableNegTrackStatus(0),
  fTreeVariablePosTrackStatus(0),

//------------------------------------------------
// HISTOGRAMS
// --- Filled on an Event-by-event basis
//------------------------------------------------
   fHistV0MultiplicityBeforeTrigSel(0),
   fHistV0MultiplicityForTrigEvt(0), 
   fHistV0MultiplicityForSelEvt(0),
   fHistV0MultiplicityForSelEvtNoTPCOnly(0),
   fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup(0),
   fHistMultiplicityBeforeTrigSel(0),
   fHistMultiplicityForTrigEvt(0),
   fHistMultiplicity(0),
   fHistMultiplicityNoTPCOnly(0),
   fHistMultiplicityNoTPCOnlyNoPileup(0),

//V0A Centrality
fHistMultiplicityV0ABeforeTrigSel(0),
fHistMultiplicityV0AForTrigEvt(0),
fHistMultiplicityV0A(0),
fHistMultiplicityV0ANoTPCOnly(0),
fHistMultiplicityV0ANoTPCOnlyNoPileup(0),

//ZNA Centrality
fHistMultiplicityZNABeforeTrigSel(0),
fHistMultiplicityZNAForTrigEvt(0),
fHistMultiplicityZNA(0),
fHistMultiplicityZNANoTPCOnly(0),
fHistMultiplicityZNANoTPCOnlyNoPileup(0),

//TRK Centrality
fHistMultiplicityTRKBeforeTrigSel(0),
fHistMultiplicityTRKForTrigEvt(0),
fHistMultiplicityTRK(0),
fHistMultiplicityTRKNoTPCOnly(0),
fHistMultiplicityTRKNoTPCOnlyNoPileup(0),

//SPD Centrality
fHistMultiplicitySPDBeforeTrigSel(0),
fHistMultiplicitySPDForTrigEvt(0),
fHistMultiplicitySPD(0),
fHistMultiplicitySPDNoTPCOnly(0),
fHistMultiplicitySPDNoTPCOnlyNoPileup(0),

	f2dHistMultiplicityVsTrueBeforeTrigSel(0),
	f2dHistMultiplicityVsTrueForTrigEvt(0),
	f2dHistMultiplicityVsTrue(0),
	f2dHistMultiplicityVsTrueNoTPCOnly(0),
	f2dHistMultiplicityVsTrueNoTPCOnlyNoPileup(0),

  //Raw Data for Vertex Z position estimator change
	f2dHistMultiplicityVsVertexZBeforeTrigSel(0),
	f2dHistMultiplicityVsVertexZForTrigEvt(0),
	f2dHistMultiplicityVsVertexZ(0),
	f2dHistMultiplicityVsVertexZNoTPCOnly(0),
	f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup(0),

  fHistGenVertexZBeforeTrigSel(0),     
  fHistGenVertexZForTrigEvt(0),
  fHistGenVertexZ(0),
  fHistGenVertexZNoTPCOnly(0),
  fHistGenVertexZNoTPCOnlyNoPileup(0),

//------------------------------------------------
// PARTICLE HISTOGRAMS
// --- Filled on a Particle-by-Particle basis
//------------------------------------------------
   f3dHistPrimAnalysisPtVsYVsMultLambda(0),
   f3dHistPrimAnalysisPtVsYVsMultAntiLambda(0),
   f3dHistPrimAnalysisPtVsYVsMultK0Short(0),
   f3dHistPrimAnalysisPtVsYCMSVsMultLambda(0),
   f3dHistPrimAnalysisPtVsYCMSVsMultAntiLambda(0),
   f3dHistPrimAnalysisPtVsYCMSVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsMultLambda(0),
   f3dHistPrimRawPtVsYVsMultAntiLambda(0),
   f3dHistPrimRawPtVsYVsMultK0Short(0),
   f3dHistPrimRawPtVsYCMSVsMultLambda(0),
   f3dHistPrimRawPtVsYCMSVsMultAntiLambda(0),
   f3dHistPrimRawPtVsYCMSVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsMultNonInjLambda(0),
   f3dHistPrimRawPtVsYVsMultNonInjAntiLambda(0),
   f3dHistPrimRawPtVsYVsMultNonInjK0Short(0),
   f3dHistPrimRawPtVsYVsMultMCLambda(0),
   f3dHistPrimRawPtVsYVsMultMCAntiLambda(0),
   f3dHistPrimRawPtVsYVsMultMCK0Short(0),
   f3dHistPrimRawPtVsYVsVertexZLambda(0),
   f3dHistPrimRawPtVsYVsVertexZAntiLambda(0),
   f3dHistPrimRawPtVsYVsVertexZK0Short(0),
   f3dHistPrimCloseToPVPtVsYVsMultLambda(0),
   f3dHistPrimCloseToPVPtVsYVsMultAntiLambda(0),
   f3dHistPrimCloseToPVPtVsYVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsDecayLengthLambda(0),
   f3dHistPrimRawPtVsYVsDecayLengthAntiLambda(0),
   f3dHistPrimRawPtVsYVsDecayLengthK0Short(0),
   f3dHistGenPtVsYVsMultXiMinus(0),
   f3dHistGenPtVsYVsMultXiPlus(0),
   f3dHistGenPtVsYVsMultOmegaMinus(0),
   f3dHistGenPtVsYVsMultOmegaPlus(0),
   f3dHistGenSelectedPtVsYVsMultXiMinus(0),
   f3dHistGenSelectedPtVsYVsMultXiPlus(0),
   f3dHistGenSelectedPtVsYVsMultOmegaMinus(0),
   f3dHistGenSelectedPtVsYVsMultOmegaPlus(0),
   f3dHistGenPtVsYCMSVsMultXiMinus(0),
   f3dHistGenPtVsYCMSVsMultXiPlus(0),
   f3dHistGenPtVsYCMSVsMultOmegaMinus(0),
   f3dHistGenPtVsYCMSVsMultOmegaPlus(0),
   f3dHistGenSelectedPtVsYCMSVsMultXiMinus(0),
   f3dHistGenSelectedPtVsYCMSVsMultXiPlus(0),
   f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus(0),
   f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus(0),
   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0),
   fHistPVxAnalysisHasHighPtLambda(0),
   fHistPVyAnalysisHasHighPtLambda(0),
   fHistPVzAnalysisHasHighPtLambda(0),
   fHistSwappedV0Counter(0)
{
   // Constructor
   // Output slot #0 writes into a TList container (Cascade)
   DefineOutput(1, TList::Class());
   DefineOutput(2, TTree::Class());
}


AliAnalysisTaskExtractPerformanceV0::~AliAnalysisTaskExtractPerformanceV0()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------

   if (fListHistV0){
      delete fListHistV0;
      fListHistV0 = 0x0;
   }
   if (fTree){
      delete fTree;
      fTree = 0x0;
   }
    //cleanup esd track cuts object too...
   if (fESDtrackCuts){
    delete fESDtrackCuts;
    fESDtrackCuts = 0x0; 
  }
}

//________________________________________________________________________
void AliAnalysisTaskExtractPerformanceV0::UserCreateOutputObjects()
{

   OpenFile(2);	
   // Called once

//------------------------------------------------

   fTree = new TTree("fTree","V0Candidates");

//------------------------------------------------
// fTree Branch definitions - V0 Tree
//------------------------------------------------

//-----------BASIC-INFO---------------------------
/* 1*/   fTree->Branch("fTreeVariablePrimaryStatus",&fTreeVariablePrimaryStatus,"fTreeVariablePrimaryStatus/I");	
/* 2*/   fTree->Branch("fTreeVariablePrimaryStatusMother",&fTreeVariablePrimaryStatusMother,"fTreeVariablePrimaryStatusMother/I");	
/* 2*/   fTree->Branch("fTreeVariableChi2V0",&fTreeVariableChi2V0,"Chi2V0/F");
/* 3*/   fTree->Branch("fTreeVariableDcaV0Daughters",&fTreeVariableDcaV0Daughters,"fTreeVariableDcaV0Daughters/F");
/* 4*/   fTree->Branch("fTreeVariableDcaPosToPrimVertex",&fTreeVariableDcaPosToPrimVertex,"fTreeVariableDcaPosToPrimVertex/F");
/* 5*/   fTree->Branch("fTreeVariableDcaNegToPrimVertex",&fTreeVariableDcaNegToPrimVertex,"fTreeVariableDcaNegToPrimVertex/F");
/* 6*/   fTree->Branch("fTreeVariableV0Radius",&fTreeVariableV0Radius,"fTreeVariableV0Radius/F");
/* 7*/   fTree->Branch("fTreeVariablePt",&fTreeVariablePt,"fTreeVariablePt/F");
/* 7*/   fTree->Branch("fTreeVariablePtMC",&fTreeVariablePtMC,"fTreeVariablePtMC/F");
/* 8*/   fTree->Branch("fTreeVariableRapK0Short",&fTreeVariableRapK0Short,"fTreeVariableRapK0Short/F");
/* 9*/   fTree->Branch("fTreeVariableRapLambda",&fTreeVariableRapLambda,"fTreeVariableRapLambda/F");
/*10*/   fTree->Branch("fTreeVariableRapMC",&fTreeVariableRapMC,"fTreeVariableRapMC/F");
/*11*/   fTree->Branch("fTreeVariableInvMassK0s",&fTreeVariableInvMassK0s,"fTreeVariableInvMassK0s/F");
/*12*/   fTree->Branch("fTreeVariableInvMassLambda",&fTreeVariableInvMassLambda,"fTreeVariableInvMassLambda/F");
/*13*/   fTree->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
/*14*/   fTree->Branch("fTreeVariableAlphaV0",&fTreeVariableAlphaV0,"fTreeVariableAlphaV0/F");
/*15*/   fTree->Branch("fTreeVariablePtArmV0",&fTreeVariablePtArmV0,"fTreeVariablePtArmV0/F");
/*16*/   fTree->Branch("fTreeVariableNegTransvMomentum",&fTreeVariableNegTransvMomentum,"fTreeVariableNegTransvMomentum/F");
/*17*/   fTree->Branch("fTreeVariablePosTransvMomentum",&fTreeVariablePosTransvMomentum,"fTreeVariablePosTransvMomentum/F");
/*18*/   fTree->Branch("fTreeVariableNegTransvMomentumMC",&fTreeVariableNegTransvMomentumMC,"fTreeVariableNegTransvMomentumMC/F");
/*19*/   fTree->Branch("fTreeVariablePosTransvMomentumMC",&fTreeVariablePosTransvMomentumMC,"fTreeVariablePosTransvMomentumMC/F");
/*20*/   fTree->Branch("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastNbrCrossedRows,"fTreeVariableLeastNbrCrossedRows/I");
/*21*/   fTree->Branch("fTreeVariableLeastRatioCrossedRowsOverFindable",&fTreeVariableLeastRatioCrossedRowsOverFindable,"fTreeVariableLeastRatioCrossedRowsOverFindable/F");
/*22*/   fTree->Branch("fTreeVariablePID",&fTreeVariablePID,"fTreeVariablePID/I");
/*23*/   fTree->Branch("fTreeVariablePIDPositive",&fTreeVariablePIDPositive,"fTreeVariablePIDPositive/I");
/*24*/   fTree->Branch("fTreeVariablePIDNegative",&fTreeVariablePIDNegative,"fTreeVariablePIDNegative/I");
/*25*/   fTree->Branch("fTreeVariablePIDMother",&fTreeVariablePIDMother,"fTreeVariablePIDMother/I");
/*26*/   fTree->Branch("fTreeVariablePtXiMother",&fTreeVariablePtMother,"fTreeVariablePtMother/F");
/*27*/   fTree->Branch("fTreeVariableV0CosineOfPointingAngle",&fTreeVariableV0CosineOfPointingAngle,"fTreeVariableV0CosineOfPointingAngle/F");
//-----------MULTIPLICITY-INFO--------------------
/*28*/   fTree->Branch("fTreeVariableMultiplicity",&fTreeVariableMultiplicity,"fTreeVariableMultiplicity/I");
/*28*/   fTree->Branch("fTreeVariableMultiplicityMC",&fTreeVariableMultiplicityMC,"fTreeVariableMultiplicityMC/I");
  /*17*/	fTree->Branch("fTreeVariableMultiplicityV0A",&fTreeVariableMultiplicityV0A,"fTreeVariableMultiplicityV0A/I");
  /*17*/	fTree->Branch("fTreeVariableMultiplicityZNA",&fTreeVariableMultiplicityZNA,"fTreeVariableMultiplicityZNA/I");
  /*17*/	fTree->Branch("fTreeVariableMultiplicityTRK",&fTreeVariableMultiplicityTRK,"fTreeVariableMultiplicityTRK/I");
  /*17*/	fTree->Branch("fTreeVariableMultiplicitySPD",&fTreeVariableMultiplicitySPD,"fTreeVariableMultiplicitySPD/I");
//------------------------------------------------
/*29*/   fTree->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
/*30*/   fTree->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
/*31*/   fTree->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
/*32*/   fTree->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
/*33*/   fTree->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
//------------------------------------------------
/*34*/   fTree->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
/*35*/   fTree->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
/*36*/   fTree->Branch("fTreeVariableV0CreationRadius",&fTreeVariableV0CreationRadius,"fTreeVariableV0CreationRadius/F");
  
    if ( fkLightWeight == kFALSE ){ // these are debugging branches!
/*37*/   fTree->Branch("fTreeVariableIndexStatus",&fTreeVariableIndexStatus,"fTreeVariableIndexStatus/I");
/*38*/   fTree->Branch("fTreeVariableIndexStatusMother",&fTreeVariableIndexStatusMother,"fTreeVariableIndexStatusMother/I");
    }

/*39*/ 	 fTree->Branch("fTreeVariableRunNumber",&fTreeVariableRunNumber,"fTreeVariableRunNumber/I");
/*40*/   fTree->Branch("fTreeVariableEventNumber",&fTreeVariableEventNumber,"fTreeVariableEventNumber/l");

    if ( fkLightWeight == kFALSE ){ // these are debugging branches!
/*34*/   fTree->Branch("fTreeVariableVertexZ",&fTreeVariableVertexZ,"fTreeVariableVertexZ/F");
    }
      
    if ( fkLightWeight == kFALSE ){ // these are debugging branches!
//-----------FOR CTAU DEBUGGING: Full Phase Space + Decay Position Info 
        fTree->Branch("fTreeVariableV0x",&fTreeVariableV0x,"fTreeVariableV0x/F");
        fTree->Branch("fTreeVariableV0y",&fTreeVariableV0y,"fTreeVariableV0y/F");
        fTree->Branch("fTreeVariableV0z",&fTreeVariableV0z,"fTreeVariableV0z/F");

        fTree->Branch("fTreeVariableV0Px",&fTreeVariableV0Px,"fTreeVariableV0Px/F");
        fTree->Branch("fTreeVariableV0Py",&fTreeVariableV0Py,"fTreeVariableV0Py/F");
        fTree->Branch("fTreeVariableV0Pz",&fTreeVariableV0Pz,"fTreeVariableV0Pz/F");

//-----------FOR CTAU DEBUGGING: Full Phase Space + Decay Position Info, perfect info from MC
        fTree->Branch("fTreeVariableMCV0x",&fTreeVariableMCV0x,"fTreeVariableMCV0x/F");
        fTree->Branch("fTreeVariableMCV0y",&fTreeVariableMCV0y,"fTreeVariableMCV0y/F");
        fTree->Branch("fTreeVariableMCV0z",&fTreeVariableMCV0z,"fTreeVariableMCV0z/F");

        fTree->Branch("fTreeVariableMCV0Px",&fTreeVariableMCV0Px,"fTreeVariableMCV0Px/F");
        fTree->Branch("fTreeVariableMCV0Py",&fTreeVariableMCV0Py,"fTreeVariableMCV0Py/F");
        fTree->Branch("fTreeVariableMCV0Pz",&fTreeVariableMCV0Pz,"fTreeVariableMCV0Pz/F");

//-----------FOR CTAU DEBUGGING: Primary vertex info 
        fTree->Branch("fTreeVariablePVx",&fTreeVariablePVx,"fTreeVariablePVx/F");
        fTree->Branch("fTreeVariablePVy",&fTreeVariablePVy,"fTreeVariablePVy/F");
        fTree->Branch("fTreeVariablePVz",&fTreeVariablePVz,"fTreeVariablePVz/F");

        fTree->Branch("fTreeVariableMCPVx",&fTreeVariableMCPVx,"fTreeVariableMCPVx/F");
        fTree->Branch("fTreeVariableMCPVy",&fTreeVariableMCPVy,"fTreeVariableMCPVy/F");
        fTree->Branch("fTreeVariableMCPVz",&fTreeVariableMCPVz,"fTreeVariableMCPVz/F");
    }

        fTree->Branch("fTreeVariableIsNonInjected",&fTreeVariableIsNonInjected,"fTreeVariableIsNonInjected/O"); //O for bOOlean...
  
  if ( fkLightWeight == kFALSE ){ // these are debugging branches!
        fTree->Branch("fTreeVariableNegTrackStatus",&fTreeVariableNegTrackStatus,"fTreeVariableNegTrackStatus/l");
        fTree->Branch("fTreeVariablePosTrackStatus",&fTreeVariablePosTrackStatus,"fTreeVariablePosTrackStatus/l");
  }
  
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
   fListHistV0 = new TList();
   fListHistV0->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner


   if(! fHistV0MultiplicityBeforeTrigSel) {
      fHistV0MultiplicityBeforeTrigSel = new TH1F("fHistV0MultiplicityBeforeTrigSel", 
         "V0s per event (before Trig. Sel.);Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityBeforeTrigSel);
   }
           
   if(! fHistV0MultiplicityForTrigEvt) {
      fHistV0MultiplicityForTrigEvt = new TH1F("fHistV0MultiplicityForTrigEvt", 
         "V0s per event (for triggered evt);Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityForTrigEvt);
   }

   if(! fHistV0MultiplicityForSelEvt) {
      fHistV0MultiplicityForSelEvt = new TH1F("fHistV0MultiplicityForSelEvt", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityForSelEvt);
   }

   if(! fHistV0MultiplicityForSelEvtNoTPCOnly) {
      fHistV0MultiplicityForSelEvtNoTPCOnly = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnly", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityForSelEvtNoTPCOnly);
   }
   if(! fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup) {
      fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup);
   }

//------------------------------------------------
// Track Multiplicity Histograms
//------------------------------------------------

   if(! fHistMultiplicityBeforeTrigSel) {
      fHistMultiplicityBeforeTrigSel = new TH1F("fHistMultiplicityBeforeTrigSel", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicityBeforeTrigSel);
   }
   if(! fHistMultiplicityForTrigEvt) {
      fHistMultiplicityForTrigEvt = new TH1F("fHistMultiplicityForTrigEvt", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicityForTrigEvt);
   }
   if(! fHistMultiplicity) {
      fHistMultiplicity = new TH1F("fHistMultiplicity", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicity);
   }
   if(! fHistMultiplicityNoTPCOnly) {
      fHistMultiplicityNoTPCOnly = new TH1F("fHistMultiplicityNoTPCOnly", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicityNoTPCOnly);
   }
   if(! fHistMultiplicityNoTPCOnlyNoPileup) {
      fHistMultiplicityNoTPCOnlyNoPileup = new TH1F("fHistMultiplicityNoTPCOnlyNoPileup", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicityNoTPCOnlyNoPileup);
   }
  
  //V0A Centrality (if PbPb / pPb)
  if(! fHistMultiplicityV0ABeforeTrigSel) {
    fHistMultiplicityV0ABeforeTrigSel = new TH1F("fHistMultiplicityV0ABeforeTrigSel",
                                                 "Centrality Distribution: V0A;V0A Centrality;Events",
                                                 200, 0, 200);
    fListHistV0->Add(fHistMultiplicityV0ABeforeTrigSel);
  }
  if(! fHistMultiplicityV0AForTrigEvt) {
    fHistMultiplicityV0AForTrigEvt = new TH1F("fHistMultiplicityV0AForTrigEvt",
                                              "Centrality Distribution: V0A;V0A Centrality;Events",
                                              200, 0, 200);
    fListHistV0->Add(fHistMultiplicityV0AForTrigEvt);
  }
  if(! fHistMultiplicityV0A) {
    fHistMultiplicityV0A = new TH1F("fHistMultiplicityV0A",
                                    "Centrality Distribution: V0A;V0A Centrality;Events",
                                    200, 0, 200);
    fListHistV0->Add(fHistMultiplicityV0A);
  }
  if(! fHistMultiplicityV0ANoTPCOnly) {
    fHistMultiplicityV0ANoTPCOnly = new TH1F("fHistMultiplicityV0ANoTPCOnly",
                                             "Centrality Distribution: V0A;V0A Centrality;Events",
                                             200, 0, 200);
    fListHistV0->Add(fHistMultiplicityV0ANoTPCOnly);
  }
  if(! fHistMultiplicityV0ANoTPCOnlyNoPileup) {
    fHistMultiplicityV0ANoTPCOnlyNoPileup = new TH1F("fHistMultiplicityV0ANoTPCOnlyNoPileup",
                                                     "Centrality Distribution: V0A;V0A Centrality;Events",
                                                     200, 0, 200);
    fListHistV0->Add(fHistMultiplicityV0ANoTPCOnlyNoPileup);
  }
  
  //ZNA Centrality (if PbPb / pPb)
  if(! fHistMultiplicityZNABeforeTrigSel) {
    fHistMultiplicityZNABeforeTrigSel = new TH1F("fHistMultiplicityZNABeforeTrigSel",
                                                 "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                                 200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNABeforeTrigSel);
  }
  if(! fHistMultiplicityZNAForTrigEvt) {
    fHistMultiplicityZNAForTrigEvt = new TH1F("fHistMultiplicityZNAForTrigEvt",
                                              "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                              200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNAForTrigEvt);
  }
  if(! fHistMultiplicityZNA) {
    fHistMultiplicityZNA = new TH1F("fHistMultiplicityZNA",
                                    "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                    200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNA);
  }
  if(! fHistMultiplicityZNANoTPCOnly) {
    fHistMultiplicityZNANoTPCOnly = new TH1F("fHistMultiplicityZNANoTPCOnly",
                                             "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                             200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNANoTPCOnly);
  }
  if(! fHistMultiplicityZNANoTPCOnlyNoPileup) {
    fHistMultiplicityZNANoTPCOnlyNoPileup = new TH1F("fHistMultiplicityZNANoTPCOnlyNoPileup",
                                                     "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                                     200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNANoTPCOnlyNoPileup);
  }
  
  //TRK Centrality (if PbPb / pPb)
  if(! fHistMultiplicityTRKBeforeTrigSel) {
    fHistMultiplicityTRKBeforeTrigSel = new TH1F("fHistMultiplicityTRKBeforeTrigSel",
                                                 "Centrality Distribution: TRK;TRK Centrality;Events",
                                                 200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRKBeforeTrigSel);
  }
  if(! fHistMultiplicityTRKForTrigEvt) {
    fHistMultiplicityTRKForTrigEvt = new TH1F("fHistMultiplicityTRKForTrigEvt",
                                              "Centrality Distribution: TRK;TRK Centrality;Events",
                                              200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRKForTrigEvt);
  }
  if(! fHistMultiplicityTRK) {
    fHistMultiplicityTRK = new TH1F("fHistMultiplicityTRK",
                                    "Centrality Distribution: TRK;TRK Centrality;Events",
                                    200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRK);
  }
  if(! fHistMultiplicityTRKNoTPCOnly) {
    fHistMultiplicityTRKNoTPCOnly = new TH1F("fHistMultiplicityTRKNoTPCOnly",
                                             "Centrality Distribution: TRK;TRK Centrality;Events",
                                             200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRKNoTPCOnly);
  }
  if(! fHistMultiplicityTRKNoTPCOnlyNoPileup) {
    fHistMultiplicityTRKNoTPCOnlyNoPileup = new TH1F("fHistMultiplicityTRKNoTPCOnlyNoPileup",
                                                     "Centrality Distribution: TRK;TRK Centrality;Events",
                                                     200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRKNoTPCOnlyNoPileup);
  }
  
  //SPD Centrality (if PbPb / pPb)
  if(! fHistMultiplicitySPDBeforeTrigSel) {
    fHistMultiplicitySPDBeforeTrigSel = new TH1F("fHistMultiplicitySPDBeforeTrigSel",
                                                 "Centrality Distribution: SPD;SPD Centrality;Events",
                                                 200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPDBeforeTrigSel);
  }
  if(! fHistMultiplicitySPDForTrigEvt) {
    fHistMultiplicitySPDForTrigEvt = new TH1F("fHistMultiplicitySPDForTrigEvt",
                                              "Centrality Distribution: SPD;SPD Centrality;Events",
                                              200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPDForTrigEvt);
  }
  if(! fHistMultiplicitySPD) {
    fHistMultiplicitySPD = new TH1F("fHistMultiplicitySPD",
                                    "Centrality Distribution: SPD;SPD Centrality;Events",
                                    200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPD);
  }
  if(! fHistMultiplicitySPDNoTPCOnly) {
    fHistMultiplicitySPDNoTPCOnly = new TH1F("fHistMultiplicitySPDNoTPCOnly",
                                             "Centrality Distribution: SPD;SPD Centrality;Events",
                                             200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPDNoTPCOnly);
  }
  if(! fHistMultiplicitySPDNoTPCOnlyNoPileup) {
    fHistMultiplicitySPDNoTPCOnlyNoPileup = new TH1F("fHistMultiplicitySPDNoTPCOnlyNoPileup",
                                                     "Centrality Distribution: SPD;SPD Centrality;Events",
                                                     200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPDNoTPCOnlyNoPileup);
  }
  

  

  //Raw Data for J/Psi paper Technique
	//TH2F    *f2dHistMultiplicityVsTrueBeforeTrigSel; 	        //! multiplicity distribution    
	//TH2F    *f2dHistMultiplicityVsTrueForTrigEvt;  		        //! multiplicity distribution
	//TH2F    *f2dHistMultiplicityVsTrue;     					        //! multiplicity distribution
	//TH2F    *f2dHistMultiplicityVsTrueNoTPCOnly;			        //! multiplicity distribution
	//TH2F    *f2dHistMultiplicityVsTrueNoTPCOnlyNoPileup;			//! multiplicity distribution

   if(! f2dHistMultiplicityVsTrueBeforeTrigSel) {
      f2dHistMultiplicityVsTrueBeforeTrigSel = new TH2F("f2dHistMultiplicityVsTrueBeforeTrigSel", 
         "Tracks per event", 200, 0, 200, 200, 0, 200); 		
      fListHistV0->Add(f2dHistMultiplicityVsTrueBeforeTrigSel);
   }
   if(! f2dHistMultiplicityVsTrueForTrigEvt) {
      f2dHistMultiplicityVsTrueForTrigEvt = new TH2F("f2dHistMultiplicityVsTrueForTrigEvt", 
         "Tracks per event", 200, 0, 200, 200, 0, 200); 		
      fListHistV0->Add(f2dHistMultiplicityVsTrueForTrigEvt);
   }
   if(! f2dHistMultiplicityVsTrue) {
      f2dHistMultiplicityVsTrue = new TH2F("f2dHistMultiplicityVsTrue", 
         "Tracks per event", 200, 0, 200, 200, 0, 200); 		
      fListHistV0->Add(f2dHistMultiplicityVsTrue);
   }
   if(! f2dHistMultiplicityVsTrueNoTPCOnly) {
      f2dHistMultiplicityVsTrueNoTPCOnly = new TH2F("f2dHistMultiplicityVsTrueNoTPCOnly", 
         "Tracks per event", 200, 0, 200, 200, 0, 200); 		
      fListHistV0->Add(f2dHistMultiplicityVsTrueNoTPCOnly);
   }
   if(! f2dHistMultiplicityVsTrueNoTPCOnlyNoPileup) {
      f2dHistMultiplicityVsTrueNoTPCOnlyNoPileup = new TH2F("f2dHistMultiplicityVsTrueNoTPCOnlyNoPileup", 
         "Tracks per event", 200, 0, 200, 200, 0, 200); 		
      fListHistV0->Add(f2dHistMultiplicityVsTrueNoTPCOnlyNoPileup);
   }


  //Raw Data for Vertex Z position estimator change
	//TH2F    *f2dHistMultiplicityVsVertexZBeforeTrigSel; 	        //! multiplicity distribution    
	//TH2F    *f2dHistMultiplicityVsVertexZForTrigEvt;  		        //! multiplicity distribution
	//TH2F    *f2dHistMultiplicityVsVertexZ;     					        //! multiplicity distribution
	//TH2F    *f2dHistMultiplicityVsVertexZNoTPCOnly;			        //! multiplicity distribution
	//TH2F    *f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup;			//! multiplicity distribution

   if(! f2dHistMultiplicityVsVertexZBeforeTrigSel) {
      f2dHistMultiplicityVsVertexZBeforeTrigSel = new TH2F("f2dHistMultiplicityVsVertexZBeforeTrigSel", 
         "Tracks per event", 200, 0, 200,400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZBeforeTrigSel);
   }
   if(! f2dHistMultiplicityVsVertexZForTrigEvt) {
      f2dHistMultiplicityVsVertexZForTrigEvt = new TH2F("f2dHistMultiplicityVsVertexZForTrigEvt", 
         "Tracks per event", 200, 0, 200, 400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZForTrigEvt);
   }
   if(! f2dHistMultiplicityVsVertexZ) {
      f2dHistMultiplicityVsVertexZ = new TH2F("f2dHistMultiplicityVsVertexZ", 
         "Tracks per event", 200, 0, 200, 400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZ);
   }
   if(! f2dHistMultiplicityVsVertexZNoTPCOnly) {
      f2dHistMultiplicityVsVertexZNoTPCOnly = new TH2F("f2dHistMultiplicityVsVertexZNoTPCOnly", 
         "Tracks per event", 200, 0, 200, 400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZNoTPCOnly);
   }
   if(! f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup) {
      f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup = new TH2F("f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup", 
         "Tracks per event", 200, 0, 200, 400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup);
   }

//Generated PVz test

//   TH1F      *fHistGenVertexZBeforeTrigSel;     //! multiplicity distribution      
//   TH1F      *fHistGenVertexZForTrigEvt;        //! multiplicity distribution
//   TH1F      *fHistGenVertexZ;                  //! multiplicity distribution
//   TH1F      *fHistGenVertexZNoTPCOnly;         //! multiplicity distribution
//   TH1F      *fHistGenVertexZNoTPCOnlyNoPileup; //! multiplicity distribution

   if(! fHistGenVertexZBeforeTrigSel) {
         fHistGenVertexZBeforeTrigSel = new TH1F("fHistGenVertexZBeforeTrigSel", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHistV0->Add(fHistGenVertexZBeforeTrigSel);
   }
   if(! fHistGenVertexZForTrigEvt) {
         fHistGenVertexZForTrigEvt = new TH1F("fHistGenVertexZForTrigEvt", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHistV0->Add(fHistGenVertexZForTrigEvt);
   }
   if(! fHistGenVertexZ) {
         fHistGenVertexZ = new TH1F("fHistGenVertexZ", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHistV0->Add(fHistGenVertexZ);
   }
   if(! fHistGenVertexZNoTPCOnly) {
         fHistGenVertexZNoTPCOnly = new TH1F("fHistGenVertexZNoTPCOnly", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHistV0->Add(fHistGenVertexZNoTPCOnly);
   }
   if(! fHistGenVertexZNoTPCOnlyNoPileup) {
         fHistGenVertexZNoTPCOnlyNoPileup = new TH1F("fHistGenVertexZNoTPCOnlyNoPileup", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHistV0->Add(fHistGenVertexZNoTPCOnlyNoPileup);
   }


//------------------------------------------------
// Generated Particle Histograms
//------------------------------------------------

   Int_t lCustomNBins = 200; 
   Double_t lCustomPtUpperLimit = 20; 
   Int_t lCustomNBinsMultiplicity = 100;

//----------------------------------
// Raw Generated (Pre-physics-selection)
//----------------------------------

//--- 3D Histo (Pt, Y, Multiplicity)  

   if(! f3dHistPrimRawPtVsYVsMultLambda) {
      f3dHistPrimRawPtVsYVsMultLambda = new TH3F( "f3dHistPrimRawPtVsYVsMultLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultLambda);
   }
   if(! f3dHistPrimRawPtVsYVsMultAntiLambda) {
      f3dHistPrimRawPtVsYVsMultAntiLambda = new TH3F( "f3dHistPrimRawPtVsYVsMultAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultAntiLambda);
   }
   if(! f3dHistPrimRawPtVsYVsMultK0Short) {
      f3dHistPrimRawPtVsYVsMultK0Short = new TH3F( "f3dHistPrimRawPtVsYVsMultK0Short", "Pt_{K0S} Vs Y_{K0S} Vs Multiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultK0Short);
   }

   if(! f3dHistPrimRawPtVsYCMSVsMultLambda) {
      f3dHistPrimRawPtVsYCMSVsMultLambda = new TH3F( "f3dHistPrimRawPtVsYCMSVsMultLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYCMSVsMultLambda);
   }
   if(! f3dHistPrimRawPtVsYCMSVsMultAntiLambda) {
      f3dHistPrimRawPtVsYCMSVsMultAntiLambda = new TH3F( "f3dHistPrimRawPtVsYCMSVsMultAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYCMSVsMultAntiLambda);
   }
   if(! f3dHistPrimRawPtVsYCMSVsMultK0Short) {
      f3dHistPrimRawPtVsYCMSVsMultK0Short = new TH3F( "f3dHistPrimRawPtVsYCMSVsMultK0Short", "Pt_{K0S} Vs Y_{K0S} Vs Multiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYCMSVsMultK0Short);
   }

//---> Non-injected particles

   if(! f3dHistPrimRawPtVsYVsMultNonInjLambda) {
      f3dHistPrimRawPtVsYVsMultNonInjLambda = new TH3F( "f3dHistPrimRawPtVsYVsMultNonInjLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultNonInjLambda);
   }
   if(! f3dHistPrimRawPtVsYVsMultNonInjAntiLambda) {
      f3dHistPrimRawPtVsYVsMultNonInjAntiLambda = new TH3F( "f3dHistPrimRawPtVsYVsMultNonInjAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultNonInjAntiLambda);
   }
   if(! f3dHistPrimRawPtVsYVsMultNonInjK0Short) {
      f3dHistPrimRawPtVsYVsMultNonInjK0Short = new TH3F( "f3dHistPrimRawPtVsYVsMultNonInjK0Short", "Pt_{K0S} Vs Y_{K0S} Vs Multiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultNonInjK0Short);
   }

//--- 3D Histo (Pt, Y, MultiplicityMC)  

   if(! f3dHistPrimRawPtVsYVsMultMCLambda) {
      f3dHistPrimRawPtVsYVsMultMCLambda = new TH3F( "f3dHistPrimRawPtVsYVsMultMCLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; MultMC", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultMCLambda);
   }
   if(! f3dHistPrimRawPtVsYVsMultMCAntiLambda) {
      f3dHistPrimRawPtVsYVsMultMCAntiLambda = new TH3F( "f3dHistPrimRawPtVsYVsMultMCAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; MultMC", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultMCAntiLambda);
   }
   if(! f3dHistPrimRawPtVsYVsMultMCK0Short) {
      f3dHistPrimRawPtVsYVsMultMCK0Short = new TH3F( "f3dHistPrimRawPtVsYVsMultMCK0Short", "Pt_{K0S} Vs Y_{K0S} Vs Multiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; MultMC", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultMCK0Short);
   }

//--- 3D Histo (Pt, Y, VertexZ)  

   if(! f3dHistPrimRawPtVsYVsVertexZLambda) {
      f3dHistPrimRawPtVsYVsVertexZLambda = new TH3F( "f3dHistPrimRawPtVsYVsVertexZLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs VertexZiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; VertexZ", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,40,-10,10);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsVertexZLambda);
   }
   if(! f3dHistPrimRawPtVsYVsVertexZAntiLambda) {
      f3dHistPrimRawPtVsYVsVertexZAntiLambda = new TH3F( "f3dHistPrimRawPtVsYVsVertexZAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs VertexZiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; VertexZ", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,40,-10,10);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsVertexZAntiLambda);
   }
   if(! f3dHistPrimRawPtVsYVsVertexZK0Short) {
      f3dHistPrimRawPtVsYVsVertexZK0Short = new TH3F( "f3dHistPrimRawPtVsYVsVertexZK0Short", "Pt_{K0S} Vs Y_{K0S} Vs VertexZiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; VertexZ", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,40,-10,10);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsVertexZK0Short);
   }

//--- 3D Histo (Pt, Y, Multiplicity), close to PV criterion

   if(! f3dHistPrimCloseToPVPtVsYVsMultLambda) {
      f3dHistPrimCloseToPVPtVsYVsMultLambda = new TH3F( "f3dHistPrimCloseToPVPtVsYVsMultLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimCloseToPVPtVsYVsMultLambda);
   }
   if(! f3dHistPrimCloseToPVPtVsYVsMultAntiLambda) {
      f3dHistPrimCloseToPVPtVsYVsMultAntiLambda = new TH3F( "f3dHistPrimCloseToPVPtVsYVsMultAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimCloseToPVPtVsYVsMultAntiLambda);
   }
   if(! f3dHistPrimCloseToPVPtVsYVsMultK0Short) {
      f3dHistPrimCloseToPVPtVsYVsMultK0Short = new TH3F( "f3dHistPrimCloseToPVPtVsYVsMultK0Short", "Pt_{K0S} Vs Y_{K0S} Vs Multiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimCloseToPVPtVsYVsMultK0Short);
   }


//--- 3D Histo (Pt, Y, Proper Decay Length)

   if(! f3dHistPrimRawPtVsYVsDecayLengthLambda) {
      f3dHistPrimRawPtVsYVsDecayLengthLambda = new TH3F( "f3dHistPrimRawPtVsYVsDecayLengthLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs DecayLength; Pt_{lambda} (GeV/c); Y_{#Lambda} ; DecayLength", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,200,0,50);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsDecayLengthLambda);
   }
   if(! f3dHistPrimRawPtVsYVsDecayLengthAntiLambda) {
      f3dHistPrimRawPtVsYVsDecayLengthAntiLambda = new TH3F( "f3dHistPrimRawPtVsYVsDecayLengthAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs DecayLength; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; DecayLength", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,200,0,50);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsDecayLengthAntiLambda);
   }
   if(! f3dHistPrimRawPtVsYVsDecayLengthK0Short) {
      f3dHistPrimRawPtVsYVsDecayLengthK0Short = new TH3F( "f3dHistPrimRawPtVsYVsDecayLengthK0Short", "Pt_{K0S} Vs Y_{K0S} Vs DecayLength; Pt_{K0S} (GeV/c); Y_{K0S} ; DecayLength", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,200,0,50);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsDecayLengthK0Short);
   }

//--------------------------------------------------------------------------------------
//--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, all generated

   if(! f3dHistGenPtVsYVsMultXiMinus) {
      f3dHistGenPtVsYVsMultXiMinus = new TH3F( "f3dHistGenPtVsYVsMultXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYVsMultXiMinus);
   }
   if(! f3dHistGenPtVsYVsMultXiPlus) {
      f3dHistGenPtVsYVsMultXiPlus = new TH3F( "f3dHistGenPtVsYVsMultXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYVsMultXiPlus);
   }
//--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus

   if(! f3dHistGenPtVsYVsMultOmegaMinus) {
      f3dHistGenPtVsYVsMultOmegaMinus = new TH3F( "f3dHistGenPtVsYVsMultOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYVsMultOmegaMinus);
   }
   if(! f3dHistGenPtVsYVsMultOmegaPlus) {
      f3dHistGenPtVsYVsMultOmegaPlus = new TH3F( "f3dHistGenPtVsYVsMultOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYVsMultOmegaPlus);
   }

//CASCADEs, Y CMS 

   if(! f3dHistGenPtVsYCMSVsMultXiMinus) {
      f3dHistGenPtVsYCMSVsMultXiMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYCMSVsMultXiMinus);
   }
   if(! f3dHistGenPtVsYCMSVsMultXiPlus) {
      f3dHistGenPtVsYCMSVsMultXiPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYCMSVsMultXiPlus);
   }
//--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus

   if(! f3dHistGenPtVsYCMSVsMultOmegaMinus) {
      f3dHistGenPtVsYCMSVsMultOmegaMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYCMSVsMultOmegaMinus);
   }
   if(! f3dHistGenPtVsYCMSVsMultOmegaPlus) {
      f3dHistGenPtVsYCMSVsMultOmegaPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYCMSVsMultOmegaPlus);
   }



//--------------------------------------------------------------------------------------
//--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, at selected analysis evts

   if(! f3dHistGenSelectedPtVsYVsMultXiMinus) {
      f3dHistGenSelectedPtVsYVsMultXiMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenSelectedPtVsYVsMultXiMinus);
   }
   if(! f3dHistGenSelectedPtVsYVsMultXiPlus) {
      f3dHistGenSelectedPtVsYVsMultXiPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenSelectedPtVsYVsMultXiPlus);
   }
//--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus

   if(! f3dHistGenSelectedPtVsYVsMultOmegaMinus) {
      f3dHistGenSelectedPtVsYVsMultOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenSelectedPtVsYVsMultOmegaMinus);
   }
   if(! f3dHistGenSelectedPtVsYVsMultOmegaPlus) {
      f3dHistGenSelectedPtVsYVsMultOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenSelectedPtVsYVsMultOmegaPlus);
   }

//CASCADES, analysis level, y CMS

//--------------------------------------------------------------------------------------
//--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, at selected analysis evts

   if(! f3dHistGenSelectedPtVsYCMSVsMultXiMinus) {
      f3dHistGenSelectedPtVsYCMSVsMultXiMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenSelectedPtVsYCMSVsMultXiMinus);
   }
   if(! f3dHistGenSelectedPtVsYCMSVsMultXiPlus) {
      f3dHistGenSelectedPtVsYCMSVsMultXiPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenSelectedPtVsYCMSVsMultXiPlus);
   }
//--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus

   if(! f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus) {
      f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus);
   }
   if(! f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus) {
      f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus);
   }


//----------------------------------
// Histos at analysis level 
//----------------------------------

   if(! f3dHistPrimAnalysisPtVsYVsMultLambda) {
      f3dHistPrimAnalysisPtVsYVsMultLambda = new TH3F( "f3dHistPrimAnalysisPtVsYVsMultLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimAnalysisPtVsYVsMultLambda);
   }
   if(! f3dHistPrimAnalysisPtVsYVsMultAntiLambda) {
      f3dHistPrimAnalysisPtVsYVsMultAntiLambda = new TH3F( "f3dHistPrimAnalysisPtVsYVsMultAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimAnalysisPtVsYVsMultAntiLambda);
   }
   if(! f3dHistPrimAnalysisPtVsYVsMultK0Short) {
      f3dHistPrimAnalysisPtVsYVsMultK0Short = new TH3F( "f3dHistPrimAnalysisPtVsYVsMultK0Short", "Pt_{K0S} Vs Y_{K0S} Vs Multiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimAnalysisPtVsYVsMultK0Short);
   }

   if(! f3dHistPrimAnalysisPtVsYCMSVsMultLambda) {
      f3dHistPrimAnalysisPtVsYCMSVsMultLambda = new TH3F( "f3dHistPrimAnalysisPtVsYCMSVsMultLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimAnalysisPtVsYCMSVsMultLambda);
   }
   if(! f3dHistPrimAnalysisPtVsYCMSVsMultAntiLambda) {
      f3dHistPrimAnalysisPtVsYCMSVsMultAntiLambda = new TH3F( "f3dHistPrimAnalysisPtVsYCMSVsMultAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimAnalysisPtVsYCMSVsMultAntiLambda);
   }
   if(! f3dHistPrimAnalysisPtVsYCMSVsMultK0Short) {
      f3dHistPrimAnalysisPtVsYCMSVsMultK0Short = new TH3F( "f3dHistPrimAnalysisPtVsYCMSVsMultK0Short", "Pt_{K0S} Vs Y_{K0S} Vs Multiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimAnalysisPtVsYCMSVsMultK0Short);
   }

//----------------------------------
// Primary Vertex Position Histos
//----------------------------------

   if(! fHistPVx) {
         fHistPVx = new TH1F("fHistPVx", 
            "PV x position;Nbr of Evts;x", 
            2000, -0.5, 0.5);       
      fListHistV0->Add(fHistPVx);
   }
   if(! fHistPVy) {
         fHistPVy = new TH1F("fHistPVy", 
            "PV y position;Nbr of Evts;y", 
            2000, -0.5, 0.5);       
      fListHistV0->Add(fHistPVy);
   }
   if(! fHistPVz) {
         fHistPVz = new TH1F("fHistPVz", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHistV0->Add(fHistPVz);
   }

   if(! fHistPVxAnalysis) {
         fHistPVxAnalysis = new TH1F("fHistPVxAnalysis", 
            "PV x position;Nbr of Evts;x", 
            2000, -0.5, 0.5);       
      fListHistV0->Add(fHistPVxAnalysis);
   }
   if(! fHistPVyAnalysis) {
         fHistPVyAnalysis = new TH1F("fHistPVyAnalysis", 
            "PV y position;Nbr of Evts;y", 
            2000, -0.5, 0.5);       
      fListHistV0->Add(fHistPVyAnalysis);
   }
   if(! fHistPVzAnalysis) {
         fHistPVzAnalysis = new TH1F("fHistPVzAnalysis", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHistV0->Add(fHistPVzAnalysis);
   }

   if(! fHistPVxAnalysisHasHighPtLambda) {
         fHistPVxAnalysisHasHighPtLambda = new TH1F("fHistPVxAnalysisHasHighPtLambda", 
            "PV x position;Nbr of Evts;x", 
            2000, -0.5, 0.5);       
      fListHistV0->Add(fHistPVxAnalysisHasHighPtLambda);
   }
   if(! fHistPVyAnalysisHasHighPtLambda) {
         fHistPVyAnalysisHasHighPtLambda = new TH1F("fHistPVyAnalysisHasHighPtLambda", 
            "PV y position;Nbr of Evts;y", 
            2000, -0.5, 0.5);       
      fListHistV0->Add(fHistPVyAnalysisHasHighPtLambda);
   }
   if(! fHistPVzAnalysisHasHighPtLambda) {
         fHistPVzAnalysisHasHighPtLambda = new TH1F("fHistPVzAnalysisHasHighPtLambda", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHistV0->Add(fHistPVzAnalysisHasHighPtLambda);
   }
   if(! fHistSwappedV0Counter) {
      fHistSwappedV0Counter = new TH1F("fHistSwappedV0Counter", 
         "Swap or not histo;Swapped (1) or not (0); count", 
         2, 0, 2); 		
      fListHistV0->Add(fHistSwappedV0Counter);
   }

   //List of Histograms: Normal
   PostData(1, fListHistV0);

   //TTree Object: Saved to base directory. Should cache to disk while saving. 
   //(Important to avoid excessive memory usage, particularly when merging)
   PostData(2, fTree);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskExtractPerformanceV0::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

   AliESDEvent *lESDevent = 0x0;
   AliMCEvent  *lMCevent  = 0x0; 
   AliStack    *lMCstack  = 0x0; 

   Int_t    lNumberOfV0s                      = -1;
   Double_t lTrkgPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
   Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
   Double_t lMagneticField                 = -10.;
   
  // Connect to the InputEvent   
  // After these lines, we should have an ESD/AOD event + the number of V0s in it.

   // Appropriate for ESD analysis! 
      
   lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
   if (!lESDevent) {
      AliWarning("ERROR: lESDevent not available \n");
      return;
   }
        
   fTreeVariableRunNumber = lESDevent->GetRunNumber();
   fTreeVariableEventNumber =  
    ( ( ((ULong64_t)lESDevent->GetPeriodNumber() ) << 36 ) |
      ( ((ULong64_t)lESDevent->GetOrbitNumber () ) << 12 ) |
        ((ULong64_t)lESDevent->GetBunchCrossNumber() )  );

   lMCevent = MCEvent();
   if (!lMCevent) {
      Printf("ERROR: Could not retrieve MC event \n");
      cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;   
      return;
   }

   lMCstack = lMCevent->Stack();
   if (!lMCstack) {
      Printf("ERROR: Could not retrieve MC stack \n");
      cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
      return;
   }
   TArrayF mcPrimaryVtx;
   AliGenEventHeader* mcHeader=lMCevent->GenEventHeader();
   if(!mcHeader) return;
   mcHeader->PrimaryVertex(mcPrimaryVtx);
        
//------------------------------------------------
// Multiplicity Information Acquistion
//------------------------------------------------

   //REVISED multiplicity estimator after 'multiplicity day' (2011)
    Int_t lMultiplicity = -100;
    Int_t lMultiplicityV0A = -100;
    Int_t lMultiplicityZNA = -100;
    Int_t lMultiplicityTRK = -100;
    Int_t lMultiplicitySPD = -100;
  
   //testing purposes
   if(fkIsNuclear == kFALSE) lMultiplicity = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,0.5);

   //---> If this is a nuclear collision, then go nuclear on "multiplicity" variable...
   //---> Warning: Experimental
   if(fkIsNuclear == kTRUE){ 
      AliCentrality* centrality;
      centrality = lESDevent->GetCentrality();
      lMultiplicity = ( ( Int_t ) ( centrality->GetCentralityPercentile(   fCentralityEstimator.Data() ) ) );
      lMultiplicityV0A = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "V0A" ) ) );
      lMultiplicityZNA = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "ZNA" ) ) );
      lMultiplicityTRK = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "TRK" ) ) );
      lMultiplicitySPD = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "SPD" ) ) );
      if (centrality->GetQuality()>1) {
        PostData(1, fListHistV0);
        PostData(2, fTree);
        return;
      }
   }
  
   //Set variable for filling tree afterwards!
   //---> pp case......: GetReferenceMultiplicity
   //---> Pb-Pb case...: Centrality by V0M
  fTreeVariableMultiplicity = lMultiplicity;
  fTreeVariableMultiplicityV0A = lMultiplicityV0A;
  fTreeVariableMultiplicityZNA = lMultiplicityZNA;
  fTreeVariableMultiplicityTRK = lMultiplicityTRK;
  fTreeVariableMultiplicitySPD = lMultiplicitySPD;

   fHistV0MultiplicityBeforeTrigSel->Fill ( lESDevent->GetNumberOfV0s() );
  
  fHistMultiplicityBeforeTrigSel->Fill ( lMultiplicity );
  fHistMultiplicityV0ABeforeTrigSel->Fill ( lMultiplicityV0A );
  fHistMultiplicityZNABeforeTrigSel->Fill ( lMultiplicityZNA );
  fHistMultiplicityTRKBeforeTrigSel->Fill ( lMultiplicityTRK );
  fHistMultiplicitySPDBeforeTrigSel->Fill ( lMultiplicitySPD );
  
//------------------------------------------------
// MC Information Acquistion
//------------------------------------------------

   Int_t iNumberOfPrimaries = -1;
   iNumberOfPrimaries = lMCstack->GetNprimary();
   if(iNumberOfPrimaries < 1) return; 
   Bool_t lHasHighPtLambda = kFALSE;

//------------------------------------------------
// Variable Definition
//------------------------------------------------

   Int_t lNbMCPrimary        = 0;

   Int_t lPdgcodeCurrentPart = 0;
   Double_t lRapCurrentPart  = 0;
   Double_t lPtCurrentPart   = 0;
  
   //Int_t lComeFromSigma      = 0;

   // current mc particle 's mother
   //Int_t iCurrentMother  = 0;
   lNbMCPrimary = lMCstack->GetNprimary();

//------------------------------------------------
// Pre-Physics Selection
//------------------------------------------------

//----- Loop on primary Xi, Omega --------------------------------------------------------------
   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < lNbMCPrimary; iCurrentLabelStack++) 
   {// This is the begining of the loop on primaries
      
      TParticle* lCurrentParticlePrimary = 0x0; 
      lCurrentParticlePrimary = lMCstack->Particle( iCurrentLabelStack );
      if(!lCurrentParticlePrimary){
         Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
         continue;
      }
      if ( TMath::Abs(lCurrentParticlePrimary->GetPdgCode()) == 3312 || TMath::Abs(lCurrentParticlePrimary->GetPdgCode()) == 3334 ) { 
         Double_t lRapXiMCPrimary = -100;
         if( (lCurrentParticlePrimary->Energy() - lCurrentParticlePrimary->Pz() +1.e-13) != 0 ) { 
           if ( (lCurrentParticlePrimary->Energy() + lCurrentParticlePrimary->Pz()) / (lCurrentParticlePrimary->Energy() - lCurrentParticlePrimary->Pz() +1.e-13) !=0 ){
             lRapXiMCPrimary = 0.5*TMath::Log( (lCurrentParticlePrimary->Energy() + lCurrentParticlePrimary->Pz()) / (lCurrentParticlePrimary->Energy() - lCurrentParticlePrimary->Pz() +1.e-13) );
           }
         }

         //=================================================================================
         // Xi Histograms
         if( lCurrentParticlePrimary->GetPdgCode() == 3312 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenPtVsYVsMultXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenPtVsYCMSVsMultXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
         }
         if( lCurrentParticlePrimary->GetPdgCode() == -3312 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenPtVsYVsMultXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenPtVsYCMSVsMultXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
         }
         // Omega Histograms
         if( lCurrentParticlePrimary->GetPdgCode() == 3334 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenPtVsYVsMultOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenPtVsYCMSVsMultOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
         }
         if( lCurrentParticlePrimary->GetPdgCode() == -3334 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenPtVsYVsMultOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenPtVsYCMSVsMultOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
         }
      } 
   }
//----- End Loop on primary Xi, Omega ----------------------------------------------------------

//--------- GENERATED NUMBER OF CHARGED PARTICLES
// ---> Set Variables to Zero again
// ---> Variable Definition

  Long_t lNumberOfCharged = 0; 

//----- Loop on Stack ----------------------------------------------------------------
   for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++) 
   {// This is the begining of the loop on tracks
      TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack);
      if(!particleOne) continue;
      if(!particleOne->GetPDG()) continue;
      Double_t lThisCharge = particleOne->GetPDG()->Charge()/3.;
      if(TMath::Abs(lThisCharge)<0.001) continue;
      if(! (lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;
     
      //Double_t gpt = particleOne -> Pt();
      Double_t geta = particleOne -> Eta(); 

      if( TMath::Abs(geta) < 0.5) lNumberOfCharged++; 
   }//End of loop on tracks
//----- End Loop on Stack ------------------------------------------------------------

   //Double_t lpArapidityShift = 0.465;
   Bool_t lStackNatural = kTRUE;
//----- Loop on Lambda, K0Short ----------------------------------------------------------------
   for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++) 
   {// This is the begining of the loop on tracks
      
      TParticle* lCurrentParticleForLambdaCheck = 0x0; 
      lCurrentParticleForLambdaCheck = lMCstack->Particle( iCurrentLabelStack );
      if(!lCurrentParticleForLambdaCheck){
         Printf("V0s loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
         continue;
      }

      //=================================================================================
      //Single-Strange checks
      // Keep only K0s, Lambda and AntiLambda:
      lPdgcodeCurrentPart = lCurrentParticleForLambdaCheck->GetPdgCode();	      

      if ( (lCurrentParticleForLambdaCheck->GetPdgCode() == 310   ) ||
           (lCurrentParticleForLambdaCheck->GetPdgCode() == 3122  ) ||
           (lCurrentParticleForLambdaCheck->GetPdgCode() == -3122 ) )
	   {
         lRapCurrentPart   = MyRapidity(lCurrentParticleForLambdaCheck->Energy(),lCurrentParticleForLambdaCheck->Pz());
         lPtCurrentPart    = lCurrentParticleForLambdaCheck->Pt();

          //Use Close to PV for filling CloseToPV histograms!
         Double_t dx, dy, dz; 

         dx = ( (mcPrimaryVtx.At(0)) - (lCurrentParticleForLambdaCheck->Vx()) ); 
         dy = ( (mcPrimaryVtx.At(1)) - (lCurrentParticleForLambdaCheck->Vy()) );
         dz = ( (mcPrimaryVtx.At(2)) - (lCurrentParticleForLambdaCheck->Vz()) );
         Double_t lDistToPV = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
         if( lDistToPV <= 0.001){ 
           if( lPdgcodeCurrentPart == 3122 ){
              f3dHistPrimCloseToPVPtVsYVsMultLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
           }
           if( lPdgcodeCurrentPart == -3122 ){
              f3dHistPrimCloseToPVPtVsYVsMultAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
           }
           if( lPdgcodeCurrentPart == 310 ){
              f3dHistPrimCloseToPVPtVsYVsMultK0Short->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
           }
         }

         //Use Physical Primaries only for filling PrimRaw Histograms!
         if ( lMCstack->IsPhysicalPrimary(iCurrentLabelStack)!=kTRUE ) continue;

          lStackNatural = lMCevent->IsFromBGEvent(iCurrentLabelStack); //Is it? 
          if (!lStackNatural){
            if (!(lCurrentParticleForLambdaCheck->GetFirstMother()<0)) 
              {lStackNatural = kTRUE;} // because there are primaries (ALICE definition) not produced in the collision
          }

         if( lPdgcodeCurrentPart == 3122 ){
            f3dHistPrimRawPtVsYVsMultLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
            f3dHistPrimRawPtVsYCMSVsMultLambda->Fill(lPtCurrentPart, lRapCurrentPart+fpArapidityShift, lMultiplicity);
            if(lStackNatural){f3dHistPrimRawPtVsYVsMultNonInjLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);}
            f3dHistPrimRawPtVsYVsMultMCLambda->Fill(lPtCurrentPart, lRapCurrentPart, lNumberOfCharged);
            f3dHistPrimRawPtVsYVsVertexZLambda->Fill(lPtCurrentPart, lRapCurrentPart, mcPrimaryVtx.At(2));
            if( TMath::Abs( lCurrentParticleForLambdaCheck->Eta() )<1.2 && lPtCurrentPart>2 ){
               lHasHighPtLambda = kTRUE; //Keep track of events with Lambda within |eta|<1.2 and pt>2
            }
         }
         if( lPdgcodeCurrentPart == -3122 ){
            f3dHistPrimRawPtVsYVsMultAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
            f3dHistPrimRawPtVsYCMSVsMultAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart+fpArapidityShift, lMultiplicity);
            if(lStackNatural){f3dHistPrimRawPtVsYVsMultNonInjAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);}
            f3dHistPrimRawPtVsYVsMultMCAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart, lNumberOfCharged);
            f3dHistPrimRawPtVsYVsVertexZAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart, mcPrimaryVtx.At(2));
         }
         if( lPdgcodeCurrentPart == 310 ){
            f3dHistPrimRawPtVsYVsMultK0Short->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
            f3dHistPrimRawPtVsYCMSVsMultK0Short->Fill(lPtCurrentPart, lRapCurrentPart+fpArapidityShift, lMultiplicity);
            if(lStackNatural){f3dHistPrimRawPtVsYVsMultNonInjK0Short->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);}
            f3dHistPrimRawPtVsYVsMultMCK0Short->Fill(lPtCurrentPart, lRapCurrentPart, lNumberOfCharged);
            f3dHistPrimRawPtVsYVsVertexZK0Short->Fill(lPtCurrentPart, lRapCurrentPart, mcPrimaryVtx.At(2));
         }
         //Decay Length Acquisition=====================================================
         Double_t decaylength = -1; 
         Double_t lV0Mass = -1; 
          
         if( !(lCurrentParticleForLambdaCheck->GetDaughter(0) < 0) ) {
            TParticle* lDght0ofV0 = lMCstack->Particle(  lCurrentParticleForLambdaCheck->GetDaughter(0) ); //get first daughter
            if(lDght0ofV0){ // skip if not defined. 
               decaylength = TMath::Sqrt(
				        TMath::Power( lCurrentParticleForLambdaCheck->Vx() - lDght0ofV0->Vx() , 2) + 
				        TMath::Power( lCurrentParticleForLambdaCheck->Vy() - lDght0ofV0->Vy() , 2) + 
				        TMath::Power( lCurrentParticleForLambdaCheck->Vz() - lDght0ofV0->Vz() , 2)
               );
               //Need to correct for relativitity! Involves multiplying by mass and dividing by momentum. 
               if(TMath::Abs( lPdgcodeCurrentPart ) == 3122 ) { lV0Mass = 1.115683; }
               if(TMath::Abs( lPdgcodeCurrentPart ) == 310 ) { lV0Mass = 0.497614; }
               if( lCurrentParticleForLambdaCheck->P() + 1e-10 != 0 ) decaylength = ( lV0Mass * decaylength ) / ( lCurrentParticleForLambdaCheck->P() + 1e-10 );
               if( lCurrentParticleForLambdaCheck->P() + 1e-10 == 0 ) decaylength = 1e+5;
            }
         }
         if( lPdgcodeCurrentPart == 3122) f3dHistPrimRawPtVsYVsDecayLengthLambda ->Fill( lPtCurrentPart, lRapCurrentPart , decaylength ); 
         if( lPdgcodeCurrentPart == -3122) f3dHistPrimRawPtVsYVsDecayLengthAntiLambda ->Fill( lPtCurrentPart, lRapCurrentPart , decaylength ); 
         if( lPdgcodeCurrentPart == 310) f3dHistPrimRawPtVsYVsDecayLengthK0Short ->Fill( lPtCurrentPart, lRapCurrentPart , decaylength ); 
      }
   }//End of loop on tracks
//----- End Loop on Lambda, K0Short ------------------------------------------------------------


   f2dHistMultiplicityVsTrueBeforeTrigSel->Fill ( lMultiplicity , lNumberOfCharged );

    fTreeVariableMultiplicityMC = lNumberOfCharged;

   fHistGenVertexZBeforeTrigSel->Fill( (mcPrimaryVtx.At(2)) );

   lPdgcodeCurrentPart = 0;
   lRapCurrentPart  = 0;
   lPtCurrentPart   = 0;

//------------------------------------------------
// Physics Selection
//------------------------------------------------
  
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  Bool_t isSelectedExtra = kTRUE; //extra sel, default YES
  isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
  
  //pA triggering: CINT7
  if( fkSwitchINT7 ) isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
  
  //Extra selection applies if with/without SDD is to be dealth with
  if( fkFastOnly == "kFastOnly"){
    //If not kFastOnly, isSelectedExtra will be kFALSE; procedure will reject it
    isSelectedExtra = (maskIsSelected & AliVEvent::kFastOnly) == AliVEvent::kFastOnly;
  }
  if( fkFastOnly == "NotkFastOnly"){
    //If not kFastOnly, isSelectedExtra will be kTRUE; procedure will accept it
    isSelectedExtra = !( (maskIsSelected & AliVEvent::kFastOnly) == AliVEvent::kFastOnly );
  }
  
  //Standard Min-Bias Selection
  if ( ! isSelected ) {
    PostData(1, fListHistV0);
    PostData(2, fTree);
    return;
  }
  //Check if goes through extra selections
  //isSelectedExtra will be true in case -> fkFastOnly==""
  //isSelectedExtra will be true in case -> fkFastOnly=="kFastOnly"    && bit kFastOnly ON
  //isSelectedExtra will be true in case -> fkFastOnly=="NotkFastOnly" && bit kFastOnly OFF
  if ( !isSelectedExtra ) {
    PostData(1, fListHistV0);
    PostData(2, fTree);
    return;
  }
  
  f2dHistMultiplicityVsTrueForTrigEvt->Fill ( lMultiplicity , lNumberOfCharged );
  fHistGenVertexZForTrigEvt->Fill( mcPrimaryVtx.At(2) );
  
//------------------------------------------------
// After Trigger Selection
//------------------------------------------------

   lNumberOfV0s          = lESDevent->GetNumberOfV0s();
  
   //Set variable for filling tree afterwards!
   fHistV0MultiplicityForTrigEvt->Fill(lNumberOfV0s);
   fHistMultiplicityForTrigEvt->Fill ( lMultiplicity );
  fHistMultiplicityV0AForTrigEvt       ->Fill( lMultiplicityV0A  );
  fHistMultiplicityZNAForTrigEvt       ->Fill( lMultiplicityZNA  );
  fHistMultiplicityTRKForTrigEvt       ->Fill( lMultiplicityTRK  );
  fHistMultiplicitySPDForTrigEvt       ->Fill( lMultiplicitySPD  );

//------------------------------------------------
// Getting: Primary Vertex + MagField Info
//------------------------------------------------

   const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
   // get the vtx stored in ESD found with tracks
   lPrimaryTrackingESDVtx->GetXYZ( lTrkgPrimaryVtxPos );
        
   const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();	
   // get the best primary vertex available for the event
   // As done in AliCascadeVertexer, we keep the one which is the best one available.
   // between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
   // This one will be used for next calculations (DCA essentially)
   lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

   Double_t lPrimaryVtxPosition[3];
   const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
   lPrimaryVtxPosition[0] = primaryVtx->GetX();
   lPrimaryVtxPosition[1] = primaryVtx->GetY();
   lPrimaryVtxPosition[2] = primaryVtx->GetZ();
   fHistPVx->Fill( lPrimaryVtxPosition[0] );
   fHistPVy->Fill( lPrimaryVtxPosition[1] );
   fHistPVz->Fill( lPrimaryVtxPosition[2] );

   f2dHistMultiplicityVsVertexZForTrigEvt->Fill( lMultiplicity, lPrimaryVtxPosition[2] );

//------------------------------------------------
// Primary Vertex Z position: SKIP
//------------------------------------------------

   if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0 ) { 
      AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !"); 
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return; 
   }

   f2dHistMultiplicityVsVertexZ->Fill( lMultiplicity, lPrimaryVtxPosition[2] );

   lMagneticField = lESDevent->GetMagneticField( );
   fHistV0MultiplicityForSelEvt ->Fill( lNumberOfV0s );
   fHistMultiplicity->Fill(lMultiplicity);
   fHistMultiplicityV0A->Fill(lMultiplicityV0A);
   fHistMultiplicityZNA->Fill(lMultiplicityZNA);
   fHistMultiplicityTRK->Fill(lMultiplicityTRK);
   fHistMultiplicitySPD->Fill(lMultiplicitySPD);
   f2dHistMultiplicityVsTrue->Fill ( lMultiplicity , lNumberOfCharged );
   fHistGenVertexZ->Fill( (mcPrimaryVtx.At(2)) );
//------------------------------------------------
// SKIP: Events with well-established PVtx
//------------------------------------------------
	
   const AliESDVertex *lPrimaryTrackingESDVtxCheck = lESDevent->GetPrimaryVertexTracks();
   const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
   if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtxCheck->GetStatus() ){
      AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return;
   }

   f2dHistMultiplicityVsVertexZNoTPCOnly->Fill( lMultiplicity, lPrimaryVtxPosition[2] );
   fHistV0MultiplicityForSelEvtNoTPCOnly ->Fill( lNumberOfV0s );
   fHistMultiplicityNoTPCOnly->Fill(lMultiplicity);
   fHistMultiplicityV0ANoTPCOnly->Fill(lMultiplicityV0A);
   fHistMultiplicityZNANoTPCOnly->Fill(lMultiplicityZNA);
   fHistMultiplicityTRKNoTPCOnly->Fill(lMultiplicityTRK);
   fHistMultiplicitySPDNoTPCOnly->Fill(lMultiplicitySPD);
   f2dHistMultiplicityVsTrueNoTPCOnly->Fill ( lMultiplicity , lNumberOfCharged );
   fHistGenVertexZNoTPCOnly->Fill( (mcPrimaryVtx.At(2)) );
//------------------------------------------------
// Pileup Rejection Studies
//------------------------------------------------

   // FIXME : quality selection regarding pile-up rejection 
   if(lESDevent->IsPileupFromSPD() && !fkIsNuclear ){// minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  -> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
      AliWarning("Pb / This is tagged as Pileup from SPD... return !");
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return;
   }
   f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup->Fill( lMultiplicity, lPrimaryVtxPosition[2] );
   fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup ->Fill( lNumberOfV0s );
 
   fHistMultiplicityNoTPCOnlyNoPileup->Fill(lMultiplicity);
   fHistMultiplicityV0ANoTPCOnlyNoPileup->Fill(lMultiplicityV0A);
   fHistMultiplicityZNANoTPCOnlyNoPileup->Fill(lMultiplicityZNA);
   fHistMultiplicityTRKNoTPCOnlyNoPileup->Fill(lMultiplicityTRK);
   fHistMultiplicitySPDNoTPCOnlyNoPileup->Fill(lMultiplicitySPD);
  
   f2dHistMultiplicityVsTrueNoTPCOnlyNoPileup->Fill ( lMultiplicity , lNumberOfCharged );
   fHistGenVertexZNoTPCOnlyNoPileup->Fill( (mcPrimaryVtx.At(2)) );
   //Do control histograms without the IsFromVertexerZ events, but consider them in analysis...
   if( ! (lESDevent->GetPrimaryVertex()->IsFromVertexerZ() )	 ){ 
      fHistPVxAnalysis->Fill( lPrimaryVtxPosition[0] );
      fHistPVyAnalysis->Fill( lPrimaryVtxPosition[1] );
      fHistPVzAnalysis->Fill( lPrimaryVtxPosition[2] );
      if ( lHasHighPtLambda == kTRUE ){ 
         fHistPVxAnalysisHasHighPtLambda->Fill( lPrimaryVtxPosition[0] );
         fHistPVyAnalysisHasHighPtLambda->Fill( lPrimaryVtxPosition[1] );
         fHistPVzAnalysisHasHighPtLambda->Fill( lPrimaryVtxPosition[2] );
      }
   }

  fTreeVariableVertexZ = lPrimaryVtxPosition[2];

  fTreeVariablePVx = lPrimaryVtxPosition[0];
  fTreeVariablePVy = lPrimaryVtxPosition[1];
  fTreeVariablePVz = lPrimaryVtxPosition[2];

  fTreeVariableMCPVx = (mcPrimaryVtx.At(0));
  fTreeVariableMCPVy = (mcPrimaryVtx.At(1));
  fTreeVariableMCPVz = (mcPrimaryVtx.At(2));

//------------------------------------------------
// stack loop starts here
//------------------------------------------------

//---> Loop over ALL PARTICLES
 
   for (Int_t iMc = 0; iMc < (lMCstack->GetNtrack()); iMc++) {  
      TParticle *p0 = lMCstack->Particle(iMc); 
      if (!p0) {
         //Printf("ERROR: particle with label %d not found in lMCstack (mc loop)", iMc);
         continue;
      }
      lPdgcodeCurrentPart = p0->GetPdgCode();

      // Keep only K0s, Lambda and AntiLambda:
      if ( (lPdgcodeCurrentPart != 310 ) && (lPdgcodeCurrentPart != 3122 ) && (lPdgcodeCurrentPart != -3122 ) ) continue;
	
      lRapCurrentPart   = MyRapidity(p0->Energy(),p0->Pz());
      lPtCurrentPart    = p0->Pt();

        //Use Physical Primaries only for filling PrimRaw Histograms!
      if ( lMCstack->IsPhysicalPrimary(iMc)!=kTRUE ) continue;

      if( lPdgcodeCurrentPart == 3122 ){
         f3dHistPrimAnalysisPtVsYVsMultLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
         f3dHistPrimAnalysisPtVsYCMSVsMultLambda->Fill(lPtCurrentPart, lRapCurrentPart+fpArapidityShift, lMultiplicity);
      }
      if( lPdgcodeCurrentPart == -3122 ){
         f3dHistPrimAnalysisPtVsYVsMultAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
         f3dHistPrimAnalysisPtVsYCMSVsMultAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart+fpArapidityShift, lMultiplicity);
      }
      if( lPdgcodeCurrentPart == 310 ){
         f3dHistPrimAnalysisPtVsYVsMultK0Short->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
         f3dHistPrimAnalysisPtVsYCMSVsMultK0Short->Fill(lPtCurrentPart, lRapCurrentPart+fpArapidityShift, lMultiplicity);
      }
   }

//----- Loop on primary Xi, Omega --------------------------------------------------------------
   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < lNbMCPrimary; iCurrentLabelStack++) 
   {// This is the begining of the loop on primaries
      
      TParticle* lCurrentParticlePrimary = 0x0; 
      lCurrentParticlePrimary = lMCstack->Particle( iCurrentLabelStack );
      if(!lCurrentParticlePrimary){
         Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
         continue;
      }
      if ( TMath::Abs(lCurrentParticlePrimary->GetPdgCode()) == 3312 || TMath::Abs(lCurrentParticlePrimary->GetPdgCode()) == 3334 ) { 
         Double_t lRapXiMCPrimary = -100;
         if( (lCurrentParticlePrimary->Energy() - lCurrentParticlePrimary->Pz() +1.e-13) != 0 ) { 
           if ( (lCurrentParticlePrimary->Energy() + lCurrentParticlePrimary->Pz()) / (lCurrentParticlePrimary->Energy() - lCurrentParticlePrimary->Pz() +1.e-13) !=0 ){
             lRapXiMCPrimary = 0.5*TMath::Log( (lCurrentParticlePrimary->Energy() + lCurrentParticlePrimary->Pz()) / (lCurrentParticlePrimary->Energy() - lCurrentParticlePrimary->Pz() +1.e-13) );
           }
         }

         //=================================================================================
         // Xi Histograms
         if( lCurrentParticlePrimary->GetPdgCode() == 3312 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenSelectedPtVsYVsMultXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenSelectedPtVsYCMSVsMultXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
         }
         if( lCurrentParticlePrimary->GetPdgCode() == -3312 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenSelectedPtVsYVsMultXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenSelectedPtVsYCMSVsMultXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
         }
         // Omega Histograms
         if( lCurrentParticlePrimary->GetPdgCode() == 3334 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenSelectedPtVsYVsMultOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
         }
         if( lCurrentParticlePrimary->GetPdgCode() == -3334 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenSelectedPtVsYVsMultOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
         }
      } 
   }
//----- End Loop on primary Xi, Omega ----------------------------------------------------------

//------------------------------------------------
// MAIN LAMBDA LOOP STARTS HERE
//------------------------------------------------

   //Variable definition
   Int_t    lOnFlyStatus = 0;
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
   
   for (Int_t iV0 = 0; iV0 < nv0s; iV0++) 
	{// This is the begining of the V0 loop
      AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
      if (!v0) continue;

      //---> Fix On-the-Fly candidates
      if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ){
        fHistSwappedV0Counter -> Fill( 1 );
      }else{
        fHistSwappedV0Counter -> Fill( 0 ); 
      }
      if ( fkUseOnTheFly ) CheckChargeV0(v0); 


      Double_t tV0mom[3];
      v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] ); 
      Double_t lV0TotalMomentum = TMath::Sqrt(
         tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );

      Double_t tDecayVertexV0[3]; v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]); 
      lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);
      lPt = v0->Pt();
      lRapK0Short = v0->RapK0Short();
      lRapLambda  = v0->RapLambda();

      //Set Variables for later filling
      fTreeVariableV0x = tDecayVertexV0[0];
      fTreeVariableV0y = tDecayVertexV0[1];
      fTreeVariableV0z = tDecayVertexV0[2];

      //Set Variables for later filling
      fTreeVariableV0Px = tV0mom[0];
      fTreeVariableV0Py = tV0mom[1];
      fTreeVariableV0Pz = tV0mom[2];

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

      //Get status flags
      fTreeVariablePosTrackStatus = pTrack->GetStatus();
      fTreeVariableNegTrackStatus = nTrack->GetStatus();
    
      if ( ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) )&&(fkTakeAllTracks==kFALSE) ) continue;
	
      //GetKinkIndex condition
      if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;

      //Findable clusters > 0 condition
      if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;

      //Compute ratio Crossed Rows / Findable clusters
      //Note: above test avoids division by zero! 
      Float_t lPosTrackCrossedRowsOverFindable = -1;
      Float_t lNegTrackCrossedRowsOverFindable = -1;
      if ( ((double)(pTrack->GetTPCNclsF()) ) != 0 ) lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF())); 
      if ( ((double)(nTrack->GetTPCNclsF()) ) != 0 ) lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF())); 

      fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
      if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable )
         fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;

      //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
      if ( (fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8)&&(fkTakeAllTracks==kFALSE) ) continue;

      //End track Quality Cuts
      //________________________________________________________________________

      lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lPrimaryVtxPosition[0],
							lPrimaryVtxPosition[1],
							lMagneticField) );

      lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lPrimaryVtxPosition[0],
							lPrimaryVtxPosition[1],
							lMagneticField) );

      lOnFlyStatus = v0->GetOnFlyStatus();
      lChi2V0 = v0->GetChi2V0();
      lDcaV0Daughters = v0->GetDcaV0Daughters();
      lDcaV0ToPrimVertex = v0->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
      lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
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

      //fTreeVariableOnFlyStatus = lOnFlyStatus;
      //fHistV0OnFlyStatus->Fill(lOnFlyStatus);

//===============================================
// Monte Carlo Association starts here
//===============================================

      //---> Set Everything to "I don't know" before starting

      fTreeVariablePIDPositive = 0;
      fTreeVariablePIDNegative = 0;

      fTreeVariableIndexStatus = 0;
      fTreeVariableIndexStatusMother = 0;

      fTreeVariablePtMother = -1;
      fTreeVariablePtMC = -1;
      fTreeVariableRapMC = -100;

      fTreeVariablePID = -1; 
      fTreeVariablePIDMother = -1;

      fTreeVariablePrimaryStatus = 0; 
      fTreeVariablePrimaryStatusMother = 0; 
      fTreeVariableV0CreationRadius = -1;

      Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrack->GetLabel() );  
      Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrack->GetLabel() );
		
      TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
      TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
	    
      fTreeVariablePosTransvMomentumMC = mcPosV0Dghter->Pt();
      fTreeVariableNegTransvMomentumMC = mcNegV0Dghter->Pt();

      Int_t lPIDPositive = mcPosV0Dghter -> GetPdgCode();
      Int_t lPIDNegative = mcNegV0Dghter -> GetPdgCode();

      fTreeVariablePIDPositive = lPIDPositive;
      fTreeVariablePIDNegative = lPIDNegative;

      Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ; 
      Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();

      if( lblMotherPosV0Dghter == lblMotherNegV0Dghter && lblMotherPosV0Dghter > -1 ){
         //either label is fine, they're equal at this stage
         TParticle* pThisV0 = lMCstack->Particle( lblMotherPosV0Dghter ); 
         //Set tree variables
         fTreeVariablePID   = pThisV0->GetPdgCode(); //PDG Code
         fTreeVariablePtMC  = pThisV0->Pt(); //Perfect Pt

          fTreeVariableIsNonInjected = lMCevent->IsFromBGEvent(lblMotherPosV0Dghter); //Is it? 
          if (!fTreeVariableIsNonInjected){
            if (!(pThisV0->GetFirstMother()<0)) 
              {fTreeVariableIsNonInjected = kTRUE;} // because there are primaries (ALICE definition) not produced in the collision
          }

         //Set Variables for later filling
         //Be careful: Vx, Vy, Vz: Creation vertex. So decay position is the 
         //Creation vertex of any one of the daughters!
         fTreeVariableMCV0x = mcPosV0Dghter->Vx();
         fTreeVariableMCV0y = mcPosV0Dghter->Vy();
         fTreeVariableMCV0z = mcPosV0Dghter->Vz();

         //Set Variables for later filling
         fTreeVariableMCV0Px = pThisV0->Px();
         fTreeVariableMCV0Py = pThisV0->Py();
         fTreeVariableMCV0Pz = pThisV0->Pz();

         //Only Interested if it's a Lambda, AntiLambda or K0s 
         //Avoid the Junction Bug! PYTHIA has particles with Px=Py=Pz=E=0 occasionally, 
         //having particle code 88 (unrecognized by PDG), for documentation purposes.
         //Even ROOT's TParticle::Y() is not prepared to deal with that exception!
         //Note that TParticle::Pt() is immune (that would just return 0)...
         //Though granted that that should be extremely rare in this precise condition...
         if( TMath::Abs(fTreeVariablePID) == 3122 || fTreeVariablePID==310 ){
            fTreeVariableRapMC = pThisV0->Y(); //Perfect Y
         }
         fTreeVariableV0CreationRadius = TMath::Sqrt(
          TMath::Power(  ( (mcPrimaryVtx.At(0)) - (pThisV0->Vx()) ) , 2) + 
          TMath::Power(  ( (mcPrimaryVtx.At(1)) - (pThisV0->Vy()) ) , 2) + 
          TMath::Power(  ( (mcPrimaryVtx.At(2)) - (pThisV0->Vz()) ) , 2) 
         );
         if( lblMotherPosV0Dghter  < lNbMCPrimary ) fTreeVariableIndexStatus = 1; //looks primary
         if( lblMotherPosV0Dghter >= lNbMCPrimary ) fTreeVariableIndexStatus = 2; //looks secondary
         if( lMCstack->IsPhysicalPrimary       (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 1; //Is Primary!
         if( lMCstack->IsSecondaryFromWeakDecay(lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 2; //Weak Decay!
         if( lMCstack->IsSecondaryFromMaterial (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 3; //Material Int!
         
         //Now we try to acquire the V0 parent particle, if possible
         Int_t lblThisV0Parent = pThisV0->GetFirstMother();
         if ( lblThisV0Parent > -1 ){ //if it has a parent, get it and store specs
            TParticle* pThisV0Parent = lMCstack->Particle( lblThisV0Parent );
            fTreeVariablePIDMother   = pThisV0Parent->GetPdgCode(); //V0 Mother PDG
            fTreeVariablePtMother    = pThisV0Parent->Pt();         //V0 Mother Pt
            //Primary Status for the V0 Mother particle 
            if( lblThisV0Parent  < lNbMCPrimary ) fTreeVariableIndexStatusMother = 1; //looks primary
            if( lblThisV0Parent >= lNbMCPrimary ) fTreeVariableIndexStatusMother = 2; //looks secondary
            if( lMCstack->IsPhysicalPrimary       (lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 1; //Is Primary!
            if( lMCstack->IsSecondaryFromWeakDecay(lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 2; //Weak Decay!
            if( lMCstack->IsSecondaryFromMaterial (lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 3; //Material Int!
         }
      }

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

//tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]
      Double_t lDistanceTravelled = TMath::Sqrt(
						TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
					);
      fTreeVariableDistOverTotMom = 1e+5;
      if( lV0TotalMomentum + 1e-10 != 0 ) fTreeVariableDistOverTotMom = lDistanceTravelled / (lV0TotalMomentum + 1e-10); //avoid division by zero, to be sure

      Double_t lMomentumPosTemp[3];
      pTrack->GetPxPyPz(lMomentumPosTemp);
      Double_t lPtPosTemporary = sqrt(pow(lMomentumPosTemp[0],2) + pow(lMomentumPosTemp[1],2));

      Double_t lMomentumNegTemp[3];
      nTrack->GetPxPyPz(lMomentumNegTemp);
      Double_t lPtNegTemporary = sqrt(pow(lMomentumNegTemp[0],2) + pow(lMomentumNegTemp[1],2));

      fTreeVariablePosTransvMomentum = lPtPosTemporary;
      fTreeVariableNegTransvMomentum = lPtNegTemporary;


//------------------------------------------------
// Fill Tree! 
//------------------------------------------------

      // The conditionals are meant to decrease excessive
      // memory usage! 

      //Modified version: Keep only OnFlyStatus == 0
      //Keep only if included in a parametric InvMass Region 20 sigmas away from peak

      //First Selection: Reject OnFly
      if( (lOnFlyStatus == 0 && fkUseOnTheFly == kFALSE) || (lOnFlyStatus != 0 && fkUseOnTheFly == kTRUE ) ){
         //Second Selection: rough 20-sigma band, parametric. 
         //K0Short: Enough to parametrize peak broadening with linear function.    
         Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeVariablePt; 
         Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeVariablePt;
         //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
         //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
         Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeVariablePt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeVariablePt); 
         Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeVariablePt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeVariablePt);
         //Do Selection      
         if( (fTreeVariableInvMassLambda     < lUpperLimitLambda  && fTreeVariableInvMassLambda     > lLowerLimitLambda     ) || 
             (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda     ) || 
             (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short    ) ){
             //Pre-selection in case this is AA...
             if( fkIsNuclear == kFALSE ) fTree->Fill();
             if( fkIsNuclear == kTRUE){ 
             //If this is a nuclear collision___________________
             // ... pre-filter with daughter eta selection only (not TPC)
               if ( TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8 ) fTree->Fill();
             }//end nuclear_____________________________________
         }
      }

//------------------------------------------------
// Fill tree over.
//------------------------------------------------


   }// This is the end of the V0 loop

//------------------------------------------------

   // Post output data.
   PostData(1, fListHistV0);
   PostData(2, fTree);
}

//________________________________________________________________________
void AliAnalysisTaskExtractPerformanceV0::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   TList *cRetrievedList = 0x0;
   cRetrievedList = (TList*)GetOutputData(1);
   if(!cRetrievedList){
      Printf("ERROR - AliAnalysisTaskExtractV0 : ouput data container list not available\n");
      return;
   }	
	
   fHistV0MultiplicityForTrigEvt = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistV0MultiplicityForTrigEvt")  );
   if (!fHistV0MultiplicityForTrigEvt) {
      Printf("ERROR - AliAnalysisTaskExtractV0 : fHistV0MultiplicityForTrigEvt not available");
      return;
   }
  
   TCanvas *canCheck = new TCanvas("AliAnalysisTaskExtractV0","V0 Multiplicity",10,10,510,510);
   canCheck->cd(1)->SetLogy();

   fHistV0MultiplicityForTrigEvt->SetMarkerStyle(22);
   fHistV0MultiplicityForTrigEvt->DrawCopy("E");
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskExtractPerformanceV0::MyRapidity(Double_t rE, Double_t rPz) const
{
   // Local calculation for rapidity
   Double_t ReturnValue = -100;
   if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ){ 
      ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
   }
   return ReturnValue;
} 

//________________________________________________________________________
void AliAnalysisTaskExtractPerformanceV0::CheckChargeV0(AliESDv0 *v0)
{
   // This function checks charge of negative and positive daughter tracks. 
   // If incorrectly defined (onfly vertexer), swaps out. 
   if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ){
      //V0 daughter track swapping is required! Note: everything is swapped here... P->N, N->P
      Long_t lCorrectNidx = v0->GetPindex();
      Long_t lCorrectPidx = v0->GetNindex();
      Double32_t	lCorrectNmom[3];
      Double32_t	lCorrectPmom[3];
      v0->GetPPxPyPz( lCorrectNmom[0], lCorrectNmom[1], lCorrectNmom[2] );
      v0->GetNPxPyPz( lCorrectPmom[0], lCorrectPmom[1], lCorrectPmom[2] );

      AliExternalTrackParam	lCorrectParamN(
        v0->GetParamP()->GetX() , 
        v0->GetParamP()->GetAlpha() , 
        v0->GetParamP()->GetParameter() , 
        v0->GetParamP()->GetCovariance() 
      );
      AliExternalTrackParam	lCorrectParamP(
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
      }else{
        //AliWarning("Found Swapped Charges and fixed.");
      }
      //________________________________________________________________
   }else{
      //Don't touch it! ---
      //Printf("Ah, nice. Charges are already ordered...");
   }
   return;
} 

