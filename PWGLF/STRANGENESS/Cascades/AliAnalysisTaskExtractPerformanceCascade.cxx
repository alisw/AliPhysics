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
//   2. Loop over stack to find Cascades, fill TH3Fs "PrimRawPt"s for Efficiency
//   3. Perform Physics Selection
//   4. Perform Primary Vertex |z|<10cm selection
//   5. Perform Primary Vertex NoTPCOnly vertexing selection (>0 contrib.)
//   6. Perform Pileup Rejection
//   7. Analysis Loops: 
//    7a. Fill TH3Fs "PrimAnalysisPt" for control purposes only
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
#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"

#include "AliAnalysisTaskExtractPerformanceCascade.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskExtractPerformanceCascade)

AliAnalysisTaskExtractPerformanceCascade::AliAnalysisTaskExtractPerformanceCascade() 
  : AliAnalysisTaskSE(), fListHist(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fUtils(0),
   fkIsNuclear   ( kFALSE ), 
   fkSwitchINT7  ( kFALSE ),
   fpArapidityShift ( 0.465 ),
   fCentralityEstimator("V0M"),
   fkpAVertexSelection( kFALSE ),
   fEtaRefMult ( 0.5 ),
   fkRunVertexers ( kFALSE ),
   fkCheckSwapping( kFALSE ),
   fkSelectPeripheral ( kFALSE ),
//------------------------------------------------
// Tree Variables
//------------------------------------------------

   fTreeCascVarCharge(0),
   fTreeCascVarMassAsXi(0),
   fTreeCascVarMassAsOmega(0),
   fTreeCascVarPt(0),
   fTreeCascVarPtMC(0),
   fTreeCascVarRapMC(0),
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
   fTreeCascVarMultiplicity(0),
   fTreeCascVarMultiplicityV0A(0),
   fTreeCascVarMultiplicityZNA(0),
   fTreeCascVarMultiplicityTRK(0),
   fTreeCascVarMultiplicitySPD(0),
   fTreeCascVarMultiplicityMC(0),
   fTreeCascVarDistOverTotMom(0),
   fTreeCascVarIsPhysicalPrimary(0),
   fTreeCascVarPID(0),
   fTreeCascVarPIDSwapped(0),
   fTreeCascVarPIDBachelor(0),
   fTreeCascVarPIDNegative(0),
   fTreeCascVarPIDPositive(0),
   fTreeCascVarBachTransMom(0),
   fTreeCascVarPosTransMom(0),
   fTreeCascVarNegTransMom(0),
   fTreeCascVarPosTransMomMC(0),
   fTreeCascVarNegTransMomMC(0),
   fTreeCascVarNegNSigmaPion(0),
   fTreeCascVarNegNSigmaProton(0),
   fTreeCascVarPosNSigmaPion(0),
   fTreeCascVarPosNSigmaProton(0),
   fTreeCascVarBachNSigmaPion(0),
   fTreeCascVarBachNSigmaKaon(0),

   fTreeCascVarkITSRefitBachelor(0),
   fTreeCascVarkITSRefitNegative(0),
   fTreeCascVarkITSRefitPositive(0),

    fTreeCascVarNegClusters(0),
    fTreeCascVarPosClusters(0),
    fTreeCascVarBachClusters(0),
    fTreeCascVarNegSharedClusters(0),
    fTreeCascVarPosSharedClusters(0),
    fTreeCascVarBachSharedClusters(0),

   fTreeCascVarEvHasXiMinus(0),
   fTreeCascVarEvHasXiPlus(0),
   fTreeCascVarEvHasOmegaMinus(0),
   fTreeCascVarEvHasOmegaPlus(0),
   fTreeCascVarEvHasLambda(0),
   fTreeCascVarEvHasAntiLambda(0),

   fTreeCascVarEvHasLowPtXiMinus(0),
   fTreeCascVarEvHasLowPtXiPlus(0),
   fTreeCascVarEvHasLowPtOmegaMinus(0),
   fTreeCascVarEvHasLowPtOmegaPlus(0),
   fTreeCascVarEvHasLowPtLambda(0),
   fTreeCascVarEvHasLowPtAntiLambda(0),

   fTreeCascVarEvHasVeryLowPtXiMinus(0),
   fTreeCascVarEvHasVeryLowPtXiPlus(0),
   fTreeCascVarEvHasVeryLowPtOmegaMinus(0),
   fTreeCascVarEvHasVeryLowPtOmegaPlus(0),
   fTreeCascVarEvHasVeryLowPtLambda(0),
   fTreeCascVarEvHasVeryLowPtAntiLambda(0),

    fTreeCascVarBadCascadeJai(0),
    fTreeCascVarDeltaDCA(0),

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


//------------------------------------------------
// PARTICLE HISTOGRAMS
// --- Filled on a Particle-by-Particle basis
//------------------------------------------------
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

//MC Mult

f3dHistGenPtVsYVsMultMCXiMinus(0),
f3dHistGenPtVsYVsMultMCXiPlus(0),
f3dHistGenPtVsYVsMultMCOmegaMinus(0),
f3dHistGenPtVsYVsMultMCOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultMCXiMinus(0),
f3dHistGenSelectedPtVsYVsMultMCXiPlus(0),
f3dHistGenSelectedPtVsYVsMultMCOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultMCOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultMCXiMinus(0),
f3dHistGenPtVsYCMSVsMultMCXiPlus(0),
f3dHistGenPtVsYCMSVsMultMCOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultMCOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultMCXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultMCXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultMCOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultMCOmegaPlus(0),

//V0A

f3dHistGenPtVsYVsMultV0AXiMinus(0),
f3dHistGenPtVsYVsMultV0AXiPlus(0),
f3dHistGenPtVsYVsMultV0AOmegaMinus(0),
f3dHistGenPtVsYVsMultV0AOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultV0AXiMinus(0),
f3dHistGenSelectedPtVsYVsMultV0AXiPlus(0),
f3dHistGenSelectedPtVsYVsMultV0AOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultV0AOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultV0AXiMinus(0),
f3dHistGenPtVsYCMSVsMultV0AXiPlus(0),
f3dHistGenPtVsYCMSVsMultV0AOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultV0AOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultV0AXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultV0AXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaPlus(0),
f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiMinus(0),
f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiPlus(0),
f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaMinus(0),
f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaPlus(0),

//ZNA

f3dHistGenPtVsYVsMultZNAXiMinus(0),
f3dHistGenPtVsYVsMultZNAXiPlus(0),
f3dHistGenPtVsYVsMultZNAOmegaMinus(0),
f3dHistGenPtVsYVsMultZNAOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultZNAXiMinus(0),
f3dHistGenSelectedPtVsYVsMultZNAXiPlus(0),
f3dHistGenSelectedPtVsYVsMultZNAOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultZNAOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultZNAXiMinus(0),
f3dHistGenPtVsYCMSVsMultZNAXiPlus(0),
f3dHistGenPtVsYCMSVsMultZNAOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultZNAOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultZNAXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultZNAXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaPlus(0),

//TRK

f3dHistGenPtVsYVsMultTRKXiMinus(0),
f3dHistGenPtVsYVsMultTRKXiPlus(0),
f3dHistGenPtVsYVsMultTRKOmegaMinus(0),
f3dHistGenPtVsYVsMultTRKOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultTRKXiMinus(0),
f3dHistGenSelectedPtVsYVsMultTRKXiPlus(0),
f3dHistGenSelectedPtVsYVsMultTRKOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultTRKOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultTRKXiMinus(0),
f3dHistGenPtVsYCMSVsMultTRKXiPlus(0),
f3dHistGenPtVsYCMSVsMultTRKOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultTRKOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultTRKXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultTRKXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaPlus(0),

//SPD

f3dHistGenPtVsYVsMultSPDXiMinus(0),
f3dHistGenPtVsYVsMultSPDXiPlus(0),
f3dHistGenPtVsYVsMultSPDOmegaMinus(0),
f3dHistGenPtVsYVsMultSPDOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultSPDXiMinus(0),
f3dHistGenSelectedPtVsYVsMultSPDXiPlus(0),
f3dHistGenSelectedPtVsYVsMultSPDOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultSPDOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultSPDXiMinus(0),
f3dHistGenPtVsYCMSVsMultSPDXiPlus(0),
f3dHistGenPtVsYCMSVsMultSPDOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultSPDOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultSPDXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultSPDXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaPlus(0),


   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0)
{
  // Dummy Constructor
}

AliAnalysisTaskExtractPerformanceCascade::AliAnalysisTaskExtractPerformanceCascade(const char *name) 
  : AliAnalysisTaskSE(name), fListHist(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fUtils(0),
   fkIsNuclear   ( kFALSE ), 
   fkSwitchINT7  ( kFALSE ),
   fpArapidityShift ( 0.465 ),
   fCentralityEstimator("V0M"),
   fkpAVertexSelection( kFALSE ),
   fEtaRefMult ( 0.5 ),
   fkRunVertexers ( kFALSE ),
   fkCheckSwapping( kFALSE ),
   fkSelectPeripheral ( kFALSE ),
//------------------------------------------------
// Tree Variables
//------------------------------------------------

   fTreeCascVarCharge(0),
   fTreeCascVarMassAsXi(0),
   fTreeCascVarMassAsOmega(0),
   fTreeCascVarPt(0),
   fTreeCascVarPtMC(0),
   fTreeCascVarRapMC(0),
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
   fTreeCascVarMultiplicity(0),
   fTreeCascVarMultiplicityV0A(0),
   fTreeCascVarMultiplicityZNA(0),
   fTreeCascVarMultiplicityTRK(0),
   fTreeCascVarMultiplicitySPD(0),
   fTreeCascVarMultiplicityMC(0),
   fTreeCascVarDistOverTotMom(0),
   fTreeCascVarIsPhysicalPrimary(0),
   fTreeCascVarPID(0),
   fTreeCascVarPIDSwapped(0),
   fTreeCascVarPIDBachelor(0),
   fTreeCascVarPIDNegative(0),
   fTreeCascVarPIDPositive(0),
   fTreeCascVarBachTransMom(0),
   fTreeCascVarPosTransMom(0),
   fTreeCascVarNegTransMom(0),
   fTreeCascVarPosTransMomMC(0),
   fTreeCascVarNegTransMomMC(0),
   fTreeCascVarNegNSigmaPion(0),
   fTreeCascVarNegNSigmaProton(0),
   fTreeCascVarPosNSigmaPion(0),
   fTreeCascVarPosNSigmaProton(0),
   fTreeCascVarBachNSigmaPion(0),
   fTreeCascVarBachNSigmaKaon(0),

   fTreeCascVarkITSRefitBachelor(0),
   fTreeCascVarkITSRefitNegative(0),
   fTreeCascVarkITSRefitPositive(0),
   
    fTreeCascVarNegClusters(0),
    fTreeCascVarPosClusters(0),
    fTreeCascVarBachClusters(0),
    fTreeCascVarNegSharedClusters(0),
    fTreeCascVarPosSharedClusters(0),
    fTreeCascVarBachSharedClusters(0),

   fTreeCascVarEvHasXiMinus(0),
   fTreeCascVarEvHasXiPlus(0),
   fTreeCascVarEvHasOmegaMinus(0),
   fTreeCascVarEvHasOmegaPlus(0),
   fTreeCascVarEvHasLambda(0),
   fTreeCascVarEvHasAntiLambda(0),

   fTreeCascVarEvHasLowPtXiMinus(0),
   fTreeCascVarEvHasLowPtXiPlus(0),
   fTreeCascVarEvHasLowPtOmegaMinus(0),
   fTreeCascVarEvHasLowPtOmegaPlus(0),
   fTreeCascVarEvHasLowPtLambda(0),
   fTreeCascVarEvHasLowPtAntiLambda(0),

   fTreeCascVarEvHasVeryLowPtXiMinus(0),
   fTreeCascVarEvHasVeryLowPtXiPlus(0),
   fTreeCascVarEvHasVeryLowPtOmegaMinus(0),
   fTreeCascVarEvHasVeryLowPtOmegaPlus(0),
   fTreeCascVarEvHasVeryLowPtLambda(0),
   fTreeCascVarEvHasVeryLowPtAntiLambda(0),
   
    fTreeCascVarBadCascadeJai(0),
    fTreeCascVarDeltaDCA(0),

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

//------------------------------------------------
// PARTICLE HISTOGRAMS
// --- Filled on a Particle-by-Particle basis
//------------------------------------------------
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

//MC Mult

f3dHistGenPtVsYVsMultMCXiMinus(0),
f3dHistGenPtVsYVsMultMCXiPlus(0),
f3dHistGenPtVsYVsMultMCOmegaMinus(0),
f3dHistGenPtVsYVsMultMCOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultMCXiMinus(0),
f3dHistGenSelectedPtVsYVsMultMCXiPlus(0),
f3dHistGenSelectedPtVsYVsMultMCOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultMCOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultMCXiMinus(0),
f3dHistGenPtVsYCMSVsMultMCXiPlus(0),
f3dHistGenPtVsYCMSVsMultMCOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultMCOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultMCXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultMCXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultMCOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultMCOmegaPlus(0),

//V0A

f3dHistGenPtVsYVsMultV0AXiMinus(0),
f3dHistGenPtVsYVsMultV0AXiPlus(0),
f3dHistGenPtVsYVsMultV0AOmegaMinus(0),
f3dHistGenPtVsYVsMultV0AOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultV0AXiMinus(0),
f3dHistGenSelectedPtVsYVsMultV0AXiPlus(0),
f3dHistGenSelectedPtVsYVsMultV0AOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultV0AOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultV0AXiMinus(0),
f3dHistGenPtVsYCMSVsMultV0AXiPlus(0),
f3dHistGenPtVsYCMSVsMultV0AOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultV0AOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultV0AXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultV0AXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaPlus(0),
f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiMinus(0),
f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiPlus(0),
f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaMinus(0),
f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaPlus(0),

//ZNA

f3dHistGenPtVsYVsMultZNAXiMinus(0),
f3dHistGenPtVsYVsMultZNAXiPlus(0),
f3dHistGenPtVsYVsMultZNAOmegaMinus(0),
f3dHistGenPtVsYVsMultZNAOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultZNAXiMinus(0),
f3dHistGenSelectedPtVsYVsMultZNAXiPlus(0),
f3dHistGenSelectedPtVsYVsMultZNAOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultZNAOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultZNAXiMinus(0),
f3dHistGenPtVsYCMSVsMultZNAXiPlus(0),
f3dHistGenPtVsYCMSVsMultZNAOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultZNAOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultZNAXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultZNAXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaPlus(0),

//TRK

f3dHistGenPtVsYVsMultTRKXiMinus(0),
f3dHistGenPtVsYVsMultTRKXiPlus(0),
f3dHistGenPtVsYVsMultTRKOmegaMinus(0),
f3dHistGenPtVsYVsMultTRKOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultTRKXiMinus(0),
f3dHistGenSelectedPtVsYVsMultTRKXiPlus(0),
f3dHistGenSelectedPtVsYVsMultTRKOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultTRKOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultTRKXiMinus(0),
f3dHistGenPtVsYCMSVsMultTRKXiPlus(0),
f3dHistGenPtVsYCMSVsMultTRKOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultTRKOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultTRKXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultTRKXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaPlus(0),

//SPD

f3dHistGenPtVsYVsMultSPDXiMinus(0),
f3dHistGenPtVsYVsMultSPDXiPlus(0),
f3dHistGenPtVsYVsMultSPDOmegaMinus(0),
f3dHistGenPtVsYVsMultSPDOmegaPlus(0),
f3dHistGenSelectedPtVsYVsMultSPDXiMinus(0),
f3dHistGenSelectedPtVsYVsMultSPDXiPlus(0),
f3dHistGenSelectedPtVsYVsMultSPDOmegaMinus(0),
f3dHistGenSelectedPtVsYVsMultSPDOmegaPlus(0),
f3dHistGenPtVsYCMSVsMultSPDXiMinus(0),
f3dHistGenPtVsYCMSVsMultSPDXiPlus(0),
f3dHistGenPtVsYCMSVsMultSPDOmegaMinus(0),
f3dHistGenPtVsYCMSVsMultSPDOmegaPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultSPDXiMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultSPDXiPlus(0),
f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaMinus(0),
f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaPlus(0),


   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0)
{
   // Constructor

   //Set Variables for re-running the cascade vertexers (as done for MS paper)
        
        // New Loose : 1st step for the 7 TeV pp analysis
        
        fV0VertexerSels[0] =  33.  ;  // max allowed chi2
        fV0VertexerSels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
        fV0VertexerSels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
        fV0VertexerSels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
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
        
   // Output slot #0 writes into a TList container (Cascade)
   DefineOutput(1, TList::Class());
   DefineOutput(2, TTree::Class());
}


AliAnalysisTaskExtractPerformanceCascade::~AliAnalysisTaskExtractPerformanceCascade()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------

   if (fListHist){
      delete fListHist;
      fListHist = 0x0;
   }
   if (fTreeCascade){
      delete fTreeCascade;
      fTreeCascade = 0x0;
   }
    //cleanup esd track cuts object too...
   if (fESDtrackCuts){
    delete fESDtrackCuts;
    fESDtrackCuts = 0x0; 
  }
  if (fUtils){
    delete fUtils;
    fUtils = 0x0;
  }

}

//________________________________________________________________________
void AliAnalysisTaskExtractPerformanceCascade::UserCreateOutputObjects()
{
   OpenFile(2);	
   // Called once

//------------------------------------------------

   fTreeCascade = new TTree("fTreeCascade","CascadeCandidates");

//------------------------------------------------
// fTreeCascade Branch definitions - Cascade Tree
//------------------------------------------------

//------------------------------------------------
// fTreeCascade Branch definitions
//------------------------------------------------

//-----------BASIC-INFO---------------------------
/* 1*/		fTreeCascade->Branch("fTreeCascVarCharge",&fTreeCascVarCharge,"fTreeCascVarCharge/I");	
/* 2*/		fTreeCascade->Branch("fTreeCascVarMassAsXi",&fTreeCascVarMassAsXi,"fTreeCascVarMassAsXi/F");
/* 3*/		fTreeCascade->Branch("fTreeCascVarMassAsOmega",&fTreeCascVarMassAsOmega,"fTreeCascVarMassAsOmega/F");
/* 4*/		fTreeCascade->Branch("fTreeCascVarPt",&fTreeCascVarPt,"fTreeCascVarPt/F");
/* 5*/		fTreeCascade->Branch("fTreeCascVarPtMC",&fTreeCascVarPtMC,"fTreeCascVarPtMC/F");
/* 6*/		fTreeCascade->Branch("fTreeCascVarRapXi",&fTreeCascVarRapXi,"fTreeCascVarRapXi/F");
/* 7*/		fTreeCascade->Branch("fTreeCascVarRapMC",&fTreeCascVarRapMC,"fTreeCascVarRapMC/F");
/* 8*/		fTreeCascade->Branch("fTreeCascVarRapOmega",&fTreeCascVarRapOmega,"fTreeCascVarRapOmega/F");
/* 9*/		fTreeCascade->Branch("fTreeCascVarNegEta",&fTreeCascVarNegEta,"fTreeCascVarNegEta/F");
/*10*/		fTreeCascade->Branch("fTreeCascVarPosEta",&fTreeCascVarPosEta,"fTreeCascVarPosEta/F");
/*11*/		fTreeCascade->Branch("fTreeCascVarBachEta",&fTreeCascVarBachEta,"fTreeCascVarBachEta/F");
//-----------INFO-FOR-CUTS------------------------
/*12*/		fTreeCascade->Branch("fTreeCascVarDCACascDaughters",&fTreeCascVarDCACascDaughters,"fTreeCascVarDCACascDaughters/F");
/*13*/		fTreeCascade->Branch("fTreeCascVarDCABachToPrimVtx",&fTreeCascVarDCABachToPrimVtx,"fTreeCascVarDCABachToPrimVtx/F");
/*14*/		fTreeCascade->Branch("fTreeCascVarDCAV0Daughters",&fTreeCascVarDCAV0Daughters,"fTreeCascVarDCAV0Daughters/F");
/*15*/		fTreeCascade->Branch("fTreeCascVarDCAV0ToPrimVtx",&fTreeCascVarDCAV0ToPrimVtx,"fTreeCascVarDCAV0ToPrimVtx/F");
/*16*/		fTreeCascade->Branch("fTreeCascVarDCAPosToPrimVtx",&fTreeCascVarDCAPosToPrimVtx,"fTreeCascVarDCAPosToPrimVtx/F");
/*17*/		fTreeCascade->Branch("fTreeCascVarDCANegToPrimVtx",&fTreeCascVarDCANegToPrimVtx,"fTreeCascVarDCANegToPrimVtx/F");
/*18*/		fTreeCascade->Branch("fTreeCascVarCascCosPointingAngle",&fTreeCascVarCascCosPointingAngle,"fTreeCascVarCascCosPointingAngle/F");
/*19*/		fTreeCascade->Branch("fTreeCascVarCascRadius",&fTreeCascVarCascRadius,"fTreeCascVarCascRadius/F");
/*20*/		fTreeCascade->Branch("fTreeCascVarV0Mass",&fTreeCascVarV0Mass,"fTreeCascVarV0Mass/F");
/*21*/		fTreeCascade->Branch("fTreeCascVarV0CosPointingAngle",&fTreeCascVarV0CosPointingAngle,"fTreeCascVarV0CosPointingAngle/F");
/*21*/		fTreeCascade->Branch("fTreeCascVarV0CosPointingAngleSpecial",&fTreeCascVarV0CosPointingAngleSpecial,"fTreeCascVarV0CosPointingAngleSpecial/F");
/*22*/		fTreeCascade->Branch("fTreeCascVarV0Radius",&fTreeCascVarV0Radius,"fTreeCascVarV0Radius/F");
/*23*/		fTreeCascade->Branch("fTreeCascVarLeastNbrClusters",&fTreeCascVarLeastNbrClusters,"fTreeCascVarLeastNbrClusters/I");

/*23*/		fTreeCascade->Branch("fTreeCascVarNegClusters",&fTreeCascVarNegClusters,"fTreeCascVarNegClusters/I");
/*23*/		fTreeCascade->Branch("fTreeCascVarPosClusters",&fTreeCascVarPosClusters,"fTreeCascVarPosClusters/I");
/*23*/		fTreeCascade->Branch("fTreeCascVarBachClusters",&fTreeCascVarBachClusters,"fTreeCascVarBachClusters/I");
/*23*/		fTreeCascade->Branch("fTreeCascVarNegSharedClusters",&fTreeCascVarNegSharedClusters,"fTreeCascVarNegSharedClusters/I");
/*23*/		fTreeCascade->Branch("fTreeCascVarPosSharedClusters",&fTreeCascVarPosSharedClusters,"fTreeCascVarPosSharedClusters/I");
/*23*/		fTreeCascade->Branch("fTreeCascVarBachSharedClusters",&fTreeCascVarBachSharedClusters,"fTreeCascVarBachSharedClusters/I");
//-----------MULTIPLICITY-INFO--------------------
/*24*/		fTreeCascade->Branch("fTreeCascVarMultiplicity",&fTreeCascVarMultiplicity,"fTreeCascVarMultiplicity/I");
/*24*/		fTreeCascade->Branch("fTreeCascVarMultiplicityV0A",&fTreeCascVarMultiplicityV0A,"fTreeCascVarMultiplicityV0A/I");
/*24*/		fTreeCascade->Branch("fTreeCascVarMultiplicityZNA",&fTreeCascVarMultiplicityZNA,"fTreeCascVarMultiplicityZNA/I");
/*24*/		fTreeCascade->Branch("fTreeCascVarMultiplicityTRK",&fTreeCascVarMultiplicityTRK,"fTreeCascVarMultiplicityTRK/I");
/*24*/		fTreeCascade->Branch("fTreeCascVarMultiplicitySPD",&fTreeCascVarMultiplicitySPD,"fTreeCascVarMultiplicitySPD/I");
/*24*/		fTreeCascade->Branch("fTreeCascVarMultiplicityMC",&fTreeCascVarMultiplicityMC,"fTreeCascVarMultiplicityMC/I");
//-----------DECAY-LENGTH-INFO--------------------
/*25*/		fTreeCascade->Branch("fTreeCascVarDistOverTotMom",&fTreeCascVarDistOverTotMom,"fTreeCascVarDistOverTotMom/F");
//-----------MC-PID-------------------------------
/*25bis*/ fTreeCascade->Branch("fTreeCascVarIsPhysicalPrimary",&fTreeCascVarIsPhysicalPrimary,"fTreeCascVarIsPhysicalPrimary/I");
/*26*/		fTreeCascade->Branch("fTreeCascVarPID",&fTreeCascVarPID,"fTreeCascVarPID/I");
/*26*/		fTreeCascade->Branch("fTreeCascVarPIDSwapped",&fTreeCascVarPIDSwapped,"fTreeCascVarPIDSwapped/I");
/*27*/		fTreeCascade->Branch("fTreeCascVarPIDBachelor",&fTreeCascVarPIDBachelor,"fTreeCascVarPIDBachelor/I");
/*28*/    fTreeCascade->Branch("fTreeCascVarPIDNegative",&fTreeCascVarPIDNegative,"fTreeCascVarPIDNegative/I");
/*29*/    fTreeCascade->Branch("fTreeCascVarPIDPositive",&fTreeCascVarPIDPositive,"fTreeCascVarPIDPositive/I");
/*30*/		fTreeCascade->Branch("fTreeCascVarBachTransMom",&fTreeCascVarBachTransMom,"fTreeCascVarBachTransMom/F");
/*30*/		fTreeCascade->Branch("fTreeCascVarPosTransMom",&fTreeCascVarPosTransMom,"fTreeCascVarPosTransMom/F");
/*31*/		fTreeCascade->Branch("fTreeCascVarNegTransMom",&fTreeCascVarNegTransMom,"fTreeCascVarNegTransMom/F");
/*32*/		fTreeCascade->Branch("fTreeCascVarPosTransMomMC",&fTreeCascVarPosTransMomMC,"fTreeCascVarPosTransMomMC/F");
/*33*/		fTreeCascade->Branch("fTreeCascVarNegTransMomMC",&fTreeCascVarNegTransMomMC,"fTreeCascVarNegTransMomMC/F");
//------------------------------------------------
/*34*/		fTreeCascade->Branch("fTreeCascVarNegNSigmaPion",&fTreeCascVarNegNSigmaPion,"fTreeCascVarNegNSigmaPion/F");
/*35*/		fTreeCascade->Branch("fTreeCascVarNegNSigmaProton",&fTreeCascVarNegNSigmaProton,"fTreeCascVarNegNSigmaProton/F");
/*36*/		fTreeCascade->Branch("fTreeCascVarPosNSigmaPion",&fTreeCascVarPosNSigmaPion,"fTreeCascVarPosNSigmaPion/F");
/*37*/		fTreeCascade->Branch("fTreeCascVarPosNSigmaProton",&fTreeCascVarPosNSigmaProton,"fTreeCascVarPosNSigmaProton/F");
/*38*/		fTreeCascade->Branch("fTreeCascVarBachNSigmaPion",&fTreeCascVarBachNSigmaPion,"fTreeCascVarBachNSigmaPion/F");
/*39*/		fTreeCascade->Branch("fTreeCascVarBachNSigmaKaon",&fTreeCascVarBachNSigmaKaon,"fTreeCascVarBachNSigmaKaon/F");

/*29*/		fTreeCascade->Branch("fTreeCascVarkITSRefitBachelor",&fTreeCascVarkITSRefitBachelor,"fTreeCascVarkITSRefitBachelor/O");
/*29*/		fTreeCascade->Branch("fTreeCascVarkITSRefitNegative",&fTreeCascVarkITSRefitNegative,"fTreeCascVarkITSRefitNegative/O");
/*29*/		fTreeCascade->Branch("fTreeCascVarkITSRefitPositive",&fTreeCascVarkITSRefitPositive,"fTreeCascVarkITSRefitPositive/O");

      fTreeCascade->Branch("fTreeCascVarNegClusters",&fTreeCascVarNegClusters,"fTreeCascVarNegClusters/I");
      fTreeCascade->Branch("fTreeCascVarPosClusters",&fTreeCascVarPosClusters,"fTreeCascVarPosClusters/I");
      fTreeCascade->Branch("fTreeCascVarBachClusters",&fTreeCascVarBachClusters,"fTreeCascVarBachClusters/I");
      fTreeCascade->Branch("fTreeCascVarNegSharedClusters",&fTreeCascVarNegSharedClusters,"fTreeCascVarNegSharedClusters/I");
      fTreeCascade->Branch("fTreeCascVarPosSharedClusters",&fTreeCascVarPosSharedClusters,"fTreeCascVarPosSharedClusters/I");
      fTreeCascade->Branch("fTreeCascVarBachSharedClusters",&fTreeCascVarBachSharedClusters,"fTreeCascVarBachSharedClusters/I");

/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasXiMinus",   &fTreeCascVarEvHasXiMinus,   "fTreeCascVarEvHasXiMinus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasXiPlus",    &fTreeCascVarEvHasXiPlus,    "fTreeCascVarEvHasXiPlus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasOmegaMinus",&fTreeCascVarEvHasOmegaMinus,"fTreeCascVarEvHasOmegaMinus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasOmegaPlus", &fTreeCascVarEvHasOmegaPlus, "fTreeCascVarEvHasOmegaPlus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasLambda",    &fTreeCascVarEvHasLambda,    "fTreeCascVarEvHasLambda/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasAntiLambda",&fTreeCascVarEvHasAntiLambda,"fTreeCascVarEvHasAntiLambda/O");

/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasLowPtXiMinus",   &fTreeCascVarEvHasLowPtXiMinus,   "fTreeCascVarEvHasLowPtXiMinus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasLowPtXiPlus",    &fTreeCascVarEvHasLowPtXiPlus,    "fTreeCascVarEvHasLowPtXiPlus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasLowPtOmegaMinus",&fTreeCascVarEvHasLowPtOmegaMinus,"fTreeCascVarEvHasLowPtOmegaMinus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasLowPtOmegaPlus", &fTreeCascVarEvHasLowPtOmegaPlus, "fTreeCascVarEvHasLowPtOmegaPlus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasLowPtLambda",    &fTreeCascVarEvHasLowPtLambda,    "fTreeCascVarEvHasLowPtLambda/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasLowPtAntiLambda",&fTreeCascVarEvHasLowPtAntiLambda,"fTreeCascVarEvHasLowPtAntiLambda/O");

/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasVeryLowPtXiMinus",   &fTreeCascVarEvHasVeryLowPtXiMinus,   "fTreeCascVarEvHasVeryLowPtXiMinus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasVeryLowPtXiPlus",    &fTreeCascVarEvHasVeryLowPtXiPlus,    "fTreeCascVarEvHasVeryLowPtXiPlus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasVeryLowPtOmegaMinus",&fTreeCascVarEvHasVeryLowPtOmegaMinus,"fTreeCascVarEvHasVeryLowPtOmegaMinus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasVeryLowPtOmegaPlus", &fTreeCascVarEvHasVeryLowPtOmegaPlus, "fTreeCascVarEvHasVeryLowPtOmegaPlus/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasVeryLowPtLambda",    &fTreeCascVarEvHasVeryLowPtLambda,    "fTreeCascVarEvHasVeryLowPtLambda/O");
/*39*/		fTreeCascade->Branch("fTreeCascVarEvHasVeryLowPtAntiLambda",&fTreeCascVarEvHasVeryLowPtAntiLambda,"fTreeCascVarEvHasVeryLowPtAntiLambda/O");

      fTreeCascade->Branch("fTreeCascVarBadCascadeJai",&fTreeCascVarBadCascadeJai,"fTreeCascVarBadCascadeJai/O");
      fTreeCascade->Branch("fTreeCascVarDeltaDCA",&fTreeCascVarDeltaDCA,"fTreeCascVarDeltaDCA/F");

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
  if(! fUtils ){
    fUtils = new AliAnalysisUtils();
  }

//------------------------------------------------
// V0 Multiplicity Histograms
//------------------------------------------------

   // Create histograms
   OpenFile(1);
   fListHist = new TList();
   fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner


   if(! fHistV0MultiplicityBeforeTrigSel) {
      fHistV0MultiplicityBeforeTrigSel = new TH1F("fHistV0MultiplicityBeforeTrigSel", 
         "V0s per event (before Trig. Sel.);Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHist->Add(fHistV0MultiplicityBeforeTrigSel);
   }
           
   if(! fHistV0MultiplicityForTrigEvt) {
      fHistV0MultiplicityForTrigEvt = new TH1F("fHistV0MultiplicityForTrigEvt", 
         "V0s per event (for triggered evt);Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHist->Add(fHistV0MultiplicityForTrigEvt);
   }

   if(! fHistV0MultiplicityForSelEvt) {
      fHistV0MultiplicityForSelEvt = new TH1F("fHistV0MultiplicityForSelEvt", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHist->Add(fHistV0MultiplicityForSelEvt);
   }

   if(! fHistV0MultiplicityForSelEvtNoTPCOnly) {
      fHistV0MultiplicityForSelEvtNoTPCOnly = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnly", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHist->Add(fHistV0MultiplicityForSelEvtNoTPCOnly);
   }
   if(! fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup) {
      fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHist->Add(fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup);
   }

//------------------------------------------------
// Track Multiplicity Histograms
//------------------------------------------------

   if(! fHistMultiplicityBeforeTrigSel) {
      fHistMultiplicityBeforeTrigSel = new TH1F("fHistMultiplicityBeforeTrigSel", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHist->Add(fHistMultiplicityBeforeTrigSel);
   }
   if(! fHistMultiplicityForTrigEvt) {
      fHistMultiplicityForTrigEvt = new TH1F("fHistMultiplicityForTrigEvt", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHist->Add(fHistMultiplicityForTrigEvt);
   }
   if(! fHistMultiplicity) {
      fHistMultiplicity = new TH1F("fHistMultiplicity", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHist->Add(fHistMultiplicity);
   }
   if(! fHistMultiplicityNoTPCOnly) {
      fHistMultiplicityNoTPCOnly = new TH1F("fHistMultiplicityNoTPCOnly", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHist->Add(fHistMultiplicityNoTPCOnly);
   }
   if(! fHistMultiplicityNoTPCOnlyNoPileup) {
      fHistMultiplicityNoTPCOnlyNoPileup = new TH1F("fHistMultiplicityNoTPCOnlyNoPileup", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHist->Add(fHistMultiplicityNoTPCOnlyNoPileup);
   }
  
  
  //V0A Centrality (if PbPb / pPb)
  if(! fHistMultiplicityV0ABeforeTrigSel) {
    fHistMultiplicityV0ABeforeTrigSel = new TH1F("fHistMultiplicityV0ABeforeTrigSel",
                                                 "Centrality Distribution: V0A;V0A Centrality;Events",
                                                 200, 0, 200);
    fListHist->Add(fHistMultiplicityV0ABeforeTrigSel);
  }
  if(! fHistMultiplicityV0AForTrigEvt) {
    fHistMultiplicityV0AForTrigEvt = new TH1F("fHistMultiplicityV0AForTrigEvt",
                                              "Centrality Distribution: V0A;V0A Centrality;Events",
                                              200, 0, 200);
    fListHist->Add(fHistMultiplicityV0AForTrigEvt);
  }
  if(! fHistMultiplicityV0A) {
    fHistMultiplicityV0A = new TH1F("fHistMultiplicityV0A",
                                    "Centrality Distribution: V0A;V0A Centrality;Events",
                                    200, 0, 200);
    fListHist->Add(fHistMultiplicityV0A);
  }
  if(! fHistMultiplicityV0ANoTPCOnly) {
    fHistMultiplicityV0ANoTPCOnly = new TH1F("fHistMultiplicityV0ANoTPCOnly",
                                             "Centrality Distribution: V0A;V0A Centrality;Events",
                                             200, 0, 200);
    fListHist->Add(fHistMultiplicityV0ANoTPCOnly);
  }
  if(! fHistMultiplicityV0ANoTPCOnlyNoPileup) {
    fHistMultiplicityV0ANoTPCOnlyNoPileup = new TH1F("fHistMultiplicityV0ANoTPCOnlyNoPileup",
                                                     "Centrality Distribution: V0A;V0A Centrality;Events",
                                                     200, 0, 200);
    fListHist->Add(fHistMultiplicityV0ANoTPCOnlyNoPileup);
  }
  
  //ZNA Centrality (if PbPb / pPb)
  if(! fHistMultiplicityZNABeforeTrigSel) {
    fHistMultiplicityZNABeforeTrigSel = new TH1F("fHistMultiplicityZNABeforeTrigSel",
                                                 "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                                 200, 0, 200);
    fListHist->Add(fHistMultiplicityZNABeforeTrigSel);
  }
  if(! fHistMultiplicityZNAForTrigEvt) {
    fHistMultiplicityZNAForTrigEvt = new TH1F("fHistMultiplicityZNAForTrigEvt",
                                              "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                              200, 0, 200);
    fListHist->Add(fHistMultiplicityZNAForTrigEvt);
  }
  if(! fHistMultiplicityZNA) {
    fHistMultiplicityZNA = new TH1F("fHistMultiplicityZNA",
                                    "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                    200, 0, 200);
    fListHist->Add(fHistMultiplicityZNA);
  }
  if(! fHistMultiplicityZNANoTPCOnly) {
    fHistMultiplicityZNANoTPCOnly = new TH1F("fHistMultiplicityZNANoTPCOnly",
                                             "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                             200, 0, 200);
    fListHist->Add(fHistMultiplicityZNANoTPCOnly);
  }
  if(! fHistMultiplicityZNANoTPCOnlyNoPileup) {
    fHistMultiplicityZNANoTPCOnlyNoPileup = new TH1F("fHistMultiplicityZNANoTPCOnlyNoPileup",
                                                     "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                                     200, 0, 200);
    fListHist->Add(fHistMultiplicityZNANoTPCOnlyNoPileup);
  }
  
  //TRK Centrality (if PbPb / pPb)
  if(! fHistMultiplicityTRKBeforeTrigSel) {
    fHistMultiplicityTRKBeforeTrigSel = new TH1F("fHistMultiplicityTRKBeforeTrigSel",
                                                 "Centrality Distribution: TRK;TRK Centrality;Events",
                                                 200, 0, 200);
    fListHist->Add(fHistMultiplicityTRKBeforeTrigSel);
  }
  if(! fHistMultiplicityTRKForTrigEvt) {
    fHistMultiplicityTRKForTrigEvt = new TH1F("fHistMultiplicityTRKForTrigEvt",
                                              "Centrality Distribution: TRK;TRK Centrality;Events",
                                              200, 0, 200);
    fListHist->Add(fHistMultiplicityTRKForTrigEvt);
  }
  if(! fHistMultiplicityTRK) {
    fHistMultiplicityTRK = new TH1F("fHistMultiplicityTRK",
                                    "Centrality Distribution: TRK;TRK Centrality;Events",
                                    200, 0, 200);
    fListHist->Add(fHistMultiplicityTRK);
  }
  if(! fHistMultiplicityTRKNoTPCOnly) {
    fHistMultiplicityTRKNoTPCOnly = new TH1F("fHistMultiplicityTRKNoTPCOnly",
                                             "Centrality Distribution: TRK;TRK Centrality;Events",
                                             200, 0, 200);
    fListHist->Add(fHistMultiplicityTRKNoTPCOnly);
  }
  if(! fHistMultiplicityTRKNoTPCOnlyNoPileup) {
    fHistMultiplicityTRKNoTPCOnlyNoPileup = new TH1F("fHistMultiplicityTRKNoTPCOnlyNoPileup",
                                                     "Centrality Distribution: TRK;TRK Centrality;Events",
                                                     200, 0, 200);
    fListHist->Add(fHistMultiplicityTRKNoTPCOnlyNoPileup);
  }
  
  //SPD Centrality (if PbPb / pPb)
  if(! fHistMultiplicitySPDBeforeTrigSel) {
    fHistMultiplicitySPDBeforeTrigSel = new TH1F("fHistMultiplicitySPDBeforeTrigSel",
                                                 "Centrality Distribution: SPD;SPD Centrality;Events",
                                                 200, 0, 200);
    fListHist->Add(fHistMultiplicitySPDBeforeTrigSel);
  }
  if(! fHistMultiplicitySPDForTrigEvt) {
    fHistMultiplicitySPDForTrigEvt = new TH1F("fHistMultiplicitySPDForTrigEvt",
                                              "Centrality Distribution: SPD;SPD Centrality;Events",
                                              200, 0, 200);
    fListHist->Add(fHistMultiplicitySPDForTrigEvt);
  }
  if(! fHistMultiplicitySPD) {
    fHistMultiplicitySPD = new TH1F("fHistMultiplicitySPD",
                                    "Centrality Distribution: SPD;SPD Centrality;Events",
                                    200, 0, 200);
    fListHist->Add(fHistMultiplicitySPD);
  }
  if(! fHistMultiplicitySPDNoTPCOnly) {
    fHistMultiplicitySPDNoTPCOnly = new TH1F("fHistMultiplicitySPDNoTPCOnly",
                                             "Centrality Distribution: SPD;SPD Centrality;Events",
                                             200, 0, 200);
    fListHist->Add(fHistMultiplicitySPDNoTPCOnly);
  }
  if(! fHistMultiplicitySPDNoTPCOnlyNoPileup) {
    fHistMultiplicitySPDNoTPCOnlyNoPileup = new TH1F("fHistMultiplicitySPDNoTPCOnlyNoPileup",
                                                     "Centrality Distribution: SPD;SPD Centrality;Events",
                                                     200, 0, 200);
    fListHist->Add(fHistMultiplicitySPDNoTPCOnlyNoPileup);
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

//--------------------------------------------------------------------------------------
//--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, all generated

   if(! f3dHistGenPtVsYVsMultXiMinus) {
      f3dHistGenPtVsYVsMultXiMinus = new TH3F( "f3dHistGenPtVsYVsMultXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenPtVsYVsMultXiMinus);
   }
   if(! f3dHistGenPtVsYVsMultXiPlus) {
      f3dHistGenPtVsYVsMultXiPlus = new TH3F( "f3dHistGenPtVsYVsMultXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenPtVsYVsMultXiPlus);
   }
//--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus

   if(! f3dHistGenPtVsYVsMultOmegaMinus) {
      f3dHistGenPtVsYVsMultOmegaMinus = new TH3F( "f3dHistGenPtVsYVsMultOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenPtVsYVsMultOmegaMinus);
   }
   if(! f3dHistGenPtVsYVsMultOmegaPlus) {
      f3dHistGenPtVsYVsMultOmegaPlus = new TH3F( "f3dHistGenPtVsYVsMultOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenPtVsYVsMultOmegaPlus);
   }

//All generated cascades, YCMS

   if(! f3dHistGenPtVsYCMSVsMultXiMinus) {
      f3dHistGenPtVsYCMSVsMultXiMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenPtVsYCMSVsMultXiMinus);
   }
   if(! f3dHistGenPtVsYCMSVsMultXiPlus) {
      f3dHistGenPtVsYCMSVsMultXiPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenPtVsYCMSVsMultXiPlus);
   }
//--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus

   if(! f3dHistGenPtVsYCMSVsMultOmegaMinus) {
      f3dHistGenPtVsYCMSVsMultOmegaMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenPtVsYCMSVsMultOmegaMinus);
   }
   if(! f3dHistGenPtVsYCMSVsMultOmegaPlus) {
      f3dHistGenPtVsYCMSVsMultOmegaPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenPtVsYCMSVsMultOmegaPlus);
   }


//--------------------------------------------------------------------------------------
//--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, at selected analysis evts

   if(! f3dHistGenSelectedPtVsYVsMultXiMinus) {
      f3dHistGenSelectedPtVsYVsMultXiMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenSelectedPtVsYVsMultXiMinus);
   }
   if(! f3dHistGenSelectedPtVsYVsMultXiPlus) {
      f3dHistGenSelectedPtVsYVsMultXiPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenSelectedPtVsYVsMultXiPlus);
   }
//--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus

   if(! f3dHistGenSelectedPtVsYVsMultOmegaMinus) {
      f3dHistGenSelectedPtVsYVsMultOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenSelectedPtVsYVsMultOmegaMinus);
   }
   if(! f3dHistGenSelectedPtVsYVsMultOmegaPlus) {
      f3dHistGenSelectedPtVsYVsMultOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenSelectedPtVsYVsMultOmegaPlus);
   }

//ANALYSIS level Cascades, YCMS


   if(! f3dHistGenSelectedPtVsYCMSVsMultXiMinus) {
      f3dHistGenSelectedPtVsYCMSVsMultXiMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultXiMinus);
   }
   if(! f3dHistGenSelectedPtVsYCMSVsMultXiPlus) {
      f3dHistGenSelectedPtVsYCMSVsMultXiPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultXiPlus);
   }
//--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus

   if(! f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus) {
      f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus);
   }
   if(! f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus) {
      f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus);
   }
  
  //======================================================================================
  //--------------------------------------------------------------------------------------
  // True Generated (For Multiplicity Unfolding)
  //--------------------------------------------------------------------------------------
  //======================================================================================
  
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, all generated
  
  if(! f3dHistGenPtVsYVsMultMCXiMinus) {
    f3dHistGenPtVsYVsMultMCXiMinus = new TH3F( "f3dHistGenPtVsYVsMultMCXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenPtVsYVsMultMCXiMinus);
  }
  if(! f3dHistGenPtVsYVsMultMCXiPlus) {
    f3dHistGenPtVsYVsMultMCXiPlus = new TH3F( "f3dHistGenPtVsYVsMultMCXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenPtVsYVsMultMCXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYVsMultMCOmegaMinus) {
    f3dHistGenPtVsYVsMultMCOmegaMinus = new TH3F( "f3dHistGenPtVsYVsMultMCOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenPtVsYVsMultMCOmegaMinus);
  }
  if(! f3dHistGenPtVsYVsMultMCOmegaPlus) {
    f3dHistGenPtVsYVsMultMCOmegaPlus = new TH3F( "f3dHistGenPtVsYVsMultMCOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenPtVsYVsMultMCOmegaPlus);
  }
  
  //All generated cascades, YCMS
  
  if(! f3dHistGenPtVsYCMSVsMultMCXiMinus) {
    f3dHistGenPtVsYCMSVsMultMCXiMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultMCXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenPtVsYCMSVsMultMCXiMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultMCXiPlus) {
    f3dHistGenPtVsYCMSVsMultMCXiPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultMCXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenPtVsYCMSVsMultMCXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYCMSVsMultMCOmegaMinus) {
    f3dHistGenPtVsYCMSVsMultMCOmegaMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultMCOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenPtVsYCMSVsMultMCOmegaMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultMCOmegaPlus) {
    f3dHistGenPtVsYCMSVsMultMCOmegaPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultMCOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenPtVsYCMSVsMultMCOmegaPlus);
  }
  

  //--------------------------------------------------------------------------------------
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, at selected analysis evts
  
  if(! f3dHistGenSelectedPtVsYVsMultMCXiMinus) {
    f3dHistGenSelectedPtVsYVsMultMCXiMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultMCXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenSelectedPtVsYVsMultMCXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultMCXiPlus) {
    f3dHistGenSelectedPtVsYVsMultMCXiPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultMCXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenSelectedPtVsYVsMultMCXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYVsMultMCOmegaMinus) {
    f3dHistGenSelectedPtVsYVsMultMCOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultMCOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenSelectedPtVsYVsMultMCOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultMCOmegaPlus) {
    f3dHistGenSelectedPtVsYVsMultMCOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultMCOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenSelectedPtVsYVsMultMCOmegaPlus);
  }
  
  //ANALYSIS level Cascades, YCMS
  
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultMCXiMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultMCXiMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultMCXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultMCXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultMCXiPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultMCXiPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultMCXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultMCXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultMCOmegaMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultMCOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultMCOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultMCOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultMCOmegaPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultMCOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultMCOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultMCOmegaPlus);
  }

  //======================================================================================
  //--------------------------------------------------------------------------------------
  // V0A
  //--------------------------------------------------------------------------------------
  //======================================================================================
  
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, all generated
  
  if(! f3dHistGenPtVsYVsMultV0AXiMinus) {
    f3dHistGenPtVsYVsMultV0AXiMinus = new TH3F( "f3dHistGenPtVsYVsMultV0AXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultV0AXiMinus);
  }
  if(! f3dHistGenPtVsYVsMultV0AXiPlus) {
    f3dHistGenPtVsYVsMultV0AXiPlus = new TH3F( "f3dHistGenPtVsYVsMultV0AXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultV0AXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYVsMultV0AOmegaMinus) {
    f3dHistGenPtVsYVsMultV0AOmegaMinus = new TH3F( "f3dHistGenPtVsYVsMultV0AOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultV0AOmegaMinus);
  }
  if(! f3dHistGenPtVsYVsMultV0AOmegaPlus) {
    f3dHistGenPtVsYVsMultV0AOmegaPlus = new TH3F( "f3dHistGenPtVsYVsMultV0AOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultV0AOmegaPlus);
  }
  
  //All generated cascades, YCMS
  
  if(! f3dHistGenPtVsYCMSVsMultV0AXiMinus) {
    f3dHistGenPtVsYCMSVsMultV0AXiMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultV0AXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultV0AXiMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultV0AXiPlus) {
    f3dHistGenPtVsYCMSVsMultV0AXiPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultV0AXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultV0AXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYCMSVsMultV0AOmegaMinus) {
    f3dHistGenPtVsYCMSVsMultV0AOmegaMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultV0AOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultV0AOmegaMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultV0AOmegaPlus) {
    f3dHistGenPtVsYCMSVsMultV0AOmegaPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultV0AOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultV0AOmegaPlus);
  }
  
  
  //--------------------------------------------------------------------------------------
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, at selected analysis evts
  
  if(! f3dHistGenSelectedPtVsYVsMultV0AXiMinus) {
    f3dHistGenSelectedPtVsYVsMultV0AXiMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultV0AXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultV0AXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultV0AXiPlus) {
    f3dHistGenSelectedPtVsYVsMultV0AXiPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultV0AXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultV0AXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYVsMultV0AOmegaMinus) {
    f3dHistGenSelectedPtVsYVsMultV0AOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultV0AOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultV0AOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultV0AOmegaPlus) {
    f3dHistGenSelectedPtVsYVsMultV0AOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultV0AOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultV0AOmegaPlus);
  }
  
  //ANALYSIS level Cascades, YCMS
  
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultV0AXiMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultV0AXiMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultV0AXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultV0AXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultV0AXiPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultV0AXiPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultV0AXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultV0AXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaPlus);
  }
  
    //Cross-check primary selection criteria
    //Should not impact much but inspired by rumors of DPMJet issues with primary selection (?)
    
    if(! f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiMinus) {
        f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiMinus = new TH3F( "f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
        if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiMinus);
    }
    if(! f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiPlus) {
        f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiPlus = new TH3F( "f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
        if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiPlus);
    }
    //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
    
    if(! f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaMinus) {
        f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaMinus = new TH3F( "f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
        if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaMinus);
    }
    if(! f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaPlus) {
        f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaPlus = new TH3F( "f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
        if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaPlus);
    }
  
  //======================================================================================
  //--------------------------------------------------------------------------------------
  // ZNA
  //--------------------------------------------------------------------------------------
  //======================================================================================
  
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, all generated
  
  if(! f3dHistGenPtVsYVsMultZNAXiMinus) {
    f3dHistGenPtVsYVsMultZNAXiMinus = new TH3F( "f3dHistGenPtVsYVsMultZNAXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultZNAXiMinus);
  }
  if(! f3dHistGenPtVsYVsMultZNAXiPlus) {
    f3dHistGenPtVsYVsMultZNAXiPlus = new TH3F( "f3dHistGenPtVsYVsMultZNAXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultZNAXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYVsMultZNAOmegaMinus) {
    f3dHistGenPtVsYVsMultZNAOmegaMinus = new TH3F( "f3dHistGenPtVsYVsMultZNAOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultZNAOmegaMinus);
  }
  if(! f3dHistGenPtVsYVsMultZNAOmegaPlus) {
    f3dHistGenPtVsYVsMultZNAOmegaPlus = new TH3F( "f3dHistGenPtVsYVsMultZNAOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultZNAOmegaPlus);
  }
  
  //All generated cascades, YCMS
  
  if(! f3dHistGenPtVsYCMSVsMultZNAXiMinus) {
    f3dHistGenPtVsYCMSVsMultZNAXiMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultZNAXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultZNAXiMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultZNAXiPlus) {
    f3dHistGenPtVsYCMSVsMultZNAXiPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultZNAXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultZNAXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYCMSVsMultZNAOmegaMinus) {
    f3dHistGenPtVsYCMSVsMultZNAOmegaMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultZNAOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultZNAOmegaMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultZNAOmegaPlus) {
    f3dHistGenPtVsYCMSVsMultZNAOmegaPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultZNAOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultZNAOmegaPlus);
  }
  
  
  //--------------------------------------------------------------------------------------
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, at selected analysis evts
  
  if(! f3dHistGenSelectedPtVsYVsMultZNAXiMinus) {
    f3dHistGenSelectedPtVsYVsMultZNAXiMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultZNAXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultZNAXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultZNAXiPlus) {
    f3dHistGenSelectedPtVsYVsMultZNAXiPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultZNAXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultZNAXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYVsMultZNAOmegaMinus) {
    f3dHistGenSelectedPtVsYVsMultZNAOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultZNAOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultZNAOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultZNAOmegaPlus) {
    f3dHistGenSelectedPtVsYVsMultZNAOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultZNAOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultZNAOmegaPlus);
  }
  
  //ANALYSIS level Cascades, YCMS
  
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultZNAXiMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultZNAXiMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultZNAXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultZNAXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultZNAXiPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultZNAXiPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultZNAXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultZNAXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaPlus);
  }

  
  //======================================================================================
  //--------------------------------------------------------------------------------------
  // TRK
  //--------------------------------------------------------------------------------------
  //======================================================================================
  
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, all generated
  
  if(! f3dHistGenPtVsYVsMultTRKXiMinus) {
    f3dHistGenPtVsYVsMultTRKXiMinus = new TH3F( "f3dHistGenPtVsYVsMultTRKXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultTRKXiMinus);
  }
  if(! f3dHistGenPtVsYVsMultTRKXiPlus) {
    f3dHistGenPtVsYVsMultTRKXiPlus = new TH3F( "f3dHistGenPtVsYVsMultTRKXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultTRKXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYVsMultTRKOmegaMinus) {
    f3dHistGenPtVsYVsMultTRKOmegaMinus = new TH3F( "f3dHistGenPtVsYVsMultTRKOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultTRKOmegaMinus);
  }
  if(! f3dHistGenPtVsYVsMultTRKOmegaPlus) {
    f3dHistGenPtVsYVsMultTRKOmegaPlus = new TH3F( "f3dHistGenPtVsYVsMultTRKOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultTRKOmegaPlus);
  }
  
  //All generated cascades, YCMS
  
  if(! f3dHistGenPtVsYCMSVsMultTRKXiMinus) {
    f3dHistGenPtVsYCMSVsMultTRKXiMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultTRKXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultTRKXiMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultTRKXiPlus) {
    f3dHistGenPtVsYCMSVsMultTRKXiPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultTRKXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultTRKXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYCMSVsMultTRKOmegaMinus) {
    f3dHistGenPtVsYCMSVsMultTRKOmegaMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultTRKOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultTRKOmegaMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultTRKOmegaPlus) {
    f3dHistGenPtVsYCMSVsMultTRKOmegaPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultTRKOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultTRKOmegaPlus);
  }
  
  
  //--------------------------------------------------------------------------------------
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, at selected analysis evts
  
  if(! f3dHistGenSelectedPtVsYVsMultTRKXiMinus) {
    f3dHistGenSelectedPtVsYVsMultTRKXiMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultTRKXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultTRKXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultTRKXiPlus) {
    f3dHistGenSelectedPtVsYVsMultTRKXiPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultTRKXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultTRKXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYVsMultTRKOmegaMinus) {
    f3dHistGenSelectedPtVsYVsMultTRKOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultTRKOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultTRKOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultTRKOmegaPlus) {
    f3dHistGenSelectedPtVsYVsMultTRKOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultTRKOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultTRKOmegaPlus);
  }
  
  //ANALYSIS level Cascades, YCMS
  
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultTRKXiMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultTRKXiMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultTRKXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultTRKXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultTRKXiPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultTRKXiPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultTRKXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultTRKXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaPlus);
  }
  
  //======================================================================================
  //--------------------------------------------------------------------------------------
  // SPD
  //--------------------------------------------------------------------------------------
  //======================================================================================
  
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, all generated
  
  if(! f3dHistGenPtVsYVsMultSPDXiMinus) {
    f3dHistGenPtVsYVsMultSPDXiMinus = new TH3F( "f3dHistGenPtVsYVsMultSPDXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultSPDXiMinus);
  }
  if(! f3dHistGenPtVsYVsMultSPDXiPlus) {
    f3dHistGenPtVsYVsMultSPDXiPlus = new TH3F( "f3dHistGenPtVsYVsMultSPDXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultSPDXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYVsMultSPDOmegaMinus) {
    f3dHistGenPtVsYVsMultSPDOmegaMinus = new TH3F( "f3dHistGenPtVsYVsMultSPDOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultSPDOmegaMinus);
  }
  if(! f3dHistGenPtVsYVsMultSPDOmegaPlus) {
    f3dHistGenPtVsYVsMultSPDOmegaPlus = new TH3F( "f3dHistGenPtVsYVsMultSPDOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYVsMultSPDOmegaPlus);
  }
  
  //All generated cascades, YCMS
  
  if(! f3dHistGenPtVsYCMSVsMultSPDXiMinus) {
    f3dHistGenPtVsYCMSVsMultSPDXiMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultSPDXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultSPDXiMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultSPDXiPlus) {
    f3dHistGenPtVsYCMSVsMultSPDXiPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultSPDXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultSPDXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenPtVsYCMSVsMultSPDOmegaMinus) {
    f3dHistGenPtVsYCMSVsMultSPDOmegaMinus = new TH3F( "f3dHistGenPtVsYCMSVsMultSPDOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultSPDOmegaMinus);
  }
  if(! f3dHistGenPtVsYCMSVsMultSPDOmegaPlus) {
    f3dHistGenPtVsYCMSVsMultSPDOmegaPlus = new TH3F( "f3dHistGenPtVsYCMSVsMultSPDOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenPtVsYCMSVsMultSPDOmegaPlus);
  }
  
  
  //--------------------------------------------------------------------------------------
  //--- 3D Histo (Pt, Y, Multiplicity) for generated XiMinus/Plus, at selected analysis evts
  
  if(! f3dHistGenSelectedPtVsYVsMultSPDXiMinus) {
    f3dHistGenSelectedPtVsYVsMultSPDXiMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultSPDXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultSPDXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultSPDXiPlus) {
    f3dHistGenSelectedPtVsYVsMultSPDXiPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultSPDXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultSPDXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYVsMultSPDOmegaMinus) {
    f3dHistGenSelectedPtVsYVsMultSPDOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYVsMultSPDOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultSPDOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYVsMultSPDOmegaPlus) {
    f3dHistGenSelectedPtVsYVsMultSPDOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYVsMultSPDOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYVsMultSPDOmegaPlus);
  }
  
  //ANALYSIS level Cascades, YCMS
  
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultSPDXiMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultSPDXiMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultSPDXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultSPDXiMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultSPDXiPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultSPDXiPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultSPDXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultSPDXiPlus);
  }
  //--- 3D Histo (Pt, Y, Multiplicity) for generated OmegaMinus/Plus
  
  if(! f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaMinus) {
    f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaMinus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaMinus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaMinus);
  }
  if(! f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaPlus) {
    f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaPlus = new TH3F( "f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaPlus", "Pt_{#Omega} Vs Y_{#Omega} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Omega} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
    if(fkIsNuclear) fListHist->Add(f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaPlus);
  }
  
//----------------------------------
// Primary Vertex Position Histos
//----------------------------------

   if(! fHistPVx) {
         fHistPVx = new TH1F("fHistPVx", 
            "PV x position;Nbr of Evts;x", 
            2000, -0.5, 0.5);       
      fListHist->Add(fHistPVx);
   }
   if(! fHistPVy) {
         fHistPVy = new TH1F("fHistPVy", 
            "PV y position;Nbr of Evts;y", 
            2000, -0.5, 0.5);       
      fListHist->Add(fHistPVy);
   }
   if(! fHistPVz) {
         fHistPVz = new TH1F("fHistPVz", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHist->Add(fHistPVz);
   }

   if(! fHistPVxAnalysis) {
         fHistPVxAnalysis = new TH1F("fHistPVxAnalysis", 
            "PV x position;Nbr of Evts;x", 
            2000, -0.5, 0.5);       
      fListHist->Add(fHistPVxAnalysis);
   }
   if(! fHistPVyAnalysis) {
         fHistPVyAnalysis = new TH1F("fHistPVyAnalysis", 
            "PV y position;Nbr of Evts;y", 
            2000, -0.5, 0.5);       
      fListHist->Add(fHistPVyAnalysis);
   }
   if(! fHistPVzAnalysis) {
         fHistPVzAnalysis = new TH1F("fHistPVzAnalysis", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHist->Add(fHistPVzAnalysis);
   }

   //List of Histograms: Normal
   PostData(1, fListHist);

   //TTree Object: Saved to base directory. Should cache to disk while saving. 
   //(Important to avoid excessive memory usage, particularly when merging)
   PostData(2, fTreeCascade);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskExtractPerformanceCascade::UserExec(Option_t *) 
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

/* --- Acquisition of exact event ID
   fTreeVariableRunNumber = lESDevent->GetRunNumber();
   fTreeVariableEventNumber =  
    ( ( ((ULong64_t)lESDevent->GetPeriodNumber() ) << 36 ) |
      ( ((ULong64_t)lESDevent->GetOrbitNumber () ) << 12 ) |
        ((ULong64_t)lESDevent->GetBunchCrossNumber() )  );
*/
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
   if(fkIsNuclear == kFALSE) lMultiplicity =  fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,  fEtaRefMult );

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

  
   //---> If this is a nuclear collision, then go nuclear on "multiplicity" variable...
   //---> Warning: Experimental
   if(fkIsNuclear == kTRUE){ 
      AliCentrality* centrality;
      centrality = lESDevent->GetCentrality();
      lMultiplicity = ( ( Int_t ) ( centrality->GetCentralityPercentile( fCentralityEstimator.Data() ) ) );
      lMultiplicityV0A = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "V0A" ) ) );
      lMultiplicityZNA = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "ZNA" ) ) );
      lMultiplicityTRK = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "TRK" ) ) );
      lMultiplicitySPD = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "CL1" ) ) );
      if (centrality->GetQuality()>1) {
        PostData(1, fListHist);
        PostData(2, fTreeCascade);
        return;
      }
   }
    
    if( fkSelectPeripheral ){
        if( lMultiplicity < 60 || lMultiplicity >= 80. ){
            //Event is outside 60-80% centrality in V0M!
            PostData(1, fListHist);
            PostData(2, fTreeCascade);
            return;
        }
    }
    
  
   //Set variable for filling tree afterwards!
   //---> pp case......: GetReferenceMultiplicity
   //---> Pb-Pb case...: Centrality by V0M

   fTreeCascVarMultiplicity = lMultiplicity;
   fTreeCascVarMultiplicityMC = lNumberOfCharged;
   fTreeCascVarMultiplicityV0A = lMultiplicityV0A;
   fTreeCascVarMultiplicityZNA = lMultiplicityZNA;
   fTreeCascVarMultiplicityTRK = lMultiplicityTRK;
   fTreeCascVarMultiplicitySPD = lMultiplicitySPD;

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

//------------------------------------------------
// Variable Definition
//------------------------------------------------

   Int_t lNbMCPrimary        = 0;

   Double_t lPtCurrentPart   = 0;
  
   //Int_t lComeFromSigma      = 0;

   // current mc particle 's mother
   //Int_t iCurrentMother  = 0;
   lNbMCPrimary = lMCstack->GetNprimary();

//------------------------------------------------
// Pre-Physics Selection
//------------------------------------------------


  fTreeCascVarEvHasXiMinus    = kFALSE; 
  fTreeCascVarEvHasXiPlus     = kFALSE; 
  fTreeCascVarEvHasOmegaMinus = kFALSE; 
  fTreeCascVarEvHasOmegaPlus  = kFALSE; 
  fTreeCascVarEvHasLambda     = kFALSE; 
  fTreeCascVarEvHasAntiLambda = kFALSE; 

  fTreeCascVarEvHasLowPtXiMinus    = kFALSE; 
  fTreeCascVarEvHasLowPtXiPlus     = kFALSE; 
  fTreeCascVarEvHasLowPtOmegaMinus = kFALSE; 
  fTreeCascVarEvHasLowPtOmegaPlus  = kFALSE; 
  fTreeCascVarEvHasLowPtLambda     = kFALSE; 
  fTreeCascVarEvHasLowPtAntiLambda = kFALSE; 

  fTreeCascVarEvHasVeryLowPtXiMinus    = kFALSE; 
  fTreeCascVarEvHasVeryLowPtXiPlus     = kFALSE; 
  fTreeCascVarEvHasVeryLowPtOmegaMinus = kFALSE; 
  fTreeCascVarEvHasVeryLowPtOmegaPlus  = kFALSE; 
  fTreeCascVarEvHasVeryLowPtLambda     = kFALSE; 
  fTreeCascVarEvHasVeryLowPtAntiLambda = kFALSE; 


//----- Loop on primary Xi, Omega --------------------------------------------------------------
   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < lNbMCPrimary; iCurrentLabelStack++) 
   {// This is the begining of the loop on primaries
      
      TParticle* lCurrentParticlePrimary = 0x0; 
      lCurrentParticlePrimary = lMCstack->Particle( iCurrentLabelStack );
      if(!lCurrentParticlePrimary){
         Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
         continue;
      }

      //Event Type Acquisition
      if(lCurrentParticlePrimary->GetPdgCode() == 3312){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) ) fTreeCascVarEvHasXiMinus     = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() ==-3312){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) ) fTreeCascVarEvHasXiPlus      = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() == 3334){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) ) fTreeCascVarEvHasOmegaMinus  = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() ==-3334){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) ) fTreeCascVarEvHasOmegaPlus   = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() == 3122){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) ) fTreeCascVarEvHasLambda      = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() ==-3122){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) ) fTreeCascVarEvHasAntiLambda  = kTRUE; }

      if(lCurrentParticlePrimary->GetPdgCode() == 3312){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<2.0) fTreeCascVarEvHasLowPtXiMinus     = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() ==-3312){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<2.0) fTreeCascVarEvHasLowPtXiPlus      = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() == 3334){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<2.0) fTreeCascVarEvHasLowPtOmegaMinus  = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() ==-3334){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<2.0) fTreeCascVarEvHasLowPtOmegaPlus   = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() == 3122){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<2.0) fTreeCascVarEvHasLowPtLambda      = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() ==-3122){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<2.0) fTreeCascVarEvHasLowPtAntiLambda  = kTRUE; }

      if(lCurrentParticlePrimary->GetPdgCode() == 3312){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<1.2) fTreeCascVarEvHasVeryLowPtXiMinus     = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() ==-3312){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<1.2) fTreeCascVarEvHasVeryLowPtXiPlus      = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() == 3334){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<1.2) fTreeCascVarEvHasVeryLowPtOmegaMinus  = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() ==-3334){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<1.2) fTreeCascVarEvHasVeryLowPtOmegaPlus   = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() == 3122){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<1.2) fTreeCascVarEvHasVeryLowPtLambda      = kTRUE; }
      if(lCurrentParticlePrimary->GetPdgCode() ==-3122){ if (TMath::Abs(lCurrentParticlePrimary->Eta())<1.2 && lMCstack->IsPhysicalPrimary(iCurrentLabelStack) && lCurrentParticlePrimary->Pt()<1.2) fTreeCascVarEvHasVeryLowPtAntiLambda  = kTRUE; }

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
           //MultMC
           f3dHistGenPtVsYVsMultMCXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lNumberOfCharged);
           f3dHistGenPtVsYCMSVsMultMCXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lNumberOfCharged);
           //V0A
           f3dHistGenPtVsYVsMultV0AXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityV0A);
           f3dHistGenPtVsYCMSVsMultV0AXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
           //ZNA
           f3dHistGenPtVsYVsMultZNAXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityZNA);
           f3dHistGenPtVsYCMSVsMultZNAXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityZNA);
           //TRK
           f3dHistGenPtVsYVsMultTRKXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityTRK);
           f3dHistGenPtVsYCMSVsMultTRKXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityTRK);
           //SPD
           f3dHistGenPtVsYVsMultSPDXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicitySPD);
           f3dHistGenPtVsYCMSVsMultSPDXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicitySPD);
         }
         if( lCurrentParticlePrimary->GetPdgCode() == -3312 ){ 
           lPtCurrentPart    = lCurrentParticlePrimary->Pt();
           f3dHistGenPtVsYVsMultXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
           f3dHistGenPtVsYCMSVsMultXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
           //MultMC
           f3dHistGenPtVsYVsMultMCXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lNumberOfCharged);
           f3dHistGenPtVsYCMSVsMultMCXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lNumberOfCharged);
           //V0A
           f3dHistGenPtVsYVsMultV0AXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityV0A);
           f3dHistGenPtVsYCMSVsMultV0AXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
           //ZNA
           f3dHistGenPtVsYVsMultZNAXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityZNA);
           f3dHistGenPtVsYCMSVsMultZNAXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityZNA);
           //TRK
           f3dHistGenPtVsYVsMultTRKXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityTRK);
           f3dHistGenPtVsYCMSVsMultTRKXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityTRK);
           //SPD
           f3dHistGenPtVsYVsMultSPDXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicitySPD);
           f3dHistGenPtVsYCMSVsMultSPDXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicitySPD);
         }
         // Omega Histograms
         if( lCurrentParticlePrimary->GetPdgCode() == 3334 ){ 
           lPtCurrentPart    = lCurrentParticlePrimary->Pt();
           f3dHistGenPtVsYVsMultOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
           f3dHistGenPtVsYCMSVsMultOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
            //MultMC
           f3dHistGenPtVsYVsMultMCOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lNumberOfCharged);
           f3dHistGenPtVsYCMSVsMultMCOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lNumberOfCharged);
           //V0A
           f3dHistGenPtVsYVsMultV0AOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityV0A);
           f3dHistGenPtVsYCMSVsMultV0AOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
           //ZNA
           f3dHistGenPtVsYVsMultZNAOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityZNA);
           f3dHistGenPtVsYCMSVsMultZNAOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityZNA);
           //TRK
           f3dHistGenPtVsYVsMultTRKOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityTRK);
           f3dHistGenPtVsYCMSVsMultTRKOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityTRK);
           //SPD
           f3dHistGenPtVsYVsMultSPDOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicitySPD);
           f3dHistGenPtVsYCMSVsMultSPDOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicitySPD);
         }
         if( lCurrentParticlePrimary->GetPdgCode() == -3334 ){
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenPtVsYVsMultOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenPtVsYCMSVsMultOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
            //MultMC
            f3dHistGenPtVsYVsMultMCOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lNumberOfCharged);
            f3dHistGenPtVsYCMSVsMultMCOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lNumberOfCharged);
           //V0A
           f3dHistGenPtVsYVsMultV0AOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityV0A);
           f3dHistGenPtVsYCMSVsMultV0AOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
           //ZNA
           f3dHistGenPtVsYVsMultZNAOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityZNA);
           f3dHistGenPtVsYCMSVsMultZNAOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityZNA);
           //TRK
           f3dHistGenPtVsYVsMultTRKOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityTRK);
           f3dHistGenPtVsYCMSVsMultTRKOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityTRK);
           //SPD
           f3dHistGenPtVsYVsMultSPDOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicitySPD);
           f3dHistGenPtVsYCMSVsMultSPDOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicitySPD);
         }
      }
   }
//----- End Loop on primary Xi, Omega ----------------------------------------------------------

// ---> Set Variables to Zero again
// ---> Variable Definition

   lPtCurrentPart   = 0;

//------------------------------------------------
// Physics Selection
//------------------------------------------------

   UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
   Bool_t isSelected = 0;
   isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;

   //pA triggering: CINT7
   if( fkSwitchINT7 ) isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;

   //Standard Min-Bias Selection
   if ( ! isSelected ) { 
        PostData(1, fListHist);
        PostData(2, fTreeCascade);
      return;
   }

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

  //------------------------------------------------
  // Primary Vertex Requirements Section:
  //  ---> pp and PbPb: Only requires |z|<10cm
  //  ---> pPb: all requirements checked at this stage
  //------------------------------------------------
  
  //Roberto's PV selection criteria, implemented 17th April 2013
  
  /* vertex selection */
  Bool_t fHasVertex = kFALSE;
  
  const AliESDVertex *vertex = lESDevent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1) {
    vertex = lESDevent->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1) fHasVertex = kFALSE;
    else fHasVertex = kTRUE;
    TString vtxTyp = vertex->GetTitle();
    Double_t cov[6]={0};
    vertex->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) fHasVertex = kFALSE;
  }
  else fHasVertex = kTRUE;
  
  //Is First event in chunk rejection: Still present!
  if(fkpAVertexSelection==kTRUE && fHasVertex == kFALSE) {
    AliWarning("Pb / | PV does not satisfy selection criteria!");
    PostData(1, fListHist);
    PostData(2, fTreeCascade);
    return;
  }
  
  //Is First event in chunk rejection: Still present!
  if(fkpAVertexSelection==kTRUE && fUtils->IsFirstEventInChunk(lESDevent)) {
    AliWarning("Pb / | This is the first event in the chunk!");
    PostData(1, fListHist);
    PostData(2, fTreeCascade);
    return;
  }
  
  //17 April Fix: Always do primary vertex Z selection, after pA vertex selection from Roberto
  if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0) {
    AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
    PostData(1, fListHist);
    PostData(2, fTreeCascade);
    return;
  }
  
  lMagneticField = lESDevent->GetMagneticField( );
  fHistV0MultiplicityForSelEvt ->Fill( lNumberOfV0s );
  fHistMultiplicity->Fill(lMultiplicity);
  fHistMultiplicityV0A->Fill(lMultiplicityV0A);
  fHistMultiplicityZNA->Fill(lMultiplicityZNA);
  fHistMultiplicityTRK->Fill(lMultiplicityTRK);
  fHistMultiplicitySPD->Fill(lMultiplicitySPD);

//------------------------------------------------
// SKIP: Events with well-established PVtx
//------------------------------------------------
	
   const AliESDVertex *lPrimaryTrackingESDVtxCheck = lESDevent->GetPrimaryVertexTracks();
   const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
   if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtxCheck->GetStatus() && fkpAVertexSelection==kFALSE ){
      AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
        PostData(1, fListHist);
        PostData(2, fTreeCascade);
      return;
   }
   fHistV0MultiplicityForSelEvtNoTPCOnly ->Fill( lNumberOfV0s );
   fHistMultiplicityNoTPCOnly->Fill(lMultiplicity);
  fHistMultiplicityV0ANoTPCOnly->Fill(lMultiplicityV0A);
  fHistMultiplicityZNANoTPCOnly->Fill(lMultiplicityZNA);
  fHistMultiplicityTRKNoTPCOnly->Fill(lMultiplicityTRK);
  fHistMultiplicitySPDNoTPCOnly->Fill(lMultiplicitySPD);

//------------------------------------------------
// Pileup Rejection Studies
//------------------------------------------------

   // FIXME : quality selection regarding pile-up rejection 
   if(lESDevent->IsPileupFromSPD() && !fkIsNuclear ){// minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  -> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
      AliWarning("Pb / This is tagged as Pileup from SPD... return !");
        PostData(1, fListHist);
        PostData(2, fTreeCascade);
      return;
   }
   fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup ->Fill( lNumberOfV0s );
   fHistMultiplicityNoTPCOnlyNoPileup->Fill(lMultiplicity);
  fHistMultiplicityV0ANoTPCOnlyNoPileup->Fill(lMultiplicityV0A);
  fHistMultiplicityZNANoTPCOnlyNoPileup->Fill(lMultiplicityZNA);
  fHistMultiplicityTRKNoTPCOnlyNoPileup->Fill(lMultiplicityTRK);
  fHistMultiplicitySPDNoTPCOnlyNoPileup->Fill(lMultiplicitySPD);

   //Do control histograms without the IsFromVertexerZ events, but consider them in analysis...
   if( ! (lESDevent->GetPrimaryVertex()->IsFromVertexerZ() )	 ){ 
      fHistPVxAnalysis->Fill( lPrimaryVtxPosition[0] );
      fHistPVyAnalysis->Fill( lPrimaryVtxPosition[1] );
      fHistPVzAnalysis->Fill( lPrimaryVtxPosition[2] );
   }

//------------------------------------------------
// stack loop starts here
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
            f3dHistGenSelectedPtVsYVsMultXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenSelectedPtVsYCMSVsMultXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
            //MultMC
            f3dHistGenSelectedPtVsYVsMultMCXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lNumberOfCharged);
            f3dHistGenSelectedPtVsYCMSVsMultMCXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lNumberOfCharged);
           //V0A
           f3dHistGenSelectedPtVsYVsMultV0AXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityV0A);
           f3dHistGenSelectedPtVsYCMSVsMultV0AXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
           //ZNA
           f3dHistGenSelectedPtVsYVsMultZNAXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityZNA);
           f3dHistGenSelectedPtVsYCMSVsMultZNAXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityZNA);
           //TRK
           f3dHistGenSelectedPtVsYVsMultTRKXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityTRK);
           f3dHistGenSelectedPtVsYCMSVsMultTRKXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityTRK);
           //SPD
           f3dHistGenSelectedPtVsYVsMultSPDXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicitySPD);
           f3dHistGenSelectedPtVsYCMSVsMultSPDXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicitySPD);
         }
         if( lCurrentParticlePrimary->GetPdgCode() == -3312 ){
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenSelectedPtVsYVsMultXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenSelectedPtVsYCMSVsMultXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
            //MultMC
            f3dHistGenSelectedPtVsYVsMultMCXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lNumberOfCharged);
            f3dHistGenSelectedPtVsYCMSVsMultMCXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lNumberOfCharged);
           //V0A
           f3dHistGenSelectedPtVsYVsMultV0AXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityV0A);
           f3dHistGenSelectedPtVsYCMSVsMultV0AXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
           //ZNA
           f3dHistGenSelectedPtVsYVsMultZNAXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityZNA);
           f3dHistGenSelectedPtVsYCMSVsMultZNAXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityZNA);
           //TRK
           f3dHistGenSelectedPtVsYVsMultTRKXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityTRK);
           f3dHistGenSelectedPtVsYCMSVsMultTRKXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityTRK);
           //SPD
           f3dHistGenSelectedPtVsYVsMultSPDXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicitySPD);
           f3dHistGenSelectedPtVsYCMSVsMultSPDXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicitySPD);
         }
         // Omega Histograms
         if( lCurrentParticlePrimary->GetPdgCode() == 3334 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenSelectedPtVsYVsMultOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenSelectedPtVsYCMSVsMultOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
            //MultMC
            f3dHistGenSelectedPtVsYVsMultMCOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lNumberOfCharged);
            f3dHistGenSelectedPtVsYCMSVsMultMCOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lNumberOfCharged);
           //V0A
           f3dHistGenSelectedPtVsYVsMultV0AOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityV0A);
           f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
           //ZNA
           f3dHistGenSelectedPtVsYVsMultZNAOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityZNA);
           f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityZNA);
           //TRK
           f3dHistGenSelectedPtVsYVsMultTRKOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityTRK);
           f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityTRK);
           //SPD
           f3dHistGenSelectedPtVsYVsMultSPDOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicitySPD);
           f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicitySPD);
         }
         if( lCurrentParticlePrimary->GetPdgCode() == -3334 ){
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenSelectedPtVsYVsMultOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
            f3dHistGenSelectedPtVsYCMSVsMultOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicity);
            //MultMC
            f3dHistGenSelectedPtVsYVsMultMCOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lNumberOfCharged);
            f3dHistGenSelectedPtVsYCMSVsMultMCOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lNumberOfCharged);
           //V0A
           f3dHistGenSelectedPtVsYVsMultV0AOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityV0A);
           f3dHistGenSelectedPtVsYCMSVsMultV0AOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
           //ZNA
           f3dHistGenSelectedPtVsYVsMultZNAOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityZNA);
           f3dHistGenSelectedPtVsYCMSVsMultZNAOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityZNA);
           //TRK
           f3dHistGenSelectedPtVsYVsMultTRKOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicityTRK);
           f3dHistGenSelectedPtVsYCMSVsMultTRKOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityTRK);
           //SPD
           f3dHistGenSelectedPtVsYVsMultSPDOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicitySPD);
           f3dHistGenSelectedPtVsYCMSVsMultSPDOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicitySPD);
         }
      }
   }
//----- End Loop on primary Xi, Omega ----------------------------------------------------------

    
    //----- Loop on primary Xi, Omega --------------------------------------------------------------
    for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
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
            //IsPhysicalPrimary: Cross-check with old index test
            // (cascades -> negligible difference ?... )
            if ( lMCstack->IsPhysicalPrimary(iCurrentLabelStack)!=kTRUE ) continue;
            
            //=================================================================================
            // Xi Histograms
            if( lCurrentParticlePrimary->GetPdgCode() == 3312 ){
                lPtCurrentPart    = lCurrentParticlePrimary->Pt();
                f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
            }
            if( lCurrentParticlePrimary->GetPdgCode() == -3312 ){
                lPtCurrentPart    = lCurrentParticlePrimary->Pt();
                f3dHistGenSelectedPrimPtVsYCMSVsMultV0AXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
            }
            // Omega Histograms
            if( lCurrentParticlePrimary->GetPdgCode() == 3334 ){
                lPtCurrentPart    = lCurrentParticlePrimary->Pt();
                f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaMinus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
            }
            if( lCurrentParticlePrimary->GetPdgCode() == -3334 ){
                lPtCurrentPart    = lCurrentParticlePrimary->Pt();
                f3dHistGenSelectedPrimPtVsYCMSVsMultV0AOmegaPlus->Fill(lPtCurrentPart, lRapXiMCPrimary+fpArapidityShift, lMultiplicityV0A);
            }
        }
    }
    //----- End Loop on primary Xi, Omega ----------------------------------------------------------

    
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
	  Double_t lXiTransvMomMC= 0. ;
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
	
	  Double_t lRapXi   = -20.0, lRapOmega = -20.0, lRapMC = -20.0; //  lEta = -20.0, lTheta = 360., lPhi = 720. ;
	  //Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
	    
    // -------------------------------------
    // II.ESD - Calculation Part dedicated to Xi vertices (ESD)
    
	  AliESDcascade *xi = lESDevent->GetCascade(iXi);
	  if (!xi) continue;
	
          
                  // - II.Step 1 : around primary vertex
                  //-------------
          //lTrkgPrimaryVtxRadius3D = TMath::Sqrt(  lTrkgPrimaryVtxPos[0] * lTrkgPrimaryVtxPos[0] +
          //                                        lTrkgPrimaryVtxPos[1] * lTrkgPrimaryVtxPos[1] +
          //                                        lTrkgPrimaryVtxPos[2] * lTrkgPrimaryVtxPos[2] );

          //lBestPrimaryVtxRadius3D = TMath::Sqrt(  lBestPrimaryVtxPos[0] * lBestPrimaryVtxPos[0] +
          //                                        lBestPrimaryVtxPos[1] * lBestPrimaryVtxPos[1] +
          //                                        lBestPrimaryVtxPos[2] * lBestPrimaryVtxPos[2] );

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
      
            //Save shared clusters information
      fTreeCascVarNegSharedClusters = nTrackXi->GetTPCnclsS(0,159);
      fTreeCascVarPosSharedClusters = pTrackXi->GetTPCnclsS(0,159);
      fTreeCascVarBachSharedClusters = bachTrackXi->GetTPCnclsS(0,159);
      
      Double_t lBMom[3], lNMom[3], lPMom[3];
      xi->GetBPxPyPz( lBMom[0], lBMom[1], lBMom[2] );
      xi->GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
      xi->GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );
      
      fTreeCascVarBachTransMom = TMath::Sqrt( lBMom[0]*lBMom[0] + lBMom[1]*lBMom[1] );
      fTreeCascVarPosTransMom  = TMath::Sqrt( lPMom[0]*lPMom[0] + lPMom[1]*lPMom[1] );
      fTreeCascVarNegTransMom  = TMath::Sqrt( lNMom[0]*lNMom[0] + lNMom[1]*lNMom[1] );
      
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

      fTreeCascVarNegClusters = lNegTPCClusters;
      fTreeCascVarPosClusters = lPosTPCClusters;
      fTreeCascVarBachClusters = lBachTPCClusters;

    fTreeCascVarkITSRefitBachelor = kTRUE; 
    fTreeCascVarkITSRefitNegative = kTRUE; 
    fTreeCascVarkITSRefitPositive = kTRUE; 

    // 1 - Poor quality related to TPCrefit
	  ULong_t pStatus    = pTrackXi->GetStatus();
	  ULong_t nStatus    = nTrackXi->GetStatus();
	  ULong_t bachStatus = bachTrackXi->GetStatus();
    if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }
    if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; }
    if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }

    //Extra Debug Information: booleans for ITS refit
    if ((pStatus&AliESDtrack::kITSrefit)    == 0) { fTreeCascVarkITSRefitPositive = kFALSE; }
    if ((nStatus&AliESDtrack::kITSrefit)    == 0) { fTreeCascVarkITSRefitNegative = kFALSE; }
    if ((bachStatus&AliESDtrack::kITSrefit) == 0) { fTreeCascVarkITSRefitBachelor = kFALSE; }

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
    //Modification to check cosPA of v0 wrt to cascade vertex
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
// Associate Cascade Candidates to Monte Carlo!
//------------------------------------------------	

//Warning: Not using Continues... Need to fill tree later!
	
	Int_t lPDGCodeCascade = 0;	

	Int_t lPID_BachMother = 0;
	Int_t lPID_NegMother = 0;
	Int_t lPID_PosMother = 0;
  fTreeCascVarIsPhysicalPrimary = 0; // 0: not defined, any candidate may have this

	  fTreeCascVarPIDPositive = 0;
	  fTreeCascVarPIDNegative = 0;
	  fTreeCascVarPIDBachelor = 0;


	if(fDebug > 5)
		cout 	<< "MC EventNumber : " << lMCevent->Header()->GetEvent() 
			<< " / MC event Number in Run : " << lMCevent->Header()->GetEventNrInRun() << endl;
	
	// - Step 4.1 : level of the V0 daughters
		
//----------------------------------------
// Regular MC ASSOCIATION STARTS HERE
//----------------------------------------

	  Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrackXi->GetLabel() );  
		  // Abs value = needed ! question of quality track association ...
	  Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
	  Int_t lblBach        = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );

	  TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
	  TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
	  TParticle* mcBach        = lMCstack->Particle( lblBach );
      
    fTreeCascVarPosTransMomMC = mcPosV0Dghter->Pt();
    fTreeCascVarNegTransMomMC = mcNegV0Dghter->Pt();

	  fTreeCascVarPIDPositive = mcPosV0Dghter -> GetPdgCode();
	  fTreeCascVarPIDNegative = mcNegV0Dghter -> GetPdgCode();
	  fTreeCascVarPIDBachelor = mcBach->GetPdgCode();

	  // - Step 4.2 : level of the Xi daughters
		
	  Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ; 
	  Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
	
	  //Rather uncivilized: Open brackets for each 'continue'
	  if(! (lblMotherPosV0Dghter != lblMotherNegV0Dghter) ) { // same mother
	  if(! (lblMotherPosV0Dghter < 0) ) { // mother != primary (!= -1)
	  if(! (lblMotherNegV0Dghter < 0) ) {
					
		// mothers = Lambda candidate ... a priori
	
	  TParticle* mcMotherPosV0Dghter = lMCstack->Particle( lblMotherPosV0Dghter );
	  TParticle* mcMotherNegV0Dghter = lMCstack->Particle( lblMotherNegV0Dghter );
			
	  // - Step 4.3 : level of Xi candidate
	
	  Int_t lblGdMotherPosV0Dghter =   mcMotherPosV0Dghter->GetFirstMother() ;
	  Int_t lblGdMotherNegV0Dghter =   mcMotherNegV0Dghter->GetFirstMother() ;
				
		if(! (lblGdMotherPosV0Dghter != lblGdMotherNegV0Dghter) ) {
		if(! (lblGdMotherPosV0Dghter < 0) ) { // primary lambda ...
		if(! (lblGdMotherNegV0Dghter < 0) ) { // primary lambda ...

		  // Gd mothers = Xi candidate ... a priori
	
	  TParticle* mcGdMotherPosV0Dghter = lMCstack->Particle( lblGdMotherPosV0Dghter );
	  TParticle* mcGdMotherNegV0Dghter = lMCstack->Particle( lblGdMotherNegV0Dghter );
					
	  Int_t lblMotherBach = (Int_t) TMath::Abs( mcBach->GetFirstMother()  );
	
  //		if( lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters
		  if(!(lblMotherBach != lblGdMotherPosV0Dghter)) { //same mother for bach and V0 daughters
	
	  TParticle* mcMotherBach = lMCstack->Particle( lblMotherBach );
	
    lPID_BachMother = mcMotherBach->GetPdgCode();
	  lPID_NegMother = mcGdMotherPosV0Dghter->GetPdgCode();
	  lPID_PosMother = mcGdMotherNegV0Dghter->GetPdgCode();
   
	  if(lPID_BachMother==lPID_NegMother && lPID_BachMother==lPID_PosMother){ 
		  lPDGCodeCascade = lPID_BachMother; 
      lXiTransvMomMC = mcMotherBach->Pt();
      if( lMCstack->IsPhysicalPrimary       (lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 1; //Is Primary!
      if( lMCstack->IsSecondaryFromWeakDecay(lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 2; //Weak Decay!
      if( lMCstack->IsSecondaryFromMaterial (lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 3; //From Material!
      if ( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) !=0 ){
        lRapMC = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
      }
	  }

  }}}}}}} //Ends all conditionals above...

//----------------------------------------
// Regular MC ASSOCIATION ENDS HERE
//----------------------------------------

//----------------------------------------
// Swapped MC Association Starts Here 
//----------------------------------------

  fTreeCascVarPIDSwapped = 0; //Nothing
  if ( fkCheckSwapping ) { //on/off switch if not needed or otherwise broken

	  Int_t lPID_BachMotherSwapped = 0;
	  Int_t lPID_NegMotherSwapped = 0;
  	Int_t lPID_PosMotherSwapped = 0;

	  Int_t lblPosV0DghterSwapped = (Int_t) TMath::Abs( pTrackXi->GetLabel() );  
		  // Abs value = needed ! question of quality track association ...
	  Int_t lblNegV0DghterSwapped = (Int_t) TMath::Abs( nTrackXi->GetLabel() );
	  Int_t lblBachSwapped        = (Int_t) TMath::Abs( bachTrackXi->GetLabel() );

    //This is SWAPPED association: swap like-sign particles
    TParticle* mcPosV0DghterSwapped = 0x0;
    TParticle* mcNegV0DghterSwapped = 0x0;
    TParticle* mcBachSwapped        = 0x0;

    //Swapping Case 1: XiMinus/OmegaMinus swapped: exchange negative V0 daughter with bachelor
    if( lChargeXi == -1 ){ 
      mcPosV0DghterSwapped = lMCstack->Particle( lblPosV0DghterSwapped );
      mcNegV0DghterSwapped = lMCstack->Particle( lblBachSwapped        );
      mcBachSwapped        = lMCstack->Particle( lblNegV0DghterSwapped );	
    }
    //Swapping Case 2: XiPlus/OmegaPlus swapped: exchange positive V0 daughter with bachelor
    if( lChargeXi ==  1 ){ 
      mcPosV0DghterSwapped = lMCstack->Particle( lblBachSwapped        );
      mcNegV0DghterSwapped = lMCstack->Particle( lblNegV0DghterSwapped );
      mcBachSwapped        = lMCstack->Particle( lblPosV0DghterSwapped );	
    }
	
    //fTreeCascVarPosTransMomMC = mcPosV0Dghter->Pt();
    //fTreeCascVarNegTransMomMC = mcNegV0Dghter->Pt();

	  //fTreeCascVarPIDPositive = mcPosV0Dghter -> GetPdgCode();
	  //fTreeCascVarPIDNegative = mcNegV0Dghter -> GetPdgCode();
	  //fTreeCascVarPIDBachelor = mcBach->GetPdgCode();

	  // - Step 4.2 : level of the Xi daughters
		
	  Int_t lblMotherPosV0DghterSwapped = mcPosV0DghterSwapped->GetFirstMother() ; 
	  Int_t lblMotherNegV0DghterSwapped = mcNegV0DghterSwapped->GetFirstMother();
	
	  //Rather uncivilized: Open brackets for each 'continue'
	  if(! (lblMotherPosV0DghterSwapped != lblMotherNegV0DghterSwapped) ) { // same mother
	  if(! (lblMotherPosV0DghterSwapped < 0) ) { // mother != primary (!= -1)
	  if(! (lblMotherNegV0DghterSwapped < 0) ) {
					
		// mothers = Lambda candidate ... a priori
	
	  TParticle* mcMotherPosV0DghterSwapped = lMCstack->Particle( lblMotherPosV0DghterSwapped );
	  TParticle* mcMotherNegV0DghterSwapped = lMCstack->Particle( lblMotherNegV0DghterSwapped );
			
	  // - Step 4.3 : level of Xi candidate
	
	  Int_t lblGdMotherPosV0DghterSwapped =   mcMotherPosV0DghterSwapped->GetFirstMother() ;
	  Int_t lblGdMotherNegV0DghterSwapped =   mcMotherNegV0DghterSwapped->GetFirstMother() ;
				
		if(! (lblGdMotherPosV0DghterSwapped != lblGdMotherNegV0DghterSwapped) ) {
		if(! (lblGdMotherPosV0DghterSwapped < 0) ) { // primary lambda ...
		if(! (lblGdMotherNegV0DghterSwapped < 0) ) { // primary lambda ...

		  // Gd mothers = Xi candidate ... a priori
	
	  TParticle* mcGdMotherPosV0DghterSwapped = lMCstack->Particle( lblGdMotherPosV0DghterSwapped );
	  TParticle* mcGdMotherNegV0DghterSwapped = lMCstack->Particle( lblGdMotherNegV0DghterSwapped );
					
	  Int_t lblMotherBachSwapped = (Int_t) TMath::Abs( mcBachSwapped->GetFirstMother()  );
	
  //		if( lblMotherBach != lblGdMotherPosV0Dghter ) continue; //same mother for bach and V0 daughters
		  if(!(lblMotherBachSwapped != lblGdMotherPosV0DghterSwapped)) { //same mother for bach and V0 daughters
	
	  TParticle* mcMotherBachSwapped = lMCstack->Particle( lblMotherBachSwapped );
	
    lPID_BachMotherSwapped = mcMotherBachSwapped->GetPdgCode();
	  lPID_NegMotherSwapped = mcGdMotherPosV0DghterSwapped->GetPdgCode();
	  lPID_PosMotherSwapped = mcGdMotherNegV0DghterSwapped->GetPdgCode();
   
	  if(lPID_BachMotherSwapped==lPID_NegMotherSwapped && lPID_BachMotherSwapped==lPID_PosMotherSwapped){ 
		  fTreeCascVarPIDSwapped = lPID_BachMotherSwapped; //there!
	  }

  }}}}}}} //Ends all conditionals above...

  }
//----------------------------------------
// Swapped MC ASSOCIATION ENDS HERE
//----------------------------------------

  //------------------------------------------------
  // Jai Salzwedel's femto-cut: better V0 exists
  //------------------------------------------------			  

  fTreeCascVarBadCascadeJai = kFALSE; 
  fTreeCascVarDeltaDCA = -100;
  Float_t DCAV0DaughtersDiff = -100; 
  for (Int_t iv0=0; iv0<lESDevent->GetNumberOfV0s(); iv0++) {
    AliESDv0 *v0 = lESDevent->GetV0(iv0);
    UInt_t posV0TrackIdx = (UInt_t) v0->GetPindex();
    UInt_t negV0TrackIdx = (UInt_t) v0->GetNindex();
    if ((posV0TrackIdx == lIdxPosXi) && (negV0TrackIdx == lIdxNegXi)) continue;
      // if both tracks are the same ones as the cascades V0 daughter tracks, then the V0 belongs to the cascade being analysed; so avoid it
    if ((posV0TrackIdx == lIdxPosXi) || (negV0TrackIdx == lIdxNegXi)) {
    DCAV0DaughtersDiff = lDcaV0DaughtersXi - v0->GetDcaV0Daughters();
    if( fTreeCascVarDeltaDCA < DCAV0DaughtersDiff ) fTreeCascVarDeltaDCA = DCAV0DaughtersDiff;
      if ( lDcaV0DaughtersXi > v0->GetDcaV0Daughters() )  {    // DCA comparison criterion
        fTreeCascVarBadCascadeJai = kTRUE;
      } //end DCA comparison 
    } // end shares a daughter check 
  } //end V0 loop 

  //------------------------------------------------
  // Set Variables for adding to tree
  //------------------------------------------------		
	
/* 1*/		fTreeCascVarCharge	= lChargeXi;
/* 2*/		if(lInvMassXiMinus!=0)    fTreeCascVarMassAsXi = lInvMassXiMinus;
/* 2*/		if(lInvMassXiPlus!=0)     fTreeCascVarMassAsXi = lInvMassXiPlus;
/* 3*/		if(lInvMassOmegaMinus!=0) fTreeCascVarMassAsOmega = lInvMassOmegaMinus;
/* 3*/		if(lInvMassOmegaPlus!=0)  fTreeCascVarMassAsOmega = lInvMassOmegaPlus;
/* 4*/		fTreeCascVarPt = lXiTransvMom;
/* 4*/		fTreeCascVarPtMC = lXiTransvMomMC;
/* 5*/		fTreeCascVarRapXi = lRapXi ;
/* 5*/		fTreeCascVarRapMC = lRapMC ;
/* 6*/		fTreeCascVarRapOmega = lRapOmega ;
/* 7*/		fTreeCascVarDCACascDaughters = lDcaXiDaughters;
/* 8*/		fTreeCascVarDCABachToPrimVtx = lDcaBachToPrimVertexXi;
/* 9*/		fTreeCascVarDCAV0Daughters = lDcaV0DaughtersXi;
/*10*/		fTreeCascVarDCAV0ToPrimVtx = lDcaV0ToPrimVertexXi;
/*11*/		fTreeCascVarDCAPosToPrimVtx = lDcaPosToPrimVertexXi;
/*12*/		fTreeCascVarDCANegToPrimVtx = lDcaNegToPrimVertexXi;
/*13*/		fTreeCascVarCascCosPointingAngle = lXiCosineOfPointingAngle;
/*14*/		fTreeCascVarCascRadius = lXiRadius;
/*15*/		fTreeCascVarV0Mass = lInvMassLambdaAsCascDghter;
/*16*/		fTreeCascVarV0CosPointingAngle = lV0CosineOfPointingAngleXi;
/*16*/		fTreeCascVarV0CosPointingAngleSpecial = lV0CosineOfPointingAngleXiSpecial;
/*17*/		fTreeCascVarV0Radius = lV0RadiusXi;
/*20*/		fTreeCascVarLeastNbrClusters = leastnumberofclusters;
/*21*/		fTreeCascVarMultiplicity = lMultiplicity; //multiplicity, whatever that may be

/*23*/		fTreeCascVarDistOverTotMom = TMath::Sqrt(
						TMath::Power( lPosXi[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( lPosXi[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( lPosXi[2] - lBestPrimaryVtxPos[2] , 2)
					);
/*23*/		fTreeCascVarDistOverTotMom /= (lXiTotMom+1e-13);
/*24*/    //Not specified here, it has been set already (TRunNumber)
/*25*/		fTreeCascVarPID = lPDGCodeCascade;

//------------------------------------------------
// Fill Tree! 
//------------------------------------------------

// The conditional is meant to decrease excessive
// memory usage! Be careful when loosening the 
// cut!

  //Xi    Mass window: 150MeV wide
  //Omega mass window: 150MeV wide

  if( (fTreeCascVarMassAsXi<1.32+0.075&&fTreeCascVarMassAsXi>1.32-0.075) ||
      (fTreeCascVarMassAsOmega<1.68+0.075&&fTreeCascVarMassAsOmega>1.68-0.075) ){
      fTreeCascade->Fill();
  }

//------------------------------------------------
// Fill tree over.
//------------------------------------------------

	}// end of the Cascade loop (ESD or AOD)

   // Post output data.
   PostData(1, fListHist);
   PostData(2, fTreeCascade);
}

//________________________________________________________________________
void AliAnalysisTaskExtractPerformanceCascade::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   TList *cRetrievedList = 0x0;
   cRetrievedList = (TList*)GetOutputData(1);
   if(!cRetrievedList){
      Printf("ERROR - AliAnalysisTaskExtractCascade : ouput data container list not available\n");
      return;
   }	
	
   fHistV0MultiplicityForTrigEvt = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistV0MultiplicityForTrigEvt")  );
   if (!fHistV0MultiplicityForTrigEvt) {
      Printf("ERROR - AliAnalysisTaskExtractCascade : fHistV0MultiplicityForTrigEvt not available");
      return;
   }
  
   TCanvas *canCheck = new TCanvas("AliAnalysisTaskExtractCascade","V0 Multiplicity",10,10,510,510);
   canCheck->cd(1)->SetLogy();

   fHistV0MultiplicityForTrigEvt->SetMarkerStyle(22);
   fHistV0MultiplicityForTrigEvt->DrawCopy("E");
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskExtractPerformanceCascade::MyRapidity(Double_t rE, Double_t rPz) const
{
   // Local calculation for rapidity
   Double_t ReturnValue = -100;
   if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ){ 
      ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
   }
   return ReturnValue;
} 
