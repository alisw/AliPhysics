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
#include "AliPPVsMultUtils.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

#include "AliPWG0Helper.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskStrangenessVsMultiplicityMC.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskStrangenessVsMultiplicityMC)

AliAnalysisTaskStrangenessVsMultiplicityMC::AliAnalysisTaskStrangenessVsMultiplicityMC()
    : AliAnalysisTaskSE(), fListHist(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fPPVsMultUtils(0), fUtils(0),
      fkSaveV0Tree      ( kFALSE ),
      fkSaveCascadeTree ( kTRUE  ),
      fkMCAssociation   ( kTRUE  ),
      fkSaveLambda( kTRUE ),
      fkSaveAntiLambda( kTRUE ),
      fkSaveK0Short( kTRUE ),
      fkSaveExtendedRefMultInfo( kFALSE ),
      fkSelectTriggerByName ( kFALSE ) ,
      fkRunVertexers    ( kTRUE  ),
      fkSkipEventSelection( kFALSE ),
      fkApplyTrackletsVsClustersCut( kFALSE ),
      fMinbAnalysis(kFALSE),
      fkMultSelection ( kFALSE ),
      fTrigType(AliVEvent::kMB),
      fTrigName(""),
      //---> Variables for fTreeEvent
      fAmplitude_V0A   (0),
      fAmplitude_V0C   (0),
      fAmplitude_V0M   (0),
      fAmplitude_V0AEq (0),
      fAmplitude_V0CEq (0),
      fAmplitude_V0MEq (0),
      fCentrality_V0A(0),
      fCentrality_V0C(0),
      fCentrality_V0M(0),
fCentrality_OnlineV0A(0),
fCentrality_OnlineV0C(0),
fCentrality_OnlineV0M(0),
fCentrality_ADA(0),
fCentrality_ADC(0),
fCentrality_ADM(0),
      fCentrality_V0AEq(0),
      fCentrality_V0CEq(0),
      fCentrality_V0MEq(0),
      fCentrality_V0B(0),
      fCentrality_V0Apartial(0),
      fCentrality_V0Cpartial(0),
      fCentrality_V0S(0),
      fCentrality_V0SB(0),
      fRefMultEta5(0),
      fRefMultEta8(0),
      fTrueMultEta5(0),
      fTrueMultEta8(0),
      fTrueMultEta10(0),
      fTrueMultVZEROA(0),
      fTrueMultVZEROC(0),
      fRunNumber(0),
      fEvSel_nTrackletsEta10(0),
      fEvSel_HasAtLeastSPDVertex(0),
      fEvSel_VtxZCut(0),
      fEvSel_IsNotPileup(0),
      fEvSel_IsNotPileupMV(0),
      fEvSel_IsNotPileupInMultBins(0),
      fEvSel_HasVtxContributor(0),
      fEvSel_Triggered(0),
      fEvSel_INELgtZERO(0),
      fEvSel_INELgtZEROStackPrimaries(0),
      fEvSel_INELgtZEROtracklets(0),
      fEvSel_INELgtZERORefMult(0),
      fEvSel_INELgtZERORefMultTracklets(0),
      fEvSel_VtxZ(0),
      fEvSel_VtxZMC(0),
      fEvSel_MCType(0),
      //---> Variables for fTreeV0
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
      fTreeVariableNegEta(0),
      fTreeVariablePosEta(0),

      fTreeVariableNSigmasPosProton(0),
      fTreeVariableNSigmasPosPion(0),
      fTreeVariableNSigmasNegProton(0),
      fTreeVariableNSigmasNegPion(0),
      fTreeVariableNegTransvMomentumMC(0),
      fTreeVariablePosTransvMomentumMC(0),

      fTreeVariableDistOverTotMom(0),
      fTreeVariableLeastNbrCrossedRows(0),
      fTreeVariableLeastRatioCrossedRowsOverFindable(0),

      fTreeVariableCentV0M(0),
      fTreeVariableCentV0A(0),
      fTreeVariableCentV0C(0),
fTreeVariableCentOnlineV0M(0),
fTreeVariableCentOnlineV0A(0),
fTreeVariableCentOnlineV0C(0),
fTreeVariableCentADM(0),
fTreeVariableCentADA(0),
fTreeVariableCentADC(0),
      fTreeVariableCentV0MEq(0),
      fTreeVariableCentV0AEq(0),
      fTreeVariableCentV0CEq(0),
      fTreeVariableCentV0B(0),
      fTreeVariableCentV0Apartial(0),
      fTreeVariableCentV0Cpartial(0),
      fTreeVariableCentV0S(0),
      fTreeVariableCentV0SB(0),

      fTreeVariableRefMultEta8(0),
      fTreeVariableRefMultEta5(0),
      fTreeVariableRunNumber(0),

      fTreeVariablePtMother(0),
fTreeVariableRapMother(0),
      fTreeVariablePID(0),
      fTreeVariablePIDPositive(0),
      fTreeVariablePIDNegative(0),
      fTreeVariablePIDMother(0),
      fTreeVariablePrimaryStatus(0),
      fTreeVariablePrimaryStatusMother(0),
      //---> Variables for fTreeCascade
      fTreeCascVarCharge(0),
      fTreeCascVarMassAsXi(0),
      fTreeCascVarMassAsOmega(0),
      fTreeCascVarPt(0),
      fTreeCascVarPtMC(0),
      fTreeCascVarRapXi(0),
      fTreeCascVarRapOmega(0),
      fTreeCascVarRapMC(0),
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
      fTreeCascVarNegTransvMomentumMC(0),
      fTreeCascVarPosTransvMomentumMC(0),
      fTreeCascVarCentV0M(0),
      fTreeCascVarCentV0A(0),
      fTreeCascVarCentV0C(0),
fTreeCascVarCentOnlineV0M(0),
fTreeCascVarCentOnlineV0A(0),
fTreeCascVarCentOnlineV0C(0),
fTreeCascVarCentADM(0),
fTreeCascVarCentADA(0),
fTreeCascVarCentADC(0),
      fTreeCascVarCentV0MEq(0),
      fTreeCascVarCentV0AEq(0),
      fTreeCascVarCentV0CEq(0),
      fTreeCascVarCentV0B(0),
      fTreeCascVarCentV0Apartial(0),
      fTreeCascVarCentV0Cpartial(0),
      fTreeCascVarCentV0S(0),
      fTreeCascVarCentV0SB(0),
      fTreeCascVarRefMultEta8(0),
      fTreeCascVarRefMultEta5(0),
      fTreeCascVarTrueMultEta5(0),
      fTreeCascVarTrueMultEta8(0),
      fTreeCascVarTrueMultVZEROA(0),
      fTreeCascVarTrueMultVZEROC(0),
      fTreeCascVarIsPhysicalPrimary(0),
      fTreeCascVarPID(0),
      fTreeCascVarRunNumber(0),
      //---> Histograms
      fHistEventCounter(0),
      //---> MC Generated Histo (analysis level)
      fHistPt_GenK0Short(0),
      fHistPt_GenLambda(0),
      fHistPt_GenAntiLambda(0),
      fHistPt_GenXiMinus(0),
      fHistPt_GenXiPlus(0),
      fHistPt_GenOmegaMinus(0),
      fHistPt_GenOmegaPlus(0),

      //VsRefMult
      fHistPtVsRefMultEta5_GenK0Short(0),
      fHistPtVsRefMultEta5_GenLambda(0),
      fHistPtVsRefMultEta5_GenAntiLambda(0),
      fHistPtVsRefMultEta5_GenXiMinus(0),
      fHistPtVsRefMultEta5_GenXiPlus(0),
      fHistPtVsRefMultEta5_GenOmegaMinus(0),
      fHistPtVsRefMultEta5_GenOmegaPlus(0),
      fHistPtVsRefMultEta8_GenK0Short(0),
      fHistPtVsRefMultEta8_GenLambda(0),
      fHistPtVsRefMultEta8_GenAntiLambda(0),
      fHistPtVsRefMultEta8_GenXiMinus(0),
      fHistPtVsRefMultEta8_GenXiPlus(0),
      fHistPtVsRefMultEta8_GenOmegaMinus(0),
      fHistPtVsRefMultEta8_GenOmegaPlus(0),

      //VsCentralities
      fHistPtVsCentV0A_GenK0Short(0),
      fHistPtVsCentV0A_GenLambda(0),
      fHistPtVsCentV0A_GenAntiLambda(0),
      fHistPtVsCentV0A_GenXiMinus(0),
      fHistPtVsCentV0A_GenXiPlus(0),
      fHistPtVsCentV0A_GenOmegaMinus(0),
      fHistPtVsCentV0A_GenOmegaPlus(0),
      fHistPtVsCentV0C_GenK0Short(0),
      fHistPtVsCentV0C_GenLambda(0),
      fHistPtVsCentV0C_GenAntiLambda(0),
      fHistPtVsCentV0C_GenXiMinus(0),
      fHistPtVsCentV0C_GenXiPlus(0),
      fHistPtVsCentV0C_GenOmegaMinus(0),
      fHistPtVsCentV0C_GenOmegaPlus(0),
      fHistPtVsCentV0M_GenK0Short(0),
      fHistPtVsCentV0M_GenLambda(0),
      fHistPtVsCentV0M_GenAntiLambda(0),
      fHistPtVsCentV0M_GenXiMinus(0),
      fHistPtVsCentV0M_GenXiPlus(0),
      fHistPtVsCentV0M_GenOmegaMinus(0),
      fHistPtVsCentV0M_GenOmegaPlus(0),
      fHistPtVsCentV0MUnselected_GenK0Short(0),
      fHistPtVsCentV0MUnselected_GenLambda(0),
      fHistPtVsCentV0MUnselected_GenAntiLambda(0),
      fHistPtVsCentV0MUnselected_GenXiMinus(0),
      fHistPtVsCentV0MUnselected_GenXiPlus(0),
      fHistPtVsCentV0MUnselected_GenOmegaMinus(0),
      fHistPtVsCentV0MUnselected_GenOmegaPlus(0),

      
      //Equalized
      fHistPtVsCentV0AEq_GenXiMinus(0),
      fHistPtVsCentV0AEq_GenXiPlus(0),
      fHistPtVsCentV0AEq_GenOmegaMinus(0),
      fHistPtVsCentV0AEq_GenOmegaPlus(0),
      fHistPtVsCentV0CEq_GenXiMinus(0),
      fHistPtVsCentV0CEq_GenXiPlus(0),
      fHistPtVsCentV0CEq_GenOmegaMinus(0),
      fHistPtVsCentV0CEq_GenOmegaPlus(0),
      fHistPtVsCentV0MEq_GenXiMinus(0),
      fHistPtVsCentV0MEq_GenXiPlus(0),
      fHistPtVsCentV0MEq_GenOmegaMinus(0),
      fHistPtVsCentV0MEq_GenOmegaPlus(0),

      //VsAmp
      fHistPtVsAmpV0A_GenXiMinus(0),
      fHistPtVsAmpV0A_GenXiPlus(0),
      fHistPtVsAmpV0A_GenOmegaMinus(0),
      fHistPtVsAmpV0A_GenOmegaPlus(0),
      fHistPtVsAmpV0C_GenXiMinus(0),
      fHistPtVsAmpV0C_GenXiPlus(0),
      fHistPtVsAmpV0C_GenOmegaMinus(0),
      fHistPtVsAmpV0C_GenOmegaPlus(0),
      fHistPtVsAmpV0M_GenXiMinus(0),
      fHistPtVsAmpV0M_GenXiPlus(0),
      fHistPtVsAmpV0M_GenOmegaMinus(0),
      fHistPtVsAmpV0M_GenOmegaPlus(0),
      //Equalized Amps
      fHistPtVsAmpV0AEq_GenXiMinus(0),
      fHistPtVsAmpV0AEq_GenXiPlus(0),
      fHistPtVsAmpV0AEq_GenOmegaMinus(0),
      fHistPtVsAmpV0AEq_GenOmegaPlus(0),
      fHistPtVsAmpV0CEq_GenXiMinus(0),
      fHistPtVsAmpV0CEq_GenXiPlus(0),
      fHistPtVsAmpV0CEq_GenOmegaMinus(0),
      fHistPtVsAmpV0CEq_GenOmegaPlus(0),
      fHistPtVsAmpV0MEq_GenXiMinus(0),
      fHistPtVsAmpV0MEq_GenXiPlus(0),
      fHistPtVsAmpV0MEq_GenOmegaMinus(0),
      fHistPtVsAmpV0MEq_GenOmegaPlus(0),
      fHistVZEROResponseStudy(0),
      fHistVZEROResponseStudyTotal(0)

//------------------------------------------------
// Tree Variables
{

}

AliAnalysisTaskStrangenessVsMultiplicityMC::AliAnalysisTaskStrangenessVsMultiplicityMC(const char *name)
    : AliAnalysisTaskSE(name), fListHist(0), fTreeEvent(0), fTreeV0(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fPPVsMultUtils(0), fUtils(0),
      fkSaveV0Tree      ( kFALSE ),
      fkSaveCascadeTree ( kTRUE  ),
      fkMCAssociation   ( kTRUE  ),
      fkSaveLambda( kTRUE ),
      fkSaveAntiLambda( kTRUE ),
      fkSaveK0Short( kTRUE ),
      fkSaveExtendedRefMultInfo( kFALSE ),
      fkSelectTriggerByName ( kFALSE ) ,
      fkRunVertexers    ( kTRUE  ),
      fkSkipEventSelection( kFALSE ),
      fkApplyTrackletsVsClustersCut( kFALSE ),
      fMinbAnalysis(kFALSE),
      fkMultSelection ( kFALSE ),
      fTrigType(AliVEvent::kMB),
      fTrigName(""),
      //---> Variables for fTreeEvent
      fAmplitude_V0A (0),
      fAmplitude_V0C (0),
      fAmplitude_V0M (0),
      fAmplitude_V0AEq (0),
      fAmplitude_V0CEq (0),
      fAmplitude_V0MEq (0),
      fCentrality_V0A(0),
      fCentrality_V0C(0),
      fCentrality_V0M(0),
fCentrality_OnlineV0A(0),
fCentrality_OnlineV0C(0),
fCentrality_OnlineV0M(0),
fCentrality_ADA(0),
fCentrality_ADC(0),
fCentrality_ADM(0),
      fCentrality_V0AEq(0),
      fCentrality_V0CEq(0),
      fCentrality_V0MEq(0),
      fCentrality_V0B(0),
      fCentrality_V0Apartial(0),
      fCentrality_V0Cpartial(0),
      fCentrality_V0S(0),
      fCentrality_V0SB(0),
      fRefMultEta5(0),
      fRefMultEta8(0),
      fTrueMultEta5(0),
      fTrueMultEta8(0),
      fTrueMultEta10(0),
      fTrueMultVZEROA(0),
      fTrueMultVZEROC(0),
      fRunNumber(0),
      fEvSel_nTrackletsEta10(0),
      fEvSel_HasAtLeastSPDVertex(0),
      fEvSel_VtxZCut(0),
      fEvSel_IsNotPileup(0),
      fEvSel_IsNotPileupMV(0),
      fEvSel_IsNotPileupInMultBins(0),
      fEvSel_HasVtxContributor(0),
      fEvSel_Triggered(0),
      fEvSel_INELgtZERO(0),
      fEvSel_INELgtZEROStackPrimaries(0),
      fEvSel_INELgtZEROtracklets(0),
      fEvSel_INELgtZERORefMult(0),
      fEvSel_INELgtZERORefMultTracklets(0),
      fEvSel_VtxZ(0),
      fEvSel_VtxZMC(0),
      fEvSel_MCType(0),
      //---> Variables for fTreeV0
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
      fTreeVariableNegEta(0),
      fTreeVariablePosEta(0),

      fTreeVariableNSigmasPosProton(0),
      fTreeVariableNSigmasPosPion(0),
      fTreeVariableNSigmasNegProton(0),
      fTreeVariableNSigmasNegPion(0),
      fTreeVariableNegTransvMomentumMC(0),
      fTreeVariablePosTransvMomentumMC(0),

      fTreeVariableDistOverTotMom(0),
      fTreeVariableLeastNbrCrossedRows(0),
      fTreeVariableLeastRatioCrossedRowsOverFindable(0),

      fTreeVariableCentV0M(0),
      fTreeVariableCentV0A(0),
      fTreeVariableCentV0C(0),
fTreeVariableCentOnlineV0M(0),
fTreeVariableCentOnlineV0A(0),
fTreeVariableCentOnlineV0C(0),
fTreeVariableCentADM(0),
fTreeVariableCentADA(0),
fTreeVariableCentADC(0),
      fTreeVariableCentV0MEq(0),
      fTreeVariableCentV0AEq(0),
      fTreeVariableCentV0CEq(0),
      fTreeVariableCentV0B(0),
      fTreeVariableCentV0Apartial(0),
      fTreeVariableCentV0Cpartial(0),
      fTreeVariableCentV0S(0),
      fTreeVariableCentV0SB(0),

      fTreeVariableRefMultEta8(0),
      fTreeVariableRefMultEta5(0),
      fTreeVariableRunNumber(0),

      fTreeVariablePtMother(0),
fTreeVariableRapMother(0),
fTreeVariablePID(0),
      fTreeVariablePIDPositive(0),
      fTreeVariablePIDNegative(0),
      fTreeVariablePIDMother(0),
      fTreeVariablePrimaryStatus(0),
      fTreeVariablePrimaryStatusMother(0),
      //---> Variables for fTreeCascade
      fTreeCascVarCharge(0),
      fTreeCascVarMassAsXi(0),
      fTreeCascVarMassAsOmega(0),
      fTreeCascVarPt(0),
      fTreeCascVarPtMC(0),
      fTreeCascVarRapXi(0),
      fTreeCascVarRapOmega(0),
      fTreeCascVarRapMC(0),
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
      fTreeCascVarNegTransvMomentumMC(0),
      fTreeCascVarPosTransvMomentumMC(0),
      fTreeCascVarCentV0M(0),
      fTreeCascVarCentV0A(0),
      fTreeCascVarCentV0C(0),
fTreeCascVarCentOnlineV0M(0),
fTreeCascVarCentOnlineV0A(0),
fTreeCascVarCentOnlineV0C(0),
fTreeCascVarCentADM(0),
fTreeCascVarCentADA(0),
fTreeCascVarCentADC(0),
      fTreeCascVarCentV0MEq(0),
      fTreeCascVarCentV0AEq(0),
      fTreeCascVarCentV0CEq(0),
      fTreeCascVarCentV0B(0),
      fTreeCascVarCentV0Apartial(0),
      fTreeCascVarCentV0Cpartial(0),
      fTreeCascVarCentV0S(0),
      fTreeCascVarCentV0SB(0),
      fTreeCascVarRefMultEta8(0),
      fTreeCascVarRefMultEta5(0),
      fTreeCascVarTrueMultEta5(0),
      fTreeCascVarTrueMultEta8(0),
      fTreeCascVarTrueMultVZEROA(0),
      fTreeCascVarTrueMultVZEROC(0),
      fTreeCascVarIsPhysicalPrimary(0),
      fTreeCascVarPID(0),
      fTreeCascVarRunNumber(0),
      //---> Histograms
      fHistEventCounter(0),
      //---> MC Generated Histo (analysis level)
      fHistPt_GenK0Short(0),
      fHistPt_GenLambda(0),
      fHistPt_GenAntiLambda(0),
      fHistPt_GenXiMinus(0),
      fHistPt_GenXiPlus(0),
      fHistPt_GenOmegaMinus(0),
      fHistPt_GenOmegaPlus(0),

      //VsRefMult
      fHistPtVsRefMultEta5_GenK0Short(0),
      fHistPtVsRefMultEta5_GenLambda(0),
      fHistPtVsRefMultEta5_GenAntiLambda(0),
      fHistPtVsRefMultEta5_GenXiMinus(0),
      fHistPtVsRefMultEta5_GenXiPlus(0),
      fHistPtVsRefMultEta5_GenOmegaMinus(0),
      fHistPtVsRefMultEta5_GenOmegaPlus(0),
      fHistPtVsRefMultEta8_GenK0Short(0),
      fHistPtVsRefMultEta8_GenLambda(0),
      fHistPtVsRefMultEta8_GenAntiLambda(0),
      fHistPtVsRefMultEta8_GenXiMinus(0),
      fHistPtVsRefMultEta8_GenXiPlus(0),
      fHistPtVsRefMultEta8_GenOmegaMinus(0),
      fHistPtVsRefMultEta8_GenOmegaPlus(0),

      //VsCentralities
      fHistPtVsCentV0A_GenK0Short(0),
      fHistPtVsCentV0A_GenLambda(0),
      fHistPtVsCentV0A_GenAntiLambda(0),
      fHistPtVsCentV0A_GenXiMinus(0),
      fHistPtVsCentV0A_GenXiPlus(0),
      fHistPtVsCentV0A_GenOmegaMinus(0),
      fHistPtVsCentV0A_GenOmegaPlus(0),
      fHistPtVsCentV0C_GenK0Short(0),
      fHistPtVsCentV0C_GenLambda(0),
      fHistPtVsCentV0C_GenAntiLambda(0),
      fHistPtVsCentV0C_GenXiMinus(0),
      fHistPtVsCentV0C_GenXiPlus(0),
      fHistPtVsCentV0C_GenOmegaMinus(0),
      fHistPtVsCentV0C_GenOmegaPlus(0),
      fHistPtVsCentV0M_GenK0Short(0),
      fHistPtVsCentV0M_GenLambda(0),
      fHistPtVsCentV0M_GenAntiLambda(0),
      fHistPtVsCentV0M_GenXiMinus(0),
      fHistPtVsCentV0M_GenXiPlus(0),
      fHistPtVsCentV0M_GenOmegaMinus(0),
      fHistPtVsCentV0M_GenOmegaPlus(0),
      fHistPtVsCentV0MUnselected_GenK0Short(0),
      fHistPtVsCentV0MUnselected_GenLambda(0),
      fHistPtVsCentV0MUnselected_GenAntiLambda(0),
      fHistPtVsCentV0MUnselected_GenXiMinus(0),
      fHistPtVsCentV0MUnselected_GenXiPlus(0),
      fHistPtVsCentV0MUnselected_GenOmegaMinus(0),
      fHistPtVsCentV0MUnselected_GenOmegaPlus(0),
      
      //Equalized
      fHistPtVsCentV0AEq_GenXiMinus(0),
      fHistPtVsCentV0AEq_GenXiPlus(0),
      fHistPtVsCentV0AEq_GenOmegaMinus(0),
      fHistPtVsCentV0AEq_GenOmegaPlus(0),
      fHistPtVsCentV0CEq_GenXiMinus(0),
      fHistPtVsCentV0CEq_GenXiPlus(0),
      fHistPtVsCentV0CEq_GenOmegaMinus(0),
      fHistPtVsCentV0CEq_GenOmegaPlus(0),
      fHistPtVsCentV0MEq_GenXiMinus(0),
      fHistPtVsCentV0MEq_GenXiPlus(0),
      fHistPtVsCentV0MEq_GenOmegaMinus(0),
      fHistPtVsCentV0MEq_GenOmegaPlus(0),

      //VsAmp
      fHistPtVsAmpV0A_GenXiMinus(0),
      fHistPtVsAmpV0A_GenXiPlus(0),
      fHistPtVsAmpV0A_GenOmegaMinus(0),
      fHistPtVsAmpV0A_GenOmegaPlus(0),
      fHistPtVsAmpV0C_GenXiMinus(0),
      fHistPtVsAmpV0C_GenXiPlus(0),
      fHistPtVsAmpV0C_GenOmegaMinus(0),
      fHistPtVsAmpV0C_GenOmegaPlus(0),
      fHistPtVsAmpV0M_GenXiMinus(0),
      fHistPtVsAmpV0M_GenXiPlus(0),
      fHistPtVsAmpV0M_GenOmegaMinus(0),
      fHistPtVsAmpV0M_GenOmegaPlus(0),
      //Equalized Amps
      fHistPtVsAmpV0AEq_GenXiMinus(0),
      fHistPtVsAmpV0AEq_GenXiPlus(0),
      fHistPtVsAmpV0AEq_GenOmegaMinus(0),
      fHistPtVsAmpV0AEq_GenOmegaPlus(0),
      fHistPtVsAmpV0CEq_GenXiMinus(0),
      fHistPtVsAmpV0CEq_GenXiPlus(0),
      fHistPtVsAmpV0CEq_GenOmegaMinus(0),
      fHistPtVsAmpV0CEq_GenOmegaPlus(0),
      fHistPtVsAmpV0MEq_GenXiMinus(0),
      fHistPtVsAmpV0MEq_GenXiPlus(0),
      fHistPtVsAmpV0MEq_GenOmegaMinus(0),
      fHistPtVsAmpV0MEq_GenOmegaPlus(0),
      fHistVZEROResponseStudy(0),
      fHistVZEROResponseStudyTotal(0)
{

    //Re-vertex: Will only apply for cascade candidates

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

    for(int i=0; i<20; i++) {
        fTreeCascVarRefMultDiffEta[i] = 0;
        fTreeVariableRefMultDiffEta[i] = 0;
        fRefMultDiffEta[i] = 0;
    }

    DefineOutput(1, TList::Class()); // Event Counter Histo
    DefineOutput(2, TTree::Class()); // Event Tree
    DefineOutput(3, TTree::Class()); // V0 Tree
    DefineOutput(4, TTree::Class()); // Cascade Tree
}


AliAnalysisTaskStrangenessVsMultiplicityMC::~AliAnalysisTaskStrangenessVsMultiplicityMC()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------

    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
    if (fTreeEvent) {
        delete fTreeEvent;
        fTreeEvent = 0x0;
    }
    if (fTreeV0) {
        delete fTreeV0;
        fTreeV0 = 0x0;
    }
    if (fTreeCascade) {
        delete fTreeCascade;
        fTreeCascade = 0x0;
    }
    if (fPPVsMultUtils) {
        delete fPPVsMultUtils;
        fPPVsMultUtils = 0x0;
    }
    if (fUtils) {
        delete fUtils;
        fUtils = 0x0;
    }
}

//________________________________________________________________________
void AliAnalysisTaskStrangenessVsMultiplicityMC::UserCreateOutputObjects()
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
    fTreeEvent->Branch("fAmplitude_V0M",&fAmplitude_V0M,"fAmplitude_V0M/F");
    fTreeEvent->Branch("fAmplitude_V0AEq",&fAmplitude_V0AEq,"fAmplitude_V0AEq/F");
    fTreeEvent->Branch("fAmplitude_V0CEq",&fAmplitude_V0CEq,"fAmplitude_V0CEq/F");
    fTreeEvent->Branch("fAmplitude_V0MEq",&fAmplitude_V0MEq,"fAmplitude_V0MEq/F");

    //Info from AliCentrality (not necessarily 'centrality' per se)
    fTreeEvent->Branch("fCentrality_V0A",&fCentrality_V0A,"fCentrality_V0A/F");
    fTreeEvent->Branch("fCentrality_V0C",&fCentrality_V0C,"fCentrality_V0C/F");
    fTreeEvent->Branch("fCentrality_V0M",&fCentrality_V0M,"fCentrality_V0M/F");
    if( !fkMultSelection ){
        fTreeEvent->Branch("fCentrality_V0AEq",&fCentrality_V0AEq,"fCentrality_V0AEq/F");
        fTreeEvent->Branch("fCentrality_V0CEq",&fCentrality_V0CEq,"fCentrality_V0CEq/F");
        fTreeEvent->Branch("fCentrality_V0MEq",&fCentrality_V0MEq,"fCentrality_V0MEq/F");
        fTreeEvent->Branch("fCentrality_V0B",&fCentrality_V0B,"fCentrality_V0B/F");
        fTreeEvent->Branch("fCentrality_V0Apartial",&fCentrality_V0Apartial,"fCentrality_V0Apartial/F");
        fTreeEvent->Branch("fCentrality_V0Cpartial",&fCentrality_V0Cpartial,"fCentrality_V0Cpartial/F");
        fTreeEvent->Branch("fCentrality_V0S",&fCentrality_V0S,"fCentrality_V0S/F");
        fTreeEvent->Branch("fCentrality_V0SB",&fCentrality_V0SB,"fCentrality_V0SB/F");
    }else{
        fTreeEvent->Branch("fCentrality_OnlineV0A",&fCentrality_OnlineV0A,"fCentrality_OnlineV0A/F");
        fTreeEvent->Branch("fCentrality_OnlineV0C",&fCentrality_OnlineV0C,"fCentrality_OnlineV0C/F");
        fTreeEvent->Branch("fCentrality_OnlineV0M",&fCentrality_OnlineV0M,"fCentrality_OnlineV0M/F");
        fTreeEvent->Branch("fCentrality_ADA",&fCentrality_ADA,"fCentrality_ADA/F");
        fTreeEvent->Branch("fCentrality_ADC",&fCentrality_ADC,"fCentrality_ADC/F");
        fTreeEvent->Branch("fCentrality_ADM",&fCentrality_ADM,"fCentrality_ADM/F");
    }

    //Official GetReferenceMultiplicity
    fTreeEvent->Branch("fRefMultEta5",&fRefMultEta5,"fRefMultEta5/I");
    fTreeEvent->Branch("fRefMultEta8",&fRefMultEta8,"fRefMultEta8/I");
    //Don't do this if not explicitly requested, takes up too much space
    if ( fkSaveExtendedRefMultInfo )
        fTreeEvent->Branch("fRefMultDiffEta",fRefMultDiffEta,"fRefMultDiffEta[20]/I");

    fTreeEvent->Branch("fTrueMultEta5",&fTrueMultEta5,"fTrueMultEta5/I");
    fTreeEvent->Branch("fTrueMultEta8",&fTrueMultEta8,"fTrueMultEta8/I");
    fTreeEvent->Branch("fTrueMultEta10",&fTrueMultEta10,"fTrueMultEta10/I");
    fTreeEvent->Branch("fTrueMultVZEROA",&fTrueMultVZEROA,"fTrueMultVZEROA/I");
    fTreeEvent->Branch("fTrueMultVZEROC",&fTrueMultVZEROC,"fTrueMultVZEROC/I");

    //Run Number
    fTreeEvent->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
    fTreeEvent->Branch("fEvSel_nTrackletsEta10", &fEvSel_nTrackletsEta10, "fEvSel_nTrackletsEta10/I");

    //Booleans for Event Selection
    fTreeEvent->Branch("fEvSel_HasAtLeastSPDVertex", &fEvSel_HasAtLeastSPDVertex, "fEvSel_HasAtLeastSPDVertex/O");
    fTreeEvent->Branch("fEvSel_VtxZCut", &fEvSel_VtxZCut, "fEvSel_VtxZCut/O");
    fTreeEvent->Branch("fEvSel_IsNotPileup", &fEvSel_IsNotPileup, "fEvSel_IsNotPileup/O");
    fTreeEvent->Branch("fEvSel_IsNotPileupMV", &fEvSel_IsNotPileupMV, "fEvSel_IsNotPileupMV/O");
    fTreeEvent->Branch("fEvSel_IsNotPileupInMultBins", &fEvSel_IsNotPileupInMultBins, "fEvSel_IsNotPileupInMultBins/O");
    fTreeEvent->Branch("fEvSel_HasVtxContributor", &fEvSel_HasVtxContributor, "fEvSel_HasVtxContributor/O");
    fTreeEvent->Branch("fEvSel_Triggered", &fEvSel_Triggered, "fEvSel_Triggered/O");
    fTreeEvent->Branch("fEvSel_INELgtZERO", &fEvSel_INELgtZERO, "fEvSel_INELgtZERO/O");
    fTreeEvent->Branch("fEvSel_INELgtZEROStackPrimaries", &fEvSel_INELgtZEROStackPrimaries, "fEvSel_INELgtZEROStackPrimaries/O");
    fTreeEvent->Branch("fEvSel_INELgtZEROtracklets", &fEvSel_INELgtZEROtracklets, "fEvSel_INELgtZEROtracklets/O");
    fTreeEvent->Branch("fEvSel_INELgtZERORefMult", &fEvSel_INELgtZERORefMult, "fEvSel_INELgtZERORefMult/O");
    fTreeEvent->Branch("fEvSel_INELgtZERORefMultTracklets", &fEvSel_INELgtZERORefMultTracklets, "fEvSel_INELgtZERORefMultTracklets/O");

    fTreeEvent->Branch("fEvSel_VtxZ", &fEvSel_VtxZ, "fEvSel_VtxZ/F");
    fTreeEvent->Branch("fEvSel_VtxZMC", &fEvSel_VtxZMC, "fEvSel_VtxZMC/F");
    fTreeEvent->Branch("fEvSel_MCType", &fEvSel_MCType, "fEvSel_MCType/I");


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
    fTreeV0->Branch("fTreeVariablePtMC",&fTreeVariablePtMC,"fTreeVariablePtMC/F");
    fTreeV0->Branch("fTreeVariableRapK0Short",&fTreeVariableRapK0Short,"fTreeVariableRapK0Short/F");
    fTreeV0->Branch("fTreeVariableRapLambda",&fTreeVariableRapLambda,"fTreeVariableRapLambda/F");
    fTreeV0->Branch("fTreeVariableRapMC",&fTreeVariableRapMC,"fTreeVariableRapMC/F");
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

    fTreeV0->Branch("fTreeVariableNegTransvMomentumMC",&fTreeVariableNegTransvMomentumMC,"fTreeVariableNegTransvMomentumMC/F");
    fTreeV0->Branch("fTreeVariablePosTransvMomentumMC",&fTreeVariablePosTransvMomentumMC,"fTreeVariablePosTransvMomentumMC/F");

    fTreeV0->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
    fTreeV0->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
//-----------MULTIPLICITY-INFO--------------------
    fTreeV0->Branch("fTreeVariableCentV0M",&fTreeVariableCentV0M,"fTreeVariableCentV0M/F");
    fTreeV0->Branch("fTreeVariableCentV0A",&fTreeVariableCentV0A,"fTreeVariableCentV0A/F");
    fTreeV0->Branch("fTreeVariableCentV0C",&fTreeVariableCentV0C,"fTreeVariableCentV0C/F");
    if ( !fkMultSelection ){
        fTreeV0->Branch("fTreeVariableCentV0MEq",&fTreeVariableCentV0MEq,"fTreeVariableCentV0MEq/F");
        fTreeV0->Branch("fTreeVariableCentV0AEq",&fTreeVariableCentV0AEq,"fTreeVariableCentV0AEq/F");
        fTreeV0->Branch("fTreeVariableCentV0CEq",&fTreeVariableCentV0CEq,"fTreeVariableCentV0CEq/F");
        fTreeV0->Branch("fTreeVariableCentV0B",&fTreeVariableCentV0B,"fTreeVariableCentV0B/F");
        fTreeV0->Branch("fTreeVariableCentV0Apartial",&fTreeVariableCentV0Apartial,"fTreeVariableCentV0Apartial/F");
        fTreeV0->Branch("fTreeVariableCentV0Cpartial",&fTreeVariableCentV0Cpartial,"fTreeVariableCentV0Cpartial/F");
        fTreeV0->Branch("fTreeVariableCentV0S",&fTreeVariableCentV0S,"fTreeVariableCentV0S/F");
        fTreeV0->Branch("fTreeVariableCentV0SB",&fTreeVariableCentV0SB,"fTreeVariableCentV0SB/F");
    }else{
        fTreeV0->Branch("fTreeVariableCentOnlineV0M",&fTreeVariableCentOnlineV0M,"fTreeVariableCentOnlineV0M/F");
        fTreeV0->Branch("fTreeVariableCentOnlineV0A",&fTreeVariableCentOnlineV0A,"fTreeVariableCentOnlineV0A/F");
        fTreeV0->Branch("fTreeVariableCentOnlineV0C",&fTreeVariableCentOnlineV0C,"fTreeVariableCentOnlineV0C/F");
        fTreeV0->Branch("fTreeVariableCentADM",&fTreeVariableCentADM,"fTreeVariableCentADM/F");
        fTreeV0->Branch("fTreeVariableCentADA",&fTreeVariableCentADA,"fTreeVariableCentADA/F");
        fTreeV0->Branch("fTreeVariableCentADC",&fTreeVariableCentADC,"fTreeVariableCentADC/F");
    }
    fTreeV0->Branch("fTreeVariableRefMultEta8",&fTreeVariableRefMultEta8,"fTreeVariableRefMultEta8/I");
    fTreeV0->Branch("fTreeVariableRefMultEta5",&fTreeVariableRefMultEta5,"fTreeVariableRefMultEta5/I");
    //Don't do this if not explicitly requested, takes up too much space
    if ( fkSaveExtendedRefMultInfo )
        fTreeV0->Branch("fTreeVariableRefMultDiffEta",fTreeVariableRefMultDiffEta,"fTreeVariableRefMultDiffEta[20]/I");
    fTreeV0->Branch("fTreeVariableRunNumber",&fTreeVariableRunNumber,"fTreeVariableRunNumber/I");
//-----------MC Exclusive info--------------------
    fTreeV0->Branch("fTreeVariablePtMother",&fTreeVariablePtMother,"fTreeVariablePtMother/F");
    fTreeV0->Branch("fTreeVariableRapMother",&fTreeVariableRapMother,"fTreeVariableRapMother/F");
    fTreeV0->Branch("fTreeVariablePID",&fTreeVariablePID,"fTreeVariablePID/I");
    fTreeV0->Branch("fTreeVariablePIDPositive",&fTreeVariablePIDPositive,"fTreeVariablePIDPositive/I");
    fTreeV0->Branch("fTreeVariablePIDNegative",&fTreeVariablePIDNegative,"fTreeVariablePIDNegative/I");
    fTreeV0->Branch("fTreeVariablePIDMother",&fTreeVariablePIDMother,"fTreeVariablePIDMother/I");
    fTreeV0->Branch("fTreeVariablePrimaryStatus",&fTreeVariablePrimaryStatus,"fTreeVariablePrimaryStatus/I");
    fTreeV0->Branch("fTreeVariablePrimaryStatusMother",&fTreeVariablePrimaryStatusMother,"fTreeVariablePrimaryStatusMother/I");
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
    fTreeCascade->Branch("fTreeCascVarPtMC",&fTreeCascVarPtMC,"fTreeCascVarPtMC/F");
    fTreeCascade->Branch("fTreeCascVarRapXi",&fTreeCascVarRapXi,"fTreeCascVarRapXi/F");
    fTreeCascade->Branch("fTreeCascVarRapOmega",&fTreeCascVarRapOmega,"fTreeCascVarRapOmega/F");
    fTreeCascade->Branch("fTreeCascVarRapMC",&fTreeCascVarRapMC,"fTreeCascVarRapMC/F");
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
    fTreeCascade->Branch("fTreeCascVarCentV0M",&fTreeCascVarCentV0M,"fTreeCascVarCentV0M/F");
    fTreeCascade->Branch("fTreeCascVarCentV0A",&fTreeCascVarCentV0A,"fTreeCascVarCentV0A/F");
    fTreeCascade->Branch("fTreeCascVarCentV0C",&fTreeCascVarCentV0C,"fTreeCascVarCentV0C/F");
    if ( !fkMultSelection ){
        fTreeCascade->Branch("fTreeCascVarCentV0MEq",&fTreeCascVarCentV0MEq,"fTreeCascVarCentV0MEq/F");
        fTreeCascade->Branch("fTreeCascVarCentV0AEq",&fTreeCascVarCentV0AEq,"fTreeCascVarCentV0AEq/F");
        fTreeCascade->Branch("fTreeCascVarCentV0CEq",&fTreeCascVarCentV0CEq,"fTreeCascVarCentV0CEq/F");
        fTreeCascade->Branch("fTreeCascVarCentV0B",&fTreeCascVarCentV0B,"fTreeCascVarCentV0B/F");
        fTreeCascade->Branch("fTreeCascVarCentV0Apartial",&fTreeCascVarCentV0Apartial,"fTreeCascVarCentV0Apartial/F");
        fTreeCascade->Branch("fTreeCascVarCentV0Cpartial",&fTreeCascVarCentV0Cpartial,"fTreeCascVarCentV0Cpartial/F");
        fTreeCascade->Branch("fTreeCascVarCentV0S",&fTreeCascVarCentV0S,"fTreeCascVarCentV0S/F");
        fTreeCascade->Branch("fTreeCascVarCentV0SB",&fTreeCascVarCentV0SB,"fTreeCascVarCentV0SB/F");
    }else{
        fTreeCascade->Branch("fTreeCascVarCentOnlineV0M",&fTreeCascVarCentOnlineV0M,"fTreeCascVarCentOnlineV0M/F");
        fTreeCascade->Branch("fTreeCascVarCentOnlineV0A",&fTreeCascVarCentOnlineV0A,"fTreeCascVarCentOnlineV0A/F");
        fTreeCascade->Branch("fTreeCascVarCentOnlineV0C",&fTreeCascVarCentOnlineV0C,"fTreeCascVarCentOnlineV0C/F");
        fTreeCascade->Branch("fTreeCascVarCentADM",&fTreeCascVarCentADM,"fTreeCascVarCentADM/F");
        fTreeCascade->Branch("fTreeCascVarCentADA",&fTreeCascVarCentADA,"fTreeCascVarCentADA/F");
        fTreeCascade->Branch("fTreeCascVarCentADC",&fTreeCascVarCentADC,"fTreeCascVarCentADC/F");
    }
    fTreeCascade->Branch("fTreeCascVarRefMultEta8",&fTreeCascVarRefMultEta8,"fTreeCascVarRefMultEta8/I");
    fTreeCascade->Branch("fTreeCascVarRefMultEta5",&fTreeCascVarRefMultEta5,"fTreeCascVarRefMultEta5/I");
    //Don't do this if not explicitly requested, takes up too much space
    if ( fkSaveExtendedRefMultInfo )
        fTreeCascade->Branch("fTreeCascVarRefMultDiffEta",fTreeCascVarRefMultDiffEta,"fTreeCascVarRefMultDiffEta[20]/I");
    fTreeCascade->Branch("fTreeCascVarTrueMultEta5",&fTreeCascVarTrueMultEta5,"fTreeCascVarTrueMultEta5/I");
    fTreeCascade->Branch("fTreeCascVarTrueMultEta8",&fTreeCascVarTrueMultEta8,"fTreeCascVarTrueMultEta8/I");
    fTreeCascade->Branch("fTreeCascVarTrueMultVZEROA",&fTreeCascVarTrueMultVZEROA,"fTreeCascVarTrueMultVZEROA/I");
    fTreeCascade->Branch("fTreeCascVarTrueMultVZEROC",&fTreeCascVarTrueMultVZEROC,"fTreeCascVarTrueMultVZEROC/I");
    fTreeCascade->Branch("fTreeCascVarIsPhysicalPrimary",&fTreeCascVarIsPhysicalPrimary,"fTreeCascVarIsPhysicalPrimary/I");
    fTreeCascade->Branch("fTreeCascVarPID",&fTreeCascVarPID,"fTreeCascVarPID/I");
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

    fTreeCascade->Branch("fTreeCascVarNegTransvMomentumMC",&fTreeCascVarNegTransvMomentumMC,"fTreeCascVarNegTransvMomentumMC/F");
    fTreeCascade->Branch("fTreeCascVarPosTransvMomentumMC",&fTreeCascVarPosTransvMomentumMC,"fTreeCascVarPosTransvMomentumMC/F");
//------------------------------------------------
// Particle Identification Setup
//------------------------------------------------

    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();

    // Multiplicity
    if(! fESDtrackCuts ) {
        fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,kFALSE);
        fESDtrackCuts->SetPtRange(0.15);  // adding pt cut
        fESDtrackCuts->SetEtaRange(-1.0, 1.0);
    }
    //Helper
    if(! fPPVsMultUtils ) {
        fPPVsMultUtils = new AliPPVsMultUtils();
    }
    //Analysis Utils
    if(! fUtils ) {
        fUtils = new AliAnalysisUtils();
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
        TString eventClasses;
        if(!fMinbAnalysis) eventClasses= "Processed;IsSelectedTrigger;IsINELgtZERO;IsAcceptedVertexPosition;IsNotPileupSPDInMultBins;HasNoInconsistentSPDandTrackVertices;IsEventSelected";
        else eventClasses="Processed;IncompleteDAQ;BkgRejTrackletsVsCls;SPDpileup;IsSelectedTrigger;RejSPDorTrackVtxNotAvailable;SPDvtxResolution;SPDVtxAndTrackVtxProximity;IsAcceptedPosition";

        TObjArray *arrNames=eventClasses.Tokenize(";");
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",arrNames->GetEntries(),0,arrNames->GetEntries());
        for(int iclass=0; iclass < arrNames->GetEntries(); iclass++) fHistEventCounter->GetXaxis()->SetBinLabel(iclass+1, arrNames->At(iclass)->GetName());
        fListHist->Add(fHistEventCounter);
    }

    //Histograms for Efficiency corrections... a bunch of them
    //1D Histograms - Fine if efficiency doesn't change vs mult (expected)
    //---> Always filled for |y|<0.5
    //V0s: basic Histos
    if(! fHistPt_GenK0Short ) {
        fHistPt_GenK0Short    = new TH1D( "fHistPt_GenK0Short",    "Generated;p_{T} (GeV/c)",200,0,20);
        fListHist->Add(fHistPt_GenK0Short);
    }
    if(! fHistPt_GenLambda ) {
        fHistPt_GenLambda     = new TH1D( "fHistPt_GenLambda",     "Generated;p_{T} (GeV/c)",200,0,20);
        fListHist->Add(fHistPt_GenLambda);
    }
    if(! fHistPt_GenAntiLambda ) {
        fHistPt_GenAntiLambda = new TH1D( "fHistPt_GenAntiLambda", "Generated;p_{T} (GeV/c)",200,0,20);
        fListHist->Add(fHistPt_GenAntiLambda);
    }
    //Cascades: basic Histos
    if(! fHistPt_GenXiMinus ) {
        fHistPt_GenXiMinus    = new TH1D( "fHistPt_GenXiMinus",    "Generated;p_{T} (GeV/c)",200,0,20);
        fListHist->Add(fHistPt_GenXiMinus);
    }
    if(! fHistPt_GenXiPlus ) {
        fHistPt_GenXiPlus     = new TH1D( "fHistPt_GenXiPlus",    "Generated;p_{T} (GeV/c)",200,0,20);
        fListHist->Add(fHistPt_GenXiPlus);
    }
    if(! fHistPt_GenOmegaMinus ) {
        fHistPt_GenOmegaMinus = new TH1D( "fHistPt_GenOmegaMinus", "Generated;p_{T} (GeV/c)",200,0,20);
        fListHist->Add(fHistPt_GenOmegaMinus);
    }
    if(! fHistPt_GenOmegaPlus ) {
        fHistPt_GenOmegaPlus  = new TH1D( "fHistPt_GenOmegaPlus",  "Generated;p_{T} (GeV/c)",200,0,20);
        fListHist->Add(fHistPt_GenOmegaPlus);
    }
    //2D Histos for vs Mult calculation
    if(! fHistPtVsRefMultEta5_GenK0Short ) {
        fHistPtVsRefMultEta5_GenK0Short    = new TH2D( "fHistPtVsRefMultEta5_GenK0Short",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta5_GenK0Short);
    }
    if(! fHistPtVsRefMultEta5_GenLambda ) {
        fHistPtVsRefMultEta5_GenLambda    = new TH2D( "fHistPtVsRefMultEta5_GenLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta5_GenLambda);
    }
    if(! fHistPtVsRefMultEta5_GenAntiLambda ) {
        fHistPtVsRefMultEta5_GenAntiLambda    = new TH2D( "fHistPtVsRefMultEta5_GenAntiLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta5_GenAntiLambda);
    }
    if(! fHistPtVsRefMultEta5_GenXiMinus ) {
        fHistPtVsRefMultEta5_GenXiMinus    = new TH2D( "fHistPtVsRefMultEta5_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta5_GenXiMinus);
    }
    if(! fHistPtVsRefMultEta5_GenXiPlus ) {
        fHistPtVsRefMultEta5_GenXiPlus     = new TH2D( "fHistPtVsRefMultEta5_GenXiPlus",        "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta5_GenXiPlus);
    }
    if(! fHistPtVsRefMultEta5_GenOmegaMinus ) {
        fHistPtVsRefMultEta5_GenOmegaMinus    = new TH2D( "fHistPtVsRefMultEta5_GenOmegaMinus", "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta5_GenOmegaMinus);
    }
    if(! fHistPtVsRefMultEta5_GenOmegaPlus ) {
        fHistPtVsRefMultEta5_GenOmegaPlus     = new TH2D( "fHistPtVsRefMultEta5_GenOmegaPlus",  "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta5_GenOmegaPlus);
    }
    if(! fHistPtVsRefMultEta8_GenK0Short ) {
        fHistPtVsRefMultEta8_GenK0Short    = new TH2D( "fHistPtVsRefMultEta8_GenK0Short",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta8_GenK0Short);
    }
    if(! fHistPtVsRefMultEta8_GenLambda ) {
        fHistPtVsRefMultEta8_GenLambda    = new TH2D( "fHistPtVsRefMultEta8_GenLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta8_GenLambda);
    }
    if(! fHistPtVsRefMultEta8_GenAntiLambda ) {
        fHistPtVsRefMultEta8_GenAntiLambda    = new TH2D( "fHistPtVsRefMultEta8_GenAntiLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta8_GenAntiLambda);
    }
    if(! fHistPtVsRefMultEta8_GenXiMinus ) {
        fHistPtVsRefMultEta8_GenXiMinus    = new TH2D( "fHistPtVsRefMultEta8_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta8_GenXiMinus);
    }
    if(! fHistPtVsRefMultEta8_GenXiPlus ) {
        fHistPtVsRefMultEta8_GenXiPlus     = new TH2D( "fHistPtVsRefMultEta8_GenXiPlus",        "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta8_GenXiPlus);
    }
    if(! fHistPtVsRefMultEta8_GenOmegaMinus ) {
        fHistPtVsRefMultEta8_GenOmegaMinus    = new TH2D( "fHistPtVsRefMultEta8_GenOmegaMinus", "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta8_GenOmegaMinus);
    }
    if(! fHistPtVsRefMultEta8_GenOmegaPlus ) {
        fHistPtVsRefMultEta8_GenOmegaPlus     = new TH2D( "fHistPtVsRefMultEta8_GenOmegaPlus",  "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsRefMultEta8_GenOmegaPlus);
    }

    //Centralities: V0A, V0C, V0M, +Eq
    if(! fHistPtVsCentV0A_GenK0Short ) {
        fHistPtVsCentV0A_GenK0Short    = new TH2D(
            "fHistPtVsCentV0A_GenK0Short",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0A_GenK0Short);
    }
    if(! fHistPtVsCentV0A_GenLambda ) {
        fHistPtVsCentV0A_GenLambda    = new TH2D(
            "fHistPtVsCentV0A_GenLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0A_GenLambda);
    }
    if(! fHistPtVsCentV0A_GenAntiLambda ) {
        fHistPtVsCentV0A_GenAntiLambda    = new TH2D(
            "fHistPtVsCentV0A_GenAntiLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0A_GenAntiLambda);
    }
    if(! fHistPtVsCentV0A_GenXiMinus ) {
        fHistPtVsCentV0A_GenXiMinus    = new TH2D(
            "fHistPtVsCentV0A_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0A_GenXiMinus);
    }
    if(! fHistPtVsCentV0A_GenXiPlus ) {
        fHistPtVsCentV0A_GenXiPlus    = new TH2D(
            "fHistPtVsCentV0A_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0A_GenXiPlus);
    }
    if(! fHistPtVsCentV0A_GenOmegaMinus ) {
        fHistPtVsCentV0A_GenOmegaMinus    = new TH2D(
            "fHistPtVsCentV0A_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0A_GenOmegaMinus);
    }
    if(! fHistPtVsCentV0A_GenOmegaPlus ) {
        fHistPtVsCentV0A_GenOmegaPlus    = new TH2D(
            "fHistPtVsCentV0A_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0A_GenOmegaPlus);
    }

    if(! fHistPtVsCentV0C_GenK0Short ) {
        fHistPtVsCentV0C_GenK0Short    = new TH2D(
            "fHistPtVsCentV0C_GenK0Short",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0C_GenK0Short);
    }
    if(! fHistPtVsCentV0C_GenLambda ) {
        fHistPtVsCentV0C_GenLambda    = new TH2D(
            "fHistPtVsCentV0C_GenLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0C_GenLambda);
    }
    if(! fHistPtVsCentV0C_GenAntiLambda ) {
        fHistPtVsCentV0C_GenAntiLambda    = new TH2D(
            "fHistPtVsCentV0C_GenAntiLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0C_GenAntiLambda);
    }
    if(! fHistPtVsCentV0C_GenXiMinus ) {
        fHistPtVsCentV0C_GenXiMinus    = new TH2D(
            "fHistPtVsCentV0C_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0C_GenXiMinus);
    }
    if(! fHistPtVsCentV0C_GenXiPlus ) {
        fHistPtVsCentV0C_GenXiPlus    = new TH2D(
            "fHistPtVsCentV0C_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0C_GenXiPlus);
    }
    if(! fHistPtVsCentV0C_GenOmegaMinus ) {
        fHistPtVsCentV0C_GenOmegaMinus    = new TH2D(
            "fHistPtVsCentV0C_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0C_GenOmegaMinus);
    }
    if(! fHistPtVsCentV0C_GenOmegaPlus ) {
        fHistPtVsCentV0C_GenOmegaPlus    = new TH2D(
            "fHistPtVsCentV0C_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0C_GenOmegaPlus);
    }

    if(! fHistPtVsCentV0M_GenK0Short ) {
        fHistPtVsCentV0M_GenK0Short    = new TH2D(
            "fHistPtVsCentV0M_GenK0Short",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0M_GenK0Short);
    }
    if(! fHistPtVsCentV0M_GenLambda ) {
        fHistPtVsCentV0M_GenLambda    = new TH2D(
            "fHistPtVsCentV0M_GenLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0M_GenLambda);
    }
    if(! fHistPtVsCentV0M_GenAntiLambda) {
        fHistPtVsCentV0M_GenAntiLambda    = new TH2D(
            "fHistPtVsCentV0M_GenAntiLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0M_GenAntiLambda);
    }
    if(! fHistPtVsCentV0M_GenXiMinus ) {
        fHistPtVsCentV0M_GenXiMinus    = new TH2D(
            "fHistPtVsCentV0M_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0M_GenXiMinus);
    }
    if(! fHistPtVsCentV0M_GenXiPlus ) {
        fHistPtVsCentV0M_GenXiPlus    = new TH2D(
            "fHistPtVsCentV0M_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0M_GenXiPlus);
    }
    if(! fHistPtVsCentV0M_GenOmegaMinus ) {
        fHistPtVsCentV0M_GenOmegaMinus    = new TH2D(
            "fHistPtVsCentV0M_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0M_GenOmegaMinus);
    }
    if(! fHistPtVsCentV0M_GenOmegaPlus ) {
        fHistPtVsCentV0M_GenOmegaPlus    = new TH2D(
            "fHistPtVsCentV0M_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0M_GenOmegaPlus);
    }

    //Special Unselected
    if(! fHistPtVsCentV0MUnselected_GenK0Short ) {
        fHistPtVsCentV0MUnselected_GenK0Short    = new TH2D(
            "fHistPtVsCentV0MUnselected_GenK0Short",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MUnselected_GenK0Short);
    }
    if(! fHistPtVsCentV0MUnselected_GenLambda ) {
        fHistPtVsCentV0MUnselected_GenLambda    = new TH2D(
            "fHistPtVsCentV0MUnselected_GenLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MUnselected_GenLambda);
    }
    if(! fHistPtVsCentV0MUnselected_GenAntiLambda) {
        fHistPtVsCentV0MUnselected_GenAntiLambda    = new TH2D(
            "fHistPtVsCentV0MUnselected_GenAntiLambda",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MUnselected_GenAntiLambda);
    }
    if(! fHistPtVsCentV0MUnselected_GenXiMinus ) {
        fHistPtVsCentV0MUnselected_GenXiMinus    = new TH2D(
            "fHistPtVsCentV0MUnselected_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MUnselected_GenXiMinus);
    }
    if(! fHistPtVsCentV0MUnselected_GenXiPlus ) {
        fHistPtVsCentV0MUnselected_GenXiPlus    = new TH2D(
            "fHistPtVsCentV0MUnselected_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MUnselected_GenXiPlus);
    }
    if(! fHistPtVsCentV0MUnselected_GenOmegaMinus ) {
        fHistPtVsCentV0MUnselected_GenOmegaMinus    = new TH2D(
            "fHistPtVsCentV0MUnselected_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MUnselected_GenOmegaMinus);
    }
    if(! fHistPtVsCentV0MUnselected_GenOmegaPlus ) {
        fHistPtVsCentV0MUnselected_GenOmegaPlus    = new TH2D(
            "fHistPtVsCentV0MUnselected_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MUnselected_GenOmegaPlus);
    }

    //Equalized
    if(! fHistPtVsCentV0AEq_GenXiMinus ) {
        fHistPtVsCentV0AEq_GenXiMinus    = new TH2D(
            "fHistPtVsCentV0AEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0AEq_GenXiMinus);
    }
    if(! fHistPtVsCentV0AEq_GenXiPlus ) {
        fHistPtVsCentV0AEq_GenXiPlus    = new TH2D(
            "fHistPtVsCentV0AEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0AEq_GenXiPlus);
    }
    if(! fHistPtVsCentV0AEq_GenOmegaMinus ) {
        fHistPtVsCentV0AEq_GenOmegaMinus    = new TH2D(
            "fHistPtVsCentV0AEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0AEq_GenOmegaMinus);
    }
    if(! fHistPtVsCentV0AEq_GenOmegaPlus ) {
        fHistPtVsCentV0AEq_GenOmegaPlus    = new TH2D(
            "fHistPtVsCentV0AEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0AEq_GenOmegaPlus);
    }

    if(! fHistPtVsCentV0CEq_GenXiMinus ) {
        fHistPtVsCentV0CEq_GenXiMinus    = new TH2D(
            "fHistPtVsCentV0CEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0CEq_GenXiMinus);
    }
    if(! fHistPtVsCentV0CEq_GenXiPlus ) {
        fHistPtVsCentV0CEq_GenXiPlus    = new TH2D(
            "fHistPtVsCentV0CEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0CEq_GenXiPlus);
    }
    if(! fHistPtVsCentV0CEq_GenOmegaMinus ) {
        fHistPtVsCentV0CEq_GenOmegaMinus    = new TH2D(
            "fHistPtVsCentV0CEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0CEq_GenOmegaMinus);
    }
    if(! fHistPtVsCentV0CEq_GenOmegaPlus ) {
        fHistPtVsCentV0CEq_GenOmegaPlus    = new TH2D(
            "fHistPtVsCentV0CEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0CEq_GenOmegaPlus);
    }

    if(! fHistPtVsCentV0MEq_GenXiMinus ) {
        fHistPtVsCentV0MEq_GenXiMinus    = new TH2D(
            "fHistPtVsCentV0MEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MEq_GenXiMinus);
    }
    if(! fHistPtVsCentV0MEq_GenXiPlus ) {
        fHistPtVsCentV0MEq_GenXiPlus    = new TH2D(
            "fHistPtVsCentV0MEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MEq_GenXiPlus);
    }
    if(! fHistPtVsCentV0MEq_GenOmegaMinus ) {
        fHistPtVsCentV0MEq_GenOmegaMinus    = new TH2D(
            "fHistPtVsCentV0MEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MEq_GenOmegaMinus);
    }
    if(! fHistPtVsCentV0MEq_GenOmegaPlus ) {
        fHistPtVsCentV0MEq_GenOmegaPlus    = new TH2D(
            "fHistPtVsCentV0MEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,200,0,200);
        fListHist->Add(fHistPtVsCentV0MEq_GenOmegaPlus);
    }

    //AMPLITUDES: V0A, V0C, V0M, +Eq
    Double_t lMaxAmplitude = 2500;
    Long_t lAmplitudeBins = 10000;
    if(! fHistPtVsAmpV0A_GenXiMinus ) {
        fHistPtVsAmpV0A_GenXiMinus    = new TH2D(
            "fHistPtVsAmpV0A_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0A_GenXiMinus);
    }
    if(! fHistPtVsAmpV0A_GenXiPlus ) {
        fHistPtVsAmpV0A_GenXiPlus    = new TH2D(
            "fHistPtVsAmpV0A_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0A_GenXiPlus);
    }
    if(! fHistPtVsAmpV0A_GenOmegaMinus ) {
        fHistPtVsAmpV0A_GenOmegaMinus    = new TH2D(
            "fHistPtVsAmpV0A_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0A_GenOmegaMinus);
    }
    if(! fHistPtVsAmpV0A_GenOmegaPlus ) {
        fHistPtVsAmpV0A_GenOmegaPlus    = new TH2D(
            "fHistPtVsAmpV0A_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0A_GenOmegaPlus);
    }

    if(! fHistPtVsAmpV0C_GenXiMinus ) {
        fHistPtVsAmpV0C_GenXiMinus    = new TH2D(
            "fHistPtVsAmpV0C_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0C_GenXiMinus);
    }
    if(! fHistPtVsAmpV0C_GenXiPlus ) {
        fHistPtVsAmpV0C_GenXiPlus    = new TH2D(
            "fHistPtVsAmpV0C_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0C_GenXiPlus);
    }
    if(! fHistPtVsAmpV0C_GenOmegaMinus ) {
        fHistPtVsAmpV0C_GenOmegaMinus    = new TH2D(
            "fHistPtVsAmpV0C_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0C_GenOmegaMinus);
    }
    if(! fHistPtVsAmpV0C_GenOmegaPlus ) {
        fHistPtVsAmpV0C_GenOmegaPlus    = new TH2D(
            "fHistPtVsAmpV0C_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0C_GenOmegaPlus);
    }

    if(! fHistPtVsAmpV0M_GenXiMinus ) {
        fHistPtVsAmpV0M_GenXiMinus    = new TH2D(
            "fHistPtVsAmpV0M_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0M_GenXiMinus);
    }
    if(! fHistPtVsAmpV0M_GenXiPlus ) {
        fHistPtVsAmpV0M_GenXiPlus    = new TH2D(
            "fHistPtVsAmpV0M_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0M_GenXiPlus);
    }
    if(! fHistPtVsAmpV0M_GenOmegaMinus ) {
        fHistPtVsAmpV0M_GenOmegaMinus    = new TH2D(
            "fHistPtVsAmpV0M_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0M_GenOmegaMinus);
    }
    if(! fHistPtVsAmpV0M_GenOmegaPlus ) {
        fHistPtVsAmpV0M_GenOmegaPlus    = new TH2D(
            "fHistPtVsAmpV0M_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0M_GenOmegaPlus);
    }

    //Equalized
    if(! fHistPtVsAmpV0AEq_GenXiMinus ) {
        fHistPtVsAmpV0AEq_GenXiMinus    = new TH2D(
            "fHistPtVsAmpV0AEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0AEq_GenXiMinus);
    }
    if(! fHistPtVsAmpV0AEq_GenXiPlus ) {
        fHistPtVsAmpV0AEq_GenXiPlus    = new TH2D(
            "fHistPtVsAmpV0AEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0AEq_GenXiPlus);
    }
    if(! fHistPtVsAmpV0AEq_GenOmegaMinus ) {
        fHistPtVsAmpV0AEq_GenOmegaMinus    = new TH2D(
            "fHistPtVsAmpV0AEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0AEq_GenOmegaMinus);
    }
    if(! fHistPtVsAmpV0AEq_GenOmegaPlus ) {
        fHistPtVsAmpV0AEq_GenOmegaPlus    = new TH2D(
            "fHistPtVsAmpV0AEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0AEq_GenOmegaPlus);
    }

    if(! fHistPtVsAmpV0CEq_GenXiMinus ) {
        fHistPtVsAmpV0CEq_GenXiMinus    = new TH2D(
            "fHistPtVsAmpV0CEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0CEq_GenXiMinus);
    }
    if(! fHistPtVsAmpV0CEq_GenXiPlus ) {
        fHistPtVsAmpV0CEq_GenXiPlus    = new TH2D(
            "fHistPtVsAmpV0CEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0CEq_GenXiPlus);
    }
    if(! fHistPtVsAmpV0CEq_GenOmegaMinus ) {
        fHistPtVsAmpV0CEq_GenOmegaMinus    = new TH2D(
            "fHistPtVsAmpV0CEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0CEq_GenOmegaMinus);
    }
    if(! fHistPtVsAmpV0CEq_GenOmegaPlus ) {
        fHistPtVsAmpV0CEq_GenOmegaPlus    = new TH2D(
            "fHistPtVsAmpV0CEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0CEq_GenOmegaPlus);
    }

    if(! fHistPtVsAmpV0MEq_GenXiMinus ) {
        fHistPtVsAmpV0MEq_GenXiMinus    = new TH2D(
            "fHistPtVsAmpV0MEq_GenXiMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0MEq_GenXiMinus);
    }
    if(! fHistPtVsAmpV0MEq_GenXiPlus ) {
        fHistPtVsAmpV0MEq_GenXiPlus    = new TH2D(
            "fHistPtVsAmpV0MEq_GenXiPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0MEq_GenXiPlus);
    }
    if(! fHistPtVsAmpV0MEq_GenOmegaMinus ) {
        fHistPtVsAmpV0MEq_GenOmegaMinus    = new TH2D(
            "fHistPtVsAmpV0MEq_GenOmegaMinus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0MEq_GenOmegaMinus);
    }
    if(! fHistPtVsAmpV0MEq_GenOmegaPlus ) {
        fHistPtVsAmpV0MEq_GenOmegaPlus    = new TH2D(
            "fHistPtVsAmpV0MEq_GenOmegaPlus",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistPtVsAmpV0MEq_GenOmegaPlus);
    }

    if(! fHistVZEROResponseStudy ) {
        fHistVZEROResponseStudy    = new TH2D(
            "fHistVZEROResponseStudy",       "Generated;p_{T} (GeV/c); Mult",200,0,20,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistVZEROResponseStudy);
    }

    if(! fHistVZEROResponseStudyTotal ) {
        fHistVZEROResponseStudyTotal    = new TH2D(
            "fHistVZEROResponseStudyTotal",       "Generated;p_{T} (GeV/c); Mult",5000,0,500,lAmplitudeBins,0,lMaxAmplitude);
        fListHist->Add(fHistVZEROResponseStudyTotal);
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
void AliAnalysisTaskStrangenessVsMultiplicityMC::UserExec(Option_t *)
{
    // Main loop
    // Called for each event

    AliESDEvent *lESDevent = 0x0;
    AliMCEvent  *lMCevent  = 0x0;
    AliStack    *lMCstack  = 0x0;

    //Zero all booleans, etc
    fEvSel_HasAtLeastSPDVertex    = kFALSE;
    fEvSel_VtxZCut                = kFALSE;
    fEvSel_IsNotPileup            = kFALSE;
    fEvSel_IsNotPileupInMultBins  = kFALSE;
    fEvSel_HasVtxContributor      = kFALSE;
    fEvSel_Triggered              = kFALSE;
    fEvSel_INELgtZEROStackPrimaries = kFALSE;
    fEvSel_VtxZ = -100;
    fEvSel_MCType = -100;


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

    lMCevent = MCEvent();
    if (!lMCevent) {
        Printf("ERROR: Could not retrieve MC event \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
    }
    //Code for the acquisition of the 'perfect' primary vertex position
    TArrayF mcPrimaryVtx;
    AliGenEventHeader* mcHeader=lMCevent->GenEventHeader();
    if(!mcHeader) return;
    mcHeader->PrimaryVertex(mcPrimaryVtx);
    fEvSel_VtxZMC = mcPrimaryVtx.At(2);

    lMCstack = lMCevent->Stack();
    if (!lMCstack) {
        Printf("ERROR: Could not retrieve MC stack \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
    }

    fRunNumber = lESDevent->GetRunNumber();

    Double_t lMagneticField = -10;
    lMagneticField = lESDevent->GetMagneticField( );

    //Before selections !
    Float_t fCentrality_V0MUnselected = fPPVsMultUtils -> GetMultiplicityPercentile(lESDevent, "V0M" , kFALSE );
    Int_t    lThisPDG  = 0;
    Double_t lThisRap  = 0;
    Double_t lThisPt   = 0;

    //----- Loop on Generated CASCADES ---------------
    for (Int_t ilab = 0;  ilab < (lMCstack->GetNtrack()); ilab++)
    {   // This is the begining of the loop on tracks

        TParticle* lPart = 0x0;
        lPart = lMCstack->Particle( ilab );
        if(!lPart) {
            Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", ilab );
            continue;
        }

        lThisPDG = lPart->GetPdgCode();

        if ( (TMath::Abs(lThisPDG) == 3312) || (TMath::Abs(lThisPDG) == 3334) || (TMath::Abs(lThisPDG) == 3122) || lThisPDG == 310 )
        {
            lThisRap   = MyRapidity(lPart->Energy(),lPart->Pz());
            lThisPt    = lPart->Pt();

            //Use Physical Primaries only for filling These Histos
            if ( lMCstack->IsPhysicalPrimary(ilab)!=kTRUE ) continue;

            if( lThisPDG ==   310 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPtVsCentV0MUnselected_GenK0Short       -> Fill (lThisPt, fCentrality_V0MUnselected);
            }
            if( lThisPDG ==  3122 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPtVsCentV0MUnselected_GenLambda       -> Fill (lThisPt, fCentrality_V0MUnselected);
            }
            if( lThisPDG == -3122 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPtVsCentV0MUnselected_GenAntiLambda       -> Fill (lThisPt, fCentrality_V0MUnselected);
            }
            if( lThisPDG ==  3312 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPtVsCentV0MUnselected_GenXiMinus       -> Fill (lThisPt, fCentrality_V0MUnselected);
            }
            if( lThisPDG == -3312 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPtVsCentV0MUnselected_GenXiPlus       -> Fill (lThisPt, fCentrality_V0MUnselected);
            }
            if( lThisPDG ==  3334 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPtVsCentV0MUnselected_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0MUnselected);
            }
            if( lThisPDG == -3334 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPtVsCentV0MUnselected_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0MUnselected);
            }
        }
    }//End of loop on tracks
//----- End Loop on Cascades ------------------------------------------------------------


//------------------------------------------------
// MC type (ND, SD, DD)
//------------------------------------------------

    AliGenEventHeader * header = lMCevent->GenEventHeader();
    Int_t processtype = AliPWG0Helper::GetPythiaEventProcessType(header);
    // non diffractive
    if (processtype !=92 && processtype !=93 && processtype != 94) fEvSel_MCType = 1;
    // single diffractive
    if ((processtype == 92 || processtype == 93)) fEvSel_MCType = 2;
    // double diffractive
    if (processtype == 94) fEvSel_MCType = 3;

    //  --- Performed entirely via AliPPVsMultUtils 
    // (except removal of incomplete events and SPDClusterVsTracklets cut) 
    //------------------------------------------------

    //Copy-paste of steps done in AliAnalysisTaskSkeleton

    fHistEventCounter->Fill(0.5);
    Bool_t isEventSelected;
    if(!fMinbAnalysis) isEventSelected = SelectEventsMultiplicityDependentAnalysis();
    else isEventSelected = SelectEventsMinimumBiasAnalysis();

    if(!isEventSelected) return;

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
    fEvSel_VtxZ = lBestPrimaryVtxPos[2];

    if(TMath::Abs(lBestPrimaryVtxPos[2]) <= 10.0 ) {
        //Passed selection!
        fEvSel_VtxZCut = kTRUE;
    }

    if( !lESDevent->IsPileupFromSPD()           ) fEvSel_IsNotPileup           = kTRUE;
    if( !lESDevent->IsPileupFromSPDInMultBins() ) fEvSel_IsNotPileupInMultBins = kTRUE;

    //First implementation of pileup from multi-vertexer (simple use of analysis utils)
    //if ( !fUtils->IsPileUpMV( lESDevent ) ) fEvSel_IsNotPileupMV = kTRUE;
    fEvSel_IsNotPileupMV = kFALSE ; //dummy

//------------------------------------------------
// Multiplicity Information Acquistion
//------------------------------------------------

    //Monte Carlo Level information !
    //--------- GENERATED NUMBER OF CHARGED PARTICLES
    // ---> Variable Definition

    Long_t lNchEta5   = 0;
    Long_t lNchEta8   = 0;
    Long_t lNchEta10  = 0;
    Long_t lNchVZEROA = 0;
    Long_t lNchVZEROC = 0;

    Float_t lPtOfParticleInsideVZEROA = -1;
    Float_t lPOfParticleInsideVZEROA = -1;

    //----- Loop on Stack ----------------------------------------------------------------
    for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++)
    {   // This is the begining of the loop on tracks
        TParticle* particleOne = lMCstack->Particle(iCurrentLabelStack);
        if(!particleOne) continue;
        if(!particleOne->GetPDG()) continue;
        Double_t lThisCharge = particleOne->GetPDG()->Charge()/3.;
        if(TMath::Abs(lThisCharge)<0.001) continue;
        if(! (lMCstack->IsPhysicalPrimary(iCurrentLabelStack)) ) continue;

        //Double_t gpt = particleOne -> Pt();
        Double_t geta = particleOne -> Eta();

        if( TMath::Abs(geta) < 0.5 ) lNchEta5++;
        if( TMath::Abs(geta) < 0.8 ) lNchEta8++;
        if( TMath::Abs(geta) < 0.8 ) lNchEta10++;
        if( TMath::Abs(geta) < 1.0 ) fEvSel_INELgtZEROStackPrimaries = kTRUE;
        if( 2.8 < geta && geta < 5.1 ) lNchVZEROA++;
        if( 2.8 < geta && geta < 5.1 ) lPtOfParticleInsideVZEROA = particleOne->Pt();
        if( 2.8 < geta && geta < 5.1 ) lPOfParticleInsideVZEROA = particleOne->P();
        if(-3.7 < geta && geta <-1.7 ) lNchVZEROC++;
    }//End of loop on tracks

    //Attribution
    fTrueMultEta5  = lNchEta5;
    fTrueMultEta8  = lNchEta8;
    fTrueMultEta10 = lNchEta10;
    fTrueMultVZEROA = lNchVZEROA;
    fTrueMultVZEROC = lNchVZEROC;
    //----- End Loop on Stack ------------------------------------------------------------

    //Standard GetReferenceMultiplicity Estimator (0.5 and 0.8)
    fRefMultEta5 = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,0.5);
    fRefMultEta8 = AliPPVsMultUtils::GetStandardReferenceMultiplicity( lESDevent );

    //Differential in eta
    //binning definition
    Float_t lEtaBinning[21] = {-1.5,-1.,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.5};
    Float_t lEtaCent        = -666.;
    Float_t lEtaRange       = -666.;
    for(Int_t i=0; i<20; i++) {
        lEtaCent  = lEtaBinning[i]+(lEtaBinning[i+1]-lEtaBinning[i])/2.;
        lEtaRange = (lEtaBinning[i+1]-lEtaBinning[i])/2.;
        if(i<2 || i>17) fRefMultDiffEta[i] = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTracklets,lEtaRange,lEtaCent);
        else fRefMultDiffEta[i] = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,lEtaRange,lEtaCent);
    }

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
    fAmplitude_V0A = multV0A;
    fAmplitude_V0C = multV0C;
    fAmplitude_V0M = multV0A+multV0C;

    if( fTrueMultVZEROA == 1 ) fHistVZEROResponseStudy->Fill( lPtOfParticleInsideVZEROA , fAmplitude_V0A );
    if( fTrueMultVZEROA == 1 ) fHistVZEROResponseStudyTotal->Fill( lPOfParticleInsideVZEROA , fAmplitude_V0A );

    // Equalized signals // From AliCentralitySelectionTask
    for(Int_t iCh = 32; iCh < 64; ++iCh) {
        Double_t mult = lESDevent->GetVZEROEqMultiplicity(iCh);
        multV0AEq += mult;
    }
    for(Int_t iCh = 0; iCh < 32; ++iCh) {
        Double_t mult = lESDevent->GetVZEROEqMultiplicity(iCh);
        multV0CEq += mult;
    }
    fAmplitude_V0AEq = multV0AEq;
    fAmplitude_V0CEq = multV0CEq;
    fAmplitude_V0MEq = multV0AEq+multV0CEq;

    fCentrality_V0A   = -100;
    fCentrality_V0C   = -100;
    fCentrality_V0M   = -100;
    fCentrality_OnlineV0A   = -100;
    fCentrality_OnlineV0C   = -100;
    fCentrality_OnlineV0M   = -100;
    fCentrality_ADA   = -100;
    fCentrality_ADC   = -100;
    fCentrality_ADM   = -100;
    fCentrality_V0AEq = -100;
    fCentrality_V0CEq = -100;
    fCentrality_V0MEq = -100;
    fCentrality_V0B        = -100;
    fCentrality_V0Apartial = -100;
    fCentrality_V0Cpartial = -100;
    fCentrality_V0S        = -100;
    fCentrality_V0SB       = -100;
    
    if( !fkMultSelection ) {
        fCentrality_V0A   = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0A"   );
        fCentrality_V0C   = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0C"   );
        fCentrality_V0M   = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0M"   );
        fCentrality_V0AEq = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0AEq" );
        fCentrality_V0CEq = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0CEq" );
        fCentrality_V0MEq = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0MEq" );
        fCentrality_V0B        = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0B"        );
        fCentrality_V0Apartial = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0Apartial" );
        fCentrality_V0Cpartial = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0Cpartial" );
        fCentrality_V0S        = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0S"        );
        fCentrality_V0SB       = fPPVsMultUtils->GetMultiplicityPercentile(lESDevent, "V0SB"       );
    }else{
        //Employ MultSelection framework (check if available first!)
        AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
        //Initial test: standard estimators only
        //FIXME: Add ADM/ADA/ADC, OnlineV0x, SPD, etc if deemed necessary
        if( MultSelection ){
            fCentrality_V0M = MultSelection->GetMultiplicityPercentile("V0M");
            fCentrality_V0A = MultSelection->GetMultiplicityPercentile("V0A");
            fCentrality_V0C = MultSelection->GetMultiplicityPercentile("V0C");
            fCentrality_OnlineV0M = MultSelection->GetMultiplicityPercentile("OnlineV0M");
            fCentrality_OnlineV0A = MultSelection->GetMultiplicityPercentile("OnlineV0A");
            fCentrality_OnlineV0C = MultSelection->GetMultiplicityPercentile("OnlineV0C");
            fCentrality_ADM = MultSelection->GetMultiplicityPercentile("ADM");
            fCentrality_ADA = MultSelection->GetMultiplicityPercentile("ADA");
            fCentrality_ADC = MultSelection->GetMultiplicityPercentile("ADC");
        }
    }
    
    //INEL > 0 check
    fEvSel_INELgtZERO          = IsINELgtZERO( lESDevent , "tracks"    );
    fEvSel_INELgtZEROtracklets = IsINELgtZERO( lESDevent , "tracklets" );

    fEvSel_INELgtZERORefMult = kFALSE;
    if ( fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC, 1.0) >= 1 ) fEvSel_INELgtZERORefMult = kTRUE;

    fEvSel_INELgtZERORefMultTracklets = kFALSE;
    fEvSel_nTrackletsEta10 = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTracklets, 1.0);
    if ( fEvSel_nTrackletsEta10 >= 1 ) fEvSel_INELgtZERORefMultTracklets = kTRUE;

    //Event-level fill
    fTreeEvent->Fill();

    //STOP HERE if skipping event selections (no point in doing the rest...)
    if( fkSkipEventSelection ) {
        PostData(1, fListHist);
        PostData(2, fTreeEvent);
        PostData(3, fTreeV0);
        PostData(4, fTreeCascade);
        return;
    }


//------------------------------------------------

//------------------------------------------------
// Fill Efficiency Denominators, please
//------------------------------------------------

    lThisPDG  = 0;
    lThisRap  = 0;
    lThisPt   = 0;

//----- Loop on Generated CASCADES ---------------
    for (Int_t ilab = 0;  ilab < (lMCstack->GetNtrack()); ilab++)
    {   // This is the begining of the loop on tracks

        TParticle* lPart = 0x0;
        lPart = lMCstack->Particle( ilab );
        if(!lPart) {
            Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", ilab );
            continue;
        }

        lThisPDG = lPart->GetPdgCode();

        if ( (TMath::Abs(lThisPDG) == 3312) || (TMath::Abs(lThisPDG) == 3334) || (TMath::Abs(lThisPDG) == 3122) || lThisPDG == 310 )
        {
            lThisRap   = MyRapidity(lPart->Energy(),lPart->Pz());
            lThisPt    = lPart->Pt();

            //Use Physical Primaries only for filling These Histos
            if ( lMCstack->IsPhysicalPrimary(ilab)!=kTRUE ) continue;

            if( lThisPDG ==   310 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPt_GenK0Short       -> Fill ( lThisPt );
                fHistPtVsRefMultEta5_GenK0Short   -> Fill (lThisPt, fRefMultEta5);
                fHistPtVsRefMultEta8_GenK0Short   -> Fill (lThisPt, fRefMultEta8);
                //Centralities
                fHistPtVsCentV0A_GenK0Short       -> Fill (lThisPt, fCentrality_V0A);
                fHistPtVsCentV0C_GenK0Short       -> Fill (lThisPt, fCentrality_V0C);
                fHistPtVsCentV0M_GenK0Short       -> Fill (lThisPt, fCentrality_V0M);
            }
            if( lThisPDG ==  3122 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPt_GenLambda     -> Fill ( lThisPt );
                fHistPtVsRefMultEta5_GenLambda   -> Fill (lThisPt, fRefMultEta5);
                fHistPtVsRefMultEta8_GenLambda   -> Fill (lThisPt, fRefMultEta8);
                //Centralities
                fHistPtVsCentV0A_GenLambda       -> Fill (lThisPt, fCentrality_V0A);
                fHistPtVsCentV0C_GenLambda       -> Fill (lThisPt, fCentrality_V0C);
                fHistPtVsCentV0M_GenLambda       -> Fill (lThisPt, fCentrality_V0M);
            }
            if( lThisPDG == -3122 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPt_GenAntiLambda -> Fill ( lThisPt );
                fHistPtVsRefMultEta5_GenAntiLambda   -> Fill (lThisPt, fRefMultEta5);
                fHistPtVsRefMultEta8_GenAntiLambda   -> Fill (lThisPt, fRefMultEta8);
                //Centralities
                fHistPtVsCentV0A_GenAntiLambda       -> Fill (lThisPt, fCentrality_V0A);
                fHistPtVsCentV0C_GenAntiLambda       -> Fill (lThisPt, fCentrality_V0C);
                fHistPtVsCentV0M_GenAntiLambda       -> Fill (lThisPt, fCentrality_V0M);
            }
            if( lThisPDG ==  3312 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPt_GenXiMinus                -> Fill (lThisPt);
                fHistPtVsRefMultEta5_GenXiMinus   -> Fill (lThisPt, fRefMultEta5);
                fHistPtVsRefMultEta8_GenXiMinus   -> Fill (lThisPt, fRefMultEta8);
                //Centralities
                fHistPtVsCentV0A_GenXiMinus       -> Fill (lThisPt, fCentrality_V0A);
                fHistPtVsCentV0C_GenXiMinus       -> Fill (lThisPt, fCentrality_V0C);
                fHistPtVsCentV0M_GenXiMinus       -> Fill (lThisPt, fCentrality_V0M);
                fHistPtVsCentV0AEq_GenXiMinus       -> Fill (lThisPt, fCentrality_V0AEq);
                fHistPtVsCentV0CEq_GenXiMinus       -> Fill (lThisPt, fCentrality_V0CEq);
                fHistPtVsCentV0MEq_GenXiMinus       -> Fill (lThisPt, fCentrality_V0MEq);
                //Amplitudes
                fHistPtVsAmpV0A_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0A);
                fHistPtVsAmpV0C_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0C);
                fHistPtVsAmpV0M_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0A + fAmplitude_V0C);
                fHistPtVsAmpV0AEq_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0AEq);
                fHistPtVsAmpV0CEq_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0CEq);
                fHistPtVsAmpV0MEq_GenXiMinus       -> Fill (lThisPt, fAmplitude_V0AEq + fAmplitude_V0CEq);
            }
            if( lThisPDG == -3312 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPt_GenXiPlus                -> Fill (lThisPt);
                fHistPtVsRefMultEta5_GenXiPlus   -> Fill (lThisPt, fRefMultEta5);
                fHistPtVsRefMultEta8_GenXiPlus   -> Fill (lThisPt, fRefMultEta8);
                //Centralities
                fHistPtVsCentV0A_GenXiPlus       -> Fill (lThisPt, fCentrality_V0A);
                fHistPtVsCentV0C_GenXiPlus       -> Fill (lThisPt, fCentrality_V0C);
                fHistPtVsCentV0M_GenXiPlus       -> Fill (lThisPt, fCentrality_V0M);
                fHistPtVsCentV0AEq_GenXiPlus       -> Fill (lThisPt, fCentrality_V0AEq);
                fHistPtVsCentV0CEq_GenXiPlus       -> Fill (lThisPt, fCentrality_V0CEq);
                fHistPtVsCentV0MEq_GenXiPlus       -> Fill (lThisPt, fCentrality_V0MEq);
                //Amplitudes
                fHistPtVsAmpV0A_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0A);
                fHistPtVsAmpV0C_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0C);
                fHistPtVsAmpV0M_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0A + fAmplitude_V0C);
                fHistPtVsAmpV0AEq_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0AEq);
                fHistPtVsAmpV0CEq_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0CEq);
                fHistPtVsAmpV0MEq_GenXiPlus       -> Fill (lThisPt, fAmplitude_V0AEq + fAmplitude_V0CEq);
            }
            if( lThisPDG ==  3334 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPt_GenOmegaMinus                -> Fill (lThisPt);
                fHistPtVsRefMultEta5_GenOmegaMinus   -> Fill (lThisPt, fRefMultEta5);
                fHistPtVsRefMultEta8_GenOmegaMinus   -> Fill (lThisPt, fRefMultEta8);
                //Centralities
                fHistPtVsCentV0A_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0A);
                fHistPtVsCentV0C_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0C);
                fHistPtVsCentV0M_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0M);
                fHistPtVsCentV0AEq_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0AEq);
                fHistPtVsCentV0CEq_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0CEq);
                fHistPtVsCentV0MEq_GenOmegaMinus       -> Fill (lThisPt, fCentrality_V0MEq);
                //Amplitudes
                fHistPtVsAmpV0A_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0A);
                fHistPtVsAmpV0C_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0C);
                fHistPtVsAmpV0M_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0A + fAmplitude_V0C);
                fHistPtVsAmpV0AEq_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0AEq);
                fHistPtVsAmpV0CEq_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0CEq);
                fHistPtVsAmpV0MEq_GenOmegaMinus       -> Fill (lThisPt, fAmplitude_V0AEq + fAmplitude_V0CEq);
            }
            if( lThisPDG == -3334 && TMath::Abs(lThisRap) < 0.5 ) {
                fHistPt_GenOmegaPlus                -> Fill (lThisPt);
                fHistPtVsRefMultEta5_GenOmegaPlus   -> Fill (lThisPt, fRefMultEta5);
                fHistPtVsRefMultEta8_GenOmegaPlus   -> Fill (lThisPt, fRefMultEta8);
                //Centralities
                fHistPtVsCentV0A_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0A);
                fHistPtVsCentV0C_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0C);
                fHistPtVsCentV0M_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0M);
                fHistPtVsCentV0AEq_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0AEq);
                fHistPtVsCentV0CEq_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0CEq);
                fHistPtVsCentV0MEq_GenOmegaPlus       -> Fill (lThisPt, fCentrality_V0MEq);
                //Amplitudes
                fHistPtVsAmpV0A_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0A);
                fHistPtVsAmpV0C_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0C);
                fHistPtVsAmpV0M_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0A + fAmplitude_V0C);
                fHistPtVsAmpV0AEq_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0AEq);
                fHistPtVsAmpV0CEq_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0CEq);
                fHistPtVsAmpV0MEq_GenOmegaPlus       -> Fill (lThisPt, fAmplitude_V0AEq + fAmplitude_V0CEq);
            }
        }
    }//End of loop on tracks
//----- End Loop on Cascades ------------------------------------------------------------

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
    {   // This is the begining of the V0 loop
        AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
        if (!v0) continue;

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


//===============================================
// Monte Carlo Association starts here
//===============================================

        //---> Set Everything to "I don't know" before starting

        fTreeVariablePIDPositive = 0;
        fTreeVariablePIDNegative = 0;

        fTreeVariablePtMother = -1;
        fTreeVariableRapMother = -100;
        fTreeVariablePtMC = -1;
        fTreeVariableRapMC = -100;

        fTreeVariablePID = -1;
        fTreeVariablePIDMother = -1;

        fTreeVariablePrimaryStatus = 0;
        fTreeVariablePrimaryStatusMother = 0;

        fTreeVariablePosTransvMomentumMC = -1;
        fTreeVariableNegTransvMomentumMC = -1;

        Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrack->GetLabel() );
        Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrack->GetLabel() );

        TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
        TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );

        Int_t lPIDPositive = mcPosV0Dghter -> GetPdgCode();
        Int_t lPIDNegative = mcNegV0Dghter -> GetPdgCode();

        fTreeVariablePosTransvMomentumMC = mcPosV0Dghter->Pt();
        fTreeVariableNegTransvMomentumMC = mcNegV0Dghter->Pt();

        fTreeVariablePIDPositive = lPIDPositive;
        fTreeVariablePIDNegative = lPIDNegative;

        Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ;
        Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();

        if( lblMotherPosV0Dghter == lblMotherNegV0Dghter && lblMotherPosV0Dghter > -1 ) {
            //either label is fine, they're equal at this stage
            TParticle* pThisV0 = lMCstack->Particle( lblMotherPosV0Dghter );
            //Set tree variables
            fTreeVariablePID   = pThisV0->GetPdgCode(); //PDG Code
            fTreeVariablePtMC  = pThisV0->Pt(); //Perfect Pt

            //Only Interested if it's a Lambda, AntiLambda or K0s
            //Avoid the Junction Bug! PYTHIA has particles with Px=Py=Pz=E=0 occasionally,
            //having particle code 88 (unrecognized by PDG), for documentation purposes.
            //Even ROOT's TParticle::Y() is not prepared to deal with that exception!
            //Note that TParticle::Pt() is immune (that would just return 0)...
            //Though granted that that should be extremely rare in this precise condition...
            if( TMath::Abs(fTreeVariablePID) == 3122 || fTreeVariablePID==310 ) {
                fTreeVariableRapMC = pThisV0->Y(); //Perfect Y
            }
            if( lMCstack->IsPhysicalPrimary       (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 1; //Is Primary!
            if( lMCstack->IsSecondaryFromWeakDecay(lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 2; //Weak Decay!
            if( lMCstack->IsSecondaryFromMaterial (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 3; //Material Int!

            //Now we try to acquire the V0 parent particle, if possible
            Int_t lblThisV0Parent = pThisV0->GetFirstMother();
            if ( lblThisV0Parent > -1 ) { //if it has a parent, get it and store specs
                TParticle* pThisV0Parent = lMCstack->Particle( lblThisV0Parent );
                fTreeVariablePIDMother   = pThisV0Parent->GetPdgCode(); //V0 Mother PDG
                fTreeVariablePtMother    = pThisV0Parent->Pt();         //V0 Mother Pt
                //NOTE: Fill only for charged xi
                if ( TMath::Abs(fTreeVariablePIDMother)==3312) fTreeVariableRapMother   = pThisV0Parent->Y();         //V0 Mother Pt
                //Primary Status for the V0 Mother particle
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

//This requires an Invariant Mass Hypothesis afterwards
        fTreeVariableDistOverTotMom = TMath::Sqrt(
                                          TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
                                          TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
                                          TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
                                      );
        fTreeVariableDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure

        //Copy Multiplicity information
        fTreeVariableCentV0M = fCentrality_V0M;
        fTreeVariableCentV0A = fCentrality_V0A;
        fTreeVariableCentV0C = fCentrality_V0C;
        fTreeVariableCentOnlineV0M = fCentrality_OnlineV0M;
        fTreeVariableCentOnlineV0A = fCentrality_OnlineV0A;
        fTreeVariableCentOnlineV0C = fCentrality_OnlineV0C;
        fTreeVariableCentADM = fCentrality_ADM;
        fTreeVariableCentADA = fCentrality_ADA;
        fTreeVariableCentADC = fCentrality_ADC;
        fTreeVariableCentV0MEq = fCentrality_V0MEq;
        fTreeVariableCentV0AEq = fCentrality_V0AEq;
        fTreeVariableCentV0CEq = fCentrality_V0CEq;
        fTreeVariableRefMultEta8 = fRefMultEta8;
        fTreeVariableRefMultEta5 = fRefMultEta5;
        fTreeVariableRunNumber = fRunNumber;
        for(Int_t i=0; i<20; i++) fTreeVariableRefMultDiffEta[i] = fRefMultDiffEta[i];

//------------------------------------------------
// Fill Tree!
//------------------------------------------------

        // The conditionals are meant to decrease excessive
        // memory usage!

        //First Selection: Reject OnFly
        if( lOnFlyStatus == 0 ) {
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
                (fTreeVariableInvMassLambda    < lUpperLimitLambda  && fTreeVariableInvMassLambda     > lLowerLimitLambda && fkSaveLambda &&
                 (!fkMCAssociation || (fkMCAssociation&&fTreeVariablePID==3122) )
                )
                ||
                //Case 2: AntiLambda Selection
                (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda && fkSaveAntiLambda &&
                 (!fkMCAssociation || (fkMCAssociation&&fTreeVariablePID==-3122) )
                )
                ||
                //Case 3: K0Short Selection
                (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short && fkSaveK0Short &&
                 (!fkMCAssociation || (fkMCAssociation&&fTreeVariablePID==310) )
                ) ) {
                //Pre-selection in case this is AA...
                if ( TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8 && fkSaveV0Tree ) {
                    fTreeV0->Fill();
                }
            }
        }
    }// This is the end of the V0 loop

//------------------------------------------------
// Fill V0 tree over.
//------------------------------------------------

//------------------------------------------------
// Rerun cascade vertexer!
//------------------------------------------------

    if( fkRunVertexers ) {
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

    for (Int_t iXi = 0; iXi < ncascades; iXi++) {
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

        Double_t lRapXi   = -20.0, lRapOmega = -20.0, lRapMC = -20;//  lEta = -20.0, lTheta = 360., lPhi = 720. ;
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
            AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!");
            continue;
        }
        if(lBachIdx == lIdxPosXi) {
            AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!");
            continue;
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

        if ((pStatus&AliESDtrack::kTPCrefit)    == 0) {
            AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!");
            continue;
        }
        if ((nStatus&AliESDtrack::kTPCrefit)    == 0) {
            AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!");
            continue;
        }
        if ((bachStatus&AliESDtrack::kTPCrefit) == 0) {
            AliWarning("Pb / Bach.   track has no TPCrefit ... continue!");
            continue;
        }

        //Extra Debug Information: booleans for ITS refit
        //if ((pStatus&AliESDtrack::kITSrefit)    == 0) { fTreeCascVarkITSRefitPositive = kFALSE; }
        //if ((nStatus&AliESDtrack::kITSrefit)    == 0) { fTreeCascVarkITSRefitNegative = kFALSE; }
        //if ((bachStatus&AliESDtrack::kITSrefit) == 0) { fTreeCascVarkITSRefitBachelor = kFALSE; }

        // 2 - Poor quality related to TPC clusters: lowest cut of 70 clusters
        if(lPosTPCClusters  < 70) {
            AliWarning("Pb / V0 Pos. track has less than 70 TPC clusters ... continue!");
            continue;
        }
        if(lNegTPCClusters  < 70) {
            AliWarning("Pb / V0 Neg. track has less than 70 TPC clusters ... continue!");
            continue;
        }
        if(lBachTPCClusters < 70) {
            AliWarning("Pb / Bach.   track has less than 70 TPC clusters ... continue!");
            continue;
        }
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

        Double_t lXiTransvMomMC= 0. ;
        Int_t lPDGCodeCascade = 0;
        Int_t lPID_BachMother = 0;
        Int_t lPID_NegMother = 0;
        Int_t lPID_PosMother = 0;
        fTreeCascVarIsPhysicalPrimary = 0; // 0: not defined, any candidate may have this
        fTreeCascVarPosTransvMomentumMC = -1;
        fTreeCascVarNegTransvMomentumMC = -1;

        if(fDebug > 5)
            cout 	<< "MC EventNumber : " << lMCevent->Header()->GetEvent()
                    << " / MC event Number in Run : " << lMCevent->Header()->GetEventNrInRun() << endl;

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

        fTreeCascVarPosTransvMomentumMC = mcPosV0Dghter->Pt();
        fTreeCascVarNegTransvMomentumMC = mcNegV0Dghter->Pt();

        //fTreeCascVarPIDPositive = mcPosV0Dghter -> GetPdgCode();
        //fTreeCascVarPIDNegative = mcNegV0Dghter -> GetPdgCode();
        //fTreeCascVarPIDBachelor = mcBach->GetPdgCode();

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

                                    if(lPID_BachMother==lPID_NegMother && lPID_BachMother==lPID_PosMother) {
                                        lPDGCodeCascade = lPID_BachMother;
                                        lXiTransvMomMC = mcMotherBach->Pt();
                                        if( lMCstack->IsPhysicalPrimary       (lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 1; //Is Primary!
                                        if( lMCstack->IsSecondaryFromWeakDecay(lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 2; //Weak Decay!
                                        if( lMCstack->IsSecondaryFromMaterial (lblMotherBach) ) fTreeCascVarIsPhysicalPrimary = 3; //From Material!
                                        if ( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) !=0 ) {
                                            lRapMC = 0.5*TMath::Log( (mcMotherBach->Energy() + mcMotherBach->Pz()) / (mcMotherBach->Energy() - mcMotherBach->Pz() +1.e-13) );
                                        }
                                    }

                                }
                            }
                        }
                    }
                }
            }
        } //Ends all conditionals above...

//----------------------------------------
// Regular MC ASSOCIATION ENDS HERE
//----------------------------------------


        //------------------------------------------------
        // Set Variables for adding to tree
        //------------------------------------------------

        fTreeCascVarCharge	= lChargeXi;
        fTreeCascVarPID = lPDGCodeCascade;
        if(lInvMassXiMinus!=0)    fTreeCascVarMassAsXi = lInvMassXiMinus;
        if(lInvMassXiPlus!=0)     fTreeCascVarMassAsXi = lInvMassXiPlus;
        if(lInvMassOmegaMinus!=0) fTreeCascVarMassAsOmega = lInvMassOmegaMinus;
        if(lInvMassOmegaPlus!=0)  fTreeCascVarMassAsOmega = lInvMassOmegaPlus;
        fTreeCascVarPt = lXiTransvMom;
        fTreeCascVarPtMC = lXiTransvMomMC;
        fTreeCascVarRapXi = lRapXi ;
        fTreeCascVarRapMC = lRapMC ;
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
        fTreeCascVarCentV0M = fCentrality_V0M;
        fTreeCascVarCentV0A = fCentrality_V0A;
        fTreeCascVarCentV0C = fCentrality_V0C;
        fTreeCascVarCentOnlineV0M = fCentrality_OnlineV0M;
        fTreeCascVarCentOnlineV0A = fCentrality_OnlineV0A;
        fTreeCascVarCentOnlineV0C = fCentrality_OnlineV0C;
        fTreeCascVarCentADM = fCentrality_ADM;
        fTreeCascVarCentADA = fCentrality_ADA;
        fTreeCascVarCentADC = fCentrality_ADC;
        fTreeCascVarCentV0MEq = fCentrality_V0MEq;
        fTreeCascVarCentV0AEq = fCentrality_V0AEq;
        fTreeCascVarCentV0CEq = fCentrality_V0CEq;
        fTreeCascVarRefMultEta8 = fRefMultEta8;
        fTreeCascVarRefMultEta5 = fRefMultEta5;
        fTreeCascVarRunNumber = fRunNumber;
        for(Int_t i=0; i<20; i++) fTreeCascVarRefMultDiffEta[i] = fRefMultDiffEta[i];

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
                                   (fTreeCascVarMassAsOmega<1.68+0.075&&fTreeCascVarMassAsOmega>1.68-0.075) ) ) {
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
void AliAnalysisTaskStrangenessVsMultiplicityMC::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query

    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityMC : ouput data container list not available\n");
        return;
    }

    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityMC : fHistEventCounter not available");
        return;
    }

    TCanvas *canCheck = new TCanvas("AliAnalysisTaskStrangenessVsMultiplicityMC","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();

    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskStrangenessVsMultiplicityMC::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}

Bool_t AliAnalysisTaskStrangenessVsMultiplicityMC::IsINELgtZERO(AliESDEvent *lESDevent, TString lType) const
{
    // This function checks if there was at least a tracklet within |eta|<1.0
    // Meant to be a cross-check before a wider implementation of such a check is in place.

    Bool_t lReturnValue = kFALSE; //No track found a priori

    if ( lType == "tracks") {
        //Step 1: Retrieve array of tracks (requires AliESDTrackCuts with some configuration)
        TObjArray* list = fESDtrackCuts->GetAcceptedTracks(lESDevent, kFALSE);
        Int_t nGoodTracks = list->GetEntries();

        // loop over esd tracks
        for (Int_t i=0; i<nGoodTracks; i++)
        {
            //Acquire Track
            AliESDtrack* esdTrack = dynamic_cast<AliESDtrack*> (list->At(i));
            if (!esdTrack)
            {
                AliDebug(AliLog::kError, Form("ERROR: Could not retrieve track %d.", i));
                continue;
            }
            //Get Basic Characteristics of track
            Float_t eta = esdTrack->Eta();
            Float_t pT  = esdTrack->Pt();

            //Check if the track is in desired phase space
            if (TMath::Abs(eta) < 1.0 && pT > 0.15) {
                lReturnValue = kTRUE;
            }
        }
        //don't forget to delete list...
        delete list;
    }

    if (lType == "tracklets") //pure tracklets check only...
    {
        const AliMultiplicity* spdmult = lESDevent->GetMultiplicity();    // spd multiplicity object
        for (Int_t i=0; i<spdmult->GetNumberOfTracklets(); ++i)
        {
            if (TMath::Abs(spdmult->GetEta(i)) < 1.0)
                lReturnValue = kTRUE;
        }
    }

    return lReturnValue;
}


Bool_t AliAnalysisTaskStrangenessVsMultiplicityMC::SelectEventsMultiplicityDependentAnalysis(){
    //
    // event selection for multiplicity dependent analysis
    //   
    AliESDEvent *lESDevent = lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return kFALSE;
    }

    //
    //Test IsEventSelected (embeds all)
    if(!fkSelectTriggerByName && AliPPVsMultUtils::IsEventSelected (lESDevent, fTrigType) ) fHistEventCounter->Fill(6.5);
    else if(fkSelectTriggerByName && AliPPVsMultUtils::IsEventSelected (lESDevent, fTrigName) ) fHistEventCounter->Fill(6.5);

    //------------------------------------------------
    //Step 1: Check for selected Trigger
    //------------------------------------------------
    if(!fkSelectTriggerByName){
       if( !AliPPVsMultUtils::IsSelectedTrigger( lESDevent, fTrigType ) && !fkSkipEventSelection) {
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
        }
    }else{
       if( !AliPPVsMultUtils::IsSelectedTrigger( lESDevent, fTrigName ) && !fkSkipEventSelection) {
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
        }
    }

    //------------------------------------------------
    //Step 1a: Discard incomplete events
    //------------------------------------------------
    if(lESDevent->IsIncompleteDAQ() && !fkSkipEventSelection){
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
    }

    //------------------------------------------------
    //Step 1b: Apply tracklets vs cluster cuts
    //------------------------------------------------
    if(fkApplyTrackletsVsClustersCut && fUtils->IsSPDClusterVsTrackletBG( lESDevent) && !fkSkipEventSelection){
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;     
    } 


    fHistEventCounter->Fill(1.5);

    //------------------------------------------------
    //Step 2: Check for INEL>0
    //------------------------------------------------
    if( !AliPPVsMultUtils::IsINELgtZERO( lESDevent ) && !fkSkipEventSelection) {
        PostData(1, fListHist);
        PostData(2, fTreeEvent);
        PostData(3, fTreeV0);
        PostData(4, fTreeCascade);
        return kFALSE;
    }
    fHistEventCounter->Fill(2.5);

    //------------------------------------------------
    //Step 3: Check for Vertex-Z position
    //------------------------------------------------
    if( !AliPPVsMultUtils::IsAcceptedVertexPosition( lESDevent ) && !fkSkipEventSelection) {
        PostData(1, fListHist);
        PostData(2, fTreeEvent);
        PostData(3, fTreeV0);
        PostData(4, fTreeCascade);
        return kFALSE;
    }
    fHistEventCounter->Fill(3.5);

    //------------------------------------------------
    //Step 4: Check for SPD Pileup
    //------------------------------------------------
    if( !AliPPVsMultUtils::IsNotPileupSPDInMultBins( lESDevent ) && !fkSkipEventSelection) {
        PostData(1, fListHist);
        PostData(2, fTreeEvent);
        PostData(3, fTreeV0);
        PostData(4, fTreeCascade);
        return kFALSE;
    }
    fHistEventCounter->Fill(4.5);

    //------------------------------------------------
    //Step 5: Check for SPD / track vertex consistency
    //------------------------------------------------
    if( !AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices( lESDevent ) && !fkSkipEventSelection) {
        PostData(1, fListHist);
        PostData(2, fTreeEvent);
        PostData(3, fTreeV0);
        PostData(4, fTreeCascade);
        return kFALSE;
    }
    fHistEventCounter->Fill(5.5);

    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
    fEvSel_Triggered = isSelected;
return kTRUE;
}

Bool_t AliAnalysisTaskStrangenessVsMultiplicityMC::SelectEventsMinimumBiasAnalysis(){
    //
    // event selection for minimum bias analysis
    //   
    AliESDEvent *lESDevent = lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return kFALSE;
    }

    // 
    //------------------------------------------------
    //Step 1.1: Discard incomplete events
    //------------------------------------------------
    if(lESDevent->IsIncompleteDAQ() && !fkSkipEventSelection){
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
    }
   
    fHistEventCounter->Fill(1.5);  

    //------------------------------------------------
    //Step 1.2: Apply tracklets vs cluster cuts
    //------------------------------------------------
    if(fkApplyTrackletsVsClustersCut && fUtils->IsSPDClusterVsTrackletBG( lESDevent) && !fkSkipEventSelection){
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
    }

    fHistEventCounter->Fill(2.5);  

    //------------------------------------------------
    //Step 1.3: Check for SPD Pileup
    //------------------------------------------------
    if( lESDevent->IsPileupFromSPD(3) && !fkSkipEventSelection) {
        PostData(1, fListHist);
        PostData(2, fTreeEvent);
        PostData(3, fTreeV0);
        PostData(4, fTreeCascade);
        return kFALSE;
    }
    fHistEventCounter->Fill(3.5);

    //------------------------------------------------
    //Step 2: Check for selected Trigger
    //------------------------------------------------
    if(!fkSelectTriggerByName){
       if( !AliPPVsMultUtils::IsSelectedTrigger( lESDevent, fTrigType ) && !fkSkipEventSelection) {
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
        }
    }else{
       if( !AliPPVsMultUtils::IsSelectedTrigger( lESDevent, fTrigName ) && !fkSkipEventSelection) {
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
        }
    }
    fHistEventCounter->Fill(4.5);
    
    //------------------------------------------------
    //Step 3:  Primary Vertex quality selection
    //------------------------------------------------
    const AliESDVertex *lESDPrimaryTrackingVtx = lESDevent->GetPrimaryVertexTracks();
    const AliESDVertex *lESDPrimarySPDVtx      = lESDevent->GetPrimaryVertexSPD();
    const AliESDVertex *lPrimaryBestESDVtx     = lESDevent->GetPrimaryVertex();

    //------------------------------------------------
    //Step 3.1: reject events if SPDVtx or TrackVtx is not available
    //------------------------------------------------
    if(!(lESDPrimarySPDVtx->GetStatus() && lESDPrimaryTrackingVtx->GetStatus()) && !fkSkipEventSelection){
	   PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
          }
    fHistEventCounter->Fill(5.5);

    //------------------------------------------------
    //Step 3.2: check the spd vertex resolution and reject if not satisfied
    //------------------------------------------------
    if (lESDPrimarySPDVtx->GetStatus() && lESDPrimarySPDVtx->IsFromVertexerZ() && !(lESDPrimarySPDVtx->GetDispersion()<0.04 && lESDPrimarySPDVtx->GetZRes()<0.25) && !fkSkipEventSelection){
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
          } 
    fHistEventCounter->Fill(6.5);

    //------------------------------------------------
    //Step 3.3: check the proximity between the spd vertex and trak vertex, and reject if not satisfied 
    //------------------------------------------------
    if((TMath::Abs(lESDPrimarySPDVtx->GetZ() - lESDPrimaryTrackingVtx->GetZ()) > 0.5) && !fkSkipEventSelection){
           PostData(1, fListHist);
           PostData(2, fTreeEvent);
           PostData(3, fTreeV0);
           PostData(4, fTreeCascade);
           return kFALSE;
          }
    fHistEventCounter->Fill(7.5);

    //------------------------------------------------
    //Step 3.4: Check for Vertex-Z position
    //------------------------------------------------
    //classical Proton-proton like selection
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
 
    if( TMath::Abs(lBestPrimaryVtxPos[2])>10.0 && !fkSkipEventSelection) {
        PostData(1, fListHist);
        PostData(2, fTreeEvent);
        PostData(3, fTreeV0);
        PostData(4, fTreeCascade);
        return kFALSE;
       }
    fHistEventCounter->Fill(8.5);
return kTRUE;
}

