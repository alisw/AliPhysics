/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//	        AliAnalysisTaskPerformanceSrange class
//    This task is for a performance study of V0 identification.
//                It works with MC info and ESD tree.
//                 Author: H.Ricaud, H.Ricaud@gsi.de
//-----------------------------------------------------------------

#include <Riostream.h>

#include <stdio.h>
#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"

#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"

#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODMCHeader.h"
#include "AliAODInputHandler.h"
//#include "AliV0vertexer.h"

#include "AliAODMCParticle.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliLog.h"

#include "AliAnalysisTaskPerformanceStrange.h"


ClassImp(AliAnalysisTaskPerformanceStrange)


//________________________________________________________________________
AliAnalysisTaskPerformanceStrange::AliAnalysisTaskPerformanceStrange()
  : AliAnalysisTaskSE(), fAnalysisType("ESD"),  fCollidingSystems(0), fUsePID("infoPID"), fUseCut("infoCut"), fListHist(0), 
    fHistMCMultiplicityPrimary(0),
    fHistMCMultiplicityTracks(0),
    fHistMCtracksProdRadiusK0s(0),
    fHistMCtracksProdRadiusLambda(0),
    fHistMCtracksProdRadiusAntiLambda(0),
    fHistMCtracksDecayRadiusK0s(0),
    fHistMCtracksDecayRadiusLambda(0),
    fHistMCtracksDecayRadiusAntiLambda(0),
    fHistMCPtAllK0s(0),
    fHistMCPtAllLambda(0),
    fHistMCPtAllAntiLambda(0),
    fHistMCProdRadiusK0s(0),
    fHistMCProdRadiusLambda(0),
    fHistMCProdRadiusAntiLambda(0),
    fHistMCPtVsYK0s(0),
    fHistMCPtVsYLambda(0),
    fHistMCPtVsYAntiLambda(0),
    fHistMCPtLambdaFromSigma(0),
    fHistMCPtAntiLambdaFromSigma(0),
    fHistNTimesRecK0s(0),
    fHistNTimesRecK0sMI(0),
    fHistNTimesRecLambda(0),
    fHistNTimesRecLambdaMI(0),
    fHistNTimesRecAntiLambda(0),
    fHistNTimesRecAntiLambdaMI(0),
    fHistNTimesRecK0sVsPt(0),
    fHistNTimesRecK0sVsPtMI(0),
    fHistNTimesRecLambdaVsPt(0),
    fHistNTimesRecLambdaVsPtMI(0),
    fHistNTimesRecAntiLambdaVsPt(0),
    fHistNTimesRecAntiLambdaVsPtMI(0),
    fHistTrackPerEvent(0),
    fHistMCDaughterTrack(0),
    fHistPrimaryVertexX(0),
    fHistPrimaryVertexY(0),
    fHistPrimaryVertexZ(0),
    fHistDcaPosToPrimVertex(0),
    fHistDcaNegToPrimVertex(0),
    fHistDcaPosToPrimVertexZoom(0),
    fHistDcaNegToPrimVertexZoom(0),
    fHistRadiusV0(0),
    fHistDecayLengthV0(0),
    fHistDcaV0Daughters(0),
    fHistChi2(0),
    fHistCosPointAngle(0),
    fHistCosPointAngleZoom(0),
    fHistProdRadius(0),
    fHistProdRadiusMI(0),
    fHistV0Multiplicity(0),
    fHistV0MultiplicityMI(0),
    fHistPtVsYK0s(0),
    fHistPtVsYK0sMI(0),
    fHistPtVsYLambda(0),
    fHistPtVsYLambdaMI(0),
    fHistPtVsYAntiLambda(0),
    fHistPtVsYAntiLambdaMI(0),
    fHistMassK0(0),
    fHistMassK0MI(0),
    fHistMassLambda(0),
    fHistMassLambdaMI(0),
    fHistMassAntiLambda(0),
    fHistMassAntiLambdaMI(0),
    fHistMassVsRadiusK0(0),
    fHistMassVsRadiusK0MI(0),
    fHistMassVsRadiusLambda(0),
    fHistMassVsRadiusLambdaMI(0),
    fHistMassVsRadiusAntiLambda(0),
    fHistMassVsRadiusAntiLambdaMI(0),
    fHistPtVsMassK0(0),
    fHistPtVsMassK0MI(0),
    fHistPtVsMassLambda(0),
    fHistPtVsMassLambdaMI(0),
    fHistPtVsMassAntiLambda(0),
    fHistPtVsMassAntiLambdaMI(0),
    fHistArmenterosPodolanski(0),
    fHistArmenterosPodolanskiMI(0),
    fHistNsigmaPosPionAntiLambda(0),
    fHistNsigmaNegProtonAntiLambda(0),
    fHistNsigmaPosProtonLambda(0),
    fHistNsigmaNegPionLambda(0),
    fHistNsigmaPosPionK0(0),
    fHistNsigmaNegPionK0(0),
    fHistAsMcPtK0(0),
    fHistAsMcPtK0MI(0),
    fHistAsMcPtLambda(0),
    fHistAsMcPtLambdaMI(0),
    fHistAsMcPtAntiLambda(0),
    fHistAsMcPtAntiLambdaMI(0),
    fHistAsMcPtZoomK0(0),
    fHistAsMcPtZoomK0MI(0),
    fHistAsMcPtZoomLambda(0),
    fHistAsMcPtZoomLambdaMI(0),
    fHistAsMcProdRadiusK0(0),
    fHistAsMcProdRadiusK0MI(0),
    fHistAsMcProdRadiusLambda(0),
    fHistAsMcProdRadiusLambdaMI(0),
    fHistAsMcProdRadiusAntiLambda(0),
    fHistAsMcProdRadiusAntiLambdaMI(0),
    fHistAsMcProdRadiusXvsYK0s(0),
    fHistAsMcProdRadiusXvsYK0sMI(0),
    fHistAsMcProdRadiusXvsYLambda(0),
    fHistAsMcProdRadiusXvsYLambdaMI(0),
    fHistAsMcProdRadiusXvsYAntiLambda(0),
    fHistAsMcProdRadiusXvsYAntiLambdaMI(0),
    fHistPidMcMassK0(0),
    fHistPidMcMassK0MI(0),
    fHistPidMcMassLambda(0),
    fHistPidMcMassLambdaMI(0),
    fHistPidMcMassAntiLambda(0),
    fHistPidMcMassAntiLambdaMI(0),
    fHistAsMcMassK0(0),
    fHistAsMcMassK0MI(0),
    fHistAsMcMassLambda(0),
    fHistAsMcMassLambdaMI(0),
    fHistAsMcMassAntiLambda(0),
    fHistAsMcMassAntiLambdaMI(0),
    fHistAsMcPtVsMassK0(0),
    fHistAsMcPtVsMassK0MI(0),
    fHistAsMcPtVsMassLambda(0),
    fHistAsMcPtVsMassLambdaMI(0),
    fHistAsMcPtVsMassAntiLambda(0),
    fHistAsMcPtVsMassAntiLambdaMI(0),
    fHistAsMcMassVsRadiusK0(0),
    fHistAsMcMassVsRadiusK0MI(0),
    fHistAsMcMassVsRadiusLambda(0),
    fHistAsMcMassVsRadiusLambdaMI(0),
    fHistAsMcMassVsRadiusAntiLambda(0),
    fHistAsMcMassVsRadiusAntiLambdaMI(0),
    fHistAsMcResxK0(0),
    fHistAsMcResyK0(0),
    fHistAsMcReszK0(0),
    fHistAsMcResrVsRadiusK0(0),
    fHistAsMcReszVsRadiusK0(0),
    fHistAsMcResxK0MI(0),
    fHistAsMcResyK0MI(0),
    fHistAsMcReszK0MI(0),
    fHistAsMcResrVsRadiusK0MI(0),
    fHistAsMcReszVsRadiusK0MI(0),
    fHistAsMcResxLambda(0),
    fHistAsMcResyLambda(0),
    fHistAsMcReszLambda(0),
    fHistAsMcResrVsRadiusLambda(0),
    fHistAsMcReszVsRadiusLambda(0),
    fHistAsMcResxLambdaMI(0),
    fHistAsMcResyLambdaMI(0),
    fHistAsMcReszLambdaMI(0),
    fHistAsMcResrVsRadiusLambdaMI(0),
    fHistAsMcReszVsRadiusLambdaMI(0),
    fHistAsMcResxAntiLambda(0),
    fHistAsMcResyAntiLambda(0),
    fHistAsMcReszAntiLambda(0),
    fHistAsMcResrVsRadiusAntiLambda(0),
    fHistAsMcReszVsRadiusAntiLambda(0),
    fHistAsMcResxAntiLambdaMI(0),
    fHistAsMcResyAntiLambdaMI(0),
    fHistAsMcReszAntiLambdaMI(0),
    fHistAsMcResrVsRadiusAntiLambdaMI(0),
    fHistAsMcReszVsRadiusAntiLambdaMI(0),
    fHistAsMcResPtK0(0),
    fHistAsMcResPtK0MI(0),
    fHistAsMcResPtLambda(0),
    fHistAsMcResPtLambdaMI(0),
    fHistAsMcResPtAntiLambda(0),
    fHistAsMcResPtAntiLambdaMI(0),
    fHistAsMcResPtVsRapK0(0),
    fHistAsMcResPtVsRapK0MI(0),
    fHistAsMcResPtVsRapLambda(0),
    fHistAsMcResPtVsRapLambdaMI(0),
    fHistAsMcResPtVsRapAntiLambda(0),
    fHistAsMcResPtVsRapAntiLambdaMI(0),
    fHistAsMcResPtVsPtK0(0),
    fHistAsMcResPtVsPtK0MI(0),
    fHistAsMcResPtVsPtLambda(0),
    fHistAsMcResPtVsPtLambdaMI(0),
    fHistAsMcResPtVsPtAntiLambda(0),
    fHistAsMcResPtVsPtAntiLambdaMI(0),
    fHistAsMcMotherPdgCodeK0s(0),
    fHistAsMcMotherPdgCodeK0sMI(0),
    fHistAsMcMotherPdgCodeLambda(0),
    fHistAsMcMotherPdgCodeLambdaMI(0),
    fHistAsMcMotherPdgCodeAntiLambda(0),
    fHistAsMcMotherPdgCodeAntiLambdaMI(0),
    fHistAsMcPtLambdaFromSigma(0),
    fHistAsMcPtLambdaFromSigmaMI(0),
    fHistAsMcPtAntiLambdaFromSigma(0),
    fHistAsMcPtAntiLambdaFromSigmaMI(0),
    fHistAsMcSecondaryPtVsYK0s(0),
    fHistAsMcSecondaryPtVsYK0sMI(0),
    fHistAsMcSecondaryPtVsYLambda(0),
    fHistAsMcSecondaryPtVsYLambdaMI(0),
    fHistAsMcSecondaryPtVsYAntiLambda(0),
    fHistAsMcSecondaryPtVsYAntiLambdaMI(0),
    fHistAsMcSecondaryProdRadiusK0s(0),
    fHistAsMcSecondaryProdRadiusK0sMI(0),
    fHistAsMcSecondaryProdRadiusLambda(0),
    fHistAsMcSecondaryProdRadiusLambdaMI(0),
    fHistAsMcSecondaryProdRadiusAntiLambda(0),
    fHistAsMcSecondaryProdRadiusAntiLambdaMI(0),
    fHistAsMcSecondaryProdRadiusXvsYK0s(0),
    fHistAsMcSecondaryProdRadiusXvsYK0sMI(0),
    fHistAsMcSecondaryProdRadiusXvsYLambda(0),
    fHistAsMcSecondaryProdRadiusXvsYLambdaMI(0),
    fHistAsMcSecondaryProdRadiusXvsYAntiLambda(0),
    fHistAsMcSecondaryProdRadiusXvsYAntiLambdaMI(0),
    fHistAsMcSecondaryMotherPdgCodeK0s(0),
    fHistAsMcSecondaryMotherPdgCodeK0sMI(0),
    fHistAsMcSecondaryMotherPdgCodeLambda(0),
    fHistAsMcSecondaryMotherPdgCodeLambdaMI(0),
    fHistAsMcSecondaryMotherPdgCodeAntiLambda(0),
    fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI(0),
    fHistAsMcSecondaryPtLambdaFromSigma(0),
    fHistAsMcSecondaryPtLambdaFromSigmaMI(0),
    fHistAsMcSecondaryPtAntiLambdaFromSigma(0),
    fHistAsMcSecondaryPtAntiLambdaFromSigmaMI(0)
    
{
  // dummy Constructor
}





//________________________________________________________________________
AliAnalysisTaskPerformanceStrange::AliAnalysisTaskPerformanceStrange(const char *name)
  : AliAnalysisTaskSE(name), fAnalysisType("ESD"), fCollidingSystems(0), fUsePID("infoPID"), fUseCut("infocut"), fListHist(),
    fHistMCMultiplicityPrimary(0),
    fHistMCMultiplicityTracks(0),
    fHistMCtracksProdRadiusK0s(0),
    fHistMCtracksProdRadiusLambda(0),
    fHistMCtracksProdRadiusAntiLambda(0),
    fHistMCtracksDecayRadiusK0s(0),
    fHistMCtracksDecayRadiusLambda(0),
    fHistMCtracksDecayRadiusAntiLambda(0),
    fHistMCPtAllK0s(0),
    fHistMCPtAllLambda(0),
    fHistMCPtAllAntiLambda(0),
    fHistMCProdRadiusK0s(0),
    fHistMCProdRadiusLambda(0),
    fHistMCProdRadiusAntiLambda(0),
    fHistMCPtVsYK0s(0),
    fHistMCPtVsYLambda(0),
    fHistMCPtVsYAntiLambda(0),
    fHistMCPtLambdaFromSigma(0),
    fHistMCPtAntiLambdaFromSigma(0),
    fHistNTimesRecK0s(0),
    fHistNTimesRecK0sMI(0),
    fHistNTimesRecLambda(0),
    fHistNTimesRecLambdaMI(0),
    fHistNTimesRecAntiLambda(0),
    fHistNTimesRecAntiLambdaMI(0),
    fHistNTimesRecK0sVsPt(0),
    fHistNTimesRecK0sVsPtMI(0),
    fHistNTimesRecLambdaVsPt(0),
    fHistNTimesRecLambdaVsPtMI(0),
    fHistNTimesRecAntiLambdaVsPt(0),
    fHistNTimesRecAntiLambdaVsPtMI(0),
    fHistTrackPerEvent(0),
    fHistMCDaughterTrack(0),
    fHistPrimaryVertexX(0),
    fHistPrimaryVertexY(0),
    fHistPrimaryVertexZ(0),
    fHistDcaPosToPrimVertex(0),
    fHistDcaNegToPrimVertex(0),
    fHistDcaPosToPrimVertexZoom(0),
    fHistDcaNegToPrimVertexZoom(0),
    fHistRadiusV0(0),
    fHistDecayLengthV0(0),
    fHistDcaV0Daughters(0),
    fHistChi2(0),
    fHistCosPointAngle(0),
    fHistCosPointAngleZoom(0),
    fHistProdRadius(0),
    fHistProdRadiusMI(0),
    fHistV0Multiplicity(0),
    fHistV0MultiplicityMI(0),
    fHistPtVsYK0s(0),
    fHistPtVsYK0sMI(0),
    fHistPtVsYLambda(0),
    fHistPtVsYLambdaMI(0),
    fHistPtVsYAntiLambda(0),
    fHistPtVsYAntiLambdaMI(0),
    fHistMassK0(0),
    fHistMassK0MI(0),
    fHistMassLambda(0),
    fHistMassLambdaMI(0),
    fHistMassAntiLambda(0),
    fHistMassAntiLambdaMI(0),
    fHistMassVsRadiusK0(0),
    fHistMassVsRadiusK0MI(0),
    fHistMassVsRadiusLambda(0),
    fHistMassVsRadiusLambdaMI(0),
    fHistMassVsRadiusAntiLambda(0),
    fHistMassVsRadiusAntiLambdaMI(0),
    fHistPtVsMassK0(0),
    fHistPtVsMassK0MI(0),
    fHistPtVsMassLambda(0),
    fHistPtVsMassLambdaMI(0),
    fHistPtVsMassAntiLambda(0),
    fHistPtVsMassAntiLambdaMI(0),
    fHistArmenterosPodolanski(0),
    fHistArmenterosPodolanskiMI(0),
    fHistNsigmaPosPionAntiLambda(0),
    fHistNsigmaNegProtonAntiLambda(0),
    fHistNsigmaPosProtonLambda(0),
    fHistNsigmaNegPionLambda(0),
    fHistNsigmaPosPionK0(0),
    fHistNsigmaNegPionK0(0),
    fHistAsMcPtK0(0),
    fHistAsMcPtK0MI(0),
    fHistAsMcPtLambda(0),
    fHistAsMcPtLambdaMI(0),
    fHistAsMcPtAntiLambda(0),
    fHistAsMcPtAntiLambdaMI(0),
    fHistAsMcPtZoomK0(0),
    fHistAsMcPtZoomK0MI(0),
    fHistAsMcPtZoomLambda(0),
    fHistAsMcPtZoomLambdaMI(0),
    fHistAsMcProdRadiusK0(0),
    fHistAsMcProdRadiusK0MI(0),
    fHistAsMcProdRadiusLambda(0),
    fHistAsMcProdRadiusLambdaMI(0),
    fHistAsMcProdRadiusAntiLambda(0),
    fHistAsMcProdRadiusAntiLambdaMI(0),
    fHistAsMcProdRadiusXvsYK0s(0),
    fHistAsMcProdRadiusXvsYK0sMI(0),
    fHistAsMcProdRadiusXvsYLambda(0),
    fHistAsMcProdRadiusXvsYLambdaMI(0),
    fHistAsMcProdRadiusXvsYAntiLambda(0),
    fHistAsMcProdRadiusXvsYAntiLambdaMI(0),
    fHistPidMcMassK0(0),
    fHistPidMcMassK0MI(0),
    fHistPidMcMassLambda(0),
    fHistPidMcMassLambdaMI(0),
    fHistPidMcMassAntiLambda(0),
    fHistPidMcMassAntiLambdaMI(0),
    fHistAsMcMassK0(0),
    fHistAsMcMassK0MI(0),
    fHistAsMcMassLambda(0),
    fHistAsMcMassLambdaMI(0),
    fHistAsMcMassAntiLambda(0),
    fHistAsMcMassAntiLambdaMI(0),
    fHistAsMcPtVsMassK0(0),
    fHistAsMcPtVsMassK0MI(0),
    fHistAsMcPtVsMassLambda(0),
    fHistAsMcPtVsMassLambdaMI(0),
    fHistAsMcPtVsMassAntiLambda(0),
    fHistAsMcPtVsMassAntiLambdaMI(0),
    fHistAsMcMassVsRadiusK0(0),
    fHistAsMcMassVsRadiusK0MI(0),
    fHistAsMcMassVsRadiusLambda(0),
    fHistAsMcMassVsRadiusLambdaMI(0),
    fHistAsMcMassVsRadiusAntiLambda(0),
    fHistAsMcMassVsRadiusAntiLambdaMI(0),
    fHistAsMcResxK0(0),
    fHistAsMcResyK0(0),
    fHistAsMcReszK0(0),
    fHistAsMcResrVsRadiusK0(0),
    fHistAsMcReszVsRadiusK0(0),
    fHistAsMcResxK0MI(0),
    fHistAsMcResyK0MI(0),
    fHistAsMcReszK0MI(0),
    fHistAsMcResrVsRadiusK0MI(0),
    fHistAsMcReszVsRadiusK0MI(0),
    fHistAsMcResxLambda(0),
    fHistAsMcResyLambda(0),
    fHistAsMcReszLambda(0),
    fHistAsMcResrVsRadiusLambda(0),
    fHistAsMcReszVsRadiusLambda(0),
    fHistAsMcResxLambdaMI(0),
    fHistAsMcResyLambdaMI(0),
    fHistAsMcReszLambdaMI(0),
    fHistAsMcResrVsRadiusLambdaMI(0),
    fHistAsMcReszVsRadiusLambdaMI(0),
    fHistAsMcResxAntiLambda(0),
    fHistAsMcResyAntiLambda(0),
    fHistAsMcReszAntiLambda(0),
    fHistAsMcResrVsRadiusAntiLambda(0),
    fHistAsMcReszVsRadiusAntiLambda(0),
    fHistAsMcResxAntiLambdaMI(0),
    fHistAsMcResyAntiLambdaMI(0),
    fHistAsMcReszAntiLambdaMI(0),
    fHistAsMcResrVsRadiusAntiLambdaMI(0),
    fHistAsMcReszVsRadiusAntiLambdaMI(0),
    fHistAsMcResPtK0(0),
    fHistAsMcResPtK0MI(0),
    fHistAsMcResPtLambda(0),
    fHistAsMcResPtLambdaMI(0),
    fHistAsMcResPtAntiLambda(0),
    fHistAsMcResPtAntiLambdaMI(0),
    fHistAsMcResPtVsRapK0(0),
    fHistAsMcResPtVsRapK0MI(0),
    fHistAsMcResPtVsRapLambda(0),
    fHistAsMcResPtVsRapLambdaMI(0),
    fHistAsMcResPtVsRapAntiLambda(0),
    fHistAsMcResPtVsRapAntiLambdaMI(0),
    fHistAsMcResPtVsPtK0(0),
    fHistAsMcResPtVsPtK0MI(0),
    fHistAsMcResPtVsPtLambda(0),
    fHistAsMcResPtVsPtLambdaMI(0),
    fHistAsMcResPtVsPtAntiLambda(0),
    fHistAsMcResPtVsPtAntiLambdaMI(0),
    fHistAsMcMotherPdgCodeK0s(0),
    fHistAsMcMotherPdgCodeK0sMI(0),
    fHistAsMcMotherPdgCodeLambda(0),
    fHistAsMcMotherPdgCodeLambdaMI(0),
    fHistAsMcMotherPdgCodeAntiLambda(0),
    fHistAsMcMotherPdgCodeAntiLambdaMI(0),
    fHistAsMcPtLambdaFromSigma(0),
    fHistAsMcPtLambdaFromSigmaMI(0),
    fHistAsMcPtAntiLambdaFromSigma(0),
    fHistAsMcPtAntiLambdaFromSigmaMI(0),
    fHistAsMcSecondaryPtVsYK0s(0),
    fHistAsMcSecondaryPtVsYK0sMI(0),
    fHistAsMcSecondaryPtVsYLambda(0),
    fHistAsMcSecondaryPtVsYLambdaMI(0),
    fHistAsMcSecondaryPtVsYAntiLambda(0),
    fHistAsMcSecondaryPtVsYAntiLambdaMI(0),
    fHistAsMcSecondaryProdRadiusK0s(0),
    fHistAsMcSecondaryProdRadiusK0sMI(0),
    fHistAsMcSecondaryProdRadiusLambda(0),
    fHistAsMcSecondaryProdRadiusLambdaMI(0),
    fHistAsMcSecondaryProdRadiusAntiLambda(0),
    fHistAsMcSecondaryProdRadiusAntiLambdaMI(0),
    fHistAsMcSecondaryProdRadiusXvsYK0s(0),
    fHistAsMcSecondaryProdRadiusXvsYK0sMI(0),
    fHistAsMcSecondaryProdRadiusXvsYLambda(0),
    fHistAsMcSecondaryProdRadiusXvsYLambdaMI(0),
    fHistAsMcSecondaryProdRadiusXvsYAntiLambda(0),
    fHistAsMcSecondaryProdRadiusXvsYAntiLambdaMI(0),
    fHistAsMcSecondaryMotherPdgCodeK0s(0),
    fHistAsMcSecondaryMotherPdgCodeK0sMI(0),
    fHistAsMcSecondaryMotherPdgCodeLambda(0),
    fHistAsMcSecondaryMotherPdgCodeLambdaMI(0),
    fHistAsMcSecondaryMotherPdgCodeAntiLambda(0),
    fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI(0),
    fHistAsMcSecondaryPtLambdaFromSigma(0),
    fHistAsMcSecondaryPtLambdaFromSigmaMI(0),
    fHistAsMcSecondaryPtAntiLambdaFromSigma(0),
    fHistAsMcSecondaryPtAntiLambdaFromSigmaMI(0)
    
{
  // Constructor
  // Define output slots only here
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once

  fListHist = new TList();
  //AliLog::SetGlobalLogLevel(AliLog::kError);

  // Bo: tbd: condition before allocation (i.e. if (!fHistMCMultiplicityPrimary){...} for each histo...

  //***************
  // MC histograms
  //***************
  
  // Multiplicity
  fHistMCMultiplicityPrimary           = new TH1F("h1MCMultiplicityPrimary", "MC Primary Particles;NPrimary;Count", 201, -0.5, 200.5);
  fListHist->Add(fHistMCMultiplicityPrimary);

  fHistMCMultiplicityTracks            = new TH1F("h1MCMultiplicityTracks", "MC Tracks;Ntracks;Count", 201, -0.5, 200.5);
  fListHist->Add(fHistMCMultiplicityTracks);

  // Production Radius of non-primary particles:
  fHistMCtracksProdRadiusK0s           = new TH2F("h2MCtracksProdRadiusK0s","Non-primary MC K^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistMCtracksProdRadiusK0s);

  fHistMCtracksProdRadiusLambda        = new TH2F("h2MCtracksProdRadiusLambda","Non-primary MC #Lambda^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistMCtracksProdRadiusLambda);

  fHistMCtracksProdRadiusAntiLambda    = new TH2F("h2MCtracksProdRadiusAntiLambda","Non-primary MC #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistMCtracksProdRadiusAntiLambda);

  // Decay Radius of non-primary particles:
  fHistMCtracksDecayRadiusK0s          = new TH1F("h1MCtracksDecayRadiusK0s","Non-primary MC K^{0} Decay Radius;r (cm)",101,-1,100);
  fListHist->Add(fHistMCtracksDecayRadiusK0s);

  fHistMCtracksDecayRadiusLambda       = new TH1F("h1MCtracksDecayRadiusLambda","Non-primary MC #Lambda^{0} Decay Radius;r (cm)",101,-1,100);
  fListHist->Add(fHistMCtracksDecayRadiusLambda);

  fHistMCtracksDecayRadiusAntiLambda   = new TH1F("h1MCtracksDecayRadiusAntiLambda","Non-primary #bar{#Lambda}^{0} Decay Radius;r (cm)",100,1,101);
  fListHist->Add(fHistMCtracksDecayRadiusAntiLambda);

  // Pt Distribution of non-primary particles:
  fHistMCPtAllK0s                      = new TH1F("h1MCPtAllK0s", "Non-primary MC K^{0};p_{t} (GeV/c);Counts",150,0,15);
  fListHist->Add(fHistMCPtAllK0s);

  fHistMCPtAllLambda                   = new TH1F("h1MCPtAllLambda", "Non-primary MC #Lambda^{0};p_{t} (GeV/c);Counts",150,0,15);
  fListHist->Add(fHistMCPtAllLambda);

  fHistMCPtAllAntiLambda               = new TH1F("h1MCPtAllAntiLambda", "Non-primary MC #bar{#Lambda}^{0};p_{t} (GeV/c);Counts",150,0,15);
  fListHist->Add(fHistMCPtAllAntiLambda);

  // Production Radius
  fHistMCProdRadiusK0s                 = new TH1F("h1MCProdRadiusK0s", "MC K^{0} Production Radius;r (cm);Count", 400, -2, 2);
  fListHist->Add(fHistMCProdRadiusK0s);

  fHistMCProdRadiusLambda              = new TH1F("h1MCProdRadiusLambda", "MC #Lambda^{0} Production Radius;r (cm);Count", 400, -2, 2);
  fListHist->Add(fHistMCProdRadiusLambda);

   fHistMCProdRadiusAntiLambda         = new TH1F("h1MCProdRadiusAntiLambda", "MC #bar{#Lambda}^{0} Production Radius;r (cm);Count", 400, -2, 2);
  fListHist->Add(fHistMCProdRadiusAntiLambda);


  // Pt and rapidity distribution:
  fHistMCPtVsYK0s               = new TH2F("h2MCPtVsYK0s", "K^{0};p_{t} (GeV/c);rapidity",150,0,15,20,-10,10);
  fListHist->Add(fHistMCPtVsYK0s);

  fHistMCPtVsYLambda            = new TH2F("h2MCPtVsYLambda", "#Lambda^{0};p_{t} (GeV/c);rapidity",150,0,15,20,-10,10);
  fListHist->Add(fHistMCPtVsYLambda);

  fHistMCPtVsYAntiLambda        = new TH2F("h2MCPtVsYAntiLambda", "#bar{#Lambda}^{0};p_{t} (GeV/c);rapidity",150,0,15,20,-10,10);
  fListHist->Add(fHistMCPtVsYAntiLambda);

  // Pt distribution of Lambda coming from Sigma decay
  fHistMCPtLambdaFromSigma      = new TH1F("h1MCPtLambdaFromSigma", "#Lambda^{0};p_{t} (GeV/c)",150,0,15);
  fListHist->Add(fHistMCPtLambdaFromSigma);

  fHistMCPtAntiLambdaFromSigma  = new TH1F("h1MCPtAntiLambdaFromSigma", "#Lambda^{0};p_{t} (GeV/c)",150,0,15);
  fListHist->Add(fHistMCPtAntiLambdaFromSigma);
 
  // Multiple reconstruction studies
  fHistNTimesRecK0s             = new TH1F("h1NTimesRecK0s","number of times a K0s is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
  fListHist->Add(fHistNTimesRecK0s);
  fHistNTimesRecK0sMI           = new TH1F("h1NTimesRecK0sMI","number of times a K0s MI is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
  fListHist->Add(fHistNTimesRecK0sMI);

  fHistNTimesRecLambda          = new TH1F("h1NTimesRecLambda","number of times a Lambda is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
  fListHist->Add(fHistNTimesRecLambda);
  fHistNTimesRecLambdaMI        = new TH1F("h1NTimesRecLambdaMI","number of times a Lambda MI is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
  fListHist->Add(fHistNTimesRecLambdaMI);

  fHistNTimesRecAntiLambda      = new TH1F("h1NTimesRecAntiLambda","number of times an AntiLambda is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
  fListHist->Add(fHistNTimesRecAntiLambda);
  fHistNTimesRecAntiLambdaMI    = new TH1F("h1NTimesRecAntiLambdaMI","number of times an AntiLambda  MI is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
  fListHist->Add(fHistNTimesRecAntiLambdaMI);

  fHistNTimesRecK0sVsPt         = new TH2F("h2NTimesRecK0sVsPt","NTimes versus Pt, K^{0} in -1<y<1;p_{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
  fListHist->Add(fHistNTimesRecK0sVsPt);
  fHistNTimesRecK0sVsPtMI       = new TH2F("h2NTimesRecK0sVsPtMI","NTimes versus Pt, K^{0}, on-the-fly finder, in -1<y<1;p_{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
  fListHist->Add(fHistNTimesRecK0sVsPtMI);

  fHistNTimesRecLambdaVsPt      = new TH2F("h2NTimesRecLambdaVsPt","NTimes versus Pt, #Lambda^{0} in -1<y<1;p_{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
  fListHist->Add(fHistNTimesRecLambdaVsPt);
  fHistNTimesRecLambdaVsPtMI    = new TH2F("h2NTimesRecLambdaVsPtMI","NTimes versus Pt, #Lambda^{0} on-the-fly finder in -1<y<1;p_{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
  fListHist->Add(fHistNTimesRecLambdaVsPtMI);

  fHistNTimesRecAntiLambdaVsPt  = new TH2F("h2NTimesRecAntiLambdaVsPt","NTimes versus Pt, #bar{#Lambda}^{0} in -1<y<1;p_{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
  fListHist->Add(fHistNTimesRecAntiLambdaVsPt);
  fHistNTimesRecAntiLambdaVsPtMI= new TH2F("h2NTimesRecAntiLambdaVsPtMI","NTimes versus Pt, #bar{#Lambda}^{0}, on-the-fly finder in -1<y<1;p_{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
  fListHist->Add(fHistNTimesRecAntiLambdaVsPtMI);

  

  //***********************************
  // Reconstructed particles histograms
  //***********************************

  // multiplicity
  fHistTrackPerEvent           = new TH1F("h1TrackPerEvent", "Tracks per event;Number of Tracks;Number of Events",50,0,50);
  fListHist->Add(fHistTrackPerEvent);

  fHistMCDaughterTrack         = new TH1F("h1MCDaughterTrack","Distribution of mc id for daughters;id tags;Counts",15,0,15);
  fListHist->Add(fHistMCDaughterTrack);

  // Primary Vertex:
  fHistPrimaryVertexX          = new TH1F("h1PrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexX);

  fHistPrimaryVertexY          = new TH1F("h1PrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexY);

  fHistPrimaryVertexZ          = new TH1F("h1PrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-5.0,5.0);
  fListHist->Add(fHistPrimaryVertexZ);

  // Cut checks:
  fHistDcaPosToPrimVertex      = new TH2F("h2DcaPosToPrimVertex", "Positive V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertex);

  fHistDcaNegToPrimVertex      = new TH2F("h2DcaNegToPrimVertex", "Negative V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertex);

  fHistDcaPosToPrimVertexZoom  = new TH2F("h2DcaPosToPrimVertexZoom", "Positive V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertexZoom);

  fHistDcaNegToPrimVertexZoom  = new TH2F("h2DcaNegToPrimVertexZoom", "Negative V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertexZoom);

  fHistRadiusV0                = new TH2F("h2RadiusV0", "Radius;Radius(cm);Status",1200,0,120,2,-0.5,1.5);
  fListHist->Add(fHistRadiusV0);

  fHistDecayLengthV0           = new TH2F("h2DecayLengthV0", "V0s decay Length;decay length(cm);Status", 240, 0, 120,2,-0.5,1.5);
  fListHist->Add(fHistDecayLengthV0);

  fHistDcaV0Daughters          = new TH2F("h2DcaV0Daughters", "DCA between daughters;dca(cm);Status", 160, 0, 4,2,-0.5,1.5);
  fListHist->Add(fHistDcaV0Daughters);

  fHistChi2                    = new TH2F("h2Chi2", "V0s chi2;chi2;Status", 33, 0, 33,2,-0.5,1.5);
  fListHist->Add(fHistChi2);

  fHistCosPointAngle           = new TH2F("h2CosPointAngle", "Cosine of V0's pointing angle", 100,0,1,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngle);

  fHistCosPointAngleZoom       = new TH2F("h2CosPointAngleZoom", "Cosine of V0's pointing angle", 100,0.9,1,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngleZoom);

  fHistProdRadius              = new TH2F("h2ProdRadius", "Production position;x (cm);y (cm)", 100,-50,50,100,-50,50);
  fListHist->Add(fHistProdRadius);

  fHistProdRadiusMI            = new TH2F("h2ProdRadiusMI", "Production position, V0s MI;x (cm);y (cm)", 100,-50,50,100,-50,50);
  fListHist->Add(fHistProdRadiusMI);

  // V0 Multiplicity
  if (!fHistV0Multiplicity) {
    if (fCollidingSystems)
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of Offline V0s;Events", 200, 0, 40000);
    else
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of Offline V0s;Events", 10, 0, 10); 
    fListHist->Add(fHistV0Multiplicity);
  }

  if (!fHistV0MultiplicityMI) {
    if (fCollidingSystems)
      fHistV0MultiplicityMI = new TH1F("fHistV0MultiplicityMI", "Multiplicity distribution;Number of On-the-fly V0s;Events", 200, 0, 40000);
    else
      fHistV0MultiplicityMI = new TH1F("fHistV0MultiplicityMI", "Multiplicity distribution;Number of On-the-fly V0s;Events", 10, 0, 10); 
    fListHist->Add(fHistV0MultiplicityMI);
  }

  // Pt and rapidity distribution:
  fHistPtVsYK0s                = new TH2F("h2PtVsYK0s", "K^{0} candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYK0s);
  fHistPtVsYK0sMI              = new TH2F("h2PtVsYK0sMI", "K^{0} MI candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYK0sMI);

  fHistPtVsYLambda             = new TH2F("h2PtVsYLambda", "#Lambda^{0} candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYLambda);
  fHistPtVsYLambdaMI           = new TH2F("h2PtVsYLambdaMI", "#Lambda^{0} MI candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYLambdaMI);

  fHistPtVsYAntiLambda         = new TH2F("h2PtVsYAntiLambda", "#bar{#Lambda}^{0} candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYAntiLambda);
  fHistPtVsYAntiLambdaMI       = new TH2F("h2PtVsYAntiLambdaMI", "#bar{#Lambda}^{0} MI candidates;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYAntiLambdaMI);

  // Mass:
  fHistMassK0                   = new TH1F("h1MassK0", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistMassK0);
  fHistMassK0MI                 = new TH1F("h1MassK0MI", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistMassK0MI);

  fHistMassLambda               = new TH1F("h1MassLambda", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassLambda);
  fHistMassLambdaMI             = new TH1F("h1MassLambdaMI", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassLambdaMI);

  fHistMassAntiLambda           = new TH1F("h1MassAntiLambda", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambda);
  fHistMassAntiLambdaMI         = new TH1F("h1MassAntiLambdaMI", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambdaMI);

  // invariant mass vs radius
  const Double_t radius[10] = {0.0,2.5,2.9,3.9,7.6,15.0,23.9,37.8,42.8,100.0};
  Int_t lNbinRadius        = 9;
  Int_t lNbinInvMassLambda = 300;

  fHistMassVsRadiusK0           = new TH2F("h2MassVsRadiusK0", "K^{0} candidates;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 200, 0.4, 0.6);
  fListHist->Add(fHistMassVsRadiusK0);

  fHistMassVsRadiusK0MI         = new TH2F("h2MassVsRadiusK0MI", "K^{0} MI candidates;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 200, 0.4, 0.6);
  fListHist->Add(fHistMassVsRadiusK0MI);
  
  fHistMassVsRadiusLambda       = new TH2F("h2MassVsRadiusLambda", "#Lambda candidates;radius (cm);M(p#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusLambda);

  fHistMassVsRadiusLambdaMI     = new TH2F("h2MassVsRadiusLambdaMI", "#Lambda MI candidates;radius (cm);M(p#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusLambdaMI);

  fHistMassVsRadiusAntiLambda   = new TH2F("h2MassVsRadiusAntiLambda", "#bar{#Lambda} candidates;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",lNbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusAntiLambda);

  fHistMassVsRadiusAntiLambdaMI = new TH2F("h2MassVsRadiusAntiLambdaMI", "#bar{#Lambda} candidates;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",lNbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusAntiLambdaMI);

  // Pt Vs Mass
  fHistPtVsMassK0               = new TH2F("h2PtVsMassK0","K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,100,0,10);
  fListHist->Add(fHistPtVsMassK0);
  fHistPtVsMassK0MI             = new TH2F("h2PtVsMassK0MI","K^{0} MIcandidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,100,0,10);
  fListHist->Add(fHistPtVsMassK0MI);

  fHistPtVsMassLambda           = new TH2F("h2PtVsMassLambda","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistPtVsMassLambda);
  fHistPtVsMassLambdaMI         = new TH2F("h2PtVsMassLambdaMI","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistPtVsMassLambdaMI);

  fHistPtVsMassAntiLambda       = new TH2F("h2PtVsMassAntiLambda","#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistPtVsMassAntiLambda);
  fHistPtVsMassAntiLambdaMI     = new TH2F("h2PtVsMassAntiLambdaMI","#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistPtVsMassAntiLambdaMI);


  fHistArmenterosPodolanski     = new TH2F("h2ArmenterosPodolanski","Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);
  fHistArmenterosPodolanskiMI   = new TH2F("h2ArmenterosPodolanskiMI","Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);


  //PID
  fHistNsigmaPosPionAntiLambda   = new TH1F("h1NsigmaPosPionAntiLambda", "Positive daughter of Antilambda;NsigmaPion;Counts",25,0,5);
  fListHist->Add(fHistNsigmaPosPionAntiLambda);

  fHistNsigmaNegProtonAntiLambda = new TH1F("h1NsigmaNegProtonAntiLambda", "Negative daughter of Antilambda;NsigmaProton;Counts",25,0,5);
  fListHist->Add(fHistNsigmaNegProtonAntiLambda);
  
  fHistNsigmaPosProtonLambda     = new TH1F("h1NsigmaPosProtonLambda", "Positive daughter of Lambda;NsigmaProton;Counts",25,0,5); 
  fListHist->Add(fHistNsigmaPosProtonLambda);
  
  fHistNsigmaNegPionLambda       = new TH1F("h1NsigmaNegPionLambda", "Negative daughter of Lambda;NsigmaPion;Counts",25,0,5);
  fListHist->Add(fHistNsigmaNegPionLambda);
  
  fHistNsigmaPosPionK0           = new TH1F("h1NsigmaPosPionK0", "Positive daughter of K0s;NsigmaPion;Counts",25,0,5);
  fListHist->Add(fHistNsigmaPosPionK0);
  
  fHistNsigmaNegPionK0           = new TH1F("h1NsigmaNegPionK0", "Negative daughter of K0s;NsigmaPion;Counts",25,0,5);
  fListHist->Add(fHistNsigmaNegPionK0);


  //********************************
  // Associated particles histograms
  //********************************

  //Pt distribution
  fHistAsMcPtK0                = new TH1F("h1AsMcPtK0", "K^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtK0);
  fHistAsMcPtK0MI              = new TH1F("h1AsMcPtK0MI", "K^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtK0MI);

  fHistAsMcPtLambda            = new TH1F("h1AsMcPtLambda", "#Lambda^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtLambda);
  fHistAsMcPtLambdaMI          = new TH1F("h1AsMcPtLambdaMI", "#Lambda^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtLambdaMI);

  fHistAsMcPtAntiLambda        = new TH1F("h1AsMcPtAntiLambda", "#bar{#Lambda}^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtAntiLambda);
  fHistAsMcPtAntiLambdaMI      = new TH1F("h1AsMcPtAntiLambdaMI", "#bar{#Lambda}^{0} associated;p_{t} (GeV/c);Counts", 150, 0, 15);
  fListHist->Add(fHistAsMcPtAntiLambdaMI);

  fHistAsMcPtZoomK0            = new TH1F("h1AsMcPtZoomK0", "K^{0} candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomK0);
  fHistAsMcPtZoomK0MI          = new TH1F("h1AsMcPtZoomK0MI", "K^{0} MI candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomK0MI);

  fHistAsMcPtZoomLambda        = new TH1F("h1AsMcPtZoomLambda", "#Lambda^{0} candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomLambda);
  fHistAsMcPtZoomLambdaMI      = new TH1F("h1AsMcPtZoomLambdaMI", "#Lambda^{0} MI candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomLambdaMI);


  // Radius distribution
  fHistAsMcProdRadiusK0               = new TH1F("h1AsMcProdRadiusK0", "K^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusK0);
  fHistAsMcProdRadiusK0MI             = new TH1F("h1AsMcProdRadiusK0MI", "K^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusK0MI);

  fHistAsMcProdRadiusLambda           = new TH1F("h1AsMcProdRadiusLambda", "#Lambda^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusLambda);
  fHistAsMcProdRadiusLambdaMI         = new TH1F("h1AsMcProdRadiusLambdaMI", "#Lambda^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusLambdaMI);

  fHistAsMcProdRadiusAntiLambda       = new TH1F("h1AsMcProdRadiusAntiLambda", "#bar{#Lambda}^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusAntiLambda);
  fHistAsMcProdRadiusAntiLambdaMI     = new TH1F("h1AsMcProdRadiusAntiLambdaMI", "#bar{#Lambda}^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusAntiLambdaMI);

  fHistAsMcProdRadiusXvsYK0s          = new TH2F("h2AsMcProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYK0s);
  fHistAsMcProdRadiusXvsYK0sMI        = new TH2F("h2AsMcProdRadiusXvsYK0sMI","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYK0sMI);

  fHistAsMcProdRadiusXvsYLambda       = new TH2F("h2AsMcProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYLambda);
  fHistAsMcProdRadiusXvsYLambdaMI     = new TH2F("h2AsMcProdRadiusXvsYLambdaMI","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYLambdaMI);

  fHistAsMcProdRadiusXvsYAntiLambda   = new TH2F("h2AsMcProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYAntiLambda);
  fHistAsMcProdRadiusXvsYAntiLambdaMI = new TH2F("h2AsMcProdRadiusXvsYAntiLambdaMI","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYAntiLambdaMI);



  // Mass
  fHistPidMcMassK0             = new TH1F("h1PidMcMassK0", "K^{0} MC PId checked;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistPidMcMassK0);
  fHistPidMcMassK0MI           = new TH1F("h1PidMcMassK0MI", "K^{0} MC PId checked;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistPidMcMassK0MI);

  fHistPidMcMassLambda         = new TH1F("h1PidMcMassLambda", "#Lambda^{0} MC PId checked;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassLambda);
  fHistPidMcMassLambdaMI       = new TH1F("h1PidMcMassLambdaMI", "#Lambda^{0} MC PId checked;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassLambdaMI);
  
  fHistPidMcMassAntiLambda     = new TH1F("h1PidMcMassAntiLambda", "#bar{#Lambda}^{0} MC PId checked;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassAntiLambda);
  fHistPidMcMassAntiLambdaMI   = new TH1F("h1PidMcMassAntiLambdaMI", "#bar{#Lambda}^{0} MC PId checked;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassAntiLambdaMI);

  fHistAsMcMassK0              = new TH1F("h1AsMcMassK0", "K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistAsMcMassK0);
  fHistAsMcMassK0MI            = new TH1F("h1AsMcMassK0MI", "K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistAsMcMassK0MI);
  
  fHistAsMcMassLambda          = new TH1F("h1AsMcMassLambda", "#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassLambda);
  fHistAsMcMassLambdaMI        = new TH1F("h1AsMcMassLambdaMI", "#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassLambdaMI);

  fHistAsMcMassAntiLambda      = new TH1F("h1AsMcMassAntiLambda", "#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassAntiLambda);
  fHistAsMcMassAntiLambdaMI    = new TH1F("h1AsMcMassAntiLambdaMI", "#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassAntiLambdaMI);

  //Pt versus Mass
  fHistAsMcPtVsMassK0               = new TH2F("h2AsMcPtVsMassK0","K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,100,0,10);
  fListHist->Add(fHistAsMcPtVsMassK0);
  fHistAsMcPtVsMassK0MI             = new TH2F("h2AsMcPtVsMassK0MI","K^{0} MIassociated;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,100,0,10);
  fListHist->Add(fHistAsMcPtVsMassK0MI);

  fHistAsMcPtVsMassLambda           = new TH2F("h2AsMcPtVsMassLambda","#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistAsMcPtVsMassLambda);
  fHistAsMcPtVsMassLambdaMI         = new TH2F("h2AsMcPtVsMassLambdaMI","#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistAsMcPtVsMassLambdaMI);

  fHistAsMcPtVsMassAntiLambda       = new TH2F("h2AsMcPtVsMassAntiLambda","#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistAsMcPtVsMassAntiLambda);
  fHistAsMcPtVsMassAntiLambdaMI     = new TH2F("h2AsMcPtVsMassAntiLambdaMI","#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistAsMcPtVsMassAntiLambdaMI);


  // invariant mass vs radius
  fHistAsMcMassVsRadiusK0             = new TH2F("h2AsMcMassVsRadiusK0", "K^{0} associated;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 500, 0.47, 0.52);
  fListHist->Add(fHistAsMcMassVsRadiusK0);

  fHistAsMcMassVsRadiusK0MI           = new TH2F("h2AsMcMassVsRadiusK0MI", "K^{0} MI associated;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 500, 0.47, 0.52);
  fListHist->Add(fHistAsMcMassVsRadiusK0MI);
  
  fHistAsMcMassVsRadiusLambda         = new TH2F("h2AsMcMassVsRadiusLambda", "#Lambda associated;radius (cm);M(p#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, lNbinInvMassLambda, 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusLambda);

  fHistAsMcMassVsRadiusLambdaMI       = new TH2F("h2AsMcMassVsRadiusLambdaMI", "#Lambda MI associated;radius (cm);M(p#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, lNbinInvMassLambda, 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusLambdaMI);

  fHistAsMcMassVsRadiusAntiLambda     = new TH2F("h2AsMcMassVsRadiusAntiLambda", "#bar{#Lambda} associated;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",lNbinRadius,radius,lNbinInvMassLambda , 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusAntiLambda);
  
  fHistAsMcMassVsRadiusAntiLambdaMI   = new TH2F("h2AsMcMassVsRadiusAntiLambdaMI", "#bar{#Lambda} MI associated;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",lNbinRadius,radius,lNbinInvMassLambda , 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusAntiLambdaMI);
    

  // Position Resolution
  fHistAsMcResxK0                     = new TH1F("h1AsMcResxK0", "K^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxK0);
  fHistAsMcResyK0                     = new TH1F("h1AsMcResyK0", "K^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyK0);
  fHistAsMcReszK0                     = new TH1F("h1AsMcReszK0", "K^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszK0);
  fHistAsMcResrVsRadiusK0             = new TH2F("h2AsMcResrVsRadiusK0", "K^{0} associated;Radius (cm);#Delta r (cm)",200,0.0,50., 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusK0);
  fHistAsMcReszVsRadiusK0             = new TH2F("h2AsMcReszVsRadiusK0", "K^{0} associated;Radius (cm);#Delta z (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusK0);

  fHistAsMcResxK0MI                   = new TH1F("h1AsMcResxK0MI", "K^{0} MI associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxK0MI);
  fHistAsMcResyK0MI                   = new TH1F("h1AsMcResyK0MI", "K^{0} MI associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyK0MI);
  fHistAsMcReszK0MI                   = new TH1F("h1AsMcReszK0MI", "K^{0} MI associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszK0MI);
  fHistAsMcResrVsRadiusK0MI           = new TH2F("h2AsMcResrVsRadiusK0MI", "K^{0} MI associated;Radius (cm);#Delta r (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusK0MI);
  fHistAsMcReszVsRadiusK0MI           = new TH2F("h2AsMcReszVsRadiusK0MI", "K^{0} MI associated;Radius (cm);#Delta z (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusK0MI);

  fHistAsMcResxLambda                 = new TH1F("h1AsMcResxLambda", "#Lambda^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxLambda);
  fHistAsMcResyLambda                 = new TH1F("h1AsMcResyLambda", "#Lambda^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyLambda);
  fHistAsMcReszLambda                 = new TH1F("h1AsMcReszLambda", "#Lambda^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszLambda);
  fHistAsMcResrVsRadiusLambda         = new TH2F("h2AsMcResrVsRadiusLambda", "#Lambda^{0} associated;Radius (cm);#Delta r (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusLambda);
  fHistAsMcReszVsRadiusLambda         = new TH2F("h2AsMcReszVsRadiusLambda", "#Lambda^{0} associated;Radius (cm);#Delta z (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusLambda);

  fHistAsMcResxLambdaMI               = new TH1F("h1AsMcResxLambdaMI", "#Lambda^{0} MI associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxLambdaMI);
  fHistAsMcResyLambdaMI               = new TH1F("h1AsMcResyLambdaMI", "#Lambda^{0} MI associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyLambdaMI);
  fHistAsMcReszLambdaMI               = new TH1F("h1AsMcReszLambdaMI", "#Lambda^{0} MI associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszLambdaMI);
  fHistAsMcResrVsRadiusLambdaMI       = new TH2F("h2AsMcResrVsRadiusLambdaMI", "#Lambda^{0} MI associated;Radius (cm);#Delta r (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusLambdaMI);
  fHistAsMcReszVsRadiusLambdaMI       = new TH2F("h2AsMcReszVsRadiusLambdaMI", "#Lambda^{0} MI associated;Radius (cm);#Delta z (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusLambdaMI);

  fHistAsMcResxAntiLambda             = new TH1F("h1AsMcResxAntiLambda", "#bar{#Lambda}^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxAntiLambda);
  fHistAsMcResyAntiLambda             = new TH1F("h1AsMcResyAntiLambda", "#bar{#Lambda}^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyAntiLambda);
  fHistAsMcReszAntiLambda             = new TH1F("h1AsMcReszAntiLambda", "#bar{#Lambda}^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszAntiLambda);
  fHistAsMcResrVsRadiusAntiLambda     = new TH2F("h2AsMcResrVsRadiusAntiLambda", "#bar{#Lambda}^{0} associated;Radius (cm);#Delta r (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusAntiLambda);
  fHistAsMcReszVsRadiusAntiLambda     = new TH2F("h2AsMcReszVsRadiusAntiLambda", "#bar{#Lambda}^{0} associated;Radius (cm);#Delta z (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusAntiLambda);

  fHistAsMcResxAntiLambdaMI           = new TH1F("h1AsMcResxAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxAntiLambdaMI);
  fHistAsMcResyAntiLambdaMI           = new TH1F("h1AsMcResyAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyAntiLambdaMI);
  fHistAsMcReszAntiLambdaMI           = new TH1F("h1AsMcReszAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszAntiLambdaMI);
  fHistAsMcResrVsRadiusAntiLambdaMI   = new TH2F("h2AsMcResrVsRadiusAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;Radius (cm);#Delta r (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusAntiLambdaMI);
  fHistAsMcReszVsRadiusAntiLambdaMI   = new TH2F("h2AsMcReszVsRadiusAntiLambdaMI", "#bar{#Lambda}^{0} MI associated;Radius (cm);#Delta z (cm)",8,radius, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusAntiLambdaMI);

  // Pt Resolution
  fHistAsMcResPtK0                   = new TH1F("h1AsMcResPtK0","Pt Resolution K^{0};#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtK0);
  fHistAsMcResPtK0MI                 = new TH1F("h1AsMcResPtK0MI","Pt Resolution K^{0} MI;#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtK0MI);
  
  fHistAsMcResPtLambda               = new TH1F("h1AsMcResPtLambda","Pt Resolution #Lambda^{0};#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtLambda);
  fHistAsMcResPtLambdaMI             = new TH1F("h1AsMcResPtLambdaMI","Pt Resolution #Lambda^{0} MI;#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtLambdaMI);

  fHistAsMcResPtAntiLambda           = new TH1F("h1AsMcResPtAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtAntiLambda);
  fHistAsMcResPtAntiLambdaMI         = new TH1F("h1AsMcResPtAntiLambdaMI","Pt Resolution #bar{#Lambda}^{0} MI;#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtAntiLambdaMI);


  fHistAsMcResPtVsRapK0              = new TH2F("h2AsMcResPtVsRapK0","Pt Resolution K^{0};#Delta Pt;Rapidity",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapK0);
  fHistAsMcResPtVsRapK0MI            = new TH2F("h2AsMcResPtVsRapK0MI","Pt Resolution K^{0} MI;#Delta Pt;Rapidity",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapK0MI);
  
  fHistAsMcResPtVsRapLambda          = new TH2F("h2AsMcResPtVsRapLambda","Pt Resolution #Lambda^{0};#Delta Pt;Rapidity",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapLambda);
  fHistAsMcResPtVsRapLambdaMI        = new TH2F("h2AsMcResPtVsRapLambdaMI","Pt Resolution #Lambda^{0} MI;#Delta Pt;Rapidity",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapLambdaMI);

  fHistAsMcResPtVsRapAntiLambda      = new TH2F("h2AsMcResPtVsRapAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Rapidity",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapAntiLambda);
  fHistAsMcResPtVsRapAntiLambdaMI    = new TH2F("h2AsMcResPtVsRapAntiLambdaMI","Pt Resolution #bar{#Lambda}^{0} MI;#Delta Pt;Rapidity",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapAntiLambdaMI);


  fHistAsMcResPtVsPtK0               = new TH2F("h2AsMcResPtVsPtK0","Pt Resolution K^{0};#Delta Pt;Pt",600,-0.15,0.15,200,0,10);
  fListHist->Add(fHistAsMcResPtVsPtK0);
  fHistAsMcResPtVsPtK0MI             = new TH2F("h2AsMcResPtVsPtK0MI","Pt Resolution K^{0} MI;#Delta Pt;Pt",600,-0.15,0.15,200,0,10);
  fListHist->Add(fHistAsMcResPtVsPtK0MI);
    
  fHistAsMcResPtVsPtLambda           = new TH2F("h2AsMcResPtVsPtLambda","Pt Resolution #Lambda^{0};#Delta Pt;Pt",600,-0.15,0.15,200,0,10);
  fListHist->Add(fHistAsMcResPtVsPtLambda);
  fHistAsMcResPtVsPtLambdaMI         = new TH2F("h2AsMcResPtVsPtLambdaMI","Pt Resolution #Lambda^{0} MI;#Delta Pt;Pt",600,-0.15,0.15,200,0,10);
  fListHist->Add(fHistAsMcResPtVsPtLambdaMI);

  fHistAsMcResPtVsPtAntiLambda       = new TH2F("h2AsMcResPtVsPtAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Pt",300,-0.15,0.15,200,0,10);
  fListHist->Add(fHistAsMcResPtVsPtAntiLambda);
  fHistAsMcResPtVsPtAntiLambdaMI     = new TH2F("h2AsMcResPtVsPtAntiLambdaMI","Pt Resolution #bar{#Lambda}^{0} MI;#Delta Pt;Pt",300,-0.15,0.15,200,0,10);
  fListHist->Add(fHistAsMcResPtVsPtAntiLambdaMI);


  // pdgcode of mother
  fHistAsMcMotherPdgCodeK0s           = new TH1F("h1AsMcMotherPdgCodeK0s","Mother of Associated K^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeK0s);
  fHistAsMcMotherPdgCodeK0sMI         = new TH1F("h1AsMcMotherPdgCodeK0sMI","Mother of Associated K^{0} MI;mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeK0sMI);

  fHistAsMcMotherPdgCodeLambda        = new TH1F("h1AsMcMotherPdgCodeLambda","Mother of Associated #Lambda^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeLambda);
  fHistAsMcMotherPdgCodeLambdaMI      = new TH1F("h1AsMcMotherPdgCodeLambdaMI","Mother of Associated #Lambda^{0} MI;mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeLambdaMI);

  fHistAsMcMotherPdgCodeAntiLambda    = new TH1F("h1AsMcMotherPdgCodeAntiLambda","Mother of Associated #bar{#Lambda}^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeAntiLambda);
  fHistAsMcMotherPdgCodeAntiLambdaMI  = new TH1F("h1AsMcMotherPdgCodeAntiLambdaMI","Mother of Associated #bar{Lambda}^{0} MI;mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeAntiLambdaMI);


  // Pt distribution Lambda from Sigma
  fHistAsMcPtLambdaFromSigma          = new TH1F("h1AsMcPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",150,0,15);
  fListHist->Add(fHistAsMcPtLambdaFromSigma);
  fHistAsMcPtLambdaFromSigmaMI        = new TH1F("h1AsMcPtLambdaFromSigmaMI","#Lambda^{0} MI associated from Sigma;p_{t} (GeV/c);Count",150,0,15);
  fListHist->Add(fHistAsMcPtLambdaFromSigmaMI);

  fHistAsMcPtAntiLambdaFromSigma      = new TH1F("h1AsMcPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",150,0,15);
  fListHist->Add(fHistAsMcPtAntiLambdaFromSigma);
  fHistAsMcPtAntiLambdaFromSigmaMI    = new TH1F("h1AsMcPtAntiLambdaFromSigmaMI","#bar{#Lambda}^{0} MI associated from Sigma;p_{t} (GeV/c);Count",150,0,15);
  fListHist->Add(fHistAsMcPtAntiLambdaFromSigmaMI);


  // Associated secondary particles:
  // Pt and rapidity distribution
  fHistAsMcSecondaryPtVsYK0s          = new TH2F("h2AsMcSecondaryPtVsYK0s", "K^{0} associated secondary;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistAsMcSecondaryPtVsYK0s);
  fHistAsMcSecondaryPtVsYK0sMI        = new TH2F("h2AsMcSecondaryPtVsYK0sMI", "K^{0} MI associated secondary;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistAsMcSecondaryPtVsYK0sMI);

  fHistAsMcSecondaryPtVsYLambda       = new TH2F("h2AsMcSecondaryPtVsYLambda", "#Lambda^{0} associated secondary;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistAsMcSecondaryPtVsYLambda);
  fHistAsMcSecondaryPtVsYLambdaMI     = new TH2F("h2AsMcSecondaryPtVsYLambdaMI", "#Lambda^{0} MI associated secondary;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistAsMcSecondaryPtVsYLambdaMI);

  fHistAsMcSecondaryPtVsYAntiLambda   = new TH2F("h2AsMcSecondaryPtVsYAntiLambda", "#bar{#Lambda}^{0} associated secondary;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistAsMcSecondaryPtVsYAntiLambda);
  fHistAsMcSecondaryPtVsYAntiLambdaMI = new TH2F("h2AsMcSecondaryPtVsYAntiLambdaMI", "#bar{#Lambda}^{0} MI associated secondary;p_{t} (GeV/c);rapidity",150,0,15,30,-1.5,1.5);
  fListHist->Add(fHistAsMcSecondaryPtVsYAntiLambdaMI);

  // Production radius
  fHistAsMcSecondaryProdRadiusK0s              = new TH1F("h1AsMcSecondaryProdRadiusK0s", "K^{0} Production Radius;r (cm);Count", 170, -2, 15);
  fListHist->Add(fHistAsMcSecondaryProdRadiusK0s);
  fHistAsMcSecondaryProdRadiusK0sMI            = new TH1F("h1AsMcSecondaryProdRadiusK0sMI", "K^{0} MI Production Radius;r (cm);Count", 170, -2, 15);
  fListHist->Add(fHistAsMcSecondaryProdRadiusK0sMI);

  fHistAsMcSecondaryProdRadiusLambda           = new TH1F("h1AsMcSecondaryProdRadiusLambda", "#Lambda^{0} Production Radius;r (cm);Count", 170, -2, 15);
  fListHist->Add(fHistAsMcSecondaryProdRadiusLambda);
  fHistAsMcSecondaryProdRadiusLambdaMI         = new TH1F("h1AsMcSecondaryProdRadiusLambdaMI", "#Lambda^{0} MI Production Radius;r (cm);Count", 170, -2, 15);
  fListHist->Add(fHistAsMcSecondaryProdRadiusLambdaMI);

  fHistAsMcSecondaryProdRadiusAntiLambda       = new TH1F("h1AsMcSecondaryProdRadiusAntiLambda", "#bar{#Lambda}^{0} Production Radius;r (cm);Count", 170, -2, 15);
  fListHist->Add(fHistAsMcSecondaryProdRadiusAntiLambda);  
  fHistAsMcSecondaryProdRadiusAntiLambdaMI     = new TH1F("h1AsMcSecondaryProdRadiusAntiLambdaMI", "#bar{#Lambda}^{0} MI Production Radius;r (cm);Count", 170, -2, 15);
  fListHist->Add(fHistAsMcSecondaryProdRadiusAntiLambdaMI);

  fHistAsMcSecondaryProdRadiusXvsYK0s          = new TH2F("h2AsMcSecondaryProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYK0s);
  fHistAsMcSecondaryProdRadiusXvsYK0sMI        = new TH2F("h2AsMcSecondaryProdRadiusXvsYK0sMI","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYK0sMI);

  fHistAsMcSecondaryProdRadiusXvsYLambda       = new TH2F("h2AsMcSecondaryProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYLambda);
  fHistAsMcSecondaryProdRadiusXvsYLambdaMI     = new TH2F("h2AsMcSecondaryProdRadiusXvsYLambdaMI","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYLambdaMI);

  fHistAsMcSecondaryProdRadiusXvsYAntiLambda   = new TH2F("h2AsMcSecondaryProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYAntiLambda);
  fHistAsMcSecondaryProdRadiusXvsYAntiLambdaMI = new TH2F("h2AsMcSecondaryProdRadiusXvsYAntiLambdaMI","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYAntiLambdaMI);

  fHistAsMcSecondaryMotherPdgCodeK0s           = new TH1F("h1AsMcSecondaryMotherPdgCodeK0s","Mother of Associated Secondary K^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeK0s);
  fHistAsMcSecondaryMotherPdgCodeK0sMI         = new TH1F("h1AsMcSecondaryMotherPdgCodeK0sMI","Mother of Associated Secondary K^{0} MI;mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeK0sMI);

  fHistAsMcSecondaryMotherPdgCodeLambda        = new TH1F("h1AsMcSecondaryMotherPdgCodeLambda","Mother of Associated Secondary #Lambda^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeLambda);
  fHistAsMcSecondaryMotherPdgCodeLambdaMI      = new TH1F("h1AsMcSecondaryMotherPdgCodeLambdaMI","Mother of Associated Secondary #Lambda^{0} MI;mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeLambdaMI);

  fHistAsMcSecondaryMotherPdgCodeAntiLambda    = new TH1F("h1AsMcSecondaryMotherPdgCodeAntiLambda","Mother of Associated Secondary #bar{#Lambda}^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeAntiLambda);
  fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI  = new TH1F("h1AsMcSecondaryMotherPdgCodeAntiLambdaMI","Mother of Associated Secondary #bar{Lambda}^{0} MI;mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI);

  // Pt distribution Lambda from Sigma
  fHistAsMcSecondaryPtLambdaFromSigma          = new TH1F("h1AsMcSecondaryPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",150,0,15);
  fListHist->Add(fHistAsMcSecondaryPtLambdaFromSigma);
  fHistAsMcSecondaryPtLambdaFromSigmaMI        = new TH1F("h1AsMcSecondaryPtLambdaFromSigmaMI","#Lambda^{0} MI associated from Sigma;p_{t} (GeV/c);Count",150,0,15);
  fListHist->Add(fHistAsMcSecondaryPtLambdaFromSigmaMI);

  fHistAsMcSecondaryPtAntiLambdaFromSigma      = new TH1F("h1AsMcSecondaryPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",150,0,15);
  fListHist->Add(fHistAsMcSecondaryPtAntiLambdaFromSigma);
  fHistAsMcSecondaryPtAntiLambdaFromSigmaMI    = new TH1F("h1AsMcSecondaryPtAntiLambdaFromSigmaMI","#bar{#Lambda}^{0} MI associated from Sigma;p_{t} (GeV/c);Count",150,0,15);
  fListHist->Add(fHistAsMcSecondaryPtAntiLambdaFromSigmaMI);

}

//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  AliStack*     stack   = 0;
  TClonesArray *mcArray = 0;

  AliVEvent* lEvent = InputEvent();
 
  if (!lEvent) {
    Printf("ERROR: Event not available");
    return;
  }
  if (!(lEvent->GetNumberOfTracks())) {
    //Printf("Strange analysis task: There is no track in this event");
    return;
  }

  fHistTrackPerEvent->Fill(lEvent->GetNumberOfTracks());

  
  if(fAnalysisType == "ESD") {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      //Printf("ERROR: Could not retrieve MC event handler");
      return;
    }    
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      //Printf("ERROR: Could not retrieve MC event");
      return;
    }    
    stack = mcEvent->Stack();
    if (!stack) {
      //Printf("ERROR: Could not retrieve stack");
      return;
   }
  }

  else if(fAnalysisType == "AOD") {

    // load MC particles
    mcArray = (TClonesArray*)lEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      Printf("strange analysis::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    AliAODMCHeader *mcHeader = 
      (AliAODMCHeader*)lEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      Printf("strange analysis::UserExec: MC header branch not found!\n");
      return;
    }
  }

  //**********************************************
  // MC loop
  //**********************************************

  Int_t lNbMCPrimary        = 0;
  Int_t lNbMCPart           = 0;

  Int_t lPdgcodeCurrentPart = 0;
  Double_t lRapCurrentPart  = 0;
  Double_t lPtCurrentPart   = 0;
  
  Int_t lComeFromSigma      = 0;

  Double_t lMaxMcProdRadiusPrimaries = 0.2;
  
  // Production Radius
  Double_t mcPosX     = 0.0,  mcPosY      = 0.0,  mcPosZ      = 0.0;
  Double_t mcPosR     = 0.0;

  // Decay Radius
  Double_t mcDecayPosX = 0, mcDecayPosY = 0, mcDecayPosR = 0;

  // current mc particle 's mother
  Int_t iCurrentMother  = 0, lPdgCurrentMother    = 0;
  Bool_t lCurrentMotherIsPrimary;

  // variables for multiple reconstruction studies:
  Int_t id0           = 0, id1          = 0;
  Int_t lLabelTrackN  = 0, lLabelTrackP = 0;
  Int_t lPartNMother  = 0, lPartPMother = 0;
  Int_t lPartPMotherPDGcode      = 0;
  Int_t lNtimesReconstructedK0s   = 0, lNtimesReconstructedLambda   = 0, lNtimesReconstructedAntiLambda   = 0;
  Int_t lNtimesReconstructedK0sMI = 0, lNtimesReconstructedLambdaMI = 0, lNtimesReconstructedAntiLambdaMI = 0;


  if(fAnalysisType == "ESD") {

    lNbMCPrimary = stack->GetNprimary();
    lNbMCPart    = stack->GetNtrack();

    fHistMCMultiplicityPrimary->Fill(lNbMCPrimary);
    fHistMCMultiplicityTracks->Fill(lNbMCPart);


    for (Int_t iMc = 0; iMc < (stack->GetNtrack()); iMc++) {  
      TParticle *p0 = stack->Particle(iMc);
      if (!p0) {
	//Printf("ERROR: particle with label %d not found in stack (mc loop)", iMc);
	continue;
      }
      lPdgcodeCurrentPart = p0->GetPdgCode();

      // Keep only K0s, Lambda and AntiLambda:
      if ( (lPdgcodeCurrentPart != 310 ) && (lPdgcodeCurrentPart != 3122 ) && (lPdgcodeCurrentPart != -3122 ) ) continue;

      lRapCurrentPart   = MyRapidity(p0->Energy(),p0->Pz());
      lPtCurrentPart    = p0->Pt();
      iCurrentMother    = p0->GetFirstMother();
      lPdgCurrentMother = stack->Particle(iCurrentMother)->GetPdgCode();

      mcPosX = p0->Vx();
      mcPosY = p0->Vy();
      mcPosZ = p0->Vz();
      mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);

      id0  = p0->GetDaughter(0);
      id1  = p0->GetDaughter(1);

       // Decay Radius and Production Radius
      if ( id0 <= lNbMCPart && id0 > 0 && id1 <= lNbMCPart && id1 > 0) {
	TParticle *pDaughter0 = stack->Particle(id0);
	mcDecayPosX = pDaughter0->Vx();
	mcDecayPosY = pDaughter0->Vy();
	mcDecayPosR = TMath::Sqrt(mcDecayPosX*mcDecayPosX+mcDecayPosY*mcDecayPosY);
      }
      else  {
	//Printf("ERROR: particle with label %d and/or %d not found in stack (mc loop)", id0,id1);
	mcDecayPosR = -1.0;
      }

      if (lPdgcodeCurrentPart==310)   {
	fHistMCtracksProdRadiusK0s->Fill(mcPosX,mcPosY);
	fHistMCtracksDecayRadiusK0s->Fill(mcDecayPosR);
	if (TMath::Abs(lRapCurrentPart) < 1) fHistMCPtAllK0s->Fill(lPtCurrentPart);
      }
      else if (lPdgcodeCurrentPart==3122)  {
	fHistMCtracksProdRadiusLambda->Fill(mcPosX,mcPosY);
	fHistMCtracksDecayRadiusLambda->Fill(mcDecayPosR);
	if (TMath::Abs(lRapCurrentPart) < 1) fHistMCPtAllLambda->Fill(lPtCurrentPart);
      }
      else if (lPdgcodeCurrentPart==-3122) {
	fHistMCtracksProdRadiusAntiLambda->Fill(mcPosX,mcPosY);
	fHistMCtracksDecayRadiusAntiLambda->Fill(mcDecayPosR);
	if (TMath::Abs(lRapCurrentPart) < 1) fHistMCPtAllAntiLambda->Fill(lPtCurrentPart);
      }

      if ( ( ( TMath::Abs(lPdgCurrentMother) == 3212)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3222)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3112)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3224)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3214)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3114) )
	   && ( iCurrentMother <= lNbMCPrimary )
	   ) lComeFromSigma = 1;
      else lComeFromSigma = 0;

      //*********************************************
      // Now keep only primary particles      
      if ( mcPosR > lMaxMcProdRadiusPrimaries ) continue;

      //********************************************
      // check if V0 is reconstructed several times  
     
      lNtimesReconstructedK0s   = 0; lNtimesReconstructedLambda   = 0; lNtimesReconstructedAntiLambda   = 0;
      lNtimesReconstructedK0sMI = 0; lNtimesReconstructedLambdaMI = 0; lNtimesReconstructedAntiLambdaMI = 0;

      for (Int_t jV0 = 0; jV0 < lEvent->GetNumberOfV0s(); jV0++) {

	lLabelTrackN  = 0; lLabelTrackP = 0;
	lPartNMother  = 0; lPartPMother = 0;
	
	AliESDv0    *vertexESD = ((AliESDEvent*)lEvent)->GetV0(jV0);
	if (!vertexESD) continue;
	
	AliESDtrack *trackNESD = ((AliESDEvent*)lEvent)->GetTrack(TMath::Abs(vertexESD->GetNindex()));
	lLabelTrackN = (UInt_t)TMath::Abs(trackNESD->GetLabel());
	if (lLabelTrackN!=id0 && lLabelTrackN!=id1) continue;
	
	AliESDtrack *trackPESD = ((AliESDEvent*)lEvent)->GetTrack(TMath::Abs(vertexESD->GetPindex()));
	lLabelTrackP = (UInt_t)TMath::Abs(trackPESD->GetLabel());
	if (lLabelTrackP!=id0 && lLabelTrackP!=id1) continue;
	
	TParticle   *lPartNESD = stack->Particle(lLabelTrackN);
	TParticle   *lPartPESD = stack->Particle(lLabelTrackP);
	lPartNMother = lPartNESD->GetFirstMother();
	lPartPMother = lPartPESD->GetFirstMother();

	lPartPMotherPDGcode = stack->Particle(lPartPMother)->GetPdgCode();
	
	switch (vertexESD->GetOnFlyStatus()){
	  
	case 0 : 
	  if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==310) ) lNtimesReconstructedK0s++;
	  else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==3122) ) lNtimesReconstructedLambda++;
	  else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==-3122) ) lNtimesReconstructedAntiLambda++;
	  break;
	  
	case 1 :
	  if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==310) ) lNtimesReconstructedK0sMI++;
	  else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==3122) ) lNtimesReconstructedLambdaMI++;
	  else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==-3122) ) lNtimesReconstructedAntiLambdaMI++;
	  break;
	  
	}	
      } // end loop over reconstructed V0s inside MC loop
      
      if (lPdgcodeCurrentPart==310) {
	fHistMCProdRadiusK0s->Fill(mcPosR);
	fHistMCPtVsYK0s->Fill(lPtCurrentPart,lRapCurrentPart);
	if (TMath::Abs(lRapCurrentPart) > 1) continue;
	fHistNTimesRecK0s->Fill(lNtimesReconstructedK0s);
	fHistNTimesRecK0sMI->Fill(lNtimesReconstructedK0s);
	fHistNTimesRecK0sVsPt->Fill(lPtCurrentPart,lNtimesReconstructedK0s);
	fHistNTimesRecK0sVsPtMI->Fill(lPtCurrentPart,lNtimesReconstructedK0sMI);
      }
      else if (lPdgcodeCurrentPart==3122) {
	fHistMCProdRadiusLambda->Fill(mcPosR);
	fHistMCPtVsYLambda->Fill(lPtCurrentPart,lRapCurrentPart);
	if (TMath::Abs(lRapCurrentPart) > 1) continue;
	fHistNTimesRecLambda->Fill(lNtimesReconstructedLambda);
	fHistNTimesRecLambdaMI->Fill(lNtimesReconstructedLambdaMI);
	fHistNTimesRecLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedLambda);
	fHistNTimesRecLambdaVsPtMI->Fill(lPtCurrentPart,lNtimesReconstructedLambdaMI);
	if (lComeFromSigma) fHistMCPtLambdaFromSigma->Fill(lPtCurrentPart);
      }
      else if (lPdgcodeCurrentPart==-3122) {
	fHistMCProdRadiusAntiLambda->Fill(mcPosR);
	fHistMCPtVsYAntiLambda->Fill(lPtCurrentPart,lRapCurrentPart);
	if (TMath::Abs(lRapCurrentPart) > 1) continue;
	fHistNTimesRecAntiLambda->Fill(lNtimesReconstructedAntiLambda);
	fHistNTimesRecAntiLambdaMI->Fill(lNtimesReconstructedAntiLambdaMI);
	fHistNTimesRecAntiLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedAntiLambda);
	fHistNTimesRecAntiLambdaVsPtMI->Fill(lPtCurrentPart,lNtimesReconstructedAntiLambdaMI);
	if (lComeFromSigma) fHistMCPtAntiLambdaFromSigma->Fill(lPtCurrentPart);
      }

    } // end loop ESD MC

  } // end ESD condition

  else if(fAnalysisType == "AOD") {
    lNbMCPart = mcArray->GetEntriesFast();
    lNbMCPrimary = 0;

    fHistMCMultiplicityTracks->Fill(lNbMCPart);

    for (Int_t iMc = 0; iMc < lNbMCPart; iMc++) {  

      AliAODMCParticle *mcAODPart = (AliAODMCParticle*)mcArray->At(iMc);
      if (!mcAODPart) {
	//Printf("Strange analysis task (mc loop): particle with label %d not found", iMc);
	continue;
      }
      lPdgcodeCurrentPart = mcAODPart->GetPdgCode();
      if (mcAODPart->IsPhysicalPrimary()) {lNbMCPrimary = lNbMCPrimary +1;}

      // Keep only K0s, Lambda and AntiLambda:
      if ( (lPdgcodeCurrentPart != 310 ) && (lPdgcodeCurrentPart != 3122 ) && (lPdgcodeCurrentPart != -3122 ) ) continue;

      lRapCurrentPart   = mcAODPart->Y();
      lPtCurrentPart    = mcAODPart->Pt();
      iCurrentMother    = mcAODPart->GetMother();
      lPdgCurrentMother = ((AliAODMCParticle*)mcArray->At(iCurrentMother))->GetPdgCode();
      lCurrentMotherIsPrimary = ((AliAODMCParticle*)mcArray->At(iCurrentMother))->IsPhysicalPrimary();

      mcPosX = mcAODPart->Xv();
      mcPosY = mcAODPart->Yv();
      mcPosZ = mcAODPart->Zv();
      mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
      
      id0  = mcAODPart->GetDaughter(0);
      id1  = mcAODPart->GetDaughter(1);

      // Decay Radius and Production Radius
      if ( id0 <= lNbMCPart && id0 > 0 && id1 <= lNbMCPart && id1 > 0) {
	AliAODMCParticle *mcAODDaughter1 = (AliAODMCParticle*)mcArray->At(id1);
	if (!mcAODPart) {
	  //Printf("Strange analysis task (mc loop): daughter not found");
	  continue;
	}
	mcDecayPosX = mcAODDaughter1->Xv();
	mcDecayPosY = mcAODDaughter1->Yv();
	mcDecayPosR = TMath::Sqrt(mcDecayPosX*mcDecayPosX+mcDecayPosY*mcDecayPosY);
      }
      else  {
	//Printf("ERROR: particle with label %d and/or %d not found in stack (mc loop)", id0,id1);
	mcDecayPosR = -1.0;
      }

      if (lPdgcodeCurrentPart==310)   {
	fHistMCtracksProdRadiusK0s->Fill(mcPosX,mcPosY);
	fHistMCtracksDecayRadiusK0s->Fill(mcDecayPosR);
	if (TMath::Abs(lRapCurrentPart) < 1) fHistMCPtAllK0s->Fill(lPtCurrentPart);
      }
      else if (lPdgcodeCurrentPart==3122)  {
	fHistMCtracksProdRadiusLambda->Fill(mcPosX,mcPosY);
	fHistMCtracksDecayRadiusLambda->Fill(mcDecayPosR);
	if (TMath::Abs(lRapCurrentPart) < 1) fHistMCPtAllLambda->Fill(lPtCurrentPart);
      }
      else if (lPdgcodeCurrentPart==-3122) {
	fHistMCtracksProdRadiusAntiLambda->Fill(mcPosX,mcPosY);
	fHistMCtracksDecayRadiusAntiLambda->Fill(mcDecayPosR);
	if (TMath::Abs(lRapCurrentPart) < 1) fHistMCPtAllAntiLambda->Fill(lPtCurrentPart);
      }

      if ( ( ( TMath::Abs(lPdgCurrentMother) == 3212)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3222)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3112)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3224)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3214)  ||
	     ( TMath::Abs(lPdgCurrentMother) == 3114) )
	   && (lCurrentMotherIsPrimary)
	   ) lComeFromSigma = 1;
      else lComeFromSigma = 0;
      
      //*********************************************
      // Now keep only primary particles      
      if ( mcPosR > lMaxMcProdRadiusPrimaries ) continue;

      //********************************************
      // check if V0 is reconstructed several times  
     
      lNtimesReconstructedK0s   = 0; lNtimesReconstructedLambda   = 0; lNtimesReconstructedAntiLambda   = 0;
      lNtimesReconstructedK0sMI = 0; lNtimesReconstructedLambdaMI = 0; lNtimesReconstructedAntiLambdaMI = 0;
      
      for (Int_t jV0 = 0; jV0 < lEvent->GetNumberOfV0s(); jV0++) {
	
	lLabelTrackN  = 0; lLabelTrackP = 0;
	lPartNMother  = 0; lPartPMother = 0;

	AliAODv0    *vertexAOD= ((AliAODEvent*)lEvent)->GetV0(jV0);
	if (!vertexAOD) continue;
	printf("enter!!");
	AliVParticle  *trackP  = ((AliVEvent*)lEvent)->GetTrack(vertexAOD->GetPosID());
	if (!trackP) continue;
	lLabelTrackP = TMath::Abs(trackP->GetLabel());
	if (lLabelTrackP!=id0 && lLabelTrackP!=id1) continue;
       
	AliVParticle  *trackN  = ((AliVEvent*)lEvent)->GetTrack(vertexAOD->GetNegID());
	if (!trackN) continue;
	lLabelTrackN = TMath::Abs(trackN->GetLabel());
	if (lLabelTrackN!=id0 && lLabelTrackN!=id1) continue;
	
	AliAODMCParticle *lPartNAOD = (AliAODMCParticle*)mcArray->At(lLabelTrackN);
	if (!lPartNAOD) continue;
	AliAODMCParticle *lPartPAOD = (AliAODMCParticle*)mcArray->At(lLabelTrackP);
	if (!lPartPAOD) continue;
	
	lPartNMother = lPartNAOD->GetMother();
	lPartPMother = lPartPAOD->GetMother();

	lPartPMotherPDGcode = ((AliAODMCParticle*)mcArray->At(lPartPMother))->GetPdgCode();
	
	switch (vertexAOD->GetOnFlyStatus()){
	  
	case 0 : 
	  if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==310) ) lNtimesReconstructedK0s++;
	  else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==3122) ) lNtimesReconstructedLambda++;
	  else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==-3122) ) lNtimesReconstructedAntiLambda++;
	  break;
	  
	case 1 :
	  if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==310) ) lNtimesReconstructedK0sMI++;
	  else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==3122) ) lNtimesReconstructedLambdaMI++;
	  else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==-3122) ) lNtimesReconstructedAntiLambdaMI++;
	  break;
	  
	}	
      } // end loop over reconstructed V0s inside MC loop


      if (lPdgcodeCurrentPart==310) {
	fHistMCProdRadiusK0s->Fill(mcPosR);
	fHistMCPtVsYK0s->Fill(lPtCurrentPart,lRapCurrentPart);
	if (TMath::Abs(lRapCurrentPart) > 1) continue;
	fHistNTimesRecK0s->Fill(lNtimesReconstructedK0s);
	fHistNTimesRecK0sMI->Fill(lNtimesReconstructedK0s);
	fHistNTimesRecK0sVsPt->Fill(lPtCurrentPart,lNtimesReconstructedK0s);
	fHistNTimesRecK0sVsPtMI->Fill(lPtCurrentPart,lNtimesReconstructedK0sMI);
      }
      else if (lPdgcodeCurrentPart==3122) {
	fHistMCProdRadiusLambda->Fill(mcPosR);
	fHistMCPtVsYLambda->Fill(lPtCurrentPart,lRapCurrentPart);
	if (TMath::Abs(lRapCurrentPart) > 1) continue;
	fHistNTimesRecLambda->Fill(lNtimesReconstructedLambda);
	fHistNTimesRecLambdaMI->Fill(lNtimesReconstructedLambdaMI);
	fHistNTimesRecLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedLambda);
	fHistNTimesRecLambdaVsPtMI->Fill(lPtCurrentPart,lNtimesReconstructedLambdaMI);
	if (lComeFromSigma) fHistMCPtLambdaFromSigma->Fill(lPtCurrentPart);
      }
      else if (lPdgcodeCurrentPart==-3122) {
	fHistMCProdRadiusAntiLambda->Fill(mcPosR);
	fHistMCPtVsYAntiLambda->Fill(lPtCurrentPart,lRapCurrentPart);
	if (TMath::Abs(lRapCurrentPart) > 1) continue;
	fHistNTimesRecAntiLambda->Fill(lNtimesReconstructedAntiLambda);
	fHistNTimesRecAntiLambdaMI->Fill(lNtimesReconstructedAntiLambdaMI);
	fHistNTimesRecAntiLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedAntiLambda);
	fHistNTimesRecAntiLambdaVsPtMI->Fill(lPtCurrentPart,lNtimesReconstructedAntiLambdaMI);
	if (lComeFromSigma) fHistMCPtAntiLambdaFromSigma->Fill(lPtCurrentPart);
      }

    } // end loop over AODMC particles 
    fHistMCMultiplicityPrimary->Fill(lNbMCPrimary);

  } // end AOD condition

  


  //************************************
  //************************************
  // ESD or AOD loop 
  //************************************

  Double_t lMagneticField = 999;

  //Multiplcity:
  Int_t    nv0sTot= 0, nv0s = 0, nv0sMI = 0;
  
  // Variables:
  
  Double_t  lMomPos[3];
  Double_t  lMomNeg[3];
  Double_t  lV0Position[3];
 
  Double_t lDcaPosToPrimVertex = 0;
  Double_t lDcaNegToPrimVertex = 0;
  Double_t lDcaV0Daughters     = 0;
  Double_t lV0cosPointAngle    = 0;
  Double_t lChi2V0             = 0;
  Double_t lV0DecayLength      = 0;
  Double_t lV0Radius           = 0;
  Double_t lDcaV0ToPrimVertex  = 0;
  
  Int_t    lOnFlyStatus        = 0;
  //Float_t   tdcaPosToPrimVertexXYZ[2], tdcaNegToPrimVertexXYZ[2]; // ..[0] = Impact parameter in XY plane and ..[1] = Impact parameter in Z            
  //Double_t  tdcaDaughterToPrimVertex[2];                          // ..[0] = Pos and ..[1] = Neg

  

  Double_t lInvMassK0s   = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
  Double_t lPt           = 0;
  Double_t lRapK0s       = 0, lRapLambda     = 0;
  Double_t lAlphaV0      = 0, lPtArmV0       = 0;

 
  Float_t lTPCsignalPos = 0;
  Float_t lTPCsignalNeg = 0;

  Double_t lMomentumTrackInPos =0;
  Double_t lMomentumTrackInNeg =0;
  
  // to study Associated V0s:
  Int_t    lIndexTrackPos       = 0, lIndexTrackNeg         = 0;
  UInt_t   lLabelTrackPos       = 0, lLabelTrackNeg         = 0;
  Int_t    lCheckPIdK0Short     = 0, lCheckMcK0Short        = 0;
  Int_t    lCheckPIdLambda      = 0, lCheckMcLambda         = 0;
  Int_t    lCheckPIdAntiLambda  = 0, lCheckMcAntiLambda     = 0;
  Int_t    lCheckSecondaryK0s   = 0, lCheckSecondaryLambda  = 0, lCheckSecondaryAntiLambda  = 0;
  Int_t    lCheckGamma          = 0;
  Double_t mcPosMotherX         = 0, mcPosMotherY           = 0, mcPosMotherZ  = 0;
  Double_t mcPosMotherR         = 0;
  Double_t mcMotherPt           = 0;

  Int_t lIndexPosMother        = 0;
  Int_t lIndexNegMother        = 0;
  Int_t lIndexMotherOfMother   = 0;
  Int_t lPDGCodePosDaughter    = 0;
  Int_t lPDGCodeNegDaughter    = 0;
  Int_t lPdgcodeMother         = 0;
  Int_t lPdgcodeMotherOfMother = 0;

 

  
  Double_t rcPosX        = 0,  rcPosY  = 0, rcPosZ  = 0;
  Double_t rcPosR        = 0;
  
  Double_t deltaPt       = 0;
  
  AliAODTrack  *myPosAodTrack  = new AliAODTrack();
  AliAODTrack  *myNegAodTrack  = new AliAODTrack();

  AliAODv0     *myAODv0        = 0;


  Double_t  lPrimaryVtxPosition[3];
  Double_t  lPrimaryVtxCov[6];
  Double_t  lPrimaryVtxChi2 = 999;

  // Bo: please use external methods for these TF1 (no TFormula like this)...

  // PID - Check Parameters depending on the AliRoot version ! Method to be improved !
  TF1 foPion("foPion", "[0]*([1]*TMath::Power(TMath::Sqrt(1 + (x/0.13957)*(x/0.13957))/(x/0.13957) , [3]) - 1 - TMath::Power(TMath::Sqrt(1 + (x/0.13957)*(x/0.13957))/(x/0.13957) , [3])*TMath::Log([2] + 1/TMath::Power((x/0.13957), [4])))",0.05,20);

  // paramters extracted from calibration....
  foPion.SetParameters(4.23232575531564326e+00,8.68482806165147636e+00,1.34000000000000005e-05,2.30445734159456084e+00,2.25624744086878559e+00);

  TF1 foProton("foProton", "[0]*([1]*TMath::Power(TMath::Sqrt(1 + (x/0.93827)*(x/0.93827))/(x/0.93827) , [3]) - 1 - TMath::Power(TMath::Sqrt(1 + (x/0.93827)*(x/0.93827))/(x/0.93827) , [3])*TMath::Log([2] + 1/TMath::Power((x/0.93827), [4])))",0.05,20);

  // paramters extracted from calibration....
  foProton.SetParameters(4.23232575531564326e+00,8.68482806165147636e+00,1.34000000000000005e-05,2.30445734159456084e+00,2.25624744086878559e+00);



  Float_t nSigmaPosPion   = 0;
  Float_t nSigmaNegPion   = 0;

  Float_t nSigmaPosProton = 0;
  Float_t nSigmaNegProton = 0;

  Float_t cutNSigma = 4.0;


  // ***********************
  // Primary Vertex 
  //************************
     
  AliAODVertex *myPrimaryVertex = 0;

  if(fAnalysisType == "ESD") {   

    const AliESDVertex *myBestPrimaryVertex = ((AliESDEvent*)lEvent)->GetPrimaryVertex();
    if (!myBestPrimaryVertex) return;
    myBestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);
    //lPrimaryVtxPosition[0] = myBestPrimaryVertex->GetXv();
    //lPrimaryVtxPosition[1] = myBestPrimaryVertex->GetYv();
    //lPrimaryVtxPosition[2] = myBestPrimaryVertex->GetZv();
    myBestPrimaryVertex->GetCovMatrix(lPrimaryVtxCov);
    lPrimaryVtxChi2 = myBestPrimaryVertex->GetChi2toNDF();
    
    myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);
    if (!myPrimaryVertex) return;
  }
  
  else if(fAnalysisType == "AOD") {
    printf("enter 1!!");
    myPrimaryVertex = ((AliAODEvent*)lEvent)->GetPrimaryVertex();
    if (!myPrimaryVertex) return;
    lPrimaryVtxPosition[0] = myPrimaryVertex->GetX();
    lPrimaryVtxPosition[1] = myPrimaryVertex->GetY();
    lPrimaryVtxPosition[2] = myPrimaryVertex->GetZ();
    
  }
  fHistPrimaryVertexX->Fill(lPrimaryVtxPosition[0]);
  fHistPrimaryVertexY->Fill(lPrimaryVtxPosition[1]);
  fHistPrimaryVertexZ->Fill(lPrimaryVtxPosition[2]);

  
  //*************************
  // V0 loop
      
  nv0sTot = lEvent->GetNumberOfV0s();

  for (Int_t iV0 = 0; iV0 < nv0sTot; iV0++) {


    lIndexPosMother     = 0; lIndexNegMother     = 0; lIndexMotherOfMother       = 0;
    lCheckPIdK0Short    = 0; lCheckMcK0Short     = 0; lCheckSecondaryK0s         = 0;
    lCheckPIdLambda     = 0; lCheckMcLambda      = 0; lCheckSecondaryLambda      = 0;
    lCheckPIdAntiLambda = 0; lCheckMcAntiLambda  = 0; lCheckSecondaryAntiLambda  = 0;       
    lComeFromSigma      = -1;
    
    
    if(fAnalysisType == "ESD") {
     
      lMagneticField = ((AliESDEvent*)lEvent)->GetMagneticField();

      AliESDv0 *v0 = ((AliESDEvent*)lEvent)->GetV0(iV0);
      if (!v0) {
	continue;
	cout<<"no access to V0 !!"<<endl;
      }
      
      // V0's Daughters
      lIndexTrackPos = TMath::Abs(v0->GetPindex());
      lIndexTrackNeg = TMath::Abs(v0->GetNindex());
      AliESDtrack *myTrackPos = ((AliESDEvent*)lEvent)->GetTrack(lIndexTrackPos);
      AliESDtrack *myTrackNeg = ((AliESDEvent*)lEvent)->GetTrack(lIndexTrackNeg);
      if (!myTrackPos || !myTrackNeg) {
	Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
	continue;
      }
      // Remove like-sign
      if ( myTrackPos->GetSign() == myTrackNeg->GetSign()){
	continue;
      } 

      lLabelTrackPos = (UInt_t)TMath::Abs(myTrackPos->GetLabel());
      lLabelTrackNeg = (UInt_t)TMath::Abs(myTrackNeg->GetLabel());


      // Tracks quality cuts 
      if ( ( (myTrackPos->GetTPCNcls()) < 80 ) || ( (myTrackNeg->GetTPCNcls()) < 80 ) ) continue;

      // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
      if( !(myTrackPos->GetStatus() & AliESDtrack::kTPCrefit)) continue;      
      if( !(myTrackNeg->GetStatus() & AliESDtrack::kTPCrefit)) continue;

      // DCA between daughter and Primary Vertex:
      lDcaPosToPrimVertex = TMath::Abs(myTrackPos->GetD(lPrimaryVtxPosition[0],
						      lPrimaryVtxPosition[1],
						      lMagneticField) );

      lDcaNegToPrimVertex = TMath::Abs(myTrackNeg->GetD(lPrimaryVtxPosition[0],
						      lPrimaryVtxPosition[1],
						      lMagneticField) );
      //if (myTrackPos) myTrackPos->GetImpactParameters(tdcaPosToPrimVertexXYZ[0],tdcaPosToPrimVertexXYZ[1]);
      //if (myTrackNeg) myTrackNeg->GetImpactParameters(tdcaNegToPrimVertexXYZ[0],tdcaNegToPrimVertexXYZ[1]);
      //tdcaDaughterToPrimVertex[0] = TMath::Sqrt(tdcaPosToPrimVertexXYZ[0]*tdcaPosToPrimVertexXYZ[0]+tdcaPosToPrimVertexXYZ[1]*tdcaPosToPrimVertexXYZ[1]);
      //tdcaDaughterToPrimVertex[1] = TMath::Sqrt(tdcaNegToPrimVertexXYZ[0]*tdcaNegToPrimVertexXYZ[0]+tdcaNegToPrimVertexXYZ[1]*tdcaNegToPrimVertexXYZ[1]);
      
      
      // VO's main characteristics:
      lOnFlyStatus       = v0->GetOnFlyStatus();
      lChi2V0            = v0->GetChi2V0();
      lDcaV0Daughters    = v0->GetDcaV0Daughters();
      lDcaV0ToPrimVertex = v0->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
      lV0cosPointAngle   = v0->GetV0CosineOfPointingAngle(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1], lPrimaryVtxPosition[2]);

      v0->GetXYZ(lV0Position[0], lV0Position[1], lV0Position[2]);

      lV0Radius      = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
      lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
		                   TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2) +
		                   TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2],2 ));

      rcPosX   = lV0Position[0];
      rcPosY   = lV0Position[1];
      rcPosZ   = lV0Position[2];
      rcPosR   = lV0Radius;

      //Multiplicity:
      if(!lOnFlyStatus) nv0s++;
      else  if(lOnFlyStatus) nv0sMI++;

      // Invariant mass
      v0->ChangeMassHypothesis(310);
      lInvMassK0s = v0->GetEffMass();
      v0->ChangeMassHypothesis(3122);
      lInvMassLambda = v0->GetEffMass();
      v0->ChangeMassHypothesis(-3122);
      lInvMassAntiLambda = v0->GetEffMass();

      // Rapidity:
      lRapK0s    = v0->Y(310);
      lRapLambda = v0->Y(3122);

      // Daughter momentum:
      v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
      v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);

      // Pt:
      lPt = v0->Pt();

      // Armenteros variables:
      lAlphaV0      = 0; 
      lPtArmV0      = 0;
      
      // PID
      lTPCsignalPos = myTrackPos->GetTPCsignal();
      lTPCsignalNeg = myTrackNeg->GetTPCsignal();

      const AliExternalTrackParam *posInnerParam = myTrackPos->GetInnerParam();
      if (!posInnerParam) continue;
      
      const AliExternalTrackParam *negInnerParam = myTrackNeg->GetInnerParam();
      if (!negInnerParam) continue;
    
      AliExternalTrackParam trackInPos(*posInnerParam);
      AliExternalTrackParam trackInNeg(*negInnerParam);

      //AliExternalTrackParam trackInPos(*myTrackPos->GetInnerParam());
      //AliExternalTrackParam trackInNeg(*myTrackNeg->GetInnerParam());


      lMomentumTrackInPos = trackInPos.GetP(); // momentum for dEdx determination
      lMomentumTrackInNeg = trackInNeg.GetP(); 
      
      if (fUsePID.Contains("withPID")) {
	nSigmaPosPion = TMath::Abs((lTPCsignalPos - foPion.Eval(lMomentumTrackInPos))/foPion.Eval(lMomentumTrackInPos))/0.06;
	
	nSigmaNegPion = TMath::Abs((lTPCsignalNeg - foPion.Eval(lMomentumTrackInNeg))/foPion.Eval(lMomentumTrackInNeg))/0.06; // 0.06 comes from momentum resolution !
	
	nSigmaPosProton = TMath::Abs((lTPCsignalPos - foProton.Eval(lMomentumTrackInPos))/foProton.Eval(lMomentumTrackInPos))/0.06;
	
	nSigmaNegProton = TMath::Abs((lTPCsignalNeg - foProton.Eval(lMomentumTrackInNeg))/foProton.Eval(lMomentumTrackInNeg))/0.06; // 0.06 comes from momentum resolution !
      }
      else {
	nSigmaPosPion = 0; nSigmaNegPion =0; nSigmaPosProton = 0; nSigmaNegProton= 0;}
      
      
      
      // Monte-Carlo particle associated to reconstructed particles:      
      //if (lLabelTrackPos < 0 || lLabelTrackNeg < 0) continue;
      TParticle  *lMCESDPartPos  = stack->Particle(lLabelTrackPos);
      if(!lMCESDPartPos) { 
	Printf("no MC particle for positive and/or negative daughter\n");
	continue;
      }
      TParticle  *lMCESDPartNeg  = stack->Particle(lLabelTrackNeg);
      if (!lMCESDPartNeg) continue;
      lPDGCodePosDaughter = lMCESDPartPos->GetPdgCode();
      lPDGCodeNegDaughter = lMCESDPartNeg->GetPdgCode();
      lIndexPosMother = lMCESDPartPos->GetFirstMother();
      lIndexNegMother = lMCESDPartNeg->GetFirstMother();
      TParticle  *lMCESDMother    = stack->Particle(lIndexPosMother);
      if (!lMCESDMother) continue;
      lPdgcodeMother         = lMCESDMother->GetPdgCode();
      lIndexMotherOfMother   = lMCESDMother->GetFirstMother();
      if (lIndexMotherOfMother ==-1) lPdgcodeMotherOfMother = 0;
      else {
	TParticle  *lMCESDMotherOfMother    = stack->Particle(lIndexMotherOfMother);
	if (!lMCESDMotherOfMother) continue;
	lPdgcodeMotherOfMother = lMCESDMotherOfMother->GetPdgCode();
      }
      
      mcPosX = lMCESDPartPos->Vx();
      mcPosY = lMCESDPartPos->Vy();
      mcPosZ = lMCESDPartPos->Vz();
      mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
      mcPosMotherX = lMCESDMother->Vx();
      mcPosMotherY = lMCESDMother->Vy();
      mcPosMotherZ = lMCESDMother->Vz();
      mcPosMotherR = TMath::Sqrt(mcPosMotherX*mcPosMotherX+mcPosMotherY*mcPosMotherY);
      
      mcMotherPt   = lMCESDMother->Pt();

    } // end ESD condition


    
    else if(fAnalysisType == "AOD") { 

      myAODv0= ((AliAODEvent*)lEvent)->GetV0(iV0);
      if (!myAODv0) {printf("no v0");continue;}

      lDcaPosToPrimVertex = myAODv0->DcaPosToPrimVertex();	
      lDcaNegToPrimVertex = myAODv0->DcaNegToPrimVertex();
      lOnFlyStatus        = myAODv0->GetOnFlyStatus();
      lChi2V0             = myAODv0->Chi2V0();
      lDcaV0Daughters     = myAODv0->DcaV0Daughters();
      lDcaV0ToPrimVertex  = myAODv0->DcaV0ToPrimVertex();
      lV0DecayLength      = myAODv0->DecayLengthV0(lPrimaryVtxPosition);
      lV0cosPointAngle    = myAODv0->CosPointingAngle(lPrimaryVtxPosition);
      lV0Radius           = myAODv0->RadiusV0();
      rcPosX              = myAODv0->DecayVertexV0X();
      rcPosY              = myAODv0->DecayVertexV0Y();
      rcPosZ              = myAODv0->DecayVertexV0Z();

      //Multiplicity:
      if(!lOnFlyStatus) nv0s++;
      else  if(lOnFlyStatus) nv0sMI++;

      lInvMassK0s = myAODv0->MassK0Short();
      lInvMassLambda = myAODv0->MassLambda();
      lInvMassAntiLambda = myAODv0->MassAntiLambda();

      lRapK0s    = myAODv0->RapK0Short();
      lRapLambda = myAODv0->RapLambda();

      lPt          = TMath::Sqrt(myAODv0->Pt2V0());

      lAlphaV0   = myAODv0->AlphaV0();
      lPtArmV0   = myAODv0->PtArmV0();
   
      // V0's Daughters
      lIndexTrackPos = TMath::Abs(myAODv0->GetPosID());
      lIndexTrackNeg = TMath::Abs(myAODv0->GetNegID());
      
      AliVParticle  *lVPartPos  = ((AliVEvent*)lEvent)->GetTrack(lIndexTrackPos);
      AliVParticle  *lVPartNeg  = ((AliVEvent*)lEvent)->GetTrack(lIndexTrackNeg);
      //AliAODTrack  *lVPartPos  = ((AliAODEvent*)lEvent)->GetTrack(lIndexTrackPos);
      //AliAODTrack  *lVPartNeg  = ((AliAODEvent*)lEvent)->GetTrack(lIndexTrackNeg);

      if (!lVPartPos ||(!lVPartNeg )) {
	Printf("strange analysis::UserExec:: Could not retreive one of the daughter track\n");
	continue;
      }

      // Quality cuts:
      // TO DO !!!!!!!

      // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
      //if( !(lVPartPos->GetStatus() & AliAODTrack::kTPCrefit)) continue;      
      //if( !(lVPartNeg->GetStatus() & AliAODTrack::kTPCrefit)) continue;
      
      lLabelTrackPos  = TMath::Abs(lVPartPos->GetLabel());
      lLabelTrackNeg  = TMath::Abs(lVPartNeg->GetLabel());

    
      // Armenteros variables:
      lAlphaV0   = myAODv0->AlphaV0();
      lPtArmV0   = myAODv0->PtArmV0();
      
      // PID not accessible with AOD !
      nSigmaPosPion = 0; nSigmaNegPion =0; nSigmaPosProton = 0; nSigmaNegProton= 0;

      
      // Monte-Carlo particle associated to reconstructed particles:        
      AliAODMCParticle *lMCAODPartPos = (AliAODMCParticle*)mcArray->At(lLabelTrackPos);
      if (!lMCAODPartPos) continue;
      AliAODMCParticle *lMCAODPartNeg = (AliAODMCParticle*)mcArray->At(lLabelTrackNeg);
      if(!lMCAODPartNeg) { 
	Printf("strange analysis::UserExec:no MC particle for negative daughter\n");
	continue;
      }
      lPDGCodePosDaughter = lMCAODPartPos->GetPdgCode();
      lPDGCodeNegDaughter = lMCAODPartNeg->GetPdgCode();
      lIndexPosMother = lMCAODPartPos->GetMother();
      lIndexNegMother = lMCAODPartNeg->GetMother();
      
      AliAODMCParticle *lMCAODMother = (AliAODMCParticle*)mcArray->At(lIndexPosMother);
      lPdgcodeMother = lMCAODMother->GetPdgCode();
      lIndexMotherOfMother  = lMCAODMother->GetMother();
      if (lIndexMotherOfMother ==-1) lPdgcodeMotherOfMother = 0;
      else {
	lPdgcodeMotherOfMother   = ((AliAODMCParticle*)mcArray->At(lIndexMotherOfMother))->GetPdgCode();
      }
      
      mcPosX = lMCAODPartPos->Xv();
      mcPosY = lMCAODPartPos->Yv();
      mcPosZ = lMCAODPartPos->Zv();
      mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
      mcPosMotherX = lMCAODMother->Xv();
      mcPosMotherY = lMCAODMother->Yv();
      mcPosMotherZ = lMCAODMother->Zv();
      mcPosMotherR = TMath::Sqrt(mcPosMotherX*mcPosMotherX+mcPosMotherY*mcPosMotherY);
      mcMotherPt   = lMCAODMother->Pt();
    
            
    } // end AOD condition
    
    
    // Look for associated particles
    if( (lIndexPosMother==-1) || (lIndexNegMother==-1) ) {
      fHistMCDaughterTrack->Fill(1);
    }
    
    else if( ( (lPDGCodePosDaughter==+211) && (lPDGCodeNegDaughter==-211) ) || 
	     ( (lPDGCodePosDaughter==-211) && (lPDGCodeNegDaughter==+211) )    
	     ) {
      lCheckPIdK0Short    = 1;
      fHistMCDaughterTrack->Fill(3);
      if ( (lIndexPosMother==lIndexNegMother) &&
	   (lPdgcodeMother==310) ) {
	if (mcPosMotherR <= lMaxMcProdRadiusPrimaries) lCheckMcK0Short  = 1;
	else lCheckSecondaryK0s = 1;
      }
    }
    else if( ( (lPDGCodePosDaughter==+2212) && (lPDGCodeNegDaughter==-211)  ) ||
	     ( (lPDGCodePosDaughter==-211)  && (lPDGCodeNegDaughter==+2212) )   
	     ) {
      lCheckPIdLambda     = 1;
      fHistMCDaughterTrack->Fill(5);
      if ( (lIndexPosMother==lIndexNegMother) &&
	   (lPdgcodeMother==3122)  ){
	if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3222) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3112) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)
	     ) lComeFromSigma = 1;
	else lComeFromSigma = 0;  
	if (mcPosMotherR <= lMaxMcProdRadiusPrimaries) lCheckMcLambda  = 1; 
	else lCheckSecondaryLambda    = 1;
      }
    }
    else if( ( (lPDGCodePosDaughter==211)   && (lPDGCodeNegDaughter==-2212) ) ||
	     ( (lPDGCodePosDaughter==-2212) && (lPDGCodeNegDaughter==211)   )	     
	     ) {
      lCheckPIdAntiLambda = 1;
      fHistMCDaughterTrack->Fill(7);
      if ( (lIndexPosMother==lIndexNegMother) &&
	   (lPdgcodeMother==-3122) ) {
	if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3222) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3112) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
	     ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)
	     ) lComeFromSigma = 1;
	else lComeFromSigma = 0;  
	if (mcPosMotherR <= lMaxMcProdRadiusPrimaries) lCheckMcAntiLambda  = 1;
	else lCheckSecondaryAntiLambda = 1;
      }
    }
    
    // Gamma conversion
    else if ( (lPDGCodePosDaughter==11) &&
	      (lPDGCodeNegDaughter==-11) &&
	      (lPdgcodeMother==22 ) )
      lCheckGamma = 1;
    

    
    
    // Cuts:
    if (fUseCut.Contains("yes")) {
      if ( (lDcaPosToPrimVertex < 0.036 ) ||
	   (lDcaNegToPrimVertex < 0.036 ) ||
	   (lDcaV0Daughters     > 0.5   ) ||
	   (lV0cosPointAngle    < 0.99 ) 
	   )	
	continue;
    }

    // Pt Resolution:
    deltaPt  = lPt - mcMotherPt;


    // filling histograms
    //cout<<"dca="<<myAODv0->DcaPosToPrimVertex()<<endl;
    fHistDcaPosToPrimVertex->Fill(lDcaPosToPrimVertex,lOnFlyStatus);
    fHistDcaNegToPrimVertex->Fill(lDcaNegToPrimVertex,lOnFlyStatus);
    fHistDcaPosToPrimVertexZoom->Fill(lDcaPosToPrimVertex,lOnFlyStatus);
    fHistDcaNegToPrimVertexZoom->Fill(lDcaNegToPrimVertex,lOnFlyStatus);
    fHistRadiusV0->Fill(lV0Radius,lOnFlyStatus);
    fHistDecayLengthV0->Fill(lV0DecayLength,lOnFlyStatus);
    fHistDcaV0Daughters->Fill(lDcaV0Daughters,lOnFlyStatus);
    fHistChi2->Fill(lChi2V0,lOnFlyStatus);
    fHistCosPointAngle->Fill(lV0cosPointAngle,lOnFlyStatus);
    if (lV0cosPointAngle >= 0.9) fHistCosPointAngleZoom->Fill(lV0cosPointAngle,lOnFlyStatus);
    if (!lOnFlyStatus) {
      fHistProdRadius->Fill(rcPosX,rcPosY);
      fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0);   
    }
    else {
      fHistProdRadiusMI->Fill(rcPosX,rcPosY);
      fHistArmenterosPodolanskiMI->Fill(lAlphaV0,lPtArmV0);    
    }
    
    // K0s associated histograms in |y| < 1:
    if ( (nSigmaPosPion < cutNSigma) && (nSigmaNegPion < cutNSigma) ) {
      if (!lOnFlyStatus) fHistPtVsYK0s->Fill(lPt,lRapK0s);
      else fHistPtVsYK0sMI->Fill(lPt,lRapK0s);

      if (TMath::Abs(lRapK0s) < 1) {

	fHistNsigmaPosPionK0->Fill(nSigmaPosPion);
	fHistNsigmaNegPionK0->Fill(nSigmaNegPion);
	
	switch (lOnFlyStatus){
	case 0 : 
	  fHistMassK0->Fill(lInvMassK0s);
	  fHistMassVsRadiusK0->Fill(lV0Radius,lInvMassK0s);
	  fHistPtVsMassK0->Fill(lInvMassK0s,lPt);
	  if(lCheckPIdK0Short) fHistPidMcMassK0->Fill(lInvMassK0s);
	  if(lCheckMcK0Short) {
	    fHistAsMcMassK0->Fill(lInvMassK0s);
	    fHistAsMcPtK0->Fill(lPt);
	    fHistAsMcPtVsMassK0->Fill(lInvMassK0s,lPt);
	    if (lPt <= 1) fHistAsMcPtZoomK0->Fill(lPt);
	    fHistAsMcMassVsRadiusK0->Fill(lV0Radius,lInvMassK0s);
	    fHistAsMcResxK0->Fill(rcPosX-mcPosX);
	    fHistAsMcResyK0->Fill(rcPosY-mcPosY);
	    fHistAsMcReszK0->Fill(rcPosZ-mcPosZ);
	    fHistAsMcResrVsRadiusK0->Fill(rcPosR,rcPosR-mcPosR);
	    fHistAsMcReszVsRadiusK0->Fill(rcPosR,rcPosZ-mcPosZ);
	    fHistAsMcProdRadiusK0->Fill(mcPosMotherR);
	    fHistAsMcProdRadiusXvsYK0s->Fill(mcPosMotherX,mcPosMotherY);
	    fHistAsMcResPtK0->Fill(deltaPt);
	    fHistAsMcResPtVsRapK0->Fill(deltaPt,lRapK0s);
	    fHistAsMcResPtVsPtK0->Fill(deltaPt,lPt);
	    // Here !
	  }
	  if (lCheckSecondaryK0s) {
	    fHistAsMcSecondaryPtVsYK0s->Fill(lPt,lRapK0s);
	    fHistAsMcSecondaryProdRadiusK0s->Fill(mcPosMotherR);
	    fHistAsMcSecondaryProdRadiusXvsYK0s->Fill(mcPosMotherX,mcPosMotherY);
	    switch (lPdgcodeMotherOfMother) {
	    case 130   : fHistAsMcSecondaryMotherPdgCodeK0s->Fill(0.5);break; // K0L
	    case 321   : fHistAsMcSecondaryMotherPdgCodeK0s->Fill(1.5);break; // K+
	    case -321  : fHistAsMcSecondaryMotherPdgCodeK0s->Fill(2.5);break; // K-
	    case -3122 : fHistAsMcSecondaryMotherPdgCodeK0s->Fill(3.5);break; //AntiLambda
	    default    : fHistAsMcSecondaryMotherPdgCodeK0s->Fill(6.5);break;
	    }
	  }
	  break;
	  
	case 1 :
	  fHistMassK0MI->Fill(lInvMassK0s);
	  fHistMassVsRadiusK0MI->Fill(lV0Radius,lInvMassK0s);
	  fHistPtVsMassK0MI->Fill(lInvMassK0s,lPt);
	  if(lCheckPIdK0Short) fHistPidMcMassK0MI->Fill(lInvMassK0s);
	  if(lCheckMcK0Short) {
	    fHistAsMcMassK0MI->Fill(lInvMassK0s);
	    fHistAsMcPtK0MI->Fill(lPt);
	    fHistAsMcPtVsMassK0MI->Fill(lInvMassK0s,lPt);
	    if (lPt <= 1) fHistAsMcPtZoomK0MI->Fill(lPt);
	    fHistAsMcMassVsRadiusK0MI->Fill(lV0Radius,lInvMassK0s);
	    fHistAsMcResxK0MI->Fill(rcPosX-mcPosX);
	    fHistAsMcResyK0MI->Fill(rcPosY-mcPosY);
	    fHistAsMcReszK0MI->Fill(rcPosZ-mcPosZ);
	    fHistAsMcResrVsRadiusK0MI->Fill(rcPosR,rcPosR-mcPosR);
	    fHistAsMcReszVsRadiusK0MI->Fill(rcPosR,rcPosZ-mcPosZ);
	    fHistAsMcProdRadiusK0MI->Fill(mcPosMotherR);
	    fHistAsMcProdRadiusXvsYK0sMI->Fill(mcPosMotherX,mcPosMotherY);
	    fHistAsMcResPtK0MI->Fill(deltaPt);
	    fHistAsMcResPtVsRapK0MI->Fill(deltaPt,lRapK0s);
	    fHistAsMcResPtVsPtK0MI->Fill(deltaPt,lPt);
	  }
	  else if (lCheckSecondaryK0s) {
	    fHistAsMcSecondaryPtVsYK0sMI->Fill(lPt,lRapK0s);
	    fHistAsMcSecondaryProdRadiusK0sMI->Fill(mcPosMotherR); 
	    fHistAsMcSecondaryProdRadiusXvsYK0sMI->Fill(mcPosMotherX,mcPosMotherY);
	    switch (lPdgcodeMotherOfMother) {
	    case 130   : fHistAsMcSecondaryMotherPdgCodeK0sMI->Fill(0.5);break; // K0L
	    case 321   : fHistAsMcSecondaryMotherPdgCodeK0sMI->Fill(1.5);break; // K+
	    case -321  : fHistAsMcSecondaryMotherPdgCodeK0sMI->Fill(2.5);break; // K-
	    case -3122 : fHistAsMcSecondaryMotherPdgCodeK0sMI->Fill(3.5);break; //AntiLambda
	    default    : fHistAsMcSecondaryMotherPdgCodeK0sMI->Fill(6.5);break;
	    }
	  }
	  break;	
	}
      } // end rapidity condition
    } // end nsigma condition
    

    // Associated Lambda histograms in |y| < 1
    if ( ( (nSigmaPosProton < cutNSigma) && (nSigmaNegPion < cutNSigma) ) || 
	 ( (nSigmaPosPion < cutNSigma)   && (nSigmaNegProton < cutNSigma) ) 
	 ) {

      if (!lOnFlyStatus) fHistPtVsYLambda->Fill(lPt,lRapLambda);
      else fHistPtVsYLambdaMI->Fill(lPt,lRapLambda);

      if (TMath::Abs(lRapLambda) < 1) {

	fHistNsigmaPosProtonLambda->Fill(nSigmaPosProton);
	fHistNsigmaNegPionLambda->Fill(nSigmaNegPion);
	switch (lOnFlyStatus){
	case 0 : 
	  fHistMassLambda->Fill(lInvMassLambda);
	  fHistMassVsRadiusLambda->Fill(lV0Radius,lInvMassLambda);
	  fHistPtVsMassLambda->Fill(lInvMassLambda,lPt);
	  if(lCheckPIdLambda) fHistPidMcMassLambda->Fill(lInvMassLambda);
	  
	  if(lCheckMcLambda) {
	    fHistAsMcMassLambda->Fill(lInvMassLambda);
	    fHistAsMcPtLambda->Fill(lPt);
	    fHistAsMcPtVsMassLambda->Fill(lInvMassLambda,lPt);
	    if (lPt <= 1) fHistAsMcPtZoomLambda->Fill(lPt);
	    fHistAsMcMassVsRadiusLambda->Fill(lV0Radius,lInvMassLambda);
	    fHistAsMcResxLambda->Fill(rcPosX-mcPosX);
	    fHistAsMcResyLambda->Fill(rcPosY-mcPosY);
	    fHistAsMcReszLambda->Fill(rcPosZ-mcPosZ);
	    fHistAsMcResrVsRadiusLambda->Fill(rcPosR,rcPosR-mcPosR);
	    fHistAsMcReszVsRadiusLambda->Fill(rcPosR,rcPosZ-mcPosZ);
	    fHistAsMcProdRadiusLambda->Fill(mcPosMotherR);
	    fHistAsMcProdRadiusXvsYLambda->Fill(mcPosMotherX,mcPosMotherY);
	    fHistAsMcResPtLambda->Fill(deltaPt);
	    fHistAsMcResPtVsRapLambda->Fill(deltaPt,lRapLambda);
	    fHistAsMcResPtVsPtLambda->Fill(deltaPt,lPt);
	    if (lComeFromSigma) fHistAsMcPtLambdaFromSigma->Fill(lPt);
	    switch (lPdgcodeMotherOfMother) {
	    case 3222 : fHistAsMcMotherPdgCodeLambda->Fill(0.5); break; // Sigma +
	    case 3212 : fHistAsMcMotherPdgCodeLambda->Fill(1.5); break; // Sigma 0
	    case 3112 : fHistAsMcMotherPdgCodeLambda->Fill(2.5); break;// Sigma -
	    case 3224 : fHistAsMcMotherPdgCodeLambda->Fill(3.5); break;// Sigma(1385) +
	    case 3214 : fHistAsMcMotherPdgCodeLambda->Fill(4.5); break;// Sigma(1385) 0
	    case 3114 : fHistAsMcMotherPdgCodeLambda->Fill(5.5); break;// Sigma(1385) -
	    case 3322 : fHistAsMcMotherPdgCodeLambda->Fill(6.5);break; // Xi 0
	    case 3312 : fHistAsMcMotherPdgCodeLambda->Fill(7.5);break; // Xi -
	    case 3334 : fHistAsMcMotherPdgCodeLambda->Fill(8.5);break; // Omega
	    case -1   : fHistAsMcMotherPdgCodeLambda->Fill(9.5);break;
	    default   : fHistAsMcMotherPdgCodeLambda->Fill(10.5);break;
	    }
	  }
	  
	  else if (lCheckSecondaryLambda) {
	    fHistAsMcSecondaryPtVsYLambda->Fill(lPt,lRapLambda);
	    fHistAsMcSecondaryProdRadiusLambda->Fill(mcPosMotherR); 
	    fHistAsMcSecondaryProdRadiusXvsYLambda->Fill(mcPosMotherX,mcPosMotherY);
	    if (lComeFromSigma) fHistAsMcSecondaryPtLambdaFromSigma->Fill(lPt);
	    switch (lPdgcodeMotherOfMother) {
	    case 3222 : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(0.5); break;// Sigma +
	    case 3212 : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(1.5); break;// Sigma 0
	    case 3112 : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(2.5); break;// Sigma -
	    case 3224 : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(3.5); break;// Sigma(1385) +
	    case 3214 : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(4.5); break;// Sigma(1385) 0
	    case 3114 : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(5.5); break;// Sigma(1385) -
	    case 3322 : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(6.5); break; // Xi 0
	    case 3312 : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(7.5); break; // Xi -
	    case 3334 : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(8.5); break; // Omega
	    case -1   : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(9.5); break;
	    default   : fHistAsMcSecondaryMotherPdgCodeLambda->Fill(10.5);break;
	    }
	  }
	  break;
	  
	case 1 :
	  fHistMassLambdaMI->Fill(lInvMassLambda);
	  fHistMassVsRadiusLambdaMI->Fill(lV0Radius,lInvMassLambda);
	  fHistPtVsMassLambdaMI->Fill(lInvMassLambda,lPt);
	  if(lCheckPIdLambda) fHistPidMcMassLambdaMI->Fill(lInvMassLambda);
	  
	  if(lCheckMcLambda) {
	    fHistAsMcMassLambdaMI->Fill(lInvMassLambda);
	    fHistAsMcPtLambdaMI->Fill(lPt);
	    fHistAsMcPtVsMassLambdaMI->Fill(lInvMassLambda,lPt);
	    fHistAsMcMassVsRadiusLambdaMI->Fill(lV0Radius,lInvMassLambda);
	    fHistAsMcResxLambdaMI->Fill(rcPosX-mcPosX);
	    fHistAsMcResyLambdaMI->Fill(rcPosY-mcPosY);
	    fHistAsMcReszLambdaMI->Fill(rcPosZ-mcPosZ);
	    fHistAsMcResrVsRadiusLambdaMI->Fill(rcPosR,rcPosR-mcPosR);
	    fHistAsMcReszVsRadiusLambdaMI->Fill(rcPosR,rcPosZ-mcPosZ);
	    fHistAsMcProdRadiusLambdaMI->Fill(mcPosMotherR);
	    fHistAsMcProdRadiusXvsYLambdaMI->Fill(mcPosMotherX,mcPosMotherY);
	    fHistAsMcResPtLambdaMI->Fill(deltaPt);
	    fHistAsMcResPtVsRapLambdaMI->Fill(deltaPt,lRapLambda);
	    fHistAsMcResPtVsPtLambdaMI->Fill(deltaPt,lPt);
	    if (lComeFromSigma) fHistAsMcPtLambdaFromSigmaMI->Fill(lPt);
	    switch (lPdgcodeMotherOfMother) {
	    case 3222 : fHistAsMcMotherPdgCodeLambdaMI->Fill(0.5); break; // Sigma +
	    case 3212 : fHistAsMcMotherPdgCodeLambdaMI->Fill(1.5); break; // Sigma 0
	    case 3112 : fHistAsMcMotherPdgCodeLambdaMI->Fill(2.5); break;// Sigma -
	    case 3224 : fHistAsMcMotherPdgCodeLambdaMI->Fill(3.5); break;// Sigma(1385) +
	    case 3214 : fHistAsMcMotherPdgCodeLambdaMI->Fill(4.5); break;// Sigma(1385) 0
	    case 3114 : fHistAsMcMotherPdgCodeLambdaMI->Fill(5.5); break;// Sigma(1385) -
	    case 3322 : fHistAsMcMotherPdgCodeLambdaMI->Fill(6.5);break; // Xi 0
	    case 3312 : fHistAsMcMotherPdgCodeLambdaMI->Fill(7.5);break; // Xi -
	    case 3334 : fHistAsMcMotherPdgCodeLambdaMI->Fill(8.5);break; // Omega
	    case -1   : fHistAsMcMotherPdgCodeLambdaMI->Fill(9.5);break;
	    default   : fHistAsMcMotherPdgCodeLambdaMI->Fill(10.5);break;
	    }
	  }
	  else if (lCheckSecondaryLambda) {
	    fHistAsMcSecondaryPtVsYLambdaMI->Fill(lPt,lRapLambda);
	    fHistAsMcSecondaryProdRadiusLambdaMI->Fill(mcPosMotherR); 
	    fHistAsMcSecondaryProdRadiusXvsYLambdaMI->Fill(mcPosMotherX,mcPosMotherY);
	    if (lComeFromSigma) fHistAsMcSecondaryPtLambdaFromSigmaMI->Fill(lPt);
	    switch (lPdgcodeMotherOfMother) {
	    case 3222 : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(0.5); break;// Sigma +
	    case 3212 : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(1.5); break;// Sigma 0
	    case 3112 : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(2.5); break;// Sigma -
	    case 3224 : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(3.5); break;// Sigma(1385) +
	    case 3214 : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(4.5); break;// Sigma(1385) 0
	    case 3114 : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(5.5); break;// Sigma(1385) -
	    case 3322 : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(6.5); break; // Xi 0
	    case 3312 : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(7.5); break; // Xi -
	    case 3334 : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(8.5); break; // Omega
	    case -1   : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(9.5); break;
	    default   : fHistAsMcSecondaryMotherPdgCodeLambdaMI->Fill(10.5);break;
	    }
	  }
	  break;	
	}
      } // end rapidity condition
    } //end nsigma condition - lambda



    // AntiLambda associated histograms in |y| < 1
    if ( ( (nSigmaPosPion < cutNSigma)   && (nSigmaNegProton < cutNSigma) ) ||
	 ( (nSigmaPosProton < cutNSigma) && (nSigmaNegPion < cutNSigma) )
	 ) {

      if (!lOnFlyStatus) fHistPtVsYAntiLambda->Fill(lPt,lRapLambda);
      else fHistPtVsYAntiLambdaMI->Fill(lPt,lRapLambda);
    
      if (TMath::Abs(lRapLambda) < 1) {

	fHistNsigmaPosPionAntiLambda->Fill(nSigmaPosPion);
	fHistNsigmaNegProtonAntiLambda->Fill(nSigmaNegProton);

	switch (lOnFlyStatus){
	case 0 : 
	  fHistMassAntiLambda->Fill(lInvMassAntiLambda);
	  fHistMassVsRadiusAntiLambda->Fill(lV0Radius,lInvMassAntiLambda);
	  fHistPtVsMassAntiLambda->Fill(lInvMassAntiLambda,lPt);
	  if (lCheckPIdAntiLambda) fHistPidMcMassAntiLambda->Fill(lInvMassAntiLambda);
	  
	  if(lCheckMcAntiLambda) {
	    fHistAsMcMassAntiLambda->Fill(lInvMassAntiLambda);
	    fHistAsMcPtAntiLambda->Fill(lPt);
	    fHistAsMcPtVsMassAntiLambda->Fill(lInvMassAntiLambda,lPt);
	    fHistAsMcMassVsRadiusAntiLambda->Fill(lV0Radius,lInvMassAntiLambda);
	    fHistAsMcResxAntiLambda->Fill(rcPosX-mcPosX);
	    fHistAsMcResyAntiLambda->Fill(rcPosY-mcPosY);
	    fHistAsMcReszAntiLambda->Fill(rcPosZ-mcPosZ);
	    fHistAsMcResrVsRadiusAntiLambda->Fill(rcPosR,rcPosR-mcPosR);
	    fHistAsMcReszVsRadiusAntiLambda->Fill(rcPosR,rcPosZ-mcPosZ);
	    fHistAsMcProdRadiusAntiLambda->Fill(mcPosMotherR);
	    fHistAsMcProdRadiusXvsYAntiLambda->Fill(mcPosMotherX,mcPosMotherY);
	    fHistAsMcResPtAntiLambda->Fill(deltaPt);
	    fHistAsMcResPtVsRapAntiLambda->Fill(deltaPt,lRapLambda);
	    fHistAsMcResPtVsPtAntiLambda->Fill(deltaPt,lPt);
	    if (lComeFromSigma) fHistAsMcPtAntiLambdaFromSigma->Fill(lPt);
	    switch (lPdgcodeMotherOfMother) {
	    case -3222 : fHistAsMcMotherPdgCodeAntiLambda->Fill(0.5); break;// Sigma +
	    case -3212 : fHistAsMcMotherPdgCodeAntiLambda->Fill(1.5); break;// Sigma 0
	    case -3112 : fHistAsMcMotherPdgCodeAntiLambda->Fill(2.5); break;// Sigma -
	    case -3224 : fHistAsMcMotherPdgCodeAntiLambda->Fill(3.5); break;// Sigma(1385) +
	    case -3214 : fHistAsMcMotherPdgCodeAntiLambda->Fill(4.5); break;// Sigma(1385) 0
	    case -3114 : fHistAsMcMotherPdgCodeAntiLambda->Fill(5.5); break;// Sigma(1385) -
	    case -3322 : fHistAsMcMotherPdgCodeAntiLambda->Fill(6.5); break; // Xi 0
	    case -3312 : fHistAsMcMotherPdgCodeAntiLambda->Fill(7.5); break; // Xi -
	    case -3334 : fHistAsMcMotherPdgCodeAntiLambda->Fill(8.5); break; // Omega
	    case -1    : fHistAsMcMotherPdgCodeAntiLambda->Fill(9.5); break;
	    default    : fHistAsMcMotherPdgCodeAntiLambda->Fill(10.5);break;
	    }
	  }
	  else if (lCheckSecondaryAntiLambda) {
	    fHistAsMcSecondaryPtVsYAntiLambda->Fill(lPt,lRapLambda);
	    fHistAsMcSecondaryProdRadiusAntiLambda->Fill(mcPosMotherR);
	    fHistAsMcSecondaryProdRadiusXvsYAntiLambda->Fill(mcPosMotherX,mcPosMotherY);
	    if (lComeFromSigma) fHistAsMcSecondaryPtAntiLambdaFromSigma->Fill(lPt);
	    switch (lPdgcodeMotherOfMother) {
	    case -3222 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(0.5); break;// Sigma +
	    case -3212 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(1.5); break;// Sigma 0
	    case -3112 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(2.5); break;// Sigma -
	    case -3224 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(3.5); break;// Sigma(1385) +
	    case -3214 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(4.5); break;// Sigma(1385) 0
	    case -3114 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(5.5); break;// Sigma(1385) -
	    case -3322 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(6.5); break; // Xi 0
	    case -3312 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(7.5); break; // Xi -
	    case -3334 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(8.5); break; // Omega
	    case -1    : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(9.5); break;
	    default    : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(10.5);break;
	    }
	  }
	  break;
	  
	  
	case 1 :
	  fHistMassAntiLambdaMI->Fill(lInvMassAntiLambda);
	  fHistMassVsRadiusAntiLambdaMI->Fill(lV0Radius,lInvMassAntiLambda);
	  fHistPtVsMassAntiLambdaMI->Fill(lInvMassAntiLambda,lPt);
	  if (lCheckPIdAntiLambda) fHistPidMcMassAntiLambdaMI->Fill(lInvMassAntiLambda);
	  
	  if(lCheckMcAntiLambda) {
	    fHistAsMcMassAntiLambdaMI->Fill(lInvMassAntiLambda);
	    fHistAsMcPtAntiLambdaMI->Fill(lPt);
	    fHistAsMcPtVsMassAntiLambdaMI->Fill(lInvMassAntiLambda,lPt);
	    fHistAsMcMassVsRadiusAntiLambdaMI->Fill(lV0Radius,lInvMassAntiLambda);
	    fHistAsMcResxAntiLambdaMI->Fill(rcPosX-mcPosX);
	    fHistAsMcResyAntiLambdaMI->Fill(rcPosY-mcPosY);
	    fHistAsMcReszAntiLambdaMI->Fill(rcPosZ-mcPosZ);
	    fHistAsMcResrVsRadiusAntiLambdaMI->Fill(rcPosR,rcPosR-mcPosR);
	    fHistAsMcReszVsRadiusAntiLambdaMI->Fill(rcPosR,rcPosZ-mcPosZ);
	    fHistAsMcProdRadiusAntiLambdaMI->Fill(mcPosMotherR);
	    fHistAsMcProdRadiusXvsYAntiLambdaMI->Fill(mcPosMotherX,mcPosMotherY);
	    fHistAsMcResPtAntiLambdaMI->Fill(deltaPt);
	    fHistAsMcResPtVsRapAntiLambdaMI->Fill(deltaPt,lRapLambda);
	    fHistAsMcResPtVsPtAntiLambdaMI->Fill(deltaPt,lPt);
	    if (lComeFromSigma) fHistAsMcPtAntiLambdaFromSigmaMI->Fill(lPt);
	    switch (lPdgcodeMotherOfMother) {
	    case -3222 : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(0.5); break;// Sigma +
	    case -3212 : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(1.5); break;// Sigma 0
	    case -3112 : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(2.5); break;// Sigma -
	    case -3224 : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(3.5); break;// Sigma(1385) +
	    case -3214 : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(4.5); break;// Sigma(1385) 0
	    case -3114 : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(5.5); break;// Sigma(1385) -
	    case -3322 : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(6.5); break; // Xi 0
	    case -3312 : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(7.5); break; // Xi -
	    case -3334 : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(8.5); break; // Omega
	    case -1    : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(9.5); break;
	    default    : fHistAsMcMotherPdgCodeAntiLambdaMI->Fill(10.5);break;
	    }
	  }
	  else if (lCheckSecondaryAntiLambda) {
	    fHistAsMcSecondaryPtVsYAntiLambdaMI->Fill(lPt,lRapLambda);
	    fHistAsMcSecondaryProdRadiusAntiLambdaMI->Fill(mcPosMotherR); 
	    fHistAsMcSecondaryProdRadiusXvsYAntiLambdaMI->Fill(mcPosMotherX,mcPosMotherY);
	    if (lComeFromSigma) fHistAsMcSecondaryPtAntiLambdaFromSigmaMI->Fill(lPt);
	    switch (lPdgcodeMotherOfMother) {
	    case -3222 : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(0.5); break;// Sigma +
	    case -3212 : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(1.5); break;// Sigma 0
	    case -3112 : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(2.5); break;// Sigma -
	    case -3224 : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(3.5); break;// Sigma(1385) +
	    case -3214 : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(4.5); break;// Sigma(1385) 0
	    case -3114 : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(5.5); break;// Sigma(1385) -
	    case -3322 : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(6.5); break; // Xi 0
	    case -3312 : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(7.5); break; // Xi -
	    case -3334 : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(8.5); break; // Omega
	    case -1    : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(9.5); break;
	    default    : fHistAsMcSecondaryMotherPdgCodeAntiLambdaMI->Fill(10.5);break;
	    }
	  }
	  break;	
	}
      } // end rapidity condition
    } // end nsigma condition - antilambda

 
  } // end V0 loop

  fHistV0Multiplicity->Fill(nv0s);
  fHistV0MultiplicityMI->Fill(nv0sMI);

  if (fAnalysisType == "ESD") { if(myPrimaryVertex) delete myPrimaryVertex;}
  if (myPosAodTrack) delete myPosAodTrack;
  if (myNegAodTrack) delete myNegAodTrack;



  
  // Post output data
  PostData(1, fListHist);
}      

//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fHistV0Multiplicity = dynamic_cast<TH1F*> (((TList*)GetOutputData(1))->FindObject("fHistV0Multiplicity"));
  if (!fHistV0Multiplicity) {
    Printf("ERROR: fHistV0Multiplicity not available");
    return;
  }

  fHistV0MultiplicityMI = dynamic_cast<TH1F*> (((TList*)GetOutputData(1))->FindObject("fHistV0MultiplicityMI"));
  if (!fHistV0MultiplicityMI) {
    Printf("ERROR: fHistV0MultiplicityMI not available");
    return;
  }

  TCanvas *canPerformanceStrange = new TCanvas("AliAnalysisTaskPerformanceStrange","Multiplicity",10,10,510,510);
  canPerformanceStrange->Divide(2,1);
  if (fHistV0Multiplicity->GetMaximum() > 0.) canPerformanceStrange->cd(1)->SetLogy();
  fHistV0Multiplicity->SetMarkerStyle(25);
  fHistV0Multiplicity->DrawCopy("E");
  if (fHistV0MultiplicityMI->GetMaximum() > 0.) canPerformanceStrange->cd(2)->SetLogy();
  fHistV0MultiplicityMI->SetMarkerStyle(24);
  fHistV0MultiplicityMI->DrawCopy("E");
  


  
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskPerformanceStrange::MyRapidity(Double_t rE, Double_t rPz) const
{
  // Local calculation for rapidity
  return 0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
} 
//----------------------------------------------------------------------------

