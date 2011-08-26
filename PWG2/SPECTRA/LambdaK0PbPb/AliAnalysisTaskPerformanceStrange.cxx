
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//	        AliAnalysisTaskPerformanceStrange class
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
//#include "TH3F.h"
#include "TF1.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"

#include "AliAnalysisManager.h"

#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliMultiplicity.h"

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
#include "AliGenEventHeader.h"

#include "AliLog.h"

#include "AliKFVertex.h"
#include "AliVertexerTracks.h"

#include "AliAnalysisTaskPerformanceStrange.h"
#include "AliAnalysisCentralitySelector.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"



ClassImp(AliAnalysisTaskPerformanceStrange)


//________________________________________________________________________
AliAnalysisTaskPerformanceStrange::AliAnalysisTaskPerformanceStrange()
: AliAnalysisTaskSE(), fAnalysisMC(0), fAnalysisType("infoType"),  fCollidingSystems(0), fUsePID("infoPID"), fUseCut("infoCut"),fDown(0),fUp(0), fESD(0), fListHist(0),fCentrSelector(0),fTracksCuts(0),fPIDResponse(0), 

  fHistMCPrimaryVertexX(0),
  fHistMCPrimaryVertexY(0),
  fHistMCPrimaryVertexZ(0),
  fHistMCMultiplicityPrimary(0),
  fHistMCMultiplicityTracks(0),
  fHistTPCTracks(0),
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
  fHistMCRapK0s(0),
  fHistMCRapInPtRangeK0s(0),
  fHistMCRapLambda(0),
  fHistMCRapInPtRangeLambda(0),
  fHistMCRapAntiLambda(0),
  fHistMCRapInPtRangeAntiLambda(0),
  fHistMCRapXi(0),
  fHistMCRapInPtRangeXi(0),
  fHistMCRapPhi(0),
  fHistMCRapInPtRangePhi(0),
////////////////////////////////////////
  fHistMCPtK0s(0),
  fHistMCPtLambda(0),
  fHistMCPtAntiLambda(0),
///////////////////////////////////////////

  fHistMCPtLambdaFromSigma(0),
  fHistMCPtAntiLambdaFromSigma(0),
  fHistNTimesRecK0s(0),
  fHistNTimesRecLambda(0),
  fHistNTimesRecAntiLambda(0),
  fHistNTimesRecK0sVsPt(0),
  fHistNTimesRecLambdaVsPt(0),
  fHistNTimesRecAntiLambdaVsPt(0),
  fHistNumberEvents(0),
  fHistTrackPerEvent(0),
  fHistTPCMult(0),
  fHistTrackletPerEvent(0),
  fHistMCDaughterTrack(0),
  fHistSPDPrimaryVertexZ(0),
  fHistPrimaryVertexX(0),
  fHistPrimaryVertexY(0),
  fHistPrimaryVertexZ(0),
  fHistPrimaryVertexResX(0),
  fHistPrimaryVertexResY(0),
  fHistPrimaryVertexResZ(0),
  fHistPrimaryVertexPosXV0events(0), 
  fHistPrimaryVertexPosYV0events(0), 
  fHistPrimaryVertexPosZV0events(0),
  fHistDaughterPt(0),
////////////////////////////////////////

  fHistDcaPosToPrimVertexK0(0),
  fHistDcaNegToPrimVertexK0(0),
  fHistRadiusV0K0(0),
  fHistDecayLengthV0K0(0),
  fHistDcaV0DaughtersK0(0),
  fHistChi2K0(0),
  fHistCosPointAngleK0(0),

  fHistDcaPosToPrimVertexK0vsMassK0(0),
  fHistDcaNegToPrimVertexK0vsMassK0(0),
  fHistRadiusV0K0vsMassK0(0),
  fHistDecayLengthV0K0vsMassK0(0),
  fHistDcaV0DaughtersK0vsMassK0(0),
  fHistCosPointAngleK0vsMassK0(0),
 
   fHistDcaPosToPrimVertexK0vsMassK0pt1(0),
   fHistDcaNegToPrimVertexK0vsMassK0pt1(0),
   fHistRadiusV0K0vsMassK0pt1(0),
   fHistDecayLengthV0K0vsMassK0pt1(0),
   fHistDcaV0DaughtersK0vsMassK0pt1(0),
   fHistCosPointAngleK0vsMassK0pt1(0),

   fHistDcaPosToPrimVertexK0vsMassK0pt2(0),
   fHistDcaNegToPrimVertexK0vsMassK0pt2(0),
   fHistRadiusV0K0vsMassK0pt2(0),
   fHistDecayLengthV0K0vsMassK0pt2(0),
   fHistDcaV0DaughtersK0vsMassK0pt2(0),
   fHistCosPointAngleK0vsMassK0pt2(0),

   fHistDcaPosToPrimVertexK0vsMassK0pt3(0),
   fHistDcaNegToPrimVertexK0vsMassK0pt3(0),
   fHistRadiusV0K0vsMassK0pt3(0),
   fHistDecayLengthV0K0vsMassK0pt3(0),
   fHistDcaV0DaughtersK0vsMassK0pt3(0),
   fHistCosPointAngleK0vsMassK0pt3(0),

/////////// Lambda ///////////////////////////

    fHistDcaPosToPrimVertexL(0),
    fHistDcaNegToPrimVertexL(0),
    fHistRadiusV0L(0),
    fHistDecayLengthV0L(0),
    fHistDcaV0DaughtersL(0),
    fHistChi2L(0),
    fHistCosPointAngleL(0),

    fHistDcaPosToPrimVertexLvsMassL(0),
    fHistDcaNegToPrimVertexLvsMassL(0),
    fHistRadiusV0LvsMassL(0),
    fHistDecayLengthV0LvsMassL(0),
    fHistDcaV0DaughtersLvsMassL(0),
    fHistCosPointAngleLvsMassL(0),


    
      fHistDcaPosToPrimVertexLambdaVsMasspt1(0),
      fHistDcaNegToPrimVertexLambdaVsMasspt1(0),
      fHistRadiusV0LambdaVsMasspt1(0),
      fHistDecayLengthV0LambdaVsMasspt1(0),
      fHistDcaV0DaughtersLambdaVsMasspt1(0),
      fHistCosPointAngleLambdaVsMasspt1(0),

      fHistDcaPosToPrimVertexLambdaVsMasspt2(0),
      fHistDcaNegToPrimVertexLambdaVsMasspt2(0),
      fHistRadiusV0LambdaVsMasspt2(0),
      fHistDecayLengthV0LambdaVsMasspt2(0),
      fHistDcaV0DaughtersLambdaVsMasspt2(0),
      fHistCosPointAngleLambdaVsMasspt2(0),

      fHistDcaPosToPrimVertexLambdaVsMasspt3(0),
      fHistDcaNegToPrimVertexLambdaVsMasspt3(0),
      fHistRadiusV0LambdaVsMasspt3(0),
      fHistDecayLengthV0LambdaVsMasspt3(0),
      fHistDcaV0DaughtersLambdaVsMasspt3(0),
      fHistCosPointAngleLambdaVsMasspt3(0),



/////////Antilambda ///////////////////
    fHistDcaPosToPrimVertexAntiL(0),
    fHistDcaNegToPrimVertexAntiL(0),
    fHistRadiusV0AntiL(0),
    fHistDecayLengthV0AntiL(0),
    fHistDcaV0DaughtersAntiL(0),
    fHistChi2AntiL(0),
    fHistCosPointAngleAntiL(0),

    fHistDcaPosToPrimVertexAntiLvsMass(0),
    fHistDcaNegToPrimVertexAntiLvsMass(0),
    fHistRadiusV0AntiLvsMass(0),
    fHistDecayLengthV0AntiLvsMass(0),
    fHistDcaV0DaughtersAntiLvsMass(0),
    fHistCosPointAngleAntiLvsMass(0),

    
    
      fHistDcaPosToPrimVertexAntiLVsMasspt1(0),
      fHistDcaNegToPrimVertexAntiLVsMasspt1(0),
      fHistRadiusV0AntiLVsMasspt1(0),
      fHistDecayLengthV0AntiLVsMasspt1(0),
      fHistDcaV0DaughtersAntiLVsMasspt1(0),
      fHistCosPointAngleAntiLVsMasspt1(0),

      fHistDcaPosToPrimVertexAntiLVsMasspt2(0),
      fHistDcaNegToPrimVertexAntiLVsMasspt2(0),
      fHistRadiusV0AntiLVsMasspt2(0),
      fHistDecayLengthV0AntiLVsMasspt2(0),
      fHistDcaV0DaughtersAntiLVsMasspt2(0),
      fHistCosPointAngleAntiLVsMasspt2(0),

      fHistDcaPosToPrimVertexAntiLVsMasspt3(0),
      fHistDcaNegToPrimVertexAntiLVsMasspt3(0),
      fHistRadiusV0AntiLVsMasspt3(0),
      fHistDecayLengthV0AntiLVsMasspt3(0),
      fHistDcaV0DaughtersAntiLVsMasspt3(0),
      fHistCosPointAngleAntiLVsMasspt3(0),

/////////////////////////////////////////
  fHistV0Multiplicity(0),
  fHistMassK0(0),
  fHistMassLambda(0),
  fHistMassAntiLambda(0),
  fHistMassVsRadiusK0(0),
  fHistMassVsRadiusLambda(0),
  fHistMassVsRadiusAntiLambda(0),

///////////////////////////////////////
  fHistPtVsMassK0(0),
  fHistPtVsMassLambda(0),
 fHistPtVsMassAntiLambda(0),
////////////////////////////////////////

  fHistArmenterosPodolanski(0),
  fHistK0sMassVsLambdaMass(0),
   fHistTPCsignalPt(0),
  fHistNSigmaProton(0),    
  fHistNsigmaPosPionAntiLambda(0),
  fHistNsigmaNegProtonAntiLambda(0),
  fHistNsigmaPosProtonLambda(0),
  fHistNsigmaNegPionLambda(0),
  fHistNsigmaPosProtonAntiLambda(0),
  fHistNsigmaNegPionAntiLambda(0),
  fHistNsigmaPosPionK0(0),
  fHistNsigmaNegPionK0(0),
  fHistAsMcRapK0(0),
  fHistAsMcRapLambda(0),
  fHistAsMcRapAntiLambda(0),
  fHistAsMcPtK0(0),
  fHistAsMcPtLambda(0),
  fHistAsMcPtAntiLambda(0),
  fHistAsMcPtZoomK0(0),
  fHistAsMcPtZoomLambda(0),
  fHistAsMcPtZoomAntiLambda(0),
  fHistAsMcProdRadiusK0(0),
  fHistAsMcProdRadiusLambda(0),
  fHistAsMcProdRadiusAntiLambda(0),
  fHistAsMcProdRadiusXvsYK0s(0),
  fHistAsMcProdRadiusXvsYLambda(0),
  fHistAsMcProdRadiusXvsYAntiLambda(0),
  fHistPidMcMassK0(0),
  fHistPidMcMassLambda(0),
  fHistPidMcMassAntiLambda(0),
  fHistAsMcMassK0(0),
  fHistAsMcMassLambda(0),
  fHistAsMcMassAntiLambda(0),
  fHistAsMcPtVsMassK0(0),
  fHistAsMcPtVsMassLambda(0),
  fHistAsMcPtVsMassAntiLambda(0),
  fHistAsMcMassVsRadiusK0(0),
  fHistAsMcMassVsRadiusLambda(0),
  fHistAsMcMassVsRadiusAntiLambda(0),
  fHistAsMcResxK0(0),
  fHistAsMcResyK0(0),
  fHistAsMcReszK0(0),
  fHistAsMcResrVsRadiusK0(0),
  fHistAsMcReszVsRadiusK0(0),
  fHistAsMcResxLambda(0),
  fHistAsMcResyLambda(0),
  fHistAsMcReszLambda(0),
  fHistAsMcResrVsRadiusLambda(0),
  fHistAsMcReszVsRadiusLambda(0),
  fHistAsMcResxAntiLambda(0),
  fHistAsMcResyAntiLambda(0),
  fHistAsMcReszAntiLambda(0),
  fHistAsMcResrVsRadiusAntiLambda(0),
  fHistAsMcReszVsRadiusAntiLambda(0),
  fHistAsMcResPtK0(0),
  fHistAsMcResPtLambda(0),
  fHistAsMcResPtAntiLambda(0),
  fHistAsMcResPtVsRapK0(0),
  fHistAsMcResPtVsRapLambda(0),
  fHistAsMcResPtVsRapAntiLambda(0),
  fHistAsMcResPtVsPtK0(0),
  fHistAsMcResPtVsPtLambda(0),
  fHistAsMcResPtVsPtAntiLambda(0),
  fHistAsMcMotherPdgCodeK0s(0),
  fHistAsMcMotherPdgCodeLambda(0),
  fHistAsMcMotherPdgCodeAntiLambda(0),
  fHistAsMcPtLambdaFromSigma(0),
  fHistAsMcPtAntiLambdaFromSigma(0),
  fHistAsMcSecondaryPtVsRapK0s(0),
  fHistAsMcSecondaryPtVsRapLambda(0),
  fHistAsMcSecondaryPtVsRapAntiLambda(0),
  fHistAsMcSecondaryProdRadiusK0s(0),
  fHistAsMcSecondaryProdRadiusLambda(0),
  fHistAsMcSecondaryProdRadiusAntiLambda(0),
  fHistAsMcSecondaryProdRadiusXvsYK0s(0),
  fHistAsMcSecondaryProdRadiusXvsYLambda(0),
  fHistAsMcSecondaryProdRadiusXvsYAntiLambda(0),
  fHistAsMcSecondaryMotherPdgCodeK0s(0),
  fHistAsMcSecondaryMotherPdgCodeLambda(0),
  fHistAsMcSecondaryMotherPdgCodeAntiLambda(0),
  fHistAsMcSecondaryPtLambdaFromSigma(0),
  fHistAsMcSecondaryPtAntiLambdaFromSigma(0)
 
{
  // Constructor
}





//________________________________________________________________________
AliAnalysisTaskPerformanceStrange::AliAnalysisTaskPerformanceStrange(const char *name)
  : AliAnalysisTaskSE(name), fAnalysisMC(0), fAnalysisType("infoType"), fCollidingSystems(0), fUsePID("infoPID"), fUseCut("infocut"),fDown(0),fUp(0), fESD(0), fListHist(),fCentrSelector(0), fTracksCuts(0),fPIDResponse(0),

    fHistMCPrimaryVertexX(0),
    fHistMCPrimaryVertexY(0),
    fHistMCPrimaryVertexZ(0),
    fHistMCMultiplicityPrimary(0),
    fHistMCMultiplicityTracks(0),
    fHistTPCTracks(0),
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
    fHistMCRapK0s(0),
    fHistMCRapInPtRangeK0s(0),
    fHistMCRapLambda(0),
    fHistMCRapInPtRangeLambda(0),
    fHistMCRapAntiLambda(0),
    fHistMCRapInPtRangeAntiLambda(0),
    fHistMCRapXi(0),
    fHistMCRapInPtRangeXi(0),
    fHistMCRapPhi(0),
    fHistMCRapInPtRangePhi(0),
    ////////////////////////////////////////////////
    fHistMCPtK0s(0),
    fHistMCPtLambda(0),
    fHistMCPtAntiLambda(0),
    /////////////////////////////////////////////////
    fHistMCPtLambdaFromSigma(0),
    fHistMCPtAntiLambdaFromSigma(0),
    fHistNTimesRecK0s(0),
    fHistNTimesRecLambda(0),
    fHistNTimesRecAntiLambda(0),
    fHistNTimesRecK0sVsPt(0),
    fHistNTimesRecLambdaVsPt(0),
    fHistNTimesRecAntiLambdaVsPt(0),
    fHistNumberEvents(0),
    fHistTrackPerEvent(0),
    fHistTPCMult(0),
    fHistTrackletPerEvent(0),
    fHistMCDaughterTrack(0),
    fHistSPDPrimaryVertexZ(0),
    fHistPrimaryVertexX(0),
    fHistPrimaryVertexY(0),
    fHistPrimaryVertexZ(0),
    fHistPrimaryVertexResX(0),
    fHistPrimaryVertexResY(0),
    fHistPrimaryVertexResZ(0),
    fHistPrimaryVertexPosXV0events(0), 
    fHistPrimaryVertexPosYV0events(0), 
    fHistPrimaryVertexPosZV0events(0),
    fHistDaughterPt(0),
    ////////////////////////////////////////

    fHistDcaPosToPrimVertexK0(0),
    fHistDcaNegToPrimVertexK0(0),
    fHistRadiusV0K0(0),
    fHistDecayLengthV0K0(0),
    fHistDcaV0DaughtersK0(0),
    fHistChi2K0(0),
    fHistCosPointAngleK0(0),

    fHistDcaPosToPrimVertexK0vsMassK0(0),
    fHistDcaNegToPrimVertexK0vsMassK0(0),
    fHistRadiusV0K0vsMassK0(0),
    fHistDecayLengthV0K0vsMassK0(0),
    fHistDcaV0DaughtersK0vsMassK0(0),
    fHistCosPointAngleK0vsMassK0(0),

    
      fHistDcaPosToPrimVertexK0vsMassK0pt1(0),
      fHistDcaNegToPrimVertexK0vsMassK0pt1(0),
      fHistRadiusV0K0vsMassK0pt1(0),
      fHistDecayLengthV0K0vsMassK0pt1(0),
      fHistDcaV0DaughtersK0vsMassK0pt1(0),
      fHistCosPointAngleK0vsMassK0pt1(0),

      fHistDcaPosToPrimVertexK0vsMassK0pt2(0),
      fHistDcaNegToPrimVertexK0vsMassK0pt2(0),
      fHistRadiusV0K0vsMassK0pt2(0),
      fHistDecayLengthV0K0vsMassK0pt2(0),
      fHistDcaV0DaughtersK0vsMassK0pt2(0),
      fHistCosPointAngleK0vsMassK0pt2(0),

      fHistDcaPosToPrimVertexK0vsMassK0pt3(0),
      fHistDcaNegToPrimVertexK0vsMassK0pt3(0),
      fHistRadiusV0K0vsMassK0pt3(0),
      fHistDecayLengthV0K0vsMassK0pt3(0),
      fHistDcaV0DaughtersK0vsMassK0pt3(0),
      fHistCosPointAngleK0vsMassK0pt3(0),
    
    /////////////////////////////////////////

    fHistDcaPosToPrimVertexL(0),
    fHistDcaNegToPrimVertexL(0),
    fHistRadiusV0L(0),
    fHistDecayLengthV0L(0),
    fHistDcaV0DaughtersL(0),
    fHistChi2L(0),
    fHistCosPointAngleL(0),

    fHistDcaPosToPrimVertexLvsMassL(0),
    fHistDcaNegToPrimVertexLvsMassL(0),
    fHistRadiusV0LvsMassL(0),
    fHistDecayLengthV0LvsMassL(0),
    fHistDcaV0DaughtersLvsMassL(0),
    fHistCosPointAngleLvsMassL(0),


    
      fHistDcaPosToPrimVertexLambdaVsMasspt1(0),
      fHistDcaNegToPrimVertexLambdaVsMasspt1(0),
      fHistRadiusV0LambdaVsMasspt1(0),
      fHistDecayLengthV0LambdaVsMasspt1(0),
      fHistDcaV0DaughtersLambdaVsMasspt1(0),
      fHistCosPointAngleLambdaVsMasspt1(0),

      fHistDcaPosToPrimVertexLambdaVsMasspt2(0),
      fHistDcaNegToPrimVertexLambdaVsMasspt2(0),
      fHistRadiusV0LambdaVsMasspt2(0),
      fHistDecayLengthV0LambdaVsMasspt2(0),
      fHistDcaV0DaughtersLambdaVsMasspt2(0),
      fHistCosPointAngleLambdaVsMasspt2(0),

      fHistDcaPosToPrimVertexLambdaVsMasspt3(0),
      fHistDcaNegToPrimVertexLambdaVsMasspt3(0),
      fHistRadiusV0LambdaVsMasspt3(0),
      fHistDecayLengthV0LambdaVsMasspt3(0),
      fHistDcaV0DaughtersLambdaVsMasspt3(0),
      fHistCosPointAngleLambdaVsMasspt3(0),

    ////////////AntiLambda
    fHistDcaPosToPrimVertexAntiL(0),
    fHistDcaNegToPrimVertexAntiL(0),
    fHistRadiusV0AntiL(0),
    fHistDecayLengthV0AntiL(0),
    fHistDcaV0DaughtersAntiL(0),
    fHistChi2AntiL(0),
    fHistCosPointAngleAntiL(0),

    fHistDcaPosToPrimVertexAntiLvsMass(0),
    fHistDcaNegToPrimVertexAntiLvsMass(0),
    fHistRadiusV0AntiLvsMass(0),
    fHistDecayLengthV0AntiLvsMass(0),
    fHistDcaV0DaughtersAntiLvsMass(0),
    fHistCosPointAngleAntiLvsMass(0),

    
    
      fHistDcaPosToPrimVertexAntiLVsMasspt1(0),
      fHistDcaNegToPrimVertexAntiLVsMasspt1(0),
      fHistRadiusV0AntiLVsMasspt1(0),
      fHistDecayLengthV0AntiLVsMasspt1(0),
      fHistDcaV0DaughtersAntiLVsMasspt1(0),
      fHistCosPointAngleAntiLVsMasspt1(0),

      fHistDcaPosToPrimVertexAntiLVsMasspt2(0),
      fHistDcaNegToPrimVertexAntiLVsMasspt2(0),
      fHistRadiusV0AntiLVsMasspt2(0),
      fHistDecayLengthV0AntiLVsMasspt2(0),
      fHistDcaV0DaughtersAntiLVsMasspt2(0),
      fHistCosPointAngleAntiLVsMasspt2(0),

      fHistDcaPosToPrimVertexAntiLVsMasspt3(0),
      fHistDcaNegToPrimVertexAntiLVsMasspt3(0),
      fHistRadiusV0AntiLVsMasspt3(0),
      fHistDecayLengthV0AntiLVsMasspt3(0),
      fHistDcaV0DaughtersAntiLVsMasspt3(0),
      fHistCosPointAngleAntiLVsMasspt3(0),


    /////////////////////////////////////////

    fHistV0Multiplicity(0),
    fHistMassK0(0),
    fHistMassLambda(0),
    fHistMassAntiLambda(0),
    fHistMassVsRadiusK0(0),
    fHistMassVsRadiusLambda(0),
    fHistMassVsRadiusAntiLambda(0),
    /////////////////////////////////////////////
    fHistPtVsMassK0(0),
    fHistPtVsMassLambda(0),
    fHistPtVsMassAntiLambda(0),
    ///////////////////////////////////////////////////
    fHistArmenterosPodolanski(0),
    fHistK0sMassVsLambdaMass(0),
    fHistTPCsignalPt(0),
    fHistNSigmaProton(0),
    fHistNsigmaPosPionAntiLambda(0),
    fHistNsigmaNegProtonAntiLambda(0),
    fHistNsigmaPosProtonLambda(0),
    fHistNsigmaNegPionLambda(0),
    fHistNsigmaPosProtonAntiLambda(0),
    fHistNsigmaNegPionAntiLambda(0),
    fHistNsigmaPosPionK0(0),
    fHistNsigmaNegPionK0(0),
    fHistAsMcRapK0(0),
    fHistAsMcRapLambda(0),
    fHistAsMcRapAntiLambda(0),
    ///////////////////////////////////
    fHistAsMcPtK0(0),
    fHistAsMcPtLambda(0),
    fHistAsMcPtAntiLambda(0),
    /////////////////////////////////////
    fHistAsMcPtZoomK0(0),
    fHistAsMcPtZoomLambda(0),
    fHistAsMcPtZoomAntiLambda(0),
    fHistAsMcProdRadiusK0(0),
    fHistAsMcProdRadiusLambda(0),
    fHistAsMcProdRadiusAntiLambda(0),
    fHistAsMcProdRadiusXvsYK0s(0),
    fHistAsMcProdRadiusXvsYLambda(0),
    fHistAsMcProdRadiusXvsYAntiLambda(0),
    fHistPidMcMassK0(0),
    fHistPidMcMassLambda(0),
    fHistPidMcMassAntiLambda(0),
    fHistAsMcMassK0(0),
    fHistAsMcMassLambda(0),
    fHistAsMcMassAntiLambda(0),
    fHistAsMcPtVsMassK0(0),
    fHistAsMcPtVsMassLambda(0),
    fHistAsMcPtVsMassAntiLambda(0),
    fHistAsMcMassVsRadiusK0(0),
    fHistAsMcMassVsRadiusLambda(0),
    fHistAsMcMassVsRadiusAntiLambda(0),
    fHistAsMcResxK0(0),
    fHistAsMcResyK0(0),
    fHistAsMcReszK0(0),
    fHistAsMcResrVsRadiusK0(0),
    fHistAsMcReszVsRadiusK0(0),
    fHistAsMcResxLambda(0),
    fHistAsMcResyLambda(0),
    fHistAsMcReszLambda(0),
    fHistAsMcResrVsRadiusLambda(0),
    fHistAsMcReszVsRadiusLambda(0),
    fHistAsMcResxAntiLambda(0),
    fHistAsMcResyAntiLambda(0),
    fHistAsMcReszAntiLambda(0),
    fHistAsMcResrVsRadiusAntiLambda(0),
    fHistAsMcReszVsRadiusAntiLambda(0),
    fHistAsMcResPtK0(0),
    fHistAsMcResPtLambda(0),
    fHistAsMcResPtAntiLambda(0),
    fHistAsMcResPtVsRapK0(0),
    fHistAsMcResPtVsRapLambda(0),
    fHistAsMcResPtVsRapAntiLambda(0),
    fHistAsMcResPtVsPtK0(0),
    fHistAsMcResPtVsPtLambda(0),
    fHistAsMcResPtVsPtAntiLambda(0),
    fHistAsMcMotherPdgCodeK0s(0),
    fHistAsMcMotherPdgCodeLambda(0),
    fHistAsMcMotherPdgCodeAntiLambda(0),
    fHistAsMcPtLambdaFromSigma(0),
    fHistAsMcPtAntiLambdaFromSigma(0),
    fHistAsMcSecondaryPtVsRapK0s(0),
    fHistAsMcSecondaryPtVsRapLambda(0),
    fHistAsMcSecondaryPtVsRapAntiLambda(0),
    fHistAsMcSecondaryProdRadiusK0s(0),
    fHistAsMcSecondaryProdRadiusLambda(0),
    fHistAsMcSecondaryProdRadiusAntiLambda(0),
    fHistAsMcSecondaryProdRadiusXvsYK0s(0),
    fHistAsMcSecondaryProdRadiusXvsYLambda(0),
    fHistAsMcSecondaryProdRadiusXvsYAntiLambda(0),
    fHistAsMcSecondaryMotherPdgCodeK0s(0),
    fHistAsMcSecondaryMotherPdgCodeLambda(0),
    fHistAsMcSecondaryMotherPdgCodeAntiLambda(0),
    fHistAsMcSecondaryPtLambdaFromSigma(0),
    fHistAsMcSecondaryPtAntiLambdaFromSigma(0)
    
{
  // Constructor

  // Define output slots only here
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliAnalysisCentralitySelector::Class());
  DefineOutput(3, AliESDtrackCuts::Class());
}
AliAnalysisTaskPerformanceStrange::~AliAnalysisTaskPerformanceStrange() {
  //
  // Destructor
  //
  if (fListHist && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())  { delete fListHist;     fListHist = 0x0;    }
  if (fCentrSelector && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())  { delete fCentrSelector;    fCentrSelector = 0x0;    }
  if (fTracksCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())  { delete fTracksCuts;     fTracksCuts = 0x0;    }


}
//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::UserCreateOutputObjects() 
{

  //******************
  // Create histograms
  //*******************
  fListHist = new TList();
  fListHist->SetOwner();
  //fListHistCuts = new TList();
  //fListHistCuts->SetOwner();

  // Bo: tbd: condition before allocation (i.e. if (!fHistMCMultiplicityPrimary){...} for each histo...

  //***************
  // MC histograms
  //***************
 
  // Primary Vertex:
  fHistMCPrimaryVertexX          = new TH1F("h1MCPrimaryVertexX", "MC Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistMCPrimaryVertexX);

  fHistMCPrimaryVertexY          = new TH1F("h1MCPrimaryVertexY", "MC Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistMCPrimaryVertexY);

  fHistMCPrimaryVertexZ          = new TH1F("h1MCPrimaryVertexZ", "MC Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fListHist->Add(fHistMCPrimaryVertexZ);
  
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
  fHistMCPtAllK0s                      = new TH1F("h1MCPtAllK0s", "Non-primary MC K^{0};p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllK0s);

  fHistMCPtAllLambda                   = new TH1F("h1MCPtAllLambda", "Non-primary MC #Lambda^{0};p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllLambda);

  fHistMCPtAllAntiLambda               = new TH1F("h1MCPtAllAntiLambda", "Non-primary MC #bar{#Lambda}^{0};p_{t} (GeV/c);Counts",240,0,12);
  fListHist->Add(fHistMCPtAllAntiLambda);

  // Production Radius
  fHistMCProdRadiusK0s                 = new TH1F("h1MCProdRadiusK0s", "MC K^{0} Production Radius;r (cm);Count", 400, -2, 2);
  fListHist->Add(fHistMCProdRadiusK0s);

  fHistMCProdRadiusLambda              = new TH1F("h1MCProdRadiusLambda", "MC #Lambda^{0} Production Radius;r (cm);Count", 400, -2, 2);
  fListHist->Add(fHistMCProdRadiusLambda);

  fHistMCProdRadiusAntiLambda         = new TH1F("h1MCProdRadiusAntiLambda", "MC #bar{#Lambda}^{0} Production Radius;r (cm);Count", 400, -2, 2);
  fListHist->Add(fHistMCProdRadiusAntiLambda);

  // Rapidity distribution:
  fHistMCRapK0s                 = new TH1F("h1MCRapK0s", "K^{0};y",160,-4,4);
  fListHist->Add(fHistMCRapK0s);

  fHistMCRapInPtRangeK0s        = new TH1F("h1MCRapInPtRangeK0s", "K^{0};y",160,-4,4);
  fListHist->Add(fHistMCRapInPtRangeK0s);

  fHistMCRapLambda              = new TH1F("h1MCRapLambda", "#Lambda;y",160,-4,4);
  fListHist->Add(fHistMCRapLambda);

  fHistMCRapInPtRangeLambda     = new TH1F("h1MCRapInPtRangeLambda", "#Lambda;y",160,-4,4);
  fListHist->Add(fHistMCRapInPtRangeLambda);

  fHistMCRapAntiLambda          = new TH1F("h1MCRapAntiLambda", "#bar{#Lambda};y",160,-4,4);
  fListHist->Add(fHistMCRapAntiLambda);

  fHistMCRapInPtRangeAntiLambda = new TH1F("h1MCRapInPtRangeAntiLambda", "#bar{#Lambda};y",160,-4,4);
  fListHist->Add(fHistMCRapInPtRangeAntiLambda);

  fHistMCRapXi                  = new TH1F("h1MCRapXi", "Xi;y",160,-4,4);
  fListHist->Add(fHistMCRapXi);

  fHistMCRapInPtRangeXi         = new TH1F("h1MCRapInPtRangeXi", "Xi;y",160,-4,4);
  fListHist->Add(fHistMCRapInPtRangeXi);

  fHistMCRapPhi                  = new TH1F("h1MCRapPhi", "Phi;y",160,-4,4);
  fListHist->Add(fHistMCRapPhi);

  fHistMCRapInPtRangePhi         = new TH1F("h1MCRapInPtRangePhi", "Phi;y",160,-4,4);
  fListHist->Add(fHistMCRapInPtRangePhi);

  // Pt distribution:
  fHistMCPtK0s               = new TH1F("h1MCPtK0s", "K^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtK0s);

  fHistMCPtLambda            = new TH1F("h1MCPtLambda", "#Lambda^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtLambda);

fHistMCPtAntiLambda            = new TH1F("h1MCPtAntiLambda", "#AntiLambda^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtAntiLambda);

  // Pt distribution of Lambda coming from Sigma decay
  fHistMCPtLambdaFromSigma      = new TH1F("h1MCPtLambdaFromSigma", "#Lambda^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtLambdaFromSigma);

  fHistMCPtAntiLambdaFromSigma  = new TH1F("h1MCPtAntiLambdaFromSigma", "#Lambda^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtAntiLambdaFromSigma);
 
  // Multiple reconstruction studies
  fHistNTimesRecK0s             = new TH1F("h1NTimesRecK0s","number of times a K0s is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
  fListHist->Add(fHistNTimesRecK0s);

  fHistNTimesRecLambda          = new TH1F("h1NTimesRecLambda","number of times a Lambda is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
  fListHist->Add(fHistNTimesRecLambda);

  fHistNTimesRecAntiLambda      = new TH1F("h1NTimesRecAntiLambda","number of times an AntiLambda is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
  fListHist->Add(fHistNTimesRecAntiLambda);

  fHistNTimesRecK0sVsPt         = new TH2F("h2NTimesRecK0sVsPt","NTimes versus Pt, K^{0} in -1<y<1;p_{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
  fListHist->Add(fHistNTimesRecK0sVsPt);

  fHistNTimesRecLambdaVsPt      = new TH2F("h2NTimesRecLambdaVsPt","NTimes versus Pt, #Lambda^{0} in -1<y<1;p_{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
  fListHist->Add(fHistNTimesRecLambdaVsPt);

  fHistNTimesRecAntiLambdaVsPt  = new TH2F("h2NTimesRecAntiLambdaVsPt","NTimes versus Pt, #bar{#Lambda}^{0} in -1<y<1;p_{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
  fListHist->Add(fHistNTimesRecAntiLambdaVsPt);

  //***********************************
  // Reconstructed particles histograms
  //***********************************

  // Number of events;
  fHistNumberEvents           = new TH1F("h1NumberEvents", "Number of events; index;Number of Events",10,0,10);
  fListHist->Add(fHistNumberEvents);

  // multiplicity
  fHistTrackPerEvent           = new TH1F("h1TrackPerEvent", "Tracks per event;Number of Tracks;Number of Events",10000,0,10000);
  fListHist->Add(fHistTrackPerEvent);

  fHistTPCMult           = new TH1F("h1HistTPCMult", "TPC tracks per event;Number of Tracks;Number of Events",10000,0,10000);
  fListHist->Add(fHistTPCMult);


  fHistTrackletPerEvent       = new TH1F("h1TrackletPerEvent", "Number of tracklets;Number of tracklets per events;Number of events",1000,0,1000);
  fListHist->Add(fHistTrackletPerEvent);

  fHistMCDaughterTrack         = new TH1F("h1MCDaughterTrack","Distribution of mc id for daughters;id tags;Counts",15,0,15);
  fListHist->Add(fHistMCDaughterTrack);

   fHistTPCTracks               = new TH1F("h1TPCTracks","Distribution of TPC tracks;Number of TPC tracks:Number of events",1000,0,10000);
  fListHist->Add(fHistTPCTracks);

  // Primary Vertex:
  fHistSPDPrimaryVertexZ          = new TH1F("h1SPDPrimaryVertexZ", "SPD Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fListHist->Add(fHistSPDPrimaryVertexZ);

  fHistPrimaryVertexX          = new TH1F("h1PrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexX);

  fHistPrimaryVertexY          = new TH1F("h1PrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexY);

  fHistPrimaryVertexZ          = new TH1F("h1PrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fListHist->Add(fHistPrimaryVertexZ);


  // Primary vertex resolution
  fHistPrimaryVertexResX          = new TH1F("h1PrimaryVertexResX", "Primary Vertex Resolution X;Primary Vertex Resolution X (cm);Events",100,-0.25,0.25);
  fListHist->Add(fHistPrimaryVertexResX);

  fHistPrimaryVertexResY          = new TH1F("h1PrimaryVertexResY", "Primary Vertex Resolution Y;Primary Vertex Resolution Y (cm);Events",100,-0.25,0.25);
  fListHist->Add(fHistPrimaryVertexResY);

  fHistPrimaryVertexResZ          = new TH1F("h1PrimaryVertexResZ", "Primary Vertex Resolution Z;Primary Vertex Resolution Z (cm);Events",200,-0.25,0.25);
  fListHist->Add(fHistPrimaryVertexResZ);
  

  // Primary Vertex in events with V0 candidates:
  fHistPrimaryVertexPosXV0events       = new TH1F("h1PrimaryVertexPosXV0events", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexPosXV0events);
  fHistPrimaryVertexPosYV0events       = new TH1F("h1PrimaryVertexPosYV0events", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexPosYV0events);
  fHistPrimaryVertexPosZV0events       = new TH1F("h1PrimaryVertexPosZV0events", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20.0,20.0);
  fListHist->Add(fHistPrimaryVertexPosZV0events);

  // Daughters Pt:
  fHistDaughterPt              = new TH2F("h2DaughterPt", "Daughter Pt;Positive Daughter Pt; Negative Daughter Pt",200,0,2,200,0,2);
  fListHist->Add(fHistDaughterPt);

  //////K0s///////////////// 2D histos: cut vs on fly status////

  fHistDcaPosToPrimVertexK0      = new TH2F("h2DcaPosToPrimVertexK0", "Positive V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertexK0);

  fHistDcaNegToPrimVertexK0      = new TH2F("h2DcaNegToPrimVertexK0", "Negative V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertexK0);


  fHistRadiusV0K0                = new TH2F("h2RadiusV0K0", "Radius;Radius(cm);Status",500,0,500,2,-0.5,1.5);
  fListHist->Add(fHistRadiusV0K0);

  fHistDecayLengthV0K0           = new TH2F("h2DecayLengthV0K0", "V0s decay Length;decay length(cm);Status", 500, 0, 500,2,-0.5,1.5);
  fListHist->Add(fHistDecayLengthV0K0);

  fHistDcaV0DaughtersK0          = new TH2F("h2DcaV0DaughtersK0", "DCA between daughters;dca(cm);Status", 300, 0, 3.0,2,-0.5,1.5);
  fListHist->Add(fHistDcaV0DaughtersK0);

  fHistChi2K0                    = new TH2F("h2Chi2K0", "V0s chi2;chi2;Status", 1000, 0, 0.1,2,-0.5,1.5);
  fListHist->Add(fHistChi2K0);

  fHistCosPointAngleK0           = new TH2F("h2CosPointAngleK0", "Cosine of V0's pointing angle", 200,0.99,1.01,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngleK0);


  ////////////K0s///////////////// 2D histos: cut vs mass////


  fHistDcaPosToPrimVertexK0vsMassK0 = new TH2F("h2DcaPosToPrimVertexK0vsMassK0", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
  fListHist->Add(fHistDcaPosToPrimVertexK0vsMassK0);

  fHistDcaNegToPrimVertexK0vsMassK0 = new TH2F("h2DcaNegToPrimVertexK0vsMassK0", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
  fListHist->Add(fHistDcaNegToPrimVertexK0vsMassK0);


  fHistRadiusV0K0vsMassK0           = new TH2F("h2RadiusV0K0vsMassK0", "Radius;Radius(cm);K0s inv. mass",110,0,110,200,0.4,0.6);
  fListHist->Add(fHistRadiusV0K0vsMassK0);

  fHistDecayLengthV0K0vsMassK0      = new TH2F("h2DecayLengthV0K0vsMassK0", "V0s decay Length;decay length(cm);K0s inv. mass", 100, 0, 100,200,0.4,0.6);
  fListHist->Add(fHistDecayLengthV0K0vsMassK0);

  fHistDcaV0DaughtersK0vsMassK0     = new TH2F("h2DcaV0DaughtersK0vsMassK0", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,200,0.4,0.6);
  fListHist->Add(fHistDcaV0DaughtersK0vsMassK0);


  fHistCosPointAngleK0vsMassK0      = new TH2F("h2CosPointAngleK0vsMassK0", "Cosine of V0's pointing angle", 200,0.997,1.007,200,0.4,0.6);
  fListHist->Add(fHistCosPointAngleK0vsMassK0);
    //// pt1
      fHistDcaPosToPrimVertexK0vsMassK0pt1 = new TH2F("h2DcaPosToPrimVertexK0vsMassK0pt1", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
      fListHist->Add(fHistDcaPosToPrimVertexK0vsMassK0pt1);

      fHistDcaNegToPrimVertexK0vsMassK0pt1 = new TH2F("h2DcaNegToPrimVertexK0vsMassK0pt1", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
      fListHist->Add(fHistDcaNegToPrimVertexK0vsMassK0pt1);

      fHistRadiusV0K0vsMassK0pt1           = new TH2F("h2RadiusV0K0vsMassK0pt1", "Radius;Radius(cm);K0s inv. mass",110,0,110,200,0.4,0.6);
      fListHist->Add(fHistRadiusV0K0vsMassK0pt1);

      fHistDecayLengthV0K0vsMassK0pt1      = new TH2F("h2DecayLengthV0K0vsMassK0pt1", "V0s decay Length;decay length(cm);K0s inv. mass",100,0,100,200,0.4,0.6);
      fListHist->Add(fHistDecayLengthV0K0vsMassK0pt1);

      fHistDcaV0DaughtersK0vsMassK0pt1     = new TH2F("h2DcaV0DaughtersK0vsMassK0pt1", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,200,0.4,0.6);
      fListHist->Add(fHistDcaV0DaughtersK0vsMassK0pt1);

      fHistCosPointAngleK0vsMassK0pt1      = new TH2F("h2CosPointAngleK0vsMassK0pt1", "Cosine of V0's pointing angle", 200,0.997,1.007,200,0.4,0.6);
      fListHist->Add(fHistCosPointAngleK0vsMassK0pt1);

      /// pt2
      fHistDcaPosToPrimVertexK0vsMassK0pt2 = new TH2F("h2DcaPosToPrimVertexK0vsMassK0pt2", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
      fListHist->Add(fHistDcaPosToPrimVertexK0vsMassK0pt2);

      fHistDcaNegToPrimVertexK0vsMassK0pt2 = new TH2F("h2DcaNegToPrimVertexK0vsMassK0pt2", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
      fListHist->Add(fHistDcaNegToPrimVertexK0vsMassK0pt2);

      fHistRadiusV0K0vsMassK0pt2           = new TH2F("h2RadiusV0K0vsMassK0pt2", "Radius;Radius(cm);K0s inv. mass",110,0,110,200,0.4,0.6);
      fListHist->Add(fHistRadiusV0K0vsMassK0pt2);

      fHistDecayLengthV0K0vsMassK0pt2      = new TH2F("h2DecayLengthV0K0vsMassK0pt2", "V0s decay Length;decay length(cm);K0s inv. mass", 100, 0, 100,200,0.4,0.6);
      fListHist->Add(fHistDecayLengthV0K0vsMassK0pt2);

      fHistDcaV0DaughtersK0vsMassK0pt2     = new TH2F("h2DcaV0DaughtersK0vsMassK0pt2", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,200,0.4,0.6);
      fListHist->Add(fHistDcaV0DaughtersK0vsMassK0pt2);
  
      fHistCosPointAngleK0vsMassK0pt2      = new TH2F("h2CosPointAngleK0vsMassK0pt2", "Cosine of V0's pointing angle", 200,0.997,1.007,200,0.4,0.6);
      fListHist->Add(fHistCosPointAngleK0vsMassK0pt2);

      /// pt3
      fHistDcaPosToPrimVertexK0vsMassK0pt3 = new TH2F("h2DcaPosToPrimVertexK0vsMassK0pt3", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
      fListHist->Add(fHistDcaPosToPrimVertexK0vsMassK0pt3);

      fHistDcaNegToPrimVertexK0vsMassK0pt3 = new TH2F("h2DcaNegToPrimVertexK0vsMassK0pt3", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
      fListHist->Add(fHistDcaNegToPrimVertexK0vsMassK0pt3);

      fHistRadiusV0K0vsMassK0pt3           = new TH2F("h2RadiusV0K0vsMassK0pt3", "Radius;Radius(cm);K0s inv. mass",110,0,110,200,0.4,0.6);
      fListHist->Add(fHistRadiusV0K0vsMassK0pt3);

      fHistDecayLengthV0K0vsMassK0pt3      = new TH2F("h2DecayLengthV0K0vsMassK0pt3", "V0s decay Length;decay length(cm);K0s inv. mass", 100, 0, 100,200,0.4,0.6);
      fListHist->Add(fHistDecayLengthV0K0vsMassK0pt3);

      fHistDcaV0DaughtersK0vsMassK0pt3     = new TH2F("h2DcaV0DaughtersK0vsMassK0pt3", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,200,0.4,0.6);
      fListHist->Add(fHistDcaV0DaughtersK0vsMassK0pt3);
  
      fHistCosPointAngleK0vsMassK0pt3      = new TH2F("h2CosPointAngleK0vsMassK0pt3", "Cosine of V0's pointing angle", 200,0.997,1.007,200,0.4,0.6);
      fListHist->Add(fHistCosPointAngleK0vsMassK0pt3);
  
  //////////Lambda////////////// 2D histos: cut vs on fly status////

  fHistDcaPosToPrimVertexL      = new TH2F("h2DcaPosToPrimVertexL", "Positive V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertexL);

  fHistDcaNegToPrimVertexL      = new TH2F("h2DcaNegToPrimVertexL", "Negative V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertexL);


  fHistRadiusV0L                = new TH2F("h2RadiusV0L", "Radius;Radius(cm);Status",100,0,110,2,-0.5,1.5);
  fListHist->Add(fHistRadiusV0L);

  fHistDecayLengthV0L           = new TH2F("h2DecayLengthV0L", "V0s decay Length;decay length(cm);Status", 500, 0, 500,2,-0.5,1.5);
  fListHist->Add(fHistDecayLengthV0L);

  fHistDcaV0DaughtersL          = new TH2F("h2DcaV0DaughtersL", "DCA between daughters;dca(cm);Status", 300, 0, 3.0,2,-0.5,1.5);
  fListHist->Add(fHistDcaV0DaughtersL);

  fHistChi2L                    = new TH2F("h2Chi2L", "V0s chi2;chi2;Status", 100, 0, 0.10,2,-0.5,1.5);
  fListHist->Add(fHistChi2L);

  fHistCosPointAngleL           = new TH2F("h2CosPointAngleL", "Cosine of V0's pointing angle", 200,0.99,1.01,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngleL);

  //////////Lambda////////////// 2D histos: cut vs mass////

  fHistDcaPosToPrimVertexLvsMassL      = new TH2F("h2DcaPosToPrimVertexLvsMassL", "Positive V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
  fListHist->Add(fHistDcaPosToPrimVertexLvsMassL);

  fHistDcaNegToPrimVertexLvsMassL      = new TH2F("h2DcaNegToPrimVertexLvsMassL", "Negative V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
  fListHist->Add(fHistDcaNegToPrimVertexLvsMassL);


  fHistRadiusV0LvsMassL                = new TH2F("h2RadiusV0LvsMassL", "Radius;Radius(cm);Status",110,0,110,140, 1.06, 1.2);
  fListHist->Add(fHistRadiusV0LvsMassL);

  fHistDecayLengthV0LvsMassL           = new TH2F("h2DecayLengthV0LvsMassL", "V0s decay Length;decay length(cm);Status", 120, 0, 120,140, 1.06, 1.2);
  fListHist->Add(fHistDecayLengthV0LvsMassL);

  fHistDcaV0DaughtersLvsMassL          = new TH2F("h2DcaV0DaughtersLvsMassL", "DCA between daughters;dca(cm);Status", 110, 0, 1.1,140, 1.06, 1.2);
  fListHist->Add(fHistDcaV0DaughtersLvsMassL);

  fHistCosPointAngleLvsMassL           = new TH2F("h2CosPointAngleLvsMassL", "Cosine of V0's pointing angle", 200,0.997,1.007,140, 1.06, 1.2);
  fListHist->Add(fHistCosPointAngleLvsMassL);

    //// pt1
      fHistDcaPosToPrimVertexLambdaVsMasspt1 = new TH2F("h2DcaPosToPrimVertexLambdaVsMasspt1", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaPosToPrimVertexLambdaVsMasspt1);

      fHistDcaNegToPrimVertexLambdaVsMasspt1 = new TH2F("h2DcaNegToPrimVertexLambdaVsMasspt1", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaNegToPrimVertexLambdaVsMasspt1);

      fHistRadiusV0LambdaVsMasspt1           = new TH2F("h2RadiusV0LambdaVsMasspt1", "Radius;Radius(cm);K0s inv. mass",110,0,110,140,1.06,1.2);
      fListHist->Add(fHistRadiusV0LambdaVsMasspt1);

      fHistDecayLengthV0LambdaVsMasspt1      = new TH2F("h2DecayLengthV0LambdaVsMasspt1", "V0s decay Length;decay length(cm);K0s inv. mass",100,0,100,140,1.06,1.2);
      fListHist->Add(fHistDecayLengthV0LambdaVsMasspt1);

      fHistDcaV0DaughtersLambdaVsMasspt1     = new TH2F("h2DcaV0DaughtersLambdaVsMasspt1", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,140,1.06,1.2);
      fListHist->Add(fHistDcaV0DaughtersLambdaVsMasspt1);

      fHistCosPointAngleLambdaVsMasspt1      = new TH2F("h2CosPointAngleLambdaVsMasspt1", "Cosine of V0's pointing angle", 200,0.997,1.007,140,1.06,1.2);
      fListHist->Add(fHistCosPointAngleLambdaVsMasspt1);

      /// pt2
      fHistDcaPosToPrimVertexLambdaVsMasspt2 = new TH2F("h2DcaPosToPrimVertexLambdaVsMasspt2", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaPosToPrimVertexLambdaVsMasspt2);

      fHistDcaNegToPrimVertexLambdaVsMasspt2 = new TH2F("h2DcaNegToPrimVertexLambdaVsMasspt2", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaNegToPrimVertexLambdaVsMasspt2);

      fHistRadiusV0LambdaVsMasspt2           = new TH2F("h2RadiusV0LambdaVsMasspt2", "Radius;Radius(cm);K0s inv. mass",110,0,110,140,1.06,1.2);
      fListHist->Add(fHistRadiusV0LambdaVsMasspt2);

      fHistDecayLengthV0LambdaVsMasspt2      = new TH2F("h2DecayLengthV0LambdaVsMasspt2", "V0s decay Length;decay length(cm);K0s inv. mass", 100, 0, 100,140,1.06,1.2);
      fListHist->Add(fHistDecayLengthV0LambdaVsMasspt2);

      fHistDcaV0DaughtersLambdaVsMasspt2     = new TH2F("h2DcaV0DaughtersLambdaVsMasspt2", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,140,1.06,1.2);
      fListHist->Add(fHistDcaV0DaughtersLambdaVsMasspt2);
  
      fHistCosPointAngleLambdaVsMasspt2      = new TH2F("h2CosPointAngleLambdaVsMasspt2", "Cosine of V0's pointing angle", 200,0.997,1.007,140,1.06,1.2);
      fListHist->Add(fHistCosPointAngleLambdaVsMasspt2);

      /// pt3
      fHistDcaPosToPrimVertexLambdaVsMasspt3 = new TH2F("h2DcaPosToPrimVertexLambdaVsMasspt3", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaPosToPrimVertexLambdaVsMasspt3);

      fHistDcaNegToPrimVertexLambdaVsMasspt3 = new TH2F("h2DcaNegToPrimVertexLambdaVsMasspt3", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaNegToPrimVertexLambdaVsMasspt3);

      fHistRadiusV0LambdaVsMasspt3           = new TH2F("h2RadiusV0LambdaVsMasspt3", "Radius;Radius(cm);K0s inv. mass",110,0,110,140,1.06,1.2);
      fListHist->Add(fHistRadiusV0LambdaVsMasspt3);

      fHistDecayLengthV0LambdaVsMasspt3      = new TH2F("h2DecayLengthV0LambdaVsMasspt3", "V0s decay Length;decay length(cm);K0s inv. mass", 100, 0, 100,140,1.06,1.2);
      fListHist->Add(fHistDecayLengthV0LambdaVsMasspt3);

      fHistDcaV0DaughtersLambdaVsMasspt3     = new TH2F("h2DcaV0DaughtersLambdaVsMasspt3", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,140,1.06,1.2);
      fListHist->Add(fHistDcaV0DaughtersLambdaVsMasspt3);
  
      fHistCosPointAngleLambdaVsMasspt3      = new TH2F("h2CosPointAngleLambdaVsMasspt3", "Cosine of V0's pointing angle", 200,0.997,1.007,140,1.06,1.2);
      fListHist->Add(fHistCosPointAngleLambdaVsMasspt3);

  //////////AntiLambda////////////// 2D histos: cut vs on fly status////

  fHistDcaPosToPrimVertexAntiL      = new TH2F("h2DcaPosToPrimVertexAntiL", "Positive V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertexAntiL);

  fHistDcaNegToPrimVertexAntiL      = new TH2F("h2DcaNegToPrimVertexAntiL", "Negative V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertexAntiL);


  fHistRadiusV0AntiL                = new TH2F("h2RadiusV0AntiL", "Radius;Radius(cm);Status",100,0,110,2,-0.5,1.5);
  fListHist->Add(fHistRadiusV0AntiL);

  fHistDecayLengthV0AntiL           = new TH2F("h2DecayLengthV0AntiL", "V0s decay Length;decay length(cm);Status", 500, 0, 500,2,-0.5,1.5);
  fListHist->Add(fHistDecayLengthV0AntiL);

  fHistDcaV0DaughtersAntiL          = new TH2F("h2DcaV0DaughtersAntiL", "DCA between daughters;dca(cm);Status", 300, 0, 3.0,2,-0.5,1.5);
  fListHist->Add(fHistDcaV0DaughtersAntiL);

  fHistChi2AntiL                    = new TH2F("h2Chi2AntiL", "V0s chi2;chi2;Status", 100, 0, 0.10,2,-0.5,1.5);
  fListHist->Add(fHistChi2AntiL);

  fHistCosPointAngleAntiL           = new TH2F("h2CosPointAngleAntiL", "Cosine of V0's pointing angle", 200,0.99,1.01,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngleAntiL);

  //////////AntiLambda////////////// 2D histos: cut vs mass////

  fHistDcaPosToPrimVertexAntiLvsMass      = new TH2F("h2DcaPosToPrimVertexAntiLvsMass", "Positive V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
  fListHist->Add(fHistDcaPosToPrimVertexAntiLvsMass);

  fHistDcaNegToPrimVertexAntiLvsMass      = new TH2F("h2DcaNegToPrimVertexAntiLvsMass", "Negative V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
  fListHist->Add(fHistDcaNegToPrimVertexAntiLvsMass);


  fHistRadiusV0AntiLvsMass                = new TH2F("h2RadiusV0AntiLvsMass", "Radius;Radius(cm);Status",110,0,110,140, 1.06, 1.2);
  fListHist->Add(fHistRadiusV0AntiLvsMass);

  fHistDecayLengthV0AntiLvsMass           = new TH2F("h2DecayLengthV0AntiLvsMass", "V0s decay Length;decay length(cm);Status", 120, 0, 120,140, 1.06, 1.2);
  fListHist->Add(fHistDecayLengthV0AntiLvsMass);

  fHistDcaV0DaughtersAntiLvsMass          = new TH2F("h2DcaV0DaughtersAntiLvsMass", "DCA between daughters;dca(cm);Status", 110, 0, 1.1,140, 1.06, 1.2);
  fListHist->Add(fHistDcaV0DaughtersAntiLvsMass);

  fHistCosPointAngleAntiLvsMass           = new TH2F("h2CosPointAngleAntiLvsMass", "Cosine of V0's pointing angle", 200,0.997,1.007,140, 1.06, 1.2);
  fListHist->Add(fHistCosPointAngleAntiLvsMass);

    //// pt1
      fHistDcaPosToPrimVertexAntiLVsMasspt1 = new TH2F("h2DcaPosToPrimVertexAntiLVsMasspt1", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaPosToPrimVertexAntiLVsMasspt1);

      fHistDcaNegToPrimVertexAntiLVsMasspt1 = new TH2F("h2DcaNegToPrimVertexAntiLVsMasspt1", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaNegToPrimVertexAntiLVsMasspt1);

      fHistRadiusV0AntiLVsMasspt1           = new TH2F("h2RadiusV0AntiLVsMasspt1", "Radius;Radius(cm);K0s inv. mass",110,0,110,140,1.06,1.2);
      fListHist->Add(fHistRadiusV0AntiLVsMasspt1);

      fHistDecayLengthV0AntiLVsMasspt1      = new TH2F("h2DecayLengthV0AntiLVsMasspt1", "V0s decay Length;decay length(cm);K0s inv. mass",100,0,100,140,1.06,1.2);
      fListHist->Add(fHistDecayLengthV0AntiLVsMasspt1);

      fHistDcaV0DaughtersAntiLVsMasspt1     = new TH2F("h2DcaV0DaughtersAntiLVsMasspt1", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,140,1.06,1.2);
      fListHist->Add(fHistDcaV0DaughtersAntiLVsMasspt1);

      fHistCosPointAngleAntiLVsMasspt1      = new TH2F("h2CosPointAngleAntiLVsMasspt1", "Cosine of V0's pointing angle", 200,0.997,1.007,140,1.06,1.2);
      fListHist->Add(fHistCosPointAngleAntiLVsMasspt1);

      /// pt2
      fHistDcaPosToPrimVertexAntiLVsMasspt2 = new TH2F("h2DcaPosToPrimVertexAntiLVsMasspt2", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaPosToPrimVertexAntiLVsMasspt2);

      fHistDcaNegToPrimVertexAntiLVsMasspt2 = new TH2F("h2DcaNegToPrimVertexAntiLVsMasspt2", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaNegToPrimVertexAntiLVsMasspt2);

      fHistRadiusV0AntiLVsMasspt2           = new TH2F("h2RadiusV0AntiLVsMasspt2", "Radius;Radius(cm);K0s inv. mass",110,0,110,140,1.06,1.2);
      fListHist->Add(fHistRadiusV0AntiLVsMasspt2);

      fHistDecayLengthV0AntiLVsMasspt2      = new TH2F("h2DecayLengthV0AntiLVsMasspt2", "V0s decay Length;decay length(cm);K0s inv. mass", 100, 0, 100,140,1.06,1.2);
      fListHist->Add(fHistDecayLengthV0AntiLVsMasspt2);

      fHistDcaV0DaughtersAntiLVsMasspt2     = new TH2F("h2DcaV0DaughtersAntiLVsMasspt2", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,140,1.06,1.2);
      fListHist->Add(fHistDcaV0DaughtersAntiLVsMasspt2);
  
      fHistCosPointAngleAntiLVsMasspt2      = new TH2F("h2CosPointAngleAntiLVsMasspt2", "Cosine of V0's pointing angle", 200,0.997,1.007,140,1.06,1.2);
      fListHist->Add(fHistCosPointAngleAntiLVsMasspt2);

      /// pt3
      fHistDcaPosToPrimVertexAntiLVsMasspt3 = new TH2F("h2DcaPosToPrimVertexAntiLVsMasspt3", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaPosToPrimVertexAntiLVsMasspt3);

      fHistDcaNegToPrimVertexAntiLVsMasspt3 = new TH2F("h2DcaNegToPrimVertexAntiLVsMasspt3", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,140,1.06,1.2);
      fListHist->Add(fHistDcaNegToPrimVertexAntiLVsMasspt3);

      fHistRadiusV0AntiLVsMasspt3           = new TH2F("h2RadiusV0AntiLVsMasspt3", "Radius;Radius(cm);K0s inv. mass",110,0,110,140,1.06,1.2);
      fListHist->Add(fHistRadiusV0AntiLVsMasspt3);

      fHistDecayLengthV0AntiLVsMasspt3      = new TH2F("h2DecayLengthV0AntiLVsMasspt3", "V0s decay Length;decay length(cm);K0s inv. mass", 100, 0, 100,140,1.06,1.2);
      fListHist->Add(fHistDecayLengthV0AntiLVsMasspt3);

      fHistDcaV0DaughtersAntiLVsMasspt3     = new TH2F("h2DcaV0DaughtersAntiLVsMasspt3", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,140,1.06,1.2);
      fListHist->Add(fHistDcaV0DaughtersAntiLVsMasspt3);
  
      fHistCosPointAngleAntiLVsMasspt3      = new TH2F("h2CosPointAngleAntiLVsMasspt3", "Cosine of V0's pointing angle", 200,0.997,1.007,140,1.06,1.2);
      fListHist->Add(fHistCosPointAngleAntiLVsMasspt3);


  // V0 Multiplicity
  if (!fHistV0Multiplicity) {
    if (fCollidingSystems)
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of Offline V0s;Events", 200, 0, 40000);
    else
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of Offline V0s;Events", 10, 0, 10); 
    fListHist->Add(fHistV0Multiplicity);
  }


  // Mass:
  fHistMassK0                   = new TH1F("h1MassK0", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 200, 0.4, 0.6);
  fListHist->Add(fHistMassK0);

  fHistMassLambda               = new TH1F("h1MassLambda", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 150, 1.05, 1.2);
  fListHist->Add(fHistMassLambda);

  fHistMassAntiLambda           = new TH1F("h1MassAntiLambda", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 150, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambda);

  /*  // invariant mass vs radius
      const Double_t radius[10] = {0.0,2.5,2.9,3.9,7.6,15.0,23.9,37.8,42.8,100.0};
      Int_t lNbinRadius        = 9;
      Int_t lNbinInvMassLambda = 300;
  */
  fHistMassVsRadiusK0           = new TH2F("h2MassVsRadiusK0", "K^{0} candidates;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",200,0,200, 200, 0.4, 0.6);
  fListHist->Add(fHistMassVsRadiusK0);

  //fHistMassVsRadiusLambda       = new TH2F("h2MassVsRadiusLambda", "#Lambda candidates;radius (cm);M(p#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 140, 1.06, 1.2);
  fHistMassVsRadiusLambda       = new TH2F("h2MassVsRadiusLambda", "#Lambda candidates;radius (cm);M(p#pi^{-}) (GeV/c^{2})",200,0,200, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusLambda);


  //fHistMassVsRadiusAntiLambda   = new TH2F("h2MassVsRadiusAntiLambda", "#bar{#Lambda} candidates;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",lNbinRadius,radius, 140, 1.06, 1.2);
  fHistMassVsRadiusAntiLambda   = new TH2F("h2MassVsRadiusAntiLambda", "#bar{#Lambda} candidates;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",200,0,200, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusAntiLambda);


  // Pt Vs Mass
  fHistPtVsMassK0               = new TH2F("h2PtVsMassK0","K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",400, 0.4, 0.6,480,0,12);
  fListHist->Add(fHistPtVsMassK0);

  fHistPtVsMassLambda           = new TH2F("h2PtVsMassLambda","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",280, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistPtVsMassLambda);
  
  fHistPtVsMassAntiLambda           = new TH2F("h2PtVsMassAntiLambda","#AntiLambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",280, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistPtVsMassAntiLambda);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///Armenteros Podolansky
  fHistArmenterosPodolanski     = new TH2F("h2ArmenterosPodolanski","Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);
  fListHist->Add(fHistArmenterosPodolanski);

  ///Inv. Mass K0s vs Inv. Mass. Lambda
  fHistK0sMassVsLambdaMass      = new TH2F("h2HistK0sMassVsLambdaMass","K^{0} vs #Lambda^{0} candidates; M(#pi^{+}#pi^{-}) (GeV/c^{2}); M(p#pi^{-}) (GeV/c^{2})",200, 0.4, 0.6,140, 1.06, 1.2);
  fListHist->Add(fHistK0sMassVsLambdaMass);

  //dE/dx vs P daughters
  fHistTPCsignalPt                            = new TH2F("h2TPCsignalVsPdaughters","TPC signal Vs p_{t} daughters;  p (GeV/c);TPC signal",1000,0,2,1000,0,1000);
  fListHist->Add(fHistTPCsignalPt);
  fHistNSigmaProton                          =new TH1F("h1NSigmaProton","Number of sigmas for proton;;Count",600,0,6);
  fListHist->Add(fHistNSigmaProton);


  //PID
  fHistNsigmaPosPionAntiLambda   = new TH1F("h1NsigmaPosPionAntiLambda", "Positive daughter of Antilambda;NsigmaPion;Counts",25,0,5);
  fListHist->Add(fHistNsigmaPosPionAntiLambda);

  fHistNsigmaNegProtonAntiLambda = new TH1F("h1NsigmaNegProtonAntiLambda", "Negative daughter of Antilambda;NsigmaProton;Counts",25,0,5);
  fListHist->Add(fHistNsigmaNegProtonAntiLambda);
  
  fHistNsigmaPosProtonLambda     = new TH1F("h1NsigmaPosProtonLambda", "Positive daughter of Lambda;NsigmaProton;Counts",25,0,5); 
  fListHist->Add(fHistNsigmaPosProtonLambda);
  
  fHistNsigmaNegPionLambda       = new TH1F("h1NsigmaNegPionLambda", "Negative daughter of Lambda;NsigmaPion;Counts",25,0,5);
  fListHist->Add(fHistNsigmaNegPionLambda);
  
  fHistNsigmaPosProtonAntiLambda     = new TH1F("h1NsigmaPosProtonAntiLambda", "Positive daughter of AntiLambda;NsigmaProton;Counts",25,0,5); 
  fListHist->Add(fHistNsigmaPosProtonAntiLambda);
  
  fHistNsigmaNegPionAntiLambda       = new TH1F("h1NsigmaNegPionAntiLambda", "Negative daughter of AntiLambda;NsigmaPion;Counts",25,0,5);
  fListHist->Add(fHistNsigmaNegPionAntiLambda);
  
  fHistNsigmaPosPionK0           = new TH1F("h1NsigmaPosPionK0", "Positive daughter of K0s;NsigmaPion;Counts",25,0,5);
  fListHist->Add(fHistNsigmaPosPionK0);
  
  fHistNsigmaNegPionK0           = new TH1F("h1NsigmaNegPionK0", "Negative daughter of K0s;NsigmaPion;Counts",25,0,5);
  fListHist->Add(fHistNsigmaNegPionK0);


  //********************************
  // Associated particles histograms
  //********************************

  // Rap distribution
  fHistAsMcRapK0                = new TH1F("h1AsMcRapK0", "K^{0} associated;eta;Counts", 60, -1.5, 1.5);
  fListHist->Add(fHistAsMcRapK0);

  fHistAsMcRapLambda            = new TH1F("h1AsMcRapLambda", "#Lambda^{0} associated;eta;Counts", 60, -1.5, 1.5);
  fListHist->Add(fHistAsMcRapLambda);

  fHistAsMcRapAntiLambda        = new TH1F("h1AsMcRapAntiLambda", "#bar{#Lambda}^{0} associated;eta;Counts", 60, -1.5, 1.5);
  fListHist->Add(fHistAsMcRapAntiLambda);

  //Pt distribution
  fHistAsMcPtK0                = new TH1F("h1AsMcPtK0", "K^{0} associated;p_{t} (GeV/c);Counts", 240,0,12);
  fListHist->Add(fHistAsMcPtK0);

  fHistAsMcPtLambda            = new TH1F("h1AsMcPtLambda", "#Lambda^{0} associated;p_{t} (GeV/c);Counts", 240,0,12);
  fListHist->Add(fHistAsMcPtLambda);

 fHistAsMcPtAntiLambda            = new TH1F("h1AsMcPtAntiLambda", "#AntiLambda^{0} associated;p_{t} (GeV/c);Counts", 240,0,12);
  fListHist->Add(fHistAsMcPtAntiLambda);


  fHistAsMcPtZoomK0            = new TH1F("h1AsMcPtZoomK0", "K^{0} candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomK0);

  fHistAsMcPtZoomLambda        = new TH1F("h1AsMcPtZoomLambda", "#Lambda^{0} candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomLambda);

 fHistAsMcPtZoomAntiLambda        = new TH1F("h1AsMcPtZoomAntiLambda", "#AntiLambda^{0} candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomAntiLambda);

  // Radius distribution
  fHistAsMcProdRadiusK0               = new TH1F("h1AsMcProdRadiusK0", "K^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusK0);

  fHistAsMcProdRadiusLambda           = new TH1F("h1AsMcProdRadiusLambda", "#Lambda^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusLambda);

  fHistAsMcProdRadiusAntiLambda       = new TH1F("h1AsMcProdRadiusAntiLambda", "#bar{#Lambda}^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusAntiLambda);

  fHistAsMcProdRadiusXvsYK0s          = new TH2F("h2AsMcProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYK0s);

  fHistAsMcProdRadiusXvsYLambda       = new TH2F("h2AsMcProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYLambda);

  fHistAsMcProdRadiusXvsYAntiLambda   = new TH2F("h2AsMcProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYAntiLambda);

  // Mass
  fHistPidMcMassK0             = new TH1F("h1PidMcMassK0", "K^{0} MC PId checked;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistPidMcMassK0);

  fHistPidMcMassLambda         = new TH1F("h1PidMcMassLambda", "#Lambda^{0} MC PId checked;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassLambda);
  
  fHistPidMcMassAntiLambda     = new TH1F("h1PidMcMassAntiLambda", "#bar{#Lambda}^{0} MC PId checked;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassAntiLambda);

  fHistAsMcMassK0              = new TH1F("h1AsMcMassK0", "K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistAsMcMassK0);
  
  fHistAsMcMassLambda          = new TH1F("h1AsMcMassLambda", "#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassLambda);

  fHistAsMcMassAntiLambda      = new TH1F("h1AsMcMassAntiLambda", "#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassAntiLambda);

  //Pt versus Mass
  fHistAsMcPtVsMassK0               = new TH2F("h2AsMcPtVsMassK0","K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,240,0,12);
  fListHist->Add(fHistAsMcPtVsMassK0);

  fHistAsMcPtVsMassLambda           = new TH2F("h2AsMcPtVsMassLambda","#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistAsMcPtVsMassLambda);

  fHistAsMcPtVsMassAntiLambda       = new TH2F("h2AsMcPtVsMassAntiLambda","#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistAsMcPtVsMassAntiLambda);


  // invariant mass vs radius
  //fHistAsMcMassVsRadiusK0             = new TH2F("h2AsMcMassVsRadiusK0", "K^{0} associated;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 500, 0.47, 0.52);
  fHistAsMcMassVsRadiusK0             = new TH2F("h2AsMcMassVsRadiusK0", "K^{0} associated;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",200,0,200, 500, 0.47, 0.52);
  fListHist->Add(fHistAsMcMassVsRadiusK0);

  //fHistAsMcMassVsRadiusLambda         = new TH2F("h2AsMcMassVsRadiusLambda", "#Lambda associated;radius (cm);M(p#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, lNbinInvMassLambda, 1.10, 1.13);
  fHistAsMcMassVsRadiusLambda         = new TH2F("h2AsMcMassVsRadiusLambda", "#Lambda associated;radius (cm);M(p#pi^{-}) (GeV/c^{2})",200,0,200, 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusLambda);

  //fHistAsMcMassVsRadiusAntiLambda     = new TH2F("h2AsMcMassVsRadiusAntiLambda", "#bar{#Lambda} associated;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",lNbinRadius,radius,lNbinInvMassLambda , 1.10, 1.13);
  fHistAsMcMassVsRadiusAntiLambda     = new TH2F("h2AsMcMassVsRadiusAntiLambda", "#bar{#Lambda} associated;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",200,0,200 , 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusAntiLambda);
  
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

  // Pt Resolution
  fHistAsMcResPtK0                   = new TH1F("h1AsMcResPtK0","Pt Resolution K^{0};#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtK0);
  
  fHistAsMcResPtLambda               = new TH1F("h1AsMcResPtLambda","Pt Resolution #Lambda^{0};#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtLambda);

  fHistAsMcResPtAntiLambda           = new TH1F("h1AsMcResPtAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtAntiLambda);


  fHistAsMcResPtVsRapK0              = new TH2F("h2AsMcResPtVsRapK0","Pt Resolution K^{0};#Delta Pt;Rap",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapK0);
  
  fHistAsMcResPtVsRapLambda          = new TH2F("h2AsMcResPtVsRapLambda","Pt Resolution #Lambda^{0};#Delta Pt;Rap",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapLambda);

  fHistAsMcResPtVsRapAntiLambda      = new TH2F("h2AsMcResPtVsRapAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Rap",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapAntiLambda);

  fHistAsMcResPtVsPtK0               = new TH2F("h2AsMcResPtVsPtK0","Pt Resolution K^{0};#Delta Pt;Pt",600,-0.15,0.15,240,0,12);
  fListHist->Add(fHistAsMcResPtVsPtK0);
    
  fHistAsMcResPtVsPtLambda           = new TH2F("h2AsMcResPtVsPtLambda","Pt Resolution #Lambda^{0};#Delta Pt;Pt",600,-0.15,0.15,240,0,12);
  fListHist->Add(fHistAsMcResPtVsPtLambda);

  fHistAsMcResPtVsPtAntiLambda       = new TH2F("h2AsMcResPtVsPtAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Pt",300,-0.15,0.15,240,0,12);
  fListHist->Add(fHistAsMcResPtVsPtAntiLambda);

  // pdgcode of mother
  fHistAsMcMotherPdgCodeK0s           = new TH1F("h1AsMcMotherPdgCodeK0s","Mother of Associated K^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeK0s);

  fHistAsMcMotherPdgCodeLambda        = new TH1F("h1AsMcMotherPdgCodeLambda","Mother of Associated #Lambda^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeLambda);

  fHistAsMcMotherPdgCodeAntiLambda    = new TH1F("h1AsMcMotherPdgCodeAntiLambda","Mother of Associated #bar{#Lambda}^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeAntiLambda);

  // Pt distribution Lambda from Sigma
  fHistAsMcPtLambdaFromSigma          = new TH1F("h1AsMcPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",240,0,12);
  fListHist->Add(fHistAsMcPtLambdaFromSigma);

  fHistAsMcPtAntiLambdaFromSigma      = new TH1F("h1AsMcPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",240,0,12);
  fListHist->Add(fHistAsMcPtAntiLambdaFromSigma);

  // Associated secondary particles:
  // Pt and rapidity distribution
  fHistAsMcSecondaryPtVsRapK0s          = new TH2F("h2AsMcSecondaryPtVsRapK0s", "K^{0} associated secondary;p_{t} (GeV/c);rapidity",240,0,12,30,-1.5,1.5);
  fListHist->Add(fHistAsMcSecondaryPtVsRapK0s);

  fHistAsMcSecondaryPtVsRapLambda       = new TH2F("h2AsMcSecondaryPtVsRapLambda", "#Lambda^{0} associated secondary;p_{t} (GeV/c);rapidity",240,0,12,30,-1.5,1.5);
  fListHist->Add(fHistAsMcSecondaryPtVsRapLambda);

  fHistAsMcSecondaryPtVsRapAntiLambda   = new TH2F("h2AsMcSecondaryPtVsRapAntiLambda", "#bar{#Lambda}^{0} associated secondary;p_{t} (GeV/c);rapidity",240,0,12,30,-1.5,1.5);
  fListHist->Add(fHistAsMcSecondaryPtVsRapAntiLambda);

  // Production radius
  fHistAsMcSecondaryProdRadiusK0s              = new TH1F("h1AsMcSecondaryProdRadiusK0s", "K^{0} Production Radius;r (cm);Count", 170, -2, 15);
  fListHist->Add(fHistAsMcSecondaryProdRadiusK0s);

  fHistAsMcSecondaryProdRadiusLambda           = new TH1F("h1AsMcSecondaryProdRadiusLambda", "#Lambda^{0} Production Radius;r (cm);Count", 170, -2, 15);
  fListHist->Add(fHistAsMcSecondaryProdRadiusLambda);

  fHistAsMcSecondaryProdRadiusAntiLambda       = new TH1F("h1AsMcSecondaryProdRadiusAntiLambda", "#bar{#Lambda}^{0} Production Radius;r (cm);Count", 170, -2, 15);
  fListHist->Add(fHistAsMcSecondaryProdRadiusAntiLambda);  

  fHistAsMcSecondaryProdRadiusXvsYK0s          = new TH2F("h2AsMcSecondaryProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYK0s);

  fHistAsMcSecondaryProdRadiusXvsYLambda       = new TH2F("h2AsMcSecondaryProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYLambda);

  fHistAsMcSecondaryProdRadiusXvsYAntiLambda   = new TH2F("h2AsMcSecondaryProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYAntiLambda);

  fHistAsMcSecondaryMotherPdgCodeK0s           = new TH1F("h1AsMcSecondaryMotherPdgCodeK0s","Mother of Associated Secondary K^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeK0s);

  fHistAsMcSecondaryMotherPdgCodeLambda        = new TH1F("h1AsMcSecondaryMotherPdgCodeLambda","Mother of Associated Secondary #Lambda^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeLambda);

  fHistAsMcSecondaryMotherPdgCodeAntiLambda    = new TH1F("h1AsMcSecondaryMotherPdgCodeAntiLambda","Mother of Associated Secondary #bar{#Lambda}^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeAntiLambda);

  // Pt distribution Lambda from Sigma
  fHistAsMcSecondaryPtLambdaFromSigma          = new TH1F("h1AsMcSecondaryPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",240,0,12);
  fListHist->Add(fHistAsMcSecondaryPtLambdaFromSigma);

  fHistAsMcSecondaryPtAntiLambdaFromSigma      = new TH1F("h1AsMcSecondaryPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",240,0,12);
  fListHist->Add(fHistAsMcSecondaryPtAntiLambdaFromSigma);
  PostData(1, fListHist);
  PostData(2, fCentrSelector);
  PostData(3, fTracksCuts);
}

//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  AliStack* stack = NULL;
  //  TClonesArray *mcArray = NULL;
  TArrayF mcPrimaryVtx;

  fESD=(AliESDEvent *)InputEvent();

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }

  // FIXME: levent not used. Can I remove it?
  AliVEvent* lEvent = InputEvent();

  if (!lEvent) {
    Printf("ERROR: Event not available");
    return;
  }
  // PID
        if (fUsePID.Contains("withPID")) {
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
	}


  //******************
  // Trigger Selection ! Warning Works only for ESD, add protection in case of AOD loop
  //******************

  Bool_t isSelected =
    (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()
     & AliVEvent::kMB);
  if (!isSelected) return;

  // Centrality selection
  static AliESDtrackCuts * trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1); // FIXME: make it a data member
  //   printf("%x, %x", fCentrSelector, trackCuts);
  Bool_t isCentralitySelected = fCentrSelector->IsCentralityBinSelected(fESD,trackCuts);
  if(!isCentralitySelected) return;
  // FIXME: add to hist number events another entry for centrality.

  // Done by the AliPhysicsSelection Task ! Only the selected events are passed to this task


  fHistNumberEvents->Fill(1.5);  // FIXME: use enum here

  //*************************
  //End track multiplicity
  //*************************

  // Remove Events with no tracks
  //if (!(fESD->GetNumberOfTracks()))  return;

  fHistNumberEvents->Fill(2.5);
  //  fHistTrackPerEvent->Fill(fESD->GetNumberOfTracks());

  //*************************************
  // Cut used:
  //*************************************
      
  // Cut Rapidity:
  Double_t lCutRap  = 0.5;

  // Cut AliKF Chi2 for Reconstructed particles
  //  Double_t cutChi2KF  = 1E3;

  // If PID is used:
  Double_t lLimitPPID    = 0.7;
  Float_t cutNSigmaLowP  = 1E3;
  Float_t cutNSigmaHighP = 1E3;
  if (fUsePID.Contains("withPID")) {
    cutNSigmaLowP  = 3.0;
    cutNSigmaHighP = 3.0;
  }


  // Cut Daughters pt (GeV/c):
  Double_t cutMinPtDaughter = 0.160;

  // Cut primary vertex:
  Double_t cutPrimVertex = 10.0;

  // Min number of TPC clusters:
  // Int_t nbMinTPCclusters = 80;

  
  //
  // PID flags:
  Int_t LambdaPID = 0;
  Int_t AntiLambdaPID = 0;
      

  //
  //  // Access MC:
  //

  if (fAnalysisMC) {
    if(fAnalysisType == "ESD") {
      AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if (!eventHandler) {
	Printf("ERROR: Could not retrieve MC event handler");
	return;
      }    
      AliMCEvent* mcEvent = eventHandler->MCEvent();
      if (!mcEvent) {
	Printf("ERROR: Could not retrieve MC event");
	return;
      }    
      stack = mcEvent->Stack();
      if (!stack) {
	Printf("ERROR: Could not retrieve stack");
	return;
      }
      
      AliGenEventHeader* mcHeader=mcEvent->GenEventHeader();
      if(!mcHeader) return;
      mcHeader->PrimaryVertex(mcPrimaryVtx);
      
    }
    
    /*
    // PID parameters for MC simulations:
    fAlephParameters[0] = 2.15898e+00/50.;
    fAlephParameters[1] = 1.75295e+01;
    fAlephParameters[2] = 3.40030e-09;
    fAlephParameters[3] = 1.96178e+00;
    fAlephParameters[4] = 3.91720e+00; 
    */
  }


  //**********************************************
  // MC loop
  //**********************************************

  Double_t lmcPrimVtxR      = 0;

  Int_t lNbMCPrimary        = 0;
  Int_t lNbMCPart           = 0;

  Int_t lPdgcodeCurrentPart = 0;
  Double_t lRapCurrentPart  = 0;
  Double_t lPtCurrentPart   = 0;
  
  Int_t lComeFromSigma      = 0;

  
  // Production Radius
  Double_t mcPosX     = 0.0,  mcPosY      = 0.0,  mcPosZ      = 0.0;
  Double_t mcPosR     = 0.0;

  // Decay Radius
  Double_t mcDecayPosX = 0, mcDecayPosY = 0, mcDecayPosR = 0;

  // current mc particle 's mother
  Int_t iCurrentMother  = 0, lPdgCurrentMother    = 0;
  //  Bool_t lCurrentMotherIsPrimary;

  // current mc particles 's daughter:
  Int_t lPdgCurrentDaughter0 = 0, lPdgCurrentDaughter1 = 0; 

  // variables for multiple reconstruction studies:
  Int_t id0           = 0, id1          = 0;
  //Int_t lLabelTrackN  = 0, lLabelTrackP = 0;
  //Int_t lPartNMother  = 0, lPartPMother = 0;
  //Int_t lPartPMotherPDGcode      = 0;
  Int_t lNtimesReconstructedK0s   = 0, lNtimesReconstructedLambda   = 0, lNtimesReconstructedAntiLambda   = 0;
  // Int_t lNtimesReconstructedK0sMI = 0, lNtimesReconstructedLambdaMI = 0, lNtimesReconstructedAntiLambdaMI = 0;

  //****************************
  // Start loop over MC particles
  if (fAnalysisMC) {

    // Primary vertex
    fHistMCPrimaryVertexX->Fill(mcPrimaryVtx.At(0));
    fHistMCPrimaryVertexY->Fill(mcPrimaryVtx.At(1));
    fHistMCPrimaryVertexZ->Fill(mcPrimaryVtx.At(2));
    
    lmcPrimVtxR = TMath::Sqrt(mcPrimaryVtx.At(0)*mcPrimaryVtx.At(0)+mcPrimaryVtx.At(1)*mcPrimaryVtx.At(1));
  

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
	
	// Keep only K0s, Lambda and AntiLambda, Xi and Phi:
	if ( (lPdgcodeCurrentPart != 310 ) && (lPdgcodeCurrentPart != 3122 ) && (lPdgcodeCurrentPart != -3122 ) && (lPdgcodeCurrentPart != 3312 ) && (lPdgcodeCurrentPart != -3312) && (lPdgcodeCurrentPart != -333) ) continue;
	
	lRapCurrentPart   = MyRapidity(p0->Energy(),p0->Pz());
	//lEtaCurrentPart   = p0->Eta();
	lPtCurrentPart    = p0->Pt();
	iCurrentMother    = p0->GetFirstMother();

	//	lPdgCurrentMother = stack->Particle(iCurrentMother)->GetPdgCode();
	if (iCurrentMother == -1){lPdgCurrentMother=0; } else {lPdgCurrentMother = stack->Particle(iCurrentMother)->GetPdgCode();}	

	mcPosX = p0->Vx();
	mcPosY = p0->Vy();
	mcPosZ = p0->Vz();
	mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
	
	id0  = p0->GetDaughter(0);
	id1  = p0->GetDaughter(1);

	// Decay Radius and Production Radius
	if ( id0 <= lNbMCPart && id0 > 0 && id1 <= lNbMCPart && id1 > 0) {
	  TParticle *pDaughter0 = stack->Particle(id0);
	  TParticle *pDaughter1 = stack->Particle(id1);
	  lPdgCurrentDaughter0 = pDaughter0->GetPdgCode();
	  lPdgCurrentDaughter1 = pDaughter1->GetPdgCode();
	  
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
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllK0s->Fill(lPtCurrentPart);
	}
	else if (lPdgcodeCurrentPart==3122)  {
	  fHistMCtracksProdRadiusLambda->Fill(mcPosX,mcPosY);
	  fHistMCtracksDecayRadiusLambda->Fill(mcDecayPosR);
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllLambda->Fill(lPtCurrentPart);
	}
	else if (lPdgcodeCurrentPart==-3122) {
	  fHistMCtracksProdRadiusAntiLambda->Fill(mcPosX,mcPosY);
	  fHistMCtracksDecayRadiusAntiLambda->Fill(mcDecayPosR);
	  if (TMath::Abs(lRapCurrentPart) < lCutRap) fHistMCPtAllAntiLambda->Fill(lPtCurrentPart);
	}
	
	if ( ( ( TMath::Abs(lPdgCurrentMother) == 3212)  ||
	       ( TMath::Abs(lPdgCurrentMother) == 3224)  ||
	       ( TMath::Abs(lPdgCurrentMother) == 3214)  ||
	       ( TMath::Abs(lPdgCurrentMother) == 3114) )
	     && ( iCurrentMother <= lNbMCPrimary )
	     ) lComeFromSigma = 1;
	else lComeFromSigma = 0;
	
	//*********************************************
	// Now keep only primary particles   
	if ( ( iMc > lNbMCPrimary ) && (!lComeFromSigma) ) continue;

	//********************************************
     
	lNtimesReconstructedK0s   = 0; lNtimesReconstructedLambda   = 0; lNtimesReconstructedAntiLambda   = 0;

 
        // Rap distribution
        if (lPdgcodeCurrentPart==310) {
	  fHistMCRapK0s->Fill(lRapCurrentPart);
	  if (lPtCurrentPart < 0.2 && lPtCurrentPart < 3.0)
	    fHistMCRapInPtRangeK0s->Fill(lRapCurrentPart);
	}

	if (lPdgcodeCurrentPart==3122) {
	  fHistMCRapLambda->Fill(lRapCurrentPart);
	  if (lPtCurrentPart < 0.6 && lPtCurrentPart < 3.5)
	    fHistMCRapInPtRangeLambda->Fill(lRapCurrentPart);
	}

	if (lPdgcodeCurrentPart==-3122) {
	  fHistMCRapAntiLambda->Fill(lRapCurrentPart);
	  if (lPtCurrentPart < 0.6 && lPtCurrentPart < 3.5)
	    fHistMCRapInPtRangeAntiLambda->Fill(lRapCurrentPart);
	}

	if (lPdgcodeCurrentPart==3312 || lPdgcodeCurrentPart==-3312) {
	  fHistMCRapXi->Fill(lRapCurrentPart);
	  if (lPtCurrentPart < 0.6 && lPtCurrentPart < 3.0)
	    fHistMCRapInPtRangeXi->Fill(lRapCurrentPart);
	}

	if (lPdgcodeCurrentPart==333) {
	  fHistMCRapPhi->Fill(lRapCurrentPart);
	  if (lPtCurrentPart < 0.7 && lPtCurrentPart < 3.0)
	    fHistMCRapInPtRangePhi->Fill(lRapCurrentPart);
	}
 
	// Rapidity Cut
	if (TMath::Abs(lRapCurrentPart) > lCutRap) continue;
 
	if (lPdgcodeCurrentPart==310) {
	  fHistMCProdRadiusK0s->Fill(mcPosR);

	  fHistMCPtK0s->Fill(lPtCurrentPart);



	  fHistNTimesRecK0s->Fill(lNtimesReconstructedK0s);
	  fHistNTimesRecK0sVsPt->Fill(lPtCurrentPart,lNtimesReconstructedK0s);
	}
	else 
	  if (lPdgcodeCurrentPart==3122) {
	    fHistMCProdRadiusLambda->Fill(mcPosR);

	    fHistMCPtLambda->Fill(lPtCurrentPart);	  


	    fHistNTimesRecLambda->Fill(lNtimesReconstructedLambda);
	    fHistNTimesRecLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedLambda);
	    if (lComeFromSigma) fHistMCPtLambdaFromSigma->Fill(lPtCurrentPart);

	    //printf("found Lambda MC pT=%e\n",lPtCurrentPart);
	    //printf("found Lambda MC Plabel=%d PPDGcode=%d Nlabel=%d NPDGcode=%d\n\n",id0,lPdgCurrentDaughter0,id1,lPdgCurrentDaughter1); 
	  
	  }
	else 
	  if (lPdgcodeCurrentPart==-3122) {
	    fHistMCProdRadiusAntiLambda->Fill(mcPosR);

	    fHistMCPtAntiLambda->Fill(lPtCurrentPart);	  


	    fHistNTimesRecAntiLambda->Fill(lNtimesReconstructedAntiLambda);
	    fHistNTimesRecAntiLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedAntiLambda);
	    if (lComeFromSigma) fHistMCPtAntiLambdaFromSigma->Fill(lPtCurrentPart);

	    //printf("found Lambda MC pT=%e\n",lPtCurrentPart);
	    //printf("found Lambda MC Plabel=%d PPDGcode=%d Nlabel=%d NPDGcode=%d\n\n",id0,lPdgCurrentDaughter0,id1,lPdgCurrentDaughter1); 
	  
	  }

	
      } // end loop ESD MC
      
    } // end ESD condition


  } // End Loop over MC condition

  



  //************************************
  // ESD loop 
  //************************************

  Double_t lMagneticField = 999;

  //Multiplcity:
  Int_t    nv0sTot= 0, nv0s = 0;
  //  Int_t nv0sMI =0;   
  // Variables:
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

  

  Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
  Double_t lPtK0s      = 0, lPtLambda      = 0, lPtAntiLambda      = 0;
  Double_t lRapK0s     = 0, lRapLambda     = 0, lRapAntiLambda     = 0;
  //  Double_t lEtaK0s     = 0, lEtaLambda     = 0, lEtaAntiLambda     = 0;
  Double_t lAlphaV0      = 0, lPtArmV0       = 0;

  Double_t lPzK0s      = 0, lPzLambda      = 0,  lPzAntiLambda      = 0;


  Double_t lV0Eta = 999;
  
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

  // Reconstructed position
  Double_t rcPosXK0s        = 0,  rcPosYK0s        = 0, rcPosZK0s        = 0;
  Double_t rcPosRK0s        = 0;
  Double_t rcPosXLambda     = 0,  rcPosYLambda     = 0, rcPosZLambda     = 0;
  Double_t rcPosRLambda     = 0;
    Double_t rcPosXAntiLambda = 0,  rcPosYAntiLambda = 0, rcPosZAntiLambda = 0;
   Double_t rcPosRAntiLambda = 0;

  // Pt resolution
   Double_t deltaPtK0s  = 0, deltaPtLambda  = 0, deltaPtAntiLambda  = 0;

  // Daughters
  AliESDtrack  *myTrackPos  = NULL;
  AliESDtrack  *myTrackNeg  = NULL;
  //  AliVParticle *lVPartPos   = NULL;
  // AliVParticle *lVPartNeg   = NULL;

  //  V0 momentum      
  Double_t V0mom[3] = {999,999,999};
  Double_t V0momentum = 0;

  // Daughters' momentum:
  Double_t  lMomPos[3] = {999,999,999};
  Double_t  lMomNeg[3] = {999,999,999};
  Double_t  lPtPos = 999, lPtNeg = 999;
  Double_t  lPPos = 999, lPNeg = 999;

  // Inner Wall parameters:
  Double_t  lMomInnerWallPos =999, lMomInnerWallNeg = 999;

  // AliKF Chi2 and Armenteros variables
  //  Double_t lChi2KFK0s  = 0, lChi2KFLambda = 0,  lChi2KFAntiLambda = 0;
  //  Double_t lAlphaV0K0s = 0, lAlphaV0Lambda = 0,  lAlphaV0AntiLambda = 0;
  //Double_t lPtArmV0K0s = 0, lPtArmV0Lambda = 0,  lPtArmV0AntiLambda = 0;
  //  Double_t lQlPos   = 0, lQlNeg   = 0;


  // PID
  Float_t nSigmaPosPion   = 0;
  Float_t nSigmaNegPion   = 0;

  Float_t nSigmaPosProton = 0;
  Float_t nSigmaNegProton = 0;
  

  Int_t lCheckPIDK0sPosDaughter        = 0, lCheckPIDK0sNegDaughter        = 0;
  Int_t lCheckPIDLambdaPosDaughter     = 0, lCheckPIDLambdaNegDaughter     = 0;
  Int_t lCheckPIDAntiLambdaPosDaughter = 0, lCheckPIDAntiLambdaNegDaughter = 0;

  
  
  //***********************
  // Primary Vertex cuts &
  // Magnetic field and Quality tracks cuts 

  Double_t  lPrimaryVtxPosition[3];
  Double_t  lPrimaryVtxCov[6];
  Double_t  lPrimaryVtxChi2 = 999;
  Double_t  lResPrimaryVtxX = 999;
  Double_t  lResPrimaryVtxY = 999;
  Double_t  lResPrimaryVtxZ = 999;
     
  AliAODVertex *myPrimaryVertex = NULL;
  //const AliVVertex *mySPDPrimaryVertex = NULL;

     
  const AliMultiplicity *myMultiplicty = ((AliESDEvent*)fESD)->GetMultiplicity();

  if(fAnalysisType == "ESD") {  
  
    // Best Primary Vertex:  
    const AliESDVertex *myBestPrimaryVertex = ((AliESDEvent*)fESD)->GetPrimaryVertex();
    myBestPrimaryVertex = ((AliESDEvent*)fESD)->GetPrimaryVertex();
    if (!myBestPrimaryVertex) return;
    if (!myBestPrimaryVertex->GetStatus()) return;
    fHistNumberEvents->Fill(3.5);
    myBestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);
    myBestPrimaryVertex->GetCovMatrix(lPrimaryVtxCov);
    if ( ( TMath::Abs(lPrimaryVtxPosition[2]) ) > cutPrimVertex) return ;
    fHistNumberEvents->Fill(4.5);    
    lPrimaryVtxChi2 = myBestPrimaryVertex->GetChi2toNDF();
    lResPrimaryVtxX = myBestPrimaryVertex->GetXRes();
    lResPrimaryVtxY = myBestPrimaryVertex->GetYRes();
    lResPrimaryVtxZ = myBestPrimaryVertex->GetZRes();

    // remove TPC-only primary vertex : retain only events with tracking + SPD vertex
    const AliESDVertex *mySPDPrimaryVertex = ((AliESDEvent*)fESD)->GetPrimaryVertexSPD();
    if (!mySPDPrimaryVertex) return;
    fHistSPDPrimaryVertexZ->Fill(mySPDPrimaryVertex->GetZ());
    const AliESDVertex *myPrimaryVertexTracking = ((AliESDEvent*)fESD)->GetPrimaryVertexTracks();
    if (!myPrimaryVertexTracking) return;
    if (!mySPDPrimaryVertex->GetStatus() && !myPrimaryVertexTracking->GetStatus() ) return;
    fHistNumberEvents->Fill(5.5);

    fHistTrackPerEvent->Fill(fESD->GetNumberOfTracks());   
    myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);
    if (!myPrimaryVertex) return;


    // Number of Tracklets:
    //const AliMultiplicity *myMultiplicty = ((AliESDEvent*)fESD)->GetMultiplicity();
    //if (myMultiplicty->GetNumberOfTracklets() < 10) return;
    fHistTrackletPerEvent->Fill(myMultiplicty->GetNumberOfTracklets());

    lMagneticField = ((AliESDEvent*)fESD)->GetMagneticField();

        fHistTPCTracks->Fill(AliESDtrackCuts::GetReferenceMultiplicity((AliESDEvent*)fESD, kTRUE));

    ////////////////////////////////////////////////////////////////////////////////////////
    Int_t i =0;


    for (Int_t jTracks=0;jTracks<fESD->GetNumberOfTracks();jTracks++){
		
      AliESDtrack* tPCtrack=fESD->GetTrack(jTracks);
      Float_t xy=0;
      Float_t z=0;
      tPCtrack->GetImpactParameters(xy,z);
      if ((fTracksCuts->IsSelected(tPCtrack))&&(xy<1.0)&&(z<1.0)) {i=i+1;}
			
    }

    fHistTPCMult->Fill(i);

    /////////////////////////////////////////////////////////////////////////////////////////

    // minimum number of clusters in TPC
    // fTracksCuts->SetMinNClustersTPC(nbMinTPCclusters);
  }
 
  fHistPrimaryVertexX->Fill(lPrimaryVtxPosition[0]);
  fHistPrimaryVertexY->Fill(lPrimaryVtxPosition[1]);
  fHistPrimaryVertexZ->Fill(lPrimaryVtxPosition[2]);
  //Double_t lrcPrimVtxR = TMath::Sqrt(lPrimaryVtxPosition[0]*lPrimaryVtxPosition[0]+lPrimaryVtxPosition[0]*lPrimaryVtxPosition[0]);

  fHistPrimaryVertexResX->Fill(lResPrimaryVtxX);
  fHistPrimaryVertexResY->Fill(lResPrimaryVtxY);
  fHistPrimaryVertexResZ->Fill(lResPrimaryVtxZ);

  //***********************
  // AliKF Primary Vertex

  AliKFVertex primaryVtxKF( *myPrimaryVertex );
  AliKFParticle::SetField(lMagneticField);




  //***Rerun the V0 finder

  //  fESD->ResetV0s();
  //  AliV0vertexer v0Vertexer;
  //  v0Vertexer.SetCuts(fCuts);
  //  v0Vertexer.Tracks2V0vertices(fESD);
  
  //*************************
  // V0 loop
      
  nv0sTot = fESD->GetNumberOfV0s();
  if (!nv0sTot) fHistNumberEvents->Fill(6.5);

  for (Int_t iV0 = 0; iV0 < nv0sTot; iV0++) {
    
    // ALiKF
    AliKFParticle* negPiKF = NULL;
    AliKFParticle* posPiKF = NULL;
    AliKFParticle* posPKF  = NULL;
    AliKFParticle* negAPKF = NULL;

    
    lIndexPosMother     = 0; lIndexNegMother     = 0; lIndexMotherOfMother       = 0;
    lCheckPIdK0Short    = 0; lCheckMcK0Short     = 0; lCheckSecondaryK0s         = 0;
    lCheckPIdLambda     = 0; lCheckMcLambda      = 0; lCheckSecondaryLambda      = 0;
    lCheckPIdAntiLambda = 0; lCheckMcAntiLambda  = 0; lCheckSecondaryAntiLambda  = 0;       
    lComeFromSigma      = -1;lCheckGamma = 0;
    
    
    if(fAnalysisType == "ESD") {


      AliESDv0 *v0 = ((AliESDEvent*)fESD)->GetV0(iV0);
      if (!v0) continue;

      //      if ((v0->Pt())<0.6) continue;
      
      // Primary vertex:
      fHistPrimaryVertexPosXV0events->Fill(lPrimaryVtxPosition[0]);
      fHistPrimaryVertexPosYV0events->Fill(lPrimaryVtxPosition[1]);
      fHistPrimaryVertexPosZV0events->Fill(lPrimaryVtxPosition[2]);
      
      // V0's Daughters
      lIndexTrackPos = TMath::Abs(v0->GetPindex());
      lIndexTrackNeg = TMath::Abs(v0->GetNindex());
      AliESDtrack *myTrackPosTest = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackPos);
      AliESDtrack *myTrackNegTest = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackNeg);
      if (!myTrackPosTest || !myTrackNegTest) {
	Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
	continue;
      }
      // Remove like-sign
      if ((Int_t)myTrackPosTest->GetSign() == (Int_t)myTrackNegTest->GetSign()){
	continue;
      } 
     
      // VO's main characteristics to check the reconstruction cuts
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


      if( myTrackPosTest->GetSign() ==1){
	
	myTrackPos = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackPos);
	myTrackNeg = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackNeg);

	// Daughters' momentum;
	v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
	v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
        
        // V0 momentum
        v0->GetXYZ(V0mom[0],V0mom[1],V0mom[2]);
        V0momentum = TMath::Sqrt(V0mom[0]*V0mom[0] +  V0mom[1]*V0mom[1] +  V0mom[2]*V0mom[2]);

	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;
	
	negPiKF = new AliKFParticle( *(v0->GetParamN()) ,-211);
	posPiKF = new AliKFParticle( *(v0->GetParamP()) ,211);
	posPKF  = new AliKFParticle( *(v0->GetParamP()) ,2212);
	negAPKF = new AliKFParticle( *(v0->GetParamN()) ,-2212);
	
      }
           
      if( myTrackPosTest->GetSign() ==-1){
	
	myTrackPos = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackNeg);
	myTrackNeg = ((AliESDEvent*)fESD)->GetTrack(lIndexTrackPos);

	// Daughters' momentum;
	v0->GetPPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
	v0->GetNPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);

	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;
	
	negPiKF = new AliKFParticle( *(v0->GetParamP()) ,-211);
	posPiKF = new AliKFParticle( *(v0->GetParamN()) ,211);
	posPKF  = new AliKFParticle( *(v0->GetParamN()) ,2212);
	negAPKF = new AliKFParticle( *(v0->GetParamP()) ,-2212);

      }

      lLabelTrackPos = (UInt_t)TMath::Abs(myTrackPos->GetLabel());
      lLabelTrackNeg = (UInt_t)TMath::Abs(myTrackNeg->GetLabel());

      // Daughters Pt and P:
      lPtPos = TMath::Sqrt(lMomPos[0]*lMomPos[0] + lMomPos[1]*lMomPos[1]);
      lPtNeg = TMath::Sqrt(lMomNeg[0]*lMomNeg[0] + lMomNeg[1]*lMomNeg[1]);

      lPPos = TMath::Sqrt(lMomPos[0]*lMomPos[0] + lMomPos[1]*lMomPos[1] + lMomPos[2]*lMomPos[2]);
      lPNeg = TMath::Sqrt(lMomNeg[0]*lMomNeg[0] + lMomNeg[1]*lMomNeg[1] + lMomNeg[2]*lMomNeg[2]);

      // Inner Wall parameter:
      const AliExternalTrackParam *myInnerWallTrackPos = myTrackPos->GetInnerParam(); 
      if(myInnerWallTrackPos) lMomInnerWallPos = myInnerWallTrackPos->GetP(); 
      const AliExternalTrackParam *myInnerWallTrackNeg = myTrackNeg->GetInnerParam(); 
      if(myInnerWallTrackNeg) lMomInnerWallNeg = myInnerWallTrackNeg->GetP(); 
	      
      // DCA between daughter and Primary Vertex:
      if (myTrackPos) lDcaPosToPrimVertex = TMath::Abs(myTrackPos->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
      
      if (myTrackNeg) lDcaNegToPrimVertex = TMath::Abs(myTrackNeg->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
      
      // Quality tracks cuts:
      if ( !(fTracksCuts->IsSelected(myTrackPos)) || !(fTracksCuts->IsSelected(myTrackNeg)) ) 


	{
	  if (negPiKF) delete negPiKF; negPiKF=NULL;
	  if (posPiKF) delete posPiKF; posPiKF=NULL;
	  if (posPKF)  delete posPKF;  posPKF=NULL;
	  if (negAPKF) delete negAPKF; negAPKF=NULL;

	  continue;}

      // Armenteros variables:
      lAlphaV0      =  v0->AlphaV0();
      lPtArmV0      =  v0->PtArmV0();

      // Pseudorapidity:
      lV0Eta = v0->Eta();
      //////////////////////////////////////////////////////////////////////////
      // Invariant mass
      v0->ChangeMassHypothesis(310);
      lInvMassK0s = v0->GetEffMass();
      lPtK0s = v0->Pt();
      lPzK0s = v0->Pz();

      v0->ChangeMassHypothesis(3122);
      lInvMassLambda = v0->GetEffMass();
      lPtLambda = v0->Pt();
      lPzLambda = v0->Pz();

      v0->ChangeMassHypothesis(-3122);
      lInvMassAntiLambda = v0->GetEffMass();
      lPtAntiLambda = v0->Pt();
      lPzAntiLambda = v0->Pz();
  

      // Rapidity:
      lRapK0s    = v0->Y(310);
      lRapLambda = v0->Y(3122);
      lRapAntiLambda = v0->Y(-3122);
	
      if (lPtK0s==0) 	{
	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;

	continue;}
      if (lPtLambda==0) 	{
	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;

	continue;}

      if (lPtAntiLambda==0) 	{
	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;

	continue;}

      ///////////////////////////////////////////////////////////////////////      

      // PID  new method July 2011
      if (fUsePID.Contains("withPID")) {
	//	nSigmaPosPion   = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackPos,AliPID::kPion));
	nSigmaPosPion =	TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion));
		//	nSigmaNegPion   = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackNeg,AliPID::kPion));
	nSigmaNegPion =	TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kPion));
	//	nSigmaPosProton = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackPos,AliPID::kProton));
	nSigmaPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton));
	//	nSigmaNegProton = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackNeg,AliPID::kProton));
        nSigmaNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kProton));
	
      }
      else {
	nSigmaPosPion = 0; nSigmaNegPion =0; nSigmaPosProton = 0; nSigmaNegProton= 0;
      }
      
      
      
      // Monte-Carlo particle associated to reconstructed particles: 
      if (fAnalysisMC) {
	//if (lLabelTrackPos < 0 || lLabelTrackNeg < 0) continue;
	TParticle  *lMCESDPartPos  = stack->Particle(lLabelTrackPos);
	if(!lMCESDPartPos) { 
	  //  Printf("no MC particle for positive and/or negative daughter\n");
	  	
	  if (negPiKF) delete negPiKF; negPiKF=NULL;
	  if (posPiKF) delete posPiKF; posPiKF=NULL;
	  if (posPKF)  delete posPKF;  posPKF=NULL;
	  if (negAPKF) delete negAPKF; negAPKF=NULL;

	  continue;
	}
	TParticle  *lMCESDPartNeg  = stack->Particle(lLabelTrackNeg);
	if (!lMCESDPartNeg) 	{
	  if (negPiKF) delete negPiKF; negPiKF=NULL;
	  if (posPiKF) delete posPiKF; posPiKF=NULL;
	  if (posPKF)  delete posPKF;  posPKF=NULL;
	  if (negAPKF) delete negAPKF; negAPKF=NULL;

	  continue;}
	lPDGCodePosDaughter = lMCESDPartPos->GetPdgCode();
	lPDGCodeNegDaughter = lMCESDPartNeg->GetPdgCode();
	lIndexPosMother = lMCESDPartPos->GetFirstMother();
	lIndexNegMother = lMCESDPartNeg->GetFirstMother();

	//////////////////////////////
	//	if (lIndexPosMother == -1) 	{
	//	if (negPiKF) delete negPiKF; negPiKF=NULL;
	//	if (posPiKF) delete posPiKF; posPiKF=NULL;
	//	if (posPKF)  delete posPKF;  posPKF=NULL;
	//	if (negAPKF) delete negAPKF; negAPKF=NULL;

	//	continue;}

	if (lIndexPosMother == -1) {


	  lPdgcodeMother = 0;
	  lIndexMotherOfMother = 0;
	  mcPosX = 0;
	  mcPosY = 0;
	  mcPosZ = 0;
	  mcPosR = 0;
	  mcPosMotherX = 0;
	  mcPosMotherY = 0;
	  mcPosMotherZ = 0;
	  mcPosMotherR = 0;
	  mcMotherPt = 1;
	}

	else {



	  TParticle  *lMCESDMother    = stack->Particle(lIndexPosMother);
	  if (!lMCESDMother) 	{
	    if (negPiKF) delete negPiKF; negPiKF=NULL;
	    if (posPiKF) delete posPiKF; posPiKF=NULL;
	    if (posPKF)  delete posPKF;  posPKF=NULL;
	    if (negAPKF) delete negAPKF; negAPKF=NULL;

	    continue;}
	  lPdgcodeMother         = lMCESDMother->GetPdgCode();
	  lIndexMotherOfMother   = lMCESDMother->GetFirstMother();
	  if (lIndexMotherOfMother ==-1) lPdgcodeMotherOfMother = 0;
	  else {
	    TParticle  *lMCESDMotherOfMother    = stack->Particle(lIndexMotherOfMother);
	    if (!lMCESDMotherOfMother) 	{
	      if (negPiKF) delete negPiKF; negPiKF=NULL;
	      if (posPiKF) delete posPiKF; posPiKF=NULL;
	      if (posPKF)  delete posPKF;  posPKF=NULL;
	      if (negAPKF) delete negAPKF; negAPKF=NULL;

	      continue;}
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
	}
      }
    } // end ESD condition


    

    
    
    // Multiplicity:
    if(!lOnFlyStatus) nv0s++;
    //    else  if(lOnFlyStatus) nv0sMI++;

    // Daughter momentum cut: ! FIX it in case of AOD !
    if ( (lPtPos  < cutMinPtDaughter ) ||
         (lPtNeg  < cutMinPtDaughter )
	 ) 	{
      if (negPiKF) delete negPiKF; negPiKF=NULL;
      if (posPiKF) delete posPiKF; posPiKF=NULL;
      if (posPKF)  delete posPKF;  posPKF=NULL;
      if (negAPKF) delete negAPKF; negAPKF=NULL;

      continue;}
    
    AliKFParticle v0K0sKF;
    v0K0sKF+=(*negPiKF);
    v0K0sKF+=(*posPiKF);
    v0K0sKF.SetProductionVertex(primaryVtxKF);
    
    AliKFParticle v0LambdaKF;
    v0LambdaKF+=(*negPiKF);
    v0LambdaKF+=(*posPKF);	
    v0LambdaKF.SetProductionVertex(primaryVtxKF);
    
    AliKFParticle v0AntiLambdaKF;
    v0AntiLambdaKF+=(*posPiKF);
    v0AntiLambdaKF+=(*negAPKF);
    v0AntiLambdaKF.SetProductionVertex(primaryVtxKF);
    // Look for associated particles:
    if (fAnalysisMC) {
      if( (lIndexPosMother==-1) || (lIndexNegMother==-1) ) {
	fHistMCDaughterTrack->Fill(1);
      }
      
      else if( ( (lPDGCodePosDaughter==+211) && (lPDGCodeNegDaughter==-211) )    
	       ) {
	lCheckPIdK0Short    = 1;
	fHistMCDaughterTrack->Fill(3);
	if ( (lIndexPosMother==lIndexNegMother) &&
	     (lPdgcodeMother==310) ) {
	  if (lIndexPosMother <= lNbMCPrimary) lCheckMcK0Short  = 1;
	  else lCheckSecondaryK0s = 1;
	}
      }
      else if( ( (lPDGCodePosDaughter==+2212) && (lPDGCodeNegDaughter==-211)  )  
	       ) {
	lCheckPIdLambda     = 1;
	fHistMCDaughterTrack->Fill(5);
	if ( (lIndexPosMother==lIndexNegMother) &&
	     (lPdgcodeMother==3122)  ){
	  if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)
	       ) lComeFromSigma = 1;
	  else lComeFromSigma = 0; 
	  if ( (lIndexPosMother <= lNbMCPrimary) || 
	       ( ( lIndexPosMother > lNbMCPrimary) && (lComeFromSigma) )
	       ) lCheckMcLambda  = 1; 
	  else lCheckSecondaryLambda    = 1;
	}
      }
      else if( ( (lPDGCodePosDaughter==211)   && (lPDGCodeNegDaughter==-2212) )	     
	       ) {
	lCheckPIdAntiLambda = 1;
	fHistMCDaughterTrack->Fill(7);
	if ( (lIndexPosMother==lIndexNegMother) &&
	     (lPdgcodeMother==-3122) ) {
	  if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3224) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3214) ||
	       ( TMath::Abs(lPdgcodeMotherOfMother) == 3114)
	       ) lComeFromSigma = 1;
	  else lComeFromSigma = 0;  
	  if ( (lIndexPosMother <= lNbMCPrimary) || 
	       ( ( lIndexPosMother > lNbMCPrimary) && (lComeFromSigma) )
	       ) lCheckMcAntiLambda  = 1;
	  else lCheckSecondaryAntiLambda = 1;
	}
      }
      
      // Gamma conversion
      else if ( (lPDGCodePosDaughter==-11) &&
		(lPDGCodeNegDaughter==11) &&
		(lPdgcodeMother==22 ) )
	lCheckGamma = 1;
    } // end "look for associated particles  
   
    
    // PID condition:
    lCheckPIDK0sPosDaughter        = 0, lCheckPIDK0sNegDaughter        = 0;
    lCheckPIDLambdaPosDaughter     = 0, lCheckPIDLambdaNegDaughter     = 0;
    lCheckPIDAntiLambdaPosDaughter = 0, lCheckPIDAntiLambdaNegDaughter = 0;

    if (lMomInnerWallPos < lLimitPPID) {
      if (nSigmaPosPion < cutNSigmaLowP)   {
	lCheckPIDK0sPosDaughter        = 1;
	lCheckPIDAntiLambdaPosDaughter = 1;
      }
      if (nSigmaPosProton < cutNSigmaLowP) lCheckPIDLambdaPosDaughter    = 1;      
    }

    else if (lMomInnerWallPos > lLimitPPID) {
      if (nSigmaPosPion < cutNSigmaHighP)   {
	lCheckPIDK0sPosDaughter        = 1;
	lCheckPIDAntiLambdaPosDaughter = 1;
      }
      if (nSigmaPosProton < cutNSigmaHighP) lCheckPIDLambdaPosDaughter    = 1;
    }

    if (lMomInnerWallNeg < lLimitPPID) {
      if (nSigmaNegPion < cutNSigmaLowP)    {
	lCheckPIDK0sNegDaughter       = 1;
	lCheckPIDLambdaNegDaughter    = 1;
      }
      if (nSigmaNegProton < cutNSigmaLowP)  lCheckPIDAntiLambdaNegDaughter = 1;
      
    }
    else if (lMomInnerWallNeg > lLimitPPID) {
      if (nSigmaNegPion < cutNSigmaHighP)   {
	lCheckPIDK0sNegDaughter       = 1;
	lCheckPIDLambdaNegDaughter    = 1;
      }
      if (nSigmaNegProton < cutNSigmaHighP) lCheckPIDAntiLambdaNegDaughter = 1;
    }
 

    /*    
    if ( lCheckPIDAntiLambdaNegDaughter==1 || lCheckPIDAntiLambdaPosDaughter==1 || lCheckPIDLambdaPosDaughter!=1 || lCheckPIDLambdaNegDaughter!=1)
      {
      	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;

	continue;
      }     
    */
    //*****************************
    // filling histograms
    //*****************************
    /*
    ///////////////values for cuts/////////////////////////////////////////////////////////////////////////////////////////
    if ((lDcaPosToPrimVertex<=0.1) || (lDcaNegToPrimVertex<=0.1) || (lDcaV0Daughters>=1.00) || 
	(lV0cosPointAngle<=0.998) || (lV0Radius<=0.9) || (lV0Radius>=100) ) 
	
      {
	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;

	continue;}
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    */
    
    if(fUsePID.Contains("withPID") && (lCheckPIDAntiLambdaNegDaughter==0) && lCheckPIDLambdaPosDaughter==1) LambdaPID = 1;
    else LambdaPID =0;
    if(fUsePID.Contains("withPID") && (lCheckPIDLambdaPosDaughter==0) && lCheckPIDAntiLambdaNegDaughter==1) AntiLambdaPID = 1;
    else AntiLambdaPID =0;

    V0momentum = TMath::Sqrt(lPzLambda*lPzLambda + lPtLambda*lPtLambda);

    if ((LambdaPID==1 && V0momentum <=1) || (V0momentum>1) ||  !(fUsePID.Contains("withPID"))){  
      if ((TMath::Abs(lRapK0s) < lCutRap) && lOnFlyStatus==0) {
      fHistTPCsignalPt->Fill(V0momentum,myTrackPos->GetTPCsignal());
      // fHistTPCsignalPt->Fill(V0momentum,myTrackNeg->GetTPCsignal());
      }
    }

    if ((V0momentum)<1 && lOnFlyStatus==0 ){
        fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0);
	  }





    // insert PID condition for K0s instead of "kTRUE" value
    if ((fUsePID.Contains("withPID") && kTRUE) || !(fUsePID.Contains("withPID"))){  
      if (TMath::Abs(lRapK0s) < lCutRap ) {

	//////2D histos: cut vs on fly status/////////////////////

	fHistDcaPosToPrimVertexK0->Fill(lDcaPosToPrimVertex,lOnFlyStatus);
	fHistDcaNegToPrimVertexK0->Fill(lDcaNegToPrimVertex,lOnFlyStatus);
	fHistRadiusV0K0->Fill(lV0Radius,lOnFlyStatus);
	fHistDecayLengthV0K0->Fill(lV0DecayLength,lOnFlyStatus);
	fHistDcaV0DaughtersK0->Fill(lDcaV0Daughters,lOnFlyStatus);
	fHistChi2K0->Fill(lChi2V0,lOnFlyStatus);
	fHistCosPointAngleK0->Fill(lV0cosPointAngle,lOnFlyStatus);

	//////2D histos: cut vs mass///////////////////// 

	if (lOnFlyStatus==0){
        

   	fHistMassK0->Fill(lInvMassK0s);
	fHistMassVsRadiusK0->Fill(rcPosRK0s,lInvMassK0s);
	fHistPtVsMassK0->Fill(lInvMassK0s,lPtK0s);


	  fHistDcaPosToPrimVertexK0vsMassK0->Fill(lDcaPosToPrimVertex,lInvMassK0s);
	  fHistDcaNegToPrimVertexK0vsMassK0->Fill(lDcaNegToPrimVertex,lInvMassK0s);
	  fHistRadiusV0K0vsMassK0->Fill(lV0Radius,lInvMassK0s);
	  fHistDecayLengthV0K0vsMassK0->Fill(lV0DecayLength,lInvMassK0s);
	  fHistDcaV0DaughtersK0vsMassK0->Fill(lDcaV0Daughters,lInvMassK0s);
	  fHistCosPointAngleK0vsMassK0->Fill(lV0cosPointAngle,lInvMassK0s);

	  if (lPtK0s>0 && lPtK0s <3){ 
	    fHistDcaPosToPrimVertexK0vsMassK0pt1->Fill(lDcaPosToPrimVertex,lInvMassK0s);
	    fHistDcaNegToPrimVertexK0vsMassK0pt1->Fill(lDcaNegToPrimVertex,lInvMassK0s);
	    fHistRadiusV0K0vsMassK0pt1->Fill(lV0Radius,lInvMassK0s);
	    fHistDecayLengthV0K0vsMassK0pt1->Fill(lV0DecayLength,lInvMassK0s);
	    fHistDcaV0DaughtersK0vsMassK0pt1->Fill(lDcaV0Daughters,lInvMassK0s);
	    fHistCosPointAngleK0vsMassK0pt1->Fill(lV0cosPointAngle,lInvMassK0s);
	  }   
	  if (lPtK0s > 3 && lPtK0s < 6){ 
	    fHistDcaPosToPrimVertexK0vsMassK0pt2->Fill(lDcaPosToPrimVertex,lInvMassK0s);
	    fHistDcaNegToPrimVertexK0vsMassK0pt2->Fill(lDcaNegToPrimVertex,lInvMassK0s);
	    fHistRadiusV0K0vsMassK0pt2->Fill(lV0Radius,lInvMassK0s);
	    fHistDecayLengthV0K0vsMassK0pt2->Fill(lV0DecayLength,lInvMassK0s);
	    fHistDcaV0DaughtersK0vsMassK0pt2->Fill(lDcaV0Daughters,lInvMassK0s);
	    fHistCosPointAngleK0vsMassK0pt2->Fill(lV0cosPointAngle,lInvMassK0s);
	  }   
	  if (lPtK0s > 6 && lPtK0s < 10){ 
	    fHistDcaPosToPrimVertexK0vsMassK0pt3->Fill(lDcaPosToPrimVertex,lInvMassK0s);
	    fHistDcaNegToPrimVertexK0vsMassK0pt3->Fill(lDcaNegToPrimVertex,lInvMassK0s);
	    fHistRadiusV0K0vsMassK0pt3->Fill(lV0Radius,lInvMassK0s);
	    fHistDecayLengthV0K0vsMassK0pt3->Fill(lV0DecayLength,lInvMassK0s);
	    fHistDcaV0DaughtersK0vsMassK0pt3->Fill(lDcaV0Daughters,lInvMassK0s);
	    fHistCosPointAngleK0vsMassK0pt3->Fill(lV0cosPointAngle,lInvMassK0s);
	  }   
	}
      } // if rap. condition
    } // end if withPID condition
    
    // insert PID condition for Lambda instead of "kTRUE" value
    //    if ((fUsePID.Contains("withPID") && kTRUE )|| !(fUsePID.Contains("withPID"))){  
    if ((LambdaPID==1 && V0momentum <=1) || (V0momentum>1) ||  !(fUsePID.Contains("withPID"))){  

      if (TMath::Abs(lRapLambda) < lCutRap) {

	//////2D histos: cut vs on fly status/////////////////////

	fHistDcaPosToPrimVertexL->Fill(lDcaPosToPrimVertex,lOnFlyStatus);
	fHistDcaNegToPrimVertexL->Fill(lDcaNegToPrimVertex,lOnFlyStatus);
	fHistRadiusV0L->Fill(lV0Radius,lOnFlyStatus);
	fHistDecayLengthV0L->Fill(lV0DecayLength,lOnFlyStatus);
	fHistDcaV0DaughtersL->Fill(lDcaV0Daughters,lOnFlyStatus);
	fHistChi2L->Fill(lChi2V0,lOnFlyStatus);
	fHistCosPointAngleL->Fill(lV0cosPointAngle,lOnFlyStatus);

	//////2D histos: cut vs mass/////////////////////

	if (lOnFlyStatus==0){
	fHistMassLambda->Fill(lInvMassLambda);
	fHistMassVsRadiusLambda->Fill(rcPosRLambda,lInvMassLambda);
	fHistPtVsMassLambda->Fill(lInvMassLambda,lPtLambda);
        if(V0momentum <=1) fHistNSigmaProton->Fill(nSigmaPosProton);

	  fHistDcaPosToPrimVertexLvsMassL->Fill(lDcaPosToPrimVertex,lInvMassLambda);
	  fHistDcaNegToPrimVertexLvsMassL->Fill(lDcaNegToPrimVertex,lInvMassLambda);
	  fHistRadiusV0LvsMassL->Fill(lV0Radius,lInvMassLambda);
	  fHistDecayLengthV0LvsMassL->Fill(lV0DecayLength,lInvMassLambda);
	  fHistDcaV0DaughtersLvsMassL->Fill(lDcaV0Daughters,lInvMassLambda);
	  fHistCosPointAngleLvsMassL->Fill(lV0cosPointAngle,lInvMassLambda);

            
	  if (lPtLambda>0 && lPtLambda <3){ 
	    fHistDcaPosToPrimVertexLambdaVsMasspt1->Fill(lDcaPosToPrimVertex,lInvMassLambda);
	    fHistDcaNegToPrimVertexLambdaVsMasspt1->Fill(lDcaNegToPrimVertex,lInvMassLambda);
	    fHistRadiusV0LambdaVsMasspt1->Fill(lV0Radius,lInvMassLambda);
	    fHistDecayLengthV0LambdaVsMasspt1->Fill(lV0DecayLength,lInvMassLambda);
	    fHistDcaV0DaughtersLambdaVsMasspt1->Fill(lDcaV0Daughters,lInvMassLambda);
	    fHistCosPointAngleLambdaVsMasspt1->Fill(lV0cosPointAngle,lInvMassLambda);
	  }   
	  if (lPtLambda > 3 && lPtLambda < 6){ 
	    fHistDcaPosToPrimVertexLambdaVsMasspt2->Fill(lDcaPosToPrimVertex,lInvMassLambda);
	    fHistDcaNegToPrimVertexLambdaVsMasspt2->Fill(lDcaNegToPrimVertex,lInvMassLambda);
	    fHistRadiusV0LambdaVsMasspt2->Fill(lV0Radius,lInvMassLambda);
	    fHistDecayLengthV0LambdaVsMasspt2->Fill(lV0DecayLength,lInvMassLambda);
	    fHistDcaV0DaughtersLambdaVsMasspt2->Fill(lDcaV0Daughters,lInvMassLambda);
	    fHistCosPointAngleLambdaVsMasspt2->Fill(lV0cosPointAngle,lInvMassLambda);
	  }   
	  if (lPtLambda > 6 && lPtLambda < 10){ 
	    fHistDcaPosToPrimVertexLambdaVsMasspt3->Fill(lDcaPosToPrimVertex,lInvMassLambda);
	    fHistDcaNegToPrimVertexLambdaVsMasspt3->Fill(lDcaNegToPrimVertex,lInvMassLambda);
	    fHistRadiusV0LambdaVsMasspt3->Fill(lV0Radius,lInvMassLambda);
	    fHistDecayLengthV0LambdaVsMasspt3->Fill(lV0DecayLength,lInvMassLambda);
	    fHistDcaV0DaughtersLambdaVsMasspt3->Fill(lDcaV0Daughters,lInvMassLambda);
	    fHistCosPointAngleLambdaVsMasspt3->Fill(lV0cosPointAngle,lInvMassLambda);
	  }   
	}
      } //end of Rap condition
    } // end of PID condition

    /////////// Anti Lambda ///////////////
    if ((AntiLambdaPID==1 && lPtAntiLambda <=1) || (lPtAntiLambda>1) ||  !(fUsePID.Contains("withPID"))){  
      if (TMath::Abs(lRapAntiLambda) < lCutRap) {

	//////2D histos: cut vs on fly status/////////////////////

	fHistDcaPosToPrimVertexAntiL->Fill(lDcaPosToPrimVertex,lOnFlyStatus);
	fHistDcaNegToPrimVertexAntiL->Fill(lDcaNegToPrimVertex,lOnFlyStatus);
	fHistRadiusV0AntiL->Fill(lV0Radius,lOnFlyStatus);
	fHistDecayLengthV0AntiL->Fill(lV0DecayLength,lOnFlyStatus);
	fHistDcaV0DaughtersAntiL->Fill(lDcaV0Daughters,lOnFlyStatus);
	fHistChi2AntiL->Fill(lChi2V0,lOnFlyStatus);
	fHistCosPointAngleAntiL->Fill(lV0cosPointAngle,lOnFlyStatus);

	//////2D histos: cut vs mass/////////////////////

	if (lOnFlyStatus==0){

		fHistMassAntiLambda->Fill(lInvMassAntiLambda);
		fHistMassVsRadiusAntiLambda->Fill(rcPosRAntiLambda,lInvMassAntiLambda);
		fHistPtVsMassAntiLambda->Fill(lInvMassAntiLambda,lPtAntiLambda);


	  fHistDcaPosToPrimVertexAntiLvsMass->Fill(lDcaPosToPrimVertex,lInvMassAntiLambda);
	  fHistDcaNegToPrimVertexAntiLvsMass->Fill(lDcaNegToPrimVertex,lInvMassAntiLambda);
	  fHistRadiusV0AntiLvsMass->Fill(lV0Radius,lInvMassAntiLambda);
	  fHistDecayLengthV0AntiLvsMass->Fill(lV0DecayLength,lInvMassAntiLambda);
	  fHistDcaV0DaughtersAntiLvsMass->Fill(lDcaV0Daughters,lInvMassAntiLambda);
	  fHistCosPointAngleAntiLvsMass->Fill(lV0cosPointAngle,lInvMassAntiLambda);

            
	  if (lPtAntiLambda>0 && lPtAntiLambda <3){ 
	    fHistDcaPosToPrimVertexAntiLVsMasspt1->Fill(lDcaPosToPrimVertex,lInvMassAntiLambda);
	    fHistDcaNegToPrimVertexAntiLVsMasspt1->Fill(lDcaNegToPrimVertex,lInvMassAntiLambda);
	    fHistRadiusV0AntiLVsMasspt1->Fill(lV0Radius,lInvMassAntiLambda);
	    fHistDecayLengthV0AntiLVsMasspt1->Fill(lV0DecayLength,lInvMassAntiLambda);
	    fHistDcaV0DaughtersAntiLVsMasspt1->Fill(lDcaV0Daughters,lInvMassAntiLambda);
	    fHistCosPointAngleAntiLVsMasspt1->Fill(lV0cosPointAngle,lInvMassAntiLambda);
	  }   
	  if (lPtAntiLambda > 3 && lPtAntiLambda < 6){ 
	    fHistDcaPosToPrimVertexAntiLVsMasspt2->Fill(lDcaPosToPrimVertex,lInvMassAntiLambda);
	    fHistDcaNegToPrimVertexAntiLVsMasspt2->Fill(lDcaNegToPrimVertex,lInvMassAntiLambda);
	    fHistRadiusV0AntiLVsMasspt2->Fill(lV0Radius,lInvMassAntiLambda);
	    fHistDecayLengthV0AntiLVsMasspt2->Fill(lV0DecayLength,lInvMassAntiLambda);
	    fHistDcaV0DaughtersAntiLVsMasspt2->Fill(lDcaV0Daughters,lInvMassAntiLambda);
	    fHistCosPointAngleAntiLVsMasspt2->Fill(lV0cosPointAngle,lInvMassAntiLambda);
	  }   
	  if (lPtAntiLambda > 6 && lPtAntiLambda < 10){ 
	    fHistDcaPosToPrimVertexAntiLVsMasspt3->Fill(lDcaPosToPrimVertex,lInvMassAntiLambda);
	    fHistDcaNegToPrimVertexAntiLVsMasspt3->Fill(lDcaNegToPrimVertex,lInvMassAntiLambda);
	    fHistRadiusV0AntiLVsMasspt3->Fill(lV0Radius,lInvMassAntiLambda);
	    fHistDecayLengthV0AntiLVsMasspt3->Fill(lV0DecayLength,lInvMassAntiLambda);
	    fHistDcaV0DaughtersAntiLVsMasspt3->Fill(lDcaV0Daughters,lInvMassAntiLambda);
	    fHistCosPointAngleAntiLVsMasspt3->Fill(lV0cosPointAngle,lInvMassAntiLambda);
	  }   
	}
      } //end of Rap condition
    } // end of PID condition


    ///////////////values for cuts end////////////////////////////////////////////////////////////////////////


    // Histo versus Rap and armenteros plot
    if (!lOnFlyStatus){
      if (lCheckMcK0Short) fHistAsMcRapK0->Fill(lRapK0s);
      if (lCheckMcLambda) fHistAsMcRapLambda->Fill(lRapLambda);
      if (lCheckMcAntiLambda) fHistAsMcRapLambda->Fill(lRapAntiLambda);
      //      fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0);
      if ((TMath::Abs(lRapK0s) < lCutRap)&&(TMath::Abs(lRapLambda) < lCutRap)) fHistK0sMassVsLambdaMass->Fill(lInvMassK0s,lInvMassLambda);
    }


    
    // K0s associated histograms in |rap| < lCutRap:


    if (TMath::Abs(lRapK0s) < lCutRap) {

      fHistNsigmaPosPionK0->Fill(nSigmaPosPion);
      fHistNsigmaNegPionK0->Fill(nSigmaNegPion);
	
      switch (lOnFlyStatus){
      case 0 : 


	//	  fHistMultVsPtVsMassK0->Fill(multiplicity ,lInvMassK0s,lPtK0s);
	if(lCheckPIdK0Short) fHistPidMcMassK0->Fill(lInvMassK0s);
	if(lCheckMcK0Short) {
	  fHistAsMcMassK0->Fill(lInvMassK0s);
	  fHistAsMcPtK0->Fill(lPtK0s);


	  fHistAsMcPtVsMassK0->Fill(lInvMassK0s,lPtK0s);
	  if (lPtK0s <= 1) fHistAsMcPtZoomK0->Fill(lPtK0s);
	  fHistAsMcMassVsRadiusK0->Fill(rcPosRK0s,lInvMassK0s);
	  fHistAsMcResxK0->Fill(rcPosXK0s-mcPosX);
	  fHistAsMcResyK0->Fill(rcPosYK0s-mcPosY);
	  fHistAsMcReszK0->Fill(rcPosZK0s-mcPosZ);
	  fHistAsMcResrVsRadiusK0->Fill(rcPosRK0s,rcPosRK0s-mcPosR);
	  fHistAsMcReszVsRadiusK0->Fill(rcPosZK0s,rcPosZK0s-mcPosZ);
	  fHistAsMcProdRadiusK0->Fill(mcPosMotherR);
	  fHistAsMcProdRadiusXvsYK0s->Fill(mcPosMotherX,mcPosMotherY);
	  fHistAsMcResPtK0->Fill(deltaPtK0s);
	  fHistAsMcResPtVsRapK0->Fill(deltaPtK0s,lRapK0s);
	  fHistAsMcResPtVsPtK0->Fill(deltaPtK0s,lPtK0s);
	}
	else if (lCheckSecondaryK0s) {
	  fHistAsMcSecondaryPtVsRapK0s->Fill(lPtK0s,lRapK0s);
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
	  
      }
    } // end rapidity condition

    

    // Associated Lambda histograms in |rap| < lCutRap
    if ((LambdaPID==1 && V0momentum <=1) || (V0momentum>1) ||  !(fUsePID.Contains("withPID"))){  

    if (TMath::Abs(lRapLambda) < lCutRap) {

      fHistNsigmaPosProtonLambda->Fill(nSigmaPosProton);
      fHistNsigmaNegPionLambda->Fill(nSigmaNegPion);
      switch (lOnFlyStatus){
      case 0 : 

	//          fHistMultVsPtVsMassLambda->Fill(multiplicity ,lInvMassLambda,lPtLambda);
	if(lCheckPIdLambda) fHistPidMcMassLambda->Fill(lInvMassLambda);
	  
	if(lCheckMcLambda) {
	  fHistAsMcMassLambda->Fill(lInvMassLambda);
	  fHistAsMcPtLambda->Fill(lPtLambda);

	  fHistAsMcPtVsMassLambda->Fill(lInvMassLambda,lPtLambda);
	  if (lPtLambda <= 1) fHistAsMcPtZoomLambda->Fill(lPtLambda);
	  fHistAsMcMassVsRadiusLambda->Fill(rcPosRLambda,lInvMassLambda);
	  fHistAsMcResxLambda->Fill(rcPosXLambda-mcPosX);
	  fHistAsMcResyLambda->Fill(rcPosYLambda-mcPosY);
	  fHistAsMcReszLambda->Fill(rcPosZLambda-mcPosZ);
	  fHistAsMcResrVsRadiusLambda->Fill(rcPosRLambda,rcPosRLambda-mcPosR);
	  fHistAsMcReszVsRadiusLambda->Fill(rcPosZLambda,rcPosZLambda-mcPosZ);
	  fHistAsMcProdRadiusLambda->Fill(mcPosMotherR);
	  fHistAsMcProdRadiusXvsYLambda->Fill(mcPosMotherX,mcPosMotherY);
	  fHistAsMcResPtLambda->Fill(deltaPtLambda);
	  fHistAsMcResPtVsRapLambda->Fill(deltaPtLambda,lRapLambda);
	  fHistAsMcResPtVsPtLambda->Fill(deltaPtLambda,lPtLambda);
	  if (lComeFromSigma) fHistAsMcPtLambdaFromSigma->Fill(lPtLambda);
	  switch (lPdgcodeMotherOfMother) {
	  case 3222 : fHistAsMcMotherPdgCodeLambda->Fill(0.5); break; // Sigma +
	  case 3212 : fHistAsMcMotherPdgCodeLambda->Fill(1.5); break; // Sigma 0
	  case 3112 : fHistAsMcMotherPdgCodeLambda->Fill(2.5); break;// Sigma -
	  case 3224 : fHistAsMcMotherPdgCodeLambda->Fill(3.5); break;// Sigma(1385) +
	  case 3214 : fHistAsMcMotherPdgCodeLambda->Fill(4.5); break;// Sigma(1385) 0
	  case 3114 : fHistAsMcMotherPdgCodeLambda->Fill(5.5); break;// Sigma(1385) -
	  case 3322 : fHistAsMcMotherPdgCodeLambda->Fill(6.5); break; // Xi 0
	  case 3312 : fHistAsMcMotherPdgCodeLambda->Fill(7.5); break; // Xi -
	  case 3334 : fHistAsMcMotherPdgCodeLambda->Fill(8.5); break; // Omega
	  case -1   : fHistAsMcMotherPdgCodeLambda->Fill(9.5); break;
	  default   : fHistAsMcMotherPdgCodeLambda->Fill(10.5);break; 
	  }
   
	  //printf("found Lambda RC dcaPos=%e dcaNeg=%e dcaDau=%e cosP=%e pT=%e mass=%e\n",lDcaPosToPrimVertex ,lDcaNegToPrimVertex ,lDcaV0Daughters,lV0cosPointAngle,lPtLambda,lInvMassLambda);
	  //printf("found Lambda RC Pindex=%d  Nindex=%d  Plabel=%d  Nlabel=%d\n\n",lIndexTrackPos,lIndexTrackNeg,lLabelTrackPos,lLabelTrackNeg);
	    
	}
	  
	else if (lCheckSecondaryLambda) {
	  fHistAsMcSecondaryPtVsRapLambda->Fill(lPtLambda,lRapLambda);
	  fHistAsMcSecondaryProdRadiusLambda->Fill(mcPosMotherR); 
	  fHistAsMcSecondaryProdRadiusXvsYLambda->Fill(mcPosMotherX,mcPosMotherY);
	  if (lComeFromSigma) fHistAsMcSecondaryPtLambdaFromSigma->Fill(lPtLambda);
	  printf(" lPdgcodeMotherOfMother= %d",lPdgcodeMotherOfMother);
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
	  
      }
    } // end rapidity condition
    }// end PID condition
    // Associated AntiLambda histograms in |rap| < lCutRap
    if ((AntiLambdaPID==1 && lPtAntiLambda <=1) || (lPtAntiLambda>1) ||  !(fUsePID.Contains("withPID"))){          

    if (TMath::Abs(lRapAntiLambda) < lCutRap) {

         fHistNsigmaPosProtonAntiLambda->Fill(nSigmaPosProton);
         fHistNsigmaNegPionAntiLambda->Fill(nSigmaNegPion);
      switch (lOnFlyStatus){
      case 0 : 

	//          fHistMultVsPtVsMassAntiLambda->Fill(multiplicity ,lInvMassAntiLambda,lPtAntiLambda);
	if(lCheckPIdAntiLambda) fHistPidMcMassAntiLambda->Fill(lInvMassAntiLambda);
	  
		if(lCheckMcAntiLambda) {
	  fHistAsMcMassAntiLambda->Fill(lInvMassAntiLambda);
	  fHistAsMcPtAntiLambda->Fill(lPtAntiLambda);


	  fHistAsMcPtVsMassAntiLambda->Fill(lInvMassAntiLambda,lPtAntiLambda);
	  if (lPtAntiLambda <= 1) fHistAsMcPtZoomAntiLambda->Fill(lPtAntiLambda);
	  fHistAsMcMassVsRadiusAntiLambda->Fill(rcPosRAntiLambda,lInvMassAntiLambda);
	  fHistAsMcResxAntiLambda->Fill(rcPosXAntiLambda-mcPosX);
	  fHistAsMcResyAntiLambda->Fill(rcPosYAntiLambda-mcPosY);
	  fHistAsMcReszAntiLambda->Fill(rcPosZAntiLambda-mcPosZ);
	  fHistAsMcResrVsRadiusAntiLambda->Fill(rcPosRAntiLambda,rcPosRAntiLambda-mcPosR);
	  fHistAsMcReszVsRadiusAntiLambda->Fill(rcPosZAntiLambda,rcPosZAntiLambda-mcPosZ);
	  fHistAsMcProdRadiusAntiLambda->Fill(mcPosMotherR);
	  fHistAsMcProdRadiusXvsYAntiLambda->Fill(mcPosMotherX,mcPosMotherY);
	  fHistAsMcResPtAntiLambda->Fill(deltaPtAntiLambda);
	  fHistAsMcResPtVsRapAntiLambda->Fill(deltaPtAntiLambda,lRapAntiLambda);
	  fHistAsMcResPtVsPtAntiLambda->Fill(deltaPtAntiLambda,lPtAntiLambda);
	  if (lComeFromSigma) fHistAsMcPtAntiLambdaFromSigma->Fill(lPtAntiLambda);
	  switch (lPdgcodeMotherOfMother) {
	  case 3222 : fHistAsMcMotherPdgCodeAntiLambda->Fill(0.5); break; // Sigma +
	  case 3212 : fHistAsMcMotherPdgCodeAntiLambda->Fill(1.5); break; // Sigma 0
	  case 3112 : fHistAsMcMotherPdgCodeAntiLambda->Fill(2.5); break;// Sigma -
	  case 3224 : fHistAsMcMotherPdgCodeAntiLambda->Fill(3.5); break;// Sigma(1385) +
	  case 3214 : fHistAsMcMotherPdgCodeAntiLambda->Fill(4.5); break;// Sigma(1385) 0
	  case 3114 : fHistAsMcMotherPdgCodeAntiLambda->Fill(5.5); break;// Sigma(1385) -
	  case 3322 : fHistAsMcMotherPdgCodeAntiLambda->Fill(6.5); break; // Xi 0
	  case 3312 : fHistAsMcMotherPdgCodeAntiLambda->Fill(7.5); break; // Xi -
	  case 3334 : fHistAsMcMotherPdgCodeAntiLambda->Fill(8.5); break; // Omega

	  case -1   : fHistAsMcMotherPdgCodeAntiLambda->Fill(9.5); break;
	  default   : fHistAsMcMotherPdgCodeAntiLambda->Fill(10.5);break; 
	  }
    
	  //printf("found AntiLambda RC dcaPos=%e dcaNeg=%e dcaDau=%e cosP=%e pT=%e mass=%e\n",lDcaPosToPrimVertex ,lDcaNegToPrimVertex ,lDcaV0Daughters,lV0cosPointAngle,lPtAntiLambda,lInvMassAntiLambda);
	  //printf("found AntiLambda RC Pindex=%d  Nindex=%d  Plabel=%d  Nlabel=%d\n\n",lIndexTrackPos,lIndexTrackNeg,lLabelTrackPos,lLabelTrackNeg);
	    
	}
	  
		else if (lCheckSecondaryAntiLambda) {
	  fHistAsMcSecondaryPtVsRapAntiLambda->Fill(lPtAntiLambda,lRapAntiLambda);
	  fHistAsMcSecondaryProdRadiusAntiLambda->Fill(mcPosMotherR); 
	  fHistAsMcSecondaryProdRadiusXvsYAntiLambda->Fill(mcPosMotherX,mcPosMotherY);
	  if (lComeFromSigma) fHistAsMcSecondaryPtAntiLambdaFromSigma->Fill(lPtAntiLambda);
	  printf(" lPdgcodeMotherOfMother= %d",lPdgcodeMotherOfMother);
	  switch (lPdgcodeMotherOfMother) {
	  case 3222 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(0.5); break;// Sigma +
	  case 3212 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(1.5); break;// Sigma 0
	  case 3112 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(2.5); break;// Sigma -
	  case 3224 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(3.5); break;// Sigma(1385) +
	  case 3214 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(4.5); break;// Sigma(1385) 0
	  case 3114 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(5.5); break;// Sigma(1385) -
	  case 3322 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(6.5); break; // Xi 0
	  case 3312 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(7.5); break; // Xi -
	  case 3334 : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(8.5); break; // Omega
	  case -1   : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(9.5); break;
	  default   : fHistAsMcSecondaryMotherPdgCodeAntiLambda->Fill(10.5);break;
	  }
	}
    	break;
	  
      }
    } // end rapidity condition
    }// end PID condition       

    if (negPiKF) delete negPiKF; negPiKF= NULL;
    if (posPiKF) delete posPiKF; posPiKF= NULL;
    if (posPKF)  delete posPKF;  posPKF = NULL;
    if (negAPKF) delete negAPKF; negAPKF= NULL;
    
  } // end V0 loop

  //  if (primaryVtxKF) delete primaryVtxKF;primaryVtxKF=NULL ;


  fHistV0Multiplicity->Fill(nv0s);
  //  fHistV0MultiplicityMI->Fill(nv0sMI);

  if (fAnalysisType == "ESD") { if(myPrimaryVertex) delete myPrimaryVertex; }

  
  // Post output data
}      

//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::Terminate(Option_t *) 
{/*
 // Draw result to the screen
 // Called once at the end of the query

 TList *cRetrievedList = 0x0;
 cRetrievedList = (TList*)GetOutputData(1);
  
 if(!cRetrievedList){
 AliWarning("ERROR - AliAnalysisTaskPerformanceStrange: output data container list not available\n"); return;
 }
  
  
 fHistV0Multiplicity = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistV0Multiplicity"));
 if (!fHistV0Multiplicity) {
 Printf("ERROR: fHistV0Multiplicity not available");
 return;
 }

 fHistV0MultiplicityMI = dynamic_cast<TH1F*> ( cRetrievedList->FindObject("fHistV0MultiplicityMI"));
 if (!fHistV0MultiplicityMI) {
 Printf("ERROR: fHistV0MultiplicityMI not available");
 return;
 }

 TCanvas *canPerformanceStrange = new TCanvas("AliAnalysisTaskCheckV0","Multiplicity",10,10,510,510);
 canPerformanceStrange->Divide(2,1);
 if (fHistV0Multiplicity->GetMaximum() > 0.) canPerformanceStrange->cd(1)->SetLogy();
 fHistV0Multiplicity->SetMarkerStyle(25);
 fHistV0Multiplicity->DrawCopy("E");
 if (fHistV0MultiplicityMI->GetMaximum() > 0.) canPerformanceStrange->cd(2)->SetLogy();
 fHistV0MultiplicityMI->SetMarkerStyle(24);
 fHistV0MultiplicityMI->DrawCopy("E");
  

 */ 
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskPerformanceStrange::MyRapidity(Double_t rE, Double_t rPz) const
{
  // Local calculation for rapidity
  return 0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
} 
//----------------------------------------------------------------------------

