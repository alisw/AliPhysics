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


ClassImp(AliAnalysisTaskPerformanceStrange)


//________________________________________________________________________
AliAnalysisTaskPerformanceStrange::AliAnalysisTaskPerformanceStrange()
  : AliAnalysisTaskSE("TaskStrange"), fAnalysisMC(0), fAnalysisType("infoType"),  fCollidingSystems(0), fUsePID("infoPID"), fUseCut("infoCut"), fUseOnTheFly(0), fESD(0), fListHist(0), fCentrSelector(0),

    // MC histograms  ---------------------------------------
    fHistMCPrimaryVertexX(0),
    fHistMCPrimaryVertexY(0),
    fHistMCPrimaryVertexZ(0),

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

    fHistMCPtK0s(0),
    fHistMCPtLambda(0),

    fHistMCPtLambdaFromSigma(0),
    fHistMCPtAntiLambdaFromSigma(0),

    fHistNTimesRecK0s(0),
    fHistNTimesRecLambda(0),
    fHistNTimesRecAntiLambda(0),
    fHistNTimesRecK0sVsPt(0),
    fHistNTimesRecLambdaVsPt(0),
    fHistNTimesRecAntiLambdaVsPt(0),
    // ------------------------------------------------------

    // Reconstructed particle histograms  -------------------
    fHistNumberEvents(0),
    fHistTrackPerEvent(0),
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

    fHistV0Multiplicity(0),

    fHistChi2KFBeforeCutK0s(0), 
    fHistChi2KFBeforeCutLambda(0), 
    fHistChi2KFBeforeCutAntiLambda(0),
    fHistChi2KFAfterCutK0s(0), 
    fHistChi2KFAfterCutLambda(0), 
    fHistChi2KFAfterCutAntiLambda(0),

    fHistMassK0(0),
    fHistMassLambda(0),
    fHistMassAntiLambda(0),
    fHistMassVsRadiusK0(0),
    fHistMassVsRadiusLambda(0),
    fHistMassVsRadiusAntiLambda(0),

    fHistPtVsMassK0(0),
    fHistPtVsMassLambda(0),
    fHistArmenterosPodolanski(0),
    // ------------------------------------------------------

    // kontrola pre |pz/pt| histograms  ---------------------
    fHistPzPtBeforeK0s(0),
    fHistPzPtAfterK0s(0),
    fHistPzPtBeforeLambda(0),
    fHistPzPtAfterLambda(0),
    // ------------------------------------------------------

    // PID histograms  --------------------------------------
    fHistNsigmaPosPionAntiLambda(0),
    fHistNsigmaNegProtonAntiLambda(0),
    fHistNsigmaPosProtonLambda(0),
    fHistNsigmaNegPionLambda(0),
    fHistNsigmaPosPionK0(0),
    fHistNsigmaNegPionK0(0),
    // ------------------------------------------------------

    // Associated particles ---------------------------------
    fHistAsMcRapK0(0),
    fHistAsMcRapLambda(0),
    fHistAsMcRapAntiLambda(0),

    fHistAsMcPtK0(0),
    fHistAsMcPtLambda(0),

    fHistAsMcPtZoomK0(0),
    fHistAsMcPtZoomLambda(0),

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
    // ------------------------------------------------------

    // Associated secondary particle histograms -------------
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

  // New V0 cuts
/*  fCuts[0]=33;    // max allowed chi2
  fCuts[1]=0.05;  // min allowed impact parameter for the 1st daughter
  fCuts[2]=0.05;  // min allowed impact parameter for the 2nd daughter
  fCuts[3]=0.5;   // max allowed DCA between the daughter tracks
  fCuts[4]=0.00;  // max allowed cosine of V0's pointing angle
  fCuts[5]=0.2;   // min radius of the fiducial volume
  fCuts[6]=100;   // max radius of the fiducial volume
*/
}

//________________________________________________________________________
AliAnalysisTaskPerformanceStrange::AliAnalysisTaskPerformanceStrange(const char *name)
  : AliAnalysisTaskSE(name), fAnalysisMC(0), fAnalysisType("infoType"), fCollidingSystems(0), fUsePID("infoPID"), fUseCut("infoCut"), fUseOnTheFly(0), fESD(0), fListHist(),fCentrSelector(0),

    // MC histograms  ---------------------------------------
    fHistMCPrimaryVertexX(0),
    fHistMCPrimaryVertexY(0),
    fHistMCPrimaryVertexZ(0),

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

    fHistMCPtK0s(0),
    fHistMCPtLambda(0),

    fHistMCPtLambdaFromSigma(0),
    fHistMCPtAntiLambdaFromSigma(0),

    fHistNTimesRecK0s(0),
    fHistNTimesRecLambda(0),
    fHistNTimesRecAntiLambda(0),
    fHistNTimesRecK0sVsPt(0),
    fHistNTimesRecLambdaVsPt(0),
    fHistNTimesRecAntiLambdaVsPt(0),
    // ------------------------------------------------------

    // Reconstructed particle histograms  -------------------
    fHistNumberEvents(0),
    fHistTrackPerEvent(0),
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

    fHistV0Multiplicity(0),

    fHistChi2KFBeforeCutK0s(0), 
    fHistChi2KFBeforeCutLambda(0), 
    fHistChi2KFBeforeCutAntiLambda(0),
    fHistChi2KFAfterCutK0s(0), 
    fHistChi2KFAfterCutLambda(0), 
    fHistChi2KFAfterCutAntiLambda(0),

    fHistMassK0(0),
    fHistMassLambda(0),
    fHistMassAntiLambda(0),
    fHistMassVsRadiusK0(0),
    fHistMassVsRadiusLambda(0),
    fHistMassVsRadiusAntiLambda(0),

    fHistPtVsMassK0(0),
    fHistPtVsMassLambda(0),
    fHistArmenterosPodolanski(0),
    // ------------------------------------------------------

    // kontrola pre |pz/pt| histograms  ---------------------
    fHistPzPtBeforeK0s(0),
    fHistPzPtAfterK0s(0),
    fHistPzPtBeforeLambda(0),
    fHistPzPtAfterLambda(0),
    // ------------------------------------------------------

    // PID histograms  --------------------------------------
    fHistNsigmaPosPionAntiLambda(0),
    fHistNsigmaNegProtonAntiLambda(0),
    fHistNsigmaPosProtonLambda(0),
    fHistNsigmaNegPionLambda(0),
    fHistNsigmaPosPionK0(0),
    fHistNsigmaNegPionK0(0),
    // ------------------------------------------------------

    // Associated particles ---------------------------------
    fHistAsMcRapK0(0),
    fHistAsMcRapLambda(0),
    fHistAsMcRapAntiLambda(0),

    fHistAsMcPtK0(0),
    fHistAsMcPtLambda(0),

    fHistAsMcPtZoomK0(0),
    fHistAsMcPtZoomLambda(0),

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
    // ------------------------------------------------------

    // Associated secondary particle histograms -------------
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

  //New V0 cuts
/*  fCuts[0]=33;    // max allowed chi2
  fCuts[1]=0.05;  // min allowed impact parameter for the 1st daughter
  fCuts[2]=0.05;  // min allowed impact parameter for the 2nd daughter
  fCuts[3]=0.5;   // max allowed DCA between the daughter tracks
  fCuts[4]=0.00;  // max allowed cosine of V0's pointing angle
  fCuts[5]=0.2;   // min radius of the fiducial volume
  fCuts[6]=100;   // max radius of the fiducial volume
*/
  // Define output slots only here
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliAnalysisCentralitySelector::Class());

}

//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::UserCreateOutputObjects() 
{
  //******************
  // Create histograms
  //*******************
  fListHist = new TList();
  fListHist->SetOwner();

  // Bo: tbd: condition before allocation (i.e. if (!fHistMCMultiplicityPrimary){...} for each histo...

  //***************
  // MC histograms
  //***************
 
  // Primary Vertex X,Y,Z:
  fHistMCPrimaryVertexX          = new TH1F("h1MCPrimaryVertexX", "MC Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistMCPrimaryVertexX);

  fHistMCPrimaryVertexY          = new TH1F("h1MCPrimaryVertexY", "MC Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistMCPrimaryVertexY);

  fHistMCPrimaryVertexZ          = new TH1F("h1MCPrimaryVertexZ", "MC Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fListHist->Add(fHistMCPrimaryVertexZ);
  
  // Multiplicity:
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

  fHistMCRapXi                  = new TH1F("h1MCRapXi", "#Xi;y",160,-4,4);
  fListHist->Add(fHistMCRapXi);

  fHistMCRapInPtRangeXi         = new TH1F("h1MCRapInPtRangeXi", "#Xi;y",160,-4,4);
  fListHist->Add(fHistMCRapInPtRangeXi);

  fHistMCRapPhi                  = new TH1F("h1MCRapPhi", "#phi;y",160,-4,4);
  fListHist->Add(fHistMCRapPhi);

  fHistMCRapInPtRangePhi         = new TH1F("h1MCRapInPtRangePhi", "#phi;y",160,-4,4);
  fListHist->Add(fHistMCRapInPtRangePhi);


  // Pt distribution:
  fHistMCPtK0s               = new TH1F("h1MCPtK0s", "K^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtK0s);

  fHistMCPtLambda            = new TH1F("h1MCPtLambda", "#Lambda^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtLambda);

  // Pt distribution of Lambda coming from Sigma decay:
  fHistMCPtLambdaFromSigma      = new TH1F("h1MCPtLambdaFromSigma", "#Lambda^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtLambdaFromSigma);

  fHistMCPtAntiLambdaFromSigma  = new TH1F("h1MCPtAntiLambdaFromSigma", "#Lambda^{0};p_{t} (GeV/c)",240,0,12);
  fListHist->Add(fHistMCPtAntiLambdaFromSigma);
 
  // Multiple reconstruction studies:
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

  // Number of events:
  fHistNumberEvents           = new TH1F("h1NumberEvents", "Number of events; index;Number of Events",10,0,10);
  fListHist->Add(fHistNumberEvents);

  // Multiplicity:
  fHistTrackPerEvent           = new TH1F("h1TrackPerEvent", "Tracks per event;Number of Tracks;Number of Events",1000,0,1000);
  fListHist->Add(fHistTrackPerEvent);

  fHistTrackletPerEvent       = new TH1F("h1TrackletPerEvent", "Number of tracklets;Number of tracklets per events;Number of events",1000,0,1000);
  fListHist->Add(fHistTrackletPerEvent);

  fHistMCDaughterTrack         = new TH1F("h1MCDaughterTrack","Distribution of mc id for daughters;id tags;Counts",15,0,15);
  fListHist->Add(fHistMCDaughterTrack);

  // Primary Vertex:
  fHistSPDPrimaryVertexZ          = new TH1F("h1SPDPrimaryVertexZ", "SPD Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fListHist->Add(fHistSPDPrimaryVertexZ);

  fHistPrimaryVertexX          = new TH1F("h1PrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexX);

  fHistPrimaryVertexY          = new TH1F("h1PrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexY);

  fHistPrimaryVertexZ          = new TH1F("h1PrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
  fListHist->Add(fHistPrimaryVertexZ);


  // Primary vertex resolution:
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

  // V0 Multiplicity:
  if (!fHistV0Multiplicity) {
    if (fCollidingSystems)
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of Offline V0s;Events", 200, 0, 40000);
    else
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of Offline V0s;Events", 10, 0, 10); 
    fListHist->Add(fHistV0Multiplicity);
  }

  // Kalman Filter Chi2:
  fHistChi2KFBeforeCutK0s               = new TH2F("h1Chi2KFBeforeCutK0s", "K^{0}  candidates;#Chi^{2});Counts", 250, 0, 50, 2,-0.5,1.5);
  fListHist->Add(fHistChi2KFBeforeCutK0s);
  fHistChi2KFBeforeCutLambda            = new TH2F("h1Chi2KFBeforeCutLambda", "#Lambda^{0}  candidates;#Chi^{2};Counts", 250, 0, 50, 2,-0.5,1.5);
  fListHist->Add(fHistChi2KFBeforeCutLambda);
  fHistChi2KFBeforeCutAntiLambda        = new TH2F("h1Chi2KFBeforeCutAntiLambda", "#bar{#Lambda}^{0}  candidates;#Chi^{2};Counts", 250, 0, 50, 2,-0.5,1.5);
  fListHist->Add(fHistChi2KFBeforeCutAntiLambda);

  fHistChi2KFAfterCutK0s               = new TH2F("h1Chi2KFAfterCutK0s", "K^{0}  candidates;#Chi^{2});Counts", 250, 0, 50, 2,-0.5,1.5);
  fListHist->Add(fHistChi2KFAfterCutK0s);
  fHistChi2KFAfterCutLambda            = new TH2F("h1Chi2KFAfterCutLambda", "#Lambda^{0}  candidates;#Chi^{2};Counts", 250, 0, 50, 2,-0.5,1.5);
  fListHist->Add(fHistChi2KFAfterCutLambda);
  fHistChi2KFAfterCutAntiLambda        = new TH2F("h1Chi2KFAfterCutAntiLambda", "#bar{#Lambda}^{0} candidates;#Chi^{2};Counts", 250, 0, 50, 2,-0.5,1.5);
  fListHist->Add(fHistChi2KFAfterCutAntiLambda);

  // Invariant mass:
  fHistMassK0                   = new TH1F("h1MassK0", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistMassK0);

  fHistMassLambda               = new TH1F("h1MassLambda", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassLambda);

  fHistMassAntiLambda           = new TH1F("h1MassAntiLambda", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambda);

  // Invariant mass vs. radius:
  const Double_t radius[10] = {0.0,2.5,2.9,3.9,7.6,15.0,23.9,37.8,42.8,100.0};
  Int_t lNbinRadius        = 9;
  Int_t lNbinInvMassLambda = 300;

  fHistMassVsRadiusK0           = new TH2F("h2MassVsRadiusK0", "K^{0} candidates;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 200, 0.4, 0.6);
  fListHist->Add(fHistMassVsRadiusK0);
  
  fHistMassVsRadiusLambda       = new TH2F("h2MassVsRadiusLambda", "#Lambda candidates;radius (cm);M(p#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusLambda);

  fHistMassVsRadiusAntiLambda   = new TH2F("h2MassVsRadiusAntiLambda", "#bar{#Lambda} candidates;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",lNbinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusAntiLambda);

  // Pt vs. mass:
  fHistPtVsMassK0               = new TH2F("h2PtVsMassK0","K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,240,0,12);
  fListHist->Add(fHistPtVsMassK0);

  fHistPtVsMassLambda           = new TH2F("h2PtVsMassLambda","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistPtVsMassLambda);

  fHistArmenterosPodolanski     = new TH2F("h2ArmenterosPodolanski","Armenteros-Podolanski phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);
  fListHist->Add(fHistArmenterosPodolanski);

  // Kontrola pre |pz/pt| histograms:
  fHistPzPtBeforeK0s	= new TH1F("h1PzPtBeforeK0s","K0s; Abs(pz/pt); count",1000,0,10);
  fListHist->Add(fHistPzPtBeforeK0s);

  fHistPzPtAfterK0s	= new TH1F("h1PzPtAfterK0s","K0s; Abs(pz/pt); count",1000,0,10);
  fListHist->Add(fHistPzPtAfterK0s);

  fHistPzPtBeforeLambda	= new TH1F("h1PzPtBeforeLambda","#Lambda^{0}; Abs(pz/pt); count",1000,0,10);
  fListHist->Add(fHistPzPtBeforeLambda);

  fHistPzPtAfterLambda	= new TH1F("h1PzPtAfterLambda","#Lambda^{0}; Abs(pz/pt); count",1000,0,10);
  fListHist->Add(fHistPzPtAfterLambda);

  // PID histograms:
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
  // Associated particle histograms
  //********************************

  // Rapidity distribution:
  fHistAsMcRapK0                = new TH1F("h1AsMcRapK0", "K^{0} associated;eta;Counts", 60, -1.5, 1.5);
  fListHist->Add(fHistAsMcRapK0);

  fHistAsMcRapLambda            = new TH1F("h1AsMcRapLambda", "#Lambda^{0} associated;eta;Counts", 60, -1.5, 1.5);
  fListHist->Add(fHistAsMcRapLambda);

  fHistAsMcRapAntiLambda        = new TH1F("h1AsMcRapAntiLambda", "#bar{#Lambda}^{0} associated;eta;Counts", 60, -1.5, 1.5);
  fListHist->Add(fHistAsMcRapAntiLambda);

  //Pt distribution:
  fHistAsMcPtK0                = new TH1F("h1AsMcPtK0", "K^{0} associated;p_{t} (GeV/c);Counts", 240,0,12);
  fListHist->Add(fHistAsMcPtK0);

  fHistAsMcPtLambda            = new TH1F("h1AsMcPtLambda", "#Lambda^{0} associated;p_{t} (GeV/c);Counts", 240,0,12);
  fListHist->Add(fHistAsMcPtLambda);

  fHistAsMcPtZoomK0            = new TH1F("h1AsMcPtZoomK0", "K^{0} candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomK0);

  fHistAsMcPtZoomLambda        = new TH1F("h1AsMcPtZoomLambda", "#Lambda^{0} candidates in -1 <y<1;p_{t} (GeV/c);Counts",20,0,1);
  fListHist->Add(fHistAsMcPtZoomLambda);

  // Radius distribution:
  fHistAsMcProdRadiusK0               = new TH1F("h1AsMcProdRadiusK0", "K^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusK0);

  fHistAsMcProdRadiusLambda           = new TH1F("h1AsMcProdRadiusLambda", "#Lambda^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusLambda);

  fHistAsMcProdRadiusAntiLambda       = new TH1F("h1AsMcProdRadiusAntiLambda", "#bar{#Lambda}^{0} associated;r (cm);Counts", 500, 0, 100);
  fListHist->Add(fHistAsMcProdRadiusAntiLambda);

  // Radius distribution vs. rapidity:
  fHistAsMcProdRadiusXvsYK0s          = new TH2F("h2AsMcProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYK0s);

  fHistAsMcProdRadiusXvsYLambda       = new TH2F("h2AsMcProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYLambda);

  fHistAsMcProdRadiusXvsYAntiLambda   = new TH2F("h2AsMcProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
  fListHist->Add(fHistAsMcProdRadiusXvsYAntiLambda);

  // Invariant mass distribution with PID checked:
  fHistPidMcMassK0             = new TH1F("h1PidMcMassK0", "K^{0} MC PId checked;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistPidMcMassK0);

  fHistPidMcMassLambda         = new TH1F("h1PidMcMassLambda", "#Lambda^{0} MC PId checked;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassLambda);
  
  fHistPidMcMassAntiLambda     = new TH1F("h1PidMcMassAntiLambda", "#bar{#Lambda}^{0} MC PId checked;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistPidMcMassAntiLambda);

  // Invariant mass distribution:
  fHistAsMcMassK0              = new TH1F("h1AsMcMassK0", "K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistAsMcMassK0);
  
  fHistAsMcMassLambda          = new TH1F("h1AsMcMassLambda", "#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassLambda);

  fHistAsMcMassAntiLambda      = new TH1F("h1AsMcMassAntiLambda", "#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistAsMcMassAntiLambda);

  // Pt vs. invariant mass:
  fHistAsMcPtVsMassK0               = new TH2F("h2AsMcPtVsMassK0","K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,240,0,12);
  fListHist->Add(fHistAsMcPtVsMassK0);

  fHistAsMcPtVsMassLambda           = new TH2F("h2AsMcPtVsMassLambda","#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistAsMcPtVsMassLambda);

  fHistAsMcPtVsMassAntiLambda       = new TH2F("h2AsMcPtVsMassAntiLambda","#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,240,0,12);
  fListHist->Add(fHistAsMcPtVsMassAntiLambda);


  // Invariant mass vs. radius:
  fHistAsMcMassVsRadiusK0             = new TH2F("h2AsMcMassVsRadiusK0", "K^{0} associated;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, 500, 0.47, 0.52);
  fListHist->Add(fHistAsMcMassVsRadiusK0);
  
  fHistAsMcMassVsRadiusLambda         = new TH2F("h2AsMcMassVsRadiusLambda", "#Lambda associated;radius (cm);M(p#pi^{-}) (GeV/c^{2})",lNbinRadius,radius, lNbinInvMassLambda, 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusLambda);

  fHistAsMcMassVsRadiusAntiLambda     = new TH2F("h2AsMcMassVsRadiusAntiLambda", "#bar{#Lambda} associated;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",lNbinRadius,radius,lNbinInvMassLambda , 1.10, 1.13);
  fListHist->Add(fHistAsMcMassVsRadiusAntiLambda);
  
  // Position resolution for K0s:
  fHistAsMcResxK0                     = new TH1F("h1AsMcResxK0", "K^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxK0);
  fHistAsMcResyK0                     = new TH1F("h1AsMcResyK0", "K^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyK0);
  fHistAsMcReszK0                     = new TH1F("h1AsMcReszK0", "K^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszK0);

  // Position resolution vs. radius for K0s:
  fHistAsMcResrVsRadiusK0             = new TH2F("h2AsMcResrVsRadiusK0", "K^{0} associated;Radius (cm);#Delta r (cm)",200,0.0,50., 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusK0);
  fHistAsMcReszVsRadiusK0             = new TH2F("h2AsMcReszVsRadiusK0", "K^{0} associated;Radius (cm);#Delta z (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusK0);

  // Position resolution for Lambda:
  fHistAsMcResxLambda                 = new TH1F("h1AsMcResxLambda", "#Lambda^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxLambda);
  fHistAsMcResyLambda                 = new TH1F("h1AsMcResyLambda", "#Lambda^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyLambda);
  fHistAsMcReszLambda                 = new TH1F("h1AsMcReszLambda", "#Lambda^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszLambda);

  // Position resolution vs. radius for Lambda:
  fHistAsMcResrVsRadiusLambda         = new TH2F("h2AsMcResrVsRadiusLambda", "#Lambda^{0} associated;Radius (cm);#Delta r (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusLambda);
  fHistAsMcReszVsRadiusLambda         = new TH2F("h2AsMcReszVsRadiusLambda", "#Lambda^{0} associated;Radius (cm);#Delta z (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusLambda);

  // Position resolution for anti-Lambda:
  fHistAsMcResxAntiLambda             = new TH1F("h1AsMcResxAntiLambda", "#bar{#Lambda}^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResxAntiLambda);
  fHistAsMcResyAntiLambda             = new TH1F("h1AsMcResyAntiLambda", "#bar{#Lambda}^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResyAntiLambda);
  fHistAsMcReszAntiLambda             = new TH1F("h1AsMcReszAntiLambda", "#bar{#Lambda}^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszAntiLambda);

  // Position resolution vs. radius for anti-Lambda:
  fHistAsMcResrVsRadiusAntiLambda     = new TH2F("h2AsMcResrVsRadiusAntiLambda", "#bar{#Lambda}^{0} associated;Radius (cm);#Delta r (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcResrVsRadiusAntiLambda);
  fHistAsMcReszVsRadiusAntiLambda     = new TH2F("h2AsMcReszVsRadiusAntiLambda", "#bar{#Lambda}^{0} associated;Radius (cm);#Delta z (cm)",200,0.0,50.0, 50, -0.25, 0.25);
  fListHist->Add(fHistAsMcReszVsRadiusAntiLambda);

  // Pt Resolution:
  fHistAsMcResPtK0                   = new TH1F("h1AsMcResPtK0","Pt Resolution K^{0};#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtK0);
  fHistAsMcResPtLambda               = new TH1F("h1AsMcResPtLambda","Pt Resolution #Lambda^{0};#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtLambda);
  fHistAsMcResPtAntiLambda           = new TH1F("h1AsMcResPtAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Counts",200,-1,1);
  fListHist->Add(fHistAsMcResPtAntiLambda);

  // Pt Resolution vs. rapidity:
  fHistAsMcResPtVsRapK0              = new TH2F("h2AsMcResPtVsRapK0","Pt Resolution K^{0};#Delta Pt;Rap",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapK0);
  fHistAsMcResPtVsRapLambda          = new TH2F("h2AsMcResPtVsRapLambda","Pt Resolution #Lambda^{0};#Delta Pt;Rap",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapLambda);
  fHistAsMcResPtVsRapAntiLambda      = new TH2F("h2AsMcResPtVsRapAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Rap",200,-1,1,20,-1,1);
  fListHist->Add(fHistAsMcResPtVsRapAntiLambda);

  // Pt Resolution vs. Pt:
  fHistAsMcResPtVsPtK0               = new TH2F("h2AsMcResPtVsPtK0","Pt Resolution K^{0};#Delta Pt;Pt",600,-0.15,0.15,240,0,12);
  fListHist->Add(fHistAsMcResPtVsPtK0);
  fHistAsMcResPtVsPtLambda           = new TH2F("h2AsMcResPtVsPtLambda","Pt Resolution #Lambda^{0};#Delta Pt;Pt",600,-0.15,0.15,240,0,12);
  fListHist->Add(fHistAsMcResPtVsPtLambda);
  fHistAsMcResPtVsPtAntiLambda       = new TH2F("h2AsMcResPtVsPtAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Pt",300,-0.15,0.15,240,0,12);
  fListHist->Add(fHistAsMcResPtVsPtAntiLambda);

  // Pdg code of mother particle:
  fHistAsMcMotherPdgCodeK0s           = new TH1F("h1AsMcMotherPdgCodeK0s","Mother of Associated K^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeK0s);
  fHistAsMcMotherPdgCodeLambda        = new TH1F("h1AsMcMotherPdgCodeLambda","Mother of Associated #Lambda^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeLambda);
  fHistAsMcMotherPdgCodeAntiLambda    = new TH1F("h1AsMcMotherPdgCodeAntiLambda","Mother of Associated #bar{#Lambda}^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcMotherPdgCodeAntiLambda);

  // Pt distribution of Lambda <- Sigma decay
  fHistAsMcPtLambdaFromSigma          = new TH1F("h1AsMcPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",240,0,12);
  fListHist->Add(fHistAsMcPtLambdaFromSigma);
  fHistAsMcPtAntiLambdaFromSigma      = new TH1F("h1AsMcPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",240,0,12);
  fListHist->Add(fHistAsMcPtAntiLambdaFromSigma);

  //*******************************************
  // Associated secondary particles histograms
  //*******************************************

  // Pt vs. rapidity distribution:
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

  // Production radius vs. rapidity:
  fHistAsMcSecondaryProdRadiusXvsYK0s          = new TH2F("h2AsMcSecondaryProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYK0s);
  fHistAsMcSecondaryProdRadiusXvsYLambda       = new TH2F("h2AsMcSecondaryProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYLambda);
  fHistAsMcSecondaryProdRadiusXvsYAntiLambda   = new TH2F("h2AsMcSecondaryProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
  fListHist->Add(fHistAsMcSecondaryProdRadiusXvsYAntiLambda);

  // Pdg code of mother particle for secondary V0s:
  fHistAsMcSecondaryMotherPdgCodeK0s           = new TH1F("h1AsMcSecondaryMotherPdgCodeK0s","Mother of Associated Secondary K^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeK0s);
  fHistAsMcSecondaryMotherPdgCodeLambda        = new TH1F("h1AsMcSecondaryMotherPdgCodeLambda","Mother of Associated Secondary #Lambda^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeLambda);
  fHistAsMcSecondaryMotherPdgCodeAntiLambda    = new TH1F("h1AsMcSecondaryMotherPdgCodeAntiLambda","Mother of Associated Secondary #bar{#Lambda}^{0};mother;counts",11,0,11);
  fListHist->Add(fHistAsMcSecondaryMotherPdgCodeAntiLambda);

  // Pt distribution of secondary Lambda <- Sigma decay:
  fHistAsMcSecondaryPtLambdaFromSigma          = new TH1F("h1AsMcSecondaryPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",240,0,12);
  fListHist->Add(fHistAsMcSecondaryPtLambdaFromSigma);
  fHistAsMcSecondaryPtAntiLambdaFromSigma      = new TH1F("h1AsMcSecondaryPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p_{t} (GeV/c);Count",240,0,12);
  fListHist->Add(fHistAsMcSecondaryPtAntiLambdaFromSigma);

  PostData(1, fListHist);
  PostData(2, fCentrSelector);
}

//________________________________________________________________________
void AliAnalysisTaskPerformanceStrange::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  AliStack* stack = NULL;
  TClonesArray *mcArray = NULL;
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

  fHistNumberEvents->Fill(0.5);

  //******************
  // Trigger Selection ! Warning Works only for ESD, add protection in case of AOD loop
  //******************

  Bool_t isSelected = 
    (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() 
     & AliVEvent::kMB);
  if (!isSelected) return;
  
  // Centrality selection
  static AliESDtrackCuts * trackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(); // FIXME: make it a data member
  Bool_t isCentralitySelected = fCentrSelector->IsCentralityBinSelected(fESD,trackCuts); 
  if(!isCentralitySelected) return;
  // FIXME: add to hist number events another entry for centrality.

 // Done by the AliPhysicsSelection Task ! Only the selected events are passed to this task

  fHistNumberEvents->Fill(1.5); // FIXME: use enum here


  //*************************
  //End track multiplicity
  //*************************


  // Remove Events with no tracks
  //if (!(fESD->GetNumberOfTracks()))  return;

  fHistNumberEvents->Fill(2.5);
  fHistTrackPerEvent->Fill(fESD->GetNumberOfTracks());

  //*************************************
  // Cut used:
  //*************************************
  // FIXME: Create a cut object, to be configured in the steering macro and to be streamed in the output to reference those cuts
  // Cut Rapidity:
  Double_t lCutRap  = 0.75;

  // Cut AliKF Chi2 for Reconstructed particles
  Double_t cutChi2KF  = 1E3;

  // If PID is used:
  Double_t lLimitPPID    = 0.7;
  Float_t cutNSigmaLowP  = 1E3;
  Float_t cutNSigmaHighP = 1E3;
  if (fUsePID.Contains("withPID")) {
    cutNSigmaLowP  = 5.0;
    cutNSigmaHighP = 3.0;
  }


  // Cut Daughters pt (GeV/c):
  Double_t cutMinPtDaughter = 0.160;

  // Cut primary vertex:
  Double_t cutPrimVertex = 10.0;

  // Min number of TPC clusters:
  Int_t nbMinTPCclusters = 80;

  //*******************
  // PID parameters:
  //*******************
  // FIXME: OADB or momber TFormula?
  Double_t fAlephParameters[5] = {0,0,0,0,0,};

  fAlephParameters[0] = 0.0283086;
  fAlephParameters[1] = 2.63394e+01;
  fAlephParameters[2] = 5.04114e-11;
  fAlephParameters[3] = 2.12543e+00;
  fAlephParameters[4] = 4.88663e+00; 


  //*******************
  // Access MC:
  //*******************
  // FIXME:: move this two branches directly in the loops below
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
    
    else if(fAnalysisType == "AOD") {
      
      // load MC particles
      mcArray = (TClonesArray*)fESD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if(!mcArray) {
	Printf("strange analysis::UserExec: MC particles branch not found!\n");
	return;
      }
      
      // load MC header
      AliAODMCHeader *mcHeader = 
	(AliAODMCHeader*)fESD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
      if(!mcHeader) {
	Printf("strange analysis::UserExec: MC header branch not found!\n");
	return;
      }
    }

    // PID parameters for MC simulations:
    // FIXME: set above, with the others
    fAlephParameters[0] = 2.15898e+00/50.;
    fAlephParameters[1] = 1.75295e+01;
    fAlephParameters[2] = 3.40030e-09;
    fAlephParameters[3] = 1.96178e+00;
    fAlephParameters[4] = 3.91720e+00; 
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
  Bool_t lCurrentMotherIsPrimary;

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
  
    // FIXME: move these loops to other functions, for better readibility?
    if(fAnalysisType == "ESD") {
      
      lNbMCPrimary = stack->GetNprimary(); // FIXME: This does not correspond to our definition of primaries, but maybe it is ok for strange particles
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

	// Decay Radius
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
	  //FIXME: shouldn't this be a fatal?
	  //Printf("ERROR: particle with label %d and/or %d not found in stack (mc loop)", id0,id1);
	  mcDecayPosR = -1.0;
	}
	// FIXME using array of histos and conversion PDGCode -> enum would make this much easier to read
	// We could also have a function FillMcHistos (pos, radius, rap...) which we call from both the AOD and ESD loops
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
	
	  // FIXME: not sure if I understand this: is it correct? (definition of primaries)
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
	//check if V0 is reconstructed several times  
     
	lNtimesReconstructedK0s   = 0; lNtimesReconstructedLambda   = 0; lNtimesReconstructedAntiLambda   = 0;

	//for (Int_t jV0 = 0; jV0 < fESD->GetNumberOfV0s(); jV0++) {
	
	//lLabelTrackN  = 0; lLabelTrackP = 0;
	//lPartNMother  = 0; lPartPMother = 0;
	
	//AliESDv0    *vertexESD = ((AliESDEvent*)fESD)->GetV0(jV0);
	//if (!vertexESD) continue;
	
	//AliESDtrack *trackNESD = ((AliESDEvent*)fESD)->GetTrack(TMath::Abs(vertexESD->GetNindex()));
	//lLabelTrackN = (UInt_t)TMath::Abs(trackNESD->GetLabel());
	//if (lLabelTrackN!=id0 && lLabelTrackN!=id1) continue;
	
	//AliESDtrack *trackPESD = ((AliESDEvent*)fESD)->GetTrack(TMath::Abs(vertexESD->GetPindex()));
	//lLabelTrackP = (UInt_t)TMath::Abs(trackPESD->GetLabel());
	//if (lLabelTrackP!=id0 && lLabelTrackP!=id1) continue;
	
	//TParticle   *lPartNESD = stack->Particle(lLabelTrackN);
	//TParticle   *lPartPESD = stack->Particle(lLabelTrackP);
	//lPartNMother = lPartNESD->GetFirstMother();
	//lPartPMother = lPartPESD->GetFirstMother();
	
	//lPartPMotherPDGcode = stack->Particle(lPartPMother)->GetPdgCode();
	
	//switch (vertexESD->GetOnFlyStatus()){
	
	//case 0 : 
	//if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==310) ) lNtimesReconstructedK0s++;
	//else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==3122) ) lNtimesReconstructedLambda++;
	//else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==-3122) ) lNtimesReconstructedAntiLambda++;
	//break;
	
	//case 1 :
	//if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==310) ) lNtimesReconstructedK0sMI++;
	//else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==3122) ) lNtimesReconstructedLambdaMI++;
	//else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==-3122) ) lNtimesReconstructedAntiLambdaMI++;
	//break;
	
	//}	
	//} // end loop over reconstructed V0s inside MC loop
	
	// FIXME: same comemtn for array of histos
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
	
      } // end loop ESD MC
      
    } // end ESD condition

    // FIXME: I skipped the AOD loop
    else if(fAnalysisType == "AOD") {
      lNbMCPart = mcArray->GetEntriesFast();
      lNbMCPrimary = 0;
      
      fHistMCMultiplicityTracks->Fill(lNbMCPart);
      
      for (Int_t iMc = 0; iMc < lNbMCPart; iMc++) {  
	
	// Primary vertex TO DO !!
	//
	
	AliAODMCParticle *mcAODPart = (AliAODMCParticle*)mcArray->At(iMc);
	if (!mcAODPart) {
	  //Printf("Strange analysis task (mc loop): particle with label %d not found", iMc);
	  continue;
	}
	lPdgcodeCurrentPart = mcAODPart->GetPdgCode();
	if (mcAODPart->IsPhysicalPrimary()) {lNbMCPrimary = lNbMCPrimary +1;}
	
	// Keep only K0s, Lambda and AntiLambda:
	if ( (lPdgcodeCurrentPart != 310 ) && (lPdgcodeCurrentPart != 3122 ) && (lPdgcodeCurrentPart != -3122 ) ) continue;
	
	//lEtaCurrentPart   = mcAODPart->Eta();
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
	     && (lCurrentMotherIsPrimary)
	     ) lComeFromSigma = 1;
	else lComeFromSigma = 0;
      
	//*********************************************
        // Now keep only primary particles 
	 
	// FIX IT !!!!    iMC is not defined !!!! FIX IT also in ESD/AOD loop !!
	if ( ( iMc > lNbMCPrimary ) && (!lComeFromSigma) ) continue;

	//********************************************
	// check if V0 is reconstructed several times  
	
	//lNtimesReconstructedK0s   = 0; lNtimesReconstructedLambda   = 0; lNtimesReconstructedAntiLambda   = 0;
        //lNtimesReconstructedK0sMI = 0; lNtimesReconstructedLambdaMI = 0; lNtimesReconstructedAntiLambdaMI = 0;
	    
        //for (Int_t jV0 = 0; jV0 < fESD->GetNumberOfV0s(); jV0++) {
	
	//lLabelTrackN  = 0; lLabelTrackP = 0;
	//lPartNMother  = 0; lPartPMother = 0;

	//AliAODv0    *vertexAOD= ((AliAODEvent*)fESD)->GetV0(jV0);
	//if (!vertexAOD) continue;
	//printf("enter!!");
	//AliVParticle  *trackP  = ((AliVEvent*)fESD)->GetTrack(vertexAOD->GetPosID());
	//if (!trackP) continue;
	//lLabelTrackP = TMath::Abs(trackP->GetLabel());
	//if (lLabelTrackP!=id0 && lLabelTrackP!=id1) continue;
       
	//AliVParticle  *trackN  = ((AliVEvent*)fESD)->GetTrack(vertexAOD->GetNegID());
	//if (!trackN) continue;
	//lLabelTrackN = TMath::Abs(trackN->GetLabel());
	//if (lLabelTrackN!=id0 && lLabelTrackN!=id1) continue;
	
	//AliAODMCParticle *lPartNAOD = (AliAODMCParticle*)mcArray->At(lLabelTrackN);
	//if (!lPartNAOD) continue;
	//AliAODMCParticle *lPartPAOD = (AliAODMCParticle*)mcArray->At(lLabelTrackP);
	//if (!lPartPAOD) continue;
	
	//lPartNMother = lPartNAOD->GetMother();
	//lPartPMother = lPartPAOD->GetMother();

	//lPartPMotherPDGcode = ((AliAODMCParticle*)mcArray->At(lPartPMother))->GetPdgCode();
	
	//switch (vertexAOD->GetOnFlyStatus()){
	  
	//case 0 : 
	  //if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==310) ) lNtimesReconstructedK0s++;
	  //else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==3122) ) lNtimesReconstructedLambda++;
	  //else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==-3122) ) lNtimesReconstructedAntiLambda++;
	  //break;
	  
	//case 1 :
	  //if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==310) ) lNtimesReconstructedK0sMI++;
	  //else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==3122) ) lNtimesReconstructedLambdaMI++;
	  //else if ( (lPartPMother==lPartNMother) && (lPartPMotherPDGcode==-3122) ) lNtimesReconstructedAntiLambdaMI++;
	  //break;
	  
       ///}	
      //} // end loop over reconstructed V0s inside MC loop
    
	if (TMath::Abs(lRapCurrentPart) > lCutRap) continue;
 
	if (lPdgcodeCurrentPart==310) {
	  fHistMCProdRadiusK0s->Fill(mcPosR);
	  fHistMCPtK0s->Fill(lPtCurrentPart);
	  fHistNTimesRecK0s->Fill(lNtimesReconstructedK0s);
	  fHistNTimesRecK0sVsPt->Fill(lPtCurrentPart,lNtimesReconstructedK0s);
	}
	else if (lPdgcodeCurrentPart==3122) {
	  fHistMCProdRadiusLambda->Fill(mcPosR);
	  fHistMCPtLambda->Fill(lPtCurrentPart);
	  fHistNTimesRecLambda->Fill(lNtimesReconstructedLambda);
	  fHistNTimesRecLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedLambda);
	  if (lComeFromSigma) fHistMCPtLambdaFromSigma->Fill(lPtCurrentPart);
	}
	else if (lPdgcodeCurrentPart==-3122) {
	  fHistMCProdRadiusAntiLambda->Fill(mcPosR);
	  fHistNTimesRecAntiLambda->Fill(lNtimesReconstructedAntiLambda);
	  fHistNTimesRecAntiLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedAntiLambda);
	  if (lComeFromSigma) fHistMCPtAntiLambdaFromSigma->Fill(lPtCurrentPart);
	}

      } // end loop over AODMC particles 
      fHistMCMultiplicityPrimary->Fill(lNbMCPrimary);
      
    } // end AOD condition

  } // End Loop over MC condition

  


  //************************************
  // ESD or AOD loop 
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
  Double_t lEtaK0s     = 0, lEtaLambda     = 0, lEtaAntiLambda     = 0;
  Double_t lAlphaV0      = 0, lPtArmV0       = 0;

  Double_t lPzK0s      = 0, lPzLambda      = 0, lPzAntiLambda      = 0;


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
  AliVParticle *lVPartPos   = NULL;
  AliVParticle *lVPartNeg   = NULL;

  // Daughters' momentum:
  Double_t  lMomPos[3] = {999,999,999};
  Double_t  lMomNeg[3] = {999,999,999};
  Double_t  lPtPos = 999, lPtNeg = 999;
  Double_t  lPPos = 999, lPNeg = 999;

  // Inner Wall parameters:
  Double_t  lMomInnerWallPos =999, lMomInnerWallNeg = 999;

  // AliKF Chi2 and Armenteros variables
  Double_t lChi2KFK0s  = 0, lChi2KFLambda = 0,  lChi2KFAntiLambda = 0;
  Double_t lAlphaV0K0s = 0, lAlphaV0Lambda = 0,  lAlphaV0AntiLambda = 0;
  Double_t lPtArmV0K0s = 0, lPtArmV0Lambda = 0,  lPtArmV0AntiLambda = 0;
  Double_t lQlPos   = 0, lQlNeg   = 0;


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

  AliESDtrackCuts *myTracksCuts = NULL;
     
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

    
    myPrimaryVertex = new AliAODVertex(lPrimaryVtxPosition, lPrimaryVtxCov, lPrimaryVtxChi2, NULL, -1, AliAODVertex::kPrimary);
    if (!myPrimaryVertex) return;


     // Number of Tracklets:
    //const AliMultiplicity *myMultiplicty = ((AliESDEvent*)fESD)->GetMultiplicity();
    //if (myMultiplicty->GetNumberOfTracklets() < 10) return;
    fHistTrackletPerEvent->Fill(myMultiplicty->GetNumberOfTracklets());

    lMagneticField = ((AliESDEvent*)fESD)->GetMagneticField();

    myTracksCuts = new AliESDtrackCuts();
    // require TPC refit
    myTracksCuts->SetRequireTPCRefit(kTRUE);
    // minimum number of clusters in TPC
    myTracksCuts->SetMinNClustersTPC(nbMinTPCclusters);

  }
  
  else if(fAnalysisType == "AOD") {
    printf("enter AOD!!");
    myPrimaryVertex = ((AliAODEvent*)fESD)->GetPrimaryVertex();
    if (!myPrimaryVertex) return;

    lPrimaryVtxPosition[0] = myPrimaryVertex->GetX();
    lPrimaryVtxPosition[1] = myPrimaryVertex->GetY();
    lPrimaryVtxPosition[2] = myPrimaryVertex->GetZ();

    // Cut on SPD vertex and fill histo Nevents: FIX it !

     // Tracks cuts FIX IT !

    // FIX it !!!
    lMagneticField = 999;   

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


  //************************************
  // PID

  AliESDpid *fESDpid = new AliESDpid(); // FIXME delete
  fESDpid->GetTPCResponse().SetBetheBlochParameters(fAlephParameters[0],fAlephParameters[1],fAlephParameters[2],fAlephParameters[3],fAlephParameters[4]); 
      


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
    lComeFromSigma      = -1;
    
    
    if(fAnalysisType == "ESD") {


      AliESDv0 *v0 = ((AliESDEvent*)fESD)->GetV0(iV0);
      if (!v0) continue;
      
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
	// FIXME: shouldn't this be fatal?
	Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
	continue;
      }
      // Remove like-sign
      if ( myTrackPosTest->GetSign() == myTrackNegTest->GetSign()){
	// FIXME: how can this happen?
	continue;
      } 
     
      //    FIXME: are the GetParamN/GetParamP reliable? If so, why do you need to check the sign below?
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

      // Inner Wall parameter (used for pid):
      const AliExternalTrackParam *myInnerWallTrackPos = myTrackPos->GetInnerParam(); 
      if(myInnerWallTrackPos) lMomInnerWallPos = myInnerWallTrackPos->GetP(); 
      const AliExternalTrackParam *myInnerWallTrackNeg = myTrackNeg->GetInnerParam(); 
      if(myInnerWallTrackNeg) lMomInnerWallNeg = myInnerWallTrackNeg->GetP(); 
	      
      // DCA between daughter and Primary Vertex:
      if (myTrackPos) lDcaPosToPrimVertex = TMath::Abs(myTrackPos->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
      
      if (myTrackNeg) lDcaNegToPrimVertex = TMath::Abs(myTrackNeg->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lMagneticField) );
      
      // Quality tracks cuts:
      if ( !(myTracksCuts->IsSelected(myTrackPos)) || !(myTracksCuts->IsSelected(myTrackNeg)) ) {
	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;  
	continue;
      }
      // Armenteros variables:
      lAlphaV0      =  v0->AlphaV0();
      lPtArmV0      =  v0->PtArmV0();

      // Pseudorapidity:
      lV0Eta = v0->Eta();
      
      // PID
      if (fUsePID.Contains("withPID")) {
	nSigmaPosPion   = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackPos,AliPID::kPion));
	
	nSigmaNegPion   = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackNeg,AliPID::kPion));

	nSigmaPosProton = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackPos,AliPID::kProton));
	
	nSigmaNegProton = TMath::Abs(fESDpid->NumberOfSigmasTPC(myTrackNeg,AliPID::kProton));
      }
      else {
	nSigmaPosPion = 0; nSigmaNegPion =0; nSigmaPosProton = 0; nSigmaNegProton= 0;
      }
      
      // Monte-Carlo particle associated to reconstructed particles: 
      if (fAnalysisMC) {
	//if (lLabelTrackPos < 0 || lLabelTrackNeg < 0) continue;
	TParticle  *lMCESDPartPos  = stack->Particle(lLabelTrackPos);
	if(!lMCESDPartPos) { 
	  Printf("no MC particle for positive and/or negative daughter\n");
	  
        if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;  
	continue;
    
	}
	TParticle  *lMCESDPartNeg  = stack->Particle(lLabelTrackNeg);
	if (!lMCESDPartNeg) {
        if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;  
	continue;
    }
	lPDGCodePosDaughter = lMCESDPartPos->GetPdgCode();
	lPDGCodeNegDaughter = lMCESDPartNeg->GetPdgCode();
	lIndexPosMother = lMCESDPartPos->GetFirstMother();
	lIndexNegMother = lMCESDPartNeg->GetFirstMother();
	
	if (lIndexPosMother == -1) {
        if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;  
	continue;
    }
	TParticle  *lMCESDMother    = stack->Particle(lIndexPosMother);
	if (!lMCESDMother) {
        if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;  
	continue;
    }
	lPdgcodeMother         = lMCESDMother->GetPdgCode();
	lIndexMotherOfMother   = lMCESDMother->GetFirstMother();
	if (lIndexMotherOfMother ==-1) lPdgcodeMotherOfMother = 0;
	else {
	  TParticle  *lMCESDMotherOfMother    = stack->Particle(lIndexMotherOfMother);
	  if (!lMCESDMotherOfMother) {
        if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;  
	continue;
    }
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

    } // end ESD condition


    
    else if(fAnalysisType == "AOD") { 

      AliAODv0     *myAODv0 = ((AliAODEvent*)fESD)->GetV0(iV0);
      if (!myAODv0) continue;

      // Primary vertex:
      fHistPrimaryVertexPosXV0events->Fill(lPrimaryVtxPosition[0]);
      fHistPrimaryVertexPosYV0events->Fill(lPrimaryVtxPosition[1]);
      fHistPrimaryVertexPosZV0events->Fill(lPrimaryVtxPosition[2]);
      

      //Multiplicity:
      if(lOnFlyStatus == fUseOnTheFly) nv0s++;

      // V0's Daughters
      lIndexTrackPos = TMath::Abs(myAODv0->GetPosID());
      lIndexTrackNeg = TMath::Abs(myAODv0->GetNegID());
      
      AliVParticle  *lVPartPosTest  = ((AliVEvent*)fESD)->GetTrack(lIndexTrackPos);
      AliVParticle  *lVPartNegTest  = ((AliVEvent*)fESD)->GetTrack(lIndexTrackNeg);
      //AliAODTrack  *lVPartPos  = ((AliAODEvent*)fESD)->GetTrack(lIndexTrackPos);
      //AliAODTrack  *lVPartNeg  = ((AliAODEvent*)fESD)->GetTrack(lIndexTrackNeg);

      if (!lVPartPosTest ||(!lVPartNegTest )) {
	Printf("strange analysis::UserExec:: Could not retreive one of the daughter track\n");
	continue;
      }

      // Quality cuts:
      // TO DO !!!!!!!

      // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
      //if( !(lVPartPosTest->GetStatus() & AliAODTrack::kTPCrefit)) continue;      
      //if( !(lVPartNegTest->GetStatus() & AliAODTrack::kTPCrefit)) continue;
      

      lDcaPosToPrimVertex = myAODv0->DcaPosToPrimVertex();	
      lDcaNegToPrimVertex = myAODv0->DcaNegToPrimVertex();
      lOnFlyStatus        = myAODv0->GetOnFlyStatus();
      lChi2V0             = myAODv0->Chi2V0();
      lDcaV0Daughters     = myAODv0->DcaV0Daughters();
      lDcaV0ToPrimVertex  = myAODv0->DcaV0ToPrimVertex();
      lV0DecayLength      = myAODv0->DecayLengthV0(lPrimaryVtxPosition);
      lV0cosPointAngle    = myAODv0->CosPointingAngle(lPrimaryVtxPosition);
      lV0Radius           = myAODv0->RadiusV0();

      if( lVPartPosTest->Charge() ==1){
	
	lVPartPos = ((AliVEvent*)fESD)->GetTrack(lIndexTrackPos);
	lVPartNeg = ((AliVEvent*)fESD)->GetTrack(lIndexTrackNeg);
	
	
	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF; posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;
	
	//negPiKF = new AliKFParticle( *(myAODv0->GetParamN()) ,-211);
	//posPiKF = new AliKFParticle( *(myAODv0->GetParamP()) ,211);
	//posPKF  = new AliKFParticle( *(myAODv0->GetParamP()) ,2212);
	//negAPKF = new AliKFParticle( *(myAODv0->GetParamN()) ,-2212);
	// TO DO !!!!!!
	negPiKF = NULL;
	posPiKF = NULL;
	posPKF  = NULL;
	negAPKF = NULL;
	
      }
           
      if( lVPartPosTest->Charge() ==-1){
	
	lVPartPos = ((AliVEvent*)fESD)->GetTrack(lIndexTrackNeg);
	lVPartNeg = ((AliVEvent*)fESD)->GetTrack(lIndexTrackPos);
	
	if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF; posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;
	
	//negPiKF = new AliKFParticle( *(myAODv0->GetParamP()) ,-211);
	//posPiKF = new AliKFParticle( *(myAODv0->GetParamN()) ,211);
	//posPKF  = new AliKFParticle( *(myAODv0->GetParamN()) ,2212);
	//negAPKF = new AliKFParticle( *(myAODv0->GetParamP()) ,-2212);
	negPiKF = NULL;
	posPiKF = NULL;
	posPKF  = NULL;
	negAPKF = NULL;
      }

      lLabelTrackPos  = TMath::Abs(lVPartPos->GetLabel());
      lLabelTrackNeg  = TMath::Abs(lVPartNeg->GetLabel());
      
      // Armenteros variables:
      lAlphaV0   = myAODv0->AlphaV0();
      lPtArmV0   = myAODv0->PtArmV0();

      // Pseudorapidity:
      lV0Eta = myAODv0->PseudoRapV0();
      
      // PID not accessible with AOD !
      nSigmaPosPion = 0; nSigmaNegPion =0; nSigmaPosProton = 0; nSigmaNegProton= 0;

      
      // Monte-Carlo particle associated to reconstructed particles:  
      if (fAnalysisMC) {
	AliAODMCParticle *lMCAODPartPos = (AliAODMCParticle*)mcArray->At(lLabelTrackPos);
	if (!lMCAODPartPos) {
        if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;  
	continue;
    }
	AliAODMCParticle *lMCAODPartNeg = (AliAODMCParticle*)mcArray->At(lLabelTrackNeg);
	if(!lMCAODPartNeg)  
	 // Printf("strange analysis::UserExec:no MC particle for negative daughter\n");
	  {
        if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;  
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
      }
            
    } // end AOD condition
    
    
    // Multiplicity:
    if(lOnFlyStatus == fUseOnTheFly) nv0s++;

    // Daughter momentum cut: ! FIX it in case of AOD !
    if ( (lPtPos  < cutMinPtDaughter ) ||
         (lPtNeg  < cutMinPtDaughter )
	 ) {
      if (negPiKF) delete negPiKF; negPiKF=NULL;
      if (posPiKF) delete posPiKF; posPiKF=NULL;
      if (posPKF)  delete posPKF;  posPKF=NULL;
      if (negAPKF) delete negAPKF; negAPKF=NULL;        
      continue;
    }
    
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
    
    // Invariant mass
    lInvMassK0s        = v0K0sKF.GetMass();
    lInvMassLambda     = v0LambdaKF.GetMass();
    lInvMassAntiLambda = v0AntiLambdaKF.GetMass();
    
    // Rapidity:
    lRapK0s        = 0.5*TMath::Log((v0K0sKF.GetE()+v0K0sKF.GetPz())/(v0K0sKF.GetE()-v0K0sKF.GetPz()+1.e-13));
    lRapLambda     = 0.5*TMath::Log((v0LambdaKF.GetE()+v0LambdaKF.GetPz())/(v0LambdaKF.GetE()-v0LambdaKF.GetPz()+1.e-13));
    lRapAntiLambda = 0.5*TMath::Log((v0AntiLambdaKF.GetE()+v0AntiLambdaKF.GetPz())/(v0AntiLambdaKF.GetE()-v0AntiLambdaKF.GetPz()+1.e-13));

    // Pseudo-rapidity
    lEtaK0s     = v0K0sKF.GetEta();
    lEtaLambda  = v0LambdaKF.GetEta();
    lEtaAntiLambda  = v0AntiLambdaKF.GetEta();

    // Pz:
    lPzK0s        = v0K0sKF.GetPz();
    lPzLambda     = v0LambdaKF.GetPz();
    lPzAntiLambda = v0AntiLambdaKF.GetPz();
    
    // Pt:
    lPtK0s        = v0K0sKF.GetPt();
    lPtLambda     = v0LambdaKF.GetPt();
    lPtAntiLambda = v0AntiLambdaKF.GetPt();

    if (lPtK0s==0) {
        if (negPiKF) delete negPiKF; negPiKF=NULL;
	if (posPiKF) delete posPiKF; posPiKF=NULL;
	if (posPKF)  delete posPKF;  posPKF=NULL;
	if (negAPKF) delete negAPKF; negAPKF=NULL;  
	// FIXME: should you really continue here and below?
	continue;
    }
    if (lPtLambda==0) {
      if (negPiKF) delete negPiKF; negPiKF=NULL;
      if (posPiKF) delete posPiKF; posPiKF=NULL;
      if (posPKF)  delete posPKF;  posPKF=NULL;
      if (negAPKF) delete negAPKF; negAPKF=NULL;  
      continue;
    }
    // Pt Resolution
    deltaPtK0s        = (lPtK0s - mcMotherPt)/mcMotherPt;
    deltaPtLambda     = (lPtLambda - mcMotherPt)/mcMotherPt;
    deltaPtAntiLambda = (lPtAntiLambda - mcMotherPt)/mcMotherPt;

    // KF Chi2
    lChi2KFK0s        = v0K0sKF.GetChi2();
    lChi2KFLambda     = v0LambdaKF.GetChi2();
    lChi2KFAntiLambda = v0AntiLambdaKF.GetChi2();
    
    // Reconstructed Position
    rcPosXK0s = v0K0sKF.GetX();
    rcPosYK0s = v0K0sKF.GetY(); 
    rcPosZK0s = v0K0sKF.GetZ();
    rcPosRK0s = TMath::Sqrt(rcPosXK0s*rcPosXK0s+rcPosYK0s*rcPosYK0s);

    rcPosXLambda = v0LambdaKF.GetX(); 
    rcPosYLambda = v0LambdaKF.GetY(); 
    rcPosZLambda = v0LambdaKF.GetZ();
    rcPosRLambda = TMath::Sqrt(rcPosXLambda*rcPosXLambda+rcPosYLambda*rcPosYLambda); 

    rcPosXAntiLambda = v0AntiLambdaKF.GetX();
    rcPosYAntiLambda = v0AntiLambdaKF.GetY(); 
    rcPosZAntiLambda = v0AntiLambdaKF.GetZ();
    rcPosRAntiLambda = TMath::Sqrt(rcPosXAntiLambda*rcPosXAntiLambda+rcPosYAntiLambda*rcPosYAntiLambda); 

    TVector3 momPos(lMomPos[0],lMomPos[1],lMomPos[2]);
    TVector3 momNeg(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
    TVector3 momTotK0s(v0K0sKF.GetPx(),v0K0sKF.GetPy(),v0K0sKF.GetPz());
    TVector3 momTotLambda(v0LambdaKF.GetPx(),v0LambdaKF.GetPy(),v0LambdaKF.GetPz());
    TVector3 momTotAntiLambda(v0AntiLambdaKF.GetPx(),v0AntiLambdaKF.GetPy(),v0AntiLambdaKF.GetPz());
    
    lQlPos = momPos.Dot(momTotK0s)/momTotK0s.Mag();
    lQlNeg = momNeg.Dot(momTotK0s)/momTotK0s.Mag();
    lAlphaV0K0s = 1.-2./(1.+lQlPos/lQlNeg);
    lQlPos = momPos.Dot(momTotLambda)/momTotLambda.Mag();
    lQlNeg = momNeg.Dot(momTotLambda)/momTotLambda.Mag();
    lAlphaV0Lambda = 1.-2./(1.+lQlPos/lQlNeg);
    lQlPos = momPos.Dot(momTotAntiLambda)/momTotAntiLambda.Mag();
    lQlNeg = momNeg.Dot(momTotAntiLambda)/momTotAntiLambda.Mag();
    lAlphaV0AntiLambda = 1.-2./(1.+lQlPos/lQlNeg);
    
    lPtArmV0K0s = momPos.Perp(momTotK0s);
    lPtArmV0Lambda = momPos.Perp(momTotLambda);
    lPtArmV0AntiLambda = momPos.Perp(momTotAntiLambda);
    
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
      else if ( (lPDGCodePosDaughter==11) &&
		(lPDGCodeNegDaughter==-11) &&
		(lPdgcodeMother==22 ) )
	lCheckGamma = 1;
    } // end "look for associated particles  
   
    
    // Cuts:
/*    if (fUseCut.Contains("yes")) {
      if ( (lDcaPosToPrimVertex < 0.036 ) ||
	   (lDcaNegToPrimVertex < 0.036 ) ||
	   (lDcaV0Daughters     > 0.5   ) ||
	   (lV0cosPointAngle    < 0.999 ) 
	   )	
	continue;
    }
*/

/*
      if ( (lDcaV0Daughters     > 0.3   ) ||
	   (lV0cosPointAngle    < 0.998 )

	   )	continue;
*/
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
    
    
    //*****************************
    // filling histograms

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
    fHistChi2KFBeforeCutK0s->Fill(lChi2KFK0s,lOnFlyStatus); 
    fHistChi2KFBeforeCutLambda->Fill(lChi2KFLambda,lOnFlyStatus);
    fHistChi2KFBeforeCutAntiLambda->Fill(lChi2KFAntiLambda,lOnFlyStatus);


    // Histo versus Rap and armenteros plot
    if (lOnFlyStatus == fUseOnTheFly){
      if (lCheckMcK0Short) fHistAsMcRapK0->Fill(lRapK0s);
      if (lCheckMcLambda) fHistAsMcRapLambda->Fill(lRapLambda);
      if (lCheckMcAntiLambda) fHistAsMcRapLambda->Fill(lRapAntiLambda);
//      fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0);
//      fHistDaughterPt->Fill(lPtPos,lPtNeg);
    }
    
    // FIXME: associated histos, what are they used for?

    // K0s associated histograms in |rap| < lCutRap:

////////////////////////////
    if ( lCheckPIDK0sPosDaughter && lCheckPIDK0sNegDaughter
	 && (lChi2KFK0s < cutChi2KF)) fHistPzPtBeforeK0s->Fill(TMath::Abs(lPzK0s/lPtK0s));
/////////////////////////////

    if ( lCheckPIDK0sPosDaughter && lCheckPIDK0sNegDaughter
	 && (lChi2KFK0s < cutChi2KF) && (TMath::Abs(lPzK0s/lPtK0s)<0.7) )     {



	fHistPzPtAfterK0s->Fill(TMath::Abs(lPzK0s/lPtK0s));


      
      fHistChi2KFAfterCutK0s->Fill(lChi2KFK0s,lOnFlyStatus);

      if (TMath::Abs(lRapK0s) < lCutRap) {

	fHistNsigmaPosPionK0->Fill(nSigmaPosPion);
	fHistNsigmaNegPionK0->Fill(nSigmaNegPion);
	
	if (lOnFlyStatus == fUseOnTheFly){
	  fHistMassK0->Fill(lInvMassK0s);
	  fHistMassVsRadiusK0->Fill(rcPosRK0s,lInvMassK0s);
	  fHistPtVsMassK0->Fill(lInvMassK0s,lPtK0s);


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
	}
      } // end rapidity condition
    } // end nsigma condition
    

    // Associated Lambda histograms in |rap| < lCutRap

////////////////////////len koly kontrole Abs(Pz/Pt)
    if ( lCheckPIDLambdaPosDaughter && lCheckPIDLambdaNegDaughter
	 && (lChi2KFLambda < cutChi2KF)) fHistPzPtBeforeLambda->Fill(TMath::Abs(lPzLambda/lPtLambda)); 
////////////////////////

    if ( lCheckPIDLambdaPosDaughter && lCheckPIDLambdaNegDaughter
	 && (lChi2KFLambda < cutChi2KF) && (TMath::Abs(lPzLambda/lPtLambda)<0.7) )  {

 


	fHistPzPtAfterLambda->Fill(TMath::Abs(lPzLambda/lPtLambda));      

      fHistChi2KFAfterCutLambda->Fill(lChi2KFLambda,lOnFlyStatus);

      if (TMath::Abs(lRapLambda) < lCutRap) {

	fHistNsigmaPosProtonLambda->Fill(nSigmaPosProton);
	fHistNsigmaNegPionLambda->Fill(nSigmaNegPion);
	if (lOnFlyStatus == fUseOnTheFly){
	  fHistMassLambda->Fill(lInvMassLambda);
	  fHistMassVsRadiusLambda->Fill(rcPosRLambda,lInvMassLambda);
	  fHistPtVsMassLambda->Fill(lInvMassLambda,lPtLambda);

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
	}
      } // end rapidity condition
    } //end nsigma condition - lambda


    if (negPiKF) delete negPiKF; negPiKF= NULL;
    if (posPiKF) delete posPiKF; posPiKF= NULL;
    if (posPKF)  delete posPKF;  posPKF = NULL;
    if (negAPKF) delete negAPKF; negAPKF= NULL;
    
  } // end V0 loop

  fHistV0Multiplicity->Fill(nv0s);
  if (fAnalysisType == "ESD") { if(myPrimaryVertex) delete myPrimaryVertex; }

  if(myTracksCuts) delete myTracksCuts;
  
  // Post output data
  PostData(1, fListHist);
  PostData(2, fCentrSelector);
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

