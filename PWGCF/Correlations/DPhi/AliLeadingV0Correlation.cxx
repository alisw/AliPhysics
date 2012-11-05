//
//            Class for Leading Charged Track+V0 Correlations Analysis

#include <TROOT.h>
#include <TList.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TRandom.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>


#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODVertex.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliEventPoolManager.h"
#include "AliCentrality.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliPhysicsSelectionTask.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliVParticle.h"
#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliStack.h"
#include "AliAODMCHeader.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliExternalTrackParam.h"

#include "AliLeadingV0Correlation.h"

#define CorrBinsX 26
#define CorrBinsY 24

ClassImp(AliLeadingV0Correlation)
ClassImp(AliLeadingBasicParticle)

//---------------------------------------------------------------------------------------
AliLeadingV0Correlation::AliLeadingV0Correlation()
   : AliAnalysisTaskSE(),
  fAODEvent(0x0),
  fPoolMgr(0x0),         
  fPoolK0(0x0),
  fPoolLambda(0x0),
  fPIDResponse(0x0),
  fPoolMaxNEvents(0), 
  fPoolMinNTracks(0), 
  fMinEventsToMix(0),
  fNzVtxBins(0), 
  fNCentBins(0),
  fcollidingSys("PbPb"),
  fpvzcut(0),
  fTrackEtaCut(0),
  fFilterBit(128),
  ftrigPtLow(0),
  ftrigPtHigh(0),
  fassocPtLow(0),
  fassocPtHigh(0),
  fAnalysisMC(0),
  fUsePID(""),
  fRapidityCut(0.),
  fCutV0Radius(0.),
  fCutDCANegToPV(0.),
  fCutDCAPosToPV(0.),
  fCutDCAV0Daughters(0.),
  fCutV0CosPA(0.),
  fSpecialArmenterosCutK0s(0.),
  fPerformMultiplicityStudy(0), 
  fLoMultBound(0),
  fHiMultBound(0),
  fOutputList(0),
  fHistMCPrimaryVertexX(0),
  fHistMCPrimaryVertexY(0),
  fHistMCPrimaryVertexZ(0),
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
  fHistMCRapLambda(0),
  fHistMCRapAntiLambda(0),
  fHistMCPtK0s(0),
  fHistMCPtLambda(0),
  fHistMCPtAntiLambda(0),
  fHistMCPtLambdaFromSigma(0),
  fHistMCPtAntiLambdaFromSigma(0),
  fHistNTimesRecK0s(0),
  fHistNTimesRecLambda(0),
  fHistNTimesRecAntiLambda(0),
  fHistNTimesRecK0sVsPt(0),
  fHistNTimesRecLambdaVsPt(0),
  fHistNTimesRecAntiLambdaVsPt(0),
  fHistMCDaughterTrack(0),
  fHistPrimRawPtVsYK0s(0),
  fHistPrimRawPtVsYLambda(0),
  fHistPrimRawPtVsYAntiLambda(0),
  fHistPrimaryVertexX(0),
  fHistPrimaryVertexY(0),
  fHistPrimaryVertexZ(0),
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
  fHistDcaPosToPrimVertexL(0),
  fHistDcaNegToPrimVertexL(0),
  fHistRadiusV0L(0),
  fHistDecayLengthV0L(0),
  fHistDcaV0DaughtersL(0),
  fHistChi2L(0),
  fHistCosPointAngleL(0),
  fHistcTauL(0),
  fHistDcaPosToPrimVertexLvsMassL(0),
  fHistDcaNegToPrimVertexLvsMassL(0),
  fHistRadiusV0LvsMassL(0),
  fHistDecayLengthV0LvsMassL(0),
  fHistDcaV0DaughtersLvsMassL(0),
  fHistCosPointAngleLvsMassL(0),
  fHistCosPointAngleLVsMassVsPtsigL(0),
  fHistCosPointAngleLVsMassVsPtbackL(0),
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
  fHistMassK0(0),
  fHistMassLambda(0),
  fHistMassAntiLambda(0),
  fHistMassVsRadiusK0(0),
  fHistMassVsRadiusLambda(0),
  fHistMassVsRadiusAntiLambda(0),
  fHistPtVsMassK0(0),
  fHistPtVsMassLambda(0),
  fHistPtVsMassAntiLambda(0),
  fHistArmenterosPodolanski(0),
  fHistNsigmaPosProtonLambda(0),
  fHistNsigmaNegPionLambda(0),
  fHistNsigmaPosProtonAntiLambda(0),
  fHistNsigmaNegPionAntiLambda(0),
  fHistNsigmaPosPionK0(0),
	fHistNsigmaNegPionK0(0),
	fHistAsMcPtK0(0),
	fHistAsMcPtLambda(0),
	fHistAsMcPtAntiLambda(0),
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
	fHistAsMcSecondaryPtLambdaFromSigma(0),
	fHistAsMcSecondaryPtAntiLambdaFromSigma(0),
	fHistSibK0(0),
	fHistMixK0(0),
	fHistSibLambda(0),
	fHistMixLambda(0),
	fHistSibK0MC(0),
	fHistMixK0MC(0),
	fHistSibLambdaMC(0),
	fHistMixLambdaMC(0)


{
	fRapidityCut=0.75;
	fUsePID="withPID";
	fCutV0Radius=0.5;
	fCutDCANegToPV=0.06;
	fCutDCAPosToPV=0.06;
	fCutDCAV0Daughters=1.0;
	fCutV0CosPA=0.995;
	fSpecialArmenterosCutK0s=5;

}
//---------------------------------------------------------------------------------------
AliLeadingV0Correlation::AliLeadingV0Correlation(const char *name)
   : AliAnalysisTaskSE(name),
	fAODEvent(0x0),
	fPoolMgr(0x0),         
	fPoolK0(0x0),
	fPoolLambda(0x0),
    fPIDResponse(0x0),
	fPoolMaxNEvents(0), 
	fPoolMinNTracks(0), 
	fMinEventsToMix(0),
	fNzVtxBins(0), 
	fNCentBins(0),
	fcollidingSys("PbPb"),
	fpvzcut(0),
    fTrackEtaCut(0),
    fFilterBit(128),
  ftrigPtLow(0),
  ftrigPtHigh(0),
  fassocPtLow(0),
  fassocPtHigh(0),
  fAnalysisMC(0),
  fUsePID(""),
  fRapidityCut(0.),
  fCutV0Radius(0.),
  fCutDCANegToPV(0.),
  fCutDCAPosToPV(0.),
  fCutDCAV0Daughters(0.),
  fCutV0CosPA(0.),
  fSpecialArmenterosCutK0s(0.),
  fPerformMultiplicityStudy(0), 
  fLoMultBound(0),
  fHiMultBound(0),
  fOutputList(0),
	fHistMCPrimaryVertexX(0),
	fHistMCPrimaryVertexY(0),
	fHistMCPrimaryVertexZ(0),
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
	fHistMCRapLambda(0),
	fHistMCRapAntiLambda(0),
	fHistMCPtK0s(0),
	fHistMCPtLambda(0),
	fHistMCPtAntiLambda(0),
	fHistMCPtLambdaFromSigma(0),
	fHistMCPtAntiLambdaFromSigma(0),
	fHistNTimesRecK0s(0),
	fHistNTimesRecLambda(0),
	fHistNTimesRecAntiLambda(0),
	fHistNTimesRecK0sVsPt(0),
	fHistNTimesRecLambdaVsPt(0),
	fHistNTimesRecAntiLambdaVsPt(0),
	fHistMCDaughterTrack(0),
	fHistPrimRawPtVsYK0s(0),
	fHistPrimRawPtVsYLambda(0),
	fHistPrimRawPtVsYAntiLambda(0),
	fHistPrimaryVertexX(0),
	fHistPrimaryVertexY(0),
	fHistPrimaryVertexZ(0),
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
	fHistDcaPosToPrimVertexL(0),
	fHistDcaNegToPrimVertexL(0),
	fHistRadiusV0L(0),
	fHistDecayLengthV0L(0),
	fHistDcaV0DaughtersL(0),
	fHistChi2L(0),
	fHistCosPointAngleL(0),
	fHistcTauL(0),
	fHistDcaPosToPrimVertexLvsMassL(0),
	fHistDcaNegToPrimVertexLvsMassL(0),
	fHistRadiusV0LvsMassL(0),
	fHistDecayLengthV0LvsMassL(0),
	fHistDcaV0DaughtersLvsMassL(0),
	fHistCosPointAngleLvsMassL(0),
	fHistCosPointAngleLVsMassVsPtsigL(0),
	fHistCosPointAngleLVsMassVsPtbackL(0),
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
	fHistMassK0(0),
	fHistMassLambda(0),
	fHistMassAntiLambda(0),
	fHistMassVsRadiusK0(0),
	fHistMassVsRadiusLambda(0),
	fHistMassVsRadiusAntiLambda(0),
	fHistPtVsMassK0(0),
	fHistPtVsMassLambda(0),
	fHistPtVsMassAntiLambda(0),
	fHistArmenterosPodolanski(0),
	fHistNsigmaPosProtonLambda(0),
	fHistNsigmaNegPionLambda(0),
	fHistNsigmaPosProtonAntiLambda(0),
	fHistNsigmaNegPionAntiLambda(0),
	fHistNsigmaPosPionK0(0),
	fHistNsigmaNegPionK0(0),
	fHistAsMcPtK0(0),
	fHistAsMcPtLambda(0),
	fHistAsMcPtAntiLambda(0),
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
	fHistAsMcSecondaryPtLambdaFromSigma(0),
	fHistAsMcSecondaryPtAntiLambdaFromSigma(0),
	fHistSibK0(0),
	fHistMixK0(0),
	fHistSibLambda(0),
	fHistMixLambda(0),
	fHistSibK0MC(0),
	fHistMixK0MC(0),
	fHistSibLambdaMC(0),
	fHistMixLambdaMC(0)
{	
	fRapidityCut=0.75;
	fUsePID="withPID";
	fCutV0Radius=0.5;
	fCutDCANegToPV=0.06;
	fCutDCAPosToPV=0.06;
	fCutDCAV0Daughters=1.0;
	fCutV0CosPA=0.995;
	fSpecialArmenterosCutK0s=5;
	
   DefineOutput(1, TList::Class());                                            
}

//---------------------------------------------------------------------------------------
AliLeadingV0Correlation::~AliLeadingV0Correlation()
{
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if (fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutputList;
   }
}
//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::UserCreateOutputObjects()
{	
   fOutputList = new TList();
   fOutputList->SetOwner();
	
	Int_t lCustomNBins = 200; 
	Double_t lCustomPtUpperLimit = 20; 
	//Int_t lCustomNBinsMultiplicity = 100;
	
	
	//---------------------------------------------- MC histograms -----------------------------------------------------//
	
	// Primary Vertex:
	fHistMCPrimaryVertexX          = new TH1F("h1MCPrimaryVertexX", "MC Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
	fOutputList->Add(fHistMCPrimaryVertexX);
	
	fHistMCPrimaryVertexY          = new TH1F("h1MCPrimaryVertexY", "MC Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
	fOutputList->Add(fHistMCPrimaryVertexY);
	
	fHistMCPrimaryVertexZ          = new TH1F("h1MCPrimaryVertexZ", "MC Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
	fOutputList->Add(fHistMCPrimaryVertexZ);
	
	
	// Production Radius of non-primary particles:
	fHistMCtracksProdRadiusK0s           = new TH2F("h2MCtracksProdRadiusK0s","Non-primary MC K^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
	fOutputList->Add(fHistMCtracksProdRadiusK0s);
	
	fHistMCtracksProdRadiusLambda        = new TH2F("h2MCtracksProdRadiusLambda","Non-primary MC #Lambda^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
	fOutputList->Add(fHistMCtracksProdRadiusLambda);
	
	fHistMCtracksProdRadiusAntiLambda    = new TH2F("h2MCtracksProdRadiusAntiLambda","Non-primary MC #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
	fOutputList->Add(fHistMCtracksProdRadiusAntiLambda);
	
	// Decay Radius of non-primary particles:
	fHistMCtracksDecayRadiusK0s          = new TH1F("h1MCtracksDecayRadiusK0s","Non-primary MC K^{0} Decay Radius;r (cm)",101,-1,100);
	fOutputList->Add(fHistMCtracksDecayRadiusK0s);
	
	fHistMCtracksDecayRadiusLambda       = new TH1F("h1MCtracksDecayRadiusLambda","Non-primary MC #Lambda^{0} Decay Radius;r (cm)",101,-1,100);
	fOutputList->Add(fHistMCtracksDecayRadiusLambda);
	
	fHistMCtracksDecayRadiusAntiLambda   = new TH1F("h1MCtracksDecayRadiusAntiLambda","Non-primary #bar{#Lambda}^{0} Decay Radius;r (cm)",100,1,101);
	fOutputList->Add(fHistMCtracksDecayRadiusAntiLambda);
	
	// Rapidity distribution:
	fHistMCRapK0s                 = new TH1F("h1MCRapK0s", "K^{0};y",160,-4,4);
	fOutputList->Add(fHistMCRapK0s);
	
	fHistMCRapLambda              = new TH1F("h1MCRapLambda", "#Lambda;y",160,-4,4);
	fOutputList->Add(fHistMCRapLambda);
	
	
	fHistMCRapAntiLambda          = new TH1F("h1MCRapAntiLambda", "#bar{#Lambda};y",160,-4,4);
	fOutputList->Add(fHistMCRapAntiLambda);
	
	
	// Production Radius
	fHistMCProdRadiusK0s                 = new TH1F("h1MCProdRadiusK0s", "MC K^{0} Production Radius;r (cm);Count", 400, -2, 2);
	fOutputList->Add(fHistMCProdRadiusK0s);
	
	fHistMCProdRadiusLambda              = new TH1F("h1MCProdRadiusLambda", "MC #Lambda^{0} Production Radius;r (cm);Count", 400, -2, 2);
	fOutputList->Add(fHistMCProdRadiusLambda);
	
	fHistMCProdRadiusAntiLambda         = new TH1F("h1MCProdRadiusAntiLambda", "MC #bar{#Lambda}^{0} Production Radius;r (cm);Count", 400, -2, 2);
	fOutputList->Add(fHistMCProdRadiusAntiLambda);
	
	
	// Pt distribution:
	fHistMCPtK0s               = new TH1F("h1MCPtK0s", "K^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtK0s);
	
	fHistMCPtLambda            = new TH1F("h1MCPtLambda", "#Lambda^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtLambda);
	
	fHistMCPtAntiLambda            = new TH1F("h1MCPtAntiLambda", "#AntiLambda^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtAntiLambda);
	
	// Pt distribution of Lambda coming from Sigma decay
	fHistMCPtLambdaFromSigma      = new TH1F("h1MCPtLambdaFromSigma", "#Lambda^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtLambdaFromSigma);
	
	fHistMCPtAntiLambdaFromSigma  = new TH1F("h1MCPtAntiLambdaFromSigma", "#Lambda^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtAntiLambdaFromSigma);
	
	// Multiple reconstruction studies
	fHistNTimesRecK0s             = new TH1F("h1NTimesRecK0s","number of times a K0s is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
	fOutputList->Add(fHistNTimesRecK0s);
	
	fHistNTimesRecLambda          = new TH1F("h1NTimesRecLambda","number of times a Lambda is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
	fOutputList->Add(fHistNTimesRecLambda);
	
	fHistNTimesRecAntiLambda      = new TH1F("h1NTimesRecAntiLambda","number of times an AntiLambda is reconstructed in -1<y<1;number of times;counts",500,-0.5,4.5);
	fOutputList->Add(fHistNTimesRecAntiLambda);
	
	fHistNTimesRecK0sVsPt         = new TH2F("h2NTimesRecK0sVsPt","NTimes versus Pt, K^{0} in -1<y<1;p{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
	fOutputList->Add(fHistNTimesRecK0sVsPt);
	
	fHistNTimesRecLambdaVsPt      = new TH2F("h2NTimesRecLambdaVsPt","NTimes versus Pt, #Lambda^{0} in -1<y<1;p{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
	fOutputList->Add(fHistNTimesRecLambdaVsPt);
	
	fHistNTimesRecAntiLambdaVsPt  = new TH2F("h2NTimesRecAntiLambdaVsPt","NTimes versus Pt, #bar{#Lambda}^{0} in -1<y<1;p{t} (GeV/c);number of times",75,0,15,5,-0.5,4.5);
	fOutputList->Add(fHistNTimesRecAntiLambdaVsPt);
	
	// Pt Distribution of non-primary particles:
	fHistMCPtAllK0s                      = new TH1F("h1MCPtAllK0s", "Non-primary MC K^{0};p{t} (GeV/c);Counts",240,0,12);
	fOutputList->Add(fHistMCPtAllK0s);
	
	fHistMCPtAllLambda                   = new TH1F("h1MCPtAllLambda", "Non-primary MC #Lambda^{0};p{t} (GeV/c);Counts",240,0,12);
	fOutputList->Add(fHistMCPtAllLambda);
	
	fHistMCPtAllAntiLambda               = new TH1F("h1MCPtAllAntiLambda", "Non-primary MC #bar{#Lambda}^{0};p{t} (GeV/c);Counts",240,0,12);
	fOutputList->Add(fHistMCPtAllAntiLambda);

	fHistMCDaughterTrack         = new TH1F("h1MCDaughterTrack","Distribution of mc id for daughters;id tags;Counts",15,0,15);
	fOutputList->Add(fHistMCDaughterTrack);
	
	fHistPrimRawPtVsYK0s = new TH2F( "f3dHistPrimRawPtVsYVsMultK0Short", "Pt{K0S} Vs Y{K0S} Vs Multiplicity; Pt{K0S} (GeV/c); Y{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2);
	fOutputList->Add(fHistPrimRawPtVsYK0s);

	
	fHistPrimRawPtVsYLambda = new TH2F( "f3dHistPrimRawPtVsYVsMultLambda", "Pt{lambda} Vs Y{#Lambda} Vs Multiplicity; Pt{lambda} (GeV/c); Y{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2);
	fOutputList->Add(fHistPrimRawPtVsYLambda);
	

	fHistPrimRawPtVsYAntiLambda = new TH2F( "f3dHistPrimRawPtVsYVsMultAntiLambda", "Pt{antilambda} Vs Y{#Lambda} Vs Multiplicity; Pt{antilambda} (GeV/c); Y{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2);
	fOutputList->Add(fHistPrimRawPtVsYAntiLambda);
	
	//---------------------------------------------End Of MC Histos-----------------------------------------------------//
	
	// Primary Vertex:
	
	fHistPrimaryVertexX          = new TH1F("h1PrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
	fOutputList->Add(fHistPrimaryVertexX);
	
	fHistPrimaryVertexY          = new TH1F("h1PrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
	fOutputList->Add(fHistPrimaryVertexY);
	
	fHistPrimaryVertexZ          = new TH1F("h1PrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
	fOutputList->Add(fHistPrimaryVertexZ);
	
	//////K0s///////////////// 2D histos: cut vs on fly status////
	
	fHistDcaPosToPrimVertexK0      = new TH2F("h2DcaPosToPrimVertexK0", "Positive V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
	fOutputList->Add(fHistDcaPosToPrimVertexK0);
	
	fHistDcaNegToPrimVertexK0      = new TH2F("h2DcaNegToPrimVertexK0", "Negative V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
	fOutputList->Add(fHistDcaNegToPrimVertexK0);
	
	
	fHistRadiusV0K0                = new TH2F("h2RadiusV0K0", "Radius;Radius(cm);Status",500,0,500,2,-0.5,1.5);
	fOutputList->Add(fHistRadiusV0K0);
	
	fHistDecayLengthV0K0           = new TH2F("h2DecayLengthV0K0", "V0s decay Length;decay length(cm);Status", 500, 0, 500,2,-0.5,1.5);
	fOutputList->Add(fHistDecayLengthV0K0);
	
	fHistDcaV0DaughtersK0          = new TH2F("h2DcaV0DaughtersK0", "DCA between daughters;dca(cm);Status", 300, 0, 3.0,2,-0.5,1.5);
	fOutputList->Add(fHistDcaV0DaughtersK0);
	
	fHistChi2K0                    = new TH2F("h2Chi2K0", "V0s chi2;chi2;Status", 1000, 0, 0.1,2,-0.5,1.5);
	fOutputList->Add(fHistChi2K0);
	
	fHistCosPointAngleK0           = new TH2F("h2CosPointAngleK0", "Cosine of V0's pointing angle", 200,0.99,1.01,2,-0.5,1.5);
	fOutputList->Add(fHistCosPointAngleK0);
	
	
	////////////K0s///////////////// 2D histos: cut vs mass////
	
	
	fHistDcaPosToPrimVertexK0vsMassK0 = new TH2F("h2DcaPosToPrimVertexK0vsMassK0", "Positive V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
	fOutputList->Add(fHistDcaPosToPrimVertexK0vsMassK0);
	
	fHistDcaNegToPrimVertexK0vsMassK0 = new TH2F("h2DcaNegToPrimVertexK0vsMassK0", "Negative V0 daughter;dca(cm);K0s inv. mass",500,0,10,200,0.4,0.6);
	fOutputList->Add(fHistDcaNegToPrimVertexK0vsMassK0);
	
	
	fHistRadiusV0K0vsMassK0           = new TH2F("h2RadiusV0K0vsMassK0", "Radius;Radius(cm);K0s inv. mass",110,0,110,200,0.4,0.6);
	fOutputList->Add(fHistRadiusV0K0vsMassK0);
	
	fHistDecayLengthV0K0vsMassK0      = new TH2F("h2DecayLengthV0K0vsMassK0", "V0s decay Length;decay length(cm);K0s inv. mass", 100, 0, 100,200,0.4,0.6);
	fOutputList->Add(fHistDecayLengthV0K0vsMassK0);
	
	fHistDcaV0DaughtersK0vsMassK0     = new TH2F("h2DcaV0DaughtersK0vsMassK0", "DCA between daughters;dca(cm);K0s inv. mass", 110, 0, 1.1,200,0.4,0.6);
	fOutputList->Add(fHistDcaV0DaughtersK0vsMassK0);
	
	
	fHistCosPointAngleK0vsMassK0      = new TH2F("h2CosPointAngleK0vsMassK0", "Cosine of V0's pointing angle", 200,0.997,1.007,200,0.4,0.6);
	fOutputList->Add(fHistCosPointAngleK0vsMassK0);
	//////////Lambda////////////// 2D histos: cut vs on fly status////
	
	fHistDcaPosToPrimVertexL      = new TH2F("h2DcaPosToPrimVertexL", "Positive V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
	fOutputList->Add(fHistDcaPosToPrimVertexL);
	
	fHistDcaNegToPrimVertexL      = new TH2F("h2DcaNegToPrimVertexL", "Negative V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
	fOutputList->Add(fHistDcaNegToPrimVertexL);
	
	fHistRadiusV0L                = new TH2F("h2RadiusV0L", "Radius;Radius(cm);Status",100,0,110,2,-0.5,1.5);
	fOutputList->Add(fHistRadiusV0L);
	
	fHistDecayLengthV0L           = new TH2F("h2DecayLengthV0L", "V0s decay Length;decay length(cm);Status", 500, 0, 500,2,-0.5,1.5);
	fOutputList->Add(fHistDecayLengthV0L);
	
	fHistDcaV0DaughtersL          = new TH2F("h2DcaV0DaughtersL", "DCA between daughters;dca(cm);Status", 300, 0, 3.0,2,-0.5,1.5);
	fOutputList->Add(fHistDcaV0DaughtersL);
	
	fHistChi2L                    = new TH2F("h2Chi2L", "V0s chi2;chi2;Status", 100, 0, 0.10,2,-0.5,1.5);
	fOutputList->Add(fHistChi2L);
	
	fHistCosPointAngleL           = new TH2F("h2CosPointAngleL", "Cosine of V0's pointing angle", 200,0.99,1.01,2,-0.5,1.5);
	fOutputList->Add(fHistCosPointAngleL);
	
	fHistcTauL                    = new TH1F("h1cTauL","cTaou of Lambdas",100,0,100);
	fOutputList->Add(fHistcTauL);
	//////////Lambda////////////// 2D histos: cut vs mass////
	fHistDcaPosToPrimVertexLvsMassL      = new TH2F("h2DcaPosToPrimVertexLvsMassL", "Positive V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaPosToPrimVertexLvsMassL);
	
	fHistDcaNegToPrimVertexLvsMassL      = new TH2F("h2DcaNegToPrimVertexLvsMassL", "Negative V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaNegToPrimVertexLvsMassL);
	
	
	fHistRadiusV0LvsMassL                = new TH2F("h2RadiusV0LvsMassL", "Radius;Radius(cm);Status",110,0,110,140, 1.06, 1.2);
	fOutputList->Add(fHistRadiusV0LvsMassL);
	
	fHistDecayLengthV0LvsMassL           = new TH2F("h2DecayLengthV0LvsMassL", "V0s decay Length;decay length(cm);Status", 120, 0, 120,140, 1.06, 1.2);
	fOutputList->Add(fHistDecayLengthV0LvsMassL);
	
	fHistDcaV0DaughtersLvsMassL          = new TH2F("h2DcaV0DaughtersLvsMassL", "DCA between daughters;dca(cm);Status", 110, 0, 1.1,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaV0DaughtersLvsMassL);
	
	fHistCosPointAngleLvsMassL           = new TH2F("h2CosPointAngleLvsMassL", "Cosine of V0's pointing angle", 200,0.997,1.007,140, 1.06, 1.2);
	fOutputList->Add(fHistCosPointAngleLvsMassL);
	
	fHistCosPointAngleLVsMassVsPtsigL           = new TH3F("h3McCosPointAngleLVsMassVsPtsigL", "Cosine of V0's pointing angle",3,0,12, 2,00.997,1.007,140, 1.06, 1.2);
	fOutputList->Add(fHistCosPointAngleLVsMassVsPtsigL);
	
	fHistCosPointAngleLVsMassVsPtbackL           = new TH3F("h3McCosPointAngleLVsMassVsPtbackL", "Cosine of V0's pointing angle",3,0,12, 20,0.997,1.007,140, 1.06, 1.2);
	fOutputList->Add(fHistCosPointAngleLVsMassVsPtbackL);
	
	
	//////////AntiLambda////////////// 2D histos: cut vs on fly status////
	
	fHistDcaPosToPrimVertexAntiL      = new TH2F("h2DcaPosToPrimVertexAntiL", "Positive V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
	fOutputList->Add(fHistDcaPosToPrimVertexAntiL);
	
	fHistDcaNegToPrimVertexAntiL      = new TH2F("h2DcaNegToPrimVertexAntiL", "Negative V0 daughter;dca(cm);Status",100,0,10,2,-0.5,1.5);
	fOutputList->Add(fHistDcaNegToPrimVertexAntiL);
	
	
	fHistRadiusV0AntiL                = new TH2F("h2RadiusV0AntiL", "Radius;Radius(cm);Status",100,0,110,2,-0.5,1.5);
	fOutputList->Add(fHistRadiusV0AntiL);
	
	fHistDecayLengthV0AntiL           = new TH2F("h2DecayLengthV0AntiL", "V0s decay Length;decay length(cm);Status", 500, 0, 500,2,-0.5,1.5);
	fOutputList->Add(fHistDecayLengthV0AntiL);
	
	fHistDcaV0DaughtersAntiL          = new TH2F("h2DcaV0DaughtersAntiL", "DCA between daughters;dca(cm);Status", 300, 0, 3.0,2,-0.5,1.5);
	fOutputList->Add(fHistDcaV0DaughtersAntiL);
	
	fHistChi2AntiL                    = new TH2F("h2Chi2AntiL", "V0s chi2;chi2;Status", 100, 0, 0.10,2,-0.5,1.5);
	fOutputList->Add(fHistChi2AntiL);
	
	fHistCosPointAngleAntiL           = new TH2F("h2CosPointAngleAntiL", "Cosine of V0's pointing angle", 200,0.99,1.01,2,-0.5,1.5);
	fOutputList->Add(fHistCosPointAngleAntiL);
	
	//////////AntiLambda////////////// 2D histos: cut vs mass////
	
	fHistDcaPosToPrimVertexAntiLvsMass      = new TH2F("h2DcaPosToPrimVertexAntiLvsMass", "Positive V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaPosToPrimVertexAntiLvsMass);
	
	fHistDcaNegToPrimVertexAntiLvsMass      = new TH2F("h2DcaNegToPrimVertexAntiLvsMass", "Negative V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaNegToPrimVertexAntiLvsMass);
	
	
	fHistRadiusV0AntiLvsMass                = new TH2F("h2RadiusV0AntiLvsMass", "Radius;Radius(cm);Status",110,0,110,140, 1.06, 1.2);
	fOutputList->Add(fHistRadiusV0AntiLvsMass);
	
	fHistDecayLengthV0AntiLvsMass           = new TH2F("h2DecayLengthV0AntiLvsMass", "V0s decay Length;decay length(cm);Status", 120, 0, 120,140, 1.06, 1.2);
	fOutputList->Add(fHistDecayLengthV0AntiLvsMass);
	
	fHistDcaV0DaughtersAntiLvsMass          = new TH2F("h2DcaV0DaughtersAntiLvsMass", "DCA between daughters;dca(cm);Status", 110, 0, 1.1,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaV0DaughtersAntiLvsMass);
	
	fHistCosPointAngleAntiLvsMass           = new TH2F("h2CosPointAngleAntiLvsMass", "Cosine of V0's pointing angle", 200,0.997,1.007,140, 1.06, 1.2);
	fOutputList->Add(fHistCosPointAngleAntiLvsMass);
	
	// Mass:
	fHistMassK0                   = new TH1F("h1MassK0", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 200, 0.4, 0.6);
	fOutputList->Add(fHistMassK0);
	
	fHistMassLambda               = new TH1F("h1MassLambda", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 150, 1.05, 1.2);
	fOutputList->Add(fHistMassLambda);
	
	fHistMassAntiLambda           = new TH1F("h1MassAntiLambda", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 150, 1.05, 1.2);
	fOutputList->Add(fHistMassAntiLambda);
	
	fHistMassVsRadiusK0           = new TH2F("h2MassVsRadiusK0", "K^{0} candidates;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",200,0,200, 200, 0.4, 0.6);
	fOutputList->Add(fHistMassVsRadiusK0);
	
	fHistMassVsRadiusLambda       = new TH2F("h2MassVsRadiusLambda", "#Lambda candidates;radius (cm);M(p#pi^{-}) (GeV/c^{2})",200,0,200, 140, 1.06, 1.2);
	fOutputList->Add(fHistMassVsRadiusLambda);
	
	fHistMassVsRadiusAntiLambda   = new TH2F("h2MassVsRadiusAntiLambda", "#bar{#Lambda} candidates;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",200,0,200, 140, 1.06, 1.2);
	fOutputList->Add(fHistMassVsRadiusAntiLambda);
	
	// Pt Vs Mass
	fHistPtVsMassK0               = new TH2F("h2PtVsMassK0","K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});p{t} (GeV/c)",400, 0.4, 0.6,240,0,12);
	fOutputList->Add(fHistPtVsMassK0);
	
	fHistPtVsMassLambda           = new TH2F("h2PtVsMassLambda","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p{t} (GeV/c)",280, 1.06, 1.2,240,0,12);
	fOutputList->Add(fHistPtVsMassLambda);
	
	fHistPtVsMassAntiLambda           = new TH2F("h2PtVsMassAntiLambda","#AntiLambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p{t} (GeV/c)",280, 1.06, 1.2,240,0,12);
	fOutputList->Add(fHistPtVsMassAntiLambda);
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	///Armenteros Podolansky
	fHistArmenterosPodolanski     = new TH2F("h2ArmenterosPodolanski","Armenteros-Podolanski phase space;#alpha;p{t} arm",100,-1.0,1.0,50,0,0.5);
	fOutputList->Add(fHistArmenterosPodolanski);
	
	//PID
	
	fHistNsigmaPosProtonLambda     = new TH1F("h1NsigmaPosProtonLambda", "Positive daughter of Lambda;NsigmaProton;Counts",25,0,5); 
	fOutputList->Add(fHistNsigmaPosProtonLambda);
	
	fHistNsigmaNegPionLambda       = new TH1F("h1NsigmaNegPionLambda", "Negative daughter of Lambda;NsigmaPion;Counts",25,0,5);
	fOutputList->Add(fHistNsigmaNegPionLambda);
	
	fHistNsigmaPosProtonAntiLambda     = new TH1F("h1NsigmaPosProtonAntiLambda", "Positive daughter of AntiLambda;NsigmaProton;Counts",25,0,5); 
	fOutputList->Add(fHistNsigmaPosProtonAntiLambda);
	
	fHistNsigmaNegPionAntiLambda       = new TH1F("h1NsigmaNegPionAntiLambda", "Negative daughter of AntiLambda;NsigmaPion;Counts",25,0,5);
	fOutputList->Add(fHistNsigmaNegPionAntiLambda);
	
	fHistNsigmaPosPionK0           = new TH1F("h1NsigmaPosPionK0", "Positive daughter of K0s;NsigmaPion;Counts",25,0,5);
	fOutputList->Add(fHistNsigmaPosPionK0);
	
	fHistNsigmaNegPionK0           = new TH1F("h1NsigmaNegPionK0", "Negative daughter of K0s;NsigmaPion;Counts",25,0,5);
	fOutputList->Add(fHistNsigmaNegPionK0);
	
	
//--------------------------------------------MC Associated histograms -----------------------------------------------------//
	
	//Pt distribution
	fHistAsMcPtK0                = new TH1F("h1AsMcPtK0", "K^{0} associated;p{t} (GeV/c);Counts", 240,0,12);
	fOutputList->Add(fHistAsMcPtK0);
	
	fHistAsMcPtLambda            = new TH1F("h1AsMcPtLambda", "#Lambda^{0} associated;p{t} (GeV/c);Counts", 240,0,12);
	fOutputList->Add(fHistAsMcPtLambda);
	
	fHistAsMcPtAntiLambda            = new TH1F("h1AsMcPtAntiLambda", "#AntiLambda^{0} associated;p{t} (GeV/c);Counts", 240,0,12);
	fOutputList->Add(fHistAsMcPtAntiLambda);
	
	// Radius distribution
	fHistAsMcProdRadiusK0               = new TH1F("h1AsMcProdRadiusK0", "K^{0} associated;r (cm);Counts", 500, 0, 100);
	fOutputList->Add(fHistAsMcProdRadiusK0);
	
	fHistAsMcProdRadiusLambda           = new TH1F("h1AsMcProdRadiusLambda", "#Lambda^{0} associated;r (cm);Counts", 500, 0, 100);
	fOutputList->Add(fHistAsMcProdRadiusLambda);
	
	fHistAsMcProdRadiusAntiLambda       = new TH1F("h1AsMcProdRadiusAntiLambda", "#bar{#Lambda}^{0} associated;r (cm);Counts", 500, 0, 100);
	fOutputList->Add(fHistAsMcProdRadiusAntiLambda);
	
	fHistAsMcProdRadiusXvsYK0s          = new TH2F("h2AsMcProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
	fOutputList->Add(fHistAsMcProdRadiusXvsYK0s);
	
	fHistAsMcProdRadiusXvsYLambda       = new TH2F("h2AsMcProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
	fOutputList->Add(fHistAsMcProdRadiusXvsYLambda);
	
	fHistAsMcProdRadiusXvsYAntiLambda   = new TH2F("h2AsMcProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
	fOutputList->Add(fHistAsMcProdRadiusXvsYAntiLambda);
	
	// Mass
	fHistPidMcMassK0             = new TH1F("h1PidMcMassK0", "K^{0} MC PId checked;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
	fOutputList->Add(fHistPidMcMassK0);
	
	fHistPidMcMassLambda         = new TH1F("h1PidMcMassLambda", "#Lambda^{0} MC PId checked;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
	fOutputList->Add(fHistPidMcMassLambda);
	
	fHistPidMcMassAntiLambda     = new TH1F("h1PidMcMassAntiLambda", "#bar{#Lambda}^{0} MC PId checked;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
	fOutputList->Add(fHistPidMcMassAntiLambda);
	
	fHistAsMcMassK0              = new TH1F("h1AsMcMassK0", "K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
	fOutputList->Add(fHistAsMcMassK0);
	
	fHistAsMcMassLambda          = new TH1F("h1AsMcMassLambda", "#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
	fOutputList->Add(fHistAsMcMassLambda);
	
	fHistAsMcMassAntiLambda      = new TH1F("h1AsMcMassAntiLambda", "#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
	fOutputList->Add(fHistAsMcMassAntiLambda);
	
	//Pt versus Mass
	fHistAsMcPtVsMassK0               = new TH2F("h2AsMcPtVsMassK0","K^{0} associated;M(#pi^{+}#pi^{-}) (GeV/c^{2});p{t} (GeV/c)",200, 0.4, 0.6,240,0,12);
	fOutputList->Add(fHistAsMcPtVsMassK0);
	
	fHistAsMcPtVsMassLambda           = new TH2F("h2AsMcPtVsMassLambda","#Lambda^{0} associated;M(p#pi^{-}) (GeV/c^{2});p{t} (GeV/c)",140, 1.06, 1.2,240,0,12);
	fOutputList->Add(fHistAsMcPtVsMassLambda);
	
	fHistAsMcPtVsMassAntiLambda       = new TH2F("h2AsMcPtVsMassAntiLambda","#bar{#Lambda}^{0} associated;M(#bar{p}#pi^{+}) (GeV/c^{2});p{t} (GeV/c)",140, 1.06, 1.2,240,0,12);
	fOutputList->Add(fHistAsMcPtVsMassAntiLambda);
	
	
	// invariant mass vs radius
	fHistAsMcMassVsRadiusK0             = new TH2F("h2AsMcMassVsRadiusK0", "K^{0} associated;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",200,0,200, 500, 0.47, 0.52);
	fOutputList->Add(fHistAsMcMassVsRadiusK0);
	
	fHistAsMcMassVsRadiusLambda         = new TH2F("h2AsMcMassVsRadiusLambda", "#Lambda associated;radius (cm);M(p#pi^{-}) (GeV/c^{2})",200,0,200, 1.10, 1.13);
	fOutputList->Add(fHistAsMcMassVsRadiusLambda);
	
	fHistAsMcMassVsRadiusAntiLambda     = new TH2F("h2AsMcMassVsRadiusAntiLambda", "#bar{#Lambda} associated;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",200,0,200 , 1.10, 1.13);
	fOutputList->Add(fHistAsMcMassVsRadiusAntiLambda);
	
	// Position Resolution
	fHistAsMcResxK0                     = new TH1F("h1AsMcResxK0", "K^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcResxK0);
	fHistAsMcResyK0                     = new TH1F("h1AsMcResyK0", "K^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcResyK0);
	fHistAsMcReszK0                     = new TH1F("h1AsMcReszK0", "K^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcReszK0);
	fHistAsMcResrVsRadiusK0             = new TH2F("h2AsMcResrVsRadiusK0", "K^{0} associated;Radius (cm);#Delta r (cm)",200,0.0,50., 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcResrVsRadiusK0);
	fHistAsMcReszVsRadiusK0             = new TH2F("h2AsMcReszVsRadiusK0", "K^{0} associated;Radius (cm);#Delta z (cm)",200,0.0,50.0, 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcReszVsRadiusK0);
	
	fHistAsMcResxLambda                 = new TH1F("h1AsMcResxLambda", "#Lambda^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcResxLambda);
	fHistAsMcResyLambda                 = new TH1F("h1AsMcResyLambda", "#Lambda^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcResyLambda);
	fHistAsMcReszLambda                 = new TH1F("h1AsMcReszLambda", "#Lambda^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcReszLambda);
	fHistAsMcResrVsRadiusLambda         = new TH2F("h2AsMcResrVsRadiusLambda", "#Lambda^{0} associated;Radius (cm);#Delta r (cm)",200,0.0,50.0, 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcResrVsRadiusLambda);
	fHistAsMcReszVsRadiusLambda         = new TH2F("h2AsMcReszVsRadiusLambda", "#Lambda^{0} associated;Radius (cm);#Delta z (cm)",200,0.0,50.0, 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcReszVsRadiusLambda);
	
	fHistAsMcResxAntiLambda             = new TH1F("h1AsMcResxAntiLambda", "#bar{#Lambda}^{0} associated;#Delta x (cm);Counts", 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcResxAntiLambda);
	fHistAsMcResyAntiLambda             = new TH1F("h1AsMcResyAntiLambda", "#bar{#Lambda}^{0} associated;#Delta y (cm);Counts", 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcResyAntiLambda);
	fHistAsMcReszAntiLambda             = new TH1F("h1AsMcReszAntiLambda", "#bar{#Lambda}^{0} associated;#Delta z (cm);Counts", 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcReszAntiLambda);
	fHistAsMcResrVsRadiusAntiLambda     = new TH2F("h2AsMcResrVsRadiusAntiLambda", "#bar{#Lambda}^{0} associated;Radius (cm);#Delta r (cm)",200,0.0,50.0, 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcResrVsRadiusAntiLambda);
	fHistAsMcReszVsRadiusAntiLambda     = new TH2F("h2AsMcReszVsRadiusAntiLambda", "#bar{#Lambda}^{0} associated;Radius (cm);#Delta z (cm)",200,0.0,50.0, 50, -0.25, 0.25);
	fOutputList->Add(fHistAsMcReszVsRadiusAntiLambda);
	
	// Pt Resolution
	fHistAsMcResPtK0                   = new TH1F("h1AsMcResPtK0","Pt Resolution K^{0};#Delta Pt;Counts",200,-1,1);
	fOutputList->Add(fHistAsMcResPtK0);
	
	fHistAsMcResPtLambda               = new TH1F("h1AsMcResPtLambda","Pt Resolution #Lambda^{0};#Delta Pt;Counts",200,-1,1);
	fOutputList->Add(fHistAsMcResPtLambda);
	
	fHistAsMcResPtAntiLambda           = new TH1F("h1AsMcResPtAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Counts",200,-1,1);
	fOutputList->Add(fHistAsMcResPtAntiLambda);
	
	
	fHistAsMcResPtVsRapK0              = new TH2F("h2AsMcResPtVsRapK0","Pt Resolution K^{0};#Delta Pt;Rap",200,-1,1,20,-1,1);
	fOutputList->Add(fHistAsMcResPtVsRapK0);
	
	fHistAsMcResPtVsRapLambda          = new TH2F("h2AsMcResPtVsRapLambda","Pt Resolution #Lambda^{0};#Delta Pt;Rap",200,-1,1,20,-1,1);
	fOutputList->Add(fHistAsMcResPtVsRapLambda);
	
	fHistAsMcResPtVsRapAntiLambda      = new TH2F("h2AsMcResPtVsRapAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Rap",200,-1,1,20,-1,1);
	fOutputList->Add(fHistAsMcResPtVsRapAntiLambda);
	
	fHistAsMcResPtVsPtK0               = new TH2F("h2AsMcResPtVsPtK0","Pt Resolution K^{0};#Delta Pt;Pt",600,-0.15,0.15,240,0,12);
	fOutputList->Add(fHistAsMcResPtVsPtK0);
    
	fHistAsMcResPtVsPtLambda           = new TH2F("h2AsMcResPtVsPtLambda","Pt Resolution #Lambda^{0};#Delta Pt;Pt",600,-0.15,0.15,240,0,12);
	fOutputList->Add(fHistAsMcResPtVsPtLambda);
	
	fHistAsMcResPtVsPtAntiLambda       = new TH2F("h2AsMcResPtVsPtAntiLambda","Pt Resolution #bar{#Lambda}^{0};#Delta Pt;Pt",300,-0.15,0.15,240,0,12);
	fOutputList->Add(fHistAsMcResPtVsPtAntiLambda);
	
	// Pt distribution Lambda from Sigma
	fHistAsMcPtLambdaFromSigma          = new TH1F("h1AsMcPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p{t} (GeV/c);Count",240,0,12);
	fOutputList->Add(fHistAsMcPtLambdaFromSigma);
	
	fHistAsMcPtAntiLambdaFromSigma      = new TH1F("h1AsMcPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p{t} (GeV/c);Count",240,0,12);
	fOutputList->Add(fHistAsMcPtAntiLambdaFromSigma);
	
	// Associated secondary particles:
	// Pt and rapidity distribution
	fHistAsMcSecondaryPtVsRapK0s          = new TH2F("h2AsMcSecondaryPtVsRapK0s", "K^{0} associated secondary;p{t} (GeV/c);rapidity",240,0,12,30,-1.5,1.5);
	fOutputList->Add(fHistAsMcSecondaryPtVsRapK0s);
	
	fHistAsMcSecondaryPtVsRapLambda       = new TH2F("h2AsMcSecondaryPtVsRapLambda", "#Lambda^{0} associated secondary;p{t} (GeV/c);rapidity",240,0,12,30,-1.5,1.5);
	fOutputList->Add(fHistAsMcSecondaryPtVsRapLambda);
	
	fHistAsMcSecondaryPtVsRapAntiLambda   = new TH2F("h2AsMcSecondaryPtVsRapAntiLambda", "#bar{#Lambda}^{0} associated secondary;p{t} (GeV/c);rapidity",240,0,12,30,-1.5,1.5);
	fOutputList->Add(fHistAsMcSecondaryPtVsRapAntiLambda);
	
	// Production radius
	fHistAsMcSecondaryProdRadiusK0s              = new TH1F("h1AsMcSecondaryProdRadiusK0s", "K^{0} Production Radius;r (cm);Count", 170, -2, 15);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusK0s);
	
	fHistAsMcSecondaryProdRadiusLambda           = new TH1F("h1AsMcSecondaryProdRadiusLambda", "#Lambda^{0} Production Radius;r (cm);Count", 170, -2, 15);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusLambda);
	
	fHistAsMcSecondaryProdRadiusAntiLambda       = new TH1F("h1AsMcSecondaryProdRadiusAntiLambda", "#bar{#Lambda}^{0} Production Radius;r (cm);Count", 170, -2, 15);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusAntiLambda);  
	
	fHistAsMcSecondaryProdRadiusXvsYK0s          = new TH2F("h2AsMcSecondaryProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusXvsYK0s);
	
	fHistAsMcSecondaryProdRadiusXvsYLambda       = new TH2F("h2AsMcSecondaryProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusXvsYLambda);
	
	fHistAsMcSecondaryProdRadiusXvsYAntiLambda   = new TH2F("h2AsMcSecondaryProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusXvsYAntiLambda);
	
	// Pt distribution Lambda from Sigma
	fHistAsMcSecondaryPtLambdaFromSigma          = new TH1F("h1AsMcSecondaryPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p{t} (GeV/c);Count",240,0,12);
	fOutputList->Add(fHistAsMcSecondaryPtLambdaFromSigma);
	
	fHistAsMcSecondaryPtAntiLambdaFromSigma      = new TH1F("h1AsMcSecondaryPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p{t} (GeV/c);Count",240,0,12);
	fOutputList->Add(fHistAsMcSecondaryPtAntiLambdaFromSigma);
	
	
	//----------------------------------------------Correlation histograms -----------------------------------------------------//
	
	fHistSibK0 = new TH2F("hfHistSibK0","",CorrBinsX,-2*fTrackEtaCut,2*fTrackEtaCut,CorrBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
	fHistSibK0->SetXTitle("#Delta #Phi");
	fHistSibK0->SetStats(0);
	fHistSibK0->Sumw2();
	fOutputList->Add(fHistSibK0);
	
	fHistMixK0 = new TH2F("hfHistMixK0","",CorrBinsX,-2*fTrackEtaCut,2*fTrackEtaCut,CorrBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
	fHistMixK0->SetXTitle("#Delta #Phi");
	fHistMixK0->SetStats(0);
	fHistMixK0->Sumw2();
	fOutputList->Add(fHistMixK0);
	
	fHistSibLambda = new TH2F("hfHistSibLambda","",CorrBinsX,-2*fTrackEtaCut,2*fTrackEtaCut,CorrBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
	fHistSibLambda->SetXTitle("#Delta #Phi");
	fHistSibLambda->SetStats(0);
	fHistSibLambda->Sumw2();
	fOutputList->Add(fHistSibLambda);
	
	fHistMixLambda= new TH2F("hfHistMixLambda","",CorrBinsX,-2*fTrackEtaCut,2*fTrackEtaCut,CorrBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
	fHistMixLambda->SetXTitle("#Delta #Phi");
	fHistMixLambda->SetStats(0);
	fHistMixLambda->Sumw2();
	fOutputList->Add(fHistMixLambda);
	
	fHistSibK0MC = new TH2F("hfHistSibK0MC","",CorrBinsX,-2*fTrackEtaCut,2*fTrackEtaCut,CorrBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
	fHistSibK0MC->SetXTitle("#Delta #Phi");
	fHistSibK0MC->SetStats(0);
	fHistSibK0MC->Sumw2();
	fOutputList->Add(fHistSibK0MC);
	
	fHistMixK0MC = new TH2F("hfHistMixK0MC","",CorrBinsX,-2*fTrackEtaCut,2*fTrackEtaCut,CorrBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
	fHistMixK0MC->SetXTitle("#Delta #Phi");
	fHistMixK0MC->SetStats(0);
	fHistMixK0MC->Sumw2();
	fOutputList->Add(fHistMixK0MC);
	
	fHistSibLambdaMC = new TH2F("hfHistSibLambdaMC","",CorrBinsX,-2*fTrackEtaCut,2*fTrackEtaCut,CorrBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
	fHistSibLambdaMC->SetXTitle("#Delta #Phi");
	fHistSibLambdaMC->SetStats(0);
	fHistSibLambdaMC->Sumw2();
	fOutputList->Add(fHistSibLambdaMC);
	
	fHistMixLambdaMC= new TH2F("hfHistMixLambdaMC","",CorrBinsX,-2*fTrackEtaCut,2*fTrackEtaCut,CorrBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
	fHistMixLambdaMC->SetXTitle("#Delta #Phi");
	fHistMixLambdaMC->SetStats(0);
	fHistMixLambdaMC->Sumw2();
	fOutputList->Add(fHistMixLambdaMC);
	
	//----------------------------------------------Event Pool-----------------------------------------------------//
	
	fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins, fCentBins, fNzVtxBins, fZvtxBins);
	if(!fPoolMgr) return;
	
	PostData(1, fOutputList);

}
//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::UserExec(Option_t *)
{
	
    AliAnalysisManager   *mgr      = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inEvMain = (AliInputEventHandler*)(mgr->GetInputEventHandler());
	if (!inEvMain) return;
	
	// Pointers to PID Response objects.	
	fPIDResponse = inEvMain->GetPIDResponse(); 
	//cout << "PID Response object: " << fPIDResponse << endl;

    fAODEvent = dynamic_cast<AliAODEvent*>(inEvMain->GetEvent());
	if(!fAODEvent) return;
	
	
	// physics selection
	UInt_t maskIsSelected = inEvMain->IsEventSelected();
    Bool_t isSelected = ((maskIsSelected & AliVEvent::kMB) || (maskIsSelected & AliVEvent::kCentral) || (maskIsSelected & AliVEvent::kSemiCentral));
    if (!isSelected) return;
	
	//-----------------------------MC Accsess------------------------------------------------
	TClonesArray *stack = 0x0;
	Double_t mcXv=0., mcYv=0., mcZv=0.;
	Int_t ntrk =0, ntrk0=0;
	
	TObjArray *selectedTracksLeadingK0MC =0x0;
    TObjArray *selectedTracksLeadingLambdaMC = 0x0;
	
	TObjArray * selectedK0MC =new TObjArray;
	selectedK0MC->SetOwner(kTRUE);
	
	TObjArray * selectedLambdaMC =new TObjArray;
	selectedLambdaMC->SetOwner(kTRUE);
	
	if (fAnalysisMC) {
		TList *lst = fAODEvent->GetList();
		stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
		if (!stack) {Printf("ERROR: stack not available");return;}
		
		AliAODMCHeader *mcHdr=(AliAODMCHeader*)lst->FindObject(AliAODMCHeader::StdBranchName());
		mcXv=mcHdr->GetVtxX(); mcYv=mcHdr->GetVtxY(); mcZv=mcHdr->GetVtxZ();
		ntrk=stack->GetEntriesFast(), ntrk0=ntrk;
		if(TMath::Abs(mcZv)>fpvzcut)return;
		
		selectedTracksLeadingK0MC = FindLeadingObjectsK0MC(stack);
		if(!selectedTracksLeadingK0MC) return;
		selectedTracksLeadingK0MC->SetOwner(kTRUE);
		
		selectedTracksLeadingLambdaMC = FindLeadingObjectsLambdaMC(stack);
		if(!selectedTracksLeadingLambdaMC) return;
		selectedTracksLeadingLambdaMC->SetOwner(kTRUE);
		
	}
	
	// If PID is used:
	Double_t lLimitPPID    = 0.7;
	Float_t cutNSigmaLowP  = 1E3;
	Float_t cutNSigmaHighP = 1E3;
	if (fUsePID=="withPID") {
		cutNSigmaLowP  = 3.0;
		cutNSigmaHighP = 3.0;
	}
	
	Double_t lmcPrimVtxR      = 0;
	
	Int_t lPdgcodeCurrentPart = 0;
	Double_t lRapCurrentPart  = 0;
	Double_t lPtCurrentPart   = 0;
	Double_t lPhiCurrentPart  = 0;
	Double_t lEtaCurrentPart  = 0;
	Int_t lComeFromSigma      = 0;
	
	// PID flags:
	Int_t LambdaPID = 0;
	Int_t AntiLambdaPID = 0;
	
	
	// Production Radius
	Double_t mcPosX     = 0.0,  mcPosY      = 0.0,  mcPosZ      = 0.0;
	Double_t mcPosR     = 0.0;
	
	// Decay Radius
	Double_t mcDecayPosX = 0, mcDecayPosY = 0, mcDecayPosR = 0;
	
	// current mc particle 's mother
	Int_t iCurrentMother  = 0, lPdgCurrentMother    = 0;
	
	// current mc particles 's daughter:
	Int_t lPdgCurrentDaughter0 = 0, lPdgCurrentDaughter1 = 0; 
	
	// variables for multiple reconstruction studies:
	Int_t id0           = 0, id1          = 0;
	Int_t lNtimesReconstructedK0s   = 0, lNtimesReconstructedLambda   = 0, lNtimesReconstructedAntiLambda   = 0;

	
	// Start loop over MC particles
	if (fAnalysisMC) {
		
		// Primary vertex
		fHistMCPrimaryVertexX->Fill(mcXv);
		fHistMCPrimaryVertexY->Fill(mcYv);
		fHistMCPrimaryVertexZ->Fill(mcZv);
		
		lmcPrimVtxR = TMath::Sqrt(mcXv*mcXv+mcYv*mcYv);
     
		
		for (Int_t iMc = 0; iMc < (ntrk); iMc++) {  
			AliAODMCParticle *p0=(AliAODMCParticle*)stack->UncheckedAt(iMc);
			if (!p0) continue;
			
			lPdgcodeCurrentPart = p0->GetPdgCode();
			
			// Keep only K0s, Lambda and AntiLambda, Xi and Phi:
			if ( (lPdgcodeCurrentPart != 310 ) && (lPdgcodeCurrentPart != 3122 ) && (lPdgcodeCurrentPart != -3122 ) && (lPdgcodeCurrentPart != 3312 ) && (lPdgcodeCurrentPart != -3312) && (lPdgcodeCurrentPart != -333) ) continue;
			
			lRapCurrentPart   = p0->Y();
			lPtCurrentPart    = p0->Pt();
			lPhiCurrentPart   = p0->Phi();
			lEtaCurrentPart   = p0->Eta();
			iCurrentMother    = p0->GetMother();
			
			AliAODMCParticle *Mother = (AliAODMCParticle*)stack->UncheckedAt(iCurrentMother);
			if (iCurrentMother == -1){lPdgCurrentMother=0; } else {lPdgCurrentMother = Mother->GetPdgCode();} 
			
			mcPosX = p0->Xv();
			mcPosY = p0->Yv();
			mcPosZ = p0->Zv();
			mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
			
			id0  = p0->GetDaughter(0);
			id1  = p0->GetDaughter(1);
			
			// Decay Radius and Production Radius
			if ( id0 <= ntrk && id0 > 0 && id1 <= ntrk && id1 > 0) {
				AliAODMCParticle *pDaughter0 = (AliAODMCParticle*)stack->UncheckedAt(id0);
				AliAODMCParticle *pDaughter1 = (AliAODMCParticle*)stack->UncheckedAt(id1);
				lPdgCurrentDaughter0 = pDaughter0->GetPdgCode();
				lPdgCurrentDaughter1 = pDaughter1->GetPdgCode();
				
				mcDecayPosX = pDaughter0->Xv();
				mcDecayPosY = pDaughter0->Yv();
				mcDecayPosR = TMath::Sqrt(mcDecayPosX*mcDecayPosX+mcDecayPosY*mcDecayPosY);
			}
			else  {mcDecayPosR = -1.0;}
			
			if (lPdgcodeCurrentPart==310)   {
				fHistMCtracksProdRadiusK0s->Fill(mcPosX,mcPosY);
				fHistMCtracksDecayRadiusK0s->Fill(mcDecayPosR);
				if (TMath::Abs(lRapCurrentPart) < fRapidityCut) fHistMCPtAllK0s->Fill(lPtCurrentPart);
			}
			else if (lPdgcodeCurrentPart==3122)  {
				fHistMCtracksProdRadiusLambda->Fill(mcPosX,mcPosY);
				fHistMCtracksDecayRadiusLambda->Fill(mcDecayPosR);
				if (TMath::Abs(lRapCurrentPart) < fRapidityCut) fHistMCPtAllLambda->Fill(lPtCurrentPart);
			}
			else if (lPdgcodeCurrentPart==-3122) {
				fHistMCtracksProdRadiusAntiLambda->Fill(mcPosX,mcPosY);
				fHistMCtracksDecayRadiusAntiLambda->Fill(mcDecayPosR);
				if (TMath::Abs(lRapCurrentPart) < fRapidityCut) fHistMCPtAllAntiLambda->Fill(lPtCurrentPart);
			}
			
			if ( ( ( TMath::Abs(lPdgCurrentMother) == 3212)  ||
				  ( TMath::Abs(lPdgCurrentMother) == 3224)  ||
				  ( TMath::Abs(lPdgCurrentMother) == 3214)  ||
				  ( TMath::Abs(lPdgCurrentMother) == 3114) )
			   && ( Mother->GetMother() == -1)
				) lComeFromSigma = 1;
			else lComeFromSigma = 0;                                                                                                                                                              
			
			Double_t dx = 0;
			Double_t dy = 0;
			Double_t dz = 0;
			Double_t ProdDistance = 0;
			
			
			dx = ( ( mcXv) - (mcPosX) );
			dy = ( ( mcYv) - (mcPosY) );
			dz = ( ( mcZv) - (mcPosZ) );
			
			ProdDistance = TMath::Sqrt(dx*dx + dy*dy + dz*dz);
			if (ProdDistance > 0.001) continue; // secondary V0      
			
			//********************************************
			
			lNtimesReconstructedK0s   = 0; lNtimesReconstructedLambda   = 0; lNtimesReconstructedAntiLambda   = 0;
			
			// Rapidity Cut
			if (TMath::Abs(lRapCurrentPart) > fRapidityCut) continue;
			
			if (lPdgcodeCurrentPart==310) {
				fHistMCRapK0s->Fill(lRapCurrentPart);
				fHistMCProdRadiusK0s->Fill(mcPosR);
				fHistMCPtK0s->Fill(lPtCurrentPart);
				fHistNTimesRecK0s->Fill(lNtimesReconstructedK0s);
				fHistNTimesRecK0sVsPt->Fill(lPtCurrentPart,lNtimesReconstructedK0s);
				fHistPrimRawPtVsYK0s->Fill(lPtCurrentPart,lRapCurrentPart);
				selectedK0MC->Add(new AliLeadingBasicParticle(lEtaCurrentPart,lPhiCurrentPart,lPtCurrentPart));
			}
			else 
			if (lPdgcodeCurrentPart==3122) {
				    fHistMCRapLambda->Fill(lRapCurrentPart);
					fHistMCProdRadiusLambda->Fill(mcPosR);
					fHistMCPtLambda->Fill(lPtCurrentPart);	  
					fHistNTimesRecLambda->Fill(lNtimesReconstructedLambda);
					fHistNTimesRecLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedLambda);
				    fHistPrimRawPtVsYLambda->Fill(lPtCurrentPart,lRapCurrentPart);
					if (lComeFromSigma) fHistMCPtLambdaFromSigma->Fill(lPtCurrentPart);
				    selectedLambdaMC->Add(new AliLeadingBasicParticle(lEtaCurrentPart,lPhiCurrentPart,lPtCurrentPart));
					
				}
			else 
			if (lPdgcodeCurrentPart==-3122) {
				        fHistMCRapAntiLambda->Fill(lRapCurrentPart);
						fHistMCProdRadiusAntiLambda->Fill(mcPosR);
						fHistMCPtAntiLambda->Fill(lPtCurrentPart);	  
						fHistNTimesRecAntiLambda->Fill(lNtimesReconstructedAntiLambda);
						fHistNTimesRecAntiLambdaVsPt->Fill(lPtCurrentPart,lNtimesReconstructedAntiLambda);
						fHistPrimRawPtVsYAntiLambda->Fill(lPtCurrentPart,lRapCurrentPart);
						if (lComeFromSigma) fHistMCPtAntiLambdaFromSigma->Fill(lPtCurrentPart);
					}
			
		} // end loop AOD MC
		
	} // End Loop over MC condition
	
	//-----------------------------------------------------------------------------------------
	
	// Vertex cut
	Double_t  lPrimaryVtxPosition[3];
	AliAODVertex *myPrimVertex = fAODEvent->GetPrimaryVertex();
	if (!myPrimVertex) return;
	myPrimVertex->GetXYZ(lPrimaryVtxPosition);
	
	Double_t lPVx = lPrimaryVtxPosition[0];
	Double_t lPVy = lPrimaryVtxPosition[1];
	Double_t lPVz = lPrimaryVtxPosition[2];
	if ((TMath::Abs(lPVz)) >= fpvzcut) return ;
	
	if (TMath::Abs(lPVx)<10e-5 && TMath::Abs(lPVy)<10e-5 && TMath::Abs(lPVz)<10e-5) return;
	
	
	fHistPrimaryVertexX->Fill(lPVx);
    fHistPrimaryVertexY->Fill(lPVy);
    fHistPrimaryVertexZ->Fill(lPVz);
	
	// Centrality definition
	AliCentrality *centralityObj = 0;
	Int_t multiplicity = -1;
	Double_t MultipOrCent = -1;
	
	// initialize the pool for event mixing
	if(fcollidingSys=="PP"){ 
		multiplicity = fAODEvent->GetNTracks();
		MultipOrCent = multiplicity; // convert from Int_t to Double_t
	}
	if(fcollidingSys=="PbPb"){ 
		centralityObj = fAODEvent->GetHeader()->GetCentralityP();
		MultipOrCent  = centralityObj->GetCentralityPercentileUnchecked("V0M");
	}
	
	Double_t * CentBins = fCentBins;
	Double_t poolmin=CentBins[0];
	Double_t poolmax=CentBins[fNCentBins];
	
	TObjArray *selectedTracksLeadingK0 = FindLeadingObjectsK0(fAODEvent);
	if(!selectedTracksLeadingK0) return;
	selectedTracksLeadingK0->SetOwner(kTRUE);
	
	TObjArray *selectedTracksLeadingLambda = FindLeadingObjectsLambda(fAODEvent);
	if(!selectedTracksLeadingLambda) return;
	selectedTracksLeadingLambda->SetOwner(kTRUE);
	
	// -------------------------------------V0 loop for reconstructed event------------------------
	
	Double_t cutcTauL   = 3*7.89;
	Double_t cutcTauK0  = 3*2.68;

	Double_t lPLambda = 0;
	Double_t lPAntiLambda = 0;
	Double_t lPK0s = 0;
	
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
	Double_t lcTauLambda         = 0;   
	Double_t lcTauAntiLambda     = 0;   
	Double_t lcTauK0s            = 0;   
	Int_t    lOnFlyStatus        = 0;
	
	Double_t lInvMassK0   = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
	Double_t lPtK0s       = 0, lPtLambda      = 0, lPtAntiLambda      = 0;
	Double_t lPhiK0s      = 0, lPhiLambda     = 0, lPhiAntiLambda     = 0;
	Double_t lEtaK0s      = 0, lEtaLambda     = 0, lEtaAntiLambda     = 0;
	Double_t lRapK0s      = 0, lRapLambda     = 0, lRapAntiLambda     = 0;
	Double_t lPzK0s       = 0, lPzLambda      = 0, lPzAntiLambda      = 0;
	Double_t lAlphaV0     = 0, lPtArmV0       = 0;
	
	
	
	Double_t lV0Eta = 999;
	
	//Associated V0s:
	
	UInt_t   lLabelTrackPos       = 0, lLabelTrackNeg         = 0;
	Int_t    lCheckPIdK0Short     = 0, lCheckMcK0Short        = 0;
	Int_t    lCheckPIdLambda      = 0, lCheckMcLambda         = 0;
	Int_t    lCheckPIdAntiLambda  = 0, lCheckMcAntiLambda     = 0;
	Int_t    lCheckSecondaryK0s   = 0, lCheckSecondaryLambda  = 0, lCheckSecondaryAntiLambda  = 0;
	Int_t    lCheckGamma          = 0;
	Double_t mcPosMotherX         = 0, mcPosMotherY           = 0, mcPosMotherZ  = 0;
	Double_t mcPosMotherR         = 0;
	Double_t mcMotherPt           = 0;
	
	Int_t lIndexPosMother         = 0;
	Int_t lIndexNegMother         = 0;
	Int_t lIndexMotherOfMother    = 0;
	Int_t lPDGCodePosDaughter     = 0;
	Int_t lPDGCodeNegDaughter     = 0;
	Int_t lPdgcodeMother          = 0;
	Int_t lPdgcodeMotherOfMother  = 0;
	
	// Reconstructed position
	Double_t rcPosXK0s        = 0,  rcPosYK0s        = 0, rcPosZK0s        = 0;
	Double_t rcPosRK0s        = 0;
	Double_t rcPosXLambda     = 0,  rcPosYLambda     = 0, rcPosZLambda     = 0;
	Double_t rcPosRLambda     = 0;
	Double_t rcPosXAntiLambda = 0,  rcPosYAntiLambda = 0, rcPosZAntiLambda = 0;
	Double_t rcPosRAntiLambda = 0;
	
	// Pt resolution
	Double_t deltaPtK0s  = 0, deltaPtLambda  = 0, deltaPtAntiLambda  = 0;
	
	//  V0 momentum      
	//  Double_t V0mom[3] = {999,999,999};
	Double_t lPosMom = 0;
	Double_t lNegMom = 0;
	
	// Daughters' momentum:
	Double_t  lMomPos[3] = {999,999,999};
	Double_t  lMomNeg[3] = {999,999,999};
	Double_t  lPtPos = 999, lPtNeg = 999;
	Double_t  lPPos = 999, lPNeg = 999;
	
	// Inner Wall parameters:
	Double_t  lMomInnerWallPos =999, lMomInnerWallNeg = 999;
	Double_t 	ldEdxPos =0.0,       ldEdxNeg=0.0;
	
	// PID
	Float_t nSigmaPosPion   = 0;
	Float_t nSigmaNegPion   = 0;
	
	Float_t nSigmaPosProton = 0;
	Float_t nSigmaNegProton = 0;
	
	
	Int_t lCheckPIDK0sPosDaughter        = 0, lCheckPIDK0sNegDaughter        = 0;
	Int_t lCheckPIDLambdaPosDaughter     = 0, lCheckPIDLambdaNegDaughter     = 0;
	Int_t lCheckPIDAntiLambdaPosDaughter = 0, lCheckPIDAntiLambdaNegDaughter = 0;

	//---------------------------------------------------------------------------------------------
	TObjArray * selectedK0 = new TObjArray;
	selectedK0->SetOwner(kTRUE);
	
	TObjArray * selectedLambda = new TObjArray;
	selectedLambda->SetOwner(kTRUE);
	
	Int_t nV0s = fAODEvent->GetNumberOfV0s();
	
	for (Int_t i = 0; i < nV0s; i++)
	{ // start of V0 slection loop
		AliAODv0* aodV0 = dynamic_cast<AliAODv0 *>(fAODEvent->GetV0(i));
		if (!aodV0) {AliError(Form("ERROR: Could not retrieve aodaodV0 %d", i));continue;}
		
		lIndexPosMother     = 0; lIndexNegMother     = 0; lIndexMotherOfMother       = 0;
		lCheckPIdK0Short    = 0; lCheckMcK0Short     = 0; lCheckSecondaryK0s         = 0;
		lCheckPIdLambda     = 0; lCheckMcLambda      = 0; lCheckSecondaryLambda      = 0;
		lCheckPIdAntiLambda = 0; lCheckMcAntiLambda  = 0; lCheckSecondaryAntiLambda  = 0;       
		lComeFromSigma      = -1;lCheckGamma = 0;
		
        // get daughters
		
   	    AliAODTrack *myTrackPos=(AliAODTrack *)(aodV0->GetDaughter(0));
        AliAODTrack *myTrackNeg=(AliAODTrack *)(aodV0->GetDaughter(1));
		
		if (!myTrackPos || !myTrackNeg) {Printf("ERROR: Could not retreive one of the daughter track");continue;}
		
        if (!IsAcseptedV0(fAODEvent,aodV0,myTrackPos,myTrackNeg)) continue;
		
		// VO's main characteristics to check the reconstruction cuts
		lOnFlyStatus       = aodV0->GetOnFlyStatus();
		lChi2V0            = aodV0->Chi2V0();
		lDcaV0Daughters    = aodV0->DcaV0Daughters();
		lDcaV0ToPrimVertex = aodV0->DcaV0ToPrimVertex();
		lV0cosPointAngle   = aodV0->CosPointingAngle(lPrimaryVtxPosition);
		
		aodV0->GetXYZ(lV0Position);
		
		lV0Radius      = TMath::Sqrt(lV0Position[0]*lV0Position[0]+lV0Position[1]*lV0Position[1]);
		lV0DecayLength = TMath::Sqrt(TMath::Power(lV0Position[0] - lPrimaryVtxPosition[0],2) +
									 TMath::Power(lV0Position[1] - lPrimaryVtxPosition[1],2) +
									 TMath::Power(lV0Position[2] - lPrimaryVtxPosition[2],2 ));
		
		lLabelTrackPos = (UInt_t)TMath::Abs(myTrackPos->GetLabel());
		lLabelTrackNeg = (UInt_t)TMath::Abs(myTrackNeg->GetLabel());
		
		// Daughters Pt and P:
		lPtPos = TMath::Sqrt(lMomPos[0]*lMomPos[0] + lMomPos[1]*lMomPos[1]);
		lPtNeg = TMath::Sqrt(lMomNeg[0]*lMomNeg[0] + lMomNeg[1]*lMomNeg[1]);
		
		lPPos = TMath::Sqrt(lMomPos[0]*lMomPos[0] + lMomPos[1]*lMomPos[1] + lMomPos[2]*lMomPos[2]);
		lPNeg = TMath::Sqrt(lMomNeg[0]*lMomNeg[0] + lMomNeg[1]*lMomNeg[1] + lMomNeg[2]*lMomNeg[2]);
		
		// V0 momentum
		lPosMom = lPPos;
		lNegMom = lPNeg;
		
		// Inner Wall parameter:
		const AliAODPid *pidPos=myTrackPos->GetDetPid();      
		const AliAODPid *pidNeg=myTrackNeg->GetDetPid();
		
		// innerWall momentum
		lMomInnerWallPos = pidPos->GetTPCmomentum(); 
		lMomInnerWallNeg = pidNeg->GetTPCmomentum();
		
		ldEdxPos = pidPos->GetTPCsignal();
		ldEdxNeg = pidNeg->GetTPCsignal();
		
		// DCA between daughter and Primary Vertex:
		if (myTrackPos) lDcaPosToPrimVertex = aodV0->DcaPosToPrimVertex();
		if (myTrackNeg) lDcaNegToPrimVertex = aodV0->DcaNegToPrimVertex();      
		
		// Quality tracks cuts:
		if ( !(IsAcseptedDaughterTrack(myTrackPos)) || !(IsAcseptedDaughterTrack(myTrackNeg)) ) { continue;}
		
		// Armenteros variables:
		lAlphaV0      =  aodV0->AlphaV0();
		lPtArmV0      =  aodV0->PtArmV0();
		
		// Pseudorapidity:
		lV0Eta = aodV0->PseudoRapV0();
		//////////////////////////////////////////////////////////////////////////
		// Invariant mass
		lInvMassK0 = aodV0->MassK0Short();
		lPtK0s = aodV0->Pt();
		lPhiK0s= aodV0->Phi();
		lEtaK0s= aodV0->Eta();
		lPzK0s = aodV0->Pz();
		
		lInvMassLambda = aodV0->MassLambda();
		lPtLambda = aodV0->Pt();
		lPhiLambda= aodV0->Phi();
		lEtaLambda= aodV0->Eta();
		lPzLambda = aodV0->Pz();
		
		lInvMassAntiLambda = aodV0->MassAntiLambda();
		lPtAntiLambda = aodV0->Pt();
		lPhiAntiLambda= aodV0->Phi();
		lEtaAntiLambda= aodV0->Eta();
		lPzAntiLambda = aodV0->Pz();
		
		// Rapidity:
		lRapK0s    = aodV0->RapK0Short();
		lRapLambda = aodV0->RapLambda();
		lRapAntiLambda = aodV0->Y(-3122);
		
		if (lPtK0s==0) {continue;}
		if (lPtLambda==0) {continue;}
		if (lPtAntiLambda==0) {continue;}
		
		
		// PID  new method July 2011
		if (fUsePID=="withPID") {
			nSigmaPosPion =	TMath::Abs(IsAccepteddEdx(lPPos,ldEdxPos,AliPID::kPion));
			nSigmaNegPion =	TMath::Abs(IsAccepteddEdx(lPNeg,ldEdxNeg,AliPID::kPion));                              
			nSigmaPosProton = TMath::Abs(IsAccepteddEdx(lPPos,ldEdxPos,AliPID::kProton));
			nSigmaNegProton = TMath::Abs(IsAccepteddEdx(lPNeg,ldEdxNeg,AliPID::kProton));
		}
		else {nSigmaPosPion = 0; nSigmaNegPion =0; nSigmaPosProton = 0; nSigmaNegProton= 0;}
		
		// Monte-Carlo particle associated to reconstructed particles: 
		if (fAnalysisMC) {
			
			AliAODMCParticle *pp=(AliAODMCParticle*)stack->UncheckedAt(lLabelTrackPos);
			if(!pp) { continue;}
			AliAODMCParticle *np=(AliAODMCParticle*)stack->UncheckedAt(lLabelTrackNeg);
			if (!np) 	{ continue;}
			
			lPDGCodePosDaughter = pp->GetPdgCode();
			lPDGCodeNegDaughter = np->GetPdgCode();
			lIndexPosMother = pp->GetMother(); 
			lIndexNegMother = np->GetMother(); 
			
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
				AliAODMCParticle *lMCAODMother=(AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother);
				if (!lMCAODMother) 	{ continue;}
				lPdgcodeMother         = lMCAODMother->GetPdgCode();
				lIndexMotherOfMother   = lMCAODMother->GetMother();
				if (lIndexMotherOfMother ==-1) lPdgcodeMotherOfMother = 0;
				else {
					AliAODMCParticle *lMCAODMotherOfMother=(AliAODMCParticle*)stack->UncheckedAt(lIndexMotherOfMother);
					if (!lMCAODMotherOfMother) 	{continue;}
					lPdgcodeMotherOfMother = lMCAODMotherOfMother->GetPdgCode();
				}
				
				mcPosX = pp->Xv();
				mcPosY = pp->Yv();
				mcPosZ = pp->Zv();
				mcPosR = TMath::Sqrt(mcPosX*mcPosX+mcPosY*mcPosY);
				mcPosMotherX = lMCAODMother->Xv();
				mcPosMotherY = lMCAODMother->Yv();
				mcPosMotherZ = lMCAODMother->Zv();
				mcPosMotherR = TMath::Sqrt(mcPosMotherX*mcPosMotherX+mcPosMotherY*mcPosMotherY);
				
				mcMotherPt   = lMCAODMother->Pt();
			}
		}
		
	}
	
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
				if (((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary()) lCheckMcK0Short  = 1;
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
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3224) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3214) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3114)
					) lComeFromSigma = 1;
				else lComeFromSigma = 0; 
				if ( ((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary() || 
					(!(((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary() ) && (lComeFromSigma) )
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
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3224) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3214) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3114)
					) lComeFromSigma = 1;
				else lComeFromSigma = 0;  
				if ( ((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary() || 
					( (!((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary()) && (lComeFromSigma) )
					) lCheckMcAntiLambda  = 1;
				else lCheckSecondaryAntiLambda = 1;
			}
		}
		
		// Gamma conversion
		else if ((lPDGCodePosDaughter==-11) &&
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
	
	//************************************filling histograms********************************//
    
	if((fUsePID=="withPID") && (lCheckPIDAntiLambdaNegDaughter==0) && lCheckPIDLambdaPosDaughter==1) LambdaPID = 1;
	else LambdaPID =0;
	if((fUsePID=="withPID") && (lCheckPIDLambdaPosDaughter==0) && lCheckPIDAntiLambdaNegDaughter==1) AntiLambdaPID = 1;
	else AntiLambdaPID =0;
	
	lPLambda = TMath::Sqrt(lPzLambda*lPzLambda + lPtLambda*lPtLambda);
	lPAntiLambda = TMath::Sqrt(lPzAntiLambda*lPzAntiLambda + lPtAntiLambda*lPtAntiLambda);
	lPK0s = TMath::Sqrt(lPzK0s*lPzK0s + lPtK0s*lPtK0s);
	
	 
	lcTauLambda     = (lV0DecayLength*lInvMassLambda)/lPLambda;
	lcTauAntiLambda = (lV0DecayLength*lInvMassAntiLambda)/lPAntiLambda; 
	lcTauK0s        = (lV0DecayLength*lInvMassK0)/lPK0s;
	
	if (lPLambda <1 && lOnFlyStatus==0 ){
        fHistcTauL->Fill(lcTauLambda);
	}
	
	//--------------------------------------------K0s---------------------------------//
    
	if (lcTauK0s< cutcTauK0){
		if (TMath::Abs(lRapK0s) < fRapidityCut ){
			if(lPtArmV0*fSpecialArmenterosCutK0s>(TMath::Abs(lAlphaV0))){
			
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
				fHistMassK0->Fill(lInvMassK0);
				fHistMassVsRadiusK0->Fill(rcPosRK0s,lInvMassK0);
				fHistPtVsMassK0->Fill(lInvMassK0,lPtK0s);
				fHistDcaPosToPrimVertexK0vsMassK0->Fill(lDcaPosToPrimVertex,lInvMassK0);
				fHistDcaNegToPrimVertexK0vsMassK0->Fill(lDcaNegToPrimVertex,lInvMassK0);
				fHistRadiusV0K0vsMassK0->Fill(lV0Radius,lInvMassK0);
				fHistDecayLengthV0K0vsMassK0->Fill(lV0DecayLength,lInvMassK0);
				fHistDcaV0DaughtersK0vsMassK0->Fill(lDcaV0Daughters,lInvMassK0);
				fHistCosPointAngleK0vsMassK0->Fill(lV0cosPointAngle,lInvMassK0);
				fHistArmenterosPodolanski->Fill(lAlphaV0,lPtArmV0);
				if(IsK0InvMass(lInvMassK0)){selectedK0->Add(new AliLeadingBasicParticle(lEtaK0s,lPhiK0s,lPtK0s));}
			}
		  }//Special Armesto Cut for K0	
		} // if rap. condition
	} // end cTau condition
	
	//-----------------------------------------Lambda---------------------------------//
	if (lcTauLambda < cutcTauL){
		if ((LambdaPID==1 && lMomInnerWallPos  <=1 ) || (lMomInnerWallPos >1 ) ||  !(fUsePID=="withPID")){ 	
			if (TMath::Abs(lRapLambda) <fRapidityCut) {
				
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
					
					fHistDcaPosToPrimVertexLvsMassL->Fill(lDcaPosToPrimVertex,lInvMassLambda);
					fHistDcaNegToPrimVertexLvsMassL->Fill(lDcaNegToPrimVertex,lInvMassLambda);
					fHistRadiusV0LvsMassL->Fill(lV0Radius,lInvMassLambda);
					fHistDecayLengthV0LvsMassL->Fill(lV0DecayLength,lInvMassLambda);
					fHistDcaV0DaughtersLvsMassL->Fill(lDcaV0Daughters,lInvMassLambda);
					fHistCosPointAngleLvsMassL->Fill(lV0cosPointAngle,lInvMassLambda);
					if(IsLambdaInvMass(lInvMassLambda)){selectedLambda->Add(new AliLeadingBasicParticle(lEtaLambda,lPhiLambda,lPtLambda));}
				}
			} //end of Rap condition
		}
	}  //end cTau condition
	
	//--------------------------------------AntiLambda---------------------------------//
	
	if (lcTauAntiLambda < cutcTauL){
		
		if ((AntiLambdaPID==1 && lMomInnerWallNeg <=1) || (lMomInnerWallNeg >1) ||  !(fUsePID=="withPID")){  
			if (TMath::Abs(lRapAntiLambda) < fRapidityCut) {
				
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
				}
			} //end of Rap condition
		} // end of PID condition
	} //end cTau condition

	//--------------------------------------K0s Associated---------------------------------//

	if (lcTauK0s< cutcTauK0) {
		if (TMath::Abs(lRapK0s) < fRapidityCut) {
			fHistNsigmaPosPionK0->Fill(nSigmaPosPion);
			fHistNsigmaNegPionK0->Fill(nSigmaNegPion);
			if(lOnFlyStatus==0){
					if(lCheckPIdK0Short) fHistPidMcMassK0->Fill(lInvMassK0);
					if(lCheckMcK0Short) {
						fHistAsMcMassK0->Fill(lInvMassK0);
						fHistAsMcPtK0->Fill(lPtK0s);
						fHistAsMcPtVsMassK0->Fill(lInvMassK0,lPtK0s);
						fHistAsMcMassVsRadiusK0->Fill(rcPosRK0s,lInvMassK0);
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
					if (lCheckSecondaryK0s) {
						fHistAsMcSecondaryPtVsRapK0s->Fill(lPtK0s,lRapK0s);
						fHistAsMcSecondaryProdRadiusK0s->Fill(mcPosMotherR);
						fHistAsMcSecondaryProdRadiusXvsYK0s->Fill(mcPosMotherX,mcPosMotherY);
					}
			}
		} // end rapidity condition
	} //end cTau condition
    
	
	//-----------------------------------Lambda Associated---------------------------------//
	if (lcTauLambda < cutcTauL){                                                                                                   
		if (TMath::Abs(lRapLambda) < fRapidityCut) {
			fHistNsigmaPosProtonLambda->Fill(nSigmaPosProton);
			fHistNsigmaNegPionLambda->Fill(nSigmaNegPion);
			if(lOnFlyStatus==0){
					if(lCheckPIdLambda) fHistPidMcMassLambda->Fill(lInvMassLambda);
					if(lCheckMcLambda) {
						fHistAsMcMassLambda->Fill(lInvMassLambda);
						fHistAsMcPtLambda->Fill(lPtLambda);
						fHistCosPointAngleLVsMassVsPtsigL->Fill(lPtLambda,lV0cosPointAngle,lInvMassLambda);
						fHistAsMcPtVsMassLambda->Fill(lInvMassLambda,lPtLambda);
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
					}
					
					if (lCheckSecondaryLambda) {
						fHistAsMcSecondaryPtVsRapLambda->Fill(lPtLambda,lRapLambda);
						fHistAsMcSecondaryProdRadiusLambda->Fill(mcPosMotherR); 
						fHistAsMcSecondaryProdRadiusXvsYLambda->Fill(mcPosMotherX,mcPosMotherY);
						if (lComeFromSigma) fHistAsMcSecondaryPtLambdaFromSigma->Fill(lPtLambda);
					}
					if(!lCheckMcLambda)fHistCosPointAngleLVsMassVsPtbackL->Fill(lPtLambda,lV0cosPointAngle,lInvMassLambda);
			}
		} // end rapidity condition
	}//end cTau condition

//------------------------------AntiLambda Associated---------------------------------//	
	if (lcTauAntiLambda < cutcTauL){
		if (TMath::Abs(lRapAntiLambda) < fRapidityCut) {
			fHistNsigmaPosProtonAntiLambda->Fill(nSigmaPosProton);
			fHistNsigmaNegPionAntiLambda->Fill(nSigmaNegPion);
			if(lOnFlyStatus==0){
					if(lCheckPIdAntiLambda) fHistPidMcMassAntiLambda->Fill(lInvMassAntiLambda);
					if(lCheckMcAntiLambda) {
						fHistAsMcMassAntiLambda->Fill(lInvMassAntiLambda);
						fHistAsMcPtAntiLambda->Fill(lPtAntiLambda);
						fHistAsMcPtVsMassAntiLambda->Fill(lInvMassAntiLambda,lPtAntiLambda);
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
	
					}
					
					if (lCheckSecondaryAntiLambda) {
						fHistAsMcSecondaryPtVsRapAntiLambda->Fill(lPtAntiLambda,lRapAntiLambda);
						fHistAsMcSecondaryProdRadiusAntiLambda->Fill(mcPosMotherR); 
						fHistAsMcSecondaryProdRadiusXvsYAntiLambda->Fill(mcPosMotherX,mcPosMotherY);
						if (lComeFromSigma) fHistAsMcSecondaryPtAntiLambdaFromSigma->Fill(lPtAntiLambda);
					}
			}
		} // end rapidity condition      
	}//end cTau condition
	
	// Correlation part
	if(fAnalysisMC){
	FillCorrelations(selectedTracksLeadingK0MC,selectedK0MC,fHistSibK0MC);
	FillCorrelations(selectedTracksLeadingLambdaMC,selectedLambdaMC,fHistSibLambdaMC);
	}
	
	FillCorrelations(selectedTracksLeadingK0,selectedK0,fHistSibK0);
	FillCorrelations(selectedTracksLeadingLambda,selectedLambda,fHistSibLambda);
	
	// Mixing part
	if(TMath::Abs(lPVz)>=10 || MultipOrCent>poolmax || MultipOrCent < poolmin) {
		if(fcollidingSys=="PP")AliInfo(Form("pp Event with Zvertex = %.2f cm and multiplicity = %.0f out of pool bounds, SKIPPING",lPVz,MultipOrCent));
		if(fcollidingSys=="PbPb") AliInfo(Form("PbPb Event with Zvertex = %.2f cm and centrality = %.1f  out of pool bounds, SKIPPING",lPVz,MultipOrCent));
		return;
	}
	
	fPoolK0 = fPoolMgr->GetEventPool(MultipOrCent, lPVz);
	if (!fPoolK0){AliInfo(Form("No pool found for multiplicity = %f, zVtx = %f cm", MultipOrCent, lPVz));return;}
	
    if (fPoolK0->IsReady() || fPoolK0->NTracksInPool() > fPoolMinNTracks  || fPoolK0->GetCurrentNEvents() >= fMinEventsToMix)
	{
		for (Int_t jMix=0; jMix<fPoolK0->GetCurrentNEvents(); jMix++) 
			FillCorrelations(selectedTracksLeadingK0, fPoolK0->GetEvent(jMix),fHistMixK0);
		fPoolK0->UpdatePool(CloneAndReduceTrackList(selectedK0));
	}
	
	fPoolLambda = fPoolMgr->GetEventPool(MultipOrCent, lPVz);
	if (!fPoolLambda){AliInfo(Form("No pool found for multiplicity = %f, zVtx = %f cm", MultipOrCent, lPVz));return;}
	
    if (fPoolLambda->IsReady() || fPoolLambda->NTracksInPool() > fPoolMinNTracks  || fPoolLambda->GetCurrentNEvents() >= fMinEventsToMix)
	{
		for (Int_t jMix=0; jMix<fPoolLambda->GetCurrentNEvents(); jMix++) 
			FillCorrelations(selectedTracksLeadingLambda, fPoolLambda->GetEvent(jMix),fHistMixLambda);
		fPoolLambda->UpdatePool(CloneAndReduceTrackList(selectedLambda));
	}
	
	PostData(1, fOutputList);

}		
//---------------------------------------------------------------------------------------
TObjArray* AliLeadingV0Correlation::CloneAndReduceTrackList(TObjArray* tracks)
{
	// clones a track list for mixing
	TObjArray* tracksClone = new TObjArray;
	tracksClone->SetOwner(kTRUE);
	
	for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
	{
		AliVParticle* particle = (AliVParticle*) tracks->At(i);
		tracksClone->Add(new AliLeadingBasicParticle(particle->Eta(), particle->Phi(), particle->Pt()));
	}
	
	return tracksClone;
}
//---------------------------------------------------------------------------------------
Int_t  AliLeadingV0Correlation::NParticles(TObject* obj)
{
	Int_t nTracks;
	
	if (obj->InheritsFrom("TClonesArray")){ // MC particles
		TClonesArray *arrayMC = static_cast<TClonesArray*>(obj);
        nTracks = arrayMC->GetEntriesFast();
	}else if (obj->InheritsFrom("TObjArray")){ // list of AliVParticle
		TObjArray *array = static_cast<TObjArray*>(obj);
        nTracks = array->GetEntriesFast();
	}else if (obj->InheritsFrom("AliAODEvent")){  // RECO AOD tracks
		AliAODEvent *aodEvent = static_cast<AliAODEvent*>(obj);
        nTracks = aodEvent->GetNTracks();
	}else if (obj->InheritsFrom("AliMCEvent")){  // RECO ESD tracks
		AliMCEvent *mcEvent = static_cast<AliMCEvent*>(obj);
        nTracks = mcEvent->GetNumberOfTracks();
	}else {
		if (fDebug > 1) AliFatal(" Analysis type not defined !!! ");
		return 0;
	}
	
	return nTracks;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsAcseptedPrimaryTrack(const AliAODTrack *itrack)
{
	if (TMath::Abs(itrack->Eta())>fTrackEtaCut) return kFALSE;
    if (!itrack->TestFilterBit(fFilterBit)) return kFALSE;	
	return kTRUE;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsAcseptedDaughterTrack(const AliAODTrack *itrack)
{
	if(TMath::Abs(itrack->Eta())>fTrackEtaCut)return kFALSE;
	
	if (!itrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
	
	Float_t nCrossedRowsTPC = itrack->GetTPCClusterInfo(2,1);
	if (nCrossedRowsTPC < 70) return kFALSE;
	
	Int_t findable=itrack->GetTPCNclsF();
	if (findable <= 0) return kFALSE;
	
	if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
	
	return kTRUE;
}
//---------------------------------------------------------------------------------------
Double_t AliLeadingV0Correlation::RangePhi(Double_t DPhi)
{
	if (DPhi < -TMath::Pi()/2)  DPhi += 2*TMath::Pi();
	if (DPhi > 3*TMath::Pi()/2) DPhi -= 2*TMath::Pi();	
	return DPhi;	
}
//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::FillCorrelations(TObjArray* particles, TObjArray* mixed,TH2F*histo)
{
 TObjArray* input = (mixed) ? mixed : particles;
 TArrayF eta(input->GetEntriesFast());
   for (Int_t i=0; i<input->GetEntriesFast(); i++) eta[i] = ((AliVParticle*) input->At(i))->Eta();
		if (particles)
		 {
		   Int_t jMax = particles->GetEntriesFast();
		   if (mixed) jMax = mixed->GetEntriesFast();
			for (Int_t i=0; i<particles->GetEntriesFast(); i++)
			{
				AliVParticle* triggerParticle = (AliVParticle*) particles->At(0); //For leading Track
				// some optimization
				Float_t triggerEta = triggerParticle->Eta();
				
				for (Int_t j=0; j<jMax; j++)
				{
				  if (!mixed && i == j)continue;
					AliVParticle* particle = 0;
					if (!mixed)particle = (AliVParticle*) particles->At(0); //For leading Track
					else particle = (AliVParticle*) mixed->At(j);
					
					// check if both particles point to the same element (does not occur for mixed events, but if subsets are mixed within the same event for cross-checks)
					if (mixed && triggerParticle == particle)continue;
					
					if (particle->Pt() > triggerParticle->Pt())continue;
					
					if (((particle->Pt())<fassocPtLow)||((particle->Pt())>fassocPtHigh)) continue;
					if (((triggerParticle->Pt())<ftrigPtLow)||((triggerParticle->Pt())>ftrigPtHigh)) continue;
					
					Double_t vars[4];
					vars[0] = triggerEta - eta[j];    
					vars[1] = particle->Pt();         //Associated Pt
					vars[2] = triggerParticle->Pt();  //Triger Pt
					vars[3] = RangePhi(triggerParticle->Phi() - particle->Phi());
					
					histo->Fill(vars[0],vars[3]);

	           }
	        }
        }
}
//---------------------------------------------------------------------------------------
TObjArray*  AliLeadingV0Correlation::FindLeadingObjectsK0(TObject *obj)
{

	Int_t nTracks = NParticles(obj);
	if( !nTracks ) return 0;
	
	// Define array of AliVParticle objects
	TObjArray* tracks = new TObjArray(nTracks);
	
	// Loop over tracks or jets
	for (Int_t ipart=0; ipart<nTracks; ++ipart) {
        AliVParticle* part = ParticleWithCuts( obj,ipart);
        if (!part) continue;
		
		if(!IsAcseptedPrimaryTrack(((AliAODTrack*)part)))continue;
		if(!IsTrackNotFromK0(ipart))continue;
		
		tracks->AddLast( part );
  	}
	// Order tracks by pT	
	QSortTracks( *tracks, 0, tracks->GetEntriesFast() );
	
	nTracks = tracks->GetEntriesFast();
	if( !nTracks ) return 0;
	
	return tracks;
}
//---------------------------------------------------------------------------------------
TObjArray*  AliLeadingV0Correlation::FindLeadingObjectsLambda(TObject *obj)
{
	
	Int_t nTracks = NParticles(obj);
	if( !nTracks ) return 0;
	
	// Define array of AliVParticle objects
	TObjArray* tracks = new TObjArray(nTracks);
	
	// Loop over tracks or jets
	for (Int_t ipart=0; ipart<nTracks; ++ipart) {
        AliVParticle* part = ParticleWithCuts( obj,ipart);
        if (!part) continue;
		
		if(!IsAcseptedPrimaryTrack(((AliAODTrack*)part)))continue;
		if(!IsTrackNotFromLambda(ipart))continue;
		
		tracks->AddLast( part );
  	}
	// Order tracks by pT	
	QSortTracks( *tracks, 0, tracks->GetEntriesFast() );
	
	nTracks = tracks->GetEntriesFast();
	if( !nTracks ) return 0;
	
	return tracks;
}
//---------------------------------------------------------------------------------------
TObjArray*  AliLeadingV0Correlation::FindLeadingObjectsK0MC(TObject *obj)
{
	
	Int_t nTracks = NParticles(obj);
	if( !nTracks ) return 0;
	
	// Define array of AliVParticle objects
	TObjArray* tracks = new TObjArray(nTracks);
	
	// Loop over tracks or jets
	for (Int_t ipart=0; ipart<nTracks; ++ipart) {
        AliVParticle* part = ParticleWithCuts( obj,ipart);
        if (!part) continue;
		Int_t pdgCodeCurrent = ((AliAODMCParticle*)part)->GetPdgCode();
		if(pdgCodeCurrent==310)continue;
		tracks->AddLast( part );
  	}
	// Order tracks by pT	
	QSortTracks( *tracks, 0, tracks->GetEntriesFast() );
	
	nTracks = tracks->GetEntriesFast();
	if( !nTracks ) return 0;
	
	return tracks;
}
//---------------------------------------------------------------------------------------
TObjArray*  AliLeadingV0Correlation::FindLeadingObjectsLambdaMC(TObject *obj)
{
	
	Int_t nTracks = NParticles(obj);
	if( !nTracks ) return 0;
	
	// Define array of AliVParticle objects
	TObjArray* tracks = new TObjArray(nTracks);
	
	// Loop over tracks or jets
	for (Int_t ipart=0; ipart<nTracks; ++ipart) {
        AliVParticle* part = ParticleWithCuts( obj,ipart);
        if (!part) continue;
		Int_t pdgCodeCurrent = ((AliAODMCParticle*)part)->GetPdgCode();
		if(pdgCodeCurrent==3122)continue;
		tracks->AddLast( part );
  	}
	// Order tracks by pT	
	QSortTracks( *tracks, 0, tracks->GetEntriesFast() );
	
	nTracks = tracks->GetEntriesFast();
	if( !nTracks ) return 0;
	
	return tracks;
}
//---------------------------------------------------------------------------------------
void  AliLeadingV0Correlation::QSortTracks(TObjArray &a, Int_t first, Int_t last)
{
	// Sort array of TObjArray of tracks by Pt using a quicksort algorithm.
	
	static TObject *tmp;
	static int i;           // "static" to save stack space
	int j;
	
	while (last - first > 1) {
		i = first;
		j = last;
		for (;;) {
			while (++i < last && ((AliVParticle*)a[i])->Pt() > ((AliVParticle*)a[first])->Pt() )
				;
			while (--j > first && ((AliVParticle*)a[j])->Pt() < ((AliVParticle*)a[first])->Pt() )
				;
			if (i >= j)
				break;
			
			tmp  = a[i];
			a[i] = a[j];
			a[j] = tmp;
		}
		if (j == first) {
			++first;
			continue;
		}
		tmp = a[first];
		a[first] = a[j];
		a[j] = tmp;
		if (j - first < last - (j + 1)) {
			QSortTracks(a, first, j);
			first = j + 1;   // QSortTracks(j + 1, last);
		} else {
			QSortTracks(a, j + 1, last);
			last = j;        // QSortTracks(first, j);
		}
	}
}
//---------------------------------------------------------------------------------------
AliVParticle*  AliLeadingV0Correlation::ParticleWithCuts(TObject* obj, Int_t ipart)
{   
	AliVParticle *part=0;
	if (obj->InheritsFrom("AliAODEvent")){ // RECO AOD TRACKS
		AliAODEvent *aodEvent = static_cast<AliAODEvent*>(obj);
        part = aodEvent->GetTrack(ipart);
		// track selection cuts
		if ( !(((AliAODTrack*)part)->TestFilterBit(fFilterBit))) return 0; 
	}
	else if(obj->InheritsFrom("TClonesArray")){ // AOD-MC PARTICLE
		TClonesArray *arrayMC = static_cast<TClonesArray*>(obj);
        part = (AliVParticle*)arrayMC->At( ipart );
		if (!part)return 0;
		// eventually only primaries
		if (!( ((AliAODMCParticle*)part)->IsPhysicalPrimary()) )return 0;
        // eventually only hadrons
	        Int_t pdgCode = ((AliAODMCParticle*)part)->GetPdgCode();
			Bool_t isHadron = TMath::Abs(pdgCode)==211 ||  // Pion
			TMath::Abs(pdgCode)==2212 || // Proton
			TMath::Abs(pdgCode)==321;    // Kaon
			if (!isHadron) return 0;				  
	}
	else {return 0;}
	// only charged
	if (!part->Charge())return 0;
	return part;
}
//---------------------------------------------------------------------------------------
Bool_t   AliLeadingV0Correlation::IsK0InvMass(const Double_t mass) const {
	
	const Float_t massK0            = 0.497;
	const Float_t sigmaK0           = 0.003;
	const Float_t nSigmaSignal      = 3.5; 
	
	return ((massK0-nSigmaSignal*sigmaK0)<=mass && mass<=(massK0 + nSigmaSignal*sigmaK0))?1:0;
}
//---------------------------------------------------------------------------------------
Bool_t   AliLeadingV0Correlation::IsLambdaInvMass(const Double_t mass) const {
	
	const Float_t massLambda        = 1.116;
	const Float_t sigmaLambda       = 0.003;
	const Float_t nSigmaSignal      = 3.5; 
	
	return ((massLambda-nSigmaSignal*sigmaLambda)<=mass && mass<=(massLambda + nSigmaSignal*sigmaLambda))?1:0;
}
//---------------------------------------------------------------------------------------
Bool_t   AliLeadingV0Correlation::IsTrackNotFromLambda(const Int_t indexTrack){
	// check wether track with index is from Lambda (= V0 with proper inv mass)
	TClonesArray* tracks = fAODEvent->GetTracks();
	for(int i=0; i<fAODEvent->GetNumberOfV0s(); i++){ // loop over V0s
		AliAODv0* aodV0 = fAODEvent->GetV0(i);
		Float_t massLambda         = aodV0->MassLambda();
		if(IsLambdaInvMass(massLambda)){
			AliAODTrack *trackPos = (AliAODTrack *) (aodV0->GetSecondaryVtx()->GetDaughter(0));
			AliAODTrack *trackNeg = (AliAODTrack *) (aodV0->GetSecondaryVtx()->GetDaughter(1));
			
			if ( !(IsAcseptedDaughterTrack(trackPos)) || !(IsAcseptedDaughterTrack(trackNeg)) ) continue;
			
			Int_t indexPos = tracks->IndexOf(trackPos);
			Int_t indexNeg = tracks->IndexOf(trackNeg);
			if(indexPos == indexTrack){ return kFALSE;}
			if(indexNeg == indexTrack){ return kFALSE;}
		}
	}
	return kTRUE;
}
//---------------------------------------------------------------------------------------
Bool_t   AliLeadingV0Correlation::IsTrackNotFromK0(const Int_t indexTrack){
	// check wether track with index is from K0 (= V0 with proper inv mass)
	TClonesArray* tracks = fAODEvent->GetTracks();
	for(int i=0; i<fAODEvent->GetNumberOfV0s(); i++){ // loop over V0s
		AliAODv0* aodV0 = fAODEvent->GetV0(i);
		Float_t massK0         = aodV0->MassK0Short();
		if(IsK0InvMass(massK0)){
			AliAODTrack *trackPos = (AliAODTrack *) (aodV0->GetSecondaryVtx()->GetDaughter(0));
			AliAODTrack *trackNeg = (AliAODTrack *) (aodV0->GetSecondaryVtx()->GetDaughter(1));
			
			if ( !(IsAcseptedDaughterTrack(trackPos)) || !(IsAcseptedDaughterTrack(trackNeg)) ) continue;
			
			Int_t indexPos = tracks->IndexOf(trackPos);
			Int_t indexNeg = tracks->IndexOf(trackNeg);
			if(indexPos == indexTrack){ return kFALSE;}
		    if(indexNeg == indexTrack){ return kFALSE;}
		}
	}
	return kTRUE;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsAcseptedV0(const AliAODEvent*aod, const AliAODv0* aodV0, const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg)
{
	if (!aodV0) return kFALSE;
	
	Double_t lRapK0s = aodV0->Y(310);
	Double_t lRapLambda = aodV0->Y(3122);
	Double_t lRapAntiLambda = aodV0->Y(-3122);
	
	if (TMath::Abs(lRapK0s)>=fRapidityCut) return kFALSE;
	if (TMath::Abs(lRapLambda)>=fRapidityCut) return kFALSE;
	if (TMath::Abs(lRapAntiLambda)>=fRapidityCut) return kFALSE;
    
	// Offline reconstructed V0 only
    if (aodV0->GetOnFlyStatus()) return kFALSE;
	
    // DCA of daughter track to Primary Vertex
    Float_t dcaNP=aodV0->DcaNegToPrimVertex();
    if (TMath::Abs(dcaNP)<fCutDCANegToPV) return kFALSE;
    Float_t dcaPP=aodV0->DcaPosToPrimVertex();
    if (TMath::Abs(dcaPP)<fCutDCAPosToPV) return kFALSE;
	
	// DCA of daughter tracks 
    Double_t dcaD=aodV0->DcaV0Daughters();
    if (dcaD>fCutDCAV0Daughters) return kFALSE;
	
	// Cosinus of pointing angle
    Double_t cpa=aodV0->CosPointingAngle(aod->GetPrimaryVertex());
    if (cpa<fCutV0CosPA) return kFALSE;
	
	// Decay Radius Cut
    Double_t xyz[3]; aodV0->GetSecondaryVtx(xyz);
    Double_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];
    if (r2<fCutV0Radius) return kFALSE;
	
    // Get daughters and check them
	myTrackPos=(AliAODTrack *)(aodV0->GetDaughter(0));
	myTrackNeg=(AliAODTrack *)(aodV0->GetDaughter(1));
	
	if (!myTrackPos||!myTrackNeg) return kFALSE;
	// Unlike signs of daughters
    if (myTrackPos->Charge() == myTrackNeg->Charge()) return kFALSE;
	
	// Track cuts for daughers
    if ( !(IsAcseptedDaughterTrack(myTrackPos)) || !(IsAcseptedDaughterTrack(myTrackNeg)) ) return kFALSE;
	
	// Minimum pt of daughters
    Double_t  lMomPos[3] = {999,999,999};
	Double_t  lMomNeg[3] = {999,999,999};
	
	lMomPos[0] = aodV0->MomPosX();
    lMomPos[1] = aodV0->MomPosY();
    lMomPos[2] = aodV0->MomPosZ();
	
    lMomNeg[0] = aodV0->MomNegX();
    lMomNeg[1] = aodV0->MomNegY();
    lMomNeg[2] = aodV0->MomNegZ();
	
    Double_t lPtPos = TMath::Sqrt(lMomPos[0]*lMomPos[0] + lMomPos[1]*lMomPos[1]);
    Double_t lPtNeg = TMath::Sqrt(lMomNeg[0]*lMomNeg[0] + lMomNeg[1]*lMomNeg[1]);
	
	Double_t cutMinPtDaughter = 0.150;
	if (lPtPos<cutMinPtDaughter || lPtNeg<cutMinPtDaughter) return kFALSE;
	
	return kTRUE;
}
//---------------------------------------------------------------------------------------
Double_t AliLeadingV0Correlation::IsAccepteddEdx(const Double_t mom,const Double_t signal, AliPID::EParticleType n){
	
	const Double_t kBBMIP(50.);
	const Double_t kBBRes(0.07);
	const Double_t kBBp1(0.76176e-1);
	const Double_t kBBp2(10.632);
	const Double_t kBBp3(0.13279e-4);
	const Double_t kBBp4(1.8631);
	const Double_t kBBp5(1.9479);
	
	Double_t mass=AliPID::ParticleMass(n); 
	Double_t betaGamma = mom/mass;
	
	const Float_t kmeanCorrection =0.1;
	Double_t bb = AliExternalTrackParam::BetheBlochAleph(betaGamma,kBBp1,kBBp2,kBBp3,kBBp4,kBBp5);
	Double_t meanCorrection =(1+(bb-1)*kmeanCorrection);
	Double_t bethe = bb * meanCorrection;
	Double_t sigma = bethe * kBBRes;
	
	Double_t dedx = signal/kBBMIP;
	Double_t nSig = (TMath::Abs(dedx - bethe))/sigma;
	
	return nSig;
}

//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::Terminate(Option_t *)
{
	//No need in the grid
}
//---------------------------------------------------------------------------------------
