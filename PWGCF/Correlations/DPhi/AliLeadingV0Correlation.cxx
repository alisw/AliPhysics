/* Leading Charged Track+V0 Correlation.(Works for Real,Monte Carlo Data)
 *                            Sandun Jayarathna
 *                          University of Houston
 *                      sandun.pahula.hewage@cern.ch
 *****************************************************************************************/
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
#include "AliAODVertex.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliEventPoolManager.h"
#include "AliCentrality.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"
#include "AliVParticle.h"
#include "AliMultiplicity.h"
#include "AliAODMCHeader.h"
#include "AliPID.h"
#include "AliExternalTrackParam.h"

#include "AliLeadingV0Correlation.h"

#define CorrBinsX 24
#define CorrBinsY 26

ClassImp(AliLeadingV0Correlation)
ClassImp(AliLeadingBasicParticle)

//---------------------------------------------------------------------------------------
AliLeadingV0Correlation::AliLeadingV0Correlation()
   : AliAnalysisTaskSE(),
	fAODEvent(0x0),
	fPoolMgr(0x0),
	fPoolMgrMC(0x0),
	fPoolK0(0x0),
	fPoolLambda(0x0),
	fPoolK0MC(0x0),
	fPoolLambdaMC(0x0),
	fPIDResponse(0x0),
	fPoolMaxNEvents(0), 
	fPoolMinNTracks(0), 
	fMinEventsToMix(0),
	fNzVtxBins(0), 
	fNCentBins(0),
	fcollidingSys(""),
	fTriggerMask(""),
	fpvzcut(0),
	fTrackEtaCut(0),
	fFilterBit(128),
	ftrigPtLow(0),
	ftrigPtHigh(0),
	fassocPtLow(0),
	fassocPtHigh(0),
	fAnalysisMC(0),
	fUsePID(""),
	fRapidityCut(0),
	fCutV0Radius(0.),
	fCutDCANegToPV(0.),
	fCutDCAPosToPV(0.),
	fCutDCAV0Daughters(0.),
	fCutV0CosPA(0.),
	fSpecialArmenterosCutK0s(0.),
	fCTauK0(0.),
	fCTauLambda(0.),
	fOutputList(0),
	fHistMCPrimaryVertexX(0),
	fHistMCPrimaryVertexY(0),
	fHistMCPrimaryVertexZ(0),
	fHistMCPtAllK0s(0),
	fHistMCPtAllLambda(0),
	fHistMCPtAllAntiLambda(0),
	fHistMCRapK0s(0),
	fHistMCRapLambda(0),
	fHistMCRapAntiLambda(0),
	fHistMCPtK0s(0),
	fHistMCPtLambda(0),
	fHistMCPtAntiLambda(0),
	fHistMCPtLambdaFromSigma(0),
	fHistMCPtAntiLambdaFromSigma(0),
	fHistPrimRawPtVsYK0s(0),
	fHistPrimRawPtVsYLambda(0),
	fHistPrimRawPtVsYAntiLambda(0),
	fHistPrimaryVertexX(0),
	fHistPrimaryVertexY(0),
	fHistPrimaryVertexZ(0),
	fHistDcaPosToPrimVertexK0vsMassK0(0),
	fHistDcaNegToPrimVertexK0vsMassK0(0),
	fHistRadiusV0K0vsMassK0(0),
	fHistDecayLengthV0K0vsMassK0(0),
	fHistDcaV0DaughtersK0vsMassK0(0),
	fHistCosPointAngleK0vsMassK0(0),
	fHistDcaPosToPrimVertexLvsMassL(0),
	fHistDcaNegToPrimVertexLvsMassL(0),
	fHistRadiusV0LvsMassL(0),
	fHistDecayLengthV0LvsMassL(0),
	fHistDcaV0DaughtersLvsMassL(0),
	fHistCosPointAngleLvsMassL(0),
	fHistDcaPosToPrimVertexAntiLvsMass(0),
	fHistDcaNegToPrimVertexAntiLvsMass(0),
	fHistRadiusV0AntiLvsMass(0),
	fHistDecayLengthV0AntiLvsMass(0),
	fHistDcaV0DaughtersAntiLvsMass(0),
	fHistCosPointAngleAntiLvsMass(0),
	fHistMassK0(0),
	fHistMassLambda(0),
	fHistMassAntiLambda(0),
	fHistPtVsMassK0(0),
	fHistPtVsMassLambda(0),
	fHistPtVsMassAntiLambda(0),
	fHistArmenterosPodolanskiK0(0),
	fHistArmenterosPodolanskiLambda(0),
	fHistArmenterosPodolanskiAntiLambda(0),
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
	fHistMixLambdaMC(0),
    fHistLeadInfo(0),
    fHistLeadInfoMC(0),
	fHistLeadInfoMix(0),
	fHistLeadInfoMixMC(0)
{

  for(Int_t iBin = 0; iBin < 100; iBin++){
    fZvtxBins[iBin] = 0.;
    fCentBins[iBin] = 0.;
  }

}
//---------------------------------------------------------------------------------------
AliLeadingV0Correlation::AliLeadingV0Correlation(const char *name)
   : AliAnalysisTaskSE(name),
	fAODEvent(0x0),
	fPoolMgr(0x0),
	fPoolMgrMC(0x0),
	fPoolK0(0x0),
	fPoolLambda(0x0),
	fPoolK0MC(0x0),
	fPoolLambdaMC(0x0),
    fPIDResponse(0x0),
	fPoolMaxNEvents(0), 
	fPoolMinNTracks(0), 
	fMinEventsToMix(0),
	fNzVtxBins(0), 
	fNCentBins(0),
	fcollidingSys(""),
	fTriggerMask(""),
	fpvzcut(0),
    fTrackEtaCut(0),
    fFilterBit(128),
	ftrigPtLow(0),
	ftrigPtHigh(0),
	fassocPtLow(0),
	fassocPtHigh(0),
	fAnalysisMC(0),
	fUsePID(""),
    fRapidityCut(0),
	fCutV0Radius(0.),
	fCutDCANegToPV(0.),
	fCutDCAPosToPV(0.),
	fCutDCAV0Daughters(0.),
	fCutV0CosPA(0.),
	fSpecialArmenterosCutK0s(0.),
    fCTauK0(0.),
    fCTauLambda(0.),
	fOutputList(0),
	fHistMCPrimaryVertexX(0),
	fHistMCPrimaryVertexY(0),
	fHistMCPrimaryVertexZ(0),
	fHistMCPtAllK0s(0),
	fHistMCPtAllLambda(0),
	fHistMCPtAllAntiLambda(0),
	fHistMCRapK0s(0),
	fHistMCRapLambda(0),
	fHistMCRapAntiLambda(0),
	fHistMCPtK0s(0),
	fHistMCPtLambda(0),
	fHistMCPtAntiLambda(0),
	fHistMCPtLambdaFromSigma(0),
	fHistMCPtAntiLambdaFromSigma(0),
	fHistPrimRawPtVsYK0s(0),
	fHistPrimRawPtVsYLambda(0),
	fHistPrimRawPtVsYAntiLambda(0),
	fHistPrimaryVertexX(0),
	fHistPrimaryVertexY(0),
	fHistPrimaryVertexZ(0),
	fHistDcaPosToPrimVertexK0vsMassK0(0),
	fHistDcaNegToPrimVertexK0vsMassK0(0),
	fHistRadiusV0K0vsMassK0(0),
	fHistDecayLengthV0K0vsMassK0(0),
	fHistDcaV0DaughtersK0vsMassK0(0),
	fHistCosPointAngleK0vsMassK0(0),
	fHistDcaPosToPrimVertexLvsMassL(0),
	fHistDcaNegToPrimVertexLvsMassL(0),
	fHistRadiusV0LvsMassL(0),
	fHistDecayLengthV0LvsMassL(0),
	fHistDcaV0DaughtersLvsMassL(0),
	fHistCosPointAngleLvsMassL(0),
	fHistDcaPosToPrimVertexAntiLvsMass(0),
	fHistDcaNegToPrimVertexAntiLvsMass(0),
	fHistRadiusV0AntiLvsMass(0),
	fHistDecayLengthV0AntiLvsMass(0),
	fHistDcaV0DaughtersAntiLvsMass(0),
	fHistCosPointAngleAntiLvsMass(0),
	fHistMassK0(0),
	fHistMassLambda(0),
	fHistMassAntiLambda(0),
	fHistPtVsMassK0(0),
	fHistPtVsMassLambda(0),
	fHistPtVsMassAntiLambda(0),
	fHistArmenterosPodolanskiK0(0),
	fHistArmenterosPodolanskiLambda(0),
	fHistArmenterosPodolanskiAntiLambda(0),
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
	fHistMixLambdaMC(0),
	fHistLeadInfo(0),
	fHistLeadInfoMC(0),
	fHistLeadInfoMix(0),
	fHistLeadInfoMixMC(0)
{	

  for(Int_t iBin = 0; iBin < 100; iBin++){
    fZvtxBins[iBin] = 0.;
    fCentBins[iBin] = 0.;
  }
  
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
	
	//---------------------------------------------- MC histograms -----------------------------------------------------//
	fHistMCPrimaryVertexX         = new TH1F("h1MCPrimaryVertexX", "MC Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
	fOutputList->Add(fHistMCPrimaryVertexX);
	
	fHistMCPrimaryVertexY         = new TH1F("h1MCPrimaryVertexY", "MC Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
	fOutputList->Add(fHistMCPrimaryVertexY);
	
	fHistMCPrimaryVertexZ         = new TH1F("h1MCPrimaryVertexZ", "MC Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
	fOutputList->Add(fHistMCPrimaryVertexZ);

	fHistMCRapK0s                 = new TH1F("h1MCRapK0s", "K^{0};y",160,-4,4);
	fOutputList->Add(fHistMCRapK0s);
	
	fHistMCRapLambda              = new TH1F("h1MCRapLambda", "#Lambda;y",160,-4,4);
	fOutputList->Add(fHistMCRapLambda);
	
	fHistMCRapAntiLambda          = new TH1F("h1MCRapAntiLambda", "#bar{#Lambda};y",160,-4,4);
	fOutputList->Add(fHistMCRapAntiLambda);

	fHistMCPtK0s                  = new TH1F("h1MCPtK0s", "K^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtK0s);
	
	fHistMCPtLambda               = new TH1F("h1MCPtLambda", "#Lambda^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtLambda);
	
	fHistMCPtAntiLambda           = new TH1F("h1MCPtAntiLambda", "#AntiLambda^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtAntiLambda);

	fHistMCPtLambdaFromSigma      = new TH1F("h1MCPtLambdaFromSigma", "#Lambda^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtLambdaFromSigma);
	
	fHistMCPtAntiLambdaFromSigma  = new TH1F("h1MCPtAntiLambdaFromSigma", "#Lambda^{0};p{t} (GeV/c)",240,0,12);
	fOutputList->Add(fHistMCPtAntiLambdaFromSigma);
	
	fHistMCPtAllK0s               = new TH1F("h1MCPtAllK0s", "Non-primary MC K^{0};p{t} (GeV/c);Counts",240,0,12);
	fOutputList->Add(fHistMCPtAllK0s);
	
	fHistMCPtAllLambda            = new TH1F("h1MCPtAllLambda", "Non-primary MC #Lambda^{0};p{t} (GeV/c);Counts",240,0,12);
	fOutputList->Add(fHistMCPtAllLambda);
	
	fHistMCPtAllAntiLambda        = new TH1F("h1MCPtAllAntiLambda", "Non-primary MC #bar{#Lambda}^{0};p{t} (GeV/c);Counts",240,0,12);
	fOutputList->Add(fHistMCPtAllAntiLambda);
	
	fHistPrimRawPtVsYK0s          = new TH2F( "f3dHistPrimRawPtVsYVsMultK0Short", "Pt{K0S} Vs Y{K0S} Vs Multiplicity; Pt{K0S} (GeV/c); Y{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2);
	fOutputList->Add(fHistPrimRawPtVsYK0s);

	
	fHistPrimRawPtVsYLambda       = new TH2F( "f3dHistPrimRawPtVsYVsMultLambda", "Pt{lambda} Vs Y{#Lambda} Vs Multiplicity; Pt{lambda} (GeV/c); Y{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2);
	fOutputList->Add(fHistPrimRawPtVsYLambda);
	
	fHistPrimRawPtVsYAntiLambda   = new TH2F( "f3dHistPrimRawPtVsYVsMultAntiLambda", "Pt{antilambda} Vs Y{#Lambda} Vs Multiplicity; Pt{antilambda} (GeV/c); Y{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2);
	fOutputList->Add(fHistPrimRawPtVsYAntiLambda);
	
	//---------------------------------------------End Of MC Histos-----------------------------------------------------//

	fHistPrimaryVertexX           = new TH1F("h1PrimaryVertexX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
	fOutputList->Add(fHistPrimaryVertexX);
	
	fHistPrimaryVertexY           = new TH1F("h1PrimaryVertexY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
	fOutputList->Add(fHistPrimaryVertexY);
	
	fHistPrimaryVertexZ           = new TH1F("h1PrimaryVertexZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-20,20);
	fOutputList->Add(fHistPrimaryVertexZ);

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
	
	fHistDcaPosToPrimVertexLvsMassL   = new TH2F("h2DcaPosToPrimVertexLvsMassL", "Positive V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaPosToPrimVertexLvsMassL);
	
	fHistDcaNegToPrimVertexLvsMassL   = new TH2F("h2DcaNegToPrimVertexLvsMassL", "Negative V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaNegToPrimVertexLvsMassL);
	
	fHistRadiusV0LvsMassL             = new TH2F("h2RadiusV0LvsMassL", "Radius;Radius(cm);Status",110,0,110,140, 1.06, 1.2);
	fOutputList->Add(fHistRadiusV0LvsMassL);
	
	fHistDecayLengthV0LvsMassL        = new TH2F("h2DecayLengthV0LvsMassL", "V0s decay Length;decay length(cm);Status", 120, 0, 120,140, 1.06, 1.2);
	fOutputList->Add(fHistDecayLengthV0LvsMassL);
	
	fHistDcaV0DaughtersLvsMassL        = new TH2F("h2DcaV0DaughtersLvsMassL", "DCA between daughters;dca(cm);Status", 110, 0, 1.1,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaV0DaughtersLvsMassL);
	
	fHistCosPointAngleLvsMassL         = new TH2F("h2CosPointAngleLvsMassL", "Cosine of V0's pointing angle", 200,0.997,1.007,140, 1.06, 1.2);
	fOutputList->Add(fHistCosPointAngleLvsMassL);

	fHistDcaPosToPrimVertexAntiLvsMass   = new TH2F("h2DcaPosToPrimVertexAntiLvsMass", "Positive V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaPosToPrimVertexAntiLvsMass);
	
	fHistDcaNegToPrimVertexAntiLvsMass   = new TH2F("h2DcaNegToPrimVertexAntiLvsMass", "Negative V0 daughter;dca(cm);Status",100,0,10,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaNegToPrimVertexAntiLvsMass);

	fHistRadiusV0AntiLvsMass             = new TH2F("h2RadiusV0AntiLvsMass", "Radius;Radius(cm);Status",110,0,110,140, 1.06, 1.2);
	fOutputList->Add(fHistRadiusV0AntiLvsMass);
	
	fHistDecayLengthV0AntiLvsMass        = new TH2F("h2DecayLengthV0AntiLvsMass", "V0s decay Length;decay length(cm);Status", 120, 0, 120,140, 1.06, 1.2);
	fOutputList->Add(fHistDecayLengthV0AntiLvsMass);
	
	fHistDcaV0DaughtersAntiLvsMass       = new TH2F("h2DcaV0DaughtersAntiLvsMass", "DCA between daughters;dca(cm);Status", 110, 0, 1.1,140, 1.06, 1.2);
	fOutputList->Add(fHistDcaV0DaughtersAntiLvsMass);
	
	fHistCosPointAngleAntiLvsMass        = new TH2F("h2CosPointAngleAntiLvsMass", "Cosine of V0's pointing angle", 200,0.997,1.007,140, 1.06, 1.2);
	fOutputList->Add(fHistCosPointAngleAntiLvsMass);
	
	fHistMassK0                          = new TH1F("h1MassK0", "K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 200, 0.4, 0.6);
	fOutputList->Add(fHistMassK0);
	
	fHistMassLambda                      = new TH1F("h1MassLambda", "#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 150, 1.05, 1.2);
	fOutputList->Add(fHistMassLambda);
	
	fHistMassAntiLambda                  = new TH1F("h1MassAntiLambda", "#bar{#Lambda}^{0} candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 150, 1.05, 1.2);
	fOutputList->Add(fHistMassAntiLambda);
	
	fHistPtVsMassK0                      = new TH2F("h2PtVsMassK0","K^{0} candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});p{t} (GeV/c)",400, 0.4, 0.6,240,0,12);
	fOutputList->Add(fHistPtVsMassK0);
	
	fHistPtVsMassLambda                  = new TH2F("h2PtVsMassLambda","#Lambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p{t} (GeV/c)",280, 1.06, 1.2,240,0,12);
	fOutputList->Add(fHistPtVsMassLambda);
	
	fHistPtVsMassAntiLambda              = new TH2F("h2PtVsMassAntiLambda","#AntiLambda^{0} candidates;M(p#pi^{-}) (GeV/c^{2});p{t} (GeV/c)",280, 1.06, 1.2,240,0,12);
	fOutputList->Add(fHistPtVsMassAntiLambda);
	
	fHistArmenterosPodolanskiK0          = new TH2F("h2ArmenterosPodolanskiK0","Armenteros-Podolanski phase space;#alpha;p{t} arm",100,-1.0,1.0,50,0,0.5);
	fOutputList->Add(fHistArmenterosPodolanskiK0);
	
	fHistArmenterosPodolanskiLambda      = new TH2F("h2ArmenterosPodolanskiLambda","Armenteros-Podolanski phase space;#alpha;p{t} arm",100,-1.0,1.0,50,0,0.5);
	fOutputList->Add(fHistArmenterosPodolanskiLambda);
	
	fHistArmenterosPodolanskiAntiLambda  = new TH2F("h2ArmenterosPodolanskiAntiLambda","Armenteros-Podolanski phase space;#alpha;p{t} arm",100,-1.0,1.0,50,0,0.5);
	fOutputList->Add(fHistArmenterosPodolanskiAntiLambda);
	
//--------------------------------------------MC Associated histograms -----------------------------------------------------//
	
	fHistAsMcPtK0                        = new TH1F("h1AsMcPtK0", "K^{0} associated;p{t} (GeV/c);Counts", 240,0,12);
	fOutputList->Add(fHistAsMcPtK0);
	
	fHistAsMcPtLambda                    = new TH1F("h1AsMcPtLambda", "#Lambda^{0} associated;p{t} (GeV/c);Counts", 240,0,12);
	fOutputList->Add(fHistAsMcPtLambda);
	
	fHistAsMcPtAntiLambda                = new TH1F("h1AsMcPtAntiLambda", "#AntiLambda^{0} associated;p{t} (GeV/c);Counts", 240,0,12);
	fOutputList->Add(fHistAsMcPtAntiLambda);
	
	fHistAsMcProdRadiusK0                = new TH1F("h1AsMcProdRadiusK0", "K^{0} associated;r (cm);Counts", 500, 0, 100);
	fOutputList->Add(fHistAsMcProdRadiusK0);
	
	fHistAsMcProdRadiusLambda            = new TH1F("h1AsMcProdRadiusLambda", "#Lambda^{0} associated;r (cm);Counts", 500, 0, 100);
	fOutputList->Add(fHistAsMcProdRadiusLambda);
	
	fHistAsMcProdRadiusAntiLambda        = new TH1F("h1AsMcProdRadiusAntiLambda", "#bar{#Lambda}^{0} associated;r (cm);Counts", 500, 0, 100);
	fOutputList->Add(fHistAsMcProdRadiusAntiLambda);
	
	fHistAsMcProdRadiusXvsYK0s           = new TH2F("h2AsMcProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
	fOutputList->Add(fHistAsMcProdRadiusXvsYK0s);
	
	fHistAsMcProdRadiusXvsYLambda        = new TH2F("h2AsMcProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
	fOutputList->Add(fHistAsMcProdRadiusXvsYLambda);
	
	fHistAsMcProdRadiusXvsYAntiLambda    = new TH2F("h2AsMcProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-50,50,200,-50,50);
	fOutputList->Add(fHistAsMcProdRadiusXvsYAntiLambda);
	
	fHistPidMcMassK0                     = new TH1F("h1PidMcMassK0", "K^{0} MC PId checked;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
	fOutputList->Add(fHistPidMcMassK0);
	
	fHistPidMcMassLambda                 = new TH1F("h1PidMcMassLambda", "#Lambda^{0} MC PId checked;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
	fOutputList->Add(fHistPidMcMassLambda);
	
	fHistPidMcMassAntiLambda             = new TH1F("h1PidMcMassAntiLambda", "#bar{#Lambda}^{0} MC PId checked;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
	fOutputList->Add(fHistPidMcMassAntiLambda);
	
	fHistAsMcPtLambdaFromSigma           = new TH1F("h1AsMcPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p{t} (GeV/c);Count",240,0,12);
	fOutputList->Add(fHistAsMcPtLambdaFromSigma);
	
	fHistAsMcPtAntiLambdaFromSigma       = new TH1F("h1AsMcPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p{t} (GeV/c);Count",240,0,12);
	fOutputList->Add(fHistAsMcPtAntiLambdaFromSigma);
	
	fHistAsMcSecondaryPtVsRapK0s         = new TH2F("h2AsMcSecondaryPtVsRapK0s", "K^{0} associated secondary;p{t} (GeV/c);rapidity",240,0,12,30,-1.5,1.5);
	fOutputList->Add(fHistAsMcSecondaryPtVsRapK0s);
	
	fHistAsMcSecondaryPtVsRapLambda      = new TH2F("h2AsMcSecondaryPtVsRapLambda", "#Lambda^{0} associated secondary;p{t} (GeV/c);rapidity",240,0,12,30,-1.5,1.5);
	fOutputList->Add(fHistAsMcSecondaryPtVsRapLambda);
	
	fHistAsMcSecondaryPtVsRapAntiLambda  = new TH2F("h2AsMcSecondaryPtVsRapAntiLambda", "#bar{#Lambda}^{0} associated secondary;p{t} (GeV/c);rapidity",240,0,12,30,-1.5,1.5);
	fOutputList->Add(fHistAsMcSecondaryPtVsRapAntiLambda);

	fHistAsMcSecondaryProdRadiusK0s      = new TH1F("h1AsMcSecondaryProdRadiusK0s", "K^{0} Production Radius;r (cm);Count", 170, -2, 15);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusK0s);
	
	fHistAsMcSecondaryProdRadiusLambda   = new TH1F("h1AsMcSecondaryProdRadiusLambda", "#Lambda^{0} Production Radius;r (cm);Count", 170, -2, 15);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusLambda);
	
	fHistAsMcSecondaryProdRadiusAntiLambda       = new TH1F("h1AsMcSecondaryProdRadiusAntiLambda", "#bar{#Lambda}^{0} Production Radius;r (cm);Count", 170, -2, 15);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusAntiLambda);  
	
	fHistAsMcSecondaryProdRadiusXvsYK0s          = new TH2F("h2AsMcSecondaryProdRadiusXvsYK0s","Associated Secondary K^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusXvsYK0s);
	
	fHistAsMcSecondaryProdRadiusXvsYLambda       = new TH2F("h2AsMcSecondaryProdRadiusXvsYLambda","Associated Secondary #Lambda^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusXvsYLambda);
	
	fHistAsMcSecondaryProdRadiusXvsYAntiLambda   = new TH2F("h2AsMcSecondaryProdRadiusXvsYAntiLambda","Associated Secondary #bar{#Lambda}^{0} Production Radius;x (cm); y (cm)",200,-20,20,200,-20,20);
	fOutputList->Add(fHistAsMcSecondaryProdRadiusXvsYAntiLambda);

	fHistAsMcSecondaryPtLambdaFromSigma          = new TH1F("h1AsMcSecondaryPtLambdaFromSigma","#Lambda}^{0} associated from Sigma;p{t} (GeV/c);Count",240,0,12);
	fOutputList->Add(fHistAsMcSecondaryPtLambdaFromSigma);
	
	fHistAsMcSecondaryPtAntiLambdaFromSigma      = new TH1F("h1AsMcSecondaryPtAntiLambdaFromSigma","#bar{#Lambda}^{0} associated from Sigma;p{t} (GeV/c);Count",240,0,12);
	fOutputList->Add(fHistAsMcSecondaryPtAntiLambdaFromSigma);
	
	//----------------------------------------------Correlation histograms -----------------------------------------------------//
	
	fHistSibK0       = new TH3F("hfHistSibK0","",CorrBinsX,-TMath::Pi()/2,3*TMath::Pi()/2,CorrBinsY,-2*fTrackEtaCut,2*fTrackEtaCut,30,0,6);
	fHistSibK0->SetStats(0);
	fHistSibK0->Sumw2();
	fOutputList->Add(fHistSibK0);
	
	fHistMixK0       = new TH3F("hfHistMixK0","",CorrBinsX,-TMath::Pi()/2,3*TMath::Pi()/2,CorrBinsY,-2*fTrackEtaCut,2*fTrackEtaCut,30,0,6);
	fHistMixK0->SetStats(0);
	fHistMixK0->Sumw2();
	fOutputList->Add(fHistMixK0);
	
	fHistSibLambda   = new TH3F("hfHistSibLambda","",CorrBinsX,-TMath::Pi()/2,3*TMath::Pi()/2,CorrBinsY,-2*fTrackEtaCut,2*fTrackEtaCut,30,0,6);
	fHistSibLambda->SetStats(0);
	fHistSibLambda->Sumw2();
	fOutputList->Add(fHistSibLambda);
	
	fHistMixLambda   = new TH3F("hfHistMixLambda","",CorrBinsX,-TMath::Pi()/2,3*TMath::Pi()/2,CorrBinsY,-2*fTrackEtaCut,2*fTrackEtaCut,30,0,6);
	fHistMixLambda->SetXTitle("#Delta #Phi");
	fHistMixLambda->SetStats(0);
	fHistMixLambda->Sumw2();
	fOutputList->Add(fHistMixLambda);
	
	fHistSibK0MC     = new TH3F("hfHistSibK0MC","",CorrBinsX,-TMath::Pi()/2,3*TMath::Pi()/2,CorrBinsY,-2*fTrackEtaCut,2*fTrackEtaCut,30,0,6);
	fHistSibK0MC->SetStats(0);
	fHistSibK0MC->Sumw2();
	fOutputList->Add(fHistSibK0MC);
	
	fHistMixK0MC     = new TH3F("hfHistMixK0MC","",CorrBinsX,-TMath::Pi()/2,3*TMath::Pi()/2,CorrBinsY,-2*fTrackEtaCut,2*fTrackEtaCut,30,0,6);
	fHistMixK0MC->SetStats(0);
	fHistMixK0MC->Sumw2();
	fOutputList->Add(fHistMixK0MC);
	
	fHistSibLambdaMC = new TH3F("hfHistSibLambdaMC","",CorrBinsX,-TMath::Pi()/2,3*TMath::Pi()/2,CorrBinsY,-2*fTrackEtaCut,2*fTrackEtaCut,30,0,6);
	fHistSibLambdaMC->SetStats(0);
	fHistSibLambdaMC->Sumw2();
	fOutputList->Add(fHistSibLambdaMC);
	
	fHistMixLambdaMC = new TH3F("hfHistMixLambdaMC","",CorrBinsX,-TMath::Pi()/2,3*TMath::Pi()/2,CorrBinsY,-2*fTrackEtaCut,2*fTrackEtaCut,30,0,6);
	fHistMixLambdaMC->SetStats(0);
	fHistMixLambdaMC->Sumw2();
	fOutputList->Add(fHistMixLambdaMC);
	
	fHistLeadInfo = new TH3F("hfHistLeadInfo","",60,0,12,CorrBinsY/2,-fTrackEtaCut,fTrackEtaCut,CorrBinsX/2,0,3*TMath::Pi());
	fOutputList->Add(fHistLeadInfo);
	
	fHistLeadInfoMC = new TH3F("hfHistLeadInfoMC","",60,0,12,CorrBinsY/2,-fTrackEtaCut,fTrackEtaCut,CorrBinsX/2,0,3*TMath::Pi());
	fOutputList->Add(fHistLeadInfoMC);
	
	fHistLeadInfoMix = new TH3F("hfHistLeadInfoMix","",60,0,12,CorrBinsY/2,-fTrackEtaCut,fTrackEtaCut,CorrBinsX/2,0,3*TMath::Pi());
	fOutputList->Add(fHistLeadInfoMix);
	
	fHistLeadInfoMixMC = new TH3F("hfHistLeadInfoMixMC","",60,0,12,CorrBinsY/2,-fTrackEtaCut,fTrackEtaCut,CorrBinsX/2,0,3*TMath::Pi());
	fOutputList->Add(fHistLeadInfoMixMC);
	
	//----------------------------------------------Event Pool-----------------------------------------------------//
	
	fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins, fCentBins, fNzVtxBins, fZvtxBins);
	if(!fPoolMgr) return;
	
	if(fAnalysisMC){
	fPoolMgrMC = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins, fCentBins, fNzVtxBins, fZvtxBins);
	if(!fPoolMgr) return;
	}
	
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
	if(!fPIDResponse) return;
	
    fAODEvent = dynamic_cast<AliAODEvent*>(inEvMain->GetEvent());
	if(!fAODEvent) return;
	
	// physics selection
	UInt_t maskIsSelected = inEvMain->IsEventSelected();
	Bool_t isSelected = 0;
	
	if( fTriggerMask == "kMB" )
	    isSelected = (maskIsSelected & AliVEvent::kMB)     == AliVEvent::kMB;
	if( fTriggerMask == "kINT7" )
		isSelected = (maskIsSelected & AliVEvent::kINT7)   == AliVEvent::kINT7;
	if( fTriggerMask == "kINT8" )
		isSelected = (maskIsSelected & AliVEvent::kINT8)   == AliVEvent::kINT8;
	if( fTriggerMask == "kAnyINT" )
		isSelected = (maskIsSelected & AliVEvent::kAnyINT) == AliVEvent::kAnyINT;
	if ( ! isSelected )return;
	
	//-----------------------------MC Accsess------------------------------------------------
	TClonesArray *stack = 0x0;
	Double_t mcXv=0., mcYv=0., mcZv=0.;
	Int_t ntrk =0, ntrk0=0;
	
	TObjArray *selectedTracksLeadingMC =0x0;
	
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
		
		selectedTracksLeadingMC = FindLeadingObjectsMC(stack);
		if(!selectedTracksLeadingMC) return;
		selectedTracksLeadingMC->SetOwner(kTRUE);
	}
	
	// If PID is used:
	Double_t lLimitPPID          = 1.0;
	Float_t  cutNSigmaLowP       = 3.0;
	Float_t  cutNSigmaHighP      = 3.0;	
	Int_t    lPdgcodeCurrentPart = 0;
	Double_t lRapCurrentPart     = 0;
	Double_t lPtCurrentPart		 = 0;
	Double_t lPhiCurrentPart	 = 0;
	Double_t lEtaCurrentPart	 = 0;
	Int_t    lComeFromSigma      = 0;
	
	// PID flags:
	Int_t LambdaPID = 0;
	Int_t AntiLambdaPID = 0;
	
	// current mc particle 's mother
	Int_t iCurrentMother  = 0, lPdgCurrentMother    = 0;

	// Start loop over MC particles
	if (fAnalysisMC) {
		
		fHistMCPrimaryVertexX->Fill(mcXv);
		fHistMCPrimaryVertexY->Fill(mcYv);
		fHistMCPrimaryVertexZ->Fill(mcZv);
		
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
			
			if (( ( TMath::Abs(lPdgCurrentMother)  == 3212)  ||
				  ( TMath::Abs(lPdgCurrentMother)  == 3224)  ||
				  ( TMath::Abs(lPdgCurrentMother)  == 3214)  ||
				  ( TMath::Abs(lPdgCurrentMother)  == 3114) )
			   && ( Mother->GetMother() == -1)
				) lComeFromSigma = 1;
			else lComeFromSigma  = 0;
			
			Bool_t isPrimary=p0->IsPhysicalPrimary();
			if (!isPrimary)continue; //keep only primary particles
			
			// Rapidity Cut
			if (TMath::Abs(lRapCurrentPart) > fRapidityCut) continue;
			
			if (lPdgcodeCurrentPart==310) {
				fHistMCRapK0s->Fill(lRapCurrentPart);
				fHistMCPtK0s->Fill(lPtCurrentPart);
				fHistPrimRawPtVsYK0s->Fill(lPtCurrentPart,lRapCurrentPart);
				fHistMCPtAllK0s->Fill(lPtCurrentPart);
				selectedK0MC->Add(new AliLeadingBasicParticle(lEtaCurrentPart,lPhiCurrentPart,lPtCurrentPart));
			}
			else 
			if (lPdgcodeCurrentPart==3122) {
				    fHistMCRapLambda->Fill(lRapCurrentPart);
					fHistMCPtLambda->Fill(lPtCurrentPart);	  
				    fHistPrimRawPtVsYLambda->Fill(lPtCurrentPart,lRapCurrentPart);
					fHistMCPtAllLambda->Fill(lPtCurrentPart);
					if (lComeFromSigma) fHistMCPtLambdaFromSigma->Fill(lPtCurrentPart);
				    selectedLambdaMC->Add(new AliLeadingBasicParticle(lEtaCurrentPart,lPhiCurrentPart,lPtCurrentPart));
					
				}
			else 
			if (lPdgcodeCurrentPart==-3122) {
				        fHistMCRapAntiLambda->Fill(lRapCurrentPart);
						fHistMCPtAntiLambda->Fill(lPtCurrentPart);	  
						fHistPrimRawPtVsYAntiLambda->Fill(lPtCurrentPart,lRapCurrentPart);
						fHistMCPtAllAntiLambda->Fill(lPtCurrentPart);
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
		if ((MultipOrCent < 0.)||(MultipOrCent > 90.)) return;
	}
	
	Double_t * CentBins = fCentBins;
	Double_t poolmin    = CentBins[0];
	Double_t poolmax    = CentBins[fNCentBins];
	
	TObjArray *selectedTracksLeading = FindLeadingObjects(fAODEvent);
	if(!selectedTracksLeading) return;
	selectedTracksLeading->SetOwner(kTRUE);
	
	// -------------------------------------V0 loop for reconstructed event------------------------
	Double_t lPLambda     = 0;
	Double_t lPAntiLambda = 0;
	Double_t lPK0s        = 0;
	
	// Variables:
	Double_t  lV0Position[3];
	
	Double_t lDcaPosToPrimVertex = 0;
	Double_t lDcaNegToPrimVertex = 0;
	Double_t lDcaV0Daughters     = 0;
	Double_t lV0cosPointAngle    = 0;
	Double_t lV0DecayLength      = 0;
	Double_t lV0Radius           = 0;
	Double_t lDcaV0ToPrimVertex  = 0;
	Double_t lcTauLambda         = 0;   
	Double_t lcTauAntiLambda     = 0;   
	Double_t lcTauK0s            = 0;   
	
	Double_t lInvMassK0   = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
	Double_t lPtK0s       = 0, lPtLambda      = 0, lPtAntiLambda      = 0;
	Double_t lPhiK0s      = 0, lPhiLambda     = 0, lPhiAntiLambda     = 0;
	Double_t lEtaK0s      = 0, lEtaLambda     = 0, lEtaAntiLambda     = 0;
	Double_t lRapK0s      = 0, lRapLambda     = 0, lRapAntiLambda     = 0;
	Double_t lPzK0s       = 0, lPzLambda      = 0, lPzAntiLambda      = 0;
	Double_t lAlphaV0     = 0, lPtArmV0       = 0;
	
	//Associated V0s:
	
	UInt_t   lLabelTrackPos       = 0, lLabelTrackNeg         = 0;
	Int_t    lCheckPIdK0Short     = 0, lCheckMcK0Short        = 0;
	Int_t    lCheckPIdLambda      = 0, lCheckMcLambda         = 0;
	Int_t    lCheckPIdAntiLambda  = 0, lCheckMcAntiLambda     = 0;
	Int_t    lCheckSecondaryK0s   = 0, lCheckSecondaryLambda  = 0, lCheckSecondaryAntiLambda  = 0;
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
	
	// Daughters' momentum:
	Double_t  lMomPos[3] = {999,999,999};
	Double_t  lMomNeg[3] = {999,999,999};
	Double_t  lPtPos = 999, lPtNeg = 999;
	
	// Inner Wall parameters:
	Double_t  lMomInnerWallPos =999, lMomInnerWallNeg = 999;
	Double_t  ldEdxPos         =0.0, ldEdxNeg         = 0.0;
	
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
		lComeFromSigma      = -1;
		
        // get daughters
		
   	    AliAODTrack *myTrackPos=(AliAODTrack *)(aodV0->GetDaughter(0));
        AliAODTrack *myTrackNeg=(AliAODTrack *)(aodV0->GetDaughter(1));
		
		if (!myTrackPos || !myTrackNeg) {Printf("ERROR: Could not retreive one of the daughter track");continue;}
        if (!IsAcseptedV0(fAODEvent,aodV0,myTrackPos,myTrackNeg)) continue;
		
		// VO's main characteristics to check the reconstruction cuts
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
		lPtAntiLambda      = aodV0->Pt();
		lPhiAntiLambda     = aodV0->Phi();
		lEtaAntiLambda     = aodV0->Eta();
		lPzAntiLambda      = aodV0->Pz();
		
		// Rapidity:
		lRapK0s    = aodV0->RapK0Short();
		lRapLambda = aodV0->RapLambda();
		lRapAntiLambda = aodV0->Y(-3122);
		
		if (lPtK0s==0) {continue;}
		if (lPtLambda==0) {continue;}
		if (lPtAntiLambda==0) {continue;}
		
		if (fUsePID=="withPID") {
			nSigmaPosPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion));
			nSigmaNegPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kPion));
			nSigmaPosProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton));
			nSigmaNegProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kProton));
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
				mcPosMotherX = lMCAODMother->Xv();
				mcPosMotherY = lMCAODMother->Yv();
				mcPosMotherZ = lMCAODMother->Zv();
				mcPosMotherR = TMath::Sqrt(mcPosMotherX*mcPosMotherX+mcPosMotherY*mcPosMotherY);
				
				mcMotherPt   = lMCAODMother->Pt();
			}
		}
		
	}
	
	if (fAnalysisMC) {		
		if(((lPDGCodePosDaughter==+211) && (lPDGCodeNegDaughter==-211))){lCheckPIdK0Short = 1;
			if ( (lIndexPosMother==lIndexNegMother) &&(lPdgcodeMother==310) ){
				if (((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary()) lCheckMcK0Short  = 1;
				else lCheckSecondaryK0s = 1;
			}
		}
		else if(((lPDGCodePosDaughter==+2212) && (lPDGCodeNegDaughter==-211))){lCheckPIdLambda = 1;
			if ((lIndexPosMother==lIndexNegMother) &&(lPdgcodeMother==3122)  ){
				if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3224) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3214) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3114))lComeFromSigma = 1;
				else lComeFromSigma = 0; 
				if (((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary()||
				   (!(((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary() )&&
				   (lComeFromSigma)))lCheckMcLambda  = 1; 
				else lCheckSecondaryLambda    = 1;
			}
		}
		else if(((lPDGCodePosDaughter==211) && (lPDGCodeNegDaughter==-2212))){lCheckPIdAntiLambda = 1;
			if ( (lIndexPosMother==lIndexNegMother) &&
				(lPdgcodeMother==-3122) ) {
				if ( ( TMath::Abs(lPdgcodeMotherOfMother) == 3212) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3224) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3214) ||
					( TMath::Abs(lPdgcodeMotherOfMother)  == 3114)
					)lComeFromSigma = 1;
				else lComeFromSigma = 0;  
				if (((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary() || 
				   ((!((AliAODMCParticle*)stack->UncheckedAt(lIndexPosMother))->IsPrimary())&&
				   (lComeFromSigma)))lCheckMcAntiLambda  = 1;
				else lCheckSecondaryAntiLambda = 1;
			}
		}
	} // end "look for associated particles 
	
	// PID condition:
	lCheckPIDK0sPosDaughter        = 0, lCheckPIDK0sNegDaughter        = 0;
	lCheckPIDLambdaPosDaughter     = 0, lCheckPIDLambdaNegDaughter     = 0;
	lCheckPIDAntiLambdaPosDaughter = 0, lCheckPIDAntiLambdaNegDaughter = 0;
	
	if (lMomInnerWallPos < lLimitPPID) {
		if (nSigmaPosPion < cutNSigmaLowP){
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
		if (nSigmaNegPion < cutNSigmaLowP){
			lCheckPIDK0sNegDaughter       = 1;
			lCheckPIDLambdaNegDaughter    = 1;
		}
		if (nSigmaNegProton < cutNSigmaLowP)  lCheckPIDAntiLambdaNegDaughter = 1;
		
	}
	else if (lMomInnerWallNeg > lLimitPPID) {
		if (nSigmaNegPion < cutNSigmaHighP){
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
	
	if(lPLambda > 0)     lcTauLambda     = (lV0DecayLength*lInvMassLambda)/lPLambda;
	if(lPAntiLambda > 0) lcTauAntiLambda = (lV0DecayLength*lInvMassAntiLambda)/lPAntiLambda; 
	if(lPK0s > 0)        lcTauK0s        = (lV0DecayLength*lInvMassK0)/lPK0s;	
	//--------------------------------------------K0s---------------------------------//
    if(!fAnalysisMC){
	if (lcTauK0s< fCTauK0){
		if (TMath::Abs(lRapK0s) < fRapidityCut ){
			if(lPtArmV0*fSpecialArmenterosCutK0s>(TMath::Abs(lAlphaV0))){
				fHistMassK0->Fill(lInvMassK0);
				fHistPtVsMassK0->Fill(lInvMassK0,lPtK0s);
				fHistDcaPosToPrimVertexK0vsMassK0->Fill(lDcaPosToPrimVertex,lInvMassK0);
				fHistDcaNegToPrimVertexK0vsMassK0->Fill(lDcaNegToPrimVertex,lInvMassK0);
				fHistRadiusV0K0vsMassK0->Fill(lV0Radius,lInvMassK0);
				fHistDecayLengthV0K0vsMassK0->Fill(lV0DecayLength,lInvMassK0);
				fHistDcaV0DaughtersK0vsMassK0->Fill(lDcaV0Daughters,lInvMassK0);
				fHistCosPointAngleK0vsMassK0->Fill(lV0cosPointAngle,lInvMassK0);
				fHistArmenterosPodolanskiK0->Fill(lAlphaV0,lPtArmV0);
			if(IsK0InvMass(lInvMassK0)){selectedK0->Add(new AliLeadingBasicParticle(lEtaK0s,lPhiK0s,lPtK0s));}
		  }//Special Armesto Cut for K0	
		} // if rap. condition
	} // end cTau condition
	//-----------------------------------------Lambda---------------------------------//
	if (lcTauLambda < fCTauLambda){
		if ((LambdaPID==1 && lMomInnerWallPos  <=1 ) || (lMomInnerWallPos >1 ) ||  !(fUsePID=="withPID")){ 	
			if (TMath::Abs(lRapLambda) <fRapidityCut) {
					fHistMassLambda->Fill(lInvMassLambda);
					fHistPtVsMassLambda->Fill(lInvMassLambda,lPtLambda);
					fHistDcaPosToPrimVertexLvsMassL->Fill(lDcaPosToPrimVertex,lInvMassLambda);
					fHistDcaNegToPrimVertexLvsMassL->Fill(lDcaNegToPrimVertex,lInvMassLambda);
					fHistRadiusV0LvsMassL->Fill(lV0Radius,lInvMassLambda);
					fHistDecayLengthV0LvsMassL->Fill(lV0DecayLength,lInvMassLambda);
					fHistDcaV0DaughtersLvsMassL->Fill(lDcaV0Daughters,lInvMassLambda);
					fHistCosPointAngleLvsMassL->Fill(lV0cosPointAngle,lInvMassLambda);
					fHistArmenterosPodolanskiLambda->Fill(lAlphaV0,lPtArmV0);
				if(IsLambdaInvMass(lInvMassLambda)){selectedLambda->Add(new AliLeadingBasicParticle(lEtaLambda,lPhiLambda,lPtLambda));}
			} //end of Rap condition
		}//End PID
	}  //end cTau condition
 
	//--------------------------------------AntiLambda---------------------------------//
	if (lcTauAntiLambda < fCTauLambda){
		if ((AntiLambdaPID==1 && lMomInnerWallNeg <=1) || (lMomInnerWallNeg >1) ||  !(fUsePID=="withPID")){  
			if (TMath::Abs(lRapAntiLambda) < fRapidityCut) {
					fHistMassAntiLambda->Fill(lInvMassAntiLambda);
					fHistPtVsMassAntiLambda->Fill(lInvMassAntiLambda,lPtAntiLambda);
					fHistDcaPosToPrimVertexAntiLvsMass->Fill(lDcaPosToPrimVertex,lInvMassAntiLambda);
					fHistDcaNegToPrimVertexAntiLvsMass->Fill(lDcaNegToPrimVertex,lInvMassAntiLambda);
					fHistRadiusV0AntiLvsMass->Fill(lV0Radius,lInvMassAntiLambda);
					fHistDecayLengthV0AntiLvsMass->Fill(lV0DecayLength,lInvMassAntiLambda);
					fHistDcaV0DaughtersAntiLvsMass->Fill(lDcaV0Daughters,lInvMassAntiLambda);
					fHistCosPointAngleAntiLvsMass->Fill(lV0cosPointAngle,lInvMassAntiLambda);
					fHistArmenterosPodolanskiAntiLambda->Fill(lAlphaV0,lPtArmV0);
			} //end of Rap condition
		} // end of PID condition
	} //end cTau condition
}
	if(fAnalysisMC) {
	//--------------------------------------K0s Associated---------------------------------//

	if (lcTauK0s< fCTauK0) {
		if (TMath::Abs(lRapK0s) < fRapidityCut) {
					if(lCheckPIdK0Short) fHistPidMcMassK0->Fill(lInvMassK0);
					if(lCheckMcK0Short) {
						fHistAsMcPtK0->Fill(lPtK0s);
						fHistAsMcProdRadiusK0->Fill(mcPosMotherR);
						fHistAsMcProdRadiusXvsYK0s->Fill(mcPosMotherX,mcPosMotherY);
					if(lPtArmV0*fSpecialArmenterosCutK0s>(TMath::Abs(lAlphaV0))){
						fHistMassK0->Fill(lInvMassK0);
						fHistPtVsMassK0->Fill(lInvMassK0,lPtK0s);
						fHistDcaPosToPrimVertexK0vsMassK0->Fill(lDcaPosToPrimVertex,lInvMassK0);
						fHistDcaNegToPrimVertexK0vsMassK0->Fill(lDcaNegToPrimVertex,lInvMassK0);
						fHistRadiusV0K0vsMassK0->Fill(lV0Radius,lInvMassK0);
						fHistDecayLengthV0K0vsMassK0->Fill(lV0DecayLength,lInvMassK0);
						fHistDcaV0DaughtersK0vsMassK0->Fill(lDcaV0Daughters,lInvMassK0);
						fHistCosPointAngleK0vsMassK0->Fill(lV0cosPointAngle,lInvMassK0);
						fHistArmenterosPodolanskiK0->Fill(lAlphaV0,lPtArmV0);
					if(IsK0InvMass(lInvMassK0)){selectedK0->Add(new AliLeadingBasicParticle(lEtaK0s,lPhiK0s,lPtK0s));}
							}
					}
					if (lCheckSecondaryK0s) {
						fHistAsMcSecondaryPtVsRapK0s->Fill(lPtK0s,lRapK0s);
						fHistAsMcSecondaryProdRadiusK0s->Fill(mcPosMotherR);
						fHistAsMcSecondaryProdRadiusXvsYK0s->Fill(mcPosMotherX,mcPosMotherY);
					}
		} // end rapidity condition
	} //end cTau condition
	//-----------------------------------Lambda Associated---------------------------------//
	if (lcTauLambda < fCTauLambda){                                                                                                   
		if (TMath::Abs(lRapLambda) < fRapidityCut) {
					if(lCheckPIdLambda) fHistPidMcMassLambda->Fill(lInvMassLambda);
					if(lCheckMcLambda) {
						fHistAsMcPtLambda->Fill(lPtLambda);
						fHistAsMcProdRadiusLambda->Fill(mcPosMotherR);
						fHistAsMcProdRadiusXvsYLambda->Fill(mcPosMotherX,mcPosMotherY);
					if (lComeFromSigma) fHistAsMcPtLambdaFromSigma->Fill(lPtLambda);
						fHistMassLambda->Fill(lInvMassLambda);
						fHistPtVsMassLambda->Fill(lInvMassLambda,lPtLambda);
						fHistDcaPosToPrimVertexLvsMassL->Fill(lDcaPosToPrimVertex,lInvMassLambda);
						fHistDcaNegToPrimVertexLvsMassL->Fill(lDcaNegToPrimVertex,lInvMassLambda);
						fHistRadiusV0LvsMassL->Fill(lV0Radius,lInvMassLambda);
						fHistDecayLengthV0LvsMassL->Fill(lV0DecayLength,lInvMassLambda);
						fHistDcaV0DaughtersLvsMassL->Fill(lDcaV0Daughters,lInvMassLambda);
						fHistCosPointAngleLvsMassL->Fill(lV0cosPointAngle,lInvMassLambda);
						fHistArmenterosPodolanskiLambda->Fill(lAlphaV0,lPtArmV0);
					if(IsLambdaInvMass(lInvMassLambda)){selectedLambda->Add(new AliLeadingBasicParticle(lEtaLambda,lPhiLambda,lPtLambda));}
					}
					if (lCheckSecondaryLambda) {
						fHistAsMcSecondaryPtVsRapLambda->Fill(lPtLambda,lRapLambda);
						fHistAsMcSecondaryProdRadiusLambda->Fill(mcPosMotherR); 
						fHistAsMcSecondaryProdRadiusXvsYLambda->Fill(mcPosMotherX,mcPosMotherY);
					if (lComeFromSigma) fHistAsMcSecondaryPtLambdaFromSigma->Fill(lPtLambda);
					}
		} // end rapidity condition
	}//end cTau condition
//------------------------------AntiLambda Associated---------------------------------//	
	if (lcTauAntiLambda < fCTauLambda){
		if (TMath::Abs(lRapAntiLambda) < fRapidityCut) {
					if(lCheckPIdAntiLambda) fHistPidMcMassAntiLambda->Fill(lInvMassAntiLambda);
					if(lCheckMcAntiLambda) {
						fHistAsMcPtAntiLambda->Fill(lPtAntiLambda);
						fHistAsMcProdRadiusAntiLambda->Fill(mcPosMotherR);
						fHistAsMcProdRadiusXvsYAntiLambda->Fill(mcPosMotherX,mcPosMotherY);
					if (lComeFromSigma) fHistAsMcPtAntiLambdaFromSigma->Fill(lPtAntiLambda);
							fHistMassAntiLambda->Fill(lInvMassAntiLambda);
							fHistPtVsMassAntiLambda->Fill(lInvMassAntiLambda,lPtAntiLambda);
							fHistDcaPosToPrimVertexAntiLvsMass->Fill(lDcaPosToPrimVertex,lInvMassAntiLambda);
							fHistDcaNegToPrimVertexAntiLvsMass->Fill(lDcaNegToPrimVertex,lInvMassAntiLambda);
							fHistRadiusV0AntiLvsMass->Fill(lV0Radius,lInvMassAntiLambda);
							fHistDecayLengthV0AntiLvsMass->Fill(lV0DecayLength,lInvMassAntiLambda);
							fHistDcaV0DaughtersAntiLvsMass->Fill(lDcaV0Daughters,lInvMassAntiLambda);
							fHistCosPointAngleAntiLvsMass->Fill(lV0cosPointAngle,lInvMassAntiLambda);
							fHistArmenterosPodolanskiAntiLambda->Fill(lAlphaV0,lPtArmV0);
					}
					if (lCheckSecondaryAntiLambda) {
						fHistAsMcSecondaryPtVsRapAntiLambda->Fill(lPtAntiLambda,lRapAntiLambda);
						fHistAsMcSecondaryProdRadiusAntiLambda->Fill(mcPosMotherR); 
						fHistAsMcSecondaryProdRadiusXvsYAntiLambda->Fill(mcPosMotherX,mcPosMotherY);
					if (lComeFromSigma) fHistAsMcSecondaryPtAntiLambdaFromSigma->Fill(lPtAntiLambda);
					}
		} // end rapidity condition      
	}//end cTau condition
	
	// Correlation part
	FillCorrelations(selectedTracksLeadingMC,selectedK0MC,fHistSibK0MC,fHistLeadInfoMC);
	FillCorrelations(selectedTracksLeadingMC,selectedLambdaMC,fHistSibLambdaMC,0);
		
	// Mixing part
	if(TMath::Abs(lPVz)>=10 || MultipOrCent>poolmax || MultipOrCent < poolmin) {
		if(fcollidingSys=="PP")AliInfo(Form("pp Event with Zvertex = %.2f cm and multiplicity = %.0f out of pool bounds, SKIPPING",mcZv,MultipOrCent));
		if(fcollidingSys=="PbPb") AliInfo(Form("PbPb Event with Zvertex = %.2f cm and centrality = %.1f  out of pool bounds, SKIPPING",mcZv,MultipOrCent));
		return;
	}
		
	fPoolK0MC = fPoolMgrMC->GetEventPool(MultipOrCent, mcZv);
	if (!fPoolK0MC){AliInfo(Form("No pool found for multiplicity = %f, zVtx = %f cm", MultipOrCent, mcZv));return;}
		
	if (fPoolK0MC->IsReady() || fPoolK0MC->NTracksInPool() > fPoolMinNTracks  || fPoolK0MC->GetCurrentNEvents() >= fMinEventsToMix)
	{
		for (Int_t jMix=0; jMix<fPoolK0MC->GetCurrentNEvents(); jMix++) 
		FillCorrelations(selectedTracksLeadingMC, fPoolK0MC->GetEvent(jMix),fHistMixK0MC,fHistLeadInfoMixMC);
		fPoolK0MC->UpdatePool(CloneAndReduceTrackList(selectedK0MC));
	}
		
	fPoolLambdaMC = fPoolMgrMC->GetEventPool(MultipOrCent, mcZv);
	if (!fPoolLambdaMC){AliInfo(Form("No pool found for multiplicity = %f, zVtx = %f cm", MultipOrCent, mcZv));return;}
		
	if (fPoolLambdaMC->IsReady() || fPoolLambdaMC->NTracksInPool() > fPoolMinNTracks  || fPoolLambdaMC->GetCurrentNEvents() >= fMinEventsToMix)
	{
		for (Int_t jMix=0; jMix<fPoolLambdaMC->GetCurrentNEvents(); jMix++) 
		FillCorrelations(selectedTracksLeadingMC, fPoolLambdaMC->GetEvent(jMix),fHistMixLambdaMC,0);
		fPoolLambdaMC->UpdatePool(CloneAndReduceTrackList(selectedLambdaMC));
	}
}
	
	FillCorrelations(selectedTracksLeading,selectedK0,fHistSibK0,fHistLeadInfo);
	FillCorrelations(selectedTracksLeading,selectedLambda,fHistSibLambda,0);
	
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
			FillCorrelations(selectedTracksLeading, fPoolK0->GetEvent(jMix),fHistMixK0,fHistLeadInfoMix);
		fPoolK0->UpdatePool(CloneAndReduceTrackList(selectedK0));
	}
	
	fPoolLambda = fPoolMgr->GetEventPool(MultipOrCent, lPVz);
	if (!fPoolLambda){AliInfo(Form("No pool found for multiplicity = %f, zVtx = %f cm", MultipOrCent, lPVz));return;}
	
    if (fPoolLambda->IsReady() || fPoolLambda->NTracksInPool() > fPoolMinNTracks  || fPoolLambda->GetCurrentNEvents() >= fMinEventsToMix)
	{
		for (Int_t jMix=0; jMix<fPoolLambda->GetCurrentNEvents(); jMix++) 
			FillCorrelations(selectedTracksLeading, fPoolLambda->GetEvent(jMix),fHistMixLambda,0);
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
void AliLeadingV0Correlation::FillCorrelations(TObjArray* particles, TObjArray* mixed,TH3F*histo,TH3F*leadinfo)
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
					
					histo->Fill(vars[3],vars[0],vars[1]);
					if(leadinfo)leadinfo->Fill(vars[2],triggerEta,triggerParticle->Phi());

	           }
	        }
        }
}
//---------------------------------------------------------------------------------------
TObjArray*  AliLeadingV0Correlation::FindLeadingObjects(TObject *obj)
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
		if(!IsTrackNotFromV0(((AliAODTrack*)part)))continue;
		
		tracks->AddLast( part );
  	}
	// Order tracks by pT	
	QSortTracks( *tracks, 0, tracks->GetEntriesFast() );
	
	nTracks = tracks->GetEntriesFast();
	if( !nTracks ) return 0;
	
	return tracks;
}
//---------------------------------------------------------------------------------------
TObjArray*  AliLeadingV0Correlation::FindLeadingObjectsMC(TObject *obj)
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
		Bool_t isHadron = TMath::Abs(pdgCodeCurrent)==211 ||  // Pion
		TMath::Abs(pdgCodeCurrent)==2212 || // Proton
		TMath::Abs(pdgCodeCurrent)==321;    // Kaon
		if (!isHadron) continue;
		tracks->AddLast( part );
  	}
	// Order tracks by pT	
	QSortTracks( *tracks, 0, tracks->GetEntriesFast() );
	
	nTracks = tracks->GetEntriesFast();
	if( !nTracks ) return 0;
	
	return tracks;
}
//---------------------------------------------------------------------------------------
void AliLeadingV0Correlation::QSortTracks(TObjArray &a, Int_t first, Int_t last)
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
Bool_t AliLeadingV0Correlation::IsK0InvMass(const Double_t mass) const 
{
	
	const Float_t massK0            = 0.497;
	const Float_t sigmaK0           = 0.003;
	const Float_t nSigmaSignal      = 3.5; 
	
	return ((massK0-nSigmaSignal*sigmaK0)<=mass && mass<=(massK0 + nSigmaSignal*sigmaK0))?1:0;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsLambdaInvMass(const Double_t mass) const 
{
	
	const Float_t massLambda        = 1.116;
	const Float_t sigmaLambda       = 0.003;
	const Float_t nSigmaSignal      = 3.5; 
	
	return ((massLambda-nSigmaSignal*sigmaLambda)<=mass && mass<=(massLambda + nSigmaSignal*sigmaLambda))?1:0;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsTrackNotFromV0(AliAODTrack* track)
{
	Int_t atrID = track->GetID();

	for(int i=0; i<fAODEvent->GetNumberOfV0s(); i++){ // loop over V0s
		AliAODv0* aodV0 = fAODEvent->GetV0(i);
		
		AliAODTrack *trackPos=(AliAODTrack *)(aodV0->GetDaughter(0));
        AliAODTrack *trackNeg=(AliAODTrack *)(aodV0->GetDaughter(1));
			
		if ( !(IsAcseptedDaughterTrack(trackPos)) || !(IsAcseptedDaughterTrack(trackNeg)) ) continue;
		//----------------------------------
		Int_t negID = trackNeg->GetID();
		Int_t posID = trackPos->GetID();
		
		if ((TMath::Abs(negID)+1)==(TMath::Abs(atrID))){ return kFALSE;}
		if ((TMath::Abs(posID)+1)==(TMath::Abs(atrID))){ return kFALSE;}
		//----------------------------------
	}
	return kTRUE;
}
//---------------------------------------------------------------------------------------
Bool_t AliLeadingV0Correlation::IsAcseptedV0(const AliAODEvent*aod, const AliAODv0* aodV0, const AliAODTrack* myTrackPos, const AliAODTrack* myTrackNeg)
{
	if (!aodV0) return kFALSE;
	
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
void AliLeadingV0Correlation::Terminate(Option_t *)
{
	//No need in the grid
}
//---------------------------------------------------------------------------------------
