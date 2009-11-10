
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

//-----------------------------------------------------------------
//                 AliAnalysisTaskStrange class
//       This task is for single strange study from ESD/AOD
//          Origin: H.Ricaud, Helene.Ricaud@IReS.in2p3.fr
//-----------------------------------------------------------------
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"

#include "AliAnalysisTaskSE.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"

#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliAODv0.h"
#include "AliAODTrack.h"

#include "AliLog.h"

#include "AliAnalysisTaskStrange.h"

ClassImp(AliAnalysisTaskStrange)

//________________________________________________________________________
AliAnalysisTaskStrange::AliAnalysisTaskStrange() 
  : AliAnalysisTaskSE(), fAnalysisType("ESD"), fCollidingSystems(0), fUseCut("infoCut"), fListHist(),
    fHistPrimaryVertexPosX(0), fHistPrimaryVertexPosY(0), fHistPrimaryVertexPosZ(0),
    fHistTrackMultiplicity(0), fHistV0Multiplicity(0),
    fHistDcaPosToPrimVertex(0), fHistDcaNegToPrimVertex(0),
    fHistDcaPosToPrimVertexZoom(0), fHistDcaNegToPrimVertexZoom(0),
    fHistRadiusV0(0), fHistDecayLengthV0(0), fHistDcaV0Daughters(0), fHistChi2(0),
    fHistCosPointAngle(0), fHistCosPointAngleZoom(0),
    fHistV0MultiplicityOff(0),
    fHistPtVsYK0sOff(0), fHistPtVsYLambdaOff(0), fHistPtVsYAntiLambdaOff(0),
    fHistMassK0sOff(0), fHistMassLambdaOff(0), fHistMassAntiLambdaOff(0),
    fHistMassVsRadiusK0sOff(0), fHistMassVsRadiusLambdaOff(0), fHistMassVsRadiusAntiLambdaOff(0),
    fHistPtVsMassK0sOff(0), fHistPtVsMassLambdaOff(0), fHistPtVsMassAntiLambdaOff(0),
    fHistArmenterosPodolanskiOff(0),
    fHistV0MultiplicityOn(0),
    fHistPtVsYK0sOn(0), fHistPtVsYLambdaOn(0), fHistPtVsYAntiLambdaOn(0),
    fHistMassK0sOn(0), fHistMassLambdaOn(0), fHistMassAntiLambdaOn(0),
    fHistMassVsRadiusK0sOn(0), fHistMassVsRadiusLambdaOn(0), fHistMassVsRadiusAntiLambdaOn(0),
    fHistPtVsMassK0sOn(0), fHistPtVsMassLambdaOn(0), fHistPtVsMassAntiLambdaOn(0),
    fHistArmenterosPodolanskiOn(0)
{
  // Dummy constructor
}
//________________________________________________________________________
AliAnalysisTaskStrange::AliAnalysisTaskStrange(const char *name) 
  : AliAnalysisTaskSE(name), fAnalysisType("ESD"), fCollidingSystems(0), fUseCut("infocut"), fListHist(),
    fHistPrimaryVertexPosX(0), fHistPrimaryVertexPosY(0), fHistPrimaryVertexPosZ(0),
    fHistTrackMultiplicity(0), fHistV0Multiplicity(0),
    fHistDcaPosToPrimVertex(0), fHistDcaNegToPrimVertex(0),
    fHistDcaPosToPrimVertexZoom(0), fHistDcaNegToPrimVertexZoom(0),
    fHistRadiusV0(0), fHistDecayLengthV0(0), fHistDcaV0Daughters(0), fHistChi2(0),
    fHistCosPointAngle(0), fHistCosPointAngleZoom(0),
    fHistV0MultiplicityOff(0),
    fHistPtVsYK0sOff(0), fHistPtVsYLambdaOff(0), fHistPtVsYAntiLambdaOff(0),
    fHistMassK0sOff(0), fHistMassLambdaOff(0), fHistMassAntiLambdaOff(0),
    fHistMassVsRadiusK0sOff(0), fHistMassVsRadiusLambdaOff(0), fHistMassVsRadiusAntiLambdaOff(0),
    fHistPtVsMassK0sOff(0), fHistPtVsMassLambdaOff(0), fHistPtVsMassAntiLambdaOff(0),
    fHistArmenterosPodolanskiOff(0),
    fHistV0MultiplicityOn(0),
    fHistPtVsYK0sOn(0), fHistPtVsYLambdaOn(0), fHistPtVsYAntiLambdaOn(0),
    fHistMassK0sOn(0), fHistMassLambdaOn(0), fHistMassAntiLambdaOn(0),
    fHistMassVsRadiusK0sOn(0), fHistMassVsRadiusLambdaOn(0), fHistMassVsRadiusAntiLambdaOn(0),
    fHistPtVsMassK0sOn(0), fHistPtVsMassLambdaOn(0), fHistPtVsMassAntiLambdaOn(0),
    fHistArmenterosPodolanskiOn(0)
{
  // Constructor
  // Define output slots only here
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskStrange::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once

  fListHist = new TList();

  // Primary Vertex:
  fHistPrimaryVertexPosX       = new TH1F("h1PrimaryVertexPosX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexPosX);
  fHistPrimaryVertexPosY       = new TH1F("h1PrimaryVertexPosY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fListHist->Add(fHistPrimaryVertexPosY);
  fHistPrimaryVertexPosZ       = new TH1F("h1PrimaryVertexPosZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-2.0,2.0);
  fListHist->Add(fHistPrimaryVertexPosZ);

  // Multiplicity:
  if (!fHistTrackMultiplicity) {
    if (fCollidingSystems)
      fHistTrackMultiplicity = new TH1F("fHistTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 200, 0, 40000);
    else
      fHistTrackMultiplicity = new TH1F("fHistTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 250, 0, 250);
    fListHist->Add(fHistTrackMultiplicity);
  }
  if (!fHistV0Multiplicity) {
    if (fCollidingSystems)
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of V0s;Events", 200, 0, 40000);
    else
      fHistV0Multiplicity = new TH1F("fHistV0Multiplicity", "Multiplicity distribution;Number of V0s;Events", 50, 0, 50);
    fListHist->Add(fHistV0Multiplicity);
  }

  // Selection checks:
  fHistDcaPosToPrimVertex      = new TH2F("h2DcaPosToPrimVertex", "Positive V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertex);
  fHistDcaNegToPrimVertex      = new TH2F("h2DcaNegToPrimVertex", "Negative V0 daughter;dca(cm);Status",500,0,5,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertex);
  fHistDcaPosToPrimVertexZoom  = new TH2F("h2DcaPosToPrimVertexZoom", "Positive V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  fListHist->Add(fHistDcaPosToPrimVertexZoom);
  fHistDcaNegToPrimVertexZoom  = new TH2F("h2DcaNegToPrimVertexZoom", "Negative V0 daughter;dca(cm);Status",100,0,0.1,2,-0.5,1.5);
  fListHist->Add(fHistDcaNegToPrimVertexZoom);
  fHistRadiusV0                = new TH2F("h2RadiusV0", "Radius;Radius(cm);Status",1000,0,100,2,-0.5,1.5);
  fListHist->Add(fHistRadiusV0);
  fHistDecayLengthV0           = new TH2F("h2DecayLengthV0", "V0s decay Length;decay length(cm);Status", 200, 0, 100,2,-0.5,1.5);
  fListHist->Add(fHistDecayLengthV0);
  fHistDcaV0Daughters          = new TH2F("h2DcaV0Daughters", "DCA between daughters;dca(cm);Status", 160, 0, 4,2,-0.5,1.5);
  fListHist->Add(fHistDcaV0Daughters);
  fHistChi2                    = new TH2F("h2Chi2", "V0s chi2;chi2;Status", 33, 0, 33,2,-0.5,1.5);
  fListHist->Add(fHistChi2);
  fHistCosPointAngle           = new TH2F("h2CosPointAngle", "Cosine of V0's pointing angle", 100,0,1,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngle);
  fHistCosPointAngleZoom       = new TH2F("h2CosPointAngleZoom", "Cosine of V0's pointing angle", 100,0.9,1,2,-0.5,1.5);
  fListHist->Add(fHistCosPointAngleZoom);

  // bounds of histograms:
  // Radius
  const Double_t radius[10] = {0.0,2.5,2.9,3.9,7.6,15.0,23.9,37.8,42.8,100.0};
  Int_t nBinRadius        = 9;

  // V0 offline distributions
  if (!fHistV0MultiplicityOff) {
    if (fCollidingSystems)
      fHistV0MultiplicityOff = new TH1F("fHistV0MultiplicityOff", "Multiplicity distribution;Number of V0s;Events", 200, 0, 40000);
    else
      fHistV0MultiplicityOff = new TH1F("fHistV0MultiplicityOff", "Multiplicity distribution;Number of V0s;Events", 50, 0, 50); 
    fListHist->Add(fHistV0MultiplicityOff);
  }
  // Pt vs rapidity:
  fHistPtVsYK0sOff             = new TH2F("h2PtVsYK0sOff", "K^{0} Offline candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYK0sOff);
  fHistPtVsYLambdaOff          = new TH2F("h2PtVsYLambdaOff", "#Lambda^{0} Offline candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYLambdaOff);
  fHistPtVsYAntiLambdaOff      = new TH2F("h2PtVsYAntiLambdaOff", "#bar{#Lambda}^{0} Offline candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYAntiLambdaOff);
  // Mass:
  fHistMassK0sOff               = new TH1F("h1MassK0sOff", "K^{0} Offline candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistMassK0sOff);
  fHistMassLambdaOff            = new TH1F("h1MassLambdaOff", "#Lambda^{0} Offline candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassLambdaOff);
  fHistMassAntiLambdaOff          = new TH1F("h1MassAntiLambdaOff", "#bar{#Lambda}^{0} Offline candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambdaOff);
  // Mass vs radius:
  fHistMassVsRadiusK0sOff           = new TH2F("h2MassVsRadiusK0sOff", "K^{0} Offline candidates;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",nBinRadius,radius, 200, 0.4, 0.6);
  fListHist->Add(fHistMassVsRadiusK0sOff);
  fHistMassVsRadiusLambdaOff       = new TH2F("h2MassVsRadiusLambdaOff", "#Lambda Offline candidates;radius (cm);M(p#pi^{-}) (GeV/c^{2})",nBinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusLambdaOff);
  fHistMassVsRadiusAntiLambdaOff = new TH2F("h2MassVsRadiusAntiLambdaOff", "#bar{#Lambda} Offline candidates;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",nBinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusAntiLambdaOff);
  // Pt Vs Mass:
  fHistPtVsMassK0sOff             = new TH2F("h2PtVsMassK0sOff","K^{0} Offline candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,100,0,10);
  fListHist->Add(fHistPtVsMassK0sOff);
  fHistPtVsMassLambdaOff         = new TH2F("h2PtVsMassLambdaOff","#Lambda^{0} Offline candidates;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistPtVsMassLambdaOff);
  fHistPtVsMassAntiLambdaOff     = new TH2F("h2PtVsMassAntiLambdaOff","#bar{#Lambda}^{0} Offline candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistPtVsMassAntiLambdaOff);
  //ArmenterosPodolanski:
  fHistArmenterosPodolanskiOff   = new TH2F("h2ArmenterosPodolanskiOff","Armenteros-Podolanski Offline phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);
  fListHist->Add(fHistArmenterosPodolanskiOff);

  // V0 on-the-fly distributions
  if (!fHistV0MultiplicityOn) {
    if (fCollidingSystems)
      fHistV0MultiplicityOn = new TH1F("fHistV0MultiplicityOn", "Multiplicity distribution;Number of V0s;Events", 200, 0, 40000);
    else
      fHistV0MultiplicityOn = new TH1F("fHistV0MultiplicityOn", "Multiplicity distribution;Number of V0s;Events", 50, 0, 50);
    fListHist->Add(fHistV0MultiplicityOn);
  }
  // Pt vs rapidity:
  fHistPtVsYK0sOn              = new TH2F("h2PtVsYK0sOn", "K^{0} Onthefly candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYK0sOn);
  fHistPtVsYLambdaOn           = new TH2F("h2PtVsYLambdaOn", "#Lambda^{0} Onthefly candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYLambdaOn);
  fHistPtVsYAntiLambdaOn       = new TH2F("h2PtVsYAntiLambdaOn", "#bar{#Lambda}^{0} Onthefly candidates;p_{t} (GeV/c);rapidity",30,0,15,30,-1.5,1.5);
  fListHist->Add(fHistPtVsYAntiLambdaOn);
  // Mass:
  fHistMassK0sOn                = new TH1F("h1MassK0sOn", "K^{0} Onthefly candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 100, 0.4, 0.6);
  fListHist->Add(fHistMassK0sOn);
  fHistMassLambdaOn            = new TH1F("h1MassLambdaOn", "#Lambda^{0} Onthefly candidates;M(p#pi^{-}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassLambdaOn);
  fHistMassAntiLambdaOn        = new TH1F("h1MassAntiLambdaOn", "#bar{#Lambda}^{0} Onthefly candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});Counts", 75, 1.05, 1.2);
  fListHist->Add(fHistMassAntiLambdaOn);
  // Mass vs radius:
  fHistMassVsRadiusK0sOn         = new TH2F("h2MassVsRadiusK0sOn", "K^{0} Onthefly candidates;radius (cm);M(#pi^{+}#pi^{-}) (GeV/c^{2})",nBinRadius,radius, 200, 0.4, 0.6);
  fListHist->Add(fHistMassVsRadiusK0sOn);
  fHistMassVsRadiusLambdaOn     = new TH2F("h2MassVsRadiusLambdaOn", "#Lambda Onthefly candidates;radius (cm);M(p#pi^{-}) (GeV/c^{2})",nBinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusLambdaOn);
  fHistMassVsRadiusAntiLambdaOn  = new TH2F("h2MassVsRadiusAntiLambdaOn", "#bar{#Lambda} Onthefly candidates;radius (cm);M(#bar{p}#pi^{+}) (GeV/c^{2})",nBinRadius,radius, 140, 1.06, 1.2);
  fListHist->Add(fHistMassVsRadiusAntiLambdaOn);
  // Pt Vs Mass:
  fHistPtVsMassK0sOn              = new TH2F("h2PtVsMassK0sOn","K^{0} Onthefly candidates;M(#pi^{+}#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",200, 0.4, 0.6,100,0,10);
  fListHist->Add(fHistPtVsMassK0sOn);
  fHistPtVsMassLambdaOn          = new TH2F("h2PtVsMassLambdaOn","#Lambda^{0} Onthefly candidates;M(p#pi^{-}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistPtVsMassLambdaOn);
  fHistPtVsMassAntiLambdaOn      = new TH2F("h2PtVsMassAntiLambdaOn","#bar{#Lambda}^{0} Onthefly candidates;M(#bar{p}#pi^{+}) (GeV/c^{2});p_{t} (GeV/c)",140, 1.06, 1.2,100,0,10);
  fListHist->Add(fHistPtVsMassAntiLambdaOn);
  //ArmenterosPodolanski:
  fHistArmenterosPodolanskiOn    = new TH2F("h2ArmenterosPodolanskiOn","Armenteros-Podolanski Onthefly phase space;#alpha;p_{t} arm",100,-1.0,1.0,50,0,0.5);
  fListHist->Add(fHistArmenterosPodolanskiOn);
}

//________________________________________________________________________
void AliAnalysisTaskStrange::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  AliVEvent* lEvent = InputEvent();
  if (!lEvent) {
    Printf("ERROR: Event not available");
    return;
  }

  if (!(lEvent->GetNumberOfTracks())) {
    //Printf("Strange analysis task: There is no track in this event");
    return;
  }
  fHistTrackMultiplicity->Fill(lEvent->GetNumberOfTracks());

  Double_t tPrimaryVtxPosition[3];

  Int_t nv0s = 0;
  nv0s = lEvent->GetNumberOfV0s();
  //Printf("Strange analysis task: There are %d v0s in this event",nv0s);

  Int_t    lOnFlyStatus = 0, nv0sOn = 0, nv0sOff = 0;
  Double_t lChi2V0 = 0;
  Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
  Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
  Double_t lV0CosineOfPointingAngle = 0;
  Double_t lV0Radius = 0;
  Double_t lV0DecayLength = 0;
  Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
  Double_t lPt       = 0, lRapK0s = 0, lRapLambda = 0;
  Double_t lAlphaV0  = 0, lPtArmV0 = 0;

  Double_t  tV0Position[3];

  Double_t lMagneticField      = 999;


  //***********************
  // ESD loop
  //***********************

  if(fAnalysisType == "ESD") {

    const AliESDVertex *primaryVtx = ((AliESDEvent*)lEvent)->GetPrimaryVertex();
    tPrimaryVtxPosition[0] = primaryVtx->GetXv();
    tPrimaryVtxPosition[1] = primaryVtx->GetYv();
    tPrimaryVtxPosition[2] = primaryVtx->GetZv();

    fHistPrimaryVertexPosX->Fill(tPrimaryVtxPosition[0]);
    fHistPrimaryVertexPosY->Fill(tPrimaryVtxPosition[1]);
    fHistPrimaryVertexPosZ->Fill(tPrimaryVtxPosition[2]);

    lMagneticField = ((AliESDEvent*)lEvent)->GetMagneticField();
  

    for (Int_t iV0 = 0; iV0 < nv0s; iV0++)
      {// This is the begining of the V0 loop  
	AliESDv0 *v0 = ((AliESDEvent*)lEvent)->GetV0(iV0);
	if (!v0) continue;

	// AliAODtrack (V0 Daughters)
	UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
	UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

	AliESDtrack *pTrack = ((AliESDEvent*)lEvent)->GetTrack(lKeyPos);
	AliESDtrack *nTrack = ((AliESDEvent*)lEvent)->GetTrack(lKeyNeg);
	if (!pTrack || !nTrack) {
	  Printf("ERROR: Could not retreive one of the daughter track");
	  continue;
	}

	// Remove like-sign
	if ( pTrack->GetSign() == nTrack->GetSign()){
	  //cout<< "like sign, continue"<< endl;
	  continue;
	} 

	// Tracks quality cuts 
	if ( ( (pTrack->GetTPCNcls()) < 80 ) || ( (nTrack->GetTPCNcls()) < 80 ) ) continue;
	
	// TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
	if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;      
	if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;

	// DCA between daughter and Primary Vertex:
	if (pTrack) lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],lMagneticField) );

	if (nTrack) lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],lMagneticField) );

	// VO's main characteristics:
	lOnFlyStatus             = v0->GetOnFlyStatus();
	lChi2V0                  = v0->GetChi2V0();
	lDcaV0Daughters          = v0->GetDcaV0Daughters();
	lDcaV0ToPrimVertex       = v0->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
	lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1], tPrimaryVtxPosition[2]);
	v0->GetXYZ(tV0Position[0], tV0Position[1], tV0Position[2]);
	lV0Radius      = TMath::Sqrt(tV0Position[0]*tV0Position[0]+tV0Position[1]*tV0Position[1]);
	lV0DecayLength = TMath::Sqrt(TMath::Power(tV0Position[0] - tPrimaryVtxPosition[0],2) +
		                     TMath::Power(tV0Position[1] - tPrimaryVtxPosition[1],2) +
		                     TMath::Power(tV0Position[2] - tPrimaryVtxPosition[2],2 ));

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
	
	// Pt:
	lPt = v0->Pt();
	
	// Armenteros variables: !!
	lAlphaV0      = v0->AlphaV0();
	lPtArmV0      = v0->PtArmV0();
	
	// Selections:
	if (fUseCut.Contains("yes")) {
	  if ( (lDcaPosToPrimVertex      < 0.036 )||
	       (lDcaNegToPrimVertex      < 0.036 )||
	       (lDcaV0Daughters          > 0.5 )  ||
	       (lV0CosineOfPointingAngle < 0.99)     
	       ) continue;
	}
    
	// Filling histograms
	fHistDcaPosToPrimVertex->Fill(lDcaPosToPrimVertex,lOnFlyStatus);
	fHistDcaNegToPrimVertex->Fill(lDcaNegToPrimVertex,lOnFlyStatus);
	fHistDcaPosToPrimVertexZoom->Fill(lDcaPosToPrimVertex,lOnFlyStatus);
	fHistDcaNegToPrimVertexZoom->Fill(lDcaNegToPrimVertex,lOnFlyStatus);
	fHistRadiusV0->Fill(lV0Radius,lOnFlyStatus);
	fHistDecayLengthV0->Fill(lV0DecayLength,lOnFlyStatus);
	fHistDcaV0Daughters->Fill(lDcaV0Daughters,lOnFlyStatus);
	fHistChi2->Fill(lChi2V0,lOnFlyStatus);
	fHistCosPointAngle->Fill(lV0CosineOfPointingAngle,lOnFlyStatus);
	if (lV0CosineOfPointingAngle >= 0.9) fHistCosPointAngleZoom->Fill(lV0CosineOfPointingAngle,lOnFlyStatus);
	if(!lOnFlyStatus){
	  nv0sOff++;
	  fHistPtVsYK0sOff->Fill(lPt,lRapK0s);
	  fHistPtVsYLambdaOff->Fill(lPt,lRapLambda);
	  fHistPtVsYAntiLambdaOff->Fill(lPt,lRapLambda);
	  fHistArmenterosPodolanskiOff->Fill(lAlphaV0,lPtArmV0);
	}
	else {
	  nv0sOn++;
	  fHistPtVsYK0sOn->Fill(lPt,lRapK0s);
	  fHistPtVsYLambdaOn->Fill(lPt,lRapLambda);
	  fHistPtVsYAntiLambdaOn->Fill(lPt,lRapLambda);
	  fHistArmenterosPodolanskiOn->Fill(lAlphaV0,lPtArmV0);
	}
	// K0s invariant mass histograms:
	if (TMath::Abs(lRapK0s) < 1) {  
	  if(!lOnFlyStatus){
	    fHistMassK0sOff->Fill(lInvMassK0s);
	    fHistMassVsRadiusK0sOff->Fill(lV0Radius,lInvMassK0s);
	    fHistPtVsMassK0sOff->Fill(lInvMassK0s,lPt);
	  }
	  else {
	    fHistMassK0sOn->Fill(lInvMassK0s);
	    fHistMassVsRadiusK0sOn->Fill(lV0Radius,lInvMassK0s);
	    fHistPtVsMassK0sOn->Fill(lInvMassK0s,lPt);
	  }
	}
	// Lambda and AntiLambda invariant mass histograms:
	if (TMath::Abs(lRapLambda) < 1) {
	  if(!lOnFlyStatus){
	    fHistMassLambdaOff->Fill(lInvMassLambda);
	    fHistMassAntiLambdaOff->Fill(lInvMassAntiLambda);
	    fHistMassVsRadiusLambdaOff->Fill(lV0Radius,lInvMassLambda);
	    fHistMassVsRadiusAntiLambdaOff->Fill(lV0Radius,lInvMassAntiLambda);
	    fHistPtVsMassLambdaOff->Fill(lInvMassLambda,lPt);
	    fHistPtVsMassAntiLambdaOff->Fill(lInvMassAntiLambda,lPt);
	  }
	  else {
	    fHistMassLambdaOn->Fill(lInvMassLambda);
	    fHistMassAntiLambdaOn->Fill(lInvMassAntiLambda);
	    fHistMassVsRadiusLambdaOn->Fill(lV0Radius,lInvMassLambda);
	    fHistMassVsRadiusAntiLambdaOn->Fill(lV0Radius,lInvMassAntiLambda);
	    fHistPtVsMassLambdaOn->Fill(lInvMassLambda,lPt);
	    fHistPtVsMassAntiLambdaOn->Fill(lInvMassAntiLambda,lPt);
	  }
	}
      } // end V0 loop

  }

  //***********************
  // AOD loop
  //***********************

  else if(fAnalysisType == "AOD") {

    const AliAODVertex *primaryVtx = ((AliAODEvent*)lEvent)->GetPrimaryVertex();
    tPrimaryVtxPosition[0] = primaryVtx->GetX();
    tPrimaryVtxPosition[1] = primaryVtx->GetY();
    tPrimaryVtxPosition[2] = primaryVtx->GetZ();

    fHistPrimaryVertexPosX->Fill(tPrimaryVtxPosition[0]);
    fHistPrimaryVertexPosY->Fill(tPrimaryVtxPosition[1]);
    fHistPrimaryVertexPosZ->Fill(tPrimaryVtxPosition[2]);
  
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) 
      {// This is the begining of the V0 loop
	AliAODv0 *myAODv0 = ((AliAODEvent*)lEvent)->GetV0(iV0);
	if (!myAODv0) continue;

	// common part
	lV0Radius                = myAODv0->RadiusV0();
	lDcaPosToPrimVertex      = myAODv0->DcaPosToPrimVertex();
	lDcaNegToPrimVertex      = myAODv0->DcaNegToPrimVertex();
	lOnFlyStatus             = myAODv0->GetOnFlyStatus();
	lChi2V0                  = myAODv0->Chi2V0();
	lDcaV0Daughters          = myAODv0->DcaV0Daughters();
	lDcaV0ToPrimVertex       = myAODv0->DcaV0ToPrimVertex();
	lV0DecayLength           = myAODv0->DecayLengthV0(tPrimaryVtxPosition);
	lV0CosineOfPointingAngle = myAODv0->CosPointingAngle(tPrimaryVtxPosition);

	lInvMassK0s        = myAODv0->MassK0Short();
	lInvMassLambda     = myAODv0->MassLambda();
	lInvMassAntiLambda = myAODv0->MassAntiLambda();

	lPt        = TMath::Sqrt(myAODv0->Pt2V0());
	lRapK0s    = myAODv0->RapK0Short();
	lRapLambda = myAODv0->RapLambda();
	lAlphaV0   = myAODv0->AlphaV0();
	lPtArmV0   = myAODv0->PtArmV0();


	// Selections:
	if (fUseCut.Contains("yes")) {
	  if ( (lDcaPosToPrimVertex      < 0.036 )||
	       (lDcaNegToPrimVertex      < 0.036 )||
	       (lDcaV0Daughters          > 0.5 )  ||
	       (lV0CosineOfPointingAngle < 0.99)     
	       ) continue;
	}
    
	// Filling histograms
	fHistDcaPosToPrimVertex->Fill(lDcaPosToPrimVertex,lOnFlyStatus);
	fHistDcaNegToPrimVertex->Fill(lDcaNegToPrimVertex,lOnFlyStatus);
	fHistDcaPosToPrimVertexZoom->Fill(lDcaPosToPrimVertex,lOnFlyStatus);
	fHistDcaNegToPrimVertexZoom->Fill(lDcaNegToPrimVertex,lOnFlyStatus);
	fHistRadiusV0->Fill(lV0Radius,lOnFlyStatus);
	fHistDecayLengthV0->Fill(lV0DecayLength,lOnFlyStatus);
	fHistDcaV0Daughters->Fill(lDcaV0Daughters,lOnFlyStatus);
	fHistChi2->Fill(lChi2V0,lOnFlyStatus);
	fHistCosPointAngle->Fill(lV0CosineOfPointingAngle,lOnFlyStatus);
	if (lV0CosineOfPointingAngle >= 0.9) fHistCosPointAngleZoom->Fill(lV0CosineOfPointingAngle,lOnFlyStatus);
	if(!lOnFlyStatus){
	  nv0sOff++;
	  fHistPtVsYK0sOff->Fill(lPt,lRapK0s);
	  fHistPtVsYLambdaOff->Fill(lPt,lRapLambda);
	  fHistPtVsYAntiLambdaOff->Fill(lPt,lRapLambda);
	  fHistArmenterosPodolanskiOff->Fill(lAlphaV0,lPtArmV0);
	}
	else {
	  nv0sOn++;
	  fHistPtVsYK0sOn->Fill(lPt,lRapK0s);
	  fHistPtVsYLambdaOn->Fill(lPt,lRapLambda);
	  fHistPtVsYAntiLambdaOn->Fill(lPt,lRapLambda);
	  fHistArmenterosPodolanskiOn->Fill(lAlphaV0,lPtArmV0);
	}
	// K0s invariant mass histograms:
	if (TMath::Abs(lRapK0s) < 1) {  
	  if(!lOnFlyStatus){
	    fHistMassK0sOff->Fill(lInvMassK0s);
	    fHistMassVsRadiusK0sOff->Fill(lV0Radius,lInvMassK0s);
	    fHistPtVsMassK0sOff->Fill(lInvMassK0s,lPt);
	  }
	  else {
	    fHistMassK0sOn->Fill(lInvMassK0s);
	    fHistMassVsRadiusK0sOn->Fill(lV0Radius,lInvMassK0s);
	    fHistPtVsMassK0sOn->Fill(lInvMassK0s,lPt);
	  }
	}
	// Lambda and AntiLambda invariant mass histograms:
	if (TMath::Abs(lRapLambda) < 1) {
	  if(!lOnFlyStatus){
	    fHistMassLambdaOff->Fill(lInvMassLambda);
	    fHistMassAntiLambdaOff->Fill(lInvMassAntiLambda);
	    fHistMassVsRadiusLambdaOff->Fill(lV0Radius,lInvMassLambda);
	    fHistMassVsRadiusAntiLambdaOff->Fill(lV0Radius,lInvMassAntiLambda);
	    fHistPtVsMassLambdaOff->Fill(lInvMassLambda,lPt);
	    fHistPtVsMassAntiLambdaOff->Fill(lInvMassAntiLambda,lPt);
	  }
	  else {
	    fHistMassLambdaOn->Fill(lInvMassLambda);
	    fHistMassAntiLambdaOn->Fill(lInvMassAntiLambda);
	    fHistMassVsRadiusLambdaOn->Fill(lV0Radius,lInvMassLambda);
	    fHistMassVsRadiusAntiLambdaOn->Fill(lV0Radius,lInvMassAntiLambda);
	    fHistPtVsMassLambdaOn->Fill(lInvMassLambda,lPt);
	    fHistPtVsMassAntiLambdaOn->Fill(lInvMassAntiLambda,lPt);
	  }
	}
      } // end V0 loop
  }

  fHistV0Multiplicity->Fill(nv0s);
  fHistV0MultiplicityOff->Fill(nv0sOff);
  fHistV0MultiplicityOn->Fill(nv0sOn);

  // Post output data.
  PostData(1, fListHist);
}    

//________________________________________________________________________
void AliAnalysisTaskStrange::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
}
