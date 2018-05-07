/**************************************************************************
* Copyright(c) 1998-2008,ALICE Experiment at CERN,All rights reserved.    *
*                                                                         *
* Author: The ALICE Off-line Project.                                     *
* Contributors are mentioned in the code where appropriate.               *
*                                                                         *
* Permission to use,copy,modify and distribute this software and its      *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee,provided that the above copyright notice appears in all     *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

/////////////////////////////////////////////////////
// AliAnalysisTaskFlowQnSPCascade:
// Analysis task to select Xi and Omega candidates for flow analysis.
//
// Author: ya.zhu@cern.ch
//////////////////////////////////////////////////////

#include <iostream>
#include "AliQnCorrectionsCutsSet.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsHistos.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliAnalysisTaskFlowVectorCorrections.h"
#include <TTimeStamp.h>
#include <TStopwatch.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include "AliMultSelection.h"
#include <AliCentrality.h>
#include <AliESDEvent.h>
#include <TList.h>
#include <TGraphErrors.h>
#include "TFile.h"
#include "TChain.h"
#include "AliEventCuts.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliVVertex.h"
#include "AliAODcascade.h"
#include "AliAODEvent.h"
#include "AliMultiplicity.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowTrackCuts.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliPIDResponse.h"
#include "AliOADBContainer.h"
#include "AliAnalysisTaskFlowQnSPCascade.h"
#include "TMath.h"
#include "AliLog.h"
#include "AliESDVertex.h"
#include "AliVEvent.h"
#include "TGrid.h"
#include  "TH1D.h"
#include "AliVEventHandler.h"
#include "AliProdInfo.h"
using namespace std;

ClassImp(AliAnalysisTaskFlowQnSPCascade)

//=======================================================================================
  AliAnalysisTaskFlowQnSPCascade::AliAnalysisTaskFlowQnSPCascade():
    AliAnalysisTaskSE(),
    fEventCuts(0x0),
    fPIDResponse(0x0),
    fFlowQnVectorMgr(0x0),
    fIsMultiStrange(0),
    fIsStrange(0),
    fHistList(0x0),
    fMinCent(0),
    fMaxCent(0),
    fHarmonicOrder(0),
    fWeightCorrection(0),
    hOmegaWeight(0x0),
    hXiWeight(0x0),

//-----------------------------para-----------------------------------------------
    XiPse(0),
    V0RadiusXi(0),
    XiRadius(0),
    dcaXiDaughters(0),
    XiCosOfPointingAngle(0),
    dcaV0ToPrimaryVtxXi(0),
    dcaBachToPrimaryVtxXi(0),
    invMassLambdaAsCascDghter(0),
    dcaV0DaughtersXi(0),
    V0CosOfPointingAngleXi(0),
    dcaPosToPrimaryVtxXi(0),
    dcaNegToPrimaryVtxXi(0),
    fV0DCAdaughters(0),
    fV0CosinePointingAngle(0),
    fDecayRad(0),
    fV0DCAToPrimVertex(0),
    V0LifeTime(0),            
    fOnline(kFALSE),
    fSepAnalysis(kFALSE),
    fUseHybridGlobalTrack(kFALSE),

//-----------------------Xi--------------------------------------------------
    fXiPseMin(0),
    fXiPseMax(0),
    fV0RadiusXiMin(0),
    fV0RadiusXiMax(0),
    fXiRadiusMin(0),
    fXiRadiusMax(0),
    fdcaXiDaughtersMax(0),
    fXiCosOfPointingAngleMin(0),
    fdcaV0ToPrimaryVtxXiMin(0),
    fdcaBachToPrimaryVtxXiMin(0),
    fLambdaMassWind(0),
    fdcaV0DaughtersXi(0),
    fV0CosOfPointingAngleXiMin(0),
    fdcaPosToPrimaryVtxXiMin(0),
    fdcaNegToPrimaryVtxXiMin(0),
//---------V0---------------------------------------------------
   // tV0Eta(0),
   // tV0Pt(0),
    tDecayRapidity(0),
    tDecayRad(0),
    tV0DCAdaughtersMax(0),
    tV0CosinePointingAngleMin(0),
    tV0DCAToPrimVertexMin(0),
    tV0LifeTimeMax(0),
    tDecayLengthV0(0),
//--------track-------------------------------------------------  
    fPrimaryTrackEta(0),
    fTrackEta(0),
    fTrackPtMin(0),
    fTPCNcls(0),
    fRPFromTPC(0),

//-------Calibration----------------
    fRemChV0A(0),
    fMultV0(0),
    fNHarm(2.),
    fQxnmV0A(0),
    fQynmV0A(0),
    fQxnsV0A(0),
    fQynsV0A(0),
    fQxnmV0C(0),
    fQynmV0C(0),
    fQxnsV0C(0),
    fQynsV0C(0),
    fRun(-1),

//--------PID-------------------------------------------------        
    fXiPIDsigma(0),
    fV0PIDsigma(0),
    fExpectedCorrectionPass(0),
    fAlternativeCorrectionPass(0){
        for(Int_t i=0;i<20;i++){
            weightXi[i]= 0;
            weightOmega[i]= 0;
        }
        
    }

//================================================================================

//=======================================================================================
AliAnalysisTaskFlowQnSPCascade::
AliAnalysisTaskFlowQnSPCascade(const char *name, double centMin, double centMax,bool WeightCorrection) :
    AliAnalysisTaskSE(name),
    fEventCuts(0x0),
    fPIDResponse(0x0),
    fFlowQnVectorMgr(0x0),
    fIsMultiStrange(0),
    fIsStrange(0),
    fHistList(0x0),
    fMinCent(centMin),
    fMaxCent(centMax),
    fHarmonicOrder(0),
    fWeightCorrection(WeightCorrection),
    hOmegaWeight(0x0),
    hXiWeight(0x0), 
   
//-----------------------------para-----------------------------------------------
    XiPse(0),
    V0RadiusXi(0),
    XiRadius(0),
    dcaXiDaughters(0),
    XiCosOfPointingAngle(0),
    dcaV0ToPrimaryVtxXi(0),
    dcaBachToPrimaryVtxXi(0),
    invMassLambdaAsCascDghter(0),
    dcaV0DaughtersXi(0),
    V0CosOfPointingAngleXi(0),
    dcaPosToPrimaryVtxXi(0),
    dcaNegToPrimaryVtxXi(0),
    fV0DCAdaughters(0),
    fV0CosinePointingAngle(0),
    fDecayRad(0),
    fV0DCAToPrimVertex(0),
    V0LifeTime(0), 
    fOnline(kFALSE),
    fSepAnalysis(kFALSE),
    fUseHybridGlobalTrack(kFALSE),
//-----------------------Xi--------------------------------------------------
    fXiPseMin(0),
    fXiPseMax(0),
    fV0RadiusXiMin(0),
    fV0RadiusXiMax(0),
    fXiRadiusMin(0),
    fXiRadiusMax(0),
    fdcaXiDaughtersMax(0),
    fXiCosOfPointingAngleMin(0),
    fdcaV0ToPrimaryVtxXiMin(0),
    fdcaBachToPrimaryVtxXiMin(0),
    fLambdaMassWind(0),
    fdcaV0DaughtersXi(0),
    fV0CosOfPointingAngleXiMin(0),
    fdcaPosToPrimaryVtxXiMin(0),
    fdcaNegToPrimaryVtxXiMin(0),
//---------V0---------------------------------------------------
  //  tV0Eta(0),
  //  tV0Pt(0),
    tDecayRapidity(0),
    tDecayRad(0),
    tV0DCAdaughtersMax(0),
    tV0CosinePointingAngleMin(0),
    tV0DCAToPrimVertexMin(0),
    tV0LifeTimeMax(0),
    tDecayLengthV0(0),
//--------track-------------------------------------------------  
    fPrimaryTrackEta(0),
    fTrackEta(0),
    fTrackPtMin(0),
    fTPCNcls(0),
    fRPFromTPC(0),

//-------Calibration----------------
    fRemChV0A(0),
    fMultV0(0),
    fNHarm(2.),
    fQxnmV0A(0),
    fQynmV0A(0),
    fQxnsV0A(0),
    fQynsV0A(0),
    fQxnmV0C(0),
    fQynmV0C(0),
    fQxnsV0C(0),
    fQynsV0C(0),

    fRun(-1),
//--------PID-------------------------------------------------        
    fXiPIDsigma(0),
    fV0PIDsigma(0),

    fExpectedCorrectionPass(0),
    fAlternativeCorrectionPass(0){
        
        for(Int_t i=0;i<20;i++){
            weightXi[i]= 0;
            weightOmega[i]= 0;
        }

  DefineInput( 0,TChain::Class());
  
  DefineInput(1, TList::Class());

if(fWeightCorrection){
  DefineInput(2, TList::Class());
}
  DefineOutput(1, TList::Class());

}

//========================================================================================
AliAnalysisTaskFlowQnSPCascade::~AliAnalysisTaskFlowQnSPCascade()
{
  if(fHistList)  delete fHistList;
  if (fEventCuts) {
	         delete fEventCuts;
  fEventCuts = 0x0;
   }
}

//==================================Event QA==============================================

void AliAnalysisTaskFlowQnSPCascade::QAEventInput()
{


  TList *tQAEvents = new TList();
  tQAEvents->SetName("Events");
  tQAEvents->SetOwner();
//----------------------------------------------------------------------------------------  
  TH1F *fhEventCent
    = new TH1F( "EventCent",
                "Event Number ; count ; Number of Event",
                1, 0.,1);

  tQAEvents->Add(fhEventCent);


   TH1F *fhEventCentAfterPilp
    = new TH1F( "EventCentAfterPilp",
                "Event distribution to centrality after pile up remove; Centrality ; Number of Event",
                70, 0., 70);
  tQAEvents->Add(fhEventCentAfterPilp);

  fHistList->Add(tQAEvents);

}
//======================================Track QA=======================================

void AliAnalysisTaskFlowQnSPCascade::QATrack()
{


  TList *tQATrack = new TList();
  tQATrack->SetName("Track");
  tQATrack->SetOwner();
//----------------------------------------------------------------------------------------  


  TH1F *fh1PhiGlobleTrack
    = new TH1F("PhiGlobleTrack",
               "#phi of Globle tracks (around the mass peak); #phi (deg); Counts", 64, 0., 6.4);
  tQATrack->Add(fh1PhiGlobleTrack);

  TH1F *fh1PhiHyBirdsTrack
    = new TH1F("PhiHyBirdsTrack",
               "#phi of HyBirds tracks (around the mass peak); #phi (deg); Counts", 64, 0., 6.4);
  tQATrack->Add(fh1PhiHyBirdsTrack);


  fHistList->Add(tQATrack);


}




//======================================Event plane====================================================

void AliAnalysisTaskFlowQnSPCascade::QAEventPlane()
{
  TList *tQAEP = new TList();
  tQAEP->SetName("EventPlane");
  tQAEP->SetOwner();

    TH1F *fhEPangleV0A = new TH1F("hEPangleV0A",
                                  "EP from V0A; #Psi; Number of Events",
                                  640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0A);
    
    TH1F *fhEPangleV0C = new TH1F("hEPangleV0C",
                                  "EP from V0C; #Psi; Number of Events",
                                  640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0C);
    
    TH1F *fhEPangleTPC = new TH1F("hEPangleTPC",
                                  "EP from TPC; #Psi; Number of Events",
                                  640, -3.2, 3.2);
    tQAEP->Add(fhEPangleTPC);
    
    
    TH1F *fhEPangleV0AAfterRescale = new TH1F("hEPangleV0AAfterRescale",
                                              "EP from V0A; #Psi; Number of Events",
                                              640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0AAfterRescale);
    
    TH1F *fhEPangleV0CAfterRescale = new TH1F("hEPangleV0CAfterRescale",
                                              "EP from V0C; #Psi; Number of Events",
                                              640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0CAfterRescale);
    
    TH1F *fhEPangleV0AAfterRaw = new TH1F("hEPangleV0AAfterRaw",
                                          "EP from V0A; #Psi; Number of Events",
                                          640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0AAfterRaw);
    
    TH1F *fhEPangleV0CAfterRaw = new TH1F("hEPangleV0CAfterRaw",
                                          "EP from V0C; #Psi; Number of Events",
                                          640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0CAfterRaw);
    
    TH1F *fhEPangleTPCAfterPlain = new TH1F("hEPangleTPCAfterplain",
                                            "EP from TPC; #Psi; Number of Events",
                                            640, -3.2, 3.2);
    tQAEP->Add(fhEPangleTPCAfterPlain);
    //---------------------------step 1-------------------------
    
    TH1F *fhEPangleV0AAfterPlain = new TH1F("hEPangleV0AAfterPlain",
                                            "EP from V0A; #Psi; Number of Events",
                                            640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0AAfterPlain);
    
    TH1F *fhEPangleV0CAfterPlain = new TH1F("hEPangleV0CAfterPlain",
                                            "EP from V0C; #Psi; Number of Events",
                                            640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0CAfterPlain);
    
    TH1F *fhEPangleTPCAfterRecent = new TH1F("hEPangleTPCAfterRecent",
                                             "EP from TPC; #Psi; Number of Events",
                                             640, -3.2, 3.2);
    tQAEP->Add(fhEPangleTPCAfterRecent);
    
    //-------------------------step 2------------------------------
    TH1F *fhEPangleV0AAfterRecent = new TH1F("hEPangleV0AAfterRecent",
                                             "EP from V0A; #Psi; Number of Events",
                                             640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0AAfterRecent);
    
    TH1F *fhEPangleV0CAfterRecent = new TH1F("hEPangleV0CAfterRecent",
                                             "EP from V0C; #Psi; Number of Events",
                                             640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0CAfterRecent);
    
    
    //-------------------------step 3------------------------------
    TH1F *fhEPangleV0AAfterAlign = new TH1F("hEPangleV0AAfterAlign",
                                            "EP from V0A; #Psi; Number of Events",
                                            640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0AAfterAlign);
    
    TH1F *fhEPangleV0CAfterAlign = new TH1F("hEPangleV0CAfterAlign",
                                            "EP from V0C; #Psi; Number of Events",
                                            640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0CAfterAlign);
    
    TH1F *fhEPangleTPCAfterTwist = new TH1F("hEPangleTPCAfterTwist",
                                            "EP from TPC; #Psi; Number of Events",
                                            640, -3.2, 3.2);
    tQAEP->Add(fhEPangleTPCAfterTwist);
    
    
    
    TH1F *fhEPangleV0AAfterTwist = new TH1F("hEPangleV0AAfterTwist",
                                            "EP from V0A; #Psi; Number of Events",
                                            640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0AAfterTwist);
    
    TH1F *fhEPangleV0CAfterTwist = new TH1F("hEPangleV0CAfterTwist",
                                            "EP from V0C; #Psi; Number of Events",
                                            640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0CAfterTwist);
    
    
    
    //-------------------------step 4--------------------------------
    TH1F *fhEPangleV0ALatest = new TH1F("hEPangleV0ALatest",
                                        "EP from V0A; #Psi; Number of Events",
                                        640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0ALatest);
    
    TH1F *fhEPangleV0CLatest = new TH1F("hEPangleV0CLatest",
                                        "EP from V0C; #Psi; Number of Events",
                                        640, -3.2, 3.2);
    tQAEP->Add(fhEPangleV0CLatest);
    
    TH1F *fhEPangleTPCLatest = new TH1F("hEPangleTPCLatest",
                                        "EP from TPC; #Psi; Number of Events",
                                        640, -3.2, 3.2);
    tQAEP->Add(fhEPangleTPCLatest);
//-------------------------------------------------------

  TProfile *fProfQnProduct = new TProfile("hProfQnProduct",
                                    "Qn product for different sub event",
				                                     3, 0.5, 3.5);
  (fProfQnProduct->GetXaxis())
  ->SetBinLabel(1, "Qn_{B}*Qn_{C}^{*}");
  (fProfQnProduct->GetXaxis())
  ->SetBinLabel(2, "Qn_{A}*Qn_{B}^{*}");
  (fProfQnProduct->GetXaxis())
  ->SetBinLabel(3, "Qn_{A}*Qn_{C}^{*}");
  tQAEP->Add(fProfQnProduct);

  fHistList->Add(tQAEP);
}

//====================================Cascade candidates============================================

void AliAnalysisTaskFlowQnSPCascade::QACascadeCandidates()
{
  
  TList *tQACascade = new TList();
  tQACascade->SetName("Cascade");
  tQACascade->SetOwner();

  TH1F *fh1DCAXiDaughters
    = new TH1F( "DcaXiDaughters",
                "DCA between Xi Daughters; DCA (cm); Number of Cascades",
                1000, 0., 0.5);
  tQACascade->Add(fh1DCAXiDaughters);

  TH1F *fh1DCABachToPrimVertex
    = new TH1F("DcaBachToPrimVertex",
               "DCA of Bach. to Prim. Vertex; DCA (cm);Number of Cascades",
               2500, 0., 2.5);
  tQACascade->Add(fh1DCABachToPrimVertex);

  TH1F *fh1XiCosOfPointingAngle
    = new TH1F("XiCosineOfPointingAngle",
               "Cos of Xi Pointing Angle; Cos (Xi Point.Angl);Number of Xis",
               200, 0.99, 1.0);
  tQACascade->Add(fh1XiCosOfPointingAngle);

  TH1F *fh1V0CosOfPointingAngle
    = new TH1F("V0CosOfPointingAngleXi",
               "Cos of V0 Pointing Angle, in cascade;Cos(V0 Point. Angl); Counts",
               200, 0.98, 1.0);
  tQACascade->Add(fh1V0CosOfPointingAngle);

  TH1F *fh1V0Radius  = new TH1F("V0RadiusXi",
                          "V0 decay radius, in cascade; radius (cm); Counts",
                          1050, 0., 105.0);
  tQACascade->Add(fh1V0Radius);

  TH1F *fh1DcaV0DaughtersXi = new TH1F("DcaV0DaughtersXi",
                                 "DCA between V0 daughters, in cascade;DCA (cm);Number of V0s", 120, 0., 6);
  tQACascade->Add(fh1DcaV0DaughtersXi);

  TH1F *fh1DcaV0ToPrimVertex = new TH1F("DcaV0ToPrimVertexXi",
                                  "DCA of V0 to Prim. Vertex, in cascade;DCA (cm);Number of Cascades", 200, 0., 1.);
  tQACascade->Add(fh1DcaV0ToPrimVertex);

  TH1F *fh1DCAPosToPrimVertex =
     new TH1F("DcaPosToPrimVertexXi",
              "DCA of V0 pos daughter to Prim. Vertex;DCA (cm);Counts",
              3000, 0, 3);
  tQACascade->Add(fh1DCAPosToPrimVertex);

  TH1F *fh1DCANegToPrimVertex
    =  new TH1F("DcaNegToPrimVertexXi",
                "DCA of V0 neg daughter to Prim. Vertex;DCA (cm);Counts",
                3000, 0, 3);
  tQACascade->Add(fh1DCANegToPrimVertex);
/*
  TH1F *fh1V0toXiCosOfPointingAngle
    = new TH1F("V0toXiCosineOfPointingAngle",
               "Cos. of V0 Ptng Angl Xi vtx; Cos(V0 Point. Angl / Xi vtx); Counts",
               100, 0.99, 1.0);
  tQACascade->Add(fh1V0toXiCosOfPointingAngle);
*/
  TH1F *fh1XiPt
    = new TH1F("XiPt" ,
               "#Xi Pt (cand. around the mass peak);p_{t}(#Xi)(GeV/c);Counts",
               100, 0.0, 10.0);
  tQACascade->Add(fh1XiPt);

  TH1F *fh1OmegaPt
    = new TH1F("OmegaPt" ,
               "#Omega Pt (cand. around the mass peak);p_{t}(#Xi)(GeV/c);Counts",
               100, 0.0, 10.0);
  tQACascade->Add(fh1OmegaPt);

  TH1F *fh1PhiXi
    = new TH1F("PhiXi",
               "#phi of #Xi candidates (around the mass peak); #phi (deg); Counts", 64, 0., 6.4);
  tQACascade->Add(fh1PhiXi);


  TH1F *fh1PhiOmega
    = new TH1F("PhiOmega",
               "#phi of #Omega candidates (around the mass peak); #phi (deg); Counts", 64, 0., 6.4);
  tQACascade->Add(fh1PhiOmega);

  TH1F *fhXiPesRap = new TH1F("hXiPesrapidity",
                          "#Xi Pesudorapidity distribution before rap. cut; y; Number of counts",
                                                    200, -1., 1.);
  tQACascade->Add(fhXiPesRap);


  TH1F *fhOmegaPesRap = new TH1F("hOmegaPesrapidity",
                          "#Omega Pesudorapidity distribution before rap. cut; y; Number of counts",
                                                    200, -1., 1.);
  tQACascade->Add(fhOmegaPesRap);

  TH2F *fh2TPCdEdxOfPion
    = new TH2F( "TPCdEdxOfPion",
                " TPC dE/dx of Pion; P(GeV/c); TPC signal (ADC) ",
                2000, 0, 10.0, 450, 0., 900.);
  tQACascade->Add(fh2TPCdEdxOfPion);

  TH2F *fh2TPCdEdxOfProton = new TH2F( "TPCdEdxOfProton",
                " TPC dE/dx of Proton; P(GeV/c); TPC signal (ADC) ",
                2000, 0, 10.0, 450, 0., 900.);
  tQACascade->Add(fh2TPCdEdxOfProton);

  TH2F *fh2TPCdEdxOfBachelorKaon = new TH2F( "TPCdEdxOfBachelorKaon",
                " TPC dE/dx of Bachelor Kaon; P(GeV/c); TPC signal (ADC) ",
                2000, 0, 10.0, 450, 0., 900.);
  tQACascade->Add(fh2TPCdEdxOfBachelorKaon);

  TH2F *fh2TPCdEdxOfBachelorPion = new TH2F( "TPCdEdxOfBachelorPion",
                " TPC dE/dx of Bachelor Pion; P(GeV/c); TPC signal (ADC) ",
                2000, 0, 10.0, 450, 0., 900.);
  tQACascade->Add(fh2TPCdEdxOfBachelorPion);

  TH2F *th2Armenteros
    = new TH2F("Armenteros",
               "#alpha_{Arm}(casc. cand.) Vs Pt_{Arm}(casc. cand.); #alpha_{Arm} ; Pt_{Arm} (GeV/c)",
               140, -1.2, 1.2, 300, 0., 0.3);
  tQACascade->Add(th2Armenteros);

  fHistList->Add(tQACascade);

}

//==========================================V0 candidate===========================================================
void AliAnalysisTaskFlowQnSPCascade::QAV0Candidates()
{

   TList *tQAV0 = new TList();
   tQAV0->SetName("V0");
   tQAV0->SetOwner();
  
   TH2F *fh2PhiV0
    = new TH2F("PhiV0com",
               "#phi of #V0 candidates (around the mass peak); Counts of phi; Counts of phi after", 64, 0., 3.2, 64, 0. , 3.2);
   tQAV0->Add(fh2PhiV0);

   TH1F *hDecayAlpha  = new TH1F( "DecayAlpha",
                "hDecayAlpha; Alpha; Number of Lambda",
		                100, -2.0, 2.0);
   tQAV0->Add(hDecayAlpha);
   
   TH2F *hQtVsAlpha
    = new TH2F("QtVsAlpha",
               "QtVsAlpha; #alpha; P_{T}^{ARM}", 100, -1, 1, 300, 0. ,0.3);
   tQAV0->Add(hQtVsAlpha);


   TH1F *hDecayDCAdaughters  = new TH1F( "V0DecayDCAdaughters",
                "DCA between V0 Daughters; DCA (cm); Number of V0",
                400, 0., 2);
   tQAV0->Add(hDecayDCAdaughters);

   TH1F *hDecayCosinePointingAngleXY  = new TH1F("V0CosineOfPointingAngle",
               "Cos of V0 Pointing Angle; Cos (V0 Point.Angl);Number of V0",
               200, 0.99, 1.0);
   tQAV0->Add(hDecayCosinePointingAngleXY);

   TH1F *hfDecayRadXY = new TH1F("V0Radius",
                         "Casc. decay transv. radius; r (cm); Counts" ,
                         1050, 0., 105.0 );
   tQAV0->Add(hfDecayRadXY);
   TH1F *hfDecayProductIPXY  = new TH1F("V0DcaBachToPrimVertex",
               "DCA of Bach. to Prim. Vertex; DCA (cm);Number of V0s",
               250, 0., 2.5);
   tQAV0->Add(hfDecayProductIPXY);
   
   TH2F *fh2TPCdEdxOfPion
    = new TH2F( "TPCdEdxOfPion",
                " TPC dE/dx of Pion; P(GeV/c); TPC signal (ADC) ",
                2000, 0, 10.0, 450, 0., 900.);
   tQAV0->Add(fh2TPCdEdxOfPion);

   TH2F *fh2TPCdEdxOfProton = new TH2F( "TPCdEdxOfProton",
                " TPC dE/dx of Proton; P(GeV/c); TPC signal (ADC) ",
                2000, 0, 10.0, 450, 0., 900.);
   tQAV0->Add(fh2TPCdEdxOfProton);


   TH2F *RadPk
    = new TH2F ("RadiuPK",
                    " Raius test; RadiusforCout; RadiusforPack ",
		                    1000, 0, 100.0, 100, 0., 100.);
   tQAV0->Add(RadPk);

   TH1F *hV0LifeTime  = new TH1F("V0LifeTime",
               "LifeTime; LifeTime;Number of V0s",
               1000, 0., 1000);
   tQAV0->Add(hV0LifeTime);
 
   TH1F *hDecaylength  = new TH1F("DecayLengthV0",
               "DecayLength; DecayLength;Number of V0s",
               200, 0., 20);
   tQAV0->Add(hDecaylength);
   
   TH1F *hDecaylengthLambda  = new TH1F("DecayLengthLambda",
               "DecayLength; DecayLength;Number of V0s",
               400, 0., 40);
   tQAV0->Add(hDecaylengthLambda);
   
   TH1F *hDecaylengthKs0  = new TH1F("DecayLengthKs0",
               "DecayLength; DecayLength;Number of V0s",
               200, 0., 20);
   tQAV0->Add(hDecaylengthKs0);


  fHistList->Add(tQAV0);

}
//
void AliAnalysisTaskFlowQnSPCascade::QACascadeFlow()
{
   TList *tFlowQACascade = new TList();
   tFlowQACascade->SetName("FlowQACascade");
   tFlowQACascade->SetOwner();

   TProfile *fProfXiUxQy[14];
   TProfile *fProfOmegaUxQy[14];
   TProfile *fProfXiUyQx[14];
   TProfile *fProfOmegaUyQx[14];
   TProfile *fProfXiUxQx[14];
   TProfile *fProfOmegaUxQx[14];
   TProfile *fProfXiUyQy[14];
   TProfile *fProfOmegaUyQy[14];

   TString PtBins[]={"0.5_1.0", "1.0_1.5", "1.5_2.0","2.0_2.5", "2.5_3.0", "3.0_3.5",
                      "3.5_4.0","4.0_4.5","4.5_5.0","5.0_6.0","6.0_7.0","7.0_8.0","8.0_9.0","9.0_10.0"};
   TString PhiBins[]={"0_0.1","0.1_0.2","0.2_0.3","0.3_0.4","0.4_0.5","0.5_0.6","0.6_0.7","0.7_0.8",
                     "0.8_0.9","0.9_1.0","1.0_1.1","1.1_1.2","1.2_1.3","1.3_1.4","1.4_1.5","1.5_1.6",
                     "1.6_1.7","1.7_1.8","1.8_1.9","1.9_2.0"};

 for(Int_t i = 0; i < 14; i++){

   
   fProfXiUxQy[i] =new TProfile(Form("hProfXiUxQyPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; u_{x}* Q_{y}",
                                      20 ,1.27 ,1.37);

    tFlowQACascade->Add(fProfXiUxQy[i]);


    fProfOmegaUxQy[i] =new TProfile(Form("hProfOmegaUxQyPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; u_{x}* Q_{y}",
                                      20, 1.62, 1.72);
    tFlowQACascade->Add(fProfOmegaUxQy[i]);

    fProfXiUyQx[i] =new TProfile(Form("hProfXiUyQxPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; u_{y}* Q_{x}",
                                      20 ,1.27 ,1.37);

    tFlowQACascade->Add(fProfXiUyQx[i]);


    fProfOmegaUyQx[i] =new TProfile(Form("hProfOmegaUyQxPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; u_{y}* Q_{x}",
                                      20, 1.62, 1.72);
    tFlowQACascade->Add(fProfOmegaUyQx[i]);

    fProfXiUxQx[i] =new TProfile(Form("hProfXiUxQxPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; u_{x}* Q_{x}",
                                      20 ,1.27 ,1.37);

    tFlowQACascade->Add(fProfXiUxQx[i]);


    fProfOmegaUxQx[i] =new TProfile(Form("hProfOmegaUxQxPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; u_{x}* Q_{x}",
                                      20, 1.62, 1.72);
    tFlowQACascade->Add(fProfOmegaUxQx[i]);


    fProfXiUyQy[i] =new TProfile(Form("hProfXiUyQyPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; u_{y}* Q_{y}",
                                      20 ,1.27 ,1.37);

    tFlowQACascade->Add(fProfXiUyQy[i]);


    fProfOmegaUyQy[i] =new TProfile(Form("hProfOmegaUyQyPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; u_{y}* Q_{y}",
                                      20, 1.62, 1.72);
    tFlowQACascade->Add(fProfOmegaUyQy[i]);

   
  }


   fHistList->Add(tFlowQACascade);

}
//==============================================Cascade Flow Analysis=======================================

void AliAnalysisTaskFlowQnSPCascade::FlowAnaCascade()

{
 
   TList *tFlowAnaCascade = new TList();
   tFlowAnaCascade->SetName("FlowAnalysisCascade");
   tFlowAnaCascade->SetOwner();

   TString PtBins[]={"0.5_1.0", "1.0_1.5", "1.5_2.0","2.0_2.5", "2.5_3.0", "3.0_3.5",
                      "3.5_4.0","4.0_4.5","4.5_5.0","5.0_6.0","6.0_7.0","7.0_8.0","8.0_9.0","9.0_10.0"};
   TString PhiBins[]={"0_0.1","0.1_0.2","0.2_0.3","0.3_0.4","0.4_0.5","0.5_0.6","0.6_0.7","0.7_0.8",
                     "0.8_0.9","0.9_1.0","1.0_1.1","1.1_1.2","1.2_1.3","1.3_1.4","1.4_1.5","1.5_1.6",
                     "1.6_1.7","1.7_1.8","1.8_1.9","1.9_2.0"}; 

  TH2F *fh2MassVsPtXiAll
    = new TH2F("MassVsPtXiAll",
               "M_{#Xi candidates} vs Pt; Pt (GeV/c); M(#Lambda, #pi) (GeV/c^{2})",
               100, 0., 10., 1600, 1.28, 1.37);
   tFlowAnaCascade->Add(fh2MassVsPtXiAll);

   TH2F *fh2MassVsPtOmegaAll
    = new TH2F("MassVsPtOmegaAll",
               "M_{#Omega candidates} vs Pt; Pt (GeV/c); M(#Lambda, K^{+}) (GeV/c^{2})",
               100, 0., 10., 2000, 1.63, 1.72);
   tFlowAnaCascade->Add(fh2MassVsPtOmegaAll);


   TH2F *fh2MassVsEtaXiAll
     = new TH2F("MassVsEtaXiAll",
                    "M_{#Xi candidates} vs #eta;  #eta ; M(#Lambda, #pi) (GeV/c^{2})",
		                   1800, -0.9, 0.9, 1600, 1.28, 1.37);
   tFlowAnaCascade->Add(fh2MassVsEtaXiAll);

   TH2F *fh2MassVsEtaOmegaAll
     = new TH2F("MassVsEtaOmegaAll",
                         "M_{#Xi candidates} vs #eta;  #eta; M(#Lambda, K^{+}) (GeV/c^{2})" ,
			                                    1800, -0.9, 0.9, 1600, 1.63, 1.72);
   tFlowAnaCascade->Add(fh2MassVsEtaOmegaAll);

   TH2F *fh2MassVsPhiXiAll
     = new TH2F("MassVsPhiXiAll",
                    "M_{#Xi candidates} vs #phi;  #phi; M(#Lambda, #pi) (GeV/c^{2})",
                                   1280, -6.4, 6.4,1600, 1.28, 1.37);
   tFlowAnaCascade->Add(fh2MassVsPhiXiAll);

   TH2F *fh2MassVsPhiOmegaAll
     = new TH2F("MassVsPhiOmegaAll",
                         "M_{#Xi candidates} vs #phi;  #Phi ; M(#Lambda, K^{+}) (GeV/c^{2})",
                                                            1280, -6.4, 6.4, 1600, 1.63, 1.72);
   tFlowAnaCascade->Add(fh2MassVsPhiOmegaAll);

  /*TH1F *fh3MassVsPhiXiAll
     = new TH1F("MassVsPtVsPhiXiAll",
                     "M_{#Xi candidates} vs #phi and pt;  M(#Lambda, K^{+}) (GeV/c^{2})",
                       1600, 1.28, 1.37);
   tFlowAnaCascade->Add(fh3MassVsPhiXiAll);
*/

   TH3F *fh3MassVsPhiXiAll
     = new TH3F("MassVsPtVsPhiXiAll",
                     "M_{#Xi candidates} vs #phi and pt;  #Phi; Pt; M(#Lambda, K^{+}) (GeV/c^{2})",
                     128, 0, 6.4, 100, 0., 10., 160, 1.28, 1.37);
   tFlowAnaCascade->Add(fh3MassVsPhiXiAll);


   TH3F *fh3MassVsPhiOmegaAll
     = new TH3F("MassVsPtVsPhiOmegaAll",
                         "M_{#Omega candidates} vs #phi;  #Phi ;M(#Lambda, K^{+}) (GeV/c^{2})",
                     128, 0, 6.4,100, 0., 10., 160, 1.63, 1.72);
   tFlowAnaCascade->Add(fh3MassVsPhiOmegaAll);


//"MassVsPtVsPhiOmegaAll"

   TProfile *fProfXiVnV0A[14];
   TProfile *fProfOmegaVnV0A[14];
   TProfile *fProfXiVnV0C[14];
   TProfile *fProfOmegaVnV0C[14];
   TProfile *fProfXiVnTPC[14];
   TProfile *fProfOmegaVnTPC[14];
   TProfile *fProfXiUn[14];
   TProfile *fProfOmegaUn[14];

 for(Int_t i = 0; i < 14; i++){

    fProfXiVnV0A[i] =new TProfile(Form("hProfXiVnV0APt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tFlowAnaCascade->Add(fProfXiVnV0A[i]);


    fProfOmegaVnV0A[i] =new TProfile(Form("hProfOmegaVnV0APt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tFlowAnaCascade->Add(fProfOmegaVnV0A[i]);

    fProfXiVnV0C[i] =new TProfile(Form("hProfXiVnV0CPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tFlowAnaCascade->Add(fProfXiVnV0C[i]);


    fProfOmegaVnV0C[i] =new TProfile(Form("hProfOmegaVnV0CPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tFlowAnaCascade->Add(fProfOmegaVnV0C[i]);

    fProfXiVnTPC[i] =new TProfile(Form("hProfXiVnTPCPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tFlowAnaCascade->Add(fProfXiVnTPC[i]);


    fProfOmegaVnTPC[i] =new TProfile(Form("hProfOmegaVnTPCPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tFlowAnaCascade->Add(fProfOmegaVnTPC[i]);
 

    fProfXiUn[i] =new TProfile(Form("hProfXiUnPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tFlowAnaCascade->Add(fProfXiUn[i]);
    

    fProfOmegaUn[i] =new TProfile(Form("hProfOmegaUnPt%s", PtBins[i].Data()),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tFlowAnaCascade->Add(fProfOmegaUn[i]);

  }
    fHistList->Add(tFlowAnaCascade);

}

//==============================================V0 Flow Analysis========================================
void AliAnalysisTaskFlowQnSPCascade::FlowAnaV0()
{

    TList *tFlowAnaV0 = new TList();
    tFlowAnaV0->SetName("FlowAnalysisV0");
    tFlowAnaV0->SetOwner();

    TH2F *fh2MassVsPtLambdaAll
    = new TH2F("MassVsPtLambdaAll",
                  "M_{#Lambda candidates} vs Pt; Pt (GeV/c); M() (GeV/c^{2})",
                                 100, 0., 10., 2000, 1.065, 1.165);
    tFlowAnaV0->Add(fh2MassVsPtLambdaAll);
    
  
   TH2F *fh2MassVsEtaLambdaAll
     = new TH2F("MassVsEtaLambdaAll",
                    "M_{#Lambda candidates} vs #eta;  #eta",
                                   1800, -0.9, 0.9, 100, 1.065, 1.165);
   tFlowAnaV0->Add(fh2MassVsEtaLambdaAll);


   TH2F *fh2MassVsPhiLambdaAll
     = new TH2F("MassVsPhiLambdaAll",
                         "M_{#Lambda candidates} vs #Phi;  #Phi",
                                                            1280, -6.4, 6.4, 1600, 1.065, 1.165);
   tFlowAnaV0->Add(fh2MassVsPhiLambdaAll);



    TProfile *fProfLambdaVnV0A[14];
    TProfile *fProfLambdaVnV0C[14];
    TProfile *fProfLambdaVnTPC[14];

    for(Int_t i = 0; i < 14; i++){

    fProfLambdaVnV0A[i] =new TProfile(Form("hProfLambdaVnV0A%d", i),
                                      "; Mass[GeV/c^{2}]; v_{n}",
                                             20, 1.065, 1.165);

    tFlowAnaV0->Add(fProfLambdaVnV0A[i]);
   
    fProfLambdaVnV0C[i] =new TProfile(Form("hProfLambdaVnV0C%d", i),
                                          "; Mass[GeV/c^{2}]; v_{n}",
					     20, 1.065, 1.165);

    tFlowAnaV0->Add(fProfLambdaVnV0C[i]);
    fProfLambdaVnTPC[i] =new TProfile(Form("hProfLambdaVnTPC%d", i),
	                                      "; Mass[GeV/c^{2}]; v_{n}",
					    20, 1.065, 1.165);
    tFlowAnaV0->Add(fProfLambdaVnTPC[i]);


 }


    TH2F *fh2MassVsPtKs0All
    = new TH2F("MassVsPtKs0All",
                  "M_{#K_{S}^{0} candidates} vs Pt; Pt (GeV/c); M() (GeV/c^{2})",
                                 100, 0., 10., 1800, 0.4, 0.58);
    tFlowAnaV0->Add(fh2MassVsPtKs0All);
    
 
   
    TH2F *fh2MassVsEtaKs0All
     = new TH2F("MassVsEtaKs0All",
                    "M_{#Ks0 candidates} vs #eta;  #eta",
                                   1800, -0.9, 0.9, 1800, 0.4, 0.58);
    tFlowAnaV0->Add(fh2MassVsEtaKs0All);
 
 

   TH2F *fh2MassVsPhiKs0All
     = new TH2F("MassVsPhiKs0All",
                         "M_{K_{S}^{0} candidates} vs #Phi;  #Phi",
                                                            1280, -6.4, 6.4, 1800, 0.4, 0.58);
   tFlowAnaV0->Add(fh2MassVsPhiKs0All);



    TProfile *fProfKs0VnV0A[14];
    TProfile *fProfKs0VnV0C[14];
    TProfile *fProfKs0VnTPC[14];

    for(Int_t i = 0; i < 14; i++){

    fProfKs0VnV0A[i] =new TProfile(Form("hProfKs0VnV0A%d", i),
                                      "; Mass[GeV/c^{2}]; v_{n}",
                                             20, 0.4, 0.58);

    tFlowAnaV0->Add(fProfKs0VnV0A[i]);
   
    fProfKs0VnV0C[i] =new TProfile(Form("hProfKs0VnV0C%d", i),
                                          "; Mass[GeV/c^{2}]; v_{n}",
					     20, 0.4, 0.58);

    tFlowAnaV0->Add(fProfKs0VnV0C[i]);
    fProfKs0VnTPC[i] =new TProfile(Form("hProfKs0VnTPC%d", i),
	                                      "; Mass[GeV/c^{2}]; v_{n}",
					    20, 0.4, 0.58);
    tFlowAnaV0->Add(fProfKs0VnTPC[i]);


 }


   fHistList->Add(tFlowAnaV0);

}

//=======================Cascade charge Separation for Flow Analysis ============================

void AliAnalysisTaskFlowQnSPCascade::SepFlowAnaCascade()

{
 //-------------------------------------------Mass------------------------------------------------
   TList *tSepFlowAnaCascade = new TList();
   tSepFlowAnaCascade->SetName("SepFlowAnalysisCascade");
   tSepFlowAnaCascade->SetOwner();

   TH2F *fh2MassVsPtXiPositive
    = new TH2F("MassVsPtXiPositive",
               "M_{#Xi candidates} vs Pt; Pt (GeV/c); M(#Lambda, #pi) (GeV/c^{2})",
               100, 0., 10., 1600, 1.28, 1.37);
   tSepFlowAnaCascade->Add(fh2MassVsPtXiPositive);

   TH2F *fh2MassVsPtXiNegative
    = new TH2F("MassVsPtXiNegative",
                      "M_{#Xi candidates} vs Pt; Pt (GeV/c); M(#Lambda, #pi) (GeV/c^{2})",
		                     100, 0., 10., 1600, 1.28, 1.37);
   tSepFlowAnaCascade->Add(fh2MassVsPtXiNegative);


   TH2F *fh2MassVsPtOmegaPositive
    = new TH2F("MassVsPtOmegaPositive",
               "M_{#Omega candidates} vs Pt; Pt (GeV/c); M(#Lambda, K^{+}) (GeV/c^{2})",
               100, 0., 10., 2000, 1.63, 1.72);
   tSepFlowAnaCascade->Add(fh2MassVsPtOmegaPositive);

   TH2F *fh2MassVsPtOmegaNegative
    = new TH2F("MassVsPtOmegaNegative",
                      "M_{#Omega candidates} vs Pt; Pt (GeV/c); M(#Lambda, K^{+}) (GeV/c^{2})",
		                     100, 0., 10., 2000, 1.63, 1.72);
   tSepFlowAnaCascade->Add(fh2MassVsPtOmegaNegative);
//---------------------------------------------Flow-------------------------------------
   TProfile *fProfXiNegativeVnV0A[14];
   TProfile *fProfOmegaNegativeVnV0A[14];
   TProfile *fProfXiNegativeVnV0C[14];
   TProfile *fProfOmegaNegativeVnV0C[14];
   TProfile *fProfXiNegativeVnTPC[14];
   TProfile *fProfOmegaNegativeVnTPC[14];

   TProfile *fProfXiPositiveVnV0A[14];
   TProfile *fProfOmegaPositiveVnV0A[14];
   TProfile *fProfXiPositiveVnV0C[14];
   TProfile *fProfOmegaPositiveVnV0C[14];
   TProfile *fProfXiPositiveVnTPC[14];
   TProfile *fProfOmegaPositiveVnTPC[14];
 for(Int_t i = 0; i < 14; i++){
//---------------------positive--------
    fProfXiPositiveVnV0A[i] =new TProfile(Form("hProfXiPositiveVnV0A%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tSepFlowAnaCascade->Add(fProfXiPositiveVnV0A[i]);


    fProfOmegaPositiveVnV0A[i] =new TProfile(Form("hProfOmegaPositiveVnV0A%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tSepFlowAnaCascade->Add(fProfOmegaPositiveVnV0A[i]);

    fProfXiPositiveVnV0C[i] =new TProfile(Form("hProfXiPositiveVnV0C%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tSepFlowAnaCascade->Add(fProfXiPositiveVnV0C[i]);


    fProfOmegaPositiveVnV0C[i] =new TProfile(Form("hProfOmegaPositiveVnV0C%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tSepFlowAnaCascade->Add(fProfOmegaPositiveVnV0C[i]);

    fProfXiPositiveVnTPC[i] =new TProfile(Form("hProfXiPositiveVnTPC%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tSepFlowAnaCascade->Add(fProfXiPositiveVnTPC[i]);


    fProfOmegaPositiveVnTPC[i] =new TProfile(Form("hProfOmegaPositiveVnTPC%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tSepFlowAnaCascade->Add(fProfOmegaPositiveVnTPC[i]);
//--------------------negative--------

    fProfXiNegativeVnV0A[i] =new TProfile(Form("hProfXiNegativeVnV0A%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tSepFlowAnaCascade->Add(fProfXiNegativeVnV0A[i]);


    fProfOmegaNegativeVnV0A[i] =new TProfile(Form("hProfOmegaNegativeVnV0A%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tSepFlowAnaCascade->Add(fProfOmegaNegativeVnV0A[i]);

    fProfXiNegativeVnV0C[i] =new TProfile(Form("hProfXiNegativeVnV0C%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tSepFlowAnaCascade->Add(fProfXiNegativeVnV0C[i]);


    fProfOmegaNegativeVnV0C[i] =new TProfile(Form("hProfOmegaNegativeVnV0C%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tSepFlowAnaCascade->Add(fProfOmegaNegativeVnV0C[i]);

    fProfXiNegativeVnTPC[i] =new TProfile(Form("hProfXiNegativeVnTPC%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.27 ,1.37);

    tSepFlowAnaCascade->Add(fProfXiNegativeVnTPC[i]);


    fProfOmegaNegativeVnTPC[i] =new TProfile(Form("hProfOmegaNegativeVnTPC%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20, 1.62, 1.72);
    tSepFlowAnaCascade->Add(fProfOmegaNegativeVnTPC[i]);


 }
    fHistList->Add(tSepFlowAnaCascade);

}

//=======================Lambda and Anti-Lambda Separation for Flow Analysis ============================

void AliAnalysisTaskFlowQnSPCascade::SepFlowAnaV0()

{
 //-------------------------------------------Mass------------------------------------------------
   TList *tSepFlowAnaV0 = new TList();
   tSepFlowAnaV0->SetName("SepFlowAnalysisV0");
   tSepFlowAnaV0->SetOwner();

   TH2F *fh2MassVsPtGeneralV0
    = new TH2F("MassVsPtGeneralV0",
              "M_{#Lambda candidates} vs Pt; Pt (GeV/c); M() (GeV/c^{2})",
                                    100, 0., 10., 2000, 1.065, 1.165);
   tSepFlowAnaV0->Add(fh2MassVsPtGeneralV0);

   TH2F *fh2MassVsPtAntiV0
    = new TH2F("MassVsPtAntiV0",
              "M_{#Lambda candidates} vs Pt; Pt (GeV/c); M() (GeV/c^{2})",
                                    100, 0., 10., 2000, 1.065, 1.165);
   tSepFlowAnaV0->Add(fh2MassVsPtAntiV0);

//---------------------------------------------Flow-------------------------------------
   TProfile *fProfGeneralV0VnV0A[16];
   TProfile *fProfGeneralV0VnV0C[16];
   TProfile *fProfGeneralV0VnTPC[16];

   TProfile *fProfAntiV0VnV0A[16];
   TProfile *fProfAntiV0VnV0C[16];
   TProfile *fProfAntiV0VnTPC[16];
 for(Int_t i = 0; i < 14; i++){
//---------------------Lambda--------
    fProfGeneralV0VnV0A[i] =new TProfile(Form("hProfGeneralV0VnV0A%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.065 ,1.165);

    tSepFlowAnaV0->Add(fProfGeneralV0VnV0A[i]);


    fProfGeneralV0VnV0C[i] =new TProfile(Form("hProfGeneralV0VnV0C%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.065 ,1.165);

    tSepFlowAnaV0->Add(fProfGeneralV0VnV0C[i]);


    fProfGeneralV0VnTPC[i] =new TProfile(Form("hProfGeneralV0VnTPC%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.065 ,1.165);

    tSepFlowAnaV0->Add(fProfGeneralV0VnTPC[i]);


//--------------------Anti-Lambda--------

    fProfAntiV0VnV0A[i] =new TProfile(Form("hProfAntiV0VnV0A%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.065 ,1.165);

    tSepFlowAnaV0->Add(fProfAntiV0VnV0A[i]);


    fProfAntiV0VnV0C[i] =new TProfile(Form("hProfAntiV0VnV0C%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.065 ,1.165);

    tSepFlowAnaV0->Add(fProfAntiV0VnV0C[i]);



    fProfAntiV0VnTPC[i] =new TProfile(Form("hProfAntiV0VnTPC%d", i),
                                     "; Mass[GeV/c^{2}]; v_{n}",
                                      20 ,1.065 ,1.165);

    tSepFlowAnaV0->Add(fProfAntiV0VnTPC[i]);



 }
    fHistList->Add(tSepFlowAnaV0);

}

//===========================================================================================================
void AliAnalysisTaskFlowQnSPCascade::UserCreateOutputObjects()
{
 //----------------------------------FlowQnVectorCorrections----------------------------------------------- 
    AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask =
       dynamic_cast<AliAnalysisTaskFlowVectorCorrections *>(AliAnalysisManager::GetAnalysisManager()->GetTask("FlowQnVectorCorrections"));
    if (flowQnVectorTask != NULL) {
    fFlowQnVectorMgr = flowQnVectorTask->GetAliQnCorrectionsManager();
  }
    else {
    AliFatal("This task needs the Flow Qn vector corrections framework and it is not present. Aborting!!!");
  }

   fFlowQnVectorMgr->GetQnVectorList()->Print("",-1);

    fEventCuts = new AliEventCuts();
//----------------------------------------------------------------------------------------------------------
    fHistList = new TList();
    fHistList->SetOwner();

//   QAEventInput();

   QAEventPlane();
   QATrack();

   TList *tEventCutQA = new TList();
   tEventCutQA->SetName("EventCutQA");
   tEventCutQA->SetOwner();


   fEventCuts->AddQAplotsToList(tEventCutQA); /// fList is your output TList

   fHistList->Add(tEventCutQA);


   if(fIsMultiStrange){
   QACascadeFlow();
   QACascadeCandidates();
   FlowAnaCascade();
   if (fSepAnalysis){
      SepFlowAnaCascade();
    }  
  }
   else if (fIsStrange){
   QAV0Candidates();
   FlowAnaV0();
   if(fSepAnalysis){
      SepFlowAnaV0();
    }

  }	  
 
  PostData(1, fHistList);
 }

//===============================================================================================================


void AliAnalysisTaskFlowQnSPCascade::UserExec(Option_t *) 
{
//----------------------------------read event---------------------------------------------
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler
      = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    if (!fPIDResponse)
        AliError("PIDResponse object was not created"); 
   
    UInt_t fSelectMask= inputHandler->IsEventSelected();
    Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7;
    if(!isINT7selected) return;

    AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    AliAODEvent *fAOD=dynamic_cast<AliAODEvent*>(InputEvent());

    if(fESD){
    
    Info("UserExec", "This task doesn't work with ESD!");
  }
    else if(fAOD){

 AliVEventHandler *inputHandler=man->GetInputEventHandler();
      if (!inputHandler) return;

      TList *uiList = inputHandler->GetUserInfo();
      AliProdInfo prodInfo(uiList);
      prodInfo.List();

      TString fRecoPassName = prodInfo.GetRecoPassName();

      if (fRecoPassName.IsNull()) return;

//-------------------------------------pile-up remove--------------------------------------

      if (!fEventCuts->AcceptEvent(fAOD)) {
	  PostData(1, fHistList);
      return;
}

   float fcent = fEventCuts->GetCentrality(); /// Centrality calculated with the default estimator (V0M for LHC15o)
   const AliVVertex* vtx = fEventCuts->GetPrimaryVertex(); /// Best primary vertex available


   if(fcent<fMinCent||fcent>=fMaxCent) return; //centrality cut

   Int_t nTrk = fAOD->GetNumberOfTracks();
   
    for (Int_t iTrk=0; iTrk<nTrk; iTrk++) {

   AliAODTrack *track = (AliAODTrack*)fAOD->GetTrack(iTrk);
     if (!track) {
      AliError(Form("%s: Could not get Track", GetName()));
      continue;
    }
      Float_t Phi = track->Phi();

      ((TH1F*)((TList*)fHistList->FindObject("Track"))->FindObject("PhiGlobleTrack"))->Fill(Phi);


   if (!track->TestFilterBit(768)) continue;
     // Float_t HybirdsPhi = track->Phi();

     ((TH1F*)((TList*)fHistList->FindObject("Track"))->FindObject("PhiHyBirdsTrack"))->Fill(Phi);

     }
//--------------------------------------------weight------------------

 TString centType[]={"05", "510", "1020","2030", "3040", "4050",
"5060"};
        Int_t cent = 100;
        if (fMinCent == 0){cent = 0;}
        if (fMinCent == 5){cent = 1;}
        if (fMinCent == 10){cent = 2;}
        if (fMinCent == 20){cent = 3;}
        if (fMinCent ==30){cent = 4;}
        if (fMinCent ==40){cent = 5;}
        if (fMinCent ==50){cent = 6;}

 if(fWeightCorrection){
     
     
    TList* weightList = dynamic_cast<TList*>(GetInputData(2));
    if(weightList){
        
    hXiWeight = (TH1F*)(weightList->FindObject("Xi")->FindObject(Form("weightCent%s",centType[cent].Data())));
    hOmegaWeight = (TH1F*)(weightList->FindObject("Omega")->FindObject(Form("weightCent%s",centType[cent].Data())));
        
     if(!hOmegaWeight || !hXiWeight){
    std::cout<<"Efficiency histograms are not available!"<<std::endl;
   }
 }
      if(!weightList){
   std::cout<<"weight histograms are not available!"<<std::endl;
      }

   
for(int s=0;s!= 20;++s){
weightXi[s]= hXiWeight->GetBinContent(s+1);
weightXi[s]= 1/weightXi[s];  

  weightOmega[s]= hOmegaWeight->GetBinContent(s+1);
  weightOmega[s]= 1/weightOmega[s];


  }
}

//------------------------------------EventPlane---------------------------------------------
        Double_t Qn_Ax=0;
        Double_t Qn_Ay=0;
        Double_t Qn_Bx=0;
        Double_t Qn_By=0;
        Double_t Qn_Cx=0;
        Double_t Qn_Cy=0;
        
        Int_t myHarmonic = fHarmonicOrder;
        
        const AliQnCorrectionsQnVector *TPCQnVectorAfterPlain;
        const AliQnCorrectionsQnVector *V0AQnVectorAfterRaw;
        const AliQnCorrectionsQnVector *V0CQnVectorAfterRaw;
        
        const AliQnCorrectionsQnVector *TPCQnVectorAfterRecent;
        const AliQnCorrectionsQnVector *V0AQnVectorAfterPlain;
        const AliQnCorrectionsQnVector *V0CQnVectorAfterPlain;
        const AliQnCorrectionsQnVector *TPCQnVectorAfterAlign;
        const AliQnCorrectionsQnVector *V0AQnVectorAfterRecent;
        const AliQnCorrectionsQnVector *V0CQnVectorAfterRecent;
        const AliQnCorrectionsQnVector *TPCQnVectorAfterTwist;
        const AliQnCorrectionsQnVector *V0AQnVectorAfterAlign;
        const AliQnCorrectionsQnVector *V0CQnVectorAfterAlign;
        const AliQnCorrectionsQnVector *TPCQnVectorLatest;
        const AliQnCorrectionsQnVector *V0AQnVectorLatest;
        const AliQnCorrectionsQnVector *V0CQnVectorLatest;
        const AliQnCorrectionsQnVector *V0AQnVectorAfterRescale;
        const AliQnCorrectionsQnVector *V0CQnVectorAfterRescale;
        const AliQnCorrectionsQnVector *V0AQnVectorAfterTwist;
        const AliQnCorrectionsQnVector *V0CQnVectorAfterTwist;
        
        Double_t psiTPCAfterPlain;
        Double_t psiV0AAfterRaw ;
        Double_t psiV0CAfterRaw ;
        
        
        Double_t psiV0AAfterRescale ;
        Double_t psiV0CAfterRescale ;
        Double_t psiV0AAfterTwist ;
        Double_t psiV0CAfterTwist ;
        
        Double_t psiTPCAfterRecent;
        Double_t psiV0AAfterPlain ;
        Double_t psiV0CAfterPlain ;
        Double_t psiTPCAfterAlign ;
        Double_t psiV0AAfterRecent ;
        Double_t psiV0CAfterRecent ;
        Double_t psiTPCAfterTwist ;
        Double_t psiV0AAfterAlign ;
        Double_t psiV0CAfterAlign ;
        Double_t psiTPCLatest ;
        Double_t psiV0ALatest ;
        Double_t psiV0CLatest ;
        
        //-----------------------------------TPC-----------------------------------
        
        TPCQnVectorAfterPlain = fFlowQnVectorMgr->GetDetectorQnVector("TPC","plain","plain");
        if (TPCQnVectorAfterPlain != NULL){
            psiTPCAfterPlain = TPCQnVectorAfterPlain->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleTPCAfterplain"))->Fill(psiTPCAfterPlain);
            
        }
        
        TPCQnVectorAfterRecent = fFlowQnVectorMgr->GetDetectorQnVector("TPC","rec","rec");
        if (TPCQnVectorAfterRecent != NULL){
            psiTPCAfterRecent = TPCQnVectorAfterRecent->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleTPCAfterRecent"))->Fill(psiTPCAfterRecent);
            
        }
        
        TPCQnVectorAfterTwist = fFlowQnVectorMgr->GetDetectorQnVector("TPC","twist","twist");
        if (TPCQnVectorAfterTwist != NULL){
            psiTPCAfterTwist = TPCQnVectorAfterTwist->EventPlane(myHarmonic);
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleTPCAfterTwist"))->Fill(psiTPCAfterTwist);
            
        }

        TPCQnVectorLatest = fFlowQnVectorMgr->GetDetectorQnVector("TPC","latest","latest");
        if (TPCQnVectorLatest != NULL){
            psiTPCLatest = TPCQnVectorLatest->EventPlane(myHarmonic);
          ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleTPCLatest"))->Fill(psiTPCLatest);
            
        }

        //----------------------------------V0A------------------------------------
        
        V0AQnVectorAfterRescale = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","rescale","rescale");
        if (V0AQnVectorAfterRescale != NULL){
            psiV0AAfterRescale = V0AQnVectorAfterRescale->EventPlane(myHarmonic);
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0AAfterRescale"))->Fill(psiV0AAfterRescale);
            
        }
        
        
        V0AQnVectorAfterRaw = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","raw","raw");
        if (V0AQnVectorAfterRaw != NULL){
            psiV0AAfterRaw = V0AQnVectorAfterRaw->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0AAfterRaw"))->Fill(psiV0AAfterRaw);
            
        }
        
        
        V0AQnVectorAfterPlain = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","plain","plain");
        if (V0AQnVectorAfterPlain != NULL){
            psiV0AAfterPlain = V0AQnVectorAfterPlain->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0AAfterPlain"))->Fill(psiV0AAfterPlain);
            
        }
        V0AQnVectorAfterRecent = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","rec","rec");
        if (V0AQnVectorAfterRecent != NULL){
            psiV0AAfterRecent = V0AQnVectorAfterRecent->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0AAfterRecent"))->Fill(psiV0AAfterRecent);
            
        }
        V0AQnVectorAfterAlign = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","align","align");
        if (V0AQnVectorAfterAlign != NULL){
            psiV0AAfterAlign = V0AQnVectorAfterAlign->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0AAfterAlign"))->Fill(psiV0AAfterAlign);
            
        }
        
        V0AQnVectorAfterTwist = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","twist","twist");
        if (V0AQnVectorAfterTwist != NULL){
            psiV0AAfterTwist = V0AQnVectorAfterTwist->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0AAfterTwist"))->Fill(psiV0AAfterTwist);
            
        }
        
        
        
        
        V0AQnVectorLatest = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","latest","latest");
        if (V0AQnVectorLatest != NULL){
            psiV0ALatest = V0AQnVectorLatest->EventPlane(myHarmonic);
            
           ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0ALatest"))->Fill(psiV0ALatest);
            
        }
        //-----------------------------------V0C----------------------------------------
        
        V0CQnVectorAfterRescale = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","rescale","rescale");
        if (V0CQnVectorAfterRescale != NULL){
            psiV0CAfterRescale = V0CQnVectorAfterRescale->EventPlane(myHarmonic);
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0CAfterRescale"))->Fill(psiV0CAfterRescale);
            
        }
        
        
        
        V0CQnVectorAfterRaw = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","raw","raw");
        if (V0CQnVectorAfterRaw != NULL){
            psiV0CAfterRaw = V0CQnVectorAfterRaw->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0CAfterRaw"))->Fill(psiV0CAfterRaw);
            
        }
        
        V0CQnVectorAfterPlain = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","plain","plain");
        if (V0CQnVectorAfterPlain != NULL){
            psiV0CAfterPlain = V0CQnVectorAfterPlain->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0CAfterPlain"))->Fill(psiV0CAfterPlain);
            
        }
        V0CQnVectorAfterRecent = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","rec","rec");
        if (V0CQnVectorAfterRecent != NULL){
            psiV0CAfterRecent = V0CQnVectorAfterRecent->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0CAfterRecent"))->Fill(psiV0CAfterRecent);
            
        }

        V0CQnVectorAfterAlign = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","align","align");
        if (V0CQnVectorAfterAlign != NULL){
            psiV0CAfterAlign = V0CQnVectorAfterAlign->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0CAfterAlign"))->Fill(psiV0CAfterAlign);
            
        }
        
        V0CQnVectorAfterTwist = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","twist","twist");
        if (V0CQnVectorAfterTwist != NULL){
            psiV0CAfterTwist = V0CQnVectorAfterTwist->EventPlane(myHarmonic);
            
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0CAfterTwist"))->Fill(psiV0CAfterTwist);
            
        }
        
        V0CQnVectorLatest = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","latest","latest");
        if (V0CQnVectorLatest != NULL){
            psiV0CLatest = V0CQnVectorLatest->EventPlane(myHarmonic);
            
           ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0CLatest"))->Fill(psiV0CLatest);
            
        }
        
        //--------------------------------------------------------------------------------
        
        const AliQnCorrectionsQnVector *TPCQnVector;
        const AliQnCorrectionsQnVector *V0AQnVector;
        const AliQnCorrectionsQnVector *V0CQnVector;
        
        Double_t psiTPC;
        Double_t psiV0A;
        Double_t psiV0C;
        
        V0AQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROA","rec","rec");
        if (V0AQnVector != NULL){
            psiV0A = V0AQnVector->EventPlane(myHarmonic);
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0A"))->Fill(psiV0A);
        }
        
        V0CQnVector = fFlowQnVectorMgr->GetDetectorQnVector("VZEROC","rec","rec");
        if (V0CQnVector != NULL){
            psiV0C = V0CQnVector->EventPlane(myHarmonic);
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleV0C"))->Fill(psiV0C);
        }
        
        TPCQnVector = fFlowQnVectorMgr->GetDetectorQnVector("TPC");
        if (TPCQnVector != NULL){
            psiTPC = TPCQnVector->EventPlane(myHarmonic);
            ((TH1F*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hEPangleTPC"))->Fill(psiTPC);
            
        }
        
        fFlowQnVectorMgr->GetQnVectorList()->Print("",-1);
        
        Qn_Ax=TMath::Cos(fHarmonicOrder*psiV0A);
        Qn_Ay=TMath::Sin(fHarmonicOrder*psiV0A);
        Qn_Cx=TMath::Cos(fHarmonicOrder*psiV0C);
        Qn_Cy=TMath::Sin(fHarmonicOrder*psiV0C);
        Qn_Bx=TMath::Cos(fHarmonicOrder*psiTPC);
        Qn_By=TMath::Sin(fHarmonicOrder*psiTPC);
        
        ((TProfile*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hProfQnProduct"))->Fill(1,  Qn_Bx * Qn_Cx + Qn_By * Qn_Cy);
        ((TProfile*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hProfQnProduct"))->Fill(2,  Qn_Ax * Qn_Bx + Qn_Ay * Qn_By);
        ((TProfile*)((TList*)fHistList->FindObject("EventPlane"))->FindObject("hProfQnProduct"))->Fill(3,  Qn_Ax * Qn_Cx + Qn_Ay * Qn_Cy);

   if(fIsMultiStrange){
     ReadFromAODCascade(fAOD, psiV0A, psiV0C, psiTPC, Qn_Ax, Qn_Ay);//用于cascade和omega的分析
 
 }
   else if (fIsStrange){

    ReadFromAODv0(fAOD, psiV0A, psiV0C, psiTPC);//用于lambda的分析
}


 }
   return;
}

//===================================================================================================

//Flow analysis of multi-strange particle xi and omega

//==================================================================================================

void AliAnalysisTaskFlowQnSPCascade::ReadFromAODCascade(AliAODEvent *fAOD, Double_t psiV0A, Double_t psiV0C, Double_t psiTPC,Double_t Qn_Ax, Double_t Qn_Ay){

   Double_t bestPrimaryVtxPos[3] = {-100., -100., -100.};

  Double_t b = fAOD->GetMagneticField();
  int nCascades=fAOD->GetNumberOfCascades();

  const AliAODVertex *primaryBestAODVtx = fAOD->GetPrimaryVertex();
  primaryBestAODVtx->GetXYZ(bestPrimaryVtxPos);
  
  for(Int_t iXi = 0; iXi < nCascades; iXi++){
//cout<<"测试"<<nCascades<<endl;  
    const Int_t nPtB = 14;
    const Int_t nPhiB =20;
    Double_t ptBins[nPtB+1] = {0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7.,8.,9.,10.};
    Double_t phiBins[nPhiB+1] = {0,0.1*TMath::Pi(),0.2*TMath::Pi(),0.3*TMath::Pi(),0.4*TMath::Pi(),
                                 0.5*TMath::Pi(),0.6*TMath::Pi(),0.7*TMath::Pi(),0.8*TMath::Pi(),0.9*TMath::Pi(),
                                 TMath::Pi(),1.1*TMath::Pi(),1.2*TMath::Pi(),1.3*TMath::Pi(),
                                 1.4*TMath::Pi(),1.5*TMath::Pi(),1.6*TMath::Pi(),1.7*TMath::Pi(),
                                 1.8*TMath::Pi(),1.9*TMath::Pi(),2*TMath::Pi()};
    TString PtBins[]={"0.5_1.0", "1.0_1.5", "1.5_2.0","2.0_2.5", "2.5_3.0", "3.0_3.5",
                      "3.5_4.0","4.0_4.5","4.5_5.0","5.0_6.0","6.0_7.0","7.0_8.0","8.0_9.0","9.0_10.0"};
    TString PhiBins[]={"0_0.1","0.1_0.2","0.2_0.3","0.3_0.4","0.4_0.5","0.5_0.6","0.6_0.7","0.7_0.8",
                     "0.8_0.9","0.9_1.0","1.0_1.1","1.1_1.2","1.2_1.3","1.3_1.4","1.4_1.5","1.5_1.6",
                     "1.6_1.7","1.7_1.8","1.8_1.9","1.9_2.0"};   

 
    Double_t dcaXiDaughters = -1.;
    Double_t XiCosOfPointingAngle = -1.;
    Double_t posXi[3] = {-1000., -1000., -1000.};
    Double_t XiRadius = -1000.;
    
    Double_t invMassLambdaAsCascDghter = 0.;
    Double_t V0Chi2Xi = -1.;
    Double_t dcaV0DaughtersXi = -1.;
    
    Double_t dcaBachToPrimaryVtxXi = -1.;
    Double_t dcaV0ToPrimaryVtxXi = -1.;
    Double_t dcaPosToPrimaryVtxXi = -1.;
    Double_t dcaNegToPrimaryVtxXi = -1.;
    Double_t V0CosOfPointingAngleXi = -1.;
    Double_t posV0Xi[3] = {-1000., -1000., -1000.};
    Double_t V0RadiusXi = -1000.;
    
    Double_t invMassXiMinus = 0.;
    Double_t invMassXiPlus = 0.;
    Double_t invMassOmegaMinus = 0.;
    Double_t invMassOmegaPlus = 0.;
    Double_t invMassXiAll =0.;
    Double_t invMassOmegaAll =0.;

    Bool_t isBachelorKaonForTPC = kFALSE;
    Bool_t isBachelorPionForTPC = kFALSE;
    Bool_t isNegPionForTPC = kFALSE;
    Bool_t isPosPionForTPC = kFALSE;
    Bool_t isNegProtonForTPC = kFALSE;
    Bool_t isPosProtonForTPC = kFALSE;
    Double_t XiPse =0.;
    Double_t XiPx = 0., XiPy = 0., XiPz = 0.;
    Double_t XiPt = 0.;
    Double_t XiPtot = 0.;
    Double_t rapXi = -20.;
    Double_t rapOmega = -20.;   
 
    Double_t bachPx = 0., bachPy = 0., bachPz = 0.;
    Double_t bachPt = 0.;
    Double_t bachPtot = 0.;
    
    Double_t V0toXiCosOfPointingAngle = 0.;
    
    Double_t phi = 6.3;
    Double_t alphaXi = -200.;
    Double_t ptArmXi = -200.;

    Double_t distToVtxZBefore = -999.;
    Double_t distToVtxZAfter = -999.;
    Double_t distToVtxXYBefore = -999.;
    Double_t distToVtxXYAfter = -999.;
    Double_t XiPAfter[3] = {-999., -999., -999.};
    Double_t phiAfter = -999.;

    const AliAODcascade *xi = fAOD->GetCascade(iXi);
    if (!xi) continue;
//----------------------------------------track selection----------------------------------------    
    AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(0) );
    AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( xi->GetDaughter(1) );
    AliAODTrack *bTrkXi 
      = dynamic_cast<AliAODTrack*>( xi->GetDecayVertexXi()->GetDaughter(0) );

    if(!pTrkXi || !nTrkXi || !bTrkXi) continue;
    
    UInt_t idxPosXi  = (UInt_t) TMath::Abs( pTrkXi->GetID() );
    UInt_t idxNegXi  = (UInt_t) TMath::Abs( nTrkXi->GetID() );
    UInt_t idxBach   = (UInt_t) TMath::Abs( bTrkXi->GetID() );

    if(idxBach == idxNegXi || idxBach == idxPosXi) continue;

    if( !IsSelected(pTrkXi)|| !IsSelected(nTrkXi)
        || !IsSelected(bTrkXi) ) continue;
//---------------------------------------candidate selection---------------------------------------------------
//----------------------------------------cascade vertex------------------------------------ 
     XiPx = xi->MomXiX();
     XiPy = xi->MomXiY();
     XiPz = xi->MomXiZ();
     XiPt = TMath::Sqrt(XiPx*XiPx + XiPy*XiPy);
     XiPtot= TMath::Sqrt(XiPx*XiPx + XiPy*XiPy + XiPz*XiPz);
     XiPse = 0.5*TMath::Log((XiPtot+XiPz)/(XiPtot-XiPz));

     alphaXi = xi->AlphaXi();
     ptArmXi = xi->PtArmXi();
     rapXi = xi->RapXi();
     rapOmega = xi->RapOmega();

    if(xi->ChargeXi()<0){
       invMassXiMinus = xi->MassXi();
       invMassOmegaMinus = xi->MassOmega();
       invMassXiAll=xi->MassXi();
       invMassOmegaAll=xi->MassOmega();
      }else{
       invMassXiPlus = xi->MassXi();
       invMassOmegaPlus = xi->MassOmega();
       invMassXiAll=xi->MassXi();
       invMassOmegaAll=xi->MassOmega();
     }

    if(xi->ChargeXi() < 0)
      invMassLambdaAsCascDghter = xi->MassLambda();
    else
      invMassLambdaAsCascDghter = xi->MassAntiLambda();//1

    dcaXiDaughters = xi->DcaXiDaughters();//2
    XiCosOfPointingAngle = xi->CosPointingAngleXi(bestPrimaryVtxPos[0],
                                                  bestPrimaryVtxPos[1],
                                                  bestPrimaryVtxPos[2]);//3
    posXi[0] = xi->DecayVertexXiX();
    posXi[1] = xi->DecayVertexXiY();
    posXi[2] = xi->DecayVertexXiZ();
    XiRadius = TMath::Sqrt(posXi[0]*posXi[0]
                           +posXi[1]*posXi[1]);

    posV0Xi[0] = xi->DecayVertexV0X();
    posV0Xi[1] = xi->DecayVertexV0Y();
    posV0Xi[2] = xi->DecayVertexV0Z();
    V0RadiusXi = TMath::Sqrt(posV0Xi[0]*posV0Xi[0]
                             +posV0Xi[1]*posV0Xi[1]);


    Double_t XiLength = TMath::Sqrt((posXi[0]-bestPrimaryVtxPos[0])*(posXi[0]-bestPrimaryVtxPos[0])
                           +(posXi[1]-bestPrimaryVtxPos[1])*(posXi[1]-bestPrimaryVtxPos[1])
                           +(posXi[2]-bestPrimaryVtxPos[2])*(posXi[2]-bestPrimaryVtxPos[2])); 


    dcaBachToPrimaryVtxXi = xi->DcaBachToPrimVertex();//4
    dcaV0ToPrimaryVtxXi = xi->DcaV0ToPrimVertex();//5

//-----------------------------------------V0 vertex------------------------------------------
     dcaV0DaughtersXi = xi->DcaV0Daughters(); //1    
     V0CosOfPointingAngleXi 
       = xi->CosPointingAngle(bestPrimaryVtxPos);//2
     dcaPosToPrimaryVtxXi = xi->DcaPosToPrimVertex();//3
     dcaNegToPrimaryVtxXi = xi->DcaNegToPrimVertex();//3
//-------------------------------------------------------------------------------------------
//cascade
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("DcaXiDaughters"))->Fill(dcaXiDaughters);
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("DcaBachToPrimVertex"))->Fill(dcaBachToPrimaryVtxXi);
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("XiCosineOfPointingAngle"))->Fill(XiCosOfPointingAngle);
//V0  
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("DcaV0DaughtersXi"))->Fill(dcaV0DaughtersXi);
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("V0CosOfPointingAngleXi"))->Fill(V0CosOfPointingAngleXi);
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("V0RadiusXi"))->Fill(V0RadiusXi);
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("DcaV0ToPrimVertexXi"))->Fill(dcaV0ToPrimaryVtxXi);
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("DcaPosToPrimVertexXi"))->Fill(dcaPosToPrimaryVtxXi);
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("DcaNegToPrimVertexXi"))->Fill(dcaNegToPrimaryVtxXi);
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("Armenteros"))->Fill(alphaXi, ptArmXi);

//----------------------------------------candidate selection cut------------------------------
//    if(XiPse < fXiPseMin || XiPse > fXiPseMax ) continue ;
    if(V0RadiusXi < fV0RadiusXiMin || V0RadiusXi > fV0RadiusXiMax ) continue;   
    if(XiRadius < fXiRadiusMin || XiRadius > fXiRadiusMax) continue;
    if(dcaXiDaughters > fdcaXiDaughtersMax) continue;
    if(XiCosOfPointingAngle < fXiCosOfPointingAngleMin) continue;
    if(dcaV0ToPrimaryVtxXi < fdcaV0ToPrimaryVtxXiMin) continue;
    if(dcaBachToPrimaryVtxXi < fdcaBachToPrimaryVtxXiMin) continue;

    if(TMath::Abs(invMassLambdaAsCascDghter-1.11568) > fLambdaMassWind) continue;

    if(dcaV0DaughtersXi > fdcaV0DaughtersXi) continue;
    if(V0CosOfPointingAngleXi <fV0CosOfPointingAngleXiMin) continue;
    if(dcaPosToPrimaryVtxXi < fdcaPosToPrimaryVtxXiMin) continue;
    if(dcaNegToPrimaryVtxXi < fdcaNegToPrimaryVtxXiMin) continue;
   // if(XiPse > 0.8) continue; 

    Bool_t PassXiRapCut=kTRUE;
      if (rapXi < fXiPseMin || rapXi > fXiPseMax) PassXiRapCut=kFALSE;
    Bool_t PassOmegaRapCut=kTRUE;
      if (rapOmega < fXiPseMin || rapOmega > fXiPseMax) PassOmegaRapCut=kFALSE;

//------------------------------------------------------------------------------------------- 
//----------------------------------------cascade QA-----------------------------------------

     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kKaon))<fXiPIDsigma)
       isBachelorKaonForTPC = kTRUE;
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(bTrkXi, AliPID::kPion))<fXiPIDsigma)
       isBachelorPionForTPC = kTRUE;

     //Negative V0 daughter
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kPion))<fXiPIDsigma)
       isNegPionForTPC = kTRUE;
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrkXi, AliPID::kProton))<fXiPIDsigma)
       isNegProtonForTPC = kTRUE;

     //Positive V0 daughter
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kPion))< fXiPIDsigma)
       isPosPionForTPC = kTRUE;
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrkXi, AliPID::kProton))<fXiPIDsigma)
       isPosProtonForTPC = kTRUE;


     if(isPosPionForTPC && pTrkXi->IsOn(AliESDtrack::kTPCin)){
      ((TH2F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("TPCdEdxOfPion"))->Fill(pTrkXi->P(),pTrkXi->GetTPCsignal());
      
    }
     if(isNegPionForTPC && nTrkXi->IsOn(AliESDtrack::kTPCin)){
    
      ((TH2F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("TPCdEdxOfPion"))->Fill(nTrkXi->P(),nTrkXi->GetTPCsignal());
 
    }
     if(isPosProtonForTPC && pTrkXi->IsOn(AliESDtrack::kTPCin)){

      ((TH2F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("TPCdEdxOfProton"))->Fill(pTrkXi->P(),pTrkXi->GetTPCsignal());

    }
     if(isNegProtonForTPC && nTrkXi->IsOn(AliESDtrack::kTPCin)){
      ((TH2F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("TPCdEdxOfProton"))->Fill(nTrkXi->P(),nTrkXi->GetTPCsignal());

    }

    if(isBachelorKaonForTPC && bTrkXi->IsOn(AliESDtrack::kTPCin)){

      ((TH2F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("TPCdEdxOfBachelorKaon"))->Fill(bTrkXi->P(),bTrkXi->GetTPCsignal());

    }
     if(isBachelorPionForTPC && bTrkXi->IsOn(AliESDtrack::kTPCin)){

      ((TH2F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("TPCdEdxOfBachelorPion"))->Fill(bTrkXi->P(),bTrkXi->GetTPCsignal());

    }

     bachPx = xi->MomBachX();
     bachPy = xi->MomBachY();
     bachPz = xi->MomBachZ();
     bachPt = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy);
     bachPtot = TMath::Sqrt(bachPx*bachPx + bachPy*bachPy + bachPz*bachPz);
     V0toXiCosOfPointingAngle = xi->CosPointingAngle( xi->GetDecayVertexXi() );


     distToVtxZBefore = posXi[2]-bestPrimaryVtxPos[2];
     distToVtxXYBefore
       = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
                     *(posXi[0] - bestPrimaryVtxPos[0])
                     +(posXi[1] - bestPrimaryVtxPos[1])
                     *(posXi[1] - bestPrimaryVtxPos[1]));

     XiPAfter[0] = XiPx;
     XiPAfter[1] = XiPy;
     XiPAfter[2] = XiPz;
     Propagate(bestPrimaryVtxPos, posXi, XiPAfter, b, xi->ChargeXi());
     distToVtxZAfter = posXi[2] - bestPrimaryVtxPos[2];
     distToVtxXYAfter = TMath::Sqrt((posXi[0] - bestPrimaryVtxPos[0])
                                    *(posXi[0] - bestPrimaryVtxPos[0])
                                    +(posXi[1] - bestPrimaryVtxPos[1])
                                    *(posXi[1] - bestPrimaryVtxPos[1]));
     phiAfter = TMath::Pi() + TMath::ATan2(-XiPAfter[1],-XiPAfter[0]);
    
     Double_t phiV0A = phiAfter;
              phiV0A -= psiV0A;
     Double_t phiTPC = phiAfter;
              phiTPC -= psiTPC;
     Double_t phiV0C = phiAfter;
	      phiV0C -= psiV0C;
//---------------------------------Flow Analysis------------------------------------------


 if(xi->ChargeXi() < 0){
        if(isPosProtonForTPC
           && isNegPionForTPC){
          if(isBachelorPionForTPC && PassXiRapCut){
            //Xi
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPtXiAll"))->Fill(XiPt, invMassXiMinus);
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsEtaXiAll"))->Fill(XiPse, invMassXiMinus);
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPhiXiAll"))->Fill(phiAfter, invMassXiMinus);

    ((TH3F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPtVsPhiXiAll"))->Fill(XiPt,phiAfter, invMassXiMinus);

           for(int r=0; r!=14; ++r) {
              if(XiPt > ptBins[r]
                 && XiPt < ptBins[r+1]){

         //   for(int s=0;s!= 20;++s){
         //     if(phiAfter >phiBins[s]&& phiAfter <phiBins[s+1]){

    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfXiVnV0APt%s", 
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiV0A)/*,weightXi[s]*/);
    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfXiVnV0CPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiV0C)/*,weightXi[s]*/);
    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfXiVnTPCPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiTPC)/*,weightXi[s]*/);
   
    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfXiUnPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)/*,weightXi[s]*/);

     ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfXiUxQyPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)*Qn_Ay/*,weightXi[s]*/);
     ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfXiUxQxPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)*Qn_Ax/*,weightXi[s]*/);
     ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfXiUyQxPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Sin(fHarmonicOrder*phiAfter)*Qn_Ax/*,weightXi[s]*/);
     ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfXiUyQyPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Sin(fHarmonicOrder*phiAfter)*Qn_Ay/*,weightXi[s]*/);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("XiPt"))->Fill(XiPt);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("PhiXi"))->Fill(phiAfter);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("hXiPesrapidity"))->Fill(XiPse);
              // }
 
            // }
            }
          }
        } 

if(isBachelorKaonForTPC && PassOmegaRapCut){

    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPtOmegaAll"))->Fill(XiPt, invMassOmegaMinus);
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsEtaOmegaAll"))->Fill(XiPse, invMassOmegaMinus);
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPhiOmegaAll"))->Fill(phiAfter, invMassOmegaMinus);



((TH3F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPtVsPhiOmegaAll"))->Fill(XiPt,phiAfter, invMassOmegaMinus);

       for(int r=0; r!=14; ++r) {
              if(XiPt > ptBins[r]
                 && XiPt < ptBins[r+1]){
  

//for(int s=0;s!= 20;++s){
//              if(phiAfter >phiBins[s]&& phiAfter <phiBins[s+1]){
 
  ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfOmegaVnV0APt%s",
                                               PtBins[r].Data())))
                                                ->Fill(invMassOmegaAll, TMath::Cos(fHarmonicOrder*phiV0A)/*,weightOmega[s]*/);
   ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfOmegaVnV0CPt%s",
                                               PtBins[r].Data())))
                                                ->Fill(invMassOmegaAll, TMath::Cos(fHarmonicOrder*phiV0C)/*,weightOmega[s]*/);
   ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfOmegaVnTPCPt%s",
                                               PtBins[r].Data())))
                                                ->Fill(invMassOmegaAll, TMath::Cos(fHarmonicOrder*phiTPC)/*,weightOmega[s]*/);
   ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfOmegaUnPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)/*,weightOmega[s]*/);

   ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfOmegaUxQyPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)*Qn_Ay/*,weightOmega[s]*/);
   ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfOmegaUxQxPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)*Qn_Ax/*,weightOmega[s]*/);
   ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfOmegaUyQxPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Sin(fHarmonicOrder*phiAfter)*Qn_Ax/*,weightOmega[s]*/);
   ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfOmegaUyQyPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Sin(fHarmonicOrder*phiAfter)*Qn_Ay/*,weightOmega[s]*/);

     //Omega
    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("OmegaPt"))->Fill(XiPt);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("PhiOmega"))->Fill(phiAfter);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("hOmegaPesrapidity"))->Fill(XiPse);

                //}
  
               //}           
              }
            }
          }
        
        }
      }//if charge < 0


 if(xi->ChargeXi() > 0){
        if(isNegProtonForTPC
           && isPosPionForTPC){
          if(isBachelorPionForTPC && PassXiRapCut){
            //Xi
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPtXiAll"))->Fill(XiPt, invMassXiPlus);

    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsEtaXiAll"))->Fill(XiPse, invMassXiPlus);
    
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPhiXiAll"))->Fill(phiAfter, invMassXiPlus);


    ((TH3F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPtVsPhiXiAll"))->Fill(XiPt,phiAfter, invMassXiPlus);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("XiPt"))->Fill(XiPt);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("PhiXi"))->Fill(phiAfter);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("hXiPesrapidity"))->Fill(XiPse);

            for(int r=0; r!=14; ++r) {
              if(XiPt > ptBins[r]
                 && XiPt < ptBins[r+1]){
 

          //  for(int s=0;s!= 20;++s){
           //   if(phiAfter >phiBins[s]&& phiAfter <phiBins[s+1]){

   ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfXiVnV0APt%s",
                                               PtBins[r].Data())))
                                  ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiV0A)/*,weightXi[s]*/);
    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfXiVnV0CPt%s",
                                               PtBins[r].Data())))
                                  ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiV0C)/*,weightXi[s]*/);
    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfXiVnTPCPt%s",
                                               PtBins[r].Data())))
                                  ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiTPC)/*,weightXi[s]*/);
    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfXiUnPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)/*,weightXi[s]*/);
  
     ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfXiUxQyPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)*Qn_Ay/*,weightXi[s]*/);
     ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfXiUxQxPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)*Qn_Ax/*,weightXi[s]*/);
     ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfXiUyQxPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Sin(fHarmonicOrder*phiAfter)*Qn_Ax/*,weightXi[s]*/);
     ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfXiUyQyPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Sin(fHarmonicOrder*phiAfter)*Qn_Ay/*,weightXi[s]*/);

                //}
               //}
              }
            }
          }
        if(isBachelorKaonForTPC && PassOmegaRapCut){
            //Omega
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPtOmegaAll"))->Fill(XiPt, invMassOmegaPlus);
             
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsEtaOmegaAll"))->Fill(XiPse, invMassOmegaPlus);
  
    ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPhiOmegaAll"))->Fill(phiAfter, invMassOmegaPlus);

((TH3F*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject("MassVsPtVsPhiOmegaAll"))->Fill(XiPt,phiAfter, invMassOmegaPlus);
    


    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("OmegaPt"))->Fill(XiPt);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("PhiOmega"))->Fill(phiAfter);

    ((TH1F*)((TList*)fHistList->FindObject("Cascade"))->FindObject("hOmegaPesrapidity"))->Fill(XiPse);
for(int r=0; r!=14; ++r) {
              if(XiPt > ptBins[r]
                 && XiPt < ptBins[r+1]){

//for(int s=0;s!= 20;++s){
//              if(phiAfter >phiBins[s]&& phiAfter <phiBins[s+1]){

    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfOmegaVnV0APt%s",
                                               PtBins[r].Data())))
                                  ->Fill(invMassOmegaAll, TMath::Cos(fHarmonicOrder*phiV0A)/*,weightOmega[s]*/);
    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfOmegaVnV0CPt%s",
                                               PtBins[r].Data())))
                                  ->Fill(invMassOmegaAll, TMath::Cos(fHarmonicOrder*phiV0C)/*,weightOmega[s]*/);
    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfOmegaVnTPCPt%s",
                                               PtBins[r].Data())))
                                  ->Fill(invMassOmegaAll, TMath::Cos(fHarmonicOrder*phiTPC)/*,weightOmega[s]*/);
    ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisCascade"))->FindObject(Form("hProfOmegaUnPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)/*,weightOmega[s]*/);
    ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfOmegaUxQyPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)*Qn_Ay/*,weightOmega[s]*/);
   ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfOmegaUxQxPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Cos(fHarmonicOrder*phiAfter)*Qn_Ax/*,weightOmega[s]*/);
   ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfOmegaUyQxPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Sin(fHarmonicOrder*phiAfter)*Qn_Ax/*,weightOmega[s]*/);
   ((TProfile*)((TList*)fHistList->FindObject("FlowQACascade"))->FindObject(Form("hProfOmegaUyQyPt%s",
                                               PtBins[r].Data())))
                                   ->Fill(invMassXiAll, TMath::Sin(fHarmonicOrder*phiAfter)*Qn_Ay/*,weightOmega[s]*/);


  //               }
   //             }
              }
            }
          }
        }
      }// if charge > 0
 
 }//end of cascade loop

  PostData(1, fHistList);

}

//==========================================================================================

//Flow Analysis of strange bayon particle lambda

//==========================================================================================

void AliAnalysisTaskFlowQnSPCascade::ReadFromAODv0(AliAODEvent *fAOD, Double_t psiV0A, Double_t psiV0C, Double_t psiTPC) {
/*
  AliFlowTrackCuts *cutsRPTPC
    = AliFlowTrackCuts::GetAODTrackCutsForFilterBit(1);
  cutsRPTPC->SetPtRange(0.15, 5.);
  cutsRPTPC->SetEtaRange(-0.8, 0.8);
  cutsRPTPC->SetMinNClustersTPC(70); 

  AliFlowTrackCuts *cutsRPTPC 
    = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts();
  cutsRPTPC->SetParamType( AliFlowTrackCuts::kGlobal );
  cutsRPTPC->SetAODfilterBit(1); // for AOD compatibility
  
   AliFlowTrackCuts *cutsRPVZE;
   cutsRPVZE = AliFlowTrackCuts::GetBetaVZEROOnlyTrackCuts();
   
   if(fRPFromTPC){

     cutsRPTPC->SetEvent(fAOD);
     }
     else{

     cutsRPVZE->SetEvent(fAOD);
     }
*/
     Int_t nV0s = fAOD->GetNumberOfV0s();
     for(Int_t iV0 = 0; iV0 < nV0s; iV0++){

         AliAODv0 *V0 = fAOD->GetV0(iV0);
         if(!V0) continue;
         if(!fOnline) if(V0->GetOnFlyStatus() ) continue;
         if(fOnline)  if(!V0->GetOnFlyStatus() ) continue;


     Double_t bestPrimaryVtxPos[3] = {-100., -100., -100.};

     Double_t b = fAOD->GetMagneticField();
//----------------------------------------------------daughter track selection-----------------------------------     
     const AliAODVertex *primaryBestAODVtx = fAOD->GetPrimaryVertex();
     primaryBestAODVtx->GetXYZ(bestPrimaryVtxPos);

     AliAODTrack *V0TrackNeg=dynamic_cast<AliAODTrack *>(V0->GetDaughter(1));
     AliAODTrack *V0TrackPos=dynamic_cast<AliAODTrack *>(V0->GetDaughter(0));


     if( !IsSelected(V0TrackNeg)|| !IsSelected(V0TrackPos) ) continue;

    if (V0TrackNeg->Charge() == V0TrackPos->Charge()) continue;


    Double_t lEta  = V0->PseudoRapV0();
    if(TMath::Abs(lEta) > fPrimaryTrackEta) continue;

//---------------------------------------------------candidate selection--------------------------------
     Double_t fV0Px = V0->MomV0X(); 
     Double_t fV0Py = V0->MomV0Y(); 
     Double_t fV0Pz = V0->MomV0Z(); 
     Double_t fV0Ptot = TMath::Sqrt(fV0Px*fV0Px+fV0Py*fV0Py+fV0Pz*fV0Pz);
     
     Double_t fLambdaMass =0;
//     Bool_t isAntiLambda =kFALSE;
//     Bool_t isGeLambda =kFALSE;
     Bool_t isNegPionForTPC = kFALSE;
     Bool_t isPosPionForTPC = kFALSE;
     Bool_t isNegProtonForTPC = kFALSE;
     Bool_t isPosProtonForTPC = kFALSE;
    
     Double_t posV0[3];
     Double_t fDecayAlpha = V0->AlphaV0(); 

 //    if( fDecayAlpha<0 ) isAntiLambda=kTRUE; 
 //    if( fDecayAlpha>0 ) isGeLambda=kTRUE;
     if(fDecayAlpha>0) fLambdaMass = V0->MassLambda();
     else fLambdaMass = V0->MassAntiLambda();
     Double_t fMassK0s = V0->MassK0Short();
     posV0[0] = V0->DecayVertexV0X();
     posV0[1] = V0->DecayVertexV0Y();
     posV0[2] = V0->DecayVertexV0Z();
     Double_t V0Radius = TMath::Sqrt(posV0[0]*posV0[0]
                           +posV0[1]*posV0[1]);

     Double_t V0Length = TMath::Sqrt((posV0[0]-bestPrimaryVtxPos[0])*(posV0[0]-bestPrimaryVtxPos[0])
                           +(posV0[1]-bestPrimaryVtxPos[1])*(posV0[1]-bestPrimaryVtxPos[1])
                           +(posV0[2]-bestPrimaryVtxPos[2])*(posV0[2]-bestPrimaryVtxPos[2])); 

     Double_t dx = V0->Xv()-bestPrimaryVtxPos[0];
     Double_t dy = V0->Yv()-bestPrimaryVtxPos[1];

     Double_t fDecayRad = TMath::Sqrt( dx*dx + dy*dy );    
     Double_t lPt=TMath::Sqrt(V0->Pt2V0());
     Double_t V0LifeTime =fLambdaMass*fDecayRad/lPt;

     Double_t fV0Rapidity = V0->RapLambda();
     Double_t fV0DCAdaughters = V0->DcaV0Daughters();
     Double_t fV0CosinePointingAngle = V0->CosPointingAngle(bestPrimaryVtxPos);
     Double_t fV0DCAToPrimVertex = V0->DcaV0ToPrimVertex(); 
     Double_t fDecayLengthV0= V0->DecayLengthV0(bestPrimaryVtxPos);

     Double_t LambdaEnergy = TMath::Sqrt( fLambdaMass*fLambdaMass + fV0Ptot*fV0Ptot );
     Double_t LambdaGamma =  LambdaEnergy/fLambdaMass;
     Double_t fDecayLengthLambda = fDecayLengthV0/LambdaGamma;

     Double_t Ks0Energy = TMath::Sqrt( fMassK0s*fMassK0s + fV0Ptot*fV0Ptot );
     Double_t Ks0Gamma =  Ks0Energy/fMassK0s;
     Double_t fDecayLengthKs0 = fDecayLengthV0/Ks0Gamma;


     Double_t fV0Pt = V0->Pt();
     Double_t fV0PhiNormal = V0->Phi();
     Double_t fV0Eta = V0->Eta();
     Double_t dQT= V0->PtArmV0();
     Double_t dALPHA= V0->AlphaV0();

     Double_t fV0Phi = TMath::ATan2(fV0Py,fV0Px);

     ((TH1F*)((TList*)fHistList->FindObject("V0"))->FindObject("DecayAlpha"))->Fill(fDecayAlpha);
//     ((TH2F*)((TList*)fHistList->FindObject("V0"))->FindObject("QtVsAlpha"))->Fill(dALPHA,dQT);
     ((TH2F*)((TList*)fHistList->FindObject("V0"))->FindObject("PhiV0com"))->Fill(fV0Phi,fV0PhiNormal);


     ((TH1F*)((TList*)fHistList->FindObject("V0"))->FindObject("V0DecayDCAdaughters"))->Fill(fV0DCAdaughters);
     ((TH1F*)((TList*)fHistList->FindObject("V0"))->FindObject("V0CosineOfPointingAngle"))->Fill(fV0CosinePointingAngle);
     ((TH1F*)((TList*)fHistList->FindObject("V0"))->FindObject("V0Radius"))->Fill(fDecayRad);
     ((TH1F*)((TList*)fHistList->FindObject("V0"))->FindObject("V0DcaBachToPrimVertex"))->Fill(fV0DCAToPrimVertex);
     ((TH2F*)((TList*)fHistList->FindObject("V0"))->FindObject("RadiuPK"))->Fill(fDecayRad,fDecayLengthV0);
     ((TH1F*)((TList*)fHistList->FindObject("V0"))->FindObject("V0LifeTime"))->Fill(V0LifeTime);
     ((TH1F*)((TList*)fHistList->FindObject("V0"))->FindObject("DecayLengthV0"))->Fill(fDecayLengthV0);
     ((TH1F*)((TList*)fHistList->FindObject("V0"))->FindObject("DecayLengthLambda"))->Fill(fDecayLengthLambda);
     ((TH1F*)((TList*)fHistList->FindObject("V0"))->FindObject("DecayLengthKs0"))->Fill(fDecayLengthKs0);



//---------------------------------------candidate selection cut------------------------
   // if(TMath::Abs(fV0Eta) > tV0Eta) continue;
   // if(fV0Pt < tV0Pt) continue;
   // if(dQT<0.2*TMath::Abs(dALPHA)) continue;
    if(fV0DCAdaughters > tV0DCAdaughtersMax) continue;
    if(fV0CosinePointingAngle < tV0CosinePointingAngleMin) continue;
    if(fDecayRad <= tDecayRad) continue;
    if(TMath::Abs(fV0Rapidity) > tDecayRapidity) continue;
    if(fV0DCAToPrimVertex < tV0DCAToPrimVertexMin) continue;

    Bool_t PassLambdaDLCut=kTRUE;
      if (fDecayLengthLambda >tDecayLengthV0*7.89) PassLambdaDLCut=kFALSE;
    Bool_t PassKs0DLCut=kTRUE;
      if ((fDecayLengthKs0 >tDecayLengthV0*2.6844)||(dQT<0.2*TMath::Abs(dALPHA))) PassKs0DLCut=kFALSE;

   ((TH2F*)((TList*)fHistList->FindObject("V0"))->FindObject("QtVsAlpha"))->Fill(dALPHA,dQT);
//----------------------------------------------------V0 QA---------------------------------------------

     //Negative V0 daughter
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(V0TrackNeg, AliPID::kPion))<fV0PIDsigma)
       isNegPionForTPC = kTRUE;
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(V0TrackNeg, AliPID::kProton))<fV0PIDsigma)
       isNegProtonForTPC = kTRUE;

     //Positive V0 daughter
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(V0TrackPos, AliPID::kPion))< fV0PIDsigma)
       isPosPionForTPC = kTRUE;
     if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(V0TrackPos, AliPID::kProton))<fV0PIDsigma)
       isPosProtonForTPC = kTRUE;

    if(PassLambdaDLCut||PassKs0DLCut){
     if(isNegPionForTPC && V0TrackNeg->IsOn(AliESDtrack::kTPCin)){
     ((TH2F*)((TList*)fHistList->FindObject("V0"))->FindObject("TPCdEdxOfPion"))->Fill(V0TrackNeg->P(),V0TrackNeg->GetTPCsignal());
     }

     if(isNegProtonForTPC && V0TrackNeg->IsOn(AliESDtrack::kTPCin)){
     ((TH2F*)((TList*)fHistList->FindObject("V0"))->FindObject("TPCdEdxOfProton"))->Fill(V0TrackNeg->P(),V0TrackNeg->GetTPCsignal());
     }

     if(isPosPionForTPC && V0TrackPos->IsOn(AliESDtrack::kTPCin)){
     ((TH2F*)((TList*)fHistList->FindObject("V0"))->FindObject("TPCdEdxOfPion"))->Fill(V0TrackPos->P(),V0TrackPos->GetTPCsignal());
     }

     if(isPosProtonForTPC && V0TrackPos->IsOn(AliESDtrack::kTPCin)){
     ((TH2F*)((TList*)fHistList->FindObject("V0"))->FindObject("TPCdEdxOfProton"))->Fill(V0TrackPos->P(),V0TrackPos->GetTPCsignal());
     }

   }

     Double_t phiV0A = fV0Phi;
                   phiV0A -= psiV0A;
     Double_t phiV0C = fV0Phi;
                   phiV0C -= psiV0C;
     Double_t phiTPC = fV0Phi;
                   phiTPC -= psiTPC;

     const Int_t nPtB = 14;
     Double_t ptBins[nPtB+1] = {0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7., 8., 9., 10};
          for(int r=0; r!=14; ++r) {
           if(fV0Pt > ptBins[r] && fV0Pt < ptBins[r+1]){
   
        if((isNegPionForTPC
           && isPosProtonForTPC && PassLambdaDLCut)||(isPosPionForTPC
           && isNegProtonForTPC && PassLambdaDLCut)){



     ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject("MassVsPtLambdaAll"))->Fill(fV0Pt, fLambdaMass);
     
      ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject("MassVsEtaLambdaAll"))->Fill(fV0Eta, fLambdaMass);

     ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject("MassVsPhiLambdaAll"))->Fill(fV0Phi, fLambdaMass);

     ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject(Form("hProfLambdaVnV0A%d", r)))
                               ->Fill(fLambdaMass, TMath::Cos(fHarmonicOrder*phiV0A));
     ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject(Form("hProfLambdaVnV0C%d", r)))
                               ->Fill(fLambdaMass, TMath::Cos(fHarmonicOrder*phiV0C));
     ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject(Form("hProfLambdaVnTPC%d", r)))
                               ->Fill(fLambdaMass, TMath::Cos(fHarmonicOrder*phiTPC));


    }
    
     if (isNegPionForTPC && isPosPionForTPC && PassKs0DLCut){

     ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject("MassVsPtKs0All"))->Fill(fV0Pt, fMassK0s);

     ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject("MassVsEtaKs0All"))->Fill(fV0Eta, fMassK0s);

     ((TH2F*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject("MassVsPhiKs0All"))->Fill(fV0Phi, fMassK0s);

     ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject(Form("hProfKs0VnV0A%d", r)))
                               ->Fill(fMassK0s, TMath::Cos(fHarmonicOrder*phiV0A));
     ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject(Form("hProfKs0VnV0C%d", r)))
                               ->Fill(fMassK0s, TMath::Cos(fHarmonicOrder*phiV0C));
     ((TProfile*)((TList*)fHistList->FindObject("FlowAnalysisV0"))->FindObject(Form("hProfKs0VnTPC%d", r)))
                               ->Fill(fMassK0s, TMath::Cos(fHarmonicOrder*phiTPC));

     } 


if (fSepAnalysis){

    if(isNegPionForTPC
           && isPosProtonForTPC && PassLambdaDLCut){

    ((TH2F*)((TList*)fHistList->FindObject("SepFlowAnalysisV0"))->FindObject("MassVsPtGeneralV0"))->Fill(fV0Pt, fLambdaMass);

   ((TProfile*)((TList*)fHistList->FindObject("SepFlowAnalysisV0"))->FindObject(Form("hProfGeneralV0VnV0A%d", r)))
                                  ->Fill(fLambdaMass, TMath::Cos(fHarmonicOrder*phiV0A));

    ((TProfile*)((TList*)fHistList->FindObject("SepFlowAnalysisV0"))->FindObject(Form("hProfGeneralV0VnV0C%d", r)))
                                  ->Fill(fLambdaMass, TMath::Cos(fHarmonicOrder*phiV0C));
    ((TProfile*)((TList*)fHistList->FindObject("SepFlowAnalysisV0"))->FindObject(Form("hProfGeneralV0VnTPC%d", r)))
                                  ->Fill(fLambdaMass, TMath::Cos(fHarmonicOrder*phiTPC));

          }

        if(isPosPionForTPC
           && isNegProtonForTPC && PassLambdaDLCut){

    ((TH2F*)((TList*)fHistList->FindObject("SepFlowAnalysisV0"))->FindObject("MassVsPtAntiV0"))->Fill(fV0Pt, fLambdaMass);
    ((TProfile*)((TList*)fHistList->FindObject("SepFlowAnalysisV0"))->FindObject(Form("hProfAntiV0VnV0A%d", r)))
                                  ->Fill(fLambdaMass, TMath::Cos(fHarmonicOrder*phiV0A));
    ((TProfile*)((TList*)fHistList->FindObject("SepFlowAnalysisV0"))->FindObject(Form("hProfAntiV0VnV0C%d", r)))
                                  ->Fill(fLambdaMass, TMath::Cos(fHarmonicOrder*phiV0C));
    ((TProfile*)((TList*)fHistList->FindObject("SepFlowAnalysisV0"))->FindObject(Form("hProfAntiV0VnTPC%d", r)))
                                  ->Fill(fLambdaMass, TMath::Cos(fHarmonicOrder*phiTPC));

           }
        }
    }

  }

   PostData(1, fHistList);

  }
  
}

//=======================================================================

Bool_t AliAnalysisTaskFlowQnSPCascade::IsSelected(AliAODTrack *t){
   // Pseudorapidity cut                                
   if (TMath::Abs(t->Eta()) > fTrackEta) return kFALSE;
   //pt cut                        
   if(t->Pt() < fTrackPtMin) return kFALSE;
   // TPC refit                   
   if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
   // Minimum number of clusters         
   Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1);
   if (nCrossedRowsTPC < fTPCNcls) return kFALSE;
   Int_t findable=t->GetTPCNclsF();
   if (findable <= 0) return kFALSE;
   if (nCrossedRowsTPC/findable < 0.8) return kFALSE;

   if (t->Chi2perNDF()> 4.0)return kFALSE;
   return kTRUE;
}



void AliAnalysisTaskFlowQnSPCascade::Terminate(Option_t *)
{

}

void AliAnalysisTaskFlowQnSPCascade::Propagate(Double_t vv[3],
                                           Double_t x[3],
                                           Double_t p[3],
                                           Double_t bz,
                                           Short_t sign){
  //Propagation to the primary vertex to determine the px and py
  //x, p are the position and momentum as input and output
  //bz is the magnetic field along z direction
  //sign is the charge of particle for propagation

  Double_t pp = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  Double_t len = (vv[2]-x[2])*pp/p[2];
  Double_t a = -kB2C*bz*sign;

  Double_t rho = a/pp;
  x[0] += p[0]*TMath::Sin(rho*len)/a - p[1]*(1-TMath::Cos(rho*len))/a;
  x[1] += p[1]*TMath::Sin(rho*len)/a + p[0]*(1-TMath::Cos(rho*len))/a;
  x[2] += p[2]*len/pp;

  Double_t p0=p[0];
  p[0] = p0  *TMath::Cos(rho*len) - p[1]*TMath::Sin(rho*len);
  p[1] = p[1]*TMath::Cos(rho*len) + p0  *TMath::Sin(rho*len);
}

