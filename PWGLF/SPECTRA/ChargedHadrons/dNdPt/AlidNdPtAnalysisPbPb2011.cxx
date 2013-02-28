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
//------------------------------------------------------------------------------
// AlidNdPtAnalysisPbPb2011 class. 
// 
// a. functionality:
// - fills analysis control histograms
// - fills generic correction matrices 
// - generates correction matrices 
//
// b. data members:
// - generic correction matrices
// - control histograms
//
// Author: J.Otwinowski 04/11/2008 
// last change: 2011-04-04 by M.Knichel
//------------------------------------------------------------------------------

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "THnSparse.h"

#include "AliHeader.h"  
#include "AliGenEventHeader.h"  
#include "AliInputEventHandler.h"  
#include "AliAnalysisManager.h"  
#include "AliStack.h"  
#include "AliESDEvent.h"  
#include "AliMCEvent.h"  
#include "AliESDtrackCuts.h"  
#include "AliLog.h" 
#include "AliMultiplicity.h"
#include "AliTracker.h"

#include "AliCentrality.h"

#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"
#include "AlidNdPtAnalysisPbPb2011.h"

#include "TSystem.h"
#include "TROOT.h"

using namespace std;

ClassImp(AlidNdPtAnalysisPbPb2011)

//_____________________________________________________________________________
AlidNdPtAnalysisPbPb2011::AlidNdPtAnalysisPbPb2011(): AlidNdPt(),
fAnalysisFolder(0),
fHistogramsOn(kFALSE),

// rec. track pt vs true track pt correlation matrix 
fTrackPtCorrelationMatrix(0),

// event level correction
fGenEventMatrix(0),
fTriggerEventMatrix(0),
fRecEventMatrix(0),

//
// track-event level correction 
//
fGenTrackEventMatrix(0),

fTriggerTrackEventMatrix(0),

fRecTrackEventMatrix(0),

// track rec. efficiency correction (fRecPrimTrackMatrix / fGenPrimTrackMatrix)
fGenTrackMatrix(0),
fGenPrimTrackMatrix(0),
fRecPrimTrackMatrix(0),

// secondary track contamination correction (fRecSecTrackMatrix / fRecTrackMatrix)
fRecTrackMatrix(0),
fRecSecTrackMatrix(0),

// multiple rec. track contamination corrections (fRecMultTrackMatrix / fRecTrackMatrix)
fRecMultTrackMatrix(0),

// event control histograms
fMCEventHist1(0),
fRecEventHist1(0),
fRecEventHist2(0),
fRecMCEventHist1(0),
fRecMCEventHist2(0),

// rec. pt and eta resolution w.r.t MC
fRecMCTrackHist1(0),

//multple reconstructed tracks
fMCMultRecTrackHist1(0), 

// rec. track control histograms
fRecTrackHist3(0),

fTriggerAnalysis(0),

fCentralityEstimator(0),

fMultNbins(0),
fPtNbins(0),
fPtCorrNbins(0),
fEtaNbins(0),
fZvNbins(0),
fCentralityNbins(0),
fBinsMult(0),
fBinsPt(0),
fBinsPtCorr(0),
fBinsEta(0),
fBinsZv(0),
fBinsCentrality(0),

fIsInit(kFALSE)
{
  // default constructor
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    fMCTrackHist1[i]=0;     
    fMCPrimTrackHist1[i]=0;     
    fMCPrimTrackHist2[i]=0;     
    fMCSecTrackHist1[i]=0;     
    fRecTrackHist1[i]=0;     
    fRecTrackHist2[i]=0;     
    fRecTrackMultHist1[i]=0;     
  }
  
  fMultNbins = 0;
  fPtNbins = 0;
  fPtCorrNbins = 0;
  fEtaNbins = 0;
  fZvNbins = 0;
  fCentralityNbins = 0;
  fBinsMult = 0;
  fBinsPt = 0;
  fBinsPtCorr = 0;
  fBinsEta = 0;
  fBinsEta = 0;
  fBinsZv = 0;
  fBinsCentrality = 0;
  
  
  //Init();
  SetCentralityEstimator();
}

//_____________________________________________________________________________
AlidNdPtAnalysisPbPb2011::AlidNdPtAnalysisPbPb2011(Char_t* name, Char_t* title): AlidNdPt(name,title),
fAnalysisFolder(0),
fHistogramsOn(kFALSE),

// rec. track pt vs true track pt correlation matrix 
fTrackPtCorrelationMatrix(0),

// event level correction
fGenEventMatrix(0),
fTriggerEventMatrix(0),
fRecEventMatrix(0),

//
// track-event level correction 
//
fGenTrackEventMatrix(0),

fTriggerTrackEventMatrix(0),

fRecTrackEventMatrix(0),

// track rec. efficiency correction (fRecPrimTrackMatrix / fGenPrimTrackMatrix)
fGenTrackMatrix(0),
fGenPrimTrackMatrix(0),
fRecPrimTrackMatrix(0),

// secondary track contamination correction (fRecTrackMatrix - fRecSecTrackMatrix)
fRecTrackMatrix(0),
fRecSecTrackMatrix(0),

// multiple rec. track contamination corrections (fRecTrackMatrix - fRecMultTrackMatrix)
fRecMultTrackMatrix(0),




// event control histograms
fMCEventHist1(0),
fRecEventHist1(0),
fRecEventHist2(0),
fRecMCEventHist1(0),
fRecMCEventHist2(0),

// rec. pt and eta resolution w.r.t MC
fRecMCTrackHist1(0),

//multiple reconstructed tracks
fMCMultRecTrackHist1(0), 

// rec. track control histograms
fRecTrackHist3(0),

fTriggerAnalysis(0),

fCentralityEstimator(0),

fMultNbins(0),
fPtNbins(0),
fPtCorrNbins(0),
fEtaNbins(0),
fZvNbins(0),
fCentralityNbins(0),
fBinsMult(0),
fBinsPt(0),
fBinsPtCorr(0),
fBinsEta(0),
fBinsZv(0),
fBinsCentrality(0),

fIsInit(kFALSE)  

{
  //
  // constructor
  //
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    fMCTrackHist1[i]=0;     
    fMCPrimTrackHist1[i]=0;     
    fMCPrimTrackHist2[i]=0;     
    fMCSecTrackHist1[i]=0;     
    fRecTrackHist1[i]=0;     
    fRecTrackHist2[i]=0;     
    fRecTrackMultHist1[i]=0; 
  }
  
  fMultNbins = 0;
  fPtNbins = 0;
  fPtCorrNbins = 0;
  fEtaNbins = 0;
  fZvNbins = 0;
  fCentralityNbins = 0;
  fBinsMult = 0;
  fBinsPt = 0;
  fBinsPtCorr = 0;
  fBinsEta = 0;
  fBinsEta = 0;
  fBinsZv = 0;
  fBinsCentrality = 0;
  
  // Init();
  SetCentralityEstimator();
}

//_____________________________________________________________________________
AlidNdPtAnalysisPbPb2011::~AlidNdPtAnalysisPbPb2011() {
  //
  // destructor
  //
  if(fTrackPtCorrelationMatrix) delete fTrackPtCorrelationMatrix; fTrackPtCorrelationMatrix=0;
  //
  if(fGenEventMatrix) delete fGenEventMatrix; fGenEventMatrix=0;
  if(fTriggerEventMatrix) delete fTriggerEventMatrix; fTriggerEventMatrix=0;
  if(fRecEventMatrix) delete fRecEventMatrix; fRecEventMatrix=0;
  
  //
  if(fGenTrackEventMatrix) delete fGenTrackEventMatrix; fGenTrackEventMatrix=0;
  if(fTriggerEventMatrix) delete fTriggerEventMatrix; fTriggerEventMatrix=0;
  if(fRecTrackEventMatrix) delete fRecTrackEventMatrix; fRecTrackEventMatrix=0;
  
  //
  if(fGenTrackMatrix) delete fGenTrackMatrix; fGenTrackMatrix=0;
  if(fGenPrimTrackMatrix) delete fGenPrimTrackMatrix; fGenPrimTrackMatrix=0;
  if(fRecPrimTrackMatrix) delete fRecPrimTrackMatrix; fRecPrimTrackMatrix=0;
  //
  if(fRecTrackMatrix) delete fRecTrackMatrix; fRecTrackMatrix=0;
  if(fRecSecTrackMatrix) delete fRecSecTrackMatrix; fRecSecTrackMatrix=0;
  // 
  if(fRecMultTrackMatrix) delete fRecMultTrackMatrix; fRecMultTrackMatrix=0;
  //
  // Control histograms
  //
  if(fMCEventHist1) delete fMCEventHist1; fMCEventHist1=0;
  if(fRecEventHist1) delete fRecEventHist1; fRecEventHist1=0;
  if(fRecEventHist2) delete fRecEventHist2; fRecEventHist2=0;
  if(fRecMCEventHist1) delete fRecMCEventHist1; fRecMCEventHist1=0;
  if(fRecMCEventHist2) delete fRecMCEventHist2; fRecMCEventHist2=0;
  
  //
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    if(fMCTrackHist1[i]) delete fMCTrackHist1[i]; fMCTrackHist1[i]=0;
    if(fMCPrimTrackHist1[i]) delete fMCPrimTrackHist1[i]; fMCPrimTrackHist1[i]=0;
    if(fMCPrimTrackHist2[i]) delete fMCPrimTrackHist2[i]; fMCPrimTrackHist2[i]=0;
    if(fMCSecTrackHist1[i]) delete fMCSecTrackHist1[i]; fMCSecTrackHist1[i]=0;
    if(fRecTrackHist1[i]) delete fRecTrackHist1[i]; fRecTrackHist1[i]=0;
    if(fRecTrackHist2[i]) delete fRecTrackHist2[i]; fRecTrackHist2[i]=0;
    if(fRecTrackMultHist1[i]) delete fRecTrackMultHist1[i]; fRecTrackMultHist1[i]=0;
  }
  if(fRecMCTrackHist1) delete fRecMCTrackHist1; fRecMCTrackHist1=0;
  if(fMCMultRecTrackHist1) delete fMCMultRecTrackHist1; fMCMultRecTrackHist1=0; 
  if(fRecTrackHist3) delete fRecTrackHist3; fRecTrackHist3=0; 
  //
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
  
  if (fTriggerAnalysis) delete fTriggerAnalysis;  fTriggerAnalysis = 0;
  
  if (fBinsMult) delete[] fBinsMult; fBinsMult=0;
  if (fBinsPt) delete[] fBinsPt; fBinsPt=0;
  if (fBinsPtCorr) delete[] fBinsPtCorr; fBinsPtCorr=0;
  if (fBinsEta) delete[] fBinsEta; fBinsEta=0;
  if (fBinsZv) delete[] fBinsZv; fBinsZv=0;
  if (fBinsCentrality) delete[] fBinsCentrality; fBinsCentrality=0;  
}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb2011::Init() {
  
  //define default binning
  Double_t binsMultDefault[48] = {-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,19.5, 20.5, 30.5, 40.5 , 50.5 , 60.5 , 70.5 , 80.5 , 90.5 , 100.5,200.5, 300.5, 400.5, 500.5, 600.5, 700.5, 800.5, 900.5, 1000.5, 2000.5, 3000.5, 4000.5, 5000.5, 6000.5, 7000.5, 8000.5, 9000.5, 10000.5 }; 
  Double_t binsPtDefault[82] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0};
  Double_t binsPtCorrDefault[37] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,3.0,4.0,200.0};    
  Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
  Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
  Double_t binsCentralityDefault[12] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};  
  
  // if no binning is set, use the default
  if (!fBinsMult)			{ SetBinsMult(48,binsMultDefault); }
  if (!fBinsPt)				{ SetBinsPt(82,binsPtDefault); }
  if (!fBinsPtCorr)			{ SetBinsPtCorr(37,binsPtCorrDefault); }
  if (!fBinsEta)				{ SetBinsEta(31,binsEtaDefault); }
  if (!fBinsZv)				{ SetBinsZv(13,binsZvDefault); }  
  if (!fBinsCentrality)	{ SetBinsCentrality(12,binsCentralityDefault); }
  
  for(Int_t i = 0; i < fPtNbins; i++)				{ AliDebug(AliLog::kDebug, Form("ptBin %i: %.2f", i, fBinsPt[i])); }
  for(Int_t i = 0; i < fPtCorrNbins; i++)		{ AliDebug(AliLog::kDebug, Form("ptCorrBin %i: %.2f", i, fBinsPtCorr[i])); }
  for(Int_t i = 0; i < fEtaNbins; i++)			{ AliDebug(AliLog::kDebug, Form("etaBin %i: %.2f", i, fBinsEta[i])); }
  for(Int_t i = 0; i < fZvNbins; i++)				{ AliDebug(AliLog::kDebug, Form("zvBin %i: %.2f", i, fBinsZv[i])); }
  for(Int_t i = 0; i < fMultNbins; i++)			{ AliDebug(AliLog::kDebug, Form("multBin %i: %.2f", i, fBinsMult[i])); }
  for(Int_t i = 0; i < fCentralityNbins; i++)	{ AliDebug(AliLog::kDebug, Form("centBin %i: %.2f", i, fBinsCentrality[i])); }
  
  Int_t binsTrackEventCorrMatrix[4]={fZvNbins-1,fPtNbins-1,fEtaNbins-1,fCentralityNbins-1};
  Int_t binsTrackEvent[4]={fZvNbins-1,fPtNbins-1,fEtaNbins-1,fCentralityNbins-1};
  Int_t binsTrackPtCorrelationMatrix[4]={fPtNbins-1,fPtNbins-1,fEtaNbins-1,fCentralityNbins-1};
  
  fTrackPtCorrelationMatrix = new THnSparseF("fTrackPtCorrelationMatrix","Pt:mcPt:mcEta:Centrality",4,binsTrackPtCorrelationMatrix);
  fTrackPtCorrelationMatrix->SetBinEdges(0,fBinsPt);
  fTrackPtCorrelationMatrix->SetBinEdges(1,fBinsPt);
  fTrackPtCorrelationMatrix->SetBinEdges(2,fBinsEta);
  fTrackPtCorrelationMatrix->SetBinEdges(3,fBinsCentrality);
  fTrackPtCorrelationMatrix->GetAxis(0)->SetTitle("Pt (GeV/c)");
  fTrackPtCorrelationMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fTrackPtCorrelationMatrix->GetAxis(2)->SetTitle("mcEta");
  fTrackPtCorrelationMatrix->GetAxis(3)->SetTitle("Centrality");
  fTrackPtCorrelationMatrix->Sumw2();
  
  //
  // Efficiency and contamination correction matrices
  //
  Int_t binsEventMatrix[3]={fZvNbins-1,fMultNbins-1,fCentralityNbins-1};
  //   Double_t minEventMatrix[3]={-30.,-0.5,0.}; 
  //   Double_t maxEventMatrix[3]={30.,10000.5,100.}; 
  
  fGenEventMatrix = new THnSparseF("fGenEventMatrix","mcZv:multMB:Centrality",3,binsEventMatrix); //,minEventMatrix,maxEventMatrix);
  fGenEventMatrix->SetBinEdges(0,fBinsZv);
  fGenEventMatrix->SetBinEdges(1,fBinsMult);
  fGenEventMatrix->SetBinEdges(2,fBinsCentrality);
  fGenEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fGenEventMatrix->GetAxis(2)->SetTitle("Centrality");
  fGenEventMatrix->Sumw2();
  //
  fTriggerEventMatrix = new THnSparseF("fTriggerEventMatrix","mcZv:multMB:Centrality",3,binsEventMatrix); //,minEventMatrix,maxEventMatrix);
  fTriggerEventMatrix->SetBinEdges(0,fBinsZv);
  fTriggerEventMatrix->SetBinEdges(1,fBinsMult);
  fTriggerEventMatrix->SetBinEdges(2,fBinsCentrality);
  fTriggerEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fTriggerEventMatrix->GetAxis(2)->SetTitle("Centrality");
  fTriggerEventMatrix->Sumw2();
  //
  fRecEventMatrix = new THnSparseF("fRecEventMatrix","mcZv:multMB:Centrality",3,binsEventMatrix); //,minEventMatrix,maxEventMatrix);
  fRecEventMatrix->SetBinEdges(0,fBinsZv);
  fRecEventMatrix->SetBinEdges(1,fBinsMult);
  fRecEventMatrix->SetBinEdges(2,fBinsCentrality);    
  fRecEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fRecEventMatrix->GetAxis(2)->SetTitle("Centrality");
  fRecEventMatrix->Sumw2();
  
  // 
  // track to event corrections
  //
  fGenTrackEventMatrix = new THnSparseF("fGenTrackEventMatrix","mcZv:mcPt:mcEta:Centrality",4,binsTrackEventCorrMatrix);
  fGenTrackEventMatrix->SetBinEdges(0,fBinsZv);
  fGenTrackEventMatrix->SetBinEdges(1,fBinsPtCorr);
  fGenTrackEventMatrix->SetBinEdges(2,fBinsEta);
  fGenTrackEventMatrix->SetBinEdges(3,fBinsCentrality);
  fGenTrackEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackEventMatrix->GetAxis(3)->SetTitle("Centrality");
  fGenTrackEventMatrix->Sumw2();
  //
  fTriggerTrackEventMatrix = new THnSparseF("fTriggerTrackEventMatrix","mcZv:mcPt:mcEta:Centrality",4,binsTrackEventCorrMatrix);
  fTriggerTrackEventMatrix->SetBinEdges(0,fBinsZv);
  fTriggerTrackEventMatrix->SetBinEdges(1,fBinsPtCorr);
  fTriggerTrackEventMatrix->SetBinEdges(2,fBinsEta);
  fTriggerTrackEventMatrix->SetBinEdges(3,fBinsCentrality);
  fTriggerTrackEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerTrackEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fTriggerTrackEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fTriggerTrackEventMatrix->GetAxis(3)->SetTitle("Centrality");
  fTriggerTrackEventMatrix->Sumw2();
  //
  fRecTrackEventMatrix = new THnSparseF("fRecTrackEventMatrix","mcZv:mcPt:mcEta:Centrality",4,binsTrackEventCorrMatrix);
  fRecTrackEventMatrix->SetBinEdges(0,fBinsZv);
  fRecTrackEventMatrix->SetBinEdges(1,fBinsPtCorr);
  fRecTrackEventMatrix->SetBinEdges(2,fBinsEta);
  fRecTrackEventMatrix->SetBinEdges(3,fBinsCentrality);
  fRecTrackEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecTrackEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecTrackEventMatrix->GetAxis(3)->SetTitle("Centrality");
  fRecTrackEventMatrix->Sumw2();
  
  //
  // tracks correction matrices
  //
  fGenTrackMatrix = new THnSparseF("fGenTrackMatrix","mcZv:mcPt:mcEta:Centrality",4,binsTrackEventCorrMatrix);
  fGenTrackMatrix->SetBinEdges(0,fBinsZv);
  fGenTrackMatrix->SetBinEdges(1,fBinsPtCorr);
  fGenTrackMatrix->SetBinEdges(2,fBinsEta);
  fGenTrackMatrix->SetBinEdges(3,fBinsCentrality);
  fGenTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackMatrix->GetAxis(3)->SetTitle("Centrality");
  fGenTrackMatrix->Sumw2();
  
  fGenPrimTrackMatrix = new THnSparseF("fGenPrimTrackMatrix","mcZv:mcPt:mcEta:Centrality",4,binsTrackEventCorrMatrix);
  fGenPrimTrackMatrix->SetBinEdges(0,fBinsZv);
  fGenPrimTrackMatrix->SetBinEdges(1,fBinsPtCorr);
  fGenPrimTrackMatrix->SetBinEdges(2,fBinsEta);
  fGenPrimTrackMatrix->SetBinEdges(3,fBinsCentrality);
  fGenPrimTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenPrimTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenPrimTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenPrimTrackMatrix->GetAxis(3)->SetTitle("Centrality");
  fGenPrimTrackMatrix->Sumw2();
  
  
  fRecPrimTrackMatrix = new THnSparseF("fRecPrimTrackMatrix","mcZv:mcPt:mcEta:Centrality",4,binsTrackEventCorrMatrix);
  fRecPrimTrackMatrix->SetBinEdges(0,fBinsZv);
  fRecPrimTrackMatrix->SetBinEdges(1,fBinsPtCorr);
  fRecPrimTrackMatrix->SetBinEdges(2,fBinsEta);
  fRecPrimTrackMatrix->SetBinEdges(3,fBinsCentrality);
  fRecPrimTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecPrimTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecPrimTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecPrimTrackMatrix->GetAxis(3)->SetTitle("Centrality");
  fRecPrimTrackMatrix->Sumw2();
  
  //
  fRecTrackMatrix = new THnSparseF("fRecTrackMatrix","mcZv:mcPt:mcEta:Centrality",4,binsTrackEventCorrMatrix);
  fRecTrackMatrix->SetBinEdges(0,fBinsZv);
  fRecTrackMatrix->SetBinEdges(1,fBinsPtCorr);
  fRecTrackMatrix->SetBinEdges(2,fBinsEta);
  fRecTrackMatrix->SetBinEdges(3,fBinsCentrality);
  fRecTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecTrackMatrix->GetAxis(3)->SetTitle("Centrality");
  fRecTrackMatrix->Sumw2();
  
  fRecSecTrackMatrix = new THnSparseF("fRecSecTrackMatrix","mcZv:mcPt:mcEta:Centrality",4,binsTrackEventCorrMatrix);
  fRecSecTrackMatrix->SetBinEdges(0,fBinsZv);
  fRecSecTrackMatrix->SetBinEdges(1,fBinsPtCorr);
  fRecSecTrackMatrix->SetBinEdges(2,fBinsEta);
  fRecSecTrackMatrix->SetBinEdges(3,fBinsCentrality);
  fRecSecTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecSecTrackMatrix->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fRecSecTrackMatrix->GetAxis(2)->SetTitle("Eta");
  fRecSecTrackMatrix->GetAxis(3)->SetTitle("Centrality");
  fRecSecTrackMatrix->Sumw2();
  
  //
  fRecMultTrackMatrix = new THnSparseF("fRecMultTrackMatrix","mcZv:mcPt:mcEta:Centrality",4,binsTrackEventCorrMatrix);
  fRecMultTrackMatrix->SetBinEdges(0,fBinsZv);
  fRecMultTrackMatrix->SetBinEdges(1,fBinsPtCorr);
  fRecMultTrackMatrix->SetBinEdges(2,fBinsEta);
  fRecMultTrackMatrix->SetBinEdges(3,fBinsCentrality);
  fRecMultTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecMultTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecMultTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecMultTrackMatrix->GetAxis(3)->SetTitle("Centrality");
  fRecMultTrackMatrix->Sumw2();
  
  //
  // Control analysis histograms
  //
  Int_t binsMCEventHist1[4]={100,100,fZvNbins-1,fCentralityNbins-1};
  Double_t minMCEventHist1[4]={-0.1,-0.1,-30.,0.}; 
  Double_t maxMCEventHist1[4]={0.1,0.1,30.,100.}; 
  fMCEventHist1 = new THnSparseF("fMCEventHist1","mcXv:mcYv:mcZv:Centrality",4,binsMCEventHist1,minMCEventHist1,maxMCEventHist1);
  fMCEventHist1->SetBinEdges(2,fBinsZv);
  fMCEventHist1->SetBinEdges(3,fBinsCentrality);
  fMCEventHist1->GetAxis(0)->SetTitle("mcXv (cm)");
  fMCEventHist1->GetAxis(1)->SetTitle("mcYv (cm)");
  fMCEventHist1->GetAxis(2)->SetTitle("mcZv (cm)");
  fMCEventHist1->GetAxis(3)->SetTitle("Centrality");  
  fMCEventHist1->Sumw2();
  
  //
  Int_t binsRecEventHist1[4]={100,100,fZvNbins-1,fCentralityNbins-1};
  Double_t minRecEventHist1[4]={-3.,-3.,-30.,0.}; 
  Double_t maxRecEventHist1[4]={3.,3.,30.,100.}; 
  
  fRecEventHist1 = new THnSparseF("fRecEventHist1","Xv:Yv:Zv:Centrality",4,binsRecEventHist1,minRecEventHist1,maxRecEventHist1);
  fRecEventHist1->SetBinEdges(2,fBinsZv);
  fRecEventHist1->SetBinEdges(3,fBinsCentrality);
  fRecEventHist1->GetAxis(0)->SetTitle("Xv (cm)");
  fRecEventHist1->GetAxis(1)->SetTitle("Yv (cm)");
  fRecEventHist1->GetAxis(2)->SetTitle("Zv (cm)");
  fRecEventHist1->GetAxis(3)->SetTitle("Centrality");  
  fRecEventHist1->Sumw2();
  
  //
  Int_t binsRecEventHist2[3]={fZvNbins-1, fMultNbins-1, fCentralityNbins-1};
  Double_t minRecEventHist2[3]={-30., -0.5, 0.}; 
  Double_t maxRecEventHist2[3]={30., 10000.5, 100.}; 
  
  fRecEventHist2 = new THnSparseF("fRecEventHist2","Zv:multMB:Centrality",3,binsRecEventHist2,minRecEventHist2,maxRecEventHist2);
  fRecEventHist2->SetBinEdges(0,fBinsZv);
  fRecEventHist2->SetBinEdges(1,fBinsMult);
  fRecEventHist2->SetBinEdges(2,fBinsCentrality);
  fRecEventHist2->GetAxis(0)->SetTitle("Zv (cm)");
  fRecEventHist2->GetAxis(1)->SetTitle("multiplicity MB");
  fRecEventHist2->GetAxis(2)->SetTitle("Centrality"); 
  fRecEventHist2->Sumw2();
  
  //
  Double_t kFact = 0.1;
  Int_t binsRecMCEventHist1[4]={100,100,100, fCentralityNbins-1};
  Double_t minRecMCEventHist1[4]={-10.0*kFact,-10.0*kFact,-10.0*kFact, 0.}; 
  Double_t maxRecMCEventHist1[4]={10.0*kFact,10.0*kFact,10.0*kFact, 100.}; 
  
  fRecMCEventHist1 = new THnSparseF("fRecMCEventHist1","Xv-mcXv:Yv-mcYv:Zv-mcZv:Centrality",4,binsRecMCEventHist1,minRecMCEventHist1,maxRecMCEventHist1);
  fRecMCEventHist1->SetBinEdges(3,fBinsCentrality);
  fRecMCEventHist1->GetAxis(0)->SetTitle("Xv-mcXv (cm)");
  fRecMCEventHist1->GetAxis(1)->SetTitle("Yv-mcYv (cm)");
  fRecMCEventHist1->GetAxis(2)->SetTitle("Zv-mcZv (cm)");
  fRecMCEventHist1->GetAxis(3)->SetTitle("Centrality"); 
  
  fRecMCEventHist1->Sumw2();
  
  //
  Int_t binsRecMCEventHist2[4]={100,100,fMultNbins-1, fCentralityNbins-1};
  Double_t minRecMCEventHist2[4]={-10.0*kFact,-10.0*kFact,-0.5, 0.}; 
  Double_t maxRecMCEventHist2[4]={10.0*kFact,10.0*kFact,10000.5, 100.}; 
  
  fRecMCEventHist2 = new THnSparseF("fRecMCEventHist2","Xv-mcXv:Zv-mcZv:mult:Centrality",4,binsRecMCEventHist2,minRecMCEventHist2,maxRecMCEventHist2);
  fRecMCEventHist2->SetBinEdges(2,fBinsMult);
  fRecMCEventHist2->SetBinEdges(3,fBinsCentrality);
  fRecMCEventHist2->GetAxis(0)->SetTitle("Xv-mcXv (cm)");
  fRecMCEventHist2->GetAxis(1)->SetTitle("Zv-mcZv (cm)");
  fRecMCEventHist2->GetAxis(2)->SetTitle("multiplicity");
  fRecMCEventHist2->GetAxis(3)->SetTitle("Centrality"); 
  fRecMCEventHist2->Sumw2();
  //
  char name[256];
  char title[256];
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) 
  {
    // THnSparse track histograms
    
    Int_t binsMCTrackHist1[4]=  {fPtNbins-1, fEtaNbins-1, 90, fCentralityNbins-1};
    Double_t minMCTrackHist1[4]={0.,-1.5,0., 0.}; 
    Double_t maxMCTrackHist1[4]={50,1.5,2.*TMath::Pi(), 100.}; 
    snprintf(name,256,"fMCTrackHist1_%d",i);
    snprintf(title,256,"mcPt:mcEta:mcPhi:Centrality");
    
    fMCTrackHist1[i] = new THnSparseF(name,title,4,binsMCTrackHist1,minMCTrackHist1,maxMCTrackHist1);
    fMCTrackHist1[i]->SetBinEdges(0,fBinsPt);
    fMCTrackHist1[i]->SetBinEdges(1,fBinsEta);
    fMCTrackHist1[i]->SetBinEdges(3,fBinsCentrality);
    fMCTrackHist1[i]->GetAxis(0)->SetTitle("mcPt (GeV/c)");
    fMCTrackHist1[i]->GetAxis(1)->SetTitle("mcEta");
    fMCTrackHist1[i]->GetAxis(2)->SetTitle("mcPhi (rad)");
    fMCTrackHist1[i]->GetAxis(3)->SetTitle("Centrality"); 
    fMCTrackHist1[i]->Sumw2();
    
    Int_t binsMCPrimTrackHist1[6]=  {fPtNbins-1,fEtaNbins-1,6,20,4000, fCentralityNbins-1};
    Double_t minMCPrimTrackHist1[6]={0.,-1.5,0.,0.,0., 0.}; 
    Double_t maxMCPrimTrackHist1[6]={50.,1.5,6.,20.,4000., 100.}; 
    snprintf(name,256,"fMCPrimTrackHist1_%d",i);
    snprintf(title,256,"mcPt:mcEta:pid:mech:mother:Centrality");
    
    fMCPrimTrackHist1[i] = new THnSparseF(name,title,6,binsMCPrimTrackHist1,minMCPrimTrackHist1,maxMCPrimTrackHist1);
    fMCPrimTrackHist1[i]->SetBinEdges(0,fBinsPt);
    fMCPrimTrackHist1[i]->SetBinEdges(1,fBinsEta);
    fMCPrimTrackHist1[i]->SetBinEdges(5,fBinsCentrality);
    fMCPrimTrackHist1[i]->GetAxis(0)->SetTitle("mcPt (GeV/c)");
    fMCPrimTrackHist1[i]->GetAxis(1)->SetTitle("mcEta");
    fMCPrimTrackHist1[i]->GetAxis(2)->SetTitle("pid");
    fMCPrimTrackHist1[i]->GetAxis(3)->SetTitle("mech");
    fMCPrimTrackHist1[i]->GetAxis(4)->SetTitle("mother");
    fMCPrimTrackHist1[i]->GetAxis(5)->SetTitle("Centrality");
    fMCPrimTrackHist1[i]->Sumw2();
    
    Int_t binsMCPrimTrackHist2[4]=  {4000,20,4000,fCentralityNbins-1};
    Double_t minMCPrimTrackHist2[4]={0.,0.,0., 0.}; 
    Double_t maxMCPrimTrackHist2[4]={4000.,20.,4000., 100.}; 
    snprintf(name,256,"fMCPrimTrackHist2_%d",i);
    snprintf(title,256,"pdg:mech:mother:Centrality");
    
    fMCPrimTrackHist2[i] = new THnSparseF(name,title,4,binsMCPrimTrackHist2,minMCPrimTrackHist2,maxMCPrimTrackHist2);
    fMCPrimTrackHist2[i]->SetBinEdges(3,fBinsCentrality);
    fMCPrimTrackHist2[i]->GetAxis(0)->SetTitle("pdg");
    fMCPrimTrackHist2[i]->GetAxis(1)->SetTitle("mech");
    fMCPrimTrackHist2[i]->GetAxis(2)->SetTitle("mother");
    fMCPrimTrackHist2[i]->GetAxis(3)->SetTitle("Centrality");
    fMCPrimTrackHist2[i]->Sumw2();
    
    Int_t binsMCSecTrackHist1[6]=  {fPtNbins-1,fEtaNbins-1,6,20,4000, fCentralityNbins-1};
    Double_t minMCSecTrackHist1[6]={0.,-1.5,0.,0.,0., 0.}; 
    Double_t maxMCSecTrackHist1[6]={50.,1.5,6.,20.,4000., 100.}; 
    snprintf(name,256,"fMCSecTrackHist1_%d",i);
    snprintf(title,256,"mcPt:mcEta:pid:mech:mother:Centrality");
    
    fMCSecTrackHist1[i] = new THnSparseF(name,title,6,binsMCSecTrackHist1,minMCSecTrackHist1,maxMCSecTrackHist1);
    fMCSecTrackHist1[i]->SetBinEdges(0,fBinsPt);
    fMCSecTrackHist1[i]->SetBinEdges(1,fBinsEta);
    fMCSecTrackHist1[i]->SetBinEdges(5,fBinsCentrality);
    fMCSecTrackHist1[i]->GetAxis(0)->SetTitle("mcPt (GeV/c)");
    fMCSecTrackHist1[i]->GetAxis(1)->SetTitle("mcEta");
    fMCSecTrackHist1[i]->GetAxis(2)->SetTitle("pid");
    fMCSecTrackHist1[i]->GetAxis(3)->SetTitle("mech");
    fMCSecTrackHist1[i]->GetAxis(4)->SetTitle("mother");
    fMCSecTrackHist1[i]->GetAxis(5)->SetTitle("Centrality");
    fMCSecTrackHist1[i]->Sumw2();
    
    Int_t binsRecTrackHist1[4]={fPtNbins-1,fEtaNbins-1,90, fCentralityNbins-1};
    Double_t minRecTrackHist1[4]={0.,-1.5,0., 0.}; 
    Double_t maxRecTrackHist1[4]={50.,1.5,2.*TMath::Pi(), 100.};
    snprintf(name,256,"fRecTrackHist1_%d",i);
    snprintf(title,256,"Pt:Eta:Phi:Centrality");
    fRecTrackHist1[i] = new THnSparseF(name,title,4,binsRecTrackHist1,minRecTrackHist1,maxRecTrackHist1);
    fRecTrackHist1[i]->SetBinEdges(0,fBinsPt);
    fRecTrackHist1[i]->SetBinEdges(1,fBinsEta);
    fRecTrackHist1[i]->SetBinEdges(3,fBinsCentrality);
    fRecTrackHist1[i]->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fRecTrackHist1[i]->GetAxis(1)->SetTitle("#eta");
    fRecTrackHist1[i]->GetAxis(2)->SetTitle("#phi (rad)");
    fRecTrackHist1[i]->GetAxis(3)->SetTitle("Centrality");
    fRecTrackHist1[i]->Sumw2();
    
    snprintf(name,256,"fRecTrackHist2_%d",i);
    snprintf(title,256,"Zv:Pt:Eta:Centrality");
    fRecTrackHist2[i] = new THnSparseF(name,title,4,binsTrackEvent);
    fRecTrackHist2[i]->SetBinEdges(0,fBinsZv);
    fRecTrackHist2[i]->SetBinEdges(1,fBinsPt);
    fRecTrackHist2[i]->SetBinEdges(2,fBinsEta);
    fRecTrackHist2[i]->SetBinEdges(3,fBinsCentrality); 
    fRecTrackHist2[i]->GetAxis(0)->SetTitle("Zv (cm)");
    fRecTrackHist2[i]->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
    fRecTrackHist2[i]->GetAxis(2)->SetTitle("#eta");
    fRecTrackHist2[i]->GetAxis(3)->SetTitle("Centrality");
    fRecTrackHist2[i]->Sumw2();
    
    // 
    Int_t binsRecTrackMultHist1[3]={fPtNbins-1,fMultNbins-1, fCentralityNbins-1};
    Double_t minRecTrackMultHist1[3]={0.,-0.5, -0.}; 
    Double_t maxRecTrackMultHist1[3]={50.,10000.5, 100.};
    snprintf(name,256,"fRecTrackMultHist_%d",i);
    snprintf(title,256,"Pt:Mult:Centrality");
    fRecTrackMultHist1[i] = new THnSparseF(name,title,3,binsRecTrackMultHist1,minRecTrackMultHist1,maxRecTrackMultHist1);
    fRecTrackMultHist1[i]->SetBinEdges(0,fBinsPt);
    fRecTrackMultHist1[i]->SetBinEdges(1,fBinsMult);
    fRecTrackMultHist1[i]->SetBinEdges(2,fBinsCentrality);   
    fRecTrackMultHist1[i]->GetAxis(0)->SetTitle("Pt (GeV/c)");
    fRecTrackMultHist1[i]->GetAxis(1)->SetTitle("multiplicity");
    fRecTrackMultHist1[i]->GetAxis(2)->SetTitle("Centrality");  
    fRecTrackMultHist1[i]->Sumw2();
  }
  
  Int_t binsRecMCTrackHist1[5] = {fPtNbins-1,fEtaNbins-1,100,100, fCentralityNbins-1};
  Double_t minRecMCTrackHist1[5]={0.,-1.5,-0.5,-0.5, 0.}; 
  Double_t maxRecMCTrackHist1[5]={50.,1.5,0.5,0.5, 100.}; 
  
  snprintf(name,256,"fRecMCTrackHist1");
  snprintf(title,256,"mcPt:mcEta:(Pt-mcPt)/mcPt:(Eta-mcEta):Centrality");
  fRecMCTrackHist1 = new THnSparseF(name,title,5,binsRecMCTrackHist1,minRecMCTrackHist1,maxRecMCTrackHist1);
  fRecMCTrackHist1->SetBinEdges(0,fBinsPt);
  fRecMCTrackHist1->SetBinEdges(1,fBinsEta);
  fRecMCTrackHist1->SetBinEdges(4,fBinsCentrality);  
  fRecMCTrackHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fRecMCTrackHist1->GetAxis(1)->SetTitle("mcEta");
  fRecMCTrackHist1->GetAxis(2)->SetTitle("(Pt-mcPt)/mcPt");
  fRecMCTrackHist1->GetAxis(3)->SetTitle("Eta-mcEta");
  fRecMCTrackHist1->GetAxis(4)->SetTitle("Centrality"); 
  
  Int_t binsMCMultRecTrackHist1[4] = {fPtNbins-1,fEtaNbins-1,6, fCentralityNbins-1};
  Double_t minMCMultRecTrackHist1[4]={0.,-1.5,0., 0.}; 
  Double_t maxMCMultRecTrackHist1[4]={50.,1.5,6., 100.}; 
  snprintf(name,256,"fMCMultRecTrackHist1");
  snprintf(title,256,"mcPt:mcEta:pid:Centrality");
  fMCMultRecTrackHist1 = new THnSparseF(name,title,4,binsMCMultRecTrackHist1,minMCMultRecTrackHist1,maxMCMultRecTrackHist1);
  fMCMultRecTrackHist1->SetBinEdges(0,fBinsPt);
  fMCMultRecTrackHist1->SetBinEdges(1,fBinsEta);
  fMCMultRecTrackHist1->SetBinEdges(3,fBinsCentrality);  
  fMCMultRecTrackHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCMultRecTrackHist1->GetAxis(1)->SetTitle("mcEta");
  fMCMultRecTrackHist1->GetAxis(2)->SetTitle("pid");
  fMCMultRecTrackHist1->GetAxis(3)->SetTitle("Centrality"); 
  
  //nClust:chi2PerClust:pt:eta:phi:Centrality
  Int_t binsRecTrackHist3[6]={160,100,fPtNbins-1,fEtaNbins-1,90, fCentralityNbins-1};
  Double_t minRecTrackHist3[6]={0., 0., 0., -1.5, 0., 0.};
  Double_t maxRecRecTrackHist3[6]={160.,10., 50., 1.5, 2.*TMath::Pi(), 100.};
  
  fRecTrackHist3 = new THnSparseF("fRecTrackHist3","nClust:chi2PerClust:pt:eta:phi:Centrality",6,binsRecTrackHist3,minRecTrackHist3,maxRecRecTrackHist3);
  fRecTrackHist3->SetBinEdges(2,fBinsPt);
  fRecTrackHist3->SetBinEdges(3,fBinsEta);
  fRecTrackHist3->SetBinEdges(5,fBinsCentrality);  
  fRecTrackHist3->GetAxis(0)->SetTitle("nClust");
  fRecTrackHist3->GetAxis(1)->SetTitle("chi2PerClust");
  fRecTrackHist3->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fRecTrackHist3->GetAxis(3)->SetTitle("#eta");
  fRecTrackHist3->GetAxis(4)->SetTitle("#phi (rad)");
  fRecTrackHist3->GetAxis(5)->SetTitle("Centrality");
  fRecTrackHist3->Sumw2();
  
  
  // init folder
  fAnalysisFolder = CreateFolder("folderdNdPt","Analysis dNdPt Folder");
  
  // init trigger analysis (for zdc cut)
  fTriggerAnalysis = new AliTriggerAnalysis;
  
  // set init flag
  fIsInit = kTRUE;
}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb2011::Process(AliESDEvent *const esdEvent, AliMCEvent *const mcEvent)
{
  //  init if not done already
  if (!fIsInit) { Init(); }
  //
  // Process real and/or simulated events
  //
  if(!esdEvent) {
    AliDebug(AliLog::kError, "esdEvent not available");
    return;
  }
  
  // get selection cuts
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AlidNdPtAcceptanceCuts *recCuts = GetRecAcceptanceCuts();   
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 
  
  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }
  if (0 == recCuts) { recCuts = accCuts;}
  
  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;
  
  // 
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler)
  {
    Printf("ERROR: Could not receive input handler");
    return;
  }
  
  if(evtCuts->IsTriggerRequired())  
  {
    // always MB
    //isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;
    isEventTriggered = inputHandler->IsEventSelected() & GetTriggerMask();
    
    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;
    //SetPhysicsTriggerSelection(physicsSelection);
    
    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }
  
  
  // centrality determination
  Float_t centralityF = -2.;
  AliCentrality *esdCentrality = esdEvent->GetCentrality();
  centralityF = esdCentrality->GetCentralityPercentile((const char*)GetCentralityEstimator());
  
  // if centrality percentile is larger then the highest bin, return
  if(centralityF > fBinsCentrality[fCentralityNbins-1]) 
  {
    AliDebug(AliLog::kError, Form("Event out of centrality binning (event: %.2f, maxallowed: %.2f)",centralityF, fBinsCentrality[fCentralityNbins-1]));
    return;
  }
  
  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);
  
  Int_t multMCTrueTracks = 0;
  if(IsUseMCInfo())
  {
    //
    if(!mcEvent) {
      AliDebug(AliLog::kError, "mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }
    
    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);
    
    Double_t vMCEventHist1[4]={vtxMC[0],vtxMC[1],vtxMC[2],centralityF};
    fMCEventHist1->Fill(vMCEventHist1);
    
    // multipliticy of all MC primary tracks
    // in Zv, pt and eta ranges)
    multMCTrueTracks = AlidNdPtHelper::GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);
    
  } // end bUseMC
  
  // get reconstructed vertex  
  const AliESDVertex* vtxESD = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
    vtxESD = esdEvent->GetPrimaryVertexTPC();
  }
  else if(GetAnalysisMode() == AlidNdPtHelper::kTPCITS) {
    vtxESD = esdEvent->GetPrimaryVertexTracks();
  }
  else {
    return;
  }
  
  if(!vtxESD) return;
  
  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD); 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());
  
  // vertex contributors
  Int_t multMBTracks = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC || GetAnalysisMode() == AlidNdPtHelper::kTPCITS) 
  {
    if(vtxESD->GetStatus()) {
      multMBTracks = vtxESD->GetNContributors();
    }
  } 
  else {
    AliDebug(AliLog::kError, Form("Found analysis type %d", GetAnalysisMode()));
    return; 
  }
  
  TObjArray *allChargedTracks=0;
  //Int_t multAll=0, multAcc=0, multRec=0;
  Int_t multAll=0, multRec=0;
  Int_t *labelsAll=0, *labelsAcc=0, *labelsRec=0;
  
  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    
    // get all charged tracks
    allChargedTracks = AlidNdPtHelper::GetAllChargedTracks(esdEvent,GetAnalysisMode());
    if(!allChargedTracks) return;
    
    Int_t entries = allChargedTracks->GetEntries();
    
    labelsAll = new Int_t[entries];
    labelsAcc = new Int_t[entries];
    labelsRec = new Int_t[entries];
    for(Int_t i=0; i<entries;++i) 
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(i);
      
      if(!track) continue;
      if(track->Charge()==0) continue;
      
      if(IsUseMCInfo() && ( mcEvent->IsFromBGEvent(track->GetLabel()) ) && (track->GetLabel() >= 0) ) 
      {
	// Printf("i = %i, label = %i", i, track->GetLabel());
	continue; 
      }
      
      // only postive charged 
      if(GetParticleMode() == AlidNdPtHelper::kPlus && track->Charge() < 0) 
	continue;
      
      // only negative charged 
      if(GetParticleMode() == AlidNdPtHelper::kMinus && track->Charge() > 0) 
	continue;
      
      //
      Double_t values[4] = {vtxESD->GetZv(),track->Pt(),track->Eta(), centralityF};	  
      
      fRecTrackHist2[AlidNdPtHelper::kAllTracks]->Fill(values);
      FillHistograms(track,stack,AlidNdPtHelper::kAllTracks,centralityF); 
      labelsAll[multAll] = TMath::Abs(track->GetLabel());
      
      multAll++;      
      if(esdTrackCuts->AcceptTrack(track) && accCuts->AcceptTrack(track) && recCuts->AcceptTrackLocalTPC(track)) {
	
	fRecTrackHist2[AlidNdPtHelper::kRecTracks]->Fill(values);
	FillHistograms(track,stack,AlidNdPtHelper::kRecTracks,centralityF); 
	labelsRec[multRec] = TMath::Abs(track->GetLabel());
	
	multRec++;
	
      }
    }//loop over entries
    
    
    // fill track multiplicity histograms
    //FillHistograms(allChargedTracks,labelsAll,multAll,labelsAcc,multAcc,labelsRec,multRec);
    
    Double_t vRecEventHist1[4] = {vtxESD->GetXv(),vtxESD->GetYv(),vtxESD->GetZv(),centralityF};
    fRecEventHist1->Fill(vRecEventHist1);
    
    Double_t vRecEventHist2[3] = {vtxESD->GetZv(),Double_t(multMBTracks),centralityF};
    fRecEventHist2->Fill(vRecEventHist2);
    
  } // triggered and event vertex
  
  if(IsUseMCInfo())  
  {
    // 
    // event level corrections (zv,N_MB)
    //
    // all inelastic
    Double_t vEventMatrix[3] = {vtxMC[2],Double_t(multMBTracks),centralityF};
    fGenEventMatrix->Fill(vEventMatrix); 
    if(isEventTriggered) fTriggerEventMatrix->Fill(vEventMatrix);
    if(isEventOK && isEventTriggered) fRecEventMatrix->Fill(vEventMatrix);
    
    //
    // track-event level corrections (zv,pt,eta)
    //
    for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc) 
    {
      TParticle* particle = stack->Particle(iMc);
      if (!particle)
	continue;
      
      // only charged particles
      if(!particle->GetPDG()) continue;
      Double_t charge = particle->GetPDG()->Charge()/3.;
      if ( TMath::Abs(charge) < 0.001 )
	continue;
      
      // only postive charged 
      if(GetParticleMode() == AlidNdPtHelper::kPlus && charge < 0.) 
	continue;
      
      // only negative charged 
      if(GetParticleMode() == AlidNdPtHelper::kMinus && charge > 0.) 
	continue;
      
      // physical primary
      Bool_t prim = stack->IsPhysicalPrimary(iMc);
      if(!prim) continue;
      
      // injected signal
      Bool_t isInjected = mcEvent->IsFromBGEvent(iMc);
      if(isInjected) continue;
      
      // checked accepted
      if(accCuts->AcceptTrack(particle)) 
      {
	Double_t vTrackEventMatrix[4] = {vtxMC[2], particle->Pt(), particle->Eta(), centralityF}; 
	fGenTrackEventMatrix->Fill(vTrackEventMatrix);
	
	if(!isEventTriggered) continue;  
	fTriggerTrackEventMatrix->Fill(vTrackEventMatrix);
	
	if(!isEventOK) continue;
	fRecTrackEventMatrix->Fill(vTrackEventMatrix);
	
      }// if(accCuts->AcceptTrack(particle))
    }// for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc)
    
    // 
    // track-level corrections (zv,pt,eta)
    //
    if(isEventOK && isEventTriggered)
    {
      
      // fill MC and rec event control histograms
      if(fHistogramsOn) {
	Double_t vRecMCEventHist1[4] = {vtxESD->GetXv()-vtxMC[0],vtxESD->GetYv()-vtxMC[1],vtxESD->GetZv()-vtxMC[2], centralityF};
	fRecMCEventHist1->Fill(vRecMCEventHist1);//
	
	Double_t vRecMCEventHist2[4] = {vtxESD->GetXv()-vtxMC[0],vtxESD->GetZv()-vtxMC[2],Double_t(multMBTracks), centralityF};
	fRecMCEventHist2->Fill(vRecMCEventHist2);
	
      }//
      
      //
      // MC histograms for track efficiency studies
      //
      for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc) 
      {
	TParticle* particle = stack->Particle(iMc);
	if (!particle)
	  continue;
	
	Double_t vTrackMatrix[4] = {vtxMC[2],particle->Pt(),particle->Eta(),centralityF}; 
	
	// only charged particles
	if(!particle->GetPDG()) continue;
	Double_t charge = particle->GetPDG()->Charge()/3.;
	if (TMath::Abs(charge) < 0.001)
	  continue;
	
	// only postive charged 
	if(GetParticleMode() == AlidNdPtHelper::kPlus && charge < 0.) 
	  continue;
	
	// only negative charged 
	if(GetParticleMode() == AlidNdPtHelper::kMinus && charge > 0.) 
	  continue;
	
	// physical primary
	Bool_t prim = stack->IsPhysicalPrimary(iMc);
	
	// injected signal
	Bool_t isInjected = mcEvent->IsFromBGEvent(iMc);
	if(isInjected) continue;
	
	// check accepted
	if(accCuts->AcceptTrack(particle)) 
	{
	  
	  if( AlidNdPtHelper::IsPrimaryParticle(stack, iMc, GetParticleMode()) ) 
	    fGenPrimTrackMatrix->Fill(vTrackMatrix);
	  
	  // fill control histograms
	  if(fHistogramsOn) 
	    FillHistograms(stack,iMc,AlidNdPtHelper::kAccTracks, centralityF); 
	  
	  // check multiple found tracks
	  Int_t multCount = 0;
	  for(Int_t iRec=0; iRec<multRec; ++iRec)
	  {
	    if(iMc == labelsRec[iRec]) 
	    {
	      multCount++;
	      if(multCount>1)
	      {  
		fRecMultTrackMatrix->Fill(vTrackMatrix);
		
		// fill control histogram
		if(fHistogramsOn) {
		  Int_t pid = AlidNdPtHelper::ConvertPdgToPid(particle);
		  Double_t vMCMultRecTrackHist1[4] = {particle->Pt(), particle->Eta(), Double_t(pid), centralityF};
		  fMCMultRecTrackHist1->Fill(vMCMultRecTrackHist1);
		}
	      }
	    }
	  }
	  
	  // check reconstructed
	  for(Int_t iRec=0; iRec<multRec; ++iRec)
	  {
	    if(iMc == labelsRec[iRec]) 
	    {
	      fRecTrackMatrix->Fill(vTrackMatrix);
	      
	      if( AlidNdPtHelper::IsPrimaryParticle(stack, iMc, GetParticleMode()) ) {
		fRecPrimTrackMatrix->Fill(vTrackMatrix);
	      }
	      if(!prim) fRecSecTrackMatrix->Fill(vTrackMatrix);
	      
	      // fill control histograms
	      if(fHistogramsOn) 
		FillHistograms(stack,iMc,AlidNdPtHelper::kRecTracks, centralityF); 
	      
	      break;
	    }//if(iMc == labelsRec[iRec])
	  }//reco tracks
	}//accepted tracks
      }//stack loop
    }//is triggered
  } // end bUseMC
  
  if(allChargedTracks) delete allChargedTracks; allChargedTracks = 0;
  if(labelsAll) delete [] labelsAll; labelsAll = 0;
  if(labelsAcc) delete [] labelsAcc; labelsAcc = 0;
  if(labelsRec) delete [] labelsRec; labelsRec = 0;
  
}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb2011::FillHistograms(TObjArray *const allChargedTracks,Int_t *const labelsAll,Int_t multAll,Int_t *const labelsAcc,Int_t multAcc,Int_t *const labelsRec,Int_t multRec, Float_t centralityF) {
  // multiplicity  histograms
  
  
  if(!allChargedTracks) return;
  if(!labelsAll) return;
  if(!labelsAcc) return;
  if(!labelsRec) return;
  
  Int_t entries = allChargedTracks->GetEntries();
  for(Int_t i=0; i<entries; ++i)
  {
    AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(i);
    if(!track) continue;
    if(track->Charge() == 0) continue;
    
    Int_t label = TMath::Abs(track->GetLabel());
    for(Int_t iAll=0; iAll<multAll; ++iAll) {
      if(label == labelsAll[iAll]) {
	Double_t v1[3] = {track->Pt(), Double_t(multAll), centralityF}; 
	fRecTrackMultHist1[AlidNdPtHelper::kAllTracks]->Fill(v1);
      }
    }
    for(Int_t iAcc=0; iAcc<multAcc; ++iAcc) {
      if(label == labelsAcc[iAcc]) {
	Double_t v2[3] = {track->Pt(), Double_t(multAcc), centralityF}; 
	fRecTrackMultHist1[AlidNdPtHelper::kAccTracks]->Fill(v2);
      }
    }
    for(Int_t iRec=0; iRec<multRec; ++iRec) {
      if(label == labelsRec[iRec]) {
	Double_t v3[3] = {track->Pt(), Double_t(multRec), centralityF}; 
	fRecTrackMultHist1[AlidNdPtHelper::kRecTracks]->Fill(v3);
      }//out
    }
  }
}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb2011::FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, AlidNdPtHelper::TrackObject trackObj, Float_t centralityF)
{
  
  //
  // Fill ESD track and MC histograms 
  //
  if(!esdTrack) return;
  
  Float_t q = esdTrack->Charge();
  if(TMath::Abs(q) < 0.001) return;
  
  Float_t pt = esdTrack->Pt();
  Float_t eta = esdTrack->Eta();
  Float_t phi = esdTrack->Phi();
  
  Float_t dca[2], bCov[3];
  esdTrack->GetImpactParameters(dca,bCov);
  
  Int_t nClust = esdTrack->GetTPCclusters(0);
  Float_t chi2PerCluster = 0.;
  if(nClust>0.) chi2PerCluster = esdTrack->GetTPCchi2()/Float_t(nClust);
  
  
  // fill histograms
  Double_t values[4] = {pt,eta,phi,centralityF};	  
  fRecTrackHist1[trackObj]->Fill(values);
  
  //
  // Fill rec vs MC information
  //
  if(!stack) return;
  
  Int_t label = TMath::Abs(esdTrack->GetLabel()); 
  //if(label == 0) return;
  
  if(label > stack->GetNtrack()) return;
  TParticle* particle = stack->Particle(label);
  if(!particle) return;
  
  Int_t motherPdg = -1;
  TParticle* mother = 0;
  
  //TParticle* prim_mother = AlidNdPtHelper::FindPrimaryMother(stack,label);
  Int_t motherLabel = particle->GetMother(0); 
  if(motherLabel>0) mother = stack->Particle(motherLabel);
  if(mother) motherPdg = TMath::Abs(mother->GetPdgCode()); // take abs for visualisation only
    //Int_t mech = particle->GetUniqueID(); // production mechanism
    
    if(!particle->GetPDG()) return;
    Double_t gq = particle->GetPDG()->Charge()/3.0; // Charge units |e|/3 
    if(TMath::Abs(gq)<0.001) return;
    Float_t gpt = particle->Pt();
  Float_t geta = particle->Eta();
  
  Double_t dpt=0;
  //printf("pt %f, gpt %f \n",pt,gpt);
  if(gpt) dpt = (pt-gpt)/gpt;
  Double_t deta = (eta-geta);
  
  // fill histograms
  if(trackObj == AlidNdPtHelper::kRecTracks)  //RecTracks???
  {
    Double_t vTrackPtCorrelationMatrix[4]={pt,gpt,geta,centralityF};
    fTrackPtCorrelationMatrix->Fill(vTrackPtCorrelationMatrix);
    
    Double_t vRecMCTrackHist1[5]={gpt,geta,dpt,deta,centralityF};
    fRecMCTrackHist1->Fill(vRecMCTrackHist1);
  }
}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb2011::FillHistograms(AliStack *const stack, Int_t label, AlidNdPtHelper::TrackObject trackObj, Float_t centralityF)
{
  
  // Fill MC histograms
  if(!stack) return;
  
  if(label > stack->GetNtrack()) return;
  TParticle* particle = stack->Particle(label);
  if(!particle) return;
  
  Int_t motherPdg = -1;
  TParticle* mother = 0;
  
  //TParticle* prim_mother = AlidNdPtHelper::FindPrimaryMother(stack,label);
  Int_t motherLabel = particle->GetMother(0); 
  if(motherLabel>0) mother = stack->Particle(motherLabel);
  if(mother) motherPdg = TMath::Abs(mother->GetPdgCode()); // take abs for visualisation only
    Int_t mech = particle->GetUniqueID(); // production mechanism
    
    if(!particle->GetPDG()) return;
    Double_t gq = particle->GetPDG()->Charge()/3.0; // Charge units |e|/3 
    if(TMath::Abs(gq) < 0.001) return;
    
    Float_t gpt = particle->Pt();
  //Float_t qgpt = particle->Pt() * gq;
  Float_t geta = particle->Eta();
  Float_t gphi = particle->Phi();
  //Float_t gpz = particle->Pz();
  
  Bool_t prim = stack->IsPhysicalPrimary(label);
  //Float_t vx = particle->Vx(); Float_t vy = particle->Vy(); Float_t vz = particle->Vz();
  
  Int_t pid = AlidNdPtHelper::ConvertPdgToPid(particle);
  
  //if(prim&&pid==5) printf("pdgcode %d, production mech %d \n",particle->GetPdgCode(),mech);
  //if(!prim) printf("motherPdg %d, particle %d, production mech %d\n",motherPdg, particle->GetPdgCode(),mech);
  
  //
  // fill histogram
  //
  Double_t vMCTrackHist1[4] = {gpt,geta,gphi,centralityF};
  fMCTrackHist1[trackObj]->Fill(vMCTrackHist1);
  
  Double_t vMCPrimTrackHist1[6] = {gpt,geta,Double_t(pid),Double_t(mech),Double_t(motherPdg),centralityF};
  Double_t vMCPrimTrackHist2[4] = {Double_t(TMath::Abs(particle->GetPdgCode())),Double_t(mech),Double_t(motherPdg),centralityF};
  
  //if(prim && AliPWG0Helper::IsPrimaryCharged(particle, label)) fMCPrimTrackHist1[trackObj]->Fill(vMCPrimTrackHist1);
  
  if(prim) { 
    fMCPrimTrackHist1[trackObj]->Fill(vMCPrimTrackHist1);
    if(pid == 5) fMCPrimTrackHist2[trackObj]->Fill(vMCPrimTrackHist2);
  }
  else { 
    fMCSecTrackHist1[trackObj]->Fill(vMCPrimTrackHist1);
  }
  
}

//_____________________________________________________________________________
Long64_t AlidNdPtAnalysisPbPb2011::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)
  
  //  init if not done already
  if (!fIsInit) { Init(); }
  
  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;
  
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  
  //
  //TList *collPhysSelection = new TList;
  
  // collection of generated histograms
  
  Int_t count=0;
  while((obj = iter->Next()) != 0) {
    AlidNdPtAnalysisPbPb2011* entry = dynamic_cast<AlidNdPtAnalysisPbPb2011*>(obj);
    if (entry == 0) continue; 
    
    // physics selection
    //printf("entry->GetPhysicsTriggerSelection() %p \n", entry->GetPhysicsTriggerSelection());
    //AliPhysicsSelection *physSel = entry->GetPhysicsTriggerSelection();
    //if( physSel ) collPhysSelection->Add(physSel); 
    
    //
    fTrackPtCorrelationMatrix->Add(entry->fTrackPtCorrelationMatrix);
    
    //
    fGenEventMatrix->Add(entry->fGenEventMatrix);
    
    fTriggerEventMatrix->Add(entry->fTriggerEventMatrix);
    
    fRecEventMatrix->Add(entry->fRecEventMatrix);
    //
    fGenTrackEventMatrix->Add(entry->fGenTrackEventMatrix);
    
    fTriggerTrackEventMatrix->Add(entry->fTriggerTrackEventMatrix);
    
    fRecTrackEventMatrix->Add(entry->fRecTrackEventMatrix);
    
    //
    fGenTrackMatrix->Add(entry->fGenTrackMatrix);
    fGenPrimTrackMatrix->Add(entry->fGenPrimTrackMatrix);
    fRecPrimTrackMatrix->Add(entry->fRecPrimTrackMatrix);
    //
    fRecTrackMatrix->Add(entry->fRecTrackMatrix);
    fRecSecTrackMatrix->Add(entry->fRecSecTrackMatrix);
    //
    fRecMultTrackMatrix->Add(entry->fRecMultTrackMatrix);
    
    //
    // control analysis histograms
    //
    fMCEventHist1->Add(entry->fMCEventHist1);
    fRecEventHist1->Add(entry->fRecEventHist1);
    fRecEventHist2->Add(entry->fRecEventHist2);
    fRecMCEventHist1->Add(entry->fRecMCEventHist1);
    fRecMCEventHist2->Add(entry->fRecMCEventHist2);
    
    
    for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) {
      fMCTrackHist1[i]->Add(entry->fMCTrackHist1[i]);
      
      fMCPrimTrackHist1[i]->Add(entry->fMCPrimTrackHist1[i]);
      fMCPrimTrackHist2[i]->Add(entry->fMCPrimTrackHist2[i]);
      fMCSecTrackHist1[i]->Add(entry->fMCSecTrackHist1[i]);
      
      fRecTrackHist1[i]->Add(entry->fRecTrackHist1[i]);
      fRecTrackHist2[i]->Add(entry->fRecTrackHist2[i]);
      fRecTrackMultHist1[i]->Add(entry->fRecTrackMultHist1[i]);
    }
    fRecMCTrackHist1->Add(entry->fRecMCTrackHist1);
    fMCMultRecTrackHist1->Add(entry->fMCMultRecTrackHist1);
    fRecTrackHist3->Add(entry->fRecTrackHist3);
    
    count++;
  }
  
  //AliPhysicsSelection *trigSelection = GetPhysicsTriggerSelection();
  //if( trigSelection ) trigSelection->Merge(collPhysSelection);
  //if(collPhysSelection) delete collPhysSelection;
  
  return count;
}



//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb2011::Analyse() 
{
  
  //  init if not done already
  if (!fIsInit) { Init(); }
  
  // Analyse histograms
  //
  TH1::AddDirectory(kFALSE);
  TH1 *h=0, *h1=0, *h2=0, *h2c = 0; 
  THnSparse *hs=0; 
  TH2 *h2D=0; 
  
  char name[256];
  TObjArray *aFolderObj = new TObjArray;
  if(!aFolderObj) return;
  
  //
  // LHC backgraund in all and 0-bins
  //
  //AliPhysicsSelection *trigSel = GetPhysicsTriggerSelection();
  //trigSel->SaveHistograms("physics_selection");
  
  //
  // Reconstructed event vertex
  //
  h = fRecEventHist1->Projection(0);
  h->SetName("Xv");
  aFolderObj->Add(h);
  
  h = fRecEventHist1->Projection(1);
  h->SetName("Yv");
  aFolderObj->Add(h);
  
  h = fRecEventHist1->Projection(2);
  h->SetName("Zv");
  aFolderObj->Add(h);
  
  //
  // multiplicity
  //
  h = fRecEventHist2->Projection(1);
  h->SetName("multMB");
  aFolderObj->Add(h);
  
  h2D = fRecEventHist2->Projection(0,1); 
  h2D->SetName("Zv_vs_multiplicity_MB");
  aFolderObj->Add(h2D);
  
  //
  // reconstructed pt histograms
  //
  h = fRecTrackHist1[0]->Projection(0);
  h->Scale(1.,"width");
  h->SetName("pt_all_ch");
  aFolderObj->Add(h);
  
  h = fRecTrackHist1[1]->Projection(0);
  h->Scale(1.,"width");
  h->SetName("pt_acc");
  aFolderObj->Add(h);
  
  h = fRecTrackHist1[2]->Projection(0);
  h->Scale(1.,"width");
  h->SetName("pt_rec");
  aFolderObj->Add(h);
  
  //
  // reconstructed eta histograms
  //
  h = fRecTrackHist1[0]->Projection(1);
  h->SetName("eta_all_ch");
  aFolderObj->Add(h);
  
  h = fRecTrackHist1[1]->Projection(1);
  h->SetName("eta_acc");
  aFolderObj->Add(h);
  
  h = fRecTrackHist1[2]->Projection(1);
  h->SetName("eta_rec");
  aFolderObj->Add(h);
  
  //
  // reconstructed phi histograms
  //
  h = fRecTrackHist1[0]->Projection(2);
  h->SetName("phi_all_ch");
  aFolderObj->Add(h);
  
  h = fRecTrackHist1[1]->Projection(2);
  h->SetName("phi_acc");
  aFolderObj->Add(h);
  
  h = fRecTrackHist1[2]->Projection(2);
  h->SetName("phi_rec");
  aFolderObj->Add(h);
  
  //
  // reconstructed eta:pt histograms
  //
  h2D = fRecTrackHist1[0]->Projection(1,0);
  h2D->SetName("pt_eta_all_ch");
  aFolderObj->Add(h2D);
  
  h2D = fRecTrackHist1[1]->Projection(1,0);
  h2D->SetName("pt_eta_acc");
  aFolderObj->Add(h2D);
  
  h2D = fRecTrackHist1[2]->Projection(1,0);
  h2D->SetName("pt_eta_rec");
  aFolderObj->Add(h2D);
  
  //
  // reconstructed phi:pt histograms
  //
  h2D = fRecTrackHist1[0]->Projection(2,0);
  h2D->SetName("pt_phi_all_ch");
  aFolderObj->Add(h2D);
  
  h2D = fRecTrackHist1[1]->Projection(2,0);
  h2D->SetName("pt_phi_acc");
  aFolderObj->Add(h2D);
  
  h2D = fRecTrackHist1[2]->Projection(2,0);
  h2D->SetName("pt_phi_rec");
  aFolderObj->Add(h2D);
  
  //
  // reconstructed phi:eta histograms
  //
  h2D = fRecTrackHist1[0]->Projection(2,1);
  h2D->SetName("eta_phi_all_ch");
  aFolderObj->Add(h2D);
  
  h2D = fRecTrackHist1[1]->Projection(2,1);
  h2D->SetName("eta_phi_acc");
  aFolderObj->Add(h2D);
  
  h2D = fRecTrackHist1[2]->Projection(2,1);
  h2D->SetName("eta_phi_rec");
  aFolderObj->Add(h2D);
  
  //
  // reconstructed nClust, chi2 vs pt, eta, phi
  //
  if(fHistogramsOn) {
    
    h2D = fRecTrackHist3->Projection(0,1);
    h2D->SetName("nClust_chi2_rec");
    aFolderObj->Add(h2D);
    
    h2D = fRecTrackHist3->Projection(0,2);
    h2D->SetName("nClust_pt_rec");
    aFolderObj->Add(h2D);
    
    h2D = fRecTrackHist3->Projection(0,3);
    h2D->SetName("nClust_eta_rec");
    aFolderObj->Add(h2D);
    
    h2D = fRecTrackHist3->Projection(0,4);
    h2D->SetName("nClust_phi_rec");
    aFolderObj->Add(h2D);
    
    h2D = fRecTrackHist3->Projection(1,2);
    h2D->SetName("chi2_pt_rec");
    aFolderObj->Add(h2D);
    
    h2D = fRecTrackHist3->Projection(1,3);
    h2D->SetName("chi2_eta_rec");
    aFolderObj->Add(h2D);
    
    h2D = fRecTrackHist3->Projection(1,4);
    h2D->SetName("chi2_phi_rec");
    aFolderObj->Add(h2D);
    
  }
  
  //
  // calculate corrections for empty events
  // with multMB==0 
  //
  
  //
  // normalised zv to generate zv for triggered events
  //
  h = fRecEventHist2->Projection(0);
  if( h->Integral() ) h->Scale(1./h->Integral());
  h->SetName("zv_distribution_norm");
  aFolderObj->Add(h);
  
  //
  // MC available
  //
  if(IsUseMCInfo()) {
    
    //
    // Event vertex resolution
    //
    h2D = fRecMCEventHist2->Projection(0,2);
    h2D->SetName("DeltaXv_vs_mult");
    aFolderObj->Add(h2D);
    
    h2D = fRecMCEventHist2->Projection(1,2);
    h2D->SetName("DeltaZv_vs_mult");
    aFolderObj->Add(h2D);
    
    //
    // normalised zv to get trigger/trigger+vertex event differences
    // F(zv) = E_trig(zv,0)/Int(E_trig(zv,0) / Sum(E_trigvtx(zv,n))/Sum(Int(E_trigvtx(zv,n))dzv)
    //
    fTriggerEventMatrix->GetAxis(1)->SetRangeUser(0.,0.);
    h = fTriggerEventMatrix->Projection(0);
    h2D = fTriggerEventMatrix->Projection(0,1);
    if(h2D->Integral()) h->Scale(1./h2D->Integral());
    
    h1 = fRecEventMatrix->Projection(0);
    h2D = fRecEventMatrix->Projection(0,1);
    if(h2D->Integral()) h1->Scale(1./h2D->Integral());
    
    h->Divide(h1);
    h->SetName("zv_empty_events_norm");
    aFolderObj->Add(h);
    
    fTriggerEventMatrix->GetAxis(1)->SetRange(1,fTriggerEventMatrix->GetAxis(1)->GetNbins());
    
    //
    // rec. vs true track pt correlation matrix
    //
    hs = (THnSparse*)fTrackPtCorrelationMatrix->Clone("track_pt_correlation_matrix");
    aFolderObj->Add(hs);
    
    //
    // trigger efficiency for INEL
    //
    h = AlidNdPtHelper::GenerateCorrMatrix(fTriggerEventMatrix->Projection(0),fGenEventMatrix->Projection(0),"zv_trig_INEL_eff_matrix");
    aFolderObj->Add(h);
    
    
    //
    // trigger bias correction (MB to INEL)
    //
    hs = AlidNdPtHelper::GenerateCorrMatrix(fGenEventMatrix,fTriggerEventMatrix,"zv_mult_trig_MBtoInel_corr_matrix");
    aFolderObj->Add(hs);
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fGenEventMatrix->Projection(0),fTriggerEventMatrix->Projection(0),"zv_trig_MBtoInel_corr_matrix");
    aFolderObj->Add(h);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenEventMatrix->Projection(0,1),fTriggerEventMatrix->Projection(0,1),"zv_mult_trig_MBtoInel_corr_matrix_2D");
    aFolderObj->Add(h2D);
    
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fGenEventMatrix->Projection(1),fTriggerEventMatrix->Projection(1),"mult_trig_MBtoInel_corr_matrix");
    aFolderObj->Add(h);
    
    
    //
    // event vertex reconstruction correction (MB)
    //
    hs = AlidNdPtHelper::GenerateCorrMatrix(fTriggerEventMatrix,fRecEventMatrix,"zv_mult_event_corr_matrix");
    aFolderObj->Add(hs);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fTriggerEventMatrix->Projection(0,1),fRecEventMatrix->Projection(0,1),"zv_mult_event_corr_matrix_2D");
    aFolderObj->Add(h2D);
    
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fTriggerEventMatrix->Projection(1),fRecEventMatrix->Projection(1),"mult_event_corr_matrix");
    aFolderObj->Add(h);
    
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fTriggerEventMatrix->Projection(0),fRecEventMatrix->Projection(0),"zv_event_corr_matrix");
    aFolderObj->Add(h);
    
    
    
    //
    // track-event trigger bias correction (MB to INEL)
    //
    hs = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackEventMatrix,fTriggerTrackEventMatrix,"zv_pt_eta_track_trig_MBtoInel_corr_matrix");
    aFolderObj->Add(hs);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackEventMatrix->Projection(1,2),fTriggerTrackEventMatrix->Projection(1,2),"eta_pt_track_trig_MBtoInel_corr_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackEventMatrix->Projection(1,0),fTriggerTrackEventMatrix->Projection(1,0),"pt_zv_track_trig_MBtoInel_corr_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackEventMatrix->Projection(2,0),fTriggerTrackEventMatrix->Projection(2,0),"zv_eta_track_trig_MBtoInel_corr_matrix");
    aFolderObj->Add(h2D);
    
    // efficiency
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fTriggerTrackEventMatrix->Projection(1),fGenTrackEventMatrix->Projection(1),"pt_track_trig_MBtoInel_eff_matrix");
    aFolderObj->Add(h);
    
    
    //
    // track-event vertex reconstruction correction (MB)
    //
    hs = AlidNdPtHelper::GenerateCorrMatrix(fTriggerTrackEventMatrix,fRecTrackEventMatrix,"zv_pt_eta_track_event_corr_matrix");
    aFolderObj->Add(hs);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fTriggerTrackEventMatrix->Projection(1,2),fRecTrackEventMatrix->Projection(1,2),"eta_pt_track_event_corr_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fTriggerTrackEventMatrix->Projection(1,0),fRecTrackEventMatrix->Projection(1,0),"pt_zv_track_event_corr_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fTriggerTrackEventMatrix->Projection(2,0),fRecTrackEventMatrix->Projection(2,0),"zv_eta_track_event_corr_matrix");
    aFolderObj->Add(h2D);
    
    // efficiency
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fRecTrackEventMatrix->Projection(1),fTriggerTrackEventMatrix->Projection(1),"pt_track_event_eff_matrix");
    aFolderObj->Add(h);
    
    
    //
    // track rec. efficiency correction
    //
    hs = AlidNdPtHelper::GenerateCorrMatrix(fGenPrimTrackMatrix,fRecPrimTrackMatrix,"zv_pt_eta_track_corr_matrix");
    aFolderObj->Add(hs);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenPrimTrackMatrix->Projection(1,2),fRecPrimTrackMatrix->Projection(1,2),"eta_pt_track_corr_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenPrimTrackMatrix->Projection(1,0),fRecPrimTrackMatrix->Projection(1,0),"pt_zv_track_corr_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenPrimTrackMatrix->Projection(2,0),fRecPrimTrackMatrix->Projection(2,0),"zv_eta_track_corr_matrix");
    aFolderObj->Add(h2D);
    
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fGenPrimTrackMatrix->Projection(0),fRecPrimTrackMatrix->Projection(0),"zv_track_corr_matrix");
    aFolderObj->Add(h);
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fGenPrimTrackMatrix->Projection(1),fRecPrimTrackMatrix->Projection(1),"pt_track_corr_matrix");
    aFolderObj->Add(h);
    
    // efficiency
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fRecPrimTrackMatrix->Projection(1), fGenPrimTrackMatrix->Projection(1),"pt_track_eff_matrix");
    aFolderObj->Add(h);
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fGenPrimTrackMatrix->Projection(2),fRecPrimTrackMatrix->Projection(2),"eta_track_corr_matrix");
    aFolderObj->Add(h);
    
    //
    // secondary track contamination correction
    //
    //hs = AlidNdPtHelper::GenerateContCorrMatrix(fRecSecTrackMatrix,fRecTrackMatrix,"zv_pt_eta_track_cont_matrix");
    hs = AlidNdPtHelper::GenerateCorrMatrix(fRecSecTrackMatrix,fRecTrackMatrix,"zv_pt_eta_track_cont_matrix");
    aFolderObj->Add(hs);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fRecSecTrackMatrix->Projection(1,2),fRecTrackMatrix->Projection(1,2),"eta_pt_track_cont_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fRecSecTrackMatrix->Projection(1,0),fRecTrackMatrix->Projection(1,0),"pt_zv_track_cont_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fRecSecTrackMatrix->Projection(2,0),fRecTrackMatrix->Projection(2,0),"zv_eta_track_cont_matrix");
    aFolderObj->Add(h2D);
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fRecSecTrackMatrix->Projection(0),fRecTrackMatrix->Projection(0),"zv_track_cont_matrix");
    aFolderObj->Add(h);
    
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fRecSecTrackMatrix->Projection(1),fRecTrackMatrix->Projection(1),"pt_track_cont_matrix");
    aFolderObj->Add(h);
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fRecSecTrackMatrix->Projection(2),fRecTrackMatrix->Projection(2),"eta_track_cont_matrix");
    aFolderObj->Add(h);
    
    //
    // multiple track reconstruction correction
    //
    //hs = AlidNdPtHelper::GenerateContCorrMatrix(fRecMultTrackMatrix,fRecTrackMatrix,"zv_pt_eta_mult_track_cont_matrix");
    hs = AlidNdPtHelper::GenerateCorrMatrix(fRecMultTrackMatrix,fRecTrackMatrix,"zv_pt_eta_mult_track_cont_matrix");
    aFolderObj->Add(hs);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fRecMultTrackMatrix->Projection(1,2),fRecTrackMatrix->Projection(1,2),"eta_pt_mult_track_cont_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fRecMultTrackMatrix->Projection(1,0),fRecTrackMatrix->Projection(1,0),"pt_zv_mult_track_cont_matrix");
    aFolderObj->Add(h2D);
    
    h2D = AlidNdPtHelper::GenerateCorrMatrix(fRecMultTrackMatrix->Projection(2,0),fRecTrackMatrix->Projection(2,0),"zv_eta_mult_track_cont_matrix");
    aFolderObj->Add(h2D);
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fRecMultTrackMatrix->Projection(0),fRecTrackMatrix->Projection(0),"zv_mult_track_cont_matrix");
    aFolderObj->Add(h);
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fRecMultTrackMatrix->Projection(1),fRecTrackMatrix->Projection(1),"pt_mult_track_cont_matrix");
    aFolderObj->Add(h);
    
    h = AlidNdPtHelper::GenerateCorrMatrix(fRecMultTrackMatrix->Projection(2),fRecTrackMatrix->Projection(2),"eta_mult_track_cont_matrix");
    aFolderObj->Add(h);
    
    //
    // Control histograms
    //
    
    if(fHistogramsOn) {
      
      // Efficiency electrons, muons, pions, kaons, protons, all
      fMCPrimTrackHist1[1]->GetAxis(2)->SetRange(1,1); 
      fMCPrimTrackHist1[2]->GetAxis(2)->SetRange(1,1); 
      h1 = fMCPrimTrackHist1[1]->Projection(0);
      h2 = fMCPrimTrackHist1[2]->Projection(0);
      h2c = (TH1D *)h2->Clone();
      h2c->Divide(h1);
      h2c->SetName("eff_pt_electrons");
      aFolderObj->Add(h2c);
      
      fMCPrimTrackHist1[1]->GetAxis(2)->SetRange(2,2); 
      fMCPrimTrackHist1[2]->GetAxis(2)->SetRange(2,2); 
      h1 = fMCPrimTrackHist1[1]->Projection(0);
      h2 = fMCPrimTrackHist1[2]->Projection(0);
      h2c = (TH1D *)h2->Clone();
      h2c->Divide(h1);
      h2c->SetName("eff_pt_muons");
      aFolderObj->Add(h2c);
      
      fMCPrimTrackHist1[1]->GetAxis(2)->SetRange(3,3); 
      fMCPrimTrackHist1[2]->GetAxis(2)->SetRange(3,3); 
      h1 = fMCPrimTrackHist1[1]->Projection(0);
      h2 = fMCPrimTrackHist1[2]->Projection(0);
      h2c = (TH1D *)h2->Clone();
      h2c->Divide(h1);
      h2c->SetName("eff_pt_pions");
      aFolderObj->Add(h2c);
      
      fMCPrimTrackHist1[1]->GetAxis(2)->SetRange(4,4); 
      fMCPrimTrackHist1[2]->GetAxis(2)->SetRange(4,4); 
      h1 = fMCPrimTrackHist1[1]->Projection(0);
      h2 = fMCPrimTrackHist1[2]->Projection(0);
      h2c = (TH1D *)h2->Clone();
      h2c->Divide(h1);
      h2c->SetName("eff_pt_kaons");
      aFolderObj->Add(h2c);
      
      fMCPrimTrackHist1[1]->GetAxis(2)->SetRange(5,5); 
      fMCPrimTrackHist1[2]->GetAxis(2)->SetRange(5,5); 
      h1 = fMCPrimTrackHist1[1]->Projection(0);
      h2 = fMCPrimTrackHist1[2]->Projection(0);
      h2c = (TH1D *)h2->Clone();
      h2c->Divide(h1);
      h2c->SetName("eff_pt_protons");
      aFolderObj->Add(h2c);
      
      fMCPrimTrackHist1[1]->GetAxis(2)->SetRange(1,5); 
      fMCPrimTrackHist1[2]->GetAxis(2)->SetRange(1,5); 
      h1 = fMCPrimTrackHist1[1]->Projection(0);
      h2 = fMCPrimTrackHist1[2]->Projection(0);
      h2c = (TH1D *)h2->Clone();
      h2c->Divide(h1);
      h2c->SetName("eff_pt_selected");
      aFolderObj->Add(h2c);
      
      fMCPrimTrackHist1[1]->GetAxis(2)->SetRange(1,6); 
      fMCPrimTrackHist1[2]->GetAxis(2)->SetRange(1,6); 
      h1 = fMCPrimTrackHist1[1]->Projection(0);
      h2 = fMCPrimTrackHist1[2]->Projection(0);
      h2c = (TH1D *)h2->Clone();
      h2c->Divide(h1);
      h2c->SetName("eff_pt_all");
      aFolderObj->Add(h2c);
      
      fMCPrimTrackHist1[1]->GetAxis(1)->SetRange(1,fMCPrimTrackHist1[1]->GetAxis(1)->GetNbins()); 
      fMCPrimTrackHist1[2]->GetAxis(1)->SetRange(1,fMCPrimTrackHist1[2]->GetAxis(1)->GetNbins());
      
      //  pt spetra
      // - rec, primaries, secondaries
      // - primaries (pid) 
      // - secondaries (pid)
      // - secondaries (mech)
      // - secondaries (mother)
      //
      
      TH1D *mcPtAccall = fMCTrackHist1[1]->Projection(0);
      mcPtAccall->SetName("mc_pt_acc_all");
      aFolderObj->Add(mcPtAccall);
      
      TH1D *mcPtAccprim = fMCPrimTrackHist1[1]->Projection(0);
      mcPtAccprim->SetName("mc_pt_acc_prim");
      aFolderObj->Add(mcPtAccprim);
      
      TH1D *mcPtRecall = fMCTrackHist1[2]->Projection(0);
      mcPtRecall->SetName("mc_pt_rec_all");
      aFolderObj->Add(mcPtRecall);
      
      TH1D *mcPtRecprim = fMCPrimTrackHist1[2]->Projection(0);
      mcPtRecprim->SetName("mc_pt_rec_prim");
      aFolderObj->Add(mcPtRecprim);
      
      TH1D *mcPtRecsec = fMCSecTrackHist1[2]->Projection(0);
      mcPtRecsec->SetName("mc_pt_rec_sec");
      aFolderObj->Add(mcPtRecsec);
      
      for(Int_t i = 0; i<6; i++) 
      { 
	snprintf(name,256,"mc_pt_rec_prim_pid_%d",i); 
	fMCPrimTrackHist1[2]->GetAxis(2)->SetRange(i+1,i+1);
	h = fMCPrimTrackHist1[2]->Projection(0);
	h->SetName(name);
	aFolderObj->Add(h);
	
	snprintf(name,256,"mc_pt_rec_sec_pid_%d",i); 
	fMCSecTrackHist1[2]->GetAxis(2)->SetRange(i+1,i+1);
	h = fMCSecTrackHist1[2]->Projection(0);
	h->SetName(name);
	aFolderObj->Add(h);
	
	// production mechanisms for given pid
	fMCSecTrackHist1[2]->GetAxis(2)->SetRange(i+1,i+1);
	
	for(Int_t j=0; j<20; j++) {
	  if(j == 4) {
	    // decay
	    
	    snprintf(name,256,"mc_pt_rec_sec_pid_%d_decay",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(0);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	    snprintf(name,256,"mc_eta_rec_sec_pid_%d_decay",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(1);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	    snprintf(name,256,"mc_mother_rec_sec_pid_%d_decay",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(4);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	  } else if (j == 5) {
	    // conversion
	    
	    snprintf(name,256,"mc_pt_rec_sec_pid_%d_conv",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(0);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	    snprintf(name,256,"mc_eta_rec_sec_pid_%d_conv",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(1);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	    snprintf(name,256,"mc_mother_rec_sec_pid_%d_conv",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(4);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	  } else if (j == 13) {
	    // mat
	    
	    snprintf(name,256,"mc_pt_rec_sec_pid_%d_mat",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(0);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	    snprintf(name,256,"mc_eta_rec_sec_pid_%d_mat",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(1);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	    snprintf(name,256,"mc_eta_mother_rec_sec_pid_%d_mat",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(4,1);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	    snprintf(name,256,"mc_mother_rec_sec_pid_%d_mat",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(4);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	    snprintf(name,256,"mc_pt_mother_rec_sec_pid_%d_mat",i); 
	    fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
	    h = fMCSecTrackHist1[2]->Projection(4,0);
	    h->SetName(name);
	    aFolderObj->Add(h);
	    
	  } else {
	    continue;
	  }
	}
	
      }
    } // end fHistogramOn
    
    //
    //  resolution histograms
    //  only for reconstructed tracks
    //
    
    TH2F *h2F=0;
    TCanvas * c = new TCanvas("resol","resol");
    c->cd();
    
    //
    fRecMCTrackHist1->GetAxis(1)->SetRangeUser(-0.8,0.79); 
    
    h2F = (TH2F*)fRecMCTrackHist1->Projection(2,0);
    h = AlidNdPtHelper::MakeResol(h2F,1,0,kTRUE,10);
    h->SetXTitle("p_{tmc} (GeV/c)");
    h->SetYTitle("(p_{t}-p_{tmc})/p_{tmc} resolution");
    h->Draw();
    h->SetName("pt_resolution_vs_mcpt");
    aFolderObj->Add(h);
    
    h2F = (TH2F*)fRecMCTrackHist1->Projection(2,0);
    h = AlidNdPtHelper::MakeResol(h2F,1,1,kTRUE,10);
    h->SetXTitle("p_{tmc} (GeV/c)");
    h->SetYTitle("(p_{t}-p_{tmc})/p_{tmc} mean");
    h->Draw();
    h->SetName("dpt_mean_vs_mcpt");
    aFolderObj->Add(h);
    
    //
    h2F = (TH2F*)fRecMCTrackHist1->Projection(3,0);
    h = AlidNdPtHelper::MakeResol(h2F,1,0,kTRUE,10);
    h->SetXTitle("p_{tmc} (GeV/c)");
    h->SetYTitle("(#eta-#eta_{mc}) resolution");
    h->Draw();
    h->SetName("eta_resolution_vs_mcpt");
    aFolderObj->Add(h);
    
    h2F = (TH2F*)fRecMCTrackHist1->Projection(3,0);
    h = AlidNdPtHelper::MakeResol(h2F,1,1,kTRUE,10);
    h->SetXTitle("p_{tmc} (GeV/c)");
    h->SetYTitle("(#eta-mc#eta) mean");
    h->Draw();
    h->SetName("deta_mean_vs_mcpt");
    aFolderObj->Add(h);
    
    // 
    fRecMCTrackHist1->GetAxis(1)->SetRange(1,fRecMCTrackHist1->GetAxis(1)->GetNbins()); 
    
    h2F = (TH2F*)fRecMCTrackHist1->Projection(2,1);
    h = AlidNdPtHelper::MakeResol(h2F,1,0,kTRUE,10);
    h->SetXTitle("#eta_{mc}");
    h->SetYTitle("(p_{t}-p_{tmc})/p_{tmc} resolution");
    h->Draw();
    h->SetName("pt_resolution_vs_mceta");
    aFolderObj->Add(h);
    
    h2F = (TH2F*)fRecMCTrackHist1->Projection(2,1);
    h = AlidNdPtHelper::MakeResol(h2F,1,1,kTRUE,10);
    h->SetXTitle("#eta_{mc}");
    h->SetYTitle("(p_{t}-p_{tmc})/p_{tmc} mean");
    h->Draw();
    h->SetName("dpt_mean_vs_mceta");
    aFolderObj->Add(h);
    
    //
    h2F = (TH2F*)fRecMCTrackHist1->Projection(3,1);
    h = AlidNdPtHelper::MakeResol(h2F,1,0,kTRUE,10);
    h->SetXTitle("#eta_{mc}");
    h->SetYTitle("(#eta-#eta_{mc}) resolution");
    h->Draw();
    h->SetName("eta_resolution_vs_mceta");
    aFolderObj->Add(h);
    
    h2F = (TH2F*)fRecMCTrackHist1->Projection(3,1);
    h = AlidNdPtHelper::MakeResol(h2F,1,1,kTRUE,10);
    h->SetXTitle("#eta_{mc}");
    h->SetYTitle("(#eta-mc#eta) mean");
    h->Draw();
    h->SetName("deta_mean_vs_mceta");
    aFolderObj->Add(h);
    
    fRecMCTrackHist1->GetAxis(0)->SetRange(1,fRecMCTrackHist1->GetAxis(0)->GetNbins()); 
    
  } // end use MC info
  
  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);
  
  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AlidNdPtAnalysisPbPb2011::ExportToFolder(TObjArray * const array) 
{
  // recreate folder avery time and export objects to new one
  //
  AlidNdPtAnalysisPbPb2011 * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();
  
  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();
  
  if(folder) { 
    // get name and title from old folder
    name = folder->GetName();  
    title = folder->GetTitle();  
    
    // delete old one
    delete folder;
    
    // create new one
    newFolder = CreateFolder(name.Data(),title.Data());
    newFolder->SetOwner();
    
    // add objects to folder
    while(i < size) {
      newFolder->Add(array->At(i));
      i++;
    }
  }
  
  return newFolder;
}

//_____________________________________________________________________________
TFolder* AlidNdPtAnalysisPbPb2011::CreateFolder(TString name,TString title) { 
  // create folder for analysed histograms
  //
  TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());
  
  return folder;
}
