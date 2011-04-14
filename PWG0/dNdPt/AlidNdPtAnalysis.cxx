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
// AlidNdPtAnalysis class. 
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
#include "AliInputEventHandler.h"  
#include "AliAnalysisManager.h"  
#include "AliGenEventHeader.h"  
#include "AliStack.h"  
#include "AliESDEvent.h"  
#include "AliMCEvent.h"  
#include "AliESDtrackCuts.h"  
#include "AliLog.h" 
#include "AliMultiplicity.h"
#include "AliTracker.h"

#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AliPhysicsSelection.h"
#include "AliTriggerAnalysis.h"

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"
#include "AlidNdPtAnalysis.h"

using namespace std;

ClassImp(AlidNdPtAnalysis)

//_____________________________________________________________________________
  AlidNdPtAnalysis::AlidNdPtAnalysis(): AlidNdPt(),
  fAnalysisFolder(0),
  fHistogramsOn(kFALSE),

  // event multiplicity correlation matrix 
  fEventMultCorrelationMatrix(0),

  // rec. track pt vs true track pt correlation matrix 
  fTrackPtCorrelationMatrix(0),

  // event level correction
  fGenEventMatrix(0),
  fGenSDEventMatrix(0),
  fGenDDEventMatrix(0),
  fGenNDEventMatrix(0),
  fGenNSDEventMatrix(0),

  fTriggerEventMatrix(0),
  fTriggerSDEventMatrix(0),
  fTriggerDDEventMatrix(0),
  fTriggerNDEventMatrix(0),
  fTriggerNSDEventMatrix(0),

  fRecEventMatrix(0),
  fRecSDEventMatrix(0),
  fRecDDEventMatrix(0),
  fRecNDEventMatrix(0),
  fRecNSDEventMatrix(0),

  //
  // track-event level correction 
  //
  fGenTrackEventMatrix(0),
  fGenTrackSDEventMatrix(0),
  fGenTrackDDEventMatrix(0),
  fGenTrackNDEventMatrix(0),
  fGenTrackNSDEventMatrix(0),

  fTriggerTrackEventMatrix(0),
  fTriggerTrackSDEventMatrix(0),
  fTriggerTrackDDEventMatrix(0),
  fTriggerTrackNDEventMatrix(0),
  fTriggerTrackNSDEventMatrix(0),

  fRecTrackEventMatrix(0),
  fRecTrackSDEventMatrix(0),
  fRecTrackDDEventMatrix(0),
  fRecTrackNDEventMatrix(0),
  fRecTrackNSDEventMatrix(0),

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
  fRecMCEventHist3(0),

  // rec. pt and eta resolution w.r.t MC
  fRecMCTrackHist1(0),

  //multple reconstructed tracks
  fMCMultRecTrackHist1(0), 

  // rec. track control histograms
  fRecTrackHist2(0),

  // Generic histograms to be corrected
  fRecEventHist(0),
  fRecTrackHist(0),

  // Candle event histogram
  fRecCandleEventMatrix(0),
  
  fMultNbins(0),
  fPtNbins(0),
  fPtCorrNbins(0),
  fEtaNbins(0),
  fZvNbins(0),
  fBinsMult(0),
  fBinsPt(0),
  fBinsPtCorr(0),
  fBinsEta(0),
  fBinsZv(0),

  fIsInit(kFALSE)  
  
{
  // default constructor
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    fMCTrackHist1[i]=0;     
    fMCPrimTrackHist1[i]=0;     
    fMCPrimTrackHist2[i]=0;     
    fMCSecTrackHist1[i]=0;     
    fRecTrackHist1[i]=0;     
    fRecTrackMultHist1[i]=0;     
  }
  //Init();
}

//_____________________________________________________________________________
AlidNdPtAnalysis::AlidNdPtAnalysis(Char_t* name, Char_t* title): AlidNdPt(name,title),
  fAnalysisFolder(0),
  fHistogramsOn(kFALSE),

  // event multiplicity correlation matrix 
  fEventMultCorrelationMatrix(0),

  // rec. track pt vs true track pt correlation matrix 
  fTrackPtCorrelationMatrix(0),

  // event level correction
  fGenEventMatrix(0),
  fGenSDEventMatrix(0),
  fGenDDEventMatrix(0),
  fGenNDEventMatrix(0),
  fGenNSDEventMatrix(0),

  fTriggerEventMatrix(0),
  fTriggerSDEventMatrix(0),
  fTriggerDDEventMatrix(0),
  fTriggerNDEventMatrix(0),
  fTriggerNSDEventMatrix(0),

  fRecEventMatrix(0),
  fRecSDEventMatrix(0),
  fRecDDEventMatrix(0),
  fRecNDEventMatrix(0),
  fRecNSDEventMatrix(0),

  //
  // track-event level correction 
  //
  fGenTrackEventMatrix(0),
  fGenTrackSDEventMatrix(0),
  fGenTrackDDEventMatrix(0),
  fGenTrackNDEventMatrix(0),
  fGenTrackNSDEventMatrix(0),

  fTriggerTrackEventMatrix(0),
  fTriggerTrackSDEventMatrix(0),
  fTriggerTrackDDEventMatrix(0),
  fTriggerTrackNDEventMatrix(0),
  fTriggerTrackNSDEventMatrix(0),

  fRecTrackEventMatrix(0),
  fRecTrackSDEventMatrix(0),
  fRecTrackDDEventMatrix(0),
  fRecTrackNDEventMatrix(0),
  fRecTrackNSDEventMatrix(0),

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
  fRecMCEventHist3(0),
 
  // rec. pt and eta resolution w.r.t MC
  fRecMCTrackHist1(0),

  //multple reconstructed tracks
  fMCMultRecTrackHist1(0), 

  // rec. track control histograms
  fRecTrackHist2(0),

  // Generic histograms to be corrected
  fRecEventHist(0),
  fRecTrackHist(0),

  // Candle event histogram
  fRecCandleEventMatrix(0),
  
  fMultNbins(0),
  fPtNbins(0),
  fPtCorrNbins(0),
  fEtaNbins(0),
  fZvNbins(0),
  fBinsMult(0),
  fBinsPt(0),
  fBinsPtCorr(0),
  fBinsEta(0),
  fBinsZv(0),

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
    fRecTrackMultHist1[i]=0; 
  }

  //Init();
}

//_____________________________________________________________________________
AlidNdPtAnalysis::~AlidNdPtAnalysis() {
  //
  // destructor
  //
  if(fEventMultCorrelationMatrix) delete fEventMultCorrelationMatrix; fEventMultCorrelationMatrix=0;
  //
  if(fTrackPtCorrelationMatrix) delete fTrackPtCorrelationMatrix; fTrackPtCorrelationMatrix=0;
  //
  if(fGenEventMatrix) delete fGenEventMatrix; fGenEventMatrix=0;
  if(fGenSDEventMatrix) delete fGenSDEventMatrix; fGenSDEventMatrix=0;
  if(fGenDDEventMatrix) delete fGenDDEventMatrix; fGenDDEventMatrix=0;
  if(fGenNDEventMatrix) delete fGenNDEventMatrix; fGenNDEventMatrix=0;
  if(fGenNSDEventMatrix) delete fGenNSDEventMatrix; fGenNSDEventMatrix=0;

  if(fTriggerEventMatrix) delete fTriggerEventMatrix; fTriggerEventMatrix=0;
  if(fTriggerSDEventMatrix) delete fTriggerSDEventMatrix; fTriggerSDEventMatrix=0;
  if(fTriggerDDEventMatrix) delete fTriggerDDEventMatrix; fTriggerDDEventMatrix=0;
  if(fTriggerNDEventMatrix) delete fTriggerNDEventMatrix; fTriggerNDEventMatrix=0;
  if(fTriggerNSDEventMatrix) delete fTriggerNSDEventMatrix; fTriggerNSDEventMatrix=0;

  if(fRecEventMatrix) delete fRecEventMatrix; fRecEventMatrix=0;
  if(fRecSDEventMatrix) delete fRecSDEventMatrix; fRecSDEventMatrix=0;
  if(fRecDDEventMatrix) delete fRecDDEventMatrix; fRecDDEventMatrix=0;
  if(fRecNDEventMatrix) delete fRecNDEventMatrix; fRecNDEventMatrix=0;
  if(fRecNSDEventMatrix) delete fRecNSDEventMatrix; fRecNSDEventMatrix=0;

  if(fRecCandleEventMatrix) delete fRecCandleEventMatrix; fRecCandleEventMatrix=0;
  //
  if(fGenTrackEventMatrix) delete fGenTrackEventMatrix; fGenTrackEventMatrix=0;
  if(fGenTrackSDEventMatrix) delete fGenTrackSDEventMatrix; fGenTrackSDEventMatrix=0;
  if(fGenTrackDDEventMatrix) delete fGenTrackDDEventMatrix; fGenTrackDDEventMatrix=0;
  if(fGenTrackNDEventMatrix) delete fGenTrackNDEventMatrix; fGenTrackNDEventMatrix=0;
  if(fGenTrackNSDEventMatrix) delete fGenTrackNSDEventMatrix; fGenTrackNSDEventMatrix=0;

  if(fTriggerTrackEventMatrix) delete fTriggerTrackEventMatrix; fTriggerTrackEventMatrix=0;
  if(fTriggerTrackSDEventMatrix) delete fTriggerTrackSDEventMatrix; fTriggerTrackSDEventMatrix=0;
  if(fTriggerTrackDDEventMatrix) delete fTriggerTrackDDEventMatrix; fTriggerTrackDDEventMatrix=0;
  if(fTriggerTrackNDEventMatrix) delete fTriggerTrackNDEventMatrix; fTriggerTrackNDEventMatrix=0;
  if(fTriggerTrackNSDEventMatrix) delete fTriggerTrackNSDEventMatrix; fTriggerTrackNSDEventMatrix=0;

  if(fRecTrackEventMatrix) delete fRecTrackEventMatrix; fRecTrackEventMatrix=0;
  if(fRecTrackSDEventMatrix) delete fRecTrackSDEventMatrix; fRecTrackSDEventMatrix=0;
  if(fRecTrackDDEventMatrix) delete fRecTrackDDEventMatrix; fRecTrackDDEventMatrix=0;
  if(fRecTrackNDEventMatrix) delete fRecTrackNDEventMatrix; fRecTrackNDEventMatrix=0;
  if(fRecTrackNSDEventMatrix) delete fRecTrackNSDEventMatrix; fRecTrackNSDEventMatrix=0;

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
  if(fRecMCEventHist3) delete fRecMCEventHist3; fRecMCEventHist3=0;
  //
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    if(fMCTrackHist1[i]) delete fMCTrackHist1[i]; fMCTrackHist1[i]=0;
    if(fMCPrimTrackHist1[i]) delete fMCPrimTrackHist1[i]; fMCPrimTrackHist1[i]=0;
    if(fMCPrimTrackHist2[i]) delete fMCPrimTrackHist2[i]; fMCPrimTrackHist2[i]=0;
    if(fMCSecTrackHist1[i]) delete fMCSecTrackHist1[i]; fMCSecTrackHist1[i]=0;
    if(fRecTrackHist1[i]) delete fRecTrackHist1[i]; fRecTrackHist1[i]=0;
    if(fRecTrackMultHist1[i]) delete fRecTrackMultHist1[i]; fRecTrackMultHist1[i]=0;
  }
  if(fRecMCTrackHist1) delete fRecMCTrackHist1; fRecMCTrackHist1=0;
  if(fMCMultRecTrackHist1) delete fMCMultRecTrackHist1; fMCMultRecTrackHist1=0; 
  if(fRecTrackHist2) delete fRecTrackHist2; fRecTrackHist2=0; 

  //
  if(fRecEventHist) delete fRecEventHist; fRecEventHist=0; 
  if(fRecTrackHist) delete fRecTrackHist; fRecTrackHist=0; 

  //
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AlidNdPtAnalysis::Init(){
  
  Int_t multNbins = 27;
  Int_t ptNbinsTrackEventCorr = 36;
  Int_t ptNbins = 68;
  Int_t etaNbins = 30;
  Int_t zvNbins = 12;

  Double_t binsMultDefault[28] = {-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,
                                     9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,
				     19.5,20.5, 21.5, 22.5, 23.5, 24.5, 29.5, 149.5};

  Double_t binsPtTrackEventCorrDefault[37] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,3.0,4.0,50.0};

  Double_t binsPtDefault[69] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};

  Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};

  Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
  
  Double_t* binsMult = binsMultDefault;
  Double_t* binsPtTrackEventCorr = binsPtTrackEventCorrDefault;
  Double_t* binsPt = binsPtDefault;
  Double_t* binsEta = binsEtaDefault;
  Double_t* binsZv = binsZvDefault;  

  if (fMultNbins > 0) { multNbins = fMultNbins; binsMult = fBinsMult; }
  if (fPtNbins > 0) { ptNbins  = fPtNbins; binsPt = fBinsPt; }
  if (fPtCorrNbins > 0) { ptNbinsTrackEventCorr  = fPtCorrNbins; binsPtTrackEventCorr = fBinsPtCorr;  }
  if (fEtaNbins > 0) { etaNbins  = fEtaNbins; binsEta = fBinsEta; }
  if (fZvNbins > 0) { zvNbins  = fZvNbins; binsZv = fBinsZv; }  

  //Int_t binsTrackMatrix[3]={zvNbins,ptNbins,etaNbins};
  Int_t binsTrackEventCorrMatrix[3]={zvNbins,ptNbinsTrackEventCorr,etaNbins};


  //
  // Generic histograms to be corrected
  //
  Int_t binsEventHist[2]={100,150};
  Double_t minEventHist[2]={-25.,-0.5}; 
  Double_t maxEventHist[2]={25.,149.5}; 

  fRecEventHist = new THnSparseF("fRecEventHist","Zv:multMB",2,binsEventHist,minEventHist,maxEventHist);
  fRecEventHist->GetAxis(0)->SetTitle("Zv (cm)");
  fRecEventHist->GetAxis(1)->SetTitle("multiplicity MB");
  fRecEventHist->Sumw2();

  //
  Int_t binsTrackHist[4]={100,ptNbins,etaNbins,150};
  Double_t minTrackHist[4]={-25.,0.,-1.5,-0.5}; 
  Double_t maxTrackHist[4]={25.,50.,1.5,149.5}; 

  fRecTrackHist = new THnSparseF("fRecTrackHist","Zv:pT:eta:multRec",4,binsTrackHist,minTrackHist,maxTrackHist);
  fRecTrackHist->SetBinEdges(1,binsPt);
  fRecTrackHist->SetBinEdges(2,binsEta);
  fRecTrackHist->GetAxis(0)->SetTitle("Zv (cm)");
  fRecTrackHist->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
  fRecTrackHist->GetAxis(2)->SetTitle("#eta");
  fRecTrackHist->GetAxis(3)->SetTitle("multiplicity MB");
  fRecTrackHist->Sumw2();

  //
  // rec. vs MC correlation matrices
  //
  Int_t binsMultTrueEventMatrix[3]={150,150,150};
  Double_t minMultTrueEventMatrix[3]={-0.5,-0.5,-0.5}; 
  Double_t maxMultTrueEventMatrix[3]={149.5,149.5,149.5}; 
  fEventMultCorrelationMatrix = new THnSparseF("fEventMultCorrelationMatrix","mult:true_mult:multMB",3,binsMultTrueEventMatrix,minMultTrueEventMatrix,maxMultTrueEventMatrix);
  fEventMultCorrelationMatrix->GetAxis(0)->SetTitle("track multiplicity");
  fEventMultCorrelationMatrix->GetAxis(1)->SetTitle("true multiplicity");
  fEventMultCorrelationMatrix->GetAxis(2)->SetTitle("MB multiplicity");
  fEventMultCorrelationMatrix->Sumw2();
  
  Int_t binsTrackPtCorrelationMatrix[3]={ptNbins,ptNbins,etaNbins};
  fTrackPtCorrelationMatrix = new THnSparseF("fTrackPtCorrelationMatrix","Pt:mcPt:mcEta",3,binsTrackPtCorrelationMatrix);
  fTrackPtCorrelationMatrix->SetBinEdges(0,binsPt);
  fTrackPtCorrelationMatrix->SetBinEdges(1,binsPt);
  fTrackPtCorrelationMatrix->SetBinEdges(2,binsEta);
  fTrackPtCorrelationMatrix->GetAxis(0)->SetTitle("Pt (GeV/c)");
  fTrackPtCorrelationMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fTrackPtCorrelationMatrix->GetAxis(2)->SetTitle("mcEta");
  fTrackPtCorrelationMatrix->Sumw2();

  //
  // Efficiency and contamination correction matrices
  //
  Int_t binsEventMatrix[2]={zvNbins,multNbins};
  Double_t minEventMatrix[2]={-25.,-0.5}; 
  Double_t maxEventMatrix[2]={25.,149.5}; 

  fGenEventMatrix = new THnSparseF("fGenEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fGenEventMatrix->SetBinEdges(0,binsZv);
  fGenEventMatrix->SetBinEdges(1,binsMult);
  fGenEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fGenEventMatrix->Sumw2();
  
  fGenSDEventMatrix = new THnSparseF("fGenSDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fGenSDEventMatrix->SetBinEdges(0,binsZv);
  fGenSDEventMatrix->SetBinEdges(1,binsMult);
  fGenSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenSDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fGenSDEventMatrix->Sumw2();
  
  fGenDDEventMatrix = new THnSparseF("fGenDDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fGenDDEventMatrix->SetBinEdges(0,binsZv);
  fGenDDEventMatrix->SetBinEdges(1,binsMult);
  fGenDDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenDDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fGenDDEventMatrix->Sumw2();
  
  fGenNDEventMatrix = new THnSparseF("fGenNDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fGenNDEventMatrix->SetBinEdges(0,binsZv);
  fGenNDEventMatrix->SetBinEdges(1,binsMult);
  fGenNDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenNDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fGenNDEventMatrix->Sumw2();

  fGenNSDEventMatrix = new THnSparseF("fGenNSDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fGenNSDEventMatrix->SetBinEdges(0,binsZv);
  fGenNSDEventMatrix->SetBinEdges(1,binsMult);
  fGenNSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenNSDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fGenNSDEventMatrix->Sumw2();

  //
  fTriggerEventMatrix = new THnSparseF("fTriggerEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fTriggerEventMatrix->SetBinEdges(0,binsZv);
  fTriggerEventMatrix->SetBinEdges(1,binsMult);
  fTriggerEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fTriggerEventMatrix->Sumw2();

  fTriggerSDEventMatrix = new THnSparseF("fTriggerSDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fTriggerSDEventMatrix->SetBinEdges(0,binsZv);
  fTriggerSDEventMatrix->SetBinEdges(1,binsMult);
  fTriggerSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerSDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fTriggerSDEventMatrix->Sumw2();
  
  fTriggerDDEventMatrix = new THnSparseF("fTriggerDDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fTriggerDDEventMatrix->SetBinEdges(0,binsZv);
  fTriggerDDEventMatrix->SetBinEdges(1,binsMult);
  fTriggerDDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerDDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fTriggerDDEventMatrix->Sumw2();
  
  fTriggerNDEventMatrix = new THnSparseF("fTriggerNDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fTriggerNDEventMatrix->SetBinEdges(0,binsZv);
  fTriggerNDEventMatrix->SetBinEdges(1,binsMult);
  fTriggerNDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerNDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fTriggerNDEventMatrix->Sumw2();
 
  fTriggerNSDEventMatrix = new THnSparseF("fTriggerNSDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fTriggerNSDEventMatrix->SetBinEdges(0,binsZv);
  fTriggerNSDEventMatrix->SetBinEdges(1,binsMult);
  fTriggerNSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerNSDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fTriggerNSDEventMatrix->Sumw2();
 
  //
  fRecEventMatrix = new THnSparseF("fRecEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fRecEventMatrix->SetBinEdges(0,binsZv);
  fRecEventMatrix->SetBinEdges(1,binsMult);
  fRecEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fRecEventMatrix->Sumw2();

  fRecSDEventMatrix = new THnSparseF("fRecSDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fRecSDEventMatrix->SetBinEdges(0,binsZv);
  fRecSDEventMatrix->SetBinEdges(1,binsMult);
  fRecSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecSDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fRecSDEventMatrix->Sumw2();
  
  fRecDDEventMatrix = new THnSparseF("fRecDDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fRecDDEventMatrix->SetBinEdges(0,binsZv);
  fRecDDEventMatrix->SetBinEdges(1,binsMult);
  fRecDDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecDDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fRecDDEventMatrix->Sumw2();
  
  fRecNDEventMatrix = new THnSparseF("fRecNDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fRecNDEventMatrix->SetBinEdges(0,binsZv);
  fRecNDEventMatrix->SetBinEdges(1,binsMult);
  fRecNDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecNDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fRecNDEventMatrix->Sumw2();
 
  fRecNSDEventMatrix = new THnSparseF("fRecNSDEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fRecNSDEventMatrix->SetBinEdges(0,binsZv);
  fRecNSDEventMatrix->SetBinEdges(1,binsMult);
  fRecNSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecNSDEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fRecNSDEventMatrix->Sumw2();

  fRecCandleEventMatrix = new THnSparseF("fRecCandleEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fRecCandleEventMatrix->SetBinEdges(0,binsZv);
  fRecCandleEventMatrix->SetBinEdges(1,binsMult);
  fRecCandleEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecCandleEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fRecCandleEventMatrix->Sumw2();

  // 
  // track to event corrections
  //

  fGenTrackEventMatrix = new THnSparseF("fGenTrackEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenTrackEventMatrix->SetBinEdges(0,binsZv);
  fGenTrackEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fGenTrackEventMatrix->SetBinEdges(2,binsEta);
  fGenTrackEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackEventMatrix->Sumw2();

  fGenTrackSDEventMatrix = new THnSparseF("fGenTrackSDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenTrackSDEventMatrix->SetBinEdges(0,binsZv);
  fGenTrackSDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fGenTrackSDEventMatrix->SetBinEdges(2,binsEta);
  fGenTrackSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackSDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackSDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackSDEventMatrix->Sumw2();

  fGenTrackDDEventMatrix = new THnSparseF("fGenTrackDDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenTrackDDEventMatrix->SetBinEdges(0,binsZv);
  fGenTrackDDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fGenTrackDDEventMatrix->SetBinEdges(2,binsEta);
  fGenTrackDDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackDDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackDDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackDDEventMatrix->Sumw2();

  fGenTrackNDEventMatrix = new THnSparseF("fGenTrackNDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenTrackNDEventMatrix->SetBinEdges(0,binsZv);
  fGenTrackNDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fGenTrackNDEventMatrix->SetBinEdges(2,binsEta);
  fGenTrackNDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackNDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackNDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackNDEventMatrix->Sumw2();

  fGenTrackNSDEventMatrix = new THnSparseF("fGenTrackNSDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenTrackNSDEventMatrix->SetBinEdges(0,binsZv);
  fGenTrackNSDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fGenTrackNSDEventMatrix->SetBinEdges(2,binsEta);
  fGenTrackNSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackNSDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackNSDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackNSDEventMatrix->Sumw2();


  //
  fTriggerTrackEventMatrix = new THnSparseF("fTriggerTrackEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fTriggerTrackEventMatrix->SetBinEdges(0,binsZv);
  fTriggerTrackEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fTriggerTrackEventMatrix->SetBinEdges(2,binsEta);
  fTriggerTrackEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerTrackEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fTriggerTrackEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fTriggerTrackEventMatrix->Sumw2();

  fTriggerTrackSDEventMatrix = new THnSparseF("fTriggerTrackSDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fTriggerTrackSDEventMatrix->SetBinEdges(0,binsZv);
  fTriggerTrackSDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fTriggerTrackSDEventMatrix->SetBinEdges(2,binsEta);
  fTriggerTrackSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerTrackSDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fTriggerTrackSDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fTriggerTrackSDEventMatrix->Sumw2();

  fTriggerTrackDDEventMatrix = new THnSparseF("fTriggerTrackDDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fTriggerTrackDDEventMatrix->SetBinEdges(0,binsZv);
  fTriggerTrackDDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fTriggerTrackDDEventMatrix->SetBinEdges(2,binsEta);
  fTriggerTrackDDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerTrackDDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fTriggerTrackDDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fTriggerTrackDDEventMatrix->Sumw2();

  fTriggerTrackNDEventMatrix = new THnSparseF("fTriggerTrackNDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fTriggerTrackNDEventMatrix->SetBinEdges(0,binsZv);
  fTriggerTrackNDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fTriggerTrackNDEventMatrix->SetBinEdges(2,binsEta);
  fTriggerTrackNDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerTrackNDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fTriggerTrackNDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fTriggerTrackNDEventMatrix->Sumw2();

  fTriggerTrackNSDEventMatrix = new THnSparseF("fTriggerTrackNSDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fTriggerTrackNSDEventMatrix->SetBinEdges(0,binsZv);
  fTriggerTrackNSDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fTriggerTrackNSDEventMatrix->SetBinEdges(2,binsEta);
  fTriggerTrackNSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerTrackNSDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fTriggerTrackNSDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fTriggerTrackNSDEventMatrix->Sumw2();

  //
  fRecTrackEventMatrix = new THnSparseF("fRecTrackEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecTrackEventMatrix->SetBinEdges(0,binsZv);
  fRecTrackEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fRecTrackEventMatrix->SetBinEdges(2,binsEta);
  fRecTrackEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecTrackEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecTrackEventMatrix->Sumw2();

  fRecTrackSDEventMatrix = new THnSparseF("fRecTrackSDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecTrackSDEventMatrix->SetBinEdges(0,binsZv);
  fRecTrackSDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fRecTrackSDEventMatrix->SetBinEdges(2,binsEta);
  fRecTrackSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackSDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecTrackSDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecTrackSDEventMatrix->Sumw2();

  fRecTrackDDEventMatrix = new THnSparseF("fRecTrackDDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecTrackDDEventMatrix->SetBinEdges(0,binsZv);
  fRecTrackDDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fRecTrackDDEventMatrix->SetBinEdges(2,binsEta);
  fRecTrackDDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackDDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecTrackDDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecTrackDDEventMatrix->Sumw2();

  fRecTrackNDEventMatrix = new THnSparseF("fRecTrackNDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecTrackNDEventMatrix->SetBinEdges(0,binsZv);
  fRecTrackNDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fRecTrackNDEventMatrix->SetBinEdges(2,binsEta);
  fRecTrackNDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackNDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecTrackNDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecTrackNDEventMatrix->Sumw2();

  fRecTrackNSDEventMatrix = new THnSparseF("fRecTrackNSDEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecTrackNSDEventMatrix->SetBinEdges(0,binsZv);
  fRecTrackNSDEventMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fRecTrackNSDEventMatrix->SetBinEdges(2,binsEta);
  fRecTrackNSDEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackNSDEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecTrackNSDEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecTrackNSDEventMatrix->Sumw2();

  //
  // tracks correction matrices
  //
  //fGenTrackMatrix = new THnSparseF("fGenTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackMatrix);
  fGenTrackMatrix = new THnSparseF("fGenTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenTrackMatrix->SetBinEdges(0,binsZv);
  //fGenTrackMatrix->SetBinEdges(1,binsPt);
  fGenTrackMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fGenTrackMatrix->SetBinEdges(2,binsEta);
  fGenTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackMatrix->Sumw2();

  //fGenPrimTrackMatrix = new THnSparseF("fGenPrimTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackMatrix);
  fGenPrimTrackMatrix = new THnSparseF("fGenPrimTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenPrimTrackMatrix->SetBinEdges(0,binsZv);
  //fGenPrimTrackMatrix->SetBinEdges(1,binsPt);
  fGenPrimTrackMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fGenPrimTrackMatrix->SetBinEdges(2,binsEta);
  fGenPrimTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenPrimTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenPrimTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenPrimTrackMatrix->Sumw2();

  //fRecPrimTrackMatrix = new THnSparseF("fRecPrimTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackMatrix);
  fRecPrimTrackMatrix = new THnSparseF("fRecPrimTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecPrimTrackMatrix->SetBinEdges(0,binsZv);
  //fRecPrimTrackMatrix->SetBinEdges(1,binsPt);
  fRecPrimTrackMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fRecPrimTrackMatrix->SetBinEdges(2,binsEta);
  fRecPrimTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecPrimTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecPrimTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecPrimTrackMatrix->Sumw2();

  //
  //fRecTrackMatrix = new THnSparseF("fRecTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackMatrix);
  fRecTrackMatrix = new THnSparseF("fRecTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecTrackMatrix->SetBinEdges(0,binsZv);
  //fRecTrackMatrix->SetBinEdges(1,binsPt);
  fRecTrackMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fRecTrackMatrix->SetBinEdges(2,binsEta);
  fRecTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackMatrix->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fRecTrackMatrix->GetAxis(2)->SetTitle("Eta");
  fRecTrackMatrix->Sumw2();

  //fRecSecTrackMatrix = new THnSparseF("fRecSecTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackMatrix);
  fRecSecTrackMatrix = new THnSparseF("fRecSecTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecSecTrackMatrix->SetBinEdges(0,binsZv);
  //fRecSecTrackMatrix->SetBinEdges(1,binsPt);
  fRecSecTrackMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fRecSecTrackMatrix->SetBinEdges(2,binsEta);
  fRecSecTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecSecTrackMatrix->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fRecSecTrackMatrix->GetAxis(2)->SetTitle("Eta");
  fRecSecTrackMatrix->Sumw2();

  //
  //fRecMultTrackMatrix = new THnSparseF("fRecMultTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackMatrix);
  fRecMultTrackMatrix = new THnSparseF("fRecMultTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecMultTrackMatrix->SetBinEdges(0,binsZv);
  //fRecMultTrackMatrix->SetBinEdges(1,binsPt);
  fRecMultTrackMatrix->SetBinEdges(1,binsPtTrackEventCorr);
  fRecMultTrackMatrix->SetBinEdges(2,binsEta);
  fRecMultTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecMultTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecMultTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecMultTrackMatrix->Sumw2();

  //
  // Control analysis histograms
  //

  Int_t binsMCEventHist1[3]={100,100,140};
  Double_t minMCEventHist1[3]={-0.1,-0.1,-35.}; 
  Double_t maxMCEventHist1[3]={0.1,0.1,35.}; 
  fMCEventHist1 = new THnSparseF("fMCEventHist1","mcXv:mcYv:mcZv",3,binsMCEventHist1,minMCEventHist1,maxMCEventHist1);
  fMCEventHist1->GetAxis(0)->SetTitle("mcXv (cm)");
  fMCEventHist1->GetAxis(1)->SetTitle("mcYv (cm)");
  fMCEventHist1->GetAxis(2)->SetTitle("mcZv (cm)");
  fMCEventHist1->Sumw2();

  //
  Int_t binsRecEventHist1[3]={100,100,140};
  Double_t minRecEventHist1[3]={-3.,-3.,-35.}; 
  Double_t maxRecEventHist1[3]={3.,3.,35.}; 
  
  fRecEventHist1 = new THnSparseF("fRecEventHist1","Xv:Yv:Zv",3,binsRecEventHist1,minRecEventHist1,maxRecEventHist1);
  fRecEventHist1->GetAxis(0)->SetTitle("Xv (cm)");
  fRecEventHist1->GetAxis(1)->SetTitle("Yv (cm)");
  fRecEventHist1->GetAxis(2)->SetTitle("Zv (cm)");
  fRecEventHist1->Sumw2();

  //
  Int_t binsRecEventHist2[3]={zvNbins,150,150};
  Double_t minRecEventHist2[3]={-25.,-0.5,-0.5}; 
  Double_t maxRecEventHist2[3]={25.,149.5,149.5}; 
  
  fRecEventHist2 = new THnSparseF("fRecEventHist2","Zv:multMB:mult",3,binsRecEventHist2,minRecEventHist2,maxRecEventHist2);
  fRecEventHist2->SetBinEdges(0,binsZv);
  fRecEventHist2->GetAxis(0)->SetTitle("Zv (cm)");
  fRecEventHist2->GetAxis(1)->SetTitle("multiplicity MB");
  fRecEventHist2->GetAxis(2)->SetTitle("multiplicity");
  fRecEventHist2->Sumw2();

  //
  Double_t kFact = 0.1;
  Int_t binsRecMCEventHist1[3]={100,100,100};
  Double_t minRecMCEventHist1[3]={-10.0*kFact,-10.0*kFact,-10.0*kFact}; 
  Double_t maxRecMCEventHist1[3]={10.0*kFact,10.0*kFact,10.0*kFact}; 
   
  fRecMCEventHist1 = new THnSparseF("fRecMCEventHist1","Xv-mcXv:Yv-mcYv:Zv-mcZv",3,binsRecMCEventHist1,minRecMCEventHist1,maxRecMCEventHist1);
  fRecMCEventHist1->GetAxis(0)->SetTitle("Xv-mcXv (cm)");
  fRecMCEventHist1->GetAxis(1)->SetTitle("Yv-mcYv (cm)");
  fRecMCEventHist1->GetAxis(2)->SetTitle("Zv-mcZv (cm)");
  fRecMCEventHist1->Sumw2();

  //
  Int_t binsRecMCEventHist2[3]={100,100,150};
  Double_t minRecMCEventHist2[3]={-10.0*kFact,-10.0*kFact,0.0}; 
  Double_t maxRecMCEventHist2[3]={10.0*kFact,10.0*kFact,149.50}; 

  fRecMCEventHist2 = new THnSparseF("fRecMCEventHist2","Xv-mcXv:Zv-mcZv:mult",3,binsRecMCEventHist2,minRecMCEventHist2,maxRecMCEventHist2);
  fRecMCEventHist2->GetAxis(0)->SetTitle("Xv-mcXv (cm)");
  fRecMCEventHist2->GetAxis(1)->SetTitle("Zv-mcZv (cm)");
  fRecMCEventHist2->GetAxis(2)->SetTitle("multiplicity");
  fRecMCEventHist2->Sumw2();

  Int_t binsRecMCEventHist3[2]={150,5};
  Double_t minRecMCEventHist3[2]={-0.5,0.0}; 
  Double_t maxRecMCEventHist3[2]={149.50,5.0}; 
  fRecMCEventHist3 = new THnSparseF("fRecMCEventHist3","mult:EventType (ND, DD, SD)",2,binsRecMCEventHist3,minRecMCEventHist3,maxRecMCEventHist3);
  fRecMCEventHist3->GetAxis(0)->SetTitle("multiplicity");
  fRecMCEventHist3->GetAxis(1)->SetTitle("EventType");
  fRecMCEventHist3->Sumw2();

  //
  char name[256];
  char title[256];
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) 
  {
  // THnSparse track histograms
 
  Int_t binsMCTrackHist1[3]=  {ptNbins, etaNbins, 90};
  Double_t minMCTrackHist1[3]={0.,-1.,0.}; 
  Double_t maxMCTrackHist1[3]={10.,1.,2.*TMath::Pi()}; 
  sprintf(name,"fMCTrackHist1_%d",i);
  sprintf(title,"mcPt:mcEta:mcPhi");
  
  fMCTrackHist1[i] = new THnSparseF(name,title,3,binsMCTrackHist1,minMCTrackHist1,maxMCTrackHist1);
  fMCTrackHist1[i]->SetBinEdges(0,binsPt);
  fMCTrackHist1[i]->SetBinEdges(1,binsEta);
  fMCTrackHist1[i]->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCTrackHist1[i]->GetAxis(1)->SetTitle("mcEta");
  fMCTrackHist1[i]->GetAxis(2)->SetTitle("mcPhi (rad)");
  fMCTrackHist1[i]->Sumw2();

  Int_t binsMCPrimTrackHist1[5]=  {ptNbins,etaNbins,6,20,4000};
  Double_t minMCPrimTrackHist1[5]={0.,-1.,0.,0.,0.}; 
  Double_t maxMCPrimTrackHist1[5]={10.,1.,6.,20.,4000.}; 
  sprintf(name,"fMCPrimTrackHist1_%d",i);
  sprintf(title,"mcPt:mcEta:pid:mech:mother");
  
  fMCPrimTrackHist1[i] = new THnSparseF(name,title,5,binsMCPrimTrackHist1,minMCPrimTrackHist1,maxMCPrimTrackHist1);
  fMCPrimTrackHist1[i]->SetBinEdges(0,binsPt);
  fMCPrimTrackHist1[i]->SetBinEdges(1,binsEta);
  fMCPrimTrackHist1[i]->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCPrimTrackHist1[i]->GetAxis(1)->SetTitle("mcEta");
  fMCPrimTrackHist1[i]->GetAxis(2)->SetTitle("pid");
  fMCPrimTrackHist1[i]->GetAxis(3)->SetTitle("mech");
  fMCPrimTrackHist1[i]->GetAxis(4)->SetTitle("mother");
  fMCPrimTrackHist1[i]->Sumw2();

  Int_t binsMCPrimTrackHist2[3]=  {4000,20,4000};
  Double_t minMCPrimTrackHist2[3]={0.,0.,0.}; 
  Double_t maxMCPrimTrackHist2[3]={4000.,20.,4000.}; 
  sprintf(name,"fMCPrimTrackHist2_%d",i);
  sprintf(title,"pdg:mech:mother");
  
  fMCPrimTrackHist2[i] = new THnSparseF(name,title,3,binsMCPrimTrackHist2,minMCPrimTrackHist2,maxMCPrimTrackHist2);
  fMCPrimTrackHist2[i]->GetAxis(0)->SetTitle("pdg");
  fMCPrimTrackHist2[i]->GetAxis(1)->SetTitle("mech");
  fMCPrimTrackHist2[i]->GetAxis(2)->SetTitle("mother");
  fMCPrimTrackHist2[i]->Sumw2();

  Int_t binsMCSecTrackHist1[5]=  {ptNbins,etaNbins,6,20,4000};
  Double_t minMCSecTrackHist1[5]={0.,-1.,0.,0.,0.}; 
  Double_t maxMCSecTrackHist1[5]={10.,1.,6.,20.,4000.}; 
  sprintf(name,"fMCSecTrackHist1_%d",i);
  sprintf(title,"mcPt:mcEta:pid:mech:mother");
  
  fMCSecTrackHist1[i] = new THnSparseF(name,title,5,binsMCSecTrackHist1,minMCSecTrackHist1,maxMCSecTrackHist1);
  fMCSecTrackHist1[i]->SetBinEdges(0,binsPt);
  fMCSecTrackHist1[i]->SetBinEdges(1,binsEta);
  fMCSecTrackHist1[i]->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCSecTrackHist1[i]->GetAxis(1)->SetTitle("mcEta");
  fMCSecTrackHist1[i]->GetAxis(2)->SetTitle("pid");
  fMCSecTrackHist1[i]->GetAxis(3)->SetTitle("mech");
  fMCSecTrackHist1[i]->GetAxis(4)->SetTitle("mother");
  fMCSecTrackHist1[i]->Sumw2();

  //

  // 
  
  Int_t binsRecTrackHist1[3]={ptNbins,etaNbins,90};
  Double_t minRecTrackHist1[3]={0.,-1.,0.}; 
  Double_t maxRecTrackHist1[3]={10.,1.,2.*TMath::Pi()};
  sprintf(name,"fRecTrackHist1_%d",i);
  sprintf(title,"Pt:Eta:Phi");
  fRecTrackHist1[i] = new THnSparseF(name,title,3,binsRecTrackHist1,minRecTrackHist1,maxRecTrackHist1);
  fRecTrackHist1[i]->SetBinEdges(0,binsPt);
  fRecTrackHist1[i]->SetBinEdges(1,binsEta);
  fRecTrackHist1[i]->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fRecTrackHist1[i]->GetAxis(1)->SetTitle("#eta");
  fRecTrackHist1[i]->GetAxis(2)->SetTitle("#phi (rad)");
  fRecTrackHist1[i]->Sumw2();

  // 
  Int_t binsRecTrackMultHist1[2]={ptNbins,150};
  Double_t minRecTrackMultHist1[2]={0.,-0.5}; 
  Double_t maxRecTrackMultHist1[2]={10.,149.5};
  sprintf(name,"fRecTrackMultHist_%d",i);
  sprintf(title,"Pt:Mult");
  fRecTrackMultHist1[i] = new THnSparseF(name,title,2,binsRecTrackMultHist1,minRecTrackMultHist1,maxRecTrackMultHist1);
  fRecTrackMultHist1[i]->SetBinEdges(0,binsPt);
  fRecTrackMultHist1[i]->GetAxis(0)->SetTitle("Pt (GeV/c)");
  fRecTrackMultHist1[i]->GetAxis(1)->SetTitle("multiplicity");
  fRecTrackMultHist1[i]->Sumw2();
  }

  Int_t binsRecMCTrackHist1[4] = {ptNbins,etaNbins,100,100};
  Double_t minRecMCTrackHist1[4]={0.,-1.,-0.5,-0.5}; 
  Double_t maxRecMCTrackHist1[4]={20.,1.,0.5,0.5}; 
  sprintf(name,"fRecMCTrackHist1");
  sprintf(title,"mcPt:mcEta:(Pt-mcPt)/mcPt:(Eta-mcEta)");
  fRecMCTrackHist1 = new THnSparseF(name,title,4,binsRecMCTrackHist1,minRecMCTrackHist1,maxRecMCTrackHist1);
  fRecMCTrackHist1->SetBinEdges(0,binsPt);
  fRecMCTrackHist1->SetBinEdges(1,binsEta);
  fRecMCTrackHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fRecMCTrackHist1->GetAxis(1)->SetTitle("mcEta");
  fRecMCTrackHist1->GetAxis(2)->SetTitle("(Pt-mcPt)/mcPt");
  fRecMCTrackHist1->GetAxis(3)->SetTitle("Eta-mcEta");

  Int_t binsMCMultRecTrackHist1[3] = {ptNbins,etaNbins,6};
  Double_t minMCMultRecTrackHist1[3]={0.,-1.,0.}; 
  Double_t maxMCMultRecTrackHist1[3]={20.,1.,6.}; 
  sprintf(name,"fMCMultRecTrackHist1");
  sprintf(title,"mcPt:mcEta:pid");
  fMCMultRecTrackHist1 = new THnSparseF(name,title,3,binsMCMultRecTrackHist1,minMCMultRecTrackHist1,maxMCMultRecTrackHist1);
  fMCMultRecTrackHist1->SetBinEdges(0,binsPt);
  fMCMultRecTrackHist1->SetBinEdges(1,binsEta);
  fMCMultRecTrackHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCMultRecTrackHist1->GetAxis(1)->SetTitle("mcEta");
  fMCMultRecTrackHist1->GetAxis(2)->SetTitle("pid");

  //nClust:chi2PerClust:pt:eta:phi
  Int_t binsRecTrackHist2[5]={160,100,ptNbins,etaNbins,90};
  Double_t minRecTrackHist2[5]={0., 0., 0, -1.5, 0.};
  Double_t maxRecRecTrackHist2[5]={160.,10., 16, 1.5, 2.*TMath::Pi()};

  fRecTrackHist2 = new THnSparseF("fRecTrackHist2","nClust:chi2PerClust:pt:eta:phi",5,binsRecTrackHist2,minRecTrackHist2,maxRecRecTrackHist2);
  fRecTrackHist2->SetBinEdges(2,binsPt);
  fRecTrackHist2->SetBinEdges(3,binsEta);
  fRecTrackHist2->GetAxis(0)->SetTitle("nClust");
  fRecTrackHist2->GetAxis(1)->SetTitle("chi2PerClust");
  fRecTrackHist2->GetAxis(2)->SetTitle("p_{T} (GeV/c)");
  fRecTrackHist2->GetAxis(3)->SetTitle("#eta");
  fRecTrackHist2->GetAxis(4)->SetTitle("#phi (rad)");
  fRecTrackHist2->Sumw2();

  // init folder
  fAnalysisFolder = CreateFolder("folderdNdPt","Analysis dNdPt Folder");
  
  // set init flag
  fIsInit = kTRUE;
}

//_____________________________________________________________________________
void AlidNdPtAnalysis::Process(AliESDEvent *const esdEvent, AliMCEvent *const mcEvent)
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
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts  || !esdTrackCuts) {
    AliDebug(AliLog::kError, "cuts not available");
    return;
  }

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
    isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;
    //SetPhysicsTriggerSelection(physicsSelection);

    if (isEventTriggered) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }

    // calculate LHC background
    if(!IsUseMCInfo()) 
    { 
      //
      // 0-multiplicity bin for LHC background correction
      //
      if( GetAnalysisMode() == AlidNdPtHelper::kTPCITS || 
          GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtx || 
	  GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtxUpdate || 
          GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || 
	  GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt ) 
      {
        physicsSelection->SetBin0CallbackViaPointer(&AlidNdPtAnalysis::IsBinZeroTrackSPDvtx);
      } else {
        physicsSelection->SetBin0CallbackViaPointer(&AlidNdPtAnalysis::IsBinZeroSPDvtx);
      }
    }
  }


  /*
  if(evtCuts->IsTriggerRequired())  
  {
    //
    trigSel = GetPhysicsTriggerSelection();
    if(!trigSel) {
      AliDebug(AliLog::kError, "cannot get trigSel");
      return;
    }
    
    if(IsUseMCInfo()) 
    { 
      trigSel->SetAnalyzeMC();

      if(GetParticleMode() == AlidNdPtHelper::kVZEROCase1)
      {
        // check V0 systematics (case1)
	// Initialization done in the macro
        trigAna = trigSel->GetTriggerAnalysis();
        if(!trigAna) 
          return;

        //trigAna->SetV0HwPars(15, 61.5, 86.5);
        //trigAna->SetV0AdcThr(15);

        isEventTriggered = trigSel->IsCollisionCandidate(esdEvent);
        //printf("MB1 & kVZEROCase1  %d \n",isEventTriggered);
        //isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
	
        if(GetTrigger() == AliTriggerAnalysis::kV0AND) 
	{
          isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
          //printf("V0AND %d \n",isEventTriggered);
        }
      }
      else if(GetParticleMode() == AlidNdPtHelper::kVZEROCase2)
      {
        // check V0 systematics (case2 only in MC)
	// Initialization done in the macro

        trigAna = trigSel->GetTriggerAnalysis();
        if(!trigAna) 
          return;

        //trigAna->SetV0HwPars(0, 0, 125);
        //trigAna->SetV0AdcThr(0);

        isEventTriggered = trigSel->IsCollisionCandidate(esdEvent);
        //isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
	
	if(GetTrigger() == AliTriggerAnalysis::kV0AND) 
	{
          isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
          //printf("V0AND %d \n",isEventTriggered);
        }
      }
      else {
        isEventTriggered = trigSel->IsCollisionCandidate(esdEvent);
        //printf("MB1 %d \n",isEventTriggered);
	
        if(GetTrigger() == AliTriggerAnalysis::kV0AND) 
	{
          trigAna = trigSel->GetTriggerAnalysis();
          if(!trigAna) 
            return;

          isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
          //printf("V0AND %d \n",isEventTriggered);
        }
      }
    }
    else {
      //
      // 0-multiplicity bin for LHC background correction
      //
      if( GetAnalysisMode() == AlidNdPtHelper::kTPCITS || 
          GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtx || 
	  GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtxUpdate || 
          GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || 
	  GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt ) 
      {
        trigSel->SetBin0CallbackViaPointer(&AlidNdPtAnalysis::IsBinZeroTrackSPDvtx);
      } else {
        trigSel->SetBin0CallbackViaPointer(&AlidNdPtAnalysis::IsBinZeroSPDvtx);
      }

      if(GetParticleMode() == AlidNdPtHelper::kVZEROCase1)
      {
        // check V0 systematics (case1)
	// Initialization done in the macro
        trigAna = trigSel->GetTriggerAnalysis();
        if(!trigAna) 
          return;

        //trigAna->SetV0HwPars(15, 61.5, 86.5);
        //trigAna->SetV0AdcThr(15);

        isEventTriggered = trigSel->IsCollisionCandidate(esdEvent);
        //isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
	
        if(GetTrigger() == AliTriggerAnalysis::kV0AND) 
	{
          isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
          //printf("V0AND %d \n",isEventTriggered);
        }
      }
      else {
        isEventTriggered = trigSel->IsCollisionCandidate(esdEvent);
        //printf("MB1 %d \n",isEventTriggered);
	
        if(GetTrigger() == AliTriggerAnalysis::kV0AND) 
	{
          trigAna = trigSel->GetTriggerAnalysis();
          if(!trigAna) 
            return;

          isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
          //printf("V0AND %d \n",isEventTriggered);
        }
      }
    }
  }
  */


  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);
  AliPWG0Helper::MCProcessType evtType = AliPWG0Helper::kInvalidProcess;

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
    // get event type (ND=0x1, DD=0x2, SD=0x4)
    evtType = AliPWG0Helper::GetEventProcessType(header);
    AliDebug(AliLog::kDebug+1, Form("Found process type %d", evtType));

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    Double_t vMCEventHist1[3]={vtxMC[0],vtxMC[1],vtxMC[2]};
    fMCEventHist1->Fill(vMCEventHist1);

    // multipliticy of all MC primary tracks
    // in Zv, pt and eta ranges)
    multMCTrueTracks = AlidNdPtHelper::GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);

  } // end bUseMC


  // get reconstructed vertex  
  const AliESDVertex* vtxESD = 0; 
  Bool_t isRecVertex = kFALSE;
  if(evtCuts->IsRecVertexRequired()) 
  {
    Bool_t bRedoTPCVertex = evtCuts->IsRedoTPCVertex();
    Bool_t bUseConstraints = evtCuts->IsUseBeamSpotConstraint();
    vtxESD = AlidNdPtHelper::GetVertex(esdEvent,evtCuts,accCuts,esdTrackCuts,GetAnalysisMode(),kFALSE,bRedoTPCVertex,bUseConstraints); 
    isRecVertex = AlidNdPtHelper::TestRecVertex(vtxESD, esdEvent->GetPrimaryVertexSPD(), GetAnalysisMode(), kFALSE);
  }

  if( IsUseMCInfo() && !evtCuts->IsRecVertexRequired() ) {
    vtxESD = new AliESDVertex(vtxMC[2],10.,genHeader->NProduced(),"smearMC");
    isRecVertex = kTRUE;
  }

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD) && isRecVertex; 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());

  // vertex contributors
  Int_t multMBTracks = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) 
  {  
     if(vtxESD->GetStatus() && isRecVertex)
       multMBTracks = AlidNdPtHelper::GetTPCMBTrackMult(esdEvent,evtCuts,accCuts,esdTrackCuts);
  } 
  else if( GetAnalysisMode() == AlidNdPtHelper::kTPCSPDvtx || 
           GetAnalysisMode()==AlidNdPtHelper::kTPCSPDvtxUpdate ) 
  {
     const AliMultiplicity* mult = esdEvent->GetMultiplicity();
     if(mult && vtxESD->GetStatus() && isRecVertex)
       multMBTracks = mult->GetNumberOfTracklets();
  } 
  else if( GetAnalysisMode() == AlidNdPtHelper::kTPCITS || 
           GetAnalysisMode()==AlidNdPtHelper::kTPCITSHybrid || 
	   GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtx || 
	   GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtxUpdate || 
           GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || 
	   GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt || 
	   GetAnalysisMode() == AlidNdPtHelper::kITSStandAloneTrackSPDvtx || 
	   GetAnalysisMode() == AlidNdPtHelper::kITSStandAloneTPCTrackSPDvtx
	   )
  {
     if(vtxESD->GetStatus() && isRecVertex)
       multMBTracks = vtxESD->GetNContributors();
  }
  else {
    AliDebug(AliLog::kError, Form("Found analysis type %d", GetAnalysisMode()));
    return; 
  }

  Bool_t isEventSelected = kTRUE;
  if(evtCuts->IsEventSelectedRequired()) 
  { 
     // select events with at least 
     // one prompt track in acceptance
     // pT>0.5 GeV/c, |eta|<0.8 for the Cross Section studies

     isEventSelected = AlidNdPtHelper::SelectEvent(esdEvent,esdTrackCuts);
     //printf("isEventSelected %d \n", isEventSelected);
  }

  TObjArray *allChargedTracks=0;
  Int_t multAll=0, multAcc=0, multRec=0;
  Int_t *labelsAll=0, *labelsAcc=0, *labelsRec=0;

  // cosmics analysis
  Int_t cosmicCount = 0;
  // high-pt tracks
  Int_t highPtCount = 0;

  // check event cuts
  if(isEventOK && isEventTriggered && isEventSelected)
  {
    // get all charged tracks
    allChargedTracks = AlidNdPtHelper::GetAllChargedTracks(esdEvent,GetAnalysisMode());
    if(!allChargedTracks) return;

    Int_t entries = allChargedTracks->GetEntries();
    //printf("entries %d \n",entries);
    

    // calculate mult of reconstructed tracks
    Int_t multRecTemp=0;
    for(Int_t i=0; i<entries;++i) 
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(i);
      if(!track) continue;
      if(track->Charge()==0) continue;


      // only postive charged 
      if(GetParticleMode() == AlidNdPtHelper::kPlus && track->Charge() < 0) 
        continue;
      
      // only negative charged 
      if(GetParticleMode() == AlidNdPtHelper::kMinus && track->Charge() > 0) 
        continue;

      if(esdTrackCuts->AcceptTrack(track)) 
      {
          if(accCuts->AcceptTrack(track)) multRecTemp++;
      }  
    }

    

    labelsAll = new Int_t[entries];
    labelsAcc = new Int_t[entries];
    labelsRec = new Int_t[entries];
    for(Int_t i=0; i<entries;++i) 
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(i);
      if(!track) continue;
      if(track->Charge()==0) continue;

      // only postive charged 
      if(GetParticleMode() == AlidNdPtHelper::kPlus && track->Charge() < 0) 
        continue;
      
      // only negative charged 
      if(GetParticleMode() == AlidNdPtHelper::kMinus && track->Charge() > 0) 
        continue;


      // esd track selection 
      // ITS stand alone
      if(GetAnalysisMode() == AlidNdPtHelper::kITSStandAloneTrackSPDvtx) 
      {
        if(!(track->GetStatus() & AliESDtrack::kITSpureSA)) continue;
        if(!(track->GetStatus() & AliESDtrack::kITSrefit)) continue;
        if(track->GetNcls(0)<4) continue;
        if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1)) continue;
      } 
      else if(GetAnalysisMode() == AlidNdPtHelper::kITSStandAloneTPCTrackSPDvtx) 
      {
        //
        // ITS and TPC stand alone tracks
	//
	if(!(track->GetStatus() & AliESDtrack::kITSpureSA)) continue;
	if(!(track->GetStatus() & AliESDtrack::kITSrefit)) continue;
	if(track->GetNcls(0)<4) continue;
	if(!track->HasPointOnITSLayer(0) && !track->HasPointOnITSLayer(1)) continue;

        // Check matching with TPC only track
	Bool_t hasMatch = kFALSE;
	for(Int_t j=0; j<entries;++j) 
        {
	  if(i==j) continue;
          AliESDtrack *track2 = (AliESDtrack*)allChargedTracks->At(j);
          if(!track2) continue;
          if(track2->Charge()==0) continue;
	  if (!track2->GetTPCInnerParam()) continue;

          // check loose cuts for TPC tracks
          if(!esdTrackCuts->AcceptTrack(track2)) continue; 

	  AliExternalTrackParam *innerTPC = new AliExternalTrackParam(*(track2->GetTPCInnerParam()));
	  if(!innerTPC) continue;
          Double_t x2[3]; track2->GetXYZ(x2);
          Double_t b2[3]; AliTracker::GetBxByBz(x2,b2);
	  Double_t dz[2],cov[3];
	  Bool_t isPropOK = innerTPC->PropagateToDCABxByBz(vtxESD,b2,kVeryBig,dz,cov);
          if(!isPropOK && innerTPC) {delete innerTPC; continue;}

          // check matching
          if (TMath::Abs(track->GetY() - innerTPC->GetY()) > 3) { delete innerTPC; continue; }
          if (TMath::Abs(track->GetSnp() - innerTPC->GetSnp()) > 0.2) { delete innerTPC; continue; }
          if (TMath::Abs(track->GetTgl() - innerTPC->GetTgl()) > 0.2) { delete innerTPC; continue; }

	  hasMatch = kTRUE;
	  if(innerTPC) delete innerTPC;
	}

        if(!hasMatch) continue;
      }
      else {
        // check track cuts
        if(!esdTrackCuts->AcceptTrack(track)) 
	  continue;
      }

      //
      Bool_t isOK = kFALSE;
      Double_t x[3]; track->GetXYZ(x);
      Double_t b[3]; AliTracker::GetBxByBz(x,b);

      //
      // if TPC-ITS hybrid tracking (kTPCITSHybrid)
      // replace track parameters with TPC-ony track parameters
      //
      if( GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybrid || 
          GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || 
	  GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt) 
      {
        // Relate TPC-only tracks to Track or SPD vertex
        isOK = track->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig);
        if(!isOK) continue;

	// replace esd track parameters with TPCinner
        AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(track->GetTPCInnerParam()));
	if (!tpcTrack) return;
        track->Set(tpcTrack->GetX(),tpcTrack->GetAlpha(),tpcTrack->GetParameter(),tpcTrack->GetCovariance());

        if(tpcTrack) delete tpcTrack; 
      } 

      FillHistograms(track,stack,vtxESD->GetZv(),AlidNdPtHelper::kAllTracks, multRecTemp); 
      labelsAll[multAll] = TMath::Abs(track->GetLabel());
      multAll++;

         //FillHistograms(track,stack,AlidNdPtHelper::kAccTracks); 
	 //labelsAcc[multAcc] = TMath::Abs(track->GetLabel());
	 //multAcc++;
      
         // check high-pt tracks
         if(accCuts->AcceptTrack(track) && track->Pt() > 6.)
	 {
	     //printf(" high pt: pt %f \n",track->Pt());
	     highPtCount++;
	 }

         // check cosmics tracks
         if( GetParticleMode()==AlidNdPtHelper::kCosmic )
         {
           if(accCuts->AcceptTrack(track)) 
	   { 
             for(Int_t j=0; j<entries;++j) 
             {
               AliESDtrack *track1 = (AliESDtrack*)allChargedTracks->At(j);
               if(!track1) continue;
               if(track1->Charge()==0) continue;

               if( esdTrackCuts->AcceptTrack(track1) && accCuts->AcceptTrack(track1) ) 
               { 
		 if ( AlidNdPtHelper::IsCosmicTrack(track,track1) ) { 
		   cosmicCount++;
		   break;
		 }
	       }
	     }
	   }
          // if(!isCosmic) continue;
         }

         // update track parameters using vertex point 
	 if( GetAnalysisMode() == AlidNdPtHelper::kTPCSPDvtxUpdate || 
	     GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtxUpdate ) 
	 {
	   // update track parameters
             AliExternalTrackParam cParam;
             isOK = track->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig, &cParam);
             if(!isOK) continue;
	     track->Set(cParam.GetX(),cParam.GetAlpha(),cParam.GetParameter(),cParam.GetCovariance());

             if(accCuts->AcceptTrack(track)) {
               FillHistograms(track,stack,vtxESD->GetZv(),AlidNdPtHelper::kRecTracks,multRecTemp); 
	       labelsRec[multRec] = TMath::Abs(track->GetLabel());
	       multRec++;
	     }  
	  }
          else {
             if(accCuts->AcceptTrack(track)) 
	     {
               FillHistograms(track,stack,vtxESD->GetZv(),AlidNdPtHelper::kRecTracks,multRecTemp); 
	       labelsRec[multRec] = TMath::Abs(track->GetLabel());
	       multRec++;
	     }
          }
     }
    // if(cosmicCount) 
       //printf("COSMIC EVENT: number %d , mult %d \n", esdEvent->GetEventNumberInFile(), multRec);

    // if(highPtCount) 
       //printf("HIGH PT EVENT: number %d , mult %d \n", esdEvent->GetEventNumberInFile(), multRec);

    // if(multRec > 30) 
       //printf("HIGH MULT EVENT: number %d , mult %d \n", esdEvent->GetEventNumberInFile(), multRec);

     // fill track multiplicity histograms
     FillHistograms(allChargedTracks,labelsAll,multAll,labelsAcc,multAcc,labelsRec,multRec);

     Double_t vRecEventHist1[3] = {vtxESD->GetXv(),vtxESD->GetYv(),vtxESD->GetZv()};
     fRecEventHist1->Fill(vRecEventHist1);

     Double_t vRecEventHist2[3] = {vtxESD->GetZv(),multMBTracks,multRec};
     fRecEventHist2->Fill(vRecEventHist2);

     // 
     Double_t vRecEventHist[2] = {vtxESD->GetZv(),multMBTracks};
     fRecEventHist->Fill(vRecEventHist);
   } 

   if(IsUseMCInfo())  
   {
     if(!mcEvent) return; 

     if(evtCuts->IsEventSelectedRequired()) 
     { 
       // select events with at least 
       // one MC primary track in acceptance
       // pT>0.5 GeV/c, |eta|<0.8 for the Cross Section studies

       Bool_t isMCEventSelected = AlidNdPtHelper::SelectMCEvent(mcEvent);
       //printf("isMCEventSelected %d \n", isMCEventSelected);
       if(!isMCEventSelected) return;  
     }

     Double_t vMultTrueEventMatrix[3] = { multRec, multMCTrueTracks, multMBTracks};
     if(isEventOK && isEventTriggered) {   
       if(TMath::Abs(vtxMC[2]) < 10.0) // both Rec. and corresponding MC events must be accepted
         fEventMultCorrelationMatrix->Fill(vMultTrueEventMatrix);
     }

     // 
     // event level corrections (zv,N_MB)
     //

     // all inelastic
     Double_t vEventMatrix[2] = {vtxMC[2],multMBTracks};
     fGenEventMatrix->Fill(vEventMatrix); 
     if(isEventTriggered) fTriggerEventMatrix->Fill(vEventMatrix);
     if(isEventOK && isEventTriggered) fRecEventMatrix->Fill(vEventMatrix);

     // single diffractive
     if(evtType == AliPWG0Helper::kSD) {
       fGenSDEventMatrix->Fill(vEventMatrix); 
       if(isEventTriggered) fTriggerSDEventMatrix->Fill(vEventMatrix);
       if(isEventOK && isEventTriggered) fRecSDEventMatrix->Fill(vEventMatrix);
     }

     // double diffractive
     if(evtType == AliPWG0Helper::kDD) {
       fGenDDEventMatrix->Fill(vEventMatrix); 
       if(isEventTriggered) fTriggerDDEventMatrix->Fill(vEventMatrix);
       if(isEventOK && isEventTriggered)  fRecDDEventMatrix->Fill(vEventMatrix);
     }

     // non diffractive
     if(evtType == AliPWG0Helper::kND) {
       fGenNDEventMatrix->Fill(vEventMatrix); 
       if(isEventTriggered) fTriggerNDEventMatrix->Fill(vEventMatrix);
       if(isEventOK && isEventTriggered) fRecNDEventMatrix->Fill(vEventMatrix);
     }

     // non single diffractive
     if(evtType != AliPWG0Helper::kSD) {
       fGenNSDEventMatrix->Fill(vEventMatrix); 
       if(isEventTriggered) fTriggerNSDEventMatrix->Fill(vEventMatrix);
       if(isEventOK && isEventTriggered) fRecNSDEventMatrix->Fill(vEventMatrix);
     }

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

       // checked accepted
       if(accCuts->AcceptTrack(particle)) 
       {
         Double_t vTrackEventMatrix[3] = {vtxMC[2], particle->Pt(), particle->Eta()}; 
         fGenTrackEventMatrix->Fill(vTrackEventMatrix);

         if(evtType == AliPWG0Helper::kSD) {
           fGenTrackSDEventMatrix->Fill(vTrackEventMatrix);
	 }
         if(evtType == AliPWG0Helper::kDD) {
           fGenTrackDDEventMatrix->Fill(vTrackEventMatrix);
	 }
         if(evtType == AliPWG0Helper::kND) {
           fGenTrackNDEventMatrix->Fill(vTrackEventMatrix);
	 }
         if(evtType != AliPWG0Helper::kSD) {
           fGenTrackNSDEventMatrix->Fill(vTrackEventMatrix);
	 }

         //
         if(!isEventTriggered) continue;  

         fTriggerTrackEventMatrix->Fill(vTrackEventMatrix);
         if(evtType == AliPWG0Helper::kSD) {
           fTriggerTrackSDEventMatrix->Fill(vTrackEventMatrix);
	 }
         if(evtType == AliPWG0Helper::kDD) {
           fTriggerTrackDDEventMatrix->Fill(vTrackEventMatrix);
	 }
         if(evtType == AliPWG0Helper::kND) {
           fTriggerTrackNDEventMatrix->Fill(vTrackEventMatrix);
	 }
         if(evtType != AliPWG0Helper::kSD) {
           fTriggerTrackNSDEventMatrix->Fill(vTrackEventMatrix);
	 }

         //
    	 if(!isEventOK) continue;  

         fRecTrackEventMatrix->Fill(vTrackEventMatrix);
         if(evtType == AliPWG0Helper::kSD) {
           fRecTrackSDEventMatrix->Fill(vTrackEventMatrix);
	 }
         if(evtType == AliPWG0Helper::kDD) {
           fRecTrackDDEventMatrix->Fill(vTrackEventMatrix);
	 }
         if(evtType == AliPWG0Helper::kND) {
           fRecTrackNDEventMatrix->Fill(vTrackEventMatrix);
	 }
         if(evtType != AliPWG0Helper::kSD) {
           fRecTrackNSDEventMatrix->Fill(vTrackEventMatrix);
	 }
       }
     }

     // 
     // track-level corrections (zv,pt,eta)
     //
     if(isEventOK && isEventTriggered)
     {

       // fill MC and rec event control histograms
       if(fHistogramsOn) {
         Double_t vRecMCEventHist1[3] = {vtxESD->GetXv()-vtxMC[0],vtxESD->GetYv()-vtxMC[1],vtxESD->GetZv()-vtxMC[2]};
         fRecMCEventHist1->Fill(vRecMCEventHist1);

         Double_t vRecMCEventHist2[3] = {vtxESD->GetXv()-vtxMC[0],vtxESD->GetZv()-vtxMC[2],multMBTracks};
         fRecMCEventHist2->Fill(vRecMCEventHist2);

         Double_t vRecMCEventHist3[2] = {multRec,evtType};
         fRecMCEventHist3->Fill(vRecMCEventHist3);
       }

       //
       // MC histograms for track efficiency studies
       //
       Int_t countRecCandle = 0;
       for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc) 
       {
         TParticle* particle = stack->Particle(iMc);
         if (!particle)
         continue;

         Double_t vTrackMatrix[3] = {vtxMC[2],particle->Pt(),particle->Eta()}; 

	 // all genertated primaries including neutral
         if( iMc < stack->GetNprimary() ) {
           //fGenTrackMatrix->Fill(vTrackMatrix);
	 }

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

         // check accepted
         if(accCuts->AcceptTrack(particle)) 
	 {

           if( AlidNdPtHelper::IsPrimaryParticle(stack, iMc, GetParticleMode()) ) 
	     fGenPrimTrackMatrix->Fill(vTrackMatrix);

	   // fill control histograms
           if(fHistogramsOn) 
	     FillHistograms(stack,iMc,AlidNdPtHelper::kAccTracks); 

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
                   Double_t vMCMultRecTrackHist1[3] = {particle->Pt(), particle->Eta(), pid};
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
	         //AliESDtrack *track = esdEvent->GetTrack(iRec);
                 //if(track && track->GetKinkIndex(0) < 0) 
		 //  printf("prim kink \n");
		 countRecCandle++;
	       }

               if(!prim) fRecSecTrackMatrix->Fill(vTrackMatrix);

	       // fill control histograms
               if(fHistogramsOn) 
                 FillHistograms(stack,iMc,AlidNdPtHelper::kRecTracks); 
	      
               break;
             }
           }
         }
       }

       if(countRecCandle>0) fRecCandleEventMatrix->Fill(vEventMatrix);
     }


   } // end bUseMC

  if(allChargedTracks) delete allChargedTracks; allChargedTracks = 0;
  if(labelsAll) delete [] labelsAll; labelsAll = 0;
  if(labelsAcc) delete [] labelsAcc; labelsAcc = 0;
  if(labelsRec) delete [] labelsRec; labelsRec = 0;

  if(!evtCuts->IsRecVertexRequired() && vtxESD != NULL) delete vtxESD;
  //if(trigAna) delete trigAna;

}

//_____________________________________________________________________________
void AlidNdPtAnalysis::FillHistograms(TObjArray *const allChargedTracks,Int_t *const labelsAll,Int_t multAll,Int_t *const labelsAcc,Int_t multAcc,Int_t *const labelsRec,Int_t multRec) {
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
         Double_t v1[2] = {track->Pt(), multAll}; 
         fRecTrackMultHist1[AlidNdPtHelper::kAllTracks]->Fill(v1);
       }
     }
     for(Int_t iAcc=0; iAcc<multAcc; ++iAcc) {
       if(label == labelsAcc[iAcc]) {
         Double_t v2[2] = {track->Pt(), multAcc}; 
         fRecTrackMultHist1[AlidNdPtHelper::kAccTracks]->Fill(v2);
       }
     }
     for(Int_t iRec=0; iRec<multRec; ++iRec) {
       if(label == labelsRec[iRec]) {
         Double_t v3[2] = {track->Pt(), multRec}; 
         fRecTrackMultHist1[AlidNdPtHelper::kRecTracks]->Fill(v3);
       }
     }
  }
}

//_____________________________________________________________________________
void AlidNdPtAnalysis::FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, const Double_t zv, AlidNdPtHelper::TrackObject trackObj, Int_t multMB)
{
  //
  // Fill ESD track and MC histograms 
  //
  if(!esdTrack) return;

  Float_t q = esdTrack->Charge();
  if(TMath::Abs(q) < 0.001) return;

  Float_t pt = esdTrack->Pt();
  //Float_t qpt = esdTrack->Pt() * q;
  Float_t eta = esdTrack->Eta();
  Float_t phi = esdTrack->Phi();

  Float_t dca[2], bCov[3];
  esdTrack->GetImpactParameters(dca,bCov);

  Int_t nClust = esdTrack->GetTPCclusters(0);
  Float_t chi2PerCluster = 0.;
  if(nClust>0.) chi2PerCluster = esdTrack->GetTPCchi2()/Float_t(nClust);


  // fill histograms
  Double_t values1[3] = {pt,eta,phi};	  
  fRecTrackHist1[trackObj]->Fill(values1);

  Double_t values[4] = {zv, pt,eta, multMB};	  
  if(trackObj == AlidNdPtHelper::kRecTracks) {
    fRecTrackHist->Fill(values);
  }

  /*
  Double_t values2[5] = {nClust,chi2PerCluster,pt,eta,phi};	  
  if(trackObj == AlidNdPtHelper::kRecTracks)  
  {
    if(fHistogramsOn)
      fRecTrackHist2->Fill(values2);
  }
  */
 
  //
  // Fill rec vs MC information
  //
  if(!stack) return;

  Int_t label = TMath::Abs(esdTrack->GetLabel()); 
  //if(label == 0) return;

  TParticle* particle = stack->Particle(label);
  if(!particle) return;

  //Bool_t prim = stack->IsPhysicalPrimary(label);
  //Int_t pid = AlidNdPtHelper::ConvertPdgToPid(particle);

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
  //Float_t qgpt = particle->Pt() * gq;
  //Float_t gphi = particle->Phi();

  Double_t dpt=0;
  //printf("pt %f, gpt %f \n",pt,gpt);
  if(gpt) dpt = (pt-gpt)/gpt;
  Double_t deta = (eta-geta);
 
  // fill histograms
  if(trackObj == AlidNdPtHelper::kRecTracks)  
  {
    Double_t vTrackPtCorrelationMatrix[3]={pt,gpt,geta};
    fTrackPtCorrelationMatrix->Fill(vTrackPtCorrelationMatrix);

    Double_t vRecMCTrackHist1[4]={gpt,geta,dpt,deta};
    fRecMCTrackHist1->Fill(vRecMCTrackHist1);
  }
}

//_____________________________________________________________________________
void AlidNdPtAnalysis::FillHistograms(AliStack *const stack, Int_t label, AlidNdPtHelper::TrackObject trackObj)
{
  // Fill MC histograms
  if(!stack) return;

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
  Double_t vMCTrackHist1[3] = {gpt,geta,gphi};
  fMCTrackHist1[trackObj]->Fill(vMCTrackHist1);

  Double_t vMCPrimTrackHist1[5] = {gpt,geta,pid,mech,motherPdg};
  Double_t vMCPrimTrackHist2[5] = {TMath::Abs(particle->GetPdgCode()),mech,motherPdg};
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
Long64_t AlidNdPtAnalysis::Merge(TCollection* const list) 
{
  //  init if not done already
  if (!fIsInit) { Init(); }
  
  // Merge list of objects (needed by PROOF)

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
    AlidNdPtAnalysis* entry = dynamic_cast<AlidNdPtAnalysis*>(obj);
    if (entry == 0) continue; 

    // physics selection
    //printf("entry->GetPhysicsTriggerSelection() %p \n", entry->GetPhysicsTriggerSelection());
    //collPhysSelection->Add(entry->GetPhysicsTriggerSelection());
    
    //
    fRecEventHist->Add(entry->fRecEventHist);
    fRecTrackHist->Add(entry->fRecTrackHist);

    //
    fEventMultCorrelationMatrix->Add(entry->fEventMultCorrelationMatrix);
    fTrackPtCorrelationMatrix->Add(entry->fTrackPtCorrelationMatrix);

    //
    fGenEventMatrix->Add(entry->fGenEventMatrix);
    fGenSDEventMatrix->Add(entry->fGenSDEventMatrix);
    fGenDDEventMatrix->Add(entry->fGenDDEventMatrix);
    fGenNDEventMatrix->Add(entry->fGenNDEventMatrix);
    fGenNSDEventMatrix->Add(entry->fGenNSDEventMatrix);

    fTriggerEventMatrix->Add(entry->fTriggerEventMatrix);
    fTriggerSDEventMatrix->Add(entry->fTriggerSDEventMatrix);
    fTriggerDDEventMatrix->Add(entry->fTriggerDDEventMatrix);
    fTriggerNDEventMatrix->Add(entry->fTriggerNDEventMatrix);
    fTriggerNSDEventMatrix->Add(entry->fTriggerNSDEventMatrix);

    fRecEventMatrix->Add(entry->fRecEventMatrix);
    fRecSDEventMatrix->Add(entry->fRecSDEventMatrix);
    fRecDDEventMatrix->Add(entry->fRecDDEventMatrix);
    fRecNDEventMatrix->Add(entry->fRecNDEventMatrix);
    fRecNSDEventMatrix->Add(entry->fRecNSDEventMatrix);

    fRecCandleEventMatrix->Add(entry->fRecCandleEventMatrix);
    //
    fGenTrackEventMatrix->Add(entry->fGenTrackEventMatrix);
    fGenTrackSDEventMatrix->Add(entry->fGenTrackSDEventMatrix);
    fGenTrackDDEventMatrix->Add(entry->fGenTrackDDEventMatrix);
    fGenTrackNDEventMatrix->Add(entry->fGenTrackNDEventMatrix);
    fGenTrackNSDEventMatrix->Add(entry->fGenTrackNSDEventMatrix);

    fTriggerTrackEventMatrix->Add(entry->fTriggerTrackEventMatrix);
    fTriggerTrackSDEventMatrix->Add(entry->fTriggerTrackSDEventMatrix);
    fTriggerTrackDDEventMatrix->Add(entry->fTriggerTrackDDEventMatrix);
    fTriggerTrackNDEventMatrix->Add(entry->fTriggerTrackNDEventMatrix);
    fTriggerTrackNSDEventMatrix->Add(entry->fTriggerTrackNSDEventMatrix);

    fRecTrackEventMatrix->Add(entry->fRecTrackEventMatrix);
    fRecTrackSDEventMatrix->Add(entry->fRecTrackSDEventMatrix);
    fRecTrackDDEventMatrix->Add(entry->fRecTrackDDEventMatrix);
    fRecTrackNDEventMatrix->Add(entry->fRecTrackNDEventMatrix);
    fRecTrackNSDEventMatrix->Add(entry->fRecTrackNSDEventMatrix);

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
    fRecMCEventHist3->Add(entry->fRecMCEventHist3);

    for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) {
      fMCTrackHist1[i]->Add(entry->fMCTrackHist1[i]);

      fMCPrimTrackHist1[i]->Add(entry->fMCPrimTrackHist1[i]);
      fMCPrimTrackHist2[i]->Add(entry->fMCPrimTrackHist2[i]);
      fMCSecTrackHist1[i]->Add(entry->fMCSecTrackHist1[i]);

      fRecTrackHist1[i]->Add(entry->fRecTrackHist1[i]);
      fRecTrackMultHist1[i]->Add(entry->fRecTrackMultHist1[i]);
    }
    fRecMCTrackHist1->Add(entry->fRecMCTrackHist1);
    fMCMultRecTrackHist1->Add(entry->fMCMultRecTrackHist1);
    fRecTrackHist2->Add(entry->fRecTrackHist2);

  count++;
  }

  //AliPhysicsSelection *trigSelection = GetPhysicsTriggerSelection();
  //trigSelection->Merge(collPhysSelection);
  //if(collPhysSelection) delete collPhysSelection;

return count;
}

//____________________________________________________________________
Bool_t AlidNdPtAnalysis::IsBinZeroTrackSPDvtx(const AliESDEvent* esdEvent) {
//
// check 0-bin
// for LHC background calculation
// return kTRUE if vertex not reconstructed or
// track multiplicity == 0
//
if(!esdEvent) return kFALSE;

  // check vertex
  const AliESDVertex *vertex = esdEvent->GetPrimaryVertexTracks();
  if(!vertex) return kTRUE;

  if(vertex->GetNContributors() < 1) {
    // SPD vertex
    vertex = esdEvent->GetPrimaryVertexSPD();
    if(!vertex) return kTRUE;
  }
  if(TMath::Abs(vertex->GetZv()) > 15.0) return kTRUE; 
  if( !AlidNdPtHelper::TestRecVertex(vertex, esdEvent->GetPrimaryVertexSPD(), AlidNdPtHelper::kTPCITSHybridTrackSPDvtx) ) 
    return kTRUE;
 
return kFALSE;
}

//____________________________________________________________________
Bool_t AlidNdPtAnalysis::IsBinZeroSPDvtx(const AliESDEvent* esdEvent) {
//
// check 0-bin
// for LHC background calculation
// return kTRUE if vertex not reconstructed or
// tracklet multiplicity == 0
//
if(!esdEvent) return kFALSE;

  // check vertex
  const AliESDVertex* vertex = esdEvent->GetPrimaryVertexSPD();
  if(!vertex) return kTRUE;
  if(TMath::Abs(vertex->GetZv()) > 15.0) return kTRUE; 
  if( !AlidNdPtHelper::TestRecVertex(vertex, esdEvent->GetPrimaryVertexSPD(), AlidNdPtHelper::kTPCSPDvtx) ) return kTRUE;
 
return kFALSE;
}

//_____________________________________________________________________________
void AlidNdPtAnalysis::Analyse() 
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

  h = fRecEventHist2->Projection(2);
  h->SetName("multiplicity");
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

    h2D = fRecTrackHist2->Projection(0,1);
    h2D->SetName("nClust_chi2_rec");
    aFolderObj->Add(h2D);

    h2D = fRecTrackHist2->Projection(0,2);
    h2D->SetName("nClust_pt_rec");
    aFolderObj->Add(h2D);

    h2D = fRecTrackHist2->Projection(0,3);
    h2D->SetName("nClust_eta_rec");
    aFolderObj->Add(h2D);

    h2D = fRecTrackHist2->Projection(0,4);
    h2D->SetName("nClust_phi_rec");
    aFolderObj->Add(h2D);

    h2D = fRecTrackHist2->Projection(1,2);
    h2D->SetName("chi2_pt_rec");
    aFolderObj->Add(h2D);

    h2D = fRecTrackHist2->Projection(1,3);
    h2D->SetName("chi2_eta_rec");
    aFolderObj->Add(h2D);

    h2D = fRecTrackHist2->Projection(1,4);
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
  // rec. vs true multiplicity correlation matrix
  //
  hs = (THnSparse*)fEventMultCorrelationMatrix->Clone("event_mult_correlation_matrix");
  aFolderObj->Add(hs);
 
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
  // trigger efficiency for NSD
  //
  h = AlidNdPtHelper::GenerateCorrMatrix(fTriggerNSDEventMatrix->Projection(0),fGenNSDEventMatrix->Projection(0),"zv_trig_NSD_eff_matrix");
  aFolderObj->Add(h);

  //
  // trigger bias correction (MB to ND)
  //
  hs = AlidNdPtHelper::GenerateCorrMatrix(fGenNDEventMatrix,fTriggerEventMatrix,"zv_mult_trig_MBtoND_corr_matrix");
  aFolderObj->Add(hs);

  h = AlidNdPtHelper::GenerateCorrMatrix(fGenNDEventMatrix->Projection(0),fTriggerEventMatrix->Projection(0),"zv_trig_MBtoND_corr_matrix");
  aFolderObj->Add(h);


  h = AlidNdPtHelper::GenerateCorrMatrix(fGenNDEventMatrix->Projection(1),fTriggerEventMatrix->Projection(1),"mult_trig_MBtoND_corr_matrix");

  aFolderObj->Add(h);

  //
  // trigger bias correction (MB to NSD)
  //
  hs = AlidNdPtHelper::GenerateCorrMatrix(fGenNSDEventMatrix,fTriggerEventMatrix,"zv_mult_trig_MBtoNSD_corr_matrix");
  aFolderObj->Add(hs);

  h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenNSDEventMatrix->Projection(0,1),fTriggerEventMatrix->Projection(0,1),"zv_mult_trig_MBtoNSD_corr_matrix_2D");
  aFolderObj->Add(h2D);

  h = AlidNdPtHelper::GenerateCorrMatrix(fGenNSDEventMatrix->Projection(0),fTriggerEventMatrix->Projection(0),"zv_trig_MBtoNSD_corr_matrix");
  aFolderObj->Add(h);


  h = AlidNdPtHelper::GenerateCorrMatrix(fGenNSDEventMatrix->Projection(1),fTriggerEventMatrix->Projection(1),"mult_trig_MBtoNSD_corr_matrix");
  aFolderObj->Add(h);


  //
  // trigger bias correction (MB to INEL)
  //
  hs = AlidNdPtHelper::GenerateCorrMatrix(fGenEventMatrix,fTriggerEventMatrix,"zv_mult_trig_MBtoInel_corr_matrix");
  aFolderObj->Add(hs);

  h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenEventMatrix->Projection(0,1),fTriggerEventMatrix->Projection(0,1),"zv_mult_trig_MBtoInel_corr_matrix_2D");
  aFolderObj->Add(h2D);

  h = AlidNdPtHelper::GenerateCorrMatrix(fGenEventMatrix->Projection(0),fTriggerEventMatrix->Projection(0),"zv_trig_MBtoInel_corr_matrix");
  aFolderObj->Add(h);

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
  // track-event trigger bias correction (MB to ND)
  //

  hs = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackNDEventMatrix,fTriggerTrackEventMatrix,"zv_pt_eta_track_trig_MBtoND_corr_matrix");
  aFolderObj->Add(hs);

  h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackNDEventMatrix->Projection(1,2),fTriggerTrackEventMatrix->Projection(1,2),"eta_pt_track_trig_MBtoND_corr_matrix");
  aFolderObj->Add(h2D);

  h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackNDEventMatrix->Projection(1,0),fTriggerTrackEventMatrix->Projection(1,0),"pt_zv_track_trig_MBtoND_corr_matrix");
  aFolderObj->Add(h2D);

  h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackNDEventMatrix->Projection(2,0),fTriggerTrackEventMatrix->Projection(2,0),"zv_eta_track_trig_MBtoND_corr_matrix");
  aFolderObj->Add(h2D);

  //
  // track-event trigger bias correction (MB to NSD)
  //
  hs = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackNSDEventMatrix,fTriggerTrackEventMatrix,"zv_pt_eta_track_trig_MBtoNSD_corr_matrix");
  aFolderObj->Add(hs);

  h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackNSDEventMatrix->Projection(1,2),fTriggerTrackEventMatrix->Projection(1,2),"eta_pt_track_trig_MBtoNSD_corr_matrix");
  aFolderObj->Add(h2D);

  h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackNSDEventMatrix->Projection(1,0),fTriggerTrackEventMatrix->Projection(1,0),"pt_zv_track_trig_MBtoNSD_corr_matrix");
  aFolderObj->Add(h2D);

  h2D = AlidNdPtHelper::GenerateCorrMatrix(fGenTrackNSDEventMatrix->Projection(2,0),fTriggerTrackEventMatrix->Projection(2,0),"zv_eta_track_trig_MBtoNSD_corr_matrix");
  aFolderObj->Add(h2D);


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
    sprintf(name,"mc_pt_rec_prim_pid_%d",i); 
    fMCPrimTrackHist1[2]->GetAxis(2)->SetRange(i+1,i+1);
    h = fMCPrimTrackHist1[2]->Projection(0);
    h->SetName(name);
    aFolderObj->Add(h);

    sprintf(name,"mc_pt_rec_sec_pid_%d",i); 
    fMCSecTrackHist1[2]->GetAxis(2)->SetRange(i+1,i+1);
    h = fMCSecTrackHist1[2]->Projection(0);
    h->SetName(name);
    aFolderObj->Add(h);

    // production mechanisms for given pid
    fMCSecTrackHist1[2]->GetAxis(2)->SetRange(i+1,i+1);

    for(Int_t j=0; j<20; j++) {
      if(j == 4) {
        // decay
	
        sprintf(name,"mc_pt_rec_sec_pid_%d_decay",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(0);
        h->SetName(name);
        aFolderObj->Add(h);

        sprintf(name,"mc_eta_rec_sec_pid_%d_decay",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(1);
        h->SetName(name);
        aFolderObj->Add(h);

        sprintf(name,"mc_mother_rec_sec_pid_%d_decay",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(4);
        h->SetName(name);
        aFolderObj->Add(h);

      } else if (j == 5) {
       // conversion

        sprintf(name,"mc_pt_rec_sec_pid_%d_conv",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(0);
        h->SetName(name);
        aFolderObj->Add(h);

        sprintf(name,"mc_eta_rec_sec_pid_%d_conv",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(1);
        h->SetName(name);
        aFolderObj->Add(h);

        sprintf(name,"mc_mother_rec_sec_pid_%d_conv",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(4);
        h->SetName(name);
        aFolderObj->Add(h);

     } else if (j == 13) {
       // mat
       
        sprintf(name,"mc_pt_rec_sec_pid_%d_mat",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(0);
        h->SetName(name);
        aFolderObj->Add(h);

        sprintf(name,"mc_eta_rec_sec_pid_%d_mat",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(1);
        h->SetName(name);
        aFolderObj->Add(h);

        sprintf(name,"mc_eta_mother_rec_sec_pid_%d_mat",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(4,1);
        h->SetName(name);
        aFolderObj->Add(h);

        sprintf(name,"mc_mother_rec_sec_pid_%d_mat",i); 
        fMCSecTrackHist1[2]->GetAxis(3)->SetRange(j+1,j+1);
        h = fMCSecTrackHist1[2]->Projection(4);
        h->SetName(name);
        aFolderObj->Add(h);

        sprintf(name,"mc_pt_mother_rec_sec_pid_%d_mat",i); 
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
TFolder* AlidNdPtAnalysis::ExportToFolder(TObjArray * const array) 
{
  // recreate folder avery time and export objects to new one
  //
  AlidNdPtAnalysis * comp=this;
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
TFolder* AlidNdPtAnalysis::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
