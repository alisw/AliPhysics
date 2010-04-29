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
// AlidNdPtCorrection class:
//
// a. functionality:
// - applies corrections on dNdPt spectra
// - fills corrected dNdPt histograms
// - fills correction control histograms 
//
// b. data members:
// - dNdPt spectra before and after correction procedure
// - control histograms
// - correction matrices (must be loaded)
// 
// Author: J.Otwinowski 04/11/2008 
//------------------------------------------------------------------------------

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "AliHeader.h"  
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

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"
#include "AlidNdPtAnalysis.h"
#include "AlidNdPtCorrection.h"

using namespace std;

ClassImp(AlidNdPtCorrection)

//_____________________________________________________________________________
//AlidNdPtCorrection::AlidNdPtCorrection(): TNamed(),
  AlidNdPtCorrection::AlidNdPtCorrection(): AlidNdPt(),
  fCorrectionFolder(0),
  fMCEventHist1(0),
  fRecEventHist1(0),
  fRecEventMultHist1(0),
  fMCAllEventMultHist1(0),
  fMCAllNDEventMultHist1(0),
  fMCAllNSDEventMultHist1(0),
  fMCTriggerMultHist1(0),
  fMCEventMultHist1(0),
  fMCAllPrimTrackMultHist1(0),
  fMCNDEventAllPrimTrackMultHist1(0),
  fMCNSDEventAllPrimTrackMultHist1(0),
  fMCTriggerPrimTrackMultHist1(0),
  fMCEventPrimTrackMultHist1(0),
  fMCAllPrimTrackTrueMultHist1(0),
  fMCNDEventAllPrimTrackTrueMultHist1(0),
  fMCNSDEventAllPrimTrackTrueMultHist1(0),
  fMCTriggerPrimTrackTrueMultHist1(0),
  fMCEventPrimTrackTrueMultHist1(0),
  fMCAllPrimTrackTrueMultHist2(0),
  fMCNDEventAllPrimTrackTrueMultHist2(0),
  fMCNSDEventAllPrimTrackTrueMultHist2(0),
  fMCTriggerPrimTrackTrueMultHist2(0),
  fMCEventPrimTrackTrueMultHist2(0),
  fMCAllPrimTrackMeanPtMult1(0),
  fMCNDEventAllPrimTrackMeanPtMult1(0),
  fMCNSDEventAllPrimTrackMeanPtMult1(0),
  fMCTriggerPrimTrackMeanPtMult1(0),
  fMCEventPrimTrackMeanPtMult1(0),
  fMCAllPrimTrackMeanPtTrueMult1(0),
  fMCNDEventAllPrimTrackMeanPtTrueMult1(0),
  fMCNSDEventAllPrimTrackMeanPtTrueMult1(0),
  fMCTriggerPrimTrackMeanPtTrueMult1(0),
  fMCEventPrimTrackMeanPtTrueMult1(0),
  fEventMultCorrelationMatrix(0),
  fZvNorm(0),
  fZvEmptyEventsNorm(0),
  fLHCBin0Background(0),
  fCorrTriggerMBtoInelEventMatrix(0),
  fCorrTriggerMBtoNDEventMatrix(0),
  fCorrTriggerMBtoNSDEventMatrix(0),
  fCorrEventMatrix(0),
  fCorrTriggerMBtoInelTrackEventMatrix(0),
  fCorrTriggerMBtoNDTrackEventMatrix(0),
  fCorrTriggerMBtoNSDTrackEventMatrix(0),
  fCorrTrackEventMatrix(0),
  fCorrTrackMatrix(0),
  fCorrHighPtTrackMatrix(0),
  fContTrackMatrix(0),
  fContMultTrackMatrix(0),
  fCorrMatrixFileName(""),
  fCosmicsHisto(0),
  fEventCount(0)
{
  // default constructor
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    fRecTrackHist1[i]=0;     
  }

  for(Int_t i=0; i<8; i++) { 
    fCorrRecTrackMultHist1[i] = 0;
    fCorrRecTrackTrueMultHist1[i] = 0;
    fCorrRecTrackTrueMultHist2[i] = 0;
    fCorrRecTrackMeanPtMultHist1[i] = 0;
    fCorrRecTrackMeanPtTrueMultHist1[i] = 0;
    fCorrRecTrackPt1[i] = 0;
  }

  for(Int_t i=0; i<5; i++) { 
    fCorrRecEventHist1[i] = 0;
    fCorrRecEventHist2[i] = 0;
  }

  Init();
}

//_____________________________________________________________________________
AlidNdPtCorrection::AlidNdPtCorrection(Char_t* name, Char_t* title, TString corrMatrixFileName): AlidNdPt(name,title),
  fCorrectionFolder(0),
  fMCEventHist1(0),
  fRecEventHist1(0),
  fRecEventMultHist1(0),
  fMCAllEventMultHist1(0),
  fMCAllNDEventMultHist1(0),
  fMCAllNSDEventMultHist1(0),
  fMCTriggerMultHist1(0),
  fMCEventMultHist1(0),
  fMCAllPrimTrackMultHist1(0),
  fMCNDEventAllPrimTrackMultHist1(0),
  fMCNSDEventAllPrimTrackMultHist1(0),
  fMCTriggerPrimTrackMultHist1(0),
  fMCEventPrimTrackMultHist1(0),
  fMCAllPrimTrackTrueMultHist1(0),
  fMCNDEventAllPrimTrackTrueMultHist1(0),
  fMCNSDEventAllPrimTrackTrueMultHist1(0),
  fMCTriggerPrimTrackTrueMultHist1(0),
  fMCEventPrimTrackTrueMultHist1(0),
  fMCAllPrimTrackTrueMultHist2(0),
  fMCNDEventAllPrimTrackTrueMultHist2(0),
  fMCNSDEventAllPrimTrackTrueMultHist2(0),
  fMCTriggerPrimTrackTrueMultHist2(0),
  fMCEventPrimTrackTrueMultHist2(0),
  fMCAllPrimTrackMeanPtMult1(0),
  fMCNDEventAllPrimTrackMeanPtMult1(0),
  fMCNSDEventAllPrimTrackMeanPtMult1(0),
  fMCTriggerPrimTrackMeanPtMult1(0),
  fMCEventPrimTrackMeanPtMult1(0),
  fMCAllPrimTrackMeanPtTrueMult1(0),
  fMCNDEventAllPrimTrackMeanPtTrueMult1(0),
  fMCNSDEventAllPrimTrackMeanPtTrueMult1(0),
  fMCTriggerPrimTrackMeanPtTrueMult1(0),
  fMCEventPrimTrackMeanPtTrueMult1(0),
  fEventMultCorrelationMatrix(0),
  fZvNorm(0),
  fZvEmptyEventsNorm(0),
  fLHCBin0Background(0),
  fCorrTriggerMBtoInelEventMatrix(0),
  fCorrTriggerMBtoNDEventMatrix(0),
  fCorrTriggerMBtoNSDEventMatrix(0),
  fCorrEventMatrix(0),
  fCorrTriggerMBtoInelTrackEventMatrix(0),
  fCorrTriggerMBtoNDTrackEventMatrix(0),
  fCorrTriggerMBtoNSDTrackEventMatrix(0),
  fCorrTrackEventMatrix(0),
  fCorrTrackMatrix(0),
  fCorrHighPtTrackMatrix(0),
  fContTrackMatrix(0),
  fContMultTrackMatrix(0),
  fCorrMatrixFileName(corrMatrixFileName),
  fCosmicsHisto(0),
  fEventCount(0)
{
  // constructor
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    fRecTrackHist1[i]=0;     
  }

  for(Int_t i=0; i<8; i++) { 
    fCorrRecTrackMultHist1[i] = 0;
    fCorrRecTrackTrueMultHist1[i] = 0;
    fCorrRecTrackTrueMultHist2[i] = 0;
    fCorrRecTrackMeanPtMultHist1[i] = 0;
    fCorrRecTrackMeanPtTrueMultHist1[i] = 0;
    fCorrRecTrackPt1[i] = 0;
  }

  for(Int_t i=0; i<5; i++) { 
    fCorrRecEventHist1[i] = 0;
    fCorrRecEventHist2[i] = 0;
  }

  Init();
}

//_____________________________________________________________________________
AlidNdPtCorrection::~AlidNdPtCorrection() {
  // 
  // destructor
  //
  if(fMCEventHist1) delete fMCEventHist1; fMCEventHist1=0;
  if(fRecEventHist1) delete fRecEventHist1; fRecEventHist1=0;
  if(fRecEventMultHist1) delete fRecEventMultHist1; fRecEventMultHist1=0;

  if(fMCAllEventMultHist1) delete fMCAllEventMultHist1; fMCAllEventMultHist1=0;
  if(fMCAllNDEventMultHist1) delete fMCAllNDEventMultHist1; fMCAllNDEventMultHist1=0;
  if(fMCAllNSDEventMultHist1) delete fMCAllNSDEventMultHist1; fMCAllNSDEventMultHist1=0;
  if(fMCTriggerMultHist1) delete fMCTriggerMultHist1; fMCTriggerMultHist1=0;
  if(fMCEventMultHist1) delete fMCEventMultHist1; fMCEventMultHist1=0;

  if(fMCAllPrimTrackMultHist1) delete fMCAllPrimTrackMultHist1; fMCAllPrimTrackMultHist1=0;
  if(fMCNDEventAllPrimTrackMultHist1) delete fMCNDEventAllPrimTrackMultHist1; fMCNDEventAllPrimTrackMultHist1=0;
  if(fMCNSDEventAllPrimTrackMultHist1) delete fMCNSDEventAllPrimTrackMultHist1; fMCNSDEventAllPrimTrackMultHist1=0;
  if(fMCTriggerPrimTrackMultHist1) delete fMCTriggerPrimTrackMultHist1; fMCTriggerPrimTrackMultHist1=0;
  if(fMCEventPrimTrackMultHist1) delete fMCEventPrimTrackMultHist1; fMCEventPrimTrackMultHist1=0;

  if(fMCAllPrimTrackTrueMultHist1) delete fMCAllPrimTrackTrueMultHist1; fMCAllPrimTrackTrueMultHist1=0;
  if(fMCNDEventAllPrimTrackTrueMultHist1) delete fMCNDEventAllPrimTrackTrueMultHist1; fMCNDEventAllPrimTrackTrueMultHist1=0;
  if(fMCNSDEventAllPrimTrackTrueMultHist1) delete fMCNSDEventAllPrimTrackTrueMultHist1; fMCNSDEventAllPrimTrackTrueMultHist1=0;
  if(fMCTriggerPrimTrackTrueMultHist1) delete fMCTriggerPrimTrackTrueMultHist1; fMCTriggerPrimTrackTrueMultHist1=0;
  if(fMCEventPrimTrackTrueMultHist1) delete fMCEventPrimTrackTrueMultHist1; fMCEventPrimTrackTrueMultHist1=0;

  if(fMCAllPrimTrackTrueMultHist2) delete fMCAllPrimTrackTrueMultHist2; fMCAllPrimTrackTrueMultHist2=0;
  if(fMCNDEventAllPrimTrackTrueMultHist2) delete fMCNDEventAllPrimTrackTrueMultHist2; fMCNDEventAllPrimTrackTrueMultHist2=0;
  if(fMCNSDEventAllPrimTrackTrueMultHist2) delete fMCNSDEventAllPrimTrackTrueMultHist2; fMCNSDEventAllPrimTrackTrueMultHist2=0;
  if(fMCTriggerPrimTrackTrueMultHist2) delete fMCTriggerPrimTrackTrueMultHist2; fMCTriggerPrimTrackTrueMultHist2=0;
  if(fMCEventPrimTrackTrueMultHist2) delete fMCEventPrimTrackTrueMultHist2; fMCEventPrimTrackTrueMultHist2=0;




  if(fMCAllPrimTrackMeanPtMult1) delete fMCAllPrimTrackMeanPtMult1; fMCAllPrimTrackMeanPtMult1=0;
  if(fMCNDEventAllPrimTrackMeanPtMult1) delete fMCNDEventAllPrimTrackMeanPtMult1; fMCNDEventAllPrimTrackMeanPtMult1=0;
  if(fMCNSDEventAllPrimTrackMeanPtMult1) delete fMCNSDEventAllPrimTrackMeanPtMult1; fMCNSDEventAllPrimTrackMeanPtMult1=0;
  if(fMCTriggerPrimTrackMeanPtMult1) delete fMCTriggerPrimTrackMeanPtMult1; fMCTriggerPrimTrackMeanPtMult1=0;
  if(fMCEventPrimTrackMeanPtMult1) delete fMCEventPrimTrackMeanPtMult1; fMCEventPrimTrackMeanPtMult1=0;

  if(fMCAllPrimTrackMeanPtTrueMult1) delete fMCAllPrimTrackMeanPtTrueMult1; fMCAllPrimTrackMeanPtTrueMult1=0;
  if(fMCNDEventAllPrimTrackMeanPtTrueMult1) delete fMCNDEventAllPrimTrackMeanPtTrueMult1; fMCNDEventAllPrimTrackMeanPtTrueMult1=0;
  if(fMCNSDEventAllPrimTrackMeanPtTrueMult1) delete fMCNSDEventAllPrimTrackMeanPtTrueMult1; fMCNSDEventAllPrimTrackMeanPtTrueMult1=0;
  if(fMCTriggerPrimTrackMeanPtTrueMult1) delete fMCTriggerPrimTrackMeanPtTrueMult1; fMCTriggerPrimTrackMeanPtTrueMult1=0;
  if(fMCEventPrimTrackMeanPtTrueMult1) delete fMCEventPrimTrackMeanPtTrueMult1; fMCEventPrimTrackMeanPtTrueMult1=0;

  if(fCosmicsHisto) delete fCosmicsHisto; fCosmicsHisto=0;
  if(fEventCount) delete fEventCount; fEventCount=0;

  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    if(fRecTrackHist1[i]) delete fRecTrackHist1[i]; fRecTrackHist1[i]=0;
  }

  for(Int_t i=0; i<8; i++) { 
    if(fCorrRecTrackMultHist1[i]) delete fCorrRecTrackMultHist1[i]; fCorrRecTrackMultHist1[i]=0;
    if(fCorrRecTrackTrueMultHist1[i]) delete fCorrRecTrackTrueMultHist1[i]; fCorrRecTrackTrueMultHist1[i]=0;
    if(fCorrRecTrackTrueMultHist2[i]) delete fCorrRecTrackTrueMultHist2[i]; fCorrRecTrackTrueMultHist2[i]=0;
    if(fCorrRecTrackMeanPtMultHist1[i]) delete fCorrRecTrackMeanPtMultHist1[i]; fCorrRecTrackMeanPtMultHist1[i]=0;
    if(fCorrRecTrackMeanPtTrueMultHist1[i]) delete fCorrRecTrackMeanPtTrueMultHist1[i]; fCorrRecTrackMeanPtTrueMultHist1[i]=0;
    if(fCorrRecTrackPt1[i]) delete fCorrRecTrackPt1[i]; fCorrRecTrackPt1[i]=0;
  }

  for(Int_t i=0; i<5; i++) { 
    if(fCorrRecEventHist1[i]) delete fCorrRecEventHist1[i]; fCorrRecEventHist1[i]=0;
    if(fCorrRecEventHist2[i]) delete fCorrRecEventHist2[i]; fCorrRecEventHist2[i]=0;
  }

  if(fCorrectionFolder) delete fCorrectionFolder; fCorrectionFolder=0;
}

//_____________________________________________________________________________
void AlidNdPtCorrection::Init(){
  //
  // Init histograms
  //
  const Int_t etaNbins = 30; 
  const Int_t zvNbins = 12;

  // UA1 bining
  //const Int_t ptNbins = 52; 
  //Double_t binsPt[ptNbins+1] = { 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.20, 4.40, 4.60, 4.80, 5.00, 5.20, 5.40, 5.60, 5.80, 6.00, 7.00, 7.60, 8.80, 9.60 }; 

 /*
  const Int_t ptNbins = 56; 
  Double_t binsPt[ptNbins+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};
  Double_t binsEta[etaNbins+1] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
  Double_t binsZv[zvNbins+1] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
  */

  const Int_t ptNbins = 55; 
  Double_t binsPt[ptNbins+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};
  Double_t binsEta[etaNbins+1] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
  Double_t binsZv[zvNbins+1] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};


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
  Int_t binsRecEventMultHist1[2]={150,150};
  Double_t minRecEventMultHist1[2]={-0.5,-0.5}; 
  Double_t maxRecEventMultHist1[2]={149.5,149.5}; 
  fRecEventMultHist1 = new THnSparseF("fRecEventMultHist1","track multiplicity:tracklet multiplicity",2,binsRecEventMultHist1,minRecEventMultHist1,maxRecEventMultHist1);
  fRecEventMultHist1->GetAxis(0)->SetTitle("track_mult");
  fRecEventMultHist1->GetAxis(1)->SetTitle("tracklet_mult");
  fRecEventMultHist1->Sumw2();

  //
  char name[256];
  char title[256];

  Int_t binsMCAllPrimTrackMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCAllPrimTrackMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCAllPrimTrackMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCAllPrimTrackMultHist1");
  sprintf(title,"mcPt:mcEta:multiplicity");
  
  fMCAllPrimTrackMultHist1 = new THnSparseF(name,title,3,binsMCAllPrimTrackMultHist1,minMCAllPrimTrackMultHist1,maxMCAllPrimTrackMultHist1);
  fMCAllPrimTrackMultHist1->SetBinEdges(0,binsPt);
  fMCAllPrimTrackMultHist1->SetBinEdges(1,binsEta);
  fMCAllPrimTrackMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCAllPrimTrackMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCAllPrimTrackMultHist1->GetAxis(2)->SetTitle("multiplicity");
  fMCAllPrimTrackMultHist1->Sumw2();

  Int_t binsMCNDEventAllPrimTrackMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCNDEventAllPrimTrackMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCNDEventAllPrimTrackMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCNDEventAllPrimTrackMultHist1");
  sprintf(title,"mcPt:mcEta:multiplicity");
  
  fMCNDEventAllPrimTrackMultHist1 = new THnSparseF(name,title,3,binsMCNDEventAllPrimTrackMultHist1,minMCNDEventAllPrimTrackMultHist1,maxMCNDEventAllPrimTrackMultHist1);
  fMCNDEventAllPrimTrackMultHist1->SetBinEdges(0,binsPt);
  fMCNDEventAllPrimTrackMultHist1->SetBinEdges(1,binsEta);
  fMCNDEventAllPrimTrackMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCNDEventAllPrimTrackMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCNDEventAllPrimTrackMultHist1->GetAxis(2)->SetTitle("multiplicity");
  fMCNDEventAllPrimTrackMultHist1->Sumw2();

  Int_t binsMCNSDEventAllPrimTrackMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCNSDEventAllPrimTrackMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCNSDEventAllPrimTrackMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCNSDEventAllPrimTrackMultHist1");
  sprintf(title,"mcPt:mcEta:multiplicity");
  
  fMCNSDEventAllPrimTrackMultHist1 = new THnSparseF(name,title,3,binsMCNSDEventAllPrimTrackMultHist1,minMCNSDEventAllPrimTrackMultHist1,maxMCNSDEventAllPrimTrackMultHist1);
  fMCNSDEventAllPrimTrackMultHist1->SetBinEdges(0,binsPt);
  fMCNSDEventAllPrimTrackMultHist1->SetBinEdges(1,binsEta);
  fMCNSDEventAllPrimTrackMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCNSDEventAllPrimTrackMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCNSDEventAllPrimTrackMultHist1->GetAxis(2)->SetTitle("multiplicity");
  fMCNSDEventAllPrimTrackMultHist1->Sumw2();

  Int_t binsMCEventTriggerPrimTrackMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCEventTriggerPrimTrackMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCEventTriggerPrimTrackMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCTriggerPrimTrackMultHist1");
  sprintf(title,"mcPt:mcEta:multiplicity");
  
  fMCTriggerPrimTrackMultHist1 = new THnSparseF(name,title,3,binsMCEventTriggerPrimTrackMultHist1,minMCEventTriggerPrimTrackMultHist1,maxMCEventTriggerPrimTrackMultHist1);
  fMCTriggerPrimTrackMultHist1->SetBinEdges(0,binsPt);
  fMCTriggerPrimTrackMultHist1->SetBinEdges(1,binsEta);
  fMCTriggerPrimTrackMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCTriggerPrimTrackMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCTriggerPrimTrackMultHist1->GetAxis(2)->SetTitle("multiplicity");
  fMCTriggerPrimTrackMultHist1->Sumw2();

  Int_t binsMCEventPrimTrackMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCEventPrimTrackMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCEventPrimTrackMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCEventPrimTrackMultHist1");
  sprintf(title,"mcPt:mcEta:multiplicity");
  
  fMCEventPrimTrackMultHist1 = new THnSparseF(name,title,3,binsMCEventPrimTrackMultHist1,minMCEventPrimTrackMultHist1,maxMCEventPrimTrackMultHist1);
  fMCEventPrimTrackMultHist1->SetBinEdges(0,binsPt);
  fMCEventPrimTrackMultHist1->SetBinEdges(1,binsEta);
  fMCEventPrimTrackMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCEventPrimTrackMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCEventPrimTrackMultHist1->GetAxis(2)->SetTitle("multiplicity");
  fMCEventPrimTrackMultHist1->Sumw2();

  //
  // true multiplicity
  //

  Int_t binsMCAllPrimTrackTrueMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCAllPrimTrackTrueMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCAllPrimTrackTrueMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCAllPrimTrackTrueMultHist1");
  sprintf(title,"mcPt:mcEta:true_mult");
  
  fMCAllPrimTrackTrueMultHist1 = new THnSparseF(name,title,3,binsMCAllPrimTrackTrueMultHist1,minMCAllPrimTrackTrueMultHist1,maxMCAllPrimTrackTrueMultHist1);
  fMCAllPrimTrackTrueMultHist1->SetBinEdges(0,binsPt);
  fMCAllPrimTrackTrueMultHist1->SetBinEdges(1,binsEta);
  fMCAllPrimTrackTrueMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCAllPrimTrackTrueMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCAllPrimTrackTrueMultHist1->GetAxis(2)->SetTitle("true_mult");
  fMCAllPrimTrackTrueMultHist1->Sumw2();

  Int_t binsMCNDEventAllPrimTrackTrueMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCNDEventAllPrimTrackTrueMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCNDEventAllPrimTrackTrueMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCNDEventAllPrimTrackTrueMultHist1");
  sprintf(title,"mcPt:mcEta:true_mult");
  
  fMCNDEventAllPrimTrackTrueMultHist1 = new THnSparseF(name,title,3,binsMCNDEventAllPrimTrackTrueMultHist1,minMCNDEventAllPrimTrackTrueMultHist1,maxMCNDEventAllPrimTrackTrueMultHist1);
  fMCNDEventAllPrimTrackTrueMultHist1->SetBinEdges(0,binsPt);
  fMCNDEventAllPrimTrackTrueMultHist1->SetBinEdges(1,binsEta);
  fMCNDEventAllPrimTrackTrueMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCNDEventAllPrimTrackTrueMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCNDEventAllPrimTrackTrueMultHist1->GetAxis(2)->SetTitle("true_mult");
  fMCNDEventAllPrimTrackTrueMultHist1->Sumw2();

  Int_t binsMCNSDEventAllPrimTrackTrueMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCNSDEventAllPrimTrackTrueMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCNSDEventAllPrimTrackTrueMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCNSDEventAllPrimTrackTrueMultHist1");
  sprintf(title,"mcPt:mcEta:true_mult");
  
  fMCNSDEventAllPrimTrackTrueMultHist1 = new THnSparseF(name,title,3,binsMCNSDEventAllPrimTrackTrueMultHist1,minMCNSDEventAllPrimTrackTrueMultHist1,maxMCNSDEventAllPrimTrackTrueMultHist1);
  fMCNSDEventAllPrimTrackTrueMultHist1->SetBinEdges(0,binsPt);
  fMCNSDEventAllPrimTrackTrueMultHist1->SetBinEdges(1,binsEta);
  fMCNSDEventAllPrimTrackTrueMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCNSDEventAllPrimTrackTrueMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCNSDEventAllPrimTrackTrueMultHist1->GetAxis(2)->SetTitle("true_mult");
  fMCNSDEventAllPrimTrackTrueMultHist1->Sumw2();

  Int_t binsMCEventTriggerPrimTrackTrueMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCEventTriggerPrimTrackTrueMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCEventTriggerPrimTrackTrueMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCTriggerPrimTrackTrueMultHist1");
  sprintf(title,"mcPt:mcEta:true_mult");
  
  fMCTriggerPrimTrackTrueMultHist1 = new THnSparseF(name,title,3,binsMCEventTriggerPrimTrackTrueMultHist1,minMCEventTriggerPrimTrackTrueMultHist1,maxMCEventTriggerPrimTrackTrueMultHist1);
  fMCTriggerPrimTrackTrueMultHist1->SetBinEdges(0,binsPt);
  fMCTriggerPrimTrackTrueMultHist1->SetBinEdges(1,binsEta);
  fMCTriggerPrimTrackTrueMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCTriggerPrimTrackTrueMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCTriggerPrimTrackTrueMultHist1->GetAxis(2)->SetTitle("true_mult");
  fMCTriggerPrimTrackTrueMultHist1->Sumw2();

  Int_t binsMCEventPrimTrackTrueMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minMCEventPrimTrackTrueMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxMCEventPrimTrackTrueMultHist1[3]={20.,1.,149.5}; 
  sprintf(name,"fMCEventPrimTrackTrueMultHist1");
  sprintf(title,"mcPt:mcEta:true_mult");
  
  fMCEventPrimTrackTrueMultHist1 = new THnSparseF(name,title,3,binsMCEventPrimTrackTrueMultHist1,minMCEventPrimTrackTrueMultHist1,maxMCEventPrimTrackTrueMultHist1);
  fMCEventPrimTrackTrueMultHist1->SetBinEdges(0,binsPt);
  fMCEventPrimTrackTrueMultHist1->SetBinEdges(1,binsEta);
  fMCEventPrimTrackTrueMultHist1->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCEventPrimTrackTrueMultHist1->GetAxis(1)->SetTitle("mcEta");
  fMCEventPrimTrackTrueMultHist1->GetAxis(2)->SetTitle("true_mult");
  fMCEventPrimTrackTrueMultHist1->Sumw2();

  //
  // mcPT vs multiplicity vs true multiplicity
  //

  Int_t binsMCAllPrimTrackTrueMultHist2[3]={ptNbins,150,150};
  Double_t minMCAllPrimTrackTrueMultHist2[3]={0.,-0.5,-0.5}; 
  Double_t maxMCAllPrimTrackTrueMultHist2[3]={20.,149.5,149.5}; 
  sprintf(name,"fMCAllPrimTrackTrueMultHist2");
  sprintf(title,"mcPt:mult:true_mult");
  
  fMCAllPrimTrackTrueMultHist2 = new THnSparseF(name,title,3,binsMCAllPrimTrackTrueMultHist2,minMCAllPrimTrackTrueMultHist2,maxMCAllPrimTrackTrueMultHist2);
  fMCAllPrimTrackTrueMultHist2->SetBinEdges(0,binsPt);
  fMCAllPrimTrackTrueMultHist2->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCAllPrimTrackTrueMultHist2->GetAxis(1)->SetTitle("mult");
  fMCAllPrimTrackTrueMultHist2->GetAxis(2)->SetTitle("true_mult");
  fMCAllPrimTrackTrueMultHist2->Sumw2();

  Int_t binsMCNDEventAllPrimTrackTrueMultHist2[3]={ptNbins,150,150};
  Double_t minMCNDEventAllPrimTrackTrueMultHist2[3]={0.,-0.5,-0.5}; 
  Double_t maxMCNDEventAllPrimTrackTrueMultHist2[3]={20.,149.5,149.5}; 
  sprintf(name,"fMCNDEventAllPrimTrackTrueMultHist2");
  sprintf(title,"mcPt:mult:true_mult");
  
  fMCNDEventAllPrimTrackTrueMultHist2 = new THnSparseF(name,title,3,binsMCNDEventAllPrimTrackTrueMultHist2,minMCNDEventAllPrimTrackTrueMultHist2,maxMCNDEventAllPrimTrackTrueMultHist2);
  fMCNDEventAllPrimTrackTrueMultHist2->SetBinEdges(0,binsPt);
  fMCNDEventAllPrimTrackTrueMultHist2->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCNDEventAllPrimTrackTrueMultHist2->GetAxis(1)->SetTitle("mult");
  fMCNDEventAllPrimTrackTrueMultHist2->GetAxis(2)->SetTitle("true_mult");
  fMCNDEventAllPrimTrackTrueMultHist2->Sumw2();

  Int_t binsMCNSDEventAllPrimTrackTrueMultHist2[3]={ptNbins,150,150};
  Double_t minMCNSDEventAllPrimTrackTrueMultHist2[3]={0.,-0.5,-0.5}; 
  Double_t maxMCNSDEventAllPrimTrackTrueMultHist2[3]={20.,149.5,149.5}; 
  sprintf(name,"fMCNSDEventAllPrimTrackTrueMultHist2");
  sprintf(title,"mcPt:mult:true_mult");
  
  fMCNSDEventAllPrimTrackTrueMultHist2 = new THnSparseF(name,title,3,binsMCNSDEventAllPrimTrackTrueMultHist2,minMCNSDEventAllPrimTrackTrueMultHist2,maxMCNSDEventAllPrimTrackTrueMultHist2);
  fMCNSDEventAllPrimTrackTrueMultHist2->SetBinEdges(0,binsPt);
  fMCNSDEventAllPrimTrackTrueMultHist2->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCNSDEventAllPrimTrackTrueMultHist2->GetAxis(1)->SetTitle("mult");
  fMCNSDEventAllPrimTrackTrueMultHist2->GetAxis(2)->SetTitle("true_mult");
  fMCNSDEventAllPrimTrackTrueMultHist2->Sumw2();

  Int_t binsMCEventTriggerPrimTrackTrueMultHist2[3]={ptNbins,150,150};
  Double_t minMCEventTriggerPrimTrackTrueMultHist2[3]={0.,-0.5,-0.5}; 
  Double_t maxMCEventTriggerPrimTrackTrueMultHist2[3]={20.,149.5,149.5}; 
  sprintf(name,"fMCTriggerPrimTrackTrueMultHist2");
  sprintf(title,"mcPt:mult:true_mult");
  
  fMCTriggerPrimTrackTrueMultHist2 = new THnSparseF(name,title,3,binsMCEventTriggerPrimTrackTrueMultHist2,minMCEventTriggerPrimTrackTrueMultHist2,maxMCEventTriggerPrimTrackTrueMultHist2);
  fMCTriggerPrimTrackTrueMultHist2->SetBinEdges(0,binsPt);
  fMCTriggerPrimTrackTrueMultHist2->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCTriggerPrimTrackTrueMultHist2->GetAxis(1)->SetTitle("mult");
  fMCTriggerPrimTrackTrueMultHist2->GetAxis(2)->SetTitle("true_mult");
  fMCTriggerPrimTrackTrueMultHist2->Sumw2();

  Int_t binsMCEventPrimTrackTrueMultHist2[3]={ptNbins,150,150};
  Double_t minMCEventPrimTrackTrueMultHist2[3]={0.,-0.5,-0.5}; 
  Double_t maxMCEventPrimTrackTrueMultHist2[3]={20.,149.5,149.5}; 
  sprintf(name,"fMCEventPrimTrackTrueMultHist2");
  sprintf(title,"mcPt:mult:true_mult");
  
  fMCEventPrimTrackTrueMultHist2 = new THnSparseF(name,title,3,binsMCEventPrimTrackTrueMultHist2,minMCEventPrimTrackTrueMultHist2,maxMCEventPrimTrackTrueMultHist2);
  fMCEventPrimTrackTrueMultHist2->SetBinEdges(0,binsPt);
  fMCEventPrimTrackTrueMultHist2->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCEventPrimTrackTrueMultHist2->GetAxis(1)->SetTitle("mult");
  fMCEventPrimTrackTrueMultHist2->GetAxis(2)->SetTitle("true_mult");
  fMCEventPrimTrackTrueMultHist2->Sumw2();


  //
  // mean pt
  //
  Int_t binsMCAllPrimTrackMeanPtTrueMult1[2]={100,150};
  Double_t minMCAllPrimTrackMeanPtTrueMult1[2]={0.,-0.5}; 
  Double_t maxMCAllPrimTrackMeanPtTrueMult1[2]={10.,149.5}; 
  sprintf(name,"fMCAllPrimTrackMeanPtTrueMult1");
  sprintf(title,"event <mcPt>:true_mult");
  
  fMCAllPrimTrackMeanPtTrueMult1 = new THnSparseF(name,title,2,binsMCAllPrimTrackMeanPtTrueMult1,minMCAllPrimTrackMeanPtTrueMult1,maxMCAllPrimTrackMeanPtTrueMult1);
  fMCAllPrimTrackMeanPtTrueMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCAllPrimTrackMeanPtTrueMult1->GetAxis(1)->SetTitle("true_mult");
  fMCAllPrimTrackMeanPtTrueMult1->Sumw2();

  Int_t binsMCAllPrimTrackMeanPtMult1[2]={100,150};
  Double_t minMCAllPrimTrackMeanPtMult1[2]={0.,-0.5}; 
  Double_t maxMCAllPrimTrackMeanPtMult1[2]={10.,149.5}; 
  sprintf(name,"fMCAllPrimTrackMeanPtMult1");
  sprintf(title,"event <mcPt>:mult");
  
  fMCAllPrimTrackMeanPtMult1 = new THnSparseF(name,title,2,binsMCAllPrimTrackMeanPtMult1,minMCAllPrimTrackMeanPtMult1,maxMCAllPrimTrackMeanPtMult1);
  fMCAllPrimTrackMeanPtMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCAllPrimTrackMeanPtMult1->GetAxis(1)->SetTitle("multiplicity");
  fMCAllPrimTrackMeanPtMult1->Sumw2();

  //
  Int_t binsMCNDEventAllPrimTrackMeanPtTrueMult1[2]={100,150};
  Double_t minMCNDEventAllPrimTrackMeanPtTrueMult1[2]={0.,-0.5}; 
  Double_t maxMCNDEventAllPrimTrackMeanPtTrueMult1[2]={10.,149.5}; 
  sprintf(name,"fMCNDEventAllPrimTrackMeanPtTrueMult1");
  sprintf(title,"event <mcPt>:true_mult");
  
  fMCNDEventAllPrimTrackMeanPtTrueMult1 = new THnSparseF(name,title,2,binsMCNDEventAllPrimTrackMeanPtTrueMult1,minMCNDEventAllPrimTrackMeanPtTrueMult1,maxMCNDEventAllPrimTrackMeanPtTrueMult1);
  fMCNDEventAllPrimTrackMeanPtTrueMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCNDEventAllPrimTrackMeanPtTrueMult1->GetAxis(1)->SetTitle("true_mult");
  fMCNDEventAllPrimTrackMeanPtTrueMult1->Sumw2();

  Int_t binsMCNDEventAllPrimTrackMeanPtMult1[2]={100,150};
  Double_t minMCNDEventAllPrimTrackMeanPtMult1[2]={0.,-0.5}; 
  Double_t maxMCNDEventAllPrimTrackMeanPtMult1[2]={10.,149.5}; 
  sprintf(name,"fMCNDEventAllPrimTrackMeanPtMult1");
  sprintf(title,"event <mcPt>:mult");
  
  fMCNDEventAllPrimTrackMeanPtMult1 = new THnSparseF(name,title,2,binsMCNDEventAllPrimTrackMeanPtMult1,minMCNDEventAllPrimTrackMeanPtMult1,maxMCNDEventAllPrimTrackMeanPtMult1);
  fMCNDEventAllPrimTrackMeanPtMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCNDEventAllPrimTrackMeanPtMult1->GetAxis(1)->SetTitle("multiplicity");
  fMCNDEventAllPrimTrackMeanPtMult1->Sumw2();

  //
  Int_t binsMCNSDEventAllPrimTrackMeanPtTrueMult1[2]={100,150};
  Double_t minMCNSDEventAllPrimTrackMeanPtTrueMult1[2]={0.,-0.5}; 
  Double_t maxMCNSDEventAllPrimTrackMeanPtTrueMult1[2]={10.,149.5}; 
  sprintf(name,"fMCNSDEventAllPrimTrackMeanPtTrueMult1");
  sprintf(title,"event <mcPt>:true_mult");
  
  fMCNSDEventAllPrimTrackMeanPtTrueMult1 = new THnSparseF(name,title,2,binsMCNSDEventAllPrimTrackMeanPtTrueMult1,minMCNSDEventAllPrimTrackMeanPtTrueMult1,maxMCNSDEventAllPrimTrackMeanPtTrueMult1);
  fMCNSDEventAllPrimTrackMeanPtTrueMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCNSDEventAllPrimTrackMeanPtTrueMult1->GetAxis(1)->SetTitle("true_mult");
  fMCNSDEventAllPrimTrackMeanPtTrueMult1->Sumw2();

  Int_t binsMCNSDEventAllPrimTrackMeanPtMult1[2]={100,150};
  Double_t minMCNSDEventAllPrimTrackMeanPtMult1[2]={0.,-0.5}; 
  Double_t maxMCNSDEventAllPrimTrackMeanPtMult1[2]={10.,149.5}; 
  sprintf(name,"fMCNSDEventAllPrimTrackMeanPtMult1");
  sprintf(title,"event <mcPt>:mult");
  
  fMCNSDEventAllPrimTrackMeanPtMult1 = new THnSparseF(name,title,2,binsMCNSDEventAllPrimTrackMeanPtMult1,minMCNSDEventAllPrimTrackMeanPtMult1,maxMCNSDEventAllPrimTrackMeanPtMult1);
  fMCNSDEventAllPrimTrackMeanPtMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCNSDEventAllPrimTrackMeanPtMult1->GetAxis(1)->SetTitle("multiplicity");
  fMCNSDEventAllPrimTrackMeanPtMult1->Sumw2();

  //
  Int_t binsMCTriggerPrimTrackMeanPtTrueMult1[2]={100,150};
  Double_t minMCTriggerPrimTrackMeanPtTrueMult1[2]={0.,-0.5}; 
  Double_t maxMCTriggerPrimTrackMeanPtTrueMult1[2]={10.,149.5}; 
  sprintf(name,"fMCTriggerPrimTrackMeanPtTrueMult1");
  sprintf(title,"event <mcPt>:true_mult");
  
  fMCTriggerPrimTrackMeanPtTrueMult1 = new THnSparseF(name,title,2,binsMCTriggerPrimTrackMeanPtTrueMult1,minMCTriggerPrimTrackMeanPtTrueMult1,maxMCTriggerPrimTrackMeanPtTrueMult1);
  fMCTriggerPrimTrackMeanPtTrueMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCTriggerPrimTrackMeanPtTrueMult1->GetAxis(1)->SetTitle("true_mult");
  fMCTriggerPrimTrackMeanPtTrueMult1->Sumw2();

  Int_t binsMCTriggerPrimTrackMeanPtMult1[2]={100,150};
  Double_t minMCTriggerPrimTrackMeanPtMult1[2]={0.,-0.5}; 
  Double_t maxMCTriggerPrimTrackMeanPtMult1[2]={10.,149.5}; 
  sprintf(name,"fMCTriggerPrimTrackMeanPtMult1");
  sprintf(title,"event <mcPt>:mult");
  
  fMCTriggerPrimTrackMeanPtMult1 = new THnSparseF(name,title,2,binsMCTriggerPrimTrackMeanPtMult1,minMCTriggerPrimTrackMeanPtMult1,maxMCTriggerPrimTrackMeanPtMult1);
  fMCTriggerPrimTrackMeanPtMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCTriggerPrimTrackMeanPtMult1->GetAxis(1)->SetTitle("multiplicity");
  fMCTriggerPrimTrackMeanPtMult1->Sumw2();

  //
  Int_t binsMCEventPrimTrackMeanPtTrueMult1[2]={100,150};
  Double_t minMCEventPrimTrackMeanPtTrueMult1[2]={0.,-0.5}; 
  Double_t maxMCEventPrimTrackMeanPtTrueMult1[2]={10.,149.5}; 
  sprintf(name,"fMCEventPrimTrackMeanPtTrueMult1");
  sprintf(title,"event <mcPt>:true_mult");
  
  fMCEventPrimTrackMeanPtTrueMult1 = new THnSparseF(name,title,2,binsMCEventPrimTrackMeanPtTrueMult1,minMCEventPrimTrackMeanPtTrueMult1,maxMCEventPrimTrackMeanPtTrueMult1);
  fMCEventPrimTrackMeanPtTrueMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCEventPrimTrackMeanPtTrueMult1->GetAxis(1)->SetTitle("true_mult");
  fMCEventPrimTrackMeanPtTrueMult1->Sumw2();

  Int_t binsMCEventPrimTrackMeanPtMult1[2]={100,150};
  Double_t minMCEventPrimTrackMeanPtMult1[2]={0.,-0.5}; 
  Double_t maxMCEventPrimTrackMeanPtMult1[2]={10.,149.5}; 
  sprintf(name,"fMCEventPrimTrackMeanPtMult1");
  sprintf(title,"event <mcPt>:mult");
  
  fMCEventPrimTrackMeanPtMult1 = new THnSparseF(name,title,2,binsMCEventPrimTrackMeanPtMult1,minMCEventPrimTrackMeanPtMult1,maxMCEventPrimTrackMeanPtMult1);
  fMCEventPrimTrackMeanPtMult1->GetAxis(0)->SetTitle("<mcPt> (GeV/c)");
  fMCEventPrimTrackMeanPtMult1->GetAxis(1)->SetTitle("multiplicity");
  fMCEventPrimTrackMeanPtMult1->Sumw2();






  //
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) 
  {
    // THnSparse track histograms
    //
    Int_t binsRecTrackHist1[3]={ptNbins,etaNbins,90};
    Double_t minRecTrackHist1[3]={0.,-1.,0.}; 
    Double_t maxRecTrackHist1[3]={10.,1.,2.*TMath::Pi()};
    sprintf(name,"fRecTrackHist1_%d",i);
    sprintf(title,"Pt:Eta:Phi");
  
    fRecTrackHist1[i] = new THnSparseF(name,title,3,binsRecTrackHist1,minRecTrackHist1,maxRecTrackHist1);
    fRecTrackHist1[i]->SetBinEdges(0,binsPt);
    fRecTrackHist1[i]->SetBinEdges(1,binsEta);
    fRecTrackHist1[i]->GetAxis(0)->SetTitle("Pt (GeV/c)");
    fRecTrackHist1[i]->GetAxis(1)->SetTitle("Eta");
    fRecTrackHist1[i]->GetAxis(2)->SetTitle("Phi (rad)");
    fRecTrackHist1[i]->Sumw2();
  }

  //
  Int_t binsCorrRecTrackMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minCorrRecTrackMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxCorrRecTrackMultHist1[3]={20.,1.,149.5};

  Int_t binsCorrRecTrackTrueMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minCorrRecTrackTrueMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxCorrRecTrackTrueMultHist1[3]={20.,1.,149.5};

  Int_t binsCorrRecTrackTrueMultHist2[3]={ptNbins,150,150};
  Double_t minCorrRecTrackTrueMultHist2[3]={0.,-0.5,-0.5}; 
  Double_t maxCorrRecTrackTrueMultHist2[3]={20.,149.5,149.5};

  //
  Int_t binsCorrRecTrackMeanPtMultHist1[2]={100,150};
  Double_t minCorrRecTrackMeanPtMultHist1[2]={0.,-0.5}; 
  Double_t maxCorrRecTrackMeanPtMultHist1[2]={10.,149.5};

  Int_t binsCorrRecTrackMeanPtTrueMultHist1[2]={100,150};
  Double_t minCorrRecTrackMeanPtTrueMultHist1[2]={0.,-0.5}; 
  Double_t maxCorrRecTrackMeanPtTrueMultHist1[2]={10.,149.5};

  Int_t binsCorrRecTrackPt1[1]={200};
  Double_t minCorrRecTrackPt1[1]={0.}; 
  Double_t maxCorrRecTrackPt1[1]={10.};

  for(Int_t i=0; i<8; i++) 
  {
    // THnSparse track histograms
    sprintf(name,"fCorrRecTrackMultHist1_%d",i);
    sprintf(title,"Pt:Eta:mult");
    fCorrRecTrackMultHist1[i] = new THnSparseF(name,title,3,binsCorrRecTrackMultHist1,minCorrRecTrackMultHist1,maxCorrRecTrackMultHist1);
    fCorrRecTrackMultHist1[i]->SetBinEdges(0,binsPt);
    fCorrRecTrackMultHist1[i]->SetBinEdges(1,binsEta);
    fCorrRecTrackMultHist1[i]->GetAxis(0)->SetTitle("Pt (GeV/c)");
    fCorrRecTrackMultHist1[i]->GetAxis(1)->SetTitle("Eta");
    fCorrRecTrackMultHist1[i]->GetAxis(2)->SetTitle("multiplicity");
    fCorrRecTrackMultHist1[i]->Sumw2();

    // THnSparse track histograms
    sprintf(name,"fCorrRecTrackTrueMultHist1_%d",i);
    sprintf(title,"Pt:Eta:true_mult");
    fCorrRecTrackTrueMultHist1[i] = new THnSparseF(name,title,3,binsCorrRecTrackTrueMultHist1,minCorrRecTrackTrueMultHist1,maxCorrRecTrackTrueMultHist1);
    fCorrRecTrackTrueMultHist1[i]->SetBinEdges(0,binsPt);
    fCorrRecTrackTrueMultHist1[i]->SetBinEdges(1,binsEta);
    fCorrRecTrackTrueMultHist1[i]->GetAxis(0)->SetTitle("Pt (GeV/c)");
    fCorrRecTrackTrueMultHist1[i]->GetAxis(1)->SetTitle("Eta");
    fCorrRecTrackTrueMultHist1[i]->GetAxis(2)->SetTitle("true multiplicity");
    fCorrRecTrackTrueMultHist1[i]->Sumw2();

    //
    sprintf(name,"fCorrRecTrackTrueMultHist2_%d",i);
    sprintf(title,"Pt:mult:true_mult");
    fCorrRecTrackTrueMultHist2[i] = new THnSparseF(name,title,3,binsCorrRecTrackTrueMultHist2,minCorrRecTrackTrueMultHist2,maxCorrRecTrackTrueMultHist2);
    fCorrRecTrackTrueMultHist2[i]->SetBinEdges(0,binsPt);
    fCorrRecTrackTrueMultHist2[i]->GetAxis(0)->SetTitle("Pt (GeV/c)");
    fCorrRecTrackTrueMultHist2[i]->GetAxis(1)->SetTitle("multiplicity");
    fCorrRecTrackTrueMultHist2[i]->GetAxis(2)->SetTitle("true multiplicity");
    fCorrRecTrackTrueMultHist2[i]->Sumw2();

    // THnSparse track histograms
    sprintf(name,"fCorrRecTrackMeanPtMultHist1_%d",i);
    sprintf(title,"<Pt>:mult");
    fCorrRecTrackMeanPtMultHist1[i] = new THnSparseF(name,title,2,binsCorrRecTrackMeanPtMultHist1,minCorrRecTrackMeanPtMultHist1,maxCorrRecTrackMeanPtMultHist1);
    fCorrRecTrackMeanPtMultHist1[i]->GetAxis(0)->SetTitle("<Pt> (GeV/c)");
    fCorrRecTrackMeanPtMultHist1[i]->GetAxis(1)->SetTitle("multiplicity");
    fCorrRecTrackMeanPtMultHist1[i]->Sumw2();

    // THnSparse track histograms
    sprintf(name,"fCorrRecTrackMeanPtTrueMultHist1_%d",i);
    sprintf(title,"<Pt>:true_mult");
    fCorrRecTrackMeanPtTrueMultHist1[i] = new THnSparseF(name,title,2,binsCorrRecTrackMeanPtTrueMultHist1,minCorrRecTrackMeanPtTrueMultHist1,maxCorrRecTrackMeanPtTrueMultHist1);
    fCorrRecTrackMeanPtTrueMultHist1[i]->GetAxis(0)->SetTitle("<Pt> (GeV/c)");
    fCorrRecTrackMeanPtTrueMultHist1[i]->GetAxis(1)->SetTitle("true multiplicity");
    fCorrRecTrackMeanPtTrueMultHist1[i]->Sumw2();

    sprintf(name,"fCorrRecTrackPt1_%d",i);
    sprintf(title,"pt small bining");
    fCorrRecTrackPt1[i] = new THnSparseF(name,title,1,binsCorrRecTrackPt1,minCorrRecTrackPt1,maxCorrRecTrackPt1);
    fCorrRecTrackPt1[i]->GetAxis(0)->SetTitle("Pt (GeV/c)");
    fCorrRecTrackPt1[i]->Sumw2();

  }

  Int_t binsEventMatrix[2]={zvNbins,150};
  Double_t minEventMatrix[2]={-25.,-0.5};
  Double_t maxEventMatrix[2]={25.,149.5};

  fMCAllEventMultHist1 = new THnSparseF("fMCAllEventMultHist1","mcZv:mult",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fMCAllEventMultHist1->SetBinEdges(0,binsZv);
  fMCAllEventMultHist1->GetAxis(0)->SetTitle("mcZv (cm)");
  fMCAllEventMultHist1->GetAxis(1)->SetTitle("multiplicity");
  fMCAllEventMultHist1->Sumw2();

  fMCAllNDEventMultHist1 = new THnSparseF("fMCAllNDEventMultHist1","mcZv:mult",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fMCAllNDEventMultHist1->SetBinEdges(0,binsZv);
  fMCAllNDEventMultHist1->GetAxis(0)->SetTitle("mcZv (cm)");
  fMCAllNDEventMultHist1->GetAxis(1)->SetTitle("multiplicity");
  fMCAllNDEventMultHist1->Sumw2();

  fMCAllNSDEventMultHist1 = new THnSparseF("fMCAllNSDEventMultHist1","mcZv:mult",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fMCAllNSDEventMultHist1->SetBinEdges(0,binsZv);
  fMCAllNSDEventMultHist1->GetAxis(0)->SetTitle("mcZv (cm)");
  fMCAllNSDEventMultHist1->GetAxis(1)->SetTitle("multiplicity");
  fMCAllNSDEventMultHist1->Sumw2();

  fMCTriggerMultHist1 = new THnSparseF("fMCTriggerMultHist1","mcZv:mult",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fMCTriggerMultHist1->SetBinEdges(0,binsZv);
  fMCTriggerMultHist1->GetAxis(0)->SetTitle("mcZv (cm)");
  fMCTriggerMultHist1->GetAxis(1)->SetTitle("multiplicity");
  fMCTriggerMultHist1->Sumw2();

  fMCEventMultHist1 = new THnSparseF("fMCEventMultHist1","mcZv:mult",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fMCEventMultHist1->SetBinEdges(0,binsZv);
  fMCEventMultHist1->GetAxis(0)->SetTitle("mcZv (cm)");
  fMCEventMultHist1->GetAxis(1)->SetTitle("multiplicity");
  fMCEventMultHist1->Sumw2();

  for(Int_t i=0; i<5; i++) 
  {
    // event corrected histograms
    sprintf(name,"fCorrRecEventHist1_%d",i);
    sprintf(title,"mcZv:mult");
    fCorrRecEventHist1[i] = new THnSparseF("fCorrRecEventHist1","mcZv:mult",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
    fCorrRecEventHist1[i]->SetBinEdges(0,binsZv);
    fCorrRecEventHist1[i]->GetAxis(0)->SetTitle("Zv (cm)");
    fCorrRecEventHist1[i]->GetAxis(1)->SetTitle("multiplicity");
    fCorrRecEventHist1[i]->Sumw2();

    // empty event corrected histograms
    sprintf(name,"fCorrRecEventHist2_%d",i);
    sprintf(title,"mcZv:mult");
    fCorrRecEventHist2[i] = new THnSparseF("fCorrRecEventHist2","mcZv:mult",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
    fCorrRecEventHist2[i]->SetBinEdges(0,binsZv);
    fCorrRecEventHist2[i]->GetAxis(0)->SetTitle("Zv (cm)");
    fCorrRecEventHist2[i]->GetAxis(1)->SetTitle("multiplicity");
    fCorrRecEventHist2[i]->Sumw2();
  }
  
  //
  // cosmics histo
  //
  Int_t binsCosmicsHisto[3]=  {151, 300, ptNbins};
  Double_t minCosmicsHisto[3]={-1.5, -2.*TMath::Pi(), 0.0}; 
  Double_t maxCosmicsHisto[3]={ 1.5, 2.*TMath::Pi(), 16.0}; 
  sprintf(name,"fCosmicsHisto");
  sprintf(title,"deta:dphi:pt");
  
  fCosmicsHisto = new THnSparseF(name,title,3,binsCosmicsHisto,minCosmicsHisto,maxCosmicsHisto);
  fCosmicsHisto->SetBinEdges(2,binsPt);
  fCosmicsHisto->GetAxis(0)->SetTitle("#Delta#eta");
  fCosmicsHisto->GetAxis(1)->SetTitle("#Delta#phi (rad)");
  fCosmicsHisto->GetAxis(2)->SetTitle("pt (GV/c)");
  fCosmicsHisto->Sumw2();

  //
  Int_t binsEventCount[2]={2,2};
  Double_t minEventCount[2]={0,0}; 
  Double_t maxEventCount[2]={2,2}; 
  fEventCount = new THnSparseF("fEventCount","trig vs trig+vertex",2,binsEventCount,minEventCount,maxEventCount);
  fEventCount->GetAxis(0)->SetTitle("trig");
  fEventCount->GetAxis(1)->SetTitle("trig+vert");
  fEventCount->Sumw2();


  // init output folder
  fCorrectionFolder = CreateFolder("folderdNdPt","Correction dNdPt Folder");

  // init correction matrices
  TFile *file = (TFile *)TFile::Open(fCorrMatrixFileName.Data()); 
  if(!file) { 
    AliDebug(AliLog::kError, "file with efficiency matrices not available");
    printf("file with efficiency matrices not available \n");
  } else {
    TFolder *folder = (TFolder *)file->FindObjectAny("folderdNdPt");
    if(!folder) { 
      AliDebug(AliLog::kError, "file without folderdNdPt");
      printf("file without folderdNdPt \n");
    } else {
      // rec. event mult vs true multiplicity 
      fEventMultCorrelationMatrix = (THnSparseF*)folder->FindObject("event_mult_correlation_matrix");
      if(!fEventMultCorrelationMatrix) {
         Printf("No %s matrix \n", "event_mult_correlation_matrix");
	 return;
      }

      //
      // event level corrections (zv,mult_MB)
      //
 
      // trigger bias correction (MBtoND) 
      fCorrTriggerMBtoNDEventMatrix = (THnSparseF*)folder->FindObject("zv_mult_trig_MBtoND_corr_matrix");
      if(!fCorrTriggerMBtoNDEventMatrix) {
         Printf("No %s matrix \n", "zv_mult_trig_MBtoND_corr_matrix");
	 return;
      }

      // trigger bias correction (MBtoNSD)
      fCorrTriggerMBtoNSDEventMatrix = (THnSparseF*)folder->FindObject("zv_mult_trig_MBtoNSD_corr_matrix");
      if(!fCorrTriggerMBtoNSDEventMatrix) {
         Printf("No %s matrix \n", "zv_mult_trig_MBtoNSD_corr_matrix");
	 return;
      }

      // trigger bias correction (MBtoInel)
      fCorrTriggerMBtoInelEventMatrix = (THnSparseF*)folder->FindObject("zv_mult_trig_MBtoInel_corr_matrix");
      if(!fCorrTriggerMBtoInelEventMatrix) {
         Printf("No %s matrix \n", "zv_mult_trig_MBtoInel_corr_matrix"); 
	 return;
      }
     
      // vertex reconstruction efficiency correction
      fCorrEventMatrix = (THnSparseF*)folder->FindObject("zv_mult_event_corr_matrix");
      if(!fCorrEventMatrix) {
         Printf("No %s matrix \n", "zv_mult_event_corr_matrix");
	 return;
      }

      //
      // histogram needed for empty events corrections
      //
      fZvNorm = (TH1D*)folder->FindObject("zv_distribution_norm");
      if(!fZvNorm) {
         Printf("No %s matrix \n", "fZvNorm");
	 return;
      }

      fZvEmptyEventsNorm = (TH1D*)folder->FindObject("zv_empty_events_norm");
      if(!fZvEmptyEventsNorm) {
         Printf("No %s matrix \n", "fZvEmptyEventsNorm");
	 return;
      }

      fLHCBin0Background = (TH1D*)folder->FindObject("hLHCBin0Background");
      if(!fLHCBin0Background) {
         Printf("No %s matrix \n", "fLHCBin0Background");
	 return;
      }

      //
      // track-event level corrections (zv,pt,eta)
      //

      // trigger bias correction (MBtoND) 
      fCorrTriggerMBtoNDTrackEventMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_trig_MBtoND_corr_matrix");
      if(!fCorrTriggerMBtoNDTrackEventMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_trig_MBtoND_corr_matrix");
	 return;
      }

      // trigger bias correction (MBtoNSD)
      fCorrTriggerMBtoNSDTrackEventMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_trig_MBtoNSD_corr_matrix");
      if(!fCorrTriggerMBtoNSDTrackEventMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_trig_MBtoNSD_corr_matrix");
	 return;
      }

      // trigger bias correction (MBtoInel) 
      fCorrTriggerMBtoInelTrackEventMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_trig_MBtoInel_corr_matrix");
      if(!fCorrTriggerMBtoInelTrackEventMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_trig_MBtoInel_corr_matrix");
	 return;
      }
    
      // vertex reconstruction efficiency correction (zv,pt,eta)
      fCorrTrackEventMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_event_corr_matrix");
      if(!fCorrTrackEventMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_event_corr_matrix");
	 return;
      }

      // track reconstruction efficiency correction (zv,pt,eta)
      fCorrTrackMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_corr_matrix");
      if(!fCorrTrackMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_corr_matrix");
	 return;
      }

      // high pt track reconstruction efficiency correction (zv,pt,eta)
      fCorrHighPtTrackMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_highPt_track_corr_matrix");
      if(!fCorrHighPtTrackMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_highPt_track_corr_matrix");
	 return;
      }

      // secondary tracks contamination correction (zv,pt,eta)
      fContTrackMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_cont_matrix");
      if(!fContTrackMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_cont_matrix");
	 return;
      }

      // multiply reconstructed tracks correction
      fContMultTrackMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_mult_track_cont_matrix");
      if(!fContMultTrackMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_mult_track_cont_matrix");
	 return;
      }
    }
  }

}

//_____________________________________________________________________________
void AlidNdPtCorrection::Process(AliESDEvent *esdEvent, AliMCEvent *mcEvent)
{
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
  AliPhysicsSelection *trigSel = NULL;
  AliTriggerAnalysis *trigAna = NULL; // needed for andV0

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
      if(GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtxUpdate || 
         GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt) 
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

  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);
  AliPWG0Helper::MCProcessType evtType = AliPWG0Helper::kInvalidProcess;
  Int_t multMCTrueTracks = 0;

  if(IsUseMCInfo())
  {
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
    //Printf("evtType %d \n", evtType);
    AliDebug(AliLog::kDebug+1, Form("Found process type %d", evtType));

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    // Fill MC event histogram
    Double_t vMCEventHist1[3]={vtxMC[0],vtxMC[1],vtxMC[2]};
    fMCEventHist1->Fill(vMCEventHist1);

    // multipliticy of all MC primary tracks
    // in Zvtx, eta ranges
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
  else if( GetAnalysisMode() == AlidNdPtHelper::kTPCITS ||  GetAnalysisMode() == AlidNdPtHelper::kTPCSPDvtx || 
           GetAnalysisMode()==AlidNdPtHelper::kTPCSPDvtxUpdate || GetAnalysisMode()==AlidNdPtHelper::kTPCITSHybrid ) 
  {
     const AliMultiplicity* mult = esdEvent->GetMultiplicity();
     //if(mult && vtxESD->GetStatus() && isRecVertex)
     if(mult)
       multMBTracks = mult->GetNumberOfTracklets();
    
  } 
  else if( GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtxUpdate || 
           GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt)
  {
     if(vtxESD->GetStatus() && isRecVertex)
       multMBTracks = vtxESD->GetNContributors();

  }
  else {
    AliDebug(AliLog::kError, Form("Found analysis type %s", GetAnalysisMode()));
    return; 
  }

  Bool_t isTrigAndVertex = isEventTriggered && isEventOK;
  Double_t vEventCount[2] = { isEventTriggered, isTrigAndVertex};
  fEventCount->Fill(vEventCount);

  //
  // correct event and track histograms
  //
  TObjArray *allChargedTracks=0;
  Int_t multRec=0, multRecTemp=0;
  Int_t *labelsRec=0;
  Bool_t isCosmic = kFALSE;


  if(isEventOK && isEventTriggered)
  {
    // get all charged tracks
    allChargedTracks = AlidNdPtHelper::GetAllChargedTracks(esdEvent,GetAnalysisMode());
    if(!allChargedTracks) return;

    Int_t entries = allChargedTracks->GetEntries();
    labelsRec = new Int_t[entries];

    // calculate mult of reconstructed tracks

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

      // cosmics analysis
      isCosmic = kFALSE;
      if( GetParticleMode()==AlidNdPtHelper::kCosmic )
      {
          for(Int_t j=0; j<entries;++j) 
          {
            AliESDtrack *track1 = (AliESDtrack*)allChargedTracks->At(j);
            if(!track1) continue;
            if(track1->Charge()==0) continue;

            if( esdTrackCuts->AcceptTrack(track) && accCuts->AcceptTrack(track) && 
	        esdTrackCuts->AcceptTrack(track1) && accCuts->AcceptTrack(track1) ) 
            { 
              isCosmic = AlidNdPtHelper::IsCosmicTrack(track, track1);
	    }
            if(isCosmic) 
	    {
	      Double_t vCosmicsHisto[3] = { track->Eta()+track1->Eta(), track->Phi()-track1->Phi(), track1->Pt() };
	      fCosmicsHisto->Fill(vCosmicsHisto);
	    }
	  }
         
        if(!isCosmic) continue;
      }

      if(esdTrackCuts->AcceptTrack(track)) 
      {
          if(accCuts->AcceptTrack(track)) multRecTemp++;
        /*
        if(GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt) {
          if(AlidNdPtHelper::IsGoodImpPar(track) && accCuts->AcceptTrack(track)) multRecTemp++;
        }
	else {
          if(accCuts->AcceptTrack(track)) multRecTemp++;
        }
	*/
      }  
    }

    /*
    // check multiplicity
    const AliMultiplicity* mult = esdEvent->GetMultiplicity();
    Int_t trackletMult = 0;
    if (mult) {
       for(Int_t i=0; i<mult->GetNumberOfTracklets(); i++) {
          if(TMath::Abs(mult->GetEta(i)) < accCuts->GetMaxEta() )
	    trackletMult++;
       }
    }
    // use tracklet multiplicity
    multRecTemp = trackletMult;
    */

    //
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
        
      // track-level corrections
      if(!esdTrackCuts->AcceptTrack(track))  continue;
      //if(GetAnalysisMode()==AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt && !AlidNdPtHelper::IsGoodImpPar(track)) continue;

        // cosmics analysis
        isCosmic = kFALSE;
        if( GetParticleMode()==AlidNdPtHelper::kCosmic )
        {
          for(Int_t j=0; j<entries;++j) 
          {
            AliESDtrack *track1 = (AliESDtrack*)allChargedTracks->At(j);
            if(!track1) continue;
            if(track1->Charge()==0) continue;

            if( esdTrackCuts->AcceptTrack(track) && accCuts->AcceptTrack(track) && 
	        esdTrackCuts->AcceptTrack(track1) && accCuts->AcceptTrack(track1) ) 
            { 
              isCosmic = AlidNdPtHelper::IsCosmicTrack(track, track1);
	    }
	  }
          if(!isCosmic) continue;
        }

        Bool_t isOK = kFALSE;
        Double_t x[3]; track->GetXYZ(x);
        Double_t b[3]; AliTracker::GetBxByBz(x,b);

        //
        // if TPC-ITS hybrid tracking (kTPCITSHybrid)
        // replace track parameters with TPC-ony track parameters
        //
        if( GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybrid || GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt) 
        {
          // Relate TPC-only tracks to Tracks or SPD vertex
          isOK = track->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig);
          if(!isOK) continue;

	  // replace esd track parameters with TPCinner
          AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(track->GetTPCInnerParam()));
	  if (!tpcTrack) return;
          track->Set(tpcTrack->GetX(),tpcTrack->GetAlpha(),tpcTrack->GetParameter(),tpcTrack->GetCovariance());

          if(tpcTrack) delete tpcTrack; 
        } 

        //
        if (GetAnalysisMode()==AlidNdPtHelper::kTPCSPDvtxUpdate || GetAnalysisMode()==AlidNdPtHelper::kTPCTrackSPDvtxUpdate) 
        {
	   //
	   // update track parameters
	   //
           AliExternalTrackParam cParam;
	   isOK = track->RelateToVertexTPC(vtxESD,esdEvent->GetMagneticField(),kVeryBig,&cParam);
	   if(!isOK) continue;
	   track->Set(cParam.GetX(),cParam.GetAlpha(),cParam.GetParameter(),cParam.GetCovariance());

           if(accCuts->AcceptTrack(track)) { 
             FillHistograms(track,stack,AlidNdPtHelper::kRecTracks,vtxESD->GetZv(),multRecTemp,multMCTrueTracks); 
	     labelsRec[multRec] = TMath::Abs(track->GetLabel());
	     multRec++;
	   }
         }
         else if (GetAnalysisMode()==AlidNdPtHelper::kMCRec && IsUseMCInfo()) 
         { 
	   //
	   // Replace rec with MC
	   //
           if(accCuts->AcceptTrack(track)) { 
	     FillHistograms(track,stack,AlidNdPtHelper::kRecTracks,vtxMC[2],multRecTemp, multMCTrueTracks); 
	     labelsRec[multRec] = TMath::Abs(track->GetLabel());
	     multRec++;
           }
	 } 
	 else  {
	   //
	   // all the rest tracking scenarios
	   //
           if(accCuts->AcceptTrack(track)) { 
             FillHistograms(track,stack,AlidNdPtHelper::kRecTracks,vtxESD->GetZv(),multRecTemp, multMCTrueTracks); 
	     labelsRec[multRec] = TMath::Abs(track->GetLabel());
	     multRec++;
	   }
         }
      }

    // event-level corrections
    if(GetAnalysisMode()==AlidNdPtHelper::kMCRec && IsUseMCInfo()) 
    { 
      FillHistograms(AlidNdPtHelper::kRecEvents,vtxMC[2],multMBTracks);
    }
    else {
      FillHistograms(AlidNdPtHelper::kRecEvents,vtxESD->GetZv(),multMBTracks);
    }

    // calculate meanPt from the event
    Double_t meanPtMult[8] = {0};  
    for (Int_t i = 0; i<8; i++) {
      if(!fCorrRecTrackMultHist1[i]) continue;
      TH1D *hp = (TH1D *)fCorrRecTrackPt1[i]->Projection(0);
      if(!hp) continue;
      meanPtMult[i] = hp->GetMean();    
      Double_t vCorrRecTrackMeanPtMultHist1[2] = {meanPtMult[i],multRecTemp};
      fCorrRecTrackMeanPtMultHist1[i]->Fill(vCorrRecTrackMeanPtMultHist1); 
      
      if( IsUseMCInfo() ) {
        Double_t vCorrRecTrackMeanPtTrueMultHist1[2] = {meanPtMult[i],multMCTrueTracks};
        fCorrRecTrackMeanPtTrueMultHist1[i]->Fill(vCorrRecTrackMeanPtTrueMultHist1); 
      }

      // reset pt histo for the next event
      if(fCorrRecTrackPt1[i])  fCorrRecTrackPt1[i]->Reset();
      if(hp) delete hp;
    }

    // control event histograms
    Double_t vRecEventHist1[3] = {vtxESD->GetXv(),vtxESD->GetYv(),vtxESD->GetZv()};
    fRecEventHist1->Fill(vRecEventHist1);

    // correlation track multiplicity vs MB track multiplicity
    Double_t vRecEventMultHist1[3] = {multRec, multMBTracks};
    fRecEventMultHist1->Fill(vRecEventMultHist1);
  }

  // empty events corrections
  // no reconstructed zv
  if( isEventTriggered && multMBTracks==0 ) 
  {
    if(GetAnalysisMode()==AlidNdPtHelper::kMCRec && IsUseMCInfo()) 
    {
      FillHistograms(AlidNdPtHelper::kTriggeredEvents,vtxMC[2],multMBTracks);
    }
    else {
      Double_t zv = fZvNorm->GetRandom();
      if(zv>evtCuts->GetMinZv() && zv<evtCuts->GetMaxZv())
        FillHistograms(AlidNdPtHelper::kTriggeredEvents,zv,multMBTracks);
    }
  }

  if(IsUseMCInfo())  
  {
    // select MC events 
    if(evtCuts->AcceptMCEvent(mcEvent))
    {
      //
      // event histograms
      //
      Double_t vMCEventMatrix[2] = {vtxMC[2],multMBTracks};
      fMCAllEventMultHist1->Fill(vMCEventMatrix);

      if(evtType == AliPWG0Helper::kND) {
        fMCAllNDEventMultHist1->Fill(vMCEventMatrix);
      }
      if(evtType != AliPWG0Helper::kSD) {
        fMCAllNSDEventMultHist1->Fill(vMCEventMatrix);
      }
      if(isEventTriggered) {
        fMCTriggerMultHist1->Fill(vMCEventMatrix);
      }
      if(isEventTriggered && isEventOK) {
        fMCEventMultHist1->Fill(vMCEventMatrix);
      }

      //
      // MC histograms for efficiency studies
      //
      Double_t sumPtMC = 0;
      Int_t nPart  = stack->GetNtrack();
      for (Int_t iMc = 0; iMc < nPart; ++iMc) 
      {
        // print MC stack info
        //AlidNdPtHelper::PrintMCInfo(stack,iMc);

        TParticle* particle = stack->Particle(iMc);
        if (!particle)
        continue;

        // only charged particles
        Double_t charge = particle->GetPDG()->Charge()/3.;
        if (TMath::Abs(charge) < 0.001)
          continue;

        // only postive charged 
        if(GetParticleMode() == AlidNdPtHelper::kPlus && charge  < 0.) 
        continue;
      
        // only negative charged 
        if(GetParticleMode() == AlidNdPtHelper::kMinus && charge > 0.) 
        continue;
      
        // physical primary
        Bool_t prim = stack->IsPhysicalPrimary(iMc);
        if(!prim) continue;

        // all primaries in acceptance
        if(!accCuts->AcceptTrack(particle)) continue;

        Double_t gpt = particle->Pt();
        Double_t geta = particle->Eta();

        // sum up pt in the event
	sumPtMC +=gpt; 

        Double_t valMCAllTrackMultHist1[3] = {gpt,geta,multRecTemp};	  
        Double_t valMCAllTrackTrueMultHist1[3] = {gpt,geta,multMCTrueTracks};	  
        Double_t valMCAllTrackTrueMultHist2[3] = {gpt,multRecTemp,multMCTrueTracks};	  

        fMCAllPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
        fMCAllPrimTrackTrueMultHist1->Fill(valMCAllTrackTrueMultHist1);
        fMCAllPrimTrackTrueMultHist2->Fill(valMCAllTrackTrueMultHist2);

        if(evtType == AliPWG0Helper::kND) {
          fMCNDEventAllPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
          fMCNDEventAllPrimTrackTrueMultHist1->Fill(valMCAllTrackTrueMultHist1);
          fMCNDEventAllPrimTrackTrueMultHist2->Fill(valMCAllTrackTrueMultHist2);
        }
        if(evtType != AliPWG0Helper::kSD) {
          fMCNSDEventAllPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
          fMCNSDEventAllPrimTrackTrueMultHist1->Fill(valMCAllTrackTrueMultHist1);
          fMCNSDEventAllPrimTrackTrueMultHist2->Fill(valMCAllTrackTrueMultHist2);
        }
        if(isEventTriggered) { 
	  fMCTriggerPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
          fMCTriggerPrimTrackTrueMultHist1->Fill(valMCAllTrackTrueMultHist1);
          fMCTriggerPrimTrackTrueMultHist2->Fill(valMCAllTrackTrueMultHist2);
	}
        if(isEventTriggered && isEventOK) { 
	  fMCEventPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
          fMCEventPrimTrackTrueMultHist1->Fill(valMCAllTrackTrueMultHist1);
          fMCEventPrimTrackTrueMultHist2->Fill(valMCAllTrackTrueMultHist2);
	}
      }

      //
      // calculate <pt> in the event
      //
      Double_t meanPtMCMult = 0;
      Double_t meanPtMCTrueMult = 0;
      if(multRecTemp) { 
        meanPtMCMult = sumPtMC/multRecTemp; 
      }
      if(multMCTrueTracks) { 
        meanPtMCTrueMult = sumPtMC/multMCTrueTracks; 
      }

      Double_t valMCMeanPtMult[2] = {meanPtMCMult,multRecTemp};	  
      Double_t valMCMeanPtTrueMult[2] = {meanPtMCTrueMult,multMCTrueTracks};	  

      fMCAllPrimTrackMeanPtMult1->Fill(valMCMeanPtMult);
      fMCAllPrimTrackMeanPtTrueMult1->Fill(valMCMeanPtTrueMult);

      if(evtType == AliPWG0Helper::kND) {
          fMCNDEventAllPrimTrackMeanPtMult1->Fill(valMCMeanPtMult);
          fMCNDEventAllPrimTrackMeanPtTrueMult1->Fill(valMCMeanPtTrueMult);
      }
      if(evtType != AliPWG0Helper::kSD) {
          fMCNSDEventAllPrimTrackMeanPtMult1->Fill(valMCMeanPtMult);
          fMCNSDEventAllPrimTrackMeanPtTrueMult1->Fill(valMCMeanPtTrueMult);
      }
      if(isEventTriggered) { 
	  fMCTriggerPrimTrackMeanPtMult1->Fill(valMCMeanPtMult);
          fMCTriggerPrimTrackMeanPtTrueMult1->Fill(valMCMeanPtTrueMult);
      }
      if(isEventTriggered && isEventOK) { 
          fMCEventPrimTrackMeanPtMult1->Fill(valMCMeanPtMult);
          fMCEventPrimTrackMeanPtTrueMult1->Fill(valMCMeanPtTrueMult);
      }
    }
  } // end bUseMC

  if(allChargedTracks) delete allChargedTracks; allChargedTracks = 0;
  if(labelsRec) delete [] labelsRec; labelsRec = 0;

  if(!evtCuts->IsRecVertexRequired() && vtxESD != NULL) delete vtxESD;
}

//_____________________________________________________________________________
void AlidNdPtCorrection::FillHistograms(AlidNdPtHelper::EventObject eventObj, Double_t zv, Int_t multMBTracks) const
{
  //
  // Fill corrected histograms
  //

  Double_t vEventMatrix[2] = {zv,multMBTracks};
  //
  // Correct for efficiency 
  //
  if(eventObj == AlidNdPtHelper::kRecEvents && multMBTracks>0)  
  {
    Double_t corrToMBF = GetCorrFactZvMult(fCorrEventMatrix,zv,multMBTracks);
    Double_t corrToInelF = GetCorrFactZvMult(fCorrTriggerMBtoInelEventMatrix,zv,multMBTracks);
    Double_t corrToNDF = GetCorrFactZvMult(fCorrTriggerMBtoNDEventMatrix,zv,multMBTracks);
    Double_t corrToNSDF = GetCorrFactZvMult(fCorrTriggerMBtoNSDEventMatrix,zv,multMBTracks);
    //printf("corrToMBF %f, corrToInelF %f, corrToNDF %f corrToNSDF %f \n",corrToMBF,corrToInelF,corrToNDF,corrToNSDF);

    fCorrRecEventHist1[0]->Fill(vEventMatrix);
    fCorrRecEventHist1[1]->Fill(vEventMatrix,corrToMBF);
    fCorrRecEventHist1[2]->Fill(vEventMatrix,corrToMBF*corrToInelF);
    fCorrRecEventHist1[3]->Fill(vEventMatrix,corrToMBF*corrToNDF);
    fCorrRecEventHist1[4]->Fill(vEventMatrix,corrToMBF*corrToNSDF);
  }

  if(eventObj==AlidNdPtHelper::kTriggeredEvents && multMBTracks==0) // empty triggered events
  {
    Double_t factLHCBack = 1.;
    if(!IsUseMCInfo()) factLHCBack = fLHCBin0Background->GetBinContent(1); 


    Int_t bin = fZvEmptyEventsNorm->FindBin(zv);
    Double_t factZ = fZvEmptyEventsNorm->GetBinContent(bin);

    Double_t corrToInelF0 = GetCorrFactZvMult(fCorrTriggerMBtoInelEventMatrix,zv,multMBTracks);
    Double_t corrToNDF0 = GetCorrFactZvMult(fCorrTriggerMBtoNDEventMatrix,zv,multMBTracks);
    Double_t corrToNSDF0 = GetCorrFactZvMult(fCorrTriggerMBtoNSDEventMatrix,zv,multMBTracks);
    //printf("factLHCBack %f, factZ %f, corrToInelF0 %f, corrToNDF0 %f, corrToNSDF0 %f \n",factLHCBack,factZ,corrToInelF0,corrToNDF0,corrToNSDF0);

    fCorrRecEventHist2[0]->Fill(vEventMatrix);
    fCorrRecEventHist2[1]->Fill(vEventMatrix,factLHCBack*factZ);
    fCorrRecEventHist2[2]->Fill(vEventMatrix,factLHCBack*factZ*corrToInelF0);
    fCorrRecEventHist2[3]->Fill(vEventMatrix,factLHCBack*factZ*corrToNDF0);
    fCorrRecEventHist2[4]->Fill(vEventMatrix,factLHCBack*factZ*corrToNSDF0);
  }
}

//_____________________________________________________________________________
void AlidNdPtCorrection::FillHistograms(AliESDtrack * const esdTrack, AliStack * const stack, AlidNdPtHelper::TrackObject trackObj, Double_t zv, Int_t mult, Int_t trueMult) const
{
  //
  // Fill ESD track and MC histograms 
  //
  if(!esdTrack) return;

  //Float_t q = esdTrack->Charge();
  Float_t pt = esdTrack->Pt();
  Float_t eta = esdTrack->Eta();
  Float_t phi = esdTrack->Phi();

  if(stack && GetAnalysisMode() == AlidNdPtHelper::kMCRec) 
  {
    Int_t label = TMath::Abs(esdTrack->GetLabel());
   
    TParticle* particle = stack->Particle(label);
    if(!particle) return;
   
    Double_t gq = particle->GetPDG()->Charge()/3.0; // Charge units |e|/3
    if(TMath::Abs(gq)<0.001) return;
    Float_t gpt = particle->Pt();
    Float_t geta = particle->Eta();
    Float_t gphi = particle->Phi();

    // replace reconstructed values with MC
    pt = gpt;
    eta = geta;
    phi = gphi;
  }

  //
  // Fill histograms
  //
  Double_t values[3] = {pt,eta,phi};	  
  fRecTrackHist1[trackObj]->Fill(values);

  //
  // Correct for contamination and efficiency 
  //
  if(trackObj == AlidNdPtHelper::kRecTracks || GetAnalysisMode() == AlidNdPtHelper::kMCRec)  
  {
    // track level corrections
    Double_t trackEffF = 1.0;  
    if(pt < 2.6) trackEffF = GetCorrFactZvPtEta(fCorrTrackMatrix,zv,pt,eta);
    else trackEffF = GetCorrFactZvPtEta(fCorrHighPtTrackMatrix,zv,pt,eta);

    Double_t trackContF = GetContFactZvPtEta(fContTrackMatrix,zv,pt,eta);
    Double_t multTrackContF = GetContFactZvPtEta(fContMultTrackMatrix,zv,pt,eta);
    //printf("zv %f, pt %f, eta %f \n",zv,pt,eta);
    //printf("trackEffF %f, trackContF %f, multTrackContF %f \n", trackEffF, trackContF, multTrackContF);
   
    // track-event level corrections
    Double_t vertexEffF = GetCorrFactZvPtEta(fCorrTrackEventMatrix,zv,pt,eta);
    Double_t trigMBToInel = GetCorrFactZvPtEta(fCorrTriggerMBtoInelTrackEventMatrix,zv,pt,eta);  
    Double_t trigMBToND = GetCorrFactZvPtEta(fCorrTriggerMBtoNDTrackEventMatrix,zv,pt,eta);
    Double_t trigMBToNSD = GetCorrFactZvPtEta(fCorrTriggerMBtoNSDTrackEventMatrix,zv,pt,eta);
    //printf("vertexEffF %f, trigMBToInel %f, trigMBToNSD %f \n", vertexEffF, trigMBToInel, trigMBToNSD);
    
    Double_t corrF[8] = { 1.0, 
                          trackContF,
			  trackContF*trackEffF,
			  trackContF*trackEffF*multTrackContF,
			  trackContF*trackEffF*multTrackContF*vertexEffF,
			  trackContF*trackEffF*multTrackContF*vertexEffF*trigMBToInel,
                          trackContF*trackEffF*multTrackContF*vertexEffF*trigMBToND,
                          trackContF*trackEffF*multTrackContF*vertexEffF*trigMBToNSD
                         }; 
 
    // Fill histograms
    Double_t valCorrRecTrackMultHist1[3] = {pt,eta,mult};	  
    Double_t valCorrRecTrackPt1[1] = {pt};	  
    for(Int_t i=0; i<8; i++) {
      fCorrRecTrackMultHist1[i]->Fill(valCorrRecTrackMultHist1,corrF[i]);
      fCorrRecTrackPt1[i]->Fill(valCorrRecTrackPt1,corrF[i]);

      if( IsUseMCInfo() ) {
        Double_t valCorrRecTrackTrueMultHist1[3] = {pt,eta,trueMult};	  
        Double_t valCorrRecTrackTrueMultHist2[3] = {pt,mult,trueMult};	  

        fCorrRecTrackTrueMultHist1[i]->Fill(valCorrRecTrackTrueMultHist1,corrF[i]);
        fCorrRecTrackTrueMultHist2[i]->Fill(valCorrRecTrackTrueMultHist2,corrF[i]);
      } 
    }
  }
}

void AlidNdPtCorrection::FillHistograms(AliStack * const stack, Int_t /*label*/, AlidNdPtHelper::TrackObject /*trackObj*/, Int_t /*mult*/) const
{
  // Fill MC histograms
  if(!stack) return;

  /*
  TParticle* particle = stack->Particle(label);
  if(!particle) return;

  Int_t mother_pdg = -1;
  TParticle* mother = 0;

  //TParticle* prim_mother = AlidNdPtHelper::FindPrimaryMother(stack,label);
  Int_t motherLabel = particle->GetMother(0); 
  if(motherLabel>0) mother = stack->Particle(motherLabel);
  if(mother) mother_pdg = TMath::Abs(mother->GetPdgCode()); // take abs for visualisation only
  Int_t mech = particle->GetUniqueID(); // production mechanism

  Double_t gq = particle->GetPDG()->Charge()/3.0; // Charge units |e|/3 
  Float_t gpt = particle->Pt();
  Float_t qgpt = particle->Pt() * gq;
  Float_t geta = particle->Eta();
  Float_t gphi = particle->Phi();
  Float_t gpz = particle->Pz();

  Bool_t prim = stack->IsPhysicalPrimary(label);
  Float_t vx = particle->Vx(); Float_t vy = particle->Vy(); Float_t vz = particle->Vz();

  Int_t pid=-1;
  if (TMath::Abs(particle->GetPdgCode()) == kElectron)         { pid = 0; }
    else if (TMath::Abs(particle->GetPdgCode()) == kMuonMinus) { pid = 1; }
    else if (TMath::Abs(particle->GetPdgCode()) == kPiPlus)    { pid = 2; }
    else if (TMath::Abs(particle->GetPdgCode()) == kKPlus)     { pid = 3; }
    else if (TMath::Abs(particle->GetPdgCode()) == kProton)    { pid = 4; }
    else                                                       { pid = 5; }
    */

  //if(!prim) printf("prim_mother %d, mother %d, particle %d, production mech %d\n",prim_mother->GetPdgCode(),mother->GetPdgCode(), particle->GetPdgCode(),mech);
  
}

//_____________________________________________________________________________
Double_t AlidNdPtCorrection::GetCorrFactZvPtEta(THnSparse * const hist, Double_t zv, Double_t pt, Double_t eta) const {
// return correction factor F(zv,pt,eta)

 if(!hist) return 1.;

 //
 TAxis *ax = hist->GetAxis(0);
 TAxis *ay = hist->GetAxis(1);
 TAxis *az = hist->GetAxis(2);

 Int_t binx = ax->FindBin(zv);
 Int_t biny = ay->FindBin(pt);
 Int_t binz = az->FindBin(eta);
 Int_t dim[3] = {binx,biny,binz};

 Double_t fact  = hist->GetBinContent(dim);  

return fact;
}

//_____________________________________________________________________________
Double_t AlidNdPtCorrection::GetContFactZvPtEta(THnSparse * const hist, Double_t zv, Double_t pt, Double_t eta) const {
// return contamination correction factor F(zv,pt,eta)

 if(!hist) return 1.0;

 //
 TAxis *ax = hist->GetAxis(0);
 TAxis *ay = hist->GetAxis(1);
 TAxis *az = hist->GetAxis(2);

 Int_t binx = ax->FindBin(zv);
 Int_t biny = ay->FindBin(pt);
 Int_t binz = az->FindBin(eta);
 Int_t dim[3] = {binx,biny,binz};

 //
 //  additional correction for secondary 
 //  particles with strangeness (data driven)
 // 
 Double_t corrFact = 1.;
 if(!IsUseMCInfo()) corrFact = AlidNdPtHelper::GetStrangenessCorrFactor(pt);
 //printf("pt %f, corrFact %f \n", pt, corrFact);

 Double_t fact  = 1.0 - corrFact*hist->GetBinContent(dim);  
 //Double_t fact  = hist->GetBinContent(dim);  

return fact;
}

//_____________________________________________________________________________
Double_t AlidNdPtCorrection::GetCorrFactZvMult(THnSparse * const hist, Double_t zv, Int_t mult) const {
// return correction factor F(zv,mult)

 if(!hist) return 1.;

 TAxis *ax = hist->GetAxis(0);
 TAxis *ay = hist->GetAxis(1);
 Int_t binx = ax->FindBin(zv);
 Int_t biny = ay->FindBin(mult);
 Int_t dim[2] = {binx,biny};

 Double_t fact  = hist->GetBinContent(dim);  


return fact;
}

//_____________________________________________________________________________
Double_t AlidNdPtCorrection::GetContFactZvMult(THnSparse * const hist, Double_t zv, Int_t mult) const {
// return contamination correction factor F(zv,mult)

 if(!hist) return 1.;

 TAxis *ax = hist->GetAxis(0);
 TAxis *ay = hist->GetAxis(1);
 Int_t binx = ax->FindBin(zv);
 Int_t biny = ay->FindBin(mult);
 Int_t dim[2] = {binx,biny};
 Double_t fact  = 1.0-hist->GetBinContent(dim);  

return fact;
}

//_____________________________________________________________________________
Long64_t AlidNdPtCorrection::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  // collection of generated histograms

  // physics selection
  TList *collPhysSelection = new TList;

  Int_t count=0;
  while((obj = iter->Next()) != 0) {
    AlidNdPtCorrection* entry = dynamic_cast<AlidNdPtCorrection*>(obj);
    if (entry == 0) continue; 
  
    collPhysSelection->Add(entry->GetPhysicsTriggerSelection());

    fEventCount->Add(entry->fEventCount);

    fMCEventHist1->Add(entry->fMCEventHist1);
    fRecEventHist1->Add(entry->fRecEventHist1);
    fRecEventMultHist1->Add(entry->fRecEventMultHist1);

    fMCAllEventMultHist1->Add(entry->fMCAllEventMultHist1);
    fMCAllNDEventMultHist1->Add(entry->fMCAllNDEventMultHist1);
    fMCAllNSDEventMultHist1->Add(entry->fMCAllNSDEventMultHist1);
    fMCTriggerMultHist1->Add(entry->fMCTriggerMultHist1);
    fMCEventMultHist1->Add(entry->fMCEventMultHist1);

    fMCAllPrimTrackMultHist1->Add(entry->fMCAllPrimTrackMultHist1);
    fMCNDEventAllPrimTrackMultHist1->Add(entry->fMCNDEventAllPrimTrackMultHist1);
    fMCNSDEventAllPrimTrackMultHist1->Add(entry->fMCNSDEventAllPrimTrackMultHist1);
    fMCTriggerPrimTrackMultHist1->Add(entry->fMCTriggerPrimTrackMultHist1);
    fMCEventPrimTrackMultHist1->Add(entry->fMCEventPrimTrackMultHist1);

    fMCAllPrimTrackTrueMultHist1->Add(entry->fMCAllPrimTrackTrueMultHist1);
    fMCNDEventAllPrimTrackTrueMultHist1->Add(entry->fMCNDEventAllPrimTrackTrueMultHist1);
    fMCNSDEventAllPrimTrackTrueMultHist1->Add(entry->fMCNSDEventAllPrimTrackTrueMultHist1);
    fMCTriggerPrimTrackTrueMultHist1->Add(entry->fMCTriggerPrimTrackTrueMultHist1);
    fMCEventPrimTrackTrueMultHist1->Add(entry->fMCEventPrimTrackTrueMultHist1);

    fMCAllPrimTrackTrueMultHist2->Add(entry->fMCAllPrimTrackTrueMultHist2);
    fMCNDEventAllPrimTrackTrueMultHist2->Add(entry->fMCNDEventAllPrimTrackTrueMultHist2);
    fMCNSDEventAllPrimTrackTrueMultHist2->Add(entry->fMCNSDEventAllPrimTrackTrueMultHist2);
    fMCTriggerPrimTrackTrueMultHist2->Add(entry->fMCTriggerPrimTrackTrueMultHist2);
    fMCEventPrimTrackTrueMultHist2->Add(entry->fMCEventPrimTrackTrueMultHist2);

    fMCAllPrimTrackMeanPtMult1->Add(entry->fMCAllPrimTrackMeanPtMult1);
    fMCNDEventAllPrimTrackMeanPtMult1->Add(entry->fMCNDEventAllPrimTrackMeanPtMult1);
    fMCNSDEventAllPrimTrackMeanPtMult1->Add(entry->fMCNSDEventAllPrimTrackMeanPtMult1);
    fMCTriggerPrimTrackMeanPtMult1->Add(entry->fMCTriggerPrimTrackMeanPtMult1);
    fMCEventPrimTrackMeanPtMult1->Add(entry->fMCEventPrimTrackMeanPtMult1);

    fMCAllPrimTrackMeanPtTrueMult1->Add(entry->fMCAllPrimTrackMeanPtTrueMult1);
    fMCNDEventAllPrimTrackMeanPtTrueMult1->Add(entry->fMCNDEventAllPrimTrackMeanPtTrueMult1);
    fMCNSDEventAllPrimTrackMeanPtTrueMult1->Add(entry->fMCNSDEventAllPrimTrackMeanPtTrueMult1);
    fMCTriggerPrimTrackMeanPtTrueMult1->Add(entry->fMCTriggerPrimTrackMeanPtTrueMult1);
    fMCEventPrimTrackMeanPtTrueMult1->Add(entry->fMCEventPrimTrackMeanPtTrueMult1);

    fCosmicsHisto->Add(entry->fCosmicsHisto);

    for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) {
      fRecTrackHist1[i]->Add(entry->fRecTrackHist1[i]);
    }

    for(Int_t i=0; i<8; i++) {
      fCorrRecTrackMultHist1[i]->Add(entry->fCorrRecTrackMultHist1[i]);
      fCorrRecTrackTrueMultHist1[i]->Add(entry->fCorrRecTrackTrueMultHist1[i]);
      fCorrRecTrackTrueMultHist2[i]->Add(entry->fCorrRecTrackTrueMultHist2[i]);

      fCorrRecTrackMeanPtMultHist1[i]->Add(entry->fCorrRecTrackMeanPtMultHist1[i]);
      fCorrRecTrackMeanPtTrueMultHist1[i]->Add(entry->fCorrRecTrackMeanPtTrueMultHist1[i]);

      fCorrRecTrackPt1[i]->Add(entry->fCorrRecTrackPt1[i]);
    }

    for(Int_t i=0; i<5; i++) {
      fCorrRecEventHist1[i]->Add(entry->fCorrRecEventHist1[i]);
      fCorrRecEventHist2[i]->Add(entry->fCorrRecEventHist2[i]);
    }

  count++;
  }

  //
  AliPhysicsSelection *trigSelection = GetPhysicsTriggerSelection();
  trigSelection->Merge(collPhysSelection);
  if(collPhysSelection) delete collPhysSelection;

return count;
}
 
//____________________________________________________________________________
Int_t AlidNdPtCorrection::GetTrueMult(THnSparse * const hist, Int_t mult) const
{
//
// get multiplicity of primary particles
//
 if(!hist) return 0;
 Int_t trueMult = 0;

 // 0 bins exluded
 TAxis *ax = hist->GetAxis(0);
 TAxis *ay = hist->GetAxis(1);
 ax->SetRange(1,ax->GetNbins());
 ay->SetRange(1,ay->GetNbins());

 // measured mult
 ax->SetRangeUser((Float_t)mult,(Float_t)mult); 

 // get true multiplicity
 TH1D *h1 = (TH1D *)hist->Projection(1);
 trueMult = (Int_t)h1->GetMean();

 return trueMult;
}

//_____________________________________________________________________________
void AlidNdPtCorrection::Analyse() 
{
  // Analyse histograms
  //
  TH1::AddDirectory(kFALSE);
  TH1 *h = 0, *hs=0, *hsc=0; 
  TH2 *h2D = 0; 

  TObjArray *aFolderObj = new TObjArray;

  //
  // get cuts
  //
  AlidNdPtEventCuts *evtCuts = GetEventCuts(); 
  AlidNdPtAcceptanceCuts *accCuts = GetAcceptanceCuts(); 
  AliESDtrackCuts *esdTrackCuts = GetTrackCuts(); 

  if(!evtCuts || !accCuts || !esdTrackCuts) {
    Error("AlidNdPtCutAnalysis::Analyse()", "cuts not available");
    return;
  }

  //
  // set min and max values
  //
  //Double_t minPt = accCuts->GetMinPt();
  //Double_t maxPt = accCuts->GetMaxPt();
  Double_t minEta = accCuts->GetMinEta();
  Double_t maxEta = accCuts->GetMaxEta()-0.00001;
 
  printf("minEta %f, maxEta %f \n",minEta, maxEta);

  //
  // LHC backgraund in all and 0-bins
  //
  AliPhysicsSelection *trigSel = GetPhysicsTriggerSelection();
  trigSel->SaveHistograms("physics_selection");

  //
  // cosmics background histo
  //
  h2D = fCosmicsHisto->Projection(0,1);
  h2D->SetName("deta_vs_dphi_cosmics");
  aFolderObj->Add(h2D);

  //
  // event level 
  //
  h = fCorrRecEventHist1[0]->Projection(1);
  h->SetName("mult_event_not_corrected");
  aFolderObj->Add(h);

  h = fCorrRecEventHist1[1]->Projection(1);
  h->SetName("mult_event_vertex_corrected");
  aFolderObj->Add(h);

  h = fCorrRecEventHist1[2]->Projection(1);
  h->SetName("mult_trigger_vertex_corrected");
  aFolderObj->Add(h);

  h = fCorrRecEventHist1[3]->Projection(1);
  h->SetName("mult_ND_trigger_vertex_corrected");
  aFolderObj->Add(h);

  h = fCorrRecEventHist1[4]->Projection(1);
  h->SetName("mult_NSD_trigger_vertex_corrected");
  aFolderObj->Add(h);

  // empty events
  h = fCorrRecEventHist2[0]->Projection(1);
  h->SetName("mult_empty_event_not_corrected");
  aFolderObj->Add(h);

  h = fCorrRecEventHist2[1]->Projection(1);
  h->SetName("mult_empty_event_vertex_corrected");
  aFolderObj->Add(h);

  h = fCorrRecEventHist2[2]->Projection(1);
  h->SetName("mult_empty_trigger_vertex_corrected");
  aFolderObj->Add(h);

  h = fCorrRecEventHist2[3]->Projection(1);
  h->SetName("mult_empty_ND_trigger_vertex_corrected");
  aFolderObj->Add(h);

  h = fCorrRecEventHist2[4]->Projection(1);
  h->SetName("mult_empty_NSD_trigger_vertex_corrected");
  aFolderObj->Add(h);
 
  //
  // MC available
  //
  if(IsUseMCInfo()) {

  // mc 
  h = fMCAllEventMultHist1->Projection(1);
  h->SetName("mc_mult_event_acc_prim");
  aFolderObj->Add(h);

  h = fMCAllNDEventMultHist1->Projection(1);
  h->SetName("mc_mult_ND_event_acc_prim");
  aFolderObj->Add(h);

  h = fMCAllNSDEventMultHist1->Projection(1);
  h->SetName("mc_mult_NSD_event_acc_prim");
  aFolderObj->Add(h);

  h = fMCTriggerMultHist1->Projection(1);
  h->SetName("mc_mult_trigger_acc_prim");
  aFolderObj->Add(h);

  h = fMCEventMultHist1->Projection(1);
  h->SetName("mc_mult_trigger_event_acc_prim");
  aFolderObj->Add(h);


  //
  // track level
  //
  
  // limit eta range
  for(Int_t i=0;i<8;i++) { 
      //fCorrRecTrackMultHist1[i]->GetAxis(0)->SetRangeUser(minPt,maxPt);
      //fCorrRecTrackMultHist1[i]->GetAxis(1)->SetRangeUser(minEta,maxEta);
  }
  //fMCAllPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  //fMCAllPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

  //fMCNDEventAllPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  //fMCNDEventAllPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

  //fMCNSDEventAllPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  //fMCNSDEventAllPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

  //fMCTriggerPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  //fMCTriggerPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

  //fMCEventPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  //fMCEventPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

  } // end use MC info 
  
  //
  h2D = fCorrRecTrackMultHist1[3]->Projection(1,0);
  h2D->SetName("pt_eta_rec_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);

  h2D = fCorrRecTrackMultHist1[4]->Projection(1,0);
  h2D->SetName("pt_eta_rec_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);

  h2D = fCorrRecTrackMultHist1[5]->Projection(1,0);
  h2D->SetName("pt_eta_rec_trig_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);

  h2D = fCorrRecTrackMultHist1[6]->Projection(1,0);
  h2D->SetName("pt_eta_rec_ND_trig_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);

  h2D = fCorrRecTrackMultHist1[7]->Projection(1,0);
  h2D->SetName("pt_eta_rec_NSD_trig_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);


  //
  h2D = fCorrRecTrackMultHist1[3]->Projection(2,0);
  h2D->SetName("pt_mult_rec_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);

  h2D = fCorrRecTrackMultHist1[4]->Projection(2,0);
  h2D->SetName("pt_mult_rec_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);

  h2D = fCorrRecTrackMultHist1[5]->Projection(2,0);
  h2D->SetName("pt_mult_rec_trig_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);

  h2D = fCorrRecTrackMultHist1[6]->Projection(2,0);
  h2D->SetName("pt_mult_rec_ND_trig_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);

  h2D = fCorrRecTrackMultHist1[7]->Projection(2,0);
  h2D->SetName("pt_mult_rec_NSD_trig_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h2D);

  // pt axis

  h = fCorrRecTrackMultHist1[0]->Projection(0);
  h->SetName("pt_rec_track_not_corrected");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_track_not_corrected_s");
  aFolderObj->Add(hs);

  //
  h = fCorrRecTrackMultHist1[1]->Projection(0);
  h->SetName("pt_rec_track_cont_corrected");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_track_cont_corrected_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone("pt_rec_track_cont_corr_fact");
  hsc->Divide((TH1D*)aFolderObj->FindObject("pt_rec_track_not_corrected_s"));
  aFolderObj->Add(hsc);

  //
  h = fCorrRecTrackMultHist1[2]->Projection(0);
  h->SetName("pt_rec_track_eff_cont_corrected");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_track_eff_cont_corrected_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone("pt_rec_track_eff_corr_fact");
  hsc->Divide((TH1D*)aFolderObj->FindObject("pt_rec_track_cont_corrected_s"));
  aFolderObj->Add(hsc);

  //
  h = fCorrRecTrackMultHist1[3]->Projection(0);
  h->SetName("pt_rec_track_mult_eff_cont_corrected");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_track_mult_eff_cont_corrected_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone("pt_rec_track_mult_corr_fact");
  hsc->Divide((TH1D*)aFolderObj->FindObject("pt_rec_track_eff_cont_corrected_s"));
  aFolderObj->Add(hsc);

  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hsc = (TH1D*)hs->Clone("pt_rec_track_mult_eff_cont_corr_fact");
  hsc->Divide((TH1D*)aFolderObj->FindObject("pt_rec_track_not_corrected_s"));
  aFolderObj->Add(hsc);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("pt_rec_track_mult_eff_cont_corrected_s_norm");
  hsc->Scale(1./(fCorrRecEventHist1[0]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

  //
  h = fCorrRecTrackMultHist1[4]->Projection(0);
  h->SetName("pt_rec_event_track_mult_eff_cont_corrected");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_event_track_mult_eff_cont_corrected_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("pt_rec_event_track_mult_eff_cont_corrected_s_norm");
  hsc->Scale(1./(fCorrRecEventHist1[1]->Projection(1)->Integral()+fCorrRecEventHist2[1]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

  //
  h = fCorrRecTrackMultHist1[5]->Projection(0);
  h->SetName("pt_rec_trig_event_track_mult_eff_cont_corrected");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_trig_event_track_mult_eff_cont_corrected_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("pt_rec_trig_event_track_mult_eff_cont_corrected_s_norm");
  hsc->Scale(1./(fCorrRecEventHist1[2]->Projection(1)->Integral() + fCorrRecEventHist2[2]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

   // positive eta
  fCorrRecTrackMultHist1[5]->GetAxis(1)->SetRangeUser(0., maxEta);

  h = fCorrRecTrackMultHist1[5]->Projection(0);
  h->SetName("pt_rec_trig_event_track_mult_eff_cont_corrected_posEta");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_trig_event_track_mult_eff_cont_corrected_s_posEta");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("pt_rec_trig_event_track_mult_eff_cont_corrected_s_norm_posEta");
  hsc->Scale(1./(fCorrRecEventHist1[2]->Projection(1)->Integral()+fCorrRecEventHist2[2]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

  // negative eta
  fCorrRecTrackMultHist1[5]->GetAxis(1)->SetRangeUser(minEta, -0.00001);

  h = fCorrRecTrackMultHist1[5]->Projection(0);
  h->SetName("pt_rec_trig_event_track_mult_eff_cont_corrected_negEta");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_trig_event_track_mult_eff_cont_corrected_s_negEta");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("pt_rec_trig_event_track_mult_eff_cont_corrected_s_norm_negEta");
  hsc->Scale(1./(fCorrRecEventHist1[2]->Projection(1)->Integral()+fCorrRecEventHist2[2]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

  fCorrRecTrackMultHist1[5]->GetAxis(1)->SetRange(1, fCorrRecTrackMultHist1[5]->GetAxis(1)->GetNbins());

  //
  h = fCorrRecTrackMultHist1[6]->Projection(0);
  h->SetName("pt_rec_ND_trig_event_track_mult_eff_cont_corrected");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_ND_trig_event_track_mult_eff_cont_corrected_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("pt_rec_ND_trig_event_track_mult_eff_cont_corrected_s_norm");
  hsc->Scale(1./(fCorrRecEventHist1[3]->Projection(1)->Integral()+fCorrRecEventHist2[3]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

  //
  h = fCorrRecTrackMultHist1[7]->Projection(0);
  h->SetName("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected_s_norm");
  hsc->Scale(1./(fCorrRecEventHist1[4]->Projection(1)->Integral() + fCorrRecEventHist2[4]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

  //
  // positive eta
  //
  fCorrRecTrackMultHist1[7]->GetAxis(1)->SetRangeUser(0., maxEta);

  h = fCorrRecTrackMultHist1[7]->Projection(0);
  h->SetName("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected_posEta");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected_s_posEta");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected_s_norm_posEta");
  hsc->Scale(1./(fCorrRecEventHist1[4]->Projection(1)->Integral()+fCorrRecEventHist2[4]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

  //
  // negative eta
  //
  fCorrRecTrackMultHist1[7]->GetAxis(1)->SetRangeUser(minEta, -0.00001);

  h = fCorrRecTrackMultHist1[7]->Projection(0);
  h->SetName("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected_negEta");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected_s_negEta");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected_s_norm_negEta");
  hsc->Scale(1./(fCorrRecEventHist1[4]->Projection(1)->Integral()+fCorrRecEventHist2[4]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

  fCorrRecTrackMultHist1[7]->GetAxis(1)->SetRange(1, fCorrRecTrackMultHist1[7]->GetAxis(1)->GetNbins());

  // eta axis
  h = fCorrRecTrackMultHist1[0]->Projection(1);
  h->SetName("eta_rec_track_not_corrected");
  aFolderObj->Add(h);
  
  h = fCorrRecTrackMultHist1[1]->Projection(1);
  h->SetName("eta_rec_track_cont_corrected");
  aFolderObj->Add(h);

  h = fCorrRecTrackMultHist1[2]->Projection(1);
  h->SetName("eta_rec_track_eff_cont_corrected");
  aFolderObj->Add(h);

  h = fCorrRecTrackMultHist1[3]->Projection(1);
  h->SetName("eta_rec_track_mult_eff_cont_corrected");
  aFolderObj->Add(h);

  h = fCorrRecTrackMultHist1[4]->Projection(1);
  h->SetName("eta_rec_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h);

  h = fCorrRecTrackMultHist1[5]->Projection(1);
  h->SetName("eta_rec_trig_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h);

  h = fCorrRecTrackMultHist1[6]->Projection(1);
  h->SetName("eta_rec_ND_trig_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h);

  h = fCorrRecTrackMultHist1[7]->Projection(1);
  h->SetName("eta_rec_NSD_trig_event_track_mult_eff_cont_corrected");
  aFolderObj->Add(h);


  //
  // MC available
  //
  if(IsUseMCInfo()) {

  //
  h2D = fMCAllPrimTrackMultHist1->Projection(2,0);
  h2D->SetName("mc_all_pt_mult_acc_prim");
  aFolderObj->Add(h2D);

  h2D = fMCAllPrimTrackMultHist1->Projection(1,0);
  h2D->SetName("mc_all_eta_pt_acc_prim");
  aFolderObj->Add(h2D);

  h2D = fMCNDEventAllPrimTrackMultHist1->Projection(2,0);
  h2D->SetName("mc_ND_all_pt_mult_acc_prim");
  aFolderObj->Add(h2D);

  h2D = fMCNDEventAllPrimTrackMultHist1->Projection(1,0);
  h2D->SetName("mc_ND_all_eta_pt_acc_prim");
  aFolderObj->Add(h2D);

  h2D = fMCNSDEventAllPrimTrackMultHist1->Projection(2,0);
  h2D->SetName("mc_NSD_all_pt_mult_acc_prim");
  aFolderObj->Add(h2D);

  h2D = fMCNSDEventAllPrimTrackMultHist1->Projection(1,0);
  h2D->SetName("mc_NSD_all_eta_pt_acc_prim");
  aFolderObj->Add(h2D);

  //

  h = fMCAllPrimTrackMultHist1->Projection(0);
  h->SetName("mc_all_pt_acc_prim");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("mc_all_pt_acc_prim_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("mc_all_pt_acc_prim_s_norm");
  hsc->Scale(1./fMCAllEventMultHist1->Projection(1)->Integral());
  aFolderObj->Add(hsc);

  h = fMCNDEventAllPrimTrackMultHist1->Projection(0);
  h->SetName("mc_ND_all_pt_acc_prim");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("mc_ND_all_pt_acc_prim_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("mc_ND_all_pt_acc_prim_s_norm");
  hsc->Scale(1./fMCAllNDEventMultHist1->Projection(1)->Integral());
  aFolderObj->Add(hsc);

  h = fMCNSDEventAllPrimTrackMultHist1->Projection(0);
  h->SetName("mc_NSD_all_pt_acc_prim");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("mc_NSD_all_pt_acc_prim_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("mc_NSD_all_pt_acc_prim_s_norm");
  hsc->Scale(1./fMCAllNSDEventMultHist1->Projection(1)->Integral());
  aFolderObj->Add(hsc);

  h = fMCTriggerPrimTrackMultHist1->Projection(0);
  h->SetName("mc_trigger_all_pt_acc_prim");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("mc_trigger_all_pt_acc_prim_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("mc_trigger_all_pt_acc_prim_s_norm");
  hsc->Scale(1./fMCTriggerMultHist1->Projection(1)->Integral());
  aFolderObj->Add(hsc);

  h = fMCEventPrimTrackMultHist1->Projection(0);
  h->SetName("mc_all_pt_acc_trig_event_prim");
  hs = AlidNdPtHelper::ScaleByBinWidth(h);
  hs->SetName("mc_all_pt_acc_trig_event_prim_s");
  aFolderObj->Add(hs);

  hsc = (TH1D*)hs->Clone();
  hsc->SetName("mc_all_pt_acc_trig_event_prim_s_norm");
  hsc->Scale(1./fMCEventMultHist1->Projection(1)->Integral());
  aFolderObj->Add(hsc);

  //

  h = fMCAllPrimTrackMultHist1->Projection(1);
  h->SetName("mc_all_eta_acc_prim");
  aFolderObj->Add(h);

  h = fMCNDEventAllPrimTrackMultHist1->Projection(1);
  h->SetName("mc_ND_all_eta_acc_prim");
  aFolderObj->Add(h);

  h = fMCNSDEventAllPrimTrackMultHist1->Projection(1);
  h->SetName("mc_NSD_all_eta_acc_prim");
  aFolderObj->Add(h);

  h = fMCTriggerPrimTrackMultHist1->Projection(1);
  h->SetName("mc_trigger_all_eta_acc_prim");
  aFolderObj->Add(h);

  h = fMCEventPrimTrackMultHist1->Projection(1);
  h->SetName("mc_all_eta_acc_trig_event_prim");
  aFolderObj->Add(h);

  //
  // calculate ratios (rec / mc)
  //
  
  hs = (TH1*)aFolderObj->FindObject("pt_rec_trig_event_track_mult_eff_cont_corrected_s_norm");
  hsc = (TH1D*)hs->Clone("ratio_pt_rec_to_mc_trig_event_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_all_pt_acc_prim_s_norm"));
  aFolderObj->Add(hsc);

  hs = (TH1*)aFolderObj->FindObject("eta_rec_trig_event_track_mult_eff_cont_corrected");
  hsc = (TH1D*)hs->Clone("ratio_eta_rec_to_mc_trig_event_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_all_eta_acc_prim"));
  aFolderObj->Add(hsc);

  //
  hs = (TH1*)aFolderObj->FindObject("pt_rec_ND_trig_event_track_mult_eff_cont_corrected_s_norm");
  hsc = (TH1D*)hs->Clone("ratio_pt_rec_to_mc_ND_trig_event_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_ND_all_pt_acc_prim_s_norm"));
  aFolderObj->Add(hsc);

  hs = (TH1*)aFolderObj->FindObject("eta_rec_ND_trig_event_track_mult_eff_cont_corrected");
  hsc = (TH1D*)hs->Clone("ratio_eta_rec_to_mc_ND_trig_event_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_ND_all_eta_acc_prim"));
  aFolderObj->Add(hsc);

  //
  hs = (TH1*)aFolderObj->FindObject("pt_rec_NSD_trig_event_track_mult_eff_cont_corrected_s_norm");
  hsc = (TH1D*)hs->Clone("ratio_pt_rec_to_mc_NSD_trig_event_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_NSD_all_pt_acc_prim_s_norm"));
  aFolderObj->Add(hsc);

  hs = (TH1*)aFolderObj->FindObject("eta_rec_NSD_trig_event_track_mult_eff_cont_corrected");
  hsc = (TH1D*)hs->Clone("ratio_eta_rec_to_mc_NSD_trig_event_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_NSD_all_eta_acc_prim"));
  aFolderObj->Add(hsc);

  //
  hs = (TH1*)aFolderObj->FindObject("pt_rec_event_track_mult_eff_cont_corrected_s_norm");
  hsc = (TH1D*)hs->Clone("ratio_pt_rec_to_mc_event_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_trigger_all_pt_acc_prim_s_norm"));
  aFolderObj->Add(hsc);

  hs = (TH1*)aFolderObj->FindObject("eta_rec_event_track_mult_eff_cont_corrected");
  hsc = (TH1D*)hs->Clone("ratio_eta_rec_to_mc_event_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_trigger_all_eta_acc_prim"));
  aFolderObj->Add(hsc);

  // track level
  hs = (TH1*)aFolderObj->FindObject("pt_rec_track_mult_eff_cont_corrected_s_norm");
  hsc = (TH1D*)hs->Clone("ratio_pt_rec_to_mc_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_all_pt_acc_trig_event_prim_s_norm"));
  aFolderObj->Add(hsc);

  hs = (TH1*)aFolderObj->FindObject("eta_rec_track_mult_eff_cont_corrected");
  hsc = (TH1D*)hs->Clone("ratio_eta_rec_to_mc_track_mult_eff_cont_corrected"); 
  hsc->Sumw2();
  hsc->Divide((TH1*)aFolderObj->FindObject("mc_all_eta_acc_trig_event_prim"));
  aFolderObj->Add(hsc);

  } // end MC infor available

  // export objects to analysis folder
  fCorrectionFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AlidNdPtCorrection::ExportToFolder(TObjArray * const array) 
{
  // recreate folder avery time and export objects to new one
  //
  AlidNdPtCorrection * comp=this;
  TFolder *folder = comp->GetCorrectionFolder();

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
TFolder* AlidNdPtCorrection::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
