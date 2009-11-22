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

#include <iostream>

#include "TFile.h"
#include "TCint.h"
#include "TH1.h"
#include "TH2.h"

#include "AliHeader.h"  
#include "AliGenEventHeader.h"  
#include "AliStack.h"  
#include "AliESDEvent.h"  
#include "AliMCEvent.h"  
#include "AliESDtrackCuts.h"  
#include "AliLog.h" 

#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"
#include "AlidNdPtCorrection.h"

using namespace std;

ClassImp(AlidNdPtCorrection)

//_____________________________________________________________________________
//AlidNdPtCorrection::AlidNdPtCorrection(): TNamed(),
  AlidNdPtCorrection::AlidNdPtCorrection(): AlidNdPt(),
  fCorrectionFolder(0),
  fMCEventHist1(0),
  fRecEventHist1(0),
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
  fEventMultCorrelationMatrix(0),
  fZvNorm(0),
  fZvEmptyEventsNorm(0),
  fCorrTriggerMBtoInelEventMatrix(0),
  fCorrTriggerMBtoNDEventMatrix(0),
  fCorrTriggerMBtoNSDEventMatrix(0),
  fCorrEventMatrix(0),
  fCorrTriggerMBtoInelTrackEventMatrix(0),
  fCorrTriggerMBtoNDTrackEventMatrix(0),
  fCorrTriggerMBtoNSDTrackEventMatrix(0),
  fCorrTrackEventMatrix(0),
  fCorrTrackMatrix(0),
  fContTrackMatrix(0),
  fContMultTrackMatrix(0),
  fCorrMatrixFileName("")
{
  // default constructor
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    fRecTrackHist1[i]=0;     
  }

  for(Int_t i=0; i<8; i++) { 
    fCorrRecTrackMultHist1[i] = 0;
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
  fEventMultCorrelationMatrix(0),
  fZvNorm(0),
  fZvEmptyEventsNorm(0),
  fCorrTriggerMBtoInelEventMatrix(0),
  fCorrTriggerMBtoNDEventMatrix(0),
  fCorrTriggerMBtoNSDEventMatrix(0),
  fCorrEventMatrix(0),
  fCorrTriggerMBtoInelTrackEventMatrix(0),
  fCorrTriggerMBtoNDTrackEventMatrix(0),
  fCorrTriggerMBtoNSDTrackEventMatrix(0),
  fCorrTrackEventMatrix(0),
  fCorrTrackMatrix(0),
  fContTrackMatrix(0),
  fContMultTrackMatrix(0),
  fCorrMatrixFileName(corrMatrixFileName)
{
  // constructor
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    fRecTrackHist1[i]=0;     
  }

  for(Int_t i=0; i<8; i++) { 
    fCorrRecTrackMultHist1[i] = 0;
    fPtvsPt[i] = 0;
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

  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) { 
    if(fRecTrackHist1[i]) delete fRecTrackHist1[i]; fRecTrackHist1[i]=0;
  }

  for(Int_t i=0; i<8; i++) { 
    if(fCorrRecTrackMultHist1[i]) delete fCorrRecTrackMultHist1[i]; fCorrRecTrackMultHist1[i]=0;
    if(fPtvsPt[i]) delete fPtvsPt[i]; fPtvsPt[i]=0;
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

  const Int_t ptNbins = 56; 
  Double_t binsPt[ptNbins+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};
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
  Double_t maxMCNSDEventAllPrimTrackMultHist1[3]={20.,1.,140.9}; 
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

  Int_t binsCorrRecTrackMultHist1[3]={ptNbins,etaNbins,150};
  Double_t minCorrRecTrackMultHist1[3]={0.,-1.,-0.5}; 
  Double_t maxCorrRecTrackMultHist1[3]={20.,1.,149.5};

  Int_t binsPtvsPt[3]={ptNbins,320};
  Double_t minPtvsPt[3]={0.,0.}; 
  Double_t maxPtvsPt[3]={20.,16.};

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

    sprintf(name,"fPtvsPt_%d",i);
    sprintf(title,"Pt:Pt");
    fPtvsPt[i] = new THnSparseF(name,title,2,binsPtvsPt,minPtvsPt,maxPtvsPt);
    fPtvsPt[i]->SetBinEdges(0,binsPt);
    fPtvsPt[i]->GetAxis(0)->SetTitle("Pt (GeV/c)"); 
    fPtvsPt[i]->GetAxis(1)->SetTitle("Pt (GeV/c)"); 
    fPtvsPt[i]->Sumw2();
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
      }

      //
      // event level corrections (zv,mult_MB)
      //
 
      // trigger bias correction (MBtoND) 
      fCorrTriggerMBtoNDEventMatrix = (THnSparseF*)folder->FindObject("zv_mult_trig_MBtoND_corr_matrix");
      if(!fCorrTriggerMBtoNDEventMatrix) {
         Printf("No %s matrix \n", "zv_mult_trig_MBtoND_corr_matrix");
      }

      // trigger bias correction (MBtoNSD)
      fCorrTriggerMBtoNSDEventMatrix = (THnSparseF*)folder->FindObject("zv_mult_trig_MBtoNSD_corr_matrix");
      if(!fCorrTriggerMBtoNSDEventMatrix) {
         Printf("No %s matrix \n", "zv_mult_trig_MBtoNSD_corr_matrix");
      }

      // trigger bias correction (MBtoInel)
      fCorrTriggerMBtoInelEventMatrix = (THnSparseF*)folder->FindObject("zv_mult_trig_MBtoInel_corr_matrix");
      if(!fCorrTriggerMBtoInelEventMatrix) {
         Printf("No %s matrix \n", "zv_mult_trig_MBtoInel_corr_matrix"); 
      }
     
      // vertex reconstruction efficiency correction
      fCorrEventMatrix = (THnSparseF*)folder->FindObject("zv_mult_event_corr_matrix");
      if(!fCorrEventMatrix) {
         Printf("No %s matrix \n", "zv_mult_event_corr_matrix");
      }

      //
      // histogram needed for empty events corrections
      //
      fZvNorm = (TH1D*)folder->FindObject("zv_distribution_norm");
      if(!fZvNorm) {
         Printf("No %s matrix \n", "fZvNorm");
      }

      fZvEmptyEventsNorm = (TH1D*)folder->FindObject("zv_empty_events_norm");
      if(!fZvEmptyEventsNorm) {
         Printf("No %s matrix \n", "fZvEmptyEventsNorm");
      }

      //
      // track-event level corrections (zv,pt,eta)
      //

      // trigger bias correction (MBtoND) 
      fCorrTriggerMBtoNDTrackEventMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_trig_MBtoND_corr_matrix");
      if(!fCorrTriggerMBtoNDTrackEventMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_trig_MBtoND_corr_matrix");
      }

      // trigger bias correction (MBtoNSD)
      fCorrTriggerMBtoNSDTrackEventMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_trig_MBtoNSD_corr_matrix");
      if(!fCorrTriggerMBtoNSDTrackEventMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_trig_MBtoNSD_corr_matrix");
      }

      // trigger bias correction (MBtoInel) 
      fCorrTriggerMBtoInelTrackEventMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_trig_MBtoInel_corr_matrix");
      if(!fCorrTriggerMBtoInelTrackEventMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_trig_MBtoInel_corr_matrix");
      }
    
      // vertex reconstruction efficiency correction (zv,pt,eta)
      fCorrTrackEventMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_event_corr_matrix");
      if(!fCorrTrackEventMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_event_corr_matrix");
      }

      // track reconstruction efficiency correction (zv,pt,eta)
      fCorrTrackMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_corr_matrix");
      if(!fCorrTrackMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_corr_matrix");
      }

      // secondary tracks contamination correction (zv,pt,eta)
      fContTrackMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_track_cont_matrix");
      if(!fContTrackMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_track_cont_matrix");
      }

      // multiply reconstructed tracks correction
      fContMultTrackMatrix = (THnSparseF*)folder->FindObject("zv_pt_eta_mult_track_cont_matrix");
      if(!fContMultTrackMatrix) {
         Printf("No %s matrix \n", "zv_pt_eta_mult_track_cont_matrix");
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
  if(evtCuts->IsTriggerRequired())  {
    static AliTriggerAnalysis* triggerAnalysis = new AliTriggerAnalysis;
    isEventTriggered = triggerAnalysis->IsTriggerFired(esdEvent, GetTrigger());
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
    // in Zvtx, pt and eta ranges
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
    isRecVertex = AlidNdPtHelper::TestRecVertex(vtxESD, GetAnalysisMode(), kFALSE); // should be moved to AcceptEvent
  }
  if( IsUseMCInfo() &&  !evtCuts->IsRecVertexRequired() ) {
    vtxESD = new AliESDVertex(vtxMC[2],10.,genHeader->NProduced());
    isRecVertex = kTRUE;
  }
  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD) && isRecVertex; 
  //printf("isEventOK %d \n",isEventOK);

  //
  // get multiplicity of min. bias tracks
  //
  Int_t multMBTracks = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC || GetAnalysisMode() == AlidNdPtHelper::kMCRec) {  
    multMBTracks = AlidNdPtHelper::GetTPCMBTrackMult(esdEvent,evtCuts,accCuts,esdTrackCuts);
  } 
  else if(GetAnalysisMode() == AlidNdPtHelper::kTPCSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kMCRec) {
    multMBTracks = AlidNdPtHelper::GetSPDMBTrackMult(esdEvent,0.0);
  } 
  else {
    AliDebug(AliLog::kError, Form("Found analysis type %s", GetAnalysisMode()));
    return; 
  }

  //
  // correct event and track histograms
  //
  TObjArray *allChargedTracks=0;
  Int_t multRec=0, multRecTemp=0;
  Int_t *labelsRec=0;

  if(isEventOK && isEventTriggered)
  {
    // get all charged tracks
    //allChargedTracks = AlidNdPtHelper::GetAllChargedTracks(esdEvent,vtxESD,GetAnalysisMode());
    allChargedTracks = AlidNdPtHelper::GetAllChargedTracks(esdEvent,GetAnalysisMode());
    if(!allChargedTracks) return;

    Int_t entries = allChargedTracks->GetEntries();
    labelsRec = new Int_t[entries];

    // calculate mult of reconstructed tracks
    for(Int_t i=0; i<entries;++i) 
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(i);
      if(!track) continue;
      if(accCuts->AcceptTrack(track)) 
      {
        if(esdTrackCuts->AcceptTrack(track)) 
	{
	  multRecTemp++;
        }
      }
    }  

    for(Int_t i=0; i<entries;++i) 
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(i);
      if(!track) continue;
        
      if(accCuts->AcceptTrack(track)) 
      {
        if(esdTrackCuts->AcceptTrack(track)) 
	{
	  // track-level corrections
          if(GetAnalysisMode()==AlidNdPtHelper::kMCRec && IsUseMCInfo()) { 
            FillHistograms(track,stack,AlidNdPtHelper::kRecTracks,vtxMC[2],multRecTemp); 
	  } else {
            FillHistograms(track,stack,AlidNdPtHelper::kRecTracks,vtxESD->GetZv(),multRecTemp); 
	  }
	  labelsRec[multRec] = TMath::Abs(track->GetLabel());
	  multRec++;
        }
      }
    }
    // event-level corrections
    if(GetAnalysisMode()==AlidNdPtHelper::kMCRec && IsUseMCInfo()) { 
      FillHistograms(AlidNdPtHelper::kRecEvents,vtxMC[2],multMBTracks);
    }
    else {
      FillHistograms(AlidNdPtHelper::kRecEvents,vtxESD->GetZv(),multMBTracks);
    }

    // control event histograms
    Double_t vRecEventHist1[3] = {vtxESD->GetXv(),vtxESD->GetYv(),vtxESD->GetZv()};
    fRecEventHist1->Fill(vRecEventHist1);
  } 

  // empty events corrections
  // no reconstructed zv
  if(isEventTriggered && multMBTracks==0) {
    if(GetAnalysisMode()==AlidNdPtHelper::kMCRec && IsUseMCInfo()) {
      FillHistograms(AlidNdPtHelper::kTriggeredEvents,vtxMC[2],multMBTracks);
    }
    else {
      Double_t zv = fZvNorm->GetRandom();
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
      if(isEventTriggered) fMCTriggerMultHist1->Fill(vMCEventMatrix);
      if(isEventTriggered && isEventOK) fMCEventMultHist1->Fill(vMCEventMatrix);

      //
      // MC histograms for efficiency studies
      //
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
        if (charge == 0.0)
          continue;
      
        // physical primary
        Bool_t prim = stack->IsPhysicalPrimary(iMc);

        // all primaries in acceptance
        if(!accCuts->AcceptTrack(particle)) continue;
        if(!prim) continue;

        Double_t gpt = particle->Pt();
        Double_t geta = particle->Eta();

        Double_t valMCAllTrackMultHist1[3] = {gpt,geta,multRec};	  
        fMCAllPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
        if(evtType == AliPWG0Helper::kND) {
          fMCNDEventAllPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
        }
        if(evtType != AliPWG0Helper::kSD) {
          fMCNSDEventAllPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
        }
        if(isEventTriggered) fMCTriggerPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
        if(isEventTriggered && isEventOK) fMCEventPrimTrackMultHist1->Fill(valMCAllTrackMultHist1);
      }
    }
  } // end bUseMC

  if(allChargedTracks) delete allChargedTracks; allChargedTracks = 0;
  if(labelsRec) delete [] labelsRec; labelsRec = 0;

  if(!evtCuts->IsRecVertexRequired() && vtxESD != NULL) 
  delete vtxESD;
}

//_____________________________________________________________________________
void AlidNdPtCorrection::FillHistograms(AlidNdPtHelper::EventObject eventObj, Double_t zv, Int_t multMBTracks)
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
    Int_t bin = fZvEmptyEventsNorm->FindBin(zv);
    Double_t Fz = fZvEmptyEventsNorm->GetBinContent(bin);
    Double_t corrToInelF0 = GetCorrFactZvMult(fCorrTriggerMBtoInelEventMatrix,zv,multMBTracks);
    Double_t corrToNDF0 = GetCorrFactZvMult(fCorrTriggerMBtoNDEventMatrix,zv,multMBTracks);
    Double_t corrToNSDF0 = GetCorrFactZvMult(fCorrTriggerMBtoNSDEventMatrix,zv,multMBTracks);
    //printf("Fz %f, corrToInelF0 %f, corrToNDF0 %f, corrToNSDF0 %f \n",Fz,corrToInelF0,corrToNDF0,corrToNSDF0);

    fCorrRecEventHist2[0]->Fill(vEventMatrix);
    fCorrRecEventHist2[1]->Fill(vEventMatrix,Fz);
    fCorrRecEventHist2[2]->Fill(vEventMatrix,Fz*corrToInelF0);
    fCorrRecEventHist2[3]->Fill(vEventMatrix,Fz*corrToNDF0);
    fCorrRecEventHist2[4]->Fill(vEventMatrix,Fz*corrToNSDF0);
  }
}

//_____________________________________________________________________________
void AlidNdPtCorrection::FillHistograms(AliESDtrack *esdTrack, AliStack *stack, AlidNdPtHelper::TrackObject trackObj, Double_t zv, Int_t mult)
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
    if(gq==0) return;
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
    Double_t trackEffF = GetCorrFactZvPtEta(fCorrTrackMatrix,zv,pt,eta);
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
    Double_t valPtvsPt[2] = {pt,pt};	  
    for(Int_t i=0; i<8; i++) {
      fCorrRecTrackMultHist1[i]->Fill(valCorrRecTrackMultHist1,corrF[i]);
      fPtvsPt[i]->Fill(valPtvsPt,corrF[i]);
    }
  }
}

//_____________________________________________________________________________
void AlidNdPtCorrection::FillHistograms(AliStack *stack, Int_t /*label*/, AlidNdPtHelper::TrackObject /*trackObj*/, Int_t /*mult*/)
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
Double_t AlidNdPtCorrection::GetCorrFactZvPtEta(THnSparse *hist, Double_t zv, Double_t pt, Double_t eta) const {
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
Double_t AlidNdPtCorrection::GetContFactZvPtEta(THnSparse *hist, Double_t zv, Double_t pt, Double_t eta) const {
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

 Double_t fact  = 1.0-hist->GetBinContent(dim);  
 //Double_t fact  = hist->GetBinContent(dim);  

return fact;
}

//_____________________________________________________________________________
Double_t AlidNdPtCorrection::GetCorrFactZvMult(THnSparse *hist, Double_t zv, Int_t mult) const {
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
Double_t AlidNdPtCorrection::GetContFactZvMult(THnSparse *hist, Double_t zv, Int_t mult) const {
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
Long64_t AlidNdPtCorrection::Merge(TCollection* list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  // collection of generated histograms

  Int_t count=0;
  while((obj = iter->Next()) != 0) {
    AlidNdPtCorrection* entry = dynamic_cast<AlidNdPtCorrection*>(obj);
    if (entry == 0) continue; 
  
    fMCEventHist1->Add(entry->fMCEventHist1);
    fRecEventHist1->Add(entry->fRecEventHist1);

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

    for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) {
      fRecTrackHist1[i]->Add(entry->fRecTrackHist1[i]);
    }

    for(Int_t i=0; i<8; i++) {
      fCorrRecTrackMultHist1[i]->Add(entry->fCorrRecTrackMultHist1[i]);
      fPtvsPt[i]->Add(entry->fPtvsPt[i]);
    }

    for(Int_t i=0; i<5; i++) {
      fCorrRecEventHist1[i]->Add(entry->fCorrRecEventHist1[i]);
      fCorrRecEventHist2[i]->Add(entry->fCorrRecEventHist2[i]);
    }

  count++;
  }

return count;
}
 
Int_t AlidNdPtCorrection::GetTrueMult(THnSparse *hist, Int_t mult)
{
 if(!hist) return 0;
 Int_t true_mult = 0;

 // 0 bins exluded
 TAxis *ax = hist->GetAxis(0);
 TAxis *ay = hist->GetAxis(1);
 ax->SetRange(1,ax->GetNbins());
 ay->SetRange(1,ay->GetNbins());

 // measured mult
 ax->SetRangeUser((Float_t)mult,(Float_t)mult); 

 // get true multiplicity
 TH1D *h1 = (TH1D *)hist->Projection(1);
 true_mult = (Int_t)h1->GetMean();

 return true_mult;
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
  Double_t minPt = accCuts->GetMinPt();
  Double_t maxPt = accCuts->GetMaxPt();
  Double_t minEta = accCuts->GetMinEta();
  Double_t maxEta = accCuts->GetMaxEta()-0.00001;

  //
  // pt profile
  //
  char name[256];
  for(Int_t i=0; i<8; i++) {
    h2D = fPtvsPt[i]->Projection(1,0);
    sprintf(name,"PtvsMeanPt_%d",i);
    aFolderObj->Add(h2D);
  }

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
      fCorrRecTrackMultHist1[i]->GetAxis(0)->SetRangeUser(minPt,maxPt);
      fCorrRecTrackMultHist1[i]->GetAxis(1)->SetRangeUser(minEta,maxEta);
  }
  fMCAllPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  fMCAllPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

  fMCNDEventAllPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  fMCNDEventAllPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

  fMCNSDEventAllPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  fMCNSDEventAllPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

  fMCTriggerPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  fMCTriggerPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

  fMCEventPrimTrackMultHist1->GetAxis(0)->SetRangeUser(minPt,maxPt);
  fMCEventPrimTrackMultHist1->GetAxis(1)->SetRangeUser(minEta,maxEta);

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
  hsc->Scale(1./(fCorrRecEventHist1[4]->Projection(1)->Integral()+fCorrRecEventHist2[4]->Projection(1)->Integral()));
  aFolderObj->Add(hsc);

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
TFolder* AlidNdPtCorrection::ExportToFolder(TObjArray * array) 
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
