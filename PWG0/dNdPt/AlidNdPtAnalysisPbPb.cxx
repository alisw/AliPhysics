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
// AlidNdPtAnalysisPbPb class. 
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
//------------------------------------------------------------------------------

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "THnSparse.h"

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
#include "AliTriggerAnalysis.h"

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"
#include "AlidNdPtAnalysisPbPb.h"


using namespace std;

ClassImp(AlidNdPtAnalysisPbPb)

//_____________________________________________________________________________
  AlidNdPtAnalysisPbPb::AlidNdPtAnalysisPbPb(): AlidNdPt(),
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
  fRecTrackHist2(0)
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
  Init();
}

//_____________________________________________________________________________
AlidNdPtAnalysisPbPb::AlidNdPtAnalysisPbPb(Char_t* name, Char_t* title): AlidNdPt(name,title),
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
  fRecTrackHist2(0)
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

  Init();
}

//_____________________________________________________________________________
AlidNdPtAnalysisPbPb::~AlidNdPtAnalysisPbPb() {
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
    if(fRecTrackMultHist1[i]) delete fRecTrackMultHist1[i]; fRecTrackMultHist1[i]=0;
  }
  if(fRecMCTrackHist1) delete fRecMCTrackHist1; fRecMCTrackHist1=0;
  if(fMCMultRecTrackHist1) delete fMCMultRecTrackHist1; fMCMultRecTrackHist1=0; 
  if(fRecTrackHist2) delete fRecTrackHist2; fRecTrackHist2=0; 
  //
  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb::Init(){
  //
  // Init histograms
  //

  const Int_t multNbins = 47;
//  const Int_t ptNbinsTrackEventCorr = 38;
  const Int_t ptNbins = 55;
  const Int_t etaNbins = 30;
  const Int_t zvNbins = 12;

//  Double_t binsMult[multNbins+1] = {-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,
//                                     9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,
// 				     19.5,20.5, 21.5, 22.5, 23.5, 24.5, 29.5, 149.5};//for pp

  Double_t binsMult[multNbins+1] = {-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,
                                     9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,
 				     19.5,20.5, 30.5, 40.5 , 50.5 , 60.5 , 70.5 , 80.5 , 90.5 , 100.5 , 
					200.5, 300.5 , 400.5 , 500.5 ,600.5, 700.5 , 800.5 , 900.5 , 1000.5 ,
					2000.5 , 3000.5 , 4000.5 , 5000.5 , 6000.5 , 7000.5 , 8000.5 , 9000.5 , 10000.5};//forPbPb



//   Double_t binsPtTrackEventCorr[ptNbinsTrackEventCorr+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,3.0,4.0,6.0,10.0,16.0};

  Double_t binsPt[ptNbins+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};

  Double_t binsEta[etaNbins+1] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};

  Double_t binsZv[zvNbins+1] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};

  Int_t binsTrackEventCorrMatrix[3]={zvNbins,ptNbins,etaNbins};

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
  Double_t minEventMatrix[2]={-30.,-0.5}; 
  Double_t maxEventMatrix[2]={30.,10000.5 }; 

  fGenEventMatrix = new THnSparseF("fGenEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fGenEventMatrix->SetBinEdges(0,binsZv);
  fGenEventMatrix->SetBinEdges(1,binsMult);
  fGenEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fGenEventMatrix->Sumw2();
  //
  fTriggerEventMatrix = new THnSparseF("fTriggerEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fTriggerEventMatrix->SetBinEdges(0,binsZv);
  fTriggerEventMatrix->SetBinEdges(1,binsMult);
  fTriggerEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fTriggerEventMatrix->Sumw2();
  //
  fRecEventMatrix = new THnSparseF("fRecEventMatrix","mcZv:multMB",2,binsEventMatrix,minEventMatrix,maxEventMatrix);
  fRecEventMatrix->SetBinEdges(0,binsZv);
  fRecEventMatrix->SetBinEdges(1,binsMult);
  fRecEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecEventMatrix->GetAxis(1)->SetTitle("multiplicity MB");
  fRecEventMatrix->Sumw2();

  // 
  // track to event corrections
  //
  fGenTrackEventMatrix = new THnSparseF("fGenTrackEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenTrackEventMatrix->SetBinEdges(0,binsZv);
  fGenTrackEventMatrix->SetBinEdges(1,binsPt);
  fGenTrackEventMatrix->SetBinEdges(2,binsEta);
  fGenTrackEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackEventMatrix->Sumw2();
  //
  fTriggerTrackEventMatrix = new THnSparseF("fTriggerTrackEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fTriggerTrackEventMatrix->SetBinEdges(0,binsZv);
  fTriggerTrackEventMatrix->SetBinEdges(1,binsPt);
  fTriggerTrackEventMatrix->SetBinEdges(2,binsEta);
  fTriggerTrackEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fTriggerTrackEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fTriggerTrackEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fTriggerTrackEventMatrix->Sumw2();
  //
  fRecTrackEventMatrix = new THnSparseF("fRecTrackEventMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecTrackEventMatrix->SetBinEdges(0,binsZv);
  fRecTrackEventMatrix->SetBinEdges(1,binsPt);
  fRecTrackEventMatrix->SetBinEdges(2,binsEta);
  fRecTrackEventMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackEventMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecTrackEventMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecTrackEventMatrix->Sumw2();

  //
  // tracks correction matrices
  //
  fGenTrackMatrix = new THnSparseF("fGenTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenTrackMatrix->SetBinEdges(0,binsZv);
  fGenTrackMatrix->SetBinEdges(1,binsPt);
  fGenTrackMatrix->SetBinEdges(2,binsEta);
  fGenTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenTrackMatrix->Sumw2();

  fGenPrimTrackMatrix = new THnSparseF("fGenPrimTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fGenPrimTrackMatrix->SetBinEdges(0,binsZv);
  fGenPrimTrackMatrix->SetBinEdges(1,binsPt);
  fGenPrimTrackMatrix->SetBinEdges(2,binsEta);
  fGenPrimTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fGenPrimTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fGenPrimTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fGenPrimTrackMatrix->Sumw2();


  fRecPrimTrackMatrix = new THnSparseF("fRecPrimTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecPrimTrackMatrix->SetBinEdges(0,binsZv);
  fRecPrimTrackMatrix->SetBinEdges(1,binsPt);
  fRecPrimTrackMatrix->SetBinEdges(2,binsEta);
  fRecPrimTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecPrimTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecPrimTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecPrimTrackMatrix->Sumw2();

  //
  fRecTrackMatrix = new THnSparseF("fRecTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecTrackMatrix->SetBinEdges(0,binsZv);
  fRecTrackMatrix->SetBinEdges(1,binsPt);
  fRecTrackMatrix->SetBinEdges(2,binsEta);
  fRecTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecTrackMatrix->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fRecTrackMatrix->GetAxis(2)->SetTitle("Eta");
  fRecTrackMatrix->Sumw2();

  fRecSecTrackMatrix = new THnSparseF("fRecSecTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecSecTrackMatrix->SetBinEdges(0,binsZv);
  fRecSecTrackMatrix->SetBinEdges(1,binsPt);
  fRecSecTrackMatrix->SetBinEdges(2,binsEta);
  fRecSecTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecSecTrackMatrix->GetAxis(1)->SetTitle("Pt (GeV/c)");
  fRecSecTrackMatrix->GetAxis(2)->SetTitle("Eta");
  fRecSecTrackMatrix->Sumw2();

  //
  fRecMultTrackMatrix = new THnSparseF("fRecMultTrackMatrix","mcZv:mcPt:mcEta",3,binsTrackEventCorrMatrix);
  fRecMultTrackMatrix->SetBinEdges(0,binsZv);
  fRecMultTrackMatrix->SetBinEdges(1,binsPt);
  fRecMultTrackMatrix->SetBinEdges(2,binsEta);
  fRecMultTrackMatrix->GetAxis(0)->SetTitle("mcZv (cm)");
  fRecMultTrackMatrix->GetAxis(1)->SetTitle("mcPt (GeV/c)");
  fRecMultTrackMatrix->GetAxis(2)->SetTitle("mcEta");
  fRecMultTrackMatrix->Sumw2();

  //
  // Control analysis histograms
  //
  Int_t binsMCEventHist1[3]={100,100,zvNbins};
  Double_t minMCEventHist1[3]={-0.1,-0.1,-30.}; 
  Double_t maxMCEventHist1[3]={0.1,0.1,30.}; 
  fMCEventHist1 = new THnSparseF("fMCEventHist1","mcXv:mcYv:mcZv",3,binsMCEventHist1,minMCEventHist1,maxMCEventHist1);
  fMCEventHist1->SetBinEdges(2,binsZv);
  fMCEventHist1->GetAxis(0)->SetTitle("mcXv (cm)");
  fMCEventHist1->GetAxis(1)->SetTitle("mcYv (cm)");
  fMCEventHist1->GetAxis(2)->SetTitle("mcZv (cm)");
  fMCEventHist1->Sumw2();

  //
  Int_t binsRecEventHist1[3]={100,100,zvNbins};
  Double_t minRecEventHist1[3]={-3.,-3.,-30.}; 
  Double_t maxRecEventHist1[3]={3.,3.,30.}; 
  
  fRecEventHist1 = new THnSparseF("fRecEventHist1","Xv:Yv:Zv",3,binsRecEventHist1,minRecEventHist1,maxRecEventHist1);
  fRecEventHist1->SetBinEdges(2,binsZv);
  fRecEventHist1->GetAxis(0)->SetTitle("Xv (cm)");
  fRecEventHist1->GetAxis(1)->SetTitle("Yv (cm)");
  fRecEventHist1->GetAxis(2)->SetTitle("Zv (cm)");
  fRecEventHist1->Sumw2();

  //
  Int_t binsRecEventHist2[3]={zvNbins,multNbins,multNbins};
  Double_t minRecEventHist2[3]={-30.,-0.5,-0.5}; 
  Double_t maxRecEventHist2[3]={30.,10000.5,10000.5}; 
  
  fRecEventHist2 = new THnSparseF("fRecEventHist2","Zv:multMB:mult",3,binsRecEventHist2,minRecEventHist2,maxRecEventHist2);
  fRecEventHist2->SetBinEdges(0,binsZv);
  fRecEventHist2->SetBinEdges(1,binsMult);
  fRecEventHist2->SetBinEdges(2,binsMult);
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
  fRecMCEventHist1->SetBinEdges(2,binsZv);
  fRecMCEventHist1->GetAxis(0)->SetTitle("Xv-mcXv (cm)");
  fRecMCEventHist1->GetAxis(1)->SetTitle("Yv-mcYv (cm)");
  fRecMCEventHist1->GetAxis(2)->SetTitle("Zv-mcZv (cm)");
  fRecMCEventHist1->Sumw2();

  //
  Int_t binsRecMCEventHist2[3]={100,100,multNbins};
  Double_t minRecMCEventHist2[3]={-10.0*kFact,-10.0*kFact,-0.5}; 
  Double_t maxRecMCEventHist2[3]={10.0*kFact,10.0*kFact,10000.5}; 

  fRecMCEventHist2 = new THnSparseF("fRecMCEventHist2","Xv-mcXv:Zv-mcZv:mult",3,binsRecMCEventHist2,minRecMCEventHist2,maxRecMCEventHist2);
  fRecMCEventHist2->SetBinEdges(2,binsMult);
  fRecMCEventHist2->GetAxis(0)->SetTitle("Xv-mcXv (cm)");
  fRecMCEventHist2->GetAxis(1)->SetTitle("Zv-mcZv (cm)");
  fRecMCEventHist2->GetAxis(2)->SetTitle("multiplicity");
  fRecMCEventHist2->Sumw2();
  //
  char name[256];
  char title[256];
  for(Int_t i=0; i<AlidNdPtHelper::kCutSteps; i++) 
  {
  // THnSparse track histograms
 
  Int_t binsMCTrackHist1[3]=  {ptNbins, etaNbins, 90};
  Double_t minMCTrackHist1[3]={0.,-1.5,0.}; 
  Double_t maxMCTrackHist1[3]={16.,1.5,2.*TMath::Pi()}; 
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
  Double_t minMCPrimTrackHist1[5]={0.,-1.5,0.,0.,0.}; 
  Double_t maxMCPrimTrackHist1[5]={16.,1.5,6.,20.,4000.}; 
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

  Int_t binsMCPrimTrackHist2[5]=  {4000,20,4000};
  Double_t minMCPrimTrackHist2[5]={0.,0.,0.}; 
  Double_t maxMCPrimTrackHist2[5]={4000.,20.,4000.}; 
  sprintf(name,"fMCPrimTrackHist2_%d",i);
  sprintf(title,"pdg:mech:mother");
  
  fMCPrimTrackHist2[i] = new THnSparseF(name,title,5,binsMCPrimTrackHist2,minMCPrimTrackHist2,maxMCPrimTrackHist2);
  fMCPrimTrackHist2[i]->GetAxis(0)->SetTitle("pdg");
  fMCPrimTrackHist2[i]->GetAxis(1)->SetTitle("mech");
  fMCPrimTrackHist2[i]->GetAxis(2)->SetTitle("mother");
  fMCPrimTrackHist2[i]->Sumw2();

  Int_t binsMCSecTrackHist1[5]=  {ptNbins,etaNbins,6,20,4000};
  Double_t minMCSecTrackHist1[5]={0.,-1.5,0.,0.,0.}; 
  Double_t maxMCSecTrackHist1[5]={16.,1.5,6.,20.,4000.}; 
  sprintf(name,"fMCSecTrackHist1_%d",i);
  sprintf(title,"mcPt:mcEta:mcPhi:pid:mech:mother");
  
  fMCSecTrackHist1[i] = new THnSparseF(name,title,5,binsMCSecTrackHist1,minMCSecTrackHist1,maxMCSecTrackHist1);
  fMCSecTrackHist1[i]->SetBinEdges(0,binsPt);
  fMCSecTrackHist1[i]->SetBinEdges(1,binsEta);
  fMCSecTrackHist1[i]->GetAxis(0)->SetTitle("mcPt (GeV/c)");
  fMCSecTrackHist1[i]->GetAxis(1)->SetTitle("mcEta");
  fMCSecTrackHist1[i]->GetAxis(2)->SetTitle("pid");
  fMCSecTrackHist1[i]->GetAxis(3)->SetTitle("mech");
  fMCSecTrackHist1[i]->GetAxis(4)->SetTitle("mother");
  fMCSecTrackHist1[i]->Sumw2();

  Int_t binsRecTrackHist1[3]={ptNbins,etaNbins,90};
  Double_t minRecTrackHist1[3]={0.,-1.5,0.}; 
  Double_t maxRecTrackHist1[3]={16.,1.5,2.*TMath::Pi()};
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
  Int_t binsRecTrackMultHist1[2]={ptNbins,multNbins};
  Double_t minRecTrackMultHist1[2]={0.,-0.5}; 
  Double_t maxRecTrackMultHist1[2]={16.,10000.5};
  sprintf(name,"fRecTrackMultHist_%d",i);
  sprintf(title,"Pt:Mult");
  fRecTrackMultHist1[i] = new THnSparseF(name,title,2,binsRecTrackMultHist1,minRecTrackMultHist1,maxRecTrackMultHist1);
  fRecTrackMultHist1[i]->SetBinEdges(0,binsPt);
  fRecTrackMultHist1[i]->SetBinEdges(1,binsMult);
  fRecTrackMultHist1[i]->GetAxis(0)->SetTitle("Pt (GeV/c)");
  fRecTrackMultHist1[i]->GetAxis(1)->SetTitle("multiplicity");
  fRecTrackMultHist1[i]->Sumw2();
  }

  Int_t binsRecMCTrackHist1[4] = {ptNbins,etaNbins,100,100};
  Double_t minRecMCTrackHist1[4]={0.,-1.5,-0.5,-0.5}; 
  Double_t maxRecMCTrackHist1[4]={16.,1.5,0.5,0.5}; 

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
  Double_t minMCMultRecTrackHist1[3]={0.,-1.5,0.}; 
  Double_t maxMCMultRecTrackHist1[3]={16.,1.5,6.}; 
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
  Double_t minRecTrackHist2[5]={0., 0., 0., -1.5, 0.};
  Double_t maxRecRecTrackHist2[5]={160.,10., 16., 1.5, 2.*TMath::Pi()};

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
}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb::Process(AliESDEvent *const esdEvent, AliMCEvent *const mcEvent)
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
  AliTriggerAnalysis *trigAna = NULL;

  if(evtCuts->IsTriggerRequired())  
  {
    //
    trigSel = GetPhysicsTriggerSelection();
    if(!trigSel) {
      printf("cannot get trigSel \n");
      return;
    }

    //
    if(IsUseMCInfo()) 
    { 
      trigSel->SetAnalyzeMC();

      
        isEventTriggered = trigSel->IsCollisionCandidate(esdEvent);
	
        if(GetTrigger() == AliTriggerAnalysis::kV0AND) 
	{
          trigAna = trigSel->GetTriggerAnalysis();
          if(!trigAna) 
            return;

          isEventTriggered = trigAna->IsOfflineTriggerFired(esdEvent, GetTrigger());
        }//if(GetTrigger() == AliTriggerAnalysis::kV0AND)
     }//if(IsUseMCInfo())
  }//if(evtCuts->IsTriggerRequired())



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

    Double_t vMCEventHist1[3]={vtxMC[0],vtxMC[1],vtxMC[2]};
    fMCEventHist1->Fill(vMCEventHist1);

    // multipliticy of all MC primary tracks
    // in Zv, pt and eta ranges)
    multMCTrueTracks = AlidNdPtHelper::GetMCTrueTrackMult(mcEvent,evtCuts,accCuts);

  } // end bUseMC

  // get reconstructed vertex  
  const AliESDVertex* vtxESD = 0; 
  Bool_t isRecVertex = kTRUE;
  if(evtCuts->IsRecVertexRequired()) 
  {
//    Bool_t bRedoTPCVertex = evtCuts->IsRedoTPCVertex();
    //Bool_t bUseConstraints = evtCuts->IsUseBeamSpotConstraint();
    //vtxESD = AlidNdPtHelper::GetVertex(esdEvent,evtCuts,accCuts,esdTrackCuts,GetAnalysisMode(),kFALSE,bRedoTPCVertex,bUseConstraints); 
    //isRecVertex = AlidNdPtHelper::TestRecVertex(vtxESD, esdEvent->GetPrimaryVertexSPD(), GetAnalysisMode(), kFALSE);
    if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
      vtxESD = esdEvent->GetPrimaryVertexTPC();
      //isRecVertex = AlidNdPtHelper::TestRecVertex(vtxESD, esdEvent->GetPrimaryVertexSPD(), GetAnalysisMode(), kFALSE);
    }
  }

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD) && isRecVertex; 
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK, isEventTriggered);
  //printf("GetAnalysisMode() %d \n",GetAnalysisMode());

  // vertex contributors
  Int_t multMBTracks = 0; 
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) 
  {
     if(vtxESD->GetStatus() && isRecVertex) {
         //multMBTracks = AlidNdPtHelper::GetTPCMBTrackMult(esdEvent,evtCuts,accCuts,esdTrackCuts);
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



      // only postive charged 
      if(GetParticleMode() == AlidNdPtHelper::kPlus && track->Charge() < 0) 
        continue;
      
      // only negative charged 
      if(GetParticleMode() == AlidNdPtHelper::kMinus && track->Charge() > 0) 
        continue;

      //
      FillHistograms(track,stack,AlidNdPtHelper::kAllTracks); 
      labelsAll[multAll] = TMath::Abs(track->GetLabel());

      multAll++;

      if(esdTrackCuts->AcceptTrack(track) && accCuts->AcceptTrack(track)) {

         FillHistograms(track,stack,AlidNdPtHelper::kRecTracks); 
         labelsRec[multRec] = TMath::Abs(track->GetLabel());

         multRec++;

      }
     }//loop over entries


     // fill track multiplicity histograms
     //FillHistograms(allChargedTracks,labelsAll,multAll,labelsAcc,multAcc,labelsRec,multRec);

     Double_t vRecEventHist1[3] = {vtxESD->GetXv(),vtxESD->GetYv(),vtxESD->GetZv()};
     fRecEventHist1->Fill(vRecEventHist1);

     Double_t vRecEventHist2[3] = {vtxESD->GetZv(),multMBTracks,multRec};
     fRecEventHist2->Fill(vRecEventHist2);

   } // triggered and event vertex

   if(IsUseMCInfo())  
   {

     //Double_t vMultTrueEventMatrix[2] = { multRec, multMCTrueTracks };

     // 
     // event level corrections (zv,N_MB)
     //
     // all inelastic
     Double_t vEventMatrix[2] = {vtxMC[2],multMBTracks};
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

       // checked accepted
       if(accCuts->AcceptTrack(particle)) 
       {
        Double_t vTrackEventMatrix[3] = {vtxMC[2], particle->Pt(), particle->Eta()}; 
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
         Double_t vRecMCEventHist1[3] = {vtxESD->GetXv()-vtxMC[0],vtxESD->GetYv()-vtxMC[1],vtxESD->GetZv()-vtxMC[2]};
         fRecMCEventHist1->Fill(vRecMCEventHist1);//

         Double_t vRecMCEventHist2[3] = {vtxESD->GetXv()-vtxMC[0],vtxESD->GetZv()-vtxMC[2],multMBTracks};
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

         Double_t vTrackMatrix[3] = {vtxMC[2],particle->Pt(),particle->Eta()}; 

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
	       }
               if(!prim) fRecSecTrackMatrix->Fill(vTrackMatrix);

	       // fill control histograms
               if(fHistogramsOn) 
                 FillHistograms(stack,iMc,AlidNdPtHelper::kRecTracks); 
	      
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

  if(!evtCuts->IsRecVertexRequired() && vtxESD != NULL) delete vtxESD;

}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb::FillHistograms(TObjArray *const allChargedTracks,Int_t *const labelsAll,Int_t multAll,Int_t *const labelsAcc,Int_t multAcc,Int_t *const labelsRec,Int_t multRec) {
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
       }//out
     }
  }
}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb::FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, AlidNdPtHelper::TrackObject trackObj)
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
  Double_t values[3] = {pt,eta,phi};	  
  fRecTrackHist1[trackObj]->Fill(values);

  //
  // Fill rec vs MC information
  //
  if(!stack) return;

  Int_t label = TMath::Abs(esdTrack->GetLabel()); 
  //if(label == 0) return;

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
    Double_t vTrackPtCorrelationMatrix[3]={pt,gpt,geta};
    fTrackPtCorrelationMatrix->Fill(vTrackPtCorrelationMatrix);

    Double_t vRecMCTrackHist1[4]={gpt,geta,dpt,deta};
    fRecMCTrackHist1->Fill(vRecMCTrackHist1);
  }
}

//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb::FillHistograms(AliStack *const stack, Int_t label, AlidNdPtHelper::TrackObject trackObj)
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
Long64_t AlidNdPtAnalysisPbPb::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  //
  TList *collPhysSelection = new TList;

  // collection of generated histograms

  Int_t count=0;
  while((obj = iter->Next()) != 0) {
    AlidNdPtAnalysisPbPb* entry = dynamic_cast<AlidNdPtAnalysisPbPb*>(obj);
    if (entry == 0) continue; 

    // physics selection
    //printf("entry->GetPhysicsTriggerSelection() %p \n", entry->GetPhysicsTriggerSelection());
    collPhysSelection->Add(entry->GetPhysicsTriggerSelection());
    
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
      fRecTrackMultHist1[i]->Add(entry->fRecTrackMultHist1[i]);
    }
    fRecMCTrackHist1->Add(entry->fRecMCTrackHist1);
    fMCMultRecTrackHist1->Add(entry->fMCMultRecTrackHist1);
    fRecTrackHist2->Add(entry->fRecTrackHist2);

  count++;
  }

  AliPhysicsSelection *trigSelection = GetPhysicsTriggerSelection();
  trigSelection->Merge(collPhysSelection);
  if(collPhysSelection) delete collPhysSelection;

return count;
}



//_____________________________________________________________________________
void AlidNdPtAnalysisPbPb::Analyse() 
{
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
  AliPhysicsSelection *trigSel = GetPhysicsTriggerSelection();
  trigSel->SaveHistograms("physics_selection");

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
TFolder* AlidNdPtAnalysisPbPb::ExportToFolder(TObjArray * const array) 
{
  // recreate folder avery time and export objects to new one
  //
  AlidNdPtAnalysisPbPb * comp=this;
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
TFolder* AlidNdPtAnalysisPbPb::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
