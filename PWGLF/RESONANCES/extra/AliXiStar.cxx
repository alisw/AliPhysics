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

////////////////////////////////////////////////////////////////////////////////
//
//  This class is used to reconstruct the neutral Xi(1530) resonance.  
//  This class essentially combines charged Xi candidates from the Xi Vertexer 
//  with primary charged pions. 
//
//  authors: Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////



#include <iostream>
#include <math.h>
#include "TChain.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjString.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"


#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODcascade.h"
#include "AliESDcascade.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"

#include "AliXiStar.h"

#define PI 3.1415927


// Author: Dhevan Gangadharan

ClassImp(AliXiStar)

//________________________________________________________________________
AliXiStar::AliXiStar():
AliAnalysisTaskSE(),
  fname(0),
  fAOD(0x0), 
  fESD(0x0), 
  fOutputList(0x0),
  fTrackCut(0x0),
  fPIDResponse(0x0),
  fEC(0x0),
  fEvt(0x0),
  fTempStruct(0x0),
  fZvertexBins(0),
  fEventsToMix(0),
  fMultBins(0),
  fMultLimits(),
  fMCcase(0),
  fAODcase(0),
  fEventCounter(0),
  fSigmaCutProton(0),
  fSigmaCutPionFirst(0),
  fSigmaCutPionSecond(0),
  fSigmaCutPionThird(0),
  fDCAVtxProton(0),
  fDCAVtxPionFirst(0),
  fDCAVtxPionSecond(0),
  fDCAVtxLambda(0),
  fDCAProtonPion(0),
  fDCALambdaPion(0),
  fLambdaDecayLengthXY(0),
  fXiDecayLengthXY(0),
  fMaxDecayLength(0),
  fLamCosTheta(0),
  fXiCosTheta(0),
  fXiStarCosTheta(0),
  fMassWindow(0),
  fCovMatrix(),
  fTrueMassPr(0), 
  fTrueMassPi(0), 
  fTrueMassK(0), 
  fTrueMassLam(0), 
  fTrueMassXi(0),
  fESDTrack4(0x0), 
  fXiTrack(0x0),
  fCutList(0)
{
}
//________________________________________________________________________
AliXiStar::AliXiStar(const char *name, Bool_t AODdecision, Bool_t MCdecision, Int_t CutListOption) 
  : AliAnalysisTaskSE(name), 
    fname(name),
    fAOD(0x0), 
    fESD(0x0), 
    fOutputList(0x0),
    fTrackCut(0x0),
    fPIDResponse(0x0),
    fEC(0x0),
    fEvt(0x0),
    fTempStruct(0x0),
    fZvertexBins(0),
    fEventsToMix(0),
    fMultBins(0),
    fMultLimits(),
    fMCcase(MCdecision),
    fAODcase(AODdecision),
    fEventCounter(0),
    fSigmaCutProton(0),
    fSigmaCutPionFirst(0),
    fSigmaCutPionSecond(0),
    fSigmaCutPionThird(0),
    fDCAVtxProton(0),
    fDCAVtxPionFirst(0),
    fDCAVtxPionSecond(0),
    fDCAVtxLambda(0),
    fDCAProtonPion(0),
    fDCALambdaPion(0),
    fLambdaDecayLengthXY(0),
    fXiDecayLengthXY(0),
    fMaxDecayLength(0),
    fLamCosTheta(0),
    fXiCosTheta(0),
    fXiStarCosTheta(0),
    fMassWindow(0),
    fCovMatrix(),
    fTrueMassPr(0), 
    fTrueMassPi(0), 
    fTrueMassK(0), 
    fTrueMassLam(0), 
    fTrueMassXi(0),
    fESDTrack4(0x0), 
    fXiTrack(0x0),
    fCutList(CutListOption)
{
  // Main Constructor
  
  // Define output slots here 
  // Output slot #1
  DefineOutput(1, TList::Class());
  
}
//________________________________________________________________________
AliXiStar::AliXiStar(const AliXiStar &obj) 
  : AliAnalysisTaskSE(obj.fname),
    fname(obj.fname),
    fAOD(obj.fAOD), 
    fESD(obj.fESD), 
    fOutputList(obj.fOutputList),
    fTrackCut(obj.fTrackCut),
    fPIDResponse(obj.fPIDResponse),
    fEC(obj.fEC),
    fEvt(obj.fEvt),
    fTempStruct(obj.fTempStruct),
    fZvertexBins(obj.fZvertexBins),
    fEventsToMix(obj.fEventsToMix),
    fMultBins(obj.fMultBins),
    fMultLimits(),
    fMCcase(obj.fMCcase),
    fAODcase(obj.fAODcase),
    fEventCounter(obj.fEventCounter),
    fSigmaCutProton(obj.fSigmaCutProton),
    fSigmaCutPionFirst(obj.fSigmaCutPionFirst),
    fSigmaCutPionSecond(obj.fSigmaCutPionSecond),
    fSigmaCutPionThird(obj.fSigmaCutPionThird),
    fDCAVtxProton(obj.fDCAVtxProton),
    fDCAVtxPionFirst(obj.fDCAVtxPionFirst),
    fDCAVtxPionSecond(obj.fDCAVtxPionSecond),
    fDCAVtxLambda(obj.fDCAVtxLambda),
    fDCAProtonPion(obj.fDCAProtonPion),
    fDCALambdaPion(obj.fDCALambdaPion),
    fLambdaDecayLengthXY(obj.fLambdaDecayLengthXY),
    fXiDecayLengthXY(obj.fXiDecayLengthXY),
    fMaxDecayLength(obj.fMaxDecayLength),
    fLamCosTheta(obj.fLamCosTheta),
    fXiCosTheta(obj.fXiCosTheta),
    fXiStarCosTheta(obj.fXiStarCosTheta),
    fMassWindow(obj.fMassWindow),
    fCovMatrix(),
    fTrueMassPr(obj.fTrueMassPr), 
    fTrueMassPi(obj.fTrueMassPi), 
    fTrueMassK(obj.fTrueMassK), 
    fTrueMassLam(obj.fTrueMassLam), 
    fTrueMassXi(obj.fTrueMassXi),
    fESDTrack4(obj.fESDTrack4), 
    fXiTrack(obj.fXiTrack),
    fCutList(obj.fCutList)
    
{
  // Copy constructor  
}
//________________________________________________________________________
AliXiStar &AliXiStar::operator=(const AliXiStar &obj) 
{
  // Assignment operator  
  if (this == &obj)
    return *this;
  
  fAOD = obj.fAOD;
  fESD = obj.fESD; 
  fOutputList = obj.fOutputList;
  fTrackCut = obj.fTrackCut;
  fPIDResponse = obj.fPIDResponse;
  fEC = obj.fEC;
  fEvt = obj.fEvt;
  fTempStruct = obj.fTempStruct;
  fZvertexBins = obj.fZvertexBins;
  fEventsToMix = obj.fEventsToMix;
  fMultBins = obj.fMultBins;
  fMultLimits[0] = obj.fMultLimits[0]; fMultLimits[1] = obj.fMultLimits[1]; fMultLimits[2] = obj.fMultLimits[2];
  fMultLimits[3] = obj.fMultLimits[3]; fMultLimits[4] = obj.fMultLimits[4]; fMultLimits[5] = obj.fMultLimits[5];
  fMultLimits[6] = obj.fMultLimits[6]; fMultLimits[7] = obj.fMultLimits[7]; fMultLimits[8] = obj.fMultLimits[8];
  fMultLimits[9] = obj.fMultLimits[9]; fMultLimits[10] = obj.fMultLimits[10]; fMultLimits[11] = obj.fMultLimits[11];
  fMCcase = obj.fMCcase;
  fAODcase = obj.fAODcase;
  fEventCounter = obj.fEventCounter;
  fSigmaCutProton = obj.fSigmaCutProton;
  fSigmaCutPionFirst = obj.fSigmaCutPionFirst;
  fSigmaCutPionSecond = obj.fSigmaCutPionSecond;
  fSigmaCutPionThird = obj.fSigmaCutPionThird;
  fDCAVtxProton = obj.fDCAVtxProton;
  fDCAVtxPionFirst = obj.fDCAVtxPionFirst;
  fDCAVtxPionSecond = obj.fDCAVtxPionSecond;
  fDCAVtxLambda = obj.fDCAVtxLambda;
  fDCAProtonPion = obj.fDCAProtonPion;
  fDCALambdaPion = obj.fDCALambdaPion;
  fLambdaDecayLengthXY = obj.fLambdaDecayLengthXY;
  fXiDecayLengthXY = obj.fXiDecayLengthXY;
  fMaxDecayLength = obj.fMaxDecayLength;
  fLamCosTheta = obj.fLamCosTheta;
  fXiCosTheta = obj.fXiCosTheta;
  fXiStarCosTheta = obj.fXiStarCosTheta;
  fMassWindow = obj.fMassWindow;
  fCovMatrix[0] = obj.fCovMatrix[0]; fCovMatrix[1] = obj.fCovMatrix[1]; fCovMatrix[2] = obj.fCovMatrix[2];
  fCovMatrix[3] = obj.fCovMatrix[3]; fCovMatrix[4] = obj.fCovMatrix[4]; fCovMatrix[5] = obj.fCovMatrix[5];
  fCovMatrix[6] = obj.fCovMatrix[6]; fCovMatrix[7] = obj.fCovMatrix[7]; fCovMatrix[8] = obj.fCovMatrix[8];
  fCovMatrix[9] = obj.fCovMatrix[9]; fCovMatrix[10] = obj.fCovMatrix[10]; fCovMatrix[11] = obj.fCovMatrix[11];
  fCovMatrix[12] = obj.fCovMatrix[12]; fCovMatrix[13] = obj.fCovMatrix[13]; fCovMatrix[14] = obj.fCovMatrix[14];
  fCovMatrix[15] = obj.fCovMatrix[15]; fCovMatrix[16] = obj.fCovMatrix[16]; fCovMatrix[17] = obj.fCovMatrix[17];
  fCovMatrix[18] = obj.fCovMatrix[18]; fCovMatrix[19] = obj.fCovMatrix[19]; fCovMatrix[20] = obj.fCovMatrix[20];
  fTrueMassPr = obj.fTrueMassPr; 
  fTrueMassPi = obj.fTrueMassPi; 
  fTrueMassK = obj.fTrueMassK; 
  fTrueMassLam = obj.fTrueMassLam; 
  fTrueMassXi = obj.fTrueMassXi;
  fESDTrack4 = obj.fESDTrack4; 
  fXiTrack = obj.fXiTrack; 
  fCutList = obj.fCutList;
    
  return (*this);
}
//________________________________________________________________________
AliXiStar::~AliXiStar()
{
  // Destructor
  if(fname) delete fname;
  if(fAOD) delete fAOD; 
  if(fESD) delete fESD; 
  if(fOutputList) delete fOutputList;
  if(fTrackCut) delete fTrackCut;
  if(fPIDResponse) delete fPIDResponse;
  if(fEC) delete fEC;
  if(fEvt) delete fEvt;
  if(fTempStruct) delete fTempStruct;
  if(fESDTrack4) delete fESDTrack4; 
  if(fXiTrack) delete fXiTrack; 
}
//________________________________________________________________________
void AliXiStar::XiStarInit()
{
  //
  //Inits cuts and analysis settings
  //
  
  fEventCounter=0;// event counter initialization
  cout<<"AliXiStar XiStarInit() call"<<endl;
   
  
  ///////////////////////////////////////////////
  // Track Cuts for ESD analysis
  fTrackCut = new AliESDtrackCuts();
  fTrackCut->SetPtRange(.15,1000);
  fTrackCut->SetAcceptKinkDaughters(kFALSE);
  fTrackCut->SetMinNClustersTPC(70);
  fTrackCut->SetRequireTPCRefit(kTRUE);
  ////////////////////////////////////////////////
  
  fZvertexBins = 20;
  fMultBins = 11;// This must also be set in AliXiStar.h
  if(fMCcase) fEventsToMix = 0;
  else fEventsToMix = 40;

  // multiplicity edges for event mixing bins
  fMultLimits[0]=0, fMultLimits[1]=5, fMultLimits[2]=10, fMultLimits[3]=15, fMultLimits[4]=20, fMultLimits[5]=25;
  fMultLimits[6]=30, fMultLimits[7]=35, fMultLimits[8]=40, fMultLimits[9]=45, fMultLimits[10]=50, fMultLimits[11]=150;
  
  
  fEC = new AliXiStarEventCollection **[fZvertexBins];
  for(unsigned short i=0; i<fZvertexBins; i++){
    
    fEC[i] = new AliXiStarEventCollection *[fMultBins];
    
    for(unsigned short j=0; j<fMultBins; j++){
      
      fEC[i][j] = new AliXiStarEventCollection(fEventsToMix+1);
    }
  }
  
  
  fTempStruct = new AliXiStarTrackStruct[kNbinsM];

  fESDTrack4 = new AliESDtrack();
  fXiTrack = new AliESDtrack();
   

  ///////////////////////////////////////////////////////
  // Reconstruction Cuts
  //
 
  
  // Xi cuts are from the 7 TeV Multi-strangeness paper values 
  if(fCutList == 0){// mean value cuts
    fSigmaCutProton = 1000;// min for protons with P > .7 GeV/c
    fSigmaCutPionFirst = 1000;// min
    fSigmaCutPionSecond = 1000;// min
    fSigmaCutPionThird =1000;// min
    fDCAVtxProton = 0.04;// min, 
    fDCAVtxPionFirst = 0.04;// min
    fDCAVtxPionSecond = 0.05;// min
    fDCAVtxLambda = 0.07;// min
    fDCAProtonPion = 1.6;// max
    fDCALambdaPion = 1.6;// max
    fLambdaDecayLengthXY = 1.4;// min
    fXiDecayLengthXY = 0.8;// min
    fMaxDecayLength = 100.;// max for Lambdas and Cascades
    fLamCosTheta = 0.97;// min; used for Xi reconstruction
    fXiCosTheta = 0.97;// min; used for Xi(1530) reconstruction
    fXiStarCosTheta = -1;// min;
    fMassWindow = .006;// window half-width
  }
  
  // Systematics Part 1: PV dca variation (stricter)
  if(fCutList == 1){
    fSigmaCutProton = 1000;// min for protons with P > .7 GeV/c
    fSigmaCutPionFirst = 1000;// min
    fSigmaCutPionSecond = 1000;// min
    fSigmaCutPionThird =1000;// min
    fDCAVtxProton = 0.2;// min, 
    fDCAVtxPionFirst = 0.5;// min
    fDCAVtxPionSecond = 0.5;// min
    fDCAVtxLambda = 0.2;// min
    fDCAProtonPion = 1.6;// max
    fDCALambdaPion = 1.6;// max
    fLambdaDecayLengthXY = 1.4;// min
    fXiDecayLengthXY = 0.8;// min
    fMaxDecayLength = 100.;// max for Lambdas and Cascades
    fLamCosTheta = 0.97;// min; used for Xi reconstruction
    fXiCosTheta = 0.97;// min; used for Xi(1530) reconstruction
    fXiStarCosTheta = -1;// min;
    fMassWindow = .006;// window half-width
  }
  
  // Systematics Part 2: Topological cut variation (looser)
  if(fCutList == 2){
    fSigmaCutProton = 1000;// min for protons with P > .7 GeV/c
    fSigmaCutPionFirst = 1000;// min
    fSigmaCutPionSecond = 1000;// min
    fSigmaCutPionThird =1000;// min
    fDCAVtxProton = 0.04;// min, 
    fDCAVtxPionFirst = 0.04;// min
    fDCAVtxPionSecond = 0.05;// min
    fDCAVtxLambda = 0.07;// min
    fDCAProtonPion = 2.4;// max
    fDCALambdaPion = 2.4;// max
    fLambdaDecayLengthXY = 0.7;// min
    fXiDecayLengthXY = 0.4;// min
    fMaxDecayLength = 100.;// max for Lambdas and Cascades
    fLamCosTheta = 0.9;// min; used for Xi reconstruction
    fXiCosTheta = 0.9;// min; used for Xi(1530) reconstruction
    fXiStarCosTheta = -1;// min;
    fMassWindow = .006;// window half-width
  }

  // Systematics Part 3: Topological cut variation (stricter)
  if(fCutList == 3){
    fSigmaCutProton = 1000;// min for protons with P > .7 GeV/c
    fSigmaCutPionFirst = 1000;// min
    fSigmaCutPionSecond = 1000;// min
    fSigmaCutPionThird =1000;// min
    fDCAVtxProton = 0.04;// min, 
    fDCAVtxPionFirst = 0.04;// min
    fDCAVtxPionSecond = 0.05;// min
    fDCAVtxLambda = 0.07;// min
    fDCAProtonPion = 1.0;// max
    fDCALambdaPion = 1.0;// max
    fLambdaDecayLengthXY = 1.8;// min
    fXiDecayLengthXY = 1.0;// min
    fMaxDecayLength = 100.;// max for Lambdas and Cascades
    fLamCosTheta = 0.99;// min; used for Xi reconstruction
    fXiCosTheta = 0.99;// min; used for Xi(1530) reconstruction
    fXiStarCosTheta = -1;// min;
    fMassWindow = .006;// window half-width
  }

  // PDG mass values
  fTrueMassPr=.93827, fTrueMassPi=.13957, fTrueMassK=.493677, fTrueMassLam=1.11568, fTrueMassXi=1.32171;
  
  // The following CovMatrix is set so that PropogateToDCA() ignores track errors. Only used to propagate Xi to third pion for XiStar reconstruction 
  for(Int_t i=0; i<21; i++) fCovMatrix[i]=0;
  fCovMatrix[0]=1, fCovMatrix[2]=1, fCovMatrix[5]=1, fCovMatrix[9]=1, fCovMatrix[14]=1, fCovMatrix[20]=1;
  
  
}
//________________________________________________________________________
void AliXiStar::UserCreateOutputObjects()
{
   
  XiStarInit();// Initialize settings
  
  // Create histograms
  fOutputList = new TList();
  fOutputList->SetOwner();
  
  TH3F *fVertexDist1 = new TH3F("fVertexDist1","Vertex Distribution",20,-1,1, 20,-1,1, 600,-30,30);
  fVertexDist1->GetXaxis()->SetTitle("X Vertex (cm)");
  fVertexDist1->GetYaxis()->SetTitle("Y Vertex (cm)");
  fVertexDist1->GetZaxis()->SetTitle("Z Vertex (cm)");
  fOutputList->Add(fVertexDist1);
  
  TH3F *fVertexDist3 = new TH3F("fVertexDist3","Vertex Distribution",20,-1,1, 20,-1,1, 600,-30,30);
  fVertexDist3->GetXaxis()->SetTitle("X Vertex (cm)");
  fVertexDist3->GetYaxis()->SetTitle("Y Vertex (cm)");
  fVertexDist3->GetZaxis()->SetTitle("Z Vertex (cm)");
  fOutputList->Add(fVertexDist3);

  TH2F *fDCADist = new TH2F("fDCADist","DCA distribution",kNbinsM,-.5,kNbinsM-.5, 100,0,10);
  fOutputList->Add(fDCADist);
  
  
  TH3F *fMultDist3d = new TH3F("fMultDist3d","Multiplicity Distribution",kNbinsM,-.5,kNbinsM-.5, kNbinsM,-.5,kNbinsM-.5, kNbinsM,-.5,kNbinsM-.5);
  fMultDist3d->GetXaxis()->SetTitle("Multiplicity");
  fMultDist3d->GetYaxis()->SetTitle("Positive Multiplicity");
  fMultDist3d->GetZaxis()->SetTitle("Negative Multiplicity");
  fMultDist3d->SetMarkerStyle(kFullCircle);
  fOutputList->Add(fMultDist3d);
  
 
  TH1F *fMultDist1 = new TH1F("fMultDist1","Multiplicity Distribution",kNbinsM,-.5,kNbinsM-.5);
  fMultDist1->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist1);
  
  TH1F *fMultDist2 = new TH1F("fMultDist2","Multiplicity Distribution",kNbinsM,-.5,kNbinsM-.5);
  fMultDist2->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist2);
  
  TH1F *fMultDist3 = new TH1F("fMultDist3","Multiplicity Distribution",kNbinsM,-.5,kNbinsM-.5);
  fMultDist3->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist3);

  TH1F *fMultDist4 = new TH1F("fMultDist4","Multiplicity Distribution",kNbinsM,-.5,kNbinsM-.5);
  fMultDist4->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist4);
  
  TH1F *fMultDist5 = new TH1F("fMultDist5","Multiplicity Distribution",kNbinsM,-.5,kNbinsM-.5);
  fMultDist5->GetXaxis()->SetTitle("Multiplicity");
  fOutputList->Add(fMultDist5);
  

  TH3F *fPtEtaDist = new TH3F("fPtEtaDist","PtEtaDist",2,-1.1,1.1, 300,0,3., 28,-1.4,1.4);
  fOutputList->Add(fPtEtaDist);

  TH3F *fPhiPtDist = new TH3F("fPhiPtDist","PhiPtDist",2,-1.1,1.1, 120,0,2*PI, 300,0,3.);
  fOutputList->Add(fPhiPtDist);
  
  TH3F *fXi = new TH3F("fXi","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  fOutputList->Add(fXi);
  TH3F *fXibar = new TH3F("fXibar","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  fOutputList->Add(fXibar);

 
  TH3F *fXiMinusPiPlus  = new TH3F("fXiMinusPiPlus","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fXiMinusPiMinus = new TH3F("fXiMinusPiMinus","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fXiPlusPiPlus   = new TH3F("fXiPlusPiPlus","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fXiPlusPiMinus  = new TH3F("fXiPlusPiMinus","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);

  TH3F *fXiMinusPiPlusbkg  = new TH3F("fXiMinusPiPlusbkg","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fXiMinusPiMinusbkg = new TH3F("fXiMinusPiMinusbkg","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fXiPlusPiPlusbkg   = new TH3F("fXiPlusPiPlusbkg","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fXiPlusPiMinusbkg  = new TH3F("fXiPlusPiMinusbkg","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);

  fOutputList->Add(fXiMinusPiPlus);
  fOutputList->Add(fXiMinusPiMinus);
  fOutputList->Add(fXiPlusPiPlus);
  fOutputList->Add(fXiPlusPiMinus);

  fOutputList->Add(fXiMinusPiPlusbkg);
  fOutputList->Add(fXiMinusPiMinusbkg);
  fOutputList->Add(fXiPlusPiPlusbkg);
  fOutputList->Add(fXiPlusPiMinusbkg);


  
  TH3F *fMCinputXiStar = new TH3F("fMCinputXiStar","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fMCinputXiStarbar = new TH3F("fMCinputXiStarbar","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  fOutputList->Add(fMCinputXiStar);
  fOutputList->Add(fMCinputXiStarbar);

  TH3F *fMCinputXi = new TH3F("fMCinputXi","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  TH3F *fMCinputXibar = new TH3F("fMCinputXibar","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  fOutputList->Add(fMCinputXi);
  fOutputList->Add(fMCinputXibar);
  
  //

  TH3F *fMCinputTotalXiStar1 = new TH3F("fMCinputTotalXiStar1","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fMCinputTotalXiStarbar1 = new TH3F("fMCinputTotalXiStarbar1","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  fOutputList->Add(fMCinputTotalXiStar1);
  fOutputList->Add(fMCinputTotalXiStarbar1);

  TH3F *fMCinputTotalXi1 = new TH3F("fMCinputTotalXi1","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  TH3F *fMCinputTotalXibar1 = new TH3F("fMCinputTotalXibar1","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  fOutputList->Add(fMCinputTotalXi1);
  fOutputList->Add(fMCinputTotalXibar1);

  //

  TH3F *fMCinputTotalXiStar3 = new TH3F("fMCinputTotalXiStar3","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fMCinputTotalXiStarbar3 = new TH3F("fMCinputTotalXiStarbar3","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  fOutputList->Add(fMCinputTotalXiStar3);
  fOutputList->Add(fMCinputTotalXiStarbar3);

  TH3F *fMCinputTotalXi3 = new TH3F("fMCinputTotalXi3","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  TH3F *fMCinputTotalXibar3 = new TH3F("fMCinputTotalXibar3","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  fOutputList->Add(fMCinputTotalXi3);
  fOutputList->Add(fMCinputTotalXibar3);

  // 

  TH3F *fMCrecXi = new TH3F("fMCrecXi","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  TH3F *fMCrecXibar = new TH3F("fMCrecXibar","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
  fOutputList->Add(fMCrecXi);
  fOutputList->Add(fMCrecXibar);

  TH3F *fMCrecXiMinusPiPlus  = new TH3F("fMCrecXiMinusPiPlus","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  TH3F *fMCrecXiPlusPiMinus  = new TH3F("fMCrecXiPlusPiMinus","Invariant Mass Distribution", 100,0,10, 40,-2,2, 1100,1.4,2.5);
  fOutputList->Add(fMCrecXiMinusPiPlus);
  fOutputList->Add(fMCrecXiPlusPiMinus);
  
  TH2F *fMCinputXiStarxiy = new TH2F("fMCinputXiStarxiy","y distribution",80,-2,2, 80,-2,2);
  TH2F *fMCinputXiStarpiony = new TH2F("fMCinputXiStarpiony","y distribution",80,-2,2, 80,-2,2);
  fOutputList->Add(fMCinputXiStarxiy);
  fOutputList->Add(fMCinputXiStarpiony);
  TH2F *fMCinputXilambday = new TH2F("fMCinputXilambday","y distribution",80,-2,2, 80,-2,2);
  TH2F *fMCinputXipiony = new TH2F("fMCinputXipiony","y distribution",80,-2,2, 80,-2,2);
  fOutputList->Add(fMCinputXilambday);
  fOutputList->Add(fMCinputXipiony);
  TH2F *fMCinputLamprotony = new TH2F("fMCinputLamprotony","y distribution",80,-2,2, 80,-2,2);
  TH2F *fMCinputLampiony = new TH2F("fMCinputLampiony","y distribution",80,-2,2, 80,-2,2);
  fOutputList->Add(fMCinputLamprotony);
  fOutputList->Add(fMCinputLampiony);
  TH2F *fMCrecXiStarxiy = new TH2F("fMCrecXiStarxiy","y distribution",80,-2,2, 80,-2,2);
  TH2F *fMCrecXiStarpiony = new TH2F("fMCrecXiStarpiony","y distribution",80,-2,2, 80,-2,2);
  fOutputList->Add(fMCrecXiStarxiy);
  fOutputList->Add(fMCrecXiStarpiony);
  TH2F *fMCrecXilambday = new TH2F("fMCrecXilambday","y distribution",80,-2,2, 80,-2,2);
  TH2F *fMCrecXipiony = new TH2F("fMCrecXipiony","y distribution",80,-2,2, 80,-2,2);
  fOutputList->Add(fMCrecXilambday);
  fOutputList->Add(fMCrecXipiony);
  TH2F *fMCrecLamprotony = new TH2F("fMCrecLamprotony","y distribution",80,-2,2, 80,-2,2);
  TH2F *fMCrecLampiony = new TH2F("fMCrecLampiony","y distribution",80,-2,2, 80,-2,2);
  fOutputList->Add(fMCrecLamprotony);
  fOutputList->Add(fMCrecLampiony);

  
  ///////////////////////////////////  
  
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliXiStar::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
  cout<<"===========  Event # "<<fEventCounter+1<<"  ==========="<<endl;
  fEventCounter++;

  if(fAODcase) {cout<<"AODs not fully supported! Exiting event."<<endl; return;}

  if(fAODcase) fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
  else fESD = dynamic_cast<AliESDEvent*> (InputEvent());
  
  if(fAODcase) {if (!fAOD) {Printf("ERROR: fAOD not available"); return;}}
  else {if (!fESD) {Printf("ERROR: fESD not available"); return;}}
  

  // ESD Trigger Cut
  if(!fAODcase){
    if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected())) {cout<<"Event Rejected"<<endl; return;}
  }
  
  ///////////////////////////////////////////////////////////
  const AliAODVertex *PrimaryVertexAOD;
  const AliESDVertex *PrimaryVertexESD;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  

  TClonesArray *mcArray       = 0x0;
  //AliAODMCParticle *mcXi;
  //AliAODMCParticle *mcXiStarD2;
  //AliAODMCParticle *mcXiStar;
  AliMCEvent  *mcEvent        = 0x0; 
  AliStack    *mcstack        = 0x0;
  TParticle   *MCLamD1esd     = 0x0;
  TParticle   *MCLamD2esd     = 0x0;
  TParticle   *MCLamesd       = 0x0;
  TParticle   *MCXiD2esd      = 0x0;
  TParticle   *MCXiesd        = 0x0;
  TParticle   *MCXiStarD2esd  = 0x0;
  TParticle   *MCXiStaresd    = 0x0;

  Double_t px1,py1,pz1, px2,py2,pz2;
  Double_t p1sq,p2sq,e1,e2,xiStarE,angle;
  Double_t dca3d;
  Float_t dca2[2];
  Double_t xiVtx[3], xiStarVtx[3];
  Double_t xiP[3], xiStarP[3];
  Double_t xiStarMom;
  Double_t xiMass, xiStarMass;
  Double_t xiPt, xiStarPt;
  Double_t xiY, xiStarY;
  Double_t xiCharge;
  Double_t decayLength, decayLengthXY;
  Double_t pDaughter1[3];
  Double_t pDaughter2[3];
  Double_t xDaughter1[3];
  Double_t xDaughter2[3];
  //
  Double_t bField=0;
  UInt_t status=0;
  Int_t positiveTracks=0, negativeTracks=0;
  Int_t myTracks=0;
  //
  Double_t primaryVtx[3]={0};
  Int_t mBin=0;
  Int_t zBin=0;
  Double_t zStep=2*10/Double_t(fZvertexBins), zStart=-10.;
  //
  Bool_t mcXiFilled=kFALSE;// So that mctracks are never used uninitialized

  if(fMCcase){
    if(fAODcase){ 
      mcArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
      if(!mcArray){
	cout<<"No MC particle branch found"<<endl;
	return;
      }
    }else {
      mcEvent = MCEvent();
      if (!mcEvent) {Printf("ERROR: Could not retrieve MC event"); return;}
      
      mcstack = mcEvent->Stack();
      if (!mcstack) {Printf("ERROR: Could not retrieve the stack"); return;}
    }
  }

 
  /////////////////////////////////////////////////
  
  
  
  if(fAODcase){// AOD case
    ////////////////////////////////
    // Vertexing
    ((TH1F*)fOutputList->FindObject("fMultDist1"))->Fill(fAOD->GetNumberOfTracks());
    PrimaryVertexAOD = fAOD->GetPrimaryVertex();
    if(!PrimaryVertexAOD) return;
    
    if(fMCcase){
      for (Int_t it = 0; it < mcArray->GetEntriesFast(); it++) {
	AliAODMCParticle *mcInputTrack = (AliAODMCParticle*)mcArray->At(it);
	if (!mcInputTrack) {
	  Error("UserExec", "Could not receive track %d", it);
	  continue;
	}
	
	if(!mcInputTrack->IsPhysicalPrimary()) continue;
	
	// Xi
	if(mcInputTrack->GetPdgCode() == +kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXi1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
	if(mcInputTrack->GetPdgCode() == -kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXibar1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
	
	// XiStar
	if(mcInputTrack->GetPdgCode() == +kXiStarCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXiStar1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
	if(mcInputTrack->GetPdgCode() == -kXiStarCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXiStarbar1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
      }
    }



    primaryVtx[0]=PrimaryVertexAOD->GetX(); primaryVtx[1]=PrimaryVertexAOD->GetY(); primaryVtx[2]=PrimaryVertexAOD->GetZ();
    ((TH3F*)fOutputList->FindObject("fVertexDist1"))->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);
    if(fabs(primaryVtx[2]) > 10.) return; // Z-Vertex Cut  
    ((TH1F*)fOutputList->FindObject("fMultDist2"))->Fill(fAOD->GetNumberOfTracks());
    

    if(fAOD->IsPileupFromSPD()) return; // Reject Pile-up events
    
    ((TH1F*)fOutputList->FindObject("fMultDist3"))->Fill(fAOD->GetNumberOfTracks());
    ((TH3F*)fOutputList->FindObject("fVertexDist3"))->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);

    if(PrimaryVertexAOD->GetNContributors() >= 1) ((TH1F*)fOutputList->FindObject("fMultDist4"))->Fill(fAOD->GetNumberOfTracks());
    
       

    Printf("There are %d tracks in this event", fAOD->GetNumberOfTracks());
    // fNtracks Cut
    if(fAOD->GetNumberOfTracks() > kNbinsM) {cout<<"More tracks than limit"<<endl; return;}
    bField = fAOD->GetMagneticField();
    
    // Track loop
    for (Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
      AliAODTrack* aodtrack = fAOD->GetTrack(i);
      if (!aodtrack) continue;
      
      status=aodtrack->GetStatus();
     
      
      if( (status&AliESDtrack::kTPCrefit)==0) continue;// Require tpcrefit
      if( aodtrack->GetTPCNcls() < 70) continue;// TPC Ncluster cut
      if(aodtrack->Pt() < 0.15) continue;
      AliAODVertex *VtxType=aodtrack->GetProdVertex();
      if((VtxType->GetType())==AliAODVertex::kKink) continue;// Reject Kinks

      Bool_t goodMomentum = aodtrack->GetPxPyPz( fTempStruct[myTracks].fP);
      if(!goodMomentum) continue;    
      aodtrack->GetXYZ( fTempStruct[myTracks].fX);
      
      
      aodtrack->GetCovarianceXYZPxPyPz( fTempStruct[myTracks].fCov);
      
      dca2[0] = sqrt( pow(fTempStruct[myTracks].fX[0] - primaryVtx[0],2) + pow(fTempStruct[myTracks].fX[1] - primaryVtx[1],2));
      dca2[1] = sqrt( pow(fTempStruct[myTracks].fX[2] - primaryVtx[2],2));
      dca3d = sqrt( pow(dca2[0],2) + pow(dca2[1],2));
      
      ((TH1F*)fOutputList->FindObject("fDCADist"))->Fill(fAOD->GetNumberOfTracks(), dca3d);
      ((TH3F*)fOutputList->FindObject("fPhiPtDist"))->Fill(aodtrack->Charge(), aodtrack->Phi(), aodtrack->Pt());
      ((TH3F*)fOutputList->FindObject("fPtEtaDist"))->Fill(aodtrack->Charge(), aodtrack->Pt(), aodtrack->Eta());
      
      
      
      
      fTempStruct[myTracks].fStatus = status;
      fTempStruct[myTracks].fFilterMap = aodtrack->GetFilterMap();
      fTempStruct[myTracks].fID = aodtrack->GetID();
      fTempStruct[myTracks].fLabel = aodtrack->GetLabel();
      fTempStruct[myTracks].fPhi = atan2(fTempStruct[myTracks].fP[1], fTempStruct[myTracks].fP[0]);
      if(fTempStruct[myTracks].fPhi < 0) fTempStruct[myTracks].fPhi += 2*PI;
      fTempStruct[myTracks].fPt = sqrt(pow(fTempStruct[myTracks].fP[0],2) + pow(fTempStruct[myTracks].fP[1],2));
      fTempStruct[myTracks].fMom = sqrt( pow(fTempStruct[myTracks].fPt,2) + pow(fTempStruct[myTracks].fP[2],2) );
      fTempStruct[myTracks].fEta = aodtrack->Eta();
      fTempStruct[myTracks].fCharge = aodtrack->Charge();
      fTempStruct[myTracks].fDCAXY = dca2[0];
      fTempStruct[myTracks].fDCAZ = dca2[1];
      fTempStruct[myTracks].fDCA = dca3d;
      fTempStruct[myTracks].fNSigmaPi = fabs(fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kPion));
      fTempStruct[myTracks].fNSigmaK = fabs(fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kKaon));
      fTempStruct[myTracks].fNSigmaPr = fabs(fPIDResponse->NumberOfSigmasTPC(aodtrack,AliPID::kProton));
         
      
      if(aodtrack->Charge() > 0) positiveTracks++;
      else negativeTracks++;
      
      
      myTracks++;
    }
  }else {// ESDs

    ((TH1F*)fOutputList->FindObject("fMultDist1"))->Fill(fESD->GetNumberOfTracks());
    PrimaryVertexESD = fESD->GetPrimaryVertex();
    if(!PrimaryVertexESD) return;

    primaryVtx[0]=PrimaryVertexESD->GetX(); primaryVtx[1]=PrimaryVertexESD->GetY(); primaryVtx[2]=PrimaryVertexESD->GetZ();
    ((TH3F*)fOutputList->FindObject("fVertexDist1"))->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);

    if(fMCcase){
      /////////////////////////////////////////////////
      // Lam mc input
      /////////////////////////////////////////////////
      for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
	TParticle *mcInputTrack = (TParticle*)mcstack->Particle(it);    
	if (!mcInputTrack) {
	  Error("UserExec", "Could not receive track %d", it);
	  continue;
	}
	
	
	// Xi
	if(mcInputTrack->GetPdgCode() == +kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXi1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
	if(mcInputTrack->GetPdgCode() == -kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXibar1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());

	// XiStar
	if(mcInputTrack->GetPdgCode() == +kXiStarCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXiStar1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
	if(mcInputTrack->GetPdgCode() == -kXiStarCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXiStarbar1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());

      }
      
      
    }
    

    
    if(fabs(primaryVtx[2]) > 10.) return; // Z-Vertex Cut  
    ((TH1F*)fOutputList->FindObject("fMultDist2"))->Fill(fESD->GetNumberOfTracks());
    

    if(fESD->IsPileupFromSPD()) return; // Reject Pile-up events
    
    ((TH1F*)fOutputList->FindObject("fMultDist3"))->Fill(fESD->GetNumberOfTracks());
    ((TH3F*)fOutputList->FindObject("fVertexDist3"))->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);

    if(PrimaryVertexESD->GetNContributors() >= 1) ((TH1F*)fOutputList->FindObject("fMultDist4"))->Fill(fESD->GetNumberOfTracks());
    
    Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
    // fNtracks Cut
    if(fESD->GetNumberOfTracks() > kNbinsM) {cout<<"More tracks than limit"<<endl; return;}
    bField = fESD->GetMagneticField();
    
    
    // Track loop
    for (Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
      AliESDtrack* esdtrack = fESD->GetTrack(i);
      if (!esdtrack) continue;
      
      status=esdtrack->GetStatus();
      
      if(!fTrackCut->AcceptTrack(esdtrack)) continue;
       
      Bool_t goodMomentum = esdtrack->GetPxPyPz( fTempStruct[myTracks].fP);
      if(!goodMomentum) continue;    
      esdtrack->GetXYZ( fTempStruct[myTracks].fX);
      
   
      esdtrack->GetCovarianceXYZPxPyPz( fTempStruct[myTracks].fCov);
      //esdtrack->GetImpactParameters(dca2, cov);
      dca2[0] = sqrt( pow(fTempStruct[myTracks].fX[0] - primaryVtx[0],2) + pow(fTempStruct[myTracks].fX[1] - primaryVtx[1],2));
      dca2[1] = sqrt( pow(fTempStruct[myTracks].fX[2] - primaryVtx[2],2));
      dca3d = sqrt( pow(dca2[0],2) + pow(dca2[1],2));
   
      ((TH1F*)fOutputList->FindObject("fDCADist"))->Fill(fESD->GetNumberOfTracks(), dca3d);
      ((TH3F*)fOutputList->FindObject("fPhiPtDist"))->Fill(esdtrack->Charge(), esdtrack->Phi(), esdtrack->Pt());
      ((TH3F*)fOutputList->FindObject("fPtEtaDist"))->Fill(esdtrack->Charge(), esdtrack->Pt(), esdtrack->Eta());
      
         
      
      fTempStruct[myTracks].fStatus = status;
      fTempStruct[myTracks].fID = esdtrack->GetID();
      fTempStruct[myTracks].fLabel = esdtrack->GetLabel();
      fTempStruct[myTracks].fPhi = atan2(fTempStruct[myTracks].fP[1], fTempStruct[myTracks].fP[0]);
      if(fTempStruct[myTracks].fPhi < 0) fTempStruct[myTracks].fPhi += 2*PI;
      fTempStruct[myTracks].fPt = sqrt(pow(fTempStruct[myTracks].fP[0],2) + pow(fTempStruct[myTracks].fP[1],2));
      fTempStruct[myTracks].fMom = sqrt( pow(fTempStruct[myTracks].fPt,2) + pow(fTempStruct[myTracks].fP[2],2) );
      fTempStruct[myTracks].fEta = esdtrack->Eta();
      fTempStruct[myTracks].fCharge = esdtrack->Charge();
      fTempStruct[myTracks].fDCAXY = dca2[0];
      fTempStruct[myTracks].fDCAZ = dca2[1];
      fTempStruct[myTracks].fDCA = dca3d;
      fTempStruct[myTracks].fNSigmaPi = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack,AliPID::kPion));
      fTempStruct[myTracks].fNSigmaK = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack,AliPID::kKaon));
      fTempStruct[myTracks].fNSigmaPr = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack,AliPID::kProton));
          
      
      if(esdtrack->Charge() > 0) positiveTracks++;
      else negativeTracks++;
      
      myTracks++;
    }

  }// end of ESD case
  
  
  
  if(myTracks >= 1) {
    ((TH1F*)fOutputList->FindObject("fMultDist5"))->Fill(myTracks);
    ((TH3F*)fOutputList->FindObject("fMultDist3d"))->Fill(positiveTracks+negativeTracks, positiveTracks, negativeTracks);
  }


  cout<<"There are "<<myTracks<<"  myTracks"<<endl;
  
  // set Z Vertex bin
  for(Int_t i=0; i<fZvertexBins; i++){
    if( (primaryVtx[2] > zStart+i*zStep) && (primaryVtx[2] < zStart+(i+1)*zStep) ){
      zBin=i;
      break;
    }
  }
  
  // set Multiplicity bin
  for(Int_t i=0; i<fMultBins; i++){
    if( ( myTracks > fMultLimits[i]) && ( myTracks <= fMultLimits[i+1]) ) { mBin=i; break;}
  }

  
    
  ////////////////////////////////////
  // Add event to buffer if > 0 tracks
  if(myTracks > 0){
    fEC[zBin][mBin]->FIFOShift();
    (fEvt) = fEC[zBin][mBin]->fEvtStr;
    (fEvt)->fNTracks = myTracks;
    for(Int_t i=0; i<myTracks; i++) (fEvt)->fTracks[i] = fTempStruct[i];
  }

  
   
  
  if(fMCcase && fAODcase){// get Input MC information for AOD case
        
    /////////////////////////////////////////////////
    // Xi mc input
    /////////////////////////////////////////////////
    for (Int_t it = 0; it < mcArray->GetEntriesFast(); it++) {
      AliAODMCParticle *mcInputTrackXi = (AliAODMCParticle*)mcArray->At(it);
          
      if (!mcInputTrackXi) {
	Error("UserExec", "Could not receive track %d", it);
	continue;
      }
      
      if(abs(mcInputTrackXi->GetPdgCode())!=kXiCode) continue;
      if(!mcInputTrackXi->IsPhysicalPrimary()) continue;


      if(mcInputTrackXi->GetPdgCode() == +kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXi3"))->Fill(mcInputTrackXi->Pt(), mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());
      else ((TH3F*)fOutputList->FindObject("fMCinputTotalXibar3"))->Fill(mcInputTrackXi->Pt(), mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());
      
      
      
      AliAODMCParticle *mcInputTrackXiD1 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXi->GetDaughter(0)));
      AliAODMCParticle *mcInputTrackXiD2 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXi->GetDaughter(1)));
      
      if(abs(mcInputTrackXiD1->GetPdgCode())!=kLambdaCode && abs(mcInputTrackXiD2->GetPdgCode())!=kLambdaCode) continue;
      if(abs(mcInputTrackXiD1->GetPdgCode())!=kPionCode && abs(mcInputTrackXiD2->GetPdgCode())!=kPionCode) continue;
      
      
      // At this point we have the right Xi decay channel
      
      AliAODMCParticle *mcInputTrackLamD1;
      AliAODMCParticle *mcInputTrackLamD2;
      
      if(abs(mcInputTrackXiD1->GetPdgCode())==kLambdaCode) {
	mcInputTrackLamD1 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiD1->GetDaughter(0)));
	mcInputTrackLamD2 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiD1->GetDaughter(1)));
      }
      else{
	mcInputTrackLamD1 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiD2->GetDaughter(0)));
	mcInputTrackLamD2 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiD2->GetDaughter(1)));
      }


      if(abs(mcInputTrackLamD1->GetPdgCode())!=kProtonCode && abs(mcInputTrackLamD2->GetPdgCode())!=kProtonCode) continue;
      if(abs(mcInputTrackLamD1->GetPdgCode())!=kPionCode && abs(mcInputTrackLamD2->GetPdgCode())!=kPionCode) continue;

      // At this point we have the right Lambda decay channel

      if(mcInputTrackXi->GetPdgCode() == +kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputXi"))->Fill(mcInputTrackXi->Pt(), mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());
      else ((TH3F*)fOutputList->FindObject("fMCinputXibar"))->Fill(mcInputTrackXi->Pt(), mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());

    }



    /////////////////////////////////////////////////
    // XiStar mc input
    /////////////////////////////////////////////////
    for (Int_t it = 0; it < mcArray->GetEntriesFast(); it++) {
      AliAODMCParticle *mcInputTrackXiStar = (AliAODMCParticle*)mcArray->At(it);
      if (!mcInputTrackXiStar) {
	Error("UserExec", "Could not receive track %d", it);
	continue;
      }
   
      if(abs(mcInputTrackXiStar->GetPdgCode())!=kXiStarCode) continue;
      if(!mcInputTrackXiStar->IsPhysicalPrimary()) continue;

      if(mcInputTrackXiStar->GetPdgCode() == +kXiStarCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXiStar3"))->Fill(mcInputTrackXiStar->Pt(), mcInputTrackXiStar->Y(), mcInputTrackXiStar->GetCalcMass());
      else ((TH3F*)fOutputList->FindObject("fMCinputTotalXiStarbar3"))->Fill(mcInputTrackXiStar->Pt(), mcInputTrackXiStar->Y(), mcInputTrackXiStar->GetCalcMass());

      
      AliAODMCParticle *mcInputTrackXiStarD1 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiStar->GetDaughter(0)));
      AliAODMCParticle *mcInputTrackXiStarD2 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiStar->GetDaughter(1)));
      
      if(abs(mcInputTrackXiStarD1->GetPdgCode())!=kXiCode && abs(mcInputTrackXiStarD2->GetPdgCode())!=kXiCode) continue;
      if(abs(mcInputTrackXiStarD1->GetPdgCode())!=kPionCode && abs(mcInputTrackXiStarD2->GetPdgCode())!=kPionCode) continue;
      

      // At this point we have the right Xi(1530) decay channel

      AliAODMCParticle *mcInputTrackXiD1;
      AliAODMCParticle *mcInputTrackXiD2;
      if(abs(mcInputTrackXiStarD1->GetPdgCode())==kXiCode) {
	mcInputTrackXiD1 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiStarD1->GetDaughter(0)));
	mcInputTrackXiD2 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiStarD1->GetDaughter(1)));
      }
      else{
	mcInputTrackXiD1 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiStarD2->GetDaughter(0)));
	mcInputTrackXiD2 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiStarD2->GetDaughter(1)));
      }
      
      if(abs(mcInputTrackXiD1->GetPdgCode())!=kLambdaCode && abs(mcInputTrackXiD2->GetPdgCode())!=kLambdaCode) continue;
      if(abs(mcInputTrackXiD1->GetPdgCode())!=kPionCode && abs(mcInputTrackXiD2->GetPdgCode())!=kPionCode) continue;
      

      // At this point we have the right Xi decay channel

      AliAODMCParticle *mcInputTrackLamD1;
      AliAODMCParticle *mcInputTrackLamD2;

      if(abs(mcInputTrackXiD1->GetPdgCode())==kLambdaCode) {
	mcInputTrackLamD1 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiD1->GetDaughter(0)));
	mcInputTrackLamD2 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiD1->GetDaughter(1)));
      }
      else{
	mcInputTrackLamD1 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiD2->GetDaughter(0)));
	mcInputTrackLamD2 = (AliAODMCParticle*)mcArray->At(abs(mcInputTrackXiD2->GetDaughter(1)));
      }

      if(abs(mcInputTrackLamD1->GetPdgCode())!=kProtonCode && abs(mcInputTrackLamD2->GetPdgCode())!=kProtonCode) continue;
      if(abs(mcInputTrackLamD1->GetPdgCode())!=kPionCode && abs(mcInputTrackLamD2->GetPdgCode())!=kPionCode) continue;

      // At this point we the right Lambda decay channel
      
      if(mcInputTrackXiStar->GetPdgCode() == +kXiStarCode) ((TH3F*)fOutputList->FindObject("fMCinputXiStar"))->Fill(mcInputTrackXiStar->Pt(), mcInputTrackXiStar->Y(), mcInputTrackXiStar->GetCalcMass());
      else ((TH3F*)fOutputList->FindObject("fMCinputXiStarbar"))->Fill(mcInputTrackXiStar->Pt(), mcInputTrackXiStar->Y(), mcInputTrackXiStar->GetCalcMass());
      
      if(abs(mcInputTrackXiStarD1->GetPdgCode())==kXiCode) {
	((TH2F*)fOutputList->FindObject("fMCinputXiStarxiy"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackXiStarD1->Y());
	((TH2F*)fOutputList->FindObject("fMCinputXiStarpiony"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackXiStarD2->Y());
      }else{
	((TH2F*)fOutputList->FindObject("fMCinputXiStarxiy"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackXiStarD2->Y());
	((TH2F*)fOutputList->FindObject("fMCinputXiStarpiony"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackXiStarD1->Y());
      }
      if(abs(mcInputTrackXiD1->GetPdgCode())==kLambdaCode){
	((TH2F*)fOutputList->FindObject("fMCinputXilambday"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackXiD1->Y());
	((TH2F*)fOutputList->FindObject("fMCinputXipiony"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackXiD2->Y());
      }else{
	((TH2F*)fOutputList->FindObject("fMCinputXilambday"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackXiD2->Y());
	((TH2F*)fOutputList->FindObject("fMCinputXipiony"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackXiD1->Y());
      }
      if(abs(mcInputTrackLamD1->GetPdgCode())==kProtonCode){
	((TH2F*)fOutputList->FindObject("fMCinputLamprotony"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackLamD1->Y());
	((TH2F*)fOutputList->FindObject("fMCinputLampiony"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackLamD2->Y());
      }else {
	((TH2F*)fOutputList->FindObject("fMCinputLamprotony"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackLamD2->Y());
	((TH2F*)fOutputList->FindObject("fMCinputLampiony"))->Fill(mcInputTrackXiStar->Y(), mcInputTrackLamD1->Y());
      }
      
      
    }
  }
  

  
  if(fMCcase && !fAODcase){// get Input MC information for ESD case

    /////////////////////////////////////////////////
    // Xi mc input
    /////////////////////////////////////////////////
    for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
      TParticle *mcInputTrackXi = (TParticle*)mcstack->Particle(it);    
      if (!mcInputTrackXi) {
	Error("UserExec", "Could not receive track %d", it);
	continue;
      }
        
      //if(!mcstack->IsPhysicalPrimary(it)) continue;
      if(abs(mcInputTrackXi->GetPdgCode())!=kXiCode) continue;
   

      if(mcInputTrackXi->GetPdgCode() == +kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXi3"))->Fill(mcInputTrackXi->Pt(), mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());
      else ((TH3F*)fOutputList->FindObject("fMCinputTotalXibar3"))->Fill(mcInputTrackXi->Pt(), mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());
    
    }
    
    
    /////////////////////////////////////////////////
    // XiStar mc input
    /////////////////////////////////////////////////
    for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
      TParticle *mcInputTrackXiStar = (TParticle*)mcstack->Particle(it);
      if (!mcInputTrackXiStar) {
	Error("UserExec", "Could not receive track %d", it);
	continue;
      }
      
      //if(!mcstack->IsPhysicalPrimary(it)) continue;
      if(abs(mcInputTrackXiStar->GetPdgCode())!=kXiStarCode) continue;
      
      if(mcInputTrackXiStar->GetPdgCode() == +kXiStarCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXiStar3"))->Fill(mcInputTrackXiStar->Pt(), mcInputTrackXiStar->Y(), mcInputTrackXiStar->GetCalcMass());
      else ((TH3F*)fOutputList->FindObject("fMCinputTotalXiStarbar3"))->Fill(mcInputTrackXiStar->Pt(), mcInputTrackXiStar->Y(), mcInputTrackXiStar->GetCalcMass());

    }
  }
  
  

  if(fAODcase) {cout<<"AOD XiVertexer not currently supported! Exiting event"<<endl; return;}

  ////////////////////////////////////////////////
  // Reconstruction
  for(Int_t i=0; i<fESD->GetNumberOfCascades(); i++){
    
    AliESDcascade *Xicandidate = fESD->GetCascade(i);
    
    if(TMath::Abs( Xicandidate->GetPindex()) == TMath::Abs( Xicandidate->GetNindex())) continue;
    if(TMath::Abs( Xicandidate->GetPindex()) == TMath::Abs( Xicandidate->GetBindex())) continue;
    if(TMath::Abs( Xicandidate->GetNindex()) == TMath::Abs( Xicandidate->GetBindex())) continue;

    AliESDtrack *pTrackXi	= fESD->GetTrack(TMath::Abs( Xicandidate->GetPindex()));
    AliESDtrack *nTrackXi	= fESD->GetTrack(TMath::Abs( Xicandidate->GetNindex()));
    AliESDtrack *bTrackXi	= fESD->GetTrack(TMath::Abs( Xicandidate->GetBindex()));
    
    // Standard track QA cuts
    if(!fTrackCut->AcceptTrack(pTrackXi)) continue;
    if(!fTrackCut->AcceptTrack(nTrackXi)) continue;
    if(!fTrackCut->AcceptTrack(bTrackXi)) continue;
    
  
    if(fabs(Xicandidate->GetDcaV0Daughters()) > fDCAProtonPion) continue; 
    if(fabs(Xicandidate->GetD(primaryVtx[0],primaryVtx[1],primaryVtx[2])) < fDCAVtxLambda) continue;
  
    if(Xicandidate->Charge() == -1){
      if(fabs(pTrackXi->GetD(primaryVtx[0],primaryVtx[1],bField)) < fDCAVtxProton) continue;
      if(fabs(nTrackXi->GetD(primaryVtx[0],primaryVtx[1],bField)) < fDCAVtxPionFirst) continue;
    }else{
      if(fabs(pTrackXi->GetD(primaryVtx[0],primaryVtx[1],bField)) < fDCAVtxPionFirst) continue;
      if(fabs(nTrackXi->GetD(primaryVtx[0],primaryVtx[1],bField)) < fDCAVtxProton) continue;
    }
  
  
    if(fabs(bTrackXi->GetD(primaryVtx[0],primaryVtx[1],bField)) < fDCAVtxPionSecond) continue;
    if(fabs(Xicandidate->GetDcaXiDaughters()) > fDCALambdaPion) continue;
    
    Double_t tempX[3]={0};
    Xicandidate->GetXYZ(tempX[0], tempX[1], tempX[2]);
    if(sqrt( pow(tempX[0],2) + pow(tempX[1],2) ) < fLambdaDecayLengthXY) continue;
    if(sqrt( pow(tempX[0],2) + pow(tempX[1],2) ) > fMaxDecayLength) continue;
    
    
    if(Xicandidate->GetCascadeCosineOfPointingAngle(primaryVtx[0],primaryVtx[1],primaryVtx[2]) < fXiCosTheta) continue;
    if(Xicandidate->GetV0CosineOfPointingAngle(primaryVtx[0],primaryVtx[1],primaryVtx[2]) < fLamCosTheta) continue;
    

    xiP[0] = Xicandidate->Px();
    xiP[1] = Xicandidate->Py();
    xiP[2] = Xicandidate->Pz();
    xiVtx[0] = Xicandidate->Xv();
    xiVtx[1] = Xicandidate->Yv();
    xiVtx[2] = Xicandidate->Zv();
    xiPt = Xicandidate->Pt();
    xiY = Xicandidate->RapXi();
    xiMass = Xicandidate->M();
    xiCharge = Xicandidate->Charge();

    decayLengthXY = sqrt( pow(xiVtx[0]-primaryVtx[0],2) + pow(xiVtx[1]-primaryVtx[1],2) );
    if(decayLengthXY < fXiDecayLengthXY) continue;// 2d version
    if(decayLengthXY > fMaxDecayLength) continue;// 2d version
    
    

    if(xiCharge == -1) ((TH3F*)fOutputList->FindObject("fXi"))->Fill(xiPt, xiY, xiMass);
    else ((TH3F*)fOutputList->FindObject("fXibar"))->Fill(xiPt, xiY, xiMass);
    
           
    // MC associaton
    mcXiFilled = kFALSE;
    if(fMCcase && !fAODcase){
      
      MCXiD2esd = (TParticle*)mcstack->Particle(abs(bTrackXi->GetLabel()));
      
      if(abs(MCXiD2esd->GetPdgCode())==kPionCode){
	
	MCLamD1esd = (TParticle*)mcstack->Particle(abs(pTrackXi->GetLabel()));
	MCLamD2esd = (TParticle*)mcstack->Particle(abs(nTrackXi->GetLabel()));
	
	if(MCLamD1esd->GetMother(0) == MCLamD2esd->GetMother(0)){
	  if(abs(MCLamD1esd->GetPdgCode())==kProtonCode || abs(MCLamD2esd->GetPdgCode())==kProtonCode) {
	    if(abs(MCLamD1esd->GetPdgCode())==kPionCode || abs(MCLamD2esd->GetPdgCode())==kPionCode) {
	      
	      MCLamesd = (TParticle*)mcstack->Particle(abs(MCLamD1esd->GetMother(0)));
	      if(abs(MCLamesd->GetPdgCode())==kLambdaCode) {
		
		if(MCLamesd->GetMother(0) == MCXiD2esd->GetMother(0)){
		  MCXiesd = (TParticle*)mcstack->Particle(abs(MCLamesd->GetMother(0)));
		  if(abs(MCXiesd->GetPdgCode())==kXiCode) {
		    
		    if(Xicandidate->Charge() == -1) {
		      ((TH3F*)fOutputList->FindObject("fMCrecXi"))->Fill(xiPt, xiY, xiMass);
		    }else {
		      ((TH3F*)fOutputList->FindObject("fMCrecXibar"))->Fill(xiPt, xiY, xiMass);
		      
		    }
		    mcXiFilled = kTRUE;
		  }
		}
	      }
	    }
	  }
	}
      }
    }// MC association
    
    
    if(fabs(xiMass-fTrueMassXi) > fMassWindow) continue;
    
    
  
    fXiTrack->Set(xiVtx, xiP, fCovMatrix, Short_t(xiCharge));
    
    
    //////////////////////////////////////////////////////////
    // Reconstruct Xi(1530)
    for(Int_t EN=0; EN<fEventsToMix+1; EN++){// Event buffer loop
      
      for(Int_t l=0; l<(fEvt+EN)->fNTracks; l++){// Present(EN=0) and Past(EN from 1 to fEventsToMix) event track loop
	
	if(EN==0) {
	  if((fEvt+EN)->fTracks[l].fID == pTrackXi->GetID()) continue;
	  if((fEvt+EN)->fTracks[l].fID == nTrackXi->GetID()) continue;
	  if((fEvt+EN)->fTracks[l].fID == bTrackXi->GetID()) continue;
	}
	
	fXiTrack->Set(xiVtx, xiP, fCovMatrix, Short_t(xiCharge));
	//cout<<(fEvt+EN)->fNTracks<<"  "<<(fEvt+EN)->fTracks[l].fCharge<<endl;
	fESDTrack4->Set((fEvt+EN)->fTracks[l].fX, (fEvt+EN)->fTracks[l].fP, (fEvt+EN)->fTracks[l].fCov, (fEvt+EN)->fTracks[l].fCharge);
	if(!fESDTrack4) continue;
	if(fAODcase){
	  if((Bool_t)(((1<<5) & (fEvt+EN)->fTracks[l].fFilterMap) == 0)) continue;// AOD filterbit cut, "Standard cuts with tight dca"
	}else{
	  if((fEvt+EN)->fTracks[l].fDCAXY > (0.0182 + 0.035/pow((fEvt+EN)->fTracks[l].fPt,1.01))) continue; 
	  if((fEvt+EN)->fTracks[l].fDCAZ > 2) continue;
	  if( (((fEvt+EN)->fTracks[l].fStatus)&AliESDtrack::kITSrefit)==0) continue;// Require itsrefit
	  // no Chi^2 cut applied for ESDs.  Info not available in my track structure.
	}
	
	if((fEvt+EN)->fTracks[l].fNSigmaPi > fSigmaCutPionThird) continue;
	if(fabs((fEvt+EN)->fTracks[l].fEta) > 0.8) continue;
	
	
	AliVertex *XiStarVtx = new AliVertex((fEvt+EN)->fTracks[l].fX,0,0);
	//fESDTrack4->PropagateToDCA(fXiTrack, bField);// Propagate tracks to dca, both tracks are budged
	if(!(fXiTrack->PropagateToDCA(XiStarVtx, bField, 3))) continue;// Propagate tracks to dca, version which assumes fESDTrack4 is already primary
	/////////////
	fXiTrack->GetPxPyPz(pDaughter1);
	fXiTrack->GetXYZ(xDaughter1);
	fESDTrack4->GetPxPyPz(pDaughter2);
	fESDTrack4->GetXYZ(xDaughter2);
	//////////////////////////
	
	
	
	xiStarVtx[0] = (xDaughter1[0]+xDaughter2[0])/2.;
	xiStarVtx[1] = (xDaughter1[1]+xDaughter2[1])/2.;
	xiStarVtx[2] = (xDaughter1[2]+xDaughter2[2])/2.;
	decayLength = sqrt(pow(xiStarVtx[0]-primaryVtx[0],2)+pow(xiStarVtx[1]-primaryVtx[1],2)+pow(xiStarVtx[2]-primaryVtx[2],2));
	
	px1=pDaughter1[0];
	py1=pDaughter1[1];
	pz1=pDaughter1[2];
	px2=pDaughter2[0];
	py2=pDaughter2[1];
	pz2=pDaughter2[2];
	
	p1sq=px1*px1+py1*py1+pz1*pz1;
	p2sq=px2*px2+py2*py2+pz2*pz2;
	if(p1sq <=0 || p2sq <=0) continue;
	
	e1=sqrt(p1sq+fTrueMassXi*fTrueMassXi);
	e2=sqrt(p2sq+fTrueMassPi*fTrueMassPi);
	angle=px1*px2+py1*py2+pz1*pz2;
	xiStarMass=fTrueMassXi*fTrueMassXi+fTrueMassPi*fTrueMassPi+2.*e1*e2-2.*angle;
	if(xiStarMass<0.) xiStarMass=1.e-8;
	xiStarMass=sqrt(xiStarMass);
	
	
	xiStarP[0] = px1+px2;
	xiStarP[1] = py1+py2;
	xiStarP[2] = pz1+pz2;
	xiStarMom = sqrt(pow(xiStarP[0],2)+pow(xiStarP[1],2)+pow(xiStarP[2],2));
	if(xiStarMom==0) continue; // So one of the following lines doesnt break
	xiStarPt = sqrt(xiStarP[0]*xiStarP[0] + xiStarP[1]*xiStarP[1]);
	xiStarY = .5*log( ((e1+e2) + xiStarP[2])/((e1+e2) - xiStarP[2]));
	xiStarE = e1 + e2;
	
	
	if( (xiStarP[0]*(xiStarVtx[0]-primaryVtx[0]) + xiStarP[1]*(xiStarVtx[1]-primaryVtx[1]) + xiStarP[2]*(xiStarVtx[2]-primaryVtx[2]))/xiStarMom/decayLength < fXiStarCosTheta) continue;

	
	if(EN==0){
	  if(fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == -1) ((TH3F*)fOutputList->FindObject("fXiMinusPiMinus"))->Fill(xiStarPt, xiStarY, xiStarMass);
	  else if(fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) ((TH3F*)fOutputList->FindObject("fXiMinusPiPlus"))->Fill(xiStarPt, xiStarY, xiStarMass);
	  else if(fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) ((TH3F*)fOutputList->FindObject("fXiPlusPiMinus"))->Fill(xiStarPt, xiStarY, xiStarMass);
	  else ((TH3F*)fOutputList->FindObject("fXiPlusPiPlus"))->Fill(xiStarPt, xiStarY, xiStarMass);
	}else {
	  if(fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == -1) ((TH3F*)fOutputList->FindObject("fXiMinusPiMinusbkg"))->Fill(xiStarPt, xiStarY, xiStarMass);
	  else if(fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) ((TH3F*)fOutputList->FindObject("fXiMinusPiPlusbkg"))->Fill(xiStarPt, xiStarY, xiStarMass);
	  else if(fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) ((TH3F*)fOutputList->FindObject("fXiPlusPiMinusbkg"))->Fill(xiStarPt, xiStarY, xiStarMass);
	  else ((TH3F*)fOutputList->FindObject("fXiPlusPiPlusbkg"))->Fill(xiStarPt, xiStarY, xiStarMass);
	}
	
	
	/*
	// MC associaton
	if(fMCcase && mcXiFilled && EN==0 && fAODcase){// AOD MC's
	  
	  MCXiStarD2 = (AliAODMCParticle*)mcArray->At(abs((fEvt)->fTracks[l].fLabel));
	  
	  if(abs(MCXiStarD2->GetPdgCode())==kPionCode){ 
	    if(MCXi->GetMother() == MCXiStarD2->GetMother()){
	      MCXiStar = (AliAODMCParticle*)mcArray->At(MCXi->GetMother());
	      if(abs(MCXiStar->GetPdgCode())==kXiStarCode) {
		
		if(fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) fMCrecXiMinusPiPlus->Fill(xiStarPt, xiStarY, xiStarMass);
		if(fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) fMCrecXiPlusPiMinus->Fill(xiStarPt, xiStarY, xiStarMass);
					
	      }
	    }
	  }
	}
	*/
	
	// MC associaton
	if(fMCcase && mcXiFilled && EN==0 && !fAODcase){// ESD MC's
	  MCXiStarD2esd = (TParticle*)mcstack->Particle(abs((fEvt)->fTracks[l].fLabel));
	  
	  if(abs(MCXiStarD2esd->GetPdgCode())==kPionCode){ 
	    if(MCXiesd->GetMother(0) == MCXiStarD2esd->GetMother(0)){
	          
	      MCXiStaresd = (TParticle*)mcstack->Particle(abs(MCXiesd->GetMother(0)));
	      if(abs(MCXiStaresd->GetPdgCode())==kXiStarCode) {
		
		if(fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) ((TH3F*)fOutputList->FindObject("fMCrecXiMinusPiPlus"))->Fill(xiStarPt, xiStarY, xiStarMass);
		if(fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) ((TH3F*)fOutputList->FindObject("fMCrecXiPlusPiMinus"))->Fill(xiStarPt, xiStarY, xiStarMass);
	      
	      }
	    }
	  }
	}
	
	
      }// 3rd pion loop
    }// Event mixing loop
    
  }// Xi loop
  

  
  // Post output data.
  PostData(1, fOutputList);
  
}
//________________________________________________________________________
void AliXiStar::Terminate(Option_t *) 
{
  cout<<"Done"<<endl;
}
//________________________________________________________________________
Double_t AliXiStar::LinearPropagateToDCA(AliESDtrack *v, AliESDtrack *t, Double_t b) {// Adapted from AliCascadeVertexer.cxx
 //--------------------------------------------------------------------
  // This function returns the DCA between the V0 and the track
  //--------------------------------------------------------------------

  Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
  Double_t r[3]; t->GetXYZ(r);
  Double_t x1=r[0], y1=r[1], z1=r[2];
  Double_t p[3]; t->GetPxPyPz(p);
  Double_t px1=p[0], py1=p[1], pz1=p[2];

  Double_t x2[3]={0};
  Double_t p2[3]={0};
  Double_t vx2,vy2,vz2;     // position and momentum of V0
  Double_t px2,py2,pz2;
  
  v->GetXYZ(x2);
  v->GetPxPyPz(p2);
  vx2=x2[0], vy2=x2[1], vz2=x2[2];
  px2=p2[0], py2=p2[1], pz2=p2[2];

// calculation dca
   
  Double_t dd= Det(vx2-x1,vy2-y1,vz2-z1,px1,py1,pz1,px2,py2,pz2);
  Double_t ax= Det(py1,pz1,py2,pz2);
  Double_t ay=-Det(px1,pz1,px2,pz2);
  Double_t az= Det(px1,py1,px2,py2);

  Double_t dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);

//points of the DCA
  Double_t t1 = Det(vx2-x1,vy2-y1,vz2-z1,px2,py2,pz2,ax,ay,az)/
                Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
  
  x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;
  

  //propagate track to the points of DCA

  x1=x1*cs1 + y1*sn1;
  if (!t->PropagateTo(x1,b)) {
    Error("PropagateToDCA","Propagation failed !");
    return 1.e+33;
  }  

  return dca;
}


//________________________________________________________________________
Double_t AliXiStar::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const {// Taken from AliCascadeVertexer
  //--------------------------------------------------------------------
  // This function calculates locally a 2x2 determinant
  //--------------------------------------------------------------------
  return a00*a11 - a01*a10;
}
//________________________________________________________________________
Double_t AliXiStar::Det(Double_t a00,Double_t a01,Double_t a02,
				 Double_t a10,Double_t a11,Double_t a12,
			      Double_t a20,Double_t a21,Double_t a22) const {// Taken from AliCascadeVertexer
  //--------------------------------------------------------------------
  // This function calculates locally a 3x3 determinant
  //--------------------------------------------------------------------
  return  a00*Det(a11,a12,a21,a22)-a01*Det(a10,a12,a20,a22)+a02*Det(a10,a11,a20,a21);
}



