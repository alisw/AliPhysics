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
// AlidNdPtCutAnalysis class.
//
// a. functionality:
// - fills generic cut histograms
// - generates cuts (selection criteria)
//
// b. data members:
// - generic cut histograms
// - control histograms
//
// Author: J.Otwinowski 04/11/2008
//
// c. flag to switch on/off THnSparse filling
// d. store information in less dimension THx histograms
//
// Modified: J.Otwinowski 21/10/2015
//------------------------------------------------------------------------------
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliTracker.h"

#include "AlidNdPtEventCuts.h"
#include "AlidNdPtAcceptanceCuts.h"
#include "AlidNdPtBackgroundCuts.h"
#include "AlidNdPtAnalysis.h"
#include "AliPhysicsSelection.h"

#include "AliPWG0Helper.h"
#include "AlidNdPtHelper.h"
#include "AlidNdPtCutAnalysis.h"

using namespace std;

ClassImp(AlidNdPtCutAnalysis)

//_____________________________________________________________________________
  AlidNdPtCutAnalysis::AlidNdPtCutAnalysis(): AlidNdPt(),
  fFillSparseHisto(kFALSE),
  fAnalysisFolder(0),
  fEventCount(0),
  fRecEventHist(0),
  fMCEventHist(0),
  fRecMCEventHist(0),
  fRecMCTrackHist(0),
  //
  fTrigVertCount(0),
  fRecEventXYZ(0),
  fRecEventXYMult(0),
  fRecEventXZMult(0),
  fRecEventYZMult(0),
  fRecEventZResZMult(0),
  fMCEventXYZ(0),
  fRecMCEventDXDYDZ(0),
  fRecMCEventDXDYMult(0),
  fRecMCEventDXDZMult(0),
  fRecMCEventDYDZMult(0),
  //
  fCrossRowsEtaPt(0),
  fChi2PerClustEtaPt(0),
  fCrossRowsOverFindEtaPt(0),
  fFracSharedClustEtaPt(0),
  fDCAyEtaPt(0),
  fDCAzEtaPt(0),
  //
  fCrossRowsPhiPt(0),
  fChi2PerClustPhiPt(0),
  fCrossRowsOverFindPhiPt(0),
  fFracSharedClustPhiPt(0),
  fDCAyPhiPt(0),
  fDCAzPhiPt(0),
  //
  fCrossRowsEtaPhi(0),
  fChi2PerClustEtaPhi(0),
  fCrossRowsOverFindEtaPhi(0),
  fFracSharedClustEtaPhi(0),
  fDCAyEtaPhi(0),
  fDCAzEtaPhi(0),
  //
  fDCAyEtaPtMCPrim(0),
  fDCAzEtaPtMCPrim(0),
  fDCAyPhiPtMCPrim(0),
  fDCAzPhiPtMCPrim(0),
  fDCAyEtaPhiMCPrim(0),
  fDCAzEtaPhiMCPrim(0),
  //
  fDCAyEtaPtMCSec(0),
  fDCAzEtaPtMCSec(0),
  fDCAyPhiPtMCSec(0),
  fDCAzPhiPtMCSec(0),
  fDCAyEtaPhiMCSec(0),
  fDCAzEtaPhiMCSec(0),
  //
  fDCAyEtaPtSecMCDecays(0),
  fDCAzEtaPtSecMCDecays(0),
  fDCAyPhiPtSecMCDecays(0),
  fDCAzPhiPtSecMCDecays(0),
  fDCAyEtaPhiSecMCDecays(0),
  fDCAzEtaPhiSecMCDecays(0),
  //
  fDCAyEtaPtSecMCDecaysK0s(0),
  fDCAzEtaPtSecMCDecaysK0s(0),
  fDCAyPhiPtSecMCDecaysK0s(0),
  fDCAzPhiPtSecMCDecaysK0s(0),
  fDCAyEtaPhiSecMCDecaysK0s(0),
  fDCAzEtaPhiSecMCDecaysK0s(0),
  //
  fDCAyEtaPtSecMCDecaysLambda(0),
  fDCAzEtaPtSecMCDecaysLambda(0),
  fDCAyPhiPtSecMCDecaysLambda(0),
  fDCAzPhiPtSecMCDecaysLambda(0),
  fDCAyEtaPhiSecMCDecaysLambda(0),
  fDCAzEtaPhiSecMCDecaysLambda(0),
  //
  fDCAyEtaPtSecMCMaterial(0),
  fDCAzEtaPtSecMCMaterial(0),
  fDCAyPhiPtSecMCMaterial(0),
  fDCAzPhiPtSecMCMaterial(0),
  fDCAyEtaPhiSecMCMaterial(0),
  fDCAzEtaPhiSecMCMaterial(0)
{
  // default constructor
  Init();
}

//_____________________________________________________________________________
AlidNdPtCutAnalysis::AlidNdPtCutAnalysis(Char_t* name, Char_t* title): AlidNdPt(name,title),
  fFillSparseHisto(kFALSE),
  fAnalysisFolder(0),
  fEventCount(0),
  fRecEventHist(0),
  fMCEventHist(0),
  fRecMCEventHist(0),
  fRecMCTrackHist(0),
  //
  fTrigVertCount(0),
  fRecEventXYZ(0),
  fRecEventXYMult(0),
  fRecEventXZMult(0),
  fRecEventYZMult(0),
  fRecEventZResZMult(0),
  fMCEventXYZ(0),
  fRecMCEventDXDYDZ(0),
  fRecMCEventDXDYMult(0),
  fRecMCEventDXDZMult(0),
  fRecMCEventDYDZMult(0),
  //
  fCrossRowsEtaPt(0),
  fChi2PerClustEtaPt(0),
  fCrossRowsOverFindEtaPt(0),
  fFracSharedClustEtaPt(0),
  fDCAyEtaPt(0),
  fDCAzEtaPt(0),
  //
  fCrossRowsPhiPt(0),
  fChi2PerClustPhiPt(0),
  fCrossRowsOverFindPhiPt(0),
  fFracSharedClustPhiPt(0),
  fDCAyPhiPt(0),
  fDCAzPhiPt(0),
  //
  fCrossRowsEtaPhi(0),
  fChi2PerClustEtaPhi(0),
  fCrossRowsOverFindEtaPhi(0),
  fFracSharedClustEtaPhi(0),
  fDCAyEtaPhi(0),
  fDCAzEtaPhi(0),
  //
  fDCAyEtaPtMCPrim(0),
  fDCAzEtaPtMCPrim(0),
  fDCAyPhiPtMCPrim(0),
  fDCAzPhiPtMCPrim(0),
  fDCAyEtaPhiMCPrim(0),
  fDCAzEtaPhiMCPrim(0),
  //
  fDCAyEtaPtMCSec(0),
  fDCAzEtaPtMCSec(0),
  fDCAyPhiPtMCSec(0),
  fDCAzPhiPtMCSec(0),
  fDCAyEtaPhiMCSec(0),
  fDCAzEtaPhiMCSec(0),
  //
  fDCAyEtaPtSecMCDecays(0),
  fDCAzEtaPtSecMCDecays(0),
  fDCAyPhiPtSecMCDecays(0),
  fDCAzPhiPtSecMCDecays(0),
  fDCAyEtaPhiSecMCDecays(0),
  fDCAzEtaPhiSecMCDecays(0),
  //
  fDCAyEtaPtSecMCDecaysK0s(0),
  fDCAzEtaPtSecMCDecaysK0s(0),
  fDCAyPhiPtSecMCDecaysK0s(0),
  fDCAzPhiPtSecMCDecaysK0s(0),
  fDCAyEtaPhiSecMCDecaysK0s(0),
  fDCAzEtaPhiSecMCDecaysK0s(0),
  //
  fDCAyEtaPtSecMCDecaysLambda(0),
  fDCAzEtaPtSecMCDecaysLambda(0),
  fDCAyPhiPtSecMCDecaysLambda(0),
  fDCAzPhiPtSecMCDecaysLambda(0),
  fDCAyEtaPhiSecMCDecaysLambda(0),
  fDCAzEtaPhiSecMCDecaysLambda(0),
  //
  fDCAyEtaPtSecMCMaterial(0),
  fDCAzEtaPtSecMCMaterial(0),
  fDCAyPhiPtSecMCMaterial(0),
  fDCAzPhiPtSecMCMaterial(0),
  fDCAyEtaPhiSecMCMaterial(0),
  fDCAzEtaPhiSecMCMaterial(0)
{
  // constructor
  Init();
}

//_____________________________________________________________________________
AlidNdPtCutAnalysis::~AlidNdPtCutAnalysis() {
    //
    if(fEventCount) delete fEventCount; fEventCount=0;
    if(fRecEventHist) delete fRecEventHist; fRecEventHist=0;
    if(fMCEventHist) delete fMCEventHist; fMCEventHist=0;
    if(fRecMCEventHist) delete fRecMCEventHist; fRecMCEventHist=0;
    if(fRecMCTrackHist) delete fRecMCTrackHist; fRecMCTrackHist=0;

    if(fTrigVertCount) delete fTrigVertCount;
    if(fRecEventXYZ) delete fRecEventXYZ;
    if(fRecEventXYMult) delete fRecEventXYMult ;
    if(fRecEventXZMult) delete fRecEventXZMult;
    if(fRecEventYZMult) delete fRecEventYZMult;
    if(fRecEventZResZMult) delete fRecEventZResZMult;
    if(fMCEventXYZ) delete fMCEventXYZ;
    if(fRecMCEventDXDYDZ) delete fRecMCEventDXDYDZ;
    if(fRecMCEventDXDYMult) delete fRecMCEventDXDYMult;
    if(fRecMCEventDXDZMult) delete fRecMCEventDXDZMult;
    if(fRecMCEventDYDZMult) delete fRecMCEventDYDZMult;

    if(fCrossRowsEtaPt) delete fCrossRowsEtaPt;
    if(fChi2PerClustEtaPt) delete fChi2PerClustEtaPt;
    if(fCrossRowsOverFindEtaPt) delete fCrossRowsOverFindEtaPt;
    if(fFracSharedClustEtaPt) delete fFracSharedClustEtaPt;
    if(fDCAyEtaPt) delete fDCAyEtaPt;
    if(fDCAzEtaPt) delete fDCAzEtaPt;

    if(fCrossRowsPhiPt) delete fCrossRowsPhiPt;
    if(fChi2PerClustPhiPt) delete fChi2PerClustPhiPt;
    if(fCrossRowsOverFindPhiPt) delete fCrossRowsOverFindPhiPt;
    if(fFracSharedClustPhiPt) delete fFracSharedClustPhiPt;
    if(fDCAyPhiPt) delete fDCAyPhiPt;
    if(fDCAzPhiPt) delete fDCAzPhiPt;

    if(fCrossRowsEtaPhi) delete fCrossRowsEtaPhi;
    if(fChi2PerClustEtaPhi) delete fChi2PerClustEtaPhi;
    if(fCrossRowsOverFindEtaPhi) delete fCrossRowsOverFindEtaPhi;
    if(fFracSharedClustEtaPhi) delete fFracSharedClustEtaPhi;
    if(fDCAyEtaPhi) delete fDCAyEtaPhi;
    if(fDCAzEtaPhi) delete fDCAzEtaPhi;

    if(fDCAyEtaPtMCPrim) delete fDCAyEtaPtMCPrim;
    if(fDCAzEtaPtMCPrim) delete fDCAzEtaPtMCPrim;
    if(fDCAyPhiPtMCPrim) delete fDCAyPhiPtMCPrim;
    if(fDCAzPhiPtMCPrim) delete fDCAzPhiPtMCPrim;
    if(fDCAyEtaPhiMCPrim) delete fDCAyEtaPhiMCPrim;
    if(fDCAzEtaPhiMCPrim) delete fDCAzEtaPhiMCPrim;

    if(fDCAyEtaPtMCSec) delete fDCAyEtaPtMCSec;
    if(fDCAzEtaPtMCSec) delete fDCAzEtaPtMCSec;
    if(fDCAyPhiPtMCSec) delete fDCAyPhiPtMCSec;
    if(fDCAzPhiPtMCSec) delete fDCAzPhiPtMCSec;
    if(fDCAyEtaPhiMCSec) delete fDCAyEtaPhiMCSec;
    if(fDCAzEtaPhiMCSec) delete fDCAzEtaPhiMCSec;

    if(fDCAyEtaPtSecMCDecays) delete fDCAyEtaPtSecMCDecays;
    if(fDCAzEtaPtSecMCDecays) delete fDCAzEtaPtSecMCDecays;
    if(fDCAyPhiPtSecMCDecays) delete fDCAyPhiPtSecMCDecays;
    if(fDCAzPhiPtSecMCDecays) delete fDCAzPhiPtSecMCDecays;
    if(fDCAyEtaPhiSecMCDecays) delete fDCAyEtaPhiSecMCDecays;
    if(fDCAzEtaPhiSecMCDecays) delete fDCAzEtaPhiSecMCDecays;

    if(fDCAyEtaPtSecMCDecaysK0s) delete fDCAyEtaPtSecMCDecaysK0s;
    if(fDCAzEtaPtSecMCDecaysK0s) delete fDCAzEtaPtSecMCDecaysK0s;
    if(fDCAyPhiPtSecMCDecaysK0s) delete fDCAyPhiPtSecMCDecaysK0s;
    if(fDCAzPhiPtSecMCDecaysK0s) delete fDCAzPhiPtSecMCDecaysK0s;
    if(fDCAyEtaPhiSecMCDecaysK0s) delete fDCAyEtaPhiSecMCDecaysK0s;
    if(fDCAzEtaPhiSecMCDecaysK0s) delete fDCAzEtaPhiSecMCDecaysK0s;

    if(fDCAyEtaPtSecMCDecaysLambda) delete fDCAyEtaPtSecMCDecaysLambda;
    if(fDCAzEtaPtSecMCDecaysLambda) delete fDCAzEtaPtSecMCDecaysLambda;
    if(fDCAyPhiPtSecMCDecaysLambda) delete fDCAyPhiPtSecMCDecaysLambda;
    if(fDCAzPhiPtSecMCDecaysLambda) delete fDCAzPhiPtSecMCDecaysLambda;
    if(fDCAyEtaPhiSecMCDecaysLambda) delete fDCAyEtaPhiSecMCDecaysLambda;
    if(fDCAzEtaPhiSecMCDecaysLambda) delete fDCAzEtaPhiSecMCDecaysLambda;

    if(fDCAyEtaPtSecMCMaterial) delete fDCAyEtaPtSecMCMaterial;
    if(fDCAzEtaPtSecMCMaterial) delete fDCAzEtaPtSecMCMaterial;
    if(fDCAyPhiPtSecMCMaterial) delete fDCAyPhiPtSecMCMaterial;
    if(fDCAzPhiPtSecMCMaterial) delete fDCAzPhiPtSecMCMaterial;
    if(fDCAyEtaPhiSecMCMaterial) delete fDCAyEtaPhiSecMCMaterial;
    if(fDCAzEtaPhiSecMCMaterial) delete fDCAzEtaPhiSecMCMaterial;

    if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AlidNdPtCutAnalysis::Init(){
  //
  // Init histograms
  //
  /*
  const Int_t ptNbins = 56;
  const Double_t ptMin = 0.;
  const Double_t ptMax = 16.;
  Double_t binsPt[ptNbins+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};
  */
  // set pt bins
  const Int_t ptNbins = 50;
  const Double_t ptMin = 1.e-2, ptMax = 50.;
  Double_t *binsPt = CreateLogAxis(ptNbins,ptMin,ptMax);

  //
  Int_t binsEventCount[2]={2,2};
  Double_t minEventCount[2]={0,0};
  Double_t maxEventCount[2]={2,2};
  fEventCount = new THnSparseF("fEventCount","trig vs trig+vertex",2,binsEventCount,minEventCount,maxEventCount);
  fEventCount->GetAxis(0)->SetTitle("trig");
  fEventCount->GetAxis(1)->SetTitle("trig+vert");
  fEventCount->Sumw2();

  fTrigVertCount = new TH2D("fTrigVertCount","trig:trigger&&vertex; trigger; trigger and vertex",2,0,2,2,0,2);

  //Xv:Yv:Zv:ResZv:Mult
  Double_t kFact = 1.0;

  Int_t binsRecEventHist[5]={80,80,100,80,150};
  Double_t minRecEventHist[5]={-3.*kFact,-3.*kFact,-35.,0.,0.};
  Double_t maxRecEventHist[5]={3.*kFact,3.*kFact,35.,10.,150.};
  fRecEventHist = new THnSparseF("fRecEventHist","Xv:Yv:Zv:ResZv:Mult",5,binsRecEventHist,minRecEventHist,maxRecEventHist);
  fRecEventHist->GetAxis(0)->SetTitle("Xv (cm)");
  fRecEventHist->GetAxis(1)->SetTitle("Yv (cm)");
  fRecEventHist->GetAxis(2)->SetTitle("Zv (cm)");
  fRecEventHist->GetAxis(3)->SetTitle("ResZv (cm)");
  fRecEventHist->GetAxis(4)->SetTitle("Mult");
  fRecEventHist->Sumw2();

  fRecEventXYZ = new TH3D("fRecEventXYZ","Xv:Yv:Zv; Xv (cm); Yv (cm); Zv (cm)", 80,-3.*kFact,3.*kFact,80,-3.*kFact,3.*kFact,100, -35,35);
  fRecEventXYMult = new TH3D("fRecEventXYMult","Xv:Yv:Mult; Xv (cm); Yv (cm); Mult", 80,-3.*kFact,3.*kFact,80,-3.*kFact,3.*kFact,150,-0.5,149.5);
  fRecEventXZMult = new TH3D("fRecEventXZMult","Xv:Zv:Mult; Xv (cm); Zv (cm); Mult", 80,-3.*kFact,3.*kFact,100,-35,35,150,-0.5,149.5);
  fRecEventYZMult = new TH3D("fRecEventYZMult","Yv:Zv:Mult; Yv (cm); Zv (cm); Mult", 80,-3.*kFact,3.*kFact,100,-35,35,150,-0.5,149.5);
  fRecEventZResZMult = new TH3D("fRecEventZResZMult","Zv:ZRes:Mult; Zv (cm); ZRes (cm); Mult",100,-35,35,80,0,10,150,-0.5,149.5);


  //Xv:Yv:Zv
  Int_t binsMCEventHist[3]={80,80,100};
  Double_t minMCEventHist[3]={-0.1,-0.1,-35.};
  Double_t maxMCEventHist[3]={0.1,0.1,35.};
  fMCEventHist = new THnSparseF("fMCEventHist","mcXv:mcYv:mcZv",3,binsMCEventHist,minMCEventHist,maxMCEventHist);
  fMCEventHist->GetAxis(0)->SetTitle("mcXv (cm)");
  fMCEventHist->GetAxis(1)->SetTitle("mcYv (cm)");
  fMCEventHist->GetAxis(2)->SetTitle("mcZv (cm)");
  fMCEventHist->Sumw2();

  fMCEventXYZ = new TH3D("fMCEventXYZ","mcXv:mcYv:mcZv; mcXv (cm); mcYv (cm); mcZv (cm)", 80,-0.1,0.1,80,-0.1,0.1,100, -35,35);


  //Xv-mcXv:Yv-mcYv:Zv-mcZv:Mult
  Int_t binsRecMCEventHist[4]={100,100,100,150};
  Double_t minRecMCEventHist[4]={-1.0*kFact,-1.0*kFact,-1.0*kFact,0.};
  Double_t maxRecMCEventHist[4]={1.0*kFact,1.0*kFact,1.0*kFact,150.};
  fRecMCEventHist = new THnSparseF("fRecMCEventHist","Xv-mcXv:Yv-mcYv:Zv-mcZv:Mult",4,binsRecMCEventHist,minRecMCEventHist,maxRecMCEventHist);
  fRecMCEventHist->GetAxis(0)->SetTitle("Xv-mcXv (cm)");
  fRecMCEventHist->GetAxis(1)->SetTitle("Yv-mcYv (cm)");
  fRecMCEventHist->GetAxis(2)->SetTitle("Zv-mcZv (cm)");
  fRecMCEventHist->GetAxis(3)->SetTitle("Mult");
  fRecMCEventHist->Sumw2();


  fRecMCEventDXDYDZ = new TH3D("fRecMCEventDXDYDZ","Xv-mcXv:Yv-mcYv:Zv-mcZv; Xv-mcXv (cm); Yv-mcYv (cm); Zv-mcZv (cm)",80,-3.*kFact,3.*kFact,80,-3.*kFact,3.*kFact,100, -35,35);
  fRecMCEventDXDYMult = new TH3D("fRecMCEventDXDYMult","Xv-mcXv:Yv-mcYv:Mult; Xv-mcXv (cm); Yv-mcYv (cm); Mult",80,-3.*kFact,3.*kFact,80,-3.*kFact,3.*kFact,150,-0.5,149.5);
  fRecMCEventDXDZMult = new TH3D("fRecMCEventDXDZMult","Xv-mcXv:Zv-mcZv:Mult; Xv-mcXv (cm); Zv-mcZv (cm); Mult",80,-3.*kFact,3.*kFact,80,-3.*kFact,3.*kFact,150,-0.5,149.5);
  fRecMCEventDYDZMult = new TH3D("fRecMCEventDYDZMult","Yv-mcYv:Zv-mcZv:Mult; Yv-mcYv (cm); Zv-mcZv (cm); Mult",80,-3.*kFact,3.*kFact,80,-3.*kFact,3.*kFact,150,-0.5,149.5);

  //
  // THnSparse track histograms
  //

 //nCrossRows:chi2PerClust:nCrossRows/nFindableClust:fracSharedClust:DCAy:DCAz:eta:phi:pt:hasStrangeMother:isFromMaterial:isPrim:charge
  Int_t binsRecMCTrackHist[13]=  { 160,  10,  20,  20, 50,  50,   20,  90,             ptNbins, 2,  2,  2,  3  };
  Double_t minRecMCTrackHist[13]={ 0.,   0.,  0.,  0., -0.5,-0.5,-1.0, 0.,             ptMin,   0., 0., 0.,-1. };
  Double_t maxRecMCTrackHist[13]={ 160., 10., 1.,  1., 0.5, 0.5,  1.0, 2.*TMath::Pi(), ptMax,   2., 2., 2., 2. };

  fRecMCTrackHist = new THnSparseF("fRecMCTrackHist","nCrossRows:chi2PerClust:nCrossRows/nFindableClust:fracSharedClust:DCAy:DCAz:eta:phi:pt:hasStrangeMother:isFromMaterial:isPrim:charge",13,binsRecMCTrackHist,minRecMCTrackHist,maxRecMCTrackHist);
  fRecMCTrackHist->SetBinEdges(8,binsPt);
  fRecMCTrackHist->GetAxis(0)->SetTitle("nCrossRows");
  fRecMCTrackHist->GetAxis(1)->SetTitle("chi2PerClust");
  fRecMCTrackHist->GetAxis(2)->SetTitle("nCrossRows/nFindableClust");
  fRecMCTrackHist->GetAxis(3)->SetTitle("fracSharedClust");
  fRecMCTrackHist->GetAxis(4)->SetTitle("DCAy (cm)");
  fRecMCTrackHist->GetAxis(5)->SetTitle("DCAz (cm)");
  fRecMCTrackHist->GetAxis(6)->SetTitle("#eta");
  fRecMCTrackHist->GetAxis(7)->SetTitle("#phi (rad)");
  fRecMCTrackHist->GetAxis(8)->SetTitle("p_{T} (GeV/c)");
  fRecMCTrackHist->GetAxis(9)->SetTitle("hasStrangeMother");
  fRecMCTrackHist->GetAxis(10)->SetTitle("isFromMaterial");
  fRecMCTrackHist->GetAxis(11)->SetTitle("isPrim");
  fRecMCTrackHist->GetAxis(12)->SetTitle("charge");
  fRecMCTrackHist->Sumw2();



  //
  fCrossRowsEtaPt = new TH3D("fCrossRowsEtaPt","nCrossRows:eta:pt; nCrossRows; eta; pt (GeV/c)", 160,0,160,20,-1,1,ptNbins,ptMin,ptMax);
  fChi2PerClustEtaPt = new TH3D("fChi2PerClustEtaPt","chi2PerClust:eta:pt; chi2PerClust; eta; pt (GeV/c)", 100,0,10,20,-1,1,ptNbins,ptMin,ptMax);
  fCrossRowsOverFindEtaPt = new TH3D("fCrossRowsOverFindEtaPt","nCrossRows/nFindableClust:eta:pt; nCrossRows/nFindableClust; eta; pt (GeV/c)", 100,0,2,20,-1,1,ptNbins,ptMin,ptMax);
  fFracSharedClustEtaPt = new TH3D("fFracSharedClustEtaPt","fracSharedClust:eta:pt;fracSharedClust; eta; pt (GeV/c)", 100,0,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyEtaPt = new TH3D("fDCAyEtaPt","DCAy:eta:pt;DCAy (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPt = new TH3D("fDCAzEtaPt","DCAz:eta:pt;DCAz (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);

  //
  fCrossRowsPhiPt = new TH3D("fCrossRowsPhiPt","nCrossRows:eta:pt; nCrossRows; phi (rad); pt (GeV/c)", 160,0,160,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fChi2PerClustPhiPt = new TH3D("fChi2PerClustPhiPt","chi2PerClust:eta:pt; chi2PerClust; phi (rad); pt (GeV/c)", 100,0,10,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fCrossRowsOverFindPhiPt = new TH3D("fCrossRowsOverFindPhiPt","nCrossRows/nFindableClust:eta:pt; nCrossRows/nFindableClust; phi (rad); pt (GeV/c)", 100,0,2,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fFracSharedClustPhiPt = new TH3D("fFracSharedClustPhiPt","fracSharedClust:eta:pt;fracSharedClust; phi (rad); pt (GeV/c)", 100,0,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyPhiPt = new TH3D("fDCAyPhiPt","DCAy:eta:pt;DCAy (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPt = new TH3D("fDCAzPhiPt","DCAz:eta:pt;DCAz (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);

  //
  fCrossRowsEtaPhi = new TH3D("fCrossRowsEtaPhi","nCrossRows:eta:phi; nCrossRows; eta; phi (rad)", 160,0,160,20,-1,1,90,0,2.*TMath::Pi());
  fChi2PerClustEtaPhi = new TH3D("fChi2PerClustEtaPhi","chi2PerClust:eta:phi; chi2PerClust; eta; phi (rad)", 100,0,10,20,-1,1,90,0,2.*TMath::Pi());
  fCrossRowsOverFindEtaPhi = new TH3D("fCrossRowsOverFindEtaPhi","nCrossRows/nFindableClust:eta:phi; nCrossRows/nFindableClust; eta; phi (rad)", 100,0,2,20,-1,1,90,0,2.*TMath::Pi());
  fFracSharedClustEtaPhi = new TH3D("fFracSharedClustEtaPhi","fracSharedClust:eta:phi;fracSharedClust; eta; phi (rad)", 100,0,1,20,-1,1,90,0,2.*TMath::Pi());
  fDCAyEtaPhi = new TH3D("fDCAyEtaPhi","DCAy:eta:phi;DCAy (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhi = new TH3D("fDCAzEtaPhi","DCAz:eta:phi;DCAz (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());

  //
  fDCAyEtaPtMCPrim = new TH3D("fDCAyEtaPtMCPrim","DCAy:eta:pt primary;DCAy (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtMCPrim = new TH3D("fDCAzEtaPtMCPrim","DCAz:eta:pt primary;DCAz (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtMCPrim = new TH3D("fDCAyPhiPtMCPrim","DCAy primary:eta:pt primary;DCAy (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtMCPrim = new TH3D("fDCAzPhiPtMCPrim","DCAz:eta:pt primary;DCAz (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiMCPrim = new TH3D("fDCAyEtaPhiMCPrim","DCAy:eta:phi primary;DCAy (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiMCPrim = new TH3D("fDCAzEtaPhiMCPrim","DCAz:eta:phi primary;DCAz (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());

  //
  fDCAyEtaPtMCSec = new TH3D("fDCAyEtaPtMCSec","DCAy:eta:pt secondary;DCAy (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtMCSec = new TH3D("fDCAzEtaPtMCSec","DCAz:eta:pt secondary;DCAz (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtMCSec = new TH3D("fDCAyPhiPtMCSec","DCAy secondary:eta:pt secondary;DCAy (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtMCSec = new TH3D("fDCAzPhiPtMCSec","DCAz:eta:pt secondary;DCAz (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiMCSec = new TH3D("fDCAyEtaPhiMCSec","DCAy:eta:phi secondary;DCAy (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiMCSec = new TH3D("fDCAzEtaPhiMCSec","DCAz:eta:phi secondary;DCAz (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());

  fDCAyEtaPtSecMCDecays = new TH3D("fDCAyEtaPtSecMCDecays","DCAy:eta:pt secondary decays;DCAy (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtSecMCDecays = new TH3D("fDCAzEtaPtSecMCDecays","DCAz:eta:pt secondary decays;DCAz (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtSecMCDecays = new TH3D("fDCAyPhiPtSecMCDecays","DCAy secondary decays:eta:pt secondary decays;DCAy (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtSecMCDecays = new TH3D("fDCAzPhiPtSecMCDecays","DCAz:eta:pt secondary decays;DCAz (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiSecMCDecays = new TH3D("fDCAyEtaPhiSecMCDecays","DCAy:eta:phi secondary decays;DCAy (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiSecMCDecays = new TH3D("fDCAzEtaPhiSecMCDecays","DCAz:eta:phi secondary decays;DCAz (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());

  fDCAyEtaPtSecMCDecaysK0s = new TH3D("fDCAyEtaPtSecMCDecaysK0s","DCAy:eta:pt secondary decays from K0s;DCAy (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtSecMCDecaysK0s = new TH3D("fDCAzEtaPtSecMCDecaysK0s","DCAz:eta:pt secondary decays from K0s;DCAz (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtSecMCDecaysK0s = new TH3D("fDCAyPhiPtSecMCDecaysK0s","DCAy secondary decays from K0s:eta:pt secondary decays from K0s;DCAy (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtSecMCDecaysK0s = new TH3D("fDCAzPhiPtSecMCDecaysK0s","DCAz:eta:pt secondary decays from K0s;DCAz (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiSecMCDecaysK0s = new TH3D("fDCAyEtaPhiSecMCDecaysK0s","DCAy:eta:phi secondary decays from K0s;DCAy (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiSecMCDecaysK0s = new TH3D("fDCAzEtaPhiSecMCDecaysK0s","DCAz:eta:phi secondary decays from K0s;DCAz (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());

  fDCAyEtaPtSecMCDecaysLambda = new TH3D("fDCAyEtaPtSecMCDecaysLambda","DCAy:eta:pt secondary decays from Lambda;DCAy (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtSecMCDecaysLambda = new TH3D("fDCAzEtaPtSecMCDecaysLambda","DCAz:eta:pt secondary decays from Lambda;DCAz (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtSecMCDecaysLambda = new TH3D("fDCAyPhiPtSecMCDecaysLambda","DCAy secondary decays from Lambda:eta:pt secondary decays from Lambda;DCAy (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtSecMCDecaysLambda = new TH3D("fDCAzPhiPtSecMCDecaysLambda","DCAz:eta:pt secondary decays from Lambda;DCAz (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiSecMCDecaysLambda = new TH3D("fDCAyEtaPhiSecMCDecaysLambda","DCAy:eta:phi secondary decays from Lambda;DCAy (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiSecMCDecaysLambda = new TH3D("fDCAzEtaPhiSecMCDecaysLambda","DCAz:eta:phi secondary decays from Lambda;DCAz (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());

  fDCAyEtaPtSecMCMaterial = new TH3D("fDCAyEtaPtSecMCMaterial","DCAy:eta:pt secondary material;DCAy (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtSecMCMaterial = new TH3D("fDCAzEtaPtSecMCMaterial","DCAz:eta:pt secondary material;DCAz (cm); eta; pt (GeV/c)", 100,-1,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtSecMCMaterial = new TH3D("fDCAyPhiPtSecMCMaterial","DCAy secondary material:eta:pt secondary material;DCAy (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtSecMCMaterial = new TH3D("fDCAzPhiPtSecMCMaterial","DCAz:eta:pt secondary material;DCAz (cm); phi (rad); pt (GeV/c)", 100,-1,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiSecMCMaterial = new TH3D("fDCAyEtaPhiSecMCMaterial","DCAy:eta:phi secondary material;DCAy (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiSecMCMaterial = new TH3D("fDCAzEtaPhiSecMCMaterial","DCAz:eta:phi secondary material;DCAz (cm); eta; phi (rad)", 100,-1,1,20,-1,1,90,0,2.*TMath::Pi());

  //nClust:chi2PerClust:nClust/nFindableClust:DCAy:DCAz:eta:phi:pt:kinkIdx:isPrim:polarity
  /*
  Int_t binsRecMCTrackHist[11]={160,80,80,100,100,90,90,ptNbins, 3, 2, 2};
  Double_t minRecMCTrackHist[11]={0., 0., 0., -1.,-1.,-1.5, 0., ptMin, -1., 0., 0.};
  Double_t maxRecMCTrackHist[11]={160.,10.,1.2, 1.,1.,1.5, 2.*TMath::Pi(), ptMax, 2., 2., 2.};

  fRecMCTrackHist = new THnSparseF("fRecMCTrackHist","nClust:chi2PerClust:nClust/nFindableClust:DCAy:DCAz:eta:phi:pt:kinkIdx:isPrim:polarity",11,binsRecMCTrackHist,minRecMCTrackHist,maxRecMCTrackHist);
  fRecMCTrackHist->SetBinEdges(7,binsPt);

  fRecMCTrackHist->GetAxis(0)->SetTitle("nClust");
  fRecMCTrackHist->GetAxis(1)->SetTitle("chi2PerClust");
  fRecMCTrackHist->GetAxis(2)->SetTitle("nClust/nFindableClust");
  fRecMCTrackHist->GetAxis(3)->SetTitle("DCAy (cm)");
  fRecMCTrackHist->GetAxis(4)->SetTitle("DCAz (cm)");
  fRecMCTrackHist->GetAxis(5)->SetTitle("#eta");
  fRecMCTrackHist->GetAxis(6)->SetTitle("#phi (rad)");
  fRecMCTrackHist->GetAxis(7)->SetTitle("p_{T} (GeV/c)");
  fRecMCTrackHist->GetAxis(8)->SetTitle("kinkIdx"); // 0 - no kink, -1 - kink mother, 1 - kink daugther
  fRecMCTrackHist->GetAxis(9)->SetTitle("isPrim");
  fRecMCTrackHist->GetAxis(10)->SetTitle("polarity");
  fRecMCTrackHist->Sumw2();
  */

  // init output folder
  fAnalysisFolder = CreateFolder("folderdNdPt","Analysis dNdPt Folder");

}

//_____________________________________________________________________________
void AlidNdPtCutAnalysis::Process(AliESDEvent *const esdEvent, AliMCEvent * const mcEvent)
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
    isEventTriggered = inputHandler->IsEventSelected() &  GetTriggerMask();

    physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) return;

    if (isEventTriggered && (GetTrigger() == AliTriggerAnalysis::kV0AND)) {
      // set trigger (V0AND)
      triggerAnalysis = physicsSelection->GetTriggerAnalysis();
      if(!triggerAnalysis) return;
      isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(esdEvent, GetTrigger());
    }
  }

  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);
  AliPWG0Helper::MCProcessType evtType = AliPWG0Helper::kInvalidProcess;

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
    AliDebug(AliLog::kDebug+1, Form("Found process type %d", evtType));

    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);

    // Fill MC event histogram
    Double_t vMCEventHist[3]={vtxMC[0],vtxMC[1],vtxMC[2]};
    if(fFillSparseHisto) { fMCEventHist->Fill(vMCEventHist); }

    fMCEventXYZ->Fill(vMCEventHist[0], vMCEventHist[1] ,vMCEventHist[2]);


  } // end bUseMC

  // get reconstructed vertex
  /*
  Bool_t bRedoTPCVertex = evtCuts->IsRedoTPCVertex();
  Bool_t bUseConstraints = evtCuts->IsUseBeamSpotConstraint();
  const AliESDVertex* vtxESD = AlidNdPtHelper::GetVertex(esdEvent,evtCuts,accCuts,esdTrackCuts,GetAnalysisMode(),kFALSE,bRedoTPCVertex,bUseConstraints);
  Bool_t isRecVertex = AlidNdPtHelper::TestRecVertex(vtxESD, esdEvent->GetPrimaryVertexSPD(), GetAnalysisMode(), kFALSE);
  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD) && isRecVertex;
  */

  const AliESDVertex* vtxESD = 0;
  Bool_t isRecVertex = kFALSE;

  vtxESD = esdEvent->GetPrimaryVertexTracks();
  if(!vtxESD || (vtxESD->GetNContributors()<1)) {
  // SPD vertex
    vtxESD = esdEvent->GetPrimaryVertexSPD();
    if (vtxESD->GetNContributors()>0) {
      if ( !vtxESD->IsFromVertexerZ() || (vtxESD->GetDispersion()<0.04 && vtxESD->GetZRes()<0.25)) isRecVertex = kTRUE;
    }
    } else {
       isRecVertex = kTRUE;
    }

  Bool_t isEventOK = evtCuts->AcceptEvent(esdEvent,mcEvent,vtxESD) && isRecVertex;


  TObjArray *allChargedTracks=0;
  Int_t multAll=0;

  //
  // event counter
  //
  //printf("isEventOK %d, isEventTriggered %d \n",isEventOK,isEventTriggered);
  Bool_t isTrigAndVertex = isEventTriggered && isEventOK;

  Double_t vEventCount[2] = { static_cast<Double_t>(isEventTriggered), static_cast<Double_t>(isTrigAndVertex)};
  if(fFillSparseHisto) {fEventCount->Fill(vEventCount);}

  fTrigVertCount->Fill(vEventCount[0],vEventCount[1]);


  //
  // cosmic background and splitted tracks
  //
  if(GetParticleMode() == AlidNdPtHelper::kBackgroundTrack)
  {
    AlidNdPtBackgroundCuts *backCuts = GetBackgroundCuts();
    if(!backCuts) return;

    for (Int_t iTrack = 0; iTrack < esdEvent->GetNumberOfTracks(); iTrack++)
    {
      AliESDtrack *track1 = esdEvent->GetTrack(iTrack);
      if(!track1) continue;
      if(track1->Charge()==0) continue;

      for (Int_t jTrack = iTrack+1; jTrack < esdEvent->GetNumberOfTracks(); jTrack++)
      {
        AliESDtrack *track2 = esdEvent->GetTrack(jTrack);
        if(!track2) continue;
        if(track2->Charge()==0) continue;

	//printf("track2->Charge() %d\n",track2->Charge());

        backCuts->IsBackgroundTrack(track1, track2);
      }
    }
  }

  // check event cuts
  if(isEventOK && isEventTriggered)
  {
    // get all charged tracks
    allChargedTracks = AlidNdPtHelper::GetAllChargedTracks(esdEvent,GetAnalysisMode());
    if(!allChargedTracks) return;

    Int_t entries = allChargedTracks->GetEntries();
    for(Int_t i=0; i<entries;++i)
    {
      AliESDtrack *track = (AliESDtrack*)allChargedTracks->At(i);
      if(!track) continue;

      if(!accCuts->AcceptTrack(track)) continue;
      if(!esdTrackCuts->AcceptTrack(track)) continue;

      //
      Bool_t isOK = kFALSE;
      Double_t x[3]; track->GetXYZ(x);
      Double_t b[3]; AliTracker::GetBxByBz(x,b);

      //
      // if TPC-ITS hybrid tracking (kTPCITSHybrid)
      // replace track parameters with TPC-ony track parameters
      //
      if( GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybrid || GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtx || GetAnalysisMode() == AlidNdPtHelper::kTPCITSHybridTrackSPDvtxDCArPt)
      {
        // Relate TPC-only tracks to SPD vertex
        isOK = track->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig);
        if(!isOK) continue;

	// replace esd track parameters with TPCinner
        AliExternalTrackParam  *tpcTrack  = new AliExternalTrackParam(*(track->GetTPCInnerParam()));
	if (!tpcTrack) return;
        track->Set(tpcTrack->GetX(),tpcTrack->GetAlpha(),tpcTrack->GetParameter(),tpcTrack->GetCovariance());

        if(tpcTrack) delete tpcTrack;
      }

      //
      if (GetAnalysisMode()==AlidNdPtHelper::kTPCSPDvtxUpdate || GetAnalysisMode() == AlidNdPtHelper::kTPCTrackSPDvtxUpdate)
      {
        //
        // update track parameters
	//
        AliExternalTrackParam cParam;
        isOK = track->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig, &cParam);
	if(!isOK) continue;

	track->Set(cParam.GetX(),cParam.GetAlpha(),cParam.GetParameter(),cParam.GetCovariance());
      }

      FillHistograms(track, stack);
      multAll++;
    }

    Double_t vRecEventHist[5] = {vtxESD->GetX(),vtxESD->GetY(),vtxESD->GetZ(),vtxESD->GetZRes(),static_cast<Double_t>(multAll)};
    if(fFillSparseHisto) { fRecEventHist->Fill(vRecEventHist); }

    fRecEventXYZ->Fill(vtxESD->GetX(),vtxESD->GetY(),vtxESD->GetZ());
    fRecEventXYMult->Fill(vtxESD->GetX(),vtxESD->GetY(),static_cast<Double_t>(multAll));
    fRecEventXZMult->Fill(vtxESD->GetX(),vtxESD->GetZ(),static_cast<Double_t>(multAll));
    fRecEventYZMult->Fill(vtxESD->GetY(),vtxESD->GetZ(),static_cast<Double_t>(multAll));
    fRecEventZResZMult->Fill(vtxESD->GetZ(),vtxESD->GetZRes(),static_cast<Double_t>(multAll));

    if(IsUseMCInfo()) {
      Double_t vRecMCEventHist[5] = {vtxESD->GetX()-vtxMC[0],vtxESD->GetY()-vtxMC[1],vtxESD->GetZ()-vtxMC[2],static_cast<Double_t>(multAll)};
      if(fFillSparseHisto) { fRecMCEventHist->Fill(vRecMCEventHist); }

      fRecMCEventDXDYDZ->Fill(vtxESD->GetX()-vtxMC[0],vtxESD->GetY()-vtxMC[1],vtxESD->GetZ()-vtxMC[2]);;
      fRecMCEventDXDYMult->Fill(vtxESD->GetX()-vtxMC[0],vtxESD->GetY()-vtxMC[1],static_cast<Double_t>(multAll));
      fRecMCEventDXDZMult->Fill(vtxESD->GetX()-vtxMC[0],vtxESD->GetZ()-vtxMC[2],static_cast<Double_t>(multAll));
      fRecMCEventDYDZMult->Fill(vtxESD->GetY()-vtxMC[1],vtxESD->GetZ()-vtxMC[2],static_cast<Double_t>(multAll));
    }
  }

  if(allChargedTracks) delete allChargedTracks; allChargedTracks = 0;

}

//_____________________________________________________________________________
void AlidNdPtCutAnalysis::FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack) const
{
  //
  // Fill ESD track and MC histograms
  //
  if(!esdTrack) return;
  if(esdTrack->Charge() == 0.) return;

  Float_t pt = esdTrack->Pt();
  Float_t eta = esdTrack->Eta();
  Float_t phi = esdTrack->Phi();


  Int_t nClust = 0;
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
    nClust = esdTrack->GetTPCNclsIter1();
  } else {
    nClust = esdTrack->GetTPCclusters(0);
  }

  Float_t chi2PerCluster = 0.;
  if(GetAnalysisMode() == AlidNdPtHelper::kTPC) {
    if(nClust>0.) chi2PerCluster = esdTrack->GetTPCchi2Iter1()/Float_t(nClust);
  } else {
    chi2PerCluster = esdTrack->GetTPCchi2()/Float_t(nClust);
  }

  Int_t nFindableClust = esdTrack->GetTPCNclsF();


  Float_t clustPerFindClust = 0.;
  if(nFindableClust>0.) clustPerFindClust = Float_t(nClust)/nFindableClust;

  Float_t b[2], bCov[3];
  esdTrack->GetImpactParameters(b,bCov);

  //
  Float_t nCrossedRowsTPC = esdTrack->GetTPCClusterInfo(2,1);

  Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
  if (esdTrack->GetTPCNclsF()>0) {
     ratioCrossedRowsOverFindableClustersTPC = esdTrack->GetTPCClusterInfo(2,1)/esdTrack->GetTPCNclsF();
  }

  //
  Int_t nClustersTPCShared = esdTrack->GetTPCnclsS();
  Float_t fracClustersTPCShared = -1.;
  fracClustersTPCShared = Float_t(nClustersTPCShared)/Float_t(nClust);


  // kink idx
  Int_t kinkIdx = 0;
  //if(esdTrack->GetKinkIndex(0) > 0.)   isKink  = kTRUE;
  if(esdTrack->GetKinkIndex(0) > 0)      kinkIdx = 1;   // kink daughter
  else if(esdTrack->GetKinkIndex(0) < 0) kinkIdx = -1;  // kink mother
  else kinkIdx = 0; // not kink

  //printf("esdTrack->GetKinkIndex(0) %d \n", esdTrack->GetKinkIndex(0));
  //printf("esdTrack->GetKinkIndex(1) %d \n", esdTrack->GetKinkIndex(1));
  //printf("esdTrack->GetKinkIndex(2) %d \n", esdTrack->GetKinkIndex(2));
  //printf("kinkIdx %d \n", kinkIdx);

  //
  // Fill rec vs MC information
  //
  Bool_t isPrim = kFALSE;
  Bool_t hasStrangeMother = kFALSE;
  Bool_t isFromMaterial = kFALSE;
  Bool_t isFromK0s = kFALSE;
  Bool_t isFromLambda = kFALSE;

  if(IsUseMCInfo()) {
    if(!stack) return;
    Int_t label = TMath::Abs(esdTrack->GetLabel());
    TParticle* particle = stack->Particle(label);
    if(!particle) return;
    if(particle->GetPDG() && particle->GetPDG()->Charge()==0.) return;
    isPrim = stack->IsPhysicalPrimary(label);

    // check whether has stange mother
    //
    Int_t motherPdg = -1;
    TParticle* mother = 0;

    Int_t motherLabel = particle->GetMother(0);
    if(motherLabel>=0) mother = stack->Particle(motherLabel);
    if(mother) motherPdg = TMath::Abs(mother->GetPdgCode()); // take abs for visualisation only
    Int_t mech = particle->GetUniqueID(); // production mechanism


    if(motherPdg==3122 || motherPdg==310) // lambda, antilambda, k0s
    {
      if(mech == 4 || mech == 5) {
          hasStrangeMother = kTRUE;
          if(motherPdg == 310) isFromK0s = kTRUE;
          if(motherPdg == 3122) isFromLambda = kTRUE;
      }
    }
    else {
      //if(isPrim==0 && mech == 13)
      //printf("mech %d \n", mech);
      if(!isPrim) isFromMaterial = kTRUE;
    }
    //if(isPrim && pt > 1.5 && kinkIdx == -1) printf("nClust  %d \n", nClust);
  }

  // fill histo
  Int_t charge = esdTrack->Charge();

  //Double_t vRecMCTrackHist[11] = { nClust,chi2PerCluster,clustPerFindClust,b[0],b[1],eta,phi,pt,kinkIdx,isPrim, polarity };
  //fRecMCTrackHist->Fill(vRecMCTrackHist);

  Double_t vRecMCTrackHist[13] = { static_cast<Double_t>(nCrossedRowsTPC), chi2PerCluster, ratioCrossedRowsOverFindableClustersTPC, fracClustersTPCShared , b[0], b[1], eta, phi, pt, static_cast<Double_t>(hasStrangeMother), static_cast<Double_t>(isFromMaterial), static_cast<Double_t>(isPrim), static_cast<Double_t>(charge) };
  if(fFillSparseHisto) { fRecMCTrackHist->Fill(vRecMCTrackHist); }

  fCrossRowsEtaPt->Fill(static_cast<Double_t>(nCrossedRowsTPC),eta,pt);
  fChi2PerClustEtaPt->Fill(chi2PerCluster,eta,pt);
  fCrossRowsOverFindEtaPt->Fill(ratioCrossedRowsOverFindableClustersTPC,eta,pt);
  fFracSharedClustEtaPt->Fill(fracClustersTPCShared,eta,pt);
  fDCAyEtaPt->Fill(b[0],eta,pt);
  fDCAzEtaPt->Fill(b[1],eta,pt);

  fCrossRowsPhiPt->Fill(static_cast<Double_t>(nCrossedRowsTPC),phi,pt);
  fChi2PerClustPhiPt->Fill(chi2PerCluster,phi,pt);
  fCrossRowsOverFindPhiPt->Fill(ratioCrossedRowsOverFindableClustersTPC,phi,pt);
  fFracSharedClustPhiPt->Fill(fracClustersTPCShared,phi,pt);
  fDCAyPhiPt->Fill(b[0],phi,pt);
  fDCAzPhiPt->Fill(b[1],phi,pt);

  fCrossRowsEtaPhi->Fill(static_cast<Double_t>(nCrossedRowsTPC),eta,phi);
  fChi2PerClustEtaPhi->Fill(chi2PerCluster,eta,phi);
  fCrossRowsOverFindEtaPhi->Fill(ratioCrossedRowsOverFindableClustersTPC,eta,phi);
  fFracSharedClustEtaPhi->Fill(fracClustersTPCShared,eta,phi);
  fDCAyEtaPhi->Fill(b[0],eta,phi);
  fDCAzEtaPhi->Fill(b[1],eta,phi);

  if(IsUseMCInfo())
  {
    if(isPrim) {
      fDCAyEtaPtMCPrim->Fill(b[0],eta,pt);
      fDCAzEtaPtMCPrim->Fill(b[1],eta,pt);
      fDCAyPhiPtMCPrim->Fill(b[0],phi,pt);
      fDCAzPhiPtMCPrim->Fill(b[1],phi,pt);
      fDCAyEtaPhiMCPrim->Fill(b[0],eta,phi);
      fDCAzEtaPhiMCPrim->Fill(b[1],eta,phi);
    } else {
      fDCAyEtaPtMCSec->Fill(b[0],eta,pt);
      fDCAzEtaPtMCSec->Fill(b[1],eta,pt);
      fDCAyPhiPtMCSec->Fill(b[0],phi,pt);
      fDCAzPhiPtMCSec->Fill(b[1],phi,pt);
      fDCAyEtaPhiMCSec->Fill(b[0],eta,phi);
      fDCAzEtaPhiMCSec->Fill(b[1],eta,phi);

      if(hasStrangeMother) {
        fDCAyEtaPtSecMCDecays->Fill(b[0],eta,pt);
        fDCAzEtaPtSecMCDecays->Fill(b[1],eta,pt);
        fDCAyPhiPtSecMCDecays->Fill(b[0],phi,pt);
        fDCAzPhiPtSecMCDecays->Fill(b[1],phi,pt);
        fDCAyEtaPhiSecMCDecays->Fill(b[0],eta,phi);
        fDCAzEtaPhiSecMCDecays->Fill(b[1],eta,phi);

            if(isFromK0s) {
                fDCAyEtaPtSecMCDecaysK0s->Fill(b[0],eta,pt);
                fDCAzEtaPtSecMCDecaysK0s->Fill(b[1],eta,pt);
                fDCAyPhiPtSecMCDecaysK0s->Fill(b[0],phi,pt);
                fDCAzPhiPtSecMCDecaysK0s->Fill(b[1],phi,pt);
                fDCAyEtaPhiSecMCDecaysK0s->Fill(b[0],eta,phi);
                fDCAzEtaPhiSecMCDecaysK0s->Fill(b[1],eta,phi);
            }

            if(isFromLambda) {
                fDCAyEtaPtSecMCDecaysLambda->Fill(b[0],eta,pt);
                fDCAzEtaPtSecMCDecaysLambda->Fill(b[1],eta,pt);
                fDCAyPhiPtSecMCDecaysLambda->Fill(b[0],phi,pt);
                fDCAzPhiPtSecMCDecaysLambda->Fill(b[1],phi,pt);
                fDCAyEtaPhiSecMCDecaysLambda->Fill(b[0],eta,phi);
                fDCAzEtaPhiSecMCDecaysLambda->Fill(b[1],eta,phi);
            }
        } else {
            fDCAyEtaPtSecMCMaterial->Fill(b[0],eta,pt);
            fDCAzEtaPtSecMCMaterial->Fill(b[1],eta,pt);
            fDCAyPhiPtSecMCMaterial->Fill(b[0],phi,pt);
            fDCAzPhiPtSecMCMaterial->Fill(b[1],phi,pt);
            fDCAyEtaPhiSecMCMaterial->Fill(b[0],eta,phi);
            fDCAzEtaPhiSecMCMaterial->Fill(b[1],eta,phi);
        }
      }
    }
}

//_____________________________________________________________________________
Long64_t AlidNdPtCutAnalysis::Merge(TCollection* const list)
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  //TList *collPhysSelection = new TList;

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) {
    AlidNdPtCutAnalysis* entry = dynamic_cast<AlidNdPtCutAnalysis*>(obj);
    if (entry == 0) continue;

    // event histo
    fEventCount->Add(entry->fEventCount);
    fRecEventHist->Add(entry->fRecEventHist);
    fRecMCEventHist->Add(entry->fRecMCEventHist);
    fMCEventHist->Add(entry->fMCEventHist);

    fTrigVertCount->Add(fTrigVertCount);

    fRecEventXYZ->Add(fRecEventXYZ);
    fRecEventXYMult->Add(fRecEventXYMult);
    fRecEventXZMult->Add(fRecEventYZMult);
    fRecEventYZMult->Add(fRecEventZResZMult);
    fRecEventZResZMult->Add(fRecEventZResZMult);

    fMCEventXYZ->Add(fMCEventXYZ);
    fRecMCEventDXDYDZ->Add(fRecMCEventDXDYDZ);
    fRecMCEventDXDYMult->Add(fRecMCEventDXDYMult);
    fRecMCEventDXDZMult->Add(fRecMCEventDXDZMult);
    fRecMCEventDYDZMult->Add(fRecMCEventDYDZMult);

    // track histo
    fRecMCTrackHist->Add(entry->fRecMCTrackHist);

    //
    fCrossRowsEtaPt->Add(entry->fCrossRowsEtaPt);
    fChi2PerClustEtaPt->Add(entry->fChi2PerClustEtaPt);
    fCrossRowsOverFindEtaPt->Add(entry->fCrossRowsOverFindEtaPt);
    fFracSharedClustEtaPt->Add(entry->fFracSharedClustEtaPt);
    fDCAyEtaPt->Add(entry->fDCAyEtaPt);
    fDCAzEtaPt->Add(entry->fDCAzEtaPt);

    fCrossRowsPhiPt->Add(entry->fCrossRowsPhiPt);
    fChi2PerClustPhiPt->Add(entry->fChi2PerClustPhiPt);
    fCrossRowsOverFindPhiPt->Add(entry->fCrossRowsOverFindPhiPt);
    fFracSharedClustPhiPt->Add(entry->fFracSharedClustPhiPt);
    fDCAyPhiPt->Add(entry->fDCAyPhiPt);
    fDCAzPhiPt->Add(entry->fDCAzPhiPt);

    fCrossRowsEtaPhi->Add(entry->fCrossRowsEtaPhi);
    fChi2PerClustEtaPhi->Add(entry->fChi2PerClustEtaPhi);
    fCrossRowsOverFindEtaPhi->Add(entry->fCrossRowsOverFindEtaPhi);
    fFracSharedClustEtaPhi->Add(entry->fFracSharedClustEtaPhi);
    fDCAyEtaPhi->Add(entry->fDCAyEtaPhi);
    fDCAzEtaPhi->Add(entry->fDCAzEtaPhi);

    //
    fDCAyEtaPtMCPrim->Add(entry->fDCAyEtaPtMCPrim);
    fDCAzEtaPtMCPrim->Add(entry->fDCAzEtaPtMCPrim);
    fDCAyPhiPtMCPrim->Add(entry->fDCAyPhiPtMCPrim);
    fDCAzPhiPtMCPrim->Add(entry->fDCAzPhiPtMCPrim);
    fDCAyEtaPhiMCPrim->Add(entry->fDCAyEtaPhiMCPrim);
    fDCAzEtaPhiMCPrim->Add(entry->fDCAzEtaPhiMCPrim);

    fDCAyEtaPtMCSec->Add(entry->fDCAyEtaPtMCSec);
    fDCAzEtaPtMCSec->Add(entry->fDCAzEtaPtMCSec);
    fDCAyPhiPtMCSec->Add(entry->fDCAyPhiPtMCSec);
    fDCAzPhiPtMCSec->Add(entry->fDCAzPhiPtMCSec);
    fDCAyEtaPhiMCSec->Add(entry->fDCAyEtaPhiMCSec);
    fDCAzEtaPhiMCSec->Add(entry->fDCAzEtaPhiMCSec);

    fDCAyEtaPtSecMCDecays->Add(entry->fDCAyEtaPtSecMCDecays);
    fDCAzEtaPtSecMCDecays->Add(entry->fDCAzEtaPtSecMCDecays);
    fDCAyPhiPtSecMCDecays->Add(entry->fDCAyPhiPtSecMCDecays);
    fDCAzPhiPtSecMCDecays->Add(entry->fDCAzPhiPtSecMCDecays);
    fDCAyEtaPhiSecMCDecays->Add(entry->fDCAyEtaPhiSecMCDecays);
    fDCAzEtaPhiSecMCDecays->Add(entry->fDCAzEtaPhiSecMCDecays);

    fDCAyEtaPtSecMCDecaysK0s->Add(entry->fDCAyEtaPtSecMCDecaysK0s);
    fDCAzEtaPtSecMCDecaysK0s->Add(entry->fDCAzEtaPtSecMCDecaysK0s);
    fDCAyPhiPtSecMCDecaysK0s->Add(entry->fDCAyPhiPtSecMCDecaysK0s);
    fDCAzPhiPtSecMCDecaysK0s->Add(entry->fDCAzPhiPtSecMCDecaysK0s);
    fDCAyEtaPhiSecMCDecaysK0s->Add(entry->fDCAyEtaPhiSecMCDecaysK0s);
    fDCAzEtaPhiSecMCDecaysK0s->Add(entry->fDCAzEtaPhiSecMCDecaysK0s);

    fDCAyEtaPtSecMCDecaysLambda->Add(entry->fDCAyEtaPtSecMCDecaysLambda);
    fDCAzEtaPtSecMCDecaysLambda->Add(entry->fDCAzEtaPtSecMCDecaysLambda);
    fDCAyPhiPtSecMCDecaysLambda->Add(entry->fDCAyPhiPtSecMCDecaysLambda);
    fDCAzPhiPtSecMCDecaysLambda->Add(entry->fDCAzPhiPtSecMCDecaysLambda);
    fDCAyEtaPhiSecMCDecaysLambda->Add(entry->fDCAyEtaPhiSecMCDecaysLambda);
    fDCAzEtaPhiSecMCDecaysLambda->Add(entry->fDCAzEtaPhiSecMCDecaysLambda);

    fDCAyEtaPtSecMCMaterial->Add(entry->fDCAyEtaPtSecMCMaterial);
    fDCAzEtaPtSecMCMaterial->Add(entry->fDCAzEtaPtSecMCMaterial);
    fDCAyPhiPtSecMCMaterial->Add(entry->fDCAyPhiPtSecMCMaterial);
    fDCAzPhiPtSecMCMaterial->Add(entry->fDCAzPhiPtSecMCMaterial);
    fDCAyEtaPhiSecMCMaterial->Add(entry->fDCAyEtaPhiSecMCMaterial);
    fDCAzEtaPhiSecMCMaterial->Add(entry->fDCAzEtaPhiSecMCMaterial);

    // physics selection
    //collPhysSelection->Add(entry->GetPhysicsTriggerSelection());

  count++;
  }

  //AliPhysicsSelection *trigSelection = GetPhysicsTriggerSelection();
  //trigSelection->Merge(collPhysSelection);

  //if(collPhysSelection) delete collPhysSelection;

return count;
}

//_____________________________________________________________________________
void AlidNdPtCutAnalysis::Analyse()
{
  //
  // Analyse histograms
  //
  TH1::AddDirectory(kFALSE);
  TObjArray *aFolderObj = new TObjArray;
  if(!aFolderObj) return;

  TH1D *h1D = 0;
  TH2D *h2D = 0;


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

  Double_t maxDCAr = accCuts->GetMaxDCAr();

  //
  // Event counters
  //
  h2D = (TH2D*)fEventCount->Projection(0,1);
  if(!h2D) return;
  h2D->SetName("trig_vs_trigANDvertex");
  aFolderObj->Add(h2D);

  fEventCount->GetAxis(0)->SetRange(2,2); // triggered
  h1D = (TH1D*)fEventCount->Projection(1);
  if(!h1D) return;
  h1D->SetTitle("rec. vertex for triggered events");
  h1D->SetName("trigANDvertex");
  aFolderObj->Add(h1D);

  //
  // Create rec. event histograms
  //
  h1D = (TH1D *)fRecEventHist->Projection(0);
  if(!h1D) return;
  h1D->SetName("rec_xv");
  aFolderObj->Add(h1D);

  h1D = (TH1D *)fRecEventHist->Projection(1);
  if(!h1D) return;
  h1D->SetName("rec_yv");
  aFolderObj->Add(h1D);

  h1D = (TH1D *)fRecEventHist->Projection(2);
  if(!h1D) return;
  h1D->SetName("rec_zv");
  aFolderObj->Add(h1D);

  h2D = (TH2D *)fRecEventHist->Projection(3,4);
  if(!h2D) return;
  h2D->SetName("rec_resZv_vs_Mult");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecEventHist->Projection(0,1);
  if(!h2D) return;
  h2D->SetName("rec_xv_vs_yv");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecEventHist->Projection(0,2);
  if(!h2D) return;
  h2D->SetName("rec_xv_vs_zv");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecEventHist->Projection(3,4);
  if(!h2D) return;
  h2D->SetName("rec_resZv_vs_Mult");
  aFolderObj->Add(h2D);

  //
  // MC available
  //
  if(IsUseMCInfo()) {

  //
  // Create mc event histograms
  //
  h2D = (TH2D *)fMCEventHist->Projection(0,1);
  if(!h2D) return;
  h2D->SetName("mc_xv_vs_yv");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fMCEventHist->Projection(0,2);
  if(!h2D) return;
  h2D->SetName("mc_xv_vs_zv");
  aFolderObj->Add(h2D);

  //
  // Create rec-mc event histograms
  //
  h2D = (TH2D *)fRecMCEventHist->Projection(0,3);
  if(!h2D) return;
  h2D->SetName("rec_mc_deltaXv_vs_mult");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCEventHist->Projection(1,3);
  if(!h2D) return;
  h2D->SetName("rec_mc_deltaYv_vs_mult");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCEventHist->Projection(2,3);
  if(!h2D) return;
  h2D->SetName("rec_mc_deltaZv_vs_mult");
  aFolderObj->Add(h2D);

  } // end use MC info



  //
  // Create rec-mc track track histograms
  //

  // DCA cuts
  fRecMCTrackHist->GetAxis(3)->SetRangeUser(-maxDCAr,maxDCAr);
  fRecMCTrackHist->GetAxis(4)->SetRangeUser(-maxDCAr,maxDCAr);

  h2D = (TH2D *)fRecMCTrackHist->Projection(7,5);
  if(!h2D) return;
  h2D->SetName("pt_vs_eta");
  aFolderObj->Add(h2D);

  fRecMCTrackHist->GetAxis(7)->SetRangeUser(minPt,maxPt);

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,5);
  if(!h2D) return;
  h2D->SetName("nClust_vs_eta");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(1,5);
  if(!h2D) return;
  h2D->SetName("chi2PerClust_vs_eta");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(2,5);
  if(!h2D) return;
  h2D->SetName("ratio_nClust_nFindableClust_vs_eta");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(5,6);
  if(!h2D) return;
  h2D->SetName("eta_vs_phi");
  aFolderObj->Add(h2D);

  //
  fRecMCTrackHist->GetAxis(5)->SetRangeUser(minEta,maxEta);

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,6);
  if(!h2D) return;
  h2D->SetName("nClust_vs_phi");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(1,6);
  if(!h2D) return;
  h2D->SetName("chi2PerClust_vs_phi");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(2,6);
  if(!h2D) return;
  h2D->SetName("ratio_nClust_nFindableClust_vs_phi");
  aFolderObj->Add(h2D);

  //
  fRecMCTrackHist->GetAxis(7)->SetRangeUser(0.0,maxPt);

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,7);
  if(!h2D) return;
  h2D->SetName("nClust_vs_pt");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(1,7);
  if(!h2D) return;
  h2D->SetName("chi2PerClust_vs_pt");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(2,7);
  if(!h2D) return;
  h2D->SetName("ratio_nClust_nFindableClust_vs_pt");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(6,7);
  if(!h2D) return;
  h2D->SetName("phi_vs_pt");
  aFolderObj->Add(h2D);


  // fiducial volume
  fRecMCTrackHist->GetAxis(5)->SetRangeUser(minEta,maxEta);
  fRecMCTrackHist->GetAxis(7)->SetRangeUser(minPt,maxPt);

  // DCA cuts
  fRecMCTrackHist->GetAxis(3)->SetRangeUser(-maxDCAr,maxDCAr);
  fRecMCTrackHist->GetAxis(4)->SetRangeUser(-maxDCAr,maxDCAr);

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,1);
  if(!h2D) return;
  h2D->SetName("nClust_vs_chi2PerClust");
  aFolderObj->Add(h2D);

  h2D = (TH2D *)fRecMCTrackHist->Projection(0,2);
  if(!h2D) return;
  h2D->SetName("nClust_vs_ratio_nClust_nFindableClust");
  aFolderObj->Add(h2D);

  //
  // DCAy cuts
  //
  fRecMCTrackHist->GetAxis(0)->SetRange(50,160); // nClust/track > 50
  fRecMCTrackHist->GetAxis(1)->SetRangeUser(0.,3.9999); // chi2/cluster < 4.0
  fRecMCTrackHist->GetAxis(3)->SetRange(1,fRecMCTrackHist->GetAxis(3)->GetNbins());
  //fRecMCTrackHist->GetAxis(4)->SetRangeUser(-1.0,1.0);
  fRecMCTrackHist->GetAxis(4)->SetRange(1,fRecMCTrackHist->GetAxis(4)->GetNbins());

  // sec
  fRecMCTrackHist->GetAxis(9)->SetRange(1,1);
  h1D = (TH1D *)fRecMCTrackHist->Projection(3);
  if(!h1D) return;
  h1D->SetName("dcay_sec");
  aFolderObj->Add(h1D);

  // prim
  fRecMCTrackHist->GetAxis(9)->SetRange(2,2);
  h1D = (TH1D *)fRecMCTrackHist->Projection(3);
  if(!h1D) return;
  h1D->SetName("dcay_prim");
  aFolderObj->Add(h1D);

  // DCAz cuts
  //fRecMCTrackHist->GetAxis(3)->SetRangeUser(-1.0,1.0);
  fRecMCTrackHist->GetAxis(4)->SetRange(1,fRecMCTrackHist->GetAxis(4)->GetNbins());

  // sec
  fRecMCTrackHist->GetAxis(9)->SetRange(1,1);
  h1D = (TH1D *)fRecMCTrackHist->Projection(4);
  if(!h1D) return;
  h1D->SetName("dcaz_sec");
  aFolderObj->Add(h1D);

  // prim
  fRecMCTrackHist->GetAxis(9)->SetRange(2,2);
  h1D = (TH1D *)fRecMCTrackHist->Projection(4);
  if(!h1D) return;
  h1D->SetName("dcaz_prim");
  aFolderObj->Add(h1D);


  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);
  if(!fAnalysisFolder) {
      if(aFolderObj) delete aFolderObj;
      return;
  }

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AlidNdPtCutAnalysis::ExportToFolder(TObjArray * const array)
{
  // recreate folder avery time and export objects to new one
  //
  AlidNdPtCutAnalysis * comp=this;
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
TFolder* AlidNdPtCutAnalysis::CreateFolder(TString name,TString title) {
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
