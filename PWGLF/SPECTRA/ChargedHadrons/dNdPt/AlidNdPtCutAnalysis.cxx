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
// Modified: F.Sozzi 03/2016
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
#include "AliMultSelection.h"

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
  fDCAyCenPt(0),
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
  fDCAyCenPtMCPrim(0),
 //
  fDCAyEtaPtMCSec(0),
  fDCAzEtaPtMCSec(0),
  fDCAyPhiPtMCSec(0),
  fDCAzPhiPtMCSec(0),
  fDCAyEtaPhiMCSec(0),
  fDCAzEtaPhiMCSec(0),
  fDCAyCenPtMCSec(0),
  //
  fDCAyEtaPtMCSecDecays(0),
  fDCAzEtaPtMCSecDecays(0),
  fDCAyPhiPtMCSecDecays(0),
  fDCAzPhiPtMCSecDecays(0),
  fDCAyEtaPhiMCSecDecays(0),
  fDCAzEtaPhiMCSecDecays(0),
  fDCAyCenPtMCSecDecays(0),
 //
  fDCAyEtaPtMCSecDecaysK0s(0),
  fDCAzEtaPtMCSecDecaysK0s(0),
  fDCAyPhiPtMCSecDecaysK0s(0),
  fDCAzPhiPtMCSecDecaysK0s(0),
  fDCAyEtaPhiMCSecDecaysK0s(0),
  fDCAzEtaPhiMCSecDecaysK0s(0),
  fDCAyCenPtMCSecDecaysK0s(0),
  //
  fDCAyEtaPtMCSecDecaysLambda(0),
  fDCAzEtaPtMCSecDecaysLambda(0),
  fDCAyPhiPtMCSecDecaysLambda(0),
  fDCAzPhiPtMCSecDecaysLambda(0),
  fDCAyEtaPhiMCSecDecaysLambda(0),
  fDCAzEtaPhiMCSecDecaysLambda(0),
  fDCAyCenPtMCSecDecaysLambda(0),
 //
  fDCAyEtaPtMCSecMaterial(0),
  fDCAzEtaPtMCSecMaterial(0),
  fDCAyPhiPtMCSecMaterial(0),
  fDCAzPhiPtMCSecMaterial(0),
  fDCAyEtaPhiMCSecMaterial(0),
  fDCAzEtaPhiMCSecMaterial(0),
  fDCAyCenPtMCSecMaterial(0)
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
  fDCAyCenPt(0),
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
  fDCAyCenPtMCPrim(0),
  //
  fDCAyEtaPtMCSec(0),
  fDCAzEtaPtMCSec(0),
  fDCAyPhiPtMCSec(0),
  fDCAzPhiPtMCSec(0),
  fDCAyEtaPhiMCSec(0),
  fDCAzEtaPhiMCSec(0),
  fDCAyCenPtMCSec(0),
  //
  fDCAyEtaPtMCSecDecays(0),
  fDCAzEtaPtMCSecDecays(0),
  fDCAyPhiPtMCSecDecays(0),
  fDCAzPhiPtMCSecDecays(0),
  fDCAyEtaPhiMCSecDecays(0),
  fDCAzEtaPhiMCSecDecays(0),
  fDCAyCenPtMCSecDecays(0),
  //
  fDCAyEtaPtMCSecDecaysK0s(0),
  fDCAzEtaPtMCSecDecaysK0s(0),
  fDCAyPhiPtMCSecDecaysK0s(0),
  fDCAzPhiPtMCSecDecaysK0s(0),
  fDCAyEtaPhiMCSecDecaysK0s(0),
  fDCAzEtaPhiMCSecDecaysK0s(0),
  fDCAyCenPtMCSecDecaysK0s(0),
  //
  fDCAyEtaPtMCSecDecaysLambda(0),
  fDCAzEtaPtMCSecDecaysLambda(0),
  fDCAyPhiPtMCSecDecaysLambda(0),
  fDCAzPhiPtMCSecDecaysLambda(0),
  fDCAyEtaPhiMCSecDecaysLambda(0),
  fDCAzEtaPhiMCSecDecaysLambda(0),
  fDCAyCenPtMCSecDecaysLambda(0),
  //
  fDCAyEtaPtMCSecMaterial(0),
  fDCAzEtaPtMCSecMaterial(0),
  fDCAyPhiPtMCSecMaterial(0),
  fDCAzPhiPtMCSecMaterial(0),
  fDCAyEtaPhiMCSecMaterial(0),
  fDCAzEtaPhiMCSecMaterial(0),
  fDCAyCenPtMCSecMaterial(0)
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
    if(fDCAyCenPt) delete fDCAyCenPt;

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
    if(fDCAyCenPtMCPrim) delete fDCAyCenPtMCPrim;

    if(fDCAyEtaPtMCSec) delete fDCAyEtaPtMCSec;
    if(fDCAzEtaPtMCSec) delete fDCAzEtaPtMCSec;
    if(fDCAyPhiPtMCSec) delete fDCAyPhiPtMCSec;
    if(fDCAzPhiPtMCSec) delete fDCAzPhiPtMCSec;
    if(fDCAyEtaPhiMCSec) delete fDCAyEtaPhiMCSec;
    if(fDCAzEtaPhiMCSec) delete fDCAzEtaPhiMCSec;
    if(fDCAyCenPtMCSec) delete fDCAyCenPtMCSec;

    if(fDCAyEtaPtMCSecDecays) delete fDCAyEtaPtMCSecDecays;
    if(fDCAzEtaPtMCSecDecays) delete fDCAzEtaPtMCSecDecays;
    if(fDCAyPhiPtMCSecDecays) delete fDCAyPhiPtMCSecDecays;
    if(fDCAzPhiPtMCSecDecays) delete fDCAzPhiPtMCSecDecays;
    if(fDCAyEtaPhiMCSecDecays) delete fDCAyEtaPhiMCSecDecays;
    if(fDCAzEtaPhiMCSecDecays) delete fDCAzEtaPhiMCSecDecays;
    if(fDCAyCenPtMCSecDecays) delete fDCAyCenPtMCSecDecays;

    if(fDCAyEtaPtMCSecDecaysK0s) delete fDCAyEtaPtMCSecDecaysK0s;
    if(fDCAzEtaPtMCSecDecaysK0s) delete fDCAzEtaPtMCSecDecaysK0s;
    if(fDCAyPhiPtMCSecDecaysK0s) delete fDCAyPhiPtMCSecDecaysK0s;
    if(fDCAzPhiPtMCSecDecaysK0s) delete fDCAzPhiPtMCSecDecaysK0s;
    if(fDCAyEtaPhiMCSecDecaysK0s) delete fDCAyEtaPhiMCSecDecaysK0s;
    if(fDCAzEtaPhiMCSecDecaysK0s) delete fDCAzEtaPhiMCSecDecaysK0s;
    if(fDCAyCenPtMCSecDecaysK0s) delete fDCAyCenPtMCSecDecaysK0s;

    if(fDCAyEtaPtMCSecDecaysLambda) delete fDCAyEtaPtMCSecDecaysLambda;
    if(fDCAzEtaPtMCSecDecaysLambda) delete fDCAzEtaPtMCSecDecaysLambda;
    if(fDCAyPhiPtMCSecDecaysLambda) delete fDCAyPhiPtMCSecDecaysLambda;
    if(fDCAzPhiPtMCSecDecaysLambda) delete fDCAzPhiPtMCSecDecaysLambda;
    if(fDCAyEtaPhiMCSecDecaysLambda) delete fDCAyEtaPhiMCSecDecaysLambda;
    if(fDCAzEtaPhiMCSecDecaysLambda) delete fDCAzEtaPhiMCSecDecaysLambda;
    if(fDCAyCenPtMCSecDecaysLambda) delete fDCAyCenPtMCSecDecaysLambda;

    if(fDCAyEtaPtMCSecMaterial) delete fDCAyEtaPtMCSecMaterial;
    if(fDCAzEtaPtMCSecMaterial) delete fDCAzEtaPtMCSecMaterial;
    if(fDCAyPhiPtMCSecMaterial) delete fDCAyPhiPtMCSecMaterial;
    if(fDCAzPhiPtMCSecMaterial) delete fDCAzPhiPtMCSecMaterial;
    if(fDCAyEtaPhiMCSecMaterial) delete fDCAyEtaPhiMCSecMaterial;
    if(fDCAzEtaPhiMCSecMaterial) delete fDCAzEtaPhiMCSecMaterial;
    if(fDCAyCenPtMCSecMaterial) delete fDCAyCenPtMCSecMaterial;

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
  const Int_t ptNbins = 100;
  const Double_t ptMin = 0, ptMax = 10.;
  Double_t *binsPt = CreateLogAxis(ptNbins,ptMin+0.001,ptMax);

  // set DCAy bins
  // const Int_t DCAybins = 300;
  // const Double_t DCAyMin = -1, DCAyMax = 1;
   const Int_t DCAybins = 100;
   const Double_t DCAyMin = -3, DCAyMax = 3;

  // set centrality bins
  const Int_t cenNbins = 20;
  const Double_t cenMin = 0, cenMax = 100;
  
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

 //nCrossRows:chi2PerClust:nCrossRows/nFindableClust:fracSharedClust:DCAy:DCAz:eta:phi:pt:isWeakDecay:isFromMaterial:isPrim:charge
  Int_t binsRecMCTrackHist[13]=  { 160,  10,  20,  20, 50,  50,   20,  90,             ptNbins, 2,  2,  2,  3  };
  Double_t minRecMCTrackHist[13]={ 0.,   0.,  0.,  0., -0.5,-0.5,-1.0, 0.,             ptMin+0.001,   0., 0., 0.,-1. };
  Double_t maxRecMCTrackHist[13]={ 160., 10., 1.,  1., 0.5, 0.5,  1.0, 2.*TMath::Pi(), ptMax,   2., 2., 2., 2. };

  fRecMCTrackHist = new THnSparseF("fRecMCTrackHist","nCrossRows:chi2PerClust:nCrossRows/nFindableClust:fracSharedClust:DCAy:DCAz:eta:phi:pt:isWeakDecay:isFromMaterial:isPrim:charge",13,binsRecMCTrackHist,minRecMCTrackHist,maxRecMCTrackHist);
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
  fRecMCTrackHist->GetAxis(9)->SetTitle("isWeakDecay");
  fRecMCTrackHist->GetAxis(10)->SetTitle("isFromMaterial");
  fRecMCTrackHist->GetAxis(11)->SetTitle("isPrim");
  fRecMCTrackHist->GetAxis(12)->SetTitle("charge");
  fRecMCTrackHist->Sumw2();



  //
  fCrossRowsEtaPt = new TH3D("fCrossRowsEtaPt","nCrossRows:eta:pt; nCrossRows; eta; pt (GeV/c)", 160,0,160,20,-1,1,ptNbins,ptMin,ptMax);
  fChi2PerClustEtaPt = new TH3D("fChi2PerClustEtaPt","chi2PerClust:eta:pt; chi2PerClust; eta; pt (GeV/c)", 100,0,10,20,-1,1,ptNbins,ptMin,ptMax);
  fCrossRowsOverFindEtaPt = new TH3D("fCrossRowsOverFindEtaPt","nCrossRows/nFindableClust:eta:pt; nCrossRows/nFindableClust; eta; pt (GeV/c)", 100,0,2,20,-1,1,ptNbins,ptMin,ptMax);
  fFracSharedClustEtaPt = new TH3D("fFracSharedClustEtaPt","fracSharedClust:eta:pt;fracSharedClust; eta; pt (GeV/c)", 100,0,1,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyEtaPt = new TH3D("fDCAyEtaPt","DCAy:eta:pt;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPt = new TH3D("fDCAzEtaPt","DCAz:eta:pt;DCAz (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyCenPt = new TH3D("fDCAyCenPt","DCAy:cen:pt;DCAy (cm); cen; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,cenNbins,cenMin,cenMax,ptNbins,ptMin,ptMax);

  //
  fCrossRowsPhiPt = new TH3D("fCrossRowsPhiPt","nCrossRows:eta:pt; nCrossRows; phi (rad); pt (GeV/c)", 160,0,160,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fChi2PerClustPhiPt = new TH3D("fChi2PerClustPhiPt","chi2PerClust:eta:pt; chi2PerClust; phi (rad); pt (GeV/c)", 100,0,10,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fCrossRowsOverFindPhiPt = new TH3D("fCrossRowsOverFindPhiPt","nCrossRows/nFindableClust:eta:pt; nCrossRows/nFindableClust; phi (rad); pt (GeV/c)", 100,0,2,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fFracSharedClustPhiPt = new TH3D("fFracSharedClustPhiPt","fracSharedClust:eta:pt;fracSharedClust; phi (rad); pt (GeV/c)", 100,0,1,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyPhiPt = new TH3D("fDCAyPhiPt","DCAy:phi:pt;DCAy (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPt = new TH3D("fDCAzPhiPt","DCAz:phi:pt;DCAz (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);

  //
  fCrossRowsEtaPhi = new TH3D("fCrossRowsEtaPhi","nCrossRows:eta:phi; nCrossRows; eta; phi (rad)", 160,0,160,20,-1,1,90,0,2.*TMath::Pi());
  fChi2PerClustEtaPhi = new TH3D("fChi2PerClustEtaPhi","chi2PerClust:eta:phi; chi2PerClust; eta; phi (rad)", 100,0,10,20,-1,1,90,0,2.*TMath::Pi());
  fCrossRowsOverFindEtaPhi = new TH3D("fCrossRowsOverFindEtaPhi","nCrossRows/nFindableClust:eta:phi; nCrossRows/nFindableClust; eta; phi (rad)", 100,0,2,20,-1,1,90,0,2.*TMath::Pi());
  fFracSharedClustEtaPhi = new TH3D("fFracSharedClustEtaPhi","fracSharedClust:eta:phi;fracSharedClust; eta; phi (rad)", 100,0,1,20,-1,1,90,0,2.*TMath::Pi());
  fDCAyEtaPhi = new TH3D("fDCAyEtaPhi","DCAy:eta:phi;DCAy (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhi = new TH3D("fDCAzEtaPhi","DCAz:eta:phi;DCAz (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());

  //
  fDCAyEtaPtMCPrim = new TH3D("fDCAyEtaPtMCPrim","DCAy:eta:pt primary;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtMCPrim = new TH3D("fDCAzEtaPtMCPrim","DCAz:eta:pt primary;DCAz (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtMCPrim = new TH3D("fDCAyPhiPtMCPrim","DCAy primary:phi:pt primary;DCAy (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtMCPrim = new TH3D("fDCAzPhiPtMCPrim","DCAz:phi:pt primary;DCAz (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiMCPrim = new TH3D("fDCAyEtaPhiMCPrim","DCAy:eta:phi primary;DCAy (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiMCPrim = new TH3D("fDCAzEtaPhiMCPrim","DCAz:eta:phi primary;DCAz (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAyCenPtMCPrim = new TH3D("fDCAyCenPtMCPrim","DCAy:cen:pt primary;DCAy (cm); cen; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,cenNbins,cenMin,cenMax,ptNbins,ptMin,ptMax);

  //
  fDCAyEtaPtMCSec = new TH3D("fDCAyEtaPtMCSec","DCAy:eta:pt secondary;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtMCSec = new TH3D("fDCAzEtaPtMCSec","DCAz:eta:pt secondary;DCAz (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtMCSec = new TH3D("fDCAyPhiPtMCSec","DCAy secondary:phi:pt secondary;DCAy (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtMCSec = new TH3D("fDCAzPhiPtMCSec","DCAz:eta:pt secondary;DCAz (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiMCSec = new TH3D("fDCAyEtaPhiMCSec","DCAy:eta:phi secondary;DCAy (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiMCSec = new TH3D("fDCAzEtaPhiMCSec","DCAz:eta:phi secondary;DCAz (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAyCenPtMCSec = new TH3D("fDCAyCenPtMCSec","DCAy:cen:pt secondary;DCAy (cm); cen; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,cenNbins,cenMin,cenMax,ptNbins,ptMin,ptMax);

  fDCAyEtaPtMCSecDecays = new TH3D("fDCAyEtaPtMCSecDecays","DCAy:eta:pt secondary decays;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtMCSecDecays = new TH3D("fDCAzEtaPtMCSecDecays","DCAz:eta:pt secondary decays;DCAz (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtMCSecDecays = new TH3D("fDCAyPhiPtMCSecDecays","DCAy secondary decays:phi:pt secondary decays;DCAy (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtMCSecDecays = new TH3D("fDCAzPhiPtMCSecDecays","DCAz:eta:pt secondary decays;DCAz (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiMCSecDecays = new TH3D("fDCAyEtaPhiMCSecDecays","DCAy:eta:phi secondary decays;DCAy (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiMCSecDecays = new TH3D("fDCAzEtaPhiMCSecDecays","DCAz:eta:phi secondary decays;DCAz (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAyCenPtMCSecDecays = new TH3D("fDCAyCenPtMCSecDecays","DCAy:cen:pt secondary decays;DCAy (cm); cen; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,cenNbins,cenMin,cenMax,ptNbins,ptMin,ptMax);

  fDCAyEtaPtMCSecDecaysK0s = new TH3D("fDCAyEtaPtMCSecDecaysK0s","DCAy:eta:pt secondary decays from K0s;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtMCSecDecaysK0s = new TH3D("fDCAzEtaPtMCSecDecaysK0s","DCAz:eta:pt secondary decays from K0s;DCAz (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtMCSecDecaysK0s = new TH3D("fDCAyPhiPtMCSecDecaysK0s","DCAy secondary decays from K0s:phi:pt secondary decays from K0s;DCAy (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtMCSecDecaysK0s = new TH3D("fDCAzPhiPtMCSecDecaysK0s","DCAz:eta:pt secondary decays from K0s;DCAz (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiMCSecDecaysK0s = new TH3D("fDCAyEtaPhiMCSecDecaysK0s","DCAy:eta:phi secondary decays from K0s;DCAy (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiMCSecDecaysK0s = new TH3D("fDCAzEtaPhiMCSecDecaysK0s","DCAz:eta:phi secondary decays from K0s;DCAz (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAyCenPtMCSecDecaysK0s = new TH3D("fDCAyCenPtMCSecDecaysK0s","DCAy:cen:pt secondary decays from K0s;DCAy (cm); cen; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,cenNbins,cenMin,cenMax,ptNbins,ptMin,ptMax);

  fDCAyEtaPtMCSecDecaysLambda = new TH3D("fDCAyEtaPtMCSecDecaysLambda","DCAy:eta:pt secondary decays from Lambda;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtMCSecDecaysLambda = new TH3D("fDCAzEtaPtMCSecDecaysLambda","DCAz:eta:pt secondary decays from Lambda;DCAz (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtMCSecDecaysLambda = new TH3D("fDCAyPhiPtMCSecDecaysLambda","DCAy secondary decays from Lambda:phi:pt secondary decays from Lambda;DCAy (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtMCSecDecaysLambda = new TH3D("fDCAzPhiPtMCSecDecaysLambda","DCAz:eta:pt secondary decays from Lambda;DCAz (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiMCSecDecaysLambda = new TH3D("fDCAyEtaPhiMCSecDecaysLambda","DCAy:eta:phi secondary decays from Lambda;DCAy (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiMCSecDecaysLambda = new TH3D("fDCAzEtaPhiMCSecDecaysLambda","DCAz:eta:phi secondary decays from Lambda;DCAz (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAyCenPtMCSecDecaysLambda = new TH3D("fDCAyCenPtMCSecDecaysLambda","DCAy:cen:pt secondary decays from Lambda;DCAy (cm); cen; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,cenNbins,cenMin,cenMax,ptNbins,ptMin,ptMax);

  fDCAyEtaPtMCSecMaterial = new TH3D("fDCAyEtaPtMCSecMaterial","DCAy:eta:pt secondary material;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAzEtaPtMCSecMaterial = new TH3D("fDCAzEtaPtMCSecMaterial","DCAz:eta:pt secondary material;DCAz (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
  fDCAyPhiPtMCSecMaterial = new TH3D("fDCAyPhiPtMCSecMaterial","DCAy secondary material:phi:pt secondary material;DCAy (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAzPhiPtMCSecMaterial = new TH3D("fDCAzPhiPtMCSecMaterial","DCAz:eta:pt secondary material;DCAz (cm); phi (rad); pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,90,0,2.*TMath::Pi(),ptNbins,ptMin,ptMax);
  fDCAyEtaPhiMCSecMaterial = new TH3D("fDCAyEtaPhiMCSecMaterial","DCAy:phi:phi secondary material;DCAy (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAzEtaPhiMCSecMaterial = new TH3D("fDCAzEtaPhiMCSecMaterial","DCAz:phi:phi secondary material;DCAz (cm); eta; phi (rad)", DCAybins, DCAyMin, DCAyMax,20,-1,1,90,0,2.*TMath::Pi());
  fDCAyCenPtMCSecMaterial = new TH3D("fDCAyCenPtMCSecMaterial","DCAy:cen:pt secondary material;DCAy (cm); cen; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,cenNbins,cenMin,cenMax,ptNbins,ptMin,ptMax);

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

   // centrality determination 
   // In case centrality is not found, put the value at -1
   //NB This piece of code is running also in case of pp, better would be disable it 

   Float_t centralityF = -1.;
   AliMultSelection *MultSelection = (AliMultSelection*) esdEvent->FindListObject("MultSelection");
  
   if ( MultSelection ){
    centralityF = MultSelection->GetMultiplicityPercentile("V0M",kFALSE);
  }
  
  
  // use MC information
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);

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

      FillHistograms(track, stack,centralityF);
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
void AlidNdPtCutAnalysis::FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, const Float_t centralityF) const
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
  Bool_t isWeakDecay = kFALSE;
  Bool_t isFromMaterial = kFALSE;
  Bool_t isFromK0s = kFALSE;
  Bool_t isFromLambda = kFALSE;

  if(IsUseMCInfo()) {
    if(!stack) return;
    Int_t label = TMath::Abs(esdTrack->GetLabel());
    TParticle* particle = stack->Particle(label);
    if(!particle) return;
    if(particle->GetPDG() && particle->GetPDG()->Charge()==0.) return;

    //
    isPrim = stack->IsPhysicalPrimary(label);
    isWeakDecay = stack->IsSecondaryFromWeakDecay(label);
    isFromMaterial = stack->IsSecondaryFromMaterial(label);

    // check whether has stange mother
    //
    Int_t motherPdg = -1;
    TParticle* mother = 0;

    Int_t motherLabel = particle->GetFirstMother();
    if(motherLabel>=0) mother = stack->Particle(motherLabel);
    if(mother) motherPdg = TMath::Abs(mother->GetPdgCode()); // take abs for visualisation only
    Int_t mech = particle->GetUniqueID(); // production mechanism

    if(isWeakDecay && motherPdg == 310) isFromK0s = kTRUE;
    if(isWeakDecay && motherPdg == 3122) isFromLambda = kTRUE;
  }

  // fill histo
  Int_t charge = esdTrack->Charge();

  //Double_t vRecMCTrackHist[11] = { nClust,chi2PerCluster,clustPerFindClust,b[0],b[1],eta,phi,pt,kinkIdx,isPrim, polarity };
  //fRecMCTrackHist->Fill(vRecMCTrackHist);

  Double_t vRecMCTrackHist[13] = { static_cast<Double_t>(nCrossedRowsTPC), chi2PerCluster, ratioCrossedRowsOverFindableClustersTPC, fracClustersTPCShared , b[0], b[1], eta, phi, pt, static_cast<Double_t>(isWeakDecay), static_cast<Double_t>(isFromMaterial), static_cast<Double_t>(isPrim), static_cast<Double_t>(charge) };
  if(fFillSparseHisto) { fRecMCTrackHist->Fill(vRecMCTrackHist); }

  fCrossRowsEtaPt->Fill(static_cast<Double_t>(nCrossedRowsTPC),eta,pt);
  fChi2PerClustEtaPt->Fill(chi2PerCluster,eta,pt);
  fCrossRowsOverFindEtaPt->Fill(ratioCrossedRowsOverFindableClustersTPC,eta,pt);
  fFracSharedClustEtaPt->Fill(fracClustersTPCShared,eta,pt);
  fDCAyEtaPt->Fill(b[0],eta,pt);
  fDCAzEtaPt->Fill(b[1],eta,pt);
  fDCAyCenPt->Fill(b[0],centralityF,pt);

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
      fDCAyCenPtMCPrim->Fill(b[0],centralityF,pt);
    } else {
      fDCAyEtaPtMCSec->Fill(b[0],eta,pt);
      fDCAzEtaPtMCSec->Fill(b[1],eta,pt);
      fDCAyPhiPtMCSec->Fill(b[0],phi,pt);
      fDCAzPhiPtMCSec->Fill(b[1],phi,pt);
      fDCAyEtaPhiMCSec->Fill(b[0],eta,phi);
      fDCAzEtaPhiMCSec->Fill(b[1],eta,phi);
      fDCAyCenPtMCSec->Fill(b[0],centralityF,pt);

      if(isWeakDecay) {
        fDCAyEtaPtMCSecDecays->Fill(b[0],eta,pt);
        fDCAzEtaPtMCSecDecays->Fill(b[1],eta,pt);
        fDCAyPhiPtMCSecDecays->Fill(b[0],phi,pt);
        fDCAzPhiPtMCSecDecays->Fill(b[1],phi,pt);
        fDCAyEtaPhiMCSecDecays->Fill(b[0],eta,phi);
        fDCAzEtaPhiMCSecDecays->Fill(b[1],eta,phi);
        fDCAyCenPtMCSecDecays->Fill(b[0],centralityF,pt);

            if(isFromK0s) {
                fDCAyEtaPtMCSecDecaysK0s->Fill(b[0],eta,pt);
                fDCAzEtaPtMCSecDecaysK0s->Fill(b[1],eta,pt);
                fDCAyPhiPtMCSecDecaysK0s->Fill(b[0],phi,pt);
                fDCAzPhiPtMCSecDecaysK0s->Fill(b[1],phi,pt);
                fDCAyEtaPhiMCSecDecaysK0s->Fill(b[0],eta,phi);
                fDCAzEtaPhiMCSecDecaysK0s->Fill(b[1],eta,phi);
		fDCAyCenPtMCSecDecaysK0s->Fill(b[0],centralityF,pt);
    }

            if(isFromLambda) {
                fDCAyEtaPtMCSecDecaysLambda->Fill(b[0],eta,pt);
                fDCAzEtaPtMCSecDecaysLambda->Fill(b[1],eta,pt);
                fDCAyPhiPtMCSecDecaysLambda->Fill(b[0],phi,pt);
                fDCAzPhiPtMCSecDecaysLambda->Fill(b[1],phi,pt);
                fDCAyEtaPhiMCSecDecaysLambda->Fill(b[0],eta,phi);
                fDCAzEtaPhiMCSecDecaysLambda->Fill(b[1],eta,phi);
		fDCAyCenPtMCSecDecaysLambda->Fill(b[0],centralityF,pt);
   }
      }
      if(isFromMaterial) {
            fDCAyEtaPtMCSecMaterial->Fill(b[0],eta,pt);
            fDCAzEtaPtMCSecMaterial->Fill(b[1],eta,pt);
            fDCAyPhiPtMCSecMaterial->Fill(b[0],phi,pt);
            fDCAzPhiPtMCSecMaterial->Fill(b[1],phi,pt);
            fDCAyEtaPhiMCSecMaterial->Fill(b[0],eta,phi);
            fDCAzEtaPhiMCSecMaterial->Fill(b[1],eta,phi);
	    fDCAyCenPtMCSecMaterial->Fill(b[0],centralityF,pt);
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
    fRecEventXZMult->Add(fRecEventXZMult);
    fRecEventYZMult->Add(fRecEventYZMult);
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
    fDCAyCenPt->Add(entry->fDCAyCenPt);

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
    fDCAyCenPtMCPrim->Add(entry->fDCAyCenPtMCPrim);

    fDCAyEtaPtMCSec->Add(entry->fDCAyEtaPtMCSec);
    fDCAzEtaPtMCSec->Add(entry->fDCAzEtaPtMCSec);
    fDCAyPhiPtMCSec->Add(entry->fDCAyPhiPtMCSec);
    fDCAzPhiPtMCSec->Add(entry->fDCAzPhiPtMCSec);
    fDCAyEtaPhiMCSec->Add(entry->fDCAyEtaPhiMCSec);
    fDCAzEtaPhiMCSec->Add(entry->fDCAzEtaPhiMCSec);
    fDCAyCenPtMCSec->Add(entry->fDCAyCenPtMCSec);

    fDCAyEtaPtMCSecDecays->Add(entry->fDCAyEtaPtMCSecDecays);
    fDCAzEtaPtMCSecDecays->Add(entry->fDCAzEtaPtMCSecDecays);
    fDCAyPhiPtMCSecDecays->Add(entry->fDCAyPhiPtMCSecDecays);
    fDCAzPhiPtMCSecDecays->Add(entry->fDCAzPhiPtMCSecDecays);
    fDCAyEtaPhiMCSecDecays->Add(entry->fDCAyEtaPhiMCSecDecays);
    fDCAzEtaPhiMCSecDecays->Add(entry->fDCAzEtaPhiMCSecDecays);
    fDCAyCenPtMCSecDecays->Add(entry->fDCAyCenPtMCSecDecays);

    fDCAyEtaPtMCSecDecaysK0s->Add(entry->fDCAyEtaPtMCSecDecaysK0s);
    fDCAzEtaPtMCSecDecaysK0s->Add(entry->fDCAzEtaPtMCSecDecaysK0s);
    fDCAyPhiPtMCSecDecaysK0s->Add(entry->fDCAyPhiPtMCSecDecaysK0s);
    fDCAzPhiPtMCSecDecaysK0s->Add(entry->fDCAzPhiPtMCSecDecaysK0s);
    fDCAyEtaPhiMCSecDecaysK0s->Add(entry->fDCAyEtaPhiMCSecDecaysK0s);
    fDCAzEtaPhiMCSecDecaysK0s->Add(entry->fDCAzEtaPhiMCSecDecaysK0s);
    fDCAyCenPtMCSecDecaysK0s->Add(entry->fDCAyCenPtMCSecDecaysK0s);

    fDCAyEtaPtMCSecDecaysLambda->Add(entry->fDCAyEtaPtMCSecDecaysLambda);
    fDCAzEtaPtMCSecDecaysLambda->Add(entry->fDCAzEtaPtMCSecDecaysLambda);
    fDCAyPhiPtMCSecDecaysLambda->Add(entry->fDCAyPhiPtMCSecDecaysLambda);
    fDCAzPhiPtMCSecDecaysLambda->Add(entry->fDCAzPhiPtMCSecDecaysLambda);
    fDCAyEtaPhiMCSecDecaysLambda->Add(entry->fDCAyEtaPhiMCSecDecaysLambda);
    fDCAzEtaPhiMCSecDecaysLambda->Add(entry->fDCAzEtaPhiMCSecDecaysLambda);
    fDCAyCenPtMCSecDecaysLambda->Add(entry->fDCAyCenPtMCSecDecaysLambda);

    fDCAyEtaPtMCSecMaterial->Add(entry->fDCAyEtaPtMCSecMaterial);
    fDCAzEtaPtMCSecMaterial->Add(entry->fDCAzEtaPtMCSecMaterial);
    fDCAyPhiPtMCSecMaterial->Add(entry->fDCAyPhiPtMCSecMaterial);
    fDCAzPhiPtMCSecMaterial->Add(entry->fDCAzPhiPtMCSecMaterial);
    fDCAyEtaPhiMCSecMaterial->Add(entry->fDCAyEtaPhiMCSecMaterial);
    fDCAzEtaPhiMCSecMaterial->Add(entry->fDCAzEtaPhiMCSecMaterial);
    fDCAyCenPtMCSecMaterial->Add(entry->fDCAyCenPtMCSecMaterial);

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
