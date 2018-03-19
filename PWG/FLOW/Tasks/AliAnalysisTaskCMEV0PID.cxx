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

/* $Id: AliAnalysisTaskCMEV0PID.cxx  Rihan Haque 14/02/2018 $ */

//-- general include---
#include "TChain.h"
#include "TTree.h"
#include "TGrid.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TMatrixDSym.h"

#include "TMath.h"
#include "stdio.h"
#include "Riostream.h"

//---- manager and handler---
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

//---V0 and ZDC info---
#include "AliAODZDC.h"
#include "AliAODVZERO.h"
#include "AliAODVertex.h"

//---AOD,ESD event--
#include "AliESDEvent.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"

//----for PID-----
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"

//----- Vevent and tracks
#include "AliVEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"

//----- must include-------
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliPhysicsSelection.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskCMEV0PID.h"


using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskCMEV0PID)


AliAnalysisTaskCMEV0PID::AliAnalysisTaskCMEV0PID(const char *name): AliAnalysisTaskSE(name),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  mfileFBHijing(NULL),
  fListFBHijing(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),
  fHistTaskConfigParameters(NULL),
  fHistPileUpCount(NULL),
  fHistMultSelPUCount(NULL),
  fHistEtaPtBefore(NULL),
  fHistEtaPtAfter(NULL),
  fHistTPCvsGlobalMultBefore(NULL),
  fHistTPCvsGlobalMultAfter(NULL),
  fHistTPCdEdxvsPBefore(NULL),
  fHistTPCdEdxvsPAfter(NULL),
  fHistTOFBetavsPBefore(NULL),
  fHistTOFBetavsPAfter(NULL),
  fHistTOFMassvsPtBefore(NULL),
  fHistTOFMatchCount(NULL),
  fHistTPCVsESDTrkBefore(NULL),
  fHistTPCVsESDTrkAfter(NULL),
  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsV0MBefore(NULL),
  fHistTPConlyVsCL1After(NULL),
  fHistTPConlyVsV0MAfter(NULL),
  fHistGlobalVsV0MBefore(NULL),
  fHistGlobalVsV0MAfter(NULL),
  fHistRawVsCorrMultFB(NULL),
  hCentvsTPCmultCuts(NULL),
  fHV0AEventPlaneVsCent(NULL),
  fHV0CEventPlaneVsCent(NULL),
  fHTPCAEventPlaneVsCent(NULL),
  fHTPCCEventPlaneVsCent(NULL),
  fV0MultChVsRun(NULL),
  fCentDistBefore(NULL),
  fCentDistAfter(NULL),
  fHCorrectV0M(NULL),
  fHAvgerageQnV0A(NULL),  
  fHAvgerageQnV0C(NULL),
  fV0AQ2xVsCentRun(NULL),
  fV0AQ2yVsCentRun(NULL),
  fV0CQ2xVsCentRun(NULL),
  fV0CQ2yVsCentRun(NULL),
  fV0AQ3xVsCentRun(NULL),
  fV0AQ3yVsCentRun(NULL),
  fV0CQ3xVsCentRun(NULL),
  fV0CQ3yVsCentRun(NULL),
  fTPCAQ2xVsCentRun(NULL),
  fTPCAQ2yVsCentRun(NULL),
  fTPCCQ2xVsCentRun(NULL),
  fTPCCQ2yVsCentRun(NULL),
  fTPCAQ3xVsCentRun(NULL),
  fTPCAQ3yVsCentRun(NULL),
  fTPCCQ3xVsCentRun(NULL),
  fTPCCQ3yVsCentRun(NULL),
  fTPCAQ4xVsCentRun(NULL),
  fTPCAQ4yVsCentRun(NULL),
  fTPCCQ4xVsCentRun(NULL),
  fTPCCQ4yVsCentRun(NULL),
  fFilterBit(1),
  gN(1),
  gM(1),
  gPsiN(2),
  fOldRunNum(111),
  fEventCount(0),
  fNSigmaCut(2),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fCentralityPercentMin(0),
  fCentralityPercentMax(90),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  bApplyMCcorr(kFALSE),
  bV0MGainCorr(kFALSE),
  bSkipPileUpCut(kFALSE), 
  bFillNUAHistPID(kFALSE),
  sPathOfMCFile("/alien"),
  sNucleiTP("PbPb"),
  fHistEventCount(NULL)
{
  for(int i=0;i<3;i++){
    fHistPtwithTPCNsigma[i]=NULL;
    fHistPtwithTOFmasscut[i]=NULL;
    fHistPtwithTOFSignal[i]=NULL;
    fHistTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCnSigmavsPtAfter[i]=NULL;
    fHistTPCTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCdEdxvsPtPIDAfter[i]=NULL;
  }
  for(int i=0;i<5;i++){
    fHCorrectNUApos[i] = NULL;
    fHCorrectNUAneg[i] = NULL;
  }
  for(int i=0;i<5;i++){  // for PID 
    fHCorrectNUAposPion[i] = NULL;
    fHCorrectNUAnegPion[i] = NULL;
    fHCorrectNUAposKaon[i] = NULL;
    fHCorrectNUAnegKaon[i] = NULL;
    fHCorrectNUAposProton[i] = NULL;
    fHCorrectNUAnegProton[i] = NULL;
  }

  for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr3p_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_NN[i][j]  =  NULL;
    }
    for(int j=0;j<4;j++) {
      fHist_Reso2n_EP_Norm_Det[i][j] =  NULL;
    }
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr3p_Pion_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_NN[i][j]  =  NULL;
    }
  }
  for(int i=0;i<4;i++){
    for(int j=0;j<5;j++){
      fHist3DEtaPhiVz_Pos_Run[i][j]=NULL;
      fHist3DEtaPhiVz_Neg_Run[i][j]=NULL;
    }
  }
  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
  }


  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}
//______________________________empty constructor_______________________
AliAnalysisTaskCMEV0PID::AliAnalysisTaskCMEV0PID():
  AliAnalysisTaskSE(),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  mfileFBHijing(NULL),
  fListFBHijing(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),
  fHistTaskConfigParameters(NULL),
  fHistPileUpCount(NULL),
  fHistMultSelPUCount(NULL),
  fHistEtaPtBefore(NULL),
  fHistEtaPtAfter(NULL),
  fHistTPCvsGlobalMultBefore(NULL),
  fHistTPCvsGlobalMultAfter(NULL),
  fHistTPCdEdxvsPBefore(NULL),
  fHistTPCdEdxvsPAfter(NULL),
  fHistTOFBetavsPBefore(NULL),
  fHistTOFBetavsPAfter(NULL),
  fHistTOFMassvsPtBefore(NULL),
  fHistTOFMatchCount(NULL),
  fHistTPCVsESDTrkBefore(NULL),
  fHistTPCVsESDTrkAfter(NULL),
  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsV0MBefore(NULL),
  fHistTPConlyVsCL1After(NULL),
  fHistTPConlyVsV0MAfter(NULL),
  fHistGlobalVsV0MBefore(NULL),
  fHistGlobalVsV0MAfter(NULL),
  fHistRawVsCorrMultFB(NULL),
  hCentvsTPCmultCuts(NULL),
  fHV0AEventPlaneVsCent(NULL),
  fHV0CEventPlaneVsCent(NULL),
  fHTPCAEventPlaneVsCent(NULL),
  fHTPCCEventPlaneVsCent(NULL),
  fV0MultChVsRun(NULL),
  fCentDistBefore(NULL),
  fCentDistAfter(NULL),
  fHCorrectV0M(NULL),
  fHAvgerageQnV0A(NULL),  
  fHAvgerageQnV0C(NULL),
  fV0AQ2xVsCentRun(NULL),
  fV0AQ2yVsCentRun(NULL),
  fV0CQ2xVsCentRun(NULL),
  fV0CQ2yVsCentRun(NULL),
  fV0AQ3xVsCentRun(NULL),
  fV0AQ3yVsCentRun(NULL),
  fV0CQ3xVsCentRun(NULL),
  fV0CQ3yVsCentRun(NULL),  
  fTPCAQ2xVsCentRun(NULL),
  fTPCAQ2yVsCentRun(NULL),
  fTPCCQ2xVsCentRun(NULL),
  fTPCCQ2yVsCentRun(NULL),
  fTPCAQ3xVsCentRun(NULL),
  fTPCAQ3yVsCentRun(NULL),
  fTPCCQ3xVsCentRun(NULL),
  fTPCCQ3yVsCentRun(NULL),
  fTPCAQ4xVsCentRun(NULL),
  fTPCAQ4yVsCentRun(NULL),
  fTPCCQ4xVsCentRun(NULL),
  fTPCCQ4yVsCentRun(NULL),
  fFilterBit(1),
  gN(1),
  gM(1),
  gPsiN(2),
  fOldRunNum(111),
  fEventCount(0),
  fNSigmaCut(2),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fCentralityPercentMin(0),
  fCentralityPercentMax(90),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  bApplyMCcorr(kFALSE),
  bV0MGainCorr(kFALSE),
  bSkipPileUpCut(kFALSE), 
  bFillNUAHistPID(kFALSE),
  sPathOfMCFile("/alien"),
  sNucleiTP("PbPb"),
  fHistEventCount(NULL)
{
  for(int i=0;i<3;i++){
    fHistPtwithTPCNsigma[i]=NULL;
    fHistPtwithTOFmasscut[i]=NULL;
    fHistPtwithTOFSignal[i]=NULL;
    fHistTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCnSigmavsPtAfter[i]=NULL;
    fHistTPCTOFnSigmavsPtAfter[i]=NULL;
    fHistTPCdEdxvsPtPIDAfter[i]=NULL;
  }
  for(int i=0;i<5;i++){
    fHCorrectNUApos[i] = NULL;
    fHCorrectNUAneg[i] = NULL;
  }
  for(int i=0;i<5;i++){  // for PID 
    fHCorrectNUAposPion[i] = NULL;
    fHCorrectNUAnegPion[i] = NULL;
    fHCorrectNUAposKaon[i] = NULL;
    fHCorrectNUAnegKaon[i] = NULL;
    fHCorrectNUAposProton[i] = NULL;
    fHCorrectNUAnegProton[i] = NULL;
  }

  for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
      fHist_Corr3p_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_EP_Norm_NN[i][j]  =  NULL;
    }
    for(int j=0;j<4;j++) {
      fHist_Reso2n_EP_Norm_Det[i][j] =  NULL;
    }
    //PID:
    for(int j=0;j<4;j++){
      fHist_Corr3p_Pion_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Pion_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Kaon_EP_Norm_NN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_PN[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_PP[i][j]  =  NULL;
      fHist_Corr3p_Proton_EP_Norm_NN[i][j]  =  NULL;
    }
  }
  for(int i=0;i<4;i++){
    for(int j=0;j<5;j++){
      fHist3DEtaPhiVz_Pos_Run[i][j]=NULL;
      fHist3DEtaPhiVz_Neg_Run[i][j]=NULL;
    }
  }
  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
  }
}

//___________________________ destructor ___________________________
AliAnalysisTaskCMEV0PID::~AliAnalysisTaskCMEV0PID()
{
  //Destructor
  //if(fPIDResponse)   delete fPIDResponse;
  //if(fMultSelection) delete fMultSelection;

  if(fListHist){
    delete fListHist;
  }  

  if(fHCorrectV0M)    delete fHCorrectV0M;
  if(fHAvgerageQnV0A) delete fHAvgerageQnV0A;
  if(fHAvgerageQnV0C) delete fHAvgerageQnV0C;

  if(mfileFBHijing->IsOpen()){
     mfileFBHijing->Close();
     if(fListFBHijing) delete fListFBHijing;
  }
  for(int i=0;i<10;i++){
    if(fFB_Efficiency_Cent[i])
      delete fFB_Efficiency_Cent[i];
  }
  for(int i=0;i<5;i++){
    if(fHCorrectNUApos[i]) delete fHCorrectNUApos[i];
    if(fHCorrectNUAneg[i]) delete fHCorrectNUAneg[i];
  }
  for(int i=0;i<5;i++){  // for PID 
    if(fHCorrectNUAposPion[i])   delete fHCorrectNUAposPion[i];
    if(fHCorrectNUAnegPion[i])   delete fHCorrectNUAnegPion[i];
    if(fHCorrectNUAposKaon[i])   delete fHCorrectNUAposKaon[i];
    if(fHCorrectNUAnegKaon[i])   delete fHCorrectNUAnegKaon[i];
    if(fHCorrectNUAposProton[i]) delete fHCorrectNUAposProton[i];
    if(fHCorrectNUAnegProton[i]) delete fHCorrectNUAnegProton[i];
  }


  if(fAnalysisUtil)  delete fAnalysisUtil; // its 'new' !!
  
}//---------------- sanity ------------------------



void AliAnalysisTaskCMEV0PID::UserCreateOutputObjects()
{

  //std::cout<<"\n..UserCreateOutputObject called.. with isCorr = "<<isCorr<<"\n ....check if succeeded...\n"<<endl; 

  //input hander
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  if (!inputHandler) {  printf("\n ...Input handler missing!!!...\n");    return; }

  //PileUp Multi-Vertex
  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);

  //pid response object
  fPIDResponse=inputHandler->GetPIDResponse();
    
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);


  SetupEventAndTaskConfigInfo();

  if(!gGrid){
    TGrid::Connect("alien://");
  }

  SetupMCcorrectionMap(sPathOfMCFile);
 


  hCentvsTPCmultCuts = new TH2F("hCentvsTPCmultCuts","TPCmult Low,high",100,0,100,5,0,5);
  hCentvsTPCmultCuts->GetYaxis()->SetBinLabel(1,"mean"); 
  hCentvsTPCmultCuts->GetYaxis()->SetBinLabel(2,"sigma");   
  fListHist->Add(hCentvsTPCmultCuts);

  SetUpCentralityOutlierCut();

  fHistTOFMatchCount = new TH2F("fHistTOFMatchCount","TofMatchFlag vs mismatch Prob",10,0,10,200,-5,5);
  fListHist->Add(fHistTOFMatchCount);

  fHistEtaPtBefore = new TH2F("fHistEtaPtBefore","Eta vs pT",100,-1.25,1.25,100,0,10);
  fListHist->Add(fHistEtaPtBefore);

  fHistEtaPtAfter  = new TH2F("fHistEtaPtAfter","Eta vs pT",100,-1.25,1.25,100,0,10);
  fListHist->Add(fHistEtaPtAfter);


  Int_t gMaxTPCFB1mult  = 0;
  Int_t gMaxGlobalmult  = 0;
  Int_t gMaxTPCcorrmult = 0;
  Int_t gMaxESDtracks   = 0;

  if(sNucleiTP=="pp"||sNucleiTP=="PP"){  
    gMaxGlobalmult  = 200;
    gMaxTPCFB1mult  = 200;
    gMaxTPCcorrmult = 500;
    gMaxESDtracks   = 1000;
    //fSkipOutlierCut = 1;
  }
  else if(sNucleiTP=="pPb"||sNucleiTP=="Pbp"||sNucleiTP=="PbP"||sNucleiTP=="PPb"){  
    gMaxGlobalmult  = 400;
    gMaxTPCFB1mult  = 400;
    gMaxTPCcorrmult = 500;
    gMaxESDtracks   = 2000;
    //fSkipOutlierCut = 1;
  }
  else{
    gMaxGlobalmult  = 4000;
    gMaxTPCFB1mult  = 4000;
    gMaxTPCcorrmult = 5000;
    gMaxESDtracks   = 20000;
    //fSkipOutlierCut =  0;
  }

  //if(bSkipPileUpCut) { fSkipOutlierCut = 1;}


  fHistTPCVsESDTrkBefore = new TH2F("fHistTPCVsESDTrkBefore","Before; TPC1; ESD trk",500,0,gMaxTPCcorrmult,200,0,gMaxESDtracks);
  fListHist->Add(fHistTPCVsESDTrkBefore);
  fHistTPCVsESDTrkAfter  = new TH2F("fHistTPCVsESDTrkAfter"," After;  TPC1; ESD trk",500,0,gMaxTPCcorrmult,200,0,gMaxESDtracks);
  fListHist->Add(fHistTPCVsESDTrkAfter);

  fHistTPCvsGlobalMultBefore = new TH2F("fHistTPCvsGlobalMultBefore","Before; Global; TPC(fb1) ",200,0,gMaxGlobalmult,200,0,gMaxTPCFB1mult);
  fListHist->Add(fHistTPCvsGlobalMultBefore);
  fHistTPCvsGlobalMultAfter  = new TH2F("fHistTPCvsGlobalMultAfter"," After;  Global; TPC(fb1) ",200,0,gMaxGlobalmult,200,0,gMaxTPCFB1mult);
  fListHist->Add(fHistTPCvsGlobalMultAfter);

  fHistGlobalVsV0MBefore = new TH2F("fHistGlobalVsV0MBefore","Before;Cent(V0M);Global",100,0,100,500,0,gMaxGlobalmult);
  fListHist->Add(fHistGlobalVsV0MBefore);
  fHistGlobalVsV0MAfter  = new TH2F("fHistGlobalVsV0MAfter"," After; Cent(V0M);Global",100,0,100,500,0,gMaxGlobalmult);
  fListHist->Add(fHistGlobalVsV0MAfter);

  fHistTPConlyVsCL1Before = new TH2F("fHistTPConlyVsCL1Before","Before;Cent(CL1); TPC(FB1)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1Before);
  fHistTPConlyVsCL1After  = new TH2F("fHistTPConlyVsCL1After","After; Cent(CL1); TPC(FB1) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1After);

  fHistTPConlyVsV0MBefore = new TH2F("fHistTPConlyVsV0MBefore","Before;Cent(V0M); TPC(FB1)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MBefore);
  fHistTPConlyVsV0MAfter  = new TH2F("fHistTPConlyVsV0MAfter","After; Cent(V0M); TPC(FB1) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MAfter);

  fHistRawVsCorrMultFB = new TH2F("fHistRawVsCorrMultFB",Form("FB%d;Mult_{raw};Mult_{corr}",fFilterBit),gMaxTPCFB1mult,0,gMaxTPCFB1mult,gMaxTPCcorrmult,0,gMaxTPCcorrmult);
  fListHist->Add(fHistRawVsCorrMultFB);


  fHistTPCdEdxvsPBefore = new TH2F("fHistTPCdEdxvsPBefore","Before; p (GeV/c); dEdx (arb)",200,-5,5,200,0,250);
  fListHist->Add(fHistTPCdEdxvsPBefore);
  fHistTPCdEdxvsPAfter  = new TH2F("fHistTPCdEdxvsPAfter"," After;  p (GeV/c); dEdx (arb)",200,-5,5, 200,0,250);
  fListHist->Add(fHistTPCdEdxvsPAfter);


  fHistTOFBetavsPBefore = new TH2F("fHistTOFBetavsPBefore","Before; p (GeV/c); beta ",200,-5,5,100,0.0,1.2);
  fListHist->Add(fHistTOFBetavsPBefore);
  fHistTOFBetavsPAfter  = new TH2F("fHistTOFBetavsPAfter"," After;  p (GeV/c); beta ",200,-5,5,100,0.0,1.2);
  fListHist->Add(fHistTOFBetavsPAfter);


  fHistTOFMassvsPtBefore = new TH2F("fHistTOFMassvsPtBefore","Before; p_{T}(GeV/c); m^{2}(GeV^{2}/c^{4})",200,-5,5,500,-0.5,4.5);
  fListHist->Add(fHistTOFMassvsPtBefore);

  fCentDistBefore = new TH1F("fCentDistBefore","no Cut; Cent (%); Events ",100,0,100);
  fListHist->Add(fCentDistBefore);

  fCentDistAfter = new TH1F("fCentDistAfter","with Cut; Cent (%); Events ",100,0,100);
  fListHist->Add(fCentDistAfter);

  
  //---------------- PID Histograms ---------------------
  //Dont forget to add histograms in the List. !!!


  char const *gSpecies[4] = {"Pion","Kaon","proton","Charge"};

  for(int i=0;i<3;i++){
    fHistPtwithTPCNsigma[i]  = new TH1F(Form("fHistPtwithTPCNsigma_%s",gSpecies[i]), Form("%s;p_{T}(GeV/c))",gSpecies[i]),200,-5,5);
    fListHist->Add(fHistPtwithTPCNsigma[i]);
    fHistPtwithTOFmasscut[i] = new TH1F(Form("fHistPtwithTOFmasscut_%s",gSpecies[i]),Form("%s;p_{T}(GeV/c))",gSpecies[i]),200,-5,5);
    fListHist->Add(fHistPtwithTOFmasscut[i]);
    fHistPtwithTOFSignal[i]  = new TH1F(Form("fHistPtwithTOFSignal_%s", gSpecies[i]),Form("%s;p_{T}(GeV/c))",gSpecies[i]),200,-5,5);
    fListHist->Add(fHistPtwithTOFSignal[i]);

    fHistTOFnSigmavsPtAfter[i] = new TH2F(Form("fHistTOFnSigmavsPtAfter_%s",gSpecies[i]),Form("%s;p_{T}(GeV/c);n#sigma_{TOF}",gSpecies[i]),200,-5,5,400,-10.0,10.0);
    fListHist->Add(fHistTOFnSigmavsPtAfter[i]);
    fHistTPCnSigmavsPtAfter[i] = new TH2F(Form("fHistTPCnSigmavsPtAfter_%s",gSpecies[i]),Form("%s;p_{T}(GeV/c);n#sigma_{TPC}",gSpecies[i]),200,-5,5,400,-10.0,10.0);
    fListHist->Add(fHistTPCnSigmavsPtAfter[i]);

    fHistTPCTOFnSigmavsPtAfter[i] = new TH3F(Form("fHistTPCTOFnSigmavsPtAfter_%s",gSpecies[i]),Form("%s; p_{T}(GeV/c); n#sigma_{TPC}; n#sigma_{TOF}",gSpecies[i]),100,0,5,400,-10,10,400,-10,10);
    fListHist->Add(fHistTPCTOFnSigmavsPtAfter[i]);

    fHistTPCdEdxvsPtPIDAfter[i] = new TH2F(Form("fHistTPCdEdxvsPtAfter_%s",gSpecies[i]),"AfterCut; p_{T} (GeV/c); dEdx (arb)",400,0,10,200,0,250);
    fListHist->Add(fHistTPCdEdxvsPtPIDAfter[i]);
  }// PID histograms done



  Double_t centRange[11]   = {0,5,10,20,30,40,50,60,70,80,90};
  char const *gDetForEP[4] = {"V0A","V0C","TPC-A","TPC-C"};

 //------------------- CME 3p correlator Charged hadrons (EP method) ------------------
  for(int i=0;i<2;i++){
    for(int j=0;j<4;j++){
     //Detector: 0 = V0A, 1 = V0C, 3 = TPCA, 4 = TPCC 
      fHist_Corr3p_EP_Norm_PN[i][j] = new TProfile(Form("fHist_Corr3p_EP_Norm_PosNeg_Mag%d_Det%d",i,j+1),Form("US, #Psi_{2} %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_EP_Norm_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EP_Norm_PN[i][j]);
      fHist_Corr3p_EP_Norm_PP[i][j] = new TProfile(Form("fHist_Corr3p_EP_Norm_PosPos_Mag%d_Det%d",i,j+1),Form("P-P, #Psi_{2} %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_EP_Norm_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EP_Norm_PP[i][j]);
      fHist_Corr3p_EP_Norm_NN[i][j] = new TProfile(Form("fHist_Corr3p_EP_Norm_NegNeg_Mag%d_Det%d",i,j+1),Form("N-N, #Psi_{2}, %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_EP_Norm_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_EP_Norm_NN[i][j]);
    }
    //EP Resolution:
    for(int j=0;j<4;j++){
    //Det: 0 = v0c-v0a, 1 = v0a-TPC, 2 = v0c-TPC, 3 =TPC-A TPC-C
      fHist_Reso2n_EP_Norm_Det[i][j]  = new TProfile(Form("fHist_Reso2n_EP_Norm_Mag%d_DetComb%d",i,j+1),"Event plane Resolution",10,centRange,"");
      //fHist_Reso2n_EP_Norm_Det[i][j]->Sumw2();
      fListHist->Add(fHist_Reso2n_EP_Norm_Det[i][j]);
    }
    //----------- PID -------------------
    for(int j=0;j<4;j++){       //Detector: 0 = V0A, 1 = V0C, 3 = TPCA, 4 = TPCC 
      //----------> Pion:
      fHist_Corr3p_Pion_EP_Norm_PN[i][j] = new TProfile(Form("fHist_Corr3p_Pion_EP_Norm_PosNeg_Mag%d_Det%d",i,j+1),Form("US, #Psi_{2} %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_Pion_EP_Norm_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_Pion_EP_Norm_PN[i][j]);
      fHist_Corr3p_Pion_EP_Norm_PP[i][j] = new TProfile(Form("fHist_Corr3p_Pion_EP_Norm_PosPos_Mag%d_Det%d",i,j+1),Form("P-P, #Psi_{2} %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_Pion_EP_Norm_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_Pion_EP_Norm_PP[i][j]);
      fHist_Corr3p_Pion_EP_Norm_NN[i][j] = new TProfile(Form("fHist_Corr3p_Pion_EP_Norm_NegNeg_Mag%d_Det%d",i,j+1),Form("N-N, #Psi_{2}, %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_Pion_EP_Norm_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_Pion_EP_Norm_NN[i][j]);
      //----------> Kaon:
      fHist_Corr3p_Kaon_EP_Norm_PN[i][j] = new TProfile(Form("fHist_Corr3p_Kaon_EP_Norm_PosNeg_Mag%d_Det%d",i,j+1),Form("US, #Psi_{2} %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_Kaon_EP_Norm_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_Kaon_EP_Norm_PN[i][j]);
      fHist_Corr3p_Kaon_EP_Norm_PP[i][j] = new TProfile(Form("fHist_Corr3p_Kaon_EP_Norm_PosPos_Mag%d_Det%d",i,j+1),Form("P-P, #Psi_{2} %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_Kaon_EP_Norm_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_Kaon_EP_Norm_PP[i][j]);
      fHist_Corr3p_Kaon_EP_Norm_NN[i][j] = new TProfile(Form("fHist_Corr3p_Kaon_EP_Norm_NegNeg_Mag%d_Det%d",i,j+1),Form("N-N, #Psi_{2}, %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_Kaon_EP_Norm_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_Kaon_EP_Norm_NN[i][j]);
      //----------> Proton:
      fHist_Corr3p_Proton_EP_Norm_PN[i][j] = new TProfile(Form("fHist_Corr3p_Proton_EP_Norm_PosNeg_Mag%d_Det%d",i,j+1),Form("US, #Psi_{2} %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_Proton_EP_Norm_PN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_Proton_EP_Norm_PN[i][j]);
      fHist_Corr3p_Proton_EP_Norm_PP[i][j] = new TProfile(Form("fHist_Corr3p_Proton_EP_Norm_PosPos_Mag%d_Det%d",i,j+1),Form("P-P, #Psi_{2} %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_Proton_EP_Norm_PP[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_Proton_EP_Norm_PP[i][j]);
      fHist_Corr3p_Proton_EP_Norm_NN[i][j] = new TProfile(Form("fHist_Corr3p_Proton_EP_Norm_NegNeg_Mag%d_Det%d",i,j+1),Form("N-N, #Psi_{2}, %s",gDetForEP[j]),10,centRange,"");
      //fHist_Corr3p_Proton_EP_Norm_NN[i][j]->Sumw2();
      fListHist->Add(fHist_Corr3p_Proton_EP_Norm_NN[i][j]);
    }//Det loop
  }//magfield loop


  Double_t truncPi = 3.1416;

  fHV0AEventPlaneVsCent = new TH2F("fHV0AEventPlaneVsCent",Form("Psi %d from V0A", gPsiN), 10,centRange,50,-0.0,truncPi);
  fListHist->Add(fHV0AEventPlaneVsCent);
  
  fHV0CEventPlaneVsCent = new TH2F("fHV0CEventPlaneVsCent",Form("Psi %d from V0C", gPsiN), 10,centRange,50,-0.0,truncPi);
  fListHist->Add(fHV0CEventPlaneVsCent);

  fHTPCAEventPlaneVsCent = new TH2F("fHTPCAEventPlaneVsCent",Form("Psi %d from V0A",gPsiN),10,centRange,50,-0.0,truncPi);
  fListHist->Add(fHTPCAEventPlaneVsCent);

  fHTPCCEventPlaneVsCent = new TH2F("fHTPCCEventPlaneVsCent",Form("Psi %d from V0A",gPsiN),10,centRange,50,-0.0,truncPi);
  fListHist->Add(fHTPCCEventPlaneVsCent);
  //-------------------------------------------------------------------------------


  //---- to store NUA and calib histograms -----
  TList *fListNUACalib = new TList();
  fListNUACalib->SetName("fListNUACalib");
  fListNUACalib->SetOwner(kTRUE);


  //----------------- V0 Calibration hist: ---------------------
  fV0MultChVsRun = new TH2F("fV0MultChVsRun","1-32 V0C, 33-64 V0A",64,0,64,20,0,100);
  fListNUACalib->Add(fV0MultChVsRun);

  const char *sCorrect[2]={"w Corr","w/o Corr"};
  Int_t isCorr = 1;

  if(fListV0MCorr){
    isCorr = 0;
  }


  fV0AQ2xVsCentRun = new TProfile("fV0ACos2nVsCentRun",Form("<Cos2> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fV0AQ2xVsCentRun);
  fV0AQ2yVsCentRun = new TProfile("fV0ASin2nVsCentRun",Form("<Sin2> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fV0AQ2yVsCentRun);
  fV0CQ2xVsCentRun = new TProfile("fV0CCos2nVsCentRun",Form("<Cos2> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fV0CQ2xVsCentRun);
  fV0CQ2yVsCentRun = new TProfile("fV0CSin2nVsCentRun",Form("<Sin2> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fV0CQ2yVsCentRun);

  fV0AQ3xVsCentRun = new TProfile("fV0ACos3nVsCentRun",Form("<Cos3> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fV0AQ3xVsCentRun);
  fV0AQ3yVsCentRun = new TProfile("fV0ASin3nVsCentRun",Form("<Sin3> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fV0AQ3yVsCentRun);
  fV0CQ3xVsCentRun = new TProfile("fV0CCos3nVsCentRun",Form("<Cos3> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fV0CQ3xVsCentRun);
  fV0CQ3yVsCentRun = new TProfile("fV0CSin3nVsCentRun",Form("<Sin3> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fV0CQ3yVsCentRun);

  isCorr = 1;
  if(fListNUACorr){
    isCorr = 0;
  }
  //------------------- TPC Qvector Recentering Histograms --------------
  fTPCAQ2xVsCentRun = new TProfile("fTPCACos2nVsCentRun",Form("<Cos2> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCAQ2xVsCentRun);
  fTPCAQ2yVsCentRun = new TProfile("fTPCASin2nVsCentRun",Form("<Sin2> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCAQ2yVsCentRun);
  fTPCCQ2xVsCentRun = new TProfile("fTPCCCos2nVsCentRun",Form("<Cos2> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCCQ2xVsCentRun);
  fTPCCQ2yVsCentRun = new TProfile("fTPCCSin2nVsCentRun",Form("<Sin2> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCCQ2yVsCentRun);

  fTPCAQ3xVsCentRun = new TProfile("fTPCACos3nVsCentRun",Form("<Cos3> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCAQ3xVsCentRun);
  fTPCAQ3yVsCentRun = new TProfile("fTPCASin3nVsCentRun",Form("<Sin3> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCAQ3yVsCentRun);
  fTPCCQ3xVsCentRun = new TProfile("fTPCCCos3nVsCentRun",Form("<Cos3> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCCQ3xVsCentRun);
  fTPCCQ3yVsCentRun = new TProfile("fTPCCSin3nVsCentRun",Form("<Sin3> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCCQ3yVsCentRun);

  fTPCAQ4xVsCentRun = new TProfile("fTPCACos4nVsCentRun",Form("<Cos4> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCAQ4xVsCentRun);
  fTPCAQ4yVsCentRun = new TProfile("fTPCASin4nVsCentRun",Form("<Sin4> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCAQ4yVsCentRun);
  fTPCCQ4xVsCentRun = new TProfile("fTPCCCos4nVsCentRun",Form("<Cos4> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCCQ4xVsCentRun);
  fTPCCQ4yVsCentRun = new TProfile("fTPCCSin4nVsCentRun",Form("<Sin4> vs cent (%s)",sCorrect[isCorr]),90,0,90,"");
  fListNUACalib->Add(fTPCCQ4yVsCentRun);

  //-------------------------------------------------------------------------------

  
  
  
 //-------------------------- Define NUA Hist for PID -----------------------------
  Int_t gCentForNUA[6] = {0,5,10,20,40,90};
  Char_t  name[100];
  Char_t title[100];

  for(int i=0;i<4;i++){
    for(int j=0;j<5;j++){
      sprintf(name,"fHistEtaPhiVz_%s_Pos_Cent%d_Run%d",gSpecies[i],j,1);
      sprintf(title,"eta,phi,Vz %sPos, Cent%d-%d, FB %d",gSpecies[i],gCentForNUA[j],gCentForNUA[j+1],fFilterBit);
      fHist3DEtaPhiVz_Pos_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
      fListNUACalib->Add(fHist3DEtaPhiVz_Pos_Run[i][j]);

      sprintf(name,"fHistEtaPhiVz_%s_Neg_Cent%d_Run%d",gSpecies[i],j,1);
      sprintf(title,"eta,phi,Vz %sNeg, Cent%d-%d, FB %d",gSpecies[i],gCentForNUA[j],gCentForNUA[j+1],fFilterBit);
      fHist3DEtaPhiVz_Neg_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
      fListNUACalib->Add(fHist3DEtaPhiVz_Neg_Run[i][j]);
    }
  }
 //---------------------------------------------------------------------------------


  fListHist->Add(fListNUACalib);

  PostData(1,fListHist);
  std::cout<<"\n.........UserCreateOutputObject called.........\n fFilterBit = "<<fFilterBit<<" CentMax = "<<fCentralityPercentMax;
  std::cout<<" PU C = "<<fPileUpConstParm<<" PsiN = "<<gPsiN<<"\n\n"<<std::endl;
}

















//______________________________________________________________________
void AliAnalysisTaskCMEV0PID::UserExec(Option_t*){
  //printf("info: UserExec is called.... 1  \n");

  //if(fEventCount==1200) return;

  Float_t stepCount = 0.5;

  fHistEventCount->Fill(stepCount); //1
  stepCount++;

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    
  if(!(fESD || fAOD)){ printf("ERROR: fESD & fAOD not available\n"); return;  }
    
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());

  if (!fVevent) { printf("ERROR: fVevent not available\n");  return;  }
   
  fHistEventCount->Fill(stepCount); //2
  stepCount++;







  //--------- Check if I have PID response object --------
  if(!fPIDResponse){
    printf("\n\n...... PIDResponse object not found..... \n\n");
    return;
  }








  //-------------- Vtx cuts ---------------
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    
  Double_t pVtxZ = -999;
  pVtxZ = pVtx->GetZ();
    
  if(TMath::Abs(pVtxZ)>10.) return;

  fHistEventCount->Fill(stepCount); //3
  stepCount++;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");

  if(!fMultSelection) { printf("\n\n **WARNING** \n::UserExec() AliMultSelection object not found.\n\n"); exit(1); }

  Float_t centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  Float_t centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
//Float_t centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
//Float_t centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");

  Float_t centrality = centrV0M;

  if(centrality<fCentralityPercentMin || centrality>fCentralityPercentMax){ 
    return;
  }

  fHistEventCount->Fill(stepCount); //4
  stepCount++;

  fCentDistBefore->Fill(centrality);
  

  





  Int_t ntracks=fAOD->GetNumberOfTracks();
  if(ntracks<2) return;              // Check this cut....!!!
    
  fHistEventCount->Fill(stepCount); //5
  stepCount++;






  Int_t cent10bin = -1;

  cent10bin = GetCentralityScaled0to10(centrality); //Centrality in 0-10 scale


//Centrality array index for NUA correcion
  Int_t cForNUA = 0;  

  if(centrality<5.0) {
    cForNUA = 0;
  }
  else if(centrality>=5.0 && centrality<10){
    cForNUA = 1; // 1=5-10,
  }
  else if(centrality>=10.0 && centrality<20) {
    cForNUA = 2; // 2 = 10-20,
  }
  else if(centrality>=20 && centrality<40){ 
    cForNUA = 3; // 3=20-40
  }
  else if(centrality>=40){
    cForNUA = 4; // 4=40-90
  }








  //---------- Magnetic field --------
  Double_t fMagField = fAOD->GetMagneticField();

  const Int_t QAindex = (fMagField > 0) ? 1 : 0;
  //---------------------------------






  //Load NUA and V0M correction map run by run:
  Int_t runNumber = fAOD->GetRunNumber();
 
  if(runNumber!=fOldRunNum) {
   
    GetNUACorrectionHist(runNumber);
    
    if(bV0MGainCorr) {
      GetV0MCorrectionHist(runNumber);
    }

    fOldRunNum = runNumber;
  } 
  //------------------------------------------
 






  //----- Event Plane variables:-------
  Double_t PsiNV0A = 0;
  Double_t PsiNV0C = 0;

  Double_t PsiNTPCA = 0; // eta <0 
  Double_t PsiNTPCC = 0; // eta >0
  Double_t PsiNTPCF = 0; // Full TPC

  Double_t sumTPCQn2x[5] = {0,0,0};  //[0]= eta<0; [1]= eta>0; [2]= -0.8 < eta < 0.8
  Double_t sumTPCQn2y[5] = {0,0,0};
  Double_t sumTPCQn3x[5] = {0,0,0};  //[0]= eta<0; [1]= eta>0; [2]= -0.8 < eta < 0.8
  Double_t sumTPCQn3y[5] = {0,0,0};
  Double_t sumTPCQn4x[5] = {0,0,0};  //[0]= eta<0; [1]= eta>0; [2]= -0.8 < eta < 0.8
  Double_t sumTPCQn4y[5] = {0,0,0};
  //------------------------------------





 

 //Variables for MC tracking correction 
  Int_t   ptBinMC = 1;
  Int_t   iBinNUA = 1;
  Float_t ptWgtMC = 1.0;
  Float_t WgtNUA  = 1.0;
  Float_t ptTrk   = 0.1;


 //-------------- Track loop for outlier and PileUp cut -------------------

 //---------------- a dobrin --------------

  Bool_t bIsPileup=kFALSE;

  Int_t isPileup = fAOD->IsPileupFromSPD(3);

  if(isPileup != 0) {
    fHistPileUpCount->Fill(0.5);
    bIsPileup=kTRUE;          
  }
  else if(PileUpMultiVertex(fAOD)) {
    fHistPileUpCount->Fill(1.5);
    bIsPileup=kTRUE;
  }
  else if(((AliAODHeader*)fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) {
    fHistPileUpCount->Fill(2.5);
    bIsPileup=kTRUE;
  }
  else if(fAOD->IsIncompleteDAQ())  {
    fHistPileUpCount->Fill(3.5);
    bIsPileup=kTRUE;
  }
  else if(fabs(centrV0M-centrCL1)> 5.0)  {//default: 7.5
//else if(fabs(centrV0M-centrCL1)> 7.5)  {//default: 7.5
    fHistPileUpCount->Fill(4.5);
    bIsPileup=kTRUE;
  }

  // check vertex consistency
  const AliAODVertex* vtTrc = fAOD->GetPrimaryVertex();
  const AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();

  if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
    fHistPileUpCount->Fill(5.5);
    bIsPileup=kTRUE;
  }

  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);

  double dz = vtTrc->GetZ() - vtSPD->GetZ();

  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = dz/errTot;
  double nsigTrc = dz/errTrc;

  if(TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
    fHistPileUpCount->Fill(6.5);
    bIsPileup=kTRUE;
  }

  Float_t multTPC     = 0;    // tpc mult estimate
  Float_t RefMultRaw  = 0;    // tpc mult estimate
  Float_t RefMultCorr = 0;    // tpc mult estimate
  Float_t RefMultRawFB = 0;
  Float_t RefMultCorrFB= 0;

  Float_t multTPCAll  = 0;    // tpc mult estimate
  Float_t multGlobal  = 0; // global multiplicity
  Float_t etaTrk, phiTrk;
  Int_t   ChTrk;
  Float_t multEtaNeg, multEtaPos, multEtaFull;

  multEtaNeg=0;
  multEtaPos=0;
  multEtaFull=0;


//const Int_t nGoodTracks = fVevent->GetNumberOfTracks();
  const Int_t nGoodTracks = fAOD->GetNumberOfTracks();

  for(Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++) { //-------------------------

    AliAODTrack* AODtrack =dynamic_cast<AliAODTrack*>(fVevent->GetTrack(iTrack));
    if(!AODtrack) continue;

    if(AODtrack->TestFilterBit(128)) multTPCAll++; //A. Dobrin, no track cuts

    ptTrk  = AODtrack->Pt();
    etaTrk = AODtrack->Eta();
    phiTrk = AODtrack->Phi();
    ChTrk  = AODtrack->Charge();
    
    if(!(AODtrack->TestFilterBit(1))) continue;
    if((ptTrk < .2) || (TMath::Abs(etaTrk) > .8) || (AODtrack->GetTPCNcls() < 70)  || (AODtrack->GetDetPid()->GetTPCsignal() < 10.0) || (AODtrack->Chi2perNDF() < 0.2)) continue;
    multTPC++;



    if(ptTrk >= 0.2 && ptTrk <= 5.0) {
      if(fFB_Efficiency_Cent[cent10bin]){
	ptBinMC = fFB_Efficiency_Cent[cent10bin]->FindBin(ptTrk);
	ptWgtMC = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBinMC);
        //cout<<iTrack<<"cent = "<<cent10bin<<" pt = "<<ptTrk<<"\t bin = "<<ptBinMC<<"\t Wgt = "<<ptWgtMC<<endl;    
	RefMultRaw++;
	RefMultCorr += ptWgtMC;
	
        if((AODtrack->TestFilterBit(fFilterBit))){
	  
	  RefMultRawFB++;
	  RefMultCorrFB += ptWgtMC;

	  if(ptTrk <= 5.0){
	    //Get NUA weights for EP:
	    if(ChTrk>0){
	      if(fHCorrectNUApos[cForNUA]){
		iBinNUA = fHCorrectNUApos[cForNUA]->FindBin(pVtxZ,phiTrk,etaTrk);
		WgtNUA  = fHCorrectNUApos[cForNUA]->GetBinContent(iBinNUA);
	      }
	      else{ WgtNUA = 1.0; }
	    }
	    else{
	      if(fHCorrectNUAneg[cForNUA]){
		iBinNUA = fHCorrectNUAneg[cForNUA]->FindBin(pVtxZ,phiTrk,etaTrk);
		WgtNUA  = fHCorrectNUAneg[cForNUA]->GetBinContent(iBinNUA);  
	      }
	      else{ WgtNUA = 1.0; }
	    }
	    
	    if(etaTrk < -0.05){
	      sumTPCQn2x[0] += WgtNUA*TMath::Cos(gPsiN*phiTrk);
	      sumTPCQn2y[0] += WgtNUA*TMath::Sin(gPsiN*phiTrk);
	      sumTPCQn3x[0] += WgtNUA*TMath::Cos(3*phiTrk);
	      sumTPCQn3y[0] += WgtNUA*TMath::Sin(3*phiTrk);
	      sumTPCQn4x[0] += WgtNUA*TMath::Cos(4*phiTrk);
	      sumTPCQn4y[0] += WgtNUA*TMath::Sin(4*phiTrk);
	      multEtaNeg++;
	    }
	    else if(etaTrk > 0.05){
	      sumTPCQn2x[1] += WgtNUA*TMath::Cos(gPsiN*phiTrk);
	      sumTPCQn2y[1] += WgtNUA*TMath::Sin(gPsiN*phiTrk);
	      sumTPCQn3x[1] += WgtNUA*TMath::Cos(3*phiTrk);
	      sumTPCQn3y[1] += WgtNUA*TMath::Sin(3*phiTrk);
	      sumTPCQn4x[1] += WgtNUA*TMath::Cos(4*phiTrk);
	      sumTPCQn4y[1] += WgtNUA*TMath::Sin(4*phiTrk);
	      multEtaPos++;
	    }
	    sumTPCQn2x[3] += WgtNUA*TMath::Cos(gPsiN*phiTrk);
	    sumTPCQn2y[3] += WgtNUA*TMath::Sin(gPsiN*phiTrk);
	    sumTPCQn3x[3] += WgtNUA*TMath::Cos(3*phiTrk);
	    sumTPCQn3y[3] += WgtNUA*TMath::Sin(3*phiTrk);
	    sumTPCQn4x[3] += WgtNUA*TMath::Cos(4*phiTrk);
	    sumTPCQn4y[3] += WgtNUA*TMath::Sin(4*phiTrk);
	    multEtaFull++;
	  }
	}
      }
    }  

    Double_t b[2] = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
    AliAODTrack copy(*AODtrack);
    if (!(copy.PropagateToDCA(fVevent->GetPrimaryVertex(), fVevent->GetMagneticField(), 100., b, bCov))) continue;
    if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3)) continue;
    multGlobal++;

  }//--- track loop outlier/PileUp ----

  Int_t multEsd = ((AliAODHeader*)fAOD->GetHeader())->GetNumberOfESDTracks();

  Float_t multESDTPCDiff = (Float_t) multEsd - fPileUpSlopeParm*multTPCAll;

  if(multESDTPCDiff > fPileUpConstParm) { 
    fHistPileUpCount->Fill(7.5);
    bIsPileup=kTRUE;
  }
  else if(bIsPileup==kFALSE) {
    if(!fMultSelection->GetThisEventIsNotPileup()){
      fHistMultSelPUCount->Fill(0.5);
      bIsPileup=kTRUE;
    }
    if(!fMultSelection->GetThisEventIsNotPileupMV()){
      fHistMultSelPUCount->Fill(1.5);
      bIsPileup=kTRUE;
    }
    if(!fMultSelection->GetThisEventIsNotPileupInMultBins()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    }
    if(!fMultSelection->GetThisEventHasNoInconsistentVertices()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    }      
    if(!fMultSelection->GetThisEventPassesTrackletVsCluster()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    } 
    if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    }      
    if(!fMultSelection->GetThisEventHasGoodVertex2016()){
      fHistMultSelPUCount->Fill(2.5);
      bIsPileup=kTRUE;
    }      
    if(bIsPileup) fHistPileUpCount->Fill(9.5);
  }  
  //-----------------------------------------------------------------



  fHistTPCVsESDTrkBefore->Fill(multTPCAll,multEsd);   //A. Dobrin

  fHistTPCvsGlobalMultBefore->Fill(multGlobal,multTPC);
     

  Bool_t  bIsOutLier=kFALSE;
  
  if(multTPC < (-20.0+1.15*multGlobal) || multTPC > (200.+1.40*multGlobal)) { bIsOutLier = kTRUE;}
       
  fHistEventCount->Fill(stepCount); //6
  stepCount++;


  fHistTPConlyVsCL1Before->Fill(centrCL1,multTPCAll);
  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTPCAll);
  fHistGlobalVsV0MBefore->Fill(centrV0M, multGlobal);



  //if bSkipPileUpCut is kTRUE then don't apply PileUp removal.
  if(!bSkipPileUpCut && bIsOutLier) return; //outlier TPC vs Global

  fHistTPCvsGlobalMultAfter->Fill(multGlobal,multTPC);

  fHistEventCount->Fill(stepCount); //7
  stepCount++;

  //cout<<"After Outlier  = "<<fSkipOutlierCut<<" multTPC = "<<multTPC<<" multGlobal = "<<multGlobal<<endl;


  if(!bSkipPileUpCut && bIsPileup) return;         //PileUp A. Dobrin

  fHistTPCVsESDTrkAfter->Fill(multTPCAll,multEsd);  

  fHistEventCount->Fill(stepCount); //8
  stepCount++;

  //cout<<"After PU cut = "<<fSkipOutlierCut<<" multTPC = "<<multTPC<<" multGlobal = "<<multGlobal<<endl;









  Int_t icentBin = centrality;
  icentBin++;

  Float_t TPCmultLowLimit  =  hCentvsTPCmultCuts->GetBinContent(icentBin,1);
  Float_t TPCmultHighLimit =  hCentvsTPCmultCuts->GetBinContent(icentBin,1);

  TPCmultLowLimit  -=  5.0 * hCentvsTPCmultCuts->GetBinContent(icentBin,2); //mean - 5sigma
  TPCmultHighLimit +=  5.0 * hCentvsTPCmultCuts->GetBinContent(icentBin,2); //mean + 5sigma
  //std::cout<<" Cent = "<<centrality<<"\t icent = "<<icentBin<<" low = "<<TPCmultLowLimit<<"\t high = "<<TPCmultHighLimit<<std::endl;

  if(!bSkipPileUpCut){
    if(multTPC<TPCmultLowLimit || multTPC>TPCmultHighLimit) return; //centrality outlier
  }

  fHistEventCount->Fill(stepCount); //9
  stepCount++;



  fHistTPConlyVsCL1After->Fill(centrCL1,multTPCAll);
  fHistTPConlyVsV0MAfter->Fill(centrV0M,multTPCAll);
  fHistGlobalVsV0MAfter->Fill(centrV0M, multGlobal);


  // MC corrected Refmult:
  fHistRawVsCorrMultFB->Fill(RefMultRawFB,RefMultCorrFB); // FB set by AddTask.. 


  Float_t EvtCent = centrality;


  if(multEtaNeg<2 || multEtaPos<2) return;  //Minimum 2 tracks in each eta

  fHistEventCount->Fill(stepCount); //10
  stepCount++;

  //--------------------------------------------------------


















  //--------------- cent CL1 <= 90 cut ------------------
  Int_t icentV0Qn = centrCL1; // cent CL1 used for V0 calibration.
  icentV0Qn += 1;

  if(icentV0Qn>90)     return;
  
  fHistEventCount->Fill(stepCount); //11
  stepCount++;


  

  //-------- V0M info ---------------
  const AliAODVZERO *fAODV0 = fAOD->GetVZEROData();

  //do v0m recentering
  Double_t QxanCor = 0, QyanCor = 0;
  Double_t QxcnCor = 0, QycnCor = 0;
   
  Double_t Qxan3  = 0., Qyan3 = 0.;
  Double_t Qxcn3  = 0., Qycn3 = 0.;
  Double_t Qxan2  = 0., Qyan2 = 0.;
  Double_t Qxcn2  = 0., Qycn2 = 0.;

  Float_t fMultv0, phiV0;
  Float_t sumMa=0;
  Float_t sumMc=0;

  for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA

    fMultv0 = fAODV0->GetMultiplicity(iV0);
 
    if(fHCorrectV0M){
      fMultv0 = fMultv0 * fHCorrectV0M->GetBinContent(iV0+1);  // Gain Correction
    //cout<<"info: run = "<<runNumber<<" cent = "<<centrCL1<<"\t channel = "<<iV0<<" gain = "<<fHCorrectV0M->GetBinContent(iV0+1)<<endl;
    }
     
    fV0MultChVsRun->Fill(iV0+0.5,centrCL1,fMultv0);
  //fV0MultChVsRun->Fill(iV0+0.5,runindex,fMultv0);

    phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);
 
    if(iV0 < 32){
      Qxcn2  += TMath::Cos(2*phiV0) * fMultv0;
      Qycn2  += TMath::Sin(2*phiV0) * fMultv0;
      Qxcn3  += TMath::Cos(3*phiV0) * fMultv0;
      Qycn3  += TMath::Sin(3*phiV0) * fMultv0;
      sumMc += fMultv0;
    }
    else if(iV0 >= 32){
      Qxan2  += TMath::Cos(2*phiV0) * fMultv0;
      Qyan2  += TMath::Sin(2*phiV0) * fMultv0;
      Qxan3  += TMath::Cos(3*phiV0) * fMultv0;
      Qyan3  += TMath::Sin(3*phiV0) * fMultv0;
      sumMa += fMultv0;
    }
  }//----- channel loop ----------

   
  if(gPsiN==3){
    QxanCor = Qxan3/sumMa;  //3rd order event plane
    QyanCor = Qyan3/sumMa;
    QxcnCor = Qxcn3/sumMc;
    QycnCor = Qycn3/sumMc;

    if(fHAvgerageQnV0C && fHAvgerageQnV0A && icentV0Qn < 91){
      QxanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,3);  //x = Cos
      QxcnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,3);  //x = Cos
      QyanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,4);  //y = Sin
      QycnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,4);  //y = Sin
    }
    //printf("\n .... I am using my own V0 gain correction for Psi3...\n");
  }
  else{
    QxanCor = Qxan2/sumMa;  //2nd order Event plane 
    QyanCor = Qyan2/sumMa;
    QxcnCor = Qxcn2/sumMc;
    QycnCor = Qycn2/sumMc;

    if(fHAvgerageQnV0C && fHAvgerageQnV0A && icentV0Qn < 91){
      QxanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,1);  //x = Cos
      QxcnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,1);  //x = Cos
      QyanCor -= fHAvgerageQnV0A->GetBinContent(icentV0Qn,2);  //y = Sin
      QycnCor -= fHAvgerageQnV0C->GetBinContent(icentV0Qn,2);  //y = Sin  
    }
    //printf("\n .... I am using my own V0 gain correction for Psi2...\n");
  }
   

 //------------ For V0-Qn Recenter and Event plane ----------
  fV0CQ2xVsCentRun->Fill(centrCL1,Qxcn2/sumMc);
  fV0CQ2yVsCentRun->Fill(centrCL1,Qycn2/sumMc);
  fV0AQ2xVsCentRun->Fill(centrCL1,Qxan2/sumMa); 
  fV0AQ2yVsCentRun->Fill(centrCL1,Qyan2/sumMa); 

  fV0CQ3xVsCentRun->Fill(centrCL1,Qxcn3/sumMc);
  fV0CQ3yVsCentRun->Fill(centrCL1,Qycn3/sumMc);
  fV0AQ3xVsCentRun->Fill(centrCL1,Qxan3/sumMa); 
  fV0AQ3yVsCentRun->Fill(centrCL1,Qyan3/sumMa); 
  
  PsiNV0C = 1.0/gPsiN*( TMath::ATan2(QycnCor,QxcnCor) ); //+ TMath::Pi() );
  if(PsiNV0C<0.) PsiNV0C += 2*TMath::Pi()/gPsiN;

  PsiNV0A = 1.0/gPsiN*( TMath::ATan2(QyanCor,QxanCor) ); //+ TMath::Pi() );
  if(PsiNV0A<0.) PsiNV0A += 2*TMath::Pi()/gPsiN;

  //fHV0CEventPlaneVsCent->Fill(EvtCent,TMath::Cos(PsiNV0C-PsiNV0A));
  fHV0CEventPlaneVsCent->Fill(EvtCent,PsiNV0C);
  fHV0AEventPlaneVsCent->Fill(EvtCent,PsiNV0A);
  //-------------------------------------------------------------












  //----- Psi_N from TPC sub & full event plane -------
  /*
  PsiNTPCA = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[0],sumTPCQn2x[0] )) ; // negetive eta
  if(PsiNTPCA<0.) PsiNTPCA += 2*TMath::Pi()/gPsiN;
  PsiNTPCC = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[1],sumTPCQn2x[1] )) ; // positive eta
  if(PsiNTPCC<0.) PsiNTPCC += 2*TMath::Pi()/gPsiN;  */

  // Enable periodicity PsiN directly:
  PsiNTPCA = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[0],sumTPCQn2x[0]) + TMath::Pi() ) ; // negetive eta
  //if(PsiNTPCA<0.) PsiNTPCA += 2*TMath::Pi()/gPsiN;
  PsiNTPCC = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[1],sumTPCQn2x[1]) + TMath::Pi() ) ; // positive eta
  //if(PsiNTPCC<0.) PsiNTPCC += 2*TMath::Pi()/gPsiN;


  //fHTPCAEventPlaneVsCent->Fill(EvtCent,TMath::Cos(PsiNTPCA-PsiNTPCC));
  //fHTPCCEventPlaneVsCent->Fill(EvtCent,TMath::Cos(PsiNTPCC-PsiNV0C));

  //fHTPCAEventPlaneVsCent->Fill(EvtCent,PsiNTPCA);
  //fHTPCCEventPlaneVsCent->Fill(EvtCent,PsiNTPCC);

  fEventCount++;



  PsiNTPCF = (1.0/gPsiN)*( TMath::ATan2(sumTPCQn2y[3],sumTPCQn2x[3] )) ; // FUll TPC eta 
  if(PsiNTPCF<0.) PsiNTPCF += 2*TMath::Pi()/gPsiN;

  //Psi3TPCA = (1.0/3)*( TMath::ATan2(sumTPCQn3y[0],sumTPCQn3x[0] )) ; // negetive eta ; Psi3 not needed now.
  //if(Psi3TPCA<0.) Psi3TPCA += 2*TMath::Pi()/3;
  //Psi3TPCC = (1.0/3)*( TMath::ATan2(sumTPCQn3y[1],sumTPCQn3x[1] )) ; // positive eta 
  //if(Psi3TPCC<0.) Psi3TPCC += 2*TMath::Pi()/3;
  //----------------------------------------------------------------











 //---- Copies of TPC-Q vectors to remove track -by- track auto-correlation -----

  Double_t sumQxTPCneg;
  Double_t sumQyTPCneg;
  Double_t sumQxTPCpos;
  Double_t sumQyTPCpos;
  Double_t sumQxTPCneg2;
  Double_t sumQyTPCneg2;
  Double_t sumQxTPCpos2;
  Double_t sumQyTPCpos2;


// Fill Centrality for run-by-run Event weight correction. 
  fCentDistAfter->Fill(centrality); 
  





  //--------- Track variable for PID/Charge studies ----------------
  Double_t   PDGmassPion   = 0.13957;
  Double_t   PDGmassKaon   = 0.49368;
  Double_t   PDGmassProton = 0.93827;

  PDGmassProton *= PDGmassProton;
  PDGmassPion   *= PDGmassPion;
  PDGmassKaon   *= PDGmassKaon;
 
  Double_t mass=0,mom = -999, pT = -999, phi = -999, eta = -999, dEdx =-999;
  Double_t length = -999., beta =-999, tofTime = -999., tof = -999.;
  Double_t dPhi1,dPhi2,dPt1,dPt2,dEta1,dEta2;
  Double_t ptw1, ptw2, w1NUA, w2NUA;
  Double_t nSigTOFpion,  nSigTPCpion;
  Double_t nSigTOFkaon,  nSigTPCkaon;
  Double_t nSigTOFproton,nSigTPCproton;
  Double_t nSigTOFpion2,  nSigTPCpion2;
  Double_t nSigTOFkaon2,  nSigTPCkaon2;
  Double_t nSigTOFproton2,nSigTPCproton2;

  Double_t c = TMath::C()*1.E-9;// velocity of light m/ns 
  Double_t  dcaXY, dcaZ, WgtEP ;
  Double_t  probMis; 
  Int_t    TOFmatch=0; 
  Int_t    charge,ptBin,dChrg1,dChrg2;



  //----------- Set the desired Harmonic ------------
  Int_t n = gN;
  Int_t m = gM;
  Int_t p =n+m;
  //------------------------------------------------


  Int_t skipPairHBT = 0;






  Double_t ptwPion1,ptwKaon1,ptwProton1; 
  Double_t ptwPion2,ptwKaon2,ptwProton2;
  Double_t wNUAPion1,wNUAKaon1,wNUAProton1;
  Double_t wNUAPion2,wNUAKaon2,wNUAProton2;

  Double_t    WgtEPPion = 1.0;
  Double_t    WgtEPKaon = 1.0;
  Double_t    WgtEPProton = 1.0;


  Bool_t isPion1 = kFALSE;
  Bool_t isKaon1  = kFALSE;
  Bool_t isProton1 = kFALSE;
  Bool_t isPion2 = kFALSE;
  Bool_t isKaon2  = kFALSE;
  Bool_t isProton2 = kFALSE;






  for(Int_t itrack = 0; itrack < ntracks; itrack++) {

    AliAODTrack *track=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(itrack));
    if(!track) continue;

    dcaXY = track->DCA();
    dcaZ  = track->ZAtDCA();

    if(!track->TestFilterBit(fFilterBit)) continue;
        
    mom=track->P();
    pT=track->Pt();
    phi=track->Phi();
    eta=track->Eta();
    dEdx=track->GetDetPid()->GetTPCsignal();
    charge = track->Charge();
    



    
    //-------------- Check TOF status ------------------
    AliPIDResponse::EDetPidStatus status;
    status = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
    TOFmatch = 0;
    if(status==AliPIDResponse::kDetPidOk){
      TOFmatch++;
    }
    probMis = fPIDResponse->GetTOFMismatchProbability(track);
    fHistTOFMatchCount->Fill(TOFmatch,probMis);

    mass = -9.9;
    beta = -0.5;

    if(TOFmatch>0 && probMis < 0.01) { 
      //This conditions are called when detector status is checked above :  
      //if((track->IsOn(AliAODTrack::kTOFin)) && (track->IsOn(AliAODTrack::kTIME)) &&  (track->IsOn(AliAODTrack::kTOFout))) {
      //if((track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid))) { //Naghmeh used it        
	tofTime = track->GetTOFsignal();  // in pico seconds
	length  = track->GetIntegratedLength();   
	tof = tofTime*1.E-3; // ns
	if (tof <= 0) tof = 9999;            
	length = length*0.01; // in meters
	tof = tof*c;
	beta = length/tof;
	mass = mom*mom*(1./(beta*beta) - 1);             
    }//------------ TOF signal -------------------------

    
    //QA histograms:
    fHistEtaPtBefore->Fill(eta,pT);
    fHistTPCdEdxvsPBefore->Fill(mom*charge,dEdx);
    fHistTOFBetavsPBefore->Fill(mom*charge,beta);
    fHistTOFMassvsPtBefore->Fill(pT*charge,mass);



    
    //-------- Apply Default track cuts for analysis: ---------
    if(pT < fMinPtCut  ||  pT > fMaxPtCut)    continue;
    if(eta <fMinEtaCut || eta > fMaxEtaCut)   continue;
    //if(!(track->TestFilterBit(fFilterBit)))   continue;
    //-----------------------------------------------------





 //============= charged hadron analysis: ============

    dPt1  =   pT;
    dPhi1 =  phi;
    dEta1 =  eta;
    dChrg1=charge;

      
    sumQxTPCneg = sumTPCQn2x[0];   // first copy to remove 1st track fron Q vector (AutoCorrelation)
    sumQyTPCneg = sumTPCQn2y[0];
    sumQxTPCpos = sumTPCQn2x[1];
    sumQyTPCpos = sumTPCQn2y[1];




    
    //--------------------- PID signals 1st track-------------------------
    nSigTOFpion   = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    nSigTOFkaon   = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    nSigTOFproton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

    nSigTPCpion   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    nSigTPCkaon   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    nSigTPCproton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);


    //cout<<"Trk "<<itrack<<" pt1 = "<<dPt1<<"\tnSigPion = "<<nSigTPCpion<<"\tnSigKaon = "<<nSigTPCkaon <<"\tnSigprot = "<<nSigTPCproton<<endl;

    isPion1 = kFALSE;
    isKaon1  = kFALSE;
    isProton1 = kFALSE;

    //------> Pion
    if(dPt1<0.6 && TMath::Abs(nSigTPCpion)<=2.5){
      isPion1 = kTRUE;
    }
    else if(dPt1>=0.6 && dPt1<=2.0 && TMath::Abs(nSigTPCpion)<=2.5 && TMath::Abs(nSigTOFpion)<=2.0 ){
      isPion1 = kTRUE;
    }
    //------> Kaon
    if(dPt1<0.6 && TMath::Abs(nSigTPCkaon)<=2.5){
      isKaon1 = kTRUE;
    }
    else if(dPt1>=0.6 && dPt1<=2.0 && TMath::Abs(nSigTPCkaon)<=2.5 && TMath::Abs(nSigTOFkaon)<=2.0){
      isKaon1 = kTRUE;
    }
    //------> Proton 
    if(dPt1<0.8 && TMath::Abs(nSigTPCproton)<=2.5){
      isProton1 = kTRUE;
    }
    else if(dPt1>=0.8 && dPt1<=3.5 && TMath::Abs(nSigTPCproton)<=2.5 && TMath::Abs(nSigTOFproton)<=2.5){
      isProton1 = kTRUE;
    }

    //-----------------------------------------------------------------


    //=================  MC wgt and NUA wgt for PID =================
    ptwPion1 = 1.0;
    ptwKaon1  = 1.0;
    ptwProton1 = 1.0;
    wNUAPion1 = 1.0;
    wNUAKaon1  = 1.0;
    wNUAProton1 = 1.0;

    //------ get MC weight and NUA for Pion track1 --------------
    if(isPion1){
      if(fFB_Efficiency_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Pion Efficiency file when available.
	ptBin    = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt1);
	ptwPion1 = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
      }
      else{ ptwPion1 = 1.0; }

      if(dChrg1>0){
	if(fHCorrectNUAposPion[cForNUA]){
	  iBinNUA   = fHCorrectNUAposPion[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAPion1 = fHCorrectNUAposPion[cForNUA]->GetBinContent(iBinNUA);
	}
	else{ wNUAPion1 = 1.0; }
      }
      else{
	if(fHCorrectNUAnegPion[cForNUA]){
	  iBinNUA   = fHCorrectNUAnegPion[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAPion1 = fHCorrectNUAnegPion[cForNUA]->GetBinContent(iBinNUA);  
	}
	else{ wNUAPion1 = 1.0; }
      }
    }


    //------ get MC weight and NUA for Kaon track1 --------------
    if(isKaon1){
      if(fFB_Efficiency_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Kaon Efficiency file when available.
	ptBin    = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt1);
	ptwKaon1 = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
      }
      else{ ptwKaon1 = 1.0; }

      if(dChrg1>0){
	if(fHCorrectNUAposKaon[cForNUA]){
	  iBinNUA   = fHCorrectNUAposKaon[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAKaon1 = fHCorrectNUAposKaon[cForNUA]->GetBinContent(iBinNUA);
	}
	else{ wNUAKaon1 = 1.0; }
      }
      else{
	if(fHCorrectNUAnegKaon[cForNUA]){
	  iBinNUA   = fHCorrectNUAnegKaon[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAKaon1 = fHCorrectNUAnegKaon[cForNUA]->GetBinContent(iBinNUA);  
	}
	else{ wNUAKaon1 = 1.0; }
      }
    }
    //------ get MC weight and NUA for Proton track1 --------------
    if(isProton1){
      if(fFB_Efficiency_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Proton Efficiency file when available.
	ptBin    = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt1);
	ptwProton1 = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
      }
      else{ ptwProton1 = 1.0; }

      if(dChrg1>0){
	if(fHCorrectNUAposProton[cForNUA]){
	  iBinNUA   = fHCorrectNUAposProton[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAProton1 = fHCorrectNUAposProton[cForNUA]->GetBinContent(iBinNUA);
	}
	else{ wNUAProton1 = 1.0; }
      }
      else{
	if(fHCorrectNUAnegProton[cForNUA]){
	  iBinNUA   = fHCorrectNUAnegProton[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	  wNUAProton1 = fHCorrectNUAnegProton[cForNUA]->GetBinContent(iBinNUA);  
	}
	else{ wNUAProton1 = 1.0; }
      }
    }
    //=========================== X ===============================



      
    //------ get MC weight and NUA for Charged  track 1--------------
    ptw1  = 1.0;

    if(fFB_Efficiency_Cent[cent10bin]){
      ptBin = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt1);
      ptw1  = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
    }
    else{ ptw1 = 1.0; }

 
    if(dChrg1>0){
      if(fHCorrectNUApos[cForNUA]){
	iBinNUA = fHCorrectNUApos[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
	w1NUA = fHCorrectNUApos[cForNUA]->GetBinContent(iBinNUA);
      }
      else{ w1NUA = 1.0; }
    }
    else{
      if(fHCorrectNUAneg[cForNUA]){
        iBinNUA = fHCorrectNUAneg[cForNUA]->FindBin(pVtxZ,dPhi1,dEta1);
        w1NUA = fHCorrectNUAneg[cForNUA]->GetBinContent(iBinNUA);  
      }
      else{ w1NUA = 1.0; }
    }



    //-------- Remove track 1 from EP calculation ----------
    if(dEta1 < -0.05){
      sumQxTPCneg -= w1NUA*TMath::Cos(gPsiN*dPhi1);
      sumQyTPCneg -= w1NUA*TMath::Sin(gPsiN*dPhi1); // [0] = eta <-0.05
    }
    else if(dEta1 > 0.05){
      sumQxTPCpos -= w1NUA*TMath::Cos(gPsiN*dPhi1);
      sumQyTPCpos -= w1NUA*TMath::Sin(gPsiN*dPhi1); // [1] = eta > 0.05
    }
    //-----------------------------------------------------






   
    //---2nd track loop (nested)---
    for(Int_t jtrack = 0; jtrack < ntracks; jtrack++) {

      if(jtrack==itrack) continue;
    
      AliAODTrack *track2=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(jtrack));
      if(!track2) continue;
      
      dPt2  = track2->Pt();
      dPhi2 = track2->Phi();
      dEta2 = track2->Eta();
      dChrg2= track2->Charge();



      //------- Apply Default track cuts for analysis: -------
      if(dPt2 < fMinPtCut  ||  dPt2 > fMaxPtCut)     continue;
      if(dEta2 <fMinEtaCut || dEta2 > fMaxEtaCut)    continue;
      if(!(track2->TestFilterBit(fFilterBit)))       continue;
      //-----------------------------------------------------
			  
      //cout<<"info: passes 1 ";


      //--------------------- PID signals 2nd track-------------------------
      nSigTOFpion2   = fPIDResponse->NumberOfSigmasTOF(track2, AliPID::kPion);
      nSigTOFkaon2   = fPIDResponse->NumberOfSigmasTOF(track2, AliPID::kKaon);
      nSigTOFproton2 = fPIDResponse->NumberOfSigmasTOF(track2, AliPID::kProton);

      nSigTPCpion2   = fPIDResponse->NumberOfSigmasTPC(track2, AliPID::kPion);
      nSigTPCkaon2   = fPIDResponse->NumberOfSigmasTPC(track2, AliPID::kKaon);
      nSigTPCproton2 = fPIDResponse->NumberOfSigmasTPC(track2, AliPID::kProton);

      isPion2 = kFALSE;
      isKaon2  = kFALSE;
      isProton2 = kFALSE;

      //------> Pion 2
      if(dPt2<0.6 && TMath::Abs(nSigTPCpion2)<=2.5){
	isPion2 = kTRUE;
      }
      else if(dPt2>=0.6 && dPt2<=2.0 && TMath::Abs(nSigTPCpion2)<=2.5 && TMath::Abs(nSigTOFpion2)<=2.0 ){
	isPion2 = kTRUE;
      }
      //------> Kaon 2
      if(dPt2<0.6 && TMath::Abs(nSigTPCkaon2)<=2.5){
	isKaon2 = kTRUE;
      }
      else if(dPt2>=0.6 && dPt2<=2.0 && TMath::Abs(nSigTPCkaon2)<=2.5 && TMath::Abs(nSigTOFkaon2)<=2.0){
	isKaon2 = kTRUE;
      }
      //------> Proton 2
      if(dPt2<0.8 && TMath::Abs(nSigTPCproton2)<=2.5){
	isProton2 = kTRUE;
      }
      else if(dPt2>=0.8 && dPt2<=3.5 && TMath::Abs(nSigTPCproton2)<=2.5 && TMath::Abs(nSigTOFproton2)<=2.5){
	isProton2 = kTRUE;
      }
      //----------------------------------------------------------------


 


      //=================  MC wgt and NUA wgt for PID =================
      WgtEPPion = 1.0;
      WgtEPKaon = 1.0;
      WgtEPProton = 1.0;

      ptwPion2 = 1.0;
      ptwKaon2  = 1.0;
      ptwProton2 = 1.0;
      wNUAPion2 = 1.0;
      wNUAKaon2  = 1.0;
      wNUAProton2 = 1.0;

      //------ get MC weight and NUA for Pion track2 --------------
      if(isPion2){
	if(fFB_Efficiency_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Pion Efficiency file when available.
	  ptBin    = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt2);
	  ptwPion2 = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
	}
	else{ ptwPion2 = 1.0; }

	if(dChrg2>0){
	  if(fHCorrectNUAposPion[cForNUA]){
	    iBinNUA   = fHCorrectNUAposPion[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAPion2 = fHCorrectNUAposPion[cForNUA]->GetBinContent(iBinNUA);
	  }
	  else{ wNUAPion2 = 1.0; }
	}
	else{
	  if(fHCorrectNUAnegPion[cForNUA]){
	    iBinNUA   = fHCorrectNUAnegPion[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAPion2 = fHCorrectNUAnegPion[cForNUA]->GetBinContent(iBinNUA);  
	  }
	  else{ wNUAPion2 = 1.0; }
	}

	WgtEPPion = ptwPion1*ptwPion2*wNUAPion1*wNUAPion2;
      }

      //------ get MC weight and NUA for Kaon track2 --------------
      if(isKaon2){
	if(fFB_Efficiency_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Kaon Efficiency file when available.
	  ptBin    = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt2);
	  ptwKaon2 = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
	}
	else{ ptwKaon2 = 1.0; }

	if(dChrg2>0){
	  if(fHCorrectNUAposKaon[cForNUA]){
	    iBinNUA   = fHCorrectNUAposKaon[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAKaon2 = fHCorrectNUAposKaon[cForNUA]->GetBinContent(iBinNUA);
	  }
	  else{ wNUAKaon2 = 1.0; }
	}
	else{
	  if(fHCorrectNUAnegKaon[cForNUA]){
	    iBinNUA   = fHCorrectNUAnegKaon[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAKaon2 = fHCorrectNUAnegKaon[cForNUA]->GetBinContent(iBinNUA);  
	  }
	  else{ wNUAKaon2 = 1.0; }
	}

	WgtEPKaon = ptwKaon1*ptwKaon2*wNUAKaon1*wNUAKaon2;
      }
      //------ get MC weight and NUA for Proton track2 --------------
      if(isProton2){
	if(fFB_Efficiency_Cent[cent10bin]){ // <-------------------- !!!! WARNING: use Proton Efficiency file when available.
	  ptBin    = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt2);
	  ptwProton2 = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
	}
	else{ ptwProton2 = 1.0; }

	if(dChrg2>0){
	  if(fHCorrectNUAposProton[cForNUA]){
	    iBinNUA   = fHCorrectNUAposProton[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAProton2 = fHCorrectNUAposProton[cForNUA]->GetBinContent(iBinNUA);
	  }
	  else{ wNUAProton2 = 1.0; }
	}
	else{
	  if(fHCorrectNUAnegProton[cForNUA]){
	    iBinNUA   = fHCorrectNUAnegProton[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	    wNUAProton2 = fHCorrectNUAnegProton[cForNUA]->GetBinContent(iBinNUA);  
	  }
	  else{ wNUAProton2 = 1.0; }
	}

	WgtEPProton = ptwProton1*ptwProton2*wNUAProton1*wNUAProton2;
      }
      //========================== X ================================




 
      //------ get MC weight and NUA (Charged) for track 2--------------
      WgtEP = 1.0;
      ptw2  = 1.0;

      if(fFB_Efficiency_Cent[cent10bin]){
	ptBin = fFB_Efficiency_Cent[cent10bin]->FindBin(dPt2);
	ptw2  = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBin);
      }
      else{ ptw2 = 1.0; }

      if(dChrg2>0){
	if(fHCorrectNUApos[cForNUA]){
	  iBinNUA = fHCorrectNUApos[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	  w2NUA   = fHCorrectNUApos[cForNUA]->GetBinContent(iBinNUA);
	}
	else{ w2NUA = 1.0; }
      }
      else{
	if(fHCorrectNUAneg[cForNUA]){
	  iBinNUA = fHCorrectNUAneg[cForNUA]->FindBin(pVtxZ,dPhi2,dEta2);
	  w2NUA   = fHCorrectNUAneg[cForNUA]->GetBinContent(iBinNUA);  
	}
	else{ w2NUA = 1.0; }
      }

      //---------- Remove track2 from EP calculation ---------
      sumQxTPCneg2 = sumQxTPCneg; //second copy to remove 2nd track from EP
      sumQyTPCneg2 = sumQyTPCneg;
      sumQxTPCpos2 = sumQxTPCpos;
      sumQyTPCpos2 = sumQyTPCpos;
      //------------------------------------------------------


      if(dEta2 < -0.05){
	sumQxTPCneg2 -= w2NUA*TMath::Cos(gPsiN*dPhi2);
	sumQyTPCneg2 -= w2NUA*TMath::Sin(gPsiN*dPhi2); // [0] = eta <-0.05
      }
      else if(dEta2 > 0.05){
	sumQxTPCpos2 -= w2NUA*TMath::Cos(gPsiN*dPhi2);
	sumQyTPCpos2 -= w2NUA*TMath::Sin(gPsiN*dPhi2); // [1] = eta > 0.05
      }

      /*
      // track by track EP:
      PsiNTPCA = (1.0/gPsiN)*( TMath::ATan2(sumQyTPCneg2,sumQxTPCneg2) );
      if(PsiNTPCA<0.) PsiNTPCA += 2*TMath::Pi()/gPsiN;
      PsiNTPCC = (1.0/gPsiN)*( TMath::ATan2(sumQyTPCpos2,sumQxTPCpos2) );
      if(PsiNTPCC<0.) PsiNTPCC += 2*TMath::Pi()/gPsiN; */

      PsiNTPCA = (1.0/gPsiN)*( TMath::ATan2(sumQyTPCneg2,sumQxTPCneg2) + TMath::Pi() );
      //if(PsiNTPCA<0.) PsiNTPCA += 2*TMath::Pi()/gPsiN;
      PsiNTPCC = (1.0/gPsiN)*( TMath::ATan2(sumQyTPCpos2,sumQxTPCpos2) + TMath::Pi() );
      //if(PsiNTPCC<0.) PsiNTPCC += 2*TMath::Pi()/gPsiN;


      fHTPCAEventPlaneVsCent->Fill(EvtCent,PsiNTPCA);
      fHTPCCEventPlaneVsCent->Fill(EvtCent,PsiNTPCC);


      // combined weight for EP:
      WgtEP = ptw1*ptw2*w1NUA*w2NUA;
      //-----------------------------------------------------------


      //cout<<", passes 2 ";
    
      if(dChrg1!=dChrg2) {
	fHist_Corr3p_EP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEP);
	fHist_Corr3p_EP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEP);
	fHist_Corr3p_EP_Norm_PN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP);
	fHist_Corr3p_EP_Norm_PN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP);	

	//-------------> PID CME ---------------
	//Pion:
	if(isPion1 && isPion2){
	  fHist_Corr3p_Pion_EP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPPion);
	  fHist_Corr3p_Pion_EP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPPion);
	  fHist_Corr3p_Pion_EP_Norm_PN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPPion);
	  fHist_Corr3p_Pion_EP_Norm_PN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPPion);	
	}
	//Kaon:
	if(isKaon1 && isKaon2){
	  fHist_Corr3p_Kaon_EP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPKaon);
	  fHist_Corr3p_Kaon_EP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPKaon);
	  fHist_Corr3p_Kaon_EP_Norm_PN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPKaon);
	  fHist_Corr3p_Kaon_EP_Norm_PN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPKaon);	
	}
	//Proton:
	if(isProton1 && isProton2){
          //cout<<"#pair pt1 = "<<dPt1<<"\tpt2 = "<<dPt2<<"\tnSigp1 = "<<nSigTPCproton<<"\tnSigp2 = "<<nSigTPCproton2<<endl;
	  fHist_Corr3p_Proton_EP_Norm_PN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPProton);
	  fHist_Corr3p_Proton_EP_Norm_PN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPProton);
	  fHist_Corr3p_Proton_EP_Norm_PN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPProton);
	  fHist_Corr3p_Proton_EP_Norm_PN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPProton);	
	}
	//------------------------------------
      }
      else if(dChrg1>0 && dChrg2>0 && skipPairHBT==0) {
	fHist_Corr3p_EP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEP);
	fHist_Corr3p_EP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEP);
	fHist_Corr3p_EP_Norm_PP[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP);	
	fHist_Corr3p_EP_Norm_PP[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP);	

	//-------------> PID CME ---------------
	//Pion:
	if(isPion1 && isPion2){
	  fHist_Corr3p_Pion_EP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPPion);
	  fHist_Corr3p_Pion_EP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPPion);
	  fHist_Corr3p_Pion_EP_Norm_PP[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPPion);
	  fHist_Corr3p_Pion_EP_Norm_PP[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPPion);	
	}
	//Kaon:
	if(isKaon1 && isKaon2){
	  fHist_Corr3p_Kaon_EP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPKaon);
	  fHist_Corr3p_Kaon_EP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPKaon);
	  fHist_Corr3p_Kaon_EP_Norm_PP[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPKaon);
	  fHist_Corr3p_Kaon_EP_Norm_PP[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPKaon);	
	}
	//Proton:
	if(isProton1 && isProton2){
	  fHist_Corr3p_Proton_EP_Norm_PP[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPProton);
	  fHist_Corr3p_Proton_EP_Norm_PP[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPProton);
	  fHist_Corr3p_Proton_EP_Norm_PP[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPProton);
	  fHist_Corr3p_Proton_EP_Norm_PP[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPProton);	
	}
	//------------------------------------
      }

      else if(dChrg1<0 && dChrg2<0 && skipPairHBT==0){
	fHist_Corr3p_EP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEP);
	fHist_Corr3p_EP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEP);
	fHist_Corr3p_EP_Norm_NN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEP);	
	fHist_Corr3p_EP_Norm_NN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEP);	

	//-------------> PID CME ---------------
	//Pion:
	if(isPion1 && isPion2){
	  fHist_Corr3p_Pion_EP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPPion);
	  fHist_Corr3p_Pion_EP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPPion);
	  fHist_Corr3p_Pion_EP_Norm_NN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPPion);
	  fHist_Corr3p_Pion_EP_Norm_NN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPPion);	
	}
	//Kaon:
	if(isKaon1 && isKaon2){
	  fHist_Corr3p_Kaon_EP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPKaon);
	  fHist_Corr3p_Kaon_EP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPKaon);
	  fHist_Corr3p_Kaon_EP_Norm_NN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPKaon);
	  fHist_Corr3p_Kaon_EP_Norm_NN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPKaon);	
	}
	//Proton:
	if(isProton1 && isProton2){
	  fHist_Corr3p_Proton_EP_Norm_NN[QAindex][0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0A), WgtEPProton);
	  fHist_Corr3p_Proton_EP_Norm_NN[QAindex][1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNV0C), WgtEPProton);
	  fHist_Corr3p_Proton_EP_Norm_NN[QAindex][2]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCA),WgtEPProton);
	  fHist_Corr3p_Proton_EP_Norm_NN[QAindex][3]->Fill(EvtCent, TMath::Cos(n*dPhi1 + m*dPhi2 - p*PsiNTPCC),WgtEPProton);	
	}
	//------------------------------------
      } 

      //cout<<", passes 3 WgtEP = "<<WgtEP<<endl; 
    }//-------- nested track loop ends ------------------









    //============ PID business starts here =============
    
    fHistEtaPtAfter->Fill(eta,pT);
    fHistTPCdEdxvsPAfter->Fill(mom*charge,dEdx);

    if(TOFmatch>0 && probMis < 0.01){
      fHistTOFBetavsPAfter->Fill(mom*charge,beta);
    }

    //nSigTOFpion=fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    //nSigTOFkaon=fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    //nSigTOFproton=fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

    fHistTOFMatchCount->Fill(TOFmatch+2,nSigTOFpion);
    fHistTOFMatchCount->Fill(TOFmatch+4,nSigTOFkaon);
    fHistTOFMatchCount->Fill(TOFmatch+6,nSigTOFproton);

    if(!TOFmatch || probMis > 0.01){  // I dont want mismatched track in my signal distribution
      nSigTOFpion   = -9.99;
      nSigTOFkaon   = -9.99;
      nSigTOFproton = -9.99;
    }

    //nSigTPCpion   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    //nSigTPCkaon   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    //nSigTPCproton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

    //0=pi, 1=K, 2=Proton
    fHistTPCTOFnSigmavsPtAfter[0]->Fill(pT,nSigTPCpion,nSigTOFpion);
    fHistTPCTOFnSigmavsPtAfter[1]->Fill(pT,nSigTPCkaon,nSigTOFkaon);
    fHistTPCTOFnSigmavsPtAfter[2]->Fill(pT,nSigTPCproton,nSigTOFproton);





    
    if(TMath::Abs(nSigTPCpion)<=fNSigmaCut){
      fHistTOFnSigmavsPtAfter[0]->Fill(pT*charge,nSigTOFpion);
      fHistPtwithTPCNsigma[0]->Fill(pT*charge);
      fHistTPCdEdxvsPtPIDAfter[0]->Fill(pT,dEdx);
      if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        fHistPtwithTOFSignal[0]->Fill(pT*charge);
      }
    }
    if(TMath::Abs(nSigTPCkaon)<=fNSigmaCut){
      fHistPtwithTPCNsigma[1]->Fill(pT*charge);
      fHistTOFnSigmavsPtAfter[1]->Fill(pT*charge,nSigTOFkaon);
      fHistTPCdEdxvsPtPIDAfter[1]->Fill(pT,dEdx);
      if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        fHistPtwithTOFSignal[1]->Fill(pT*charge);
      }
    }
    if(TMath::Abs(nSigTPCproton)<=fNSigmaCut){
      fHistPtwithTPCNsigma[2]->Fill(pT*charge);
      fHistTOFnSigmavsPtAfter[2]->Fill(pT*charge,nSigTOFproton);
      fHistTPCdEdxvsPtPIDAfter[2]->Fill(pT,dEdx);
      if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        fHistPtwithTOFSignal[2]->Fill(pT*charge);
      }
    }
   



    //========> nSigmaTOF distribution for circular cut <===============
    //if(TMath::Sqrt(nSigTPCpion*nSigTPCpion+nSigTOFpion*nSigTOFpion)<=fNSigmaCut){
      //fHistTOFnSigmavsPtAfter[0]->Fill(pT*charge,nSigTOFpion);
      //if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        //fHistPtwithTOFSignal[0]->Fill(pT*charge);
      //}
    //}
    //if(TMath::Sqrt(nSigTPCkaon*nSigTPCkaon+nSigTOFkaon*nSigTOFkaon)<=fNSigmaCut){
      //fHistTOFnSigmavsPtAfter[1]->Fill(pT*charge,nSigTOFkaon);
      //if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        //fHistPtwithTOFSignal[1]->Fill(pT*charge);
      //}
    //}
    //if(TMath::Sqrt(nSigTPCproton*nSigTPCproton+nSigTOFproton*nSigTOFproton)<=fNSigmaCut){
      //fHistTOFnSigmavsPtAfter[2]->Fill(pT*charge,nSigTOFproton);
      //if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        //fHistPtwithTOFSignal[2]->Fill(pT*charge);
      //}
    //} 


    

    if(TMath::Abs(nSigTPCpion)<=fNSigmaCut){
      fHistPtwithTPCNsigma[0]->Fill(pT*charge);
      fHistTPCdEdxvsPtPIDAfter[0]->Fill(pT,dEdx);
    }
    if(TMath::Abs(nSigTPCkaon)<=fNSigmaCut){
      fHistPtwithTPCNsigma[1]->Fill(pT*charge);
      fHistTPCdEdxvsPtPIDAfter[1]->Fill(pT,dEdx);
    }
    if(TMath::Abs(nSigTPCproton)<=fNSigmaCut){
      fHistPtwithTPCNsigma[2]->Fill(pT*charge);
      fHistTPCdEdxvsPtPIDAfter[2]->Fill(pT,dEdx);
    }


 



    //-------------- Fill NUA for Charged tracks ----------------
    if(charge>0){
      fHist3DEtaPhiVz_Pos_Run[3][cForNUA]->Fill(pVtxZ,phi,eta);
    }
    else if(charge<0){
      fHist3DEtaPhiVz_Neg_Run[3][cForNUA]->Fill(pVtxZ,phi,eta);
    }
     
    if(bFillNUAHistPID){
      //============== Fill NUA Histograms for Pion ---------------------
      if(pT<0.6 && TMath::Abs(nSigTPCpion)<=2.5){
	if(charge>0){
	  fHist3DEtaPhiVz_Pos_Run[0][cForNUA]->Fill(pVtxZ,phi,eta);
	}
	else if(charge<0){
	  fHist3DEtaPhiVz_Neg_Run[0][cForNUA]->Fill(pVtxZ,phi,eta);
	}
      }
      else if(pT>=0.6 && pT<=2.0 && TMath::Abs(nSigTPCpion)<=2.5 && TMath::Abs(nSigTOFpion)<=2.0 ){
	if(charge>0){
	  fHist3DEtaPhiVz_Pos_Run[0][cForNUA]->Fill(pVtxZ,phi,eta);
	}
	else if(charge<0){
	  fHist3DEtaPhiVz_Neg_Run[0][cForNUA]->Fill(pVtxZ,phi,eta);
	}
      }
     
      //============== Fill NUA Histograms for Kaon ---------------------
      if(pT<0.6 && TMath::Abs(nSigTPCkaon)<=2.5){
	if(charge>0){
	  fHist3DEtaPhiVz_Pos_Run[1][cForNUA]->Fill(pVtxZ,phi,eta);
	}
	else if(charge<0){
	  fHist3DEtaPhiVz_Neg_Run[1][cForNUA]->Fill(pVtxZ,phi,eta);
	}
      }
      else if(pT>=0.6 && pT<=2.0 && TMath::Abs(nSigTPCkaon)<=2.5 && TMath::Abs(nSigTOFkaon)<=2.0){
	if(charge>0){
	  fHist3DEtaPhiVz_Pos_Run[1][cForNUA]->Fill(pVtxZ,phi,eta);
	}
	else if(charge<0){
	  fHist3DEtaPhiVz_Neg_Run[1][cForNUA]->Fill(pVtxZ,phi,eta);
	}
      }

      //============== Fill NUA Histograms for proton ---------------------
      if(pT<0.8 && TMath::Abs(nSigTPCproton)<=2.5){
	if(charge>0){
	  fHist3DEtaPhiVz_Pos_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
	}
	else if(charge<0){
	  fHist3DEtaPhiVz_Neg_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
	}
      }
      else if(pT>=0.8 && pT<=3.5 && TMath::Abs(nSigTPCproton)<=2.5 && TMath::Abs(nSigTOFproton)<=2.5){
	if(charge>0){
	  fHist3DEtaPhiVz_Pos_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
	}
	else if(charge<0){
	  fHist3DEtaPhiVz_Neg_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
	}
      }

    }// Fill NUA for PID or not?


    //------------------------------------------------------------
    if(!TOFmatch || probMis > 0.01 || beta>0.2) continue;

    // nSigmaTPC distribution for Fixed nSigmaTOF cut
   
    if(TMath::Abs(nSigTOFpion)<=fNSigmaCut){
      fHistTPCnSigmavsPtAfter[0]->Fill(pT*charge,nSigTPCpion);
      fHistPtwithTOFmasscut[0]->Fill(pT*charge);
    }
    if(TMath::Abs(nSigTOFkaon)<=fNSigmaCut){
      fHistTPCnSigmavsPtAfter[1]->Fill(pT*charge,nSigTPCkaon);
      fHistPtwithTOFmasscut[1]->Fill(pT*charge);
    }
    if(TMath::Abs(nSigTOFproton)<=fNSigmaCut){
      fHistTPCnSigmavsPtAfter[2]->Fill(pT*charge,nSigTPCproton);
      fHistPtwithTOFmasscut[2]->Fill(pT*charge);
    }





  }
 //===================== track loop ends ============================








  //---------- Do event by event business here  -------





  // ----------------- Store Q-Vectors --------------
  //V0A-V0C 
  fHist_Reso2n_EP_Norm_Det[QAindex][0]->Fill(EvtCent, TMath::Cos(gPsiN*(PsiNV0A-PsiNV0C)));
  //V0A-TPC 
  fHist_Reso2n_EP_Norm_Det[QAindex][1]->Fill(EvtCent, TMath::Cos(gPsiN*(PsiNV0A-PsiNTPCF)));
  //V0C-TPC 
  fHist_Reso2n_EP_Norm_Det[QAindex][2]->Fill(EvtCent, TMath::Cos(gPsiN*(PsiNV0C-PsiNTPCF)));
  //TPCa -TPCc 
  fHist_Reso2n_EP_Norm_Det[QAindex][3]->Fill(EvtCent, TMath::Cos(gPsiN*(PsiNTPCA-PsiNTPCC)));




  //---------  Store TPC-Qn for Recenter ---------
  sumTPCQn2x[0] = sumTPCQn2x[0]/multEtaNeg;
  sumTPCQn2y[0] = sumTPCQn2y[0]/multEtaNeg;
  sumTPCQn3x[0] = sumTPCQn3x[0]/multEtaNeg;
  sumTPCQn3y[0] = sumTPCQn3y[0]/multEtaNeg;
  sumTPCQn4x[0] = sumTPCQn4x[0]/multEtaNeg;
  sumTPCQn4y[0] = sumTPCQn4y[0]/multEtaNeg;

  sumTPCQn2x[1] = sumTPCQn2x[1]/multEtaPos;
  sumTPCQn2y[1] = sumTPCQn2y[1]/multEtaPos;
  sumTPCQn3x[1] = sumTPCQn3x[1]/multEtaPos;
  sumTPCQn3y[1] = sumTPCQn3y[1]/multEtaPos;
  sumTPCQn4x[1] = sumTPCQn4x[1]/multEtaPos;
  sumTPCQn4y[1] = sumTPCQn4y[1]/multEtaPos;

  fTPCAQ2xVsCentRun->Fill(EvtCent,sumTPCQn2x[0]); 
  fTPCAQ2yVsCentRun->Fill(EvtCent,sumTPCQn2y[0]);
  fTPCCQ2xVsCentRun->Fill(EvtCent,sumTPCQn2x[1]);
  fTPCCQ2yVsCentRun->Fill(EvtCent,sumTPCQn2y[1]);

  fTPCAQ3xVsCentRun->Fill(EvtCent,sumTPCQn3x[0]); 
  fTPCAQ3yVsCentRun->Fill(EvtCent,sumTPCQn3y[0]);
  fTPCCQ3xVsCentRun->Fill(EvtCent,sumTPCQn3x[1]);
  fTPCCQ3yVsCentRun->Fill(EvtCent,sumTPCQn3y[1]);

  fTPCAQ4xVsCentRun->Fill(EvtCent,sumTPCQn4x[0]); 
  fTPCAQ4yVsCentRun->Fill(EvtCent,sumTPCQn4y[0]);
  fTPCCQ4xVsCentRun->Fill(EvtCent,sumTPCQn4x[1]);
  fTPCCQ4yVsCentRun->Fill(EvtCent,sumTPCQn4y[1]);


    



  PostData(1,fListHist);

  fHistEventCount->Fill(14.5); //15th bin is last one
  stepCount++;




}//================ UserExec ==============

























//////////////////// FUNCTIONS //////////////////////




double AliAnalysisTaskCMEV0PID::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    AliDebug(2,"\n\n ::GetWDist => One of vertices is not valid\n\n");
    return 0;
  }
  static TMatrixDSym vVb(3);
  double dist = -1;
  double dx = v0->GetX()-v1->GetX();
  double dy = v0->GetY()-v1->GetY();
  double dz = v0->GetZ()-v1->GetZ();
  double cov0[6],cov1[6];
  v0->GetCovarianceMatrix(cov0);
  v1->GetCovarianceMatrix(cov1);
  vVb(0,0) = cov0[0]+cov1[0];
  vVb(1,1) = cov0[2]+cov1[2];
  vVb(2,2) = cov0[5]+cov1[5];
  vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
  vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
  vVb.InvertFast();
  if (!vVb.IsValid()) {
    AliDebug(2,"Singular Matrix\n");
    return dist;
  }
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
  +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;
}


Bool_t AliAnalysisTaskCMEV0PID::PileUpMultiVertex(const AliAODEvent* faod)
 {  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if(!(nPlp=faod->GetNumberOfPileupVerticesTracks()))
  return kFALSE;

  vtPrm = faod->GetPrimaryVertex();
  if(vtPrm == faod->GetPrimaryVertexSPD())
  return kTRUE;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for(int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)faod->GetPileupVertexTracks(ipl);
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF()    > kMaxPlpChi2)    continue;
    //int bcPlp = vtPlp->GetBC();
    //if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2)
    // return kTRUE; // pile-up from other BC

    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist)        continue;

    return kTRUE; // pile-up: well separated vertices
    }
  return kFALSE;
}



void AliAnalysisTaskCMEV0PID::SetupMCcorrectionMap(TString sMCfilePath){

  if(bApplyMCcorr){
    if(!gGrid){
      TGrid::Connect("alien://");
    }

    if(!mfileFBHijing){

      mfileFBHijing = TFile::Open(sMCfilePath,"READ");

      fListFBHijing = dynamic_cast<TList*>(mfileFBHijing->FindObjectAny("fMcEffiHij"));

      if(!fListFBHijing){
	std::cout<<"\n\n !!!!**** Warning: FB Efficiency File/List not found *****\n\n"<<std::endl;
	//exit(1);
      }
      else if(fListFBHijing) {
        cout<<"\n =========> Info: FB Efficiency is Used from file = "<<sMCfilePath.Data()<<endl;
	for(int i=0;i<10;i++) {
	  fFB_Efficiency_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_%d",i));
	  //std::cout<<" input MC hist"<<i<<" = "<<fFB_Efficiency_Cent[i]->GetName()<<std::endl;
	}
      }
    }
  }
  else{ // if MC efficiency Not used/ file not found, then use weight = 1.
    for(int i=0;i<10;i++){
      fFB_Efficiency_Cent[i] = new TH1D(Form("eff_unbiased_%d",i),"",1,0,50.); 
      fFB_Efficiency_Cent[i]->SetBinContent(1,1.0);
    }
    if(bApplyMCcorr){ printf("\n!!*****  !!!! WARNING !!!! *****!!\n MC correction File not found, using Efficiency = 1.0 !!\n\n");}
  }   
}




Int_t AliAnalysisTaskCMEV0PID::GetCentralityScaled0to10(Float_t fCent){

 Int_t cIndex = 0;

 if(fCent<5.0) {
   cIndex  = 0; 
 }
 else if(fCent>=5.0 && fCent<10){
   cIndex  = 1;
 }
 else if(fCent>=10.0) {
   cIndex = abs(fCent/10.0)+1;
 }
 return cIndex;
}

void AliAnalysisTaskCMEV0PID::SetUpCentralityOutlierCut(){
  //std::cout<<" centrality outlier function called "<<std::endl;
 Float_t fMeanTPC[100] = {2902.95,2758.33,2642.78,2536.67,2435.37,2340.06,2248.44,2163.71,2080.49,2001.54,1925.86,1852.64,1781.97,1715.56,1650.53,1587.23,1527.51,1468.19,1412.73,1357.86,1305.35,1254.33,1205.57,1157.28,1111.53,1066.42,1023.15,981.594,940.795,901.766,863.651,826.183,790.53,756.358,722.654,690.513,659.443,628.807,599.748,571.664,544.446,518.042,492.369,468.072,444.694,422.487,400.104,379.129,359.147,339.62,320.817,302.788,285.791,269.015,253.688,238.671,224.039,209.932,196.915,184.647,172.76,161.381,150.395,140.288,131.033,121.58,113.112,104.938,97.3078,90.2178,83.5974,77.2645,70.7126,65.4424,60.1404,55.5644,50.8314,46.3761,43.024,38.625,35.3435,32.2304,29.4192,26.821,24.3303,21.9332,19.4215,16.7163,14.9414,13.1092,0.};

 Float_t fSigmaTPC[100] = {122.209,107.901,103.452,100.498,97.7403,94.7845,93.2543,90.0548,88.1106,85.7382,84.0812,82.2978,80.3817,78.6002,77.3448,75.5086,73.6842,71.9733,70.3447,69.1999,67.878,66.3511,65.0406,63.4866,62.4409,60.7899,59.1328,58.426,56.8618,55.8871,54.1031,53.4959,52.0482,51.0441,49.6218,48.7646,47.5166,46.5247,45.0727,44.4311,43.4531,42.0404,41.0238,40.1384,39.2588,38.2461,36.5951,36.0552,35.3727,33.7883,32.7167,32.4486,31.3709,30.3444,29.505,28.5139,27.4471,26.5359,25.9506,25.127,24.3797,23.2985,22.279,21.4698,20.781,20.8193,19.9509,18.8036,17.9145,16.961,16.7375,15.852,14.9324,14.7663,13.5969,13.4533,12.3067,12.7835,11.7283,10.6758,10.6676,10.6492,9.04614,8.89065,8.66093,8.50997,7.98812,6.91087,7.12045,7.29593,0.};

 for(int i=0;i<90;i++) {
   hCentvsTPCmultCuts->SetBinContent(i+1,1,fMeanTPC[i]);
   hCentvsTPCmultCuts->SetBinContent(i+1,2,fSigmaTPC[i]);
 }
}






void AliAnalysisTaskCMEV0PID::SetupEventAndTaskConfigInfo(){

  fHistTaskConfigParameters = new TH1F("fHistTaskConfigParameters","Task parameters",20,0,20);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(1,"FilterBit");
  fHistTaskConfigParameters->SetBinContent(1,fFilterBit);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(2,"n#sigmaTPC");
  fHistTaskConfigParameters->SetBinContent(2,fNSigmaCut);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(3,"MinPt");
  fHistTaskConfigParameters->SetBinContent(3,fMinPtCut);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(4,"MaxPt");
  fHistTaskConfigParameters->SetBinContent(4,fMaxPtCut);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(5,"MinEta");
  fHistTaskConfigParameters->SetBinContent(5,fMinEtaCut);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(6,"MaxEta");
  fHistTaskConfigParameters->SetBinContent(6,fMaxEtaCut);

  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(11,"CentralityMin");
  fHistTaskConfigParameters->SetBinContent(11,fCentralityPercentMin);
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(12,"CentralityMax");
  fHistTaskConfigParameters->SetBinContent(12,fCentralityPercentMax);

  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(13,"VertexMin(cm)");
  fHistTaskConfigParameters->GetXaxis()->SetBinLabel(14,"VertexMax(cm)");

  fListHist->Add(fHistTaskConfigParameters);



  fHistPileUpCount = new TH1F("fHistPileUpCount", "fHistPileUpCount", 15, 0., 15.);
  fHistPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fHistPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fHistPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultComb08");
  fHistPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fHistPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>5.0");
  fHistPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fHistPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  Int_t puConst = fPileUpConstParm;
  fHistPileUpCount->GetXaxis()->SetBinLabel(8,Form("multESDTPCDif>%d",puConst));
  fHistPileUpCount->GetXaxis()->SetBinLabel(9,Form("multGlobTPCDif>%d",puConst));
  fHistPileUpCount->GetXaxis()->SetBinLabel(10,"PileUpMultSelTask");
  fListHist->Add(fHistPileUpCount);


  fHistMultSelPUCount = new TH1F("fHistMultSelPileUpCount", "no of PU Event from MultSelTask", 10, 0., 10);
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(1,"PileUp");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(2,"PileUpMV");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(3,"PileUpMultBins");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(4,"InconsistentVtx");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(5,"TrackletVsCluster");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(6,"IncompleteDAQ");
  fHistMultSelPUCount->GetXaxis()->SetBinLabel(7,"NotGoodVertex2016");
  fListHist->Add(fHistMultSelPUCount);


  fHistEventCount = new TH1F("fHistEventCount","Event counts",15,0,15);
  fHistEventCount->GetXaxis()->SetBinLabel(1,"Called UserExec()");
  fHistEventCount->GetXaxis()->SetBinLabel(2,"Called Exec()");
  fHistEventCount->GetXaxis()->SetBinLabel(3,"AOD Exist");
  fHistEventCount->GetXaxis()->SetBinLabel(4,"Vz < 10");
  fHistEventCount->GetXaxis()->SetBinLabel(5,Form("%2.0f<Cent<%2.0f",fCentralityPercentMin,fCentralityPercentMax));
  fHistEventCount->GetXaxis()->SetBinLabel(6,"noAODtrack > 2 ");
  fHistEventCount->GetXaxis()->SetBinLabel(7,"TPC vs Global");
  fHistEventCount->GetXaxis()->SetBinLabel(8,"TPC128 vs ESD");
  fHistEventCount->GetXaxis()->SetBinLabel(9,"Cent vs TPC");
  fHistEventCount->GetXaxis()->SetBinLabel(10,"mult eta+/- > 2");
  fHistEventCount->GetXaxis()->SetBinLabel(11,"centCL1 < 90");
  fHistEventCount->GetXaxis()->SetBinLabel(15,"Survived Events");
  fListHist->Add(fHistEventCount);

  //fHistEventCount->Fill(1);

}





void AliAnalysisTaskCMEV0PID::GetV0MCorrectionHist(Int_t run)
{
  if(fListV0MCorr){
    fHCorrectV0M    = (TH1D *) fListV0MCorr->FindObject(Form("fHistV0Gain_Run%d",run));
    fHAvgerageQnV0A = (TH2D *) fListV0MCorr->FindObject(Form("fHistAvgQnV0A_Run%d",run));
    fHAvgerageQnV0C = (TH2D *) fListV0MCorr->FindObject(Form("fHistAvgQnV0C_Run%d",run));

    if(fHCorrectV0M){
      cout<<"\n =========== Info:: Setting up V0 gain correction for run = "<<run<<"============"<<endl;
    }
    else{
      cout<<"\n =========== Info:: No V0 gain correction..!!! for run = "<<run<<"============"<<endl;
    }
  }
  else{
    cout<<"\n ======== Error:: List of V0 gain correction not found for run "<<run<<"============"<<endl;
    fHCorrectV0M  = new TH1D("fHCorrectV0M","",64,0,64);
    for(int i=1;i<=64;i++){
      fHCorrectV0M->SetBinContent(i,1.0);
    }
    fHAvgerageQnV0A = new TH2D("fHAvgerageQnV0A_empty","<Cos2>,<Sin2>,<Cos3>,<Sin3> V0A",90,0,90,8,0,8);
    fHAvgerageQnV0C = new TH2D("fHAvgerageQnV0C_empty","<Cos2>,<Sin2>,<Cos3>,<Sin3> V0A",90,0,90,8,0,8);
    for(int i=1;i<=90;i++){
      for(int j=1;j<=8;j++){
        fHAvgerageQnV0A->SetBinContent(i,j,0.0);
        fHAvgerageQnV0C->SetBinContent(i,j,0.0);
      }
    }
  }
}







void AliAnalysisTaskCMEV0PID::GetNUACorrectionHist(Int_t run)
{

  if(fListNUACorr){
    for(int i=0;i<5;i++){
      fHCorrectNUApos[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pos_Cent%d_Run%d",i,run)); 
      fHCorrectNUAneg[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Neg_Cent%d_Run%d",i,run));
    }
    if(fHCorrectNUApos[0] && fHCorrectNUApos[0]){
      cout<<"\n=========== Info:: Setting up NUA corrections for run "<<run<<"============"<<endl;
    }
  }
  else {
    printf("\n ******** Warning: No NUA Correction for Charge Particle in run %d, Use Wgt = 1.0 ********* \n",run);
    for(int i=0;i<5;i++){
      fHCorrectNUApos[i] = new TH3D(Form("fHCorrectNUApos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAneg[i] = new TH3D(Form("fHCorrectNUAneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUApos[i]->SetBinContent(1,1,1,1.0);
      fHCorrectNUAneg[i]->SetBinContent(1,1,1,1.0);
      //exit(1);
    }
  }

  //===================   PID:  ==========================
  if(fListNUACorr){
    for(int i=0;i<5;i++){
      fHCorrectNUAposPion[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Pos_Cent%d_Run%d",i,run)); //
      fHCorrectNUAnegPion[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Neg_Cent%d_Run%d",i,run));

      fHCorrectNUAposKaon[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Pos_Cent%d_Run%d",i,run)); //
      fHCorrectNUAnegKaon[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Neg_Cent%d_Run%d",i,run));

      fHCorrectNUAposProton[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Pos_Cent%d_Run%d",i,run)); 
      fHCorrectNUAnegProton[i] = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Neg_Cent%d_Run%d",i,run));
    }
    if(fHCorrectNUAposPion[0] && fHCorrectNUAposKaon[0] && fHCorrectNUAposProton[0]) {
      cout<<"\n=========== Info:: Setting up --> PID NUA corrections for run =  "<<run<<"============"<<endl;
    }
    else{
      cout<<"\n=========== WARNING :: PID NUA corrections NOT found for run =  "<<run<<"============"<<endl;
    }
  }
  else {
    printf("\n ******** Warning: No NUA Correction for PID for run %d, Use Wgt = 1.0 ********* \n",run);
    for(int i=0;i<5;i++){
      fHCorrectNUAposPion[i] = new TH3D(Form("fHCorrectNUAPionpos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposPion[i] = new TH3D(Form("fHCorrectNUAPionneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposPion[i]->SetBinContent(1,1,1,1.0);
      fHCorrectNUAposPion[i]->SetBinContent(1,1,1,1.0);

      fHCorrectNUAposKaon[i] = new TH3D(Form("fHCorrectNUAKaonpos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposKaon[i] = new TH3D(Form("fHCorrectNUAKaonneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposKaon[i]->SetBinContent(1,1,1,1.0);
      fHCorrectNUAposKaon[i]->SetBinContent(1,1,1,1.0);

      fHCorrectNUAposProton[i] = new TH3D(Form("fHCorrectNUAProtonpos_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposProton[i] = new TH3D(Form("fHCorrectNUAProtonneg_cent%d",i),"",1,-10,10,1,0,6.284,1,-0.9,0.9); 
      fHCorrectNUAposProton[i]->SetBinContent(1,1,1,1.0);
      fHCorrectNUAposProton[i]->SetBinContent(1,1,1,1.0);
      //exit(1);
    }
  }

}

