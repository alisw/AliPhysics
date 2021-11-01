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
/* $Id: AliAnalysisTaskCMW.cxx  Rihan Haque, 18/09/2019 (ver1) $ */
/* This class is derived from AliAnalysisTaskCVE           */
/* This analysis Task contains trivial corrections to CMW */

//-- general include---
#include "TChain.h"
#include "TTree.h"
#include "TGrid.h"
#include "TROOT.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TMatrixDSym.h"
#include "TParticle.h"
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
#include "AliTimeRangeCut.h"

//----- For MC event------
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

//----for PID-----
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"

//----- Vevent and tracks
#include "AliVEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliCentrality.h"


//----- must include-------
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliPhysicsSelection.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskCMWPU2018dca.h"

using std::cout;
using std::endl;
using std::vector;


ClassImp(AliAnalysisTaskCMWPU2018dca)

AliAnalysisTaskCMWPU2018dca::AliAnalysisTaskCMWPU2018dca(const char *name): AliAnalysisTaskSE(name),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fListTRKCorr(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),  
  fCentralityMin(0),
  fCentralityMax(90),
  gHarmonic(2),
  fParticle(1),
  fFilterBit(1),
  fTPCclustMin(70),
  bUseKinkTracks(kFALSE),
  fNSigmaTPCCut(2.0),
  fNSigmaTOFCut(2.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fDCAxyMax(2.4),
  fDCAzMax(3.2),
  fChi2(4.0),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  bSkipPileUpCut(kFALSE), 
  fEtaGapNeg(-0.1),
  fEtaGapPos(0.1),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),    
  fdEdxMin(10.0),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  sCentrEstimator("V0M"),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fHistAChrgVsCent(NULL),
  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsV0MBefore(NULL),
  fHistCL0VsV0MBefore(NULL),  
  fHistTPConlyVsCL1After(NULL),
  fHistTPConlyVsV0MAfter(NULL),
  fHistCL0VsV0MAfter(NULL),
  fHistGlobalVsV0MBefore(NULL),
  fHistGlobalVsV0MAfter(NULL),
  fTPCvsGlobalTrkBefore(NULL),
  fTPCvsGlobalTrkAfter(NULL),
  fHistTPCVsESDTrkBefore(NULL),
  fHistTPCVsESDTrkAfter(NULL),
  fHdcaxy(NULL),
  fHdcaz(NULL),

  fHistPileUpCount(NULL),
  fHistEventCount(NULL),
  fHCorrectEVNTWGTChrg(NULL),
  fHFillNUAPosPID(NULL),
  fHFillNUANegPID(NULL),
  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL)   
{

  //for(int i=0; i<5; i++){
  //fHFillNUAPosPID[i]  = NULL;
  //fHFillNUANegPID[i]  = NULL; 
  //}

  
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskCMWPU2018dca::AliAnalysisTaskCMWPU2018dca():
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fListTRKCorr(NULL),
  fListNUACorr(NULL),
  fListV0MCorr(NULL),
  fCentralityMin(0),
  fCentralityMax(90),
  gHarmonic(2),
  fParticle(1),
  fFilterBit(1),
  fTPCclustMin(70),
  bUseKinkTracks(kFALSE),
  fNSigmaTPCCut(2.0),
  fNSigmaTOFCut(2.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fDCAxyMax(2.4),
  fDCAzMax(3.2),
  fChi2(4.0),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  bSkipPileUpCut(kFALSE), 
  fEtaGapNeg(-0.1),
  fEtaGapPos(0.1),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),    
  fdEdxMin(10.0),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  sCentrEstimator("V0M"),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),  
  fHistAChrgVsCent(NULL),

  fHistTPConlyVsCL1Before(NULL),
  fHistTPConlyVsV0MBefore(NULL),
  fHistCL0VsV0MBefore(NULL),
  fHistTPConlyVsCL1After(NULL),
  fHistTPConlyVsV0MAfter(NULL),
  fHistCL0VsV0MAfter(NULL),
  fHistGlobalVsV0MBefore(NULL),
  fHistGlobalVsV0MAfter(NULL),
  fTPCvsGlobalTrkBefore(NULL),
  fTPCvsGlobalTrkAfter(NULL),
  fHistTPCVsESDTrkBefore(NULL),
  fHistTPCVsESDTrkAfter(NULL),
  fHdcaxy(NULL),
  fHdcaz(NULL),

  fHistPileUpCount(NULL),
  fHistEventCount(NULL),
  fHCorrectEVNTWGTChrg(NULL),
  fHFillNUAPosPID(NULL),
  fHFillNUANegPID(NULL),
  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL)    
{

  //for(int i=0; i<5; i++){
  //fHFillNUAPosPID[i]  = NULL;
  //fHFillNUANegPID[i]  = NULL; 
  //}
  
  //Not needed for Empty Constructor:
  //DefineInput(0,TChain::Class());
  //DefineOutput(1,TList::Class());
}
  
//__________________ destructor ___________________
AliAnalysisTaskCMWPU2018dca::~AliAnalysisTaskCMWPU2018dca()
{
  if(fListHist)      delete fListHist;  
  if(fAnalysisUtil)  delete fAnalysisUtil;   // because its 'new' !!

  //Delete the clones
  if(fListTRKCorr)  delete fListTRKCorr;
  if(fListNUACorr)  delete fListNUACorr;
  if(fListV0MCorr)  delete fListV0MCorr;

  //delete the functions:
  if(fSPDCutPU)     delete fSPDCutPU;
  if(fV0CutPU)      delete fV0CutPU;
  if(fMultCutPU)    delete fMultCutPU;
  if(fCenCutLowPU)  delete fCenCutLowPU; 
  if(fCenCutHighPU) delete fCenCutHighPU; 
      
}










//________________ Define Histograms _______________
void AliAnalysisTaskCMWPU2018dca::UserCreateOutputObjects()
{
  //std::cout<<"\n UserCreateOutputObject: function begins...\n"<<endl; 
  //Get The Input Hander:
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  
  if (!inputHandler) {  printf("\n***** ERROR *****\n Input handler missing, Status:QUIT!\n");    exit(1);}




  //PileUp Multi-Vertex
  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);
  
  //// obtain the PID response object if needed:
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) {  printf("\n***** ERROR *****\n fPIDResponse missing, Status:QUIT!\n");    exit(1);}    
    
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);

  SetupEventAndTaskConfigInfo();


  /// Bunch of QA histograms:  
  fCentDistBeforCut = new TH1F("fCentDistBeforCut","Cent w/o any Cuts; Cent (%); no.Events ",100,0,100);
  fListHist->Add(fCentDistBeforCut);
  fCentDistAfterCut = new TH1F("fCentDistAfterCut","Cent with all Cut; Cent (%); no.Events ",100,0,100);
  fListHist->Add(fCentDistAfterCut);

  fHdcaxy = new TH1F("fHdcaxy","fHdcaxy",1000,-3,3);
  fListHist->Add(fHdcaxy);
  fHdcaz = new TH1F("fHdcaz","fHdcaz",1000,-3,3);
  fListHist->Add(fHdcaz);

  /*  
  fHistEtaPtBeforCut = new TH2F("fHistEtaPtBeforCut","#eta vs p_{T} (wFB, w/o  cut)",24,-1.2,1.2,50,0,5);
  fListHist->Add(fHistEtaPtBeforCut);
  fHistEtaPhiBeforCut = new TH2F("fHistPhiEtaBeforCut","#phi vs #eta (wFB, w/o cut)",50,0,6.2835,24,-1.2,1.2);
  fListHist->Add(fHistEtaPhiBeforCut);  
  fHistEtaPhiAfterCut = new TH2F("fHistPhiEtaAfterCut","#phi vs #eta (with Wgts)",50,0,6.2835,24,-1.2,1.2);
  fListHist->Add(fHistEtaPhiAfterCut);
  */



  Int_t gMaxGlobalmult  = 4000;
  Int_t gMaxTPCcorrmult = 5000;
  Int_t gMaxESDtracks   = 20000;





  fHistGlobalVsV0MBefore = new TH2F("fHistGlobalVsV0MBefore","Before;Cent(V0M);Global",100,0,100,250,0,gMaxGlobalmult);
  fListHist->Add(fHistGlobalVsV0MBefore);
  fHistGlobalVsV0MAfter  = new TH2F("fHistGlobalVsV0MAfter"," After; Cent(V0M);Global",100,0,100,250,0,gMaxGlobalmult);
  fListHist->Add(fHistGlobalVsV0MAfter);

  fHistTPConlyVsCL1Before = new TH2F("fHistTPConlyVsCL1Before","Before;Cent(CL1); TPC(FB128)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1Before);
  fHistTPConlyVsCL1After  = new TH2F("fHistTPConlyVsCL1After","After; Cent(CL1); TPC(FB128) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsCL1After);

  fHistTPConlyVsV0MBefore = new TH2F("fHistTPConlyVsV0MBefore","Before;Cent(V0M); TPC(FB128)",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MBefore);
  fHistTPConlyVsV0MAfter  = new TH2F("fHistTPConlyVsV0MAfter","After; Cent(V0M); TPC(FB128) ",100,0,100,250,0,gMaxTPCcorrmult);
  fListHist->Add(fHistTPConlyVsV0MAfter);


  fHistCL0VsV0MBefore = new TH2F("fHistCL0VsV0MBefore","Before;Cent(V0M); Cent(CL0)",100,0,100,100,0,100);
  fListHist->Add(fHistCL0VsV0MBefore);
  fHistCL0VsV0MAfter  = new TH2F("fHistCL0VsV0MAfter","After; Cent(V0M); Cent(CL0) ",100,0,100,100,0,100);
  fListHist->Add(fHistCL0VsV0MAfter);


  fTPCvsGlobalTrkBefore = new TH2F("fTPCvsGlobalTrkBefore","Global(FB32) vs TPC(FB128)",250,0,5000,250,0,5000);
  fListHist->Add(fTPCvsGlobalTrkBefore);
  fTPCvsGlobalTrkAfter = new TH2F("fTPCvsGlobalTrkAfter","Global(FB32) vs TPC(FB128)",250,0,5000,250,0,5000);
  fListHist->Add(fTPCvsGlobalTrkAfter);
  
  //fHistTPCVsESDTrkBefore = new TH2F("fHistTPCVsESDTrkBefore","Before; TPC1; ESD trk",250,0,gMaxTPCcorrmult,250,0,gMaxESDtracks);
  //fListHist->Add(fHistTPCVsESDTrkBefore);
  //fHistTPCVsESDTrkAfter  = new TH2F("fHistTPCVsESDTrkAfter"," After;  TPC1; ESD trk",250,0,gMaxTPCcorrmult,250,0,gMaxESDtracks);
  //fListHist->Add(fHistTPCVsESDTrkAfter);






  
  //----------- User's histograms: --------------
  Char_t  name[1000];
  Char_t title[1000];
 

  Double_t centRange[11] = {0,5,10,20,30,40,50,60,70,80,90};

  fHistAChrgVsCent = new TH2F("fHistAChrgVsCent","Ach vs Cent;Cent;Ach",18,0,90,1000,-1.0,1.0);
  fListHist->Add(fHistAChrgVsCent);



  Int_t gCentForNUA[6] = {0,5,10,20,40,90};
  Char_t cpid[10];

  if(fParticle==1)      sprintf(cpid,"Pion,Id %d",fParticle);
  else if(fParticle==2) sprintf(cpid,"Kaon,Id %d",fParticle);
  else if(fParticle==3) sprintf(cpid,"Prot,Id %d",fParticle);
  else  sprintf(cpid,"Charge,Id %d",fParticle);


  sprintf(name,"fHistEtaPhiVz_%d_Pos_Cent%d_Run%d",fParticle,0,1); 
  sprintf(title,"%s Pos, Cent%d-%d, FB %d",cpid,0,90,fFilterBit);
  fHFillNUAPosPID = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
  fListHist->Add(fHFillNUAPosPID);

  sprintf(name,"fHistEtaPhiVz_%d_Neg_Cent%d_Run%d",fParticle,0,1); 
  sprintf(title,"%s Neg, Cent%d-%d, FB %d",cpid,0,90,fFilterBit);
  fHFillNUANegPID = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
  fListHist->Add(fHFillNUANegPID);

  
  //// PileUp Removal Functions:
  fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);
  
  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);



  
    


  if(fListTRKCorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist for MC tracking Efficiency Found.!!\n"<<std::endl;
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList for Trk Efficiency Correction!!\n using TrkWgt = 1.0 \n "<<std::endl;
  }

  if(fListNUACorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist for NUA Correction Found.!!\n"<<std::endl;
    //fListNUACorr->ls(); 
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList NUA Correction!!\n using NUAWgt = 1.0 \n "<<std::endl;
  }

  if(fListV0MCorr){
    std::cout<<"\n UserCreateOutputObject::Info() Tlist for EVNTWGT Correction Found.!!\n"<<std::endl;
    //fListV0MCorr->ls();
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList EVNTWGT Correction!!\n using EVNTWGTWgt = 1.0 \n "<<std::endl;
  }
  
  //fParticle = 3;
  
  std::cout<<"\n UserCreateOutputObject; PID = "<<fParticle<<" FB = "<<fFilterBit<<" harmonic = "<<gHarmonic<<"...\n"<<endl;
  
  
  PostData(1,fListHist);
  
}










//____________________________ Call Event by Event ___________________________________
void AliAnalysisTaskCMWPU2018dca::UserExec(Option_t*) {
 
  //std::cout<<" Info:UserExec() called ..!!!\n";






  
  Float_t stepCount = 0.5;

  fHistEventCount->Fill(stepCount);   //1
  stepCount++;





  fAOD = dynamic_cast <AliAODEvent*> (InputEvent());
  fESD = dynamic_cast <AliESDEvent*> (InputEvent());
  if(!(fESD || fAOD)) {
    printf("ERROR: fESD & fAOD not available\n");
    return;
  }

  //std::cout<<" Info:UserExec()  AOD Exist ..!!!\n";

  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent){
    printf("ERROR: fVevent not available\n");
    return;
  }





  fHistEventCount->Fill(stepCount); //2
  stepCount++;

  
  //fHistEventCount->Fill(stepCount); //2
  //stepCount++;


  

    

  //-------------- Vtx cuts ---------------
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t pVtxZ = -999;
  pVtxZ  = pVtx->GetZ();

  //std::cout<<" Info:UserExec()  AliVEvent Exist pVtxZ = "<<pVtxZ<<"..!!!\n";
    
  if(pVtxZ < fMinVzCut || pVtxZ > fMaxVzCut){
    return;
  }



  //std::cout<<" Info:UserExec()  Vz<10cm check ..!!!\n";
  
  fHistEventCount->Fill(stepCount); //3
  stepCount++;

  Float_t centrality = -99.0;
  Float_t centrV0M   = -99.0;
  //Float_t centrCL1   = -99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");  // Must never comment this
  if(!fMultSelection) {
    printf("\n...**ERROR**...\n UserExec() AliMultSelection object not found\n Status:Quit!! \n");
    exit(1);
  }


  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  //centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");

 
  centrality = centrV0M;  // This Is Always Default, changes below for other options:
  
  if(sCentrEstimator=="CL0"){
    centrality = fMultSelection->GetMultiplicityPercentile("CL0");;
  }
  else if(sCentrEstimator=="CL1"){
    centrality = fMultSelection->GetMultiplicityPercentile("CL1");;
  }
  else if(sCentrEstimator=="V0C"){
    centrality = fMultSelection->GetMultiplicityPercentile("V0C");
  }
  else if(sCentrEstimator=="V0A"){
    centrality = fMultSelection->GetMultiplicityPercentile("V0A");
  }
  else if(sCentrEstimator=="TRK"){
    centrality = fMultSelection->GetMultiplicityPercentile("TRK");
  }
  else{
    centrality = centrV0M;
  }

  
  fCentDistBeforCut->Fill(centrality);
  




  //Pile up cut function called--------------------------

  Bool_t kPileupEvent = kFALSE;

  kPileupEvent = CheckEventIsPileUp2018(fAOD);
  //if (!bSkipPileUpCut)
  if(kPileupEvent)  return;
    
  
  /*
  AliTimeRangeCut  *fTimeRangeCut;
  fTimeRangeCut = new AliTimeRangeCut;
  fTimeRangeCut->InitFromEvent(fAOD);
  if (fTimeRangeCut->CutEvent(fAOD))
    return;
  */



  //////////////////////////////////////////////////



  if(centrality<fCentralityMin || centrality>fCentralityMax){ 
    return;
  }


  //std::cout<<" Info:UserExec()  Centrality checkd..!!!\n";
  
  ///----> Get 0-10 index of centrality:
  Int_t cent10bin = -99, iCent = -99;

  if(centrality<5.0) {
    cent10bin  = 0; 
  }
  else if(centrality>=5.0 && centrality<10){
    cent10bin  = 1;
  }
  else if(centrality>=10.0) {
    cent10bin = abs(centrality/10.0)+1;
  }

  iCent = cent10bin;


  //Centrality array index for NUA correcion
  // Remove Centrality dependence..
  
  /*Int_t cForNUA = 0;  
  
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
  }*/









  

  Int_t ntracks = fAOD->GetNumberOfTracks();
  if(ntracks < 4) return;                        // Minimum 4 tracks per event. 

  fHistEventCount->Fill(stepCount); //7
  stepCount++;

  //std::cout<<" Info:UserExec()  minimum ntracks checkd..!!!\n";
  
  

  //////----> Get Magnetic field and RunNo.---------
  // Float_t fMagField = fAOD->GetMagneticField();
  // const Int_t QAindex = (fMagField > 0) ? 1 : 0;
  Int_t runNumber = fAOD->GetRunNumber();
  //------------------------------------------------

  /*
  if(fListTRKCorr) GetMCCorrectionHist(runNumber,centrality);  //use centrality dependent MC efficiency (Temporary!!)

  

  if(fListNUACorr){
     GetNUACorrectionHist(runNumber,0);         //Charge
     GetNUACorrectionHist(runNumber,fParticle); //1=pion, 2=kaon, 3=proton
  }
  
  if(fListV0MCorr){
     GetV0MCorrectionHist(runNumber,0);         //Charge
     //GetEVNTWGTCorrectionHist(runNumber,fParticle); //1=pion, 2=kaon, 3=proton
  }
  */




  
  Float_t fMultTPCFull = 0;  // TPC mult estimate
  Float_t fMultGlobal  = 0;  // global track multiplicity
  

  Float_t trkPt=0,trkPhi=0,trkEta=0,trkDCAxy=0.0, trkDCAz=0.0, trkDCAxyA=0.0, trkDCAzA=0.0;
  Float_t trkChi2=0,trkdEdx=0,trkWgt=1.0;
  Double_t dcaxys=0.0;
  Int_t   trkChrg=0, trkTpcNC=0;
  ////PID variables:
  Double_t nSigTOFpion=-99, nSigTPCpion=-99;
  Double_t nSigTOFkaon=-99, nSigTPCkaon=-99;
  Double_t nSigTOFprot=-99, nSigTPCprot=-99;  
  //Bool_t   bTOFmatch= kFALSE;
  Bool_t   isItPion = kFALSE, isItKaon= kFALSE, isItProt= kFALSE;
	  
  ////User's variable:
  Double_t fSumTPCQn2xNeg=0, fSumTPCQn2yNeg=0, fSumTPCQn2xPos=0, fSumTPCQn2yPos=0;
  Double_t fSumTPCQn2xNegWhole=0, fSumTPCQn2yNegWhole=0, fSumTPCQn2xPosWhole=0, fSumTPCQn2yPosWhole=0;
  Double_t fSumTPCQn2xNegChNeg=0, fSumTPCQn2yNegChNeg=0, fSumTPCQn2xPosChNeg=0, fSumTPCQn2yPosChNeg=0;
  Double_t fSumWgtEtaNeg=0, fSumWgtEtaPos=0;
  Double_t fSumWgtEtaNegWhole=0, fSumWgtEtaPosWhole=0;
  Double_t fSumWgtEtaNegChNeg=0, fSumWgtEtaPosChNeg=0;
  Double_t fNumOfPos = 0;
  Double_t fNumOfNeg = 0;
  Double_t ptWgtMC = 1.0, WgtNUA = 1.0;
  Int_t ptBinMC = 1, iBinNUA = 1;
  Double_t binCont = 1.0;
  
  Float_t fWgtEvent = 1.0; //Event Weight if Any.

  /// TO be passed as Argument:
  Double_t gPsiN = gHarmonic;

  
  //std::cout<<" Info:UserExec()  Starting track loop 1..!!!\n";
  
  ////////-----> Starting 1st track Loop -----------
  
  for(Int_t iTrack = 0; iTrack < ntracks; iTrack++) { 

    AliAODTrack* AODtrack = dynamic_cast <AliAODTrack*> (fVevent->GetTrack(iTrack));
    if(!AODtrack) continue;
    
    if(AODtrack->TestFilterBit(fFilterBit)) {  //// Only use FB tracks. 

      trkPt    = AODtrack->Pt();
      trkPhi   = AODtrack->Phi();
      trkEta   = AODtrack->Eta();
      trkChrg  = AODtrack->Charge();
      trkChi2  = AODtrack->Chi2perNDF();
      trkTpcNC = AODtrack->GetTPCNcls();
      trkDCAxy=  AODtrack->DCA();
      trkDCAz=   AODtrack->ZAtDCA();

            
      /// This Next function is called After Filter bit is validated!! (Otherwise code breaks!)
      trkdEdx  = AODtrack->GetDetPid()->GetTPCsignal();

      if((trkPt <= 10) && (trkPt >= 0.2) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= fChi2) && TMath::Abs(trkChrg)) {


	if ((TMath::Abs(trkDCAxy)==999)||(TMath::Abs(trkDCAz)==999))
	  {
	    Double_t    bval[2] = {-99., -99.};
	    Double_t    bCov[3] = {-99., -99., -99.};
	    AliAODTrack copy(*AODtrack);
	    if(copy.PropagateToDCA(fVevent->GetPrimaryVertex(), fVevent->GetMagneticField(), 100., bval, bCov) && TMath::Abs(bval[0]) < 0.3 &&  TMath::Abs(bval[1]) < 0.3){
	      fMultGlobal++;
	    }

	    trkDCAzA = bval[1];
	    trkDCAxyA = bval[0];
          }
	else
	  {
            trkDCAzA=trkDCAz;
            trkDCAxyA=trkDCAxy;
	  }

	dcaxys=7*(0.0026+(0.005/(TMath::Power(trkPt,1.01))));


	if (TMath::Abs(trkDCAzA)<fDCAzMax && TMath::Abs(trkDCAxyA)<dcaxys)
	  {


	    fHdcaz->Fill(trkDCAzA);
	    fHdcaxy->Fill(trkDCAxyA);

	          
	//---------->  Here I do All my track level analysis:
    	
	if(trkChrg > 0){	  
	  fNumOfPos += 1;
	}
	else{
	  fNumOfNeg += 1;
	}
	
	  }	
		//<---------- User track analysis Done---------------

      }// when all trackCuts applied      
    }//-------> when FB is validated


  
  }///------> 1st Track loop Ends here.<--------




  if(iCent > 9) return;

  fHistEventCount->Fill(stepCount); //8
  stepCount++;

  
  fHistEventCount->Fill(stepCount); //9
  stepCount++;
  

  
  
  Float_t fAchrgNet = (fNumOfPos - fNumOfNeg)/(fNumOfPos + fNumOfNeg); // Efficiency & NUA Corrected!
  fHistAChrgVsCent->Fill(centrality, fAchrgNet, fWgtEvent);
  













  
  ////////-----> Starting 2nd track Loop -----------
  /// Track loop variable (not used in previous loop)
 


  
  for(Int_t iTrack = 0; iTrack < ntracks; iTrack++) { 

    AliAODTrack* AODtrack = dynamic_cast <AliAODTrack*> (fVevent->GetTrack(iTrack));
    if(!AODtrack) continue;
    
    if(AODtrack->TestFilterBit(fFilterBit)) {  //// Only use FB tracks. 

      trkPt    = AODtrack->Pt();
      trkPhi   = AODtrack->Phi();
      trkEta   = AODtrack->Eta();
      trkChrg  = AODtrack->Charge();
      trkChi2  = AODtrack->Chi2perNDF();
      trkTpcNC = AODtrack->GetTPCNcls();
      trkDCAxy=  AODtrack->DCA();
      trkDCAz=   AODtrack->ZAtDCA();
            
      /// This Next function is called After Filter bit is validated!! (Otherwise code breaks!)
      trkdEdx  = AODtrack->GetDetPid()->GetTPCsignal();  

      

      //Apply track cuts here:
      if((trkPt <= fMaxPtCut) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= fChi2) && TMath::Abs(trkChrg)) {

   

	if ((TMath::Abs(trkDCAxy)==999)||(TMath::Abs(trkDCAz)==999))
	  {
	    Double_t    bval[2] = {-99., -99.};
	    Double_t    bCov[3] = {-99., -99., -99.};
	    AliAODTrack copy(*AODtrack);
	    if(copy.PropagateToDCA(fVevent->GetPrimaryVertex(), fVevent->GetMagneticField(), 100., bval, bCov) && TMath::Abs(bval[0]) < 0.3 &&  TMath::Abs(bval[1]) < 0.3){
	      fMultGlobal++;
	    }

	    trkDCAzA = bval[1];
	    trkDCAxyA = bval[0];
          }
	else
	  {
            trkDCAzA=trkDCAz;
            trkDCAxyA=trkDCAxy;
	  }

	dcaxys=7*(0.0026+(0.005/(TMath::Power(trkPt,1.01))));


	if (TMath::Abs(trkDCAzA)<fDCAzMax && TMath::Abs(trkDCAxyA)<dcaxys)
	  
	  {


    
	//---------->  Here I do All my track level analysis:

	
	isItPion = kFALSE;
	isItKaon = kFALSE;
	isItProt = kFALSE;
     
	///=========> Get TPC/TOF nSigma for PID




	//----- Pion
	if(fParticle==1){
	  nSigTPCpion = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kPion);
	  nSigTOFpion = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kPion);
          Double_t nSigRMSTPCTOFpion = TMath::Sqrt(nSigTPCpion*nSigTPCpion + nSigTOFpion*nSigTOFpion);

	
	  if(trkPt<=0.5 && TMath::Abs(nSigTPCpion)<=fNSigmaTPCCut){
	    isItPion = kTRUE;
	  }
	  // Using Circular cut for Pion: 	  
	  else if(trkPt>0.5 && TMath::Abs(nSigRMSTPCTOFpion)<=fNSigmaTOFCut){
	  isItPion = kTRUE;
	  }
	}

	//----- Kaon
	else if(fParticle==2){
	  nSigTPCkaon = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kKaon);	
	  nSigTOFkaon = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kKaon);
          Double_t nSigRMSTPCTOFkaon = TMath::Sqrt(nSigTPCkaon*nSigTPCkaon + nSigTOFkaon*nSigTOFkaon);


	  if(trkPt<=0.45 && TMath::Abs(nSigTPCkaon)<=fNSigmaTPCCut){
	    isItKaon = kTRUE;
	  }
	  else if(trkPt>0.45 && TMath::Abs(nSigRMSTPCTOFkaon)<=fNSigmaTOFCut){
	    isItKaon = kTRUE;
	  }
	}


	else if(fParticle==3){
	  nSigTPCprot = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kProton);    
	  nSigTOFprot = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kProton);
          Double_t nSigRMSTPCTOFprot = TMath::Sqrt(nSigTPCprot*nSigTPCprot + nSigTOFprot*nSigTOFprot);
      
	  //----- Proton 
	  if(trkPt<=0.6 && TMath::Abs(nSigTPCprot)<=fNSigmaTPCCut){
	    isItProt = kTRUE;
	    if(trkChrg>0 && trkPt<0.4) isItProt = kFALSE;  
	  }
	  else if(trkPt>0.6 && TMath::Abs(nSigRMSTPCTOFprot)<=fNSigmaTOFCut){
	    isItProt = kTRUE;
	  }
	}



    	//-------- PID selection is done ---------

	ptWgtMC = 1.0;
	WgtNUA  = 1.0;
        
	
	//fHistEtaPhiAfterCut->Fill(trkPhi,trkEta,trkWgt);

	
	
	if(trkChrg > 0){	  
	  
	  if(fParticle==0)
	    fHFillNUAPosPID->Fill(pVtxZ,trkPhi,trkEta);
	    //fHFillNUAPosPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);	 
	    	  
	  if(fParticle==1 && isItPion){
	    //fHFillNUAPosPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	    fHFillNUAPosPID->Fill(pVtxZ,trkPhi,trkEta);
	  }
	  if(fParticle==2 && isItKaon){
	    //fHFillNUAPosPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	    fHFillNUAPosPID->Fill(pVtxZ,trkPhi,trkEta);
	  }
	  if(fParticle==3 && isItProt){
	    //fHFillNUAPosPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	    fHFillNUAPosPID->Fill(pVtxZ,trkPhi,trkEta);
	  }

	  
	  
	}///+ve Ch done	
	else{  //-Ve charge
		  
	  if(fParticle==0)
	    fHFillNUANegPID->Fill(pVtxZ,trkPhi,trkEta);
	    //fHFillNUANegPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	    
	  if(fParticle==1 && isItPion){
	    //fHFillNUANegPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	    fHFillNUANegPID->Fill(pVtxZ,trkPhi,trkEta);
	  }
	  if(fParticle==2 && isItKaon){
	    //fHFillNUANegPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	    fHFillNUANegPID->Fill(pVtxZ,trkPhi,trkEta);
	  }
	  if(fParticle==3 && isItProt){
	    //fHFillNUANegPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	    fHFillNUANegPID->Fill(pVtxZ,trkPhi,trkEta);
	  }
	  
	}/// if -ve Particle
	//----------- v2 vs Ach filled ---------
         
	  }
      }//----> with all trackCuts applied      

    }//-----> if FB is validated
    
  }///------> 2nd Track loop Ends here.<--------
 

  /// For cumulant method:



  
  //Last lines in Event loop
  fCentDistAfterCut->Fill(centrality);
  fHistEventCount->Fill(14.5); //final Event count.
  PostData(1,fListHist);

  // std::cout<<" Info:UserExec()  Call Finished ..!!!\n";
}//---------------- UserExec ----------------------








Bool_t AliAnalysisTaskCMWPU2018dca::CheckEventIsPileUp2018(AliAODEvent *faod) {


  /*
  TF1 *fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);
    

  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  TF1 *fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);
    

  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  TF1 *fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  TF1 *fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);
    

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  TF1 *fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);




  Bool_t BisPileup=kFALSE;


  Double_t centrV0M=300;
  Double_t centrCL1=300;
  Double_t centrCL0=300;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
 
  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }


  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");



  Int_t nITSClsLy0 = faod->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = faod->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;
    

  AliAODTracklets* aodTrkl = (AliAODTracklets*)faod->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();


  const Int_t nTracks = faod->GetNumberOfTracks();
    

  Int_t multTrk = 0;
    

  for (Int_t it = 0; it < nTracks; it++) {
        

    AliAODTrack* aodTrk = (AliAODTrack*)faod->GetTrack(it);
        

    if (!aodTrk){
      delete aodTrk;
      continue;
    }


        

    if (aodTrk->TestFilterBit(32)){
      if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
	multTrk++;
            
    }     

  }




    

  AliAODVZERO* aodV0 = faod->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  Int_t tpcClsTot = faod->GetNumberOfTPCClusters();
  Float_t nclsDif = Float_t(tpcClsTot) - (60932.9 + 69.2897*multV0Tot - 0.000217837*multV0Tot*multV0Tot);


  if (centrCL0 < fCenCutLowPU->Eval(centrV0M))
    BisPileup=kTRUE;

  if (centrCL0 > fCenCutHighPU->Eval(centrV0M))
    BisPileup=kTRUE;        


  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls))
    BisPileup=kTRUE;
        

  if (multV0On < fV0CutPU->Eval(multV0Tot))
    BisPileup=kTRUE;


  if (Float_t(multTrk) < fMultCutPU->Eval(centrV0M))
    BisPileup=kTRUE;


        

  if (((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0)
    BisPileup=kTRUE;


  if (faod->IsIncompleteDAQ())
    BisPileup=kTRUE;
        
  //if (nclsDif > 200000)//can be increased to 200000
  // BisPileup=kTRUE;



  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTrk);
  fHistCL0VsV0MBefore->Fill(centrV0M,centrCL0);


  if (!BisPileup)
    {
      fHistTPConlyVsV0MAfter->Fill(centrV0M,multTrk);
      fHistCL0VsV0MAfter->Fill(centrV0M,centrCL0);
    }

  */


  /*
  TF1 *fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  TF1 *fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);

  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};

  TF1 *fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  TF1 *fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  TF1 *fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);
 */



  Bool_t BisPileup=kFALSE;


  Double_t centrV0M=-99.0;
  Double_t centrCL1=-99.0;
  Double_t centrCL0=-99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");

  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }


  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");



  Int_t nITSClsLy0 = faod->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = faod->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;


  AliAODTracklets* aodTrkl = (AliAODTracklets*)faod->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();


  const Int_t nTracks = faod->GetNumberOfTracks();

    

  Int_t multTrk = 0;


  for (Int_t it = 0; it < nTracks; it++) {


    AliAODTrack* aodTrk = (AliAODTrack*)faod->GetTrack(it);

    if (!aodTrk){
      delete aodTrk;
      continue;
    }


    if (aodTrk->TestFilterBit(32)){
      if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;

    }

  }



  AliAODVZERO* aodV0 = faod->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  Int_t tpcClsTot = faod->GetNumberOfTPCClusters();
  Float_t nclsDif = Float_t(tpcClsTot) - (60932.9 + 69.2897*multV0Tot - 0.000217837*multV0Tot*multV0Tot);


  if (centrCL0 < fCenCutLowPU->Eval(centrV0M))
    {
      //cout<<"*****************hi i am in 1st**************************:"<<endl;
      BisPileup=kTRUE;
    }

  if (centrCL0 > fCenCutHighPU->Eval(centrV0M))
    {
      //cout<<"*****************hi i am in 2nd**************************:"<<endl;
      BisPileup=kTRUE;
    }


  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls))
    {
      //cout<<"*****************hi i am in 3rd**************************:"<<endl;
      BisPileup=kTRUE;
    }
        

  if (multV0On < fV0CutPU->Eval(multV0Tot))
    {
      //cout<<"*****************hi i am in 4th**************************:"<<endl;
      BisPileup=kTRUE;
    }


  if (Float_t(multTrk) < fMultCutPU->Eval(centrV0M))
    {
      //cout<<"*****************hi i am in 5th**************************:"<<endl;
      BisPileup=kTRUE;
    }


  if (((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0)
    BisPileup=kTRUE;


  if (faod->IsIncompleteDAQ())
    BisPileup=kTRUE;
        
  //if (nclsDif > 200000)//can be increased to 200000
  // BisPileup=kTRUE;



  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTrk);
  fHistCL0VsV0MBefore->Fill(centrV0M,centrCL0);


  if (!BisPileup)
    {
      fHistTPConlyVsV0MAfter->Fill(centrV0M,multTrk);
      fHistCL0VsV0MAfter->Fill(centrV0M,centrCL0);
    }

 return BisPileup; 

}














Bool_t AliAnalysisTaskCMWPU2018dca::CheckEventIsPileUp(AliAODEvent *faod) {




  cout<<"---------------------------hi i am in 2015 dataset-------------------------------------"<<endl;


  Bool_t BisPileup=kFALSE;


  Double_t centrV0M=300;
  Double_t centrCL1=300;
  Double_t centrCL0=300;
  Double_t centrTRK=300;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
 
  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }


  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
  centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");


  //-- pile-up a la Dobrin for LHC15o -----
  if(PileUpMultiVertex(faod)) {
    fHistPileUpCount->Fill(0.5);    
    BisPileup=kTRUE;
  }
  Int_t isPileup = faod->IsPileupFromSPD(3);
  if(isPileup != 0) {
    fHistPileUpCount->Fill(1.5);
    BisPileup=kTRUE;          
  }
  if(((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
    fHistPileUpCount->Fill(2.5);
    BisPileup=kTRUE;
  }
  if(faod->IsIncompleteDAQ())  {
    fHistPileUpCount->Fill(3.5);
    BisPileup=kTRUE;
  }
  if(fabs(centrV0M-centrCL1)> 5.0)  {//default: 7.5
    fHistPileUpCount->Fill(4.5);
    BisPileup=kTRUE;
  }



  // check vertex consistency
  const AliAODVertex* vtTrc = faod->GetPrimaryVertex();
  const AliAODVertex* vtSPD = faod->GetPrimaryVertexSPD();

  if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
    fHistPileUpCount->Fill(5.5);
    BisPileup=kTRUE;
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
    BisPileup=kTRUE;
  }




  Int_t multTPC = 0;
  Int_t multTPCAll = 0;
  Int_t multITSfb96 = 0;
  Int_t multITSfb32 = 0;

  Int_t multTPCFE = 0;
  Int_t multGlobal = 0;
  Int_t multTPCuncut = 0;

  Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks();

  const Int_t nTracks = faod->GetNumberOfTracks();

  for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
    //AliNanoAODTrack* track = dynamic_cast<AliNanoAODTrack*>(faod->GetTrack(iTracks));
    AliAODTrack* track = (AliAODTrack*)faod->GetTrack(iTracks);
    if(!track)  continue;
  
    //---------- old method -----------
    if(track->TestFilterBit(128))
      multTPCAll++;
    //if(track->TestFilterBit(96))
    //multITSfb96++;


    //----------------------------------
    if(track->TestFilterBit(1))  multTPCuncut++;
    if(track->TestFilterBit(32)) multITSfb32++;


    if(track->Pt()<0.2 || track->Pt()>10.0 || TMath::Abs(track->Eta())>0.8 || track->GetTPCNcls()<fTPCclustMin || track->GetTPCsignal()<10.0)
      continue;
    if(track->GetDetPid() && track->Chi2perNDF() > 0.2) multTPC++;
    if(track->TestFilterBit(1) && track->Chi2perNDF()>0.2)  multTPCFE++;
    if(!track->TestFilterBit(16) || track->Chi2perNDF()<0.1)   continue;
                
    Double_t b[2]    = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
                
    AliAODTrack copy(*track);
    Double_t magField = faod->GetMagneticField();
                
    if(magField!=0){     
      if(track->PropagateToDCA(faod->GetPrimaryVertex(), magField, 100., b, bCov) && TMath::Abs(b[0]) < 0.3 && TMath::Abs(b[1]) < 0.3) 
	multGlobal++;    
    }
  }




  //fixed for test:
  // Double_t fPileUpSlopeParm = 3.43;
  //Double_t fPileUpConstParm = 43;

  Double_t multESDTPCDif  = multEsd  - fPileUpSlopeParm*multTPCAll;
  

  Bool_t  bIsOutLier=kFALSE;
  if(multTPC < (-20.0+1.15*multGlobal) || multTPC > (200.+1.45*multGlobal)) { bIsOutLier = kTRUE;}

  //  fHistEventCount->Fill(stepCount); //4
  //stepCount++;





  fTPCvsGlobalTrkBefore->Fill(multGlobal,multTPC);
  //fHistTPCVsESDTrkBefore->Fill(multTPCAll,multEsd);   //A. Dobrin
  fHistTPConlyVsCL1Before->Fill(centrCL1,multTPCAll);
  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTPCAll);
  fHistGlobalVsV0MBefore->Fill(centrV0M, multGlobal);



  if(multESDTPCDif > fPileUpConstParm) { 
    fHistPileUpCount->Fill(7.5);
    BisPileup=kTRUE;
  }
  if(BisPileup==kFALSE) {
    if(!fMultSelection->GetThisEventIsNotPileup())  BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
    if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
    if(BisPileup)     fHistPileUpCount->Fill(9.5);
  }  


  if (!bIsOutLier){
    fTPCvsGlobalTrkAfter->Fill(multGlobal,multTPC);

    //fHistEventCount->Fill(stepCount); //5
    //stepCount++;

  }
    if(!BisPileup) {
      //fHistTPCVsESDTrkAfter->Fill(multTPCAll,multEsd);  
      fHistTPConlyVsCL1After->Fill(centrCL1,multTPCAll);
      fHistTPConlyVsV0MAfter->Fill(centrV0M,multTPCAll);
      fHistGlobalVsV0MAfter->Fill(centrV0M, multGlobal);  
    //fHistEventCount->Fill(stepCount); //6
    //stepCount++
    }


  return BisPileup; 
} //-------pile up function ------









Bool_t AliAnalysisTaskCMWPU2018dca::PileUpMultiVertex(const AliAODEvent* faod)
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




double AliAnalysisTaskCMWPU2018dca::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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











void AliAnalysisTaskCMWPU2018dca::SetupEventAndTaskConfigInfo(){




  fHistPileUpCount = new TH1F("fHistPileUpCount", "fHistPileUpCount", 15, 0., 15.);
  fHistPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fHistPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fHistPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultComb08");
  fHistPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fHistPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>5.0");
  fHistPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fHistPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  Int_t puConst = fPileUpConstParm;
  //Int_t puConst = 3.0;
  fHistPileUpCount->GetXaxis()->SetBinLabel(8,Form("multESDTPCDif>%d",puConst));
  fHistPileUpCount->GetXaxis()->SetBinLabel(9,Form("multGlobTPCDif>%d",puConst));
  fHistPileUpCount->GetXaxis()->SetBinLabel(10,"PileUpMultSelTask");
  fListHist->Add(fHistPileUpCount);







  fHistEventCount = new TH1F("fHistEventCount","Event counts",15,0,15);
  fHistEventCount->GetXaxis()->SetBinLabel(1,"Called UserExec()");
  fHistEventCount->GetXaxis()->SetBinLabel(2,"Called Exec()");
  //fHistEventCount->GetXaxis()->SetBinLabel(3,"AOD Exist");
  fHistEventCount->GetXaxis()->SetBinLabel(3,"Vz < 10");
  fHistEventCount->GetXaxis()->SetBinLabel(4,"TPC vs Global");
  fHistEventCount->GetXaxis()->SetBinLabel(5,"TPC128 vs ESD");
  fHistEventCount->GetXaxis()->SetBinLabel(6,"Cent vs TPC");
  //fHistEventCount->GetXaxis()->SetBinLabel(8,Form("%2.0f<Cent<%2.0f",fCentralityMin,fCentralityMax));
  fHistEventCount->GetXaxis()->SetBinLabel(7,"noAODtrack > 2 ");  
  fHistEventCount->GetXaxis()->SetBinLabel(8,"cent < 90");  
  fHistEventCount->GetXaxis()->SetBinLabel(9,"mult eta+/- > 2");
  fHistEventCount->GetXaxis()->SetBinLabel(15,"Survived Events");
  fListHist->Add(fHistEventCount);
  //fHistEventCount->Fill(1);
}//----------SetupEventAndTaskConfigInfo-----------




Int_t AliAnalysisTaskCMWPU2018dca::GetCentralityScaled0to10(Double_t fCent){

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
 
}//------------GetCentralityScaled0to10------------

void AliAnalysisTaskCMWPU2018dca::GetNUACorrectionHist(Int_t run, Int_t kParticleID)
{

  //// Do not allocate memory for Reading NUA correction histograms.
  //// Because this macro only fills NUA histograms. 

  //// Removed the unnecessary Histograms!!
  
  // else {
  //   printf("\n ******** Warning: No NUA Correction File/List...!! \n Run= %d, Use NUA Wgt = 1.0 ********* \n",run);
  // }
}







void AliAnalysisTaskCMWPU2018dca::GetV0MCorrectionHist(Int_t run, Int_t kParticleID)
{

  if(fListV0MCorr){
    //cout<<"*****************AT LAST I GOT INTO EVNTWGT CORRECTION FUNCTION******************************"<<endl;
    if(kParticleID==0){ //charge
    fHCorrectEVNTWGTChrg = (TH1F *) fListV0MCorr->FindObject(Form("hwgtCharge_Run%d",run)); 	
    

    //cout<<"Run no:"<<run<<endl;

}
  }
    else{
      // cout<<"*****************AT LAST I GOT INTO EVNTWGT CORRECTION FUNCTION BUT LIST NOT PRESENT******************************"<<endl;
      fHCorrectEVNTWGTChrg=NULL;
    }
}
    
 



////---------- SetUp Tracking Efficiency Correction Map ---------------

void AliAnalysisTaskCMWPU2018dca::GetMCCorrectionHist(Int_t run,Float_t centr){

  if(fListTRKCorr) {
    //cout<<"\n =========> Info: Found TList with MC Tracking Corr Histograms <=========== "<<endl;
    /// Default: Centrality Independent MC efficiency:


    /// Do Nothing, as this Task only Fills NUA histrograms.
    /*
    fHCorrectMCposChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgPos");
    fHCorrectMCposPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionPos");
    fHCorrectMCposKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonPos");
    fHCorrectMCposProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtPos");

    fHCorrectMCnegChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgNeg");
    fHCorrectMCnegPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionNeg");
    fHCorrectMCnegKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonNeg");
    fHCorrectMCnegProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtNeg");
    */
    

    //for(int i=0;i<10;i++) {
    //fFB_Efficiency_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_%d",i));
    //}
  }
  // else if(!fListTRKCorr){
  //   std::cout<<"\n\n !!!!**** Warning : No FB Efficiency Correction, run = "<<run<<"..!!!**** \n using MC TrkWgt = 1.0 \n"<<std::endl;
  // }
}



















