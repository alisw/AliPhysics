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

/* $Id: AliAnalysisTaskCVE.cxx  Rihan Haque, 18/09/2019 (ver1) $ */

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
#include "AliAnalysisTaskv2pt.h"

using std::cout;
using std::endl;
using std::vector;


ClassImp(AliAnalysisTaskv2pt)

AliAnalysisTaskv2pt::AliAnalysisTaskv2pt(const char *name): AliAnalysisTaskSE(name),
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
  fMaxevpt(2.0),
  fEtaGapNeg(-0.1),
  fEtaGapPos(0.1),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),    
  fdEdxMin(10.0),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  puswitch(0),
  sCentrEstimator("V0M"),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fHCorrectMCposChrg(NULL),
  fHCorrectMCposPion(NULL),
  fHCorrectMCposKaon(NULL),
  fHCorrectMCposProt(NULL),
  fHCorrectMCnegChrg(NULL),
  fHCorrectMCnegPion(NULL),
  fHCorrectMCnegKaon(NULL),
  fHCorrectMCnegProt(NULL),      

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
  fHistv2cumCentChrgAll(NULL),
  fHistv2cumCentChrgAllNeg(NULL),
  fHistPileUpCount(NULL),  

  fHCorrectNUAposChrg(NULL),  
  fHCorrectNUAnegChrg(NULL),  
  fHCorrectNUAposPion(NULL),  
  fHCorrectNUAnegPion(NULL),  
  fHCorrectNUAposKaon(NULL),  
  fHCorrectNUAnegKaon(NULL),  
  fHCorrectNUAposProt(NULL),  
  fHCorrectNUAnegProt(NULL),

  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),   

  fHistEventCount(NULL)
{

  
  for(int i=0;i<1;i++){
    for(int j=0;j<10;j++){

      fHistv2AchChrgPos[i][j] = NULL;
      fHistv2AchPionPos[i][j] = NULL;
      fHistv2AchKaonPos[i][j] = NULL;
      fHistv2AchProtPos[i][j] = NULL;

      fHistv2AchChrgNeg[i][j] = NULL;
      fHistv2AchPionNeg[i][j] = NULL;
      fHistv2AchKaonNeg[i][j] = NULL;
      fHistv2AchProtNeg[i][j] = NULL;
      
    }
  }
  
  
 
   for(int i=0; i<10; i++){
     fHistv2cumAchChrgAll[i] = NULL;
   }
  
  
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskv2pt::AliAnalysisTaskv2pt():
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
  fMaxevpt(2.0),  
  fEtaGapNeg(-0.1),
  fEtaGapPos(0.1),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),    
  fdEdxMin(10.0),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  fPileUpSlopeParm(3.43),
  fPileUpConstParm(43),
  puswitch(0),
  sCentrEstimator("V0M"),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fHCorrectMCposChrg(NULL),
  fHCorrectMCposPion(NULL),
  fHCorrectMCposKaon(NULL),
  fHCorrectMCposProt(NULL),
  fHCorrectMCnegChrg(NULL),
  fHCorrectMCnegPion(NULL),
  fHCorrectMCnegKaon(NULL),
  fHCorrectMCnegProt(NULL),    
 
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
  fHistv2cumCentChrgAll(NULL),
  fHistv2cumCentChrgAllNeg(NULL),
  fHistPileUpCount(NULL),  

  fHCorrectNUAposChrg(NULL),  
  fHCorrectNUAnegChrg(NULL),  
  fHCorrectNUAposPion(NULL),  
  fHCorrectNUAnegPion(NULL),  
  fHCorrectNUAposKaon(NULL),  
  fHCorrectNUAnegKaon(NULL),  
  fHCorrectNUAposProt(NULL),  
  fHCorrectNUAnegProt(NULL),


  fSPDCutPU(NULL),
  fV0CutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),   
        

  fHistEventCount(NULL)
{
  
  for(int i=0;i<1;i++){
    for(int j=0;j<10;j++){
      //for (int k=0;k<1;k++){
      fHistv2AchChrgPos[i][j] = NULL;
      fHistv2AchPionPos[i][j] = NULL;
      fHistv2AchKaonPos[i][j] = NULL;
      fHistv2AchProtPos[i][j] = NULL;

      fHistv2AchChrgNeg[i][j] = NULL;
      fHistv2AchPionNeg[i][j] = NULL;
      fHistv2AchKaonNeg[i][j] = NULL;
      fHistv2AchProtNeg[i][j] = NULL;
      
    }
  }
    
  for(int i=0; i<10; i++){
    fHistv2cumAchChrgAll[i] = NULL;
  }
  
  
  //Not needed for Empty Constructor:
  //DefineInput(0,TChain::Class());
  //DefineOutput(1,TList::Class());
}
  
//__________________ destructor ___________________
AliAnalysisTaskv2pt::~AliAnalysisTaskv2pt()
{
  if(fListHist)      delete fListHist;  
  if(fAnalysisUtil)  delete fAnalysisUtil;   // because its 'new' !!

  //Delete the clones
  if(fListTRKCorr)  delete fListTRKCorr;
  if(fListNUACorr)  delete fListNUACorr;
  if(fListV0MCorr)  delete fListV0MCorr;
      

  if(fSPDCutPU)     delete fSPDCutPU;
  if(fV0CutPU)      delete fV0CutPU;
  if(fMultCutPU)    delete fMultCutPU;
  if(fCenCutLowPU)  delete fCenCutLowPU; 
  if(fCenCutHighPU) delete fCenCutHighPU; 
  

}










//________________ Define Histograms _______________
void AliAnalysisTaskv2pt::UserCreateOutputObjects()
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
  
  

  Int_t gMaxGlobalmult  = 4000;
  Int_t gMaxTPCcorrmult = 5000;
  Int_t gMaxESDtracks   = 20000;





  fHistGlobalVsV0MBefore = new TH2F("fHistGlobalVsV0MBefore","Before;Cent(V0M);Global",100,0,100,500,0,gMaxGlobalmult);
  fListHist->Add(fHistGlobalVsV0MBefore);
  fHistGlobalVsV0MAfter  = new TH2F("fHistGlobalVsV0MAfter"," After; Cent(V0M);Global",100,0,100,500,0,gMaxGlobalmult);
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
  
  fHistTPCVsESDTrkBefore = new TH2F("fHistTPCVsESDTrkBefore","Before; TPC1; ESD trk",500,0,gMaxTPCcorrmult,200,0,gMaxESDtracks);
  fListHist->Add(fHistTPCVsESDTrkBefore);
  fHistTPCVsESDTrkAfter  = new TH2F("fHistTPCVsESDTrkAfter"," After;  TPC1; ESD trk",500,0,gMaxTPCcorrmult,200,0,gMaxESDtracks);
  fListHist->Add(fHistTPCVsESDTrkAfter);

  
  //----------- User's histograms: --------------
  Char_t  name[100];
  Char_t title[100];
 

  Double_t centRange[11] = {0,5,10,20,30,40,50,60,70,80,90};

  
  const int NPT=21;
  Double_t PT[NPT+1]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.2,1.4,1.6,1.8,2.0};

  
 
  // v2 vs Ach
  for(int i=0;i<1;i++){
    for(int j=0;j<10;j++){
      //for (int k=0;k<1;k++){
      ////Charge:
      sprintf(name,"fHistv2AchChrgPos_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchChrgPos[i][j] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchChrgPos[i][j]->Sumw2();
      fListHist->Add(fHistv2AchChrgPos[i][j]);
      sprintf(name,"fHistv2AchChrgNeg_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchChrgNeg[i][j] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchChrgNeg[i][j]->Sumw2();
      fListHist->Add(fHistv2AchChrgNeg[i][j]);      

      //// Pion:
      sprintf(name,"fHistv2AchPionPos_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchPionPos[i][j] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchPionPos[i][j]->Sumw2();
      fListHist->Add(fHistv2AchPionPos[i][j]);
      sprintf(name,"fHistv2AchPionNeg_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchPionNeg[i][j] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchPionNeg[i][j]->Sumw2();
      fListHist->Add(fHistv2AchPionNeg[i][j]);      
 
      //// Kaon:
      sprintf(name,"fHistv2AchKaonPos_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchKaonPos[i][j] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchKaonPos[i][j]->Sumw2();
      fListHist->Add(fHistv2AchKaonPos[i][j]);
      sprintf(name,"fHistv2AchKaonNeg_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchKaonNeg[i][j] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchKaonNeg[i][j]->Sumw2();
      fListHist->Add(fHistv2AchKaonNeg[i][j]);      

      //// Proton:
      sprintf(name,"fHistv2AchProtPos_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchProtPos[i][j] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchProtPos[i][j]->Sumw2();
      fListHist->Add(fHistv2AchProtPos[i][j]);
      sprintf(name,"fHistv2AchProtNeg_Method%d_Cent%d",i,j);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchProtNeg[i][j] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchProtNeg[i][j]->Sumw2();
      fListHist->Add(fHistv2AchProtNeg[i][j]);  

         
      }
  }
  
  
 

  Int_t gCentForNUA[6] = {0,5,10,20,40,90};
  Char_t cpid[10];

 
  if(fParticle==1) sprintf(cpid,"Pion,Id %d",fParticle);
  else if(fParticle==2) sprintf(cpid,"Kaon,Id %d",fParticle);
  else if(fParticle==3) sprintf(cpid,"Prot,Id %d",fParticle);
  else  sprintf(cpid,"Charge,Id %d",fParticle);

 
  
  for(int i=0; i<10; i++){
    ////Charge:
    sprintf(name,"fHistv2cumAchChrgAllQcumCent%d",i);
    sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[i],centRange[i+1]);
    fHistv2cumAchChrgAll[i] = new TProfile(name,title,1,-10,10,"");
    fHistv2cumAchChrgAll[i] ->Sumw2();
    fListHist->Add(fHistv2cumAchChrgAll[i]);
  }


  sprintf(name,"fHistv2cumCentChrgAllQcumCent%d",0);
  sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[0],centRange[1]);
  fHistv2cumCentChrgAll = new TProfile(name,title,18,0,90,"");
  fHistv2cumCentChrgAll ->Sumw2();
  fListHist->Add(fHistv2cumCentChrgAll);

  sprintf(name,"fHistv2cumCentChrgAllQcumNegCent%d",0);
  sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[0],centRange[1]);
  fHistv2cumCentChrgAllNeg = new TProfile(name,title,18,0,90,"");
  fHistv2cumCentChrgAllNeg ->Sumw2();
  fListHist->Add(fHistv2cumCentChrgAllNeg);
  



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
  }
  else{
    std::cout<<"\n\n ******* WARNING No TList NUA Correction!!\n using NUAWgt = 1.0 \n "<<std::endl;
  }

  //fParticle = 3;
  
  std::cout<<"\n UserCreateOutputObject; PID = "<<fParticle<<" FB = "<<fFilterBit<<" harmonic = "<<gHarmonic<<"...\n"<<endl;
  
  
  PostData(1,fListHist);
  
}










//____________________________ Call Event by Event ___________________________________
void AliAnalysisTaskv2pt::UserExec(Option_t*) {
 
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


  fHistEventCount->Fill(stepCount);
  stepCount++;


  

  //-------------- Vtx cuts ---------------
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t pVtxZ = -999;
  pVtxZ  = pVtx->GetZ();
  

    
  if(pVtxZ < fMinVzCut || pVtxZ > fMaxVzCut){
    return;
  }

  
  fHistEventCount->Fill(stepCount); //3
  stepCount++;

  Float_t centrality = -99.0;
  Float_t centerV0M   = -99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");  // Must never comment this
  if(!fMultSelection) {
    printf("\n...**ERROR**...\n UserExec() AliMultSelection object not found\n Status:Quit!! \n");
    exit(1);
  }


  centerV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  //centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");

 
  centrality = centerV0M;  // This Is Always Default, changes below for other options:
  
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
    centrality = centerV0M;
  }

  
  fCentDistBeforCut->Fill(centrality);


  Bool_t kPileupEvent = kFALSE;

  if(puswitch==0)
  kPileupEvent = CheckEventIsPileUp2018(fAOD);
  else if(puswitch==1) 
  kPileupEvent = CheckEventIsPileUp(fAOD);
  //if (!bSkipPileUpCut)
  if(kPileupEvent)  return;


  
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
  








  

  Int_t ntracks = fAOD->GetNumberOfTracks();
  if(ntracks < 4) return;                        // Minimum 4 tracks per event. 

  fHistEventCount->Fill(stepCount); //4
  stepCount++;

  //std::cout<<" Info:UserExec()  minimum ntracks checkd..!!!\n";
  
  

  //////----> Get Magnetic field and RunNo.---------
  // Float_t fMagField = fAOD->GetMagneticField();
  // const Int_t QAindex = (fMagField > 0) ? 1 : 0;
  Int_t runNumber = fAOD->GetRunNumber();
  //------------------------------------------------

  if(fListTRKCorr) GetMCCorrectionHist(runNumber);

  if(fListNUACorr){
     GetNUACorrectionHist(runNumber,0);         //Charge
     GetNUACorrectionHist(runNumber,fParticle); //1=pion, 2=kaon, 3=proton
    //GetNUACorrectionHist(runNumber,1);
    //GetNUACorrectionHist(runNumber,2);
    //GetNUACorrectionHist(runNumber,3); 
  }





  
  Float_t fMultTPCFull = 0;  // TPC mult estimate
  Float_t fMultGlobal  = 0;  // global track multiplicity
  
  //Float_t fMultwRawFB  = 0;  // Uncorrected Multiplicity
  //Float_t fMultCorrFB  = 0;  // Corrected Multiplicity  
  //Int_t   gMultEtaNeg  = 0;
  //Int_t   gMultEtaPos  = 0;
  //Int_t   gMultEtaAll  = 0;

  Float_t trkPt=0,trkPhi=0,trkEta=0;
  Float_t trkChi2=0,trkdEdx=0,trkWgt=1.0;
  Int_t   trkChrg=0, trkTpcNC=0;
  ////PID variables:
  Double_t nSigTOFpion=-99, nSigTPCpion=-99;
  Double_t nSigTOFkaon=-99, nSigTPCkaon=-99;
  Double_t nSigTOFprot=-99, nSigTPCprot=-99;  
  //Bool_t   bTOFmatch= kFALSE;
  Bool_t   isItPion = kFALSE, isItKaon= kFALSE, isItProt= kFALSE;
	  
  ////User's variable:
  Double_t fSumTPCQn2xNeg=0, fSumTPCQn2yNeg=0, fSumTPCQn2xPos=0, fSumTPCQn2yPos=0;
  Double_t fSumTPCQn2xNegChNeg=0, fSumTPCQn2yNegChNeg=0, fSumTPCQn2xPosChNeg=0, fSumTPCQn2yPosChNeg=0;
  Double_t fSumWgtEtaNeg=0, fSumWgtEtaPos=0;
  Double_t fSumWgtEtaNegChNeg=0, fSumWgtEtaPosChNeg=0;
  Double_t fNumOfPos = 0;
  Double_t fNumOfNeg = 0;
  Double_t ptWgtMC = 1.0, WgtNUA = 1.0;
  Int_t    ptBinMC = 1, iBinNUA = 1;
  Double_t binCont = 1.0;
  
  Float_t fWgtEvent = 1.0; //Event Weight if Any.

  /// TO be passed as Argument:
  //Double_t fEtaGapNeg = -0.100;
  //Double_t fEtaGapPos =  0.100;
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

      
      /// This Next function is called After Filter bit is validated!! (Otherwise code breaks!)
      trkdEdx  = AODtrack->GetDetPid()->GetTPCsignal();  

      //Apply track cuts here:
      if((trkPt <= 10.0) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= 4.0) && TMath::Abs(trkChrg)) {

	//---------->  Here I do All my track level analysis:
    



	  

	//------> Get NUA weights for EP <----------

	//if(trkPt > 2.0) continue; ///// *********  Trk cut for Event Plane: 0.2 < pT < 2.0; *********
	
	WgtNUA = 1.0;
	ptWgtMC = 1.0;
	
	if(trkChrg>0){

	  if(fHCorrectMCposChrg){
	    ptBinMC = fHCorrectMCposChrg->FindBin(trkPt);    //Charge independent MC correction atm.
	    binCont = fHCorrectMCposChrg->GetBinContent(ptBinMC);
	    if(binCont!=0) ptWgtMC = 1.0/binCont;      
	  }
	  
	  if(fHCorrectNUAposChrg){
	    //iBinNUA = fHCorrectNUAposChrg[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUA  = fHCorrectNUAposChrg[cForNUA]->GetBinContent(iBinNUA);	    
	    iBinNUA = fHCorrectNUAposChrg->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUA  = fHCorrectNUAposChrg->GetBinContent(iBinNUA);
	  
	  }
	}
	else{
	  
	  if(fHCorrectMCnegChrg){
	    ptBinMC = fHCorrectMCnegChrg->FindBin(trkPt);    //Charge independent MC correction atm.	    
	    binCont = fHCorrectMCnegChrg->GetBinContent(ptBinMC);
	    if(binCont!=0) ptWgtMC = 1.0/binCont;  
	  }
	  
	  if(fHCorrectNUAnegChrg){
	    //iBinNUA = fHCorrectNUAnegChrg[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUA  = fHCorrectNUAnegChrg[cForNUA]->GetBinContent(iBinNUA);  
	    iBinNUA = fHCorrectNUAnegChrg->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUA  = fHCorrectNUAnegChrg->GetBinContent(iBinNUA);  
	  }
	}

	//RefMultRawFB++;
	//RefMultCorrFB += ptWgtMC;
	if(WgtNUA>1e3) WgtNUA = 1;   // For empty bins which cause Wgt = infinity !!

	trkWgt = WgtNUA*ptWgtMC;


	if(trkChrg > 0){	  
	  fNumOfPos += 1;
	}
	else{
	  fNumOfNeg += 1;
	}

	
	
	
	if (trkPt<fMaxevpt)
	  {
	    if(trkChrg > 0){	  
	    if(trkEta < fEtaGapNeg){
	  fSumTPCQn2xNeg += trkWgt*TMath::Cos(gPsiN*trkPhi);
	  fSumTPCQn2yNeg += trkWgt*TMath::Sin(gPsiN*trkPhi);
	  fSumWgtEtaNeg  += trkWgt;
	}
	    else if(trkEta > fEtaGapPos){
	  fSumTPCQn2xPos += trkWgt*TMath::Cos(gPsiN*trkPhi);
	  fSumTPCQn2yPos += trkWgt*TMath::Sin(gPsiN*trkPhi);
	  fSumWgtEtaPos  += trkWgt;
	}
	    }
	    else
	      {
	    if(trkEta < fEtaGapNeg){
	  fSumTPCQn2xNegChNeg += trkWgt*TMath::Cos(gPsiN*trkPhi);
	  fSumTPCQn2yNegChNeg += trkWgt*TMath::Sin(gPsiN*trkPhi);
	  fSumWgtEtaNegChNeg  += trkWgt;
	}
	    else if(trkEta > fEtaGapPos){
	  fSumTPCQn2xPosChNeg += trkWgt*TMath::Cos(gPsiN*trkPhi);
	  fSumTPCQn2yPosChNeg += trkWgt*TMath::Sin(gPsiN*trkPhi);
	  fSumWgtEtaPosChNeg  += trkWgt;
	}
	      }


	  }
	


	
	
	//<---------- User track analysis Done---------------

      }// when all trackCuts applied      
    }//-------> when FB is validated



    
  }///------> 1st Track loop Ends here.<--------




  if(iCent > 9) return;

  fHistEventCount->Fill(stepCount); //5
  stepCount++;


  if(fSumWgtEtaNeg <= 0 || fSumWgtEtaPos <= 0) return;
  
  fHistEventCount->Fill(stepCount); //6
  stepCount++;
  
  
  





  const int noofptbins=21;


  Double_t nSigRMSTPCTOFpion;
  Double_t nSigRMSTPCTOFkaon;
  Double_t nSigRMSTPCTOFprot;



  
  ////////-----> Starting 2nd track Loop -----------
  /// Track loop variable (not used in previous loop)
  Double_t sumQxTPCneg=0, sumQyTPCneg=0, sumQxTPCpos=0, sumQyTPCpos=0;
  Double_t sumWgtneg=0, sumWgtpos=0;
  Double_t uqRe=0, uqIm=0;

  //// Variable for Eta-Pos:
  Double_t sumQ2xChrgPosEtaPos[noofptbins]={0.0}, sumQ2yChrgNegEtaPos[noofptbins]={0.0};
  Double_t sumQ2xPionPosEtaPos[noofptbins]={0.0}, sumQ2yPionNegEtaPos[noofptbins]={0.0};
  Double_t sumQ2xKaonPosEtaPos[noofptbins]={0.0}, sumQ2yKaonNegEtaPos[noofptbins]={0.0};
  Double_t sumQ2xProtPosEtaPos[noofptbins]={0.0}, sumQ2yProtNegEtaPos[noofptbins]={0.0};
  Double_t sumQ2xChrgNegEtaPos[noofptbins]={0.0}, sumQ2yChrgPosEtaPos[noofptbins]={0.0};
  Double_t sumQ2xPionNegEtaPos[noofptbins]={0.0}, sumQ2yPionPosEtaPos[noofptbins]={0.0};
  Double_t sumQ2xKaonNegEtaPos[noofptbins]={0.0}, sumQ2yKaonPosEtaPos[noofptbins]={0.0};
  Double_t sumQ2xProtNegEtaPos[noofptbins]={0.0}, sumQ2yProtPosEtaPos[noofptbins]={0.0};
  
  Double_t  NumOfChrgPosEtaPos[noofptbins]={0.0},  NumOfChrgNegEtaPos[noofptbins]={0.0};    
  Double_t  NumOfPionPosEtaPos[noofptbins]={0.0},  NumOfPionNegEtaPos[noofptbins]={0.0};
  Double_t  NumOfKaonPosEtaPos[noofptbins]={0.0},  NumOfKaonNegEtaPos[noofptbins]={0.0};
  Double_t  NumOfProtPosEtaPos[noofptbins]={0.0},  NumOfProtNegEtaPos[noofptbins]={0.0};

  //// Now for Eta-Neg:
  Double_t sumQ2xChrgPosEtaNeg[noofptbins]={0.0}, sumQ2yChrgNegEtaNeg[noofptbins]={0.0};
  Double_t sumQ2xPionPosEtaNeg[noofptbins]={0.0}, sumQ2yPionNegEtaNeg[noofptbins]={0.0};
  Double_t sumQ2xKaonPosEtaNeg[noofptbins]={0.0}, sumQ2yKaonNegEtaNeg[noofptbins]={0.0};
  Double_t sumQ2xProtPosEtaNeg[noofptbins]={0.0}, sumQ2yProtNegEtaNeg[noofptbins]={0.0};
  Double_t sumQ2xChrgNegEtaNeg[noofptbins]={0.0}, sumQ2yChrgPosEtaNeg[noofptbins]={0.0};
  Double_t sumQ2xPionNegEtaNeg[noofptbins]={0.0}, sumQ2yPionPosEtaNeg[noofptbins]={0.0};
  Double_t sumQ2xKaonNegEtaNeg[noofptbins]={0.0}, sumQ2yKaonPosEtaNeg[noofptbins]={0.0};
  Double_t sumQ2xProtNegEtaNeg[noofptbins]={0.0}, sumQ2yProtPosEtaNeg[noofptbins]={0.0};

  Double_t  NumOfChrgPosEtaNeg[noofptbins]={0.0},  NumOfChrgNegEtaNeg[noofptbins]={0.0};   
  Double_t  NumOfPionPosEtaNeg[noofptbins]={0.0},  NumOfPionNegEtaNeg[noofptbins]={0.0};
  Double_t  NumOfKaonPosEtaNeg[noofptbins]={0.0},  NumOfKaonNegEtaNeg[noofptbins]={0.0};
  Double_t  NumOfProtPosEtaNeg[noofptbins]={0.0},  NumOfProtNegEtaNeg[noofptbins]={0.0};



  Double_t ptcenter[noofptbins]={0.0};
  Double_t lowpt[noofptbins]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.2,1.4,1.6,1.8};
  Double_t highpt[noofptbins]={0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.2,1.4,1.6,1.8,2.0};

  
  Double_t trkWgtPion=1.0, trkWgtKaon=1.0, trkWgtProt=1.0;
  Double_t WgtNUAPion=1.0, WgtNUAKaon=1.0, WgtNUAProt=1.0;
  Double_t ptWgtMCPion=1.0, ptWgtMCKaon=1.0, ptWgtMCProt=1.0;

  
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

      
      /// This Next function is called After Filter bit is validated!! (Otherwise code breaks!)
      trkdEdx  = AODtrack->GetDetPid()->GetTPCsignal();  

      //Apply track cuts here:
      if((trkPt <= fMaxPtCut) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= 4.0) && TMath::Abs(trkChrg)) {

	
	//---------->  Here I do All my track level analysis:

	
	isItPion = kFALSE;
	isItKaon = kFALSE;
	isItProt = kFALSE;
     
	///=========> Get TPC/TOF nSigma for PID


	//----- Pion
	if(fParticle==1){
	  nSigTPCpion = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kPion);
	  nSigTOFpion = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kPion);
	

	  nSigRMSTPCTOFpion = TMath::Sqrt(nSigTPCpion*nSigTPCpion + nSigTOFpion*nSigTOFpion);


	  if(trkPt<=0.5 && TMath::Abs(nSigTPCpion)<=fNSigmaTPCCut){
	    isItPion = kTRUE;
	  }
	  else if(trkPt>0.5 && trkPt<=10.0 && TMath::Abs(nSigRMSTPCTOFpion)<=fNSigmaTOFCut){
	    isItPion = kTRUE;
	  }
	}
	//----- Kaon
	else if(fParticle==2){
	  nSigTPCkaon = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kKaon);	
	  nSigTOFkaon = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kKaon);

	  nSigRMSTPCTOFkaon = TMath::Sqrt(nSigTPCkaon*nSigTPCkaon + nSigTOFkaon*nSigTOFkaon);

	  if(trkPt<=0.45 && TMath::Abs(nSigTPCkaon)<=fNSigmaTPCCut){
	    isItKaon = kTRUE;
	  }
	  else if(trkPt>0.45 && trkPt<=10.0 && TMath::Abs(nSigRMSTPCTOFkaon)<=fNSigmaTOFCut){
	    //else if(trkPt>0.45 && trkPt<=10.0 && TMath::Abs(nSigTPCkaon)<=fNSigmaTPCCut && TMath::Abs(nSigTOFkaon)<=fNSigmaTOFCut){
	    isItKaon = kTRUE;
	  }
	}
	else if(fParticle==3){
	  nSigTPCprot = fPIDResponse->NumberOfSigmasTPC(AODtrack, AliPID::kProton);    
	  nSigTOFprot = fPIDResponse->NumberOfSigmasTOF(AODtrack, AliPID::kProton);
      
	  //----- Proton 
	  if(trkPt<=0.8 && TMath::Abs(nSigTPCprot)<=fNSigmaTPCCut){
	    isItProt = kTRUE;
	    if(trkChrg>0 && trkPt<0.4) isItProt = kFALSE;  //Proton below 0.4 GeV has beam Pipe Contamination
	  }
	  else if(trkPt>0.8 && trkPt<=10.0 && TMath::Abs(nSigTPCprot)<=fNSigmaTPCCut && TMath::Abs(nSigTOFprot)<=fNSigmaTOFCut){  
	    isItProt = kTRUE;
	  }
	}
    	//-------- PID selection is done ---------
      
	ptWgtMC = 1.0;
	WgtNUA  = 1.0;
	
	trkWgtPion = 1.0;  trkWgtKaon = 1.0;  trkWgtProt = 1.0;
	WgtNUAPion = 1.0;  WgtNUAKaon = 1.0;  WgtNUAProt = 1.0; 
	ptWgtMCPion = 1.0; ptWgtMCKaon = 1.0; ptWgtMCProt = 1.0;
	

	
	if(trkChrg>0){
	  
	  ///Tracking:
	  if(fHCorrectMCposChrg){
	    ptBinMC = fHCorrectMCposChrg->FindBin(trkPt);    	    
	    binCont = fHCorrectMCposChrg->GetBinContent(ptBinMC);
	    if(binCont!=0) ptWgtMC = 1.0/binCont;  
	  }
	  if(fParticle==1 && isItPion && fHCorrectMCposPion){
	    ptBinMC = fHCorrectMCposPion->FindBin(trkPt);    
	    binCont = fHCorrectMCposPion->GetBinContent(ptBinMC);
	    if(binCont!=0)  ptWgtMCPion = 1.0/binCont; 	      
	  }
	  if(fParticle==2 && isItKaon && fHCorrectMCposKaon){
	    ptBinMC = fHCorrectMCposKaon->FindBin(trkPt);
	    binCont = fHCorrectMCposKaon->GetBinContent(ptBinMC);
	    if(binCont!=0)  ptWgtMCKaon = 1.0/binCont;
	  }
	  if(fParticle==3 && isItProt && fHCorrectMCposProt){
	    ptBinMC     = fHCorrectMCposProt->FindBin(trkPt);
	    binCont = fHCorrectMCposProt->GetBinContent(ptBinMC);
	    if(binCont!=0)  ptWgtMCProt = 1.0/binCont; 	    
	  }

	  ///NUA:
	  if(fHCorrectNUAposChrg){
	    //iBinNUA = fHCorrectNUAposChrg[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUA  = fHCorrectNUAposChrg[cForNUA]->GetBinContent(iBinNUA);
	    iBinNUA = fHCorrectNUAposChrg->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUA  = fHCorrectNUAposChrg->GetBinContent(iBinNUA);

}
	  if(fParticle==1 && isItPion && fHCorrectNUAposPion){
	    //iBinNUA     = fHCorrectNUAposPion[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUAPion  = fHCorrectNUAposPion[cForNUA]->GetBinContent(iBinNUA);
	    iBinNUA     = fHCorrectNUAposPion->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAPion  = fHCorrectNUAposPion->GetBinContent(iBinNUA);
	  }
	  if(fParticle==2 && isItKaon && fHCorrectNUAposKaon){
	    //iBinNUA     = fHCorrectNUAposKaon[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUAKaon  = fHCorrectNUAposKaon[cForNUA]->GetBinContent(iBinNUA);
	    iBinNUA     = fHCorrectNUAposKaon->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAKaon  = fHCorrectNUAposKaon->GetBinContent(iBinNUA);
	  }
	  if(fParticle==3 && isItProt && fHCorrectNUAposProt){
	    //iBinNUA     = fHCorrectNUAposProt[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUAProt  = fHCorrectNUAposProt[cForNUA]->GetBinContent(iBinNUA);
	    iBinNUA     = fHCorrectNUAposProt->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAProt  = fHCorrectNUAposProt->GetBinContent(iBinNUA);
	  }       	  
	}
	else{
	  ///Tracking:
	  if(fHCorrectMCnegChrg){
	    ptBinMC = fHCorrectMCnegChrg->FindBin(trkPt);    	   
	    binCont = fHCorrectMCnegChrg->GetBinContent(ptBinMC);
	    if(binCont!=0) ptWgtMC = 1.0/binCont; 
	  }
	  if(fParticle==1 && isItPion && fHCorrectMCnegPion){
	    ptBinMC = fHCorrectMCnegPion->FindBin(trkPt);    
	    binCont = fHCorrectMCnegPion->GetBinContent(ptBinMC);
	    if(binCont!=0)  ptWgtMCPion = 1.0/binCont; 	    
	  }
	  if(fParticle==2 && isItKaon && fHCorrectMCnegKaon){
	    ptBinMC = fHCorrectMCnegKaon->FindBin(trkPt);
	    binCont = fHCorrectMCnegKaon->GetBinContent(ptBinMC);
	    if(binCont!=0)  ptWgtMCKaon = 1.0/binCont;
	  }
	  if(fParticle==3 && isItProt && fHCorrectMCnegProt){
	    ptBinMC     = fHCorrectMCnegProt->FindBin(trkPt);
	    binCont = fHCorrectMCnegProt->GetBinContent(ptBinMC);
	    if(binCont!=0)  ptWgtMCProt = 1.0/binCont; 	    
	  }
	  
	  ////NUA:
	  if(fHCorrectNUAnegChrg){
	    //iBinNUA = fHCorrectNUAnegChrg[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUA  = fHCorrectNUAnegChrg[cForNUA]->GetBinContent(iBinNUA);  
	    iBinNUA = fHCorrectNUAnegChrg->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUA  = fHCorrectNUAnegChrg->GetBinContent(iBinNUA);  
	  }
	  if(fParticle==1 && isItPion && fHCorrectNUAnegPion){
	    //iBinNUA     = fHCorrectNUAnegPion[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUAPion  = fHCorrectNUAnegPion[cForNUA]->GetBinContent(iBinNUA);
	    iBinNUA     = fHCorrectNUAnegPion->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAPion  = fHCorrectNUAnegPion->GetBinContent(iBinNUA);
	  }
	  if(fParticle==2 && isItKaon && fHCorrectNUAnegKaon){
	    //iBinNUA     = fHCorrectNUAnegKaon[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUAKaon  = fHCorrectNUAnegKaon[cForNUA]->GetBinContent(iBinNUA);
	    iBinNUA     = fHCorrectNUAnegKaon->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAKaon  = fHCorrectNUAnegKaon->GetBinContent(iBinNUA);
	  }
	  if(fParticle==3 && isItProt && fHCorrectNUAnegProt){
	    //iBinNUA     = fHCorrectNUAnegProt[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    //WgtNUAProt  = fHCorrectNUAnegProt[cForNUA]->GetBinContent(iBinNUA);
	    iBinNUA     = fHCorrectNUAnegProt->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAProt  = fHCorrectNUAnegProt->GetBinContent(iBinNUA);
	  }	  
	}

	/// infinity Weight protection:
	if(WgtNUA>1e3)     WgtNUA = 1.0;
	if(WgtNUAPion>1e3) WgtNUAPion = 1.0;
	if(WgtNUAKaon>1e3) WgtNUAKaon = 1.0;
	if(WgtNUAProt>1e3) WgtNUAProt = 1.0;

	
	
	trkWgt     = WgtNUA*ptWgtMC;
	trkWgtPion = WgtNUAPion*ptWgtMCPion;
	trkWgtKaon = WgtNUAKaon*ptWgtMCKaon;
	trkWgtProt = WgtNUAProt*ptWgtMCProt;


	///----------------------------------------------------------------	




    
	uqRe = TMath::Cos(gPsiN*trkPhi);
	uqIm = TMath::Sin(gPsiN*trkPhi);
	



	///---------- pT cut for v2 vs Ach -------------
	
      
	for (int j=0;j<noofptbins;j++) //pt loop
	  
	  {

	    if (trkPt>lowpt[j] && trkPt<highpt[j])
	      {
		
		ptcenter[j]=(lowpt[j]+highpt[j])/2.0;
		
	
	if(trkChrg > 0){
	
	  
	  if(trkEta > fEtaGapPos){
	    sumQ2xChrgPosEtaPos[j] += trkWgt*uqRe;
	    sumQ2yChrgPosEtaPos[j] += trkWgt*uqIm;
	    NumOfChrgPosEtaPos[j]  += trkWgt;
	  }
	  else if(trkEta < fEtaGapNeg){
	    sumQ2xChrgPosEtaNeg[j] += trkWgt*uqRe;
	    sumQ2yChrgPosEtaNeg[j] += trkWgt*uqIm;
	    NumOfChrgPosEtaNeg[j]  += trkWgt;
	  }

	 
	  
	  if(fParticle==1 && isItPion){
	
	    if(trkEta > fEtaGapPos){
	      sumQ2xPionPosEtaPos[j] += trkWgtPion*uqRe;
	      sumQ2yPionPosEtaPos[j] += trkWgtPion*uqIm;
	      NumOfPionPosEtaPos[j]  += trkWgtPion;
	    }
	    else if(trkEta < fEtaGapNeg){
	      sumQ2xPionPosEtaNeg[j] += trkWgtPion*uqRe;
	      sumQ2yPionPosEtaNeg[j] += trkWgtPion*uqIm;
	      NumOfPionPosEtaNeg[j]  += trkWgtPion;
	    }
	    //fHFillNUAPosPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	  }
	  if(fParticle==2 && isItKaon){
	    // fHistv2AchKaonPos[0][iCent]->Fill(fAchrgNet, (uqRe*sumQxTPCneg + uqIm*sumQyTPCneg)/sumWgtneg, trkWgtKaon);
	      
	    if(trkEta > fEtaGapPos){
	      sumQ2xKaonPosEtaPos[j] += trkWgtKaon*uqRe;
	      sumQ2yKaonPosEtaPos[j] += trkWgtKaon*uqIm;
	      NumOfKaonPosEtaPos[j]  += trkWgtKaon;
	    }
	    else if(trkEta < fEtaGapNeg){
	      sumQ2xKaonPosEtaNeg[j] += trkWgtKaon*uqRe;
	      sumQ2yKaonPosEtaNeg[j] += trkWgtKaon*uqIm;
	      NumOfKaonPosEtaNeg[j]  += trkWgtKaon;
	    }
	
	  }
	  if(fParticle==3 && isItProt){
	    //fHistv2AchProtPos[0][iCent]->Fill(fAchrgNet, (uqRe*sumQxTPCneg + uqIm*sumQyTPCneg)/sumWgtneg, trkWgtProt);
	      
	    if(trkEta > fEtaGapPos){
	      sumQ2xProtPosEtaPos[j] += trkWgtProt*uqRe;
	      sumQ2yProtPosEtaPos[j] += trkWgtProt*uqIm;
	      NumOfProtPosEtaPos[j]  += trkWgtProt;
	    }
	    else if(trkEta < fEtaGapNeg){
	      sumQ2xProtPosEtaNeg[j] += trkWgtProt*uqRe;
	      sumQ2yProtPosEtaNeg[j] += trkWgtProt*uqIm;
	      NumOfProtPosEtaNeg[j]  += trkWgtProt;
	    }
	
	  }

	  
	  
	}///+ve Ch done	
	else{  //-Ve charge
	  
	  
	  if(trkEta > fEtaGapPos){
	    sumQ2xChrgNegEtaPos[j] += trkWgt*uqRe;
	    sumQ2yChrgNegEtaPos[j] += trkWgt*uqIm;
	    NumOfChrgNegEtaPos[j]  += trkWgt;
	  }
	  else if(trkEta < fEtaGapNeg){
	    sumQ2xChrgNegEtaNeg[j] += trkWgt*uqRe;
	    sumQ2yChrgNegEtaNeg[j] += trkWgt*uqIm;
	    NumOfChrgNegEtaNeg[j]  += trkWgt;
	  }
	   
	  if(fParticle==1 && isItPion){
	    // fHistv2AchPionNeg[0][iCent]->Fill(fAchrgNet, (uqRe*sumQxTPCneg + uqIm*sumQyTPCneg)/sumWgtneg, trkWgtPion);

	    if(trkEta > fEtaGapPos){
	      sumQ2xPionNegEtaPos[j] += trkWgtPion*uqRe;
	      sumQ2yPionNegEtaPos[j] += trkWgtPion*uqIm;
	      NumOfPionNegEtaPos[j]  += trkWgtPion;
	    }
	    else if(trkEta < fEtaGapNeg){
	      sumQ2xPionNegEtaNeg[j] += trkWgtPion*uqRe;
	      sumQ2yPionNegEtaNeg[j] += trkWgtPion*uqIm;
	      NumOfPionNegEtaNeg[j]  += trkWgtPion;
	    }	   
	    //fHFillNUANegPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);	  
	  }
	  if(fParticle==2 && isItKaon){
	    //fHistv2AchKaonNeg[0][iCent]->Fill(fAchrgNet, (uqRe*sumQxTPCneg + uqIm*sumQyTPCneg)/sumWgtneg, trkWgtKaon);

	    if(trkEta > fEtaGapPos){
	      sumQ2xKaonNegEtaPos[j] += trkWgtKaon*uqRe;
	      sumQ2yKaonNegEtaPos[j] += trkWgtKaon*uqIm;
	      NumOfKaonNegEtaPos[j]  += trkWgtKaon;
	    }
	    else if(trkEta < fEtaGapNeg){
	      sumQ2xKaonNegEtaNeg[j] += trkWgtKaon*uqRe;
	      sumQ2yKaonNegEtaNeg[j] += trkWgtKaon*uqIm;
	      NumOfKaonNegEtaNeg[j]  += trkWgtKaon;
	    }	  
	
	  }
	  if(fParticle==3 && isItProt){
	    //fHistv2AchProtNeg[0][iCent]->Fill(fAchrgNet, (uqRe*sumQxTPCneg + uqIm*sumQyTPCneg)/sumWgtneg, trkWgtProt);

	    if(trkEta > fEtaGapPos){
	      sumQ2xProtNegEtaPos[j] += trkWgtProt*uqRe;
	      sumQ2yProtNegEtaPos[j] += trkWgtProt*uqIm;
	      NumOfProtNegEtaPos[j]  += trkWgtProt;
	    }
	    else if(trkEta < fEtaGapNeg){
	      sumQ2xProtNegEtaNeg[j] += trkWgtProt*uqRe;
	      sumQ2yProtNegEtaNeg[j] += trkWgtProt*uqIm;
	      NumOfProtNegEtaNeg[j]  += trkWgtProt;
	    }
	
	  }

	

	}/// if -ve Particle
	//----------- v2 vs Ach filled ---------

	      }
	  }
       
       

      }//with all trackCuts applied      
    }//-------> if FB is validated
    
  }///------> 2nd Track loop Ends here.<--------
 

    
  /// For cumulant method:


  /// Charge All(+-):


  Double_t c2WeightChrg     = fSumWgtEtaPos*fSumWgtEtaNeg;
  Double_t c2WeightChrgNeg     = fSumWgtEtaPosChNeg*fSumWgtEtaNegChNeg;
  if (c2WeightChrg!=0)
    {
  Double_t c2cumulantChrg   = (fSumTPCQn2xPos*fSumTPCQn2xNeg + fSumTPCQn2yPos*fSumTPCQn2yNeg)/c2WeightChrg;
  fHistv2cumAchChrgAll[iCent]->Fill(1, c2cumulantChrg, c2WeightChrg);   /// for denominator
  fHistv2cumCentChrgAll->Fill(centrality, c2cumulantChrg, c2WeightChrg);   /// for denominator
    }

  if (c2WeightChrgNeg!=0)
    {
  Double_t c2cumulantChrgNeg   = (fSumTPCQn2xPosChNeg*fSumTPCQn2xNegChNeg + fSumTPCQn2yPosChNeg*fSumTPCQn2yNegChNeg)/c2WeightChrgNeg;
  fHistv2cumCentChrgAllNeg->Fill(centrality, c2cumulantChrgNeg, c2WeightChrgNeg);   /// for denominator
    }



  if (fSumWgtEtaNeg!=0)
    {
  for (int j=0;j<noofptbins;j++)  //for v2 vs pt filling loop

    {

  
  if(NumOfChrgPosEtaPos[j]>0){
    ///Chrg+:  
    Double_t c2WeightChrgPos   =  NumOfChrgPosEtaPos[j]*fSumWgtEtaNeg;
    Double_t c2cumulantChrgPos =  (sumQ2xChrgPosEtaPos[j]*fSumTPCQn2xNeg + sumQ2yChrgPosEtaPos[j]*fSumTPCQn2yNeg)/c2WeightChrgPos;
    fHistv2AchChrgPos[0][iCent]->Fill(ptcenter[j], c2cumulantChrgPos, c2WeightChrgPos);   /// for denominator
  }
  
  if(NumOfChrgNegEtaPos[j]>0){  
    ///Chrg-:  
    Double_t c2WeightChrgNeg   =  NumOfChrgNegEtaPos[j]*fSumWgtEtaNeg;
    Double_t c2cumulantChrgNeg =  (sumQ2xChrgNegEtaPos[j]*fSumTPCQn2xNeg + sumQ2yChrgNegEtaPos[j]*fSumTPCQn2yNeg)/c2WeightChrgNeg;
    fHistv2AchChrgNeg[0][iCent]->Fill(ptcenter[j], c2cumulantChrgNeg, c2WeightChrgNeg);   /// for denominator
  }
  






  
  if(fParticle==1){
    if(NumOfPionPosEtaPos[j]>0){
      ///Pion+:  
      Double_t c2WeightPionPos   =  NumOfPionPosEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantPionPos =  (sumQ2xPionPosEtaPos[j]*fSumTPCQn2xNeg + sumQ2yPionPosEtaPos[j]*fSumTPCQn2yNeg)/c2WeightPionPos;
      fHistv2AchPionPos[0][iCent]->Fill(ptcenter[j], c2cumulantPionPos, c2WeightPionPos);   /// for denominator
    }
    if(NumOfPionNegEtaPos[j]>0){  
      ///Pion-:  
      Double_t c2WeightPionNeg   =  NumOfPionNegEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantPionNeg =  (sumQ2xPionNegEtaPos[j]*fSumTPCQn2xNeg + sumQ2yPionNegEtaPos[j]*fSumTPCQn2yNeg)/c2WeightPionNeg;
      fHistv2AchPionNeg[0][iCent]->Fill(ptcenter[j], c2cumulantPionNeg, c2WeightPionNeg);   /// for denominator
    }
  }

  if(fParticle==2){
    if(NumOfKaonPosEtaPos[j]>0){
      ///Kaon+:  
      Double_t c2WeightKaonPos   =  NumOfKaonPosEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantKaonPos =  (sumQ2xKaonPosEtaPos[j]*fSumTPCQn2xNeg + sumQ2yKaonPosEtaPos[j]*fSumTPCQn2yNeg)/c2WeightKaonPos;
      fHistv2AchKaonPos[0][iCent]->Fill(ptcenter[j], c2cumulantKaonPos, c2WeightKaonPos);   /// for denominator
    }
    if(NumOfKaonNegEtaPos[j]>0){
      ///Kaon-:  
      Double_t c2WeightKaonNeg   =  NumOfKaonNegEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantKaonNeg =  (sumQ2xKaonNegEtaPos[j]*fSumTPCQn2xNeg + sumQ2yKaonNegEtaPos[j]*fSumTPCQn2yNeg)/c2WeightKaonNeg;
      fHistv2AchKaonNeg[0][iCent]->Fill(ptcenter[j], c2cumulantKaonNeg, c2WeightKaonNeg);   /// for denominator
    }
  }

  if(fParticle==3){
    if(NumOfProtPosEtaPos[j]>0){
      ///Prot+:  
      Double_t c2WeightProtPos   =  NumOfProtPosEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantProtPos =  (sumQ2xProtPosEtaPos[j]*fSumTPCQn2xNeg + sumQ2yProtPosEtaPos[j]*fSumTPCQn2yNeg)/c2WeightProtPos;
      fHistv2AchProtPos[0][iCent]->Fill(ptcenter[j], c2cumulantProtPos, c2WeightProtPos);   /// for denominator
    }

    if(NumOfProtNegEtaPos[j]>0){
      ///Prot-:  
      Double_t c2WeightProtNeg   =  NumOfProtNegEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantProtNeg =  (sumQ2xProtNegEtaPos[j]*fSumTPCQn2xNeg + sumQ2yProtNegEtaPos[j]*fSumTPCQn2yNeg)/c2WeightProtNeg;
      fHistv2AchProtNeg[0][iCent]->Fill(ptcenter[j], c2cumulantProtNeg, c2WeightProtNeg);   /// for denominator
    }
  }
  
    }
  
    }

  //Last lines in Event loop
  fCentDistAfterCut->Fill(centrality);
  fHistEventCount->Fill(14.5); //final Event count.
  PostData(1,fListHist);

  // std::cout<<" Info:UserExec()  Call Finished ..!!!\n";
  

}//---------------- UserExec ----------------------




Bool_t AliAnalysisTaskv2pt::CheckEventIsPileUp(AliAODEvent *faod) {
  
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
  fHistTPCVsESDTrkBefore->Fill(multTPCAll,multEsd);   //A. Dobrin
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
  if(!BisPileup){
    fHistTPCVsESDTrkAfter->Fill(multTPCAll,multEsd);  
    fHistTPConlyVsCL1After->Fill(centrCL1,multTPCAll);
    fHistTPConlyVsV0MAfter->Fill(centrV0M,multTPCAll);
    fHistGlobalVsV0MAfter->Fill(centrV0M, multGlobal);
  
    //fHistEventCount->Fill(stepCount); //6
    //stepCount++
  }

  
  return BisPileup; 

  
} //-------pile up function ------




Bool_t AliAnalysisTaskv2pt::PileUpMultiVertex(const AliAODEvent* faod)
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




double AliAnalysisTaskv2pt::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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












Bool_t AliAnalysisTaskv2pt::CheckEventIsPileUp2018(AliAODEvent *faod) {


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
    {
      //cout<<"*****************hi i am in 6th**************************:"<<endl;
      BisPileup=kTRUE;
    }
  

  if (faod->IsIncompleteDAQ())
    {
      //cout<<"*****************hi i am in 7th**************************:"<<endl;
      BisPileup=kTRUE;
    }
        
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
















void AliAnalysisTaskv2pt::SetupEventAndTaskConfigInfo(){
 
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
  fHistEventCount->GetXaxis()->SetBinLabel(3,"AOD Exist");
  fHistEventCount->GetXaxis()->SetBinLabel(4,"Vz < 10");
  fHistEventCount->GetXaxis()->SetBinLabel(5,Form("%2.0f<Cent<%2.0f",fCentralityMin,fCentralityMax));
  fHistEventCount->GetXaxis()->SetBinLabel(6,"noAODtrack > 2 ");
  fHistEventCount->GetXaxis()->SetBinLabel(7,"TPC vs Global");
  fHistEventCount->GetXaxis()->SetBinLabel(8,"TPC128 vs ESD");
  fHistEventCount->GetXaxis()->SetBinLabel(9,"Cent vs TPC");
  fHistEventCount->GetXaxis()->SetBinLabel(10,"mult eta+/- > 2");
  fHistEventCount->GetXaxis()->SetBinLabel(11,"centCL1 < 90");
  fHistEventCount->GetXaxis()->SetBinLabel(15,"Survived Events");
  fListHist->Add(fHistEventCount);
  //fHistEventCount->Fill(1);
}//----------SetupEventAndTaskConfigInfo-----------




Int_t AliAnalysisTaskv2pt::GetCentralityScaled0to10(Double_t fCent){

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

void AliAnalysisTaskv2pt::GetNUACorrectionHist(Int_t run, Int_t kParticleID)
{
  /*
  if(fListNUACorr){

    if(kParticleID==0){ //charge
      for(int i=0;i<5;i++){
	fHCorrectNUAposChrg[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Pos_Cent%d_Run%d",i,run)); 
	fHCorrectNUAnegChrg[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Neg_Cent%d_Run%d",i,run));
      }
    }
    else if(kParticleID==1){ //Pion
      for(int i=0;i<5;i++){
	fHCorrectNUAposPion[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Pos_Cent%d_Run%d",i,run)); 
	fHCorrectNUAnegPion[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Neg_Cent%d_Run%d",i,run));
      }
    }
    else if(kParticleID==2){ //Kaon
      for(int i=0;i<5;i++){
	fHCorrectNUAposKaon[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Pos_Cent%d_Run%d",i,run)); 
	fHCorrectNUAnegKaon[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Neg_Cent%d_Run%d",i,run));
      }
    }
    else if(kParticleID==3){ //Proton
      for(int i=0;i<5;i++){
	fHCorrectNUAposProt[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Pos_Cent%d_Run%d",i,run)); 
	fHCorrectNUAnegProt[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Neg_Cent%d_Run%d",i,run));
      }
    }
    else{
      for(int i=0;i<5;i++){
	fHCorrectNUAposChrg[i]=NULL;
	fHCorrectNUAnegChrg[i]=NULL;
	fHCorrectNUAposPion[i]=NULL;
	fHCorrectNUAnegPion[i]=NULL;
	fHCorrectNUAposKaon[i]=NULL;
	fHCorrectNUAnegKaon[i]=NULL;
	fHCorrectNUAposProt[i]=NULL;
	fHCorrectNUAnegProt[i]=NULL;
      }
    }
    
  }//------> if list Exist

  // else {
  //   printf("\n ******** Warning: No NUA Correction File/List...!! \n Run= %d, Use NUA Wgt = 1.0 ********* \n",run);
  // }
  */



  if(fListNUACorr){

    if(kParticleID==0){ //charge
      //for(int i=0;i<5;i++){
      //fHCorrectNUAposChrg[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Pos_Cent%d_Run%d",i,run)); 
      //fHCorrectNUAnegChrg[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Neg_Cent%d_Run%d",i,run));

      fHCorrectNUAposChrg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Pos_Cent%d_Run%d",0,run)); 
      fHCorrectNUAnegChrg = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Charge_Neg_Cent%d_Run%d",0,run));
      //}
    }
    else if(kParticleID==1){ //Pion
      //for(int i=0;i<5;i++){
      //fHCorrectNUAposPion[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Pos_Cent%d_Run%d",i,run)); 
      //fHCorrectNUAnegPion[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Neg_Cent%d_Run%d",i,run));
      fHCorrectNUAposPion = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Pos_Cent%d_Run%d",0,run)); 
      fHCorrectNUAnegPion = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pion_Neg_Cent%d_Run%d",0,run));
      //}
    }
    else if(kParticleID==2){ //Kaon
      //for(int i=0;i<5;i++){
      //fHCorrectNUAposKaon[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Pos_Cent%d_Run%d",i,run)); 
      //fHCorrectNUAnegKaon[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Neg_Cent%d_Run%d",i,run));

      fHCorrectNUAposKaon = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Pos_Cent%d_Run%d",0,run)); 
      fHCorrectNUAnegKaon = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Kaon_Neg_Cent%d_Run%d",0,run));
      //}
    }
    else if(kParticleID==3){ //Proton
      //for(int i=0;i<5;i++){
      //fHCorrectNUAposProt[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Pos_Cent%d_Run%d",i,run)); 
      //fHCorrectNUAnegProt[i] = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Neg_Cent%d_Run%d",i,run));

      fHCorrectNUAposProt = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Pos_Cent%d_Run%d",0,run)); 
      fHCorrectNUAnegProt = (TH3F *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Proton_Neg_Cent%d_Run%d",0,run));
      //}
    }
    else{
      //for(int i=0;i<5;i++){
      /*fHCorrectNUAposChrg[i]=NULL;
	fHCorrectNUAnegChrg[i]=NULL;
	fHCorrectNUAposPion[i]=NULL;
	fHCorrectNUAnegPion[i]=NULL;
	fHCorrectNUAposKaon[i]=NULL;
	fHCorrectNUAnegKaon[i]=NULL;
	fHCorrectNUAposProt[i]=NULL;
	fHCorrectNUAnegProt[i]=NULL;
      */
      fHCorrectNUAposChrg=NULL;
      fHCorrectNUAnegChrg=NULL;
      fHCorrectNUAposPion=NULL;
      fHCorrectNUAnegPion=NULL;
      fHCorrectNUAposKaon=NULL;
      fHCorrectNUAnegKaon=NULL;
      fHCorrectNUAposProt=NULL;
      fHCorrectNUAnegProt=NULL;
      //}
    }
    
  }//------> if list Exist

  // else {
  //   printf("\n ******** Warning: No NUA Correction File/List...!! \n Run= %d, Use NUA Wgt = 1.0 ********* \n",run);
  // }
}






////---------- SetUp Tracking Efficiency Correction Map ---------------

void AliAnalysisTaskv2pt::GetMCCorrectionHist(Int_t run){

  if(fListTRKCorr) {
    //cout<<"\n =========> Info: Found TList with MC Tracking Corr Histograms <=========== "<<endl;
    fHCorrectMCposChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgPos");
    fHCorrectMCposPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionPos");
    fHCorrectMCposKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonPos");
    fHCorrectMCposProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtPos");

    fHCorrectMCnegChrg =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyChrgNeg");
    fHCorrectMCnegPion =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyPionNeg");
    fHCorrectMCnegKaon =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyKaonNeg");
    fHCorrectMCnegProt =  (TH1D *) fListTRKCorr->FindObject("trkEfficiencyProtNeg");

    //for(int i=0;i<10;i++) {
    //fFB_Efficiency_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_%d",i));
    //}
  }
  // else if(!fListTRKCorr){
  //   std::cout<<"\n\n !!!!**** Warning : No FB Efficiency Correction, run = "<<run<<"..!!!**** \n using MC TrkWgt = 1.0 \n"<<std::endl;
  // }
}
