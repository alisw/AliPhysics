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
#include "AliAnalysisTaskCVEUNIpt.h"

using std::cout;
using std::endl;
using std::vector;


ClassImp(AliAnalysisTaskCVEUNIpt)

AliAnalysisTaskCVEUNIpt::AliAnalysisTaskCVEUNIpt(const char *name): AliAnalysisTaskSE(name),
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
  sCentrEstimator("V0M"),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fHistEtaPtBeforCut(NULL),
  fHistEtaPhiBeforCut(NULL),
  fHistEtaPhiAfterCut(NULL),  
  fHCorrectMCposChrg(NULL),
  fHCorrectMCposPion(NULL),
  fHCorrectMCposKaon(NULL),
  fHCorrectMCposProt(NULL),
  fHCorrectMCnegChrg(NULL),
  fHCorrectMCnegPion(NULL),
  fHCorrectMCnegKaon(NULL),
  fHCorrectMCnegProt(NULL),      
  fHistAChrgVsCent(NULL),
  
  fHistEventCount(NULL)
{

  for(int i=0;i<2;i++){
    for(int j=0;j<10;j++){
  for (int k=0;k<10;k++){
      fHistv2AchChrgPos[i][j][k] = NULL;
      fHistv2AchKaonPos[i][j][k] = NULL;    
      fHistv2AchProtPos[i][j][k] = NULL;      
      fHistv2AchPionPos[i][j][k] = NULL;
      
      fHistv2AchChrgNeg[i][j][k] = NULL;
      fHistv2AchPionNeg[i][j][k] = NULL;      
      fHistv2AchKaonNeg[i][j][k] = NULL;
      fHistv2AchProtNeg[i][j][k] = NULL;      
    }
  }
  }
  
  for(int i=0; i<5; i++){
    fHCorrectNUAposChrg[i] = NULL;  
    fHCorrectNUAnegChrg[i] = NULL;  
    fHCorrectNUAposPion[i] = NULL;  
    fHCorrectNUAnegPion[i] = NULL;  
    fHCorrectNUAposKaon[i] = NULL;  
    fHCorrectNUAnegKaon[i] = NULL;  
    fHCorrectNUAposProt[i] = NULL;  
    fHCorrectNUAnegProt[i] = NULL;    
  }

  for(int i=0; i<5; i++){
    fHFillNUAPosPID[i]  = NULL;
    fHFillNUANegPID[i]  = NULL; 
  }

  for(int i=0; i<2; i++){
    fHistEPResolution[i] = NULL;
  }

  for(int i=0; i<10; i++){
    fHistEPResolutionAch[i] = NULL;
  }

   for(int i=0; i<10; i++){
     for(int k=0; k<10; k++){ 
   fHistv2cumAchChrgAll[i][k] = NULL;
  }
  }
  
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskCVEUNIpt::AliAnalysisTaskCVEUNIpt():
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
  sCentrEstimator("V0M"),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fHistEtaPtBeforCut(NULL),
  fHistEtaPhiBeforCut(NULL),
  fHistEtaPhiAfterCut(NULL),  
  fHCorrectMCposChrg(NULL),
  fHCorrectMCposPion(NULL),
  fHCorrectMCposKaon(NULL),
  fHCorrectMCposProt(NULL),
  fHCorrectMCnegChrg(NULL),
  fHCorrectMCnegPion(NULL),
  fHCorrectMCnegKaon(NULL),
  fHCorrectMCnegProt(NULL),    
  fHistAChrgVsCent(NULL),
  
  fHistEventCount(NULL)
{

  for(int i=0;i<2;i++){
    for(int j=0;j<10;j++){
  for(int k=0;k<10;k++){
      fHistv2AchChrgPos[i][j][k] = NULL;
      fHistv2AchKaonPos[i][j][k] = NULL;    
      fHistv2AchProtPos[i][j][k] = NULL;      
      fHistv2AchPionPos[i][j][k] = NULL;
      
      fHistv2AchChrgNeg[i][j][k] = NULL;
      fHistv2AchPionNeg[i][j][k] = NULL;      
      fHistv2AchKaonNeg[i][j][k] = NULL;
      fHistv2AchProtNeg[i][j][k] = NULL;      
    }
  }
  }

  for(int i=0; i<5; i++){
    fHCorrectNUAposChrg[i] = NULL;  
    fHCorrectNUAnegChrg[i] = NULL;  
    fHCorrectNUAposPion[i] = NULL;  
    fHCorrectNUAnegPion[i] = NULL;  
    fHCorrectNUAposKaon[i] = NULL;  
    fHCorrectNUAnegKaon[i] = NULL;  
    fHCorrectNUAposProt[i] = NULL;  
    fHCorrectNUAnegProt[i] = NULL;  
  }

  for(int i=0; i<5; i++){
    fHFillNUAPosPID[i]  = NULL;
    fHFillNUANegPID[i]  = NULL; 
  }

  for(int i=0; i<2; i++){
    fHistEPResolution[i] = NULL;
  }

  for(int i=0; i<10; i++){
    fHistEPResolutionAch[i] = NULL;
  }
    
  for(int i=0; i<10; i++){
    for(int k=0; k<10; k++){  
  fHistv2cumAchChrgAll[i][k] = NULL;
  }
  }
  
  //Not needed for Empty Constructor:
  //DefineInput(0,TChain::Class());
  //DefineOutput(1,TList::Class());
}
  
//__________________ destructor ___________________
AliAnalysisTaskCVEUNIpt::~AliAnalysisTaskCVEUNIpt()
{
  if(fListHist)      delete fListHist;  
  if(fAnalysisUtil)  delete fAnalysisUtil;   // because its 'new' !!

  //Delete the clones
  if(fListTRKCorr)  delete fListTRKCorr;
  if(fListNUACorr)  delete fListNUACorr;
  if(fListV0MCorr)  delete fListV0MCorr;
      
}










//________________ Define Histograms _______________
void AliAnalysisTaskCVEUNIpt::UserCreateOutputObjects()
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
  
  fHistEtaPtBeforCut = new TH2F("fHistEtaPtBeforCut","#eta vs p_{T} (wFB, w/o  cut)",24,-1.2,1.2,50,0,5);
  fListHist->Add(fHistEtaPtBeforCut);
  fHistEtaPhiBeforCut = new TH2F("fHistPhiEtaBeforCut","#phi vs #eta (wFB, w/o cut)",50,0,6.2835,24,-1.2,1.2);
  fListHist->Add(fHistEtaPhiBeforCut);  
  fHistEtaPhiAfterCut = new TH2F("fHistPhiEtaAfterCut","#phi vs #eta (with Wgts)",50,0,6.2835,24,-1.2,1.2);
  fListHist->Add(fHistEtaPhiAfterCut);
  

  
  //----------- User's histograms: --------------
  Char_t  name[100];
  Char_t title[100];
 

  Double_t centRange[11] = {0,5,10,20,30,40,50,60,70,80,90};

  fHistAChrgVsCent = new TH2F("fHistAChrgVsCent","Ach vs Cent;Cent;Ach",10,centRange,500,-1.0,1.0);
  fListHist->Add(fHistAChrgVsCent);

  const int NPT=21;
  //Double_t PT[NPT+1]={0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,3.0,4.0,5.0};
  //Double_t PT[NPT+1]={0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0};
  Double_t PT[NPT+1]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.2,1.4,1.6,1.8,2.0};

  
  // v2 vs Ach


  for(int i=0;i<2;i++){
    for(int j=0;j<10;j++){
  for (int k=0;k<10;k++){
      ////Charge:
      sprintf(name,"fHistv2AchChrgPos_Method%d_Cent%d_Ach%d",i,j,k);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchChrgPos[i][j][k] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchChrgPos[i][j][k]->Sumw2();
      fListHist->Add(fHistv2AchChrgPos[i][j][k]);
      sprintf(name,"fHistv2AchChrgNeg_Method%d_Cent%d_Ach%d",i,j,k);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchChrgNeg[i][j][k] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchChrgNeg[i][j][k]->Sumw2();
      fListHist->Add(fHistv2AchChrgNeg[i][j][k]);      

      //// Pion:
      sprintf(name,"fHistv2AchPionPos_Method%d_Cent%d_Ach%d",i,j,k);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchPionPos[i][j][k] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchPionPos[i][j][k]->Sumw2();
      fListHist->Add(fHistv2AchPionPos[i][j][k]);
      sprintf(name,"fHistv2AchPionNeg_Method%d_Cent%d_Ach%d",i,j,k);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchPionNeg[i][j][k] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchPionNeg[i][j][k]->Sumw2();
      fListHist->Add(fHistv2AchPionNeg[i][j][k]);      
 
      //// Kaon:
      sprintf(name,"fHistv2AchKaonPos_Method%d_Cent%d_Ach%d",i,j,k);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchKaonPos[i][j][k] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchKaonPos[i][j][k]->Sumw2();
      fListHist->Add(fHistv2AchKaonPos[i][j][k]);
      sprintf(name,"fHistv2AchKaonNeg_Method%d_Cent%d_Ach%d",i,j,k);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchKaonNeg[i][j][k] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchKaonNeg[i][j][k]->Sumw2();
      fListHist->Add(fHistv2AchKaonNeg[i][j][k]);      

      //// Proton:
      sprintf(name,"fHistv2AchProtPos_Method%d_Cent%d_Ach%d",i,j,k);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchProtPos[i][j][k] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchProtPos[i][j][k]->Sumw2();
      fListHist->Add(fHistv2AchProtPos[i][j][k]);
      sprintf(name,"fHistv2AchProtNeg_Method%d_Cent%d_Ach%d",i,j,k);
      sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[j],centRange[j+1]);
      fHistv2AchProtNeg[i][j][k] = new TProfile(name,title,NPT,PT,"");
      fHistv2AchProtNeg[i][j][k]->Sumw2();
      fListHist->Add(fHistv2AchProtNeg[i][j][k]);      
    }
  }
  }


  Int_t gCentForNUA[6] = {0,5,10,20,40,90};
  Char_t cpid[10];

 
  if(fParticle==1) sprintf(cpid,"Pion,Id %d",fParticle);
  else if(fParticle==2) sprintf(cpid,"Kaon,Id %d",fParticle);
  else if(fParticle==3) sprintf(cpid,"Prot,Id %d",fParticle);
  else  sprintf(cpid,"Charge,Id %d",fParticle);

 
  for(int i=0; i<5; i++){  //centrality bins
    sprintf(name,"fHistEtaPhiVz_%d_Pos_Cent%d_Run%d",fParticle,i,1); 
    sprintf(title,"%s Pos, Cent%d-%d, FB %d",cpid,gCentForNUA[i],gCentForNUA[i+1],fFilterBit);
    fHFillNUAPosPID[i] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
    fListHist->Add(fHFillNUAPosPID[i]);

    sprintf(name,"fHistEtaPhiVz_%d_Neg_Cent%d_Run%d",fParticle,i,1); 
    sprintf(title,"%s Neg, Cent%d-%d, FB %d",cpid,gCentForNUA[i],gCentForNUA[i+1],fFilterBit);
    fHFillNUANegPID[i] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
    fListHist->Add(fHFillNUANegPID[i]);    
  }



  
  for(int i=0; i<10; i++){
    ////Charge:
    sprintf(name,"fHistResolutionvsAch_Cent%d",i);
    sprintf(title,"Cent %2.0f-%2.0f; A_{ch}; Resolution",centRange[i],centRange[i+1]);
    fHistEPResolutionAch[i] = new TProfile(name,title,10,-0.1,0.1,"");
    fHistEPResolutionAch[i]->Sumw2();
    fListHist->Add(fHistEPResolutionAch[i]);
  }

    for(int i=0; i<10; i++){
      for(int k=0; k<10; k++){
    ////Charge:
    sprintf(name,"fHistv2cumAchChrgAllQcumCent%d_Ach%d",i,k);
    sprintf(title,"Cent %2.0f-%2.0f; p_{T}; v_{2}",centRange[i],centRange[i+1]);
    fHistv2cumAchChrgAll[i][k] = new TProfile(name,title,1,-10,10,"");
    fHistv2cumAchChrgAll[i][k]->Sumw2();
    fListHist->Add(fHistv2cumAchChrgAll[i][k]);
  }
  }

  
  for(int i=0; i<2; i++){
    sprintf(name,"fHistEPResolution_Method%d",i);
    sprintf(title,"EP Resolution; centrality; Reso");
    fHistEPResolution[i] = new TProfile(name,title,10,centRange,"");
    fHistEPResolution[i]->Sumw2();
    fListHist->Add(fHistEPResolution[i]);
  }



  



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
void AliAnalysisTaskCVEUNIpt::UserExec(Option_t*) {
 
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
  Double_t fSumWgtEtaNeg=0, fSumWgtEtaPos=0;
  Double_t fNumOfPos = 0;
  Double_t fNumOfNeg = 0;
  Double_t ptWgtMC = 1.0, WgtNUA = 1.0;
  Int_t ptBinMC = 1, iBinNUA = 1;
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

      fHistEtaPtBeforCut->Fill(trkEta, trkPt);
      fHistEtaPhiBeforCut->Fill(trkPhi,trkEta);
      
      /// This Next function is called After Filter bit is validated!! (Otherwise code breaks!)
      trkdEdx  = AODtrack->GetDetPid()->GetTPCsignal();  

      //Apply track cuts here:
      //if((trkPt <= fMaxPtCut) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= 4.0) && TMath::Abs(trkChrg)) {

      if((trkPt <= 10.0) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= 4.0) && TMath::Abs(trkChrg)) {

	//dcaXY  = track->DCA();
	//dcaZ   = track->ZAtDCA();
        
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
	  
	  if(fHCorrectNUAposChrg[cForNUA]){
	    iBinNUA = fHCorrectNUAposChrg[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUA  = fHCorrectNUAposChrg[cForNUA]->GetBinContent(iBinNUA);	    
	  }
	}
	else{
	  
	  if(fHCorrectMCnegChrg){
	    ptBinMC = fHCorrectMCnegChrg->FindBin(trkPt);    //Charge independent MC correction atm.	    
	    binCont = fHCorrectMCnegChrg->GetBinContent(ptBinMC);
	    if(binCont!=0) ptWgtMC = 1.0/binCont;  
	  }
	  
	  if(fHCorrectNUAnegChrg[cForNUA]){
	    iBinNUA = fHCorrectNUAnegChrg[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUA  = fHCorrectNUAnegChrg[cForNUA]->GetBinContent(iBinNUA);  
	  }
	}

	//RefMultRawFB++;
	//RefMultCorrFB += ptWgtMC;
	if(WgtNUA>1e3) WgtNUA = 1;   // For empty bins which cause Wgt = infinity !!

	trkWgt = WgtNUA*ptWgtMC;


	if(trkChrg > 0){	  
	  fNumOfPos += trkWgt;
	}
	else{
	  fNumOfNeg += trkWgt;
	}

	
	//if(iTrack%10==0){
	//std::cout<<" pT = "<<trkPt<<"\t MCWgt = "<<ptWgtMC<<"\t Eta = "<<trkEta<<"\t NUAwgt = "<<WgtNUA<<"\t TotalWgt = "<<trkWgt<<endl;
	//}
       

	if (trkPt<fMaxevpt)
	  {

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
	


	
	
	//<---------- User track analysis Done---------------

      }// when all trackCuts applied      
    }//-------> when FB is validated




    
    ///------- For Pile-UP removal Purpose only-----
    if(AODtrack->TestFilterBit(128))      fMultTPCFull++; // A. Dobrin TPC vs ESD PileUp Cut.
    if(!AODtrack->TestFilterBit(16) || AODtrack->Chi2perNDF() < 0.1) continue;
    Double_t    bval[2] = {-99., -99.};
    Double_t    bCov[3] = {-99., -99., -99.};
    AliAODTrack copy(*AODtrack);
    if(copy.PropagateToDCA(fVevent->GetPrimaryVertex(), fVevent->GetMagneticField(), 100., bval, bCov) && TMath::Abs(bval[0]) < 0.3 && TMath::Abs(bval[1]) < 0.3){
      fMultGlobal++;
    }///MultGlobal Condition
    
  }///------> 1st Track loop Ends here.<--------




  if(iCent > 9) return;

  fHistEventCount->Fill(stepCount); //5
  stepCount++;


  if(fSumWgtEtaNeg <= 0 || fSumWgtEtaPos <= 0) return;
  
  fHistEventCount->Fill(stepCount); //6
  stepCount++;
  

  //// Decide If Pile-Up cut to be applied.
  //fHistEventCount->Fill(stepCount); //7
  //stepCount++;

  
  
  Float_t fAchrgNet = (fNumOfPos - fNumOfNeg)/(fNumOfPos + fNumOfNeg); // Efficiency & NUA Corrected!
  /*
  fHistAChrgVsCent->Fill(centrality, fAchrgNet, fWgtEvent);
  
  // Do Event plane defination, Fill Resolution Etc.
  // Double_t Psi2EtaPos =  Formula for PsiN:
  
  Double_t ResolutionEP = fSumTPCQn2xNeg*fSumTPCQn2xPos + fSumTPCQn2yNeg*fSumTPCQn2yPos;
  ResolutionEP  = ResolutionEP/(fSumWgtEtaPos*fSumWgtEtaNeg);
  Double_t ResWgt = TMath::Sqrt(fSumWgtEtaPos*fSumWgtEtaNeg);

  /// Resolution vs cent:
  fHistEPResolution[0]->Fill(centrality,ResolutionEP,ResWgt*fWgtEvent);

  /// Resolution vs Ach for each cent:
  fHistEPResolutionAch[iCent]->Fill(fAchrgNet,ResolutionEP,ResWgt*fWgtEvent);
  */


  
  // if (TMath::Abs(fAchrgNet)>0.1) return;

 Int_t Ach=0;
 
  if (fAchrgNet>-0.1 && fAchrgNet<-0.08)
    Ach=0;
  else if (fAchrgNet>-0.08 && fAchrgNet<-0.06)
   Ach=1;
  else if (fAchrgNet>-0.06 && fAchrgNet<-0.04)
  Ach=2;
  else if (fAchrgNet>-0.04 && fAchrgNet<-0.02)
   Ach=3;
  else if (fAchrgNet>-0.02 && fAchrgNet<0.0)
   Ach=4;
  else if (fAchrgNet>0.0 && fAchrgNet<0.02)
   Ach=5;
  else if (fAchrgNet>0.02 && fAchrgNet<0.04)
   Ach=6;
  else if (fAchrgNet>0.04 && fAchrgNet<0.06)
   Ach=7;
  else if (fAchrgNet>0.06 && fAchrgNet<0.08)
   Ach=8;
  else if (fAchrgNet>0.08 && fAchrgNet<0.1)
   Ach=9;





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



  // Double_t lowpt[noofptbins]={0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0,3,4};
  //Double_t highpt[noofptbins]={0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,3,4,5};
  Double_t ptcenter[noofptbins]={0.0};

  // Double_t lowpt[noofptbins]={0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0};
  // Double_t highpt[noofptbins]={0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0};

  // Double_t PT[NPT+1]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.2,1.4,1.6,1.8,2.0};
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

      fHistEtaPtBeforCut->Fill(trkEta, trkPt);
      fHistEtaPhiBeforCut->Fill(trkPhi,trkEta);
      
      /// This Next function is called After Filter bit is validated!! (Otherwise code breaks!)
      trkdEdx  = AODtrack->GetDetPid()->GetTPCsignal();  

      //Apply track cuts here:
      if((trkPt <= fMaxPtCut) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= 4.0) && TMath::Abs(trkChrg)) {

	//dcaXY  = track->DCA();
	//dcaZ   = track->ZAtDCA();
        
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


	  if(trkPt<=0.6 && TMath::Abs(nSigTPCpion)<=fNSigmaTPCCut){
	    isItPion = kTRUE;
	  }
	  else if(trkPt>0.6 && trkPt<=10.0 && TMath::Abs(nSigRMSTPCTOFpion)<=fNSigmaTOFCut){
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
	  if(fHCorrectNUAposChrg[cForNUA]){
	    iBinNUA = fHCorrectNUAposChrg[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUA  = fHCorrectNUAposChrg[cForNUA]->GetBinContent(iBinNUA);
	  }
	  if(fParticle==1 && isItPion && fHCorrectNUAposPion[cForNUA]){
	    iBinNUA     = fHCorrectNUAposPion[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAPion  = fHCorrectNUAposPion[cForNUA]->GetBinContent(iBinNUA);
	  }
	  if(fParticle==2 && isItKaon && fHCorrectNUAposKaon[cForNUA]){
	    iBinNUA     = fHCorrectNUAposKaon[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAKaon  = fHCorrectNUAposKaon[cForNUA]->GetBinContent(iBinNUA);
	  }
	  if(fParticle==3 && isItProt && fHCorrectNUAposProt[cForNUA]){
	    iBinNUA     = fHCorrectNUAposProt[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAProt  = fHCorrectNUAposProt[cForNUA]->GetBinContent(iBinNUA);
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
	  if(fHCorrectNUAnegChrg[cForNUA]){
	    iBinNUA = fHCorrectNUAnegChrg[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUA  = fHCorrectNUAnegChrg[cForNUA]->GetBinContent(iBinNUA);  
	  }
	  if(fParticle==1 && isItPion && fHCorrectNUAnegPion[cForNUA]){
	    iBinNUA     = fHCorrectNUAnegPion[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAPion  = fHCorrectNUAnegPion[cForNUA]->GetBinContent(iBinNUA);
	  }
	  if(fParticle==2 && isItKaon && fHCorrectNUAnegKaon[cForNUA]){
	    iBinNUA     = fHCorrectNUAnegKaon[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAKaon  = fHCorrectNUAnegKaon[cForNUA]->GetBinContent(iBinNUA);
	  }
	  if(fParticle==3 && isItProt && fHCorrectNUAnegProt[cForNUA]){
	    iBinNUA     = fHCorrectNUAnegProt[cForNUA]->FindBin(pVtxZ,trkPhi,trkEta);
	    WgtNUAProt  = fHCorrectNUAnegProt[cForNUA]->GetBinContent(iBinNUA);
	  }	  
	}

	/// infinity Weight protection:
	if(WgtNUA>1e3)     WgtNUA = 1.0;
	if(WgtNUAPion>1e3) WgtNUAPion = 1.0;
	if(WgtNUAKaon>1e3) WgtNUAPion = 1.0;
	if(WgtNUAProt>1e3) WgtNUAPion = 1.0;

	
	
	trkWgt     = WgtNUA*ptWgtMC;
	trkWgtPion = WgtNUAPion*ptWgtMCPion;
	trkWgtKaon = WgtNUAKaon*ptWgtMCKaon;
	trkWgtProt = WgtNUAProt*ptWgtMCProt;


	///----------------------------------------------------------------	




    
	uqRe = TMath::Cos(gPsiN*trkPhi);
	uqIm = TMath::Sin(gPsiN*trkPhi);
	
	//----> Remove Track from EP calculation ------
	sumQxTPCneg = fSumTPCQn2xNeg;   
	sumQyTPCneg = fSumTPCQn2yNeg;
	sumQxTPCpos = fSumTPCQn2xPos;
	sumQyTPCpos = fSumTPCQn2yPos;
    
	sumWgtneg = fSumWgtEtaNeg;
	sumWgtpos = fSumWgtEtaPos;

	//// remove AutoCorrelation:
	if(trkEta < fEtaGapNeg){
	  sumQxTPCneg -= trkWgt*uqRe;
	  sumQyTPCneg -= trkWgt*uqIm; 
	  sumWgtneg   -= trkWgt;
	}
	else if(trkEta > fEtaGapPos){
	  sumQxTPCpos -= trkWgt*uqRe;
	  sumQyTPCpos -= trkWgt*uqIm; 
	  sumWgtpos   -= trkWgt;
	}
	//---------------------------------------------

	
	fHistEtaPhiAfterCut->Fill(trkPhi,trkEta,trkWgt);




	///---------- pT cut for v2 vs Ach -------------
	
	// if(trkPt > 2.0) continue;

      
	for (int j=0;j<noofptbins;j++) //pt loop
	  
	  {

	    if (trkPt>lowpt[j] && trkPt<highpt[j])
	      {
		
		ptcenter[j]=(lowpt[j]+highpt[j])/2.0;
		
	
	if(trkChrg > 0){
	
	  //fHistv2AchChrgPos[0][iCent]->Fill(fAchrgNet,  (uqRe*sumQxTPCneg + uqIm*sumQyTPCneg)/sumWgtneg, trkWgt); //This Trk weigth is for uQ

	  // if(fParticle==0)
	  //fHFillNUAPosPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	  
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
	    //fHistv2AchPionPos[0][iCent]->Fill(fAchrgNet, (uqRe*sumQxTPCneg + uqIm*sumQyTPCneg)/sumWgtneg, trkWgtPion);

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
	    fHFillNUAPosPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
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
	    fHFillNUAPosPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	  }

	  
	  
	}///+ve Ch done	
	else{  //-Ve charge
	  
	  //fHistv2AchChrgNeg[0][iCent]->Fill(fAchrgNet,   (uqRe*sumQxTPCneg + uqIm*sumQyTPCneg)/sumWgtneg, trkWgt);

	  //if(fParticle==0)
	  //fHFillNUANegPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);
	  
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
	    fHFillNUANegPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);	  
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
	    fHFillNUANegPID[cForNUA]->Fill(pVtxZ,trkPhi,trkEta);	  	    
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
  Double_t c2cumulantChrg   = (fSumTPCQn2xPos*fSumTPCQn2xNeg + fSumTPCQn2yPos*fSumTPCQn2yNeg)/c2WeightChrg;
  fHistv2cumAchChrgAll[iCent][Ach]->Fill(1, c2cumulantChrg, c2WeightChrg);   /// for denominator

  for (int j=0;j<noofptbins;j++)  //for v2 vs pt filling loop

    {

  
  if(NumOfChrgPosEtaPos[j]>0){
    ///Chrg+:  
    Double_t c2WeightChrgPos   =  NumOfChrgPosEtaPos[j]*fSumWgtEtaNeg;
    Double_t c2cumulantChrgPos =  (sumQ2xChrgPosEtaPos[j]*fSumTPCQn2xNeg + sumQ2yChrgPosEtaPos[j]*fSumTPCQn2yNeg)/c2WeightChrgPos;
    fHistv2AchChrgPos[1][iCent][Ach]->Fill(ptcenter[j], c2cumulantChrgPos, c2WeightChrgPos);   /// for denominator
  }
  
  if(NumOfChrgNegEtaPos[j]>0){  
    ///Chrg-:  
    Double_t c2WeightChrgNeg   =  NumOfChrgNegEtaPos[j]*fSumWgtEtaNeg;
    Double_t c2cumulantChrgNeg =  (sumQ2xChrgNegEtaPos[j]*fSumTPCQn2xNeg + sumQ2yChrgNegEtaPos[j]*fSumTPCQn2yNeg)/c2WeightChrgNeg;
    fHistv2AchChrgNeg[1][iCent][Ach]->Fill(ptcenter[j], c2cumulantChrgNeg, c2WeightChrgNeg);   /// for denominator
  }
  






  
  if(fParticle==1){
    if(NumOfPionPosEtaPos[j]>0){
      ///Pion+:  
      Double_t c2WeightPionPos   =  NumOfPionPosEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantPionPos =  (sumQ2xPionPosEtaPos[j]*fSumTPCQn2xNeg + sumQ2yPionPosEtaPos[j]*fSumTPCQn2yNeg)/c2WeightPionPos;
      fHistv2AchPionPos[1][iCent][Ach]->Fill(ptcenter[j], c2cumulantPionPos, c2WeightPionPos);   /// for denominator
    }
    if(NumOfPionNegEtaPos[j]>0){  
      ///Pion-:  
      Double_t c2WeightPionNeg   =  NumOfPionNegEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantPionNeg =  (sumQ2xPionNegEtaPos[j]*fSumTPCQn2xNeg + sumQ2yPionNegEtaPos[j]*fSumTPCQn2yNeg)/c2WeightPionNeg;
      fHistv2AchPionNeg[1][iCent][Ach]->Fill(ptcenter[j], c2cumulantPionNeg, c2WeightPionNeg);   /// for denominator
    }
  }

  if(fParticle==2){
    if(NumOfKaonPosEtaPos[j]>0){
      ///Kaon+:  
      Double_t c2WeightKaonPos   =  NumOfKaonPosEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantKaonPos =  (sumQ2xKaonPosEtaPos[j]*fSumTPCQn2xNeg + sumQ2yKaonPosEtaPos[j]*fSumTPCQn2yNeg)/c2WeightKaonPos;
      fHistv2AchKaonPos[1][iCent][Ach]->Fill(ptcenter[j], c2cumulantKaonPos, c2WeightKaonPos);   /// for denominator
    }
    if(NumOfKaonNegEtaPos[j]>0){
      ///Kaon-:  
      Double_t c2WeightKaonNeg   =  NumOfKaonNegEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantKaonNeg =  (sumQ2xKaonNegEtaPos[j]*fSumTPCQn2xNeg + sumQ2yKaonNegEtaPos[j]*fSumTPCQn2yNeg)/c2WeightKaonNeg;
      fHistv2AchKaonNeg[1][iCent][Ach]->Fill(ptcenter[j], c2cumulantKaonNeg, c2WeightKaonNeg);   /// for denominator
    }
  }

  if(fParticle==3){
    if(NumOfProtPosEtaPos[j]>0){
      ///Prot+:  
      Double_t c2WeightProtPos   =  NumOfProtPosEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantProtPos =  (sumQ2xProtPosEtaPos[j]*fSumTPCQn2xNeg + sumQ2yProtPosEtaPos[j]*fSumTPCQn2yNeg)/c2WeightProtPos;
      fHistv2AchProtPos[1][iCent][Ach]->Fill(ptcenter[j], c2cumulantProtPos, c2WeightProtPos);   /// for denominator
    }

    if(NumOfProtNegEtaPos[j]>0){
      ///Prot-:  
      Double_t c2WeightProtNeg   =  NumOfProtNegEtaPos[j]*fSumWgtEtaNeg;
      Double_t c2cumulantProtNeg =  (sumQ2xProtNegEtaPos[j]*fSumTPCQn2xNeg + sumQ2yProtNegEtaPos[j]*fSumTPCQn2yNeg)/c2WeightProtNeg;
      fHistv2AchProtNeg[1][iCent][Ach]->Fill(ptcenter[j], c2cumulantProtNeg, c2WeightProtNeg);   /// for denominator
    }
  }
  
    }
  
  //Last lines in Event loop
  fCentDistAfterCut->Fill(centrality);
  fHistEventCount->Fill(14.5); //final Event count.
  PostData(1,fListHist);

  // std::cout<<" Info:UserExec()  Call Finished ..!!!\n";
}//---------------- UserExec ----------------------














void AliAnalysisTaskCVEUNIpt::SetupEventAndTaskConfigInfo(){
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




Int_t AliAnalysisTaskCVEUNIpt::GetCentralityScaled0to10(Double_t fCent){

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

void AliAnalysisTaskCVEUNIpt::GetNUACorrectionHist(Int_t run, Int_t kParticleID)
{

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
}



////---------- SetUp Tracking Efficiency Correction Map ---------------

void AliAnalysisTaskCVEUNIpt::GetMCCorrectionHist(Int_t run){

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
