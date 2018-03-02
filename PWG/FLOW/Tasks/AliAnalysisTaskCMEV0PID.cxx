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
  fHistTaskConfigParameters(NULL),
  fHistPileUpCount(NULL),
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
  fSkipOutlierCut(0),
  fFilterBit(1),
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
  for(int i=0;i<3;i++){
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
AliAnalysisTaskCMEV0PID::AliAnalysisTaskCMEV0PID():AliAnalysisTaskSE(),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  mfileFBHijing(NULL),
  fListFBHijing(NULL),
  fHistTaskConfigParameters(NULL),
  fHistPileUpCount(NULL),
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
  fSkipOutlierCut(0),
  fFilterBit(1),
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
  for(int i=0;i<3;i++){
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
  if(mfileFBHijing->IsOpen()){
     mfileFBHijing->Close();
     if(fListFBHijing) delete fListFBHijing;
  }
  for(int i=0;i<10;i++){
    if(fFB_Efficiency_Cent[i])
      delete fFB_Efficiency_Cent[i];
  }
  if(fAnalysisUtil){
    delete fAnalysisUtil; // its 'new' !!
  }
}//---------------- sanity ------------------------



void AliAnalysisTaskCMEV0PID::UserCreateOutputObjects()
{
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
    fSkipOutlierCut = 1;
  }
  else if(sNucleiTP=="pPb"||sNucleiTP=="Pbp"||sNucleiTP=="PbP"||sNucleiTP=="PPb"){  
    gMaxGlobalmult  = 400;
    gMaxTPCFB1mult  = 400;
    gMaxTPCcorrmult = 500;
    gMaxESDtracks   = 2000;
    fSkipOutlierCut = 1;
  }
  else{
    gMaxGlobalmult  = 4000;
    gMaxTPCFB1mult  = 4000;
    gMaxTPCcorrmult = 5000;
    gMaxESDtracks   = 20000;
    fSkipOutlierCut =  0;
  }



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


  //---------------- PID Histograms ---------------------
  //Dont forget to add histograms in the List. !!!


  char const *gSpecies[3] = {"Pion","Kaon","proton"};

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



  
  
 //-------------------------- Define NUA Hist for PID ------------------------------------
  Int_t gCentForNUA[6] = {0,5,10,20,40,90};
  Char_t  name[100];
  Char_t title[100];

  for(int i=0;i<3;i++){
    for(int j=0;j<5;j++){
      sprintf(name,"fHistEtaPhiVz_%s_Pos_Cent%d_Run%d",gSpecies[i],j,1);
      sprintf(title,"eta,phi,Vz %sPos, Cent%d-%d, FB %d",gSpecies[i],gCentForNUA[j],gCentForNUA[j+1],fFilterBit);
      fHist3DEtaPhiVz_Pos_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
      fListHist->Add(fHist3DEtaPhiVz_Pos_Run[i][j]);

      sprintf(name,"fHistEtaPhiVz_%s_Neg_Cent%d_Run%d",gSpecies[i],j,1);
      sprintf(title,"eta,phi,Vz %sNeg, Cent%d-%d, FB %d",gSpecies[i],gCentForNUA[j],gCentForNUA[j+1],fFilterBit);
      fHist3DEtaPhiVz_Neg_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.283185,16,-0.8,0.8); 
      fListHist->Add(fHist3DEtaPhiVz_Neg_Run[i][j]);
    }
  }
 //---------------------------------------------------------------------------------




  PostData(1,fListHist);
  std::cout<<"........ UserCreateOutputObject called ... fFilterBit = "<<fFilterBit<<" Cmax = "<<fCentralityPercentMax <<"\n"<<std::endl;
}







//______________________________________________________________________
void AliAnalysisTaskCMEV0PID::UserExec(Option_t*){
  //printf("info: UserExec is called.... 1  \n");

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


  Int_t ntracks=fAOD->GetNumberOfTracks();

  if(ntracks<2) return;              // Check this cut....!!!
    

  fHistEventCount->Fill(stepCount); //5
  stepCount++;


  Int_t cent10bin = -1;

  cent10bin = GetCentralityScaled0to10(centrality); //Centrality in 0-10 scale



  Int_t cForNUA = 0;  //Centrality array index for NUA correcion

 if(centrality<5.0) {
   cForNUA = 0;
 }
 else if(centrality>=5.0 && centrality<10){
   cForNUA = 1; // 1=5-10,
 }
 else if(centrality>=10.0 && centrality<20) {
   cForNUA = 2; // 2 = 10-20,
 }
 else if(centrality>=20 && centrality<40.0){ 
   cForNUA = 3; // 3=20-40
 }
 else if(centrality>=40.0){
   cForNUA = 4; // 4=40-90
 }
 


 //Variables for MC tracking correction 
  Int_t   ptBinMC = 1;
  Float_t ptWgtMC = 1.0;
  Float_t ptTrk  = 0.1;


 //-------------- Track loop for outlier and PileUp cut -------------------

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

  //---------------- a dobrin --------------

  Bool_t  bIsOutLier=kTRUE;

  Float_t multTPC     = 0;    // tpc mult estimate
  Float_t RefMultRaw  = 0;    // tpc mult estimate
  Float_t RefMultCorr = 0;    // tpc mult estimate
  Float_t RefMultRawFB = 0;
  Float_t RefMultCorrFB= 0;

  Float_t multTPCAll  = 0;    // tpc mult estimate
  Float_t multGlobal  = 0; // global multiplicity
    

//const Int_t nGoodTracks = fVevent->GetNumberOfTracks();
  const Int_t nGoodTracks = fAOD->GetNumberOfTracks();

  for(Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++) { //-------------------------

    AliAODTrack* AODtrack =dynamic_cast<AliAODTrack*>(fVevent->GetTrack(iTrack));
    if(!AODtrack) continue;

    if(AODtrack->TestFilterBit(128)) multTPCAll++; //A. Dobrin, no track cuts

    if(!(AODtrack->TestFilterBit(1))) continue;
    if((AODtrack->Pt() < .2) || (TMath::Abs(AODtrack->Eta()) > .8) || (AODtrack->GetTPCNcls() < 70)  || (AODtrack->GetDetPid()->GetTPCsignal() < 10.0) || (AODtrack->Chi2perNDF() < 0.2)) continue;
    multTPC++;

    ptTrk   = AODtrack->Pt();

    if(ptTrk >= 0.2 && ptTrk <= 5.0) {
      if(fFB_Efficiency_Cent[cent10bin]){
	ptBinMC = fFB_Efficiency_Cent[cent10bin]->FindBin(ptTrk);
	ptWgtMC    = 1.0/fFB_Efficiency_Cent[cent10bin]->GetBinContent(ptBinMC);
        //cout<<iTrack<<"cent = "<<cent10bin<<" pt = "<<ptTrk<<"\t bin = "<<ptBinMC<<"\t Wgt = "<<ptWgtMC<<endl;    
	RefMultRaw++;
	RefMultCorr += ptWgtMC;
        if((AODtrack->TestFilterBit(fFilterBit))){
	  RefMultRawFB++;
	  RefMultCorrFB += ptWgtMC;
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
    if(!fMultSelection->GetThisEventIsNotPileup()) bIsPileup=kTRUE;
    if(!fMultSelection->GetThisEventIsNotPileupMV()) bIsPileup=kTRUE;
    if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) bIsPileup=kTRUE;
    if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) bIsPileup=kTRUE;
    if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) bIsPileup=kTRUE;
    if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) bIsPileup=kTRUE;
    if(!fMultSelection->GetThisEventHasGoodVertex2016()) bIsPileup=kTRUE;
    if(bIsPileup) fHistPileUpCount->Fill(9.5);
  }  
  //-----------------------------------------------------------------



  fHistTPCVsESDTrkBefore->Fill(multTPCAll,multEsd);   //A. Dobrin

  fHistTPCvsGlobalMultBefore->Fill(multGlobal,multTPC);
     

        
//if(multTPC < (-40.3+1.22*multGlobal) || multTPC > (32.1+1.59*multGlobal)) { bIsOutLier = kFALSE;} //Run-1 or 2011
  if(multTPC < (-20.0+1.15*multGlobal) || multTPC > (200.+1.40*multGlobal)) { bIsOutLier = kFALSE;}
       
  fHistEventCount->Fill(stepCount); //6
  stepCount++;


  fHistTPConlyVsCL1Before->Fill(centrCL1,multTPCAll);
  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTPCAll);
  fHistGlobalVsV0MBefore->Fill(centrV0M, multGlobal);







  if(!fSkipOutlierCut && bIsOutLier) return; //outlier TPC vs Global

  fHistTPCvsGlobalMultAfter->Fill(multGlobal,multTPC);

  fHistEventCount->Fill(stepCount); //7
  stepCount++;


  if(!fSkipOutlierCut && bIsPileup) return;         //PileUp A. Dobrin

  fHistTPCVsESDTrkAfter->Fill(multTPCAll,multEsd);  

  fHistEventCount->Fill(stepCount); //8
  stepCount++;










  Int_t icentBin = centrality;
  icentBin++;

  Float_t TPCmultLowLimit  =  hCentvsTPCmultCuts->GetBinContent(icentBin,1);
  Float_t TPCmultHighLimit =  hCentvsTPCmultCuts->GetBinContent(icentBin,1);

  TPCmultLowLimit  -=  5.0 * hCentvsTPCmultCuts->GetBinContent(icentBin,2); //mean - 5sigma
  TPCmultHighLimit +=  5.0 * hCentvsTPCmultCuts->GetBinContent(icentBin,2); //mean + 5sigma
  //std::cout<<" Cent = "<<centrality<<"\t icent = "<<icentBin<<" low = "<<TPCmultLowLimit<<"\t high = "<<TPCmultHighLimit<<std::endl;

  if(!fSkipOutlierCut){
    if(multTPC<TPCmultLowLimit || multTPC>TPCmultHighLimit) return;
  }

  fHistEventCount->Fill(stepCount); //9
  stepCount++;



  fHistTPConlyVsCL1After->Fill(centrCL1,multTPCAll);
  fHistTPConlyVsV0MAfter->Fill(centrV0M,multTPCAll);
  fHistGlobalVsV0MAfter->Fill(centrV0M, multGlobal);


  // MC corrected Refmult:
  fHistRawVsCorrMultFB->Fill(RefMultRawFB,RefMultCorrFB); // FB set by AddTask.. 

  //--------------------------------------------------------
















  //if(nGoodTracks!=ntracks) std::cout<<" Event with ntracks = "<<ntracks<<"\t good Tracks = "<<nGoodTracks<<std::endl;

            






















  //--------- PID works begin ----------
  Double_t   PDGmassPion   = 0.13957;
  Double_t   PDGmassKaon   = 0.49368;
  Double_t   PDGmassProton = 0.93827;

  PDGmassProton *= PDGmassProton;
  PDGmassPion   *= PDGmassPion;
  PDGmassKaon   *= PDGmassKaon;
 
  Double_t mass=0,mom = -999, pT = -999, phi = -999, eta = -999, dEdx =-999;
  Double_t length = -999., beta =-999, tofTime = -999., tof = -999.;
  Double_t nSigTOFpion, nSigTPCpion;
  Double_t nSigTOFkaon, nSigTPCkaon;
  Double_t nSigTOFproton, nSigTPCproton;
  Double_t c = TMath::C()*1.E-9;// velocity of light m/ns 
  Float_t  dcaXY, dcaZ;
  Float_t  probMis; 
  Int_t    TOFmatch=0; 
  Int_t    charge;
        
  //--------- Track Loop for PID studies ----------------

  for(Int_t itrack = 0; itrack < ntracks; itrack++) {

    AliAODTrack *track=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(itrack));
    if(!track) continue;

    dcaXY = track->DCA();
    dcaZ  = track->ZAtDCA();
        
    //cout<<"track->GetFilterMap()= "<<track->GetFilterMap()<<endl;
    if(!track->TestFilterBit(fFilterBit)) continue;
        
    mom=track->P();
    pT=track->Pt();
    phi=track->Phi();
    eta=track->Eta();
    dEdx=track->GetDetPid()->GetTPCsignal();
    charge = track->Charge();
    fHistEtaPtBefore->Fill(eta,pT);





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

    //if(probMis < 0.20)  // if u want to reduce mismatch using also TPC
    if(TOFmatch>0 && probMis < 0.01) { 
      //This conditions are called when detector status is checked above :  
      //if((track->IsOn(AliAODTrack::kTOFin)) && (track->IsOn(AliAODTrack::kTIME)) &&  (track->IsOn(AliAODTrack::kTOFout))) {
      //if((track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid))) { //Naghmeh used it
                
	tofTime = track->GetTOFsignal();       // in pico seconds
	length  = track->GetIntegratedLength();
                
	tof = tofTime*1.E-3; // ns
	if (tof <= 0) continue;            
	//if (length <= 0) continue;
            
	length = length*0.01; // in meters
	tof = tof*c;
	beta = length/tof;
	mass = mom*mom*(1./(beta*beta) - 1);
	//cout<<"tof = "<<tof<<"\t length = "<<length<<"\t beta = "<<beta<<endl;                    
    }//TOF signal
    

    fHistTPCdEdxvsPBefore->Fill(mom*charge,dEdx);
    fHistTOFBetavsPBefore->Fill(mom*charge,beta);
    fHistTOFMassvsPtBefore->Fill(pT*charge,mass);




    //-------- Apply all track cuts for analsys: ---------
    if(pT < fMinPtCut  || pT > fMaxPtCut)     continue;
    if(eta < fMinEtaCut || eta > fMaxEtaCut)  continue;
    if(!(track->TestFilterBit(fFilterBit)))   continue;
    //-----------------------------------------------------




    fHistEtaPtAfter->Fill(eta,pT);
    fHistTPCdEdxvsPAfter->Fill(mom*charge,dEdx);

    if(TOFmatch>0 && probMis < 0.01){
      fHistTOFBetavsPAfter->Fill(mom*charge,beta);
    }

    nSigTOFpion=fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
    nSigTOFkaon=fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
    nSigTOFproton=fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

    fHistTOFMatchCount->Fill(TOFmatch+2,nSigTOFpion);
    fHistTOFMatchCount->Fill(TOFmatch+4,nSigTOFkaon);
    fHistTOFMatchCount->Fill(TOFmatch+6,nSigTOFproton);

    if(!TOFmatch || probMis > 0.01){  // I dont want mismatched track in my signal distribution
      nSigTOFpion   = -9.99;
      nSigTOFkaon   = -9.99;
      nSigTOFproton = -9.99;
    }

    nSigTPCpion   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    nSigTPCkaon   = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    nSigTPCproton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

    //0=pi, 1=K, 2=Proton
    fHistTPCTOFnSigmavsPtAfter[0]->Fill(pT,nSigTPCpion,nSigTOFpion);
    fHistTPCTOFnSigmavsPtAfter[1]->Fill(pT,nSigTPCkaon,nSigTOFkaon);
    fHistTPCTOFnSigmavsPtAfter[2]->Fill(pT,nSigTPCproton,nSigTOFproton);





    /*
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
    */



    // nSigmaTOF distribution for circular cut

    if(TMath::Sqrt(nSigTPCpion*nSigTPCpion+nSigTOFpion*nSigTOFpion)<=fNSigmaCut){
      fHistTOFnSigmavsPtAfter[0]->Fill(pT*charge,nSigTOFpion);
      if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        fHistPtwithTOFSignal[0]->Fill(pT*charge);
      }
    }

    if(TMath::Sqrt(nSigTPCkaon*nSigTPCkaon+nSigTOFkaon*nSigTOFkaon)<=fNSigmaCut){
      fHistTOFnSigmavsPtAfter[1]->Fill(pT*charge,nSigTOFkaon);
      if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        fHistPtwithTOFSignal[1]->Fill(pT*charge);
      }
    }

    if(TMath::Sqrt(nSigTPCproton*nSigTPCproton+nSigTOFproton*nSigTOFproton)<=fNSigmaCut){
      fHistTOFnSigmavsPtAfter[2]->Fill(pT*charge,nSigTOFproton);
      if(TOFmatch>0 && probMis < 0.01 && beta>0.2){
        fHistPtwithTOFSignal[2]->Fill(pT*charge);
      }
    }


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



    //============== Fill NUA Histograms for Pion ---------------------
     if(pT<0.6 && TMath::Abs(nSigTPCpion)<=2.5){
       if(charge>0){
	 fHist3DEtaPhiVz_Pos_Run[0][cForNUA]->Fill(pVtxZ,phi,eta);
       }
       else if(charge<0){
	 fHist3DEtaPhiVz_Neg_Run[0][cForNUA]->Fill(pVtxZ,phi,eta);
       }
     }
     else if(pT>=0.6 && pT<=2.0 && TMath::Abs(nSigTPCpion)<=2.0 && TMath::Abs(nSigTOFpion)<=2.0 ){
       if(charge>0){
	 fHist3DEtaPhiVz_Pos_Run[0][cForNUA]->Fill(pVtxZ,phi,eta);
       }
       else if(charge<0){
	 fHist3DEtaPhiVz_Neg_Run[0][cForNUA]->Fill(pVtxZ,phi,eta);
       }
     }
     else if(pT>=2.0 && pT<=4.0  && TMath::Abs(nSigTPCpion)<=1.0 && TMath::Abs(nSigTOFpion)<=1.0){
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
     else if(pT>=0.6 && pT<=3.0 && TMath::Abs(nSigTPCkaon)<=1.5 && TMath::Abs(nSigTOFkaon)<=1.5){
       if(charge>0){
	 fHist3DEtaPhiVz_Pos_Run[1][cForNUA]->Fill(pVtxZ,phi,eta);
       }
       else if(charge<0){
	 fHist3DEtaPhiVz_Neg_Run[1][cForNUA]->Fill(pVtxZ,phi,eta);
       }
     }
     else if(pT>=3.0 && pT<=4.0 && TMath::Abs(nSigTPCkaon)<=1.0 && TMath::Abs(nSigTOFkaon)<=1.0){
       if(charge>0){
	 fHist3DEtaPhiVz_Pos_Run[1][cForNUA]->Fill(pVtxZ,phi,eta);
       }
       else if(charge<0){
	 fHist3DEtaPhiVz_Neg_Run[1][cForNUA]->Fill(pVtxZ,phi,eta);
       }
     }



    //============== Fill NUA Histograms for proton ---------------------
     if(pT<0.9 && TMath::Abs(nSigTPCproton)<=2.5){
       if(charge>0){
	 fHist3DEtaPhiVz_Pos_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
       }
       else if(charge<0){
	 fHist3DEtaPhiVz_Neg_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
       }
     }
     else if(pT>=0.9 && pT<=3.5 && TMath::Abs(nSigTPCproton)<=2.5 && TMath::Abs(nSigTOFproton)<=2.0){
       if(charge>0){
	 fHist3DEtaPhiVz_Pos_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
       }
       else if(charge<0){
	 fHist3DEtaPhiVz_Neg_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
       }
     }
     else if(pT>=3.5 && pT<=4.0 && TMath::Abs(nSigTPCproton)<=1.5 && TMath::Abs(nSigTOFproton)<=1.0){
       if(charge>0){
	 fHist3DEtaPhiVz_Pos_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
       }
       else if(charge<0){
	 fHist3DEtaPhiVz_Neg_Run[2][cForNUA]->Fill(pVtxZ,phi,eta);
       }
     }
 


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

  }//track loop ends






  PostData(1,fListHist);

  fHistEventCount->Fill(stepCount); //10
  stepCount++;

}//================ UserExec ==============







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
    if(bApplyMCcorr){ printf("\n\n!!*****  Warning *****!!\n MC correction File not found, using Efficiency = 1.0 !!\n\n");}
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




  fHistEventCount = new TH1F("fHistEventCount","Event counts",15,0,15);
  fHistEventCount->GetXaxis()->SetBinLabel(1,"Called UserExec()");
  fHistEventCount->GetXaxis()->SetBinLabel(2,"Called Exec()");
  fHistEventCount->GetXaxis()->SetBinLabel(3,"AOD Exist");
  fHistEventCount->GetXaxis()->SetBinLabel(4,"Vz<10");
  fHistEventCount->GetXaxis()->SetBinLabel(5,Form("%2.0f<Cent<%2.0f",fCentralityPercentMin,fCentralityPercentMax));
  fHistEventCount->GetXaxis()->SetBinLabel(6,"noAODtrack>2 ");
  fHistEventCount->GetXaxis()->SetBinLabel(7,"TPC vs Global");
  fHistEventCount->GetXaxis()->SetBinLabel(8,"TPC128 vs ESD");
  fHistEventCount->GetXaxis()->SetBinLabel(9,"Cent vs TPC");
  fHistEventCount->GetXaxis()->SetBinLabel(10,"Survived Events");
  fListHist->Add(fHistEventCount);

  //fHistEventCount->Fill(1);

}



