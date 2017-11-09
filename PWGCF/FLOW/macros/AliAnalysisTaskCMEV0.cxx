/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////
// CME correlator study using VO.
// Author: Rihan Haque (mhaque@cern.ch)
///////////////////////////////////////////////



#include <TGrid.h>
#include <TFile.h>
#include <TList.h>
#include "TMatrixDSym.h"
#include "Riostream.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"
#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliVVertex.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliFlowEvent.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskCMEV0.h"

using std::endl;
using std::cout;

ClassImp(AliAnalysisTaskCMEV0)

AliAnalysisTaskCMEV0::AliAnalysisTaskCMEV0(const TString name):AliAnalysisTaskSE(name),
  fEvent(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHistos(NULL),
  fListCalibs(NULL),
  fListFBHijing(NULL),
  fListNUACorr(NULL),
  fListZDNCorr(NULL),
  fRejectPileUp(kTRUE),
  fRejectPileUpTight(kTRUE),
  bFillAvgTPCQn(kFALSE),
  bFillEtaPhiNUA(kFALSE),
  bApplyNUACorr(kFALSE),
  bApplyZDCCorr(kFALSE),
  bApplyNUAforEP(kFALSE),
  sDataSet("2015"),
  sAnalysisSet("DoGainEq"),
  sCentEstimator("V0"),
  fRunFlag(0),
  fOldRunNum(0),
  fievent(0),
  EvtCent(0.0),
  fHarmonic(2.0),
  fHist_Event_count(NULL),
  fPileUpMultSelCount(NULL),
  fPileUpCount(NULL),
  fMultV0(NULL),
  fQxnmV0A(NULL), 
  fQynmV0A(NULL),     
  fQxnsV0A(NULL),       
  fQynsV0A(NULL),       
  fQxnmV0C(NULL),        
  fQynmV0C(NULL),    
  fQxnsV0C(NULL),     
  fQynsV0C(NULL),
  fHCos1nPosChEtaPosVzPos(NULL),
  fHCos1nPosChEtaNegVzPos(NULL),
  fHCos1nPosChEtaPosVzNeg(NULL),
  fHCos1nPosChEtaNegVzNeg(NULL),
  fHSin1nPosChEtaPosVzPos(NULL),
  fHSin1nPosChEtaNegVzPos(NULL),
  fHSin1nPosChEtaPosVzNeg(NULL),
  fHSin1nPosChEtaNegVzNeg(NULL),
  fHCos1nNegChEtaPosVzPos(NULL),
  fHCos1nNegChEtaNegVzPos(NULL),
  fHCos1nNegChEtaPosVzNeg(NULL),
  fHCos1nNegChEtaNegVzNeg(NULL),
  fHSin1nNegChEtaPosVzPos(NULL),
  fHSin1nNegChEtaNegVzPos(NULL),
  fHSin1nNegChEtaPosVzNeg(NULL),
  fHSin1nNegChEtaNegVzNeg(NULL),
  fHV0AEventPlaneVsCent(NULL),
  fHV0CEventPlaneVsCent(NULL),
  fHTPCEventPlaneVsCent(NULL),
  fHCorrectNUApos(NULL),
  fHCorrectNUAneg(NULL),
  fHEnergyZNCvsCent(NULL),
  fHEnergyZNAvsCent(NULL),
  fHEnergyZPCvsCent(NULL),
  fHEnergyZPAvsCent(NULL),
  fHEnergyZNCvsCentRun(NULL),
  fHEnergyZNAvsCentRun(NULL),
  fHEnergyZPCvsCentRun(NULL),
  fHEnergyZPAvsCentRun(NULL),
  fHEnergyZPCvsZPA(NULL),
  fHEnergyZNCvsZNA(NULL),
  hUnderOverBinNUApos(NULL),
  hUnderOverBinNUAneg(NULL),
  fHCentBinTrkRecenter(NULL),
  fHCorrectZDNP(NULL),
  fV0AQnxVsCentRun(NULL),
  fV0AQnyVsCentRun(NULL),
  fV0CQnxVsCentRun(NULL),
  fV0CQnyVsCentRun(NULL),
  fTPCQnxVsCentRun(NULL),
  fTPCQnyVsCentRun(NULL),
  mListNUAPos(NULL),
  mListNUANeg(NULL),
  fileNUApos(NULL),
  fileNUAneg(NULL)
{
  for(int i=0;i<90;i++){
    runNums[i] = 0;
    for(int j=0;j<4;j++){
     fHist3DEtaPhiVz_Pos_Run[j][i] = NULL;
     fHist3DEtaPhiVz_Neg_Run[j][i] = NULL;
    }
  }
  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
  }
  for(int i=0;i<3;i++){
    fHist_Corr3p_SP_Norm_PN[i] =  NULL;
    fHist_Corr3p_SP_Norm_PP[i] =  NULL;
    fHist_Corr3p_SP_Norm_NN[i] =  NULL;
    fHist_Reso2n_SP_Norm_Det[i] =  NULL;
  }
  for(int i=0;i<3;i++){
    fHist_Corr3p_EP_Norm_PN[i] =  NULL;
    fHist_Corr3p_EP_Norm_PP[i] =  NULL;
    fHist_Corr3p_EP_Norm_NN[i] =  NULL;
    fHist_Reso2n_EP_Norm_Det[i] =  NULL;
  }
  for(int i=0;i<3;i++){
    fHist_Corr3p_ZDN_SP_PN[i]  =  NULL;
    fHist_Corr3p_ZDN_SP_PP[i]  =  NULL;
    fHist_Corr3p_ZDN_SP_NN[i]  =  NULL;
    fHist_Reso2n_ZDN_SP_Det[i] = NULL;
  }
  for(int i=0;i<6;i++){
    fHist_Corr3p_pTSum_EP_V0A_PN[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0A_PP[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0A_NN[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0C_PN[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0C_PP[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0C_NN[i] = NULL;
  }
  for(int i=0;i<6;i++){
    fHist_Corr3p_pTDiff_EP_V0A_PN[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0A_PP[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0A_NN[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0C_PN[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0C_PP[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0C_NN[i] = NULL;
  }
  for(int i=0;i<6;i++){
    fHist_Corr3p_EtaDiff_EP_V0A_PN[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0A_PP[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0A_NN[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0C_PN[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0C_PP[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0C_NN[i] = NULL;
  }
  for(int i=0;i<10;i++){
    fHistChPosvsEtaPtRun[i] = NULL;
    fHistChNegvsEtaPtRun[i] = NULL;
  }


  DefineInput(1, AliFlowEventSimple::Class()); // Input slot #1: AliFlowEventSimple

  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());
}//-------- real ---------

AliAnalysisTaskCMEV0::AliAnalysisTaskCMEV0(): AliAnalysisTaskSE(),
  fEvent(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHistos(NULL),
  fListCalibs(NULL),
  fListFBHijing(NULL),
  fListNUACorr(NULL),
  fListZDNCorr(NULL),
  fRejectPileUp(kTRUE),
  fRejectPileUpTight(kTRUE),
  bFillAvgTPCQn(kFALSE),
  bFillEtaPhiNUA(kFALSE),
  bApplyNUACorr(kFALSE),
  bApplyZDCCorr(kFALSE),
  bApplyNUAforEP(kFALSE),
  sDataSet("2015"),
  sAnalysisSet("DoGainEq"),
  sCentEstimator("V0"),
  fRunFlag(0),
  fOldRunNum(0),
  fievent(0),
  EvtCent(0.0),
  fHarmonic(2.0),
  fHist_Event_count(NULL),
  fPileUpMultSelCount(NULL),
  fPileUpCount(NULL),
  fMultV0(NULL),
  fQxnmV0A(NULL), 
  fQynmV0A(NULL),     
  fQxnsV0A(NULL),       
  fQynsV0A(NULL),       
  fQxnmV0C(NULL),        
  fQynmV0C(NULL),    
  fQxnsV0C(NULL),     
  fQynsV0C(NULL),
  fHCos1nPosChEtaPosVzPos(NULL),
  fHCos1nPosChEtaNegVzPos(NULL),
  fHCos1nPosChEtaPosVzNeg(NULL),
  fHCos1nPosChEtaNegVzNeg(NULL),
  fHSin1nPosChEtaPosVzPos(NULL),
  fHSin1nPosChEtaNegVzPos(NULL),
  fHSin1nPosChEtaPosVzNeg(NULL),
  fHSin1nPosChEtaNegVzNeg(NULL),
  fHCos1nNegChEtaPosVzPos(NULL),
  fHCos1nNegChEtaNegVzPos(NULL),
  fHCos1nNegChEtaPosVzNeg(NULL),
  fHCos1nNegChEtaNegVzNeg(NULL),
  fHSin1nNegChEtaPosVzPos(NULL),
  fHSin1nNegChEtaNegVzPos(NULL),
  fHSin1nNegChEtaPosVzNeg(NULL),
  fHSin1nNegChEtaNegVzNeg(NULL),
  fHV0AEventPlaneVsCent(NULL),
  fHV0CEventPlaneVsCent(NULL),
  fHTPCEventPlaneVsCent(NULL),
  fHCorrectNUApos(NULL),
  fHCorrectNUAneg(NULL),
  fHEnergyZNCvsCent(NULL),
  fHEnergyZNAvsCent(NULL),
  fHEnergyZPCvsCent(NULL),
  fHEnergyZPAvsCent(NULL),
  fHEnergyZNCvsCentRun(NULL),
  fHEnergyZNAvsCentRun(NULL),
  fHEnergyZPCvsCentRun(NULL),
  fHEnergyZPAvsCentRun(NULL),
  fHEnergyZPCvsZPA(NULL),
  fHEnergyZNCvsZNA(NULL),
  hUnderOverBinNUApos(NULL),
  hUnderOverBinNUAneg(NULL),
  fHCentBinTrkRecenter(NULL),
  fHCorrectZDNP(NULL),
  fV0AQnxVsCentRun(NULL),
  fV0AQnyVsCentRun(NULL),
  fV0CQnxVsCentRun(NULL),
  fV0CQnyVsCentRun(NULL),
  fTPCQnxVsCentRun(NULL),
  fTPCQnyVsCentRun(NULL),
  mListNUAPos(NULL),
  mListNUANeg(NULL),
  fileNUApos(NULL),
  fileNUAneg(NULL)
{
  for(int i=0;i<90;i++){
    runNums[i] = 0;
    for(int j=0;j<4;j++){
     fHist3DEtaPhiVz_Pos_Run[j][i] = NULL;
     fHist3DEtaPhiVz_Neg_Run[j][i] = NULL;
    }
  }
  for(int i=0;i<10;i++){
    fFB_Efficiency_Cent[i] = NULL;
  }
  for(int i=0;i<3;i++){
    fHist_Corr3p_SP_Norm_PN[i] =  NULL;
    fHist_Corr3p_SP_Norm_PP[i] =  NULL;
    fHist_Corr3p_SP_Norm_NN[i] =  NULL;
    fHist_Reso2n_SP_Norm_Det[i] =  NULL;
  }
  for(int i=0;i<3;i++){
    fHist_Corr3p_EP_Norm_PN[i] =  NULL;
    fHist_Corr3p_EP_Norm_PP[i] =  NULL;
    fHist_Corr3p_EP_Norm_NN[i] =  NULL;
    fHist_Reso2n_EP_Norm_Det[i] =  NULL;
  }
  for(int i=0;i<3;i++){
    fHist_Corr3p_ZDN_SP_PN[i]  =  NULL;
    fHist_Corr3p_ZDN_SP_PP[i]  =  NULL;
    fHist_Corr3p_ZDN_SP_NN[i]  =  NULL;
    fHist_Reso2n_ZDN_SP_Det[i] = NULL;
  }
  for(int i=0;i<6;i++){
    fHist_Corr3p_pTSum_EP_V0A_PN[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0A_PP[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0A_NN[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0C_PN[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0C_PP[i] = NULL;
    fHist_Corr3p_pTSum_EP_V0C_NN[i] = NULL;
  }
  for(int i=0;i<6;i++){
    fHist_Corr3p_pTDiff_EP_V0A_PN[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0A_PP[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0A_NN[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0C_PN[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0C_PP[i] = NULL;
    fHist_Corr3p_pTDiff_EP_V0C_NN[i] = NULL;
  }
  for(int i=0;i<6;i++){
    fHist_Corr3p_EtaDiff_EP_V0A_PN[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0A_PP[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0A_NN[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0C_PN[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0C_PP[i] = NULL;
    fHist_Corr3p_EtaDiff_EP_V0C_NN[i] = NULL;
  }
  for(int i=0;i<10;i++){
    fHistChPosvsEtaPtRun[i] = NULL;
    fHistChNegvsEtaPtRun[i] = NULL;
  }


}//-----------




void AliAnalysisTaskCMEV0::UserCreateOutputObjects()
{
 this->InitializeRunArray(sDataSet);

 fAnalysisUtil = new AliAnalysisUtils();
 fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
 fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);


 if(fListFBHijing) {
  for(int i=0;i<10;i++) {
    fFB_Efficiency_Cent[i] = (TH1D *) fListFBHijing->FindObject(Form("eff_unbiased_%d",i));
  }
 }
 else{ // if MC efficiency not found then use weight = 1.
  printf("\n\n!!*****  Warning *****!!\n FilterBit efficiency not found, use = 1.0 !!\n\n");
   for(int i=0;i<10;i++){
     fFB_Efficiency_Cent[i] = new TH1D(Form("eff_unbiased_%d",i),"",1000,0,50.); 
     for(int j=1;j<=fFB_Efficiency_Cent[i]->GetNbinsX();j++){
       fFB_Efficiency_Cent[i]->SetBinContent(j,1.000);
     }
   } 
 }



 fListHistos = new TList();
 fListHistos->SetOwner(kTRUE);

 fListCalibs = new TList();
 fListCalibs->SetOwner(kTRUE);

 this->DefineHistograms();

 PostData(1,fListHistos); 
 PostData(2,fListCalibs); 

}

AliAnalysisTaskCMEV0::~AliAnalysisTaskCMEV0()
{
  delete                 fListHistos;        
  delete                 fListCalibs;         

  delete              fMultSelection; 
  delete               fAnalysisUtil; // it is '= new' !!!

  delete        fHCentBinTrkRecenter;

}







void AliAnalysisTaskCMEV0::UserExec(Option_t *)
{

 Float_t stepCount = 0.5;
 fHist_Event_count->Fill(stepCount); //1
 stepCount++;

 AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(1));

 if(!aod || !fEvent){
   printf("\n ... ::UserExec = no AOD or Flow Event, \n.... EXIT .....  \n");
   return;
 }

 fHist_Event_count->Fill(stepCount); //2
 stepCount++;



 //---------pileup rejection: --------
 Bool_t kPileupEvent =  CheckEventIsPileUp(aod);

 if(kPileupEvent)    return;

 fHist_Event_count->Fill(stepCount); //3
 stepCount++;



 Double_t centrV0M = -1;
 Double_t centrCL1 = -1;
 Double_t centrCL0 = -1;
 Double_t centrTRK = -1;
 EvtCent = -1;

 if(sDataSet=="2010"||sDataSet=="2011"){
    centrV0M = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
    centrCL1 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");
    centrCL0 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL0");
    centrTRK = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("TRK");
 }
 else{
  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
   if(!fMultSelection) {
     printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
     exit(1);
   }
    centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
    centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
    centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");
 }



 if(sCentEstimator=="V0" || sCentEstimator=="V0M"){
   EvtCent = centrV0M;
 }
 else{
   EvtCent = centrCL1;
 }



 Int_t iCentSPD = centrCL1;

 if(iCentSPD >= 90) return;

 fHist_Event_count->Fill(stepCount); //4
 stepCount++;


 Int_t cIndex = 0;
 Int_t cForNUA = 0;

 if(EvtCent<5.0) {
   cIndex  = 0; 
   cForNUA = 0;
 }
 else if(EvtCent>=5.0 && EvtCent<10){
   cIndex  = 1;
   cForNUA = 1;
 }
 else if(EvtCent>=10.0 && EvtCent<40.0) {
   cForNUA = 2;
 }
 else if(EvtCent>=40.0) {
     cForNUA = 3;  
 }





  //------- load v0 calib file of Alex --------
 Int_t runindex = -111;
 Int_t runNumber = aod->GetRunNumber();
 runindex = GetCurrentRunIndex(runNumber);

 if(runNumber!=fOldRunNum){ 
   OpenInfoCalbration(runNumber);
   if(bApplyNUACorr){
     GetNUACorrectionHist(runNumber,EvtCent);
   }
   if(bApplyZDCCorr){
     GetZDCCorrectionHist(runNumber);
   }
   fOldRunNum = runNumber;
 }


//-------- V0 info ---------------
 const AliAODVZERO *fAODV0 = aod->GetVZEROData();

 Double_t Qxan  = 0, Qyan  = 0;
 Double_t Qxcn  = 0, Qycn  = 0;
 Double_t sumMa = 0, sumMc = 0;

 GetV0QvectAndMult(fAODV0, Qxan, Qyan, sumMa, Qxcn, Qycn, sumMc);


 Double_t QyanCor = (Qyan - fQynmV0A->GetBinContent(iCentSPD+1))/fQynsV0A->GetBinContent(iCentSPD+1);
 Double_t QycnCor = (Qycn - fQynmV0C->GetBinContent(iCentSPD+1))/fQynsV0C->GetBinContent(iCentSPD+1);

 Double_t QxanCor = Qxan;
 Double_t QxcnCor = Qxcn;

 if(fHarmonic != 4.){
   QxanCor = (Qxan - fQxnmV0A->GetBinContent(iCentSPD+1))/fQxnsV0A->GetBinContent(iCentSPD+1);
   QxcnCor = (Qxcn - fQxnmV0C->GetBinContent(iCentSPD+1))/fQxnsV0C->GetBinContent(iCentSPD+1);
 }

 if(sumMa < 0 || sumMc < 0)   return;
 fHist_Event_count->Fill(stepCount); //5
 stepCount++;

 QxanCor = QxanCor;
 QyanCor = QyanCor;
 QxcnCor = QxcnCor;
 QycnCor = QycnCor;


 AliAODVertex *pVertex = aod->GetPrimaryVertex();
 Double_t VtxZ    =   pVertex->GetZ();



 //------------- Load ZDC info: -----------
 Double_t energyZNC=0.,energyZNA=0.,energyZPC=0.,energyZPA=0.;

 energyZNC = ((AliVAODHeader*)aod->GetHeader())->GetZDCN1Energy();
 energyZNA = ((AliVAODHeader*)aod->GetHeader())->GetZDCN2Energy();
 energyZPC = ((AliVAODHeader*)aod->GetHeader())->GetZDCP1Energy();
 energyZPA = ((AliVAODHeader*)aod->GetHeader())->GetZDCP2Energy();


 //Apply Gain Correction:
 if(bApplyZDCCorr && fListZDNCorr){
   energyZNC = energyZNC * fHCorrectZDNP->GetBinContent((abs(EvtCent)+1),1);
   energyZNA = energyZNA * fHCorrectZDNP->GetBinContent((abs(EvtCent)+1),2);
   energyZPC = energyZPC * fHCorrectZDNP->GetBinContent((abs(EvtCent)+1),3);
   energyZPA = energyZPA * fHCorrectZDNP->GetBinContent((abs(EvtCent)+1),4);
 }

 fHEnergyZNCvsCent->Fill(EvtCent,energyZNC);
 fHEnergyZNAvsCent->Fill(EvtCent,energyZNA);
 fHEnergyZPCvsCent->Fill(EvtCent,energyZPC);
 fHEnergyZPAvsCent->Fill(EvtCent,energyZPA);
 
 /*
 fHEnergyZNCvsCentRun->Fill(EvtCent,runindex,energyZNC);
 fHEnergyZNAvsCentRun->Fill(EvtCent,runindex,energyZNA);
 fHEnergyZPCvsCentRun->Fill(EvtCent,runindex,energyZPC);
 fHEnergyZPAvsCentRun->Fill(EvtCent,runindex,energyZPA);

 if(EvtCent>30. && EvtCent<40.){
   fHEnergyZPCvsZPA->Fill(energyZPC,energyZPA);
   fHEnergyZNCvsZNA->Fill(energyZNC,energyZNA);
 }
 */




 //Store Qvector for tracks:
 Double_t QxPos[2] = {0.,};  //[0]=n*phi,[1]=2*n*phi
 Double_t QyPos[2] = {0.,};
 Double_t QxNeg[2] = {0.,};  
 Double_t QyNeg[2] = {0.,};

 //Store Qvector for Auto-correlations
 Double_t QxAutoPos[2] = {0.,};  
 Double_t QyAutoPos[2] = {0.,};
 Double_t QxAutoNeg[2] = {0.,};  
 Double_t QyAutoNeg[2] = {0.,};

 //Multiplicity of POIs:
 Double_t MPOIpos = 0;
 Double_t MPOIneg = 0;

 //TPC Event plane Qvectors:
 Double_t QxTPC[2] = {0.,};
 Double_t QyTPC[2] = {0.,};

 Int_t ptBin;
 Int_t nRefMult = 0;
 Int_t n = fHarmonic/2.0;  // Because, cos(nphi1 + nphi2 - 2nPsi)//
 
 Int_t iCentRec = 0;  //centrality bin for track recenter 
 iCentRec = fHCentBinTrkRecenter->FindBin(EvtCent);

 //Int_t nBinVz=1,nBinPhi=1,nBinEta=1;

 //CME using Event plane:
 Double_t pi =  TMath::Pi();
 Double_t Psi2V0C = 0.5*TMath::ATan2(QycnCor,QxcnCor);
 if(Psi2V0C<0.) Psi2V0C += pi;

 Double_t Psi2V0A = 0.5*TMath::ATan2(QyanCor,QxanCor);
 if(Psi2V0A<0.) Psi2V0A += pi;

 fV0AQnxVsCentRun->Fill(centrCL1,runindex,QycnCor); // For V0-Q Recenter QA.
 fV0AQnyVsCentRun->Fill(centrCL1,runindex,QyanCor); 
 fV0CQnxVsCentRun->Fill(centrCL1,runindex,QxcnCor);
 fV0CQnyVsCentRun->Fill(centrCL1,runindex,QycnCor);



 Int_t iBinNUA = 1;

 Int_t iTracks = fEvent->NumberOfTracks();

 AliFlowTrackSimple*   pTrack1 = NULL;
 AliFlowTrackSimple*   pTrack2 = NULL;

 Double_t dPhi1,dPt1,dEta1,dChrg1;
 Double_t dPhi2,dPt2,dEta2,dChrg2;
 Double_t WgtEP = 1.0;

 Double_t ptw1 = 1.0, w1NUA = 1.0;
 Double_t ptw2 = 1.0, w2NUA = 1.0;

 for(int i=0; i<iTracks; i++) {
   pTrack1    =     fEvent->GetTrack(i);
   if(!pTrack1)                continue;
   dPhi1      =          pTrack1->Phi(); //0-2pi range
   dPt1       =          pTrack1-> Pt();
   dEta1      =          pTrack1->Eta();
   dChrg1     =       pTrack1->Charge();

   if(!pTrack1->InPOISelection())
   continue;


   if(bFillEtaPhiNUA){
     if(dChrg1>0){
       fHist3DEtaPhiVz_Pos_Run[cForNUA][runindex]->Fill(VtxZ,dPhi1,dEta1);
     }
     else if(dChrg1<0){
       fHist3DEtaPhiVz_Neg_Run[cForNUA][runindex]->Fill(VtxZ,dPhi1,dEta1);
     }
   }

   if(bApplyNUACorr){
     //get NUA weights: 
     if(dChrg1>0){
        iBinNUA = fHCorrectNUApos->FindBin(VtxZ,dPhi1,dEta1);
        w1NUA = 1.0/fHCorrectNUApos->GetBinContent(iBinNUA);  //Jacopo's 
        //w1NUA = fHCorrectNUApos->GetBinContent(iBinNUA);      //Rihan's
	if(w1NUA > 1e5 ) w1NUA = 1.0;
      }
      else if(dChrg1<0){
        iBinNUA = fHCorrectNUAneg->FindBin(VtxZ,dPhi1,dEta1);
        w1NUA = 1.0/fHCorrectNUAneg->GetBinContent(iBinNUA); //Jacopo: bincontent = 1/Wgts 
        //w1NUA = fHCorrectNUAneg->GetBinContent(iBinNUA);     //Rihan: bincontent = Wgts
	if(w1NUA > 1e5 ) w1NUA = 1.0;
      }
   }

   if(dChrg1>0){
     fHistChPosvsEtaPtRun[cIndex]->Fill(dPt1,dEta1,runindex,w1NUA);
   }
   else if(dChrg1<0){
     fHistChNegvsEtaPtRun[cIndex]->Fill(dPt1,dEta1,runindex,w1NUA);
   }


   ptBin = fFB_Efficiency_Cent[cIndex]->FindBin(dPt1);
   ptw1   = 1.0/fFB_Efficiency_Cent[cIndex]->GetBinContent(ptBin);
  
   if(dChrg1>0){
     QxPos[0] += ptw1*w1NUA*TMath::Cos(n*dPhi1);
     QyPos[0] += ptw1*w1NUA*TMath::Sin(n*dPhi1);
     QxAutoPos[0] += ptw1*ptw1*w1NUA*w1NUA*TMath::Cos(2.*n*dPhi1);
     QyAutoPos[0] += ptw1*ptw1*w1NUA*w1NUA*TMath::Sin(2.*n*dPhi1);
     MPOIpos += ptw1*w1NUA;
   }
   else if(dChrg1<0){
     QxNeg[0] += ptw1*w1NUA*TMath::Cos(n*dPhi1);
     QyNeg[0] += ptw1*w1NUA*TMath::Sin(n*dPhi1);
     QxAutoNeg[0] += ptw1*ptw1*w1NUA*w1NUA*TMath::Cos(2.*n*dPhi1);
     QyAutoNeg[0] += ptw1*ptw1*w1NUA*w1NUA*TMath::Sin(2.*n*dPhi1);
     MPOIneg += ptw1*w1NUA;
   }

   QxTPC[0] += ptw1*w1NUA*TMath::Cos(n*dPhi1);
   QyTPC[0] += ptw1*w1NUA*TMath::Sin(n*dPhi1);
   QxTPC[1] += ptw1*w1NUA*TMath::Cos(2.*n*dPhi1);
   QyTPC[1] += ptw1*w1NUA*TMath::Sin(2.*n*dPhi1);

   nRefMult++;

   if(bFillAvgTPCQn) { //fill TPC Q-vect: 
     if(VtxZ>0 && dEta1>=0){
       if(dChrg1>0){
         fHCos1nPosChEtaPosVzPos->Fill(EvtCent,runindex,ptw1*TMath::Cos(n*dPhi1));
         fHSin1nPosChEtaPosVzPos->Fill(EvtCent,runindex,ptw1*TMath::Sin(n*dPhi1));
       }
       else if(dChrg1<0){
         fHCos1nNegChEtaPosVzPos->Fill(EvtCent,runindex,ptw1*TMath::Cos(n*dPhi1));
         fHSin1nNegChEtaPosVzPos->Fill(EvtCent,runindex,ptw1*TMath::Sin(n*dPhi1));
       }
     }
     else if(VtxZ>0 && dEta1<0){
       if(dChrg1>0){
         fHCos1nPosChEtaNegVzPos->Fill(EvtCent,runindex,ptw1*TMath::Cos(n*dPhi1));
         fHSin1nPosChEtaNegVzPos->Fill(EvtCent,runindex,ptw1*TMath::Sin(n*dPhi1));
       }
       else if(dChrg1<0){
         fHCos1nNegChEtaNegVzPos->Fill(EvtCent,runindex,ptw1*TMath::Cos(n*dPhi1));
         fHSin1nNegChEtaNegVzPos->Fill(EvtCent,runindex,ptw1*TMath::Sin(n*dPhi1));
       }
     }
     else if(VtxZ<0 && dEta1>=0){
       if(dChrg1>0){
         fHCos1nPosChEtaPosVzNeg->Fill(EvtCent,runindex,ptw1*TMath::Cos(n*dPhi1));
         fHSin1nPosChEtaPosVzNeg->Fill(EvtCent,runindex,ptw1*TMath::Sin(n*dPhi1));
       }
       else if(dChrg1<0){
         fHCos1nNegChEtaPosVzNeg->Fill(EvtCent,runindex,ptw1*TMath::Cos(n*dPhi1));
         fHSin1nNegChEtaPosVzNeg->Fill(EvtCent,runindex,ptw1*TMath::Sin(n*dPhi1));
       }

     }
     else if(VtxZ<0 && dEta1<0){
       if(dChrg1>0){
         fHCos1nPosChEtaNegVzNeg->Fill(EvtCent,runindex,ptw1*TMath::Cos(n*dPhi1));
         fHSin1nPosChEtaNegVzNeg->Fill(EvtCent,runindex,ptw1*TMath::Sin(n*dPhi1));
       }
       else if(dChrg1<0){
         fHCos1nNegChEtaNegVzNeg->Fill(EvtCent,runindex,ptw1*TMath::Cos(n*dPhi1));
         fHSin1nNegChEtaNegVzNeg->Fill(EvtCent,runindex,ptw1*TMath::Sin(n*dPhi1));
       }
     }
   }

   //if(fabs(w1NUA)>2.0){
   //  std::cout<<" track = "<<i<<" Ch = "<<dChrg1<<"\tPt = "<<dPt1<<"\tptw = "<<ptw1<<"\tw1NUA = "<<w1NUA<<"\tVtxZ = "<<VtxZ<<std::endl;
   //}


   //2nd track loop:
   for(int j=0; j<iTracks; j++) {

     if(j==i) continue;  //Auto-correlation removed.

     pTrack2    =     fEvent->GetTrack(j);
     if(!pTrack2)                continue;
     dPhi2      =          pTrack2->Phi(); //0-2pi range
     dPt2       =          pTrack2-> Pt();
     dEta2      =          pTrack2->Eta();
     dChrg2     =       pTrack2->Charge();

     if(!pTrack2->InPOISelection())
     continue;
   
     if(bApplyNUAforEP){
       ptBin = fFB_Efficiency_Cent[cIndex]->FindBin(dPt2);
       ptw2   = 1.0/fFB_Efficiency_Cent[cIndex]->GetBinContent(ptBin);

       if(bApplyNUACorr){
         if(dChrg2>0){
           iBinNUA = fHCorrectNUApos->FindBin(VtxZ,dPhi2,dEta2);
           w2NUA   = 1.0/fHCorrectNUApos->GetBinContent(iBinNUA);    
  	   if(w2NUA > 1e5 ) w2NUA = 1.0;
	 }
         else if(dChrg2<0){
           iBinNUA = fHCorrectNUApos->FindBin(VtxZ,dPhi2,dEta2);
           w2NUA   = 1.0/fHCorrectNUApos->GetBinContent(iBinNUA);  
  	   if(w2NUA > 1e5 ) w2NUA = 1.0;  
	 }
       }
       WgtEP = ptw1*ptw2*w1NUA*w2NUA;
     }

     if(dChrg1!=dChrg2){
       fHist_Corr3p_EP_Norm_PN[0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),WgtEP);
       fHist_Corr3p_EP_Norm_PN[1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),WgtEP);
       if(cIndex<6){
         fHist_Corr3p_pTSum_EP_V0A_PN[cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),1.0);
         fHist_Corr3p_pTSum_EP_V0C_PN[cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),1.0);
         fHist_Corr3p_pTDiff_EP_V0A_PN[cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),1.0);
         fHist_Corr3p_pTDiff_EP_V0C_PN[cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),1.0);
         fHist_Corr3p_EtaDiff_EP_V0A_PN[cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),WgtEP);
         fHist_Corr3p_EtaDiff_EP_V0C_PN[cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),WgtEP);
       }
     }
     else if(dChrg1>0 && dChrg2>0){
       fHist_Corr3p_EP_Norm_PP[0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),WgtEP);
       fHist_Corr3p_EP_Norm_PP[1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),WgtEP);
       if(cIndex<6){
         fHist_Corr3p_pTSum_EP_V0A_PP[cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),1.0);
         fHist_Corr3p_pTSum_EP_V0C_PP[cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),1.0);
         fHist_Corr3p_pTDiff_EP_V0A_PP[cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),1.0);
         fHist_Corr3p_pTDiff_EP_V0C_PP[cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),1.0);
         fHist_Corr3p_EtaDiff_EP_V0A_PP[cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),WgtEP);
         fHist_Corr3p_EtaDiff_EP_V0C_PP[cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),WgtEP);
       }
     }
     else if(dChrg1<0 && dChrg2<0){
       fHist_Corr3p_EP_Norm_NN[0]->Fill(EvtCent, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),WgtEP);
       fHist_Corr3p_EP_Norm_NN[1]->Fill(EvtCent, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),WgtEP);
       if(cIndex<6){
         fHist_Corr3p_pTSum_EP_V0A_NN[cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),1.0);
         fHist_Corr3p_pTSum_EP_V0C_NN[cIndex]->Fill((dPt1+dPt2)*0.5, TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),1.0);
         fHist_Corr3p_pTDiff_EP_V0A_NN[cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),1.0);
         fHist_Corr3p_pTDiff_EP_V0C_NN[cIndex]->Fill(TMath::Abs(dPt1-dPt2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),1.0);
         fHist_Corr3p_EtaDiff_EP_V0A_NN[cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0A),WgtEP);
         fHist_Corr3p_EtaDiff_EP_V0C_NN[cIndex]->Fill(TMath::Abs(dEta1-dEta2), TMath::Cos(n*dPhi1 + n*dPhi2 - 2.*n*Psi2V0C),WgtEP);
       }
     }
   }//nested loop ends
 }//------ track loop ends ------


 Double_t QTPCRe = 0.;
 Double_t QTPCIm = 0.;

 QTPCRe = QxTPC[1];
 QTPCIm = QyTPC[1];

 fTPCQnxVsCentRun->Fill(centrCL1,runindex,QxcnCor);
 fTPCQnyVsCentRun->Fill(centrCL1,runindex,QycnCor);


 Double_t Psi2TPC = 0.5*(TMath::ATan2(QTPCIm,QTPCRe));
 if(Psi2TPC < 0.) Psi2TPC += pi;


//V0A-V0C EP resolution:
 fHist_Reso2n_EP_Norm_Det[0]->Fill(EvtCent, TMath::Cos(2.*(Psi2V0A-Psi2V0C)));
//V0A-TPC EP resolution:
 fHist_Reso2n_EP_Norm_Det[1]->Fill(EvtCent, TMath::Cos(2.*(Psi2V0A-Psi2TPC)));
//V0C-TPC EP resolution:
 fHist_Reso2n_EP_Norm_Det[2]->Fill(EvtCent, TMath::Cos(2.*(Psi2V0C-Psi2TPC)));


 fHV0AEventPlaneVsCent->Fill(EvtCent,Psi2V0A);
 fHV0CEventPlaneVsCent->Fill(EvtCent,Psi2V0C);
 fHTPCEventPlaneVsCent->Fill(EvtCent,Psi2TPC);


 Double_t uPRe=0.,uNRe=0.,uPIm=0.,uNIm=0.,uN2Re=0.,uN2Im=0.,uP2Re=0.,uP2Im=0.;
 Double_t uPM=0.,uNM=0.;
 //Double_t VCIm=0., VPRe=0.,VPIm=0.;


 uPM = MPOIpos; uNM = MPOIneg;

 uPRe =  QxPos[0];  uNRe =  QxNeg[0]; 
 uPIm =  QyPos[0];  uNIm =  QyNeg[0];

 uP2Re = QxAutoPos[0]; uP2Im =  QyAutoPos[0];
 uN2Re = QxAutoNeg[0]; uN2Im =  QyAutoNeg[0];

 Double_t TwoQpQnV = 0.,TwoQpQpV=0.,TwoQnQnV=0.;

 if(uPM > 1 && uNM > 1){
   //CME w.r.t V0A EP:
   TwoQpQnV = ((uPRe*uNRe-uPIm*uNIm)*QxanCor + (uPRe*uNIm+uPIm*uNRe)*QyanCor) / (uPM*uNM) ;
   TwoQpQpV = ((uPRe*uPRe-uPIm*uPIm-uP2Re)*QxanCor + (2.*uPRe*uPIm-uP2Im)*QyanCor) / (uPM*(uPM-1.)) ;
   TwoQnQnV = ((uNRe*uNRe-uNIm*uNIm-uN2Re)*QxanCor + (2.*uNRe*uNIm-uN2Im)*QyanCor) / (uNM*(uNM-1.)) ;
   //Fill profiles:
   fHist_Corr3p_SP_Norm_PN[0]->Fill(EvtCent, TwoQpQnV, uPM*uNM);
   fHist_Corr3p_SP_Norm_PP[0]->Fill(EvtCent, TwoQpQpV, (uPM*(uPM-1.)));
   fHist_Corr3p_SP_Norm_NN[0]->Fill(EvtCent, TwoQnQnV, (uNM*(uNM-1.)));

   //-------- CME - ZDN correlators:---------
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  0., TwoQpQnV, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  1., energyZPA, 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  2., energyZPC, 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  4., energyZNA, 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  5., energyZNC, 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  7., TwoQpQnV*energyZPA, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  8., TwoQpQnV*energyZPC, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent,  9., TwoQpQnV*(energyZPC+energyZPA), uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent, 10., TwoQpQnV*energyZNA, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent, 11., TwoQpQnV*energyZNC, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[0]->Fill(EvtCent, 12., TwoQpQnV*(energyZNC+energyZNA), uPM*uNM);

   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  0., TwoQpQpV, (uPM*(uPM-1.)));
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  1., energyZPA, 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  2., energyZPC, 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  4., energyZNA, 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  5., energyZNC, 1.);
   //fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  7., TwoQpQpV*energyZPA, (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  8., TwoQpQpV*energyZPC, (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent,  9., TwoQpQpV*(energyZPC+energyZPA), (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent, 10., TwoQpQpV*energyZNA, (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent, 11., TwoQpQpV*energyZNC, (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[0]->Fill(EvtCent, 12., TwoQpQpV*(energyZNC+energyZNA), (uPM*(uPM-1.)));

   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  0., TwoQnQnV, (uNM*(uNM-1.)));
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  1., energyZPA, 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  2., energyZPC, 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  4., energyZNA, 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  5., energyZNC, 1.);
   //fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  7., TwoQnQnV*energyZPA, (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  8., TwoQnQnV*energyZPC, (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent,  9., TwoQnQnV*(energyZPC+energyZPA), (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent, 10., TwoQnQnV*energyZNA, (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent, 11., TwoQnQnV*energyZNC, (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[0]->Fill(EvtCent, 12., TwoQnQnV*(energyZNC+energyZNA), (uNM*(uNM-1.)));
   //-----------------------------------------




   //CME w.r.t V0C EP:
   TwoQpQnV = 0.; TwoQpQpV = 0.; TwoQnQnV = 0.;
   TwoQpQnV = ((uPRe*uNRe-uPIm*uNIm)*QxcnCor + (uPRe*uNIm+uPIm*uNRe)*QycnCor) / (uPM*uNM) ;
   TwoQpQpV = ((uPRe*uPRe-uPIm*uPIm-uP2Re)*QxcnCor + (2.*uPRe*uPIm-uP2Im)*QycnCor) / (uPM*(uPM-1.)) ;
   TwoQnQnV = ((uNRe*uNRe-uNIm*uNIm-uN2Re)*QxcnCor + (2.*uNRe*uNIm-uN2Im)*QycnCor) / (uNM*(uNM-1.)) ;

   //Fill profiles:
   fHist_Corr3p_SP_Norm_PN[1]->Fill(EvtCent, TwoQpQnV, uPM*uNM);
   fHist_Corr3p_SP_Norm_PP[1]->Fill(EvtCent, TwoQpQpV, (uPM*(uPM-1.)));
   fHist_Corr3p_SP_Norm_NN[1]->Fill(EvtCent, TwoQnQnV, (uNM*(uNM-1.)));
  
   //-------- CME - ZDN correlators:---------
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  0., TwoQpQnV, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  1., energyZPA, 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  2., energyZPC, 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  4., energyZNA, 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  5., energyZNC, 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  7., TwoQpQnV*energyZPA, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  8., TwoQpQnV*energyZPC, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent,  9., TwoQpQnV*(energyZPC+energyZPA), uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent, 10., TwoQpQnV*energyZNA, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent, 11., TwoQpQnV*energyZNC, uPM*uNM);
   fHist_Corr3p_ZDN_SP_PN[1]->Fill(EvtCent, 12., TwoQpQnV*(energyZNC+energyZNA), uPM*uNM);

   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  0., TwoQpQpV, (uPM*(uPM-1.)));
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  1., energyZPA, 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  2., energyZPC, 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  4., energyZNA, 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  5., energyZNC, 1.);
   //fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  7., TwoQpQpV*energyZPA, (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  8., TwoQpQpV*energyZPC, (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent,  9., TwoQpQpV*(energyZPC+energyZPA), (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent, 10., TwoQpQpV*energyZNA, (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent, 11., TwoQpQpV*energyZNC, (uPM*(uPM-1.)));
   fHist_Corr3p_ZDN_SP_PP[1]->Fill(EvtCent, 12., TwoQpQpV*(energyZNC+energyZNA), (uPM*(uPM-1.)));

   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  0., TwoQnQnV, (uNM*(uNM-1.)));
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  1., energyZPA, 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  2., energyZPC, 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  3., (energyZPA+energyZPC), 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  4., energyZNA, 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  5., energyZNC, 1.);
   //fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  6., (energyZNA+energyZNC), 1.);
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  7., TwoQnQnV*energyZPA, (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  8., TwoQnQnV*energyZPC, (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent,  9., TwoQnQnV*(energyZPC+energyZPA), (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent, 10., TwoQnQnV*energyZNA, (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent, 11., TwoQnQnV*energyZNC, (uNM*(uNM-1.)));
   fHist_Corr3p_ZDN_SP_NN[1]->Fill(EvtCent, 12., TwoQnQnV*(energyZNC+energyZNA), (uNM*(uNM-1.)));
   //-----------------------------------------



   //V0A-V0C SP resolution:
   fHist_Reso2n_SP_Norm_Det[0]->Fill(EvtCent, (QxcnCor*QxanCor+QycnCor*QyanCor));
   //V0A-TPC SP resolution:
   fHist_Reso2n_SP_Norm_Det[1]->Fill(EvtCent, (QTPCRe*QxanCor+QTPCIm*QyanCor));
   //V0C-TPC SP resolution:
   fHist_Reso2n_SP_Norm_Det[2]->Fill(EvtCent, (QTPCRe*QxcnCor+QTPCIm*QycnCor));

   //V0A-V0C:
   fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 0., (QxcnCor*QxanCor+QycnCor*QyanCor)); 
   fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 1., (QxcnCor*QxanCor+QycnCor*QyanCor)*energyZPA); 
   fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 2., (QxcnCor*QxanCor+QycnCor*QyanCor)*energyZPC); 
   fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 3., (QxcnCor*QxanCor+QycnCor*QyanCor)*(energyZPA+energyZPC)); 
   fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 4., (QxcnCor*QxanCor+QycnCor*QyanCor)*energyZNA); 
   fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 5., (QxcnCor*QxanCor+QycnCor*QyanCor)*energyZNC); 
   fHist_Reso2n_ZDN_SP_Det[0]->Fill(EvtCent, 6., (QxcnCor*QxanCor+QycnCor*QyanCor)*(energyZNA+energyZNC)); 

   //V0A-TPC 
   fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 0., (QTPCRe*QxanCor+QTPCIm*QyanCor));
   fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 1., (QTPCRe*QxanCor+QTPCIm*QyanCor)*energyZPA);
   fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 2., (QTPCRe*QxanCor+QTPCIm*QyanCor)*energyZPC);
   fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 3., (QTPCRe*QxanCor+QTPCIm*QyanCor)*(energyZPA+energyZPC));
   fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 4., (QTPCRe*QxanCor+QTPCIm*QyanCor)*energyZNA);
   fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 5., (QTPCRe*QxanCor+QTPCIm*QyanCor)*energyZNC);
   fHist_Reso2n_ZDN_SP_Det[1]->Fill(EvtCent, 6., (QTPCRe*QxanCor+QTPCIm*QyanCor)*(energyZNA+energyZNC));

  //V0C-TPC
   fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 0., (QTPCRe*QxcnCor+QTPCIm*QycnCor));
   fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 1., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*energyZPA);
   fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 2., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*energyZPC);
   fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 3., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*(energyZPA+energyZPC));
   fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 1., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*energyZNA);
   fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 2., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*energyZNC);
   fHist_Reso2n_ZDN_SP_Det[2]->Fill(EvtCent, 3., (QTPCRe*QxcnCor+QTPCIm*QycnCor)*(energyZNA+energyZNC));

 }


 PostData(1,fListHistos);
 PostData(2,fListCalibs); 
 
 fHist_Event_count->Fill(9.5);

 //if(fievent%20==0) {
   //std::cout<<"irun = "<<runindex<<" n "<<n<<" cent= "<<EvtCent<<"\tsumMa= "<<sumMa<<"\tQxA = "<<QxanCor<<std::endl;
   //std::cout<<" cent= "<<EvtCent<<"\teZNC= "<<energyZNC<<"\teZPC = "<<energyZPC<<"\teZNA= "<<energyZNA<<"\teZPA = "<<energyZPA<<std::endl;
   //}

 fievent++;

}
//======================= UserExec done =========================



void AliAnalysisTaskCMEV0::Terminate(Option_t *)
{
 AliDebug(2,"\n ... AliAnalysisTaskCMEV0::Terminate() is being called ...  \n");
}













void AliAnalysisTaskCMEV0::GetV0QvectAndMult(const AliAODVZERO *aodV0, Double_t& Qxan,Double_t& Qyan,Double_t& sumMa,Double_t& Qxcn,Double_t& Qycn,Double_t& sumMc) 
{
  for(Int_t iV0 = 0; iV0 < 64; iV0++) {
  /*if(fRemChV0A){
     if(iV0 == 46)
        continue;
    }*/
    Double_t phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);
    Float_t multv0 = aodV0->GetMultiplicity(iV0);

    if(iV0 < 32) {
      Double_t multCorC = -10;

      if(iV0 < 8)
         multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(1);
      else if(iV0 >= 8 && iV0 < 16)
         multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(9);
      else if(iV0 >= 16 && iV0 < 24)
         multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(17);
      else if(iV0 >= 24 && iV0 < 32)
         multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(25);

      if(multCorC < 0){
         cout<<"Problem with multiplicity in V0C"<<endl;
         continue;
      }
      Qxcn += TMath::Cos(fHarmonic*phiV0) * multCorC;
      Qycn += TMath::Sin(fHarmonic*phiV0) * multCorC;

      sumMc = sumMc + multCorC;
    } 
    else{
      Double_t multCorA = -10;

      if(iV0 >= 32 && iV0 < 40)
         multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(33);
      else if(iV0 >= 40 && iV0 < 48)
         multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(41);
      else if(iV0 >= 48 && iV0 < 56)
         multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(49);
      else if(iV0 >= 56 && iV0 < 64)
         multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(57);

      if(multCorA < 0){
         cout<<"Problem with multiplicity in V0A"<<endl;
         continue;
      }
      Qxan += TMath::Cos(fHarmonic*phiV0) * multCorA;
      Qyan += TMath::Sin(fHarmonic*phiV0) * multCorA;

      sumMa = sumMa + multCorA;
    }
  }
}//-------- Get V0 QVect and Multiplicity -----

void AliAnalysisTaskCMEV0::GetZDCCorrectionHist(Int_t run)
{
  if(fListZDNCorr){
    fHCorrectZDNP = (TH2D *) fListZDNCorr->FindObject(Form("fHist_ZDC_All_Wgt_VsCent_Run%d",run));
  }
  else{
    fHCorrectZDNP = NULL;
    printf("\n\n ********** ZDC Wgt Histograms NotFound ***************\n\n");
    exit(1);
  }
}

void AliAnalysisTaskCMEV0::GetNUACorrectionHist(Int_t run, Float_t cent)
{
  /*if(fListNUACorr){
    fHCorrectNUApos = (TH3F *) fListNUACorr->FindObject(Form("fHist_WgtPos_EtaPhiVz_Run%d",run));
    fHCorrectNUAneg = (TH3F *) fListNUACorr->FindObject(Form("fHist_WgtNeg_EtaPhiVz_Run%d",run));
  }
  else{
    printf("\n\n ********** NUA Correction Histograms NotFound ***************\n\n");
    exit(1);
  }
  */

  Int_t centBin = -1;

  if(cent<5.0){ centBin = 0;}
  else if(cent>=5.0 && cent<10.0){centBin = 1;}
  else if(cent>=10.0 && cent<20.0){centBin = 2;}
  else if(cent>=20.0) {centBin = 3;}

  /*
  if(fListNUACorr){
    fHCorrectNUApos = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Pos_Cent%d_Run%d",centBin,run));
    fHCorrectNUAneg = (TH3D *) fListNUACorr->FindObject(Form("fHist_NUA_VzPhiEta_Neg_Cent%d_Run%d",centBin,run));
  }
  else{
    printf("\n\n ********** NUA Correction File/mList Not Found EXIT..!! ***************\n\n");
    fHCorrectNUApos = NULL;
    fHCorrectNUAneg = NULL;
    exit(1);
  }
  */

  if(!gGrid){
    TGrid::Connect("alien://");
  }

 //from Rihan's file:
  /*
  TFile *fNUApn = NULL;

  fNUApn = TFile::Open("alien:///alice/cern.ch/user/m/mhaque/calib_files/Run2015o_NUA_PosNeg_RunbyRun_Oct5.root");
  TList  *mListPN = dynamic_cast<TList*> (fNUApn->FindObjectAny("fListTrkRecnter"));

  if(mListPN){
    fHCorrectNUApos = (TH3F *) mListPN->FindObject(Form("fHist_NUA_VzPhiEta_Pos_Cent%d_Run%d",centBin,run));
    fHCorrectNUAneg = (TH3F *) mListPN->FindObject(Form("fHist_NUA_VzPhiEta_Neg_Cent%d_Run%d",centBin,run));
  }
  else{
    printf("\n\n ********** NUA Correction File/mList Not Found EXIT(1) ***************\n\n");
    fHCorrectNUApos = NULL;
    fHCorrectNUAneg = NULL;
    exit(1);
  }
  */

  //from Jacopo's NUA file:
  if(!fileNUApos){
    fileNUApos = TFile::Open("alien:///alice/cern.ch/user/m/mhaque/calib_files/15oHI_FB768_PosCh_CenPhiEtaWeights_VtxRbR.root");
    mListNUAPos = dynamic_cast<TList*> (fileNUApos->FindObjectAny("CenPhiEta Weights"));
  }
  if(!fileNUAneg){
    fileNUAneg = TFile::Open("alien:///alice/cern.ch/user/m/mhaque/calib_files/15oHI_FB768_NegCh_CenPhiEtaWeights_VtxRbR.root");
    mListNUANeg = dynamic_cast<TList*> (fileNUAneg->FindObjectAny("CenPhiEta Weights"));
  }

  if(mListNUAPos){
    fHCorrectNUApos = (TH3D *) mListNUAPos->FindObject(Form("CRCQVecPhiHistVtx[%d][%d]",centBin,run));
  }
  if(mListNUANeg){
    fHCorrectNUAneg = (TH3D *) mListNUANeg->FindObject(Form("CRCQVecPhiHistVtx[%d][%d]",centBin,run));
  }
  //here:
  if(!fHCorrectNUApos || !fHCorrectNUAneg){
    printf("\n\n ******** could not open NUA Histograms for run %d, Use Wgt = 1.0 *********\n\n",run);
    fHCorrectNUApos = new TH3D("fHCorrectNUApos","",1,-10,10,1,0,6.2832,1,-0.8,0.8); 
    fHCorrectNUAneg = new TH3D("fHCorrectNUAneg","",1,-10,10,1,0,6.2832,1,-0.8,0.8); 
    fHCorrectNUApos->SetBinContent(1,1,1,1.0);
    fHCorrectNUAneg->SetBinContent(1,1,1,1.0);
    //exit(1);
  }
}



void AliAnalysisTaskCMEV0::OpenInfoCalbration(Int_t run)
{
  if(!gGrid){
    TGrid::Connect("alien://");
  }

  TFile* foadb = 0;
  //if(!fRemChV0A)
     foadb = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibV0HIR.root");
  //else
   //foadb = TFile::Open("alien:///alice/cern.ch/user/a/adobrin/calibV0HIRNoCh46V0A.root");

    if(!foadb){
        printf("OADB V0 calibration file cannot be opened\n");
        return;
    }

    AliOADBContainer* cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorPfpx");
    if(!cont){
        printf("OADB object hMultV0BefCorr is not available in the file\n");
        return;
    }
    if(!(cont->GetObject(run))){
        printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
        return;
    }
    fMultV0 = ((TH1D*) cont->GetObject(run));

    AliOADBContainer* contQxnam = 0;
    if (fHarmonic == 2.)
        contQxnam = (AliOADBContainer*) foadb->Get("fqxa2m");
    else
        contQxnam = (AliOADBContainer*) foadb->Get("fqxa3m");

    if(!contQxnam){
        printf("OADB object fqxanm is not available in the file\n");
        return;
    }
    if(!(contQxnam->GetObject(run))){
        printf("OADB object fqxanm is not available for run %i\n", run);
        return;
    }
    fQxnmV0A = ((TH1D*) contQxnam->GetObject(run));

    AliOADBContainer* contQynam = 0;
    if (fHarmonic == 2.)
        contQynam = (AliOADBContainer*) foadb->Get("fqya2m");
    else if (fHarmonic == 3.)
        contQynam = (AliOADBContainer*) foadb->Get("fqya3m");
    else if (fHarmonic == 4.)
        contQynam = (AliOADBContainer*) foadb->Get("fqya4m");

    if(!contQynam){
        printf("OADB object fqyanm is not available in the file\n");
        return;
    }
    if(!(contQynam->GetObject(run))){
        printf("OADB object fqyanm is not available for run %i\n", run);
        return;
    }
    fQynmV0A = ((TH1D*) contQynam->GetObject(run));

    AliOADBContainer* contQxnas = 0;
    if (fHarmonic == 2.)
        contQxnas = (AliOADBContainer*) foadb->Get("fqxa2s");
    else
        contQxnas = (AliOADBContainer*) foadb->Get("fqxa3s");

    if(!contQxnas){
        printf("OADB object fqxans is not available in the file\n");
        return;
    }
    if(!(contQxnas->GetObject(run))){
        printf("OADB object fqxans is not available for run %i\n", run);
        return;
    }
    fQxnsV0A = ((TH1D*) contQxnas->GetObject(run));

    AliOADBContainer* contQynas = 0;
    if (fHarmonic == 2.)
        contQynas = (AliOADBContainer*) foadb->Get("fqya2s");
    else if (fHarmonic == 3.)
        contQynas = (AliOADBContainer*) foadb->Get("fqya3s");
    else if (fHarmonic == 4.)
        contQynas = (AliOADBContainer*) foadb->Get("fqya4s");

    if(!contQynas){
        printf("OADB object fqyans is not available in the file\n");
        return;
    }
    if(!(contQynas->GetObject(run))){
        printf("OADB object fqyans is not available for run %i\n", run);
        return;
    }
    fQynsV0A = ((TH1D*) contQynas->GetObject(run));



    AliOADBContainer* contQxncm = 0;
    if (fHarmonic == 2.)
        contQxncm = (AliOADBContainer*) foadb->Get("fqxc2m");
    else
        contQxncm = (AliOADBContainer*) foadb->Get("fqxc3m");

    if(!contQxncm){
        printf("OADB object fqxcnm is not available in the file\n");
        return;
    }
    if(!(contQxncm->GetObject(run))){
        printf("OADB object fqxcnm is not available for run %i\n", run);
        return;
    }
    fQxnmV0C = ((TH1D*) contQxncm->GetObject(run));



    AliOADBContainer* contQyncm = 0;
    if (fHarmonic == 2.)
        contQyncm = (AliOADBContainer*) foadb->Get("fqyc2m");
    else if (fHarmonic == 3.)
        contQyncm = (AliOADBContainer*) foadb->Get("fqyc3m");
    else if (fHarmonic == 4.)
        contQyncm = (AliOADBContainer*) foadb->Get("fqyc4m");

    if(!contQyncm){
        printf("OADB object fqyc2m is not available in the file\n");
        return;
    }
    if(!(contQyncm->GetObject(run))){
        printf("OADB object fqyc2m is not available for run %i\n", run);
        return;
    }
    fQynmV0C = ((TH1D*) contQyncm->GetObject(run));


    AliOADBContainer* contQxncs = 0;
    if (fHarmonic == 2.)
        contQxncs = (AliOADBContainer*) foadb->Get("fqxc2s");
    else
        contQxncs = (AliOADBContainer*) foadb->Get("fqxc3s");

    if(!contQxncs){
        printf("OADB object fqxc2s is not available in the file\n");
        return;
    }
    if(!(contQxncs->GetObject(run))){
        printf("OADB object fqxc2s is not available for run %i\n", run);
        return;
    }
    fQxnsV0C = ((TH1D*) contQxncs->GetObject(run));


    AliOADBContainer* contQyncs = 0;
    if (fHarmonic == 2.)
        contQyncs = (AliOADBContainer*) foadb->Get("fqyc2s");
    else if (fHarmonic == 3.)
        contQyncs = (AliOADBContainer*) foadb->Get("fqyc3s");
    else if (fHarmonic == 4.)
        contQyncs = (AliOADBContainer*) foadb->Get("fqyc4s");

    if(!contQyncs){
        printf("OADB object fqycnm is not available in the file\n");
        return;
    }
    if(!(contQyncs->GetObject(run))){
        printf("OADB object fqycns is not available for run %i\n", run);
        return;
    }
    fQynsV0C = ((TH1D*) contQyncs->GetObject(run));
}//------- OADB container ------




double AliAnalysisTaskCMEV0::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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

Bool_t AliAnalysisTaskCMEV0::PileUpMultiVertex(const AliAODEvent* faod)
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






Bool_t AliAnalysisTaskCMEV0::CheckEventIsPileUp(AliAODEvent *faod) {

 Bool_t BisPileup=kFALSE;

 Double_t centrV0M=300;
 Double_t centrCL1=300;
 Double_t centrCL0=300;
 Double_t centrTRK=300;

 if(sDataSet=="2010"||sDataSet=="2011"){
   centrV0M = ((AliVAODHeader*)faod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
   centrCL1 = ((AliVAODHeader*)faod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");
   centrCL0 = ((AliVAODHeader*)faod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL0");
   centrTRK = ((AliVAODHeader*)faod->GetHeader())->GetCentralityP()->GetCentralityPercentile("TRK");
 }
 else{
   fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
   if(!fMultSelection) {
     printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
     exit(1);
   }
   centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
   centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
   centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
   centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");
 }// 2015


 if(fRejectPileUp && InputEvent()) {
  //if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
    if(sDataSet!="2015") {
      if(PileUpMultiVertex(faod)) {
         fPileUpCount->Fill(0.5);
         BisPileup=kTRUE;
      }
      Int_t isPileup = faod->IsPileupFromSPD(3);
      if(isPileup != 0) {
         fPileUpCount->Fill(1.5);
       //BisPileup=kTRUE; //
      }
      if(((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
         fPileUpCount->Fill(2.5);
         BisPileup=kTRUE;
      }
      if(faod->IsIncompleteDAQ())  {
         fPileUpCount->Fill(3.5);
         BisPileup=kTRUE;
      }

    //check vertex consistency
      const AliAODVertex* vtTrc = faod->GetPrimaryVertex();
      const AliAODVertex* vtSPD = faod->GetPrimaryVertexSPD();

      if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
        fPileUpCount->Fill(5.5);
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
         fPileUpCount->Fill(6.5);
         BisPileup=kTRUE;
      }
      if(fAnalysisUtil->IsPileUpEvent(InputEvent())) {
         fPileUpCount->Fill(7.5);
         BisPileup=kTRUE;
      }
    } ////------ dataset 2010,2011 pile up ---------

    else{ //------------ pileup for 2015 data ----------------- 
      if(!fMultSelection->GetThisEventIsNotPileup())
         fPileUpMultSelCount->Fill(0.5);
      if(!fMultSelection->GetThisEventIsNotPileupMV())
         fPileUpMultSelCount->Fill(1.5);
      if(!fMultSelection->GetThisEventIsNotPileupInMultBins())
         fPileUpMultSelCount->Fill(2.5);
      if(!fMultSelection->GetThisEventHasNoInconsistentVertices())
         fPileUpMultSelCount->Fill(3.5);
      if(!fMultSelection->GetThisEventPassesTrackletVsCluster())
         fPileUpMultSelCount->Fill(4.5);
      if(!fMultSelection->GetThisEventIsNotAsymmetricInVZERO())
         fPileUpMultSelCount->Fill(5.5);
      if(!fMultSelection->GetThisEventIsNotIncompleteDAQ())
         fPileUpMultSelCount->Fill(6.5);
      if(!fMultSelection->GetThisEventHasGoodVertex2016())
         fPileUpMultSelCount->Fill(7.5);

         BisPileup=kFALSE;

     //-- pile-up a la Dobrin for LHC15o -----
      if(PileUpMultiVertex(faod)) {
         fPileUpCount->Fill(0.5);
         BisPileup=kTRUE;
      }
      Int_t isPileup = faod->IsPileupFromSPD(3);
      if(isPileup != 0) {
         fPileUpCount->Fill(1.5);
         BisPileup=kTRUE;          
      }
      if(((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
         fPileUpCount->Fill(2.5);
         BisPileup=kTRUE;
      }
      if(faod->IsIncompleteDAQ())  {
         fPileUpCount->Fill(3.5);
         BisPileup=kTRUE;
      }
      if(fabs(centrV0M-centrCL1)>7.5)  {
         fPileUpCount->Fill(4.5);
         BisPileup=kTRUE;
      }

     // check vertex consistency
      const AliAODVertex* vtTrc = faod->GetPrimaryVertex();
      const AliAODVertex* vtSPD = faod->GetPrimaryVertexSPD();

      if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
        fPileUpCount->Fill(5.5);
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
        fPileUpCount->Fill(6.5);
        BisPileup=kTRUE;
      }

      //cuts on tracks
      const Int_t nTracks = faod->GetNumberOfTracks();
      Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks();

      //Int_t multTrk = 0;
      //Int_t multTrkBefC = 0;
      //Int_t multTrkTOFBefC = 0;
      Int_t multTPC = 0;

      for(Int_t it = 0; it < nTracks; it++) {
        AliAODTrack* aodTrk = (AliAODTrack*)faod->GetTrack(it);
        if(!aodTrk) {
          delete aodTrk;
          continue;
        }
       //if(aodTrk->TestFilterBit(32)){
       //   multTrkBefC++;
       //   if(TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
       //     multTrkTOFBefC++;
       //     if((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2) && (aodTrk->Pt() < 20.))
       //       multTrk++;
       //}
        if(aodTrk->TestFilterBit(128))
           multTPC++;
      } // end of for AOD track loop

      Double_t multTPCn      = multTPC;
      Double_t multEsdn      = multEsd;
      Double_t multESDTPCDif = multEsdn - multTPCn*3.38;

      /*if(multESDTPCDif > (fRejectPileUpTight?700.:15000.)) {
        fPileUpCount->Fill(7.5);
        BisPileup=kTRUE;
        }*/
      if(multESDTPCDif > 15000.){
        fPileUpCount->Fill(7.5);
        BisPileup=kTRUE;
      }
      else if(fRejectPileUpTight) {
        if(multESDTPCDif > 700.) {
          fPileUpCount->Fill(8.5);
          BisPileup=kTRUE;
        }
        if(BisPileup==kFALSE) {
          if(!fMultSelection->GetThisEventIsNotPileup()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
          if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
          if(BisPileup) fPileUpCount->Fill(9.5);
        }  
      }
    }
  }

 return BisPileup; 
}//-------pile up function ------




Int_t AliAnalysisTaskCMEV0::GetCurrentRunIndex(Int_t  run) {

 Int_t irun = -1;

 for(int i=0;i<fRunFlag;i++){
  if(run==runNums[i])
   {
    irun = i;
    break;
   }
 }
 if(irun<0) {
   printf("\n ... **WARNING** \n::UserExec() runnumber not listed.\n EXIT..\n");
 }

 return irun;
}//------ GetCurrentRunIndex --------


void AliAnalysisTaskCMEV0::InitializeRunArray(TString sPeriod){

 Int_t runArray_2010[89] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137243, 137236, 137235, 137232, 137231, 137162, 137161};

 Int_t runArray_2011[68] = {167915, 168115, 168460, 169035, 169238, 169859, 170228, 167920, 168310, 168464, 169091, 169411, 169923, 170230, 167985, 168311, 168467, 169094, 169415, 170027, 170268, 167987, 168322, 168511, 169138, 169417, 170081, 170269, 167988, 168325, 168512, 169144, 169835, 170155, 170270, 168069, 168341, 168514, 169145, 169837, 170159, 170306, 168076, 168342, 168777, 169148, 169838, 170163, 170308, 168105, 168361, 168826, 169156, 169846, 170193, 170309, 168107, 168362, 168988, 169160, 169855, 170203, 168108, 168458, 168992, 169167, 169858, 170204};

 Int_t runArray_2015[90] = {246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246676, 246675, 246540, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};

 if(sPeriod=="2010"){
  fRunFlag = 89;
  for(int i=0;i<fRunFlag;i++)
    runNums[i] = runArray_2010[i];
 }
 else if(sPeriod=="2011"){
   fRunFlag = 68;                 //<-- 2011 check
  for(int i=0;i<fRunFlag;i++)
    runNums[i] = runArray_2011[i];
 }
 else if(sPeriod=="2015"){
  fRunFlag = 90;
  for(int i=0;i<fRunFlag;i++)
    runNums[i] = runArray_2015[i];
 }
 else{
   printf("\n\n ***** Run Number not defined for this data set. *******\n\n Please modify code..\n\n");
   exit(1);
 }

}//------- InitializeRunArray ---------



void AliAnalysisTaskCMEV0::DefineHistograms(){

  fHist_Event_count = new TH1F("fHist_Event_count"," ",20,0,20);
  fHist_Event_count->GetXaxis()->SetBinLabel(1,"Called Exec()");
  fHist_Event_count->GetXaxis()->SetBinLabel(2,"AOD Exist");
  fHist_Event_count->GetXaxis()->SetBinLabel(3,"PileUp");
  fHist_Event_count->GetXaxis()->SetBinLabel(4,"iCentSPD<90");
  fHist_Event_count->GetXaxis()->SetBinLabel(5,"V0 Mult>0");
  fHist_Event_count->GetXaxis()->SetBinLabel(6,"..TBA..");

  fHist_Event_count->GetXaxis()->SetBinLabel(10,"Final Event");
  fListHistos->Add(fHist_Event_count);

  fPileUpCount = new TH1F("fPileUpCount", "fPileUpCount", 12, 0., 12.);
  fPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultiplicityComb08");
  fPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>7.5");
  fPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(8,"multESDTPCDif=15000");
  fPileUpCount->GetXaxis()->SetBinLabel(9,"multESDTPCDif=700");
  fPileUpCount->GetXaxis()->SetBinLabel(10,"extraPileUpMultSel");
  fListHistos->Add(fPileUpCount);

  fPileUpMultSelCount = new TH1F("fPileUpMultSelCount", "fPileUpMultSelCount", 10, 0., 10.);
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(1,"IsNotPileup");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(2,"IsNotPileupMV");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(3,"IsNotPileupInMultBins");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(4,"InconsistentVertices");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(5,"TrackletVsCluster");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(6,"AsymmetricInVZERO");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(7,"IncompleteDAQ");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(8,"GoodVertex2016");
  fListHistos->Add(fPileUpMultSelCount);





  //----- CME-ZDN correlator SP method---------
  for(int i=0;i<3;i++){ //ZDN_SP
   //Detector: 0 = V0A, 1 = V0C, 3 = Q-cumulant
    fHist_Corr3p_ZDN_SP_PN[i] = new TProfile2D(Form("fHist_Corr3p_ZDN_SP_PosNeg_Det%d",i+1),"opposit charge correlator",90,0,90,15,0,15,"");
    fHist_Corr3p_ZDN_SP_PN[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_ZDN_SP_PN[i]);
    fHist_Corr3p_ZDN_SP_PP[i] = new TProfile2D(Form("fHist_Corr3p_ZDN_SP_PosPos_Det%d",i+1),"pos-pos charge correlator",90,0,90,15,0,15,"");
    fHist_Corr3p_ZDN_SP_PP[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_ZDN_SP_PP[i]);
    fHist_Corr3p_ZDN_SP_NN[i] = new TProfile2D(Form("fHist_Corr3p_ZDN_SP_NegNeg_Det%d",i+1),"neg-neg charge correlator",90,0,90,15,0,15,"");
    fHist_Corr3p_ZDN_SP_NN[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_ZDN_SP_NN[i]);
  }
  //SP Resolution:
  for(int i=0;i<3;i++){
  //Det: 0 = v0c-v0a, 1 = v0a-TPC, 2 = v0c-TPC, 
    fHist_Reso2n_ZDN_SP_Det[i]  = new TProfile2D(Form("fHist_Reso2n_ZDN_SP_DetComb%d",i+1),"Event plane Resolution",90,0,90,8,0,8,"");
    fHist_Reso2n_ZDN_SP_Det[i]->Sumw2();
    fListHistos->Add(fHist_Reso2n_ZDN_SP_Det[i]);
  }
  //--------------------------------------------

  Double_t centRange[11]    = {0,5,10,20,30,40,50,60,70,80,90};

  //----- CME SP method histograms ---------
  for(int i=0;i<3;i++){
   //Detector: 0 = V0A, 1 = V0C, 3 = Q-cumulant
    fHist_Corr3p_SP_Norm_PN[i] = new TProfile(Form("fHist_Corr3p_SP_Norm_PosNeg_Det%d",i+1),"opposit charge correlator",10,centRange,"");
    fHist_Corr3p_SP_Norm_PN[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_SP_Norm_PN[i]);
    fHist_Corr3p_SP_Norm_PP[i] = new TProfile(Form("fHist_Corr3p_SP_Norm_PosPos_Det%d",i+1),"pos-pos charge correlator",10,centRange,"");
    fHist_Corr3p_SP_Norm_PP[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_SP_Norm_PP[i]);
    fHist_Corr3p_SP_Norm_NN[i] = new TProfile(Form("fHist_Corr3p_SP_Norm_NegNeg_Det%d",i+1),"neg-neg charge correlator",10,centRange,"");
    fHist_Corr3p_SP_Norm_NN[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_SP_Norm_NN[i]);
  }
  //SP Resolution:
  for(int i=0;i<3;i++){
  //Det: 0 = v0c-v0a, 1 = v0a-TPC, 2 = v0c-TPC, 
    fHist_Reso2n_SP_Norm_Det[i]  = new TProfile(Form("fHist_Reso2n_SP_Norm_DetComb%d",i+1),"Event plane Resolution",10,centRange,"");
    fHist_Reso2n_SP_Norm_Det[i]->Sumw2();
    fListHistos->Add(fHist_Reso2n_SP_Norm_Det[i]);
  }

  //----- CME EP method histograms ---------
  for(int i=0;i<3;i++){
   //Detector: 0 = V0A, 1 = V0C, 3 = Q-cumulant
    fHist_Corr3p_EP_Norm_PN[i] = new TProfile(Form("fHist_Corr3p_EP_Norm_PosNeg_Det%d",i+1),"opposit charge correlator",10,centRange,"");
    fHist_Corr3p_EP_Norm_PN[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_EP_Norm_PN[i]);
    fHist_Corr3p_EP_Norm_PP[i] = new TProfile(Form("fHist_Corr3p_EP_Norm_PosPos_Det%d",i+1),"pos-pos charge correlator",10,centRange,"");
    fHist_Corr3p_EP_Norm_PP[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_EP_Norm_PP[i]);
    fHist_Corr3p_EP_Norm_NN[i] = new TProfile(Form("fHist_Corr3p_EP_Norm_NegNeg_Det%d",i+1),"neg-neg charge correlator",10,centRange,"");
    fHist_Corr3p_EP_Norm_NN[i]->Sumw2();
    fListHistos->Add(fHist_Corr3p_EP_Norm_NN[i]);
  }
  //EP Resolution:
  for(int i=0;i<3;i++){
  //Det: 0 = v0c-v0a, 1 = v0a-TPC, 2 = v0c-TPC, 
    fHist_Reso2n_EP_Norm_Det[i]  = new TProfile(Form("fHist_Reso2n_EP_Norm_DetComb%d",i+1),"Event plane Resolution",10,centRange,"");
    fHist_Reso2n_EP_Norm_Det[i]->Sumw2();
    fListHistos->Add(fHist_Reso2n_EP_Norm_Det[i]);
  }

  fHV0AEventPlaneVsCent = new TH2F("fHV0AEventPlaneVsCent","Psi2 from V0A",10,centRange,50,0,3.1415);
  fListHistos->Add(fHV0AEventPlaneVsCent);
  fHV0CEventPlaneVsCent = new TH2F("fHV0CEventPlaneVsCent","Psi2 from V0C",10,centRange,50,0,3.1415);
  fListHistos->Add(fHV0CEventPlaneVsCent);
  fHTPCEventPlaneVsCent = new TH2F("fHTPCEventPlaneVsCent","Psi2 from TPC",10,centRange,50,0,3.1415);
  fListHistos->Add(fHTPCEventPlaneVsCent);


  Char_t name[100], title[100];
 
  //TH1::SetDefaultSumw2();

  //Differential in pT:
  //Double_t pTRange[24] = {0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.,2.33,2.66,3.,3.5,4.,5.,6.,8.,10.};

  Double_t pTRange[21] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0}; 
 
  for(int i=0;i<6;i++){
    sprintf(name,"fHist_Corr3p_pTSum_EP_V0A_PN_Cent%d",i);
    sprintf(title,"PN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTSum_EP_V0A_PN[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTSum_EP_V0A_PN[i]);

    sprintf(name,"fHist_Corr3p_pTSum_EP_V0A_PP_Cent%d",i);
    sprintf(title,"PP 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTSum_EP_V0A_PP[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTSum_EP_V0A_PP[i]);

    sprintf(name,"fHist_Corr3p_pTSum_EP_V0A_NN_Cent%d",i);
    sprintf(title,"NN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTSum_EP_V0A_NN[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTSum_EP_V0A_NN[i]);
    //-----v0c----
    sprintf(name,"fHist_Corr3p_pTSum_EP_V0C_PN_Cent%d",i);
    sprintf(title,"PN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTSum_EP_V0C_PN[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTSum_EP_V0C_PN[i]);

    sprintf(name,"fHist_Corr3p_pTSum_EP_V0C_PP_Cent%d",i);
    sprintf(title,"PP 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTSum_EP_V0C_PP[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTSum_EP_V0C_PP[i]);

    sprintf(name,"fHist_Corr3p_pTSum_EP_V0C_NN_Cent%d",i);
    sprintf(title,"NN 3p vs (pT1+pT2)/2, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTSum_EP_V0C_NN[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTSum_EP_V0C_NN[i]);
  }


  for(int i=0;i<6;i++){
    sprintf(name,"fHist_Corr3p_pTDiff_EP_V0A_PN_Cent%d",i);
    sprintf(title,"PN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTDiff_EP_V0A_PN[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0A_PN[i]);

    sprintf(name,"fHist_Corr3p_pTDiff_EP_V0A_PP_Cent%d",i);
    sprintf(title,"PP 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTDiff_EP_V0A_PP[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0A_PP[i]);

    sprintf(name,"fHist_Corr3p_pTDiff_EP_V0A_NN_Cent%d",i);
    sprintf(title,"NN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTDiff_EP_V0A_NN[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0A_NN[i]);
    //-----v0c----
    sprintf(name,"fHist_Corr3p_pTDiff_EP_V0C_PN_Cent%d",i);
    sprintf(title,"PN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTDiff_EP_V0C_PN[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0C_PN[i]);

    sprintf(name,"fHist_Corr3p_pTDiff_EP_V0C_PP_Cent%d",i);
    sprintf(title,"PP 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTDiff_EP_V0C_PP[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0C_PP[i]);

    sprintf(name,"fHist_Corr3p_pTDiff_EP_V0C_NN_Cent%d",i);
    sprintf(title,"NN 3p vs |pT1-pT2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_pTDiff_EP_V0C_NN[i] = new TProfile(name,title,20,pTRange,"");
    fListHistos->Add(fHist_Corr3p_pTDiff_EP_V0C_NN[i]);
  }

 
  Double_t EtaRange[9] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6};

  for(int i=0;i<6;i++){
    sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0A_PN_Cent%d",i);
    sprintf(title,"PN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_EtaDiff_EP_V0A_PN[i] = new TProfile(name,title,8,EtaRange,"");
    fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0A_PN[i]);

    sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0A_PP_Cent%d",i);
    sprintf(title,"PP 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_EtaDiff_EP_V0A_PP[i] = new TProfile(name,title,8,EtaRange,"");
    fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0A_PP[i]);

    sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0A_NN_Cent%d",i);
    sprintf(title,"NN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_EtaDiff_EP_V0A_NN[i] = new TProfile(name,title,8,EtaRange,"");
    fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0A_NN[i]);
    //-----v0c----
    sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0C_PN_Cent%d",i);
    sprintf(title,"PN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_EtaDiff_EP_V0C_PN[i] = new TProfile(name,title,8,EtaRange,"");
    fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0C_PN[i]);

    sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0C_PP_Cent%d",i);
    sprintf(title,"PP 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_EtaDiff_EP_V0C_PP[i] = new TProfile(name,title,8,EtaRange,"");
    fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0C_PP[i]);

    sprintf(name,"fHist_Corr3p_EtaDiff_EP_V0C_NN_Cent%d",i);
    sprintf(title,"NN 3p vs |Eta1-Eta2|, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHist_Corr3p_EtaDiff_EP_V0C_NN[i] = new TProfile(name,title,8,EtaRange,"");
    fListHistos->Add(fHist_Corr3p_EtaDiff_EP_V0C_NN[i]);
  }





  //=========== Calibtation Histograms etc ==================
  fListCalibs->Add(fHist_Event_count);

  //for debug only, remove after stable code
  hUnderOverBinNUApos = new TH1F("hUnderOverBinNUApos","",90,0,90);
  fListCalibs->Add(hUnderOverBinNUApos);
  hUnderOverBinNUAneg = new TH1F("hUnderOverBinNUAneg","",90,0,90); 
  fListCalibs->Add(hUnderOverBinNUAneg);

  Double_t fCentBinQvect[16] = {0.,2.5,5,10,15,20,25,30,35,40,45,50,60,70,80,90};

  fHEnergyZNCvsCent = new TH2F("fHEnergyZNCvsCent","ZNC Energy vs cent",15,fCentBinQvect,4000,0,200000);
  fListCalibs->Add(fHEnergyZNCvsCent);
  fHEnergyZNAvsCent = new TH2F("fHEnergyZNAvsCent","ZNA Energy vs cent",15,fCentBinQvect,4000,0,200000);
  fListCalibs->Add(fHEnergyZNAvsCent);
  fHEnergyZPCvsCent = new TH2F("fHEnergyZPCvsCent","ZPC Energy vs cent",15,fCentBinQvect,2000,0,50000);
  fListCalibs->Add(fHEnergyZPCvsCent);
  fHEnergyZPAvsCent = new TH2F("fHEnergyZPAvsCent","ZPA Energy vs cent",15,fCentBinQvect,2000,0,50000);
  fListCalibs->Add(fHEnergyZPAvsCent);

  fHEnergyZNCvsCentRun = new TProfile2D("fHEnergyZNCvsCentRun","",90,0,90,90,0,90);
  fListCalibs->Add(fHEnergyZNCvsCentRun);
  fHEnergyZNAvsCentRun = new TProfile2D("fHEnergyZNAvsCentRun","",90,0,90,90,0,90);
  fListCalibs->Add(fHEnergyZNAvsCentRun);
  fHEnergyZPCvsCentRun = new TProfile2D("fHEnergyZPCvsCentRun","",90,0,90,90,0,90);
  fListCalibs->Add(fHEnergyZPCvsCentRun);
  fHEnergyZPAvsCentRun = new TProfile2D("fHEnergyZPAvsCentRun","",90,0,90,90,0,90);
  fListCalibs->Add(fHEnergyZPAvsCentRun);

  fHEnergyZPCvsZPA = new TH2F("fHEnergyZPCvsZPA","ZNC Energy vs cent",500,0,50000,500,0,50000);
  fListCalibs->Add(fHEnergyZPCvsZPA);
  fHEnergyZNCvsZNA = new TH2F("fHEnergyZNCvsZNA","ZNC Energy vs cent",500,50000,150000,500,50000,150000);
  fListCalibs->Add(fHEnergyZNCvsZNA);

  fHCentBinTrkRecenter = new TH1F("fHCentBinTrkRecenter","centrality Binning",15,fCentBinQvect);

  if(bFillAvgTPCQn){
    fHCos1nPosChEtaPosVzPos = new TProfile2D("fHCos1nPosChEtaPosVzPos","Cos(n*phi),Eta>0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fHCos1nPosChEtaPosVzPos->GetXaxis()->SetTitle("Centrality");
    fHCos1nPosChEtaPosVzPos->GetYaxis()->SetTitle("Run");
    fListCalibs->Add(fHCos1nPosChEtaPosVzPos);
    fHCos1nPosChEtaNegVzPos = new TProfile2D("fHCos1nPosChEtaNegVzPos","Cos(n*phi),Eta<0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHCos1nPosChEtaNegVzPos);
    fHCos1nPosChEtaPosVzNeg = new TProfile2D("fHCos1nPosChEtaPosVzNeg","Cos(n*phi),Eta>0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHCos1nPosChEtaPosVzNeg);
    fHCos1nPosChEtaNegVzNeg = new TProfile2D("fHCos1nPosChEtaNegVzNeg","Cos(n*phi),Eta<0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHCos1nPosChEtaNegVzNeg);
    fHSin1nPosChEtaPosVzPos = new TProfile2D("fHSin1nPosChEtaPosVzPos","Sin(n*phi),Eta>0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHSin1nPosChEtaPosVzPos);
    fHSin1nPosChEtaNegVzPos = new TProfile2D("fHSin1nPosChEtaNegVzPos","Sin(n*phi),Eta<0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHSin1nPosChEtaNegVzPos);
    fHSin1nPosChEtaPosVzNeg = new TProfile2D("fHSin1nPosChEtaPosVzNeg","Sin(n*phi),Eta>0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHSin1nPosChEtaPosVzNeg);
    fHSin1nPosChEtaNegVzNeg = new TProfile2D("fHSin1nPosChEtaNegVzNeg","Sin(n*phi),Eta<0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHSin1nPosChEtaNegVzNeg);

    fHCos1nNegChEtaPosVzPos = new TProfile2D("fHCos1nNegChEtaPosVzPos","Cos(n*phi),Eta>0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fHCos1nNegChEtaPosVzPos->GetXaxis()->SetTitle("Centrality");
    fHCos1nNegChEtaPosVzPos->GetYaxis()->SetTitle("Run");
    fListCalibs->Add(fHCos1nNegChEtaPosVzPos);
    fHCos1nNegChEtaNegVzPos = new TProfile2D("fHCos1nNegChEtaNegVzPos","Cos(n*phi),Eta<0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHCos1nNegChEtaNegVzPos);
    fHCos1nNegChEtaPosVzNeg = new TProfile2D("fHCos1nNegChEtaPosVzNeg","Cos(n*phi),Eta>0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHCos1nNegChEtaPosVzNeg);
    fHCos1nNegChEtaNegVzNeg = new TProfile2D("fHCos1nNegChEtaNegVzNeg","Cos(n*phi),Eta<0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHCos1nNegChEtaNegVzNeg);
    fHSin1nNegChEtaPosVzPos = new TProfile2D("fHSin1nNegChEtaPosVzPos","Sin(n*phi),Eta>0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHSin1nNegChEtaPosVzPos);
    fHSin1nNegChEtaNegVzPos = new TProfile2D("fHSin1nNegChEtaNegVzPos","Sin(n*phi),Eta<0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHSin1nNegChEtaNegVzPos);
    fHSin1nNegChEtaPosVzNeg = new TProfile2D("fHSin1nNegChEtaPosVzNeg","Sin(n*phi),Eta>0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHSin1nNegChEtaPosVzNeg);
    fHSin1nNegChEtaNegVzNeg = new TProfile2D("fHSin1nNegChEtaNegVzNeg","Sin(n*phi),Eta<0,Vz>0",90,0,90,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHSin1nNegChEtaNegVzNeg);
  }


  fV0AQnxVsCentRun = new TProfile2D("fV0ACos2nVsCentRun","<Cos2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"s");
  fListCalibs->Add(fV0AQnxVsCentRun);
  fV0AQnyVsCentRun = new TProfile2D("fV0ASin2nVsCentRun","<Sin2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"s");
  fListCalibs->Add(fV0AQnyVsCentRun);

  fV0CQnxVsCentRun = new TProfile2D("fV0CCos2nVsCentRun","<Cos2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"s");
  fListCalibs->Add(fV0CQnxVsCentRun);
  fV0CQnyVsCentRun = new TProfile2D("fV0CSin2nVsCentRun","<Sin2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"s");
  fListCalibs->Add(fV0CQnyVsCentRun);

  fTPCQnxVsCentRun = new TProfile2D("fTPCCos2nVsCentRun","<Cos2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"s");
  fListCalibs->Add(fTPCQnxVsCentRun);
  fTPCQnyVsCentRun = new TProfile2D("fTPCSin2nVsCentRun","<Sin2> vs cent,Run",90,0,90,fRunFlag,0,fRunFlag,"s");
  fListCalibs->Add(fTPCQnyVsCentRun);


  for(int i=0;i<10;i++){
    sprintf(name,"fHistChPosvsEtaPtRun_Cent%d",i);
    sprintf(title,"Pos Ch, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHistChPosvsEtaPtRun[i]  = new TH3F(name,title,25,0.2,5.2,16,-0.8,0.8,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHistChPosvsEtaPtRun[i]);
    sprintf(name,"fHistChNegvsEtaPtRun_Cent%d",i);
    sprintf(title,"Neg Ch, Cent %2.0f-%2.0f",centRange[i],centRange[i+1]);
    fHistChNegvsEtaPtRun[i]  = new TH3F(name,title,25,0.2,5.2,16,-0.8,0.8,fRunFlag,0,fRunFlag);
    fListCalibs->Add(fHistChNegvsEtaPtRun[i]);
  }


  Int_t gCentForNUA[5] = {0,5,10,40,90};

  if(bFillEtaPhiNUA){
   for(int i=0;i<4;i++){
    for(int j=0;j<fRunFlag;j++){
      sprintf(name,"fHistEtaPhiVz_Pos_Cent%d_Run%d",i,runNums[j]);
      sprintf(title,"eta,phi,Vz Pos Cent%d-%d%%",gCentForNUA[i],gCentForNUA[i+1]);
      fHist3DEtaPhiVz_Pos_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.2832,10,-0.8,0.8); 
      fListCalibs->Add(fHist3DEtaPhiVz_Pos_Run[i][j]);

      sprintf(name,"fHistEtaPhiVz_Neg_Cent%d_Run%d",i,runNums[j]);
      sprintf(title,"eta,phi,Vz Pos Cent%d-%d%%",gCentForNUA[i],gCentForNUA[i+1]);
      fHist3DEtaPhiVz_Neg_Run[i][j] = new TH3F(name,title,10,-10,10,50,0,6.2832,10,-0.8,0.8); 
      fListCalibs->Add(fHist3DEtaPhiVz_Neg_Run[i][j]);
    }
   }
  }



}


