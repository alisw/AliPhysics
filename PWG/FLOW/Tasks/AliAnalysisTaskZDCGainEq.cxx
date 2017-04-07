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
// 
// Author: Rihan Haque (mhaque@cern.ch)
///////////////////////////////////////////////


#include "Riostream.h" //needed as include

class TFile;
class TList;
class AliAnalysisTaskSE;

#include "TProfile.h"  //needed as include
#include "TProfile2D.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "AliMultSelection.h"
#include "AliVVertex.h"

#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "AliAnalysisTaskZDCGainEq.h"

using std::endl;
using std::cout;


ClassImp(AliAnalysisTaskZDCGainEq)

//________________________________________________________________________
AliAnalysisTaskZDCGainEq::AliAnalysisTaskZDCGainEq(const char *name) :
  AliAnalysisTaskSE(name),
  fEvent(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHistos(NULL),
  fListZDCQxy(NULL),
  fListZDCWgt(NULL),
  fListDummy1(NULL),
  fListHijing(NULL),
  fListSubRun(NULL),
  fRejectPileUpTight(kFALSE),
  fRejectPileUp(kFALSE),
  bFillCosSin(kFALSE),
  bFillZDCQAon(kFALSE),
  bRunAveragedQn(kFALSE),
  fHarmonic(2),
  frunflag(0),
  fievent(0),
  vxBin(0),
  vyBin(0),
  vzBin(0),
  fcheckOnce(0),
  fOldRunNum(0),
  fHist_Event_count(NULL),
  fPileUpMultSelCount(NULL),
  fPileUpCount(NULL),
  fHist_ChanWgt_ZDCC(NULL),
  fHist_ChanWgt_ZDCA(NULL),
  fHist_Vx_ArrayFinder(NULL),
  fHist_Vy_ArrayFinder(NULL),
  fHist_Vz_ArrayFinder(NULL),
  fHist_Task_config(NULL),
  fHist_Cent_woZDCcut(NULL),
  fHist_Cent_wiZDCcut(NULL),
  fHist_CutParameters(NULL),
  fHist_Psi1_ZDCC_wGainCorr(NULL),
  fHist_Psi1_ZDCA_wGainCorr(NULL)
{
  for(int i=0;i<90;i++){
   runNums[i] = 0;
    for(int j=0;j<10;j++){
     fHist_znCx_V0_VxVy[i][j] = NULL;
     fHist_znCy_V0_VxVy[i][j] = NULL;
     fHist_znAx_V0_VxVy[i][j] = NULL;
     fHist_znAy_V0_VxVy[i][j] = NULL;
    }     
  }

  for(int i=0;i<4;i++){
    for(int j=0;j<5;j++){
     fHist_Qx_vs_Obs_woCorr[i][j] = NULL;
     fHist_XX_vs_Obs_woCorr[i][j] = NULL;
    }
  }

  for(int i=0;i<2;i++){
    VxCut[i] = 0;
    VyCut[i] = 0;
    VzCut[i] = 0;
  }


  DefineInput(1, AliFlowEventSimple::Class()); // Input slot #1 works with an AliFlowEventSimple
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());

  fDataSet="2010";
  fAnalysisSet="DoGainEq";
  sCentEstimator="V0";

 //fTotalQvector = new TString("QaQb");         // "QaQb" (means Qa+Qb), "Qa"  or "Qb"

}//-------------constructor-----------------

//________________________________________________
AliAnalysisTaskZDCGainEq::AliAnalysisTaskZDCGainEq() :
  AliAnalysisTaskSE(),
  fEvent(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHistos(NULL),
  fListZDCQxy(NULL),
  fListZDCWgt(NULL),
  fListDummy1(NULL),
  fListHijing(NULL),
  fListSubRun(NULL),
  fRejectPileUpTight(kFALSE),
  fRejectPileUp(kFALSE),
  bFillCosSin(kFALSE),
  bFillZDCQAon(kFALSE),
  bRunAveragedQn(kFALSE),
  fHarmonic(2),
  frunflag(0),
  fievent(0),
  vxBin(0),
  vyBin(0),
  vzBin(0),
  fcheckOnce(0),
  fOldRunNum(0),
  fHist_Event_count(NULL),
  fPileUpMultSelCount(NULL),
  fPileUpCount(NULL),
  fHist_ChanWgt_ZDCC(NULL),
  fHist_ChanWgt_ZDCA(NULL),
  fHist_Vx_ArrayFinder(NULL),
  fHist_Vy_ArrayFinder(NULL),
  fHist_Vz_ArrayFinder(NULL),
  fHist_Task_config(NULL),
  fHist_Cent_woZDCcut(NULL),
  fHist_Cent_wiZDCcut(NULL),
  fHist_CutParameters(NULL),
  fHist_Psi1_ZDCC_wGainCorr(NULL),
  fHist_Psi1_ZDCA_wGainCorr(NULL)
{
  for(int i=0;i<90;i++){
   runNums[i] = 0;
    for(int j=0;j<10;j++){
     fHist_znCx_V0_VxVy[i][j] = NULL;
     fHist_znCy_V0_VxVy[i][j] = NULL;
     fHist_znAx_V0_VxVy[i][j] = NULL;
     fHist_znAy_V0_VxVy[i][j] = NULL;
    }     
  }

  for(int i=0;i<4;i++){
    for(int j=0;j<5;j++){
     fHist_Qx_vs_Obs_woCorr[i][j] = NULL;
     fHist_XX_vs_Obs_woCorr[i][j] = NULL;
    }
  }

  for(int i=0;i<2;i++){
    VxCut[i] = 0;
    VyCut[i] = 0;
    VzCut[i] = 0;
  }

  fDataSet="2010";
  fAnalysisSet="DoGainEq";
  sCentEstimator="V0";
}





//________________________________________________________________________
AliAnalysisTaskZDCGainEq::~AliAnalysisTaskZDCGainEq()
{
  delete                 fListHistos;        
  delete                 fListZDCQxy;         
  delete                 fListZDCWgt;         
  delete                 fListDummy1;        
  delete                 fListHijing;        

  delete              fMultSelection; 
  delete               fAnalysisUtil; // it is '= new' !!!


 //these histograms are not in any list:
 //therefore, deleted manually.
 //delete         *fHist_ChanWgt_ZDCC;  //can't delete, Execption thrown.!! why?
 //delete         *fHist_ChanWgt_ZDCA;  
 //delete       *fHist_Vx_ArrayFinder; 
 //delete       *fHist_Vy_ArrayFinder; 
 //delete       *fHist_Vz_ArrayFinder; 


 //printf("\n\n ~AliAnalysisTaskZDCGainEq::Destructor is called..!!\n\n");
}

//________________________________________________________________________
void AliAnalysisTaskZDCGainEq::UserCreateOutputObjects()
{

 Int_t runArray_2010[89] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137243, 137236, 137235, 137232, 137231, 137162, 137161};

 Int_t runArray_2011[68] = {167915, 168115, 168460, 169035, 169238, 169859, 170228, 167920, 168310, 168464, 169091, 169411, 169923, 170230, 167985, 168311, 168467, 169094, 169415, 170027, 170268, 167987, 168322, 168511, 169138, 169417, 170081, 170269, 167988, 168325, 168512, 169144, 169835, 170155, 170270, 168069, 168341, 168514, 169145, 169837, 170159, 170306, 168076, 168342, 168777, 169148, 169838, 170163, 170308, 168105, 168361, 168826, 169156, 169846, 170193, 170309, 168107, 168362, 168988, 169160, 169855, 170203, 168108, 168458, 168992, 169167, 169858, 170204};

 Int_t runArray_2015[90] = {246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246676, 246675, 246540, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};


 if(fDataSet=="2010"){
  frunflag = 89;
  for(int i=0;i<frunflag;i++)
    runNums[i] = runArray_2010[i];
 }
 if(fDataSet=="2011"){
   frunflag = 10;                 //<--------- 2011 runnumbers have to be checked and put...
  for(int i=0;i<frunflag;i++)
    runNums[i] = runArray_2011[i];
 }
 if(fDataSet=="2015"){
  frunflag = 90;
  for(int i=0;i<frunflag;i++)
    runNums[i] = runArray_2015[i];
 }

  fListHistos = new TList();
  fListHistos->SetOwner(kTRUE);

  fHist_Event_count    = new TH1F("fHist_Event_count"," ",20,0,20);
  fListHistos->Add(fHist_Event_count);

  fPileUpCount = new TH1F("fPileUpCount", "fPileUpCount", 9, 0., 9.);
  fPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultiplicityComb08");
  fPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>7.5");
  fPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(8,"multESDTPCDif");
  fPileUpCount->GetXaxis()->SetBinLabel(9,"extraPileUpMultSel");
  fListHistos->Add(fPileUpCount);

  fPileUpMultSelCount = new TH1F("fPileUpMultSelCount", "fPileUpMultSelCount", 8, 0., 8.);
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(1,"IsNotPileup");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(2,"IsNotPileupMV");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(3,"IsNotPileupInMultBins");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(4,"InconsistentVertices");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(5,"TrackletVsCluster");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(6,"AsymmetricInVZERO");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(7,"IncompleteDAQ");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(8,"GoodVertex2016");
  fListHistos->Add(fPileUpMultSelCount);

  fHist_Task_config = new TH1F("fTask_Configuration", "Task Configuration", 20, 0., 20.);
  fHist_Task_config->GetXaxis()->SetBinLabel(1,"IsZDCQAon");
  fHist_Task_config->GetXaxis()->SetBinLabel(2,"IsCentV0M");
  fHist_Task_config->GetXaxis()->SetBinLabel(3,"IsCentTPC");
  fHist_Task_config->GetXaxis()->SetBinLabel(4,"IsCentCL1");
  fHist_Task_config->GetXaxis()->SetBinLabel(5,"IsPileUpOn");
  fHist_Task_config->GetXaxis()->SetBinLabel(6,"IsPileUpTightOn");
  fHist_Task_config->GetXaxis()->SetBinLabel(7,"IsFillGainEq");
  fHist_Task_config->GetXaxis()->SetBinLabel(8,"IsDoGainEq");
  fHist_Task_config->GetXaxis()->SetBinLabel(9,"IsAvailGainFile");
  fHist_Task_config->GetXaxis()->SetBinLabel(10,"IsFillCosSin");
  fHist_Task_config->GetXaxis()->SetBinLabel(11,"IsDoRecenter");
  fHist_Task_config->GetXaxis()->SetBinLabel(12,"IsAvailRecntFile");
  fHist_Task_config->GetXaxis()->SetBinLabel(13,"IsRunByRun");
  fListHistos->Add(fHist_Task_config);


  fHist_CutParameters = new TH1F("fTask_CutParameters", "Variable Cut Values", 10, 0., 10.);
  fHist_CutParameters->GetXaxis()->SetBinLabel(1,"VzCutLowValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(2,"VzCutHighValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(3,"NBins_Vz");
  fHist_CutParameters->GetXaxis()->SetBinLabel(4,"VxCutHighValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(5,"VxCutLowValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(6,"NBins_Vx");
  fHist_CutParameters->GetXaxis()->SetBinLabel(7,"VyCutHighValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(8,"VyCutLowValue");
  fHist_CutParameters->GetXaxis()->SetBinLabel(9,"NBins_Vy");
  fListHistos->Add(fHist_CutParameters);

  fHist_Cent_woZDCcut =  new TH1F("fHist_Cent_before_ZDCcut"," ",100,0,100);
  fListHistos->Add(fHist_Cent_woZDCcut);
  fHist_Cent_wiZDCcut =  new TH1F("fHist_Cent_afterr_ZDCcut"," ",100,0,100);
  fListHistos->Add(fHist_Cent_wiZDCcut);



  Double_t centRange[12]    = {0,5,10,20,30,40,50,60,70,80,90,100};


 //-------- define and add the recenter histograms ----------------
  if(fDataSet=="2010") {
    //vxBin = 9; 
    //vyBin = 9;
    //vzBin = 10;   

    VxCut[0] =  -0.030;
    VxCut[1] =   0.015;
    VyCut[0] =   0.150;
    VyCut[1] =   0.204;
    VzCut[0] =   -10.0;
    VzCut[1] =     9.4;
  }

  if(fDataSet=="2015"){
    //vxBin = 8;
    //vyBin = 8;
    //vzBin = 10;  

    VxCut[0] =   0.060;
    VxCut[1] =   0.086;
    VyCut[0] =   0.321;
    VyCut[1] =   0.345;
    VzCut[0] =   -10.0;
    VzCut[1] =    10.0;
  }

  if(fDataSet=="2011"){
     AliDebug(2,"\n\n !!** WARNING ***!!  \n cuts not defined for DATASET: %s\n ...EXIT...\n\n)");
     exit(1);
  }



  if(bFillZDCQAon){

   fHist_Task_config->Fill(0.5);

   fHist_Psi1_ZDCC_wGainCorr      = new TH1F("fHist_Psi1_ZDCC_wGainCorr","", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_Psi1_ZDCC_wGainCorr);
   fHist_Psi1_ZDCA_wGainCorr      = new TH1F("fHist_Psi1_ZDCA_wGainCorr","", 200, 0,2.*TMath::Pi());
   fListHistos->Add(fHist_Psi1_ZDCA_wGainCorr);

   TString     sNameQn[4] = {"Xa","Xc","Ya","Yc"};
   TString    sNameQn2[4] = {"XaXc","YaYc","XcYa","YcXa"};
   TString     sNameVs[5] = {"Cent","Mult","Vx","Vy","Vz"};      
   Int_t     nBinNumVs[5] = {100,  100,       50,       50, 400}; // number of bins for "Cent", "Mult", "Vx","Vy","Vz"
   Double_t  lBinLowVs[5] = {  0,    0, VxCut[0], VyCut[0], -10}; // bin low  value: "Cent", "Mult", "Vx", "Vy", "Vz"
   Double_t lBinHighVs[5] = {100, 4000, VxCut[1], VyCut[1],  10}; // bin high value: "Cent", "Mult", "Vx", "Vy", "Vz"

   char name[100];

   for(int i=0;i<4;i++) {
     for(int j=0;j<5;j++) {//fHist_Qx_vs_Obs_woCorr
      //store: X,Y position for recenter:
       sprintf(name,"fHist_%s_vs_%s_woCorr",static_cast<const char*>(sNameQn[i]),static_cast<const char*>(sNameVs[j]));
       fHist_Qx_vs_Obs_woCorr[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
       fListHistos->Add(fHist_Qx_vs_Obs_woCorr[i][j]);

       sprintf(name,"fHist_%s_vs_%s_woCorr",static_cast<const char*>(sNameQn2[i]),static_cast<const char*>(sNameVs[j]));
       fHist_XX_vs_Obs_woCorr[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
       fListHistos->Add(fHist_XX_vs_Obs_woCorr[i][j]);

       //sprintf(name,"fHist_%s_vs_%s_withCor",static_cast<const char*>(sNameQn[i]),static_cast<const char*>(sNameVs[j]));
       //fHist_Qx_vs_Obs_withCor[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
       //fListHistos->Add(fHist_Qx_vs_Obs_withCor[i][j]);

       //sprintf(name,"fHist_%s_vs_%s_withCor",static_cast<const char*>(sNameQn2[i]),static_cast<const char*>(sNameVs[j]));
       //fHist_XX_vs_Obs_withCor[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
       //fListHistos->Add(fHist_XX_vs_Obs_withCor[i][j]);
     }
   }
  }


  const int NbinVt =  (vyBin-1)*vxBin + vxBin;

  fHist_Vx_ArrayFinder = new TH1F("fHist_Vx_ArrayFinder","",vxBin,VxCut[0],VxCut[1]);
  fHist_Vy_ArrayFinder = new TH1F("fHist_Vy_ArrayFinder","",vyBin,VyCut[0],VyCut[1]);
  fHist_Vz_ArrayFinder = new TH1F("fHist_Vz_ArrayFinder","",vzBin,VzCut[0],VzCut[1]);


  if(vzBin>10){
    printf("\n\n::UserCreateOutPutObject()\n Vz Binning more than 10 not allowed\n Exit \n\n");
    exit(1);
  }

  const int VzBinIter = (const int) vzBin;


  //filling the Q vectors for recenter:
  if(fAnalysisSet=="DoGainEq") {

    fListZDCQxy = new TList();
    fListZDCQxy->SetOwner(kTRUE);

    if(!bRunAveragedQn) {
     for(int i=0;i<frunflag;i++) {
      for(int j=0;j<VzBinIter;j++) {
         fHist_znCx_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znCx_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90);
         fHist_znCy_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znCy_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90);
         fHist_znAx_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znAx_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90);
         fHist_znAy_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znAy_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90);

         fListZDCQxy->Add(fHist_znCx_V0_VxVy[i][j]);
         fListZDCQxy->Add(fHist_znCy_V0_VxVy[i][j]);
         fListZDCQxy->Add(fHist_znAx_V0_VxVy[i][j]);
         fListZDCQxy->Add(fHist_znAy_V0_VxVy[i][j]);
      }
     }
    }
    else{
     for(int j=0;j<VzBinIter;j++) {
         fHist_znCx_V0_VxVy[0][j] = new TProfile2D(Form("fHist_znCx_V0_Run%d_Vz%d",0,j+1),"",NbinVt,0,NbinVt,90,0,90);
         fHist_znCy_V0_VxVy[0][j] = new TProfile2D(Form("fHist_znCy_V0_Run%d_Vz%d",0,j+1),"",NbinVt,0,NbinVt,90,0,90);
         fHist_znAx_V0_VxVy[0][j] = new TProfile2D(Form("fHist_znAx_V0_Run%d_Vz%d",0,j+1),"",NbinVt,0,NbinVt,90,0,90);
         fHist_znAy_V0_VxVy[0][j] = new TProfile2D(Form("fHist_znAy_V0_Run%d_Vz%d",0,j+1),"",NbinVt,0,NbinVt,90,0,90);

         fListZDCQxy->Add(fHist_znCx_V0_VxVy[0][j]);
         fListZDCQxy->Add(fHist_znCy_V0_VxVy[0][j]);
         fListZDCQxy->Add(fHist_znAx_V0_VxVy[0][j]);
         fListZDCQxy->Add(fHist_znAy_V0_VxVy[0][j]);
      }
    }
  }


  if(sCentEstimator=="V0"){
    fHist_Task_config->Fill(1.5);
  }
  if(sCentEstimator=="TPC"){
    fHist_Task_config->Fill(2.5);
  }
  if(sCentEstimator=="CL1"){
    fHist_Task_config->Fill(3.5);
  }
  if(fRejectPileUp){
    fHist_Task_config->Fill(4.5);
  }
  if(fRejectPileUpTight){
    fHist_Task_config->Fill(5.5);
  }
  if(fAnalysisSet=="FillGainEq"){
    fHist_Task_config->Fill(6.5);
  }
  if(fAnalysisSet=="DoGainEq"){
    fHist_Task_config->Fill(7.5);
  }
  if(fListZDCWgt){
    fHist_Task_config->Fill(8.5);
  }
  if(bFillCosSin){
    fHist_Task_config->Fill(9.5);
  }
  if(!bRunAveragedQn){
    fHist_Task_config->Fill(12.5);
  }



  fHist_CutParameters->SetBinContent(1,VzCut[0]);
  fHist_CutParameters->SetBinContent(2,VzCut[1]);
  fHist_CutParameters->SetBinContent(3,vzBin);
  fHist_CutParameters->SetBinContent(4,VxCut[0]);
  fHist_CutParameters->SetBinContent(5,VxCut[1]);
  fHist_CutParameters->SetBinContent(6,vxBin);
  fHist_CutParameters->SetBinContent(7,VyCut[0]);
  fHist_CutParameters->SetBinContent(8,VyCut[1]);
  fHist_CutParameters->SetBinContent(9,vyBin);



  PostData(1,fListHistos); 

  if(fAnalysisSet=="DoGainEq") {
    PostData(2,fListZDCQxy); 
  }
  else{
    fListDummy1 = new TList();
    fListDummy1->SetOwner(kTRUE);
    fListDummy1->Add(fHist_Event_count);
    PostData(2,fListDummy1); 
  }


  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);


  fcheckOnce = 0;
  fOldRunNum = 0;
  printf("\n\n::UserCreateOutPutObject()\n Runflag= %d, dataset: %s, VxyBin = %d, Analysis= %s\n\n",frunflag,fDataSet.Data(),NbinVt,fAnalysisSet.Data());

}

//________________________________________________________________________
void AliAnalysisTaskZDCGainEq::UserExec(Option_t *)
{

  int stepCount = 0;

  //printf("\n ... ::UserExec() is being called. 1 Step %d...  \n",stepCount);

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
  fEvent           = dynamic_cast<AliFlowEventSimple*>(GetInputData(1));


  if(!aod){
    printf("\n ... ::UserExec = no aod found.....  \n");
    return;
  }

  fHist_Event_count->Fill(stepCount);
  stepCount++;


  AliAODVertex *pVertex = aod->GetPrimaryVertex();
  Double_t Vxyz[3]    =         {0,0,0};
  Vxyz[0]             = pVertex->GetX();
  Vxyz[1]             = pVertex->GetY();
  Vxyz[2]             = pVertex->GetZ();

  //------- Apply Necessary event cuts ---------
  if(Vxyz[2] >= VzCut[1]  || Vxyz[2] <= VzCut[0])  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(Vxyz[0] >= VxCut[1]  || Vxyz[0] <= VxCut[0])  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(Vxyz[1] >= VyCut[1]  || Vxyz[1] <= VyCut[0])  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;


  Double_t EvtCent = fEvent->GetCentrality();


  //---------- get runindex: --------------
  Int_t runNumber = aod->GetRunNumber();
  Int_t runindex = -111;

 for(int i=0;i<frunflag;i++){
   if(runNumber==runNums[i])
    {
      runindex = i;
      break;
    }
  }
 if(runindex<0) {
    printf("\n ... **WARNING** \n::UserExec() runnumber not listed.\n EXIT..\n");
    //exit(1);
  }
 //-----------------------------------------











 //--------- starting pileup rejection work: --------
   Double_t centrV0M=300;
   Double_t centrCL1=300;
   Double_t centrCL0=300;
   Double_t centrTRK=300;

  if(fDataSet=="2010"||fDataSet=="2011"){
    centrV0M = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
    centrCL1 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");
    centrCL0 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL0");
    centrTRK = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("TRK");
  }
  else{
    fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
     if(!fMultSelection) {
       printf("\n **WARNING** ::UserExec() AliMultSelection object not found. Step# %d\n",stepCount);
       exit(1);
     }
    centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
    centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
    centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");
  }// 2015

  Bool_t BisPileup=kFALSE;

  if(fRejectPileUp && InputEvent()) {
    //if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
    if(fDataSet!="2015") {
          if(plpMV(aod)) {
            fPileUpCount->Fill(0.5);
            BisPileup=kTRUE;
          }
          Int_t isPileup = aod->IsPileupFromSPD(3);
          if(isPileup != 0) {
            fPileUpCount->Fill(1.5);
            //BisPileup=kTRUE;       //  ----- Rihan: pileup from SPD not used for 2010
          }
          if(((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
            fPileUpCount->Fill(2.5);
            BisPileup=kTRUE;
          }
          if(aod->IsIncompleteDAQ())  {
            fPileUpCount->Fill(3.5);
            BisPileup=kTRUE;
          }

    //check vertex consistency
     const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
     const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();

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
        }

        else {
          // pileup from AliMultSelection for 2015
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

      // pile-up a la Dobrin for LHC15o
          if(plpMV(aod)) {
            fPileUpCount->Fill(0.5);
            BisPileup=kTRUE;
          }

          Int_t isPileup = aod->IsPileupFromSPD(3);
          if(isPileup != 0) {
            fPileUpCount->Fill(1.5);
            BisPileup=kTRUE;          
          }

          if(((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
            fPileUpCount->Fill(2.5);
            BisPileup=kTRUE;
          }

          if(aod->IsIncompleteDAQ())  {
            fPileUpCount->Fill(3.5);
            BisPileup=kTRUE;
          }

          if(fabs(centrV0M-centrCL1)>7.5)  {
            fPileUpCount->Fill(4.5);
            BisPileup=kTRUE;
          }

     // check vertex consistency
     const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
     const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();

     if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
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

       // cuts on tracks
          const Int_t nTracks = aod->GetNumberOfTracks();
          Int_t multEsd = ((AliAODHeader*)aod->GetHeader())->GetNumberOfESDTracks();

          Int_t multTrk = 0;
          Int_t multTrkBefC = 0;
          Int_t multTrkTOFBefC = 0;
          Int_t multTPC = 0;

    for(Int_t it = 0; it < nTracks; it++) {
       AliAODTrack* aodTrk = (AliAODTrack*)aod->GetTrack(it);
        if(!aodTrk) {
          delete aodTrk;
          continue;
         }
//      if(aodTrk->TestFilterBit(32)){
//         multTrkBefC++;
//       if(TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
//         multTrkTOFBefC++;
//         if((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2) && (aodTrk->Pt() < 20.))
//            multTrk++;
//       }
         if(aodTrk->TestFilterBit(128))
              multTPC++;
        } // end of for (Int_t it = 0; it < nTracks; it++)

      Double_t multTPCn = multTPC;
      Double_t multEsdn = multEsd;
      Double_t multESDTPCDif = multEsdn - multTPCn*3.38;

     if(multESDTPCDif > (fRejectPileUpTight?700.:15000.)) {
        fPileUpCount->Fill(7.5);
        BisPileup=kTRUE;
      }

     if(fRejectPileUpTight) {
       if(BisPileup==kFALSE) {
              if(!fMultSelection->GetThisEventIsNotPileup()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
              if(BisPileup) fPileUpCount->Fill(8.5);
            }
          }
	}
      }

 //----------- pile up rejection done ---------







 //read Histograms for channel Wgt:
  if(runNumber!=fOldRunNum){  //if runNumber changed, re-read new list.
    fcheckOnce = 0;
  }

  if(!fcheckOnce && fAnalysisSet=="DoGainEq") {
    fHist_ChanWgt_ZDCC = (TH1F *) fListZDCWgt->FindObject(Form("fHist1F_ZDCC_ChannelWgt_Run%d",runNumber));
    fHist_ChanWgt_ZDCA = (TH1F *) fListZDCWgt->FindObject(Form("fHist1F_ZDCA_ChannelWgt_Run%d",runNumber));

    fcheckOnce++;
    fOldRunNum = runNumber;
  }


  Double_t ChanWgtZDCC[4] = {1.,1.,1.,1.};
  Double_t ChanWgtZDCA[4] = {1.,1.,1.,1.};

  Int_t iCentBin = abs(EvtCent) + 1;
  Int_t iWgtBin = -1;

  if(fAnalysisSet=="DoGainEq") {
    if(fHist_ChanWgt_ZDCC && fHist_ChanWgt_ZDCA){
      for(int ich=1; ich<=4;  ich++){
        iWgtBin = 4*(iCentBin-1) + ich;
        ChanWgtZDCC[ich-1] = fHist_ChanWgt_ZDCC->GetBinContent(iWgtBin);
        ChanWgtZDCA[ich-1] = fHist_ChanWgt_ZDCA->GetBinContent(iWgtBin);
      }
    }
    else{
        printf("\n\n **WARNING**\n ZDC Channel Weights not found. Using weights = 1.0 \n\n");
        //exit(1);
    }
  }



  fHist_Cent_woZDCcut->Fill(EvtCent);


  //----------- Read ZDC information ----------
  AliAODZDC *aodZDC  = aod->GetZDCData();
  Float_t                            fZDCGainAlpha = 0.500;  // fZDCGainAlpha : Jacopo uses 0.35 ??
  Float_t energyZNC  = (Float_t)  (aodZDC->GetZNCEnergy());
  Float_t energyZPC  = (Float_t)  (aodZDC->GetZPCEnergy());
  Float_t energyZNA  = (Float_t)  (aodZDC->GetZNAEnergy());
  Float_t energyZPA  = (Float_t)  (aodZDC->GetZPAEnergy());
  Float_t energyZEM1 = (Float_t) (aodZDC->GetZEM1Energy());
  Float_t energyZEM2 = (Float_t) (aodZDC->GetZEM2Energy());

  const Double_t * towZNC   = aodZDC->GetZNCTowerEnergy();
  const Double_t * towZPC   = aodZDC->GetZPCTowerEnergy();
  const Double_t * towZNA   = aodZDC->GetZNATowerEnergy();
  const Double_t * towZPA   = aodZDC->GetZPATowerEnergy();

  const Double_t * towZNClg = aodZDC->GetZNCTowerEnergyLR(); // Low gain something, should not be used.
  const Double_t * towZNAlg = aodZDC->GetZNATowerEnergyLR();

  Double_t towZPClg[5] = {0.,};
  Double_t towZPAlg[5] = {0.,};

  for(Int_t it=0; it<5; it++) {
    towZPClg[it] = 8*towZPC[it];
    towZPAlg[it] = 8*towZNA[it];
  }

  //sanity: remove if any of ZDC_C_A has negetive Energy:
  if(towZNC[1]<0 || towZNC[2]<0 || towZNC[3]<0 || towZNC[4]<0)  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(towZNA[1]<0 || towZNA[2]<0 || towZNA[3]<0 || towZNA[4]<0)  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;


//********** Get centroid from ZDCs **************

  Double_t xyZNC[2]={999.,999.};
  Double_t xyZNA[2]={999.,999.};

  Float_t zncEnergy=0., znaEnergy=0.;

  for(Int_t i=0; i<5; i++){
    zncEnergy += towZNC[i];
    znaEnergy += towZNA[i];
  }



  

/*-----------------------Not used---------------------------------
  Int_t CenBin = GetCenBin(centrperc);
  Double_t zvtxpos[3]={0.,0.,0.};
  fFlowEvent->GetVertexPosition(zvtxpos);
  Int_t RunNum=fFlowEvent->GetRun();
  if(fTowerEqList) {
   if(RunNum!=fCachedRunNum) {
      for(Int_t i=0; i<8; i++) {
      fTowerGainEq[i] = (TH1D*)(fTowerEqList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fhnTowerGainEqFactor[%d][%d]",RunNum,i)));
     }
   }
 }
  Bool_t fUseMCCen = kFALSE; //rihan:hardcoded
  if (fUseMCCen) {
   if(aod->GetRunNumber() < 209122)
    aodZDC->GetZNCentroidInPbPb(1380., xyZNC, xyZNA);
    else
    aodZDC->GetZNCentroidInPbPb(2510., xyZNC, xyZNA);
  }
  else {
  //set tower gain equalization, if available
  if(fTowerEqList) {
   for(Int_t i=0; i<8; i++)
   {
    if(fTowerGainEq[i])
    AvTowerGain[i] = fTowerGainEq[i]->GetBinContent(fTowerGainEq[i]->FindBin(centrperc));
   }
 }//---------------------------------------------------------------  */



  Double_t  towCalibZNC[5] = {0,}; 
  Double_t  towCalibZNA[5] = {0,};

// towZNC[] is constant; so Need to make a copy
  for(int i=0;i<5;i++){
    towCalibZNC[i] = towZNC[i];
    towCalibZNA[i] = towZNA[i];
  }

  // Now calibrate the energy of channel 1-4:
  for(int i=0;i<4;i++){
    towCalibZNC[i+1] = ChanWgtZDCC[i]*towCalibZNC[i+1];
    if(ChanWgtZDCA[i] < 4.){
      towCalibZNA[i+1] = ChanWgtZDCA[i]*towCalibZNA[i+1];
    }
  }
 //manually get Energy in ZDC-A channel 2:
  if(ChanWgtZDCA[1] >= 4.0){
    towCalibZNA[2] = towCalibZNA[0] - towCalibZNA[1] - towCalibZNA[3] - towCalibZNA[4];
  }



  Double_t AvTowerGain[8] = {1., 1., 1., 1., 1., 1., 1., 1.};
  
  const Float_t x[4] = {-1.75,  1.75,-1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};

  Float_t numXZNC=0., numYZNC=0., denZNC=0., wZNC;
  Float_t numXZNA=0., numYZNA=0., denZNA=0., wZNA;

  for(Int_t i=0; i<4; i++) {
    if(towCalibZNC[i+1]>0.) {
       wZNC = TMath::Power(towCalibZNC[i+1], fZDCGainAlpha)*AvTowerGain[i];
       numXZNC += x[i]*wZNC;
       numYZNC += y[i]*wZNC;
       denZNC  += wZNC;
    }

    if(towCalibZNA[i+1]>0.) {
       wZNA = TMath::Power(towCalibZNA[i+1], fZDCGainAlpha)*AvTowerGain[i+4];
       numXZNA += x[i]*wZNA;
       numYZNA += y[i]*wZNA;
       denZNA  += wZNA;
    }
  }

  if(denZNC!=0) {
    xyZNC[0] = numXZNC/denZNC;
    xyZNC[1] = numYZNC/denZNC;
  }
  else{
     xyZNC[0]  = 999.;
     xyZNC[1]  = 999.;
     zncEnergy =   0.;
  }
  if(denZNA!=0) {
     xyZNA[0] = numXZNA/denZNA;
     xyZNA[1] = numYZNA/denZNA;
  }
  else{
     xyZNA[0]  = 999.;
     xyZNA[1]  = 999.;
     znaEnergy =   0.;
  }


  
  xyZNA[0] = -1.*xyZNA[0]; //----- Important: zdcA_X = -zdcA_X ---------



  if(sqrt(xyZNC[0]*xyZNC[0] + xyZNC[1]*xyZNC[1])>1.5)  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(sqrt(xyZNA[0]*xyZNA[0] + xyZNA[1]*xyZNA[1])>1.5)  return;

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  fHist_Event_count->Fill(stepCount);
  stepCount++;



  Int_t indexVx = fHist_Vx_ArrayFinder->FindBin(Vxyz[0]);
  Int_t indexVy = fHist_Vy_ArrayFinder->FindBin(Vxyz[1]);
  Int_t indexVz = fHist_Vz_ArrayFinder->FindBin(Vxyz[2]);
  Double_t tVertexBin1 = 0;

  if(fAnalysisSet=="DoGainEq"){

    tVertexBin1 = (Double_t) (indexVy-1)*vxBin + (Double_t)indexVx - 0.5 ; 

    if(bFillCosSin) {
       double psi1C = TMath::ATan2(xyZNC[1],xyZNC[0]);
       if(psi1C<0){
         psi1C += 2.*TMath::Pi();
       }
       fHist_Psi1_ZDCC_wGainCorr->Fill(psi1C);

       double psi1A = TMath::ATan2(xyZNA[1],xyZNA[0]);
       if(psi1A<0){
         psi1A += 2.*TMath::Pi();
       }
       fHist_Psi1_ZDCA_wGainCorr->Fill(psi1A);

       if(!bRunAveragedQn){
         fHist_znCx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Cos(psi1C)); 
         fHist_znCy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Sin(psi1C));
         fHist_znAx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Cos(psi1A));
         fHist_znAy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Sin(psi1A));
       }
       else{
         fHist_znCx_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Cos(psi1C)); 
         fHist_znCy_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Sin(psi1C));
         fHist_znAx_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Cos(psi1A));
         fHist_znAy_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,TMath::Sin(psi1A));
       }
    }
    else{

       if(!bRunAveragedQn){
         fHist_znCx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[0]); 
         fHist_znCy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[1]);
         fHist_znAx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[0]);
         fHist_znAy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[1]);
       }
       else{
         fHist_znCx_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[0]); 
         fHist_znCy_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[1]);
         fHist_znAx_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[0]);
         fHist_znAy_V0_VxVy[0][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[1]);
       }
    }
  }


  Int_t nRefMult = fEvent->GetReferenceMultiplicity();

  if(bFillZDCQAon){

    Double_t  FillVsWith[5]  = {EvtCent,static_cast<Double_t>(nRefMult), Vxyz[0], Vxyz[1], Vxyz[2]};
    Double_t  FillValueQx[4] = {xyZNA[0],xyZNC[0],xyZNA[1],xyZNC[1]};
    Double_t  FillValueXX[4] = {xyZNA[0]*xyZNC[0],xyZNA[1]*xyZNC[1],xyZNC[0]*xyZNA[1],xyZNC[1]*xyZNA[0]}; //XaXc,YaYc,XcYa,YcXa

    //fill the uncorrected Qx,XX.  /* Corrected to be added later */
    for(int i=0;i<4;i++){
     for(int j=0;j<5;j++){
        fHist_Qx_vs_Obs_woCorr[i][j]->Fill(FillVsWith[j],FillValueQx[i]);
        fHist_XX_vs_Obs_woCorr[i][j]->Fill(FillVsWith[j],FillValueXX[i]);
     }
    }
  }

  fHist_Cent_wiZDCcut->Fill(EvtCent);


  Int_t    iTracks = fEvent->NumberOfTracks();


  //if(fievent%20==0){
  //std::cout<<fievent<<" run = "<<runNumber<<" cTPC= "<<EvtCent<<"\tRefMult = "<<nRefMult
  //         <<"\t towZNA[2] = "<<towZNA[2]<<"\t towCalibZNA[2] = "<<towCalibZNA[2]<<std::endl;    } 
  

  PostData(1,fListHistos);

  if(fAnalysisSet=="DoGainEq"){
    PostData(2,fListZDCQxy); 
  }
  else{
    PostData(2,fListDummy1);  //default if no 'fAnalysisSet' found.
  }


  fievent++;

  //printf("\n ... ::UserExec() is being called. Step last %d... Event %d  \n",stepCount,fievent);
  
}// UserExec ends




void AliAnalysisTaskZDCGainEq::Terminate(Option_t *)
{
      AliDebug(2,"\n ... AliAnalysisTaskZDCGainEq::Terminate() is being called ...  \n");
}


double AliAnalysisTaskZDCGainEq::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
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

 Bool_t AliAnalysisTaskZDCGainEq::plpMV(const AliAODEvent* aod)
 {  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if(!(nPlp=aod->GetNumberOfPileupVerticesTracks()))
  return kFALSE;

  vtPrm = aod->GetPrimaryVertex();
  if(vtPrm == aod->GetPrimaryVertexSPD())
  return kTRUE;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for(int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)aod->GetPileupVertexTracks(ipl);
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




