/*
  Unit test for some functions classes  used in the $ALICE_ROOT/TPC/Base directory:
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER\
   -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW -I$ALICE_ROOT/PWG1 -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/Base -I$ALICE_ROOT/TPC/Calib");
   
  .L $ALICE_ROOT/TPC/Base/test/UnitTest.C+ 
  UnitTestAliTPCCalPadTree();
  TestCorrection_AliTPCExBTwistAddCorrectionCompact();
  TestCorrection_AliTPCFCVoltError3DAddCorrectionCompact();
  TestCorrection_AliTPCRocVoltError3DAddCorrectionCompact();

*/

#include "TF1.h"
#include "TMath.h"
#include "TLinearFitter.h"
#include "TFile.h"
#include "AliSysInfo.h"
#include "TTree.h"
#include "AliLog.h"
#include "THn.h"
#include "TRandom.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalibViewer.h"
#include "AliTPCcalibDButil.h"
#include "AliTPCCorrection.h"
#include "AliTPCComposedCorrection.h"
#include "AliTPCExBTwist.h"
#include "AliTPCFCVoltError3D.h"
#include "AliTPCROCVoltError3D.h"
//#include "AliTPCBoundaryVoltError.h"

//
// PARAMETERS to set from outside:
//
TString baseDir="/hera/alice/wiechula/calib/guiTrees";  // TO  FIX specification of inout data
//
//


void  UnitTestAliTPCCalPadTree(){
  //
  //  Make a UnitTest of the AliTPCCalPad 
  //   a.) TTree functionaility
  //   b.) MedianFilterFunctionality
  //   c.) LTMFilterFunctionality
  // 
  TObjArray *fArray = new TObjArray(100);
  TTree * treePad=AliTPCcalibDButil::ConnectGainTrees(baseDir);
  for (Int_t i=0; i<5; i+=2){
    AliTPCCalPad * padLx = AliTPCCalPad::MakePadFromTree(treePad,"lx.fElements","Lx",kTRUE);
    AliTPCCalPad * padLy = AliTPCCalPad::MakePadFromTree(treePad,"ly.fElements","Ly",kTRUE);
    AliTPCCalPad * padLLx = AliTPCCalPad::MakePadFromTree(treePad,"lx.fElements","LLx",kTRUE);
    AliTPCCalPad * padLLy = AliTPCCalPad::MakePadFromTree(treePad,"ly.fElements","LLy",kTRUE);
    AliTPCCalPad * padMax = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MaxCharge.fElements","QMax",kTRUE);
    AliTPCCalPad * padMean = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MeanCharge.fElements","QTot",kTRUE);
    AliTPCCalPad * padMaxL = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MaxCharge.fElements","QMax",kTRUE);
    AliTPCCalPad * padMeanL = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MeanCharge.fElements","QTot",kTRUE);
    if (i>0) {
      padLx->MedianFilter(i,2*i);
      padLy->MedianFilter(i,2*i);
      padLLx->LTMFilter(i,2*i,1.00, 0);
      padLLy->LTMFilter(i,2*i,1.00, 0);
      padMax->MedianFilter(i,2*i);
      padMean->MedianFilter(i,2*i);
      padMaxL->LTMFilter(i,2*i,0.8,0);
      padMeanL->LTMFilter(i,2*i,0.8,0);
    }
    padLx->SetName(TString::Format("Lx%d",i).Data());
    padLy->SetName(TString::Format("Ly%d",i).Data());
    padLLx->SetName(TString::Format("LLx%d",i).Data());
    padLLy->SetName(TString::Format("LLy%d",i).Data());
    padMax->SetName(TString::Format("QMax%d",i).Data());
    padMean->SetName(TString::Format("QTot%d",i).Data());
    padMaxL->SetName(TString::Format("QMaxL%d",i).Data());
    padMeanL->SetName(TString::Format("QTotL%d",i).Data());
    fArray->AddLast(padLx);
    fArray->AddLast(padLy);
    fArray->AddLast(padLLx);
    fArray->AddLast(padLLy);
    fArray->AddLast(padMax);
    fArray->AddLast(padMean);
    fArray->AddLast(padMaxL);
    fArray->AddLast(padMeanL);
  }
  AliTPCCalibViewer::MakeTree("QAtest.root", fArray,0);
  //
  // 2.) Check invariants
  //
  TFile*fout= TFile::Open("QAtest.root");
  TTree * tree  = (TTree*)fout->Get("calPads");
  Int_t isOutM0 = tree->Draw("(Ly2.fElements-Ly0.fElements)>>his0(100,-10,10)","abs((Ly2.fElements-Ly0.fElements))>2","goff");
  Int_t isOutM1=tree->Draw("(Lx2.fElements-Lx0.fElements)/0.75>>his1(100,-10,10)","abs((Lx2.fElements-Lx0.fElements))>0","goff");
  printf("IsOut=%d\t%d\n",isOutM0,isOutM1);
  if ((isOutM0+isOutM1)==0) ::Info("UnitTestAliTPCCalPadTree","MedianTest OK");
  if (isOutM0||isOutM1) ::Fatal("UnitTestAliTPCCalPadTree","MedianTest FAILED");
  //
  Int_t isOutL0 = tree->Draw("(LLy2.fElements-Ly0.fElements)>>his0(100,-10,10)","abs((LLy2.fElements-LLy0.fElements))>0","goff");
  Int_t isOutL1=tree->Draw("(LLx2.fElements-Lx0.fElements)/0.75>>his1(100,-10,10)","abs((LLx2.fElements-LLx0.fElements))>0","goff");
  printf("IsOut=%d\t%d\n",isOutL0,isOutL1);
  if ((isOutL0+isOutL1)==0) ::Info("UnitTestAliTPCCalPadTree","LTMTest OK");
  if (isOutL0||isOutL1) ::Fatal("UnitTestAliTPCCalPadTree","LTMTest FAILED");
}


Bool_t  TestCorrection_AliTPCExBTwistAddCorrectionCompact(){
  //
  // 
  // 1.) Test ExB twist AddCorrectionCompact
  //
  Bool_t isOK[10]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  AliTPCComposedCorrection *compCorrTwist = new AliTPCComposedCorrection();
  AliTPCExBTwist  *twistX    = new  AliTPCExBTwist;
  AliTPCExBTwist  *twistY    = new  AliTPCExBTwist;
  twistX->SetXTwist(0.001);  // 1 mrad twist in x
  twistY->SetYTwist(0.001);  // 1 mrad twist in x
  isOK[0]&=compCorrTwist->AddCorrectionCompact(twistX,1);
  isOK[0]&=compCorrTwist->AddCorrectionCompact(twistY,1);
  isOK[0]&=compCorrTwist->AddCorrectionCompact(twistY,-1);
  isOK[0]&=compCorrTwist->AddCorrectionCompact(twistX,-1);
  isOK[1]=compCorrTwist->GetCorrections()->GetEntries()==1;
  AliTPCExBTwist  *twistRes=0;
  if (isOK[1]==kFALSE){
    isOK[2]=kFALSE;
    isOK[3]=kFALSE;
    isOK[4]=kFALSE;
  }else{
    twistRes=  dynamic_cast<AliTPCExBTwist *>(compCorrTwist->GetSubCorrection(0));
    if (twistRes==NULL){
      isOK[2]=kFALSE;
      isOK[3]=kFALSE;
      isOK[4]=kFALSE;
    }else{
      isOK[3] &= (twistRes->GetXTwist()==0);
      isOK[4] &= (twistRes->GetYTwist()==0);
    }
  }
  Bool_t res=kTRUE;
  for (Int_t i=0; i<5; i++) res&=isOK[i];
  {
    if (isOK[0]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCExBTwist -ADD FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCExBTwist -ADD OK");
    }
    if (isOK[1]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCExBTwist - wrong entries  FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCExBTwist - entries  OK");
    }
    if (isOK[2]==kFALSE || isOK[3]==kFALSE ||isOK[4]==kFALSE ){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCExBTwist - inconsitent entries  FAILED");    
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCExBTwist - consistent entries  OK");    
    }
  }    
  return res;
} 


Bool_t  TestCorrection_AliTPCFCVoltError3DAddCorrectionCompact(){
  //
  // TestCorrection_AliTPCFCVoltError3DAddCorrectionCompact
  //
  const Float_t kEpsilon=0.000001;
  Bool_t isOK[10]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  AliTPCComposedCorrection *compCorrComp = new AliTPCComposedCorrection();
  AliTPCFCVoltError3D  *corr0    = new  AliTPCFCVoltError3D;
  AliTPCFCVoltError3D  *corr1    = new  AliTPCFCVoltError3D;
  for (Int_t isec=0; isec<36; isec++){
    corr0->SetRodVoltShiftA(isec,TMath::Cos(TMath::Pi()*isec/36),kFALSE);
    corr0->SetRodVoltShiftC(isec,TMath::Cos(TMath::Pi()*isec/36),kFALSE);
    corr1->SetRodVoltShiftA(isec,TMath::Sin(TMath::Pi()*isec/36),kFALSE);
    corr1->SetRodVoltShiftC(isec,TMath::Sin(TMath::Pi()*isec/36),kFALSE);
    corr1->SetCopperRodShiftA(isec,TMath::Sin(TMath::Pi()*isec/36),kFALSE);
    corr1->SetCopperRodShiftC(isec,TMath::Sin(TMath::Pi()*isec/36),kFALSE);
  }
  //
  isOK[0]&=compCorrComp->AddCorrectionCompact(corr0,1);
  isOK[0]&=compCorrComp->AddCorrectionCompact(corr1,1);
  isOK[0]&=compCorrComp->AddCorrectionCompact(corr1,-1);
  isOK[0]&=compCorrComp->AddCorrectionCompact(corr0,-1);
  isOK[1]=compCorrComp->GetCorrections()->GetEntries()==1;
  AliTPCFCVoltError3D  *corrRes=0;
  if (isOK[1]==kFALSE){
    isOK[2]=kFALSE;
    isOK[3]=kFALSE;
    isOK[4]=kFALSE;
  }else{
    corrRes=  dynamic_cast<AliTPCFCVoltError3D *>(compCorrComp->GetSubCorrection(0));
    if (corrRes==NULL){
      isOK[2]=kFALSE;
      isOK[3]=kFALSE;
      isOK[4]=kFALSE;
    }else{
      for (Int_t isec=0; isec<36; isec++){
	isOK[3] &=( TMath::Abs(corrRes->GetRodVoltShiftA(isec))<kEpsilon);
	isOK[4] &=( TMath::Abs(corrRes->GetRodVoltShiftC(isec))<kEpsilon);
	isOK[5] &=( TMath::Abs(corrRes->GetCopperRodShiftA(isec))<kEpsilon);
	isOK[6] &=( TMath::Abs(corrRes->GetCopperRodShiftC(isec))<kEpsilon);
      }
    }
  }
  Bool_t res=kTRUE;
  for (Int_t i=0; i<5; i++) res&=isOK[i];
  {
    if (isOK[0]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCFCVoltError3D -ADD FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCFCVoltError3D -ADD OK");
    }
    if (isOK[1]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCFCVoltError3D - wrong entries  FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCFCVoltError3D - entries  OK");
    }
    if (isOK[2]==kFALSE || isOK[3]==kFALSE ||isOK[4]==kFALSE ){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCFCVoltError3D - inconsitent entries  FAILED");    
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCFCVoltError3D - consistent entries  OK");    
    }
  }    
  return res;
}



Bool_t  TestCorrection_AliTPCRocVoltError3DAddCorrectionCompact(){
  //
  // AliTPCRocVoltError3DAddCorrectionCompact
  //
  const Float_t kEpsilon=0.00000001;
  Bool_t isOK[10]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  AliTPCComposedCorrection *compCorrROCVoltError3D = new AliTPCComposedCorrection();
  AliTPCROCVoltError3D  *corr0    = new  AliTPCROCVoltError3D;
  AliTPCROCVoltError3D  *corr1    = new  AliTPCROCVoltError3D;
  TMatrixD matrixDz(72,3);
  for (Int_t isec=0; isec<72; isec++){
    matrixDz(isec,0)=gRandom->Rndm()*0.1;
    matrixDz(isec,1)=gRandom->Rndm()*0.001;
    matrixDz(isec,2)=gRandom->Rndm()*0.001;
  }
  corr0->SetROCData(&matrixDz);
  matrixDz*=0.5;
  corr1->SetROCData(&matrixDz);
  //
  isOK[0]&=compCorrROCVoltError3D->AddCorrectionCompact(corr0,1);
  isOK[0]&=compCorrROCVoltError3D->AddCorrectionCompact(corr1,1);
  isOK[0]&=compCorrROCVoltError3D->AddCorrectionCompact(corr1,-1);
  isOK[0]&=compCorrROCVoltError3D->AddCorrectionCompact(corr0,-1);
  isOK[1]=compCorrROCVoltError3D->GetCorrections()->GetEntries()==1;
  AliTPCROCVoltError3D  *corrRes=0;
  if (isOK[1]==kFALSE){
    isOK[2]=kFALSE;
    isOK[3]=kFALSE;
    isOK[4]=kFALSE;
  }else{
    corrRes=  dynamic_cast<AliTPCROCVoltError3D *>(compCorrROCVoltError3D->GetSubCorrection(0));
    if (corrRes==NULL){
      isOK[2]=kFALSE;
      isOK[3]=kFALSE;
      isOK[4]=kFALSE;
    }else{
      isOK[3]=TMath::Abs(corrRes->GetMatrix()->Sum())<kEpsilon;
    }
  }
  Bool_t res=kTRUE;
  for (Int_t i=0; i<5; i++) res&=isOK[i];
  {
    if (isOK[0]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCROCVoltError3D -ADD FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCROCVoltError3D -ADD OK");
    }
    if (isOK[1]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCROCVoltError3D - wrong entries  FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCROCVoltError3D - entries  OK");
    }
    if (isOK[2]==kFALSE || isOK[3]==kFALSE ||isOK[4]==kFALSE ){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCROCVoltError3D - inconsitent entries  FAILED");    
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCROCVoltError3D - consistent entries  OK");    
    }
  }    
  return res;
} 
