/// \file UnitTest.C
///
/// Unit test for some functions classes  used in the $ALICE_ROOT/TPC/Base directory:
///
/// ~~~{.cpp}
/// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/install/include -I$ALICE_ROOT/STEER   -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW  -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/TPCbase -I$ALICE_ROOT/TPCcalib");
///
/// .L $ALICE_ROOT/TPC/Base/test/UnitTest.C+
/// UnitTestAliTPCCalPadTree();
/// TestCorrection_AliTPCCorrection_AddCorrectionCompact();
/// ~~~

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
#include "AliTPCBoundaryVoltError.h"
#include "AliTPCCalibGlobalMisalignment.h"
#include "AliCDBEntry.h"
#include "TStopwatch.h"
#include "TGeoMatrix.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
// PARAMETERS to set from outside:
TString baseDir="/hera/alice/wiechula/calib/guiTrees";  // TO  FIX specification of inout data
//
//

Bool_t  TestCorrection_AliTPCCorrection_AddCorrectionCompact();
Bool_t  TestCorrection_AliTPCExBTwistAddCorrectionCompact();
Bool_t  TestCorrection_AliTPCFCVoltError3DAddCorrectionCompact();
Bool_t  TestCorrection_AliTPCRocVoltError3DAddCorrectionCompact();
Bool_t  TestCorrection_AliTPCBoundaryVoltErrorAddCorrectionCompact();
Bool_t  TestCorrection_AliTPCCalibGlobalMisalignmentAddCorrectionCompact();
Bool_t  TestCorrection_AliTPCComposedCorrectionAddCorrectionCompact();
Bool_t TestCorrection_AliTPCComposedCorrectionAddCorrectionCompact_TPCCalibCorrection(Bool_t fast=kFALSE);

void  UnitTestAliTPCCalPadTree(){
  ///  Make a UnitTest of the AliTPCCalPad
  ///   a.) TTree functionaility
  ///   b.) MedianFilterFunctionality
  ///   c.) LTMFilterFunctionality

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
  /// 1.) Test ExB twist AddCorrectionCompact

  Bool_t isOK[10]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  AliTPCComposedCorrection *compCorrTwist = new AliTPCComposedCorrection();
  AliTPCExBTwist  *twistX    = new  AliTPCExBTwist;
  AliTPCExBTwist  *twistY    = new  AliTPCExBTwist;
  twistX->SetXTwist(0.001);  // 1 mrad twist in x
  twistY->SetYTwist(0.001);  // 1 mrad twist in x
  isOK[0]&=compCorrTwist->AddCorrectionCompact(twistX,0.5);
  isOK[0]&=compCorrTwist->AddCorrectionCompact(twistY,0.5);
  isOK[0]&=compCorrTwist->AddCorrectionCompact(twistX,0.5);
  isOK[0]&=compCorrTwist->AddCorrectionCompact(twistY,0.5);
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
  /// TestCorrection_AliTPCFCVoltError3DAddCorrectionCompact

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
  isOK[0]&=compCorrComp->AddCorrectionCompact(corr0,0.5);
  isOK[0]&=compCorrComp->AddCorrectionCompact(corr1,0.5);
  isOK[0]&=compCorrComp->AddCorrectionCompact(corr0,0.5);
  isOK[0]&=compCorrComp->AddCorrectionCompact(corr1,0.5);
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
  /// AliTPCRocVoltError3DAddCorrectionCompact

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
  isOK[0]&=compCorrROCVoltError3D->AddCorrectionCompact(corr0,0.5);
  isOK[0]&=compCorrROCVoltError3D->AddCorrectionCompact(corr1,0.5);
  isOK[0]&=compCorrROCVoltError3D->AddCorrectionCompact(corr0,0.5);
  isOK[0]&=compCorrROCVoltError3D->AddCorrectionCompact(corr1,0.5);
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





Bool_t  TestCorrection_AliTPCBoundaryVoltErrorAddCorrectionCompact(){
  /// AliTPCBoundaryVoltErrorAddCorrectionCompact

  const Float_t kEpsilon=0.00000001;
  Bool_t isOK[10]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  AliTPCComposedCorrection *compCorrBoundaryVoltError = new AliTPCComposedCorrection();
  AliTPCBoundaryVoltError  *corr0    = new  AliTPCBoundaryVoltError;
  AliTPCBoundaryVoltError  *corr1    = new  AliTPCBoundaryVoltError;
  Float_t boundaries[8];
  for (Int_t ibound=0; ibound<8; ibound++){ 
    boundaries[ibound]=gRandom->Rndm()-0.5;
  }
  corr0->SetBoundariesA(boundaries);
  corr1->SetBoundariesA(boundaries);  
  //
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr0,0.5);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr1,0.5);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr0,0.5);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr1,0.5);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr1,-1);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr0,-1);
  isOK[1]=compCorrBoundaryVoltError->GetCorrections()->GetEntries()==1;
  AliTPCBoundaryVoltError  *corrRes=0;
  if (isOK[1]==kFALSE){
    isOK[2]=kFALSE;
    isOK[3]=kFALSE;
    isOK[4]=kFALSE;
  }else{
    corrRes=  dynamic_cast<AliTPCBoundaryVoltError *>(compCorrBoundaryVoltError->GetSubCorrection(0));
    if (corrRes==NULL){
      isOK[2]=kFALSE;
      isOK[3]=kFALSE;
      isOK[4]=kFALSE;
    }else{
      for (Int_t ibound=0; ibound<8; ibound++){ 
	isOK[3]&=TMath::Abs(corrRes->GetBoundariesA(ibound))<kEpsilon;
	isOK[3]&=TMath::Abs(corrRes->GetBoundariesC(ibound))<kEpsilon;
      }
    }
  }
  Bool_t res=kTRUE;
  for (Int_t i=0; i<5; i++) res&=isOK[i];
  {
    if (isOK[0]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCBoundaryVoltError -ADD FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCBoundaryVoltError -ADD OK");
    }
    if (isOK[1]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCBoundaryVoltError - wrong entries  FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCBoundaryVoltError - entries  OK");
    }
    if (isOK[2]==kFALSE || isOK[3]==kFALSE ||isOK[4]==kFALSE ){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCBoundaryVoltError - inconsitent entries  FAILED");    
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCBoundaryVoltError - consistent entries  OK");    
    }
  }    
  return res;
} 





Bool_t  TestCorrection_AliTPCCalibGlobalMisalignmentAddCorrectionCompact(){
  /// AliTPCCalibGlobalMisalignmentAddCorrectionCompact
  /// Invariant used in test is not exact it is only approximate - as matrix multiplication is not comulative
  ///  !!!! BUG FOUND ????
  ///  hmatrix1->GetTranslation()[idelta]=xxx; // does not work as expected Translation is set,  visible in Print but not used  later

  const Float_t kEpsilon=0.0001;
  Bool_t isOK[10]={kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE};
  Double_t delta[3]={0.01,0.02,0.03};
  
  AliTPCComposedCorrection *compCorrBoundaryVoltError = new AliTPCComposedCorrection();
  AliTPCCalibGlobalMisalignment  *corr0    = new  AliTPCCalibGlobalMisalignment;
  AliTPCCalibGlobalMisalignment  *corr1    = new  AliTPCCalibGlobalMisalignment;
  AliTPCCalibGlobalMisalignment  *corr2    = new  AliTPCCalibGlobalMisalignment;
  TObjArray sectorAlign(72);
  
  TGeoHMatrix *hmatrix0 = new TGeoHMatrix;
  TGeoHMatrix *hmatrix1 = new TGeoHMatrix;
  hmatrix0->RotateX(TMath::RadToDeg()*0.0001);
  hmatrix0->RotateY(TMath::RadToDeg()*0.0002);
  hmatrix0->RotateZ(TMath::RadToDeg()*0.0003);
  hmatrix1->SetTranslation(delta);
  for (Int_t isec=0; isec<72; isec++){
    if ((isec%2)==0) sectorAlign.AddAt(hmatrix0,isec);
    if ((isec%2)==1) sectorAlign.AddAt(hmatrix1,isec);
  }
  corr0->SetAlignGlobal(hmatrix0);
  corr1->SetAlignGlobal(hmatrix1);
  corr0->SetAlignGlobalDelta(hmatrix1);
  corr1->SetAlignGlobalDelta(hmatrix0);
  corr2->SetAlignSectors(&sectorAlign);
  //
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr0,0.5);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr1,0.5);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr0,0.5);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr1,0.5);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr0,-1);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr1,-1);
  //
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr2,1);
  isOK[0]&=compCorrBoundaryVoltError->AddCorrectionCompact(corr2,-1);
  //
  //
  isOK[1]=compCorrBoundaryVoltError->GetCorrections()->GetEntries()==1;
  AliTPCCalibGlobalMisalignment  *corrRes=0;
  if (isOK[1]==kFALSE){
    isOK[2]=kFALSE;
    isOK[3]=kFALSE;
    isOK[4]=kFALSE;
  }else{
    corrRes=  dynamic_cast<AliTPCCalibGlobalMisalignment *>(compCorrBoundaryVoltError->GetSubCorrection(0));
    if (corrRes==NULL){
      isOK[2]=kFALSE;
      isOK[3]=kFALSE;
      isOK[4]=kFALSE;
    }else{
      for (Int_t itrans=0; itrans<3; itrans++){
	isOK[2+itrans]&=TMath::Abs(corrRes->GetAlignGlobal()->GetTranslation()[itrans])<kEpsilon;
	for (Int_t isec=0; isec<72; isec++){
	  isOK[2+itrans]&=TMath::Abs(((TGeoHMatrix*)(corrRes->GetAlignSectors()->At(isec)))->GetTranslation()[itrans])<kEpsilon;
	}
      }
      corrRes->GetAlignGlobal()->Print();
      corrRes->GetAlignSectors()->At(0)->Print();
      corrRes->GetAlignSectors()->At(1)->Print();      
    }
  }
  Bool_t res=kTRUE; 
  for (Int_t i=0; i<5; i++) res&=isOK[i];
  {
    if (isOK[0]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCCalibGlobalMisalignment -ADD FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCCalibGlobalMisalignment -ADD OK");
    }
    if (isOK[1]==kFALSE){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCCalibGlobalMisalignment - wrong entries  FAILED");
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCCalibGlobalMisalignment - entries  OK");
    }
    if (isOK[2]==kFALSE || isOK[3]==kFALSE ||isOK[4]==kFALSE ){
      ::Error("TestCorrection_AddCorrectionCompact","AliTPCCalibGlobalMisalignment - inconsitent entries  FAILED");    
    }else{
      ::Info("TestCorrection_AddCorrectionCompact","AliTPCCalibGlobalMisalignment - consistent entries  OK");    
    }
  }    
  return res;
} 

Bool_t  TestCorrection_AliTPCCorrection_AddCorrectionCompact(){
  ///

  TestCorrection_AliTPCExBTwistAddCorrectionCompact();
  TestCorrection_AliTPCFCVoltError3DAddCorrectionCompact();
  TestCorrection_AliTPCRocVoltError3DAddCorrectionCompact();
  TestCorrection_AliTPCBoundaryVoltErrorAddCorrectionCompact();
  TestCorrection_AliTPCCalibGlobalMisalignmentAddCorrectionCompact();
}

Bool_t TestCorrection_AliTPCComposedCorrectionAddCorrectionCompact_TPCCalibCorrection(Bool_t fast){
  /// Test the

  const Int_t npointsTest=10000;
  const Float_t kEpsilon=0.001;  //10 microns
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG));
  //
  // 0.) Read an input OCDB entry 
  //
  TFile * f = TFile::Open("$ALICE_OCDB/alice/data/2010/OCDB/TPC/Calib/Correction/Run0_999999999_v8_s0.root");
  AliCDBEntry * entry=(AliCDBEntry*)f->Get("AliCDBEntry");
  TObjArray * corrArray  = (TObjArray *)entry->GetObject();
  AliTPCComposedCorrection *compInput = (AliTPCComposedCorrection *)corrArray->At(0);  
  AliTPCComposedCorrection *compInputFast = new AliTPCComposedCorrection;
  Int_t ncorrs = compInput->GetCorrections()->GetEntries();
  TObjArray arrayInputFast(ncorrs);
  //
  // 1.) Test each individual correction
  //
  for (Int_t icorr=0; icorr<ncorrs; icorr++){
    TString clName=compInput->GetSubCorrection(icorr)->IsA()->GetName();
    if (fast){ // skip slow correction
      if ( clName.Contains("AliTPCFCVoltError3D"))continue;
      if ( clName.Contains("AliTPCROCVoltError3D"))continue;
    }
    //    if ( clName.Contains("AliTPCExBBShape"))continue;
    AliTPCCorrection *corrInput=compInput->GetSubCorrection(icorr);
    //
    ::Info("TestCorrection_AliTPCComposedCorrectionAddCorrectionCompact_TPCCalibCorrection",TString::Format("%s\t%s",corrInput->IsA()->GetName(),corrInput->GetName()).Data());
    AliTPCComposedCorrection *compTest0= new  AliTPCComposedCorrection;
    AliTPCComposedCorrection *compTest1= new  AliTPCComposedCorrection;
    compTest0->AddCorrectionCompact(corrInput,0.5);
    compTest0->AddCorrectionCompact(corrInput,0.5);
    compTest1->AddCorrectionCompact(corrInput,1);
    compTest1->AddCorrectionCompact(corrInput,-1);
    corrInput->AddVisualCorrection(corrInput,10);
    compTest0->AddVisualCorrection(compTest0,11);
    compTest1->AddVisualCorrection(compTest1,12);
    compTest0->SetOmegaTauT1T2(0.35,1,1);
    compTest1->SetOmegaTauT1T2(0.35,1,1);
    corrInput->SetOmegaTauT1T2(0.35,1,1);
    for (Int_t icoord=0; icoord<3; icoord++){
      TVectorD dvecTest0(npointsTest);
      TVectorD dvecTest1(npointsTest);
      for (Int_t ipoint=0; ipoint<npointsTest; ipoint++){
	Double_t r= 85.+gRandom->Rndm()*150;
	Double_t phi= gRandom->Rndm()*TMath::TwoPi();
	Double_t z=500*(gRandom->Rndm()-0.5);
	dvecTest0[ipoint]=AliTPCCorrection::GetCorrXYZ(r*TMath::Cos(phi),r*TMath::Sin(phi),z, icoord, 11)-AliTPCCorrection::GetCorrXYZ(r*TMath::Cos(phi),r*TMath::Sin(phi),z, icoord, 10);
	dvecTest1[ipoint]=AliTPCCorrection::GetCorrXYZ(r*TMath::Cos(phi),r*TMath::Sin(phi),z, icoord, 12);
      }
      Double_t mean0 = TMath::Mean(npointsTest, dvecTest0.GetMatrixArray());
      Double_t rms0  = TMath::RMS(npointsTest, 	dvecTest0.GetMatrixArray());
      Double_t mean1 = TMath::Mean(npointsTest, dvecTest1.GetMatrixArray());
      Double_t rms1  = TMath::RMS(npointsTest, 	dvecTest1.GetMatrixArray());
      if (TMath::Abs(rms0)>kEpsilon){
	::Error("TestCorrection_AliTPCComposedCorrectionAddCorrectionCompact_TPCCalibCorrection",TString::Format("Test0:\t%s\t%3.5f\t%3.5f FAILED",clName.Data(),mean0,rms0).Data());
      }else{
	::Info("TestCorrection_AliTPCComposedCorrectionAddCorrectionCompact_TPCCalibCorrection",TString::Format("Test0:\t%s\t%3.5f\t%3.5f OK",clName.Data(),mean0,rms0).Data());
      }
      if (TMath::Abs(rms1)>kEpsilon){
	::Error("TestCorrection_AliTPCComposedCorrectionAddCorrectionCompact_TPCCalibCorrection",TString::Format("Test1:\t%s\t%3.5f\t%3.5f FAILED",clName.Data(),mean1,rms1).Data());
      }else{
	::Info("TestCorrection_AliTPCComposedCorrectionAddCorrectionCompact_TPCCalibCorrection",TString::Format("Test1:\t%s\t%3.5f\t%3.5f OK",clName.Data(),mean1,rms1).Data());
      }      
    }
  }
}


Bool_t TestCorrection_AliTPCComposedCorrectionAddCorrectionCompact(){
  /// Tests of AliTPCComposedCorrection
  ///  1.) Make linear combination  correction example using weights.
  ///      Test correction  checking invariant inverse x orig  (there are simpler way to do inversion using AliTPCInverseCorrection)
  ///
  ///  2.) Make compact for of the Composed correction. Test correction  checking invariant inverse x orig

  const Int_t npointsTest=10000;
  const Float_t kEpsilon=0.0001;  // using Floating point precission
  //
  // 0.) Read an input OCDB entry 
  //
  TFile * f = TFile::Open("$ALICE_OCDB/alice/data/2010/OCDB/TPC/Calib/Correction/Run0_999999999_v8_s0.root");
  AliCDBEntry * entry=(AliCDBEntry*)f->Get("AliCDBEntry");
  TObjArray * corrArray  = (TObjArray *)entry->GetObject();
  AliTPCComposedCorrection *compInput = (AliTPCComposedCorrection *)corrArray->At(0);  
  AliTPCComposedCorrection *compInputFast = new AliTPCComposedCorrection;
  Int_t ncorrs = compInput->GetCorrections()->GetEntries();
  TObjArray arrayInputFast(ncorrs);
  for (Int_t icorr=0; icorr<ncorrs; icorr++){
    TString clName=compInput->GetSubCorrection(icorr)->IsA()->GetName();
    if ( clName.Contains("AliTPCFCVoltError3D"))continue;
    if ( clName.Contains("AliTPCROCVoltError3D"))continue;
    if ( clName.Contains("AliTPCExBBShape"))continue;
    arrayInputFast.AddLast(compInput->GetSubCorrection(icorr));
  }
  compInputFast->SetCorrections(&arrayInputFast);
  


  //
  // 1.) Make linear combination  correction example using weights.
  //     Test correction  checking invariant inverse x orig  (there are simpler way to do inversion using AliTPCInverseCorrection)
  AliTPCComposedCorrection *compInverse = new AliTPCComposedCorrection();
  Bool_t isOK1[10]={kTRUE};
  TObjArray * collection=   dynamic_cast<TObjArray*>(compInput->GetCorrections());
  Int_t entries = collection->GetEntries();
  TVectorD weights(entries+1);
  for (Int_t i=0; i<entries+1; i++) weights[i]=-1.0;
  weights[0]=1.;
  TObjArray * arrayInvariant = new TObjArray(entries+1);
  arrayInvariant->AddLast(compInput);
  for (Int_t i=0; i<entries; i++) arrayInvariant->AddLast( collection->At(i));
  compInverse->SetCorrections( arrayInvariant);
  compInverse->SetWeights(&weights);
  compInverse->AddVisualCorrection(compInverse,1);
  compInput->AddVisualCorrection(compInput,2);
  TF1 finv1("finv1","AliTPCCorrection::GetCorrXYZ(x,x,100,0,1)",85,245);
  //
  TVectorD vecCompInverse(npointsTest);
  for (Int_t icoord=0; icoord<3; icoord++){
    for (Int_t ipoint=0; ipoint<npointsTest; ipoint++){
      Double_t r= 85.+gRandom->Rndm()*150;
      Double_t phi= gRandom->Rndm()*TMath::TwoPi();
      Double_t z=500*(gRandom->Rndm()-0.5);
    vecCompInverse[ipoint]=AliTPCCorrection::GetCorrXYZ(r*TMath::Cos(phi),r*TMath::Sin(phi),z, icoord, 1);
    }
    Double_t rms=TMath::RMS(npointsTest,vecCompInverse.GetMatrixArray());
    Double_t mean=TMath::Mean(npointsTest,vecCompInverse.GetMatrixArray());
    isOK1[icoord]=TMath::Abs(rms)<kEpsilon;
    isOK1[icoord]&=TMath::Abs(mean)<kEpsilon;
  }
  if (isOK1[0]==kFALSE || isOK1[1]==kFALSE ||isOK1[2]==kFALSE ){
    ::Error("TestCorrection_AddCorrectionCompact",TString::Format("AliTPCComposedCorrection - Test1 (%d,%d,%d) FAILED",isOK1[0], isOK1[1],isOK1[2]).Data());    
  }else{
    ::Info("TestCorrection_AddCorrectionCompact","AliTPCComposedCorrection - Test1  OK");    
  } 
  //
  //  2.) Make compact for of the Composed correction. Test correction  checking invariant inverse x orig
  //      This take time - dostortion has to be recalculated
  AliTPCComposedCorrection *compOutInverseCompact = new AliTPCComposedCorrection();
  compOutInverseCompact->AddCorrectionCompact(compInputFast,1);
  compOutInverseCompact->AddCorrectionCompact(compInputFast,-1);
  compOutInverseCompact->SetOmegaTauT1T2(0,1,1);
  compInputFast->SetOmegaTauT1T2(0,1,1);
  compOutInverseCompact->AddVisualCorrection(compOutInverseCompact,10);  
  compInputFast->AddVisualCorrection(compInput,3);
  TStopwatch timer;
  //
  TF1 fcomp("fcomp","AliTPCCorrection::GetCorrXYZ(x,x,100,0,10)",85,245);
  TF1 forig("forig","-AliTPCCorrection::GetCorrXYZ(x,x,100,0,3)",85,245);
  TF1 fdiff("fdiff","AliTPCCorrection::GetCorrXYZ(x,x,100,0,10)+AliTPCCorrection::GetCorrXYZ(x,x,100,0,2)",85,245);
  timer.Print();


  return kTRUE;
}
