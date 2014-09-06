/*
  Unit test for fucntions used in the Base directory:
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER\
   -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW -I$ALICE_ROOT/PWG1 -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/Base -I$ALICE_ROOT/TPC/Calib");
   
  .L $ALICE_ROOT/TPC/Base/test/UnitTest.C+ 
  UnitTestAliTPCCalPadTree();

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

//
// PARAMETERS to set from outside:
//
TString baseDir="/hera/alice/wiechula/calib/guiTrees";
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
    AliTPCCalPad * padLx = AliTPCCalPad::MakePadFromTree(treePad,"lx.fElements","Lx");
    AliTPCCalPad * padLy = AliTPCCalPad::MakePadFromTree(treePad,"ly.fElements","Ly");
    AliTPCCalPad * padLLx = AliTPCCalPad::MakePadFromTree(treePad,"lx.fElements","LLx");
    AliTPCCalPad * padLLy = AliTPCCalPad::MakePadFromTree(treePad,"ly.fElements","LLy");
    AliTPCCalPad * padMax = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MaxCharge.fElements","QMax");
    AliTPCCalPad * padMean = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MeanCharge.fElements","QTot");
    AliTPCCalPad * padMaxL = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MaxCharge.fElements","QMax");
    AliTPCCalPad * padMeanL = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MeanCharge.fElements","QTot");
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
