/*
  Unit test for fucntions used in the calibrations
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER\
   -I$ALICE_ROOT/TPC -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/TOF -I$ALICE_ROOT/RAW -I$ALICE_ROOT/PWG1 -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/Base -I$ALICE_ROOT/TPC/Calib");
   
  .L $ALICE_ROOT/TPC/Base/test/UnitTest.C+ 
 
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
// PARAMETERS:
//
TString baseDir="/hera/alice/wiechula/calib/guiTrees";
//
//


void  UnitTestAliTPCCalPadTree(){
  //
  //
  //
  TObjArray *fArray = new TObjArray(100);
  TTree * treePad=AliTPCcalibDButil::ConnectGainTrees(baseDir);
  for (Int_t i=0; i<3; i++){
    AliTPCCalPad * padLx = AliTPCCalPad::MakePadFromTree(treePad,"lx.fElements","QMax");
    AliTPCCalPad * padLy = AliTPCCalPad::MakePadFromTree(treePad,"ly.fElements","QMax");
    AliTPCCalPad * padMax = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MaxCharge.fElements","QMax");
    AliTPCCalPad * padMean = AliTPCCalPad::MakePadFromTree(treePad,"QA.2010.LHC10d.MeanCharge.fElements","QTot");
    if (i>0) {
      padLx->MedianFilter(i,2*i);
      padLy->MedianFilter(i,2*i);
      padMax->MedianFilter(i,2*i);
      padMean->MedianFilter(i,2*i);
    }
    padLx->SetName(TString::Format("QLx%d",i).Data());
    padLy->SetName(TString::Format("QLy%d",i).Data());
    padMax->SetName(TString::Format("QMax%d",i).Data());
    padMean->SetName(TString::Format("QTot%d",i).Data());
    fArray->AddLast(padLx);
    fArray->AddLast(padLy);
    fArray->AddLast(padMax);
    fArray->AddLast(padMean);
  }
  AliTPCCalibViewer::MakeTree("QAtest.root", fArray,0);
  //
  TFile*fout= TFile::Open("QAtest.root");
  TTree * tree  = (TTree*)fout->Get("calPads");
  
}
