/*
  Unit test for fucntions used in the calibrations
  .L $ALICE_ROOT/TPC/Calib/test/UnitTest.C+

  
*/

#include "TF1.h"
#include "TMath.h"
#include "TLinearFitter.h"
#include "TFile.h"
#include "AliTPCcalibAlign.h"
#include "AliTPCcalibLaser.h"
#include "AliSysInfo.h"
#include "TTree.h"

const Double_t kAlmost0=1.e-13;

void UnitTestPlaneFitter();
void UnitTestF1Plane();


void UnitTest(){ 
  UnitTestPlaneFitter();
  UnitTestF1Plane();
}

void UnitTestPlaneFitter(){
  //
  // Unit test plane fitter
  //
  TLinearFitter fitter(6,"x[0]++x[1]++x[2]++x[3]++x[4]++x[5]");
  Double_t x[6]={0};
  TVectorD vec(6);
  for (Int_t i=0; i<10000; i++) { 
    x[0]=i; 
    x[1]=i*i; 
    x[2]=i*i*i; 
    x[3]=TMath::Power(i,1/2.); 
    x[4]=TMath::Power(i,1/3.); 
    x[5]=TMath::Power(i,1/4.); 
    fitter.AddPoint(x,i*i*i);
  }
  fitter.Eval();
  fitter.GetParameters(vec);
  vec.Print();
  if (TMath::Abs(vec[2]-1)>kAlmost0){
    ::Fatal("UnitTestPlaneFitter","Wrong value\n");    
  }else{
    ::Info("UnitTestPlaneFitter","OK");
  }
}


void UnitTestF1Plane(){
  //
  // Test multiplane interface
  //
  TF1 f1("f1","x[0]++x[1]++x[2]");
  Double_t xxx[6]={1,2,3,4,5,6}; 
  Double_t par[6]={1,2,3,4,5,6}; 
  f1.EvalPar(xxx,xxx);  //shoul be 14  = 1*1+2*2+3*3
  //
  if (TMath::Abs(f1.EvalPar(xxx,par)-14.)>kAlmost0){
    ::Fatal("UnitTestF1Plane","Wrong value\n");
  }else{
    ::Info("UnitTestF1Plane", Form("OK: Diff=%f\n",f1.EvalPar(xxx,par)-14 ));
  }
}



void UnitTestAliTPCcalibAlignStreamer(const char *fname="/hera/alice/local/benchmark/vAN-20140518/000128503/cpass1/CalibObjects.root"){
  //
  // 0.) ReadPart
  //
  TFile *fin= TFile::Open(fname);
  AliSysInfo::AddStamp("LoadFile");
  AliTPCcalibAlign * align = (AliTPCcalibAlign * )fin->Get("TPCAlign/alignTPC");
  AliSysInfo::AddStamp("LoadAlign");
  AliTPCcalibLaser * laser = (AliTPCcalibLaser * )fin->Get("TPCAlign/laserTPC");
  AliSysInfo::AddStamp("LoadLaser");
  TTree * tree =AliSysInfo::MakeTree("syswatch.log");
  tree->Scan("sname:deltaVM:VM","","colsize=30");
  //
  // 1.) Write part
  //
  TFile * fout=new TFile("testAliTPCcalibAlignStreamer.root","recreate");
  align->Write();
  fout->ls();
  delete align;
  AliSysInfo::AddStamp("deleteAlign");
  delete laser;
  AliSysInfo::AddStamp("deleteLaser");
  delete fout;
  //
  // 2.) Check memory
  //
  tree =AliSysInfo::MakeTree("syswatch.log");
  tree->Scan("sname:deltaVM:VM","","colsize=30");

}

