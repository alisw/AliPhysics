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
#include "AliLog.h"
#include "THn.h"
#include "TRandom.h"


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
  // test streamer of the AliTPCcalibAlign::Streamer
  //   0.) Read old data part
  //   1.) Fill part 
  //   2.) Write part
  //   3.) Read back - consistency check
  //   4.) Destructor check
  //   5.) Memory usage print
  //
  AliLog::SetClassDebugLevel("AliTPCcalibAlign",1);
  AliTPCcalibAlign * align=0;
  Int_t nPoints=1000000;
  //
  //  0.) Read old data part
  //
  TFile *fin= TFile::Open(fname);
  if (fin){
    AliSysInfo::AddStamp("LoadFile");
    align = (AliTPCcalibAlign * )fin->Get("TPCAlign/alignTPC");
    AliSysInfo::AddStamp("LoadAlign");
    fin->Close();
    delete fin;
    if (align->GetClusterDelta(0)==NULL){
      ::Error("UnitTestAliTPCcalibAlignStreamer","Not back compatible class- GetClusterDelta");
      align->MakeResidualHistos();    
    }
    if (align->GetTrackletDelta(0)==NULL){
      ::Error("UnitTestAliTPCcalibAlignStreamer","Not back compatible class- GetTrackletDelta");
      align->MakeResidualHistosTracklet();
    }
  }else{
  }


  //
  // 1.) Fill part test
  //
  for (Int_t ipoint=0; ipoint<nPoints; ipoint++){
    Double_t xxx[10]={0};
    for (Int_t ihis=0; ihis<2; ihis++){
      THn* his = align->GetClusterDelta(ihis);
      for (Int_t iaxis=0; iaxis<his->GetNdimensions(); iaxis++) {
	xxx[iaxis]=his->GetAxis(iaxis)->GetXmin()+gRandom->Rndm()*(his->GetAxis(iaxis)->GetXmax()-his->GetAxis(iaxis)->GetXmin());
      }
      his->Fill(xxx);
    }
    for (Int_t ihis=0; ihis<4; ihis++){
      THnSparse* his = align->GetTrackletDelta(ihis);
      for (Int_t iaxis=0; iaxis<his->GetNdimensions(); iaxis++) {
	xxx[iaxis]=his->GetAxis(iaxis)->GetXmin()+gRandom->Rndm()*(his->GetAxis(iaxis)->GetXmax()-his->GetAxis(iaxis)->GetXmin());
      }
      his->Fill(xxx);
    }
  } 
  AliSysInfo::AddStamp("FillTrees");
  //
  // 2.) Write part
  //
  TFile * fout=new TFile("testAliTPCcalibAlignStreamer.root","recreate");
  AliSysInfo::AddStamp("WriteAlignStart");
  align->Write("alignTPC"); 
  AliSysInfo::AddStamp("WriteAlignEnd");
  fout->ls();
  fout->Close();
  delete fout;
  //
  // 3.) Read back - consistency check
  //
  fin=new TFile("testAliTPCcalibAlignStreamer.root");
  AliTPCcalibAlign * align2 = (AliTPCcalibAlign *)fin->Get("alignTPC");  
  AliSysInfo::AddStamp("ReadAlign2");
  if (align2==NULL){
    ::Fatal("UnitTestAliTPCcalibAlignStreamer","Alignemnt not read");
  }else{
    ::Info("UnitTestAliTPCcalibAlignStreamer","Alignemnt read-OK");
  }
  if (align2->GetClusterDelta(0)==NULL){
    ::Fatal("UnitTestAliTPCcalibAlignStreamer","histogram GetClusterDelta(0) not read");
  }else{
    ::Info("UnitTestAliTPCcalibAlignStreamer","histogram read GetClusterDelta(0) -OK");
  }
  if (align2->GetTrackletDelta(0)==NULL){
    ::Fatal("UnitTestAliTPCcalibAlignStreamer","histogram GetTrackletDelta(0)not read");
  }else{
    ::Info("UnitTestAliTPCcalibAlignStreamer","histogram read GetTrackletDelta(0) -OK");
  }

  if (align2->GetClusterDelta(0)->GetEntries()!=align->GetClusterDelta(0)->GetEntries()){
    ::Fatal("UnitTestAliTPCcalibAlignStreamer","histogram with different entries");
  }else{
    ::Info("UnitTestAliTPCcalibAlignStreamer","histogram cont. GettrackletDelta(0) -OK");
  }
  if (align2->GetTrackletDelta(0)->GetEntries()!=align->GetTrackletDelta(0)->GetEntries()){
    ::Fatal("UnitTestAliTPCcalibAlignStreamer","histogram with different entries");
  }
  //
  // 4.) Destructor check
  //
  delete align2;
  AliSysInfo::AddStamp("deleteAlign2");
  delete align;
  AliSysInfo::AddStamp("deleteAlign");
  //
  // 5.) Memory usage print
  //
  TTree * treeSys =AliSysInfo::MakeTree("syswatch.log");
  treeSys->Scan("sname:deltaVM:VM:pI.fMemResident","","colsize=30:15:15:20");

}

