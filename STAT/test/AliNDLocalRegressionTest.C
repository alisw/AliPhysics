/*
  Unit test for the  AliNDLocalRegression class:

  Tests:
  1.) Check consistency of the error estimates and bias for smaple with gausian noise
  2.) Check persistancy streamer
  3.) Check outlier handling


  gSystem->AddIncludePath("-I$ALICE_ROOT/../src/STAT/");
  .L $ALICE_ROOT/../src/STAT/AliNDLocalRegression.cxx+
  .L $ALICE_ROOT/../src/STAT/test/AliNDLocalRegressionTest.C+  
  AliNDLocalRegressionTest(5000);
 
*/


#include "THn.h"
#include "AliNDLocalRegression.h"
void UnitTestGaussNoise();
void UnitTestStreamer();



TTree * treeIn=0;
AliNDLocalRegression *pfitNDIdeal=0;       // ideal fit without noise
AliNDLocalRegression *pfitNDGaus0=0;       // fit with noisy data  - sample 0
AliNDLocalRegression *pfitNDGaus1=0;       // fit with noisy data  - sample 1
AliNDLocalRegression *pfitNDBreit0=0;       // fit with noisy data  - BreitWigner 0
AliNDLocalRegression *pfitNDBreit1=0;       // fit with noisy data  - BreitWigner 1

THn   *hN=0; 
TFormula *pformula= 0;


void UnitTestGaussNoise(){
  //
  // Unit Test pulls Gauss
  // Compare 2 Regression objects obtained from 2 indpendent training sample
  //
  // Test:
  //   Bias < 10  *error
  //   Pull < 1+ 3*error
  //   
  treeIn->Draw("(AliNDLocalRegression::GetCorrND(3,xyz0,xyz1)-AliNDLocalRegression::GetCorrND(2,xyz0,xyz1))/sqrt(AliNDLocalRegression::GetCorrNDError(3,xyz0,xyz1)**2+AliNDLocalRegression::GetCorrNDError(2,xyz0,xyz1)**2)>>pullsGaus01(200,-20,20)","","");
  Double_t meanPullsGaus01 = treeIn->GetHistogram()->GetMean();
  Double_t meanPullsGaus01Err = treeIn->GetHistogram()->GetMeanError();
  Double_t rmsPullsGaus01 = treeIn->GetHistogram()->GetRMS();
  Double_t rmsPullsGaus01Err = treeIn->GetHistogram()->GetRMSError();
  if (TMath::Abs(meanPullsGaus01)<10*meanPullsGaus01Err) {
    ::Info( "AliNDLocalRegressionTest","mean pull OK %3.3f\t+-%3.3f", meanPullsGaus01, meanPullsGaus01Err);
  }else{
    ::Error( "AliNDLocalRegressionTest","mean pull OK %3.3f\t+-%3.3f", meanPullsGaus01, meanPullsGaus01Err);
  }
  if (rmsPullsGaus01<1+rmsPullsGaus01Err) {
    ::Info( "AliNDLocalRegressionTest"," rms pull OK %3.3f\t+-%3.3f", rmsPullsGaus01, rmsPullsGaus01Err);
  }
}

void UnitTestStreamer(){
  //
  // Check consistency of streamer
  //
  TFile * f = TFile::Open("fitNDLocalTestStreamer.root","recreate");
  pfitNDGaus0->Write("pfitNDGaus0");
  f->Close();
  f = TFile::Open("fitNDLocalTestStreamer.root");
  AliNDLocalRegression *streamTest1 = (AliNDLocalRegression *)f->Get("pfitNDGaus0");
  streamTest1->AddVisualCorrection(streamTest1,100);
  Int_t entries = treeIn->Draw("AliNDLocalRegression::GetCorrND(100,xyz0,xyz1)-AliNDLocalRegression::GetCorrND(2,xyz0,xyz1)","");
  Int_t rms = TMath::RMS(entries,treeIn->GetV1());
  if (rms!=0){
    ::Error( "AliNDLocalRegressionTest","Streamer problem");
  }else{
    ::Error( "AliNDLocalRegressionTest","Streamer OK");
  }
}



void AliNDLocalRegressionTest(Int_t npoints=10000, Int_t ndim=2, const char *sfromula="cos(7*x[0]/pi)*sin(19*x[1]/pi)", Double_t err=0.01){
  //
  // Local regression test method
  //
  // Int_t npoints=100000; Int_t ndim=2; const char *sfromula="cos(10*x[0])*cos(15*x[1])"; 
  //
  pformula=new TFormula("pformula", sfromula);
  pfitNDIdeal  = new  AliNDLocalRegression;
  pfitNDGaus0  = new  AliNDLocalRegression;
  pfitNDGaus1  = new  AliNDLocalRegression;
  pfitNDBreit0 = new  AliNDLocalRegression;
  pfitNDBreit1 = new  AliNDLocalRegression;
  //
  // 0.) Initialization of variables and THn
  // 
  TTreeSRedirector *pcstreamIn         = new TTreeSRedirector("fitNDLocalTestInput.root","recreate");
  TTreeSRedirector *pcstreamOutIdeal   = new TTreeSRedirector("fitNDLocalTestOutputIdeal.root","recreate");  
  TTreeSRedirector *pcstreamOutGaus0   = new TTreeSRedirector("fitNDLocalTestOutputGaus0.root","recreate");  
  TTreeSRedirector *pcstreamOutGaus1   = new TTreeSRedirector("fitNDLocalTestOutputGaus1.root","recreate");  
  TTreeSRedirector *pcstreamOutBreit0  = new TTreeSRedirector("fitNDLocalTestOutputBreit0.root","recreate");  
  TTreeSRedirector *pcstreamOutBreit1  = new TTreeSRedirector("fitNDLocalTestOutputBreit1.root","recreate");  
  
  Double_t *xyz     = new Double_t[ndim];
  Double_t *sxyz    = new Double_t[ndim];
  Int_t    *nbins   = new Int_t[ndim];
  Double_t *xmin    = new Double_t[ndim];
  Double_t *xmax    = new Double_t[ndim];
  TString **chxyz = new TString*[ndim];
  for (Int_t idim=0; idim<ndim; idim++) {
    chxyz[idim]=new TString(TString::Format("xyz%d=",idim).Data());
    nbins[idim]=40;
    xmin[idim]=0;
    xmax[idim]=1;
  }  
  hN= new THnF("exampleFit","exampleFit", ndim, nbins, xmin,xmax);
  //
  // 1.) generate random input points
  //
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    for (Int_t idim=0; idim<ndim; idim++){
      xyz[idim]=gRandom->Rndm();
    }
    Double_t value=pformula->EvalPar(xyz,0);
    Double_t noise = gRandom->Gaus()*err;
    Double_t noiseBreit = gRandom->BreitWigner()*err;
    (*pcstreamIn)<<"testInput"<<
      "val="<<value<<
      "err="<<err<<
      "noise="<<noise<<
      "noiseBreit="<<noiseBreit;      
    for (Int_t idim=0; idim<ndim; idim++){
      (*pcstreamIn)<<"testInput"<<chxyz[idim]->Data()<<xyz[idim];
    }
    (*pcstreamIn)<<"testInput"<<"\n";
  }
  delete pcstreamIn;
  pcstreamIn  = new TTreeSRedirector("fitNDLocalTestInput.root");
  treeIn= (TTree*)(pcstreamIn->GetFile()->Get("testInput"));
  //   treeIn->Draw("val:xyz0:xyz1>>his(20,0,1,20,0,1)","","profsurf2",10000); // visualization of input
  //
  // 2.) Make local fits 
  //
  pfitNDIdeal->SetStreamer(pcstreamOutIdeal);
  pfitNDGaus0->SetStreamer(pcstreamOutGaus0);
  pfitNDGaus1->SetStreamer(pcstreamOutGaus1);
  pfitNDBreit0->SetStreamer(pcstreamOutBreit0);
  pfitNDBreit1->SetStreamer(pcstreamOutBreit1);
  //
  pfitNDIdeal->SetHistogram((THn*)(hN->Clone()));
  pfitNDGaus0->SetHistogram((THn*)(hN->Clone()));
  pfitNDGaus1->SetHistogram((THn*)(hN->Clone()));
  pfitNDBreit0->SetHistogram((THn*)(hN->Clone()));
  pfitNDBreit1->SetHistogram((THn*)(hN->Clone()));
  //
  pfitNDIdeal->MakeFit(treeIn, "val:err", "xyz0:xyz1","Entry$%2==1", "0.05:0.05","2:2",0.001);
  pfitNDGaus0->MakeFit(treeIn, "val+noise:err", "xyz0:xyz1","Entry$%2==0", "0.05:0.05","2:2",0.001);
  pfitNDGaus1->MakeFit(treeIn, "val+noise:err", "xyz0:xyz1","Entry$%2==1", "0.05:0.05","2:2",0.001);
  pfitNDBreit0->MakeFit(treeIn, "val+noiseBreit:err", "xyz0:xyz1","Entry$%2==0", "0.05:0.05","2:2",0.001);
  pfitNDBreit1->MakeFit(treeIn, "val+noiseBreit:err", "xyz0:xyz1","Entry$%2==1", "0.05:0.05","2:2",0.001);
  //
  pfitNDIdeal->AddVisualCorrection(pfitNDIdeal,1);
  pfitNDGaus0->AddVisualCorrection(pfitNDGaus0,2);
  pfitNDGaus1->AddVisualCorrection(pfitNDGaus1,3);
  pfitNDBreit0->AddVisualCorrection(pfitNDBreit0,4);
  pfitNDBreit1->AddVisualCorrection(pfitNDBreit1,5);

  UnitTestGaussNoise();
  UnitTestStreamer();
}


