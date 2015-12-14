/*
  Unit test for the  AliNDLocalRegression class:

  Tests:
  1.) Check consistency of the error estimates and bias for smaple with gausian noise
  2.) Check persistancy streamer
  3.) Check outlier handling
  .x $NOTES/aux/rootlogon.C
  gSystem->AddIncludePath("-I$ALICE_ROOT/../src/STAT/");
  gSystem->SetIncludePath("-I$ALICE_PHYSICS/include -I$ALICE_ROOT/include");
  .L $ALICE_ROOT/../src/STAT/test/AliNDLocalRegressionTest.C+  
  AliNDLocalRegressionTest(10000,2,"cos(7*x[0]/pi)*sin(11*x[1]/pi)",0.1);
  //
  AliNDLocalRegressionTest(10000,2,"x[0]+x[1]",0.1);

*/


#include "THn.h"
#include "TCanvas.h"
#include "AliNDLocalRegression.h" 
#include "TStyle.h"
#include "TPaveText.h"

void UnitTestGaussNoise();
void UnitTestStreamer();
Bool_t UnitTestContrain();


TTree * treeIn=0;
AliNDLocalRegression *pfitNDIdeal=0;       // ideal fit without noise
AliNDLocalRegression *pfitNDGaus0=0;       // fit with noisy data  - sample 0
AliNDLocalRegression *pfitNDGaus1=0;       // fit with noisy data  - sample 1
AliNDLocalRegression *pfitNDGaus2=0;       // fit with noisy data  - sample 2 - 2 gaussian signal+ backgorund
AliNDLocalRegression *pfitNDGaus3=0;       // fit with noisy data  - sample 2 - 2 gaussian signal+ backgorund

AliNDLocalRegression *pfitNDBreit0=0;       // fit with noisy data  - BreitWigner 0
AliNDLocalRegression *pfitNDBreit1=0;       // fit with noisy data  - BreitWigner 1
AliNDLocalRegression *pfitNDBreit2=0;       // fit with noisy data  - BreitWigner 2 - without outlier removal

THn   *hN=0; 
TFormula *pformula= 0;

void PlotData(TH1F *hData,TString xTitle="xTitle",TString yTitle="yTitle",Color_t color=kBlack,TString zTitle="zTitle",Double_t rms=999999., Double_t eRms=0.,Double_t mean=999999.,Double_t eMean=0.){
  //
  //
  //
  if(color==(kRed+2)){hData->SetMarkerStyle(20);}
  if(color==(kBlue+2)){hData->SetMarkerStyle(21);}
  if(color==(kGreen+2)){hData->SetMarkerStyle(22);hData->SetMarkerSize(1.3);}
  
  hData->SetMarkerColor(color);
  hData->SetLineColor(color);
  hData->GetXaxis()->SetTitle(xTitle.Data());
  hData->GetYaxis()->SetTitle(yTitle.Data());
  hData->GetZaxis()->SetTitle(zTitle.Data());
  hData->GetXaxis()->SetTitleOffset(1.2);
  hData->GetXaxis()->SetTitleSize(0.05);
  hData->GetYaxis()->SetTitleOffset(1.3);
  hData->GetYaxis()->SetTitleSize(0.05);
  hData->GetXaxis()->SetLabelSize(0.035);
  hData->GetYaxis()->SetLabelSize(0.035);
  hData->GetXaxis()->SetDecimals();
  hData->GetYaxis()->SetDecimals();
  hData->GetZaxis()->SetDecimals();
  hData->Sumw2();
  hData->Draw("pe1");
  
  if(mean!=999999.){
    TPaveText *text1 = new TPaveText(0.21,0.82,0.51,0.92,"NDC");
    text1->SetTextFont(43);
    text1->SetTextSize(30.);
    text1->SetBorderSize(1);
    text1->SetFillColor(kWhite);
    text1->AddText(Form("Mean: %0.2f #pm %0.2f",mean,eMean));
    text1->AddText(Form("RMS: %0.2f #pm %0.2f",rms,eRms));
    text1->Draw();
  }
  if(rms!=999999. && mean==999999.){
    TPaveText *text1 = new TPaveText(0.21,0.87,0.51,0.92,"NDC");
    text1->SetTextFont(43);
    text1->SetTextSize(30.);
    text1->SetBorderSize(1);
    text1->SetFillColor(kWhite);
    text1->AddText(Form("RMS: %0.2f",rms));
    text1->Draw();
  }
}

void UnitTestGaussNoise(){
  //
  // Unit Test pulls Gauss
  // Compare 2 Regression objects obtained from 2 indpendent training sample
  //
  // Test:
  //   Bias < 10  *error
  //   Pull < 1+ 3*error
  //   
  TCanvas * canvasUnitGausNoise = new TCanvas("canvasUnitGausNoise","canvasUnitGausNoise",800,800);
  treeIn->Draw("(AliNDLocalRegression::GetCorrND(3,xyz0,xyz1)-AliNDLocalRegression::GetCorrND(2,xyz0,xyz1))/sqrt(AliNDLocalRegression::GetCorrNDError(3,xyz0,xyz1)**2+AliNDLocalRegression::GetCorrNDError(2,xyz0,xyz1)**2)>>pullsGaus32(401,-20.5,20.5)","","");
	TH1F   *pullsGaus32 = (TH1F*)gPad->GetPrimitive("pullsGaus32");
  Double_t meanPullsGaus01 = treeIn->GetHistogram()->GetMean();
  Double_t meanPullsGaus01Err = treeIn->GetHistogram()->GetMeanError();
  Double_t rmsPullsGaus01 = treeIn->GetHistogram()->GetRMS();
  Double_t rmsPullsGaus01Err = treeIn->GetHistogram()->GetRMSError();
  if (TMath::Abs(meanPullsGaus01)<10*meanPullsGaus01Err) {
    ::Info( "AliNDLocalRegressionTest","mean pull OK %3.3f\t+-%3.3f", meanPullsGaus01, meanPullsGaus01Err);
  }else{
    ::Error( "AliNDLocalRegressionTest","mean pull NOT OK %3.3f\t+-%3.3f", meanPullsGaus01, meanPullsGaus01Err);
  }
  if (rmsPullsGaus01<1+rmsPullsGaus01Err) {
    ::Info( "AliNDLocalRegressionTest"," rms pull OK %3.3f\t+-%3.3f", rmsPullsGaus01, rmsPullsGaus01Err);
  }else{
    ::Error( "AliNDLocalRegressionTest"," rms pull NOT OK %3.3f\t+-%3.3f", rmsPullsGaus01, rmsPullsGaus01Err);
  }
  //
  //
  //
	PlotData(pullsGaus32,"Gaus pulls","counts (arb. units)",kRed+2,"zTitle",rmsPullsGaus01,rmsPullsGaus01Err,meanPullsGaus01,meanPullsGaus01Err);
  canvasUnitGausNoise->SaveAs("AliNDLocalRegressionTest.canvasUnitTestGaussNoise.png");
}

void UnitTestGaussNoisePlusOutliers(){
  //
  // Unit Test pulls Gauss
  // Compare 2 Regression objects obtained from 2 indpendent training sample
  //
  // Test:
  //   Bias < 10  *error
  //   Pull < 1+ 3*error
  //   
  TCanvas * canvasUnitGausNoisesOutliers = new TCanvas("canvasUnitGausNoisesOutliers","canvasUnitGausNoisesOutliers",800,800);
  treeIn->Draw("(AliNDLocalRegression::GetCorrND(7,xyz0,xyz1)-AliNDLocalRegression::GetCorrND(2,xyz0,xyz1))/sqrt(AliNDLocalRegression::GetCorrNDError(7,xyz0,xyz1)**2+AliNDLocalRegression::GetCorrNDError(2,xyz0,xyz1)**2)>>pullsGaus72(401,-20.5,20.5)","","");
  TH1F   *pullsGaus72 = (TH1F*)gPad->GetPrimitive("pullsGaus72");
  Double_t meanPullsGaus72 = treeIn->GetHistogram()->GetMean();
  Double_t meanPullsGaus72Err = treeIn->GetHistogram()->GetMeanError();
  Double_t rmsPullsGaus72 = treeIn->GetHistogram()->GetRMS();
  Double_t rmsPullsGaus72Err = treeIn->GetHistogram()->GetRMSError();
  if (TMath::Abs(meanPullsGaus72)<10*meanPullsGaus72Err) {
    ::Info( "AliNDLocalRegressionTest::UnitTestGaussNoisePluOutliers","mean pull OK %3.3f\t+-%3.3f", meanPullsGaus72, meanPullsGaus72Err);
  }else{
    ::Error( "AliNDLocalRegressionTest::UnitTestGaussNoisePluOutliers","mean pull NOT OK %3.3f\t+-%3.3f", meanPullsGaus72, meanPullsGaus72Err);
  }
  if (rmsPullsGaus72<1+rmsPullsGaus72Err) {
    ::Info( "AliNDLocalRegressionTest::UnitTestGaussNoisePluOutliers"," rms pull OK %3.3f\t+-%3.3f", rmsPullsGaus72, rmsPullsGaus72Err);
  }else{
    ::Error( "AliNDLocalRegressionTest::UnitTestGaussNoisePluOutliers"," rms pull NOT OK %3.3f\t+-%3.3f", rmsPullsGaus72, rmsPullsGaus72Err);
  }
  //
  //
  //
  PlotData(pullsGaus72,"Gaus pulls","counts (arb. units)",kRed+2,"zTitle",rmsPullsGaus72,rmsPullsGaus72Err,meanPullsGaus72,meanPullsGaus72Err);
  canvasUnitGausNoisesOutliers->SaveAs("AliNDLocalRegressionTest.UnitTestGaussNoisePluOutliers.png");
}



void UnitTestBreitWignerNoise(){
  //
  // Unit Test pulls BreitWigner
  // Compare 2 Regression objects obtained from 2 independent training sample
  //
  // Test:
  //   Bias < 10  *error
  //   Pull < 1+ 3*error
  //   
  TCanvas * canvasUnitBreitWigner = new TCanvas("canvasUnitBretWigner","canvasUnitBreitWigner",800,800);
  
  treeIn->Draw("(AliNDLocalRegression::GetCorrND(5,xyz0,xyz1)-AliNDLocalRegression::GetCorrND(4,xyz0,xyz1))/sqrt(AliNDLocalRegression::GetCorrNDError(4,xyz0,xyz1)**2+AliNDLocalRegression::GetCorrNDError(5,xyz0,xyz1)**2)>>pullsBreiWigner54(401,-20.5,20.5)","","");
	TH1F   *pullsBreiWigner54 = (TH1F*)gPad->GetPrimitive("pullsBreiWigner54");
  Double_t meanPullsBreitWigner = treeIn->GetHistogram()->GetMean();
  Double_t meanPullsBreitWignerErr = treeIn->GetHistogram()->GetMeanError();
  Double_t rmsPullsBreitWigner = treeIn->GetHistogram()->GetRMS();
  Double_t rmsPullsBreitWignerErr = treeIn->GetHistogram()->GetRMSError();
  if (TMath::Abs(meanPullsBreitWigner)<10*meanPullsBreitWignerErr) {
    ::Info( "AliNDLocalRegressionTest::UnitTestBreitWignerNoise","mean pull OK %3.3f\t+-%3.3f", meanPullsBreitWigner, meanPullsBreitWignerErr);
  }else{
    ::Error( "AliNDLocalRegressionTest::UnitTestBreitWignerNoise","mean pull NOT OK %3.3f\t+-%3.3f", meanPullsBreitWigner, meanPullsBreitWignerErr);
  }
  if (rmsPullsBreitWigner<1+rmsPullsBreitWignerErr) {
    ::Info( "AliNDLocalRegressionTest::UnitTestBreitWignerNoise"," rms pull OK %3.3f\t+-%3.3f", rmsPullsBreitWigner, rmsPullsBreitWignerErr);
  }else{
    ::Error( "AliNDLocalRegressionTest::UnitTestBreitWignerNoise"," rms pull NOT OK %3.3f\t+-%3.3f", rmsPullsBreitWigner, rmsPullsBreitWignerErr);
  }
	PlotData(pullsBreiWigner54,"BreitWigner pulls","counts (arb. units)",kBlue+2,"zTitle",rmsPullsBreitWigner,rmsPullsBreitWignerErr,meanPullsBreitWigner,meanPullsBreitWignerErr);
  canvasUnitBreitWigner->SaveAs("AliNDLocalRegressionTest.canvasUnitBreitWigner.png");
}



void UnitTestStreamer(){
  //
  // Check consistency of streamer
  //
  TCanvas * canvasUnitTestStreamer = new TCanvas("canvasUnitTestStreamer","canvasUnitTestStreamer",800,800);
  TFile * f = TFile::Open("fitNDLocalTestStreamer.root","recreate");
  pfitNDGaus0->Write("pfitNDGaus0");
  f->Close();
  f = TFile::Open("fitNDLocalTestStreamer.root");
  AliNDLocalRegression *streamTest1 = (AliNDLocalRegression *)f->Get("pfitNDGaus0");
  streamTest1->AddVisualCorrection(streamTest1,100);
  Int_t entries = treeIn->Draw("AliNDLocalRegression::GetCorrND(100,xyz0,xyz1)-AliNDLocalRegression::GetCorrND(2,xyz0,xyz1)>>streamerTest(201,-1.05,1.05)","");
  TH1F   *streamerTest = (TH1F*)gPad->GetPrimitive("streamerTest");
  Int_t rms = TMath::RMS(entries,treeIn->GetV1());
  if (rms!=0){
    ::Error( "AliNDLocalRegressionTest::UnitTestStreamer","Streamer problem");
  }else{
    ::Info( "AliNDLocalRegressionTest::UnitTestStreamer","Streamer OK");
  }
  PlotData(streamerTest,"streamer test","counts (arb. units)",kGreen+2,"zTitle",rms);
  canvasUnitTestStreamer->SaveAs(" AliNDLocalRegressionTest.UnitTestStreamer.png");
}



void AliNDLocalRegressionTest(Int_t npoints=10000, Int_t ndim=2, const char *sfromula="cos(7*x[0]/pi)*sin(11*x[1]/pi)", Double_t err=0.1){
  //
  // Local regression test method
  //
  // Int_t npoints=100000; Int_t ndim=2; const char *sfromula="cos(10*x[0])*cos(15*x[1])"; 
  //
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetOptStat(0);
  //
  pformula=new TFormula("pformula", sfromula);
  pfitNDIdeal  = new  AliNDLocalRegression("pfitNDIdeal","pfitNDIdeal");
  pfitNDGaus0  = new  AliNDLocalRegression("pfitNDGaus0","pfitNDGaus0");
  pfitNDGaus1  = new  AliNDLocalRegression("pfitNDGaus1","pfitNDGaus1");
  pfitNDGaus2  = new  AliNDLocalRegression("pfitNDGaus2","pfitNDGaus2");
  pfitNDGaus3  = new  AliNDLocalRegression("pfitNDGaus3","pfitNDGaus3");
  pfitNDBreit0 = new  AliNDLocalRegression("pfitNDBreit0","pfitNDBreit0");
  pfitNDBreit1 = new  AliNDLocalRegression("pfitNDBreit1","pfitNDBreit1");
  pfitNDBreit2 = new  AliNDLocalRegression("pfitNDBreit2","pfitNDBreit2");
  //
  // 0.) Initialization of variables and THn
  // 
  TTreeSRedirector *pcstreamIn         = new TTreeSRedirector("fitNDLocalTestInput.root","recreate");
  TTreeSRedirector *pcstreamOutIdeal   = new TTreeSRedirector("fitNDLocalTestOutputIdeal.root","recreate");  
  TTreeSRedirector *pcstreamOutGaus0   = new TTreeSRedirector("fitNDLocalTestOutputGaus0.root","recreate");  
  TTreeSRedirector *pcstreamOutGaus1   = new TTreeSRedirector("fitNDLocalTestOutputGaus1.root","recreate");  
  TTreeSRedirector *pcstreamOutGaus2   = new TTreeSRedirector("fitNDLocalTestOutputGaus2.root","recreate");  
  TTreeSRedirector *pcstreamOutGaus3   = new TTreeSRedirector("fitNDLocalTestOutputGaus3.root","recreate");  
  TTreeSRedirector *pcstreamOutBreit0  = new TTreeSRedirector("fitNDLocalTestOutputBreit0.root","recreate");  
  TTreeSRedirector *pcstreamOutBreit1  = new TTreeSRedirector("fitNDLocalTestOutputBreit1.root","recreate");  
  TTreeSRedirector *pcstreamOutBreit2  = new TTreeSRedirector("fitNDLocalTestOutputBreit2.root","recreate");  
  
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
  for (Int_t ipoint=0; ipoint<TMath::Abs(npoints); ipoint++){
    for (Int_t idim=0; idim<ndim; idim++){
      xyz[idim]=gRandom->Rndm();
    }
    Double_t value=pformula->EvalPar(xyz,0);
    Double_t noise = gRandom->Gaus()*err;
    Double_t noise2 = noise*(1+(gRandom->Rndm()<0.1)*100);  // noise with 10 percent of outliers
    Double_t noiseBreit = gRandom->BreitWigner()*err;
    (*pcstreamIn)<<"testInput"<<
      "val="<<value<<
      "err="<<err<<
      "noise="<<noise<<            // gausian noise 
      "noise2="<<noise2<<          // gausian noise + 10% of outliers
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
  if (npoints>0){
    pfitNDIdeal->SetStreamer(pcstreamOutIdeal);
    pfitNDGaus0->SetStreamer(pcstreamOutGaus0);
    pfitNDGaus1->SetStreamer(pcstreamOutGaus1);
    pfitNDGaus2->SetStreamer(pcstreamOutGaus2);
    pfitNDGaus3->SetStreamer(pcstreamOutGaus3);
    pfitNDBreit0->SetStreamer(pcstreamOutBreit0);
    pfitNDBreit1->SetStreamer(pcstreamOutBreit1);
    pfitNDBreit2->SetStreamer(pcstreamOutBreit2);
  }
  //
  pfitNDIdeal->SetHistogram((THn*)(hN->Clone()));
  pfitNDGaus0->SetHistogram((THn*)(hN->Clone()));
  pfitNDGaus1->SetHistogram((THn*)(hN->Clone()));
  pfitNDGaus2->SetHistogram((THn*)(hN->Clone()));
  pfitNDGaus3->SetHistogram((THn*)(hN->Clone()));
  pfitNDGaus2->SetCuts(3,0.8,1);
  pfitNDGaus3->SetCuts(0,0.8,0);

  pfitNDBreit0->SetCuts(3,0.8,1);
  pfitNDBreit1->SetCuts(3,0.8,1);
  pfitNDBreit2->SetCuts(0,0.8,1);
  pfitNDBreit0->SetHistogram((THn*)(hN->Clone()));
  pfitNDBreit1->SetHistogram((THn*)(hN->Clone()));
  pfitNDBreit2->SetHistogram((THn*)(hN->Clone()));
  //
  pfitNDIdeal->MakeFit(treeIn, "val:err", "xyz0:xyz1","Entry$%2==1", "0.02:0.02","2:2",0.0001);
  pfitNDGaus0->MakeFit(treeIn, "val+noise:err", "xyz0:xyz1","Entry$%2==0", "0.02:0.02","2:2",0.0001);  // sample Gaussian1
  pfitNDGaus1->MakeFit(treeIn, "val+noise:err", "xyz0:xyz1","Entry$%2==1", "0.02:0.02","2:2",0.0001);  // sample Gaussian2
  pfitNDGaus2->MakeFit(treeIn, "val+noise2:err", "xyz0:xyz1","Entry$%2==1", "0.02:0.02","2:2",0.0001);  // sample Gaussian2 - with tails robust
  pfitNDGaus3->MakeFit(treeIn, "val+noise2:err", "xyz0:xyz1","Entry$%2==1", "0.02:0.02","2:2",0.0001);  // sample Gaussian2 - with tails non robust
  pfitNDBreit0->MakeFit(treeIn, "val+noiseBreit:err", "xyz0:xyz1","Entry$%2==0", "0.02:0.02","2:2",0.0001);  // sample Breit0
  pfitNDBreit1->MakeFit(treeIn, "val+noiseBreit:err", "xyz0:xyz1","Entry$%2==1", "0.02:0.02","2:2",0.0001);  // sample Breit1
  pfitNDBreit2->MakeFit(treeIn, "val+noiseBreit:err", "xyz0:xyz1","Entry$%2==1", "0.02:0.02","2:2",0.0001);  // sample Breit2 without outlier filtering
  //
  if (npoints<0) return;  // callgrind mode of operation
  pfitNDIdeal->AddVisualCorrection(pfitNDIdeal,1);
  pfitNDGaus0->AddVisualCorrection(pfitNDGaus0,2);
  pfitNDGaus1->AddVisualCorrection(pfitNDGaus1,3);
  pfitNDGaus2->AddVisualCorrection(pfitNDGaus2,7);
  pfitNDGaus3->AddVisualCorrection(pfitNDGaus3,8);
  pfitNDBreit0->AddVisualCorrection(pfitNDBreit0,4);
  pfitNDBreit1->AddVisualCorrection(pfitNDBreit1,5);
  pfitNDBreit2->AddVisualCorrection(pfitNDBreit2,6);
  //

  TObjArray * array = AliNDLocalRegression::GetVisualCorrections();
  for (Int_t i=0; i<array->GetEntries(); i++){
    AliNDLocalRegression * regression = ( AliNDLocalRegression *)array->At(i);
    if (regression==NULL) continue;
    regression->AddVisualCorrection(regression);
    Int_t hashIndex = regression->GetVisualCorrectionIndex();
    treeIn->SetAlias( regression->GetName(), TString::Format("AliNDLocalRegression::GetCorrND(%d,xyz0,xyz1+0)",hashIndex).Data());
  }

  
  delete pcstreamOutGaus0;
  delete pcstreamOutGaus1;
  delete pcstreamOutGaus2;
  delete pcstreamOutGaus3;
  delete pcstreamOutBreit0;
  delete pcstreamOutBreit1;
  delete pcstreamOutBreit2;
  UnitTestGaussNoise(); 
  UnitTestGaussNoisePlusOutliers();
  UnitTestBreitWignerNoise();
  UnitTestStreamer();
  UnitTestContrain();
}


Bool_t UnitTestContrain(){ 
  //
  // To test:
  //    1.) Improvement using the constraint with/without
  //    2.) Compare 2 different constraints
  //    3.) Repeating Kalman update data are smoother but the error estimate is underestimated
  //             
  Int_t nDims=2;
  Int_t indexes[2]={0,1};
  Double_t relWeight0[6]={1,1,1,1,1,1};
  Double_t relWeight1[6]={1,1,10,1,1,10};
  TTreeSRedirector*pcstream=new TTreeSRedirector("constrainStream.root","recreate");
  
  AliNDLocalRegression * regression0 = ( AliNDLocalRegression *)AliNDLocalRegression::GetVisualCorrections()->FindObject("pfitNDGaus0");
  AliNDLocalRegression * regression1 = ( AliNDLocalRegression *)AliNDLocalRegression::GetVisualCorrections()->FindObject("pfitNDGaus1");
  AliNDLocalRegression *regressionUpdate0 = (AliNDLocalRegression *)regression0->Clone();
  AliNDLocalRegression *regressionUpdate1 = (AliNDLocalRegression *)regression1->Clone();
 
  for (Int_t iter=0; iter<5; iter++){
    regressionUpdate0->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeight0, pcstream);
    regressionUpdate1->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeight1, pcstream);
  }
  regressionUpdate1->AddWeekConstrainsAtBoundaries(nDims, indexes,relWeight1, pcstream,kTRUE);


  regressionUpdate0->SetName("pfitNDGaus0_Updated");
  regressionUpdate1->SetName("pfitNDGaus1_Updated");
  AliNDLocalRegression::AddVisualCorrection(regressionUpdate0);
  AliNDLocalRegression::AddVisualCorrection(regressionUpdate1);
  treeIn->SetAlias( regressionUpdate0->GetName(), TString::Format("AliNDLocalRegression::GetCorrND(%d,xyz0,xyz1+0)", regressionUpdate0->GetVisualCorrectionIndex()).Data());
  treeIn->SetAlias( regressionUpdate1->GetName(), TString::Format("AliNDLocalRegression::GetCorrND(%d,xyz0,xyz1+0)", regressionUpdate1->GetVisualCorrectionIndex()).Data());
  treeIn->SetAlias( TString::Format("%sErr",regressionUpdate0->GetName()).Data(), TString::Format("AliNDLocalRegression::GetCorrNDError(%d,xyz0,xyz1+0)", regressionUpdate0->GetVisualCorrectionIndex()).Data());
  treeIn->SetAlias( TString::Format("%sError",regressionUpdate1->GetName()).Data(), TString::Format("AliNDLocalRegression::GetCorrNDError(%d,xyz0,xyz1+0)", regressionUpdate1->GetVisualCorrectionIndex()).Data());
  delete pcstream;

  TH1 * his[10]={0};
  TCanvas *canvasUnitTestConstrain = new TCanvas("canvasUnitTestConstrain","canvasUnitTestConstrain");
  treeIn->SetLineColor(4);
  treeIn->Draw("(pfitNDGaus1_Updated-pfitNDIdeal)>>hisSmoothedC(100,-0.2,0.2)","abs(xyz0-0.5)<0.45&&abs(xyz1-0.05)<0.45","");
  his[0]=treeIn->GetHistogram();
  treeIn->SetLineColor(2);
  treeIn->Draw("(pfitNDGaus0_Updated-pfitNDIdeal)>>hisSmoothed(100,-0.2,0.2)","abs(xyz0-0.5)<0.45&&abs(xyz1-0.05)<0.45","");
  his[1]=treeIn->GetHistogram();
  treeIn->SetLineColor(1);
  treeIn->SetLineColor(1);
  treeIn->Draw("(pfitNDGaus0-pfitNDIdeal)>>hisNotSmoothed(100,-0.2,0.2)","abs(xyz0-0.5)<0.45&&abs(xyz1-0.05)<0.45","same");
  his[2]=treeIn->GetHistogram();
  his[0]->Fit("gaus","+");
  his[1]->Fit("gaus","+");
  his[2]->Fit("gaus","+");
  his[0]->Draw();
  his[1]->Draw("same");
  his[2]->Draw("same");
  canvasUnitTestConstrain->SaveAs("canvasUnitTestConstrain.png");
  
}

