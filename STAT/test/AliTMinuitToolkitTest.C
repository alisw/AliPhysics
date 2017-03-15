/*
   .L  $ALICE_ROOT/../src/STAT/test/AliTMinuitToolkitTest.C+
   Demonstrate performance of the AliTMinuitToolkitTest.C
   later also alarms base on invaraints should be implemented

*/

#include "AliTMinuitToolkit.h"
#include "TTree.h"
#include "TTreeStream.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "AliSysInfo.h"
#include "TGraphErrors.h"
#include "TStatToolkit.h"

Int_t fitEntries=400;
const Int_t kLineColors[5]={1,2,4,3,6};
const Int_t kLineStyle[5]={1,7,10,1,2};
TF1 likeGausCachy("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike,-10,10,2);
TF1 likeAbs("likeAbs", "abs(x)",-10,10);

TTree *inputTree=0;
void GenerateInput();
void Test1D(Int_t bootStrapIter=400);
void TestHistogram();

void AliTMinuitToolkitTest(Int_t iters, Int_t nPoints=1000){
  gStyle->SetLabelSize(0.065,"XYZ");
  gStyle->SetTitleSize(0.065,"XYZ");
  gStyle->SetLabelSize(0.07,"Y");
  gStyle->SetTitleSize(0.07,"Y");

  fitEntries=nPoints;
  TestHistogram();
  GenerateInput();
  Test1D(iters);
}

void GenerateInput(){
  //
  //  0. Generate input data
  TTreeSRedirector *pcstream = new TTreeSRedirector("AliTMinuiToolkitTestInput.root","recreate");
  for (Int_t i=0; i<fitEntries; i++){
    Double_t x0=2*(gRandom->Rndm()-0.5), x1=2*(gRandom->Rndm()-0.5),x2=2*(gRandom->Rndm()-0.5);
    if (inputTree==NULL) {
      (*pcstream)<<"data"<<"testx0="<<x0<<"testx1="<<x1<<"testx2="<<x2<<"\n";
      inputTree=((*pcstream)<<"data").GetTree();
    }else{
      inputTree->Fill();
    }
  }
}

 
void TestHistogram() {
  //
  // This test function shows the basic working principles of this class 
  // and illustrates how a robust fit can improve the results
  //
  
  // 1. provide some example histogram
  gStyle->SetOptStat(0);
  TH1F * hist = new TH1F("test", "AliTMinuitToolkit: Test histogram fit (e^{-x}) with outliers", 20,0,4);
  hist->SetMarkerStyle(25);
  TRandom * rand = new TRandom();
  for (Int_t i = 0; i < 10000; i++) {
    hist->Fill(rand->Exp(1));
    if (i < 1000) hist->Fill(3); //"outliers"
    if (i < 1070) hist->Fill(3.5);
    if (i < 670) hist->Fill(2);
    if (i < 770) hist->Fill(1.5);//"outliers"
    if (i < 740) hist->Fill(1);
  }
  TCanvas * canv = new TCanvas();
  canv->cd(1);hist->Draw("error"); ;hist->Draw("same LHist"); 
  // declare fit functions
  TF1 *finput= new TF1("ffit1", "[0]*TMath::Exp(-[1]*x)", 0, 6);
  TF1 *ffit1 = new TF1("ffit1", "[0]*TMath::Exp(-[1]*x)", 0, 6);
  TF1 *ffit2 = new TF1("ffit2", "[0]*TMath::Exp(-[1]*x)", 0, 6);
  TF1 *ffit3 = new TF1("ffit3", "[0]*TMath::Exp(-[1]*x)", 0, 6); 
  TVectorD oParam(2);
  TMatrixD initParam(2,4); // param,error,min,max
  initParam(0,0)=1; initParam(0,1)=1; initParam(0,2)=0; initParam(0,3)=100000;
  initParam(1,0)=1; initParam(1,1)=1; initParam(1,2)=0; initParam(1,3)=10;

  // 1.) example fit without robust option
  AliTMinuitToolkit * tool = new AliTMinuitToolkit();
  TFormula *aFormExp = new TFormula("formExp", "[0]*TMath::Exp(-[1]*x)");
  tool->SetFitFunction(aFormExp,0);
  tool->SetInitialParam(&initParam);
  tool->FitHistogram(hist);
  ffit1->SetLineColor(1); ffit1->SetParameters(tool->GetParameters()->GetMatrixArray()); ffit1->Draw("same");
  
  // 2.) Use "robust" custom user defined log likelihood function e.g gauss+bckg. cachy
  TF1 * fcost = new TF1("1","abs(x)<10?-log(0.95*exp(-x**2)+0.05/(1+x**2)):-log(0.05/(1+x**2))",-20,20); // 95 % gaus + 5% cachy
  tool->SetLogLikelihoodFunction(fcost);
  tool->FitHistogram(hist);
  ffit2->SetLineColor(2); ffit2->SetParameters(tool->GetParameters()->GetMatrixArray());  ffit2->Draw("same");
  
  // 3.) Use predefined huber cost function 
  tool->SetInitialParam(&initParam); 
  tool->SetLogLikelihoodFunction(NULL);
  tool->EnableRobust(true);
  tool->FitHistogram(hist);
  ffit3->SetLineColor(4); ffit3->SetParameters(tool->GetParameters()->GetMatrixArray());  ffit3->Draw("same");
  // 
  TLegend *legend = new TLegend(0.4,0.7,0.89,0.89,"Test AliTMinuitTolkit - Expontetial fit with outliers");
  legend->SetBorderSize(0);
  legend->AddEntry(hist,"Histogram");
  legend->AddEntry(ffit1,"Default chi2 minimization");
  legend->AddEntry(ffit2,"User defined likelihood (0.95*gaus+0.05*cachy)");
  legend->AddEntry(ffit3,"Huber likelihood");
  legend->Draw();
}




void Test1D(Int_t bootStrapIter){
  //
  // 1D fit example
  TVectorD oParam(2); oParam[0]=0; oParam[1]=5;
  TMatrixD initParam(2,4); // param,error,min,max
  initParam(0,0)=1; initParam(0,1)=1; initParam(0,0)=0; initParam(0,1)=100000;
  initParam(1,0)=1; initParam(1,1)=1; initParam(1,0)=0; initParam(1,1)=20;

  TF1 *fitFunctions[5];
  TH1 *resHistograms[5];
  TF1 formula1D("formula1D","[0]+[1]*x[0]",-1,1);
  inputTree->SetAlias("test1D","5*testx0");
  inputTree->SetAlias("noise","(rndm<0.7)?AliTMinuitToolkit::RrndmGaus():4*AliTMinuitToolkit::RrndmLandau()");
  TString  selection="1";
  AliTMinuitToolkit * tool1D = new AliTMinuitToolkit("AliTMinutiTookitTest1D.root");
  tool1D->SetVerbose(0x1);
  tool1D->SetFitFunction(&formula1D,kTRUE);
  tool1D->SetInitialParam(&initParam);
  tool1D->FillFitter(inputTree,"test1D+noise:1/sqrt(12.+0)","testx0", "", 0,fitEntries);
  //
  formula1D.SetParameters(oParam.GetMatrixArray());
  fitFunctions[0]= (TF1*)formula1D.DrawClone("same");
 
  // 1.1)  Standard fit
  tool1D->EnableRobust(kFALSE); 
  tool1D->SetInitialParam(&initParam);
  tool1D->Fit(kTRUE); 
  inputTree->SetAlias("fitChi2Norm",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[1]= (TF1*)formula1D.DrawClone("same");
  tool1D->Bootstrap(bootStrapIter,"bootstrapChi2Norm");
  tool1D->TwoFoldCrossValidation(bootStrapIter,"cvChi2Norm");
  // 1.2)   Huber cost function 
  tool1D->EnableRobust(kTRUE); 
  tool1D->SetInitialParam(&initParam);
  tool1D->Fit(kTRUE);
  inputTree->SetAlias("fitHuberNorm",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[2]= (TF1*)formula1D.DrawClone("same");
  tool1D->Bootstrap(bootStrapIter,"bootstrapHuberNorm");
  tool1D->TwoFoldCrossValidation(bootStrapIter,"cvHuberNorm");
  // 1.3)  User cost function
 tool1D->SetLogLikelihoodFunction(&likeAbs);
  tool1D->Fit(kTRUE);
  inputTree->SetAlias("fitLikeAbs",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[3]= (TF1*)formula1D.DrawClone("same");
  tool1D->Bootstrap(bootStrapIter,"bootstrapLikeAbs");
  tool1D->TwoFoldCrossValidation(bootStrapIter,"cvLikeAbs");
  // 1.3)  User cost function
  likeGausCachy.SetParameters(0.8,1);
  tool1D->SetLogLikelihoodFunction(&likeGausCachy);
  tool1D->Fit();
  inputTree->SetAlias("fitUserGausAndCachy",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[4]= (TF1*)formula1D.DrawClone("same");
  tool1D->Bootstrap(bootStrapIter,"bootstrapUserGausAndCachy");
  tool1D->TwoFoldCrossValidation(bootStrapIter,"cvUserGausAndCachy");
  delete tool1D;
  //
  // Draw Results
  //
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *canvasTest1D = new TCanvas("Test1D","Test1D",1400,1000);
  canvasTest1D->Divide(1,3);
  TLatex latex;
  latex.SetTextSize(0.055);

  canvasTest1D->cd(1);
  TGraphErrors *gr0 = TStatToolkit::MakeGraphErrors(inputTree,"test1D+noise:testx0:1","",25,1,0.5,0,TMath::Min(fitEntries,100));
  gr0->SetMaximum(30); gr0->SetMinimum(-30);
  TGraphErrors *gr1 = TStatToolkit::MakeGraphErrors(inputTree,"test1D+noise:testx0:1","",25,1,0.5,0,TMath::Min(fitEntries,100));
  gr1->SetMaximum(5); gr1->SetMinimum(-5);

  gr0->Draw("ap");
  for (Int_t ifit=0; ifit<5; ifit++){
    fitFunctions[ifit]->SetLineColor(kLineColors[ifit]);
    fitFunctions[ifit]->SetLineStyle(kLineStyle[ifit]);
    fitFunctions[ifit]->SetLineWidth(3);
    fitFunctions[ifit]->Draw("same");     
  }
  latex.SetTextSize(0.07);
  latex.DrawLatexNDC(0.11,0.8,"Input y=0+5*x+#epsilon ");
  latex.DrawLatexNDC(0.11,0.7,"      #epsilon=(0.8 Gaus +0.2 Landau)");
  latex.SetTextSize(0.055);
  canvasTest1D->cd(2);

  gr1->Draw("ap");
  for (Int_t ifit=0; ifit<5; ifit++){
    fitFunctions[ifit]->SetLineColor(kLineColors[ifit]);
    fitFunctions[ifit]->SetLineStyle(kLineStyle[ifit]);
    fitFunctions[ifit]->SetLineWidth(3);
    fitFunctions[ifit]->Draw("same");     
  }
  TLegend * legend = new TLegend(0.11,0.7,0.6,0.89,"Unbinned1D fit for different cost function (ZOOM)");
  legend->SetBorderSize(0); legend->SetNColumns(2);  
  legend->AddEntry(gr0, "Input data+noise","pl");
  legend->AddEntry(fitFunctions[0],"MC input","l");
  legend->AddEntry(fitFunctions[1],"Fit: Chi2 cost","l");
  legend->AddEntry(fitFunctions[2],"Fit: Huber norm","l");
  legend->AddEntry(fitFunctions[3],"Fit: LogLikeAbs","l");
  legend->AddEntry(fitFunctions[4],"Fit: Log(80%Gaus+20%Cauchy)","l");
  legend->Draw();
  //
  //for (Int_t ifit=0; ifit<3; ifit++) resHistograms[ifit]->Fit("gaus");
  // TTre
  TVirtualPad *pad =canvasTest1D->cd(3);
  gStyle->SetOptTitle();
  TFile *fReport = TFile::Open("AliTMinutiTookitTest1D.root");
  TTree *tBootstrap = (TTree*)fReport->Get("bootstrap");
  TTree *tCrossValidation = (TTree*)fReport->Get("crossValidation");  
  //
  pad->Divide(2,1);
  TH2 *h2D=0;
  pad->cd(1);
  tCrossValidation->SetAlias("dP0","(param0.fElements[0]-param1.fElements[0])");
  tCrossValidation->SetAlias("dP1","(param0.fElements[1]-param1.fElements[1])");
  tCrossValidation->Draw("dP0:reportName.GetName()>>hCrosValidDPar0(4,0,4,200,-1,1)","","colz");
  h2D=(TH2*)(tCrossValidation->GetHistogram());
  TGraph *grRMS0=TStatToolkit::MakeStat1D(h2D, 0, 1., 1, 21, 2);
  grRMS0->Draw("p");
  pad->cd(2);
  tCrossValidation->Draw("dP1:reportName.GetName()>>hCrosValidDPar1(4,0,4,200,-1,1)","","colz");
  h2D=(TH2*)(tCrossValidation->GetHistogram());
  TGraph *grRMS1=TStatToolkit::MakeStat1D(h2D, 0, 1., 1, 21, 2);
  grRMS1->Draw("p");
  pad->cd(1);
  latex.DrawLatexNDC(0.3,0.85,"Two fold cross validation. Param0");
  pad->cd(2);
  latex.DrawLatexNDC(0.3,0.85,"Two fold cross validation. Param1");
  for (Int_t ifit=0; ifit<4; ifit++){ 
    pad->cd(1);
    latex.SetTextColor(kLineColors[ifit+1]);
    latex.DrawLatexNDC(0.35,0.85-0.05*(ifit+1),TString::Format("%s:  RMS %.2f ",h2D->GetXaxis()->GetBinLabel(1+ifit), grRMS0->GetY()[ifit]).Data());
    pad->cd(2);
    latex.SetTextColor(kLineColors[ifit+1]);
    latex.DrawLatexNDC(0.35,0.85-0.05*(ifit+1),TString::Format("%s:  RMS %.2f ",h2D->GetXaxis()->GetBinLabel(1+ifit), grRMS1->GetY()[ifit]).Data());
  }
  canvasTest1D->SaveAs("AliTMinuitToolkitTest.Test1D.png");
  canvasTest1D->SaveAs("AliTMinuitToolkitTest.Test1D.pdf");

}
