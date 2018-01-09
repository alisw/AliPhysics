/// \ingroup STAT/test
/// \brief  test of AliTMinuiToolkit class - log likelihood  fits


/*!
* Example usage:
\code
   AliDrawStyle::SetDefaults();
   AliDrawStyle::ApplyStyle("figTemplate");
   .L  $AliRoot_SRC/STAT/test/AliTMinuitToolkitTest.C+
\endcode

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
#include "AliDrawStyle.h"

Int_t fitEntries=400;
const Int_t kNDim=10;
const Int_t kLineColors[5]={1,2,4,3,6};
const Int_t kLineStyle[5]={1,7,10,1,2};
TF1 likeGausCachy("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike,-10,10,2);
TF1 likeAbs("likeAbs", "abs(x)",-10,10);

TTree *inputTree=0;
void GenerateInput();
void Test1D(Int_t bootStrapIter=400);
void AliTMinuitToolkit_TestHistogram(Int_t nIter);   // can be used

/// Run AliTMinuitToolkitTest
/// \param nIter
/// \param nPoints
void AliTMinuitToolkitTest(Int_t nIter, Int_t nPoints=1000){
  AliDrawStyle::SetDefaults();
  AliDrawStyle::ApplyStyle("figTemplate");
  fitEntries=nPoints;
  AliTMinuitToolkit_TestHistogram(100);
  GenerateInput();
  Test1D(nIter);
}


///\brief Generate input data for the test
void GenerateInput(){
  //
  //  0. Generate input data n dimensional array 1,1
  TTreeSRedirector *pcstream = new TTreeSRedirector("AliTMinuiToolkitTestInput.root","recreate");
  TVectorD xxx(kNDim);
  for (Int_t i=0; i<fitEntries; i++){
    for (Int_t iDim=0; iDim<kNDim; iDim++){
      xxx[iDim]=2.*(gRandom->Rndm()-0.5);
    }
    Double_t noiseG=gRandom->Gaus();
    Double_t noiseL=gRandom->Landau();
    if (inputTree==NULL) {
      (*pcstream)<<"data"<<
                 "x.="<<&xxx<<
                 "noiseG="<<noiseG<<
                 "noiseL="<<noiseL<<
                 "\n";
      inputTree=((*pcstream)<<"data").GetTree();
    }else{
      inputTree->Fill();
    }
  }
  delete pcstream;
  TFile *f = TFile::Open("AliTMinuiToolkitTestInput.root");
  inputTree= (TTree*)f->Get("data");
}


/// \brief This test function shows the basic working principles of AliTMinuitToolkit
///  * illustrates how a robust fit can improve the results for histogram fitting
///  * low level interface shown
///    * user create fitter
///    * define log likelihood function
///  * Algorithm:
///   * 0.)  define input:
///     * 0.1) create example histogram with outliers
///     * 0.1) create fit functions
///   * 1.) Make fits
///     * 1.1) Example fits with chi2 minimization
///     * 1.2) Use "robust" custom user defined log likelihood function e.g gauss+background cachy
///     * 1.3) Use predefined huber cost function
///   * 2.) Use different fit strategies to avoid local minimum
///
///  \image html AliTMinuitToolkit_TestHistogram.png
///  TODO - add exponential fitter to the list of predefined fitters
/// \param nIter
void AliTMinuitToolkit_TestHistogram(Int_t nIter) {
  //
  gStyle->SetOptStat(0);
  TCanvas * canvas = new TCanvas("AliTMinuitToolkit_TestHistogram","AliTMinuitToolkit_TestHistogram",1200,1000);
  canvas->Divide(3,3,0,0);
  TMatrixD initParam(2, 4); // param,error,min,max
  initParam(0, 0) = 20000;
  initParam(0, 1) = 100;
  initParam(0, 2) = 0;
  initParam(0, 3) = 100000;
  initParam(1, 0) = 1;
  initParam(1, 1) = 1;
  initParam(1, 2) = 0;
  initParam(1, 3) = 10;
  //
  TTreeSRedirector *pcstream=new TTreeSRedirector("AliTMinuitToolkit_TestHistogram.root","recreate");
  std::map<string,TF1*> fitMap;
  TF1 *likeGausCachy = new TF1("likeGausCachy", AliTMinuitToolkit::GaussCachyLogLike, -10, 10, 2);
  likeGausCachy->SetParameters(0.9, 1);
  TF1 *likePseudoHuber = new TF1("likePseudoHuber", AliTMinuitToolkit::PseudoHuberLogLike, -10, 10, 2);
  likePseudoHuber->SetParameter(0,3);

  for (Int_t iter=0; iter<nIter; iter++) {
    Double_t slope=1+gRandom->Gaus(0,0.1);
    if (iter <= 9) canvas->cd(iter+1);
    // 0.1 provide some example histogram
    TH1F *hist = new TH1F("test", "AliTMinuitToolkit: Test histogram fit (e^{-x}) with outliers", 50, 0, 4);
    hist->SetMarkerStyle(25);
    TRandom *random = new TRandom();
    for (Int_t i = 0; i < 20000; i++) {
      hist->Fill(random->Exp(slope));
    }
    for (Int_t iOutlier = 0; iOutlier < 10; iOutlier++) {
      Double_t position = gRandom->Rndm() * 4;
      Double_t value = gRandom->Rndm() * hist->GetEntries();
      hist->Fill(position, value);
    }
    hist->SetMinimum(0);hist->SetMaximum(2500);
    hist->Draw("error");
    hist->Draw("same LHist");

    // 0.2 declare fit functions
    TF1 *finput = new TF1("funFit1", "[0]*TMath::Exp(-[1]*x)", 0, 6);
    TF1 *funFit1 = new TF1("funFit1", "[0]*TMath::Exp(-[1]*x)", 0, 6);
    TF1 *funFit2 = new TF1("funFit2", "[0]*TMath::Exp(-[1]*x)", 0, 6);
    TF1 *funFit3 = new TF1("funFit3", "[0]*TMath::Exp(-[1]*x)", 0, 6);
    TVectorD oParam(2);
    // 1.1) example fit without robust option
    AliTMinuitToolkit *fitter = new AliTMinuitToolkit();
    TF1 *aFormExp = new TF1("formExp", "[0]*TMath::Exp(-[1]*x)");
    fitter->SetFitFunction(aFormExp, 0);
    fitter->SetInitialParam(&initParam);
    fitter->FitHistogram(hist);
    funFit1->SetLineColor(1);
    funFit1->SetParameters(fitter->GetParameters()->GetMatrixArray());
    funFit1->Draw("same");
    fitMap["chi2"]=(TF1*)fitter->GetFormula()->Clone();
    // 1.2) Use "robust" custom user defined log likelihood function e.g gauss+background cachy
    fitter->SetLogLikelihoodFunction(likeGausCachy);
    fitter->SetInitialParam(&initParam);
    fitter->FitHistogram(hist);
    funFit2->SetLineColor(2);
    funFit2->SetParameters(fitter->GetParameters()->GetMatrixArray());
    funFit2->Draw("same");
    // 1.3) Use predefined huber cost function
    fitter->SetInitialParam(&initParam);
    fitter->SetLogLikelihoodFunction(NULL);
    fitter->EnableRobust(true);
    fitter->FitHistogram(hist);
    funFit3->SetLineColor(4);
    funFit3->SetParameters(fitter->GetParameters()->GetMatrixArray());
    funFit3->Draw("same");
    fitMap["huber"]=(TF1*)fitter->GetFormula()->Clone();
    //
    // 2.) Test different fit strategies
    fitter->SetLogLikelihoodFunction(likePseudoHuber);
    AliTMinuitToolkit::SetPredefinedFitter("ExpFit", fitter);
    AliTMinuitToolkit::Fit(hist, "ExpFit", "", NULL, "funOption(2,2,1)");
    fitMap["pseudoHuber"]=(TF1*)fitter->GetFormula()->Clone();
    AliTMinuitToolkit::Fit(hist, "ExpFit", "misac(10,50)", NULL, "funOption(2,2,1)");
    fitMap["misacH(10,50)"]=(TF1*)fitter->GetFormula()->Clone();
    AliTMinuitToolkit::Fit(hist, "ExpFit", "misac(10,100)", NULL, "funOption(4,2,2)");
    fitMap["misacH(10,100)"]=(TF1*)fitter->GetFormula()->Clone();
    AliTMinuitToolkit::Fit(hist, "ExpFit", "bootstrap50", NULL, "funOption(6,2,3)");
    fitMap["bootstrapH50"]=(TF1*)fitter->GetFormula()->Clone();
    //
    fitter->SetLogLikelihoodFunction(likeGausCachy);
    AliTMinuitToolkit::Fit(hist, "ExpFit", "misac(10,50)", NULL, "funOption(2,2,1)");
    fitMap["misacGC(10,50)"]=(TF1*)fitter->GetFormula()->Clone();
    AliTMinuitToolkit::Fit(hist, "ExpFit", "misac(10,100)", NULL, "funOption(4,2,2)");
    fitMap["misacGC(10,100)"]=(TF1*)fitter->GetFormula()->Clone();
    AliTMinuitToolkit::Fit(hist, "ExpFit", "bootstrap50", NULL, "funOption(6,2,3)");
    fitMap["bootstrapGC50"]=(TF1*)fitter->GetFormula()->Clone();
    //
    //
    // 2. Draw and save results
    TLegend *legend = new TLegend(0.4, 0.7, 0.89, 0.89, "Test AliTMinuitToolkit - Exponential fit with outliers");
    legend->SetBorderSize(0);
    legend->AddEntry(hist, "Histogram");
    legend->AddEntry(funFit1, "Default chi2 minimization");
    legend->AddEntry(funFit2, "User defined likelihood (0.95*gaus+0.05*cachy)");
    legend->AddEntry(funFit3, "Huber likelihood");
    legend->Draw();
    //
    (*pcstream)<<"test"<<
               "slope="<<slope<<
               "his.="<<hist<<
               "chi2.="<<fitMap["chi2"]<<
               "huber.="<<fitMap["huber"]<<
               "pseudoHuber.="<<fitMap["pseudoHuber"]<<
               "misacH1050.="<<fitMap["misacH(10,50)"]<<
               "misacH10100.="<<fitMap["misacH(10,100)"]<<
               "bootstrapH50.="<<fitMap["bootstrapH50"]<<
               "misacGC1050.="<<fitMap["misacGC(10,50)"]<<
               "misacGC10100.="<<fitMap["misacGC(10,100)"]<<
               "bootstrapGC50.="<<fitMap["bootstrapGC50"]<<

              "\n";
  }
  canvas->SaveAs("AliTMinuitToolkit_TestHistogram.png");
  delete pcstream;
}




void Test1D(Int_t bootStrapIter){
  //
  // 1D fit example
  //    chi2, huber norm, abs(norm) and gaus+cachy log likelihood function
  // Input: 
  //    y=p[0]+p[1]*x;   // p[0]=0; p[1]=2 
  //
  TVectorD oParam(2); oParam[0]=0; oParam[1]=2;
  TMatrixD initParam(2,4); // param,error,min,max
  initParam(0,0)=1; initParam(0,1)=1; initParam(0,2)=0; initParam(0,3)=100000;
  initParam(1,0)=1; initParam(1,1)=1; initParam(1,2)=0; initParam(1,3)=20;

  TF1 *fitFunctions[5];
  TH1 *resHistograms[5];
  TF1 formula1D("formula1D","[0]+[1]*x[0]",-1,1);
  inputTree->SetAlias("test1D","2*x.fElements[0]");
  inputTree->SetAlias("noise","(rndm<0.7)?noiseG:4*noiseL");
  inputTree->SetAlias("X0","x.fElements[0]");

  TString  selection="1";
  AliTMinuitToolkit * tool1D = new AliTMinuitToolkit("AliTMinuitToolkitTest1D.root");
  tool1D->SetVerbose(0x1);
  tool1D->SetFitFunction(&formula1D,kTRUE);
  tool1D->SetInitialParam(&initParam);
  tool1D->FillFitter(inputTree,"test1D+noise:1/sqrt(12.+0)","X0", "", 0,fitEntries);
  //
  formula1D.SetParameters(oParam.GetMatrixArray());
  fitFunctions[0]= (TF1*)formula1D.DrawClone("same");
 
  // 1.1)  Standard fit
  tool1D->EnableRobust(kFALSE); 
  tool1D->SetInitialParam(&initParam);
  tool1D->Fit(""); 
  inputTree->SetAlias("fitChi2Norm",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[1]= (TF1*)formula1D.DrawClone("same");
  tool1D->Bootstrap(bootStrapIter,"bootstrapChi2Norm");
  tool1D->TwoFoldCrossValidation(bootStrapIter,"cvChi2Norm");
  // 1.2)   Huber cost function 
  tool1D->EnableRobust(kTRUE); 
  tool1D->SetInitialParam(&initParam);
  tool1D->Fit("");
  inputTree->SetAlias("fitHuberNorm",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[2]= (TF1*)formula1D.DrawClone("same");
  tool1D->Bootstrap(bootStrapIter,"bootstrapHuberNorm");
  tool1D->TwoFoldCrossValidation(bootStrapIter,"cvHuberNorm");
  // 1.3)  User cost function
 tool1D->SetLogLikelihoodFunction(&likeAbs);
  tool1D->Fit("");
  inputTree->SetAlias("fitLikeAbs",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[3]= (TF1*)formula1D.DrawClone("same");
  tool1D->Bootstrap(bootStrapIter,"bootstrapLikeAbs");
  tool1D->TwoFoldCrossValidation(bootStrapIter,"cvLikeAbs");
  // 1.3)  User cost function
  likeGausCachy.SetParameters(0.8,1);
  tool1D->SetLogLikelihoodFunction(&likeGausCachy);
  tool1D->Fit("");
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
  TGraphErrors *gr0 = TStatToolkit::MakeGraphErrors(inputTree,"test1D+noise:X0:1","",25,1,0.5,0,TMath::Min(fitEntries,100));
  gr0->SetMaximum(30); gr0->SetMinimum(-30);
  TGraphErrors *gr1 = TStatToolkit::MakeGraphErrors(inputTree,"test1D+noise:X0:1","",25,1,0.5,0,TMath::Min(fitEntries,100));
  gr1->SetMaximum(5); gr1->SetMinimum(-5);

  gr0->Draw("ap");
  for (Int_t iFit=0; iFit<5; iFit++){
    fitFunctions[iFit]->SetLineColor(kLineColors[iFit]);
    fitFunctions[iFit]->SetLineStyle(kLineStyle[iFit]);
    fitFunctions[iFit]->SetLineWidth(3);
    fitFunctions[iFit]->Draw("same");
  }
  latex.SetTextSize(0.07);
  latex.DrawLatexNDC(0.11,0.8,"Input y=0+2*x+#epsilon ");
  latex.DrawLatexNDC(0.11,0.7,"      #epsilon=(0.8 Gaus +0.2 Landau)");
  latex.SetTextSize(0.055);
  canvasTest1D->cd(2);

  gr1->Draw("ap");
  for (Int_t iFit=0; iFit<5; iFit++){
    fitFunctions[iFit]->SetLineColor(kLineColors[iFit]);
    fitFunctions[iFit]->SetLineStyle(kLineStyle[iFit]);
    fitFunctions[iFit]->SetLineWidth(3);
    fitFunctions[iFit]->Draw("same");
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
  //for (Int_t iFit=0; iFit<3; iFit++) resHistograms[iFit]->Fit("gaus");
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
  for (Int_t iFit=0; iFit<4; iFit++){
    pad->cd(1);
    latex.SetTextColor(kLineColors[iFit+1]);
    latex.DrawLatexNDC(0.35,0.85-0.05*(iFit+1),TString::Format("%s:  RMS %.2f ",h2D->GetXaxis()->GetBinLabel(1+iFit), grRMS0->GetY()[iFit]).Data());
    pad->cd(2);
    latex.SetTextColor(kLineColors[iFit+1]);
    latex.DrawLatexNDC(0.35,0.85-0.05*(iFit+1),TString::Format("%s:  RMS %.2f ",h2D->GetXaxis()->GetBinLabel(1+iFit), grRMS1->GetY()[iFit]).Data());
  }
  canvasTest1D->SaveAs("AliTMinuitToolkitTest.Test1D.png");
  canvasTest1D->SaveAs("AliTMinuitToolkitTest.Test1D.pdf");

}

