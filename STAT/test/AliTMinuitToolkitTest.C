/*
   .x $ALICE_ROOT/../src/STAT/test/AliTMinuitToolkitTest.C+
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


const Int_t fitEntries=20000;
const Int_t kLineColors[4]={1,2,4,6};
const Int_t kLineStyle[4]={1,7,10,3};
TTree *inputTree=0;
void GenerateInput();
void Test1D();


void AliTMinuitToolkitTest(){
  GenerateInput();
  Test1D();
}

void GenerateInput(){
  //
  //  0. Generate input data
  TTreeSRedirector *pcstream = new TTreeSRedirector("AliTMinuiToolkitTest0.root","recreate");
  for (Int_t i=0; i<fitEntries; i++){
    Double_t x0=gRandom->Rndm(), x1=gRandom->Rndm(),x2=gRandom->Rndm();
    if (inputTree==NULL) {
      (*pcstream)<<"data"<<"testx0="<<x0<<"testx1="<<x1<<"testx2="<<x2<<"\n";
      inputTree=((*pcstream)<<"data").GetTree();
    }else{
      inputTree->Fill();
    }
  }
}




void Test1D(){
  //
  // 1D fit example
  TVectorD vParam(2); vParam[0]=0; vParam[1]=0;
  TF1 *fitFunctions[3];
  TH1 *resHistograms[3];
  TF1 formula1D("formula1D","[0]+[1]*x[0]");
  inputTree->SetAlias("test1D","2*testx0");
  inputTree->SetAlias("noise","AliTMinuitToolkit::RrndmGaus()+(rndm<0.2)*AliTMinuitToolkit::RrndmLandau()");
  TString  selection="1";
  AliTMinuitToolkit * tool1D = new AliTMinuitToolkit();
  tool1D->SetFitFunction(&formula1D,kTRUE);
  tool1D->SetInitialParam(&vParam);
  tool1D->FillFitter(inputTree,"test1D+noise:1/sqrt(12.+0)","testx0", "", 0,fitEntries);
 
  // 1.1)  Standard fit
  tool1D->EnableRobust(kFALSE); 
  tool1D->SetInitialParam(&vParam);
  tool1D->Fit();
  inputTree->SetAlias("fitChi2Norm",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[0]= (TF1*)formula1D.DrawClone("same");

  // 1.2)   Huber cost function 
  tool1D->EnableRobust(kTRUE); 
  tool1D->SetInitialParam(&vParam);
  tool1D->Fit();
  inputTree->SetAlias("fitHuberNorm",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[1]= (TF1*)formula1D.DrawClone("same");

  // 1.3)  User cost function
  TF1 * fcost = new TF1("1","abs(x)<10?-log(0.8*exp(-x**2)+0.2/(1+x**2)):-log(0.2/(1+x**2))",-20,20); // 80 % gaus + 20% cachy
  tool1D->SetLogLikelihoodFunction(fcost);
  tool1D->Fit();
  inputTree->SetAlias("fitUserGausAndCachy",tool1D->GetFitFunctionAsAlias().Data());
  formula1D.SetParameters(tool1D->GetParameters()->GetMatrixArray());
  fitFunctions[2]= (TF1*)formula1D.DrawClone("same");
  //
  // Draw Results
  //
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *canvasTest1D = new TCanvas("Test1D","Test1D",1000,800);
  canvasTest1D->Divide(1,2);
  canvasTest1D->cd(1);
  inputTree->Draw("test1D+noise:testx0>>his(100,0,1,100,-10,20)","abs(test1D+noise)<1","colz");
  for (Int_t ifit=0; ifit<3; ifit++){
    fitFunctions[ifit]->SetLineColor(kLineColors[ifit]);
    fitFunctions[ifit]->SetLineStyle(kLineStyle[ifit]);
    fitFunctions[ifit]->SetLineWidth(3);
    fitFunctions[ifit]->Draw("same");     
  }
  TLegend * legend = new TLegend(0.11,0.7,0.6,0.89,"Unbinned1D fit");
  legend->SetBorderSize(0);
  legend->AddEntry(fitFunctions[0],"Chi2 minimization","l");
  legend->AddEntry(fitFunctions[1],"Huber norm","l");
  legend->AddEntry(fitFunctions[2],"Custom log likelihood (80%Gaus+20%Cauchy)","l");
  legend->Draw();
  //
  canvasTest1D->cd(2);
  inputTree->Draw("test1D+noise-fitChi2Norm>>hisChi2Norm(100,-10,20)"); 
  resHistograms[0]=inputTree->GetHistogram();
  inputTree->Draw("test1D+noise-fitHuberNorm>>hisHuberNorm(100,-10,20)","","");
  resHistograms[1]=inputTree->GetHistogram();
  inputTree->Draw("test1D+noise-fitUserGausAndCachy>>hisUserGausAndCachy(100,-10,20)","","");
  resHistograms[2]=inputTree->GetHistogram();
  //for (Int_t ifit=0; ifit<3; ifit++) resHistograms[ifit]->Fit("gaus");
  for (Int_t ifit=0; ifit<3; ifit++){
    resHistograms[ifit]->SetLineColor(kLineColors[ifit]);
    resHistograms[ifit]->SetLineStyle(kLineStyle[ifit]);
    resHistograms[ifit]->SetLineWidth(3);
    resHistograms[ifit]->Draw("same");     
  }
  canvasTest1D->SaveAs("AliTMinuitToolkitTest.Test1D.png");


}
