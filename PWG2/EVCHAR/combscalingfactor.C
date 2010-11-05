#include "Riostream.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h" 
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCanvas.h"
void combscalingfactor(Char_t *fileall, Char_t *filecomb, Bool_t scalebkg = kFALSE, Float_t scfac=1.) {

 //  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);
//  gStyle->SetOptLogy(kFALSE);
  gStyle->SetFrameLineWidth(2.5);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetHistLineWidth(2.5);
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleSize(0.05, "x");
  gStyle->SetTitleSize(0.05, "y");
  gStyle->SetTitleOffset(1.1, "x");
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFillColor(0);
//  gROOT->ForceStyle();



 TFile *fall = new TFile(fileall);
 TFile *fcomb = new TFile(filecomb);

 TList *lall = (TList*) fall->Get("cOutput");
 TList *lcomb = (TList*) fcomb->Get("cOutput");

 TH2F *hall = (TH2F*) lall->FindObject("fHistSPDdePhideTheta");
 TH2F *hcomb = (TH2F*) lcomb->FindObject("fHistSPDdePhideTheta");

 TH1D* hproall = new TH1D("_xall","",2000,-1.,1.);
 hproall->Sumw2(); 
 hproall = hall->ProjectionX("_xall",0,-1,"e");

 TH1D* hprocomb = new TH1D("_xcomb","",2000,-1.,1.);
 hprocomb->Sumw2(); 
 hprocomb = hcomb->ProjectionX("_xcomb",0,-1,"e");


 TH1D* hratio = new TH1D("ratio","",2000,-1.,1.); 
 hratio->Divide(hproall,hprocomb,1.,1.); 

 TF1 *fleft = new TF1("fleft","pol0",-.8,-0.6); //check values
 TF1 *fright = new TF1("fright","pol0",0.6,.8);

 hratio->Fit("fleft","R");
 hratio->Fit("fright","R");

/*
 // Projection on DeltaTheta 
 TH1D* hproall = new TH1D("_xall","",1000,-.25,.25);
 hproall->Sumw2(); 
 hproall = hall->ProjectionY("_xall",0,-1,"e");

 TH1D* hprocomb = new TH1D("_xcomb","",1000,-.25,.25);
 hprocomb->Sumw2(); 
 hprocomb = hcomb->ProjectionY("_xcomb",0,-1,"e");

 TH1D* hratio = new TH1D("ratio","",1000,-.25,.25); 
 hratio->Divide(hproall,hprocomb,1.,1.);
 
 TF1 *fleft = new TF1("fleft","pol0",-.8,-0.6); //check values
 TF1 *fright = new TF1("fright","pol0",0.6,.8);
 hratio->Fit("fleft","R");
 hratio->Fit("fright","R");
*/

  TFile *fout= new TFile("scaling.root","RECREATE");

  new TCanvas();
  hproall->Draw("histo");
  
  if (scalebkg) hprocomb->Scale(scfac);
  hprocomb->Draw("histo,same");

  cout<<"Percentage of background"<< hprocomb->Integral(1,1000)/hproall->Integral(1,1000)<<endl;
  cout<<"Percentage of background in |#D#phi|<0.08"<< hprocomb->Integral(921,1080)/hproall->Integral(921,1080)<<endl; 
  hproall->Write();
  hprocomb->SetLineColor(kBlue);
  hprocomb->Write(); 
  fout->Write();
  fout->Close();

}
