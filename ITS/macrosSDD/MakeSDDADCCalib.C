#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#endif

// For each module, finds the dEdx distribution in 8 drift time intervals. Fits each distribution with LanGaus. For each module builds the histogram of
// Most Probable Value vs drift time and fits it with pol1. Stores the fit parameters vs module number in histograms hmpvModpar0 and hmpvModpar1.

Double_t LangausFun(Double_t *x, Double_t *par);


void MakeSDDADCCalib(TString foldname = ".",TString filename ="QAresults_barrel.root")
{
  
  TFile *fin = new TFile(Form("%s/%s",foldname.Data(),filename.Data()));
  TString lname=Form("clistITSAlignQA");
  TDirectory *dirFile=(TDirectory*)fin->Get("ITSAlignQA");
  TList *cOutput=0x0;
  cOutput = (TList*)dirFile->Get(lname.Data());
  
  const Int_t nDrTimeBin=8;
  Float_t drTimeLim[nDrTimeBin+1]={0,800,1600,2400,3200,4000,4800,5600,6400};
  
  if(!cOutput){
    Printf("E: Cannot open TList %s",lname.Data());
    return;
  }
  
  TH2F *hdEdxvsDrTime[260];
  TF1 *lfun[260]; 
  TH1F *hmpv = new TH1F("hmpv","hmpv",nDrTimeBin,drTimeLim);
  TH1F *hsig = new TH1F("hsig","hsig",nDrTimeBin,drTimeLim);
  TH1F *hsigl = new TH1F("hsigl","hsigl",nDrTimeBin,drTimeLim);
  TH1F *hmpvModpar0 = new TH1F("hmpvModpar0","hmpvModpar0",280,229.5,509.5);
  TH1F *hmpvModpar1 = new TH1F("hmpvModpar1","hmpvModpar1",280,229.5,509.5);
  TH1F *hsigModpar0 = new TH1F("hsigModpar0","hsigModpar0",280,229.5,509.5);
  TH1F *hsigModpar1 = new TH1F("hsigModpar1","hsigModpar1",280,229.5,509.5);
  TH1F *hsiglModpar0 = new TH1F("hsiglModpar0","hsiglModpar0",280,229.5,509.5);
  TH1F *hsiglModpar1 = new TH1F("hsiglModpar1","hsiglModpar1",280,229.5,509.5);

  TString fout="SDDADCCalibResults.root";
  TFile *out=new TFile(Form("%s/%s",foldname.Data(),fout.Data()),"recreate");

  TCanvas *chdEdxproj=new TCanvas("chdEdxproj","chdEdxproj",1000,800);
  chdEdxproj->Divide(4,2);
  TCanvas* cmod=new TCanvas("cmod","module",600,900);
  cmod->Divide(1,3);
  for(Int_t ihist=0;ihist<260;ihist++){//loop on modules
  // for(Int_t ihist=0;ihist<10;ihist++){//only 10 to test
    Int_t imod=ihist+240;
    hdEdxvsDrTime[ihist]=(TH2F*)cOutput->FindObject(Form("hSDDdEdxvsDrTime%i",imod));
    chdEdxproj->Clear("D");
    //    chdEdxproj->Divide(4,2);
    printf("Mod. # %i \n",imod);
    hmpv->SetTitle(Form("MPV Mod. # %i",imod));
    hmpv->SetName(Form("MPVModule%i",imod));
    hsig->SetTitle(Form("Gauss sigma Mod. # %i",imod));
    hsigl->SetTitle(Form("Landau sigma Mod. # %i",imod));
    hmpv->GetXaxis()->SetTitle("Drift time (ns)");
    hsig->GetXaxis()->SetTitle("Drift time (ns)");
    hsigl->GetXaxis()->SetTitle("Drift time (ns)");
    hmpv->GetYaxis()->SetTitle("MPV (keV)");
    hsig->GetYaxis()->SetTitle("Gaussian #sigma (keV)");
    hsigl->GetYaxis()->SetTitle("Landau #sigma (keV)");
    TH1F *hdEdxproj[nDrTimeBin];
    for(Int_t idEdx=0;idEdx<nDrTimeBin;idEdx++){//loop on DrTime Bin
      hdEdxproj[idEdx]=(TH1F*)hdEdxvsDrTime[ihist]->ProjectionY(Form("%i",idEdx),hdEdxvsDrTime[ihist]->GetXaxis()->FindBin(drTimeLim[idEdx]),hdEdxvsDrTime[ihist]->GetXaxis()->FindBin(drTimeLim[idEdx+1])); 
      hdEdxproj[idEdx]->SetTitle(hdEdxvsDrTime[ihist]->GetTitle());
      lfun[ihist] = new TF1(Form("LangausFun%d",imod),LangausFun,50.,300.,4); //Langaus fit on a DrTime slice
      lfun[ihist]->SetLineWidth(2);
      lfun[ihist]->SetParameter(0,5.);
      lfun[ihist]->SetParameter(1,80.);
      lfun[ihist]->SetParameter(2,hdEdxproj[idEdx]->GetEntries()/10.);
      lfun[ihist]->SetParameter(3,10.);
      lfun[ihist]->SetParLimits(3,0.,20);
      hdEdxproj[idEdx]->Fit(lfun[ihist],"NQLR");
      hdEdxproj[idEdx]->GetXaxis()->SetTitle(Form("dE/dx, time interval %d",idEdx));
      hdEdxproj[idEdx]->GetYaxis()->SetTitle("Events");
      Float_t mpv=lfun[ihist]->GetParameter(1);
      Float_t empv=lfun[ihist]->GetParError(1);
      Float_t sig=lfun[ihist]->GetParameter(3);
      Float_t esig=lfun[ihist]->GetParError(3);
      Float_t sigl=lfun[ihist]->GetParameter(0);
      Float_t esigl=lfun[ihist]->GetParError(0);
      //filling histos mpv vs DrTime for each module
      hmpv->Fill(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1]),mpv);
      hmpv->SetBinError(hmpv->FindBin(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1])),empv);
      hsig->Fill(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1]),sig);
      hsig->SetBinError(hsig->FindBin(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1])),esig);
      hsigl->Fill(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1]),sigl);
      hsigl->SetBinError(hsigl->FindBin(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1])),esigl);
      chdEdxproj->cd(idEdx+1);
      hdEdxproj[idEdx]->Draw();
      lfun[ihist]->Draw("same");
      chdEdxproj->Update();
    }//end loop on DrTime slices
    TF1 *pol1mpv = new TF1("pol1mpv","pol1mpv",0,6400);
    if(imod!=469)hmpv->Fit(pol1mpv,"NQLR","",1000,5400);
    else hmpv->Fit(pol1mpv,"NQLR","",1000,3000); //Mod 469 only one part of the dr Time region is full
    if(imod==376){//Hard coded, bad module
      pol1mpv->SetParameter(0,84);
      pol1mpv->SetParameter(1,0);
    }
    hmpv->Write();
    hmpvModpar1->Fill(imod,pol1mpv->GetParameter(1));
    hmpvModpar1->SetBinError(hmpvModpar1->FindBin(imod),pol1mpv->GetParError(1));
    Printf("Par0_mpv:%f Err_Par0_mpv:%f    Par1_mpv:%f Err_Par1_mpv:%f",pol1mpv->GetParameter(0),pol1mpv->GetParError(0),pol1mpv->GetParameter(1),pol1mpv->GetParError(1));
    hmpvModpar0->Fill(imod,pol1mpv->GetParameter(0));
    hmpvModpar0->SetBinError(hmpvModpar0->FindBin(imod),pol1mpv->GetParError(0));
    cmod->cd(1);
    hmpv->Draw();
    pol1mpv->Draw("SAME");    
    cmod->Update();
    hmpv->Reset("M");
    
    TF1 *pol1sig = new TF1("pol1sig","pol1sig",0,6400);
    hsig->Fit(pol1sig,"NQLR");
    hsigModpar0->Fill(imod,pol1sig->GetParameter(0));
    hsigModpar0->SetBinError(hsigModpar0->FindBin(imod),pol1sig->GetParError(0));
    hsigModpar1->Fill(imod,pol1sig->GetParameter(1));
    hsigModpar1->SetBinError(hsigModpar1->FindBin(imod),pol1sig->GetParError(1));
    cmod->cd(2);
    hsig->Draw();
    pol1sig->Draw("SAME");
    cmod->Update();
    hsig->Reset("M");
    
    TF1 *pol1sigl = new TF1("pol1sigl","pol1sigl",0,6400);
    hsigl->Fit(pol1sigl,"NQLR");
    hsiglModpar0->Fill(imod,pol1sigl->GetParameter(0));
    hsiglModpar0->SetBinError(hsiglModpar0->FindBin(imod),pol1sigl->GetParError(0));
    hsiglModpar1->Fill(imod,pol1sigl->GetParameter(1));
    hsiglModpar1->SetBinError(hsiglModpar1->FindBin(imod),pol1sigl->GetParError(1));
    cmod->cd(3);
    hsigl->Draw();
    pol1sigl->Draw("SAME");
    cmod->Update();
    hsigl->Reset("M");

  }//end loop on modules
    
  hmpvModpar0->GetXaxis()->SetTitle("SDD Module");
  hmpvModpar1->GetXaxis()->SetTitle("SDD Module");
  hsigModpar0->GetXaxis()->SetTitle("SDD Module");
  hsigModpar1->GetXaxis()->SetTitle("SDD Module");
  hsiglModpar0->GetXaxis()->SetTitle("SDD Module");
  hsiglModpar1->GetXaxis()->SetTitle("SDD Module");
  hmpvModpar0->GetYaxis()->SetTitle("par0 MPV");
  hmpvModpar1->GetYaxis()->SetTitle("par1 MPV");
  hsigModpar0->GetYaxis()->SetTitle("par0 sigma Gaus");
  hsigModpar1->GetYaxis()->SetTitle("par1 sigma Gaus");
  hsiglModpar0->GetYaxis()->SetTitle("par0 sigma Landau");
  hsiglModpar1->GetYaxis()->SetTitle("par1 sigma Landau");
  Double_t offset=1.3;
  hmpvModpar0->GetYaxis()->SetTitleOffset(offset);
  hmpvModpar1->GetYaxis()->SetTitleOffset(offset);
  hsigModpar0->GetYaxis()->SetTitleOffset(offset);
  hsigModpar1->GetYaxis()->SetTitleOffset(offset);
  hsiglModpar0->GetYaxis()->SetTitleOffset(offset);
  hsiglModpar1->GetYaxis()->SetTitleOffset(offset);
    
  TCanvas *chmpvModpar0=new TCanvas("chmpvModpar0","chmpvModpar0",1000,800);
  chmpvModpar0->cd();
  hmpvModpar0->Draw();
  TCanvas *chmpvModpar1=new TCanvas("chmpvModpar1","chmpvModpar1",1000,800);
  chmpvModpar1->cd();
  hmpvModpar1->Draw();
  TCanvas *chsigModpar0=new TCanvas("chsigModpar0","chsigModpar0",1000,800);
  chsigModpar0->cd();
  hsigModpar0->Draw();
  TCanvas *chsigModpar1=new TCanvas("chsigModpar1","chsigModpar1",1000,800);
  chsigModpar1->cd();
  hsigModpar1->Draw();
  TCanvas *chsiglModpar0=new TCanvas("chsiglModpar0","chsiglModpar0",1000,800);
  chsiglModpar0->cd();
  hsiglModpar0->SetMinimum(0);
  hsiglModpar0->SetMaximum(20);
  hsiglModpar0->Draw();
  TCanvas *chsiglModpar1=new TCanvas("chsiglModpar1","chsiglModpar1",1000,800);
  chsiglModpar1->cd();
  hsiglModpar1->Draw();
  
  hmpvModpar0->Write();
  hmpvModpar1->Write();
  hsigModpar0->Write();
  hsigModpar1->Write();
  hsiglModpar0->Write();
  hsiglModpar1->Write();
  out->Close();
  delete out;
  
} //end main


Double_t LangausFun(Double_t *x, Double_t *par) {
  
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}
