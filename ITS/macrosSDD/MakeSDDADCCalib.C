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

//________________________________________________________________________________________

Double_t LangausFun(Double_t *x, Double_t *par) {
    
  //Fit parameters:
  //par[0] = Width (scale) parameter of Landau density
  //par[1] = Most Probable (MP, location) parameter of Landau density
  //par[2] = Total area (integral -inf to inf, normalization constant)
  //par[3] = Width (sigma) of convoluted Gaussian function
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

//________________________________________________________________________________________

void MakeSDDADCCalib(Int_t run = 245705,TString foldname = "LHC15o",TString filename ="QAresults_barrel.root"){
    
  TFile *fin = new TFile(Form("%s/%s",foldname.Data(),filename.Data()));
  TDirectory *dirFile = (TDirectory*)fin->Get("ITSAlignQA");
  TString lname = Form("clistITSAlignQA");
  TList *cOutput = (TList*)dirFile->Get(lname.Data());
    
  const Int_t nDrTimeBin = 8;
  const Int_t nModules = 260;
  Float_t drTimeLim[nDrTimeBin+1] = {0,800,1600,2400,3200,4000,4800,5600,6400};
    
  if(!cOutput) {
    Printf("E: Cannot open TList %s",lname.Data());
    return;
  }
    
  TH2F *hdEdxvsDrTime;
  TF1  *lfun = new TF1("LangausFun",LangausFun,50.,300.,4); //Langaus fit on a DrTime slice
    
  TH1F *hmpv  = new TH1F("hmpv","hmpv",nDrTimeBin,drTimeLim);
  TH1F *hsig  = new TH1F("hsig","hsig",nDrTimeBin,drTimeLim);
  TH1F *hsigl = new TH1F("hsigl","hsigl",nDrTimeBin,drTimeLim);
  TH1F *hmpvModpar0  = new TH1F("hmpvModpar0","hmpvModpar0; SDD Module; par0 MPV",260,239.5,499.5);
  TH1F *hmpvModpar1  = new TH1F("hmpvModpar1","hmpvModpar1; SDD Module; par1 MPV",260,239.5,499.5);
  TH1F *hsigModpar0  = new TH1F("hsigModpar0","hsigModpar0; SDD Module; par0 sigma Gaus",260,239.5,499.5);
  TH1F *hsigModpar1  = new TH1F("hsigModpar1","hsigModpar1; SDD Module; par1 sigma Gaus",260,239.5,499.5);
  TH1F *hsiglModpar0 = new TH1F("hsiglModpar0","hsiglModpar0; SDD Module; par0 sigma Landau",260,239.5,499.5);
  TH1F *hsiglModpar1 = new TH1F("hsiglModpar1","hsiglModpar1; SDD Module; par1 sigma Landau",260,239.5,499.5);
    
  TString fout = Form("SDDADCCalibResults_%d.root",run);
  TFile *out = new TFile(Form("%s/%s",foldname.Data(),fout.Data()),"recreate");
    
  //    TCanvas *chdEdxproj = new TCanvas("chdEdxproj","chdEdxproj",1000,800);
  //    TCanvas *cmod       = new TCanvas("cmod","module",600,900);
  //    chdEdxproj->Divide(4,2);
  //    cmod->Divide(1,3);
    
  for(Int_t ihist = 0; ihist < nModules; ihist++){//loop on modules
        
    Int_t imod = ihist+240;
    hdEdxvsDrTime = (TH2F*)cOutput->FindObject(Form("hSDDdEdxvsDrTime%i",imod));
    //        chdEdxproj->Clear("D");
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
        
    for(Int_t idEdx = 0; idEdx < nDrTimeBin; idEdx++){//loop on DrTime Bin
            
      int minTimeBin = hdEdxvsDrTime->GetXaxis()->FindBin(drTimeLim[idEdx]);
      int maxTimeBin = hdEdxvsDrTime->GetXaxis()->FindBin(drTimeLim[idEdx+1]);
      hdEdxproj[idEdx] = (TH1F*)hdEdxvsDrTime->ProjectionY(Form("%i",idEdx),minTimeBin,maxTimeBin);
      hdEdxproj[idEdx]->SetTitle(hdEdxvsDrTime->GetTitle());
            
      lfun = new TF1(Form("LangausFun%d",imod),LangausFun,50.,300.,4); //Langaus fit on a DrTime slice
      lfun->SetLineWidth(2);
      lfun->SetParameter(0,5.);
      lfun->SetParameter(1,80.);
      lfun->SetParameter(2,hdEdxproj[idEdx]->GetEntries()/10.);
      lfun->SetParameter(3,10.);
      lfun->SetParLimits(3,0.,20);
            
      hdEdxproj[idEdx]->Fit(lfun,"0NQLR");
      hdEdxproj[idEdx]->GetXaxis()->SetTitle(Form("dE/dx, time interval %d",idEdx));
      hdEdxproj[idEdx]->GetYaxis()->SetTitle("Events");
            
      Float_t mpv   = lfun->GetParameter(1);
      Float_t empv  = lfun->GetParError(1);
      Float_t sig   = lfun->GetParameter(3);
      Float_t esig  = lfun->GetParError(3);
      Float_t sigl  = lfun->GetParameter(0);
      Float_t esigl = lfun->GetParError(0);
            
      //filling histos mpv vs DrTime for each module
      hmpv->SetBinContent(idEdx+1,mpv);
      hmpv->SetBinError(hmpv->FindBin(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1])),empv);
      hsig->SetBinContent(idEdx+1,sig);
      hsig->SetBinError(hsig->FindBin(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1])),esig);
      hsigl->SetBinContent(idEdx+1,sigl);
      hsigl->SetBinError(hsigl->FindBin(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1])),esigl);
            
      //            chdEdxproj->cd(idEdx+1);
      //            hdEdxproj[idEdx]->Draw();
      //            lfun->Draw("same");
      //            chdEdxproj->Update();
    }//end loop on DrTime slices
        
    TF1 *pol1mpv = new TF1("pol1mpv","pol1mpv",0,6400);
    if(imod!=469) hmpv->Fit(pol1mpv,"0NQLR","",1000,5400);
    else hmpv->Fit(pol1mpv,"0NQLR","",1000,3000); //Mod 469 only one part of the dr Time region is full
    if(imod==376) {//Hard coded, bad module
      pol1mpv->SetParameter(0,84);
      pol1mpv->SetParameter(1,0);
    }
    hmpv->Write();
    hmpvModpar1->Fill(imod,pol1mpv->GetParameter(1));
    hmpvModpar1->SetBinError(hmpvModpar1->FindBin(imod),pol1mpv->GetParError(1));
    hmpvModpar0->Fill(imod,pol1mpv->GetParameter(0));
    hmpvModpar0->SetBinError(hmpvModpar0->FindBin(imod),pol1mpv->GetParError(0));
    //Printf("Par0_mpv:%f Err_Par0_mpv:%f    Par1_mpv:%f Err_Par1_mpv:%f",pol1mpv->GetParameter(0),pol1mpv->GetParError(0),pol1mpv->GetParameter(1),pol1mpv->GetParError(1));
        
    //        cmod->cd(1);
    //        hmpv->Draw();
    //        pol1mpv->Draw("SAME");
    //        cmod->Update();
    //        hmpv->Reset("M");
        
    TF1 *pol1sig = new TF1("pol1sig","pol1sig",0,6400);
    hsig->Fit(pol1sig,"0NQLR");
    hsigModpar0->Fill(imod,pol1sig->GetParameter(0));
    hsigModpar0->SetBinError(hsigModpar0->FindBin(imod),pol1sig->GetParError(0));
    hsigModpar1->Fill(imod,pol1sig->GetParameter(1));
    hsigModpar1->SetBinError(hsigModpar1->FindBin(imod),pol1sig->GetParError(1));
        
    //        cmod->cd(2);
    //        hsig->Draw();
    //        pol1sig->Draw("SAME");
    //        cmod->Update();
    //        hsig->Reset("M");
        
    TF1 *pol1sigl = new TF1("pol1sigl","pol1sigl",0,6400);
    hsigl->Fit(pol1sigl,"0NQLR");
    hsiglModpar0->Fill(imod,pol1sigl->GetParameter(0));
    hsiglModpar0->SetBinError(hsiglModpar0->FindBin(imod),pol1sigl->GetParError(0));
    hsiglModpar1->Fill(imod,pol1sigl->GetParameter(1));
    hsiglModpar1->SetBinError(hsiglModpar1->FindBin(imod),pol1sigl->GetParError(1));
        
    //        cmod->cd(3);
    //        hsigl->Draw();
    //        pol1sigl->Draw("SAME");
    //        cmod->Update();
    //        hsigl->Reset("M");
        
  }//end loop on modules
    
  Double_t offset=1.3;
  hmpvModpar0->GetYaxis()->SetTitleOffset(offset);
  hmpvModpar1->GetYaxis()->SetTitleOffset(offset);
  hsigModpar0->GetYaxis()->SetTitleOffset(offset);
  hsigModpar1->GetYaxis()->SetTitleOffset(offset);
  hsiglModpar0->GetYaxis()->SetTitleOffset(offset);
  hsiglModpar1->GetYaxis()->SetTitleOffset(offset);
    
  TCanvas *chmpvModpar0  = new TCanvas("chmpvModpar0","chmpvModpar0",1000,800);
  hmpvModpar0->Draw();
  TCanvas *chmpvModpar1  = new TCanvas("chmpvModpar1","chmpvModpar1",1000,800);
  hmpvModpar1->Draw();
  TCanvas *chsigModpar0  = new TCanvas("chsigModpar0","chsigModpar0",1000,800);
  hsigModpar0->Draw();
  TCanvas *chsigModpar1  = new TCanvas("chsigModpar1","chsigModpar1",1000,800);
  hsigModpar1->Draw();
  TCanvas *chsiglModpar0 = new TCanvas("chsiglModpar0","chsiglModpar0",1000,800);
  hsiglModpar0->SetMinimum(0);
  hsiglModpar0->SetMaximum(20);
  hsiglModpar0->Draw();
  TCanvas *chsiglModpar1 = new TCanvas("chsiglModpar1","chsiglModpar1",1000,800);
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


