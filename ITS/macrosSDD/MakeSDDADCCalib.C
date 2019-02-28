#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLine.h>
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


void MakeSDDADCCalib(Int_t run = 245705,TString foldname = "15o_Bunch4",TString filename ="CalibObjects.root", Int_t year = 2015, TString period = "LHC15o", Bool_t readLocal=kFALSE){
    
  //****************** Connection to alien *****************************************
    
  if(!readLocal){
    TGrid::Connect("alien://",0,0,"t");
    //TGrid *gGrid = TGrid::Connect("alien");
    if(!gGrid||!gGrid->IsConnected()) {
      printf("gGrid not found! exit macro\n");
      return;
    }
  }
    
  const Int_t nDrTimeBin = 8;
  const Int_t nModules = 260;
  Float_t drTimeLim[nDrTimeBin+1] = {0,800,1600,2400,3200,4000,4800,5600,6400};
    
    
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
    
  TCanvas *chdEdxproj = new TCanvas("chdEdxproj","chdEdxproj",1600,800);
  //    TCanvas *cmod       = new TCanvas("cmod","module",600,900);
  chdEdxproj->Divide(5,2);
  //    cmod->Divide(1,3);
  Bool_t firstpage=kTRUE;

  Char_t path[200];
    
  //    sprintf(path,"alien:///alice/data/%04i/%s/%09i/cpass1_pass1/QAresults_barrel.root",year,period.Data(),run);
  //    TFile *fin = TFile::Open(path);
  //    if(!fin)return;
  //    TDirectory *dirFile = (TDirectory*)fin->Get("ITSAlignQA");
  //    TString lname = Form("clistITSAlignQA");
  //    TList *cOutput = (TList*)dirFile->Get(lname.Data());
  //    if(!cOutput) {
  //        Printf("E: Cannot open TList %s",lname.Data());
  //        return;
  //    }
    
  TFile *fin;
  if(!readLocal) {
    sprintf(path,"alien:///alice/data/%04i/%s/%09i/cpass1_pass1/OCDB/CalibObjects.root",year,period.Data(),run);
    //        sprintf(path,"alien:///alice/data/%04i/%s/%09i/zdc_special_wTPC_cpass1/OCDB/CalibObjects.root",year,period.Data(),run);
    fin = TFile::Open(path);
  }
  else fin = new TFile(Form("%s/%s",foldname.Data(),filename.Data()));
    
  if(!fin)return;
  TString lname = Form("clistSDDCalib");
  TList *cOutput= 0x0;
  if(filename.Contains("CalibObjects")) cOutput = (TList*)fin->Get(lname.Data());
  else{
    TDirectoryFile* drf=(TDirectoryFile*)fin->Get("ITSAlignQA");
    lname="clistITSAlignQA";
    if(drf) cOutput = (TList*)drf->Get(lname.Data());
  }
  
  if(!cOutput) {
    Printf("E: Cannot open TList %s",lname.Data());
    return;
  }
  TLatex* textTooFew=new TLatex(0.2,0.8,"Too few entries");
  textTooFew->SetNDC();
  textTooFew->SetTextFont(63);
  textTooFew->SetTextColor(2);
  textTooFew->SetTextSize(28);
  
  TLatex* textfldFit=new TLatex(0.2,0.8,"Fit did not converge");
  textfldFit->SetNDC();
  textfldFit->SetTextFont(63);
  textfldFit->SetTextColor(2);
  textfldFit->SetTextSize(28);
    
  TLatex* textbadFit=new TLatex(0.2,0.8,"Bad fit params");
  textbadFit->SetNDC();
  textbadFit->SetTextFont(63);
  textbadFit->SetTextColor(2);
  textbadFit->SetTextSize(28);

  TLine* refLine=new TLine(0.,84.,6400.,84.);
  refLine->SetLineColor(kGreen+1);
  refLine->SetLineWidth(2);
  refLine->SetLineStyle(2);
  TLine* refLineL=new TLine(0.,79.,6400.,79.);
  refLineL->SetLineColor(kRed-7);
  refLineL->SetLineWidth(2);
  refLineL->SetLineStyle(3);
  TLine* refLineH=new TLine(0.,89.,6400.,89.);
  refLineH->SetLineColor(kRed-7);
  refLineH->SetLineWidth(2);
  refLineH->SetLineStyle(3);

  //    for(Int_t ihist = 6; ihist < 7; ihist++){//loop on modules
  for(Int_t ihist = 0; ihist < nModules; ihist++){//loop on modules
    Int_t imod = ihist+240;
    hdEdxvsDrTime = (TH2F*)cOutput->FindObject(Form("hSDDdEdxvsDrTime%i",imod));
    chdEdxproj->Clear("D");
    printf("Mod. # %i \n",imod);
    hmpv->Reset("ICESM");
    hsig->Reset("ICESM");
    hsigl->Reset("ICESM");
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
    Double_t minMPV=999999;
    Double_t maxMPV=1;
    for(Int_t idEdx = 0; idEdx < nDrTimeBin; idEdx++){//loop on DrTime Bin
            
      int minTimeBin = hdEdxvsDrTime->GetXaxis()->FindBin(drTimeLim[idEdx]);
      int maxTimeBin = hdEdxvsDrTime->GetXaxis()->FindBin(drTimeLim[idEdx+1]);
      hdEdxproj[idEdx] = (TH1F*)hdEdxvsDrTime->ProjectionY(Form("%i",idEdx),minTimeBin,maxTimeBin);
      hdEdxproj[idEdx]->SetTitle(hdEdxvsDrTime->GetTitle());
      hdEdxproj[idEdx]->GetXaxis()->SetTitle(Form("dE/dx, time interval %d",idEdx));
      hdEdxproj[idEdx]->GetYaxis()->SetTitle("Events");

      chdEdxproj->cd(idEdx+1);
      hdEdxproj[idEdx]->Draw();

      if(hdEdxproj[idEdx]->GetEntries()<50){
	printf("Module %d, time bin %d, too few entries, SKIP\n",imod,idEdx);
	textTooFew->Draw();
	chdEdxproj->Update();
	continue;
      }
      lfun = new TF1(Form("LangausFun%d",imod),LangausFun,50.,300.,4); //Langaus fit on a DrTime slice
      lfun->SetLineWidth(2);
      lfun->SetParameter(0,5.);
      lfun->SetParameter(1,80.);
      lfun->SetParameter(2,hdEdxproj[idEdx]->GetEntries()/10.);
      lfun->SetParameter(3,10.);
      lfun->SetParLimits(3,1,20);

      Int_t fitstat = hdEdxproj[idEdx]->Fit(lfun,"0NQLR");
      lfun->Draw("same");
      chdEdxproj->Update();

      if(fitstat!=0){
	printf("   Fit failed in time interval %d\n",idEdx);
	textfldFit->Draw();
	chdEdxproj->Update();
	continue;
      }
      Float_t mpv   = lfun->GetParameter(1);
      Float_t empv  = lfun->GetParError(1);
      Float_t sig   = lfun->GetParameter(3);
      Float_t esig  = lfun->GetParError(3);
      Float_t sigl  = lfun->GetParameter(0);
      Float_t esigl = lfun->GetParError(0);
      if(sigl<1. || mpv<40. || mpv>200.){
	printf("   Bad fit parameters in time interval %d: mpv=%f sig=%f sigl=%f\n",idEdx,mpv,sig,sigl);
	textbadFit->Draw();
	chdEdxproj->Update();
	continue;
      }
      printf("   TimeBin %d Entries =%.0f  mpv=%f  sig=%f sigl=%f\n",idEdx,hdEdxproj[idEdx]->GetEntries(),mpv,sig,sigl);
      //filling histos mpv vs DrTime for each module
      if((mpv+empv)>maxMPV) maxMPV=mpv+empv;
      if((mpv-empv)<minMPV) minMPV=mpv-empv;
      hmpv->SetBinContent(idEdx+1,mpv);
      hmpv->SetBinError(hmpv->FindBin(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1])),empv);
      hsig->SetBinContent(idEdx+1,sig);
      hsig->SetBinError(hsig->FindBin(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1])),esig);
      hsigl->SetBinContent(idEdx+1,sigl);
      hsigl->SetBinError(hsigl->FindBin(0.5*(drTimeLim[idEdx]+drTimeLim[idEdx+1])),esigl);

    }//end loop on DrTime slices
    Double_t maxFit=5400;
    Double_t minFit=1000;
    if(minMPV<maxMPV){
      hmpv->SetMinimum(0.95*minMPV);
      hmpv->SetMaximum(1.02*maxMPV);
    }
    TF1 *pol1mpv = new TF1("pol1mpv","pol1",0,6400);
    //Mod 469 only one part of the dr Time region is full 
    if(imod==469) maxFit=3000;
    if(hmpv->GetEntries()<=3) pol1mpv->FixParameter(1,0);
    if(hmpv->GetEntries()>2){
      hmpv->Fit(pol1mpv,"0NQ","",minFit,maxFit);
    }else{
      pol1mpv->SetParameter(0,84);
      pol1mpv->SetParameter(1,0);
    }
    if(imod==376) {//Hard coded, bad module
      pol1mpv->SetParameter(0,84);
      pol1mpv->SetParameter(1,0);
    }
    chdEdxproj->cd(nDrTimeBin+1);
    if(hmpv->GetMaximum()<85) hmpv->SetMaximum(85);
    if(hmpv->GetMinimum()>83) hmpv->SetMinimum(83);
    hmpv->Draw();
    pol1mpv->Draw("same");
    refLine->Draw("same");
    if(hmpv->GetMinimum()<79) refLineL->Draw("same");
    if(hmpv->GetMaximum()>89)refLineH->Draw("same");
    chdEdxproj->Update();
    if(firstpage){
      chdEdxproj->Print("LanGausFits.pdf(");
      firstpage=kFALSE;
    }else{
      chdEdxproj->Print("LanGausFits.pdf");
    }
    out->cd();
    hmpv->Write();
    hmpvModpar1->Fill(imod,pol1mpv->GetParameter(1));
    hmpvModpar1->SetBinError(hmpvModpar1->FindBin(imod),pol1mpv->GetParError(1));
    hmpvModpar0->Fill(imod,pol1mpv->GetParameter(0));
    hmpvModpar0->SetBinError(hmpvModpar0->FindBin(imod),pol1mpv->GetParError(0));
    Printf("Par0_mpv:%f Err_Par0_mpv:%f    Par1_mpv:%f Err_Par1_mpv:%f",pol1mpv->GetParameter(0),pol1mpv->GetParError(0),pol1mpv->GetParameter(1),pol1mpv->GetParError(1));
        
    //        cmod->cd(1);
    //        hmpv->Draw();
    //        pol1mpv->Draw("SAME");
    //        cmod->Update();
    //        hmpv->Reset("M");
        
    TF1 *pol1sig = new TF1("pol1sig","pol1",0,6400);
    if(hsig->GetEntries()>2){
      hsig->Fit(pol1sig,"0NQ","",minFit,maxFit);
    }
    hsigModpar0->Fill(imod,pol1sig->GetParameter(0));
    hsigModpar0->SetBinError(hsigModpar0->FindBin(imod),pol1sig->GetParError(0));
    hsigModpar1->Fill(imod,pol1sig->GetParameter(1));
    hsigModpar1->SetBinError(hsigModpar1->FindBin(imod),pol1sig->GetParError(1));
        
    //        cmod->cd(2);
    //        hsig->Draw();
    //        pol1sig->Draw("SAME");
    //        cmod->Update();
    //        hsig->Reset("M");
        
    TF1 *pol1sigl = new TF1("pol1sigl","pol1",0,6400);
    if(hsigl->GetEntries()>2){
      hsigl->Fit(pol1sigl,"0NQ","",minFit,maxFit);
    }
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
    
  TCanvas *chmpv  = new TCanvas("chmpv","chmpv",1600,800);
  chmpv->Divide(2,1);
  chmpv->cd(1);
  hmpvModpar0->Draw();
  chmpv->cd(2);
  hmpvModpar1->Draw();
  chmpv->Print("LanGausFits.pdf");


  TCanvas *chsig  = new TCanvas("chsig","chsig",1600,800);
  chsig->Divide(2,1);
  chsig->cd(1);
  hsigModpar0->Draw();
  chsig->cd(2);
  hsigModpar1->Draw();
  chsig->Print("LanGausFits.pdf");

  TCanvas *chsigl = new TCanvas("chsigl","chsigl",1600,800);
  chsigl->Divide(2,1);
  chsigl->cd(1);
  hsiglModpar0->SetMinimum(0);
  hsiglModpar0->SetMaximum(20);
  hsiglModpar0->Draw();
  chsigl->cd(2);
  hsiglModpar1->Draw();
  chsigl->Print("LanGausFits.pdf)");
   
  out->cd();
  hmpvModpar0->Write();
  hmpvModpar1->Write();
  hsigModpar0->Write();
  hsigModpar1->Write();
  hsiglModpar0->Write();
  hsiglModpar1->Write();
  out->Close();
  delete out;
    
} //end main


