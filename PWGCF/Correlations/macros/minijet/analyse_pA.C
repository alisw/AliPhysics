/* $Id:  $ */
//--------------------------------------------------
//
// macro to do the final analysis step 
// uses input of analysis class AliAnalysisTaskPhiCorrelation
//
//  Author : Emilia Leogrande (University of Utrecht)
//
//-------------------------------------------------

#include <TChain.h>
#include <TList.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TCanvas.h>
#include "TRandom.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TF1.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TStyle.h" 
#include "TROOT.h"
#include "TColor.h"

#include <iostream>
using namespace std;

void analyseEmy2(Bool_t zyam = kTRUE);  // if zyam = kFALSE, fit is used
Double_t fitFunction(Double_t *x ,Double_t *par); // fit function using constant + 3 gaussians
Double_t fitFunction2Gaus(Double_t *x ,Double_t *par); // fit function using constant + 2 gaussians

//input file and mixed event removed file
TFile *fileData=0x0;
TFile *fileDataEMremoved = 0x0;

const int multclass = 20;

TH1D *fDeltaPhiNch[multclass];
TH1D *fDeltaEtaNch[multclass];
TH1D *fSignalDPhi[multclass];
TH1D *fSignalNSDPhi[multclass];
TH1D *fSignalASDPhi[multclass];
TH1D *fRidge1DPhi[multclass];
TH1D *fRidge2DPhi[multclass];
TH1D *fRidgeDPhi[multclass];
TH1D *fSymmRidgeNotScaled[multclass];
TH1D *fSymmRidge[multclass];
TH1D *fFinal1DPhi[multclass];
TH1D *fFinalDPhi[multclass];

TString flag = "R";
TF1 *fTotal2Gaus[multclass];       // fit with 2 gaussians + const
TF1 *fTotal[multclass];            // fit with 3 gaussians + const

//properties of histogram
const int bins = 72; //
Double_t binWidth=2*TMath::Pi()/bins;

const int binsDeta = 48;


Double_t max_bin_for_etagap =   1.2;
Double_t min_bin_for_etagap =   -1.2;
Double_t max_eta = 1.8;
Double_t min_eta = -1.8;

//________________________________________________________________________________________________________________
//
Double_t fitFunction(Double_t *x ,Double_t *par)
{
  // fit function for 3 gaus + constant  
  
  // parameters for Gaussian
  Double_t A1     = par[0];
  Double_t sigma1 = par[1];
  Double_t A2     = par[2];
  Double_t sigma2 = par[3];
  Double_t A3     = par[4];
  Double_t sigma3 = par[5];
  Double_t integral = par[6];

  Double_t constante = (integral-
			TMath::Sqrt(TMath::Pi()*2)/ binWidth*
			(A1 * sigma1 + A2 * sigma2 + A3*sigma3))/bins;
  Double_t q  = x[0];
  
  //fit value
  Double_t fitval = constante +
    (q>-0.5*TMath::Pi()&&q<0.5*TMath::Pi())*(
					     A1 * exp(- q * q / (2 * sigma1 *sigma1)) +
					     A1 * exp(-((q - TMath::TwoPi())) * ((q - TMath::TwoPi())) / ( 2 * sigma1 * sigma1))
					     )
    +
    (q>-0.2*TMath::Pi()&&q<0.2*TMath::Pi())*(
					     A2 * exp(- q * q / (2 * sigma2 *sigma2)) +
					     A2 * exp(-((q - TMath::TwoPi())) * ((q - TMath::TwoPi())) / ( 2 * sigma2 * sigma2))
					     )
    +
    (q>0.5*TMath::Pi()&&q<1.5*TMath::Pi())*(
					    A3 * exp(-((q - TMath::Pi())) * ((q - TMath::Pi())) / ( 2 * sigma3 * sigma3)) +
					    A3 * exp(-((q + TMath::Pi())) * ((q + TMath::Pi())) / (2 * sigma3 * sigma3))
					    );
  return fitval;
}

//________________________________________________________________________________________________________________
//
Double_t fitFunction2Gaus(Double_t *x ,Double_t *par)
{
  // fit function for 2 gaus + constant  

  // parameters for Gaussian
  Double_t A1     = par[0];
  Double_t sigma1 = par[1];
  Double_t A3     = par[2];
  Double_t sigma3 = par[3];
  Double_t integral = par[4];

  Double_t constante = (integral -
			TMath::Sqrt(TMath::Pi()*2)/ binWidth*
			(A1 * sigma1 + A3*sigma3))/bins;
  Double_t q  = x[0];
  
  //fit value
  Double_t fitval = constante +
    (q>-0.5*TMath::Pi()&&q<0.5*TMath::Pi())*(
					     A1 * exp(- q * q / (2 * sigma1 *sigma1)) +
					     A1 * exp(-((q - TMath::TwoPi())) * ((q - TMath::TwoPi())) / ( 2 * sigma1 * sigma1)) 
					     )
    +
    (q>0.5*TMath::Pi()&&q<1.5*TMath::Pi())*(
					    A3 * exp(-((q - TMath::Pi())) * ((q - TMath::Pi())) / ( 2 * sigma3 * sigma3)) +
					    A3 * exp(-((q + TMath::Pi())) * ((q + TMath::Pi())) / (2 * sigma3 * sigma3))
					    );
  return fitval;
}

//_______________________________________________________________________________________________________________
//
Double_t fline(Double_t *x, Double_t *par){
    
    if(x[0]>-1.8 && x[0]<=0){
        return par[0]+par[1]*x[0];
    }
    else if(x[0]>0 && x[0]<1.8){
        return par[2]+par[3]*x[0];
    }
    else
        return 0;
}


//________________________________________________________________________________________________________________
//
void analyseEmy2(Bool_t zyam){


  // plot style
  gStyle->SetOptStat(0);
  const Int_t NRGBs = 5;
  const Int_t NCont = 500;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
    
  //style
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
    
  //-------------- TRIGGERS AND EVENTS
  
  TH2D *dphideta[multclass];
  TH1D * trigger = 0x0;
  TH1D * event = 0x0;
  
  fileData = TFile::Open("dphi_corr.root");
  trigger = (TH1D*)fileData->Get("triggers_0");
  event = (TH1D*)fileData->Get("events");
    
  // get average trigger particles per event
  TProfile *p0 = (TProfile*)trigger->Clone();
  TProfile *p1 = (TProfile*)event->Clone();
  p0->Sumw2();
  p1->Sumw2();
  p0->Divide(p0,p1,1,1,"B");
    
  // copy triggers and events in the new dphi_corr with the Mixed Event removed
  TH1D *triggerCopy = 0x0;
  TH1D *eventCopy = 0x0;
    
  triggerCopy = (TH1D*)trigger->Clone();
  eventCopy = (TH1D*)event->Clone();
    
  fileDataEMremoved = TFile::Open("dphi_corr_MEremoved.root","RECREATE");
  triggerCopy->SetName("triggers_0");
  triggerCopy->Write();
  eventCopy->SetName("events");
  eventCopy->Write();
  fileDataEMremoved->Close();
  
    
  //-------------- MIXED EVENT REMOVAL: restores the right number of particles in the detector acceptance but keeps the detector azimuthal unefficiencies corrections and cures the dip in (0,0) from two-trak cuts
  // Removing the event mixing: S/M (from dphi_corr) * M (from the triangle)
    
    Double_t triangle_factor[binsDeta]={0};

    TH2D *s_over_m[multclass];
    TH1D *s_m_deta[multclass];
    TH2D *s_over_m_x_m[multclass];
    
    for(Int_t i=0;i<multclass;i++){
        s_over_m[i] = (TH2D*)fileData->Get(Form("dphi_0_0_%d",i));
        s_m_deta[i] = (TH1D*)s_over_m[i]->ProjectionY()->Clone();
        s_over_m_x_m[i] = (TH2D*)s_over_m[i]->Clone();
        s_over_m_x_m[i]->Reset();
    }
    
    
    TF1 *f2 = new TF1("f2",fline,min_eta,max_eta,4);
    
    f2->FixParameter(0,1);
    f2->FixParameter(1,1/max_eta);
    f2->FixParameter(2,1);
    f2->FixParameter(3,-1/max_eta);
    
    for(Int_t i=0;i<binsDeta;i++){
        
        triangle_factor[i] = f2->Eval(s_m_deta[0]->GetBinCenter(i+1));

    }
    


    //--scale each deta bin of the old TH2 with the triangle_factor[deta]
    
    for(Int_t i=0;i<multclass;i++){
        for(Int_t j=0;j<binsDeta;j++){
            for(Int_t k=0;k<bins;k++){
                    s_over_m_x_m[i] -> SetBinContent(k+1,j+1,(s_over_m[i]->GetBinContent(k+1,j+1))*triangle_factor[j]);
                    s_over_m_x_m[i]->SetBinError(k+1,j+1,(s_over_m[i]->GetBinError(k+1,j+1))*triangle_factor[j]);
            }
        }
    }
    
    fileDataEMremoved = TFile::Open("dphi_corr_MEremoved.root","UPDATE");
    
    for(Int_t i=0;i<multclass;i++){
        
        s_over_m_x_m[i]->SetName(Form("dphiNoMixed_%d",i));
        s_over_m_x_m[i]->Write();
        
    }
    
    

    //-------------- DOUBLE RIDGE SUBTRACTION: gets rid of no-jet related components (v3 is still kept => effect added to the systematics) 
    
    // the ridge, estimated via an etagap, has to be scaled since it sits on the triangle 
    Double_t scale_for_ridge_NS = 0, scale_for_ridge_AS = 0;
    
        
    scale_for_ridge_NS = f2->Integral(min_bin_for_etagap,max_bin_for_etagap)/(f2->Integral(min_eta,min_bin_for_etagap)+f2->Integral(max_bin_for_etagap,max_eta)); //there is etagap in the NS
    cout<<"scaling NS:"<<scale_for_ridge_NS<<endl;
        
    scale_for_ridge_AS = f2->Integral(min_eta,max_eta)/(f2->Integral(min_eta,min_bin_for_etagap)+f2->Integral(max_bin_for_etagap,max_eta)); // there is no etagap in the AS
    cout<<"scaling AS:"<<scale_for_ridge_AS<<endl;
    
  // Double ridge subtraction
    
  TCanvas *c = new TCanvas();
  c->Divide(5,4);
  
  for(Int_t i=0;i<multclass;i++){
  c->cd(i+1);
        
        
      dphideta[i] = (TH2D*)fileDataEMremoved->Get(Form("dphiNoMixed_%d",i));

        
      // phi and eta projections
      fDeltaPhiNch[i] = (TH1D*)dphideta[i]->ProjectionX()->Clone();
      if(!zyam)
          fDeltaPhiNch[i]->Scale(binWidth);    //gaussians include the binwidth, so when using the fit, the histograms must be scaled first
      fDeltaPhiNch[i]->Draw();
    
      fDeltaEtaNch[i] = (TH1D*)dphideta[i]->ProjectionY()->Clone();
    
      // signal NS: |DEta|<max_bin_for_etagap; signal AS: |DEta|<max_eta
      fSignalNSDPhi[i] = (TH1D*)dphideta[i]->ProjectionX(Form("|DEta|<%f",max_bin_for_etagap),fDeltaEtaNch[i]->FindBin(min_bin_for_etagap+0.0001),fDeltaEtaNch[i]->FindBin(max_bin_for_etagap-0.0001))->Clone();
      fSignalASDPhi[i] = (TH1D*)dphideta[i]->ProjectionX(Form("|DEta|<%f",max_eta))->Clone();
      
      fSignalDPhi[i] = (TH1D*)fSignalASDPhi[i]->Clone();
      fSignalDPhi[i]->Reset();
      fSignalDPhi[i]->Sumw2();
      
      for(Int_t k=0;k<bins/2;k++){
          fSignalDPhi[i]->SetBinContent(k+1,fSignalNSDPhi[i]->GetBinContent(k+1));
          fSignalDPhi[i]->SetBinError(k+1, fSignalNSDPhi[i]->GetBinError(k+1));
      }
      for(Int_t k=bins/2;k<bins;k++){
          fSignalDPhi[i]->SetBinContent(k+1,fSignalASDPhi[i]->GetBinContent(k+1));
          fSignalDPhi[i]->SetBinError(k+1, fSignalASDPhi[i]->GetBinError(k+1));
      }
      if(!zyam)
          fSignalDPhi[i]->Scale(binWidth);
        
      // ridge1 DEta<min_bin_for_etagap
      fRidge1DPhi[i] = (TH1D*)dphideta[i]->ProjectionX(Form("DEta<%f",min_bin_for_etagap),1,fDeltaEtaNch[i]->FindBin(min_bin_for_etagap-0.0001))->Clone();
      if(!zyam)
          fRidge1DPhi[i]->Scale(binWidth);
      fRidge1DPhi[i]->SetMarkerColor(kRed);

      // ridge2 DEta>max_bin_for_etagap
      fRidge2DPhi[i] = (TH1D*)dphideta[i]->ProjectionX(Form("DEta>%f",max_bin_for_etagap),fDeltaEtaNch[i]->FindBin(max_bin_for_etagap+0.0001),fDeltaEtaNch[i]->GetNbinsX())->Clone();
      if(!zyam)
          fRidge2DPhi[i]->Scale(binWidth);
      fRidge2DPhi[i]->SetMarkerColor(kBlue);

      // ridge = ridge1 + ridge2
      fRidgeDPhi[i] = (TH1D*)fRidge1DPhi[i]->Clone("fRidge");
      fRidgeDPhi[i]->Reset();
      fRidgeDPhi[i]->Sumw2();
      fRidgeDPhi[i]->Add(fRidge1DPhi[i],fRidge2DPhi[i],1,1);
      //fRidgeDPhi[i]->Scale(scale_for_ridge);

      // symmetrize NS ridge in the AS
      fSymmRidgeNotScaled[i] = (TH1D*)fRidgeDPhi[i]->Clone("fSymmRidgeNotScaled");
      
      for(Int_t k=fSymmRidgeNotScaled[i]->GetNbinsX()/2+1;k<=fSymmRidgeNotScaled[i]->GetNbinsX();k++){
 
          fSymmRidgeNotScaled[i]->SetBinContent(k,fSymmRidgeNotScaled[i]->GetBinContent(fSymmRidgeNotScaled[i]->GetNbinsX()+1-k));

      }
      
      // scale the symmetrized ridge according to NS or AS
      fSymmRidge[i] = (TH1D*)fSymmRidgeNotScaled[i]->Clone("fSymmRidge");

      for(Int_t k=0;k<bins/2;k++){
          fSymmRidge[i]->SetBinContent(k+1,(fSymmRidgeNotScaled[i]->GetBinContent(k+1))*scale_for_ridge_NS);
      }
      for(Int_t k=bins/2;k<bins;k++){
          fSymmRidge[i]->SetBinContent(k+1,(fSymmRidgeNotScaled[i]->GetBinContent(k+1))*scale_for_ridge_AS);
      }

      
      // signal - symmetric ridge
      
      if(zyam){
          fFinal1DPhi[i] = new TH1D(Form("fFinal1DPhi[%d]",i),Form("fFinal1DPhi[%d]",i),bins,-0.5*TMath::Pi(),1.5*TMath::Pi());
          fFinal1DPhi[i]->Add(fSignalDPhi[i],fSymmRidge[i],1,-1);
          fFinal1DPhi[i]->Sumw2();
          fFinalDPhi[i] = (TH1D*)fFinal1DPhi[i]->Clone("fFinal"); // zyam: average between the two min values => sum first half of NS in the second half and second half of AS in the first half, so zyam = min/2
          fFinalDPhi[i]->Reset();
          fFinalDPhi[i]->Sumw2();
      
          for(Int_t k=1;k<=bins/4;k++){
              fFinalDPhi[i]->SetBinContent(k,0.);
              fFinalDPhi[i]->SetBinContent(k+bins/4,fFinal1DPhi[i]->GetBinContent(k+bins/4)+fFinal1DPhi[i]->GetBinContent(bins/4+1-k));
              fFinalDPhi[i]->SetBinError(k+bins/4,TMath::Sqrt(pow(fFinal1DPhi[i]->GetBinError(k+bins/4),2)+pow(fFinal1DPhi[i]->GetBinError(bins/4+1-k),2)));
              fFinalDPhi[i]->SetBinContent(k+bins/2,fFinal1DPhi[i]->GetBinContent(k+bins/2)+fFinal1DPhi[i]->GetBinContent(bins+1-k));
              fFinalDPhi[i]->SetBinError(k+bins/2,TMath::Sqrt(pow(fFinal1DPhi[i]->GetBinError(k+bins/2),2)+pow(fFinal1DPhi[i]->GetBinError(bins+1-k),2)));
              fFinalDPhi[i]->SetBinContent(k+bins/4*3,0.);
          
          }
      }
      
      else{

          fFinalDPhi[i] = (TH1D*)fSignalDPhi[i]->Clone();
          fFinalDPhi[i]->Reset();
          fFinalDPhi[i]->Sumw2();
          fFinalDPhi[i]->Add(fSignalDPhi[i],fSymmRidge[i],1,-1);
      }
      
  }

  // store the pair yields in a file (the yields are *not* normalized to the Ntriggers)
    
  TFile* file_yields = 0x0;
  if(zyam)
      file_yields = TFile::Open("PairYields_zyam.root","RECREATE");
  else
      file_yields = TFile::Open("PairYields_fit.root","RECREATE");


  for(Int_t i=0;i<multclass;i++){
      fDeltaEtaNch[i]->SetName(Form("DeltaEta_0_0_%d",i));
      fDeltaEtaNch[i]->Write();
      fDeltaPhiNch[i]->SetName(Form("Correlation bin %d in dphi",i));
      fDeltaPhiNch[i]->Write();
      fSignalDPhi[i]->SetName(Form("Signal_0_0_%d",i));
      fSignalDPhi[i]->Write();
      fRidgeDPhi[i]->SetName(Form("Ridge_0_0_%d",i));
      fRidgeDPhi[i]->Write();
      fSymmRidgeNotScaled[i]->SetName(Form("Symmetric_Ridge_NotScaled_0_0_%d",i));
      fSymmRidgeNotScaled[i]->Write();
      fSymmRidge[i]->SetName(Form("Symmetric_Ridge_0_0_%d",i));
      fSymmRidge[i]->Write();
      fFinalDPhi[i]->SetName(Form("Pure_Signal_0_0_%d",i));
      fFinalDPhi[i]->Write();
  }
  file_yields->Close();

  //-------------- CORRELATION OBSERVABLES: per-trigger yields, triggers and uncorrelated seeds
    
  Float_t baseline[multclass]={0};
  
  TGraphErrors *fNearSideIntegral = new TGraphErrors();
  fNearSideIntegral->SetName("fNearSideIntegral");
  fNearSideIntegral->SetMarkerColor(kGreen+2);
  fNearSideIntegral->SetLineColor(kGreen+2);
  fNearSideIntegral->SetLineWidth(1);
  fNearSideIntegral->SetMarkerStyle(4);

  TGraphErrors *fAwaySideIntegral = new TGraphErrors();
  fAwaySideIntegral->SetName("fAwaySideIntegral");
  fAwaySideIntegral->SetMarkerColor(kBlue);
  fAwaySideIntegral->SetLineColor(kBlue);
  fAwaySideIntegral->SetLineWidth(1);
  fAwaySideIntegral->SetMarkerStyle(4);

  TGraphErrors *fBothSideIntegral = new TGraphErrors();
  fBothSideIntegral->SetName("fBothSideIntegral");
  fBothSideIntegral->SetMarkerColor(kMagenta);
  fBothSideIntegral->SetLineColor(kMagenta);
  fBothSideIntegral->SetLineWidth(1);
  fBothSideIntegral->SetMarkerStyle(4);

    
  TGraphErrors *fNjets = new TGraphErrors();
  fNjets->SetName("fNjets");
  fNjets->SetMarkerColor(kCyan+2);
  fNjets->SetLineColor(kCyan+2);
  fNjets->SetLineWidth(1);
  fNjets->SetMarkerStyle(4);

  TGraphErrors *fTriggerAverage = new TGraphErrors();
  fTriggerAverage->SetName("fTriggerAverage");
  fTriggerAverage->SetMarkerColor(kBlack);
  fTriggerAverage->SetLineColor(kBlack);
  fTriggerAverage->SetLineWidth(1);
  fTriggerAverage->SetMarkerStyle(4);

  Int_t points=0;
  Double_t minbin[multclass] = {0};
  
  //  extract information out of dphi histograms
  TCanvas * cYields= new TCanvas("cYields", "cYields", 150, 150, 820, 620);
  cYields->Divide(5,4);
    
  for(Int_t i=0;i<multclass;i++){
  cYields->cd(i+1);
      

  if(zyam) {
      
      if(fFinalDPhi[i]->Integral()>0){
          fFinalDPhi[i]->GetXaxis()->SetRange(bins/4+1,bins/4*3);
          baseline[i]=fFinalDPhi[i]->GetMinimum()/2;
          minbin[i] = fFinalDPhi[i]->GetMinimumBin();
          fFinalDPhi[i]->GetXaxis()->UnZoom();
          
          for(Int_t k=0;k<bins;k++){
              if(fFinalDPhi[i]->GetBinContent(k+1)!=0)
                  fFinalDPhi[i]->SetBinContent(k+1,fFinalDPhi[i]->GetBinContent(k+1)-baseline[i]);
              else
                  fFinalDPhi[i]->SetBinContent(k+1,0.);
          }
          
          fFinalDPhi[i]->DrawClone("");
          
          fFinalDPhi[i]->SetTitle(Form("0.7<p_{T,trig}<5.0 - 0.7<p_{T,assoc}<5.0 - %d-%d %",i*5,(i+1)*5));
          fFinalDPhi[i]->SetTitle("1/N_{trig} dN_{assoc}/d#Delta#varphi (rad^{-1})");          
          //-
          Double_t errorNS = 0;
          Double_t nearSideResult = (fFinalDPhi[i]->IntegralAndError(0,minbin[i],errorNS,"width"))/trigger->GetBinContent(i+1);
          Double_t nearSideError = errorNS/trigger->GetBinContent(i+1); 
          fNearSideIntegral->SetPoint(points,i, nearSideResult);
          fNearSideIntegral->SetPointError(points,0.5,errorNS/trigger->GetBinContent(i+1));
          //-
          
          //--
          Double_t errorAS = 0;
          Double_t awaySideResult = (fFinalDPhi[i]->IntegralAndError(minbin[i],bins,errorAS,"width"))/trigger->GetBinContent(i+1);
          Double_t awaySideError = errorAS/trigger->GetBinContent(i+1); 
          fAwaySideIntegral->SetPoint(points,i, awaySideResult );
          fAwaySideIntegral->SetPointError(points,0.5, errorAS/trigger->GetBinContent(i+1));
          //--
          
          //---
          Double_t bothSideResult = nearSideResult + awaySideResult;
          Double_t bothSideError = bothSideResult * TMath::Sqrt(pow(errorNS,2)+pow(errorAS,2))/trigger->GetBinContent(i+1);
          fBothSideIntegral->SetPoint(points,i, bothSideResult );
          fBothSideIntegral->SetPointError(points,0.5, bothSideError );      
          //---
          

          
      }
      else{
          fNearSideIntegral->SetPoint(points,i, 0);
          fAwaySideIntegral->SetPoint(points,i, 0);
          fBothSideIntegral->SetPoint(points,i,0);
      }
      Double_t p0BinContent=p0->GetBinContent(i+1);
      Double_t p0BinError=p0->GetBinError(i+1);
      
      //--------
      Double_t njets =  p0BinContent/(1+bothSideResult); 
      Double_t njetsError = njets*TMath::Sqrt(bothSideError*bothSideError/(1+bothSideResult)/(1+bothSideResult)+p0BinError*p0BinError/p0BinContent/p0BinContent);
      fNjets->SetPoint(points,i, njets );
      fNjets->SetPointError(points,0.5,njetsError );
      
      //-------
      
      fTriggerAverage->SetPoint(points,i, p0BinContent);
      fTriggerAverage->SetPointError(points,0.5, p0BinError);
      
  }
      
  else if (!zyam){ 

      if(fFinalDPhi[i]->Integral()>0){

          //first fit function: 2 gauss + const
          fTotal2Gaus[i] = new TF1(Form("gaus3and2_%d",i), fitFunction2Gaus , -0.5*TMath::Pi(), 1.5*TMath::Pi(), 5);
          fTotal2Gaus[i]->SetName(Form("gaus3_%d",i));
          fTotal2Gaus[i]->SetParNames ("A1","sigma1","A3", "sigma3");
          fTotal2Gaus[i]->SetLineColor(kRed);
          fTotal2Gaus[i]->SetLineWidth(2);
    
          baseline[i]=fFinalDPhi[i]->GetMinimum();
          Double_t integr_for_const_2 = fFinalDPhi[i]->Integral();
        
          fTotal2Gaus[i]->FixParameter(4,integr_for_const_2);
          fTotal2Gaus[i]->SetParameters( fFinalDPhi[i]->GetBinContent(fFinalDPhi[i]->GetXaxis()->FindFixBin(0)) - baseline[i] , 0.6 , fFinalDPhi[i]->GetBinContent(fFinalDPhi[i]->GetXaxis()->FindFixBin(TMath::Pi()))-baseline[i] , 0.6);
      
          fTotal2Gaus[i]->SetParLimits(0, 0, (fFinalDPhi[i]->GetBinContent(fFinalDPhi[i]->GetXaxis()->FindFixBin(0))-baseline[i])*2);
          fTotal2Gaus[i]->SetParLimits(1, 0.01, 10);
          fTotal2Gaus[i]->SetParLimits(2, 0, (fFinalDPhi[i]->GetBinContent(fFinalDPhi[i]->GetXaxis()->FindFixBin(TMath::Pi()))-baseline[i])*2);
          fTotal2Gaus[i]->SetParLimits(3, 0.01, 10);
      
          fTotal2Gaus[i]->SetLineColor(kRed);
          fTotal2Gaus[i]->SetLineWidth(2);

          fFinalDPhi[i]->Fit(fTotal2Gaus[i],flag);
          fFinalDPhi[i]->SetMinimum(0);
          fFinalDPhi[i]->DrawClone("");
          fTotal2Gaus[i] ->DrawClone("same");
        
          Double_t A11     = fTotal2Gaus[i]->GetParameter(0);
          Double_t sigma11 = fTotal2Gaus[i]->GetParameter(1);
          Double_t A31     = fTotal2Gaus[i]->GetParameter(2);
          Double_t sigma31 = fTotal2Gaus[i]->GetParameter(3);

          Double_t a1e1 = fTotal2Gaus[i]->GetParError(0);
          Double_t s1e1 = fTotal2Gaus[i]->GetParError(1);
          Double_t a3e1 = fTotal2Gaus[i]->GetParError(2);
          Double_t s3e1 = fTotal2Gaus[i]->GetParError(3);
        
      
          Double_t T11 = A11*sigma11; 
          Double_t T31 = A31*sigma31;
          Double_t t11 = T11*TMath::Sqrt(a1e1*a1e1/A11/A11 + s1e1*s1e1/sigma11/sigma11); 
          Double_t t31 = T31*TMath::Sqrt(a3e1*a3e1/A31/A31 + s3e1*s3e1/sigma31/sigma31);
  

          //second fit: 3 gauss + const
          fTotal[i] = new TF1(Form("gaus3_%d",i), fitFunction , -0.5*TMath::Pi(), 1.5*TMath::Pi(), 7);
          fTotal[i]->SetName(Form("gaus3_%d",i));
          fTotal[i]->SetParNames ("A1","sigma1","A2","sigma2", "A3", "sigma3","integral");
          fTotal[i]->SetLineColor(kRed);
          fTotal[i]->SetLineWidth(2);
    
          Double_t integr_for_const = fFinalDPhi[i]->Integral();
    
        
          fTotal[i]->FixParameter(0,A11);
          fTotal[i]->FixParameter(1,sigma11*1.2);
          fTotal[i]->FixParameter(2,A11);
          fTotal[i]->FixParameter(3,sigma11*0.7);
          fTotal[i]->FixParameter(4,A31);
          fTotal[i]->FixParameter(5,sigma31);
          fTotal[i]->FixParameter(6,integr_for_const);

          fTotal[i]->SetParLimits(0, 0, (fFinalDPhi[i]->GetBinContent(fFinalDPhi[i]->GetXaxis()->FindFixBin(0))-baseline[i])*2);
          fTotal[i]->SetParLimits(1, 0.3, 10); 
          fTotal[i]->SetParLimits(2, 0, (fFinalDPhi[i]->GetBinContent(fFinalDPhi[i]->GetXaxis()->FindFixBin(0))-baseline[i])*2);
          fTotal[i]->SetParLimits(3, 0.12, 0.4);
          fTotal[i]->SetParLimits(4, 0, (fFinalDPhi[i]->GetBinContent(fFinalDPhi[i]->
							      GetXaxis()->FindFixBin(TMath::Pi()))-baseline[i])*2);
          fTotal[i]->SetParLimits(5, 0.01, 10);
    
          fTotal[i]->SetLineColor(kRed);
          fTotal[i]->SetLineWidth(2);


          fFinalDPhi[i]->Fit(fTotal[i],flag);
          fFinalDPhi[i]->SetMinimum(0);
          fFinalDPhi[i]->DrawClone("");
          fFinalDPhi[i]->SetTitle(Form("0.7<p_{T,trig}<5.0 - 0.7<p_{T,assoc}<5.0 - %d-%d %",i*5,(i+1)*5));
          fFinalDPhi[i]->SetTitle("1/N_{trig} dN_{assoc}/d#Delta#varphi (rad^{-1})");
          fTotal[i]->DrawClone("same");
        
          Double_t A1     = fTotal[i]->GetParameter(0);
          Double_t sigma1 = fTotal[i]->GetParameter(1);
          Double_t A2     = fTotal[i]->GetParameter(2);
          Double_t sigma2 = fTotal[i]->GetParameter(3);
          Double_t A3     = fTotal[i]->GetParameter(4);
          Double_t sigma3 = fTotal[i]->GetParameter(5);

      
          //define each gaussian and constant to be drawn with different colors on top of each other
        
          TF1 * fConstant = new TF1("konst", "pol0(0)",-0.5*TMath::Pi(), 1.5*TMath::Pi());
          fConstant->SetParameter(0,(integr_for_const - TMath::Sqrt(TMath::Pi()*2)/binWidth*(A1*sigma1+A2*sigma2+A3*sigma3))/bins);
          fConstant->SetLineColor(kBlue);
          fConstant->Draw("same");
      
          //gaus 1 NS
          TF1 * fGaussian1 = new TF1("fGaussian1", "[0]*exp(-x*x/(2*[1]*[1])) +[0] * exp(-(x-TMath::TwoPi())*(x-TMath::TwoPi())/(2*[1]*[1]))",-0.5*TMath::Pi(), 1.5*TMath::Pi());
          fGaussian1->SetParameters(fTotal[i]->GetParameter(0),fTotal[i]->GetParameter(1));
          fGaussian1->SetLineColor(kMagenta);
          fGaussian1->SetLineStyle(1);
          fGaussian1->Draw("same");
      
          //gaus 2 NS
          TF1 * fGaussian2 = new TF1("fGaussian2", "[0]*exp(-x*x/(2*[1]*[1])) +[0] * exp(-(x-TMath::TwoPi())*(x-TMath::TwoPi())/(2*[1]*[1]))",-0.5*TMath::Pi(), 1.5*TMath::Pi());
          fGaussian2->SetLineColor(kGreen+2);
          fGaussian2->SetParameters(fTotal[i]->GetParameter(2),fTotal[i]->GetParameter(3));
          fGaussian2->Draw("same");
      
          //gaus 3 AS
          TF1 * fGaussian3 = new TF1("fGaussian3", "[0] * exp(-((x-TMath::Pi()))*((x-TMath::Pi()))/(2*[1]*[1]))+[0] * exp(-((x+TMath::Pi()))*((x+TMath::Pi()))/(2*[1]*[1]))",-0.5*TMath::Pi(), 1.5*TMath::Pi());
          fGaussian3->SetLineColor(kCyan);
          fGaussian3->SetParameters(fTotal[i]->GetParameter(4), fTotal[i]->GetParameter(5));
          fGaussian3->Draw("same");
        
        
          Double_t a1e = fTotal[i]->GetParError(0);
          Double_t s1e = fTotal[i]->GetParError(1);
          Double_t a2e = fTotal[i]->GetParError(2);
          Double_t s2e = fTotal[i]->GetParError(3);
          Double_t a3e = fTotal[i]->GetParError(4);
          Double_t s3e = fTotal[i]->GetParError(5);

          Double_t T1 = A1*sigma1;
          Double_t T2 = A2*sigma2;
          Double_t T3 = A3*sigma3;
          Double_t t1 = T1*TMath::Sqrt(a1e*a1e/A1/A1 + s1e*s1e/sigma1/sigma1);
          Double_t t2 = T2*TMath::Sqrt(a2e*a2e/A2/A2 + s2e*s2e/sigma2/sigma2);
          Double_t t3 = T3*TMath::Sqrt(a3e*a3e/A3/A3 + s3e*s3e/sigma3/sigma3);
              
          //-
          Double_t nearSideResult = TMath::Sqrt(TMath::Pi()*2)/ binWidth* (A1 * sigma1 + A2 * sigma2)/trigger->GetBinContent(i+1);
          Double_t nearSideError = nearSideResult * TMath::Sqrt((t1*t1 + t2*t2)/(T1+T2)/(T1+T2)+ 1./trigger->GetBinContent(i+1));
          fNearSideIntegral->SetPoint(points,i, nearSideResult);
          fNearSideIntegral->SetPointError(points,0.5,nearSideError);
        
          //-

          //--
          Double_t awaySideResult = TMath::Sqrt(TMath::Pi()*2)/ binWidth* 
        (A3 * sigma3)/trigger->GetBinContent(i+1);
          Double_t awaySideError = awaySideResult*TMath::Sqrt(a3e*a3e/A3/A3 + s3e*s3e/sigma3/sigma3 + 1/trigger->GetBinContent(i+1));
          fAwaySideIntegral->SetPoint(points,i, awaySideResult );
          fAwaySideIntegral->SetPointError(points,0.5, awaySideError );         
          //--

          //---
          bothSideResult = TMath::Sqrt(TMath::Pi()*2)/ binWidth* (A1 * sigma1 + A2 * sigma2 + A3 * sigma3 )/trigger->GetBinContent(i+1); 
          bothSideError = nearSideResult *  TMath::Sqrt((t1*t1 + t2*t2 + t3*t3)/(T1+T2+T3)/(T1+T2+T3)+ 1./trigger->GetBinContent(i+1));
          fBothSideIntegral->SetPoint(points,i, bothSideResult );
          fBothSideIntegral->SetPointError(points,0.5, bothSideError );      
          //---
            
    }
    else{
        
        fNearSideIntegral->SetPoint(points,i, 0);
        fAwaySideIntegral->SetPoint(points,i, 0);
        fBothSideIntegral->SetPoint(points,i,0);
    
    }
    Double_t p0BinContent=p0->GetBinContent(i+1);
    Double_t p0BinError=p0->GetBinError(i+1);
    
    //--------
    Double_t njets =  p0BinContent/(1+bothSideResult); 
    Double_t njetsError = njets*TMath::Sqrt(bothSideError*bothSideError/(1+bothSideResult)/(1+bothSideResult) + p0BinError*p0BinError/p0BinContent/p0BinContent);
    fNjets->SetPoint(points,i, njets );
    fNjets->SetPointError(points,0.5,njetsError );
    //-------
      
    fTriggerAverage->SetPoint(points,i, p0BinContent);
    fTriggerAverage->SetPointError(points,0.5, p0BinError);
      
    
  }
      points++;
  }


  TFile* file = 0x0;
  if(zyam)
      file = TFile::Open("njet_zyam.root","RECREATE");
  else
      file = TFile::Open("njet_fit.root","RECREATE");

  fNearSideIntegral->Write();
  fAwaySideIntegral->Write();
  fBothSideIntegral->Write();
  fNjets->Write();
  fTriggerAverage->Write();

  file->Close();



}


