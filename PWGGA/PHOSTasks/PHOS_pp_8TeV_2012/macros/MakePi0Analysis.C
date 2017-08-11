#include <TFile.h>
#include <TList.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <Stdio.h>
#include <stdlib.h>
#include <TGaxis.h>
#include <THnSparse.h>

using namespace std;

TFile *tfout = new TFile("pi0.root","recreate");

void MakePi0Analysis(){

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetPadBottomMargin(0.11);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(1);
  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);
  TGaxis::SetMaxDigits(3);
  
  Double_t bins[] = {1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0, 3.2, 3.4, 3.6,3.8,4.0,4.5, 5.0, 5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,10.0,11,12}; //Nicolas LHC12d PHI
  const Int_t binnum = sizeof(bins)/sizeof(Double_t) - 1;
  
  TString cwd = gSystem->WorkingDirectory();

  TH1F* fHistPeakMean  = new TH1F("fHistPeakMean","",binnum,bins);
  TH1F* fHistPeakWidth = new TH1F("fHistPeakWidth","",binnum,bins);
  TH1F* fHistRawYield  = new TH1F("fHistRawYield","",binnum,bins);
  TH1F* fHistRawYieldPerEvent  = new TH1F("fHistRawYieldPerEvent","",binnum,bins);
  
  TCanvas *c1[binnum];
  TCanvas *c2[binnum];
  TCanvas *c3;

  TFile* ef     = TFile::Open("../AnalysisResults.root");
  TList* list   = (TList*)ef->Get("list_kINT7_Pi0");
  TList *listEv = (TList*)ef->Get("list_kINT7_Event");
  
  TH1F* fHistEvents = (TH1F*)listEv->FindObject("fHistAnalyzedEvents");
  Double_t nEvt = fHistEvents->GetEntries();
  
  TString tofName[] = {"","_TOFcut1","_TOFcut2"};
  TString fitName[] = {"CrystalBall","AsymmGauss"};
  TString modName[] = {"","_M1","_M3"};
  
  Int_t ntof = 3;
  
  Double_t signal_range_min=0.11;
  Double_t signal_range_max=0.15;
  Double_t bg_range_min=0.050;
  Double_t bg_range_max=0.250;
  
  for(Int_t itof=0; itof<ntof; ++itof){
    for(Int_t iMod=0; iMod<1; iMod++){
      for(Int_t iFit=0; iFit<2; iFit++){
	for(Int_t iPol=0; iPol<3; iPol++){
	  
	  Int_t pol = iPol;

	  TString fNameSame = "fHistMassTwoGammas"+tofName[itof]+modName[iMod];
	  TString fNameMix  = "fHistMixMassTwoGammas"+tofName[itof]+modName[iMod];
	  
	  THnSparse *fTHnSparseRecPi0;
	  THnSparse *fTHnSparseRecMixPi0;
	  THnSparse *fTHnSparseGeneratePi0;
	  
	  TH2* hPi0A08    = 0x0;
	  TH2* hMixPi0A08 = 0x0;
	  hPi0A08    = (TH2*)list->FindObject(fNameSame)->Clone();
	  hMixPi0A08 = (TH2*)list->FindObject(fNameMix)->Clone();

	  Int_t sbin = 0;
	  Int_t lbin = 0;
      
	  TH1F* fHistFinalRatio[binnum];
	  TH1F* fHistFinalSame[binnum];
	  TH1F* fHistFinalBG[binnum];
	  TH1F* fHistFinalSignal[binnum];
	  TH1F* fHistOnlyFitSignal[binnum];
	  
	  TF1* fFitFinalRatio[binnum];
	  TF1* fFitFinalSignal[binnum];
	  TF1* fFitOnlyFitSignal[binnum];

	  for(Int_t i=0; i<binnum; i++){    
	    
	    Double_t bin_width = (bins[i+1]-bins[i]);
	    Double_t bin_mean  = (bins[i+1]+bins[i])/2;
	    
	    Int_t rebin = 1;
	    
	    if(bin_mean>5.0){
	      rebin=10;
	    }
	    else{
	      rebin=5;
	    }
	    
	    Int_t sproj_bin = hPi0A08->GetYaxis()->FindBin(bins[i]);
	    Int_t lproj_bin = hPi0A08->GetYaxis()->FindBin(bins[i+1])-1;
	    
	    TH1* fHistBasicSame = (TH1*)hPi0A08   ->ProjectionX(Form("fHistBasicSame_No%d",i+1),sproj_bin,lproj_bin,""); 
	    TH1* fHistBasicMix  = (TH1*)hMixPi0A08->ProjectionX(Form("fHistBasicMix_No%d",i+1),sproj_bin,lproj_bin,""); 
	    fHistBasicSame->Rebin(rebin);
	    fHistBasicMix->Rebin(rebin);
	    fHistBasicSame->Sumw2();
	    fHistBasicMix->Sumw2();
	    fHistBasicSame->GetYaxis()->SetTitle(Form("dN/dM per %.0f MeV/c^{2}",fHistBasicSame->GetBinWidth(1)*1000));
	    fHistBasicMix->GetYaxis()->SetTitle(Form("dN/dM per %.0f MeV/c^{2}",fHistBasicSame->GetBinWidth(1)*1000));
	    fHistBasicSame->GetXaxis()->SetRange(fHistBasicSame->GetXaxis()->FindBin(0.),fHistBasicSame->GetXaxis()->FindBin(0.3)-1);
	    fHistBasicMix->GetXaxis()->SetRange(fHistBasicSame->GetXaxis()->FindBin(0.),fHistBasicSame->GetXaxis()->FindBin(0.3)-1);
	    fHistBasicSame->SetMarkerStyle(20);
	    fHistBasicMix->SetMarkerStyle(24);
	    fHistBasicSame->SetMarkerColor(kBlack);
	    fHistBasicMix->SetMarkerColor(kBlue);
	    fHistBasicSame->SetLineColor(kBlack);
	    fHistBasicMix->SetLineColor(kBlue);
	    
	    fHistFinalSame[i]     = (TH1F*)fHistBasicSame->Clone();
	    fHistOnlyFitSignal[i] = (TH1F*)fHistBasicSame->Clone();
	    fHistOnlyFitSignal[i]->SetName(Form("fHistOnlyFitSignal_Pol%d_No%d",pol,i+1)+fitName[iFit]+tofName[itof]+modName[iMod]);


	    TH1F* fHistOnlyFit_Signal = (TH1F*)fHistOnlyFitSignal[i]->Clone();
	    TH1F* fHistOnlyFit_BG     = (TH1F*)fHistOnlyFitSignal[i]->Clone();
	    
	    TH1F* fHistRatio     = (TH1F*)fHistBasicSame->Clone();
	    fHistRatio->Divide(fHistBasicMix);
	    fHistRatio->SetName(Form("cRatioSameBG_Pol%d_No%d",pol,i+1)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    
	    TH1F* fHistRatio_Signal = (TH1F*)fHistRatio->Clone();
	    TH1F* fHistRatio_BG     = (TH1F*)fHistRatio->Clone();
	    
	    Int_t ssignal_bin = fHistRatio_Signal->GetXaxis()->FindBin(signal_range_min);
	    Int_t lsignal_bin = fHistRatio_Signal->GetXaxis()->FindBin(signal_range_max);
	    Int_t sbg_bin = fHistRatio_BG->GetXaxis()->FindBin(bg_range_min);
	    Int_t lbg_bin = fHistRatio_BG->GetXaxis()->FindBin(bg_range_max);
	    
	    for(Int_t j=ssignal_bin; j<lsignal_bin; ++j){
	      fHistRatio_BG->SetBinContent(j,0);
	      fHistRatio_BG->SetBinError(j,0);
	      fHistOnlyFit_BG->SetBinContent(j,0);
	      fHistOnlyFit_BG->SetBinError(j,0);
	    }
	    for(Int_t j=sbg_bin; j<ssignal_bin; ++j){
	      fHistRatio_Signal->SetBinContent(j,0);
	      fHistRatio_Signal->SetBinError(j,0);
	      fHistOnlyFit_Signal->SetBinContent(j,0);
	      fHistOnlyFit_Signal->SetBinError(j,0);
	    }
	    for(Int_t j=lsignal_bin; j<lbg_bin; ++j){
	      fHistRatio_Signal->SetBinContent(j,0);
	      fHistRatio_Signal->SetBinError(j,0);
	      fHistOnlyFit_Signal->SetBinContent(j,0);
	      fHistOnlyFit_Signal->SetBinError(j,0);
	    }
	    
	    ///////////////////////////////////////////////////////////////////////////////////////////////////
	    // Combinatrial background analysis
	    ///////////////////////////////////////////////////////////////////////////////////////////////////
	    
	    TF1 *fFitRatio_Signal = NULL;
	    if(iFit==0){
	      fFitRatio_Signal = new TF1("fFitRatio_Signal",
					 "[4]*((x-[2])/[3] > -[0] ? 1.0:0.0 )*exp(-pow(x-[2],2)/(2*pow([3],2))) + ((x-[2])/[3] <= -[0] ? 1.0:0.0 )*[4]*pow([1]/[0],[1])*exp(-[0]*[0]/2) * pow([1]/[0]-[0]-(x-[2])/[3],-[1])+[5]",
					 signal_range_min,signal_range_max);
	      fFitRatio_Signal->SetParameters(1.,6.,0.135,0.005,0.01);
	      fFitRatio_Signal->FixParameter(0,1.16053);
	      fFitRatio_Signal->FixParameter(1,6.);
	      fFitRatio_Signal->SetParLimits(2,0.130,0.140);
	      fFitRatio_Signal->SetParLimits(3,0.005,0.03);
	    }
	    else if(iFit==1){//[0]=A, [1]=M_pi0, [2]=sigma, [3]=lamda
	      fFitRatio_Signal = new TF1("fFitRatio_Signal",
					 "[0] *( ([1]>x ? 1.0:0.0 ) * (1 - exp(-0.5*pow((x-[1])/[2],2))) * exp((x-[1])/[3]) + exp(-0.5*pow((x-[1])/[2],2)) )",
					 signal_range_min,signal_range_max);
	      fFitRatio_Signal->SetParameters(1,0.135,6.46624e-03,7.47626e-03);
	      fFitRatio_Signal->SetParLimits(1,0.130,0.140);
	      fFitRatio_Signal->SetParLimits(2,0.001,0.03);
	      fFitRatio_Signal->SetParLimits(3,0.001,0.03);
	    }
	    
	    fHistRatio_Signal->Fit(fFitRatio_Signal,"RNQ");
	    fHistRatio_Signal->Fit(fFitRatio_Signal,"RNQ");
	    fHistRatio_Signal->Fit(fFitRatio_Signal,"RNQ");
	    
	    TF1 *fFitRatio_BG = 0x0;
	    if(pol==0){
	      fFitRatio_BG = new TF1("fFitRatio_BG","[0]",0,1);
	    }
	    else if(pol==1){
	      fFitRatio_BG = new TF1("fFitRatio_BG","[0]+[1]*x",0,1);
	    }
	    else if(pol==2){
	      fFitRatio_BG = new TF1("fFitRatio_BG","[0]+[1]*x+[2]*x*x",0,1);
	    }
	    fHistRatio_BG->Fit(fFitRatio_BG,"RNQ","",bg_range_min,bg_range_max);
	    fHistRatio_BG->Fit(fFitRatio_BG,"RNQ","",bg_range_min,bg_range_max);
	    fHistRatio_BG->Fit(fFitRatio_BG,"RNQ","",bg_range_min,bg_range_max);
	    
	    TF1 *fFitRatio_SignalBG = 0x0;
	    TF1 *fFitScale = 0x0;
	    
	    if(iFit==0){
	      if(pol==0){
		fFitRatio_SignalBG= new TF1("fFitRatio_SignalBG",
					    "[4]*((x-[2])/[3] > -[0] ? 1.0:0.0 )*exp(-pow(x-[2],2)/(2*pow([3],2))) + ((x-[2])/[3] <= -[0] ? 1.0:0.0 )*[4]*pow([1]/[0],[1])*exp(-[0]*[0]/2) * pow([1]/[0]-[0]-(x-[2])/[3],-[1])+[5]",0,1);
		fFitRatio_SignalBG->SetParameter(0,fFitRatio_Signal->GetParameter(0));
		fFitRatio_SignalBG->SetParameter(1,fFitRatio_Signal->GetParameter(1));
		fFitRatio_SignalBG->SetParameter(2,fFitRatio_Signal->GetParameter(2));
		fFitRatio_SignalBG->SetParameter(3,fFitRatio_Signal->GetParameter(3));
		fFitRatio_SignalBG->SetParameter(4,fFitRatio_Signal->GetParameter(4));
		fFitRatio_SignalBG->SetParameter(5,fFitRatio_BG->GetParameter(0));
		fFitScale = new TF1("fFitScale","pol0",0.,1);
	      }
	      else if(pol==1){
		fFitRatio_SignalBG = new TF1("fFitRatio_SignalBG",
					     "[4]*((x-[2])/[3] > -[0] ? 1.0:0.0 )*exp(-pow(x-[2],2)/(2*pow([3],2))) + ((x-[2])/[3] <= -[0] ? 1.0:0.0 )*[4]*pow([1]/[0],[1])*exp(-[0]*[0]/2) * pow([1]/[0]-[0]-(x-[2])/[3],-[1])+[5]+[6]*x",
					     0,1);
		fFitRatio_SignalBG->SetParameter(0,fFitRatio_Signal->GetParameter(0));
		fFitRatio_SignalBG->SetParameter(1,fFitRatio_Signal->GetParameter(1));
		fFitRatio_SignalBG->SetParameter(2,fFitRatio_Signal->GetParameter(2));
		fFitRatio_SignalBG->SetParameter(3,fFitRatio_Signal->GetParameter(3));
		fFitRatio_SignalBG->SetParameter(4,fFitRatio_Signal->GetParameter(4));
		fFitRatio_SignalBG->SetParameter(5,fFitRatio_BG->GetParameter(0));
		fFitScale = new TF1("fFitScale","pol1",0,1);
	      }
	      else if(pol==2){
		fFitRatio_SignalBG = new TF1("fFitRatio_SignalBG",
					     "[4]*((x-[2])/[3] > -[0] ? 1.0:0.0 )*exp(-pow(x-[2],2)/(2*pow([3],2))) + ((x-[2])/[3] <= -[0] ? 1.0:0.0 )*[4]*pow([1]/[0],[1])*exp(-[0]*[0]/2) * pow([1]/[0]-[0]-(x-[2])/[3],-[1])+[5]+[6]*x+[7]*x*x",0,1);
		fFitRatio_SignalBG->SetParameter(0,fFitRatio_Signal->GetParameter(0));
		fFitRatio_SignalBG->SetParameter(1,fFitRatio_Signal->GetParameter(1));
		fFitRatio_SignalBG->SetParameter(2,fFitRatio_Signal->GetParameter(2));
		fFitRatio_SignalBG->SetParameter(3,fFitRatio_Signal->GetParameter(3));
		fFitRatio_SignalBG->SetParameter(4,fFitRatio_Signal->GetParameter(4));
		fFitRatio_SignalBG->SetParameter(5,fFitRatio_BG->GetParameter(0));
		fFitRatio_SignalBG->SetParameter(6,fFitRatio_BG->GetParameter(1));
		fFitRatio_SignalBG->SetParameter(7,fFitRatio_BG->GetParameter(2));
		fFitScale = new TF1("fFitScale","pol2",0,1);
	      }
	      
	      fFitRatio_SignalBG->FixParameter(0,1.16053);
	      fFitRatio_SignalBG->FixParameter(1,6.);
	      fFitRatio_SignalBG->SetParLimits(2,0.130,0.140);
	      fFitRatio_SignalBG->SetParLimits(3,0.005,0.03);
	    }
	    
	    if(iFit==1){
	      if(pol==0){
		fFitRatio_SignalBG= new TF1("fFitRatio_SignalBG",
					    "[0] *( ([1]>x ? 1.0:0.0 ) * (1 - exp(-0.5*pow((x-[1])/[2],2))) * exp((x-[1])/[3]) + exp(-0.5*pow((x-[1])/[2],2)) ) + [4]",
					    0,1);
		fFitRatio_SignalBG->SetParameter(0,fFitRatio_Signal->GetParameter(0));
		fFitRatio_SignalBG->SetParameter(1,fFitRatio_Signal->GetParameter(1));
		fFitRatio_SignalBG->SetParameter(2,fFitRatio_Signal->GetParameter(2));
		fFitRatio_SignalBG->SetParameter(3,fFitRatio_Signal->GetParameter(3));
		fFitRatio_SignalBG->SetParameter(4,fFitRatio_BG->GetParameter(0));
		fFitScale = new TF1("fFitScale","pol0",0,1);
	      }
	      else if(pol==1){
		fFitRatio_SignalBG = new TF1("fFitRatio_SignalBG",
					     "[0] *( ([1]>x ? 1.0:0.0 ) * (1 - exp(-0.5*pow((x-[1])/[2],2))) * exp((x-[1])/[3]) + exp(-0.5*pow((x-[1])/[2],2)) ) + [4]+[5]*x",
					     0,1);
		fFitRatio_SignalBG->SetParameter(0,fFitRatio_Signal->GetParameter(0));
		fFitRatio_SignalBG->SetParameter(1,fFitRatio_Signal->GetParameter(1));
		fFitRatio_SignalBG->SetParameter(2,fFitRatio_Signal->GetParameter(2));
		fFitRatio_SignalBG->SetParameter(3,fFitRatio_Signal->GetParameter(3));
		fFitRatio_SignalBG->SetParameter(4,fFitRatio_BG->GetParameter(0));
		fFitRatio_SignalBG->SetParameter(5,fFitRatio_BG->GetParameter(1));
		fFitScale = new TF1("fFitScale","pol1",0,1);
	      }
	      else if(pol==2){
		fFitRatio_SignalBG = new TF1("fFitRatio_SignalBG",
					     "[0] *( ([1]>x ? 1.0:0.0 ) * (1 - exp(-0.5*pow((x-[1])/[2],2))) * exp((x-[1])/[3]) + exp(-0.5*pow((x-[1])/[2],2)) ) + [4]+[5]*x+[6]*x*x",
					     0,1);
		fFitRatio_SignalBG->SetParameter(0,fFitRatio_Signal->GetParameter(0));
		fFitRatio_SignalBG->SetParameter(1,fFitRatio_Signal->GetParameter(1));
		fFitRatio_SignalBG->SetParameter(2,fFitRatio_Signal->GetParameter(2));
		fFitRatio_SignalBG->SetParameter(3,fFitRatio_Signal->GetParameter(3));
		fFitRatio_SignalBG->SetParameter(4,fFitRatio_BG->GetParameter(0));
		fFitRatio_SignalBG->SetParameter(5,fFitRatio_BG->GetParameter(1));
		fFitRatio_SignalBG->SetParameter(6,fFitRatio_BG->GetParameter(2));
		fFitScale = new TF1("fFitScale","pol2",0,1);
	      }
	      fFitRatio_SignalBG->SetParLimits(1,0.130,0.140);
	      fFitRatio_SignalBG->SetParLimits(2,0.005,0.03);
	      fFitRatio_SignalBG->SetParLimits(3,0.005,0.03);
	    }

	    fHistRatio->Fit(fFitRatio_SignalBG,"NRQ","",bg_range_min,bg_range_max);
	    fHistRatio->Fit(fFitRatio_SignalBG,"NRQ","",bg_range_min,bg_range_max);
	    fHistRatio->Fit(fFitRatio_SignalBG,"NRQ","",bg_range_min,bg_range_max);
	    
	    fHistFinalRatio[i] = (TH1F*)fHistRatio->Clone();;
	    fFitFinalRatio[i]  = (TF1*)fFitRatio_SignalBG->Clone();
	    
	    if(iFit==0){
	      if(pol==0){
		fFitScale->SetParameter(0,fFitRatio_SignalBG->GetParameter(5));
	      }
	      else if(pol==1){
		fFitScale->SetParameters(fFitRatio_SignalBG->GetParameter(5),fFitRatio_SignalBG->GetParameter(6));
	      }
	      else if(pol==2){
		fFitScale->SetParameters(fFitRatio_SignalBG->GetParameter(5),fFitRatio_SignalBG->GetParameter(6),fFitRatio_SignalBG->GetParameter(7));
	      }
	    }	
	    if(iFit==1){
	      if(pol==0){
		fFitScale->SetParameter(0,fFitRatio_SignalBG->GetParameter(4));
	      }
	      else if(pol==1){
		fFitScale->SetParameters(fFitRatio_SignalBG->GetParameter(4),fFitRatio_SignalBG->GetParameter(5));
	      }
	      else if(pol==2){
		fFitScale->SetParameters(fFitRatio_SignalBG->GetParameter(4),fFitRatio_SignalBG->GetParameter(5),fFitRatio_SignalBG->GetParameter(6));
	      }
	    }	
	    
	    c3 = new TCanvas("c3","c3",600,600);
	    c3->cd();
	    fHistRatio->Draw();
	    fFitRatio_SignalBG->Draw("same");
	    fFitScale->Draw("same");
	    c3->SetName(Form("cRatioSameBG_Pol%d_No%d",pol,i+1)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    tfout->WriteTObject(c3);

	    TH1F* fHistScaledMix = (TH1F*)fHistBasicMix->Clone();
	    fHistScaledMix->Multiply(fFitScale);
	    fHistFinalBG[i] = (TH1F*)fHistScaledMix->Clone();;
	    
	    TH1F* fHistSignal = (TH1F*)fHistBasicSame->Clone();
	    fHistSignal->SetName(Form("fHistSignal_Pol%d_No%d",pol,i+1)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    fHistScaledMix->SetName(Form("fHistScaledMix_Pol%d_No%d",pol,i+1)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    tfout->WriteTObject(fHistSignal);
	    tfout->WriteTObject(fHistScaledMix);
	    
	    for(Int_t j=0; j<fHistBasicSame->GetXaxis()->FindBin(0.3); ++j){
	      Double_t same   = fHistSignal->GetBinContent(j);
	      Double_t e_same = fHistSignal->GetBinError(j);
	      Double_t mix    = fHistScaledMix->GetBinContent(j);
	      Double_t e_mix  = fHistScaledMix->GetBinError(j);
	      
	      Double_t signal   = same - mix;
	      Double_t e_signal = sqrt(pow(e_same,2)+pow(e_mix,2));
	      
	      if(same>0){
		signal   = same - mix;
		e_signal = sqrt(pow(e_same,2)+pow(e_mix,2));
	      }
	      else{
		signal   = same;
		e_signal = e_same;
	      }
	      
	      fHistSignal->SetBinContent(j,signal);
	      fHistSignal->SetBinError(j,e_signal);
	    }
	    
	    fHistFinalSignal[i] = (TH1F*)fHistSignal->Clone();
	    fHistFinalSignal[i]->SetTitle(Form("%.2f < #it{p}_{T} %.2f (GeV/c)", bins[i], bins[i+1]));

	    TF1 *fFitSignal = NULL;
	    
	    if(iFit==0){
	      fFitSignal = new TF1("fFitSignal", "[4]*((x-[2])/[3] > -[0] ? 1.0:0.0 )*exp(-pow(x-[2],2)/(2*pow([3],2))) + ((x-[2])/[3] <= -[0] ? 1.0:0.0 )*[4]*pow([1]/[0],[1])*exp(-[0]*[0]/2) * pow([1]/[0]-[0]-(x-[2])/[3],-[1])",
				   0,0.3);
	      fFitSignal->SetParameters(1.,6.,0.135,0.005,0.01);

	      fFitSignal->FixParameter(0,1.16053);
	      fFitSignal->FixParameter(1,6.);

	      fFitSignal->SetParLimits(2,0.120,0.140);
	      fFitSignal->SetParLimits(3,0.005,0.03);
	    
	    }	
	    else if(iFit==1){
	      fFitSignal = new TF1("fFitSignal",
				   "[0] *( ([1]>x ? 1.0:0.0 ) * (1 - exp(-0.5*pow((x-[1])/[2],2))) * exp((x-[1])/[3]) + exp(-0.5*pow((x-[1])/[2],2)) )",
				   0,0.3);
	      fFitSignal->SetParameters(1,0.135,6.46624e-03,7.47626e-03);
	      fFitSignal->SetParLimits(1,0.125,0.140);
	      fFitSignal->SetParLimits(2,0.001,0.01);
	      fFitSignal->SetParLimits(3,0.001,0.01);
	    }
	    
	    fHistSignal->Fit(fFitSignal,"RNQ","",signal_range_min,signal_range_max);
	    fHistSignal->Fit(fFitSignal,"RNQ","",0.12,signal_range_max);
	    
	    
	    if(iFit==0){
	      fHistSignal->Fit(fFitSignal,"RNQ","",0.12,0.150);
	    }	
	    else{
	      fHistSignal->Fit(fFitSignal,"RNQ","",0.12,0.150);
	    }

	    fFitFinalSignal[i]  = (TF1*)fFitSignal->Clone();
	    
	    Double_t mean   = 0;
	    Double_t e_mean = 0;
	    Double_t signal_sigma   = 0;
	    Double_t e_signal_sigma = 0;
	    Double_t signal_window_min = 0;
	    Double_t signal_window_max = 0;
	    
	    if(iFit==0){
	      mean   = fFitSignal->GetParameter(2);
	      e_mean = fFitSignal->GetParError(2);
	      signal_sigma   = fabs(fFitSignal->GetParameter(3));
	      e_signal_sigma = fabs(fFitSignal->GetParError(3));
	      
	      signal_window_min = mean - signal_sigma*5;
	      signal_window_max = mean + signal_sigma*3;
	      
	    }
	    else if(iFit==1){
	      mean   = fFitSignal->GetParameter(1);
	      e_mean = fFitSignal->GetParError(1);
	      signal_sigma   = fabs(fFitSignal->GetParameter(2));
	      e_signal_sigma = fabs(fFitSignal->GetParError(2));
	      signal_window_min = mean - signal_sigma*5;
	      signal_window_max = mean + signal_sigma*3;
	    }
	  
	    fHistPeakMean->SetBinContent(i+1,mean);
	    fHistPeakMean->SetBinError(i+1,e_mean);
	    fHistPeakWidth->SetBinContent(i+1,signal_sigma);
	    fHistPeakWidth->SetBinError(i+1,e_signal_sigma);
	    
	    Int_t signal_window_bin_min = fHistSignal->GetXaxis()->FindBin(signal_window_min);
	    Int_t signal_window_bin_max = fHistSignal->GetXaxis()->FindBin(signal_window_max);
	    
	    Double_t num_pi0 = 0;
	    Double_t e_num_pi0 = 0;
	    
	    for(Int_t j=signal_window_bin_min; j<signal_window_bin_max; ++j){
	      num_pi0   += fHistSignal->GetBinContent(j);
	      e_num_pi0 += pow(fHistSignal->GetBinError(j),2);
	    }
	    e_num_pi0 = sqrt(e_num_pi0);
	    
	    fHistRawYield->SetBinContent(i+1,num_pi0/bin_width);
	    fHistRawYield->SetBinError(i+1,e_num_pi0/bin_width);
	    
	  }

	  fHistRawYield->SetName(Form("fHistRawYieldPol%d",pol)+fitName[iFit]+tofName[itof]+modName[iMod]);
	  tfout->WriteTObject(fHistRawYield);
	  
	  fHistRawYieldPerEvent = (TH1F*)fHistRawYield->Clone();
	  fHistRawYieldPerEvent->Scale(1./nEvt);
	  fHistRawYieldPerEvent->SetName(Form("fHistRawYieldPerEventPol%d",pol)+fitName[iFit]+tofName[itof]+modName[iMod]);
	  tfout->WriteTObject(fHistRawYieldPerEvent);
	  
	  fHistPeakMean->SetName(Form("fHistPeakMeanPol%d",pol)+fitName[iFit]+tofName[itof]+modName[iMod]);
	  fHistPeakWidth->SetName(Form("fHistPeakWidthPol%d",pol)+fitName[iFit]+tofName[itof]+modName[iMod]);
	  tfout->WriteTObject(fHistPeakMean);
	  tfout->WriteTObject(fHistPeakWidth);
	  
	  c1[itof] = new TCanvas("c1"+fitName[iFit]+tofName[itof],"",1200,1800);
	  c2[itof] = new TCanvas("c2"+fitName[iFit]+tofName[itof],"",1200,1800);
	  
	  c1[itof]->Divide(4,6);
	  for(Int_t i=0; i<binnum; ++i){
	    c1[itof]->cd(i+1);
	    fHistFinalSignal[i]->Draw();
	    fFitFinalSignal[i]->Draw("same");
	  }
	  
	  c1[itof]->SaveAs(Form("cInvariantMassSpectrumPol%d",pol)+fitName[iFit]+tofName[itof]+modName[iMod]+".eps");
	  
	  c2[itof]->Divide(4,6);
	  for(Int_t i=0; i<binnum; ++i){
	    c2[itof]->cd(i+1);
	    fHistFinalSame[i]->Draw();
	    fHistFinalBG[i]->Draw("same");
	  }
	  c2[itof]->SaveAs(Form("cSameScaledMixPol%d",pol)+fitName[iFit]+tofName[itof]+modName[iMod]+".eps");
	  
	  for(Int_t i=0; i<binnum; ++i){
	    
	    fHistFinalSignal[i]->SetName(Form("fHistFinalSignal_No%d_Pol%d",i+1,iPol)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    fFitFinalSignal[i]->SetName(Form("fFitFinalSignal_No%d_Pol%d",i+1,iPol)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    fHistFinalSame[i]->SetName(Form("fHistFinalSame_No%d_Pol%d",i+1,iPol)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    fHistFinalBG[i]->SetName(Form("fHistFinalBG_No%d_Pol%d",i+1,iPol)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    
	    tfout->WriteTObject(fHistFinalSignal[i]);
	    tfout->WriteTObject(fFitFinalSignal[i]);
	    tfout->WriteTObject(fHistFinalSame[i]);
	    tfout->WriteTObject(fHistFinalBG[i]);
	    
	    c3 = new TCanvas("c3","",600,600);
	    c3->cd(1);
	    fHistFinalSame[i]->Draw();
	    fHistFinalBG[i]->Draw("same");
	    c3->SetName(Form("cRawSignalBG_Pol%d_No%d",pol,i+1)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    //tfout->WriteTObject(c3);
	  }
	  for(Int_t i=0; i<binnum; ++i){
	    c3 = new TCanvas("c3","",600,600);
	    c3->cd(1);
	    fHistFinalSignal[i]->Draw("");
	    fFitFinalSignal[i]->Draw("same");
	    c3->SetName(Form("cRawSignal_Pol%d_No%d",pol,i+1)+fitName[iFit]+tofName[itof]+modName[iMod]);
	    //tfout->WriteTObject(c3);
	  }
	  
	}	
      }
    }
  }  
  
}
