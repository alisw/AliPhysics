#include <TH1.h>
#include <TMath.h>
#include <TLegend.h>


namespace RawProduction {
  class Output;
}

const int nFrom=4;
const int nTo=4;
float fromRanges[nFrom] = {0.04, 0.05, 0.07, 0.10};
float toRanges[nTo] = {0.20, 0.25, 0.30, 0.40};
RawProduction::Output* output[nFrom][nTo] = {0};

const Double_t maxyPlot = 5., minyPlot = 1.e-3;
const Double_t minReasonableY=1.e-7;

TH1* GetHist(RawProduction::Output* out, const char* trigger, const char* pid, int cent, const char* methode)
{
  //TH1* hist[4][4] = {0};
  char name[256];
  sprintf(name, "%s/c%03i/%s/%s", trigger, cent, pid, methode);
  return out->GetHistogram(name);
}

void SetVariation(TH1* hist, const char* trigger, const char* pid, int cent, const char* methode)
{
  for(int ptb=1; ptb <= hist->GetNbinsX(); ++ptb) {
    int rejected =0;
    Double_t y[nFrom][nTo] = {0};
    for(int fidx=0; fidx<nFrom; fidx++) {
      for(int tidx=0; tidx<nTo; tidx++) {
	TH1* fthist = GetHist(output[fidx][tidx], trigger, pid, cent, methode);
	if( 0. == fthist->GetBinContent(ptb) )
	  rejected++;
	y[fidx][tidx] = fthist->GetBinContent(ptb);
      }
    }
    if( 0 == rejected ) {
      Double_t N = (nFrom*nTo);
      Double_t m = TMath::Mean(N, (Double_t*)y);
      Double_t s = TMath::RMS(N, (Double_t*)y)*TMath::Sqrt(N/(N-1));
      Double_t s_rel = s/m;
      Double_t e_s_rel = s_rel / TMath::Sqrt(2*(N-1));
      hist->SetBinContent(ptb, s_rel);
      hist->SetBinError(ptb, e_s_rel);
      //Printf("%f %f", s_rel, e_s_rel);
    }
    else {
      hist->SetBinContent(ptb, 0.);
      hist->SetBinError(ptb, 0.);
    }
  }
  if( hist->GetMaximum() > maxyPlot )
    Printf("Warning, maximum of %s is %e, larger range max: %e", hist->GetName(), hist->GetMaximum(), maxyPlot);
  if( hist->GetMinimum(minReasonableY) < minyPlot)
    Printf("Warning, minimum of %s is %e, larger range min: %e", hist->GetName(), hist->GetMinimum(minReasonableY), minyPlot);
}

void DrawVar(const char* trigger, const char* pid, int cent)
{
  gStyle->SetOptStat(0);
  TString name = Form("%s_%s_c%03i", trigger, pid, cent);
  TCanvas* canv = new TCanvas(name.Data(), name.Data());
  
  TLegend* leg = new TLegend(0.87, 0.12, 0.99,0.32);
  
  // yr1
  TH1* hist = (TH1*)GetHist(output[0][0], trigger, pid, cent, "yr1")->Clone(Form("%s_%s_c%03i_yr1", trigger, pid, cent));
  hist->SetTitle(Form("Variation of peak position with fit Range, %s, %s, %s", trigger, pid, RawProduction::GetCentString(cent)));
  SetVariation(hist, trigger, pid, cent, "yr1");
  hist->GetXaxis()->SetTitle("p_T [GeV/c]");
  hist->GetYaxis()->SetTitle("#frac{s}{#bar #hat #mu}");
  canv->SetLogy();
  hist->GetYaxis()->SetRangeUser(minyPlot, maxyPlot);
  hist->SetMarkerStyle(21);
  leg->AddEntry(hist, "yr1", "lep");
  hist->DrawCopy("E");
    
  // yr2
  hist = (TH1*)GetHist(output[0][0], trigger, pid, cent, "yr2")->Clone(Form("%s_%s_c%03i_yr2", trigger, pid, cent));
  SetVariation(hist, trigger, pid, cent, "yr2");
  hist->SetMarkerStyle(23);
  hist->SetMarkerColor(kRed);
  hist->SetLineColor(kRed);
  leg->AddEntry(hist, "yr2", "lep");
  hist->Draw("same");
  
  // yr1int
  hist = (TH1*)GetHist(output[0][0], trigger, pid, cent, "yr1int")->Clone(Form("%s_%s_c%03i_yr1int", trigger, pid, cent));
  SetVariation(hist, trigger, pid, cent, "yr1int");
  hist->SetMarkerStyle(24);
  hist->SetMarkerColor(kGreen);
  hist->SetLineColor(kGreen);
  leg->AddEntry(hist, "yr1int", "lep");
  hist->Draw("same");

  // yr2int
  hist = (TH1*)GetHist(output[0][0], trigger, pid, cent, "yr2int")->Clone(Form("%s_%s_c%03i_yr2int", trigger, pid, cent));
  SetVariation(hist, trigger, pid, cent, "yr2int");
  hist->SetMarkerStyle(25);
  hist->SetMarkerColor(kBlue);
  hist->SetLineColor(kBlue);
  leg->AddEntry(hist, "yr2int", "lep");
  hist->Draw("same");

  leg->Draw();

  canv->SaveAs(Form("imgs/RV_c%03i_%s_%s.png", cent, trigger, pid));
  canv->SaveAs(Form("imgs/RV_c%03i_%s_%s.pdf", cent, trigger, pid));
}


void DrawRangeVar()
{
  for(int fidx=0; fidx<nFrom; fidx++) {
    for(int tidx=0; tidx<nTo; tidx++) {
      output[fidx][tidx] = new RawProduction::Output(Form("RawProduction_%.2f_%.2f.root", fromRanges[fidx], toRanges[tidx]));
    }
  }
  
  //TStringToken triggers("kMB kCentral kSemiCentral kPHOSPb", " ");
  TStringToken triggers("kPHOSPb kCentral kSemiCentral kMB", " ");
  while(triggers.NextToken()) {
    //TStringToken pids("All Allcore Allwou Disp Disp2 Dispcore Disp2core Dispwou CPV CPVcore CPV2 CPV2core Both Bothcore Both2 Both2core", " ");
    TStringToken pids("All Allcore Disp CPV Both", " ");
    while(pids.NextToken()) {
      for(int cent = -11; cent < 0; ++cent) {
	if(triggers.EqualTo("kMB") || triggers.EqualTo("kPHOSPb")) {
	  if( -1 == cent || -11 == cent || -10 == cent || -6 == cent )
	    DrawVar(triggers.Data(), pids.Data(), cent);
	}
	if(triggers.EqualTo("kCentral") ) {
	  if( -1 == cent )
	    DrawVar(triggers.Data(), pids.Data(), cent);
	}
	if(triggers.EqualTo("kSemiCentral") ) {
	  if( -11 == cent )
	    DrawVar(triggers.Data(), pids.Data(), cent);
	}
      } // cent
    } // pid
  } // triggers
}


void DrawRangeVariation()
{
  gROOT->LoadMacro("MakeRawProduction.C+g");

  DrawRangeVar();
}
