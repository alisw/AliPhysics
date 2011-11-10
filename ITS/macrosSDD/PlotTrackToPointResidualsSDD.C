#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TString.h>
#endif


// Macro to plot the output of the ITS aligment task for SDD detector
// For the moment just plot quantities summed over modules
// Orgin: F. Prino


void PlotTrackToPointResidualsSDD(){
  TFile *fil=new TFile("AnalysisResults.root");
  TDirectoryFile* df=(TDirectoryFile*)fil->Get("ITSAlignQA");
  TList* l=(TList*)df->Get("clistITSAlignQA");
  TH2F* hSDDResidXvsX=0x0;
  TH2F* hSDDResidX=0x0;
  for(Int_t iMod=240; iMod<500; iMod++){
    TString hname=Form("hSDDResidXvsX%d",iMod);
    TH2F* h=(TH2F*)l->FindObject(hname);
    if(h){
      if(hSDDResidXvsX==0x0){
	hSDDResidXvsX=(TH2F*)h->Clone("hSDDResidXvsXAll");
      }else{
	hSDDResidXvsX->Add(h);
      }
    }
    TString hname2=Form("hSDDResidX%d",iMod);
    TH2F* h2=(TH2F*)l->FindObject(hname2);
    if(h2){
      if(hSDDResidX==0x0){
	hSDDResidX=(TH2F*)h2->Clone("hSDDResidXvsPtAll");
      }else{
	hSDDResidX->Add(h2);
      }
    }
  }

  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  TCanvas* c1=new TCanvas("c1","",700,1000);
  c1->Divide(1,2);
  c1->cd(1);
  hSDDResidXvsX->Draw("colz");
  hSDDResidXvsX->GetYaxis()->SetTitle("Track-to-point residual (cm)");
  hSDDResidXvsX->GetXaxis()->SetTitle("X local (cm)");
  TProfile* hpfx=hSDDResidXvsX->ProfileX();
  hpfx->SetLineWidth(2);
  hpfx->Draw("same");
  c1->cd(2);
  Int_t midbin=hSDDResidXvsX->GetXaxis()->FindBin(0.);
  TH1D* hleft=hSDDResidXvsX->ProjectionY("hleft",1,midbin-2);
  TH1D* hright=hSDDResidXvsX->ProjectionY("hright",midbin+2,hSDDResidXvsX->GetNbinsX());
  hleft->Draw();
  hleft->GetYaxis()->SetTitle("Track-to-point residual (cm)");
  c1->Update();
  TPaveStats *st1=(TPaveStats*)hleft->GetListOfFunctions()->FindObject("stats");
  st1->SetY1NDC(0.71);
  st1->SetY2NDC(0.9);
  hright->SetLineColor(2);
  hright->Draw("sames");
  c1->Update();
  TPaveStats *st2=(TPaveStats*)hright->GetListOfFunctions()->FindObject("stats");
  st2->SetY1NDC(0.51);
  st2->SetY2NDC(0.7);
  st2->SetTextColor(2);
  c1->Modified();
  c1->Update();
  TLatex *t2=new TLatex(0.14,0.7,"Left sides");
  t2->SetTextSize(0.048);
  t2->SetNDC();
  TLatex *t3=new TLatex(0.14,0.5,"Right sides");
  t3->SetNDC();
  t3->SetTextSize(0.048);
  t3->SetTextColor(2);
  t2->Draw();
  t3->Draw();

  TCanvas* c2=new TCanvas("c2","",700,1000);
  c2->Divide(1,2);
  c2->cd(1);
  hSDDResidX->GetXaxis()->SetRangeUser(0.,10.);
  hSDDResidX->Draw("colz");
  hSDDResidX->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hSDDResidX->GetYaxis()->SetTitle("Track-to-point residual (cm)");
  c2->cd(2);
  TGraphErrors* grms=new TGraphErrors(0);
  for(Int_t iBin=1; iBin<=hSDDResidX->GetNbinsX(); iBin++){
    TH1D* htmp=hSDDResidX->ProjectionY("htmp",iBin,iBin);
    if(htmp->GetEntries()>10.){
      Int_t npt=grms->GetN();
      grms->SetPoint(npt,hSDDResidX->GetXaxis()->GetBinCenter(iBin),htmp->GetRMS());
      grms->SetPointError(npt,0.5*hSDDResidX->GetXaxis()->GetBinWidth(iBin),htmp->GetRMSError());
    }
  }
  grms->SetMarkerStyle(20);
  grms->GetXaxis()->SetLimits(0.,10.);
  grms->Draw("AP");
  grms->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  grms->GetYaxis()->SetTitle("RMS of Track-to-point residual (cm)");
}
