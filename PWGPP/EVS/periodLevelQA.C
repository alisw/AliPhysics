#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "map"
using namespace std;
// TODO read number of bits from AliVEvent?
#define NBITS 29
TString bitNames[NBITS] = {
"kMB",
"kINT7",
"kMUON",
"kHighMult",
"kEMC1",
"kCINT5",
"kCMUS5",
"kMUSH7",
"kMUL7",
"kMUU7",
"kEMC7",
"kMUS7",
"kPHI1",
"kPHI7/kPHI8",
"kEMCEJE",
"kEMCEGA",
"kCentral",
"kSemiCentral",
"kDG5",
"kZED",
"kSPI7/kSPI8",
"kINT8",
"kMuonSingleLowPt8",
"kMuonSingleHighPt8",
"kMuonLikeLowPt8",
"kMuonUnlikeLowPt8",
"kMuonUnlikeLowPt0",
"kUserDefined",
"kTRD"
};

void SetHisto(TH1D* h);
void AddFillSeparationLines(TH1D* h, map<Int_t,Int_t> &fills);

void periodLevelQA(TString inputFileName ="trending.root"){
  TFile* f = new TFile(inputFileName.Data());
  TTree* t = (TTree*) f->Get("trending");
  Int_t run;
  Int_t fill = 0;
  Int_t nBCsPerOrbit = 0;
  Double_t duration = 0;
  Double_t mu = 0;
  Double_t lumi_seen = 0;
  UInt_t l0b = 0;
  Int_t all[NBITS];
  Int_t accepted[NBITS];
  
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("fill",&fill);
  t->SetBranchAddress("bcs",&nBCsPerOrbit);
  t->SetBranchAddress("duration",&duration);
  t->SetBranchAddress("mu",&mu);
  t->SetBranchAddress("l0b",&l0b);
  t->SetBranchAddress("all",&all);
  t->SetBranchAddress("accepted",&accepted);
  t->SetBranchAddress("lumi_seen",&lumi_seen);
  Int_t nRuns = t->GetEntries();
  map<Int_t,Int_t> fills;
  TH1D* hAll[NBITS];
  TH1D* hAccepted[NBITS] = {0x0};
  TH1D* hRejected[NBITS] = {0x0};
  TH1D* hAcceptedFraction[NBITS] = {0x0};
  TH1D* hRejectedFraction[NBITS] = {0x0};
  
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.5);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(1);
  TCanvas* dummy = new TCanvas("dummy","dummy",1800,500);
  gPad->SetMargin(0.05,0.01,0.18,0.06);
  gPad->Print("global_properties.pdf[");
  gPad->Print("accepted_event_statistics.pdf[");
  gPad->Print("accepted_fraction.pdf[");
  gPad->Print("rejected_fraction.pdf[");

  TH1D* hMu       = new TH1D("hMu","Average number of minimum bias collisions per BC",nRuns,0,nRuns);
  TH1D* hBCs      = new TH1D("hBCs","Number of colliding bunches",nRuns,0,nRuns);
  TH1D* hDuration = new TH1D("hDuration","Duration, s",nRuns,0,nRuns);
  TH1D* hLumiSeen = new TH1D("hLumiSeen","Luminosity seen, nb-1",nRuns,0,nRuns);

  for (Int_t r=0;r<nRuns;r++){
    t->GetEntry(r);
    hMu->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hBCs->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hDuration->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hLumiSeen->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hMu->SetBinContent(r+1,mu);
    hBCs->SetBinContent(r+1,nBCsPerOrbit);
    hDuration->SetBinContent(r+1,duration);
    hLumiSeen->SetBinContent(r+1,lumi_seen/1000.);
    fills[run]=fill;
  }
  SetHisto(hMu);
  SetHisto(hBCs);
  SetHisto(hDuration);
  SetHisto(hLumiSeen);
  TCanvas* cMu = new TCanvas("mu","mu",1800,500);
  cMu->SetMargin(0.05,0.01,0.18,0.06);
  hMu->Draw("h");
  AddFillSeparationLines(hMu,fills);
  gPad->Print("global_properties.pdf");

  TCanvas* cBCs = new TCanvas("bcs","bcs",1800,500);
  cBCs->SetMargin(0.05,0.01,0.18,0.06);
  hBCs->Draw("h");
  AddFillSeparationLines(hBCs,fills);
  gPad->Print("global_properties.pdf");
  
  TCanvas* cDuration = new TCanvas("duration","duration",1800,500);
  cDuration->SetMargin(0.05,0.01,0.18,0.06);
  hDuration->SetTitle(Form("Duration in seconds: total= %.0f s = %.0f h",hDuration->Integral(),hDuration->Integral()/3600));
  hDuration->Draw("h");
  AddFillSeparationLines(hDuration,fills);
  gPad->Print("global_properties.pdf");

  TCanvas* cLumiSeen = new TCanvas("lumiseen","lumi seen",1800,500);
  cLumiSeen->SetMargin(0.05,0.01,0.18,0.06);
  hLumiSeen->SetTitle(Form("Luminosity seen [1/ub]: total= %.0f",hLumiSeen->Integral()));
  hLumiSeen->Draw("h");
  AddFillSeparationLines(hLumiSeen,fills);
  gPad->Print("global_properties.pdf");


  dummy->Print("global_properties.pdf]");

  for (Int_t ibit=0;ibit<NBITS;ibit++) {
    printf("bit=%i\n",ibit);
    const char* bitName = bitNames[ibit];
    hAll[ibit]      = new TH1D(Form("hAll%02i"     ,ibit),Form("All: %s"     ,bitName),nRuns,0,nRuns);
    hAccepted[ibit] = new TH1D(Form("hAccepted%02i",ibit),Form("Accepted: %s",bitName),nRuns,0,nRuns);
    hRejected[ibit] = new TH1D(Form("hRejected%02i",ibit),Form("Rejected: %s",bitName),nRuns,0,nRuns);
    for (Int_t r=0;r<nRuns;r++){
      t->GetEntry(r);
      hAll[ibit]     ->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
      hAccepted[ibit]->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
      hRejected[ibit]->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
      hAll[ibit]->SetBinContent(r+1,all[ibit]);
      if (all[ibit]) hAccepted[ibit]->SetBinContent(r+1,accepted[ibit]);
      if (all[ibit]) hRejected[ibit]->SetBinContent(r+1,all[ibit]-accepted[ibit]);
    }
    if (hAll[ibit]->Integral()<1) continue;
    SetHisto(hAll[ibit]);
    SetHisto(hAccepted[ibit]);
    SetHisto(hRejected[ibit]);
    hAll[ibit]->Sumw2();
    hAccepted[ibit]->Sumw2();
    hRejected[ibit]->Sumw2();
    hAll[ibit]->SetLineColor(kBlue);
    hAll[ibit]->SetFillColor(kBlue);
    hAccepted[ibit]->SetLineColor(kGreen+2);
    hAccepted[ibit]->SetFillColor(kGreen+2);
    hAcceptedFraction[ibit] = (TH1D*) hAll[ibit]->Clone(Form("hAcceptedFraction%02i",ibit));
    hAcceptedFraction[ibit]->SetTitle(Form("Accepted fraction: %s",bitName));
    hAcceptedFraction[ibit]->Divide(hAccepted[ibit],hAll[ibit],1,1,"B");
    hAcceptedFraction[ibit]->SetFillColor(0);
    hAcceptedFraction[ibit]->SetLineWidth(2);
    hRejectedFraction[ibit] = (TH1D*) hAll[ibit]->Clone(Form("hRejectedFraction%02i",ibit));
    hRejectedFraction[ibit]->SetTitle(Form("Rejected fraction: %s",bitName));
    hRejectedFraction[ibit]->Divide(hRejected[ibit],hAll[ibit],1,1,"B");
    hRejectedFraction[ibit]->SetFillColor(0);
    hRejectedFraction[ibit]->SetLineWidth(2);
    
    TCanvas* cAll = new TCanvas(Form("all_%s",bitName),Form("All: %s",bitName),1800,500);
    cAll->SetMargin(0.05,0.01,0.18,0.06);
    hAll[ibit]->SetTitle(Form("%s: total=%.0f accepted=%.0f",bitName,hAll[ibit]->Integral(),hAccepted[ibit]->Integral()));
    hAll[ibit]->Draw("h");
    hAccepted[ibit]->Draw("h same");
    AddFillSeparationLines(hAccepted[ibit],fills);
    gPad->RedrawAxis();
    gPad->Print("accepted_event_statistics.pdf");
    
    TCanvas* cAcceptedFraction = new TCanvas(Form("accepted_fraction_%s",bitName),Form("Accepted Fraction: %s",bitName),1800,500);
    cAcceptedFraction->SetMargin(0.05,0.01,0.18,0.06);
    hAcceptedFraction[ibit]->SetTitle(Form("%s: average accepted fraction = %.3f",bitName,hAccepted[ibit]->Integral()/hAll[ibit]->Integral()));
    hAcceptedFraction[ibit]->Draw();
    AddFillSeparationLines(hAcceptedFraction[ibit],fills);
    gPad->Print("accepted_fraction.pdf");
    
    TCanvas* cRejectedFraction = new TCanvas(Form("rejected_fraction_%s",bitName),Form("Rejected Fraction: %s",bitName),1800,500);
    cRejectedFraction->SetMargin(0.05,0.01,0.18,0.06);
    hRejectedFraction[ibit]->SetTitle(Form("%s: average rejected fraction = %.3f",bitName,hRejected[ibit]->Integral()/hAll[ibit]->Integral()));
    hRejectedFraction[ibit]->Draw();
    AddFillSeparationLines(hRejectedFraction[ibit],fills);
    gPad->Print("rejected_fraction.pdf");

  }
  dummy->Print("accepted_event_statistics.pdf]");
  dummy->Print("accepted_fraction.pdf]");
  dummy->Print("rejected_fraction.pdf]");

  TFile* fglobal = new TFile("global_properties.root","recreate");
  hMu->Write();
  hBCs->Write();
  hDuration->Write();
  hLumiSeen->Write();
  fglobal->Close();

  TFile* fstat = new TFile("accepted_event_statistics.root","recreate");
  for (Int_t ibit=0;ibit<NBITS;ibit++) {
    if (hAll[ibit]) hAll[ibit]->Write();
    if (hAccepted[ibit]) hAccepted[ibit]->Write();
  }
  fstat->Close();

  TFile* faccepted = new TFile("accepted_fraction.root","recreate");
  for (Int_t ibit=0;ibit<NBITS;ibit++) {
    if (hRejectedFraction[ibit]) hRejectedFraction[ibit]->Write();
  }
  faccepted->Close();

  TFile* frejected = new TFile("rejected_fraction.root","recreate");
  for (Int_t ibit=0;ibit<NBITS;ibit++) {
    if (hRejectedFraction[ibit]) hRejectedFraction[ibit]->Write();
  }
  frejected->Close();

}

void SetHisto(TH1D* h){
  h->SetTitleFont(43);
  h->SetTitleSize(25);
  h->GetYaxis()->SetTitleFont(43);
  h->GetXaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetTitleSize(25);
  h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetLabelSize(25);
  h->GetYaxis()->SetTickLength(0.01);
  h->GetYaxis()->SetTitleOffset(0.6);
  h->GetYaxis()->SetDecimals(1);
  h->LabelsOption("av");
  h->SetLineWidth(2);
  h->SetLineColor(kBlue);
  h->SetMinimum(0);
}

void AddFillSeparationLines(TH1D* h, map<Int_t,Int_t> &fills){
  gPad->Update();
  Double_t ymin = gPad->GetUymin();
  Double_t ymax = gPad->GetUymax();
  TLine * fillSeparationLine = new TLine();
  fillSeparationLine->SetLineColor(kRed);
  fillSeparationLine->SetLineWidth(1);
  for(Int_t iBin = 1; iBin < h->GetXaxis()->GetNbins(); iBin++) {
    UInt_t runOld = atoi(h->GetXaxis()->GetBinLabel(iBin));
    UInt_t runNew = atoi(h->GetXaxis()->GetBinLabel(iBin + 1));
    if (fills[runOld]==fills[runNew]) continue;
    fillSeparationLine->DrawLine(iBin,ymin,iBin,ymax);
  }
}
