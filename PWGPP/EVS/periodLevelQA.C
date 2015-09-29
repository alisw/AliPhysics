#ifndef __CINT__
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "AliDAQ.h"
#endif
#include "map"
using namespace std;
// TODO read number of bits from AliVEvent?
#define NBITS 29
#define NMAXCLASSES 100

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
"kPHI78",
"kEMCEJE",
"kEMCEGA",
"kCentral",
"kSemiCentral",
"kDG5",
"kZED",
"kSPI78",
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
void SetHisto(TH2D* h);
void AddFillSeparationLines(TH1D* h, map<Int_t,Int_t> &fills);

//void periodLevelQA(TString inputFileName ="/afs/cern.ch/work/a/aliqaevs/www/data/2012/LHC12h/pass1/trending.root"){
void periodLevelQA(TString inputFileName ="trending.root"){
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.5);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(1);

  TFile* f = new TFile(inputFileName.Data());
  TTree* t = (TTree*) f->Get("trending");
  TObjArray* classes = new TObjArray();
  TObjString* activeDetectors = new TObjString();
  Int_t run                = 0;
  Int_t fill               = 0;
  Double_t run_duration    = 0;
  Int_t nBCsPerOrbit       = 0;
  Double_t refCounts       = 0;
  Double_t mu              = 0;
  Double_t lumi_seen       = 0;
  Double_t interactionRate = 0;
  ULong64_t class_l0b[NMAXCLASSES]         = {0};
  ULong64_t class_l0a[NMAXCLASSES]         = {0};
  ULong64_t class_l1b[NMAXCLASSES]         = {0};
  ULong64_t class_l1a[NMAXCLASSES]         = {0};
  ULong64_t class_l2b[NMAXCLASSES]         = {0};
  ULong64_t class_l2a[NMAXCLASSES]         = {0};
  Double_t  class_lifetime[NMAXCLASSES]    = {0};
  Double_t  class_lumi[NMAXCLASSES]        = {0};
  ULong64_t alias_recorded[NBITS]          = {0};
  ULong64_t alias_reconstructed[NBITS]     = {0};
  ULong64_t alias_accepted[NBITS]          = {0};
  Double_t alias_l0b_rate[NBITS]           = {0};
  Double_t alias_lifetime[NBITS]           = {0};
  Double_t alias_lumi_recorded[NBITS]      = {0};
  Double_t alias_lumi_reconstructed[NBITS] = {0};
  Double_t alias_lumi_accepted[NBITS]      = {0};
  Int_t timeStart = 0;
  Int_t timeEnd = 0;
  
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("fill",&fill);
  t->SetBranchAddress("bcs",&nBCsPerOrbit);
  t->SetBranchAddress("run_duration",&run_duration);
  t->SetBranchAddress("mu",&mu);
  t->SetBranchAddress("interactionRate",&interactionRate);
  t->SetBranchAddress("refCounts",&refCounts);
  t->SetBranchAddress("lumi_seen",&lumi_seen);
  t->SetBranchAddress("classes",&classes);
  t->SetBranchAddress("class_l0b",&class_l0b);
  t->SetBranchAddress("class_l0a",&class_l0a);
  t->SetBranchAddress("class_l1b",&class_l1b);
  t->SetBranchAddress("class_l1a",&class_l1a);
  t->SetBranchAddress("class_l2b",&class_l2b);
  t->SetBranchAddress("class_l2a",&class_l2a);
  t->SetBranchAddress("class_lifetime",&class_lifetime);
  t->SetBranchAddress("class_lumi",&class_lumi);
  t->SetBranchAddress("alias_recorded",&alias_recorded);
  t->SetBranchAddress("alias_reconstructed",&alias_reconstructed);
  t->SetBranchAddress("alias_accepted",&alias_accepted);
  t->SetBranchAddress("alias_l0b_rate",&alias_lifetime);
  t->SetBranchAddress("alias_lifetime",&alias_lifetime);
  t->SetBranchAddress("alias_lumi_recorded",&alias_lumi_recorded);
  t->SetBranchAddress("alias_lumi_reconstructed",&alias_lumi_reconstructed);
  t->SetBranchAddress("alias_lumi_accepted",&alias_lumi_accepted);
  t->SetBranchAddress("activeDetectors",&activeDetectors);
  t->SetBranchAddress("timeStart",&timeStart);
  t->SetBranchAddress("timeEnd",&timeEnd);

  Int_t nRuns = t->GetEntries();
  TH2D* hClassL0BvsRun      = new TH2D("hClassL0BVsRun","Class L0B vs run",nRuns,0,nRuns,1,0,1);
  TH2D* hClassL2AvsRun      = new TH2D("hClassL2AVsRun","Class L2A vs run",nRuns,0,nRuns,1,0,1);
  TH2D* hClassLifetimeVsRun = new TH2D("hClassLifetimeVsRun","Lifetime class-by-class vs run",nRuns,0,nRuns,1,0,1);
  TH2D* hClassLumiVsRun     = new TH2D("hClassLumiVsRun","Luminosity class-by-class vs run",nRuns,0,nRuns,1,0,1);
  hClassL0BvsRun->SetBit(TH1::kCanRebin);
  hClassL2AvsRun->SetBit(TH1::kCanRebin);
  hClassLifetimeVsRun->SetBit(TH1::kCanRebin);
  hClassLumiVsRun->SetBit(TH1::kCanRebin);
  for (Int_t r=0;r<nRuns;r++){
    t->GetEntry(r);
    printf("run=%i\n",run);
    hClassL0BvsRun->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hClassL2AvsRun->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hClassLifetimeVsRun->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hClassLumiVsRun->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    for (Int_t i=0;i<classes->GetEntriesFast();i++){
      TString className = classes->At(i)->GetName();
      if      (className.Contains("-A-"))      continue;
      else if (className.Contains("-C-"))      continue;
      else if (className.Contains("-E-"))      continue;
      else if (className.Contains("-AC-"))     continue;
      else if (className.Contains("-ACE-"))    continue;
      else if (className.Contains("-GA-"))     continue;
      else if (className.Contains("-GC-"))     continue;
      else if (className.Contains("1A-ABCE-")) continue;
      else if (className.Contains("1C-ABCE-")) continue;
      else if (className.Contains("C0LSR-ABCE-")) continue;
//      printf("%30s %12lli %10lli %10lli %10lli %10lli %10lli\n",className.Data(),L0B[i],L0A[i],L1B[i],L1A[i],L2B[i],L2A[i]);
      hClassL0BvsRun->Fill(Form("%i",run),className.Data(),Double_t(class_l0b[i]));
      hClassL2AvsRun->Fill(Form("%i",run),className.Data(),Double_t(class_l2a[i]));
      hClassLifetimeVsRun->Fill(Form("%i",run),className.Data(),class_lifetime[i]);
      hClassLumiVsRun    ->Fill(Form("%i",run),className.Data(),class_lumi[i]);
    }
  }
  
  hClassL0BvsRun->LabelsDeflate("Y");
  hClassL2AvsRun->LabelsDeflate("Y");
  hClassLifetimeVsRun->LabelsDeflate("Y");
  hClassLumiVsRun->LabelsDeflate("Y");

  map<Int_t,Int_t> fills;
  TH1D* hRecorded[NBITS]          = {0x0};
  TH1D* hReconstructed[NBITS]     = {0x0};
  TH1D* hAccepted[NBITS]          = {0x0};
  TH1D* hRejected[NBITS]          = {0x0};
  TH1D* hAcceptedFraction[NBITS]  = {0x0};
  TH1D* hRejectedFraction[NBITS]  = {0x0};
  TH1D* hLumiRecorded[NBITS]      = {0x0};
  TH1D* hLumiReconstructed[NBITS] = {0x0};
  TH1D* hLumiAccepted[NBITS]      = {0x0};
  

  const Int_t nDetectors=19;

  TH2D* hActiveDetectors = new TH2D("hActiveDetectors","Active detectors",nRuns,0,nRuns,nDetectors,0,nDetectors);
  TH1D* hInteractionRate = new TH1D("hInteractionRate","INEL interaction rate [Hz]",nRuns,0,nRuns);
  TH1D* hMu       = new TH1D("hMu","Average number of INEL collisions per BC",nRuns,0,nRuns);
  TH1D* hBCs      = new TH1D("hBCs","Number of colliding bunches",nRuns,0,nRuns);
  TH1D* hDuration = new TH1D("hDuration","Duration, s",nRuns,0,nRuns);
  TH1D* hLumiSeen = new TH1D("hLumiSeen","Luminosity seen, nb-1",nRuns,0,nRuns);

  TString detName[nDetectors]={"ACORDE","AD","CPV","EMCAL","FMD","HMPID","ITSSPD","ITSSDD","ITSSSD","MUONTRK","MUONTRG","PHOS","PMD","T0","TOF","TPC","TRD","VZERO","ZDC"};
  for (Int_t iDet=0;iDet<nDetectors;iDet++) hActiveDetectors->GetYaxis()->SetBinLabel(iDet+1,detName[iDet].Data());
  
  for (Int_t r=0;r<nRuns;r++){
    t->GetEntry(r);
    hActiveDetectors->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hInteractionRate->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hMu->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hBCs->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hDuration->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    hLumiSeen->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
    for (Int_t iDet=0;iDet<nDetectors;iDet++) {
     hActiveDetectors->Fill(Form("%i",run),detName[iDet].Data(),activeDetectors->String().Contains(detName[iDet].Data()));
    }
    hInteractionRate->SetBinContent(r+1,interactionRate);
    hMu->SetBinContent(r+1,mu);
    hBCs->SetBinContent(r+1,nBCsPerOrbit);
    hDuration->SetBinContent(r+1,run_duration);
    hLumiSeen->SetBinContent(r+1,lumi_seen);
    fills[run]=fill;
  }
  
  
  
  TFile* fclassL0B = new TFile("class_L0B_counts.root","recreate");
  TCanvas* cClassL0B = new TCanvas("cClassL0B","Class L0B vs run",1800,900);
  gPad->SetMargin(0.15,0.01,0.08,0.06);
  SetHisto(hClassL0BvsRun);
  hClassL0BvsRun->Draw("col");
  gPad->Print("class_L0B_counts.pdf(");
  hClassL0BvsRun->Write();
  TCanvas* cL0B = new TCanvas("cL0B","L0B vs run",1800,500);
  cL0B->SetMargin(0.05,0.01,0.18,0.06);
  for (Int_t i=1;i<=hClassL0BvsRun->GetNbinsY();i++){
    TH1D* h = (TH1D*) hClassL0BvsRun->ProjectionX(Form("hClassL0BvsRun_%02i",i),i,i);
    h->SetTitle(Form("%s L0B counts: %.0f",hClassL0BvsRun->GetYaxis()->GetBinLabel(i),h->Integral()));
    SetHisto(h);
    h->Draw();
    AddFillSeparationLines(h,fills);
    gPad->Print("class_L0B_counts.pdf");
    h->Write();
  }
  gPad->Print("class_L0B_counts.pdf]");
  fclassL0B->Close();
  
  TFile* fclassL2A = new TFile("class_L2A_counts.root","recreate");
  TCanvas* cClassL2A = new TCanvas("cClassL2A","Class L2A vs run",1800,900);
  gPad->SetMargin(0.15,0.01,0.08,0.06);
  SetHisto(hClassL2AvsRun);
  hClassL2AvsRun->Draw("col");
  gPad->Print("class_L2A_counts.pdf(");
  hClassL2AvsRun->Write();
  TCanvas* cL2A = new TCanvas("cL2A","L2A vs run",1800,500);
  cL2A->SetMargin(0.05,0.01,0.18,0.06);
  for (Int_t i=1;i<=hClassL2AvsRun->GetNbinsY();i++){
    TH1D* h = (TH1D*) hClassL2AvsRun->ProjectionX(Form("hClassL2AvsRun_%02i",i),i,i);
    h->SetTitle(Form("%s L2A counts: %.0f",hClassL2AvsRun->GetYaxis()->GetBinLabel(i),h->Integral()));
    SetHisto(h);
    h->Draw();
    AddFillSeparationLines(h,fills);
    gPad->Print("class_L2A_counts.pdf");
    h->Write();
  }
  gPad->Print("class_L2A_counts.pdf]");
  fclassL2A->Close();
  
  TFile* fclassLifetime = new TFile("class_lifetime.root","recreate");
  TCanvas* cClassLifetime = new TCanvas("cClassLifetime","Lifetime class-by-class vs run",1800,900);
  gPad->SetMargin(0.15,0.01,0.08,0.06);
  SetHisto(hClassLifetimeVsRun);
  hClassLifetimeVsRun->Draw("col");
  gPad->Print("class_lifetime.pdf(");
  hClassLifetimeVsRun->Write();
  TCanvas* cL2A = new TCanvas("cLifetime","Lifetime vs run",1800,500);
  cL2A->SetMargin(0.05,0.01,0.18,0.06);
  for (Int_t i=1;i<=hClassLifetimeVsRun->GetNbinsY();i++){
    TH1D* h = (TH1D*) hClassLifetimeVsRun->ProjectionX(Form("hClassLifetimeVsRun_%02i",i),i,i);
    h->SetTitle(Form("%s lifetime",hClassLifetimeVsRun->GetYaxis()->GetBinLabel(i)));
    SetHisto(h);
    h->SetFillColor(0);
    h->Draw("");
    AddFillSeparationLines(h,fills);
    gPad->Print("class_lifetime.pdf");
    h->Write();
  }
  gPad->Print("class_lifetime.pdf]");
  fclassLifetime->Close();
  
  TFile* fclassLumi = new TFile("class_lumi.root","recreate");
  TCanvas* cClassLumi = new TCanvas("cClassLumi","Luminosity class-by-class vs run",1800,900);
  gPad->SetMargin(0.15,0.01,0.08,0.06);
  SetHisto(hClassLumiVsRun);
  hClassLumiVsRun->Draw("col");
  gPad->Print("class_lumi.pdf(");
//  TCanvas* clumi = new TCanvas("clumi","lumi vs run",1800,500);
//  clumi->SetMargin(0.05,0.01,0.18,0.06);
//
//  for (Int_t i=1;i<=hClassLumiVsRun->GetNbinsY();i++){
//    TH1D* h = (TH1D*) hClassLumiVsRun->ProjectionX(Form("hClassLumiVsRun_%02i",i),i,i);
//    h->SetTitle(Form("%s luminosity [ub-1]: %.0g",hClassLumiVsRun->GetYaxis()->GetBinLabel(i),h->Integral()));
//    SetHisto(h);
//    h->Draw();
//    AddFillSeparationLines(h,fills);
//    gPad->Print("class_lumi.pdf");
//  }
  gPad->Print("class_lumi.pdf]");
  fclassLumi->Close();
  
  SetHisto(hActiveDetectors);
  SetHisto(hInteractionRate);
  SetHisto(hMu);
  SetHisto(hBCs);
  SetHisto(hDuration);
  SetHisto(hLumiSeen);

  TCanvas* cActiveDetectors = new TCanvas("active_detectors","Active Detectors",1800,500);
  cActiveDetectors->SetMargin(0.05,0.01,0.18,0.06);
  hActiveDetectors->GetYaxis()->SetLabelOffset(0.001);
  hActiveDetectors->SetMaximum(2);
  hActiveDetectors->Draw("col");
  AddFillSeparationLines(hMu,fills);
  gPad->Print("global_properties.pdf(");

  TCanvas* cInteractionRate = new TCanvas("cInteractionRate","Interaction Rate",1800,500);
  cInteractionRate->SetMargin(0.05,0.01,0.18,0.06);
  hInteractionRate->SetFillColor(0);
  hInteractionRate->Draw();
  AddFillSeparationLines(hMu,fills);
  gPad->Print("global_properties.pdf");

  TCanvas* cMu = new TCanvas("mu","mu",1800,500);
  cMu->SetMargin(0.05,0.01,0.18,0.06);
  hMu->SetFillColor(0);
  hMu->Draw("h");
  AddFillSeparationLines(hMu,fills);
  gPad->Print("global_properties.pdf");

  TCanvas* cBCs = new TCanvas("bcs","bcs",1800,500);
  cBCs->SetMargin(0.05,0.01,0.18,0.06);
  hBCs->SetFillColor(0);
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
  hLumiSeen->SetTitle(Form("Luminosity seen [1/ub]: total= %.3f",hLumiSeen->Integral()));
  hLumiSeen->Draw("h");
  AddFillSeparationLines(hLumiSeen,fills);
  gPad->Print("global_properties.pdf)");

  TCanvas* dummy = new TCanvas("dummy","dummy",1800,500);
  gPad->SetMargin(0.05,0.01,0.18,0.06);
  gPad->Print("alias_event_statistics.pdf[");
  gPad->Print("alias_lumi_statistics.pdf[");
  gPad->Print("accepted_fraction.pdf[");
  gPad->Print("rejected_fraction.pdf[");

  TCanvas* cCounts           = new TCanvas("c_alias_counts"   ,"c_alias_counts"   ,1800,500);
  TCanvas* cLumi             = new TCanvas("c_lumi"           ,"c_lumi"           ,1800,500);
  TCanvas* cAcceptedFraction = new TCanvas("accepted_fraction","accepted fraction",1800,500);
  TCanvas* cRejectedFraction = new TCanvas("rejected_fraction","rejected fraction",1800,500);
  cCounts          ->SetMargin(0.05,0.01,0.18,0.06);
  cLumi            ->SetMargin(0.05,0.01,0.18,0.06);
  cAcceptedFraction->SetMargin(0.05,0.01,0.18,0.06);
  cRejectedFraction->SetMargin(0.05,0.01,0.18,0.06);
  cRejectedFraction->SetLogy();

  for (Int_t ibit=0;ibit<NBITS;ibit++) {
    const char* bitName = bitNames[ibit];

    hRecorded[ibit]          = new TH1D(Form("hRecorded%02i"         ,ibit),Form("Recorded: %s"          ,bitName),nRuns,0,nRuns);
    hReconstructed[ibit]     = new TH1D(Form("hReconstructed%02i"    ,ibit),Form("Reconstructed: %s"     ,bitName),nRuns,0,nRuns);
    hAccepted[ibit]          = new TH1D(Form("hAccepted%02i"         ,ibit),Form("Accepted: %s"          ,bitName),nRuns,0,nRuns);
    hRejected[ibit]          = new TH1D(Form("hRejected%02i"         ,ibit),Form("Rejected: %s"          ,bitName),nRuns,0,nRuns);
    hLumiRecorded[ibit]      = new TH1D(Form("hLumiRecorded%02i"     ,ibit),Form("Lumi recorded: %s"     ,bitName),nRuns,0,nRuns);
    hLumiReconstructed[ibit] = new TH1D(Form("hLumiReconstructed%02i",ibit),Form("Lumi reconstructed: %s",bitName),nRuns,0,nRuns);
    hLumiAccepted[ibit]      = new TH1D(Form("hLumiAccepted%02i"     ,ibit),Form("Lumi accepted: %s"     ,bitName),nRuns,0,nRuns);
    
    for (Int_t r=0;r<nRuns;r++){
      t->GetEntry(r);
      hRecorded[ibit]         ->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
      hReconstructed[ibit]    ->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
      hAccepted[ibit]         ->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
      hRejected[ibit]         ->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
      hLumiRecorded[ibit]     ->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
      hLumiReconstructed[ibit]->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
      hLumiAccepted[ibit]     ->GetXaxis()->SetBinLabel(r+1,Form("%i",run));
//      printf("recorded = %i",alias_recorded[ibit]);
//      printf("reconstructed = %i",alias_reconstructed[ibit]);
      hRecorded[ibit]         ->SetBinContent(r+1,alias_recorded[ibit]);
      hReconstructed[ibit]    ->SetBinContent(r+1,alias_reconstructed[ibit]);
      hAccepted[ibit]         ->SetBinContent(r+1,alias_accepted[ibit]);
      hRejected[ibit]         ->SetBinContent(r+1,alias_reconstructed[ibit]-alias_accepted[ibit]);
      hLumiRecorded[ibit]     ->SetBinContent(r+1,alias_lumi_recorded[ibit]);
      hLumiReconstructed[ibit]->SetBinContent(r+1,alias_lumi_reconstructed[ibit]);
      hLumiAccepted[ibit]     ->SetBinContent(r+1,alias_lumi_accepted[ibit]);
    }
    if (hRecorded[ibit]->Integral()<1) continue;
    printf("bit=%i %s\n",ibit,bitName);

    SetHisto(hRecorded[ibit]);
    SetHisto(hReconstructed[ibit]);
    SetHisto(hAccepted[ibit]);
    SetHisto(hRejected[ibit]);
    SetHisto(hLumiRecorded[ibit]);
    SetHisto(hLumiReconstructed[ibit]);
    SetHisto(hLumiAccepted[ibit]);
    hRecorded[ibit]->Sumw2();
    hReconstructed[ibit]->Sumw2();
    hAccepted[ibit]->Sumw2();
    hRejected[ibit]->Sumw2();
    hRecorded[ibit]->SetLineColor(kRed+2);
    hRecorded[ibit]->SetFillColor(kRed+2);
    hReconstructed[ibit]->SetLineColor(kBlue);
    hReconstructed[ibit]->SetFillColor(kBlue);
    hAccepted[ibit]->SetLineColor(kGreen+2);
    hAccepted[ibit]->SetFillColor(kGreen+2);
    hLumiRecorded[ibit]->SetLineColor(kRed+2);
    hLumiRecorded[ibit]->SetFillColor(kRed+2);
    hLumiReconstructed[ibit]->SetLineColor(kBlue);
    hLumiReconstructed[ibit]->SetFillColor(kBlue);
    hLumiAccepted[ibit]->SetLineColor(kGreen+2);
    hLumiAccepted[ibit]->SetFillColor(kGreen+2);
    
    cCounts->cd();
    hRecorded[ibit]->SetTitle(Form("%s trigger counts: recorded=%0.f, total=%.0f, accepted=%.0f",bitName,hRecorded[ibit]->Integral(),hReconstructed[ibit]->Integral(),hAccepted[ibit]->Integral()));
    hRecorded[ibit]->Draw("h");
    hReconstructed[ibit]->Draw("h same");
    hAccepted[ibit]->Draw("h same");
    AddFillSeparationLines(hAccepted[ibit],fills);
    gPad->RedrawAxis();
    gPad->Print("alias_event_statistics.pdf");

    cLumi->cd();
    hLumiRecorded[ibit]->SetTitle(Form("%s luminosity [ub-1]: recorded=%.0g, total=%.0g, accepted=%.0g",bitName,hLumiRecorded[ibit]->Integral(),hLumiReconstructed[ibit]->Integral(),hLumiAccepted[ibit]->Integral()));
    hLumiRecorded[ibit]->Draw("h");
    hLumiReconstructed[ibit]->Draw("h same");
    hLumiAccepted[ibit]->Draw("h same");
    AddFillSeparationLines(hLumiAccepted[ibit],fills);
    gPad->RedrawAxis();
    gPad->Print("alias_lumi_statistics.pdf");

    if (hReconstructed[ibit]->Integral()<1) continue;
    hAcceptedFraction[ibit] = (TH1D*) hReconstructed[ibit]->Clone(Form("hAcceptedFraction%02i",ibit));
    hAcceptedFraction[ibit]->SetTitle(Form("Accepted fraction: %s",bitName));
    hAcceptedFraction[ibit]->Divide(hAccepted[ibit],hReconstructed[ibit],1,1,"B");
    hAcceptedFraction[ibit]->SetFillColor(0);
    hAcceptedFraction[ibit]->SetLineWidth(2);
    hRejectedFraction[ibit] = (TH1D*) hReconstructed[ibit]->Clone(Form("hRejectedFraction%02i",ibit));
    hRejectedFraction[ibit]->SetTitle(Form("Rejected fraction: %s",bitName));
    hRejectedFraction[ibit]->Divide(hRejected[ibit],hReconstructed[ibit],1,1,"B");
    hRejectedFraction[ibit]->SetFillColor(0);
    hRejectedFraction[ibit]->SetLineWidth(2);

    cAcceptedFraction->cd();
    hAcceptedFraction[ibit]->SetTitle(Form("%s: average accepted fraction = %.3f",bitName,hAccepted[ibit]->Integral()/hReconstructed[ibit]->Integral()));
    hAcceptedFraction[ibit]->Draw();
    AddFillSeparationLines(hAcceptedFraction[ibit],fills);
    gPad->Print("accepted_fraction.pdf");
    
    cRejectedFraction->cd();
    hRejectedFraction[ibit]->SetMaximum(1);
    hRejectedFraction[ibit]->SetTitle(Form("%s: average rejected fraction = %.3f",bitName,hRejected[ibit]->Integral()/hReconstructed[ibit]->Integral()));
    hRejectedFraction[ibit]->Draw();
    AddFillSeparationLines(hRejectedFraction[ibit],fills);
    gPad->Print("rejected_fraction.pdf");

  }
  dummy->Print("alias_event_statistics.pdf]");
  dummy->Print("alias_lumi_statistics.pdf]");
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
    if (hRecorded[ibit])      hRecorded[ibit]->Write();
    if (hReconstructed[ibit]) hReconstructed[ibit]->Write();
    if (hAccepted[ibit])      hAccepted[ibit]->Write();
  }
  fstat->Close();

  TFile* faccepted = new TFile("accepted_fraction.root","recreate");
  for (Int_t ibit=0;ibit<NBITS;ibit++) {
    if (hAcceptedFraction[ibit]) hAcceptedFraction[ibit]->Write();
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
  h->SetLineColor(kBlue);
  h->SetFillColor(kBlue);
  h->LabelsOption("av");
  h->SetLineWidth(2);
  h->SetMinimum(0);
}

void SetHisto(TH2D* h){
  h->SetTitleFont(43);
  h->SetTitleSize(15);
  h->GetYaxis()->SetTitleFont(43);
  h->GetXaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelFont(43);
  h->GetZaxis()->SetLabelFont(43);
  h->GetYaxis()->SetTitleSize(25);
  h->GetXaxis()->SetLabelSize(15);
  h->GetYaxis()->SetLabelSize(15);
  h->GetZaxis()->SetLabelSize(15);
  h->GetYaxis()->SetTickLength(0.01);
  h->GetYaxis()->SetDecimals(1);
  h->GetXaxis()->LabelsOption("av");
  h->GetYaxis()->LabelsOption("a");
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
