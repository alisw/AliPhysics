#include "AliHighPtDeDxBase.h"

#include <TMath.h>
#include <TStyle.h>
#include <TROOT.h>

#ifndef __IOSTREAM__
#include <iostream>
#endif

using namespace std;

ClassImp(AliHighPtDeDxBase);

//
// AliHighPtDeDxBase class
//
// This class contains the AliHighPtDeDxBase information 
//

//_________________________________________________________
AliHighPtDeDxBase::AliHighPtDeDxBase():
  TNamed(),
  fIsMc(kFALSE),
  fUseRunCut(kFALSE),
  fRun(0),
  fUseEtaCut(kFALSE),
  fUseEtaCutAbs(kFALSE),
  fEtaLow(-0.8),
  fEtaHigh(0.8),
  fUseFilterCut(kFALSE),
  fFilter(1),
  fUsePhiCut(kFALSE),
  fPhiCutLow(0x0),
  fPhiCutHigh(0x0),
  fEventVtxStatus(-999),
  fEventVtxStatusMc(-999),
  fEventRun(0),
  fEventMag(0),
  fEventTrigger(-1),
  fTrackCharge(0),
  fTrackEta(0),
  fTrackP(0),
  fTrackPt(0),
  fTrackFilter(-999),
  fTrackPhi(-999),
  fTrackPhiCutVariable(-999),
  fTrackDeDx(-999),
  fTrackNcl(-999),
  fTrackPidMc(0),
  fTrackPrimaryMc(-1),
  hVtxStatus(0x0),
  hNevents(0x0),
  hPhi(0x0),
  hEta(0x0),
  hPt(0x0),
  hMeanPt(0x0),
  hNclVsPhiVsPtBefore(0x0),
  hNclVsPhiVsPtAfter(0x0)
{
  // default constructor - do not use
}

//_________________________________________________________
AliHighPtDeDxBase::AliHighPtDeDxBase(const char* name, const char* title):
  TNamed(name, title),
  fIsMc(kFALSE),
  fUseRunCut(kFALSE),
  fRun(0),
  fUseEtaCut(kFALSE),
  fUseEtaCutAbs(kFALSE),
  fEtaLow(-0.8),
  fEtaHigh(0.8),
  fUseFilterCut(kFALSE),
  fFilter(1),
  fUsePhiCut(kFALSE),
  fPhiCutLow(0x0),
  fPhiCutHigh(0x0),
  fEventVtxStatus(-999),
  fEventVtxStatusMc(-999),
  fEventRun(0),
  fEventMag(0),
  fEventTrigger(-1),
  fTrackCharge(0),
  fTrackEta(0),
  fTrackP(0),
  fTrackPt(0),
  fTrackFilter(-999),
  fTrackPhi(-999),
  fTrackPhiCutVariable(-999),
  fTrackDeDx(-999),
  fTrackNcl(-999),
  fTrackPidMc(0),
  fTrackPrimaryMc(-1),
  hVtxStatus(0x0),
  hNevents(0x0),
  hPhi(0x0),
  hEta(0x0),
  hPt(0x0),
  hMeanPt(0x0),
  hNclVsPhiVsPtBefore(0x0),
  hNclVsPhiVsPtAfter(0x0)
{
  // named constructor
}

//_________________________________________________________
AliHighPtDeDxBase::~AliHighPtDeDxBase()
{
  delete hVtxStatus;
  delete hNevents;
  delete hEta;
  delete hPhi;
  delete hPt;
  delete hMeanPt;
  delete hNclVsPhiVsPtBefore;
  delete hNclVsPhiVsPtAfter;
}

//_________________________________________________________
void AliHighPtDeDxBase::Init(Int_t nPtBins, Double_t* ptBins)
{
  //
  // Create histograms
  //
  hVtxStatus = new TH1D("hVtxStatus", "Vtx status - No Vtx = -1, Vtx outside cut = 0, Vtx inside = 1",
			3, -1.5, 1.5);
  hVtxStatus->Sumw2();
  hVtxStatus->SetDirectory(0);

  hNevents = new TH1D("hNevents", "N events - No Vtx = 0, Vtx OK = 1",
		      2, 0, 2);
  hNevents->Sumw2();
  hNevents->SetDirectory(0);
  
  hEta = new TH1D("hEta", "#eta distribution; #eta; Counts",
		  100, -1.0, 1.0);
  hEta->Sumw2();
  hEta->SetDirectory(0);
  
  hPhi = new TH1D("hPhi", "#varphi distribution; #varphi; Counts",
		  360, 0, TMath::TwoPi());
  hPhi->Sumw2();
  hPhi->SetDirectory(0);
  
  hPt = new TH1D("hPt", "p_{T} spectrum; p_{T} [GeV/c]; Counts",
		 nPtBins, ptBins);
  hPt->Sumw2();
  hPt->SetDirectory(0);
  
  hMeanPt = new TProfile("hMeanPt", "mean p_{T}; p_{T} [GeV/c]; mean p_{T}",
			 nPtBins, ptBins);
  hMeanPt->SetDirectory(0);

  const Int_t nPhiBins = 50;
  const Double_t phiBinSize = TMath::Pi()/9.0/nPhiBins;
  Double_t phiBins[nPhiBins+1];

  for(Int_t i = 0; i <= nPhiBins; i++) {

    phiBins[i] = phiBinSize*i;
  }

  const Int_t nNclBins = 45;
  const Double_t nclBinSize = 90.0/nNclBins;
  Double_t nclBins[nNclBins+1];

  for(Int_t i = 0; i <= nNclBins; i++) {

    nclBins[i] = 69.5 + nclBinSize*i;
  }
  
  hNclVsPhiVsPtBefore = new TH3F("hNclVsPhiVsPtBefore", "<Ncl> vs Pt and #phi (before #phi cut); p_{T} [GeV/c]; fmod(#phi, #pi/9.0)",
				 nPtBins, ptBins, nPhiBins, phiBins, nNclBins, nclBins); 
  hNclVsPhiVsPtBefore->SetDirectory(0);
  hNclVsPhiVsPtAfter = new TH3F("hNclVsPhiVsPtAfter", "<Ncl> vs Pt and #phi (after #phi cut); p_{T} [GeV/c]; fmod(#phi, #pi/9.0)",
				nPtBins, ptBins, nPhiBins, phiBins, nNclBins, nclBins); 
  hNclVsPhiVsPtAfter->SetDirectory(0);
}

//_________________________________________________________
Bool_t AliHighPtDeDxBase::EventAccepted()
{
  if(fUseRunCut && fRun!=fEventRun)
    return kFALSE;

  return kTRUE;
}

//_________________________________________________________
Bool_t AliHighPtDeDxBase::TrackAccepted()
{
  if(fUseFilterCut && !(fTrackFilter&fFilter))
    return kFALSE;

  if(fUseEtaCut && (fTrackEta<fEtaLow || fTrackEta>fEtaHigh))
    return kFALSE;
  
  if(fUseEtaCutAbs && (TMath::Abs(fTrackEta)<fEtaLow || TMath::Abs(fTrackEta)>fEtaHigh))
    return kFALSE;
  
  if(fUsePhiCut) {

    fTrackPhiCutVariable = fTrackPhi;
    if(fEventMag<0)    // for negatve polarity field
      fTrackPhiCutVariable = TMath::TwoPi() - fTrackPhiCutVariable;
    if(fTrackCharge<0) // for negatve charge
      fTrackPhiCutVariable = TMath::TwoPi()-fTrackPhiCutVariable;
    if(fTrackPhiCutVariable < 0)
      cout << "Warning!!!!! phi < 0: " << fTrackPhiCutVariable << endl;
    
    fTrackPhiCutVariable += TMath::Pi()/18.0; // to center gap in the middle
    fTrackPhiCutVariable = fmod(fTrackPhiCutVariable, TMath::Pi()/9.0);

    hNclVsPhiVsPtBefore->Fill(fTrackPt, fTrackPhiCutVariable, fTrackNcl);

    if(fTrackPt>2.0 && fTrackPhiCutVariable<fPhiCutHigh->Eval(fTrackPt) 
       && fTrackPhiCutVariable>fPhiCutLow->Eval(fTrackPt))
      return kFALSE; // reject track

    hNclVsPhiVsPtAfter->Fill(fTrackPt, fTrackPhiCutVariable, fTrackNcl);
  }
  
  return kTRUE;
}

//_________________________________________________________
void AliHighPtDeDxBase::FillEventInfo()
{
  if(fEventTrigger==1) {
    hVtxStatus->Fill(fEventVtxStatus);
    if(fEventVtxStatus != 0) {
      hNevents->Fill(1.0+0.5*fEventVtxStatus);
    }
  }
}

//_________________________________________________________
void AliHighPtDeDxBase::FillTrackInfo(Float_t weight)
{
  hEta->Fill(fTrackEta, weight);
  hPhi->Fill(fTrackPhi, weight);
  hPt->Fill(fTrackPt, weight);
  hMeanPt->Fill(fTrackPt, fTrackPt);
}

//_________________________________________________________
void AliHighPtDeDxBase::Print(Option_t* option) const
{
  cout << ClassName() << " : " << GetName() << endl  
       << "Event cuts: " << endl;
  if(fUseRunCut)
    cout << "   Run cut: " << fRun << endl;
  else
    cout << "   Run cut is diabled " << endl;
  
  cout << "Track cuts: " << endl;
  if(fUseFilterCut)
    cout << "   Filter cut: " << fFilter << endl;
  else
    cout << "   Filter cut is diabled " << endl;
  if(fUseEtaCut)
    cout << "   Eta range: " << fEtaLow << " - " << fEtaHigh << endl;
  if(fUseEtaCutAbs)
    cout << "   |Eta| range: " << fEtaLow << " - " << fEtaHigh << endl;
  else
    cout << "   Eta cut is diabled " << endl;
  if(fUsePhiCut)
    cout << "   Phi cut is ENABLED" << endl;
  else
    cout << "   Phi cut is diabled " << endl;
}

//_________________________________________________________
void AliHighPtDeDxBase::SetUseEtaCut(Bool_t  value)
{
  fUseEtaCut = value;

  if(value == kTRUE)
    fUseEtaCutAbs = kFALSE;
}

//_________________________________________________________
void AliHighPtDeDxBase::SetUseEtaCutAbs(Bool_t  value)
{
  fUseEtaCutAbs = value;

  if(value == kTRUE)
    fUseEtaCut = kFALSE;
}

//_________________________________________________________
TF1* AliHighPtDeDxBase::GetStandardPhiCutLow()
{
  //  TF1* cutLow  = new TF1("StandardPhiCutLow",  "-0.01/x+pi/18.0-0.015", 0, 50);
  TF1* cutLow  = new TF1("StandardPhiCutLow",  "0.1/x/x+pi/18.0-0.025", 0, 50);
  return cutLow;
}

//_________________________________________________________
TF1* AliHighPtDeDxBase::GetStandardPhiCutHigh()
{
  //  TF1* cutHigh = new TF1("StandardPhiCutHigh", "0.55/x/x+pi/18.0+0.03", 0, 50);
  TF1* cutHigh = new TF1("StandardPhiCutHigh", "0.12/x+pi/18.0+0.035", 0, 50);
  return cutHigh;
}


//___________________________________________________________________________
TCanvas* AliHighPtDeDxBase::DrawPhiCutHistograms()
{
  gStyle->SetOptStat(0);

  TCanvas* cPhiCut = FindCanvas("cPhiCut", 1200, 800);
  cPhiCut->SetTitle("phi cut histograms");
  cPhiCut->Clear();
  cPhiCut->Divide(3,2);

  // cPhiCut->cd(1);
  // TH2D* hPhiVsPtBefore = (TH2D*)hNclVsPhiVsPtBefore->Project3D("yx");
  // hPhiVsPtBefore->SetName("hPhiVsPtBefore");
  // hPhiVsPtBefore->SetTitle("Phi vs p_{T} (before cuts)");
  // MakeNice2dHisto(hPhiVsPtBefore, gPad, kTRUE);
  // hPhiVsPtBefore->Draw("COL");
  // fPhiCutHigh->SetRange(2.0, 50.0);
  // fPhiCutHigh->Draw("SAME");
  // fPhiCutLow->SetRange(2.0, 50.0);
  // fPhiCutLow->Draw("SAME");
  
  cPhiCut->cd(1);
  TProfile2D* hNclBefore = hNclVsPhiVsPtBefore->Project3DProfile("yx");
  hNclBefore->SetName("hNclBefore");
  hNclBefore->SetTitle("<Ncl> (before cuts); p_{T} [GeV/c]; #varphi'");
  MakeNice2dHisto(hNclBefore, gPad, kTRUE);
  DrawNice(hNclBefore, 0, 0, 0, "COLZ");
  fPhiCutHigh->Draw("SAME");
  fPhiCutLow->Draw("SAME");
  
  cPhiCut->cd(2);
  gPad->SetGridy();
  TH2D* hNclVsPtBefore = (TH2D*)hNclVsPhiVsPtBefore->Project3D("zx");
  hNclVsPtBefore->SetName("hNclVsPtBefore");
  hNclVsPtBefore->SetTitle("; p_{T} [GeV/c]; Ncl;");
  MakeNice2dHisto(hNclVsPtBefore, gPad, kTRUE);
  hNclVsPtBefore->Draw("COL");
  TProfile* hNclVsPtBeforeProf = hNclVsPtBefore->ProfileX();
  hNclVsPtBeforeProf->SetMarkerStyle(29);
  hNclVsPtBeforeProf->SetMarkerColor(2);
  hNclVsPtBeforeProf->Draw("SAME");

  // cPhiCut->cd(2);
  // TH2D* hPhiVsPtAfter = (TH2D*)hNclVsPhiVsPtAfter->Project3D("yx");
  // MakeNice2dHisto(hNclVsPtAfter, gPad, kTRUE);
  // hPhiVsPtAfter->SetName("hPhiVsPtAfter");
  // hPhiVsPtAfter->SetTitle("Phi vs p_{T} (after cuts)");
  // DrawNice(hPhiVsPtAfter, 0, 0, 0, "COL");
			  
  cPhiCut->cd(4);
  TProfile2D* hNclAfter = (TProfile2D*)hNclVsPhiVsPtAfter->Project3DProfile("yx");
  hNclAfter->SetName("hNclAfter");
  hNclAfter->SetTitle("<Ncl> (after cuts); p_{T} [GeV/c]; #varphi'");
  MakeNice2dHisto(hNclAfter, gPad, kTRUE);
  DrawNice(hNclAfter, 0, 0, 0, "COLZ");
  
  cPhiCut->cd(5);
  gPad->SetGridy();
  TH2D* hNclVsPtAfter = (TH2D*)hNclVsPhiVsPtAfter->Project3D("zx");
  hNclVsPtAfter->SetName("hNclVsPtAfter");
  hNclVsPtAfter->SetTitle("; p_{T} [GeV/c]; Ncl;");
  MakeNice2dHisto(hNclVsPtAfter, gPad, kTRUE);
  hNclVsPtAfter->Draw("COL");
  TProfile* hNclVsPtAfterProf = hNclVsPtAfter->ProfileX();
  hNclVsPtAfterProf->SetMarkerStyle(29);
  hNclVsPtAfterProf->SetMarkerColor(2);
  hNclVsPtAfterProf->Draw("SAME");

  cPhiCut->cd(3);
  gPad->SetGridy();
  TH1D* hEfficiency = (TH1D*)hNclVsPhiVsPtAfter->Project3D("x");
  hEfficiency->SetName("hEfficiency");
  hEfficiency->SetTitle("; p_{T} [GeV/c]; Cut efficiency;");  
  TH1D* hHelp = (TH1D*)hNclVsPhiVsPtBefore->Project3D("x");
  hEfficiency->Divide(hEfficiency, hHelp, 1, 1, "B");
  delete hHelp;
  MakeNice1dHisto(hEfficiency, gPad);
  hEfficiency->SetMarkerStyle(20);
  hEfficiency->Draw("P");

  return cPhiCut;
}

/*
  These methods does not fit in very well here, but for now I put them
  here. Maybe they should go in some tool class.
 */

//_______________________________________________________________________
void AliHighPtDeDxBase::MakeNice1dHisto(TH1* hist, TVirtualPad* c1)
{
  if(c1) {
    
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.02);
    c1->SetBottomMargin(0.12);
    c1->SetTopMargin(0.06);
  }

  gStyle->SetTitleH(0.08);
  gStyle->SetTitleW(0.46);
  gStyle->SetTitleX(0.29);
  gStyle->SetTitleY(0.99);
  hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetTitleOffset(0.7);
  hist->GetYaxis()->SetTitleOffset(0.8);
}

//_______________________________________________________________________
void AliHighPtDeDxBase::MakeNice2dHisto(TH2* hist, TVirtualPad* c1, Bool_t colz)
{
  gStyle->SetTitleH(0.08);
  gStyle->SetTitleW(0.76);
  gStyle->SetTitleX(0.14);
  gStyle->SetTitleY(0.98);
  if(colz)
    gStyle->SetTitleW(0.68);

  hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetTitleOffset(0.7);
  hist->GetYaxis()->SetTitleOffset(0.8);
  if(c1) {
    
    c1->SetLeftMargin(0.12);
    if(!colz)
      c1->SetRightMargin(0.02);
    else
      c1->SetRightMargin(0.12);      
    c1->SetBottomMargin(0.12);
    c1->SetTopMargin(0.06);
  }
}
 
//______________________________________________________________________
TCanvas* AliHighPtDeDxBase::FindCanvas(const Char_t* canvasName, 
				       Int_t xwidth, Int_t ywidth)
{
  
  TCanvas *c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(canvasName);
  if(!c1)
    c1 = new TCanvas(canvasName, canvasName, xwidth, ywidth);
  else
    c1->SetWindowSize(xwidth, ywidth);
  
  return c1;
}

//_______________________________________________________________________
TCanvas* AliHighPtDeDxBase::DrawNice(TH1* hist, const Char_t* canvasName,
				     Int_t xwidth, Int_t ywidth, const Char_t* option)
{
  TCanvas* canv = 0;

  if(canvasName) {
    canv = FindCanvas(canvasName, xwidth, ywidth);
    canv->Clear();
  }

  Bool_t is2dHist = hist->InheritsFrom("TH2");
  
  if(is2dHist) {
    
    TString opt(option);
    opt.ToUpper();
    
    if (opt.Contains("COL"))
      MakeNice2dHisto((TH2*)hist, canv, kTRUE);
    else
      MakeNice2dHisto((TH2*)hist, canv);
    
  } else
    MakeNice1dHisto(hist, canv);
  
  hist->Draw(option);

  return canv;
}
