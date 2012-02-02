
#include "TROOT.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TKey.h"
#include "TRandom.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TBox.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TPolyLine.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TObjArray.h"
#include "TGraph2D.h"
#include "TGaxis.h"
#include "TArrow.h"
#include "TASImage.h"
#include "TFitResult.h"
#include "TStopwatch.h"
#include "TProfile2D.h"
#include "TRegexp.h"
#include <iostream>
#include <fstream>
#include <vector>

// Plot utilities class - loaded by rootlogon.C + .rootrc
#include "Utils.h"
PlotUtils utils;
const char* inputFileName = ""; // Assign in macro
const char* cfDef = "cA"; // s, m, cA or cB. Only for phi projections.
double minRidgeDEta = 0.8;
double maxRidgeDEta = 1.799;
const int gNMax = 12; // Global const -- largest n in VnDelta

// Don't do global fit to points outside these pt bin
// indices. Globally scoped for availability in vntvna(). Reassigned
// within individual functions.
double minFitPtBin = 0, maxFitPtBin = 999;

// Switches
//============================================================
bool usePtBinCenter = 1;
bool ispp = 0;

// Functions
//============================================================
void SetProtonProtonFlag(int val) { ispp = val; }
void SetMinRidgeDEta(int val) { minRidgeDEta = val; }
bool IsProtonProton() { return ispp; }
void SetCFDef(char* def) { cfDef = def; }
void SaveGraphs(const char* outputFileName, TString opt = "");
TF1* DoGlobalFit(int n, int k, int includedPtBin1 = 0, int includedPtBin2 = 999,
		 const char* region = "RIDGE", const char* cfDef = "cA", TString opt = "");
void VnDelta(int n, const TH1 &h, double &vnd, double &vnd_err, TString opt="");
void VnDelta(int n, const TF1 &f, double &vnd, double &vnd_err, TString opt="");
void MakeVnDVsQGraphs(int n1, int n2, int k, const char* region, const char* corrtype, TString opt="");
TF1* Harmonic(TH1* h, int n, TString opt);
TF1* HarmonicSum(TH1* h, int n1, int n2, TString opt="");
TH1* Hist(const char* region, const char* type, int i, int j, int k, TString opt="");
TH2* EtaPhiHist(int cent, int ptt, int pta);
TH2* PttPtaHist();
TH2F* Agreement2DHist(int k, int n);
void ijkFromHistName(TH1* h, int& i, int& j, int& k);

TGraphErrors* VnDeltaVsN(int i, int j, int k, int nMax, TString opt);
TGraphErrors* VnDeltaNFVsN(int i, int j, int k, int nMax, TString opt);
TGraphErrors* GlobalvnVsN(int i, int k, int nMax, TString opt);
TGraphErrors* VnDeltaVsPtAssoc(int i, int k, int n, TString opt="");
TGraphErrors* VnDeltaVsPtTrig(int j, int k, int n, TString opt="");
TGraphErrors* AgreementVsPtSum(int k, int n);
TGraphErrors* VnDVsQ(int n, int k, const char* region = "RIDGE", const char* cfDef = "cA", TString opt="");
TGraphErrors* GlobalFitVsPtAssoc(int i, int k, int n, TString opt=""); // v_n(ptt)v_n(pta) vs. assoc pt
TGraphErrors* vnGFvsPt(int n, int k, int ipt1, int ipt2,
		       const char* region = "RIDGE", 
		       const char* corrtype = "cA", TString opt = "");
TF1* GlobalFitVsQ(int n, int k, int ipt1=0, int ipt2=999, const char* region="RIDGE",
		  const char* corrtype="cA", TString opt=""); // The real fit function
//TF1* GlobalFitVsQSubset(int n, int k, int i, TString opt = "");
TGraph* Luzumv1(int cent, const char* pid);
double MomCorrection(int i, int j, int k);

// Global fit function - give to TF1 constructor
double vntvna(double *x, double *par);
double Chi2(TH1* h, TF1* f, TString opt = ""); // "", "ndf" or "prob"
double ReducedChi2(int i, int k, int n);
double SineError(int i, int j, int k); // RMS of <sin n D phi> vs. n about zero
double vnGF(int n, int k, int ptBin, int includedPtBin1 = 0, int includedPtBin2 = 999,
	    const char* region = "RIDGE", const char* cfDef = "cA", TString opt = "");

// Show +/- solutions in chi2 space by varying one fixed parameter (par)
TCanvas* FitStudy(int icent, int n, int par);
TCanvas* DrawQ(int k, int n, TString opt = "");
TCanvas* SingleDraw(int cent, int ptt, int pta, TString opt="");
TCanvas* SingleDrawEta(int cent, int ptt, int pta, TString opt);
TCanvas* SingleDrawPlain(int cent, int ptt, int pta, TString opt);
TCanvas* DrawVnVsPtTrig(int k, int npta, int ptabins[], TString opt = "");
TCanvas* SingleDraw2D(int cent, int ptt, int pta, TString opt="");
TCanvas* MultiDraw(int cent, TString opt);
TCanvas* DrawVnFromGlobalFit(int n, int ipt1=0, int ipt2=999, int ncb = 999, 
			     int centbins[] = 0, TString opt=""); // little v_n
TCanvas* TwoPanelVert(double ydiv, const char* canvName, const char* canvTitle);
TCanvas* DrawChi2vsPtTrig(int npta, int ptabins[], int ncb, int centbins[], TString opt);
TCanvas* DrawAgreement(int icent, TString opt = "");
TCanvas* Drawv1to5(int ncb, int centbins[], int ipt1=0, int ipt2=999, TString opt="");
TCanvas* Drawv2to5(int ncb, int centbins[], int ipt1=0, int ipt2=999, TString opt="");
TCanvas* DrawGlobalvnVsN(int npt, int ptbins[], int ncb, int centbins[], TString opt="");

// functions needed for above routine
void initialize(const char* ifName = inputFileName);
void AddPrelimStampToCurrentPad(float x1, float y1, float x2, float y2, TString opt="");

void SaveCanvases(TObjArray* canvases, const char* fileName);
void SaveCanvasesFromFile(const char* rootFile, const char* targetDir, const char* tag, const char* fileType);
TObjArray* GetObjectsFromFile(TFile& file, TString clname, TString dir="");


// Global variables
//============================================================
TFile *fin=0;
int maxcent;
int centlow[100];  // low and high percentile for each selection
int centhigh[100];
int maxpt;
double ptmean[100];
double ptlow[100];
double pthigh[100];

int colorsPale[] = {kOrange-3, kGreen-6, kAzure-9, kGray+1, kRed-9, kRed-7, kBlue-9, kViolet-9};
int colorsDark[] = {kOrange-8, kGreen+1, kAzure+2, kBlack,  kCyan+1, kRed+0, kBlue+1, kViolet+1};
int centColors[] = {kBlack, kRed, kOrange-3, kGreen+1, kCyan+2, kBlue-1, kViolet, kGreen+1, kOrange-3, kRed+0,  kOrange-8, kRed+1, kBlue+1};
int centColorsPale[] = {kGray, kRed-9, kOrange-4, kGreen-7, kCyan-7, kBlue-9, kViolet-9, kGreen, kOrange, kRed-4,  kOrange-9, kRed-1, kBlue-1};

TLatex ltx;  // NDC
TLatex ltx2; // absolute

// Global container for graphs, TF1s, etc.
TObjArray* gList = new TObjArray();
TObjArray* gFitList = new TObjArray();
TF1* axisFn = 0;

TGraphErrors *gi, *gj, *gk;
char* trigLabel(int i);
char* asscLabel(int j);
char* centLabel(int k);
double MeanPt(int i, int j, int k, TString t_or_a = "t" /* "t" or "a" */, 
	      TString opt="" /* "" or "squared" */);
int CentBin(double cntLo, double cntHi); 
int PtBin(double ptLo, double ptHi);
int CentColor(double cen1, double cen2);
int CentColorPale(double cen1, double cen2);
int PtColor(double p1, double p2);
int PtColorPale(double p1, double p2);
TStopwatch timer;

// The business
//============================================================

void SaveGraphs(const char* outputFileName, TString opt)
{
  const char* mode = "recreate";
  if (opt.Contains("update"))
    mode = "update";

  TFile* outputFile = new TFile(outputFileName, mode);
  gList->Write();
  gFitList->Write();
  outputFile->Close();
  return;
}

void initialize(const char* ifName) {

  // TH1F files with histograms for:
  // NSJET_y_PTTRIG_PTASSOC_CENT (tracks with |delta-eta| < 0.8)
  // RIDGE_y_PTTRIG_PTASSOC_CENT (tracks with |delta-eta| = 0.8 - 1.6)

  if (fin) // don't keep re-initializing
    return;

  fin = new TFile(ifName);
  if (!fin) {
    Error("FourierPlus - initialize()", "%s not found", ifName);
    gSystem->Exit(1);
  }
  
  TList* info = (TList*)fin->Get("FileInfo");
  if (info)
    info->Print();

  gi = (TGraphErrors*) fin->Get("TrigPtBins");
  gj = (TGraphErrors*) fin->Get("AsscPtBins");
  gk = (TGraphErrors*) fin->Get("EvCentBins");

  maxpt = gi->GetN();
  maxcent = gk->GetN();
  gList->Add(gi);   gList->Add(gj);   gList->Add(gk);

  if (gj->GetN() != maxpt)
    Warning("initialize()", 
	    "Trigger and assc. binning are not symmetric: %d vs. %d", 
	    maxpt, gj->GetN() );
  
  for (int k=0; k<maxcent; k++) {
    centlow[k] = gk->GetX()[k] - gk->GetEX()[k];
    centhigh[k] = gk->GetX()[k] + gk->GetEX()[k];
    //    cout << centlow[k] << " " << centhigh[k] << endl;
  }

  for (int i=0; i<maxpt; i++) {
    ptlow[i] = gi->GetX()[i] - gi->GetEX()[i];
    pthigh[i] = gi->GetX()[i] + gi->GetEX()[i];
    //    cout << ptlow[i] << " " << pthigh[i] << endl;
  }

  // OK, leave hard-coded for now....
  ptmean[0] = 0.66; // 0.5-1 GeV
  ptmean[1] = 1.30; // 1-2 GeV
  ptmean[2] = 2.30; // 2-3 GeV
  ptmean[3] = 3.30; // 3-4 GeV
  ptmean[4] = 4.60; // 4-6 GeV
  ptmean[5] = 6.60; // 6-8 GeV
  ptmean[6] = 11.7; // 8-15 GeV

  if (usePtBinCenter)
    for (int i=0; i<maxpt; i++) {
      ptmean[i] = gi->GetX()[i];
    }

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  ltx.SetNDC();
}

TH1* Hist(const char* region, const char* type, int i, int j, int k, TString opt)
{
  // Options: 
  // "sys": return measured pts. w/systematic errors. 
  // "hi_sys": meas. pts + sys, unmodified stat. errors.
  // "lo_sys": meas. pts - sys, unmodified stat. errors.

  TH1* h_in=0, *h_out=0;
  const char* name = Form("PbPb/%s_%s_%d_%d_%d",region,type,i,j,k);
  h_in = (TH1*)fin->Get(name);
  if (!h_in) {
    Error("Hist()","%s not found, exiting.", name);
    gSystem->Exit(-1);
  }

  h_out = (TH1*)h_in->Clone(Form("%s_%s", name, opt.Data()));
  
  // Stat --> sys error bars. Error estimated as difference in Fourier
  // series (0<n<6) between same-event dist. and CF.
  if (opt.Contains("sys")) {
    const char* s = Form("PbPb/%s_s_%d_%d_%d",region,i,j,k);
    TH1* htmp = (TH1*)fin->Get(s);
    TH1* hs = 0;
    if (htmp)
      hs = (TH1*)htmp->Clone();
    else
      Error("Hist()", "No hist %s", s);

    double c = double(hs->GetNbinsX())/hs->Integral();
    hs->Scale(c);

    TF1* fs = HarmonicSum(hs,  1,5);
    TF1* fc = HarmonicSum(h_in,1,5);

    for (int n=1; n<=hs->GetNbinsX(); n++) {
      double x = hs->GetBinCenter(n);
      double err = TMath::Abs(fs->Eval(x) - fc->Eval(x));
      double bc = h_in->GetBinContent(n);
      
      if (opt.Contains("hi_sys"))
	h_out->SetBinContent(n, bc+err);
      else if (opt.Contains("lo_sys"))
	h_out->SetBinContent(n, bc-err);
      else {
	h_out->SetBinError(n, err);
	h_out->SetFillColor(kGray);
      }
    }
    
    delete htmp;
    delete hs;
    delete fs;
    delete fc;
  }

  return h_out;
}

void ijkFromHistName(TH1* h, int& i, int& j, int& k)
{
  TString tok, name(h->GetName());
  Ssiz_t from = 0;
  int counter = 0;
  while (name.Tokenize(tok, from, "_")) {
    if (counter==2) i = tok.Atoi();
    if (counter==3) j = tok.Atoi();
    if (counter==4) k = tok.Atoi();
    counter++;
  }
  if (0)
    cout << i<<j<<k << endl;

  return;
}

double PtBinCenterSum(int ptt, int pta)
{
  if ((ptt<0 || ptt>=maxpt) || (pta>ptt))
    Error("PtBinCenterSum()","Bad arg(s): ptt %d pta %d", ptt, pta);
 
  double pti = gi->GetX()[ptt];
  double ptj = gj->GetX()[pta];
  return pti+ptj;
}

int GlobalIndex(int ptt, int pta)
{
  // Find the global index from the trig and assc. pt bins.
  // It is useful that TMath::Binomial(n,k) returns 0 when n < k.
  if (pta > ptt) return -1;
  return TMath::Binomial(ptt+1, 2) + pta;
}

void ijFromGlobalIndex(int q, int& ptt, int& pta)
{
  // Assign ptt and pta from the global index q. This function is
  // effectively the inverse of GlobalIndex(i,j).
  for (int i=0; i<maxpt; i++) {
    for (int j=0; j<=i; j++) {
      if (GlobalIndex(i,j)==q) {
	ptt = i; pta = j;
	return;
      }
    }
  }
  return;
}

// Needs to be updated! Has hard-coded numbers like 27.9. -AMA 11/15/2011
TCanvas* FitStudy(int icent, int n, int par)
{
  // First get the VnDelta vs q graph...
  TGraphErrors* gVn = VnDVsQ(n, icent);
  
  // Then setup fit function
  TF1 *fn = new TF1(Form("fn_cent%dto%d_n%d", centlow[icent],centhigh[icent], n), 
		    vntvna, 0.0, 27.9, maxpt);
  fn->SetNpx(1000);
  fn->SetLineColor(kRed);

  TCanvas* c = new TCanvas(Form("chi2_cent%dto%d_n%d", centlow[icent],centhigh[icent], n), 
			   Form("chi2 test for %s, varying v_{%d}(%.1f)",
				centLabel(icent), n, ptmean[par]), 
			   700, 1200);
  c->Divide(1, 3, 0.001, 0.001);
  
  TGraphErrors *g = new TGraphErrors();
  TGraphErrors *chi2 = new TGraphErrors();
  chi2->SetTitle(Form("Total #chi^{2} for n=%d;v_{%d}(p_{T}=%.2f);total #chi^{2}", n,n, ptmean[par]));
  utils.set_tgraph_props(g, kGray, kGray, kFullCircle, 1.3);
  utils.set_tgraph_props(chi2, kBlue, kBlue, kFullCircle, 0.8);
  g->SetLineWidth(3);

  double bestvn     = vnGF(n,icent,par,0,999,"RIDGE","cA","");
  double bestvn_err = vnGF(n,icent,par,0,999,"RIDGE","cA","err");

  int nsteps = 500;

  double parmax = TMath::Abs(2*(bestvn + bestvn_err));
  double parmin = -parmax;

  // First panel: step thru one fixed vn
  c->cd(1);
  TH1F* hf1 = gPad->DrawFrame(0.0, 1.2*parmin, 12.5, 1.2*parmax);
  hf1->SetTitle(Form("title;p_{T};v_{%d}(p_{T})", n));

  for (int step=0; step<=nsteps; step++) {
    double stepsize =  (parmax-parmin)/nsteps;
    double val = parmin + step*stepsize;
    
    for (int ip=0;ip<maxpt;ip++) fn->SetParameter(ip,0.0);
    fn->FixParameter(par, val);

    gVn->Fit(fn, "QR");
    for (int i=0;i<maxpt;i++) {
      g->SetPoint(i, ptmean[i], fn->GetParameter(i));
      g->SetPointError(i, 0.0, fn->GetParError(i));
    }
    if (step%50==0)
      g->DrawClone("elpsame");

    // store chi2 at this val
    chi2->SetPoint(step, val, fn->GetChisquare());
  }
  // Finally, draw once with no fixed parameters...
  utils.set_tgraph_props(g, kBlack, kBlack, kFullCircle, 1.3);
  for (int ip=0;ip<maxpt;ip++) fn->ReleaseParameter(ip);
  gVn->Fit(fn,"QR");
  for (int i=0;i<maxpt;i++) {
    g->SetPoint(i, ptmean[i], fn->GetParameter(i));
    g->SetPointError(i, 0.0, fn->GetParError(i));
  }
  g->DrawClone("elpsame");
  ltx.SetTextSize(0.08);
  ltx.DrawLatex(0.6, 0.8, Form("v_{%d}, %s", n, centLabel(icent)));

  // Plot chi2 over full range
  c->cd(2);
  c->Update();
  chi2->DrawClone("alp");

  // and zoomed to see 1 sigma interval
  c->cd(3);
  c->Update();
  double ymin = 1e9;//chi2->GetHistogram()->GetMinimum();
  for (int j=0; j<chi2->GetN(); j++)
    if (chi2->GetY()[j] < ymin)
      ymin = chi2->GetY()[j];
  double nsigma = n==2? 3 : 3;
  TH1F* hf3 = gPad->DrawFrame(parmin, ymin - 1, parmax, ymin + nsigma);
  hf3->SetTitle(Form("Total #chi^{2} for n=%d;v_{%d}(p_{T}=%.2f);total #chi^{2}", n,n, ptmean[par]));
  chi2->DrawClone("lpsame");
  TLine best, onesig;
  best.SetLineWidth(3);
  onesig.SetLineWidth(3);
  onesig.SetLineColor(kGray);
  best.DrawLine(parmin, ymin, parmax, ymin);
  onesig.DrawLine(parmin, ymin+1, parmax, ymin+1);
  TLatex lbl;
  lbl.DrawLatex(parmin + 0.05*(parmax-parmin), ymin+0.1, Form("Best v_{%d}(%.2f)",n,ptmean[par]));
  lbl.DrawLatex(parmin + 0.05*(parmax-parmin), ymin+1.1, "+ 1#sigma");

  return c;
}

TCanvas* MultiDraw(int cent, TString opt)
{
  // TODO: make this function much more flexible. Pass in an array of
  //  bins or histos and m,n ints to plot in an m x n grid.
  TCanvas* c = new TCanvas(Form("cfmulti_%s_%d",opt.Data(),cent),
			   Form("Corr. Functions - %s", centLabel(cent)),
			   1200, 900);
  c->Divide(maxpt, maxpt-2, 0, 0);

  int ipad=1;

  // Combine first 3 rows (i=0,1,2) to save space
  // Then draw the rest as usual
  for (int i=0; i<maxpt; i++) {
    for (int j=0; j<maxpt; j++) {
      if(j<=i) {
	TCanvas* cs = 0;
	if (opt.Contains("2D")) 
	  cs = SingleDraw2D(cent, i, j, "");
	else if (opt.Contains("deta"))
 	  cs = SingleDrawEta(cent, i, j, "");
	else
	  cs = SingleDraw(cent, i, j, "");
	c->cd(ipad);
	cs->DrawClonePad();
	ipad++;
	if (ipad==maxpt) ipad++;
      }
      else if (i > 2) ipad++;
    }
  }

  // Add some labeling in empty pads
  ltx.SetTextSize(0.15);
  c->cd(1*maxpt);
  ltx.DrawLatex(0.1, 0.8, "p_{T}^{trig}: 0.5-1,");
  ltx.DrawLatex(0.1, 0.6, "1-2, 2-3 GeV/c");  
  c->cd(2*maxpt);
  ltx.DrawLatex(0.1, 0.8, "p_{T}^{trig}:");
  ltx.DrawLatex(0.1, 0.6, "3-4 GeV/c");  
  c->cd(3*maxpt);
  ltx.DrawLatex(0.1, 0.8, "p_{T}^{trig}:");
  ltx.DrawLatex(0.1, 0.6, "4-6 GeV/c");  
  c->cd(4*maxpt);
  ltx.DrawLatex(0.1, 0.8, "p_{T}^{trig}:");
  ltx.DrawLatex(0.1, 0.6, "6-8 GeV/c");  

  c->cd(2*maxpt-2);
  ltx.SetTextSize(0.2);
  ltx.DrawLatex(0.1, 0.7, "Pb+Pb");
  ltx.DrawLatex(0.1, 0.4, "LHC10h");
  ltx.DrawLatex(0.1, 0.1, centLabel(cent));  

  return c;
}

TH2* EtaPhiHist(int cent, int ptt, int pta)
{
 initialize();
  TH2* htmp = (TH2*)fin->Get(Form("PbPb/ETAPHI_c_%1d_%1d_%1d",ptt,pta,cent));
  if (!htmp) {
    Error("EtaPhiHist()",
	  "Problem getting pointer to histogram. cent %d ptt %d pta %d", 
	  cent, ptt, pta);
    return 0;
  }
  TH2* h = (TH2*)htmp->Clone(Form("EtaPhi_c_%1d_%1d_%1d",ptt,pta,cent));
  TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis(), *az = h->GetZaxis();
  ax->SetTitle("#Delta#phi");
  ay->SetTitle("#Delta#eta");
  az->SetTitle("C(#Delta#phi, #Delta#eta)");
  ay->SetRangeUser(-maxRidgeDEta, maxRidgeDEta);
  ax->SetTitleOffset(1.3);
  ay->SetTitleOffset(1.3);
  az->SetTitleOffset(1.4);

  return h;
}

TCanvas* SingleDraw2D(int cent, int ptt, int pta, TString opt)
{
  initialize();
  TCanvas* c = 0; // The canvas to be returned
  double labelScale = 1.3;
  double nDiv = 5;
  
  // See if we have already made this once
  const char* name;
  if (opt=="") 
    name = Form("etaphi_cent%dto%d_%d_%d",
		centlow[cent],centhigh[cent],ptt,pta );
  else
    name = Form("etaphi_cent%dto%d_%d_%d_%s",
		centlow[cent],centhigh[cent],ptt,pta, opt.Data() );

  const char* title = Form("dphideta_cent%dto%d_pTtrig%.2gto%.2g_pTassoc%.2gto%.2g_%s",
			   centlow[cent],centhigh[cent],
			   ptlow[ptt], pthigh[ptt],
			   ptlow[pta], pthigh[pta],
			   opt.Data());

  if (gDirectory->FindObject(name)) {
    c = (TCanvas*) gDirectory->FindObject(name);
    return c;
  }
  else
    c = new TCanvas(name, title, 1);

  TH2* h = EtaPhiHist(cent, ptt,pta);
  if (!h) {
    Error("SingleDraw2D()",
	  "Problem getting pointer to histogram. cent %d ptt %d pta %d", 
	  cent, ptt, pta);
    return 0;
  }
  utils.make_nice_axes(c, h, labelScale, nDiv, nDiv);
  TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis(), *az = h->GetZaxis();
  ax->SetTitle("#Delta#phi [rad]");
  ay->SetTitle("#Delta#eta");
  az->SetTitle("C(#Delta#phi, #Delta#eta)");
  ay->SetRangeUser(-maxRidgeDEta, maxRidgeDEta);
  ax->SetTitleOffset(1.3);
  ay->SetTitleOffset(1.3);
  az->SetTitleOffset(1.4);

  double czlo=0, czhi=5.0;
  if (opt.Contains("zoom")) 
    az->SetRangeUser(czlo, czhi);

  if (opt.Contains("colz"))
    h->Draw("colz");
  else
    h->Draw("surf1");

  ltx.SetTextSize(0.05);
  if (opt.Contains("colz")) {
    if (opt.Contains("zoom")) 
      ltx.DrawLatex(0.6, 0.9, Form("C(#Delta#phi) (C(#Delta#phi) < %.3g)", czhi));
    else
      ltx.DrawLatex(0.6, 0.9, Form("C(#Delta#phi) (full range)"));
  }
  
  // Label pt, centrality, energy
  TString ptStr(Form("#splitline"
		     "{%.3g < p_{T}^{t} < %.3g GeV/c}"
		     "{%.3g < p_{T}^{a} < %.3g GeV/c}",
		     ptlow[ptt], pthigh[ptt],
		     ptlow[pta], pthigh[pta]));
  ltx.SetTextSize(0.07);
  ltx.DrawLatex(0.01, 0.86, ptStr.Data());
  if (ispp)
    ltx.DrawLatex(0.65, 0.9, Form("#splitline{pp 2.76 TeV}{}"));
  else
    ltx.DrawLatex(0.65, 0.9, Form("#splitline{Pb-Pb 2.76 TeV}{%s}", centLabel(cent)));

  return c;
}

TCanvas* SingleDrawEta(int cent, int ptt, int pta, TString opt)
{
  TCanvas* c = 0;
  int lineWidth  = opt.Contains("small") ? 1 : 4;
  double mrkSize = opt.Contains("small") ? 0.5 : 1.5;

  const char* name = Form("deta_cent%dto%d_%d_%d", 
			  centlow[cent],centhigh[cent],ptt,pta);

  // See if we have already made this once
  if (gDirectory->FindObject(name)) {
    c = (TCanvas*) gDirectory->FindObject(name);
    Info("SingleDrawEta()", "%s already created", c->GetName());
    return c;
  }
  else
    c = new TCanvas(name, name, 600, 600);

  TH1 *hNS = 0, *hAS = 0;    // divide-then-project
  TH1 *hNS2 = 0, *hAS2 = 0;  // project-then-divide

  bool includeProjThenDiv = 0;

  if (1) { // Recreate CFs using divide-then-project method

    TH2* hc = (TH2*)fin->Get(Form("PbPb/ETAPHI_c_%1d_%1d_%1d",ptt,pta,cent));
    if (!hc) {
      Error("SingleDrawEta()", "Problem getting hc\n");
      gSystem->Exit(-1);
    }

    int bin1 = hc->GetXaxis()->FindBin(-TMath::PiOver2()+0.001);
    int bin2 = hc->GetXaxis()->FindBin(TMath::PiOver2()-0.001);
    int bin3 = hc->GetXaxis()->FindBin(TMath::PiOver2()+0.001);
    int bin4 = hc->GetXaxis()->FindBin(3*TMath::PiOver2()-0.001);
    
    bool smallEtaRange = false;
    if (smallEtaRange) {
      bin1 = hc->GetXaxis()->FindBin(-TMath::PiOver2()+0.001);
      bin2 = hc->GetXaxis()->FindBin(TMath::PiOver2()-0.001);
      bin3 = hc->GetXaxis()->FindBin(TMath::Pi()-0.5);
      bin4 = hc->GetXaxis()->FindBin(TMath::Pi()+0.5);
    }

    const char* ns_name = Form("ns_deta_cent%dto%d_%d_%d", 
			       centlow[cent],centhigh[cent],ptt,pta);
    const char* as_name = Form("as_deta_cent%dto%d_%d_%d", 
			       centlow[cent],centhigh[cent],ptt,pta);
    
    hNS = hc->ProjectionY(ns_name, bin1, bin2, "e");
    hAS = hc->ProjectionY(as_name, bin3, bin4, "e");

    hNS->Scale(1./(bin2-bin1+1));
    hAS->Scale(1./(bin4-bin3+1));

  }

  if (includeProjThenDiv) { // As created by Project.C (project-then-divide)
    const char* ns = Form("PbPb/ETA_NS_c_%1d_%1d_%1d",ptt,pta,cent);
    const char* as = Form("PbPb/ETA_AS_c_%1d_%1d_%1d",ptt,pta,cent);
    hNS2 = (TH1*)fin->Get(ns);
    hAS2 = (TH1*)fin->Get(as);
    if (!hNS2) {
      Error("SingleDrawEta()", "Problem getting %s\n", ns);
      gSystem->Exit(-1);
    }
    if (!hAS2) {
      Error("SingleDrawEta()", "Problem getting %s\n", as);
      gSystem->Exit(-1);
    }
  }

  utils.set_hist_props(hNS, kBlack, kNone, kBlack, kFullDotLarge, mrkSize);
  utils.set_hist_props(hAS, kRed, kNone, kRed, kFullDotLarge, mrkSize);
  utils.set_ylimits(hNS, hAS);
  hNS->SetLineWidth(lineWidth);
  hAS->SetLineWidth(lineWidth);

  // hNS->Draw("hist");
  // hAS->Draw("histsame");
  hNS->DrawClone("ep");
  hAS->DrawClone("epsame");

  if (includeProjThenDiv) {
    utils.set_hist_props(hNS2, kGray+2, kNone, kGray+2, kFullDotLarge, mrkSize);
    utils.set_hist_props(hAS2, kCyan, kNone, kCyan, kFullDotLarge, mrkSize);
    hNS2->SetLineWidth(lineWidth);
    hAS2->SetLineWidth(lineWidth);

    hNS2->DrawClone("epsame");
    hAS2->DrawClone("epsame");
  }

  return c;
}

TCanvas* SingleDrawPlain(int cent, int ptt, int pta, TString opt)
{
  initialize();
  const char* region = "RIDGE";
  const char* yTitle = Form("%.2g < |#Delta#eta| < %.2g", 
			    minRidgeDEta, maxRidgeDEta);
  if (opt.Contains("ALL")) {
    region = "ALL";
    yTitle = Form("%.2g < |#Delta#eta| < %.2g", 0.0, maxRidgeDEta);
  }
  if (opt.Contains("NSJET")) {
    region = "NSJET";
    yTitle = Form("%.2g < |#Delta#eta| < %.2g", 0.0, minRidgeDEta);
  }
  
  const char* name = Form("%s_cent%dto%d_%d_%d_%s", region,
			  centlow[cent],centhigh[cent],ptt,pta, opt.Data());
  const char* title = Form("dphiPlain_%s_cent%dto%d_pTtrig%.2gto%.2g_pTassoc%.2gto%.2g_%s",
			   region, 
			   centlow[cent],centhigh[cent],
			   ptlow[ptt], pthigh[ptt],
			   ptlow[pta], pthigh[pta],
			   opt.Data());
  TCanvas* c = 0;
  if (gDirectory->FindObject(name)) {

    c = (TCanvas*) gDirectory->FindObject(name);
    Info("SingleDrawPlain()", "%s already created", c->GetName());
    return c;
  }
  else
    c = new TCanvas(name, title, 700, 700);
  c->cd();

  TH1* h_sys = Hist(region, cfDef, ptt,pta,cent, "sys");
  TH1* h     = Hist(region, cfDef, ptt,pta,cent);
  h->SetName("h");

  TAxis* ax = h->GetXaxis();
  TAxis* ay = h->GetYaxis();

  h->SetTitle(Form(";#Delta#phi [rad];C(#Delta#phi),  %s", yTitle));
  h->SetLineWidth(2);
  ax->SetNdivisions(208);
  ay->SetNdivisions(208);
  ax->SetLabelSize(0.04);
  ay->SetLabelSize(0.04);
  ax->SetTitleSize(0.05);
  ay->SetTitleSize(0.05);
  ax->SetTitleOffset(1.5);
  ay->SetTitleOffset(-1.9);
  c->SetFrameLineWidth(2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.15);
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);
  ax->CenterTitle();
  ay->CenterTitle();
  ax->SetTicks("+-");
  ay->SetTicks("+-");

  utils.set_ylimits(h, h, 0.05);
  utils.set_hist_props(h, kBlack, kNone, kBlack, kFullDotLarge, 1.4);
  h->Draw();
  h_sys->Draw("e2psame");
  h->Draw("same");
  if (opt.Contains("harmonics"))
    for (int n=1; n<6; n++)
      Harmonic(h, n, "draw");
  if (opt.Contains("sum"))
    HarmonicSum(h,1,5,"draw");
  
  ltx.SetTextSize(0.045);
  double top = 0.92, gap = 0.08;
  double x1 = (TString(region)=="RIDGE" && (ptlow[ptt] > 4 || ispp))? 0.25 : 0.5;
  ltx.DrawLatex(x1, top-0*gap, Form("p_{T}^{t} %s GeV/c", trigLabel(ptt)));  
  ltx.DrawLatex(x1, top-1*gap, Form("p_{T}^{a} %s GeV/c", asscLabel(pta)));  
  if (ispp)
  ltx.DrawLatex(x1, top-2*gap, Form("pp 2.76 TeV"));  
  else
  ltx.DrawLatex(x1, top-2*gap, Form("Pb-Pb 2.76 TeV, %s", centLabel(cent)));  
  if (0)
    ltx.DrawLatex(0.22, 0.92, Form("|#Delta#eta| > %.1f", minRidgeDEta));  

  // ltx.SetTextColor(kBlack);
  // ltx.SetTextSize(0.03);
  // ltx.DrawLatex(0.7, 0.18, "statistical error only");  
 
  return c;
}

TF1* Harmonic(TH1* h, int n, TString opt)
{
  TF1* f = 0;
  double vndelta=0, vndelta_err=0;
  double phiMin = -0.5*TMath::Pi();
  double phiMax = +1.5*TMath::Pi();
  double b0 = h->Integral() / h->GetNbinsX();
  int color = opt.Contains("gray")? kGray : colorsPale[n-1];

  VnDelta(n, *h, vndelta, vndelta_err);

  f = new TF1(Form("%scos%ddphi", h->GetName(), n), 
	      "[0]*(1. + 2.*[1]*cos([2]*x))",
	      phiMin, phiMax);
    f->SetParameters(b0, vndelta, n);
    f->SetLineColor(color);
    f->SetLineWidth(3);
    if (opt.Contains("draw")) 
      f->Draw("same");

    return f;
}

TF1* HarmonicSum(TH1* h, int n1, int n2, TString opt)
{
  if (n1<1)
    Warning("HarmonicSum()", "n1=%d but should be > 0", n1);

  double phiMin = -0.5*TMath::Pi();
  double phiMax = +1.5*TMath::Pi();

  // Function definition string
  TString fdef("[0]*(1.0 + 2.0*(");
  for (int n=n1; n<=n2; n++) {
    fdef.Append(Form(" [%d]*cos(%d*x) ", n, n));
    fdef.Append( n<n2 ? "+" : "))");
  }

  // Normalization
  TF1 *ff = new TF1("ff", fdef.Data(), phiMin,phiMax);
  ff->SetParameter(0, h->Integral()/h->GetNbinsX());

  // Set parameters 
  for (int n=n1; n<=n2; n++) {

    // Sum of vnt x vna from global fit
    if (opt.Contains("global") ) {
      int i=-1, j=-1, k=-1;
      ijkFromHistName(h,i,j,k);
      double vnt = vnGF(n,k,i,0,999,"RIDGE","cA","");
      double vna = vnGF(n,k,j,0,999,"RIDGE","cA","");
      ff->SetParameter(n, vnt*vna); // TODO: consider including errors
      ff->SetLineColor(kRed+1);
    }
    else { // sum of directly-extracted pair moments
      double Vn=0, VnErr=0;
      VnDelta(n, *h, Vn, VnErr);
      ff->SetParameter(n, Vn);
      ff->SetLineColor(kBlack); // ff->SetLineColor(kAzure);
      //      ff->SetLineWidth(4);
      ff->SetLineStyle(kDashed);
      if (n2==10) {
	ff->SetLineColor(kGreen);
      }
    } 
  }

  if (opt.Contains("draw")) 
    ff->Draw("same");
  
  return ff;
}

TCanvas* SingleDraw(int cent, int ptt, int pta, TString opt)
{
  int lineWidth  = opt.Contains("small") ? 1 : 4;
  double mrkSize = opt.Contains("small") ? 0.5 : 1.5;
  TString dEtaLabel = "";

  const char* region = "RIDGE";
  if (opt.Contains("ALL")) {
    region = "ALL";
    dEtaLabel = Form("%.2g < |#Delta#eta| < %.2g", 0.0, maxRidgeDEta);
  }
  if (opt.Contains("NSJET")) {
    region = "NSJET";
    dEtaLabel = Form("%.2g < |#Delta#eta| < %.2g", 0.0, minRidgeDEta);
  }
  if (TString(region)=="RIDGE")
    dEtaLabel = Form("%.2g < |#Delta#eta| < %.2g", minRidgeDEta, maxRidgeDEta);
  
  initialize();

  TCanvas* c1 = 0; // The canvas to be returned

  // See if we have already made this once
  TString name = Form("dphi_etamin%02d_cent%dto%d_%d_%d", (int)(10*minRidgeDEta), 
			  centlow[cent],centhigh[cent],ptt,pta);
  TString title = Form("dphiFourier_%s_cent%dto%d_pTtrig%.2gto%.2g_pTassoc%.2gto%.2g",
			   region, 
			   centlow[cent],centhigh[cent],
			   ptlow[ptt], pthigh[ptt],
			   ptlow[pta], pthigh[pta]);

  if (!opt.IsNull()) {
    name.Append("_"); name.Append(opt);
    title.Append("_"); title.Append(opt);
  }

  if (gDirectory->FindObject(name.Data())) {

    c1 = (TCanvas*) gDirectory->FindObject(name.Data());
    Info("SingleDraw()", "%s already created", c1->GetName());
    return c1;
  }
  else
    c1 = TwoPanelVert(0.3, name.Data(), title.Data());

  c1->cd(1);
  gPad->SetTopMargin(0.1);
  c1->cd();

  TH1* h = Hist(region, cfDef, ptt,pta,cent);
  TH1* h_sys = Hist(region, cfDef, ptt,pta,cent, "sys");
  utils.set_hist_props(h, kBlack, kNone, kBlack, kFullDotLarge, mrkSize);
  h->SetLineWidth(opt.Contains("small") ? 1 : 2);
  h->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h->GetYaxis()->SetTitle("C(#Delta#phi)");

  TF1* ff5 = HarmonicSum(h, 1, 5, "");
  ff5->SetLineWidth(lineWidth);
  TF1* fflowg = HarmonicSum(h, 1, 5, "global");
  fflowg->SetLineWidth(lineWidth);
  TF1* ff10 = HarmonicSum(h, 1, 10, "");
  ff10->SetLineWidth(lineWidth);
  TF1* ff12 = HarmonicSum(h, 1, 12, "");
  ff12->SetLineWidth(lineWidth);

  c1->cd(1);

  int nXticks = 404, nYticks = 210;
  double pad1Scale = 1.6;
  double pad2Scale = 2.3;
  utils.make_nice_axes(c1, h, pad1Scale, nXticks, nYticks, 0.05, 0.01);
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleOffset(1.5);
  utils.set_ylimits(h,h);

  h->Draw("ep");
  h_sys->Draw("e2psame");
  h->Draw("epsame");

  if (opt.Contains("gray"))
    for (int n=1; n<6; n++)
      Harmonic(h, n, "gray, draw");
  else if (opt.Contains("color"))
    for (int n=1; n<6; n++)
      Harmonic(h, n, "draw");
  
  if (opt.Contains("upto10")) {
    ff10->SetLineWidth(lineWidth+2);
    ff10->SetLineColor(kBlack);
    ff10->DrawClone("same");
    ff10->SetLineWidth(lineWidth);
    ff10->SetLineColor(kGreen);
    ff10->SetLineStyle(kDashed);
    ff10->Draw("same");
  }

  //  ff10->Draw("same");
  ff5->Draw("same");
  h->Draw("p,e,same");
  
  if (opt.Contains("global")) {
    fflowg->Draw("same");
  }

  ltx.SetTextSize(0.06);
  double top = 0.82, gap = 0.08;
  double x1 = (TString(region)=="RIDGE" && ptlow[ptt] > 6)? 0.25 : 0.6;
  if (ispp)
    ltx.DrawLatex(0.2, 0.92, Form("pp 2.76 TeV"));  
  else
    ltx.DrawLatex(0.2, 0.92, Form("Pb-Pb 2.76 TeV, %s central", centLabel(cent)));  
  ltx.DrawLatex(x1, top-0*gap, Form("%.3g < p_{T}^{t} < %.3g GeV/c",
				    ptlow[ptt], pthigh[ptt]));  
  ltx.DrawLatex(x1, top-1*gap, Form("%.3g < p_{T}^{a} < %.3g GeV/c", 
				    ptlow[pta], pthigh[pta]));  
  // ltx.DrawLatex(x1, top-0*gap, Form("p_{T}^{trig} %s GeV/c", trigLabel(ptt)));  
  // ltx.DrawLatex(x1, top-1*gap, Form("p_{T}^{assoc} %s GeV/c", asscLabel(pta)));  
  ltx.DrawLatex(x1, top-2*gap, dEtaLabel.Data());  


  c1->cd(2);

  TH1F *hsub = static_cast <TH1F *> (h->Clone());
  TH1F *hsub10 = static_cast <TH1F *> (h->Clone());
  TH1F *hsubg = (TH1F*) (h->Clone());

  utils.set_hist_props(hsub10,  kBlack, kNone, kBlack, kOpenCircle, mrkSize);
  utils.set_hist_props(hsub,  kAzure, kNone, kAzure,   kFullSquare, mrkSize);
  utils.set_hist_props(hsubg, kRed+1, kNone, kRed+1, kFullSquare, mrkSize);

  for (int i=1;i<=hsub->GetNbinsX();i++) {
    double val, err;
    val = h->GetBinContent(i) * h->GetBinWidth(i);
    err = h->GetBinError(i)   * h->GetBinWidth(i);
    double dphi1 = h->GetXaxis()->GetBinLowEdge(i);
    double dphi2 = h->GetXaxis()->GetBinUpEdge(i);
    double fSum5, fSum10, fSumGlobal;
    fSum5      = ff5->Integral(dphi1, dphi2);
    fSum10     = ff10->Integral(dphi1, dphi2);
    fSumGlobal = fflowg->Integral(dphi1, dphi2);
    hsub->SetBinContent(i,   val/fSum5);
    hsub->SetBinError(i,     err/fSum5);
    hsub10->SetBinContent(i, val/fSum10);
    hsub10->SetBinError(i,   err/fSum10);
    hsubg->SetBinContent(i,  val/fSumGlobal);
    hsubg->SetBinError(i,    err/fSumGlobal);
  }
  utils.make_nice_axes(c1, hsub, pad2Scale, nXticks, 5, 0.05, 0.01);
  if (opt.Contains("global"))
    utils.set_ylimits(hsub, hsubg, 0.2);
  else 
    utils.set_ylimits(hsub, hsub10, 0.2);
  
  hsub->SetTitle(Form(";#Delta#phi [rad];ratio"));
  //  hsub->SetTitle(Form(";#Delta#phi [rad];data/#Sigma_{n=1}^{%s}", opt.Contains("global")?"5":"N"));
  hsub->GetXaxis()->SetTitleOffset(1.25);
  hsub->GetYaxis()->SetTitleOffset(0.7);
  hsub->GetXaxis()->SetTicks("+-");   hsub->GetXaxis()->SetTickLength(0.08);
  hsub10->SetTitle(";#Delta#phi [rad];data/#Sigma_{n=1}^{N}");
  hsub10->GetXaxis()->SetTitleOffset(1.25);
  hsub10->GetYaxis()->SetTitleOffset(0.7);
  hsubg->SetTitle(";#Delta#phi [rad];data/calc");
  hsubg->GetXaxis()->SetTitleOffset(1.25);
  hsubg->GetYaxis()->SetTitleOffset(0.7);

  TLegend* leg = new TLegend(0.2, 0.01, 0.38, 0.38);
  leg->SetFillColor(kNone);
  leg->SetBorderSize(1);
  if (opt.Contains("global")) {
    leg->AddEntry(hsub, "#LTcos(n#Delta#phi)#GT", "epl");
    leg->SetBorderSize(0);
  }
  else
    leg->AddEntry(hsub, "1#leqn#leq5", "epl");

  TLine unity;
  hsub->Draw("ep");

  if (opt.Contains("global")) {
    hsubg->Draw("epsame");
    leg->AddEntry(hsubg, "v(p_{Tt})v(p_{Ta})", "epl");
  }
  else if (opt.Contains("upto10")) {
    TH1* hsub10green = (TH1*)hsub10->Clone();
    utils.set_hist_props(hsub10green,  kBlack, kNone, kGreen, kFullCircle, mrkSize);
    hsub10green->Draw("epsame");
    hsub10->Draw("epsame");
    leg->AddEntry(hsub10green, "1#leqn#leq10", "epl");
  }

  unity.DrawLine(-TMath::PiOver2(), 1.0, 3*TMath::PiOver2(), 1.0);

  c1->cd(2);
  if (opt.Contains("global"))
    leg->Draw();

  if (1) {
    // calculate chi^2
    TF1* flat = new TF1("flat", "1.0", -TMath::PiOver2(), 3*TMath::PiOver2());
    double chi2 = Chi2(hsub, flat);
    double chi2_10 = Chi2(hsub10, flat);
    int ndof = hsub->GetNbinsX();
    delete flat;
    c1->cd(1);
    ltx.DrawLatex(0.65, 0.05, Form("#chi^{2}/ndf = %.3g / %d", chi2, ndof-1));
    if (0)
      cout << "Chi2[5] = " << chi2 << " Chi2[10] = " << chi2_10 
	   << " ndof = " << ndof-1 << " (Prob = " << TMath::Prob(chi2,ndof-6) 
	   << " , " << TMath::Prob(chi2_10,ndof-11) << " )" << endl;
  }
  return c1;
}

double Chi2(TH1* h, TF1* f, TString opt) {
  double chi2 = 0;
  double ndf = 0;
  for (int i=1;i<=h->GetNbinsX();i++) {
    double fy = f->Eval(h->GetBinCenter(i));
    chi2 += pow( (h->GetBinContent(i) - fy)/h->GetBinError(i) , 2.0);
    ndf++;
  }
  ndf -= 1;

  if (opt.Contains("ndf") && ndf > 0)
    chi2 /= ndf;

  if (opt.Contains("prob") )
    return TMath::Prob(chi2, ndf);
 
  return chi2; 
}

double vntvna(double *x, double *par) 
{
  // This function is called many times during a fit. For a 1-d fit,
  // x[0] is the x-position of the nth point,
  // i.e. graph->GetX()[n]. For a single-valued TGraph, it can be used
  // to indicate which point the fitter is currently working on.  The
  // par array contains maxpt free parameters that are continually
  // updated during the fit process.

  // q is the global index
  int q = (int)x[0]; // Assume points are at q + 0.5.
  int i = 999; 
  int j = 999;
  int nPoints = (maxpt*(maxpt+1))/2;
  TF1::fgRejectPoint = kFALSE;

  // Assign i and j from q
  ijFromGlobalIndex(q, i, j);

  if (q < 0 || q >= nPoints) {
    Error("vntvna()", 
	  "Problem finding global index. x[0]=%f, q=%d, i=%d, j=%d",
	  x[0], q, i, j);
    
    return 999;
  }
  
  // Exclude points from fit. Warning: TF1::fgRejectPoint is a global
  // flag checked by all TF1s, beware of side-effects.
  if (j < minFitPtBin || j > maxFitPtBin) {
    TF1::RejectPoint();
    if (0)
      cout << Form("Rejecting x[0],q,i,j = %f %d %d %d", 
		   x[0], q, i, j)
	   << endl;
    return 0;
  }
  
  double v2t = par[i];
  double v2a = par[j];
  int n = par[maxpt];                       // Fourier coeff. index
  int k = par[maxpt+1];                     // cent. bin index
  double fitval = v2t * v2a;

  // Add extra fit term for momentum conservation
  // It is 1/<\sum pt2> * <ptt><pta>
  // From ML & JYO, PRL 106, 102301 (2011) eqs 3-7
  double mptt=0, mpta=0, v1fac=0;
  if (n==1) {
    mptt = MeanPt(i, j, k, "t");
    mpta = MeanPt(i, j, k, "a");
    v1fac = par[maxpt+2];

    fitval += v1fac * mptt * mpta;  
  }
  
  // Odd vnd's < 0 at high pt. // can I replace the below by saying if
  // (vndelta < 0) then make fitval < 0?
  if (j >= PtBin(5., 6.) && (n%2) ) 
    fitval = -fitval;

  return fitval;
}

void VnDelta(int n, const TH1 &h, double &vnd, double &vnd_err, TString opt)
{
  // Extract fourier coefficient n directly from h. If opt is "hi" or
  // "lo", the coefficient is calculated for the points + or - their
  // errors. "hi" and "lo" are meant for an h with systematic (not
  // statistical) errors, so vnd_err is not calculated.
  double VN = 0.0;
  double norm = 0.0;

  // Calculate and set coeff. vnd
  for (int ib=1; ib<=h.GetNbinsX(); ib++) {
    double deltaphi = h.GetBinCenter(ib);
    double weight   = h.GetBinContent(ib);

    if (opt.Contains("hi"))
      weight += h.GetBinError(ib);
    else if (opt.Contains("lo"))
      weight -= h.GetBinError(ib);

    if (opt.Contains("sin") )
      VN += sin(n * deltaphi ) * weight;
    else
      VN += cos(n * deltaphi ) * weight;
    
    norm += weight;
  }
  
  if (norm==0) 
    ::Error("VnDelta()", "Div/0 error");

  VN = VN / norm;

  // // Mom. cons correction?????????????????????????????????
  // if (n==1) VN += 0.002*ptmean[i]*ptmean[j];

  vnd = VN;
  if (opt.Contains("hi") || opt.Contains("lo")) {
    vnd_err = 0;
    return;
  }

  // Statistical uncertainty vnd_err
  double quad_sum_uncertainty = 0.0;
  for (int ib=1; ib<=h.GetNbinsX(); ib++) {
    double deltaphi = h.GetBinCenter(ib);
    double dfdwi = cos(n * deltaphi)/norm - VN/norm;

    if (opt.Contains("sin") )
      dfdwi = sin(n * deltaphi)/norm - VN/norm;

    double sigma = h.GetBinError(ib);
    quad_sum_uncertainty += dfdwi*dfdwi * sigma*sigma;
  }

  vnd_err = TMath::Sqrt(quad_sum_uncertainty);
  return;
}

void VnDelta(int n, const TF1 &f, double &vnd, double &vnd_err, TString opt)
{
  // Extract fourier coefficient n from f.
  double VN = 0.0;
  double norm = 0.0;
  int nSteps = 100;
  double x1=0, x2=0;
  f.GetRange(x1, x2);
  double dx = (x2-x1)/nSteps;

  if (!opt.IsNull() ) cout<<opt.Data()<<endl;

  // Calculate and set coeff. vnd
  for (int ib=0; ib<nSteps; ib++) {

    double deltaphi = x1 + ib*dx; 
    double weight   = f.Eval(deltaphi);

    VN += cos(n * deltaphi ) * weight;
    norm += weight;
  }
  
  if (norm==0) 
    Error("VnDelta()", "Div/0 error");

  VN = VN / norm;
  vnd = VN;
  vnd_err = 0; // for now
  return;
}



// The argument ydiv is the division between pads as a fraction of the
// height.
TCanvas* TwoPanelVert(double ydiv, const char* canvName, const char* canvTitle)
{
  TCanvas* c = new TCanvas(canvName, canvTitle, 800, 800);
  // Extents of pads inside canvas
  double x1=0.01, x2=0.99, y1=0.01, y2=0.99;
  double lm=0.2, rm=0.02;

  // Division between top and bottom pads
  double ysplit = y1 + ydiv*(y2-y1);

  // Divide first, adjust afterward.
  c->Divide(1,2);
  TPad* ptop = (TPad*)(c->GetPad(1));
  TPad* pbot = (TPad*)(c->GetPad(2));

  ptop->SetPad(x1, ysplit, x2, y2);
  ptop->SetMargin(lm, rm, 0.01, 0.01); // left, right, bottom, top.
  ptop->SetFrameLineWidth(2);

  pbot->SetPad(x1, y1, x2, ysplit);
  pbot->SetMargin(lm, rm, 0.4, 0.01); // left, right, bottom, top.
  pbot->SetFrameLineWidth(2);
  pbot->SetTicks();

  c->cd();
  return c;
}

char* trigLabel(int i)
{
  if (i<0 || i>maxpt)
    return Form("Error: no trig pt bin %d", i);
  
  double x1 = gi->GetX()[i] - gi->GetEX()[i];
  double x2 = gi->GetX()[i] + gi->GetEX()[i];

  // int prec1 = 0, prec2 = 0;
  // if (x1 - (int)x1 > 0)
  //   prec1 = 1;
  // if (x2 - (int)x2 > 0)
  //   prec2 = 1;

  return Form("%.3g-%.3g", x1, x2);
}

char* asscLabel(int j)
{
  if (j<0 || j>maxpt)
    return Form("Error: no assc pt bin %d", j);
  
  double x1 = gj->GetX()[j] - gj->GetEX()[j];
  double x2 = gj->GetX()[j] + gj->GetEX()[j];

  /* Stupid coding...
  int prec1 = 0, prec2 = 0;

  if (x1 - (int)x1 > 0)
    prec1 = 1;
  if (x2 - (int)x2 > 0)
    prec2 = 1;
  */

  return Form("%.3g-%.3g", x1, x2);
}

char* centLabel(int k)
{
  if (k<0 || k>maxcent)
    return Form("Error: no cent bin %d", k);
  
  double c1 = gk->GetX()[k] - gk->GetEX()[k];
  double c2 = gk->GetX()[k] + gk->GetEX()[k];

  return Form("%.0f-%.0f%%", c1, c2);
}


void VnLimits(double &minVn_k, double &maxVn_k, int k, int imin, int imax)
{
  // Find limits for each centrality bin k between imin (inclusive)
  // and imax (also inclusive) to plot a common y axis.
  //  for (int k=0; k < maxcent; k++) {
    for (int n=1; n<=5; n++) {
      for (int i=imin; i<=imax; i++) {
	TGraphErrors* gc = VnDeltaVsPtAssoc(i, k, n);
	if (!gc) {
	  Error("VnLimits()", "No graph i %d k %d n %d", i,k,n );
	  continue;
	}
	for (int j=0; j<gc->GetN(); j++) {
	  double vn = gc->GetY()[j];
	  double evn = gc->GetEY()[j];
	  if (vn+evn > maxVn_k) maxVn_k = vn+evn;
	  if (vn-evn < minVn_k) minVn_k = vn-evn;
	}
      }
    }

  return;
}

TCanvas* DrawQ(int k, int n, TString opt)
{
  TStopwatch ts1, ts2;
  int fitCurveColor = kRed-4;
  int divColor = kBlue;
  const char* region = "RIDGE";
  if (opt.Contains("ALL"))
    region = "ALL";

  TString title = Form("globalfit_%d_cent%dto%d",
			   n, centlow[k],centhigh[k]);
  TString name = Form("V%dvsQ_etamin%02d_cent%dto%d",
		      n, (int)(10*minRidgeDEta), 
		      centlow[k],centhigh[k]);
  
  if (!opt.IsNull()) {
    name.Append("_");
    name.Append(opt.Data());
    title.Append("_");
    title.Append(opt.Data());
  }

  TCanvas* c = TwoPanelVert(0.4, name.Data(), title.Data());
  int cwidth = opt.Contains("highptfit") ? 500 : 1200; 
  c->SetCanvasSize(cwidth, 450);
  c->SetWindowSize(cwidth+50, 500);

  if (opt.Contains("highptfit")) {
    c->GetPad(1)->SetLeftMargin(0.2);
    c->GetPad(2)->SetLeftMargin(0.2);
    c->GetPad(1)->SetRightMargin(0.05);
    c->GetPad(2)->SetRightMargin(0.05);

  }

  TLegend* lq = new TLegend(0.14, 0.56, 0.25, 0.83);
  lq->SetFillColor(kNone);
  lq->SetBorderSize(0);

  // g:  global VnDelta vs global index (stat errors).
  // gs: global VnDelta vs global index (sys errors).
  // fn: global fit TF1.
  // r:  ratio of g/fit.
  TString opt1 = "";
  if (opt.Contains("ALL"))
    opt1.Append("_ALL");
  if (opt.Contains("ptcons"))
    opt1.Append("_ptcons");

  TGraphErrors* g  = VnDVsQ(n, k, region, cfDef, opt1);
  opt1.Append("_sys");
  TGraphErrors* gs = VnDVsQ(n, k, region, cfDef, opt1);

  int ptaLo = 0, ptaHi = 999;
  if (opt.Contains("lowptfit")) {
    ptaHi = PtBin(3., 4.);
  }
  if (opt.Contains("highptfit")) {
    ptaLo = PtBin(5., 6.);
  }
  // If custom pta range was passed in, e.g.
  // Form("ptabin%dto%d", PtBin(1.5, 2.0), PtBin(8, 15))
  if (opt.Contains("ptabin")) {
    TRegexp re("ptabin[0-9]+to[0-9]+");
    TRegexp re1("ptabin[0-9]+");
    TRegexp re2("to[0-9]+");
    TString rx = opt(re), pta1s = rx(re1), pta2s = rx(re2);
    pta1s.ReplaceAll("ptabin","");
    pta2s.ReplaceAll("to","");
    ptaLo = pta1s.Atoi();
    ptaHi = pta2s.Atoi();
  }

  TF1* fn    = GlobalFitVsQ(n, k, ptaLo, ptaHi, region, cfDef, "");
  TF1* fn_hi = GlobalFitVsQ(n, k, ptaLo, ptaHi, region, cfDef, "hi_sys");
  TF1* fn_lo = GlobalFitVsQ(n, k, ptaLo, ptaHi, region, cfDef, "lo_sys");

  // Systematic band on ratio line at unity
  TH1F *r_sys = new TH1F(Form("r_sys_%s", title.Data()), "r_sys", 
			 g->GetN()-0.5, 0, g->GetN()-0.5);
  utils.set_hist_props(r_sys, kBlack, kGray, kBlack, kDot, 1.0);

  TGraphErrors* r = new TGraphErrors();
  TString rname(g->GetName());
  rname.Append("_ratio");
  r->SetName(rname.Data());

  if (!fn) {
    Error("DrawQ()","Problem getting/creating global fit function");
    gSystem->Exit(-1);
  }

  // Set visual properties
  utils.set_tgraph_props(g, kBlack, kBlack, kFullSquare, 1.0);
  utils.set_tgraph_props(gs, kBlack, kBlack, kDot, 0.5);
  gs->SetFillColor(kGray);
  utils.set_tgraph_props(r, kBlack, kBlack, kFullCircle, 1.0);
  g->SetLineWidth(1);
  r->SetLineWidth(1);
  fn->SetLineWidth(1);
  fn->SetLineColor(fitCurveColor);
  fn_hi->SetLineStyle(kDashed);
  fn_lo->SetLineStyle(kDashed);
  fn_hi->SetLineWidth(1);
  fn_lo->SetLineWidth(1);

  // Compute ratio for r graph
  for (int q=0; q<g->GetN(); q++) {
    int ii=0, jj=0;
    ijFromGlobalIndex(q, ii, jj);

    double fy = fn->Eval(q+0.5);
    double y = fy==0? 1e9 : g->GetY()[q]/fy;
    double rerr = fn_hi->Eval(q+0.5) / fy - 1;

    if (opt.Contains("highptfit") && n%2 && g->GetY()[q]!=0)
      y = fy/g->GetY()[q];

    // No error band when no fit point!
    if (jj < ptaLo || jj > ptaHi)
      rerr = 0;

    r->SetPoint(q, g->GetX()[q], y);
    r->SetPointError(q, 0.0, TMath::Abs(g->GetEY()[q]/fy));

    r_sys->SetBinContent(q+1, 1.0);
    r_sys->SetBinError(q+1, rerr);
  }

  // Set up plot ranges for top panel (g and f)
  int nx = g->GetN();
  double x1,y1,x2,y2;
  double sf = 1; // scale factor
  fn->GetRange(x1,x2);
  if (opt.Contains("highptfit")) {
    x1 = 45;
  }
  y1 = sf*fn->GetMinimum(x1,x2);
  y2 = sf*fn->GetMaximum(x1,x2);
  double marg = 0.4*(y2-y1);
  //  if (opt.Contains("lowptfit") || opt.Contains("highptfit") ) marg = 2*(y2-y1);
  x1 -= 0.5;
  x2 += 0.5;
  y1 -= marg;
  y2 += marg;

  // Set up plot ranges for lower panel (r)
  double ry1=0.21, ry2=1.79;
  if (opt.Contains("highptfit")) {
    ry1 = 0.01;
    ry2 = 1.99;
  }

  // Get / make frame histos for top and bottom
  TH2F* hf1 = 0; // (TH2F*)gDirectory->Get("qhf1");
  TH2F* hf2 = 0; // (TH2F*)gDirectory->Get("qhf2");
  TString name1 = Form("qhf1_%s", title.Data());
  TString name2 = Form("qhf2_%s", title.Data());
  // TString sv = n==1 ? "#LTcos(#Delta#phi)#GT" :
  // Form("#LTcos(%d#Delta#phi)#GT", n);
  TString sv = Form("V_{%d#Delta}", n);
  TString sgf = Form("(v_{%d}^{t}v_{%d}^{a})_{fit}", n, n);

  if (!hf1)
    hf1 = new TH2F(name1.Data(), "hf1", nx+1, x1, x2, 1000, y1, y2);
  if (!hf2)
    hf2 = new TH2F(name2.Data(), "hf2", nx+1, x1, x2, 1000, ry1, ry2);

  // Axes
  TAxis* ax1 = hf1->GetXaxis();
  TAxis* ay1 = hf1->GetYaxis();
  TAxis* ax2 = hf2->GetXaxis();
  TAxis* ay2 = hf2->GetYaxis();
  hf1->SetBins(nx+1, x1, x2, 100, y1, y2);
  hf2->SetBins(nx+1, x1, x2, 100, ry1, ry2);
  hf1->SetAxisRange(x1, x2, "X");  
  hf1->SetAxisRange(y1, y2, "Y");  
  hf2->SetAxisRange(x1, x2, "X");  
  hf2->SetAxisRange(ry1, ry2, "Y");  
  ax1->SetLabelOffset(0.02);
  ax1->SetLabelSize(0.04);  ay1->SetLabelSize(0.08);
  ax2->SetLabelSize(0.18);   ay2->SetLabelSize(0.12);
  ax1->SetNdivisions(1);   ax2->SetNdivisions(1);
  ay1->SetNdivisions(208);  ay2->SetNdivisions(opt.Contains("highptfit")?104:210);
  ax1->SetTickLength(0.02);

  //ay1->SetTitle(Form("%s and #color[%d]{%s}",sv.Data(), fitCurveColor, sgf.Data()));
  ay1->SetTitle(Form("%s", sv.Data()));
  ay2->SetTitle(Form("#frac{%s}{fit}", sv.Data()));
  //  ay2->SetTitle(Form("#frac{%s}{%s}", sv.Data(), sgf.Data()));
  // ax2->SetTitle(Form("#splitline{associated p_{T}^{a}}{#color[%d]{   trigger p_{T}^{t}}} [GeV/c]", divColor));
  ax2->SetTitle(Form("#color[%d]{p_{T}^{t}}, p_{T}^{a} [GeV/c]", divColor));
  ay1->SetTitleOffset(opt.Contains("highptfit") ? 1. : 0.45);
  ay1->SetTitleSize(0.09);
  ay1->CenterTitle();
  ay2->SetTitleOffset(opt.Contains("highptfit") ? 0.6 : 0.26);
  ay2->SetTitleSize(0.14);
  ay2->CenterTitle();
  ax2->SetLabelOffset(2.1);
  ax2->SetTitleOffset(1.3);
  ax2->SetTitleSize(0.14);
  ax2->CenterTitle();

  double leftMarg = opt.Contains("highptfit") ? 0.2 : 0.08;
  c->cd(1);
  gPad->SetFrameLineWidth(1);
  gPad->SetLeftMargin(leftMarg);
  hf1->Draw();
  TLine zero;
  zero.DrawLine(x1, 0., x2, 0.);

  
  // NEW========================================
  for (int i=0; i<gi->GetN(); i++) {
    int q1 = GlobalIndex(i, 1), q2 = GlobalIndex(i, i);
    TF1* ifn = 0, *ifn_hi = 0, *ifn_lo = 0;
    //    TF1* _ifn = 0, *_ifn_hi = 0, *_ifn_lo = 0; // temp, unscaled.

    // This works, but has lines between \ptt groups :(
    // TH1F *hfn=0, *hfn_hi=0, *hfn_lo=0;
    // int nq = g->GetN();
    // hfn = new TH1F(Form("%s_%d", fn->GetName(), i), Form("%s_%d", fn->GetName(), i),
    // 		   nq, 0, nq);
    // hfn->SetLineColor(kMagenta);
    // for (int q=0; q<g->GetN(); q++) {
    //   hfn->SetBinContent(q+1, sf*fn->Eval(q+0.5));
    // }

    ifn = (TF1*)       fn->Clone(Form("%s_%d", fn->GetName(), i));
    ifn_hi = (TF1*) fn_hi->Clone(Form("%s_%d", fn_hi->GetName(), i));
    ifn_lo = (TF1*) fn_lo->Clone(Form("%s_%d", fn_lo->GetName(), i));

    ifn->SetName(Form("%s_%d", fn->GetName(), i));
    ifn_hi->SetName(Form("%s_%d", fn_hi->GetName(), i));
    ifn_lo->SetName(Form("%s_%d", fn_lo->GetName(), i));

    double r1 = q1-0.9, r2 = q2+0.92;
    ifn->SetRange(r1, r2);
    ifn_hi->SetRange(r1, r2);
    ifn_lo->SetRange(r1, r2);


    // for (int ip=0; ip<ifn->GetNpar(); ip++) {
    //   double par = ifn->GetParameter(ip);
    //   ifn->SetParameter(ip, 1000*par);
    //   cout<<par<< " " << ifn->GetParameter(ip) << endl;
    // }


    // ifn    = new TF1(Form("%s_%d",    fn->GetName(), i), Form("100*%s",   _ifn->GetName()));
    // ifn_hi = new TF1(Form("%s_%d", fn_hi->GetName(), i), Form("100*%s", _ifn_hi->GetName()));
    // ifn_lo = new TF1(Form("%s_%d", fn_lo->GetName(), i), Form("100*%s", _ifn_lo->GetName()));

    ifn_hi->Draw("same");
    ifn_lo->Draw("same");
    ifn->Draw("same");
    //    hfn->Draw("same");
  }
  // NEW========================================

  for (int i=0; i<gi->GetN(); i++) {
    char* newName = Form("%s_ptt%d", g->GetName(), i);
    TGraphErrors* gcol = (TGraphErrors*)g->Clone(newName);
    int col = PtColor(ptlow[i], pthigh[i]);
    utils.set_tgraph_props(gcol, kBlack, col, kFullSquare, 1.0);
    gcol->SetLineWidth(1);

    newName = Form("%s_ptt%d", gs->GetName(), i);
    TGraphErrors* g_sy = (TGraphErrors*)gs->Clone(newName);
    int col_sy = PtColorPale(ptlow[i], pthigh[i]);
    utils.set_tgraph_props(g_sy, kBlack, col_sy, kDot, 0.5);
    g_sy->SetFillColor(col_sy);
    
    for (int q=0; q<gcol->GetN(); q++) {
      int i2=0, j=0;
      ijFromGlobalIndex(q,i2,j);
      if (i2!=i) {
	gcol->SetPoint(q, gcol->GetX()[q], -99.);
	g_sy->SetPoint(q, gcol->GetX()[q], -99.);
      }

      // Scale points by sf
      gcol->SetPoint(q, gcol->GetX()[q], sf*gcol->GetY()[q]);
      g_sy->SetPoint(q, g_sy->GetX()[q], sf*g_sy->GetY()[q]);
      gcol->SetPointError(q, gcol->GetEX()[q], sf*gcol->GetEY()[q]);
      g_sy->SetPointError(q, g_sy->GetEX()[q], sf*g_sy->GetEY()[q]);

    }
    g_sy->Draw("e2psame");
    gcol->Draw("ezpsame");
  }
  
  c->cd(2);
  gPad->SetFrameLineWidth(1);
  gPad->SetLeftMargin(leftMarg);
  hf2->Draw();
  r_sys->Draw("e2psame");
  
  // Overlay lines, arrows, text, etc.
  TLine div, div_a; // division ticks for trigger and assoc
  TArrow a;
  double tipSize = 0.008; // size of arrowhead
  a.SetLineWidth(1);
  div.SetLineWidth(2);
  div.SetLineColor(fitCurveColor);
  div.DrawLine(x1, 1.0, x2, 1.0);
  r->Draw("ezpsame");
  gList->Add(r);

  for (int i=0; i<gi->GetN(); i++) {
    char* newName = Form("%s_ptt%d", r->GetName(), i);
    TGraphErrors* rcol = (TGraphErrors*)r->Clone(newName);
    int col = PtColor(ptlow[i], pthigh[i]);
    utils.set_tgraph_props(rcol, kBlack, col, kFullCircle, 1.0);
    for (int q=0; q<rcol->GetN(); q++) {
      int i2=0, j=0;
      ijFromGlobalIndex(q,i2,j);
      if (i2!=i) rcol->SetPoint(q, rcol->GetX()[q], -99.);
    }
    rcol->SetLineWidth(1);
    rcol->Draw("ezpsame");
  }
  
  div.SetLineColor(divColor);
  div.SetLineStyle(kSolid);

  // Label centrality and fourier index
  c->cd(1);
  ltx.SetTextSize(0.14);
  ltx.DrawLatex(opt.Contains("highptfit")?0.25:0.14, 0.85, Form("n = %d", n));

  if (!opt.Contains("highptfit")) {
  if (n<=2) {
    ltx.SetTextSize(0.1);
    if (ispp)
    ltx.DrawLatex(0.14, 0.70, Form("pp 2.76 TeV"));
    else {
      ltx.DrawLatex(0.14, 0.70, Form("Pb-Pb 2.76 TeV"));
      ltx.DrawLatex(0.14, 0.61, Form("%d-%d%%", centlow[k],centhigh[k]));
    }
    if (opt.Contains("ptcons")) {
    ltx.DrawLatex(0.14, 0.58, Form("#LT#sum p_{T}^{2}#GT = 4.76 GeV^{2}"));
    }

  }
  if (n==3) {
    lq->AddEntry(g, Form("V_{n#Delta}"), "epl");
    lq->AddEntry(fn, Form("fit"), "l");
    lq->Draw();
  }
  // Always show
    ltx.SetTextSize(0.07);
  double ptaFitLow = ptlow[ptaLo];
  double ptaFitHigh = ptaHi<999? pthigh[ptaHi] : pthigh[maxpt-1];
  //  if (!opt.Contains("highptfit"))
    ltx.DrawLatex(0.14, 0.09, 
		  Form("Global fit range: %3.2g < p_{T}^{a} < min(p_{T}^{t}, %3.2g GeV/c)",
		       ptaFitLow, ptaFitHigh));
  }

  for (int q=0; q<g->GetN(); q++) {
    int i=0, j=0;
    ijFromGlobalIndex(q,i,j);
    double drop = 0.35;
    double ptt = gi->GetX()[i] + gi->GetEX()[i];
    double pta = gj->GetX()[j] + gj->GetEX()[j];
    double tickLength = 0.03*(y2-y1);

    // Draw divisions between trigger pt intervals...
    bool skipTicks = opt.Contains("highptfit") && (pta < 5.);
    if (i==j && !skipTicks) {
      //      if (ptt > 5) {
	c->cd(1);
	div.DrawLine(q+1, y1, q+1, y1+tickLength);
	c->cd(2);
	tickLength = drop-0.1;
	div.DrawLine(q+1, ry2-tickLength, q+1, ry2+tickLength);
	div.DrawLine(q+1, ry1-tickLength, q+1, ry1+tickLength);
	//      }

      // Alignment: align = 10*HorizontalAlign + VerticalAlign
      // For horizontal alignment the following convention applies:
      // 1=left adjusted, 2=centered, 3=right adjusted
      // For vertical alignment the following convention applies:
      // 1=bottom adjusted, 2=centered, 3=top adjusted
      ltx2.SetTextAlign(23);
      ltx2.SetTextSize(0.12);
      ltx2.SetTextColor(divColor);

      if (ptt != 0.75) // too cluttered otherwise 
	ltx2.DrawLatex(q+1, ry1-drop, Form("%.2g", ptt) ); // pt_trig labels
      }
    
    if (pta < ptt && (pta == (int)pta || pta == 0.5) && !skipTicks) {
	div_a.DrawLine(q+1, ry1-drop/8, q+1, ry1+drop/4);
	ltx2.SetTextColor(kBlack);
	ltx2.SetTextFont(52); // italic
	ltx2.DrawLatex(q+1, ry1-drop/4, Form("%.2g", pta) ); // pt_assoc labels
	ltx2.SetTextColor(divColor);
	ltx2.SetTextFont(62); // 62 is default
    }

    // Draw arrows to off-scale points
    c->cd(1);
    double y = g->GetY()[q];
    double ymax = hf1->GetYaxis()->GetXmax();
    double ymin = hf1->GetYaxis()->GetXmin();
    double len = 0.06*(ymax-ymin);
    double gap = 0.01*(ymax-ymin);

    if (y > ymax)
      a.DrawArrow(q+0.5, ymax-len-gap, q+0.5, ymax-gap, tipSize, ">");
    if (y < ymin)
      a.DrawArrow(q+0.5, ymin+len+gap, q+0.5, ymin+gap, tipSize, ">");

    // Do the same for the ratio.
    // If fit has not been applied to some points, don't draw the arrow.
    c->cd(2);
    //    gPad->SetTickx();
    y = r->GetY()[q];
    ymax = hf2->GetYaxis()->GetXmax();
    ymin = hf2->GetYaxis()->GetXmin();
    len = 0.06/0.4*(ymax-ymin); // 0.4 from TwoPanelVert ctor
    gap = 0.02*(ymax-ymin);

    if (0) { // Draw at upper or lower edge of panel
      if (y > ymax)
	a.DrawArrow(q+0.5, ymax-len-gap, q+0.5, ymax-gap, tipSize, ">");
      if (y < ymin)
	a.DrawArrow(q+0.5, ymin+len+gap, q+0.5, ymin+gap, tipSize, ">");
    }
    // Draw with arrow based from ref. line.
    // No arrow when no fit point!
    if (j < ptaLo || j > ptaHi)
      continue;
    if (y > ymax)
      a.DrawArrow(q+0.5, 1.0, q+0.5, 1.0+len, tipSize, ">");
    if (y < ymin)
      a.DrawArrow(q+0.5, 1.0, q+0.5, 1.0-len, tipSize, ">");
  } // end q loop
  
  //  ax2->LabelsOption("h");
  c->cd();
  c->Update();
  return c;
}

TCanvas* DrawVn(int k, int ptt1, int ptt2)
{
  // Draw Fourier coefficients on one canvas for each centrality bin k
  // for trigger pt bins from ptt1 to (and including) ptt2.
  // The x-axis is associated pt.
  const char* name = Form("vndelta_etamin%02d_cent%dto%d_trig%dto%d", 				 
			  (int)(10*minRidgeDEta), centlow[k],centhigh[k], ptt1, ptt2);
  TCanvas* cv = new TCanvas(name, name, 1200, 500);
  TLegend* lv = new TLegend(0.09, 0.72, 0.27, 0.88); // w.r.t. cv
  lv->SetFillColor(kNone);
  lv->SetBorderSize(0);

  double lv2bottom = ptt2<5? 0.65 : 0.8;
  TLegend* lv2 = new TLegend(0.82, lv2bottom, 0.94, 0.95, "p_{T}^{t} (GeV/c)");
  lv2->SetFillColor(kNone);
  lv2->SetBorderSize(0);
  utils.padsetup(cv, 5, 1, "", 0.12, 0.01, 0.03, 0.2);

  // Find limits
  double minVn_k=1000, maxVn_k=-1000;	  
  VnLimits(minVn_k, maxVn_k, k, ptt1, ptt2);

    for (int n=1; n<=5; n++) {
      cv->cd(n);
      double vertMargin = 0.1*(maxVn_k - minVn_k);
      double xmax = gj->GetX()[ptt2] + gj->GetEX()[ptt2];
      double xmin =  -0.43;
      TH1F* hf = gPad->DrawFrame(xmin, minVn_k - vertMargin,
				 xmax, maxVn_k + vertMargin);
      int nXticks = 404;
      int nYticks = 210;
      if (n==1) {
	utils.make_nice_axes(cv, hf, 2., nXticks, nYticks, -0.02, 0.02);
	hf->SetLabelOffset(-0.01, "X");
      }
      else {
	utils.make_nice_axes(cv, hf, 3.4, nXticks, nYticks, -0.075, 0.0);
	hf->SetLabelOffset(-0.057, "X");
	hf->SetTickLength(0.02, "X");
	hf->SetTickLength(0.05, "Y");
      }
      hf->GetXaxis()->SetTicks("+-");
      hf->GetYaxis()->SetTicks("+-");
      
      // Loop backwards to improve visibility of low-pt points
      for (int i=ptt2; i>=ptt1; i--) {
	
	TGraphErrors *gc=0;
	gc = VnDeltaVsPtAssoc(i, k, n);
	if (!gc) {
	  Error("DrawVn()", "Problem finding graph n%d, i%d, k%d" ,n,i,k);
	  continue;
	}
	
	utils.set_tgraph_props(gc, colorsDark[i], colorsDark[i], kFullSquare, 1.6);
	
	TGraphErrors* gfit = GlobalFitVsPtAssoc(i, k, n);
	
	gfit->Draw("plsame");
	gc->Draw("epsame");

	if (n==1)
	  lv2->AddEntry(gc,trigLabel(i), "p");

      }
    }

    cv->cd();
    TPad* overlay = new TPad("overlay", "overlay", 0,0,1,1);
    overlay->SetFillStyle(4000); // transparent
    overlay->Draw();
    overlay->cd();

    ltx.SetTextSize(0.07);
    lv2->Draw();
    ltx.DrawLatex(0.14, 0.76, Form("%s", centLabel(k)));

    for (int n=1; n<=4; n++) {
      cv->cd(n);
      ltx.SetTextSize(n==1? 0.1 : 0.15);
      ltx.DrawLatex((n==1)? 0.45 : 0.1, 0.86, Form("V_{%d#Delta}", n));
    }
    cv->cd(5);
    ltx.DrawLatex(0.72, 0.86, Form("V_{%d#Delta}", 5));
    overlay->cd();
    ltx.SetTextSize(0.075);
    ltx.DrawLatex(0.4, 0.05, "associated p_{T} [GeV/c]");
    ltx.SetTextAngle(90);
    ltx.DrawLatex(0.03, 0.4, "V_{n#Delta}(p_{T}^{assoc.})");
    ltx.SetTextAngle(0);

  return cv;
}

TCanvas* DrawVnVsPtTrig(int k, int npta, int ptabins[], TString opt)
{
  TString title = Form("VnDeltaVsPtTrig_cent%dto%d_pTassoc%.2gto%.2g",
		       centlow[k],centhigh[k],
		       ptlow[ptabins[0]], pthigh[ptabins[npta-1]]);
  if (!opt.IsNull()) {
    title.Append("_");
    title.Append(opt);
  }
  TCanvas* cv = new TCanvas(title.Data(), title.Data(), 1400, 500);
  double lv2bot = 0.9 - 0.07*npta;
  TLegend* lv2 = new TLegend(0.82, lv2bot, 0.98, 0.95, "p_{T}^{a} (GeV/c)");
  lv2->SetFillColor(kNone);
  lv2->SetBorderSize(0);
  utils.padsetup(cv, 5, 1, "", 0.12, 0.01, 0.03, 0.2);

  // Start without knowing y-limits, adjust below.
  double minVn_k=-1., maxVn_k=1.;	  
  TObjArray* graphs = new TObjArray();
  
  for (int n=1; n<=5; n++) {
    cv->cd(n);
    double vertMargin = 0.1*(maxVn_k - minVn_k);
    double xmax = 11.1;
    double xmin = -0.43;

    TH1F* hf = gPad->DrawFrame(xmin, minVn_k - 4*vertMargin,
			       xmax, maxVn_k + vertMargin);
    TLine l;
    l.DrawLine(xmin, 0, xmax, 0);
    
    int nXticks = 404;
    int nYticks = 210;
    if (n==1) {
      utils.make_nice_axes(cv, hf, 2., nXticks, nYticks, -0.02, 0.02);
      hf->SetLabelOffset(-0.01, "X");
    }
    else {
      utils.make_nice_axes(cv, hf, 3.4, nXticks, nYticks, -0.075, 0.0);
      hf->SetLabelOffset(-0.057, "X");
      hf->SetTickLength(0.02, "X");
      hf->SetTickLength(0.05, "Y");
    }
    hf->GetXaxis()->SetTicks("+-");
    hf->GetYaxis()->SetTicks("+-");

    for (int jpta=0; jpta<npta; jpta++) {
      int j = ptabins[jpta];
      TGraphErrors *gc = VnDeltaVsPtTrig(j, k, n);
      if (!gc) {
	Error("DrawVnVsPtTrig()", "Problem finding gc n%d, i%d, k%d" ,n,j,k);
	continue;
      }
      TGraphErrors *gs = VnDeltaVsPtTrig(j, k, n, "sys");
      if (!gs) {
	Error("DrawVnVsPtTrig()", "Problem finding gs n%d, i%d, k%d" ,n,j,k);
	continue;
      }

      int col = PtColor(ptlow[j], pthigh[j]);
      utils.set_tgraph_props(gc, col, col, kFullSquare, 1.2);
      gs->SetFillColor(PtColorPale(ptlow[j], pthigh[j]));
      TGraphErrors* open = (TGraphErrors*)gc->Clone();
      utils.set_tgraph_props(open, kBlack, kBlack, kOpenSquare, 1.2);
	
      gs->Draw("e2psame");
      gc->Draw("epsame");
      open->Draw("epsame");

      if (n==1)
	lv2->AddEntry(gc,asscLabel(j), "elp");

      graphs->Add(gc); // store graph array to find y limits for plot
    }
  }

  // Adjust y limits
  double newYmin = 0, newYmax = 0, newSpace = 0;
  for (int m=0; m<graphs->GetEntries(); m++) {
    TGraphErrors* tg = (TGraphErrors*)graphs->At(m);
    double yhi = TMath::MaxElement(tg->GetN()-1, tg->GetY());
    double ylo = TMath::MinElement(tg->GetN()-1, tg->GetY());
    double ey = TMath::MaxElement(tg->GetN()-1, tg->GetEY());
    ylo -= ey;
    if (yhi > newYmax) newYmax = yhi;
    if (ylo < newYmin) newYmin = ylo;
  }
  newSpace = 0.1*(newYmax - newYmin);
  for (int n=1; n<=5; n++) {
    cv->cd(n);
    TH1F* hf = (TH1F*)gPad->FindObject("hframe");
    hf->GetYaxis()->SetRangeUser(newYmin-3*newSpace, newYmax+newSpace);
    gPad->Modified(); // Signal to redraw
  }
 
  cv->cd();
  TPad* overlay = new TPad("overlay", "overlay", 0,0,1,1);
  overlay->SetFillStyle(4000); // transparent
  overlay->Draw();
  overlay->cd();
  lv2->Draw();
  ltx.SetTextSize(0.05);
  ltx.DrawLatex(0.14, 0.88, Form("%s 2.76 TeV", ispp?"pp":"Pb-Pb"));
  ltx.SetTextSize(0.07);
  if (!ispp)
    ltx.DrawLatex(0.14, 0.8, Form("%s", centLabel(k)));

  for (int n=1; n<=5; n++) {
    cv->cd(n);
    ltx.SetTextSize(n==1? 0.1 : 0.15);
    if (n==1)
      ltx.DrawLatex(0.44, 0.24, Form("n=%d", n));
    else
      ltx.DrawLatex(0.04, 0.24, Form("n=%d", n));

  }
  overlay->cd();
  ltx.SetTextSize(0.075);
  ltx.DrawLatex(0.48, 0.05, "trigger p_{T}^{t} [GeV/c]");
  ltx.SetTextAngle(90);
  ltx.DrawLatex(0.06, 0.48, "V_{n#Delta}  [10^{-2}]");
  ltx.SetTextAngle(0);
  return cv;
}

TGraphErrors* AgreementVsPtSum(int k, int n)
{
  TGraphErrors* gVn = VnDVsQ(n, k);
  TF1* fVn = GlobalFitVsQ(n,k);
  TGraphErrors* g = new TGraphErrors();
  int N = gVn->GetN(); // And f better have the same N.
  int i = 0, j = 0;
  if(N<=0) {
    Error("AgreementVsPtSum()", "Invalid number of points: %d", N);
    gSystem->Exit(1);
  }

  for (int q=0; q<N; q++) {
    ijFromGlobalIndex(q,i,j);
    double gy = gVn->GetY()[q];
    double fy = fVn->Eval(q+0.5);
    double r = TMath::Abs((gy-fy)/gy); 
    g->SetPoint(q, PtBinCenterSum(i, j), r);
  }
  g->SetMarkerStyle(kFullCircle);
  return g;
}

TGraph2D* Agreement2D(int k, int n)
{
  TGraphErrors* gVn = VnDVsQ(n, k);
  TF1* fVn = GlobalFitVsQ(n,k);
  TGraph2D* g = new TGraph2D();
  int N = gVn->GetN(); // And f better have the same N.
  int i = 0, j = 0;
  if(N<=0) {
    Error("AgreementVsPtSum()", "Invalid number of points: %d", N);
    gSystem->Exit(1);
  }

  for (int q=0; q<N; q++) {
    ijFromGlobalIndex(q,i,j);
    double gy = gVn->GetY()[q];
    double fy = fVn->Eval(q+0.5);
    double r = TMath::Abs((gy-fy)/fy); 
    g->SetPoint(q, ptmean[j], ptmean[i], r);
  }
  return g;
}

TH2F* Agreement2DHist(int k, int n)
{
  TGraphErrors* gVn = VnDVsQ(n, k);
  TF1* fVn = GlobalFitVsQ(n,k);

  // Pt binning array
  double ptarr[100];
  for (int i=0; i<maxpt; i++) ptarr[i] = gi->GetX()[i] - gi->GetEX()[i];
  ptarr[maxpt] = gi->GetX()[maxpt-1] + gi->GetEX()[maxpt-1];

  int N = gVn->GetN(); // And f better have the same N.
  int i = 0, j = 0;
  if(N<=0) {
    Error("AgreementVsPtSum()", "Invalid number of points: %d", N);
    gSystem->Exit(1);
  }
  const char* name = Form("agmt_cent%dto%d_n%d", 
			  centlow[k],centhigh[k], n);
  TH2F* h = new TH2F(name, name, maxpt, ptarr, maxpt, ptarr);
  h->SetTitle("title;p_{T}^{a};p_{T}^{t};");

  for (int q=0; q<N; q++) {
    ijFromGlobalIndex(q,i,j);
    double gy = gVn->GetY()[q];
    double fy = fVn->Eval(q+0.5);
    
    // double r = fy==0? 0 : TMath::Abs(gy-fy) / fy; 
    double r = gy==0? 0 : TMath::Abs((gy-fy) / gy); 
    //    double r = fy==0? 0 : TMath::Abs(gy / fy); 
    h->SetBinContent(j+1, i+1, r);
  }
  return h;
}

TH2* PttPtaHist()
{
  // Pt binning array
  double ptarr[100];
  for (int i=0; i<maxpt; i++) ptarr[i] = gi->GetX()[i] - gi->GetEX()[i];
  ptarr[maxpt] = gi->GetX()[maxpt-1] + gi->GetEX()[maxpt-1];

  int i = 0, j = 0;
  const char* name = Form("PttPtaHist");
  TH2F* h = new TH2F(name, name, maxpt, ptarr, maxpt, ptarr);
  h->SetTitle("title;p_{T}^{a} (GeV/c);p_{T}^{t} (GeV/c);");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  int N = maxpt * (maxpt + 1) / 2;
  for (int q=0; q<N; q++) {
    ijFromGlobalIndex(q,i,j);
    h->SetBinContent(j+1, i+1, 1.);
  }
  return h;
}

TGraphErrors* Chi2VsTrigPt(int k, int n)
{
  TGraphErrors* g = new TGraphErrors();
  for (int i=0; i<maxpt; i++) {
    g->SetPoint(i, ptmean[i], ReducedChi2(i, k, n));
  }
  g->SetMarkerStyle(kFullCircle);
  return g;
}

TGraphErrors* 
CorrFnChi2VsTrigPt(int j, int k, int n1, int n2, TString opt)
{
  TGraphErrors* g = new TGraphErrors();
  int ip=0;
  for (int i=j; i<maxpt; i++) {
    TH1* h = Hist("RIDGE", "cA", i, j, k);
    double x = MeanPt(i, j, k, "t");
    double y = Chi2(h, HarmonicSum(h, n1, n2, "global"), opt);
    g->SetPoint(ip, x, y);
    ip++;
  }
  g->SetMarkerStyle(kFullCircle);
  return g;
}

double ReducedChi2(int i, int k, int n)
{
  TGraphErrors* g = VnDeltaVsPtAssoc(i, k, n);
  TGraphErrors* f = GlobalFitVsPtAssoc(i, k, n);
  int N = g->GetN(); // And f better have the same N.
  double chi2 = 0;

  // Value (or one-sigma error if opt=="err") of v_n
  //  double vnerr = Bestvn(k,n,i,"err");
  double vnerr = vnGF(n,k,i,0,999,"RIDGE","cA","err");

  if(N<=0 || N!= f->GetN()) {
    Error("ReducedChi2()", "Invalid number of points: %d", N);
    gSystem->Exit(1);
  }

  for (int j=0; j<N; j++) {
    double errj = g->GetEY()[j];

    chi2 += 
      (g->GetY()[j] - f->GetY()[j]) * (g->GetY()[j] - f->GetY()[j]) / 
      pow(errj + vnerr, 2);
  }
  chi2 /= N;
  return chi2;
}

TGraphErrors* VnDeltaVsPtAssoc(int i, int k, int n, TString opt)
{
  // Create a TGraph from VnDelta points
  TGraphErrors* gALL = VnDVsQ(n, k, "RIDGE", "cA", opt);
  TGraphErrors* g = new TGraphErrors();
  TString gname(Form("VnDeltaVsPtAssoc_ptt%d_cent%d_n%d",i,k,n));
  gname.Append(opt);
  g->SetName(gname.Data());

  for (int j=0; j<=i; j++) {
    int q = GlobalIndex(i, j);
    g->SetPoint(j, gj->GetX()[j], gALL->GetY()[q] );
    g->SetPointError(j, 0,  gALL->GetEY()[q] );
    utils.set_tgraph_props(g, colorsDark[i], colorsDark[i], kFullSquare, 1.0); 
    g->SetLineWidth(2);
  }
  
  return g;
}

TGraphErrors* VnDeltaVsPtTrig(int j, int k, int n, TString opt)
{
  // Create a TGraph from VnDelta points
  TGraphErrors* gALL = VnDVsQ(n, k, "RIDGE", "cA", opt);
  TGraphErrors* g = new TGraphErrors();
  TString gname(Form("VnDeltaVsPtTrig_pta%d_cent%dto%d_n%d",
		     j,centlow[k],centhigh[k],n));
  gname.Append(opt);
  g->SetName(gname.Data());

  int ip=0;
  for (int i=j; i<maxpt; i++) {
    int q = GlobalIndex(i, j);
    
    double xerr = 0.;
    if (opt.Contains("sys")) {
      if (gi->GetEX()[i] < 0.126)
	xerr = gi->GetEX()[i];
      else
	xerr = 0.250;
    }
    double mpt = MeanPt(i, j, k, "t"); 
    double sf = 100; // Scale y-value; make sure to label accordingly!!
    g->SetPoint(ip, mpt, sf*gALL->GetY()[q] ); 
    g->SetPointError(ip, xerr,  sf*gALL->GetEY()[q] );

    utils.set_tgraph_props(g, colorsDark[j], colorsDark[j], kFullSquare, 1.0); 
    g->SetLineWidth(2);
    ip++;
  }
  
  return g;
}

TGraphErrors* VnDeltaVsN(int i, int j, int k, int nMax, TString opt)
{
  // Create a TGraph from VnDelta points
  // Stat errors by default. Call with opt "sys" to get sys errors.
  // Call with opt "sine" to get residual <sin n Delta phi> for sys. est.
  int q = GlobalIndex(i, j);
  TString gname(Form("VnDeltaVsN_ptt%.2gto%.2g_pta%.2gto%.2g_cent%dto%d",
		     ptlow[i],pthigh[i],ptlow[j],pthigh[j],
		     centlow[k],centhigh[k]));
  if (!opt.IsNull())  {
    gname.Append("_");
    gname.Append(opt);
  }
  TGraphErrors* g = (TGraphErrors*) gList->FindObject(gname.Data());
  if (g) {
    return g;
  } 
  else {
    g = new TGraphErrors();
    g->SetName(gname.Data());
  }

  int ip=0;
  gStyle->SetEndErrorSize(5);
  for (int n=1; n<=nMax; n++) {

    TGraphErrors* gALL = VnDVsQ(n, k, "RIDGE", "cA", opt);
    double xerr = opt.Contains("sys") ? 0. : 0.0; // 0.2
    g->SetPoint(ip, n, 100*gALL->GetY()[q] );
    g->SetPointError(ip, xerr,  100*gALL->GetEY()[q] );
    if (0)
      cout<<gALL->GetName()<< " " << g->GetY()[ip] << endl;

    ip++;
  }

  utils.set_tgraph_props(g, kBlack, kBlack, kFullCircle, 1.2); 
  g->SetLineWidth(2);
  gList->Add(g);
  return g;
}

TGraphErrors* GlobalvnVsN(int i, int k, int nMax, TString opt)
{
  // Create a TGraph from vn{GF} points
  // Stat errors by default. Call with opt "sys" to get sys errors.

  TString gname(Form("GlobalvnVsN_pt%.2gto%.2g_cent%dto%d",
		     ptlow[i],pthigh[i], centlow[k],centhigh[k]));
  if (!opt.IsNull())  {
    gname.Append("_");
    gname.Append(opt);
  }
  TGraphErrors* g = (TGraphErrors*) gList->FindObject(gname.Data());
  if (g) {
    return g;
  } 
  else {
    g = new TGraphErrors();
    g->SetName(gname.Data());
  }
  g->SetName(gname.Data());

  int ip=0;
  for (int n=1; n<=nMax; n++) {

    // point, stat err, sys err.
    double vnp = vnGF(n, k, i, 0, 999, "RIDGE", "cA", "");
    double vne = vnGF(n, k, i, 0, 999, "RIDGE", "cA", "err");
    double vns = vnGF(n, k, i, 0, 999, "RIDGE", "cA", "sys");
    double xerr = opt.Contains("sys") ? 0. : 0.0; // 0.2
    double yerr = opt.Contains("sys") ? vns : vne;

    g->SetPoint(ip, n, vnp );
    g->SetPointError(ip, xerr, yerr);
    ip++;
  }

  utils.set_tgraph_props(g, kBlack, kBlack, kFullCircle, 1.2); 
  g->SetLineWidth(2);

  return g;
}

TGraphErrors* VnDeltaNFVsN(int i, int j, int k, int nMax, TString opt)
{

  TGraphErrors* VnD = VnDeltaVsN(i,j,k,nMax,opt);
  TGraphErrors* vnt = GlobalvnVsN(i,k,nMax,opt);
  TGraphErrors* vna = GlobalvnVsN(j,k,nMax,opt);
  TString gname(Form("VnDeltaNFVsN_ptt%.2gto%.2g_pta%.2gto%.2g_cent%dto%d",
		     ptlow[i],pthigh[i],ptlow[j],pthigh[j],
		     centlow[k],centhigh[k]));
  if (!opt.IsNull())  {
    gname.Append("_");
    gname.Append(opt);
  }
  TGraphErrors* g = (TGraphErrors*) gList->FindObject(gname.Data());
  if (g) {
    return g;
  } 
  else {
    g = new TGraphErrors();
    g->SetName(gname.Data());
  }
  g->SetName(gname.Data());

  int ip=0;
  for (int n=1; n<=nMax; n++) {

    // Uncertainty vne is stat or sys depending on opt.
    double vnp = VnD->GetY()[ip] - 100*vnt->GetY()[ip]*vna->GetY()[ip];
    double vne = VnD->GetEY()[ip] - 100*vnt->GetEY()[ip]*vna->GetEY()[ip];
    g->SetPoint(ip, n, vnp );
    g->SetPointError(ip, 0.0, vne);
    ip++;
  }

  utils.set_tgraph_props(g, kBlack, kBlack, kOpenSquare, 1.2); 
  g->SetLineWidth(2);
  return g;
}


TF1* GlobalFitVsQ(int n, int k, int ipt1, int ipt2, const char* region, 
		  const char* corrtype, TString opt)
{

  // vn{GF}(trig) x vn{GF}(assc) vs. q, where the points between (and
  // including) ipt1 and ipt2 were used in the global fit.
  // Note: gfname must be identical format to fnName in GlobalFitVsQ().
  TF1* fitfn = 0;
  TString errType("meas");
  if (opt.Contains("hi_sys"))
    errType = "hi";
  if (opt.Contains("lo_sys"))
    errType = "lo";
  TString gfname = 
    Form("globalFitFn_cent%dto%d_n%d_ptbins%dto%d_%s_corrtype_%s_%s", 
	 centlow[k],centhigh[k], n, ipt1, ipt2, 
	 region, corrtype, errType.Data());
  fitfn = (TF1*)gFitList->FindObject(gfname.Data());
  if (fitfn) {
    if (0)
      Info("GlobalFitVsQ()", "Returning cached TF1 %s", gfname.Data());
    return fitfn;
  }

  // If nonexistent, create (& store) a new one
  fitfn = DoGlobalFit(n,k,ipt1,ipt2,region,corrtype,opt);
  //  fitfn = (TF1*)gFitList->FindObject(gfname.Data());
  
  if (!fitfn)
    Warning("GlobalFitVsQ()", 
	    "no fitfn %s,\neven after calling DoGlobalFit()", 
	    gfname.Data());

  return fitfn;
}

TF1* DoGlobalFit(int n, int k, int ipt1, int ipt2, 
		 const char* region, const char* corrtype, TString opt)
{
  // fnName must be identical format to gfname in GlobalFitVsQ().
  TF1 *globalFitFn = 0;
  TGraphErrors* g = 0;
  TString errType("meas");
  if (opt.Contains("hi_sys"))
    errType = "hi";
  if (opt.Contains("lo_sys"))
    errType = "lo";
  TString fnName = 
    Form("globalFitFn_cent%dto%d_n%d_ptbins%dto%d_%s_corrtype_%s_%s", 
	 centlow[k],centhigh[k], n, ipt1, ipt2, 
	 region, corrtype, errType.Data());

  g = VnDVsQ(n, k, region, corrtype, opt);

  if (!g) {
    Error("DoGlobalFit()", "No graph found");
    return 0;
  }

  // Assign global variables that will get checked in the vntvna
  // function during the fit.
  minFitPtBin = ipt1;
  maxFitPtBin = ipt2;

  // Set up the global fit function...  
  double qmax = g->GetN() - 0.4999;
  int nPars = maxpt + 3;
  globalFitFn = new TF1(fnName.Data(), vntvna, 0.0, qmax, nPars);
  
  /* Not needed now, I just flip the v1 points later...AMA 11/15/2011 */
  // Each parameter has an equally good minimum at + or - v_nbest. The
  // starting value is chosen to steer the fit to the positive "root".
  // v_1 is special because it is the only coefficient that crosses
  // zero at low pt, so it gets a negative starting point. Pars are
  // vn(ipt), then last one is n.
  //  double startVal = (n==1)? -0.02 : 0.01;

  double startVal = 0.01;
  for (int ip=0; ip<maxpt; ip++) { 
    globalFitFn->SetParameter(ip, startVal);
    globalFitFn->SetParName(ip, Form("v_%d(%.3g)",n,ptmean[ip]));
  }

  // Pass n into vntvna(). Access as par[maxpt].
  globalFitFn->FixParameter(maxpt, double(n));

  // Pass k into vntvna(). Access as par[maxpt+1].
  globalFitFn->FixParameter(maxpt+1, double(k));

  // Mom. conservation term prefactor. Access as par[maxpt+2].
  if (n==1) {
    globalFitFn->FixParameter(maxpt+2, 0.);
    //    globalFitFn->SetParLimits(maxpt+2, 0., 1e-3);
  }
  else
    globalFitFn->FixParameter(maxpt+2, 0.);
  
  globalFitFn->SetNpx(1000);

  TFitResultPtr result = g->Fit(globalFitFn,"QRNS"); // fast, takes ~15 ms

  if (0) {
    cout<<globalFitFn->GetName();
    result->Print();
  }
  // TMinuit status: 
  // migradResult + 10*minosResult + 100*hesseResult + 1000*improveResult
  // 0:ok, <0:user error, >0: see above.
  int fitStatus = result;
  if (fitStatus)
    Warning("DoGlobalFit()", "TMinuit status %d", fitStatus);

  gFitList->Add(globalFitFn);

  return globalFitFn;
}

double vnGF(int n, int k, int ptBin, int ipt1, int ipt2,
	    const char* region, const char* corrtype, TString opt)
{

  double vn, err;
  TF1* gfn = GlobalFitVsQ(n,k,ipt1,ipt2,region,corrtype,opt);
  if (!gfn) 
    Error("vnGF()", "!gfn");  
  vn = gfn->GetParameter(ptBin); // ptBin is the q index
  err = gfn->GetParError(ptBin);
  
  if (opt.Contains("sys")) {
    TF1* hi = GlobalFitVsQ(n,k,ipt1,ipt2,region,corrtype,"hi_sys");
    TF1* lo = GlobalFitVsQ(n,k,ipt1,ipt2,region,corrtype,"lo_sys");
    double vhi = hi->GetParameter(ptBin);
    double vlo = lo->GetParameter(ptBin);
    return 0.5*TMath::Abs(vhi-vlo);
  }

  // v1 can be negative, otherwise ensure the positive solution is chosen
  // if (n>1 && vn < 0)
  //   vn *= -1;

  // if (ptlow[ipt1] > 4.) // !!!!!!!!!!!!!!!!!!!!!!!!!
  //   vn = TMath::Abs(vn);

  // if (n==1) 
  //   vn += MomCorrection(i, j, k);

  return opt.Contains("err") ? err : vn;
}

TGraphErrors* vnGFvsPt(int n, int k, int ipt1, int ipt2, 
		       const char* region, const char* corrtype, TString opt)
{

  // vn{GF} vs. pt, where the points between (and including) ipt1 and
  // ipt2 were used in the global fit.
  bool isHighPtFit = ipt1 > 0 && ipt2==999;
  bool flipSign = n==1;  // Flip v1?
  if (isHighPtFit) 
    flipSign = 0;
  if (opt.Contains("keepsign"))
    flipSign = 0;
  TString name(Form("v%dGFvsPt_fitpt%dto%d_cent%dto%d",
		  n, ipt1, ipt2, centlow[k], centhigh[k]));
  if (opt.Contains("sys"))
    name.Append("_sys");
  if (opt.Contains("ptcons"))
    name.Append("_ptcons");

  TGraphErrors* g = (TGraphErrors*)gList->FindObject(name.Data());
  if (g)
    return g;
  
  g = new TGraphErrors();
  g->SetName(name);
  
  for (int i=0; i<maxpt; i++) {
    TString opt1 = "";
    if (opt.Contains("ptcons"))
      opt1.Append("_ptcons");
    double vn = vnGF(n,k,i,ipt1,ipt2,region,corrtype,opt1);
    
    if (flipSign)
      vn = -vn;
    
    //    if (isHighPtFit && i<ipt1)
    if (i<ipt1)
      vn = 0.;

    double mpt = MeanPt(i, 0, k, "t"); 
    if (i<ipt1) 
      mpt = -10.;
    g->SetPoint(i, mpt, vn);
    
    if (opt.Contains("sys")) {
      opt1.Append("_sys");
      double ex  = (gi->GetEX()[i] < 0.126) ? gi->GetEX()[i] : 0.250;
      double vns = vnGF(n,k,i,ipt1,ipt2,region,corrtype, opt1);
      g->SetPointError(i, ex, vns);
      //      g->SetPointError(i, gi->GetEX()[i], vns); // too wide 
    }
    else {
      opt1.Append("_err");
      double vne = vnGF(n,k,i,ipt1,ipt2,region,corrtype,opt1);
      g->SetPointError(i, 0.0, vne);
    } 
  }

  int col = CentColor(centlow[k], centhigh[k]);
  int col_pale = CentColorPale(centlow[k], centhigh[k]);
  g->SetLineWidth(2);
  g->SetFillColor(col_pale); // for sys
  utils.set_tgraph_props(g, col, col, kFullCircle, 1.2);

  if (0) cout<<opt.Data()<<endl;
  gList->Add(g);
  return g;
}

TGraph* Luzumv1(int cent, const char* pid) 
{
  // v1 for pions vs pt from viscous hydro calc by Matt Luzum
  // valid cent args 0-4 - 10% bins.
  // pid = "pi", "K", "p".
  const char* dat = Form("v1LHCGlauber%ssv08.dat", pid);
  TString fmt("%lg "); // pt
  for (int ky=0; ky<5; ky++) 
    fmt.Append( ky==cent? "%lg " : "%*lg "); // centrality - skip if *lg
  return new TGraph(dat, fmt.Data());
}

TGraphErrors* GlobalFitVsPtAssoc(int i, int k, int n, TString opt)
{
  // Create a TGraph of vnt{gf} x vna{gf} vs. pt_assoc from gfit
  TF1* gfit = GlobalFitVsQ(n, k, 0, 999, "RIDGE", cfDef, opt);
  TGraphErrors* g = new TGraphErrors();
  TString gname(Form("GlobalFitVsPtAssoc_ptt%d_cent%d_n%d",i,k,n));
  if (!opt.IsNull())
    gname.Append(opt);
  g->SetName(gname);

  for (int j=0; j<=i; j++) {
    int q = GlobalIndex(i, j);
    g->SetPoint(j, gj->GetX()[j], gfit->Eval(q + 0.5) );
    g->SetPointError(j, 0, 0 );
    utils.set_tgraph_props(g, colorsPale[i], colorsPale[i], kFullCircle, 1.0); 
    g->SetLineWidth(4);
  }
  return g;
}

double SineError(int i, int j, int k)
{
  double rms = 0;
  double resid_sin=0, rs_err=0;
  TH1* hst  = Hist("RIDGE", "cA", i, j, k, ""); 	
  double arr[100] = {0};
  //  arr[0] = 0;
  for (int n=1; n<=gNMax; n++) {
    VnDelta(n, *hst,  resid_sin, rs_err, "sin");
    arr[n-1] = resid_sin;
  }

  rms = TMath::RMS(gNMax, arr) / 2;
  if (0)
    cout << " SineError: " << rms << endl;

  return rms;
}

TCanvas* DrawVnDeltaVsN(int nptt, int pttbins[], int npta, 
			int ptabins[], int ncb, int centbins[], TString opt)
{

  double t1=ptlow[pttbins[0]], t2=pthigh[pttbins[nptt-1]];
  double a1=ptlow[ptabins[0]], a2=pthigh[ptabins[npta-1]];
  TString cname = Form("VnDeltaVsN_ptt%.2gto%.2g_pta%.2gto%.2g",t1,t2,a1,a2);
  if (!opt.IsNull()) {
    cname.Append("_");
    cname.Append(opt);
  }
  TCanvas* c = new TCanvas(cname.Data(), cname.Data(), 650,700);
  double legbottom = 0.9 - 0.07*ncb;
  TLegend* leg = new TLegend(0.62, legbottom, 0.98, 0.98, "Centrality");
  leg->SetFillColor(kNone);
  leg->SetBorderSize(0);
  leg->SetName("leg");

  TObjArray* graphs = new TObjArray();
  double xmin = 0.4;
  double xmax = opt.Contains("ext") ? gNMax + 0.5 : 8.5;
  int nXticks =  10, nYticks = 210;
  TLine l;
  TH1F* hf = gPad->DrawFrame(xmin, 0.0001, xmax, 0.26);
  //  bool scale1000 = 0;

  hf->SetTitle(Form(";n;V_{n#Delta} [10^{-2}]"));
  if (opt.Contains("sine"))
    hf->SetTitle(Form(";n;#LTsin(n#Delta#phi)#GT [10^{-2}]"));

  // if (opt.Contains("ext") )  
  //   if (!opt.Contains("sine")) {
  //     scale1000 = 0;
  //     hf->SetTitle(Form(";n;V_{n#Delta}"));
  // }

  l.SetLineStyle(kDashed);

  utils.make_nice_axes(c, hf, 1.3, nXticks, nYticks, 1.4,-1.4);
  hf->SetLabelOffset(0.0, "X");
  hf->SetTitleOffset(-1.96, "Y");

  hf->GetXaxis()->SetTicks("+-");
  hf->GetYaxis()->SetTicks("+-");
  l.DrawLine(hf->GetXaxis()->GetXmin(), 0, hf->GetXaxis()->GetXmax(),0);
  gPad->SetLeftMargin(0.22);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);

  for (int i=0; i<nptt; i++) {
    for (int j=0; j<npta; j++) {
      for(int k=0; k<ncb; k++) {

	// Stat errors on g, sys errors on gs.
	TString opt1 = opt.Contains("sine") ? "sine" : "";
	TGraphErrors* g = 
	  VnDeltaVsN(pttbins[i], ptabins[j], centbins[k], gNMax, opt1);
	TGraphErrors* gs = 
	  VnDeltaVsN(pttbins[i], ptabins[j], centbins[k], gNMax, "sys");
	TGraphErrors* gc = (TGraphErrors*)g->Clone();
	int dark = CentColor(centlow[centbins[k]], centhigh[centbins[k]]);
	int pale = CentColorPale(centlow[centbins[k]], centhigh[centbins[k]]);
	utils.set_tgraph_props(g, dark, dark, kFullCircle, 1.5); 
	utils.set_tgraph_props(gs, dark, dark, kDot, 2.); 
	utils.set_tgraph_props(gc, kBlack, kBlack, kOpenCircle, 1.5);
	g->SetFillColor(pale);
	gs->SetFillColor(dark);

	TGraphErrors* gnf = 0;
	TGraphErrors* gnf_sys = 0;
	TGraphErrors* gnf_solid = 0;
	if (opt.Contains("vnf")) {  // Non-factorizing part
	  gnf = 
	    VnDeltaNFVsN(pttbins[i], ptabins[j], centbins[k], gNMax, "");
	  gnf_sys = 
	    VnDeltaNFVsN(pttbins[i], ptabins[j], centbins[k], gNMax, "sys");

	  // Just for appearance...
	  gnf_solid = (TGraphErrors*) gnf->Clone(Form("%s_solid", gnf->GetName()));
	  utils.set_tgraph_props(gnf_solid, pale, pale, kFullSquare, gnf->GetMarkerSize());
	}
	// if (scale1000) {
	//   for (int ip=0; ip<g->GetN(); ip++) {
	//     g->SetPoint(ip,   g->GetX()[ip], 1000* g->GetY()[ip]);
	//     gs->SetPoint(ip, gs->GetX()[ip], 1000*gs->GetY()[ip]);
	//     gc->SetPoint(ip, gc->GetX()[ip], 1000*gc->GetY()[ip]);

	//     g->SetPointError(ip,   g->GetEX()[ip], 1000* g->GetEY()[ip]);
	//     gs->SetPointError(ip, gs->GetEX()[ip], 1000*gs->GetEY()[ip]);
	//     gc->SetPointError(ip, gc->GetEX()[ip], 1000*gc->GetEY()[ip]);
	//   }
	// }

	if (opt.Contains("sine")) {
	  g->Fit("pol0", "Q");
	  TF1* gfit = g->GetFunction("pol0");
	  gfit->SetLineColor(dark);
	}

	else if (!opt.Contains("nobars"))
	  g->Draw("Bsame");  // bars

	if (opt.Contains("vnf")) {  // Non-factorizing part
	  utils.draw_errorbox(gnf_sys, pale, 0.3);
	  gnf_solid->Draw("epsame");
	  gnf->Draw("epsame");
	}

	g->Draw("psame");  // colored markers
	gc->Draw("epsame");  // open markers w/black stat error bars

	if (!opt.Contains("sine"))
	  utils.draw_errorbox(gs, dark, 0.3);

	graphs->Add(g);
	leg->AddEntry(g,Form("%d-%d%%",centlow[centbins[k]],centhigh[centbins[k]]),"l,p");

	if (opt.Contains("vnf")) {  // Non-factorizing part
	  leg->AddEntry(gnf_solid,"NF V_{n#Delta}","l,p");
	  //	  leg->AddEntry(gnf_solid,"[v_{n}(p_{T}^{t}) #times v_{n}(p_{T}^{a})]_{GF}","l,p");
	  leg->SetHeader(""); 
	}

      }
    }
  }

  // Adjust y limits
  double newYmin = 0, newYmax = 0, newSpace = 0;
  for (int m=0; m<graphs->GetEntries(); m++) {
    TGraphErrors* tg = (TGraphErrors*)graphs->At(m);
    double yhi = TMath::MaxElement(tg->GetN()-1, tg->GetY());
    double ylo = TMath::MinElement(tg->GetN()-1, tg->GetY());
    if (yhi > newYmax) newYmax = yhi;
    if (ylo < newYmin) newYmin = ylo;
  }

  newSpace = 0.1*(newYmax - newYmin);

  if (opt.Contains("ext"))
    hf->GetYaxis()->SetRangeUser(newYmin-newSpace, newYmax+4*newSpace);
  else
    hf->GetYaxis()->SetRangeUser(newYmin-newSpace, newYmax+newSpace);

  // if (opt.Contains("sine")) {
  //   hf->GetYaxis()->SetRangeUser(0.0001, newYmax+4*newSpace);
  //   gPad->SetLogy();
  // }

  if (!opt.Contains("sine")) 
    leg->Draw();
 
  // Draw text
  ltx.SetTextSize(0.04); // 0.05

  if ( opt.Contains("sine") )  
    legbottom = 0.99;
  else if (opt.Contains("ext") && !opt.Contains("vnf") )  
    legbottom -= 0.4;
  
  double xstart = 0.6;
  ltx.DrawLatex(xstart, legbottom-1*0.07, Form("%.2g < p_{T}^{t} < %.2g GeV/c", t1, t2));
  ltx.DrawLatex(xstart, legbottom-2*0.07, Form("%.2g < p_{T}^{a} < %.2g GeV/c", a1, a2));
  
  gPad->Update();
  gPad->Modified();
  return c; 
}

TCanvas* DrawGlobalvnVsN(int npt, int ptbins[], int ncb, int centbins[], TString opt)
{
  double t1=ptlow[ptbins[0]], t2=pthigh[ptbins[npt-1]];
  TString cname = Form("GlobalvnVsN_pt%.2gto%.2g",t1,t2);
  if (!opt.IsNull()) {
    cname.Append("_");
    cname.Append(opt);
  }
  TCanvas* c = new TCanvas(cname.Data(), cname.Data(), 600,700);
  double legbottom = 0.9 - 0.07*ncb;
  TLegend* leg = new TLegend(0.62, legbottom, 0.98, 0.98, "Centrality");
  leg->SetFillColor(kNone);
  leg->SetBorderSize(0);
  leg->SetName("leg");

  TObjArray* graphs = new TObjArray();
  double xmin = opt.Contains("zoom") ? 2.6 : 0.4;
  double xmax = 8.5;
  int nXticks =  10, nYticks = 210;
  TLine l;
  TH1F* hf = gPad->DrawFrame(xmin, -0.06, xmax, 0.26);
  hf->SetTitle(Form(";n;v_{n}{GF}"));
  l.SetLineStyle(kDashed);

  utils.make_nice_axes(c, hf, 1.3, nXticks, nYticks, 1.4,-1.4);
  hf->SetLabelOffset(0.0, "X");
  if (opt.Contains("zoom"))
    hf->SetTitleOffset(-2.3, "Y");
  else
    hf->SetTitleOffset(-1.9, "Y");

  hf->GetXaxis()->SetTicks("+-");
  hf->GetYaxis()->SetTicks("+-");
  l.DrawLine(hf->GetXaxis()->GetXmin(), 0, hf->GetXaxis()->GetXmax(),0);
  if (opt.Contains("zoom"))
    gPad->SetLeftMargin(0.25);
  else
    gPad->SetLeftMargin(0.2);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);

  for (int i=0; i<npt; i++) {
    for(int k=0; k<ncb; k++) {

      // Stat errors on g, sys errors on gs.
      TGraphErrors* g  = GlobalvnVsN(ptbins[i], centbins[k], gNMax, "");
      TGraphErrors* gs = GlobalvnVsN(ptbins[i], centbins[k], gNMax, "sys");
      TGraphErrors* gc = (TGraphErrors*)g->Clone();
      int dark = CentColor(centlow[centbins[k]], centhigh[centbins[k]]);
      int pale = CentColorPale(centlow[centbins[k]], centhigh[centbins[k]]);
      utils.set_tgraph_props(g, dark, dark, kFullCircle, 1.5); 
      //      utils.set_tgraph_props(gs, dark, dark, kDot, 2.); 
      utils.set_tgraph_props(gc, kBlack, kBlack, kOpenCircle, 1.5);
      g->SetFillColor(pale);
      //      gs->SetFillColor(dark);
      g->Draw("Bsame");  // colored bars
      utils.draw_errorbox(gs, dark, 0.3);
      g->Draw("psame");  // colored markers
      //      gs->Draw("e[]psame"); // sys error filled rectangles
      gc->Draw("epsame");  // open markers w/black stat error bars
      graphs->Add(g);
      leg->AddEntry(g,Form("%d-%d%%",centlow[centbins[k]],centhigh[centbins[k]]),"l,p");
    }
  }

  // Adjust y limits
  double newYmin = 0, newYmax = 0, newSpace = 0;
  for (int m=0; m<graphs->GetEntries(); m++) {
    TGraphErrors* tg = (TGraphErrors*)graphs->At(m);
    double yhi = TMath::MaxElement(tg->GetN()-1, tg->GetY());
    double ylo = TMath::MinElement(tg->GetN()-1, tg->GetY());
    if (yhi > newYmax) newYmax = yhi;
    if (ylo < newYmin) newYmin = ylo;
  }

  newSpace = 0.1*(newYmax - newYmin);
  hf->GetYaxis()->SetRangeUser(newYmin-newSpace, newYmax+newSpace);
  if (opt.Contains("zoom"))
    hf->GetYaxis()->SetRangeUser(-0.0008, 0.0085);

  leg->Draw();
  c->Update();
 
  // Draw text
  // if (opt.Contains("zoom"))
  //   ltx.SetTextSize(0.03); // 0.04
  // else

  // ltx.SetTextSize(0.04); // 0.05

  double xstart = (opt.Contains("zoom")) ? 0.7 : 0.6;
  ltx.DrawLatex(xstart, legbottom-1*0.07, Form("%.2g < p_{T} < %.2g GeV/c", t1, t2));
  // ltx.DrawLatex(xstart, legbottom-2*0.07, Form("%.2g < p_{T}^{assoc} < %.2g GeV/c", a1, a2));

  return c; 
}


TCanvas* DrawVnFromGlobalFit(int n, int ipt1, int ipt2, int ncb, 
			     int centbins[], TString opt)
{
  // Draw vn for one n on one canvas - many centralities
  initialize();
  double maxPtFit = ipt2==999 ? pthigh[maxpt-1] : pthigh[ipt2];
  const char* cname = Form("v%d_global%s_ptfit%.2gto%.2g",
			   n, opt.Data(), ptlow[ipt1], maxPtFit);
  const char* title = Form("v%d_global%s_ptfit%.2gto%.2g",
			   n, opt.Data(), ptlow[ipt1], maxPtFit);
  double lvleft = opt.Contains("multi")? 0.74 : 0.82;
  double lvbot = opt.Contains("multi")? 0.25 : 0.56;
  int cvHeight = opt.Contains("multi")? 500 : 500;
  TCanvas* cv = new TCanvas(cname, title, 700, cvHeight);
  TLegend* lv = new TLegend(lvleft, lvbot, 0.99, 0.97);
  lv->SetFillColor(kNone); lv->SetName("lv");
  lv->SetMargin(0.3);
  lv->SetBorderSize(1);
  lv->SetHeader("#splitline{Fit range}{(GeV/c):}");
  utils.padsetup(cv, 1, 1, "", 0.12, 0.01, 0.03, 0.2); // ######
  //  gPad->SetRightMargin(0.15);


  TPad* inset = 0;

  TObjArray* graphs = new TObjArray();
  double xmin = 0., xmax = 14.8; // 8.3;
  double ymin = -0.06, ymax = 0.26;
  if (opt.Contains("multi")) {
    ymin = -1.19;
    ymax = 0.52;
  }
  TH1F* hf = gPad->DrawFrame(xmin, ymin, xmax, ymax);
  hf->SetTitle(Form(";p_{T}^{t} (GeV/c);v_{%d}{GF}",n));
  int nXticks =  210;
  int nYticks = 210;
  TLine l;
  l.SetLineStyle(kDashed);
  if (n==1) {
    cv->cd();
    inset = new TPad("inset", "inset", 0.18, 0.18, 0.45, 0.59, kNone, 1, 0);
    inset->Draw();
    inset->cd();
    TH1F* hfi = gPad->DrawFrame(-0.05, -0.04, 4.2, 0.16);
    TAxis* hx = hfi->GetXaxis();     
    TAxis* hy = hfi->GetYaxis();
    hx->SetLabelSize(0.1);
    hy->SetLabelSize(0.1);
    hx->SetNdivisions(104);
    hy->SetNdivisions(104);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.2);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.01);
    //    utils.make_nice_axes(inset, hfi, 3, nXticks, nYticks, 1.4,-1.4);
    l.DrawLine(-0.05, 0, 4.2, 0);
    cv->cd();
  }

  utils.make_nice_axes(cv, hf, 1.3, nXticks, nYticks, 1.4,-1.4);
  hf->SetLabelOffset(0.0, "X");
  hf->SetTitleOffset(1.5, "Y");
  hf->GetXaxis()->SetTicks("+");
  //  hf->GetYaxis()->SetTicks("+-");

  TGraph* v1h[5][3];  
  if (opt.Contains("luzum")) {
    int lSty[] = {7, kSolid, kSolid}; 
    int lCol[] = {kBlack, kCyan, kBlue, kGreen+1, kRed}; 
    TString pidString[] = {"pi", "K", "p"};
    for (int m=0; m<5; m++) {
      for (int pid=0; pid<3; pid++) {
	v1h[m][pid] = Luzumv1(m, pidString[pid].Data()); 	
	TGraph* g1 = v1h[m][pid];
	g1->SetLineStyle(lSty[pid]);
	g1->SetLineColor(lCol[m]);
	g1->SetLineWidth(2);
      }
    }
    hf->GetXaxis()->SetRangeUser(0, 4.8);
    l.DrawLine(0, 0, 4.8, 0);
    v1h[0][0]->Draw("plsame");
    v1h[0][2]->Draw("plsame");
    v1h[4][0]->Draw("plsame");
    v1h[4][2]->Draw("plsame");
  }
  else
    l.DrawLine(hf->GetXaxis()->GetXmin(), 0, hf->GetXaxis()->GetXmax(),0);

  for (int cb=0; cb<ncb; cb++) {
    int k = centbins[cb];

    if (!opt.Contains("multi")) {    
      // last 2 args are included pt bins
      TGraphErrors* gtmp = vnGFvsPt(n, k, ipt1, ipt2);

      // clone to modify without ripples
      TGraphErrors* gc = (TGraphErrors*) gtmp->Clone();
      if (!gc) {
	Error("DrawVnFromGlobalFit()", "Problem making graph" );
	continue;
      }
      TGraphErrors* gtmp_sys = vnGFvsPt(n, k, ipt1, ipt2, "RIDGE", cfDef, "sys");
      TGraphErrors* gs = (TGraphErrors*) gtmp_sys->Clone();

      if (opt.Contains("luzum")) {
	if (centlow[k]==20)
	  utils.set_tgraph_props(gc, kBlue, kBlue, kFullCircle, 1.2);
	if (centlow[k]==40)
	  utils.set_tgraph_props(gc, kRed, kRed, kFullCircle, 1.2);
      }

      utils.draw_errorbox(gs, gs->GetFillColor());
      gc->Draw("epsame");
      graphs->Add(gc);

      TGraphErrors* open = (TGraphErrors*) gc->Clone(Form("%s_open", gc->GetName()));
      utils.set_tgraph_props(open, kBlack, kBlack, kOpenCircle, 1.2);
      open->Draw("epsame");

      lv->AddEntry(gc,Form("%d-%d%%",centlow[k],centhigh[k]),"l,p");
    }
    // NEW
    if (opt.Contains("multi")) {
      int b1,b2,b3,b4;
      b1=PtBin(0.25,0.5);
      b2=PtBin(0.75,1.0);
      b3=PtBin(2.0, 2.5);
      b4=PtBin(3.0, 4.0);
      TGraphErrors* gg1 = vnGFvsPt(n, k, b1, b2, "RIDGE", cfDef, "");
      TGraphErrors* gg2 = vnGFvsPt(n, k, b3, b4, "RIDGE", cfDef, "keepsign");
      TGraphErrors* gs1 = vnGFvsPt(n, k, b1, b2, "RIDGE", cfDef, "sys");
      TGraphErrors* gs2 = vnGFvsPt(n, k, b3, b4, "RIDGE", cfDef, "sys_keepsign");
      utils.set_tgraph_props(gg1, gg1->GetLineColor(), gg1->GetMarkerColor(), kOpenCircle, 1.0);
      utils.set_tgraph_props(gg2, gg2->GetLineColor(), gg2->GetMarkerColor(), kFullSquare, 1.0);
      
      utils.draw_errorbox(gs1, gs1->GetFillColor());
      utils.draw_errorbox(gs2, gs2->GetFillColor());
      gg1->Draw("epsame");
      gg2->Draw("epsame");
      graphs->Add(gs1);
      graphs->Add(gs2);
      const char* s1 = Form("#splitline{%.2g < p_{T}^{a} < %.2g  }{%d-%d%%}",
			    ptlow[b1], pthigh[b2], centlow[k],centhigh[k]);
      const char* s2 = Form("#splitline{%.2g < p_{T}^{a} < %.2g  }{%d-%d%%}",
			    ptlow[b3], pthigh[b4], centlow[k],centhigh[k]);
      lv->AddEntry(gg1, s1, "l,p");
      lv->AddEntry(gg2, s2, "l,p");

      
      if (n==1) {
	inset->cd();
	gg1->Draw("epsame");
	gg2->Draw("epsame");
	cv->cd();
      }
      
    }

  } // centrality bin loop

  if (opt.Contains("luzum")) {
    lv->AddEntry(v1h[0][0], Form("v.h. (#pi)"),"l,p");
    lv->AddEntry(v1h[0][2], Form("v.h. (p)"), "l,p");
  }
  lv->Draw();
  ltx.SetTextSize(0.05);

  if (opt.Contains("multi")) {
    gPad->SetTopMargin(0.02);
    if (ispp)
    ltx.DrawLatex(0.19, 0.92, Form("%s", "#splitline{pp}{2.76 TeV}"));
    else
    ltx.DrawLatex(0.19, 0.92, Form("%s", "#splitline{Pb-Pb}{2.76 TeV}"));
    ltx.DrawLatex(0.49, 0.94, Form("%.2g < |#Delta#eta| < 1.8", minRidgeDEta) ); 

  }
  else {
    ltx.DrawLatex(0.17, 0.88, Form("%.2g < |#Delta#eta| < 1.8", minRidgeDEta) ); 
    ltx.DrawLatex(0.17, 0.94, Form("%.2g < fit p_{T}^{a} < %.2g GeV/c",
				   ptlow[ipt1], maxPtFit));
    ltx.DrawLatex(0.19, 0.78, Form("%s", "#splitline{Pb-Pb}{2.76 TeV}"));
  }
  double newYmin = 0, newYmax = 0, newSpace = 0;
  for (int m=0; m<graphs->GetEntries(); m++) {
    TGraphErrors* tg = (TGraphErrors*)graphs->At(m);
    double yhi = TMath::MaxElement(tg->GetN()-1, tg->GetY());
    double ylo = TMath::MinElement(tg->GetN()-1, tg->GetY());
    if (yhi > newYmax) newYmax = yhi;
    if (ylo < newYmin) newYmin = ylo;
  }

  if (opt.Contains("multi")) {
    if(n==1) 
      newYmin *= 2.2;
    if(n==2) 
      newYmin = 0.02;
  }

  newSpace = 0.2*(newYmax - newYmin);
  //  if (opt.Contains("multi")) newSpace /= 2;
  if (n==5) newSpace *= 2;
  hf->GetYaxis()->SetRangeUser(newYmin-newSpace, newYmax+newSpace);
  cv->Modified();
  cv->Update();
  cv->Draw();
  return cv;
}

TCanvas* Drawv1to5(int ncb, int centbins[], int ipt1, int ipt2, TString opt)
{
  if (!opt.IsNull())
    cout<<opt.Data()<<endl;

  initialize();
  int cent1 = centlow[0];
  int cent2 = centhigh[centbins[ncb-1]];
  double pmax = 15.;
  if (ipt2 < 999.)
    pmax = pthigh[ipt2];
  const char* cname = Form("vn_etamin%02d_cent%dto%d_fitptbin%dto%d", 
			   (int)(10*minRidgeDEta), cent1, cent2, ipt1, ipt2);
  const char* title = Form("v1to5gf_etamin%02d_cent%dto%d_fitpt%.2gto%.2g", 
			   (int)(10*minRidgeDEta), cent1, cent2, ptlow[ipt1], pmax);
  TCanvas* cv = new TCanvas(cname, title, 1400, 500);
  TLegend* lv = new TLegend(0.09, 0.72, 0.27, 0.88);
  lv->SetFillColor(kNone);
  lv->SetBorderSize(0);
  TLine zero;
  double lv2bottom = 0.62;
  TLegend* lv2 = new TLegend(0.82, lv2bottom, 0.93, 0.95, "Centrality");
  lv2->SetFillColor(kNone);
  lv2->SetBorderSize(0);
  utils.padsetup(cv, 5, 1, "", 0.12, 0.01, 0.03, 0.2);

  for (int n=1; n<=5; n++) {
    Printf("Drawv1to5() n = %d / 5", n);
    cv->cd(n);
    double xmin = -0.43, xmax =  8.3;
    double ymin = -0.09, ymax =  0.26;
    if (ptlow[ipt1] > 4.) {
      xmin = 4.9;
      ymin = -0.1;
    }
    if (pthigh[ipt2] > 6.) {
      xmax = 11.8;
      ymax = 1;
    }
    if (opt.Contains("split")) {
      xmin = -0.43;
      xmax = 11.8;
      ymin = -0.2;
      ymax = 1;
    }

    TH1F* hf = gPad->DrawFrame(xmin, ymin, xmax, ymax);
    int nXticks = 206;
    int nYticks = 208;
    if (n==1) {
      utils.make_nice_axes(cv, hf, 2., nXticks, nYticks, -0.02, 0.02);
      hf->SetLabelOffset(-0.01, "X");
    }
    else {
      utils.make_nice_axes(cv, hf, 3.4, nXticks, nYticks, -0.075, 0.0);
      hf->SetLabelOffset(-0.057, "X");
      hf->SetTickLength(0.02, "X");
      hf->SetTickLength(0.05, "Y");
    }
    hf->GetXaxis()->SetTicks("+-");
    hf->GetYaxis()->SetTicks("+-");
    zero.SetLineStyle(kDashed);
    zero.DrawLine(xmin, 0., xmax, 0.);

    // Loop thru every cent bin. Skip non-requested ones.
    // for (int k=0; k<maxcent; k++) {
    //   bool centOK = 1;
    //   if (ncb != 999) {
    // 	centOK = 0;
    // 	for (int cb=0; cb<ncb; cb++) 
    // 	  if(k==centbins[cb])
    // 	    centOK = 1;
    //   }
    //   if (!centOK)
    // 	continue;

    for (int cb=0; cb<ncb; cb++) {
      int k = centbins[cb];
      TGraphErrors* gc = vnGFvsPt(n, k, ipt1, ipt2);
      TGraphErrors* gs = vnGFvsPt(n, k, ipt1, ipt2, "RIDGE", cfDef, "sys");

      if (!gc) {
	Error("Drawv1to5()", "Problem finding graph");
	continue;
      }

      TGraphErrors* gc_hipt=0;
      TGraphErrors* gs_hipt=0;
      if (opt.Contains("highptfit") || opt.Contains("split")) {
	gc_hipt = vnGFvsPt(n, k, ipt2+1, 999, "highptfit");
	gs_hipt = vnGFvsPt(n, k, ipt2+1, 999, "RIDGE", cfDef, "sys");
	gc_hipt->SetMarkerStyle(kOpenCircle);
      }

      int pale = gs->GetFillColor(); // it is CentColor(centlow[k], centhigh[k]).
      utils.draw_errorbox(gs, pale);
      gc->Draw("epsame");

      if (opt.Contains("split")) {
	utils.draw_errorbox(gs_hipt, pale);
	gc_hipt->Draw("epsame");
      }
      
      TGraphErrors* open = (TGraphErrors*) gc->Clone(Form("%s_open", gc->GetName()));
      utils.set_tgraph_props(open, kBlack, kBlack, kOpenCircle, 1.2);
      open->Draw("epsame");

      if (n==1)
	lv2->AddEntry(gc,Form("%d-%d%%",centlow[k],centhigh[k]),"l,p");

    } // cent bin k
  } // Fourier moment n
  
  cv->cd();
  TPad* overlay = new TPad("overlay", "overlay", 0,0,1,1);
  overlay->SetFillStyle(4000); // transparent
  overlay->Draw();
  overlay->cd();
  
  ltx.SetTextSize(0.07);
  lv2->Draw();
  
  for (int n=1; n<=5; n++) {
    cv->cd(n);
    ltx.SetTextSize(n==1? 0.1 : 0.15);
    ltx.DrawLatex((n==1)? 0.85 : 0.7, 0.86, Form("v_{%d}", n));

    if (0 && n==1)
      AddPrelimStampToCurrentPad(0.45,0.7,0.75,0.95, "");
  }
  overlay->cd();
  ltx.SetTextSize(0.075);
  ltx.DrawLatex(0.5, 0.05, "p_{T} [GeV/c]");
  ltx.SetTextAngle(90);
  ltx.DrawLatex(0.05, 0.5, "v_{n }(p_{T})");
  ltx.SetTextAngle(0);
  ltx.SetTextSize(0.05);
  ltx.DrawLatex(0.14, 0.88, Form("%s", "#splitline{Pb-Pb}{2.76 TeV}"));
  ltx.SetTextSize(0.05);
  ltx.DrawLatex(0.13, 0.24, Form("%.2g < |#Delta#eta| < 1.8", minRidgeDEta) ); 
  
  const char* frange = Form("%.2g < fit p_{T}^{a} < %.2g GeV",
  			    ptlow[ipt1], pmax);  
  const char* frange1 = Form("#splitline{solid:}{%.2g < fit p_{T}^{a} < %.2g GeV}",
			    ptlow[ipt1], pthigh[ipt2]);  
  const char* frange2 = Form("#splitline{open:}{%.2g < fit p_{T}^{a} < %.2g GeV}",
			    ptlow[ipt2+1], pthigh[maxpt-1]);  

  if (opt.Contains("split")) {
    ltx.SetTextSize(0.04);
    ltx.DrawLatex(0.31, 0.8, frange1);
    ltx.DrawLatex(0.31, 0.7, frange2);
  }
  else
    ltx.DrawLatex(0.3, 0.27, frange);

  return cv;
}

TCanvas* Drawv2to5(int ncb, int centbins[], int ipt1, int ipt2, TString opt)
{
  if (!opt.IsNull())
    cout<<opt.Data()<<endl;

  initialize();
  int cent1 = centlow[0];
  int cent2 = centhigh[centbins[ncb-1]];
  double pmax = 15.;
  if (ipt2 < 999.)
    pmax = pthigh[ipt2];
  const char* cname = Form("v2to5_etamin%02d_cent%dto%d_fitptbin%dto%d", 
			   (int)(10*minRidgeDEta), cent1, cent2, ipt1, ipt2);
  const char* title = Form("v2to5gf_etamin%02d_cent%dto%d_fitpt%.2gto%.2g", 
			   (int)(10*minRidgeDEta), cent1, cent2, ptlow[ipt1], pmax);
  TCanvas* cv = new TCanvas(cname, title, 1400, 500);
  TLine zero;
  double lv2bottom = 0.58; // 0.5;
  TLegend* lv2 = new TLegend(0.77, lv2bottom, 0.9, 0.95, "Centrality");
  lv2->SetFillColor(kNone);
  lv2->SetBorderSize(0);
  utils.padsetup(cv, 4, 1, "", 0.08, 0.01, 0.03, 0.2);

  for (int n=2; n<=5; n++) {
    Printf("Drawv2to5() n = %d", n);
    cv->cd(n-1);
    double xmin = -0.2;
    double xmax =  7.6;
    TH1F* hf = gPad->DrawFrame(xmin, -0.12, xmax, 0.26); //was y = -0.03, 0.26
    int nXticks = 206;
    int nYticks = 208;
    if (n==2) {
      utils.make_nice_axes(cv, hf, 2., nXticks, nYticks, -0.02, 0.02);
      hf->SetLabelOffset(0.0, "X"); //       hf->SetLabelOffset(-0.01, "X");
    }
    else {
      utils.make_nice_axes(cv, hf, 2.6, nXticks, nYticks, -0.075, 0.0);
      hf->SetLabelOffset(-0.02, "X");
      hf->SetTickLength(0.02, "X");
      hf->SetTickLength(0.05, "Y");
    }
    hf->GetXaxis()->SetTicks("+-");
    hf->GetYaxis()->SetTicks("+-");
    zero.SetLineStyle(kDashed);
    zero.DrawLine(xmin, 0., xmax, 0.);

    for (int cb=0; cb<ncb; cb++) {
      int k = centbins[cb];
      
      TGraphErrors* gc = vnGFvsPt(n, k, ipt1, ipt2);
      TGraphErrors* gs = vnGFvsPt(n, k, ipt1, ipt2, "RIDGE", cfDef, "sys");

      if (!gc) {
	Error("Drawv2to5()", "Problem finding graph");
	continue;
      }

      int pale = gs->GetFillColor(); // it is CentColor(centlow[k], centhigh[k]).
      utils.draw_errorbox(gs, pale);
      //      gs->Draw("e2psame");
      gc->Draw("epsame");
      TGraphErrors* open = (TGraphErrors*) gc->Clone(Form("%s_open", gc->GetName()));
      utils.set_tgraph_props(open, kBlack, kBlack, kOpenCircle, 1.2);
      open->Draw("epsame");

      if (n==2)
	lv2->AddEntry(gc,Form("%d-%d%%",centlow[k],centhigh[k]),"l,p");

    } // cent bin k
  } // Fourier moment n
  
  cv->cd();
  TPad* overlay = new TPad("overlay", "overlay", 0,0,1,1);
  overlay->SetFillStyle(4000); // transparent
  overlay->Draw();
  overlay->cd();
  
  ltx.SetTextSize(0.07);
  lv2->Draw();
  
  for (int n=2; n<=5; n++) {
    cv->cd(n-1);
    ltx.SetTextSize(n==2? 0.1 : 0.15);
    ltx.DrawLatex((n==2)? 0.85 : 0.7, 0.86, Form("v_{%d}", n));
  }
  overlay->cd();
  ltx.SetTextSize(0.075);
  ltx.DrawLatex(0.48, 0.04, "p_{T} [GeV/c]");
  ltx.SetTextAngle(90);
  ltx.DrawLatex(0.02, 0.5, "v_{n }{GF}");
  ltx.SetTextAngle(0);
  ltx.SetTextSize(0.05);
  ltx.DrawLatex(0.33, 0.88, Form("%s", "#splitline{Pb-Pb}{2.76 TeV}"));
  ltx.SetTextSize(0.05);
  ltx.DrawLatex(0.33, 0.75, Form("%.2g < |#Delta#eta| < 1.8", minRidgeDEta) ); 
  
  const char* frange = Form("#splitline{Global fit range}{%.2g < p_{T} < %.2g GeV}", 
			    ptlow[ipt1], pmax);  
  ltx.DrawLatex(0.56, 0.75, frange);
  
  return cv;
}

TCanvas* DrawChi2vsPtTrig(int npta, int ptabins[], int ncb, int centbins[], TString opt)
{
    
  // TODO: options-- global or 1to5, chi2/N or prob
  double a1=ptlow[ptabins[0]], a2=pthigh[ptabins[npta-1]];
  int c1=centlow[centbins[0]], c2=centhigh[centbins[ncb-1]];
  TString cname = Form("chi2_etamin%02d_pta%.2gto%.2g_cent%dto%d",
		       (int)(10*minRidgeDEta), a1, a2, c1, c2 );
  
  if (!opt.IsNull() ) {
    cname.Append("_");
    cname.Append(opt.Data());
  }

  TCanvas* cv = new TCanvas(cname, cname, 1);
  TLegend* lv = new TLegend(0.09, 0.72, 0.27, 0.88); // w.r.t. cv
  lv->SetFillColor(kNone);
  lv->SetBorderSize(0);

  double lv2bottom = 0.62;
  TLegend* lv2 = new TLegend(0.82, lv2bottom, 0.93, 0.95, "Centrality"); // w.r.t. cv
  lv2->SetFillColor(kNone);
  lv2->SetBorderSize(0);
  utils.padsetup(cv, 1, 1, "", 0.12, 0.01, 0.03, 0.2);
  
  double xmin = -0.43;
  double xmax =  12.5;
  TH1F* hf = gPad->DrawFrame(xmin, -0.09, xmax, 26);
  // int nXticks = 210;
  // int nYticks = 210;
  hf->GetXaxis()->SetTicks("+-");
  hf->GetYaxis()->SetTicks("+-");
    
// IN PROGRESS---------------
    for (int j=0; j<npta; j++) {
      for(int k=0; k<ncb; k++) {

	TGraphErrors *gc = CorrFnChi2VsTrigPt(ptabins[j], centbins[k], 1, 5, "ndf");
	if (!gc) {
	  Error("DrawChi2()", "!gc");
	  continue;
	}    
	int col = CentColor(centlow[centbins[k]], centhigh[centbins[k]]);
	utils.set_tgraph_props(gc, col, col, kFullCircle, 1.6);
	gc->Draw("epsame");
	lv2->AddEntry(gc, centLabel(centbins[k]), "l,p");
      }
    }
  
  // cv->cd();
  // TPad* overlay = new TPad("overlay", "overlay", 0,0,1,1);
  // overlay->SetFillStyle(4000); // transparent
  // overlay->Draw();
  // overlay->cd();
  
  // ltx.SetTextSize(0.07);
  // lv2->Draw();
  // overlay->cd();
  // ltx.SetTextSize(0.075);
  // ltx.DrawLatex(0.5, 0.05, "trigger p_{T} [GeV/c]");
  // ltx.SetTextAngle(90);
  // ltx.DrawLatex(0.05, 0.5, "#chi^{2} / NDF");
  // ltx.SetTextAngle(0);

  return cv;
}


TCanvas* DrawAgreement(int icent, TString opt)
{
  // Draw Chi2 on one canvas for each
  // centrality bin k.

  initialize();
  int cent1 = centlow[icent];
  int cent2 = centhigh[icent];

  if (icent==999) {
    cent1 = 0;
    cent2 = 90;
  }
  int k = icent;
  TString cname = Form("agmt_etamin%02d_cent%dto%d_%s", 
		       (int)(10*minRidgeDEta), cent1, cent2, opt.Data());
  if (!opt.IsNull()) {
    cname.Append("_");
    cname.Append(opt.Data());
  }

  TCanvas* cv = new TCanvas(cname.Data(), cname.Data(), 1400, 500);
  TLegend* lv = new TLegend(0.09, 0.72, 0.27, 0.88);
  lv->SetFillColor(kNone);
  lv->SetBorderSize(0);

  double lv2bottom = 0.62;
  TLegend* lv2 = new TLegend(0.82, lv2bottom, 0.93, 0.95, "Centrality");
  lv2->SetFillColor(kNone);
  lv2->SetBorderSize(0);
  cv->Divide(5,1);

  for (int n=1; n<=5; n++) {
    cv->cd(n);

      TH2F *gc = Agreement2DHist(k,n);
      if (!gc) {
	Error("DrawAgreement()", "Problem finding graph" );
	continue;
      }

      //      gc->GetZaxis()->SetRangeUser(1.e-1, 1000);

      gPad->SetLogx(); gPad->SetLogy(); gPad->SetLogz();

      if (opt.Contains("lego")) {
	gc->SetTitleOffset(0.8, "X");
	gc->GetXaxis()->SetLabelSize(0.07);
	gc->GetYaxis()->SetLabelSize(0.07);
	gc->GetZaxis()->SetLabelSize(0.07);
	gc->GetXaxis()->SetTitleSize(0.1);
	gc->GetYaxis()->SetTitleSize(0.1);
	gc->GetZaxis()->SetTitleSize(0.1);
	gc->GetXaxis()->SetNdivisions(8);
	gc->GetYaxis()->SetNdivisions(8);
	
	gPad->SetTopMargin(0.2);
	gPad->SetLeftMargin(0.12);
	gPad->SetRightMargin(0.01);
	gc->Draw("lego2");
      }
      else if (opt.Contains("contour"))
	gc->Draw("colz");

      if (n==1) {
	lv2->AddEntry(gc, Form("%d-%d%%",centlow[k],centhigh[k]),"l,p");
      }
      
  } // Fourier moment n
  
  cv->cd();
  TPad* overlay = new TPad("overlay", "overlay", 0,0,1,1);
  overlay->SetFillStyle(4000); // transparent
  overlay->Draw();
  overlay->cd();
  
  ltx.SetTextSize(0.07);
  const char* str = Form("#left|V_{n#Delta} / (v_{n}^{t} v_{n}^{a})#right|");
  ltx.DrawLatex(0.42, 0.92, str);
  ltx.DrawLatex(0.91, 0.92, Form("%s", centLabel(k)));
  
  for (int n=1; n<=5; n++) {
    cv->cd(n);
    ltx.SetTextSize(0.10);
    ltx.DrawLatex(0.05, 0.76, Form("n = %d", n));
  }
  return cv;
}

void AddPrelimStampToCurrentPad(float x1, float y1, float x2, float y2, TString opt)
{
  TASImage* logo = new TASImage("../commonfigs/alice-prelim.png"); // 287x323
  TPad* logopad = new TPad("logopad", "logo", x1,y1,x2,y2);
  if (opt.Contains("border")) {
    logopad->SetFillColor(kRed);
    logopad->SetBorderSize(2);
  }
  logopad->Draw();  
  logopad->cd(); 
  logo->Draw("");
  return;
}

int CentBin(double cntLo, double cntHi) 
{
  int bin = -123;
  for (int k=0; k<maxcent; k++) {
    if (centlow[k]==cntLo && centhigh[k]==cntHi)
      bin = k;
  }
  if (bin < 0 && 0)
    Warning("CentBin","No bin for %.2g - %.2g", cntLo, cntHi);

  return bin;
}

int PtBin(double ptLo, double ptHi) 
{
  int bin = -123;
  for (int i=0; i<maxpt; i++) {
    if (ptlow[i]==ptLo && pthigh[i]==ptHi)
      bin = i;
  }
  if (bin < 0)
    Warning("PtBin","No bin for %.2g - %.2g", ptLo, ptHi);

  return bin;
}

double MomCorrection(int i, int j, int k)
{
  if(0)cout<<i<<j<<k<<endl;
  // For momentum (non)conservation in v1.
  // Multiplicitave weight: w = w_t*w_a.
  // w_t = ptt - <pt^2>/<pt>
  // w_a = pta - <pt^2>/<pt>

  //  double c = 1./3 * -0.00185;

  //  (ptmean[i]-)* ptmean[j]

  double correction =  1./3 * -0.00185 * ptmean[i]* ptmean[j];
  return correction;
}

void MakeVnDVsQGraphs(int n1, int n2, int k, const char* region, const char* corrtype, TString opt)
{
  // The first two errType constants have points set at the central measured value of
  // VnDelta, with statistical and systematic error bars
  // respectively. The last two have points set at the upper and lower
  // systematic value, with statistical error bars set to be the same
  // as for meas_stat.
  const int nErrTypes = 4;
  enum errType {MEAS_STAT, MEAS_SYS, HI, LO};
  //  const char* errTypeStr[] = {"meas_stat","meas_sys","hi","lo"};
  
  // Fill a huge array with all the coefficients, then use it to set the tgraphs.
  const int nN = gNMax + 1; // 1-10, don't use 0.
  const int nI = maxpt+1;
  const int nJ = maxpt+1;
  const int nT = nErrTypes;
  double vVal[nN][nI][nJ][nT];
  double vErr[nN][nI][nJ][nT];
  double vMix[nN][nI][nJ][nT]; // Just the acceptance correction unc.

  double vSin[nN][nI][nJ][nT]; // <sin n delta phi> for sys. estimation
  double vSinErr[nN][nI][nJ][nT];

  for (int i=0; i<maxpt; i++) {
    for (int j=0; j<=i; j++) {

      double esin  = 0;  // residual sine error
      double ebw   = 0;  // bin width error \propto n / (# dphi bins)
      double ecnt  = 0;  // centrality determination
      double emix  = 0;  // acc. correction error
      double etrk  = 0;  // tracking resolution
      double ev1   = 0;  // v1
      double esys  = 0;  // quadrature sum

      TH1* hst=0, *hsys=0; // stat and sys error, but same points

      hst  = Hist(region, corrtype, i, j, k, "");
      hsys = Hist(region, corrtype, i, j, k, "sys"); 

      if (!hst) Error("MakeVnDVsQGraphs()","!hst");
      if (!hsys) Error("MakeVnDVsQGraphs()","!hsys");

      esin = SineError(i,j,k); // resid_sin / 2;

      for (int n=n1; n<=n2; n++) { // Order of coeff.
	
	// VnDelta on measured points, with statistical errors.
	double val=0, err=0, resid_sin=0, rs_err=0;
	VnDelta(n, *hst,  val, err, "");
	VnDelta(n, *hst,  resid_sin, rs_err, "sin");
	
	// VnDelta on high and low sys values. Stat err. not used.
	double val_hi=0, xxx=0;
	VnDelta(n, *hsys, val_hi, xxx, "hi");
	double val_lo=0, yyy=0;
	VnDelta(n, *hsys, val_lo, yyy, "lo");

	// Momentum conservation correction
	if (n==1 && opt.Contains("ptcons")) { // eqs. 3-7, PRL 106, 102301 (2011)
	  double factor =  ispp ? 0.212 : 0; // 7e-5; <------ trial & error: 7e-5 must be close for 0-20%
	  val    += factor*ptmean[i]*ptmean[j];
	  val_hi += factor*ptmean[i]*ptmean[j];
	  val_lo += factor*ptmean[i]*ptmean[j];
	}

	int nbins = hst->GetNbinsX();
	double mean = (val_hi + val_lo + val) / 3.;
	double ms = 
	  (val_hi-mean)*(val_hi-mean) + 
	  (val_lo-mean)*(val_lo-mean) + 
	  (val   -mean)*(val   -mean);
	double rms = TMath::Sqrt(ms/2); // sqrt(1/(N-1) sum (x_i - <x>)^2 )


	ebw = n/nbins/TMath::Sqrt(12)*val; // This is about 0.008 * n * val
	etrk = 0.01*ptmean[j]*val; // see comments in comparisons/VnDcomparison.C 
	ecnt = 0.01*val;
	emix = rms/2; // take half not to double-count stat err. // 0.25*TMath::Abs(val_hi - val_lo);
	ev1 = (n==1) ? 0.001*ptmean[i]*ptmean[j]*val : 0.;

	esys = TMath::Sqrt(ecnt*ecnt + emix*emix + ebw*ebw + etrk*etrk + ev1*ev1 + esin*esin);
	
	double pointVal = val;
	double pointErr = err;

	vVal[n][i][j][MEAS_STAT] = pointVal;
	vErr[n][i][j][MEAS_STAT] = pointErr;
	vVal[n][i][j][MEAS_SYS] = pointVal;
	vErr[n][i][j][MEAS_SYS] = esys;
	vVal[n][i][j][HI] = pointVal + esys;
	vErr[n][i][j][HI] = pointErr;
	vVal[n][i][j][LO] = pointVal - esys;
	vErr[n][i][j][LO] = pointErr;
	vSin[n][i][j][MEAS_STAT] = resid_sin;
	vSinErr[n][i][j][MEAS_STAT] = rs_err;
	vMix[n][i][j][MEAS_SYS] = emix;

      } // n
    } // j
  } // i

  // Now set the graphs from the arrays.
  for (int n=n1; n<=n2; n++) {
    for (int t=0; t<nErrTypes; t++) {
  
      TString gname = Form("V%dDeltaVsGlobalIndex_cent%dto%d",
			   n,centlow[k],centhigh[k]);
      if (t==HI)
	gname.Append("hi");
      else if (t==LO)
	gname.Append("lo");
      else if (t==MEAS_SYS)
	gname.Append("meas_sys");
      else
	gname.Append("meas_stat");

      if (opt.Contains("sine"))
	gname.Append("_sine");

      if (opt.Contains("ALL"))
	gname.Append("_ALL");
      
      TGraphErrors* g = (TGraphErrors*) gList->FindObject(gname.Data());
      if (!g) {
	g = new TGraphErrors();
	g->SetName(gname.Data());
      }

      TGraphErrors* gmix = 0;
      if (t==MEAS_SYS) {
	gmix = new TGraphErrors();
	gname.Append("_mix");
	gmix->SetName(gname.Data());
      }
      
      double xer = 0.4; // was 0.5
      for (int i=0; i<maxpt; i++) {
	for (int j=0; j<=i; j++) {
	  int q = GlobalIndex(i, j);
	  g->SetPoint(q, q+0.5, vVal[n][i][j][t]);
	  g->SetPointError(q, xer, vErr[n][i][j][t]);

	  if (opt.Contains("sine")) {
	    g->SetPoint(q, q+0.5, vSin[n][i][j][MEAS_STAT]);
	    g->SetPointError(q, xer, vSinErr[n][i][j][MEAS_STAT]);
	  }

	  if (t==MEAS_SYS) {
	    gmix->SetPoint(q, q+0.5, vVal[n][i][j][t]);
	    gmix->SetPointError(q, xer, vMix[n][i][j][t]);
	  }

	}
      }
      
      gList->Add(g);
      if (gmix) gList->Add(gmix);
    }
  }
  
  return;
}

TGraphErrors* VnDVsQ(int n, int k, const char* region, const char* corrtype, TString opt)
{

  // Graph of pair fourier coefficients VnDelta vs. the linear
  // global index q. Arguments:
  // 
  // n, k:   fourier index, centrality bin
  // region: "RIDGE", "NSJET" or "ALL"
  // corrtype:  "s", "m", "cA", "cB"
  //
  // Select opt "hi_sys" or "lo_sys" for hi and low pts
  // with central stat errors and "sys" for the central points
  // with sys. errors. Just like Hist().

  TString gname = Form("V%dDeltaVsGlobalIndex_cent%dto%d",
		       n,centlow[k],centhigh[k]);

  if (opt.Contains("hi"))
    gname.Append("hi");
  else if (opt.Contains("lo"))
    gname.Append("lo");
  else if (opt.Contains("sys"))
    gname.Append("meas_sys");
  else
    gname.Append("meas_stat");

  if (opt.Contains("sine"))
    gname.Append("_sine");

  if (opt.Contains("ALL"))
    gname.Append("_ALL");
  
  TGraphErrors* g = (TGraphErrors*) gList->FindObject(gname.Data()); 
  if (g) {
    if (0)
      Info("VnDVsQ()", "Returning preexisting graph");
    return g;
  }

  TString opt1 = (opt.Contains("sine")) ? "sine" : "";
  if (opt.Contains("ALL"))
    opt1.Append("_ALL");
  if (opt.Contains("ptcons"))
    opt1.Append("_ptcons");

  MakeVnDVsQGraphs(1, gNMax, k, region, corrtype, opt1);
  g = (TGraphErrors*) gList->FindObject(gname.Data()); 
  if (!g) {
    Error("VnDVsQ()", "No graph %s", gname.Data());
    gSystem->Exit(-1);
  }
  return g;
}

double MeanPt(int i, int j, int k, TString t_or_a, TString opt)
{
  // For combined centralities, map k's to the closest available bin
  // double cb1[] =  {0, 0, 0,2,2, 1,3,0,    0, 1, 2, 3, 4,  5, 10, 20, 30, 40, 60 };
  // double cb2[] =  {10,20,2,5,10,3,5,5,    1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 90 };
// OBJ: TProfile2D	hMeanPtTrgCen1	cen=1 (0.5) : 0 at: 0x102af06e0
//  OBJ: TProfile2D	hMeanPtAssCen1	cen=1 (0.5) : 0 at: 0x108308e60
//  OBJ: TProfile2D	hMean2PtTrgCen1	cen=1 (0.5) : 0 at: 0x108309480
//  OBJ: TProfile2D	hMean2PtAssCen1	cen=1 (0.5) : 0 at: 0x108309ad0
//  OBJ: TProfile2D	hMeanPtTrgCen2	cen=2 (1.5) : 0 at: 0x10830a120
//  OBJ: TProfile2D	hMeanPtAssCen2	cen=2 (1.5) : 0 at: 0x10830a770
//  OBJ: TProfile2D	hMean2PtTrgCen2	cen=2 (1.5) : 0 at: 0x10830adc0
//  OBJ: TProfile2D	hMean2PtAssCen2	cen=2 (1.5) : 0 at: 0x10830b410
//  OBJ: TProfile2D	hMeanPtTrgCen3	cen=3 (2.5) : 0 at: 0x10830ba60
//  OBJ: TProfile2D	hMeanPtAssCen3	cen=3 (2.5) : 0 at: 0x10830c0b0
//  OBJ: TProfile2D	hMean2PtTrgCen3	cen=3 (2.5) : 0 at: 0x10830c700
//  OBJ: TProfile2D	hMean2PtAssCen3	cen=3 (2.5) : 0 at: 0x10830cd50
//  OBJ: TProfile2D	hMeanPtTrgCen4	cen=4 (3.5) : 0 at: 0x10830d3a0
//  OBJ: TProfile2D	hMeanPtAssCen4	cen=4 (3.5) : 0 at: 0x10830d9f0
//  OBJ: TProfile2D	hMean2PtTrgCen4	cen=4 (3.5) : 0 at: 0x10830e040
//  OBJ: TProfile2D	hMean2PtAssCen4	cen=4 (3.5) : 0 at: 0x10830e690
//  OBJ: TProfile2D	hMeanPtTrgCen5	cen=5 (4.5) : 0 at: 0x10830ece0
//  OBJ: TProfile2D	hMeanPtAssCen5	cen=5 (4.5) : 0 at: 0x10830f330
//  OBJ: TProfile2D	hMean2PtTrgCen5	cen=5 (4.5) : 0 at: 0x10830f980
//  OBJ: TProfile2D	hMean2PtAssCen5	cen=5 (4.5) : 0 at: 0x10830ffd0
//  OBJ: TProfile2D	hMeanPtTrgCen6	cen=6 (7.5) : 0 at: 0x108310620
//  OBJ: TProfile2D	hMeanPtAssCen6	cen=6 (7.5) : 0 at: 0x108310c70
//  OBJ: TProfile2D	hMean2PtTrgCen6	cen=6 (7.5) : 0 at: 0x1083112c0
//  OBJ: TProfile2D	hMean2PtAssCen6	cen=6 (7.5) : 0 at: 0x108311910
//  OBJ: TProfile2D	hMeanPtTrgCen7	cen=7 (15.0) : 0 at: 0x108311f60
//  OBJ: TProfile2D	hMeanPtAssCen7	cen=7 (15.0) : 0 at: 0x1083125b0
//  OBJ: TProfile2D	hMean2PtTrgCen7	cen=7 (15.0) : 0 at: 0x108312c00
//  OBJ: TProfile2D	hMean2PtAssCen7	cen=7 (15.0) : 0 at: 0x108313250
//  OBJ: TProfile2D	hMeanPtTrgCen8	cen=8 (25.0) : 0 at: 0x1083138a0
//  OBJ: TProfile2D	hMeanPtAssCen8	cen=8 (25.0) : 0 at: 0x108313ef0
//  OBJ: TProfile2D	hMean2PtTrgCen8	cen=8 (25.0) : 0 at: 0x108314540
//  OBJ: TProfile2D	hMean2PtAssCen8	cen=8 (25.0) : 0 at: 0x108314b90
//  OBJ: TProfile2D	hMeanPtTrgCen9	cen=9 (35.0) : 0 at: 0x1083151e0
//  OBJ: TProfile2D	hMeanPtAssCen9	cen=9 (35.0) : 0 at: 0x108315830
//  OBJ: TProfile2D	hMean2PtTrgCen9	cen=9 (35.0) : 0 at: 0x108315e80
//  OBJ: TProfile2D	hMean2PtAssCen9	cen=9 (35.0) : 0 at: 0x1083164d0
//  OBJ: TProfile2D	hMeanPtTrgCen10	cen=10 (45.0) : 0 at: 0x108316b20
//  OBJ: TProfile2D	hMeanPtAssCen10	cen=10 (45.0) : 0 at: 0x108317170
//  OBJ: TProfile2D	hMean2PtTrgCen10	cen=10 (45.0) : 0 at: 0x1083177c0
//  OBJ: TProfile2D	hMean2PtAssCen10	cen=10 (45.0) : 0 at: 0x108317e10
//  OBJ: TProfile2D	hMeanPtTrgCen11	cen=11 (55.0) : 0 at: 0x108318460
//  OBJ: TProfile2D	hMeanPtAssCen11	cen=11 (55.0) : 0 at: 0x108318ab0
//  OBJ: TProfile2D	hMean2PtTrgCen11	cen=11 (55.0) : 0 at: 0x108319100
//  OBJ: TProfile2D	hMean2PtAssCen11	cen=11 (55.0) : 0 at: 0x108319750
//  OBJ: TProfile2D	hMeanPtTrgCen12	cen=12 (75.0) : 0 at: 0x108319da0
//  OBJ: TProfile2D	hMeanPtAssCen12	cen=12 (75.0) : 0 at: 0x10831a3f0
//  OBJ: TProfile2D	hMean2PtTrgCen12	cen=12 (75.0) : 0 at: 0x10831aa40
//  OBJ: TProfile2D	hMean2PtAssCen12	cen=12 (75.0) : 0 at: 0x10831b090

  int k2 = -1;
  if (k == CentBin(0,  1))   k2 = 1;  
  else if (k == CentBin(1,  2))   k2 = 2;  
  else if (k == CentBin(2,  3))   k2 = 3;  
  else if (k == CentBin(3,  4))   k2 = 4;  
  else if (k == CentBin(4,  5))   k2 = 5;  
  else if (k == CentBin(5, 10))   k2 = 6;  
  else if (k == CentBin(10,20))   k2 = 7;  
  else if (k == CentBin(20,30))   k2 = 8;  
  else if (k == CentBin(30,40))   k2 = 9;  
  else if (k == CentBin(40,50))   k2 = 10; 
  else if (k == CentBin(50,60))   k2 = 11; 
  else if (k == CentBin(60,90))   k2 = 12; 
  else if (k == CentBin(0, 10))   k2 = 5;    //CentBin(4, 5);   
  else if (k == CentBin(0, 20))   k2 = 6;    //CentBin(5, 10);  
  else if (k == CentBin(0,  2))   k2 = 1;    //CentBin(0, 1);   
  else if (k == CentBin(2,  5))   k2 = 4;    //CentBin(3, 4);   
  else if (k == CentBin(2, 10))   k2 = 6;    //CentBin(5, 10);  
  else if (k == CentBin(1,  3))   k2 = 2;    //CentBin(1, 2);   
  else if (k == CentBin(3,  5))   k2 = 4;    //CentBin(3, 4);   
  else if (k == CentBin(0,  5))   k2 = 3;    //CentBin(2, 3);   

  else if (k == CentBin(0, 0))   k2 = 1;    //pp
  else {
    Error("MeanPt()", "cent bin %d (%d-%d%%) not assigned", k, centlow[k], centhigh[k]);
  }

  const char* m = opt.Contains("squared") ? "Mean2" : "Mean";

  TString prName_t = Form("h%sPtTrgCen%d", m, k2);
  TString prName_a = Form("h%sPtAssCen%d", m, k2);
  TProfile2D* hmpt_t = (TProfile2D*)fin->Get(prName_t.Data());
  TProfile2D* hmpt_a = (TProfile2D*)fin->Get(prName_a.Data());

  if (!hmpt_t) {
    Error("MeanPt", "Problem finding TProfile %s. Input %d", prName_t.Data(), k);
    gSystem->Exit(-1);
  }
  //  cout << hmpt_t->GetTitle() << " " << centlow[k] << " " << centhigh[k] << endl;

  double mptt = hmpt_t->GetBinContent(i+1, j+1);
  double mpta = hmpt_a->GetBinContent(i+1, j+1);

  if (0) {
    // x axis is trig pt,  y axis assc pt
    TAxis* ax = hmpt_t->GetXaxis();
    TAxis* ay = hmpt_t->GetYaxis();
    Printf("trig %.3g - %.3g, assc: %.3g - %.3g, <pt_t,a> %.3g, %.3g",
	   ax->GetBinLowEdge(i+1), 
	   ax->GetBinUpEdge(i+1), 
	   ay->GetBinLowEdge(j+1), 
	   ay->GetBinUpEdge(j+1),
	   mptt, mpta );
  }

  if (t_or_a.Contains("t"))
    return mptt;
  if (t_or_a.Contains("a"))
    return mpta;
  else
    return -1;
}


int CentColor(double cen1, double cen2) 
{
  if (cen1==0 && cen2==10)
    return kGray+3;
  if (cen1==0 && cen2==20)
    return kMagenta+3;
  if (cen1==0 && cen2==2)
    return kBlack;
  if (cen1==2 && cen2==10)
    return kRed;
  if (cen1==10 && cen2==20)
    return kOrange-3;
  if (cen1==20 && cen2==30)
    return kGreen+1;
  if (cen1==30 && cen2==40)
    return kYellow+2;
  if (cen1==40 && cen2==50)
    return kAzure;
  if (cen1==50 && cen2==60)
    return kBlue-1;
  if (cen1==60 && cen2==90)
    return kViolet;

  return kBlack;
}

int CentColorPale(double cen1, double cen2) 
{
  if (cen1==0 && cen2==10)
    return kGray+1;
  if (cen1==0 && cen2==20)
    return kPink-1; // kPink-6;
  if (cen1==0 && cen2==2)
    return kGray+1;
  if (cen1==2 && cen2==10)
    return kRed-9;
  if (cen1==10 && cen2==20)
    return kOrange-4;
  if (cen1==20 && cen2==30)
    return kGreen-7;
  if (cen1==30 && cen2==40)
    return kYellow-8;
  if (cen1==40 && cen2==50)
    return kAzure-4;
  if (cen1==50 && cen2==60)
    return kBlue-9;
  if (cen1==60 && cen2==90)
    return kViolet-9;

  return kBlack;
}

int PtColor(double p1, double p2)
{
 if (p1==0.25 && p2==0.5)
    return kGray+1;
 if (p1==0.5 && p2==0.75)
    return kYellow+3;
 if (p1==0.75 && p2==1.0)
   return kPink-6;
 if (p1==1.0 && p2==1.5)
   return kBlack;
 if (p1==1.5 && p2==2.0)
   return kRed;
 if (p1==2.0 && p2==2.5)
   return kOrange-3;
 if (p1==2.5 && p2==3)
   return kGreen+2;
 if (p1==3 && p2==4)
   return kBlue;
 if (p1==4 && p2==5)
   return kViolet + 2;
 if (p1==5 && p2==6)
   return kViolet+10;
 if (p1==4 && p2==6)
   return kViolet + 2;
 if (p1==6 && p2==8)
   return kOrange-7;
 if (p1==8 && p2==15)
   return kOrange+9;

 return kBlack; 
}

int PtColorPale(double p1, double p2)
{
 if (p1==0.25 && p2==0.5)
    return kGray;
 if (p1==0.5 && p2==0.75)
    return kYellow+1;
 if (p1==0.75 && p2==1.0)
   return kPink-4;
 if (p1==1.0 && p2==1.5)
   return kGray+1;
 if (p1==1.5 && p2==2.0)
   return kRed-7;
 if (p1==2.0 && p2==2.5)
   return kOrange-4;
 if (p1==2.5 && p2==3)
   return kSpring-4;
 if (p1==3 && p2==4)
   return kAzure+6;
 if (p1==4 && p2==5)
   return kViolet+1;
 if (p1==5 && p2==6)
   return kViolet+8;
 if (p1==4 && p2==6)
   return kViolet+1;
 if (p1==6 && p2==8)
   return kOrange-9;
 if (p1==8 && p2==15)
   return kOrange+6;

 return kGray; 
}

int PtaColor(double p1, double p2)
{
 if (p1==0.25 && p2==0.5)
    return kGray+1;
 if (p1==0.5 && p2==0.75)
   return kBlack;
 if (p1==0.75 && p2==1.0)
   return kRed+1;
 if (p1==1.0 && p2==1.5)
   return kOrange-3;//kBlack;
 if (p1==1.5 && p2==2.0)
   return kGreen+2;//kRed;
 if (p1==2.0 && p2==2.5)
   return kAzure-3;//kOrange-3;
 if (p1==2.5 && p2==3)
   return kViolet + 2;//kBlue;
 if (p1==3 && p2==4)
   return kViolet + 10;//kBlue;
 if (p1==4 && p2==5)
   return kMagenta+3; //kViolet + 2;
 if (p1==5 && p2==6)
   return kViolet;
 if (p1==4 && p2==6)
   return kMagenta-5;
 if (p1==6 && p2==8)
   return kGray+2;
 if (p1==8 && p2==15)
   return kGray+3;

 return kBlack; 
}

int PtaColorPale(double p1, double p2)
{
 if (p1==0.25 && p2==0.5)
    return kGray;
 if (p1==0.5 && p2==0.75)
   return kGray+2;
 if (p1==0.75 && p2==1.0)
   return kPink-2;
 if (p1==1.0 && p2==1.5)
   return kOrange-4;
 if (p1==1.5 && p2==2.0)
   return kGreen-3;
 if (p1==2.0 && p2==2.5)
   return kAzure-4;
 if (p1==2.5 && p2==3)
   return kViolet;
 if (p1==3 && p2==4)
   return kViolet + 2;
 if (p1==4 && p2==5)
   return kMagenta+3;
 if (p1==5 && p2==6)
   return kViolet+10;
 if (p1==4 && p2==6)
   return kMagenta-5;
 if (p1==6 && p2==8)
   return kGray+2;
 if (p1==8 && p2==15)
   return kGray+3;

 return kGray; 
}

void SaveCanvases(TObjArray* canvases, const char* fileName)
{
  TFile* f = new TFile(fileName, "recreate");

  for (int n=0; n<canvases->GetEntries(); n++) {
    TCanvas* c = (TCanvas*)canvases->At(n);
    c->Write(c->GetTitle());
  }
  f->Close();
  return;
}

void SaveCanvasesFromFile(const char* rootFile, const char* targetDir, const char* tag, const char* fileType)
{
  // Get a list of canvases from rootFile into array, then save each
  // to its own file in targetDir/. fileType = "eps", "pdf", "C",
  // "png", etc. Not all formats have been tested.
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  TString name = "";
  TString base(targetDir);
  TFile *cFile = new TFile(rootFile, "read");
  cout << cFile->GetName() << endl;
  TObjArray* cList = GetObjectsFromFile(*cFile, "TCanvas");

  for (int n=0; n<cList->GetEntries(); n++) {
    TCanvas* c = (TCanvas*)cList->At(n);
    if (c) {
      name = "";
      name = base;
      name += TString("/");
      name += TString(fileType);
      name += TString("/");
      name += TString(c->GetTitle());
      name += TString(".");
      name += TString(fileType);
      cout<<name.Data()<<endl;

      c->Draw();
      c->Modified();
      c->Update();
      c->SaveAs(name.Data());
    }
    else
      Error("SaveCanvasesFromFile()", "!c");
  }

  if (1) {
    utils.print_pdf(cList, Form("%s/all-figs%s", targetDir, tag), "pdf");
  }
  
  return;
}

//void SaveCanvases(TObjArray* canvases, const char* fileName)
//void SaveCanvasesFromFile(const char* rootFile, const char* targetDir, const char* fileType)
TObjArray* GetObjectsFromFile(TFile& file, TString clname, TString dir)
{
  file.cd(dir.Data());

  TObjArray* objList = new TObjArray();
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  
  while ((key=(TKey*)next())) {
    TString className(key->GetClassName());
    TString keyName(key->GetName());
    if (1) 
      printf("%10s %20s\n", className.Data(), keyName.Data());
    
    if (className.Contains(clname)) {
      objList->Add(gDirectory->Get(keyName.Data()));
    }
  }

  cout << objList->GetEntries() << " objects retrieved from "
       << file.GetName()  << "/" << gDirectory->GetName() 
       << endl;

  return objList;
}
