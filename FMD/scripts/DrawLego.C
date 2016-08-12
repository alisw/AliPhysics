//____________________________________________________________________
//
// Script to plot the lego information 
// 
#ifndef __CINT__
# include <TFile.h>
# include <TCanvas.h>
# include <TPad.h>
# include <TH1D.h>
# include <TString.h>
# include <TLatex.h>
# include <TLegend.h>
# include <THStack.h>
# include <TStyle.h>
# include <iostream>
using namespace std;
#endif
/** @defgroup FMD_lego_script LEGO script
    @ingroup FMD_script 
*/
//____________________________________________________________________
/** @ingroup FMD_lego_script
    @param ang 
    @return  */
Float_t 
deg2eta(Float_t ang) 
{
  if (ang == 180) ang -= .001;
  if (ang == 0)   ang += .001;
  Float_t theta = ang * TMath::Pi() / 180;
  Float_t eta   = - TMath::Log(TMath::Tan(theta / 2));
  return eta;
}


//____________________________________________________________________
/** @ingroup FMD_lego_script
    @param which 
    @param what 
    @param back 
    @return  */
TH1* 
getHisto(const Char_t* which, const Char_t* what, TH1* back=0) 
{
  TFile* file = TFile::Open(Form("Lego_%s.root", which), "READ");
  if (!file) {
    cerr << "Couldn't open the file 'Lego_" << which << ".root" 
	 << endl;
    return 0;
  }
  TH1D*  h1d  = static_cast<TH1D*>(file->Get(Form("h%s_py", what)));
  if (!h1d) {
    cerr << "Couldn't find h" << what << "_py in " 
	 <<  "Lego_" << which << ".root" << endl;
    return 0;
  }

  TAxis* xaxis = h1d->GetXaxis();
  Int_t n = xaxis->GetNbins();
  TArrayF bins(n-1);
  for (Int_t i = n-1; i > 1; i--) {
    Float_t ang = xaxis->GetBinUpEdge(i);
    Float_t eta = deg2eta(ang);
    bins[n-i-1] = deg2eta(xaxis->GetBinUpEdge(i));
  }
  bins[n-2] = deg2eta(xaxis->GetBinLowEdge(2));
  
  TH1F* heta = new TH1F(Form("%s_eta", what), h1d->GetTitle(),
			n-2, bins.fArray);
  heta->SetXTitle("#eta");
  heta->GetXaxis()->SetTitleSize(0);
  heta->GetXaxis()->SetTitleOffset(1.5);
  heta->GetXaxis()->SetTitleFont(132);
  heta->GetXaxis()->SetLabelFont(132);
  heta->SetYTitle(Form("%s per degree", h1d->GetTitle()));
  heta->GetYaxis()->SetTitleOffset(1.5);
  heta->GetYaxis()->SetTitleFont(132);
  heta->GetYaxis()->SetLabelFont(132);
  heta->GetYaxis()->SetTitleSize(0);
  heta->GetZaxis()->SetTitle(heta->GetTitle());
  heta->GetZaxis()->SetTitleOffset(1.5);
  heta->GetZaxis()->SetTitleFont(132);
  heta->GetZaxis()->SetLabelFont(132);
  heta->GetYaxis()->SetTitleOffset(1.5);
  for (Int_t i = 2; i < n; i++) {
    Float_t ang = xaxis->GetBinUpEdge(i);
    Float_t eta = deg2eta(ang);
    Float_t y   = h1d->GetBinContent(i);
    Float_t j   = heta->FindBin(eta);
    if (back) {
      Float_t b = back->GetBinContent(j);
      if (y - b <= 0) y = .000001;
      else            y = y - b;
    }
    // cout << i << ": " << ang << " -> " << eta << " = " << y << endl;    
    heta->SetBinContent(j, y);
  }
  return heta;
}

//____________________________________________________________________
/** @ingroup FMD_lego_script
    @param what 
*/
void
drawLego(const char* what="abso") 
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetLabelFont(132, "xyz");
  gStyle->SetTitleFont(132, "xyz");
  gStyle->SetTitleOffset(1.5, "y");

  TH1* nothing = getHisto("Nothing", what);
  TH1* its     = getHisto("ITS",     what, nothing);
  TH1* fmd     = getHisto("FMD",     what, nothing);
  TH1* pipe    = getHisto("PIPE",    what, nothing);
  TH1* inner   = getHisto("Inner",   what);

  if (!inner || !pipe || !fmd || !its || !nothing) {
    cerr << "Failed to get a histogram!" << endl;
    return;
  }
  TCanvas* c = new TCanvas(Form("single_%s", what), 
			   Form("Single %s", what), 800, 800);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(0);
  TPad* p1 = new TPad("p1", "p1", 0.0, 0.5, 0.5, 1.0, 0, 0, 0);
  TPad* p2 = new TPad("p2", "p2", 0.5, 0.5, 1.0, 1.0, 0, 0, 0);
  TPad* p3 = new TPad("p3", "p3", 0.0, 0.0, 0.5, 0.5, 0, 0, 0);
  TPad* p4 = new TPad("p4", "p4", 0.5, 0.0, 1.0, 0.5, 0, 0, 0);
  
  TLatex* latex = new TLatex(0,0,"");
  latex->SetTextFont(132);

  float logmax = inner->GetMaximum();
  float logmin = .001 * logmax;
  float latexy = 6 * logmin;
  
  c->cd();
  p1->SetLogy();
  p1->SetGridy();
  p1->SetTopMargin(0.15);
  p1->SetBottomMargin(0);
  p1->SetRightMargin(0);
  p1->Draw();
  p1->cd();
  fmd->SetFillColor(4);
  fmd->GetYaxis()->SetRangeUser(logmin, logmax);
  fmd->GetYaxis()->SetTitleSize(.04);
  fmd->Draw();
  latex->DrawLatex(-1, latexy, "FMD only");
  
  c->cd();
  p2->SetLogy();
  p2->SetGridy();
  p2->SetTopMargin(0.15);
  p2->SetBottomMargin(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(0.15);
  p2->Draw();
  p2->cd();
  its->SetFillColor(2);
  its->GetYaxis()->SetRangeUser(logmin, logmax);
  its->Draw();
  latex->DrawLatex(-1, latexy, "ITS only");

  c->cd();
  p3->SetLogy();
  p3->SetGridy();
  p3->SetTopMargin(0);
  p3->SetLeftMargin(0.15);
  p3->SetRightMargin(0);
  p3->Draw();
  p3->cd();
  pipe->SetFillColor(3);
  pipe->GetYaxis()->SetRangeUser(logmin, logmax);
  pipe->Draw();
  latex->DrawLatex(-1, latexy, "PIPE only");
  
  c->cd();
  p4->SetLogy();
  p4->SetGridy();
  p4->SetTopMargin(0);
  p4->SetLeftMargin(0.);
  p4->SetRightMargin(0.15);
  p4->Draw();
  p4->cd();
  inner->GetYaxis()->SetRangeUser(logmin, logmax);
  inner->GetXaxis()->SetTitleSize(.04);
  inner->SetFillColor(5);
  inner->Draw();
  latex->DrawLatex(-1, latexy, "PIPE, ITS, FMD and Air");
  
  c->Modified();
  c->cd();
  c->Print(Form("%s_single.png", what));

  TCanvas* accum = new TCanvas(Form("accum_%s", what), 
			       Form("Accumalted %s",what), 
			       800, 500);
  accum->SetLogy();
  accum->SetFillColor(0);
  accum->SetBorderMode(0);
  accum->SetBorderSize(0);

  THStack* stack = new THStack("stack", "Stack");
  nothing->SetFillColor(6);
  stack->Add(nothing);
  stack->Add(pipe);
  stack->Add(fmd);
  stack->Add(its);

  TLegend* legend = new TLegend(.15, .65, .27, .95);
  legend->SetFillColor(0);
  legend->SetBorderSize(1);
  legend->AddEntry(its, "ITS", "f");
  legend->AddEntry(fmd, "FMD", "f");
  legend->AddEntry(pipe, "PIPE", "f");
  legend->AddEntry(nothing, "Air", "f");
  
  stack->SetMinimum(nothing->GetMinimum());
  stack->SetMaximum(logmax);
  stack->Draw();
  stack->GetXaxis()->SetTitle("#eta");
  // stack->GetYaxis()->SetRangeUser(, );
  stack->GetYaxis()->SetTitle(fmd->GetTitle());
  legend->Draw();
  accum->Modified();
  
  accum->Modified(); 
  accum->cd();
  // accum->Print(Form("%s_accum.eps", what));
  accum->Print(Form("%s_accum.png", what));
  
}

  
  
//____________________________________________________________________
/** @ingroup FMD_lego_script
 */
void 
DrawLego() 
{
  drawLego("abso");
  drawLego("radl");
  drawLego("gcm2");
}

//____________________________________________________________________
//
// EOF
//

  
