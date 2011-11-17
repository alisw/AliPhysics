/**
 * @file   Draw123.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Nov 17 11:07:15 2011
 * 
 * @brief  This scripts draws the energy loss distribution for single,
 * double, and triple hits in the FMD as resolved by the sharing
 * filter 
 * 
 * @ingroup pwg2_forward_scripts_qa
 */
#ifndef __CINT__
# include <TH1.h>
# include <TH2.h>
# include <TList.h>
# include <TFile.h>
# include <TString.h>
# include <TError.h>
# include <TPad.h>
# include <TCanvas.h>
# include <TLine.h>
# include <TLatex.h>
# include <TStyle.h>
# include <TLegend.h>
#else
class TList;
#endif

/** 
 * Draw the energy loss spectra of single, double, and triple hits in
 * a particular ring 
 * 
 * @param p Parent list 
 * @param d Detector 
 * @param r Ring 
 *
 * @deprecated Use the QATrender instead
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawRing123(TList* p, UShort_t d, Char_t r)
{
  if (!p) return;

  TList* ring = static_cast<TList*>(p->FindObject(Form("FMD%d%c",d,r)));
  if (!ring) { 
    Error("Draw123", "List FMD%d%c not found in %s",d,r,p->GetName());
    return;
  }
  
  TH1* one   = static_cast<TH1*>(ring->FindObject("singleEloss"));
  TH1* two   = static_cast<TH1*>(ring->FindObject("doubleEloss"));
  TH1* three = static_cast<TH1*>(ring->FindObject("tripleEloss"));
  if (!one || !two || !three) { 
    Error("DrawRing123", "Histograms of Eloss not found in FMD%d%c", d, r);
    return;
  }
  one->SetStats(0);
  one->SetTitle(Form("FMD%d%c", d, r));
  one->GetXaxis()->SetRangeUser(0, 8);

  gPad->SetLogy();
  gPad->SetFillColor(0);

  one->Draw();
  if (two)   two->Draw("same");
  if (three) three->Draw("same");

  TLegend* l = new TLegend(.6, .6, .95, 1);
  l->SetFillColor(0);
  l->SetBorderSize(0);
  l->AddEntry(one);
  if (two)   l->AddEntry(two);
  if (three) l->AddEntry(three);
  l->Draw();

  gPad->cd();
}


/** 
 * Draw the energy loss distribution of singles, doubles, and triples
 * as given by the sharing filter 
 * 
 * @param filename Input file name
 * @param folder   Input folder (TList) in input file
 *
 * @deprecated Use the QATrender instead
 * @ingroup pwg2_forward_scripts_qa
 */
void
Draw123(const char* filename="forward.root", 
	const char* folder="ForwardResults")
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(.6);
  
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("Draw123", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get(folder));
  if (!forward) { 
    Error("Draw123", "List %s not found in %s", folder, filename);
    return;
  }

  TList* sf = static_cast<TList*>(forward->FindObject("fmdSharingFilter"));
  if (!sf) { 
    Error("Draw123", "List fmdSharingFilter not found in Forward");
    return;
  }
  
  TCanvas* c = new TCanvas("123", 
			   "singles, doubles, and tripples", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  c->Divide(3, 2, 0, 0);
  
  c->cd(1); DrawRing123(sf, 1, 'I');
  c->cd(2); DrawRing123(sf, 2, 'I');
  c->cd(5); DrawRing123(sf, 2, 'O');
  c->cd(3); DrawRing123(sf, 3, 'I');
  c->cd(6); DrawRing123(sf, 3, 'O');
  TVirtualPad* p = c->cd(4);
  // p->SetTopMargin(0.05);
  p->SetRightMargin(0.15);
  p->SetFillColor(0);
  TH2D* highCuts = static_cast<TH2D*>(sf->FindObject("highCuts"));
  if (highCuts) highCuts->Draw("colz");
  c->cd();
  c->SaveAs("123.png");
}
 
//
// EOF
//
