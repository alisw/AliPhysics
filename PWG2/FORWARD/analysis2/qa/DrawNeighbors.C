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
#else
class TList;
#endif
/** 
 * Draw the correlation of neighboring strips before/after merging 
 * 
 * @param p  List
 * @param d  Detector
 * @param r  Ring
 *
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawRingNeighbors(TList* p, UShort_t d, Char_t r)
{
  if (!p) return;

  TList* ring = static_cast<TList*>(p->FindObject(Form("FMD%d%c",d,r)));
  if (!ring) { 
    Error("DrawNeighbors", "List FMD%d%c not found in %s",d,r,p->GetName());
    return;
  }
  
  TH2* before = static_cast<TH2D*>(ring->FindObject("neighborsBefore"));
  if (!before) { 
    Error("DrawRingNeighbors", "Histogram neighborsBefore not found in FMD%d%c",
	  d, r);
    return;
  }
  TH2* after = static_cast<TH2D*>(ring->FindObject("neighborsAfter"));
  if (!after) { 
    Error("DrawRingNeighbors", "Histogram neighborsAfter not found in FMD%d%c",
	  d, r);
    return;
  }
  gPad->SetLogz();
  gPad->SetFillColor(0);
  TPad* pad = (TPad*)gPad;
  if (d == 3) { 
    pad->SetPad(pad->GetXlowNDC(), pad->GetYlowNDC(), .99, 
		 pad->GetYlowNDC()+pad->GetHNDC());
    pad->SetRightMargin(0.15);
  }
  // gStyle->SetTitleY(gPad->GetBottomMargin());

  before->SetTitle(Form("FMD%d%c",d,r));
  before->Draw("colz");
  after->Draw("same box");

  before->GetXaxis()->SetRangeUser(-.5, 2);
  before->GetYaxis()->SetRangeUser(-.5, 2);

  TLatex* ltx = new TLatex(gPad->GetLeftMargin()+.01, 
			   gPad->GetBottomMargin()+.01, 
			   before->GetTitle());
  ltx->SetNDC();
  ltx->SetTextSize(.07);
  ltx->Draw();

  gPad->cd();
}

/** 
 * Draw the correlation of neighboring strips before/after merging 
 * 
 * @param filename 
 *
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawNeighbors(const char* filename="forward.root", 
	      const char* folder="ForwardResults")
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleX(.1);
  gStyle->SetTitleY(.1);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptTitle(0);

  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawNeighbors", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get(folder));
  if (!forward) { 
    Error("DrawNeighbors", "List %s not found in %s", folder, filename);
    return;
  }

  TList* sf = static_cast<TList*>(forward->FindObject("fmdSharingFilter"));
  if (!sf) { 
    Error("DrawNeighbors", "List fmdSharingFilter not found in Forward");
    return;
  }
  
  TCanvas* c = new TCanvas("neighbors", "Correlation of Neighbor strips", 
			   900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->Divide(3, 2, 0, 0);
  
  c->cd(1); DrawRingNeighbors(sf, 1, 'I');
  c->cd(2); DrawRingNeighbors(sf, 2, 'I');
  c->cd(5); DrawRingNeighbors(sf, 2, 'O');
  c->cd(3); DrawRingNeighbors(sf, 3, 'I');
  c->cd(6); DrawRingNeighbors(sf, 3, 'O');
  c->cd(4)->SetFillColor(0);
  c->cd();
  c->SaveAs("neighbors.png");
}

  
  
 
//
// EOF
// 
