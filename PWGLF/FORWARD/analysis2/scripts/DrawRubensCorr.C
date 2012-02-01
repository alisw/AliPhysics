/** 
 * Draw Rubens corrections 
 * 
 * @param fname 
 * @param hname 
 *
 * @ingroup pwg2_forward_scripts
 */
void
DrawRubensCorr(const char* fname="rubensRatio.root",
	       const char* hname = "dNdEtaCor1D_cls")
{
  TFile* f = TFile::Open(fname, "READ");
  if (!f) { 
    Error("DrawRubensCorr", "%s not found", fname);
    return;
  }
  TH1D* h = static_cast<TH1D*>(f->Get(hname));
  if (!h) {
    Error("DrawRubensCorr", "%s not found in %s", hname, fname);
    return;
  }

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  TCanvas* c = new TCanvas("c", "c", 800, 800);
  c->SetFillColor(0);
  c->SetFillStyle(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetRightMargin(0.03);
  c->SetLeftMargin(0.15);
  c->SetTopMargin(0.03);
  
  h->GetXaxis()->SetTitleFont(132);
  h->GetXaxis()->SetLabelFont(132);
  h->GetXaxis()->SetDecimals();
  h->GetYaxis()->SetTitleFont(132);
  h->GetYaxis()->SetLabelFont(132);
  h->GetYaxis()->SetDecimals();
  h->GetYaxis()->SetTitleOffset(1.8);
  h->SetYTitle("Clusters/Tracklets");
  h->Draw();
  
  TF1* g = new TF1("g", "pol2", -2.2, 2.2);
  g->SetLineWidth(2);
  h->Fit(g, "NQ");
  g->Draw("same");
  Info("DrawRubensCorr", "Result of Fit:\n"
       "  chi^2/NDF = %7.4f / %2d = %7.4f\n"
       "  A =   %+4.2f +/-   %5.3f\n"
       "  B =  %+5.3f +/-  %6.4f\n"
       "  C = %+6.4f +/- %7.5f\n"
       "  f(x)=%4.2f%+5.3f*x%+6.4f*x*x",
       g->GetChisquare(), g->GetNDF(), g->GetChisquare()/g->GetNDF(), 
       g->GetParameter(0), g->GetParError(0), 
       g->GetParameter(1), g->GetParError(1), 
       g->GetParameter(2), g->GetParError(2), 
       g->GetParameter(0), g->GetParameter(1), g->GetParameter(2));

  TLatex* ltx = new TLatex(.2, .9, Form("f(x)=%4.2f%+5.3fx%+6.4fx^{2}",
					g->GetParameter(0), 
					g->GetParameter(1), 
					g->GetParameter(2)));
  ltx->SetNDC();
  ltx->SetTextFont(132);
  ltx->Draw();

  c->SaveAs("rubens_corr.png");
}

//
// EOF
//
