void
ClearCanvas(TCanvas* c)
{
  c->SetLeftMargin(.1);
  c->SetRightMargin(.05);
  c->SetBottomMargin(.1);
  c->SetTopMargin(.05);
  c->Clear();
}

void
DrawELossFits(const char* fname, const char* option="err")
{
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  TFile* file = TFile::Open(fname, "READ");
  if (!file) { 
    Error("DrawELossFits", "Failed to open %s", fname);
    return;
  }
  TString pname(fname);
  pname.ReplaceAll(".root", ".pdf");

  AliFMDCorrELossFit* fits = 
    static_cast<AliFMDCorrELossFit*>(file->Get("elossfits"));
  if (!fits) { 
    Error("DrawELossFits", "Object 'elossfits' not found in %s", fname);
    return;
  }

  TCanvas* c = new TCanvas("c", "c", 800 / TMath::Sqrt(2), 800);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->Print(Form("%s[", pname.Data()));
  
  gStyle->SetOptStat(0);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(.7);
  gStyle->SetTitleY(1);
  gStyle->SetTitleW(.3);
  gStyle->SetTitleH(.1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameBorderMode(1);

  ClearCanvas(c);

  TLatex* ll = new TLatex(.5,.8, fname);
  ll->SetTextAlign(22);
  ll->SetTextSize(0.05);
  ll->SetNDC();
  ll->Draw();

  TLatex* l = new TLatex(.5,.8, fname);
  l->SetNDC();
  l->SetTextSize(0.03);
  l->SetTextFont(132);
  l->SetTextAlign(12);
  l->DrawLatex(0.2, 0.70, "1^{st} page is a summary of fit parameters");
  l->DrawLatex(0.2, 0.67, "2^{nd} page is a summary of relative errors");
  l->DrawLatex(0.2, 0.64, "Subsequent pages shows the fitted functions");
  l->DrawLatex(0.3, 0.60, "Black line is the full fitted function");
  l->DrawLatex(0.3, 0.57, "Coloured lines are the individual N-mip components");
  l->DrawLatex(0.3, 0.54, "Full drawn lines correspond to used components");
  l->DrawLatex(0.3, 0.51, "Dashed lines correspond to ignored components");
  l->DrawLatex(0.2, 0.47, "Each component has the form");
  l->DrawLatex(0.3, 0.42, "f_{n}(x; #Delta, #xi, #sigma') = "
	       "#int_{-#infty}^{+#infty}d#Delta' "
	       "landau(x; #Delta, #xi)gaus(x; #Delta', #sigma')");
  l->DrawLatex(0.2, 0.37, "The full function is given by");
  l->DrawLatex(0.3, 0.32, "f_{N}(x; #Delta, #xi, #sigma', #bf{a}) = "
	       "#sum_{i=1}^{N} a_{i} "
	       "f_{i}(x; #Delta_{i}, #xi_{i}, #sigma_{i}')");
  l->DrawLatex(0.3, 0.26, "#Delta_{i} = i (#Delta_{1} + #xi_{1} log(i))");
  l->DrawLatex(0.3, 0.23, "#xi_{i} = i #xi_{1}");
  l->DrawLatex(0.3, 0.20, "#sigma_{i} = #sqrt{i} #sigma_{1}");
  l->DrawLatex(0.3, 0.17, "#sigma_{n} #dot{=} 0");
  l->DrawLatex(0.3, 0.14, "#sigma'^{2} = #sigma^{2}_{n} + #sigma^{2}");
  l->DrawLatex(0.3, 0.11, "a_{1} = 1");
  c->Print(pname.Data(), "Title:Title page");

  ClearCanvas(c);

  fits->Draw(option);
  c->Print(pname.Data(), "Title:Fit overview");

  ClearCanvas(c);

  fits->Draw("rel");
  c->Print(pname.Data(), "Title:Relative parameter errors");

  Int_t nPad = 6;
  for (UShort_t d=1; d<=3; d++) { 
    UShort_t nQ = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nQ; q++) { 
      Char_t r = (q == 0 ? 'I' : 'O');

      TObjArray*  ra = fits->GetRingArray(d, r);
      if (!ra) continue;

      AliFMDCorrELossFit::ELossFit* fit = 0;
      TIter next(ra);
      Int_t i = 0;
      while ((fit = static_cast<AliFMDCorrELossFit::ELossFit*>(next()))) {
	if ((i % nPad) == 0) { 
	  ClearCanvas(c);
	  c->Divide(2,nPad/2,0,0);
	}
	TVirtualPad* p = c->cd((i % nPad) + 1);
	p->SetLogy();
	p->SetGridx();
	p->SetGridy();
	p->SetFillColor(kWhite);
	fit->Draw("comp");

	if ((i % nPad) == (nPad-1)) 
	  c->Print(pname.Data(), 
		   Form("Title:FMD%d%c page %d", d, r, (i/nPad)+1));
	i++;
      }
      if (i % nPad != 0) 
	c->Print(pname.Data(), 
		 Form("Title:FMD%d%c page %d", d, r, (i/nPad)+1));
    }
  }

  c->Print(Form("%s]", pname.Data()));
}
