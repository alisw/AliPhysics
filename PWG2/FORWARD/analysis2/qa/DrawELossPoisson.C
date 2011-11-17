/**
 * @file   DrawELossPoisson.C
 * @author Christian Holm Christensen <cholm@dalsgaard.hehi.nbi.dk>
 * @date   Thu Jul  7 10:24:58 2011
 * 
 * @brief  A script to draw the Poisson vs Energy Loss correlation 
 * 
 *
 * @deprecated Use QATrender instead 
 * @ingroup pwg2_forward_scripts_qa
 * 
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
# include <TMath.h>
# include <TF1.h>
# include <TLine.h>
# include <TLatex.h>
# include <TLinearFitter.h>
# include <TStyle.h>
#else
class TList;
#endif

/** 
 * Draw the poisson @f$N_{ch}@f$ estimate against the @f$\Delta@f$
 * @f$N_{ch}@f$ estimate and do a regression line between the two 
 * 
 * @param p            List 
 * @param d            Detector
 * @param r            Ring
 * @param xmin         Minimum
 * @param xmax         Maximum
 * 
 * @return The regression coefficient 
 *
 * @deprecated Use QATrender instead 
 * @ingroup pwg2_forward_scripts_qa
 */
Double_t
DrawRingELossPoisson(TList* p, UShort_t d, Char_t r, 
		     Double_t xmin, Double_t xmax)
{
  if (!p) return 0;

  TList* ring = static_cast<TList*>(p->FindObject(Form("FMD%d%c",d,r)));
  if (!ring) { 
    Error("DrawELossPoisson", "List FMD%d%c not found in %s",d,r,p->GetName());
    return 0;
  }
  
  TH2* corr = static_cast<TH2D*>(ring->FindObject("elossVsPoisson"));
  if (!corr) { 
    Error("DrawRingELossPoisson", "Histogram esdEloss not found in FMD%d%c",
	  d, r);
    return 0;
  }
  TPad* pad = static_cast<TPad*>(gPad);
  pad->SetGridy();
  pad->SetGridx();
  pad->SetLogz();
  if (d == 3) { 
    pad->SetPad(pad->GetXlowNDC(), pad->GetYlowNDC(), .99, 
		 pad->GetYlowNDC()+pad->GetHNDC());
    pad->SetRightMargin(0.15);
  }
  pad->SetFillColor(0);
  if (xmax < 0) xmax = corr->GetXaxis()->GetXmax();
  corr->GetXaxis()->SetRangeUser(xmin,xmax);
  corr->GetYaxis()->SetRangeUser(xmin,xmax);
  corr->SetTitle(Form("FMD%d%c",d,r));
  // Info("", "Entries: %d, integral: %f", 
  //      int(corr->GetEntries()), corr->Integral()); 
  corr->Draw("colz");
  if (corr->GetEntries() <= 0) return 0;

  // Calculate the linear regression 
  // Double_t dx    = (xmax-xmin);
  Double_t rxy   = corr->GetCorrelationFactor();
  Double_t sx    = corr->GetRMS(1);
  Double_t sy    = corr->GetRMS(2);
  Double_t sx2   = sx*sx;
  Double_t sy2   = sy*sy;
  Double_t cxy   = corr->GetCovariance();
  Double_t mx    = corr->GetMean(1);
  Double_t my    = corr->GetMean(2);
  Double_t delta = 1;
#if 0
  Double_t beta  = rxy * sy / sx;
#else
  if (TMath::Abs(cxy) < 1e-6) return 0;
  Double_t beta  = ((sy2 - delta*sx2 + 
		     TMath::Sqrt(TMath::Power(sy2-delta*sx2,2) + 
				 4*delta*cxy*cxy))
		    / 2 / cxy);
#endif
  Double_t alpha = my - beta * mx;
  TF1* f = new TF1("f", "[0]+[1]*x", xmin, xmax);
  f->SetParameters(alpha, beta);
  f->SetLineWidth(1);
  f->SetLineStyle(1);
  f->Draw("same");

  TLine* l = new TLine(xmin,xmin,xmax,xmax);
  l->SetLineWidth(1);
  l->SetLineStyle(2);
  l->SetLineColor(kBlack);
  l->Draw();
  // corr->GetXaxis()->SetRangeUser(0, 2);

#if 0
  Info("", "FMD%d%c correlation coefficient: %9.5f%% "
       "line y = %f + %f * x", d, r, 100*rxy, alpha, beta);
#endif

  Double_t x = pad->GetLeftMargin()+.01;
  Double_t y = 1-pad->GetTopMargin()-.01;
  TLatex* ltx = new TLatex(x, y, Form("FMD%d%c", d, r));
  ltx->SetNDC();
  ltx->SetTextAlign(13);
  ltx->SetTextSize(0.08);
  ltx->Draw();
  y -= 0.12;
  ltx->SetTextSize(0.06);
  ltx->DrawLatex(x,y,"Deming regression: y=#alpha+#beta x");
  x += .02;
  y -= .06;
  ltx->SetTextSize(0.05);
  ltx->DrawLatex(x, y, Form("#alpha=%5.3f", alpha));
  y -= .06;
  ltx->DrawLatex(x, y, Form("#beta= %5.3f", beta));

  TLinearFitter* fitter = new TLinearFitter(1);
  fitter->SetFormula("1 ++ x");
  for (Int_t i = 1; i <= corr->GetNbinsX(); i++) { 
    x = corr->GetXaxis()->GetBinCenter(i);
    if (x < xmin || x > xmax) continue;
    for (Int_t j = 1; j <= corr->GetNbinsY(); j++) {
      y = corr->GetYaxis()->GetBinCenter(j);
      if (y < xmin || y > xmax) continue;
      Double_t c = corr->GetBinContent(i,j);
      if (c < .1) continue;
      fitter->AddPoint(&x, y, c);
    }
  }
  Bool_t robust = false;
  if (robust)
    fitter->EvalRobust();
  else 
    fitter->Eval();
#if 0
  for (Int_t i = 0; i < 2; i++) { 
    std::cout << i << "\t" 
	      << fitter->GetParName(i) << "\t"
	      << fitter->GetParameter(i) << "\t";
    if (!robust) 
      std::cout << fitter->GetParError(i); 
    std::cout << std::endl;
  }
  Double_t chi2 = fitter->GetChisquare();
  Int_t    ndf  = (fitter->GetNpoints() - 
		   fitter->GetNumberFreeParameters() );
  std::cout << "chi2/ndf: " << chi2 << '/' << ndf 
	    << '=' << chi2 / ndf << std::endl;
#endif

  // corr->Scale(1. / corr->GetMaximum());
  pad->cd();
  return corr->GetCorrelationFactor();
}

/** 
 * Draw the correlation between the Poisson @f$N_{ch}@f$ estimate
 * and the @f$\Delta@f$ @f$N_{ch}@f$ estimate and do a regression
 * line between the two for each ring
 * 
 * @param filename File to read
 * @param xmax     Minimum X
 * @param xmin     Maximum X 
 *
 * @deprecated Use QATrender instead 
 * @ingroup pwg2_forward_scripts_qa
 */
void
DrawELossPoisson(const char* filename="forward.root", 
		 const char* folder="ForwardResults",
		 Double_t xmax=-1,
		 Double_t xmin=-1)
{
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetTitleX(.4);
  // gStyle->SetTitleY(.1);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  
  TFile* file = TFile::Open(filename, "READ");
  if (!file) { 
    Error("DrawELossPoisson", "failed to open %s", filename);
    return;
  }

  TList* forward = static_cast<TList*>(file->Get(folder));
  if (!forward) { 
    Error("DrawELossPoisson", "List %s not found in %s", folder, filename);
    return;
  }

  TList* dc = static_cast<TList*>(forward->FindObject("fmdDensityCalculator"));
  if (!dc) { 
    Error("DrawELossPoisson", "List fmdDensityCalculator not found in Forward");
    return;
  }
  
  TCanvas* c = new TCanvas("elossVsPoisson", 
			   "N_ch from ELoss versus from Poisson", 900, 700);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetHighLightColor(0);
  c->Divide(3, 2, 0, 0);

  Double_t corrs[5];
  c->cd(1); corrs[0] = DrawRingELossPoisson(dc, 1, 'I', xmin, xmax);
  c->cd(2); corrs[1] = DrawRingELossPoisson(dc, 2, 'I', xmin, xmax);
  c->cd(5); corrs[2] = DrawRingELossPoisson(dc, 2, 'O', xmin, xmax);
  c->cd(3); corrs[3] = DrawRingELossPoisson(dc, 3, 'I', xmin, xmax);
  c->cd(6); corrs[4] = DrawRingELossPoisson(dc, 3, 'O', xmin, xmax);

  TVirtualPad* p = c->cd(4);
  p->SetTopMargin(0.05);
  p->SetRightMargin(0.10);
  p->SetLeftMargin(0.15);
  p->SetBottomMargin(0.15);
  p->SetFillColor(0);

  TH1D* hc = new TH1D("corrs", "Correlation factors", 5, .5, 5.5);
  hc->SetFillColor(kRed+1);
  hc->SetFillStyle(3001);
  hc->SetMinimum(0.0);
  hc->SetMaximum(1.3);
  hc->GetXaxis()->SetBinLabel(1,"FMD1i"); hc->SetBinContent(1,corrs[0]);
  hc->GetXaxis()->SetBinLabel(2,"FMD2i"); hc->SetBinContent(2,corrs[1]);
  hc->GetXaxis()->SetBinLabel(3,"FMD2o"); hc->SetBinContent(3,corrs[2]);
  hc->GetXaxis()->SetBinLabel(4,"FMD3i"); hc->SetBinContent(4,corrs[3]);
  hc->GetXaxis()->SetBinLabel(5,"FMD3o"); hc->SetBinContent(5,corrs[4]);
  hc->GetXaxis()->SetLabelSize(0.08);
  hc->GetYaxis()->SetTitle("r");
  hc->SetMarkerSize(1.5);
  hc->Draw("text hist");

  // TH2D* highCuts = static_cast<TH2D*>(dc->FindObject("highCuts"));
  // if (highCuts) highCuts->Draw("colz");
  c->cd();
  c->SaveAs("elossVsPoisson.png");
}
//
// EOF
//
