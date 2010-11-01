//____________________________________________________________________
//
// $Id$
//
// Script that contains a class to draw eloss from hits, versus ADC
// counts from digits, using the AliFMDInputHits class in the util library. 
//
// It draws the energy loss versus the p/(mq^2).  It can be overlayed
// with the Bethe-Bloc curve to show how the simulation behaves
// relative to the expected. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH1D.h>
#include <AliFMDHit.h>
#include <AliFMDDigit.h>
#include <AliFMDInput.h>
#include <AliFMDUShortMap.h>
#include <AliFMDFloatMap.h>
#include <AliFMDRecPoint.h>
#include <AliESDFMD.h>
#include <AliLog.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPolyMarker.h>
#include <TSpectrum.h>
#include <TList.h>
#include <TGraph.h>

Double_t landau(Double_t* xp, Double_t* pp)
{
  Double_t x     = xp[0];
  Double_t A     = pp[0];
  Double_t mpv   = pp[1];
  Double_t w     = pp[2];

  Double_t v     = mpv; //  + w * 0.22278;
  return A * TMath::Landau(x, v, w); 
}

Double_t foldLandau(Double_t* xp, Double_t* pp)
{
  Double_t x     = xp[0];
  Double_t A     = pp[0];
  Double_t mpv   = pp[1];
  Double_t w     = pp[2];
  Double_t B     = pp[3];
  Double_t sigma = pp[4];
  
  return A * (TMath::Landau(x,mpv,w) + B * TMath::Gaus(x, mpv, sigma));
}

/** @class DrawESD
    @brief Draw digit ADC versus Rec point mult
    @code 
    Root> .L Compile.C
    Root> Compile("DrawESD.C")
    Root> DrawESD c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawESD : public AliFMDInput
{
private:
  TH1D* fMult; // Histogram 
  const Double_t fCorr;
  TList fCleanup;
public:
  //__________________________________________________________________
  DrawESD(Int_t n=2000, Double_t mmin=-0.5, Double_t mmax=15.5) 
    : fCorr(1) // 0.68377 / 1.1)
  { 
    AddLoad(kESD);
    fMult = new TH1D("mult", "#DeltaE/#DeltaE_{MIP)", n, mmin, mmax);
    fMult->Sumw2();
    fMult->SetFillColor(kRed+1);
    fMult->SetFillStyle(3001);
    fMult->SetXTitle("#DeltaE/#DeltaE_{MIP}");
    fCleanup.Add(fMult);
  }
  ~DrawESD() 
  {
    fCleanup.Delete();
  }
  //__________________________________________________________________
  /** Begining of event
      @param ev Event number
      @return @c false on error */
  Bool_t Begin(Int_t ev) 
  {
    return AliFMDInput::Begin(ev);
  }
  //__________________________________________________________________
  Bool_t ProcessESD(UShort_t /* det */, 
		    Char_t   /* rng */, 
		    UShort_t /* sec */, 
		    UShort_t /* str */, 
		    Float_t  eta, 
		    Float_t  mult)
  {
    // Cache the energy loss 
    // if (mult > 0) 
    //   Info("ProcessESD", "FMD%d%c[%2d,%3d]=%f", det, rng, sec, str, mult);
    if (mult <= 0 || mult == AliESDFMD::kInvalidMult) return kTRUE;
    Double_t x = mult;
    if (!fESD->IsAngleCorrected()) {
      Double_t theta = 2 * TMath::ATan(TMath::Exp(-eta));
      Double_t corr  = TMath::Abs(TMath::Cos(theta));
      Double_t cmult = corr * mult;
      x = cmult;
    }
    if (x > 0.001) fMult->Fill(x);
    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessESDs()
  {
    if (!AliFMDInput::ProcessESDs()) return kFALSE;
    // Info("ProcessESDs", "ESD is %sangle corrected", 
    //      fESD->IsAngleCorrected() ? "" : "not ");
    return kTRUE;
  }
  //__________________________________________________________________
  TF1* FitPeak(Int_t n, TH1D* hist, Double_t min, Double_t& max)
  {
    if (TMath::Abs(max-min) < .25) return 0;
    std::cout << "Fit peak in range " << min << " to " << max << std::endl;
    TF1* l = new TF1(Form("l%02d", n), "landau", min, max);
    hist->GetXaxis()->SetRangeUser(0, 4);
    hist->Fit(l, "0Q", "", min, max);
    Double_t mpv   = l->GetParameter(1);
    Double_t empv  = l->GetParError(1);
    Double_t sigma = l->GetParameter(2);
    l->SetRange(mpv-empv, mpv+3*sigma);
    hist->Fit(l, "EMQ0", "", mpv-3*empv, mpv+3*sigma);
    std::cout << "Peak # " << n << " [" << min << "," << max << "]\n"
	      << " MPV: " << l->GetParameter(1) 
	      << " +/- "  << l->GetParError(1) 
	      << " Var: " << l->GetParameter(2) 
	      << " +/- "  << l->GetParError(2)
	      << " Chi^2/NDF: " << l->GetChisquare() / l->GetNDF() 
	      << std::endl;
    mpv   = l->GetParameter(1);
    sigma = l->GetParameter(2);
    min   = mpv - sigma * 2; // (n==1 ? 3 : 2);
    max   = mpv + sigma * 3;
    // l->SetRange(min, max);
    l->Draw("same");
    return  l;
  }
  //__________________________________________________________________
  void MaxInRange(TH1D* hist, Double_t min, Double_t& mean, Double_t& var)
  {
    hist->GetXaxis()->SetRangeUser(min, 4);
    mean = hist->GetMean();
    var  = hist->GetRMS();
  }

  //__________________________________________________________________
  const char* PrettyFloat(float x)
  {
    if (x == 0) return Form("%9.4f", x);
    float e = TMath::Floor(TMath::Log10(TMath::Abs(x)));
    if (TMath::Abs(e) < 4) {
      return Form("%9.4f", x);
    }
    float f = x / TMath::Power(10,e);
    return Form("%4.2f#times10^{%d}", f, int(e));
  }
  //__________________________________________________________________
  void ShowFit(Double_t x1, Double_t y1, const char* title, 
	       TF1* f, Double_t dx=0, Double_t dy=0.05)
  {
    Double_t x = x1, y = y1;
    TLatex* latex = new TLatex(x, y, title);
    latex->SetTextFont(132);
    latex->SetTextSize(0.8*dy);
    latex->SetNDC();
    latex->Draw();
    x -= dx;
    y -= dy;
    const Double_t eqDx=0.1;
    Double_t chi2 = f->GetChisquare();
    Int_t    ndf  = f->GetNDF();
    Double_t prob = f->GetProb();
    latex->DrawLatex(x, y, "#chi^{2}/NDF");
    latex->DrawLatex(x+eqDx, y, Form("= %7.4f/%3d=%5.2f (%7.3f%%)", 
				     chi2, ndf, chi2/ndf, 100*prob));
    Int_t     n = f->GetNpar();
    Double_t* p = f->GetParameters();
    Double_t* e = f->GetParErrors();
    for (int i = 0; i < n; i++) { 
      x -= dx;
      y -= dy;
      latex->DrawLatex(x, y, f->GetParName(i));
      latex->DrawLatex(x+eqDx, y, Form("= %s", PrettyFloat(p[i])));
      latex->DrawLatex(x+2.2*eqDx, y, Form("#pm %s", PrettyFloat(e[i])));
    }
  }      
  //__________________________________________________________________
  Bool_t Finish()
  {
    Info("Finish", "Will draw results");
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit(1111111);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    
    if (fMult->GetEntries() <= 0) return kFALSE;
    
    TCanvas* c = new TCanvas("c", "C", 1200, 800);
    fCleanup.Add(c);
    c->cd();
    c->SetLogy();
    c->SetTopMargin(0.05);
    c->SetRightMargin(0.05);
    c->SetFillColor(0);
    c->SetBorderMode(0);

    TLegend* leg = new TLegend(.1, .1, .4, .2);
    leg->SetFillColor(0);
    leg->SetBorderSize(1);

    DrawMult(c, leg);

    Double_t xmax=0, xmin=0, ymax=0;
    FindMinMax(xmin, xmax, ymax);

    FitLandau(xmin, xmax, ymax, leg);
    DrawResponse(xmax, ymax, leg);
    
    // TF1* f = FitMultiLandau(xmin, xmax, leg);
    // DrawLandaus(f);

    c->cd();
    leg->Draw();

    c->Modified();
    c->Update();
    c->cd();
    c->SaveAs("esd_eloss.png");

    fCleanup.Add(leg);

    return kTRUE;
    
  }
  //__________________________________________________________________
  /** 
   * Draw the multiplicity distribution
   * 
   */
  void DrawMult(TCanvas* /* c */, TLegend* leg) 
  {
    fMult->GetXaxis()->SetRangeUser(0.2,20);
    Double_t integral = fMult->Integral();
    Info("DrawMult", "Integral in range [0.2,20]=%f (%f entries)", 
	 integral, fMult->GetEntries());
    fMult->Scale(1. / integral);
    Info("DrawMult", "Integral in range [0.2,20]=%f (%f entries)", 
	 fMult->Integral(), fMult->GetEntries());
    Double_t max = 1.5 * fMult->GetMaximum();
    fMult->GetXaxis()->SetRangeUser(0,4);
    fMult->SetMaximum(max);
    fMult->SetStats(kFALSE);
    fMult->Draw("e same");
    leg->AddEntry(fMult, "Strip signal", "lf");

  }
  //__________________________________________________________________
  /** 
   * Find the minimum and maximum values in range
   * 
   * @param xmin 
   * @param xmax 
   * @param ymax 
   */
  void FindMinMax(Double_t& xmin, Double_t& xmax, Double_t& ymax)
  {
    fMult->GetXaxis()->SetRangeUser(0.1,4); 
    TSpectrum    s(10);
    Int_t        nPeak  = s.Search(fMult);
    fMult->GetXaxis()->SetRangeUser(0.1, 4);
    Int_t        bmax   = fMult->GetMaximumBin();
    xmax                = fMult->GetBinCenter(bmax);
    fMult->GetXaxis()->SetRangeUser(0.1, xmax);
    Double_t     bmin   = fMult->GetMinimumBin();
    xmin                = fMult->GetBinCenter(bmin);
    ymax                = fMult->GetBinContent(bmax);
    Info("Finish", "Simple peak finding found x_max=%f, x_min=%f y_max=%g", 
	 xmax, xmin, ymax);

    if (nPeak > 1) {
      TPolyMarker* pm = static_cast<TPolyMarker*>(fMult->GetListOfFunctions()
						  ->FindObject("TPolyMarker"));

      // Peaks are ordered by size 
      Double_t* peakX = pm->GetX();
      Double_t* peakY = pm->GetY();
      Double_t  xlmax = peakX[1];
      xmax            = peakX[0];
      ymax            = peakY[0];
      if (xmax < xlmax) { 
	xmax  = xlmax; 
	xlmax = peakX[0]; 
	ymax  = peakY[1];
      }

      fMult->GetXaxis()->SetRangeUser(xlmax, xmax);
      bmin  = fMult->GetMinimumBin();
      xmin  = fMult->GetBinCenter(bmin);
      Info("Finish", "Spectrum peak finding found x_max=%f, x_min=%f y_max=%g", 
	   xmax, xmin, ymax);
    }
  }
  //__________________________________________________________________
  /** 
   * Fit a landau and a landau+gaussian to the data
   * 
   * @param xmin  Minimum of peak range
   * @param xmax  Maximum of peak range
   * @param ymax  Y value in MIP peak
   * @param leg   Legend
   */
  void FitLandau(Double_t xmin, Double_t xmax, Double_t& ymax, TLegend* leg)
  {
    fMult->GetXaxis()->SetRangeUser(xmin, xmax+(xmax-xmin));
    Double_t mean = fMult->GetMean();
    Double_t var  = fMult->GetRMS();
    Info("Finish", "Peak range [%f,%f] mean=%f, var=%f", 
	 xmin, xmax+(xmax-xmin), mean, var);
    Double_t lowCut = mean-var;
    Double_t hiCut  = 4; // 2*mean;
    Info("Finish", "Low cut set to %f", lowCut);
    fMult->GetXaxis()->SetRangeUser(0, hiCut);
    
    TF1* pl = MakeLandau(lowCut, hiCut, xmax, var/2);
    fMult->Fit(pl, "NEM", "", lowCut, hiCut);
    ymax = pl->GetMaximum();
    Info("Finish", "Maximum of landau is at %f (y_max=%g)", 
	 pl->GetMaximumX(), ymax);

    TF1* gl  = MakeFoldLandau(lowCut, hiCut, pl, xmax, var/2);
    gl->SetLineColor(kRed+1);
    // gl->SetParLimits(1, xmax-var, xmax+var);
    fMult->Fit(gl, "NEM", "", lowCut, hiCut);
    TF1* l  = MakeLandau(lowCut, hiCut, xmax, var/2);
    l->SetLineColor(kGreen+1);
    l->SetParameters(gl->GetParameter(0),
		     gl->GetParameter(1),
		     gl->GetParameter(2));
    TF1* g  = new TF1("g", "gaus", lowCut, hiCut);
    g->SetParNames("A", "#mu", "#sigma");
    g->SetLineColor(kBlue+1);
    g->SetParameters(gl->GetParameter(3)*gl->GetParameter(0),
		     gl->GetParameter(1),
		     gl->GetParameter(2));
    fMult->GetXaxis()->SetRangeUser(0,4);
    fMult->DrawCopy("E");
    fMult->DrawCopy("HIST SAME");
    pl->Draw("same");
    gl->Draw("same");
    g->Draw("Same");
    l->Draw("Same");
    fCleanup.Add(pl);
    fCleanup.Add(l);
    fCleanup.Add(g);
    fCleanup.Add(gl);

    ShowFit(.5, .90, "Landau", pl);
    ShowFit(.5, .65, "Landau+Gaussian", gl);

    leg->AddEntry(pl, Form("Landau fit - #chi^{2}/NDF=%f", 
			   pl->GetChisquare()/pl->GetNDF()), "l");
    leg->AddEntry(gl, Form("Landau+Gaussian fit - #chi^{2}/NDF=%f", 
			   gl->GetChisquare()/gl->GetNDF()), "l");
    leg->AddEntry(l, "Landau part", "l");
    leg->AddEntry(g, "Gaussian part", "l");
  }
  
  //__________________________________________________________________
  /** 
   * Superimpose the response graph from the RPP 
   * 
   * @param ymax Y value of multiplicity spectra in the landau peak
   * @param leg  Legend
   */  
  void DrawResponse(Double_t xmax, Double_t ymax, TLegend* leg) 
  {
    TGraph*   resp = GetResp();
    // TGraph*   corr = GetCorr();
    Double_t* x    = resp->GetX();
    Double_t* y    = resp->GetY(); // [MeV cm^2/g]
    TGraph*   gr   = new TGraph(resp->GetN());
    gr->SetName(Form("%sCorr", resp->GetName()));
    gr->SetTitle(resp->GetTitle());
    gr->SetLineStyle(resp->GetLineStyle());
    gr->SetLineColor(kMagenta+1);
    gr->SetLineWidth(2);
    TGraph*   gr2   = new TGraph(resp->GetN());
    gr2->SetName(Form("%sCorr", resp->GetName()));
    gr2->SetTitle(resp->GetTitle());
    gr2->SetLineStyle(resp->GetLineStyle());
    gr2->SetLineColor(kCyan+1);
    gr2->SetLineWidth(2);
    // const Double_t rho = 2.33;   // [g/cm^3] 
    // const Double_t thk = 0.320;  // [cm]
    const Double_t mip = 1.664;  // [MeV cm^2/g]
    // const Double_t bgm = 3.4601; // beta*gamma of a MIP
    // Double_t  xs2 = corr->Eval(bgm); // [1]
    // Double_t  xss = 1.1;
    Double_t  xs  = 1/mip;
    Double_t xxmax = 0;
    Double_t yymax = 0;
    for (Int_t i = 0; i < gr->GetN(); i++) {
      if (y[i] > yymax) { 
	yymax = y[i];
	xxmax = xs * x[i];
      }
      gr->SetPoint(i, x[i] * xs, y[i] * ymax);
    }
    Info("DrawResponse", "Maximum at x=%f (xmax=%f)", xxmax, xmax);
    Double_t xs2 = xmax / xxmax;
    Info("DrawResponse", "Correction factor: %f", xs2);
    for (Int_t i = 0; i < gr->GetN(); i++) {
      gr2->SetPoint(i, x[i] * xs * xs2, y[i] * ymax);
    }
    gr->Draw("C same");
    gr2->Draw("C same");

    leg->AddEntry(gr, "Response", "l");
    leg->AddEntry(gr2, "Response", "l");
  }
  //__________________________________________________________________
  /** 
   * Fit sum of landaus to the multiplicity distribution
   * 
   * @param xmin Minimum of range 
   * @param xmax Maximum of range 
   * @param leg  Legend
   * 
   * @return Fitted function 
   */    
  TF1* FitMultiLandau(Double_t xmin, Double_t xmax, TLegend* leg)
  {
    fMult->GetXaxis()->SetRangeUser(xmin, xmax+(xmax-xmin));
    Double_t mean = fMult->GetMean();
    Double_t rms  = fMult->GetRMS();

    Double_t x1   = mean-rms/2; // .75; // .8;  // .65 / fCorr;
    Double_t x2   = mean+rms/2; // 1.3; // 1.7; // fCorr;
    TF1*     l1   = FitPeak(1, fMult, x1, x2);
    x2            = TMath::Max(mean+rms/2, x2);
    MaxInRange(fMult, x2, mean, rms);
    Double_t x3   = mean + rms;
    TF1*     l2   = FitPeak(2, fMult, x2, x3);
    TF1*     f    = 0;
    Double_t diff = 0;
    if (l2) {
      diff          = l2->GetParameter(1)-l1->GetParameter(1);
      f             = new TF1("user", "landau(0)+landau(3)", x1, x3);
      f->SetParNames("A_{1}", "Mpv_{1}", "#sigma_{1}",
		     "A_{2}", "Mpv_{2}", "#sigma_{2}",
		     "A_{3}", "Mpv_{3}", "#sigma_{3}");
      f->SetParameters(l1->GetParameter(0), 
		       l1->GetParameter(1), 
		       l1->GetParameter(2), 
		       l2 ? l2->GetParameter(0) : 0, 
		       l2 ? l2->GetParameter(1) : 0, 
		       l2 ? l2->GetParameter(2) : 0,
		       l2->GetParameter(0)/10, 
		       l2->GetParameter(1) + diff, 
		       l2->GetParameter(2));
    }
    else { 
      x3 = x2;
      f  = new TF1("usr", "landau", x1, x3);
    }
    
    std::cout << "Range: " << x1 << "-" << x3 << std::endl;
    
    fMult->GetXaxis()->SetRangeUser(0, 4);
    fMult->Fit(f, "NQR", "", x1, x3);
    fMult->Fit(f, "NMER", "E1", x1, x3);

    l1->SetLineColor(3);
    l1->SetLineWidth(2);
    l1->SetRange(0, 4);
    l1->Draw("same");
    if (l2) {
      l2->SetLineColor(4);
      l2->SetLineWidth(2);
      l2->SetRange(0, 4);
      l2->Draw("same");
    }
    f->SetLineWidth(2);
    f->SetRange(0, 4);
    f->Draw("same");

    fCleanup.Add(l1);
    if (l2) fCleanup.Add(l2);
    fCleanup.Add(f);


    leg->AddEntry(l1, "1 particle Landau", "l");
    if (l2) leg->AddEntry(l2, "2 particle Landau", "l");
    leg->AddEntry(f,  "1+2 particle Landau", "l");

    return f;
  }
  //__________________________________________________________________
  /** 
   * Draw landau functions in a separate canvas 
   * 
   * @param f Multi-landau function 
   */  
  void DrawLandaus(TF1* f)
  {
    TCanvas* c = new TCanvas("c2", "Landaus");
    // c->SetLogy();
    fMult->DrawClone("axis");
    f->Draw("same");
    Double_t* p1 = f->GetParameters();
    Double_t* p2 = &(p1[3]);
    TF1* ll1 = new TF1("ll1", "landau", 0, 4);
    ll1->SetParameters(p1);
    ll1->SetLineColor(3);
    ll1->Draw("same");
    TF1* ll2 = new TF1("ll2", "landau", 0, 4);
    ll2->SetParameters(p2);
    ll2->SetLineColor(4);
    ll2->Draw("same");

    Double_t diff = ll2->GetParameter(1)-ll1->GetParameter(1);
    Double_t y1   = fMult->GetMinimum() * 1.1;
    Double_t y2   = fMult->GetMaximum() * .9;
    Double_t xc1  = p1[1]-3*p1[2];
    Double_t xc2  = p2[1]-2*p2[2];
    Double_t xc3  = p2[1]-2*p2[2]+diff;
    TLine* c1 = new TLine(xc1, y1, xc1, y2);
    c1->Draw("same");
    TLine* c2 = new TLine(xc2, y1, xc2, y2);
    c2->Draw("same");
    TLine* c3 = new TLine(xc3, y1, xc3, y2);
    c3->Draw("same");

    TLegend* l = new TLegend(0.6, 0.6, .89, .89);
    l->AddEntry(ll1, "1 particle Landau", "l");
    l->AddEntry(ll2, "2 particle Landau", "l");
    l->AddEntry(f,  "1+2 particle Landau", "l");
    l->SetFillColor(0);
    l->Draw("same");

    c->Modified();
    c->Update();
    c->cd();
  }

  //__________________________________________________________________
  /** 
   * Get the response functin @f$ f(\Delta_p/x)@f$ from Review of
   * Particle Physics (fig. 27.7).  It is scaled to the value at MPV.
   *
   * @return Graph of response 
   */ 
  TGraph* GetResp()
  {
    static TGraph*  graph = 0;
    if (!graph) {
      graph = new TGraph;
      graph->SetName("si_resp");
      graph->SetTitle("f(#Delta/x) scaled to the MPV value ");
      graph->GetHistogram()->SetXTitle("#Delta/x (MeVcm^{2}/g)");
      graph->GetHistogram()->SetYTitle("f(#Delta/x)");
      graph->SetLineColor(kBlue+1);
      graph->SetLineWidth(2);
      graph->SetFillColor(kBlue+1);
      graph->SetMarkerStyle(21);
      graph->SetMarkerSize(0.6);
#if 0
      // Figure 27.7 or Review of Particle physics - Straggeling function in 
      // silicon of 500 MeV Pions, normalised to unity at the most probable 
      // value.   
      graph->SetPoint(0,0.808094,0.00377358);
      graph->SetPoint(1,0.860313,0.0566038);
      graph->SetPoint(2,0.891645,0.116981);
      graph->SetPoint(3,0.912533,0.181132);
      graph->SetPoint(4,0.928198,0.260377);
      graph->SetPoint(5,0.938642,0.320755);
      graph->SetPoint(6,0.954308,0.377358);
      graph->SetPoint(7,0.964752,0.433962);
      graph->SetPoint(8,0.975196,0.490566);
      graph->SetPoint(9,0.98564,0.550943);
      graph->SetPoint(10,0.996084,0.611321);
      graph->SetPoint(11,1.00653,0.667925);
      graph->SetPoint(12,1.02219,0.732075);
      graph->SetPoint(13,1.03264,0.784906);
      graph->SetPoint(14,1.0483,0.845283);
      graph->SetPoint(15,1.06397,0.901887);
      graph->SetPoint(16,1.09008,0.958491);
      graph->SetPoint(17,1.10574,0.984906);
      graph->SetPoint(18,1.13708,1);
      graph->SetPoint(19,1.13708,1);
      graph->SetPoint(20,1.15796,0.988679);
      graph->SetPoint(21,1.17363,0.966038);
      graph->SetPoint(22,1.19974,0.916981);
      graph->SetPoint(23,1.2154,0.89434);
      graph->SetPoint(24,1.23629,0.837736);
      graph->SetPoint(25,1.2624,0.784906);
      graph->SetPoint(26,1.28329,0.724528);
      graph->SetPoint(27,1.3094,0.664151);
      graph->SetPoint(28,1.32507,0.611321);
      graph->SetPoint(29,1.3564,0.550943);
      graph->SetPoint(30,1.41384,0.445283);
      graph->SetPoint(31,1.44517,0.392453);
      graph->SetPoint(32,1.48695,0.335849);
      graph->SetPoint(33,1.52872,0.286792);
      graph->SetPoint(34,1.58094,0.237736);
      graph->SetPoint(35,1.63838,0.196226);
      graph->SetPoint(36,1.68016,0.169811);
      graph->SetPoint(37,1.75326,0.135849);
      graph->SetPoint(38,1.81593,0.113208);
      graph->SetPoint(39,1.89426,0.0981132);
      graph->SetPoint(40,1.96214,0.0830189);
      graph->SetPoint(41,2.0718,0.0641509);
      graph->SetPoint(42,2.19191,0.0490566);
      graph->SetPoint(43,2.31723,0.0415094);
      graph->SetPoint(44,2.453,0.0301887);
      graph->SetPoint(45,2.53133,0.0264151);
      graph->SetPoint(46,2.57833,0.0264151);
#else
      graph->SetPoint(0,0.8115124,0.009771987);
      graph->SetPoint(1,0.9198646,0.228013);
      graph->SetPoint(2,0.996614,0.5895765);
      graph->SetPoint(3,1.041761,0.8241042);
      graph->SetPoint(4,1.059819,0.8794788);
      graph->SetPoint(5,1.077878,0.9348534);
      graph->SetPoint(6,1.100451,0.980456);
      graph->SetPoint(7,1.141084,0.9967427);
      graph->SetPoint(8,1.204289,0.9153094);
      graph->SetPoint(9,1.276524,0.742671);
      graph->SetPoint(10,1.402935,0.465798);
      graph->SetPoint(11,1.515801,0.3029316);
      graph->SetPoint(12,1.73702,0.1465798);
      graph->SetPoint(13,1.985327,0.08143322);
      graph->SetPoint(14,2.301354,0.04234528);
      graph->SetPoint(15,2.56772,0.02931596);
#endif
    }
    return graph;
  }
  //__________________________________________________________________
  /** 
   * Get the correction to Bethe-Bloc from Review of Particle Physics
   * (fig 27.8).
   *
   * @return correction graph
   */
  TGraph* GetCorr() 
  {
    static TGraph* graph = 0;
    if (!graph) {
      graph = new TGraph(14);
      graph->SetName("graph");
      graph->SetTitle("(#Delta_{p}/x)/(dE/dx)|_{mip} for 320#mu Si");
      graph->GetHistogram()->SetXTitle("#beta#gamma = p/m");
      graph->GetHistogram()->SetYTitle("#frac{#Delta_{p}/x)}{(dE/dx)|_{mip}}");
      graph->SetFillColor(1);
      graph->SetLineColor(7);
      graph->SetMarkerStyle(21);
      graph->SetMarkerSize(0.6);
      graph->SetPoint(0,1.196058,0.9944915);
      graph->SetPoint(1,1.28502,0.9411017);
      graph->SetPoint(2,1.484334,0.8559322);
      graph->SetPoint(3,1.984617,0.7491525);
      graph->SetPoint(4,2.658367,0.6983051);
      graph->SetPoint(5,3.780227,0.6779661);
      graph->SetPoint(6,4.997358,0.6741525);
      graph->SetPoint(7,8.611026,0.684322);
      graph->SetPoint(8,15.28296,0.6995763);
      graph->SetPoint(9,41.54516,0.7186441);
      graph->SetPoint(10,98.91461,0.7288136);
      graph->SetPoint(11,203.2734,0.7326271);
      graph->SetPoint(12,505.6421,0.7338983);
      graph->SetPoint(13,896.973,0.7338983);
    }
    return graph;
  }
  //__________________________________________________________________
  /** 
   * Make a Landau function object. 
   * 
   * @param min  Minimum of fit range 
   * @param max  Maximum of fit range
   * @param p    Peak position
   * @param v    Variance around peak
   * 
   * @return Landau function object
   */
  TF1* MakeLandau(Double_t min, Double_t max, Double_t p, Double_t v)
  {
    TF1* f = new TF1("l", "landau", min, max);
    f->SetParNames("A", "#delta", "#xi");
    f->SetParameters(1, p, p/10);
    if      (false) f->FixParameter(1,p);
    else if (false) f->SetParLimits(1, p-v, p+v);
    return f;
  }
  //__________________________________________________________________
  /** 
   * Make a Landau, folded with a gaussian, function object
   * 
   * @param min  Minimum of fit range 
   * @param max  Maximum of fit range
   * @param l    Seed Landau function object
   * @param p    Peak position
   * @param v    Variance around peak
   * 
   * @return Landau+Gaus function object
   */
  TF1* MakeFoldLandau(Double_t min, Double_t max, TF1* l, 
		      Double_t p, Double_t v)
  {
    TF1* f = new TF1("gl", &foldLandau, min, max, 5);
    f->SetParNames("A", "#delta", "#xi", "B", "#sigma");
    f->SetParameters(l->GetParameter(0), 
		     l->GetParameter(1),
		     l->GetParameter(2),
		     l->GetParameter(0)/1000,
		     l->GetParameter(2));
    f->SetParLimits(3, 1e-7, 1e-1);
    if      (false) f->FixParameter(1,p);
    else if (true)  f->SetParLimits(1, p-2*v, p+2*v);
    return f;
  }

  ClassDef(DrawESD,0);
  
};

//____________________________________________________________________
//
// EOF
//
