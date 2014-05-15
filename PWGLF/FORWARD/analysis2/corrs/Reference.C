#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TLatex.h>

struct Reference 
{
  static Double_t function(Double_t* x, Double_t* par)
  {
    const Double_t invSq2Pi = 1. / TMath::Sqrt(TMath::TwoPi());
    const Double_t mpShift  = -0.22278298;
    const Int_t    np       = 100;
    const Double_t sc       = 5;

    Double_t       sum      = 0;
    Double_t       x0       = x[0];
    Double_t       xi       = par[0];
    Double_t       mpc      = par[1] - mpShift * xi;
    Double_t       area     = par[2];
    Double_t       sigma    = par[3];
    Double_t       xLow     = x0 - sc * sigma;
    Double_t       xUpp     = x0 + sc * sigma;
    Double_t       step     = (xUpp - xLow) / np;
    
    for (Int_t i = 1; i <= np/2; i++) { 
      Double_t x1 = xLow + (i - .5) * step;
      Double_t x2 = xUpp - (i - .5) * step;
      
      sum += TMath::Landau(x1, mpc, xi, true) * TMath::Gaus(x0,x1,sigma);
      sum += TMath::Landau(x2, mpc, xi, true) * TMath::Gaus(x0,x2,sigma);
    }

    return area * step * sum * invSq2Pi / sigma;
  }

  static TF1* fit(TH1*      hist, 
		  Double_t* range, 
		  Double_t* guess, 
		  Double_t* lowLimits, 
		  Double_t* uppLimits, 
		  Double_t* params, 
		  Double_t* errors, 
		  Double_t& chi2, 
		  Int_t&    ndf)
  {
    TF1* func = new TF1(Form("func%s", hist->GetName()), &function, 
			range[0], range[1], 4);
    func->SetParameters(guess);
    func->SetParNames("#xi", "#Delta_{p}", "C", "#sigma");
    for (Int_t i = 0; i < 4; i++) 
      func->SetParLimits(i, lowLimits[i], uppLimits[i]);
      
    hist->Fit(func, "RB0");
    
    for (Int_t i = 0; i < 4; i++) {
      params[i] = func->GetParameter(i);
      errors[i] = func->GetParError(i);
    }
    chi2 = func->GetChisquare();
    ndf  = func->GetNDF();
    
    return func;
  }
  static Bool_t search(TF1*   f,    Double_t start, Double_t step,
		       Bool_t peak, Double_t& res)
  {
    const Int_t maxCalls = 10000;
    Double_t    lOld     = -2;
    Double_t    l        = peak ? -1 : -1e300;
    Int_t       i        = 0;
    Double_t    x        = 0;
    Double_t    cur      = start;
    
    while ((i < maxCalls) && 
	   TMath::Abs(l - lOld) > 1e-6 && 
	   TMath::Abs(step) > 1e-8) {
      i++;
      lOld = l;
      x    = cur + step;
      l    = f->Eval(x);
      if (!peak) l = TMath::Abs(l-res);
      // Printf("mode=%d x=%g step=%g l=%g lOld=%g", mode, x, step, l, lOld);
      
      if ((peak && l < lOld) || (!peak && l > lOld)) 
	step = -step/10; // Go the other way in smaller steps

      cur += step;
    }
    if (i >= maxCalls) return false;
    res = x;
    return true;
  }
  static Int_t find(TF1* f, Int_t iXi, Int_t iMpc, 
		    Double_t& maxX, Double_t& fwhm)
  {

    Double_t    xi       = f->GetParameter(iXi); // params[0];
    Double_t    mpc      = f->GetParameter(iMpc); // params[1];
    if (!search(f, mpc-0.1*xi,  0.05 * xi, true, maxX))  return false;
    Double_t    left     = f->Eval(maxX) / 2;
    Double_t    right    = left;
    if (!search(f, maxX-.5*xi, xi,         false, left))  return false;
    if (!search(f, maxX+xi,    -xi,        false, right)) return false;
    fwhm = right - left;

    return true;
  }
  static TH1* histo()
  {
    Int_t data[100] = {0,    0,  0,  0,  0,  0,  2,  6, 11, 18,
		       18,  55, 90,141,255,323,454,563,681,737,
		       821,796,832,720,637,558,519,460,357,291,
		       279,241,212,153,164,139,106, 95, 91, 76,
		       80,  80, 59, 58, 51, 30, 49, 23, 35, 28,
		       23,  22, 27, 27, 24, 20, 16, 17, 14, 20,
		       12,  12, 13, 10, 17,  7,  6, 12,  6, 12,
		       4,    9,  9, 10,  3,  4,  5,  2,  4,  1,
		       5,    5,  1,  7,  1,  6,  3,  3,  3,  4,
		       5,    4,  4,  2,  2,  7,  2,  4};
    TH1D *ret = new TH1D("snr","Signal-to-noise",100,0,100);
    for (Int_t i = 0; i < 100; i++) ret->Fill(i, data[i]);
    
    ret->SetFillColor(kRed+1);
    ret->SetFillStyle(3001);
    ret->SetXTitle("#Delta");
    ret->SetDirectory(0);

    return ret;
  }
  static void maxFwhm(TF1* f, Int_t iXi, Int_t iMpc)
  {
    printf("Find max and FWHM ...");
    Double_t maxX = 0;
    Double_t fwhm = 0;
    if (!find(f, iXi, iMpc, maxX, fwhm)) 
      Printf(" failed");
    else {
      Printf(" Max: %f FWHM: %f", maxX, fwhm);
      Double_t y = f->Eval(maxX);
      TGraph* maxG = new TGraph(1);
      maxG->SetPoint(0, maxX, y);
      maxG->SetMarkerStyle(20);
      maxG->SetMarkerColor(kBlue+2);
      maxG->Draw("p");

      TLatex* maxT = new TLatex(1.2*maxX, y, Form("Max: %f", maxX));
      maxT->SetTextAlign(11);
      maxT->Draw();

      TGraph* fwhmG = new TGraph(2);
      fwhmG->SetPoint(0, maxX-fwhm/2, y/2);
      fwhmG->SetPoint(1, maxX+fwhm/2, y/2);
      fwhmG->SetMarkerStyle(21);
      fwhmG->SetMarkerColor(kMagenta+2);
      fwhmG->Draw("pl");

      TLatex* fwhmT = new TLatex(1.2*(maxX+fwhm/2), y/2, 
				 Form("FWHM: %f", fwhm));
      fwhmT->SetTextAlign(11);
      fwhmT->Draw();
    }      

  }
  static void test()
  {
    TH1* hist = histo();
    
    Double_t range[] = { .3 * hist->GetMean(), 3 * hist->GetMean() };
    Double_t guess[] = { 1.8, 20., 5e4, 3.0 };
    Double_t low[]   = { 0.5,  5.,  1., 0.4 };
    Double_t upp[]   = { 5.0, 50., 1e6, 5.0 };
    Double_t chi2    = 0;
    Int_t    ndf     = 0;
    Double_t params[4];
    Double_t errors[4];
    
    printf("Fitting ...");
    TF1* f = fit(hist, range, guess, low, upp, params, errors, chi2, ndf);
    Printf(" done");
    f->Print();
    
    printf("Drawing ...");
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);
    hist->Draw();
    f->Draw("same");
    Printf(" done");

    maxFwhm(f, 0, 1);
  }
};
