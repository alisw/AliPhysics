#include <Riostream.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>

Double_t bw(Double_t *x, Double_t *par)
{
	// Fit parameters:
	// par[0] = normalization factor
	// par[1] = peak position
	// par[2] = FWHM
	
	return par[0] * TMath::BreitWigner(x[0], par[1], par[2]);
}

Double_t bwgaus(Double_t *x, Double_t *par) 
{
	// Fit parameters:
	// par[0] = global normalization constant
	// par[1] = BW peak position
	// par[2] = FWHM
	// par[3] = sigma of convoluted gaussian
	
	// Numeric constants
    Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
	
	// Control constants
	Double_t nsteps = 100.0;      // number of convolution steps
	Double_t nsigma =   3.0;      // convolution extends to +-sc Gaussian sigmas
	
	// Range of convolution integral
    Double_t x1 = x[0] - nsigma * par[3];
	Double_t x2 = x[0] + nsigma * par[3];
	Double_t dx = (x2 - x1) / nsteps;
	
	// Variables
	Double_t i, xx, fbw, fgaus, sum = 0.0;
	
	// Convolution integral of Breit-Wigner and Gaussian by sum
	for(i = 1.0; i <= 0.5 * nsteps; i++) {

		xx = x1 + (i - 0.5) * dx;
		fbw = TMath::BreitWigner(xx, par[1], par[2]);
		fgaus = TMath::Gaus(x[0], xx, par[3]);
        sum += fbw * fgaus;

		xx = x2 - (i - 0.5) * dx;
		fbw = TMath::BreitWigner(xx, par[1], par[2]);
		fgaus = TMath::Gaus(x[0], xx, par[3]);
        sum += fbw * fgaus;
		
      }

      return (par[0] * dx * sum * invsq2pi / par[3]);
}

Double_t bwgaus_line(Double_t *x, Double_t *par)
{
	Double_t bwg = bwgaus(x, par);
	return bwg + par[4] + par[5] * x[0];
}

void view 
(Int_t rebin = 8, Double_t mult = 4.0, Double_t fitmin = 0.7, Double_t fitmax = 0.8,
 const char *filename = "kstar.invmass.root") 
{
	TFile *f = TFile::Open(filename);
	
	TH1D * hsign   = (TH1D*)f->Get("h_Pi(+)_K(-)");
	TH1D * hsign_2 = (TH1D*)f->Get("h_Pi(-)_K(+)");
	
	TH1D * hmix    = (TH1D*)f->Get("hmix_Pi(+)_K(-)");
	TH1D * hmix_2  = (TH1D*)f->Get("hmix_Pi(-)_K(+)");
	TH1D * hmix_3  = (TH1D*)f->Get("hmix_K(+)_Pi(-)");
	TH1D * hmix_4  = (TH1D*)f->Get("hmix_K(-)_Pi(+)");
	
	TH1D * htrue   = (TH1D*)f->Get("h_Pi(+)_K(-)_true");
	TH1D * htrue_2 = (TH1D*)f->Get("h_Pi(-)_K(+)_true");
	
	hsign->Rebin(rebin);
	hsign_2->Rebin(rebin);
	
	htrue->Rebin(rebin);
	htrue_2->Rebin(rebin);
	
	hmix->Rebin(rebin);
	hmix_2->Rebin(rebin);
	hmix_3->Rebin(rebin);
	hmix_4->Rebin(rebin);
		
	hsign->Add(hsign_2);
	htrue->Add(htrue_2);
	
	hmix->Add(hmix_2);
	hmix->Add(hmix_3);
	hmix->Add(hmix_4);
	
	hsign->Sumw2();
	hmix->Sumw2();
	htrue->Sumw2();
	
	cout << "# signal  entries: " << hsign->Integral() << endl;
	cout << "# backgr. entries: " << hmix->Integral() << endl;
	cout << "# true    pairs:   " << htrue->Integral() << endl;
	
	// mean and width of lambda*
	Double_t mean = 0.896;
	Double_t gamma = 0.05;
	
	// find bin limits for fit
	Int_t b1 = hsign->GetXaxis()->FindFixBin(fitmin);
	Int_t b2 = hsign->GetXaxis()->FindFixBin(fitmax);
	hmix->Scale( hsign->Integral(b1, b2) / hmix->Integral(b1, b2) );
	
	// drawing options
	gStyle->SetOptStat(0);
	gROOT->SetStyle("Plain");
	
	// draw signal and background together
	TCanvas *c1 = new TCanvas("c1", "", 0, 0, 640, 480);
	TH1D *hsign1 = (TH1D*)hsign->Clone();
	TH1D *hmix1 = (TH1D*)hmix->Clone();
	hsign1->GetXaxis()->SetTitle("Inv. mass (GeV/c^{2})");
	hsign1->GetYaxis()->SetTitle("counts");	
	hsign1->SetTitle("");
	hsign1->SetMarkerStyle(8);
	hsign1->SetMarkerSize(0.8);
	hsign1->SetStats(0);
	hsign1->Draw("PE");
	//hmix1->SetLineColor(kRed);
	hmix1->SetMarkerStyle(4);
	hmix1->SetStats(0);
	hmix1->GetXaxis()->SetTitle("Inv. mass (GeV/c^{2})");
	hmix1->GetYaxis()->SetTitle("counts");
	hmix1->Draw("CESAME");
	c1->Update();
	
	// draw subtraction and fit with BW+gaus
	TCanvas *c2 = new TCanvas("c2", "", 50, 50, 640, 480);
	TH1D *hdiff = (TH1D*)hsign1->Clone();
	hdiff->Add(hmix, -1.0);
	hdiff->GetXaxis()->SetRangeUser(mean - mult*gamma, mean + mult*gamma);
//	TF1 *fcn = new TF1("fcn", bwgaus_line, mean - mult*gamma, mean + mult*gamma, 6);
	TF1 *fcn = new TF1("fcn", bwgaus, mean - mult*gamma, mean + mult*gamma, 4);
	fcn->SetParameter(0, hdiff->GetMaximum() * 0.6);
	fcn->SetParameter(1, mean);
	fcn->SetParameter(2, gamma);
	fcn->SetParameter(3, 0.001);
//	fcn->SetParameter(4, 20);
//	fcn->SetParameter(5, -50);
//	fcn->SetParNames("Constant", "BW_peak", "BW_gamma", "Gaus_sigma", "Line_A", "Line_B");
	fcn->SetParNames("Constant", "BW_peak", "BW_gamma", "Gaus_sigma");
	fcn->SetLineColor(kGreen);
//	hdiff->Fit(fcn, "RE");
	hdiff->SetStats(0);
	hdiff->GetXaxis()->SetTitle("Inv. mass (GeV/c^{2})");
	hdiff->GetYaxis()->SetTitle("counts");
	hdiff->Draw("PE");
	cout << "Fit results: " << endl;
	cout << "Peak center: " << fcn->GetParameter(1) << " +/- " << fcn->GetParError(1) << endl;
	cout << "Peak width : " << fcn->GetParameter(2) << " +/- " << fcn->GetParError(2) << endl;
	cout << "Gaus sigma : " << fcn->GetParameter(3) << " +/- " << fcn->GetParError(3) << endl;
	c2->Update();
	
	// draw true pairs
	TCanvas *c3 = new TCanvas("c3", "", 100, 100, 640, 480);
	TH1D *htrue1 = (TH1D*)htrue->Clone();
	htrue1->SetMarkerStyle(21);
	htrue1->SetMarkerColor(kRed);
	htrue1->GetXaxis()->SetRangeUser(mean - mult*gamma, mean + mult*gamma);
	htrue1->SetStats(0);
	htrue1->GetXaxis()->SetTitle("Inv. mass (GeV/c^{2})");
	htrue1->GetYaxis()->SetTitle("counts");
	htrue1->Draw("PE");
	c3->Update();
	
	// draw true pairs with diff
	TCanvas *c4 = new TCanvas("c4", "", 150, 150, 640, 480);
	hdiff->Draw("PE");
	htrue1->SetMarkerStyle(4);
	htrue1->SetMarkerColor(1);
	htrue1->Draw("PEsame");
	c4->Update();
	
	// print true and reconstructed pairs
	b1 = htrue->GetXaxis()->FindFixBin(mean - 1.5*gamma);
	b2 = htrue->GetXaxis()->FindFixBin(mean + 1.5*gamma);
	TF1 *fbw = new TF1("fbw", bw, mean - mult*gamma, mean + mult*gamma, 3);
	fbw->SetParameter(0, fcn->GetParameter(0));
	fbw->SetParameter(1, fcn->GetParameter(1));
	fbw->SetParameter(2, fcn->GetParameter(2));
	Double_t S = htrue->Integral(b1, b2);
	Double_t B = hsign->Integral(b1, b2);
	cout << "True pairs:        " << S << endl;
	cout << "Reconstructed:     " << fbw->Integral(mean - 1.5*gamma, mean + 1.5*gamma) / htrue->GetBinWidth(1) << endl;
	cout << "Background pairs:  " << B << endl;
	cout << "S/B (+/- 1.5 G):   " << S / B << endl;
	cout << "sign. (+/- 1.5 G): " << S / TMath::Sqrt(S + B) << endl;
	
	// save pictures
	Text_t answer;
	cout << "Save canvases as eps files (y/n)? ";
	cin >> answer;
	if (answer == 'y' || answer == 'Y') {
		c1->SaveAs("lambda_sign_bg.eps");
		c2->SaveAs("lambda_diff.eps");
		c3->SaveAs("lambda_true.eps");
		c4->SaveAs("lambda_comparison.eps");
	}
}
	
