/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// d/p ratio
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TVirtualFitter.h>
#include <TPaveText.h>
#include <TAxis.h>
#endif

#include "B2.h"
#include "Config.h"

void GetRatio(Double_t x, Double_t y, Double_t errx, Double_t erry, Double_t& r, Double_t& rerr);
void Draw(TGraph* gr, Int_t marker, Int_t color, const TString& xtitle, const TString& ytitle);

void RatioMult(  const TString&  pSpectra   = "~/alice/output/Proton-lhc10d-Mult-Spectra.root"
               , const TString&  dSpectra   = "~/alice/output/Deuteron-lhc10d-hilow-Mult-Spectra.root"
               , const TString&  ptag       = "lhc10d"
               , const TString&  dtag       = "lhc10d"
               , const Int_t     nx         = B2Mult::kNmult
               , const TString*  xtag       = B2Mult::kMultTag
               , const Double_t* x          = B2Mult::kKNOmult
               , const Double_t* xerr       = B2Mult::kKNOmultErr
               , const TString&  xname      = B2Mult::kKNOmultName
               , const Bool_t    qTsallis   = 0
               , const TString&  outputfile = "~/alice/output/Particle-Ratios-lhc10d-Tsallis.root"
               , const TString&  otag       = "Tsallis")
{

//
// d/p ratio as a function of X
//
	const Int_t kNspec = 2;
	const Int_t kNpart = 2;
	
	const TString kParticle[kNspec][kNpart] = { {"Proton", "AntiProton"}, { "Deuteron", "AntiDeuteron"} };
	
	const Double_t xmin[kNspec] = {0., 0.};
	const Double_t xmax[kNspec] = {100, 100};
	
	const Double_t xminf[kNspec] = {0.4, 0.8};
	const Double_t xmaxf[kNspec] = {2., 2.2};
	
	TVirtualFitter::SetDefaultFitter("Minuit2");
	
	TFile* fproton = new TFile(pSpectra.Data(), "read");
	if(fproton->IsZombie()) exit(1);
	
	TFile* fdeuteron = new TFile(dSpectra.Data(), "read");
	if(fdeuteron->IsZombie()) exit(1);
	
	TFile* finput[kNspec] = { fproton, fdeuteron };
	TString tag[kNspec] = { ptag, dtag };
	
	// particle ratios
	TGraphErrors* grRatio[kNspec]           = { new TGraphErrors(), new TGraphErrors() };
	TGraphErrors* grMixRatio[kNspec]        = { new TGraphErrors(), new TGraphErrors() };
	
	// model parameters
	TGraphErrors* grP0[kNspec][kNpart]      = {{new TGraphErrors(), new TGraphErrors()},
	                                           {new TGraphErrors(), new TGraphErrors()}};
	TGraphErrors* grP1[kNspec][kNpart]      = {{new TGraphErrors(), new TGraphErrors()},
	                                           {new TGraphErrors(), new TGraphErrors()}};
	TGraphErrors* grP2[kNspec][kNpart]      = {{new TGraphErrors(), new TGraphErrors()},
	                                           {new TGraphErrors(), new TGraphErrors()}};
	
	TGraphErrors* grdNdy[kNspec][kNpart]    = {{new TGraphErrors(), new TGraphErrors()},
	                                           {new TGraphErrors(), new TGraphErrors()}};
	
	TGraphErrors* grMeanPt[kNspec][kNpart]  = {{new TGraphErrors(), new TGraphErrors()},
	                                           {new TGraphErrors(), new TGraphErrors()}};
	
	using std::cout;
	using std::endl;
	
	TGraphErrors* grDYieldPt[nx][kNspec][kNpart];
	TF1* fncTsallis[nx][kNspec][kNpart];
	
	for(Int_t i=0; i<nx; ++i)
	{
		Double_t dNdy[kNspec][kNpart];
		Double_t dNdyErr[kNspec][kNpart];
		
		cout << endl << xtag[i] << endl << endl;
		
		for(Int_t j=0; j<kNspec; ++j)
		{
			for(Int_t k=0; k<kNpart; ++k)
			{
				// data
				
				grDYieldPt[i][j][k] = FindObj<TGraphErrors>(finput[j], tag[j] + xtag[i], kParticle[j][k] + "_StatSystErr_DiffYield_Pt");
				
				// Tsallis distribution
				
				fncTsallis[i][j][k] = (qTsallis) ? QTsallisDYield(GetMass(kParticle[j][k]), kParticle[j][k] + "_Tsallis_DiffYield_Pt",xmin[j],xmax[j])
				                                 : TsallisDYield(GetMass(kParticle[j][k]), kParticle[j][k] + "_Tsallis_DiffYield_Pt",xmin[j],xmax[j]);
				
				// fit to data
				
				TFitResultPtr r = grDYieldPt[i][j][k]->Fit(fncTsallis[i][j][k], "RENSQ");
				Int_t status = r;
				
				Int_t n=0;
				while(status != 0 && n<7)
				{
					r = grDYieldPt[i][j][k]->Fit(fncTsallis[i][j][k], "RENSQ");
					status = r;
					++n;
					printf("\t*** %s, status: %d (%d)\n", kParticle[j][k].Data(), status,n);
				}
				
				// integral
				
				Double_t integral = fncTsallis[i][j][k]->Integral(xmin[j], xmax[j]);
				Double_t intErr   = fncTsallis[i][j][k]->IntegralError(xmin[j], xmax[j], r->GetParams(),  r->GetCovarianceMatrix().GetMatrixArray());
				
				// mean pt (dy cancels out)
				
				TF1* fncPtTsallis = (qTsallis) ? PtQTsallisDYield(GetMass(kParticle[j][k]), kParticle[j][k] + "_Tsallis_MeanPt",xmin[j],xmax[j])
				                            : PtTsallisDYield(GetMass(kParticle[j][k]), kParticle[j][k] + "_Tsallis_MeanPt",xmin[j],xmax[j]) ;
				
				Double_t par[]   = {  fncTsallis[i][j][k]->GetParameter(0)
				                    , fncTsallis[i][j][k]->GetParameter(1)
					            , fncTsallis[i][j][k]->GetParameter(2)
					           };
					
				Double_t parErr[] = {  fncTsallis[i][j][k]->GetParError(0)
					             , fncTsallis[i][j][k]->GetParError(1)
					             , fncTsallis[i][j][k]->GetParError(2)
					            };
				
				fncPtTsallis->SetParameters(par);
				fncPtTsallis->SetParErrors(parErr);
				
				Double_t intPtTsallis    = fncPtTsallis->Integral(xmin[j], xmax[j]);
				Double_t intPtTsallisErr = fncPtTsallis->IntegralError(xmin[j], xmax[j], r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray());
				
				Double_t meanPt, meanPtErr;
				GetRatio(intPtTsallis, integral, intPtTsallisErr, intErr, meanPt, meanPtErr);
				
				grMeanPt[j][k]->SetPoint(i, x[i], meanPt);
				grMeanPt[j][k]->SetPointError(i, xerr[i], meanPtErr);
				
				delete fncPtTsallis;
				
				// extrapolation fraction
				
				Double_t intFrac    = fncTsallis[i][j][k]->Integral(xminf[j], xmaxf[j]);
				Double_t intFracErr = fncTsallis[i][j][k]->IntegralError(xminf[j], xmaxf[j], r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray());
				
				Double_t extrFrac, extrFracErr;
				GetRatio(intFrac, integral, intFracErr, intErr, extrFrac, extrFracErr);
				
				// dN/dy
				
				dNdy[j][k]    = (qTsallis) ? fncTsallis[i][j][k]->GetParameter(0) : integral;
				dNdyErr[j][k] = (qTsallis) ? fncTsallis[i][j][k]->GetParError(0)  : intErr;
				
				grdNdy[j][k]->SetPoint(i, x[i], dNdy[j][k]);
				grdNdy[j][k]->SetPointError(i, xerr[i], dNdyErr[j][k]);
				
				// debug
				printf("\t*** %s, status: %d\tdN/dy: %g +/- %g\t extr: %g +/- %g\t<pt>: %g +/- %g\n", kParticle[j][k].Data(), status, dNdy[j][k], dNdyErr[j][k], 1.-extrFrac, extrFracErr, meanPt, meanPtErr);
			}
		}
		
		// ratios
		Double_t ratio[kNspec], ratioErr[kNspec];
		for(Int_t j=0; j<kNspec; ++j)
		{
			GetRatio(dNdy[j][1], dNdy[j][0], dNdyErr[j][1], dNdyErr[j][0], ratio[j], ratioErr[j]);
			
			grRatio[j]->SetPoint(i, x[i], ratio[j]);
			grRatio[j]->SetPointError(i, 0, ratioErr[j]);
		}
		
		// mixed ratios
		Double_t mixRatio[kNpart], mixRatioErr[kNpart];
		for(Int_t j=0; j<kNpart; ++j)
		{
			GetRatio(dNdy[1][j], dNdy[0][j], dNdyErr[1][j], dNdyErr[0][j], mixRatio[j], mixRatioErr[j]);
			
			grMixRatio[j]->SetPoint(i, x[i], mixRatio[j]);
			grMixRatio[j]->SetPointError(i, xerr[i], mixRatioErr[j]);
		}
		
		// model parameters
		for(Int_t j=0; j<kNspec; ++j)
		{
			for(Int_t k=0; k<kNpart; ++k)
			{
				grP0[j][k]->SetPoint(i, x[i], fncTsallis[i][j][k]->GetParameter(0));
				grP0[j][k]->SetPointError(i, xerr[i], fncTsallis[i][j][k]->GetParError(0));
				
				grP1[j][k]->SetPoint(i, x[i], fncTsallis[i][j][k]->GetParameter(1));
				grP1[j][k]->SetPointError(i, xerr[i], fncTsallis[i][j][k]->GetParError(1));
				
				grP2[j][k]->SetPoint(i, x[i], fncTsallis[i][j][k]->GetParameter(2));
				grP2[j][k]->SetPointError(i, xerr[i], fncTsallis[i][j][k]->GetParError(2));
			}
		}
	}
	
	// save
	
	TFile* foutput = new TFile(outputfile.Data(), "recreate");
	if (foutput->IsZombie()) exit(1);
	
	foutput->mkdir(otag.Data());
	foutput->cd(otag.Data());
	
	for(Int_t i=0; i<kNspec; ++i)
	{
		grRatio[i]->SetName(Form("%s%s_Ratio", kParticle[i][1].Data(), kParticle[i][0].Data()));
		grRatio[i]->Write();
		
		grMixRatio[i]->SetName(Form("%s%s_Ratio", kParticle[1][i].Data(), kParticle[0][i].Data()));
		grMixRatio[i]->Write();
	}
	
	for(Int_t i=0; i<kNspec; ++i)
	{
		for(Int_t j=0; j<kNpart; ++j)
		{
			grP0[i][j]->SetName(Form("%s_P0", kParticle[i][j].Data()));
			grP1[i][j]->SetName(Form("%s_P1", kParticle[i][j].Data()));
			grP2[i][j]->SetName(Form("%s_P2", kParticle[i][j].Data()));
			
			grP0[i][j]->Write();
			grP1[i][j]->Write();
			grP2[i][j]->Write();
			
			grdNdy[i][j]->SetName(Form("%s_dNdy", kParticle[i][j].Data()));
			grdNdy[i][j]->Write();
			
			grMeanPt[i][j]->SetName(Form("%s_MeanPt", kParticle[i][j].Data()));
			grMeanPt[i][j]->Write();
		}
	}
	
	// draw
	
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetOptStat(0);
	
	const Int_t kNCol = 3;
	
	TCanvas* c0[kNspec];
	
	for(Int_t i=0; i<kNspec; ++i)
	{
		c0[i] = new TCanvas(Form("c0.%s",kParticle[i][0].Data()), Form("%s model parameters",kParticle[i][0].Data()));
		c0[i]->Divide(kNCol,2);
		
		for(Int_t j=0; j<kNpart; ++j)
		{
			c0[i]->cd(kNCol*j+1);
			Draw(grP0[i][j], kFullCircle, kBlue, xname, fncTsallis[0][0][0]->GetParName(0));
			
			c0[i]->cd(kNCol*j+2);
			Draw(grP1[i][j], kFullCircle, kBlue, xname, fncTsallis[0][0][0]->GetParName(1));
			
			c0[i]->cd(kNCol*j+3);
			Draw(grP2[i][j], kFullCircle, kBlue, xname, fncTsallis[0][0][0]->GetParName(2));
		}
	}
	
	TCanvas* c1 = new TCanvas("c1", "Particle ratios");
	c1->Divide(2,2);
	
	c1->cd(1);
	Draw(grRatio[0], kFullCircle, kRed, xname, "#bar{p}/p");
	
	c1->cd(2);
	Draw(grRatio[1], kFullCircle, kRed, xname, "#bar{d}/d");
	
	c1->cd(3);
	Draw(grMixRatio[0], kFullCircle, kRed, xname, "d/p");
	
	c1->cd(4);
	Draw(grMixRatio[1], kFullCircle, kRed, xname, "#bar{d}/#bar{p}");
	
	// spectra
	
	TPaveText* label[nx][kNspec][kNpart];
	
	TCanvas* c2[kNspec][kNpart];
	
	for(Int_t j=0; j<kNspec; ++j)
	{
		for(Int_t k=0; k<kNpart; ++k)
		{
			c2[j][k] = new TCanvas(Form("c2.%s",kParticle[j][k].Data()), Form("%s data",kParticle[j][k].Data()));
			c2[j][k]->Divide(3,2);
			for(Int_t i=0; i<nx && i<6; ++i)
			{
				gPad->SetLogy(0);
				c2[j][k]->cd(i+1);
				Draw(grDYieldPt[i][j][k], kFullCircle, kBlue, "p_{T} (GeV/c)", "1/N_{ev} d^{2}N/dp_{T}dy (GeV^{-1})");
				fncTsallis[i][j][k]->SetLineWidth(1);
				fncTsallis[i][j][k]->Draw("same");
				
				label[i][j][k] = new TPaveText(0.5658525,0.7471751,0.814296,0.835452,"brNDC");
				label[i][j][k]->AddText(Form("%s = %.2f",xname.Data(),x[i]));
				label[i][j][k]->SetBorderSize(0);
				label[i][j][k]->SetFillColor(0);
				label[i][j][k]->SetTextSize(0.06);
				label[i][j][k]->Draw();
			}
		}
	}
	
	delete foutput;
	delete fproton;
	delete fdeuteron;
	
	// free memory
}

void GetRatio(Double_t x, Double_t y, Double_t errx, Double_t erry, Double_t& r, Double_t& rerr)
{
//
// get the ratio of x/y, its error and save in r and rerr
//
	if(x == 0 || y == 0)
	{
		r=0;
		rerr=1e+16;
		
		return;
	}
	
	r = x/y;
	rerr = r*TMath::Sqrt(TMath::Power(errx/x,2)+TMath::Power(erry/y,2));
}

void Draw(TGraph* gr, Int_t marker, Int_t color, const TString& xtitle, const TString& ytitle)
{
//
// draw a graph in current canvas
//
	gr->SetMarkerStyle(marker);
	gr->SetMarkerColor(color);
	gr->SetLineColor(color);
	gr->GetXaxis()->SetTitle(xtitle.Data());
	gr->GetYaxis()->SetTitle(ytitle.Data());
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->Draw("zAP");
}
