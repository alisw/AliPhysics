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

#include "B2.h"
#include "Config.h"

void GetRatio(Double_t x, Double_t y, Double_t errx, Double_t erry, Double_t& r, Double_t& rerr);
void Draw(TGraph* gr, Int_t marker, Int_t color, const TString& xtitle, const TString& ytitle);

void RatioMult(  const TString&  pSpectra   = "~/alice/output/Proton-lhc10d-Mult-Spectra.root"
               , const TString&  dSpectra   = "~/alice/output/Deuteron-lhc10d-Mult-Spectra.root"
               , const TString&  ptag       = "lhc10d"
               , const TString&  dtag       = "lhc10d"
               , const Int_t     nx         = B2mult::kNmult
               , const TString*  xtag       = B2mult::kMultTag
               , const Double_t* x          = B2mult::kKNOmult
               , const Double_t* xerr       = B2mult::kKNOmultErr
               , const TString&  xname      = B2mult::kKNOmultName
               , const Bool_t    tsallis    = 1
               , const TString&  outputfile = "~/alice/output/Ratio-Mult.root"
               , const TString&  otag       = "Tsallis")
{
//
// d/p ratio as a function of X
//
	const Int_t kNspec = 2;
	const Int_t kNpart = 2;
	
	const TString kParticle[kNspec][kNpart] = { {"Proton", "AntiProton"}, { "Deuteron", "AntiDeuteron"} };
	
	//TVirtualFitter::SetDefaultFitter("Minuit2");
	
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
	TGraphErrors* grChi2Ndf[kNspec][kNpart] = {{new TGraphErrors(), new TGraphErrors()},
	                                           {new TGraphErrors(), new TGraphErrors()}};
	
	TGraphErrors* grdNdy[kNspec][kNpart]    = {{new TGraphErrors(), new TGraphErrors()},
	                                           {new TGraphErrors(), new TGraphErrors()}};
	
	TGraphErrors* grDYieldPt[nx][kNspec][kNpart];
	TF1* fncTsallis[nx][kNspec][kNpart];
	
	for(Int_t i=0; i<nx; ++i)
	{
		Double_t dNdy[kNspec][kNpart];
		Double_t dNdyErr[kNspec][kNpart];
		
		for(Int_t j=0; j<kNspec; ++j)
		{
			for(Int_t k=0; k<kNpart; ++k)
			{
				grDYieldPt[i][j][k] = (TGraphErrors*)FindObj(finput[j], tag[j] + "-" + xtag[i], kParticle[j][k] + "_SysErr_DiffYield_Pt");
				
				if(tsallis)
				{
					fncTsallis[i][j][k] = TsallisDYield(GetMass(kParticle[j][k]), kParticle[j][k] + "_Tsallis_DiffYield_Pt",0,10);
				}
				else
				{
					fncTsallis[i][j][k] = TsallisParetoDYield(GetMass(kParticle[j][k]), kParticle[j][k] + "_Tsallis_DiffYield_Pt",0,10);
				}
				
				TFitResultPtr r = grDYieldPt[i][j][k]->Fit(fncTsallis[i][j][k], "RENSQ");
				
				if(tsallis)
				{
					dNdy[j][k] = fncTsallis[i][j][k]->Integral(0,100);
					dNdyErr[j][k] = fncTsallis[i][j][k]->IntegralError(0,100,r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray());
				}
				else
				{
					dNdy[j][k] = fncTsallis[i][j][k]->GetParameter(0);
					dNdyErr[j][k] = fncTsallis[i][j][k]->GetParError(0);
				}
				
				Int_t status = r;
				printf("status: %d\tdN/dy: %g +/- %g\n",status,dNdy[j][k],dNdyErr[j][k]);
				
				grdNdy[j][k]->SetPoint(i, x[i], dNdy[j][k]);
				grdNdy[j][k]->SetPointError(i, xerr[i], dNdyErr[j][k]);
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
				
				grChi2Ndf[j][k]->SetPoint(i, x[i], fncTsallis[i][j][k]->GetChisquare()/fncTsallis[i][j][k]->GetNDF());
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
			grChi2Ndf[i][j]->SetName(Form("%s_Chi2Ndf", kParticle[i][j].Data()));
			
			grdNdy[i][j]->Write();
			grChi2Ndf[i][j]->Write();
		}
	}
	
	// draw
	
	TStyle* st = new TStyle();
	
	st->SetPadTickX(1);
	st->SetPadTickY(1);
	st->SetPadGridX(1);
	st->SetPadGridY(1);
	st->SetCanvasColor(0);
	st->SetFrameBorderMode(0);
	st->SetFrameFillColor(0);
	st->SetLabelFont(62,"XYZ");
	st->SetTitleFont(62,"XYZ");
	
	st->cd();
	gROOT->ForceStyle();
	
	const Int_t kNCol = 4;
	
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
			
			c0[i]->cd(kNCol*j+4);
			Draw(grChi2Ndf[i][j], kFullCircle, kBlue, xname, "#chi^{2}/ndf");
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
	
	TCanvas* c2[kNspec];
	
	for(Int_t j=0; j<kNspec; ++j)
	{
		c2[j] = new TCanvas(Form("c2.%s",kParticle[j][0].Data()), Form("%s data",kParticle[j][0].Data()));
		c2[j]->Divide(nx,2);
		
		for(Int_t k=0; k<kNpart; ++k)
		{
			for(Int_t i=0; i<nx; ++i)
			{
				gPad->SetLogy(0);
				c2[j]->cd(nx*k+i+1);
				Draw(grDYieldPt[i][j][k], kFullCircle, kBlue, "p_{T} (GeV/c)", "DYield");
				fncTsallis[i][j][k]->Draw("same");
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
	gr->Draw("AP");
}
