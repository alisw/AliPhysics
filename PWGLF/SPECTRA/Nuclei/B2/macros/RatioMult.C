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

// particle ratios as a function of multiplicity
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStyle.h>

#include "B2.h"
#include "Config.h"

void GetRatio(Double_t x, Double_t y, Double_t errx, Double_t erry, Double_t& r, Double_t& rerr);
void Draw(TGraph* gr, Int_t marker, Int_t color, const TString& xtitle, const TString& ytitle);

void RatioMult(const TString& pSpectra = "~/alice/output/Proton-lhc10d-Mult-Spectra.root",
               const TString& ptag     = "lhc10d",
               const TString& dSpectra = "~/alice/output/Deuteron-lhc10d-Mult-Spectra.root",
               const TString& dtag     = "lhc10d")
{
//
// particle ratios as a function of multiplicity (from fitting distributions)
//
	using namespace B2mult;
	
	const Int_t kNpart = 2;
	
	const TString kProton[kNpart] = { "Proton", "AntiProton" };
	const TString kDeuteron[kNpart] = { "Deuteron", "AntiDeuteron" };
	
	// open files
	
	TFile* fproton = new TFile(pSpectra.Data());
	if(fproton->IsZombie()) exit(1);
	
	TFile* fdeuteron = new TFile(dSpectra.Data());
	if(fdeuteron->IsZombie()) exit(1);
	
	// particle ratios
	TGraphErrors*  grRatio[kNpart]    = { new TGraphErrors(), new TGraphErrors() };
	TGraphErrors*  grMixRatio[kNpart] = { new TGraphErrors(), new TGraphErrors() };
	
	// model parameters
	TGraphErrors*  grProtondNdy[kNpart]    = { new TGraphErrors(), new TGraphErrors() };
	TGraphErrors*  grProtonN[kNpart]       = { new TGraphErrors(), new TGraphErrors() };
	TGraphErrors*  grProtonC[kNpart]       = { new TGraphErrors(), new TGraphErrors() };
	TGraphErrors*  grProtonChi2Ndf[kNpart] = { new TGraphErrors(), new TGraphErrors() };
	//TGraphErrors*  grProtonProb[kNpart]    = { new TGraphErrors(), new TGraphErrors() };
	
	TGraphErrors*  grDeuterondNdy[kNpart]    = { new TGraphErrors(), new TGraphErrors() };
	TGraphErrors*  grDeuteronN[kNpart]       = { new TGraphErrors(), new TGraphErrors() };
	TGraphErrors*  grDeuteronC[kNpart]       = { new TGraphErrors(), new TGraphErrors() };
	TGraphErrors*  grDeuteronChi2Ndf[kNpart] = { new TGraphErrors(), new TGraphErrors() };
	//TGraphErrors*  grDeuteronProb[kNpart]    = { new TGraphErrors(), new TGraphErrors() };
	
	for(Int_t i=0; i<kNmult; ++i)
	{
		TF1* fncProton[kNpart];
		TF1* fncDeuteron[kNpart];
		
		for(Int_t j=0; j<kNpart; ++j)
		{
			fncProton[j] = (TF1*)FindObj(fproton, ptag + "-" + kMultClass[i], kProton[j] + "_Fit_InvDiffYield_Pt");
			fncDeuteron[j] = (TF1*)FindObj(fdeuteron, dtag + "-" + kMultClass[i], kDeuteron[j] + "_Fit_InvDiffYield_Pt");
		}
		
		// integrated ratios
		
		Double_t dNdy[kNpart] = {fncProton[0]->GetParameter(0), fncDeuteron[0]->GetParameter(0)};
		Double_t dNdyErr[kNpart] = {fncProton[0]->GetParError(0), fncDeuteron[0]->GetParError(0)};
		
		Double_t dNdyBar[kNpart] = {fncProton[1]->GetParameter(0), fncDeuteron[1]->GetParameter(0)};
		Double_t dNdyBarErr[kNpart] = {fncProton[1]->GetParError(0), fncDeuteron[1]->GetParError(0)};
		
		// ratios
		Double_t ratio[kNpart], ratioErr[kNpart];
		for(Int_t j=0; j<kNpart; ++j)
		{
			GetRatio(dNdyBar[j], dNdy[j], dNdyBarErr[j], dNdyErr[j], ratio[j], ratioErr[j]);
			grRatio[j]->SetPoint(i, kKNOmult[i], ratio[j]);
			grRatio[j]->SetPointError(i, 0, ratioErr[j]);
		}
		
		// mixed ratios
		Double_t mixRatio[kNpart], mixRatioErr[kNpart];
		
		GetRatio(dNdy[1], dNdy[0], dNdyErr[1], dNdyErr[0], mixRatio[0], mixRatioErr[0]);
		GetRatio(dNdyBar[1], dNdyBar[0], dNdyBarErr[1], dNdyBarErr[0], mixRatio[1], mixRatioErr[1]);
		
		for(Int_t j=0; j<kNpart; ++j)
		{
			grMixRatio[j]->SetPoint(i, kKNOmult[i], mixRatio[j]);
			grMixRatio[j]->SetPointError(i, 0, mixRatioErr[j]);
		}
		
		// model parameters
		for(Int_t j=0; j<kNpart; ++j)
		{
			grProtondNdy[j]->SetPoint(i, kKNOmult[i], fncProton[j]->GetParameter(0));
			grProtondNdy[j]->SetPointError(i, 0, fncProton[j]->GetParError(0));
			
			grProtonN[j]->SetPoint(i, kKNOmult[i], fncProton[j]->GetParameter(1));
			grProtonN[j]->SetPointError(i, 0, fncProton[j]->GetParError(1));
			
			grProtonC[j]->SetPoint(i, kKNOmult[i], fncProton[j]->GetParameter(2));
			grProtonC[j]->SetPointError(i, 0, fncProton[j]->GetParError(2));
			
			grProtonChi2Ndf[j]->SetPoint(i, kKNOmult[i], fncProton[j]->GetChisquare()/fncProton[j]->GetNDF());
			//grProtonProb[j]->SetPoint(i, kKNOmult[i], fncProton[j]->GetProb());
			
			// deuteron
			grDeuterondNdy[j]->SetPoint(i, kKNOmult[i], fncDeuteron[j]->GetParameter(0));
			grDeuterondNdy[j]->SetPointError(i, 0, fncDeuteron[j]->GetParError(0));
			
			grDeuteronN[j]->SetPoint(i, kKNOmult[i], fncDeuteron[j]->GetParameter(1));
			grDeuteronN[j]->SetPointError(i, 0, fncDeuteron[j]->GetParError(1));
			
			grDeuteronC[j]->SetPoint(i, kKNOmult[i], fncDeuteron[j]->GetParameter(2));
			grDeuteronC[j]->SetPointError(i, 0, fncDeuteron[j]->GetParError(2));
			
			grDeuteronChi2Ndf[j]->SetPoint(i, kKNOmult[i], fncDeuteron[j]->GetChisquare()/fncDeuteron[j]->GetNDF());
			//grDeuteronProb[j]->SetPoint(i, kKNOmult[i], fncDeuteron[j]->GetProb());
		}
	}
	
	delete fproton;
	delete fdeuteron;
	
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
	
	TCanvas* c0 = new TCanvas("c0", "proton model parameters");
	c0->Divide(kNCol,2);
	
	for(Int_t j=0; j<kNpart; ++j)
	{
		c0->cd(kNCol*j+1);
		Draw(grProtondNdy[j], kFullCircle, kBlue, "z", "dN/dy");
		
		c0->cd(kNCol*j+2);
		Draw(grProtonN[j], kFullCircle, kBlue, "z", "n");
		
		c0->cd(kNCol*j+3);
		Draw(grProtonC[j], kFullCircle, kBlue, "z", "C");
		
		c0->cd(kNCol*j+4);
		Draw(grProtonChi2Ndf[j], kFullCircle, kBlue, "z", "#chi^{2}/ndf");
		
	/*	c0->cd(kNCol*j+5);
		Draw(grProtonProb[j], kFullCircle, kBlue, "z", "Prob");
	*/
	}
	
	TCanvas* c1 = new TCanvas("c1", "deuteron model parameters");
	c1->Divide(kNCol,2);
	
	for(Int_t j=0; j<kNpart; ++j)
	{
		c1->cd(kNCol*j+1);
		Draw(grDeuterondNdy[j], kFullCircle, kBlue, "z", "dN/dy");
		
		c1->cd(kNCol*j+2);
		Draw(grDeuteronN[j], kFullCircle, kBlue, "z", "n");
		
		c1->cd(kNCol*j+3);
		Draw(grDeuteronC[j], kFullCircle, kBlue, "z", "C");
		
		c1->cd(kNCol*j+4);
		Draw(grDeuteronChi2Ndf[j], kFullCircle, kBlue, "z", "#chi^{2}/ndf");
		
	/*	c1->cd(kNCol*j+5);
		Draw(grDeuteronProb[j], kFullCircle, kBlue, "z", "Prob");
	*/
	}
	
	TCanvas* c2 = new TCanvas("c2","particle ratios");
	c2->Divide(2,2);
	
	c2->cd(1);
	Draw(grRatio[0], kFullCircle, kRed, "z", "#bar{p}/p");
	
	c2->cd(2);
	Draw(grRatio[1], kFullCircle, kRed, "z", "#bar{d}/d");
	
	c2->cd(3);
	Draw(grMixRatio[0], kFullCircle, kRed, "z", "d/p");
	
	c2->cd(4);
	Draw(grMixRatio[1], kFullCircle, kRed, "z", "#bar{d}/#bar{p}");
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
