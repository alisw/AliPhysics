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

// Draw pt corrections for debugging
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TF1.h>
#endif

#include "B2.h"

void DrawPair(TH1* hX, TH1* hY, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, const TString& title="", const char* option="E", Int_t xMarker=kFullCircle, Int_t yMarker=kFullCircle, Int_t xColor=kBlue, Int_t yColor=kRed);

void DrawCorr(const TString& species="Deuteron", const TString& inputFile="corrections.root", const TString& tag="")
{
//
// Draw pt corrections for debugging
//
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	
	Double_t xmin = 0;
	Double_t xmax = 3.5;
	
	const Int_t kNpart = 2;
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	const TString kPrefix[] = {"", "Anti"};
	
	// Reconstruction efficiency
	
	TCanvas* c2 = new TCanvas(Form("%s.Efficiency",species.Data()), Form("Reconstruction Efficiency for (Anti)%ss",species.Data()));
	c2->Divide(2,2);
	
	TH1D* hEffTrigPt[kNpart];
	TH1D* hEffVtxPt[kNpart];
	TH1D* hEffAccPt[kNpart];
	TH1D* hEffAccTrkPt[kNpart];
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		hEffTrigPt[i] = FindObj<TH1D>(finput, tag, kPrefix[i] + species + "_Eff_Trig_Pt");
		hEffVtxPt[i] = FindObj<TH1D>(finput, tag, kPrefix[i] + species + "_Eff_Vtx_Pt");
		hEffAccPt[i] = FindObj<TH1D>(finput, tag, kPrefix[i] + species + "_Eff_Acc_Pt");
		hEffAccTrkPt[i] = FindObj<TH1D>(finput, tag, kPrefix[i] + species + "_Eff_AccTrk_Pt");
	}
	
	c2->cd(1);
	DrawPair(hEffTrigPt[0], hEffTrigPt[1], xmin, xmax, 0, 1.1);
	
	c2->cd(2);
	DrawPair(hEffVtxPt[0], hEffVtxPt[1], xmin, xmax, 0, 1.1);
	
	c2->cd(3);
	DrawPair(hEffAccPt[0], hEffAccPt[1], xmin, xmax, 0, 1.1);
	
	c2->cd(4);
	DrawPair(hEffAccTrkPt[0], hEffAccTrkPt[1], xmin, xmax, 0, 1.1);
	
	// Secondaries
	
	TCanvas* c3 = new TCanvas(Form("%s.Secondaries",species.Data()), Form("Fraction of secondaries for (Anti)%ss",species.Data()));
	
	c3->Divide(2,2);
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		TH1D* hFracFdwnPt = FindObj<TH1D>(finput, tag, kPrefix[i] + species + "_Frac_Fdwn_Pt");
		TH1D* hFracMatPt  = FindObj<TH1D>(finput, tag, kPrefix[i] + species + "_Frac_Mat_Pt");
		
		TF1* fncFracFdwnPt = FindObj<TF1>(finput, tag, kPrefix[i] + species + "_Frac_Fdwn_Fit_Pt");
		TF1* fncFracMatPt  = FindObj<TF1>(finput, tag, kPrefix[i] + species + "_Frac_Mat_Fit_Pt");
		
		c3->cd(2*i+1);
		hFracFdwnPt->SetAxisRange(xmin, xmax, "X");
		hFracFdwnPt->SetAxisRange(0., 0.6, "Y");
		hFracFdwnPt->GetYaxis()->SetTitleOffset(1.4);
		
		hFracFdwnPt->SetMarkerColor(kBlue);
		hFracFdwnPt->SetLineColor(kBlue);
		hFracFdwnPt->SetMarkerStyle(kFullCircle);
		hFracFdwnPt->DrawCopy("E");
		
		fncFracFdwnPt->SetLineColor(kRed);
		fncFracFdwnPt->SetLineWidth(1);
		fncFracFdwnPt->Draw("same");
		
		c3->cd(2*i+2);
		hFracMatPt->SetAxisRange(xmin, xmax, "X");
		hFracMatPt->SetAxisRange(0.,0.6, "Y");
		hFracMatPt->GetYaxis()->SetTitleOffset(1.4);
		
		hFracMatPt->SetMarkerColor(kBlue);
		hFracMatPt->SetLineColor(kBlue);
		hFracMatPt->SetMarkerStyle(kFullCircle);
		hFracMatPt->DrawCopy("E");
		
		fncFracMatPt->SetLineColor(kRed);
		fncFracMatPt->SetLineWidth(1);
		fncFracMatPt->Draw("same");
	}
}

void DrawPair(TH1* hX, TH1* hY, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, const TString& title, const char* option, Int_t xMarker, Int_t yMarker, Int_t xColor, Int_t yColor)
{
//
// Draw a pair of histograms in the current pad
//
	hX->SetTitle(title.Data());
	hX->SetAxisRange(xmin, xmax, "X");
	hX->SetAxisRange(ymin, ymax, "Y");
	hX->SetMarkerColor(xColor);
	hX->SetLineColor(xColor);
	hX->SetMarkerStyle(xMarker);
	
	hY->SetMarkerColor(yColor);
	hY->SetLineColor(yColor);
	hY->SetMarkerStyle(yMarker);
	
	hX->DrawCopy(option);
	hY->DrawCopy(Form("same%s",option));
}
