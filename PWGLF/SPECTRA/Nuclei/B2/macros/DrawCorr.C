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

#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TF1.h>

#include "B2.h"

void DrawPair(TH1* hX, TH1* hY, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, const char* title="", const char* option="E", Int_t xMarker=kFullCircle, Int_t yMarker=kFullCircle, Int_t xColor=kBlue, Int_t yColor=kRed);

void DrawCorr(const TString& species="Deuteron", const TString& inputFile="corrections.root", const TString& tag="")
{
//
// Draw pt corrections for debugging
//
	Double_t xmin = 0;
	Double_t xmax = 3.5;
	
	const Int_t kNpart = 2;
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	const TString kPrefix[] = {"", "Anti"};
	
	// Response matrix
	
	TCanvas* c0 = new TCanvas(Form("%s.Unfolding", species.Data()), Form("Response matrix for (Anti)%ss", species.Data()));
	
	c0->Divide(2,1);
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		c0->cd(i+1);
		
		TH2D* hResponseMtx = (TH2D*)FindObj(finput, tag, Form("%s%s_Response_Matrix",kPrefix[i].Data(), species.Data()));
		hResponseMtx->SetAxisRange(xmin, xmax, "X");
		hResponseMtx->SetAxisRange(xmin, xmax, "Y");
		hResponseMtx->GetYaxis()->SetTitleOffset(1.3);
		hResponseMtx->DrawCopy("cont0colZ");
	}
	
	// Reconstruction efficiency
	
	TCanvas* c2 = new TCanvas(Form("%s.Efficiency",species.Data()), Form("Reconstruction Efficiency for (Anti)%ss",species.Data()));
	c2->Divide(2,2);
	
	TH1D* hEffVtxPt[kNpart];
	TH1D* hEffAccTrkPt[kNpart];
	TH1D* hEffTrigPt[kNpart];
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		hEffTrigPt[i] = (TH1D*)FindObj(finput, tag, kPrefix[i] + species + "_Eff_Trig_Pt");
		hEffVtxPt[i] = (TH1D*)FindObj(finput, tag, kPrefix[i] + species + "_Eff_Vtx_Pt");
		hEffAccTrkPt[i] = (TH1D*)FindObj(finput, tag, kPrefix[i] + species + "_Eff_AccTrk_Pt");
	}
	
	c2->cd(1);
	DrawPair(hEffTrigPt[0], hEffTrigPt[1], xmin, xmax, 0, 1.1);
	
	c2->cd(2);
	DrawPair(hEffVtxPt[0], hEffVtxPt[1], xmin, xmax, 0, 1.1);
	
	c2->cd(3);
	DrawPair(hEffAccTrkPt[0], hEffAccTrkPt[1], xmin, xmax, 0, 1.1);
	
	// Secondaries
	
	TCanvas* c3 = new TCanvas(Form("%s.Secondaries",species.Data()), Form("Fraction of secondaries for (Anti)%ss",species.Data()));
	
	c3->Divide(2,2);
	
	for(Int_t i=0; i<kNpart; ++i)
	{
		TH1D* hFracFdwnPt = (TH1D*)FindObj(finput, tag, kPrefix[i] + species + "_Frac_Fdwn_Pt");
		TH1D* hFracMatPt  = (TH1D*)FindObj(finput, tag, kPrefix[i] + species + "_Frac_Mat_Pt");
		
		TF1* fncFracFdwnPt = (TF1*)FindObj(finput, tag, kPrefix[i] + species + "_Frac_Fdwn_Fit_Pt");
		TF1* fncFracMatPt  = (TF1*)FindObj(finput, tag, kPrefix[i] + species + "_Frac_Mat_Fit_Pt");
		
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

void DrawPair(TH1* hX, TH1* hY, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax, const char* title, const char* option, Int_t xMarker, Int_t yMarker, Int_t xColor, Int_t yColor)
{
//
// Draw a pair of histograms in the current pad
//
	hX->SetTitle(title);
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
