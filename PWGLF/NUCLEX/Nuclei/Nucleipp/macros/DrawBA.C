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

// Draw BA
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#endif

#include "B2.h"

void DrawBA(const TString& inputFile="b2.root", const TString& tag="test", const TString& nucleus="Deuteron")
{
//
// Draw B2, B3
//
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	Int_t A = (nucleus.Contains("Deuteron")) ? 2 : 3;
	
	TCanvas* c0 = new TCanvas(Form("%s.B%d", nucleus.Data(),A), Form("Coalescence parameter for %s", nucleus.Data()));
	
	if(A==2) c0->Divide(2,1);
	
	// BA
	
	TGraphErrors* grBAPt = FindObj<TGraphErrors>(finput, tag, nucleus + Form("_B%d_Pt",A));
	
	if(A==2) c0->cd(1);
	gPad->SetLogy(1);
	
	grBAPt->GetXaxis()->SetTitle("#it{p}_{T}/A (GeV/#it{c})");
	grBAPt->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}/#it{c}^{3})");
	if(A==3) grBAPt->GetYaxis()->SetTitle("#it{B}_{3} (GeV^{4}/#it{c}^{6})");
	grBAPt->SetLineColor(kRed);
	grBAPt->SetMarkerStyle(kFullCircle);
	grBAPt->SetMarkerColor(kRed);
	grBAPt->Draw("zAP");
	
	// homogeneity volume
	if(A==2)
	{
		TGraphErrors* grR3Pt = FindObj<TGraphErrors>(finput, tag, nucleus + "_R3_Pt");
		
		c0->cd(2);
		
		grR3Pt->GetXaxis()->SetTitle("#it{p}_{T}/A (GeV/#it{c})");
		grR3Pt->GetYaxis()->SetTitle("#it{R}_{side}^{2} #it{R}_{long} (fm^{3})");
		grR3Pt->GetYaxis()->SetTitleOffset(1.3);
		grR3Pt->SetLineColor(kRed);
		grR3Pt->SetMarkerStyle(kFullCircle);
		grR3Pt->SetMarkerColor(kRed);
		grR3Pt->Draw("zAP");
	}
}
