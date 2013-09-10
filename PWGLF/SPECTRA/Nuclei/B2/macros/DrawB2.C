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

// Draw B2
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#endif

#include "B2.h"

void DrawB2(const TString& inputFile="b2.root", const TString& tag="test", const TString& prefix="", const TString& species="Deuteron")
{
//
// Draw B2
//
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	TString nucleus = prefix + species;
	TString suffix = "";
	if(prefix=="Anti") suffix="bar";
	
	TCanvas* c0 = new TCanvas(Form("%s.B2", nucleus.Data()), Form("Coalescence parameter for %s", nucleus.Data()));
	
	c0->Divide(2,1);
	
	// B2
	
	TGraphErrors* grB2Pt = FindObj<TGraphErrors>(finput, tag, Form("B2%s_Pt",suffix.Data()));
	
	c0->cd(1);
	gPad->SetLogy(1);
	
	grB2Pt->GetXaxis()->SetTitle("#it{p}_{T}/A (GeV/#it{c})");
	grB2Pt->GetYaxis()->SetTitle("#it{B}_{2} (GeV^{2}#it{c}^{-3})");
	grB2Pt->SetLineColor(kRed);
	grB2Pt->SetMarkerStyle(kFullCircle);
	grB2Pt->SetMarkerColor(kRed);
	grB2Pt->Draw("zAP");
	
	// homogeneity volume
	
	TGraphErrors* grR3Pt = FindObj<TGraphErrors>(finput, tag, Form("R3%s_Pt",suffix.Data()));
	
	c0->cd(2);
	
	grR3Pt->GetXaxis()->SetTitle("#it{p}_{T}/A (GeV/#it{c})");
	grR3Pt->GetYaxis()->SetTitle("#it{R}_{side}^{2} #it{R}_{long} (fm^{3})");
	grR3Pt->GetYaxis()->SetTitleOffset(1.3);
	grR3Pt->SetLineColor(kRed);
	grR3Pt->SetMarkerStyle(kFullCircle);
	grR3Pt->SetMarkerColor(kRed);
	grR3Pt->Draw("zAP");
}
