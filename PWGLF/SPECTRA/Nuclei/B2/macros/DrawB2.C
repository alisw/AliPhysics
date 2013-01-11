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

#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TStyle.h>

#include "B2.h"

void DrawB2(const TString& inputFile="b2.root", const TString& tag="test", const TString& prefix="", const TString& species="Deuteron")
{
//
// Draw B2
//
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
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	TString nucleus = prefix + species;
	TString suffix = "";
	if(prefix=="Anti") suffix="bar";
	
	TCanvas* c0 = new TCanvas(Form("%s.B2", nucleus.Data()), Form("Coalescence parameter for %s", nucleus.Data()));
	
	c0->Divide(2,1);
	
	// B2
	
	TGraphErrors* grB2Pt = (TGraphErrors*)FindObj(finput, tag, Form("B2%s_Pt",suffix.Data()));
	TGraphErrors* grSysErrB2Pt = (TGraphErrors*)FindObj(finput, tag, Form("B2%s_SysErr_Pt",suffix.Data()));
	
	c0->cd(1);
	gPad->SetLogy(1);
	
	grSysErrB2Pt->GetXaxis()->SetTitle("p_{T}/A (GeV/c)");
	grSysErrB2Pt->GetYaxis()->SetTitle("B_{2} (GeV^{2}c^{-3})");
	grSysErrB2Pt->SetLineColor(kRed);
	grSysErrB2Pt->SetFillStyle(1001);
	grSysErrB2Pt->SetFillColor(kRed-10);
	grSysErrB2Pt->Draw("A3");
	
	grB2Pt->SetLineColor(kRed);
	grB2Pt->SetMarkerStyle(kFullCircle);
	grB2Pt->SetMarkerColor(kRed);
	grB2Pt->Draw("P");
	
	// homogeneity volume
	
	TGraphErrors* grR3Pt = (TGraphErrors*)FindObj(finput, tag, Form("R3%s_Pt",suffix.Data()));
	TGraphErrors* grSysErrR3Pt = (TGraphErrors*)FindObj(finput, tag, Form("R3%s_SysErr_Pt",suffix.Data()));
	
	c0->cd(2);
	
	grSysErrR3Pt->GetXaxis()->SetTitle("p_{T}/A (GeV/c)");
	grSysErrR3Pt->GetYaxis()->SetTitle("R_{side}^{2} R_{long} (fm^{3})");
	grSysErrR3Pt->GetYaxis()->SetTitleOffset(1.3);
	grSysErrR3Pt->SetLineColor(kRed);
	grSysErrR3Pt->SetFillStyle(1001);
	grSysErrR3Pt->SetFillColor(kRed-10);
	grSysErrR3Pt->Draw("A3");
	
	grR3Pt->SetLineColor(kRed);
	grR3Pt->SetMarkerStyle(kFullCircle);
	grR3Pt->SetMarkerColor(kRed);
	grR3Pt->Draw("P");
}
