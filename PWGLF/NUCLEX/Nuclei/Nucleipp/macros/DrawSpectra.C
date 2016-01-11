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

// Draw differential yields
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#endif

#include "B2.h"

void DrawSpectra(const TString& inputFile="spectra.root", const TString& tag="lhc10d", const TString& particle="Deuteron")
{
//
// Draw differential yields
//
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptLogy(1);
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	TCanvas* c0 = new TCanvas(Form("%s.Spectra", particle.Data()), Form("Spectra for %ss", particle.Data()));
	c0->Divide(2,1);
	
	// differential yield
	
	c0->cd(1);
	
	TGraphErrors* grDYieldPt = FindObj<TGraphErrors>(finput, tag, particle + "_DiffYield_Pt");
	
	grDYieldPt->SetMarkerStyle(kFullCircle);
	grDYieldPt->SetMarkerColor(kRed);
	grDYieldPt->SetLineColor(kRed);
	grDYieldPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	grDYieldPt->GetYaxis()->SetTitle("1/#it{N}_{ev} d^{2}#it{N}/(d#it{p}_{T}d#it{y}) (GeV^{-1}#it{c}^{2})");
	grDYieldPt->GetYaxis()->SetTitleOffset(1.30);
	grDYieldPt->Draw("zAP");
	
	// invariant differential yield
	
	c0->cd(2);
	
	TGraphErrors* grInvDYieldPt = FindObj<TGraphErrors>(finput, tag, particle + "_InvDiffYield_Pt");
	
	grInvDYieldPt->SetMarkerStyle(kFullCircle);
	grInvDYieldPt->SetMarkerColor(kRed);
	grInvDYieldPt->SetLineColor(kRed);
	grInvDYieldPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	grInvDYieldPt->GetYaxis()->SetTitle("1/(2#pi#it{N}_{ev}) d^{2}#it{N}/(#it{p}_{T}d#it{p}_{T}d#it{y}) (GeV^{-2}#it{c}^{3})");
	grInvDYieldPt->GetYaxis()->SetTitleOffset(1.30);
	grInvDYieldPt->Draw("zAP");
}
