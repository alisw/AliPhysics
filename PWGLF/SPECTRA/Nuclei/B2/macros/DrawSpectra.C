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

#include <Riostream.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "B2.h"

void DrawSpectra(const TString& inputFile="spectra.root", const TString& tag="lhc10d", const TString& particle="Deuteron")
{
//
// Draw differential yields
//
	gStyle->SetOptLogy(1);
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	TCanvas* c0 = new TCanvas(Form("%s.Spectra", particle.Data()), Form("Spectra for %ss", particle.Data()));
	c0->Divide(2,2);
	
	// differential yields
	
	c0->cd(1);
	
	TGraphErrors* grDYieldPt = (TGraphErrors*)FindObj(finput, tag, particle + "_DiffYield_Pt");
	TGraphErrors* grSysErrDYieldPt = (TGraphErrors*)FindObj(finput, tag, particle + "_SysErr_DiffYield_Pt");
	
	grSysErrDYieldPt->SetLineColor(kRed);
	grSysErrDYieldPt->SetFillStyle(0);
	grSysErrDYieldPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	grSysErrDYieldPt->GetYaxis()->SetTitle("#frac{1}{N_{inel}} #frac{d^{2}N}{dp_{T}dy} (GeV^{-1}c^{2})");
	grSysErrDYieldPt->Draw("A5");
	
	grDYieldPt->SetMarkerStyle(kFullCircle);
	grDYieldPt->SetMarkerColor(kRed);
	grDYieldPt->SetLineColor(kRed);
	grDYieldPt->Draw("P");
	
	// invariant differential yields
	
	c0->cd(2);
	
	TGraphErrors* grInvDYieldPt = (TGraphErrors*)FindObj(finput, tag, particle + "_InvDiffYield_Pt");
	TGraphErrors* grSysErrInvDYieldPt = (TGraphErrors*)FindObj(finput, tag, particle + "_SysErr_InvDiffYield_Pt");
	
	grSysErrInvDYieldPt->SetLineColor(kRed);
	grSysErrInvDYieldPt->SetFillStyle(0);
	grSysErrInvDYieldPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	grSysErrInvDYieldPt->GetYaxis()->SetTitle("#frac{1}{2#piN_{inel}} #frac{d^{2}N}{p_{T}dp_{T}dy} (GeV^{-2}c^{3})");
	grSysErrInvDYieldPt->Draw("A5");
	
	grInvDYieldPt->SetMarkerStyle(kFullCircle);
	grInvDYieldPt->SetMarkerColor(kRed);
	grInvDYieldPt->SetLineColor(kRed);
	grInvDYieldPt->Draw("P");
	
	// invariant differential cross section
	
	c0->cd(3);
	
	TGraphErrors* grInvDXsectPt = (TGraphErrors*)FindObj(finput, tag, particle + "_InvDiffXSection_Pt");
	TGraphErrors* grSysErrInvDXsectPt = (TGraphErrors*)FindObj(finput, tag, particle + "_SysErr_InvDiffXSection_Pt");
	
	grSysErrInvDXsectPt->SetLineColor(kRed);
	grSysErrInvDXsectPt->SetFillStyle(0);
	grSysErrInvDXsectPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	grSysErrInvDXsectPt->GetYaxis()->SetTitle("E #frac{d^{3}#sigma}{dp^{3}} (mb GeV^{-2}c^{3})");
	grSysErrInvDXsectPt->Draw("A5");
	
	grInvDXsectPt->SetMarkerStyle(kFullCircle);
	grInvDXsectPt->SetMarkerColor(kRed);
	grInvDXsectPt->SetLineColor(kRed);
	grInvDXsectPt->Draw("P");
}
