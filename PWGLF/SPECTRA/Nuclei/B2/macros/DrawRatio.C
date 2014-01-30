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

// Draw antiparticle/particle ratio
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#endif

#include "B2.h"

void DrawRatio(const TString& inputFile="pt.root", const TString& tag="test", const TString& species="Deuteron")
{
//
// Draw antiparticle/particle ratio
//
	Double_t xmin = 0;
	Double_t xmax = 3.5;
	Double_t ymin = 0.73;
	Double_t ymax = 1.14;
	
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	// antiparticle/particle ratio
	
	TString ytitle = "Negative/Positive";
	if(species=="Proton")        ytitle = "#bar{p} / p";
	else if(species=="Deuteron") ytitle = "#bar{d} / d";
	else if(species=="Triton")   ytitle = "#bar{t} / t";
	else if(species=="He3")      ytitle = "{}^{3}#bar{He} / {}^{3}He";
	else if(species=="Alpha")    ytitle = "#bar{#alpha} / #alpha";
	
	TH1D* hRatioPt = FindObj<TH1D>(finput, tag, Form("Anti%s%s_Ratio_Pt", species.Data(), species.Data()));
	
	TCanvas* c0 = new TCanvas(Form("c0.Ratio%s",species.Data()), Form("Anti%s/%s ratio", species.Data(), species.Data()));
	c0->cd();
	
	hRatioPt->SetTitle("");
	hRatioPt->SetAxisRange(xmin,xmax,"X");
	hRatioPt->SetAxisRange(ymin,ymax,"Y");
	hRatioPt->SetLineColor(kRed);
	hRatioPt->SetMarkerStyle(kFullCircle);
	hRatioPt->SetMarkerColor(kRed);
	hRatioPt->SetYTitle(ytitle);
	hRatioPt->DrawCopy();
}
