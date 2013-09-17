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

// Draw corrections of secondaries
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#endif

#include "B2.h"

Int_t GetNumberCols(Int_t lowbin, Int_t hibin);

TCanvas* DrawMC(TH1D** hData, TH1D** hPrim, TH1D** hMat, TH1D** hFdwn, Int_t lowbin, Int_t hibin, const TString& name, const TString& title, const TH1D* hPt, Double_t dcaxyMin, Double_t dcaxyMax);

TCanvas* DrawFit(TH1D** hData, TH1D** hDataFit, TH1D** hPrim, TH1D** hMat, TH1D** hFdwn, Int_t lowbin, Int_t hibin, const TString& name, const TString& title, const TH1D* hPt, Double_t dcaxyMin, Double_t dcaxyMax);

TCanvas* DrawGoF(TH1D** hData, TH1D** hDataFit, Int_t lowbin, Int_t hibin, const TString& name, const TString& title, const TH1D* hPt, Double_t dcaxyMin, Double_t dcaxyMax);

void DrawSec(const TString& inputFile="debug.root", const TString& tag="test", const TString& particle="Proton", Double_t ptmin=0.5, Double_t ptmax=3., Double_t dcaxyMin=-1.5, Double_t dcaxyMax=1.5)
{
//
// draw corrections of secondaries
//
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(1);
	
	const Int_t kNBin = 50;
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	// data and templates
	
	TH1D* hData[kNBin];
	TH1D* hDataMC[kNBin];
	TH1D* hPrim[kNBin];
	TH1D* hMat[kNBin];
	TH1D* hFdwn[kNBin];
	
	// fitted data
	
	TH1D* hDataFit[kNBin];
	TH1D* hPrimFit[kNBin];
	TH1D* hMatFit[kNBin];
	TH1D* hFdwnFit[kNBin];
	
	TH1D* hPidPt = FindObj<TH1D>(finput, tag, particle + "_PID_Pt");
	
	Int_t lowbin = hPidPt->GetXaxis()->FindFixBin(ptmin);
	Int_t hibin  = hPidPt->GetXaxis()->FindFixBin(ptmax);
	
	for(Int_t i=lowbin; i<hibin && i<kNBin; ++i)
	{
		// MC
		hData[i]   = FindObj<TH1D>(finput, tag, particle + Form("_Data_DCAxy_%02d",i));
		hDataMC[i] = FindObj<TH1D>(finput, tag, particle + Form("_SimData_DCAxy_%02d",i));
		hPrim[i]   = FindObj<TH1D>(finput, tag, particle + Form("_Prim_DCAxy_%02d",i));
		hMat[i]    = FindObj<TH1D>(finput, tag, particle + Form("_Mat_DCAxy_%02d",i));
		
		if(particle.Contains("Proton"))
		{
			hFdwn[i] = FindObj<TH1D>(finput, tag, particle + Form("_Fdwn_DCAxy_%02d",i));
		}
		else
		{
			hFdwn[i] = 0;
		}
		
		// fitted data
		hDataFit[i] = FindObj<TH1D>(finput, tag, particle + Form("_Fit_Data_DCAxy_%02d",i));
		hPrimFit[i] = FindObj<TH1D>(finput, tag, particle + Form("_Fit_Prim_DCAxy_%02d",i));
		hMatFit[i]  = FindObj<TH1D>(finput, tag, particle + Form("_Fit_Mat_DCAxy_%02d",i));
		
		if(particle.Contains("Proton"))
		{
			hFdwnFit[i] = FindObj<TH1D>(finput, tag, particle + Form("_Fit_Fdwn_DCAxy_%02d",i));
		}
		else
		{
			hFdwnFit[i] = 0;
		}
	}
	
	// draw
	
	TCanvas* c0 = 0;
	
	if(hDataMC[lowbin]!=0)
	{
		c0 = DrawMC(hDataMC, hPrim, hMat, hFdwn, lowbin, hibin, particle + ".MC", Form("MC simulation for %s",particle.Data()), hPidPt, dcaxyMin, dcaxyMax);
		c0->Update();
	}
	
	TCanvas* c1;
	
	if(hDataFit[lowbin]!=0)
	{
		c1 = DrawFit(hData, hDataFit, hPrimFit, hMatFit, hFdwnFit, lowbin, hibin, particle + ".Fit", Form("Fitted DCA for %s",particle.Data()), hPidPt, dcaxyMin, dcaxyMax);
		c1->Update();
	}
	
	TCanvas* c2;
	
	if(hDataFit[lowbin]!=0)
	{
		c2 = DrawGoF(hData, hDataFit, lowbin, hibin, particle + ".GoF", Form("GoF for %s",particle.Data()), hPidPt, dcaxyMin, dcaxyMax);
		c2->Update();
	}
}

Int_t GetNumberCols(Int_t lowbin, Int_t hibin)
{
//
// number of rows and cols for plotting the slices
// square matrix (cols = rows)
//
	return TMath::CeilNint(TMath::Sqrt(TMath::Abs(hibin - lowbin)));
}

TCanvas* DrawMC(TH1D** hData, TH1D** hPrim, TH1D** hMat, TH1D** hFdwn, Int_t lowbin, Int_t hibin, const TString& name, const TString& title, const TH1D* hPt, Double_t dcaxyMin, Double_t dcaxyMax)
{
//
// Draw MC DCA distributions for each pt bin
//
	Int_t ncol = GetNumberCols(lowbin, hibin);
	
	TCanvas* c0 = new TCanvas(name.Data(), title.Data());
	c0->Divide(ncol, ncol);
	
	for(Int_t i=lowbin; i < hibin; ++i)
	{
		c0->cd(i+1-lowbin);
		
		gPad->SetLogy();
		
		hData[i]->SetTitle(Form("%0.2f < p_{T} < %0.2f GeV/c", hPt->GetBinLowEdge(i), hPt->GetBinLowEdge(i)+hPt->GetBinWidth(i)));
		hData[i]->SetMinimum(0.2);
		
		hData[i]->SetAxisRange(dcaxyMin, dcaxyMax, "X");
		
		hData[i]->SetMarkerStyle(kFullCircle);
		hData[i]->Draw("E");
		
		if(hPrim[i] != 0)
		{
			hPrim[i]->SetLineColor(9);
			hPrim[i]->Draw("sameHist");
		}
		
		if(hMat[i] != 0)
		{
			hMat[i]->SetLineColor(46);
			hMat[i]->Draw("sameHist");
		}
		
		if(hFdwn != 0)
		if(hFdwn[i] != 0)
		{
			hFdwn[i]->SetLineColor(8);
			hFdwn[i]->Draw("sameHist");
		}
	}
	
	return c0;
}

TCanvas* DrawFit(TH1D** hData, TH1D** hDataFit, TH1D** hPrim, TH1D** hMat, TH1D** hFdwn, Int_t lowbin, Int_t hibin, const TString& name, const TString& title, const TH1D* hPt, Double_t dcaxyMin, Double_t dcaxyMax)
{
//
// Draw DCA fitted model for each pt bin in current canvas
//
	Int_t ncol = GetNumberCols(lowbin, hibin);
	
	TCanvas* c0 = new TCanvas(name.Data(), title.Data());
	c0->Divide(ncol, ncol);
	
	for(Int_t i=lowbin; i < hibin && i < ncol*ncol; ++i)
	{
		c0->cd(i+1-lowbin);
		
		gPad->SetLogy();
		
		hData[i]->SetTitle(Form("%0.2f < #it{p}_{T} < %0.2f GeV/c", hPt->GetBinLowEdge(i), hPt->GetBinLowEdge(i)+hPt->GetBinWidth(i)));
		hData[i]->SetMinimum(0.2);
		
		hData[i]->SetAxisRange(dcaxyMin, dcaxyMax, "X");
		
		hData[i]->SetMarkerStyle(kFullCircle);
		hData[i]->Draw("E");
		
		hDataFit[i]->SetLineColor(2);
		hDataFit[i]->Draw("sameHist");
		
		if(hPrim[i] != 0)
		{
			hPrim[i]->SetLineColor(9);
			hPrim[i]->Draw("sameHist");
		}
		
		if(hMat[i] != 0)
		{
			hMat[i]->SetLineColor(46);
			hMat[i]->Draw("sameHist");
		}
		
		if(hFdwn != 0)
		if(hFdwn[i] != 0)
		{
			hFdwn[i]->SetLineColor(8);
			hFdwn[i]->Draw("sameHist");
		}
	}
	
	return c0;
}

TCanvas* DrawGoF(TH1D** hData, TH1D** hDataFit, Int_t lowbin, Int_t hibin, const TString& name, const TString& title, const TH1D* hPt, Double_t dcaxyMin, Double_t dcaxyMax)
{
//
// Draw a goodness of fit for each pt bin consisting of data/fit
//
	Int_t ncol = GetNumberCols(lowbin, hibin);
	
	TCanvas* c0 = new TCanvas(name.Data(), title.Data());
	c0->Divide(ncol, ncol);
	
	for(Int_t i=lowbin; i < hibin; ++i)
	{
		c0->cd(i+1-lowbin);
		
		TH1D* hDiv = (TH1D*)hData[i]->Clone("_Div_Data_Fit_");
		hDiv->Sumw2();
		hDiv->Divide(hDataFit[i]);
		
		hDiv->SetTitle(Form("%0.2f < p_{T} < %0.2f GeV/c", hPt->GetBinLowEdge(i), hPt->GetBinLowEdge(i)+hPt->GetBinWidth(i)));
		hDiv->SetAxisRange(dcaxyMin, dcaxyMax, "X");
		hDiv->SetAxisRange(0,3,"Y");
		hDiv->SetMarkerStyle(kFullCircle);
		hDiv->DrawCopy("E");
		
		delete hDiv;
	}
	
	return c0;
}
