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

// Draw corrected pt
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <RooWorkspace.h>
#include <RooMsgService.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooAbsData.h>
#include <RooAbsPdf.h>
#endif

#include "B2.h"

void DrawPt(const TString& inputFile="debug.root", const TString& tag="test", const TString& particle="AntiDeuteron", Double_t ptmax=3., Bool_t m2pid=0, Double_t ptpid=1.2)
{
//
// Draw corrected pt for debugging
//
	Double_t xmin = 0;
	Double_t xmax = 3.5;
	
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(1);
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	TH1D* hPidPt = FindObj<TH1D>(finput, tag, Form("%s_PID_Pt",particle.Data()));
	
	Int_t hiptbin  = hPidPt->GetXaxis()->FindFixBin(ptmax);
	Int_t lowm2bin = hPidPt->GetXaxis()->FindFixBin(ptpid);
	Int_t him2bin  = hPidPt->GetXaxis()->FindFixBin(ptmax);
	
	// m2 data fitted models
	
	if(m2pid && (hiptbin>lowm2bin))
	{
		using namespace RooFit;
		
		// disable verbose in RooFit
		RooMsgService::instance().setSilentMode(1);
		RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
		
		TCanvas* c0 = new TCanvas(Form("%s.M2",particle.Data()), Form("M2 for %ss",particle.Data()));
		c0->Divide(4,5);
		
		TGraphErrors* grFitSigmaPt = new TGraphErrors();
		grFitSigmaPt->SetName(Form("%s_Fit_Sigma_Pt", particle.Data()));
		
		TGraphErrors* grFitM2Pt = new TGraphErrors();
		grFitM2Pt->SetName(Form("%s_Fit_M2_Pt", particle.Data()));
		
		for(Int_t i=lowm2bin, j=0; i<him2bin && i-lowm2bin < 20 && i < hiptbin; ++i)
		{
			c0->cd(i-lowm2bin+1);
			
			RooWorkspace* w= FindObj<RooWorkspace>(finput, tag, Form("%s_M2_%02d",particle.Data(),i));
			
			RooPlot* m2frame = w->var("x")->frame();
			
			w->data("data")->plotOn(m2frame);
			
			w->pdf("model")->plotOn(m2frame, Components(*(w->pdf("Sd"))),LineWidth(1), LineColor(8));
			w->pdf("model")->plotOn(m2frame, Components(*(w->pdf("Bkg"))),LineWidth(1), LineColor(46),LineStyle(kDashed));
			w->pdf("model")->plotOn(m2frame, LineWidth(1));
			
			m2frame->SetTitle(Form("%0.2f < #it{p}_{T} < %0.2f GeV/#it{c}", hPidPt->GetBinLowEdge(i), hPidPt->GetBinLowEdge(i)+hPidPt->GetBinWidth(i)));
			m2frame->SetMinimum(0.2);
			
			m2frame->Draw();
			
			Double_t pt = hPidPt->GetBinCenter(i);
			
			grFitSigmaPt->SetPoint(j, pt, w->var("sigma")->getVal());
			grFitSigmaPt->SetPointError(j, 0, w->var("sigma")->getError());
			
			grFitM2Pt->SetPoint(j, pt, w->var("mu")->getVal());
			grFitM2Pt->SetPointError(j++, 0, w->var("mu")->getError());
		}
		
		c0->Update();
		
		// model parameters
		
		TCanvas* c1 = new TCanvas(Form("%s.M2.FitParameters",particle.Data()), Form("M2 model parameters for %ss",particle.Data()));
		
		c1->Divide(2,1);
		
		c1->cd(1);
		
		grFitSigmaPt->SetMarkerStyle(kFullCircle);
		grFitSigmaPt->SetMarkerColor(kBlue);
		grFitSigmaPt->SetLineColor(kBlue);
		grFitSigmaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		grFitSigmaPt->GetYaxis()->SetTitle("#sigma");
		grFitSigmaPt->Draw("ALP");
		
		c1->cd(2);
		
		grFitM2Pt->SetMarkerStyle(kFullCircle);
		grFitM2Pt->SetMarkerColor(kBlue);
		grFitM2Pt->SetLineColor(kBlue);
		grFitM2Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		grFitM2Pt->GetYaxis()->SetTitle("#mu");
		grFitM2Pt->Draw("ALP");
		
		//delete grFitSigmaPt;
		//delete grFitM2Pt;
		
		c1->Update();
	}
	
	// remaining corrections
	
	TCanvas* c2 = new TCanvas(Form("%s.Pt",particle.Data()), Form("Pt for %s",particle.Data()));
	c2->SetLogy();
	
	const Int_t kNum = 4;
	const TString kCorr[kNum]  = { "PID", "PidCorr", "SecCorr", "EffCorr"};
	const TString kLabel[kNum] = { "Raw", "PID", "Secondaries","Efficiency" };
	const Int_t kColor[]   = { kRed, kAzure, kOrange+1, kGreen-3, kGreen-2};
	const Int_t kMarker[]  = { kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCircle, kOpenTriangleUp};
	
	TLegend* legend = new TLegend(0.5689655,0.6355932,0.8362069,0.8326271,0,"brNDC");
	legend->SetTextSize(0.03);
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	
	TH1D* hPt[kNum];
	
	for(Int_t i=0; i<kNum; ++i)
	{
		hPt[i] = FindObj<TH1D>(finput, tag, particle + "_" + kCorr[i] + "_Pt");
		hPt[i]->SetLineColor(kColor[i]);
		hPt[i]->SetMarkerColor(kColor[i]);
		hPt[i]->SetMarkerStyle(kMarker[i]);
		legend->AddEntry(hPt[i], kLabel[i], "lp");
	}
	
	hPt[kNum-1]->SetTitle(particle.Data());
	hPt[kNum-1]->SetAxisRange(xmin,xmax,"X");
	hPt[kNum-1]->Draw("E");
	
	for(Int_t i=0; i<kNum-1; ++i) hPt[i]->Draw("sameE");
	legend->Draw();
	
	c2->Update();
}
