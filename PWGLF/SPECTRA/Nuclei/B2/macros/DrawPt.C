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

#include <Riostream.h>
#include <TROOT.h>
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
#include <RooAbsPdf.h>
#include <RooAbsData.h>

#include "B2.h"

void DrawPt(const TString& inputFile="debug.root", const TString& tag="test", const TString& particle="AntiDeuteron", Int_t hiptbin=17, Bool_t m2pid=0, Int_t lowm2bin=9, Int_t him2bin=17)
{
//
// Draw corrected pt for debugging
//
	Double_t xmin = 0;
	Double_t xmax = 3.5;
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	// m2 data fitted models
	
	if(m2pid && (hiptbin>lowm2bin))
	{
		TH1D* hPidPt = (TH1D*)FindObj(finput, tag, Form("%s_PID_Pt",particle.Data()));
		Double_t binwidth = hPidPt->GetBinWidth(0);
		
		using namespace RooFit;
		
		// disable verbose in RooFit
		RooMsgService::instance().setSilentMode(1);
		RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
		
		TCanvas* c0 = new TCanvas(Form("%s.M2",particle.Data()), Form("M2 for %ss",particle.Data()));
		c0->Divide(3,3);
		
		TGraphErrors* grFitSigmaPt = new TGraphErrors();
		grFitSigmaPt->SetName(Form("%s_Fit_Sigma_Pt", particle.Data()));
		
		TGraphErrors* grFitM2Pt = new TGraphErrors();
		grFitM2Pt->SetName(Form("%s_Fit_M2_Pt", particle.Data()));
		
		for(Int_t i=lowm2bin, j=0; i<him2bin && i-lowm2bin < 9 && i < hiptbin; ++i)
		{
			c0->cd(i-lowm2bin+1);
			//gPad->SetLogy(0);
			
			RooWorkspace* w= (RooWorkspace*)FindObj(finput, tag, Form("%s_M2_%02d",particle.Data(),i));
			
			RooPlot* m2frame = w->var("m2")->frame();
			
			w->data("data")->plotOn(m2frame);
			
			w->pdf("model")->plotOn(m2frame, LineWidth(2));
			w->pdf("model")->plotOn(m2frame, Components(*(w->pdf("Sd"))),LineWidth(1), LineColor(9));
			w->pdf("model")->plotOn(m2frame, Components(*(w->pdf("Bkg"))),LineWidth(1), LineColor(46),LineStyle(kDashed));
			
			Double_t ptMin = (i-1)*binwidth;
			Double_t ptMax = i*binwidth;
			Double_t pt = (ptMin+ptMax)/2.;
			
			m2frame->SetTitle(Form("%0.2f < p_{T} < %0.2f GeV/c", ptMin, ptMax));
			m2frame->SetMinimum(0.2);
		
			m2frame->Draw();
			
			grFitSigmaPt->SetPoint(j, pt, w->var("width")->getVal());
			grFitSigmaPt->SetPointError(j, 0, w->var("width")->getError());
			
			grFitM2Pt->SetPoint(j, pt, w->var("mean")->getVal());
			grFitM2Pt->SetPointError(j++, 0, w->var("mean")->getError());
		}
		
		// model parameters
		
		TCanvas* c1 = new TCanvas(Form("%s.M2.FitParameters",particle.Data()), Form("M2 model parameters for %ss",particle.Data()));
		
		c1->Divide(2,1);
		
		c1->cd(1);
		gPad->SetGridx();
		gPad->SetGridy();
		
		grFitSigmaPt->SetMarkerStyle(kFullCircle);
		grFitSigmaPt->SetMarkerColor(kBlue);
		grFitSigmaPt->SetLineColor(kBlue);
		grFitSigmaPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		grFitSigmaPt->GetYaxis()->SetTitle("#sigma_{m^{2}}");
		grFitSigmaPt->Draw("ALP");
		
		c1->cd(2);
		gPad->SetGridx();
		gPad->SetGridy();
		
		grFitM2Pt->SetMarkerStyle(kFullCircle);
		grFitM2Pt->SetMarkerColor(kBlue);
		grFitM2Pt->SetLineColor(kBlue);
		grFitM2Pt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		grFitM2Pt->GetYaxis()->SetTitle("m^{2} (GeV^{2}/c^{4})");
		grFitM2Pt->Draw("ALP");
		
		//delete grFitSigmaPt;
		//delete grFitM2Pt;
	}
	
	// remaining corrections
	
	const Int_t kNum = 5;
	const TString kCorr[kNum]  = { "PID", "M2Corr", "SecCor", "Unfolded", "EffCor"};
	const TString kLabel[kNum] = { "Raw", "PID contamination", "Secondaries", "Unfolding","Efficiency" };
	const Int_t kColor[kNum]   = { kRed, kAzure, kOrange+1, kGreen-2, kGreen-3};
	const Int_t kMarker[kNum]  = { kFullCircle, kOpenCircle, kFullTriangleUp, kOpenTriangleUp, kFullCircle};
	
	TLegend* legend = new TLegend(0.5689655,0.6355932,0.8362069,0.8326271,0,"brNDC");
	legend->SetTextSize(0.03);
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	
	TH1D* hPt[kNum];
	
	for(Int_t i=0; i<kNum; ++i)
	{
		hPt[i] = (TH1D*)FindObj(finput, tag, particle + "_" + kCorr[i] + "_Pt");
		hPt[i]->SetLineColor(kColor[i]);
		hPt[i]->SetMarkerColor(kColor[i]);
		hPt[i]->SetMarkerStyle(kMarker[i]);
		legend->AddEntry(hPt[i], kLabel[i], "lp");
	}
	
	TCanvas* c2 = new TCanvas(Form("%s.Pt",particle.Data()), Form("Pt for %s",particle.Data()));
	c2->SetLogy();
	
	hPt[kNum-1]->SetTitle(particle.Data());
	hPt[kNum-1]->SetAxisRange(xmin,xmax,"X");
	hPt[kNum-1]->Draw("E");
	for(Int_t i=0; i<kNum-1; ++i) hPt[i]->Draw("sameE");
	
	legend->Draw();
}
