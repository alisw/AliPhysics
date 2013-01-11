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

// Draw histograms/graphs with same name located in different directories
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TLine.h>
#include <TStyle.h>

TGraphErrors* Divide(const TGraphErrors* grX, const TGraphErrors* grY, const TString& name);

void DrawDir(const TString& inputFile,
             const TString& graphName,
             const TString& refdir="",
             Double_t xmin = 0,
             Double_t xmax = 3.5,
             Double_t ymin = 0,
             Double_t ymax = 1.,
             const TString& xtitle = "",
             const TString& ytitle = "",
             Int_t option = 1,
             const TString& canvasName = "c0",
             const TString& canvasTitle = "DrawDir")
{
//
// draw histograms/graphs with same name located in different directories
// if option = 0 no comparison,
// if option = 1 draw a bottom pad with the comparison, and
// if option = 2 draw the comparasion in a different canvas
//
	const Int_t kMax = 10; // maximum number of subdirectories
	
	const Int_t kColor[kMax]  = { kRed, kBlue, kOrange+1, kGreen-2, kGreen+2, kAzure, kViolet+10, kAzure+2, kOrange+2, kSpring-7 };
	
	const Int_t kMarker[kMax] = { kFullCircle, kFullCircle, kFullCircle, kFullCircle, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar };
	
	TFile* finput = new TFile(inputFile.Data());
	if (finput->IsZombie()) exit(1);
	
	// iterate over the list of keys
	
	TObject* obj[kMax];
	TString* subdir[kMax];
	obj[0] = 0; // reference object
	Int_t nDir = 0;
	if(option>0) nDir = 1;
	
	TIter next(finput->GetListOfKeys());
	TKey* key = 0;
	while( (key = (TKey*)next()) && (nDir < kMax) )
	{
		TString classname = key->GetClassName();
		if(classname == "TDirectoryFile")
		{
			TObject* i = 0;
			finput->GetObject(Form("%s/%s;1", key->GetName(), graphName.Data()),i);
			if(i == 0)
			{
				finput->Error("GetObject","%s/%s not found", key->GetName(), graphName.Data());
				exit(1);
			}
			else if(i->InheritsFrom("TH1") || i->InheritsFrom("TGraph"))
			{
				Int_t j = 0;
				
				if(option==0) j = nDir++;
				else if(TString(key->GetName()) != refdir) j = nDir++;
				
				obj[j] = i;
				(dynamic_cast<TAttLine*>(i))->SetLineColor(kColor[j]);
				(dynamic_cast<TAttMarker*>(i))->SetMarkerColor(kColor[j]);
				(dynamic_cast<TAttMarker*>(i))->SetMarkerStyle(kMarker[j]);
				subdir[j] = new TString(key->GetName());
			}
			else
			{
				finput->Warning("GetObject", "%s does not contain %s", key->GetName(), graphName.Data());
			}
		}
	}
	
	// compare w.r.t. reference data
	
	TGraphErrors* grDiv[kMax];
	
	if(option > 0)
	{
		if(obj[0] == 0)
		{
			finput->Error("GetObject", "reference directory %s not found", refdir.Data());
			for(Int_t i=0; i < nDir; ++i) delete subdir[i];
			exit(1);
		}
		
		TGraphErrors* grY = 0;
		
		if(obj[0]->InheritsFrom("TH1"))
		{
			grY = new TGraphErrors(dynamic_cast<TH1*>(obj[0]));
		}
		else if(obj[0]->InheritsFrom("TGraph"))
		{
			grY = new TGraphErrors(*dynamic_cast<TGraphErrors*>(obj[0]));
		}
		
		TGraphErrors* grX = 0;
		
		for(Int_t i=1; i < nDir; ++i)
		{
			if(obj[i]->InheritsFrom("TH1"))
			{
				grX = new TGraphErrors(dynamic_cast<TH1*>(obj[i]));
			}
			else if(obj[i]->InheritsFrom("TGraph"))
			{
				grX = new TGraphErrors(*dynamic_cast<TGraphErrors*>(obj[i]));
			}
			
			grDiv[i] = Divide(grX, grY, Form("%s_Ratio", obj[i]->GetName()));
			
			grDiv[i]->SetLineColor(kColor[i]);
			grDiv[i]->SetMarkerColor(kColor[i]);
			grDiv[i]->SetMarkerStyle(kMarker[i]);
			
			delete grX;
		}
		
		delete grY;
	}
	
	// draw
	
	TStyle* st = new TStyle();
	
	st->SetOptTitle(0);
	st->SetOptStat(0);
	
	st->SetPadTickX(1);
	st->SetPadTickY(1);
	st->SetPadGridX(1);
	st->SetPadGridY(1);
	
	st->SetCanvasColor(0);
	st->SetFrameBorderMode(0);
	st->SetFrameFillColor(0);
	st->SetTitleFillColor(0);
	st->SetLabelFont(62,"XYZ");
	st->SetTitleFont(62,"XYZ");
	
	st->cd();
	gROOT->ForceStyle();
	
	TCanvas* c0 = new TCanvas(canvasName.Data(), canvasTitle.Data());
	
	if(option == 1) // create a top pad
	{
		TPad* top = new TPad("top", "Variable", 0, 0.25, 1, 1, 0, 0, 0);
		
		top->SetBottomMargin(0.);
		top->Draw();
		
		top->cd();
		
		TH1F* frm = top->DrawFrame(xmin, ymin , xmax, ymax);
		frm->GetYaxis()->SetTitle(ytitle);
	}
	else // only draw the frame
	{
		TH1F* frm = c0->DrawFrame(xmin, ymin , xmax, ymax);
		frm->GetXaxis()->SetTitle(xtitle);
		frm->GetYaxis()->SetTitle(ytitle);
	}
	
	// draw objects in current pad
	for(Int_t i=0; i < nDir; ++i)
	{
		if(obj[i]->InheritsFrom("TH1"))
		{
			obj[i]->Draw("same");
		}
		else if(obj[i]->InheritsFrom("TGraph"))
		{
			obj[i]->Draw("PZ");
		}
	}
	
	// build a legend
	TLegend* legend = new TLegend(0.5718391,0.6991525,0.8390805,0.8368644,0,"brNDC");
	legend->SetTextSize(0.03);
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	
	for(Int_t i=0; i < nDir; ++i)
	{
		legend->AddEntry(obj[i], subdir[i]->Data(), "lp");
	}
	
	legend->Draw();
	
	if(option == 1) // create a bottom pad
	{
		c0->cd();
		
		TPad* bottom = new TPad("bottom", "ratio", 0, 0, 1, 0.25, 0, 0, 0);
		
		bottom->SetBottomMargin(0.3);
		bottom->SetTopMargin(0.);
		bottom->Draw();
		
		bottom->cd();
		
		TH1F* frm = bottom->DrawFrame(xmin, 0.7 , xmax, 1.3);
		
		frm->GetXaxis()->SetLabelSize(0.13);
		frm->GetXaxis()->SetTitle(xtitle);
		frm->GetXaxis()->SetTitleSize(0.13);
		frm->GetYaxis()->SetTitle("Ratio");
		frm->GetYaxis()->SetLabelSize(0.1);
		frm->GetYaxis()->SetTitleSize(0.12);
		frm->GetYaxis()->SetTitleOffset(0.3);
	}
	else if(option == 2) // create a new canvas
	{
		TCanvas* c1 = new TCanvas(Form("%s.Ratio",canvasName.Data()), Form("%s ratio",canvasTitle.Data()));
		
		TH1F* frm = c1->DrawFrame(xmin, 0.5 ,xmax, 1.5);
		frm->GetXaxis()->SetTitle(xtitle);
		frm->GetYaxis()->SetTitle("Ratio");
	}
	
	// draw comparison
	if(nDir>0 && option > 0)
	{
		for(Int_t i=1; i<nDir; ++i)
		{
			grDiv[i]->Draw("PZ");
		}
		
		// draw a red line for the reference
		TLine* ref = new TLine(xmin,1,xmax,1);
		ref->SetLineWidth(1);
		ref->SetLineColor(kColor[0]);
		ref->Draw();
		
		if(option == 2)
		{
			TLegend* legendRatio = new TLegend(0.5718391,0.6991525,0.8390805,0.8368644,0,"brNDC");
			legendRatio->SetTextSize(0.03);
			legendRatio->SetFillColor(0);
			legendRatio->SetBorderSize(0);
			
			legendRatio->AddEntry(ref, subdir[0]->Data(), "l");
			for(Int_t i=1; i < nDir; ++i)
			{
				legendRatio->AddEntry(grDiv[i], subdir[i]->Data(), "lp");
			}
			
			legendRatio->Draw();
		}
	}
	
	// release memory
}

TGraphErrors* Divide(const TGraphErrors* grX, const TGraphErrors* grY, const TString& name)
{
//
// Divide two TGraphErrors using first one as reference
//
	TGraphErrors* grQ = new TGraphErrors();
	grQ->SetName(name.Data());
	
	for(Int_t i=0; i < grX->GetN(); ++i)
	{
		Double_t x, y1;
		grX->GetPoint(i, x, y1);
		Double_t y2 = grY->Eval(x);
		
		if(y1 == 0 || y2 == 0) continue;
		
		Double_t r = y1/y2;
		
		Double_t erry1 = grX->GetErrorY(i);
		Double_t erry2 = grY->GetErrorY(i); // guess error
		Double_t err = r*TMath::Sqrt(TMath::Power(erry1/y1,2)+TMath::Power(erry2/y2,2));
		
		grQ->SetPoint(i, x, r);
		grQ->SetPointError(i, grX->GetErrorX(i), err);
	}
	
	return grQ;
}
