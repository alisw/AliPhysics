/************************************************************************* 
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. * 
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

// -----------------------------------------------------------------------
//  Tools for drawing/ manipulating histograms. (NEEDS CLEANUP!)
// -----------------------------------------------------------------------
//  Author: Misha Veldhoen (misha.veldhoen@cern.ch)

#include <iostream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TArray.h"

#include "AliHistToolsDiHadronPID.h"

using namespace std;

// -----------------------------------------------------------------------
TH1F* AliHistToolsDiHadronPID::RebinVariableBinning(
	const TH1F* histIn, const Double_t* binsx, const Int_t Nbinsx, const Bool_t density) {

	// Rebins a histogram (hin) with a variable binning to a histogram
	// with another variable binning (binsx). If the "density" flag is set,
	// then it is expected that the bins are per unit x-axis, otherwise an
	// absolute count is assumed.
	//
	// TODO: determine over/under-flow bins.

	// Gather info from the original histogram.
	const TArrayD* binsxIn = histIn->GetXaxis()->GetXbins();
	TArrayD* binsxOut = new TArrayD(Nbinsx + 1,binsx);
	Int_t NbinsxIn = binsxIn->GetSize() - 1;
	const char* nameIn = histIn->GetName();
	const char* titleIn = histIn->GetTitle();
	const char* xAxisTitleIn = histIn->GetXaxis()->GetTitle(); 
	const char* yAxisTitleIn = histIn->GetYaxis()->GetTitle();

	// Gather info from the output histogram and create it.
	Int_t NbinsxOut = binsxOut->GetSize() - 1;
	const char* nameOut = Form("%s (rebinned)",nameIn);
	const char* titleOut = Form("%s_rebinned;%s;%s",titleIn,xAxisTitleIn,yAxisTitleIn);
	TH1F* histOut = new TH1F(nameOut,titleOut,NbinsxOut,binsxOut->GetArray());
	//histOut->Sumw2();

	// Filling the bins.
	for (Int_t iBinOut = 0; iBinOut < NbinsxOut; iBinOut++) {

		Double_t binOutLowEdge = binsxOut->At(iBinOut);
		Double_t binOutHighEdge = binsxOut->At(iBinOut+1);
		Double_t binOutWidth = binOutHighEdge - binOutLowEdge;
		//Double_t binOutCenter = binOutHighEdge - binOutWidth/2.;

		// Setting all errors of the outgoing histogram to zero (needed later):
		histOut->SetBinError(iBinOut+1,0.);

		for (Int_t iBinIn = 0; iBinIn < NbinsxIn; iBinIn++) {
			
			Double_t binInLowEdge = binsxIn->At(iBinIn);
			Double_t binInHighEdge = binsxIn->At(iBinIn+1);
			Double_t binInWidth = binInHighEdge - binInLowEdge;
			//Double_t binInCenter = binInHighEdge - binInWidth/2.;
			Double_t binInContent = histIn->GetBinContent(iBinIn+1);
			Double_t binInError2 = histIn->GetBinError(iBinIn+1)*histIn->GetBinError(iBinIn+1);

			/* -- Determining what percentage of the in-bin is included 
			   in the current out-bin -- */
			Double_t PercentageIncluded = 0.;

			// Position of the low edge of the in-bin:
			// In  1   2   3            
			// Out   |   |
			Int_t binInLowEdgePos = 0;
			if (binInLowEdge < binOutLowEdge) binInLowEdgePos = 1;
			else if (binInLowEdge > binOutHighEdge) binInLowEdgePos = 3;
			else binInLowEdgePos = 2;

			// Same for the high edge of the in-bin:
			Int_t binInHighEdgePos = 0;
			if (binInHighEdge < binOutLowEdge) binInHighEdgePos = 1;
			else if (binInHighEdge > binOutHighEdge) binInHighEdgePos = 3;
			else binInHighEdgePos = 2;

			// In-bin not included.
			// In        | |     or     | |
			// Out |   |                    |    |
			if ( binInLowEdgePos == 3 || binInHighEdgePos == 1 ) {
				//cout<<"In-bin not included."<<endl;
				PercentageIncluded = 0.;
			}

			// In-bin partially included (left side).
			// In    |    |
			// Out |  |
			else if ( binInLowEdgePos == 2 && binInHighEdgePos == 3 ) {
				//cout<<"In-bin partially included (left side)."<<endl;
				PercentageIncluded = (binOutHighEdge - binInLowEdge)/(binInHighEdge - binInLowEdge);
			}

			// In-bin partially included (middle).
			// In    |     |
			// Out     |  |
			else if ( binInLowEdgePos == 1 && binInHighEdgePos == 3 ) {
				//cout<<"In-bin partially included (middle)."<<endl;
				PercentageIncluded = (binOutHighEdge - binOutLowEdge)/(binInHighEdge - binInLowEdge);
			}

			// In-bin partially included (right side).
			// In    |     |
			// Out        |  |
			else if ( binInLowEdgePos == 1 && binInHighEdgePos == 2 ) {
				//cout<<"In-bin partially included (right side)."<<endl;
				PercentageIncluded = (binInHighEdge - binOutLowEdge)/(binInHighEdge - binInLowEdge);
			}

			// In-bin completely included.
			// In    |  |
			// Out  |    |
			else if ( binInLowEdgePos == 2 && binInHighEdgePos == 2 ) {
				//cout<<"In-bin completely included."<<endl;
				PercentageIncluded = 1.;
			}

			// Give some status:
			//cout<<"Bin Out: "<<iBinOut+1<<" ("<<binOutLowEdge<<","<<binOutHighEdge<<")  Bin In: "<<iBinIn+1<<" ("<<binInLowEdge<<","<<binInHighEdge<<"), content: "<<binInContent<<" +/- "<<TMath::Sqrt(binInError2)<<". Bin in Pos: "<<binInLowEdgePos<<" "<<binInHighEdgePos<<" pct: "<<PercentageIncluded<<endl; 

			//else cout<<"Something weird just happened."<<endl;

			Double_t binOutAddedError2 = 0.;
			binOutAddedError2 = binInError2 * PercentageIncluded * PercentageIncluded;
			Double_t binOutAddedContent = 0.;
			binOutAddedContent = binInContent * PercentageIncluded;

			if (density) {
				binOutAddedContent*=binInWidth;
			}
			
			//histOut->Fill(binOutCenter,binOutAddedContent);
			histOut->SetBinContent(iBinOut+1,histOut->GetBinContent(iBinOut+1)+binOutAddedContent);

			//cout<<histOut->GetBinError(iBinOut+1)<<endl;
			//cout<<binOutAddedError2<<endl;
			histOut->SetBinError(iBinOut+1,histOut->GetBinError(iBinOut+1)+binOutAddedError2);
			//cout<<histOut->GetBinError(iBinOut+1)<<endl;
		}

		// Once a certain binOut has been filled all the way we have to scale
		// it if the "density" flag is set.

		histOut->SetBinContent(iBinOut + 1,histOut->GetBinContent(iBinOut+1)/binOutWidth);
		histOut->SetBinError(iBinOut + 1,TMath::Sqrt(histOut->GetBinError(iBinOut+1)));
		//cout<<histOut->GetBinError(iBinOut+1)<<endl;
	}

	return histOut;

}

// -----------------------------------------------------------------------
TH1F* AliHistToolsDiHadronPID::RebinVariableBinning(const TH1F* histIn, const TH1F* histAxis, const Bool_t density) {

	// Rebins histogram histIn to the x-axis of histAxis
	TAxis* xaxis = histAxis->GetXaxis();
	Int_t nbinsx = xaxis->GetNbins();
	const Double_t* binsx = (xaxis->GetXbins())->GetArray();
	return RebinVariableBinning(histIn, const_cast<Double_t*>(binsx), nbinsx, density);

}

// -----------------------------------------------------------------------
TH1F* AliHistToolsDiHadronPID::TrimHisto(const TH1F* histo, const Int_t firstbin, const Int_t lastbin) {

	const char* name = histo->GetName();
	const char* title = histo->GetTitle();
	if (firstbin < 0) return 0x0;
	if (lastbin > histo->GetNbinsX()) return 0x0;
	if (firstbin > lastbin) return 0x0;

	Int_t Nbins = lastbin - firstbin + 1;
	cout<<firstbin<<" "<<lastbin<<" "<<Nbins<<endl;
	//const Int_t NbinsC = Nbins;
	Double_t newaxis[41];

	for (Int_t ii = firstbin; ii < lastbin+1; ii++) {
		newaxis[ii - firstbin] = histo->GetBinLowEdge(ii); 
		//cout<<ii-firstbin<<" "<<newaxis[ii - firstbin]<<endl;
	}
	newaxis[Nbins] = histo->GetBinLowEdge(lastbin+1);

	TH1F* hout = new TH1F("hout","hout",Nbins,newaxis);
	hout->SetDirectory(0);
	hout->SetTitle(title);
	hout->SetName(name);

	for (Int_t ii = 0; ii < Nbins+1; ii++) {
		hout->SetBinContent(ii,histo->GetBinContent(firstbin+ii));
		hout->SetBinError(ii,histo->GetBinError(firstbin+ii));
	} 

	return hout;

}

// -----------------------------------------------------------------------
void AliHistToolsDiHadronPID::ConstMinusHist(TH1F* histo, const Float_t cc) {

	// h -> (c-h)
	Int_t nbins = histo->GetNbinsX();
	for (Int_t iBin = 0; iBin < (nbins + 1); iBin++) {
		Float_t bincontent = histo->GetBinContent(iBin);
		histo->SetBinContent(iBin,(cc - bincontent));
	}

}

// -----------------------------------------------------------------------
TH3F* AliHistToolsDiHadronPID::MakeHist3D(const char* name, const char* title, 
	const Int_t nbinsX, const Double_t minX, const Double_t maxX,
	const Int_t nbinsY, const Double_t minY, const Double_t maxY,
	const Int_t nbinsZ, const Double_t* zaxis) {

	const Double_t* xaxis = const_cast<Double_t*>(CreateAxis(nbinsX,minX,maxX));
	const Double_t* yaxis = const_cast<Double_t*>(CreateAxis(nbinsY,minY,maxY));

	TH3F* hout = new TH3F(name,title,nbinsX,xaxis,nbinsY,yaxis,nbinsZ,zaxis);

	return hout;

}

// -----------------------------------------------------------------------
TH2F* AliHistToolsDiHadronPID::Function2DToHist2D(const TF2* function, const TH2* grid) {

	// Creates a 2D histogram of "function" using the binning of 
	// the histogram "grid". 
	//  - We assume "grid" to have fixed binwidth.
	//  - Bins are only filled if they are within the domain of the function.
	//  - Histogram takes over the color of the function.

	// Gathering info about the axes.
	TAxis* Xaxis = grid->GetXaxis();
	Int_t NbinsX = Xaxis->GetNbins();
	TAxis* Yaxis = grid->GetYaxis();
	Int_t NbinsY = Yaxis->GetNbins();

	// Determining function range.
	Double_t Xmin = 0.; 
	Double_t Xmax = 0.; 
	Double_t Ymin = 0.; 
	Double_t Ymax = 0.;
	function->GetRange(Xmin, Ymin, Xmax, Ymax);

	// Creating the histogram.
	TH2F* hout = new TH2F(Form("%s_hist", function->GetName()),
		Form("%s (Hist);%s;%s", function->GetTitle(), Xaxis->GetTitle(), Yaxis->GetTitle()),
		NbinsX, Xaxis->GetBinLowEdge(1), Xaxis->GetBinUpEdge(NbinsX),
		NbinsY, Yaxis->GetBinLowEdge(1), Yaxis->GetBinUpEdge(NbinsY));
	hout->SetLineColor(function->GetLineColor());
	hout->SetMarkerColor(function->GetLineColor());
	hout->SetDirectory(0);

	// Filling the histogram.
	for (Int_t iBinX = 1; iBinX < (NbinsX + 1); iBinX++) {
		for (Int_t iBinY = 1; iBinY < (NbinsY + 1); iBinY++) {
			Double_t CoordX = Xaxis->GetBinCenter(iBinX);
			Double_t CoordY = Yaxis->GetBinCenter(iBinY);

			if ( (CoordX < Xmin) || (CoordX > Xmax) || (CoordY < Ymin) || (CoordY > Ymax) ) {continue;}

			hout->SetBinContent(iBinX, iBinY, function->Eval(CoordX, CoordY));
		}		
	}

	return hout;

}

// -----------------------------------------------------------------------
TObjArray* AliHistToolsDiHadronPID::CreateSpectraComparison(const char* name, 
	const char* title, const TH1F* h1, const TH1F* h2, const Int_t markerstyle, const Bool_t logy) {

	// - Creates a window comparing two histograms h1, and h2.
	// - Returns an array of pointers to the objects created
	//   in this function.

	//Int_t optstat = gStyle->GetOptStat();
	gStyle->SetOptStat(0);
	//c->UseCurrentStyle();
	//gStyle->SetOptStat(optstat);

	TH1F* spectra[2];
	TH1F* division;
	Int_t spectracolor[2] = {1,2};
	Int_t divisioncolor = 4;

	// Cloning and configuring spectra.
	spectra[0] = (TH1F*)h1->Clone();
	//spectra[0]->SetDirectory(0);
	spectra[1] = (TH1F*)h2->Clone();
	//spactra[1]->SetDirectory(0);	
	for (Int_t iSpectra = 0; iSpectra < 2; iSpectra++) {
		spectra[iSpectra]->Sumw2();
		spectra[iSpectra]->SetMarkerStyle(markerstyle);
		spectra[iSpectra]->SetLineColor(spectracolor[iSpectra]);
		spectra[iSpectra]->SetMarkerColor(spectracolor[iSpectra]);
	}

	// Creating the division.
	division = (TH1F*)spectra[0]->Clone("hDivision");
	division->Divide(spectra[1]);
	division->SetLineColor(divisioncolor);
	division->SetMarkerColor(divisioncolor);
	division->GetYaxis()->SetTitle("Ratio");

	// Creating window
	//TCanvas* c = new TCanvas(name,title,0,0,400,400);
	TCanvas* c = TCanvas::MakeDefCanvas();
	c->SetName(name);
	c->SetTitle(title);
	TPad* p1 = new TPad("p1","pad1",0,0.25,1,1);
	p1->SetRightMargin(0.05);
	p1->SetBottomMargin(0.);
	TPad* p2 = new TPad("p2","pad2",0,0,1,0.25);
	p2->SetTopMargin(0.);
	p2->SetRightMargin(0.05);
	p2->SetBottomMargin(0.3);
	p2->SetFillStyle(4000);
	p1->Draw();
	p2->Draw();

	// Determining plotting range.
	Double_t max = TMath::Max(spectra[0]->GetMaximum(),spectra[1]->GetMaximum());
	Double_t min = TMath::Min(spectra[0]->GetMinimum(),spectra[1]->GetMinimum());
	Double_t range = max-min;
	spectra[0]->SetMinimum(min-0.05*range);
	spectra[0]->SetMaximum(max+0.05*range);

	// Drawing
	p1->cd();
	if (logy) {p1->SetLogy();}
	spectra[0]->Draw("e");
	spectra[0]->GetXaxis()->SetLabelSize(0);
	spectra[0]->GetYaxis()->SetLabelSize(0.05);
	spectra[0]->GetYaxis()->SetTitleSize(0.05);
	spectra[0]->GetYaxis()->SetTitleOffset(1.1);
	//Double_t labelsize = spectra[0]->GetYaxis()->GetLabelSize();
	spectra[1]->Draw("same e");
	p2->cd();
	division->SetTitle("");
	division->GetXaxis()->SetLabelSize(0.1);
	division->GetXaxis()->SetTitleSize(0.15);
	division->GetXaxis()->SetTitleOffset(0.6);
	division->GetYaxis()->SetLabelSize(0.1);
	division->GetYaxis()->SetTitleSize(0.12);
	division->GetYaxis()->SetTitleOffset(0.3);
	division->Draw("e");

	//c->UseCurrentStyle();

	//gStyle->SetOptStat(optstat);

	TLegend* legend = new TLegend(0.75,0.8,0.95,0.95);
	legend->AddEntry(spectra[0],h1->GetTitle(),"lpe");
	legend->AddEntry(spectra[1],h2->GetTitle(),"lpe");
	p1->cd();
	legend->Draw("same");

	// returning the created objects.
	TObjArray* returnobjects = new TObjArray(6);
	returnobjects->AddLast(c);
	returnobjects->AddLast(p1);
	returnobjects->AddLast(p2);
	returnobjects->AddLast(spectra[0]);
	returnobjects->AddLast(spectra[1]);
	returnobjects->AddLast(division);

	return returnobjects;

}

// -----------------------------------------------------------------------
Double_t* AliHistToolsDiHadronPID::CreateAxis(const Int_t nbins, const Double_t min, const Double_t max) {

	if (nbins <= 0) return 0x0;
	if (max < min) return 0x0;

	Double_t* axis = new Double_t[nbins + 1];
	Double_t binsize = (max - min)/((Double_t)nbins);
	for (Int_t iBin = 0; iBin < nbins + 1; iBin++) {
		axis[iBin] = min + ((Double_t)iBin) * binsize;
	}

	return axis;
}






