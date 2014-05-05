#ifndef ALIHISTTOOLSDIHADRONPID_H
#define ALIHISTTOOLSDIHADRONPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
* See cxx source for full Copyright notice */ 
/* $Id$ */

/*
class TF2;
class TH1F;
class TH2;
class TH2F;
class TH3F;
class TCanvas;
*/

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF2.h"
#include "TCanvas.h"

class AliHistToolsDiHadronPID {

public:
	AliHistToolsDiHadronPID() {};

protected:
	~AliHistToolsDiHadronPID() {};

public:

	// Histogram Manipulation.
	static TH1F* RebinVariableBinning(const TH1F* histIn, const Double_t* binsx, Int_t Nbinsx, Bool_t density = kTRUE);
	static TH1F* RebinVariableBinning(const TH1F* histIn, const TH1F* histAxis, Bool_t density = kTRUE);
	static TH1F* RebinVariableBinning(const TH1F* histIn, const TAxis* xaxis, Bool_t density = kTRUE);
	static TH1F* TrimHisto(const TH1F* histo, Int_t firstbin, Int_t lastbin);
	static void ConstMinusHist(TH1F* histo, Float_t cc = 1);
	static TH3F* MakeHist3D(const char* name, const char* title, 
		Int_t nbinsX, Double_t minX, Double_t maxX,
		Int_t nbinsY, Double_t minY, Double_t maxY,
		Int_t nbinsZ, const Double_t* zaxis);

	// Creating histograms from other histograms or functions.
	static TH2F* Function2DToHist2D(const TF2* function, const TH2* grid);

	// Histogram Visualization.
	static TCanvas* CreateSpectraComparison(const char* name, const char* title, const TH1F* h1, const TH1F* h2, Int_t markerstyle = 8, Bool_t logy = kTRUE);

//private:
	static Double_t* CreateAxis(Int_t nbins, Double_t min, Double_t max);

};

#endif
