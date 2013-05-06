#ifndef ALIHISTTOOLSDIHADRONPID_H
#define ALIHISTTOOLSDIHADRONPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. * 
* See cxx source for full Copyright notice */ 
/* $Id$ */

class TF2;
class TH1F;
class TH2F;
class TH3F;

class AliHistToolsDiHadronPID {

public:
	AliHistToolsDiHadronPID() {};

protected:
	~AliHistToolsDiHadronPID() {};

public:

	// Histogram Manipulation.
	static TH1F* RebinVariableBinning(const TH1F* histIn, const Double_t* binsx, const Int_t Nbinsx, const Bool_t density = kTRUE);
	static TH1F* RebinVariableBinning(const TH1F* histIn, const TH1F* histAxis, const Bool_t density = kTRUE);
	static TH1F* TrimHisto(const TH1F* histo, const Int_t firstbin, const Int_t lastbin);
	static void ConstMinusHist(TH1F* histo, const Float_t cc = 1);
	static TH3F* MakeHist3D(const char* name, const char* title, 
		const Int_t nbinsX, const Double_t minX, const Double_t maxX,
		const Int_t nbinsY, const Double_t minY, const Double_t maxY,
		const Int_t nbinsZ, const Double_t* zaxis);

	// Creating histograms from other histograms or functions.
	static TH2F* Function2DToHist2D(const TF2* function, const TH2* grid);

	// Histogram Visualization.
	static TObjArray* CreateSpectraComparison(const char* name, const char* title, const TH1F* h1, const TH1F* h2, const Int_t markerstyle = 8, const Bool_t logy = kTRUE);

private:
	static Double_t* CreateAxis(const Int_t nbins, const Double_t min, const Double_t max);

};

#endif
