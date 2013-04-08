#ifndef ALIHISTTOOLSDIHADRONPID_H
#define ALIHISTTOOLSDIHADRONPID_H

class AliHistToolsDiHadronPID {

public:
	AliHistToolsDiHadronPID() {};

protected:
	~AliHistToolsDiHadronPID() {};

public:

	// Histogram Manipulation.
	static TH1F* RebinVariableBinning(TH1F* histIn, Double_t* binsx, Int_t Nbinsx, Bool_t density = kTRUE);
	static TH1F* RebinVariableBinning(TH1F* histIn, TH1F* histAxis, Bool_t density = kTRUE) {

		// Rebins histogram histIn to the x-axis of histAxis
		TAxis* xaxis = histAxis->GetXaxis();
		Int_t nbinsx = xaxis->GetNbins();
		const Double_t* binsx = (xaxis->GetXbins())->GetArray();
		return RebinVariableBinning(histIn, const_cast<Double_t*>(binsx), nbinsx, density);

	}
	static TH1F* TrimHisto(TH1F* histo, Int_t firstbin, Int_t lastbin);
	static void ConstMinusHist(TH1F* histo, const Float_t cc = 1) {

		// h -> (c-h)
		Int_t nbins = histo->GetNbinsX();
		for (Int_t iBin = 0; iBin < (nbins + 1); iBin++) {
			Float_t bincontent = histo->GetBinContent(iBin);
			histo->SetBinContent(iBin,(cc - bincontent));
		}

	}
	static TH3F* MakeHist3D(const char* name, const char* title, 
		Int_t nbinsX, Double_t minX, Double_t maxX,
		Int_t nbinsY, Double_t minY, Double_t maxY,
		Int_t nbinsZ, const Float_t* zaxis) {

		const Float_t* xaxis = const_cast<Float_t*>(CreateAxis(nbinsX,minX,maxX));
		const Float_t* yaxis = const_cast<Float_t*>(CreateAxis(nbinsY,minY,maxY));

		TH3F* hout = new TH3F(name,title,nbinsX,xaxis,nbinsY,yaxis,nbinsZ,zaxis);
	
		return hout;
	}

	// Histogram Visualization.
	static TObjArray* CreateSpectraComparison(const char* name, const char* title, TH1F* h1, TH1F* h2, Int_t markerstyle = 8, Bool_t logy = kTRUE);

private:
	static Float_t* CreateAxis(const Int_t nbins, Double_t min, Double_t max) {
		if (nbins <= 0) return 0x0;
		if (max < min) return 0x0;

		Float_t* axis = new Float_t[nbins + 1];
		Float_t binsize = (max - min)/((Float_t)nbins);
		for (Int_t iBin = 0; iBin < nbins + 1; iBin++) {
			axis[iBin] = min + ((Float_t)iBin) * binsize;
		}

		return axis;
	}


};

#endif
