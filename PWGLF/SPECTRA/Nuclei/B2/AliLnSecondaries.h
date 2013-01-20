#ifndef ALILNSECONDARIES_H
#define ALILNSECONDARIES_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// removal of secondaries using DCA distributions
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <TString.h>

class TString;
class TH1D;
class TH2D;
class RooWorkspace;
class RooDataSet;
class RooRealVar;
class TF1;
class TFractionFitter;

class AliLnSecondaries: public TObject
{
  public:
	
	AliLnSecondaries(const TString& particle, const TString& dataFilename, const TString& simuFilename, const TString& outputFilename, const TString& otag);
	virtual ~AliLnSecondaries();
	
	const TString* GetOutputFilename() const { return &fOutputFilename; }
	
	Int_t Exec();
	
	void SetParticle(const TString& particle) { fParticle = particle; }
	
	void SetOutputTag(const TString& tag) { fOutputTag = tag; }
	
	void SetCorBins(Int_t lowbin, Int_t hibin) { fLowBin = lowbin; fHiBin = hibin; }
	void SetDCAxyInterval(Double_t lowdca, Double_t hidca) { fMinDCAxy = lowdca; fMaxDCAxy = hidca; }
	
	void SetNBin(Int_t nbin) { fNbin = nbin; }
	
	void SetProcedure(Int_t prod) { fFracProc=prod; }
	void SetAntiNucleusAsTemplate(Bool_t flag=1) { fANucTemplate = flag; }
	void SetMatDCAxyModel(Int_t model) { fMatDCAxyMod = model; }
	void SetScalingFactors(Double_t mat, Double_t fd) { fScMat=mat; fScFd=fd; }
	
	enum { kTFractionFitter=0, kRooFit, kMonteCarlo };
	enum { kGeantDCAxy=0, kFlatDCAxy};
	
  private:
 
	AliLnSecondaries(const AliLnSecondaries& other);
	AliLnSecondaries& operator=(const AliLnSecondaries& other);
	
	void GetFraction(TH1D* hPrimPt) const;
	void GetFraction(TH1D* hPrimFracPt, TH1D* hSecFracPt, const TH2D* hDCAxyPt, const TH2D* hMCDCAxyPt, const TH2D* hPrimDCAxyPt, const TH2D* hSecDCAxyPt, const TString& secName) const;
	void GetFraction(TH1D* hFracPt[3], const TH2D* hDCAxyPt, const TH2D* hMCDCAxyPt, const TH2D* hPrimDCAxyPt, const TH2D* hMatDCAxyPt, const TH2D* hFdwnDCAxyPt) const;
	
	Int_t GetTFFfractions(Double_t* frac, Double_t* err, TH1D* hData, TH1D* hPrim, TH1D* hMat, TH1D* hFdwn, Int_t ibin) const;
	Int_t GetTFFfractions(Double_t* frac, Double_t* err, TH1D* hData, TH1D* hPrim, TH1D* hSec, Int_t ibin, const TString& secName) const;
	
	void GetRooFitFractions(Double_t* frac, Double_t* err, const TH1D* hData, const TH1D* hPrim, const TH1D* hSec, Int_t ibin, const TString& secName) const;
	void GetRooFitFractions(Double_t* frac, Double_t* err, const TH1D* hData, const TH1D* hPrim, const TH1D* hMat, const TH1D* hFdwn, Int_t ibin) const;
	
	TH2D* GetFlatDCAxyPt(Double_t norm, const TH2D* hDCAxyPt, const TString& name) const;
	
	TF1* GetMatFraction(const TString& name) const;
	TF1* GetFdwnFraction(const TString& name) const;
	
	TH1D* ZeroClone(const TH1D* h, const TString& name) const;
	void WriteTFFdebug(TH1D* hData, TFractionFitter* fit, Int_t status, Int_t ibin, const char* contrib[], Double_t* frac, Int_t kmax) const;
	
  private:
  
	TString fParticle; // particle
	
	TString fDataFilename; // data filename
	TString fSimuFilename; // simulation filename
	TString fOutputFilename; // output filename
	TString fOutputTag; // tag for the ouput
	
	Int_t fLowBin; // low pt bin for the corrections
	Int_t fHiBin ; // high pt bin for the corrections
	
	Int_t fNbin; // for rebinning DCA distributions
	
	Double_t fMinDCAxy; // low DCAxy value
	Double_t fMaxDCAxy; // high DCAxy value
	
	Int_t fFracProc;  // procedure for estimating the fractions
	Int_t fMatDCAxyMod; // DCAxy model for secondaries from materials
	Bool_t fANucTemplate; // enable antinucleus as template for primaries
	
	Double_t fScMat; // scaling factor for material fraction
	Double_t fScFd; // scaling factor for feed-down fraction
	
	ClassDef(AliLnSecondaries,1)
};

#endif // ALILNSECONDARIES_H
