#ifndef ALILNSPECTRA_H
#define ALILNSPECTRA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// invariant differential yields and cross sections
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>

class TGraphErrors;
class TF1;

class AliLnSpectra: public TObject
{
  public:
	
	AliLnSpectra(const TString& particle, const TString& ptFilename, const TString& tag, const TString& outputFilename, const TString& otag, const Double_t xsec[3]);
	
	virtual ~AliLnSpectra();
	
	TGraphErrors* GetDiffYieldPt(const TH1D* hPt, Int_t nEvent, const TString& name) const;
	TGraphErrors* GetInvDiffYieldPt(const TH1D* hPt, Int_t nEvent, const TString& name) const;
	TGraphErrors* GetInvDiffYieldPt(const TGraphErrors* grDYieldPt, const TString& name) const;
	TGraphErrors* GetInvDiffXsectionPt(const TGraphErrors* grInvDYieldPt, const Double_t* sigma, const TString& name) const;
	
	TGraphErrors* AddSystError(const TGraphErrors* gr, Double_t percent, const TString& name) const;
	
	Int_t Exec();
	
	void SetRapidityInterval(Double_t ymin, Double_t ymax) { fYMin = ymin; fYMax = ymax; }
	void SetExtrapolateToINEL(Bool_t flag=1) { fINEL = flag; }
	
	void SetOnlyGeneration(Bool_t flag=1) { fIsOnlyGen = flag; }
	
	void SetScalingFactor(Double_t syserr=1) { fSysErr = syserr; }
	
	void SetInelXSection(Double_t xsec, Double_t statErr, Double_t systErr) { fInelXsec[0] = xsec; fInelXsec[1] = statErr; fInelXsec[2] = systErr; }
	
	TF1* Tsallis(Double_t m0, const TString& name, Double_t xmin=0., Double_t xmax=10.) const;
	TF1* Tsallis(Double_t m0, Double_t xsect, const TString& name, Double_t xmin=0., Double_t xmax=10.) const;
	TF1* TsallisDiffYield(Double_t m0, const TString& name, Double_t xmin=0., Double_t xmax=10.) const;
	
  private:
 
	AliLnSpectra(const AliLnSpectra& other);
	AliLnSpectra& operator=(const AliLnSpectra& other);
	
  private:
 
	TString fParticle; // particle name
	
	TString fPtFilename; // pt filename
	TString fTag; // tag for pt file
	
	TString fOutputFilename; // output filename
	TString fOutputTag; // tag for output file
	
	Double_t fYMin; // rapidity low limit
	Double_t fYMax; // rapidity high limit
	
	Bool_t fINEL; // extrapolate to inelastic events
	Bool_t fIsOnlyGen; // if no need for correction
	
	Double_t fSysErr; // variation for systematics
	Double_t fInelXsec[3]; // total inelastic cross section, syst. and stat. errors (mb)
	
	ClassDef(AliLnSpectra,1)
};

#endif // ALILNSPECTRA_H
