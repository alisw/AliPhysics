#ifndef ALILNSPECTRA_H
#define ALILNSPECTRA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// (invariant) differential yield
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>

class TGraphErrors;
class TString;
class TF1;

class AliLnSpectra: public TObject
{
  public:
	
	AliLnSpectra(const TString& particle, const TString& ptFilename, const TString& tag, const TString& outputFilename, const TString& otag);
	
	virtual ~AliLnSpectra();
	
	TGraphErrors* GetDiffYieldPt(const TH1D* hPt, Int_t nEvent, const TString& name) const;
	TGraphErrors* GetInvDiffYieldPt(const TH1D* hPt, Int_t nEvent, const TString& name) const;
	TGraphErrors* GetInvDiffYieldPt(const TGraphErrors* grDYieldPt, const TString& name) const;
	
	Int_t Exec();
	
	void SetRapidityInterval(Double_t ymin, Double_t ymax) { fYMin = ymin; fYMax = ymax; }
	void SetExtrapolateToINEL(Bool_t flag=1) { fINEL = flag; }
	
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
	
	ClassDef(AliLnSpectra,1)
};

#endif // ALILNSPECTRA_H
