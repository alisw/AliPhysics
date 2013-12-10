#ifndef ALILNBA_H
#define ALILNBA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// coalescence parameter
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>

class TGraphErrors;

class AliLnBA: public TObject
{
  public:
	
	AliLnBA(const TString& protonSpectra, const TString& nucleusSpectra, const TString& protonTag, const TString& nucleusTag, const TString& outputFilename, const TString& otag, Int_t a, Int_t z);
	
	virtual ~AliLnBA();
	
	TGraphErrors* GetBAPt(const TGraphErrors* grPrtInvDYieldPt, const TGraphErrors* grNucInvDYieldPt, const TString& name) const;
	
	TString GetNucleusName() const { return fNucleusName; }
	Int_t GetA() const { return fA; }
	Int_t GetZ() const { return fZ; }
	
	Double_t Rside2Rlong(Double_t pt, Double_t B2, Double_t Cd) const;
	TGraphErrors* Rside2Rlong(const TGraphErrors* grB2, const TString& name, Double_t Cd) const;
	
	Int_t Run();
	
	void SetNucleus(Int_t a, Int_t z);
	void SetCd(Double_t Cd) { fCd = Cd; }
	
  private:
 
	AliLnBA(const AliLnBA& other);
	AliLnBA& operator=(const AliLnBA& other);
	
	Double_t GetErrorY(const TGraphErrors* gr, Double_t x0) const;
	
  private:
 
	TString fProtonSpectra; // proton spectra filename
	TString fProtonTag; // tag for proton file
	TString fNucleusSpectra; // nucleus spectra filename
	TString fNucleusTag; // tag for nucleus file
	
	TString fOutputFilename; // output filename
	TString fOutputTag; // tag for output file
	
	Int_t fA; // nucleus mass
	Int_t fZ; // nucleus charge
	TString fNucleusName; // nucleus name
	
	Double_t fCd; // correction factor for homogeneity volume
	
	ClassDef(AliLnBA,2)
};

#endif // ALILNBA_H
