#ifndef ALILNB2_H
#define ALILNB2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// coalescence parameter
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>

class TGraphErrors;

class AliLnB2: public TObject
{
  public:
	
	AliLnB2(const TString& protonSpectra, const TString& nucleusSpectra, const TString& protonTag, const TString& nucleusTag, const TString& outputFilename, const TString& otag, Int_t a, Int_t z);
	
	virtual ~AliLnB2();
	
	TGraphErrors* GetBAPt(const TGraphErrors* grPrtInvDYieldPt, const TGraphErrors* grNucInvDYieldPt, const TString& name, Double_t errPt=0) const;
	
	Double_t Rside2Rlong(Double_t pt, Double_t B2, Double_t Cd) const;
	TGraphErrors* Rside2Rlong(const TGraphErrors* grB2, const TString& name, Double_t Cd, Double_t errPt=0) const;
	
	Int_t Run();
	
	void SetNucleus(Int_t a, Int_t z);
	void SetCd(Double_t Cd) { fCd = Cd; }
	
  private:
 
	AliLnB2(const AliLnB2& other);
	AliLnB2& operator=(const AliLnB2& other);
	
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
	
	ClassDef(AliLnB2,1)
};

#endif // ALILNB2_H
