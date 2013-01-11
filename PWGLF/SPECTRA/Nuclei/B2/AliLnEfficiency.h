#ifndef ALILNEFFICIENCY_H
#define ALILNEFFICIENCY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// efficiency correction
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <TString.h>

class TString;
class TF1;

class AliLnEfficiency: public TObject
{
  public:
	
	AliLnEfficiency(const TString& particle, const TString& simuFilename, const TString& outputFilename, const TString& otag);
	virtual ~AliLnEfficiency();
	
	TF1* GetProtonG3FlukaCor(const TString& name) const;
	TF1* GetAntiProtonG3FlukaCor(const TString& name) const;
	
	const TString* GetOutputFilename() const { return &fOutputFilename; }
	
	Int_t Exec();
	
	void SetParticle(const TString& particle) { fParticle = particle; }
	
	void SetOutputTag(const TString& tag) { fOutputTag = tag; }
	
	void SetG3Fluka(Bool_t flag=1) { fG3Fluka = flag; }
	
  private:
 
	AliLnEfficiency(const AliLnEfficiency& other);
	AliLnEfficiency& operator=(const AliLnEfficiency& other);
	
  private:
 
	TString fParticle; // particle name
	
	TString fSimuFilename; // simulation filename
	TString fOutputFilename; // output filename
	TString fOutputTag; // tag for the ouput file
	
	Bool_t fG3Fluka; // use G3/FLUKA correction for (anti)protons
	
	ClassDef(AliLnEfficiency,1)
};

#endif // ALILNEFFICIENCY_H
