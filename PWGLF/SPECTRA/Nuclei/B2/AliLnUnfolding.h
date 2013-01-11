#ifndef ALILNUNFOLDING_H
#define ALILNUNFOLDING_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// unfolding correction for the pt
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <TString.h>

class TString;

class AliLnUnfolding: public TObject
{
  public:
	
	AliLnUnfolding(const TString& particle, const TString& simuFilename, const TString& outputFilename, const TString& otag);
	virtual ~AliLnUnfolding();
	
	const TString* GetOutputFilename() const { return &fOutputFilename; }
	
	Int_t Exec();
	
	void SetParticle(const TString& particle) { fParticle = particle; }
	
	void SetOutputTag(const TString& tag) { fOutputTag = tag; }
	
  private:
 
	AliLnUnfolding(const AliLnUnfolding& other);
	AliLnUnfolding& operator=(const AliLnUnfolding& other);
	
  private:

	TString fParticle; // particle
	
	TString fSimuFilename; // simulation filename
	TString fOutputFilename; // output filename
	TString fOutputTag; // tag for the ouput
	
	ClassDef(AliLnUnfolding,1)
};

#endif // ALILNUNFOLDING_H
