#ifndef ALILNFAKETRACKS_H
#define ALILNFAKETRACKS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// fake track fraction
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <TString.h>

class TString;

class AliLnFakeTracks: public TObject
{
  public:
	
	AliLnFakeTracks(const TString& particle, const TString& simuFilename, const TString& outputFilename, const TString& otag);
	virtual ~AliLnFakeTracks();
	
	const TString* GetOutputFilename() const { return &fOutputFilename; }
	
	Int_t Exec();
	
	void SetParticle(const TString& particle) { fParticle = particle; }
	
	void SetOutputTag(const TString& tag) { fOutputTag = tag; }
	
  private:
 
	AliLnFakeTracks(const AliLnFakeTracks& other);
	AliLnFakeTracks& operator=(const AliLnFakeTracks& other);
	
  private:
 
	TString fParticle; // particle
	
	TString fSimuFilename; // simulation filename
	TString fOutputFilename; // output filename
	TString fOutputTag; // tag for the ouput
	
	
	ClassDef(AliLnFakeTracks,1)
};

#endif // ALILNFAKETRACKS_H
