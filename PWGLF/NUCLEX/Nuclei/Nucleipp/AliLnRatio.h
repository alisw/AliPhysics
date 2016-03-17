#ifndef ALILNRATIO_H
#define ALILNRATIO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// antiparticle / particle ratio
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>

class AliLnRatio: public TObject
{
  public:
	
	AliLnRatio(const TString& species, const TString& ptFilename, const TString& ptTag, const TString& outputFilename, const TString& otag);
	
	virtual ~AliLnRatio();
	
	Int_t Exec();
	
  private:
 
	AliLnRatio(const AliLnRatio& other);
	AliLnRatio& operator=(const AliLnRatio& other);
	
  private:
 
	TString fSpecies; // particle species
	TString fPtFilename; // pt filename
	TString fPtTag; // tag for pt
	
	TString fOutputFilename; // output filename
	TString fOutputTag; // tag for the output
	
	ClassDef(AliLnRatio,1)
};

#endif // ALILNRATIO_H
