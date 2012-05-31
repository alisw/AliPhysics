#ifndef ALI_T0_PREPROCESSOR_COSMIC_H
#define ALI_T0_PREPRECESSOR_COSMIC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliT0PreprocessorCosmic.h 24021 2008-02-19 18:53:25Z alla $ */


// T0 preprocessor. 
// Takes data from DCS and passes it to the class AliTOFDataDCS for processing and writes the result to the Reference DB.
// Takes data form DAQ (both from Laser Calibration and Physics runs), processes it, and stores either to OCDB or to Reference DB.

#include "AliPreprocessor.h"

class AliT0DataDCS;

class AliT0PreprocessorCosmic: public AliPreprocessor 
{
  public:
        AliT0PreprocessorCosmic(): AliPreprocessor("T00",0),  
	  fData(0)
 { }
        AliT0PreprocessorCosmic(AliShuttleInterface* shuttle);
	virtual ~AliT0PreprocessorCosmic();
  
  protected:
        virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
	virtual UInt_t Process(TMap* dcsAliasMap);
	virtual Bool_t ProcessDCS();

  private:
	AliT0PreprocessorCosmic(const AliT0PreprocessorCosmic & proc); // copy constructor	
	AliT0PreprocessorCosmic& operator=(const AliT0PreprocessorCosmic&); //operator
	UInt_t ProcessDCSDataPoints(TMap* dcsAliasMap);
 	UInt_t ProcessLaser();
 	UInt_t ProcessPhysics();
 	AliT0DataDCS *fData;			// Data member to process DCS data	
 
	ClassDef(AliT0PreprocessorCosmic, 1)
};

#endif
