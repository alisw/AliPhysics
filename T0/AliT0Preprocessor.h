#ifndef ALI_T0_PREPROCESSOR_H
#define ALI_T0_PREPRECESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */


// T0 preprocessor. 
// Takes data from DCS and passes it to the class AliT0DataDCS for processing and writes the result to the Reference DB.
// Takes data form DAQ (both from Laser Calibration and Physics runs), processes it, and stores either to OCDB or to Reference DB.

#include "AliPreprocessor.h"

class AliT0DataDCS;

class AliT0Preprocessor: public AliPreprocessor 
{
  public:
        AliT0Preprocessor(): AliPreprocessor("T00",0),  
	  fData(0)
 { }
        AliT0Preprocessor(AliShuttleInterface* shuttle);
	virtual ~AliT0Preprocessor();
  
  protected:
        virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
	virtual UInt_t Process(TMap* dcsAliasMap);
	virtual Bool_t ProcessDCS();

  private:
	AliT0Preprocessor(const AliT0Preprocessor & proc); // copy constructor	
	AliT0Preprocessor& operator=(const AliT0Preprocessor&); //operator
	UInt_t ProcessDCSDataPoints(TMap* dcsAliasMap);
	UInt_t ProcessLaser();
 	UInt_t ProcessPhysics();
 	AliT0DataDCS *fData;			// Data member to process DCS data	
 
	ClassDef(AliT0Preprocessor, 2)
};

typedef AliT0Preprocessor AliSTARTPreprocessor; // for backward compatibility

#endif
