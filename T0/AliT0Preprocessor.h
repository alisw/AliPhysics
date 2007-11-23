#ifndef ALI_T0_PREPROCESSOR_H
#define ALI_T0_PREPRECESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliPreprocessor.h"

class AliT0DataDCS;

class AliT0Preprocessor: public AliPreprocessor 
{
  public:
        AliT0Preprocessor(): AliPreprocessor("T00",0) { }
        AliT0Preprocessor(AliShuttleInterface* shuttle);
	virtual ~AliT0Preprocessor();
  
  protected:
        virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
	virtual UInt_t Process(TMap* dcsAliasMap);

  private:
	//AliT0Calc *fData;
	AliT0DataDCS *fData;	
	ClassDef(AliT0Preprocessor, 1)
};

typedef AliT0Preprocessor AliSTARTPreprocessor; // for backward compatibility

#endif
