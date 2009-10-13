#ifndef ALICTPTIMEPARAMS_H
#define ALICTPTIMEPARAMS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#include<TNamed.h>
#include<TObjArray.h>

class TNamed;

class AliCTPTimeParams : public TNamed {

public:
			AliCTPTimeParams();
              virtual   ~AliCTPTimeParams();
	
                        AliCTPTimeParams(const AliCTPTimeParams &timeparams);
      AliCTPTimeParams& operator=(const AliCTPTimeParams& timeparams);		

             
	static AliCTPTimeParams* LoadCTPTimeParams(TString filename);	
	static AliCTPTimeParams* LoadCTPTimeParamsFromString(const char* timeparams);
		Bool_t ProcessCTPTimeParamsLine(const char* line);
		  void AddInput( TString& inputName, UInt_t& inputLevel, UInt_t inputDelay, TString inputEdge );
		  void AddDelayL0L1L2(UInt_t delayL1L0, UInt_t delayL2L0);
	  virtual void Print(const Option_t* opt="") const;
		//Setters

		//Getters
		UInt_t  GetDelayL1L0()   const { return fDelayL1L0; }
	        UInt_t  GetDelayL2L0()   const { return fDelayL2L0; }                  const TObjArray* GetInputTimeParams() const { return &fCTPInputTimeParams; }
                enum {kNMaxInputs = 60}; //CTP can manage up to 60 trigger detector inputs
private:
			UInt_t fDelayL1L0;
			UInt_t fDelayL2L0;
			TObjArray fCTPInputTimeParams;

  ClassDef( AliCTPTimeParams, 1 ) 
};

#endif
