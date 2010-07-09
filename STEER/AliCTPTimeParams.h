#ifndef ALICTPTIMEPARAMS_H
#define ALICTPTIMEPARAMS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#include<TNamed.h>
#include<TObjArray.h>
#include<AliCTPInputTimeParams.h>

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
		  void AddInput( TString& inputName, UInt_t& inputLevel, UInt_t inputDelay, TString inputEdge, UInt_t deltamin, UInt_t deltamax );
		  void AddDelayL0L1L2(Int_t delayL1L0, UInt_t delayL2L0);
	  virtual void Print(const Option_t* opt="") const;
		//Setters

		//Getters
		Int_t   GetDelayL1L0()   const { return fDelayL1L0; }
	        UInt_t  GetDelayL2L0()   const { return fDelayL2L0; }                  const TObjArray* GetInputTimeParams() const { return &fCTPInputTimeParams; }
      AliCTPInputTimeParams* GetTimeParamsForInput(TString inputname);
      Int_t GetDeltasforClass(TString classname,Int_t& deltamin,Int_t& deltamax);

                enum {kNMaxInputs = 60}; //CTP can manage up to 60 trigger detector inputs
private:
			Int_t  fDelayL1L0;
			UInt_t fDelayL2L0;
			TObjArray fCTPInputTimeParams;

  ClassDef( AliCTPTimeParams, 3 ) 
};

#endif
