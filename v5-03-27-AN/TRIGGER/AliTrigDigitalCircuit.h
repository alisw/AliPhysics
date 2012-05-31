#ifndef ALITRIGDIGITALCIRCUIT_H
#define ALITRIGDIGITALCIRCUIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 27/07/2009

//==============================================================================
//   AliTrigDigitalCircuit - Base class for digital circuits having N Boolean 
//      inputs and one Boolean output. Derived classes must implement the pure 
//      virtual method Trigger() that will return the Boolean response 
//      of the circuit as function of the inputs.
//==============================================================================

#ifndef ALITRIGDEVICE_H
#include "AliTrigDevice.h"
#endif

#ifndef ROOT_TBits
#include "TBits.h"
#endif

class AliTrigEvent;
class AliTrigConnector;

class AliTrigDigitalCircuit : public AliTrigDevice {

public:
  AliTrigDigitalCircuit() : AliTrigDevice(), fLastOutput(kFALSE), fConnector(0), fInputs() {}
  AliTrigDigitalCircuit(const char *name, UInt_t ninputs) : AliTrigDevice(name, ninputs, 1), fLastOutput(kFALSE), fConnector(0), fInputs(ninputs) {}
  virtual ~AliTrigDigitalCircuit();

  virtual Bool_t            Connect(Int_t output, AliTrigDevice *other, Int_t at_input);
  virtual Bool_t            Response(Int_t output);
  // Get/Set inputs
  Bool_t                    GetInputValue(Int_t input) const {return fInputs.TestBitNumber(input);}
  virtual void              ResetInputs() {fInputs.ResetAllBits();}
  virtual Bool_t            SetInputType(Int_t /*input*/, const char */*classname*/) {return kFALSE;}
  virtual Bool_t            SetInputValue(Int_t input, Bool_t value) {fInputs.SetBitNumber(input,value); return kTRUE;}
  virtual Bool_t            SetInputValue(Int_t /*input*/, AliTrigEvent */*signal*/) {return  kFALSE;}
private:
   // Circuit response function. 
  AliTrigDigitalCircuit(const AliTrigDigitalCircuit &other);
  AliTrigDigitalCircuit &operator=(const AliTrigDigitalCircuit &other);
  virtual Bool_t            Trigger() = 0;
   
protected:
  Bool_t                    fLastOutput; // Output recorded after the last Response() call.
  AliTrigConnector         *fConnector;  // Connector for the circuit output
  TBits                     fInputs;
   
  ClassDef(AliTrigDigitalCircuit,1)  // Base class for digital circuits
};
#endif
