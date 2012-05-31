#ifndef ALITRIGMODULE_H
#define ALITRIGMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Andrei Gheata, 27/07/2009

//==============================================================================
//   AliTrigModule - Base class for trigger devices handling generic events. 
//      A module has arbitrary number of inputs and outputs. Derived classes must 
//      implement CreateDevice() and Trigger() metods.
//==============================================================================

#ifndef ALITRIGDEVICE_H
#include "AliTrigDevice.h"
#endif

class TObjArray;
class AliTrigEvent;

class AliTrigModule : public AliTrigDevice {

public:
  AliTrigModule() : AliTrigDevice(), fInputs(0), fOutputs(0), fOutputConnectors(0) {}
  AliTrigModule(const char *name, Int_t ninputs, Int_t noutputs) : AliTrigDevice(name, ninputs, noutputs), fInputs(0), fOutputs(0), fOutputConnectors(0) {}
  virtual ~AliTrigModule();

  virtual Bool_t            Connect(Int_t output, AliTrigDevice *other, Int_t at_input);
  virtual Bool_t            CreateDevice() = 0;
  void                      DefineInput(Int_t islot, AliTrigEvent *event);
  void                      DefineOutput(Int_t islot, AliTrigEvent *event);
  virtual Bool_t            Response(Int_t output);
  // Get/Set inputs
  AliTrigEvent             *GetInputValue(Int_t input) const;
  AliTrigEvent             *GetOutputValue(Int_t output) const;
  virtual void              ResetInputs();
  virtual Bool_t            SetInputType(Int_t /*input*/, const char */*classname*/) {return kFALSE;}
  virtual Bool_t            SetInputValue(Int_t /*input*/, Bool_t /*value*/)         {return kFALSE;}
  virtual Bool_t            SetInputValue(Int_t input, AliTrigEvent *event);
private:
   // Circuit response function. 
  AliTrigModule(const AliTrigModule &other);
  AliTrigModule &operator=(const AliTrigModule &other);
  virtual Bool_t            Trigger(Int_t ioutput) = 0;
   
protected:
  TObjArray                *fInputs;           // Array of input events
  TObjArray                *fOutputs;          // Array of output events
  TObjArray                *fOutputConnectors; // Array of output connectors
   
  ClassDef(AliTrigModule,1)  // Base class for a trigger module handling events
};
#endif
