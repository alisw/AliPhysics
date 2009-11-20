#ifndef ALITRIGDEVICE_H
#define ALITRIGDEVICE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 27/07/2009

//==============================================================================
//   AliTrigDevice - Base class for a generic device.
//==============================================================================

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class AliTrigScheduler;

class AliTrigDevice : public TNamed {

public:
  AliTrigDevice() : TNamed(), fNinputs(0), fNoutputs(0) {}
  AliTrigDevice(const char *name, UInt_t ninputs, UInt_t noutputs) 
              : TNamed(name,""), fNinputs(ninputs), fNoutputs(noutputs) {}
  AliTrigDevice(const AliTrigDevice &other);
  virtual ~AliTrigDevice();
  AliTrigDevice &operator=(const AliTrigDevice &other);

  virtual void              AddDevice(AliTrigDevice *other);
  Int_t                     GetNcomponents() const;
  AliTrigDevice            *GetComponent(Int_t n);
  //____________________________________________________________________________
  // Connectivity to other devices. The method will create a connector between
  // an output of this device to one input of the other.
  virtual Bool_t            Connect(UInt_t output, AliTrigDevice *other, UInt_t at_input) = 0;

  //____________________________________________________________________________
  // Response functions to be overloaded by specific devices. Has to propagate
  // the response to all connected devices. Representing the output #n of the device.
  virtual Bool_t            Response(UInt_t output = 0) = 0;

  //____________________________________________________________________________
  // Register the response functions to be replayed by the global scheduler.
  // The delay argument is in arbitrary time units with respect to the startup
  // reference.
  Bool_t                    RegisterResponseFunction(AliTrigScheduler *calendar, UInt_t output, Int_t delay, Int_t group) const;

  //____________________________________________________________________________
  // Setting the value for a given input for digital devices of general ones
  // that are handling generic signals.
  virtual const char       *GetOutputType(UInt_t output) {return 0;}
  virtual Bool_t            SetInputType(UInt_t input, const char *classname) = 0;
  virtual Bool_t            SetInputValue(UInt_t input, Bool_t value) = 0;
  virtual Bool_t            SetInputValue(UInt_t input, AliTrigSignal *signal) = 0;

  //____________________________________________________________________________
  // Device-dependent inputs reset method
  virtual void              ResetInputs() = 0;

  UInt_t                    GetNinputs() const {return fNinputs;}
  UInt_t                    GetNoutputs() const {return fNoutputs;}
   
protected:
  UInt_t                    fNinputs;  // Number of inputs.
  UInt_t                    fNoutputs; // Number of outputs.
  TObjArray                *fComponents; // Component devices if any
   
  ClassDef(AliTrigDevice,1)  // Base class for trigger devices
};
#endif
