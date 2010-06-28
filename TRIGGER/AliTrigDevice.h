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

class AliTrigEvent;
class AliTrigScheduler;
class AliTrigScheduledResponse;

class AliTrigDevice : public TNamed {

private:
  AliTrigDevice(const AliTrigDevice &other);
  AliTrigDevice &operator=(const AliTrigDevice &other);

public:
  AliTrigDevice();
  AliTrigDevice(const char *name, Int_t ninputs, Int_t noutputs);
  virtual ~AliTrigDevice();

  virtual void              AddDevice(AliTrigDevice *other);
  Int_t                     GetNcomponents() const;
  AliTrigDevice            *GetComponent(Int_t n);
  AliTrigScheduledResponse *GetResponseFunction(const char *name);
  AliTrigScheduler         *GetScheduler() const {return fScheduler;}
  
  //____________________________________________________________________________
  // Device creation method to be implemented by derived classes. The response
  // functions are registered here. Connections between component devices should
  // also be handled in this method.
  virtual Bool_t            CreateDevice() {return kTRUE;}
  //____________________________________________________________________________
  // Connectivity to other devices. The method will create a connector between
  // an output of this device to one input of the other.
  virtual Bool_t            Connect(Int_t /*output*/, AliTrigDevice */*other*/, Int_t /*at_input*/) {return kTRUE;}

  //____________________________________________________________________________
  // Response functions to be implemented by specific devices. Has to propagate
  // the response to all connected devices. Representing the output #n of the device.
  virtual Bool_t            Response(Int_t output = 0) = 0;

  //____________________________________________________________________________
  // Create the response functions of the device.
  // The delay argument is in arbitrary time units with respect to the startup
  // reference. Note that the created scheduled entry must be registered to the
  // device scheduler via: fDevice->AddScheduledEntry() method
  AliTrigScheduledResponse *RegisterResponseFunction(const char *name, Int_t output, Int_t delay);

  //____________________________________________________________________________
  // Setting the value for a given input for digital devices of general ones
  // that are handling generic signals.
  virtual const char       *GetOutputType(Int_t /*output*/) {return 0;}
  virtual Bool_t            SetInputType(Int_t input, const char *classname) = 0;
  virtual Bool_t            SetInputValue(Int_t input, Bool_t value) = 0;
  virtual Bool_t            SetInputValue(Int_t input, AliTrigEvent *event) = 0;

  //____________________________________________________________________________
  // Device-dependent inputs reset method
  virtual void              ResetInputs() = 0;

  void                      SetNinputs(Int_t ninputs)   {fNinputs = ninputs;}
  void                      SetNoutputs(Int_t noutputs) {fNoutputs = noutputs;}
  Int_t                     GetNinputs() const {return fNinputs;}
  Int_t                     GetNoutputs() const {return fNoutputs;}
   
protected:
  Int_t                     fNinputs;            // Number of inputs
  Int_t                     fNoutputs;           // Number of outputs
  AliTrigScheduler         *fScheduler;         // Device scheduler
  TObjArray                *fComponents;        // Component devices
  TObjArray                *fResponseFunctions; // List of response functions
   
  ClassDef(AliTrigDevice,1)  // Base class for trigger devices
};
#endif
