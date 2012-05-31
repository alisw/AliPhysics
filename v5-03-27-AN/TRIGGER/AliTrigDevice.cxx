/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */
// Author: Andrei Gheata, 28/07/2009


#include "AliTrigDevice.h"

#include <TObjArray.h>
#include "AliTrigScheduler.h"
#include "AliTrigScheduledEntry.h"

ClassImp(AliTrigDevice)
//==============================================================================
//   AliTrigDevice - Generic device class. A device has a number of inputs and
// outputs. The data handled by the device can be either Boolean (digital
// devices) or arbitrary (wrapped by the class AliTrigSignal). A device must
// provide a response function that may depend on the output id. To replay the 
// device response function for a given output id, the device MUST register the
// output via the RegisterResponseFunction() method providing the delay in arbitrary time
// units. After the execution of the response for some output, the result will
// be propagated to all devices connected to this output. The method CreateDevice() 
// must be implemented by all devices and should connect all component devices
// and register all response functions. 
// The ResetInputs() method is called during simulation after the execution of 
// all response functions.
//==============================================================================

//______________________________________________________________________________
AliTrigDevice::AliTrigDevice()
              :TNamed(),
               fNinputs(0),
               fNoutputs(0),
               fScheduler(NULL),
               fComponents(NULL),
               fResponseFunctions(NULL) 
{
// I/O constructor.
}

//______________________________________________________________________________
AliTrigDevice::AliTrigDevice(const char *name, Int_t ninputs, Int_t noutputs)
              :TNamed(name, ""),
               fNinputs(ninputs),
               fNoutputs(noutputs),
               fScheduler(new AliTrigScheduler(name)),
               fComponents(NULL),
               fResponseFunctions(NULL) 
{
// Constructor.
}

//______________________________________________________________________________
AliTrigDevice::~AliTrigDevice()
{
// Destructor.
  delete fScheduler;
  if (fComponents) {fComponents->Delete(); delete fComponents;}
  if (fResponseFunctions) {fResponseFunctions->Delete(); delete fResponseFunctions;}
}   

//______________________________________________________________________________
void AliTrigDevice::AddDevice(AliTrigDevice *other)
{
// Add another device as component of this device.
  if (!fComponents) fComponents = new TObjArray();
  fComponents->Add(other);
}

//______________________________________________________________________________
Int_t AliTrigDevice::GetNcomponents() const
{
// Returns number of components.
  if (!fComponents) return 0;
  return fComponents->GetEntriesFast();
}

//______________________________________________________________________________
AliTrigDevice *AliTrigDevice::GetComponent(Int_t n)
{
// Get component at index n.
  if (!fComponents) return NULL;
  return (AliTrigDevice*)fComponents->At(n);
}

//______________________________________________________________________________
AliTrigScheduledResponse *AliTrigDevice::GetResponseFunction(const char *name)
{
// Get a response function by name.
  if (!fResponseFunctions) return NULL;
  return (AliTrigScheduledResponse*)fResponseFunctions->FindObject(name);
}  

//______________________________________________________________________________
AliTrigScheduledResponse *AliTrigDevice::RegisterResponseFunction(const char *name, Int_t output, Int_t delay)
{
// Creates a response function of the device. The delay argument is in arbitrary 
// time units with respect to the startup reference. Note that the created 
// scheduled entry must be registered to the device scheduler via: 
//    fDevice->AddScheduledEntry() method, otherwise it will not be replayed.
  if (!fResponseFunctions) fResponseFunctions = new TObjArray();
  if (fResponseFunctions->FindObject(name)) {
    Error("RegisterResponseFunction", "A response function named %s was already registered for device %s",
          name, GetName());
    return NULL;
  }        
  AliTrigScheduledResponse *response = new AliTrigScheduledResponse(name, (AliTrigDevice*)this, output, delay);
  fResponseFunctions->Add(response);
  return response;
}
