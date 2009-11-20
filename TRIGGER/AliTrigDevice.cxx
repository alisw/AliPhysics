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

//==============================================================================
//   AliTrigDevice - Generic device class. A device has a number of inputs and
// outputs. The data handled by the device can be either Boolean (digital
// devices) or arbitrary (wrapped by the class AliTrigSignal). A device must
// provide a response function that may depend on the output id. To replay the 
// device response function for a given output id, the device MUST register the
// output via the RegisterOutput() method providing the delay in arbitrary time
// units. After the execution of the response for some output, the Emit() method
// will be invoked to propagate the computed result to all devices connected to 
// this output. The method Connect() must be implemented by all devices and should
// create a connector specific to the device types that are linked. 
// The ResetInputs() method is called during simulation after the execution of 
// all response functions.
//==============================================================================

#include "AliTrigDevice.h"

ClassImp(AliTrigDevice)

//______________________________________________________________________________
AliTrigDevice::AliTrigDevice()
{
// Destructor.
   if (fComponents) delete fComponents;
}   

//______________________________________________________________________________
AliTrigDevice::AliTrigDevice(const AliTrigDevice &other)
                    :TNamed(other), 
                     fNinputs(other.fNinputs), 
                     fNoutputs(other.fNoutputs),
                     fComponents(0)
{
// Copy ctor.
  if (other.fComponents) {
     fComponents = new TObjArray();
     TIter next(other.fComponents);
     AliTrigDevice *dev;
     while ((dev=(AliTrigDevice*)next())) fComponents->Add(dev);
  }   
}                        

//______________________________________________________________________________
AliTrigDevice& AliTrigDevice::operator=(const AliTrigDevice &other)
{
// Assignment
  if (&other == this) return *this;
  TNamed::operator=(other);
  fNinputs  = other.fNinputs;
  fNoutputs = other.fNoutputs;
  fComponents = 0;
  if (other.fComponents) {
     fComponents = new TObjArray();
     TIter next(other.fComponents);
     AliTrigDevice *dev;
     while ((dev=(AliTrigDevice*)next())) fComponents->Add(dev);
  }   
  return *this;
}

//______________________________________________________________________________
void AliTrigDevice::AddDevice(AliTrigDevice *other)
{
// Add another device as component of this device.
  if (!fComponents) fComponents = new TObjArray();
  fComponents->Add(other);
}

//______________________________________________________________________________
Int_t AliTrigDevice::GetNcomponents()
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
   return fComponents->At(n);
}
   
//______________________________________________________________________________
Bool_t AliTrigDevice::RegisterResponseFunction(AliTrigScheduler *calendar, UInt_t output, Int_t delay) const
{
// Register the response functions to be replayed by the provided scheduler.
// The delay argument is in arbitrary time units with respect to the startup
// reference of the simulation.
// CALLING SEQUENCE:
//   The delay behaves like a BUSY gate for the device and MUST be called when
//   configuring the whole setup of devices. The method may fail in case the
//   determinism is not satisfied with other connected devices providing inputs.
   return calendar->RegisterResponseFunction(this, output, delay);
}
