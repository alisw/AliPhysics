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
// Author: Andrei Gheata, 04/01/2010

#include "AliTrigScheduledEntry.h"

#include "AliTrigScheduler.h"
#include "AliTrigDevice.h"

ClassImp(AliTrigScheduledEntry)
//==============================================================================
//
//   AliTrigScheduledEntry - ABC for scheduled responses of a device that is
//                           able to fire-up single response functions or the 
//                           full device scheduled sequence. The start time is
//                           in arbitrary units and in case it is 0 will not be
//                           considered when ordering by time by schedulers.
// 
//==============================================================================

//______________________________________________________________________________
AliTrigScheduledEntry::AliTrigScheduledEntry(const char *name, AliTrigDevice *device, Int_t start)
                      :TNamed(name,""), 
                       fStartTime(start), 
                       fDevice(device) 
{
// Default constructor. The only way to set the device.
}


ClassImp(AliTrigScheduledResponse)
//==============================================================================
//
//   AliTrigScheduledResponse - Scheduled device response function. Fires-up a
//                              single response function at a time.
//
//==============================================================================

//______________________________________________________________________________
AliTrigScheduledResponse::AliTrigScheduledResponse(const char *name, AliTrigDevice *device, Int_t output, Int_t start)
                         :AliTrigScheduledEntry(name, device, start),
                          fOutputID(output)
{
// Default constructor. The only way to set the device and response function.
}

//______________________________________________________________________________
void AliTrigScheduledResponse::FireUp(Int_t /*time*/)
{
// Virtual fire-up method. Calls the device response function.
   fDevice->Response(fOutputID);
}

ClassImp(AliTrigScheduledDevice)
//==============================================================================
//
//   AliTrigScheduledDevice - Scheduled entry for a full device sequence. Invokes
//                            the device scheduler when firing-up.
//
//==============================================================================

//______________________________________________________________________________
AliTrigScheduledDevice::AliTrigScheduledDevice(const char *name, AliTrigDevice *device, Int_t start)
                         :AliTrigScheduledEntry(name, device, start)
{
// Default constructor. The only way to set the device. Device scheduler must be set up.
}

//______________________________________________________________________________
void AliTrigScheduledDevice::FireUp(Int_t time)
{
// Virtual fire-up method. Calls the device response function.
   fDevice->GetScheduler()->FireUp(time);
}
