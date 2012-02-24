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
//   AliTrigConnector - General connector class. A connector links a feeder
// device output to a number of other device inputs (clients). It transmits the 
// signal to all clients and calls their SetInput() method.
//==============================================================================

#include "AliTrigConnector.h"

#include <TObjArray.h>
#include <TClass.h>

#include "AliTrigDevice.h"
#include "AliTrigEvent.h"


ClassImp(AliTrigConnector)

//______________________________________________________________________________
AliTrigConnector::~AliTrigConnector()
{
// Destructor.
  if (fInputs) delete [] fInputs;
  if (fDevices) delete fDevices;
}

//______________________________________________________________________________
void AliTrigConnector::Connect(AliTrigDevice *client, Int_t input)
{
// Adds the device and its input to the list of clients.
  // Loop array of inputs to check if this input is already connected.
  for (Int_t i=0; i<fNclients; i++) {
    if (fInputs[i]==input && fDevices->At(i)==client) {
      Info("Connect", "Output #%d of device %s already connected to input #%d of device%s",
           fOutput, fFeeder->GetName(), input, client->GetName());
      return;
    } 
  }
//  if (strcmp(client->GetInputType(fFeeder->GetOutputType(fOutput))) {
//    Fatal("Cannot connect output slot #%d (type %s) of device %s to input slot #%d of device %s. Aborting",
//            fOutput, fFeeder->GetInputType(fOutput), fFeeder->GetName(), input, client->GetName());
//  }          
  if (!fArraySize) {
    fArraySize = 8;
    fInputs = new Int_t[fArraySize];
    fDevices = new TObjArray(fArraySize);
  }
  if (fNclients >= fArraySize) {
    fArraySize *= 2;
    Int_t *array = new Int_t[fArraySize];
    memcpy(array, fInputs, fNclients*sizeof(Int_t));
    delete [] fInputs;
    fInputs = array;
  }
  fInputs[fNclients] = input;
  fDevices->Add(client);
  fNclients++;
}    

//______________________________________________________________________________
void AliTrigConnector::Print(Option_t */*option*/) const
{
// Print info about this connector.
  Printf("   feeder: output #%d of device %s\n", fOutput, fFeeder->GetName());
  Printf("   client devices:\n");
  for (Int_t i=0; i<fNclients; i++) Printf("      #%d %s\n", fInputs[i], fDevices->At(i)->GetName());
}

//______________________________________________________________________________
Bool_t AliTrigConnector::Transmit(Bool_t value)
{
// Transmit Boolean signal from feeder to all clients.
  AliTrigDevice *nextclient;
  Bool_t transmit = kTRUE;
  for (Int_t i=0; i<fNclients; i++) {
    nextclient = (AliTrigDevice*)fDevices->At(i);
    Bool_t done = nextclient->SetInputValue(fInputs[i], value);
    if (!done) {
      Error("Transmit", "Connector %s: Boolean value cannot be transmitted to input %d of device %s",
            GetName(), i,  nextclient->GetName());
      transmit = kFALSE;      
    }
  }
  return transmit;
}

//______________________________________________________________________________
Bool_t AliTrigConnector::Transmit(AliTrigEvent *event)
{
// Transmit Boolean signal from feeder to all clients.
  AliTrigDevice *nextclient;
  Bool_t transmit = kTRUE;
  for (Int_t i=0; i<fNclients; i++) {
    nextclient = (AliTrigDevice*)fDevices->At(i);
    Bool_t done = nextclient->SetInputValue(fInputs[i], event);
    if (!done) {
      Error("Transmit", "Connector %s: Event cannot be transmitted to input %d of device %s",
            GetName(), i,  nextclient->GetName());
      transmit = kFALSE;      
    }
  }
  return transmit;    
}
