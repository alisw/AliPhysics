#ifndef ALITRIGCONNECTOR_H
#define ALITRIGCONNECTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 28/07/2009

//==============================================================================
//   AliTrigConnector - Class representing a connector between an output of a 
// device (feeder) and an arbitrary number of inputs of other devices.
//==============================================================================

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TObjArray;
class AliTrigDevice;
class AliTrigEvent;

class AliTrigConnector : public TNamed {

public:
  AliTrigConnector() : TNamed(), fFeeder(0), fOutput(0), fNclients(0), fArraySize(0), fInputs(0), fDevices(0) {}
  AliTrigConnector(const char *name, AliTrigDevice *feeder, Int_t output) : TNamed(name, ""), fFeeder(feeder), fOutput(output), fNclients(0), fArraySize(0), fInputs(0), fDevices(0) {}
  virtual ~AliTrigConnector();

  
  // Connect a client input.
  void                      Connect(AliTrigDevice *client, Int_t input);

  virtual void              Print(Option_t *option="") const;  
  
  // Transmit the feeder signal to all connected inputs. Different device types
  // call different Transmit() methods.
  Bool_t                    Transmit(Bool_t value);
  Bool_t                    Transmit(AliTrigEvent *event);
  
private:
  AliTrigConnector(const AliTrigConnector &other);
  AliTrigConnector &operator=(const AliTrigConnector &other);

  AliTrigDevice            *fFeeder;    // Feeder device
  Int_t                     fOutput;    // Output slot index for the feeder
  Int_t                     fNclients;  // Number of clients
  Int_t                     fArraySize; // Size of the clients array
  Int_t                    *fInputs;    //[fArraySize] Array of input slot indices
  TObjArray                *fDevices;   // Array of client devices
   
  ClassDef(AliTrigConnector,1)  // Class representing a connector between devices.
};
#endif
