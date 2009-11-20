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
//   AliTrigDigitalCircuit - Device that has N Boolean inputs and one Boolean
// output. This is a base class and derived digital circuits must implement the
// response function Trigger()
//==============================================================================

#include "AliTrigDigitalCircuit.h"

ClassImp(AliTrigDigitalCircuit)

//______________________________________________________________________________
AliTrigDigitalCircuit::AliTrigDigitalCircuit(const AliTrigDigitalCircuit &other)
                    :AliTrigDevice(other),
                     fLastOutput(other.fLastOutput),
                     fConnector(0),
                     fInputs(other.fInputs)
{
// Copy ctor.
  if (other.fConnector) fConnector = new AliTrigConnector(*other.fConnector);
}                        

//______________________________________________________________________________
AliTrigDigitalCircuit::~AliTrigDigitalCircuit()
{
// Destructor
  if (fConnector) delete fConnector;
}  

//______________________________________________________________________________
AliTrigDigitalCircuit& AliTrigDigitalCircuit::operator=(const AliTrigDigitalCircuit &other)
{
// Assignment
  if (&other == this) return *this;
  AliTrigDevice::operator=(other);
  fLastOutput = other.fLastOutput;
  if (other.fConnector) fConnector = new AliTrigConnector(*other.fConnector);
  else fConnector = 0;
  fInputs = other.fInputs;
  return *this;
}

//______________________________________________________________________________
Bool_t AliTrigDigitalCircuit::Connect(UInt_t output, AliTrigDevice *other, UInt_t at_input)
{
// Connect to an input of another device.
  if (!fConnector) fConnector = new AliTrigConnector(this, 0);
  fConnector->Connect(other, at_input);
}

//______________________________________________________________________________
Bool_t AliTrigDigitalCircuit::Response(UInt_t /*output*/)
{
// Response function of the digital circuit. Calling user-defined one.
  fLastOutput = CircuitResponse();
  if (fConnector) fConnector->Transmit(fLastOutput);
  return fLastOutput;
}   
