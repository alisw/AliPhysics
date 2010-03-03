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
#include "AliTrigConnector.h"

ClassImp(AliTrigDigitalCircuit)

//______________________________________________________________________________
AliTrigDigitalCircuit::~AliTrigDigitalCircuit()
{
// Destructor
  if (fConnector) delete fConnector;
}  

//______________________________________________________________________________
Bool_t AliTrigDigitalCircuit::Connect(Int_t output, AliTrigDevice *other, Int_t at_input)
{
// Connect to an input of another device.
  if (!fConnector) fConnector = new AliTrigConnector(Form("wire_%s_%d", GetName(), output), (AliTrigDevice*)this, 0);
  fConnector->Connect(other, at_input);
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliTrigDigitalCircuit::Response(Int_t /*output*/)
{
// Response function of the digital circuit. Calling user-defined one.
  fLastOutput = Trigger();
  if (fConnector) fConnector->Transmit(fLastOutput);
  return fLastOutput;
}   
