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

// Author: Andrei Gheata, 28/07/2009

//==============================================================================
//   AliTrigModule - Base class for trigger devices handling generic events. 
//      A module has arbitrary number of inputs and outputs. Derived classes must 
//      implement CreateDevice() and Trigger() metods.
//==============================================================================

#include "AliTrigModule.h"

#include "TObjArray.h"
#include "AliTrigConnector.h"
#include "AliTrigEvent.h"

ClassImp(AliTrigModule)

//______________________________________________________________________________
AliTrigModule::~AliTrigModule()
{
// Destructor
  if (fInputs) {fInputs->Delete(); delete fInputs;}
  if (fOutputs) {fOutputs->Delete(); delete fOutputs;}
  if (fOutputConnectors) {fOutputConnectors->Delete(); delete fOutputConnectors;}
}  

//______________________________________________________________________________
Bool_t AliTrigModule::Connect(Int_t output, AliTrigDevice *other, Int_t at_input)
{
// Connect to an input of another device.
  if (!fOutputConnectors) fOutputConnectors = new TObjArray(fNoutputs);
  AliTrigConnector *connector = new AliTrigConnector(Form("wire_%s_%d", GetName(), output), (AliTrigDevice*)this, 0);
  connector->Connect(other, at_input);
  fOutputConnectors->AddAt(connector, output);
  return kTRUE;
}

//______________________________________________________________________________
void AliTrigModule::DefineInput(Int_t islot, AliTrigEvent *event)
{
// Define an input slot and provide an event derived from AliTrigEvent.
  if (!fInputs) fInputs = new TObjArray(fNinputs);
  fInputs->AddAt(event, islot);
}   

//______________________________________________________________________________
void AliTrigModule::DefineOutput(Int_t islot, AliTrigEvent *event)
{
// Define an output slot and provide an event derived from AliTrigEvent.
  if (!fOutputs) fOutputs = new TObjArray(fNoutputs);
  fOutputs->AddAt(event, islot);
}   

//______________________________________________________________________________
AliTrigEvent *AliTrigModule::GetInputValue(Int_t input) const
{
// Get current input value for a slot.
  if (!fInputs) return 0;
  return (AliTrigEvent*)fInputs->At(input);
}   

//______________________________________________________________________________
AliTrigEvent *AliTrigModule::GetOutputValue(Int_t output) const
{
// Get current input value for a slot.
  if (!fOutputs) return 0;
  return (AliTrigEvent*)fOutputs->At(output);
}   

//______________________________________________________________________________
Bool_t AliTrigModule::Response(Int_t ioutput)
{
// Response function of the digital circuit. Calling user-defined one.
  Bool_t response = Trigger(ioutput);
  AliTrigConnector *connector = (AliTrigConnector*)fOutputConnectors->At(ioutput);
  if (connector) connector->Transmit(GetOutputValue(ioutput));
  return response;
}   

//______________________________________________________________________________
void AliTrigModule::ResetInputs()
{
// Reset all inputs
}

//______________________________________________________________________________
Bool_t AliTrigModule::SetInputValue(Int_t input, AliTrigEvent *event)
{
// A way to set directly the input value.
  if (!fInputs) return kFALSE;
  AliTrigEvent *input_event = GetInputValue(input);
  if (!input_event) return kFALSE;
  input_event->ImportData(event);
  return kTRUE;
}
