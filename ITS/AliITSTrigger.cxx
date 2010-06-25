/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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

////////////////////////////////////////////////////////////////////////
//                                                                    //
// Simulates generation of Fast-OR signals from SPD (if needed).      //
// Processes the Fast-OR signals generated in AliITSsimulationSPD.    //
// Provides inputs for AliCentralTrigger.                             //
//                                                                    //
// Version 2, Henrik Tydesjo, Feb 2009                                //
// Version 1, D. Elia, C. Jorgensen, Mar 2006                         //
// Version 0, J. Conrad, E. Lopez Torres, Oct 2005                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "AliITSTrigger.h"
#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliITSLoader.h"
#include "AliTriggerInput.h"

ClassImp(AliITSTrigger)

//______________________________________________________________________
AliITSTrigger::AliITSTrigger() : 
  AliTriggerDetector(),
  fPITprocessor()
{
  //standard constructor
  SetName("ITS");
}
//______________________________________________________________________
AliITSTrigger::AliITSTrigger(AliITSTriggerConditions* cond) : 
  AliTriggerDetector(),
  fPITprocessor(cond)
{
  // optional constructor
  SetName("ITS");
}
//______________________________________________________________________
void AliITSTrigger::SetTriggerConditions(AliITSTriggerConditions* cond) {
  // Sets the trigger conditions, normally coming from OCDB
  fPITprocessor.SetTriggerConditions(cond);
}
//______________________________________________________________________
void AliITSTrigger::CreateInputs() {
  // Create inputs, based on OCDB Pixel Trigger Conditions
  if( fInputs.GetEntriesFast() > 0 ) return; // Inputs already created, no need to proceed
  
  // Load trigger conditions from OCDB if needed
  if (! fPITprocessor.TriggerConditionsSet() ) {
    AliError("Trigger conditions not set. No inputs created.");
    return;
  }

  UInt_t numInputs = fPITprocessor.GetNumOutputs();
  AliInfo(Form("Number of trigger inputs: %d",numInputs));
  for (UInt_t inp=0; inp<numInputs; inp++) {
    fInputs.AddLast( new AliTriggerInput(fPITprocessor.GetOutputLabel(inp), "SPD", 0) );
  }
}
//______________________________________________________________________
void AliITSTrigger::Trigger() {
  // Performs Pixel Trigger processing of the simulated fast-or signals

  // Get the FO signals for this event
  AliITSFOSignalsSPD* foSignals = NULL;
  AliRunLoader* runLoader = AliRunLoader::Instance();
  AliITSLoader* itsLoader = (AliITSLoader*) runLoader->GetLoader("ITSLoader");
  if (!itsLoader) {
    AliError("ITS loader is NULL.");
  }

   else {
      itsLoader->LoadDigits();
      TTree *tree = itsLoader->TreeD();
      if(!tree) {
        AliError("TreeD not available");
	itsLoader->UnloadDigits();
        return;
      }
      foSignals = (AliITSFOSignalsSPD*)tree->GetUserInfo()->FindObject("AliITSFOSignalsSPD");
      if(!foSignals) AliError("FO signals not retrieved");
     }

  // Process the FO signals
  if (foSignals) {
    fPITprocessor.PreprocessFOSignals(foSignals);
    UInt_t numInputs = fPITprocessor.GetNumOutputs();
    for (UInt_t inp=0; inp<numInputs; inp++) {
      if (fPITprocessor.ProcessFOSignalsIndex(inp, foSignals)) {
	SetInput(fPITprocessor.GetOutputLabel(inp));
      }
    }
  }
  else {
    AliError("Fast-OR signals not available. No trigger processing done.");
  }
  if (itsLoader) itsLoader->UnloadDigits();
}
