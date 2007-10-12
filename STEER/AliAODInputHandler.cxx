/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//     Event handler for AOD input 
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TChain.h>

#include "AliAODInputHandler.h"
#include "AliAODEvent.h"

ClassImp(AliAODInputHandler)

//______________________________________________________________________________
AliAODInputHandler::AliAODInputHandler() :
  AliInputEventHandler()
{
  // default constructor
}

//______________________________________________________________________________
AliAODInputHandler::~AliAODInputHandler() 
{
// destructor
}

//______________________________________________________________________________
AliAODInputHandler::AliAODInputHandler(const char* name, const char* title):
    AliInputEventHandler(name, title)
{
}

Bool_t AliAODInputHandler::InitIO(Option_t* /*opt*/)
{
    // Get pointer to AOD event
    fEvent = new AliAODEvent();
    fEvent->ReadFromTree(fChain);
    return kTRUE;
}

