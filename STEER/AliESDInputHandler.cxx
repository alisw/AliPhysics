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
//     Event handler for ESD input 
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TChain.h>

#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESD.h"

ClassImp(AliESDInputHandler)

//______________________________________________________________________________
AliESDInputHandler::AliESDInputHandler() :
  AliInputEventHandler()
{
  // default constructor
}

//______________________________________________________________________________
AliESDInputHandler::~AliESDInputHandler() 
{
// destructor
}

//______________________________________________________________________________
AliESDInputHandler::AliESDInputHandler(const char* name, const char* title):
    AliInputEventHandler(name, title)
{
}

Bool_t AliESDInputHandler::InitIO(Option_t* /*opt*/)
{
    // Get pointer to ESD event
    fEvent = new AliESDEvent();
    fEvent->ReadFromTree(fChain);
    return kTRUE;
}

Bool_t AliESDInputHandler::BeginEvent()
{
    // Copy from old to new format if necessary
    AliESD* old = ((AliESDEvent*) fEvent)->GetAliESDOld();
    if (old) {
	printf("It's an old ESD \n");
	
	((AliESDEvent*)fEvent)->CopyFromOldESD();
	old->Reset();
    }
    return kTRUE;
}

