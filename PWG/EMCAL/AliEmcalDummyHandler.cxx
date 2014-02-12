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

//-------------------------------------------------------------------------
//     Dummy implementation of the input handler e.g. for the case of event-loop steered analysis
//     Author: M .Verweij
//-------------------------------------------------------------------------

#include "AliEmcalDummyHandler.h"
#include <TTree.h>

ClassImp(AliEmcalDummyHandler)

//______________________________________________________________________________
AliEmcalDummyHandler::AliEmcalDummyHandler() :
  AliDummyHandler(),
  fEvent(0)

{
  // Default constructor
}

//______________________________________________________________________________
AliEmcalDummyHandler::AliEmcalDummyHandler(const char* name, const char* title):
  AliDummyHandler(name, title),
  fEvent(0)
{
   // Constructor
}

//______________________________________________________________________________
AliEmcalDummyHandler::~AliEmcalDummyHandler() 
{
// Destructor
}
