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
//     Event IO handler base class
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include "AliVEventHandler.h"

ClassImp(AliVEventHandler)

//______________________________________________________________________________
AliVEventHandler::AliVEventHandler() :
  TNamed()
{
  // default constructor
}

//______________________________________________________________________________
AliVEventHandler::~AliVEventHandler() 
{
// destructor
}

//______________________________________________________________________________
AliVEventHandler::AliVEventHandler(const char* name, const char* title):
    TNamed(name, title)
{
}
//______________________________________________________________________________
void AliVEventHandler::Lock()
{
// Security lock. This is to detect NORMAL user errors and not really to
// protect against intentional hacks.
   if (IsLocked()) return;
   TObject::SetBit(kHandlerLocked, kTRUE);
}

//______________________________________________________________________________
void AliVEventHandler::UnLock()
{
// Verbose unlocking. Hackers will be punished ;-) ... 
   if (!IsLocked()) return;
   TObject::SetBit(kHandlerLocked, kFALSE);
}

//______________________________________________________________________________
void AliVEventHandler::Changed()
{
// All critical setters pass through the Changed method that throws an exception 
// in case the lock was set.
   if (IsLocked()) Fatal("Changed","Critical setter of a handler called in locked mode");
}

//______________________________________________________________________________
const char *AliVEventHandler::GetExtraOutputs(Bool_t) const
{
// Returns extra outputs. If merge is requested, returns only files to be
// merged. Implementation in derived classes.
   return 0;
}
