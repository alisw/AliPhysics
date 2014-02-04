// @(#)root/eve:$Id$
// Author: Matevz Tadel 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <TError.h>
#include <TSystem.h>

#include "AliDimIntNotifier.h"

//______________________________________________________________________________
// Full description of AliDimIntNotifier
//

ClassImp(AliDimIntNotifier)

AliDimIntNotifier::AliDimIntNotifier(const TString& service) :
  DimUpdatedInfo(service, -1),
  fLastMessage(-1)
{

}

void AliDimIntNotifier::infoHandler()
{
	// Handle DIM message
	fLastMessage = getData() ? getInt() : -1;
	DimMessage(fLastMessage);
}

void AliDimIntNotifier::DimMessage(Int_t)
{

  if (fLastMessage != -1)
  {
    Emit("DimMessage(Int_t)", fLastMessage);
  }
  gSystem->ProcessEvents();
}
