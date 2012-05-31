// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveHOMERSource.h"

//______________________________________________________________________________
// AliEveHOMERSource
//

ClassImp(AliEveHOMERSource)

AliEveHOMERSource::AliEveHOMERSource(const Text_t* n, const Text_t* t) :
  TEveElement (),
  TNamed      (n, t),
  fSrcId      (0),
  fSrcState   (0)
{}

/******************************************************************************/

//______________________________________________________________________________
Bool_t AliEveHOMERSource::SetRnrState(Bool_t rnr)
{
   // Set render state of this element and of its children to the same
   // value.

  if (fSrcState)
    fSrcState->fState = rnr;

  return TEveElement::SetRnrState(rnr);
}
