// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveMagField.h"
#include "AliEveEventManager.h"

#include <AliMagF.h>

//______________________________________________________________________________
// Full description of AliEveMagField
//

ClassImp(AliEveMagField)

//______________________________________________________________________________
AliEveMagField::AliEveMagField(AliMagF* mf) :
  TEveMagField(),
  fField(mf)
{
  // Constructor.

  if (fField == 0)
  {
    fField = AliEveEventManager::AssertMagField();
  }
}

//______________________________________________________________________________
TEveVector AliEveMagField::GetField(Float_t x, Float_t y, Float_t z) const
{
  // Return magnetic field at requested point.

  Double_t rb[3] = { x, y, z };
  Double_t bb[3];

  fField->Field(rb, bb);

  TEveVector b(bb);
  b *= -0.1f;
  return b;
}
