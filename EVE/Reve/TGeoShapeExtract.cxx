// $Header$

// Copyright (C) 1999-2005, Matevz Tadel. All rights reserved.
// This file is part of GLED, released under GNU General Public License version 2.
// For the licensing terms see $GLEDSYS/LICENSE or http://www.gnu.org/.

//__________________________________________________________________________
// TGeoShapeExtract
//
// Vessel to carry hand-picked geometry from gled to reve.
// This class exists in both frameworks.

#include "TGeoShapeExtract.h"

#include <TList.h>
#include <TGeoShape.h>

ClassImp(TGeoShapeExtract)

/**************************************************************************/

TGeoShapeExtract::TGeoShapeExtract(const Text_t* n, const Text_t* t) :
  TNamed(n,t),
  mRnrSelf     (true),
  mRnrElements (true),
  mShape       (0),
  mElements    (0)
{
  memset(mTrans, 0, sizeof(mTrans));
  mTrans[0] = mTrans[5] = mTrans[10] = mTrans[15] = 1;
  mRGBA [0] = mRGBA [1] = mRGBA [2]  = mRGBA [3]  = 1;
}

TGeoShapeExtract::~TGeoShapeExtract()
{
  delete mShape;
  delete mElements;
}

/**************************************************************************/

Bool_t TGeoShapeExtract::HasElements()
{
  return mElements != 0 && mElements->GetSize() > 0;
}

void TGeoShapeExtract::AddElement(TGeoShapeExtract* gse)
{
  if (mElements == 0)
    mElements = new TList;

  mElements->Add(gse);
}

/**************************************************************************/

void TGeoShapeExtract::SetTrans(const Double_t arr[16])
{
  for(Int_t i=0; i<16; ++i)
    mTrans[i] = arr[i];
}

void TGeoShapeExtract::SetRGBA (const Float_t  arr[4])
{
  for(Int_t i=0; i<4; ++i)
    mRGBA[i] = arr[i];
}
