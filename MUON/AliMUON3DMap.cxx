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

// $Id$

#include "AliMUON3DMap.h"

ClassImp(AliMUON3DMap)

#include "AliMUON1DMap.h"
#include "AliMUON2DMap.h"

#include "Riostream.h"

//_____________________________________________________________________________
AliMUON3DMap::AliMUON3DMap() : AliMUONV3DStore(), fStore(new AliMUON1DMap)
{
}

//_____________________________________________________________________________
AliMUON3DMap::~AliMUON3DMap()
{
  delete fStore;
}

//_____________________________________________________________________________
TObject* 
AliMUON3DMap::Get(Int_t i, Int_t j, Int_t k) const
{
  AliMUONV2DStore* m = static_cast<AliMUONV2DStore*>(fStore->Get(i));
  if (!m) return 0x0;
  return m->Get(j,k);
}

//_____________________________________________________________________________
Bool_t 
AliMUON3DMap::IsOwner() const
{
  return kTRUE;
}

//_____________________________________________________________________________
void
AliMUON3DMap::Print(Option_t*) const
{
  cout << "Would need an iterator here to be able to print !" << endl;
}

//_____________________________________________________________________________
Bool_t 
AliMUON3DMap::Set(Int_t i, Int_t j, Int_t k, TObject* object, Bool_t replace)
{
  AliMUONV2DStore* m = static_cast<AliMUONV2DStore*>(fStore->Get(i));
  if (!m)
  {
    fStore->Set(i,new AliMUON2DMap,replace);
    m = static_cast<AliMUONV2DStore*>(fStore->Get(i));
  }
  return m->Set(j,k,object,replace);
}



