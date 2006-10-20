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

#include "AliMUONCheckItemIterator.h"
#include "TExMap.h"
#include "AliMpExMap.h"
#include "AliMUONCheckItem.h"

/// \class AliMUONCheckItemIterator
///
/// Iterator on AliMUONCheckItem objects
/// 
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONCheckItemIterator)
/// \endcond

//_____________________________________________________________________________
AliMUONCheckItemIterator::AliMUONCheckItemIterator() : TObject(), fIter(0x0)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONCheckItemIterator::AliMUONCheckItemIterator(const AliMUONCheckItem& item)
: TObject(),
fIter(0x0)
{
  /// ctor
  AliMpExMap* m = item.fMissing;
  fIter = new TExMapIter(m->GetIterator());
}

//_____________________________________________________________________________
AliMUONCheckItemIterator::~AliMUONCheckItemIterator()
{
  /// dtor
  delete fIter;
}

//_____________________________________________________________________________
void
AliMUONCheckItemIterator::First()
{
  /// Rewind the iterator
  if ( fIter) fIter->Reset();
}

//_____________________________________________________________________________
TObject*
AliMUONCheckItemIterator::Next()
{
  /// Advance one object. Return 0 if ended.
  if (!fIter) return 0x0;
  Long_t key, value;
  Bool_t ok = fIter->Next(key,value);
  if (ok)
  {
    return reinterpret_cast<TObject*>(value);
  }
  else
  {
    return 0x0;
  }
}
