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

#include "AliMUON2DMapIterator.h"

//-----------------------------------------------------------------------------
/// \class AliMUON2DMapIterator
/// Implementation of TIterator for 2Dmaps
/// 
/// A simple implementation of VDataIterator for 2Dmaps.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMpExMap.h"
#include "AliMpExMapIterator.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUON2DMapIterator)
/// \endcond

//_____________________________________________________________________________
AliMUON2DMapIterator::AliMUON2DMapIterator(const AliMpExMap& theMap)
: TIterator(),
fkMap(&theMap),
fIter1(theMap.CreateIterator()),
fIter2(NextIterator())
{
  /// default ctor
  Reset();
}

//_____________________________________________________________________________
AliMUON2DMapIterator& 
AliMUON2DMapIterator::operator=(const TIterator& /*rhs*/)
{
  // overriden operator= (imposed by Root's definition of TIterator::operator= ?)
  
  AliFatalGeneral("operator=(TIterator&)",""); // as in copy ctor
  return *this;
}

//_____________________________________________________________________________
AliMUON2DMapIterator::~AliMUON2DMapIterator()
{
  /// dtor
  delete fIter1;
  delete fIter2;
}

//_____________________________________________________________________________
const TCollection* 
AliMUON2DMapIterator::GetCollection() const
{
  /// Return 0 as we're not really dealing with a Root TCollection...
  return 0x0;
}

//_____________________________________________________________________________
TIterator*
AliMUON2DMapIterator::NextIterator()
{
  /// Get next map (from fIter1) and create an iterator to it
  
  AliMpExMap* m = static_cast<AliMpExMap*>(fIter1->Next());

  if (!m) return 0x0;
  
  return m->CreateIterator();
}

//_____________________________________________________________________________
TObject*
AliMUON2DMapIterator::Next()
{
  /// return next object
  
  if (!fIter2) return 0x0;

  TObject* o = fIter2->Next();
  
  if (!o)
  {
    delete fIter2;
    fIter2 = NextIterator();
    return Next();
  }
  
  return o;
}

//_____________________________________________________________________________
void
AliMUON2DMapIterator::Reset()
{
  /// rewind the iterator
  
  delete fIter2;
  fIter1->Reset();
  fIter2 = NextIterator();
}

