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

//-----------------------------------------------------------------------------
/// \class AliMUON2DMapIteratorByI
///
/// Implementation of TIterator for 2D maps
/// 
/// An implementation of TIterator for 2D maps, which can iterate
/// on a range of i values (i being the first element of the couple
/// (i,j) used to index values in the map).
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUON2DMapIteratorByI.h"
#include "AliMpExMapIterator.h"
#include "AliMpExMap.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUON2DMapIteratorByI)
/// \endcond

//_____________________________________________________________________________
AliMUON2DMapIteratorByI::AliMUON2DMapIteratorByI(const AliMpExMap& theMap, Int_t firstI, Int_t lastI)
: TIterator(),
fkMap(&theMap),
fIter1(theMap.CreateIterator()),
fIter2(0x0),
fFirstI(firstI),
fLastI(lastI),
fCurrentI(-1)
{
  /// default ctor
  Reset();
}

//_____________________________________________________________________________
AliMUON2DMapIteratorByI& 
AliMUON2DMapIteratorByI::operator=(const TIterator& /*rhs*/)
{
  // overriden operator= (imposed by Root's definition of TIterator::operator= ?)
  
  AliFatalGeneral("operator=(TIterator&)",""); // as in copy ctor
  return *this;
}

//_____________________________________________________________________________
AliMUON2DMapIteratorByI::~AliMUON2DMapIteratorByI()
{
  /// dtor
  delete fIter1;
  delete fIter2;
}

//_____________________________________________________________________________
const TCollection* 
AliMUON2DMapIteratorByI::GetCollection() const
{
  /// Return 0 as we're not really dealing with a Root TCollection...
  return 0x0;
}

//_____________________________________________________________________________
AliMpExMapIterator*
AliMUON2DMapIteratorByI::NextIterator()
{
  /// Get next map (from fIter1) and create an iterator to it
  
  AliMpExMap* m = static_cast<AliMpExMap*>(fIter1->Next(fCurrentI));
  
  if (!m) return 0x0;
  
  if ( fCurrentI < fFirstI || fCurrentI > fLastI ) return NextIterator(); // try again

  return m->CreateIterator();
}

//_____________________________________________________________________________
TObject*
AliMUON2DMapIteratorByI::Next()
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
AliMUON2DMapIteratorByI::Reset()
{
  /// rewind the iterator
  
  delete fIter2;
  fIter1->Reset();
  fIter2 = NextIterator();
  fCurrentI = -1;
}

