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

#include "AliMpExMap.h"

/// \cond CLASSIMP
ClassImp(AliMUON2DMapIteratorByI)
/// \endcond

//_____________________________________________________________________________
AliMUON2DMapIteratorByI::AliMUON2DMapIteratorByI(const AliMpExMap& theMap,
                                                 Int_t firstI,
                                                 Int_t lastI)
: TIterator(),
fkMap(&theMap),
fIter2(0x0),
fCurrentI(-1),
fCurrentJ(-1),
fFirstI(firstI),
fLastI(lastI)
{
  /// default ctor
  Reset();
}

//_____________________________________________________________________________
AliMUON2DMapIteratorByI::AliMUON2DMapIteratorByI(const AliMUON2DMapIteratorByI& rhs)
:TIterator(rhs),
fkMap(rhs.fkMap),
fIter2(0x0),
fCurrentI(rhs.fCurrentI),
fCurrentJ(rhs.fCurrentJ),
fFirstI(rhs.fFirstI),
fLastI(rhs.fLastI)
{
  /// copy ctor
  if ( rhs.fIter2 ) fIter2 = new TExMapIter(*(rhs.fIter2));
}

//_____________________________________________________________________________
AliMUON2DMapIteratorByI&
AliMUON2DMapIteratorByI::operator=(const AliMUON2DMapIteratorByI& rhs)
{
  /// assignment operator
  if ( this != &rhs ) 
  {
    fkMap = rhs.fkMap;
    fIter2 = 0x0;
    if ( rhs.fIter2 ) fIter2 = new TExMapIter(*(rhs.fIter2));
    fCurrentI = rhs.fCurrentI;
    fCurrentJ = rhs.fCurrentJ;
    fFirstI = rhs.fFirstI;
    fLastI = rhs.fLastI;
  }
  return *this;
}

//_____________________________________________________________________________
TIterator&
AliMUON2DMapIteratorByI::operator=(const TIterator& rhs)
{
  /// overriden assigment operator (imposed by Root's declaration of TIterator ?)
  if ( this != &rhs && rhs.IsA() == AliMUON2DMapIteratorByI::Class() ) 
  {
    const AliMUON2DMapIteratorByI& rhs1 = static_cast<const AliMUON2DMapIteratorByI&>(rhs);
    fkMap = rhs1.fkMap;
    fIter2 = 0x0;
    if ( rhs1.fIter2 ) fIter2 = new TExMapIter(*(rhs1.fIter2));
    fCurrentI = rhs1.fCurrentI;
    fCurrentJ = rhs1.fCurrentJ;
    fFirstI = rhs1.fFirstI;
    fLastI = rhs1.fLastI;
  }
  return *this;
}

//_____________________________________________________________________________
AliMUON2DMapIteratorByI::~AliMUON2DMapIteratorByI()
{
  /// dtor
  delete fIter2;
}

//_____________________________________________________________________________
const TCollection* 
AliMUON2DMapIteratorByI::GetCollection() const
{
  /// Not implemented
  return 0x0;
}

//_____________________________________________________________________________
TObject*
AliMUON2DMapIteratorByI::GetValue(TExMapIter& iter, Int_t& theKey) const
{
  /// return the value corresponding to theKey in iterator iter
  theKey = -1;
  Long_t key, value;
  Bool_t ok = iter.Next(key,value);
  if (!ok) return 0x0;
  theKey = (Int_t)(key & 0xFFFF);
  return reinterpret_cast<TObject*>(value);
}

//_____________________________________________________________________________
TObject*
AliMUON2DMapIteratorByI::Next()
{
  /// logic :
  /// get TObject* from fIter2
  /// if null, increment fIter2 by getting next map from fMap
  
  if (!fIter2) return 0x0;
  
  TObject* o = GetValue(*fIter2,fCurrentJ);
  if (!o)
  {
    // fIter2 exhausted, try to get the next one
    delete fIter2;
    fIter2 = 0x0;
    AliMpExMap* m(0x0);
    while ( !m && fCurrentI < fLastI ) 
    {
      ++fCurrentI;
      m = static_cast<AliMpExMap*>(fkMap->GetValue(fCurrentI));
    }
    if (!m) return 0x0; // we are done
    fIter2 = new TExMapIter(m->GetIterator());
    o = GetValue(*fIter2,fCurrentJ);
  }
  
  return o;
}

//_____________________________________________________________________________
void
AliMUON2DMapIteratorByI::Reset()
{
  /// rewind the iterator
  delete fIter2;
  fIter2 = 0x0;
  fCurrentI = fFirstI;
  AliMpExMap* m;
  
  while ( !(  m = static_cast<AliMpExMap*>(fkMap->GetValue(fCurrentI) ) ) && 
          fCurrentI < fLastI )
  {
    ++fCurrentI;
  }
  
  if ( m ) 
  {
    fIter2 = new TExMapIter(m->GetIterator());
  }
}
