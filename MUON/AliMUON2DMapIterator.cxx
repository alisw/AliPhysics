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

/// \cond CLASSIMP
ClassImp(AliMUON2DMapIterator)
/// \endcond

//_____________________________________________________________________________
AliMUON2DMapIterator::AliMUON2DMapIterator(const AliMpExMap& theMap)
: TIterator(), 
fIter(theMap.GetIterator()),
fIter2(0x0),
fCurrentI(-1),
fCurrentJ(-1)
{
  /// default ctor
  Reset();
}

//_____________________________________________________________________________
AliMUON2DMapIterator::AliMUON2DMapIterator(const AliMUON2DMapIterator& rhs)
: TIterator(rhs),
fIter(rhs.fIter),
fIter2(0x0),
fCurrentI(rhs.fCurrentI),
fCurrentJ(rhs.fCurrentJ)
{
  /// copy ctor
  if ( rhs.fIter2 ) fIter2 = new TExMapIter(*rhs.fIter2);

}
//_____________________________________________________________________________
AliMUON2DMapIterator& 
AliMUON2DMapIterator::operator=(const AliMUON2DMapIterator& rhs)
{
  /// assignment operator
  if ( this != &rhs ) 
  {
    fIter = rhs.fIter;
    fIter2 = 0x0;
    if ( rhs.fIter2 ) fIter2 = new TExMapIter(*rhs.fIter2);
    fCurrentI = rhs.fCurrentI;
    fCurrentJ = rhs.fCurrentJ;
  }
  return *this;
}

//_____________________________________________________________________________
TIterator& 
AliMUON2DMapIterator::operator=(const TIterator& rhs)
{
  /// overriden operator= (imposed by Root's definition of TIterator::operator= ?)
  
  if ( this != &rhs && rhs.IsA() == AliMUON2DMapIterator::Class() ) 
  {
    const AliMUON2DMapIterator& rhs1 = static_cast<const AliMUON2DMapIterator&>(rhs);
    
    fIter = rhs1.fIter;
    fIter2 = 0x0;
    if ( rhs1.fIter2 ) fIter2 = new TExMapIter(*rhs1.fIter2);
    fCurrentI = rhs1.fCurrentI;
    fCurrentJ = rhs1.fCurrentJ;
  }
  return *this;
}

//_____________________________________________________________________________
AliMUON2DMapIterator::~AliMUON2DMapIterator()
{
  /// dtor
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
TObject*
AliMUON2DMapIterator::GetValue(TExMapIter& iter, Int_t& theKey) const
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
AliMpExMap*
AliMUON2DMapIterator::GetMap(TExMapIter& iter, Int_t& key)  const
{
  /// get the map corresponding to key
  AliMpExMap* rv(0x0);
  TObject* o = GetValue(iter,key);
  if (o)
  {
    rv = dynamic_cast<AliMpExMap*>(o);
  }
  return rv;
}

//_____________________________________________________________________________
TObject*
AliMUON2DMapIterator::Next()
{
  /// logic :
  /// get TObject* from fIter2
  /// if null, increment fIter2 by getting next on fIter
  
  if (!fIter2) return 0x0;
  
  TObject* o = GetValue(*fIter2,fCurrentJ);
  if (!o)
  {
    // fIter2 exhausted, try to get the next one
    AliMpExMap* m = GetMap(fIter,fCurrentI);
    if (!m)
    {
      // nothing left, we are done.
      return 0x0;
    }
    delete fIter2;
    fIter2 = new TExMapIter(m->GetIterator());
    o = GetValue(*fIter2,fCurrentJ);
  }
  
  return o;
}

//_____________________________________________________________________________
void
AliMUON2DMapIterator::Reset()
{
  /// rewind the iterator
  delete fIter2;
  fIter.Reset();
  AliMpExMap* m = GetMap(fIter,fCurrentI);
  if (m)
  {
    fIter2 = new TExMapIter(m->GetIterator());
  }  
}
