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
AliMUON2DMapIterator::AliMUON2DMapIterator(AliMpExMap* theMap)
: TIterator(), 
fMap(theMap),
fCurrentMap(0x0),
fI(-1),
fJ(-1)
{
  /// default ctor
  Reset();
}

//_____________________________________________________________________________
AliMUON2DMapIterator::AliMUON2DMapIterator(const AliMUON2DMapIterator& rhs)
: TIterator(rhs),
fMap(rhs.fMap),
fCurrentMap(rhs.fCurrentMap),
fI(rhs.fI),
fJ(rhs.fI)
{
  /// copy ctor
}

//_____________________________________________________________________________
AliMUON2DMapIterator& 
AliMUON2DMapIterator::operator=(const AliMUON2DMapIterator& rhs)
{
  /// assignment operator
  if ( this != &rhs ) 
  {
    fMap = rhs.fMap;
    fCurrentMap = rhs.fCurrentMap;
    fI = rhs.fI;
    fJ = rhs.fJ;
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
    fMap = rhs1.fMap;
    fCurrentMap = rhs1.fCurrentMap;
    fI = rhs1.fI;
    fJ = rhs1.fJ;
  }
  return *this;
}

//_____________________________________________________________________________
AliMUON2DMapIterator::~AliMUON2DMapIterator()
{
  /// dtor
}

//_____________________________________________________________________________
const TCollection* 
AliMUON2DMapIterator::GetCollection() const
{
  /// Return 0 as we're not really dealing with a Root TCollection...
  return 0x0;
}

//_____________________________________________________________________________
AliMpExMap*
AliMUON2DMapIterator::Map(Int_t i) const
{
  /// Get the map at a given index
  return static_cast<AliMpExMap*>(fMap->GetObjectFast(i));
}

//_____________________________________________________________________________
TObject*
AliMUON2DMapIterator::Next()
{
  /// return next object
  
  if (!fCurrentMap) return 0x0;
  
  ++fJ;
  
  if ( fJ < fCurrentMap->GetSize() ) 
  {
    return fCurrentMap->GetObjectFast(fJ);
  }
  else
  {
    ++fI;
    if ( fI < fMap->GetSize() )
    {
      fCurrentMap = Map(fI);
      fJ = -1;
      return Next();
    }
    return 0x0;
  }
}

//_____________________________________________________________________________
void
AliMUON2DMapIterator::Reset()
{
  /// rewind the iterator
  fI = -1;
  fJ = -1;
  fCurrentMap = 0x0;
  
  if ( fMap->GetSize() > 0 )
  {
    fI = 0;
    fCurrentMap = Map(fI);
    fJ = -1;
  }
}

