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
/// \class AliMUON1DMapIterator
/// Implementation of TIterator for 1Dmaps
/// 
/// A simple implementation of VDataIterator for 1Dmaps.
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUON1DMapIterator.h"
#include "AliMpExMap.h"
#include "AliLog.h"
#include "AliMpIntPair.h"
#include "AliMpExMap.h"

/// \cond CLASSIMP
ClassImp(AliMUON1DMapIterator)
/// \endcond

//_____________________________________________________________________________
AliMUON1DMapIterator::AliMUON1DMapIterator(AliMpExMap& theMap)
: TIterator(), 
fIter(theMap.GetIterator()),
fCurrentI(-1)
{
  /// default ctor
  Reset();
}

//_____________________________________________________________________________
AliMUON1DMapIterator::AliMUON1DMapIterator(const AliMUON1DMapIterator& rhs)
: TIterator(rhs),
  fIter(rhs.fIter),
  fCurrentI(rhs.fCurrentI)        
{
    /// copy ctor
}

//_____________________________________________________________________________
AliMUON1DMapIterator& 
AliMUON1DMapIterator::operator=(const AliMUON1DMapIterator& rhs)
{
  /// assignment operator
  if ( this != &rhs )
  {
    fIter = rhs.fIter;
    fCurrentI = rhs.fCurrentI;
  }
  return *this;
}

//_____________________________________________________________________________
TIterator& 
AliMUON1DMapIterator::operator=(const TIterator& rhs)
{
  /// overriden operator= (imposed by Root definition of TIterator::operator= ?)
  
  if ( this != &rhs && rhs.IsA() == AliMUON1DMapIterator::Class() )
  {
    const AliMUON1DMapIterator& rhs1 = static_cast<const AliMUON1DMapIterator&>(rhs);
    fIter = rhs1.fIter;
    fCurrentI = rhs1.fCurrentI;
  }
  return *this;
}


//_____________________________________________________________________________
AliMUON1DMapIterator::~AliMUON1DMapIterator()
{
  /// dtor
}

//_____________________________________________________________________________
TObject*
AliMUON1DMapIterator::Next()
{
  /// Return next object in iteration

  Long_t key, value;
  Bool_t ok = fIter.Next(key,value);
  if (!ok) return 0x0;
  return reinterpret_cast<TObject*>(value);
}

//_____________________________________________________________________________
void
AliMUON1DMapIterator::Reset()
{
  /// rewind the iterator
  fIter.Reset();
}
