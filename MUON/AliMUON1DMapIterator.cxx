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

#include "AliMUON1DMapIterator.h"
#include "AliMpExMap.h"
#include "AliLog.h"
#include "AliMpIntPair.h"
#include "AliMUONObjectPair.h"
#include "AliMpExMap.h"

/// \class AliMUON1DMapIterator
/// Implementation of AliMUONVDataIterator for 1Dmaps
/// 
/// A simple implementation of VDataIterator for 1Dmaps.
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUON1DMapIterator)
/// \endcond

//_____________________________________________________________________________
AliMUON1DMapIterator::AliMUON1DMapIterator(AliMpExMap& theMap)
: AliMUONVDataIterator(), 
fIter(theMap.GetIterator()),
fCurrentI(-1)
{
  /// default ctor
  Reset();
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
  ///  

  Long_t key, value;
  Bool_t ok = fIter.Next(key,value);
  if (!ok) return 0x0;
  Int_t i = (Int_t)(key & 0xFFFF);
  TObject* o = reinterpret_cast<TObject*>(value);
  
  return new AliMUONObjectPair(new AliMpIntPair(i,0),o,
                               kTRUE, /* owner of intpair */
                               kFALSE /* but not of o */);
}

//_____________________________________________________________________________
void
AliMUON1DMapIterator::Reset()
{
  /// rewind the iterator
  fIter.Reset();
}

//_____________________________________________________________________________
Bool_t
AliMUON1DMapIterator::Remove()
{
  /// to be implemented if needed
  AliInfo("Not supported yet");
  return kFALSE;
}
