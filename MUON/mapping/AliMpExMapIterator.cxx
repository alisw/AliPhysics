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
// Class AliMpExMapIterator
// ------------------------
// Implementation of TIterator for AliMpExMap
// Author: Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMpExMapIterator.h"
#include "AliMpExMap.h"

#include "AliLog.h"

#include <TClass.h>
#include <TExMap.h>
#include <TString.h>

/// \cond CLASSIMP
ClassImp(AliMpExMapIterator)
/// \endcond

//_____________________________________________________________________________
AliMpExMapIterator::AliMpExMapIterator(const AliMpExMap& theMap)
: TIterator(),
  fIterator(new TExMapIter(&(theMap.fMap)))
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMpExMapIterator::AliMpExMapIterator(const AliMpExMapIterator& rhs)
: TIterator(rhs),
  fIterator(rhs.fIterator)
{
/// Copy constructor
}

//_____________________________________________________________________________
AliMpExMapIterator& 
AliMpExMapIterator::operator=(const AliMpExMapIterator& rhs)
{
/// Assignment operator

  if  ( this != &rhs )
  {
    fIterator = rhs.fIterator;
  }
  return *this;
}

//_____________________________________________________________________________
AliMpExMapIterator& 
AliMpExMapIterator::operator=(const TIterator& rhs)
{
  /// Overriden operator= (imposed by Root's definition of TIterator::operator= ?)
  
  if ( this != &rhs && rhs.IsA() == AliMpExMapIterator::Class() ) 
  {
    const AliMpExMapIterator& rhs1 = static_cast<const AliMpExMapIterator&>(rhs);
    fIterator = rhs1.fIterator;
  }
  return *this;
}

//_____________________________________________________________________________
AliMpExMapIterator::~AliMpExMapIterator()
{
/// Destructor

  delete fIterator;
}

//_____________________________________________________________________________
#if ROOT_VERSION_CODE >= 334081
//#if ROOT_VERSION_CODE >= 333824  // needed with Root v5.24.00-patches 
Bool_t 
AliMpExMapIterator::Next(Long64_t& index, TObject*& object)
#else
Bool_t 
AliMpExMapIterator::Next(Long_t& index, TObject*& object)
#endif
{
/// Move to next object in iteration

#if ROOT_VERSION_CODE >= 334081
//#if ROOT_VERSION_CODE >= 333824  // needed with Root v5.24.00-patches 
  Long64_t value(0);
#else
  Long_t value(0);
#endif

  object = 0;
  
  Bool_t rv = fIterator->Next(index,value);

  if ( rv )
  {
    object = reinterpret_cast<TObject*> (value);
  }
  
  return rv;
}

//_____________________________________________________________________________
TObject* 
AliMpExMapIterator::Next()
{
/// Return the next object in iteration.
/// The returned object must not be deleted by the user.  

#if ROOT_VERSION_CODE >= 334081
//#if ROOT_VERSION_CODE >= 333824  // needed with Root v5.24.00-patches 
  Long64_t dummy;
#else
  Long_t dummy;
#endif
  TObject* o(0x0);
  Next(dummy,o);
  return o;
}

//_____________________________________________________________________________
TObject* 
AliMpExMapIterator::Next(Int_t& key)
{
/// Return the next object in iteration and fill the key.
/// The returned object must not be deleted by the user.  

  TObject* o;
#if ROOT_VERSION_CODE >= 334081
//#if ROOT_VERSION_CODE >= 333824  // needed with Root v5.24.00-patches 
  Long64_t index;
#else
  Long_t index;
#endif
  Next(index,o);
  key = (Int_t)(index);
  return o;
}

//_____________________________________________________________________________
TObject* 
AliMpExMapIterator::Next(Int_t& keyFirst, Int_t& keySecond)
{
/// Return the next object in iteration and fill the key.
/// The returned object must not be deleted by the user.  

#if ROOT_VERSION_CODE >= 334081
//#if ROOT_VERSION_CODE >= 333824  // needed with Root v5.24.00-patches 
  Long64_t index;
#else
  Long_t index;
#endif
  TObject* o(0x0);
  Next(index,o);
  keyFirst = AliMpExMap::GetPairFirst(index);
  keySecond = AliMpExMap::GetPairSecond(index);
  return o;
}

//_____________________________________________________________________________
TObject* 
AliMpExMapIterator::Next(TString& key)
{
/// Return the next object in iteration and fill the key.
/// The returned object must not be deleted by the user.  

#if ROOT_VERSION_CODE >= 334081
//#if ROOT_VERSION_CODE >= 333824  // needed with Root v5.24.00-patches 
  Long64_t index;
#else
  Long_t index;
#endif
  TObject* o(0x0);
  Next(index,o);
  key = AliMpExMap::GetString(index);
  return o;
}

//_____________________________________________________________________________
void 
AliMpExMapIterator::Reset()
{
/// Reset

  fIterator->Reset();
}

//_____________________________________________________________________________
const TCollection* 
AliMpExMapIterator::GetCollection() const 
{
/// Nothing to be returned here, AliMpExMap is not a TCollection

  return 0x0;
}
