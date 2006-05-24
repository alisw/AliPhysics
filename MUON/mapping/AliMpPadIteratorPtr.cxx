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
// $MpId: AliMpPadIteratorPtr.cxx,v 1.7 2006/05/24 13:58:29 ivana Exp $
// Category: basic
//
// Class AliMpPadIteratorPtr
// --------------------------
// Pointer to the virtual pad iterator;
// enables to allocate the virtual pad iterator on stack.
// Usage:
// MVIndexed* myIndexed = MyIndexed()
// MVIterator& it = *AliMpPadIteratorPtr(myIndexed->CreateIterator());
//
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpPadIteratorPtr.h"
#include "AliMpVPadIterator.h"

/// \cond CLASSIMP
ClassImp(AliMpPadIteratorPtr)
/// \endcond

//_____________________________________________________________________________
AliMpPadIteratorPtr::AliMpPadIteratorPtr(AliMpVPadIterator* it)
  : fIterator(it)
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMpPadIteratorPtr::AliMpPadIteratorPtr(const AliMpPadIteratorPtr& right) 
  : TObject(right) 
{
/// Protected copy constructor (not provided) 

  Fatal("AliMpPadIteratorPtr", "Copy constructor not provided.");
}

//_____________________________________________________________________________
AliMpPadIteratorPtr::~AliMpPadIteratorPtr() 
{
/// Destructor

  delete fIterator;
}

// operators

//_____________________________________________________________________________
AliMpPadIteratorPtr& 
AliMpPadIteratorPtr::operator=(const AliMpPadIteratorPtr& right)
{
/// Protected assignment operator (not provided) 

  // check assignment to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignment operator not provided.");
    
  return *this;  
}    

