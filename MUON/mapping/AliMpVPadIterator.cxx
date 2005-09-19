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
// $MpId: AliMpVPadIterator.cxx,v 1.5 2005/08/26 15:43:36 ivana Exp $
// Category: basic
//
// Class AliMpVPadIterator
// -----------------------
// Abstract base class, which defines an iterator over pads
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpVPadIterator.h"

ClassImp(AliMpVPadIterator)

//___________________________________________________________________
AliMpVPadIterator::AliMpVPadIterator():
    TObject()
{
/// Default constructor
}

//___________________________________________________________________
AliMpVPadIterator::AliMpVPadIterator(const AliMpVPadIterator& right)
  : TObject(right)
{
/// Copy constructor
}

//___________________________________________________________________
AliMpVPadIterator::~AliMpVPadIterator()
{
/// Destructor
}

//
// operators
//

//___________________________________________________________________
AliMpVPadIterator& 
AliMpVPadIterator::operator = (const AliMpVPadIterator& right)
{
/// Assignment operator

  // check assignment to self
  if (this == &right) return *this;

  // base class assignment
  TObject::operator=(right);

  return *this;
}  

