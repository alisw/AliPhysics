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
// $MpId: AliMpPadPair.cxx,v 1.6 2005/08/26 15:43:36 ivana Exp $
// Category: basic
//
// Class AliMpPadPair
// ------------------
// Wrap up for std::pair<AliMpPad, AliMpPad>
// to avoid problems with CINT.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpPadPair.h"

ClassImp(AliMpPadPair)


//_____________________________________________________________________________
AliMpPadPair::AliMpPadPair(const AliMpPad& pad1, const AliMpPad& pad2)
  : TObject(),
    fPadFirst(pad1),
    fPadSecond(pad2) 
{
/// Standard constructor
}

//_____________________________________________________________________________
AliMpPadPair::AliMpPadPair(const AliMpPadPair& right)
  : TObject(),
    fPadFirst(right.GetFirst()),
    fPadSecond(right.GetSecond()) 
{
/// Copy constructor
}

//_____________________________________________________________________________
AliMpPadPair::AliMpPadPair()
  : TObject(),
    fPadFirst(AliMpPad::Invalid()),
    fPadSecond(AliMpPad::Invalid()) 
{
/// Default constructor
}

//_____________________________________________________________________________
AliMpPadPair::~AliMpPadPair() 
{
/// Destructor
}

//_____________________________________________________________________________
Bool_t AliMpPadPair::operator == (const AliMpPadPair& right) const
{
/// Equality operator 

  return (fPadFirst == right.fPadFirst && fPadSecond == right.fPadSecond);
}

//_____________________________________________________________________________
Bool_t AliMpPadPair::operator!= (const AliMpPadPair& right) const
{
/// Non-equality operator 

  return !(*this == right);
}

//_____________________________________________________________________________
AliMpPadPair& AliMpPadPair::operator = (const AliMpPadPair& right) 
{
/// Assignment operator 

  // check assignment to self
  if (this == &right) return *this;

  // base class assignment
  TObject::operator=(right);

  // assignment operator
  fPadFirst = right.fPadFirst;
  fPadSecond = right.fPadSecond;
  
  return *this;
}


