// $Id$
// Category: basic
//
// Class AliMpPadPair
// ------------------
// Wrap up for std::pair<AliMpPad, AliMpPad>
// to avoid problems with CINT.
//
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include "AliMpPadPair.h"

ClassImp(AliMpPadPair)


//_____________________________________________________________________________
AliMpPadPair::AliMpPadPair(const AliMpPad& pad1, const AliMpPad& pad2)
  : TObject(),
    fPair(pad1, pad2) {
//
}

//_____________________________________________________________________________
AliMpPadPair::AliMpPadPair(const AliMpPadPair& right)
  : TObject(),
    fPair(right.GetFirst(), right.GetSecond()) {
//
}

//_____________________________________________________________________________
AliMpPadPair::AliMpPadPair()
  : TObject(),
    fPair(AliMpPad::Invalid(), AliMpPad::Invalid()) {
//
}

//_____________________________________________________________________________
AliMpPadPair::~AliMpPadPair() {
//
}

//_____________________________________________________________________________
Bool_t AliMpPadPair::operator == (const AliMpPadPair& right) const
{
  return fPair == right.fPair;
}

//_____________________________________________________________________________
Bool_t AliMpPadPair::operator!= (const AliMpPadPair& right) const
{
  return !(*this == right);
}

//_____________________________________________________________________________
AliMpPadPair& AliMpPadPair::operator = (const AliMpPadPair& right) 
{
  // check assignement to self
  if (this == &right) return *this;

  // base class assignement
  TObject::operator=(right);

  // assignement operator
  fPair = right.fPair;
  
  return *this;
}


