// $Id$
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
    fPadSecond(pad2) {
//
}

//_____________________________________________________________________________
AliMpPadPair::AliMpPadPair(const AliMpPadPair& right)
  : TObject(),
    fPadFirst(right.GetFirst()),
    fPadSecond(right.GetSecond()) {
//
}

//_____________________________________________________________________________
AliMpPadPair::AliMpPadPair()
  : TObject(),
    fPadFirst(AliMpPad::Invalid()),
    fPadSecond(AliMpPad::Invalid()) {
//
}

//_____________________________________________________________________________
AliMpPadPair::~AliMpPadPair() {
//
}

//_____________________________________________________________________________
Bool_t AliMpPadPair::operator == (const AliMpPadPair& right) const
{
  return (fPadFirst == right.fPadFirst && fPadSecond == right.fPadSecond);
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
  fPadFirst = right.fPadFirst;
  fPadSecond = right.fPadSecond;
  
  return *this;
}


