// $Id$
// Category: basic
//
// Class AliMpConstants
// --------------------
// Class for globally used constants definition.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay

#include <TMath.h>
#include <TVector2.h>

#include "AliMpConstants.h"

ClassImp(AliMpConstants)

// static data
const Double_t AliMpConstants::fgkLengthTolerance = 1e-04; // 0.1 mum;
const Double_t AliMpConstants::fgkLengthStep = 1.; // 1 mm;
const Int_t    AliMpConstants::fgkStartPadIndex = 1;

//_____________________________________________________________________________
AliMpConstants::AliMpConstants()
  : TObject() {
//  
}

//_____________________________________________________________________________
AliMpConstants::~AliMpConstants() {
//
}

//_____________________________________________________________________________
Bool_t  AliMpConstants::IsEqual(Double_t length1, Double_t length2)
{
// Compares lengths within the length tolerance.
// ---

  return TMath::Abs(length1 - length2) < fgkLengthTolerance;
}  


//_____________________________________________________________________________
Bool_t  AliMpConstants::IsEqual(const TVector2& v1, const TVector2& v2)
{
// Compares x, y vector coordinates within the length tolerance.
// ---

  return (  TMath::Abs(v1.X() - v2.X()) 
          + TMath::Abs(v1.Y() - v2.Y())) < 2.*fgkLengthTolerance;
}
