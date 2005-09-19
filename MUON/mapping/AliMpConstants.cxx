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
// $MpId: AliMpConstants.cxx,v 1.6 2005/08/26 15:43:36 ivana Exp $
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
  : TObject() 
{
/// Default constructor  
}

//_____________________________________________________________________________
AliMpConstants::~AliMpConstants() 
{
///Destructor
}

//_____________________________________________________________________________
Bool_t  AliMpConstants::IsEqual(Double_t length1, Double_t length2)
{
/// Compare lengths within the length tolerance.

  return TMath::Abs(length1 - length2) < fgkLengthTolerance;
}  


//_____________________________________________________________________________
Bool_t  AliMpConstants::IsEqual(const TVector2& v1, const TVector2& v2)
{
/// Compare x, y vector coordinates within the length tolerance.

  return (  TMath::Abs(v1.X() - v2.X()) 
          + TMath::Abs(v1.Y() - v2.Y())) < 2.*fgkLengthTolerance;
}
