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

/* $Id$ */

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1SpecialMotif
// ----------------------------
// Encapsulate the distance between the center of a given daughter card
// and the pad/kapton connector.
// Included in AliRoot 2003/01/28

#include "AliMUONSt1SpecialMotif.h"


//__________________________________________________________________________
AliMUONSt1SpecialMotif::AliMUONSt1SpecialMotif(const TVector2& delta, 
                                               Double_t rotAngle)
  :fDelta(delta),
   fRotAngle(rotAngle)
{
// normal constructor
}

//__________________________________________________________________________
AliMUONSt1SpecialMotif::AliMUONSt1SpecialMotif()
  :fDelta(TVector2(0.,0.)),
   fRotAngle(0.)
{
// default constructor
}

//__________________________________________________________________________
AliMUONSt1SpecialMotif::AliMUONSt1SpecialMotif(const AliMUONSt1SpecialMotif& src)
  :fDelta(src.fDelta),
   fRotAngle(src.fRotAngle)
  
{
// copy constructor
}

//__________________________________________________________________________
AliMUONSt1SpecialMotif::~AliMUONSt1SpecialMotif()
{
// destructor
}

