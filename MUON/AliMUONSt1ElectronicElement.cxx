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
// Class AliMUONSt1ElectronicElement
// ---------------------------------
// Describes a set of pads either by defining
// a range of indices, or
// a range of (centimeters) positions or
// a range of electronic channel numbers or
// a range of MANU numbers or, finally,
// a range of gassiplex/MANAS numbers, in a given range of MANU addresses
// Included in AliRoot 2003/01/28

#include "AliMUONSt1ElectronicElement.h"
#include "AliMpPad.h"

ClassImp(AliMUONSt1ElectronicElement)

//______________________________________________________________________________
AliMUONSt1ElectronicElement::AliMUONSt1ElectronicElement()
  :TObject(),fDescription(kNone)
{
// default constructor
}

//______________________________________________________________________________
AliMUONSt1ElectronicElement::AliMUONSt1ElectronicElement(TDescription descr)
  :TObject(),fDescription(descr)
{
// normal constructor
  switch (descr){
    case kXY : SetRange(0,0.,0.); SetRange(1,0.,0.); break;
    default:   SetRange(0,0,0); SetRange(1,0,0); break;
  }
}

//______________________________________________________________________________
AliMUONSt1ElectronicElement::~AliMUONSt1ElectronicElement() 
{
// destructor
}

//______________________________________________________________________________
void AliMUONSt1ElectronicElement::SetRange(Int_t numVar,Int_t i1,Int_t i2)
{
// set the range of the <numVar>th variables, in all cases but kXY
// ---

  if (fDescription==kXY) {
    fRanges[numVar][0].x = (Double_t)i1;
    fRanges[numVar][1].x = (Double_t)i2;
  } else {
    fRanges[numVar][0].i = i1;
    fRanges[numVar][1].i = i2;
  }
}

//______________________________________________________________________________
void AliMUONSt1ElectronicElement::SetRange(Int_t numVar,Double_t x1,Double_t x2)
{
// set the range of the <numVar>th variable, in cases kXY
// ---

  if (fDescription==kXY) {
    fRanges[numVar][0].x = x1;
    fRanges[numVar][1].x = x2;
  } else {
    fRanges[numVar][0].i = (Int_t)x1;
    fRanges[numVar][1].i = (Int_t)x2;
  }
}

//______________________________________________________________________________
Bool_t AliMUONSt1ElectronicElement::IsInRange(Int_t numVar,Int_t i) const
{
// is the given value in the <numVar>th variable
// ---

  return (fRanges[numVar][0].i<=i) && (fRanges[numVar][1].i>=i);
}

//______________________________________________________________________________
Bool_t AliMUONSt1ElectronicElement::IsInRange(Int_t numVar,Double_t x) const
{
// is the given value in the <numVar>th variable
// ---

  return (fRanges[numVar][0].x<=x) && (fRanges[numVar][1].x>=x);
}

//______________________________________________________________________________
Bool_t AliMUONSt1ElectronicElement::Contains(const AliMpPad& pad) const
{
// is the pad <pad> contained in this range
// ---

  switch(fDescription){
    case kNone:
      return kFALSE;
    case kIJ  :
      return (  IsInRange(0,pad.GetIndices().GetFirst())
             && IsInRange(1,pad.GetIndices().GetSecond())
	     );
    case kXY  :
      return (  IsInRange(0,pad.Position().X())
             && IsInRange(1,pad.Position().Y())
	     );
    case kMGC  :
      return (  IsInRange(0,pad.GetLocation().GetFirst())
             && IsInRange(1,pad.GetLocation().GetSecond())
	     );
    case kMG   :
      return (  IsInRange(0,pad.GetLocation().GetFirst())
             && IsInRange(1,pad.GetLocation().GetSecond() >> 4)
	     );
    case kM    :
      return (  IsInRange(0,pad.GetLocation().GetFirst()));
      
    default: return kFALSE;
  }
}
