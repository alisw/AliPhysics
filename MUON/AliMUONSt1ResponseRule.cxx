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
// Class AliMUONSt1ResponseRule
// -----------------------------
// Describes a response rule.
// A "rule" is defined as being a set of electronic filters to be applied 
// (ie. a set of AliMUONSt1ResponseParameter) and a set of cathode pads to 
// which these filters should be applied (set of AliMUONSt1ElectronicElement)
// Included in AliRoot 2003/01/28

#include "AliMpPad.h"

#include "AliMUONSt1ResponseRule.h"
#include "AliMUONSt1ElectronicElement.h"
#include "AliMUONSt1ResponseParameter.h"

ClassImp(AliMUONSt1ResponseRule);

//__________________________________________________________________________
AliMUONSt1ResponseRule::AliMUONSt1ResponseRule()
  : TObject(),
    fElementList(),
    fParameters()
{
// default constructor
}

//__________________________________________________________________________
AliMUONSt1ResponseRule::~AliMUONSt1ResponseRule()
{
// destructor
}

//__________________________________________________________________________
void AliMUONSt1ResponseRule::AddElement(AliMUONSt1ElectronicElement* element)
{
// Add an electronic element to the list
// ---

  fElementList.Add(element);
}

//__________________________________________________________________________
void AliMUONSt1ResponseRule::AddParameter(AliMUONSt1ResponseParameter* param)
{
// Add an electronics parameter for this rule
// ---

  fParameters.Add(param);
}

//__________________________________________________________________________
Bool_t AliMUONSt1ResponseRule::Contains(const AliMpPad& pad) const
{
// Is this pad is contained in this rule's list
// ---

  TIter next(&fElementList);
  AliMUONSt1ElectronicElement* el;
  while ((el = static_cast<AliMUONSt1ElectronicElement*>(next()))){
    if (el->Contains(pad)) return kTRUE;
  }
  return kFALSE;
}
