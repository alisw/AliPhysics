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

//-----------------------------------------------------------------------
//    Class to debug entry and exit from a volume
//    Used by AliLego class
//    Author: A.Morsch
//-----------------------------------------------------------------------

#include "AliDebugVolume.h"

ClassImp(AliDebugVolume)


//_______________________________________________________________________
AliDebugVolume::AliDebugVolume():
  fCopy(0),
  fStep(0),
  fX(0),
  fY(0),
  fZ(0),
  fStatus(0)
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________
AliDebugVolume::AliDebugVolume(const char *name, Int_t copy, Float_t step, 
                               Float_t x, Float_t y, Float_t z, Int_t status):
  TNamed(name, "Debug Volume"),
  fCopy(copy),
  fStep(step),
  fX(x),
  fY(y),
  fZ(z),
  fStatus(status)
{
  //
  // Normal constructor
  //
}


//_______________________________________________________________________
Bool_t  AliDebugVolume::IsVEqual(const char* name, Int_t copy) const
{
  return (copy == fCopy && strcmp(name, fName) == 0);
}

//_______________________________________________________________________
const char*   AliDebugVolume::Status() const
{
  // Returns the status of the particle with respect to 
  // the current volume (Undefined, Entering, Exiting) 
  const char* tmp;
  tmp = "Undefined";
  if (fStatus == 1) tmp = "Entering";
  if (fStatus == 2) tmp = "Exiting";   
  return tmp;
}
