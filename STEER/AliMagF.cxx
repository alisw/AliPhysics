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

//----------------------------------------------------------------------
// Basic magnetic field class
// Used in all the detectors, and also in the traking classes
// Author:
//----------------------------------------------------------------------

#include "AliLog.h"
#include "AliMagF.h"
Bool_t AliMagF::fgReadField = kTRUE;
ClassImp(AliMagF)

//_______________________________________________________________________
AliMagF::AliMagF():
  fMap(0),
  fType(0),
  fInteg(0),
  fFactor(0),
  fMax(0)
{
  //
  // Default constructor
  //
}

//_______________________________________________________________________
AliMagF::AliMagF(const char *name, const char *title, Int_t integ, 
                 Float_t factor, Float_t fmax):
  TNamed(name,title),
  fMap(0),
  fType(0),
  fInteg(0),
  fFactor(factor),
  fMax(fmax)
{
  //
  // Standard constructor
  //
    if(integ<0 || integ > 2) {
      AliWarning(Form(
              "Invalid magnetic field flag: %5d; Helix tracking chosen instead"
              ,integ));
      fInteg = 2;
    } else {
      fInteg = integ;
    }
    fType = kUndef;
    //
}

//_______________________________________________________________________
void AliMagF::Field(Float_t*, Float_t *b) const
{
  //
  // Method to return the field in one point -- dummy in this case
  //
  AliWarning("Undefined MagF Field called, returning 0");
  b[0]=b[1]=b[2]=0;
}
