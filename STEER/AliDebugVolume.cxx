/*
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

/*
$Log$
*/

#include "AliDebugVolume.h"

ClassImp(AliDebugVolume)



AliDebugVolume::AliDebugVolume()
{
  //
  // Default constructor
  //
}


AliDebugVolume::AliDebugVolume(const char *name,
		 Int_t copy, Float_t step, Float_t x, Float_t y, Float_t z, Int_t status)
    : TNamed(name, "Debug Volume")
{
//
// Constructor
//
    fCopy = copy;
    fX    = x;
    fY    = y;
    fZ    = z;
    fStep = step;
    fStatus = status;
}



Bool_t  AliDebugVolume::IsEqual(const char* name, const Int_t copy)
{
    return (copy == fCopy && strcmp(name, fName) == 0);
}

char*   AliDebugVolume::Status() const
{
    char* tmp;
    tmp = "Undefined";
    if (fStatus == 1) tmp = "Entering";
    if (fStatus == 2) tmp = "Exiting";   
    return tmp;
}


void AliDebugVolume::Copy(AliDebugVolume &volume) const
{
  //
  // Copy *this onto debug volume -- not implemented
  //
  Fatal("Copy","Not implemented!\n");
}













