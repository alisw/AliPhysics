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

/*
  $Log$
  Revision 1.1  2000/06/12 15:35:17  jbarbosa
  Cleaned up version.

*/


#include "AliRICHTransientDigit.h"
#include <TObjArray.h>

ClassImp(AliRICHTransientDigit)
    
//____________________________________________________________________________
AliRICHTransientDigit::AliRICHTransientDigit(Int_t ich, Int_t *digits): 
    AliRICHDigit(digits)
{
    //
    // Creates a RICH digit list object
    //
    
    fChamber = ich;
    fTrackList   = new TObjArray;
    
}
//_____________________________________________________________________________

AliRICHTransientDigit::~AliRICHTransientDigit()
{
  delete fTrackList;
}
