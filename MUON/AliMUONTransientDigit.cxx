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
$Log $
*/

#include "AliMUONTransientDigit.h"
#include <TObjArray.h>

ClassImp(AliMUONTransientDigit)

 
//____________________________________________________________________________

AliMUONTransientDigit::AliMUONTransientDigit(const AliMUONTransientDigit& digit)
{
// dummy copy constructor
}


AliMUONTransientDigit::AliMUONTransientDigit(Int_t ich, Int_t *digits): 
    AliMUONDigit(digits)
{
    //
    // Creates a MUON digit list object
    //

    fChamber     = ich;
    fTrackList   = new TObjArray;
    
}

AliMUONTransientDigit::~AliMUONTransientDigit() 
{
delete fTrackList;
}

AliMUONTransientDigit& AliMUONTransientDigit::operator =(const AliMUONTransientDigit& rhs)
{
// Dummy assignment operator
    return *this;
}


