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

//_________________________________________________________________________
// Digit class for PHOS that contains an absolute ID and an energy
//*-- Author : Laurent Aphecetche  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

#include <iostream>

// --- AliRoot header files ---

#include "AliPHOSDigit.h"


ClassImp(AliPHOSDigit)

//____________________________________________________________________________
AliPHOSDigit::AliPHOSDigit(Int_t id, Int_t DigEnergy) : AliDigitNew( ), fId(id),fAmp(DigEnergy)
{  
  // This part should be a true Digitization, but is not for the moment.
  //  fAmp = energy ;
  fId = id;
  fAmp = DigEnergy ;
}

//____________________________________________________________________________
Bool_t AliPHOSDigit::operator==(AliPHOSDigit const &Digit) const 
{
  if ( fId == Digit.fId ) 
    return kTRUE ;
  else 
    return kFALSE ;
}

//____________________________________________________________________________
AliPHOSDigit& AliPHOSDigit::operator+(AliPHOSDigit const &Digit) 
{
  fAmp += Digit.fAmp ;

  return *this ;
}

//____________________________________________________________________________
ostream& operator << ( ostream& out , const AliPHOSDigit& Digit)
{
  out << "ID " << Digit.fId << " Energy = " << Digit.fAmp ;

  return out ;
}

//____________________________________________________________________________
Int_t AliPHOSDigit::Compare(TObject * obj)
{
  Int_t rv ; 

  AliPHOSDigit * digit = (AliPHOSDigit *)obj ; 

  Int_t iddiff = fId - digit->GetId() ; 

  if ( iddiff > 0 ) 
    rv = 1 ;
  else if ( iddiff < 0 )
    rv = -1 ; 
  else
    rv = 0 ;
  
  return rv ; 

}
