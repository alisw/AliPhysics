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
#include <cassert> 

// --- AliRoot header files ---

#include "AliPHOSDigit.h"


ClassImp(AliPHOSDigit)

//____________________________________________________________________________
  AliPHOSDigit::AliPHOSDigit() 
{
  fNprimary = 0 ;  
  fPrimary1 = -1 ; 
  fPrimary2 = -1 ; 
  fPrimary3 = -1 ;
}

//____________________________________________________________________________
AliPHOSDigit::AliPHOSDigit(Int_t primary, Int_t id, Int_t DigEnergy) 
{  
  fId         = id ;
  fAmp        = DigEnergy ;
  fPrimary1   = primary ;
  fNprimary   = 1 ; 
}

//____________________________________________________________________________
AliPHOSDigit::AliPHOSDigit(const AliPHOSDigit & digit) 
{  
  fId       = digit.fId;
  fAmp      = digit.fAmp ;
  fNprimary = digit.fNprimary ;
  fPrimary1 = digit.fPrimary1 ;
  fPrimary2 = digit.fPrimary2 ;
  fPrimary3 = digit.fPrimary3 ;
}

//____________________________________________________________________________
AliPHOSDigit::~AliPHOSDigit()
{
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

//____________________________________________________________________________
Int_t AliPHOSDigit::GetPrimary(Int_t index) const
{
  Int_t rv = -1 ; 
  if ( index > 3 )
    cout << "AliPHOSDigit  ERROR > only 3 primaries allowed" << endl ; 
  else {
    switch (index) {  
    case 1 :
      rv = fPrimary1 ; 
      break ; 
    case 2 :
      rv = fPrimary2 ; 
      break ; 
    case 3 :
      rv = fPrimary3 ; 
      break ; 
    }
  } 

  return rv ; 

}

//____________________________________________________________________________
Bool_t AliPHOSDigit::operator==(AliPHOSDigit const & digit) const 
{
  if ( fId == digit.fId ) 
    return kTRUE ;
  else 
    return kFALSE ;
}
 
//____________________________________________________________________________
AliPHOSDigit& AliPHOSDigit::operator+(AliPHOSDigit const & digit) 
{
  fAmp += digit.fAmp ;
 
  // Here comes something crummy ... but I do not know how to stream pointers
  // because AliPHOSDigit is in a TCLonesArray
  
  Int_t tempo1[3] ; 
  tempo1[0] = fPrimary1 ; 
  tempo1[1] = fPrimary2 ; 
  tempo1[2] = fPrimary3 ; 

  Int_t tempo2[3] ; 
  tempo2[0] = digit.fPrimary1 ; 
  tempo2[1] = digit.fPrimary2 ; 
  tempo2[2] = digit.fPrimary3 ; 

  Int_t max1 = fNprimary ; 
  Int_t max2 = digit.fNprimary ; 
 
  if ( fNprimary >= 3 ) {
    cout << "AliPHOSDigit + operator  ERROR > too many primaries, modify AliPHOSDigit" << endl ; 
  } 
  else {
    fNprimary += digit.fNprimary ; 
    if ( fNprimary > 3 ) {
      cout << "AliPHOSDigit + operator  ERROR > too many primaries, modify AliPHOSDigit" << endl ; 
      fNprimary = 3 ;
    }
    
    Int_t tempo3[3] ;
    Int_t index ; 
    for (index = 0 ; index < max1 ; index++)
      tempo3[index] = tempo1[index] ;
    for (index = 0 ; index < max2 ; index++)
      tempo3[index+max1] = tempo2[index] ; 
    
    fPrimary1 = tempo3[0] ; 
    fPrimary2 = tempo3[1] ; 
    fPrimary3 = tempo3[2] ; 

  }
 // end of crummy stuff      

  return *this ;
}

//____________________________________________________________________________
ostream& operator << ( ostream& out , const AliPHOSDigit & digit)
{
  out << "ID " << digit.fId << " Energy = " << digit.fAmp ;

  return out ;
}


