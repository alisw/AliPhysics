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
  AliPHOSDigit::AliPHOSDigit() : fPrimary(0)
{
}

//____________________________________________________________________________
AliPHOSDigit::AliPHOSDigit(Int_t primary, Int_t id, Int_t DigEnergy) 
{  
  fId         = id ;
  fAmp        = DigEnergy ;
  fPrimary    = new Int_t[1] ; 
  fPrimary[0] = primary ;
  fNprimary   = 1 ; 
}

//____________________________________________________________________________
AliPHOSDigit::AliPHOSDigit(const AliPHOSDigit & digit) 
{  
  fId       = digit.fId;
  fAmp      = digit.fAmp ;
  fNprimary = digit.GetNprimary() ;
  fPrimary  = new Int_t[fNprimary] ;
  Int_t * primary = digit.GetPrimary() ;
  Int_t  index ;
  for ( index = 0 ; index < fNprimary ; index++ )
    fPrimary[index] = primary[index] ;
}

//____________________________________________________________________________
AliPHOSDigit::~AliPHOSDigit()
{
  delete fPrimary ; 
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
  
  Int_t * tempo = new Int_t[fNprimary] ; 
  Int_t index ; 
  
  Int_t oldfNprimary = fNprimary ; 

  for ( index = 0 ; index < oldfNprimary ; index++ ){
    tempo[index] = fPrimary[index] ; 
  }  
 
  delete fPrimary ; 
  fNprimary += digit.GetNprimary() ; 
  fPrimary = new Int_t[fNprimary] ; 
  
  for ( index = 0 ; index < oldfNprimary  ; index++ ) { 
    fPrimary[index] = tempo[index] ; 
  }

  Int_t jndex = 0 ; 
  for ( index = oldfNprimary ; index < fNprimary ; index++ ) { 
    fPrimary[index] = digit.fPrimary[jndex] ; 
    jndex++ ; 
  }
      
  return *this ;
}

//____________________________________________________________________________
ostream& operator << ( ostream& out , const AliPHOSDigit & digit)
{
  out << "ID " << digit.fId << " Energy = " << digit.fAmp ;

  return out ;
}

//______________________________________________________________________________
void AliPHOSDigit::Streamer(TBuffer &R__b)
{
  assert(0==1) ; 
   // Stream an object of class AliPHOSDigit.
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliDigitNew::Streamer(R__b);
      R__b.ReadArray(fPrimary);
      R__b >> fNprimary;
   } else {
      R__b.WriteVersion(AliPHOSDigit::IsA());
      AliDigitNew::Streamer(R__b);
      R__b.WriteArray(fPrimary, fNprimary);
      R__b << fNprimary;

   }
}

   
