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

//_________________________________________________________________________
//  PHOS digit: Id
//              energy
//              3 identifiers for the primary particle(s) at the origine of the digit
//  The digits are made in FinishEvent() by summing all the hits in a single PHOS crystal or PPSD gas cell
//  It would be nice to replace the 3 identifiers by an array, but, because digits are kept in a TClonesQArray,
//   it is not possible to stream such an array... (beyond my understqnding!)
//
//*-- Author: Laurent Aphecetche & Yves Schutz (SUBATECH)


// --- ROOT system ---

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSDigit.h"


ClassImp(AliPHOSDigit)

//____________________________________________________________________________
  AliPHOSDigit::AliPHOSDigit() 
{
  // default ctor 

  fIndexInList = -1 ; 
  fNprimary    = 0 ;  
  fNMaxPrimary = 5 ; 
}

//____________________________________________________________________________
AliPHOSDigit::AliPHOSDigit(Int_t primary, Int_t id, Int_t DigEnergy, Int_t index) 
{  
  // ctor with all data 

  fNMaxPrimary = 5 ; 
  fAmp         = DigEnergy ;
  fId          = id ;
  fIndexInList = index ; 
  if( primary != -1){
    fNprimary    = 1 ; 
    fPrimary[0]  = primary ;
  }
  else{  //If the contribution of this primary smaller than fDigitThreshold (AliPHOSv1)
    fNprimary    = 0 ; 
    fPrimary[0]  = -1 ;
  }
  Int_t i ;
  for ( i = 1; i < fNMaxPrimary ; i++)
    fPrimary[i]  = -1 ; 
}

//____________________________________________________________________________
AliPHOSDigit::AliPHOSDigit(const AliPHOSDigit & digit) 
{
  // copy ctor
  

  fNMaxPrimary = digit.fNMaxPrimary ;  
  Int_t i ;
  for ( i = 0; i < fNMaxPrimary ; i++)
    fPrimary[i]  = digit.fPrimary[i] ;
  fAmp         = digit.fAmp ;
  fId          = digit.fId;
  fIndexInList = digit.fIndexInList ; 
  fNprimary    = digit.fNprimary ;
}

//____________________________________________________________________________
AliPHOSDigit::~AliPHOSDigit() 
{
  // Delete array of primiries if any
  
}

//____________________________________________________________________________
Int_t AliPHOSDigit::Compare(TObject * obj)
{
  // Compares two digits with respect to its Id
  // to sort according increasing Id

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
  // retrieves the primary particle number given its index in the list 
  Int_t rv = -1 ;
  if ( index <= fNprimary ){
    rv = fPrimary[index-1] ;
  } 

  return rv ; 
  
}

//____________________________________________________________________________
Bool_t AliPHOSDigit::operator==(AliPHOSDigit const & digit) const 
{
  // Two digits are equal if they have the same Id
  
  if ( fId == digit.fId ) 
    return kTRUE ;
  else 
    return kFALSE ;
}
 
//____________________________________________________________________________
AliPHOSDigit& AliPHOSDigit::operator+(AliPHOSDigit const & digit) 
{
  // Adds the amplitude of digits and completes the list of primary particles
  // if amplitude is larger than 
  
  fAmp += digit.fAmp ;
  
  Int_t max1 = fNprimary ; 
  Int_t max2 = digit.fNprimary ; 
  
  if ( fNprimary + digit.fNprimary < fNMaxPrimary ) {
    fNprimary += digit.fNprimary ; 
  
    Int_t index ; 
    for (index = 0 ; index < max2 ; index++)
      fPrimary[index+max1] = (digit.fPrimary)[index] ; 
  }
  else{
    cout << "AliPHOSDigit >> Increase NMaxPrimary"<< endl ;
    Int_t index ; 
    for (index = max1; index < fNMaxPrimary ; index++)
      fPrimary[index] = (digit.fPrimary)[index-max1] ; 
    fNprimary = fNMaxPrimary ;
  }  
  return *this ;
}

//____________________________________________________________________________
ostream& operator << ( ostream& out , const AliPHOSDigit & digit)
{
  // Prints the data of the digit
  
  out << "ID " << digit.fId << " Energy = " << digit.fAmp << endl ; 
  Int_t i ;
  for(i=0;i<digit.fNprimary;i++)
    out << "Primary " << i+1 << " = " << digit.fPrimary[i] << endl ;
  out << "Position in list = " << digit.fIndexInList << endl ; 
  return out ;
}


