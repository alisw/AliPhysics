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
//
//*-- Author: Laurent Aphecetche & Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)


// --- ROOT system ---

// --- Standard library ---

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
AliPHOSDigit::AliPHOSDigit(Int_t primary, Int_t id, Int_t digEnergy, Float_t time, Int_t index) 
{  
  // ctor with all data 

  fNMaxPrimary = 5 ; 
  fAmp         = digEnergy ;
  fTime        = time ;
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
  fTime        = digit.fTime ;
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
Int_t AliPHOSDigit::Compare(const TObject * obj) const
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
  if ( index <= fNprimary && index > 0){
    rv = fPrimary[index-1] ;
  } 

  return rv ; 
  
}
//____________________________________________________________________________
void AliPHOSDigit::Print(Option_t *option) const
{
  printf("PHOS digit: Amp=%d, Id=%d\n",fAmp,fId);
}
//____________________________________________________________________________
void AliPHOSDigit::ShiftPrimary(Int_t shift)
{
  //shifts primary number to BIG offset, to separate primary in different TreeK
  Int_t index  ;
  for(index = 0; index <fNprimary; index ++ ){
    fPrimary[index] = fPrimary[index]+ shift ;
  } 
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
  
  Int_t toAdd = fNprimary ;
  if(digit.fNprimary>0){
    if(fAmp < digit.fAmp){//most energetic primary in second digit => first primaries in list from second digit
      for (Int_t index = 0 ; index < digit.fNprimary ; index++){
	for (Int_t old = 0 ; old < fNprimary ; old++) { //already have this primary?
	  if(fPrimary[old] == (digit.fPrimary)[index]){
	    fPrimary[old] = -1 ; //removed
	    toAdd-- ;
	    break ;
	  }
	}
      }
      Int_t nNewPrimaries = digit.fNprimary+toAdd ;
      if(nNewPrimaries >fNMaxPrimary)  //Do not change primary list
	Error("Operator +", "Increase NMaxPrimary") ;
      else{
	for(Int_t index = fNprimary-1 ; index >=0 ; index--){ //move old primaries	
	  if(fPrimary[index]>-1){
	    toAdd-- ;
	    fPrimary[fNprimary+toAdd]=fPrimary[index] ;
	  }
	}
	//copy new primaries
	for(Int_t index = 0; index < digit.fNprimary ; index++){
	  fPrimary[index] = (digit.fPrimary)[index] ;	
	}
	fNprimary = nNewPrimaries ;
      }
    }
    else{ //add new primaries to the end
      for(Int_t index = 0 ; index < digit.fNprimary ; index++){
	Bool_t deja = kTRUE ;
	for(Int_t old = 0 ; old < fNprimary; old++) { //already have this primary?
	  if(fPrimary[old] == (digit.fPrimary)[index]){
	    deja = kFALSE;
	    break ;
	  }      
	}
	if(deja){
	  fPrimary[fNprimary] = (digit.fPrimary)[index] ; 
	  fNprimary++ ;
	  if(fNprimary>fNMaxPrimary) {
	    Error("Operator +", "Increase NMaxPrimary") ;
	    break ;
	  }
	}
      }
      
    }
  }
  
  fAmp += digit.fAmp ;
  if(fTime > digit.fTime)
    fTime = digit.fTime ;
  
  return *this ;
}

//____________________________________________________________________________
ostream& operator << ( ostream& out , const AliPHOSDigit & digit)
{
  // Prints the data of the digit
  
//   out << "ID " << digit.fId << " Energy = " << digit.fAmp << " Time = " << digit.fTime << endl ; 
//   Int_t i ;
//   for(i=0;i<digit.fNprimary;i++)
//     out << "Primary " << i+1 << " = " << digit.fPrimary[i] << endl ;
//   out << "Position in list = " << digit.fIndexInList << endl ; 
  digit.Warning("operator <<", "Implement differently") ; 
  return out ;
}


