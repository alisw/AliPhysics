#ifndef ALIPHOSDIGIT_H
#define ALIPHOSDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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

#include "TObject.h" 

// --- Standard library ---

// --- AliRoot header files ---

#include "AliDigitNew.h"

class AliPHOSDigit : public AliDigitNew {

  friend class ostream& operator << ( ostream& , const AliPHOSDigit&) ;

 public:
  
  AliPHOSDigit() ;
  AliPHOSDigit(Int_t primary, Int_t id, Int_t DigEnergy, Int_t index = -1) ;
  AliPHOSDigit(const AliPHOSDigit & digit) ;
  virtual ~AliPHOSDigit(){
    // dtor 
  } 

  Bool_t operator==(const AliPHOSDigit &rValue) const;
  AliPHOSDigit& operator+(AliPHOSDigit const &rValue) ;
    
  Int_t   Compare(TObject * obj) ;  
  Int_t   GetNprimary() const { 
    // returns the number of primaries
    return fNprimary ; }
  Int_t   GetPrimary(Int_t index) const ; 
  Bool_t  IsSortable() const { 
    // says that AliPHOSDigits are sortable (needed for Sort method
    return kTRUE ; }
  void    SetAmp(Int_t Amp) { 
    // sets the amplitude data member 
    fAmp=Amp ; } 

 private:

  Int_t fNprimary ;     // Number of primaries
  Int_t fNMaxPrimary ;  //! Max Number of primaries
  Int_t * fPrimary ;    //[fNprimary]  Array of primaries       
    
  ClassDef(AliPHOSDigit,1)   // Digit in PHOS 

} ;

#endif //  ALIPHOSDIGIT_H
