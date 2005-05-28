#ifndef ALIPHOSDIGIT_H
#define ALIPHOSDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

//_________________________________________________________________________
//  PHOS digit: Id
//              energy
//              3 identifiers for the primary particle(s) at the origine of the digit
//  The digits are made in FinishEvent() by summing all the hits in a single PHOS crystal or PPSD gas cell
//  It would be nice to replace the 3 identifiers by an array, but, because digits are kept in a TClonesArray,
//   it is not possible to stream such an array... (beyond my understqnding!)
//
//*-- Author: Laurent Aphecetche & Yves Schutz (SUBATECH)

// --- ROOT system ---

//#include "TObject.h" 

// --- Standard library ---

// --- AliRoot header files ---

#include "AliDigitNew.h"

class AliPHOSDigit : public AliDigitNew {

  friend ostream& operator << ( ostream& , const AliPHOSDigit&) ;

 public:
  
  AliPHOSDigit() ;
  AliPHOSDigit(Int_t primary, Int_t id, Int_t DigEnergy, Float_t Time, Int_t index = -1) ;
  AliPHOSDigit(const AliPHOSDigit & digit) ;
  virtual ~AliPHOSDigit() ;

  Bool_t operator==(const AliPHOSDigit &rValue) const;
  AliPHOSDigit& operator+(AliPHOSDigit const &rValue) ;
  AliPHOSDigit& operator*(Float_t factor) ; 
  //  AliPHOSDigit& operator=(const AliPHOSDigit) {
  //  Fatal("operator = ", "not implemented") ; return *this ; } 
  Int_t   Compare(const TObject * obj) const ;  
  Int_t   GetNprimary() const { return fNprimary ; }
  Int_t   GetPrimary(Int_t index) const ; 
  Float_t GetTime(void) const {return fTime ;}
  Float_t GetTimeR(void) const {return fTimeR ;}
  Bool_t  IsSortable() const { return kTRUE ; }
  void    Print(const Option_t * = "") const;
  void    SetAmp(Int_t Amp) { fAmp=Amp ; } 
  void    SetTime(Float_t time) {fTime = time ;}
  void    SetTimeR(Float_t time) {fTimeR = time ;}
  void    ShiftPrimary(Int_t shift); // shift to separate different TreeK in merging

 private:

  Int_t fNprimary ;        // Number of primaries
  Int_t * fPrimary ;       //[fNprimary] Array of primaries      
  Float_t fTime ;          // Calculcated time 
  Float_t fTimeR ;         // Earliest time: to be used by Digits2Raw
    
  ClassDef(AliPHOSDigit,3)   // Digit in PHOS 

} ;

#endif //  ALIPHOSDIGIT_H
