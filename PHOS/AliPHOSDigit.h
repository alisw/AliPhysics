#ifndef ALIPHOSDIGIT_H
#define ALIPHOSDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  The digit  class: a list of abs Id, energy//
//  Version SUBATECH                          //
//  Author Laurent Aphecetche     SUBATECH    //
//      comment: added sortable YS            //  
//                                            //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TObject.h" 

// --- Standard library ---

// --- AliRoot header files ---

#include "AliDigitNew.h"

class AliPHOSDigit : public AliDigitNew {
  
public:
 
  AliPHOSDigit() ;
  AliPHOSDigit(Int_t primary, Int_t id, Int_t DigEnergy) ;
  AliPHOSDigit(const AliPHOSDigit & digit) ;
  virtual ~AliPHOSDigit() ;  

  Bool_t operator==(AliPHOSDigit const &rValue) const;
  AliPHOSDigit& operator+(AliPHOSDigit const &rValue) ;
  
  friend ostream& operator << ( ostream& , const AliPHOSDigit&) ;
  
  Int_t   Compare(TObject * obj) ;  
  Int_t   GetNprimary() const { return fNprimary ; }
  Int_t   GetPrimary(Int_t index) const ; 
  Bool_t  IsSortable() const { return kTRUE ; }
  void    SetAmp(Int_t Amp) { fAmp=Amp ; } 

private:

  Int_t fPrimary1 ;          // first primary (because I do not know how to stream *fPrimary) 
  Int_t fPrimary2 ;          // second primary (because I do not know how to stream *fPrimary) 
  Int_t fPrimary3 ;          // third primary (because I do not know how to stream *fPrimary) 
  Int_t fNprimary ;          // Number of primaries
  
  ClassDef(AliPHOSDigit,1)   // Digit in PHOS, version 1 

} ;

#endif //  ALIPHOSDIGIT_H
