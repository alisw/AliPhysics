#ifndef ALIPHOSTICK_H
#define ALIPHOSTICK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Class for PHOS time digitization
//                  
//*-- Author: Dmitri Peressounko (SUBATECH)


// --- ROOT system ---
#include "TObject.h"
// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSTick: public TObject {

public:
  AliPHOSTick() ;          
  AliPHOSTick(Float_t time, Float_t a, Float_t slope) ;
  virtual ~AliPHOSTick(){}

  Int_t   Compare(const TObject * obj) const ;  
  Bool_t  IsSortable() const { return kTRUE ; }

  Float_t CrossingTime(Float_t threshold) const
    //Calculates time, when rizing front of the signal crosses 
    {if(fB) return fTime + (threshold - fA)/fB ;
    else return 1. ;} //return very big time

  Float_t GetTime(void){return fTime ;}

  void operator+=(AliPHOSTick const &rValue) ;


private:
  Float_t fTime ;     //!time of the beginning of this tick
  Float_t fA ;        //!constant
  Float_t fB ;        //!slope        

  ClassDef(AliPHOSTick,1)  // description 

};

#endif // AliPHOSTICK_H
