#ifndef ALIEMCALTICK_H
#define ALIEMCALTICK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Class for EMCAL time digitization
//  holds info on a single
//  digit time bucket
//                  
//*-- Author: Dmitri Peressounko (SUBATECH)


// --- ROOT system ---
#include "TObject.h"
// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALTick: public TObject {

public:
  AliEMCALTick() ;          
  AliEMCALTick(Float_t time, Float_t a, Float_t slope) ;
  virtual ~AliEMCALTick(){}

  Int_t   Compare(const TObject * obj) const ;  
  Bool_t  IsSortable() const { return kTRUE ; }

  Float_t CrossingTime(Float_t threshold) const
    //Calculates time, when rizing front of the signal crosses 
    {if(fB) return fTime + (threshold - fA)/fB ;
    else return 1. ;} //return very big time

  Float_t GetTime(void) const {return fTime ;}

  void operator+=(AliEMCALTick const &rValue) ;


private:
  Float_t fTime ;     //!time of the beginning of this tick
  Float_t fA ;        //!constant
  Float_t fB ;        //!slope        

  ClassDef(AliEMCALTick,1)  // description 

};

#endif // AliEMCALTICK_H
