////////////////////////////////////////////////
//  Digitization class for set: TOF           //
//  AliTOFRawDigit  class                     //
//  Interface                                 // 
//  Description                               //
//*-- Authors: Pierella, Seganti, Vicinanza   //
//    (Bologna and Salerno University)        //
////////////////////////////////////////////////


#ifndef ALITOFRAWDIGIT_H
#define ALITOFRAWDIGIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "TObject.h"
#include "TClonesArray.h"

//_______________________________________________________
class AliTOFRawDigit : public TObject{

public:
  AliTOFRawDigit(); 
  virtual ~AliTOFRawDigit(){};
  
protected:
  Int_t fTreeD;     // class under construction
  Int_t fRawDigits; // class under construction
    
    
  ClassDef(AliTOFRawDigit,2) // TOF Digit in rawdata format
};

#endif /* ALITOFRAWDIGIT_H */
