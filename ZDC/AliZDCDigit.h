#ifndef ALIZDCDIGIT_H
#define ALIZDCDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//
//
//   ZDC digit = ADC Channels for each PM 
//
//_________________________________________________________________________

#include "AliDigitNew.h"

class AliZDCDigit : public AliDigitNew {

 public:
  
  AliZDCDigit() ;
  AliZDCDigit(Int_t *Sector, Int_t ADCValue);
  AliZDCDigit(const AliZDCDigit & digit);
  virtual ~AliZDCDigit() {}

  // Getters 
  virtual Float_t   GetSector(Int_t i) {return fSector[i];}
  virtual Float_t   GetADCValue() {return fADCValue;}

  // Operators
  Int_t operator == (AliZDCDigit &digit) {
    // Two digits are equal if they refers to the detector
    // in the same sub-volume (same procedure as for hits)
    Int_t i;
    for(i=0; i<2; i++) if(fSector[i]!=digit.GetSector(i)) return 0;
    return 1;
  }
  virtual AliZDCDigit& operator + (AliZDCDigit &digit) {
    // Adds the amplitude of digits 
    fADCValue += digit.fADCValue ;
    return *this ;
  }
  
 protected:

  //Data members
  Int_t   fSector[2];         // Detecor and tower in which light is produced
  Float_t fADCValue;          // ADC channel value

  // Print method
  virtual void Print(Option_t *) {
     printf(" -> DIGIT: Detector =  %d Quadrant =  %d ADCCh =  %f\n ",
     fSector[0], fSector[1], fADCValue);
  }
    
  ClassDef(AliZDCDigit,1)   // Digits in ZDC 

} ;

#endif //  ALIZDCDIGIT_H
