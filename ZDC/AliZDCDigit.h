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
  AliZDCDigit(Int_t Det, Int_t Quad, Float_t ADCValue);
  AliZDCDigit(const AliZDCDigit & digit);

  // Getters 
  virtual Float_t   GetDetector() {return fDetector;}
  virtual Float_t   GetQuadrant() {return fQuadrant;}
  virtual Float_t   GetADCValue() {return fADCValue;}

  virtual ~AliZDCDigit(){} 

  // Operators
  Int_t operator == (AliZDCDigit &digit) {
    // Two digits are equal if they refers to the detector
    // in the same sub-volume (same procedure as for hits)
    if (fDetector != digit.fDetector) return 0;
    if (fQuadrant != digit.fQuadrant) return 0;
    return 1;
  }
  virtual AliZDCDigit& operator + (AliZDCDigit &digit) {
    // Adds the amplitude of digits 

    fADCValue += digit.fADCValue ;
    return *this ;
  }
  
 protected:

//  Int_t   fNprimary;          // Number of primaries
  Int_t   fDetector;          // Detector
  Int_t   fQuadrant;          // Quadrant
  Float_t fADCValue;          // ADC channel value

  // Print method
  virtual void Print(Option_t *) {
     printf(" -> DIGIT: Det =  %d Quad =  %d ADCCh =  %f\n ",
     fDetector, fQuadrant, fADCValue);
  }
    
  ClassDef(AliZDCDigit,1)   // Digits in ZDC 

} ;

#endif //  ALIZDCDIGIT_H
