#ifndef ALIZDCDIGIT_H
#define ALIZDCDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//**********************************************************************
//
//   		Class for ZDC digit 
//   	      ADC Channels for each PM 
//   	   5 for hadronic ZDCs 1 for EM ZDCs
//
//**********************************************************************

#include<TObject.h>

class AliZDCDigit : public TObject {

 public:
  
  AliZDCDigit() ;
  AliZDCDigit(Int_t *Sector, Int_t *ADCValue);
  AliZDCDigit(const AliZDCDigit & digit);
  virtual ~AliZDCDigit() {}

  // Getters 
  Int_t   GetSector(Int_t i)	  {return fSector[i];}
  Int_t   GetADCValue(Int_t i)    {return fADCValue[i];}

  // Operators
  // Two digits are equal if they refers to the detector
  // in the same sub-volume (same procedure as for hits)
  Int_t operator == (AliZDCDigit &digit){
    Int_t i;
    for(i=0; i<2; i++) if(fSector[i]!=digit.GetSector(i)) return 0;
    return 1;
  }
  // Adds the amplitude of digits 
  virtual AliZDCDigit operator + (AliZDCDigit &digit){
    for(Int_t i = 0; i < 2; i++) fADCValue[i] += digit.fADCValue[i];
    return *this;
  }

  // Print method
  virtual void Print(Option_t *) const {
     printf("\t AliZDCDigit -> Detector %d Quadrant %d: ADC HighGain=  %d ADC LowGain=  %d\n ",
     fSector[0], fSector[1], fADCValue[0], fADCValue[1]);
  }
  
 protected:

  //Data members
  Int_t  fSector[2];         // Detector and tower in which light is produced
  Int_t  fADCValue[2];       // ADC channel value (0 = high gain, 1 = low gain)
    
  ClassDef(AliZDCDigit,4)   // Digits in ZDC 

} ;

#endif //  ALIZDCDIGIT_H

