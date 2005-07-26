////////////////////////////////////////////////
//                                            //
//  Digit class for TOF                       //
//  Interface                                 //
//  Getters, Setters and member variables     //
//  declared here                             //
//                                            //
////////////////////////////////////////////////

#ifndef ALITOFDIGIT_H
#define ALITOFDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliDigit.h"

class AliTOFdigit : public AliDigit {

  //overloading of the streamer << operator
friend ostream& operator << ( ostream& , const AliTOFdigit&) ;

 public:
  AliTOFdigit(){}
  AliTOFdigit(Int_t* tracks, Int_t* vol, Float_t* digit);
// new ctor for sdigits
  AliTOFdigit(Int_t sector, Int_t plate, Int_t strip, Int_t padx, Int_t padz, Float_t tdc, Float_t adc);
// copy ctor
  AliTOFdigit(const AliTOFdigit & digit) ;
  virtual ~AliTOFdigit(){}
  void            GetLocation(Int_t* Loc) const;
  Int_t           GetTotPad() const;
  void            AddTrack(Int_t track);
  // getters for AliTOFdigit object 
  Float_t GetTdc()    const     {return fTdc;}
  Float_t GetAdc()    const     {return fAdc;}
  Int_t   GetSector() const     {return fSector;}
  Int_t   GetPlate()  const     {return fPlate;}
  Int_t   GetStrip()  const     {return fStrip;}
  Int_t   GetPadx()   const     {return fPadx;}
  Int_t   GetPadz()   const     {return fPadz;}

  // setters for AliTOFdigit object
  void    SetTdc(Float_t TDC){fTdc = TDC;}
  void    SetAdc(Float_t ADC){fAdc = ADC;}

  //overloading of ==, + operators (summable digits)
  
  Bool_t operator==(const AliTOFdigit& digit) const;
  AliTOFdigit& operator+(AliTOFdigit const &digit) ;  


protected:
  Int_t   fSector;  // number of sector
  Int_t   fPlate;   // number of plate
  Int_t   fStrip;   // number of strip
  Int_t   fPadx;    // number of pad along x
  Int_t   fPadz;    // number of pad along z
  Float_t fTdc;     // tdc channel value, to be multiplied by
		    // AliTOFGeometry::TdcBinWidth() to have the
		    // time-of-flight measurement
  Float_t fAdc;     // adc channel value, to be multiplie by
		    // AliTOFSDigitizer::GetAdcBin() to have the
		    // 'charge' measurement

  ClassDef(AliTOFdigit,2)  // Digit for Time Of Flight
};

#endif /* ALITOFDIGIT_H */
