#ifndef ALITOFSDIGIT_H
#define ALITOFSDIGIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//                                            //
//  Class for TOF SDigits                     //
//                                            //
////////////////////////////////////////////////

/* $Id$ */

#include "TObject.h"
#include "TArrayF.h"
#include "TArrayI.h"

class AliTOFSDigit : public TObject {

  //overloading of the streamer << operator
  //friend ostream& operator << ( ostream& , const AliTOFSDigit&) ;

 public:
  AliTOFSDigit();
  AliTOFSDigit(Int_t tracknum, Int_t* vol, Int_t* digit);
// new ctor for sdigits
  AliTOFSDigit(Int_t sector, Int_t plate, Int_t strip, Int_t padx, Int_t padz, Int_t tdc, Int_t adc);
// copy ctor
  AliTOFSDigit(const AliTOFSDigit & digit) ;
  AliTOFSDigit& operator=(const AliTOFSDigit & digit) ;
  virtual ~AliTOFSDigit();
  void            GetLocation(Int_t* Loc) const;
  Int_t           GetTotPad() const;

  void Update(Float_t tdcbin, Int_t tdc, Int_t adc, Int_t track);
  void Update(AliTOFSDigit* sdig);

// getters for AliTOFSDigit object 
  Int_t   GetNDigits() const    {return fNDigits;}
  Int_t GetTdc(Int_t i) const {return fTdc->At(i);}
  Int_t GetAdc(Int_t i) const {return fAdc->At(i);}
//  Int_t   GetNTracks(Int_t i) const {return fTracks[i]->GetSize();}
  Int_t   GetTrack(Int_t i, Int_t j) const {return fTracks->At(i*kMAXDIGITS+j);}
  Int_t   GetSector() const     {return fSector;}
  Int_t   GetPlate()  const     {return fPlate;}
  Int_t   GetStrip()  const     {return fStrip;}
  Int_t   GetPadx()   const     {return fPadx;}
  Int_t   GetPadz()   const     {return fPadz;}

  enum {
    kMAXDIGITS = 3 // number 3 is a legacy from AliDigit object
  };

protected:

  Int_t   fSector;  // number of sector
  Int_t   fPlate;   // number of plate
  Int_t   fStrip;   // number of strip
  Int_t   fPadx;    // number of pad along x
  Int_t   fPadz;    // number of pad along z
  Int_t   fNDigits;  // dimension of fTdc array
  TArrayI *fTdc;     // tdc values for sdigit
  TArrayI *fAdc;     // adc values for sdigit
  TArrayI *fTracks;  // contributing tracks, kMAXDIGITS entries per
                     // 1 tdc value
  ClassDef(AliTOFSDigit,2)  // SDigit for Time Of Flight
};

#endif /* ALITOFSDIGIT_H */
