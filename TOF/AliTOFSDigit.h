#ifndef ALITOFSDIGIT_H
#define ALITOFSDIGIT_H

////////////////////////////////////////////////
//                                            //
//  Class for TOF SDigits                     //
//                                            //
////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TArrayF.h"
#include "TArrayI.h"
#include "AliDigit.h"

//class TArrayF;
class AliTOFGeometry;

// number 3 is a legacy from AliDigit object
const Int_t kMAXDIGITS = 3;

class AliTOFSDigit : public TObject {

  //overloading of the streamer << operator
  //friend ostream& operator << ( ostream& , const AliTOFSDigit&) ;

 public:
  AliTOFSDigit();
  AliTOFSDigit(Int_t tracknum, Int_t* vol, Float_t* digit);
// new ctor for sdigits
  AliTOFSDigit(Int_t sector, Int_t plate, Int_t strip, Int_t padx, Int_t padz, Float_t tdc, Float_t adc);
// copy ctor
  AliTOFSDigit(const AliTOFSDigit & digit) ;
  virtual ~AliTOFSDigit();
  void            GetLocation(Int_t* Loc) const;
  Int_t           GetTotPad(AliTOFGeometry *tofGeom) const;

  void Update(Float_t tdcbin, Int_t tdc, Int_t adc, Int_t track);
  void Update(AliTOFSDigit* sdig);

// getters for AliTOFSDigit object 
  Int_t   GetNDigits() const    {return fNDigits;}
  Float_t GetTdc(Int_t i) const {return fTdc->At(i);}
  Float_t GetAdc(Int_t i) const {return fAdc->At(i);}
//  Int_t   GetNTracks(Int_t i) const {return fTracks[i]->GetSize();}
  Int_t   GetTrack(Int_t i, Int_t j) const {return fTracks->At(i*kMAXDIGITS+j);}
  Int_t   GetSector() const     {return fSector;}
  Int_t   GetPlate()  const     {return fPlate;}
  Int_t   GetStrip()  const     {return fStrip;}
  Int_t   GetPadx()   const     {return fPadx;}
  Int_t   GetPadz()   const     {return fPadz;}

protected:
  Int_t   fSector;  // number of sector
  Int_t   fPlate;   // number of plate
  Int_t   fStrip;   // number of strip
  Int_t   fPadx;    // number of pad along x
  Int_t   fPadz;    // number of pad along z
  Int_t   fNDigits;  // dimension of fTdc array
  TArrayF *fTdc;     // tdc values for sdigit
  TArrayF *fAdc;     // adc values for sdigit
  TArrayI *fTracks;  // contributing tracks, kMAXDIGITS entries per
                     // 1 tdc value

//  Float_t *fTdc;    //[fNDigits] tdc values for sdigit
//  Float_t *fAdc;    //[fNDigits] adc values for sdigit
//  Int_t **fTracks;  //[fNDigits] contributing tracks, pointers to
		      //arrays with track indices

  ClassDef(AliTOFSDigit,1)  // SDigit for Time Of Flight
};

#endif /* ALITOFSDIGIT_H */
