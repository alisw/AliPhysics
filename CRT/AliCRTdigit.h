#ifndef ALICRTDIGIT_H
#define ALICRTDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//  CRT digit: Id
//
// The digits are made in FinishEvent() by summing all the hits in a 
// counter.
////////////////////////////////////////////////////////////////////////////

#include "AliDigit.h"

class TArrayF;
class TArrayI;

class AliCRTdigit: public AliDigit  {
public:
  AliCRTdigit();
  AliCRTdigit(Int_t* tracks, Int_t* vol, Float_t* digit);
  AliCRTdigit(const AliCRTdigit& digit);
  virtual ~AliCRTdigit();

  AliCRTdigit& operator= (const AliCRTdigit& digit);

protected:
  Int_t     fSector;  // number of sector
  Int_t     fPlate;   // number of plate
  Int_t     fStrip;   // number of strip
  Int_t     fPadx;    // number of pad along x
  Int_t     fPadz;    // number of pad along z
  Int_t     fNDigits;  // dimension of fTdc array
  TArrayF*  fTdc;     // tdc values for sdigit
  TArrayF*  fAdc;     // adc values for sdigit

private:
    ClassDef(AliCRTdigit,1)  //Digit (Header) object for set : CRT (ACORDE)
};
#endif // ALICRTDIGIT_H
