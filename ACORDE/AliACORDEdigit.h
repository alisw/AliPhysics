#ifndef ALIACORDEDIGIT_H
#define ALIACORDEDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//  ACORDE digit: Id
//
// The digits are made in FinishEvent() by summing all the hits in a 
// counter.
////////////////////////////////////////////////////////////////////////////

#include "AliDigit.h"

class TArrayF;
class TArrayI;

class AliACORDEdigit: public AliDigit  {
public:
  AliACORDEdigit();
  AliACORDEdigit(Int_t* tracks, Int_t* vol, Float_t* digit);
  AliACORDEdigit(const AliACORDEdigit& digit);
  virtual ~AliACORDEdigit();

  AliACORDEdigit& operator= (const AliACORDEdigit& digit);

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
    ClassDef(AliACORDEdigit,1)  //Digit (Header) object for set : ACORDE (ACORDE)
};

typedef AliACORDEdigit AliCRTdigit; // for backward compatibility

#endif // ALIACORDEDIGIT_H
