#ifndef ALITRDDIGIT_H
#define ALITRDDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdigit.h,v */

#include "AliDigitNew.h"

//_____________________________________________________________________________
class AliTRDdigit : public AliDigitNew {

 public:

  AliTRDdigit();
  AliTRDdigit(Bool_t isRaw, Int_t *digits, Int_t *amp);
  virtual ~AliTRDdigit();

  static  UInt_t RawDigit()          { return fgkRawDigit; };

          Int_t  GetAmp() const      { if (TestBit(fgkRawDigit))
                                         return DecodeAmp();
                                       else
                                         return fAmp; };
          Int_t  GetDetector() const { return fId;   };
          Int_t  GetRow() const      { return fRow;  };
          Int_t  GetCol() const      { return fCol;  };
          Int_t  GetTime() const     { return fTime; };

  virtual Int_t  DecodeAmp() const;

 protected:

  static const UInt_t fgkRawDigit; // Marks a raw digit

  UShort_t     fRow;               // Pad row number
  UShort_t     fCol;               // Pad col number
  UShort_t     fTime;              // Time bucket

  ClassDef(AliTRDdigit,2)          // Digit for the TRD

};

#endif
