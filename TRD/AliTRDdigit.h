#ifndef ALITRDDIGIT_H
#define ALITRDDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdigit.h,v */

#include "AliDigitNew.h"

const UInt_t kRawDigit = 0x00000001;

//_____________________________________________________________________________
class AliTRDdigit : public AliDigitNew {

 public:

  AliTRDdigit();
  AliTRDdigit(Bool_t isRaw, Int_t *digits, Int_t *amp);
  virtual ~AliTRDdigit();

          Int_t GetAmp() const    { if (TestBit(kRawDigit))
                                      return DecodeAmp();
                                    else
                                      return fAmp; };
          Int_t GetDetector()     { return fId;   };
          Int_t GetRow()          { return fRow;  };
          Int_t GetCol()          { return fCol;  };
          Int_t GetTime()         { return fTime; };

          Int_t DecodeAmp() const { return 0;     };

 protected:

  Int_t        fRow;              // Pad row number
  Int_t        fCol;              // Pad col number
  Int_t        fTime;             // Time bucket

  ClassDef(AliTRDdigit,1)         // Digit for the TRD

};

#endif
