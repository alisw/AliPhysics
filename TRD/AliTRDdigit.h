#ifndef ALITRDDIGIT_H
#define ALITRDDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDdigit.h,v */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  The TRD digit                                                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliDigitNew.h"

//_____________________________________________________________________________
class AliTRDdigit : public AliDigitNew {

 public:

  AliTRDdigit();
  AliTRDdigit(Int_t * const digits, const Int_t *amp);
  virtual ~AliTRDdigit();

          Int_t  GetDetector() const { return fId;   };

          Int_t  GetRow() const      { return fRow;  };
          Int_t  GetCol() const      { return fCol;  };
          Int_t  GetTime() const     { return fTime; };

 protected:

               UShort_t fRow;        // Pad row number
               UShort_t fCol;        // Pad col number
               UShort_t fTime;       // Time bucket

  ClassDef(AliTRDdigit,3)            // Digit for the TRD

};

#endif
