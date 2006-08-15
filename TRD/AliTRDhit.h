#ifndef ALITRDHIT_H
#define ALITRDHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Hit class for the TRD                                                 //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliHit.h"

class AliTRDhit : public AliHit {

 public:

  AliTRDhit();
  AliTRDhit(Int_t shunt, Int_t track, Int_t det, Float_t *hits, Int_t q);
  virtual ~AliTRDhit();

          Int_t    GetDetector() const         { return fDetector; }
          Int_t    GetCharge() const           { return fQ;        }

          void     SetX(Float_t x)             { fX        =   x;  }
          void     SetY(Float_t y)             { fY        =   y;  }
          void     SetZ(Float_t z)             { fZ        =   z;  }
          void     SetDetector(Int_t det)      { fDetector = det;  }
          void     SetQ(Int_t q)               { fQ        =   q;  }

          Bool_t   FromDrift() const           { return TestBit(kDrift);         }
          Bool_t   FromAmplification() const   { return TestBit(kAmplification); }
          Bool_t   FromTRphoton() const        { return (fQ < 0);                }
	  
          void     SetDrift()                  { SetBit(kDrift);                 }
          void     SetAmplification()          { SetBit(kAmplification);         }
          void     SetTRphoton()               { SetBit(kTRphoton);              }

 protected:

          enum {
              kDrift         = 0x00000001   //  Hit is from the drift region
            , kAmplification = 0x00000002   //  Hit is from the amplification region
            , kTRphoton      = 0x00000004   //  Hit is from a TR photon
            , kTest          = 0x00000008   //  Hit is a special test hit
          };

          UShort_t fDetector;               //  TRD detector number
          Short_t  fQ;                      //  Charge created by a hit. TR signals are negative.

  ClassDef(AliTRDhit,3)                     //  Hit for the Transition Radiation Detector

};

#endif
