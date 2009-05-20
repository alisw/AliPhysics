#ifndef ALITRDV0INFO_H
#define ALITRDV0INFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#ifndef Root_TObject
#include "TObject.h"
#endif

#ifndef ALITRDGEOMETRY_H
#include "AliTRDgeometry.h"
#endif

#ifndef ALIPID_H
#include "AliPID.h"
#endif

class AliESDv0;
class AliTRDv0Info : public TObject
{
public:
  enum ETRDv0Info{
    kNV0param = 10
   ,kNlayer   = AliTRDgeometry::kNlayer
   ,kNDetectors = 2
  };
  AliTRDv0Info();
  AliTRDv0Info(AliESDv0 *v0);
  void            Print(Option_t *opt=0x0) const;


private:
  Int_t   fStatus;              // track status
  Float_t fV0param[kNV0param];  // V0 parameters (momentum, variance, etc)
  Float_t fPplus[2*kNlayer];    // momentum and variance for the positive daughter  
  Float_t fPminus[2*kNlayer];   // momentum and variance for the negative daughter  
  Float_t fPID[kNDetectors][AliPID::kSPECIES]; // PID provided by TPC and TOF

  ClassDef(AliTRDv0Info, 0) // extracted V0 MC information
};


#endif

