#ifndef ALITRDPIDESD_H
#define ALITRDPIDESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TRD PID class
// A very naive design... And the implementation is even poorer... 
// Should be made better by the detector experts...
//-------------------------------------------------------
#include <Rtypes.h>

class AliESD;

class AliTRDpidESD {
public:
  AliTRDpidESD(Double_t *param);
  Int_t MakePID(AliESD *event);
  static Double_t Bethe(Double_t bg);
private:
  Double_t fMIP;          // dEdx for MIP
  Double_t fRes;          // relative dEdx resolution
  Double_t fRange;        // one particle type PID range (in sigmas)
  ClassDef(AliTRDpidESD,1)   // TRD PID class
};

#endif


