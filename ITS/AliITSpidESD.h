#ifndef ALIITSpIDESD_H
#define ALIITSpIDESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    ITS PID class
// A very naive design... Should be made better by the detector experts...
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------
#include <Rtypes.h>

class AliESD;

class AliITSpidESD {
public:
  AliITSpidESD(Double_t *param);
  virtual ~AliITSpidESD() {}
  Int_t MakePID(AliESD *event);
  static Double_t Bethe(Double_t bg);
private:
  Double_t fMIP;          // dEdx for MIP
  Double_t fRes;          // relative dEdx resolution
  Double_t fRange;        // one particle type PID range (in sigmas)
  ClassDef(AliITSpidESD,1)   // ITS PID class
};

#endif


