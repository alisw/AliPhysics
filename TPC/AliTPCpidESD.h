#ifndef ALITPCpIDESD_H
#define ALITPCpIDESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TPC PID class
// A very naive design... Should be made better by the detector experts...
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------
#include <Rtypes.h>

class AliESD;

class AliTPCpidESD {
public:
  AliTPCpidESD():fMIP(0.),fRes(0.),fRange(0.){}
  AliTPCpidESD(Double_t *param);
  virtual ~AliTPCpidESD() {}
  Int_t MakePID(AliESD *event);
  static Double_t Bethe(Double_t bg);
private:
  Double_t fMIP;          // dEdx for MIP
  Double_t fRes;          // relative dEdx resolution
  Double_t fRange;        // one particle type PID range (in sigmas)
  ClassDef(AliTPCpidESD,1)   // TPC PID class
};

#endif


