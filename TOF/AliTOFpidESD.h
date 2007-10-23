#ifndef ALITOFPIDESD_H
#define ALITOFPIDESD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TOF PID class
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include "TObject.h"

class AliESDEvent;

class AliTOFGeometry;

class AliTOFpidESD : public TObject {
public:
  AliTOFpidESD();
  AliTOFpidESD(Double_t *param);
 ~AliTOFpidESD(){}
 
  void     SetMaxMismatchProbability(Double_t p) {fPmax=p;}
  Double_t GetMaxMismatchProbability() {return fPmax;}

  Int_t MakePID(AliESDEvent *event);
  Int_t MakePID(AliESDEvent *event, Double_t timeZero);

private:
  Double_t GetMismatchProbability(Double_t p,Double_t mass) const;

  Double_t fSigma;        // intrinsic TOF resolution
  Double_t fRange;        // one particle type PID range (in sigmas)
  Double_t fPmax;         // "maximal" probability of mismathing (at ~0.5 GeV/c)

  ClassDef(AliTOFpidESD,2)   // TOF PID class
};

#endif
