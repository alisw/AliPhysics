#ifndef ALITOFPIDRESPONSE_H
#define ALITOFPIDRESPONSE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TOF PID class
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include "TObject.h"
#include "AliPID.h"

class AliTOFPIDResponse : public TObject {
public:
  AliTOFPIDResponse();
  AliTOFPIDResponse(Double_t *param);
 ~AliTOFPIDResponse(){}

  void     SetTimeResolution(Float_t res) { fSigma = res; }
  void     SetTimeZero(Double_t t0) { fTime0=t0; }
  Double_t GetTimeZero() const { return fTime0; }

  void     SetMaxMismatchProbability(Double_t p) {fPmax=p;}
  Double_t GetMaxMismatchProbability() const {return fPmax;}

  Double_t GetExpectedSigma(Float_t mom, Float_t tof, Float_t mass) const;


  Double_t GetMismatchProbability(Double_t p,Double_t mass) const;

 private:
  Double_t fSigma;        // intrinsic TOF resolution
  Double_t fPmax;         // "maximal" probability of mismathing (at ~0.5 GeV/c)
  Double_t fTime0;        // time zero

  ClassDef(AliTOFPIDResponse,1)   // TOF PID class
};

#endif
