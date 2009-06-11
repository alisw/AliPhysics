#ifndef ALITOFPIDESD_H
#define ALITOFPIDESD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TOF PID class
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include "TObject.h"
#include "AliPID.h"

class AliESDEvent;
class AliESDtrack;

class AliTOFpidESD : public TObject {
public:
  AliTOFpidESD();
  AliTOFpidESD(Double_t *param);
 ~AliTOFpidESD(){}
 
  void     SetTimeZero(Double_t t0) { fTime0=t0; }
  Double_t GetTimeZero() const { return fTime0; }

  void     SetMaxMismatchProbability(Double_t p) {fPmax=p;}
  Double_t GetMaxMismatchProbability() const {return fPmax;}

  Int_t MakePID(AliESDEvent *event);
  Int_t MakePID(AliESDEvent *event, Double_t timeZero);

  Bool_t ExpectedSignals(const AliESDtrack *t, 
                          Double_t s[], 
                          Int_t n=AliPID::kSPECIES) const;
  Bool_t ExpectedSigmas(const AliESDtrack *t, 
                         Double_t s[],
                         Int_t n=AliPID::kSPECIES) const;
  Bool_t NumberOfSigmas(const AliESDtrack *t, 
                         Double_t s[],
                         Int_t n=AliPID::kSPECIES) const;

  Double_t GetExpectedSignal(const AliESDtrack *t,
                     AliPID::EParticleType n=AliPID::kKaon) const;
  Double_t GetExpectedSigma(const AliESDtrack *t,
                     AliPID::EParticleType n=AliPID::kKaon) const;
  Double_t GetNumberOfSigmas(const AliESDtrack *t,
                     AliPID::EParticleType n=AliPID::kKaon) const;

private:
  Double_t GetMismatchProbability(Double_t p,Double_t mass) const;

  Double_t fSigma;        // intrinsic TOF resolution
  Double_t fRange;        // one particle type PID range (in sigmas)
  Double_t fPmax;         // "maximal" probability of mismathing (at ~0.5 GeV/c)
  Double_t fTime0;        // time zero

  ClassDef(AliTOFpidESD,3)   // TOF PID class
};

#endif
