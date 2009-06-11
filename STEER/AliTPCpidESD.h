#ifndef ALITPCpIDESD_H
#define ALITPCpIDESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TPC PID class
// A very naive design... Should be made better by the detector experts...
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
// With many additions and modifications suggested by
//      Alexander Kalweit, GSI, alexander.philipp.kalweit@cern.ch
//      Dariusz Miskowiec, GSI, D.Miskowiec@gsi.de
//-------------------------------------------------------
#include <Rtypes.h>

#include "AliPID.h"

class AliESDEvent;
class AliESDtrack;

class AliTPCpidESD {
public:
  AliTPCpidESD();
  AliTPCpidESD(Double_t *param);
  virtual ~AliTPCpidESD() {}
  void SetBetheBlochParameters(Double_t kp1,
                               Double_t kp2,
                               Double_t kp3,
                               Double_t kp4,
                               Double_t kp5
                               );
  Int_t MakePID(AliESDEvent *event);
  Double_t Bethe(Double_t bg) const;

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
  Double_t fMIP;          // dEdx for MIP
  Double_t fRes;          // relative dEdx resolution
  Double_t fRange;        // one particle type PID range (in sigmas)

  Double_t fKp1;   // Parameters
  Double_t fKp2;   //    of
  Double_t fKp3;   // the ALEPH
  Double_t fKp4;   // Bethe-Bloch
  Double_t fKp5;   // formula

  ClassDef(AliTPCpidESD,2)   // TPC PID class
};

#endif


