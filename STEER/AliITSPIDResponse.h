#ifndef ALIITSPIDRESPONSE_H
#define ALIITSPIDRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                    ITS PID response class
//
//
//-------------------------------------------------------
//#include <Rtypes.h>
#include <TObject.h>

#include "AliPID.h"

class AliITSPIDResponse : public TObject {

public:
  AliITSPIDResponse();
  AliITSPIDResponse(Double_t *param);
 ~AliITSPIDResponse() {}
  Double_t Bethe(Double_t p,Double_t mass) const;
  Double_t GetResolution(Double_t bethe) const;
  void GetITSProbabilities(Float_t mom, Double_t qclu[4], Double_t condprobfun[AliPID::kSPECIES]) const;
  Float_t GetNumberOfSigmas(Float_t mom, Float_t signal, AliPID::EParticleType type) const {
  Float_t bethe = Bethe(mom,AliPID::ParticleMass(type));
  return (signal - bethe)/GetResolution(bethe);
}

private:

  // Data members for truncated mean method
  Float_t fRes;             // relative dEdx resolution
  Double_t fKp1;             // ALEPH BB param 1
  Double_t fKp2;             // ALEPH BB param 2
  Double_t fKp3;             // ALEPH BB param 3
  Double_t fKp4;             // ALEPH BB param 4
  Double_t fKp5;             // ALEPH BB param 5

  ClassDef(AliITSPIDResponse,1)   // ITS PID class
};

#endif


