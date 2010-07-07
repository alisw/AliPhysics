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
  AliITSPIDResponse(Bool_t isMC=kFALSE);
  AliITSPIDResponse(Double_t *param);
 ~AliITSPIDResponse() {}
 Double_t BetheAleph(Double_t p,Double_t mass) const;
 Double_t Bethe(Double_t p, Double_t mass, Bool_t iSA=kFALSE) const;
 Double_t GetResolution(Double_t bethe, Int_t nPtsForPid=4, Bool_t isSA=kFALSE) const;
 void GetITSProbabilities(Float_t mom, Double_t qclu[4], Double_t condprobfun[AliPID::kSPECIES]) const;
 Float_t GetNumberOfSigmas(Float_t mom, Float_t signal, AliPID::EParticleType type, Int_t nPtsForPid=4, Bool_t isSA=kFALSE) const {
   Float_t bethe = Bethe(mom,AliPID::ParticleMass(type),isSA);
   return (signal - bethe)/GetResolution(bethe,nPtsForPid,isSA);
 }

private:


  // Data members for truncated mean method
  Float_t  fRes;             // relative dEdx resolution
  Double_t fKp1;             // ALEPH BB param 1
  Double_t fKp2;             // ALEPH BB param 2
  Double_t fKp3;             // ALEPH BB param 3
  Double_t fKp4;             // ALEPH BB param 4
  Double_t fKp5;             // ALEPH BB param 
  Double_t  fBBsa[5];         // parameters of BB for SA tracks
  Double_t  fBBtpcits[5];     // parameters of BB for TPC+ITS tracks
  Float_t  fResolSA[5];      // resolutions vs. n. of SDD/SSD points
  Float_t  fResolTPCITS[5];  // resolutions vs. n. of SDD/SSD points

  ClassDef(AliITSPIDResponse,2)   // ITS PID class
};

#endif


