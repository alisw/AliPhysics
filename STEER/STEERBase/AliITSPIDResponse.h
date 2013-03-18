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

class AliVTrack;

class AliITSPIDResponse : public TObject {

public:
  AliITSPIDResponse(Bool_t isMC=kFALSE);
  //AliITSPIDResponse(Double_t *param);
 ~AliITSPIDResponse() {}

 void SetBetheBlochParamsITSTPC(Double_t* param){
   for(Int_t iPar=0; iPar<5; iPar++) fBBtpcits[iPar]=param[iPar];
 }
 void SetBetheBlochParamsITSsa(Double_t* param){
   for(Int_t iPar=0; iPar<5; iPar++) fBBsa[iPar]=param[iPar];
 }
 
 void SetBetheBlochHybridParamsITSsa(Double_t* param){
   for(Int_t iPar=0; iPar<9; iPar++) fBBsaHybrid[iPar]=param[iPar];
 }
 void SetElectronBetheBlochParamsITSsa(Double_t* param){
   for(Int_t iPar=0; iPar<5; iPar++) fBBsaElectron[iPar]=param[iPar];
 }

 Double_t BetheAleph(Double_t p,Double_t mass) const;
 Double_t Bethe(Double_t p, Double_t mass, Bool_t iSA=kFALSE) const;
 Double_t BetheITSsaHybrid(Double_t p, Double_t mass) const;
 Double_t GetResolution(Double_t bethe, Int_t nPtsForPid=4, Bool_t isSA=kFALSE) const;
 void GetITSProbabilities(Float_t mom, Double_t qclu[4], Double_t condprobfun[AliPID::kSPECIES],Bool_t isMC=kFALSE) const;

 Double_t GetNumberOfSigmas( const AliVTrack* track, AliPID::EParticleType species) const;

 Double_t GetSignalDelta( const AliVTrack* track, AliPID::EParticleType species) const;
 
 Float_t GetNumberOfSigmas(Float_t mom, Float_t signal, AliPID::EParticleType type, Int_t nPtsForPid=4, Bool_t isSA=kFALSE) const {
   const Double_t chargeFactor = TMath::Power(AliPID::ParticleCharge(type),2.);
   Float_t bethe = Bethe(mom,AliPID::ParticleMassZ(type),isSA)*chargeFactor;
   return (signal - bethe)/GetResolution(bethe,nPtsForPid,isSA);
 }
 Int_t GetParticleIdFromdEdxVsP(Float_t mom, Float_t signal, Bool_t isSA=kFALSE) const;

private:


  // Data members for truncated mean method
  Float_t  fRes;             // relative dEdx resolution
  Double_t fKp1;             // ALEPH BB param 1
  Double_t fKp2;             // ALEPH BB param 2
  Double_t fKp3;             // ALEPH BB param 3
  Double_t fKp4;             // ALEPH BB param 4
  Double_t fKp5;             // ALEPH BB param 
  Double_t  fBBsa[5];        // parameters of BB for SA tracks
  Double_t  fBBsaHybrid[9];  // parameters of Hybrid BB for SA tracks, PHOB + Polinomial al low beta*gamma
  Double_t  fBBsaElectron[5];// parameters of electron BB for SA tracks
  Double_t  fBBtpcits[5];     // parameters of BB for TPC+ITS tracks
  Float_t  fResolSA[5];      // resolutions vs. n. of SDD/SSD points
  Float_t  fResolTPCITS[5];  // resolutions vs. n. of SDD/SSD points

  ClassDef(AliITSPIDResponse,4)   // ITS PID class
};

#endif


