#ifndef ALIAODPIDUTIL_H
#define ALIAODPIDUTIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAODpidUtil.h 38493 2010-01-26 16:33:03Z hristov $ */

//-------------------------------------------------------
//                    Combined PID class
//                    for the AOD class
//   Origin: Rosa Romita, GSI, r.romita@gsi.de 
//   Modified: Jens Wiechula, Uni Tuebingen, jens.wiechula@cern.ch
//-------------------------------------------------------
#include <Rtypes.h>
#include <TMatrixD.h>
#include "AliAODTrack.h" // Needed for inline functions
#include "AliAODPid.h" // Needed for inline functions
//#include "HMPID/AliHMPID.h"

#include "AliPIDResponse.h"

class AliAODEvent;
class AliVParticle;

class AliAODpidUtil : public AliPIDResponse  {
public:
  //TODO: isMC???
  AliAODpidUtil(Bool_t isMC = kFALSE): AliPIDResponse(isMC), fRange(5.) {;}
  virtual ~AliAODpidUtil() {;}

  Int_t MakePID(const AliAODTrack *track,Double_t *p) const;
  void MakeTPCPID(const AliAODTrack *track,Double_t *p) const;
  void MakeITSPID(const AliAODTrack *track,Double_t *p) const;
  void MakeTOFPID(const AliAODTrack *track,Double_t *p) const;
  //  void MakeHMPIDPID(AliESDtrack *track);
  void MakeTRDPID(const AliAODTrack *track,Double_t *p) const;

  virtual Float_t NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type) const;
  
private:
  Float_t           fRange;          // nSigma max in likelihood
  
  ClassDef(AliAODpidUtil,2)  // PID calculation class
};

inline Float_t AliAODpidUtil::NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type) const {
  AliAODTrack *track=(AliAODTrack*)vtrack;
  Double_t times[AliPID::kSPECIES];
  Double_t sigmaTOFPid[AliPID::kSPECIES];
  AliAODPid *pidObj = track->GetDetPid();
  if (!pidObj) return -999.;
  pidObj->GetIntegratedTimes(times);
  pidObj->GetTOFpidResolution(sigmaTOFPid);
  if (sigmaTOFPid[type]>0) return (pidObj->GetTOFsignal() - times[type])/sigmaTOFPid[type];
  else return (pidObj->GetTOFsignal() - times[type])/fTOFResponse.GetExpectedSigma(track->P(),times[type],AliPID::ParticleMass(type));
}

#endif


