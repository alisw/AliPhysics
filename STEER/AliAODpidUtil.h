#ifndef ALIAODPIDUTIL_H
#define ALIAODPIDUTIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAODpidUtil.h 38493 2010-01-26 16:33:03Z hristov $ */

//-------------------------------------------------------
//                    Combined PID class
//                    for the AOD class
//   Origin: Rosa Romita, GSI, r.romita@gsi.de 
//-------------------------------------------------------
#include <Rtypes.h>
#include <TMatrixD.h>
#include "AliAODTrack.h" // Needed for inline functions
#include "AliAODPid.h" // Needed for inline functions
#include "AliTPCPIDResponse.h"
#include "AliITSPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliTRDPIDResponse.h"
//#include "HMPID/AliHMPID.h"

class AliAODEvent;

class AliAODpidUtil {
public:
  AliAODpidUtil(): fRange(5.), fTPCResponse(), fITSResponse(), fTOFResponse(), fTRDResponse() {;}
  virtual ~AliAODpidUtil() {;}


  Int_t MakePID(AliAODTrack *track,Double_t *p) const;
  void MakeTPCPID(AliAODTrack *track,Double_t *p) const;
  void MakeITSPID(AliAODTrack *track,Double_t *p) const;
  void MakeTOFPID(AliAODTrack *track,Double_t *p) const;
  //  void MakeHMPIDPID(AliESDtrack *track);
  void MakeTRDPID(AliAODTrack *track,Double_t *p) const;

  Float_t NumberOfSigmasTPC(const AliAODTrack *track, AliPID::EParticleType type) const;
  Float_t NumberOfSigmasTOF(const AliAODTrack *track, AliPID::EParticleType type) const;
  Float_t NumberOfSigmasITS(const AliAODTrack *track, AliPID::EParticleType type) const;

  AliITSPIDResponse &GetITSResponse() {return fITSResponse;}
  AliTPCPIDResponse &GetTPCResponse() {return fTPCResponse;}
  AliTOFPIDResponse &GetTOFResponse() {return fTOFResponse;}

private:
  Float_t           fRange;          // nSigma max in likelihood
  AliTPCPIDResponse fTPCResponse;
  AliITSPIDResponse fITSResponse;
  AliTOFPIDResponse fTOFResponse;
  AliTRDPIDResponse fTRDResponse;

  ClassDef(AliAODpidUtil,1)  // PID calculation class
};

inline Float_t AliAODpidUtil::NumberOfSigmasTPC(const AliAODTrack *track, AliPID::EParticleType type) const {
  
  Double_t mom = track->P();
  AliAODPid *pidObj = track->GetDetPid();
  if (pidObj)
    mom = pidObj->GetTPCmomentum();
  return fTPCResponse.GetNumberOfSigmas(mom,pidObj->GetTPCsignal(),0,type); 
}

inline Float_t AliAODpidUtil::NumberOfSigmasTOF(const AliAODTrack *track, AliPID::EParticleType type) const {
  Double_t times[AliPID::kSPECIES];
  AliAODPid *pidObj = track->GetDetPid();
  pidObj->GetIntegratedTimes(times);
  return (pidObj->GetTOFsignal() - times[type])/fTOFResponse.GetExpectedSigma(track->P(),times[type],AliPID::ParticleMass(type));
}

inline Float_t AliAODpidUtil::NumberOfSigmasITS(const AliAODTrack *track, AliPID::EParticleType type) const {
  AliAODPid *pidObj = track->GetDetPid();
  Float_t dEdx=pidObj->GetITSsignal();
  UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
    if(clumap&(1<<i)) ++nPointsForPid;
  }
  if(track->GetTPCNcls()>0){
    Float_t mom=pidObj->GetTPCmomentum();
    return fITSResponse.GetNumberOfSigmas(mom,dEdx,type,nPointsForPid,kFALSE);
  }else{
    Float_t mom=track->P();
    return fITSResponse.GetNumberOfSigmas(mom,dEdx,type,nPointsForPid,kTRUE);
  }
}
#endif


