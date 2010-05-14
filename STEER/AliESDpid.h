#ifndef ALIESDPID_H
#define ALIESDPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                    Combined PID class
//           for the Event Summary Data class
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------
#include <Rtypes.h>
#include "AliESDtrack.h" // Needed for inline functions
#include "AliTPCPIDResponse.h"
#include "AliITSPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliTRDPIDResponse.h"
//#include "HMPID/AliHMPID.h"
//#include "TRD/AliTRDpidESD.h"

class AliESDEvent;

class AliESDpid {
public:
  AliESDpid(): fRange(5.), fITSPIDmethod(kITSTruncMean), fTPCResponse(), fITSResponse(), fTOFResponse(), fTRDResponse(){;}
  virtual ~AliESDpid() {}
  Int_t MakePID(AliESDEvent *event, Bool_t TPCOnly = kFALSE, Float_t TimeZeroTOF=0) const;
  void MakeTPCPID(AliESDtrack *track) const;
  void MakeITSPID(AliESDtrack *track) const;
  void MakeTOFPID(AliESDtrack *track, Float_t TimeZeroTOF) const;
  //  void MakeHMPIDPID(AliESDtrack *track);
  void MakeTRDPID(AliESDtrack *track) const;
  void CombinePID(AliESDtrack *track) const;

  enum ITSPIDmethod { kITSTruncMean, kITSLikelihood };
  void SetITSPIDmethod(ITSPIDmethod pmeth) { fITSPIDmethod = pmeth; }

  Float_t NumberOfSigmasTPC(const AliESDtrack *track, AliPID::EParticleType type) const;
  Float_t NumberOfSigmasTOF(const AliESDtrack *track, AliPID::EParticleType type, const Float_t TimeZeroTOF) const;
  Float_t NumberOfSigmasITS(const AliESDtrack *track, AliPID::EParticleType type) const;

  AliITSPIDResponse &GetITSResponse() {return fITSResponse;}
  AliTPCPIDResponse &GetTPCResponse() {return fTPCResponse;}
  AliTOFPIDResponse &GetTOFResponse() {return fTOFResponse;}
  AliTRDPIDResponse &GetTRDResponse() {return fTRDResponse;}


private:
  Float_t           fRange;          // nSigma max in likelihood
  ITSPIDmethod      fITSPIDmethod;   // 0 = trunc mean; 1 = likelihood 
  AliTPCPIDResponse fTPCResponse;
  AliITSPIDResponse fITSResponse;
  AliTOFPIDResponse fTOFResponse;
  // AliHMPIDPIDResponse fHMPIDResponse;
  AliTRDPIDResponse fTRDResponse;

  ClassDef(AliESDpid,4)  // PID calculation class
};

inline Float_t AliESDpid::NumberOfSigmasTPC(const AliESDtrack *track, AliPID::EParticleType type) const {
  Double_t mom = track->GetP();
  const AliExternalTrackParam *in = track->GetInnerParam();
  if (in)
    mom = in->GetP();
  return fTPCResponse.GetNumberOfSigmas(mom,track->GetTPCsignal(),track->GetTPCsignalN(),type); 
}

inline Float_t AliESDpid::NumberOfSigmasTOF(const AliESDtrack *track, AliPID::EParticleType type, const Float_t TimeZeroTOF) const {
  Double_t times[AliPID::kSPECIES];
  track->GetIntegratedTimes(times);
  return (track->GetTOFsignal() - TimeZeroTOF - times[type])/fTOFResponse.GetExpectedSigma(track->GetP(),times[type],AliPID::ParticleMass(type));
}

inline Float_t AliESDpid::NumberOfSigmasITS(const AliESDtrack *track, AliPID::EParticleType type) const {
  return fITSResponse.GetNumberOfSigmas(track->GetP(),track->GetITSsignal(),type); 
}
#endif


