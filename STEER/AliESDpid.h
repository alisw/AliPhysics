#ifndef ALIESDPID_H
#define ALIESDPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                    Combined PID class
//           for the Event Summary Data class
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//   Modified: Jens Wiechula, Uni Tuebingen, jens.wiechula@cern.ch
//-------------------------------------------------------
#include <Rtypes.h>
#include "AliESDtrack.h" // Needed for inline functions

//#include "HMPID/AliHMPID.h"
//#include "TRD/AliTRDpidESD.h"

#include "AliPIDResponse.h"

class AliESDEvent;
class AliVEvent;
class AliVParticle;

class AliESDpid : public AliPIDResponse  {
public:
  AliESDpid(Bool_t forMC=kFALSE): AliPIDResponse(forMC), fRange(5.), fRangeTOFMismatch(5.), fITSPIDmethod(kITSTruncMean) {;}
  virtual ~AliESDpid() {}
  
  Int_t MakePID(AliESDEvent *event, Bool_t TPCOnly = kFALSE, Float_t timeZeroTOF=9999) const;
  void MakeTPCPID(AliESDtrack *track) const;
  void MakeITSPID(AliESDtrack *track) const;
  void MakeTOFPID(AliESDtrack *track, Float_t /*timeZeroTOF*/) const;
  Bool_t CheckTOFMatching(AliESDtrack *track) const;
  //  void MakeHMPIDPID(AliESDtrack *track);
  void MakeTRDPID(AliESDtrack *track) const;
  void CombinePID(AliESDtrack *track) const;

  enum ITSPIDmethod { kITSTruncMean, kITSLikelihood };
  void SetITSPIDmethod(ITSPIDmethod pmeth) { fITSPIDmethod = pmeth; }
  
  virtual Float_t NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type, const Float_t timeZeroTOF) const;
  virtual Float_t NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type) const {return NumberOfSigmasTOF(vtrack,type,0); }
  
  void SetTOFResponse(AliVEvent *vevent,EStartTimeType_t option);

  void SetNMaxSigmaTOFTPCMismatch(Float_t range) {fRangeTOFMismatch=range;}
  Float_t GetNMaxSigmaTOFTPCMismatch() const {return fRangeTOFMismatch;}

private:
  Float_t           fRange;          // nSigma max in likelihood
  Float_t           fRangeTOFMismatch; // nSigma max for TOF matching with TPC
  ITSPIDmethod      fITSPIDmethod;   // 0 = trunc mean; 1 = likelihood 

  ClassDef(AliESDpid,6)  // PID calculation class
};


inline Float_t AliESDpid::NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type, const Float_t /*timeZeroTOF*/) const {
  AliESDtrack *track=(AliESDtrack*)vtrack;
  Double_t times[AliPID::kSPECIES];
  track->GetIntegratedTimes(times);
  return (track->GetTOFsignal() - fTOFResponse.GetStartTime(track->GetP()) - times[type])/fTOFResponse.GetExpectedSigma(track->GetP(),times[type],AliPID::ParticleMass(type));
}

#endif


