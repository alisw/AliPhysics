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
#include "AliMCEventHandler.h"

//#include "HMPID/AliHMPID.h"
//#include "TRD/AliTRDpidESD.h"

#include "AliPIDResponse.h"

class AliESDEvent;
class AliVEvent;
class AliVParticle;

class AliESDpid : public AliPIDResponse  {
public:
  AliESDpid(Bool_t forMC=kFALSE): AliPIDResponse(forMC), fRangeTOFMismatch(5.), fEventHandler(NULL) {;}
AliESDpid(const AliESDpid&a): AliPIDResponse(a), fRangeTOFMismatch(a.fRangeTOFMismatch), fEventHandler(NULL){;};
    AliESDpid& operator=(const AliESDpid& a){AliPIDResponse::operator=(a); fRangeTOFMismatch=a.fRangeTOFMismatch; fEventHandler=NULL; return *this;};
  virtual ~AliESDpid() {}
  
  Int_t MakePID(AliESDEvent *event, Bool_t TPCOnly = kFALSE, Float_t timeZeroTOF=9999) const;
  void MakeTPCPID(AliESDtrack *track) const;
  void MakeITSPID(AliESDtrack *track) const;
  void MakeTOFPID(AliESDtrack *track, Float_t /*timeZeroTOF*/) const;
  Bool_t CheckTOFMatching(AliESDtrack *track) const;
  //  void MakeHMPIDPID(AliESDtrack *track);
  void MakeTRDPID(AliESDtrack *track) const;
  void CombinePID(AliESDtrack *track) const;

//   Float_t NumberOfSigmasTOF(const AliVParticle *track, AliPID::EParticleType type) const {return AliPIDResponse::NumberOfSigmasTOF(track,type);}
//   Float_t GetNumberOfSigmasTOF(const AliVParticle *track, AliPID::EParticleType type, const Float_t timeZeroTOF) const;
  
  void SetNMaxSigmaTOFTPCMismatch(Float_t range) {fRangeTOFMismatch=range;}
  Float_t GetNMaxSigmaTOFTPCMismatch() const {return fRangeTOFMismatch;}

  Float_t GetTPCsignalTunedOnData(const AliVTrack *t) const;

  void SetEventHandler(AliVEventHandler *event){fEventHandler=event;};
protected:
  virtual Float_t GetSignalDeltaTOFold(const AliVParticle *track, AliPID::EParticleType type) const;
  virtual Float_t GetNumberOfSigmasTOFold(const AliVParticle *track, AliPID::EParticleType type) const;

private:

  Float_t           fRangeTOFMismatch; // nSigma max for TOF matching with TPC
  AliVEventHandler *fEventHandler; //! MC event handler
  
  ClassDef(AliESDpid,7)  // PID calculation class
};


#endif


