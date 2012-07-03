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
//   Modified: Pietro Antonioli, INFN BO, pietro.antonioli@bo.infn.it
//-------------------------------------------------------
#include <Rtypes.h>
#include <TMatrixD.h>
#include <AliLog.h>
#include "AliAODEvent.h" // Needed for inline functions
#include "AliAODTrack.h" // Needed for inline functions
#include "AliAODPid.h" // Needed for inline functions
#include "AliTOFHeader.h" //Needed for inline functions
//#include "HMPID/AliHMPID.h"

#include "AliPIDResponse.h"

class AliAODEvent;
class AliVParticle;

class AliAODpidUtil : public AliPIDResponse  {
public:
  //TODO: isMC???
  AliAODpidUtil(Bool_t isMC = kFALSE): AliPIDResponse(isMC) {;}
  virtual ~AliAODpidUtil() {;}

  Int_t MakePID(const AliAODTrack *track,Double_t *p) const;
  void MakeTPCPID(const AliAODTrack *track,Double_t *p) const;
  void MakeITSPID(const AliAODTrack *track,Double_t *p) const;
  void MakeTOFPID(const AliAODTrack *track,Double_t *p) const;
  //  void MakeHMPIDPID(AliESDtrack *track);
  void MakeTRDPID(const AliAODTrack *track,Double_t *p) const;

  virtual Float_t NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type) const;
  
private:
  
  ClassDef(AliAODpidUtil,3)  // PID calculation class
};

inline Float_t AliAODpidUtil::NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type) const {
  AliAODTrack *track=(AliAODTrack*)vtrack;
  Bool_t oldAod=kTRUE;
  Double_t sigTOF;
  AliAODPid *pidObj = track->GetDetPid();
  if (!pidObj) return -999.;
  Double_t tofTime=pidObj->GetTOFsignal();
  Double_t expTime=fTOFResponse.GetExpectedSignal((AliVTrack*)vtrack,type);
  Double_t sigmaTOFPid[AliPID::kSPECIES];
  pidObj->GetTOFpidResolution(sigmaTOFPid);
  AliAODEvent *event=(AliAODEvent*)track->GetAODEvent();
  if (event) {  // protection if the user didn't call GetTrack, which sets the internal pointer
    AliTOFHeader* tofH=(AliTOFHeader*)event->GetTOFHeader();
    if (tofH && (TMath::Abs(sigmaTOFPid[0]) <= 1.E-16) ) { // new AOD
	sigTOF=fTOFResponse.GetExpectedSigma(track->P(),expTime,AliPID::ParticleMassZ(type)); //fTOFResponse is set in InitialiseEvent
	tofTime -= fTOFResponse.GetStartTime(vtrack->P());
	oldAod=kFALSE;
    } 
  } else {
    AliError("pointer to AliAODEvent not found, please call GetTrack to set it");
    return -996.;
  }
  if (oldAod) { // old AOD
    if (type <= AliPID::kProton) {
      sigTOF=sigmaTOFPid[type];
    } else return -998.;  // light nuclei cannot be supported on old AOD because we don't have timeZero resolution
  }
  if (sigTOF>0) return (tofTime - expTime)/sigTOF;
  else return -997.;
}

#endif


