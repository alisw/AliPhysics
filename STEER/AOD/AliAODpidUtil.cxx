/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAODpidUtil.cxx 38329 2010-01-17 19:17:24Z hristov $ */

//-----------------------------------------------------------------
//           Implementation of the combined PID class
//           For the AOD Class
//           containing information on the particle identification
//      Origin: Rosa Romita, GSI, r.romita@gsi.de
//-----------------------------------------------------------------

#include "TRandom.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliAODpidUtil.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliTRDPIDResponse.h"
#include "AliESDtrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

#include <AliDetectorPID.h>

ClassImp(AliAODpidUtil)

//_________________________________________________________________________
Float_t AliAODpidUtil::GetTPCsignalTunedOnData(const AliVTrack *t) const {
    AliAODTrack *track = (AliAODTrack *) t;
    Float_t dedx = track->GetTPCsignalTunedOnData();
    if(dedx > 0) return dedx;

    dedx = t->GetTPCsignal();
    track->SetTPCsignalTunedOnData(dedx);

    if(dedx < 20) return dedx;

    
    AliPID::EParticleType type = AliPID::kPion;
    
    AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(track->GetAODEvent()->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (mcHeader) {
	TClonesArray *mcArray = (TClonesArray*)track->GetAODEvent()->GetList()->FindObject(AliAODMCParticle::StdBranchName());
	
	Bool_t kGood = kTRUE;
	
	Int_t iS = TMath::Abs(((AliAODMCParticle*)mcArray->At(TMath::Abs(t->GetLabel())))->GetPdgCode());
	if(iS==AliPID::ParticleCode(AliPID::kElectron)){
	    type = AliPID::kElectron;
	}
	else if(iS==AliPID::ParticleCode(AliPID::kMuon)){
	    type = AliPID::kMuon;
	}
	else if(iS==AliPID::ParticleCode(AliPID::kPion)){
	    type = AliPID::kPion;
	}
	else if(iS==AliPID::ParticleCode(AliPID::kKaon)){
	    type = AliPID::kKaon;
	}
	else if(iS==AliPID::ParticleCode(AliPID::kProton)){
	    type = AliPID::kProton;
	}
	else if(iS==AliPID::ParticleCode(AliPID::kDeuteron)){ // d
	    type = AliPID::kDeuteron;
	}
	else if(iS==AliPID::ParticleCode(AliPID::kTriton)){ // t
	    type = AliPID::kTriton;
	}
	else if(iS==AliPID::ParticleCode(AliPID::kHe3)){ // 3He
	    type = AliPID::kHe3;
	}
	else if(iS==AliPID::ParticleCode(AliPID::kAlpha)){ // 4He
	    type = AliPID::kAlpha;
	}
	else
	    kGood = kFALSE;

	if(kGood){
	    //TODO maybe introduce different dEdxSources?
        Double_t bethe = fTPCResponse.GetExpectedSignal(track, type, AliTPCPIDResponse::kdEdxDefault, this->UseTPCEtaCorrection());
        Double_t sigma = fTPCResponse.GetExpectedSigma(track, type, AliTPCPIDResponse::kdEdxDefault, this->UseTPCEtaCorrection());
        dedx = gRandom->Gaus(bethe,sigma);
        
// 	    if(iS == AliPID::ParticleCode(AliPID::kHe3) || iS == AliPID::ParticleCode(AliPID::kAlpha)) dedx *= 5;
	}

    }

    track->SetTPCsignalTunedOnData(dedx);
    return dedx;
}

//_________________________________________________________________________
Float_t AliAODpidUtil::GetSignalDeltaTOFold(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Number of sigma implementation for the TOF
  //
  
  AliAODTrack *track=(AliAODTrack*)vtrack;
  
  AliAODPid *pidObj = track->GetDetPid();
  if (!pidObj) return -9999.;
  Double_t tofTime=pidObj->GetTOFsignal();
  const Double_t expTime=fTOFResponse.GetExpectedSignal((AliVTrack*)vtrack,type);
  Double_t sigmaTOFPid[AliPID::kSPECIES];
  pidObj->GetTOFpidResolution(sigmaTOFPid);
  AliAODEvent *event=(AliAODEvent*)track->GetAODEvent();
  if (event) {  // protection if the user didn't call GetTrack, which sets the internal pointer
    AliTOFHeader* tofH=(AliTOFHeader*)event->GetTOFHeader();
    if (tofH && (TMath::Abs(sigmaTOFPid[0]) <= 1.E-16) ) { // new AOD
        tofTime -= fTOFResponse.GetStartTime(vtrack->P());
    }
  } else {
    AliError("pointer to AliAODEvent not found, please call GetTrack to set it");
    return -9999.;
  }
  
  return tofTime - expTime;
}

//_________________________________________________________________________
Float_t AliAODpidUtil::GetNumberOfSigmasTOFold(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Number of sigma implementation for the TOF
  //
  
  AliAODTrack *track=(AliAODTrack*)vtrack;

  Bool_t oldAod=kTRUE;
  Double_t sigTOF=0.;
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
