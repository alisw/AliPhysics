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

  Int_t AliAODpidUtil::MakePID(const AliAODTrack *track,Double_t *p) const {
  //
  //  Calculate probabilities for all detectors, except if TPConly==kTRUE
  //  and combine PID
  //  
  //   Option TPConly==kTRUE is used during reconstruction, 
  //  because ITS tracking uses TPC pid
  //  HMPID and TRD pid are done in detector reconstructors
  //

  /*
    Float_t TimeZeroTOF = 0;
    if (subtractT0) 
    TimeZeroTOF = event->GetT0();
  */
  Int_t ns=AliPID::kSPECIES;
  Double_t tpcPid[AliPID::kSPECIES];
  MakeTPCPID(track,tpcPid);
  Double_t itsPid[AliPID::kSPECIES];
  Double_t tofPid[AliPID::kSPECIES];
  Double_t trdPid[AliPID::kSPECIES];
  MakeITSPID(track,itsPid);
  MakeTOFPID(track,tofPid);
  //MakeHMPIDPID(track);
  MakeTRDPID(track,trdPid);
  for (Int_t j=0; j<ns; j++) {
    p[j]=tpcPid[j]*itsPid[j]*tofPid[j]*trdPid[j];
  }

  return 0;
}
//_________________________________________________________________________
Float_t AliAODpidUtil::GetTPCsignalTunedOnData(const AliVTrack *t) const {
    AliAODTrack *track = (AliAODTrack *) t;
    Float_t dedx = track->GetTPCsignalTunedOnData();
    if(dedx > 0) return dedx;

    Double_t mom = t->GetTPCmomentum();

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
	    Double_t bethe=fTPCResponse.GetExpectedSignal(mom,type);
	    Double_t sigma=fTPCResponse.GetExpectedSigma(mom,t->GetTPCsignalN(),type);
	    dedx = gRandom->Gaus(bethe,sigma);

	    if(iS == AliPID::ParticleCode(AliPID::kHe3) || iS == AliPID::ParticleCode(AliPID::kAlpha)) dedx *= 5;
	}

    }

    track->SetTPCsignalTunedOnData(dedx);
    return dedx;
}
//_________________________________________________________________________
void AliAODpidUtil::MakeTPCPID(const AliAODTrack *track,Double_t *p) const
{
  //
  //  TPC pid using bethe-bloch and gaussian response
  //

  if ((track->GetStatus()&AliESDtrack::kTPCin )==0) return;

  AliAODPid *pidObj = track->GetDetPid();
  Double_t mom      = track->P();
  Double_t dedx     = 0.;  
  UShort_t nTPCClus = 0;
  if (pidObj) {
      nTPCClus = pidObj->GetTPCsignalN();
      dedx     = pidObj->GetTPCsignal();
      mom      = pidObj->GetTPCmomentum();
  }
  
  Bool_t mismatch=kTRUE;

  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    AliPID::EParticleType type=AliPID::EParticleType(j);
    Double_t bethe=fTPCResponse.GetExpectedSignal(mom,type); 
    Double_t sigma=fTPCResponse.GetExpectedSigma(mom,nTPCClus,type);
    if (TMath::Abs(dedx-bethe) > fRange*sigma) {
      p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
    } else {
      p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
      mismatch=kFALSE;
    }

  }

  if (mismatch)
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;


  return;
}
//_________________________________________________________________________
void AliAODpidUtil::MakeITSPID(const AliAODTrack *track,Double_t *p) const
{
  //
  // ITS PID
  //  1) Truncated mean method
  //


  if ((track->GetStatus()&AliESDtrack::kITSin)==0) return;
  UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
   if(clumap&(1<<i)) ++nPointsForPid;
  }
  if(nPointsForPid<3) { // track not to be used for combined PID purposes
    for (Int_t j=0; j<AliPID::kSPECIES; j++) 
      p[j] = 1./AliPID::kSPECIES;
    return;
  }
  Double_t mom=track->P();  
  AliAODPid *pidObj = track->GetDetPid();

  Double_t dedx = 0.;
  if (pidObj) {
      dedx = pidObj->GetITSsignal();
  }
  
  Bool_t mismatch = kTRUE;
  Bool_t isSA = kTRUE;
  if(track->GetStatus() & AliESDtrack::kTPCin){
    isSA = kFALSE;
    if (pidObj)
      mom = pidObj->GetTPCmomentum();
  }
  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    Double_t mass = AliPID::ParticleMass(j);//GeV/c^2
    Double_t bethe = fITSResponse.Bethe(mom,mass);
    Double_t sigma = fITSResponse.GetResolution(bethe,nPointsForPid,isSA);
    if (TMath::Abs(dedx-bethe) > fRange*sigma) {
      p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
    } else {
      p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
      mismatch=kFALSE;
    }

    // Check for particles heavier than (AliPID::kSPECIES - 1)

  }

  if (mismatch)
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;

  return;

}
//_________________________________________________________________________
void AliAODpidUtil::MakeTOFPID(const AliAODTrack *track, Double_t *p) const
{
  //
  //   TOF PID using gaussian response
  //
  if ((track->GetStatus()&AliESDtrack::kTOFout )==0)    return;
  if ((track->GetStatus()&AliESDtrack::kTIME )==0)     return;
  if ((track->GetStatus()&AliESDtrack::kTOFpid )==0)   return;

  Double_t time[AliPID::kSPECIES];
  Double_t sigma[AliPID::kSPECIES];
  AliAODPid *pidObj = track->GetDetPid();
  pidObj->GetIntegratedTimes(time);

  for (Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++) {
    sigma[iPart] = fTOFResponse.GetExpectedSigma(track->P(),time[iPart],AliPID::ParticleMass(iPart));
  }

  AliDebugGeneral("AliESDpid::MakeTOFPID",2,
		  Form("Expected TOF signals [ps]: %f %f %f %f %f",
		       time[AliPID::kElectron],
		       time[AliPID::kMuon],
		       time[AliPID::kPion],
		       time[AliPID::kKaon],
		       time[AliPID::kProton]));

  AliDebugGeneral("AliESDpid::MakeTOFPID",2,
		  Form("Expected TOF std deviations [ps]: %f %f %f %f %f",
		       sigma[AliPID::kElectron],
		       sigma[AliPID::kMuon],
		       sigma[AliPID::kPion],
		       sigma[AliPID::kKaon],
		       sigma[AliPID::kProton]
		       ));

  Double_t tof = pidObj->GetTOFsignal();

  Bool_t mismatch = kTRUE;
  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    Double_t sig = sigma[j];
    if (TMath::Abs(tof-time[j]) > fRange*sig) {
      p[j] = TMath::Exp(-0.5*fRange*fRange)/sig;
    } else
      p[j] = TMath::Exp(-0.5*(tof-time[j])*(tof-time[j])/(sig*sig))/sig;

    // Check the mismatching
    Double_t mass = AliPID::ParticleMass(j);
    Double_t pm = fTOFResponse.GetMismatchProbability(track->P(),mass);
    if (p[j]>pm) mismatch = kFALSE;

    // Check for particles heavier than (AliPID::kSPECIES - 1)

  }

  if (mismatch)
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;

  return;
}
//_________________________________________________________________________
void AliAODpidUtil::MakeTRDPID(const AliAODTrack *track,Double_t *p) const
{
  ComputeTRDProbability(track, AliPID::kSPECIES, p); 
  return;
}

//_________________________________________________________________________
Float_t AliAODpidUtil::NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Number of sigma implementation for the TOF
  //
  
  AliAODTrack *track=(AliAODTrack*)vtrack;

  // look for cached value first
  if (track->GetDetectorPID()){
    return track->GetDetectorPID()->GetNumberOfSigmas(kTOF, type);
  }  
  
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
