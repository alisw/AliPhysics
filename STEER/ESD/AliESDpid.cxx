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

/* $Id$ */

//-----------------------------------------------------------------
//           Implementation of the combined PID class
//           For the Event Summary Data Class
//           produced by the reconstruction process
//           and containing information on the particle identification
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "TArrayI.h"
#include "TArrayF.h"

#include "TRandom.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliTOFHeader.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"

#include <AliDetectorPID.h>

ClassImp(AliESDpid)

Int_t AliESDpid::MakePID(AliESDEvent *event, Bool_t TPConly, Float_t timeZeroTOF) const {
  //
  //  Calculate probabilities for all detectors, except if TPConly==kTRUE
  //  and combine PID
  //  
  //   Option TPConly==kTRUE is used during reconstruction, 
  //  because ITS tracking uses TPC pid
  //  HMPID and TRD pid are done in detector reconstructors
  //

  /*
  Float_t timeZeroTOF = 0;
  if (subtractT0) 
    timeZeroTOF = event->GetT0();
  */
  Int_t nTrk=event->GetNumberOfTracks();
  for (Int_t iTrk=0; iTrk<nTrk; iTrk++) {  
    AliESDtrack *track=event->GetTrack(iTrk);
    MakeTPCPID(track);
    if (!TPConly) {
      MakeITSPID(track);
      MakeTOFPID(track, timeZeroTOF);
      //MakeHMPIDPID(track);
      //MakeTRDPID(track);
    }
    CombinePID(track);
  }
  return 0;
}
//_________________________________________________________________________
Float_t AliESDpid::GetTPCsignalTunedOnData(const AliVTrack *t) const {
    AliESDtrack *track = (AliESDtrack *) t;
    Float_t dedx = track->GetTPCsignalTunedOnData();

    if(dedx > 0) return dedx;

    dedx = t->GetTPCsignal();
    track->SetTPCsignalTunedOnData(dedx);
    if(dedx < 20) return dedx;

    AliPID::EParticleType type = AliPID::kPion;

    AliMCEventHandler* eventHandler=dynamic_cast<AliMCEventHandler*>(fEventHandler);
    if (eventHandler) {
	AliMCEvent* mcEvent = eventHandler->MCEvent();
	if(mcEvent){
	    Bool_t kGood = kTRUE;
	    AliMCParticle *MCpart = (AliMCParticle *) mcEvent->GetTrack(TMath::Abs(t->GetLabel()));
	    TParticle *part = MCpart->Particle();
	    
	    Int_t iS = TMath::Abs(part->GetPdgCode());

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
// 		if(iS == AliPID::ParticleCode(AliPID::kHe3) || iS == AliPID::ParticleCode(AliPID::kAlpha)) dedx *= 5;
	    }
	}
    }

    track->SetTPCsignalTunedOnData(dedx);
    return dedx;
}
//_________________________________________________________________________
void AliESDpid::MakeTPCPID(AliESDtrack *track) const
{
  //
  //  TPC pid using bethe-bloch and gaussian response
  //
  if ((track->GetStatus()&AliESDtrack::kTPCin )==0)
    if ((track->GetStatus()&AliESDtrack::kTPCout)==0) return;

    Double_t mom = track->GetP();
    const AliExternalTrackParam *in=track->GetInnerParam();
    if (in) mom = in->GetP();

    Double_t p[AliPID::kSPECIES];
    Double_t dedx=track->GetTPCsignal(); 
    Bool_t mismatch=kTRUE, heavy=kTRUE;

    for (Int_t j=0; j<AliPID::kSPECIES; j++) {
      AliPID::EParticleType type=AliPID::EParticleType(j);
      Double_t bethe=fTPCResponse.GetExpectedSignal(mom,type); 
      Double_t sigma=fTPCResponse.GetExpectedSigma(mom,track->GetTPCsignalN(),type);
      if (TMath::Abs(dedx-bethe) > fRange*sigma) {
	p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
      } else {
        p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
        mismatch=kFALSE;
      }

      // Check for particles heavier than (AliPID::kSPECIES - 1)
      if (dedx < (bethe + fRange*sigma)) heavy=kFALSE;

    }

    if (mismatch)
       for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;

    track->SetTPCpid(p);

    if (heavy) track->ResetStatus(AliESDtrack::kTPCpid);

}
//_________________________________________________________________________
void AliESDpid::MakeITSPID(AliESDtrack *track) const
{
  //
  // ITS PID
  // Two options, depending on fITSPIDmethod:
  //  1) Truncated mean method
  //  2) Likelihood, using charges measured in all 4 layers and 
  //     Landau+gaus response functions
  //

  if ((track->GetStatus()&AliESDtrack::kITSin)==0 &&
      (track->GetStatus()&AliESDtrack::kITSout)==0) return;

  Double_t mom=track->GetP();  
  if (fITSPIDmethod == kITSTruncMean) {
    Double_t dedx=track->GetITSsignal();
    Bool_t isSA=kTRUE;
    Double_t momITS=mom;
    ULong_t trStatus=track->GetStatus();
    if(trStatus&AliESDtrack::kTPCin) isSA=kFALSE;
    UChar_t clumap=track->GetITSClusterMap();
    Int_t nPointsForPid=0;
    for(Int_t i=2; i<6; i++){
      if(clumap&(1<<i)) ++nPointsForPid;
    }

    if(nPointsForPid<3) { // track not to be used for combined PID purposes
      track->ResetStatus(AliESDtrack::kITSpid);
      return;
    }

    Double_t p[10];

    Bool_t mismatch=kTRUE, heavy=kTRUE;
    for (Int_t j=0; j<AliPID::kSPECIES; j++) {
      Double_t mass=AliPID::ParticleMass(j);//GeV/c^2
      Double_t bethe=fITSResponse.Bethe(momITS,mass);
      Double_t sigma=fITSResponse.GetResolution(bethe,nPointsForPid,isSA);
      if (TMath::Abs(dedx-bethe) > fRange*sigma) {
	p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
      } else {
        p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
        mismatch=kFALSE;
      }

      // Check for particles heavier than (AliPID::kSPECIES - 1)
      if (dedx < (bethe + fRange*sigma)) heavy=kFALSE;

    }

    if (mismatch)
       for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;

    track->SetITSpid(p);

    if (heavy) track->ResetStatus(AliESDtrack::kITSpid);
  }
  else {  // Likelihood method
    Double_t condprobfun[AliPID::kSPECIES];
    Double_t qclu[4];
    track->GetITSdEdxSamples(qclu);
    fITSResponse.GetITSProbabilities(mom,qclu,condprobfun);
    track->SetITSpid(condprobfun);
  }

}
//_________________________________________________________________________
void AliESDpid::MakeTOFPID(AliESDtrack *track, Float_t /*timeZeroTOF*/) const
{
  //
  //   TOF PID using gaussian response
  //

  if ((track->GetStatus()&AliESDtrack::kTOFout)==0) return;
  if ((track->GetStatus()&AliESDtrack::kTIME)==0) return;
  if ((track->GetStatus()&AliESDtrack::kITSin)==0) return;

  Int_t ibin = fTOFResponse.GetMomBin(track->GetP());
  Float_t timezero = fTOFResponse.GetT0bin(ibin);

  Double_t time[AliPID::kSPECIES];
  track->GetIntegratedTimes(time);

  Double_t sigma[AliPID::kSPECIES];
  for (Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++) {
    sigma[iPart] = fTOFResponse.GetExpectedSigma(track->GetP(),time[iPart],AliPID::ParticleMass(iPart));
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

  Double_t tof = track->GetTOFsignal() - timezero;

  Double_t p[AliPID::kSPECIES];
//   Bool_t mismatch = kTRUE;
  Bool_t heavy = kTRUE;
  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    Double_t sig = sigma[j];
    if (TMath::Abs(tof-time[j]) > (fRange+2)*sig) {
	p[j] = TMath::Exp(-0.5*(fRange+2)*(fRange+2))/sig;
    } else
      p[j] = TMath::Exp(-0.5*(tof-time[j])*(tof-time[j])/(sig*sig))/sig;

    // Check the mismatching
    
//     Double_t mass = AliPID::ParticleMass(j);
//     Double_t pm = fTOFResponse.GetMismatchProbability(track->GetP(),mass);
//     if (p[j]>pm) mismatch = kFALSE;

    // Check for particles heavier than (AliPID::kSPECIES - 1)
    if (tof < (time[j] + fRange*sig)) heavy=kFALSE;

  }

  /*
    if (mismatch)
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;
  */

  track->SetTOFpid(p);

  if (heavy) track->ResetStatus(AliESDtrack::kTOFpid);
  
  // kTOFmismatch flas is not set because deprecated from 18/02/2013
  //  if (!CheckTOFMatching(track)) track->SetStatus(AliESDtrack::kTOFmismatch);
  //  else track->ResetStatus(AliESDtrack::kTOFmismatch);
}
//_________________________________________________________________________
void AliESDpid::MakeTRDPID(AliESDtrack *track) const
{
  //
  // Method to recalculate the TRD PID probabilities
  //
  Double_t prob[AliPID::kSPECIES];
  GetComputeTRDProbability(track, AliPID::kSPECIES, prob);
  track->SetTRDpid(prob);
}
//_________________________________________________________________________
void AliESDpid::CombinePID(AliESDtrack *track) const
{
  //
  // Combine the information of various detectors
  // to determine the Particle Identification
  //
  Int_t ns=AliPID::kSPECIES;
  Double_t p[10]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

  if (track->IsOn(AliESDtrack::kITSpid)) {
    Double_t d[10];
    track->GetITSpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kTPCpid)) {
    Double_t d[10];
    track->GetTPCpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kTRDpid)) {
    Double_t d[10];
    track->GetTRDpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kTOFpid)) {
    Double_t d[10];
    track->GetTOFpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kHMPIDpid)) {
    Double_t d[10];
    track->GetHMPIDpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  track->SetESDpid(p);
}
//_________________________________________________________________________
Bool_t AliESDpid::CheckTOFMatching(AliESDtrack *track) const{
  //
  // Check pid matching of TOF with TPC as reference
  //
    Bool_t status = kFALSE;
    
    Double_t exptimes[5];
    track->GetIntegratedTimes(exptimes);
    
    Float_t p = track->P();
    
    Float_t dedx = track->GetTPCsignal();
    Float_t time = track->GetTOFsignal() - fTOFResponse.GetStartTime(p);
    
    Double_t ptpc[3];
    track->GetInnerPxPyPz(ptpc);
    Float_t momtpc=TMath::Sqrt(ptpc[0]*ptpc[0] + ptpc[1]*ptpc[1] + ptpc[2]*ptpc[2]);
    
    for(Int_t i=0;i < 5;i++){
	AliPID::EParticleType type=AliPID::EParticleType(i);
	
	Float_t resolutionTOF = fTOFResponse.GetExpectedSigma(p, exptimes[i], AliPID::ParticleMass(i));
	if(TMath::Abs(exptimes[i] - time) < fRange * resolutionTOF){
	    Float_t dedxExp = fTPCResponse.GetExpectedSignal(momtpc,type);
	    Float_t resolutionTPC = fTPCResponse.GetExpectedSigma(momtpc,track->GetTPCsignalN(),type);
	    
	    if(TMath::Abs(dedx - dedxExp) < fRangeTOFMismatch * resolutionTPC){
		status = kTRUE;
	    }
	}
    }
    
    // for nuclei
    Float_t resolutionTOFpr = fTOFResponse.GetExpectedSigma(p, exptimes[4], AliPID::ParticleMass(4));
    if(!status && (exptimes[4] + fRange*resolutionTOFpr < time)) status = kTRUE;
    
    
    return status;
}

//_________________________________________________________________________
Float_t AliESDpid::GetSignalDeltaTOFold(const AliVParticle *track, AliPID::EParticleType type, Bool_t ratio/*=kFALSE*/) const
{
  //
  // TOF signal - expected
  //
  AliVTrack *vtrack=(AliVTrack*)track;
  
  const Double_t expTime = fTOFResponse.GetExpectedSignal(vtrack,type);
  const Double_t tofTime=vtrack->GetTOFsignal() - fTOFResponse.GetStartTime(vtrack->P());
  Double_t delta=-9999.;

  if (!ratio) delta=tofTime-expTime;
  else if (expTime>1.e-20) delta=tofTime/expTime;
  
  return delta;
}

//_________________________________________________________________________
Float_t AliESDpid::GetNumberOfSigmasTOFold(const AliVParticle *track, AliPID::EParticleType type) const
{
  //
  // Number of sigma implementation for the TOF
  //

  AliVTrack *vtrack=(AliVTrack*)track;
  
  Double_t expTime = fTOFResponse.GetExpectedSignal(vtrack,type);
  return (vtrack->GetTOFsignal() - fTOFResponse.GetStartTime(vtrack->P()) - expTime)/fTOFResponse.GetExpectedSigma(vtrack->P(),expTime,AliPID::ParticleMassZ(type));
}
