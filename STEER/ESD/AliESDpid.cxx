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

/* $Id: AliESDpid.cxx 64123 2013-09-05 15:09:53Z morsch $ */

//-----------------------------------------------------------------
//           Implementation of the combined PID class
//           For the Event Summary Data Class
//           produced by the reconstruction process
//           and containing information on the particle identification
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "TArrayI.h"
#include "TArrayF.h"

#include "AliLog.h"
#include "AliPID.h"
#include "AliTOFHeader.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliTOFPIDParams.h"

#include <AliDetectorPID.h>

ClassImp(AliESDpid)

Bool_t AliESDpid::fgUseElectronExclusionBands = kFALSE;
Int_t  AliESDpid::fgNSpeciesForTracking = AliPID::kSPECIESC;

Int_t AliESDpid::MakePID(AliESDEvent *event, Bool_t TPConly, Float_t /*timeZeroTOF*/) const {
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
      //MakeTOFPID(track, timeZeroTOF);
      //MakeHMPIDPID(track);
      //MakeTRDPID(track);
    }
    CombinePID(track);
  }
  return 0;
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

    Double_t p[AliPID::kSPECIES];

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

  Double_t time[AliPID::kSPECIESC];
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
  Double_t p[AliPID::kSPECIES]={1.};

  if (track->IsOn(AliESDtrack::kITSpid)) {
    Double_t d[AliPID::kSPECIES];
    track->GetITSpid(d);
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kTPCpid)) {
    Double_t d[AliPID::kSPECIES];
    track->GetTPCpid(d);
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kTRDpid)) {
    Double_t d[AliPID::kSPECIES];
    track->GetTRDpid(d);
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kTOFpid)) {
    Double_t d[AliPID::kSPECIES];
    track->GetTOFpid(d);
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kHMPIDpid)) {
    Double_t d[AliPID::kSPECIES];
    track->GetHMPIDpid(d);
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]*=d[j];
  }

  track->SetESDpid(p);
}
//_________________________________________________________________________
Bool_t AliESDpid::CheckTOFMatching(AliESDtrack *track) const{
  //
  // Check pid matching of TOF with TPC as reference
  //
    Bool_t status = kFALSE;

    Double_t exptimes[AliPID::kSPECIESC];
    track->GetIntegratedTimes(exptimes);

    Float_t p = track->P();

    Float_t dedx = track->GetTPCsignal();
    Float_t time = track->GetTOFsignal() - fTOFResponse.GetStartTime(p);

    Double_t ptpc[3];
    track->GetInnerPxPyPz(ptpc);
    Float_t momtpc=TMath::Sqrt(ptpc[0]*ptpc[0] + ptpc[1]*ptpc[1] + ptpc[2]*ptpc[2]);

    for(Int_t i=0;i < AliPID::kSPECIES;i++){
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
  Double_t tofTime = 99999;
  if (fTuneMConData && ((fTuneMConDataMask & kDetTOF) == kDetTOF) ) tofTime = (Double_t)this->GetTOFsignalTunedOnData((AliVTrack*)vtrack);
  else tofTime=vtrack->GetTOFsignal();
  tofTime = tofTime  - fTOFResponse.GetStartTime(vtrack->P());
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
  Double_t tofTime = 99999;
  if (fTuneMConData && ((fTuneMConDataMask & kDetTOF) == kDetTOF) ) tofTime = (Double_t)this->GetTOFsignalTunedOnData((AliVTrack*)vtrack);
  else tofTime=vtrack->GetTOFsignal();

  Double_t expTime = fTOFResponse.GetExpectedSignal(vtrack,type);
  return (tofTime - fTOFResponse.GetStartTime(vtrack->P()) - expTime)/fTOFResponse.GetExpectedSigma(vtrack->P(),expTime,AliPID::ParticleMassZ(type));
}

//_________________________________________________________________________
void AliESDpid::SetPIDForTracking(AliESDtrack *esdtr) const
{
  // assign mass for tracking
  //

  // in principle AliPIDCombined could be used to also set priors
  //AliPIDCombined pidProb;
  //pidProb.SetDetectorMask(kDetTPC);              // use only TPC information, couls also be changed
  //pidProb.SetSelectedSpecies(AliPID::kSPECIESC); // number of species to use
  //pidProb.SetDefaultTPCPriors();                 // load default priors

  Double_t prob[AliPID::kSPECIESC]={0.};
  EDetPidStatus pidStatus=ComputePIDProbability(kTPC, esdtr, AliPID::kSPECIESC, prob);
  // check if a valid signal was found, otherwise return pion mass
  if (pidStatus==kDetNoSignal) { //kDetPidOk) {
    esdtr->SetPIDForTracking(AliPID::kPion);
    return;
  }

  // or with AliPIDCombined
  // pidProb.ComputeProbabilities(esdtr, this, p);

  // find max probability
  Float_t max=0.,min=1.e9;
  Int_t pid=-1;
  //  for (Int_t i=0; i<AliPID::kSPECIESC; ++i) {
  for (Int_t i=0; i<fgNSpeciesForTracking; ++i) {  // ignore nuclei
    if (prob[i]>max) {pid=i; max=prob[i];}
    if (prob[i]<min) min=prob[i];
  }

  //int pid = AliPID::kPion; // this should be substituted by real most probable TPC pid (e,mu -> pion) or poin if no PID possible

  //
  if (pid>AliPID::kSPECIESC-1 || (min>=max)) pid = AliPID::kPion;
  //
  if (pid==0 && fgUseElectronExclusionBands) { // dE/dx "crossing points" in the TPC 
    Double_t p = esdtr->GetP();
    if ((p>0.38)&&(p<0.48)) {
      if (prob[0]<prob[3]*10.) pid = AliPID::kKaon;
    }
    else if ((p>0.75)&&(p<0.85)) {
      if (prob[0]<prob[4]*10.) pid = AliPID::kProton;
    } 
  }

  esdtr->SetPIDForTracking( pid );
  //
}


//_________________________________________________________________________
void AliESDpid::MakePIDForTracking(AliESDEvent *event) const
{
  // assign masses using for tracking
  Int_t nTrk=event->GetNumberOfTracks();
  for (Int_t iTrk=nTrk; iTrk--;) {
    AliESDtrack *track = event->GetTrack(iTrk);
    SetPIDForTracking(track);
  }
}

//_________________________________________________________________________
void AliESDpid::SetUseElectronExclusionBands(Bool_t val)
{
  // disable electron tag in e/K and e/p crossing (a la Run1)
  fgUseElectronExclusionBands=val;
  AliInfoClassF("Exclude electron tag in e/K and e/p crossings: %d",fgUseElectronExclusionBands);
}

//_________________________________________________________________________
void AliESDpid::SetNSpeciesForTracking(Int_t n)
{
  // set max number of species to consider for tracking
  fgNSpeciesForTracking = n>0 ? n : AliPID::kSPECIESC;
  AliInfoClassF("First %d species will be considered for tracking PID",fgNSpeciesForTracking);
}
