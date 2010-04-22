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

#include "AliLog.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"

ClassImp(AliESDpid)

Int_t AliESDpid::MakePID(AliESDEvent *event, Bool_t TPConly, Float_t TimeZeroTOF) const {
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
  Int_t nTrk=event->GetNumberOfTracks();
  for (Int_t iTrk=0; iTrk<nTrk; iTrk++) {  
    AliESDtrack *track=event->GetTrack(iTrk);
    MakeTPCPID(track);
    if (!TPConly) {
      MakeITSPID(track);
      MakeTOFPID(track, TimeZeroTOF);
      //MakeHMPIDPID(track);
      MakeTRDPID(track);
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
       for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1/AliPID::kSPECIES;

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
    Double_t p[10];
    Bool_t mismatch=kTRUE, heavy=kTRUE;
    for (Int_t j=0; j<AliPID::kSPECIES; j++) {
      Double_t mass=AliPID::ParticleMass(j);//GeV/c^2
      Double_t bethe=fITSResponse.Bethe(mom,mass);
      Double_t sigma=fITSResponse.GetResolution(bethe);
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
void AliESDpid::MakeTOFPID(AliESDtrack *track, Float_t TimeZeroTOF) const
{
  //
  //   TOF PID using gaussian response
  //
  if ((track->GetStatus()&AliESDtrack::kTOFout)==0) return;
  if ((track->GetStatus()&AliESDtrack::kTIME)==0) return;

  Double_t time[AliPID::kSPECIESN];
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

  Double_t tof = track->GetTOFsignal() - TimeZeroTOF;

  Double_t p[AliPID::kSPECIES];
  Bool_t mismatch = kTRUE, heavy = kTRUE;
  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    Double_t sig = sigma[j];
    if (TMath::Abs(tof-time[j]) > fRange*sig) {
      p[j] = TMath::Exp(-0.5*fRange*fRange)/sig;
    } else
      p[j] = TMath::Exp(-0.5*(tof-time[j])*(tof-time[j])/(sig*sig))/sig;

    // Check the mismatching
    Double_t mass = AliPID::ParticleMass(j);
    Double_t pm = fTOFResponse.GetMismatchProbability(track->GetP(),mass);
    if (p[j]>pm) mismatch = kFALSE;

    // Check for particles heavier than (AliPID::kSPECIES - 1)
    if (tof < (time[j] + fRange*sig)) heavy=kFALSE;

  }

  if (mismatch)
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1/AliPID::kSPECIES;

  track->SetTOFpid(p);

  if (heavy) track->ResetStatus(AliESDtrack::kTOFpid);    
}
//_________________________________________________________________________
void AliESDpid::MakeTRDPID(AliESDtrack *track) const
{
  //
  // Method to recalculate the TRD PID probabilities
  //
  if((track->GetStatus()&AliESDtrack::kTRDout)==0) return;
  Double_t prob[AliPID::kSPECIES]; Float_t mom[6];
  Double_t dedx[48];  // Allocate space for the maximum number of TRD slices
  for(Int_t ilayer = 0; ilayer < 6; ilayer++){
    mom[ilayer] = track->GetTRDmomentum(ilayer);
    for(Int_t islice = 0; islice < track->GetNumberOfTRDslices(); islice++){
      dedx[ilayer*track->GetNumberOfTRDslices()+islice] = track->GetTRDslice(ilayer, islice);
    }
  }
  fTRDResponse.GetResponse(track->GetNumberOfTRDslices(), dedx, mom, prob);
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
