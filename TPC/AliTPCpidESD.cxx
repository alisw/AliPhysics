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

//-----------------------------------------------------------------
//           Implementation of the TPC PID class
// Very naive one... Should be made better by the detector experts...
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
// With many additions and modifications suggested by
//      Alexander Kalweit, GSI, alexander.philipp.kalweit@cern.ch
//      Dariusz Miskowiec, GSI, D.Miskowiec@gsi.de
//-----------------------------------------------------------------

#include "AliTPCpidESD.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMathBase.h"

ClassImp(AliTPCpidESD)

//_________________________________________________________________________
AliTPCpidESD::AliTPCpidESD():
    fMIP(50.),
    fRes(0.07),
    fRange(5.),
    fKp1(0.76176e-1),
    fKp2(10.632),
    fKp3(0.13279e-4),
    fKp4(1.8631),
    fKp5(1.9479)
{
  //
  //  The default constructor
  //
}

//_________________________________________________________________________
AliTPCpidESD::AliTPCpidESD(Double_t *param):
    fMIP(param[0]),
    fRes(param[1]),
    fRange(param[2]),
    fKp1(0.76176e-1),
    fKp2(10.632),
    fKp3(0.13279e-4),
    fKp4(1.8631),
    fKp5(1.9479)
{
  //
  //  The main constructor
  //
}

Double_t AliTPCpidESD::Bethe(Double_t betaGamma) const {
  //
  // This is the Bethe-Bloch function normalised to 1 at the minimum
  // WARNING
  // Simulated and reconstructed Bethe-Bloch differs
  //           Simulated  curve is the dNprim/dx
  //           Reconstructed is proportianal dNtot/dx
  // Temporary fix for production -  Simple linear correction function
  // Future    2 Bethe Bloch formulas needed
  //           1. for simulation
  //           2. for reconstructed PID
  //
  const Float_t kmeanCorrection =0.1;
  Double_t bb=AliMathBase::BetheBlochAleph(betaGamma,fKp1,fKp2,fKp3,fKp4,fKp5);
  Double_t meanCorrection =(1+(bb-1)*kmeanCorrection);
  bb *= meanCorrection;
  return bb;
}

//_________________________________________________________________________
void AliTPCpidESD::SetBetheBlochParameters(Double_t kp1,
                             Double_t kp2,
                             Double_t kp3,
                             Double_t kp4,
                             Double_t kp5) {
  //
  // Set the parameters of the ALEPH Bethe-Bloch formula
  //
  fKp1=kp1;
  fKp2=kp2;
  fKp3=kp3;
  fKp4=kp4;
  fKp5=kp5;
}

//_________________________________________________________________________
Bool_t AliTPCpidESD::ExpectedSigmas(const AliESDtrack *t,
                                     Double_t s[],
                                     Int_t n) const {
  //
  // Calculate the expected dE/dx resolution as the function of 
  // the information stored in the track.
  //
  // At the moment, this resolution is just proportional to the expected 
  // signal. This can be improved. By taking into account the number of
  // assigned clusters, for example.  
  //
  Bool_t ok=kFALSE;
  Double_t signals[AliPID::kSPECIESN]; 
  if (ExpectedSignals(t,signals,n)) {
     for (Int_t i=0; i<n; i++) s[i] = fRes*signals[i];
     ok=kTRUE;
  }
  return ok;
}

//_________________________________________________________________________
Bool_t AliTPCpidESD::NumberOfSigmas(const AliESDtrack *t,
                                     Double_t s[],
                                     Int_t n) const {
  //
  // Calculate the deviation of the actual PID signal from the expected
  // signal, in units of expected sigmas.
  //
  Bool_t ok=kFALSE;

  Double_t dedx=t->GetTPCsignal()/fMIP;
  Double_t sigmas[AliPID::kSPECIESN];
  if (ExpectedSigmas(t,sigmas,n)) {
     Double_t signals[AliPID::kSPECIESN];
     if (ExpectedSignals(t,signals,n)) {
       for (Int_t i=0; i<n; i++) s[i] = (signals[i] - dedx)/sigmas[i];
       ok=kTRUE;
     }
  }
  return ok;
}

//_________________________________________________________________________
Bool_t AliTPCpidESD::ExpectedSignals(const AliESDtrack *t, 
                                      Double_t s[],
                                      Int_t n) const {
  //
  // Calculates the expected PID signals as the function of 
  // the information stored in the track.
  //  
  // At the moment, these signals are just the results of calling the 
  // Bethe-Bloch formula. 
  // This can be improved. By taking into account the number of
  // assigned clusters and/or the track dip angle, for example.  
  //
  
  Double_t mom=t->GetP();
  const AliExternalTrackParam *in=t->GetInnerParam();
  if (in) mom=in->GetP();

  for (Int_t i=0; i<n; i++) {
      Double_t mass=AliPID::ParticleMass(i); 
      s[i]=Bethe(mom/mass);
  }

  return kTRUE;
}

//_________________________________________________________________________
Double_t AliTPCpidESD::GetExpectedSignal(const AliESDtrack *t,
                                         AliPID::EParticleType n) const {
  //
  // Calculates the expected PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  // At the moment, these signals are just the results of calling the 
  // Bethe-Bloch formula. 
  // This can be improved. By taking into account the number of
  // assigned clusters and/or the track dip angle, for example.  
  //
  
  Double_t mom=t->GetP();
  const AliExternalTrackParam *in=t->GetInnerParam();
  if (in) mom=in->GetP();

  Double_t mass=AliPID::ParticleMass(n); 
  return Bethe(mom/mass);
}

//_________________________________________________________________________
Double_t AliTPCpidESD::GetExpectedSigma(const AliESDtrack *t,
                                         AliPID::EParticleType n) const {
  //
  // Calculates the expected sigma of the PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  //
  // At the moment, this sigma is just proportional to the expected 
  // signal. This can be improved. By taking into account the number of
  // assigned clusters, for example.  
  //
  
  return fRes*GetExpectedSignal(t,n);
}

//_________________________________________________________________________
Double_t AliTPCpidESD::GetNumberOfSigmas(const AliESDtrack *t,
                                         AliPID::EParticleType n) const {
  // 
  // Calculate the deviation of the actual PID signal from the expected
  // signal, in units of expected sigmas, for the specified particle type
  //
  
    Double_t dedx=t->GetTPCsignal()/fMIP;
    return (dedx - GetExpectedSignal(t,n))/GetExpectedSigma(t,n);
}

//_________________________________________________________________________
Int_t AliTPCpidESD::MakePID(AliESDEvent *event)
{
  //
  //  This function calculates the "detector response" PID probabilities 
  //
  Int_t ntrk=event->GetNumberOfTracks();
  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    if ((t->GetStatus()&AliESDtrack::kTPCin )==0)
      if ((t->GetStatus()&AliESDtrack::kTPCout)==0) continue;
    Double_t p[10];
    Double_t dedx=t->GetTPCsignal()/fMIP;
    Bool_t mismatch=kTRUE, heavy=kTRUE;
    for (Int_t j=0; j<AliPID::kSPECIES; j++) {
      AliPID::EParticleType type=AliPID::EParticleType(j);
      Double_t bethe=GetExpectedSignal(t,type); 
      Double_t sigma=GetExpectedSigma(t,type);
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

    t->SetTPCpid(p);

    if (heavy) t->ResetStatus(AliESDtrack::kTPCpid);

  }
  return 0;
}
