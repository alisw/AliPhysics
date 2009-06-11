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

//-----------------------------------------------------------------//
//                                                                 //
//           Implementation of the TOF PID class                   //
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch         //
//                                                                 //
//-----------------------------------------------------------------//

#include "TMath.h"
#include "AliLog.h"

#include "AliESDtrack.h"
#include "AliESDEvent.h"

#include "AliTOFpidESD.h"

ClassImp(AliTOFpidESD)

//_________________________________________________________________________
AliTOFpidESD::AliTOFpidESD(): 
  fSigma(0),
  fRange(0),
  fPmax(0),         // zero at 0.5 GeV/c for pp
  fTime0(0)
{
}
//_________________________________________________________________________
AliTOFpidESD::AliTOFpidESD(Double_t *param):
  fSigma(param[0]),
  fRange(param[1]),
  fPmax(0),          // zero at 0.5 GeV/c for pp
  fTime0(0)
{
  //
  //  The main constructor
  //
  //

  //fPmax=TMath::Exp(-0.5*3*3)/fSigma; // ~3 sigma at 0.5 GeV/c for PbPb 
}

//_________________________________________________________________________
Double_t 
AliTOFpidESD::GetMismatchProbability(Double_t p, Double_t mass) const {
  //
  // Returns the probability of mismatching 
  // assuming 1/(p*beta)^2 scaling
  //
  const Double_t km=0.5;                   // "reference" momentum (GeV/c)

  Double_t ref2=km*km*km*km/(km*km + mass*mass);// "reference" (p*beta)^2
  Double_t p2beta2=p*p*p*p/(p*p + mass*mass);

  return fPmax*ref2/p2beta2;
}

//_________________________________________________________________________
Int_t AliTOFpidESD::MakePID(AliESDEvent *event, Double_t timeZero)
{
  //
  //  This function calculates the "detector response" PID probabilities
  //                Just for a bare hint... 

  AliDebug(1,Form("TOF PID Parameters: Sigma (ps)= %f, Range= %f",fSigma,fRange));
  AliDebug(1,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n");

  fTime0=timeZero;
  return MakePID(event);  
}

//_________________________________________________________________________
Int_t AliTOFpidESD::MakePID(AliESDEvent *event)
{
  //
  //  This function calculates the "detector response" PID probabilities
  //                Just for a bare hint... 

  Int_t ntrk=event->GetNumberOfTracks();
  AliESDtrack **tracks=new AliESDtrack*[ntrk];

  Int_t i;
  for (i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    tracks[i]=t;
  }

  for (i=0; i<ntrk; i++) {
    AliESDtrack *t=tracks[i];
    if ((t->GetStatus()&AliESDtrack::kTOFout)==0) continue;
    if ((t->GetStatus()&AliESDtrack::kTIME)==0) continue;

    Double_t time[AliPID::kSPECIES];
    if (!ExpectedSignals(t,time,AliPID::kSPECIES)) continue;
    Double_t sigma[AliPID::kSPECIES];
    if (!ExpectedSigmas(t,sigma,AliPID::kSPECIES)) continue;

    AliDebug(2,Form("Expected TOF signals [ps]: %f %f %f %f %f",
		    GetExpectedSignal(t,AliPID::kElectron),
		    GetExpectedSignal(t,AliPID::kMuon),
		    GetExpectedSignal(t,AliPID::kPion),
		    GetExpectedSignal(t,AliPID::kKaon),
		    GetExpectedSignal(t,AliPID::kProton)
		    ));

    AliDebug(2,Form("Expected TOF std deviations [ps]: %f %f %f %f %f",
		    GetExpectedSigma(t,AliPID::kElectron),
		    GetExpectedSigma(t,AliPID::kMuon),
		    GetExpectedSigma(t,AliPID::kPion),
		    GetExpectedSigma(t,AliPID::kKaon),
		    GetExpectedSigma(t,AliPID::kProton)
		    ));

    AliDebug(2,Form("Expected TOF std deviations [number of expected sigmas]: %f %f %f %f %f",
		    GetNumberOfSigmas(t,AliPID::kElectron),
		    GetNumberOfSigmas(t,AliPID::kMuon),
		    GetNumberOfSigmas(t,AliPID::kPion),
		    GetNumberOfSigmas(t,AliPID::kKaon),
		    GetNumberOfSigmas(t,AliPID::kProton)
		    ));

    Double_t tof = t->GetTOFsignal() - fTime0;

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
      Double_t pm = GetMismatchProbability(t->GetP(),mass);
      if (p[j]>pm) mismatch = kFALSE;

      // Check for particles heavier than (AliPID::kSPECIES - 1)
      if (tof < (time[j] + fRange*sig)) heavy=kFALSE;

    }

    if (mismatch)
       for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1/AliPID::kSPECIES;

    t->SetTOFpid(p);

    if (heavy) t->ResetStatus(AliESDtrack::kTOFpid);    

  }

  delete[] tracks;
  
  return 0;
}

//_________________________________________________________________________
Bool_t AliTOFpidESD::ExpectedSignals(const AliESDtrack *t,
                                     Double_t s[],Int_t n) const
{
  //
  // Return the expected PID signals for the involved particle species
  //

  if (n > AliPID::kSPECIESN)          return kFALSE;
  if ( !t->IsOn(AliESDtrack::kTIME) ) return kFALSE;

  Double_t time[AliPID::kSPECIESN];
  t->GetIntegratedTimes(time);
  for (Int_t i=0; i<n; i++) s[i]=time[i];
  return kTRUE;
}

//_________________________________________________________________________
Bool_t AliTOFpidESD::ExpectedSigmas(const AliESDtrack *t,
                                     Double_t s[],Int_t n) const
{
  //
  // Return the expected sigma of PID signals for the involved
  // particle species.
  // This approximate (but reasonable) formula takes into account the
  // relative momentum resolution.
  //

  Double_t time[AliPID::kSPECIESN];
  if ( !ExpectedSignals(t,time,n) ) return kFALSE;

  Double_t mom = t->GetP();
  Double_t dpp = 0.01;      //mean relative pt resolution;
  if (mom>0.5) dpp = 0.01*mom;
  for (Int_t i=0; i<n; i++) {
    Double_t mass = AliPID::ParticleMass(i);
    Double_t sigma = dpp*time[i]/(1.+ mom*mom/(mass*mass));
    s[i] = TMath::Sqrt(sigma*sigma + fSigma*fSigma);
  }
  return kTRUE;  
}

//_________________________________________________________________________
Bool_t AliTOFpidESD::NumberOfSigmas(const AliESDtrack *t,
                                    Double_t s[],Int_t n) const
{
  //
  // Returns the deviation of the actual PID signal from the expected
  // signal, in units of expected sigmas.
  //

  Double_t time[AliPID::kSPECIESN];
  if ( !ExpectedSignals(t,time,n) ) return kFALSE;

  if ( !ExpectedSigmas(t,s,n) ) return kFALSE;
  
  Double_t tof = t->GetTOFsignal() - fTime0;
  for (Int_t i=0; i<n; i++) s[i] = (time[i]-tof)/s[i];
      
  return kTRUE;
}

//_________________________________________________________________________
 Double_t AliTOFpidESD::GetExpectedSignal(const AliESDtrack *t,
                                          AliPID::EParticleType n) const
{
  //
  // Return the expected PID signal for the specified particle type.
  // If the operation is not possible, return a negative value.
  //

  if (Int_t(n) >= AliPID::kSPECIESN)        return -1.;
  if ( !t->IsOn(AliESDtrack::kTIME) ) return -1.;

  Double_t time[AliPID::kSPECIESN];
  t->GetIntegratedTimes(time);

  return time[n];
}

//_________________________________________________________________________
 Double_t AliTOFpidESD::GetExpectedSigma(const AliESDtrack *t,
                                         AliPID::EParticleType n) const
{
  //
  // Return the expected sigma of the PID signal for the specified
  // particle type.
  // If the operation is not possible, return a negative value.
  //

  Double_t time[AliPID::kSPECIESN];
  if ( !ExpectedSignals(t,time,AliPID::kSPECIESN) ) return -1.;

  Double_t mom = t->GetP();
  Double_t dpp = 0.01;      //mean relative pt resolution;
  if (mom>0.5) dpp = 0.01*mom;

  Double_t mass = AliPID::ParticleMass(n);
  Double_t sigma = dpp*time[n]/(1.+ mom*mom/(mass*mass));

  return TMath::Sqrt(sigma*sigma + fSigma*fSigma);
}

//_________________________________________________________________________
 Double_t AliTOFpidESD::GetNumberOfSigmas(const AliESDtrack *t,
                                          AliPID::EParticleType n) const
{
  //
  // Returns the deviation of the actual PID signal from the expected
  // signal for the specified particle type, in units of expected
  // sigmas.
  // If the operation is not possible, return a negative value.
  //

  Double_t time=GetExpectedSignal(t,n);;
  if (time < 0.) return -1.;

  Double_t sigma=GetExpectedSigma(t,n);
  if (sigma < 0.) return -1;
  
  Double_t tof=t->GetTOFsignal() - fTime0;
  return (time-tof)/sigma;
}
