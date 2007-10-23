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
  fPmax(0)         // zero at 0.5 GeV/c for pp
{
}
//_________________________________________________________________________
AliTOFpidESD::AliTOFpidESD(Double_t *param):
  fSigma(param[0]),
  fRange(param[1]),
  fPmax(0)          // zero at 0.5 GeV/c for pp
{
  //
  //  The main constructor
  //
  //

  //fPmax=TMath::Exp(-0.5*3*3)/fSigma; // ~3 sigma at 0.5 GeV/c for PbPb 
}

Double_t 
AliTOFpidESD::GetMismatchProbability(Double_t p, Double_t mass) const {
  //
  // Returns the probability of mismatching 
  // assuming 1/(p*beta)^2 scaling
  //
  const Double_t m=0.5;                   // "reference" momentum (GeV/c)

  Double_t ref2=m*m*m*m/(m*m + mass*mass);// "reference" (p*beta)^2
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
    Double_t tof=t->GetTOFsignal()-timeZero;
    Double_t time[10]; t->GetIntegratedTimes(time);
    Double_t p[10];
    Double_t mom=t->GetP();
    Bool_t mismatch=kTRUE, heavy=kTRUE;
    for (Int_t j=0; j<AliPID::kSPECIES; j++) {
      Double_t mass=AliPID::ParticleMass(j);
      Double_t dpp=0.01;      //mean relative pt resolution;
      if (mom>0.5) dpp=0.01*mom;
      Double_t sigma=dpp*time[j]/(1.+ mom*mom/(mass*mass));
      sigma=TMath::Sqrt(sigma*sigma + fSigma*fSigma);
      if (TMath::Abs(tof-time[j]) > fRange*sigma) {
	p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
      } else 
        p[j]=TMath::Exp(-0.5*(tof-time[j])*(tof-time[j])/(sigma*sigma))/sigma;

      // Check the mismatching
      Double_t pm=GetMismatchProbability(mom,mass);
      if (p[j]>pm) mismatch=kFALSE;

      // Check for particles heavier than (AliPID::kSPECIES - 1)
      if (tof < (time[j] + fRange*sigma)) heavy=kFALSE;

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
    Double_t tof=t->GetTOFsignal();
    Double_t time[10]; t->GetIntegratedTimes(time);
    Double_t p[10];
    Double_t mom=t->GetP();
    Bool_t mismatch=kTRUE, heavy=kTRUE;
    for (Int_t j=0; j<AliPID::kSPECIES; j++) {
      Double_t mass=AliPID::ParticleMass(j);
      Double_t dpp=0.01;      //mean relative pt resolution;
      if (mom>0.5) dpp=0.01*mom;
      Double_t sigma=dpp*time[j]/(1.+ mom*mom/(mass*mass));
      sigma=TMath::Sqrt(sigma*sigma + fSigma*fSigma);
      if (TMath::Abs(tof-time[j]) > fRange*sigma) {
	p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
      } else
        p[j]=TMath::Exp(-0.5*(tof-time[j])*(tof-time[j])/(sigma*sigma))/sigma;

      // Check the mismatching
      Double_t pm=GetMismatchProbability(mom,mass);
      if (p[j]>pm) mismatch=kFALSE;

      // Check for particles heavier than (AliPID::kSPECIES - 1)
      if (tof < (time[j] + fRange*sigma)) heavy=kFALSE;

    }

    if (mismatch)
       for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1/AliPID::kSPECIES;

    t->SetTOFpid(p);

    if (heavy) t->ResetStatus(AliESDtrack::kTOFpid);    

  }

  delete[] tracks;
  
  return 0;
}

