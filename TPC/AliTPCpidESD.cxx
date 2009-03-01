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
//-----------------------------------------------------------------

#include "AliTPCpidESD.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMathBase.h"

ClassImp(AliTPCpidESD)

//_________________________________________________________________________
  AliTPCpidESD::AliTPCpidESD(Double_t *param):
    fMIP(0.),
    fRes(0.),
    fRange(0.)
{
  //
  //  The main constructor
  //
  fMIP=param[0];
  fRes=param[1];
  fRange=param[2];
}

Double_t AliTPCpidESD::Bethe(Double_t betaGamma) {
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
  Double_t bb = AliMathBase::BetheBlochAleph(betaGamma);
  Double_t meanCorrection =(1+(bb-1)*kmeanCorrection);
  bb *= meanCorrection;
  return bb;
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
    Double_t mom=t->GetP();
    const AliExternalTrackParam *in=t->GetInnerParam();
    if (in) mom=in->GetP();
    Double_t dedx=t->GetTPCsignal()/fMIP;
    Bool_t mismatch=kTRUE, heavy=kTRUE;
    for (Int_t j=0; j<AliPID::kSPECIES; j++) {
      Double_t mass=AliPID::ParticleMass(j);
      Double_t bethe=Bethe(mom/mass); 
      Double_t sigma=fRes*bethe;
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
