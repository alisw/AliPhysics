/**************************************************************************
 * Copyright(c) 2005-2007, ALICE Experiment at CERN, All rights reserved. *
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
// ITS PID method # 1
//           Implementation of the ITS PID class
// Very naive one... Should be made better by the detector experts...
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------
#include "AliITSpidESD.h"
#include "AliITSpidESD1.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"

ClassImp(AliITSpidESD1)

AliITSpidESD1::AliITSpidESD1(): AliITSpidESD(),
fMIP(0),
fRes(0),
fRange(0) 
{
  //Default constructor
}
//_________________________________________________________________________
AliITSpidESD1::AliITSpidESD1(Double_t *param): AliITSpidESD(),
fMIP(param[0]),
fRes(param[1]),
fRange(param[2])
{
  //
  //  The main constructor
  //
}


//_________________________________________________________________________
Int_t AliITSpidESD1::MakePID(AliESDEvent *event)
{
  //
  //  This function calculates the "detector response" PID probabilities 
  //
  Int_t ntrk=event->GetNumberOfTracks();
  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    if ((t->GetStatus()&AliESDtrack::kITSin )==0)
      if ((t->GetStatus()&AliESDtrack::kITSout)==0) continue;
    Double_t mom=t->GetP();
    Double_t dedx=t->GetITSsignal()/fMIP;
    Double_t p[10];
    Bool_t mismatch=kTRUE, heavy=kTRUE;
    for (Int_t j=0; j<AliPID::kSPECIES; j++) {
      Double_t mass=AliPID::ParticleMass(j);//GeV/c^2
      Double_t bethe=AliITSpidESD::Bethe(mom,mass);
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

    t->SetITSpid(p);

    if (heavy) t->ResetStatus(AliESDtrack::kITSpid);

  }
  return 0;
}
