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
#include "AliESD.h"
#include "AliESDtrack.h"

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

Double_t AliTPCpidESD::Bethe(Double_t bg) {
  //
  // This is the Bethe-Bloch function normalised to 1 at the minimum
  //
  Double_t bg2=bg*bg;
  Double_t bethe;
  if (bg<3.5e1) 
      bethe=(1.+ bg2)/bg2*(log(5940*bg2) - bg2/(1.+ bg2));
  else // Density effect ( approximately :) 
      bethe=1.15*(1.+ bg2)/bg2*(log(3.5*5940*bg) - bg2/(1.+ bg2));
  return bethe/11.091;
}

//_________________________________________________________________________
Int_t AliTPCpidESD::MakePID(AliESD *event)
{
  //
  //  This function calculates the "detector response" PID probabilities 
  //
  Int_t ntrk=event->GetNumberOfTracks();
  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    if ((t->GetStatus()&AliESDtrack::kTPCin )==0)
      if ((t->GetStatus()&AliESDtrack::kTPCout)==0) continue;
    Int_t ns=AliPID::kSPECIES;
    Double_t p[10];
    for (Int_t j=0; j<ns; j++) {
      Double_t mass=AliPID::ParticleMass(j);
      Double_t mom=t->GetP();
      Double_t dedx=t->GetTPCsignal()/fMIP;
      Double_t bethe=Bethe(mom/mass); 
      Double_t sigma=fRes*bethe;
      if (TMath::Abs(dedx-bethe) > fRange*sigma) {
	p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
        continue;
      }
      p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
    }
    t->SetTPCpid(p);
  }
  return 0;
}
