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
//           Implementation of the TRD PID class
// Very naive one... And the implementation is even poorer... 
// Should be made better by the detector experts...
//-----------------------------------------------------------------

#include "AliTRDpidESD.h"
#include "AliESD.h"
#include "AliESDtrack.h"

ClassImp(AliTRDpidESD)

//_________________________________________________________________________
AliTRDpidESD::AliTRDpidESD(Double_t *param)
{
  //
  //  The main constructor
  //
  fMIP=param[0];   // MIP signal
  fRes=param[1];   // relative resolution
  fRange=param[2]; // PID "range" (in sigmas)
}

Double_t AliTRDpidESD::Bethe(Double_t bg) 
{
  //
  // Parametrization of the Bethe-Bloch-curve
  // The parametrization is the same as for the TPC and is taken from Lehrhaus.
  //

  // This parameters have been adjusted to averaged values from GEANT
  const Double_t kP1 = 7.17960e-02;
  const Double_t kP2 = 8.54196;
  const Double_t kP3 = 1.38065e-06;
  const Double_t kP4 = 5.30972;
  const Double_t kP5 = 2.83798;

  // This parameters have been adjusted to Xe-data found in:
  // Allison & Cobb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kP1 = 0.76176E-1;
  //const Double_t kP2 = 10.632;
  //const Double_t kP3 = 3.17983E-6;
  //const Double_t kP4 = 1.8631;
  //const Double_t kP5 = 1.9479;

  // Lower cutoff of the Bethe-Bloch-curve to limit step sizes
  const Double_t kBgMin = 0.8;
  const Double_t kBBMax = 6.83298;
  //const Double_t kBgMin = 0.6;
  //const Double_t kBBMax = 17.2809;
  //const Double_t kBgMin = 0.4;
  //const Double_t kBBMax = 82.0;

  if (bg > kBgMin) {
    Double_t yy = bg / TMath::Sqrt(1. + bg*bg);
    Double_t aa = TMath::Power(yy,kP4);
    Double_t bb = TMath::Power((1./bg),kP5);
             bb = TMath::Log(kP3 + bb);
    return ((kP2 - aa - bb)*kP1 / aa);
  }
  else {
    return kBBMax;
  }

}

//_________________________________________________________________________
Int_t AliTRDpidESD::MakePID(AliESD *event)
{
  //
  //  This function calculates the "detector response" PID probabilities 
  //
  static const Double_t masses[]={
    0.000511, 0.105658, 0.139570, 0.493677, 0.938272, 1.875613
  };
  Int_t ntrk=event->GetNumberOfTracks();
  for (Int_t i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    if ((t->GetStatus()&AliESDtrack::kTRDin)==0)
       if ((t->GetStatus()&AliESDtrack::kTRDout)==0)
          if ((t->GetStatus()&AliESDtrack::kTRDrefit)==0) continue;
    Int_t ns=AliESDtrack::kSPECIES;
    Double_t p[10];
    for (Int_t j=0; j<ns; j++) {
      Double_t mass=masses[j];
      Double_t mom=t->GetP();
      Double_t dedx=t->GetTRDsignal()/fMIP;
      Double_t bethe=Bethe(mom/mass); 
      Double_t sigma=fRes*bethe;
      if (TMath::Abs(dedx-bethe) > fRange*sigma) {
	p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
        continue;
      }
      p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
    }
    t->SetTRDpid(p);
  }
  return 0;
}
