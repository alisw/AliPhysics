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

Double_t AliTRDpidESD::Bethe(Double_t bg) {
  //
  // This is the Bethe-Bloch function normalised to 1 at the minimum
  //
  Double_t bg2=bg*bg;
  Double_t bethe;
  if (bg<3.5e1) 
      bethe=(1.+ bg2)/bg2*(log(5940*bg2) - bg2/(1.+ bg2));
  else // Density effect ( approximately :) 
      bethe=(1.+ bg2)/bg2*(log(3.5*5940*bg) - bg2/(1.+ bg2));
  return bethe/11.091;
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
