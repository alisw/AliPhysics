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

#include "AliESDtrack.h"
#include "AliESD.h"

#include "AliTOFpidESD.h"

ClassImp(AliTOFpidESD)

//_________________________________________________________________________
AliTOFpidESD::AliTOFpidESD(Double_t *param) {
  //
  //  The main constructor
  //
  fN=0; fEventN=0;

  fSigma=param[0];
  fRange=param[1];

}

//_________________________________________________________________________
Int_t AliTOFpidESD::MakePID(AliESD *event)
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
    for (Int_t j=0; j<AliPID::kSPECIES; j++) {
      Double_t mass=AliPID::ParticleMass(j);
      Double_t dpp=0.01;      //mean relative pt resolution;
      if (mom>0.5) dpp=0.01*mom;
      Double_t sigma=dpp*time[j]/(1.+ mom*mom/(mass*mass));
      sigma=TMath::Sqrt(sigma*sigma + fSigma*fSigma);
      if (TMath::Abs(tof-time[j]) > fRange*sigma) {
	p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
        continue;
      }
      p[j]=TMath::Exp(-0.5*(tof-time[j])*(tof-time[j])/(sigma*sigma))/sigma;
    }
    t->SetTOFpid(p);
  }

  delete[] tracks;
  
  return 0;
}

