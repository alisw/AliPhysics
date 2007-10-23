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

#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"

ClassImp(AliESDpid)

//_________________________________________________________________________
Int_t AliESDpid::MakePID(AliESDEvent *event)
{
  //
  // Combine the information of various detectors
  // to determine the Particle Identification
  //
  Int_t ntrk=event->GetNumberOfTracks();
  for (Int_t i=0; i<ntrk; i++) {
    Int_t ns=AliPID::kSPECIES;
    Double_t p[10]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

    AliESDtrack *t=event->GetTrack(i);

    if (t->IsOn(AliESDtrack::kITSpid)) {
      Double_t d[10];
      t->GetITSpid(d);
      for (Int_t j=0; j<ns; j++) p[j]*=d[j];
    }

    if (t->IsOn(AliESDtrack::kTPCpid)) {
      Double_t d[10];
      t->GetTPCpid(d);
      for (Int_t j=0; j<ns; j++) p[j]*=d[j];
    }

    if (t->IsOn(AliESDtrack::kTRDpid)) {
      Double_t d[10];
      t->GetTRDpid(d);
      for (Int_t j=0; j<ns; j++) p[j]*=d[j];
    }

    if (t->IsOn(AliESDtrack::kTOFpid)) {
      Double_t d[10];
      t->GetTOFpid(d);
      for (Int_t j=0; j<ns; j++) p[j]*=d[j];
    }

    if (t->IsOn(AliESDtrack::kHMPIDpid)) {
      Double_t d[10];
      t->GetHMPIDpid(d);
      for (Int_t j=0; j<ns; j++) p[j]*=d[j];
    }

    t->SetESDpid(p);
  }

  return 0;
}
