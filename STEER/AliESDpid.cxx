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
//           Implementation of the combined PID class
//
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "AliESDpid.h"
#include "AliESD.h"
#include "AliESDtrack.h"

ClassImp(AliESDpid)

//_________________________________________________________________________
Int_t AliESDpid::MakePID(AliESD *event)
{
  Int_t ntrk=event->GetNumberOfTracks();
  for (Int_t i=0; i<ntrk; i++) {
    Int_t ns=AliESDtrack::kSPECIES;
    Double_t p[10]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

    AliESDtrack *t=event->GetTrack(i);

    if ((t->GetStatus()&AliESDtrack::kITSpid )!=0) {
      Double_t d[10];
      t->GetITSpid(d);
      for (Int_t j=0; j<ns; j++) p[j]*=d[j];
    }

    if ((t->GetStatus()&AliESDtrack::kTPCpid )!=0) {
      Double_t d[10];
      t->GetTPCpid(d);
      for (Int_t j=0; j<ns; j++) p[j]*=d[j];
    }

    if ((t->GetStatus()&AliESDtrack::kTRDpid )!=0) {
      Double_t d[10];
      t->GetTRDpid(d);
      for (Int_t j=0; j<ns; j++) p[j]*=d[j];
    }

    if (t->GetP()>0.7) // accept the TOF only for the high momenta
      if ((t->GetStatus()&AliESDtrack::kTOFpid )!=0) {
	Double_t d[10];
	t->GetTOFpid(d);
	for (Int_t j=0; j<ns; j++) p[j]*=d[j];
      }

    t->SetESDpid(p);
  }
  return 0;
}
