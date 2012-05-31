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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// HMPID tender: - recalculate pid bit using tighter cuts
//               - this is needed for all 2010 and 11a-c data
//             
//             
// Contacts: Giacomo.Volpe@ba.infn.it                                        //
//           Jens.Wiechula@cern.ch                                           //
///////////////////////////////////////////////////////////////////////////////
#include <TMath.h>
#include <TRandom.h>
#include <AliLog.h>
#include <AliESDEvent.h>
#include <AliESDtrack.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliESDpid.h>
#include <AliTender.h>

#include "AliHMPIDTenderSupply.h"

ClassImp(AliHMPIDTenderSupply)

AliHMPIDTenderSupply::AliHMPIDTenderSupply() :
AliTenderSupply()
{
  //
  // default ctor
  //
}

//_____________________________________________________
AliHMPIDTenderSupply::AliHMPIDTenderSupply(const char *name, const AliTender *tender) :
AliTenderSupply(name,tender)
{
  //
  // named ctor
  //
}

//_____________________________________________________
void AliHMPIDTenderSupply::Init()
{
  // Initialise HMPID tender
  
}

//_____________________________________________________
void AliHMPIDTenderSupply::ProcessEvent()
{
  //
  // recalculate HMPIDpid bit
  // 
  

  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;
  
  // re-evaluate the HMPIDpid bit for all tracks
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliESDtrack *track=event->GetTrack(itrack);
    if (!itrack) continue;
    //reset pid bit first
    track->ResetStatus(AliESDtrack::kHMPIDpid);

    Float_t xPc=0., yPc=0., xMip=0., yMip=0., thetaTrk=0., phiTrk=0.;
    Int_t nPhot=0, qMip=0;
 
    track->GetHMPIDtrk(xPc,yPc,thetaTrk,phiTrk);
    track->GetHMPIDmip(xMip,yMip,qMip,nPhot);
    //
    //make cuts, just an example, THIS NEEDS TO BE CHANGED
    //
    //if ((track->GetStatus()&AliESDtrack::kHMPIDout)!=AliESDtrack::kHMPIDout) continue;
   
    Float_t dist = TMath::Sqrt((xPc-xMip)*(xPc-xMip) + (yPc-yMip)*(yPc-yMip));    

    if(dist > 0.7 || nPhot> 30 || qMip < 100  ) continue;

    //set pid bit, track was accepted
    track->SetStatus(AliESDtrack::kHMPIDpid);
    
  }
  
  
}
