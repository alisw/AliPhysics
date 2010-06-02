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
// TRD tender: reapply pid on the fly                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <AliLog.h>
#include <TTree.h>
#include <TChain.h>
#include <AliPID.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliESDpid.h>
#include <AliESDtrack.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliTRDPIDResponse.h>
#include <AliTender.h>

#include "AliTRDTenderSupply.h"


AliTRDTenderSupply::AliTRDTenderSupply() :
  AliTenderSupply(),
  fESDpid(0x0),
  fPIDmethod(k1DLQpid)
{
  //
  // default ctor
  //
}

//_____________________________________________________
AliTRDTenderSupply::AliTRDTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fESDpid(0x0),
  fPIDmethod(k1DLQpid)
{
  //
  // named ctor
  //
}

AliTRDTenderSupply::~AliTRDTenderSupply()
{
  //
  // dtor
  //
}

//_____________________________________________________
void AliTRDTenderSupply::Init()
{
  //
  // Initialise TRD tender
  //

  //
  // Set event information
  //
  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  
  // 1DLQ PID implemented in the AliESD object
  fESDpid=fTender->GetESDhandler()->GetESDpid();
  if (!fESDpid) {
    fESDpid=new AliESDpid;
    fTender->GetESDhandler()->SetESDpid(fESDpid);
  }
  // Set Normalisation Factors
  if(mgr->GetMCtruthEventHandler()){
    // Assume MC
    fESDpid->GetTRDResponse().SetGainNormalisationFactor(1.8315);
  }
  else{
    // Assume Data
    fESDpid->GetTRDResponse().SetGainNormalisationFactor(1.14);
  }
  
}

//_____________________________________________________
void AliTRDTenderSupply::ProcessEvent()
{
  //
  // Reapply pid information
  //


  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;
  Int_t ntracks=event->GetNumberOfTracks();
  
  //
  // recalculate PID probabilities
  //
  switch(fPIDmethod){
    case k1DLQpid:
      for(Int_t itrack = 0; itrack < ntracks; itrack++) 
        fESDpid->MakeTRDPID(event->GetTrack(itrack));
      break;
    default:
      AliError("PID Method not implemented (yet)");
  }
}
