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
// PID tender: Do combined PID                                               //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <AliESDpid.h>
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include "AliTender.h"

#include "AliPIDTenderSupply.h"

AliPIDTenderSupply::AliPIDTenderSupply() :
  AliTenderSupply()
{
  //
  // default ctor
  //
}

//_____________________________________________________
AliPIDTenderSupply::AliPIDTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender)
{
  //
  // named ctor
  //
}

//_____________________________________________________
void AliPIDTenderSupply::ProcessEvent()
{
  //
  // Combine PID information
  //

  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;

  AliESDpid *pid=fTender->GetESDhandler()->GetESDpid();
  if (!pid) return;
  //
  // recalculate combined PID probabilities
  //
  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++)
    pid->CombinePID(event->GetTrack(itrack));
  
}
