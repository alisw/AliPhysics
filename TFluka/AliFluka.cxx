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

/*
$Log$
*/

#include "AliFluka.h"
#include "AliRun.h"
#include "AliGenerator.h"

ClassImp(AliFluka)

void AliFluka::ProcessRun(Int_t nevent)
{
  //
  // Process the run
  //
  Int_t todo = TMath::Abs(nevent);
  for (Int_t i=0; i<todo; i++) {
  // Process one run (one run = one event)
     gAlice->BeginEvent();
     ProcessEvent();
     gAlice->FinishEvent();
  }
}



void AliFluka::ProcessEvent()
{
  //
  // Process one event
  //
    gAlice->Generator()->Generate();
}
