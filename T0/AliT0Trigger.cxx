/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                         *
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

#include <Riostream.h>
#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTriggerInput.h"

#include "AliT0.h"
#include "AliT0Loader.h"
#include "AliT0digit.h"
#include "AliT0Trigger.h"

//----------------------------------------------------------------------
ClassImp(AliT0Trigger)

//----------------------------------------------------------------------
AliT0Trigger::AliT0Trigger()
  : AliTriggerDetector() 
{
   SetName("T0");
   CreateInputs();
}

//----------------------------------------------------------------------
void AliT0Trigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "START_A_L0", "Signal on T0-A",  0x0100 ) );
   fInputs.AddLast( new AliTriggerInput( "START_C_L0", "Signal on T0-C", 0x0200 ) );
   fInputs.AddLast( new AliTriggerInput( "START_Vertex_L0", " Vertex T0-C&T0-A ", 0x0400 ) );
   fInputs.AddLast( new AliTriggerInput( "START_Centr_L0", "Centrality central",  0x0800 ) );
   fInputs.AddLast( new AliTriggerInput( "START_SemiCentral_L0", "Centrality semicentral",  0x1000 ) );
}

//----------------------------------------------------------------------
void AliT0Trigger::Trigger()
{
   AliRunLoader* runLoader = gAlice->GetRunLoader();
   AliLoader * fT0Loader = runLoader->GetLoader("T0Loader");
   //   AliT0digit *fDigits; 
   fT0Loader->LoadDigits("READ");
   // Creating T0 data container

   TTree* treeD = fT0Loader->TreeD();
  if (!treeD) {
    AliError("no digits tree");
    return;
  }
  AliT0digit *fDigits = new AliT0digit();

  TBranch *brDigits = treeD->GetBranch("T0");
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    AliError("Branch T0 DIGIT not found");
    return;
  } 
  brDigits->GetEntry(0);
  Int_t   besttimeright = fDigits->BestTimeRight();
  Int_t   besttimeleft = fDigits->BestTimeLeft();
  Int_t   timeDiff = fDigits->TimeDiff();
  Int_t    sumMult=   fDigits->SumMult();

  if (besttimeright > 0 && besttimeright <99999  )  SetInput("START_A_L0");
  if (besttimeleft>0  && besttimeleft<99999)   SetInput("START_C_L0"); 
  if (timeDiff >5500 && timeDiff < 6500)       SetInput("START_Vertex_L0");
  if (sumMult > 175)                           SetInput("START_Centr_L0");
  if (sumMult>155 && sumMult <= 175)           SetInput("START_SemiCentral_L0");;

   
}
