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

#include "AliSTART.h"
#include "AliSTARTLoader.h"
#include "AliSTARTdigit.h"
#include "AliSTARTTrigger.h"

//----------------------------------------------------------------------
ClassImp(AliSTARTTrigger)

//----------------------------------------------------------------------
AliSTARTTrigger::AliSTARTTrigger()
  : AliTriggerDetector() 
{
   SetName("START");
   CreateInputs();
}

//----------------------------------------------------------------------
void AliSTARTTrigger::CreateInputs()
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
void AliSTARTTrigger::Trigger()
{
   AliRunLoader* runLoader = gAlice->GetRunLoader();
   AliLoader * fSTARTLoader = runLoader->GetLoader("STARTLoader");
   //   AliSTARTdigit *fDigits; 
   fSTARTLoader->LoadDigits("READ");
   // Creating START data container

   TTree* treeD = fSTARTLoader->TreeD();
  if (!treeD) {
    AliError("no digits tree");
    return;
  }
  AliSTARTdigit *fDigits = new AliSTARTdigit();

  TBranch *brDigits = treeD->GetBranch("START");
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    AliError("Branch START DIGIT not found");
    return;
  } 
  brDigits->GetEntry(0);
  Int_t   besttimeright = fDigits->BestTimeRight();
  Int_t   besttimeleft = fDigits->BestTimeLeft();
  Int_t   timeDiff = besttimeright-besttimeleft;
  Int_t    sumMult=   fDigits->SumMult();

  if (besttimeright>0 && besttimeright<10000)  SetInput("START_A_L0");
  if (besttimeleft>0  && besttimeleft<10000)   SetInput("START_C_L0"); 
  if (TMath::Abs(timeDiff) < 7000)                SetInput("START_Vertex_L0");
  if (sumMult>2300)                            SetInput("START_Centr_L0");
  if (sumMult>1800)                            SetInput("START_SemiCentral_L0");;

   
}
