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
   
   fInputs.AddLast( new AliTriggerInput( "START_A_L0", "Signal on T0-A",  0x01 ) );
   fInputs.AddLast( new AliTriggerInput( "START_C_L0", "Signal on T0-C", 0x02 ) );
   fInputs.AddLast( new AliTriggerInput( "START_Vertex_L0", " Vertex T0-C&T0-A ", 0x04 ) );
   fInputs.AddLast( new AliTriggerInput( "START_Centr_L0", "Centrality central",  0x08 ) );
   fInputs.AddLast( new AliTriggerInput( "START_SemiCentral_L0", "Centrality semicentral",  0x08 ) );
}

//----------------------------------------------------------------------
void AliSTARTTrigger::Trigger()
{

   AliRunLoader* runLoader = gAlice->GetRunLoader();
   AliLoader * fSTARTLoader = runLoader->GetLoader("STARTLoader");
   //   AliSTARTdigit *fDigits; 
   fSTARTLoader->LoadDigits("READ");
   // Creating START data container
   Int_t idebug = 1;
   TTree* treeD = fSTARTLoader->TreeD();
  if (!treeD) {
    AliError("no digits tree");
    return;
  }

  TBranch *brDigits = treeD->GetBranch("START");
  if (brDigits) {
    brDigits->SetAddress(&fDigits);
  }else{
    AliError("Branch START DIGIT not found");
    exit(111);
  } 
  brDigits->GetEntry(0);
  Int_t   besttimeright = fDigits->BestTimeRight();
  Int_t   besttimeleft = fDigits->BestTimeLeft();
  Int_t   timeDiff = besttimeright-besttimeleft;
  Int_t   meanTime = (besttimeright+besttimeleft)/2;
  Int_t sumMult=   fDigits->SumMult();
  cout<<besttimeright<<" "<<besttimeleft<<" "<< timeDiff<<" "<<meanTime<<endl;

  Bool_t trigg[5];
  if (besttimeright<1000) trigg[0]=kTRUE;
  if (besttimeleft<1000)  trigg[1]=kTRUE;
  if (timeDiff<5)      trigg[2]=kTRUE;
  if (sumMult>2300)    trigg[3]=kTRUE;
  if (sumMult>1800)    trigg[4]=kTRUE;


   if( trigg[0] )  SetInput("START_T0A_L0");
   if( trigg[1] )  SetInput("START_T0C_L0");
   if( trigg[2] )  SetInput("START_T0vertex_L0");
   if( trigg[3] )  SetInput("START_T0Centr_L0");
   if( trigg[4] )  SetInput("START_T0SemiCentral_L0");
   cout<<" Trigger "<<trigg[0]<<" "<<trigg[1]<<endl;
   
}
