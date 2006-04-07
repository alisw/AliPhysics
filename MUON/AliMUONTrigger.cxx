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

#include "AliMUON.h"
#include "AliMUONLoader.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONTriggerDecision.h"
#include "AliMUONTrigger.h"

//----------------------------------------------------------------------
ClassImp(AliMUONTrigger)

//----------------------------------------------------------------------
AliMUONTrigger::AliMUONTrigger()
  : AliTriggerDetector() 
{
   SetName("MUON");
   CreateInputs();
}

//----------------------------------------------------------------------
void AliMUONTrigger::CreateInputs()
{
   // inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "MUON_SPlus_LPt_L0", "Single Plus Low Pt",  0x01 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_SPlus_HPt_L0", "Single Plus High Pt", 0x02 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_SPlus_All_L0", "Single Plus All",     0x04 ) );

   fInputs.AddLast( new AliTriggerInput( "MUON_SMinus_LPt_L0", "Single Minus Low Pt",  0x08 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_SMinus_HPt_L0", "Single Minus High Pt", 0x10 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_SMinus_All_L0", "Single Minus All",     0x20 ) );

   fInputs.AddLast( new AliTriggerInput( "MUON_SUndef_LPt_L0", "Single Undefined Low Pt",  0x40 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_SUndef_HPt_L0", "Single Undefined High Pt", 0x80 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_SUndef_All_L0", "Single Undefined All",     0x100 ) );

   fInputs.AddLast( new AliTriggerInput( "MUON_Unlike_LPt_L0", "Single Unlike Sign pair Low Pt",  0x200 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_Unlike_HPt_L0", "Single Unlike Sign pair High Pt", 0x400 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_Unlike_All_L0", "Single Unlike Sign pair All",     0x800 ) );

   fInputs.AddLast( new AliTriggerInput( "MUON_Like_LPt_L0", "Single Like Sign pair Low Pt",  0x1000 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_Like_HPt_L0", "Single Like Sign pair High Pt", 0x2000 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_Like_All_L0", "Single Like Sign pair All",     0x4000 ) );
}

//----------------------------------------------------------------------
void AliMUONTrigger::Trigger()
{


   AliRunLoader* runLoader = gAlice->GetRunLoader();

   AliLoader * muonLoader = runLoader->GetLoader("MUONLoader");
   muonLoader->LoadDigits("READ");
   // Creating MUON data container
   AliMUONData* muonData = new AliMUONData(muonLoader,"MUON","MUON");
   muonData->SetTreeAddress("D");
   muonData->GetDigits();
   Int_t idebug = 1;
   // Creating MUONTriggerDecision
   AliMUONTriggerDecision* decision = new AliMUONTriggerDecision(muonLoader ,idebug,muonData);
   AliMUONDigit * mDigit;
   Int_t tracks[10];
   Int_t charges[10];
   Int_t digits[7];

   for(Int_t ichamber=10; ichamber<14; ichamber++) {
      Int_t idigit, ndigits;
      ndigits = (Int_t) muonData->Digits(ichamber)->GetEntriesFast();
//            printf(">>> Chamber Cathode ndigits %d %d %d\n",ichamber,icathode,ndigits);
      for(idigit=0; idigit<ndigits; idigit++) {
         mDigit = static_cast<AliMUONDigit*>(muonData->Digits(ichamber)->At(idigit));
         digits[0] = mDigit->PadX();
         digits[1] = mDigit->PadY();
         digits[2] = mDigit->Cathode();
         digits[3] = mDigit->Signal();
         digits[4] = mDigit->Physics();
         digits[5] = mDigit->Hit();
         digits[6] = mDigit->DetElemId();

         Int_t digitindex = 0 ;
  //       printf("ichamber ix iy %d %d %d \n",ichamber,mDigit->PadX(),mDigit->PadY());

         decision->AddDigit(ichamber, tracks, charges, digits, digitindex );
      } // loop on digits
   } // loop on chambers
   muonData->ResetDigits();
   decision->Trigger();
   decision->ClearDigits();

   // Set the trigger inputs = "global decision" 
   Int_t singlePlus[3];  // tot num of single plus
   Int_t singleMinus[3]; // tot num of single minus
   Int_t singleUndef[3]; // tot num of single undefined
   Int_t pairUnlike[3];  // tot num of unlike-sign pairs
   Int_t pairLike[3];    // tot num of like-sign pairs
   decision->GetGlobalTrigger(singlePlus, singleMinus, singleUndef, pairUnlike, pairLike);

   if( singlePlus[0] )  SetInput("MUON_SPlus_LPt_L0");
   if( singlePlus[1] )  SetInput("MUON_SPlus_HPt_L0");
   if( singlePlus[2] )  SetInput("MUON_SPlus_All_L0");
   
   if( singleMinus[0] ) SetInput("MUON_SMinus_LPt_L0");
   if( singleMinus[1] ) SetInput("MUON_SMinus_HPt_L0");
   if( singleMinus[2] ) SetInput("MUON_SMinus_All_L0");
   
   if( singleUndef[0] ) SetInput("MUON_SUndef_LPt_L0");
   if( singleUndef[1] ) SetInput("MUON_SUndef_HPt_L0");
   if( singleUndef[2] ) SetInput("MUON_SUndef_All_L0");
   
   if( pairUnlike[0] )  SetInput("MUON_Unlike_LPt_L0");
   if( pairUnlike[1] )  SetInput("MUON_Unlike_HPt_L0");
   if( pairUnlike[2] )  SetInput("MUON_Unlike_All_L0");
   
   if( pairLike[0] )    SetInput("MUON_Like_LPt_L0");
   if( pairLike[1] )    SetInput("MUON_Like_HPt_L0");
   if( pairLike[2] )    SetInput("MUON_Like_All_L0");
}
