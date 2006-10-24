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
#include "AliMUONGlobalTrigger.h"
#include "AliMUONTrigger.h"

///
/// \class AliMUONTrigger
///
/// Implementation of AliTriggerDetector for MUON detector
///
/// So far, the inputs are taken from AliMUONTriggerDecision object
/// April 06, E.L.T.
/// May 06, taken info from Global Trigger branch (Ch.F)

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
  // sets the trigger inputs
  AliMUONGlobalTrigger* globalTrigger;
  TClonesArray* globalTriggerArray;

   AliRunLoader* runLoader = gAlice->GetRunLoader();

   AliLoader * muonLoader = runLoader->GetLoader("MUONLoader");
   muonLoader->LoadDigits("READ");

   // Creating MUON data container
   AliMUONData* muonData = new AliMUONData(muonLoader,"MUON","MUON");

   // get global info
   muonData->SetTreeAddress("GLT");
   muonData->GetTriggerD();
   globalTriggerArray = muonData->GlobalTrigger(); 
   if (globalTriggerArray == 0x0) { 
     AliWarning("No Global Trigger Array available");
     return;
   }
   globalTrigger = (AliMUONGlobalTrigger*)globalTriggerArray->UncheckedAt(0);

   if (globalTrigger == 0x0) { 
     AliWarning("No Global Trigger available");
     return;
   }
   // set CTP
   if (globalTrigger->SinglePlusLpt())  SetInput("MUON_SPlus_LPt_L0");
   if (globalTrigger->SinglePlusHpt())  SetInput("MUON_SPlus_HPt_L0");
   
   if (globalTrigger->SingleMinusLpt()) SetInput("MUON_SMinus_LPt_L0");
   if (globalTrigger->SingleMinusHpt()) SetInput("MUON_SMinus_HPt_L0");
   
   if (globalTrigger->SingleUndefLpt()) SetInput("MUON_SUndef_LPt_L0");
   if (globalTrigger->SingleUndefHpt()) SetInput("MUON_SUndef_HPt_L0");
   
   if (globalTrigger->PairUnlikeLpt())  SetInput("MUON_Unlike_LPt_L0");
   if (globalTrigger->PairUnlikeHpt())  SetInput("MUON_Unlike_HPt_L0");
   
   if (globalTrigger->PairLikeLpt())    SetInput("MUON_Like_LPt_L0");
   if (globalTrigger->PairLikeHpt())    SetInput("MUON_Like_HPt_L0");

   muonData->ResetTrigger();
   muonLoader->UnloadDigits();

}
