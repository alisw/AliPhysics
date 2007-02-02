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

#include <TClonesArray.h>

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
   
   fInputs.AddLast( new AliTriggerInput( "MUON_Single_LPt_L0", "Single Low Pt",  0x01 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_Single_HPt_L0", "Single High Pt", 0x02 ) );

   fInputs.AddLast( new AliTriggerInput( "MUON_Unlike_LPt_L0", "Dimuon Unlike Sign pair Low Pt",  0x04 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_Unlike_HPt_L0", "Dimuon Unlike Sign pair High Pt", 0x08 ) );

   fInputs.AddLast( new AliTriggerInput( "MUON_Like_LPt_L0", "Dimuon Like Sign pair Low Pt",  0x10 ) );
   fInputs.AddLast( new AliTriggerInput( "MUON_Like_HPt_L0", "Dimuon Like Sign pair High Pt", 0x20 ) );
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
   if (globalTrigger->SingleLpt())      SetInput("MUON_Single_LPt_L0");
   if (globalTrigger->SingleHpt())      SetInput("MUON_Single_HPt_L0");
   
   if (globalTrigger->PairUnlikeLpt())  SetInput("MUON_Unlike_LPt_L0");
   if (globalTrigger->PairUnlikeHpt())  SetInput("MUON_Unlike_HPt_L0");
   
   if (globalTrigger->PairLikeLpt())    SetInput("MUON_Like_LPt_L0");
   if (globalTrigger->PairLikeHpt())    SetInput("MUON_Like_HPt_L0");

   muonData->ResetTrigger();
   muonLoader->UnloadDigits();

}
