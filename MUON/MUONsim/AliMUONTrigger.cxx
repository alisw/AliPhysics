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

#include <TTree.h>
#include "AliMUONTrigger.h"

#include "AliLog.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONVTriggerStore.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliTriggerInput.h"


//-----------------------------------------------------------------------------
/// \class AliMUONTrigger
///
/// Implementation of AliTriggerDetector for MUON detector
///
/// So far, the inputs are taken from AliMUONTriggerDecision object
/// April 06, E.L.T.
/// May 06, taken info from Global Trigger branch (Ch.F)
///
/// \author E. Lopez Torres
//-----------------------------------------------------------------------------

//----------------------------------------------------------------------
/// \cond CLASSIMP
ClassImp(AliMUONTrigger)
/// \endcond

//----------------------------------------------------------------------
AliMUONTrigger::AliMUONTrigger()
  : AliTriggerDetector(), fTriggerStore(0x0)
{
/// Default constructor

   SetName("MUON");
   CreateInputs();
}

//----------------------------------------------------------------------
AliMUONTrigger::~AliMUONTrigger()
{
  /// Destructor
  delete fTriggerStore;
}

//----------------------------------------------------------------------
void AliMUONTrigger::CreateInputs()
{
   /// inputs 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "0MSL", "MUONTRG",  0 ) );
   fInputs.AddLast( new AliTriggerInput( "0MSH", "MUONTRG", 0 ) );

   fInputs.AddLast( new AliTriggerInput( "0MUL", "MUONTRG",  0 ) );
   fInputs.AddLast( new AliTriggerInput( "0MUH", "MUONTRG", 0 ) );

   fInputs.AddLast( new AliTriggerInput( "0MLL", "MUONTRG",  0 ) );
   fInputs.AddLast( new AliTriggerInput( "0MLH", "MUONTRG", 0 ) );
}

//----------------------------------------------------------------------
void AliMUONTrigger::Trigger()
{
  /// sets the trigger inputs

   AliRunLoader* runLoader = AliRunLoader::Instance();
  
   AliLoader * muonLoader = runLoader->GetDetectorLoader("MUON");
   muonLoader->LoadDigits("READ");

   TTree* treeD = muonLoader->TreeD();
   
   if (!treeD)
   {
     AliError("No TreeD available. Cannot make trigger");
     return;
   }
   
   if (!fTriggerStore) 
   {  
     fTriggerStore = AliMUONVTriggerStore::Create(*treeD);
     if (!fTriggerStore)
     {
       AliError("Could not create triggerStore from treeD");
       return;
     }     
   }

   Bool_t ok = fTriggerStore->Connect(*treeD,kTRUE);

   if (!ok)
   {
     AliError("Could not read trigger from TreeD !");
     return;
   }
   
   treeD->GetEvent(0);
   
   AliMUONGlobalTrigger* globalTrigger = fTriggerStore->Global();
   if (globalTrigger == 0x0) 
   { 
     AliWarning("No Global Trigger available");
   }
   else
   {
     // set CTP
     if (globalTrigger->SingleLpt())      SetInput("0MSL");
     if (globalTrigger->SingleHpt())      SetInput("0MSH");
     
     if (globalTrigger->PairUnlikeLpt())  SetInput("0MUL");
     if (globalTrigger->PairUnlikeHpt())  SetInput("0MUH");
     
     if (globalTrigger->PairLikeLpt())    SetInput("0MLL");
     if (globalTrigger->PairLikeHpt())    SetInput("0MLH");
   }
   muonLoader->UnloadDigits();
   fTriggerStore->Clear();
}
