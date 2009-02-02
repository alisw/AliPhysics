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

/* $Id$ */

#include <TClonesArray.h>

#include "AliRun.h"
#include "AliRunLoader.h"

#include "AliVZEROTrigger.h"
#include "AliVZEROTriggerMask.h"

//______________________________________________________________________
ClassImp(AliVZEROTrigger)
////////////////////////////////////////////////////////////////////////
//
// Version 1
//
// AliVZEROTrigger: 
//
////////////////////////////////////////////////////////////////////////

//______________________________________________________________________

AliVZEROTrigger::AliVZEROTrigger()
  :AliTriggerDetector(),
   fAdcThresHold(0.0),
   fTimeWindowWidthBBA(50.0),
   fTimeWindowWidthBGA(20.0),
   fTimeWindowWidthBBC(50.0),
   fTimeWindowWidthBGC(20.0)
   
{
   SetName("VZERO");
   CreateInputs();

   SetAdcThreshold();
}
//______________________________________________________________________
void AliVZEROTrigger::CreateInputs()
{
   // inputs

   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;

   fInputs.AddLast( new AliTriggerInput( "VZERO_LEFT", "VZERO", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "VZERO_RIGHT","VZERO", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "VZERO_AND",  "VZERO", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "VZERO_OR",   "VZERO", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "VZERO_BEAMGAS", "VZERO", 0 ) );
}

//______________________________________________________________________
void AliVZEROTrigger::Trigger()
{
  
  //  ********** Get run loader for the current event **********
  AliRunLoader* runLoader = AliRunLoader::Instance();

  AliVZEROLoader* loader = 
    (AliVZEROLoader* )runLoader->GetLoader( "VZEROLoader" );

  loader->LoadDigits("READ");
  TTree* vzeroDigitsTree = loader->TreeD();
  if (!vzeroDigitsTree) return;

  TClonesArray* vzeroDigits = new TClonesArray("AliVZEROdigit",1000);
  TBranch* digitBranch = vzeroDigitsTree->GetBranch("VZERODigit");
  digitBranch->SetAddress(&vzeroDigits);

  AliVZEROTriggerMask *TriggerMask = new AliVZEROTriggerMask();
  TriggerMask->SetAdcThreshold(fAdcThresHold);
  TriggerMask->SetTimeWindowWidthBBA(fTimeWindowWidthBBA);
  TriggerMask->SetTimeWindowWidthBGA(fTimeWindowWidthBGA);
  TriggerMask->SetTimeWindowWidthBBC(fTimeWindowWidthBBC);
  TriggerMask->SetTimeWindowWidthBGC(fTimeWindowWidthBGC);
  TriggerMask->FillMasks(vzeroDigitsTree,vzeroDigits);

  if ( (TriggerMask->GetBGtriggerV0A()>0) ||
       (TriggerMask->GetBGtriggerV0C()>0)) SetInput( "VZERO_BEAMGAS" );
  if (TriggerMask->GetBBtriggerV0A()>0)  SetInput( "VZERO_LEFT" );
  if (TriggerMask->GetBBtriggerV0C()>0)  SetInput( "VZERO_RIGHT" );
  if ( (TriggerMask->GetBBtriggerV0A()>0) ||
       (TriggerMask->GetBBtriggerV0C()>0)) SetInput( "VZERO_OR" );
  if ( (TriggerMask->GetBBtriggerV0A()>0) &&
       (TriggerMask->GetBBtriggerV0C()>0)) SetInput( "VZERO_AND" );

  return;
}

