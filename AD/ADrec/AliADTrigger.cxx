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
/* $Id: AliADTrigger.cxx 49869 2011-05-18 04:49:51Z hristov $ */
// ---------------------
// Class AliADTrigger
// ---------------------
// Top class to simulate the AD trigger response
// This class is only used for interface with AliTriggerDetector
// Its create and Set  Inputs of the CTP
//


#include <TClonesArray.h>
#include <TTree.h>

#include "AliRun.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliTriggerInput.h"

#include "AliADdigit.h"
#include "AliADTrigger.h"
#include "AliADTriggerSimulator.h"

//______________________________________________________________________
ClassImp(AliADTrigger)

//______________________________________________________________________

AliADTrigger::AliADTrigger():AliTriggerDetector()
{
   SetName("AD");
   CreateInputs();
}
//______________________________________________________________________
void AliADTrigger::CreateInputs()
{
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;

   fInputs.AddLast( new AliTriggerInput( "AD_BBA_AND_BBC", "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_BBA_OR_BBC","AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_BGA_AND_BBC",  "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "0UGA",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_BGC_AND_BBA", "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "0UGC",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_CTA1_AND_CTC1",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_CTA1_OR_CTC1",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_CTA2_AND_CTC2",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_CTA2_OR_CTC2",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_MTA_AND_MTC",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_MTA_OR_MTC",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "0UBA",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "0UBC",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_BGA_OR_BGC",   "AD", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "AD_BEAMGAS",   "AD", 0 ) );
}

//______________________________________________________________________
void AliADTrigger::Trigger()
{
  //  ********** Get run loader for the current event **********
   AliRunLoader* runLoader = AliRunLoader::Instance();

   AliLoader* loader = runLoader->GetLoader( "ADLoader" );

   if(!loader) {
      AliError("Can not get AD loader");
      return;
   }
   loader->LoadDigits("update");
   TTree* ADDigitsTree = loader->TreeD();

   if (!ADDigitsTree) {
      AliError("Can not get the AD digit tree");
      return;
   }
   TClonesArray* ADDigits = NULL;
   TBranch* digitBranch = ADDigitsTree->GetBranch("ADDigit");
   digitBranch->SetAddress(&ADDigits);
   
   AliADTriggerSimulator * triggerSimulator = new AliADTriggerSimulator(ADDigitsTree,ADDigits);
   triggerSimulator->Run();
	
   loader->WriteDigits("OVERWRITE");  
   loader->UnloadDigits();     

   if(triggerSimulator->GetBBAandBBC())   SetInput( "AD_BBA_AND_BBC" );
   if(triggerSimulator->GetBBAorBBC())      SetInput( "AD_BBA_OR_BBC" );
   if(triggerSimulator->GetBGAandBBC())   SetInput( "AD_BGA_AND_BBC" );
   if(triggerSimulator->GetBGA())         SetInput( "0UGA" );
   if(triggerSimulator->GetBGCandBBA())   SetInput( "AD_BGC_AND_BBA" );
   if(triggerSimulator->GetBGC())         SetInput( "0UGC" );
   if(triggerSimulator->GetCTA1andCTC1())   SetInput( "AD_CTA1_AND_CTC1" );
   if(triggerSimulator->GetCTA1orCTC1())   SetInput( "AD_CTA1_OR_CTC1" );
   if(triggerSimulator->GetCTA2andCTC2())   SetInput( "AD_CTA2_AND_CTC2" );
   if(triggerSimulator->GetCTA1orCTC1())   SetInput( "AD_CTA1_OR_CTC1" );
   if(triggerSimulator->GetMTAandMTC())   SetInput( "AD_MTA_AND_MTC" );
   if(triggerSimulator->GetMTAorMTC())      SetInput( "AD_MTA_OR_MTC" );
   if(triggerSimulator->GetBBA())         SetInput( "0UBA" );
   if(triggerSimulator->GetBBC())         SetInput( "0UBC" );
   if(triggerSimulator->GetBGAorBGC())      SetInput( "AD_BGA_OR_BGC" );
   if(triggerSimulator->GetBeamGas())      SetInput( "AD_BEAMGAS" );

  return;
}



