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
// ---------------------
// Class AliVZEROTrigger
// ---------------------
// Top class to simulate the VZERO trigger response
// This class is only used for interface with AliTriggerDetector
// Its create and Set  Inputs of the CTP
// The Calculation of the trigger response is done into AliVZEROTriggerSimulator
//


#include <TClonesArray.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTriggerInput.h"

#include "AliVZEROdigit.h"
#include "AliVZEROTriggerSimulator.h"
#include "AliVZEROTrigger.h"

//______________________________________________________________________
ClassImp(AliVZEROTrigger)

//______________________________________________________________________

AliVZEROTrigger::AliVZEROTrigger():AliTriggerDetector()
{
   SetName("VZERO");
   CreateInputs();
}
//______________________________________________________________________
void AliVZEROTrigger::CreateInputs()
{
	// Do not create inputs again!!
	if( fInputs.GetEntriesFast() > 0 ) return;

	fInputs.AddLast( new AliTriggerInput( "VZERO_BBA_AND_BBC", "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_BBA_OR_BBC","VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_BGA_AND_BBC",  "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "0VGA",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_BGC_AND_BBA", "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "0VGC",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_CTA1_AND_CTC1",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_CTA1_OR_CTC1",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_CTA2_AND_CTC2",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_CTA2_OR_CTC2",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_MTA_AND_MTC",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_MTA_OR_MTC",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "0VBA",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "0VBC",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_BGA_OR_BGC",   "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_BEAMGAS",   "VZERO", 0 ) );

	// The following are kept for compatibility with the CTP configuration file. Will have to be removed at some point
	fInputs.AddLast( new AliTriggerInput( "VZERO_AND", "VZERO", 0 ) );
	fInputs.AddLast( new AliTriggerInput( "VZERO_OR","VZERO", 0 ) );
}

//______________________________________________________________________
void AliVZEROTrigger::Trigger()
{
  //  ********** Get run loader for the current event **********
	AliRunLoader* runLoader = AliRunLoader::Instance();

	AliLoader* loader = runLoader->GetLoader( "VZEROLoader" );

	if(!loader) {
		AliError("Can not get VZERO loader");
		return;
	}
	loader->LoadDigits("update");
	TTree* vzeroDigitsTree = loader->TreeD();

	if (!vzeroDigitsTree) {
		AliError("Can not get the VZERO digit tree");
		return;
	}
	TClonesArray* vzeroDigits = NULL;
	TBranch* digitBranch = vzeroDigitsTree->GetBranch("VZERODigit");
	digitBranch->SetAddress(&vzeroDigits);

	AliVZEROTriggerSimulator * triggerSimulator = new AliVZEROTriggerSimulator(vzeroDigitsTree,vzeroDigits);
	

	triggerSimulator->Run();
	
	loader->WriteDigits("OVERWRITE");  
	loader->UnloadDigits();     

	if(triggerSimulator->GetBBAandBBC())	SetInput( "VZERO_BBA_AND_BBC" );
	if(triggerSimulator->GetBBAorBBC())		SetInput( "VZERO_BBA_OR_BBC" );
	if(triggerSimulator->GetBGAandBBC())	SetInput( "VZERO_BGA_AND_BBC" );
	if(triggerSimulator->GetBGA())			SetInput( "0VGA" );
	if(triggerSimulator->GetBGCandBBA())	SetInput( "VZERO_BGC_AND_BBA" );
	if(triggerSimulator->GetBGC())			SetInput( "0VGC" );
	if(triggerSimulator->GetCTA1andCTC1())	SetInput( "VZERO_CTA1_AND_CTC1" );
	if(triggerSimulator->GetCTA1orCTC1())	SetInput( "VZERO_CTA1_OR_CTC1" );
	if(triggerSimulator->GetCTA2andCTC2())	SetInput( "VZERO_CTA2_AND_CTC2" );
	if(triggerSimulator->GetCTA1orCTC1())	SetInput( "VZERO_CTA1_OR_CTC1" );
 	if(triggerSimulator->GetMTAandMTC())	SetInput( "VZERO_MTA_AND_MTC" );
 	if(triggerSimulator->GetMTAorMTC())		SetInput( "VZERO_MTA_OR_MTC" );
 	if(triggerSimulator->GetBBA())			SetInput( "0VBA" );
 	if(triggerSimulator->GetBBC())			SetInput( "0VBC" );
 	if(triggerSimulator->GetBGAorBGC())		SetInput( "VZERO_BGA_OR_BGC" );
 	if(triggerSimulator->GetBeamGas())		SetInput( "VZERO_BEAMGAS" );

	// The following are kept for compatibility with the CTP configuration file. Will have to be removed at some point
	if(triggerSimulator->GetBBAandBBC())    SetInput( "VZERO_AND" );
	if(triggerSimulator->GetBBAorBBC())             SetInput( "VZERO_OR" );

  return;
}



