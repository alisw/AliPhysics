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

/////////////////////////////////////////////////////////////////////////
//
// Task to Copy ESDs
//
////////////////////////////////////////////////////////////////////////


#include <TTree.h>

#include "AliAnalysisTaskFilter.h"
#include "AliAnalysisTaskCopyESD.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliESDfriend.h"

ClassImp(AliAnalysisTaskCopyESD)

////////////////////////////////////////////////////////////////////////

AliAnalysisTaskCopyESD::AliAnalysisTaskCopyESD():
	AliAnalysisTaskFilter(),
	fESDEvent(0x0),
	fESDfriend(0x0)
{
	// Default constructor
	//DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TTree
	//DefineOutput(0,TTree::Class());  //My private output
}
//-------------------------------------------------------------------------------
AliAnalysisTaskCopyESD::AliAnalysisTaskCopyESD(const char* name):
	AliAnalysisTaskFilter(name),
	fESDEvent(0x0),
	fESDfriend(0x0)
{
	// Constructor
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskCopyESD::UserCreateOutputObjects()
{
	// Create the output container
	AliInfo("In UserCreateOuptputObject");
}

//-------------------------------------------------------------------------------
void AliAnalysisTaskCopyESD::Init()
{
	// Initialization
	if (fDebug > 1) AliInfo("Init() \n");
}


//-------------------------------------------------------------------------------
void AliAnalysisTaskCopyESD::UserExec(Option_t */*option*/)
{
	// Execute analysis for current event
	//
	
	AliInfo("Copying event");	
	AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
	fESDEvent = ESDEvent(); // get the output ESD
	fESDfriend = ESDfriend();  // get the output friend
	esd->Copy(*fESDEvent);

	// releasing tracks to avoid copying them 
	//	for (Int_t i = 0; i<fESDEvent->GetNumberOfTracks(); i++){
	//	fESDEvent->GetTrack(i)->ReleaseESDfriendTrack();
	//}
}


//-------------------------------------------------------------------------------
void AliAnalysisTaskCopyESD::Terminate(Option_t */*option*/)
{
	// Terminate analysis
	//
	if (fDebug > 1) printf("AnalysisCopyESD: Terminate() \n");
}

