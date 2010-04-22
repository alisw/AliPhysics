/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/*$Id$*/

/////////////////////////////////////////////////////////////
//
//   Test task to add an object to the new AliESDfriends file
//
// /////////////////////////////////////////////////////////////

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>

#include "AliLog.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskAddObject.h"


ClassImp(AliAnalysisTaskAddObject)


//________________________________________________________________________
AliAnalysisTaskAddObject::AliAnalysisTaskAddObject():
AliAnalysisTask(),
fESDInput(0x0),
fESDfriendInput(0x0),
fESDhandler(0x0),
fh(0x0)
{
	// Dummy Constructor
	
}

//________________________________________________________________________
AliAnalysisTaskAddObject::AliAnalysisTaskAddObject(const char* name):
AliAnalysisTask(name,"Adding an object"),
fESDInput(0),
fESDfriendInput(0),
fESDhandler(0x0),
fh(0x0)
{
	// Constructor
	
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TTree
	//	DefineOutput(0,TTree::Class());  
	// Output slot #1 writes into a TH1D
	DefineOutput(0,TH1D::Class());  
}

//________________________________________________________________________
AliAnalysisTaskAddObject::~AliAnalysisTaskAddObject()
{

	// dtor
	if (fh){
		delete fh;
		fh = 0x0;
	}
}  

//______________________________________________________________________________
void AliAnalysisTaskAddObject::ConnectInputData(Option_t* /*option*/)
{
	//
	// Connect the input data
	//

	printf("AliAnalysisTaskAddObject::ConnectInputData()\n");
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	if (!mgr) AliFatal("No analysis manager available");
	fESDhandler = dynamic_cast<AliESDInputHandler *>(mgr->GetInputEventHandler());
    
	if (fESDhandler) {
		fESDInput = fESDhandler->GetEvent();
	} else {
		AliFatal("No ESD input event handler connected") ; 
	}
}
//________________________________________________________________________
void AliAnalysisTaskAddObject::CreateOutputObjects()
{
	//
	// Create the output container
	//
	//OpenFile(0,"UPDATE");
	fh = new TH1D("fh1","Integrated Length",100,0,1000);
	return;
}

//________________________________________________________________________
void AliAnalysisTaskAddObject::Exec(Option_t */*option*/)
{

	//	if (fDebug > 1) {
	Long_t entry = fESDhandler->GetReadEntry();
	AliDebug(2,Form("AliAnalysisTaskAddObject::Exec() %s ==> processing event %lld", fESDhandler->GetTree()->GetCurrentFile()->GetName(),entry));
	//}  
	fESDInput = fESDhandler->GetEvent();
	if(!fESDInput) {
		printf("AliAnalysisTaskAddObject::Exec(): no ESD \n");
		return;
	} 
	for (Int_t i = 0; i< fESDInput->GetNumberOfTracks(); i++){
		AliESDtrack* t = fESDInput->GetTrack(i);
		Double_t l = t->GetIntegratedLength();
		fh->Fill(l);
	}
	PostData(0,fh);
	return;
}

//________________________________________________________________________
void AliAnalysisTaskAddObject::Terminate(Option_t */*option*/)
{
	// Terminate analysis
	//
	AliDebug(2,"AliAnalysisTaskAddObject: Terminate() \n");
	
	return;
}
