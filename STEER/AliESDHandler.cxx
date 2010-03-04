
/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-------------------------------------------------------------------------
//
//     Implementation of the Virtual Event Handler Interface for ESD
//
//-------------------------------------------------------------------------


#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliESDHandler.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"


ClassImp(AliESDHandler)

//______________________________________________________________________________
AliESDHandler::AliESDHandler() :
	AliVEventHandler(),
	fESDEvent(NULL),
	fesdf(NULL),
	fTreeE(NULL),
	fFileE(NULL),
	fFileEF(NULL),
	fFileName("")
{
	// default constructor
}

//______________________________________________________________________________
AliESDHandler::AliESDHandler(const char* name, const char* title):
	AliVEventHandler(name, title),
	fESDEvent(NULL),
	fesdf(NULL),
	fTreeE(NULL),
	fFileE(NULL),
	fFileEF(NULL),
	fFileName("")
{

	// constructor with name and title

}

//______________________________________________________________________________
AliESDHandler::~AliESDHandler() 
{
	// Destructor.
	delete fESDEvent;
	delete fesdf;
	if(fFileE){
		// is already handled in TerminateIO
		fFileE->Close();
		delete fFileE;
	}
	if(fFileEF){
		// is already handled in TerminateIO
		fFileEF->Close();
		delete fFileEF;
	}
	delete fTreeE;
}

//______________________________________________________________________________
Bool_t AliESDHandler::Init(Option_t* opt)
{
	//
	// Initialize IO
	//
	
	// File opening according to execution mode
	TString option(opt);
	option.ToLower();
	TDirectory *owd = gDirectory;
	if (option.Contains("proof")) {
		// proof
		// Merging via files. Need to access analysis manager via interpreter.
		gROOT->ProcessLine(Form("AliAnalysisManager::GetAnalysisManager()->OpenProofFile(\"%s\", \"RECREATE\");", fFileName.Data()));
		gROOT->ProcessLine(Form("AliAnalysisManager::GetAnalysisManager()->GetCommonOutputContainer()->SetFile((TFile*)0x%lx);", gFile));
		fFileE = gFile;
	} else {
		// local and grid
		fFileE = new TFile(fFileName.Data(), "RECREATE");
	}
	CreateTree(1);
	CreateFriends(1);
	owd->cd();
	
	return kTRUE;
}


//______________________________________________________________________________
Bool_t AliESDHandler::FinishEvent()
{
	//
	// Fill the tree 
	//

	FillTree();
	
	// resetting
	fESDEvent->Reset();
	fesdf->~AliESDfriend();
	new(fesdf) AliESDfriend();  
	return kTRUE;
}

//______________________________________________________________________________
Bool_t AliESDHandler::Terminate()
{
	//
	// Terminate 
	//

	AddESDtoTreeUserInfo();
	return kTRUE;
}

//______________________________________________________________________________
Bool_t AliESDHandler::TerminateIO()
{
	//
	// Terminate IO
	//

	if (fFileE) {
		fFileE->cd();
		fTreeE->Write();
		fFileE->Close();
		delete fFileE;
		fFileE = 0;
	}

	return kTRUE;
}


//______________________________________________________________________________
void AliESDHandler::CreateTree(Int_t /*flag*/)
{
	//
	// Creates the ESD Tree
	// 

	fTreeE = new TTree("esdTree", "AliESD tree");
	// Create the ESDevent object
	if(!fESDEvent){
		fESDEvent = new AliESDEvent();
		fESDEvent->CreateStdContent();
	}
	fESDEvent->WriteToTree(fTreeE);
}
//______________________________________________________________________________
void AliESDHandler::CreateFriends(Int_t /*flag*/)
{
	fesdf = new AliESDfriend();

      	TBranch *br=fTreeE->Branch("ESDfriend.","AliESDfriend", &fesdf);
	br->SetFile("AliESDfriends_v1.root");
	fESDEvent->AddObject(fesdf);
}

//______________________________________________________________________________
void AliESDHandler::FillTree()
{
	//
	// Fill the ESD Tree
	//

	AliDebug(2,Form("number of friend tracks = %d\n",fesdf->GetNumberOfTracks()));

	fFileE->cd();
	fTreeE->Fill();
}

//______________________________________________________________________________
void AliESDHandler::AddESDtoTreeUserInfo()
{
	//
	// Add aod event to tree user info
	//

	fTreeE->GetUserInfo()->Add(fESDEvent);
}



