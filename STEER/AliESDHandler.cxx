
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
	fesdf(NULL),
	fTreeEF(NULL),
	fFileEF(NULL),
	fFileName("AliESDfriends_v1.root"),
	fIsEventSelectedForFriends(kFALSE)
{

	// default constructor
}

//______________________________________________________________________________
AliESDHandler::AliESDHandler(const char* name, const char* title):
	AliVEventHandler(name, title),
	fesdf(NULL),
	fTreeEF(NULL),
	fFileEF(NULL),
	fFileName("AliESDfriends_v1.root"),
	fIsEventSelectedForFriends(kFALSE)
{

	// constructor with name and title

}

//______________________________________________________________________________
AliESDHandler::~AliESDHandler() 
{
	// Destructor.
	delete fesdf;
	if(fFileEF){
		// is already handled in TerminateIO
		fFileEF->Close();
		delete fFileEF;
	}
	delete fTreeEF;
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

	fesdf = new AliESDfriend();

	// Open the file with friends
	if (option.Contains("proof")) {
		// proof
		// Merging via files. Need to access analysis manager via interpreter.
		gROOT->ProcessLine(Form("AliAnalysisManager::GetAnalysisManager()->OpenProofFile(\"%s\", \"RECREATE\");", fFileName.Data()));
		gROOT->ProcessLine(Form("AliAnalysisManager::GetAnalysisManager()->GetCommonOutputContainer()->SetFile((TFile*)0x%lx);", gFile));
		fFileEF = gFile;
	} else {
		// local and grid
		fFileEF = new TFile(fFileName.Data(), "RECREATE");
	}

	// Create the friends tree
	fFileEF->cd();
	fTreeEF = new TTree("esdFriendTree", "Tree with ESD friends");
      	fTreeEF->Branch("ESDfriend.","AliESDfriend", &fesdf);

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
	if (fesdf) fesdf->~AliESDfriend();
	new(fesdf) AliESDfriend();  
	return kTRUE;
}

//______________________________________________________________________________
Bool_t AliESDHandler::Terminate()
{
	//
	// Terminate 
	//

	return kTRUE;
}

//______________________________________________________________________________
Bool_t AliESDHandler::TerminateIO()
{
	//
	// Terminate IO
	//

	if (fFileEF) {
		fFileEF->cd();
		fTreeEF->Write();
		fFileEF->Close();
		delete fFileEF;
		fFileEF = 0;
	}

	return kTRUE;
}

//______________________________________________________________________________
void AliESDHandler::FillTree()
{
	//
	// Fill the ESD Tree
	//
	if (fIsEventSelectedForFriends){
		AliDebug(2,Form("number of friend tracks = %d\n",fesdf->GetNumberOfTracks()));
	}
	else {
		fesdf->SetSkipBit(kTRUE);
	}
	AliDebug(2,Form("friend = %p",fesdf));
	fFileEF->cd();
	fTreeEF->Fill();
}
