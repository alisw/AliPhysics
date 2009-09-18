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

/* $Id$ */

//-------------------------------------------------------------------------
//     Event handler for AOD input 
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TTree.h>
#include <TList.h>
#include <TNamed.h>
#include <TFile.h>

#include "AliAODInputHandler.h"
#include "AliAODEvent.h"

ClassImp(AliAODInputHandler)

static Option_t *gAODDataType = "AOD";

//______________________________________________________________________________
AliAODInputHandler::AliAODInputHandler() :
    AliInputEventHandler(),
    fEvent(0),
    fMCEvent(new AliMCEvent()),
    fFriends(new TList()),
    fMergeEvents(kFALSE),
    fFileToMerge(0),
    fTreeToMerge(0),
    fAODEventToMerge(0)
{
  // Default constructor
}

//______________________________________________________________________________
AliAODInputHandler::AliAODInputHandler(const char* name, const char* title):
  AliInputEventHandler(name, title),
  fEvent(0),
  fMCEvent(new AliMCEvent()),
  fFriends(new TList()),
  fMergeEvents(kFALSE),
  fFileToMerge(0),
  fTreeToMerge(0),
  fAODEventToMerge(0)
{
    // Constructor
}

//______________________________________________________________________________
AliAODInputHandler::~AliAODInputHandler() 
{
// Destructor
    fFriends->Delete();
}


Bool_t AliAODInputHandler::Init(TTree* tree, Option_t* /*opt*/)
{
    // Initialisation necessary for each new tree
    if (!fMergeEvents) {
	fTree = tree;
	TIter next(fFriends);
	TNamed* obj;
	
	if (!fTree) return kFALSE;
	fTree->GetEntry(0);
	TString aodTreeFName,aodFriendTreeFName;
	
	while((obj = (TNamed*)next())) {
	    if (fTree->GetTree()) {
		aodTreeFName = (fTree->GetTree()->GetCurrentFile())->GetName();
		aodFriendTreeFName = aodTreeFName;
		aodFriendTreeFName.ReplaceAll("AliAOD.root",obj->GetName());
		aodFriendTreeFName.ReplaceAll("AliAODs.root",obj->GetName());
		(fTree->GetTree())->AddFriend("aodTree", aodFriendTreeFName.Data());
	    } else {
		aodTreeFName = (fTree->GetCurrentFile())->GetName();
		aodFriendTreeFName = aodTreeFName;
		aodFriendTreeFName.ReplaceAll("AliAOD.root",obj->GetName());
		aodFriendTreeFName.ReplaceAll("AliAODs.root",obj->GetName());
		fTree->AddFriend("aodTree", aodFriendTreeFName.Data());
	    }
	}
    } else {
	// Friends have to be merged
	TNamed* filename = (TNamed*) (fFriends->At(0));
	fFileToMerge = new TFile(filename->GetName());
	if (fFileToMerge) {
	    fFileToMerge->GetObject("aodTree", fTreeToMerge);
	    if (!fAODEventToMerge) fAODEventToMerge = new AliAODEvent();
	    fAODEventToMerge->ReadFromTree(fTreeToMerge);
	}
    }
    
    
 

    SwitchOffBranches();
    SwitchOnBranches();
    
    // Get pointer to AOD event
    if (!fEvent) fEvent = new AliAODEvent();

    fEvent->ReadFromTree(fTree);

    return kTRUE;
}

Bool_t AliAODInputHandler::BeginEvent(Long64_t entry)
{
    //
    TClonesArray* mcParticles = (TClonesArray*) (fEvent->FindListObject("mcparticles"));
    if (mcParticles) fMCEvent->SetParticleArray(mcParticles);
    if (fTreeToMerge) fTreeToMerge->GetEntry(entry);
    
    return kTRUE;
}

void AliAODInputHandler::AddFriend(char* filename)
{
    // Add a friend tree 
    TNamed* obj = new TNamed(filename, filename);
    fFriends->Add(obj);
}

Option_t *AliAODInputHandler::GetDataType() const
{
// Returns handled data type.
   return gAODDataType;
}
