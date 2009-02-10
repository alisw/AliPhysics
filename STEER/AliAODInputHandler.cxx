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

#include "AliAODInputHandler.h"
#include "AliAODEvent.h"

ClassImp(AliAODInputHandler)

static Option_t *gAODDataType = "AOD";

//______________________________________________________________________________
AliAODInputHandler::AliAODInputHandler() :
    AliInputEventHandler(),
    fEvent(0),
    fFriends(new TList())
{
  // Default constructor
}

//______________________________________________________________________________
AliAODInputHandler::AliAODInputHandler(const char* name, const char* title):
  AliInputEventHandler(name, title),
  fEvent(0),
  fFriends(new TList())
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
    fTree = tree;
    TIter next(fFriends);
    TNamed* obj;
    
    while((obj = (TNamed*)next())) {
	if (fTree->GetTree()) {
	    (fTree->GetTree())->AddFriend("aodTree", obj->GetName());
	} else {
	    fTree->AddFriend("aodTree", obj->GetName());
	}
    }
    
    if (!fTree) return kFALSE;
    // Get pointer to AOD event
    if (fEvent) {
      delete fEvent;
      fEvent = 0;
    }
    fEvent = new AliAODEvent();

    fEvent->ReadFromTree(fTree);
    return kTRUE;
}

Bool_t AliAODInputHandler::BeginEvent(Long64_t /*entry*/)
{
    //
    //if (fTree) fTree->BranchRef();
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
