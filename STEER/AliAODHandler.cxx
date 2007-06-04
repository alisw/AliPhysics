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
//     Implementation of the Virtual Event Handler Interface for AOD
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TFile.h>
#include <TTree.h>
#include "AliAODHandler.h"
#include "AliAODEvent.h"

ClassImp(AliAODHandler)

//______________________________________________________________________________
AliAODHandler::AliAODHandler() :
    AliVirtualEventHandler(),
    fAODEvent(NULL),
    fAODFile(NULL)
{
  // default constructor
}

//______________________________________________________________________________
AliAODHandler::AliAODHandler(const char* name, const char* title):
    AliVirtualEventHandler(name, title),
    fAODEvent(NULL),
    fAODFile(NULL)
{
}

//______________________________________________________________________________
AliAODHandler::~AliAODHandler() 
{
// destructor
}


Bool_t AliAODHandler::InitIO()
{
    // Initialize IO
    fAODEvent = new AliAODEvent();
    fAODEvent->CreateStdContent();
//    fAODFile  = new TFile("aod.root", "recreate");
    CreateTree();
    
    return kTRUE;
}

Bool_t AliAODHandler::Fill()
{
    // Fill data structures
    printf(">>>>>>>>>>> AliAODHandler::Fill()\n");
    
    FillTree();
    fAODEvent->ClearStd();
    
    return kTRUE;
}

Bool_t AliAODHandler::Terminate()
{
    // Terminate 
    AddAODtoTreeUserInfo();
    return kTRUE;
}

Bool_t AliAODHandler::TerminateIO()
{
    // Terminate IO

//    fAODFile->cd();
//    fTreeA->Write();
//    fAODFile->Close();
    
    return kTRUE;
}


void AliAODHandler::CreateTree()
{
    // Creates the AOD Tree
    fTreeA = new TTree("AOD", "AliAOD tree");
    fTreeA->Branch(fAODEvent->GetList());
}

void AliAODHandler::FillTree()
{
    // Fill the AOD Tree
    fTreeA->Fill();
}


void AliAODHandler::AddAODtoTreeUserInfo()
{
    // Add aod event to tree user info
    fTreeA->GetUserInfo()->Add(fAODEvent);
}
