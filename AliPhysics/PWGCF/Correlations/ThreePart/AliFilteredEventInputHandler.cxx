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
//     Event handler for AliFilteredEvent input 
//     Author: Paul Baetzing, UIO
//-------------------------------------------------------------------------

#include <TSystem.h>
#include <TTree.h>
#include <TList.h>
#include <TNamed.h>
#include <TFile.h>
#include <TH2.h>

#include "AliFilteredEventInputHandler.h"
#include "AliFilteredEvent.h"

ClassImp(AliFilteredEventInputHandler)

static Option_t *gFilteredEventDataType = "AliFilteredEvent";

//______________________________________________________________________________
AliFilteredEventInputHandler::AliFilteredEventInputHandler() :
    AliInputEventHandler(),
    fEvent(0)
{
  // Default constructor
}

//______________________________________________________________________________
AliFilteredEventInputHandler::AliFilteredEventInputHandler(const char* name, const char* title):
  AliInputEventHandler(name, title),
  fEvent(0)
 {
    // Constructor
}

//______________________________________________________________________________
AliFilteredEventInputHandler::~AliFilteredEventInputHandler() 
{
// Destructor
}

Bool_t AliFilteredEventInputHandler::Init(TTree* tree, Option_t* opt)
{
    // Initialisation necessary for each new tree
    fTree = tree;
    if (!fTree) return kFALSE;
    fTree->GetEntries();
//     SetInactiveBranches(" ");
//     SetActiveBranches("DstTree");
//     
    SwitchOffBranches();
    SwitchOnBranches();
    
    // Get pointer to AOD event
    if (!fEvent) fEvent = new AliFilteredEvent();
    
    tree->SetBranchAddress("Event",&fEvent);
//     fEvent->ReadFromTree(fTree);
    
    return kTRUE;
}
//______________________________________________________________________________
Bool_t AliFilteredEventInputHandler::BeginEvent(Long64_t entry)
{
    // Begin event
    static Int_t prevRunNumber = -1;
    if (prevRunNumber != fEvent->GetRunNumber() ) {
      prevRunNumber = fEvent->GetRunNumber();
    } 
    fTree->GetEvent(entry);
    // When merging, get current event number from GetReadEntry(), 
    // entry gives the events in the current file

    // set transient pointer to event inside tracks
//     fEvent->ConnectTracks();

    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliFilteredEventInputHandler::Notify(const char* path)
{
  // Notifaction of directory change
  fUserInfo=fTree->GetTree()->GetUserInfo();
    
  TTree *ttree = fTree->GetTree();
  if (!ttree) ttree = fTree;
  TString statFname(ttree->GetCurrentFile()->GetName());
  AliInfo(Form("Moving to file %s", statFname.Data()));
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliFilteredEventInputHandler::FinishEvent()
{
  // Finish event
  if (fEvent) fEvent->Reset();
  return kTRUE;
}

//______________________________________________________________________________
Option_t *AliFilteredEventInputHandler::GetDataType() const
{
// Returns handled data type.
   return gFilteredEventDataType;
}

