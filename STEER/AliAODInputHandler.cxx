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

#include <TSystem.h>
#include <TTree.h>
#include <TList.h>
#include <TNamed.h>
#include <TFile.h>
#include <TH2.h>

#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliVCuts.h"
#include "AliMCEvent.h"

ClassImp(AliAODInputHandler)

static Option_t *gAODDataType = "AOD";

//______________________________________________________________________________
AliAODInputHandler::AliAODInputHandler() :
    AliInputEventHandler(),
    fEvent(0),
    fMCEvent(new AliMCEvent()),
    fFriends(new TList()),
    fMergeEvents(kFALSE),
    fFriendsConnected(kFALSE),
    fFileToMerge(0),
    fTreeToMerge(0),
    fAODEventToMerge(0),
    fMergeOffset(0)
{
  // Default constructor
  fHistStatistics[0] = fHistStatistics[1] = NULL;
}

//______________________________________________________________________________
AliAODInputHandler::AliAODInputHandler(const char* name, const char* title):
  AliInputEventHandler(name, title),
  fEvent(0),
  fMCEvent(new AliMCEvent()),
  fFriends(new TList()),
  fMergeEvents(kFALSE),
  fFriendsConnected(kFALSE),
  fFileToMerge(0),
  fTreeToMerge(0),
  fAODEventToMerge(0),
  fMergeOffset(0)
{
    // Constructor
  fHistStatistics[0] = fHistStatistics[1] = NULL;
}

//______________________________________________________________________________
AliAODInputHandler::~AliAODInputHandler() 
{
// Destructor
  fFriends->Delete();
  if (fHistStatistics[0]) {
    delete fHistStatistics[0];
    fHistStatistics[0] = 0;
  }
  if (fHistStatistics[1]) {
    delete fHistStatistics[1];
    fHistStatistics[1] = 0;
  }
}

//______________________________________________________________________________
Bool_t AliAODInputHandler::Init(TTree* tree, Option_t* opt)
{
    // Initialisation necessary for each new tree
    fTree = tree;
    if (!fTree) return kFALSE;
    fTree->GetEntries();
    ConnectFriends();

    SwitchOffBranches();
    SwitchOnBranches();
    
    // Get pointer to AOD event
    if (!fEvent) fEvent = new AliAODEvent();
    
    fEvent->ReadFromTree(fTree);
    
    if (fMixingHandler) fMixingHandler->Init(tree, opt);
    
    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODInputHandler::BeginEvent(Long64_t entry)
{
    // Begin event
    TClonesArray* mcParticles = (TClonesArray*) (fEvent->FindListObject("mcparticles"));
    if (mcParticles) fMCEvent->SetParticleArray(mcParticles);
    if (fTreeToMerge) fTreeToMerge->GetEntry(entry + fMergeOffset);
    
    fIsSelectedResult = fEvent->GetHeader()->GetOfflineTrigger();

    if (fMixingHandler) fMixingHandler->BeginEvent(entry);
    
    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODInputHandler::Notify(const char* path)
{
  // Notifaction of directory change
  if (fMixingHandler) fMixingHandler->Notify(path);
  if (!fFriendsConnected) {
      ConnectFriends();
      fEvent->ReadFromTree(fTree, "reconnect");
  }
  fFriendsConnected = kFALSE;
    
  TTree *ttree = fTree->GetTree();
  if (!ttree) ttree = fTree;
  TString statFname(gSystem->DirName(ttree->GetCurrentFile()->GetName()));
  statFname += "/EventStat_temp.root";
  TFile *statFile = TFile::Open(statFname, "READ");
  if (statFile) {
     TList *list = (TList*)statFile->Get("cstatsout");
     if (list) {
        AliVCuts *physSel = (AliVCuts*)list->At(0);
        if (physSel) {
           TH2F *hAll = dynamic_cast<TH2F*>(physSel->GetStatistics("ALL"));
           TH2F *hBin0 = dynamic_cast<TH2F*>(physSel->GetStatistics("BIN0"));
           if (fHistStatistics[0] && hAll) {
              TList tmplist;
              tmplist.Add(hAll);
              fHistStatistics[0]->Merge(&tmplist);
              tmplist.Clear();
              tmplist.Add(hBin0);
              if (fHistStatistics[1] && hBin0) fHistStatistics[1]->Merge(&tmplist);
           } else {
             fHistStatistics[0] = static_cast<TH2F*>(hAll->Clone());
             fHistStatistics[1] = static_cast<TH2F*>(hBin0->Clone());
             fHistStatistics[0]->SetDirectory(0);
             fHistStatistics[1]->SetDirectory(0);
           }   
        }
        delete list;
     }
     delete statFile;
  }
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliAODInputHandler::FinishEvent()
{
  // Finish event
  if (fMixingHandler) fMixingHandler->FinishEvent();
  return kTRUE;
}

//______________________________________________________________________________
void AliAODInputHandler::AddFriend(char* filename)
{
    // Add a friend tree 
    TNamed* obj = new TNamed(filename, filename);
    fFriends->Add(obj);
}

//______________________________________________________________________________
Option_t *AliAODInputHandler::GetDataType() const
{
// Returns handled data type.
   return gAODDataType;
}

//______________________________________________________________________________
TObject *AliAODInputHandler::GetStatistics(Option_t *option) const
{
// Get the statistics histogram(s) from the physics selection object. This
// should be called during FinishTaskOutput(). Option can be empty (default
// statistics histogram) or BIN0.
   TString opt(option);
   opt.ToUpper();
   if (opt=="BIN0") return fHistStatistics[1];
   return fHistStatistics[0];
}   

void AliAODInputHandler::ConnectFriends()
{
    // Connect the friend trees 
    if (!fMergeEvents) {
	TIter next(fFriends);
	TNamed* obj;
	TString aodTreeFName,aodFriendTreeFName;
	TTree *ttree = fTree->GetTree();
	if (!ttree) ttree = fTree;
	aodTreeFName = ttree->GetCurrentFile()->GetName();
	
	while((obj = (TNamed*)next())) {
	    aodFriendTreeFName = aodTreeFName;
	    aodFriendTreeFName.ReplaceAll("AliAOD.root",obj->GetName());
	    aodFriendTreeFName.ReplaceAll("AliAODs.root",obj->GetName());
	    ttree->AddFriend("aodTree", aodFriendTreeFName.Data());
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
    fFriendsConnected = kTRUE;
}

