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
//     Event handler for ESD input 
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TArchiveFile.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjString.h>
#include <TProcessID.h>
#include <TMap.h>

#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliVCuts.h"
#include "AliESD.h"
#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliLog.h"

ClassImp(AliESDInputHandler)

static Option_t *gESDDataType = "ESD";

//______________________________________________________________________________
AliESDInputHandler::AliESDInputHandler() :
  AliInputEventHandler(),
  fEvent(0x0),
  fFriend(0x0),
  fESDpid(0x0),
  fAnalysisType(0),
  fNEvents(0),
  fHLTEvent(0x0),
  fHLTTree(0x0),
  fUseHLT(kFALSE),
  fTagCutSumm(0x0),
  fUseTags(kFALSE),
  fChainT(0),
  fTreeT(0),
  fRunTag(0),
  fEventTag(0),
  fReadFriends(0),
  fFriendFileName("AliESDfriends.root")
{
  // default constructor
}

//______________________________________________________________________________
AliESDInputHandler::~AliESDInputHandler() 
{
  //  destructor
  if (fRunTag) delete fRunTag;
}

//______________________________________________________________________________
AliESDInputHandler::AliESDInputHandler(const char* name, const char* title):
    AliInputEventHandler(name, title), fEvent(0x0), fFriend(0x0), fESDpid(0x0), fAnalysisType(0),
    fNEvents(0),  fHLTEvent(0x0), fHLTTree(0x0), fUseHLT(kFALSE), fTagCutSumm(0x0), fUseTags(kFALSE), fChainT(0), fTreeT(0), fRunTag(0), fEventTag(0), fReadFriends(0), fFriendFileName("AliESDfriends.root")
{
    // Constructor
}

//______________________________________________________________________________
Bool_t AliESDInputHandler::Init(TTree* tree,  Option_t* opt)
{
    //
    // Initialisation necessary for each new tree 
    // 
    fAnalysisType = opt;
    fTree = tree;
    
    if (!fTree) return kFALSE;
//    fTree->GetEntry(0);
    

    if (!fEvent) fEvent = new AliESDEvent();
    fEvent->ReadFromTree(fTree);
    fNEvents = fTree->GetEntries();

    if (fMixingHandler) fMixingHandler->Init(tree,  opt);

    return kTRUE;
}

//______________________________________________________________________________
Bool_t AliESDInputHandler::BeginEvent(Long64_t entry)
{
    
    // Copy from old to new format if necessary
  static Bool_t called = kFALSE;
  if (!called && fEventCuts && IsUserCallSelectionMask())
     AliInfo(Form("The ESD input handler expects that the first task calls AliESDInputHandler::CheckSelectionMask() %s", fEventCuts->ClassName()));
  AliESD* old = ((AliESDEvent*) fEvent)->GetAliESDOld();
  if (old) {
   ((AliESDEvent*)fEvent)->CopyFromOldESD();
   old->Reset();
  }

  if (fHLTTree) {
      fHLTTree->GetEntry(entry);
  }
  
  fNewEvent = kTRUE;
  //
  // Event selection
  // 
  fIsSelectedResult = 0;
  if (fEventCuts && !IsUserCallSelectionMask())
      fIsSelectedResult = fEventCuts->GetSelectionMask((AliESDEvent*)fEvent); 
  //
  // Friends
  ((AliESDEvent*)fEvent)->SetESDfriend(fFriend);
  called = kTRUE;

  if (fMixingHandler) fMixingHandler->BeginEvent(entry);
  if (fUseTags && fRunTag) {
    fEventTag = 0;
    if (entry >= fRunTag->GetNEvents()) {
      AliError(Form("Current event %d does not match max range from run tag: 0-%d", (Int_t)entry, fRunTag->GetNEvents()));
      return kTRUE;
    }
    fEventTag = fRunTag->GetEventTag(entry);   
  }      
  return kTRUE;
}

//______________________________________________________________________________
void AliESDInputHandler::CheckSelectionMask()
{
// This method can be called by a task only if IsUserCallSelectionMask is true.
   if (!fEventCuts || !IsUserCallSelectionMask()) return;
   fIsSelectedResult = fEventCuts->GetSelectionMask((AliESDEvent*)fEvent);
}
   
//______________________________________________________________________________
Bool_t  AliESDInputHandler::FinishEvent()
{
    // Finish the event 
  if(fEvent)fEvent->Reset();
  if (fMixingHandler) fMixingHandler->FinishEvent();
  return kTRUE;
} 

//______________________________________________________________________________
Bool_t AliESDInputHandler::Notify(const char* path)
{
  // Notify a directory change
  static Bool_t firsttime = kFALSE;
  AliInfo(Form("Directory change %s \n", path));
  //
  // Handle the friends first
  //
  if (!fTree->FindBranch("ESDfriend.") && fReadFriends) {
    // Try to add ESDfriend. branch as friend
    TString esdTreeFName, esdFriendTreeFName;    
    esdTreeFName = (fTree->GetCurrentFile())->GetName();
    esdFriendTreeFName = esdTreeFName;
    esdFriendTreeFName.ReplaceAll("AliESDs.root", fFriendFileName.Data());
    TTree* cTree = fTree->GetTree();
    if (!cTree) cTree = fTree;      
    cTree->AddFriend("esdFriendTree", esdFriendTreeFName.Data());
    cTree->SetBranchStatus("ESDfriend.", 1);
    fFriend = (AliESDfriend*)(fEvent->FindListObject("AliESDfriend"));
    if (fFriend) cTree->SetBranchAddress("ESDfriend.", &fFriend);
  } 
  //
  //
  SwitchOffBranches();
  SwitchOnBranches();
  fFriend = (AliESDfriend*)(fEvent->FindListObject("AliESDfriend"));
  //
  if (fUseHLT) {
    // Get HLTesdTree from current file
    TTree* cTree = fTree;
    if (fTree->GetTree()) cTree = fTree->GetTree();
    TFile* cFile = cTree->GetCurrentFile();
    cFile->GetObject("HLTesdTree", fHLTTree);
	
    if (fHLTTree) {
	   if (!fHLTEvent) fHLTEvent = new AliESDEvent();
	   fHLTEvent->ReadFromTree(fHLTTree);
    }
  }

  if (!fUseTags) {
    if (fMixingHandler) fMixingHandler->Notify(path);
    return kTRUE;
  }
    
  Bool_t zip = kFALSE;
    
  // Setup the base path
  TString pathName(path);
  Int_t index = pathName.Index("#");
  if (index>=0) {
    zip = kTRUE;
    pathName = pathName(0,index);
  } else {
    pathName = gSystem->DirName(pathName);
  }
  if (fTree->GetCurrentFile()->GetArchive()) zip = kTRUE;
  if (pathName.IsNull()) pathName = "./";  
  printf("AliESDInputHandler::Notify() Path: %s\n", pathName.Data());

  if (fRunTag) {
    fRunTag->Clear();
  } else {
    fRunTag = new AliRunTag();
  }
    
  const char* tagPattern = "ESD.tag.root";
  TString sname;
  TString tagFilename;
  if (zip) {
    TObjArray* arr = fTree->GetCurrentFile()->GetArchive()->GetMembers();
    TIter next(arr);
    TObject *objarchive;
    while ((objarchive = next())) {
      sname = objarchive->GetName();
	   if (sname.Contains(tagPattern)) { 
        tagFilename = pathName;
        if (index>=0) tagFilename += "#";
        else tagFilename += "/";
        tagFilename += sname;
        AliInfo(Form("Tag file found: %s\n", tagFilename.Data()));
        break; // There should be only one such file in the archive
      }//pattern check
    } // archive file loop
  } else {
    void * dirp = gSystem->OpenDirectory(pathName.Data());
    while(1) {
      sname = gSystem->GetDirEntry(dirp);
      if (sname.IsNull()) break;
      if (sname.Contains(tagPattern)) { 
        tagFilename = pathName;
        tagFilename += "/";
        tagFilename += sname;
        AliInfo(Form("Tag file found: %s\n", tagFilename.Data()));
        break;
      }//pattern check
    }//directory loop
    gSystem->FreeDirectory(dirp);
  }
  if (tagFilename.IsNull()) {
    if (firsttime) AliWarning(Form("Tag file not found in directory: %s", pathName.Data()));
    firsttime = kFALSE;
    delete fRunTag; fRunTag = 0;
    return kTRUE;
  }
  TFile *tagfile = TFile::Open(tagFilename);
  if (!tagfile) {
    AliError(Form("Cannot open tag file: %s", tagFilename.Data()));
    delete fRunTag; fRunTag = 0;
    return kTRUE;
  }   
  fTreeT = (TTree*)tagfile->Get("T"); // file is the owner
  if (!fTreeT) {
    AliError(Form("Cannot get tree of tags from file: %s", tagFilename.Data()));
    delete fRunTag; fRunTag = 0;
    return kTRUE;
  }
    
  fTreeT->SetBranchAddress("AliTAG",&fRunTag);
  fTreeT->GetEntry(0);
  delete tagfile;
  // Notify the mixing handler after the tags are loaded
  if (fMixingHandler) fMixingHandler->Notify(path);
  return kTRUE;
}

//______________________________________________________________________________
Option_t *AliESDInputHandler::GetDataType() const
{
// Returns handled data type.
   return gESDDataType;
}

//______________________________________________________________________________
Int_t AliESDInputHandler::GetNEventAcceptedInFile()
{
  // Get number of events in file accepted by the tag cuts
  // return -1 if no info is available
  if (!fTagCutSumm) {
    TList *luo = fTree->GetUserInfo();
    if (!luo) {
      AliInfo(Form("No user info in input tree - no tag cut summary\n"));
      return -1;
    }
    for (int iluo=0; iluo<luo->GetEntries(); iluo++) {
      fTagCutSumm = dynamic_cast<TMap *>(luo->At(iluo));
      if (fTagCutSumm) break;
    }
    if (!fTagCutSumm) {
      AliInfo(Form("No tag summary map in input tree\n"));
      return -1;
    }
  }

  TObjString *ostr = 0;
  if (fTagCutSumm->FindObject(fTree->GetCurrentFile()->GetName()))
    ostr = (TObjString *) fTagCutSumm->GetValue(fTree->GetCurrentFile()->GetName());
  else {
    AliInfo(Form("No tag cut summary for file %s\n", fTree->GetCurrentFile()->GetName()));
    return -1;
  }
  char *iTagInfo;
  iTagInfo = strdup(ostr->GetString().Data());

  Int_t iAcc = atoi(strtok(iTagInfo, ","));
  
  AliInfo(Form("Got %i accepted events for file %s", iAcc,  fTree->GetCurrentFile()->GetName()));
  
  free(iTagInfo);

  return iAcc;
}

//______________________________________________________________________________
Int_t AliESDInputHandler::GetNEventRejectedInFile()
{
  // Get number of events in file rejected by the tag cuts
  // return -1 if no info is available
  if (!fTagCutSumm) {
    TList *luo = fTree->GetUserInfo();
    if (!luo) {
      AliInfo(Form("No user info in input tree - no tag cut summary\n"));
      return -1;
    }
    for (int iluo=0; iluo<luo->GetEntries(); iluo++) {
      fTagCutSumm = dynamic_cast<TMap *>(luo->At(iluo));
      if (fTagCutSumm) break;
    }
    if (!fTagCutSumm) {
      AliInfo(Form("No tag summary map in input tree\n"));
      return -1;
    }
  }

  TObjString *ostr = 0;
  if (fTagCutSumm->FindObject(fTree->GetCurrentFile()->GetName()))
    ostr = (TObjString *) fTagCutSumm->GetValue(fTree->GetCurrentFile()->GetName());
  else {
    AliInfo(Form("No tag cut summary for file %s\n", fTree->GetCurrentFile()->GetName()));
    return -1;
  }
  char *iTagInfo;
  iTagInfo = strdup(ostr->GetString().Data());

  strtok(iTagInfo, ",");
  Int_t iRej = atoi(strtok(NULL, ","));
  
  AliInfo(Form("Got %i accepted events for file %s", iRej,  fTree->GetCurrentFile()->GetName()));
  
  free(iTagInfo);

  return iRej;
}

//______________________________________________________________________________
Bool_t AliESDInputHandler::GetCutSummaryForChain(Int_t *aTotal,  Int_t *aAccepted,  Int_t *aRejected)
{
  // Get number of events in the full chain
  // Count accepted and rejected events
  // return kFALSE if no info is available
  if (!fTagCutSumm) {
    TList *luo = fTree->GetUserInfo();
    if (!luo) {
      AliInfo(Form("No user info in input tree - no tag cut summary\n"));
      return kFALSE;
    }
    for (int iluo=0; iluo<luo->GetEntries(); iluo++) {
      fTagCutSumm = dynamic_cast<TMap *>(luo->At(iluo));
      if (fTagCutSumm) break;
    }
    if (!fTagCutSumm) {
      AliInfo(Form("No tag summary map in input tree\n"));
      return kFALSE;
    }
  }
  
  TMapIter *tIter = new TMapIter(fTagCutSumm);
  
  Int_t iTotList=0, iAccList=0, iRejList=0;

  TObject *cobj;
  while ((cobj = tIter->Next())) {
    TObjString *kstr = (TObjString *) cobj;
    TObjString *vstr = (TObjString *) fTagCutSumm->GetValue(kstr->GetString().Data());
    //    printf("Got object value %s %s\n", kstr->GetString().Data(), vstr->GetString().Data());
    char *iTagInfo;
    iTagInfo = strdup(vstr->GetString().Data());
    
    Int_t iAcc = atoi(strtok(iTagInfo, ","));
    Int_t iRej = atoi(strtok(NULL, ","));
    
    iAccList += iAcc;
    iRejList += iRej;
    iTotList += (iAcc+iRej);
  }

  *aTotal = iTotList;
  *aAccepted = iAccList;
  *aRejected = iRejList;

  return kTRUE;
}

//______________________________________________________________________________
Int_t AliESDInputHandler::GetNFilesEmpty()
{
  // Count number of files in which all events were de-selected
  // For such files Notify() will NOT be called
  // return -1 if no info is available
  if (!fTagCutSumm) {
    TList *luo = fTree->GetUserInfo();
    if (!luo) {
      AliInfo(Form("No user info in input tree - no tag cut summary\n"));
      return -1;
    }
    for (int iluo=0; iluo<luo->GetEntries(); iluo++) {
      fTagCutSumm = dynamic_cast<TMap *>(luo->At(iluo));
      if (fTagCutSumm) break;
    }
    if (!fTagCutSumm) {
      AliInfo(Form("No tag summary map in input tree\n"));
      return -1;
    }
  }
  
  TMapIter *tIter = new TMapIter(fTagCutSumm);
  
  Int_t iFilesEmpty = 0;

  TObject *cobj;
  while ((cobj = tIter->Next())) {
    TObjString *kstr = (TObjString *) cobj;
    TObjString *vstr = (TObjString *) fTagCutSumm->GetValue(kstr->GetString().Data());
    //    printf("Got object value %s %s\n", kstr->GetString().Data(), vstr->GetString().Data());
    char *iTagInfo;
    iTagInfo = strdup(vstr->GetString().Data());
    
    Int_t iAcc = atoi(strtok(iTagInfo, ","));
    Int_t iRej = atoi(strtok(NULL, ","));
    
    if ((iAcc == 0) && ((iRej+iAcc)>0))
      iFilesEmpty++;
  }

  return iFilesEmpty;
  
}

//______________________________________________________________________________
TObject *AliESDInputHandler::GetStatistics(Option_t *option) const
{
// Get the statistics histogram(s) from the physics selection object. This
// should be called during FinishTaskOutput(). Option can be empty (default
// statistics histogram) or BIN0.
   if (!fEventCuts) return NULL;
   TString opt(option);
   opt.ToUpper();
   if (opt=="BIN0") return fEventCuts->GetStatistics("BIN0");
   else return fEventCuts->GetStatistics("ALL");
}   
