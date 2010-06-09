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
  fReadFriends(0),
  fFriendFileName("AliESDfriends.root")
{
  // default constructor
}

//______________________________________________________________________________
AliESDInputHandler::~AliESDInputHandler() 
{
  //  destructor
}

//______________________________________________________________________________
AliESDInputHandler::AliESDInputHandler(const char* name, const char* title):
    AliInputEventHandler(name, title), fEvent(0x0), fFriend(0x0), fESDpid(0x0), fAnalysisType(0),
    fNEvents(0),  fHLTEvent(0x0), fHLTTree(0x0), fUseHLT(kFALSE), fTagCutSumm(0x0), fUseTags(kFALSE), fChainT(0), fTreeT(0), fRunTag(0), fReadFriends(0), fFriendFileName("AliESDfriends.root")
{
    // Constructor
}

Bool_t AliESDInputHandler::Init(TTree* tree,  Option_t* opt)
{
    //
    // Initialisation necessary for each new tree 
    // 
    fAnalysisType = opt;
    fTree         = tree;
    
    if (!fTree) return kFALSE;
    fTree->GetEntry(0);
    

    if (!fEvent) fEvent = new AliESDEvent();
    fEvent->ReadFromTree(fTree);
    fNEvents = fTree->GetEntries();
    return kTRUE;
}

Bool_t AliESDInputHandler::BeginEvent(Long64_t entry)
{
    
    // Copy from old to new format if necessary
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
  if (fEventCuts)
    fIsSelected = fEventCuts->IsSelected((AliESDEvent*)fEvent); 
  //
  // Friends
  ((AliESDEvent*)fEvent)->SetESDfriend(fFriend);
  return kTRUE;
}

Bool_t  AliESDInputHandler::FinishEvent()
{
    // Finish the event 
    if(fEvent)fEvent->Reset();
    return kTRUE;
} 

Bool_t AliESDInputHandler::Notify(const char* path)
{
    // Notify a directory change
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
      cTree->SetBranchAddress("ESDfriend.", &fFriend);
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




    if (!fUseTags) return (kTRUE);
    
    Bool_t zip = kFALSE;
    
    TString fileName(path);
    if(fileName.Contains("#AliESDs.root")){
	zip = kTRUE;
    } 
    else if (fileName.Contains("AliESDs.root")){
	fileName.ReplaceAll("AliESDs.root", "");
    }
    else if(fileName.Contains("#AliAOD.root")){
	zip = kTRUE;
    }
    else if(fileName.Contains("AliAOD.root")){
	fileName.ReplaceAll("AliAOD.root", "");
    }
    else if(fileName.Contains("#galice.root")){
	// For running with galice and kinematics alone...
	zip = kTRUE;
    }
    else if(fileName.Contains("galice.root")){
	// For running with galice and kinematics alone...
	fileName.ReplaceAll("galice.root", "");
    }

    
    TString pathName("./");
    if (fileName.Length() != 0) {
	pathName = fileName;
    }
    
    printf("AliESDInputHandler::Notify() Path: %s\n", pathName.Data());

    if (fRunTag) {
	fRunTag->Clear();
    } else {
	fRunTag = new AliRunTag();
    }
    
    delete fTreeT; fTreeT = 0;
    
    if (fChainT) {
	delete fChainT;
	fChainT = 0;
    }
    
    if (!fChainT) {
	fChainT = new TChain("T");
    }
    


    const char* tagPattern = "ESD.tag.root";
    const char* name = 0x0;
    TString tagFilename;
    if (zip) {
	TFile* file = fTree->GetCurrentFile();
	TArchiveFile* arch = file->GetArchive();
	TObjArray* arr = arch->GetMembers();
	TIter next(arr);
	
	while ((file = (TFile*) next())) {
	    name = file->GetName();
	    if (strstr(name,tagPattern)) { 
		tagFilename = pathName.Data();
		tagFilename += "#";
		tagFilename += name;
		fChainT->Add(tagFilename);  
		AliInfo(Form("Adding %s to tag chain \n", tagFilename.Data()));
	    }//pattern check
	} // archive file loop
    } else {
	void * dirp = gSystem->OpenDirectory(pathName.Data());
	while((name = gSystem->GetDirEntry(dirp))) {
	    if (strstr(name,tagPattern)) { 
		tagFilename = pathName.Data();
		tagFilename += "/";
		tagFilename += name;
		fChainT->Add(tagFilename);  
		AliInfo(Form("Adding %s to tag chain \n", tagFilename.Data()));
	    }//pattern check
	}//directory loop
    }
    fChainT->SetBranchAddress("AliTAG",&fRunTag);
    fChainT->GetEntry(0);
    return kTRUE;
}



Option_t *AliESDInputHandler::GetDataType() const
{
// Returns handled data type.
   return gESDDataType;
}

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
Bool_t AliESDInputHandler::GetCutSummaryForChain(Int_t *aTotal, Int_t *aAccepted, Int_t *aRejected)
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
