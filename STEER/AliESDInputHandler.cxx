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

#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliRunTag.h"
#include "AliEventTag.h"
#include "AliLog.h"

ClassImp(AliESDInputHandler)

//______________________________________________________________________________
AliESDInputHandler::AliESDInputHandler() :
  AliInputEventHandler(),
  fEvent(0x0),
  fBranches(""),
  fBranchesOn(""),
  fAnalysisType(0),
  fUseTags(kFALSE),
  fChainT(0),
  fTreeT(0),
  fRunTag(0)
{
  // default constructor
}

//______________________________________________________________________________
AliESDInputHandler::~AliESDInputHandler() 
{
  //  destructor
  //  delete fEvent;
}

//______________________________________________________________________________
AliESDInputHandler::AliESDInputHandler(const char* name, const char* title):
    AliInputEventHandler(name, title), fEvent(0x0), fBranches(""), fBranchesOn(""), fAnalysisType(0),
     fUseTags(kFALSE), fChainT(0), fTreeT(0), fRunTag(0)
{
    // Constructor
}

Bool_t AliESDInputHandler::Init(TTree* tree,  Option_t* opt)
{
    // Initialisation necessary for each new tree 
    fAnalysisType = opt;
    fTree         = tree;
    
    if (!fTree) return kFALSE;
    // Get pointer to ESD event
    SwitchOffBranches();
    SwitchOnBranches();
    
    if (fEvent) {
      delete fEvent;
      fEvent = 0;
    }
    fEvent = new AliESDEvent();

    fEvent->ReadFromTree(fTree);
    return kTRUE;
}

Bool_t AliESDInputHandler::BeginEvent(Long64_t /*entry*/)
{
    // Copy from old to new format if necessary
  AliESD* old = ((AliESDEvent*) fEvent)->GetAliESDOld();
  if (old) {
	((AliESDEvent*)fEvent)->CopyFromOldESD();
	old->Reset();
  }
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
    if (!fUseTags) return (kTRUE);
    
    Bool_t zip = kFALSE;
    
    TString fileName(path);
    if(fileName.Contains("#AliESDs.root")){
	fileName.ReplaceAll("#AliESDs.root", "");
	zip = kTRUE;
    } 
    else if (fileName.Contains("AliESDs.root")){
	fileName.ReplaceAll("AliESDs.root", "");
    }
    else if(fileName.Contains("#AliAOD.root")){
	fileName.ReplaceAll("#AliAOD.root", "");
	zip = kTRUE;
    }
    else if(fileName.Contains("AliAOD.root")){
	fileName.ReplaceAll("AliAOD.root", "");
    }
    else if(fileName.Contains("#galice.root")){
	// For running with galice and kinematics alone...
	fileName.ReplaceAll("#galice.root", "");
	zip = kTRUE;
    }
    else if(fileName.Contains("galice.root")){
	// For running with galice and kinematics alone...
	fileName.ReplaceAll("galice.root", "");
    }

    
   TString* pathName = new TString("./");
   if (fileName.Length() != 0) {
       *pathName = fileName;
   }

    printf("AliESDInputHandler::Notify() Path: %s\n", pathName->Data());
    
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
	TFile* file = TFile::Open(fileName.Data());
	TArchiveFile* arch = file->GetArchive();
	TObjArray* arr = arch->GetMembers();
	TIter next(arr);
	
	while ((file = (TFile*) next())) {
	    name = file->GetName();
	    if (strstr(name,tagPattern)) { 
		tagFilename = pathName->Data();
		tagFilename += "#";
		tagFilename += name;
		fChainT->Add(tagFilename);  
		AliInfo(Form("Adding %s to tag chain \n", tagFilename.Data()));
	    }//pattern check
	} // archive file loop
    } else {
	void * dirp = gSystem->OpenDirectory(pathName->Data());
	while((name = gSystem->GetDirEntry(dirp))) {
	    if (strstr(name,tagPattern)) { 
		tagFilename = pathName->Data();
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


void AliESDInputHandler::SwitchOffBranches() const {
  //
  // Switch of branches on user request
    TObjArray * tokens = fBranches.Tokenize(" ");
    Int_t ntok = tokens->GetEntries();
    for (Int_t i = 0; i < ntok; i++)  {
	TString str = ((TObjString*) tokens->At(i))->GetString();
	if (str.Length() == 0)
	    continue;
	fTree->SetBranchStatus(Form("%s%s%s","*", str.Data(), "*"), 0);
	AliInfo(Form("Branch %s switched off \n", str.Data()));
    }
}

void AliESDInputHandler::SwitchOnBranches() const {
  //
  // Switch of branches on user request
  TObjArray * tokens = fBranchesOn.Tokenize(" ");
  Int_t ntok = tokens->GetEntries();

  for (Int_t i = 0; i < ntok; i++)  {
      TString str = ((TObjString*) tokens->At(i))->GetString();
      if (str.Length() == 0)
	  continue;
      fTree->SetBranchStatus(Form("%s%s%s","*", str.Data(), "*"), 1);
      AliInfo(Form("Branch %s switched on \n", str.Data()));
  }
}
