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
//     Event handler for ESD input reading the RecPoint Trees in parallel
//     Author: Andreas Morsch, CERN
//-------------------------------------------------------------------------

#include <TTree.h>
#include <TList.h>
#include <TFile.h>
#include <TArchiveFile.h>
#include <TSystemDirectory.h>
#include <TString.h>
#include <TObjString.h>
#include <TProcessID.h>

#include "AliESDInputHandlerRP.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliLog.h"

ClassImp(AliESDInputHandlerRP)

//______________________________________________________________________________
AliESDInputHandlerRP::AliESDInputHandlerRP() :
    AliESDInputHandler(),
    fRTrees(   new TList()),
    fRFiles(   new TList()),
    fDetectors(new TList()),
    fDirR(0),
    fEventNumber(-1),
    fFileNumber(0),
    fEventsPerFile(0),
    fExtension(""),
    fPathName(new TString("./")),
    fIsArchive(kFALSE)
{
  // Default constructor
}


//______________________________________________________________________________
AliESDInputHandlerRP::AliESDInputHandlerRP(const char* name, const char* title):
    AliESDInputHandler(name, title),
    fRTrees(   new TList()),
    fRFiles(   new TList()),
    fDetectors(new TList()),
    fDirR(0),
    fEventNumber(-1),
    fFileNumber(0),
    fEventsPerFile(0),
    fExtension(""),
    fPathName(new TString("./")),
    fIsArchive(kFALSE)
{
    // Constructor
}

//______________________________________________________________________________
AliESDInputHandlerRP::~AliESDInputHandlerRP() 
{
  // Destructor
}

Bool_t AliESDInputHandlerRP::Init(Option_t* opt)
{
    //
    // Initialize input
    //
    if (!(strcmp(opt, "proof")) || !(strcmp(opt, "local"))) return kTRUE;
    //
    TIter next(fDetectors);
    TNamed* det;
    TFile* file = 0;
    while ((det = (TNamed*) next()))
    {
	if (!fIsArchive) {
	    file = TFile::Open(Form("%s%s.RecPoints.root", fPathName->Data(), det->GetName()));
	} else {
	    file = TFile::Open(Form("%s#%s.RecPoints.root", fPathName->Data(), det->GetName()));
	}
	if (!file) AliFatal(Form("AliESDInputHandlerRP: %s.RecPoints.root not found in %s ! \n", det->GetName(), fPathName->Data()));
	fRFiles->Add(file);
    }

    if (file) {
	fEventsPerFile = file->GetNkeys() - file->GetNProcessIDs();
    } else {
	AliFatal(Form("AliESDInputHandlerRP: No file with RecPoints found in %s ! \n", fPathName->Data()));
    }
    

    // Reset the event number
    fEventNumber      = -1;
    fFileNumber       =  0;
    // Get number of events from esd tree 
    printf("AliESDInputHandler::Init() %d %d\n",__LINE__, fNEvents);
    return kTRUE;
}

Bool_t AliESDInputHandlerRP::BeginEvent(Long64_t entry)
{
    // Begin the next event
    // Delegate first to base class
    AliESDInputHandler::BeginEvent(entry);
//
    if (entry == -1) {
	fEventNumber++;
	entry = fEventNumber;
    } else {
	fEventNumber = entry;
    }
    
    if (entry >= fNEvents) {
	AliWarning(Form("AliESDInputHandlerRP: Event number out of range %5d %5d\n", entry, fNEvents));
	return kFALSE;
    }
    return LoadEvent(entry);
}

Bool_t AliESDInputHandlerRP::LoadEvent(Int_t iev)
{
    // Load the event number iev
    //
    // Calculate the file number
    Int_t inew  = iev / fEventsPerFile;
    if (inew != fFileNumber) {
	fFileNumber = inew;
	if (!OpenFile(fFileNumber)){
	    return kFALSE;
	}
    }
    // Folder name
    char folder[20];
    sprintf(folder, "Event%d", iev);
    // Tree R
    TIter next(fRFiles);
    TFile* file;
    Int_t idx = 0;
    while ((file = (TFile*) next()))
    {
	file->GetObject(folder, fDirR);
	if (!fDirR) {
	    AliWarning(Form("AliESDInputHandlerRP: Event #%5d not found\n", iev));
	    return kFALSE;
	}
	TTree* tree = 0;
	fDirR ->GetObject("TreeR", tree);
	fRTrees->AddAt(tree, idx++);
    }
    return kTRUE;
}

Bool_t AliESDInputHandlerRP::OpenFile(Int_t i)
{
    // Open file i
    Bool_t ok = kTRUE;
    if (i > 0) {
	fExtension = Form("%d", i);
    } else {
	fExtension = "";
    }
    
    fRFiles->Delete();
    TIter next(fDetectors);
    TNamed* det;
    TFile* file;
    while ((det = (TNamed*) next()))
    {
	if (!fIsArchive) {
	    file = TFile::Open(Form("%s%s.RecPoints%s.root", fPathName->Data(), det->GetName(), fExtension));
	} else {
	    file = TFile::Open(Form("%s#%s.RecPoints%s.root", fPathName->Data(), det->GetName(), fExtension));
	}
	if (!file) AliFatal(Form("AliESDInputHandlerRP: RecPoints.root not found in %s ! \n", fPathName->Data()));
	fRFiles->Add(file);
    }
    return ok;
}

Bool_t AliESDInputHandlerRP::Notify(const char *path)
{
  // Notify about directory change
  // The directory is taken from the 'path' argument
  // 
    // Get path to directory
    TString fileName(path);
    if(fileName.Contains("AliESDs.root")){
	fileName.ReplaceAll("AliESDs.root", "");
    }
    // If this is an archive it will contain a # 
    if(fileName.Contains("#")){
	fileName.ReplaceAll("#", "");
    }
    //
    // At this point we have a path to the directory or to the archive
    *fPathName = fileName;
    //
    // Now filter the files containing RecPoints *.RecPoints.*
    fIsArchive = kFALSE;
    if (fPathName->Contains(".zip")) fIsArchive = kTRUE;

    TSeqCollection* members;
    
    if (fIsArchive) {
	// Archive
	TFile* file = TFile::Open(fPathName->Data());
	TArchiveFile* arch = file->GetArchive();
	members = arch->GetMembers();
    } else {
	// Directory
	TSystemDirectory dir(".", fPathName->Data());
	members = dir.GetListOfFiles();
    }
    
    TIter next(members);
    TFile* entry;
    Int_t ien = 0;
    fRTrees->Clear();
    
    while ( (entry = (TFile*) next()) )
    {
	printf("File %s \n", entry->GetName());
	TString name(entry->GetName());
	TObjArray* tokens = name.Tokenize(".");
	Int_t ntok = tokens->GetEntries();
	if (ntok <= 1) continue;
	TString str = ((TObjString*) tokens->At(1))->GetString();
	if (!(strcmp(str.Data(), "RecPoints"))){
	    TString det = ((TObjString*) tokens->At(0))->GetString();
	    printf("Name  %s \n", det.Data());
	    TNamed* ent = new TNamed(det.Data(), det.Data());
	    ent->SetUniqueID(ien++);
	    fDetectors->Add(ent);
	}
    } // loop over files


    // Now we have the path and the list of detectors
    
    printf("AliESDInputHandlerRP::Notify() Path: %s\n", fPathName->Data());
    //
    ResetIO();
    InitIO("");
    // Some clean-up
    members->Delete();

    return kTRUE;
}

Bool_t AliESDInputHandlerRP::FinishEvent()
{
    // Clean-up after each event
    delete fDirR;  fDirR = 0;
    AliESDInputHandler::FinishEvent();
    return kTRUE;
}

void AliESDInputHandlerRP::ResetIO()
{
// Delete trees and files
//    fRTrees->Clear();
    fRFiles->Clear();
    fExtension="";
}

TTree* AliESDInputHandlerRP::GetTreeR(char* det)
{
// Return pointer to RecPoint tree for detector det
    TNamed* entry = (TNamed*) (fDetectors->FindObject(det));
    if (!entry) {
	AliWarning(Form("AliESDInputHandlerRP: No RecPoints for detector %s available \n", det));
	return 0;
    } else {
	Int_t ien = entry->GetUniqueID();
	return ((TTree*) (fRTrees->At(ien)));
    }
}
