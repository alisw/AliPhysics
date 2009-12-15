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
#include <TObjArray.h>
#include <TProcessID.h>
#include <TSystem.h>

#include "AliESDInputHandlerRPITS.h"
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliLog.h"

ClassImp(AliESDInputHandlerRPITS)

//______________________________________________________________________________
AliESDInputHandlerRPITS::AliESDInputHandlerRPITS() :
    AliESDInputHandlerRP()
{
  // Default constructor
}


//______________________________________________________________________________
AliESDInputHandlerRPITS::AliESDInputHandlerRPITS(const char* name, const char* title):
    AliESDInputHandlerRP(name, title)
{
    // Constructor
}

//______________________________________________________________________________
AliESDInputHandlerRPITS::~AliESDInputHandlerRPITS() 
{
  // Destructor
}


Bool_t AliESDInputHandlerRPITS::Notify(const char *path)
{
  // Notify about directory change
  // The directory is taken from the 'path' argument
  // 
    AliInfo(Form("Directory change %s \n", path));
    // Get path to directory
    TString fileName(path);

    if(fileName.Contains("#")){
    // If this is an archive it will contain a # 
      fIsArchive = kTRUE;
    } else  if(fileName.Contains("AliESDs.root")){
      fileName.ReplaceAll("AliESDs.root", "");
    }

    //
    // At this point we have a path to the directory or to the archive anchor
    *fPathName = fileName;
    //
    // Now filter the files containing RecPoints *.RecPoints.*

    TSeqCollection* members;

    
    if (fIsArchive) {
	// Archive
      TFile* file = TFile::Open(fPathName->Data());
      TArchiveFile* arch = file->GetArchive();
      members = arch->GetMembers();
    } else {
	// Directory or alien archive
      if (fileName.BeginsWith("alien:")) {
        TFile* file = TFile::Open(Form("%s/root_archive.zip", fPathName->Data()));
        TArchiveFile* arch = file->GetArchive();
        members = arch->GetMembers();
      } else {  
        TString wd = gSystem->WorkingDirectory();
        TSystemDirectory dir(".", fPathName->Data());
        members = dir.GetListOfFiles();
        gSystem->cd(wd);
      }  
    }

    TIter next(members);
    TFile* entry;
    Int_t ien = 0;
    fDetectors->Delete();
    
    while ( (entry = (TFile*) next()) )
    {
	TString name(entry->GetName());
	TObjArray* tokens = name.Tokenize(".");
	Int_t ntok = tokens->GetEntries();
	if (ntok <= 1) continue;
	TString str = ((TObjString*) tokens->At(1))->GetString();
	if (!(strcmp(str.Data(), "RecPoints"))){
	    TString det = ((TObjString*) tokens->At(0))->GetString();
	    if (!(strcmp(det.Data(), "ITS"))) {
		printf("Found file with RecPoints for %s \n", det.Data());
		TNamed* ent = new TNamed(det.Data(), det.Data());
		fRTrees->AddAt(0, ien);
		ent->SetUniqueID(ien++);
		fDetectors->Add(ent);
	    }
	}
	if(tokens) delete tokens;
    } // loop over files


    // Now we have the path and the list of detectors
    
    printf("AliESDInputHandlerRPITS::Notify() Path: %s\n", fPathName->Data());
    //
    ResetIO();
    InitIO("");
    // Some clean-up
    members->Delete();

    return kTRUE;
}

