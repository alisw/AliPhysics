/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// access classes for a data base in a (local) file                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>
#include "AliLog.h"
#include "AliRunData.h"
#include "AliRunDataFile.h"


ClassImp(AliRunDataFile)


//_____________________________________________________________________________
AliRunDataFile::AliRunDataFile(const char* fileName, Bool_t readOnly) :
  AliRunDataStorage(),
  fFile(NULL)
{
// constructor

  if (!fileName) {
    AliError("no file name given");
    return;
  }
  TDirectory* saveDir = gDirectory;
  fFile = TFile::Open(fileName, ((readOnly) ? "READ" : "UPDATE"));
  if (saveDir) saveDir->cd(); else gROOT->cd();
  if (!fFile || !fFile->IsOpen()) {
    AliError(Form("could not open file %s", fileName));
    fFile = NULL;
  }
}

//_____________________________________________________________________________
AliRunDataFile::~AliRunDataFile()
{
// destructor

  if (fFile) {
    fFile->Close();
    delete fFile;
  }
}

//_____________________________________________________________________________
AliRunDataFile::AliRunDataFile(const AliRunDataFile& /*db*/) :
  AliRunDataStorage(),
  fFile(NULL)
{
// copy constructor

  AliFatal("not implemented");
}

//_____________________________________________________________________________
AliRunDataFile& AliRunDataFile::operator = (const AliRunDataFile& /*db*/)
{
// assignment operator

  AliFatal("not implemented");
  return *this;
}


//_____________________________________________________________________________
AliRunData* AliRunDataFile::GetEntry(AliMetaData& metaData, Int_t runNumber)
{
// get an object from the data base

  // go to the directory
  TDirectory* saveDir = gDirectory;
  TString name(metaData.GetName());
  Int_t last = name.Last('/');
  if (last < 0) {
    fFile->cd();
  } else {
    TString dirName(name(0, last));
    if (!fFile->cd(dirName)) {
      AliDebug(1, Form("no directory %s found", dirName.Data()));
      if (saveDir) saveDir->cd(); else gROOT->cd();
      return NULL;
    }
    name.Remove(0, last+1);
  }

  TKey* key = fFile->GetKey(name);
  if (!key) {
    AliDebug(1, Form("no object with name %s found", metaData.GetName()));
    if (saveDir) saveDir->cd(); else gROOT->cd();
    return NULL;
  }
  Int_t nCycles = key->GetCycle();

  // find the closest entry
  AliRunData* closestEntry = NULL;
  for (Int_t iCycle = nCycles; iCycle > 0; iCycle--) {
    key = fFile->GetKey(name, iCycle);
    if (!key) continue;
    AliRunData* entry = (AliRunData*) key->ReadObj();
    if (!entry) continue;
    if (!entry->InheritsFrom(AliRunData::Class())) {
      AliMetaData metaData;
      entry = new AliRunData(entry, metaData);
    }
    if (!entry->GetMetaData().IsValid(runNumber, &metaData) ||
	(entry->Compare(closestEntry) <= 0)) {
      delete entry;
      continue;
    }
    delete closestEntry;
    closestEntry = entry;
  }
  if (saveDir) saveDir->cd(); else gROOT->cd();
  if (!closestEntry) return NULL;
  return closestEntry;
}

//_____________________________________________________________________________
Bool_t AliRunDataFile::PutEntry(AliRunData* entry)
{
// put an object into the data base

  if (!entry || !fFile) return kFALSE;
  if (!fFile->IsWritable()) {
    AliError(Form("The data base file was opened in read only mode. "
		  "The object %s was not inserted", entry->GetName()));
    return kFALSE;
  }
  TDirectory* saveDir = gDirectory;

  // go to or create the directory
  TString name(entry->GetName());
  while (name.BeginsWith("/")) name.Remove(0);
  TDirectory* dir = fFile;
  Int_t index = -1;
  while ((index = name.Index("/")) >= 0) {
    TString dirName(name(0, index));
    if ((index > 0) && !dir->Get(dirName)) dir->mkdir(dirName);
    dir->cd(dirName);
    dir = gDirectory;
    name.Remove(0, index+1);
  } 

  // determine the version number
  Int_t version = 0;
  TKey* key = fFile->GetKey(name);
  if (key) {
    Int_t nCycles = key->GetCycle();
    for (Int_t iCycle = nCycles; iCycle > 0; iCycle--) {
      key = fFile->GetKey(entry->GetName(), iCycle);
      if (!key) continue;
      AliRunData* oldEntry = (AliRunData*) key->ReadObj();
      if (!oldEntry) continue;
      if (oldEntry->InheritsFrom(AliRunData::Class())) {
	if (version <= oldEntry->GetMetaData().GetVersion()) {
	  version = oldEntry->GetMetaData().GetVersion()+1;
	}
      }
      delete oldEntry;
    }
  }
  entry->SetVersion(version);

  Bool_t result = (entry->Write(name) != 0);
  if (saveDir) saveDir->cd(); else gROOT->cd();
  return result;
}
