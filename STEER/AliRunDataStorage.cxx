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
// base class for data base access classes                                   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TFile.h>
#include <TKey.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliMetaData.h"
#include "AliRunData.h"
#include "AliRunDataStorage.h"


ClassImp(AliRunDataStorage)


AliRunDataStorage* AliRunDataStorage::fgInstance = NULL;


//_____________________________________________________________________________
AliRunDataStorage::AliRunDataStorage() :
  TObject(),
  fSelection(),
  fEntries(),
  fRecordFile(NULL)
{
// default constructor

  if (fgInstance) delete fgInstance;
  fgInstance = this;
}

//_____________________________________________________________________________
AliRunDataStorage::~AliRunDataStorage()
{
// destructor

  fSelection.Delete();
  fEntries.Delete();
  if (fRecordFile) {
    fRecordFile->Close();
    delete fRecordFile;
  }
  fgInstance = NULL;
}

//_____________________________________________________________________________
AliRunDataStorage::AliRunDataStorage(const AliRunDataStorage& db) :
  TObject(db),
  fSelection(),
  fEntries(),
  fRecordFile(NULL)
{
// copy constructor

  AliFatal("not implemented");
}

//_____________________________________________________________________________
AliRunDataStorage& AliRunDataStorage::operator = (const AliRunDataStorage& /*db*/)
{
// assignment operator

  AliFatal("not implemented");
  return *this;
}


//_____________________________________________________________________________
const TObject* AliRunDataStorage::Get(const char* name, Int_t runNumber)
{
// get an object from the data base
// (AliRunDataStorage is the owner of the returned object)

  AliMetaData defaultMetaData;
  AliMetaData* selectedMetaData = &defaultMetaData;

  // look for a meta data selection
  for (Int_t i = 0; i < fSelection.GetEntriesFast(); i++) {
    AliMetaData* selection = (AliMetaData*) fSelection[i];
    if (!selection) continue;
    if (selection->Matches(name, runNumber)) {
      selectedMetaData = selection;
    }
  }

  // get the entry
  AliMetaData metaData(*selectedMetaData);
  metaData.SetName(name);
  AliRunData* entry = GetEntry(metaData, runNumber);
  if (entry) {
    AliDebug(2, "got the entry:");
    ToAliDebug(2, entry->Dump());
  } else {
    AliDebug(2, Form("got no entry for %s", name));
  }

  // update array of current entries
  if (!entry) return NULL;
  TObject* oldEntry = fEntries.FindObject(entry->GetName());
  if (oldEntry) {
    delete fEntries.Remove(oldEntry);
  }
  fEntries.Add(entry);

  // record entry to a file
  if (fRecordFile) {
    Bool_t isAlreadyRecorded = kFALSE;
    TDirectory* dir = gDirectory;
    fRecordFile->cd();
    TKey* key = fRecordFile->GetKey(entry->GetName());
    if (key) {
      Int_t nCycles = key->GetCycle();
      for (Int_t iCycle = nCycles; iCycle > 0; iCycle--) {
	key = fRecordFile->GetKey(entry->GetName(), iCycle);
	if (!key) continue;
	AliRunData* recEntry = (AliRunData*) key->ReadObj();
	if (!recEntry) continue;
	if (recEntry->InheritsFrom(AliRunData::Class()) && 
	    (recEntry->GetMetaData() == entry->GetMetaData())) {
	  isAlreadyRecorded = kTRUE;
	}
	delete recEntry;
	if (isAlreadyRecorded) break;
      }
    }
    if (!isAlreadyRecorded) {
      if (entry->Write() == 0) {
	AliError(Form("could not record entry %s", entry->GetName()));
      }
    }
    if (dir) dir->cd(); else gROOT->cd();
  }

  return entry->GetObject();
}


//_____________________________________________________________________________
Bool_t AliRunDataStorage::Put(const TObject* object, 
			      const AliMetaData& metaData)
{
// put an object into the data base
// (AliRunDataStorage does not adopt the object)

  if (!object) return kFALSE;
  AliRunData entry(object->Clone(), metaData);
  return PutEntry(&entry);
}

//_____________________________________________________________________________
Bool_t AliRunDataStorage::PutEntry(AliRunData* entry)
{
// put an object into the data base

  if (!entry) return kFALSE;
  AliError(Form("This is a read only data base. "
		"The object %s was not inserted", entry->GetName()));
  return kFALSE;
}


//_____________________________________________________________________________
void AliRunDataStorage::Select(const AliMetaData& metaData)
{
// add some meta data selection criteria

  fSelection.Add(new AliMetaData(metaData));
}


//_____________________________________________________________________________
Bool_t AliRunDataStorage::RecordToFile(const char* fileName)
{
// record entries retrieved from the data base to a file with the given name

  if (fRecordFile) {
    fRecordFile->Close();
    delete fRecordFile;
  }

  TDirectory* dir = gDirectory;
  fRecordFile = TFile::Open(fileName, "UPDATE");
  if (dir) dir->cd(); else gROOT->cd();
  if (!fRecordFile || !fRecordFile->IsOpen()) {
    AliError(Form("could not open file %s", fileName));
    delete fRecordFile;
    fRecordFile = NULL;
    return kFALSE;
  }
  return kTRUE;
}


//_____________________________________________________________________________
AliRunDataStorage* AliRunDataStorage::Instance()
{
// return the current instance of the DB

  return fgInstance;
}
