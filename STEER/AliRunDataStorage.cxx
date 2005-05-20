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
#include "AliSelectionMetaData.h"
#include "AliObjectMetaData.h"
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
// (AliRunDataStorage is NOT the owner of the returned object)
// name must be in the form "Detector/DBType/DetSpecType"
// es: "ZDC/Calib/Pedestals"

  AliSelectionMetaData defaultMetaData;
  AliSelectionMetaData* selectedMetaData = &defaultMetaData;

  // look for a meta data selection
  for (Int_t i = 0; i < fSelection.GetEntriesFast(); i++) {
    AliSelectionMetaData* selection = (AliSelectionMetaData*) fSelection[i];
    if (!selection) continue;
    if (selection->Matches(name, runNumber)) {
      selectedMetaData = selection;
    }
  }

  // get the entry
  AliSelectionMetaData selMetaData(*selectedMetaData);
  selMetaData.SetName(name);
  AliRunData* entry = GetEntry(selMetaData, runNumber);
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

  // record entry to a file (in the same way as AliRunDataFile::PutEntry, 
  // so that the file can be opened as a AliRunDataFile!)

  if (fRecordFile) {
    fRecordFile->cd();
    TDirectory* saveDir = gDirectory;

    // go to or create the directory
    TString strname(name);
    while (strname.BeginsWith("/")) strname.Remove(0);
    TDirectory* dir = fRecordFile;
    Int_t index = -1;
    while ((index = strname.Index("/")) >= 0) {
      TString dirName(strname(0, index));
      if ((index > 0) && !dir->Get(dirName)) dir->mkdir(dirName);
      dir->cd(dirName);
      dir = gDirectory;
      strname.Remove(0, index+1);
    } 

    entry->Write(strname);
    if (saveDir) saveDir->cd(); else gROOT->cd();
  
  }
  
  return (entry->GetObject())->Clone();

}


//_____________________________________________________________________________
Bool_t AliRunDataStorage::Put(const TObject* object, 
			      const AliObjectMetaData& objMetaData)
{
// put an object into the data base
// (AliRunDataStorage does not adopt the object)
// location of where the object is stored is defined by 
// the AliObjectMetaData's name ("Detector/DBType/DetSpecType")
// and run Range. Storage is handled by the PutEntry method
// of the current AliRunDataStorage instance. 

  if (!object) return kFALSE;
  AliRunData entry(object->Clone(), objMetaData);
  return PutEntry(&entry);
}

//_____________________________________________________________________________
Bool_t AliRunDataStorage::PutEntry(AliRunData* entry)
{
// put an object into the data base
// Refer to the specific method of the current AliRunDataStorage instance
// (AliRunDataFile, AliRunDataOrganizedFile, AliRunDataAlien)

  if (!entry) return kFALSE;
  AliError(Form("This is a read only data base. "
		"The object %s was not inserted", entry->GetName()));
  return kFALSE;
}


//_____________________________________________________________________________
void AliRunDataStorage::Select(const AliSelectionMetaData& selMetaData)
{
// add some meta data selection criteria

  fSelection.Add(new AliSelectionMetaData(selMetaData));
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
// return the current instance of the DB (AliRunDataFile, AliRunOrganizedDataFile...)
// Example of usage: after creating an istance of AliRunDataStorage:
// AliRunDataStorage::Instance()->Get(...)

  return fgInstance;
}


//_____________________________________________________________________________
const AliObjectMetaData& AliRunDataStorage::GetObjectMetaData(const char* name)
{
// Returns the object's metadata of the already retrieved object
// (useful, for example, if you want to know the format of the object you have 
// retrieved)

AliRunData *entry = (AliRunData*) fEntries.FindObject(name);
 if(!entry){
    AliError(Form("Entry %s not found! You make me crash!",name));
 }
 return entry->GetObjectMetaData(); 

}
